import unittest
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from nasa_mouse_generative.archs4_catalog import (
    classify_eligibility,
    classify_tissues,
    select_balanced,
)
from nasa_mouse_generative.config import (
    BenchmarkConfig,
    PreprocessingConfig,
    TrainingConfig,
)
from nasa_mouse_generative.effect_validation import compare_real_synthetic_effects
from nasa_mouse_generative.experiment_plan import expand_matrix
from nasa_mouse_generative.metrics import (
    _write_accession_effect_validation,
    _write_per_tissue_generation_metrics,
)
from nasa_mouse_generative.preprocessing import FittedPreprocessor, apply_stats
from nasa_mouse_generative.osdr_expression import _read_accession_block
from nasa_mouse_generative.split_plan import build_pooled_plan
from nasa_mouse_generative.tissues import canonicalize_material


def _archs4_row(**overrides):
    row = {
        "geo_accession": "GSM1",
        "series_id": "GSE1",
        "title": "normal liver",
        "source_name_ch1": "liver",
        "characteristics_ch1": "adult wild type control",
        "library_strategy": "RNA-Seq",
        "library_source": "TRANSCRIPTOMIC",
        "organism_ch1": "Mus musculus",
        "singlecellprobability": 0.01,
    }
    row.update(overrides)
    return row


class TissueTests(unittest.TestCase):
    def test_material_aliases_and_specificity(self):
        self.assertEqual(canonicalize_material("Cells, Cultured"), "cultured_cells")
        self.assertEqual(canonicalize_material("Bone-marrow"), "bone_marrow")
        self.assertEqual(canonicalize_material("Optic nerve"), "optic_nerve")
        self.assertEqual(canonicalize_material("Tibialis anterior"), "skeletal_muscle")


class Archs4CatalogTests(unittest.TestCase):
    def test_bulk_health_and_control_cohorts(self):
        frame = pd.DataFrame(
            [
                _archs4_row(),
                _archs4_row(
                    geo_accession="GSM2",
                    title="liver sample",
                    characteristics_ch1="adult sample",
                ),
                _archs4_row(
                    geo_accession="GSM3",
                    title="liver tumor",
                    characteristics_ch1="cancer treatment",
                ),
                _archs4_row(
                    geo_accession="GSM4",
                    title="WT_scRNA_batch1",
                    singlecellprobability=0.1,
                ),
                _archs4_row(
                    geo_accession="GSM5",
                    title="normal liver spaceflight sample",
                ),
                _archs4_row(
                    geo_accession="GSM6",
                    title="Sik2/3 DKO",
                ),
                _archs4_row(
                    geo_accession="GSM7",
                    title="embryonic liver E10.5",
                ),
                _archs4_row(
                    geo_accession="GSM8",
                    title="LTHSC_WT_30",
                    characteristics_ch1="adult control protocol SmartSeqV2",
                ),
                _archs4_row(
                    geo_accession="GSM9",
                    title="ECKO1",
                    characteristics_ch1="adult liver endothelial cells",
                ),
            ]
        )
        classified = classify_eligibility(classify_tissues(frame, ["liver"]))
        by_id = classified.set_index("geo_accession")

        self.assertTrue(by_id.loc["GSM1", "eligible_control_only"])
        self.assertTrue(by_id.loc["GSM2", "eligible_healthy_preferred"])
        self.assertFalse(by_id.loc["GSM2", "eligible_control_only"])
        self.assertTrue(by_id.loc["GSM3", "eligible_broad"])
        self.assertFalse(by_id.loc["GSM3", "eligible_healthy_preferred"])
        self.assertFalse(by_id.loc["GSM4", "bulk_like"])
        self.assertFalse(by_id.loc["GSM5", "eligible_broad"])
        self.assertFalse(by_id.loc["GSM6", "eligible_healthy_preferred"])
        self.assertFalse(by_id.loc["GSM7", "eligible_healthy_preferred"])
        self.assertFalse(by_id.loc["GSM8", "bulk_like"])
        self.assertFalse(by_id.loc["GSM9", "eligible_healthy_preferred"])

    def test_balancing_caps_series_and_equalizes_tissue_weight(self):
        rows = []
        for tissue in ("liver", "kidney"):
            for series in ("GSE1", "GSE2"):
                for index in range(4):
                    rows.append(
                        {
                            "geo_accession": f"{tissue}-{series}-{index}",
                            "series_id": series,
                            "canonical_tissue": tissue,
                            "health_status": "explicit_control_like",
                            "eligible": True,
                        }
                    )
        selected = select_balanced(
            pd.DataFrame(rows),
            eligibility_column="eligible",
            max_per_tissue=3,
            max_per_series=2,
            seed=7,
        )
        self.assertLessEqual(selected.groupby("canonical_tissue").size().max(), 3)
        self.assertLessEqual(
            selected.groupby(["canonical_tissue", "series_id"]).size().max(), 2
        )
        weights = selected.groupby("canonical_tissue")[
            "hierarchical_sampling_weight"
        ].sum()
        np.testing.assert_allclose(weights.to_numpy(), [0.5, 0.5])


class PreprocessingTests(unittest.TestCase):
    def test_training_fold_scaler_is_reused(self):
        spec = PreprocessingConfig(
            input_units="raw_counts",
            library_normalization="none",
            transform="none",
            scaler="zscore",
            harmonization="none",
        )
        processor = FittedPreprocessor(spec)
        train = np.asarray([[1.0, 2.0], [3.0, 8.0], [5.0, 4.0]])
        fitted = processor.fit_transform(train, ["A", "A", "B"])
        np.testing.assert_allclose(fitted.mean(axis=0), 0.0, atol=1e-6)

        test = np.asarray([[20.0, 30.0]])
        observed = processor.transform(test, ["unseen"])
        expected = apply_stats(test, processor.final_stats)
        np.testing.assert_allclose(observed, expected)

    def test_unseen_study_uses_training_global_fallback(self):
        spec = PreprocessingConfig(
            input_units="raw_counts",
            library_normalization="none",
            transform="none",
            scaler="none",
            harmonization="within_study_zscore",
            unseen_study_policy="global_train_fallback",
        )
        processor = FittedPreprocessor(spec)
        train = np.asarray([[1.0], [3.0], [7.0], [9.0]])
        processor.fit_transform(train, ["A", "A", "B", "B"])
        test = np.asarray([[20.0], [24.0]])
        observed = processor.transform(test, ["C", "C"])
        expected = apply_stats(test, processor.global_stats)
        np.testing.assert_allclose(observed, expected)
        self.assertGreater(abs(float(observed.mean())), 1.0)


class EffectValidationTests(unittest.TestCase):
    def test_unified_evaluator_writes_per_tissue_generation_metrics(self):
        rng = np.random.default_rng(8)
        real = rng.normal(size=(12, 8)).astype(np.float32)
        synthetic = (real + rng.normal(scale=0.1, size=real.shape)).astype(np.float32)
        samples = pd.DataFrame(
            {
                "tissue": ["liver"] * 6 + ["skin"] * 6,
                "condition": (["flight"] * 3 + ["ground_control"] * 3) * 2,
            }
        )
        with tempfile.TemporaryDirectory() as directory:
            table, paths = _write_per_tissue_generation_metrics(
                Path(directory),
                real=real,
                synthetic=synthetic,
                samples=samples,
                max_samples=100,
            )
            self.assertEqual(table["tissue"].tolist(), ["liver", "skin"])
            self.assertTrue(Path(paths["table"]).exists())
            self.assertTrue(Path(paths["plot"]).exists())

    def test_unified_evaluator_writes_accession_effect_artifacts(self):
        rows = []
        real = []
        synthetic = []
        for accession_index, accession in enumerate(("OSD-1", "OSD-2", "OSD-3")):
            baseline = float(accession_index + 1)
            for condition, delta in (("ground_control", 0.0), ("flight", 1.0)):
                for replicate in range(2):
                    rows.append(
                        {
                            "accession": accession,
                            "tissue": "liver",
                            "condition": condition,
                        }
                    )
                    jitter = 0.05 * replicate
                    real.append(
                        [baseline + delta + jitter, baseline + 0.5 * delta - jitter]
                    )
                    synthetic.append(
                        [
                            baseline + 0.9 * delta + jitter,
                            baseline + 0.45 * delta - 0.8 * jitter,
                        ]
                    )
        with tempfile.TemporaryDirectory() as directory:
            summary, paths = _write_accession_effect_validation(
                Path(directory),
                real_normalized=np.asarray(real),
                synthetic_normalized=np.asarray(synthetic),
                samples=pd.DataFrame(rows),
                feature_names=["gene_a", "gene_b"],
            )
            self.assertEqual(summary["accessions"], 3)
            self.assertGreater(summary["meta_effect_correlation"], 0.99)
            self.assertTrue(Path(paths["comparison"]).exists())
            self.assertTrue(Path(paths["summary"]).exists())

    def test_accession_aware_real_synthetic_effect_recovery(self):
        rows = []
        real = []
        synthetic = []
        for accession_index, accession in enumerate(("OSD-1", "OSD-2", "OSD-3")):
            baseline = float(accession_index * 2)
            for condition, delta in (("ground_control", 0.0), ("flight", 1.0)):
                for replicate in range(3):
                    jitter = 0.05 * (replicate - 1)
                    rows.append(
                        {
                            "accession": accession,
                            "tissue": "liver",
                            "condition": condition,
                        }
                    )
                    real.append(
                        [baseline + delta + jitter, baseline - 0.5 * delta + jitter]
                    )
                    synthetic.append(
                        [baseline + 0.9 * delta + jitter, baseline - 0.4 * delta + jitter]
                    )
        tables, summary = compare_real_synthetic_effects(
            np.asarray(real),
            np.asarray(synthetic),
            pd.DataFrame(rows),
            ["up", "down"],
        )
        self.assertEqual(summary["accessions"], 3)
        self.assertAlmostEqual(summary["meta_direction_agreement"], 1.0)
        self.assertGreater(summary["meta_effect_correlation"], 0.99)
        self.assertEqual(len(tables["real_per_accession"]), 6)
        self.assertEqual(len(tables["real_leave_one_out"]), 6)

    def test_effect_validation_reports_insufficient_accessions(self):
        samples = pd.DataFrame(
            {
                "accession": ["OSD-1", "OSD-1"],
                "tissue": ["liver", "liver"],
                "condition": ["flight", "ground_control"],
            }
        )
        _, summary = compare_real_synthetic_effects(
            np.ones((2, 1)), np.ones((2, 1)), samples, ["gene"]
        )
        self.assertEqual(
            summary["status"], "insufficient_accessions_with_both_conditions"
        )


class OsdrExpressionTests(unittest.TestCase):
    def test_technical_replicate_counts_are_summed(self):
        metadata = pd.DataFrame(
            [
                {
                    "id.accession": "OSD-1",
                    "id.sample name": "mouse1_techrep1",
                    "condition_inferred": "flight",
                    "tissue_canonical": "liver",
                },
                {
                    "id.accession": "OSD-1",
                    "id.sample name": "mouse1_techrep2",
                    "condition_inferred": "flight",
                    "tissue_canonical": "liver",
                },
            ]
        )
        counts = pd.DataFrame(
            {
                "gene": ["ENSMUSG0001.2", "ENSMUSG0002"],
                "folder/mouse1_techrep1": [1, 3],
                "folder/mouse1_techrep2": [2, 4],
            }
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "counts.csv"
            counts.to_csv(path, index=False)
            block, retained, missing = _read_accession_block(
                path, metadata, technical_replicate_policy="sum"
            )
        np.testing.assert_allclose(block.iloc[:, 0].to_numpy(), [3, 7])
        self.assertEqual(block.index.tolist(), ["ENSMUSG0001", "ENSMUSG0002"])
        self.assertEqual(retained[0]["api_profile_count"], 2)
        self.assertEqual(missing, [])


class PlanningTests(unittest.TestCase):
    def test_genejepa_generation_is_rejected(self):
        config = BenchmarkConfig(
            training=TrainingConfig(model="genejepa", task="conditional_generation")
        )
        with self.assertRaisesRegex(ValueError, "does not generate expression"):
            config.validate()

    def test_matrix_marks_non_generator_as_blocked(self):
        matrix = {
            "phases": [
                {
                    "name": "test",
                    "axes": {
                        "model": ["genejepa"],
                        "task": ["conditional_generation"],
                    },
                }
            ]
        }
        row = expand_matrix(matrix).iloc[0]
        self.assertEqual(row["status"], "capability_blocked")

    def test_pooled_split_keeps_every_tissue_condition_in_training(self):
        rows = []
        for accession in ("OSD-1", "OSD-2", "OSD-3", "OSD-4"):
            for tissue in ("liver", "kidney"):
                for condition in ("flight", "ground_control"):
                    rows.append(
                        {
                            "id.accession": accession,
                            "tissue_canonical": tissue,
                            "condition_inferred": condition,
                            "profile_id": f"{accession}-{tissue}-{condition}",
                        }
                    )
        metadata = pd.DataFrame(rows)
        plan = build_pooled_plan(
            metadata,
            seed=5,
            validation_fraction=0.2,
            test_fraction=0.2,
        )
        training = set(plan.loc[plan["role"].eq("training"), "id.accession"])
        retained = metadata.loc[metadata["id.accession"].isin(training)]
        observed = set(
            zip(retained["tissue_canonical"], retained["condition_inferred"])
        )
        expected = set(
            zip(metadata["tissue_canonical"], metadata["condition_inferred"])
        )
        self.assertEqual(observed, expected)


if __name__ == "__main__":
    unittest.main()
