import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
import sys
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd

from nasa_mouse_generative.conditioning import (
    CategoryEncoder,
    UNKNOWN,
    canonicalize_sex,
)
from nasa_mouse_generative.config import (
    BenchmarkConfig,
    DataConfig,
    ExecutionConfig,
    FeatureConfig,
    PreprocessingConfig,
    TrainingConfig,
    load_config_with_overrides,
)
from nasa_mouse_generative.adapters import load_adapter
from nasa_mouse_generative.adapters.diffusion import DiffusionAdapter
from nasa_mouse_generative.adapters.wgan import WGANAdapter
from nasa_mouse_generative.generate import _default_profile
from nasa_mouse_generative.evaluate import (
    _assert_preprocessing_matches,
    _restore_saved_conditioning,
)
from nasa_mouse_generative.harmonizers import CombatHarmonizer, CombatSeqHarmonizer
from nasa_mouse_generative.metrics import (
    accession_effect_selection,
    conditional_effect_selection,
    fidelity_selection,
)
from nasa_mouse_generative.paper_metrics import paper_distribution_metrics
from nasa_mouse_generative.preprocessing import FittedPreprocessor, ScaleStats
from nasa_mouse_generative.profiles import resolve_preprocessing_profile
from nasa_mouse_generative.runner import _claim_run_identity
from nasa_mouse_generative.scoreboard import _per_tissue_diagnostics, _unified_row
from nasa_mouse_generative.training_data import (
    DataPartition,
    _retain_readable_archs4_metadata,
    _split_archs4_selection,
    _single_accession_roles,
    extract_archs4_matrix,
    prepare_training_data,
)


class RuntimeConfigTests(unittest.TestCase):
    def test_paper_metrics_only_subsample_adversarial_accuracy(self):
        rng = np.random.default_rng(31)
        real = rng.normal(size=(30, 8))
        synthetic = real + rng.normal(scale=0.05, size=real.shape)
        metrics = paper_distribution_metrics(
            real,
            synthetic,
            max_samples=30,
            neighbors=3,
            adversarial_max_samples=10,
            seed=4,
        )
        self.assertEqual(metrics["metric_samples"], 30)
        self.assertEqual(metrics["adversarial_metric_samples"], 10)

    def test_accession_effect_gate_requires_replicated_effect_recovery(self):
        passing = accession_effect_selection(
            {
                "accessions": 4,
                "meta_effect_correlation": 0.7,
                "meta_direction_agreement": 0.8,
            }
        )
        confounded = accession_effect_selection(
            {
                "accessions": 12,
                "meta_effect_correlation": -0.02,
                "meta_direction_agreement": 0.55,
            }
        )
        self.assertTrue(passing["passed"])
        self.assertFalse(confounded["passed"])

    def test_scoreboard_counts_independently_passing_tissues(self):
        result = _per_tissue_diagnostics(
            [
                {
                    "all_fidelity_metrics_pass": True,
                    "flt_gc_delta_correlation": 0.4,
                    "flt_gc_direction_agreement": 0.6,
                },
                {
                    "all_fidelity_metrics_pass": False,
                    "flt_gc_delta_correlation": 0.1,
                    "flt_gc_direction_agreement": 0.7,
                },
            ]
        )
        self.assertEqual(result["per_tissue_fidelity_evaluable"], 2)
        self.assertEqual(result["per_tissue_fidelity_passes"], 1)
        self.assertEqual(result["per_tissue_condition_passes"], 1)

    def test_scoreboard_uses_run_namespace_for_model_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            summary = (
                Path(temporary)
                / "runs"
                / "lacan_diffusion"
                / "run-x"
                / "run_summary.json"
            )
            summary.parent.mkdir(parents=True)
            summary.write_text(
                '{"run_id": "run-x", "model": "/tmp/model.pt"}\n',
                encoding="utf-8",
            )
            row = _unified_row(summary)
        self.assertEqual(row["model"], "lacan_diffusion")
        self.assertEqual(row["implementation"], "")

    def test_dotted_overrides_resolve_model_parameters(self):
        config = load_config_with_overrides(
            "configs/generative/benchmark/default.yaml",
            [
                "training.regime=osdr_only",
                "training.model_parameters.epochs=3",
                "features.max_genes=64",
            ],
        )
        self.assertEqual(config.training.regime, "osdr_only")
        self.assertEqual(config.training.model_parameters["epochs"], 3)
        self.assertEqual(config.features.max_genes, 64)

    def test_archs4_only_cannot_claim_flight_conditioning(self):
        config = BenchmarkConfig(
            training=TrainingConfig(
                regime="archs4_only", condition_on_flight=True
            )
        )
        with self.assertRaisesRegex(ValueError, "ARCHS4 has no flight"):
            config.validate()

    def test_combat_requires_explicit_transductive_permission(self):
        config = BenchmarkConfig(
            preprocessing=PreprocessingConfig(harmonization="combat")
        )
        with self.assertRaisesRegex(ValueError, "transductive sensitivity"):
            config.validate()

    def test_unconditional_combat_cannot_preserve_condition(self):
        config = BenchmarkConfig(
            preprocessing=PreprocessingConfig(harmonization="combat"),
            training=replace(TrainingConfig(), condition_on_flight=False),
            validation=replace(
                BenchmarkConfig().validation,
                allow_transductive_preprocessing=True,
            ),
        )
        with self.assertRaisesRegex(ValueError, "unconditional negative-control"):
            config.validate()

    def test_shared_preprocessing_profile_is_resolved(self):
        config = load_config_with_overrides(
            "configs/generative/benchmark/default.yaml",
            ["preprocessing.profile=shared_log1p_cpm_maxabs"],
        )
        resolved = resolve_preprocessing_profile(config)
        self.assertEqual(resolved.preprocessing.library_normalization, "cpm")
        self.assertEqual(resolved.preprocessing.transform, "log1p")
        self.assertEqual(resolved.preprocessing.scaler, "maxabs")

    def test_wgan_native_and_nasa_cpm_profiles_are_distinct(self):
        native = load_config_with_overrides(
            "configs/generative/benchmark/default.yaml",
            ["preprocessing.profile=model_native"],
        )
        nasa = load_config_with_overrides(
            "configs/generative/benchmark/default.yaml",
            ["preprocessing.profile=wgan_nasa_cpm_zscore"],
        )
        native = resolve_preprocessing_profile(native)
        nasa = resolve_preprocessing_profile(nasa)
        self.assertEqual(native.preprocessing.library_normalization, "none")
        self.assertEqual(nasa.preprocessing.library_normalization, "cpm")
        self.assertEqual(native.preprocessing.transform, "log1p")
        self.assertEqual(native.preprocessing.scaler, "zscore")

    def test_run_identity_checks_legacy_summary_before_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "run_summary.json").write_text(
                '{"run_sha256": "old-digest"}\n', encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "already belongs"):
                _claim_run_identity(
                    root,
                    identifier="named-run",
                    digest="new-digest",
                    model="vinas_wgan_gp",
                )
            self.assertFalse((root / "run_identity.json").exists())


class ConditioningTests(unittest.TestCase):
    def test_sex_labels_are_canonicalized_without_strain_false_positives(self):
        self.assertEqual(canonicalize_sex("M"), "male")
        self.assertEqual(canonicalize_sex("virgin female"), "female")
        self.assertEqual(canonicalize_sex("12 weeks/female"), "female")
        self.assertEqual(canonicalize_sex("pooled male and female"), "mixed")
        self.assertEqual(canonicalize_sex("3 males and 4 females"), "mixed")
        self.assertEqual(canonicalize_sex("C57BL/6J"), "unknown_sex")

    def test_unseen_categories_use_explicit_unknown_code(self):
        encoder = CategoryEncoder.fit(
            [pd.DataFrame({"condition": ["flight", "ground_control"]})],
            ["condition"],
        )
        observed = encoder.transform(pd.DataFrame({"condition": ["new_condition"]}))
        self.assertEqual(
            int(observed[0, 0]), encoder.vocabularies["condition"].index(UNKNOWN)
        )

    def test_reconstructed_evaluation_restores_saved_encoder_and_order(self):
        saved_obs = pd.DataFrame(
            {
                "profile_id": ["p1", "p2"],
                "accession": ["OSD-1", "OSD-2"],
                "condition": ["flight", "ground_control"],
            }
        )
        current_obs = saved_obs.iloc[::-1].reset_index(drop=True)
        encoder = CategoryEncoder.fit([saved_obs], ["condition"])
        partition = DataPartition(
            name="train",
            matrix=np.asarray([[2.0], [1.0]], dtype=np.float32),
            obs=current_obs,
            categories=encoder.transform(current_obs),
            weights=np.asarray([0.6, 0.4], dtype=np.float32),
        )
        adapter = SimpleNamespace(
            covariates=encoder.covariates,
            cardinalities=encoder.cardinalities,
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            encoder.save(root / "categorical_encoder.json")
            saved_obs.to_csv(root / "train_obs.tsv.gz", sep="\t", index=False)
            restored = _restore_saved_conditioning(
                {"train": partition}, root, adapter
            )["train"]
        np.testing.assert_array_equal(restored.matrix[:, 0], [1.0, 2.0])
        self.assertEqual(restored.obs["profile_id"].tolist(), ["p1", "p2"])
        np.testing.assert_array_equal(
            restored.categories, encoder.transform(saved_obs)
        )

    def test_reconstructed_evaluation_rejects_changed_preprocessing_fit(self):
        first = FittedPreprocessor(PreprocessingConfig())
        second = FittedPreprocessor(PreprocessingConfig())
        first.final_stats = ScaleStats(
            center=np.asarray([0.0]), scale=np.asarray([1.0])
        )
        second.final_stats = ScaleStats(
            center=np.asarray([0.1]), scale=np.asarray([1.0])
        )
        with self.assertRaisesRegex(ValueError, "final_stats differs"):
            _assert_preprocessing_matches(first, second)

    def test_single_accession_fallback_preserves_each_condition_in_training(self):
        frame = pd.DataFrame(
            {
                "profile_id": [f"p{index}" for index in range(20)],
                "tissue": ["liver"] * 20,
                "condition": ["flight"] * 10 + ["ground_control"] * 10,
            }
        )
        roles = _single_accession_roles(
            frame, seed=2, validation_fraction=0.2, test_fraction=0.2
        )
        for condition in ("flight", "ground_control"):
            selected = frame.loc[roles.eq("train"), "condition"]
            self.assertIn(condition, set(selected))
        self.assertIn("validation", set(roles))
        self.assertIn("test", set(roles))

    def test_generation_default_uses_observed_joint_profile(self):
        table = pd.DataFrame(
            {
                "condition": ["flight", "ground_control"],
                "tissue": ["skin", "liver"],
                "material_type": ["dorsal skin", "Liver"],
                "study": ["OSD-1", "OSD-2"],
            }
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table.to_csv(root / "train_obs.tsv.gz", sep="\t", index=False)
            profile, matched = _default_profile(
                root,
                ("condition", "tissue", "material_type"),
                {"tissue": "skin"},
            )
        self.assertTrue(matched)
        self.assertEqual(profile["material_type"], "dorsal skin")
        self.assertEqual(profile["study"], "OSD-1")

    def test_generation_default_reports_unobserved_constraint_combination(self):
        table = pd.DataFrame(
            {
                "condition": ["flight", "ground_control"],
                "tissue": ["skin", "liver"],
                "study": ["OSD-1", "OSD-2"],
            }
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table.to_csv(root / "train_obs.tsv.gz", sep="\t", index=False)
            _, matched = _default_profile(
                root,
                ("condition", "tissue"),
                {"condition": "flight", "tissue": "liver"},
            )
        self.assertFalse(matched)


class PreprocessorSerializationTests(unittest.TestCase):
    def test_round_trip_and_inverse_transform(self):
        processor = FittedPreprocessor(
            PreprocessingConfig(
                input_units="raw_counts",
                library_normalization="none",
                transform="log1p",
                scaler="zscore",
                harmonization="none",
            )
        )
        matrix = np.asarray([[1.0, 2.0], [3.0, 8.0], [5.0, 4.0]])
        transformed = processor.fit_transform(matrix, ["A", "A", "B"])
        with tempfile.TemporaryDirectory() as directory:
            processor.save(directory)
            restored = FittedPreprocessor.load(directory)
            second = restored.transform(matrix, ["A", "A", "B"])
            inverted = restored.inverse_transform(second, ["A", "A", "B"])
        np.testing.assert_allclose(second, transformed)
        np.testing.assert_allclose(inverted, matrix, atol=1e-5)

    def test_inverse_cpm_is_compositional(self):
        processor = FittedPreprocessor(
            PreprocessingConfig(
                input_units="raw_counts",
                library_normalization="cpm",
                transform="log1p",
                scaler="zscore",
                harmonization="none",
            )
        )
        matrix = np.asarray([[1.0, 2.0, 7.0], [4.0, 1.0, 5.0]])
        transformed = processor.fit_transform(matrix, ["A", "B"])
        inverted = processor.inverse_transform(transformed, ["A", "B"])
        np.testing.assert_allclose(inverted.sum(axis=1), 1_000_000.0, rtol=1e-6)

    def test_nonzero_global_scaler_preserves_absent_gene_tokens(self):
        processor = FittedPreprocessor(
            PreprocessingConfig(
                input_units="raw_counts",
                library_normalization="none",
                transform="none",
                scaler="nonzero_global_zscore",
                harmonization="none",
            )
        )
        matrix = np.asarray([[0.0, 2.0, 4.0], [8.0, 0.0, 0.0]])
        transformed = processor.fit_transform(matrix, ["A", "B"])
        np.testing.assert_array_equal(transformed[matrix == 0], 0.0)
        np.testing.assert_allclose(transformed[matrix != 0].mean(), 0.0, atol=1e-6)
        inverted = processor.inverse_transform(transformed, ["A", "B"])
        np.testing.assert_allclose(inverted, matrix, atol=1e-6)

    @staticmethod
    def _harmonization_data():
        rng = np.random.default_rng(12)
        matrix = rng.normal(5.0, 0.5, size=(16, 8)).astype(np.float32)
        studies = np.asarray(["A"] * 8 + ["B"] * 8)
        matrix[8:] += 3.0
        metadata = pd.DataFrame(
            {
                "study": studies,
                "condition": ["flight", "ground_control"] * 8,
                "tissue": ["liver"] * 16,
                "sex": ["female"] * 16,
                "source": ["osdr"] * 16,
            }
        )
        return matrix, studies, metadata

    def test_combat_corrects_batch_location_and_round_trips(self):
        matrix, studies, metadata = self._harmonization_data()
        processor = FittedPreprocessor(
            PreprocessingConfig(
                input_units="cpm",
                library_normalization="none",
                transform="none",
                scaler="none",
                harmonization="combat",
                harmonization_parameters={"batch_key": "study"},
            ),
            device_spec="cpu",
            seed=3,
        )
        corrected = processor.fit_transform(
            matrix, studies, metadata=metadata
        )
        self.assertLess(abs(float(corrected[:8].mean() - corrected[8:].mean())), 0.2)
        with tempfile.TemporaryDirectory() as directory:
            processor.save(directory)
            restored = FittedPreprocessor.load(directory)
            second = restored.transform(
                matrix, studies, metadata=metadata
            )
        np.testing.assert_allclose(second, corrected, atol=1e-5)

    def test_mober_projects_with_frozen_serialized_model(self):
        matrix, studies, metadata = self._harmonization_data()
        processor = FittedPreprocessor(
            PreprocessingConfig(
                input_units="cpm",
                library_normalization="none",
                transform="none",
                scaler="zscore",
                harmonization="mober",
                harmonization_parameters={
                    "batch_key": "study",
                    "epochs": 1,
                    "batch_size": 4,
                    "encoding_dim": 4,
                    "projection_batch_size": 8,
                },
            ),
            device_spec="cpu",
            seed=3,
        )
        corrected = processor.fit_transform(
            matrix, studies, metadata=metadata
        )
        with tempfile.TemporaryDirectory() as directory:
            processor.save(directory)
            restored = FittedPreprocessor.load(directory)
            projected = restored.transform(
                matrix[:4], studies[:4], metadata=metadata.iloc[:4]
            )
        self.assertEqual(corrected.shape, matrix.shape)
        self.assertEqual(projected.shape, (4, matrix.shape[1]))
        self.assertTrue(np.isfinite(projected).all())


class CombatSeqPolicyTests(unittest.TestCase):
    def test_fractional_counts_require_explicit_rounding_policy(self):
        values = np.asarray([[1.2, 3.0], [2.0, 4.7]])
        strict = CombatSeqHarmonizer(
            covariates=("condition",),
            parameters={},
            device_spec="cpu",
            seed=1,
        )
        with self.assertRaisesRegex(ValueError, "fractional count-like"):
            strict._integer_counts(values)
        rounded = CombatSeqHarmonizer(
            covariates=("condition",),
            parameters={"noninteger_policy": "round"},
            device_spec="cpu",
            seed=1,
        )
        observed = rounded._integer_counts(values)
        np.testing.assert_array_equal(observed, [[1, 3], [2, 5]])
        self.assertEqual(rounded.rounding_audit[0]["entries_rounded"], 2)

    def test_combat_confounded_covariate_is_strict_by_default(self):
        values = np.arange(24, dtype=np.float32).reshape(6, 4)
        studies = ["A"] * 3 + ["B"] * 3
        metadata = pd.DataFrame(
            {
                "study": studies,
                "condition": ["flight"] * 3 + ["ground_control"] * 3,
            }
        )
        strict = CombatHarmonizer(
            covariates=("condition",),
            parameters={"batch_key": "study"},
            device_spec="cpu",
            seed=1,
        )
        with self.assertRaisesRegex(ValueError, "confounded with batch"):
            strict.fit_transform(values, studies, metadata)

        permissive = CombatHarmonizer(
            covariates=("condition",),
            parameters={
                "batch_key": "study",
                "confounded_covariate_policy": "drop",
            },
            device_spec="cpu",
            seed=1,
        )
        corrected = permissive.fit_transform(values, studies, metadata)
        self.assertEqual(corrected.shape, values.shape)
        self.assertEqual(permissive.audit()["dropped_covariates"], ["condition"])


class MBatchAdapterTests(unittest.TestCase):
    @unittest.skipUnless(
        Path("assets/model_sources/MBatch/apps/MBatch/R").exists()
        and Path(sys.executable).with_name("Rscript").exists(),
        "optional pinned MBatch source and conda Rscript are required",
    )
    def test_official_mbatch_methods_round_trip_heldout_batches(self):
        rng = np.random.default_rng(7)
        train_studies = np.repeat(["A", "B", "C"], 8)
        train = np.maximum(
            rng.normal(5.0, 1.0, (24, 12))
            + np.repeat([0.0, 2.0, -1.0], 8)[:, None],
            0.0,
        ).astype(np.float32)
        heldout_studies = np.repeat(["D", "E"], 3)
        heldout = np.maximum(
            rng.normal(5.0, 1.0, (6, 12))
            + np.repeat([3.0, -2.0], 3)[:, None],
            0.0,
        ).astype(np.float32)
        train_metadata = pd.DataFrame({"study": train_studies})
        heldout_metadata = pd.DataFrame({"study": heldout_studies})
        for method in (
            "mbatch_median_polish",
            "mbatch_empirical_bayes",
            "mbatch_anova",
        ):
            with self.subTest(method=method):
                processor = FittedPreprocessor(
                    PreprocessingConfig(
                        input_units="normalized_counts",
                        library_normalization="none",
                        transform="none",
                        scaler="none",
                        harmonization=method,
                        harmonization_covariates=(),
                        harmonization_parameters={
                            "batch_key": "study",
                            "anchor_samples": 12,
                            "nonfinite_policy": "identity_gene",
                        },
                        unseen_study_policy="transductive_unlabeled",
                    ),
                    device_spec="cpu",
                    seed=7,
                )
                corrected = processor.fit_transform(
                    train, train_studies, metadata=train_metadata
                )
                projected = processor.transform(
                    heldout,
                    heldout_studies,
                    metadata=heldout_metadata,
                    allow_transductive=True,
                )
                with tempfile.TemporaryDirectory() as directory:
                    processor.save(directory)
                    restored = FittedPreprocessor.load(directory)
                    replay = restored.transform(
                        heldout,
                        heldout_studies,
                        metadata=heldout_metadata,
                        allow_transductive=True,
                    )
                self.assertEqual(corrected.shape, train.shape)
                self.assertTrue(np.isfinite(projected).all())
                np.testing.assert_allclose(projected, replay, atol=1e-5)


class Archs4ExtractionTests(unittest.TestCase):
    def test_unreadable_archs4_profiles_are_removed_without_zero_imputation(self):
        metadata = pd.DataFrame(
            {
                "archs4_sample_index": [11, 12, 13],
                "geo_accession": ["GSM11", "GSM12", "GSM13"],
            }
        )
        matrix = np.ones((2, 4), dtype=np.float32)
        retained = _retain_readable_archs4_metadata(
            metadata,
            matrix,
            {"skipped_corrupt_sample_indices": [12]},
        )
        self.assertEqual(retained["archs4_sample_index"].tolist(), [11, 13])

    def test_archs4_series_split_never_crosses_roles(self):
        rows = []
        for tissue in ("liver", "kidney"):
            for series_index in range(8):
                for sample_index in range(2):
                    rows.append(
                        {
                            "geo_accession": f"{tissue}-{series_index}-{sample_index}",
                            "series_id": f"{tissue}-GSE{series_index}",
                            "canonical_tissue": tissue,
                        }
                    )
        metadata, audit = _split_archs4_selection(
            pd.DataFrame(rows), BenchmarkConfig()
        )
        self.assertEqual(metadata.groupby("series_id")["role"].nunique().max(), 1)
        self.assertEqual(set(metadata["role"]), {"train", "validation", "test"})
        self.assertEqual(
            set(metadata.loc[metadata["role"].eq("train"), "canonical_tissue"]),
            {"liver", "kidney"},
        )
        self.assertEqual(audit["split_unit"], "GEO series_id")

    def test_archs4_only_preparation_does_not_require_osdr_matrix(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "archs4.h5"
            genes = [f"ENSMUSG{index:011d}" for index in range(6)]
            rng = np.random.default_rng(8)
            with h5py.File(source, "w") as handle:
                handle.create_dataset(
                    "meta/genes/ensembl_gene",
                    data=np.asarray(genes, dtype="S"),
                )
                handle.create_dataset(
                    "data/expression",
                    data=rng.integers(0, 100, size=(6, 24), dtype=np.uint32),
                )
            catalog = root / "catalog"
            catalog.mkdir()
            rows = []
            for sample_index in range(24):
                tissue = "liver" if sample_index < 12 else "kidney"
                rows.append(
                    {
                        "geo_accession": f"GSM{sample_index}",
                        "series_id": f"GSE{sample_index // 2}",
                        "canonical_tissue": tissue,
                        "archs4_sample_index": sample_index,
                        "source_name_ch1": tissue,
                        "characteristics_ch1": "adult wild type control",
                        "library_strategy": "RNA-Seq",
                    }
                )
            pd.DataFrame(rows).to_csv(
                catalog / "archs4_healthy_preferred_balanced.tsv.gz",
                sep="\t",
                index=False,
                compression="gzip",
            )
            config = BenchmarkConfig(
                output_root=str(root / "output"),
                data=replace(
                    DataConfig(),
                    archs4_h5=str(source),
                    archs4_catalog_dir=str(catalog),
                    osdr_h5ad=str(root / "does-not-exist.h5ad"),
                ),
                features=FeatureConfig(space="all_shared", max_genes=4),
                preprocessing=PreprocessingConfig(
                    input_units="raw_counts",
                    library_normalization="none",
                    transform="log1p",
                    scaler="none",
                ),
                training=TrainingConfig(
                    model="lacan_diffusion",
                    task="conditional_generation",
                    regime="archs4_only",
                    condition_on_flight=False,
                    conditioning_covariates=("tissue",),
                ),
                execution=replace(
                    ExecutionConfig(), cache_archs4=True, device="cpu"
                ),
            )
            prepared = prepare_training_data(config)
            cached = prepare_training_data(config)
        self.assertEqual(sum(map(len, prepared.partitions.values())), 24)
        self.assertFalse(prepared.metadata["osdr_expression_used"])
        self.assertEqual(prepared.reference.name, "train")
        self.assertEqual(len(prepared.genes), 4)
        self.assertFalse(prepared.metadata["features"]["selection_cache_hit"])
        self.assertTrue(cached.metadata["features"]["selection_cache_hit"])
        self.assertEqual(prepared.genes, cached.genes)

    def test_selected_columns_and_gene_order_are_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "archs4.h5"
            with h5py.File(source, "w") as handle:
                handle.create_dataset(
                    "meta/genes/ensembl_gene",
                    data=np.asarray([b"g1", b"g2", b"g3"]),
                )
                handle.create_dataset(
                    "data/expression",
                    data=np.asarray([[1, 2], [3, 4], [5, 6]], dtype=np.uint32),
                )
            config = BenchmarkConfig(
                output_root=str(root / "out"),
                data=replace(DataConfig(), archs4_h5=str(source)),
                execution=replace(ExecutionConfig(), cache_archs4=True),
            )
            metadata = pd.DataFrame(
                {
                    "archs4_sample_index": [1, 0],
                    "hierarchical_sampling_weight": [0.5, 0.5],
                }
            )
            matrix, first = extract_archs4_matrix(config, metadata, ["g3", "g1"])
            cached, second = extract_archs4_matrix(config, metadata, ["g3", "g1"])
        np.testing.assert_array_equal(matrix, [[6, 2], [5, 1]])
        np.testing.assert_array_equal(cached, matrix)
        self.assertFalse(first["cache_hit"])
        self.assertTrue(second["cache_hit"])


def _partition() -> DataPartition:
    rng = np.random.default_rng(7)
    matrix = rng.normal(size=(12, 16)).astype(np.float32)
    categories = np.column_stack([np.arange(12) % 3, np.arange(12) % 2])
    obs = pd.DataFrame(
        {
            "profile_id": [f"p{index}" for index in range(12)],
            "accession": ["A"] * 12,
            "study": ["A"] * 12,
            "tissue": ["liver"] * 12,
            "condition": ["flight", "ground_control"] * 6,
        }
    )
    return DataPartition(
        name="train",
        matrix=matrix,
        obs=obs,
        categories=categories.astype(np.int64),
        weights=np.full(12, 1 / 12, dtype=np.float32),
    )


class AdapterTests(unittest.TestCase):
    def test_wgan_checkpoint_resume_and_reload(self):
        partition = _partition()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            kwargs = dict(
                genes=[f"g{index}" for index in range(16)],
                cardinalities=[3, 2],
                covariates=("condition", "tissue"),
                parameters={
                    "hidden_dims": [8],
                    "noise_dim": 4,
                    "batch_size": 6,
                    "critic_steps": 1,
                    "gradient_penalty": 1.0,
                },
                device_spec="cpu",
                output_dir=root,
                checkpoint_every=1,
                seed=3,
            )
            first = WGANAdapter(**kwargs, resume=False)
            first.fit_stage(partition, stage="osdr", epochs=1, learning_rate=1e-4)
            resumed = WGANAdapter(**kwargs, resume=True)
            resumed.fit_stage(partition, stage="osdr", epochs=2, learning_rate=1e-4)
            resumed.save_final()
            loaded = load_adapter(root, device_spec="cpu")
            generated = loaded.generate(partition.categories[:3], seed=9)
        self.assertEqual(resumed.state.completed_epochs["osdr"], 2)
        self.assertEqual(generated.shape, (3, 16))
        self.assertTrue(np.isfinite(generated).all())


class MetricTests(unittest.TestCase):
    def test_finite_sample_corr_gate_does_not_hide_absolute_paper_result(self):
        selection = fidelity_selection(
            {
                "correlation_matrix_agreement": 0.90,
                "correlation_real_bootstrap_p05": 0.88,
                "precision": 0.98,
                "recall": 0.90,
                "f1": 0.94,
                "adversarial_accuracy": 0.50,
                "frechet_ratio_to_real_split_p95": 0.8,
                "real_global_std": 1.0,
                "fake_global_std": 0.95,
            },
            {"fraction_below_training_p01": 0.0},
        )
        self.assertTrue(selection["eligible_for_model_selection"])
        self.assertFalse(selection["meets_absolute_paper_benchmark"])
        self.assertEqual(
            selection["fidelity_gate"]["requirements"][
                "correlation_matrix_agreement"
            ]["minimum"],
            0.88,
        )

    def test_fidelity_selection_requires_every_paper_metric_without_composite(self):
        selection = fidelity_selection(
            {
                "correlation_matrix_agreement": 0.99,
                "precision": 0.98,
                "recall": 0.90,
                "f1": 0.94,
                "adversarial_accuracy": 0.50,
                "frechet_ratio_to_real_split_p95": 0.8,
                "real_global_std": 1.0,
                "fake_global_std": 0.95,
            },
            {"fraction_below_training_p01": 0.0},
        )
        self.assertTrue(selection["eligible_for_model_selection"])
        self.assertNotIn("heldout_fidelity_composite", selection)
        self.assertEqual(selection["fidelity_gate"]["failed_metrics"], [])

        selection["fidelity_gate"] = fidelity_selection(
            {
                "correlation_matrix_agreement": 0.99,
                "precision": 0.98,
                "recall": 0.90,
                "f1": 0.94,
                "adversarial_accuracy": 0.96,
                "frechet_ratio_to_real_split_p95": 0.8,
                "real_global_std": 1.0,
                "fake_global_std": 0.95,
            },
            {"fraction_below_training_p01": 0.0},
        )["fidelity_gate"]
        self.assertIn(
            "adversarial_accuracy", selection["fidelity_gate"]["failed_metrics"]
        )

    def test_fidelity_gates_reject_collapse(self):
        selection = fidelity_selection(
            {
                "gene_mean_correlation": 0.9,
                "gene_std_correlation": 0.1,
                "f1": 0.0,
                "adversarial_accuracy": 1.0,
                "recall": 0.0,
                "real_global_std": 1.0,
                "fake_global_std": 0.01,
            },
            {"fraction_below_training_p01": 0.0},
        )
        self.assertFalse(selection["diversity_gate"]["passed"])
        self.assertFalse(selection["eligible_for_model_selection"])

    def test_condition_effect_gate_requires_correlation_and_direction(self):
        self.assertTrue(
            conditional_effect_selection(
                {"delta_correlation": 0.31, "direction_agreement": 0.56}
            )["passed"]
        )
        self.assertFalse(
            conditional_effect_selection(
                {"delta_correlation": 0.29, "direction_agreement": 0.80}
            )["passed"]
        )

    def test_fidelity_gate_rejects_distinguishable_low_fidelity_samples(self):
        selection = fidelity_selection(
            {
                "gene_mean_correlation": 0.32,
                "gene_std_correlation": 0.56,
                "f1": 0.53,
                "adversarial_accuracy": 0.95,
                "recall": 0.49,
                "real_global_std": 1.0,
                "fake_global_std": 0.84,
            },
            {"fraction_below_training_p01": 0.0},
        )
        self.assertFalse(selection["fidelity_gate"]["passed"])
        self.assertTrue(selection["diversity_gate"]["passed"])
        self.assertTrue(selection["memorization_gate"]["passed"])
        self.assertFalse(selection["eligible_for_model_selection"])


class DiffusionAdapterTests(unittest.TestCase):
    def test_diffusion_save_reload_and_generation(self):
        partition = _partition()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            adapter = DiffusionAdapter(
                genes=[f"g{index}" for index in range(16)],
                cardinalities=[3, 2],
                covariates=("condition", "tissue"),
                parameters={
                    "hidden_dim": 8,
                    "n_blocks": 1,
                    "batch_size": 6,
                    "diffusion_timesteps": 8,
                    "sample_steps": 2,
                    "n_landmarks": 4,
                    "landmark_strategy": "hvg",
                    "use_amp": False,
                    "reconstruction_samples": 12,
                },
                device_spec="cpu",
                output_dir=root,
                checkpoint_every=1,
                resume=False,
                seed=3,
                reconstruction_matrix=partition.matrix,
                l1000_map="",
            )
            adapter.fit_stage(
                partition, stage="osdr", epochs=1, learning_rate=1e-3
            )
            adapter.save_final()
            loaded = load_adapter(root, device_spec="cpu")
            generated = loaded.generate(partition.categories[:3], seed=9)
            trajectory = loaded.generate_trajectory(
                partition.categories[:3],
                seed=9,
                snapshot_timesteps=(8, 2, 0),
                sample_steps=8,
            )
        self.assertEqual(generated.shape, (3, 16))
        self.assertTrue(np.isfinite(generated).all())
        self.assertEqual(set(trajectory), {8, 2, 0})
        self.assertTrue(all(values.shape == (3, 16) for values in trajectory.values()))


if __name__ == "__main__":
    unittest.main()
