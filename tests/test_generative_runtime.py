import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from nasa_mouse_generative.conditioning import CategoryEncoder, UNKNOWN
from nasa_mouse_generative.config import (
    BenchmarkConfig,
    DataConfig,
    ExecutionConfig,
    PreprocessingConfig,
    TrainingConfig,
    load_config_with_overrides,
)
from nasa_mouse_generative.adapters import load_adapter
from nasa_mouse_generative.adapters.diffusion import DiffusionAdapter
from nasa_mouse_generative.adapters.wgan import WGANAdapter
from nasa_mouse_generative.generate import _default_profile
from nasa_mouse_generative.metrics import fidelity_selection
from nasa_mouse_generative.preprocessing import FittedPreprocessor
from nasa_mouse_generative.profiles import resolve_preprocessing_profile
from nasa_mouse_generative.runner import _claim_run_identity
from nasa_mouse_generative.training_data import (
    DataPartition,
    _single_accession_roles,
    extract_archs4_matrix,
)


class RuntimeConfigTests(unittest.TestCase):
    def test_dotted_overrides_resolve_model_parameters(self):
        config = load_config_with_overrides(
            "configs/generative/default.yaml",
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

    def test_shared_preprocessing_profile_is_resolved(self):
        config = load_config_with_overrides(
            "configs/generative/default.yaml",
            ["preprocessing.profile=shared_log1p_cpm_maxabs"],
        )
        resolved = resolve_preprocessing_profile(config)
        self.assertEqual(resolved.preprocessing.library_normalization, "cpm")
        self.assertEqual(resolved.preprocessing.transform, "log1p")
        self.assertEqual(resolved.preprocessing.scaler, "maxabs")

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
    def test_unseen_categories_use_explicit_unknown_code(self):
        encoder = CategoryEncoder.fit(
            [pd.DataFrame({"condition": ["flight", "ground_control"]})],
            ["condition"],
        )
        observed = encoder.transform(pd.DataFrame({"condition": ["new_condition"]}))
        self.assertEqual(
            int(observed[0, 0]), encoder.vocabularies["condition"].index(UNKNOWN)
        )

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


class Archs4ExtractionTests(unittest.TestCase):
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
        self.assertEqual(generated.shape, (3, 16))
        self.assertTrue(np.isfinite(generated).all())


if __name__ == "__main__":
    unittest.main()
