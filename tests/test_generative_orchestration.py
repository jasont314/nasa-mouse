import argparse
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd
import torch

from nasa_mouse_generative.config import BenchmarkConfig, load_config
from nasa_mouse_generative.matrix_runner import (
    _execution_mask,
    _expand_tissues,
    _initial_status,
    _resolve_row_ids,
    config_for_row,
)
from nasa_mouse_generative.metrics import (
    _score_classifier,
    _synthetic_training_profiles,
    classifier_utility,
)
from nasa_mouse_generative.paper_contracts import verify_pinned_source
from nasa_mouse_generative.preprocessing import FittedPreprocessor
from nasa_mouse_generative.profiles import load_model_parameters
from nasa_mouse_generative.runner import run as run_training
from nasa_mouse_generative.training_data import DataPartition, effective_covariates
from nasa_mouse_generative.adapters.wgan import WGANAdapter
from nasa_mouse_wgan.model import ConditionalWGANGP, embedding_dim


def _matrix_row(**updates):
    row = {
        "model": "vinas_wgan_gp",
        "model_profile": "practical_screen",
        "task": "conditional_generation",
        "preprocessing_profile": "shared_log1p_cpm_zscore",
        "feature_space": "hvg_2000",
        "harmonization": "none",
        "training_regime": "osdr_only",
        "accession_scope": "all_eligible",
        "tissue_mode": "pooled_conditioned",
        "condition_on_flight": True,
        "study_policy": "not_conditioned",
        "seed": 2020,
    }
    row.update(updates)
    return row


class PaperContractTests(unittest.TestCase):
    def test_wgan_accelerated_gamma_matches_released_numpy_definition(self):
        rng = np.random.default_rng(9)
        real = rng.normal(size=(30, 8))
        fake = 0.7 * real + rng.normal(scale=0.4, size=real.shape)
        upper = np.triu_indices(real.shape[1], k=1)
        expected = np.corrcoef(
            1.0 - np.corrcoef(real, rowvar=False)[upper],
            1.0 - np.corrcoef(fake, rowvar=False)[upper],
        )[0, 1]
        observed = WGANAdapter._gamma_coefficient(
            real, fake, device=torch.device("cpu")
        )
        self.assertAlmostEqual(observed, expected, places=12)

    def test_wgan_uses_uncapped_official_embedding_rule(self):
        self.assertEqual(embedding_dim(10_000), 101)
        self.assertEqual(embedding_dim(1), 2)

    def test_wgan_paper_topology_includes_released_numeric_placeholder(self):
        model = ConditionalWGANGP(
            expression_dim=10,
            categorical_cardinalities=[4, 9],
            noise_dim=64,
            numeric_dim=1,
            hidden_dims=(256, 256),
        )
        first = model.generator.network[0]
        self.assertEqual(first.in_features, 64 + 1 + 3 + 4)
        self.assertEqual(first.out_features, 256)
        self.assertEqual(model.generator.network[-1].out_features, 10)

    def test_wgan_paper_rmsprop_and_stopping_variants_match_source_audit(self):
        base = load_config("configs/generative/default.yaml")
        parameters = load_model_parameters(
            replace(
                base,
                training=replace(base.training, model_profile="paper_native"),
            )
        )
        with tempfile.TemporaryDirectory() as directory:
            adapter = WGANAdapter(
                genes=["g1", "g2"],
                cardinalities=[3, 2],
                covariates=("tissue", "data_source"),
                parameters=parameters,
                device_spec="cpu",
                output_dir=Path(directory),
                checkpoint_every=100,
                resume=False,
                seed=1,
                num_workers=0,
                source_path="assets/model_sources/adversarial-gene-expression",
            )
            generator_optimizer, _ = adapter._optimizers(5e-4)
            self.assertEqual(generator_optimizer.defaults["alpha"], 0.9)
            self.assertEqual(generator_optimizer.defaults["eps"], 1e-7)
            self.assertEqual(adapter._early_stopping_checks(), 10)
            self.assertEqual(
                [epoch for epoch in range(1, 13) if adapter._is_monitor_epoch(epoch)],
                [1, 6, 11],
            )
            adapter.parameters.update(
                {
                    "early_stopping_variant": "paper_text",
                    "early_stopping_patience_epochs": 30,
                    "early_stopping_evaluate_every_epochs": 5,
                }
            )
            self.assertEqual(adapter._early_stopping_checks(), 6)

    def test_paper_native_wgan_lr_cannot_be_silently_changed(self):
        base = load_config("configs/generative/default.yaml")
        changed = replace(
            base,
            training=replace(
                base.training,
                model_profile="paper_native",
                model_parameters={"learning_rate": 5e-5},
            ),
        )
        with self.assertRaisesRegex(ValueError, "paper-native contract"):
            load_model_parameters(changed)

    def test_paper_native_wgan_stopping_cannot_be_disabled(self):
        base = load_config("configs/generative/default.yaml")
        changed = replace(
            base,
            training=replace(
                base.training,
                model_profile="paper_native",
                model_parameters={"early_stopping": False},
            ),
        )
        with self.assertRaisesRegex(ValueError, "paper-native contract"):
            load_model_parameters(changed)

    def test_all_pinned_source_hashes_verify(self):
        roots = {
            "vinas_wgan_gp": "assets/model_sources/adversarial-gene-expression",
            "lacan_diffusion": "assets/model_sources/rna-diffusion",
            "genejepa": "assets/model_sources/GeneJEPA",
        }
        for model, root in roots.items():
            with self.subTest(model=model):
                manifest = verify_pinned_source(model, root)
                self.assertTrue(manifest["source_file_sha256"])


class MatrixResolutionTests(unittest.TestCase):
    def test_pooled_matrix_rows_receive_stable_resume_ids(self):
        table = pd.DataFrame([_matrix_row(tissue_mode="pooled_conditioned")])
        first = _expand_tissues(
            table, inventory_path=Path("unused.tsv"), tissue_filter=""
        )
        second = _expand_tissues(
            table, inventory_path=Path("unused.tsv"), tissue_filter=""
        )
        self.assertEqual(first.loc[0, "row_id"], second.loc[0, "row_id"])

    def test_phase_filter_does_not_remove_other_rows_from_ledger(self):
        plan = pd.DataFrame(
            [
                {**_matrix_row(), "phase": "first", "row_id": "current-a"},
                {**_matrix_row(), "phase": "second", "row_id": "current-b"},
            ]
        )
        existing = pd.DataFrame(
            [
                {
                    **_matrix_row(),
                    "phase": "retired",
                    "row_id": "old-row",
                    "status": "complete",
                    "run_summary": "old.json",
                },
                {
                    **_matrix_row(),
                    "phase": "retired",
                    "row_id": "never-started",
                    "status": "planned",
                    "run_summary": "",
                },
            ]
        )
        status = _initial_status(plan, existing)
        selected = _execution_mask(status, phases=["second"], tissue_filter="")
        self.assertEqual(set(status["row_id"]), {"current-a", "current-b", "old-row"})
        self.assertEqual(status.loc[selected, "row_id"].tolist(), ["current-b"])
        self.assertFalse(
            bool(status.loc[status["row_id"].eq("old-row"), "in_current_plan"].iloc[0])
        )

    def test_matrix_identity_uses_resolved_parameters(self):
        table = _expand_tissues(
            pd.DataFrame([{**_matrix_row(), "phase": "screen"}]),
            inventory_path=Path("unused.tsv"),
            tissue_filter="",
        )
        base = load_config("configs/generative/default.yaml")
        changed = replace(
            base,
            training=replace(
                base.training, model_parameters={"learning_rate": 0.000321}
            ),
        )
        first = _resolve_row_ids(table, base).loc[0, "row_id"]
        second = _resolve_row_ids(table, changed).loc[0, "row_id"]
        self.assertNotEqual(first, second)

    def test_study_conditioning_is_an_explicit_alternative(self):
        config = config_for_row(
            load_config("configs/generative/default.yaml"),
            _matrix_row(study_policy="conditioned"),
        )
        self.assertIn("study", config.training.conditioning_covariates)
        self.assertEqual(config.training.study_policy, "conditioned")

    def test_matrix_conditioning_profile_can_reduce_sparse_covariates(self):
        config = config_for_row(
            load_config("configs/generative/default.yaml"),
            _matrix_row(conditioning_profile="condition_tissue"),
        )
        self.assertEqual(
            config.training.conditioning_covariates,
            ("condition", "tissue"),
        )

    def test_matrix_rejects_unknown_conditioning_profile(self):
        with self.assertRaisesRegex(ValueError, "Unknown conditioning_profile"):
            config_for_row(
                load_config("configs/generative/default.yaml"),
                _matrix_row(conditioning_profile="not_a_profile"),
            )

    def test_unconditional_control_removes_condition_at_runtime(self):
        config = config_for_row(
            load_config("configs/generative/default.yaml"),
            _matrix_row(condition_on_flight=False),
        )
        self.assertNotIn("condition", effective_covariates(config))

    def test_single_and_selected_scopes_require_explicit_accessions(self):
        base = load_config("configs/generative/default.yaml")
        with self.assertRaisesRegex(ValueError, "awaiting_accession_selection"):
            config_for_row(base, _matrix_row(accession_scope="single"))
        selected_base = replace(
            base,
            data=replace(
                base.data,
                osdr_include_accessions=("OSD-1", "OSD-2"),
            ),
        )
        selected = config_for_row(
            selected_base, _matrix_row(accession_scope="selected")
        )
        self.assertEqual(selected.data.osdr_accession_scope, "selected")

    def test_every_harmonizer_has_a_valid_configured_arm(self):
        base = load_config("configs/generative/default.yaml")
        for method in (
            "none",
            "within_study_zscore",
            "within_study_then_global_zscore",
            "combat",
            "combat_seq",
            "mober",
        ):
            row = _matrix_row(
                harmonization=method,
                preprocessing_profile=(
                    "shared_raw" if method == "combat_seq" else "shared_log1p_cpm_zscore"
                ),
            )
            with self.subTest(method=method):
                config = config_for_row(base, row)
                self.assertEqual(config.preprocessing.harmonization, method)
                if method in {"combat", "combat_seq"}:
                    self.assertTrue(
                        config.validation.allow_transductive_preprocessing
                    )


class GenerationParameterTests(unittest.TestCase):
    def test_balanced_accuracy_averages_only_heldout_classes(self):
        class FixedModel:
            def predict(self, matrix):
                return np.asarray(["liver", "training_only", "skin", "skin"])

        result = _score_classifier(
            FixedModel(),
            np.zeros((4, 2)),
            np.asarray(["liver", "liver", "skin", "skin"]),
        )
        self.assertAlmostEqual(result["balanced_accuracy"], 0.75)

    def test_ratio_and_per_profile_cap_control_generated_count(self):
        categories = np.asarray([[0], [0], [1], [1], [2], [2], [3], [3], [0], [1]])
        labels = np.asarray(["flight"] * 4 + ["ground_control"] * 4 + ["flight"] * 2)
        partition = DataPartition(
            name="train",
            matrix=np.zeros((10, 3), dtype=np.float32),
            obs=pd.DataFrame({"condition": labels}),
            categories=categories,
            weights=np.full(10, 0.1, dtype=np.float32),
        )
        sampled, sampled_labels, audit = _synthetic_training_profiles(
            partition,
            ratio=2.0,
            samples_per_covariate_profile=3,
            max_samples=100,
            seed=4,
        )
        self.assertEqual(len(sampled), 12)
        self.assertEqual(len(sampled_labels), 12)
        self.assertTrue(audit["capacity_limited"])

    def test_augmentation_is_not_scored_before_quality_gates_pass(self):
        rng = np.random.default_rng(2)
        train = rng.normal(size=(20, 5))
        evaluation = rng.normal(size=(10, 5))
        train_labels = np.asarray(["flight"] * 10 + ["ground_control"] * 10)
        evaluation_labels = np.asarray(["flight"] * 5 + ["ground_control"] * 5)
        result = classifier_utility(
            train,
            train_labels,
            evaluation,
            evaluation_labels,
            synthetic_train=train.copy(),
            synthetic_labels=train_labels,
            allow_augmentation=False,
        )
        self.assertEqual(
            result["augmentation_status"], "blocked_by_generator_quality_gates"
        )
        self.assertNotIn("real_plus_synthetic_train_real_evaluation", result)


class RepeatExecutionTests(unittest.TestCase):
    def test_repeat_count_expands_to_deterministic_seed_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            args = argparse.Namespace(
                config="configs/generative/default.yaml",
                set=[f"output_root={directory}", "training.repeats=3"],
                tissue="",
                all_tissues=False,
                run_name="repeat_test",
                smoke=False,
            )
            observed = []

            def fake_train(config, **kwargs):
                observed.append((config.training.seed, kwargs["run_name"]))
                return Path(directory) / f"seed_{config.training.seed}.json"

            with mock.patch(
                "nasa_mouse_generative.runner.train_one", side_effect=fake_train
            ):
                summary = run_training(args)
        self.assertEqual([seed for seed, _ in observed], [2020, 2021, 2022])
        self.assertTrue(all(name.endswith(str(seed)) for seed, name in observed))
        self.assertEqual(summary.name, "repeat_batch_summary.json")


class PreprocessingFamilyTests(unittest.TestCase):
    def test_transform_and_scaler_families_are_finite(self):
        matrix = np.asarray(
            [[1.0, 2.0, 7.0], [4.0, 1.0, 5.0], [2.0, 8.0, 1.0]],
            dtype=np.float32,
        )
        studies = ["A", "A", "B"]
        for transform in ("none", "log1p", "log2p1"):
            for scaler in (
                "none",
                "zscore",
                "global_zscore",
                "nonzero_global_zscore",
                "robust",
                "maxabs",
            ):
                with self.subTest(transform=transform, scaler=scaler):
                    spec = replace(
                        BenchmarkConfig().preprocessing,
                        library_normalization="none",
                        transform=transform,
                        scaler=scaler,
                    )
                    values = FittedPreprocessor(spec).fit_transform(matrix, studies)
                    self.assertTrue(np.isfinite(values).all())


if __name__ == "__main__":
    unittest.main()
