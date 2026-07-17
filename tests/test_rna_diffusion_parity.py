from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import anndata as ad
import numpy as np
import pandas as pd
import torch

from nasa_mouse_rna_diffusion.conditional_config import load_conditional_config
from nasa_mouse_rna_diffusion.conditional_data import (
    _deseq2_median_of_ratios,
    _explicit_accession_roles,
    _full_transcriptome_tpm,
    _joint_class_labels,
    _within_study_roles,
)
from nasa_mouse_rna_diffusion.conditional_evaluate import (
    _class_probe,
    _per_tissue_fidelity,
)
from nasa_mouse_rna_diffusion.conditional_train import (
    _expanded_condition_state,
    _scaled_optimizer_step,
)
from nasa_mouse_rna_diffusion.config import load_config
from nasa_mouse_rna_diffusion.data import (
    _group_split_indices,
    _targeted_candidates,
)
from nasa_mouse_rna_diffusion.evaluate import (
    _nearest_neighbor_adversarial_accuracy,
)
from nasa_mouse_rna_diffusion.upstream import (
    ddim_trajectory,
    model_config,
    noise_estimation_loss,
    quadratic_beta_schedule,
    upstream_denoising_module,
    upstream_loss_module,
    upstream_model_class,
    verify_source,
)
from nasa_mouse_rna_diffusion.real_effect_ceiling import (
    analyze_tissue,
    load_development_expression,
)
from nasa_mouse_rna_diffusion.factorized_adapter import (
    FactorizedAdapterDDIM,
    build_factorized_schema,
    encode_factorized_labels,
    neutralize_group,
)
from nasa_mouse_rna_diffusion.factorized_calibrate import CovarianceCalibrator


SMALL_MODEL = {
    "hidden_dims": [16, 16],
    "dropout": 0.1,
    "time_embedding_dim": 1,
    "tissue_embedding_dim": 2,
    "sinusoidal_time": False,
    "diffusion_timesteps": 8,
    "beta_schedule": "quad",
    "beta_start": 0.0001,
    "beta_end": 0.02,
    "ema_decay": 0.999,
}


class UpstreamParityTests(unittest.TestCase):
    def _model(self):
        config = model_config(expression_dim=12, num_classes=3, model=SMALL_MODEL)
        torch.manual_seed(9)
        model = upstream_model_class()(config)
        model.eval()
        return model

    def test_vendored_source_is_pinned(self):
        manifest = verify_source()
        self.assertEqual(
            manifest["source_commit"],
            "cde890154698fcea96c924804aaff04af3351b48",
        )

    def test_noise_loss_matches_upstream_function(self):
        model = self._model()
        clean = torch.randn(4, 12)
        timesteps = torch.tensor([0, 2, 5, 7])
        noise = torch.randn_like(clean)
        labels = torch.nn.functional.one_hot(
            torch.tensor([0, 1, 2, 1]), num_classes=3
        )
        betas = quadratic_beta_schedule(
            beta_start=0.0001, beta_end=0.02, timesteps=8
        )
        observed = noise_estimation_loss(
            model, clean, timesteps, noise, betas, labels
        )
        expected = upstream_loss_module().noise_estimation_loss(
            model, clean, timesteps, noise, betas, labels
        )
        torch.testing.assert_close(observed[0], expected[0])
        torch.testing.assert_close(observed[1], expected[1])

    def test_ddim_final_state_matches_upstream_function(self):
        model = self._model()
        initial = torch.randn(4, 12)
        labels = torch.nn.functional.one_hot(
            torch.tensor([0, 1, 2, 1]), num_classes=3
        )
        betas = quadratic_beta_schedule(
            beta_start=0.0001, beta_end=0.02, timesteps=8
        )
        torch.manual_seed(15)
        observed = ddim_trajectory(
            initial.clone(),
            labels,
            model,
            betas,
            sequence=range(8),
            snapshot_timesteps=(8, 2, 0),
        )[0]
        torch.manual_seed(15)
        expected = upstream_denoising_module().generalized_steps(
            initial.clone(), range(8), model, betas, y=labels, eta=0.0
        )[0][-1]
        torch.testing.assert_close(observed, expected)

    def test_stochastic_ddim_is_deterministic_with_explicit_generator(self):
        model = self._model()
        initial = torch.randn(4, 12)
        labels = torch.nn.functional.one_hot(
            torch.tensor([0, 1, 2, 1]), num_classes=3
        )
        betas = quadratic_beta_schedule(
            beta_start=0.0001, beta_end=0.02, timesteps=8
        )
        outputs = []
        for _ in range(2):
            outputs.append(
                ddim_trajectory(
                    initial.clone(),
                    labels,
                    model,
                    betas,
                    sequence=range(8),
                    snapshot_timesteps=(0,),
                    eta=0.5,
                    generator=torch.Generator().manual_seed(41),
                )[0]
            )
        torch.testing.assert_close(outputs[0], outputs[1])


class PaperConfigurationTests(unittest.TestCase):
    def test_paper_configuration_rejects_proxy_architecture(self):
        source = Path("configs/rna_diffusion/archs4_mouse_paper_parity.yaml")
        text = source.read_text(encoding="utf-8").replace(
            "hidden_dims: [8192, 8192]", "hidden_dims: [512, 512]"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.yaml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "model.hidden_dims"):
                load_config(path)

    def test_osdr_extension_retains_exact_model_and_training_contract(self):
        config = load_conditional_config(
            "configs/rna_diffusion/osdr_pooled_flt_gc_paper_architecture.yaml"
        )
        self.assertEqual(config["model"]["hidden_dims"], [8192, 8192])
        self.assertEqual(config["training"]["epochs"], 15000)
        self.assertEqual(
            config["data"]["conditioning_covariates"], ["tissue", "condition"]
        )
        pretrained = load_conditional_config(
            "configs/rna_diffusion/osdr_archs4_pretrain_flt_gc_paper_architecture.yaml"
        )
        self.assertEqual(
            pretrained["training"]["regime"],
            "archs4_pretrain_osdr_finetune",
        )

    def test_targeted_reference_retains_exact_training_contract(self):
        config = load_config(
            "configs/rna_diffusion/archs4_liver_enriched_paper_native.yaml"
        )
        self.assertEqual(config["model"]["hidden_dims"], [8192, 8192])
        self.assertEqual(config["training"]["epochs"], 15000)
        self.assertEqual(config["data"]["profiles_per_tissue"]["liver"], 9466)


class ConditionalDataTests(unittest.TestCase):
    def test_explicit_accession_split_keeps_whole_studies(self):
        rows = pd.DataFrame(
            {
                "accession": ["train"] * 2 + ["val-a"] * 2 + ["val-b"] * 2 + ["test"] * 2,
                "condition": ["flight", "ground_control"] * 4,
            }
        )
        roles, audit = _explicit_accession_roles(
            rows,
            validation_accessions=["val-a", "val-b"],
            test_accessions=["test"],
        )
        self.assertEqual(set(roles.iloc[:2]), {"train"})
        self.assertEqual(set(roles.iloc[2:6]), {"validation"})
        self.assertEqual(set(roles.iloc[6:]), {"test"})
        self.assertEqual(
            audit["accessions_by_role"]["validation"], ["val-a", "val-b"]
        )
    def test_targeted_archs4_selection_obeys_tissue_quotas(self):
        metadata = pd.DataFrame(
            {
                "canonical_tissue": ["liver"] * 5 + ["kidney"] * 3,
                "selection_rank_within_tissue": [5, 1, 4, 2, 3, 3, 1, 2],
                "archs4_sample_index": np.arange(8),
            }
        )
        selected, available = _targeted_candidates(
            metadata, quotas={"liver": 4, "kidney": 2}, seed=7
        )
        self.assertEqual(selected["canonical_tissue"].value_counts()["liver"], 4)
        self.assertEqual(selected["canonical_tissue"].value_counts()["kidney"], 2)
        self.assertEqual(available, {"kidney": 3, "liver": 5})

    def test_series_split_has_no_group_leakage_and_keeps_classes(self):
        rows = []
        for tissue in ("liver", "kidney"):
            for series in range(20):
                rows.extend(
                    {
                        "canonical_tissue": tissue,
                        "series_id": f"{tissue}-{series}",
                    }
                    for _ in range(3)
                )
        metadata = pd.DataFrame(rows)
        partitions, groups = _group_split_indices(
            metadata,
            fractions={"train": 0.7, "validation": 0.15, "test": 0.15},
            seed=11,
            group_column="series_id",
        )
        role_groups = {role: set(groups[index]) for role, index in partitions.items()}
        self.assertTrue(role_groups["train"].isdisjoint(role_groups["validation"]))
        self.assertTrue(role_groups["train"].isdisjoint(role_groups["test"]))
        self.assertTrue(role_groups["validation"].isdisjoint(role_groups["test"]))
        for indices in partitions.values():
            self.assertEqual(
                set(metadata.loc[indices, "canonical_tissue"]), {"liver", "kidney"}
            )

    def test_within_study_split_retains_each_stratum_in_training(self):
        rows = pd.DataFrame(
            {
                "profile_id": [f"p{index}" for index in range(12)],
                "accession": ["a"] * 6 + ["b"] * 6,
                "tissue": ["liver"] * 12,
                "condition": ["flight"] * 3
                + ["ground_control"] * 3
                + ["flight"] * 3
                + ["ground_control"] * 3,
            }
        )
        roles = _within_study_roles(
            rows, seed=3, validation_fraction=0.2, test_fraction=0.2
        )
        observed = rows.assign(role=roles).groupby(
            ["accession", "condition", "role"]
        ).size()
        for accession in ("a", "b"):
            for condition in ("flight", "ground_control"):
                self.assertEqual(observed[accession, condition, "train"], 1)
                self.assertEqual(observed[accession, condition, "validation"], 1)
                self.assertEqual(observed[accession, condition, "test"], 1)

    def test_tpm_denominator_uses_genes_outside_landmark_panel(self):
        matrix = np.asarray(
            [[100.0, 100.0, 0.0], [100.0, 0.0, 900.0]], dtype=np.float32
        )
        adata = ad.AnnData(matrix)
        observed = _full_transcriptome_tpm(
            adata,
            np.asarray([0, 1]),
            np.asarray([0]),
            np.asarray([1000.0, 2000.0, 1000.0]),
            chunk_size=1,
        )
        np.testing.assert_allclose(
            observed[:, 0], [100.0 / 150.0 * 1e6, 100.0 / 1000.0 * 1e6]
        )

    def test_deseq2_median_of_ratios_removes_library_scale(self):
        base = np.asarray([5.0, 11.0, 17.0, 23.0])
        matrix = np.stack([base, base * 2.0, base * 4.0])
        normalized, size_factors, audit = _deseq2_median_of_ratios(matrix)
        np.testing.assert_allclose(normalized[0], normalized[1], rtol=1e-6)
        np.testing.assert_allclose(normalized[0], normalized[2], rtol=1e-6)
        np.testing.assert_allclose(
            size_factors / size_factors[0], [1.0, 2.0, 4.0], rtol=1e-6
        )
        self.assertEqual(audit["eligible_positive_genes"], 4)

    def test_deseq2_median_of_ratios_rounds_like_official_pipeline(self):
        matrix = np.asarray([[1.2, 4.8], [2.2, 9.7]], dtype=np.float64)
        normalized, _, audit = _deseq2_median_of_ratios(matrix)
        self.assertTrue(np.isfinite(normalized).all())
        self.assertEqual(audit["fractional_input_fraction"], 1.0)
        self.assertAlmostEqual(audit["maximum_rounding_distance"], 0.3)

    def test_joint_classes_preserve_named_covariates(self):
        rows = pd.DataFrame(
            {
                "tissue": ["liver", "liver"],
                "condition": ["flight", "ground_control"],
            }
        )
        labels = _joint_class_labels(rows, ("tissue", "condition"))
        self.assertEqual(
            labels.tolist(),
            [
                "tissue=liver||condition=flight",
                "tissue=liver||condition=ground_control",
            ],
        )

    def test_pretrained_tissue_columns_map_to_reference_classes(self):
        source = {
            "y_emb.weight": torch.arange(4, dtype=torch.float32).reshape(2, 2),
            "mid.0.w1.weight": torch.arange(7, dtype=torch.float32).reshape(1, 7),
            "mid.0.w1.bias": torch.asarray([4.0]),
        }
        template = {
            "y_emb.weight": torch.zeros(3, 2),
            "mid.0.w1.weight": torch.zeros(1, 9),
            "mid.0.w1.bias": torch.zeros(1),
        }
        expanded, audit = _expanded_condition_state(
            template,
            source,
            old_classes=["a", "b"],
            new_classes=[
                "tissue=a||condition=flight",
                "tissue=a||condition=reference",
                "tissue=b||condition=reference",
            ],
            embedding_dim=2,
        )
        torch.testing.assert_close(expanded["mid.0.w1.weight"][:, :3], source["mid.0.w1.weight"][:, :3])
        torch.testing.assert_close(expanded["mid.0.w1.weight"][:, 5:7], source["mid.0.w1.weight"][:, 3:5])
        torch.testing.assert_close(expanded["mid.0.w1.weight"][:, 7:9], source["mid.0.w1.weight"][:, 5:7])
        torch.testing.assert_close(expanded["mid.0.w1.weight"][:, 3:5], torch.zeros(1, 2))
        torch.testing.assert_close(expanded["y_emb.weight"][0], torch.zeros(2))
        torch.testing.assert_close(expanded["y_emb.weight"][1], source["y_emb.weight"][0])
        torch.testing.assert_close(expanded["y_emb.weight"][2], source["y_emb.weight"][1])
        self.assertEqual(audit["mapped_classes"], 2)

    def test_function_preserving_tissue_expansion_matches_source_outputs(self):
        old_classes = ["a", "b"]
        new_classes = [
            "tissue=a||condition=flight",
            "tissue=a||condition=ground_control",
            "tissue=a||condition=reference",
            "tissue=b||condition=flight",
            "tissue=b||condition=ground_control",
            "tissue=b||condition=reference",
        ]
        torch.manual_seed(17)
        old_model = upstream_model_class()(
            model_config(expression_dim=12, num_classes=2, model=SMALL_MODEL)
        ).eval()
        torch.manual_seed(23)
        new_model = upstream_model_class()(
            model_config(expression_dim=12, num_classes=6, model=SMALL_MODEL)
        ).eval()
        expanded, audit = _expanded_condition_state(
            new_model.state_dict(),
            old_model.state_dict(),
            old_classes=old_classes,
            new_classes=new_classes,
            embedding_dim=2,
            strategy="function_preserving_tissue",
        )
        new_model.load_state_dict(expanded)
        expression = torch.randn(4, 12)
        timesteps = torch.tensor([0, 2, 5, 7])
        for old_index, new_indices in ((0, (0, 1, 2)), (1, (3, 4, 5))):
            old_labels = torch.nn.functional.one_hot(
                torch.full((4,), old_index), num_classes=2
            )
            expected = old_model(expression, timesteps, old_labels)
            for new_index in new_indices:
                new_labels = torch.nn.functional.one_hot(
                    torch.full((4,), new_index), num_classes=6
                )
                observed = new_model(expression, timesteps, new_labels)
                torch.testing.assert_close(observed, expected, atol=2e-6, rtol=2e-6)
        self.assertEqual(audit["mapped_new_joint_classes"], 6)
        self.assertEqual(audit["copied_one_hot_embedding_rows"], 2)

    def test_amp_skip_detection_uses_scale_backoff(self):
        class FakeScaler:
            def __init__(self, before, after):
                self.values = iter((before, after))

            def get_scale(self):
                return next(self.values)

            def step(self, optimizer):
                return None

            def update(self):
                return None

        succeeded, before, after = _scaled_optimizer_step(
            FakeScaler(65536.0, 32768.0), object()
        )
        self.assertFalse(succeeded)
        self.assertEqual((before, after), (65536.0, 32768.0))

        succeeded, _, _ = _scaled_optimizer_step(
            FakeScaler(32768.0, 32768.0), object()
        )
        self.assertTrue(succeeded)


class RealEffectCeilingTests(unittest.TestCase):
    def test_coherent_cross_accession_effect_has_high_ceiling(self):
        rng = np.random.default_rng(7)
        rows = []
        expression = []
        effect = np.linspace(-1.5, 1.5, 24)
        for accession_index in range(4):
            offset = rng.normal(0.0, 0.4, len(effect))
            for condition in ("ground_control", "flight"):
                for _ in range(6):
                    rows.append(
                        {
                            "tissue": "skeletal_muscle",
                            "accession": f"OSD-{accession_index}",
                            "condition": condition,
                        }
                    )
                    expression.append(
                        offset
                        + (effect if condition == "flight" else 0.0)
                        + rng.normal(0.0, 0.08, len(effect))
                    )
        summary, detail = analyze_tissue(
            np.asarray(expression),
            pd.DataFrame(rows),
            tissue="skeletal_muscle",
            permutation_repeats=20,
            seed=3,
        )
        self.assertEqual(summary["eligible_accessions"], 4)
        self.assertGreater(summary["loo_effect_correlation_median"], 0.9)
        self.assertGreater(summary["loo_balanced_accuracy_median"], 0.9)
        self.assertEqual(len(detail), 4)

    def test_loader_rejects_locked_test(self):
        with self.assertRaisesRegex(ValueError, "locked test"):
            load_development_expression("unused.h5", "unused.tsv", roles=["test"])


class FactorizedAdapterTests(unittest.TestCase):
    def _schema(self):
        samples = pd.DataFrame(
            {
                "tissue": ["liver", "liver", "skeletal_muscle"],
                "condition": ["flight", "ground_control", "flight"],
                "sex": ["female", "male", "female"],
                "muscle_group": ["unknown", "unknown", "soleus"],
            }
        )
        return samples, build_factorized_schema(
            samples, ["liver", "skeletal_muscle"]
        )

    def test_zero_adapter_preserves_pretrained_function(self):
        samples, schema = self._schema()
        torch.manual_seed(31)
        base = upstream_model_class()(
            model_config(expression_dim=12, num_classes=2, model=SMALL_MODEL)
        ).eval()
        adapter = FactorizedAdapterDDIM(
            base, schema, domain_lora_rank=4, domain_lora_alpha=4
        ).eval()
        labels = torch.from_numpy(encode_factorized_labels(samples, schema))
        expression = torch.randn(3, 12)
        timesteps = torch.tensor([1.0, 3.0, 6.0])
        expected = base(expression, timesteps, labels[:, :2])
        observed = adapter(expression, timesteps, labels)
        torch.testing.assert_close(observed, expected)

    def test_stage_selection_exposes_only_requested_adapter_group(self):
        _, schema = self._schema()
        base = upstream_model_class()(
            model_config(expression_dim=12, num_classes=2, model=SMALL_MODEL)
        )
        adapter = FactorizedAdapterDDIM(base, schema)
        adapter.set_trainable_groups(["condition"])
        trainable = [
            name for name, value in adapter.named_parameters() if value.requires_grad
        ]
        self.assertTrue(trainable)
        self.assertTrue(all("layers.condition" in name for name in trainable))

        adapter = FactorizedAdapterDDIM(base, schema, domain_lora_rank=2)
        adapter.set_trainable_groups(["domain"])
        trainable = [
            name for name, value in adapter.named_parameters() if value.requires_grad
        ]
        self.assertTrue(any("_domain_lora." in name for name in trainable))

    def test_condition_neutralization_keeps_base_and_domain_labels(self):
        samples, schema = self._schema()
        labels = torch.from_numpy(encode_factorized_labels(samples, schema))
        neutral = neutralize_group(labels, schema, "condition")
        self.assertTrue(torch.equal(labels[:, : schema.base_width], neutral[:, : schema.base_width]))
        condition = schema.group_slices()["condition"]
        start = schema.base_width + condition.start
        stop = schema.base_width + condition.stop
        self.assertEqual(int(neutral[:, start:stop].sum()), 0)
        self.assertGreater(int(labels[:, start:stop].sum()), 0)

    def test_study_conditioning_adds_main_and_interaction_features(self):
        samples, _ = self._schema()
        samples["accession"] = ["OSD-1", "OSD-1", "OSD-2"]
        schema = build_factorized_schema(
            samples,
            ["liver", "skeletal_muscle"],
            include_study=True,
        )
        self.assertIn("study=OSD-1", schema.groups["domain"])
        self.assertIn(
            "study_condition=OSD-2::flight", schema.groups["condition"]
        )
        labels = encode_factorized_labels(samples, schema)
        self.assertEqual(labels.shape[1], schema.total_width)

    def test_covariance_calibration_matches_training_correlation_structure(self):
        rng = np.random.default_rng(91)
        real = rng.normal(size=(800, 12)) @ rng.normal(size=(12, 12))
        synthetic = rng.normal(size=(800, 12)) @ rng.normal(size=(12, 12))
        calibrated = CovarianceCalibrator(1e-6).fit(real, synthetic).apply(
            synthetic
        )
        upper = np.triu_indices(real.shape[1], 1)
        agreement = np.corrcoef(
            np.corrcoef(real, rowvar=False)[upper],
            np.corrcoef(calibrated, rowvar=False)[upper],
        )[0, 1]
        self.assertGreater(agreement, 0.999)


class EvaluationMetricTests(unittest.TestCase):
    def test_per_tissue_fidelity_reports_each_sufficient_group(self):
        rng = np.random.default_rng(9)
        real = rng.normal(size=(12, 8))
        synthetic = real + rng.normal(scale=0.1, size=real.shape)
        samples = pd.DataFrame({"tissue": ["a"] * 6 + ["b"] * 6})
        table = _per_tissue_fidelity(real, synthetic, samples)
        self.assertEqual(table["tissue"].tolist(), ["a", "b"])
        self.assertEqual(table["profiles"].tolist(), [6, 6])
        self.assertTrue(
            table["adversarial_accuracy"].between(0.0, 1.0).all()
        )

    def test_class_probe_excludes_training_classes_absent_from_evaluation(self):
        real_train = np.asarray(
            [[-2.0], [-1.0], [1.0], [2.0], [100.0], [101.0]], dtype=float
        )
        synthetic = np.asarray([[-1.5], [1.5]], dtype=float)
        result = _class_probe(
            real_train,
            ["a", "a", "b", "b", "unused", "unused"],
            synthetic,
            ["a", "b"],
            seed=0,
        )
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["train_profiles"], 4)
        self.assertEqual(result["excluded_train_profiles"], 2)
        self.assertEqual(result["evaluation_classes"], ["a", "b"])
        self.assertEqual(result["balanced_accuracy"], 1.0)

    def test_nearest_neighbor_adversarial_accuracy_extremes(self):
        real = np.asarray([[0.0], [1.0], [2.0], [3.0]])
        self.assertEqual(
            _nearest_neighbor_adversarial_accuracy(real, real.copy()), 0.0
        )
        far = real + 100.0
        self.assertEqual(_nearest_neighbor_adversarial_accuracy(real, far), 1.0)


if __name__ == "__main__":
    unittest.main()
