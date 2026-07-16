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


class ConditionalDataTests(unittest.TestCase):
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
