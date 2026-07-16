from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import numpy as np
import torch

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


class EvaluationMetricTests(unittest.TestCase):
    def test_nearest_neighbor_adversarial_accuracy_extremes(self):
        real = np.asarray([[0.0], [1.0], [2.0], [3.0]])
        self.assertEqual(
            _nearest_neighbor_adversarial_accuracy(real, real.copy()), 0.0
        )
        far = real + 100.0
        self.assertEqual(_nearest_neighbor_adversarial_accuracy(real, far), 1.0)


if __name__ == "__main__":
    unittest.main()
