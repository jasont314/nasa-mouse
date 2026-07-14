"""Vinas et al. conditional WGAN-GP adapter."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import torch

from nasa_mouse_wgan.model import ConditionalWGANGP
from nasa_mouse_wgan.training import TrainConfig, critic_features, train_epoch

from ..training_data import DataPartition
from .base import ModelAdapter, weighted_loader


class WGANAdapter(ModelAdapter):
    adapter_id = "vinas_wgan_gp"
    supports_generation = True

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        hidden_dims = tuple(
            int(value) for value in self.parameters.get("hidden_dims", [256, 256])
        )
        self.model_config = {
            "expression_dim": len(self.genes),
            "categorical_cardinalities": self.cardinalities,
            "noise_dim": int(self.parameters.get("noise_dim", 64)),
            "hidden_dims": hidden_dims,
        }
        self.model = ConditionalWGANGP(**self.model_config).to(self.device)
        self._resume_payload: dict[str, Any] | None = None
        if self.resume and self.checkpoint_path.exists():
            payload = torch.load(
                self.checkpoint_path, map_location=self.device, weights_only=False
            )
            self._validate_payload(payload)
            self.model.load_state_dict(payload["model_state_dict"])
            self._restore_common(payload)
            self._resume_payload = payload

    def _validate_payload(self, payload: dict[str, Any]) -> None:
        if payload.get("adapter_id") != self.adapter_id:
            raise ValueError("Checkpoint belongs to another adapter")
        if payload.get("genes") != self.genes:
            raise ValueError("Checkpoint gene order differs from this run")
        if list(payload.get("cardinalities", [])) != self.cardinalities:
            raise ValueError("Checkpoint categorical cardinalities differ from this run")

    def _optimizers(self, learning_rate: float):
        name = str(self.parameters.get("optimizer", "adam")).lower()
        if name == "rmsprop":
            optim_g = torch.optim.RMSprop(
                self.model.generator.parameters(), lr=float(learning_rate)
            )
            optim_d = torch.optim.RMSprop(
                self.model.critic.parameters(), lr=float(learning_rate)
            )
        elif name == "adam":
            optim_g = torch.optim.Adam(
                self.model.generator.parameters(),
                lr=float(learning_rate),
                betas=(0.5, 0.9),
            )
            optim_d = torch.optim.Adam(
                self.model.critic.parameters(),
                lr=float(learning_rate),
                betas=(0.5, 0.9),
            )
        else:
            raise ValueError(f"Unsupported WGAN optimizer: {name}")
        return optim_g, optim_d

    def _save_checkpoint(self, stage: str, optim_g, optim_d) -> None:
        payload = self._common_payload()
        payload.update(
            {
                "model_config": self.model_config,
                "model_state_dict": self.model.state_dict(),
                "active_stage": stage,
                "optimizer_g_state_dict": optim_g.state_dict(),
                "optimizer_d_state_dict": optim_d.state_dict(),
            }
        )
        self._atomic_torch_save(payload, self.checkpoint_path)

    def fit_stage(
        self, partition: DataPartition, *, stage: str, epochs: int, learning_rate: float
    ) -> list[dict[str, Any]]:
        if len(partition) < 2:
            raise ValueError(f"WGAN stage {stage} needs at least two profiles")
        completed = int(self.state.completed_epochs.get(stage, 0))
        if completed >= int(epochs):
            return [row for row in self.state.history if row.get("stage") == stage]
        optim_g, optim_d = self._optimizers(learning_rate)
        if self._resume_payload and self._resume_payload.get("active_stage") == stage:
            optim_g.load_state_dict(self._resume_payload["optimizer_g_state_dict"])
            optim_d.load_state_dict(self._resume_payload["optimizer_d_state_dict"])
        train_config = TrainConfig(
            epochs=1,
            batch_size=self.batch_size,
            learning_rate=float(learning_rate),
            critic_steps=int(self.parameters.get("critic_steps", 5)),
            gradient_penalty=float(self.parameters.get("gradient_penalty", 10.0)),
            seed=self.seed,
        )
        for epoch in range(completed + 1, int(epochs) + 1):
            loader = weighted_loader(
                partition,
                batch_size=self.batch_size,
                seed=self.seed + epoch + 10000 * len(self.state.completed_epochs),
                num_workers=self.num_workers,
            )
            metrics = train_epoch(
                self.model,
                loader,
                config=train_config,
                optim_g=optim_g,
                optim_d=optim_d,
                device=self.device,
            )
            if not all(np.isfinite(value) for value in metrics.values()):
                raise FloatingPointError(f"Non-finite WGAN metrics at {stage} epoch {epoch}")
            self.state.global_steps += len(loader)
            row = {
                "stage": stage,
                "epoch": epoch,
                "learning_rate": float(learning_rate),
                **metrics,
            }
            self.state.history.append(row)
            self.state.completed_epochs[stage] = epoch
            if epoch % self.checkpoint_every == 0 or epoch == int(epochs):
                self._save_checkpoint(stage, optim_g, optim_d)
                self.write_history()
        self._resume_payload = None
        return [row for row in self.state.history if row.get("stage") == stage]

    def encode(self, partition: DataPartition) -> np.ndarray:
        scores, features = critic_features(
            self.model,
            partition.matrix,
            partition.categories,
            batch_size=self.batch_size,
            device=self.device,
        )
        return np.column_stack([scores, features]).astype(np.float32)

    def generate(
        self, categories: np.ndarray, *, seed: int, batch_size: int | None = None
    ) -> np.ndarray:
        batch_size = int(batch_size or self.batch_size)
        categories = np.asarray(categories, dtype=np.int64)
        generator = torch.Generator(device=self.device)
        generator.manual_seed(int(seed))
        generated = []
        self.model.eval()
        with torch.no_grad():
            for start in range(0, len(categories), batch_size):
                cats = torch.as_tensor(
                    categories[start : start + batch_size],
                    dtype=torch.long,
                    device=self.device,
                )
                noise = torch.randn(
                    (len(cats), self.model.noise_dim),
                    generator=generator,
                    device=self.device,
                )
                generated.append(self.model.generator(noise, cats).cpu().numpy())
        if not generated:
            return np.empty((0, len(self.genes)), dtype=np.float32)
        result = np.concatenate(generated).astype(np.float32)
        if not np.isfinite(result).all():
            raise FloatingPointError("WGAN generated non-finite expression")
        return result

    def save_final(self) -> Path:
        path = self.output_dir / "model.pt"
        payload = self._common_payload()
        payload.update(
            {
                "model_config": self.model_config,
                "model_state_dict": self.model.state_dict(),
            }
        )
        self._atomic_torch_save(payload, path)
        self.write_history()
        self.write_adapter_summary()
        return path

    @classmethod
    def load(cls, output_dir: Path, *, device_spec: str = "auto") -> "WGANAdapter":
        payload = torch.load(
            Path(output_dir) / "model.pt", map_location="cpu", weights_only=False
        )
        adapter = cls(
            genes=list(payload["genes"]),
            cardinalities=list(payload["cardinalities"]),
            covariates=tuple(payload["covariates"]),
            parameters=dict(payload["parameters"]),
            device_spec=device_spec,
            output_dir=Path(output_dir),
            checkpoint_every=1,
            resume=False,
            seed=int(payload["seed"]),
            num_workers=0,
        )
        adapter.model.load_state_dict(payload["model_state_dict"])
        adapter._restore_common(payload)
        adapter.model.eval()
        return adapter
