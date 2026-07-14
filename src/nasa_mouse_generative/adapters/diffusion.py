"""Lacan et al. landmark-space DDPM/DDIM adapter."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Any

import numpy as np
import torch
from sklearn.linear_model import Ridge

from nasa_mouse_diffusion.data import select_landmarks
from nasa_mouse_diffusion.diffusion import beta_schedule, noise_estimation_loss, sample
from nasa_mouse_diffusion.model import ConditionalDiffusionMLP

from ..training_data import DataPartition
from .base import ModelAdapter, weighted_loader


class DiffusionAdapter(ModelAdapter):
    adapter_id = "lacan_diffusion"
    supports_generation = True

    def __init__(
        self,
        *,
        reconstruction_matrix: np.ndarray | None,
        l1000_map: str,
        serialized_payload: dict[str, Any] | None = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.l1000_map = str(l1000_map)
        checkpoint_payload = None
        if self.resume and self.checkpoint_path.exists():
            checkpoint_payload = torch.load(
                self.checkpoint_path, map_location=self.device, weights_only=False
            )
        initialization_payload = serialized_payload or checkpoint_payload
        if initialization_payload is not None:
            self.landmark_genes = list(initialization_payload["landmark_genes"])
            self.target_genes = list(initialization_payload["target_genes"])
            self.landmark_indices = np.asarray(
                initialization_payload["landmark_indices"], dtype=int
            )
            self.target_indices = np.asarray(
                initialization_payload["target_indices"], dtype=int
            )
            self.landmark_source = str(initialization_payload["landmark_source"])
            self.model_config = dict(initialization_payload["model_config"])
            self.reconstruction_coef = np.asarray(
                initialization_payload["reconstruction_coef"], dtype=np.float32
            )
            self.reconstruction_intercept = np.asarray(
                initialization_payload["reconstruction_intercept"], dtype=np.float32
            )
        else:
            if reconstruction_matrix is None:
                raise ValueError("New diffusion adapters require reconstruction data")
            strategy = str(self.parameters.get("landmark_strategy", "l1000_or_hvg"))
            n_landmarks = min(
                int(self.parameters.get("n_landmarks", 512)), len(self.genes)
            )
            (
                self.landmark_genes,
                self.target_genes,
                self.landmark_indices,
                self.target_indices,
                self.landmark_source,
            ) = select_landmarks(
                reconstruction_matrix,
                self.genes,
                n_landmarks=n_landmarks,
                strategy=strategy,
                l1000_map=l1000_map,
                min_l1000=min(300, max(1, n_landmarks)),
            )
            self.model_config = {
                "expression_dim": len(self.landmark_genes),
                "categorical_cardinalities": self.cardinalities,
                "hidden_dim": int(self.parameters.get("hidden_dim", 512)),
                "n_blocks": int(self.parameters.get("n_blocks", 2)),
                "dropout": float(self.parameters.get("dropout", 0.1)),
                "time_embedding_dim": int(
                    self.parameters.get("time_embedding_dim", 64)
                ),
                "categorical_embedding_dim": int(
                    self.parameters.get("categorical_embedding_dim", 8)
                ),
                "sinusoidal_time": bool(
                    self.parameters.get("sinusoidal_time", False)
                ),
                "num_timesteps": int(
                    self.parameters.get("diffusion_timesteps", 1000)
                ),
            }
            self.reconstruction_coef, self.reconstruction_intercept = (
                self._fit_reconstructor(reconstruction_matrix)
            )
        self.model = ConditionalDiffusionMLP(**self.model_config).to(self.device)
        self.betas = torch.as_tensor(
            beta_schedule(
                str(self.parameters.get("beta_schedule", "quad")),
                beta_start=float(self.parameters.get("beta_start", 0.0001)),
                beta_end=float(self.parameters.get("beta_end", 0.02)),
                timesteps=self.model_config["num_timesteps"],
            ),
            dtype=torch.float32,
            device=self.device,
        )
        self.ema_state = {
            key: value.detach().clone()
            for key, value in self.model.state_dict().items()
        }
        self._resume_payload: dict[str, Any] | None = None
        if checkpoint_payload is not None:
            payload = checkpoint_payload
            self._validate_payload(payload)
            self.model.load_state_dict(payload["model_state_dict"])
            self.ema_state = {
                key: value.to(self.device)
                for key, value in payload.get("ema_state_dict", {}).items()
            } or {
                key: value.detach().clone()
                for key, value in self.model.state_dict().items()
            }
            self.reconstruction_coef = np.asarray(
                payload["reconstruction_coef"], dtype=np.float32
            )
            self.reconstruction_intercept = np.asarray(
                payload["reconstruction_intercept"], dtype=np.float32
            )
            self._restore_common(payload)
            self._resume_payload = payload

    def _validate_payload(self, payload: dict[str, Any]) -> None:
        if payload.get("adapter_id") != self.adapter_id:
            raise ValueError("Checkpoint belongs to another adapter")
        if payload.get("genes") != self.genes:
            raise ValueError("Checkpoint gene order differs from this run")
        if list(payload.get("landmark_indices", [])) != self.landmark_indices.tolist():
            raise ValueError("Checkpoint landmark genes differ from this run")

    def _fit_reconstructor(self, matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if len(self.target_indices) == 0:
            return (
                np.empty((0, len(self.landmark_indices)), dtype=np.float32),
                np.empty(0, dtype=np.float32),
            )
        limit = int(self.parameters.get("reconstruction_samples", 10000))
        values = np.asarray(matrix, dtype=np.float32)
        if limit > 0 and len(values) > limit:
            rng = np.random.default_rng(self.seed)
            values = values[rng.choice(len(values), limit, replace=False)]
        model = Ridge(
            alpha=float(self.parameters.get("reconstruction_alpha", 1.0)),
            fit_intercept=True,
        )
        model.fit(values[:, self.landmark_indices], values[:, self.target_indices])
        return (
            np.asarray(model.coef_, dtype=np.float32),
            np.asarray(model.intercept_, dtype=np.float32),
        )

    def _reconstruct(self, landmarks: np.ndarray) -> np.ndarray:
        full = np.zeros((len(landmarks), len(self.genes)), dtype=np.float32)
        full[:, self.landmark_indices] = landmarks
        if len(self.target_indices):
            full[:, self.target_indices] = (
                landmarks @ self.reconstruction_coef.T
                + self.reconstruction_intercept.reshape(1, -1)
            )
        return full

    def _update_ema(self) -> None:
        if not bool(self.parameters.get("use_ema", True)):
            return
        decay = float(self.parameters.get("ema_decay", 0.999))
        current = self.model.state_dict()
        for key, value in current.items():
            if value.is_floating_point():
                self.ema_state[key].mul_(decay).add_(value.detach(), alpha=1.0 - decay)
            else:
                self.ema_state[key].copy_(value)

    @contextmanager
    def _ema_weights(self):
        if not bool(self.parameters.get("use_ema", True)):
            yield
            return
        current = {
            key: value.detach().clone()
            for key, value in self.model.state_dict().items()
        }
        self.model.load_state_dict(self.ema_state)
        try:
            yield
        finally:
            self.model.load_state_dict(current)

    def _save_checkpoint(self, stage: str, optimizer, scaler) -> None:
        payload = self._common_payload()
        payload.update(
            {
                "model_config": self.model_config,
                "model_state_dict": self.model.state_dict(),
                "ema_state_dict": self.ema_state,
                "betas": self.betas.detach().cpu(),
                "active_stage": stage,
                "optimizer_state_dict": optimizer.state_dict(),
                "scaler_state_dict": scaler.state_dict(),
                "landmark_genes": self.landmark_genes,
                "target_genes": self.target_genes,
                "landmark_indices": self.landmark_indices.tolist(),
                "target_indices": self.target_indices.tolist(),
                "landmark_source": self.landmark_source,
                "l1000_map": self.l1000_map,
                "reconstruction_coef": self.reconstruction_coef,
                "reconstruction_intercept": self.reconstruction_intercept,
            }
        )
        self._atomic_torch_save(payload, self.checkpoint_path)

    def fit_stage(
        self, partition: DataPartition, *, stage: str, epochs: int, learning_rate: float
    ) -> list[dict[str, Any]]:
        if len(partition) < 2:
            raise ValueError(f"Diffusion stage {stage} needs at least two profiles")
        completed = int(self.state.completed_epochs.get(stage, 0))
        if completed >= int(epochs):
            return [row for row in self.state.history if row.get("stage") == stage]
        optimizer = torch.optim.AdamW(
            self.model.parameters(), lr=float(learning_rate), eps=1e-8
        )
        amp = bool(self.parameters.get("use_amp", True)) and self.device.type == "cuda"
        scaler = torch.amp.GradScaler("cuda", enabled=amp)
        if self._resume_payload and self._resume_payload.get("active_stage") == stage:
            optimizer.load_state_dict(self._resume_payload["optimizer_state_dict"])
            if self._resume_payload.get("scaler_state_dict"):
                scaler.load_state_dict(self._resume_payload["scaler_state_dict"])
        for epoch in range(completed + 1, int(epochs) + 1):
            loader = weighted_loader(
                partition,
                batch_size=self.batch_size,
                seed=self.seed + epoch + 10000 * len(self.state.completed_epochs),
                num_workers=self.num_workers,
            )
            losses: list[float] = []
            errors: list[float] = []
            self.model.train()
            for full, categories in loader:
                expression = full[:, self.landmark_indices].to(self.device)
                categories = categories.to(self.device)
                optimizer.zero_grad(set_to_none=True)
                with torch.autocast(
                    device_type=self.device.type, enabled=amp, dtype=torch.float16
                ):
                    loss, error = noise_estimation_loss(
                        self.model, expression, categories, self.betas
                    )
                if not torch.isfinite(loss):
                    raise FloatingPointError(
                        f"Non-finite diffusion loss at {stage} epoch {epoch}"
                    )
                scaler.scale(loss).backward()
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
                scaler.step(optimizer)
                scaler.update()
                self._update_ema()
                losses.append(float(loss.detach().cpu()))
                errors.append(float(error.detach().cpu()))
                self.state.global_steps += 1
            row = {
                "stage": stage,
                "epoch": epoch,
                "learning_rate": float(learning_rate),
                "loss": float(np.mean(losses)),
                "noise_abs_error": float(np.mean(errors)),
            }
            self.state.history.append(row)
            self.state.completed_epochs[stage] = epoch
            if epoch % self.checkpoint_every == 0 or epoch == int(epochs):
                self._save_checkpoint(stage, optimizer, scaler)
                self.write_history()
        self._resume_payload = None
        return [row for row in self.state.history if row.get("stage") == stage]

    def encode(self, partition: DataPartition) -> np.ndarray:
        features = []
        self.model.eval()
        with self._ema_weights(), torch.no_grad():
            for start in range(0, len(partition), self.batch_size):
                end = min(start + self.batch_size, len(partition))
                expression = torch.as_tensor(
                    partition.matrix[start:end, self.landmark_indices],
                    dtype=torch.float32,
                    device=self.device,
                )
                categories = torch.as_tensor(
                    partition.categories[start:end],
                    dtype=torch.long,
                    device=self.device,
                )
                timesteps = torch.zeros(
                    end - start, dtype=torch.long, device=self.device
                )
                features.append(
                    self.model.features(expression, timesteps, categories)
                    .detach()
                    .cpu()
                    .numpy()
                )
        return np.concatenate(features).astype(np.float32)

    def generate(
        self, categories: np.ndarray, *, seed: int, batch_size: int | None = None
    ) -> np.ndarray:
        batch_size = int(batch_size or self.batch_size)
        categories = np.asarray(categories, dtype=np.int64)
        sample_steps = int(self.parameters.get("sample_steps", 50))
        eta = float(self.parameters.get("eta", 0.0))
        generator = torch.Generator(device=self.device)
        generator.manual_seed(int(seed))
        generated = []
        self.model.eval()
        with self._ema_weights(), torch.no_grad():
            for start in range(0, len(categories), batch_size):
                end = min(start + batch_size, len(categories))
                noise = torch.randn(
                    (end - start, len(self.landmark_indices)),
                    generator=generator,
                    device=self.device,
                )
                landmarks = sample(
                    self.model,
                    categories[start:end],
                    betas=self.betas,
                    sample_steps=sample_steps,
                    eta=eta,
                    noise=noise,
                    device=self.device,
                )
                generated.append(landmarks.cpu().numpy())
        if not generated:
            return np.empty((0, len(self.genes)), dtype=np.float32)
        result = self._reconstruct(np.concatenate(generated).astype(np.float32))
        if not np.isfinite(result).all():
            raise FloatingPointError("Diffusion generated non-finite expression")
        return result

    def save_final(self) -> Path:
        path = self.output_dir / "model.pt"
        payload = self._common_payload()
        payload.update(
            {
                "model_config": self.model_config,
                "model_state_dict": self.model.state_dict(),
                "ema_state_dict": self.ema_state,
                "betas": self.betas.detach().cpu(),
                "landmark_genes": self.landmark_genes,
                "target_genes": self.target_genes,
                "landmark_indices": self.landmark_indices.tolist(),
                "target_indices": self.target_indices.tolist(),
                "landmark_source": self.landmark_source,
                "l1000_map": self.l1000_map,
                "reconstruction_coef": self.reconstruction_coef,
                "reconstruction_intercept": self.reconstruction_intercept,
            }
        )
        self._atomic_torch_save(payload, path)
        self.write_history()
        self.write_adapter_summary()
        return path

    @classmethod
    def load(
        cls, output_dir: Path, *, device_spec: str = "auto"
    ) -> "DiffusionAdapter":
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
            reconstruction_matrix=None,
            l1000_map=str(payload.get("l1000_map", "")),
            serialized_payload=payload,
        )
        adapter.model.load_state_dict(payload["model_state_dict"])
        adapter.ema_state = {
            key: value.to(adapter.device)
            for key, value in payload["ema_state_dict"].items()
        }
        adapter.betas = payload["betas"].to(adapter.device)
        adapter._restore_common(payload)
        adapter.model.eval()
        return adapter
