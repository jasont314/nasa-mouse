"""Official-architecture MOBER training and deterministic projection adapter."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F
from torch.utils.data import DataLoader, TensorDataset

from .harmonizers import Harmonizer, _batch_values, _metadata_frame, _resolve_batch_key


class _Encoder(nn.Module):
    def __init__(self, genes: int, encoding_dim: int) -> None:
        super().__init__()
        self.activation = nn.SELU()
        self.fc1 = nn.Linear(genes, 256)
        self.bn1 = nn.BatchNorm1d(256, momentum=0.01, eps=0.001)
        self.dropout1 = nn.Dropout(0.1)
        self.fc2 = nn.Linear(256, 128)
        self.bn2 = nn.BatchNorm1d(128, momentum=0.01, eps=0.001)
        self.dropout2 = nn.Dropout(0.1)
        self.means = nn.Linear(128, encoding_dim)
        self.log_variances = nn.Linear(128, encoding_dim)

    def forward(self, values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        hidden = self.dropout1(self.activation(self.bn1(self.fc1(values))))
        hidden = self.dropout2(self.activation(self.bn2(self.fc2(hidden))))
        means = self.means(hidden)
        standard_deviations = torch.exp(0.5 * self.log_variances(hidden)) + 1e-4
        latent = means + standard_deviations * torch.randn_like(standard_deviations)
        return means, standard_deviations, latent


class _Decoder(nn.Module):
    def __init__(self, genes: int, encoding_dim: int, batches: int) -> None:
        super().__init__()
        self.activation = nn.SELU()
        self.batch_layer = nn.Linear(batches, batches)
        self.batch_norm = nn.BatchNorm1d(batches, momentum=0.01, eps=0.001)
        self.fc1 = nn.Linear(encoding_dim + batches, 128)
        self.bn1 = nn.BatchNorm1d(128, momentum=0.01, eps=0.001)
        self.fc2 = nn.Linear(128, 256)
        self.bn2 = nn.BatchNorm1d(256, momentum=0.01, eps=0.001)
        self.output = nn.Linear(256, genes)

    def forward(self, latent: torch.Tensor, batches: torch.Tensor) -> torch.Tensor:
        batch_embedding = self.activation(self.batch_norm(self.batch_layer(batches)))
        hidden = torch.cat([latent, batch_embedding], dim=1)
        hidden = self.activation(self.bn1(self.fc1(hidden)))
        hidden = self.activation(self.bn2(self.fc2(hidden)))
        return torch.relu(self.output(hidden))


class _BatchVAE(nn.Module):
    def __init__(self, genes: int, encoding_dim: int, batches: int) -> None:
        super().__init__()
        self.encoder = _Encoder(genes, encoding_dim)
        self.decoder = _Decoder(genes, encoding_dim, batches)

    def forward(
        self, values: torch.Tensor, batches: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        means, standard_deviations, latent = self.encoder(values)
        decoded = self.decoder(latent, batches)
        return decoded, latent, means, standard_deviations


class _Adversary(nn.Module):
    def __init__(self, encoding_dim: int, batches: int) -> None:
        super().__init__()
        self.activation = nn.SELU()
        self.fc1 = nn.Linear(encoding_dim, encoding_dim)
        self.bn1 = nn.BatchNorm1d(encoding_dim, momentum=0.01, eps=0.001)
        self.fc2 = nn.Linear(encoding_dim, encoding_dim)
        self.bn2 = nn.BatchNorm1d(encoding_dim, momentum=0.01, eps=0.001)
        self.output = nn.Linear(encoding_dim, batches)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        hidden = self.activation(self.bn1(self.fc1(values)))
        hidden = self.activation(self.bn2(self.fc2(hidden)))
        return F.log_softmax(self.output(hidden), dim=1)


class MoberHarmonizer(Harmonizer):
    method_id = "mober"
    is_transductive = False

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.batch_key = ""
        self.batch_levels: list[str] = []
        self.target_batch = ""
        self.encoding_dim = int(self.parameters.get("encoding_dim", 64))
        self.genes = 0
        self.model: _BatchVAE | None = None
        self.history: list[dict[str, float]] = []
        self.device = self._resolve_device(self.device_spec)

    @staticmethod
    def _resolve_device(device_spec: str) -> torch.device:
        if device_spec == "auto":
            return torch.device("cuda" if torch.cuda.is_available() else "cpu")
        if device_spec == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("MOBER requested CUDA, but CUDA is unavailable")
        return torch.device(device_spec)

    def _batch_encoding(self, batches: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        mapping = {level: index for index, level in enumerate(self.batch_levels)}
        unknown = sorted(set(batches).difference(mapping))
        if unknown:
            raise ValueError(
                f"MOBER cannot encode unseen batch levels for {self.batch_key}: {unknown}"
            )
        indices = np.asarray([mapping[level] for level in batches], dtype=np.int64)
        one_hot = np.eye(len(self.batch_levels), dtype=np.float32)[indices]
        return indices, one_hot

    def _choose_target(self, frame: pd.DataFrame, batches: np.ndarray) -> str:
        requested = str(self.parameters.get("target_batch", "auto"))
        if requested != "auto":
            if requested not in self.batch_levels:
                raise ValueError(
                    f"MOBER target_batch={requested!r} is not a training batch"
                )
            return requested
        if self.batch_key == "source" and "osdr" in self.batch_levels:
            return "osdr"
        counts = pd.Series(batches).value_counts()
        maximum = int(counts.max())
        return sorted(counts.loc[counts.eq(maximum)].index.astype(str))[0]

    def fit_transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
    ) -> np.ndarray:
        values = np.asarray(matrix, dtype=np.float32)
        if values.ndim != 2 or len(values) < 4:
            raise ValueError("MOBER requires at least four samples")
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(
                "MOBER requires finite non-negative expression before final scaling"
            )
        frame = _metadata_frame(studies, metadata, len(values))
        self.batch_key = _resolve_batch_key(frame, self.parameters)
        batches = _batch_values(frame, self.batch_key)
        self.batch_levels = sorted(set(batches))
        if len(self.batch_levels) < 2:
            raise ValueError("MOBER requires at least two training batches/origins")
        max_batches = int(self.parameters.get("max_batches", 512))
        if len(self.batch_levels) > max_batches:
            raise ValueError(
                f"MOBER resolved {len(self.batch_levels)} batches, above max_batches="
                f"{max_batches}. Use batch_key=source or restrict the cohort."
            )
        self.target_batch = self._choose_target(frame, batches)
        batch_indices, one_hot = self._batch_encoding(batches)
        self.genes = values.shape[1]
        torch.manual_seed(self.seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(self.seed)
        self.model = _BatchVAE(
            self.genes, self.encoding_dim, len(self.batch_levels)
        ).to(self.device)
        adversary = _Adversary(self.encoding_dim, len(self.batch_levels)).to(
            self.device
        )
        learning_rate = float(self.parameters.get("learning_rate", 1e-3))
        adversary_learning_rate = float(
            self.parameters.get("adversary_learning_rate", 1e-3)
        )
        model_optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)
        adversary_optimizer = torch.optim.Adam(
            adversary.parameters(), lr=adversary_learning_rate
        )
        counts = np.bincount(batch_indices, minlength=len(self.batch_levels))
        class_weights = len(values) / np.maximum(counts, 1) / len(self.batch_levels)
        class_weights_tensor = torch.as_tensor(
            class_weights, dtype=torch.float32, device=self.device
        )
        dataset = TensorDataset(
            torch.from_numpy(values),
            torch.from_numpy(one_hot),
            torch.from_numpy(batch_indices),
        )
        batch_size = min(max(2, int(self.parameters.get("batch_size", 256))), len(values))
        generator = torch.Generator().manual_seed(self.seed)
        loader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=True,
            drop_last=len(values) % batch_size == 1,
            generator=generator,
        )
        epochs = int(self.parameters.get("epochs", 300))
        adversary_weight = float(self.parameters.get("adversary_weight", 0.01))
        kl_weight = float(self.parameters.get("kl_weight", 1e-5))
        self.history = []
        for epoch in range(epochs):
            self.model.train()
            adversary.train()
            reconstruction_total = 0.0
            adversary_total = 0.0
            batches_seen = 0
            for expression, batch_one_hot, batch_index in loader:
                expression = expression.to(self.device)
                batch_one_hot = batch_one_hot.to(self.device)
                batch_index = batch_index.to(self.device)
                decoded, latent, means, standard_deviations = self.model(
                    expression, batch_one_hot
                )

                adversary_optimizer.zero_grad(set_to_none=True)
                adversary_loss = F.nll_loss(
                    adversary(latent.detach()),
                    batch_index,
                    weight=class_weights_tensor,
                )
                adversary_loss.backward()
                adversary_optimizer.step()

                for parameter in adversary.parameters():
                    parameter.requires_grad_(False)
                model_optimizer.zero_grad(set_to_none=True)
                reconstruction = F.mse_loss(decoded, expression)
                variance = torch.square(standard_deviations)
                kl = 0.5 * torch.mean(
                    torch.square(means) + variance - torch.log(variance) - 1.0
                )
                confusion = F.nll_loss(
                    adversary(latent), batch_index, weight=class_weights_tensor
                )
                total = reconstruction + kl_weight * kl - adversary_weight * confusion
                total.backward()
                model_optimizer.step()
                for parameter in adversary.parameters():
                    parameter.requires_grad_(True)

                reconstruction_total += float(reconstruction.detach().cpu())
                adversary_total += float(adversary_loss.detach().cpu())
                batches_seen += 1
            self.history.append(
                {
                    "epoch": float(epoch + 1),
                    "reconstruction_loss": reconstruction_total / max(batches_seen, 1),
                    "adversary_loss": adversary_total / max(batches_seen, 1),
                }
            )
        return self._project(values)

    def _project(self, values: np.ndarray) -> np.ndarray:
        if self.model is None:
            raise RuntimeError("MOBER is not fitted")
        target_index = self.batch_levels.index(self.target_batch)
        batch_size = max(2, int(self.parameters.get("projection_batch_size", 1024)))
        self.model.eval()
        projected: list[np.ndarray] = []
        with torch.no_grad():
            for start in range(0, len(values), batch_size):
                expression = torch.as_tensor(
                    values[start : start + batch_size],
                    dtype=torch.float32,
                    device=self.device,
                )
                one_hot = torch.zeros(
                    (len(expression), len(self.batch_levels)),
                    dtype=torch.float32,
                    device=self.device,
                )
                one_hot[:, target_index] = 1.0
                means, _, _ = self.model.encoder(expression)
                decoded = self.model.decoder(means, one_hot)
                projected.append(decoded.cpu().numpy())
        result = np.concatenate(projected)
        if not np.isfinite(result).all():
            raise RuntimeError("MOBER projection produced non-finite expression")
        return result.astype(np.float32)

    def transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
        *,
        allow_transductive: bool,
    ) -> np.ndarray:
        del studies, metadata, allow_transductive
        values = np.asarray(matrix, dtype=np.float32)
        if values.shape[1] != self.genes:
            raise ValueError(
                f"MOBER expected {self.genes} genes, received {values.shape[1]}"
            )
        return self._project(values)

    def save(self, directory: Path) -> dict[str, object]:
        if self.model is None:
            raise RuntimeError("Cannot save an unfitted MOBER harmonizer")
        path = directory / "mober_harmonizer.pt"
        state = {
            key: value.detach().cpu() for key, value in self.model.state_dict().items()
        }
        torch.save({"model_state_dict": state}, path)
        return {**self.audit(), "artifact": path.name}

    @classmethod
    def load(cls, directory: Path, payload: dict[str, object]) -> "MoberHarmonizer":
        harmonizer = cls(
            covariates=tuple(payload.get("covariates", ())),
            parameters=dict(payload.get("parameters", {})),
            device_spec=str(
                payload.get("load_device", payload.get("device_spec", "auto"))
            ),
            seed=int(payload.get("seed", 0)),
        )
        harmonizer.batch_key = str(payload["batch_key"])
        harmonizer.batch_levels = list(map(str, payload["batch_levels"]))
        harmonizer.target_batch = str(payload["target_batch"])
        harmonizer.encoding_dim = int(payload["encoding_dim"])
        harmonizer.genes = int(payload["genes"])
        harmonizer.history = list(payload.get("history", []))
        harmonizer.model = _BatchVAE(
            harmonizer.genes,
            harmonizer.encoding_dim,
            len(harmonizer.batch_levels),
        ).to(harmonizer.device)
        state = torch.load(
            directory / str(payload["artifact"]),
            map_location=harmonizer.device,
            weights_only=True,
        )
        harmonizer.model.load_state_dict(state["model_state_dict"])
        harmonizer.model.eval()
        return harmonizer

    def audit(self) -> dict[str, object]:
        return {
            "method": self.method_id,
            "fold_behavior": "inductive_projection_to_training_target_batch",
            "batch_key": self.batch_key,
            "batch_levels": self.batch_levels,
            "target_batch": self.target_batch,
            "covariates": list(self.covariates),
            "covariates_consumed_by_model": [],
            "outcome_informed": False,
            "parameters": self.parameters,
            "seed": self.seed,
            "encoding_dim": self.encoding_dim,
            "genes": self.genes,
            "device": str(self.device),
            "device_spec": self.device_spec,
            "cuda_device_name": (
                torch.cuda.get_device_name(self.device)
                if self.device.type == "cuda"
                else ""
            ),
            "history": self.history,
            "deterministic_projection": "encoder_mean_decoded_onto_target_batch",
        }
