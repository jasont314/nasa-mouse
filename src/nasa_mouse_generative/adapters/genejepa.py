"""Official GeneJEPA representation adapter for exploratory bulk transfer."""

from __future__ import annotations

import importlib
from pathlib import Path
import sys
from typing import Any

import numpy as np
import torch
import torch.nn.functional as functional

from ..training_data import DataPartition
from .base import ModelAdapter, weighted_loader


PINNED_COMMIT = "a2f4d7218b17f2f52cc5f1cc94420c8ef1ae3265"


def _official_classes(source_path: str | Path):
    source = Path(source_path)
    if not (source / "genejepa" / "models.py").exists():
        raise FileNotFoundError(
            f"Pinned GeneJEPA source not found at {source}. Run "
            "`python -m nasa_mouse_generative prepare-upstreams`."
        )
    source_string = str(source.resolve())
    if source_string not in sys.path:
        sys.path.insert(0, source_string)
    configs = importlib.import_module("genejepa.configs")
    models = importlib.import_module("genejepa.models")
    return configs.ModelConfig, models.GenePerceiverJEPA


def _vicreg(values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    values = values.float()
    std = torch.sqrt(values.var(dim=0, unbiased=False) + 1e-4)
    variance = functional.relu(1.0 - std).mean()
    if len(values) < 2:
        covariance = torch.zeros((), device=values.device)
    else:
        centered = values - values.mean(dim=0)
        covariance_matrix = centered.T @ centered / (len(values) - 1)
        covariance = (
            covariance_matrix.square().sum()
            - covariance_matrix.diagonal().square().sum()
        ) / values.shape[1]
    return variance, covariance


class GeneJEPAAdapter(ModelAdapter):
    adapter_id = "genejepa"
    supports_generation = False

    def __init__(self, *, source_path: str, **kwargs) -> None:
        super().__init__(**kwargs)
        self.source_path = str(source_path)
        ModelConfig, GenePerceiverJEPA = _official_classes(self.source_path)
        self.model_config = {
            "d": int(self.parameters.get("d", 128)),
            "latents_L": int(self.parameters.get("latents_L", 64)),
            "blocks_D": int(self.parameters.get("blocks_D", 4)),
            "heads_h": int(self.parameters.get("heads_h", 4)),
            "cross_attn_chunk_size": int(
                self.parameters.get("cross_attn_chunk_size", 32)
            ),
            "gene_vocab_size": len(self.genes),
            "mask_ratio": float(self.parameters.get("mask_ratio", 0.45)),
            "num_targets": int(self.parameters.get("num_targets", 1)),
            "min_context_genes": int(
                self.parameters.get("min_context_genes", 128)
            ),
            "min_target_genes_per_block": int(
                self.parameters.get("min_target_genes_per_block", 16)
            ),
            "ema_start_decay": float(
                self.parameters.get("ema_start_decay", 0.99)
            ),
            "ema_end_decay": float(self.parameters.get("ema_end_decay", 0.999)),
            "ema_warmup_steps": int(self.parameters.get("ema_warmup_steps", 0)),
            "identity_value_split_ratio": float(
                self.parameters.get("identity_value_split_ratio", 0.5)
            ),
            "fourier_num_frequencies": int(
                self.parameters.get("fourier_num_frequencies", 64)
            ),
            "fourier_min_freq": float(
                self.parameters.get("fourier_min_freq", 0.1)
            ),
            "fourier_max_freq": float(
                self.parameters.get("fourier_max_freq", 100.0)
            ),
            "fourier_freq_scale": float(
                self.parameters.get("fourier_freq_scale", 1.0)
            ),
            "predictor_depth": int(self.parameters.get("predictor_depth", 3)),
            "predictor_expansion_factor": int(
                self.parameters.get("predictor_expansion_factor", 4)
            ),
        }
        self.max_tokens = min(
            int(self.parameters.get("max_tokens", 2048)), len(self.genes)
        )
        minimum = (
            self.model_config["min_context_genes"]
            + self.model_config["num_targets"]
            * self.model_config["min_target_genes_per_block"]
        )
        if self.max_tokens < minimum:
            raise ValueError(
                f"GeneJEPA max_tokens={self.max_tokens} is below its required "
                f"context/target token count {minimum}"
            )
        if self.model_config["d"] % self.model_config["heads_h"]:
            raise ValueError("GeneJEPA d must be divisible by heads_h")
        self.model = GenePerceiverJEPA(ModelConfig(**self.model_config)).to(self.device)
        self._hard_sync_teacher()
        self._resume_payload: dict[str, Any] | None = None
        if self.resume and self.checkpoint_path.exists():
            payload = torch.load(
                self.checkpoint_path, map_location=self.device, weights_only=False
            )
            self._validate_payload(payload)
            self.model.load_state_dict(payload["model_state_dict"])
            self._restore_common(payload)
            self._resume_payload = payload

    def _hard_sync_teacher(self) -> None:
        teacher = self.model.teacher_encoder.ema_model
        teacher.load_state_dict(self.model.student_encoder.state_dict())
        teacher.eval()
        for parameter in teacher.parameters():
            parameter.requires_grad = False

    def _validate_payload(self, payload: dict[str, Any]) -> None:
        if payload.get("adapter_id") != self.adapter_id:
            raise ValueError("Checkpoint belongs to another adapter")
        if payload.get("genes") != self.genes:
            raise ValueError("Checkpoint gene order differs from this run")

    def _ragged(self, matrix: torch.Tensor) -> dict[str, torch.Tensor]:
        indices: list[torch.Tensor] = []
        values: list[torch.Tensor] = []
        lengths: list[int] = []
        for row in matrix:
            finite = torch.isfinite(row)
            available = torch.nonzero(finite & row.ne(0), as_tuple=False).flatten()
            if len(available) < self.max_tokens:
                available = torch.nonzero(finite, as_tuple=False).flatten()
            if len(available) > self.max_tokens:
                ranking = torch.topk(
                    row[available].abs(), self.max_tokens, sorted=False
                ).indices
                available = available[ranking]
            available = available.sort().values
            indices.append(available.long())
            values.append(row[available].float())
            lengths.append(len(available))
        flat_indices = torch.cat(indices)
        flat_values = torch.cat(values)
        offsets = functional.pad(
            torch.as_tensor(lengths, device=matrix.device, dtype=torch.long).cumsum(0),
            (1, 0),
        )
        return {"indices": flat_indices, "values": flat_values, "offsets": offsets}

    def _save_checkpoint(self, stage: str, optimizer) -> None:
        payload = self._common_payload()
        payload.update(
            {
                "model_config": self.model_config,
                "model_state_dict": self.model.state_dict(),
                "active_stage": stage,
                "optimizer_state_dict": optimizer.state_dict(),
                "source_path": self.source_path,
                "source_commit": PINNED_COMMIT,
                "max_tokens": self.max_tokens,
            }
        )
        self._atomic_torch_save(payload, self.checkpoint_path)

    def fit_stage(
        self, partition: DataPartition, *, stage: str, epochs: int, learning_rate: float
    ) -> list[dict[str, Any]]:
        if len(partition) < 2:
            raise ValueError(f"GeneJEPA stage {stage} needs at least two profiles")
        completed = int(self.state.completed_epochs.get(stage, 0))
        if completed >= int(epochs):
            return [row for row in self.state.history if row.get("stage") == stage]
        optimizer = torch.optim.AdamW(
            [parameter for parameter in self.model.parameters() if parameter.requires_grad],
            lr=float(learning_rate),
            weight_decay=float(self.parameters.get("weight_decay", 2e-4)),
            betas=(0.9, 0.98),
        )
        if self._resume_payload and self._resume_payload.get("active_stage") == stage:
            optimizer.load_state_dict(self._resume_payload["optimizer_state_dict"])
        sim_coeff = float(self.parameters.get("sim_coeff", 1.0))
        var_coeff = float(self.parameters.get("var_coeff", 25.0))
        cov_coeff = float(self.parameters.get("cov_coeff", 1.0))
        for epoch in range(completed + 1, int(epochs) + 1):
            loader = weighted_loader(
                partition,
                batch_size=self.batch_size,
                seed=self.seed + epoch + 10000 * len(self.state.completed_epochs),
                num_workers=self.num_workers,
            )
            losses: list[float] = []
            similarities: list[float] = []
            self.model.train()
            for expression, _ in loader:
                expression = expression.to(self.device)
                optimizer.zero_grad(set_to_none=True)
                predicted, target, student = self.model(**self._ragged(expression))
                if predicted.numel() == 0:
                    continue
                similarity = 1.0 - (
                    functional.normalize(predicted.float(), dim=1)
                    * functional.normalize(target.detach().float(), dim=1)
                ).sum(dim=1).mean()
                variance, covariance = _vicreg(predicted)
                student_variance, student_covariance = _vicreg(student)
                loss = (
                    sim_coeff * similarity
                    + var_coeff * variance
                    + cov_coeff * covariance
                    + 20.0 * student_variance
                    + student_covariance
                )
                if not torch.isfinite(loss):
                    raise FloatingPointError(
                        f"Non-finite GeneJEPA loss at {stage} epoch {epoch}"
                    )
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
                optimizer.step()
                progress = self.state.global_steps / max(1, int(epochs) * len(loader))
                start_decay = float(self.parameters.get("ema_start_decay", 0.99))
                end_decay = float(self.parameters.get("ema_end_decay", 0.999))
                self.model.teacher_encoder.beta = start_decay + (
                    end_decay - start_decay
                ) * min(progress, 1.0)
                if self.state.global_steps >= int(
                    self.parameters.get("ema_warmup_steps", 0)
                ):
                    self.model.update_teacher()
                self.state.global_steps += 1
                losses.append(float(loss.detach().cpu()))
                similarities.append(float(similarity.detach().cpu()))
            if not losses:
                raise RuntimeError(
                    "GeneJEPA produced no valid masked batches; increase selected genes "
                    "or lower min_context_genes"
                )
            row = {
                "stage": stage,
                "epoch": epoch,
                "learning_rate": float(learning_rate),
                "loss": float(np.mean(losses)),
                "cosine_loss": float(np.mean(similarities)),
            }
            self.state.history.append(row)
            self.state.completed_epochs[stage] = epoch
            if epoch % self.checkpoint_every == 0 or epoch == int(epochs):
                self._save_checkpoint(stage, optimizer)
                self.write_history()
        self._resume_payload = None
        return [row for row in self.state.history if row.get("stage") == stage]

    def encode(self, partition: DataPartition) -> np.ndarray:
        outputs = []
        self.model.eval()
        with torch.no_grad():
            for start in range(0, len(partition), self.batch_size):
                expression = torch.as_tensor(
                    partition.matrix[start : start + self.batch_size],
                    dtype=torch.float32,
                    device=self.device,
                )
                outputs.append(
                    self.model.get_embedding(**self._ragged(expression), use_teacher=True)
                    .detach()
                    .cpu()
                    .numpy()
                )
        return np.concatenate(outputs).astype(np.float32)

    def save_final(self) -> Path:
        path = self.output_dir / "model.pt"
        payload = self._common_payload()
        payload.update(
            {
                "model_config": self.model_config,
                "model_state_dict": self.model.state_dict(),
                "source_path": self.source_path,
                "source_commit": PINNED_COMMIT,
                "max_tokens": self.max_tokens,
            }
        )
        self._atomic_torch_save(payload, path)
        self.write_history()
        self.write_adapter_summary()
        return path

    @classmethod
    def load(
        cls, output_dir: Path, *, device_spec: str = "auto"
    ) -> "GeneJEPAAdapter":
        payload = torch.load(
            Path(output_dir) / "model.pt", map_location="cpu", weights_only=False
        )
        parameters = dict(payload["parameters"])
        parameters["max_tokens"] = int(payload["max_tokens"])
        adapter = cls(
            genes=list(payload["genes"]),
            cardinalities=list(payload["cardinalities"]),
            covariates=tuple(payload["covariates"]),
            parameters=parameters,
            device_spec=device_spec,
            output_dir=Path(output_dir),
            checkpoint_every=1,
            resume=False,
            seed=int(payload["seed"]),
            num_workers=0,
            source_path=str(payload["source_path"]),
        )
        adapter.model.load_state_dict(payload["model_state_dict"])
        adapter._restore_common(payload)
        adapter.model.eval()
        return adapter
