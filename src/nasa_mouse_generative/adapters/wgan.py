"""Vinas et al. conditional WGAN-GP adapter."""

from __future__ import annotations

from pathlib import Path
from typing import Any
import math

import numpy as np
import torch

from nasa_mouse_wgan.model import ConditionalWGANGP
from nasa_mouse_wgan.training import TrainConfig, critic_features, train_epoch

from ..paper_contracts import verify_pinned_source
from ..training_data import DataPartition
from .base import ModelAdapter, weighted_loader


class WGANAdapter(ModelAdapter):
    adapter_id = "vinas_wgan_gp"
    supports_generation = True

    def __init__(
        self,
        *,
        source_path: str = "",
        validation_partition: DataPartition | None = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.source_path = str(source_path)
        self.validation_partition = validation_partition
        if self.source_path and Path(self.source_path).exists():
            self.source_manifest = verify_pinned_source(
                self.adapter_id, self.source_path
            )
            self.source_manifest["implementation"] = "pytorch_equivalent_port"
        elif bool(self.parameters.get("_paper_native", False)):
            raise FileNotFoundError(
                "Paper-native WGAN requires the pinned source checkout. Run "
                "`python -m nasa_mouse_generative prepare-upstreams`."
            )
        hidden_dims = tuple(
            int(value) for value in self.parameters.get("hidden_dims", [256, 256])
        )
        self.model_config = {
            "expression_dim": len(self.genes),
            "categorical_cardinalities": self.cardinalities,
            "noise_dim": int(self.parameters.get("noise_dim", 64)),
            "numeric_dim": int(self.parameters.get("numeric_dim", 0)),
            "hidden_dims": hidden_dims,
        }
        self.model = ConditionalWGANGP(**self.model_config).to(self.device)
        self.early_stopped_stages: set[str] = set()
        self.early_stopping_state: dict[str, dict[str, float | int]] = {}
        self._resume_payload: dict[str, Any] | None = None
        if self.resume and self.checkpoint_path.exists():
            payload = torch.load(
                self.checkpoint_path, map_location=self.device, weights_only=False
            )
            self._validate_payload(payload)
            self.model.load_state_dict(payload["model_state_dict"])
            self._restore_common(payload)
            self.early_stopped_stages = set(
                map(str, payload.get("early_stopped_stages", []))
            )
            self.early_stopping_state = {
                str(key): dict(value)
                for key, value in payload.get("early_stopping_state", {}).items()
            }
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
                self.model.generator.parameters(),
                lr=float(learning_rate),
                alpha=float(self.parameters.get("rmsprop_alpha", 0.9)),
                eps=float(self.parameters.get("rmsprop_epsilon", 1e-7)),
                momentum=0.0,
                centered=False,
            )
            optim_d = torch.optim.RMSprop(
                self.model.critic.parameters(),
                lr=float(learning_rate),
                alpha=float(self.parameters.get("rmsprop_alpha", 0.9)),
                eps=float(self.parameters.get("rmsprop_epsilon", 1e-7)),
                momentum=0.0,
                centered=False,
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
                "source_path": self.source_path,
                "early_stopped_stages": sorted(self.early_stopped_stages),
                "early_stopping_state": self.early_stopping_state,
            }
        )
        self._atomic_torch_save(payload, self.checkpoint_path)

    @staticmethod
    def _gamma_coefficient(
        real: np.ndarray, fake: np.ndarray, *, device: torch.device
    ) -> float:
        if real.shape[1] < 2 or len(real) < 3 or len(fake) < 3:
            return float("nan")
        with torch.no_grad():
            first_matrix = torch.as_tensor(real, dtype=torch.float64, device=device)
            second_matrix = torch.as_tensor(fake, dtype=torch.float64, device=device)

            def correlation(matrix: torch.Tensor) -> torch.Tensor:
                centered = matrix - matrix.mean(dim=0, keepdim=True)
                norms = torch.linalg.vector_norm(centered, dim=0)
                denominator = torch.outer(norms, norms)
                return centered.T.mm(centered) / denominator

            upper = torch.triu_indices(
                real.shape[1], real.shape[1], offset=1, device=device
            )
            first = 1.0 - correlation(first_matrix)[upper[0], upper[1]]
            second = 1.0 - correlation(second_matrix)[upper[0], upper[1]]
            finite = torch.isfinite(first) & torch.isfinite(second)
            if int(finite.sum().item()) < 2:
                return float("nan")
            first = first[finite]
            second = second[finite]
            first = first - first.mean()
            second = second - second.mean()
            denominator = torch.linalg.vector_norm(first) * torch.linalg.vector_norm(
                second
            )
            if not torch.isfinite(denominator) or float(denominator.item()) == 0.0:
                return float("nan")
            return float((torch.dot(first, second) / denominator).item())

    def _monitor_score(self) -> float:
        partition = self.validation_partition
        if partition is None or len(partition) < 3:
            return float("nan")
        generated = self.generate(partition.categories, seed=self.seed + 7919)
        maximum = int(self.parameters.get("early_stopping_max_genes", 2000))
        if maximum > 0 and partition.matrix.shape[1] > maximum:
            variances = np.var(partition.matrix, axis=0, dtype=np.float64)
            selected = np.argsort(-variances, kind="stable")[:maximum]
            real = partition.matrix[:, selected]
            generated = generated[:, selected]
        else:
            real = partition.matrix
        return self._gamma_coefficient(real, generated, device=self.device)

    def _early_stopping_checks(self) -> int:
        variant = str(
            self.parameters.get("early_stopping_variant", "released_code")
        )
        every = max(
            1,
            int(self.parameters.get("early_stopping_evaluate_every_epochs", 5)),
        )
        if variant == "paper_text":
            patience_epochs = int(
                self.parameters.get("early_stopping_patience_epochs", 30)
            )
            return max(1, math.ceil(patience_epochs / every))
        if variant != "released_code":
            raise ValueError(f"Unknown WGAN early-stopping variant: {variant}")
        return max(
            1, int(self.parameters.get("early_stopping_patience_checks", 10))
        )

    def _is_monitor_epoch(self, epoch: int) -> bool:
        first = max(1, int(self.parameters.get("early_stopping_first_epoch", 1)))
        every = max(
            1,
            int(self.parameters.get("early_stopping_evaluate_every_epochs", 5)),
        )
        return epoch >= first and (epoch - first) % every == 0

    def _paper_loader(self, partition: DataPartition):
        permutation = np.random.default_rng(self.seed).permutation(len(partition))
        expression = torch.as_tensor(
            partition.matrix[permutation], dtype=torch.float32
        )
        categories = torch.as_tensor(
            partition.categories[permutation], dtype=torch.long
        )
        dataset = torch.utils.data.TensorDataset(expression, categories)
        return torch.utils.data.DataLoader(
            dataset,
            batch_size=min(self.batch_size, len(dataset)),
            shuffle=False,
            drop_last=False,
            num_workers=self.num_workers,
            pin_memory=torch.cuda.is_available(),
        )

    def fit_stage(
        self, partition: DataPartition, *, stage: str, epochs: int, learning_rate: float
    ) -> list[dict[str, Any]]:
        if len(partition) < 2:
            raise ValueError(f"WGAN stage {stage} needs at least two profiles")
        completed = int(self.state.completed_epochs.get(stage, 0))
        if completed >= int(epochs) or stage in self.early_stopped_stages:
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
        early_stopping = bool(self.parameters.get("early_stopping", False))
        monitor_compatible = bool(
            self.validation_partition is not None
            and len(self.validation_partition) >= 3
            and (
                stage != "reference"
                or set(
                    self.validation_partition.obs.get(
                        "source", np.asarray([], dtype=str)
                    ).astype(str)
                )
                == {"archs4"}
            )
        )
        state = self.early_stopping_state.setdefault(
            stage, {"best_score": float("-inf"), "checks_without_improvement": 0}
        )
        patience_checks = self._early_stopping_checks() if early_stopping else 0
        progress_every = max(
            1, int(self.parameters.get("progress_every_epochs", 10))
        )
        best_path = self.checkpoint_dir / f"best_{stage}.pt"
        for epoch in range(completed + 1, int(epochs) + 1):
            if bool(self.parameters.get("weighted_sampling", True)):
                loader = weighted_loader(
                    partition,
                    batch_size=self.batch_size,
                    seed=self.seed
                    + epoch
                    + 10000 * len(self.state.completed_epochs),
                    num_workers=self.num_workers,
                )
            else:
                loader = self._paper_loader(partition)
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
            if (
                epoch == completed + 1
                or epoch % progress_every == 0
                or epoch == int(epochs)
            ):
                print(
                    f"[wgan] stage={stage} epoch={epoch}/{epochs} "
                    f"critic={float(metrics['critic_loss']):.6f} "
                    f"generator={float(metrics['generator_loss']):.6f} "
                    f"wasserstein={float(metrics['wasserstein_estimate']):.6f} "
                    f"gradient_penalty={float(metrics['gradient_penalty']):.6f}",
                    flush=True,
                )
            should_stop = False
            if early_stopping and monitor_compatible and self._is_monitor_epoch(epoch):
                score = self._monitor_score()
                row["early_stopping_score"] = score
                if np.isfinite(score) and score > float(state["best_score"]):
                    state["best_score"] = float(score)
                    state["checks_without_improvement"] = 0
                    self._atomic_torch_save(
                        {
                            "model_state_dict": self.model.state_dict(),
                            "stage": stage,
                            "epoch": epoch,
                            "score": score,
                        },
                        best_path,
                    )
                else:
                    state["checks_without_improvement"] = int(
                        state["checks_without_improvement"]
                    ) + 1
                row["early_stopping_checks_without_improvement"] = int(
                    state["checks_without_improvement"]
                )
                self.write_history()
                print(
                    f"[wgan] monitor stage={stage} epoch={epoch} "
                    f"gamma={float(score):.6f} "
                    f"best={float(state['best_score']):.6f} "
                    f"checks_without_improvement="
                    f"{int(state['checks_without_improvement'])}/"
                    f"{patience_checks}",
                    flush=True,
                )
                should_stop = (
                    int(state["checks_without_improvement"]) >= patience_checks
                )
            if epoch % self.checkpoint_every == 0 or epoch == int(epochs):
                self._save_checkpoint(stage, optim_g, optim_d)
                self.write_history()
            if should_stop:
                self.early_stopped_stages.add(stage)
                if best_path.exists():
                    best = torch.load(
                        best_path, map_location=self.device, weights_only=False
                    )
                    self.model.load_state_dict(best["model_state_dict"])
                self._save_checkpoint(stage, optim_g, optim_d)
                self.write_history()
                print(
                    f"[wgan] early stop stage={stage} epoch={epoch} "
                    f"variant={self.parameters.get('early_stopping_variant')} "
                    f"best_gamma={float(state['best_score']):.6f}",
                    flush=True,
                )
                break
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
                "source_path": self.source_path,
                "early_stopped_stages": sorted(self.early_stopped_stages),
                "early_stopping_state": self.early_stopping_state,
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
            source_path=str(payload.get("source_path", "")),
            validation_partition=None,
        )
        adapter.model.load_state_dict(payload["model_state_dict"])
        adapter._restore_common(payload)
        adapter.model.eval()
        return adapter
