"""Common adapter contract, device handling, and checkpoint utilities."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
import json
from pathlib import Path
import random
from typing import Any

import numpy as np
import pandas as pd
import torch

from ..training_data import DataPartition


@dataclass
class AdapterState:
    completed_epochs: dict[str, int]
    history: list[dict[str, Any]]
    global_steps: int = 0


def resolve_device(spec: str) -> torch.device:
    if spec == "cpu":
        return torch.device("cpu")
    if spec == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("execution.device=cuda but CUDA is unavailable")
        return torch.device("cuda")
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def seed_everything(seed: int) -> None:
    random.seed(int(seed))
    np.random.seed(int(seed))
    torch.manual_seed(int(seed))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(seed))


def weighted_loader(
    partition: DataPartition,
    *,
    batch_size: int,
    seed: int,
    num_workers: int = 0,
    num_samples: int = 0,
):
    expression = torch.as_tensor(partition.matrix, dtype=torch.float32)
    categories = torch.as_tensor(partition.categories, dtype=torch.long)
    dataset = torch.utils.data.TensorDataset(expression, categories)
    generator = torch.Generator()
    generator.manual_seed(int(seed))
    weights = torch.as_tensor(partition.weights, dtype=torch.double)
    if len(weights) != len(dataset) or not torch.isfinite(weights).all() or weights.sum() <= 0:
        weights = torch.ones(len(dataset), dtype=torch.double)
    sampler = torch.utils.data.WeightedRandomSampler(
        weights,
        num_samples=int(num_samples) if int(num_samples) > 0 else len(dataset),
        replacement=True,
        generator=generator,
    )
    return torch.utils.data.DataLoader(
        dataset,
        batch_size=min(int(batch_size), len(dataset)),
        sampler=sampler,
        num_workers=int(num_workers),
        drop_last=False,
        pin_memory=torch.cuda.is_available(),
    )


class ModelAdapter(ABC):
    adapter_id: str
    supports_generation: bool = False

    def __init__(
        self,
        *,
        genes: list[str],
        cardinalities: list[int],
        covariates: tuple[str, ...],
        parameters: dict[str, Any],
        device_spec: str,
        output_dir: Path,
        checkpoint_every: int,
        resume: bool,
        seed: int,
        num_workers: int = 0,
    ) -> None:
        self.genes = list(genes)
        self.cardinalities = list(map(int, cardinalities))
        self.covariates = tuple(covariates)
        self.parameters = dict(parameters)
        self.device = resolve_device(device_spec)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.checkpoint_dir = self.output_dir / "checkpoints"
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.checkpoint_every = int(checkpoint_every)
        self.resume = bool(resume)
        self.seed = int(seed)
        self.num_workers = int(num_workers)
        self.state = AdapterState(completed_epochs={}, history=[])
        self.source_manifest: dict[str, Any] = {}
        seed_everything(self.seed)

    @property
    def checkpoint_path(self) -> Path:
        return self.checkpoint_dir / "latest.pt"

    @property
    def batch_size(self) -> int:
        return int(self.parameters.get("batch_size", 128))

    def device_summary(self) -> dict[str, object]:
        return {
            "torch_version": str(torch.__version__),
            "cuda_available": bool(torch.cuda.is_available()),
            "device": str(self.device),
            "cuda_device_name": (
                torch.cuda.get_device_name(self.device)
                if self.device.type == "cuda"
                else ""
            ),
        }

    def _common_payload(self) -> dict[str, Any]:
        return {
            "adapter_id": self.adapter_id,
            "genes": self.genes,
            "cardinalities": self.cardinalities,
            "covariates": self.covariates,
            "parameters": self.parameters,
            "seed": self.seed,
            "state": {
                "completed_epochs": self.state.completed_epochs,
                "history": self.state.history,
                "global_steps": self.state.global_steps,
            },
            "source_manifest": self.source_manifest,
            "rng_state": torch.get_rng_state(),
            "cuda_rng_state": (
                torch.cuda.get_rng_state_all() if torch.cuda.is_available() else []
            ),
        }

    def _restore_common(self, payload: dict[str, Any]) -> None:
        state = payload.get("state", {})
        self.state = AdapterState(
            completed_epochs={
                str(key): int(value)
                for key, value in state.get("completed_epochs", {}).items()
            },
            history=list(state.get("history", [])),
            global_steps=int(state.get("global_steps", 0)),
        )
        if "rng_state" in payload:
            torch.set_rng_state(payload["rng_state"].cpu())
        if torch.cuda.is_available() and payload.get("cuda_rng_state"):
            torch.cuda.set_rng_state_all(
                [state.cpu() for state in payload["cuda_rng_state"]]
            )
        self.source_manifest = dict(payload.get("source_manifest", {}))

    def _atomic_torch_save(self, payload: dict[str, Any], path: Path) -> None:
        temporary = path.with_suffix(path.suffix + ".tmp")
        torch.save(payload, temporary)
        temporary.replace(path)

    def write_history(self) -> Path:
        path = self.output_dir / "training_history.tsv"
        pd.DataFrame(self.state.history).to_csv(path, sep="\t", index=False)
        return path

    def write_adapter_summary(self) -> Path:
        path = self.output_dir / "adapter_summary.json"
        payload = {
            "adapter_id": self.adapter_id,
            "supports_generation": self.supports_generation,
            "genes": len(self.genes),
            "covariates": list(self.covariates),
            "cardinalities": self.cardinalities,
            "parameters": self.parameters,
            "device": self.device_summary(),
            "completed_epochs": self.state.completed_epochs,
            "global_steps": self.state.global_steps,
            "source_manifest": self.source_manifest,
        }
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        return path

    @abstractmethod
    def fit_stage(
        self, partition: DataPartition, *, stage: str, epochs: int, learning_rate: float
    ) -> list[dict[str, Any]]:
        raise NotImplementedError

    @abstractmethod
    def encode(self, partition: DataPartition) -> np.ndarray:
        raise NotImplementedError

    def generate(
        self, categories: np.ndarray, *, seed: int, batch_size: int | None = None
    ) -> np.ndarray:
        raise RuntimeError(f"{self.adapter_id} does not generate expression")

    @abstractmethod
    def save_final(self) -> Path:
        raise NotImplementedError
