"""Executable adapters for the three benchmark model families."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..config import BenchmarkConfig
from ..training_data import PreparedTrainingData
from .base import ModelAdapter


def create_adapter(
    config: BenchmarkConfig,
    data: PreparedTrainingData,
    parameters: dict[str, Any],
    output_dir: str | Path,
) -> ModelAdapter:
    common = {
        "genes": data.genes,
        "cardinalities": data.encoder.cardinalities,
        "covariates": data.covariates,
        "parameters": parameters,
        "device_spec": config.execution.device,
        "output_dir": Path(output_dir),
        "checkpoint_every": config.execution.checkpoint_every_epochs,
        "num_workers": int(
            parameters.get("num_workers", config.execution.num_workers)
        ),
        "resume": config.execution.resume,
        "seed": config.training.seed,
    }
    if config.training.model == "vinas_wgan_gp":
        from .wgan import WGANAdapter

        return WGANAdapter(
            **common,
            source_path=config.execution.wgan_source,
            validation_partition=data.partitions.get("validation"),
        )
    if config.training.model == "lacan_diffusion":
        from .diffusion import DiffusionAdapter

        fit_partition = data.reference if data.reference is not None else data.train
        return DiffusionAdapter(
            **common,
            reconstruction_matrix=fit_partition.matrix,
            l1000_map=config.features.l1000_map,
            source_path=config.execution.diffusion_source,
        )
    if config.training.model == "genejepa":
        from .genejepa import GeneJEPAAdapter

        return GeneJEPAAdapter(
            **common,
            source_path=config.execution.genejepa_source,
        )
    raise ValueError(f"Unsupported model adapter: {config.training.model}")


def load_adapter(output_dir: str | Path, *, device_spec: str = "auto") -> ModelAdapter:
    import torch

    root = Path(output_dir)
    payload = torch.load(root / "model.pt", map_location="cpu", weights_only=False)
    adapter_id = payload.get("adapter_id")
    if adapter_id == "vinas_wgan_gp":
        from .wgan import WGANAdapter

        return WGANAdapter.load(root, device_spec=device_spec)
    if adapter_id == "lacan_diffusion":
        from .diffusion import DiffusionAdapter

        return DiffusionAdapter.load(root, device_spec=device_spec)
    if adapter_id == "genejepa":
        from .genejepa import GeneJEPAAdapter

        return GeneJEPAAdapter.load(root, device_spec=device_spec)
    raise ValueError(f"Unknown adapter in {root / 'model.pt'}: {adapter_id!r}")


__all__ = ["ModelAdapter", "create_adapter", "load_adapter"]
