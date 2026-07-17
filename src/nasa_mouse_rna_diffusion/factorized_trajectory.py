"""Validation-only denoising trajectory for a factorized ModelDDIM adapter."""

from __future__ import annotations

import json
from pathlib import Path
import time
from typing import Iterable

import numpy as np
import pandas as pd
import torch

from .evaluate import _plot_trajectory
from .factorized_adapter import encode_factorized_labels, load_factorized_role
from .factorized_config import load_factorized_config
from .factorized_evaluate import (
    _GuidedModel,
    _balanced_indices,
    _load_adapter_model,
)
from .upstream import ddim_trajectory, quadratic_beta_schedule


def plot_factorized_trajectory(
    config_path: str | Path,
    *,
    sample_count: int = 128,
    batch_size: int = 64,
    snapshot_timesteps: Iterable[int] = (1000, 200, 0),
    sampling_seed: int | None = None,
    model_artifact: str = "model.pt",
) -> Path:
    """Plot exact DDIM states without loading the locked test role."""

    if sample_count <= 0 or batch_size <= 0:
        raise ValueError("sample_count and batch_size must be positive")
    if Path(model_artifact).name != model_artifact:
        raise ValueError("model_artifact must be a filename under the run directory")

    config = load_factorized_config(config_path)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("Factorized paper ModelDDIM trajectory generation requires CUDA")

    model, schema, payload = _load_adapter_model(
        config, device, artifact_name=model_artifact
    )
    data_options = config["data"]
    validation = load_factorized_role(
        data_options["prepared_h5"], data_options["samples_tsv"], "validation"
    )
    seed = (
        int(config["run"]["seed"]) + 3000
        if sampling_seed is None
        else int(sampling_seed)
    )
    selected = _balanced_indices(
        validation["samples"], min(int(sample_count), len(validation["samples"])), seed
    )
    samples = validation["samples"].iloc[selected].reset_index(drop=True)
    labels = encode_factorized_labels(samples, schema)
    timesteps = int(config["model"]["diffusion_timesteps"])
    requested = tuple(dict.fromkeys(map(int, snapshot_timesteps)))
    invalid = [value for value in requested if value < 0 or value > timesteps]
    if invalid:
        raise ValueError(f"snapshot timesteps must be in [0, {timesteps}]: {invalid}")

    betas = quadratic_beta_schedule(
        beta_start=float(config["model"]["beta_start"]),
        beta_end=float(config["model"]["beta_end"]),
        timesteps=timesteps,
    ).to(device)
    guided = _GuidedModel(model, schema, scale=1.0)
    generator = torch.Generator(device=device).manual_seed(seed)
    collected: dict[int, list[np.ndarray]] = {value: [] for value in requested}
    started = time.time()
    genes = len(validation["genes"])
    for start in range(0, len(labels), int(batch_size)):
        end = min(start + int(batch_size), len(labels))
        condition = torch.as_tensor(labels[start:end], device=device)
        noise = torch.randn(
            (end - start, genes), generator=generator, device=device
        )
        trajectory = ddim_trajectory(
            noise,
            condition,
            guided,
            betas,
            sequence=range(timesteps),
            snapshot_timesteps=requested,
            eta=0.0,
            generator=generator,
        )
        for timestep, expression in trajectory.items():
            collected[timestep].append(expression.numpy().astype(np.float32))
        print(
            f"[factorized-ddim:trajectory] sampled={end}/{len(labels)}",
            flush=True,
        )
    snapshots = {
        timestep: np.concatenate(values) for timestep, values in collected.items()
    }

    output = (
        Path(config["run"]["output_dir"])
        / "evaluation"
        / "validation_denoising_trajectory"
    )
    output.mkdir(parents=True, exist_ok=True)
    tissues = sorted(samples["tissue"].astype(str).unique())
    tissue_to_index = {value: index for index, value in enumerate(tissues)}
    tissue_labels = samples["tissue"].astype(str).map(tissue_to_index).to_numpy()
    figure_path, coordinates, pca = _plot_trajectory(
        snapshots=snapshots,
        labels=tissue_labels,
        classes=tissues,
        real_background=validation["expression"],
        output=output,
        seed=seed,
        title="Study-conditioned DDIM denoising across OSDR mouse tissues",
        background_label="real OSDR validation",
        filename="osdr_factorized_ddim_trajectory_pca.png",
    )

    np.savez_compressed(
        output / "trajectory_scaled_expression.npz",
        source_row=np.asarray(validation["source_row"])[selected],
        tissue=samples["tissue"].to_numpy(dtype=str),
        condition=samples["condition"].to_numpy(dtype=str),
        genes=np.asarray(validation["genes"]),
        **{f"t{timestep}": values for timestep, values in snapshots.items()},
    )
    metadata_columns = [
        column
        for column in (
            "_row_index",
            "accession",
            "tissue",
            "condition",
            "material_type",
            "sex",
            "muscle_group",
        )
        if column in samples.columns
    ]
    for timestep, values in coordinates.items():
        table = samples[metadata_columns].copy()
        table.insert(0, "PC2", values[:, 1])
        table.insert(0, "PC1", values[:, 0])
        table.to_csv(
            output / f"trajectory_pca_t{timestep}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    summary = {
        "status": "complete",
        "model": str(Path(config["run"]["output_dir"]) / model_artifact),
        "model_format": payload.get("format"),
        "split": "validation",
        "test_loaded": False,
        "profiles": int(len(samples)),
        "validation_background_profiles": int(len(validation["expression"])),
        "genes": int(genes),
        "tissues": tissues,
        "sampling_seed": seed,
        "sampling_steps": timesteps,
        "snapshot_timesteps": list(requested),
        "guidance_scale": 1.0,
        "eta": 0.0,
        "posthoc_calibration_applied": False,
        "posthoc_calibration_note": (
            "The accepted train-only distribution calibration is an endpoint "
            "operation and is not a DDIM denoising step."
        ),
        "pca_fit": "all real OSDR validation profiles in model-scaled expression space",
        "pca_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "device": str(device),
        "gpu": torch.cuda.get_device_name(device),
        "elapsed_seconds": time.time() - started,
        "figure": str(figure_path),
    }
    with (output / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    return output
