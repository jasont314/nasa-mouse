"""Contrastive FLT/GC guidance for a trained conditional ModelDDIM."""

from __future__ import annotations

import json
from pathlib import Path
import time
from typing import Iterable

import numpy as np
import pandas as pd
import torch

from .conditional_config import load_conditional_config
from .conditional_data import load_conditional_prepared
from .evaluate import _load_model
from .upstream import compute_alpha, quadratic_beta_schedule, verify_source


def opposite_condition_indices(
    classes: Iterable[str], class_indices: np.ndarray
) -> np.ndarray:
    """Map each joint tissue/condition class to its opposite condition class."""

    labels = list(map(str, classes))
    index_by_label = {label: index for index, label in enumerate(labels)}
    opposite: dict[int, int] = {}
    for index, label in enumerate(labels):
        if "condition=flight" in label:
            target = label.replace("condition=flight", "condition=ground_control")
        elif "condition=ground_control" in label:
            target = label.replace("condition=ground_control", "condition=flight")
        else:
            continue
        if target in index_by_label:
            opposite[index] = index_by_label[target]
    requested = np.asarray(class_indices, dtype=np.int64)
    missing = sorted(set(requested.tolist()) - set(opposite))
    if missing:
        names = [labels[index] for index in missing[:5]]
        raise ValueError(f"Classes lack an opposite FLT/GC condition: {names}")
    return np.asarray([opposite[int(index)] for index in requested], dtype=np.int64)


def contrastive_ddim_final(
    initial_noise: torch.Tensor,
    target_labels: torch.Tensor,
    opposite_labels: torch.Tensor,
    model: torch.nn.Module,
    betas: torch.Tensor,
    *,
    sequence: Iterable[int],
    guidance_scale: float,
    eta: float = 0.0,
    generator: torch.Generator | None = None,
) -> torch.Tensor:
    """Denoise while extrapolating from opposite toward requested condition."""

    scale = float(guidance_scale)
    if not np.isfinite(scale) or scale < 0:
        raise ValueError("guidance_scale must be finite and nonnegative")
    if target_labels.shape != opposite_labels.shape:
        raise ValueError("Target and opposite labels must have matching shapes")
    sequence = list(map(int, sequence))
    sequence_next = [-1] + sequence[:-1]
    current = initial_noise
    with torch.no_grad():
        for current_t, next_t in zip(reversed(sequence), reversed(sequence_next)):
            timestep = torch.full(
                (len(current),), current_t, dtype=torch.long, device=current.device
            )
            next_timestep = torch.full(
                (len(current),), next_t, dtype=torch.long, device=current.device
            )
            alpha = compute_alpha(betas, timestep)
            alpha_next = compute_alpha(betas, next_timestep)
            target_noise = model(current, timestep, target_labels)
            if scale == 1.0:
                guided_noise = target_noise
            else:
                opposite_noise = model(current, timestep, opposite_labels)
                guided_noise = opposite_noise + scale * (
                    target_noise - opposite_noise
                )
            predicted_clean = (
                current - guided_noise * (1.0 - alpha).sqrt()
            ) / alpha.sqrt()
            stochastic_scale = eta * (
                (1.0 - alpha / alpha_next)
                * (1.0 - alpha_next)
                / (1.0 - alpha)
            ).sqrt()
            residual_scale = ((1.0 - alpha_next) - stochastic_scale**2).sqrt()
            if eta:
                stochastic_noise = torch.randn(
                    current.shape,
                    generator=generator,
                    device=current.device,
                    dtype=current.dtype,
                )
            else:
                stochastic_noise = 0.0
            current = (
                alpha_next.sqrt() * predicted_clean
                + stochastic_scale * stochastic_noise
                + residual_scale * guided_noise
            )
    return current.detach().cpu()


def _scale_label(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace(".", "p")


def generate_contrastive_training(
    config_path: str | Path,
    *,
    guidance_scales: Iterable[float],
    seeds: Iterable[int],
    batch_size: int = 256,
    eta: float = 0.0,
) -> Path:
    """Generate train-condition profiles for classifier augmentation."""

    config = load_conditional_config(config_path)
    verify_source(config["run"]["source_root"])
    prepared = load_conditional_prepared(config["data"]["prepared_h5"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("Contrastive ModelDDIM generation requires CUDA")
    model, payload, model_path = _load_model(config, prepared, device)
    class_indices = prepared["train"]["class_index"]
    opposite_indices = opposite_condition_indices(
        prepared["classes"], class_indices
    )
    betas = quadratic_beta_schedule(
        beta_start=float(config["model"]["beta_start"]),
        beta_end=float(config["model"]["beta_end"]),
        timesteps=int(config["model"]["diffusion_timesteps"]),
    ).to(device)
    output = Path(config["run"]["output_dir"]) / "evaluation" / "contrastive_training"
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for scale in tuple(map(float, guidance_scales)):
        for seed in tuple(map(int, seeds)):
            path = output / f"scale_{_scale_label(scale)}_seed_{seed}.npz"
            if path.exists():
                rows.append(
                    {
                        "guidance_scale": scale,
                        "seed": seed,
                        "profiles": len(class_indices),
                        "seconds": 0.0,
                        "path": str(path),
                        "status": "existing",
                    }
                )
                continue
            generator = torch.Generator(device=device).manual_seed(seed)
            collected: list[np.ndarray] = []
            started = time.time()
            for start in range(0, len(class_indices), int(batch_size)):
                end = min(start + int(batch_size), len(class_indices))
                target_index = torch.as_tensor(
                    class_indices[start:end], dtype=torch.long, device=device
                )
                opposite_index = torch.as_tensor(
                    opposite_indices[start:end], dtype=torch.long, device=device
                )
                target = torch.nn.functional.one_hot(
                    target_index, num_classes=len(prepared["classes"])
                ).long()
                opposite = torch.nn.functional.one_hot(
                    opposite_index, num_classes=len(prepared["classes"])
                ).long()
                noise = torch.randn(
                    (end - start, len(prepared["genes"])),
                    generator=generator,
                    device=device,
                )
                generated = contrastive_ddim_final(
                    noise,
                    target,
                    opposite,
                    model,
                    betas,
                    sequence=range(int(config["model"]["diffusion_timesteps"])),
                    guidance_scale=scale,
                    eta=float(eta),
                    generator=generator,
                )
                collected.append(generated.numpy().astype(np.float32))
            matrix = np.concatenate(collected)
            np.savez_compressed(
                path,
                scaled_expression=matrix,
                class_index=class_indices,
                opposite_class_index=opposite_indices,
                source_row=prepared["train"]["source_row"],
                genes=np.asarray(prepared["genes"]),
                guidance_scale=scale,
                sampling_seed=seed,
                eta=float(eta),
            )
            rows.append(
                {
                    "guidance_scale": scale,
                    "seed": seed,
                    "profiles": len(matrix),
                    "seconds": float(time.time() - started),
                    "path": str(path),
                    "status": "generated",
                }
            )
            print(
                f"[contrastive-ddim] scale={scale:g} seed={seed} "
                f"profiles={len(matrix)}",
                flush=True,
            )
    table = pd.DataFrame(rows)
    table.to_csv(output / "generation_manifest.tsv", sep="\t", index=False)
    summary = {
        "status": "complete",
        "model": str(model_path),
        "model_epoch": int(payload.get("epoch", 0)),
        "prepared_data": str(config["data"]["prepared_h5"]),
        "split": "train",
        "profiles": int(len(class_indices)),
        "guidance_scales": sorted(set(table["guidance_scale"].astype(float))),
        "seeds": sorted(set(table["seed"].astype(int))),
        "eta": float(eta),
        "device": torch.cuda.get_device_name(device),
        "interpretation": (
            "Scale 1 is ordinary conditional DDIM. Larger scales extrapolate "
            "the learned requested-condition prediction away from the opposite "
            "condition for the same tissue."
        ),
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary_path
