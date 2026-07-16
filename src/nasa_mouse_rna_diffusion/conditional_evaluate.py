"""Held-out fidelity and FLT/GC evaluation for conditional upstream ModelDDIM."""

from __future__ import annotations

import json
from pathlib import Path
import time
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import torch

from nasa_mouse_diffusion.evaluate import generated_quality
from nasa_mouse_generative.metrics import (
    _condition_effect,
    classifier_utility,
    fidelity_selection,
    memorization_metrics,
)

from .conditional_config import load_conditional_config
from .conditional_data import load_conditional_prepared
from .evaluate import _classifier_predictions, _correlation, _load_model
from .upstream import ddim_trajectory, quadratic_beta_schedule, verify_source


def _balanced_indices(labels: np.ndarray, limit: int, seed: int) -> np.ndarray:
    labels = np.asarray(labels)
    if limit <= 0 or len(labels) <= limit:
        return np.arange(len(labels), dtype=int)
    rng = np.random.default_rng(seed)
    groups = {
        label: rng.permutation(np.flatnonzero(labels == label))
        for label in sorted(set(labels.tolist()))
    }
    selected: list[int] = []
    offset = 0
    while len(selected) < limit:
        added = False
        for label in groups:
            if offset < len(groups[label]):
                selected.append(int(groups[label][offset]))
                added = True
                if len(selected) == limit:
                    break
        if not added:
            break
        offset += 1
    return np.asarray(selected, dtype=int)


def _sample(
    *,
    model,
    class_indices: np.ndarray,
    genes: int,
    classes: int,
    betas: torch.Tensor,
    timesteps: int,
    batch_size: int,
    seed: int,
    device: torch.device,
) -> np.ndarray:
    generator = torch.Generator(device=device).manual_seed(int(seed))
    collected: list[np.ndarray] = []
    for start in range(0, len(class_indices), int(batch_size)):
        end = min(start + int(batch_size), len(class_indices))
        class_batch = torch.as_tensor(
            class_indices[start:end], dtype=torch.long, device=device
        )
        labels = torch.nn.functional.one_hot(
            class_batch, num_classes=int(classes)
        ).long()
        noise = torch.randn(
            (end - start, int(genes)), generator=generator, device=device
        )
        generated = ddim_trajectory(
            noise,
            labels,
            model,
            betas,
            sequence=range(int(timesteps)),
            snapshot_timesteps=(0,),
            eta=0.0,
        )[0]
        collected.append(generated.numpy().astype(np.float32))
        print(
            f"[conditional-ddim:evaluate] sampled {end}/{len(class_indices)}",
            flush=True,
        )
    return np.concatenate(collected)


def _plot_real_synthetic_pca(
    real: np.ndarray,
    synthetic: np.ndarray,
    samples: pd.DataFrame,
    output: Path,
) -> Path:
    coordinates = PCA(n_components=2, random_state=0).fit_transform(
        np.concatenate([real, synthetic])
    )
    real_coordinates = coordinates[: len(real)]
    synthetic_coordinates = coordinates[len(real) :]
    colors = {"flight": "#C14924", "ground_control": "#176B87"}
    figure, axes = plt.subplots(1, 2, figsize=(12.2, 5.1))
    for condition in ("ground_control", "flight"):
        mask = samples["condition"].astype(str).eq(condition).to_numpy()
        axes[0].scatter(
            real_coordinates[mask, 0],
            real_coordinates[mask, 1],
            s=25,
            color=colors[condition],
            marker="o",
            alpha=0.65,
            edgecolors="none",
            label=f"real {condition}",
        )
        axes[0].scatter(
            synthetic_coordinates[mask, 0],
            synthetic_coordinates[mask, 1],
            s=27,
            color=colors[condition],
            marker="x",
            alpha=0.75,
            label=f"synthetic {condition}",
        )
    axes[0].set_title("Condition and source")
    axes[0].legend(frameon=False, fontsize=8)
    tissues = sorted(samples["tissue"].astype(str).unique())
    palette = plt.get_cmap("tab20", len(tissues))
    for index, tissue in enumerate(tissues):
        mask = samples["tissue"].astype(str).eq(tissue).to_numpy()
        axes[1].scatter(
            real_coordinates[mask, 0],
            real_coordinates[mask, 1],
            s=20,
            color=palette(index),
            marker="o",
            alpha=0.35,
            edgecolors="none",
        )
        axes[1].scatter(
            synthetic_coordinates[mask, 0],
            synthetic_coordinates[mask, 1],
            s=24,
            color=palette(index),
            marker="x",
            alpha=0.75,
            label=tissue,
        )
    axes[1].set_title("Expected tissue class")
    axes[1].legend(frameon=False, fontsize=6, ncol=2, bbox_to_anchor=(1.01, 1.0))
    for axis in axes:
        axis.set_xlabel("PCA 1")
        axis.set_ylabel("PCA 2")
        axis.grid(alpha=0.16)
    figure.suptitle("Held-out OSDR real and conditional synthetic expression", fontweight="bold")
    figure.tight_layout()
    path = output / "real_vs_synthetic_pca.png"
    figure.savefig(path, dpi=220, bbox_inches="tight")
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return path


def _read_pathways(path: str | Path, genes: list[str]) -> dict[str, np.ndarray]:
    gene_map = {gene: index for index, gene in enumerate(genes)}
    result: dict[str, np.ndarray] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            indices = sorted({gene_map[gene] for gene in fields[2:] if gene in gene_map})
            if len(indices) >= 3:
                result[fields[0]] = np.asarray(indices, dtype=int)
    return result


def _effect_tables(
    real_tpm: np.ndarray,
    synthetic_tpm: np.ndarray,
    samples: pd.DataFrame,
    genes: list[str],
    pathway_file: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression_real = np.log1p(np.maximum(real_tpm, 0.0))
    expression_synthetic = np.log1p(np.maximum(synthetic_tpm, 0.0))
    gene_rows: list[dict[str, object]] = []
    pathway_rows: list[dict[str, object]] = []
    pathways = _read_pathways(pathway_file, genes)
    for tissue in sorted(samples["tissue"].astype(str).unique()):
        tissue_mask = samples["tissue"].astype(str).eq(tissue).to_numpy()
        flight = tissue_mask & samples["condition"].astype(str).eq("flight").to_numpy()
        ground = tissue_mask & samples["condition"].astype(str).eq(
            "ground_control"
        ).to_numpy()
        if flight.sum() < 2 or ground.sum() < 2:
            continue
        real_delta = expression_real[flight].mean(axis=0) - expression_real[
            ground
        ].mean(axis=0)
        synthetic_delta = expression_synthetic[flight].mean(axis=0) - expression_synthetic[
            ground
        ].mean(axis=0)
        for gene, first, second in zip(genes, real_delta, synthetic_delta):
            gene_rows.append(
                {
                    "tissue": tissue,
                    "gene_id": gene,
                    "real_log1p_tpm_delta": float(first),
                    "synthetic_log1p_tpm_delta": float(second),
                    "direction_agrees": bool(np.sign(first) == np.sign(second)),
                }
            )
        for pathway, indices in pathways.items():
            first = float(real_delta[indices].mean())
            second = float(synthetic_delta[indices].mean())
            pathway_rows.append(
                {
                    "tissue": tissue,
                    "pathway": pathway,
                    "landmark_genes": int(len(indices)),
                    "real_mean_gene_delta": first,
                    "synthetic_mean_gene_delta": second,
                    "absolute_error": abs(first - second),
                    "direction_agrees": bool(np.sign(first) == np.sign(second)),
                }
            )
    return pd.DataFrame(gene_rows), pd.DataFrame(pathway_rows)


def _class_probe(
    real_train: np.ndarray,
    train_labels: Iterable[str],
    synthetic: np.ndarray,
    expected_labels: Iterable[str],
    *,
    seed: int,
) -> dict[str, float]:
    return _classifier_predictions(
        real_train,
        np.asarray(list(train_labels), dtype=str),
        synthetic,
        np.asarray(list(expected_labels), dtype=str),
        seed=seed,
    )[0]


def evaluate_conditional(
    config_path: str | Path, *, unlock_test: bool = False
) -> Path:
    config = load_conditional_config(config_path)
    options = config["evaluation"]
    split = str(options["split"])
    if split == "test" and not (unlock_test and bool(options.get("unlock_test", False))):
        raise ValueError("Locked OSDR test evaluation requires both config and CLI unlock")
    if split not in {"validation", "test"}:
        raise ValueError("Conditional evaluation split must be validation or test")
    verify_source(config["run"]["source_root"])
    prepared_path = Path(config["data"]["prepared_h5"])
    prepared = load_conditional_prepared(prepared_path)
    samples = pd.read_csv(prepared_path.with_suffix(".samples.tsv.gz"), sep="\t")
    train_samples = samples.loc[samples["role"].astype(str).eq("train")].reset_index(
        drop=True
    )
    evaluation_samples = samples.loc[
        samples["role"].astype(str).eq(split)
    ].reset_index(drop=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, payload, model_path = _load_model(config, prepared, device)
    betas = quadratic_beta_schedule(
        beta_start=float(config["model"]["beta_start"]),
        beta_end=float(config["model"]["beta_end"]),
        timesteps=int(config["model"]["diffusion_timesteps"]),
    ).to(device)
    output = Path(config["run"]["output_dir"]) / "evaluation" / split
    output.mkdir(parents=True, exist_ok=True)
    metric_indices = _balanced_indices(
        prepared[split]["class_index"],
        int(options["metric_samples"]),
        int(config["run"]["seed"]) + 31,
    )
    real = prepared[split]["expression"][metric_indices]
    sample_frame = evaluation_samples.iloc[metric_indices].reset_index(drop=True)
    class_indices = prepared[split]["class_index"][metric_indices]
    started = time.time()
    synthetic = _sample(
        model=model,
        class_indices=class_indices,
        genes=len(prepared["genes"]),
        classes=len(prepared["classes"]),
        betas=betas,
        timesteps=int(config["model"]["diffusion_timesteps"]),
        batch_size=int(options["quality_batch_size"]),
        seed=int(config["run"]["seed"]) + 101,
        device=device,
    )
    sampling_seconds = time.time() - started
    synthetic_tpm = synthetic * prepared["maxabs_scale"].reshape(1, -1)
    real_tpm = prepared[split]["tpm"][metric_indices]
    fidelity = generated_quality(real, synthetic, max_pr_samples=len(real))
    memorization = memorization_metrics(
        prepared["train"]["expression"],
        synthetic,
        max_samples=max(len(real), 50),
        seed=int(config["run"]["seed"]),
    )
    selection = fidelity_selection(fidelity, memorization)
    effect = _condition_effect(
        real,
        synthetic,
        sample_frame["condition"].astype(str).to_numpy(),
    )

    train_indices = _balanced_indices(
        prepared["train"]["class_index"],
        int(options["quality_samples"]),
        int(config["run"]["seed"]) + 53,
    )
    synthetic_train = _sample(
        model=model,
        class_indices=prepared["train"]["class_index"][train_indices],
        genes=len(prepared["genes"]),
        classes=len(prepared["classes"]),
        betas=betas,
        timesteps=int(config["model"]["diffusion_timesteps"]),
        batch_size=int(options["quality_batch_size"]),
        seed=int(config["run"]["seed"]) + 211,
        device=device,
    )
    train_condition = train_samples.iloc[train_indices]["condition"].astype(str).to_numpy()
    evaluation_condition = sample_frame["condition"].astype(str).to_numpy()
    utility = classifier_utility(
        prepared["train"]["expression"],
        train_samples["condition"].astype(str).to_numpy(),
        real,
        evaluation_condition,
        synthetic_train=synthetic_train,
        synthetic_labels=train_condition,
        allow_augmentation=bool(selection["eligible_for_model_selection"]),
    )
    condition_consistency = _class_probe(
        prepared["train"]["expression"],
        train_samples["condition"].astype(str),
        synthetic,
        sample_frame["condition"].astype(str),
        seed=int(config["run"]["seed"]),
    )
    tissue_consistency = _class_probe(
        prepared["train"]["expression"],
        train_samples["tissue"].astype(str),
        synthetic,
        sample_frame["tissue"].astype(str),
        seed=int(config["run"]["seed"]),
    )
    pca_path = _plot_real_synthetic_pca(real, synthetic, sample_frame, output)
    gene_effects, pathway_effects = _effect_tables(
        real_tpm,
        synthetic_tpm,
        sample_frame,
        prepared["genes"],
        "data/pathways/reactome_current_mouse_ensembl.gmt",
    )
    gene_effect_path = output / "flt_gc_gene_effect_recovery.tsv.gz"
    pathway_effect_path = output / "flt_gc_pathway_effect_recovery.tsv.gz"
    gene_effects.to_csv(
        gene_effect_path, sep="\t", index=False, compression="gzip"
    )
    pathway_effects.to_csv(
        pathway_effect_path, sep="\t", index=False, compression="gzip"
    )
    pathway_summary = {
        "rows": int(len(pathway_effects)),
        "delta_correlation": _correlation(
            pathway_effects.get("real_mean_gene_delta", []),
            pathway_effects.get("synthetic_mean_gene_delta", []),
        ),
        "direction_agreement": (
            float(pathway_effects["direction_agrees"].mean())
            if len(pathway_effects)
            else float("nan")
        ),
    }
    np.savez_compressed(
        output / "conditional_synthetic_expression.npz",
        scaled_expression=synthetic,
        tpm_unclipped=synthetic_tpm.astype(np.float32),
        class_index=class_indices,
        source_row=prepared[split]["source_row"][metric_indices],
        genes=np.asarray(prepared["genes"]),
    )
    summary = {
        "status": "complete",
        "split": split,
        "locked_test_unlocked": bool(split == "test" and unlock_test),
        "model": str(model_path),
        "checkpoint_epoch": int(payload.get("epoch", 0)),
        "checkpoint_global_step": int(payload.get("global_step", 0)),
        "device": torch.cuda.get_device_name(device) if device.type == "cuda" else "cpu",
        "samples": int(len(real)),
        "sampling_seconds": float(sampling_seconds),
        "fidelity_transformed": fidelity,
        "memorization": memorization,
        "model_selection": selection,
        "flt_gc_effect_recovery": effect,
        "flt_gc_classifier_utility": utility,
        "condition_consistency": condition_consistency,
        "tissue_consistency": tissue_consistency,
        "pathway_effect_recovery": pathway_summary,
        "synthetic_tpm_unclipped_min": float(synthetic_tpm.min()),
        "synthetic_tpm_negative_fraction": float((synthetic_tpm < 0).mean()),
        "plots": {"real_vs_synthetic_pca": str(pca_path)},
        "tables": {
            "gene_effect_recovery": str(gene_effect_path),
            "pathway_effect_recovery": str(pathway_effect_path),
        },
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    return summary_path
