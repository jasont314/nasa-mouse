"""Held-out fidelity and FLT/GC evaluation for conditional upstream ModelDDIM."""

from __future__ import annotations

import json
from pathlib import Path
import re
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
from nasa_mouse_generative.effect_validation import compare_real_synthetic_effects
from nasa_mouse_generative.metrics import (
    _condition_effect,
    accession_effect_selection,
    classifier_utility,
    conditional_effect_selection,
    fidelity_selection,
    memorization_metrics,
)

from .conditional_config import load_conditional_config
from .conditional_data import load_conditional_prepared
from .evaluate import (
    _classifier_predictions,
    _correlation,
    _load_model,
    _plot_training_history,
)
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
    eta: float,
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
            eta=float(eta),
            generator=generator,
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


def _plot_accession_pca(
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
    accessions = sorted(samples["accession"].astype(str).unique())
    palette = plt.get_cmap("turbo", max(len(accessions), 1))
    figure, axis = plt.subplots(figsize=(8.2, 5.6))
    for index, accession in enumerate(accessions):
        mask = samples["accession"].astype(str).eq(accession).to_numpy()
        axis.scatter(
            real_coordinates[mask, 0],
            real_coordinates[mask, 1],
            s=24,
            color=palette(index),
            marker="o",
            alpha=0.42,
            edgecolors="none",
        )
        axis.scatter(
            synthetic_coordinates[mask, 0],
            synthetic_coordinates[mask, 1],
            s=27,
            color=palette(index),
            marker="x",
            alpha=0.78,
            label=accession,
        )
    axis.set_xlabel("PCA 1")
    axis.set_ylabel("PCA 2")
    axis.set_title("Held-out real (circle) and synthetic (x) by accession")
    axis.grid(alpha=0.16)
    axis.legend(frameon=False, fontsize=7, ncol=2, bbox_to_anchor=(1.02, 1.0))
    figure.tight_layout()
    path = output / "real_vs_synthetic_pca_by_accession.png"
    figure.savefig(path, dpi=220, bbox_inches="tight")
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return path


def _plot_effect_recovery(
    table: pd.DataFrame,
    *,
    real_column: str,
    synthetic_column: str,
    output: Path,
    title: str,
) -> Path | None:
    if table.empty:
        return None
    real = table[real_column].to_numpy(dtype=float)
    synthetic = table[synthetic_column].to_numpy(dtype=float)
    finite = np.isfinite(real) & np.isfinite(synthetic)
    if finite.sum() < 2:
        return None
    figure, axis = plt.subplots(figsize=(6.2, 5.4))
    if "tissue" in table:
        tissues = sorted(table["tissue"].astype(str).unique())
        palette = plt.get_cmap("turbo", max(len(tissues), 1))
        for index, tissue in enumerate(tissues):
            mask = finite & table["tissue"].astype(str).eq(tissue).to_numpy()
            axis.scatter(
                real[mask],
                synthetic[mask],
                s=13,
                alpha=0.55,
                color=palette(index),
                edgecolors="none",
                label=tissue,
            )
        axis.legend(
            frameon=False, fontsize=7, ncol=2, bbox_to_anchor=(1.02, 1.0)
        )
    else:
        axis.scatter(
            real[finite], synthetic[finite], s=13, alpha=0.5, edgecolors="none"
        )
    lower = float(min(real[finite].min(), synthetic[finite].min()))
    upper = float(max(real[finite].max(), synthetic[finite].max()))
    axis.plot([lower, upper], [lower, upper], color="#333333", linestyle="--")
    axis.set_xlabel("Real FLT - GC effect")
    axis.set_ylabel("Synthetic FLT - GC effect")
    axis.set_title(title)
    axis.grid(alpha=0.16)
    figure.tight_layout()
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)
    return output


def _plot_per_tissue_utility(table: pd.DataFrame, output: Path) -> Path | None:
    if table.empty:
        return None
    columns = [
        ("real_train_balanced_accuracy", "Real train"),
        ("synthetic_train_balanced_accuracy", "Synthetic train"),
        ("augmented_train_balanced_accuracy", "Augmented train"),
    ]
    positions = np.arange(len(table))
    width = 0.24
    figure, axis = plt.subplots(figsize=(max(7.0, len(table) * 0.75), 5.2))
    for index, (column, label) in enumerate(columns):
        if column in table and table[column].notna().any():
            axis.bar(
                positions + (index - 1) * width,
                table[column].to_numpy(dtype=float),
                width=width,
                label=label,
            )
    axis.axhline(0.5, color="#444444", linewidth=1, linestyle="--")
    axis.set_xticks(positions)
    axis.set_xticklabels(table["tissue"].astype(str), rotation=45, ha="right")
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("Held-out balanced accuracy")
    axis.set_title("FLT/GC classifier utility by tissue")
    axis.legend(frameon=False)
    axis.grid(axis="y", alpha=0.16)
    figure.tight_layout()
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)
    return output


def _per_tissue_fidelity(
    real: np.ndarray,
    synthetic: np.ndarray,
    samples: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for tissue in sorted(samples["tissue"].astype(str).unique()):
        mask = samples["tissue"].astype(str).eq(tissue).to_numpy()
        if mask.sum() < 5:
            continue
        quality = generated_quality(
            real[mask], synthetic[mask], max_pr_samples=int(mask.sum())
        )
        adversarial_accuracy = float(quality["adversarial_accuracy"])
        rows.append(
            {
                "tissue": tissue,
                "profiles": int(mask.sum()),
                "gene_mean_correlation": quality["gene_mean_correlation"],
                "gene_std_correlation": quality["gene_std_correlation"],
                "precision": quality["precision"],
                "recall": quality["recall"],
                "precision_recall_f1": quality["f1"],
                "nearest_neighbor_adversarial_accuracy": adversarial_accuracy,
                "adversarial_indistinguishability": max(
                    0.0, 1.0 - 2.0 * abs(adversarial_accuracy - 0.5)
                ),
                "synthetic_to_real_global_std_ratio": float(
                    quality["fake_global_std"]
                    / max(float(quality["real_global_std"]), 1e-8)
                ),
            }
        )
    return pd.DataFrame(rows)


def _plot_per_tissue_fidelity(table: pd.DataFrame, output: Path) -> Path | None:
    if table.empty:
        return None
    columns = [
        ("gene_mean_correlation", "Gene mean"),
        ("gene_std_correlation", "Gene SD"),
        ("precision_recall_f1", "PR F1"),
        ("adversarial_indistinguishability", "NN indistinguishability"),
    ]
    positions = np.arange(len(table))
    width = 0.19
    figure, axis = plt.subplots(figsize=(max(8.0, len(table) * 0.9), 5.4))
    for index, (column, label) in enumerate(columns):
        axis.bar(
            positions + (index - 1.5) * width,
            table[column].to_numpy(dtype=float),
            width=width,
            label=label,
        )
    axis.set_xticks(positions)
    axis.set_xticklabels(table["tissue"].astype(str), rotation=45, ha="right")
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("Held-out score")
    axis.set_title("Conditional synthetic fidelity by tissue")
    axis.legend(frameon=False, ncol=2)
    axis.grid(axis="y", alpha=0.16)
    figure.tight_layout()
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)
    return output


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


def _pathway_scores(
    expression: np.ndarray, pathways: dict[str, np.ndarray]
) -> tuple[np.ndarray, list[str]]:
    names = list(pathways)
    if not names:
        return np.empty((len(expression), 0), dtype=np.float32), []
    scores = np.column_stack(
        [expression[:, pathways[name]].mean(axis=1) for name in names]
    )
    return scores.astype(np.float32), names


def _write_accession_validation(
    output: Path,
    *,
    real_tpm: np.ndarray,
    synthetic_tpm: np.ndarray,
    samples: pd.DataFrame,
    genes: list[str],
    pathway_file: str | Path,
) -> tuple[dict[str, object], dict[str, str]]:
    directory = output / "accession_validation"
    directory.mkdir(parents=True, exist_ok=True)
    real_expression = np.log1p(np.maximum(real_tpm, 0.0))
    synthetic_expression = np.log1p(np.maximum(synthetic_tpm, 0.0))
    gene_tables, gene_summary = compare_real_synthetic_effects(
        real_expression, synthetic_expression, samples, genes
    )
    pathways = _read_pathways(pathway_file, genes)
    real_pathways, pathway_names = _pathway_scores(real_expression, pathways)
    synthetic_pathways, _ = _pathway_scores(synthetic_expression, pathways)
    pathway_tables, pathway_summary = compare_real_synthetic_effects(
        real_pathways, synthetic_pathways, samples, pathway_names
    )
    paths: dict[str, str] = {}
    for level, tables in (("gene", gene_tables), ("pathway", pathway_tables)):
        for name, table in tables.items():
            path = directory / f"{level}_{name}.tsv.gz"
            table.to_csv(path, sep="\t", index=False, compression="gzip")
            paths[f"{level}_{name}"] = str(path)
        plot_path = directory / f"{level}_meta_effect_recovery.png"
        plotted = _plot_effect_recovery(
            tables["comparison"],
            real_column="real_meta_effect",
            synthetic_column="synthetic_meta_effect",
            output=plot_path,
            title=f"Accession-aware {level} effect recovery",
        )
        if plotted is not None:
            paths[f"{level}_meta_effect_recovery_plot"] = str(plotted)
    summary = {"gene": gene_summary, "pathway": pathway_summary}
    summary_path = directory / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    paths["summary"] = str(summary_path)
    return summary, paths


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
) -> dict[str, object]:
    train_labels_array = np.asarray(list(train_labels), dtype=str)
    expected_labels_array = np.asarray(list(expected_labels), dtype=str)
    expected_classes = np.unique(expected_labels_array)
    retained = np.isin(train_labels_array, expected_classes)
    retained_classes = np.unique(train_labels_array[retained])
    if len(expected_classes) < 2 or len(retained_classes) < 2:
        return {
            "status": "insufficient_shared_classes",
            "balanced_accuracy": float("nan"),
            "accuracy": float("nan"),
            "train_profiles": int(retained.sum()),
            "excluded_train_profiles": int((~retained).sum()),
            "evaluation_classes": expected_classes.tolist(),
        }
    metrics = _classifier_predictions(
        real_train[retained],
        train_labels_array[retained],
        synthetic,
        expected_labels_array,
        seed=seed,
    )[0]
    return {
        "status": "complete",
        **metrics,
        "train_profiles": int(retained.sum()),
        "excluded_train_profiles": int((~retained).sum()),
        "evaluation_classes": expected_classes.tolist(),
    }


def _per_tissue_classifier_utility(
    *,
    real_train: np.ndarray,
    train_samples: pd.DataFrame,
    synthetic_train: np.ndarray,
    synthetic_samples: pd.DataFrame,
    real_evaluation: np.ndarray,
    evaluation_samples: pd.DataFrame,
    allow_augmentation: bool,
) -> tuple[dict[str, object], pd.DataFrame]:
    results: dict[str, object] = {}
    rows: list[dict[str, object]] = []
    for tissue in sorted(evaluation_samples["tissue"].astype(str).unique()):
        train_mask = train_samples["tissue"].astype(str).eq(tissue).to_numpy()
        synthetic_mask = (
            synthetic_samples["tissue"].astype(str).eq(tissue).to_numpy()
        )
        evaluation_mask = (
            evaluation_samples["tissue"].astype(str).eq(tissue).to_numpy()
        )
        train_labels = train_samples.loc[train_mask, "condition"].astype(str).to_numpy()
        evaluation_labels = (
            evaluation_samples.loc[evaluation_mask, "condition"].astype(str).to_numpy()
        )
        if (
            len(np.unique(train_labels)) < 2
            or len(np.unique(evaluation_labels)) < 2
            or train_mask.sum() < 4
            or evaluation_mask.sum() < 4
        ):
            results[tissue] = {"status": "insufficient_two_condition_data"}
            continue
        synthetic_labels = (
            synthetic_samples.loc[synthetic_mask, "condition"].astype(str).to_numpy()
        )
        utility = classifier_utility(
            real_train[train_mask],
            train_labels,
            real_evaluation[evaluation_mask],
            evaluation_labels,
            synthetic_train=synthetic_train[synthetic_mask],
            synthetic_labels=synthetic_labels,
            allow_augmentation=allow_augmentation,
        )
        results[tissue] = utility
        rows.append(
            {
                "tissue": tissue,
                "real_train_profiles": int(train_mask.sum()),
                "synthetic_train_profiles": int(synthetic_mask.sum()),
                "heldout_real_profiles": int(evaluation_mask.sum()),
                "real_train_balanced_accuracy": utility.get(
                    "real_train_real_evaluation", {}
                ).get("balanced_accuracy", float("nan")),
                "synthetic_train_balanced_accuracy": utility.get(
                    "synthetic_train_real_evaluation", {}
                ).get("balanced_accuracy", float("nan")),
                "augmented_train_balanced_accuracy": utility.get(
                    "real_plus_synthetic_train_real_evaluation", {}
                ).get("balanced_accuracy", float("nan")),
                "augmentation_status": utility.get("augmentation_status", ""),
            }
        )
    return results, pd.DataFrame(rows)


def evaluate_conditional(
    config_path: str | Path,
    *,
    unlock_test: bool = False,
    eta_override: float | None = None,
    evaluation_variant: str = "",
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
    eta = float(options.get("eta", 0.0) if eta_override is None else eta_override)
    if not 0.0 <= eta <= 1.0:
        raise ValueError("DDIM sampling eta must be between zero and one")
    variant = str(evaluation_variant).strip()
    if variant and not re.fullmatch(r"[A-Za-z0-9_.-]+", variant):
        raise ValueError("Evaluation variant may contain only letters, digits, ._- characters")
    output_name = split if not variant else f"{split}_{variant}"
    output = Path(config["run"]["output_dir"]) / "evaluation" / output_name
    output.mkdir(parents=True, exist_ok=True)
    history_plot = _plot_training_history(Path(config["run"]["output_dir"]), output)
    metric_indices = _balanced_indices(
        prepared[split]["class_index"],
        int(options["metric_samples"]),
        int(config["run"]["seed"]) + 31,
    )
    all_class_indices = prepared[split]["class_index"]
    started = time.time()
    synthetic_all = _sample(
        model=model,
        class_indices=all_class_indices,
        genes=len(prepared["genes"]),
        classes=len(prepared["classes"]),
        betas=betas,
        timesteps=int(config["model"]["diffusion_timesteps"]),
        batch_size=int(options["quality_batch_size"]),
        seed=int(config["run"]["seed"]) + 101,
        eta=eta,
        device=device,
    )
    sampling_seconds = time.time() - started
    real = prepared[split]["expression"][metric_indices]
    synthetic = synthetic_all[metric_indices]
    sample_frame = evaluation_samples.iloc[metric_indices].reset_index(drop=True)
    synthetic_tpm_all = synthetic_all * prepared["maxabs_scale"].reshape(1, -1)
    real_tpm_all = prepared[split]["tpm"]
    synthetic_tpm = synthetic_tpm_all[metric_indices]
    real_tpm = real_tpm_all[metric_indices]
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
    condition_effect_gate = conditional_effect_selection(effect)
    accession_validation, accession_paths = _write_accession_validation(
        output,
        real_tpm=real_tpm_all,
        synthetic_tpm=synthetic_tpm_all,
        samples=evaluation_samples.reset_index(drop=True),
        genes=prepared["genes"],
        pathway_file="data/pathways/reactome_current_mouse_ensembl.gmt",
    )
    accession_effect_gate = accession_effect_selection(
        accession_validation["gene"]
    )
    accession_effect_gate["feature_level"] = "gene"
    per_tissue_fidelity_table = _per_tissue_fidelity(
        real, synthetic, sample_frame
    )
    per_tissue_fidelity_path = output / "per_tissue_fidelity.tsv"
    per_tissue_fidelity_table.to_csv(
        per_tissue_fidelity_path, sep="\t", index=False
    )
    per_tissue_fidelity_plot = _plot_per_tissue_fidelity(
        per_tissue_fidelity_table, output / "per_tissue_fidelity.png"
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
        eta=eta,
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
        allow_augmentation=bool(
            selection["eligible_for_model_selection"]
            and condition_effect_gate["passed"]
            and accession_effect_gate["passed"]
        ),
    )
    synthetic_train_samples = train_samples.iloc[train_indices].reset_index(drop=True)
    per_tissue_utility, per_tissue_utility_table = _per_tissue_classifier_utility(
        real_train=prepared["train"]["expression"],
        train_samples=train_samples,
        synthetic_train=synthetic_train,
        synthetic_samples=synthetic_train_samples,
        real_evaluation=prepared[split]["expression"],
        evaluation_samples=evaluation_samples,
        allow_augmentation=bool(
            selection["eligible_for_model_selection"]
            and condition_effect_gate["passed"]
            and accession_effect_gate["passed"]
        ),
    )
    per_tissue_utility_path = output / "per_tissue_flt_gc_classifier_utility.tsv"
    per_tissue_utility_table.to_csv(
        per_tissue_utility_path, sep="\t", index=False
    )
    per_tissue_utility_plot = _plot_per_tissue_utility(
        per_tissue_utility_table,
        output / "per_tissue_flt_gc_classifier_utility.png",
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
    accession_pca_path = _plot_accession_pca(
        real, synthetic, sample_frame, output
    )
    gene_effects, pathway_effects = _effect_tables(
        real_tpm_all,
        synthetic_tpm_all,
        evaluation_samples.reset_index(drop=True),
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
    gene_effect_plot = _plot_effect_recovery(
        gene_effects,
        real_column="real_log1p_tpm_delta",
        synthetic_column="synthetic_log1p_tpm_delta",
        output=output / "flt_gc_gene_effect_recovery.png",
        title="Per-tissue gene effect recovery",
    )
    pathway_effect_plot = _plot_effect_recovery(
        pathway_effects,
        real_column="real_mean_gene_delta",
        synthetic_column="synthetic_mean_gene_delta",
        output=output / "flt_gc_pathway_effect_recovery.png",
        title="Per-tissue Reactome effect recovery",
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
        scaled_expression=synthetic_all,
        tpm_unclipped=synthetic_tpm_all.astype(np.float32),
        tpm_nonnegative=np.maximum(synthetic_tpm_all, 0.0).astype(np.float32),
        class_index=all_class_indices,
        source_row=prepared[split]["source_row"],
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
        "samples": int(len(synthetic_all)),
        "fidelity_metric_samples": int(len(real)),
        "sampling_seconds": float(sampling_seconds),
        "sampling_eta": eta,
        "evaluation_variant": variant,
        "fidelity_transformed": fidelity,
        "memorization": memorization,
        "model_selection": selection,
        "per_tissue_fidelity": {
            str(row["tissue"]): {
                key: value
                for key, value in row.items()
                if key != "tissue"
            }
            for row in per_tissue_fidelity_table.to_dict(orient="records")
        },
        "flt_gc_effect_recovery": effect,
        "conditional_effect_gate": condition_effect_gate,
        "accession_effect_gate": accession_effect_gate,
        "flt_gc_classifier_utility": utility,
        "per_tissue_flt_gc_classifier_utility": per_tissue_utility,
        "condition_consistency": condition_consistency,
        "tissue_consistency": tissue_consistency,
        "pathway_effect_recovery": pathway_summary,
        "accession_effect_validation": accession_validation,
        "synthetic_tpm_unclipped_min": float(synthetic_tpm_all.min()),
        "synthetic_tpm_negative_fraction": float(
            (synthetic_tpm_all < 0).mean()
        ),
        "inverse_transform_policy": {
            "audit_matrix": "tpm_unclipped",
            "downstream_export_matrix": "tpm_nonnegative",
            "negative_values": "clip_to_zero_after_quality_metrics",
        },
        "plots": {
            "real_vs_synthetic_pca": str(pca_path),
            "real_vs_synthetic_pca_by_accession": str(accession_pca_path),
            "training_history": str(history_plot),
            "gene_effect_recovery": (
                str(gene_effect_plot) if gene_effect_plot is not None else ""
            ),
            "pathway_effect_recovery": (
                str(pathway_effect_plot)
                if pathway_effect_plot is not None
                else ""
            ),
            "per_tissue_flt_gc_classifier_utility": (
                str(per_tissue_utility_plot)
                if per_tissue_utility_plot is not None
                else ""
            ),
            "per_tissue_fidelity": (
                str(per_tissue_fidelity_plot)
                if per_tissue_fidelity_plot is not None
                else ""
            ),
        },
        "tables": {
            "gene_effect_recovery": str(gene_effect_path),
            "pathway_effect_recovery": str(pathway_effect_path),
            "per_tissue_flt_gc_classifier_utility": str(
                per_tissue_utility_path
            ),
            "per_tissue_fidelity": str(per_tissue_fidelity_path),
            **accession_paths,
        },
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    return summary_path
