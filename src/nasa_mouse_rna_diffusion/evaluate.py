"""Generate paper-style trajectory plots and held-out tissue metrics."""

from __future__ import annotations

import json
from pathlib import Path
import time
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import sqrtm
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, silhouette_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import torch

from nasa_mouse_generative.metrics import fidelity_selection, memorization_metrics
from nasa_mouse_generative.paper_metrics import paper_distribution_metrics

from .config import load_config
from .data import load_prepared
from .upstream import (
    EMA,
    ddim_trajectory,
    model_config,
    quadratic_beta_schedule,
    upstream_model_class,
    verify_source,
)


def _correlation(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=np.float64).ravel()
    second = np.asarray(second, dtype=np.float64).ravel()
    valid = np.isfinite(first) & np.isfinite(second)
    if valid.sum() < 2 or np.std(first[valid]) == 0 or np.std(second[valid]) == 0:
        return float("nan")
    return float(np.corrcoef(first[valid], second[valid])[0, 1])


def _classifier_predictions(
    train_expression: np.ndarray,
    train_labels: np.ndarray,
    test_expression: np.ndarray,
    test_labels: np.ndarray,
    *,
    seed: int,
) -> tuple[dict[str, float], np.ndarray]:
    classifier = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            max_iter=2000,
            solver="lbfgs",
            random_state=int(seed),
        ),
    )
    classifier.fit(train_expression, train_labels)
    prediction = classifier.predict(test_expression)
    return (
        {
            "balanced_accuracy": float(
                balanced_accuracy_score(test_labels, prediction)
            ),
            "accuracy": float(np.mean(prediction == test_labels)),
        },
        prediction,
    )


def _plot_tissue_probe_recall(table: pd.DataFrame, output: Path) -> Path:
    positions = np.arange(len(table))
    figure, axis = plt.subplots(figsize=(9.4, 8.0))
    width = 0.38
    axis.barh(
        positions - width / 2,
        table["real_train_probe_recall"],
        height=width,
        color="#176B87",
        label="real train -> real test",
    )
    axis.barh(
        positions + width / 2,
        table["synthetic_train_probe_recall"],
        height=width,
        color="#C14924",
        label="synthetic train -> real test",
    )
    axis.set_yticks(positions, table["tissue"])
    axis.set_xlim(0, 1.0)
    axis.set_xlabel("Held-out test recall")
    axis.set_title("Per-tissue reverse validation", fontweight="bold", pad=50)
    axis.grid(axis="x", alpha=0.2, linewidth=0.7)
    axis.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
    )
    figure.tight_layout()
    path = output / "per_tissue_reverse_validation.png"
    figure.savefig(path, dpi=220, bbox_inches="tight")
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return path


def _logistic_adversarial_accuracy(
    real: np.ndarray, synthetic: np.ndarray, seed: int
) -> float:
    count = min(len(real), len(synthetic))
    expression = np.concatenate([real[:count], synthetic[:count]])
    labels = np.concatenate(
        [np.zeros(count, dtype=int), np.ones(count, dtype=int)]
    )
    train_x, test_x, train_y, test_y = train_test_split(
        expression,
        labels,
        test_size=0.3,
        random_state=int(seed),
        stratify=labels,
    )
    classifier = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=1000, random_state=int(seed)),
    )
    classifier.fit(train_x, train_y)
    return float(np.mean(classifier.predict(test_x) == test_y))


def _nearest_neighbor_adversarial_accuracy(
    real: np.ndarray, synthetic: np.ndarray
) -> float:
    """Compute the paper's 1-NN real-versus-synthetic adversarial accuracy."""

    real_distances = cdist(real, real)
    synthetic_distances = cdist(synthetic, synthetic)
    np.fill_diagonal(real_distances, np.inf)
    np.fill_diagonal(synthetic_distances, np.inf)
    cross = cdist(real, synthetic)
    real_to_real = real_distances.min(axis=1)
    real_to_synthetic = cross.min(axis=1)
    synthetic_to_real = cross.min(axis=0)
    synthetic_to_synthetic = synthetic_distances.min(axis=1)
    return float(
        0.5
        * (
            np.mean(real_to_synthetic > real_to_real)
            + np.mean(synthetic_to_real > synthetic_to_synthetic)
        )
    )


def _precision_recall(real: np.ndarray, synthetic: np.ndarray, neighbors: int = 10) -> dict[str, float]:
    real_distances = cdist(real, real)
    synthetic_distances = cdist(synthetic, synthetic)
    cross = cdist(real, synthetic)
    real_radius = np.partition(real_distances, neighbors, axis=1)[:, neighbors]
    synthetic_radius = np.partition(
        synthetic_distances, neighbors, axis=1
    )[:, neighbors]
    precision = np.mean((cross <= real_radius[:, None]).any(axis=0))
    recall = np.mean((cross <= synthetic_radius[None, :]).any(axis=1))
    return {"precision": float(precision), "recall": float(recall)}


def _frechet(first: np.ndarray, second: np.ndarray) -> float:
    first_mean = first.mean(axis=0)
    second_mean = second.mean(axis=0)
    first_covariance = np.cov(first, rowvar=False)
    second_covariance = np.cov(second, rowvar=False)
    covariance_mean = sqrtm(first_covariance @ second_covariance)
    if np.iscomplexobj(covariance_mean):
        covariance_mean = covariance_mean.real
    difference = first_mean - second_mean
    return float(
        difference @ difference
        + np.trace(first_covariance + second_covariance - 2 * covariance_mean)
    )


def _load_model(config: dict[str, Any], prepared: dict[str, Any], device: torch.device):
    output = Path(config["run"]["output_dir"])
    model_path = output / "model.pt"
    checkpoint_path = output / "checkpoints/latest.pt"
    path = model_path if model_path.exists() else checkpoint_path
    if not path.exists():
        raise FileNotFoundError(f"No trained model or checkpoint under {output}")
    payload = torch.load(path, map_location=device, weights_only=False)
    namespace = model_config(
        expression_dim=len(prepared["genes"]),
        num_classes=len(prepared["classes"]),
        model=config["model"],
    )
    model = upstream_model_class()(namespace).to(device)
    model.load_state_dict(payload["model_state_dict"])
    ema = EMA(model, float(config["model"]["ema_decay"]))
    ema.load_state_dict(payload["ema_state_dict"], device)
    ema.copy_to(model)
    model.eval()
    return model, payload, path


def _plot_trajectory(
    *,
    snapshots: dict[int, np.ndarray],
    labels: np.ndarray,
    classes: list[str],
    real_background: np.ndarray,
    output: Path,
    seed: int,
) -> tuple[Path, dict[int, np.ndarray], PCA]:
    pca = PCA(n_components=2, random_state=int(seed)).fit(real_background)
    real_coordinates = pca.transform(real_background)
    coordinates = {
        timestep: pca.transform(expression)
        for timestep, expression in snapshots.items()
    }
    colors = plt.get_cmap("tab20", len(classes))
    order = sorted(snapshots, reverse=True)
    figure, axes = plt.subplots(1, len(order), figsize=(16.2, 5.2), sharex=True, sharey=True)
    for axis, timestep in zip(axes, order):
        axis.scatter(
            real_coordinates[:, 0],
            real_coordinates[:, 1],
            s=5,
            alpha=0.10,
            color="#686868",
            edgecolors="none",
            rasterized=True,
            label="real ARCHS4",
        )
        for class_index, class_name in enumerate(classes):
            mask = labels == class_index
            axis.scatter(
                coordinates[timestep][mask, 0],
                coordinates[timestep][mask, 1],
                s=11,
                alpha=0.82,
                color=colors(class_index),
                edgecolors="none",
                rasterized=True,
                label=class_name,
            )
        axis.set_title(f"t = {timestep}", fontsize=14, fontweight="bold")
        axis.set_xlabel("PC1", fontweight="bold")
    axes[0].set_ylabel("PC2", fontweight="bold")
    handles, legend_labels = axes[-1].get_legend_handles_labels()
    figure.legend(
        handles,
        legend_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.03),
        ncol=7,
        fontsize=8,
        frameon=False,
    )
    figure.suptitle(
        "Upstream-parity DDIM generation across ARCHS4 mouse tissues",
        fontsize=15,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0.15, 1, 0.94))
    path = output / "archs4_mouse_ddim_trajectory_pca.png"
    figure.savefig(path, dpi=240, bbox_inches="tight")
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return path, coordinates, pca


def _plot_training_history(run_dir: Path, output: Path) -> Path:
    history_path = run_dir / "training_history.tsv"
    history = pd.read_csv(history_path, sep="\t")
    if history.empty:
        raise ValueError(f"Training history is empty: {history_path}")
    figure, axes = plt.subplots(2, 1, figsize=(9.2, 6.8), sharex=True)
    axes[0].plot(
        history["epoch"], history["loss"], color="#176B87", linewidth=1.2
    )
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Summed noise MSE")
    axes[0].set_title("Paper-parity DDIM training", fontweight="bold")
    error_axis = axes[0].twinx()
    error_axis.plot(
        history["epoch"],
        history["noise_absolute_error"],
        color="#C14924",
        linewidth=1.0,
        alpha=0.8,
    )
    error_axis.set_ylabel("Noise MAE")
    axes[1].plot(
        history["epoch"],
        history["learning_rate"],
        color="#4F6D3A",
        linewidth=1.2,
    )
    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Learning rate")
    for axis in (*axes, error_axis):
        axis.grid(alpha=0.18, linewidth=0.7)
    figure.tight_layout()
    path = output / "training_history.png"
    figure.savefig(path, dpi=220, bbox_inches="tight")
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return path


def evaluate(config_path: str | Path) -> Path:
    config = load_config(config_path)
    verify_source(config["run"]["source_root"])
    prepared = load_prepared(config["data"]["prepared_h5"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, payload, model_path = _load_model(config, prepared, device)
    options = config["evaluation"]
    run_dir = Path(config["run"]["output_dir"])
    output = run_dir / "evaluation"
    output.mkdir(parents=True, exist_ok=True)
    history_figure = _plot_training_history(run_dir, output)
    seed = int(config["run"]["seed"]) + 101
    generator = torch.Generator(device=device).manual_seed(seed)
    sample_count = int(options["trajectory_samples"])
    class_indices = torch.randint(
        0,
        len(prepared["classes"]),
        (sample_count,),
        generator=generator,
        device=device,
    )
    initial_noise = torch.randn(
        (sample_count, len(prepared["genes"])),
        generator=generator,
        device=device,
    )
    betas = quadratic_beta_schedule(
        beta_start=float(config["model"]["beta_start"]),
        beta_end=float(config["model"]["beta_end"]),
        timesteps=int(config["model"]["diffusion_timesteps"]),
    ).to(device)
    requested = tuple(map(int, options["trajectory_timesteps"]))
    collected: dict[int, list[np.ndarray]] = {value: [] for value in requested}
    labels_numpy: list[np.ndarray] = []
    started = time.time()
    batch_size = int(options["trajectory_batch_size"])
    for start in range(0, sample_count, batch_size):
        end = min(start + batch_size, sample_count)
        labels = torch.nn.functional.one_hot(
            class_indices[start:end], num_classes=len(prepared["classes"])
        ).long()
        trajectory = ddim_trajectory(
            initial_noise[start:end],
            labels,
            model,
            betas,
            sequence=range(int(config["model"]["diffusion_timesteps"])),
            snapshot_timesteps=requested,
            eta=0.0,
        )
        for timestep, expression in trajectory.items():
            collected[timestep].append(expression.numpy().astype(np.float32))
        labels_numpy.append(class_indices[start:end].cpu().numpy())
        print(f"[rna-diffusion:evaluate] sampled {end}/{sample_count}", flush=True)
    trajectory_seconds = time.time() - started
    snapshots = {
        timestep: np.concatenate(values) for timestep, values in collected.items()
    }
    generated_labels = np.concatenate(labels_numpy)
    np.savez_compressed(
        output / "trajectory_scaled_expression.npz",
        labels=generated_labels,
        classes=np.asarray(prepared["classes"]),
        **{f"t{timestep}": values for timestep, values in snapshots.items()},
    )

    rng = np.random.default_rng(seed)
    train_expression = prepared["train"]["expression"]
    background_count = min(int(options["pca_background_samples"]), len(train_expression))
    background_indices = rng.choice(
        len(train_expression), background_count, replace=False
    )
    figure_path, coordinates, pca = _plot_trajectory(
        snapshots=snapshots,
        labels=generated_labels,
        classes=prepared["classes"],
        real_background=train_expression[background_indices],
        output=output,
        seed=seed,
    )
    for timestep, values in coordinates.items():
        pd.DataFrame(
            {
                "PC1": values[:, 0],
                "PC2": values[:, 1],
                "class_index": generated_labels,
                "tissue": [prepared["classes"][value] for value in generated_labels],
            }
        ).to_csv(
            output / f"trajectory_pca_t{timestep}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    quality_count = min(int(options["quality_samples"]), len(train_expression))
    quality_labels = prepared["train"]["class_index"][:quality_count]
    quality_generated: list[np.ndarray] = []
    quality_started = time.time()
    quality_batch_size = int(options["quality_batch_size"])
    for start in range(0, quality_count, quality_batch_size):
        end = min(start + quality_batch_size, quality_count)
        class_batch = torch.as_tensor(
            quality_labels[start:end], dtype=torch.long, device=device
        )
        labels = torch.nn.functional.one_hot(
            class_batch, num_classes=len(prepared["classes"])
        ).long()
        noise = torch.randn(
            (end - start, len(prepared["genes"])),
            generator=generator,
            device=device,
        )
        generated = ddim_trajectory(
            noise,
            labels,
            model,
            betas,
            sequence=range(int(config["model"]["diffusion_timesteps"])),
            snapshot_timesteps=(0,),
            eta=0.0,
        )[0]
        quality_generated.append(generated.numpy().astype(np.float32))
        print(
            f"[rna-diffusion:evaluate] quality samples {end}/{quality_count}",
            flush=True,
        )
    quality_final = np.concatenate(quality_generated)
    quality_seconds = time.time() - quality_started
    quality_tpm = quality_final * prepared["maxabs_scale"].reshape(1, -1)
    np.savez_compressed(
        output / "quality_synthetic_expression.npz",
        scaled_expression=quality_final,
        tpm_unclipped=quality_tpm.astype(np.float32),
        class_index=quality_labels,
        classes=np.asarray(prepared["classes"]),
        genes=np.asarray(prepared["genes"]),
        maxabs_scale=prepared["maxabs_scale"],
    )
    pd.DataFrame(
        {
            "synthetic_profile": np.arange(quality_count, dtype=int),
            "class_index": quality_labels,
            "tissue": [prepared["classes"][value] for value in quality_labels],
        }
    ).to_csv(
        output / "quality_synthetic_samples.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    test_expression = prepared["test"]["expression"]
    test_labels = prepared["test"]["class_index"]
    real_baseline, real_prediction = _classifier_predictions(
        train_expression,
        prepared["train"]["class_index"],
        test_expression,
        test_labels,
        seed=seed,
    )
    synthetic_utility, synthetic_prediction = _classifier_predictions(
        quality_final,
        quality_labels,
        test_expression,
        test_labels,
        seed=seed,
    )
    class_values = np.arange(len(prepared["classes"]), dtype=int)
    real_confusion = confusion_matrix(
        test_labels, real_prediction, labels=class_values, normalize="true"
    )
    synthetic_confusion = confusion_matrix(
        test_labels, synthetic_prediction, labels=class_values, normalize="true"
    )
    per_tissue = pd.DataFrame(
        {
            "class_index": class_values,
            "tissue": prepared["classes"],
            "test_profiles": [int(np.sum(test_labels == value)) for value in class_values],
            "real_train_probe_recall": np.diag(real_confusion),
            "synthetic_train_probe_recall": np.diag(synthetic_confusion),
        }
    )
    per_tissue["synthetic_minus_real_recall"] = (
        per_tissue["synthetic_train_probe_recall"]
        - per_tissue["real_train_probe_recall"]
    )
    per_tissue_path = output / "per_tissue_reverse_validation.tsv"
    per_tissue.to_csv(per_tissue_path, sep="\t", index=False)
    confusion_path = output / "synthetic_to_real_test_confusion.tsv"
    pd.DataFrame(
        synthetic_confusion,
        index=prepared["classes"],
        columns=prepared["classes"],
    ).rename_axis("true_tissue").to_csv(confusion_path, sep="\t")
    tissue_probe_figure = _plot_tissue_probe_recall(per_tissue, output)
    metric_count = min(
        int(options["metric_samples"]), len(quality_final), len(train_expression)
    )
    real_metric_indices = rng.choice(
        len(train_expression), metric_count, replace=False
    )
    synthetic_metric_indices = rng.choice(
        len(quality_final), metric_count, replace=False
    )
    real_metric = train_expression[real_metric_indices]
    synthetic_metric = quality_final[synthetic_metric_indices]
    embedding = PCA(n_components=50, random_state=seed).fit(train_expression)
    real_embedding = embedding.transform(real_metric)
    synthetic_embedding = embedding.transform(synthetic_metric)
    precision_recall_pca = _precision_recall(
        real_embedding, synthetic_embedding, neighbors=10
    )
    paper_quality = paper_distribution_metrics(
        real_metric,
        synthetic_metric,
        max_samples=metric_count,
        neighbors=10,
        seed=seed,
    )
    precision_recall_scaled_expression = {
        "precision": paper_quality["precision"],
        "recall": paper_quality["recall"],
    }
    nearest_adversarial = float(paper_quality["adversarial_accuracy"])
    memorization = memorization_metrics(
        train_expression,
        quality_final,
        max_samples=metric_count,
        seed=seed,
    )
    model_selection = fidelity_selection(
        {
            **paper_quality,
            "gene_mean_correlation": _correlation(
                real_metric.mean(axis=0), synthetic_metric.mean(axis=0)
            ),
            "gene_std_correlation": _correlation(
                real_metric.std(axis=0), synthetic_metric.std(axis=0)
            ),
            "real_global_std": float(real_metric.std()),
            "fake_global_std": float(synthetic_metric.std()),
        },
        memorization,
    )
    summary = {
        "status": "complete",
        "model": str(model_path),
        "checkpoint_epoch": int(payload.get("epoch", 0)),
        "checkpoint_global_step": int(payload.get("global_step", 0)),
        "device": torch.cuda.get_device_name(device) if device.type == "cuda" else "cpu",
        "generated_profiles": sample_count,
        "quality_generated_profiles": quality_count,
        "synthetic_scaled_negative_fraction": float((quality_final < 0).mean()),
        "synthetic_scaled_above_one_fraction": float((quality_final > 1).mean()),
        "synthetic_tpm_unclipped_min": float(quality_tpm.min()),
        "synthetic_tpm_unclipped_max": float(quality_tpm.max()),
        "diffusion_steps": int(config["model"]["diffusion_timesteps"]),
        "trajectory_sampling_seconds": float(trajectory_seconds),
        "quality_sampling_seconds": float(quality_seconds),
        "figure": str(figure_path),
        "training_history_figure": str(history_figure),
        "per_tissue_reverse_validation": str(per_tissue_path),
        "per_tissue_reverse_validation_figure": str(tissue_probe_figure),
        "synthetic_to_real_test_confusion": str(confusion_path),
        "pca_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "real_train_to_test_tissue_classifier": real_baseline,
        "synthetic_to_real_test_tissue_classifier": synthetic_utility,
        "synthetic_t0_pca_silhouette": float(
            silhouette_score(coordinates[0], generated_labels)
        ),
        "gene_mean_correlation": _correlation(
            real_metric.mean(axis=0), synthetic_metric.mean(axis=0)
        ),
        "gene_standard_deviation_correlation": _correlation(
            real_metric.std(axis=0), synthetic_metric.std(axis=0)
        ),
        "gene_correlation_matrix_agreement": paper_quality[
            "correlation_matrix_agreement"
        ],
        "precision_recall_in_scaled_l974": precision_recall_scaled_expression,
        "precision_recall_in_train_pca50": precision_recall_pca,
        "memorization": memorization,
        "model_selection": model_selection,
        "frechet_distance_in_train_pca50": paper_quality["frechet_pca"],
        "frechet_real_split_p95_in_train_pca50": paper_quality[
            "frechet_real_split_p95"
        ],
        "frechet_ratio_to_real_split_p95": paper_quality[
            "frechet_ratio_to_real_split_p95"
        ],
        "nearest_neighbor_adversarial_accuracy_in_scaled_l974": (
            nearest_adversarial
        ),
        "logistic_adversarial_accuracy_in_scaled_l974": _logistic_adversarial_accuracy(
            real_metric, synthetic_metric, seed
        ),
        "metric_limitations": (
            "Tissue utility uses a fixed logistic probe and Frechet distance uses a "
            "train-fitted PCA-50 embedding because the paper's human GTEx MLP feature "
            "extractor is not transferable to mouse genes. Precision/recall and the "
            "nearest-neighbor adversarial accuracy are also reported directly in the "
            "scaled 974-gene space, following the paper's unsupervised definitions, "
            "but on a fixed 2,000-profile subset for tractable pairwise distances."
        ),
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "README.md").write_text(
        "# Paper-Style Mouse DDIM Evaluation\n\n"
        "The trajectory uses the EMA weights and all 1,000 upstream DDIM steps. "
        "One PCA basis is fitted only on real ARCHS4 training profiles and reused "
        "for `t=1000`, `t=200`, and `t=0`.\n\n"
        "`training_history.png` shows convergence and the OneCycle schedule. See "
        "`summary.json` for quantitative held-out tissue metrics. Synthetic L974 "
        "model-scale and inverse-MaxAbs TPM matrices are stored together in "
        "`quality_synthetic_expression.npz`; TPM values are not silently clipped. "
        "Per-tissue reverse-validation recall is provided as TSV and PNG/PDF.\n\n"
        f"Synthetic-to-real held-out tissue balanced accuracy is "
        f"{synthetic_utility['balanced_accuracy']:.4f}, compared with "
        f"{real_baseline['balanced_accuracy']:.4f} for real-to-real. Direct "
        f"scaled-L974 precision/recall are "
        f"{precision_recall_scaled_expression['precision']:.4f}/"
        f"{precision_recall_scaled_expression['recall']:.4f}, gene-correlation-"
        f"matrix agreement is {summary['gene_correlation_matrix_agreement']:.4f}, "
        f"and nearest-neighbor adversarial accuracy is "
        f"{summary['nearest_neighbor_adversarial_accuracy_in_scaled_l974']:.3f}. "
        f"The generated matrix has "
        f"{100 * summary['synthetic_scaled_negative_fraction']:.2f}% negative "
        "entries, retained so any downstream nonnegativity policy remains "
        "explicit.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2), flush=True)
    return summary_path
