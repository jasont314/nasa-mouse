"""Held-out fidelity, memorization, and FLT/GC utility metrics."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from nasa_mouse_diffusion.evaluate import generated_quality

from .adapters.base import ModelAdapter
from .preprocessing import FittedPreprocessor
from .training_data import DataPartition


def _safe_correlation(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if len(first) < 2 or np.std(first) == 0 or np.std(second) == 0:
        return float("nan")
    return float(np.corrcoef(first, second)[0, 1])


def _condition_effect(
    real: np.ndarray, fake: np.ndarray, labels: np.ndarray
) -> dict[str, float]:
    labels = np.asarray(labels, dtype=str)
    flight = labels == "flight"
    ground = labels == "ground_control"
    if flight.sum() < 2 or ground.sum() < 2:
        return {
            "delta_correlation": float("nan"),
            "delta_rmse": float("nan"),
            "direction_agreement": float("nan"),
        }
    real_delta = real[flight].mean(axis=0) - real[ground].mean(axis=0)
    fake_delta = fake[flight].mean(axis=0) - fake[ground].mean(axis=0)
    return {
        "delta_correlation": _safe_correlation(real_delta, fake_delta),
        "delta_rmse": float(np.sqrt(np.mean((real_delta - fake_delta) ** 2))),
        "direction_agreement": float(
            np.mean(np.sign(real_delta) == np.sign(fake_delta))
        ),
    }


def conditional_effect_selection(effect: dict[str, float]) -> dict[str, object]:
    minimums = {"delta_correlation": 0.30, "direction_agreement": 0.55}
    observed = {
        name: float(effect.get(name, float("nan"))) for name in minimums
    }
    passed = bool(
        all(
            np.isfinite(observed[name]) and observed[name] >= minimum
            for name, minimum in minimums.items()
        )
    )
    return {"passed": passed, "observed": observed, "minimums": minimums}


def _classifier() -> Any:
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced", random_state=0),
    )


def _score_classifier(model, matrix: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    prediction = model.predict(matrix)
    result = {
        "balanced_accuracy": float(balanced_accuracy_score(labels, prediction))
    }
    if len(np.unique(labels)) == 2 and hasattr(model, "predict_proba"):
        probability = model.predict_proba(matrix)
        classes = list(model.classes_)
        flight_index = classes.index("flight") if "flight" in classes else 1
        result["roc_auc"] = float(
            roc_auc_score(labels == "flight", probability[:, flight_index])
        )
    else:
        result["roc_auc"] = float("nan")
    return result


def classifier_utility(
    real_train: np.ndarray,
    train_labels: np.ndarray,
    real_evaluation: np.ndarray,
    evaluation_labels: np.ndarray,
    *,
    synthetic_train: np.ndarray | None = None,
    synthetic_labels: np.ndarray | None = None,
    allow_augmentation: bool = True,
) -> dict[str, object]:
    train_labels = np.asarray(train_labels, dtype=str)
    evaluation_labels = np.asarray(evaluation_labels, dtype=str)
    if len(np.unique(train_labels)) < 2 or len(np.unique(evaluation_labels)) < 2:
        return {"status": "insufficient_two_condition_data"}
    real_model = _classifier().fit(real_train, train_labels)
    result: dict[str, object] = {
        "status": "complete",
        "real_train_real_evaluation": _score_classifier(
            real_model, real_evaluation, evaluation_labels
        ),
    }
    if synthetic_train is not None:
        if synthetic_labels is None:
            if len(synthetic_train) != len(train_labels):
                raise ValueError(
                    "Synthetic labels are required when synthetic and real counts differ"
                )
            synthetic_labels = train_labels
        synthetic_labels = np.asarray(synthetic_labels, dtype=str)
        if len(synthetic_labels) != len(synthetic_train):
            raise ValueError("Synthetic expression and labels have different lengths")
        if len(np.unique(synthetic_labels)) < 2:
            result["synthetic_status"] = "insufficient_two_condition_data"
            return result
        synthetic_model = _classifier().fit(synthetic_train, synthetic_labels)
        result["synthetic_train_real_evaluation"] = _score_classifier(
            synthetic_model, real_evaluation, evaluation_labels
        )
        if allow_augmentation:
            augmented_model = _classifier().fit(
                np.concatenate([real_train, synthetic_train]),
                np.concatenate([train_labels, synthetic_labels]),
            )
            result["real_plus_synthetic_train_real_evaluation"] = _score_classifier(
                augmented_model, real_evaluation, evaluation_labels
            )
            result["augmentation_status"] = "evaluated_after_quality_gates_passed"
        else:
            result["augmentation_status"] = "blocked_by_generator_quality_gates"
    return result


def _synthetic_training_profiles(
    partition: DataPartition,
    *,
    ratio: float,
    samples_per_covariate_profile: int,
    max_samples: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Sample observed joint covariate profiles without inventing combinations."""

    labels = partition.obs["condition"].astype(str).to_numpy()
    columns = [f"category_{index}" for index in range(partition.categories.shape[1])]
    frame = pd.DataFrame(partition.categories, columns=columns)
    frame["condition_label"] = labels
    profiles = frame.drop_duplicates().reset_index(drop=True)
    repetitions = max(1, int(samples_per_covariate_profile))
    pool = profiles.loc[profiles.index.repeat(repetitions)].reset_index(drop=True)
    requested = max(2, int(round(len(partition) * float(ratio))))
    if max_samples > 0:
        requested = min(requested, max(2, int(round(max_samples * float(ratio)))))
    target = min(requested, len(pool))
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    pool_labels = pool["condition_label"].astype(str).to_numpy()
    for label in sorted(set(pool_labels)):
        positions = np.flatnonzero(pool_labels == label)
        share = max(1, int(round(target * len(positions) / len(pool))))
        selected.extend(
            rng.choice(positions, min(share, len(positions)), replace=False).tolist()
        )
    if len(selected) > target:
        selected = rng.choice(selected, target, replace=False).tolist()
    elif len(selected) < target:
        remaining = np.setdiff1d(np.arange(len(pool)), np.asarray(selected, dtype=int))
        extra = rng.choice(
            remaining, min(target - len(selected), len(remaining)), replace=False
        )
        selected.extend(extra.tolist())
    selected_array = np.asarray(selected, dtype=int)
    return (
        pool.iloc[selected_array][columns].to_numpy(dtype=np.int64),
        pool.iloc[selected_array]["condition_label"].astype(str).to_numpy(),
        {
            "ratio": float(ratio),
            "requested": int(requested),
            "generated": int(len(selected_array)),
            "observed_profiles": int(len(profiles)),
            "per_profile_cap": repetitions,
            "capacity_limited": bool(len(selected_array) < requested),
        },
    )


def memorization_metrics(
    train: np.ndarray, synthetic: np.ndarray, *, max_samples: int, seed: int
) -> dict[str, float]:
    n_train = min(len(train), int(max_samples))
    n_fake = min(len(synthetic), int(max_samples))
    if n_train < 5 or n_fake < 2:
        return {
            "nearest_train_distance_median": float("nan"),
            "training_leave_one_out_distance_p01": float("nan"),
            "fraction_below_training_p01": float("nan"),
        }
    rng = np.random.default_rng(seed)
    train = train[rng.choice(len(train), n_train, replace=False)]
    synthetic = synthetic[rng.choice(len(synthetic), n_fake, replace=False)]
    dimensions = min(50, train.shape[1], n_train - 1)
    pca = PCA(n_components=dimensions, random_state=seed)
    train_embedding = pca.fit_transform(train)
    synthetic_embedding = pca.transform(synthetic)
    neighbors = NearestNeighbors(n_neighbors=2).fit(train_embedding)
    training_distances = neighbors.kneighbors(train_embedding)[0][:, 1]
    synthetic_distances = neighbors.kneighbors(
        synthetic_embedding, n_neighbors=1
    )[0][:, 0]
    threshold = float(np.quantile(training_distances, 0.01))
    return {
        "nearest_train_distance_median": float(np.median(synthetic_distances)),
        "training_leave_one_out_distance_p01": threshold,
        "fraction_below_training_p01": float(
            np.mean(synthetic_distances < threshold)
        ),
    }


def fidelity_selection(
    fidelity: dict[str, float], memorization: dict[str, float]
) -> dict[str, object]:
    def bounded(value: object) -> float:
        number = float(value)
        if not np.isfinite(number):
            return 0.0
        return float(np.clip(number, 0.0, 1.0))

    adversarial = float(fidelity.get("adversarial_accuracy", float("nan")))
    adversarial_score = (
        float(np.clip(1.0 - 2.0 * abs(adversarial - 0.5), 0.0, 1.0))
        if np.isfinite(adversarial)
        else 0.0
    )
    components = {
        "gene_mean_correlation": bounded(
            fidelity.get("gene_mean_correlation", float("nan"))
        ),
        "gene_std_correlation": bounded(
            fidelity.get("gene_std_correlation", float("nan"))
        ),
        "precision_recall_f1": bounded(fidelity.get("f1", float("nan"))),
        "adversarial_indistinguishability": adversarial_score,
    }
    real_std = float(fidelity.get("real_global_std", 0.0))
    fake_std = float(fidelity.get("fake_global_std", 0.0))
    diversity_ratio = fake_std / max(real_std, 1e-8)
    composite = float(np.mean(list(components.values())))
    fidelity_thresholds = {
        "heldout_fidelity_composite": 0.70,
        "gene_mean_correlation": 0.80,
        "gene_std_correlation": 0.70,
        "precision_recall_f1": 0.50,
        "adversarial_indistinguishability": 0.50,
    }
    fidelity_observed = {
        "heldout_fidelity_composite": composite,
        **components,
    }
    fidelity_pass = bool(
        all(
            fidelity_observed[name] >= minimum
            for name, minimum in fidelity_thresholds.items()
        )
    )
    diversity_pass = bool(
        float(fidelity.get("recall", 0.0)) >= 0.1 and diversity_ratio >= 0.1
    )
    memorization_fraction = float(
        memorization.get("fraction_below_training_p01", float("nan"))
    )
    memorization_pass = bool(
        np.isfinite(memorization_fraction) and memorization_fraction <= 0.05
    )
    return {
        "heldout_fidelity_composite": composite,
        "components": components,
        "fidelity_gate": {
            "passed": fidelity_pass,
            "observed": fidelity_observed,
            "minimums": fidelity_thresholds,
        },
        "diversity_gate": {
            "passed": diversity_pass,
            "recall_minimum": 0.1,
            "global_std_ratio": diversity_ratio,
            "global_std_ratio_minimum": 0.1,
        },
        "memorization_gate": {
            "passed": memorization_pass,
            "fraction_below_training_p01_maximum": 0.05,
        },
        "eligible_for_model_selection": (
            fidelity_pass and diversity_pass and memorization_pass
        ),
    }


def _plot_representation(
    path: Path,
    embeddings: np.ndarray,
    labels: np.ndarray,
    *,
    title: str = "Held-out representation, condition input neutralized",
) -> str:
    if len(embeddings) < 3 or embeddings.shape[1] < 2:
        return ""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coordinates = PCA(n_components=2, random_state=0).fit_transform(embeddings)
    unique_labels = sorted(set(labels))
    fixed_colors = {"flight": "#c43c39", "ground_control": "#2878b5"}
    palette = plt.get_cmap("turbo", max(len(unique_labels), 1))
    colors = {
        label: fixed_colors.get(label, palette(index))
        for index, label in enumerate(unique_labels)
    }
    figure, axis = plt.subplots(figsize=(6.4, 5.0))
    for label in unique_labels:
        mask = labels == label
        axis.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            s=24,
            alpha=0.75,
            label=label,
            color=colors.get(label, "#666666"),
            edgecolors="none",
        )
    axis.set_xlabel("PCA 1")
    axis.set_ylabel("PCA 2")
    axis.set_title(title)
    axis.legend(
        frameon=False,
        fontsize=7,
        ncol=2 if len(unique_labels) > 8 else 1,
        bbox_to_anchor=(1.02, 1.0) if len(unique_labels) > 8 else None,
        loc="upper left" if len(unique_labels) > 8 else "best",
    )
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return str(path)


def _plot_representation_umap(
    path: Path,
    embeddings: np.ndarray,
    labels: np.ndarray,
    *,
    title: str,
    seed: int,
) -> str:
    if len(embeddings) < 5 or embeddings.shape[1] < 2:
        return ""
    try:
        from umap import UMAP
    except ImportError:
        return ""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coordinates = UMAP(
        n_components=2,
        n_neighbors=min(15, len(embeddings) - 1),
        min_dist=0.2,
        random_state=int(seed),
    ).fit_transform(embeddings)
    unique_labels = sorted(set(labels))
    palette = plt.get_cmap("turbo", max(len(unique_labels), 1))
    figure, axis = plt.subplots(figsize=(6.4, 5.0))
    for index, label in enumerate(unique_labels):
        mask = labels == label
        axis.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            s=22,
            alpha=0.72,
            label=label,
            color=palette(index),
            edgecolors="none",
        )
    axis.set_xlabel("UMAP 1")
    axis.set_ylabel("UMAP 2")
    axis.set_title(title)
    axis.legend(
        frameon=False,
        fontsize=7,
        ncol=2 if len(unique_labels) > 8 else 1,
        bbox_to_anchor=(1.02, 1.0) if len(unique_labels) > 8 else None,
        loc="upper left" if len(unique_labels) > 8 else "best",
    )
    figure.tight_layout()
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return str(path)


def _plot_generation(
    output: Path,
    real: np.ndarray,
    fake: np.ndarray,
    labels: np.ndarray,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths: dict[str, str] = {}
    combined = np.concatenate([real, fake])
    if len(real) >= 3 and combined.shape[1] >= 2:
        coordinates = PCA(n_components=2, random_state=0).fit_transform(combined)
        real_coordinates = coordinates[: len(real)]
        fake_coordinates = coordinates[len(real) :]
        colors = {"flight": "#c43c39", "ground_control": "#2878b5"}
        figure, axis = plt.subplots(figsize=(6.8, 5.2))
        for label in sorted(set(labels)):
            mask = labels == label
            color = colors.get(label, "#666666")
            axis.scatter(
                real_coordinates[mask, 0],
                real_coordinates[mask, 1],
                s=26,
                alpha=0.75,
                color=color,
                marker="o",
                label=f"real {label}",
                edgecolors="none",
            )
            axis.scatter(
                fake_coordinates[mask, 0],
                fake_coordinates[mask, 1],
                s=28,
                alpha=0.75,
                color=color,
                marker="x",
                label=f"synthetic {label}",
            )
        axis.set_xlabel("PCA 1")
        axis.set_ylabel("PCA 2")
        axis.set_title("Held-out real and conditioned synthetic expression")
        axis.legend(frameon=False, fontsize=8, ncol=2)
        figure.tight_layout()
        pca_path = output / "real_vs_synthetic_pca.png"
        figure.savefig(pca_path, dpi=180)
        plt.close(figure)
        paths["real_vs_synthetic_pca"] = str(pca_path)

    figure, axis = plt.subplots(figsize=(5.4, 5.2))
    real_mean = real.mean(axis=0)
    fake_mean = fake.mean(axis=0)
    axis.scatter(real_mean, fake_mean, s=12, alpha=0.55, color="#3b6f8f")
    limits = [
        min(float(real_mean.min()), float(fake_mean.min())),
        max(float(real_mean.max()), float(fake_mean.max())),
    ]
    axis.plot(limits, limits, color="#333333", linewidth=1, linestyle="--")
    axis.set_xlabel("Real gene mean")
    axis.set_ylabel("Synthetic gene mean")
    axis.set_title("Gene-mean fidelity")
    figure.tight_layout()
    mean_path = output / "gene_mean_fidelity.png"
    figure.savefig(mean_path, dpi=180)
    plt.close(figure)
    paths["gene_mean_fidelity"] = str(mean_path)
    return paths


def _subsample_partition(
    partition: DataPartition, max_samples: int, seed: int
) -> DataPartition:
    if max_samples <= 0 or len(partition) <= max_samples:
        return partition
    rng = np.random.default_rng(seed)
    labels = partition.obs["condition"].astype(str).to_numpy()
    selected: list[int] = []
    for label in sorted(set(labels)):
        positions = np.flatnonzero(labels == label)
        take = max(1, round(max_samples * len(positions) / len(labels)))
        selected.extend(rng.choice(positions, min(take, len(positions)), replace=False))
    selected = sorted(selected[:max_samples])
    return DataPartition(
        name=partition.name,
        matrix=partition.matrix[selected],
        obs=partition.obs.iloc[selected].reset_index(drop=True),
        categories=partition.categories[selected],
        weights=partition.weights[selected],
    )


def _neutralize_condition(
    partition: DataPartition, adapter: ModelAdapter
) -> DataPartition:
    categories = partition.categories.copy()
    if "condition" in adapter.covariates:
        # Category code 1 is the encoder's explicit __unknown__ token. Using one
        # constant prevents the supplied FLT/GC label from trivially leaking into
        # a representation-quality classifier.
        categories[:, adapter.covariates.index("condition")] = 1
    return DataPartition(
        name=partition.name,
        matrix=partition.matrix,
        obs=partition.obs,
        categories=categories,
        weights=partition.weights,
    )


def write_embeddings(
    path: Path, partition: DataPartition, embeddings: np.ndarray, prefix: str
) -> Path:
    columns = [f"{prefix}_{index:04d}" for index in range(embeddings.shape[1])]
    frame = pd.concat(
        [
            partition.obs.reset_index(drop=True),
            pd.DataFrame(embeddings, columns=columns),
        ],
        axis=1,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def evaluate_model(
    adapter: ModelAdapter,
    partitions: dict[str, DataPartition],
    preprocessor: FittedPreprocessor,
    *,
    split: str,
    output_dir: str | Path,
    seed: int,
    max_samples: int,
    save_generated_matrix: bool,
    samples_per_covariate_profile: int = 1,
    synthetic_to_real_ratios: tuple[float, ...] = (1.0,),
) -> Path:
    if split not in partitions:
        raise ValueError(f"Unknown evaluation split: {split}")
    evaluation = partitions[split]
    if len(evaluation) == 0:
        raise ValueError(f"Evaluation split {split!r} is empty")
    train = partitions["train"]
    metric_evaluation = _subsample_partition(evaluation, max_samples, seed)
    metric_train = _subsample_partition(train, max_samples, seed + 1)
    output = Path(output_dir) / split
    output.mkdir(parents=True, exist_ok=True)

    train_embeddings = adapter.encode(_neutralize_condition(metric_train, adapter))
    evaluation_embeddings = adapter.encode(
        _neutralize_condition(metric_evaluation, adapter)
    )
    write_embeddings(
        output / "evaluation_embeddings.tsv.gz",
        metric_evaluation,
        evaluation_embeddings,
        adapter.adapter_id.upper(),
    )
    representation_utility = classifier_utility(
        train_embeddings,
        metric_train.obs["condition"].astype(str).to_numpy(),
        evaluation_embeddings,
        metric_evaluation.obs["condition"].astype(str).to_numpy(),
    )
    representation_plot = _plot_representation(
        output / "evaluation_embedding_pca.png",
        evaluation_embeddings,
        metric_evaluation.obs["condition"].astype(str).to_numpy(),
    )
    tissue_train_labels = metric_train.obs["tissue"].astype(str).to_numpy()
    tissue_evaluation_labels = (
        metric_evaluation.obs["tissue"].astype(str).to_numpy()
    )
    representation_tissue_utility = classifier_utility(
        train_embeddings,
        tissue_train_labels,
        evaluation_embeddings,
        tissue_evaluation_labels,
    )
    expression_tissue_utility = classifier_utility(
        metric_train.matrix,
        tissue_train_labels,
        metric_evaluation.matrix,
        tissue_evaluation_labels,
    )
    tissue_pca_plot = _plot_representation(
        output / "evaluation_embedding_pca_by_tissue.png",
        evaluation_embeddings,
        tissue_evaluation_labels,
        title="Held-out representation by tissue",
    )
    tissue_umap_plot = _plot_representation_umap(
        output / "evaluation_embedding_umap_by_tissue.png",
        evaluation_embeddings,
        tissue_evaluation_labels,
        title="Held-out representation by tissue",
        seed=seed,
    )
    harmonization_audit = preprocessor.audit()
    summary: dict[str, object] = {
        "adapter_id": adapter.adapter_id,
        "split": split,
        "samples": len(metric_evaluation),
        "harmonization": harmonization_audit,
        "prediction_metric_interpretation": (
            "outcome_informed_preprocessing_not_valid_for_blind_flt_gc_prediction"
            if harmonization_audit["outcome_informed"]
            else "outcome_blind_preprocessing"
        ),
        "representation_condition_input": "neutralized_to___unknown__",
        "representation_flt_gc_utility": representation_utility,
        "representation_tissue_utility": representation_tissue_utility,
        "expression_tissue_utility": expression_tissue_utility,
        "plots": {
            "evaluation_embedding_pca": representation_plot,
            "evaluation_embedding_pca_by_tissue": tissue_pca_plot,
            "evaluation_embedding_umap_by_tissue": tissue_umap_plot,
        },
        "generation": {"status": "not_supported"},
    }

    if adapter.supports_generation:
        fake_evaluation = adapter.generate(
            metric_evaluation.categories, seed=seed + 17
        )
        fidelity = generated_quality(
            metric_evaluation.matrix,
            fake_evaluation,
            max_pr_samples=max_samples,
        )
        effect = _condition_effect(
            metric_evaluation.matrix,
            fake_evaluation,
            metric_evaluation.obs["condition"].astype(str).to_numpy(),
        )
        memorization = memorization_metrics(
            metric_train.matrix,
            fake_evaluation,
            max_samples=max_samples,
            seed=seed,
        )
        fake_normalized = preprocessor.inverse_transform(
            fake_evaluation, metric_evaluation.obs["study"]
        )
        normalized_real = preprocessor.inverse_transform(
            metric_evaluation.matrix, metric_evaluation.obs["study"]
        )
        normalized_quality = {
            "gene_mean_correlation": _safe_correlation(
                normalized_real.mean(axis=0), fake_normalized.mean(axis=0)
            ),
            "gene_std_correlation": _safe_correlation(
                normalized_real.std(axis=0), fake_normalized.std(axis=0)
            ),
            "output_units": preprocessor.output_units,
            "minimum": float(fake_normalized.min()),
            "maximum": float(fake_normalized.max()),
        }
        selection = fidelity_selection(fidelity, memorization)
        condition_effect_gate = conditional_effect_selection(effect)
        utility_by_ratio: dict[str, object] = {}
        generation_audit: dict[str, object] = {}
        for ratio_index, ratio in enumerate(synthetic_to_real_ratios):
            categories, synthetic_labels, audit = _synthetic_training_profiles(
                metric_train,
                ratio=float(ratio),
                samples_per_covariate_profile=samples_per_covariate_profile,
                max_samples=max_samples,
                seed=seed + 1000 + ratio_index,
            )
            key = f"ratio_{float(ratio):g}"
            generation_audit[key] = audit
            synthetic_train = adapter.generate(
                categories, seed=seed + 2000 + ratio_index
            )
            utility_by_ratio[key] = classifier_utility(
                metric_train.matrix,
                metric_train.obs["condition"].astype(str).to_numpy(),
                metric_evaluation.matrix,
                metric_evaluation.obs["condition"].astype(str).to_numpy(),
                synthetic_train=synthetic_train,
                synthetic_labels=synthetic_labels,
                allow_augmentation=bool(
                    selection["eligible_for_model_selection"]
                    and condition_effect_gate["passed"]
                ),
            )
        generation_plots = _plot_generation(
            output,
            metric_evaluation.matrix,
            fake_evaluation,
            metric_evaluation.obs["condition"].astype(str).to_numpy(),
        )
        summary["plots"].update(generation_plots)
        summary["generation"] = {
            "status": "complete",
            "fidelity_transformed": fidelity,
            "fidelity_normalized_units": normalized_quality,
            "flt_gc_effect_recovery": effect,
            "conditional_effect_gate": condition_effect_gate,
            "memorization": memorization,
            "model_selection": selection,
            "configured_generation": generation_audit,
            "expression_flt_gc_utility_by_ratio": utility_by_ratio,
        }
        if save_generated_matrix:
            np.savez_compressed(
                output / "generated_expression.npz",
                transformed=fake_evaluation,
                normalized=fake_normalized,
            )

    summary_path = output / "evaluation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    flat_rows = []

    def flatten(prefix: str, value: object) -> None:
        if isinstance(value, dict):
            for key, item in value.items():
                flatten(f"{prefix}.{key}" if prefix else key, item)
        elif isinstance(value, (str, int, float, bool)) or value is None:
            flat_rows.append({"metric": prefix, "value": value})

    flatten("", summary)
    pd.DataFrame(flat_rows).to_csv(
        output / "evaluation_metrics.tsv", sep="\t", index=False
    )
    return summary_path
