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
        synthetic_model = _classifier().fit(synthetic_train, train_labels)
        augmented_model = _classifier().fit(
            np.concatenate([real_train, synthetic_train]),
            np.concatenate([train_labels, train_labels]),
        )
        result["synthetic_train_real_evaluation"] = _score_classifier(
            synthetic_model, real_evaluation, evaluation_labels
        )
        result["real_plus_synthetic_train_real_evaluation"] = _score_classifier(
            augmented_model, real_evaluation, evaluation_labels
        )
    return result


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
    summary: dict[str, object] = {
        "adapter_id": adapter.adapter_id,
        "split": split,
        "samples": len(metric_evaluation),
        "representation_condition_input": "neutralized_to___unknown__",
        "representation_flt_gc_utility": representation_utility,
        "generation": {"status": "not_supported"},
    }

    if adapter.supports_generation:
        fake_evaluation = adapter.generate(
            metric_evaluation.categories, seed=seed + 17
        )
        fake_train = adapter.generate(metric_train.categories, seed=seed + 29)
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
        expression_utility = classifier_utility(
            metric_train.matrix,
            metric_train.obs["condition"].astype(str).to_numpy(),
            metric_evaluation.matrix,
            metric_evaluation.obs["condition"].astype(str).to_numpy(),
            synthetic_train=fake_train,
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
        summary["generation"] = {
            "status": "complete",
            "fidelity_transformed": fidelity,
            "fidelity_normalized_units": normalized_quality,
            "flt_gc_effect_recovery": effect,
            "memorization": memorization,
            "expression_flt_gc_utility": expression_utility,
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
