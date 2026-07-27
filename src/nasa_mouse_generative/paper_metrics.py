"""Paper-aligned distribution metrics for synthetic expression matrices."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import sqrtm
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA


@dataclass(frozen=True)
class PaperMetricThresholds:
    """Independent quality requirements derived from the DDIM paper benchmark."""

    correlation_matrix_agreement_minimum: float = 0.98
    precision_minimum: float = 0.95
    recall_minimum: float = 0.85
    f1_minimum: float = 0.90
    adversarial_accuracy_minimum: float = 0.40
    adversarial_accuracy_maximum: float = 0.60
    frechet_ratio_to_real_split_p95_maximum: float = 1.0


DEFAULT_PAPER_THRESHOLDS = PaperMetricThresholds()


def safe_correlation(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=np.float64).ravel()
    second = np.asarray(second, dtype=np.float64).ravel()
    finite = np.isfinite(first) & np.isfinite(second)
    if finite.sum() < 2:
        return float("nan")
    first = first[finite]
    second = second[finite]
    if np.std(first) == 0 or np.std(second) == 0:
        return float("nan")
    return float(np.corrcoef(first, second)[0, 1])


def correlation_matrix_agreement(real: np.ndarray, synthetic: np.ndarray) -> float:
    """Return the paper's gamma agreement between gene-correlation matrices."""

    with np.errstate(divide="ignore", invalid="ignore"):
        real_correlation = np.corrcoef(real, rowvar=False)
        synthetic_correlation = np.corrcoef(synthetic, rowvar=False)
    upper = np.triu_indices(real.shape[1], k=1)
    return safe_correlation(
        real_correlation[upper], synthetic_correlation[upper]
    )


def correlation_bootstrap_reference(
    real: np.ndarray, *, repeats: int = 20, seed: int = 2026
) -> dict[str, float | int]:
    """Estimate finite-sample gamma variation from same-size real bootstraps."""

    real = np.asarray(real, dtype=np.float32)
    if len(real) < 3 or int(repeats) < 1:
        return {
            "correlation_real_bootstrap_repeats": 0,
            "correlation_real_bootstrap_p05": float("nan"),
            "correlation_real_bootstrap_median": float("nan"),
            "correlation_real_bootstrap_p95": float("nan"),
        }
    rng = np.random.default_rng(seed)
    values = np.asarray(
        [
            correlation_matrix_agreement(
                real,
                real[rng.choice(len(real), len(real), replace=True)],
            )
            for _ in range(int(repeats))
        ],
        dtype=float,
    )
    values = values[np.isfinite(values)]
    if not len(values):
        return {
            "correlation_real_bootstrap_repeats": 0,
            "correlation_real_bootstrap_p05": float("nan"),
            "correlation_real_bootstrap_median": float("nan"),
            "correlation_real_bootstrap_p95": float("nan"),
        }
    return {
        "correlation_real_bootstrap_repeats": int(len(values)),
        "correlation_real_bootstrap_p05": float(np.quantile(values, 0.05)),
        "correlation_real_bootstrap_median": float(np.median(values)),
        "correlation_real_bootstrap_p95": float(np.quantile(values, 0.95)),
    }


def precision_recall(
    real: np.ndarray, synthetic: np.ndarray, *, neighbors: int
) -> dict[str, float]:
    """Compute the paper's k-NN hypersphere precision and recall."""

    if neighbors < 1:
        raise ValueError("neighbors must be positive")
    if len(real) <= neighbors or len(synthetic) <= neighbors:
        return {
            "precision": float("nan"),
            "recall": float("nan"),
            "f1": float("nan"),
        }
    real_distances = cdist(real, real, metric="sqeuclidean")
    synthetic_distances = cdist(synthetic, synthetic, metric="sqeuclidean")
    cross_distances = cdist(real, synthetic, metric="sqeuclidean")
    # Index zero is each sample itself, matching the upstream implementation.
    real_radius = np.partition(real_distances, neighbors, axis=1)[:, neighbors]
    synthetic_radius = np.partition(
        synthetic_distances, neighbors, axis=1
    )[:, neighbors]
    precision = float(
        np.mean((cross_distances <= real_radius[:, None]).any(axis=0))
    )
    recall = float(
        np.mean((cross_distances <= synthetic_radius[None, :]).any(axis=1))
    )
    f1 = (
        float(2.0 * precision * recall / (precision + recall))
        if precision + recall > 0
        else 0.0
    )
    return {"precision": precision, "recall": recall, "f1": f1}


def adversarial_accuracy(real: np.ndarray, synthetic: np.ndarray) -> float:
    """Compute symmetric 1-NN real-versus-synthetic adversarial accuracy."""

    if len(real) < 2 or len(synthetic) < 2:
        return float("nan")
    real_distances = cdist(real, real, metric="sqeuclidean")
    synthetic_distances = cdist(synthetic, synthetic, metric="sqeuclidean")
    np.fill_diagonal(real_distances, np.inf)
    np.fill_diagonal(synthetic_distances, np.inf)
    cross_distances = cdist(real, synthetic, metric="sqeuclidean")
    return float(
        0.5
        * (
            np.mean(cross_distances.min(axis=1) > real_distances.min(axis=1))
            + np.mean(
                cross_distances.min(axis=0) > synthetic_distances.min(axis=1)
            )
        )
    )


def frechet_distance(first: np.ndarray, second: np.ndarray) -> float:
    first_mean = first.mean(axis=0)
    second_mean = second.mean(axis=0)
    first_covariance = np.atleast_2d(np.cov(first, rowvar=False))
    second_covariance = np.atleast_2d(np.cov(second, rowvar=False))
    covariance_mean = sqrtm(first_covariance @ second_covariance)
    if np.iscomplexobj(covariance_mean):
        covariance_mean = covariance_mean.real
    difference = first_mean - second_mean
    return float(
        difference @ difference
        + np.trace(
            first_covariance + second_covariance - 2.0 * covariance_mean
        )
    )


def pca_frechet_with_real_reference(
    real: np.ndarray,
    synthetic: np.ndarray,
    *,
    dimensions: int = 50,
    reference_repeats: int = 20,
    seed: int = 2026,
) -> dict[str, float | int | str]:
    """Compare PCA FD with finite-sample FD between independent real halves."""

    half = len(real) // 2
    components = min(int(dimensions), real.shape[1], half - 1)
    if components < 1 or reference_repeats < 1:
        return {
            "frechet_pca": float("nan"),
            "frechet_real_split_median": float("nan"),
            "frechet_real_split_p95": float("nan"),
            "frechet_ratio_to_real_split_p95": float("nan"),
            "frechet_embedding_dimensions": 0,
            "frechet_embedding": "train-fitted PCA",
        }
    pca = PCA(n_components=components, random_state=int(seed)).fit(real)
    real_embedding = pca.transform(real)
    synthetic_embedding = pca.transform(synthetic)
    observed = frechet_distance(real_embedding, synthetic_embedding)
    rng = np.random.default_rng(seed)
    references: list[float] = []
    for _ in range(int(reference_repeats)):
        indices = rng.permutation(len(real))
        references.append(
            frechet_distance(
                real_embedding[indices[:half]],
                real_embedding[indices[half : 2 * half]],
            )
        )
    p95 = float(np.quantile(references, 0.95))
    return {
        "frechet_pca": observed,
        "frechet_real_split_median": float(np.median(references)),
        "frechet_real_split_p95": p95,
        "frechet_ratio_to_real_split_p95": (
            float(observed / p95) if p95 > 0 else float("nan")
        ),
        "frechet_embedding_dimensions": int(components),
        "frechet_embedding": "real-fitted PCA",
    }


def paper_distribution_metrics(
    real: np.ndarray,
    synthetic: np.ndarray,
    *,
    max_samples: int = 2000,
    neighbors: int | None = None,
    adversarial_max_samples: int = 2048,
    correlation_reference_repeats: int = 20,
    seed: int = 2026,
) -> dict[str, float | int | str]:
    """Evaluate one real/synthetic pair with the paper's unsupervised metrics."""

    real = np.asarray(real, dtype=np.float32)
    synthetic = np.asarray(synthetic, dtype=np.float32)
    if real.ndim != 2 or synthetic.ndim != 2:
        raise ValueError("real and synthetic expression must be two-dimensional")
    if real.shape[1] != synthetic.shape[1]:
        raise ValueError("real and synthetic expression must have the same genes")
    count = min(len(real), len(synthetic), int(max_samples))
    if count < 3:
        return {
            "metric_samples": int(count),
            "correlation_matrix_agreement": float("nan"),
            "precision": float("nan"),
            "recall": float("nan"),
            "f1": float("nan"),
            "adversarial_accuracy": float("nan"),
            "frechet_pca": float("nan"),
            "frechet_real_split_p95": float("nan"),
            "frechet_ratio_to_real_split_p95": float("nan"),
        }
    rng = np.random.default_rng(seed)
    if count != len(real) or count != len(synthetic):
        if len(real) == len(synthetic):
            indices = rng.choice(len(real), count, replace=False)
            real = real[indices]
            synthetic = synthetic[indices]
        else:
            real = real[rng.choice(len(real), count, replace=False)]
            synthetic = synthetic[
                rng.choice(len(synthetic), count, replace=False)
            ]
    if neighbors is None:
        neighbors = 10 if real.shape[1] <= 1000 else 50
    adversarial_count = min(count, int(adversarial_max_samples))
    if adversarial_count < count:
        adversarial_indices = rng.choice(
            count, adversarial_count, replace=False
        )
        adversarial_real = real[adversarial_indices]
        adversarial_synthetic = synthetic[adversarial_indices]
    else:
        adversarial_real = real
        adversarial_synthetic = synthetic
    result: dict[str, float | int | str] = {
        "metric_samples": int(count),
        "adversarial_metric_samples": int(adversarial_count),
        "genes": int(real.shape[1]),
        "neighbors": int(neighbors),
        "correlation_matrix_agreement": correlation_matrix_agreement(
            real, synthetic
        ),
        "adversarial_accuracy": adversarial_accuracy(
            adversarial_real, adversarial_synthetic
        ),
    }
    result.update(precision_recall(real, synthetic, neighbors=int(neighbors)))
    result.update(
        correlation_bootstrap_reference(
            real,
            repeats=int(correlation_reference_repeats),
            seed=int(seed) + 17,
        )
    )
    result.update(
        pca_frechet_with_real_reference(real, synthetic, seed=int(seed))
    )
    return result


def paper_metric_selection(
    metrics: dict[str, object],
    *,
    thresholds: PaperMetricThresholds = DEFAULT_PAPER_THRESHOLDS,
    finite_sample_calibrated: bool = False,
) -> dict[str, object]:
    """Require every paper-aligned quality metric to pass independently."""

    observed = {
        "correlation_matrix_agreement": float(
            metrics.get("correlation_matrix_agreement", float("nan"))
        ),
        "precision": float(metrics.get("precision", float("nan"))),
        "recall": float(metrics.get("recall", float("nan"))),
        "f1": float(metrics.get("f1", float("nan"))),
        "adversarial_accuracy": float(
            metrics.get("adversarial_accuracy", float("nan"))
        ),
        "frechet_ratio_to_real_split_p95": float(
            metrics.get(
                "frechet_ratio_to_real_split_p95", float("nan")
            )
        ),
    }
    correlation_minimum = thresholds.correlation_matrix_agreement_minimum
    correlation_reference = float(
        metrics.get("correlation_real_bootstrap_p05", float("nan"))
    )
    if finite_sample_calibrated and np.isfinite(correlation_reference):
        correlation_minimum = min(correlation_minimum, correlation_reference)
    requirements = {
        "correlation_matrix_agreement": {
            "minimum": correlation_minimum,
            "absolute_paper_minimum": (
                thresholds.correlation_matrix_agreement_minimum
            ),
            "real_bootstrap_p05": correlation_reference,
            "policy": (
                "minimum_of_absolute_paper_and_real_bootstrap_p05"
                if finite_sample_calibrated
                else "absolute_paper"
            ),
        },
        "precision": {"minimum": thresholds.precision_minimum},
        "recall": {"minimum": thresholds.recall_minimum},
        "f1": {"minimum": thresholds.f1_minimum},
        "adversarial_accuracy": {
            "minimum": thresholds.adversarial_accuracy_minimum,
            "maximum": thresholds.adversarial_accuracy_maximum,
        },
        "frechet_ratio_to_real_split_p95": {
            "maximum": thresholds.frechet_ratio_to_real_split_p95_maximum
        },
    }
    checks: dict[str, bool] = {}
    for name, requirement in requirements.items():
        value = observed[name]
        checks[name] = bool(
            np.isfinite(value)
            and value >= float(requirement.get("minimum", -np.inf))
            and value <= float(requirement.get("maximum", np.inf))
        )
    return {
        "passed": bool(all(checks.values())),
        "checks": checks,
        "observed": observed,
        "requirements": requirements,
        "failed_metrics": [name for name, passed in checks.items() if not passed],
        "selection_rule": (
            "all_metrics_must_pass_independently_with_finite_sample_corr"
            if finite_sample_calibrated
            else "all_absolute_paper_metrics_must_pass_independently"
        ),
    }
