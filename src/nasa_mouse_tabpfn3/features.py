"""Fold-local feature selection for high-dimensional OSDR RNA-seq matrices."""

from __future__ import annotations

from dataclasses import dataclass

from nasa_mouse_glare.io import require_import


@dataclass
class FeatureSelection:
    """Selected feature positions and diagnostics."""

    feature_mode: str
    indices: object
    genes: list[str]
    variance: object
    prevalence: object
    univariate_rank: object


def select_features(
    x_train,
    y_train,
    genes: list[str],
    *,
    feature_mode: str,
    min_expr_fraction: float = 0.05,
    hvg_top_n: int = 2000,
    max_features: int = 0,
) -> FeatureSelection:
    """Select features using only the training fold.

    `all_expressed` keeps all genes expressed in enough training samples, while
    `hvg` keeps the highest-variance genes after the same expression filter.
    The optional `max_features` cap is a runtime guard, not part of the default
    analysis.
    """

    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    x_train = np.asarray(x_train, dtype="float32")
    y_train = np.asarray(y_train, dtype="int64")
    if x_train.ndim != 2:
        raise ValueError("x_train must be a 2D matrix")

    prevalence = (x_train > 0).mean(axis=0)
    variance = x_train.var(axis=0)
    expressed = prevalence >= float(min_expr_fraction)
    expressed &= variance > 0
    candidate = np.flatnonzero(expressed)
    if candidate.size == 0:
        candidate = np.flatnonzero(variance > 0)
    if candidate.size == 0:
        candidate = np.arange(x_train.shape[1])

    mode = feature_mode.strip().lower()
    if mode in {"all", "all_expressed", "all_genes"}:
        selected = candidate
    elif mode in {"hvg", "highly_variable"}:
        order = np.argsort(-variance[candidate], kind="stable")
        selected = candidate[order[: min(int(hvg_top_n), candidate.size)]]
    else:
        raise ValueError(f"Unsupported feature mode: {feature_mode}")

    if max_features and max_features > 0 and selected.size > max_features:
        order = np.argsort(-variance[selected], kind="stable")
        selected = selected[order[: int(max_features)]]

    selected = np.asarray(selected, dtype="int64")
    ranks = univariate_feature_rank(x_train[:, selected], y_train)
    return FeatureSelection(
        feature_mode=mode,
        indices=selected,
        genes=[genes[int(index)] for index in selected.tolist()],
        variance=variance[selected],
        prevalence=prevalence[selected],
        univariate_rank=ranks,
    )


def univariate_feature_rank(x_train, y_train):
    """Rank selected features by absolute standardized FLT-GC difference."""

    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    x_train = np.asarray(x_train, dtype="float32")
    y_train = np.asarray(y_train, dtype="int64")
    flight = x_train[y_train == 1]
    ground = x_train[y_train == 0]
    if flight.size == 0 or ground.size == 0:
        return np.arange(x_train.shape[1], dtype="int64")
    mean_diff = flight.mean(axis=0) - ground.mean(axis=0)
    pooled = x_train.std(axis=0)
    pooled = np.where(pooled < 1e-6, 1.0, pooled)
    scores = np.abs(mean_diff / pooled)
    return np.argsort(-scores, kind="stable").astype("int64")


def top_univariate_candidate_indices(selection: FeatureSelection, top_n: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    ranks = np.asarray(selection.univariate_rank, dtype="int64")
    if top_n and top_n > 0:
        ranks = ranks[: min(int(top_n), ranks.size)]
    return ranks

