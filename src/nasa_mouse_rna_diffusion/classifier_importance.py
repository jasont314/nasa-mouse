"""Compare per-arm permutation importance and linear SHAP contributions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata
from sklearn.preprocessing import StandardScaler

from .generated_feature_guidance import _fit_classifier, _reactome_enrichment
from .within_study_feature_stability import (
    METRICS,
    _arm_specs,
    _labels,
    _load_data,
    _metric_set,
    _muscle_group_analysis_data,
    _scaled_views,
    _valid_nested_split,
)


ARM_ORDER = (
    "real_only",
    "generated_only",
    "real_plus_generated",
    "guided_real_only",
    "guided_low_weight",
)
ARM_LABELS = {
    "real_only": "Real only",
    "generated_only": "Synthetic only",
    "real_plus_generated": "Real + synthetic",
    "guided_real_only": "Guided; real fit",
    "guided_low_weight": "Guided; 5% synthetic",
}


def _safe_spearman(left: Iterable[float], right: Iterable[float]) -> float:
    frame = pd.DataFrame({"left": left, "right": right}).dropna()
    if len(frame) < 3:
        return float("nan")
    left_rank = frame["left"].rank(method="average").to_numpy(dtype=float)
    right_rank = frame["right"].rank(method="average").to_numpy(dtype=float)
    if np.std(left_rank) == 0.0 or np.std(right_rank) == 0.0:
        return float("nan")
    return float(np.corrcoef(left_rank, right_rank)[0, 1])


def _recenter_evaluation_draw(
    synthetic_train: np.ndarray,
    synthetic_evaluation: np.ndarray,
    real_train: np.ndarray,
    train_labels: np.ndarray,
    evaluation_labels: np.ndarray,
) -> np.ndarray:
    result = np.asarray(synthetic_evaluation, dtype=np.float64).copy()
    for condition in (0, 1):
        train_mask = train_labels == condition
        evaluation_mask = evaluation_labels == condition
        if not train_mask.any() or not evaluation_mask.any():
            continue
        offset = (
            real_train[train_mask].mean(axis=0)
            - synthetic_train[train_mask].mean(axis=0)
        )
        result[evaluation_mask] += offset
    return result


def _training_background(
    arm: dict[str, Any],
    selected: np.ndarray,
    real_scaled: np.ndarray,
    synthetic_scaled: list[np.ndarray],
    recentered_scaled: list[np.ndarray],
) -> np.ndarray:
    training = str(arm["training"])
    real = real_scaled[:, selected]
    draws = recentered_scaled if bool(arm["recenter"]) else synthetic_scaled
    synthetic = np.concatenate([draw[:, selected] for draw in draws])
    if training == "real":
        return real.mean(axis=0)
    if training == "synthetic":
        return synthetic.mean(axis=0)
    if training != "combined":
        raise ValueError(f"Unsupported training mode: {training}")
    per_row_weight = float(arm["synthetic_weight"]) / len(draws)
    expression = np.concatenate((real, synthetic))
    weights = np.concatenate(
        (
            np.ones(len(real), dtype=float),
            np.full(len(synthetic), per_row_weight, dtype=float),
        )
    )
    return np.average(expression, axis=0, weights=weights)


def _fit_selected_arm(
    arm: dict[str, Any],
    selected: np.ndarray,
    *,
    regularization_c: float,
    real_scaled: np.ndarray,
    labels: np.ndarray,
    synthetic_scaled: list[np.ndarray],
    recentered_scaled: list[np.ndarray],
    seed: int,
):
    """Refit an arm using the exact feature indices saved by the nested run."""

    training = str(arm["training"])
    real = real_scaled[:, selected]
    draws = recentered_scaled if bool(arm["recenter"]) else synthetic_scaled
    synthetic = np.concatenate([draw[:, selected] for draw in draws])
    synthetic_labels = np.tile(labels, len(draws))
    if training == "real":
        return _fit_classifier(
            real,
            labels,
            regularization_c=regularization_c,
            seed=seed,
        )
    if training == "synthetic":
        return _fit_classifier(
            synthetic,
            synthetic_labels,
            regularization_c=regularization_c,
            seed=seed,
        )
    if training != "combined":
        raise ValueError(f"Unsupported training mode: {training}")
    synthetic_weight = float(arm["synthetic_weight"]) / len(draws)
    expression = np.concatenate((real, synthetic))
    combined_labels = np.concatenate((labels, synthetic_labels))
    sample_weight = np.concatenate(
        (
            np.ones(len(real), dtype=float),
            np.full(len(synthetic), synthetic_weight, dtype=float),
        )
    )
    return _fit_classifier(
        expression,
        combined_labels,
        regularization_c=regularization_c,
        seed=seed,
        sample_weight=sample_weight,
    )


def _permutation_rows(
    classifier,
    expression: np.ndarray,
    labels: np.ndarray,
    genes: list[str],
    symbols: dict[str, str],
    *,
    permutation_repeats: int,
    seed: int,
    blocks: list[np.ndarray] | None = None,
) -> tuple[pd.DataFrame, dict[str, float]]:
    expression = np.asarray(expression, dtype=float)
    labels = np.asarray(labels, dtype=np.int64)
    if expression.ndim != 2 or expression.shape[1] != len(genes):
        raise ValueError("Permutation expression and gene dimensions differ")
    if len(expression) != len(labels):
        raise ValueError("Permutation expression and labels differ")
    probability = classifier.predict_proba(expression)[:, 1]
    baseline = _metric_set(labels, probability)
    baseline_logit = np.asarray(classifier.decision_function(expression), dtype=float)
    coefficients = np.asarray(classifier.coef_[0], dtype=float)
    shuffle_blocks = blocks or [np.arange(len(expression), dtype=np.int64)]
    if sorted(np.concatenate(shuffle_blocks).tolist()) != list(range(len(expression))):
        raise ValueError("Permutation blocks must partition evaluation rows")
    rng = np.random.default_rng(int(seed))
    rows: list[dict[str, object]] = []
    for column, gene in enumerate(genes):
        original = expression[:, column]
        permuted = np.tile(original, (int(permutation_repeats), 1))
        for block in shuffle_blocks:
            orders = np.argsort(
                rng.random((int(permutation_repeats), len(block))),
                axis=1,
                kind="stable",
            )
            permuted[:, block] = original[block][orders]
        logits = baseline_logit.reshape(1, -1) + coefficients[column] * (
            permuted - original.reshape(1, -1)
        )
        permuted_metrics = _metric_matrix_from_logits(labels, logits)
        row: dict[str, object] = {
            "gene": gene,
            "symbol": symbols.get(gene, gene),
        }
        for metric in METRICS:
            array = baseline[metric] - permuted_metrics[metric]
            row[f"baseline_{metric}"] = float(baseline[metric])
            row[f"permutation_{metric}_mean"] = float(array.mean())
            row[f"permutation_{metric}_sd"] = float(
                array.std(ddof=1) if len(array) > 1 else 0.0
            )
            row[f"permutation_{metric}_positive_fraction"] = float(
                np.mean(array > 0.0)
            )
        rows.append(row)
    return pd.DataFrame(rows), baseline


def _metric_matrix_from_logits(
    labels: np.ndarray, logits: np.ndarray
) -> dict[str, np.ndarray]:
    """Evaluate rows of logits with tie-aware vectorized binary metrics."""

    labels = np.asarray(labels, dtype=np.int64)
    logits = np.asarray(logits, dtype=float)
    if logits.ndim == 1:
        logits = logits.reshape(1, -1)
    if logits.ndim != 2 or logits.shape[1] != len(labels):
        raise ValueError("Logit rows and labels differ")
    positive = labels == 1
    negative = labels == 0
    positive_count = int(positive.sum())
    negative_count = int(negative.sum())
    if positive_count == 0 or negative_count == 0:
        nan = np.full(logits.shape[0], np.nan, dtype=float)
        return {metric: nan.copy() for metric in METRICS}

    predicted = logits >= 0.0
    true_positive_rate = predicted[:, positive].mean(axis=1)
    true_negative_rate = (~predicted[:, negative]).mean(axis=1)
    balanced_accuracy = 0.5 * (true_positive_rate + true_negative_rate)

    ranks = rankdata(logits, method="average", axis=1)
    roc_auc = (
        ranks[:, positive].sum(axis=1)
        - positive_count * (positive_count + 1) / 2.0
    ) / (positive_count * negative_count)

    order = np.argsort(-logits, axis=1, kind="stable")
    sorted_scores = np.take_along_axis(logits, order, axis=1)
    sorted_labels = labels[order]
    cumulative_positive = np.cumsum(sorted_labels, axis=1)
    group_end = np.ones_like(sorted_scores, dtype=bool)
    group_end[:, :-1] = sorted_scores[:, :-1] != sorted_scores[:, 1:]
    cumulative_at_end = np.maximum.accumulate(
        np.where(group_end, cumulative_positive, 0), axis=1
    )
    previous_end = np.concatenate(
        (
            np.zeros((len(logits), 1), dtype=int),
            cumulative_at_end[:, :-1],
        ),
        axis=1,
    )
    group_positive = cumulative_positive - previous_end
    precision = cumulative_positive / np.arange(1, logits.shape[1] + 1)
    average_precision = (
        np.where(group_end, group_positive * precision, 0.0).sum(axis=1)
        / positive_count
    )
    return {
        "balanced_accuracy": balanced_accuracy,
        "roc_auc": roc_auc,
        "average_precision": average_precision,
    }


def _linear_shap_rows(
    classifier,
    expression: np.ndarray,
    labels: np.ndarray,
    background: np.ndarray,
    genes: list[str],
    symbols: dict[str, str],
) -> tuple[pd.DataFrame, float]:
    """Exact interventional SHAP values for a linear log-odds model."""

    expression = np.asarray(expression, dtype=float)
    background = np.asarray(background, dtype=float)
    coefficients = np.asarray(classifier.coef_[0], dtype=float)
    values = (expression - background.reshape(1, -1)) * coefficients.reshape(1, -1)
    expected_logit = float(classifier.intercept_[0] + background @ coefficients)
    reconstructed = expected_logit + values.sum(axis=1)
    observed = np.asarray(classifier.decision_function(expression), dtype=float)
    reconstruction_error = float(np.max(np.abs(reconstructed - observed)))
    flight = labels == 1
    ground = labels == 0
    rows = []
    for column, gene in enumerate(genes):
        contribution = values[:, column]
        flight_mean = float(contribution[flight].mean())
        ground_mean = float(contribution[ground].mean())
        rows.append(
            {
                "gene": gene,
                "symbol": symbols.get(gene, gene),
                "linear_shap_mean_absolute": float(np.mean(np.abs(contribution))),
                "linear_shap_mean_flight": flight_mean,
                "linear_shap_mean_ground_control": ground_mean,
                "linear_shap_flight_minus_ground": flight_mean - ground_mean,
                "linear_shap_reconstruction_max_error": reconstruction_error,
            }
        )
    return pd.DataFrame(rows), reconstruction_error


def _aggregate_importance(
    rows: pd.DataFrame, completed_repeats: pd.DataFrame
) -> pd.DataFrame:
    value_columns = [
        "classifier_coefficient",
        *[
            f"baseline_{metric}"
            for metric in METRICS
        ],
        *[
            f"permutation_{metric}_{suffix}"
            for metric in METRICS
            for suffix in ("mean", "sd", "positive_fraction")
        ],
        "linear_shap_mean_absolute",
        "linear_shap_mean_flight",
        "linear_shap_mean_ground_control",
        "linear_shap_flight_minus_ground",
        "linear_shap_reconstruction_max_error",
    ]
    grouped = rows.groupby(
        ["scope", "tissue", "arm", "domain", "gene", "symbol"],
        observed=True,
    )
    aggregate_spec: dict[str, tuple[str, object]] = {
        "selected_repeats": ("repeat", "nunique"),
        "median_classifier_coefficient": ("classifier_coefficient", "median"),
        "mean_absolute_classifier_coefficient": (
            "classifier_coefficient",
            lambda values: float(np.mean(np.abs(values))),
        ),
        "coefficient_sign_agreement": (
            "classifier_coefficient",
            lambda values: float(abs(np.mean(np.sign(values)))),
        ),
    }
    for column in value_columns:
        if column == "classifier_coefficient":
            continue
        aggregate_spec[column] = (column, "mean")
    for metric in METRICS:
        column = f"permutation_{metric}_mean"
        aggregate_spec[f"permutation_{metric}_positive_repeat_fraction"] = (
            column,
            lambda values: float(np.mean(np.asarray(values) > 0.0)),
        )
    table = grouped.agg(**aggregate_spec).reset_index()
    table = table.merge(
        completed_repeats,
        on=["scope", "tissue"],
        how="left",
        validate="many_to_one",
    )
    table["selection_frequency"] = (
        table["selected_repeats"] / table["completed_repeats"]
    )
    return table.sort_values(
        ["scope", "tissue", "arm", "domain", "permutation_roc_auc_mean"],
        ascending=[True, True, True, True, False],
        kind="stable",
    )


def _comparison_pattern(
    *,
    real_stable: bool,
    arm_stable: bool,
    coefficient_match: bool,
    real_importance: float,
    arm_real_importance: float,
    arm_real_positive_fraction: float,
    arm_synthetic_importance: float,
    arm_synthetic_positive_fraction: float,
) -> str:
    real_transfer = (
        arm_real_importance > 0.0 and arm_real_positive_fraction >= 0.5
    )
    synthetic_support = (
        arm_synthetic_importance > 0.0
        and arm_synthetic_positive_fraction >= 0.5
    )
    if arm_stable and not real_stable:
        if real_transfer:
            return "synthetic_emergent_real_transfer"
        if synthetic_support:
            return "synthetic_domain_only"
        return "synthetic_emergent_no_positive_importance"
    if arm_stable and real_stable:
        if not coefficient_match:
            return "shared_direction_conflict"
        if not real_transfer:
            if synthetic_support:
                return "shared_synthetic_domain_only"
            return "shared_no_positive_importance"
        if arm_real_importance >= real_importance:
            return "shared_reinforced"
        return "shared_attenuated"
    if real_stable and not arm_stable:
        return "real_only"
    return "unstable"


def _compare_arms(
    aggregate: pd.DataFrame,
    thresholds: dict[tuple[str, str], tuple[float, float]],
) -> pd.DataFrame:
    real_domain = aggregate.loc[aggregate["domain"].eq("real")]
    synthetic_domain = aggregate.loc[aggregate["domain"].eq("synthetic")]
    rows: list[dict[str, object]] = []
    for (scope, tissue), frame in real_domain.groupby(
        ["scope", "tissue"], sort=True, observed=True
    ):
        baseline = frame.loc[frame["arm"].eq("real_only")].set_index("gene")
        frequency_threshold, sign_threshold = thresholds[(str(scope), str(tissue))]
        for arm in ARM_ORDER[1:]:
            arm_real = frame.loc[frame["arm"].eq(arm)].set_index("gene")
            arm_synthetic = synthetic_domain.loc[
                synthetic_domain["scope"].eq(scope)
                & synthetic_domain["tissue"].eq(tissue)
                & synthetic_domain["arm"].eq(arm)
            ].set_index("gene")
            genes = sorted(set(baseline.index) | set(arm_real.index) | set(arm_synthetic.index))
            for gene in genes:
                real_row = baseline.loc[gene] if gene in baseline.index else None
                arm_row = arm_real.loc[gene] if gene in arm_real.index else None
                synthetic_row = (
                    arm_synthetic.loc[gene]
                    if gene in arm_synthetic.index
                    else None
                )
                real_frequency = (
                    float(real_row["selection_frequency"])
                    if real_row is not None
                    else 0.0
                )
                arm_frequency = (
                    float(arm_row["selection_frequency"])
                    if arm_row is not None
                    else 0.0
                )
                real_sign_agreement = (
                    float(real_row["coefficient_sign_agreement"])
                    if real_row is not None
                    else 0.0
                )
                arm_sign_agreement = (
                    float(arm_row["coefficient_sign_agreement"])
                    if arm_row is not None
                    else 0.0
                )
                real_stable = (
                    real_frequency >= frequency_threshold
                    and real_sign_agreement >= sign_threshold
                )
                arm_stable = (
                    arm_frequency >= frequency_threshold
                    and arm_sign_agreement >= sign_threshold
                )
                real_coefficient = (
                    float(real_row["median_classifier_coefficient"])
                    if real_row is not None
                    else float("nan")
                )
                arm_coefficient = (
                    float(arm_row["median_classifier_coefficient"])
                    if arm_row is not None
                    else float("nan")
                )
                real_importance = (
                    float(real_row["permutation_roc_auc_mean"])
                    if real_row is not None
                    else 0.0
                )
                arm_real_importance = (
                    float(arm_row["permutation_roc_auc_mean"])
                    if arm_row is not None
                    else 0.0
                )
                arm_real_positive = (
                    float(arm_row["permutation_roc_auc_positive_repeat_fraction"])
                    if arm_row is not None
                    else 0.0
                )
                arm_synthetic_importance = (
                    float(synthetic_row["permutation_roc_auc_mean"])
                    if synthetic_row is not None
                    else 0.0
                )
                arm_synthetic_positive = (
                    float(
                        synthetic_row[
                            "permutation_roc_auc_positive_repeat_fraction"
                        ]
                    )
                    if synthetic_row is not None
                    else 0.0
                )
                coefficient_match = bool(
                    np.isfinite(real_coefficient)
                    and np.isfinite(arm_coefficient)
                    and np.sign(real_coefficient) == np.sign(arm_coefficient)
                )
                symbol = ""
                for candidate in (arm_row, real_row, synthetic_row):
                    if candidate is not None:
                        symbol = str(candidate["symbol"])
                        break
                rows.append(
                    {
                        "scope": scope,
                        "tissue": tissue,
                        "arm": arm,
                        "gene": gene,
                        "symbol": symbol,
                        "minimum_selection_frequency": frequency_threshold,
                        "minimum_coefficient_sign_agreement": sign_threshold,
                        "real_selection_frequency": real_frequency,
                        "arm_selection_frequency": arm_frequency,
                        "real_coefficient_sign_agreement": real_sign_agreement,
                        "arm_coefficient_sign_agreement": arm_sign_agreement,
                        "real_stable": real_stable,
                        "arm_stable": arm_stable,
                        "real_coefficient": real_coefficient,
                        "arm_coefficient": arm_coefficient,
                        "coefficient_direction_match": coefficient_match,
                        "real_permutation_roc_auc": real_importance,
                        "arm_real_permutation_roc_auc": arm_real_importance,
                        "arm_real_positive_repeat_fraction": arm_real_positive,
                        "arm_synthetic_permutation_roc_auc": arm_synthetic_importance,
                        "arm_synthetic_positive_repeat_fraction": arm_synthetic_positive,
                        "arm_minus_real_permutation_roc_auc": (
                            arm_real_importance - real_importance
                        ),
                        "pattern": _comparison_pattern(
                            real_stable=real_stable,
                            arm_stable=arm_stable,
                            coefficient_match=coefficient_match,
                            real_importance=real_importance,
                            arm_real_importance=arm_real_importance,
                            arm_real_positive_fraction=arm_real_positive,
                            arm_synthetic_importance=arm_synthetic_importance,
                            arm_synthetic_positive_fraction=arm_synthetic_positive,
                        ),
                    }
                )
    return pd.DataFrame(rows)


def _arm_similarity(
    aggregate: pd.DataFrame, comparison: pd.DataFrame
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    real_domain = aggregate.loc[aggregate["domain"].eq("real")]
    synthetic_domain = aggregate.loc[aggregate["domain"].eq("synthetic")]
    for (scope, tissue, arm), frame in comparison.groupby(
        ["scope", "tissue", "arm"], sort=True, observed=True
    ):
        baseline_real = real_domain.loc[
            real_domain["scope"].eq(scope)
            & real_domain["tissue"].eq(tissue)
            & real_domain["arm"].eq("real_only")
        ].set_index("gene")
        arm_real = real_domain.loc[
            real_domain["scope"].eq(scope)
            & real_domain["tissue"].eq(tissue)
            & real_domain["arm"].eq(arm)
        ].set_index("gene")
        arm_synthetic = synthetic_domain.loc[
            synthetic_domain["scope"].eq(scope)
            & synthetic_domain["tissue"].eq(tissue)
            & synthetic_domain["arm"].eq(arm)
        ].set_index("gene")
        domain_genes = sorted(set(arm_real.index) | set(arm_synthetic.index))
        arm_real_values = [
            float(arm_real.loc[gene, "permutation_roc_auc_mean"])
            if gene in arm_real.index
            else 0.0
            for gene in domain_genes
        ]
        arm_synthetic_values = [
            float(arm_synthetic.loc[gene, "permutation_roc_auc_mean"])
            if gene in arm_synthetic.index
            else 0.0
            for gene in domain_genes
        ]
        common = frame.loc[frame["real_stable"] & frame["arm_stable"]]
        same_direction = (
            float(common["coefficient_direction_match"].mean())
            if len(common)
            else float("nan")
        )
        real_top = set(
            baseline_real.sort_values("permutation_roc_auc_mean", ascending=False)
            .head(20)
            .index
        )
        arm_top = set(
            arm_real.sort_values("permutation_roc_auc_mean", ascending=False)
            .head(20)
            .index
        )
        rows.append(
            {
                "scope": scope,
                "tissue": tissue,
                "arm": arm,
                "real_vs_arm_importance_spearman": _safe_spearman(
                    frame["real_permutation_roc_auc"],
                    frame["arm_real_permutation_roc_auc"],
                ),
                "arm_real_vs_synthetic_domain_spearman": _safe_spearman(
                    arm_real_values, arm_synthetic_values
                ),
                "top20_overlap": len(real_top & arm_top),
                "shared_stable_genes": int(len(common)),
                "shared_coefficient_direction_agreement": same_direction,
                **{
                    f"pattern_{pattern}": int((frame["pattern"] == pattern).sum())
                    for pattern in sorted(frame["pattern"].unique())
                },
            }
        )
    return pd.DataFrame(rows)


def _importance_enrichment(
    comparison: pd.DataFrame,
    *,
    background: list[str],
    gmt_path: Path,
    symbols: dict[str, str],
) -> pd.DataFrame:
    target_patterns = (
        "synthetic_emergent_real_transfer",
        "synthetic_domain_only",
        "shared_reinforced",
    )
    tables: list[pd.DataFrame] = []
    for (scope, tissue, arm, pattern), frame in comparison.loc[
        comparison["pattern"].isin(target_patterns)
    ].groupby(["scope", "tissue", "arm", "pattern"], observed=True):
        genes = frame["gene"].astype(str).tolist()
        if len(genes) < 2:
            continue
        table = _reactome_enrichment(genes, background, gmt_path, symbols)
        if table.empty:
            continue
        table.insert(0, "pattern", pattern)
        table.insert(0, "arm", arm)
        table.insert(0, "tissue", tissue)
        table.insert(0, "scope", scope)
        tables.append(table)
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def _selected_arm_tables(
    aggregate: pd.DataFrame,
    comparison: pd.DataFrame,
    choices: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    choice_columns = [
        "scope",
        "tissue",
        "selected_arm",
        "generated_arm_eligible_all_metrics",
    ]
    selected_importance = aggregate.merge(
        choices[choice_columns],
        on=["scope", "tissue"],
        how="inner",
        validate="many_to_one",
    )
    selected_importance = selected_importance.loc[
        selected_importance["arm"].eq(selected_importance["selected_arm"])
    ].copy()
    selected_comparison = comparison.merge(
        choices[choice_columns],
        on=["scope", "tissue"],
        how="inner",
        validate="many_to_one",
    )
    selected_comparison = selected_comparison.loc[
        selected_comparison["arm"].eq(selected_comparison["selected_arm"])
    ].copy()
    return selected_importance, selected_comparison


def _synthetic_informed_gene_importance(
    inventory_path: Path,
    selected_importance: pd.DataFrame,
    aggregate: pd.DataFrame,
) -> pd.DataFrame:
    inventory = pd.read_csv(inventory_path, sep="\t")
    scope_mapping = {
        "canonical_tissue": "tissue",
        "skeletal_muscle_group": "muscle_group",
    }
    inventory["scope"] = inventory["analysis_scope"].map(scope_mapping)
    if inventory["scope"].isna().any():
        unknown = sorted(inventory.loc[inventory["scope"].isna(), "analysis_scope"].unique())
        raise ValueError(f"Unknown synthetic-informed analysis scopes: {unknown}")
    available_units = selected_importance[["scope", "tissue"]].drop_duplicates()
    inventory = inventory.merge(
        available_units,
        on=["scope", "tissue"],
        how="inner",
        validate="many_to_one",
    )

    selected = selected_importance.loc[
        selected_importance["domain"].eq("real")
    ].copy()
    selected["permutation_rank_within_unit"] = selected.groupby(
        ["scope", "tissue"], observed=True
    )["permutation_roc_auc_mean"].rank(method="min", ascending=False)
    selected_columns = [
        "scope",
        "tissue",
        "gene",
        "selected_arm",
        "selection_frequency",
        "coefficient_sign_agreement",
        "median_classifier_coefficient",
        "permutation_balanced_accuracy_mean",
        "permutation_balanced_accuracy_positive_repeat_fraction",
        "permutation_roc_auc_mean",
        "permutation_roc_auc_positive_repeat_fraction",
        "permutation_average_precision_mean",
        "permutation_average_precision_positive_repeat_fraction",
        "linear_shap_mean_absolute",
        "linear_shap_flight_minus_ground",
        "permutation_rank_within_unit",
    ]
    selected = selected[selected_columns].rename(
        columns={
            column: f"selected_arm_{column}"
            for column in selected_columns
            if column not in {"scope", "tissue", "gene", "selected_arm"}
        }
    )
    baseline_columns = [
        "scope",
        "tissue",
        "gene",
        "selection_frequency",
        "coefficient_sign_agreement",
        "median_classifier_coefficient",
        "permutation_roc_auc_mean",
        "permutation_roc_auc_positive_repeat_fraction",
        "linear_shap_flight_minus_ground",
    ]
    baseline = aggregate.loc[
        aggregate["domain"].eq("real") & aggregate["arm"].eq("real_only"),
        baseline_columns,
    ].rename(
        columns={
            column: f"real_only_{column}"
            for column in baseline_columns
            if column not in {"scope", "tissue", "gene"}
        }
    )
    table = inventory.merge(
        selected,
        on=["scope", "tissue", "gene"],
        how="left",
        validate="one_to_one",
    ).merge(
        baseline,
        on=["scope", "tissue", "gene"],
        how="left",
        validate="one_to_one",
    )
    if table["selected_arm"].isna().any():
        missing = table.loc[
            table["selected_arm"].isna(), ["analysis_scope", "tissue", "gene"]
        ]
        raise RuntimeError(
            "Synthetic-informed genes are missing selected-arm importance: "
            f"{missing.head().to_dict(orient='records')}"
        )
    table["positive_held_out_real_permutation_importance"] = (
        table["selected_arm_permutation_roc_auc_mean"].gt(0.0)
        & table[
            "selected_arm_permutation_roc_auc_positive_repeat_fraction"
        ].ge(0.5)
    )
    table["positive_linear_shap_flt_gc_separation"] = table[
        "selected_arm_linear_shap_flight_minus_ground"
    ].gt(0.0)
    table["permutation_and_shap_supported"] = (
        table["positive_held_out_real_permutation_importance"]
        & table["positive_linear_shap_flt_gc_separation"]
    )
    table["selected_minus_real_only_permutation_roc_auc"] = (
        table["selected_arm_permutation_roc_auc_mean"]
        - table["real_only_permutation_roc_auc_mean"]
    )
    return table


def _plot_unit_importance(
    aggregate: pd.DataFrame, output: Path, *, scope: str, tissue: str
) -> None:
    subset = aggregate.loc[
        aggregate["scope"].eq(scope) & aggregate["tissue"].eq(tissue)
    ].copy()
    if subset.empty:
        return
    real = subset.loc[subset["domain"].eq("real")]
    score = (
        real.groupby(["gene", "symbol"], observed=True)["permutation_roc_auc_mean"]
        .apply(lambda values: float(np.max(np.abs(values))))
        .sort_values(ascending=False)
    )
    top_pairs = list(score.head(20).index)
    if not top_pairs:
        return
    top_genes = [pair[0] for pair in top_pairs]
    gene_labels = [pair[1] or pair[0] for pair in top_pairs]
    panels = [
        ("real", "permutation_roc_auc_mean", "Held-out real\nAUROC permutation drop"),
        (
            "synthetic",
            "permutation_roc_auc_mean",
            "Held-out synthetic\nAUROC permutation drop",
        ),
        (
            "real",
            "linear_shap_flight_minus_ground",
            "Held-out real\nlinear SHAP FLT-GC",
        ),
    ]
    figure, axes = plt.subplots(1, 3, figsize=(15.5, 8.2), sharey=True)
    for axis, (domain, value, title) in zip(axes, panels):
        frame = subset.loc[subset["domain"].eq(domain)]
        pivot = frame.pivot_table(index="gene", columns="arm", values=value)
        matrix = pivot.reindex(index=top_genes, columns=ARM_ORDER).fillna(0.0)
        limit = float(np.nanmax(np.abs(matrix.to_numpy())))
        limit = max(limit, 1e-6)
        image = axis.imshow(
            matrix,
            aspect="auto",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
        )
        axis.set_title(title, fontsize=11, weight="bold")
        axis.set_xticks(
            np.arange(len(ARM_ORDER)),
            [ARM_LABELS[arm] for arm in ARM_ORDER],
            rotation=45,
            ha="right",
            fontsize=8,
        )
        axis.set_yticks(np.arange(len(top_genes)), gene_labels, fontsize=8)
        figure.colorbar(image, ax=axis, fraction=0.035, pad=0.02)
    axes[0].set_ylabel("Gene")
    figure.suptitle(
        f"{tissue.replace('_', ' ').title()}: classifier importance by arm",
        fontsize=14,
        weight="bold",
    )
    figure.subplots_adjust(left=0.10, right=0.95, top=0.90, bottom=0.20, wspace=0.28)
    directory = output / scope / tissue
    directory.mkdir(parents=True, exist_ok=True)
    figure.savefig(directory / "classifier_importance_heatmaps.png", dpi=220)
    figure.savefig(directory / "classifier_importance_heatmaps.pdf")
    plt.close(figure)


def _plot_domain_similarity(similarity: pd.DataFrame, output: Path) -> None:
    if similarity.empty:
        return
    labels = similarity["scope"].astype(str) + ":" + similarity["tissue"].astype(str)
    frame = similarity.assign(unit=labels).pivot_table(
        index="unit",
        columns="arm",
        values="arm_real_vs_synthetic_domain_spearman",
    )
    frame = frame.reindex(columns=ARM_ORDER[1:]).sort_index()
    figure, axis = plt.subplots(figsize=(8.5, max(5.0, 0.28 * len(frame))))
    image = axis.imshow(frame.fillna(0.0), aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    axis.set_xticks(
        np.arange(len(frame.columns)),
        [ARM_LABELS[arm] for arm in frame.columns],
        rotation=35,
        ha="right",
    )
    axis.set_yticks(np.arange(len(frame)), frame.index, fontsize=7)
    axis.set_title("Real versus synthetic-domain permutation-importance correlation")
    figure.colorbar(image, ax=axis, label="Spearman correlation")
    figure.tight_layout()
    figure.savefig(output / "arm_domain_importance_correlation.png", dpi=220)
    figure.savefig(output / "arm_domain_importance_correlation.pdf")
    plt.close(figure)


def _write_readme(output: Path, summary: dict[str, Any]) -> None:
    text = f"""# Classifier importance analysis

This analysis reconstructs the five selected tissue-specific classifiers from
their saved nested-development hyperparameters. It does not retrain the DDIM or
rerun classifier hyperparameter selection.

For each outer split and arm, gene permutation importance is evaluated on held-out
real profiles and on the corresponding held-out synthetic profiles. Synthetic
permutations are performed independently inside each DDIM draw. Balanced accuracy,
AUROC, and average precision drops are reported separately.

The `linear_shap_*` columns are exact interventional SHAP contributions for the
standardized logistic model: `coefficient * (value - training background)`. Their
sum reconstructs each sample's log-odds relative to the training-background
expectation.

## Outputs

- `importance_by_repeat.tsv.gz`: per-repeat, per-gene permutation and SHAP results.
- `importance_summary.tsv.gz`: importance aggregated across repeated outer splits.
- `arm_vs_real_gene_comparison.tsv.gz`: synthetic-arm patterns relative to real-only.
- `selected_arm_importance.tsv.gz`: importance for each unit's retained utility arm.
- `selected_arm_vs_real_gene_comparison.tsv.gz`: retained synthetic arms only.
- `synthetic_informed_bh_fdr_gene_importance.tsv`: permutation and SHAP crosswalk
  for the manuscript's promoted and reinforced BH-FDR genes.
- `tissue_arm_similarity.tsv`: arm and domain rank correlations and top-20 overlap.
- `reactome_importance_enrichment.tsv.gz`: exploratory enrichment of importance patterns.
- `top_importance_genes.tsv`: top genes by held-out-real AUROC permutation drop.
- `<scope>/<tissue>/classifier_importance_heatmaps.png`: per-tissue arm comparison.

## Scope

- Completed analysis units: {summary['completed_units']}
- Permutation repeats per fitted classifier: {summary['permutation_repeats']}
- Reconstructed fitted classifiers: {summary['fitted_classifiers']}
- Neural-network retraining: no

## Interpretation

Synthetic-only importance describes the FLT/GC pattern encoded by the generator.
Importance on held-out real profiles asks whether that pattern transfers back to
observed expression. It is not independent biological evidence because the fixed
DDIM was trained before these nested splits. Gene correlation can also distribute
or suppress marginal permutation importance. Biological claims still require
observed OSDR effect direction and FDR support.
"""
    (output / "README.md").write_text(text, encoding="utf-8")


def _run_source(
    source: dict[str, Any],
    *,
    permutation_repeats: int,
    run_seed: int,
    coefficient_tolerance: float,
    metric_tolerance: float,
    units_override: set[str] | None,
) -> tuple[
    pd.DataFrame,
    dict[str, tuple[float, float]],
    list[str],
    dict[str, str],
    Path,
]:
    scope = str(source["scope"])
    workflow_path = Path(source["workflow_config"])
    workflow = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    analysis_dir = Path(source["analysis_dir"])
    metrics = pd.read_csv(analysis_dir / "nested_arm_metrics.tsv.gz", sep="\t")
    saved_features = pd.read_csv(
        analysis_dir / "nested_selected_features.tsv.gz", sep="\t"
    )
    inventory = pd.read_csv(analysis_dir / "tissue_inventory.tsv", sep="\t")
    data = _load_data(workflow)
    if scope == "muscle_group":
        data, _ = _muscle_group_analysis_data(
            data, workflow["analysis"].get("muscle_groups")
        )
    elif scope != "tissue":
        raise ValueError(f"Unsupported scope: {scope}")
    unit_order = inventory["tissue"].astype(str).tolist()
    completed = set(metrics["tissue"].astype(str).unique())
    if units_override is not None:
        completed &= units_override
    arms = _arm_specs(workflow)
    if tuple(arms) != ARM_ORDER:
        raise ValueError(
            f"Saved workflow arm order differs from classifier-importance arms: {tuple(arms)}"
        )
    arm_offsets = {arm: index for index, arm in enumerate(arms)}
    base_seed = int(workflow["run"]["seed"])
    gene_lookup = {gene: index for index, gene in enumerate(data.genes)}
    result_tables: list[pd.DataFrame] = []
    completed_repeats: dict[str, int] = {}
    for tissue in unit_order:
        if tissue not in completed:
            continue
        print(f"[classifier-importance] {scope}:{tissue}", flush=True)
        tissue_position = unit_order.index(tissue)
        tissue_seed = base_seed + tissue_position * 100_000
        mask = data.development_samples["tissue"].astype(str).eq(tissue)
        positions = np.flatnonzero(mask.to_numpy())
        real = data.development_expression[positions]
        samples = data.development_samples.loc[mask].reset_index(drop=True)
        labels = _labels(samples)
        synthetic = [draw[positions] for draw in data.synthetic_draws.values()]
        tissue_metrics = metrics.loc[metrics["tissue"].astype(str).eq(tissue)]
        repeats = sorted(tissue_metrics["repeat"].astype(int).unique())
        completed_repeats[tissue] = len(repeats)
        for repeat in repeats:
            split = _valid_nested_split(
                samples,
                outer_fraction=float(workflow["analysis"]["outer_fraction"]),
                inner_fraction=float(workflow["analysis"]["inner_fraction"]),
                seed=tissue_seed + repeat,
            )
            if split is None:
                raise RuntimeError(f"Could not reconstruct split {scope}:{tissue}:{repeat}")
            _, _, outer_test = split
            outer_train = np.setdiff1d(np.arange(len(samples)), outer_test)
            scaler = StandardScaler().fit(real[outer_train])
            real_train_scaled, synthetic_train_scaled, recentered_train_scaled = (
                _scaled_views(
                    real[outer_train],
                    labels[outer_train],
                    [draw[outer_train] for draw in synthetic],
                    scaler,
                )
            )
            real_evaluation_scaled = scaler.transform(real[outer_test])
            raw_synthetic_evaluation = [draw[outer_test] for draw in synthetic]
            for arm in ARM_ORDER:
                candidate_rows = tissue_metrics.loc[
                    tissue_metrics["repeat"].astype(int).eq(repeat)
                    & tissue_metrics["arm"].eq(arm)
                ]
                if len(candidate_rows) != 1:
                    raise RuntimeError(
                        f"Expected one saved candidate for {scope}:{tissue}:{repeat}:{arm}"
                    )
                candidate = candidate_rows.iloc[0]
                fit_seed = tissue_seed + repeat * 100 + arm_offsets[arm]
                feature_rows = saved_features.loc[
                    saved_features["tissue"].astype(str).eq(tissue)
                    & saved_features["repeat"].astype(int).eq(repeat)
                    & saved_features["arm"].eq(arm)
                ]
                if len(feature_rows) != int(candidate["feature_count"]):
                    raise RuntimeError(
                        "Saved feature count differs from selected candidate for "
                        f"{scope}:{tissue}:{repeat}:{arm}"
                    )
                missing_genes = sorted(
                    set(feature_rows["gene"].astype(str)).difference(gene_lookup)
                )
                if missing_genes:
                    raise RuntimeError(
                        f"Saved genes are absent from the panel: {missing_genes[:5]}"
                    )
                selected = np.asarray(
                    [gene_lookup[str(gene)] for gene in feature_rows["gene"]],
                    dtype=np.int64,
                )
                classifier = _fit_selected_arm(
                    arms[arm],
                    selected,
                    regularization_c=float(candidate["regularization_c"]),
                    real_scaled=real_train_scaled,
                    labels=labels[outer_train],
                    synthetic_scaled=synthetic_train_scaled,
                    recentered_scaled=recentered_train_scaled,
                    seed=fit_seed,
                )
                saved_coefficients = feature_rows[
                    "classifier_coefficient"
                ].to_numpy(dtype=float)
                coefficient_error = float(
                    np.max(np.abs(classifier.coef_[0] - saved_coefficients))
                )
                if coefficient_error > coefficient_tolerance:
                    raise RuntimeError(
                        "Refitted coefficients do not reproduce the saved classifier for "
                        f"{scope}:{tissue}:{repeat}:{arm}; max error={coefficient_error:.3g}, "
                        f"tolerance={coefficient_tolerance:.3g}"
                    )
                # Preserve the exact fitted feature weights being interpreted. The
                # independently refitted intercept is retained and checked below.
                classifier.coef_[0] = saved_coefficients.copy()
                genes = [data.genes[index] for index in selected]
                background = _training_background(
                    arms[arm],
                    selected,
                    real_train_scaled,
                    synthetic_train_scaled,
                    recentered_train_scaled,
                )
                real_evaluation = real_evaluation_scaled[:, selected]
                real_permutation, real_baseline = _permutation_rows(
                    classifier,
                    real_evaluation,
                    labels[outer_test],
                    genes,
                    data.symbols,
                    permutation_repeats=permutation_repeats,
                    seed=run_seed + tissue_position * 1_000_000 + repeat * 10_000 + arm_offsets[arm] * 100,
                )
                metric_error = float(
                    max(
                        abs(float(real_baseline[metric]) - float(candidate[metric]))
                        for metric in METRICS
                    )
                )
                if metric_error > metric_tolerance:
                    raise RuntimeError(
                        "Refitted classifier does not reproduce saved held-out metrics for "
                        f"{scope}:{tissue}:{repeat}:{arm}; max error={metric_error:.3g}, "
                        f"tolerance={metric_tolerance:.3g}"
                    )
                real_shap, _ = _linear_shap_rows(
                    classifier,
                    real_evaluation,
                    labels[outer_test],
                    background,
                    genes,
                    data.symbols,
                )
                real_table = real_permutation.merge(
                    real_shap, on=["gene", "symbol"], validate="one_to_one"
                )
                real_table["domain"] = "real"

                synthetic_evaluation_draws: list[np.ndarray] = []
                for draw_index, draw in enumerate(raw_synthetic_evaluation):
                    values = draw
                    if bool(arms[arm]["recenter"]):
                        values = _recenter_evaluation_draw(
                            synthetic[draw_index][outer_train],
                            draw,
                            real[outer_train],
                            labels[outer_train],
                            labels[outer_test],
                        )
                    synthetic_evaluation_draws.append(
                        scaler.transform(values)[:, selected]
                    )
                synthetic_evaluation = np.concatenate(synthetic_evaluation_draws)
                synthetic_labels = np.tile(
                    labels[outer_test], len(synthetic_evaluation_draws)
                )
                block_size = len(outer_test)
                blocks = [
                    np.arange(
                        index * block_size,
                        (index + 1) * block_size,
                        dtype=np.int64,
                    )
                    for index in range(len(synthetic_evaluation_draws))
                ]
                synthetic_permutation, _ = _permutation_rows(
                    classifier,
                    synthetic_evaluation,
                    synthetic_labels,
                    genes,
                    data.symbols,
                    permutation_repeats=permutation_repeats,
                    seed=run_seed + tissue_position * 1_000_000 + repeat * 10_000 + arm_offsets[arm] * 100 + 1,
                    blocks=blocks,
                )
                synthetic_shap, _ = _linear_shap_rows(
                    classifier,
                    synthetic_evaluation,
                    synthetic_labels,
                    background,
                    genes,
                    data.symbols,
                )
                synthetic_table = synthetic_permutation.merge(
                    synthetic_shap,
                    on=["gene", "symbol"],
                    validate="one_to_one",
                )
                synthetic_table["domain"] = "synthetic"
                table = pd.concat((real_table, synthetic_table), ignore_index=True)
                table.insert(0, "metric_reconstruction_max_error", metric_error)
                table.insert(0, "coefficient_reconstruction_max_error", coefficient_error)
                table.insert(0, "classifier_coefficient", np.tile(classifier.coef_[0], 2))
                table.insert(0, "regularization_c", float(candidate["regularization_c"]))
                table.insert(0, "feature_count", int(candidate["feature_count"]))
                table.insert(0, "rank_method", str(candidate["rank_method"]))
                table.insert(0, "arm", arm)
                table.insert(0, "repeat", repeat)
                table.insert(0, "tissue", tissue)
                table.insert(0, "scope", scope)
                result_tables.append(table)
    repeat_table = pd.concat(result_tables, ignore_index=True)
    repeats_table = pd.DataFrame(
        [
            {"scope": scope, "tissue": tissue, "completed_repeats": repeats}
            for tissue, repeats in completed_repeats.items()
        ]
    )
    minimum = float(workflow["stability"]["minimum_selection_frequency"])
    minimum_sign = float(
        workflow["stability"]["minimum_coefficient_sign_agreement"]
    )
    thresholds = {
        tissue: (minimum, minimum_sign) for tissue in completed_repeats
    }
    return repeat_table, thresholds, data.genes, data.symbols, Path(
        workflow["annotations"]["reactome_gmt"]
    )


def run(
    config_path: Path,
    *,
    scopes_override: set[str] | None = None,
    units_override: set[str] | None = None,
    permutation_repeats_override: int | None = None,
    output_override: Path | None = None,
) -> Path:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    output = output_override or Path(config["run"]["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    permutation_repeats = int(
        permutation_repeats_override
        or config["analysis"]["permutation_repeats"]
    )
    run_seed = int(config["run"]["seed"])
    coefficient_tolerance = float(
        config["analysis"].get("coefficient_tolerance", 0.05)
    )
    metric_tolerance = float(
        config["analysis"].get("metric_tolerance", 1e-10)
    )
    repeat_tables: list[pd.DataFrame] = []
    choice_tables: list[pd.DataFrame] = []
    threshold_lookup: dict[tuple[str, str], tuple[float, float]] = {}
    completed_rows: list[dict[str, object]] = []
    all_genes: list[str] | None = None
    all_symbols: dict[str, str] = {}
    reactome_path: Path | None = None
    for source in config["sources"]:
        scope = str(source["scope"])
        if scopes_override is not None and scope not in scopes_override:
            continue
        choices = pd.read_csv(
            Path(source["analysis_dir"]) / "tissue_arm_choices.tsv", sep="\t"
        )
        choices.insert(0, "scope", scope)
        choice_tables.append(choices)
        table, thresholds, genes, symbols, gmt_path = _run_source(
            source,
            permutation_repeats=permutation_repeats,
            run_seed=run_seed,
            coefficient_tolerance=coefficient_tolerance,
            metric_tolerance=metric_tolerance,
            units_override=units_override,
        )
        repeat_tables.append(table)
        for tissue, threshold in thresholds.items():
            threshold_lookup[(scope, tissue)] = threshold
            completed_rows.append(
                {
                    "scope": scope,
                    "tissue": tissue,
                    "completed_repeats": int(
                        table.loc[table["tissue"].eq(tissue), "repeat"].nunique()
                    ),
                }
            )
        if all_genes is None:
            all_genes = genes
        elif all_genes != genes:
            raise ValueError("Importance sources use different gene panels")
        all_symbols.update(symbols)
        if reactome_path is None:
            reactome_path = gmt_path
        elif reactome_path != gmt_path:
            raise ValueError("Importance sources use different Reactome files")
    if not repeat_tables:
        raise RuntimeError("No classifier-importance source was selected")
    repeat_table = pd.concat(repeat_tables, ignore_index=True)
    completed = pd.DataFrame(completed_rows).drop_duplicates(["scope", "tissue"])
    aggregate = _aggregate_importance(repeat_table, completed)
    comparison = _compare_arms(aggregate, threshold_lookup)
    similarity = _arm_similarity(aggregate, comparison)
    choices = pd.concat(choice_tables, ignore_index=True)
    selected_importance, selected_comparison = _selected_arm_tables(
        aggregate, comparison, choices
    )
    inventory_path_value = config["analysis"].get("synthetic_informed_genes")
    synthetic_informed = (
        _synthetic_informed_gene_importance(
            Path(inventory_path_value), selected_importance, aggregate
        )
        if inventory_path_value
        else pd.DataFrame()
    )
    if reactome_path is None:
        raise RuntimeError("No Reactome annotation path was loaded")
    enrichment = _importance_enrichment(
        comparison,
        background=all_genes or [],
        gmt_path=reactome_path,
        symbols=all_symbols,
    )
    top = (
        aggregate.loc[aggregate["domain"].eq("real")]
        .sort_values(
            ["scope", "tissue", "arm", "permutation_roc_auc_mean"],
            ascending=[True, True, True, False],
            kind="stable",
        )
        .groupby(["scope", "tissue", "arm"], observed=True)
        .head(int(config["analysis"].get("top_genes_per_arm", 20)))
    )
    repeat_table.to_csv(
        output / "importance_by_repeat.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    aggregate.to_csv(
        output / "importance_summary.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    comparison.to_csv(
        output / "arm_vs_real_gene_comparison.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    selected_importance.to_csv(
        output / "selected_arm_importance.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    selected_comparison.to_csv(
        output / "selected_arm_vs_real_gene_comparison.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    if not synthetic_informed.empty:
        synthetic_informed.to_csv(
            output / "synthetic_informed_bh_fdr_gene_importance.tsv",
            sep="\t",
            index=False,
        )
    similarity.to_csv(output / "tissue_arm_similarity.tsv", sep="\t", index=False)
    enrichment.to_csv(
        output / "reactome_importance_enrichment.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    top.to_csv(output / "top_importance_genes.tsv", sep="\t", index=False)
    for row in completed.itertuples(index=False):
        _plot_unit_importance(
            aggregate, output, scope=str(row.scope), tissue=str(row.tissue)
        )
    _plot_domain_similarity(similarity, output)
    pattern_counts = comparison["pattern"].value_counts().to_dict()
    selected_pattern_counts = selected_comparison["pattern"].value_counts().to_dict()
    arm_medians = (
        similarity.groupby("arm", observed=True)[
            [
                "real_vs_arm_importance_spearman",
                "arm_real_vs_synthetic_domain_spearman",
                "top20_overlap",
            ]
        ]
        .median()
        .to_dict(orient="index")
    )
    summary = {
        "status": "complete",
        "config": str(config_path.resolve()),
        "output": str(output.resolve()),
        "completed_units": int(len(completed)),
        "scopes": completed["scope"].value_counts().to_dict(),
        "permutation_repeats": permutation_repeats,
        "fitted_classifiers": int(
            repeat_table[["scope", "tissue", "repeat", "arm"]]
            .drop_duplicates()
            .shape[0]
        ),
        "importance_rows": int(len(repeat_table)),
        "pattern_counts": {key: int(value) for key, value in pattern_counts.items()},
        "selected_arm_pattern_counts": {
            key: int(value) for key, value in selected_pattern_counts.items()
        },
        "synthetic_informed_bh_fdr_genes": int(len(synthetic_informed)),
        "synthetic_informed_positive_permutation": int(
            synthetic_informed.get(
                "positive_held_out_real_permutation_importance", pd.Series(dtype=bool)
            ).sum()
        ),
        "synthetic_informed_positive_shap_separation": int(
            synthetic_informed.get(
                "positive_linear_shap_flt_gc_separation", pd.Series(dtype=bool)
            ).sum()
        ),
        "median_similarity_by_arm": arm_medians,
        "linear_shap_method": "exact interventional linear SHAP on log-odds",
        "neural_network_retrained": False,
        "classifier_hyperparameters_reselected": False,
        "maximum_coefficient_reconstruction_error": float(
            repeat_table["coefficient_reconstruction_max_error"].max()
        ),
        "coefficient_reconstruction_tolerance": coefficient_tolerance,
        "maximum_metric_reconstruction_error": float(
            repeat_table["metric_reconstruction_max_error"].max()
        ),
        "metric_reconstruction_tolerance": metric_tolerance,
        "limitations": [
            "The fixed DDIM was trained before the repeated classifier splits.",
            "Outer splits retain represented accessions and measure within-study interpolation.",
            "Repeated outer splits overlap and are not independent observations.",
            "Marginal permutation importance can be diluted across correlated genes.",
            "Synthetic-domain importance is model-implied evidence, not an independent biological replicate.",
        ],
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    _write_readme(output, summary)
    print(json.dumps(summary, indent=2), flush=True)
    return output / "summary.json"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--scopes", nargs="+", choices=("tissue", "muscle_group"))
    parser.add_argument("--units", nargs="+")
    parser.add_argument("--permutation-repeats", type=int)
    parser.add_argument("--output", type=Path)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    run(
        arguments.config,
        scopes_override=set(arguments.scopes) if arguments.scopes else None,
        units_override=set(arguments.units) if arguments.units else None,
        permutation_repeats_override=arguments.permutation_repeats,
        output_override=arguments.output,
    )


if __name__ == "__main__":
    main()
