"""Matched all-gene real, generated, and augmented FLT/GC classifiers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from .classifier_importance import (
    ARM_LABELS,
    _aggregate_importance,
    _linear_shap_rows,
    _permutation_rows,
    _safe_spearman,
)
from .generated_feature_guidance import _reactome_enrichment
from .within_study_feature_stability import (
    METRICS,
    WorkflowData,
    _labels,
    _load_data,
    _metric_set,
    _muscle_group_analysis_data,
    _valid_nested_split,
)


ARMS = ("real_only", "generated_only", "real_plus_generated")
SYNTHETIC_ARMS = ARMS[1:]
MATCHED_ARM_LABELS = {arm: ARM_LABELS[arm] for arm in ARMS}


def _fit_ridge(
    expression: np.ndarray,
    labels: np.ndarray,
    *,
    regularization_c: float,
    seed: int,
    sample_weight: np.ndarray | None = None,
) -> LogisticRegression:
    classifier = LogisticRegression(
        C=float(regularization_c),
        class_weight="balanced",
        max_iter=5000,
        random_state=int(seed),
        solver="lbfgs",
    )
    classifier.fit(expression, labels, sample_weight=sample_weight)
    return classifier


def _fit_matched_arm(
    arm: str,
    *,
    real_scaled: np.ndarray,
    labels: np.ndarray,
    synthetic_scaled: list[np.ndarray],
    regularization_c: float,
    synthetic_weight: float,
    seed: int,
) -> LogisticRegression:
    if arm not in ARMS:
        raise ValueError(f"Unknown matched classifier arm: {arm}")
    if arm == "real_only":
        return _fit_ridge(
            real_scaled,
            labels,
            regularization_c=regularization_c,
            seed=seed,
        )
    synthetic = np.concatenate(synthetic_scaled)
    synthetic_labels = np.tile(labels, len(synthetic_scaled))
    if arm == "generated_only":
        return _fit_ridge(
            synthetic,
            synthetic_labels,
            regularization_c=regularization_c,
            seed=seed,
        )
    per_row_weight = float(synthetic_weight) / len(synthetic_scaled)
    expression = np.concatenate((real_scaled, synthetic))
    combined_labels = np.concatenate((labels, synthetic_labels))
    sample_weight = np.concatenate(
        (
            np.ones(len(real_scaled), dtype=float),
            np.full(len(synthetic), per_row_weight, dtype=float),
        )
    )
    return _fit_ridge(
        expression,
        combined_labels,
        regularization_c=regularization_c,
        seed=seed,
        sample_weight=sample_weight,
    )


def _accession_blocks(samples: pd.DataFrame, *, offset: int = 0) -> list[np.ndarray]:
    blocks: list[np.ndarray] = []
    for _, indices in samples.groupby("accession", sort=True, observed=True).groups.items():
        blocks.append(np.asarray(list(indices), dtype=np.int64) + int(offset))
    if not blocks:
        raise ValueError("No accession blocks were available")
    return blocks


def _synthetic_accession_blocks(
    samples: pd.DataFrame, number_of_draws: int
) -> list[np.ndarray]:
    block_size = len(samples)
    blocks: list[np.ndarray] = []
    for draw in range(int(number_of_draws)):
        blocks.extend(_accession_blocks(samples, offset=draw * block_size))
    return blocks


def _accession_macro_metrics(
    samples: pd.DataFrame, labels: np.ndarray, probabilities: np.ndarray
) -> dict[str, float]:
    rows: list[dict[str, float]] = []
    for _, indices in samples.groupby("accession", sort=True, observed=True).groups.items():
        positions = np.asarray(list(indices), dtype=np.int64)
        subset_labels = labels[positions]
        if np.unique(subset_labels).size != 2:
            continue
        rows.append(_metric_set(subset_labels, probabilities[positions]))
    if not rows:
        return {metric: float("nan") for metric in METRICS}
    return {
        metric: float(np.mean([row[metric] for row in rows])) for metric in METRICS
    }


def _evaluate_classifier(
    classifier: LogisticRegression,
    expression: np.ndarray,
    labels: np.ndarray,
    samples: pd.DataFrame,
) -> dict[str, float]:
    probability = classifier.predict_proba(expression)[:, 1]
    pooled = _metric_set(labels, probability)
    macro = _accession_macro_metrics(samples, labels, probability)
    return {
        **pooled,
        **{f"accession_macro_{metric}": value for metric, value in macro.items()},
    }


def _select_shared_regularization(
    real: np.ndarray,
    labels: np.ndarray,
    samples: pd.DataFrame,
    *,
    inner_train: np.ndarray,
    inner_validation: np.ndarray,
    regularization_grid: Iterable[float],
    seed: int,
) -> tuple[float, pd.DataFrame]:
    """Select one ridge penalty using real data and apply it to every arm."""

    scaler = StandardScaler().fit(real[inner_train])
    train = scaler.transform(real[inner_train])
    validation = scaler.transform(real[inner_validation])
    rows: list[dict[str, float]] = []
    for index, regularization_c in enumerate(regularization_grid):
        classifier = _fit_ridge(
            train,
            labels[inner_train],
            regularization_c=float(regularization_c),
            seed=int(seed) + index,
        )
        metrics = _evaluate_classifier(
            classifier,
            validation,
            labels[inner_validation],
            samples.loc[inner_validation].reset_index(drop=True),
        )
        rows.append({"regularization_c": float(regularization_c), **metrics})
    table = pd.DataFrame(rows)
    selected = table.sort_values(
        [
            "balanced_accuracy",
            "roc_auc",
            "average_precision",
            "accession_macro_balanced_accuracy",
            "accession_macro_roc_auc",
            "accession_macro_average_precision",
            "regularization_c",
        ],
        ascending=[False, False, False, False, False, False, True],
        kind="stable",
    ).iloc[0]
    return float(selected["regularization_c"]), table


def _run_unit(
    *,
    scope: str,
    tissue: str,
    data: WorkflowData,
    repeats: int,
    tissue_seed: int,
    outer_fraction: float,
    inner_fraction: float,
    regularization_grid: list[float],
    synthetic_weight: float,
    permutation_repeats: int,
    permutation_seed: int,
    metrics_only: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mask = data.development_samples["tissue"].astype(str).eq(tissue)
    positions = np.flatnonzero(mask.to_numpy())
    real = data.development_expression[positions]
    samples = data.development_samples.loc[mask].reset_index(drop=True)
    labels = _labels(samples)
    synthetic = [draw[positions] for draw in data.synthetic_draws.values()]
    metric_rows: list[dict[str, object]] = []
    candidate_tables: list[pd.DataFrame] = []
    importance_tables: list[pd.DataFrame] = []
    for repeat in range(int(repeats)):
        split = _valid_nested_split(
            samples,
            outer_fraction=outer_fraction,
            inner_fraction=inner_fraction,
            seed=tissue_seed + repeat,
        )
        if split is None:
            continue
        inner_train, inner_validation, outer_test = split
        outer_train = np.setdiff1d(np.arange(len(samples)), outer_test)
        selected_c, candidates = _select_shared_regularization(
            real,
            labels,
            samples,
            inner_train=inner_train,
            inner_validation=inner_validation,
            regularization_grid=regularization_grid,
            seed=tissue_seed + repeat * 100,
        )
        candidates.insert(0, "repeat", repeat)
        candidates.insert(0, "tissue", tissue)
        candidates.insert(0, "scope", scope)
        candidate_tables.append(candidates)

        scaler = StandardScaler().fit(real[outer_train])
        real_train = scaler.transform(real[outer_train])
        real_test = scaler.transform(real[outer_test])
        synthetic_train = [scaler.transform(draw[outer_train]) for draw in synthetic]
        synthetic_test_draws = [scaler.transform(draw[outer_test]) for draw in synthetic]
        synthetic_test = np.concatenate(synthetic_test_draws)
        synthetic_test_labels = np.tile(labels[outer_test], len(synthetic_test_draws))
        real_test_samples = samples.loc[outer_test].reset_index(drop=True)
        shared_background = real_train.mean(axis=0)
        for arm_offset, arm in enumerate(ARMS):
            classifier = _fit_matched_arm(
                arm,
                real_scaled=real_train,
                labels=labels[outer_train],
                synthetic_scaled=synthetic_train,
                regularization_c=selected_c,
                synthetic_weight=synthetic_weight,
                seed=tissue_seed + repeat * 100 + arm_offset,
            )
            real_metrics = _evaluate_classifier(
                classifier,
                real_test,
                labels[outer_test],
                real_test_samples,
            )
            metric_rows.append(
                {
                    "scope": scope,
                    "tissue": tissue,
                    "repeat": repeat,
                    "arm": arm,
                    "genes": int(real.shape[1]),
                    "regularization_c": selected_c,
                    "outer_train_profiles": int(len(outer_train)),
                    "outer_test_profiles": int(len(outer_test)),
                    **real_metrics,
                }
            )
            if metrics_only:
                continue
            genes = data.genes
            real_permutation, _ = _permutation_rows(
                classifier,
                real_test,
                labels[outer_test],
                genes,
                data.symbols,
                permutation_repeats=permutation_repeats,
                seed=(
                    permutation_seed
                    + tissue_seed
                    + repeat * 10_000
                    + arm_offset * 100
                ),
                blocks=_accession_blocks(real_test_samples),
            )
            real_shap, _ = _linear_shap_rows(
                classifier,
                real_test,
                labels[outer_test],
                shared_background,
                genes,
                data.symbols,
            )
            real_table = real_permutation.merge(
                real_shap, on=["gene", "symbol"], validate="one_to_one"
            )
            real_table["domain"] = "real"

            synthetic_permutation, _ = _permutation_rows(
                classifier,
                synthetic_test,
                synthetic_test_labels,
                genes,
                data.symbols,
                permutation_repeats=permutation_repeats,
                seed=(
                    permutation_seed
                    + tissue_seed
                    + repeat * 10_000
                    + arm_offset * 100
                    + 1
                ),
                blocks=_synthetic_accession_blocks(
                    real_test_samples, len(synthetic_test_draws)
                ),
            )
            synthetic_shap, _ = _linear_shap_rows(
                classifier,
                synthetic_test,
                synthetic_test_labels,
                shared_background,
                genes,
                data.symbols,
            )
            synthetic_table = synthetic_permutation.merge(
                synthetic_shap, on=["gene", "symbol"], validate="one_to_one"
            )
            synthetic_table["domain"] = "synthetic"
            table = pd.concat((real_table, synthetic_table), ignore_index=True)
            table.insert(
                0,
                "classifier_coefficient",
                np.tile(classifier.coef_[0], 2),
            )
            table.insert(0, "regularization_c", selected_c)
            table.insert(0, "arm", arm)
            table.insert(0, "repeat", repeat)
            table.insert(0, "tissue", tissue)
            table.insert(0, "scope", scope)
            importance_tables.append(table)
    if not metric_rows:
        raise RuntimeError(f"No valid matched splits for {scope}:{tissue}")
    return (
        pd.DataFrame(metric_rows),
        pd.concat(candidate_tables, ignore_index=True),
        (
            pd.concat(importance_tables, ignore_index=True)
            if importance_tables
            else pd.DataFrame()
        ),
    )


def _metric_summary(metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_columns = [
        *METRICS,
        *[f"accession_macro_{metric}" for metric in METRICS],
    ]
    summary_spec: dict[str, tuple[str, str]] = {
        "repeats": ("repeat", "nunique"),
    }
    for metric in metric_columns:
        summary_spec[f"mean_{metric}"] = (metric, "mean")
        summary_spec[f"sd_{metric}"] = (metric, "std")
    summary = (
        metrics.groupby(["scope", "tissue", "arm"], observed=True)
        .agg(**summary_spec)
        .reset_index()
    )
    rows: list[dict[str, object]] = []
    for (scope, tissue), frame in metrics.groupby(
        ["scope", "tissue"], sort=True, observed=True
    ):
        baseline = frame.loc[frame["arm"].eq("real_only")]
        for arm in SYNTHETIC_ARMS:
            candidate = frame.loc[frame["arm"].eq(arm)]
            paired = baseline.merge(
                candidate,
                on="repeat",
                how="inner",
                suffixes=("_real", "_arm"),
                validate="one_to_one",
            )
            record: dict[str, object] = {
                "scope": scope,
                "tissue": tissue,
                "arm": arm,
                "paired_repeats": int(len(paired)),
            }
            nonworse = []
            strict = []
            for metric in metric_columns:
                delta = paired[f"{metric}_arm"] - paired[f"{metric}_real"]
                record[f"real_mean_{metric}"] = float(
                    paired[f"{metric}_real"].mean()
                )
                record[f"arm_mean_{metric}"] = float(
                    paired[f"{metric}_arm"].mean()
                )
                record[f"mean_delta_{metric}"] = float(delta.mean())
                record[f"nonworse_rate_{metric}"] = float((delta >= 0.0).mean())
                record[f"strict_win_rate_{metric}"] = float((delta > 0.0).mean())
                nonworse.append(delta.to_numpy(dtype=float) >= 0.0)
                strict.append(delta.to_numpy(dtype=float) > 0.0)
            pooled_means_nonworse = all(
                float(record[f"mean_delta_{metric}"]) >= -1e-12 for metric in METRICS
            )
            macro_means_nonworse = all(
                float(record[f"mean_delta_accession_macro_{metric}"]) >= -1e-12
                for metric in METRICS
            )
            record["pooled_mean_all_metrics_nonworse"] = pooled_means_nonworse
            record["macro_mean_all_metrics_nonworse"] = macro_means_nonworse
            record["joint_mean_all_metrics_nonworse"] = bool(
                pooled_means_nonworse and macro_means_nonworse
            )
            pooled_matrix = np.column_stack(nonworse[: len(METRICS)])
            pooled_strict = np.column_stack(strict[: len(METRICS)])
            record["pooled_all_metrics_nonworse_rate"] = float(
                np.mean(np.all(pooled_matrix, axis=1))
            )
            record["pooled_all_metrics_strict_win_rate"] = float(
                np.mean(np.all(pooled_strict, axis=1))
            )
            rows.append(record)
    return summary, pd.DataFrame(rows)


def _positive_importance(
    frame: pd.DataFrame,
    prefix: str,
    *,
    minimum_importance: float,
    minimum_fraction: float,
) -> pd.Series:
    return frame[f"{prefix}_permutation_roc_auc_mean"].ge(minimum_importance) & frame[
        f"{prefix}_permutation_roc_auc_positive_repeat_fraction"
    ].ge(minimum_fraction)


def _importance_pattern(frame: pd.DataFrame) -> pd.Series:
    real_positive = frame["real_only_positive_importance"]
    arm_positive = frame["arm_real_positive_importance"]
    synthetic_positive = frame["arm_synthetic_positive_importance"]
    direction_match = frame["coefficient_direction_match"]
    reinforced = (
        arm_positive
        & real_positive
        & direction_match
        & frame["arm_real_permutation_roc_auc_mean"].ge(
            frame["real_only_permutation_roc_auc_mean"]
        )
    )
    attenuated = arm_positive & real_positive & direction_match & ~reinforced
    return pd.Series(
        np.select(
            [
                arm_positive & ~real_positive,
                arm_positive & real_positive & ~direction_match,
                reinforced,
                attenuated,
                ~arm_positive & real_positive,
                ~arm_positive & ~real_positive & synthetic_positive,
            ],
            [
                "synthetic_promoted_real_transfer",
                "shared_direction_conflict",
                "shared_reinforced",
                "shared_attenuated",
                "real_only",
                "synthetic_domain_only",
            ],
            default="no_marginal_importance",
        ),
        index=frame.index,
    )


def _importance_comparison(
    aggregate: pd.DataFrame,
    *,
    minimum_importance: float,
    minimum_fraction: float,
) -> pd.DataFrame:
    real_domain = aggregate.loc[aggregate["domain"].eq("real")]
    synthetic_domain = aggregate.loc[aggregate["domain"].eq("synthetic")]
    value_columns = [
        "selection_frequency",
        "median_classifier_coefficient",
        "coefficient_sign_agreement",
        "permutation_balanced_accuracy_mean",
        "permutation_balanced_accuracy_positive_repeat_fraction",
        "permutation_roc_auc_mean",
        "permutation_roc_auc_positive_repeat_fraction",
        "permutation_average_precision_mean",
        "permutation_average_precision_positive_repeat_fraction",
        "linear_shap_mean_absolute",
        "linear_shap_flight_minus_ground",
    ]
    keys = ["scope", "tissue", "gene", "symbol"]
    baseline = real_domain.loc[real_domain["arm"].eq("real_only"), keys + value_columns]
    baseline = baseline.rename(
        columns={column: f"real_only_{column}" for column in value_columns}
    )
    tables: list[pd.DataFrame] = []
    for arm in SYNTHETIC_ARMS:
        arm_real = real_domain.loc[
            real_domain["arm"].eq(arm), keys + value_columns
        ].rename(columns={column: f"arm_real_{column}" for column in value_columns})
        arm_synthetic = synthetic_domain.loc[
            synthetic_domain["arm"].eq(arm), keys + value_columns
        ].rename(
            columns={column: f"arm_synthetic_{column}" for column in value_columns}
        )
        table = baseline.merge(
            arm_real, on=keys, how="inner", validate="one_to_one"
        ).merge(arm_synthetic, on=keys, how="inner", validate="one_to_one")
        table.insert(2, "arm", arm)
        table["real_only_positive_importance"] = _positive_importance(
            table,
            "real_only",
            minimum_importance=minimum_importance,
            minimum_fraction=minimum_fraction,
        )
        table["arm_real_positive_importance"] = _positive_importance(
            table,
            "arm_real",
            minimum_importance=minimum_importance,
            minimum_fraction=minimum_fraction,
        )
        table["arm_synthetic_positive_importance"] = _positive_importance(
            table,
            "arm_synthetic",
            minimum_importance=minimum_importance,
            minimum_fraction=minimum_fraction,
        )
        table["coefficient_direction_match"] = np.sign(
            table["real_only_median_classifier_coefficient"]
        ).eq(np.sign(table["arm_real_median_classifier_coefficient"]))
        table["arm_minus_real_permutation_roc_auc"] = (
            table["arm_real_permutation_roc_auc_mean"]
            - table["real_only_permutation_roc_auc_mean"]
        )
        table["pattern"] = _importance_pattern(table)
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


def _importance_similarity(comparison: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (scope, tissue, arm), frame in comparison.groupby(
        ["scope", "tissue", "arm"], sort=True, observed=True
    ):
        real_top = set(
            frame.nlargest(20, "real_only_permutation_roc_auc_mean")["gene"]
        )
        arm_top = set(frame.nlargest(20, "arm_real_permutation_roc_auc_mean")["gene"])
        rows.append(
            {
                "scope": scope,
                "tissue": tissue,
                "arm": arm,
                "real_vs_arm_importance_spearman": _safe_spearman(
                    frame["real_only_permutation_roc_auc_mean"],
                    frame["arm_real_permutation_roc_auc_mean"],
                ),
                "arm_real_vs_synthetic_domain_spearman": _safe_spearman(
                    frame["arm_real_permutation_roc_auc_mean"],
                    frame["arm_synthetic_permutation_roc_auc_mean"],
                ),
                "top20_overlap": int(len(real_top & arm_top)),
                **{
                    f"pattern_{pattern}": int((frame["pattern"] == pattern).sum())
                    for pattern in sorted(frame["pattern"].unique())
                },
            }
        )
    return pd.DataFrame(rows)


def _bh_fdr_crosswalk(
    inventory_path: Path,
    comparison: pd.DataFrame,
    utility: pd.DataFrame,
) -> pd.DataFrame:
    inventory = pd.read_csv(inventory_path, sep="\t")
    inventory["scope"] = inventory["analysis_scope"].map(
        {
            "canonical_tissue": "tissue",
            "skeletal_muscle_group": "muscle_group",
        }
    )
    if inventory["scope"].isna().any():
        raise ValueError("BH-FDR inventory has unknown analysis scopes")
    importance = comparison.rename(columns={"symbol": "importance_symbol"})
    table = inventory.merge(
        importance,
        on=["scope", "tissue", "gene"],
        how="inner",
        validate="one_to_many",
    ).merge(
        utility[
            [
                "scope",
                "tissue",
                "arm",
                "pooled_mean_all_metrics_nonworse",
                "macro_mean_all_metrics_nonworse",
                "joint_mean_all_metrics_nonworse",
            ]
        ],
        on=["scope", "tissue", "arm"],
        how="left",
        validate="many_to_one",
    )
    table["symbol_matches_importance_annotation"] = table["symbol"].eq(
        table["importance_symbol"]
    )
    table["arm_coefficient_matches_real_effect"] = np.sign(
        table["arm_real_median_classifier_coefficient"]
    ).eq(np.sign(table["meta_effect"]))
    table["arm_shap_supports_flt_gc"] = table[
        "arm_real_linear_shap_flight_minus_ground"
    ].gt(0.0)
    table["matched_importance_interpretation"] = np.select(
        [
            table["pattern"].eq("synthetic_promoted_real_transfer"),
            table["pattern"].isin(["shared_reinforced", "shared_attenuated"]),
            table["pattern"].eq("real_only"),
            table["pattern"].eq("synthetic_domain_only"),
            table["pattern"].eq("shared_direction_conflict"),
        ],
        [
            "synthetic_promoted_importance",
            "reinforced_importance",
            "real_only_importance",
            "synthetic_domain_only",
            "direction_conflict",
        ],
        default="no_repeat_consistent_marginal_importance",
    )
    table["eligible_synthetic_biological_candidate"] = (
        table["joint_mean_all_metrics_nonworse"].fillna(False)
        & table["matched_importance_interpretation"].isin(
            ["synthetic_promoted_importance", "reinforced_importance"]
        )
        & table["arm_coefficient_matches_real_effect"]
        & table["arm_shap_supports_flt_gc"]
    )
    return table


def _candidate_reactome_enrichment(
    eligible: pd.DataFrame,
    *,
    background: list[str],
    gmt_path: Path,
    symbols: dict[str, str],
) -> pd.DataFrame:
    unique = eligible.drop_duplicates(["scope", "tissue", "gene"])
    tables: list[pd.DataFrame] = []
    for (scope, tissue), frame in unique.groupby(
        ["scope", "tissue"], sort=True, observed=True
    ):
        gene_sets = {
            "all": frame,
            "flt_higher": frame.loc[frame["flt_gc_direction"].eq("FLT_higher")],
            "flt_lower": frame.loc[frame["flt_gc_direction"].eq("FLT_lower")],
        }
        for gene_set, subset in gene_sets.items():
            if len(subset) < 2:
                continue
            table = _reactome_enrichment(
                subset["gene"].astype(str).tolist(),
                background,
                gmt_path,
                symbols,
            )
            if table.empty:
                continue
            table.insert(0, "gene_set", gene_set)
            table.insert(0, "tissue", tissue)
            table.insert(0, "scope", scope)
            tables.append(table)
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def _plot_utility(utility: pd.DataFrame, output: Path) -> None:
    units = sorted(
        utility.assign(unit=utility["scope"] + ":" + utility["tissue"])["unit"].unique()
    )
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(13.5, max(7.0, 0.31 * len(units))),
        sharey=True,
    )
    for axis, arm in zip(axes, SYNTHETIC_ARMS):
        frame = utility.loc[utility["arm"].eq(arm)].assign(
            unit=lambda values: values["scope"] + ":" + values["tissue"]
        )
        matrix = (
            frame.set_index("unit")
            .reindex(units)[[f"mean_delta_{metric}" for metric in METRICS]]
            .to_numpy(dtype=float)
        )
        limit = max(float(np.nanmax(np.abs(matrix))), 0.05)
        image = axis.imshow(
            matrix,
            aspect="auto",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
        )
        axis.set_xticks(
            np.arange(len(METRICS)),
            ["Balanced accuracy", "AUROC", "Average precision"],
            rotation=30,
            ha="right",
        )
        axis.set_yticks(np.arange(len(units)))
        if axis is axes[0]:
            axis.set_yticklabels(units, fontsize=7)
        else:
            axis.tick_params(axis="y", labelleft=False, left=False)
        axis.set_title(MATCHED_ARM_LABELS[arm], weight="bold")
        colorbar = figure.colorbar(image, ax=axis, fraction=0.035, pad=0.02)
        colorbar.set_label("Mean change")
    axes[0].set_ylabel("Analysis unit")
    figure.suptitle(
        "Matched all-gene classifier change from real-only",
        fontsize=14,
        weight="bold",
    )
    figure.subplots_adjust(
        left=0.18,
        right=0.92,
        top=0.92,
        bottom=0.15,
        wspace=0.20,
    )
    figure.savefig(output / "matched_classifier_metric_deltas.png", dpi=220)
    figure.savefig(output / "matched_classifier_metric_deltas.pdf")
    plt.close(figure)


def _plot_unit_importance(
    aggregate: pd.DataFrame,
    comparison: pd.DataFrame,
    output: Path,
    *,
    scope: str,
    tissue: str,
    minimum_importance: float,
) -> None:
    subset = aggregate.loc[
        aggregate["scope"].eq(scope)
        & aggregate["tissue"].eq(tissue)
        & aggregate["domain"].eq("real")
    ].copy()
    compare = comparison.loc[
        comparison["scope"].eq(scope) & comparison["tissue"].eq(tissue)
    ].copy()
    if subset.empty or compare.empty:
        return
    top = (
        subset.groupby(["gene", "symbol"], observed=True)[
            "permutation_roc_auc_mean"
        ]
        .max()
        .nlargest(25)
    )
    top_genes = [gene for gene, _ in top.index]
    top_labels = [symbol or gene for gene, symbol in top.index]
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 11.0))
    pattern_colors = {
        "synthetic_promoted_real_transfer": "#E76F51",
        "shared_reinforced": "#2A9D8F",
        "shared_attenuated": "#75B8AD",
        "real_only": "#457B9D",
        "synthetic_domain_only": "#E9C46A",
        "shared_direction_conflict": "#9B5DE5",
        "no_marginal_importance": "#C8C8C8",
    }
    pattern_labels = {
        "synthetic_promoted_real_transfer": "Synthetic-promoted transfer",
        "shared_reinforced": "Shared, reinforced",
        "shared_attenuated": "Shared, attenuated",
        "real_only": "Real-only importance",
        "synthetic_domain_only": "Synthetic-domain only",
        "shared_direction_conflict": "Direction conflict",
        "no_marginal_importance": "Below importance gate",
    }
    for axis, arm in zip(axes[0], SYNTHETIC_ARMS):
        frame = compare.loc[compare["arm"].eq(arm)]
        axis.scatter(
            frame["real_only_permutation_roc_auc_mean"],
            frame["arm_real_permutation_roc_auc_mean"],
            s=10,
            c=[pattern_colors[str(pattern)] for pattern in frame["pattern"]],
            alpha=0.65,
            linewidths=0,
        )
        limits = np.asarray(
            [
                axis.get_xlim()[0],
                axis.get_xlim()[1],
                axis.get_ylim()[0],
                axis.get_ylim()[1],
            ]
        )
        low, high = float(limits.min()), float(limits.max())
        axis.plot([low, high], [low, high], color="#333333", linewidth=0.8)
        axis.axhline(minimum_importance, color="#777777", linewidth=0.7, linestyle=":")
        axis.axvline(minimum_importance, color="#777777", linewidth=0.7, linestyle=":")
        axis.set_xlim(low, high)
        axis.set_ylim(low, high)
        axis.set_xlabel("Real-only AUROC permutation loss")
        axis.set_ylabel(f"{MATCHED_ARM_LABELS[arm]} loss")
        axis.set_title(MATCHED_ARM_LABELS[arm], weight="bold")
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgewidth=0,
            markersize=6,
            label=pattern_labels[pattern],
        )
        for pattern, color in pattern_colors.items()
    ]
    figure.legend(
        handles=legend_handles,
        loc="center",
        bbox_to_anchor=(0.53, 0.505),
        ncol=4,
        frameon=False,
        fontsize=8,
    )
    panels = [
        ("permutation_roc_auc_mean", "Held-out real AUROC permutation loss"),
        ("linear_shap_flight_minus_ground", "Linear SHAP FLT-GC contribution"),
    ]
    for axis, (value, title) in zip(axes[1], panels):
        pivot = subset.pivot_table(index="gene", columns="arm", values=value)
        matrix = pivot.reindex(index=top_genes, columns=ARMS).fillna(0.0)
        limit = max(float(np.nanmax(np.abs(matrix.to_numpy()))), 1e-6)
        image = axis.imshow(
            matrix,
            aspect="auto",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
        )
        axis.set_xticks(
            np.arange(len(ARMS)),
            [MATCHED_ARM_LABELS[arm] for arm in ARMS],
            rotation=30,
            ha="right",
        )
        axis.set_yticks(np.arange(len(top_genes)), top_labels, fontsize=8)
        axis.set_title(title, weight="bold")
        figure.colorbar(image, ax=axis, fraction=0.035, pad=0.02)
    figure.suptitle(
        f"{tissue.replace('_', ' ').title()}: matched all-gene classifiers",
        fontsize=15,
        weight="bold",
    )
    figure.subplots_adjust(
        left=0.10,
        right=0.96,
        top=0.93,
        bottom=0.10,
        hspace=0.52,
        wspace=0.28,
    )
    directory = output / scope / tissue
    directory.mkdir(parents=True, exist_ok=True)
    figure.savefig(directory / "matched_classifier_importance.png", dpi=220)
    figure.savefig(directory / "matched_classifier_importance.pdf")
    plt.close(figure)


def _plot_bh_candidate_counts(crosswalk: pd.DataFrame, output: Path) -> None:
    if crosswalk.empty:
        return
    supported = crosswalk.loc[crosswalk["eligible_synthetic_biological_candidate"]]
    counts = (
        supported.groupby(["scope", "tissue", "arm", "matched_importance_interpretation"], observed=True)
        .size()
        .rename("genes")
        .reset_index()
    )
    if counts.empty:
        return
    units = sorted((counts["scope"] + ":" + counts["tissue"]).unique())
    figure, axes = plt.subplots(1, 2, figsize=(11.5, max(7.0, 0.31 * len(units))), sharey=True)
    colors = {
        "synthetic_promoted_importance": "#E76F51",
        "reinforced_importance": "#2A9D8F",
    }
    for axis, arm in zip(axes, SYNTHETIC_ARMS):
        frame = counts.loc[counts["arm"].eq(arm)].assign(
            unit=lambda values: values["scope"] + ":" + values["tissue"]
        )
        pivot = frame.pivot_table(
            index="unit",
            columns="matched_importance_interpretation",
            values="genes",
            aggfunc="sum",
            fill_value=0,
        ).reindex(units, fill_value=0)
        left = np.zeros(len(units))
        for category in ("reinforced_importance", "synthetic_promoted_importance"):
            values = (
                pivot[category].to_numpy(dtype=float)
                if category in pivot
                else np.zeros(len(units))
            )
            axis.barh(
                np.arange(len(units)),
                values,
                left=left,
                color=colors[category],
                label=category.replace("_importance", "").replace("_", " "),
            )
            left += values
        axis.set_yticks(np.arange(len(units)), units, fontsize=7)
        axis.set_xlabel("BH-FDR genes")
        axis.set_title(MATCHED_ARM_LABELS[arm], weight="bold")
        axis.legend(frameon=False, fontsize=8)
    figure.suptitle(
        "BH-FDR genes with eligible matched importance support",
        fontsize=14,
        weight="bold",
    )
    figure.subplots_adjust(left=0.18, right=0.98, top=0.92, bottom=0.09, wspace=0.20)
    figure.savefig(output / "bh_fdr_matched_importance_counts.png", dpi=220)
    figure.savefig(output / "bh_fdr_matched_importance_counts.pdf")
    plt.close(figure)


def _write_readme(output: Path, summary: dict[str, Any]) -> None:
    text = f"""# Matched all-gene classifier analysis

This workflow compares real-only, generated-only, and real-plus-generated ridge
classifiers using the same 974 genes, preprocessing, outer splits, and model
family. One regularization value is selected from inner real-only data and then
applied unchanged to every arm in that split.

Permutation importance is evaluated on identical held-out real profiles while
shuffling within accession. Generated-domain importance is secondary and is used
to identify generator-specific patterns. Linear SHAP uses the same real-training
background for every arm.

## Scope

- Analysis units: {summary['completed_units']}
- Repeated outer splits: {summary['repeats']}
- Genes per classifier: {summary['genes']}
- Fitted outer classifiers: {summary['fitted_classifiers']}
- Permutations per gene and fit: {summary['permutation_repeats']}
- Minimum mean AUROC permutation loss: {summary['minimum_permutation_roc_auc']}
- Minimum positive outer-split fraction: {summary['minimum_positive_outer_fraction']}
- DDIM retraining: no

## Outputs

- `nested_metrics.tsv.gz`: held-out real performance for every split and arm.
- `arm_utility.tsv`: paired synthetic-arm changes from real-only.
- `importance_summary.tsv.gz`: all-gene permutation and SHAP summaries.
- `arm_gene_comparison.tsv.gz`: matched real-only versus synthetic-arm importance.
- `bh_fdr_matched_importance.tsv.gz`: real BH-FDR genes joined to matched importance.
- `eligible_bh_fdr_candidates.tsv`: compact synthetic-supported BH-FDR candidates.
- `eligible_bh_fdr_reactome.tsv.gz`: Reactome enrichment of retained candidates.
- `matched_classifier_metric_deltas.png`: tissue-level performance changes.
- `<scope>/<tissue>/matched_classifier_importance.png`: importance comparison.

## Limits

Outer profiles are held out from classifier fitting and regularization selection.
The DDIM was fixed before these repeated classifier splits and had seen the
original generator-training role, so this is not independent generator
validation. Generated profiles are not biological replicates and do not enter
the BH-FDR analysis.
"""
    (output / "README.md").write_text(text, encoding="utf-8")


def _run_source(
    source: dict[str, Any],
    *,
    regularization_grid: list[float],
    synthetic_weight: float,
    permutation_repeats: int,
    permutation_seed: int,
    units_override: set[str] | None,
    metrics_only: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    scope = str(source["scope"])
    workflow = yaml.safe_load(
        Path(source["workflow_config"]).read_text(encoding="utf-8")
    )
    analysis_dir = Path(source["analysis_dir"])
    inventory = pd.read_csv(analysis_dir / "tissue_inventory.tsv", sep="\t")
    completed = set(
        inventory.loc[inventory["status"].eq("complete"), "tissue"].astype(str)
    )
    if units_override is not None:
        completed &= units_override
    data = _load_data(workflow)
    if scope == "muscle_group":
        data, _ = _muscle_group_analysis_data(
            data, workflow["analysis"].get("muscle_groups")
        )
    elif scope != "tissue":
        raise ValueError(f"Unsupported matched-classifier scope: {scope}")
    unit_order = inventory["tissue"].astype(str).tolist()
    base_seed = int(workflow["run"]["seed"])
    metrics: list[pd.DataFrame] = []
    candidates: list[pd.DataFrame] = []
    importance: list[pd.DataFrame] = []
    completed_units: list[str] = []
    for position, tissue in enumerate(unit_order):
        if tissue not in completed:
            continue
        print(f"[matched-all-gene] {scope}:{tissue}", flush=True)
        unit_metrics, unit_candidates, unit_importance = _run_unit(
            scope=scope,
            tissue=tissue,
            data=data,
            repeats=int(workflow["analysis"]["repeats"]),
            tissue_seed=base_seed + position * 100_000,
            outer_fraction=float(workflow["analysis"]["outer_fraction"]),
            inner_fraction=float(workflow["analysis"]["inner_fraction"]),
            regularization_grid=regularization_grid,
            synthetic_weight=synthetic_weight,
            permutation_repeats=permutation_repeats,
            permutation_seed=permutation_seed,
            metrics_only=metrics_only,
        )
        metrics.append(unit_metrics)
        candidates.append(unit_candidates)
        if not unit_importance.empty:
            importance.append(unit_importance)
        completed_units.append(tissue)
    if not metrics:
        raise RuntimeError(f"No units completed for scope {scope}")
    return (
        pd.concat(metrics, ignore_index=True),
        pd.concat(candidates, ignore_index=True),
        pd.concat(importance, ignore_index=True) if importance else pd.DataFrame(),
        completed_units,
    )


def run(
    config_path: Path,
    *,
    scopes_override: set[str] | None = None,
    units_override: set[str] | None = None,
    permutation_repeats_override: int | None = None,
    output_override: Path | None = None,
    metrics_only: bool = False,
) -> Path:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    output = output_override or Path(config["run"]["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    options = config["analysis"]
    regularization_grid = list(map(float, options["regularization_c"]))
    permutation_repeats = int(
        permutation_repeats_override or options["permutation_repeats"]
    )
    synthetic_weight = float(options.get("synthetic_weight", 1.0))
    minimum_importance = float(options.get("minimum_permutation_roc_auc", 0.001))
    minimum_fraction = float(options.get("minimum_positive_outer_fraction", 0.5))
    permutation_seed = int(config["run"]["seed"])
    metric_tables: list[pd.DataFrame] = []
    candidate_tables: list[pd.DataFrame] = []
    importance_tables: list[pd.DataFrame] = []
    completed_units: list[tuple[str, str]] = []
    for source in config["sources"]:
        scope = str(source["scope"])
        if scopes_override is not None and scope not in scopes_override:
            continue
        metrics, candidates, importance, units = _run_source(
            source,
            regularization_grid=regularization_grid,
            synthetic_weight=synthetic_weight,
            permutation_repeats=permutation_repeats,
            permutation_seed=permutation_seed,
            units_override=units_override,
            metrics_only=metrics_only,
        )
        metric_tables.append(metrics)
        candidate_tables.append(candidates)
        importance_tables.append(importance)
        completed_units.extend((scope, tissue) for tissue in units)
    if not metric_tables:
        raise RuntimeError("No matched all-gene source was selected")
    metrics = pd.concat(metric_tables, ignore_index=True)
    candidates = pd.concat(candidate_tables, ignore_index=True)
    importance_by_repeat = (
        pd.concat(importance_tables, ignore_index=True)
        if importance_tables and any(not table.empty for table in importance_tables)
        else pd.DataFrame()
    )
    completed = (
        metrics.groupby(["scope", "tissue"], observed=True)["repeat"]
        .nunique()
        .rename("completed_repeats")
        .reset_index()
    )
    metric_summary, utility = _metric_summary(metrics)

    metrics.to_csv(output / "nested_metrics.tsv.gz", sep="\t", index=False)
    candidates.to_csv(
        output / "inner_real_regularization_metrics.tsv.gz", sep="\t", index=False
    )
    metric_summary.to_csv(output / "arm_summary.tsv", sep="\t", index=False)
    utility.to_csv(output / "arm_utility.tsv", sep="\t", index=False)
    _plot_utility(utility, output)

    if metrics_only:
        summary = {
            "status": "complete",
            "mode": "metrics_only",
            "config": str(config_path.resolve()),
            "output": str(output.resolve()),
            "completed_units": int(len(completed)),
            "scopes": completed["scope"].value_counts().to_dict(),
            "repeats": int(metrics["repeat"].nunique()),
            "genes": int(metrics["genes"].max()),
            "fitted_classifiers": int(
                metrics[["scope", "tissue", "repeat", "arm"]]
                .drop_duplicates()
                .shape[0]
            ),
            "synthetic_weight": synthetic_weight,
            "shared_regularization_selected_from": "inner real-only profiles",
            "joint_utility_units_by_arm": {
                arm: int(
                    utility.loc[
                        utility["arm"].eq(arm),
                        "joint_mean_all_metrics_nonworse",
                    ].sum()
                )
                for arm in SYNTHETIC_ARMS
            },
            "interpretation": (
                "For the combined arm, total synthetic training weight is "
                f"{synthetic_weight:g} times total real training weight."
            ),
        }
        (output / "README.md").write_text(
            "# Matched all-gene classifier metric screen\n\n"
            "This metric-only run uses the same 974 genes, preprocessing, nested "
            "splits, real-data regularization selection, and classifier family as "
            "the full matched analysis. It skips permutation importance and SHAP.\n\n"
            f"The combined arm gives all synthetic profiles {synthetic_weight:g} "
            "times the total weight of the real profiles. All arms are evaluated "
            "on the same held-out real OSDR profiles.\n",
            encoding="utf-8",
        )
        (output / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(summary, indent=2), flush=True)
        return output / "summary.json"

    aggregate = _aggregate_importance(importance_by_repeat, completed)
    comparison = _importance_comparison(
        aggregate,
        minimum_importance=minimum_importance,
        minimum_fraction=minimum_fraction,
    )
    similarity = _importance_similarity(comparison)
    crosswalk = _bh_fdr_crosswalk(
        Path(config["annotations"]["bh_fdr_inventory"]), comparison, utility
    )

    importance_by_repeat.to_csv(
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
        output / "arm_gene_comparison.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    similarity.to_csv(output / "arm_importance_similarity.tsv", sep="\t", index=False)
    crosswalk.to_csv(
        output / "bh_fdr_matched_importance.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    eligible = crosswalk.loc[
        crosswalk["eligible_synthetic_biological_candidate"]
    ].copy()
    eligible.to_csv(
        output / "eligible_bh_fdr_candidates.tsv", sep="\t", index=False
    )
    gene_symbols = (
        aggregate[["gene", "symbol"]]
        .drop_duplicates("gene")
        .set_index("gene")["symbol"]
        .astype(str)
        .to_dict()
    )
    candidate_enrichment = _candidate_reactome_enrichment(
        eligible,
        background=sorted(aggregate["gene"].astype(str).unique()),
        gmt_path=Path(config["annotations"]["reactome_gmt"]),
        symbols=gene_symbols,
    )
    candidate_enrichment.to_csv(
        output / "eligible_bh_fdr_reactome.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    top = (
        aggregate.loc[aggregate["domain"].eq("real")]
        .sort_values(
            ["scope", "tissue", "arm", "permutation_roc_auc_mean"],
            ascending=[True, True, True, False],
            kind="stable",
        )
        .groupby(["scope", "tissue", "arm"], observed=True)
        .head(int(options.get("top_genes_per_arm", 25)))
    )
    top.to_csv(output / "top_gene_importance.tsv", sep="\t", index=False)

    _plot_bh_candidate_counts(crosswalk, output)
    for scope, tissue in completed_units:
        _plot_unit_importance(
            aggregate,
            comparison,
            output,
            scope=scope,
            tissue=tissue,
            minimum_importance=minimum_importance,
        )

    pattern_counts = comparison["pattern"].value_counts().to_dict()
    eligible_counts = (
        crosswalk.loc[crosswalk["eligible_synthetic_biological_candidate"]]
        .groupby("arm", observed=True)
        .size()
        .to_dict()
    )
    summary = {
        "status": "complete",
        "config": str(config_path.resolve()),
        "output": str(output.resolve()),
        "completed_units": int(len(completed)),
        "scopes": completed["scope"].value_counts().to_dict(),
        "repeats": int(metrics["repeat"].nunique()),
        "genes": int(importance_by_repeat["gene"].nunique()),
        "fitted_classifiers": int(
            metrics[["scope", "tissue", "repeat", "arm"]].drop_duplicates().shape[0]
        ),
        "permutation_repeats": permutation_repeats,
        "minimum_permutation_roc_auc": minimum_importance,
        "minimum_positive_outer_fraction": minimum_fraction,
        "synthetic_weight": synthetic_weight,
        "shared_regularization_selected_from": "inner real-only profiles",
        "importance_permutation_blocks": "accession within domain and DDIM draw",
        "pattern_counts": {key: int(value) for key, value in pattern_counts.items()},
        "eligible_bh_fdr_candidates_by_arm": {
            key: int(value) for key, value in eligible_counts.items()
        },
        "eligible_bh_fdr_arm_rows": int(len(eligible)),
        "eligible_bh_fdr_unique_associations": int(
            eligible[["scope", "tissue", "gene"]].drop_duplicates().shape[0]
        ),
        "eligible_bh_fdr_units": sorted(
            (eligible["scope"] + ":" + eligible["tissue"]).unique().tolist()
        ),
        "significant_candidate_reactome_terms": int(
            candidate_enrichment["fdr"].lt(0.05).sum()
            if not candidate_enrichment.empty
            else 0
        ),
        "joint_utility_units_by_arm": {
            arm: int(
                utility.loc[
                    utility["arm"].eq(arm), "joint_mean_all_metrics_nonworse"
                ].sum()
            )
            for arm in SYNTHETIC_ARMS
        },
        "limitations": [
            "The DDIM was fixed before and not nested inside classifier splits.",
            "The original generator-training role can appear in classifier outer tests.",
            "Repeated outer splits overlap and are not independent observations.",
            "Generated profiles are not biological replicates.",
            "Permutation importance can be diluted among correlated genes.",
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
    parser.add_argument(
        "--metrics-only",
        action="store_true",
        help="Skip permutation importance and SHAP; write classifier metrics only.",
    )
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    run(
        arguments.config,
        scopes_override=set(arguments.scopes) if arguments.scopes else None,
        units_override=set(arguments.units) if arguments.units else None,
        permutation_repeats_override=arguments.permutation_repeats,
        output_override=arguments.output,
        metrics_only=arguments.metrics_only,
    )


if __name__ == "__main__":
    main()
