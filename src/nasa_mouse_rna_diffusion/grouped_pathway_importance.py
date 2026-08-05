"""Evaluate Reactome pathways as grouped features in matched FLT/GC classifiers."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.preprocessing import StandardScaler

from nasa_mouse_generative.effect_validation import accession_effects, random_effects_table

from .classifier_importance import _metric_matrix_from_logits
from .matched_all_gene_classifiers import (
    ARMS,
    MATCHED_ARM_LABELS,
    SYNTHETIC_ARMS,
    _accession_blocks,
    _evaluate_classifier,
    _fit_matched_arm,
    _metric_summary,
    _select_shared_regularization,
)
from .within_study_feature_stability import (
    METRICS,
    WorkflowData,
    _labels,
    _load_data,
    _metric_set,
    _muscle_group_analysis_data,
    _valid_nested_split,
)


@dataclass(frozen=True)
class PathwayGroup:
    """One Reactome pathway intersected with the classifier feature space."""

    term: str
    description: str
    url: str
    indices: tuple[int, ...]
    genes: tuple[str, ...]
    symbols: tuple[str, ...]


def _load_pathway_groups(
    gmt_path: Path,
    genes: list[str],
    symbols: dict[str, str],
    *,
    minimum_genes: int,
    maximum_genes: int,
    term_metadata_path: Path | None = None,
) -> list[PathwayGroup]:
    gene_index = {str(gene): index for index, gene in enumerate(genes)}
    term_metadata: dict[str, dict[str, str]] = {}
    if term_metadata_path is not None:
        metadata = pd.read_csv(term_metadata_path, sep="\t").fillna("")
        term_metadata = metadata.set_index("term")[["name", "url"]].to_dict(
            orient="index"
        )
    merged: dict[str, dict[str, object]] = {}
    with gmt_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            term, description, *members = fields
            record = merged.setdefault(
                term,
                {"description": description, "genes": set()},
            )
            record["genes"].update(gene for gene in members if gene in gene_index)

    groups: list[PathwayGroup] = []
    for term, record in sorted(merged.items()):
        members = sorted(record["genes"], key=gene_index.__getitem__)
        if not minimum_genes <= len(members) <= maximum_genes:
            continue
        groups.append(
            PathwayGroup(
                term=str(term),
                description=term_metadata.get(term, {}).get(
                    "name", term.split("_", 1)[-1].replace("_", " ").title()
                ),
                url=term_metadata.get(term, {}).get(
                    "url", str(record["description"])
                ),
                indices=tuple(gene_index[gene] for gene in members),
                genes=tuple(members),
                symbols=tuple(symbols.get(gene, gene) for gene in members),
            )
        )
    if not groups:
        raise ValueError(f"No eligible pathways overlap {gmt_path}")
    return groups


def _membership_matrix(
    number_of_genes: int, groups: list[PathwayGroup]
) -> np.ndarray:
    membership = np.zeros((number_of_genes, len(groups)), dtype=np.float64)
    for column, group in enumerate(groups):
        membership[np.asarray(group.indices, dtype=np.int64), column] = 1.0
    return membership


def _validate_blocks(blocks: list[np.ndarray], number_of_rows: int) -> None:
    observed = sorted(np.concatenate(blocks).tolist()) if blocks else []
    if observed != list(range(number_of_rows)):
        raise ValueError("Permutation blocks must partition evaluation rows")


def _group_permutation_rows(
    classifier,
    expression: np.ndarray,
    labels: np.ndarray,
    groups: list[PathwayGroup],
    membership: np.ndarray,
    *,
    permutation_repeats: int,
    seed: int,
    blocks: list[np.ndarray],
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Jointly permute all pathway genes with one row order per accession.

    For a linear model, jointly permuting a pathway is exactly equivalent to
    permuting its summed logit contribution. This preserves the covariance among
    pathway genes while breaking the pathway-outcome relation.
    """

    expression = np.asarray(expression, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int64)
    if membership.shape != (expression.shape[1], len(groups)):
        raise ValueError("Pathway membership and expression dimensions differ")
    if len(labels) != len(expression):
        raise ValueError("Permutation expression and labels differ")
    _validate_blocks(blocks, len(expression))

    probability = classifier.predict_proba(expression)[:, 1]
    baseline = _metric_set(labels, probability)
    baseline_logit = np.asarray(classifier.decision_function(expression), dtype=float)
    coefficients = np.asarray(classifier.coef_[0], dtype=float)
    contributions = (expression * coefficients.reshape(1, -1)) @ membership
    rng = np.random.default_rng(int(seed))
    rows: list[dict[str, object]] = []
    for column, group in enumerate(groups):
        original = contributions[:, column]
        permuted = np.tile(original, (int(permutation_repeats), 1))
        for block in blocks:
            orders = np.argsort(
                rng.random((int(permutation_repeats), len(block))),
                axis=1,
                kind="stable",
            )
            permuted[:, block] = original[block][orders]
        logits = baseline_logit.reshape(1, -1) + permuted - original.reshape(1, -1)
        permuted_metrics = _metric_matrix_from_logits(labels, logits)
        row: dict[str, object] = {
            "term": group.term,
            "description": group.description,
            "url": group.url,
            "pathway_genes": len(group.genes),
            "genes": ",".join(group.genes),
            "symbols": ",".join(group.symbols),
        }
        for metric in METRICS:
            losses = baseline[metric] - permuted_metrics[metric]
            row[f"baseline_{metric}"] = baseline[metric]
            row[f"permutation_{metric}_mean"] = float(losses.mean())
            row[f"permutation_{metric}_sd"] = float(
                losses.std(ddof=1) if len(losses) > 1 else 0.0
            )
            row[f"permutation_{metric}_positive_fraction"] = float(
                np.mean(losses > 0.0)
            )
        rows.append(row)
    return pd.DataFrame(rows), baseline


def _group_linear_shap_rows(
    classifier,
    expression: np.ndarray,
    labels: np.ndarray,
    background: np.ndarray,
    groups: list[PathwayGroup],
    membership: np.ndarray,
) -> pd.DataFrame:
    """Sum exact linear interventional SHAP values over pathway members."""

    expression = np.asarray(expression, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int64)
    background = np.asarray(background, dtype=np.float64)
    coefficients = np.asarray(classifier.coef_[0], dtype=np.float64)
    if expression.shape[1] != len(background) or membership.shape[0] != len(background):
        raise ValueError("SHAP background, expression, and pathways differ")
    gene_values = (
        expression - background.reshape(1, -1)
    ) * coefficients.reshape(1, -1)
    group_values = gene_values @ membership
    flight = labels == 1
    ground = labels == 0
    rows: list[dict[str, object]] = []
    for column, group in enumerate(groups):
        contribution = group_values[:, column]
        flight_mean = float(contribution[flight].mean())
        ground_mean = float(contribution[ground].mean())
        rows.append(
            {
                "term": group.term,
                "group_shap_mean_absolute": float(np.mean(np.abs(contribution))),
                "group_shap_mean_flight": flight_mean,
                "group_shap_mean_ground_control": ground_mean,
                "group_shap_flight_minus_ground": flight_mean - ground_mean,
            }
        )
    return pd.DataFrame(rows)


def _within_accession_pathway_scores(
    expression: np.ndarray,
    samples: pd.DataFrame,
    membership: np.ndarray,
) -> np.ndarray:
    """Average within-accession gene z scores for each Reactome pathway."""

    expression = np.asarray(expression, dtype=np.float64)
    if len(expression) != len(samples) or expression.shape[1] != membership.shape[0]:
        raise ValueError("Pathway score inputs do not align")
    standardized = np.zeros_like(expression, dtype=np.float64)
    for _, indices in samples.groupby("accession", sort=True, observed=True).groups.items():
        positions = np.asarray(list(indices), dtype=np.int64)
        block = expression[positions]
        scale = block.std(axis=0, ddof=1)
        scale[~np.isfinite(scale) | (scale == 0.0)] = 1.0
        standardized[positions] = (block - block.mean(axis=0)) / scale
    sizes = membership.sum(axis=0)
    return standardized @ (membership / sizes.reshape(1, -1))


def _pathway_associations(
    *,
    scope: str,
    tissue: str,
    data: WorkflowData,
    groups: list[PathwayGroup],
    membership: np.ndarray,
) -> pd.DataFrame:
    mask = data.all_samples["tissue"].astype(str).eq(tissue)
    positions = np.flatnonzero(mask.to_numpy())
    samples = data.all_samples.loc[mask].reset_index(drop=True)
    expression = data.all_expression[positions]
    scores = _within_accession_pathway_scores(expression, samples, membership)
    terms = [group.term for group in groups]
    effects = accession_effects(scores, samples, terms)
    meta = random_effects_table(effects)
    metadata = pd.DataFrame(
        {
            "term": terms,
            "description": [group.description for group in groups],
            "url": [group.url for group in groups],
            "pathway_genes": [len(group.genes) for group in groups],
            "genes": [",".join(group.genes) for group in groups],
            "symbols": [",".join(group.symbols) for group in groups],
        }
    )
    if meta.empty:
        result = metadata.copy()
        for column in (
            "n_accessions",
            "meta_effect",
            "meta_se",
            "meta_p",
            "tau2",
            "i2",
            "n_accession_same_direction",
            "n_accession_opposite_direction",
            "meta_fdr",
        ):
            result[column] = np.nan
    else:
        result = metadata.merge(
            meta.rename(columns={"feature": "term"}),
            on="term",
            how="left",
            validate="one_to_one",
        )
    result.insert(0, "tissue", tissue)
    result.insert(0, "scope", scope)
    result["accession_direction_fraction"] = (
        result["n_accession_same_direction"] / result["n_accessions"]
    )
    result["flt_gc_direction"] = np.where(
        result["meta_effect"].gt(0.0), "FLT_higher", "FLT_lower"
    )
    return result


def _run_unit(
    *,
    scope: str,
    tissue: str,
    data: WorkflowData,
    groups: list[PathwayGroup],
    membership: np.ndarray,
    repeats: int,
    tissue_seed: int,
    outer_fraction: float,
    inner_fraction: float,
    regularization_grid: list[float],
    synthetic_weight: float,
    permutation_repeats: int,
    permutation_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = data.development_samples["tissue"].astype(str).eq(tissue)
    positions = np.flatnonzero(mask.to_numpy())
    real = data.development_expression[positions]
    samples = data.development_samples.loc[mask].reset_index(drop=True)
    labels = _labels(samples)
    synthetic = [draw[positions] for draw in data.synthetic_draws.values()]
    pathway_tables: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []
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
        selected_c, _ = _select_shared_regularization(
            real,
            labels,
            samples,
            inner_train=inner_train,
            inner_validation=inner_validation,
            regularization_grid=regularization_grid,
            seed=tissue_seed + repeat * 100,
        )
        scaler = StandardScaler().fit(real[outer_train])
        real_train = scaler.transform(real[outer_train])
        real_test = scaler.transform(real[outer_test])
        synthetic_train = [scaler.transform(draw[outer_train]) for draw in synthetic]
        real_test_samples = samples.loc[outer_test].reset_index(drop=True)
        background = real_train.mean(axis=0)
        blocks = _accession_blocks(real_test_samples)
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
            permutation, baseline = _group_permutation_rows(
                classifier,
                real_test,
                labels[outer_test],
                groups,
                membership,
                permutation_repeats=permutation_repeats,
                seed=(
                    permutation_seed
                    + tissue_seed
                    + repeat * 10_000
                    + arm_offset * 100
                ),
                blocks=blocks,
            )
            shap = _group_linear_shap_rows(
                classifier,
                real_test,
                labels[outer_test],
                background,
                groups,
                membership,
            )
            table = permutation.merge(shap, on="term", validate="one_to_one")
            table.insert(0, "regularization_c", selected_c)
            table.insert(0, "arm", arm)
            table.insert(0, "repeat", repeat)
            table.insert(0, "tissue", tissue)
            table.insert(0, "scope", scope)
            pathway_tables.append(table)
            evaluation_metrics = _evaluate_classifier(
                classifier,
                real_test,
                labels[outer_test],
                real_test_samples,
            )
            for metric in METRICS:
                if not np.isclose(evaluation_metrics[metric], baseline[metric]):
                    raise RuntimeError(
                        f"Grouped and classifier {metric} baselines differ"
                    )
            metric_rows.append(
                {
                    "scope": scope,
                    "tissue": tissue,
                    "repeat": repeat,
                    "arm": arm,
                    "regularization_c": selected_c,
                    **evaluation_metrics,
                }
            )
    if not pathway_tables:
        raise RuntimeError(f"No valid grouped-pathway splits for {scope}:{tissue}")
    return pd.concat(pathway_tables, ignore_index=True), pd.DataFrame(metric_rows)


def _aggregate_pathways(
    rows: pd.DataFrame, completed_repeats: pd.DataFrame
) -> pd.DataFrame:
    grouped = rows.groupby(
        ["scope", "tissue", "arm", "term"], sort=True, observed=True
    )
    specification: dict[str, tuple[str, Any]] = {
        "description": ("description", "first"),
        "url": ("url", "first"),
        "pathway_genes": ("pathway_genes", "first"),
        "genes": ("genes", "first"),
        "symbols": ("symbols", "first"),
        "evaluated_repeats": ("repeat", "nunique"),
        "median_regularization_c": ("regularization_c", "median"),
    }
    for metric in METRICS:
        for suffix in ("mean", "sd", "positive_fraction"):
            column = f"permutation_{metric}_{suffix}"
            specification[column] = (column, "mean")
        specification[f"permutation_{metric}_positive_repeat_fraction"] = (
            f"permutation_{metric}_mean",
            lambda values: float(np.mean(np.asarray(values) > 0.0)),
        )
    for column in (
        "group_shap_mean_absolute",
        "group_shap_mean_flight",
        "group_shap_mean_ground_control",
        "group_shap_flight_minus_ground",
    ):
        specification[column] = (column, "mean")
    table = grouped.agg(**specification).reset_index()
    return table.merge(
        completed_repeats,
        on=["scope", "tissue"],
        how="left",
        validate="many_to_one",
    )


def _positive_group_importance(
    frame: pd.DataFrame,
    prefix: str,
    *,
    minimum_importance: float,
    minimum_fraction: float,
) -> pd.Series:
    return frame[f"{prefix}_permutation_roc_auc_mean"].ge(minimum_importance) & frame[
        f"{prefix}_permutation_roc_auc_positive_repeat_fraction"
    ].ge(minimum_fraction)


def _compare_arms(
    aggregate: pd.DataFrame,
    utility: pd.DataFrame,
    associations: pd.DataFrame,
    *,
    minimum_importance: float,
    minimum_fraction: float,
) -> pd.DataFrame:
    keys = ["scope", "tissue", "term"]
    values = [
        "permutation_balanced_accuracy_mean",
        "permutation_balanced_accuracy_positive_repeat_fraction",
        "permutation_roc_auc_mean",
        "permutation_roc_auc_positive_repeat_fraction",
        "permutation_average_precision_mean",
        "permutation_average_precision_positive_repeat_fraction",
        "group_shap_mean_absolute",
        "group_shap_flight_minus_ground",
    ]
    baseline = aggregate.loc[
        aggregate["arm"].eq("real_only"), keys + values
    ].rename(columns={column: f"real_only_{column}" for column in values})
    tables: list[pd.DataFrame] = []
    for arm in SYNTHETIC_ARMS:
        candidate = aggregate.loc[
            aggregate["arm"].eq(arm), keys + values
        ].rename(columns={column: f"arm_{column}" for column in values})
        table = baseline.merge(candidate, on=keys, validate="one_to_one")
        table.insert(2, "arm", arm)
        table["real_only_positive_group_importance"] = _positive_group_importance(
            table,
            "real_only",
            minimum_importance=minimum_importance,
            minimum_fraction=minimum_fraction,
        )
        table["arm_positive_group_importance"] = _positive_group_importance(
            table,
            "arm",
            minimum_importance=minimum_importance,
            minimum_fraction=minimum_fraction,
        )
        real_positive = table["real_only_positive_group_importance"]
        arm_positive = table["arm_positive_group_importance"]
        reinforced = arm_positive & real_positive & table[
            "arm_permutation_roc_auc_mean"
        ].ge(table["real_only_permutation_roc_auc_mean"])
        table["group_importance_pattern"] = np.select(
            [
                arm_positive & ~real_positive,
                reinforced,
                arm_positive & real_positive & ~reinforced,
                ~arm_positive & real_positive,
            ],
            [
                "synthetic_promoted_group",
                "shared_reinforced_group",
                "shared_attenuated_group",
                "real_only_group",
            ],
            default="below_group_importance_gate",
        )
        tables.append(table)
    result = pd.concat(tables, ignore_index=True)
    association_columns = [
        "scope",
        "tissue",
        "term",
        "description",
        "url",
        "pathway_genes",
        "genes",
        "symbols",
        "n_accessions",
        "meta_effect",
        "meta_se",
        "meta_p",
        "meta_fdr",
        "tau2",
        "i2",
        "n_accession_same_direction",
        "n_accession_opposite_direction",
        "accession_direction_fraction",
        "flt_gc_direction",
    ]
    result = result.merge(
        associations[association_columns],
        on=keys,
        how="left",
        validate="many_to_one",
    ).merge(
        utility[
            [
                "scope",
                "tissue",
                "arm",
                "joint_mean_all_metrics_nonworse",
            ]
        ],
        on=["scope", "tissue", "arm"],
        how="left",
        validate="many_to_one",
    )
    result["group_shap_supports_flt_gc"] = result[
        "arm_group_shap_flight_minus_ground"
    ].gt(0.0)
    result["eligible_synthetic_pathway"] = (
        result["joint_mean_all_metrics_nonworse"].fillna(False)
        & result["group_importance_pattern"].isin(
            ["synthetic_promoted_group", "shared_reinforced_group"]
        )
        & result["group_shap_supports_flt_gc"]
        & result["meta_fdr"].lt(0.05)
    )
    return result


def _nonredundant_top_pathways(
    comparison: pd.DataFrame,
    *,
    maximum_per_arm: int,
    maximum_jaccard: float,
) -> pd.DataFrame:
    eligible = comparison.loc[comparison["eligible_synthetic_pathway"]].copy()
    selected_rows: list[pd.Series] = []
    for _, frame in eligible.groupby(
        ["scope", "tissue", "arm"], sort=True, observed=True
    ):
        ordered = frame.sort_values(
            ["meta_fdr", "arm_permutation_roc_auc_mean", "term"],
            ascending=[True, False, True],
            kind="stable",
        )
        selected_sets: list[set[str]] = []
        for _, row in ordered.iterrows():
            genes = set(str(row["genes"]).split(","))
            if any(
                len(genes & prior) / len(genes | prior) > maximum_jaccard
                for prior in selected_sets
            ):
                continue
            selected_rows.append(row)
            selected_sets.append(genes)
            if len(selected_sets) >= maximum_per_arm:
                break
    return pd.DataFrame(selected_rows).reset_index(drop=True)


def _plot_unit(
    aggregate: pd.DataFrame,
    output: Path,
    *,
    scope: str,
    tissue: str,
    top_count: int,
) -> None:
    frame = aggregate.loc[
        aggregate["scope"].eq(scope) & aggregate["tissue"].eq(tissue)
    ]
    if frame.empty:
        return
    top_terms = (
        frame.groupby(["term", "description"], observed=True)[
            "permutation_roc_auc_mean"
        ]
        .max()
        .nlargest(top_count)
        .index
    )
    terms = [term for term, _ in top_terms]
    descriptions = [description for _, description in top_terms]
    labels = [
        description if len(description) <= 58 else description[:55] + "..."
        for description in descriptions
    ]
    pivot = frame.pivot_table(
        index="term",
        columns="arm",
        values="permutation_roc_auc_mean",
    ).reindex(index=terms, columns=ARMS)
    positions = np.arange(len(terms))
    width = 0.24
    colors = {
        "real_only": "#457B9D",
        "generated_only": "#E9C46A",
        "real_plus_generated": "#E76F51",
    }
    figure, axis = plt.subplots(figsize=(12.5, max(6.0, 0.42 * len(terms))))
    for index, arm in enumerate(ARMS):
        axis.barh(
            positions + (index - 1) * width,
            pivot[arm].to_numpy(dtype=float),
            height=width,
            color=colors[arm],
            label=MATCHED_ARM_LABELS[arm],
        )
    axis.axvline(0.0, color="#333333", linewidth=0.8)
    axis.set_yticks(positions, labels, fontsize=8)
    axis.invert_yaxis()
    axis.set_xlabel("Held-out real AUROC loss after joint pathway permutation")
    axis.set_title(
        f"{tissue.replace('_', ' ').title()}: grouped Reactome importance",
        weight="bold",
    )
    axis.legend(frameon=False, ncol=3, loc="lower right")
    figure.tight_layout()
    directory = output / scope / tissue
    directory.mkdir(parents=True, exist_ok=True)
    figure.savefig(directory / "grouped_pathway_importance.png", dpi=220)
    figure.savefig(directory / "grouped_pathway_importance.pdf")
    plt.close(figure)


def _plot_eligible_summary(pathways: pd.DataFrame, output: Path) -> None:
    if pathways.empty:
        return
    frame = pathways.sort_values(
        ["scope", "tissue", "meta_fdr", "arm", "description"], kind="stable"
    ).reset_index(drop=True)
    arm_abbreviation = {
        "generated_only": "synthetic",
        "real_plus_generated": "real + synthetic",
    }
    labels = [
        f"{row.tissue.replace('_', ' ').title()}: {row.description} "
        f"[{arm_abbreviation.get(row.arm, row.arm)}]"
        for row in frame.itertuples()
    ]
    labels = [label if len(label) <= 82 else label[:79] + "..." for label in labels]
    colors = [
        "#D37B00" if direction == "FLT_higher" else "#2D6496"
        for direction in frame["flt_gc_direction"]
    ]
    positions = np.arange(len(frame))
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(16.0, max(6.4, 0.56 * len(frame))),
        sharey=True,
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.0]},
    )
    panels = [
        (
            "arm_permutation_roc_auc_mean",
            "Held-out AUROC loss\nafter group permutation",
            False,
        ),
        (
            "arm_group_shap_flight_minus_ground",
            "Grouped SHAP\nFLT - GC contribution",
            False,
        ),
        (
            "meta_effect",
            "Observed pathway score\nFLT - GC effect",
            True,
        ),
    ]
    for axis, (column, xlabel, diverging) in zip(axes, panels):
        axis.barh(positions, frame[column], color=colors, alpha=0.90)
        axis.axvline(0.0, color="#333333", linewidth=0.8)
        axis.set_xlabel(xlabel)
        axis.grid(axis="x", color="#E4E8EB", linewidth=0.7)
        axis.set_axisbelow(True)
        if diverging:
            limit = max(abs(float(frame[column].min())), abs(float(frame[column].max())))
            axis.set_xlim(-1.12 * limit, 1.12 * limit)
    axes[0].set_yticks(positions, labels, fontsize=8)
    axes[0].invert_yaxis()
    figure.suptitle(
        "Synthetic-supported grouped Reactome pathways",
        fontsize=15,
        weight="bold",
    )
    figure.text(
        0.52,
        0.02,
        "Blue: lower pathway score in FLT. Orange: higher pathway score in FLT. "
        "All association FDR values are below 0.05 in observed OSDR profiles.",
        ha="center",
        fontsize=9,
    )
    figure.subplots_adjust(left=0.39, right=0.98, top=0.90, bottom=0.14, wspace=0.24)
    figure.savefig(output / "eligible_grouped_pathway_evidence.png", dpi=220)
    figure.savefig(output / "eligible_grouped_pathway_evidence.pdf")
    plt.close(figure)


def _write_readme(
    output: Path,
    summary: dict[str, Any],
    top_pathways: pd.DataFrame,
) -> None:
    result_rows: list[str] = []
    for row in top_pathways.itertuples():
        result_rows.append(
            "| "
            + " | ".join(
                [
                    str(row.tissue).replace("_", " "),
                    MATCHED_ARM_LABELS.get(str(row.arm), str(row.arm)),
                    str(row.description),
                    str(row.flt_gc_direction).replace("_", " "),
                    f"{float(row.arm_permutation_roc_auc_mean):.4f}",
                    f"{float(row.arm_group_shap_flight_minus_ground):.4f}",
                    f"{float(row.meta_fdr):.4g}",
                ]
            )
            + " |"
        )
    result_table = "\n".join(result_rows) if result_rows else "No pathway passed every gate."
    text = f"""# Grouped Reactome importance

This analysis refits the completed matched all-gene classifiers and evaluates
Reactome pathways as joint feature groups on the same held-out real profiles.
It does not retrain the DDIM.

For pathway permutation, every member gene receives the same row permutation
within each accession. This preserves within-pathway gene covariance while
breaking the pathway-outcome relationship. For the linear ridge classifier this
is computed exactly by permuting the pathway's summed logit contribution.

Grouped linear SHAP is the exact sum of member-gene interventional SHAP values.
It describes the pathway's contribution to the fitted logit; overlapping
Reactome pathways are not independent and their SHAP values must not be summed
across pathways.

Real-data pathway association is evaluated separately. Genes are standardized
within accession, pathway scores are their mean, and FLT-GC effects are combined
with accession-aware random effects. BH FDR is calculated within each analysis
unit over all eligible Reactome pathways. Generated profiles never enter this
association test.

## Scope

- Analysis units: {summary['completed_units']}
- Eligible Reactome pathways: {summary['pathways']}
- Classifier genes: {summary['genes']}
- Outer repeats: {summary['repeats']}
- Joint permutations per pathway and fit: {summary['permutation_repeats']}
- Minimum AUROC permutation loss: {summary['minimum_permutation_roc_auc']}
- Minimum positive outer-repeat fraction: {summary['minimum_positive_outer_fraction']}

## Main results

Eligible pathway groups occurred in {', '.join(unit.split(':', 1)[1].replace('_', ' ') for unit in summary['eligible_units']) or 'no analysis unit'}.
Reactome contains overlapping parent and child terms, so the rows below are a
Jaccard-filtered interpretation set rather than independent discoveries.

| Tissue | Classifier arm | Reactome pathway | Observed direction | AUROC loss | Grouped SHAP FLT-GC | Real BH FDR |
|---|---|---|---|---:|---:|---:|
{result_table}

## Outputs

- `pathway_importance_by_repeat.tsv.gz`: pathway permutation and grouped SHAP for every fit.
- `pathway_importance_summary.tsv.gz`: values aggregated over outer repeats.
- `arm_utility.tsv`: same-run synthetic-arm changes from the real-only classifier.
- `real_pathway_random_effects.tsv.gz`: observed-data pathway-score association and BH FDR.
- `pathway_arm_comparison.tsv.gz`: matched real-only versus synthetic-arm results.
- `eligible_synthetic_pathways.tsv.gz`: pathways passing utility, grouped importance, SHAP-direction, and real BH-FDR gates.
- `top_nonredundant_pathways.tsv`: a Jaccard-filtered interpretation table.
- `eligible_grouped_pathway_evidence.png`: grouped permutation, SHAP, and observed-effect summary.
- `<scope>/<tissue>/grouped_pathway_importance.png`: top pathway importance by arm.

## Interpretation limits

Grouped importance tests whether a pathway contributes collectively, not which
member is causal. Pathways overlap extensively, pathway size affects the amount
of information removed, and repeated outer splits are not independent studies.
The pathway-score FDR is real-data evidence; grouped permutation and SHAP are
predictive interpretation. Independent biological replication remains required.
"""
    (output / "README.md").write_text(text, encoding="utf-8")


def run(
    config_path: Path,
    *,
    scopes_override: set[str] | None = None,
    units_override: set[str] | None = None,
    permutation_repeats_override: int | None = None,
    output_override: Path | None = None,
) -> Path:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    matched_config_path = Path(config["matched_config"])
    matched_config = yaml.safe_load(matched_config_path.read_text(encoding="utf-8"))
    output = output_override or Path(config["run"]["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    options = config["analysis"]
    permutation_repeats = int(
        permutation_repeats_override or options["permutation_repeats"]
    )
    regularization_grid = list(
        map(float, matched_config["analysis"]["regularization_c"])
    )
    synthetic_weight = float(matched_config["analysis"].get("synthetic_weight", 1.0))
    permutation_seed = int(config["run"]["seed"])
    minimum_importance = float(options["minimum_permutation_roc_auc"])
    minimum_fraction = float(options["minimum_positive_outer_fraction"])

    importance_tables: list[pd.DataFrame] = []
    metric_tables: list[pd.DataFrame] = []
    association_tables: list[pd.DataFrame] = []
    completed_units: list[tuple[str, str]] = []
    pathway_terms: set[str] | None = None
    pathway_count = 0
    gene_count = 0
    for source in matched_config["sources"]:
        scope = str(source["scope"])
        if scopes_override is not None and scope not in scopes_override:
            continue
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
            raise ValueError(f"Unsupported grouped-pathway scope: {scope}")
        groups = _load_pathway_groups(
            Path(matched_config["annotations"]["reactome_gmt"]),
            data.genes,
            data.symbols,
            minimum_genes=int(options["minimum_pathway_genes"]),
            maximum_genes=int(options["maximum_pathway_genes"]),
            term_metadata_path=Path(config["reactome_terms"]),
        )
        terms = {group.term for group in groups}
        if pathway_terms is not None and terms != pathway_terms:
            raise ValueError("Tissue and muscle-group pathway feature spaces differ")
        pathway_terms = terms
        pathway_count = len(groups)
        gene_count = len(data.genes)
        membership = _membership_matrix(len(data.genes), groups)
        unit_order = inventory["tissue"].astype(str).tolist()
        base_seed = int(workflow["run"]["seed"])
        for position, tissue in enumerate(unit_order):
            if tissue not in completed:
                continue
            print(f"[grouped-pathway] {scope}:{tissue}", flush=True)
            importance, metrics = _run_unit(
                scope=scope,
                tissue=tissue,
                data=data,
                groups=groups,
                membership=membership,
                repeats=int(workflow["analysis"]["repeats"]),
                tissue_seed=base_seed + position * 100_000,
                outer_fraction=float(workflow["analysis"]["outer_fraction"]),
                inner_fraction=float(workflow["analysis"]["inner_fraction"]),
                regularization_grid=regularization_grid,
                synthetic_weight=synthetic_weight,
                permutation_repeats=permutation_repeats,
                permutation_seed=permutation_seed,
            )
            associations = _pathway_associations(
                scope=scope,
                tissue=tissue,
                data=data,
                groups=groups,
                membership=membership,
            )
            importance_tables.append(importance)
            metric_tables.append(metrics)
            association_tables.append(associations)
            completed_units.append((scope, tissue))
    if not importance_tables:
        raise RuntimeError("No grouped-pathway source was selected")

    rows = pd.concat(importance_tables, ignore_index=True)
    metrics = pd.concat(metric_tables, ignore_index=True)
    associations = pd.concat(association_tables, ignore_index=True)
    completed = (
        metrics.groupby(["scope", "tissue"], observed=True)["repeat"]
        .nunique()
        .rename("completed_repeats")
        .reset_index()
    )
    aggregate = _aggregate_pathways(rows, completed)
    _, utility = _metric_summary(metrics)
    comparison = _compare_arms(
        aggregate,
        utility,
        associations,
        minimum_importance=minimum_importance,
        minimum_fraction=minimum_fraction,
    )
    eligible = comparison.loc[comparison["eligible_synthetic_pathway"]].copy()
    nonredundant = _nonredundant_top_pathways(
        comparison,
        maximum_per_arm=int(options["top_nonredundant_per_arm"]),
        maximum_jaccard=float(options["maximum_gene_jaccard"]),
    )

    rows.to_csv(
        output / "pathway_importance_by_repeat.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    metrics.to_csv(output / "classifier_metrics.tsv.gz", sep="\t", index=False)
    utility.to_csv(output / "arm_utility.tsv", sep="\t", index=False)
    aggregate.to_csv(
        output / "pathway_importance_summary.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    associations.to_csv(
        output / "real_pathway_random_effects.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    comparison.to_csv(
        output / "pathway_arm_comparison.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    eligible.to_csv(
        output / "eligible_synthetic_pathways.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    nonredundant.to_csv(
        output / "top_nonredundant_pathways.tsv", sep="\t", index=False
    )
    for scope, tissue in completed_units:
        _plot_unit(
            aggregate,
            output,
            scope=scope,
            tissue=tissue,
            top_count=int(options["top_pathways_per_plot"]),
        )

    summary = {
        "status": "complete",
        "config": str(config_path.resolve()),
        "matched_config": str(matched_config_path.resolve()),
        "matched_output": str(Path(config["matched_output"]).resolve()),
        "output": str(output.resolve()),
        "completed_units": len(completed_units),
        "scopes": completed["scope"].value_counts().to_dict(),
        "repeats": int(completed["completed_repeats"].max()),
        "genes": gene_count,
        "pathways": pathway_count,
        "permutation_repeats": permutation_repeats,
        "minimum_permutation_roc_auc": minimum_importance,
        "minimum_positive_outer_fraction": minimum_fraction,
        "eligible_pathway_arm_rows": int(len(eligible)),
        "eligible_unique_pathways": int(
            eligible[["scope", "tissue", "term"]].drop_duplicates().shape[0]
        ),
        "eligible_units": sorted(
            (eligible["scope"] + ":" + eligible["tissue"]).unique().tolist()
        ),
        "eligible_by_unit": {
            f"{scope}:{tissue}": int(len(frame))
            for (scope, tissue), frame in eligible.groupby(
                ["scope", "tissue"], sort=True, observed=True
            )
        },
        "nonredundant_rows": int(len(nonredundant)),
        "ddim_retrained": False,
        "association_profiles": "observed OSDR only",
        "limitations": [
            "Reactome pathways overlap and are not independent features.",
            "Larger pathways can remove more predictive information when permuted.",
            "Grouped SHAP attributes fitted-model logit contribution, not causality.",
            "Repeated classifier splits are not independent biological studies.",
        ],
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    _plot_eligible_summary(nonredundant, output)
    _write_readme(output, summary, nonredundant)
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
