"""Aggregate cross-fitted whole-study transfer across eligible OSDR tissues."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from nasa_mouse_generative.effect_validation import random_effects_table

from .generated_feature_guidance import (
    METRICS,
    _aggregate,
    _evaluate_fold,
    _symbol_mapping,
)


def _safe_correlation(first: Iterable[float], second: Iterable[float]) -> float:
    first_values = np.asarray(list(first), dtype=float)
    second_values = np.asarray(list(second), dtype=float)
    finite = np.isfinite(first_values) & np.isfinite(second_values)
    if (
        finite.sum() < 2
        or np.std(first_values[finite]) == 0
        or np.std(second_values[finite]) == 0
    ):
        return float("nan")
    return float(np.corrcoef(first_values[finite], second_values[finite])[0, 1])


def _effect_recovery(
    real_effects: pd.DataFrame,
    synthetic_effects: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Summarize cross-fitted gene effects by accession and tissue."""

    keys = ["fold", "accession", "tissue", "feature"]
    real = real_effects.rename(
        columns={
            "flight_minus_ground": "real_effect",
            "effect_variance": "real_variance",
        }
    )
    synthetic = synthetic_effects.rename(
        columns={
            "flight_minus_ground": "synthetic_effect",
            "effect_variance": "synthetic_variance",
        }
    )
    comparison = real[
        keys + ["n_flight", "n_ground_control", "real_effect", "real_variance"]
    ].merge(
        synthetic[keys + ["synthetic_effect", "synthetic_variance"]],
        on=keys,
        how="inner",
        validate="one_to_one",
    )
    if comparison.empty:
        raise ValueError("No aligned real and synthetic accession effects")
    duplicated = comparison.duplicated(["accession", "tissue", "feature"])
    if duplicated.any():
        examples = comparison.loc[
            duplicated, ["fold", "accession", "tissue", "feature"]
        ].head()
        raise ValueError(
            "A tissue/accession was evaluated in more than one outer fold: "
            f"{examples.to_dict(orient='records')}"
        )
    comparison["direction_agrees"] = np.sign(comparison["real_effect"]).eq(
        np.sign(comparison["synthetic_effect"])
    )
    comparison["squared_error"] = np.square(
        comparison["real_effect"] - comparison["synthetic_effect"]
    )

    accession_rows: list[dict[str, object]] = []
    for (tissue, accession, fold), frame in comparison.groupby(
        ["tissue", "accession", "fold"], sort=True, observed=True
    ):
        accession_rows.append(
            {
                "tissue": str(tissue),
                "accession": str(accession),
                "fold": str(fold),
                "features": int(len(frame)),
                "effect_correlation": _safe_correlation(
                    frame["real_effect"], frame["synthetic_effect"]
                ),
                "direction_agreement": float(frame["direction_agrees"].mean()),
                "effect_rmse": float(np.sqrt(frame["squared_error"].mean())),
            }
        )
    accession_summary = pd.DataFrame(accession_rows)

    tissue_rows: list[dict[str, object]] = []
    meta_rows: list[pd.DataFrame] = []
    for tissue, frame in comparison.groupby("tissue", sort=True, observed=True):
        real_frame = frame[
            ["accession", "feature", "real_effect", "real_variance"]
        ].rename(
            columns={
                "real_effect": "flight_minus_ground",
                "real_variance": "effect_variance",
            }
        )
        synthetic_frame = frame[
            ["accession", "feature", "synthetic_effect", "synthetic_variance"]
        ].rename(
            columns={
                "synthetic_effect": "flight_minus_ground",
                "synthetic_variance": "effect_variance",
            }
        )
        real_meta = random_effects_table(real_frame).add_prefix("real_").rename(
            columns={"real_feature": "feature"}
        )
        synthetic_meta = random_effects_table(synthetic_frame).add_prefix(
            "synthetic_"
        ).rename(columns={"synthetic_feature": "feature"})
        meta = real_meta.merge(
            synthetic_meta, on="feature", how="inner", validate="one_to_one"
        )
        meta.insert(0, "tissue", str(tissue))
        meta["direction_agrees"] = np.sign(meta["real_meta_effect"]).eq(
            np.sign(meta["synthetic_meta_effect"])
        )
        meta_rows.append(meta)
        accession_subset = accession_summary.loc[
            accession_summary["tissue"].eq(str(tissue))
        ]
        tissue_rows.append(
            {
                "tissue": str(tissue),
                "accessions": int(frame["accession"].nunique()),
                "features": int(len(meta)),
                "meta_effect_correlation": _safe_correlation(
                    meta["real_meta_effect"], meta["synthetic_meta_effect"]
                ),
                "meta_direction_agreement": float(meta["direction_agrees"].mean()),
                "mean_accession_effect_correlation": float(
                    accession_subset["effect_correlation"].mean()
                ),
                "mean_accession_direction_agreement": float(
                    accession_subset["direction_agreement"].mean()
                ),
            }
        )
    return comparison, accession_summary, pd.DataFrame(tissue_rows), pd.concat(
        meta_rows, ignore_index=True
    )


def _plot_transfer_metrics(table: pd.DataFrame, output: Path) -> None:
    ordered = table.sort_values("delta_balanced_accuracy")
    positions = np.arange(len(ordered))
    figure, axes = plt.subplots(1, 3, figsize=(13.4, 6.8), sharey=True)
    labels = {
        "balanced_accuracy": "Balanced accuracy",
        "roc_auc": "AUROC",
        "average_precision": "Average precision",
    }
    for axis, metric in zip(axes, METRICS):
        values = ordered[f"delta_{metric}"].to_numpy(dtype=float)
        axis.axvline(0, color="#6B7280", linewidth=1)
        axis.scatter(
            values,
            positions,
            color=np.where(values >= 0, "#16847F", "#D96652"),
            s=46,
            zorder=3,
        )
        axis.set_title(labels[metric], fontweight="bold")
        axis.set_xlabel("Deployed minus real-only")
        axis.grid(axis="x", alpha=0.18)
    axes[0].set_yticks(positions, ordered["tissue"].str.replace("_", " "))
    figure.suptitle(
        "Whole-study transfer on real held-out OSDR accessions",
        fontweight="bold",
    )
    figure.tight_layout()
    figure.savefig(output / "whole_study_transfer_metric_deltas.png", dpi=240)
    figure.savefig(output / "whole_study_transfer_metric_deltas.pdf")
    plt.close(figure)


def _plot_effect_levels(
    pooled: pd.DataFrame,
    tissues: pd.DataFrame,
    accessions: pd.DataFrame,
    output: Path,
) -> None:
    ordered = tissues.sort_values("meta_effect_correlation")
    positions = np.arange(len(ordered))
    figure, axes = plt.subplots(1, 3, figsize=(14.2, 6.7))
    axes[0].scatter(
        pooled["delta_correlation"],
        pooled["direction_agreement"],
        color="#3D6FA3",
        s=58,
    )
    label_offsets = {
        "fold0": (-29, -12),
        "fold1": (5, -12),
        "fold2": (5, 6),
    }
    for row in pooled.itertuples(index=False):
        axes[0].annotate(
            str(row.fold),
            (row.delta_correlation, row.direction_agreement),
            xytext=label_offsets.get(str(row.fold), (4, 4)),
            textcoords="offset points",
            fontsize=8,
        )
    axes[0].set_xlabel("Effect correlation")
    axes[0].set_ylabel("Direction agreement")
    axes[0].set_title("Pooled across tissues", fontweight="bold")

    axes[1].axvline(0, color="#6B7280", linewidth=1)
    axes[1].scatter(
        ordered["meta_effect_correlation"],
        positions,
        color="#16847F",
        s=46,
    )
    axes[1].set_yticks(positions, ordered["tissue"].str.replace("_", " "))
    axes[1].set_xlabel("Cross-fitted meta-effect correlation")
    axes[1].set_title("Per tissue", fontweight="bold")

    tissue_order = ordered["tissue"].tolist()
    lookup = {tissue: index for index, tissue in enumerate(tissue_order)}
    jitter = np.linspace(-0.12, 0.12, max(1, len(accessions)))
    for index, row in enumerate(accessions.itertuples(index=False)):
        axes[2].scatter(
            row.effect_correlation,
            lookup[str(row.tissue)] + jitter[index % len(jitter)],
            color="#D69A24",
            s=24,
            alpha=0.78,
        )
    axes[2].axvline(0, color="#6B7280", linewidth=1)
    axes[2].set_yticks(positions, ordered["tissue"].str.replace("_", " "))
    axes[2].set_xlabel("Gene-effect correlation")
    axes[2].set_title("Per tissue and accession", fontweight="bold")
    for axis in axes:
        axis.grid(alpha=0.16)
    figure.suptitle(
        "Conditional FLT/GC recovery answers different aggregation questions",
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(output / "effect_recovery_levels.png", dpi=240)
    figure.savefig(output / "effect_recovery_levels.pdf")
    plt.close(figure)


def _read_fold_effects(
    folds: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pooled_rows: list[dict[str, object]] = []
    real_tables: list[pd.DataFrame] = []
    synthetic_tables: list[pd.DataFrame] = []
    for fold in folds:
        label = str(fold["label"])
        evaluation = Path(str(fold["evaluation_dir"]))
        summary = json.loads((evaluation / "summary.json").read_text(encoding="utf-8"))
        pooled_rows.append({"fold": label, **summary["flt_gc_effect_recovery"]})
        for name, destination in (
            ("gene_real_per_accession.tsv.gz", real_tables),
            ("gene_synthetic_per_accession.tsv.gz", synthetic_tables),
        ):
            table = pd.read_csv(evaluation / "accession_validation" / name, sep="\t")
            table.insert(0, "fold", label)
            destination.append(table)
    return (
        pd.DataFrame(pooled_rows),
        pd.concat(real_tables, ignore_index=True),
        pd.concat(synthetic_tables, ignore_index=True),
    )


def _write_readme(
    output: Path,
    tissue_results: pd.DataFrame,
    pooled_effects: pd.DataFrame,
    tissue_effects: pd.DataFrame,
    summary: dict[str, object],
) -> None:
    table = tissue_results[
        [
            "tissue",
            "accessions",
            "baseline_balanced_accuracy",
            "deployed_balanced_accuracy",
            "delta_balanced_accuracy",
            "delta_roc_auc",
            "delta_average_precision",
            "inner_gate_passed_folds",
        ]
    ].copy()
    effect = tissue_effects[
        [
            "tissue",
            "meta_effect_correlation",
            "meta_direction_agreement",
            "mean_accession_effect_correlation",
        ]
    ]
    table = table.merge(effect, on="tissue", how="left")
    readme = f"""# Twelve-tissue whole-study transfer

Three cross-study outer folds were fixed before evaluation. Every eligible OSDR
accession from the twelve tissues was assigned to test exactly once. For each
fold, the held-out accessions were absent from OSDR diffusion adaptation,
feature-policy selection, and classifier fitting. The ARCHS4 backbone excluded
all GEO series linked to eligible OSDR accessions.

This remains retrospective because the repository had already been used to
inspect related OSDR outcomes. It is nevertheless a uniform cross-fitted
comparison and replaces the earlier mixture of tissue-specific holdout recipes.

## Classification transfer

{_markdown_table(table)}

Global accession-macro deployed-minus-real-only changes were:

- balanced accuracy: {summary['accession_macro']['delta']['balanced_accuracy']:.3f}
- AUROC: {summary['accession_macro']['delta']['roc_auc']:.3f}
- average precision: {summary['accession_macro']['delta']['average_precision']:.3f}

## Three effect-recovery levels

1. **Pooled FLT/GC recovery** compares one real and synthetic 974-gene effect
   vector after combining all held-out tissues in a fold. It is a broad
   condition-label sanity check and can be dominated by tissue composition.
2. **Per-tissue recovery** compares cross-fitted random-effects FLT/GC vectors
   separately within each tissue. Each tissue contributes independently.
3. **Per-tissue, per-accession recovery** compares real and synthetic FLT/GC
   vectors inside each held-out study. This is the most local test and exposes
   study heterogeneity that pooled recovery can hide.

Pooled fold results:

{_markdown_table(pooled_effects)}

Generated profiles are never counted as biological replicates. Biological
association P values and BH-FDR remain calculated from real OSDR samples only.
"""
    (output / "README.md").write_text(readme, encoding="utf-8")


def _markdown_table(frame: pd.DataFrame) -> str:
    """Render a small result table without pandas' optional tabulate dependency."""

    def format_cell(value: object) -> str:
        if pd.isna(value):
            rendered = "NA"
        elif isinstance(value, (float, np.floating)):
            rendered = f"{float(value):.3f}"
        else:
            rendered = str(value)
        return rendered.replace("|", "\\|").replace("\n", " ")

    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(format_cell(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, divider, *rows])


def run(config_path: Path, *, seed: int = 6200) -> dict[str, object]:
    source = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    output = Path(str(source["output_root"]))
    output.mkdir(parents=True, exist_ok=True)
    folds = list(source["folds"])
    tissues = list(map(str, source["tissues"]))
    if len(folds) < 2 or len(tissues) != len(set(tissues)):
        raise ValueError("Whole-study transfer requires multiple folds and unique tissues")
    symbols = _symbol_mapping(Path(str(source["annotations"]["archs4_h5"])))

    tissue_rows: list[dict[str, object]] = []
    accession_tables: list[pd.DataFrame] = []
    prediction_tables: list[pd.DataFrame] = []
    for tissue_index, tissue in enumerate(tissues):
        print(f"[whole-study-transfer] evaluating {tissue}", flush=True)
        tissue_output = output / tissue
        config = {
            "tissue": tissue,
            "output_dir": str(tissue_output),
            "annotations": source["annotations"],
            "grid": source["grid"],
            "folds": folds,
        }
        for fold_index, fold in enumerate(folds):
            _evaluate_fold(
                fold,
                config,
                symbols,
                seed=seed + tissue_index * 100 + fold_index,
            )
        tissue_summary = _aggregate(config, symbols)
        row: dict[str, object] = {
            "tissue": tissue,
            "outer_folds": tissue_summary["outer_folds"],
            "accessions": tissue_summary["outer_accessions"],
            "profiles": tissue_summary["outer_profiles"],
            "inner_gate_passed_folds": tissue_summary["inner_gate_passed_folds"],
            "promising_by_predeclared_rule": tissue_summary[
                "promising_by_predeclared_rule"
            ],
        }
        for metric in METRICS:
            row[f"baseline_{metric}"] = tissue_summary["accession_macro"][
                "baseline"
            ][metric]
            row[f"deployed_{metric}"] = tissue_summary["accession_macro"][
                "deployed"
            ][metric]
            row[f"delta_{metric}"] = tissue_summary["deployed_minus_baseline"][
                metric
            ]
        tissue_rows.append(row)
        accessions = pd.read_csv(
            tissue_output / "outer_accession_results.tsv", sep="\t"
        )
        accessions.insert(0, "tissue", tissue)
        accession_tables.append(accessions)
        predictions = pd.read_csv(
            tissue_output / "outer_predictions.tsv.gz", sep="\t"
        )
        predictions.insert(0, "analysis_tissue", tissue)
        prediction_tables.append(predictions)

    tissue_results = pd.DataFrame(tissue_rows)
    accession_results = pd.concat(accession_tables, ignore_index=True)
    predictions = pd.concat(prediction_tables, ignore_index=True)
    tissue_results.to_csv(output / "tissue_results.tsv", sep="\t", index=False)
    accession_results.to_csv(
        output / "accession_results.tsv", sep="\t", index=False
    )
    predictions.to_csv(output / "predictions.tsv.gz", sep="\t", index=False)

    baseline_keys = accession_results.loc[
        accession_results["arm"].eq("baseline"), ["tissue", "accession"]
    ]
    if baseline_keys.duplicated().any():
        raise ValueError("At least one tissue/accession was tested in multiple folds")
    macro: dict[str, dict[str, float]] = {}
    for arm in ("baseline", "generated", "deployed"):
        subset = accession_results.loc[accession_results["arm"].eq(arm)]
        macro[arm] = {metric: float(subset[metric].mean()) for metric in METRICS}
    delta = {
        metric: macro["deployed"][metric] - macro["baseline"][metric]
        for metric in METRICS
    }

    pooled_effects, real_effects, synthetic_effects = _read_fold_effects(folds)
    (
        effect_comparison,
        accession_effects,
        tissue_effects,
        tissue_meta_effects,
    ) = _effect_recovery(real_effects, synthetic_effects)
    pooled_effects.to_csv(output / "pooled_effect_recovery.tsv", sep="\t", index=False)
    effect_comparison.to_csv(
        output / "cross_fitted_accession_gene_effects.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    accession_effects.to_csv(
        output / "accession_effect_recovery.tsv", sep="\t", index=False
    )
    tissue_effects.to_csv(
        output / "tissue_effect_recovery.tsv", sep="\t", index=False
    )
    tissue_meta_effects.to_csv(
        output / "tissue_meta_effects.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    summary = {
        "status": "complete",
        "design": "three_fold_cross_fitted_whole_accession_transfer",
        "tissues": tissues,
        "outer_folds": len(folds),
        "test_tissue_accession_pairs": int(len(baseline_keys)),
        "test_accessions": int(baseline_keys["accession"].nunique()),
        "test_profiles": int(len(predictions)),
        "accession_macro": {"baseline": macro["baseline"], "deployed": macro["deployed"], "delta": delta},
        "tissues_passing_predeclared_rule": int(
            tissue_results["promising_by_predeclared_rule"].sum()
        ),
        "pooled_effect_recovery": {
            "mean_correlation": float(pooled_effects["delta_correlation"].mean()),
            "mean_direction_agreement": float(
                pooled_effects["direction_agreement"].mean()
            ),
        },
        "macro_tissue_effect_recovery": {
            "mean_meta_effect_correlation": float(
                tissue_effects["meta_effect_correlation"].mean()
            ),
            "mean_meta_direction_agreement": float(
                tissue_effects["meta_direction_agreement"].mean()
            ),
            "mean_accession_effect_correlation": float(
                accession_effects["effect_correlation"].mean()
            ),
            "mean_accession_direction_agreement": float(
                accession_effects["direction_agreement"].mean()
            ),
        },
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    _plot_transfer_metrics(tissue_results, output)
    _plot_effect_levels(pooled_effects, tissue_effects, accession_effects, output)
    _write_readme(output, tissue_results, pooled_effects, tissue_effects, summary)
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--seed", default=6200, type=int)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    run(arguments.config, seed=arguments.seed)


if __name__ == "__main__":
    main()
