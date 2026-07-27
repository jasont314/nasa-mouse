"""Run a fixed DDIM-guided feature policy on fresh tissue test accessions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.stats import binomtest

from .generated_feature_guidance import (
    METRICS,
    _aggregate,
    _evaluate_fold,
    _symbol_mapping,
)


def _paired_accession_table(table: pd.DataFrame) -> pd.DataFrame:
    subset = table.loc[table["arm"].isin(["baseline", "generated"])]
    wide = subset.pivot(
        index=["tissue", "accession", "profiles"],
        columns="arm",
        values=list(METRICS),
    ).reset_index()
    columns = []
    for metric, arm in wide.columns:
        columns.append(metric if not arm else f"{arm}_{metric}")
    wide.columns = columns
    for metric in METRICS:
        wide[f"delta_{metric}"] = (
            wide[f"generated_{metric}"] - wide[f"baseline_{metric}"]
        )
    return wide


def _bootstrap_intervals(
    accessions: pd.DataFrame, *, seed: int, repeats: int = 50000
) -> dict[str, dict[str, float]]:
    rng = np.random.default_rng(seed)
    number = len(accessions)
    intervals: dict[str, dict[str, float]] = {}
    for metric in METRICS:
        delta = accessions[f"delta_{metric}"].to_numpy(dtype=float)
        draws = delta[rng.integers(0, number, size=(repeats, number))].mean(axis=1)
        intervals[metric] = {
            "low": float(np.quantile(draws, 0.025)),
            "high": float(np.quantile(draws, 0.975)),
        }
    return intervals


def _plot_tissues(table: pd.DataFrame, output: Path) -> None:
    positions = np.arange(len(table))
    width = 0.36
    figure, axis = plt.subplots(figsize=(10, 5.8))
    axis.bar(
        positions - width / 2,
        table["baseline_balanced_accuracy"],
        width,
        color="#4C78A8",
        label="Real-only baseline",
    )
    axis.bar(
        positions + width / 2,
        table["generated_balanced_accuracy"],
        width,
        color="#F58518",
        label="Generated-feature model",
    )
    axis.axhline(0.5, color="#333333", linestyle="--", linewidth=1)
    axis.set_xticks(positions, table["tissue"], rotation=25, ha="right")
    axis.set_ylim(0, 1.03)
    axis.set_ylabel("Fresh-test accession-macro balanced accuracy")
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output / "tissue_balanced_accuracy.png", dpi=220)
    figure.savefig(output / "tissue_balanced_accuracy.pdf")
    plt.close(figure)


def run(config_path: Path, *, seed: int = 5050) -> dict[str, object]:
    source = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    output = Path(str(source["output_root"]))
    output.mkdir(parents=True, exist_ok=True)
    symbols = _symbol_mapping(Path(str(source["annotations"]["archs4_h5"])))
    tissue_rows: list[dict[str, object]] = []
    accession_tables: list[pd.DataFrame] = []
    prediction_tables: list[pd.DataFrame] = []
    feature_tables: list[pd.DataFrame] = []

    for offset, tissue in enumerate(source["tissues"]):
        print(f"[generated-feature-transfer] evaluating {tissue}", flush=True)
        tissue_output = output / str(tissue)
        config = {
            "tissue": str(tissue),
            "output_dir": str(tissue_output),
            "annotations": source["annotations"],
            "grid": source["grid"],
            "folds": [source["fold"]],
        }
        _evaluate_fold(
            source["fold"],
            config,
            symbols,
            seed=seed + offset,
        )
        summary = _aggregate(config, symbols)
        baseline = summary["accession_macro"]["baseline"]
        generated = summary["accession_macro"]["generated"]
        row: dict[str, object] = {
            "tissue": str(tissue),
            "test_accessions": summary["outer_accessions"],
            "test_profiles": summary["outer_profiles"],
        }
        for metric in METRICS:
            row[f"baseline_{metric}"] = baseline[metric]
            row[f"generated_{metric}"] = generated[metric]
            row[f"delta_{metric}"] = generated[metric] - baseline[metric]
        row["improved_ba_without_auc_loss"] = bool(
            row["delta_balanced_accuracy"] > 0 and row["delta_roc_auc"] >= 0
        )
        tissue_rows.append(row)

        accessions = pd.read_csv(tissue_output / "outer_accession_results.tsv", sep="\t")
        accessions.insert(0, "tissue", str(tissue))
        accession_tables.append(accessions)
        predictions = pd.read_csv(tissue_output / "outer_predictions.tsv.gz", sep="\t")
        predictions["tissue"] = str(tissue)
        prediction_tables.append(predictions)
        features = pd.read_csv(tissue_output / "selected_features_all_folds.tsv.gz", sep="\t")
        features.insert(0, "tissue", str(tissue))
        feature_tables.append(features)

    tissues = pd.DataFrame(tissue_rows)
    accessions = pd.concat(accession_tables, ignore_index=True)
    predictions = pd.concat(prediction_tables, ignore_index=True)
    features = pd.concat(feature_tables, ignore_index=True)
    paired = _paired_accession_table(accessions)
    tissues.to_csv(output / "tissue_results.tsv", sep="\t", index=False)
    paired.to_csv(output / "accession_results.tsv", sep="\t", index=False)
    predictions.to_csv(output / "predictions.tsv.gz", sep="\t", index=False)
    features.to_csv(output / "selected_features.tsv.gz", sep="\t", index=False)

    cross_tissue_features = (
        features.groupby(["gene", "symbol"], observed=True)
        .agg(
            tissues_selected=("tissue", "nunique"),
            mean_classifier_coefficient=("classifier_coefficient", "mean"),
            mean_absolute_classifier_coefficient=(
                "classifier_coefficient",
                lambda values: float(np.mean(np.abs(values))),
            ),
            coefficient_sign_agreement=(
                "classifier_coefficient",
                lambda values: float(abs(np.mean(np.sign(values)))),
            ),
        )
        .reset_index()
        .sort_values(
            ["tissues_selected", "coefficient_sign_agreement", "mean_absolute_classifier_coefficient"],
            ascending=[False, False, False],
        )
    )
    cross_tissue_features.to_csv(
        output / "cross_tissue_feature_stability.tsv", sep="\t", index=False
    )

    baseline_macro = {
        metric: float(paired[f"baseline_{metric}"].mean()) for metric in METRICS
    }
    generated_macro = {
        metric: float(paired[f"generated_{metric}"].mean()) for metric in METRICS
    }
    delta = {
        metric: generated_macro[metric] - baseline_macro[metric] for metric in METRICS
    }
    real_correct = (predictions["baseline_probability"] >= 0.5) == predictions["label"]
    generated_correct = (predictions["generated_probability"] >= 0.5) == predictions["label"]
    gained = int((~real_correct & generated_correct).sum())
    lost = int((real_correct & ~generated_correct).sum())
    summary = {
        "status": "complete",
        "design": "one_time_fresh_tissue_transfer",
        "tissues": list(map(str, source["tissues"])),
        "test_accessions": int(len(paired)),
        "test_profiles": int(len(predictions)),
        "tissues_improved_ba_without_auc_loss": int(
            tissues["improved_ba_without_auc_loss"].sum()
        ),
        "accession_macro": {
            "baseline": baseline_macro,
            "generated": generated_macro,
            "delta": delta,
            "paired_bootstrap_95pct_ci": _bootstrap_intervals(
                paired, seed=seed
            ),
        },
        "classification_changes": {
            "baseline_correct": int(real_correct.sum()),
            "generated_correct": int(generated_correct.sum()),
            "wrong_to_correct": gained,
            "correct_to_wrong": lost,
            "net_correct": gained - lost,
            "exact_mcnemar_p": float(
                binomtest(gained, gained + lost, p=0.5).pvalue
            ),
        },
        "successful_transfer": bool(
            delta["balanced_accuracy"] > 0
            and delta["roc_auc"] >= 0
            and delta["average_precision"] >= 0
        ),
        "top_cross_tissue_features": cross_tissue_features.head(20).to_dict(
            orient="records"
        ),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    _plot_tissues(tissues, output)
    print(json.dumps(summary, indent=2))
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--seed", default=5050, type=int)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    run(arguments.config, seed=arguments.seed)


if __name__ == "__main__":
    main()
