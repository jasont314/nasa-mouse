"""Muscle-composition QC for aggregate OSDR liver runs."""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .aggregate_liver_mober import load_gene_symbols
from .osd379_tissue_qc import (
    LIVER_MARKERS,
    MUSCLE_MARKERS,
    SEVERE_MIN_MARKER_FRACTION_PERCENT,
    SEVERE_MIN_MARKERS_OVER_100,
    robust_z,
)


DEFAULT_RUN_DIR = "outputs/mober_liver_ribo6_osdr"
DEFAULT_OSDR_H5 = "assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5"
DEFAULT_OUTPUT_DIR = "outputs/mober_liver_ribo6_osdr/muscle_outlier_qc"
MIN_STRATUM_GROUP_SIZE = 4
MIN_BROAD_GROUP_SIZE = 6


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def marker_indices(
    genes: np.ndarray,
    gene_symbols: dict[str, str],
    markers: list[str],
    require_all: bool = True,
) -> tuple[list[int], dict[str, str]]:
    gene_to_index = {str(gene): index for index, gene in enumerate(genes)}
    symbol_to_gene: dict[str, str] = {}
    for gene in genes:
        symbol = gene_symbols.get(str(gene), "")
        if symbol in markers and symbol not in symbol_to_gene:
            symbol_to_gene[symbol] = str(gene)
    missing = sorted(set(markers) - set(symbol_to_gene))
    if require_all and missing:
        raise ValueError(f"Missing marker genes from aggregate matrix: {missing}")
    indices = [
        gene_to_index[symbol_to_gene[symbol]]
        for symbol in markers
        if symbol in symbol_to_gene
    ]
    return indices, symbol_to_gene


def load_run_matrix(run_dir: Path) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    target = np.load(run_dir / "controlled_target.npz")
    genes = target["genes"].astype(str)
    flt = target["flt"].astype(float, copy=False)
    gc = target["gc"].astype(float, copy=False)
    flt_features = target["flt_features"].astype(str)
    gc_features = target["gc_features"].astype(str)
    expression = np.concatenate([flt, gc], axis=1)
    features = list(flt_features) + list(gc_features)
    metadata = pd.read_csv(run_dir / "retained_profile_features.tsv", sep="\t")
    if len(metadata) != len(features):
        raise ValueError(
            "Metadata rows do not match target matrix samples: "
            f"{len(metadata)} vs {len(features)}"
        )
    metadata = metadata.copy()
    metadata.insert(0, "sample", features)
    return expression, genes, metadata


def infer_collection(row: pd.Series) -> str:
    text = " ".join(
        str(row.get(column, ""))
        for column in ["profile", "source_name", "official_source_name", "h5_sample_name"]
    ).upper()
    if "ISS-T" in text or "ISS_T" in text or "ISS TERMINAL" in text:
        return "ISS-T"
    if "LAR" in text or "LIVE ANIMAL RETURN" in text:
        return "LAR"
    return "unknown"


def infer_age_label(row: pd.Series) -> str:
    profile = str(row.get("profile", "")).upper()
    if re.search(r"(^|_)OLD($|_)", profile):
        return "OLD"
    if re.search(r"(^|_)YNG($|_)", profile):
        return "YNG"
    text = " ".join(
        str(row.get(column, ""))
        for column in ["age", "age_at_launch"]
        if str(row.get(column, "")).strip() and str(row.get(column, "")).lower() != "nan"
    ).upper()
    cleaned = re.sub(r"\s+", " ", text).strip()
    return cleaned or "unknown"


def add_inferred_strata(scores: pd.DataFrame) -> pd.DataFrame:
    scores = scores.copy()
    scores["collection_inferred"] = scores.apply(infer_collection, axis=1)
    scores["age_label_inferred"] = scores.apply(infer_age_label, axis=1)
    return scores


def add_group_zscores(scores: pd.DataFrame) -> pd.DataFrame:
    scores = scores.copy()
    scores = add_inferred_strata(scores)
    scores["accession_group_size"] = scores.groupby("h5_accession")[
        "sample"
    ].transform("size")
    scores["accession_condition_group_size"] = scores.groupby(
        ["h5_accession", "location"], dropna=False
    )["sample"].transform("size")
    scores["muscle_mean_log2_robust_z_within_accession"] = scores.groupby(
        "h5_accession", dropna=False
    )["muscle_mean_log2"].transform(robust_z)
    scores["muscle_mean_log2_robust_z_within_accession_condition"] = scores.groupby(
        ["h5_accession", "location"], dropna=False
    )["muscle_mean_log2"].transform(robust_z)
    stratum_columns = [
        "h5_accession",
        "location",
        "collection_inferred",
        "age_label_inferred",
    ]
    scores["stratum_group_size"] = scores.groupby(stratum_columns, dropna=False)[
        "sample"
    ].transform("size")
    scores["muscle_mean_log2_robust_z_within_stratum"] = scores.groupby(
        stratum_columns, dropna=False
    )["muscle_mean_log2"].transform(robust_z)
    return scores


def score_samples(
    expression: np.ndarray,
    genes: np.ndarray,
    metadata: pd.DataFrame,
    gene_symbols: dict[str, str],
) -> tuple[pd.DataFrame, dict]:
    muscle_indices, muscle_ids = marker_indices(
        genes, gene_symbols, MUSCLE_MARKERS, require_all=True
    )
    liver_indices, liver_ids = marker_indices(
        genes, gene_symbols, LIVER_MARKERS, require_all=False
    )
    library_sizes = expression.sum(axis=0)
    library_sizes[library_sizes <= 0] = np.nan
    logged = np.log2(expression + 1.0)
    muscle = expression[muscle_indices, :]
    liver = expression[liver_indices, :]

    scores = metadata.copy()
    scores["muscle_mean_log2"] = logged[muscle_indices, :].mean(axis=0)
    scores["muscle_median_log2"] = np.median(logged[muscle_indices, :], axis=0)
    scores["muscle_markers_over_100"] = (muscle > 100).sum(axis=0)
    scores["muscle_marker_fraction_percent"] = 100.0 * muscle.sum(axis=0) / library_sizes
    scores["liver_mean_log2"] = logged[liver_indices, :].mean(axis=0)
    scores["liver_markers_over_100"] = (liver > 100).sum(axis=0)
    scores["muscle_liver_mean_log2_delta"] = (
        scores["muscle_mean_log2"] - scores["liver_mean_log2"]
    )
    scores["high_muscle_abundance"] = scores["muscle_markers_over_100"].ge(
        SEVERE_MIN_MARKERS_OVER_100
    ) & scores["muscle_marker_fraction_percent"].ge(
        SEVERE_MIN_MARKER_FRACTION_PERCENT
    )
    scores = add_group_zscores(scores)
    scores["relative_stratum_outlier"] = (
        scores["stratum_group_size"].ge(MIN_STRATUM_GROUP_SIZE)
        & scores["muscle_mean_log2_robust_z_within_stratum"].ge(3.5)
    )
    scores["relative_accession_outlier"] = (
        scores["accession_group_size"].ge(MIN_BROAD_GROUP_SIZE)
        & scores["muscle_mean_log2_robust_z_within_accession"].ge(3.5)
    )
    scores["relative_accession_condition_outlier"] = (
        scores["accession_condition_group_size"].ge(MIN_BROAD_GROUP_SIZE)
        & scores["muscle_mean_log2_robust_z_within_accession_condition"].ge(3.5)
    )
    scores["candidate_muscle_outlier"] = (
        scores["high_muscle_abundance"] & scores["relative_stratum_outlier"]
    )
    scores["candidate_small_group_warning"] = (
        scores["candidate_muscle_outlier"]
        & scores["stratum_group_size"].lt(MIN_BROAD_GROUP_SIZE)
    )
    scores["broad_review_muscle_outlier"] = (
        scores["high_muscle_abundance"]
        & ~scores["candidate_muscle_outlier"]
        & (
            scores["relative_accession_outlier"]
            | scores["relative_accession_condition_outlier"]
        )
    )
    scores["review_high_abundance_not_relative"] = (
        scores["high_muscle_abundance"]
        & ~scores["candidate_muscle_outlier"]
        & ~scores["broad_review_muscle_outlier"]
    )

    summary = {
        "muscle_markers": muscle_ids,
        "liver_markers": liver_ids,
        "rules": {
            "high_muscle_abundance": {
                "muscle_markers_over_100_min": SEVERE_MIN_MARKERS_OVER_100,
                "muscle_marker_fraction_percent_min": SEVERE_MIN_MARKER_FRACTION_PERCENT,
            },
            "relative_outlier": {
                "robust_z_min": 3.5,
                "minimum_stratum_group_size": MIN_STRATUM_GROUP_SIZE,
                "minimum_broad_group_size": MIN_BROAD_GROUP_SIZE,
                "candidate_group": (
                    "h5_accession + location + collection_inferred + age_label_inferred"
                ),
                "broad_review_groups": ["h5_accession", "h5_accession + location"],
            },
        },
    }
    return scores.sort_values("muscle_marker_fraction_percent", ascending=False), summary


def accession_condition_summary(scores: pd.DataFrame) -> pd.DataFrame:
    return (
        scores.groupby(["h5_accession", "location"], dropna=False)
        .agg(
            samples=("sample", "count"),
            median_muscle_fraction_percent=("muscle_marker_fraction_percent", "median"),
            max_muscle_fraction_percent=("muscle_marker_fraction_percent", "max"),
            median_muscle_mean_log2=("muscle_mean_log2", "median"),
            max_muscle_mean_log2=("muscle_mean_log2", "max"),
            high_muscle_abundance=("high_muscle_abundance", "sum"),
            candidate_muscle_outlier=("candidate_muscle_outlier", "sum"),
            broad_review_muscle_outlier=("broad_review_muscle_outlier", "sum"),
        )
        .reset_index()
        .sort_values(
            ["candidate_muscle_outlier", "high_muscle_abundance", "max_muscle_fraction_percent"],
            ascending=[False, False, False],
        )
    )


def write_report(output_dir: Path, scores: pd.DataFrame, summary: dict) -> None:
    candidates = scores.loc[scores["candidate_muscle_outlier"]]
    broad_review = scores.loc[scores["broad_review_muscle_outlier"]]
    high_only = scores.loc[scores["review_high_abundance_not_relative"]]
    top_columns = [
        "h5_accession",
        "location",
        "collection_inferred",
        "age_label_inferred",
        "profile",
        "sample",
        "muscle_marker_fraction_percent",
        "muscle_markers_over_100",
        "muscle_mean_log2",
        "liver_mean_log2",
        "stratum_group_size",
        "muscle_mean_log2_robust_z_within_stratum",
        "muscle_mean_log2_robust_z_within_accession",
        "muscle_mean_log2_robust_z_within_accession_condition",
        "candidate_muscle_outlier",
        "candidate_small_group_warning",
        "broad_review_muscle_outlier",
        "review_high_abundance_not_relative",
    ]
    candidates_lines = (
        candidates[top_columns].to_csv(sep="\t", index=False).strip().splitlines()
        if len(candidates)
        else ["none"]
    )
    broad_review_lines = (
        broad_review[top_columns].to_csv(sep="\t", index=False).strip().splitlines()
        if len(broad_review)
        else ["none"]
    )
    high_lines = (
        high_only[top_columns].head(25).to_csv(sep="\t", index=False).strip().splitlines()
        if len(high_only)
        else ["none"]
    )
    accession_lines = (
        pd.DataFrame(summary["accession_condition_summary"])
        .to_csv(sep="\t", index=False)
        .strip()
        .splitlines()
    )
    report = [
        "# Aggregate Liver Muscle Outlier QC",
        "",
        f"- Samples scored: {summary['samples_scored']}",
        f"- Muscle markers found: {len(summary['muscle_markers'])}",
        f"- High muscle-abundance samples: {summary['high_muscle_abundance_samples']}",
        f"- Candidate muscle outliers: {summary['candidate_muscle_outliers']}",
        f"- Broad-review muscle outliers: {summary['broad_review_muscle_outliers']}",
        "",
        "A candidate requires high muscle-marker abundance and relative elevation",
        "within the finest inferred accession/condition/collection/age stratum.",
        "Broad-review rows meet the abundance rule and are accession-level",
        "relative outliers, but not stratum-level candidates.",
        "",
        "## Candidate Outliers",
        "",
        "```tsv",
        *candidates_lines,
        "```",
        "",
        "## Broad-Review Outliers",
        "",
        "```tsv",
        *broad_review_lines,
        "```",
        "",
        "## High Abundance But Not Relative Outliers",
        "",
        "```tsv",
        *high_lines,
        "```",
        "",
        "## Accession/Condition Summary",
        "",
        "```tsv",
        *accession_lines,
        "```",
    ]
    (output_dir / "AGGREGATE_LIVER_MUSCLE_QC.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score aggregate OSDR liver samples for muscle-marker outliers."
    )
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR)
    parser.add_argument("--osdr-h5", default=DEFAULT_OSDR_H5)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log(f"Loading aggregate run from {run_dir}")
    expression, genes, metadata = load_run_matrix(run_dir)
    symbols = load_gene_symbols(args.osdr_h5)
    log("Scoring muscle and liver marker panels")
    scores, summary = score_samples(expression, genes, metadata, symbols)
    grouped = accession_condition_summary(scores)

    scores.to_csv(output_dir / "sample_muscle_marker_scores.tsv", sep="\t", index=False)
    scores.loc[scores["candidate_muscle_outlier"]].to_csv(
        output_dir / "candidate_muscle_outliers.tsv", sep="\t", index=False
    )
    scores.loc[scores["broad_review_muscle_outlier"]].to_csv(
        output_dir / "broad_review_muscle_outliers.tsv", sep="\t", index=False
    )
    scores.loc[scores["high_muscle_abundance"]].to_csv(
        output_dir / "high_muscle_abundance_samples.tsv", sep="\t", index=False
    )
    grouped.to_csv(output_dir / "accession_condition_muscle_summary.tsv", sep="\t", index=False)
    suggested = scores.loc[scores["candidate_muscle_outlier"], "profile"].astype(str)
    (output_dir / "suggested_candidate_exclusion_profiles.txt").write_text(
        "\n".join(suggested.tolist()) + ("\n" if len(suggested) else ""),
        encoding="utf-8",
    )

    summary.update(
        {
            "run_dir": str(run_dir),
            "output_dir": str(output_dir),
            "samples_scored": int(len(scores)),
            "high_muscle_abundance_samples": int(scores["high_muscle_abundance"].sum()),
            "candidate_muscle_outliers": int(scores["candidate_muscle_outlier"].sum()),
            "broad_review_muscle_outliers": int(
                scores["broad_review_muscle_outlier"].sum()
            ),
            "candidate_profiles": suggested.tolist(),
            "accession_condition_summary": grouped.to_dict(orient="records"),
        }
    )
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_report(output_dir, scores, summary)
    log(f"Saved muscle outlier QC to {output_dir}")


if __name__ == "__main__":
    main()
