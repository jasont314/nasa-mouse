"""Summarize signed DGEA direction inside validated GLARE modules."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from .multi_tissue_validation import term_module_genes


DEFAULT_ROOT = "outputs/glare_multi_tissue_api"
DEFAULT_VALIDATION_DIR = "outputs/glare_multi_tissue_api/validation_stack_terms15"


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def safe_float(value: object) -> float:
    try:
        value = float(value)
    except Exception:
        return math.nan
    return value if math.isfinite(value) else math.nan


def safe_median(series: pd.Series) -> float:
    values = numeric(series).dropna()
    if values.empty:
        return math.nan
    return safe_float(values.median())


def safe_mean(series: pd.Series) -> float:
    values = numeric(series).dropna()
    if values.empty:
        return math.nan
    return safe_float(values.mean())


def direction_label(
    flt_up_sig: int,
    gc_up_sig: int,
    tested_up: int,
    tested_down: int,
    median_log2fc_tested: float,
    min_sig_genes: int,
    dominant_fraction: float,
) -> str:
    sig_total = flt_up_sig + gc_up_sig
    if sig_total >= min_sig_genes:
        flt_frac = flt_up_sig / sig_total
        gc_frac = gc_up_sig / sig_total
        if flt_frac >= dominant_fraction:
            return "FLT-up DEG-supported"
        if gc_frac >= dominant_fraction:
            return "GC-up DEG-supported"
        if flt_up_sig and gc_up_sig:
            return "mixed/reorganized DEG"
        return "ambiguous DEG"

    tested_total = tested_up + tested_down
    if tested_total:
        flt_frac = tested_up / tested_total
        gc_frac = tested_down / tested_total
        if flt_frac >= dominant_fraction and median_log2fc_tested > 0:
            return "weak FLT-up trend"
        if gc_frac >= dominant_fraction and median_log2fc_tested < 0:
            return "weak GC-up trend"
    return "no clear DGEA direction"


def meta_direction_label(
    flt_up_sig_sum: int,
    gc_up_sig_sum: int,
    mean_flight_minus_ground: float,
    score_fdr: float,
    score_empirical_p: float,
    score_direction_consistency: float,
    min_sig_genes: int,
    dominant_fraction: float,
) -> str:
    sig_total = flt_up_sig_sum + gc_up_sig_sum
    if sig_total >= min_sig_genes:
        flt_frac = flt_up_sig_sum / sig_total
        gc_frac = gc_up_sig_sum / sig_total
        if flt_frac >= dominant_fraction:
            return "FLT-up DEG-supported"
        if gc_frac >= dominant_fraction:
            return "GC-up DEG-supported"
        if flt_up_sig_sum and gc_up_sig_sum:
            return "mixed/reorganized DEG"
        return "ambiguous DEG"

    score_supported = (
        math.isfinite(score_fdr)
        and score_fdr < 0.05
        and math.isfinite(score_empirical_p)
        and score_empirical_p <= 0.05
    )
    if score_supported:
        if score_direction_consistency >= 0.67 and mean_flight_minus_ground > 0:
            return "hidden module-score FLT-up"
        if score_direction_consistency >= 0.67 and mean_flight_minus_ground < 0:
            return "hidden module-score GC-up"
        return "module-score significant, direction inconsistent"
    return "ambiguous/no validated direction"


def summarize_candidate_direction(
    root: Path,
    candidates: pd.DataFrame,
    alpha: float,
    lfc_cutoff: float,
    min_sig_genes: int,
    dominant_fraction: float,
) -> pd.DataFrame:
    rows = []
    for candidate in candidates.itertuples(index=False):
        enrich_path = root / candidate.tissue / "per_study" / "glare_cluster_reactome_enrichment.tsv"
        if not enrich_path.exists():
            continue
        enrichment = pd.read_csv(enrich_path, sep="\t")
        genes_by_study = term_module_genes(enrichment, str(candidate.term))
        for accession, module_genes in sorted(genes_by_study.items()):
            gene_path = root / candidate.tissue / "dgea_comparison" / f"{accession}_gene_level_glare_dgea.tsv"
            if not gene_path.exists() or len(module_genes) < 5:
                continue
            gene_table = pd.read_csv(gene_path, sep="\t")
            required = {"gene_id", "log2FoldChange", "padj", "tested_dgea"}
            if not required.issubset(gene_table.columns):
                continue

            module = pd.DataFrame({"gene_id": sorted(module_genes)})
            merged = module.merge(gene_table, on="gene_id", how="left")
            tested_mask = bool_series(merged["tested_dgea"]) if "tested_dgea" in merged else pd.Series(False, index=merged.index)
            log2fc = numeric(merged["log2FoldChange"])
            padj = numeric(merged["padj"])
            tested = merged.loc[tested_mask].copy()
            tested_log2fc = numeric(tested["log2FoldChange"])
            tested_padj = numeric(tested["padj"])
            sig_mask = tested_padj.lt(alpha) & tested_log2fc.abs().ge(lfc_cutoff)
            sig = tested.loc[sig_mask].copy()
            sig_log2fc = numeric(sig["log2FoldChange"])

            flt_up_sig = int(sig_log2fc.gt(0).sum())
            gc_up_sig = int(sig_log2fc.lt(0).sum())
            tested_up = int(tested_log2fc.gt(0).sum())
            tested_down = int(tested_log2fc.lt(0).sum())
            sig_total = flt_up_sig + gc_up_sig
            tested_signed_total = tested_up + tested_down
            median_tested = safe_median(tested["log2FoldChange"])
            median_sig = safe_median(sig["log2FoldChange"])
            mean_tested = safe_mean(tested["log2FoldChange"])
            mean_sig = safe_mean(sig["log2FoldChange"])

            rows.append(
                {
                    "tissue": candidate.tissue,
                    "module_class": candidate.module_class,
                    "term": candidate.term,
                    "clean_term": candidate.clean_term,
                    "accession": accession,
                    "module_genes": int(len(module_genes)),
                    "genes_with_dgea_row": int(merged["padj"].notna().sum()),
                    "tested_dgea_genes": int(len(tested)),
                    "significant_abs_lfc_genes": int(sig_total),
                    "flt_up_sig_genes": flt_up_sig,
                    "gc_up_sig_genes": gc_up_sig,
                    "flt_up_sig_fraction": flt_up_sig / sig_total if sig_total else math.nan,
                    "gc_up_sig_fraction": gc_up_sig / sig_total if sig_total else math.nan,
                    "tested_flt_up_genes": tested_up,
                    "tested_gc_up_genes": tested_down,
                    "tested_flt_up_fraction": tested_up / tested_signed_total if tested_signed_total else math.nan,
                    "tested_gc_up_fraction": tested_down / tested_signed_total if tested_signed_total else math.nan,
                    "median_log2fc_tested": median_tested,
                    "mean_log2fc_tested": mean_tested,
                    "median_log2fc_sig": median_sig,
                    "mean_log2fc_sig": mean_sig,
                    "direction_label": direction_label(
                        flt_up_sig,
                        gc_up_sig,
                        tested_up,
                        tested_down,
                        median_tested,
                        min_sig_genes,
                        dominant_fraction,
                    ),
                }
            )
    return pd.DataFrame(rows)


def meta_summarize_direction(
    by_study: pd.DataFrame,
    module_score_meta: pd.DataFrame,
    min_sig_genes: int,
    dominant_fraction: float,
) -> pd.DataFrame:
    if by_study.empty:
        return pd.DataFrame()
    rows = []
    group_cols = ["tissue", "module_class", "term", "clean_term"]
    for keys, group in by_study.groupby(group_cols, sort=True):
        flt_up_sig_sum = int(group["flt_up_sig_genes"].sum())
        gc_up_sig_sum = int(group["gc_up_sig_genes"].sum())
        sig_total = flt_up_sig_sum + gc_up_sig_sum
        study_labels = group["direction_label"].astype(str)
        rows.append(
            {
                "tissue": keys[0],
                "module_class": keys[1],
                "term": keys[2],
                "clean_term": keys[3],
                "studies_with_module": int(group["accession"].nunique()),
                "module_genes_median": safe_median(group["module_genes"]),
                "tested_dgea_genes_median": safe_median(group["tested_dgea_genes"]),
                "significant_abs_lfc_genes_sum": sig_total,
                "flt_up_sig_genes_sum": flt_up_sig_sum,
                "gc_up_sig_genes_sum": gc_up_sig_sum,
                "flt_up_sig_fraction": flt_up_sig_sum / sig_total if sig_total else math.nan,
                "gc_up_sig_fraction": gc_up_sig_sum / sig_total if sig_total else math.nan,
                "median_log2fc_tested_median": safe_median(group["median_log2fc_tested"]),
                "median_log2fc_sig_median": safe_median(group["median_log2fc_sig"]),
                "studies_flt_up_deg_supported": int(study_labels.eq("FLT-up DEG-supported").sum()),
                "studies_gc_up_deg_supported": int(study_labels.eq("GC-up DEG-supported").sum()),
                "studies_mixed_deg": int(study_labels.eq("mixed/reorganized DEG").sum()),
                "studies_no_clear_dgea": int(study_labels.eq("no clear DGEA direction").sum()),
            }
        )

    result = pd.DataFrame(rows)
    if not module_score_meta.empty:
        score_cols = [
            "tissue",
            "module_class",
            "term",
            "clean_term",
            "combined_welch_fdr_bh",
            "mean_flight_minus_ground",
            "median_flight_minus_ground",
            "direction_consistency",
            "median_empirical_abs_p",
        ]
        result = result.merge(
            module_score_meta[[col for col in score_cols if col in module_score_meta.columns]],
            on=group_cols,
            how="left",
        )
    for col in ["combined_welch_fdr_bh", "mean_flight_minus_ground", "direction_consistency", "median_empirical_abs_p"]:
        if col not in result:
            result[col] = math.nan
    result["direction_call"] = [
        meta_direction_label(
            int(row.flt_up_sig_genes_sum),
            int(row.gc_up_sig_genes_sum),
            safe_float(row.mean_flight_minus_ground),
            safe_float(row.combined_welch_fdr_bh),
            safe_float(row.median_empirical_abs_p),
            safe_float(row.direction_consistency),
            min_sig_genes,
            dominant_fraction,
        )
        for row in result.itertuples(index=False)
    ]
    return result.sort_values(
        ["module_class", "direction_call", "combined_welch_fdr_bh", "median_empirical_abs_p"],
        na_position="last",
    )


def markdown_table(frame: pd.DataFrame, columns: list[str], max_rows: int = 30) -> list[str]:
    if frame.empty:
        return ["No rows."]
    display = frame[[col for col in columns if col in frame.columns]].head(max_rows).copy()
    headers = list(display.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in display.itertuples(index=False):
        values = []
        for value in row:
            if isinstance(value, float):
                values.append("" if math.isnan(value) else f"{value:.4g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def write_report(meta: pd.DataFrame, output_dir: Path) -> None:
    lines = [
        "# Signed Direction Summary For GLARE Candidate Modules",
        "",
        "This report adds within-module DGEA direction checks to the validated GLARE candidate modules.",
        "Positive log2FC means higher in spaceflight; negative log2FC means higher in ground control.",
        "",
        "The GLARE paper used DEG proportion, enrichment, heatmaps, and follow-up feature analyses to validate modules.",
        "This signed table is a stricter extension for the multi-study mouse analysis, because enrichment alone does not say whether a pathway is activated, suppressed, or mixed.",
        "",
        "## Direction Call Counts",
        "",
    ]
    counts = (
        meta.groupby(["module_class", "direction_call"], as_index=False)
        .size()
        .rename(columns={"size": "modules"})
        .sort_values(["module_class", "modules"], ascending=[True, False])
    )
    lines.extend(markdown_table(counts, ["module_class", "direction_call", "modules"], max_rows=50))
    lines.extend(["", "## Strict GLARE-Only Hidden Direction Candidates", ""])
    hidden = meta[
        meta["module_class"].eq("glare_only")
        & meta["direction_call"].astype(str).str.startswith("hidden module-score")
    ].sort_values(["combined_welch_fdr_bh", "median_empirical_abs_p"])
    lines.extend(
        markdown_table(
            hidden,
            [
                "tissue",
                "clean_term",
                "direction_call",
                "combined_welch_fdr_bh",
                "median_empirical_abs_p",
                "direction_consistency",
                "mean_flight_minus_ground",
                "significant_abs_lfc_genes_sum",
            ],
            max_rows=30,
        )
    )
    lines.extend(["", "## Mixed/Reorganized DEG Modules", ""])
    mixed = meta[meta["direction_call"].eq("mixed/reorganized DEG")].sort_values(
        ["module_class", "combined_welch_fdr_bh"], na_position="last"
    )
    lines.extend(
        markdown_table(
            mixed,
            [
                "tissue",
                "module_class",
                "clean_term",
                "flt_up_sig_genes_sum",
                "gc_up_sig_genes_sum",
                "combined_welch_fdr_bh",
                "mean_flight_minus_ground",
            ],
            max_rows=30,
        )
    )
    lines.extend(["", "## Output Files", ""])
    lines.extend(
        [
            "- `candidate_module_signed_dgea_by_study.tsv`: per-study signed counts for each module.",
            "- `candidate_module_signed_dgea_meta.tsv`: cross-study signed summary joined to module-score validation.",
        ]
    )
    (output_dir / "SIGNED_MODULE_DIRECTION_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> Path:
    root = Path(args.root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(args.candidate_modules, sep="\t")
    score_meta_path = Path(args.module_score_meta)
    module_score_meta = pd.read_csv(score_meta_path, sep="\t") if score_meta_path.exists() else pd.DataFrame()

    by_study = summarize_candidate_direction(
        root,
        candidates,
        args.alpha,
        args.lfc_cutoff,
        args.min_sig_genes,
        args.dominant_fraction,
    )
    meta = meta_summarize_direction(by_study, module_score_meta, args.min_sig_genes, args.dominant_fraction)

    by_study.to_csv(output_dir / "candidate_module_signed_dgea_by_study.tsv", sep="\t", index=False)
    meta.to_csv(output_dir / "candidate_module_signed_dgea_meta.tsv", sep="\t", index=False)
    write_report(meta, output_dir)
    return output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=DEFAULT_ROOT)
    parser.add_argument("--validation-dir", default=DEFAULT_VALIDATION_DIR)
    parser.add_argument("--candidate-modules", default="")
    parser.add_argument("--module-score-meta", default="")
    parser.add_argument("--output-dir", default="")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--lfc-cutoff", type=float, default=0.0)
    parser.add_argument("--min-sig-genes", type=int, default=3)
    parser.add_argument("--dominant-fraction", type=float, default=0.70)
    args = parser.parse_args()
    validation_dir = Path(args.validation_dir)
    if not args.candidate_modules:
        args.candidate_modules = str(validation_dir / "candidate_modules.tsv")
    if not args.module_score_meta:
        args.module_score_meta = str(validation_dir / "candidate_module_score_meta.tsv")
    if not args.output_dir:
        args.output_dir = str(validation_dir / "module_direction")
    return args


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
