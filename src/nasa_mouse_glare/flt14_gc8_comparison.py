"""Focused comparison of FLT14 and its primary GC8 counterpart."""

from __future__ import annotations

import argparse
import math
import time
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_RUN_DIR = "outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers"


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def to_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.fillna(False).map(
        lambda value: str(value).strip().lower() in {"true", "1", "yes"}
    )


def safe_mean(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return math.nan
    return float(values.mean())


def safe_percentile(series: pd.Series, q: float) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return math.nan
    return float(np.percentile(values, q))


def labels_for(group: pd.DataFrame, n: int = 20, *, significant: bool = False) -> str:
    subset = group.copy()
    if significant:
        subset = subset.loc[subset["significant_fdr05_abs_log2fc1"]].copy()
        if subset.empty:
            return ""
        subset["abs_lfc"] = subset["meta_log2_fold_change"].abs()
        subset = subset.sort_values("abs_lfc", ascending=False)
    else:
        subset = subset.sort_values("procrustes_latent_shift", ascending=False)
    return ",".join(subset["gene_label"].head(n).astype(str))


def load_table(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "post_analysis" / "paired_cluster_report" / "paired_gene_level_table.tsv"
    table = pd.read_csv(path, sep="\t")
    table["flt_cluster"] = table["flt_cluster"].astype(int)
    table["gc_cluster"] = table["gc_cluster"].astype(int)
    for column in ("eligible_meta", "significant_fdr05_abs_log2fc1"):
        table[column] = to_bool(table[column])
    if "gene_label" not in table:
        table["gene_label"] = table["gene_symbol"].where(
            table["gene_symbol"].notna(), table["gene_id"]
        )
    return table


def module_summary(
    name: str,
    definition: str,
    group: pd.DataFrame,
    flt_gene_ids: set[str],
    gc_gene_ids: set[str],
) -> dict[str, Any]:
    sig = group.loc[group["significant_fdr05_abs_log2fc1"]]
    eligible = group.loc[group["eligible_meta"]]
    gene_ids = set(group["gene_id"].astype(str))
    flt_overlap = gene_ids & flt_gene_ids
    gc_overlap = gene_ids & gc_gene_ids
    return {
        "module": name,
        "definition": definition,
        "gene_count": int(len(group)),
        "overlap_with_FLT14_genes": int(len(flt_overlap)),
        "fraction_of_FLT14": (
            float(len(flt_overlap) / len(flt_gene_ids)) if flt_gene_ids else math.nan
        ),
        "overlap_with_GC8_genes": int(len(gc_overlap)),
        "fraction_of_GC8": (
            float(len(gc_overlap) / len(gc_gene_ids)) if gc_gene_ids else math.nan
        ),
        "mean_latent_shift": safe_mean(group["procrustes_latent_shift"]),
        "median_latent_shift": safe_percentile(group["procrustes_latent_shift"], 50),
        "p90_latent_shift": safe_percentile(group["procrustes_latent_shift"], 90),
        "eligible_dgea_genes": int(len(eligible)),
        "significant_dgea_genes": int(len(sig)),
        "significant_up_genes": int((sig["meta_log2_fold_change"] > 0).sum()),
        "significant_down_genes": int((sig["meta_log2_fold_change"] < 0).sum()),
        "mean_meta_log2fc_eligible": safe_mean(eligible["meta_log2_fold_change"]),
        "mean_meta_log2fc_significant": safe_mean(sig["meta_log2_fold_change"]),
        "top_shift_genes": labels_for(group, n=20),
        "top_significant_genes": labels_for(group, n=20, significant=True),
    }


def summarize_partition(
    group: pd.DataFrame,
    by: str,
    denominator: int,
    label_column: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for value, subgroup in group.groupby(by, sort=True):
        sig = subgroup.loc[subgroup["significant_fdr05_abs_log2fc1"]]
        rows.append(
            {
                label_column: int(value),
                "gene_count": int(len(subgroup)),
                "fraction": float(len(subgroup) / denominator) if denominator else math.nan,
                "mean_latent_shift": safe_mean(subgroup["procrustes_latent_shift"]),
                "median_latent_shift": safe_percentile(
                    subgroup["procrustes_latent_shift"], 50
                ),
                "significant_dgea_genes": int(len(sig)),
                "top_shift_genes": labels_for(subgroup, n=12),
                "top_significant_genes": labels_for(
                    subgroup, n=12, significant=True
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("gene_count", ascending=False)


def overlap_rows(modules: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    sets = {
        name: set(frame["gene_id"].astype(str))
        for name, frame in modules.items()
    }
    for left, right in combinations(sets, 2):
        intersection = sets[left] & sets[right]
        union = sets[left] | sets[right]
        rows.append(
            {
                "left": left,
                "right": right,
                "left_genes": len(sets[left]),
                "right_genes": len(sets[right]),
                "intersection_genes": len(intersection),
                "jaccard": len(intersection) / len(union) if union else math.nan,
            }
        )
    return pd.DataFrame(rows).sort_values("intersection_genes", ascending=False)


def write_gene_lists(modules: dict[str, pd.DataFrame], output_dir: Path) -> None:
    gene_dir = output_dir / "gene_lists"
    gene_dir.mkdir(parents=True, exist_ok=True)
    lists: dict[str, pd.Series] = {}
    manifest = []
    for name, frame in modules.items():
        values = frame["gene_label"].dropna().astype(str).drop_duplicates()
        lists[name] = values.reset_index(drop=True)
        (gene_dir / f"{name}.txt").write_text(
            "\n".join(values.tolist()) + "\n", encoding="utf-8"
        )
        manifest.append({"list": name, "genes": int(len(values))})

    pd.DataFrame(lists).to_csv(gene_dir / "flt14_gc8_gene_lists.csv", index=False)
    pd.DataFrame(manifest).to_csv(
        gene_dir / "flt14_gc8_gene_list_manifest.tsv", sep="\t", index=False
    )


def markdown_table(df: pd.DataFrame, columns: list[str], n: int | None = None) -> list[str]:
    subset = df.loc[:, columns].copy()
    if n is not None:
        subset = subset.head(n)
    for column in subset.select_dtypes(include=[float]).columns:
        subset[column] = subset[column].map(lambda value: f"{value:.4g}")
    return subset.to_csv(sep="\t", index=False).strip().splitlines()


def write_report(
    output_dir: Path,
    module_df: pd.DataFrame,
    flt14_by_gc: pd.DataFrame,
    gc8_by_flt: pd.DataFrame,
    overlap: pd.DataFrame,
) -> None:
    flt14 = module_df.loc[module_df["module"] == "FLT14_all"].iloc[0]
    flt14_gc8 = module_df.loc[module_df["module"] == "FLT14_and_GC8"].iloc[0]
    flt14_not_gc8 = module_df.loc[module_df["module"] == "FLT14_not_GC8"].iloc[0]
    gc8 = module_df.loc[module_df["module"] == "GC8_all"].iloc[0]

    lines = [
        "# FLT14 vs GC8 Focused Comparison",
        "",
        "This compares `FLT14_all`, `FLT14_not_GC8`, and `GC8_all` from the",
        "12-filter MOBER + GLARE liver run. `FLT14_and_GC8` and",
        "`GC8_not_FLT14` are included as controls so the split is explicit.",
        "",
        "## Main Answer",
        "",
        f"- `GC8` is the primary GC match for `FLT14`, but only "
        f"{flt14_gc8['gene_count']:,}/{flt14['gene_count']:,} "
        f"FLT14 genes map there "
        f"({float(flt14_gc8['fraction_of_FLT14']):.1%}).",
        f"- `FLT14_not_GC8` contains {flt14_not_gc8['gene_count']:,} genes "
        f"({float(flt14_not_gc8['fraction_of_FLT14']):.1%} of FLT14).",
        f"- `GC8_all` contains {gc8['gene_count']:,} genes, but only "
        f"{float(flt14_gc8['fraction_of_GC8']):.1%} of GC8 comes from FLT14.",
        "- Therefore GC8 is the closest ground counterpart, not a clean",
        "  one-to-one liver equivalent.",
        "",
        "## Module Summary",
        "",
        "```tsv",
        *markdown_table(
            module_df,
            [
                "module",
                "gene_count",
                "overlap_with_FLT14_genes",
                "fraction_of_FLT14",
                "overlap_with_GC8_genes",
                "fraction_of_GC8",
                "mean_latent_shift",
                "significant_dgea_genes",
                "significant_up_genes",
                "significant_down_genes",
                "top_significant_genes",
            ],
        ),
        "```",
        "",
        "## Where FLT14 Genes Go In GC",
        "",
        "```tsv",
        *markdown_table(
            flt14_by_gc,
            [
                "gc_cluster",
                "gene_count",
                "fraction",
                "mean_latent_shift",
                "significant_dgea_genes",
                "top_shift_genes",
                "top_significant_genes",
            ],
            12,
        ),
        "```",
        "",
        "## What GC8 Is Made Of",
        "",
        "```tsv",
        *markdown_table(
            gc8_by_flt,
            [
                "flt_cluster",
                "gene_count",
                "fraction",
                "mean_latent_shift",
                "significant_dgea_genes",
                "top_shift_genes",
                "top_significant_genes",
            ],
            12,
        ),
        "```",
        "",
        "## Set Overlaps",
        "",
        "```tsv",
        *markdown_table(
            overlap,
            ["left", "right", "left_genes", "right_genes", "intersection_genes", "jaccard"],
            12,
        ),
        "```",
        "",
        "## Interpretation",
        "",
        "- `FLT14_and_GC8` is the matched portion of the module. It carries",
        "  `Mup15`, `Mup17`, and `Cdkn1a` among the significant DGEA-overlap",
        "  genes.",
        "- `FLT14_not_GC8` is the split/reorganized portion. Its only significant",
        "  DGEA-overlap gene in this analysis is `Apoa4`, but its high-shift gene",
        "  list includes liver metabolism genes mixed with contractile-associated",
        "  genes.",
        "- `GC8_all` is broad. Most GC8 genes come from FLT2 and FLT9 rather than",
        "  FLT14, so enrichment on all GC8 should be interpreted as the GC",
        "  neighborhood around FLT14, not as FLT14's direct equivalent.",
        "",
        "## Outputs",
        "",
        "- `flt14_gc8_module_summary.tsv`",
        "- `flt14_destination_gc_clusters.tsv`",
        "- `gc8_source_flt_clusters.tsv`",
        "- `flt14_gc8_set_overlaps.tsv`",
        "- `flt14_gc8_gene_membership.tsv`",
        "- `gene_lists/flt14_gc8_gene_lists.csv`",
    ]
    (output_dir / "FLT14_GC8_COMPARISON.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare FLT14, FLT14_not_GC8, and GC8 modules."
    )
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR)
    parser.add_argument("--flt-cluster", type=int, default=14)
    parser.add_argument("--gc-cluster", type=int, default=8)
    parser.add_argument(
        "--output-dir",
        help="Defaults to <run-dir>/post_analysis/flt14_gc8_comparison.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else run_dir / "post_analysis" / "flt14_gc8_comparison"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    log("Loading paired cluster gene table")
    table = load_table(run_dir)
    flt = int(args.flt_cluster)
    gc = int(args.gc_cluster)

    modules = {
        f"FLT{flt}_all": table.loc[table["flt_cluster"] == flt].copy(),
        f"FLT{flt}_and_GC{gc}": table.loc[
            (table["flt_cluster"] == flt) & (table["gc_cluster"] == gc)
        ].copy(),
        f"FLT{flt}_not_GC{gc}": table.loc[
            (table["flt_cluster"] == flt) & (table["gc_cluster"] != gc)
        ].copy(),
        f"GC{gc}_all": table.loc[table["gc_cluster"] == gc].copy(),
        f"GC{gc}_not_FLT{flt}": table.loc[
            (table["gc_cluster"] == gc) & (table["flt_cluster"] != flt)
        ].copy(),
    }

    flt_gene_ids = set(modules[f"FLT{flt}_all"]["gene_id"].astype(str))
    gc_gene_ids = set(modules[f"GC{gc}_all"]["gene_id"].astype(str))
    rows = [
        module_summary(name, definition, modules[name], flt_gene_ids, gc_gene_ids)
        for name, definition in [
            (f"FLT{flt}_all", f"All genes in FLT cluster {flt}"),
            (f"FLT{flt}_and_GC{gc}", f"FLT {flt} genes that also map to GC {gc}"),
            (f"FLT{flt}_not_GC{gc}", f"FLT {flt} genes outside GC {gc}"),
            (f"GC{gc}_all", f"All genes in GC cluster {gc}"),
            (f"GC{gc}_not_FLT{flt}", f"GC {gc} genes outside FLT {flt}"),
        ]
    ]
    module_df = pd.DataFrame(rows)

    flt14_by_gc = summarize_partition(
        modules[f"FLT{flt}_all"], "gc_cluster", len(flt_gene_ids), "gc_cluster"
    )
    gc8_by_flt = summarize_partition(
        modules[f"GC{gc}_all"], "flt_cluster", len(gc_gene_ids), "flt_cluster"
    )
    overlap = overlap_rows(modules)

    log("Writing focused comparison outputs")
    module_df.to_csv(output_dir / "flt14_gc8_module_summary.tsv", sep="\t", index=False, na_rep="NA")
    flt14_by_gc.to_csv(output_dir / "flt14_destination_gc_clusters.tsv", sep="\t", index=False, na_rep="NA")
    gc8_by_flt.to_csv(output_dir / "gc8_source_flt_clusters.tsv", sep="\t", index=False, na_rep="NA")
    overlap.to_csv(output_dir / "flt14_gc8_set_overlaps.tsv", sep="\t", index=False, na_rep="NA")

    membership = table.assign(
        in_FLT14_all=table["gene_id"].isin(modules[f"FLT{flt}_all"]["gene_id"]),
        in_FLT14_and_GC8=table["gene_id"].isin(modules[f"FLT{flt}_and_GC{gc}"]["gene_id"]),
        in_FLT14_not_GC8=table["gene_id"].isin(modules[f"FLT{flt}_not_GC{gc}"]["gene_id"]),
        in_GC8_all=table["gene_id"].isin(modules[f"GC{gc}_all"]["gene_id"]),
        in_GC8_not_FLT14=table["gene_id"].isin(modules[f"GC{gc}_not_FLT{flt}"]["gene_id"]),
    )
    membership.loc[
        membership[
            [
                "in_FLT14_all",
                "in_FLT14_and_GC8",
                "in_FLT14_not_GC8",
                "in_GC8_all",
                "in_GC8_not_FLT14",
            ]
        ].any(axis=1)
    ].to_csv(
        output_dir / "flt14_gc8_gene_membership.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    write_gene_lists(modules, output_dir)
    write_report(output_dir, module_df, flt14_by_gc, gc8_by_flt, overlap)
    log(f"Saved FLT{flt}/GC{gc} comparison to {output_dir}")


if __name__ == "__main__":
    main()
