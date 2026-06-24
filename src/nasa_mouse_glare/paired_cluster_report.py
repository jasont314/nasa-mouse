"""Build paired FLT-vs-GC GLARE cluster reports.

The aggregate liver GLARE runs fine-tune separate SAE models for flight and
ground-control samples. This report compares those two gene representations by
using the existing Procrustes-aligned latent shift and cluster-overlap tables.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_RUN_DIR = "outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers"


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def clean_label(row: pd.Series) -> str:
    symbol = row.get("gene_symbol")
    if isinstance(symbol, str) and symbol and symbol != "NA":
        return symbol
    return str(row["gene_id"])


def to_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.fillna(False).map(
        lambda value: str(value).strip().lower() in {"true", "1", "yes"}
    )


def read_optional_tsv(path: Path) -> pd.DataFrame | None:
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


def load_gene_table(run_dir: Path) -> pd.DataFrame:
    raw_overlap = (
        run_dir
        / "post_analysis"
        / "raw_count_deseq2_glare_overlap"
        / "deseq2_genes_with_glare_clusters.tsv"
    )
    if raw_overlap.exists():
        table = pd.read_csv(raw_overlap, sep="\t")
    else:
        table = pd.read_csv(
            run_dir / "post_analysis" / "gene_cluster_comparison.tsv", sep="\t"
        )
        meta = read_optional_tsv(run_dir / "post_analysis" / "meta_dgea.tsv")
        if meta is not None:
            table = table.merge(meta, on="gene_id", how="left")

    for column in ("flt_cluster", "gc_cluster"):
        table[column] = table[column].astype(int)
    for column in ("eligible_meta", "significant_fdr05_abs_log2fc1"):
        if column in table:
            table[column] = to_bool(table[column])
        else:
            table[column] = False
    if "gene_symbol" not in table:
        table["gene_symbol"] = ""
    table["gene_label"] = table.apply(clean_label, axis=1)
    return table


def top_labels(
    group: pd.DataFrame,
    *,
    n: int = 12,
    sort_column: str = "procrustes_latent_shift",
    descending: bool = True,
    only_significant: bool = False,
) -> str:
    subset = group
    if only_significant:
        subset = subset.loc[subset["significant_fdr05_abs_log2fc1"]]
    if subset.empty:
        return ""
    if sort_column == "abs_meta_log2_fold_change":
        subset = subset.assign(
            abs_meta_log2_fold_change=subset["meta_log2_fold_change"].abs()
        )
    if sort_column in subset:
        subset = subset.sort_values(sort_column, ascending=not descending)
    return ",".join(subset["gene_label"].head(n).astype(str))


def safe_mean(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce")
    if values.notna().any():
        return float(values.mean())
    return math.nan


def safe_percentile(series: pd.Series, q: float) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return math.nan
    return float(np.percentile(values, q))


def summarize_edges(table: pd.DataFrame) -> pd.DataFrame:
    flt_totals = table.groupby("flt_cluster")["gene_id"].size()
    gc_totals = table.groupby("gc_cluster")["gene_id"].size()
    rows: list[dict[str, Any]] = []

    for (flt_cluster, gc_cluster), group in table.groupby(
        ["flt_cluster", "gc_cluster"], sort=True
    ):
        eligible = group.loc[group["eligible_meta"]]
        sig = group.loc[group["significant_fdr05_abs_log2fc1"]]
        rows.append(
            {
                "flt_cluster": int(flt_cluster),
                "gc_cluster": int(gc_cluster),
                "gene_count": int(len(group)),
                "fraction_of_flt_cluster": float(
                    len(group) / flt_totals.loc[flt_cluster]
                ),
                "fraction_of_gc_cluster": float(len(group) / gc_totals.loc[gc_cluster]),
                "mean_latent_shift": safe_mean(group["procrustes_latent_shift"]),
                "median_latent_shift": safe_percentile(
                    group["procrustes_latent_shift"], 50
                ),
                "p90_latent_shift": safe_percentile(
                    group["procrustes_latent_shift"], 90
                ),
                "eligible_dgea_genes": int(len(eligible)),
                "significant_dgea_genes": int(len(sig)),
                "mean_meta_log2fc_eligible": safe_mean(
                    eligible.get("meta_log2_fold_change", pd.Series(dtype=float))
                ),
                "mean_abs_meta_log2fc_sig": safe_mean(
                    sig.get("meta_log2_fold_change", pd.Series(dtype=float)).abs()
                ),
                "top_shift_genes": top_labels(group, n=12),
                "top_significant_genes": top_labels(
                    group,
                    n=12,
                    sort_column="abs_meta_log2_fold_change",
                    only_significant=True,
                ),
            }
        )

    edges = pd.DataFrame(rows)
    flt_top = edges.sort_values(
        ["flt_cluster", "gene_count"], ascending=[True, False]
    ).drop_duplicates("flt_cluster")
    gc_top = edges.sort_values(
        ["gc_cluster", "gene_count"], ascending=[True, False]
    ).drop_duplicates("gc_cluster")
    edges = edges.merge(
        flt_top[["flt_cluster", "gc_cluster"]].rename(
            columns={"gc_cluster": "primary_gc_for_flt"}
        ),
        on="flt_cluster",
        how="left",
    )
    edges = edges.merge(
        gc_top[["gc_cluster", "flt_cluster"]].rename(
            columns={"flt_cluster": "primary_flt_for_gc"}
        ),
        on="gc_cluster",
        how="left",
    )
    edges["is_primary_flt_to_gc"] = edges["gc_cluster"].eq(
        edges["primary_gc_for_flt"]
    )
    edges["is_primary_gc_to_flt"] = edges["flt_cluster"].eq(
        edges["primary_flt_for_gc"]
    )
    edges["edge_role"] = np.select(
        [
            edges["is_primary_flt_to_gc"] & edges["is_primary_gc_to_flt"],
            edges["is_primary_flt_to_gc"],
            edges["is_primary_gc_to_flt"],
        ],
        ["reciprocal_primary", "flt_primary_only", "gc_primary_only"],
        default="secondary",
    )
    return edges.sort_values(["flt_cluster", "gene_count"], ascending=[True, False])


def compact_destinations(edges: pd.DataFrame, flt_cluster: int, skip_gc: int) -> str:
    subset = edges.loc[
        (edges["flt_cluster"] == flt_cluster) & (edges["gc_cluster"] != skip_gc)
    ].sort_values("gene_count", ascending=False)
    parts = []
    for row in subset.head(6).itertuples(index=False):
        parts.append(
            f"GC{int(row.gc_cluster)}:{int(row.gene_count)}"
            f"({float(row.fraction_of_flt_cluster):.2f})"
        )
    return ",".join(parts)


def classify_cluster(match_fraction: float, mean_shift: float, shift_p75: float) -> str:
    if match_fraction >= 0.85 and mean_shift < shift_p75:
        return "conserved"
    if match_fraction >= 0.85:
        return "matched_but_shifted"
    if match_fraction >= 0.65 and mean_shift >= shift_p75:
        return "partially_split_shifted"
    if match_fraction < 0.65 and mean_shift >= shift_p75:
        return "split_reorganized"
    if match_fraction < 0.65:
        return "split_low_shift"
    return "partially_split"


def summarize_flt_clusters(table: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    shift_p75 = float(table["procrustes_latent_shift"].quantile(0.75))
    rows: list[dict[str, Any]] = []
    top_edges = edges.loc[edges["is_primary_flt_to_gc"]].copy()

    for edge in top_edges.sort_values("flt_cluster").itertuples(index=False):
        flt_cluster = int(edge.flt_cluster)
        group = table.loc[table["flt_cluster"] == flt_cluster]
        eligible = group.loc[group["eligible_meta"]]
        sig = group.loc[group["significant_fdr05_abs_log2fc1"]]
        gc_clusters_spanned = int(group["gc_cluster"].nunique())
        match_fraction = float(edge.fraction_of_flt_cluster)
        mean_shift = safe_mean(group["procrustes_latent_shift"])
        priority_score = (
            mean_shift
            * (1.0 - match_fraction)
            * math.log2(gc_clusters_spanned + 1)
        )
        rows.append(
            {
                "flt_cluster": flt_cluster,
                "matched_gc_cluster": int(edge.gc_cluster),
                "gene_count": int(len(group)),
                "matched_gene_count": int(edge.gene_count),
                "matched_fraction": match_fraction,
                "gc_clusters_spanned": gc_clusters_spanned,
                "secondary_gc_clusters": compact_destinations(
                    edges, flt_cluster, int(edge.gc_cluster)
                ),
                "mean_latent_shift": mean_shift,
                "median_latent_shift": safe_percentile(
                    group["procrustes_latent_shift"], 50
                ),
                "p90_latent_shift": safe_percentile(
                    group["procrustes_latent_shift"], 90
                ),
                "eligible_dgea_genes": int(len(eligible)),
                "significant_dgea_genes": int(len(sig)),
                "mean_meta_log2fc_sig": safe_mean(
                    sig.get("meta_log2_fold_change", pd.Series(dtype=float))
                ),
                "top_shift_genes": top_labels(group, n=12),
                "top_significant_genes": top_labels(
                    group,
                    n=12,
                    sort_column="abs_meta_log2_fold_change",
                    only_significant=True,
                ),
                "cluster_status": classify_cluster(
                    match_fraction, mean_shift, shift_p75
                ),
                "reorganization_priority_score": float(priority_score),
            }
        )

    summary = pd.DataFrame(rows)
    return summary.sort_values(
        ["reorganization_priority_score", "mean_latent_shift"],
        ascending=[False, False],
    )


def summarize_gc_clusters(table: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    top_edges = edges.loc[edges["is_primary_gc_to_flt"]].copy()
    for edge in top_edges.sort_values("gc_cluster").itertuples(index=False):
        gc_cluster = int(edge.gc_cluster)
        group = table.loc[table["gc_cluster"] == gc_cluster]
        rows.append(
            {
                "gc_cluster": gc_cluster,
                "matched_flt_cluster": int(edge.flt_cluster),
                "gene_count": int(len(group)),
                "matched_gene_count": int(edge.gene_count),
                "matched_fraction": float(edge.fraction_of_gc_cluster),
                "flt_clusters_spanned": int(group["flt_cluster"].nunique()),
                "mean_latent_shift": safe_mean(group["procrustes_latent_shift"]),
                "median_latent_shift": safe_percentile(
                    group["procrustes_latent_shift"], 50
                ),
                "significant_dgea_genes": int(
                    group["significant_fdr05_abs_log2fc1"].sum()
                ),
                "top_shift_genes": top_labels(group, n=12),
                "top_significant_genes": top_labels(
                    group,
                    n=12,
                    sort_column="abs_meta_log2_fold_change",
                    only_significant=True,
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["mean_latent_shift", "matched_fraction"], ascending=[False, True]
    )


def summarize_displaced_modules(table: pd.DataFrame, flt_summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in flt_summary.itertuples(index=False):
        flt_cluster = int(row.flt_cluster)
        matched_gc = int(row.matched_gc_cluster)
        group = table.loc[table["flt_cluster"] == flt_cluster]
        displaced = group.loc[group["gc_cluster"] != matched_gc]
        if displaced.empty:
            continue
        sig = displaced.loc[displaced["significant_fdr05_abs_log2fc1"]]
        rows.append(
            {
                "module": f"FLT{flt_cluster}_not_GC{matched_gc}",
                "flt_cluster": flt_cluster,
                "excluded_matched_gc_cluster": matched_gc,
                "gene_count": int(len(displaced)),
                "fraction_of_flt_cluster": float(len(displaced) / len(group)),
                "gc_clusters_spanned": int(displaced["gc_cluster"].nunique()),
                "mean_latent_shift": safe_mean(displaced["procrustes_latent_shift"]),
                "median_latent_shift": safe_percentile(
                    displaced["procrustes_latent_shift"], 50
                ),
                "significant_dgea_genes": int(len(sig)),
                "mean_meta_log2fc_sig": safe_mean(
                    sig.get("meta_log2_fold_change", pd.Series(dtype=float))
                ),
                "top_shift_genes": top_labels(displaced, n=20),
                "top_significant_genes": top_labels(
                    displaced,
                    n=20,
                    sort_column="abs_meta_log2_fold_change",
                    only_significant=True,
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["mean_latent_shift", "fraction_of_flt_cluster"], ascending=[False, False]
    )


def write_gene_lists(
    table: pd.DataFrame,
    flt_summary: pd.DataFrame,
    modules: pd.DataFrame,
    output_dir: Path,
    top_n: int,
) -> None:
    gene_list_dir = output_dir / "gene_lists"
    gene_list_dir.mkdir(parents=True, exist_ok=True)
    lists: dict[str, pd.Series] = {}

    for row in flt_summary.head(top_n).itertuples(index=False):
        cluster = int(row.flt_cluster)
        matched_gc = int(row.matched_gc_cluster)
        group = table.loc[table["flt_cluster"] == cluster]
        displaced = group.loc[group["gc_cluster"] != matched_gc]
        lists[f"FLT{cluster}_all"] = group["gene_label"].dropna().astype(str)
        if len(displaced) >= 20:
            lists[f"FLT{cluster}_not_GC{matched_gc}"] = (
                displaced["gene_label"].dropna().astype(str)
            )

    for row in modules.head(top_n).itertuples(index=False):
        module = str(row.module)
        if module in lists:
            continue
        group = table.loc[
            (table["flt_cluster"] == int(row.flt_cluster))
            & (table["gc_cluster"] != int(row.excluded_matched_gc_cluster))
        ]
        lists[module] = group["gene_label"].dropna().astype(str)

    if not lists:
        return
    pd.DataFrame({name: values.reset_index(drop=True) for name, values in lists.items()}).to_csv(
        gene_list_dir / "paired_reorganized_gene_lists.csv", index=False
    )
    manifest = pd.DataFrame(
        [{"list": name, "genes": int(values.nunique())} for name, values in lists.items()]
    )
    manifest.to_csv(gene_list_dir / "paired_reorganized_gene_list_manifest.tsv", sep="\t", index=False)


def markdown_table(df: pd.DataFrame, columns: list[str], n: int) -> list[str]:
    subset = df.loc[:, columns].head(n).copy()
    for column in subset.select_dtypes(include=[float]).columns:
        subset[column] = subset[column].map(lambda value: f"{value:.3g}")
    return subset.to_csv(sep="\t", index=False).strip().splitlines()


def write_report(
    run_dir: Path,
    output_dir: Path,
    table: pd.DataFrame,
    edges: pd.DataFrame,
    flt_summary: pd.DataFrame,
    gc_summary: pd.DataFrame,
    modules: pd.DataFrame,
) -> None:
    summary_json = read_optional_json(
        run_dir / "post_analysis" / "post_analysis_summary.json"
    )
    cluster_summary = (summary_json or {}).get("cluster_comparison", {})
    ari = cluster_summary.get("adjusted_rand_index")
    nmi = cluster_summary.get("normalized_mutual_information")

    top_flt = markdown_table(
        flt_summary,
        [
            "flt_cluster",
            "matched_gc_cluster",
            "gene_count",
            "matched_fraction",
            "gc_clusters_spanned",
            "mean_latent_shift",
            "significant_dgea_genes",
            "cluster_status",
            "top_shift_genes",
        ],
        10,
    )
    top_modules = markdown_table(
        modules,
        [
            "module",
            "gene_count",
            "fraction_of_flt_cluster",
            "gc_clusters_spanned",
            "mean_latent_shift",
            "significant_dgea_genes",
            "top_shift_genes",
        ],
        10,
    )
    top_edges = markdown_table(
        edges.sort_values(["mean_latent_shift", "gene_count"], ascending=[False, False]),
        [
            "flt_cluster",
            "gc_cluster",
            "edge_role",
            "gene_count",
            "fraction_of_flt_cluster",
            "fraction_of_gc_cluster",
            "mean_latent_shift",
            "significant_dgea_genes",
            "top_shift_genes",
        ],
        12,
    )

    lines = [
        "# Paired FLT/GC GLARE Cluster Report",
        "",
        "This report compares flight and ground-control liver GLARE gene",
        "representations after MOBER correction and removal of the 12 candidate",
        "muscle-outlier profiles. FLT and GC were fine-tuned as separate SAE",
        "models, so latent-space distances use the existing Procrustes alignment",
        "rather than direct unaligned coordinate subtraction.",
        "",
        "## Overall Structure",
        "",
        f"- Genes compared: {len(table):,}",
        f"- FLT clusters: {table['flt_cluster'].nunique()}",
        f"- GC clusters: {table['gc_cluster'].nunique()}",
        f"- Adjusted Rand index: {ari:.4f}" if ari is not None else "- Adjusted Rand index: NA",
        f"- Normalized mutual information: {nmi:.4f}" if nmi is not None else "- Normalized mutual information: NA",
        f"- Median Procrustes latent shift: {table['procrustes_latent_shift'].median():.4f}",
        "",
        "A high match fraction means most genes in a FLT cluster remain together",
        "in one GC cluster. A high latent shift means those genes move farther",
        "after aligning the two SAE spaces. The strongest candidates are therefore",
        "clusters with both split structure and high latent shift.",
        "",
        "## Top FLT Reorganization Candidates",
        "",
        "```tsv",
        *top_flt,
        "```",
        "",
        "## Displaced FLT Modules",
        "",
        "`FLT##_not_GC##` means genes from a FLT cluster that do not map to that",
        "cluster's primary GC counterpart.",
        "",
        "```tsv",
        *top_modules,
        "```",
        "",
        "## Highest-Shift FLT-to-GC Edges",
        "",
        "```tsv",
        *top_edges,
        "```",
        "",
        "## Interpretation",
        "",
        "- FLT cluster 14 is the strongest non-obvious reorganization candidate:",
        "  it has high latent shift and only about half of its genes map to its",
        "  top GC counterpart. It still needs gene-level review because some",
        "  high-shift genes are contractile or smooth-muscle associated.",
        "- FLT clusters 13 and 15 are high-shift and DGEA-enriched, but their",
        "  top genes are dominated by contractile/muscle markers, so they should",
        "  be treated as residual composition or contamination-sensitive modules.",
        "- FLT cluster 3 is mostly conserved with GC cluster 1, despite having many",
        "  DGEA-overlapping genes. That pattern is more consistent with a stable",
        "  module whose expression changes than with a strongly reorganized module.",
        "- Cluster-level reorganization is not the same thing as flight-up",
        "  expression. Direction still needs to be read from the DGEA columns and",
        "  per-study consistency.",
        "",
        "## Outputs",
        "",
        "- `flt_to_gc_paired_cluster_summary.tsv`",
        "- `gc_to_flt_paired_cluster_summary.tsv`",
        "- `flt_gc_cluster_edges.tsv`",
        "- `flt_displaced_modules.tsv`",
        "- `paired_gene_level_table.tsv`",
        "- `gene_lists/paired_reorganized_gene_lists.csv`",
    ]
    (output_dir / "PAIRED_CLUSTER_REPORT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def read_optional_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create paired FLT-vs-GC GLARE cluster reorganization reports."
    )
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR)
    parser.add_argument(
        "--output-dir",
        help="Defaults to <run-dir>/post_analysis/paired_cluster_report.",
    )
    parser.add_argument("--top-gene-list-clusters", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else run_dir / "post_analysis" / "paired_cluster_report"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    log("Loading gene-level FLT/GC cluster comparison")
    table = load_gene_table(run_dir)
    edges = summarize_edges(table)
    flt_summary = summarize_flt_clusters(table, edges)
    gc_summary = summarize_gc_clusters(table, edges)
    modules = summarize_displaced_modules(table, flt_summary)

    log("Writing paired cluster outputs")
    table.to_csv(
        output_dir / "paired_gene_level_table.tsv", sep="\t", index=False, na_rep="NA"
    )
    edges.to_csv(output_dir / "flt_gc_cluster_edges.tsv", sep="\t", index=False, na_rep="NA")
    flt_summary.to_csv(
        output_dir / "flt_to_gc_paired_cluster_summary.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    gc_summary.to_csv(
        output_dir / "gc_to_flt_paired_cluster_summary.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    modules.to_csv(
        output_dir / "flt_displaced_modules.tsv", sep="\t", index=False, na_rep="NA"
    )
    write_gene_lists(
        table, flt_summary, modules, output_dir, args.top_gene_list_clusters
    )
    write_report(run_dir, output_dir, table, edges, flt_summary, gc_summary, modules)
    log(f"Saved paired cluster report to {output_dir}")


if __name__ == "__main__":
    main()
