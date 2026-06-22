"""Compare aggregate liver DESeq2 hits with FLT/GC GLARE clusters."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


DEFAULT_RUN_DIR = "outputs/glare_tms_liver_aggregated_osdr_flt_gc"


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return adjusted
    valid_p = p[valid]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    adjusted[valid] = restored
    return adjusted


def cluster_summary(merged: pd.DataFrame, cluster_column: str) -> pd.DataFrame:
    universe = merged.loc[merged["eligible_meta"].fillna(False)].copy()
    total_sig = int(universe["significant_fdr05_abs_log2fc1"].sum())
    total_not_sig = int(len(universe) - total_sig)
    rows = []

    for cluster_id, group in universe.groupby(cluster_column, sort=True):
        sig = group.loc[group["significant_fdr05_abs_log2fc1"]]
        n_sig = int(len(sig))
        n_total = int(len(group))
        n_not_sig = n_total - n_sig
        outside_sig = total_sig - n_sig
        outside_not_sig = total_not_sig - n_not_sig
        if n_total == 0 or total_sig == 0:
            odds_ratio = np.nan
            p_value = np.nan
        else:
            odds_ratio, p_value = fisher_exact(
                [[n_sig, n_not_sig], [outside_sig, outside_not_sig]],
                alternative="greater",
            )
        top_sig = (
            sig.assign(abs_lfc=sig["meta_log2_fold_change"].abs())
            .sort_values("abs_lfc", ascending=False)
            .head(10)
        )
        top_symbols = [
            symbol if isinstance(symbol, str) and symbol else gene
            for symbol, gene in zip(top_sig["gene_symbol"], top_sig["gene_id"])
        ]
        rows.append(
            {
                cluster_column: int(cluster_id),
                "eligible_genes": n_total,
                "significant_genes": n_sig,
                "significant_fraction": n_sig / n_total if n_total else np.nan,
                "fisher_odds_ratio": odds_ratio,
                "fisher_p_value": p_value,
                "mean_meta_log2fc_sig": (
                    float(sig["meta_log2_fold_change"].mean()) if n_sig else np.nan
                ),
                "mean_abs_meta_log2fc_sig": (
                    float(sig["meta_log2_fold_change"].abs().mean()) if n_sig else np.nan
                ),
                "mean_latent_shift_sig": (
                    float(sig["procrustes_latent_shift"].mean()) if n_sig else np.nan
                ),
                "top_significant_genes": ",".join(top_symbols),
            }
        )

    summary = pd.DataFrame(rows)
    summary["fisher_fdr_bh"] = bh_adjust(summary["fisher_p_value"].to_numpy())
    return summary.sort_values(
        ["significant_genes", "fisher_fdr_bh"], ascending=[False, True]
    )


def write_report(output_dir: Path, merged: pd.DataFrame, flt: pd.DataFrame, gc: pd.DataFrame) -> None:
    eligible = merged.loc[merged["eligible_meta"].fillna(False)]
    significant = eligible.loc[eligible["significant_fdr05_abs_log2fc1"]]
    flt_top = flt.head(8).to_csv(sep="\t", index=False).strip().splitlines()
    gc_top = gc.head(8).to_csv(sep="\t", index=False).strip().splitlines()
    text = [
        "# DESeq2 to GLARE Cluster Overlap",
        "",
        "DESeq2 aggregate meta-analysis hits were joined to the FLT and GC",
        "GLARE consensus clusters and Procrustes latent shifts.",
        "",
        f"- Eligible genes with GLARE clusters: {len(eligible):,}",
        f"- Significant DESeq2 meta genes: {len(significant):,}",
        "",
        "## FLT Clusters",
        "",
        "```tsv",
        *flt_top,
        "```",
        "",
        "## GC Clusters",
        "",
        "```tsv",
        *gc_top,
        "```",
    ]
    (output_dir / "DESEQ2_GLARE_OVERLAP_SUMMARY.md").write_text(
        "\n".join(text) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Join aggregate DESeq2 meta-analysis results to GLARE clusters."
    )
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR)
    parser.add_argument(
        "--deseq2-meta",
        help="Defaults to <run-dir>/post_analysis/deseq2_meta/deseq2_meta.tsv.",
    )
    parser.add_argument(
        "--output-dir",
        help="Defaults to <run-dir>/post_analysis/deseq2_glare_overlap.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    deseq2_meta = (
        Path(args.deseq2_meta)
        if args.deseq2_meta
        else run_dir / "post_analysis" / "deseq2_meta" / "deseq2_meta.tsv"
    )
    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else run_dir / "post_analysis" / "deseq2_glare_overlap"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    log("Loading DESeq2 and GLARE cluster comparison tables")
    meta = pd.read_csv(deseq2_meta, sep="\t")
    clusters = pd.read_csv(
        run_dir / "post_analysis" / "gene_cluster_comparison.tsv", sep="\t"
    )
    merged = clusters.merge(meta, on="gene_id", how="left")
    merged["eligible_meta"] = merged["eligible_meta"].fillna(False).astype(bool)
    merged["significant_fdr05_abs_log2fc1"] = (
        merged["significant_fdr05_abs_log2fc1"].fillna(False).astype(bool)
    )

    log("Summarizing DESeq2 hits by FLT and GC clusters")
    flt = cluster_summary(merged, "flt_cluster")
    gc = cluster_summary(merged, "gc_cluster")

    merged.to_csv(
        output_dir / "deseq2_genes_with_glare_clusters.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    flt.to_csv(
        output_dir / "deseq2_by_flt_cluster.tsv", sep="\t", index=False, na_rep="NA"
    )
    gc.to_csv(
        output_dir / "deseq2_by_gc_cluster.tsv", sep="\t", index=False, na_rep="NA"
    )
    write_report(output_dir, merged, flt, gc)
    log(f"Saved DESeq2/GLARE overlap outputs to {output_dir}")


if __name__ == "__main__":
    main()
