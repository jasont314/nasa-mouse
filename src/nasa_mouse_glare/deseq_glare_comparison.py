"""Compare filtered OSD-379 DESeq2 results with GLARE representations."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.linalg import orthogonal_procrustes
from scipy.stats import fisher_exact, hypergeom, mannwhitneyu, spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

from .cluster_enrichment import read_gmt

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def cluster_enrichment(
    gene_table: pd.DataFrame, clusters: pd.DataFrame, location: str
) -> pd.DataFrame:
    merged = clusters[["gene_id", "consensus"]].merge(
        gene_table[["gene_id", "is_deg", "direction", "log2FoldChange"]],
        on="gene_id",
        validate="one_to_one",
    )
    total_deg = int(merged["is_deg"].sum())
    total_non_deg = int((~merged["is_deg"]).sum())
    records = []
    for cluster, group in merged.groupby("consensus"):
        cluster_deg = int(group["is_deg"].sum())
        cluster_non_deg = int((~group["is_deg"]).sum())
        odds_ratio, p_value = fisher_exact(
            [
                [cluster_deg, cluster_non_deg],
                [total_deg - cluster_deg, total_non_deg - cluster_non_deg],
            ],
            alternative="greater",
        )
        records.append(
            {
                "location": location,
                "cluster": int(cluster),
                "gene_count": int(len(group)),
                "deg_count": cluster_deg,
                "deg_proportion": float(group["is_deg"].mean()),
                "up_count": int((group["direction"] == "up").sum()),
                "down_count": int((group["direction"] == "down").sum()),
                "median_log2_fold_change": float(group["log2FoldChange"].median()),
                "fisher_odds_ratio": float(odds_ratio),
                "fisher_p_value": float(p_value),
            }
        )
    result = pd.DataFrame(records)
    result["fisher_fdr_bh"] = multipletests(
        result["fisher_p_value"], method="fdr_bh"
    )[1]
    return result.sort_values(
        ["fisher_fdr_bh", "deg_proportion"], ascending=[True, False]
    )


def procrustes_shift(flt: np.ndarray, gc: np.ndarray) -> np.ndarray:
    flt_scaled = StandardScaler().fit_transform(flt)
    gc_scaled = StandardScaler().fit_transform(gc)
    rotation, scale = orthogonal_procrustes(gc_scaled, flt_scaled)
    aligned_gc = gc_scaled @ rotation * scale
    aligned_gc /= max(scale, np.finfo(float).eps)
    return np.linalg.norm(flt_scaled - aligned_gc, axis=1)


def normalize_entrez(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text.split(";")[0].split(",")[0].strip()


def load_annotation(path: str | Path) -> pd.DataFrame:
    rows = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        indices = {name: header.index(name) for name in ("ENSEMBL", "SYMBOL", "ENTREZID")}
        for row in reader:
            rows.append(
                {
                    name: row[index] if index < len(row) else ""
                    for name, index in indices.items()
                }
            )
    return pd.DataFrame(rows).drop_duplicates("ENSEMBL")


def pathway_enrichment(
    query: set[str],
    universe: set[str],
    gmt_path: Path,
    method: str,
    min_overlap: int = 2,
) -> pd.DataFrame:
    records = []
    query = query & universe
    for gene_set in read_gmt(gmt_path):
        term_genes = gene_set["genes"] & universe
        overlap = query & term_genes
        if len(overlap) < min_overlap:
            continue
        p_value = hypergeom.sf(
            len(overlap) - 1,
            len(universe),
            len(term_genes),
            len(query),
        )
        records.append(
            {
                "method": method,
                "term": gene_set["term"],
                "description": gene_set["description"],
                "query_size": len(query),
                "universe_size": len(universe),
                "term_size": len(term_genes),
                "overlap": len(overlap),
                "enrichment_ratio": (
                    (len(overlap) / len(query))
                    / (len(term_genes) / len(universe))
                ),
                "p_value": float(p_value),
                "overlap_gene_ids": ",".join(sorted(overlap)),
            }
        )
    result = pd.DataFrame(records)
    if result.empty:
        return result
    result["fdr_bh"] = multipletests(result["p_value"], method="fdr_bh")[1]
    return result.sort_values(["fdr_bh", "p_value", "overlap"])


def plot_comparison(gene_table: pd.DataFrame, output_dir: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    colors = np.where(gene_table["is_deg"], "#c43c39", "#8c8c8c")
    axes[0].scatter(
        gene_table["abs_log2_fold_change"],
        gene_table["procrustes_latent_shift"],
        c=colors,
        s=7,
        alpha=0.45,
        linewidths=0,
    )
    axes[0].set_xlabel("Absolute DESeq2 log2 fold change")
    axes[0].set_ylabel("Procrustes-aligned GLARE latent shift")
    axes[0].spines[["top", "right"]].set_visible(False)

    non_deg = gene_table.loc[~gene_table["is_deg"], "procrustes_latent_shift"]
    deg = gene_table.loc[gene_table["is_deg"], "procrustes_latent_shift"]
    axes[1].boxplot(
        [non_deg, deg],
        tick_labels=["Not DEG", "DEG"],
        showfliers=False,
        patch_artist=True,
        boxprops={"facecolor": "#d0d0d0"},
        medianprops={"color": "#202020"},
    )
    axes[1].set_ylabel("Procrustes-aligned GLARE latent shift")
    axes[1].spines[["top", "right"]].set_visible(False)
    figure.savefig(output_dir / "deseq2_vs_glare.png", dpi=180)
    plt.close(figure)


def run(args: argparse.Namespace) -> dict:
    run_dir = Path(args.run_dir)
    deseq_dir = Path(args.deseq_dir)
    output_dir = run_dir / "deseq2_glare_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)

    deseq = pd.read_csv(deseq_dir / "global_flight_vs_ground.tsv", sep="\t")
    target = np.load(run_dir / "controlled_target.npz")
    genes = target["genes"].astype(str)
    flt_latent = np.load(run_dir / "FLT_FTSAE_representation.npy")
    gc_latent = np.load(run_dir / "GC_FTSAE_representation.npy")
    shifts = procrustes_shift(flt_latent, gc_latent)

    gene_table = pd.DataFrame(
        {"gene_id": genes, "procrustes_latent_shift": shifts}
    ).merge(deseq, on="gene_id", how="left", validate="one_to_one")
    annotation = load_annotation(args.annotation)
    annotation["ENSEMBL"] = annotation["ENSEMBL"].astype(str)
    annotation["entrez_id"] = annotation["ENTREZID"].map(normalize_entrez)
    gene_table = gene_table.merge(
        annotation[["ENSEMBL", "SYMBOL", "entrez_id"]],
        left_on="gene_id",
        right_on="ENSEMBL",
        how="left",
        validate="one_to_one",
    ).drop(columns="ENSEMBL")
    gene_table["is_deg"] = gene_table[
        "significant_padj05_abs_lfc1"
    ].fillna(False).astype(bool)
    gene_table["direction"] = gene_table["direction"].fillna("not_tested")
    gene_table["abs_log2_fold_change"] = gene_table["log2FoldChange"].abs()
    gene_table["latent_shift_rank"] = gene_table[
        "procrustes_latent_shift"
    ].rank(method="min", ascending=False)
    gene_table.to_csv(output_dir / "gene_level_comparison.tsv", sep="\t", index=False)

    enrichments = []
    cluster_tables = {}
    for location in ("FLT", "GC"):
        clusters = pd.read_csv(
            run_dir / "clustering" / f"{location}_gene_clusters.tsv", sep="\t"
        )
        cluster_tables[location] = clusters
        enrichments.append(cluster_enrichment(gene_table, clusters, location))
    enrichment = pd.concat(enrichments, ignore_index=True)
    enrichment.to_csv(
        output_dir / "cluster_deg_enrichment.tsv", sep="\t", index=False
    )

    cluster_comparison = cluster_tables["FLT"][["gene_id", "consensus"]].merge(
        cluster_tables["GC"][["gene_id", "consensus"]],
        on="gene_id",
        suffixes=("_flt", "_gc"),
        validate="one_to_one",
    )
    cluster_ari = adjusted_rand_score(
        cluster_comparison["consensus_flt"], cluster_comparison["consensus_gc"]
    )
    cluster_nmi = normalized_mutual_info_score(
        cluster_comparison["consensus_flt"], cluster_comparison["consensus_gc"]
    )

    tested = gene_table.dropna(
        subset=["log2FoldChange", "procrustes_latent_shift"]
    ).copy()
    correlation = spearmanr(
        tested["abs_log2_fold_change"],
        tested["procrustes_latent_shift"],
        nan_policy="omit",
    )
    deg_shift = tested.loc[tested["is_deg"], "procrustes_latent_shift"]
    non_deg_shift = tested.loc[~tested["is_deg"], "procrustes_latent_shift"]
    if len(deg_shift) and len(non_deg_shift):
        shift_test = mannwhitneyu(deg_shift, non_deg_shift, alternative="greater")
        latent_auc = roc_auc_score(
            tested["is_deg"], tested["procrustes_latent_shift"]
        )
    else:
        shift_test = None
        latent_auc = float("nan")

    significant_clusters = enrichment.loc[enrichment["fisher_fdr_bh"] < 0.05]
    plot_comparison(gene_table, output_dir)
    top_degs = (
        tested.loc[tested["is_deg"]]
        .sort_values(["padj", "abs_log2_fold_change"], ascending=[True, False])
        .head(25)
    )
    top_degs.to_csv(output_dir / "top_deseq2_genes.tsv", sep="\t", index=False)

    pathway_tested = tested.copy()
    universe_ids = set(pathway_tested["gene_id"])
    deg_ids = set(pathway_tested.loc[pathway_tested["is_deg"], "gene_id"])
    top_shift = (
        pathway_tested.sort_values("procrustes_latent_shift", ascending=False)
        .drop_duplicates("gene_id")
        .head(len(deg_ids))
    )
    shift_ids = set(top_shift["gene_id"])
    top_shift.to_csv(
        output_dir / "top_glare_latent_shift_genes.tsv", sep="\t", index=False
    )
    deseq_pathways = pathway_enrichment(
        deg_ids,
        universe_ids,
        Path(args.reactome_gmt),
        "DESeq2_DEG",
    )
    glare_pathways = pathway_enrichment(
        shift_ids,
        universe_ids,
        Path(args.reactome_gmt),
        "GLARE_top_latent_shift",
    )
    deseq_pathways.to_csv(
        output_dir / "reactome_deseq2_deg.tsv", sep="\t", index=False
    )
    glare_pathways.to_csv(
        output_dir / "reactome_glare_top_shift.tsv", sep="\t", index=False
    )
    significant_deseq_terms = (
        set(deseq_pathways.loc[deseq_pathways["fdr_bh"] < 0.05, "term"])
        if not deseq_pathways.empty
        else set()
    )
    significant_glare_terms = (
        set(glare_pathways.loc[glare_pathways["fdr_bh"] < 0.05, "term"])
        if not glare_pathways.empty
        else set()
    )
    pathway_comparison = pd.DataFrame(
        {
            "term": sorted(significant_deseq_terms | significant_glare_terms),
        }
    )
    if not pathway_comparison.empty:
        pathway_comparison["significant_deseq2"] = pathway_comparison["term"].isin(
            significant_deseq_terms
        )
        pathway_comparison["significant_glare_shift"] = pathway_comparison[
            "term"
        ].isin(significant_glare_terms)
    pathway_comparison.to_csv(
        output_dir / "significant_reactome_pathway_comparison.tsv",
        sep="\t",
        index=False,
    )

    reference_comparison = {}
    original_cluster12 = pd.DataFrame()
    reference_run = Path(args.reference_run) if args.reference_run else None
    if reference_run and reference_run.exists():
        for location in ("FLT", "GC"):
            original = pd.read_csv(
                reference_run / "clustering" / f"{location}_gene_clusters.tsv",
                sep="\t",
            )
            filtered = cluster_tables[location]
            aligned = original[["gene_id", "consensus"]].merge(
                filtered[["gene_id", "consensus"]],
                on="gene_id",
                suffixes=("_original", "_filtered"),
                validate="one_to_one",
            )
            reference_comparison[location] = {
                "ari": float(
                    adjusted_rand_score(
                        aligned["consensus_original"],
                        aligned["consensus_filtered"],
                    )
                ),
                "nmi": float(
                    normalized_mutual_info_score(
                        aligned["consensus_original"],
                        aligned["consensus_filtered"],
                    )
                ),
            }
            if location == "FLT":
                original_cluster12 = (
                    original.loc[original["consensus"].eq(12), ["gene_id"]]
                    .merge(
                        filtered[["gene_id", "consensus"]],
                        on="gene_id",
                        validate="one_to_one",
                    )
                    .merge(
                        gene_table[["gene_id", "is_deg"]],
                        on="gene_id",
                        validate="one_to_one",
                    )
                    .groupby("consensus", as_index=False)
                    .agg(gene_count=("gene_id", "size"), deg_count=("is_deg", "sum"))
                    .rename(columns={"consensus": "filtered_flt_cluster"})
                    .sort_values("gene_count", ascending=False)
                )
                original_cluster12.to_csv(
                    output_dir / "original_flt_cluster12_filtered_distribution.tsv",
                    sep="\t",
                    index=False,
                )

    summary = {
        "deseq2_genes_tested_in_glare": int(tested["log2FoldChange"].notna().sum()),
        "deseq2_significant_genes": int(tested["is_deg"].sum()),
        "deseq2_up": int((tested["direction"] == "up").sum()),
        "deseq2_down": int((tested["direction"] == "down").sum()),
        "flt_gc_cluster_ari": float(cluster_ari),
        "flt_gc_cluster_nmi": float(cluster_nmi),
        "abs_lfc_vs_latent_shift_spearman_rho": float(correlation.statistic),
        "abs_lfc_vs_latent_shift_spearman_p": float(correlation.pvalue),
        "deg_latent_shift_mannwhitney_u": (
            float(shift_test.statistic) if shift_test else None
        ),
        "deg_latent_shift_mannwhitney_p": (
            float(shift_test.pvalue) if shift_test else None
        ),
        "latent_shift_deg_roc_auc": float(latent_auc),
        "deg_enriched_glare_clusters_fdr05": int(len(significant_clusters)),
        "significant_cluster_records": significant_clusters.to_dict(
            orient="records"
        ),
        "latent_alignment_note": (
            "GC latent coordinates were standardized and orthogonally aligned "
            "to FLT before calculating gene-wise distances; separate SAE spaces "
            "are otherwise not directly coordinate-comparable."
        ),
        "reactome": {
            "background_mouse_ensembl_genes": len(universe_ids),
            "deseq2_query_mouse_ensembl_genes": len(deg_ids),
            "glare_top_shift_query_mouse_ensembl_genes": len(shift_ids),
            "deseq2_significant_terms_fdr05": len(significant_deseq_terms),
            "glare_shift_significant_terms_fdr05": len(significant_glare_terms),
            "shared_significant_terms": len(
                significant_deseq_terms & significant_glare_terms
            ),
            "shared_terms": sorted(
                significant_deseq_terms & significant_glare_terms
            ),
            "top_deseq2_terms": (
                deseq_pathways.loc[
                    deseq_pathways["fdr_bh"] < 0.05, "term"
                ].head(10).tolist()
                if not deseq_pathways.empty
                else []
            ),
            "top_glare_shift_terms": (
                glare_pathways.loc[
                    glare_pathways["fdr_bh"] < 0.05, "term"
                ].head(10).tolist()
                if not glare_pathways.empty
                else []
            ),
        },
        "filtered_vs_original_glare": reference_comparison,
        "original_flt_cluster12": {
            "genes": int(original_cluster12["gene_count"].sum())
            if not original_cluster12.empty
            else 0,
            "filtered_clusters_spanned": int(len(original_cluster12)),
            "deseq2_significant_genes": int(original_cluster12["deg_count"].sum())
            if not original_cluster12.empty
            else 0,
        },
    }
    (output_dir / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    report = f"""# Filtered DESeq2 vs GLARE

- DESeq2 tested {summary["deseq2_genes_tested_in_glare"]:,} genes represented
  in GLARE and called {summary["deseq2_significant_genes"]:,} significant at
  adjusted p < 0.05 and absolute log2 fold change >= 1
  ({summary["deseq2_up"]} up, {summary["deseq2_down"]} down).
- {summary["deg_enriched_glare_clusters_fdr05"]} FLT/GC consensus clusters are
  enriched for DEGs by one-sided Fisher test with BH FDR < 0.05.
- Absolute DESeq2 effect and aligned GLARE latent shift have Spearman rho
  {summary["abs_lfc_vs_latent_shift_spearman_rho"]:.3f}
  (p={summary["abs_lfc_vs_latent_shift_spearman_p"]:.3g}).
- GLARE latent shift separates DEGs from non-DEGs with ROC AUC
  {summary["latent_shift_deg_roc_auc"]:.3f}; the one-sided Mann-Whitney p-value
  is {summary["deg_latent_shift_mannwhitney_p"]:.3g}.
- FLT and GC consensus partitions have ARI {cluster_ari:.3f} and NMI
  {cluster_nmi:.3f}.
- Reactome ORA found {len(significant_deseq_terms)} significant terms for
  mapped DESeq2 DEGs and {len(significant_glare_terms)} for an equal-sized set
  of top GLARE-shift genes; {len(significant_deseq_terms & significant_glare_terms)}
  terms are shared.
- The original muscle-heavy FLT cluster 12 now spans
  {len(original_cluster12)} filtered FLT clusters and contains
  {int(original_cluster12["deg_count"].sum()) if not original_cluster12.empty else 0}
  composition-adjusted DEGs.
- Filtered-vs-original clustering ARI is
  {reference_comparison.get("FLT", {}).get("ari", float("nan")):.3f} for FLT
  and {reference_comparison.get("GC", {}).get("ari", float("nan")):.3f} for GC.

DESeq2 and GLARE answer different questions. DESeq2 tests sample-level mean
expression changes after adjusting for age/collection stratum. GLARE groups
genes by nonlinear expression representations and can capture coordinated
structure even when individual genes do not pass a DEG threshold. Agreement is
therefore assessed through DEG enrichment within GLARE clusters and through
the relationship between DESeq2 effect size and aligned latent-space shift.
"""
    (output_dir / "COMPARISON_SUMMARY.md").write_text(report, encoding="utf-8")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir", default="outputs/glare_filtered_tms_liver_osd379"
    )
    parser.add_argument(
        "--reference-run",
        default="outputs/glare_paper_tms_liver_osd379",
        help="Original unfiltered GLARE run used for cluster-stability comparison.",
    )
    parser.add_argument(
        "--deseq-dir",
        default="outputs/glare_filtered_tms_liver_osd379/deseq2",
    )
    parser.add_argument(
        "--annotation",
        default="assets/osdr/GLDS-379_rna_seq_differential_expression_GLbulkRNAseq.csv",
    )
    parser.add_argument(
        "--reactome-gmt",
        default=(
            "src/expiMap_reproducibility/metadata/"
            "c2.cp.reactome.v4.0_mouseEID.gmt"
        ),
    )
    return parser.parse_args()


def main() -> None:
    print(json.dumps(run(parse_args()), indent=2))


if __name__ == "__main__":
    main()
