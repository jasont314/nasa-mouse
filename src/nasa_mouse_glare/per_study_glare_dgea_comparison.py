"""Compare per-study GLARE modules with matched per-study DESeq2 results."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import orthogonal_procrustes
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    roc_auc_score,
)
from sklearn.preprocessing import StandardScaler

from .aggregate_liver_analysis import load_gene_symbols
from .cluster_enrichment import bh_fdr
from .per_study_glare import (
    DEFAULT_OUTPUT_DIR as DEFAULT_GLARE_DIR,
    DEFAULT_STUDIES,
)
from .study_specific_pathway_recurrence import (
    infer_stratum,
    read_gmt,
    run_ora_enrichment,
    run_rank_sum_enrichment,
)


DEFAULT_COMPARISON_OUTPUT_DIR = (
    "outputs/glare_per_study_liver_noercc_12filter/dgea_comparison"
)
DEFAULT_OSDR_H5 = "assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5"
DEFAULT_REACTOME_GMT = (
    "src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt"
)


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def default_rscript() -> str:
    candidate = Path(sys.executable).with_name("Rscript")
    if candidate.exists():
        return str(candidate)
    return shutil.which("Rscript") or "Rscript"


def condition_for_location(location: str) -> str:
    if location == "FLT":
        return "flight"
    if location == "GC":
        return "ground"
    raise ValueError(f"Unsupported GLARE location: {location}")


def unique_sample_name(feature: str, seen: dict[str, int]) -> str:
    count = seen.get(feature, 0) + 1
    seen[feature] = count
    if count == 1:
        return feature
    return f"{feature}__dup{count}"


def export_matched_deseq2_inputs(
    glare_dir: Path,
    output_dir: Path,
    studies: list[str],
    symbols: dict[str, str],
) -> dict[str, str | int]:
    input_dir = output_dir / "deseq2_inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    genes_reference: np.ndarray | None = None
    count_blocks = []
    metadata_rows = []
    seen_samples: dict[str, int] = {}

    for accession in studies:
        run_dir = glare_dir / accession
        target_path = run_dir / "controlled_target.npz"
        retained_path = run_dir / "retained_profile_features.tsv"
        if not target_path.exists() or not retained_path.exists():
            raise SystemExit(f"Missing per-study GLARE run files for {accession}")

        with np.load(target_path) as target:
            genes = target["genes"].astype(str)
            if genes_reference is None:
                genes_reference = genes
            elif not np.array_equal(genes_reference, genes):
                raise SystemExit(f"Gene order differs for {accession}")

            retained = pd.read_csv(retained_path, sep="\t", keep_default_na=False)
            for location, matrix_key, feature_key in [
                ("FLT", "flt", "flt_features"),
                ("GC", "gc", "gc_features"),
            ]:
                matrix = np.rint(target[matrix_key]).astype(np.int64, copy=False)
                features = target[feature_key].astype(str)
                if matrix.shape[1] != len(features):
                    raise SystemExit(
                        f"{accession} {location}: feature count does not match matrix"
                    )
                count_blocks.append(matrix)
                retained_location = retained.loc[retained["location"].eq(location)].copy()
                used_indices: set[int] = set()
                for feature in features:
                    candidates = retained_location.index[
                        retained_location["feature"].astype(str).eq(feature)
                    ].tolist()
                    unused = [index for index in candidates if index not in used_indices]
                    if not unused:
                        raise SystemExit(
                            f"{accession} {location}: no retained metadata row for {feature}"
                        )
                    retained_index = unused[0]
                    used_indices.add(retained_index)
                    row = retained_location.loc[retained_index]
                    sample = unique_sample_name(str(feature), seen_samples)
                    profile = str(row.get("profile", feature))
                    metadata_rows.append(
                        {
                            "sample": sample,
                            "accession": accession,
                            "condition": condition_for_location(location),
                            "stratum": infer_stratum(profile, accession),
                            "profile": profile,
                            "original_feature": feature,
                            "source_profile_index": row.get("source_profile_index", ""),
                            "official_sample_name": row.get("official_sample_name", ""),
                            "official_material_type": row.get("official_material_type", ""),
                            "official_tissue": row.get("official_tissue", ""),
                            "library_selection": row.get("library_selection", ""),
                            "library_layout": row.get("library_layout", ""),
                            "sex": row.get("sex", ""),
                            "strain": row.get("strain", ""),
                            "ercc_status": ercc_status(profile),
                        }
                    )

    if genes_reference is None:
        raise SystemExit("No studies were selected")

    counts = np.concatenate(count_blocks, axis=1)
    samples = [row["sample"] for row in metadata_rows]
    count_table = pd.DataFrame(counts, index=genes_reference, columns=samples)
    count_table.index.name = "gene_id"
    counts_path = input_dir / "counts.tsv"
    count_table.to_csv(counts_path, sep="\t")

    metadata = pd.DataFrame(metadata_rows)
    metadata_path = input_dir / "sample_metadata.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)

    symbol_path = input_dir / "gene_symbols.tsv"
    pd.DataFrame(
        {
            "gene_id": genes_reference,
            "gene_symbol": [symbols.get(gene, "") for gene in genes_reference],
        }
    ).to_csv(symbol_path, sep="\t", index=False)

    counts_by_study = (
        metadata.groupby(["accession", "condition"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    counts_by_study.to_csv(input_dir / "study_condition_counts.tsv", sep="\t", index=False)
    ercc_counts = (
        metadata.groupby(["accession", "condition", "ercc_status"])
        .size()
        .reset_index(name="profiles")
    )
    ercc_counts.to_csv(input_dir / "study_ercc_counts.tsv", sep="\t", index=False)

    return {
        "counts": str(counts_path),
        "metadata": str(metadata_path),
        "gene_symbols": str(symbol_path),
        "study_condition_counts": str(input_dir / "study_condition_counts.tsv"),
        "study_ercc_counts": str(input_dir / "study_ercc_counts.tsv"),
        "genes": int(len(genes_reference)),
        "samples": int(len(samples)),
    }


def ercc_status(profile: str) -> str:
    if "_noERCC_" in profile or profile.endswith("_noERCC"):
        return "noERCC"
    if "_wERCC_" in profile or profile.endswith("_wERCC"):
        return "wERCC"
    return "not_annotated"


def run_deseq2(
    rscript: str,
    input_paths: dict[str, str | int],
    output_dir: Path,
    alpha: float,
    min_count: int,
    min_samples: int,
) -> None:
    script = Path(__file__).with_name("study_specific_deseq2.R")
    command = [
        rscript,
        str(script),
        "--counts",
        str(input_paths["counts"]),
        "--metadata",
        str(input_paths["metadata"]),
        "--gene-symbols",
        str(input_paths["gene_symbols"]),
        "--output-dir",
        str(output_dir / "deseq2"),
        "--alpha",
        str(alpha),
        "--min-count",
        str(min_count),
        "--min-samples",
        str(min_samples),
    ]
    subprocess.run(command, check=True)


def procrustes_shift(flt: pd.DataFrame, gc: pd.DataFrame) -> pd.DataFrame:
    shared = flt.index.intersection(gc.index)
    flt_values = flt.loc[shared].to_numpy(dtype=float)
    gc_values = gc.loc[shared].to_numpy(dtype=float)
    flt_scaled = StandardScaler().fit_transform(flt_values)
    gc_scaled = StandardScaler().fit_transform(gc_values)
    rotation, _ = orthogonal_procrustes(gc_scaled, flt_scaled)
    aligned_gc = gc_scaled @ rotation
    shift = np.linalg.norm(flt_scaled - aligned_gc, axis=1)
    return pd.DataFrame({"gene_id": shared.astype(str), "procrustes_latent_shift": shift})


def bh_adjust_nullable(values: list[float]) -> list[float]:
    array = np.asarray(values, dtype=float)
    adjusted = np.full(array.shape, np.nan, dtype=float)
    valid = np.isfinite(array)
    if valid.any():
        adjusted[valid] = bh_fdr(array[valid])
    return adjusted.tolist()


def top_gene_names(table: pd.DataFrame, max_genes: int = 10) -> str:
    if table.empty:
        return ""
    sorted_table = table.sort_values(
        ["padj", "pvalue", "abs_log2_fold_change"],
        ascending=[True, True, False],
        na_position="last",
    ).head(max_genes)
    names = []
    for row in sorted_table.itertuples():
        symbol = getattr(row, "gene_symbol", "")
        names.append(symbol if isinstance(symbol, str) and symbol else row.gene_id)
    return ",".join(names)


def cluster_dgea_summary(
    gene_table: pd.DataFrame,
    clusters: pd.DataFrame,
    accession: str,
    location: str,
    alpha: float,
    lfc_cutoff: float,
) -> pd.DataFrame:
    merged = clusters[["gene_id", "consensus"]].merge(
        gene_table,
        on="gene_id",
        how="left",
        validate="one_to_one",
    )
    tested = merged.loc[merged["tested_dgea"]].copy()
    total_sig = int(tested["significant_padj05"].sum())
    total_not_sig = int(len(tested) - total_sig)
    rows = []
    all_gene_counts = clusters.groupby("consensus")["gene_id"].size().to_dict()

    for cluster, group in tested.groupby("consensus", sort=True):
        sig = group.loc[group["significant_padj05"]].copy()
        sig_lfc = group.loc[group["significant_padj05_abs_lfc"]].copy()
        n_sig = int(len(sig))
        n_total = int(len(group))
        n_not_sig = n_total - n_sig
        if total_sig > 0 and n_total > 0:
            outside_sig = total_sig - n_sig
            outside_not_sig = total_not_sig - n_not_sig
            odds_ratio, p_value = fisher_exact(
                [[n_sig, n_not_sig], [outside_sig, outside_not_sig]],
                alternative="greater",
            )
        else:
            odds_ratio = np.nan
            p_value = np.nan
        rows.append(
            {
                "accession": accession,
                "location": location,
                "cluster": int(cluster),
                "cluster_genes": int(all_gene_counts.get(cluster, 0)),
                "tested_dgea_genes": n_total,
                "significant_padj05_genes": n_sig,
                "significant_padj05_abs_lfc_genes": int(len(sig_lfc)),
                "significant_fraction": n_sig / n_total if n_total else np.nan,
                "mean_log2fc_tested": float(group["log2FoldChange"].mean())
                if n_total
                else np.nan,
                "mean_log2fc_sig": float(sig["log2FoldChange"].mean())
                if n_sig
                else np.nan,
                "mean_abs_stat_tested": float(group["stat"].abs().mean())
                if n_total
                else np.nan,
                "mean_latent_shift_tested": float(
                    group["procrustes_latent_shift"].mean()
                )
                if n_total
                else np.nan,
                "fisher_odds_ratio": float(odds_ratio)
                if np.isfinite(odds_ratio)
                else np.nan,
                "fisher_p_value": float(p_value) if np.isfinite(p_value) else np.nan,
                "top_significant_genes": top_gene_names(sig),
            }
        )

    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fisher_fdr_bh"] = bh_adjust_nullable(result["fisher_p_value"].tolist())
    return result.sort_values(
        [
            "significant_padj05_genes",
            "fisher_fdr_bh",
            "tested_dgea_genes",
        ],
        ascending=[False, True, False],
        na_position="last",
    )


def study_gene_comparison(
    glare_dir: Path,
    accession: str,
    deseq: pd.DataFrame,
    alpha: float,
    lfc_cutoff: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    run_dir = glare_dir / accession
    flt_latent = pd.read_csv(run_dir / "FLT_gene_latent.tsv", sep="\t").set_index(
        "gene_id"
    )
    gc_latent = pd.read_csv(run_dir / "GC_gene_latent.tsv", sep="\t").set_index(
        "gene_id"
    )
    shifts = procrustes_shift(flt_latent, gc_latent)
    flt_clusters = pd.read_csv(
        run_dir / "clustering" / "FLT_gene_clusters.tsv", sep="\t"
    ).rename(columns={"consensus": "flt_cluster"})
    gc_clusters = pd.read_csv(
        run_dir / "clustering" / "GC_gene_clusters.tsv", sep="\t"
    ).rename(columns={"consensus": "gc_cluster"})
    clusters = flt_clusters[["gene_id", "flt_cluster"]].merge(
        gc_clusters[["gene_id", "gc_cluster"]],
        on="gene_id",
        how="inner",
        validate="one_to_one",
    )
    study_deseq = deseq.loc[deseq["accession"].eq(accession)].copy()
    study_deseq["gene_id"] = study_deseq["gene_id"].astype(str)

    gene_table = (
        clusters.merge(shifts, on="gene_id", how="left", validate="one_to_one")
        .merge(study_deseq, on="gene_id", how="left", validate="one_to_one")
        .copy()
    )
    gene_table["tested_dgea"] = gene_table["stat"].notna()
    gene_table["significant_padj05"] = (
        gene_table["significant_padj05"].fillna(False).astype(bool)
    )
    gene_table["abs_log2_fold_change"] = gene_table["log2FoldChange"].abs()
    gene_table["significant_padj05_abs_lfc"] = (
        gene_table["significant_padj05"]
        & (gene_table["abs_log2_fold_change"] >= lfc_cutoff)
    )
    gene_table["direction"] = gene_table["direction"].fillna("not_tested")

    cluster_tables = []
    cluster_tables.append(
        cluster_dgea_summary(
            gene_table,
            flt_clusters.rename(columns={"flt_cluster": "consensus"}),
            accession,
            "FLT",
            alpha,
            lfc_cutoff,
        )
    )
    cluster_tables.append(
        cluster_dgea_summary(
            gene_table,
            gc_clusters.rename(columns={"gc_cluster": "consensus"}),
            accession,
            "GC",
            alpha,
            lfc_cutoff,
        )
    )
    cluster_summary = pd.concat(cluster_tables, ignore_index=True)

    tested = gene_table.loc[gene_table["tested_dgea"]].copy()
    sig = tested.loc[tested["significant_padj05"]].copy()
    if len(tested) >= 3 and tested["abs_log2_fold_change"].nunique() > 1:
        lfc_corr = spearmanr(
            tested["abs_log2_fold_change"],
            tested["procrustes_latent_shift"],
            nan_policy="omit",
        )
        stat_corr = spearmanr(
            tested["stat"].abs(),
            tested["procrustes_latent_shift"],
            nan_policy="omit",
        )
    else:
        lfc_corr = stat_corr = None

    if tested["significant_padj05"].nunique() == 2:
        auc = float(
            roc_auc_score(
                tested["significant_padj05"],
                tested["procrustes_latent_shift"],
            )
        )
        shift_test = mannwhitneyu(
            tested.loc[tested["significant_padj05"], "procrustes_latent_shift"],
            tested.loc[~tested["significant_padj05"], "procrustes_latent_shift"],
            alternative="greater",
        )
        shift_p = float(shift_test.pvalue)
    else:
        auc = np.nan
        shift_p = np.nan

    partition_table = gene_table.dropna(subset=["flt_cluster", "gc_cluster"])
    cluster_ari = adjusted_rand_score(
        partition_table["flt_cluster"], partition_table["gc_cluster"]
    )
    cluster_nmi = normalized_mutual_info_score(
        partition_table["flt_cluster"], partition_table["gc_cluster"]
    )

    significant_clusters = cluster_summary.loc[
        cluster_summary["fisher_fdr_bh"].lt(alpha)
    ]
    summary = {
        "accession": accession,
        "genes_in_glare": int(len(gene_table)),
        "genes_tested_dgea": int(gene_table["tested_dgea"].sum()),
        "significant_padj05": int(gene_table["significant_padj05"].sum()),
        "significant_padj05_abs_lfc": int(
            gene_table["significant_padj05_abs_lfc"].sum()
        ),
        "significant_up": int((sig["log2FoldChange"] > 0).sum()),
        "significant_down": int((sig["log2FoldChange"] < 0).sum()),
        "deg_enriched_glare_clusters_fdr05": int(len(significant_clusters)),
        "abs_lfc_vs_latent_shift_spearman_rho": float(lfc_corr.statistic)
        if lfc_corr is not None
        else np.nan,
        "abs_lfc_vs_latent_shift_spearman_p": float(lfc_corr.pvalue)
        if lfc_corr is not None
        else np.nan,
        "abs_stat_vs_latent_shift_spearman_rho": float(stat_corr.statistic)
        if stat_corr is not None
        else np.nan,
        "abs_stat_vs_latent_shift_spearman_p": float(stat_corr.pvalue)
        if stat_corr is not None
        else np.nan,
        "latent_shift_deg_roc_auc": auc,
        "deg_latent_shift_mannwhitney_p": shift_p,
        "flt_gc_cluster_ari": float(cluster_ari),
        "flt_gc_cluster_nmi": float(cluster_nmi),
    }
    return gene_table, cluster_summary, summary


def summarize_glare_terms(glare_enrichment: pd.DataFrame, alpha: float) -> pd.DataFrame:
    if glare_enrichment.empty:
        return pd.DataFrame()
    significant = glare_enrichment.loc[glare_enrichment["fdr_bh"].lt(alpha)].copy()
    if significant.empty:
        return pd.DataFrame()
    rows = []
    for (accession, term), group in significant.groupby(["accession", "term"]):
        rows.append(
            {
                "accession": accession,
                "term": term,
                "glare_best_fdr_bh": float(group["fdr_bh"].min()),
                "glare_cluster_count": int(
                    len(group[["location", "cluster"]].drop_duplicates())
                ),
                "glare_locations": ",".join(sorted(group["location"].unique())),
                "glare_example_clusters": ";".join(
                    f"{row.location}{int(row.cluster)}"
                    for row in group.sort_values("fdr_bh").head(8).itertuples()
                ),
            }
        )
    return pd.DataFrame(rows)


def compare_pathways(
    deseq: pd.DataFrame,
    glare_dir: Path,
    output_dir: Path,
    studies: list[str],
    reactome_gmt: Path,
    alpha: float,
    min_pathway_size: int,
    max_pathway_size: int,
    min_overlap: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gene_sets = read_gmt(reactome_gmt)
    rank = run_rank_sum_enrichment(
        deseq,
        gene_sets,
        studies,
        min_pathway_size,
        max_pathway_size,
    )
    ora = run_ora_enrichment(deseq, gene_sets, studies, alpha, min_overlap)

    glare_enrichment_path = glare_dir / "glare_cluster_reactome_enrichment.tsv"
    if glare_enrichment_path.exists():
        glare_enrichment = pd.read_csv(glare_enrichment_path, sep="\t")
    else:
        glare_enrichment = pd.DataFrame()
    glare_terms = summarize_glare_terms(glare_enrichment, alpha)

    def overlap_with_glare(dgea_terms: pd.DataFrame) -> pd.DataFrame:
        if dgea_terms.empty or glare_terms.empty:
            return pd.DataFrame()
        return dgea_terms.loc[dgea_terms["fdr_bh"].lt(alpha)].merge(
            glare_terms,
            on=["accession", "term"],
            how="inner",
            validate="many_to_one",
        )

    rank_overlap = overlap_with_glare(rank)
    ora_overlap = overlap_with_glare(ora)

    recurring = pd.DataFrame()
    if not rank_overlap.empty:
        rows = []
        for (direction, term, clean_term), group in rank_overlap.groupby(
            ["direction", "term", "clean_term"]
        ):
            rows.append(
                {
                    "direction": direction,
                    "term": term,
                    "clean_term": clean_term,
                    "study_count": int(group["accession"].nunique()),
                    "accessions": ",".join(sorted(group["accession"].unique())),
                    "best_dgea_fdr_bh": float(group["fdr_bh"].min()),
                    "best_glare_fdr_bh": float(group["glare_best_fdr_bh"].min()),
                    "mean_wald_stat_shift_mean": float(
                        group["mean_wald_stat_shift"].mean()
                    ),
                    "example_glare_clusters": ";".join(
                        group.sort_values("glare_best_fdr_bh")[
                            "glare_example_clusters"
                        ]
                        .head(4)
                        .tolist()
                    ),
                }
            )
        recurring = pd.DataFrame(rows).sort_values(
            ["study_count", "best_dgea_fdr_bh", "best_glare_fdr_bh"],
            ascending=[False, True, True],
        )

    rank.to_csv(output_dir / "dgea_reactome_rank_sum_pathways.tsv", sep="\t", index=False)
    ora.to_csv(output_dir / "dgea_reactome_ora_pathways.tsv", sep="\t", index=False)
    glare_terms.to_csv(
        output_dir / "significant_glare_reactome_terms_by_study.tsv",
        sep="\t",
        index=False,
    )
    rank_overlap.to_csv(
        output_dir / "dgea_rank_pathway_glare_overlap.tsv", sep="\t", index=False
    )
    ora_overlap.to_csv(
        output_dir / "dgea_ora_pathway_glare_overlap.tsv", sep="\t", index=False
    )
    recurring.to_csv(
        output_dir / "recurring_dgea_glare_pathway_overlap.tsv",
        sep="\t",
        index=False,
    )
    return rank, ora, rank_overlap, recurring


def markdown_table(frame: pd.DataFrame, columns: list[str], max_rows: int = 16) -> list[str]:
    if frame.empty:
        return ["No rows."]
    display = frame.loc[:, [col for col in columns if col in frame.columns]].head(max_rows)
    lines = [
        "| " + " | ".join(display.columns) + " |",
        "| " + " | ".join(["---"] * len(display.columns)) + " |",
    ]
    for _, row in display.iterrows():
        values = []
        for value in row.tolist():
            if isinstance(value, float):
                values.append(f"{value:.3g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def write_report(
    output_dir: Path,
    summary: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    recurring_overlap: pd.DataFrame,
    rank_overlap: pd.DataFrame,
    input_paths: dict[str, str | int],
    alpha: float,
    lfc_cutoff: float,
) -> None:
    enriched_clusters = cluster_summary.loc[cluster_summary["fisher_fdr_bh"].lt(alpha)]
    lines = [
        "# Per-Study GLARE vs DGEA",
        "",
        "Each study was analyzed separately. DESeq2 inputs were exported from the",
        "same per-study GLARE target matrices, so this comparison uses the 12",
        "muscle-outlier filter and the same ERCC/noERCC selections as the GLARE runs.",
        "",
        "## Inputs",
        "",
        f"- Counts: `{input_paths['counts']}`",
        f"- Metadata: `{input_paths['metadata']}`",
        f"- Gene symbols: `{input_paths['gene_symbols']}`",
        f"- DESeq2 alpha: {alpha}",
        f"- Cluster DEG strict flag: adjusted p < {alpha} and abs(log2FC) >= {lfc_cutoff}",
        "",
        "## Study Summary",
        "",
        *markdown_table(
            summary,
            [
                "accession",
                "n_flight",
                "n_ground",
                "deseq2_design",
                "genes_tested_dgea",
                "significant_padj05",
                "significant_padj05_abs_lfc",
                "deg_enriched_glare_clusters_fdr05",
                "abs_stat_vs_latent_shift_spearman_rho",
                "latent_shift_deg_roc_auc",
                "rank_pathway_glare_overlaps",
            ],
        ),
        "",
        "## DEG-Enriched GLARE Clusters",
        "",
    ]
    if enriched_clusters.empty:
        lines.append("No per-study GLARE clusters were significantly enriched for DESeq2 FDR genes.")
    else:
        lines.extend(
            markdown_table(
                enriched_clusters.sort_values(
                    ["fisher_fdr_bh", "significant_padj05_genes"],
                    ascending=[True, False],
                ),
                [
                    "accession",
                    "location",
                    "cluster",
                    "tested_dgea_genes",
                    "significant_padj05_genes",
                    "significant_padj05_abs_lfc_genes",
                    "fisher_odds_ratio",
                    "fisher_fdr_bh",
                    "top_significant_genes",
                ],
                max_rows=24,
            )
        )
    lines.extend(
        [
            "",
            "## Recurring DGEA/GLARE Pathway Concordance",
            "",
        ]
    )
    if recurring_overlap.empty:
        lines.append(
            "No Reactome pathway was significant in both DESeq2 rank-sum testing and GLARE cluster enrichment in two or more studies."
        )
    else:
        lines.extend(
            markdown_table(
                recurring_overlap,
                [
                    "direction",
                    "clean_term",
                    "study_count",
                    "accessions",
                    "best_dgea_fdr_bh",
                    "best_glare_fdr_bh",
                    "mean_wald_stat_shift_mean",
                    "example_glare_clusters",
                ],
                max_rows=24,
            )
        )
    lines.extend(
        [
            "",
            "## Per-Study DGEA/GLARE Pathway Overlaps",
            "",
        ]
    )
    if rank_overlap.empty:
        lines.append("No same-study rank-sum pathway overlaps passed FDR.")
    else:
        lines.extend(
            markdown_table(
                rank_overlap.sort_values(["accession", "fdr_bh", "glare_best_fdr_bh"]),
                [
                    "accession",
                    "direction",
                    "clean_term",
                    "fdr_bh",
                    "mean_wald_stat_shift",
                    "glare_best_fdr_bh",
                    "glare_locations",
                    "glare_example_clusters",
                ],
                max_rows=40,
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "DESeq2 and GLARE are being used for different evidence layers.",
            "DESeq2 is the per-study sample-level FLT-vs-GC statistical test.",
            "GLARE is treated as module discovery: it is most useful where its",
            "clusters are enriched for DESeq2 genes or recover Reactome pathways",
            "that are also shifted in per-study DESeq2 rankings.",
            "",
            "Small studies can show module structure but have limited DESeq2 power;",
            "interpret `OSD-48`, `OSD-137`, and `OSD-168` as support/sensitivity",
            "unless their signals recur in the larger studies.",
        ]
    )
    (output_dir / "GLARE_VS_DGEA_PER_STUDY.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glare-dir", default=DEFAULT_GLARE_DIR)
    parser.add_argument("--output-dir", default=DEFAULT_COMPARISON_OUTPUT_DIR)
    parser.add_argument("--osdr-h5", default=DEFAULT_OSDR_H5)
    parser.add_argument("--reactome-gmt", default=DEFAULT_REACTOME_GMT)
    parser.add_argument("--studies", nargs="+", default=DEFAULT_STUDIES)
    parser.add_argument("--rscript", default=default_rscript())
    parser.add_argument("--skip-deseq2", action="store_true")
    parser.add_argument("--export-only", action="store_true")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--lfc-cutoff", type=float, default=1.0)
    parser.add_argument("--min-count", type=int, default=10)
    parser.add_argument("--min-samples", type=int, default=3)
    parser.add_argument("--min-overlap", type=int, default=3)
    parser.add_argument("--min-pathway-size", type=int, default=10)
    parser.add_argument("--max-pathway-size", type=int, default=500)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    glare_dir = Path(args.glare_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log("Loading gene symbols")
    symbols = load_gene_symbols(args.osdr_h5)
    log("Exporting matched DESeq2 inputs from per-study GLARE runs")
    input_paths = export_matched_deseq2_inputs(
        glare_dir,
        output_dir,
        args.studies,
        symbols,
    )
    if args.export_only:
        log(f"Exported DESeq2 inputs to {output_dir / 'deseq2_inputs'}")
        return

    if not args.skip_deseq2:
        log("Running per-study DESeq2")
        run_deseq2(
            args.rscript,
            input_paths,
            output_dir,
            args.alpha,
            args.min_count,
            args.min_samples,
        )

    deseq_path = output_dir / "deseq2" / "per_study_deseq2.tsv"
    deseq_summary_path = output_dir / "deseq2" / "study_deseq2_summary.tsv"
    if not deseq_path.exists():
        raise SystemExit(f"Missing DESeq2 output: {deseq_path}")
    deseq = pd.read_csv(deseq_path, sep="\t")
    deseq["gene_id"] = deseq["gene_id"].astype(str)

    log("Comparing DESeq2 genes with per-study GLARE clusters")
    gene_tables = []
    cluster_tables = []
    summaries = []
    for accession in args.studies:
        gene_table, cluster_summary, study_summary = study_gene_comparison(
            glare_dir,
            accession,
            deseq,
            args.alpha,
            args.lfc_cutoff,
        )
        gene_table.to_csv(
            output_dir / f"{accession}_gene_level_glare_dgea.tsv",
            sep="\t",
            index=False,
        )
        cluster_summary.to_csv(
            output_dir / f"{accession}_cluster_dgea_enrichment.tsv",
            sep="\t",
            index=False,
        )
        gene_tables.append(gene_table.assign(accession_for_table=accession))
        cluster_tables.append(cluster_summary)
        summaries.append(study_summary)

    all_clusters = pd.concat(cluster_tables, ignore_index=True)
    all_clusters.to_csv(
        output_dir / "cluster_dgea_enrichment_all_studies.tsv",
        sep="\t",
        index=False,
    )

    log("Comparing DESeq2 Reactome pathways with GLARE cluster Reactome pathways")
    rank, ora, rank_overlap, recurring_overlap = compare_pathways(
        deseq,
        glare_dir,
        output_dir,
        args.studies,
        Path(args.reactome_gmt),
        args.alpha,
        args.min_pathway_size,
        args.max_pathway_size,
        args.min_overlap,
    )

    summary = pd.DataFrame(summaries)
    if deseq_summary_path.exists():
        deseq_summary = pd.read_csv(deseq_summary_path, sep="\t")
        summary = summary.merge(
            deseq_summary[
                [
                    "accession",
                    "n_flight",
                    "n_ground",
                    "genes_tested",
                    "design",
                    "dispersion_fit",
                ]
            ].rename(
                columns={
                    "design": "deseq2_design",
                    "dispersion_fit": "deseq2_dispersion_fit",
                }
            ),
            on="accession",
            how="left",
            validate="one_to_one",
        )
    else:
        summary["n_flight"] = np.nan
        summary["n_ground"] = np.nan
        summary["deseq2_design"] = ""
        summary["deseq2_dispersion_fit"] = ""
    if rank_overlap.empty:
        summary["rank_pathway_glare_overlaps"] = 0
    else:
        overlap_counts = rank_overlap.groupby("accession").size()
        summary["rank_pathway_glare_overlaps"] = (
            summary["accession"].map(overlap_counts).fillna(0).astype(int)
        )
    summary = summary[
        [
            "accession",
            "n_flight",
            "n_ground",
            "deseq2_design",
            "deseq2_dispersion_fit",
            "genes_in_glare",
            "genes_tested_dgea",
            "significant_padj05",
            "significant_padj05_abs_lfc",
            "significant_up",
            "significant_down",
            "deg_enriched_glare_clusters_fdr05",
            "abs_lfc_vs_latent_shift_spearman_rho",
            "abs_lfc_vs_latent_shift_spearman_p",
            "abs_stat_vs_latent_shift_spearman_rho",
            "abs_stat_vs_latent_shift_spearman_p",
            "latent_shift_deg_roc_auc",
            "deg_latent_shift_mannwhitney_p",
            "flt_gc_cluster_ari",
            "flt_gc_cluster_nmi",
            "rank_pathway_glare_overlaps",
        ]
    ]
    summary.to_csv(output_dir / "per_study_glare_dgea_summary.tsv", sep="\t", index=False)
    (output_dir / "per_study_glare_dgea_summary.json").write_text(
        json.dumps(summary.to_dict(orient="records"), indent=2) + "\n",
        encoding="utf-8",
    )
    write_report(
        output_dir,
        summary,
        all_clusters,
        recurring_overlap,
        rank_overlap,
        input_paths,
        args.alpha,
        args.lfc_cutoff,
    )
    log(f"Saved GLARE-vs-DGEA comparison to {output_dir}")


if __name__ == "__main__":
    main()
