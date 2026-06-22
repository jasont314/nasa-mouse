"""Post-analysis for aggregated OSDR liver FLT/GC GLARE results."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import orthogonal_procrustes
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.preprocessing import StandardScaler

from .io import require_import


DEFAULT_RUN_DIR = "outputs/glare_tms_liver_aggregated_osdr_flt_gc"
DEFAULT_OSDR_H5 = "assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5"
BATCH_VARIABLES = [
    "h5_accession",
    "mission",
    "library_selection",
    "library_layout",
    "sex",
    "strain",
    "sequencing_instrument",
]


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


def decode_array(values) -> list[str]:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8", "replace"))
        else:
            decoded.append(str(value))
    return decoded


def load_gene_symbols(osdr_h5: str | Path) -> dict[str, str]:
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")
    with h5py.File(osdr_h5, "r") as handle:
        genes = decode_array(handle["/meta/genes/ensembl_gene"][:])
        symbols = decode_array(handle["/meta/genes/symbol"][:])
    return dict(zip(genes, symbols))


def load_run(run_dir: Path) -> dict:
    target = np.load(run_dir / "controlled_target.npz")
    retained = pd.read_csv(run_dir / "retained_profile_features.tsv", sep="\t")
    clusters = {
        location: pd.read_csv(
            run_dir / "clustering" / f"{location}_gene_clusters.tsv", sep="\t"
        )
        for location in ("FLT", "GC")
    }
    latent = {
        location: pd.read_csv(
            run_dir / f"{location}_gene_latent.tsv", sep="\t"
        ).set_index("gene_id")
        for location in ("FLT", "GC")
    }
    return {
        "target": target,
        "retained": retained,
        "clusters": clusters,
        "latent": latent,
    }


def one_way_eta(values: np.ndarray, groups: pd.Series) -> tuple[float, float, int]:
    frame = pd.DataFrame({"value": values, "group": groups.astype(str).to_numpy()})
    frame = frame.loc[frame["group"].str.len() > 0].dropna()
    counts = frame["group"].value_counts()
    valid_groups = counts[counts >= 2].index.tolist()
    frame = frame.loc[frame["group"].isin(valid_groups)]
    if frame["group"].nunique() < 2:
        return np.nan, np.nan, int(frame["group"].nunique())
    group_arrays = [
        group["value"].to_numpy(dtype=float)
        for _, group in frame.groupby("group", sort=False)
    ]
    grand = frame["value"].mean()
    ss_total = float(((frame["value"] - grand) ** 2).sum())
    if ss_total <= 0:
        return 0.0, 1.0, int(frame["group"].nunique())
    ss_between = 0.0
    for group_values in group_arrays:
        ss_between += len(group_values) * float((group_values.mean() - grand) ** 2)
    eta_squared = ss_between / ss_total
    try:
        p_value = float(stats.f_oneway(*group_arrays).pvalue)
    except Exception:
        p_value = np.nan
    return eta_squared, p_value, int(frame["group"].nunique())


def cluster_sample_matrix(
    expression: np.ndarray,
    genes: np.ndarray,
    cluster_table: pd.DataFrame,
    metadata: pd.DataFrame,
    location: str,
) -> pd.DataFrame:
    gene_to_index = {gene: index for index, gene in enumerate(genes)}
    rows = []
    for cluster_id, cluster_rows in cluster_table.groupby("consensus", sort=True):
        indices = [
            gene_to_index[gene]
            for gene in cluster_rows["gene_id"].astype(str)
            if gene in gene_to_index
        ]
        if not indices:
            continue
        values = np.log2(expression[indices, :] + 1.0).mean(axis=0)
        for profile_row, value in zip(metadata.itertuples(index=False), values):
            rows.append(
                {
                    "location": location,
                    "cluster": int(cluster_id),
                    "n_genes": int(len(indices)),
                    "feature": profile_row.feature,
                    "profile": profile_row.profile,
                    "cluster_mean_log2_expression": float(value),
                    "h5_accession": profile_row.h5_accession,
                    "mission": getattr(
                        profile_row,
                        "mission",
                        getattr(profile_row, "project_identifier", ""),
                    ),
                    "library_selection": getattr(profile_row, "library_selection", ""),
                    "library_layout": getattr(profile_row, "library_layout", ""),
                    "sex": getattr(profile_row, "sex", ""),
                    "strain": getattr(profile_row, "strain", ""),
                    "sequencing_instrument": getattr(
                        profile_row, "sequencing_instrument", ""
                    ),
                }
            )
    return pd.DataFrame(rows)


def run_batch_qc(run: dict, output_dir: Path) -> dict:
    target = run["target"]
    genes = target["genes"].astype(str)
    retained = run["retained"]
    sample_frames = []
    qc_rows = []

    for location, matrix_key in (("FLT", "flt"), ("GC", "gc")):
        metadata = retained.loc[retained["location"].eq(location)].reset_index(drop=True)
        if "mission" not in metadata.columns:
            metadata["mission"] = metadata.get("project_identifier", "")
        expression = target[matrix_key]
        sample_expression = cluster_sample_matrix(
            expression,
            genes,
            run["clusters"][location],
            metadata,
            location,
        )
        sample_frames.append(sample_expression)
        for (cluster_id, n_genes), cluster_frame in sample_expression.groupby(
            ["cluster", "n_genes"], sort=True
        ):
            values = cluster_frame["cluster_mean_log2_expression"].to_numpy(dtype=float)
            for variable in BATCH_VARIABLES:
                eta, p_value, n_groups = one_way_eta(values, cluster_frame[variable])
                qc_rows.append(
                    {
                        "location": location,
                        "cluster": int(cluster_id),
                        "n_genes": int(n_genes),
                        "variable": variable,
                        "n_groups": n_groups,
                        "eta_squared": eta,
                        "p_value": p_value,
                    }
                )

    sample_table = pd.concat(sample_frames, ignore_index=True)
    qc = pd.DataFrame(qc_rows)
    qc["fdr_bh"] = qc.groupby(["location", "variable"])["p_value"].transform(
        lambda x: bh_adjust(x.to_numpy(dtype=float))
    )
    qc["strong_batch_driver"] = qc["eta_squared"].ge(0.25) & qc["fdr_bh"].lt(0.05)

    summary = (
        qc.groupby(["location", "variable"])
        .agg(
            clusters_tested=("cluster", "nunique"),
            median_eta_squared=("eta_squared", "median"),
            max_eta_squared=("eta_squared", "max"),
            strong_driver_clusters=("strong_batch_driver", "sum"),
        )
        .reset_index()
        .sort_values(["location", "max_eta_squared"], ascending=[True, False])
    )
    top = (
        qc.sort_values(["eta_squared", "fdr_bh"], ascending=[False, True])
        .groupby(["location", "cluster"], as_index=False)
        .head(1)
        .sort_values(["location", "eta_squared"], ascending=[True, False])
    )

    sample_table.to_csv(output_dir / "cluster_sample_expression.tsv", sep="\t", index=False)
    qc.to_csv(output_dir / "batch_driver_cluster_qc.tsv", sep="\t", index=False)
    summary.to_csv(output_dir / "batch_driver_summary.tsv", sep="\t", index=False)
    top.to_csv(output_dir / "top_batch_driver_by_cluster.tsv", sep="\t", index=False)

    return {
        "batch_qc_rows": int(len(qc)),
        "strong_driver_clusters": int(qc["strong_batch_driver"].sum()),
        "summary": summary.to_dict(orient="records"),
    }


def export_deseq2_inputs(run: dict, symbols: dict[str, str], output_dir: Path) -> dict:
    input_dir = output_dir / "deseq2_inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    target = run["target"]
    genes = target["genes"].astype(str)
    flt_features = target["flt_features"].astype(str)
    gc_features = target["gc_features"].astype(str)
    counts = np.concatenate([target["flt"], target["gc"]], axis=1)
    counts = np.rint(counts).astype(np.int64, copy=False)
    original_columns = list(flt_features) + list(gc_features)
    seen: dict[str, int] = {}
    count_columns = []
    for column in original_columns:
        duplicate_index = seen.get(column, 0)
        seen[column] = duplicate_index + 1
        if duplicate_index == 0:
            count_columns.append(column)
        else:
            count_columns.append(f"{column}__dup{duplicate_index + 1}")
    count_table = pd.DataFrame(counts, index=genes, columns=count_columns)
    count_table.index.name = "gene_id"
    count_path = input_dir / "counts.tsv"
    count_table.to_csv(count_path, sep="\t")

    retained = run["retained"].copy()
    if len(retained) != len(count_columns):
        raise ValueError(
            "Retained metadata row count does not match FLT+GC count columns: "
            f"{len(retained)} metadata rows vs {len(count_columns)} count columns"
        )
    metadata = retained[
        [
            "feature",
            "location",
            "h5_accession",
            "project_identifier",
            "library_selection",
            "library_layout",
            "sex",
            "strain",
            "sequencing_instrument",
        ]
    ].copy()
    metadata.insert(0, "sample", count_columns)
    metadata = metadata.rename(
        columns={
            "feature": "original_feature",
            "location": "condition",
            "h5_accession": "accession",
            "project_identifier": "mission",
        }
    )
    metadata_path = input_dir / "sample_metadata.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)

    symbol_path = input_dir / "gene_symbols.tsv"
    pd.DataFrame(
        {
            "gene_id": genes,
            "gene_symbol": [symbols.get(gene, "") for gene in genes],
        }
    ).to_csv(symbol_path, sep="\t", index=False)

    return {
        "counts": str(count_path),
        "metadata": str(metadata_path),
        "gene_symbols": str(symbol_path),
        "genes": int(len(genes)),
        "samples": int(len(count_columns)),
    }


def procrustes_latent_shift(flt_latent: pd.DataFrame, gc_latent: pd.DataFrame) -> pd.DataFrame:
    shared = flt_latent.index.intersection(gc_latent.index)
    flt = flt_latent.loc[shared].to_numpy(dtype=float)
    gc = gc_latent.loc[shared].to_numpy(dtype=float)
    flt_scaled = StandardScaler().fit_transform(flt)
    gc_scaled = StandardScaler().fit_transform(gc)
    rotation, _ = orthogonal_procrustes(gc_scaled, flt_scaled)
    aligned_gc = gc_scaled @ rotation
    shift = np.linalg.norm(flt_scaled - aligned_gc, axis=1)
    return pd.DataFrame({"gene_id": shared, "procrustes_latent_shift": shift})


def cluster_entropy(counts: np.ndarray) -> float:
    total = counts.sum()
    if total <= 0:
        return np.nan
    probs = counts[counts > 0] / total
    return float(-(probs * np.log2(probs)).sum())


def run_cluster_comparison(run: dict, output_dir: Path) -> dict:
    flt = run["clusters"]["FLT"][["gene_id", "consensus"]].rename(
        columns={"consensus": "flt_cluster"}
    )
    gc = run["clusters"]["GC"][["gene_id", "consensus"]].rename(
        columns={"consensus": "gc_cluster"}
    )
    merged = flt.merge(gc, on="gene_id", how="inner")
    ari = float(adjusted_rand_score(merged["flt_cluster"], merged["gc_cluster"]))
    nmi = float(normalized_mutual_info_score(merged["flt_cluster"], merged["gc_cluster"]))
    contingency = pd.crosstab(merged["flt_cluster"], merged["gc_cluster"])
    contingency.to_csv(output_dir / "flt_gc_cluster_contingency.tsv", sep="\t")

    flt_rows = []
    for cluster_id, group in merged.groupby("flt_cluster", sort=True):
        counts = group["gc_cluster"].value_counts()
        flt_rows.append(
            {
                "flt_cluster": int(cluster_id),
                "gene_count": int(len(group)),
                "gc_clusters_spanned": int(counts.size),
                "top_gc_cluster": int(counts.index[0]),
                "top_gc_overlap": int(counts.iloc[0]),
                "top_gc_overlap_fraction": float(counts.iloc[0] / len(group)),
                "gc_cluster_entropy": cluster_entropy(counts.to_numpy()),
            }
        )
    gc_rows = []
    for cluster_id, group in merged.groupby("gc_cluster", sort=True):
        counts = group["flt_cluster"].value_counts()
        gc_rows.append(
            {
                "gc_cluster": int(cluster_id),
                "gene_count": int(len(group)),
                "flt_clusters_spanned": int(counts.size),
                "top_flt_cluster": int(counts.index[0]),
                "top_flt_overlap": int(counts.iloc[0]),
                "top_flt_overlap_fraction": float(counts.iloc[0] / len(group)),
                "flt_cluster_entropy": cluster_entropy(counts.to_numpy()),
            }
        )

    shift = procrustes_latent_shift(run["latent"]["FLT"], run["latent"]["GC"])
    gene_level = merged.merge(shift, on="gene_id", how="left")
    flt_shift = (
        gene_level.groupby("flt_cluster")
        .agg(
            gene_count=("gene_id", "count"),
            mean_latent_shift=("procrustes_latent_shift", "mean"),
            median_latent_shift=("procrustes_latent_shift", "median"),
            p90_latent_shift=("procrustes_latent_shift", lambda x: np.nanpercentile(x, 90)),
        )
        .reset_index()
        .sort_values("mean_latent_shift", ascending=False)
    )
    gc_shift = (
        gene_level.groupby("gc_cluster")
        .agg(
            gene_count=("gene_id", "count"),
            mean_latent_shift=("procrustes_latent_shift", "mean"),
            median_latent_shift=("procrustes_latent_shift", "median"),
            p90_latent_shift=("procrustes_latent_shift", lambda x: np.nanpercentile(x, 90)),
        )
        .reset_index()
        .sort_values("mean_latent_shift", ascending=False)
    )

    pd.DataFrame(flt_rows).to_csv(
        output_dir / "flt_cluster_gc_overlap_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(gc_rows).to_csv(
        output_dir / "gc_cluster_flt_overlap_summary.tsv", sep="\t", index=False
    )
    gene_level.to_csv(output_dir / "gene_cluster_comparison.tsv", sep="\t", index=False)
    flt_shift.to_csv(output_dir / "flt_cluster_latent_shift_summary.tsv", sep="\t", index=False)
    gc_shift.to_csv(output_dir / "gc_cluster_latent_shift_summary.tsv", sep="\t", index=False)

    return {
        "genes_compared": int(len(gene_level)),
        "adjusted_rand_index": ari,
        "normalized_mutual_information": nmi,
        "mean_latent_shift": float(gene_level["procrustes_latent_shift"].mean()),
        "median_latent_shift": float(gene_level["procrustes_latent_shift"].median()),
        "top_flt_shift_clusters": flt_shift.head(5).to_dict(orient="records"),
        "top_gc_shift_clusters": gc_shift.head(5).to_dict(orient="records"),
    }


def welch_p_value(effect: np.ndarray, se: np.ndarray, vf: np.ndarray, vg: np.ndarray, nf: int, ng: int) -> np.ndarray:
    denom = (vf / nf + vg / ng) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        df = denom / ((vf / nf) ** 2 / max(nf - 1, 1) + (vg / ng) ** 2 / max(ng - 1, 1))
        t_stat = effect / se
        p_value = 2.0 * stats.t.sf(np.abs(t_stat), df)
    p_value[~np.isfinite(p_value)] = np.nan
    return p_value


def run_meta_dgea(run: dict, symbols: dict[str, str], output_dir: Path) -> dict:
    target = run["target"]
    genes = target["genes"].astype(str)
    flt = np.log2(target["flt"].astype(float) + 1.0)
    gc = np.log2(target["gc"].astype(float) + 1.0)
    retained = run["retained"]
    flt_meta = retained.loc[retained["location"].eq("FLT")].reset_index(drop=True)
    gc_meta = retained.loc[retained["location"].eq("GC")].reset_index(drop=True)
    accessions = sorted(set(flt_meta["h5_accession"]) & set(gc_meta["h5_accession"]))

    per_study_rows = []
    effect_stack = []
    se_stack = []
    accession_stack = []
    for accession in accessions:
        flt_idx = flt_meta.index[flt_meta["h5_accession"].eq(accession)].to_numpy()
        gc_idx = gc_meta.index[gc_meta["h5_accession"].eq(accession)].to_numpy()
        if len(flt_idx) < 2 or len(gc_idx) < 2:
            continue
        xf = flt[:, flt_idx]
        xg = gc[:, gc_idx]
        mean_f = xf.mean(axis=1)
        mean_g = xg.mean(axis=1)
        var_f = xf.var(axis=1, ddof=1)
        var_g = xg.var(axis=1, ddof=1)
        effect = mean_f - mean_g
        se = np.sqrt(var_f / len(flt_idx) + var_g / len(gc_idx))
        p_value = welch_p_value(effect, se, var_f, var_g, len(flt_idx), len(gc_idx))
        fdr = bh_adjust(p_value)
        for gene, gene_effect, gene_se, gene_p, gene_fdr in zip(
            genes, effect, se, p_value, fdr
        ):
            per_study_rows.append(
                {
                    "gene_id": gene,
                    "gene_symbol": symbols.get(gene, ""),
                    "accession": accession,
                    "n_flt": int(len(flt_idx)),
                    "n_gc": int(len(gc_idx)),
                    "log2_fold_change": float(gene_effect),
                    "se": float(gene_se) if np.isfinite(gene_se) else np.nan,
                    "p_value": float(gene_p) if np.isfinite(gene_p) else np.nan,
                    "fdr_bh": float(gene_fdr) if np.isfinite(gene_fdr) else np.nan,
                    "direction": (
                        "up" if gene_effect > 0 else "down" if gene_effect < 0 else "flat"
                    ),
                }
            )
        effect_stack.append(effect)
        se_stack.append(se)
        accession_stack.append(accession)

    effects = np.vstack(effect_stack)
    ses = np.vstack(se_stack)
    valid = np.isfinite(effects) & np.isfinite(ses) & (ses > 0)
    weights = np.zeros_like(ses, dtype=float)
    weights[valid] = 1.0 / np.square(ses[valid])
    weight_sum = weights.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        meta_effect = (weights * np.where(valid, effects, 0.0)).sum(axis=0) / weight_sum
        meta_se = np.sqrt(1.0 / weight_sum)
        meta_z = meta_effect / meta_se
        meta_p = 2.0 * stats.norm.sf(np.abs(meta_z))
    meta_effect[~np.isfinite(meta_effect)] = np.nan
    meta_se[~np.isfinite(meta_se)] = np.nan
    meta_z[~np.isfinite(meta_z)] = np.nan
    meta_p[~np.isfinite(meta_p)] = np.nan

    q = np.where(valid, weights * np.square(effects - meta_effect), 0.0).sum(axis=0)
    df = valid.sum(axis=0) - 1
    with np.errstate(divide="ignore", invalid="ignore"):
        i2 = np.maximum(0.0, (q - df) / q) * 100.0
        heterogeneity_p = stats.chi2.sf(q, df)
    i2[~np.isfinite(i2)] = np.nan
    heterogeneity_p[~np.isfinite(heterogeneity_p)] = np.nan

    direction_up = (effects > 0).sum(axis=0)
    direction_down = (effects < 0).sum(axis=0)
    studies_tested = valid.sum(axis=0)
    meta_fdr = bh_adjust(meta_p)
    meta = pd.DataFrame(
        {
            "gene_id": genes,
            "gene_symbol": [symbols.get(gene, "") for gene in genes],
            "studies_tested": studies_tested.astype(int),
            "meta_log2_fold_change": meta_effect,
            "meta_se": meta_se,
            "meta_z": meta_z,
            "meta_p_value": meta_p,
            "meta_fdr_bh": meta_fdr,
            "direction_up_studies": direction_up.astype(int),
            "direction_down_studies": direction_down.astype(int),
            "q_heterogeneity": q,
            "i2_percent": i2,
            "heterogeneity_p_value": heterogeneity_p,
        }
    )
    meta["significant_fdr05_abs_log2fc1"] = meta["meta_fdr_bh"].lt(0.05) & meta[
        "meta_log2_fold_change"
    ].abs().ge(1.0)
    meta["consistent_direction"] = np.where(
        meta["direction_up_studies"].eq(meta["studies_tested"]),
        "all_up",
        np.where(
            meta["direction_down_studies"].eq(meta["studies_tested"]),
            "all_down",
            "mixed",
        ),
    )
    meta = meta.sort_values(["meta_fdr_bh", "meta_p_value"], na_position="last")

    per_study = pd.DataFrame(per_study_rows)
    per_study.to_csv(output_dir / "per_study_dgea.tsv", sep="\t", index=False)
    meta.to_csv(output_dir / "meta_dgea.tsv", sep="\t", index=False)
    meta.head(200).to_csv(output_dir / "top_meta_dgea_genes.tsv", sep="\t", index=False)

    sig = meta.loc[meta["significant_fdr05_abs_log2fc1"]]
    return {
        "accessions_tested": accession_stack,
        "genes_tested": int(len(meta)),
        "significant_fdr05_abs_log2fc1": int(len(sig)),
        "significant_up": int((sig["meta_log2_fold_change"] > 0).sum()),
        "significant_down": int((sig["meta_log2_fold_change"] < 0).sum()),
        "top_genes": meta.head(20)[
            [
                "gene_id",
                "gene_symbol",
                "meta_log2_fold_change",
                "meta_fdr_bh",
                "consistent_direction",
            ]
        ].to_dict(orient="records"),
    }


def write_report(output_dir: Path, summary: dict) -> None:
    batch_summary = pd.DataFrame(summary["batch_qc"]["summary"])
    batch_summary = batch_summary.sort_values(["location", "variable"])
    batch_lines = batch_summary.to_csv(sep="\t", index=False).strip().splitlines()
    text = [
        "# Aggregated Liver Post-Analysis",
        "",
        "## Study/Batch QC",
        "",
        "Cluster mean expression was tested against accession, mission, library,",
        "sex, strain, and sequencing instrument within FLT and GC separately.",
        "`mission` is the OSDR project identifier field.",
        "",
        "```tsv",
        *batch_lines,
        "```",
        "",
        "## FLT vs GC Cluster Structure",
        "",
        f"- Genes compared: {summary['cluster_comparison']['genes_compared']:,}",
        (
            "- Adjusted Rand index: "
            f"{summary['cluster_comparison']['adjusted_rand_index']:.4f}"
        ),
        (
            "- Normalized mutual information: "
            f"{summary['cluster_comparison']['normalized_mutual_information']:.4f}"
        ),
        (
            "- Median Procrustes latent shift: "
            f"{summary['cluster_comparison']['median_latent_shift']:.4f}"
        ),
        "",
        "## Exploratory Meta-DGEA",
        "",
        "This is not DESeq2 from raw counts. It is an exploratory fixed-effect",
        "meta-analysis over per-study Welch tests on log2(normalized expression + 1)",
        "from the integrated HDF5/aligned target matrix.",
        "",
        f"- Accessions tested: {', '.join(summary['meta_dgea']['accessions_tested'])}",
        f"- Genes tested: {summary['meta_dgea']['genes_tested']:,}",
        (
            "- Significant genes at FDR < 0.05 and abs(log2FC) >= 1: "
            f"{summary['meta_dgea']['significant_fdr05_abs_log2fc1']}"
        ),
        f"- Up: {summary['meta_dgea']['significant_up']}",
        f"- Down: {summary['meta_dgea']['significant_down']}",
        "",
        "Top meta-analysis genes are in `top_meta_dgea_genes.tsv`.",
        "",
        "## Raw-Count DESeq2 Inputs",
        "",
        "The HDF5 expression matrix is integer count-like data. Per-study DESeq2",
        "inputs were exported for a study-aware FLT-vs-GC analysis:",
        "",
        f"- Counts: `{summary['deseq2_inputs']['counts']}`",
        f"- Metadata: `{summary['deseq2_inputs']['metadata']}`",
        f"- Gene symbols: `{summary['deseq2_inputs']['gene_symbols']}`",
    ]
    (output_dir / "POST_ANALYSIS_SUMMARY.md").write_text(
        "\n".join(text) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze aggregated OSDR liver GLARE FLT/GC results."
    )
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR)
    parser.add_argument("--osdr-h5", default=DEFAULT_OSDR_H5)
    parser.add_argument(
        "--output-dir",
        help="Defaults to <run-dir>/post_analysis.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    output_dir = Path(args.output_dir) if args.output_dir else run_dir / "post_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    log(f"Loading aggregate GLARE run from {run_dir}")
    run = load_run(run_dir)
    symbols = load_gene_symbols(args.osdr_h5)

    log("Running study/batch QC")
    batch_qc = run_batch_qc(run, output_dir)
    log("Comparing FLT and GC cluster structures")
    cluster_comparison = run_cluster_comparison(run, output_dir)
    log("Running exploratory per-study/meta DGEA")
    meta_dgea = run_meta_dgea(run, symbols, output_dir)
    log("Exporting raw-count DESeq2 inputs")
    deseq2_inputs = export_deseq2_inputs(run, symbols, output_dir)

    summary = {
        "run_dir": str(run_dir),
        "output_dir": str(output_dir),
        "batch_qc": batch_qc,
        "cluster_comparison": cluster_comparison,
        "meta_dgea": meta_dgea,
        "deseq2_inputs": deseq2_inputs,
        "notes": [
            "DGEA is exploratory normalized-expression meta-analysis, not raw-count DESeq2.",
            "Raw-count DESeq2 inputs are exported for per-study DESeq2 meta-analysis.",
            "Batch variables are partially confounded in this aggregate design.",
        ],
    }
    (output_dir / "post_analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_report(output_dir, summary)
    log(f"Saved post-analysis outputs to {output_dir}")


if __name__ == "__main__":
    main()
