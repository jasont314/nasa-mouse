"""Expression summaries and stratified significance tests for gene clusters."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from .cluster_enrichment import bh_fdr
from .io import load_matrix_bundle, require_import


DEFAULT_GENE_CLUSTERS = (
    "outputs/glare_hpt_tms_facs_osdr/post_finetune/"
    "ensemble_clustering/gene_clusters.tsv"
)
DEFAULT_TARGET_MANIFEST = "data/processed/tms_facs_osdr_aligned.target.manifest.json"
DEFAULT_PROFILE_METADATA = (
    "outputs/glare_hpt_tms_facs_osdr/post_finetune/profile_metadata.tsv"
)
DEFAULT_OUTPUT_DIR = (
    "outputs/glare_hpt_tms_facs_osdr/post_finetune/ensemble_analysis"
)


TISSUE_PATTERNS = [
    ("liver", r"(^|[_-])LVR([_-]|$)"),
    ("kidney", r"(^|[_-])(R-)?KDN([_-]|$)"),
    ("spleen", r"(^|[_-])SPL([_-]|$)"),
    ("lung", r"(^|[_-])LNG([_-]|$)"),
    ("thymus", r"(^|[_-])TMS([_-]|$)"),
    ("heart", r"(^|[_-])HRT([_-]|$)"),
    ("hippocampus", r"(^|[_-])HPC([_-]|$)"),
    ("cerebellum", r"(^|[_-])CB([_-]|$)"),
    ("brain", r"(^|[_-])BRN([_-]|$)"),
    ("retina", r"(^|[_-])(LRTN|RRTN|RTN)([_-]|$)"),
    ("optic_nerve", r"(^|[_-])ON([_-]|$)"),
    ("eye", r"(^|[_-])EYE([_-]|$)"),
    ("adrenal_gland", r"(^|[_-])ADR([_-]|$)"),
    ("colon", r"(^|[_-])CLN([_-]|$)"),
    ("cecum", r"(^|[_-])CEC([_-]|$)"),
    ("mammary_gland", r"(^|[_-])MG([_-]|$)"),
    ("skin", r"(^|[_-])(DSKN|FSKN|SKN)([_-]|$)"),
    (
        "skeletal_muscle",
        r"(^|[_-])(R-)?(EDL|QUAD|SLS|TA|GST|GST-R)([_-]|$)",
    ),
]


def infer_tissue(profile: str) -> str:
    value = str(profile).upper()
    for tissue, pattern in TISSUE_PATTERNS:
        if re.search(pattern, value):
            return tissue
    return "unknown"


def tissue_inference_rule(profile: str) -> str:
    value = str(profile).upper()
    for tissue, pattern in TISSUE_PATTERNS:
        match = re.search(pattern, value)
        if match:
            token = match.group(0).strip("_-")
            return f"{token}->{tissue}"
    return ""


def validate_alignment(bundle, gene_clusters, metadata) -> None:
    genes = gene_clusters["gene"].astype(str).tolist()
    profiles = metadata["profile"].astype(str).tolist()
    if list(map(str, bundle.genes)) != genes:
        raise SystemExit(
            "Target manifest genes do not match the supplied gene cluster order."
        )
    if list(map(str, bundle.profiles)) != profiles:
        raise SystemExit(
            "Target manifest profiles do not match profile metadata order."
        )


def compute_cluster_sample_expression(bundle, gene_clusters):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    labels = gene_clusters["gene_cluster"].astype(int).to_numpy()
    clusters = sorted(set(labels.tolist()))
    values = []
    for cluster in clusters:
        gene_mask = labels == cluster
        sample_means = np.asarray(bundle.matrix[gene_mask, :].mean(axis=0)).ravel()
        values.append(sample_means)
    matrix = np.vstack(values)
    return clusters, matrix, pd.Series(
        [(labels == cluster).sum() for cluster in clusters],
        index=clusters,
        dtype=int,
    )


def write_sample_expression(clusters, expression, n_genes, metadata, output_dir: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    metadata_cols = [
        "profile",
        "id.accession",
        "investigation.study.comment.project identifier",
        "study.source name",
        "condition_inferred",
        "flight_status_inferred",
        "tissue_inferred",
        "tissue_inference_rule",
        "tissue_inference_confidence",
    ]
    available_cols = [col for col in metadata_cols if col in metadata]
    for row_idx, cluster in enumerate(clusters):
        frame = metadata[available_cols].copy()
        frame.insert(0, "n_genes", int(n_genes[cluster]))
        frame.insert(0, "gene_cluster", cluster)
        frame["mean_expression"] = expression[row_idx]
        rows.append(frame)
    result = pd.concat(rows, ignore_index=True)
    path = output_dir / "cluster_sample_expression.tsv"
    result.to_csv(path, sep="\t", index=False)
    return result, path


def write_group_summary(
    sample_expression,
    group_col: str,
    output_dir: Path,
    filename: str,
):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    grouped = (
        sample_expression.groupby(["gene_cluster", "n_genes", group_col], dropna=False)[
            "mean_expression"
        ]
        .mean()
        .reset_index()
    )
    wide = grouped.pivot(
        index=["gene_cluster", "n_genes"],
        columns=group_col,
        values="mean_expression",
    ).reset_index()
    wide.columns.name = None
    path = output_dir / filename
    wide.to_csv(path, sep="\t", index=False)
    return path


def write_flight_ground_summary(sample_expression, output_dir: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    focus = sample_expression[
        sample_expression["condition_inferred"].isin(["flight", "ground_control"])
    ]
    grouped = (
        focus.groupby(["gene_cluster", "n_genes", "condition_inferred"])[
            "mean_expression"
        ]
        .mean()
        .unstack()
        .reset_index()
    )
    grouped["flight_mean"] = grouped["flight"]
    grouped["ground_or_control_mean"] = grouped["ground_control"]
    grouped["flight_minus_ground_or_control"] = (
        grouped["flight"] - grouped["ground_control"]
    )
    grouped["comparison"] = "flight_vs_ground_control"
    result = grouped[
        [
            "gene_cluster",
            "n_genes",
            "flight_mean",
            "ground_or_control_mean",
            "flight_minus_ground_or_control",
            "comparison",
        ]
    ]
    path = output_dir / "gene_cluster_flight_ground_summary.tsv"
    result.to_csv(path, sep="\t", index=False)
    return path


def paired_effects(sample_expression, strata_cols: list[str]):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    focus = sample_expression[
        sample_expression["condition_inferred"].isin(["flight", "ground_control"])
    ].copy()
    grouping = ["gene_cluster", "n_genes", *strata_cols, "condition_inferred"]
    grouped = (
        focus.groupby(grouping, dropna=False)["mean_expression"].agg(["mean", "count"])
        .reset_index()
    )
    means = grouped.pivot(
        index=["gene_cluster", "n_genes", *strata_cols],
        columns="condition_inferred",
        values="mean",
    )
    counts = grouped.pivot(
        index=["gene_cluster", "n_genes", *strata_cols],
        columns="condition_inferred",
        values="count",
    )
    valid = means.dropna(subset=["flight", "ground_control"]).reset_index()
    valid["n_flight_samples"] = (
        counts.loc[valid.set_index(["gene_cluster", "n_genes", *strata_cols]).index, "flight"]
        .to_numpy()
        .astype(int)
    )
    valid["n_ground_samples"] = (
        counts.loc[
            valid.set_index(["gene_cluster", "n_genes", *strata_cols]).index,
            "ground_control",
        ]
        .to_numpy()
        .astype(int)
    )
    valid["flight_minus_ground"] = valid["flight"] - valid["ground_control"]
    valid["log2_flight_ground_ratio"] = np.log2(
        (valid["flight"] + 1.0) / (valid["ground_control"] + 1.0)
    )
    return valid


def signed_rank_summary(effects, group_cols: list[str], min_pairs: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import(
        "scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt"
    )

    rows = []
    for keys, group in effects.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        values = group["log2_flight_ground_ratio"].astype(float).to_numpy()
        nonzero = values[values != 0]
        n_pairs = len(values)
        if n_pairs >= min_pairs and len(nonzero):
            wilcoxon_p = float(
                scipy_stats.wilcoxon(
                    nonzero,
                    zero_method="wilcox",
                    alternative="two-sided",
                    method="auto",
                ).pvalue
            )
            n_positive = int((values > 0).sum())
            n_negative = int((values < 0).sum())
            sign_n = n_positive + n_negative
            sign_p = (
                float(
                    scipy_stats.binomtest(
                        n_positive,
                        sign_n,
                        p=0.5,
                        alternative="two-sided",
                    ).pvalue
                )
                if sign_n
                else 1.0
            )
        else:
            wilcoxon_p = np.nan
            sign_p = np.nan
            n_positive = int((values > 0).sum())
            n_negative = int((values < 0).sum())
        row = dict(zip(group_cols, keys))
        row.update(
            {
                "n_paired_strata": n_pairs,
                "mean_log2_flight_ground_ratio": float(np.mean(values)),
                "median_log2_flight_ground_ratio": float(np.median(values)),
                "mean_flight_minus_ground": float(
                    group["flight_minus_ground"].mean()
                ),
                "positive_strata": n_positive,
                "negative_strata": n_negative,
                "wilcoxon_p_value": wilcoxon_p,
                "sign_test_p_value": sign_p,
            }
        )
        rows.append(row)
    result = pd.DataFrame(rows)
    result["wilcoxon_fdr_bh"] = bh_fdr_with_nan(result["wilcoxon_p_value"])
    result["sign_test_fdr_bh"] = bh_fdr_with_nan(result["sign_test_p_value"])
    return result


def bh_fdr_with_nan(values):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    array = np.asarray(values, dtype=float)
    result = np.full(len(array), np.nan)
    valid = np.isfinite(array)
    if valid.any():
        result[valid] = bh_fdr(array[valid])
    return result


def eta_squared(values, groups) -> float:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    frame = pd.DataFrame({"value": values, "group": groups})
    frame = frame[frame["group"].astype(str).ne("")]
    if frame.empty or frame["group"].nunique() < 2:
        return float("nan")
    grand_mean = frame["value"].mean()
    total_ss = float(((frame["value"] - grand_mean) ** 2).sum())
    if total_ss == 0:
        return 0.0
    between_ss = 0.0
    for _, group in frame.groupby("group"):
        between_ss += len(group) * float((group["value"].mean() - grand_mean) ** 2)
    return between_ss / total_ss


def write_variance_summary(sample_expression, output_dir: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    group_cols = [
        "id.accession",
        "investigation.study.comment.project identifier",
        "tissue_inferred",
        "condition_inferred",
    ]
    for cluster, group in sample_expression.groupby("gene_cluster"):
        row = {
            "gene_cluster": int(cluster),
            "n_genes": int(group["n_genes"].iloc[0]),
        }
        for column in group_cols:
            row[f"eta_squared_{column}"] = eta_squared(
                group["mean_expression"], group[column]
            )
        rows.append(row)
    result = pd.DataFrame(rows)
    path = output_dir / "stratification_variance_summary.tsv"
    result.to_csv(path, sep="\t", index=False)
    return path


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bundle = load_matrix_bundle(args.target_manifest)
    gene_clusters = pd.read_csv(args.gene_clusters, sep="\t")
    gene_clusters["gene_cluster"] = gene_clusters["gene_cluster"].astype(int)
    metadata = pd.read_csv(args.profile_metadata, sep="\t", keep_default_na=False)
    metadata["tissue_inferred"] = metadata["profile"].map(infer_tissue)
    metadata["tissue_inference_rule"] = metadata["profile"].map(
        tissue_inference_rule
    )
    metadata["tissue_inference_confidence"] = metadata["tissue_inferred"].map(
        lambda value: "high_explicit_sample_token"
        if value != "unknown"
        else "unassigned"
    )
    validate_alignment(bundle, gene_clusters, metadata)

    gene_clusters_path = output_dir / "gene_clusters.tsv"
    gene_clusters.to_csv(gene_clusters_path, sep="\t", index=False)
    metadata_path = output_dir / "profile_metadata.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)

    clusters, expression, n_genes = compute_cluster_sample_expression(
        bundle, gene_clusters
    )
    sample_expression, sample_expression_path = write_sample_expression(
        clusters,
        expression,
        n_genes,
        metadata,
        output_dir,
    )

    condition_path = write_group_summary(
        sample_expression,
        "condition_inferred",
        output_dir,
        "gene_cluster_expression_by_condition_inferred.tsv",
    )
    accession_path = write_group_summary(
        sample_expression,
        "id.accession",
        output_dir,
        "gene_cluster_expression_by_id.accession.tsv",
    )
    project_path = write_group_summary(
        sample_expression,
        "investigation.study.comment.project identifier",
        output_dir,
        "gene_cluster_expression_by_project_identifier.tsv",
    )
    tissue_path = write_group_summary(
        sample_expression,
        "tissue_inferred",
        output_dir,
        "gene_cluster_expression_by_tissue_inferred.tsv",
    )
    shift_path = write_flight_ground_summary(sample_expression, output_dir)

    accession_effects = paired_effects(sample_expression, ["id.accession"])
    accession_effects_path = output_dir / "flight_ground_effects_by_accession.tsv"
    accession_effects.to_csv(accession_effects_path, sep="\t", index=False)
    accession_tests = signed_rank_summary(
        accession_effects,
        ["gene_cluster"],
        args.min_paired_studies,
    )
    accession_tests_path = output_dir / "flight_ground_significance_by_accession.tsv"
    accession_tests.to_csv(accession_tests_path, sep="\t", index=False)

    tissue_effects = paired_effects(
        sample_expression,
        ["id.accession", "tissue_inferred"],
    )
    tissue_effects_path = output_dir / "flight_ground_effects_by_accession_tissue.tsv"
    tissue_effects.to_csv(tissue_effects_path, sep="\t", index=False)
    known_tissue_effects = tissue_effects[
        tissue_effects["tissue_inferred"].ne("unknown")
    ]
    tissue_tests = signed_rank_summary(
        known_tissue_effects,
        ["gene_cluster", "tissue_inferred"],
        args.min_paired_studies,
    )
    tissue_tests_path = output_dir / "flight_ground_significance_by_tissue.tsv"
    tissue_tests.to_csv(tissue_tests_path, sep="\t", index=False)

    variance_path = write_variance_summary(sample_expression, output_dir)

    tissue_counts = (
        metadata.groupby(["tissue_inferred", "condition_inferred"])
        .size()
        .rename("n_samples")
        .reset_index()
    )
    tissue_counts_path = output_dir / "tissue_inference_counts.tsv"
    tissue_counts.to_csv(tissue_counts_path, sep="\t", index=False)

    known_tissue = metadata[metadata["tissue_inferred"].ne("unknown")]
    accession_tissue_counts = known_tissue.groupby("id.accession")[
        "tissue_inferred"
    ].nunique()
    multi_tissue_accessions = accession_tissue_counts[
        accession_tissue_counts > 1
    ].index.tolist()
    flight_ground = metadata[
        metadata["condition_inferred"].isin(["flight", "ground_control"])
    ]
    tissue_audit = {
        "accuracy_status": (
            "not directly measurable because the source HDF5 has no explicit "
            "tissue ground-truth column"
        ),
        "method": "deterministic explicit tissue tokens in profile sample IDs",
        "all_profile_coverage": float(
            metadata["tissue_inferred"].ne("unknown").mean()
        ),
        "flight_ground_profile_coverage": float(
            flight_ground["tissue_inferred"].ne("unknown").mean()
        ),
        "accessions_with_inferred_tissue": int(len(accession_tissue_counts)),
        "single_tissue_accessions": int(accession_tissue_counts.eq(1).sum()),
        "multi_tissue_accessions": multi_tissue_accessions,
        "internal_consistency_note": (
            "OSD-164 is the only multi-tissue accession and explicitly contains "
            "liver and spleen sample IDs."
        ),
    }
    tissue_audit_path = output_dir / "tissue_inference_audit.json"
    tissue_audit_path.write_text(
        json.dumps(tissue_audit, indent=2) + "\n",
        encoding="utf-8",
    )

    summary = {
        "gene_clusters": str(args.gene_clusters),
        "target_manifest": str(args.target_manifest),
        "profile_metadata": str(args.profile_metadata),
        "n_clusters": len(clusters),
        "n_profiles": len(metadata),
        "n_accessions": int(metadata["id.accession"].nunique()),
        "tissue_inference": {
            "source": "profile sample-name tokens; not an explicit OSDR tissue field",
            "known_profiles": int(metadata["tissue_inferred"].ne("unknown").sum()),
            "unknown_profiles": int(metadata["tissue_inferred"].eq("unknown").sum()),
            "known_flight_ground_profiles": int(
                metadata[
                    metadata["condition_inferred"].isin(
                        ["flight", "ground_control"]
                    )
                ]["tissue_inferred"]
                .ne("unknown")
                .sum()
            ),
            "audit": tissue_audit,
        },
        "significance": {
            "effect": "log2((mean flight expression + 1)/(mean ground expression + 1))",
            "test": "two-sided Wilcoxon signed-rank across paired accessions",
            "secondary_test": "two-sided binomial sign test",
            "multiple_testing": "Benjamini-Hochberg FDR",
            "min_paired_studies": args.min_paired_studies,
            "accession_cluster_tests": int(len(accession_tests)),
            "accession_nominal_p_lt_0_05": int(
                accession_tests["wilcoxon_p_value"].lt(0.05).sum()
            ),
            "accession_fdr_lt_0_05": int(
                accession_tests["wilcoxon_fdr_bh"].lt(0.05).sum()
            ),
            "tissue_cluster_tests": int(len(tissue_tests)),
            "tissue_nominal_p_lt_0_05": int(
                tissue_tests["wilcoxon_p_value"].lt(0.05).sum()
            ),
            "tissue_fdr_lt_0_05": int(
                tissue_tests["wilcoxon_fdr_bh"].lt(0.05).sum()
            ),
        },
        "outputs": {
            "gene_clusters": str(gene_clusters_path),
            "profile_metadata": str(metadata_path),
            "cluster_sample_expression": str(sample_expression_path),
            "condition_summary": str(condition_path),
            "accession_summary": str(accession_path),
            "project_summary": str(project_path),
            "tissue_summary": str(tissue_path),
            "flight_ground_summary": str(shift_path),
            "accession_effects": str(accession_effects_path),
            "accession_significance": str(accession_tests_path),
            "accession_tissue_effects": str(tissue_effects_path),
            "tissue_significance": str(tissue_tests_path),
            "variance_summary": str(variance_path),
            "tissue_inference_counts": str(tissue_counts_path),
            "tissue_inference_audit": str(tissue_audit_path),
        },
    }
    summary_path = output_dir / "cluster_stratified_analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"clusters={len(clusters)}")
    print(
        "known_tissue_profiles="
        f"{summary['tissue_inference']['known_profiles']}/{len(metadata)}"
    )
    print(f"summary={summary_path}")
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute study/tissue-stratified cluster expression tests."
    )
    parser.add_argument("--gene-clusters", default=DEFAULT_GENE_CLUSTERS)
    parser.add_argument("--target-manifest", default=DEFAULT_TARGET_MANIFEST)
    parser.add_argument("--profile-metadata", default=DEFAULT_PROFILE_METADATA)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-paired-studies", type=int, default=3)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
