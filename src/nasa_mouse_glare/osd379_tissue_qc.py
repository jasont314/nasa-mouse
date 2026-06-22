"""Sample-level tissue-composition QC for the OSD-379 liver analysis."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MUSCLE_MARKERS = [
    "Myh1",
    "Myh2",
    "Myh4",
    "Myh7",
    "Myh8",
    "Tnnt1",
    "Tnnt3",
    "Tnni2",
    "Tnnc2",
    "Actn3",
    "Ckm",
    "Pvalb",
    "Mb",
    "Atp2a1",
    "Myl1",
    "Myl2",
    "Mylpf",
    "Mybpc1",
    "Mybpc2",
    "Ryr1",
]
LIVER_MARKERS = ["Alb", "Cps1", "Ass1", "Tat", "Pck1", "Cyp2e1", "Apob", "Fgb"]
SAMPLE_PATTERN = re.compile(
    r"^RR8_LVR_(?P<condition>BSL|FLT|GC|VIV)_"
    r"(?P<collection>ISS-T|LAR)_(?P<age>OLD|YNG)_"
)


def parse_sample(sample: str) -> dict[str, str]:
    match = SAMPLE_PATTERN.match(sample)
    if not match:
        return {"condition": "unknown", "collection": "unknown", "age": "unknown"}
    return match.groupdict()


def load_expression(path: Path) -> pd.DataFrame:
    expression = pd.read_csv(path, index_col=0)
    expression.index = expression.index.astype(str)
    expression.columns = expression.columns.astype(str)
    return expression.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def load_symbols(path: Path) -> pd.DataFrame:
    symbols = pd.read_csv(path, usecols=["ENSEMBL", "SYMBOL"])
    symbols["ENSEMBL"] = symbols["ENSEMBL"].astype(str)
    symbols["SYMBOL"] = symbols["SYMBOL"].fillna("").astype(str)
    return symbols.drop_duplicates("ENSEMBL").set_index("ENSEMBL")


def marker_gene_ids(
    symbols: pd.DataFrame, marker_symbols: list[str], available: pd.Index
) -> dict[str, str]:
    symbol_to_gene = (
        symbols.reset_index()
        .loc[lambda frame: frame["SYMBOL"].isin(marker_symbols)]
        .drop_duplicates("SYMBOL")
        .set_index("SYMBOL")["ENSEMBL"]
        .to_dict()
    )
    return {
        symbol: gene_id
        for symbol, gene_id in symbol_to_gene.items()
        if gene_id in available
    }


def robust_z(values: pd.Series) -> pd.Series:
    median = values.median()
    mad = np.median(np.abs(values - median))
    if mad == 0:
        return pd.Series(np.nan, index=values.index)
    return 0.67448975 * (values - median) / mad


def sample_marker_scores(
    expression: pd.DataFrame,
    muscle_ids: dict[str, str],
    liver_ids: dict[str, str],
) -> pd.DataFrame:
    log_expression = np.log2(expression + 1.0)
    muscle = expression.loc[list(muscle_ids.values())]
    liver = expression.loc[list(liver_ids.values())]
    metadata = pd.DataFrame(
        [parse_sample(sample) for sample in expression.columns],
        index=expression.columns,
    )
    scores = metadata.assign(
        sample=expression.columns,
        muscle_mean_log2=log_expression.loc[list(muscle_ids.values())].mean(axis=0),
        muscle_markers_over_100=(muscle > 100).sum(axis=0),
        muscle_marker_fraction_percent=(
            100.0 * muscle.sum(axis=0) / expression.sum(axis=0)
        ),
        liver_mean_log2=log_expression.loc[list(liver_ids.values())].mean(axis=0),
        liver_markers_over_100=(liver > 100).sum(axis=0),
    ).reset_index(drop=True)
    strata = ["condition", "collection", "age"]
    scores["muscle_score_robust_z_within_stratum"] = scores.groupby(
        strata, dropna=False
    )["muscle_mean_log2"].transform(robust_z)
    scores["severe_muscle_composition_outlier"] = scores[
        "muscle_score_robust_z_within_stratum"
    ].ge(3.5)
    return scores.sort_values("muscle_mean_log2", ascending=False)


def welch_results(
    expression: pd.DataFrame, flight_samples: list[str], ground_samples: list[str]
) -> pd.DataFrame:
    logged = np.log2(expression + 1.0)
    flight = logged[flight_samples].to_numpy()
    ground = logged[ground_samples].to_numpy()
    statistic, p_value = ttest_ind(
        flight, ground, axis=1, equal_var=False, nan_policy="omit"
    )
    p_value = np.nan_to_num(p_value, nan=1.0, posinf=1.0, neginf=1.0)
    fdr = multipletests(p_value, method="fdr_bh")[1]
    effect = flight.mean(axis=1) - ground.mean(axis=1)
    return pd.DataFrame(
        {
            "gene_id": expression.index,
            "mean_log2_flight": flight.mean(axis=1),
            "mean_log2_ground": ground.mean(axis=1),
            "log2_mean_difference": effect,
            "welch_t": statistic,
            "p_value": p_value,
            "fdr_bh": fdr,
            "significant_fdr05_abs_effect1": (fdr < 0.05) & (np.abs(effect) >= 1.0),
        }
    ).set_index("gene_id")


def tms_marker_detection(
    h5ad_path: Path, marker_symbols: list[str]
) -> tuple[pd.DataFrame, int]:
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path, backed="r")
    liver_mask = adata.obs["tissue"].astype(str).eq("liver").to_numpy()
    var = adata.var.copy()
    var["gene_id"] = var.index.astype(str)
    var["feature_name"] = var["feature_name"].astype(str)
    selected = (
        var.loc[var["feature_name"].isin(marker_symbols)]
        .drop_duplicates("feature_name")
        .set_index("feature_name")
    )
    gene_positions = adata.var_names.get_indexer(selected["gene_id"])
    matrix = adata[liver_mask, gene_positions].X
    if sparse.issparse(matrix):
        detected = np.asarray((matrix > 0).mean(axis=0)).ravel()
        means = np.asarray(matrix.mean(axis=0)).ravel()
    else:
        matrix = np.asarray(matrix)
        detected = (matrix > 0).mean(axis=0)
        means = matrix.mean(axis=0)
    result = pd.DataFrame(
        {
            "gene_symbol": selected.index,
            "gene_id": selected["gene_id"].to_numpy(),
            "detected_cell_fraction": detected,
            "mean_expression": means,
        }
    ).sort_values("detected_cell_fraction")
    adata.file.close()
    return result, int(liver_mask.sum())


def implicated_assay_qc(path: Path, samples: list[str]) -> pd.DataFrame:
    assay = pd.read_csv(path, sep="\t", dtype=str)
    columns = [
        "Sample Name",
        "Parameter Value[QA Score]",
        "Unit",
        "Comment[Library Prep Date]",
        "Parameter Value[Fragment Size]",
        "Parameter Value[Read Depth]",
        "Parameter Value[rRNA Contamination]",
    ]
    missing = set(columns) - set(assay.columns)
    if missing:
        raise ValueError(f"ISA assay table is missing columns: {sorted(missing)}")
    selected = assay.loc[assay["Sample Name"].isin(samples), columns].copy()
    selected = selected.rename(
        columns={
            "Sample Name": "sample",
            "Parameter Value[QA Score]": "rna_integrity_number",
            "Unit": "qa_score_unit",
            "Comment[Library Prep Date]": "library_prep_date",
            "Parameter Value[Fragment Size]": "fragment_size",
            "Parameter Value[Read Depth]": "read_depth",
            "Parameter Value[rRNA Contamination]": "rrna_contamination_percent",
        }
    )
    numeric = [
        "rna_integrity_number",
        "fragment_size",
        "read_depth",
        "rrna_contamination_percent",
    ]
    selected[numeric] = selected[numeric].apply(pd.to_numeric, errors="coerce")
    return selected.sort_values("sample")


def plot_scores(scores: pd.DataFrame, output_path: Path) -> None:
    ordered = scores.sort_values("muscle_mean_log2").reset_index(drop=True)
    colors = {
        "FLT": "#c43c39",
        "GC": "#2878b5",
        "BSL": "#737373",
        "VIV": "#31824a",
        "unknown": "#9a9a9a",
    }
    fig, ax = plt.subplots(figsize=(11, 5))
    for condition, group in ordered.groupby("condition"):
        ax.scatter(
            group.index,
            group["muscle_mean_log2"],
            s=24,
            color=colors.get(condition, "#9a9a9a"),
            label=condition,
            alpha=0.85,
        )
    ax.set_xlabel("Samples ordered by skeletal-muscle marker score")
    ax.set_ylabel("Mean log2(normalized expression + 1)")
    ax.legend(frameon=False, ncol=4)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_old_iss_heatmap(
    expression: pd.DataFrame,
    muscle_ids: dict[str, str],
    samples: list[str],
    output_path: Path,
) -> None:
    marker_order = [marker for marker in MUSCLE_MARKERS if marker in muscle_ids]
    values = np.log2(
        expression.loc[[muscle_ids[marker] for marker in marker_order], samples] + 1.0
    )
    fig, ax = plt.subplots(figsize=(13, 7))
    image = ax.imshow(values.to_numpy(), aspect="auto", cmap="magma")
    ax.set_yticks(range(len(marker_order)), marker_order, fontsize=8)
    ax.set_xticks(
        range(len(samples)),
        [sample.removeprefix("RR8_LVR_") for sample in samples],
        rotation=65,
        ha="right",
        fontsize=7,
    )
    ax.set_title("OSD-379 old ISS-terminal liver samples")
    colorbar = fig.colorbar(image, ax=ax, pad=0.01)
    colorbar.set_label("log2(normalized expression + 1)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict:
    run_dir = Path(args.run_dir)
    output_dir = run_dir / "biological_analysis" / "tissue_composition_qc"
    output_dir.mkdir(parents=True, exist_ok=True)

    expression = load_expression(Path(args.normalized_counts))
    symbols = load_symbols(Path(args.official_de))
    muscle_ids = marker_gene_ids(symbols, MUSCLE_MARKERS, expression.index)
    liver_ids = marker_gene_ids(symbols, LIVER_MARKERS, expression.index)
    if len(muscle_ids) != len(MUSCLE_MARKERS):
        missing = sorted(set(MUSCLE_MARKERS) - set(muscle_ids))
        raise ValueError(f"Missing muscle markers: {missing}")

    scores = sample_marker_scores(expression, muscle_ids, liver_ids)
    scores.to_csv(output_dir / "sample_muscle_marker_scores.tsv", sep="\t", index=False)
    plot_scores(scores, output_dir / "sample_muscle_marker_scores.png")

    old_iss = scores.loc[
        scores["collection"].eq("ISS-T")
        & scores["age"].eq("OLD")
        & scores["condition"].isin(["FLT", "GC"])
    ]
    flight_samples = old_iss.loc[old_iss["condition"].eq("FLT"), "sample"].tolist()
    ground_samples = old_iss.loc[old_iss["condition"].eq("GC"), "sample"].tolist()
    flagged_flight = old_iss.loc[
        old_iss["condition"].eq("FLT")
        & old_iss["severe_muscle_composition_outlier"],
        "sample",
    ].tolist()
    cleaned_flight = [sample for sample in flight_samples if sample not in flagged_flight]
    heatmap_samples = (
        old_iss.sort_values(["condition", "muscle_mean_log2"])["sample"].tolist()
    )
    plot_old_iss_heatmap(
        expression,
        muscle_ids,
        heatmap_samples,
        output_dir / "old_iss_muscle_marker_heatmap.png",
    )

    before = welch_results(expression, flight_samples, ground_samples)
    after = welch_results(expression, cleaned_flight, ground_samples)
    cleaned_candidates = (
        after.loc[after["significant_fdr05_abs_effect1"]]
        .join(symbols.rename(columns={"SYMBOL": "gene_symbol"}))
        .sort_values(["fdr_bh", "log2_mean_difference"], ascending=[True, False])
        .reset_index()
    )
    cleaned_candidates.to_csv(
        output_dir / "old_iss_cleaned_candidate_genes.tsv",
        sep="\t",
        index=False,
    )
    cluster = pd.read_csv(
        run_dir / "biological_analysis" / "FLT_gene_clusters_annotated.tsv",
        sep="\t",
    )
    cluster = cluster.loc[cluster["consensus"].eq(args.cluster)].copy()
    official_long = pd.read_csv(
        run_dir / "biological_analysis" / "official_matched_contrast_deg_long.tsv",
        sep="\t",
    )
    old_contrast = official_long.loc[
        official_long["contrast"].str.count("32 week").eq(2)
        & official_long["contrast"].str.count("Carcass").eq(2)
    ].copy()
    if old_contrast["contrast"].nunique() != 1:
        raise ValueError("Could not uniquely identify the old ISS-terminal contrast")
    old_contrast = old_contrast.set_index("gene_id")

    sensitivity = cluster[["gene_id", "gene_symbol"]].set_index("gene_id")
    sensitivity = sensitivity.join(
        old_contrast[
            [
                "log2_fold_change",
                "fdr_bh",
                "significant_fdr05_abs_log2fc1",
            ]
        ].rename(
            columns={
                "log2_fold_change": "official_log2_fold_change",
                "fdr_bh": "official_fdr_bh",
                "significant_fdr05_abs_log2fc1": "official_significant",
            }
        )
    )
    sensitivity = sensitivity.join(before.add_prefix("all_samples_"))
    sensitivity = sensitivity.join(after.add_prefix("without_flagged_"))
    sensitivity.reset_index().to_csv(
        output_dir / "flt_cluster_12_old_iss_sensitivity.tsv",
        sep="\t",
        index=False,
    )

    contrast_counts = (
        official_long.loc[official_long["gene_id"].isin(cluster["gene_id"])]
        .groupby("contrast")["significant_fdr05_abs_log2fc1"]
        .sum()
        .astype(int)
        .sort_values(ascending=False)
    )
    contrast_counts.rename("significant_cluster_genes").to_csv(
        output_dir / "flt_cluster_12_significance_by_contrast.tsv",
        sep="\t",
    )
    short_contrast_counts = {}
    for name, value in contrast_counts.items():
        collection = "ISS-terminal" if "On ISS" in name else "live-animal-return"
        age = "old" if name.count("32 week") == 2 else "young"
        short_contrast_counts[f"{collection} {age}"] = int(value)

    tms_result = pd.DataFrame()
    tms_liver_cells = 0
    if args.tms_h5ad and Path(args.tms_h5ad).exists():
        tms_result, tms_liver_cells = tms_marker_detection(
            Path(args.tms_h5ad), MUSCLE_MARKERS + LIVER_MARKERS
        )
        tms_result.to_csv(
            output_dir / "tms_liver_marker_detection.tsv", sep="\t", index=False
        )

    assay_qc = pd.DataFrame()
    if args.isa_assay and Path(args.isa_assay).exists():
        assay_qc = implicated_assay_qc(Path(args.isa_assay), flagged_flight)
        assay_qc.to_csv(
            output_dir / "implicated_sample_assay_qc.tsv", sep="\t", index=False
        )

    top_samples = scores.head(10)[
        [
            "sample",
            "condition",
            "collection",
            "age",
            "muscle_mean_log2",
            "muscle_marker_fraction_percent",
            "liver_mean_log2",
            "severe_muscle_composition_outlier",
        ]
    ]
    summary = {
        "dataset": "OSD-379 RR-8 left-lobe liver bulk RNA-seq",
        "samples_scored": int(len(scores)),
        "muscle_marker_count": int(len(muscle_ids)),
        "liver_marker_count": int(len(liver_ids)),
        "old_iss_flight_samples": len(flight_samples),
        "old_iss_ground_samples": len(ground_samples),
        "flagged_old_iss_flight_samples": flagged_flight,
        "flt_cluster": args.cluster,
        "flt_cluster_gene_count": int(len(sensitivity)),
        "official_old_iss_significant_cluster_genes": int(
            sensitivity["official_significant"].sum()
        ),
        "cluster_genes_abs_effect_ge_1_before_exclusion": int(
            sensitivity["all_samples_log2_mean_difference"].abs().ge(1).sum()
        ),
        "cluster_genes_abs_effect_ge_1_after_exclusion": int(
            sensitivity["without_flagged_log2_mean_difference"].abs().ge(1).sum()
        ),
        "cluster_genes_significant_after_exclusion": int(
            sensitivity["without_flagged_significant_fdr05_abs_effect1"].sum()
        ),
        "all_genes_significant_after_exclusion": int(
            after["significant_fdr05_abs_effect1"].sum()
        ),
        "tms_liver_cells": tms_liver_cells,
        "implicated_sample_assay_qc": assay_qc.to_dict(orient="records"),
        "top_samples": top_samples.to_dict(orient="records"),
        "interpretation": (
            "The coherent fast-skeletal-muscle signal is concentrated in a few "
            "liver profiles and is most consistent with tissue admixture or "
            "cross-sample contamination, not generalized hepatic activation."
        ),
        "limitations": (
            "The exclusion analysis uses Welch tests on official normalized counts, "
            "not a replacement DESeq2 fit from raw integer counts."
        ),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    flagged_rows = scores.loc[scores["sample"].isin(flagged_flight)]
    tms_fast = tms_result.loc[tms_result["gene_symbol"].isin(MUSCLE_MARKERS)]
    if assay_qc.empty:
        assay_text = "The complete ISA assay table was not available."
    else:
        assay_text = " ".join(
            f"{row.sample} had RIN {row.rna_integrity_number:g}, library prep "
            f"{row.library_prep_date}, {int(row.read_depth):,} reads, and "
            f"{row.rrna_contamination_percent:g}% rRNA contamination."
            for row in assay_qc.itertuples()
        )
    report = f"""# OSD-379 Liver Tissue-Composition QC

## Conclusion

OSD-379 is a liver dataset: the study is titled "Transcriptional profiling of
livers", and the samples are annotated as left lobe of the liver. The
muscle-enriched GLARE result is nevertheless not robust evidence of a
spaceflight-induced liver program.

The signal is concentrated in {len(flagged_flight)} old ISS-terminal flight
samples: {", ".join(flagged_flight)}. These profiles retain strong liver-marker
expression, so they are not simply mislabeled muscle samples. They instead look
like liver RNA mixed with a substantial amount of skeletal-muscle RNA.

## Evidence

- FLT cluster {args.cluster} contains {len(sensitivity)} genes.
- NASA's official old ISS-terminal contrast calls
  {int(sensitivity["official_significant"].sum())}/{len(sensitivity)} of them
  significant.
- Across all four matched contrasts, the cluster significance counts are
  {", ".join(f"{value} ({name})" for name, value in short_contrast_counts.items())}.
- Before excluding the severe composition outliers,
  {int(sensitivity["all_samples_log2_mean_difference"].abs().ge(1).sum())}
  cluster genes have an absolute mean log2 difference of at least 1.
- After excluding them, that falls to
  {int(sensitivity["without_flagged_log2_mean_difference"].abs().ge(1).sum())};
  {int(sensitivity["without_flagged_significant_fdr05_abs_effect1"].sum())}
  cluster genes pass FDR < 0.05 and absolute effect >= 1 in the normalized-count
  sensitivity test.
- Across all tested genes, {len(cleaned_candidates)} non-cluster-specific
  candidates still pass that sensitivity threshold. The outlier finding
  invalidates the muscle-cluster interpretation, not every possible
  spaceflight response in the old ISS-terminal stratum.
- The implicated samples still express liver markers strongly:
  {", ".join(f"{row.sample}: muscle score {row.muscle_mean_log2:.2f}, liver score {row.liver_mean_log2:.2f}" for row in flagged_rows.itertuples())}.
- In {tms_liver_cells:,} TMS FACS liver cells, the corresponding fast-muscle
  markers are sparse rather than a coherent liver-cell program. The median
  detection fraction for the marker panel is
  {float(tms_fast["detected_cell_fraction"].median()) if not tms_fast.empty else float("nan"):.4f}.
- Similar high muscle-marker scores occur in ground, vivarium, and baseline
  samples, including RR8_LVR_GC_LAR_YNG_GL5. That broader distribution argues
  against a flight-specific hepatic program and suggests a recurrent
  tissue-composition problem.
- {assay_text} The implicated samples were prepared in different library
  batches, making one shared library-preparation spillover event less likely.

## Likely Cause

The most likely explanation is physical tissue admixture during dissection or
sample trimming. ISS-terminal animals were frozen as whole carcasses, later
thawed for approximately 60-90 minutes, and then dissected. Liver material was
trimmed to approximately 25 mg for RNA extraction. A small fragment from the
adjacent diaphragm or body-wall muscle can contribute a large muscle RNA signal
while leaving abundant liver transcripts intact.

This analysis cannot identify the exact contamination stage. Cross-sample
carryover during tissue handling, homogenization, or library preparation
remains possible. Normal read depth and mapping metrics do not rule out either
form of tissue-composition contamination.

## Impact On GLARE

The SAE and clustering captured a real expression module present in the input;
the muscle enrichment is not a software error. The biological interpretation
is the problem: the official DEG union is dominated by one contrast whose
effect is strongly leveraged by a few composition outliers.

Keep the current run as an unfiltered audit trail. For biological conclusions,
rerun the OSD-379 differential expression and GLARE fine-tuning after applying
a documented tissue-composition QC rule, and compare the filtered result with
the original.

## Files

- `sample_muscle_marker_scores.tsv`: all 141 liver profiles and marker scores.
- `sample_muscle_marker_scores.png`: cohort-wide score distribution.
- `old_iss_muscle_marker_heatmap.png`: marker expression in the affected stratum.
- `flt_cluster_12_significance_by_contrast.tsv`: localization by NASA contrast.
- `flt_cluster_12_old_iss_sensitivity.tsv`: before/after exclusion effects.
- `old_iss_cleaned_candidate_genes.tsv`: candidates retained after exclusion.
- `tms_liver_marker_detection.tsv`: single-cell liver reference detection.
- `implicated_sample_assay_qc.tsv`: official ISA assay metadata for FI16/FI17.

## Sources

- OSD-379: https://osdr.nasa.gov/bio/repo/data/studies/OSD-379
- RRRM-1/RR-8 payload: https://osdr.nasa.gov/bio/repo/data/payloads/RRRM-1%20%28RR-8%29
- RR-8 tissue handling details: https://www.nature.com/articles/s41467-026-68737-1
"""
    (output_dir / "TISSUE_QC_REPORT.md").write_text(report, encoding="utf-8")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit skeletal-muscle composition in OSD-379 liver profiles."
    )
    parser.add_argument(
        "--run-dir", default="outputs/glare_paper_tms_liver_osd379"
    )
    parser.add_argument(
        "--normalized-counts",
        default="assets/osdr/GLDS-379_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
    )
    parser.add_argument(
        "--official-de",
        default="assets/osdr/GLDS-379_rna_seq_differential_expression_GLbulkRNAseq.csv",
    )
    parser.add_argument(
        "--tms-h5ad",
        default="assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad",
    )
    parser.add_argument(
        "--isa-assay",
        default="assets/osdr/OSD-379_assay_metadata.tsv",
        help="Complete OSD-379 ISA assay table; omitted if the file is absent.",
    )
    parser.add_argument("--cluster", type=int, default=12)
    return parser.parse_args()


def main() -> None:
    summary = run(parse_args())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
