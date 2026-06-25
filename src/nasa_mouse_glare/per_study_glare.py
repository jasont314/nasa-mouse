"""Run GLARE independently for selected liver FLT-vs-GC OSDR studies."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import torch

from .aggregate_liver_finetune import (
    DEFAULT_OSDR_H5,
    DEFAULT_TARGET_MANIFEST,
    load_excluded_profiles,
    select_aggregate_profiles,
)
from .cluster_enrichment import bh_fdr, read_gmt
from .io import require_import
from .paper_finetune import (
    finetune_location,
    format_elapsed,
    infer_pretrain_input_dim,
    log,
    write_outlier_audit,
)


DEFAULT_STUDIES = ["OSD-379", "OSD-245", "OSD-463", "OSD-168", "OSD-48", "OSD-137"]
DEFAULT_OUTPUT_DIR = "outputs/glare_per_study_liver_noercc_12filter"
DEFAULT_PRETRAINED_WEIGHTS = (
    "outputs/glare_paper_tms_liver_osd379/pretraining/sc_shulse_pretrained_reproduced.pth"
)
DEFAULT_EXCLUDE_PROFILES = "data/filters/aggregate_liver_12_muscle_candidate_profiles.txt"
DEFAULT_REACTOME_GMT = (
    "src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt"
)


def run_clustering(run_dir: Path, skip_tsne: bool) -> None:
    command = [
        sys.executable,
        "-m",
        "nasa_mouse_glare.paper_clustering",
        "--run-dir",
        str(run_dir),
    ]
    if skip_tsne:
        command.append("--skip-tsne")
    subprocess.run(command, check=True)


def cluster_enrichment_for_run(
    run_dir: Path,
    gmt_path: Path,
    min_overlap: int,
) -> pd.DataFrame:
    hypergeom = require_import(
        "scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt"
    ).hypergeom

    target = np.load(run_dir / "controlled_target.npz")
    genes = target["genes"].astype(str).tolist()
    universe = set(genes)
    gene_sets = read_gmt(gmt_path)
    rows = []
    for location in ("FLT", "GC"):
        clusters_path = run_dir / "clustering" / f"{location}_gene_clusters.tsv"
        if not clusters_path.exists():
            continue
        clusters = pd.read_csv(clusters_path, sep="\t")
        for cluster, cluster_table in clusters.groupby("consensus"):
            query = set(cluster_table["gene_id"].astype(str)) & universe
            if len(query) < min_overlap:
                continue
            for gene_set in gene_sets:
                term_genes = gene_set["genes"] & universe
                overlap = query & term_genes
                if len(overlap) < min_overlap:
                    continue
                pvalue = float(
                    hypergeom.sf(
                        len(overlap) - 1,
                        len(universe),
                        len(term_genes),
                        len(query),
                    )
                )
                rows.append(
                    {
                        "accession": run_dir.name,
                        "location": location,
                        "cluster": int(cluster),
                        "term": gene_set["term"],
                        "description": gene_set["description"],
                        "cluster_genes": len(query),
                        "overlap": len(overlap),
                        "term_genes_in_universe": len(term_genes),
                        "universe_genes": len(universe),
                        "p_value": pvalue,
                        "overlap_genes": ",".join(sorted(overlap)),
                    }
                )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr_bh"] = result.groupby(
        ["accession", "location", "cluster"], group_keys=False
    )["p_value"].transform(bh_fdr)
    return result.sort_values(
        ["accession", "location", "cluster", "fdr_bh", "p_value"],
        ascending=[True, True, True, True, True],
    )


def summarize_recurring_terms(enrichment: pd.DataFrame, alpha: float) -> pd.DataFrame:
    if enrichment.empty:
        return pd.DataFrame()
    significant = enrichment[enrichment["fdr_bh"] < alpha].copy()
    if significant.empty:
        return pd.DataFrame()
    rows = []
    for (location, term), group in significant.groupby(["location", "term"]):
        rows.append(
            {
                "location": location,
                "term": term,
                "study_count": group["accession"].nunique(),
                "accessions": ",".join(sorted(group["accession"].unique())),
                "cluster_count": len(group[["accession", "cluster"]].drop_duplicates()),
                "best_fdr_bh": group["fdr_bh"].min(),
                "max_overlap": int(group["overlap"].max()),
                "example_clusters": ";".join(
                    f"{row.accession}:{row.location}{int(row.cluster)}"
                    for row in group.sort_values("fdr_bh").head(8).itertuples()
                ),
            }
        )
    recurrence = pd.DataFrame(rows)
    return recurrence.sort_values(
        ["study_count", "best_fdr_bh", "term"],
        ascending=[False, True, True],
    )


def markdown_table(frame: pd.DataFrame, columns: list[str], max_rows: int = 20) -> list[str]:
    if frame.empty:
        return ["No rows."]
    display = frame.loc[:, columns].head(max_rows)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in display.iterrows():
        values = []
        for column in columns:
            value = row[column]
            if isinstance(value, float):
                value = f"{value:.3g}"
            values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def write_report(
    output_dir: Path,
    study_summaries: list[dict],
    recurrence: pd.DataFrame,
    alpha: float,
) -> None:
    study_rows = []
    for summary in study_summaries:
        condition_counts = pd.DataFrame(summary["condition_counts"])
        if {"FLT", "GC"}.issubset(condition_counts.columns):
            flt = int(condition_counts["FLT"].iloc[0])
            gc = int(condition_counts["GC"].iloc[0])
        elif "condition_label" in condition_counts.columns:
            pivot = condition_counts.pivot_table(
                index="h5_accession",
                columns="condition_label",
                values="total",
                aggfunc="sum",
                fill_value=0,
            )
            flt = int(pivot["FLT"].iloc[0]) if "FLT" in pivot.columns and len(pivot) else 0
            gc = int(pivot["GC"].iloc[0]) if "GC" in pivot.columns and len(pivot) else 0
        else:
            flt = 0
            gc = 0
        locations = {row["location"]: row for row in summary["locations"]}
        study_rows.append(
            {
                "accession": summary["accession"],
                "FLT": flt,
                "GC": gc,
                "ercc_dropped": summary["selection"]["ercc_policy"]["profiles_dropped"],
                "unique_wERCC_retained": summary["selection"]["ercc_policy"][
                    "unique_wERCC_profiles_retained"
                ],
                "FLT_best_loss": locations.get("FLT", {}).get("best_loss", ""),
                "GC_best_loss": locations.get("GC", {}).get("best_loss", ""),
            }
        )
    study_table = pd.DataFrame(study_rows)
    if recurrence.empty or "study_count" not in recurrence.columns:
        recurring = pd.DataFrame()
    else:
        recurring = recurrence[recurrence["study_count"] >= 2].copy()
    lines = [
        "# Per-Study Liver GLARE",
        "",
        "GLARE was run independently for each liver study. FLT and GC were fine-tuned separately per study.",
        "ERCC handling uses `prefer_noercc`: when both wERCC and noERCC profiles exist for the same biological sample, noERCC is retained.",
        "",
        "## Study Runs",
        "",
        *markdown_table(
            study_table,
            [
                "accession",
                "FLT",
                "GC",
                "ercc_dropped",
                "unique_wERCC_retained",
                "FLT_best_loss",
                "GC_best_loss",
            ],
            max_rows=20,
        ),
        "",
        "## Recurring GLARE Cluster Reactome Terms",
        "",
        f"Terms below are significant at cluster-level FDR < {alpha} and recur in at least two studies.",
        "",
    ]
    if recurring.empty:
        lines.append("No Reactome terms recurred across at least two per-study GLARE runs.")
    else:
        lines.extend(
            markdown_table(
                recurring,
                [
                    "location",
                    "term",
                    "study_count",
                    "accessions",
                    "cluster_count",
                    "best_fdr_bh",
                    "max_overlap",
                    "example_clusters",
                ],
                max_rows=40,
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Use these per-study GLARE terms as module-level support. The primary evidence remains per-study DESeq2 recurrence; GLARE is strongest when its recurring modules agree with recurring DESeq2 pathways.",
        ]
    )
    (output_dir / "PER_STUDY_GLARE_SUMMARY.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def run_study(args: argparse.Namespace, accession: str, input_dim: int, device) -> dict:
    study_dir = Path(args.output_dir) / accession
    prepared = select_aggregate_profiles(
        args.target_manifest,
        args.osdr_h5,
        [accession],
        study_dir,
        load_excluded_profiles(args.exclude_profiles_file, args.exclude_profile),
        args.ercc_policy,
    )
    counts = prepared["counts"].reset_index().rename(columns={"index": "h5_accession"})
    count_records = counts.reset_index(drop=True).to_dict(orient="records")
    flt_profiles = prepared["matrices"]["FLT"].shape[1]
    gc_profiles = prepared["matrices"]["GC"].shape[1]
    log(
        f"{accession}: prepared {len(prepared['genes'])} genes, "
        f"{flt_profiles} FLT and {gc_profiles} GC profiles"
    )
    if flt_profiles < args.min_profiles_per_condition or gc_profiles < args.min_profiles_per_condition:
        raise ValueError(
            f"{accession} has too few profiles after filters: "
            f"FLT={flt_profiles}, GC={gc_profiles}, "
            f"minimum={args.min_profiles_per_condition}"
        )
    if args.prepare_only:
        locations = []
    else:
        for location, matrix in prepared["matrices"].items():
            write_outlier_audit(matrix, prepared["genes"], location, study_dir)
        locations = [
            finetune_location(
                prepared["matrices"][location],
                prepared["genes"],
                location,
                Path(args.pretrained_weights),
                study_dir,
                device,
                input_dim,
                args.epochs,
                args.batch_size,
                args.seed,
            )
            for location in ("FLT", "GC")
        ]
        if not args.skip_clustering:
            run_clustering(study_dir, not args.with_tsne)

    summary = {
        "method": "per-study GLARE released 16-dimensional SAE with separate FLT/GC fine-tuning",
        "accession": accession,
        "target_manifest": prepared["target_manifest"],
        "target_expression_input": prepared["input_path"],
        "target_expression_kind": "aligned_osdr_liver_hdf5_expression",
        "osdr_h5": prepared["osdr_h5"],
        "pretrained_weights": args.pretrained_weights,
        "pretrained_input_dim": input_dim,
        "device": str(device),
        "seed_reused_for_each_location": args.seed,
        "architecture": [128, 64, 32, 16],
        "learning_rate": 1e-3,
        "weight_decay": 0,
        "sparsity_penalty": 1e-5,
        "batch_size": args.batch_size,
        "selection": {
            "accession": accession,
            "excluded_profiles_requested": prepared["excluded_profiles_requested"],
            "excluded_profiles_matched": prepared["excluded_profiles_matched"],
            "ercc_policy": prepared["ercc_policy"],
        },
        "condition_counts": count_records,
        "locations": locations,
    }
    (study_dir / "finetune_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", default=DEFAULT_TARGET_MANIFEST)
    parser.add_argument("--osdr-h5", default=DEFAULT_OSDR_H5)
    parser.add_argument("--studies", nargs="+", default=DEFAULT_STUDIES)
    parser.add_argument("--pretrained-weights", default=DEFAULT_PRETRAINED_WEIGHTS)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--exclude-profiles-file", default=DEFAULT_EXCLUDE_PROFILES)
    parser.add_argument("--exclude-profile", action="append", default=[])
    parser.add_argument(
        "--ercc-policy",
        choices=["keep_all", "prefer_noercc"],
        default="prefer_noercc",
    )
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1996)
    parser.add_argument("--min-profiles-per-condition", type=int, default=5)
    parser.add_argument("--reactome-gmt", default=DEFAULT_REACTOME_GMT)
    parser.add_argument("--enrichment-alpha", type=float, default=0.05)
    parser.add_argument("--min-overlap", type=int, default=5)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--skip-clustering", action="store_true")
    parser.add_argument(
        "--with-tsne",
        action="store_true",
        help="Also run t-SNE visualizations. PCA plots are always generated.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start = time.perf_counter()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pretrained_weights = Path(args.pretrained_weights)
    input_dim = infer_pretrain_input_dim(pretrained_weights)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    study_summaries = []
    for accession in args.studies:
        study_summaries.append(run_study(args, accession, input_dim, device))

    enrichment = pd.DataFrame()
    recurrence = pd.DataFrame()
    if not args.prepare_only and not args.skip_clustering:
        enrichment_tables = [
            cluster_enrichment_for_run(Path(args.output_dir) / accession, Path(args.reactome_gmt), args.min_overlap)
            for accession in args.studies
        ]
        enrichment_tables = [table for table in enrichment_tables if not table.empty]
        if enrichment_tables:
            enrichment = pd.concat(enrichment_tables, ignore_index=True)
            enrichment.to_csv(
                output_dir / "glare_cluster_reactome_enrichment.tsv",
                sep="\t",
                index=False,
            )
            recurrence = summarize_recurring_terms(enrichment, args.enrichment_alpha)
            recurrence.to_csv(
                output_dir / "recurring_glare_reactome_terms.tsv",
                sep="\t",
                index=False,
            )

    summary = {
        "studies": args.studies,
        "output_dir": str(output_dir),
        "pretrained_weights": args.pretrained_weights,
        "ercc_policy": args.ercc_policy,
        "exclude_profiles_file": args.exclude_profiles_file,
        "min_profiles_per_condition": args.min_profiles_per_condition,
        "study_summaries": [str(output_dir / accession / "finetune_summary.json") for accession in args.studies],
        "enrichment": str(output_dir / "glare_cluster_reactome_enrichment.tsv") if not enrichment.empty else "",
        "recurrence": str(output_dir / "recurring_glare_reactome_terms.tsv") if not recurrence.empty else "",
        "elapsed_seconds": round(time.perf_counter() - start, 3),
        "elapsed": format_elapsed(time.perf_counter() - start),
    }
    (output_dir / "per_study_glare_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    write_report(output_dir, study_summaries, recurrence, args.enrichment_alpha)
    log(f"Saved per-study GLARE summary: {output_dir / 'PER_STUDY_GLARE_SUMMARY.md'}")


if __name__ == "__main__":
    main()
