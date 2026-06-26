"""Summarize expiMap OSDR and ARCHS4 reference-query runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import require_import


DEFAULT_TISSUES = ("liver", "kidney")
DIRECT_RUNS = {
    "raw_counts": {
        "dir": "raw_counts_nb_50epoch",
        "status": "primary",
        "note": "raw NASA OSDR API counts with negative-binomial loss",
    },
    "cpm": {
        "dir": "cpm_mse_50epoch",
        "status": "sensitivity",
        "note": "CPM-normalized values with MSE loss",
    },
    "log1p_cpm": {
        "dir": "log1p_cpm_mse_50epoch",
        "status": "sensitivity",
        "note": "log1p(CPM) values with MSE loss",
    },
}


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_table(path: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def relpath(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root))
    except ValueError:
        return str(path)


def markdown_link(path: Path | str) -> str:
    return f"`{path}`"


def format_number(value) -> str:
    if value is None:
        return "NA"
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if numeric != numeric:
        return "NA"
    if abs(numeric) >= 1000:
        return f"{numeric:,.0f}"
    if abs(numeric) >= 1:
        return f"{numeric:.4g}"
    return f"{numeric:.3g}"


def fdr_summary(comparison) -> dict:
    if comparison.empty:
        return {
            "top_term": "NA",
            "top_effect": None,
            "min_welch_fdr": None,
            "min_mannwhitney_fdr": None,
            "n_welch_fdr_0_10": 0,
            "n_mannwhitney_fdr_0_10": 0,
        }
    top = comparison.iloc[0]
    return {
        "top_term": str(top.get("term", "NA")),
        "top_effect": top.get("flight_minus_ground"),
        "min_welch_fdr": comparison["welch_fdr"].astype(float).min()
        if "welch_fdr" in comparison
        else None,
        "min_mannwhitney_fdr": comparison["mannwhitney_fdr"].astype(float).min()
        if "mannwhitney_fdr" in comparison
        else None,
        "n_welch_fdr_0_10": int((comparison["welch_fdr"].astype(float) < 0.10).sum())
        if "welch_fdr" in comparison
        else 0,
        "n_mannwhitney_fdr_0_10": int(
            (comparison["mannwhitney_fdr"].astype(float) < 0.10).sum()
        )
        if "mannwhitney_fdr" in comparison
        else 0,
    }


def direct_tissue_summary(root: Path, tissue: str) -> dict:
    tissue_dir = root / f"outputs/expimap_direct_osdr_{tissue}"
    input_manifest = read_json(tissue_dir / "input/input_manifest.json")
    runs = {}
    for transform, spec in DIRECT_RUNS.items():
        run_dir = tissue_dir / spec["dir"]
        training = read_json(run_dir / "training_summary.json")
        analysis = read_json(run_dir / "analysis/analysis_summary.json")
        comparison = read_table(run_dir / "analysis/flt_vs_gc_pathway_comparison.tsv")
        study_tests = read_table(run_dir / "analysis/flight_ground_study_aware_tests.tsv")
        run_summary = {
            "path": relpath(root, run_dir),
            "transformation": transform,
            "status": spec["status"],
            "note": spec["note"],
            "exists": bool(training and analysis and not comparison.empty),
            "training": training,
            "analysis": analysis,
            **fdr_summary(comparison),
        }
        if not study_tests.empty and "wilcoxon_fdr" in study_tests:
            run_summary["min_study_aware_fdr"] = (
                study_tests["wilcoxon_fdr"].astype(float).min()
            )
            run_summary["n_study_aware_fdr_0_10"] = int(
                (study_tests["wilcoxon_fdr"].astype(float) < 0.10).sum()
            )
        else:
            run_summary["min_study_aware_fdr"] = None
            run_summary["n_study_aware_fdr_0_10"] = 0
        runs[transform] = run_summary
    preprocessing_dir = tissue_dir / "preprocessing_comparison_50epoch"
    return {
        "tissue": tissue,
        "path": relpath(root, tissue_dir),
        "input_manifest": input_manifest,
        "runs": runs,
        "preprocessing_summary": read_table(
            preprocessing_dir / "preprocessing_run_summary.tsv"
        ),
        "preprocessing_correlations": read_table(
            preprocessing_dir / "preprocessing_effect_correlations.tsv"
        ),
        "preprocessing_manifest": read_json(
            preprocessing_dir / "preprocessing_comparison_manifest.json"
        ),
    }


def reference_tissue_summary(root: Path, tissue: str) -> dict:
    tissue_dir = root / f"outputs/expimap_archs4_reference_osdr_query_{tissue}"
    reference_input = read_json(tissue_dir / "reference_input_1000/reference_input_manifest.json")
    reference_training = read_json(tissue_dir / "reference_nb_1000_50epoch/training_summary.json")
    query_mapping = read_json(tissue_dir / "query_nb_1000ref_50epoch/query_mapping_summary.json")
    analysis = read_json(tissue_dir / "query_nb_1000ref_50epoch/analysis/analysis_summary.json")
    comparison = read_table(tissue_dir / "query_nb_1000ref_50epoch/analysis/flt_vs_gc_pathway_comparison.tsv")
    study_tests = read_table(tissue_dir / "query_nb_1000ref_50epoch/analysis/flight_ground_study_aware_tests.tsv")
    summary = {
        "tissue": tissue,
        "path": relpath(root, tissue_dir),
        "reference_input": reference_input,
        "reference_training": reference_training,
        "query_mapping": query_mapping,
        "analysis": analysis,
        **fdr_summary(comparison),
    }
    if not study_tests.empty and "wilcoxon_fdr" in study_tests:
        summary["min_study_aware_fdr"] = study_tests["wilcoxon_fdr"].astype(float).min()
        summary["n_study_aware_fdr_0_10"] = int(
            (study_tests["wilcoxon_fdr"].astype(float) < 0.10).sum()
        )
    else:
        summary["min_study_aware_fdr"] = None
        summary["n_study_aware_fdr_0_10"] = 0
    return summary


def write_table(lines: list[str], headers: list[str], rows: list[list[str]]) -> None:
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")


def add_osdr_discovery(lines: list[str], root: Path) -> None:
    tissue_counts = read_table(
        root / "data/osdr_api/osdr_api_mouse_bulk_rnaseq_tissue_counts.tsv"
    )
    summary = read_json(root / "data/osdr_api/osdr_api_mouse_bulk_rnaseq_summary.json")
    lines.extend(
        [
            "## OSDR API Discovery",
            "",
            "OSDR inputs are built from the NASA OSDR Biological Data API, not from the older local integrated OSDR HDF5.",
            "",
            f"- API-selected samples: {format_number(summary.get('counts', {}).get('samples'))}",
            f"- OSD accessions: {format_number(summary.get('counts', {}).get('accessions'))}",
            f"- count files: {format_number(summary.get('counts', {}).get('files'))}",
            f"- tissues: {format_number(summary.get('counts', {}).get('tissues'))}",
            "",
        ]
    )
    if not tissue_counts.empty:
        rows = []
        for _, row in tissue_counts.head(12).iterrows():
            rows.append(
                [
                    str(row["tissue_final"]),
                    format_number(row.get("flight")),
                    format_number(row.get("ground_control")),
                    format_number(row.get("total_flt_gc")),
                ]
            )
        write_table(lines, ["tissue", "flight", "ground_control", "total"], rows)
        extra = tissue_counts.loc[
            ~tissue_counts["tissue_final"].isin(DEFAULT_TISSUES)
        ].head(8)
        if not extra.empty:
            rows = []
            for _, row in extra.iterrows():
                rows.append(
                    [
                        str(row["tissue_final"]),
                        format_number(row.get("flight")),
                        format_number(row.get("ground_control")),
                        format_number(row.get("total_flt_gc")),
                    ]
                )
            lines.extend(["Top additional tissues by available FLT/GC samples:", ""])
            write_table(lines, ["tissue", "flight", "ground_control", "total"], rows)


def add_direct_section(lines: list[str], direct: dict[str, dict]) -> None:
    lines.extend(
        [
            "## Direct OSDR expiMap",
            "",
            "These runs train expiMap directly on each tissue-specific OSDR FLT/GC dataset.",
            "The current direct summaries use matched 50-epoch runs. Raw counts with negative-binomial loss are the primary analysis; CPM and log1p(CPM) use MSE loss and are sensitivity checks.",
            "",
        ]
    )
    rows = []
    for tissue, summary in direct.items():
        counts = summary["input_manifest"].get("counts", {})
        rows.append(
            [
                tissue,
                format_number(counts.get("samples")),
                format_number(counts.get("flight")),
                format_number(counts.get("ground_control")),
                format_number(counts.get("accessions")),
                format_number(counts.get("genes")),
                format_number(counts.get("pathways")),
                markdown_link(Path(summary["path"]) / "input"),
            ]
        )
    write_table(
        lines,
        ["tissue", "samples", "flight", "ground", "accessions", "genes", "pathways", "input"],
        rows,
    )

    rows = []
    for tissue, summary in direct.items():
        for transform, run in summary["runs"].items():
            training = run["training"]
            rows.append(
                [
                    tissue,
                    transform,
                    run["status"],
                    str(training.get("recon_loss", "missing")),
                    format_number(training.get("epochs")),
                    format_number(run["min_welch_fdr"]),
                    format_number(run["n_welch_fdr_0_10"]),
                    format_number(run["min_study_aware_fdr"]),
                    format_number(run["n_study_aware_fdr_0_10"]),
                    run["top_term"],
                ]
            )
    write_table(
        lines,
        [
            "tissue",
            "transform",
            "role",
            "loss",
            "epochs",
            "min Welch FDR",
            "Welch FDR<0.10",
            "min study FDR",
            "study FDR<0.10",
            "top aggregate term",
        ],
        rows,
    )
    lines.extend(
        [
            "Current direct-run interpretation:",
            "",
            "- Each current direct analysis directory contains PCA and UMAP coordinates plus condition/accession-colored plots.",
            "- Liver has one aggregate FDR-significant pathway in the 50-epoch direct runs: RNA Polymerase II transcription elongation is lower in flight and is the top term across raw counts, CPM, and log1p(CPM).",
            "- Kidney has no aggregate pathway FDR < 0.10 in the matched 50-epoch direct runs.",
            "- 50-epoch effect rankings are highly correlated across raw-count NB and CPM/log1p(CPM) MSE sensitivity runs, unlike the earlier 3-epoch validation runs.",
            "- Study-aware accession-level tests should still be treated as exploratory unless they agree with aggregate and preprocessing-stable effects.",
            "",
        ]
    )


def add_preprocessing_section(lines: list[str], direct: dict[str, dict]) -> None:
    lines.extend(["## Preprocessing Comparison", ""])
    for tissue, summary in direct.items():
        lines.extend([f"### {tissue}", ""])
        run_table = summary["preprocessing_summary"]
        if not run_table.empty:
            rows = []
            for _, row in run_table.iterrows():
                rows.append(
                    [
                        str(row.get("transformation")),
                        str(row.get("recon_loss")),
                        str(row.get("validity")),
                        format_number(row.get("min_welch_fdr")),
                        format_number(row.get("n_welch_fdr_0_10")),
                        str(row.get("top_welch_term")),
                    ]
                )
            write_table(
                lines,
                ["transform", "loss", "validity", "min Welch FDR", "Welch FDR<0.10", "top term"],
                rows,
            )
        corr_table = summary["preprocessing_correlations"]
        if not corr_table.empty:
            rows = []
            for _, row in corr_table.iterrows():
                rows.append(
                    [
                        f"{row.get('left_transformation')} vs {row.get('right_transformation')}",
                        format_number(row.get("spearman_effect_rho")),
                        format_number(row.get("top50_overlap")),
                        format_number(row.get("top50_jaccard")),
                    ]
                )
            write_table(
                lines,
                ["comparison", "Spearman effect rho", "top50 overlap", "top50 Jaccard"],
                rows,
            )


def add_reference_section(lines: list[str], reference: dict[str, dict]) -> None:
    lines.extend(
        [
            "## ARCHS4 Reference-Query expiMap",
            "",
            "These runs train a tissue-filtered ARCHS4 mouse bulk reference and map the OSDR tissue dataset as query.",
            "The current liver/kidney reference-query summaries use 1000 ARCHS4 tissue-filtered reference samples, 50 reference-training epochs, and 50 query-mapping epochs.",
            "",
        ]
    )
    rows = []
    for tissue, summary in reference.items():
        ref_counts = summary["reference_input"].get("counts", {})
        query = summary["query_mapping"]
        reference_training = summary["reference_training"]
        rows.append(
            [
                tissue,
                format_number(ref_counts.get("samples")),
                format_number(reference_training.get("epochs")),
                format_number(ref_counts.get("genes")),
                format_number(ref_counts.get("pathways")),
                format_number(query.get("n_query_samples")),
                format_number(query.get("n_query_genes")),
                format_number(query.get("epochs")),
                format_number(summary["min_welch_fdr"]),
                format_number(summary["n_welch_fdr_0_10"]),
                summary["top_term"],
            ]
        )
    write_table(
        lines,
        [
            "tissue",
            "ARCHS4 ref samples",
            "ref epochs",
            "ref genes",
            "pathways",
            "query samples",
            "query genes",
            "query epochs",
            "min Welch FDR",
            "Welch FDR<0.10",
            "top query term",
        ],
        rows,
    )
    lines.extend(
        [
            "Reference-query interpretation:",
            "",
            "- These runs use tissue-filtered, leakage-excluded ARCHS4 references with the same Reactome architecture as the direct OSDR runs.",
            "- Reference-query preprocessing is raw-count NB only in the current workflow; CPM/log1p(CPM) comparisons are direct-workflow sensitivity analyses.",
            "- Each current query analysis directory contains PCA and UMAP coordinates plus condition/accession-colored plots.",
            "- They are still bounded 1000-sample reference runs rather than all available ARCHS4 tissue samples, so compare them with direct OSDR results before treating a pathway as robust.",
            "",
        ]
    )


def add_workflow_agreement(lines: list[str], root: Path, direct: dict[str, dict]) -> None:
    lines.extend(["## Direct vs Reference-Query Agreement", ""])
    rows = []
    for tissue, summary in direct.items():
        term = summary["runs"]["raw_counts"]["top_term"]
        direct_path = (
            root
            / f"outputs/expimap_direct_osdr_{tissue}/raw_counts_nb_50epoch/analysis/flt_vs_gc_pathway_comparison.tsv"
        )
        reference_path = (
            root
            / f"outputs/expimap_archs4_reference_osdr_query_{tissue}/query_nb_1000ref_50epoch/analysis/flt_vs_gc_pathway_comparison.tsv"
        )
        direct_table = read_table(direct_path)
        reference_table = read_table(reference_path)
        direct_row = direct_table.loc[direct_table["term"].eq(term)].head(1)
        reference_row = reference_table.loc[reference_table["term"].eq(term)].head(1)
        if direct_row.empty or reference_row.empty:
            continue
        direct_effect = float(direct_row["flight_minus_ground"].iloc[0])
        reference_effect = float(reference_row["flight_minus_ground"].iloc[0])
        direct_fdr = float(direct_row["welch_fdr"].iloc[0])
        reference_fdr = float(reference_row["welch_fdr"].iloc[0])
        rows.append(
            [
                tissue,
                term,
                format_number(direct_effect),
                format_number(direct_fdr),
                format_number(reference_effect),
                format_number(reference_fdr),
                "yes" if direct_effect * reference_effect > 0 else "no",
                "yes" if reference_fdr < 0.10 else "no",
            ]
        )
    if rows:
        write_table(
            lines,
            [
                "tissue",
                "direct raw-count top term",
                "direct effect",
                "direct Welch FDR",
                "reference-query effect",
                "reference-query Welch FDR",
                "same direction",
                "reference FDR<0.10",
            ],
            rows,
        )
    lines.extend(
        [
            "Current workflow-agreement interpretation:",
            "",
            "- The direct liver top term has the same negative direction in the bounded ARCHS4 reference-query run, but it is not FDR-significant there.",
            "- Kidney has no aggregate FDR-significant direct or reference-query pathway signal.",
            "- Treat the direct liver signal as preprocessing-stable but not yet reference-query-confirmed.",
            "",
        ]
    )


def add_archs4_section(lines: list[str], root: Path) -> None:
    archs4 = read_table(root / "data/archs4/archs4_mouse_tissue_summary.tsv")
    if archs4.empty:
        return
    lines.extend(["## ARCHS4 Tissue Candidates", ""])
    rows = []
    for _, row in archs4.iterrows():
        rows.append(
            [
                str(row.get("tissue")),
                format_number(row.get("usable_nonleakage_bulk_like_samples")),
            ]
        )
    write_table(lines, ["tissue", "usable nonleakage bulk-like samples"], rows)


def build_report(root: Path, tissues: tuple[str, ...]) -> tuple[str, dict]:
    direct = {tissue: direct_tissue_summary(root, tissue) for tissue in tissues}
    reference = {tissue: reference_tissue_summary(root, tissue) for tissue in tissues}
    lines = [
        "# expiMap Results Summary",
        "",
        "Generated from current local manifests and ignored output artifacts.",
        "This report intentionally distinguishes pipeline validation from final biological inference.",
        "",
    ]
    add_osdr_discovery(lines, root)
    add_direct_section(lines, direct)
    add_preprocessing_section(lines, direct)
    add_reference_section(lines, reference)
    add_workflow_agreement(lines, root, direct)
    add_archs4_section(lines, root)
    lines.extend(
        [
            "## Current Bottom Line",
            "",
            "- API-derived OSDR liver and kidney direct expiMap analyses are implemented and runnable.",
            "- The matched 50-epoch direct liver runs nominate lower flight RNA Polymerase II transcription elongation as the most robust current signal.",
            "- The matched 50-epoch direct kidney runs do not show aggregate FDR-significant FLT-vs-GC pathway shifts.",
            "- CPM/log1p(CPM) sensitivity results are now rank-stable against raw-count NB for the matched 50-epoch direct runs.",
            "- ARCHS4 reference-query is implemented for liver and kidney with bounded 1000-sample tissue-filtered references; it does not currently confirm the direct liver signal at aggregate FDR < 0.10.",
            "",
            "## Next Full-Run Gate",
            "",
            "Before making stronger scientific claims, run larger ARCHS4 reference subsets or all available nonleakage tissue samples, then regenerate this report.",
            "",
        ]
    )
    summary = {
        "tissues": list(tissues),
        "direct": direct,
        "reference_query": reference,
    }
    return "\n".join(lines), summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a markdown summary of expiMap direct and reference-query runs."
    )
    parser.add_argument("--root", default=".")
    parser.add_argument("--output", default="docs/expimap_results.md")
    parser.add_argument("--json-output", default="docs/expimap_results_summary.json")
    parser.add_argument("--tissue", action="append", choices=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    tissues = tuple(args.tissue or DEFAULT_TISSUES)
    report, summary = build_report(root, tissues)
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(report + "\n", encoding="utf-8")
    json_output = root / args.json_output
    json_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.write_text(json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(output), "summary": str(json_output)}, indent=2))


if __name__ == "__main__":
    main()
