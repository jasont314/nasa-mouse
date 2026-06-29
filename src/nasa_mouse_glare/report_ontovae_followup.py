"""Build a focused OntoVAE follow-up report from completed output files."""

from __future__ import annotations

import argparse
from pathlib import Path

from .io import require_import


DEFAULT_SUMMARY = Path("outputs/ontovae_pipeline/summary/ontovae_run_summary.tsv")
DEFAULT_OUTPUT_DIR = Path("outputs/ontovae_pipeline/followup")
DEFAULT_DOC = Path("docs/ontovae_followup_report.md")

PLOT_FILES = {
    "pca_condition": "pathway_score_pca.png",
    "pca_accession": "pathway_score_pca_by_accession.png",
    "umap_condition": "pathway_score_umap.png",
    "umap_accession": "pathway_score_umap_by_accession.png",
    "heatmap": "top_pathway_shift_heatmap.png",
}

PRIORITY_KEYS = [
    ("skeletal_muscle", "soleus", "archs4_pretrain_osdr_finetune", "finetuned_or_direct"),
    ("liver", "", "archs4_pretrain_osdr_finetune", "finetuned_or_direct"),
    ("skeletal_muscle", "", "archs4_pretrain_osdr_finetune", "finetuned_or_direct"),
    ("skeletal_muscle", "quadriceps", "archs4_pretrain_osdr_finetune", "finetuned_or_direct"),
    ("spleen", "", "archs4_pretrain_osdr_finetune", "finetuned_or_direct"),
    ("thymus", "", "archs4_pretrain_osdr_finetune", "finetuned_or_direct"),
]


def read_tsv(path: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    return pd.read_csv(path, sep="\t")


def clean_group(value) -> str:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    if pd.isna(value):
        return ""
    return str(value)


def score_dirs(run_dir: Path, score_set: str) -> tuple[Path, Path]:
    if score_set == "pre_finetune_projection":
        return run_dir / "pretrained_query_analysis", run_dir / "pretrained_query_accession_validation"
    return run_dir / "analysis", run_dir / "accession_validation"


def stable_terms_for_row(row) -> list[dict]:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    run_dir = Path(row["output_dir"])
    analysis_dir, validation_dir = score_dirs(run_dir, row["score_set"])
    meta_path = validation_dir / "random_effects_meta_analysis.tsv"
    loo_path = validation_dir / "leave_one_out_summary.tsv"
    if not meta_path.exists() or not loo_path.exists():
        return []

    meta = read_tsv(meta_path)
    loo = read_tsv(loo_path)
    if meta.empty or loo.empty:
        return []
    merged = meta.merge(loo, on="term", how="inner")
    if merged.empty:
        return []

    stable = merged.loc[
        (merged["meta_fdr"] < 0.05)
        & (merged["maximum_leave_one_out_fdr"] < 0.05)
        & (merged["n_same_direction"] == merged["n_leave_one_out"])
    ].copy()
    if stable.empty:
        return []

    terms_path = run_dir / "terms.tsv"
    if terms_path.exists():
        terms = read_tsv(terms_path)
        stable = stable.merge(terms, on="term", how="left")
    else:
        stable["description"] = ""
        stable["n_genes"] = pd.NA

    stable = stable.sort_values(["meta_fdr", "maximum_leave_one_out_fdr", "term"])
    rows = []
    for rank, (_, item) in enumerate(stable.iterrows(), start=1):
        plot_paths = {
            name: str(analysis_dir / filename)
            for name, filename in PLOT_FILES.items()
            if (analysis_dir / filename).exists()
        }
        rows.append(
            {
                "tissue": row["tissue"],
                "group": clean_group(row.get("group", "")),
                "mode": row["mode"],
                "score_set": row["score_set"],
                "stable_rank": rank,
                "term": item["term"],
                "description": item.get("description", ""),
                "n_genes": item.get("n_genes", ""),
                "meta_effect": item["meta_effect"],
                "meta_fdr": item["meta_fdr"],
                "meta_p": item.get("meta_p", ""),
                "tau2": item.get("tau2", ""),
                "i2": item.get("i2", ""),
                "n_accessions": item.get("n_accessions", ""),
                "n_accession_same_direction": item.get("n_accession_same_direction", ""),
                "n_accession_opposite_direction": item.get(
                    "n_accession_opposite_direction", ""
                ),
                "n_leave_one_out": item["n_leave_one_out"],
                "n_same_direction": item["n_same_direction"],
                "minimum_leave_one_out_fdr": item["minimum_leave_one_out_fdr"],
                "maximum_leave_one_out_fdr": item["maximum_leave_one_out_fdr"],
                "query_samples": row.get("query_samples", ""),
                "reference_samples": row.get("reference_samples", ""),
                "genes": row.get("genes", ""),
                "terms": row.get("terms", ""),
                "output_dir": str(run_dir),
                "analysis_dir": str(analysis_dir),
                "validation_dir": str(validation_dir),
                **plot_paths,
            }
        )
    return rows


def collect_stable_terms(summary_path: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    summary = read_tsv(summary_path)
    if "group" in summary.columns:
        summary["group"] = summary["group"].map(clean_group)
    stable_rows: list[dict] = []
    for _, row in summary.loc[summary["loo_stable_fdr_lt_005"] > 0].iterrows():
        stable_rows.extend(stable_terms_for_row(row))
    return pd.DataFrame(stable_rows), summary


def collect_top_genes(stable_terms, *, terms_per_score_set: int, genes_per_term: int):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    if stable_terms.empty:
        return pd.DataFrame()
    rows = []
    grouped = stable_terms.groupby(["output_dir", "score_set"], sort=False)
    for (_, score_set), group in grouped:
        selected = group.sort_values("stable_rank").head(terms_per_score_set)
        for _, term_row in selected.iterrows():
            weights_path = Path(term_row["output_dir"]) / "term_gene_weights_top.tsv"
            if not weights_path.exists():
                continue
            weights = read_tsv(weights_path)
            subset = weights.loc[weights["term"].eq(term_row["term"])].head(genes_per_term)
            for _, gene_row in subset.iterrows():
                rows.append(
                    {
                        "tissue": term_row["tissue"],
                        "group": term_row["group"],
                        "mode": term_row["mode"],
                        "score_set": score_set,
                        "stable_rank": int(term_row["stable_rank"]),
                        "term": term_row["term"],
                        "meta_effect": term_row["meta_effect"],
                        "meta_fdr": term_row["meta_fdr"],
                        "maximum_leave_one_out_fdr": term_row[
                            "maximum_leave_one_out_fdr"
                        ],
                        "gene_rank": int(gene_row["rank"]),
                        "gene": gene_row["gene"],
                        "decoder_weight": gene_row["decoder_weight"],
                        "output_dir": term_row["output_dir"],
                    }
                )
    return pd.DataFrame(rows)


def rel(path: str | Path, base: Path = Path(".")) -> str:
    path = Path(path)
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def md_link(path: str | Path, label: str | None = None, *, from_docs: bool = False) -> str:
    path = Path(path)
    target = Path("..") / path if from_docs else path
    return f"[{label or path.name}]({target.as_posix()})"


def format_float(value, digits: int = 3) -> str:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    if pd.isna(value):
        return ""
    value = float(value)
    if value != 0 and abs(value) < 0.001:
        return f"{value:.2e}"
    return f"{value:.{digits}g}"


def key_tuple(row) -> tuple[str, str, str, str]:
    return (
        str(row["tissue"]),
        clean_group(row.get("group", "")),
        str(row["mode"]),
        str(row["score_set"]),
    )


def write_markdown(
    *,
    stable_terms,
    top_genes,
    summary,
    output_dir: Path,
    doc_path: Path,
) -> None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    final_stable = stable_terms.loc[
        stable_terms["score_set"].eq("finetuned_or_direct")
    ].copy()
    score_set_counts = (
        stable_terms.groupby(["tissue", "group", "mode", "score_set"], dropna=False)
        .size()
        .reset_index(name="loo_stable_terms")
    )
    first_terms = (
        stable_terms.sort_values(["tissue", "group", "mode", "score_set", "stable_rank"])
        .groupby(["tissue", "group", "mode", "score_set"], dropna=False)
        .head(1)
    )
    score_set_counts = score_set_counts.merge(
        first_terms[
            [
                "tissue",
                "group",
                "mode",
                "score_set",
                "term",
                "meta_effect",
                "meta_fdr",
                "maximum_leave_one_out_fdr",
            ]
        ],
        on=["tissue", "group", "mode", "score_set"],
        how="left",
    )

    priority_frames = []
    for key in PRIORITY_KEYS:
        mask = stable_terms.apply(lambda row: key_tuple(row) == key, axis=1)
        if mask.any():
            priority_frames.append(stable_terms.loc[mask].sort_values("stable_rank").head(5))
    priority = pd.concat(priority_frames, ignore_index=True) if priority_frames else pd.DataFrame()

    lines = [
        "# OntoVAE Follow-up Report",
        "",
        "This report is generated from completed OntoVAE outputs. It focuses on",
        "strict leave-one-accession-out stable FLT vs GC pathway/program shifts,",
        "then links directly to the PCA/UMAP/heatmap visualizations for manual",
        "inspection.",
        "",
        "Primary machine-readable outputs:",
        "",
        f"- `{rel(output_dir / 'ontovae_stable_followup_terms.tsv')}`",
        f"- `{rel(output_dir / 'ontovae_stable_followup_top_genes.tsv')}`",
        "",
        "Strict LOO-stable means `meta_fdr < 0.05`,",
        "`maximum_leave_one_out_fdr < 0.05`, and all leave-one-accession-out",
        "effects keep the same direction.",
        "",
        "## Focus Runs",
        "",
    ]

    focus_rows = score_set_counts.loc[
        score_set_counts["score_set"].eq("finetuned_or_direct")
    ].sort_values(
        ["tissue", "group", "mode"]
    )
    if not focus_rows.empty:
        lines.extend(
            [
                "| tissue | group | mode | stable terms | top stable term | effect | FDR | max LOO FDR |",
                "| --- | --- | --- | ---: | --- | ---: | ---: | ---: |",
            ]
        )
        for _, row in focus_rows.iterrows():
            lines.append(
                "| {tissue} | {group} | {mode} | {n} | `{term}` | {effect} | {fdr} | {loo} |".format(
                    tissue=row["tissue"],
                    group=row["group"],
                    mode=row["mode"],
                    n=int(row["loo_stable_terms"]),
                    term=row["term"],
                    effect=format_float(row["meta_effect"]),
                    fdr=format_float(row["meta_fdr"]),
                    loo=format_float(row["maximum_leave_one_out_fdr"]),
                )
            )
    else:
        lines.append("No strict LOO-stable final OntoVAE rows were found.")

    lines.extend(["", "## Priority Pathways", ""])
    if not priority.empty:
        lines.extend(
            [
                "| tissue | group | mode | score set | rank | term | effect | FDR | max LOO FDR | top decoder genes |",
                "| --- | --- | --- | --- | ---: | --- | ---: | ---: | ---: | --- |",
            ]
        )
        for _, row in priority.iterrows():
            genes = top_genes.loc[
                (top_genes["output_dir"].eq(row["output_dir"]))
                & (top_genes["score_set"].eq(row["score_set"]))
                & (top_genes["term"].eq(row["term"]))
            ].sort_values("gene_rank")
            gene_text = ", ".join(genes["gene"].head(5).astype(str).tolist())
            lines.append(
                "| {tissue} | {group} | {mode} | {score_set} | {rank} | `{term}` | {effect} | {fdr} | {loo} | {genes} |".format(
                    tissue=row["tissue"],
                    group=row["group"],
                    mode=row["mode"],
                    score_set=row["score_set"],
                    rank=int(row["stable_rank"]),
                    term=row["term"],
                    effect=format_float(row["meta_effect"]),
                    fdr=format_float(row["meta_fdr"]),
                    loo=format_float(row["maximum_leave_one_out_fdr"]),
                    genes=gene_text,
                )
            )
    else:
        lines.append("No priority pathways were found.")

    lines.extend(["", "## Plot Review Queue", ""])
    review = []
    for key in PRIORITY_KEYS:
        mask = stable_terms.apply(lambda row: key_tuple(row) == key, axis=1)
        if mask.any():
            review.append(stable_terms.loc[mask].sort_values("stable_rank").iloc[0])
    for row in review:
        lines.extend(
            [
                f"### {row['tissue']} {row['group']} {row['mode']} {row['score_set']}".strip(),
                "",
                f"Top stable term: `{row['term']}` "
                f"(effect {format_float(row['meta_effect'])}, "
                f"FDR {format_float(row['meta_fdr'])}, "
                f"max LOO FDR {format_float(row['maximum_leave_one_out_fdr'])}).",
                "",
                "- "
                + md_link(row["heatmap"], "top pathway heatmap", from_docs=True)
                + " | "
                + md_link(row["umap_accession"], "UMAP by accession", from_docs=True)
                + " | "
                + md_link(row["pca_accession"], "PCA by accession", from_docs=True),
                "",
                f"![Top pathway heatmap](../{row['heatmap']})",
                "",
            ]
        )

    lines.extend(
        [
            "## Frozen Projection vs Fine-tuned Scores",
            "",
            "Rows with `score_set = pre_finetune_projection` are OSDR samples scored",
            "by the ARCHS4-pretrained OntoVAE before OSDR fine-tuning. Rows with",
            "`finetuned_or_direct` are final fine-tuned or direct OSDR scores.",
            "",
            "| tissue | group | mode | score set | stable terms | top stable term | effect | FDR |",
            "| --- | --- | --- | --- | ---: | --- | ---: | ---: |",
        ]
    )
    for _, row in score_set_counts.sort_values(
        ["tissue", "group", "mode", "score_set"]
    ).iterrows():
        lines.append(
            "| {tissue} | {group} | {mode} | {score_set} | {n} | `{term}` | {effect} | {fdr} |".format(
                tissue=row["tissue"],
                group=row["group"],
                mode=row["mode"],
                score_set=row["score_set"],
                n=int(row["loo_stable_terms"]),
                term=row["term"],
                effect=format_float(row["meta_effect"]),
                fdr=format_float(row["meta_fdr"]),
            )
        )

    lines.extend(
        [
            "",
            "## Caveats",
            "",
            "- Decoder genes are listed as Ensembl mouse gene IDs because the OntoVAE",
            "  AnnData inputs do not carry a project-wide gene-symbol annotation.",
            "- These are pathway/program score shifts, not direct gene-level DGEA.",
            "- OntoVAE here uses ARCHS4 pretraining plus OSDR fine-tuning; it is not",
            "  native scArches query mapping.",
            "- Skin, kidney, EDL, gastrocnemius, and tibialis anterior had exploratory",
            "  random-effects hits but no strict LOO-stable final hits.",
        ]
    )

    doc_path.parent.mkdir(parents=True, exist_ok=True)
    doc_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    readme = output_dir / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# OntoVAE Follow-up Outputs",
                "",
                "Generated focused follow-up tables for strict LOO-stable OntoVAE",
                "FLT vs GC pathway/program shifts.",
                "",
                f"- Stable pathway table: `{(output_dir / 'ontovae_stable_followup_terms.tsv').name}`",
                f"- Top decoder genes: `{(output_dir / 'ontovae_stable_followup_top_genes.tsv').name}`",
                f"- Rendered report: `{doc_path}`",
                "",
                f"Stable term rows: {len(stable_terms)}",
                f"Top-gene rows: {len(top_genes)}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def run(args) -> Path:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stable_terms, summary = collect_stable_terms(Path(args.summary))
    top_genes = collect_top_genes(
        stable_terms,
        terms_per_score_set=args.terms_per_score_set,
        genes_per_term=args.genes_per_term,
    )

    terms_path = output_dir / "ontovae_stable_followup_terms.tsv"
    genes_path = output_dir / "ontovae_stable_followup_top_genes.tsv"
    stable_terms.to_csv(terms_path, sep="\t", index=False)
    top_genes.to_csv(genes_path, sep="\t", index=False)
    write_markdown(
        stable_terms=stable_terms,
        top_genes=top_genes,
        summary=summary,
        output_dir=output_dir,
        doc_path=Path(args.doc),
    )
    print(
        {
            "stable_terms": len(stable_terms),
            "top_gene_rows": len(top_genes),
            "terms": str(terms_path),
            "top_genes": str(genes_path),
            "doc": str(args.doc),
        }
    )
    return terms_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--doc", type=Path, default=DEFAULT_DOC)
    parser.add_argument("--terms-per-score-set", type=int, default=5)
    parser.add_argument("--genes-per-term", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
