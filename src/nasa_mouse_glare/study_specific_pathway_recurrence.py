"""Study-first FLT-vs-GC pathway recurrence analysis for OSDR liver data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from .io import require_import


DEFAULT_OUTPUT_DIR = (
    "outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/"
    "post_analysis/study_specific_pathway_recurrence"
)
DEFAULT_H5 = "assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5"
DEFAULT_LIVER_METADATA = "data/processed/osdr_mouse_bulk_liver.profile_metadata.tsv"
DEFAULT_EXCLUDE_PROFILES = "data/filters/aggregate_liver_12_muscle_candidate_profiles.txt"
DEFAULT_REACTOME_GMT = (
    "src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt"
)
DEFAULT_AGGREGATE_GLARE_TERMS = (
    "outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/"
    "glare_original_style/metascape_results/t_lwm_r7q/top_heatmap_terms.tsv"
)
DEFAULT_STUDIES = ["OSD-379", "OSD-245", "OSD-463"]


def read_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def decode_h5_values(values) -> list[str]:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(str(value))
    return decoded


def infer_stratum(profile: str, accession: str) -> str:
    if accession == "OSD-379":
        match = re.search(r"_LVR_(?:FLT|GC)_([^_]+)_([^_]+)_", profile)
        if match:
            return f"{match.group(1)}_{match.group(2)}"
    if accession == "OSD-245":
        match = re.search(r"_LVR_(?:FLT|GC)_([^_]+)_", profile)
        if match:
            return match.group(1)
    return "all"


def condition_label(condition_inferred: str) -> str:
    if condition_inferred == "flight":
        return "flight"
    if condition_inferred == "ground_control":
        return "ground"
    raise ValueError(f"Unsupported condition: {condition_inferred}")


def extract_raw_count_inputs(
    h5_path: Path,
    liver_metadata_path: Path,
    exclude_profiles_path: Path,
    studies: list[str],
    output_dir: Path,
) -> dict[str, str]:
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    input_dir = output_dir / "raw_deseq2_inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(liver_metadata_path, sep="\t", keep_default_na=False)
    excluded = set(read_lines(exclude_profiles_path))
    selected = metadata[
        metadata["id.accession"].isin(studies)
        & metadata["condition_inferred"].isin(["flight", "ground_control"])
        & ~metadata["profile"].isin(excluded)
    ].copy()
    if selected.empty:
        raise SystemExit("No liver FLT/GC profiles remain after filtering")

    selected["sample"] = selected["profile"].astype(str)
    selected["accession"] = selected["id.accession"].astype(str)
    selected["condition"] = selected["condition_inferred"].map(condition_label)
    selected["stratum"] = [
        infer_stratum(profile, accession)
        for profile, accession in zip(selected["profile"], selected["accession"])
    ]
    selected["source_profile_index"] = selected["source_profile_index"].astype(int)
    selected = selected.sort_values(["accession", "condition", "stratum", "sample"])

    if selected["sample"].duplicated().any():
        duplicates = selected.loc[selected["sample"].duplicated(), "sample"].tolist()
        raise SystemExit(f"Duplicate profile names after filtering: {duplicates[:10]}")

    with h5py.File(h5_path, "r") as handle:
        gene_ids = decode_h5_values(handle["meta/genes/ensembl_gene"][:])
        gene_symbols = decode_h5_values(handle["meta/genes/symbol"][:])
        counts = handle["data/expression"][:, selected["source_profile_index"].tolist()]

    counts_table = pd.DataFrame(counts, index=gene_ids, columns=selected["sample"].tolist())
    counts_table.index.name = "gene_id"
    counts_path = input_dir / "counts.tsv"
    counts_table.to_csv(counts_path, sep="\t")

    metadata_columns = [
        "sample",
        "accession",
        "condition",
        "stratum",
        "profile",
        "source_profile_index",
        "official_sample_name",
        "official_material_type",
        "official_tissue",
        "library_selection",
        "library_layout",
        "sex",
        "strain",
    ]
    available_metadata_columns = [col for col in metadata_columns if col in selected.columns]
    sample_metadata = selected[available_metadata_columns].copy()
    metadata_path = input_dir / "sample_metadata.tsv"
    sample_metadata.to_csv(metadata_path, sep="\t", index=False)

    symbols_path = input_dir / "gene_symbols.tsv"
    pd.DataFrame({"gene_id": gene_ids, "gene_symbol": gene_symbols}).to_csv(
        symbols_path,
        sep="\t",
        index=False,
    )

    counts_by_study = (
        sample_metadata.groupby(["accession", "condition"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    counts_by_study.to_csv(input_dir / "study_condition_counts.tsv", sep="\t", index=False)

    stratum_counts = (
        sample_metadata.groupby(["accession", "stratum", "condition"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    stratum_counts.to_csv(input_dir / "study_stratum_condition_counts.tsv", sep="\t", index=False)

    return {
        "counts": str(counts_path),
        "metadata": str(metadata_path),
        "gene_symbols": str(symbols_path),
        "study_condition_counts": str(input_dir / "study_condition_counts.tsv"),
        "study_stratum_condition_counts": str(input_dir / "study_stratum_condition_counts.tsv"),
    }


def default_rscript() -> str:
    candidate = Path(sys.executable).with_name("Rscript")
    if candidate.exists():
        return str(candidate)
    return shutil.which("Rscript") or "Rscript"


def run_deseq2(
    rscript: str,
    counts: str,
    metadata: str,
    gene_symbols: str,
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
        counts,
        "--metadata",
        metadata,
        "--gene-symbols",
        gene_symbols,
        "--output-dir",
        str(output_dir / "raw_deseq2"),
        "--alpha",
        str(alpha),
        "--min-count",
        str(min_count),
        "--min-samples",
        str(min_samples),
    ]
    subprocess.run(command, check=True)


def read_gmt(path: Path) -> list[dict[str, object]]:
    gene_sets = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term, description, *genes = parts
            gene_sets.append(
                {
                    "term": term,
                    "clean_term": clean_reactome_term(term),
                    "description": description,
                    "genes": {gene for gene in genes if gene},
                }
            )
    return gene_sets


def clean_reactome_term(term: str) -> str:
    if term.startswith("REACTOME_"):
        term = term[len("REACTOME_") :]
    return term.replace("_", " ").title()


def bh_fdr(pvalues) -> list[float]:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return adjusted.tolist()
    valid_values = values[valid]
    order = np.argsort(valid_values)
    ordered = valid_values[order]
    adjusted_ordered = np.empty_like(ordered)
    running_min = 1.0
    n_tests = len(ordered)
    for idx in range(n_tests - 1, -1, -1):
        rank = idx + 1
        running_min = min(running_min, ordered[idx] * n_tests / rank)
        adjusted_ordered[idx] = running_min
    restored = np.empty_like(adjusted_ordered)
    restored[order] = np.minimum(adjusted_ordered, 1.0)
    adjusted[valid] = restored
    return adjusted.tolist()


def run_rank_sum_enrichment(
    deseq,
    gene_sets: list[dict[str, object]],
    studies: list[str],
    min_size: int,
    max_size: int,
):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    mannwhitneyu = require_import(
        "scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt"
    ).mannwhitneyu

    rows = []
    for accession in studies:
        study = deseq[deseq["accession"] == accession].dropna(subset=["stat"]).copy()
        if study.empty:
            continue
        stats = study.groupby("gene_id")["stat"].first()
        universe = set(stats.index)
        all_values = stats.to_numpy(dtype=float)
        all_mean = float(np.mean(all_values))
        for gene_set in gene_sets:
            genes = sorted(gene_set["genes"] & universe)
            set_size = len(genes)
            if set_size < min_size or set_size > max_size:
                continue
            in_values = stats.loc[genes].to_numpy(dtype=float)
            outside_count = len(all_values) - len(in_values)
            if outside_count < min_size:
                continue
            in_mean = float(np.mean(in_values))
            outside_values = stats.loc[~stats.index.isin(genes)].to_numpy(dtype=float)
            for direction, alternative in [
                ("up_in_flight", "greater"),
                ("down_in_flight", "less"),
            ]:
                test = mannwhitneyu(in_values, outside_values, alternative=alternative)
                rows.append(
                    {
                        "accession": accession,
                        "direction": direction,
                        "term": gene_set["term"],
                        "clean_term": gene_set["clean_term"],
                        "description": gene_set["description"],
                        "pathway_genes_tested": set_size,
                        "universe_genes_tested": len(all_values),
                        "mean_wald_stat_in_pathway": in_mean,
                        "mean_wald_stat_all_genes": all_mean,
                        "mean_wald_stat_shift": in_mean - all_mean,
                        "p_value": float(test.pvalue),
                        "overlap_genes": ",".join(genes[:200]),
                    }
                )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr_bh"] = result.groupby(["accession", "direction"])["p_value"].transform(bh_fdr)
    return result.sort_values(
        ["accession", "direction", "fdr_bh", "p_value"],
        ascending=[True, True, True, True],
    )


def run_ora_enrichment(
    deseq,
    gene_sets: list[dict[str, object]],
    studies: list[str],
    alpha: float,
    min_overlap: int,
):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    hypergeom = require_import(
        "scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt"
    ).hypergeom

    rows = []
    for accession in studies:
        study = deseq[deseq["accession"] == accession].dropna(subset=["padj"]).copy()
        universe = set(study["gene_id"].astype(str))
        if not universe:
            continue
        queries = {
            "up_in_flight": set(
                study.loc[
                    (study["padj"] < alpha) & (study["log2FoldChange"] > 0),
                    "gene_id",
                ].astype(str)
            ),
            "down_in_flight": set(
                study.loc[
                    (study["padj"] < alpha) & (study["log2FoldChange"] < 0),
                    "gene_id",
                ].astype(str)
            ),
        }
        for direction, query in queries.items():
            query = query & universe
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
                        "accession": accession,
                        "direction": direction,
                        "term": gene_set["term"],
                        "clean_term": gene_set["clean_term"],
                        "description": gene_set["description"],
                        "query_genes": len(query),
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
    result["fdr_bh"] = result.groupby(["accession", "direction"])["p_value"].transform(bh_fdr)
    return result.sort_values(
        ["accession", "direction", "fdr_bh", "p_value"],
        ascending=[True, True, True, True],
    )


def pathway_categories(text: str) -> set[str]:
    text = text.lower()
    categories = set()
    keyword_groups = {
        "lipid_steroid_xenobiotic_metabolism": [
            "lipid",
            "fatty acid",
            "cholesterol",
            "steroid",
            "xenobiotic",
            "biological oxidation",
            "cytochrome",
            "triglyceride",
            "monocarboxylic",
            "small molecule",
            "sulfur",
        ],
        "translation_ribosome": ["translation", "ribosom", "rrna", "trna"],
        "immune_complement_coagulation": [
            "immune",
            "immunoregulatory",
            "lymphoid",
            "inflammasome",
            "interferon",
            "complement",
            "coagulation",
            "clot",
            "clotting",
            "hemostasis",
            "platelet",
            "cytokine",
            "inflammatory",
        ],
        "mitochondrial_energy": [
            "mitochond",
            "respiratory",
            "oxidative phosphorylation",
            "electron transport",
            "tca",
            "citric acid",
        ],
        "muscle_cytoskeleton": [
            "muscle",
            "contractile",
            "cytoskeleton",
            "actin",
            "myosin",
            "sarcomere",
        ],
        "protein_processing_localization": [
            "protein localization",
            "protein processing",
            "ubiquitin",
            "proteasome",
            "quality control",
            "folding",
        ],
    }
    for category, keywords in keyword_groups.items():
        if any(keyword in text for keyword in keywords):
            categories.add(category)
    return categories


def aggregate_glare_categories(path: Path) -> set[str]:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    if not path.exists():
        return set()
    table = pd.read_csv(path, sep="\t")
    categories = set()
    for _, row in table.head(50).iterrows():
        text = " ".join(str(row.get(col, "")) for col in table.columns)
        categories.update(pathway_categories(text))
    return categories


def summarize_rank_recurrence(rank_enrichment, alpha: float, aggregate_categories: set[str]):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    if rank_enrichment.empty:
        return pd.DataFrame()
    significant = rank_enrichment[rank_enrichment["fdr_bh"] < alpha].copy()
    if significant.empty:
        return pd.DataFrame(
            columns=[
                "direction",
                "term",
                "clean_term",
                "study_count",
                "accessions",
                "best_fdr_bh",
                "mean_wald_stat_shift_mean",
                "aggregate_glare_support_categories",
            ]
        )

    rows = []
    for (direction, term, clean_term), group in significant.groupby(
        ["direction", "term", "clean_term"]
    ):
        categories = pathway_categories(clean_term)
        support_categories = sorted(categories & aggregate_categories)
        rows.append(
            {
                "direction": direction,
                "term": term,
                "clean_term": clean_term,
                "study_count": group["accession"].nunique(),
                "accessions": ",".join(sorted(group["accession"].unique())),
                "best_fdr_bh": group["fdr_bh"].min(),
                "mean_wald_stat_shift_mean": group["mean_wald_stat_shift"].mean(),
                "aggregate_glare_support_categories": ",".join(support_categories),
            }
        )
    recurrence = pd.DataFrame(rows)
    return recurrence.sort_values(
        ["study_count", "best_fdr_bh", "clean_term"],
        ascending=[False, True, True],
    )


def top_gene_table(deseq, studies: list[str], alpha: float):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    for accession in studies:
        study = deseq[deseq["accession"] == accession].copy()
        for direction, sign in [("up_in_flight", 1), ("down_in_flight", -1)]:
            if sign > 0:
                selected = study[study["log2FoldChange"] > 0].copy()
            else:
                selected = study[study["log2FoldChange"] < 0].copy()
            selected = selected.sort_values(["padj", "pvalue"], na_position="last").head(50)
            selected["rank_direction"] = direction
            selected["significant_for_query"] = selected["padj"] < alpha
            rows.append(selected)
    if not rows:
        return pd.DataFrame()
    result = pd.concat(rows, ignore_index=True)
    keep = [
        "accession",
        "rank_direction",
        "gene_id",
        "gene_symbol",
        "baseMean",
        "log2FoldChange",
        "stat",
        "pvalue",
        "padj",
        "significant_for_query",
    ]
    return result[[col for col in keep if col in result.columns]]


def markdown_table(frame, columns: list[str], max_rows: int = 12) -> list[str]:
    if frame.empty:
        return ["No rows."]
    display = frame.loc[:, columns].head(max_rows).copy()
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    lines = [header, separator]
    for _, row in display.iterrows():
        values = []
        for col in columns:
            value = row[col]
            if isinstance(value, float):
                value = f"{value:.3g}"
            values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def write_report(
    output_dir: Path,
    studies: list[str],
    input_paths: dict[str, str],
    deseq_summary,
    rank_enrichment,
    ora_enrichment,
    recurrence,
    alpha: float,
    min_size: int,
    max_size: int,
    aggregate_categories: set[str],
) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    sample_counts = pd.read_csv(input_paths["study_condition_counts"], sep="\t")
    stratum_counts = pd.read_csv(input_paths["study_stratum_condition_counts"], sep="\t")

    recurring = recurrence[recurrence["study_count"] >= 2].copy() if not recurrence.empty else recurrence
    top_rank = (
        rank_enrichment.sort_values(["accession", "direction", "fdr_bh", "p_value"])
        .groupby(["accession", "direction"], as_index=False)
        .head(5)
        if not rank_enrichment.empty
        else rank_enrichment
    )
    top_ora = (
        ora_enrichment.sort_values(["accession", "direction", "fdr_bh", "p_value"])
        .groupby(["accession", "direction"], as_index=False)
        .head(5)
        if not ora_enrichment.empty
        else ora_enrichment
    )

    lines = [
        "# Study-Specific FLT-vs-GC Pathway Recurrence",
        "",
        "This analysis treats OSD-379, OSD-245, and OSD-463 as separate primary studies.",
        "Aggregate GLARE is used only as exploratory support after study-specific raw-count DESeq2 and Reactome pathway testing.",
        "",
        "## Inputs",
        "",
        f"- Raw expression HDF5: `{DEFAULT_H5}`",
        f"- Liver metadata: `{DEFAULT_LIVER_METADATA}`",
        f"- Excluded muscle-outlier profiles: `{DEFAULT_EXCLUDE_PROFILES}`",
        f"- Reactome GMT: `{DEFAULT_REACTOME_GMT}`",
        f"- DESeq2/pathway alpha: {alpha}",
        f"- Rank-sum pathway size range: {min_size}-{max_size} tested genes",
        "",
        "## Sample Counts After 12-Outlier Filter",
        "",
        *markdown_table(sample_counts, list(sample_counts.columns), max_rows=10),
        "",
        "## Study-Internal Strata",
        "",
        *markdown_table(stratum_counts, list(stratum_counts.columns), max_rows=20),
        "",
        "## DESeq2 Summary",
        "",
        *markdown_table(
            deseq_summary,
            [
                "accession",
                "n_flight",
                "n_ground",
                "genes_tested",
                "significant_padj05",
                "significant_up",
                "significant_down",
                "design",
                "dispersion_fit",
            ],
            max_rows=10,
        ),
        "",
        "## Recurring Reactome Pathways",
        "",
    ]

    if recurring.empty:
        lines.extend(
            [
                "No same-direction Reactome pathway passed FDR < "
                f"{alpha} in at least two of the three studies.",
                "That means the strongest evidence remains study-specific rather than recurrent across OSD-379, OSD-245, and OSD-463.",
            ]
        )
    else:
        lines.extend(
            markdown_table(
                recurring,
                [
                    "direction",
                    "clean_term",
                    "study_count",
                    "accessions",
                    "best_fdr_bh",
                    "mean_wald_stat_shift_mean",
                    "aggregate_glare_support_categories",
                ],
                max_rows=20,
            )
        )

    lines.extend(
        [
            "",
            "## Top Per-Study Rank-Sum Pathways",
            "",
            "Reactome rank-sum tests ask whether all genes in a pathway are shifted up or down in the DESeq2 Wald-statistic ranking.",
            "",
            *markdown_table(
                top_rank,
                [
                    "accession",
                    "direction",
                    "clean_term",
                    "pathway_genes_tested",
                    "mean_wald_stat_shift",
                    "p_value",
                    "fdr_bh",
                ],
                max_rows=30,
            ),
            "",
            "## ORA On Significant DESeq2 Genes",
            "",
        ]
    )
    if top_ora.empty:
        lines.append(
            "ORA had no pathway rows, usually because a study/direction had too few FDR-significant DE genes for overlap testing."
        )
    else:
        lines.extend(
            markdown_table(
                top_ora,
                [
                    "accession",
                    "direction",
                    "clean_term",
                    "query_genes",
                    "overlap",
                    "p_value",
                    "fdr_bh",
                ],
                max_rows=30,
            )
        )

    lines.extend(
        [
            "",
            "## Aggregate GLARE Support",
            "",
            "Aggregate GLARE is not used as the main evidence here because prior audits showed strong study/mission/batch structure.",
            "It is only used to flag whether study-specific Reactome terms fall into broad categories also seen in aggregate GLARE/Metascape clusters.",
            f"Aggregate GLARE support categories detected: {', '.join(sorted(aggregate_categories)) or 'none'}",
            "",
            "## Interpretation Rule",
            "",
            "Use a pathway as stronger evidence only when it appears in the same direction in multiple individual studies.",
            "Use one-study pathways as study-specific findings, and use aggregate GLARE clusters only as exploratory context for gene modules or follow-up hypotheses.",
        ]
    )

    report_path = output_dir / "STUDY_SPECIFIC_PATHWAY_RECURRENCE.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", default=DEFAULT_H5)
    parser.add_argument("--liver-metadata", default=DEFAULT_LIVER_METADATA)
    parser.add_argument("--exclude-profiles", default=DEFAULT_EXCLUDE_PROFILES)
    parser.add_argument("--reactome-gmt", default=DEFAULT_REACTOME_GMT)
    parser.add_argument("--aggregate-glare-terms", default=DEFAULT_AGGREGATE_GLARE_TERMS)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--studies", nargs="+", default=DEFAULT_STUDIES)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-count", type=int, default=10)
    parser.add_argument("--min-samples", type=int, default=3)
    parser.add_argument("--min-pathway-size", type=int, default=10)
    parser.add_argument("--max-pathway-size", type=int, default=500)
    parser.add_argument("--min-ora-overlap", type=int, default=3)
    parser.add_argument("--rscript", default=default_rscript())
    parser.add_argument("--skip-deseq2", action="store_true")
    return parser.parse_args()


def main() -> None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    input_paths = extract_raw_count_inputs(
        Path(args.h5),
        Path(args.liver_metadata),
        Path(args.exclude_profiles),
        args.studies,
        output_dir,
    )
    if not args.skip_deseq2:
        run_deseq2(
            args.rscript,
            input_paths["counts"],
            input_paths["metadata"],
            input_paths["gene_symbols"],
            output_dir,
            args.alpha,
            args.min_count,
            args.min_samples,
        )

    deseq_path = output_dir / "raw_deseq2" / "per_study_deseq2.tsv"
    summary_path = output_dir / "raw_deseq2" / "study_deseq2_summary.tsv"
    if not deseq_path.exists() or not summary_path.exists():
        raise SystemExit(
            "Missing DESeq2 results. Re-run without --skip-deseq2 or check raw_deseq2 output."
        )

    deseq = pd.read_csv(deseq_path, sep="\t")
    deseq_summary = pd.read_csv(summary_path, sep="\t")
    gene_sets = read_gmt(Path(args.reactome_gmt))
    aggregate_categories = aggregate_glare_categories(Path(args.aggregate_glare_terms))

    rank_enrichment = run_rank_sum_enrichment(
        deseq,
        gene_sets,
        args.studies,
        args.min_pathway_size,
        args.max_pathway_size,
    )
    rank_path = output_dir / "reactome_rank_sum_pathways.tsv"
    rank_enrichment.to_csv(rank_path, sep="\t", index=False)

    ora_enrichment = run_ora_enrichment(
        deseq,
        gene_sets,
        args.studies,
        args.alpha,
        args.min_ora_overlap,
    )
    ora_path = output_dir / "reactome_ora_significant_deseq2_pathways.tsv"
    ora_enrichment.to_csv(ora_path, sep="\t", index=False)

    recurrence = summarize_rank_recurrence(
        rank_enrichment,
        args.alpha,
        aggregate_categories,
    )
    recurrence_path = output_dir / "recurring_reactome_pathways.tsv"
    recurrence.to_csv(recurrence_path, sep="\t", index=False)

    top_genes = top_gene_table(deseq, args.studies, args.alpha)
    top_genes.to_csv(output_dir / "top_deseq2_genes_by_study_direction.tsv", sep="\t", index=False)

    report_path = write_report(
        output_dir,
        args.studies,
        input_paths,
        deseq_summary,
        rank_enrichment,
        ora_enrichment,
        recurrence,
        args.alpha,
        args.min_pathway_size,
        args.max_pathway_size,
        aggregate_categories,
    )

    summary = {
        "studies": args.studies,
        "output_dir": str(output_dir),
        "inputs": input_paths,
        "deseq2": {
            "per_study": str(deseq_path),
            "summary": str(summary_path),
        },
        "pathways": {
            "rank_sum": str(rank_path),
            "ora": str(ora_path),
            "recurrence": str(recurrence_path),
        },
        "report": str(report_path),
        "aggregate_glare_support_categories": sorted(aggregate_categories),
    }
    (output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
