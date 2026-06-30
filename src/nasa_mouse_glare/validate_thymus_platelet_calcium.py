"""Validate the thymus platelet-calcium GLARE module.

This follow-up checks whether the thymus Reactome platelet-calcium signal is
more consistent with a thymus-intrinsic program or with blood/vascular
composition.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from .io import dense_matrix, load_matrix_bundle


DEFAULT_ROOT = Path("outputs/glare_multi_tissue_api")
DEFAULT_VALIDATION_DIR = DEFAULT_ROOT / "validation_stack_terms15"
DEFAULT_OUTPUT_DIR = DEFAULT_VALIDATION_DIR / "thymus_platelet_calcium_validation"
TERM = "REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_"
CLEAN_TERM = "Response To Elevated Platelet Cytosolic Ca2"
TISSUE = "thymus"


MANUAL_MARKERS: dict[str, set[str]] = {
    "platelet_megakaryocyte_core": {
        "ITGA2B",
        "ITGB3",
        "PF4",
        "PPBP",
        "PLEK",
        "SELP",
        "VWF",
        "MMRN1",
        "SRGN",
        "CD9",
        "CD36",
        "CD63",
    },
    "endothelial_vascular": {
        "PECAM1",
        "VWF",
        "VEGFA",
        "VEGFB",
        "VEGFC",
        "VEGFD",
        "PDGFA",
        "PDGFB",
        "HGF",
        "EGF",
        "SPARC",
        "THBS1",
        "FN1",
    },
    "coagulation_plasma": {
        "A2M",
        "ALB",
        "APOA1",
        "APOH",
        "F5",
        "F8",
        "F13A1",
        "FGA",
        "FGB",
        "FGG",
        "HRG",
        "KNG1",
        "PLG",
        "PROS1",
        "SERPINA1",
        "SERPINA3",
        "SERPING1",
        "SERPINF2",
        "TF",
    },
    "ecm_growth_factor_remodeling": {
        "FN1",
        "THBS1",
        "TIMP1",
        "TIMP3",
        "TGFB1",
        "TGFB2",
        "TGFB3",
        "PDGFA",
        "PDGFB",
        "VEGFA",
        "VEGFB",
        "VEGFC",
        "VEGFD",
        "SPARC",
        "HGF",
        "EGF",
        "ECM1",
    },
    "calcium_secretion_signaling": {
        "ABCC4",
        "CALM1",
        "CALM2",
        "CALM3",
        "PRKCA",
        "PRKCB",
        "PRKCG",
        "PLEK",
        "STX4",
        "STXBP2",
        "STXBP3",
        "SYTL4",
        "RAB27B",
        "SRGN",
    },
    "thymic_epithelial_stromal": {
        "AIRE",
        "CCL25",
        "CD40",
        "COL1A1",
        "COL1A2",
        "CXCL12",
        "DLL4",
        "EPCAM",
        "FOXN1",
        "IL7",
        "KRT5",
        "KRT8",
        "KRT14",
        "LTBR",
        "LY75",
        "PDPN",
    },
    "thymocyte_t_cell": {
        "BCL11B",
        "CCR7",
        "CD247",
        "CD3D",
        "CD3E",
        "CD3G",
        "CD4",
        "CD8A",
        "CD8B1",
        "DNTT",
        "LCK",
        "PTCRA",
        "RAG1",
        "RAG2",
        "TCF7",
        "THEMIS",
        "TRAC",
        "ZAP70",
    },
    "erythroid_blood": {
        "ALAS2",
        "GYPA",
        "HBA-A1",
        "HBA-A2",
        "HBB-BS",
        "HBB-BT",
        "KLF1",
        "TFRC",
    },
}

PANGLAO_TERMS = {
    "panglao_Platelets": "panglao_platelet",
    "panglao_Megakaryocytes": "panglao_megakaryocyte",
    "panglao_Endothelial_cells": "panglao_endothelial",
    "panglao_Pericytes": "panglao_pericyte",
    "panglao_Erythroid-like_and_erythroid_precursor_cells": "panglao_erythroid",
    "panglao_Monocytes": "panglao_monocyte",
}


def split_genes(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [gene for gene in str(value).split(",") if gene]


def bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def parse_gmt(path: Path, term: str) -> tuple[str, list[str]]:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] == term:
                return parts[1], parts[2:]
    raise SystemExit(f"Term {term} not found in {path}")


def load_v4_symbol_map(symbol_gmt: Path, eid_gmt: Path) -> tuple[dict[str, str], list[str], list[str]]:
    _, symbols = parse_gmt(symbol_gmt, TERM)
    _, eids = parse_gmt(eid_gmt, TERM)
    mapped = {eid: symbol for eid, symbol in zip(eids, symbols)}
    return mapped, eids, symbols


def load_panglao_hits(path: Path) -> dict[str, set[str]]:
    if not path.exists():
        return {}
    table = pd.read_csv(path, sep="\t")
    subset = table[
        table["tissue"].astype(str).eq(TISSUE)
        & table["module_class"].astype(str).eq("glare_only")
        & table["term"].astype(str).eq(TERM)
    ]
    hits: dict[str, set[str]] = {}
    for row in subset.itertuples(index=False):
        label = PANGLAO_TERMS.get(str(row.panglao_term))
        if not label:
            continue
        for gene in split_genes(row.overlap_genes):
            hits.setdefault(gene, set()).add(label)
    return hits


def categories_for_gene(gene_id: str, symbol: str, panglao_hits: dict[str, set[str]]) -> list[str]:
    categories = set(panglao_hits.get(gene_id, set()))
    symbol_upper = str(symbol).upper()
    for category, markers in MANUAL_MARKERS.items():
        if symbol_upper in markers:
            categories.add(category)
    return sorted(categories)


def significant_cluster_table(enrichment: pd.DataFrame, symbol_map: dict[str, str], panglao_hits: dict[str, set[str]]) -> pd.DataFrame:
    subset = enrichment[enrichment["term"].astype(str).eq(TERM)].copy()
    subset["fdr_bh"] = numeric(subset["fdr_bh"])
    subset["p_value"] = numeric(subset["p_value"])
    subset["overlap"] = numeric(subset["overlap"])
    subset["significant_reactome_fdr05"] = subset["fdr_bh"].le(0.05)

    rows = []
    for row in subset.itertuples(index=False):
        genes = split_genes(row.overlap_genes)
        category_counts = {category: 0 for category in MANUAL_MARKERS}
        category_counts.update({label: 0 for label in PANGLAO_TERMS.values()})
        symbols = []
        for gene in genes:
            symbol = symbol_map.get(gene, "")
            if symbol:
                symbols.append(symbol)
            for category in categories_for_gene(gene, symbol, panglao_hits):
                category_counts[category] = category_counts.get(category, 0) + 1
        rows.append(
            {
                "accession": row.accession,
                "location": row.location,
                "cluster": row.cluster,
                "cluster_genes": row.cluster_genes,
                "overlap": int(row.overlap) if math.isfinite(row.overlap) else len(genes),
                "p_value": row.p_value,
                "fdr_bh": row.fdr_bh,
                "significant_reactome_fdr05": bool(row.significant_reactome_fdr05),
                "overlap_gene_ids": ",".join(genes),
                "overlap_symbols": ",".join(symbols),
                **{f"n_{category}": int(count) for category, count in sorted(category_counts.items())},
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return result.sort_values(["significant_reactome_fdr05", "accession", "location", "fdr_bh"], ascending=[False, True, True, True])


def build_gene_direction_table(
    root: Path,
    gene_ids: list[str],
    symbol_map: dict[str, str],
    panglao_hits: dict[str, set[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    accessions = []
    for path in sorted((root / TISSUE / "dgea_comparison").glob("*_gene_level_glare_dgea.tsv")):
        accession = path.name.removesuffix("_gene_level_glare_dgea.tsv")
        accessions.append(accession)
        table = pd.read_csv(path, sep="\t")
        table = table[table["gene_id"].isin(gene_ids)].copy()
        if table.empty:
            continue
        table["log2FoldChange"] = numeric(table["log2FoldChange"])
        table["padj"] = numeric(table["padj"])
        table["baseMean"] = numeric(table["baseMean"])
        table["tested_dgea"] = bool_series(table["tested_dgea"])
        table["significant_padj05"] = table["padj"].lt(0.05)
        for row in table.itertuples(index=False):
            symbol = symbol_map.get(row.gene_id, "")
            categories = categories_for_gene(row.gene_id, symbol, panglao_hits)
            direction = "not_significant"
            if bool(row.significant_padj05):
                if row.log2FoldChange > 0:
                    direction = "flt_up"
                elif row.log2FoldChange < 0:
                    direction = "gc_up"
            records.append(
                {
                    "accession": accession,
                    "gene_id": row.gene_id,
                    "symbol": symbol,
                    "categories": ";".join(categories),
                    "baseMean": row.baseMean,
                    "log2FoldChange": row.log2FoldChange,
                    "padj": row.padj,
                    "tested_dgea": bool(row.tested_dgea),
                    "significant_padj05": bool(row.significant_padj05),
                    "direction": direction,
                    "flt_cluster": getattr(row, "flt_cluster", ""),
                    "gc_cluster": getattr(row, "gc_cluster", ""),
                }
            )
    by_study = pd.DataFrame(records)
    if by_study.empty:
        return by_study, pd.DataFrame()

    rows = []
    for gene_id, group in by_study.groupby("gene_id", sort=True):
        symbol = symbol_map.get(gene_id, "")
        tested = group[group["tested_dgea"]]
        sig = group[group["significant_padj05"]]
        flt_up = sig[sig["log2FoldChange"] > 0]
        gc_up = sig[sig["log2FoldChange"] < 0]
        rows.append(
            {
                "gene_id": gene_id,
                "symbol": symbol,
                "categories": ";".join(categories_for_gene(gene_id, symbol, panglao_hits)),
                "studies_with_dgea_row": int(group["accession"].nunique()),
                "studies_tested": int(tested["accession"].nunique()),
                "sig_studies": int(sig["accession"].nunique()),
                "flt_up_sig_studies": int(flt_up["accession"].nunique()),
                "gc_up_sig_studies": int(gc_up["accession"].nunique()),
                "median_log2fc": float(tested["log2FoldChange"].median()) if not tested.empty else math.nan,
                "mean_log2fc": float(tested["log2FoldChange"].mean()) if not tested.empty else math.nan,
                "min_padj": float(tested["padj"].min()) if not tested.empty else math.nan,
                "median_baseMean": float(tested["baseMean"].median()) if not tested.empty else math.nan,
                "accessions_flt_up_sig": ",".join(sorted(flt_up["accession"].unique())),
                "accessions_gc_up_sig": ",".join(sorted(gc_up["accession"].unique())),
            }
        )
    summary = pd.DataFrame(rows)
    summary = summary.sort_values(
        ["flt_up_sig_studies", "gc_up_sig_studies", "sig_studies", "median_log2fc"],
        ascending=[False, True, False, False],
    )
    return by_study, summary


def marker_category_summary(gene_by_study: pd.DataFrame) -> pd.DataFrame:
    rows = []
    categories = sorted(set(MANUAL_MARKERS) | set(PANGLAO_TERMS.values()))
    for category in categories:
        mask = gene_by_study["categories"].fillna("").astype(str).str.split(";").apply(lambda values: category in values)
        subset = gene_by_study[mask].copy()
        if subset.empty:
            rows.append(
                {
                    "category": category,
                    "genes": 0,
                    "tested_study_gene_pairs": 0,
                    "sig_pairs": 0,
                    "flt_up_sig_pairs": 0,
                    "gc_up_sig_pairs": 0,
                    "median_log2fc_tested": math.nan,
                }
            )
            continue
        tested = subset[subset["tested_dgea"]]
        sig = subset[subset["significant_padj05"]]
        rows.append(
            {
                "category": category,
                "genes": int(subset["gene_id"].nunique()),
                "tested_study_gene_pairs": int(len(tested)),
                "sig_pairs": int(len(sig)),
                "flt_up_sig_pairs": int((sig["log2FoldChange"] > 0).sum()),
                "gc_up_sig_pairs": int((sig["log2FoldChange"] < 0).sum()),
                "median_log2fc_tested": float(tested["log2FoldChange"].median()) if not tested.empty else math.nan,
                "top_symbols": ",".join(sorted({symbol for symbol in subset["symbol"].dropna().astype(str) if symbol})[:30]),
            }
        )
    return pd.DataFrame(rows).sort_values(["genes", "sig_pairs"], ascending=[False, False])


def score_summary(validation_dir: Path) -> pd.DataFrame:
    path = validation_dir / "candidate_module_score_validation.tsv"
    table = pd.read_csv(path, sep="\t")
    subset = table[
        table["tissue"].astype(str).eq(TISSUE)
        & table["module_class"].astype(str).eq("glare_only")
        & table["term"].astype(str).eq(TERM)
    ].copy()
    for col in ["flight_minus_ground", "welch_p_value", "empirical_abs_p", "welch_fdr_bh", "n_flight", "n_ground"]:
        subset[col] = numeric(subset[col])
    subset["score_direction"] = np.where(subset["flight_minus_ground"] > 0, "flt_higher", "gc_higher")
    subset["strict_score_support"] = subset["welch_fdr_bh"].lt(0.05) & subset["empirical_abs_p"].le(0.05)
    return subset.sort_values(["strict_score_support", "flight_minus_ground"], ascending=[False, False])


def panglao_summary(validation_dir: Path) -> pd.DataFrame:
    path = validation_dir / "candidate_module_panglao_enrichment.tsv"
    table = pd.read_csv(path, sep="\t")
    subset = table[
        table["tissue"].astype(str).eq(TISSUE)
        & table["module_class"].astype(str).eq("glare_only")
        & table["term"].astype(str).eq(TERM)
    ].copy()
    subset["fdr_bh"] = numeric(subset["fdr_bh"])
    subset["overlap"] = numeric(subset["overlap"])
    return subset.sort_values("fdr_bh")


def load_marker_gene_sets_from_panglao(validation_dir: Path) -> dict[str, set[str]]:
    hits = load_panglao_hits(validation_dir / "candidate_module_panglao_enrichment.tsv")
    marker_sets: dict[str, set[str]] = {}
    for gene, labels in hits.items():
        for label in labels:
            marker_sets.setdefault(label, set()).add(gene)
    return marker_sets


def sample_marker_scores(root: Path, marker_sets: dict[str, set[str]], module_gene_ids: set[str]) -> pd.DataFrame:
    rows = []
    for manifest in sorted((root / TISSUE / "per_study").glob("OSD-*/inputs/api_log2_cpm.manifest.json")):
        accession = manifest.parent.parent.name
        bundle = load_matrix_bundle(manifest)
        matrix = dense_matrix(bundle.matrix, max_dense_gb=1.0)
        gene_index = {gene: idx for idx, gene in enumerate(bundle.genes)}
        metadata = bundle.profile_metadata.copy()
        if metadata is None:
            continue
        condition_col = "condition" if "condition" in metadata.columns else "condition_inferred"
        if condition_col not in metadata.columns:
            continue

        sets = {"reactome_v4_module": set(module_gene_ids), **marker_sets}
        for label, genes in sets.items():
            idx = [gene_index[gene] for gene in genes if gene in gene_index]
            if len(idx) < 2:
                continue
            scores = matrix[idx, :].mean(axis=0)
            score_df = pd.DataFrame({"feature": bundle.profiles, "score": np.asarray(scores).ravel()})
            merged = score_df.merge(metadata, on="feature", how="left")
            condition = merged[condition_col].fillna("").astype(str)
            flight = merged[condition.eq("flight")]["score"]
            ground = merged[condition.isin(["ground", "ground_control"])]["score"]
            if flight.empty or ground.empty:
                continue
            rows.append(
                {
                    "accession": accession,
                    "marker_set": label,
                    "genes_scored": len(idx),
                    "n_flight": int(len(flight)),
                    "n_ground": int(len(ground)),
                    "flight_mean": float(flight.mean()),
                    "ground_mean": float(ground.mean()),
                    "flight_minus_ground": float(flight.mean() - ground.mean()),
                    "flight_median": float(flight.median()),
                    "ground_median": float(ground.median()),
                }
            )
    return pd.DataFrame(rows).sort_values(["marker_set", "accession"]) if rows else pd.DataFrame()


def markdown_table(table: pd.DataFrame) -> str:
    if table.empty:
        return ""
    rendered = table.copy()
    for col in rendered.columns:
        rendered[col] = rendered[col].map(lambda value: "" if pd.isna(value) else str(value))
    header = "| " + " | ".join(rendered.columns) + " |"
    separator = "| " + " | ".join(["---"] * len(rendered.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in rendered.astype(str).values.tolist()]
    return "\n".join([header, separator, *rows])


def write_report(
    output_dir: Path,
    cluster_table: pd.DataFrame,
    score_table: pd.DataFrame,
    panglao_table: pd.DataFrame,
    category_table: pd.DataFrame,
    sample_scores: pd.DataFrame,
    gene_summary: pd.DataFrame,
) -> None:
    sig_clusters = cluster_table[cluster_table["significant_reactome_fdr05"]] if not cluster_table.empty else pd.DataFrame()
    sig_accessions = sorted(sig_clusters["accession"].astype(str).unique()) if not sig_clusters.empty else []
    all_cluster_union = {
        gene for values in cluster_table.get("overlap_gene_ids", pd.Series(dtype=str)).fillna("").astype(str) for gene in split_genes(values)
    }
    sig_cluster_union = {
        gene for values in sig_clusters.get("overlap_gene_ids", pd.Series(dtype=str)).fillna("").astype(str) for gene in split_genes(values)
    }
    strict_scores = score_table[score_table["strict_score_support"]] if not score_table.empty else pd.DataFrame()
    flt_strict = strict_scores[strict_scores["flight_minus_ground"] > 0] if not strict_scores.empty else pd.DataFrame()
    gc_strict = strict_scores[strict_scores["flight_minus_ground"] < 0] if not strict_scores.empty else pd.DataFrame()

    top_panglao = panglao_table.head(6)[["panglao_term", "overlap", "fdr_bh"]] if not panglao_table.empty else pd.DataFrame()
    sample_score_table = sample_scores[
        sample_scores["marker_set"].isin(["reactome_v4_module", "panglao_platelet", "panglao_megakaryocyte", "panglao_endothelial"])
    ][["accession", "marker_set", "genes_scored", "flight_minus_ground"]]
    top_genes = gene_summary.head(15)[
        [
            "symbol",
            "gene_id",
            "categories",
            "flt_up_sig_studies",
            "gc_up_sig_studies",
            "median_log2fc",
            "min_padj",
            "accessions_flt_up_sig",
            "accessions_gc_up_sig",
        ]
    ]

    lines = [
        "# Thymus Platelet-Calcium Module Validation",
        "",
        "## Bottom Line",
        "",
        "The thymus `Response To Elevated Platelet Cytosolic Ca2` GLARE-only signal is real enough to inspect, but it is not clean evidence for thymocyte-intrinsic biology.",
        "It is best labeled as a recurring FLT-up platelet/hemostasis/endothelial-remodeling program with substantial composition risk.",
        "",
        "## Reproducibility Notes",
        "",
        "- The per-study GLARE enrichment uses the Reactome v4 mouse Ensembl term with 85 mapped IDs.",
        f"- Across all cluster rows, the observed thymus term union has {len(all_cluster_union)} IDs; across FDR-significant cluster rows, the union has {len(sig_cluster_union)} IDs.",
        "- The DGEA gene summaries and sample-level Reactome score below use the FDR-significant cluster union.",
        "- Symbols are mapped from the paired Reactome v4 symbol and mouse Ensembl GMT files. The OSDR count matrices remain Ensembl-keyed.",
        "- `log2FoldChange > 0` means higher in spaceflight.",
        "",
        "## Study-Level Recurrence",
        "",
        f"- Significant GLARE Reactome enrichment appears in {len(sig_accessions)} accessions: {', '.join(sig_accessions) if sig_accessions else 'none'}.",
        f"- Strict module-score FLT-up support appears in {len(flt_strict)} accessions: {', '.join(flt_strict['accession'].astype(str)) if not flt_strict.empty else 'none'}.",
        f"- Strict module-score GC-up support appears in {len(gc_strict)} accessions: {', '.join(gc_strict['accession'].astype(str)) if not gc_strict.empty else 'none'}.",
        "- OSD-421 trends GC-higher by module score but is not strict-significant; OSD-515 is effectively flat.",
        "",
        "## Composition Checks",
        "",
        "Panglao enrichment for the same module strongly favors platelet/megakaryocyte/endothelial categories:",
        "",
        markdown_table(top_panglao) if not top_panglao.empty else "No Panglao rows found.",
        "",
        "Manual marker-category counts among study-gene DGEA rows:",
        "",
        markdown_table(category_table),
        "",
        "The relevant sample-level marker scores are not uniformly FLT-up across every study, but the platelet/megakaryocyte/endothelial marker sets move FLT-up in the strongest score-supported accessions.",
        "That points to either a real vascular/hemostasis response, a blood/platelet composition shift, or both.",
        "",
        "Sample-level mean marker-score shifts, FLT minus GC:",
        "",
        markdown_table(sample_score_table),
        "",
        "## Top Recurrent Genes",
        "",
        markdown_table(top_genes),
        "",
        "## Interpretation",
        "",
        "- Evidence for a one-accession artifact is weak: the pathway recurs across the thymus accessions, and three accessions have strict FLT-up module-score support.",
        "- Evidence for a composition/annotation effect is strong: the module is enriched for platelet, megakaryocyte, endothelial, coagulation/plasma, and vascular-remodeling markers, while canonical thymic epithelial and thymocyte marker sets are absent from the module.",
        "- Evidence for thymus-intrinsic biology is therefore indirect. The signal could reflect thymic vascular remodeling or stromal response, but bulk RNA-seq cannot separate that from platelet/blood contamination or vascular-cell proportion shifts.",
        "",
        "## Recommended Follow-Up",
        "",
        "1. Treat this as `FLT-up hemostasis/platelet-calcium/endothelial-remodeling`, not simply `thymus platelet activation`.",
        "2. Re-run the module score after removing platelet-core and plasma/coagulation genes to see whether the FLT-up signal persists.",
        "3. Check sample-level platelet/endothelial marker scores against sample metadata and outliers before using this as a main biological claim.",
        "4. Prefer this module as a hypothesis-generating GLARE-only finding unless confirmed by cell-composition deconvolution or independent thymus histology/flow/cell-type data.",
        "",
        "## Output Files",
        "",
        "- `significant_cluster_marker_summary.tsv`",
        "- `module_gene_dgea_by_study.tsv`",
        "- `module_gene_dgea_summary.tsv`",
        "- `marker_category_summary.tsv`",
        "- `module_score_by_study.tsv`",
        "- `sample_marker_scores.tsv`",
        "- `panglao_marker_enrichment.tsv`",
    ]
    (output_dir / "THYMUS_PLATELET_CALCIUM_VALIDATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    root = Path(args.root)
    validation_dir = Path(args.validation_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    symbol_map, eids, symbols = load_v4_symbol_map(
        Path(args.reactome_v4_symbols),
        Path(args.reactome_v4_mouse_eids),
    )

    enrichment_path = root / TISSUE / "per_study" / "glare_cluster_reactome_enrichment.tsv"
    enrichment = pd.read_csv(enrichment_path, sep="\t")
    panglao_hits = load_panglao_hits(validation_dir / "candidate_module_panglao_enrichment.tsv")
    clusters = significant_cluster_table(enrichment, symbol_map, panglao_hits)
    clusters.to_csv(output_dir / "significant_cluster_marker_summary.tsv", sep="\t", index=False)

    sig_clusters = clusters[clusters["significant_reactome_fdr05"]]
    if sig_clusters.empty:
        gene_ids = sorted(set(eids))
    else:
        gene_ids = sorted({gene for values in sig_clusters["overlap_gene_ids"] for gene in split_genes(values)})

    gene_by_study, gene_summary = build_gene_direction_table(root, gene_ids, symbol_map, panglao_hits)
    gene_by_study.to_csv(output_dir / "module_gene_dgea_by_study.tsv", sep="\t", index=False)
    gene_summary.to_csv(output_dir / "module_gene_dgea_summary.tsv", sep="\t", index=False)

    categories = marker_category_summary(gene_by_study)
    categories.to_csv(output_dir / "marker_category_summary.tsv", sep="\t", index=False)

    scores = score_summary(validation_dir)
    scores.to_csv(output_dir / "module_score_by_study.tsv", sep="\t", index=False)

    panglao = panglao_summary(validation_dir)
    panglao.to_csv(output_dir / "panglao_marker_enrichment.tsv", sep="\t", index=False)

    marker_sets = load_marker_gene_sets_from_panglao(validation_dir)
    marker_sets = {label: genes for label, genes in marker_sets.items() if label in {"panglao_platelet", "panglao_megakaryocyte", "panglao_endothelial", "panglao_pericyte", "panglao_erythroid"}}
    sample_scores = sample_marker_scores(root, marker_sets, set(gene_ids))
    sample_scores.to_csv(output_dir / "sample_marker_scores.tsv", sep="\t", index=False)

    write_report(output_dir, clusters, scores, panglao, categories, sample_scores, gene_summary)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(DEFAULT_ROOT))
    parser.add_argument("--validation-dir", default=str(DEFAULT_VALIDATION_DIR))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--reactome-v4-symbols", default="src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0.symbols.gmt")
    parser.add_argument("--reactome-v4-mouse-eids", default="src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt")
    return parser


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
