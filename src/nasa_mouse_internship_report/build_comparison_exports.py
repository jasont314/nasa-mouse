"""Build comparison-ready expiMap and generative feature tables."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import zipfile

import pandas as pd
from openpyxl.styles import Font, PatternFill


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "outputs/comparison/selected_features"
EXPIMAP_SOURCE = ROOT / "paper/asgsr_expimap_hvg/source_data"
GENERATIVE_SOURCE = ROOT / "paper/synthetic_guided_spaceflight/source_data"
IMPORTANCE_SOURCE = (
    ROOT
    / "outputs/generative/benchmark/analyses/"
    "classifier_importance_osdr_disjoint_v1"
)
GENE_MAP = ROOT / "data/reference/gencode_vM39_mouse_gene_symbols.tsv.gz"
ARM_CHOICE_SOURCES = {
    "tissue": (
        ROOT
        / "outputs/generative/benchmark/analyses/"
        "within_study_generated_feature_stability_osdr_disjoint_v1/"
        "tissue_arm_choices.tsv"
    ),
    "muscle_group": (
        ROOT
        / "outputs/generative/benchmark/analyses/"
        "within_study_generated_feature_stability_muscle_groups_osdr_disjoint_v1/"
        "tissue_arm_choices.tsv"
    ),
}


def _read(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", low_memory=False)


def _write_tsv(frame: pd.DataFrame, path: Path) -> None:
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        compression=compression,
        lineterminator="\n",
        na_rep="NA",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalize_xlsx(path: Path) -> None:
    """Remove ZIP timestamps so rebuilding the workbook is deterministic."""

    normalized = path.with_suffix(".normalized.xlsx")
    with (
        zipfile.ZipFile(path) as source,
        zipfile.ZipFile(
            normalized,
            "w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=9,
        ) as target,
    ):
        for name in sorted(source.namelist()):
            original = source.getinfo(name)
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = original.create_system
            info.external_attr = original.external_attr
            contents = source.read(name)
            if name == "docProps/core.xml":
                contents = re.sub(
                    rb"(<dcterms:modified[^>]*>)[^<]+(</dcterms:modified>)",
                    rb"\g<1>2000-01-01T00:00:00Z\g<2>",
                    contents,
                )
            target.writestr(info, contents)
    normalized.replace(path)


def _symbol_map() -> pd.Series:
    table = _read(GENE_MAP)
    if table["gene_id"].duplicated().any():
        raise ValueError("GENCODE symbol map contains duplicate gene IDs")
    return table.set_index("gene_id")["gene_symbol"]


def _fill_symbols(frame: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    frame = frame.copy()
    mapped = frame["gene_id"].map(symbols)
    current = frame.get("gene_symbol", pd.Series("", index=frame.index)).astype(str)
    missing = current.eq("") | current.str.startswith("ENSMUSG") | current.eq("nan")
    frame["gene_symbol"] = current.where(~missing, mapped).fillna(frame["gene_id"])
    return frame


def _reactome_ids(terms: pd.Series) -> pd.Series:
    identifiers = terms.astype(str).str.extract(r"^(R-MMU-\d+)", expand=False)
    if identifiers.isna().any():
        unknown = sorted(terms.loc[identifiers.isna()].astype(str).unique())
        raise ValueError(f"Could not parse Reactome IDs from terms: {unknown}")
    return identifiers


def build_expimap_tables() -> dict[str, pd.DataFrame]:
    pathways = _read(EXPIMAP_SOURCE / "table_2_retained_pathway_evidence.tsv")
    pathways.insert(0, "method", "expiMap")
    pathways.insert(1, "feature_type", "retained_pathway")
    pathways.insert(5, "reactome_id", _reactome_ids(pathways["term"]))
    pathways.insert(
        6,
        "pathway_url",
        pathways["reactome_id"].map(
            lambda value: f"https://reactome.org/PathwayBrowser/#/{value}"
        ),
    )
    pathways.insert(
        7,
        "flt_gc_direction",
        pathways["seed_effect_median"].map(
            lambda value: "FLT_higher" if float(value) > 0 else "FLT_lower"
        ),
    )
    pathways = pathways.sort_values(
        ["analysis_role", "tissue", "seed_effect_median"],
        kind="stable",
    ).reset_index(drop=True)

    members = _read(
        EXPIMAP_SOURCE / "table_s35_retained_pathway_member_gene_effects.tsv.gz"
    ).rename(columns={"gene_id": "gene_id", "gene_symbol": "gene_symbol"})
    evidence = pathways[
        [
            "tissue",
            "term",
            "analysis_role",
            "evidence_role",
            "seed_effect_median",
            "robustness_support_count",
            "robustness_status",
            "gsea_fdr",
        ]
    ]
    members = members.merge(
        evidence,
        on=["tissue", "term"],
        how="left",
        validate="many_to_one",
    )
    members = _fill_symbols(members, _symbol_map())
    same = members["same_direction_as_pathway_score"].astype(bool)
    significant = members["gene_fdr"].lt(0.05)
    members.insert(0, "method", "expiMap")
    members.insert(1, "feature_type", "retained_pathway_member_gene")
    members["independently_selected_gene"] = False
    members["member_support_class"] = "discordant_or_nonsignificant"
    members.loc[same, "member_support_class"] = "concordant"
    members.loc[~same & significant, "member_support_class"] = "discordant_bh_fdr"
    members.loc[same & significant, "member_support_class"] = "concordant_bh_fdr"
    members["gene_flt_gc_direction"] = members[
        "project_balanced_gene_log2cpm_effect"
    ].map(lambda value: "FLT_higher" if float(value) > 0 else "FLT_lower")
    members = members.sort_values(
        [
            "analysis_role",
            "tissue",
            "term",
            "same_direction_as_pathway_score",
            "gene_fdr",
            "absolute_gene_effect",
        ],
        ascending=[True, True, True, False, True, False],
        kind="stable",
    ).reset_index(drop=True)

    rows: list[dict[str, object]] = []
    for (tissue, gene_id, gene_symbol), frame in members.groupby(
        ["tissue", "gene_id", "gene_symbol"],
        sort=True,
        observed=True,
    ):
        concordant = frame["same_direction_as_pathway_score"].astype(bool)
        gene_effect = float(frame["project_balanced_gene_log2cpm_effect"].iloc[0])
        rows.append(
            {
                "method": "expiMap",
                "analysis_scope": "canonical_tissue",
                "tissue": tissue,
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "comparison_role": "member_of_retained_expimap_pathway",
                "independently_selected_gene": False,
                "n_retained_pathways": int(frame["term"].nunique()),
                "retained_pathway_ids": ";".join(
                    sorted(_reactome_ids(frame["term"]).unique())
                ),
                "retained_pathway_terms": ";".join(
                    sorted(frame["term"].unique())
                ),
                "retained_pathway_labels": ";".join(
                    sorted(frame["display_label"].unique())
                ),
                "n_concordant_pathways": int(concordant.sum()),
                "n_discordant_pathways": int((~concordant).sum()),
                "any_concordant_bh_fdr": bool(
                    (concordant & frame["gene_fdr"].lt(0.05)).any()
                ),
                "minimum_gene_fdr": float(frame["gene_fdr"].min()),
                "project_balanced_gene_log2cpm_effect": gene_effect,
                "gene_flt_gc_direction": (
                    "FLT_higher" if gene_effect > 0 else "FLT_lower"
                ),
                "analysis_roles": ";".join(sorted(frame["analysis_role"].unique())),
                "pathway_robustness_statuses": ";".join(
                    sorted(frame["robustness_status"].unique())
                ),
            }
        )
    genes = pd.DataFrame(rows).sort_values(
        ["analysis_roles", "tissue", "any_concordant_bh_fdr", "minimum_gene_fdr"],
        ascending=[True, True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    return {
        "expimap_pathways": pathways,
        "expimap_members": members,
        "expimap_genes": genes,
    }


def _candidate_flags() -> tuple[pd.DataFrame, pd.DataFrame]:
    matched = _read(
        GENERATIVE_SOURCE / "table_s22_matched_gene_literature_annotations.tsv"
    ).sort_values(["analysis_scope", "tissue", "meta_fdr", "symbol"], kind="stable")
    consensus = _read(
        GENERATIVE_SOURCE / "table_s16_promoted_gene_literature_annotations.tsv"
    ).sort_values(
        ["analysis_scope", "tissue", "real_meta_fdr", "symbol"],
        kind="stable",
    )
    return matched.reset_index(drop=True), consensus.reset_index(drop=True)


def build_generative_tables() -> dict[str, pd.DataFrame]:
    all_comparison = _read(IMPORTANCE_SOURCE / "arm_vs_real_gene_comparison.tsv.gz")
    all_comparison["analysis_scope"] = all_comparison["scope"].map(
        {"tissue": "canonical_tissue", "muscle_group": "skeletal_muscle_group"}
    )
    if all_comparison["analysis_scope"].isna().any():
        raise ValueError("Unknown all-arm generative feature-importance scope")
    all_comparison.insert(0, "method", "conditional_DDIM_classifier")
    all_comparison.insert(1, "feature_type", "all_arm_stable_gene_feature")
    all_stable = all_comparison.loc[
        all_comparison["real_stable"].astype(bool)
        | all_comparison["arm_stable"].astype(bool)
    ].copy()
    all_stable = all_stable.sort_values(
        ["analysis_scope", "tissue", "arm", "symbol"],
        kind="stable",
    ).reset_index(drop=True)

    comparison = _read(
        IMPORTANCE_SOURCE / "selected_arm_vs_real_gene_comparison.tsv.gz"
    )
    importance = _read(IMPORTANCE_SOURCE / "importance_summary.tsv.gz")
    comparison["analysis_scope"] = comparison["scope"].map(
        {"tissue": "canonical_tissue", "muscle_group": "skeletal_muscle_group"}
    )
    if comparison["analysis_scope"].isna().any():
        raise ValueError("Unknown generative feature-importance scope")

    real_domain = importance.loc[importance["domain"].eq("real")].copy()
    baseline = real_domain.loc[real_domain["arm"].eq("real_only")].copy()
    baseline = baseline[
        [
            "scope",
            "tissue",
            "gene",
            "linear_shap_mean_absolute",
            "linear_shap_flight_minus_ground",
        ]
    ].rename(
        columns={
            "linear_shap_mean_absolute": "real_linear_shap_mean_absolute",
            "linear_shap_flight_minus_ground": (
                "real_linear_shap_flight_minus_ground"
            ),
        }
    )
    selected = real_domain.loc[~real_domain["arm"].eq("real_only")].copy()
    selected = selected[
        [
            "scope",
            "tissue",
            "arm",
            "gene",
            "linear_shap_mean_absolute",
            "linear_shap_flight_minus_ground",
        ]
    ].rename(
        columns={
            "linear_shap_mean_absolute": "selected_linear_shap_mean_absolute",
            "linear_shap_flight_minus_ground": (
                "selected_linear_shap_flight_minus_ground"
            ),
        }
    )
    comparison = comparison.merge(
        baseline,
        on=["scope", "tissue", "gene"],
        how="left",
        validate="one_to_one",
    ).merge(
        selected,
        on=["scope", "tissue", "arm", "gene"],
        how="left",
        validate="one_to_one",
    )
    comparison["selected_arm_importance_rank"] = comparison.groupby(
        ["scope", "tissue"], observed=True
    )["arm_real_permutation_roc_auc"].rank(
        method="min", ascending=False
    )
    comparison["real_only_importance_rank"] = comparison.groupby(
        ["scope", "tissue"], observed=True
    )["real_permutation_roc_auc"].rank(method="min", ascending=False)

    matched, consensus = _candidate_flags()
    candidate_keys = ["analysis_scope", "tissue", "gene"]
    matched_columns = candidate_keys + [
        "matched_statuses",
        "literature_classification",
        "interpretation",
    ]
    matched_flags = matched[matched_columns].rename(
        columns={
            "literature_classification": "matched_literature_classification",
            "interpretation": "matched_interpretation",
        }
    )
    matched_flags["matched_primary"] = True
    consensus_flags = consensus[
        candidate_keys
        + [
            "selection_status",
            "literature_classification",
            "interpretation",
        ]
    ].rename(
        columns={
            "selection_status": "consensus_selection_status",
            "literature_classification": "consensus_literature_classification",
            "interpretation": "consensus_interpretation",
        }
    )
    consensus_flags["consensus_secondary"] = True
    comparison = comparison.merge(
        matched_flags,
        on=candidate_keys,
        how="left",
        validate="many_to_one",
    ).merge(
        consensus_flags,
        on=candidate_keys,
        how="left",
        validate="many_to_one",
    )
    comparison["matched_primary"] = comparison["matched_primary"].fillna(False)
    comparison["consensus_secondary"] = comparison[
        "consensus_secondary"
    ].fillna(False)

    real_associations = _read(
        GENERATIVE_SOURCE / "table_s11_all_random_effects_bh_fdr_genes.tsv"
    )[
        [
            "analysis_scope",
            "tissue",
            "gene",
            "flt_gc_direction",
            "meta_effect",
            "meta_fdr",
            "n_accessions",
        ]
    ]
    comparison = comparison.merge(
        real_associations,
        on=candidate_keys,
        how="left",
        validate="many_to_one",
    )
    comparison.insert(0, "method", "conditional_DDIM_classifier")
    comparison.insert(1, "feature_type", "selected_gene_feature")
    comparison = comparison.sort_values(
        ["analysis_scope", "tissue", "selected_arm_importance_rank", "symbol"],
        kind="stable",
    ).reset_index(drop=True)

    stable = comparison.loc[
        comparison["real_stable"].astype(bool) | comparison["arm_stable"].astype(bool)
    ].copy()
    grouped = _read(
        GENERATIVE_SOURCE / "table_s23_grouped_pathway_literature_annotations.tsv"
    ).sort_values(["scope", "tissue", "meta_fdr", "term"], kind="stable")
    return {
        "generative_all_arm_stable": all_stable,
        "generative_full": comparison,
        "generative_stable": stable.reset_index(drop=True),
        "generative_matched": matched,
        "generative_consensus": consensus,
        "generative_grouped": grouped.reset_index(drop=True),
    }


def build_generative_coverage(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Summarize feature-result availability for all 27 benchmark units."""

    choice_frames = []
    for scope, path in ARM_CHOICE_SOURCES.items():
        frame = _read(path)
        frame.insert(0, "scope", scope)
        choice_frames.append(frame)
    coverage = pd.concat(choice_frames, ignore_index=True)
    coverage["analysis_scope"] = coverage["scope"].map(
        {"tissue": "canonical_tissue", "muscle_group": "skeletal_muscle_group"}
    )
    coverage["selected_arm_feature_comparison_available"] = ~coverage[
        "selected_arm"
    ].eq("real_only")
    coverage["feature_comparison_status"] = coverage[
        "selected_arm_feature_comparison_available"
    ].map(
        {
            True: "selected synthetic-supported arm compared with real-only",
            False: "real-only retained; no selected synthetic-arm comparison",
        }
    )

    keys = ["analysis_scope", "tissue"]
    all_arm = tables["generative_all_arm_stable"]
    all_arm_rows = (
        all_arm.groupby(keys, observed=True)
        .size()
        .rename("all_arm_stable_feature_row_count")
    )
    all_arm_genes = (
        all_arm.groupby(keys, observed=True)["gene"]
        .nunique()
        .rename("all_arm_stable_unique_gene_count")
    )
    coverage = coverage.merge(
        all_arm_rows,
        on=keys,
        how="left",
        validate="one_to_one",
    ).merge(
        all_arm_genes,
        on=keys,
        how="left",
        validate="one_to_one",
    )
    count_sources = {
        "full_selected_arm_feature_count": tables["generative_full"],
        "stable_selected_feature_count": tables["generative_stable"],
        "matched_primary_gene_count": tables["generative_matched"],
        "consensus_secondary_gene_count": tables["generative_consensus"],
    }
    for column, frame in count_sources.items():
        counts = frame.groupby(keys, observed=True).size().rename(column)
        coverage = coverage.merge(counts, on=keys, how="left", validate="one_to_one")

    grouped = tables["generative_grouped"].copy()
    grouped["analysis_scope"] = grouped["scope"].map(
        {"tissue": "canonical_tissue", "muscle_group": "skeletal_muscle_group"}
    )
    grouped_counts = (
        grouped.groupby(keys, observed=True)
        .size()
        .rename("grouped_pathway_count")
    )
    coverage = coverage.merge(
        grouped_counts,
        on=keys,
        how="left",
        validate="one_to_one",
    )
    count_columns = [
        "all_arm_stable_feature_row_count",
        "all_arm_stable_unique_gene_count",
        *count_sources,
        "grouped_pathway_count",
    ]
    coverage[count_columns] = coverage[count_columns].fillna(0).astype(int)
    columns = [
        "analysis_scope",
        "scope",
        "tissue",
        "selected_arm",
        "generated_arm_eligible_all_metrics",
        "selected_arm_feature_comparison_available",
        "feature_comparison_status",
        *count_columns,
    ]
    return coverage[columns].sort_values(keys, kind="stable").reset_index(drop=True)


def build_pathway_crosswalk(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    expimap = tables["expimap_pathways"]
    expimap_rows = pd.DataFrame(
        {
            "method": "expiMap",
            "analysis": "retained_pathway",
            "evidence_tier": expimap["analysis_role"],
            "tissue": expimap["tissue"],
            "pathway_id": expimap["reactome_id"],
            "pathway_term": expimap["term"],
            "pathway_name": expimap["display_label"],
            "pathway_url": expimap["pathway_url"],
            "flt_gc_direction": expimap["flt_gc_direction"],
            "effect": expimap["seed_effect_median"],
            "fdr": expimap["gsea_fdr"],
            "selection_status": expimap["robustness_status"],
            "literature_classification": expimap["evidence_role"],
        }
    )
    grouped = tables["generative_grouped"]
    grouped_ids = _reactome_ids(grouped["term"])
    grouped_rows = pd.DataFrame(
        {
            "method": "conditional_DDIM_classifier",
            "analysis": "grouped_permutation_and_shap",
            "evidence_tier": "primary_grouped_importance",
            "tissue": grouped["tissue"],
            "pathway_id": grouped_ids,
            "pathway_term": grouped["term"],
            "pathway_name": grouped["description"],
            "pathway_url": grouped["url"],
            "flt_gc_direction": grouped["flt_gc_direction"],
            "effect": grouped["meta_effect"],
            "fdr": grouped["meta_fdr"],
            "selection_status": grouped["group_importance_patterns"],
            "literature_classification": grouped["literature_classification"],
        }
    )
    return pd.concat([expimap_rows, grouped_rows], ignore_index=True).sort_values(
        ["method", "evidence_tier", "tissue", "fdr", "pathway_id"],
        kind="stable",
    ).reset_index(drop=True)


def build_gene_crosswalk(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    expimap = tables["expimap_genes"].copy()
    expimap = expimap.rename(
        columns={
            "gene_symbol": "expimap_gene_symbol",
            "project_balanced_gene_log2cpm_effect": "expimap_gene_effect",
            "gene_flt_gc_direction": "expimap_flt_gc_direction",
        }
    )
    expimap["expimap_pathway_member"] = True
    expimap_columns = [
        "analysis_scope",
        "tissue",
        "gene_id",
        "expimap_gene_symbol",
        "expimap_pathway_member",
        "n_retained_pathways",
        "retained_pathway_ids",
        "retained_pathway_terms",
        "retained_pathway_labels",
        "n_concordant_pathways",
        "any_concordant_bh_fdr",
        "minimum_gene_fdr",
        "expimap_gene_effect",
        "expimap_flt_gc_direction",
        "analysis_roles",
    ]
    expimap = expimap[expimap_columns]

    all_arm = tables["generative_all_arm_stable"].copy()
    all_arm["stable_synthetic_arm"] = all_arm["arm"].where(
        all_arm["arm_stable"].astype(bool)
    )
    all_arm = (
        all_arm.groupby(["analysis_scope", "tissue", "gene"], observed=True)
        .agg(
            any_arm_gene_symbol=("symbol", "first"),
            generative_any_arm_real_stable=("real_stable", "max"),
            generative_stable_synthetic_arms=(
                "stable_synthetic_arm",
                lambda values: ";".join(sorted(set(values.dropna()))),
            ),
            generative_all_compared_arms=(
                "arm",
                lambda values: ";".join(sorted(set(values))),
            ),
            generative_all_arm_patterns=(
                "pattern",
                lambda values: ";".join(sorted(set(values))),
            ),
            generative_max_arm_real_permutation_roc_auc=(
                "arm_real_permutation_roc_auc",
                "max",
            ),
            generative_max_arm_minus_real_permutation_roc_auc=(
                "arm_minus_real_permutation_roc_auc",
                "max",
            ),
        )
        .reset_index()
        .rename(columns={"gene": "gene_id"})
    )
    all_arm["generative_any_arm_stable_feature"] = True

    generative = tables["generative_stable"].copy().rename(
        columns={
            "gene": "gene_id",
            "symbol": "generative_gene_symbol",
            "arm": "generative_selected_arm",
            "pattern": "generative_selection_pattern",
            "meta_effect": "generative_real_meta_effect",
            "meta_fdr": "generative_real_meta_fdr",
            "flt_gc_direction": "generative_flt_gc_direction",
        }
    )
    generative["generative_selected_arm_stable_feature"] = True
    generative_columns = [
        "analysis_scope",
        "tissue",
        "gene_id",
        "generative_gene_symbol",
        "generative_selected_arm_stable_feature",
        "generative_selected_arm",
        "generative_selection_pattern",
        "real_stable",
        "arm_stable",
        "real_selection_frequency",
        "arm_selection_frequency",
        "real_permutation_roc_auc",
        "arm_real_permutation_roc_auc",
        "arm_minus_real_permutation_roc_auc",
        "selected_arm_importance_rank",
        "real_linear_shap_flight_minus_ground",
        "selected_linear_shap_flight_minus_ground",
        "matched_primary",
        "consensus_secondary",
        "generative_real_meta_effect",
        "generative_real_meta_fdr",
        "generative_flt_gc_direction",
    ]
    generative = generative[generative_columns]

    matched = tables["generative_matched"][
        ["analysis_scope", "tissue", "gene", "symbol"]
    ].rename(columns={"gene": "gene_id", "symbol": "matched_gene_symbol"})
    matched["matched_primary_candidate"] = True
    consensus = tables["generative_consensus"][
        ["analysis_scope", "tissue", "gene", "symbol"]
    ].rename(columns={"gene": "gene_id", "symbol": "consensus_gene_symbol"})
    consensus["consensus_secondary_candidate"] = True

    keys = ["analysis_scope", "tissue", "gene_id"]
    all_keys = pd.concat(
        [
            expimap[keys],
            all_arm[keys],
            generative[keys],
            matched[keys],
            consensus[keys],
        ],
        ignore_index=True,
    ).drop_duplicates()
    crosswalk = all_keys.merge(expimap, on=keys, how="left", validate="one_to_one")
    crosswalk = crosswalk.merge(
        all_arm,
        on=keys,
        how="left",
        validate="one_to_one",
    ).merge(
        generative,
        on=keys,
        how="left",
        validate="one_to_one",
        suffixes=("", "_stable"),
    ).merge(
        matched,
        on=keys,
        how="left",
        validate="one_to_one",
        suffixes=("", "_matched"),
    ).merge(
        consensus,
        on=keys,
        how="left",
        validate="one_to_one",
        suffixes=("", "_consensus"),
    )
    for column in (
        "expimap_pathway_member",
        "generative_any_arm_stable_feature",
        "generative_selected_arm_stable_feature",
        "matched_primary",
        "matched_primary_candidate",
        "consensus_secondary",
        "consensus_secondary_candidate",
    ):
        if column in crosswalk:
            crosswalk[column] = crosswalk[column].fillna(False).astype(bool)
    crosswalk["matched_primary"] = crosswalk["matched_primary"] | crosswalk[
        "matched_primary_candidate"
    ]
    crosswalk["consensus_secondary"] = crosswalk[
        "consensus_secondary"
    ] | crosswalk["consensus_secondary_candidate"]
    crosswalk = crosswalk.drop(
        columns=[
            column
            for column in (
                "matched_primary_candidate",
                "consensus_secondary_candidate",
            )
            if column in crosswalk
        ]
    )
    symbol_columns = [
        column
        for column in (
            "generative_gene_symbol",
            "any_arm_gene_symbol",
            "matched_gene_symbol",
            "consensus_gene_symbol",
            "expimap_gene_symbol",
        )
        if column in crosswalk
    ]
    crosswalk.insert(
        3,
        "gene_symbol",
        crosswalk[symbol_columns].bfill(axis=1).iloc[:, 0],
    )
    crosswalk["gene_symbol"] = crosswalk["gene_symbol"].fillna(crosswalk["gene_id"])
    crosswalk = crosswalk.drop(columns=symbol_columns)
    return crosswalk.sort_values(
        [
            "analysis_scope",
            "tissue",
            "matched_primary",
            "consensus_secondary",
            "gene_symbol",
        ],
        ascending=[True, True, False, False, True],
        kind="stable",
    ).reset_index(drop=True)


def _write_workbook(tables: dict[str, pd.DataFrame], path: Path) -> None:
    sheet_map = {
        "gene_crosswalk": "gene_crosswalk",
        "pathway_crosswalk": "pathway_crosswalk",
        "generative_coverage": "gen_analysis_coverage",
        "expimap_pathways": "expimap_pathways",
        "expimap_genes": "expimap_genes",
        "expimap_members": "expimap_members",
        "generative_all_arm_stable": "gen_all_arm_stable",
        "generative_stable": "gen_stable_features",
        "generative_matched": "gen_matched_genes",
        "generative_consensus": "gen_consensus_genes",
        "generative_grouped": "gen_grouped_pathways",
    }
    guide = pd.DataFrame(
        [
            {
                "sheet": sheet,
                "rows": len(tables[key]),
                "purpose": {
                    "gene_crosswalk": (
                        "Start here for gene-list comparisons across methods."
                    ),
                    "pathway_crosswalk": "Start here for selected pathway comparisons.",
                    "expimap_pathways": "All 16 retained expiMap pathway records.",
                    "expimap_genes": (
                        "One row per tissue and retained-pathway member gene."
                    ),
                    "expimap_members": (
                        "Every retained pathway-to-gene membership and gene effect."
                    ),
                    "generative_stable": (
                        "Stable real-only or selected-arm classifier features."
                    ),
                    "generative_all_arm_stable": (
                        "Stable features for every synthetic arm in all 27 units."
                    ),
                    "generative_matched": (
                        "Primary matched synthetic-supported gene results."
                    ),
                    "generative_consensus": (
                        "Secondary promoted or reinforced consensus genes."
                    ),
                    "generative_grouped": (
                        "Selected grouped Reactome importance results."
                    ),
                    "generative_coverage": (
                        "Availability and result counts for all 27 analysis units."
                    ),
                }[key],
                "selection_definition": {
                    "gene_crosswalk": (
                        "Filter the five boolean method and evidence flags."
                    ),
                    "pathway_crosswalk": (
                        "Join on canonical pathway_id; pathway_term keeps labels."
                    ),
                    "generative_coverage": (
                        "A zero count means no result in that named analysis."
                    ),
                    "expimap_pathways": (
                        "Pathways retained by the final tissue evidence workflow."
                    ),
                    "expimap_genes": (
                        "Members of retained pathways, not a selected gene panel."
                    ),
                    "expimap_members": (
                        "One row for each retained pathway and member-gene pair."
                    ),
                    "generative_all_arm_stable": (
                        "At least 50% selection and 75% sign agreement in real or arm."
                    ),
                    "generative_stable": (
                        "Same stability rule, restricted to chosen synthetic arms."
                    ),
                    "generative_matched": (
                        "Real BH-FDR association plus synthetic-supported importance."
                    ),
                    "generative_consensus": (
                        "Secondary promoted or reinforced compact-panel result."
                    ),
                    "generative_grouped": (
                        "Reactome group importance supported by permutation and SHAP."
                    ),
                }[key],
            }
            for key, sheet in sheet_map.items()
        ]
    )
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        guide.to_excel(writer, sheet_name="guide", index=False)
        for key, sheet in sheet_map.items():
            tables[key].to_excel(writer, sheet_name=sheet, index=False)
        workbook = writer.book
        workbook.properties.creator = "Jason Trinh"
        workbook.properties.title = "NASA mouse selected feature comparison"
        workbook.properties.created = datetime(2000, 1, 1)
        workbook.properties.modified = datetime(2000, 1, 1)
        for worksheet in workbook.worksheets:
            worksheet.freeze_panes = "A2"
            worksheet.auto_filter.ref = worksheet.dimensions
            for cell in worksheet[1]:
                cell.font = Font(bold=True, color="FFFFFF")
                cell.fill = PatternFill("solid", fgColor="1F4E78")
            for column in worksheet.columns:
                values = [str(cell.value or "") for cell in column[:200]]
                width = min(max(max(map(len, values), default=0) + 2, 10), 42)
                worksheet.column_dimensions[column[0].column_letter].width = width
    _normalize_xlsx(path)


def run(output_dir: Path = OUTPUT_DIR) -> dict[str, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = build_expimap_tables()
    tables.update(build_generative_tables())
    tables["generative_coverage"] = build_generative_coverage(tables)
    tables["gene_crosswalk"] = build_gene_crosswalk(tables)
    tables["pathway_crosswalk"] = build_pathway_crosswalk(tables)

    outputs = {
        "gene_crosswalk": "gene_crosswalk.tsv",
        "pathway_crosswalk": "pathway_crosswalk.tsv",
        "expimap_pathways": "expimap_retained_pathways.tsv",
        "expimap_genes": "expimap_retained_pathway_gene_summary.tsv",
        "expimap_members": "expimap_retained_pathway_members.tsv.gz",
        "generative_full": "generative_full_selected_feature_comparison.tsv.gz",
        "generative_stable": "generative_selected_arm_stable_features.tsv",
        "generative_all_arm_stable": "generative_all_arm_stable_features.tsv.gz",
        "generative_matched": "generative_matched_genes.tsv",
        "generative_consensus": "generative_consensus_genes.tsv",
        "generative_grouped": "generative_grouped_pathways.tsv",
        "generative_coverage": "generative_analysis_coverage.tsv",
    }
    for key, filename in outputs.items():
        _write_tsv(tables[key], output_dir / filename)
    workbook = output_dir / "selected_feature_comparison.xlsx"
    _write_workbook(tables, workbook)

    source_paths = [
        EXPIMAP_SOURCE / "table_2_retained_pathway_evidence.tsv",
        EXPIMAP_SOURCE / "table_s35_retained_pathway_member_gene_effects.tsv.gz",
        GENERATIVE_SOURCE / "table_s11_all_random_effects_bh_fdr_genes.tsv",
        GENERATIVE_SOURCE / "table_s16_promoted_gene_literature_annotations.tsv",
        GENERATIVE_SOURCE / "table_s22_matched_gene_literature_annotations.tsv",
        GENERATIVE_SOURCE / "table_s23_grouped_pathway_literature_annotations.tsv",
        IMPORTANCE_SOURCE / "importance_summary.tsv.gz",
        IMPORTANCE_SOURCE / "arm_vs_real_gene_comparison.tsv.gz",
        IMPORTANCE_SOURCE / "selected_arm_vs_real_gene_comparison.tsv.gz",
        *ARM_CHOICE_SOURCES.values(),
        GENE_MAP,
    ]
    manifest = {
        "sources": {
            str(path.relative_to(ROOT)): _sha256(path) for path in source_paths
        },
        "outputs": {
            filename: {
                "rows": int(len(tables[key])),
                "sha256": _sha256(output_dir / filename),
            }
            for key, filename in outputs.items()
        },
        "workbook": {
            "path": workbook.name,
            "sha256": _sha256(workbook),
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return tables


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(args.output_dir)


if __name__ == "__main__":
    main()
