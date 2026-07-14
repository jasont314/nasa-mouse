"""Build source tables and figures for the ASGSR expiMap HVG paper."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
PAPER_DIR = ROOT / "paper/asgsr_expimap_hvg"
SOURCE_DIR = PAPER_DIR / "source_data"
FIGURE_DIR = PAPER_DIR / "figures"


@dataclass(frozen=True)
class ModelConfig:
    tissue: str
    display_name: str
    run_dir: Path
    input_manifest: Path
    reference_summary: Path
    query_input: Path
    literature_labels: Path
    confounded_accessions: tuple[str, ...] = ()
    covariate_input: Path | None = None
    primary_excluded_accessions: tuple[str, ...] = ()
    primary_exclusion_reason: str = ""
    full_input_run_dir: Path | None = None


CONFIGS = (
    ModelConfig(
        tissue="thymus",
        display_name="Thymus",
        run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_thymus/tutorial_hvg_2000/query_nb_250epoch_seed2020",
        input_manifest=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_thymus/tutorial_hvg_2000/input/tutorial_hvg_input_manifest.json",
        reference_summary=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_thymus/tutorial_hvg_2000/reference_nb_400epoch_seed2020/training_summary.json",
        query_input=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_thymus/tutorial_hvg_2000/input/osdr_thymus_query_tutorial_hvg_raw_counts.h5ad",
        literature_labels=ROOT
        / "presentation/expimap/literature_reviewed_hvg/thymus_hvg_literature_review_labels.tsv",
        confounded_accessions=("OSD-289",),
    ),
    ModelConfig(
        tissue="skin",
        display_name="Skin",
        run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_skin/tutorial_hvg_2000/query_nb_250epoch_seed2020",
        input_manifest=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_skin/tutorial_hvg_2000/input/tutorial_hvg_input_manifest.json",
        reference_summary=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_skin/tutorial_hvg_2000/reference_nb_400epoch_seed2020/training_summary.json",
        query_input=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_skin/tutorial_hvg_2000/input/osdr_skin_query_tutorial_hvg_raw_counts.h5ad",
        literature_labels=ROOT
        / "presentation/expimap/literature_reviewed_hvg/skin_hvg_literature_review_labels.tsv",
    ),
    ModelConfig(
        tissue="liver",
        display_name="Liver",
        run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/query_nb_250epoch_seed2020_primary_deduplicated",
        input_manifest=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/input/tutorial_hvg_input_manifest.json",
        reference_summary=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/reference_nb_400epoch_seed2020/training_summary.json",
        query_input=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/input/osdr_liver_query_tutorial_hvg_primary_deduplicated_raw_counts.h5ad",
        literature_labels=ROOT
        / "presentation/expimap/literature_reviewed_hvg/liver_hvg_literature_review_labels.tsv",
        covariate_input=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/input/osdr_liver_query_tutorial_hvg_raw_counts.h5ad",
        primary_excluded_accessions=("OSD-164", "OSD-168"),
        primary_exclusion_reason=(
            "OSD-164 overlaps OSD-47 animals; OSD-168 repackages RR-1/RR-3 "
            "cohorts represented by OSD-48/OSD-137 and includes ERCC technical variants."
        ),
        full_input_run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/query_nb_250epoch_seed2020",
    ),
    ModelConfig(
        tissue="soleus",
        display_name="Soleus",
        run_dir=ROOT
        / "outputs/expimap_muscle_targeted_combined_min8/tutorial_hvg_soleus_2000/query_nb_250epoch_seed2020",
        input_manifest=ROOT
        / "outputs/expimap_muscle_targeted_combined_min8/tutorial_hvg_soleus_2000/input/tutorial_hvg_input_manifest.json",
        reference_summary=ROOT
        / "outputs/expimap_muscle_targeted_combined_min8/tutorial_hvg_edl_2000/reference_nb_400epoch_seed2020/training_summary.json",
        query_input=ROOT
        / "outputs/expimap_muscle_targeted_combined_min8/tutorial_hvg_soleus_2000/input/osdr_skeletal_muscle_soleus_query_tutorial_hvg_raw_counts.h5ad",
        literature_labels=ROOT
        / "presentation/expimap/literature_reviewed_hvg/soleus_hvg_literature_review_labels.tsv",
        confounded_accessions=("OSD-714",),
    ),
)


CURATED_PATHWAYS = {
    "thymus": (
        (
            "R-MMU-69278_CELL_CYCLE_MITOTIC",
            "Mitotic cell cycle",
            "aligned",
            "Reduced thymocyte proliferation",
            "Horie2019;Gridley2013",
        ),
        (
            "R-MMU-73894_DNA_REPAIR",
            "DNA repair",
            "complementary",
            "Reduced DNA-repair program",
            "Gridley2013;Luxton2020",
        ),
        (
            "R-MMU-202403_TCR_SIGNALING",
            "T-cell receptor signaling",
            "aligned",
            "Reduced adaptive immune signaling",
            "Lebsack2010;Gridley2009",
        ),
        (
            "R-MMU-166166_MYD88_INDEPENDENT_TLR4_CASCADE",
            "Innate TLR signaling",
            "complementary",
            "Innate inflammatory signaling",
            "Shimizu2023;Okamura2024",
        ),
        (
            "R-MMU-198933_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
            "Lymphoid-stromal interactions",
            "complementary",
            "Weaker thymocyte-niche coordination",
            "Horie2019;Grandke2026",
        ),
        (
            "R-MMU-8980692_RHOA_GTPASE_CYCLE",
            "RHOA cytoskeletal cycle",
            "complementary",
            "Reduced migration and adhesion signaling",
            "Horie2019;Grandke2026",
        ),
        (
            "R-MMU-1474244_EXTRACELLULAR_MATRIX_ORGANIZATION",
            "Extracellular matrix organization",
            "complementary",
            "Thymic stromal-scaffold response",
            "Grandke2026",
        ),
    ),
    "skin": (
        (
            "R-MMU-6805567_KERATINIZATION",
            "Keratinization",
            "aligned",
            "Lower in pooled skin and under true microgravity at both MHU-2 sites",
            "Cope2024;Neutelings2015",
        ),
        (
            "R-MMU-73894_DNA_REPAIR",
            "DNA repair",
            "context_sensitive",
            "Lower under true microgravity at both MHU-2 sites but heterogeneous across missions",
            "Cope2024;Park2024",
        ),
        (
            "R-MMU-421270_CELL_CELL_JUNCTION_ORGANIZATION",
            "Cell-cell junction organization",
            "aligned",
            "Lower pooled and true-microgravity barrier and tissue-coordination program",
            "Park2024;Cope2024",
        ),
        (
            "R-MMU-3247509_CHROMATIN_MODIFYING_ENZYMES",
            "Chromatin-modifying enzymes",
            "complementary",
            "Lower true-microgravity chromatin-regulatory program at both MHU-2 sites",
            "Cope2024;Park2024",
        ),
        (
            "R-MMU-190828_GAP_JUNCTION_TRAFFICKING",
            "Gap-junction trafficking",
            "complementary",
            "Lower true-microgravity direct cell-communication program at both MHU-2 sites",
            "Park2024;Zhao2025",
        ),
        (
            "R-MMU-156580_PHASE_II_CONJUGATION_OF_COMPOUNDS",
            "Phase II detoxification",
            "context_sensitive",
            "Lower broad phase-II node but higher fully nested glutathione-conjugation node",
            "Mao2014;Cope2024",
        ),
        (
            "R-MMU-5358351_SIGNALING_BY_HEDGEHOG",
            "Hedgehog signaling",
            "complementary",
            "Lower true-microgravity epithelial and follicular signaling at both MHU-2 sites",
            "Neutelings2015;Cope2024",
        ),
        (
            "R-MMU-428157_SPHINGOLIPID_METABOLISM",
            "Sphingolipid metabolism",
            "complementary",
            "Lower true-microgravity barrier-lipid program at both MHU-2 sites",
            "Mao2014;Cope2024",
        ),
    ),
    "liver": (
        (
            "R-MMU-211897_CYTOCHROME_P450_ARRANGED_BY_SUBSTRATE_TYPE",
            "Cytochrome P450",
            "context_sensitive",
            "Higher de-duplicated mean but evenly divided independent cohort directions",
            "Moskaleva2015",
        ),
        (
            "R-MMU-422356_REGULATION_OF_INSULIN_SECRETION",
            "Regulation of insulin secretion",
            "aligned",
            "Lower insulin-regulatory program in 9 of 10 de-duplicated cohort sources",
            "Mathyk2024",
        ),
        (
            "R-MMU-156590_GLUTATHIONE_CONJUGATION",
            "Glutathione conjugation",
            "context_sensitive",
            "Weakly higher de-duplicated mean with heterogeneous directions and metabolite conflict",
            "Kurosawa2021",
        ),
        (
            "R-MMU-2132295_MHC_CLASS_II_ANTIGEN_PRESENTATION",
            "MHC class II antigen presentation",
            "complementary",
            "Lower adaptive immune communication in 8 of 10 de-duplicated cohort sources",
            "Shimizu2023;daSilveira2020",
        ),
        (
            "R-MMU-202403_TCR_SIGNALING",
            "T-cell receptor signaling",
            "complementary",
            "Lower adaptive immune signaling in 9 of 10 de-duplicated cohort sources",
            "Shimizu2023;daSilveira2020",
        ),
        (
            "R-MMU-9012999_RHO_GTPASE_CYCLE",
            "Rho-family GTPase cycle",
            "complementary",
            "Lower mechanosensitive cytoskeletal program in 9 of 10 de-duplicated cohort sources",
            "Li2026;Jonscher2016",
        ),
        (
            "R-MMU-1474244_EXTRACELLULAR_MATRIX_ORGANIZATION",
            "Extracellular matrix organization",
            "context_sensitive",
            "Predominantly lower coordinated matrix-organization and maintenance program after de-duplication",
            "Jonscher2016;Grandke2026",
        ),
    ),
    "soleus": (
        (
            "R-MMU-390522_STRIATED_MUSCLE_CONTRACTION",
            "Striated-muscle contraction",
            "context_sensitive",
            "Contractile direction changes in the restricted sensitivity",
            "Gambara2017;Tascher2017;Tsuji2026",
        ),
        (
            "R-MMU-8978868_FATTY_ACID_METABOLISM",
            "Fatty-acid metabolism",
            "context_sensitive",
            "Metabolic direction changes in the restricted sensitivity",
            "Tascher2017;Mathyk2024",
        ),
        (
            "R-MMU-168256_IMMUNE_SYSTEM",
            "Immune-system signaling",
            "context_sensitive",
            "Inflammatory signal attenuates in the restricted sensitivity",
            "Gambara2017",
        ),
        (
            "R-MMU-1280215_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
            "Cytokine signaling",
            "context_sensitive",
            "Study-consistent decrease conflicts with prior inflammatory direction",
            "Sandona2012;Gambara2017",
        ),
        (
            "R-MMU-1474228_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
            "Extracellular matrix degradation",
            "complementary",
            "Increased matrix turnover",
            "Tascher2017;Murgia2024",
        ),
        (
            "R-MMU-1630316_GLYCOSAMINOGLYCAN_METABOLISM",
            "Glycosaminoglycan metabolism",
            "complementary",
            "Altered muscle-support matrix",
            "Tascher2017;Murgia2024",
        ),
        (
            "R-MMU-73894_DNA_REPAIR",
            "DNA repair",
            "complementary",
            "Cellular stress response",
            "Gambara2017;Tascher2017",
        ),
    ),
}


ROLE_COLORS = {
    "aligned": "#24834b",
    "complementary": "#2869a8",
    "context_sensitive": "#b36a19",
}
TISSUE_COLORS = {
    "thymus": "#7855a6",
    "skin": "#bd5b2a",
    "liver": "#1a7a73",
    "soleus": "#4c6584",
}

SKIN_PROJECT_LOOKUP = {
    "OSD-238": "MHU-2",
    "OSD-239": "MHU-2",
    "OSD-240": "RR-5",
    "OSD-241": "RR-5",
    "OSD-243": "RR-6",
    "OSD-254": "RR-7",
}


def bh_fdr(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = np.where(np.isfinite(values), values, 1.0)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def parse_gmt(path: Path) -> dict[str, set[str]]:
    result = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                result[fields[0]] = set(fields[2:])
    return result


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def latent_directions(config: ModelConfig) -> pd.DataFrame:
    import anndata as ad

    path = config.run_dir / "latent_enrich_condition/model_adata_with_condition_bf.h5ad"
    adata = ad.read_h5ad(path, backed="r")
    terms = list(map(str, adata.uns["terms"]))
    directions = np.asarray(adata.uns["directions"], dtype=float)
    if len(terms) != len(directions):
        raise RuntimeError(f"Latent directions do not match terms for {config.tissue}")
    return pd.DataFrame(
        {
            "term": terms,
            "latent_orientation": directions,
            "active_latent_program": directions != 0,
        }
    )


def model_summary(config: ModelConfig) -> dict:
    manifest = load_json(config.input_manifest)
    training = load_json(config.reference_summary)
    mapping = load_json(config.run_dir / "query_mapping_summary.json")
    scores = pd.read_csv(config.run_dir / "query_pathway_scores.tsv", sep="\t")
    source_scores = scores
    if config.full_input_run_dir is not None:
        source_scores = pd.read_csv(
            config.full_input_run_dir / "query_pathway_scores.tsv", sep="\t"
        )
    return {
        "tissue": config.tissue,
        "reference_samples": manifest["n_reference_samples"],
        "reference_series": training["n_conditions"],
        "query_samples": int(len(scores)),
        "query_flight": int(scores["condition_inferred"].eq("flight").sum()),
        "query_ground_control": int(
            scores["condition_inferred"].eq("ground_control").sum()
        ),
        "query_accessions": int(scores["id.accession"].nunique()),
        "source_query_samples_before_primary_filter": int(len(source_scores)),
        "source_query_accessions_before_primary_filter": int(
            source_scores["id.accession"].nunique()
        ),
        "primary_excluded_accessions": ";".join(
            config.primary_excluded_accessions
        ),
        "hvg_requested": manifest["n_top_genes_requested"],
        "genes_after_filter": manifest["n_genes_after_term_filter"],
        "reactome_programs": manifest["n_terms_after_hvg_filter"],
        "reference_epochs_requested": training["epochs"],
        "reference_epochs_completed": training["training"]["epochs_completed"],
        "query_epochs": mapping["epochs"],
        "reconstruction_loss": training["recon_loss"],
        "hidden_layers": "x".join(map(str, training["hidden_layers"])),
        "gpu": mapping["torch"]["cuda_device_name"],
        "posterior_mean_scores": mapping["posterior_mean_latent"],
    }


def gene_level_results(config: ModelConfig) -> pd.DataFrame:
    import anndata as ad
    from scipy import stats

    adata = ad.read_h5ad(config.query_input)
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    library_sizes = matrix.sum(axis=1, keepdims=True)
    library_sizes[library_sizes == 0] = 1.0
    log2_cpm = np.log2(matrix / library_sizes * 1_000_000 + 1.0)

    condition = adata.obs["condition_inferred"].astype(str).to_numpy()
    flight = condition == "flight"
    ground = condition == "ground_control"
    pooled_effect = log2_cpm[flight].mean(axis=0) - log2_cpm[ground].mean(axis=0)
    pooled_p = np.array(
        [
            stats.ttest_ind(
                log2_cpm[flight, index],
                log2_cpm[ground, index],
                equal_var=False,
                nan_policy="omit",
            ).pvalue
            for index in range(log2_cpm.shape[1])
        ]
    )

    accession_effects = []
    accession_names = []
    for accession, indexes in adata.obs.groupby("id.accession", observed=True).indices.items():
        indexes = np.asarray(indexes)
        local_condition = condition[indexes]
        local_flight = local_condition == "flight"
        local_ground = local_condition == "ground_control"
        if not local_flight.any() or not local_ground.any():
            continue
        accession_names.append(str(accession))
        accession_effects.append(
            log2_cpm[indexes[local_flight]].mean(axis=0)
            - log2_cpm[indexes[local_ground]].mean(axis=0)
        )
    accession_effects = np.vstack(accession_effects)
    accession_mean = accession_effects.mean(axis=0)
    accession_p = np.array(
        [
            stats.ttest_1samp(accession_effects[:, index], 0.0).pvalue
            for index in range(log2_cpm.shape[1])
        ]
    )
    signs = np.sign(accession_effects)
    result = pd.DataFrame(
        {
            "tissue": config.tissue,
            "gene_id": adata.var_names.astype(str),
            "gene_symbol": adata.var.get(
                "gene_symbol", pd.Series(adata.var_names, index=adata.var_names)
            ).astype(str).to_numpy(),
            "pooled_log2cpm_flight_minus_ground": pooled_effect,
            "pooled_p": pooled_p,
            "pooled_fdr": bh_fdr(pooled_p),
            "study_balanced_log2cpm_flight_minus_ground": accession_mean,
            "study_t_p": accession_p,
            "study_t_fdr": bh_fdr(accession_p),
            "n_accessions": len(accession_names),
            "n_accessions_positive": (signs > 0).sum(axis=0),
            "n_accessions_negative": (signs < 0).sum(axis=0),
        }
    )
    return result


def covariate_audit(config: ModelConfig) -> pd.DataFrame:
    import anndata as ad

    source = config.covariate_input or config.query_input
    obs = ad.read_h5ad(source, backed="r").obs.copy()
    rows = []
    for accession, frame in obs.groupby("id.accession", observed=True):
        record = {
            "tissue": config.tissue,
            "accession": str(accession),
            "n_samples": len(frame),
            "n_flight": int(frame["condition_inferred"].astype(str).eq("flight").sum()),
            "n_ground_control": int(
                frame["condition_inferred"].astype(str).eq("ground_control").sum()
            ),
            "sex": ";".join(
                sorted(frame["study.characteristics.sex"].astype(str).unique())
            ),
            "strain": ";".join(
                sorted(frame["study.characteristics.strain"].astype(str).unique())
            ),
            "genotype": ";".join(
                sorted(frame["study.characteristics.genotype"].astype(str).unique())
            ),
            "age_available_in_analysis_metadata": False,
        }
        for field, short in (
            ("study.characteristics.sex", "sex"),
            ("study.characteristics.strain", "strain"),
            ("study.characteristics.genotype", "genotype"),
        ):
            flight_levels = set(
                frame.loc[
                    frame["condition_inferred"].astype(str).eq("flight"), field
                ].astype(str)
            )
            ground_levels = set(
                frame.loc[
                    frame["condition_inferred"].astype(str).eq("ground_control"), field
                ].astype(str)
            )
            record[f"condition_{short}_disjoint"] = bool(
                flight_levels and ground_levels and flight_levels.isdisjoint(ground_levels)
            )
        record["excluded_in_restricted_sensitivity"] = str(accession) in set(
            config.confounded_accessions
        )
        record["included_in_primary_analysis"] = str(accession) not in set(
            config.primary_excluded_accessions
        )
        record["primary_exclusion_reason"] = (
            config.primary_exclusion_reason
            if str(accession) in set(config.primary_excluded_accessions)
            else ""
        )
        rows.append(record)
    return pd.DataFrame(rows)


def pathway_results(
    config: ModelConfig,
    gene_results: pd.DataFrame,
    gmt: dict[str, set[str]],
    directions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pooled = pd.read_csv(
        config.run_dir / "analysis/flt_vs_gc_pathway_comparison.tsv", sep="\t"
    )
    study = pd.read_csv(
        config.run_dir / "analysis/flight_ground_study_aware_tests.tsv", sep="\t"
    ).rename(columns={"n_accessions": "n_accessions_study_test"})
    meta = pd.read_csv(
        config.run_dir / "accession_validation/random_effects_meta_analysis.tsv",
        sep="\t",
    ).rename(columns={"n_accessions": "n_accessions_meta"})
    loo = pd.read_csv(
        config.run_dir / "accession_validation/leave_one_out_summary.tsv", sep="\t"
    )
    accession = pd.read_csv(
        config.run_dir / "accession_validation/per_accession_effects.tsv", sep="\t"
    )
    direction_lookup = directions.set_index("term")["latent_orientation"]
    for frame, columns in (
        (pooled, ("mean_flight", "mean_ground_control", "flight_minus_ground")),
        (study, ("mean_accession_effect", "median_accession_effect")),
        (meta, ("meta_effect",)),
        (accession, ("flight_minus_ground",)),
    ):
        orientation = frame["term"].map(direction_lookup).astype(float)
        for column in columns:
            frame[f"raw_{column}"] = frame[column]
            frame[column] = frame[column] * orientation
            frame.loc[orientation.eq(0), column] = np.nan
    labels = pd.read_csv(config.literature_labels, sep="\t")[
        [
            "term",
            "reviewed_category",
            "reviewed_category_label",
            "review_rationale",
            "citations",
        ]
    ]
    labels = labels.rename(
        columns={
            column: f"legacy_unoriented_{column}"
            for column in labels.columns
            if column != "term"
        }
    )
    merged = (
        pooled.merge(study, on="term", how="left")
        .merge(meta, on="term", how="left")
        .merge(loo, on="term", how="left")
        .merge(directions, on="term", how="left")
        .merge(labels, on="term", how="left")
    )
    merged.insert(0, "tissue", config.tissue)
    merged["absolute_study_balanced_effect"] = merged["mean_accession_effect"].abs()
    merged["within_tissue_magnitude_percentile"] = merged[
        "absolute_study_balanced_effect"
    ].rank(pct=True, method="average")

    clean = accession.loc[
        ~accession["id.accession"].astype(str).isin(config.confounded_accessions)
    ]
    clean_summary = (
        clean.groupby("term")["flight_minus_ground"]
        .agg(
            restricted_study_balanced_effect="mean",
            restricted_effect_min="min",
            restricted_effect_max="max",
            restricted_accessions="size",
            restricted_positive=lambda values: int((values > 0).sum()),
            restricted_negative=lambda values: int((values < 0).sum()),
        )
        .reset_index()
    )
    merged = merged.merge(clean_summary, on="term", how="left")
    merged["restricted_direction_matches_full"] = (
        np.sign(merged["restricted_study_balanced_effect"])
        == np.sign(merged["mean_accession_effect"])
    )

    if config.full_input_run_dir is not None:
        full_input = pd.read_csv(
            config.full_input_run_dir
            / "accession_validation/per_accession_effects.tsv",
            sep="\t",
        )
        full_orientation = full_input["term"].map(direction_lookup).astype(float)
        full_input["flight_minus_ground"] = (
            full_input["flight_minus_ground"] * full_orientation
        )
        full_summary = (
            full_input.groupby("term")["flight_minus_ground"]
            .agg(
                full_input_study_balanced_effect="mean",
                full_input_accessions="size",
                full_input_positive=lambda values: int((values > 0).sum()),
                full_input_negative=lambda values: int((values < 0).sum()),
            )
            .reset_index()
        )
        merged = merged.merge(full_summary, on="term", how="left")

    curated = pd.DataFrame(
        CURATED_PATHWAYS[config.tissue],
        columns=[
            "term",
            "short_label",
            "evidence_role",
            "paper_interpretation",
            "paper_citations",
        ],
    )
    selected = curated.merge(merged, on="term", how="left", validate="one_to_one")
    if selected["mean_accession_effect"].isna().any():
        missing = selected.loc[selected["mean_accession_effect"].isna(), "term"].tolist()
        raise RuntimeError(f"Missing curated pathways for {config.tissue}: {missing}")

    # The original heatmap labels were assigned before latent-sign orientation and
    # are retained only in the all-pathway provenance table. Excluding them here
    # prevents superseded directions from being mistaken for the paper review.
    selected = selected.drop(
        columns=[
            column
            for column in selected.columns
            if column.startswith("legacy_unoriented_")
        ]
    )

    gene_lookup = gene_results.set_index("gene_id")
    support_rows = []
    for row in selected.itertuples(index=False):
        measured = sorted(set(gmt.get(row.term, set())).intersection(gene_lookup.index))
        genes = gene_lookup.loc[measured] if measured else gene_lookup.iloc[0:0]
        significant = genes.loc[genes["pooled_fdr"].lt(0.05)]
        direction = np.sign(row.mean_accession_effect)
        same = int(
            (np.sign(significant["pooled_log2cpm_flight_minus_ground"]) == direction).sum()
        )
        opposite = int(
            (np.sign(significant["pooled_log2cpm_flight_minus_ground"]) == -direction).sum()
        )
        top = significant.assign(
            absolute_effect=significant["pooled_log2cpm_flight_minus_ground"].abs()
        ).sort_values("absolute_effect", ascending=False)
        examples = "; ".join(
            f"{item.gene_symbol}:{item.pooled_log2cpm_flight_minus_ground:+.2f}"
            for item in top.head(5).itertuples()
        )
        support_rows.append(
            {
                "term": row.term,
                "pathway_genes_in_hvg_model": len(measured),
                "pathway_genes_pooled_fdr_lt_005": len(significant),
                "fraction_pathway_genes_pooled_fdr_lt_005": (
                    len(significant) / len(measured) if measured else np.nan
                ),
                "significant_genes_same_direction": same,
                "significant_genes_opposite_direction": opposite,
                "top_gene_examples": examples,
            }
        )
    selected = selected.merge(pd.DataFrame(support_rows), on="term", how="left")
    return merged, selected


def skin_protocol_context_results(
    config: ModelConfig, directions: pd.DataFrame
) -> pd.DataFrame:
    scores = pd.read_csv(config.run_dir / "query_pathway_scores.tsv", sep="\t")
    curated = pd.DataFrame(
        CURATED_PATHWAYS["skin"],
        columns=[
            "term",
            "short_label",
            "evidence_role",
            "paper_interpretation",
            "paper_citations",
        ],
    )
    terms = curated["term"].tolist()
    orientation = directions.set_index("term")["latent_orientation"]
    for term in terms:
        scores[term] = scores[term].astype(float) * float(orientation.loc[term])

    rows = []

    def add_contrast(
        frame: pd.DataFrame,
        accession: str,
        display_label: str,
        context_group: str,
        contrast_order: int,
        site: str,
        gravity: str,
        collection: str,
    ) -> None:
        flight = frame.loc[frame["condition_inferred"].eq("flight")]
        ground = frame.loc[frame["condition_inferred"].eq("ground_control")]
        if flight.empty or ground.empty:
            raise RuntimeError(f"Incomplete skin context contrast: {display_label}")
        for pathway_order, pathway in enumerate(curated.itertuples(index=False)):
            rows.append(
                {
                    "tissue": "skin",
                    "accession": accession,
                    "display_label": display_label,
                    "context_group": context_group,
                    "contrast_order": contrast_order,
                    "site": site,
                    "gravity": gravity,
                    "collection_context": collection,
                    "term": pathway.term,
                    "short_label": pathway.short_label,
                    "pathway_order": pathway_order,
                    "n_flight": int(len(flight)),
                    "n_ground_control": int(len(ground)),
                    "flight_minus_ground": float(
                        flight[pathway.term].mean() - ground[pathway.term].mean()
                    ),
                }
            )

    # MHU-2 includes true microgravity and onboard centrifuge-generated 1 g.
    # Each flight subgroup is compared with the same Earth 1 g control animals.
    for accession, site, start_order in (
        ("OSD-238", "dorsal skin", 0),
        ("OSD-239", "femoral skin", 2),
    ):
        frame = scores.loc[scores["id.accession"].eq(accession)]
        ground = frame.loc[frame["condition_inferred"].eq("ground_control")]
        for offset, token, short_gravity in (
            (0, "FLT_uG", "microgravity"),
            (1, "FLT_1G", "onboard artificial 1 g"),
        ):
            flight = frame.loc[
                frame["id.sample name"].astype(str).str.contains(token, regex=False)
            ]
            contrast = pd.concat([flight, ground], ignore_index=True)
            add_contrast(
                contrast,
                accession,
                f"{site.split()[0].title()} {short_gravity.replace('onboard ', '')}",
                "MHU-2 gravity and skin site",
                start_order + offset,
                site,
                short_gravity,
                "live return; collected within one day",
            )

    for accession, site, order in (
        ("OSD-240", "dorsal skin", 4),
        ("OSD-241", "femoral skin", 5),
    ):
        add_contrast(
            scores.loc[scores["id.accession"].eq(accession)],
            accession,
            f"RR-5 {site.split()[0]} 30-day recovery",
            "Recovery, endpoint, duration, and strain",
            order,
            site,
            "microgravity exposure",
            "live return followed by approximately 30 days of Earth recovery",
        )

    for token, label, order, collection in (
        ("LAR", "RR-6 live return ~30 days", 6, "live return after approximately 30 days"),
        ("ISS-T", "RR-6 terminal ~60 days", 7, "ISS terminal after approximately 60 days"),
    ):
        frame = scores.loc[
            scores["id.accession"].eq("OSD-243")
            & scores["id.sample name"].astype(str).str.contains(token, regex=False)
        ]
        add_contrast(
            frame,
            "OSD-243",
            label,
            "Recovery, endpoint, duration, and strain",
            order,
            "dorsal skin",
            "microgravity exposure",
            collection,
        )

    order = 8
    for strain_token, strain_label in (("C3H-HeJ", "C3H/HeJ"), ("C57-6J", "C57BL/6J")):
        for duration_token, duration_label in (("25days", "25 days"), ("75days", "75 days")):
            frame = scores.loc[
                scores["id.accession"].eq("OSD-254")
                & scores["id.sample name"].astype(str).str.contains(
                    strain_token, regex=False
                )
                & scores["id.sample name"].astype(str).str.contains(
                    duration_token, regex=False
                )
            ]
            add_contrast(
                frame,
                "OSD-254",
                f"RR-7 {strain_label} {duration_label}",
                "Recovery, endpoint, duration, and strain",
                order,
                "dorsal skin",
                "microgravity exposure",
                f"ISS terminal after {duration_label}",
            )
            order += 1
    return pd.DataFrame(rows).sort_values(["contrast_order", "pathway_order"])


def add_context_aware_annotations(
    selected: pd.DataFrame,
    skin_context: pd.DataFrame,
    accession_effects: pd.DataFrame,
) -> pd.DataFrame:
    """Attach protocol-aware effects without replacing literature evidence roles."""
    annotated = selected.copy()
    annotated["primary_analysis_scope"] = "Study-balanced flight versus ground control"
    annotated["context_aware_annotation"] = ""
    annotated["protocol_consistency_status"] = "No additional protocol restriction"
    annotated["context_sensitive_reason"] = ""

    context_columns = {
        "Dorsal microgravity": "skin_dorsal_microgravity_effect",
        "Femoral microgravity": "skin_femoral_microgravity_effect",
        "Dorsal artificial 1 g": "skin_dorsal_artificial_1g_effect",
        "Femoral artificial 1 g": "skin_femoral_artificial_1g_effect",
    }
    for output_column in context_columns.values():
        annotated[output_column] = np.nan
    annotated["skin_true_microgravity_same_direction"] = pd.NA
    annotated["skin_project_balanced_effect"] = np.nan
    annotated["skin_projects"] = np.nan
    annotated["skin_project_positive"] = np.nan
    annotated["skin_project_negative"] = np.nan

    skin_accession = accession_effects.loc[
        accession_effects["tissue"].eq("skin")
    ].copy()
    skin_accession["project"] = skin_accession["id.accession"].map(
        SKIN_PROJECT_LOOKUP
    )
    skin_project = (
        skin_accession.groupby(["term", "project"], as_index=False)[
            "flight_minus_ground"
        ]
        .mean()
        .groupby("term")["flight_minus_ground"]
        .agg(
            skin_project_balanced_effect="mean",
            skin_projects="size",
            skin_project_positive=lambda values: int((values > 0).sum()),
            skin_project_negative=lambda values: int((values < 0).sum()),
        )
    )

    skin_lookup = skin_context.pivot(
        index="term", columns="display_label", values="flight_minus_ground"
    )
    for index, row in annotated.loc[annotated["tissue"].eq("skin")].iterrows():
        effects = skin_lookup.loc[row["term"]]
        for display_label, output_column in context_columns.items():
            annotated.at[index, output_column] = float(effects[display_label])
        dorsal_microgravity = float(effects["Dorsal microgravity"])
        femoral_microgravity = float(effects["Femoral microgravity"])
        dorsal_artificial = float(effects["Dorsal artificial 1 g"])
        femoral_artificial = float(effects["Femoral artificial 1 g"])
        same_microgravity_direction = (
            np.sign(dorsal_microgravity) == np.sign(femoral_microgravity)
        )
        project = skin_project.loc[row["term"]]
        for column in (
            "skin_project_balanced_effect",
            "skin_projects",
            "skin_project_positive",
            "skin_project_negative",
        ):
            annotated.at[index, column] = project[column]
        annotated.at[index, "skin_true_microgravity_same_direction"] = bool(
            same_microgravity_direction
        )
        annotated.at[index, "primary_analysis_scope"] = (
            "Six-accession broad spaceflight estimate; MHU-2 artificial 1 g is "
            "included in the top-level flight label"
        )
        annotated.at[index, "protocol_consistency_status"] = (
            "MHU-2 true-microgravity direction agrees between dorsal and femoral skin"
        )
        annotated.at[index, "context_aware_annotation"] = (
            f"Pooled effect {row['mean_accession_effect']:+.3f}. True microgravity "
            f"was {dorsal_microgravity:+.3f} in dorsal and "
            f"{femoral_microgravity:+.3f} in femoral skin; onboard artificial 1 g "
            f"was {dorsal_artificial:+.3f} in dorsal and "
            f"{femoral_artificial:+.3f} in femoral skin. After paired sites were "
            f"collapsed to four mission projects, the project-balanced effect was "
            f"{project['skin_project_balanced_effect']:+.3f} "
            f"({int(project['skin_project_negative'])} lower, "
            f"{int(project['skin_project_positive'])} higher). Literature role remains "
            f"{row['evidence_role'].replace('_', ' ')} because protocol consistency "
            "and prior-literature alignment are separate judgments."
        )

    skin_dna = annotated["tissue"].eq("skin") & annotated["short_label"].eq(
        "DNA repair"
    )
    annotated.loc[skin_dna, "context_sensitive_reason"] = (
        "True microgravity is lower at both MHU-2 sites, but the broad six-accession "
        "effects divide 3 lower and 3 higher. Collapsing paired sites gives 3 lower "
        "and 1 higher mission-project effects, so cross-project heterogeneity remains."
    )

    for index, row in annotated.loc[annotated["tissue"].eq("liver")].iterrows():
        primary = float(row["mean_accession_effect"])
        original = float(row["full_input_study_balanced_effect"])
        direction_status = (
            "direction preserved"
            if np.sign(primary) == np.sign(original)
            else "direction changed"
        )
        annotated.at[index, "primary_analysis_scope"] = (
            "De-duplicated 10-accession remap; overlapping OSD-164 and OSD-168 excluded"
        )
        annotated.at[index, "protocol_consistency_status"] = (
            f"{direction_status.title()} relative to the original 12-accession input"
        )
        annotated.at[index, "context_aware_annotation"] = (
            f"De-duplicated effect {primary:+.3f} versus original 12-accession "
            f"effect {original:+.3f}; {direction_status}. Evidence role is based on "
            "the de-duplicated primary result and literature review."
        )

    context_reasons = {
        ("skin", "Phase II detoxification"): (
            "The broad phase-II node is lower, but its fully nested glutathione-"
            "conjugation child is higher and the two scores are strongly anticorrelated; "
            "a uniform detoxification direction is not identifiable."
        ),
        ("liver", "Cytochrome P450"): (
            "The de-duplicated sources divide 5 higher and 5 lower, so the positive "
            "mean is not a common cohort direction."
        ),
        ("liver", "Glutathione conjugation"): (
            "The effect is small, directions divide 6 higher and 4 lower, and a higher "
            "latent score does not match reported lower hepatic glutathione pools."
        ),
        ("liver", "Extracellular matrix organization"): (
            "Seven of 10 sources are lower and the negative effects are larger than the "
            "three positive effects. This suggests reduced coordinated matrix organization "
            "or maintenance, not lower total matrix abundance or absence of remodeling."
        ),
        ("soleus", "Immune-system signaling"): (
            "The all-accession increase attenuates and changes sign after excluding "
            "condition-strain-confounded OSD-714."
        ),
        ("soleus", "Fatty-acid metabolism"): (
            "The all-accession decrease changes to an increase after excluding "
            "condition-strain-confounded OSD-714."
        ),
        ("soleus", "Striated-muscle contraction"): (
            "The all-accession increase changes to a decrease after excluding "
            "condition-strain-confounded OSD-714."
        ),
        ("soleus", "Cytokine signaling"): (
            "All three accessions are lower, but that direction conflicts with increased "
            "inflammatory signaling reported in some prior soleus studies."
        ),
    }
    for (tissue, label), reason in context_reasons.items():
        mask = annotated["tissue"].eq(tissue) & annotated["short_label"].eq(label)
        annotated.loc[mask, "context_sensitive_reason"] = reason

    for tissue, excluded in (("thymus", "OSD-289"), ("soleus", "OSD-714")):
        mask = annotated["tissue"].eq(tissue)
        preserved = annotated.loc[mask, "restricted_direction_matches_full"].map(
            {True: "preserved", False: "changed"}
        )
        annotated.loc[mask, "protocol_consistency_status"] = (
            "Direction " + preserved + f" after excluding {excluded}"
        )

    return annotated


def systematic_pathway_screen(
    pathway: pd.DataFrame,
    selected: pd.DataFrame,
    accession_effects: pd.DataFrame,
) -> pd.DataFrame:
    """Rank every pathway transparently without assigning unreviewed biology labels."""
    import re

    screen = pathway.copy()
    curated = selected[
        ["tissue", "term", "short_label", "evidence_role"]
    ].rename(columns={"short_label": "curated_short_label"})
    screen = screen.merge(curated, on=["tissue", "term"], how="left")
    screen["curated_for_main_figures"] = screen["evidence_role"].notna()
    screen["active_latent_program"] = screen["active_latent_program"].fillna(False).astype(bool)
    same_direction = np.where(
        screen["mean_accession_effect"].ge(0),
        screen["restricted_positive"],
        screen["restricted_negative"],
    )
    screen["primary_direction_accessions"] = same_direction
    screen["primary_direction_agreement"] = same_direction / screen[
        "restricted_accessions"
    ].replace(0, np.nan)
    screen["absolute_effect_rank_active"] = np.nan
    for tissue, indices in screen.loc[screen["active_latent_program"]].groupby(
        "tissue"
    ).groups.items():
        screen.loc[indices, "absolute_effect_rank_active"] = screen.loc[
            indices, "absolute_study_balanced_effect"
        ].rank(method="first", ascending=False)
    screen["top_20_absolute_effect_active"] = screen[
        "absolute_effect_rank_active"
    ].le(20)

    skin_accession = accession_effects.loc[
        accession_effects["tissue"].eq("skin")
    ].copy()
    skin_accession["project"] = skin_accession["id.accession"].map(
        SKIN_PROJECT_LOOKUP
    )
    skin_project = (
        skin_accession.groupby(["term", "project"], as_index=False)[
            "flight_minus_ground"
        ]
        .mean()
        .groupby("term")["flight_minus_ground"]
        .agg(
            skin_project_balanced_effect="mean",
            skin_projects="size",
            skin_project_positive=lambda values: int((values > 0).sum()),
            skin_project_negative=lambda values: int((values < 0).sum()),
        )
        .reset_index()
    )
    screen = screen.merge(skin_project, on="term", how="left")
    skin_mask = screen["tissue"].eq("skin")
    screen.loc[~skin_mask, list(skin_project.columns[1:])] = np.nan
    screen["skin_project_direction_matches_accession"] = pd.NA
    screen.loc[skin_mask, "skin_project_direction_matches_accession"] = (
        np.sign(screen.loc[skin_mask, "skin_project_balanced_effect"])
        == np.sign(screen.loc[skin_mask, "mean_accession_effect"])
    )
    screen["skin_project_effect_change"] = np.nan
    screen.loc[skin_mask, "skin_project_effect_change"] = (
        screen.loc[skin_mask, "skin_project_balanced_effect"]
        - screen.loc[skin_mask, "mean_accession_effect"]
    )
    screen["skin_project_absolute_effect_rank_active"] = np.nan
    skin_active = skin_mask & screen["active_latent_program"]
    screen.loc[skin_active, "skin_project_absolute_effect_rank_active"] = screen.loc[
        skin_active, "skin_project_balanced_effect"
    ].abs().rank(method="first", ascending=False)
    screen["skin_project_top_20_absolute_effect_active"] = screen[
        "skin_project_absolute_effect_rank_active"
    ].le(20)

    def display_name(term: str) -> str:
        name = re.sub(r"^R-MMU-\d+_", "", str(term)).replace("_", " ").title()
        for old, new in (
            ("Dna", "DNA"),
            ("Rna", "RNA"),
            ("Gtpase", "GTPase"),
            ("Mhc", "MHC"),
            ("Tlr", "TLR"),
            ("Er ", "ER "),
            ("Ecm", "ECM"),
            ("Cdh1", "CDH1"),
            ("Gpcr", "GPCR"),
        ):
            name = name.replace(old, new)
        return name

    screen["screen_display_label"] = screen["term"].map(display_name)
    screen.loc[screen["curated_for_main_figures"], "screen_display_label"] = screen.loc[
        screen["curated_for_main_figures"], "curated_short_label"
    ]
    screen["screening_note"] = np.where(
        screen["curated_for_main_figures"],
        "Individually literature-reviewed for the main figures",
        np.where(
            screen["top_20_absolute_effect_active"],
            "High-magnitude screening candidate; biological relevance not yet assigned",
            "Complete-result record; not individually reviewed for the main figures",
        ),
    )
    return screen.sort_values(
        ["tissue", "absolute_effect_rank_active", "term"], na_position="last"
    )


def save_figure(fig, name: str) -> None:
    fig.savefig(FIGURE_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{name}.pdf", bbox_inches="tight")


def plot_workflow(summary: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    fig, ax = plt.subplots(figsize=(13.2, 4.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = (
        (0.02, 0.20, 0.18, 0.62, "ARCHS4 reference", "Tissue-matched mouse\nbulk RNA-seq\nSpaceflight studies excluded"),
        (0.23, 0.20, 0.18, 0.62, "Reactome architecture", "Current mouse Ensembl\npathways\n2,000 reference-selected HVGs"),
        (0.44, 0.20, 0.18, 0.62, "expiMap reference", "Negative-binomial model\n3 x 300 hidden units\nA100 GPU"),
        (0.65, 0.20, 0.16, 0.62, "OSDR query mapping", "Flight and ground-control\nbulk RNA-seq\nAccession encoded as\nquery condition"),
        (0.84, 0.20, 0.14, 0.62, "Robustness", "ssGSEA and GSEA\nHeld-out projects\nThree full seeds\nComposition proxies\nLiterature audit"),
    )
    fills = ("#e7f0f5", "#edf3e8", "#e9e7f3", "#f6ece4", "#eef0f1")
    for index, (x, y, width, height, title, body) in enumerate(boxes):
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.012,rounding_size=0.015",
            linewidth=1.1,
            edgecolor="#37434a",
            facecolor=fills[index],
        )
        ax.add_patch(patch)
        ax.text(x + width / 2, y + height - 0.12, title, ha="center", va="center", fontsize=11, weight="bold")
        ax.text(x + width / 2, y + height / 2 - 0.03, body, ha="center", va="center", fontsize=7.8, linespacing=1.4)
        if index < len(boxes) - 1:
            next_x = boxes[index + 1][0]
            ax.add_patch(
                FancyArrowPatch(
                    (x + width + 0.005, 0.51),
                    (next_x - 0.005, 0.51),
                    arrowstyle="-|>",
                    mutation_scale=14,
                    linewidth=1.2,
                    color="#4b565c",
                )
            )
    counts = "  |  ".join(
        f"{row.tissue.title()}: {int(row.query_samples)} query samples, {int(row.reactome_programs)} programs"
        for row in summary.itertuples()
    )
    ax.text(0.5, 0.08, counts, ha="center", va="center", fontsize=8.8, color="#30383c")
    ax.set_title("Reference-guided, pathway-constrained analysis of NASA OSDR mouse transcriptomes", fontsize=14, weight="bold", pad=12)
    save_figure(fig, "figure_1_workflow")
    plt.close(fig)


def plot_pathway_shifts(selected: pd.DataFrame, accession_effects: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.8), constrained_layout=True)
    for panel, config in enumerate(CONFIGS):
        ax = axes.flat[panel]
        subset = selected.loc[selected["tissue"].eq(config.tissue)].copy()
        subset = subset.sort_values("mean_accession_effect")
        y = np.arange(len(subset))
        for position, row in zip(y, subset.itertuples()):
            points = accession_effects.loc[
                accession_effects["tissue"].eq(config.tissue)
                & accession_effects["term"].eq(row.term)
            ]
            clean = points.loc[
                ~points["id.accession"].astype(str).isin(config.confounded_accessions)
            ]
            flagged = points.loc[
                points["id.accession"].astype(str).isin(config.confounded_accessions)
            ]
            ax.hlines(
                position,
                points["flight_minus_ground"].min(),
                points["flight_minus_ground"].max(),
                color="#c5c9cc",
                linewidth=1.0,
                zorder=1,
            )
            ax.scatter(
                clean["flight_minus_ground"],
                np.full(len(clean), position),
                s=22,
                color="#8d969b",
                alpha=0.78,
                zorder=2,
            )
            if not flagged.empty:
                ax.scatter(
                    flagged["flight_minus_ground"],
                    np.full(len(flagged), position),
                    s=52,
                    facecolor="#f3c676",
                    edgecolor="#9b5b12",
                    linewidth=1.0,
                    zorder=3,
                )
            ax.scatter(
                row.mean_accession_effect,
                position,
                s=72,
                marker="D",
                color=ROLE_COLORS[row.evidence_role],
                edgecolor="white",
                linewidth=0.7,
                zorder=4,
            )
            if config.confounded_accessions:
                ax.scatter(
                    row.restricted_study_balanced_effect,
                    position,
                    s=62,
                    marker="D",
                    facecolor="white",
                    edgecolor="#20282c",
                    linewidth=1.2,
                    zorder=5,
                )
        ax.axvline(0, color="#30383c", linewidth=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(subset["short_label"], fontsize=9)
        ax.set_xlabel("Flight minus ground expiMap pathway score", fontsize=9)
        ax.set_title(
            f"{config.display_name}  |  {int(subset['n_accessions_meta'].iloc[0])} OSD accessions",
            loc="left",
            fontsize=12,
            weight="bold",
        )
        ax.grid(axis="x", color="#e4e6e7", linewidth=0.7)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
    fig.suptitle(
        "Representative pathway shifts and study sensitivity",
        fontsize=15,
        weight="bold",
    )
    role_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="D",
            linestyle="",
            markerfacecolor=color,
            markeredgecolor="white",
            markersize=7,
            label=label,
        )
        for label, color in (
            ("Literature-aligned", ROLE_COLORS["aligned"]),
            ("Complementary", ROLE_COLORS["complementary"]),
            ("Context-sensitive", ROLE_COLORS["context_sensitive"]),
        )
    ]
    fig.legend(handles=role_handles, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.965), fontsize=8.5)
    fig.text(
        0.5,
        -0.015,
        "Small circles: accession-specific effects. Colored diamonds: study-balanced mean. "
        "Hollow diamonds: mean after excluding condition-strain-confounded accessions. "
        "Orange circles: excluded accessions. Liver shows the 10-cohort de-duplicated primary run; "
        "skin uses the broad flight label, with gravity and site separated in Figure 7.",
        ha="center",
        fontsize=9,
    )
    save_figure(fig, "figure_2_tissue_pathway_shifts")
    plt.close(fig)


def plot_evidence_map(selected: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt
    import textwrap

    frame = selected.copy()
    frame["direction"] = np.where(frame["mean_accession_effect"].ge(0), "Higher in flight", "Lower in flight")
    frame["signed_percentile"] = np.sign(frame["mean_accession_effect"]) * frame[
        "within_tissue_magnitude_percentile"
    ]
    frame["gene_support"] = frame["fraction_pathway_genes_pooled_fdr_lt_005"].fillna(0)

    fig = plt.figure(figsize=(15.2, 8.9), constrained_layout=True)
    grid = fig.add_gridspec(2, 4, height_ratios=(1.0, 0.09))
    axes = [fig.add_subplot(grid[0, index]) for index in range(4)]
    footer_ax = fig.add_subplot(grid[1, :])
    footer_ax.axis("off")
    for ax, config in zip(axes, CONFIGS):
        subset = frame.loc[frame["tissue"].eq(config.tissue)].copy()
        subset = subset.sort_values("signed_percentile")
        y = np.arange(len(subset))
        colors = np.where(subset["signed_percentile"].ge(0), "#b53b32", "#2d68a7")
        sizes = 55 + 850 * subset["gene_support"].clip(upper=0.30)
        ax.scatter(subset["signed_percentile"], y, s=sizes, c=colors, alpha=0.88, edgecolor="white", linewidth=0.8)
        ax.axvline(0, color="#40484c", linewidth=0.8)
        # Leave enough horizontal padding for the largest gene-support bubbles.
        ax.set_xlim(-1.14, 1.14)
        ax.set_xticks([-1, -0.5, 0, 0.5, 1])
        ax.set_xticklabels(["large\nlower", "", "0", "", "large\nhigher"], fontsize=8)
        ax.set_yticks(y)
        ax.set_yticklabels(subset["short_label"], fontsize=8.5)
        for tick, role in zip(ax.get_yticklabels(), subset["evidence_role"]):
            tick.set_color(ROLE_COLORS[role])
        ax.set_title(config.display_name, fontsize=12, weight="bold")
        ax.grid(axis="x", color="#e4e6e7", linewidth=0.7)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
    fig.suptitle(
        "Within-tissue pathway magnitude, direction, and gene-level support",
        fontsize=15,
        weight="bold",
    )
    footer = (
        "Position is the signed within-model magnitude percentile; latent-score scales are not compared across tissues. "
        "Circle area reflects the fraction of measured pathway genes with pooled gene-level FDR < 0.05. "
        "Green labels are literature-aligned, blue labels are complementary hypotheses, and orange labels are context-sensitive."
    )
    footer_ax.text(
        0.5,
        0.5,
        textwrap.fill(footer, width=150),
        ha="center",
        va="center",
        fontsize=8.8,
        linespacing=1.2,
        transform=footer_ax.transAxes,
    )
    save_figure(fig, "figure_3_evidence_map")
    plt.close(fig)


def plot_sensitivity(selected: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    frame = selected.loc[
        selected["tissue"].isin(["thymus", "liver", "soleus"])
    ].copy()
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 6.8), constrained_layout=True)
    for ax, tissue in zip(axes, ("thymus", "liver", "soleus")):
        subset = frame.loc[frame["tissue"].eq(tissue)].copy()
        subset = subset.sort_values("mean_accession_effect")
        y = np.arange(len(subset))
        if tissue == "liver":
            comparison = subset["full_input_study_balanced_effect"]
            primary_label = "De-duplicated primary (10 accessions)"
            comparison_label = "Original input (12 accessions)"
            title = "Liver  |  exclude OSD-164 and OSD-168"
        else:
            comparison = subset["restricted_study_balanced_effect"]
            primary_label = "All accessions"
            comparison_label = "Restricted sensitivity"
            dropped = ", ".join(
                next(c.confounded_accessions for c in CONFIGS if c.tissue == tissue)
            )
            title = f"{tissue.title()}  |  exclude {dropped}"
        ax.hlines(
            y,
            subset["mean_accession_effect"],
            comparison,
            color="#b6bcc0",
            linewidth=1.4,
        )
        ax.scatter(
            subset["mean_accession_effect"],
            y,
            s=62,
            color=TISSUE_COLORS[tissue],
            label=primary_label,
            zorder=3,
        )
        ax.scatter(
            comparison,
            y,
            s=62,
            facecolor="white",
            edgecolor="#20282c",
            linewidth=1.1,
            label=comparison_label,
            zorder=4,
        )
        ax.axvline(0, color="#30383c", linewidth=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(subset["short_label"], fontsize=8.4)
        ax.set_title(title, loc="left", fontsize=11.2, weight="bold")
        ax.set_xlabel("Study-balanced flight minus ground score", fontsize=9)
        ax.grid(axis="x", color="#e4e6e7", linewidth=0.7)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.legend(frameon=False, fontsize=7.8, loc="best")
    fig.suptitle(
        "Primary-analysis sensitivity to condition confounding and cohort duplication",
        fontsize=15,
        weight="bold",
    )
    save_figure(fig, "figure_4_primary_analysis_sensitivity")
    plt.close(fig)


def plot_skin_protocol_context(context: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt
    import textwrap

    groups = (
        "MHU-2 gravity and skin site",
        "Recovery, endpoint, duration, and strain",
    )
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16.2, 8.2),
        gridspec_kw={"width_ratios": [1, 2]},
        constrained_layout=True,
    )
    vmax = float(np.ceil(context["flight_minus_ground"].abs().max() * 2) / 2)
    image = None
    for ax, group in zip(axes, groups):
        subset = context.loc[context["context_group"].eq(group)]
        matrix = subset.pivot(
            index="pathway_order",
            columns="contrast_order",
            values="flight_minus_ground",
        ).sort_index(axis=0).sort_index(axis=1)
        contrast_labels = (
            subset[["contrast_order", "display_label"]]
            .drop_duplicates()
            .sort_values("contrast_order")["display_label"]
            .tolist()
        )
        pathway_labels = (
            subset[["pathway_order", "short_label"]]
            .drop_duplicates()
            .sort_values("pathway_order")["short_label"]
            .tolist()
        )
        image = ax.imshow(matrix.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(np.arange(len(contrast_labels)))
        ax.set_xticklabels(
            [textwrap.fill(label, width=16) for label in contrast_labels],
            fontsize=7.8,
            rotation=38 if len(contrast_labels) > 4 else 0,
            ha="right" if len(contrast_labels) > 4 else "center",
            rotation_mode="anchor",
        )
        ax.set_yticks(np.arange(len(pathway_labels)))
        ax.set_yticklabels(pathway_labels, fontsize=9)
        ax.set_title(group, fontsize=11.5, weight="bold", loc="left")
        for row in range(matrix.shape[0]):
            for column in range(matrix.shape[1]):
                value = matrix.iloc[row, column]
                color = "white" if abs(value) > vmax * 0.55 else "#202628"
                ax.text(
                    column,
                    row,
                    f"{value:+.2f}",
                    ha="center",
                    va="center",
                    fontsize=7.7,
                    color=color,
                )
        for spine in ax.spines.values():
            spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=axes, fraction=0.025, pad=0.015)
    colorbar.set_label("Decoder-oriented flight minus matched ground score")
    fig.suptitle(
        "Skin pathway shifts vary with gravity and collection context",
        fontsize=15,
        weight="bold",
    )
    fig.text(
        0.5,
        -0.015,
        "MHU-2 microgravity and onboard artificial 1 g subgroups share Earth 1 g controls. "
        "RR-5 was collected after about 30 days of Earth recovery; RR-6 and RR-7 include live-return or ISS-terminal contrasts.",
        ha="center",
        fontsize=8.8,
    )
    save_figure(fig, "figure_7_skin_protocol_context")
    plt.close(fig)


def plot_broad_pathway_screen(screen: pd.DataFrame) -> None:
    """Show the highest-magnitude active programs without implying validation."""
    import matplotlib.pyplot as plt
    import textwrap

    fig, axes = plt.subplots(2, 2, figsize=(17.2, 15.2), constrained_layout=True)
    for ax, tissue in zip(axes.flat, ("thymus", "skin", "liver", "soleus")):
        subset = screen.loc[
            screen["tissue"].eq(tissue)
            & screen["top_20_absolute_effect_active"]
        ].copy()
        subset = subset.sort_values("mean_accession_effect")
        y = np.arange(len(subset))
        colors = np.where(
            subset["mean_accession_effect"].ge(0), "#b3483f", "#326ea8"
        )
        ax.scatter(
            subset["mean_accession_effect"],
            y,
            c=colors,
            s=np.where(subset["curated_for_main_figures"], 70, 42),
            edgecolor=np.where(
                subset["curated_for_main_figures"], "#171b1d", "white"
            ),
            linewidth=np.where(subset["curated_for_main_figures"], 1.2, 0.6),
            zorder=3,
        )
        ax.axvline(0, color="#30383c", linewidth=0.8)
        labels = [
            textwrap.fill(
                f"{label}{' *' if curated else ''}", width=44
            )
            for label, curated in zip(
                subset["screen_display_label"],
                subset["curated_for_main_figures"],
            )
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7.2)
        ax.set_title(
            f"{tissue.title()}  |  top 20 of "
            f"{int(screen.loc[screen['tissue'].eq(tissue), 'active_latent_program'].sum())} active programs",
            loc="left",
            fontsize=11.5,
            weight="bold",
        )
        ax.set_xlabel("Study-balanced flight minus ground score", fontsize=9)
        ax.grid(axis="x", color="#e4e6e7", linewidth=0.7)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
    fig.suptitle(
        "Broad pathway screen: highest absolute effects per tissue",
        fontsize=15,
        weight="bold",
    )
    fig.text(
        0.5,
        -0.01,
        "Programs are ranked only by absolute study-balanced effect among active latent dimensions. "
        "Asterisks and black outlines identify pathways individually reviewed for the main figures. "
        "Unmarked terms are screening candidates, not validated biological claims; related Reactome terms are not independent.",
        ha="center",
        fontsize=8.8,
        wrap=True,
    )
    save_figure(fig, "figure_s1_broad_pathway_screen")
    plt.close(fig)


def plot_skin_project_balance(screen: pd.DataFrame) -> None:
    """Compare six-accession and four-project effects for every active skin program."""
    import matplotlib.pyplot as plt
    import textwrap

    skin = screen.loc[
        screen["tissue"].eq("skin") & screen["active_latent_program"]
    ].copy()
    same_direction = skin["skin_project_direction_matches_accession"].astype(bool)
    correlation = skin["mean_accession_effect"].corr(
        skin["skin_project_balanced_effect"]
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16.4, 8.8),
        gridspec_kw={"width_ratios": [1.0, 1.18]},
        constrained_layout=True,
    )
    ax = axes[0]
    limit = float(
        np.ceil(
            skin[
                ["mean_accession_effect", "skin_project_balanced_effect"]
            ].abs().to_numpy().max()
            * 10
        )
        / 10
    )
    ax.plot([-limit, limit], [-limit, limit], color="#9aa2a6", linestyle="--", linewidth=1)
    ax.axhline(0, color="#30383c", linewidth=0.7)
    ax.axvline(0, color="#30383c", linewidth=0.7)
    unreviewed = skin.loc[~skin["curated_for_main_figures"]]
    reviewed = skin.loc[skin["curated_for_main_figures"]]
    ax.scatter(
        unreviewed["mean_accession_effect"],
        unreviewed["skin_project_balanced_effect"],
        s=24,
        color="#8b969b",
        alpha=0.62,
        edgecolor="white",
        linewidth=0.4,
        label="Unreviewed active program",
    )
    ax.scatter(
        reviewed["mean_accession_effect"],
        reviewed["skin_project_balanced_effect"],
        s=68,
        color=TISSUE_COLORS["skin"],
        edgecolor="#171b1d",
        linewidth=1.0,
        label="Main-figure reviewed program",
        zorder=4,
    )
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Six-accession study-balanced effect")
    ax.set_ylabel("Four-project-balanced effect")
    ax.set_title(
        "All active skin programs",
        loc="left",
        fontsize=12,
        weight="bold",
    )
    ax.text(
        0.03,
        0.97,
        f"Pearson r = {correlation:.3f}\n"
        f"Same direction: {int(same_direction.sum())}/{len(skin)}\n"
        f"Sign changes: {len(skin) - int(same_direction.sum())}, all |effect| < 0.036",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "#c7cccf", "pad": 5},
    )
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.grid(color="#e5e7e8", linewidth=0.6)
    ax.set_axisbelow(True)

    ax = axes[1]
    top = skin.loc[skin["skin_project_top_20_absolute_effect_active"]].copy()
    top = top.sort_values("skin_project_balanced_effect")
    y = np.arange(len(top))
    colors = np.where(
        top["skin_project_balanced_effect"].ge(0), "#b3483f", "#326ea8"
    )
    ax.scatter(
        top["skin_project_balanced_effect"],
        y,
        c=colors,
        s=np.where(top["curated_for_main_figures"], 72, 44),
        edgecolor=np.where(top["curated_for_main_figures"], "#171b1d", "white"),
        linewidth=np.where(top["curated_for_main_figures"], 1.2, 0.6),
        zorder=3,
    )
    ax.axvline(0, color="#30383c", linewidth=0.8)
    labels = [
        textwrap.fill(f"{label}{' *' if curated else ''}", width=46)
        for label, curated in zip(
            top["screen_display_label"], top["curated_for_main_figures"]
        )
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.6)
    ax.set_xlabel("Four-project-balanced flight minus ground score")
    ax.set_title(
        "Top 20 project-balanced absolute effects",
        loc="left",
        fontsize=12,
        weight="bold",
    )
    ax.grid(axis="x", color="#e5e7e8", linewidth=0.6)
    ax.set_axisbelow(True)
    for panel in axes:
        for spine in ("top", "right", "left"):
            panel.spines[spine].set_visible(False)
    fig.suptitle(
        "Skin paired-site sensitivity across all 319 Reactome programs",
        fontsize=15,
        weight="bold",
    )
    fig.text(
        0.5,
        -0.01,
        "MHU-2 and RR-5 dorsal/femoral effects were averaged within project before averaging four mission projects. "
        "Asterisks and black outlines mark the eight programs reviewed for the main skin narrative; all other points were included in this sensitivity analysis.",
        ha="center",
        fontsize=8.8,
        wrap=True,
    )
    save_figure(fig, "figure_s2_skin_project_balance_sensitivity")
    plt.close(fig)


def plot_conceptual_summary() -> None:
    """Draw a non-quantitative synthesis of the complementary pathway hypotheses."""
    import matplotlib.pyplot as plt
    import textwrap
    from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Polygon, Rectangle

    up_color = "#a83b32"
    down_color = "#2d68a7"
    context_color = ROLE_COLORS["context_sensitive"]
    panels = (
        {
            "tissue": "Thymus",
            "summary": "Lower repair and cytoskeletal state is triangulated; lower niche interaction is model-specific.",
            "higher": (),
            "lower": ("DNA repair", "Lymphoid-stromal interactions", "RHOA cytoskeletal cycle"),
            "context": "Innate TLR and ECM were method-supported but failed seed or composition robustness.",
            "icon": "thymus",
            "fill": "#f4e5e9",
        },
        {
            "tissue": "Skin",
            "summary": "Regulatory, repair, regenerative, and barrier-lipid programs are reproducibly lower in flight.",
            "higher": (),
            "lower": (
                "Chromatin regulation",
                "DNA repair",
                "Hedgehog signaling",
                "Sphingolipid metabolism",
                "Cell-cell junction organization*",
            ),
            "context": "*Cell-cell junction organization is internally robust but not supported by ssGSEA; keratinization and gap-junction nodes were seed-sensitive.",
            "icon": "skin",
            "fill": "#f6eee4",
        },
        {
            "tissue": "Liver",
            "summary": "A lower adaptive-immune axis is reproducible beside heterogeneous established metabolism.",
            "higher": (),
            "lower": (
                "MHC II antigen presentation",
                "T-cell receptor signaling",
            ),
            "context": "P450 is project-heterogeneous; Rho-family and ECM directions fail seed or composition robustness.",
            "icon": "liver",
            "fill": "#f2e5e3",
        },
        {
            "tissue": "Soleus",
            "summary": "No reviewed soleus pathway passes all five robustness checks.",
            "higher": (),
            "lower": (),
            "context": "Matrix, repair, metabolic, immune, and contraction programs remain sensitivity-dependent candidates.",
            "icon": "soleus",
            "fill": "#f2e7e4",
        },
    )

    def pathway_block(items: tuple[str, ...], empty_label: str) -> str:
        if not items:
            return empty_label
        return "\n".join(
            textwrap.fill(
                f"\u2022 {item}",
                width=28,
                subsequent_indent="  ",
            )
            for item in items
        )

    fig, axes = plt.subplots(2, 2, figsize=(15.6, 9.4), constrained_layout=True)
    for ax, panel in zip(axes.flat, panels):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.add_patch(
            FancyBboxPatch(
                (0.01, 0.02),
                0.98,
                0.96,
                boxstyle="round,pad=0.012,rounding_size=0.018",
                facecolor="#fafbfb",
                edgecolor="#69757a",
                linewidth=1.0,
            )
        )
        ax.text(0.04, 0.91, panel["tissue"], fontsize=15, weight="bold", va="top")
        ax.text(
            0.04,
            0.82,
            textwrap.fill(panel["summary"], width=68),
            fontsize=9.2,
            color="#30383c",
            va="top",
            linespacing=1.25,
        )

        if panel["icon"] == "thymus":
            ax.add_patch(Ellipse((0.13, 0.48), 0.13, 0.28, angle=12, facecolor=panel["fill"], edgecolor="#914f63", linewidth=1.4))
            ax.add_patch(Ellipse((0.21, 0.48), 0.13, 0.28, angle=-12, facecolor=panel["fill"], edgecolor="#914f63", linewidth=1.4))
            ax.plot([0.17, 0.17], [0.36, 0.60], color="#c58ba0", linewidth=1.0)
        elif panel["icon"] == "skin":
            colors = ("#e8c7ae", "#dba88d", "#f0d7a7")
            for index, color in enumerate(colors):
                ax.add_patch(Rectangle((0.055, 0.37 + index * 0.09), 0.22, 0.085, facecolor=color, edgecolor="white"))
            for x in (0.08, 0.13, 0.18, 0.23):
                ax.add_patch(Circle((x, 0.60), 0.015, facecolor="#9d6a78", edgecolor="none"))
            ax.plot([0.07, 0.25], [0.455, 0.455], color="#7f5966", linewidth=1.2, linestyle="--")
        elif panel["icon"] == "liver":
            ax.add_patch(
                Polygon(
                    [(0.05, 0.44), (0.09, 0.61), (0.22, 0.64), (0.28, 0.54), (0.25, 0.39), (0.12, 0.35)],
                    closed=True,
                    facecolor=panel["fill"],
                    edgecolor="#92463f",
                    linewidth=1.4,
                )
            )
            ax.add_patch(Circle((0.17, 0.50), 0.018, facecolor="#f7f8f8", edgecolor="#92463f"))
            for angle in np.linspace(0, 2 * np.pi, 7)[:-1]:
                ax.plot([0.17, 0.17 + 0.08 * np.cos(angle)], [0.50, 0.50 + 0.08 * np.sin(angle)], color="#c9857e", linewidth=0.8)
        else:
            for index, y in enumerate((0.38, 0.45, 0.52, 0.59)):
                ax.add_patch(FancyBboxPatch((0.055, y), 0.22, 0.055, boxstyle="round,pad=0.004", facecolor="#e6b4aa", edgecolor="#9d5b53", linewidth=0.9))
                for x in np.linspace(0.075, 0.25, 7):
                    ax.plot([x, x + 0.012], [y + 0.012, y + 0.043], color="#f8e5e1", linewidth=0.7)

        ax.text(
            0.33,
            0.68,
            "\u2191  Higher in flight",
            fontsize=9.2,
            weight="bold",
            color=up_color,
            va="top",
        )
        ax.text(
            0.33,
            0.61,
            pathway_block(panel["higher"], "No robust higher\nprogram"),
            fontsize=8.0,
            va="top",
            linespacing=1.28,
        )

        ax.plot([0.64, 0.64], [0.27, 0.69], color="#d7dbdd", linewidth=0.8)
        ax.text(
            0.68,
            0.68,
            "\u2193  Lower in flight",
            fontsize=9.2,
            weight="bold",
            color=down_color,
            va="top",
        )
        ax.text(
            0.68,
            0.61,
            pathway_block(panel["lower"], "No robust lower\nprogram"),
            fontsize=8.0,
            va="top",
            linespacing=1.28,
        )

        ax.add_patch(
            FancyBboxPatch(
                (0.31, 0.055),
                0.65,
                0.15,
                boxstyle="round,pad=0.009,rounding_size=0.012",
                facecolor="#fff8ef",
                edgecolor=context_color,
                linewidth=1.0,
            )
        )
        ax.text(
            0.33,
            0.13,
            textwrap.fill("Context: " + panel["context"], width=76),
            fontsize=7.8,
            color="#5a3a18",
            va="center",
            linespacing=1.2,
        )

    fig.suptitle(
        "Robustness-filtered perspective: reproducible and model-specific pathway layers",
        fontsize=16,
        weight="bold",
    )
    fig.text(
        0.5,
        -0.015,
        "Conceptual synthesis, not a quantitative or mechanistic model. Arrows report pathway-score direction; "
        "orange notes identify programs that fail at least one robustness check.",
        ha="center",
        fontsize=9,
    )
    save_figure(fig, "figure_5_complementary_process_model")
    plt.close(fig)


def plot_generated_process_illustration() -> None:
    """Pair GPT Image 2 tissue artwork with literature-grounded interpretation."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    source_path = FIGURE_DIR / "source/figure_6_gpt_image_2_process_art.png"
    if not source_path.exists():
        raise FileNotFoundError(f"Missing generated illustration source: {source_path}")

    image = plt.imread(source_path)
    height, width = image.shape[:2]
    gutter = 14
    half_height = height // 2
    half_width = width // 2
    crops = (
        image[: half_height - gutter // 2, : half_width - gutter // 2],
        image[: half_height - gutter // 2, half_width + gutter // 2 :],
        image[half_height + gutter // 2 :, : half_width - gutter // 2],
        image[half_height + gutter // 2 :, half_width + gutter // 2 :],
    )
    panels = (
        {
            "tissue": "Thymus",
            "references": "[6-11]",
            "established": "Atrophy, fewer cycling thymocytes,\nand lower adaptive immune activity.",
            "complementary": (
                "Lower repair and cytoskeletal programs are\ntriangulated; lower lymphoid-stromal interaction\n"
                "is an internally robust model-specific signal."
            ),
            "direction": "Repair and cytoskeletal state lower; niche interaction lower",
        },
        {
            "tissue": "Skin",
            "references": "[12-16]",
            "established": "Barrier injury, inflammation, oxidative stress,\nand impaired wound repair.",
            "complementary": (
                "Chromatin, repair, Hedgehog, and barrier-lipid\nprograms are reproducibly lower; cell-junction\n"
                "organization is an added model-specific layer."
            ),
            "direction": "Regulatory, regenerative, barrier-lipid, and junction state lower",
        },
        {
            "tissue": "Liver",
            "references": "[17-23]",
            "established": "Lipid, xenobiotic, insulin-related,\nand mitochondrial dysregulation.",
            "complementary": (
                "MHC II antigen-presentation and T-cell receptor\nprograms are reproducibly lower; cytoskeletal\n"
                "and matrix directions remain model-sensitive."
            ),
            "direction": "Adaptive-immune communication lower",
        },
        {
            "tissue": "Soleus",
            "references": "[24-30]",
            "established": "Atrophy, fiber-type transition, mitochondrial stress,\nand contractile change.",
            "complementary": (
                "No reviewed pathway passes all conventional,\nheld-out, seed, and composition checks; matrix\n"
                "programs remain follow-up candidates."
            ),
            "direction": "No robust complementary pathway direction",
        },
    )

    fig = plt.figure(figsize=(15.2, 9.1), constrained_layout=False)
    outer = fig.add_gridspec(
        2,
        2,
        left=0.035,
        right=0.985,
        bottom=0.085,
        top=0.89,
        wspace=0.08,
        hspace=0.15,
    )
    for index, (crop, panel) in enumerate(zip(crops, panels)):
        row, column = divmod(index, 2)
        inner = outer[row, column].subgridspec(
            1, 2, width_ratios=(1.02, 1.48), wspace=0.055
        )
        image_ax = fig.add_subplot(inner[0, 0])
        text_ax = fig.add_subplot(inner[0, 1])
        image_ax.imshow(crop)
        image_ax.axis("off")
        text_ax.set_xlim(0, 1)
        text_ax.set_ylim(0, 1)
        text_ax.axis("off")

        text_ax.text(0.0, 0.96, panel["tissue"], fontsize=15, weight="bold", va="top")
        text_ax.text(
            0.0, 0.80, f"Established phenotype {panel['references']}", fontsize=9.3, weight="bold",
            color="#4f5a60", va="top",
        )
        text_ax.text(
            0.0, 0.71, panel["established"], fontsize=8.9, color="#30383c",
            va="top", wrap=True, linespacing=1.35,
        )
        text_ax.text(
            0.0, 0.48, "Robustness-filtered expiMap perspective", fontsize=9.3,
            weight="bold", color=ROLE_COLORS["complementary"], va="top",
        )
        text_ax.text(
            0.0, 0.39, panel["complementary"], fontsize=8.9, color="#20282c",
            va="top", wrap=True, linespacing=1.35,
        )
        text_ax.text(
            0.0, 0.08, panel["direction"], fontsize=8.3, weight="bold",
            color="#244f73", va="bottom", wrap=True,
        )
        box = outer[row, column].get_position(fig)
        fig.add_artist(
            Rectangle(
                (box.x0 - 0.008, box.y0 - 0.015),
                box.width + 0.016,
                box.height + 0.03,
                transform=fig.transFigure,
                fill=False,
                edgecolor="#a6afb3",
                linewidth=0.8,
            )
        )

    fig.suptitle(
        "From established organ phenotypes to robustness-filtered pathway hypotheses",
        fontsize=16,
        weight="bold",
        y=0.965,
    )
    fig.text(
        0.5,
        0.025,
        "Conceptual artwork only. Labels summarize literature and robustness-filtered expiMap scores; "
        "the illustration does not depict measured cell states or causal mechanisms.",
        ha="center",
        fontsize=8.8,
        color="#4b555a",
    )
    save_figure(fig, "figure_6_generated_biological_processes")
    plt.close(fig)


def write_readme(summary: pd.DataFrame) -> None:
    lines = [
        "# ASGSR expiMap HVG paper package",
        "",
        "This directory contains the conference abstract, manuscript, source tables, and reproducible figures for the thymus, skin, liver, and soleus HVG expiMap analysis.",
        "",
        "The primary descriptive quantity is the mean of accession-specific flight-minus-ground pathway-score differences after orienting each latent dimension with the sign of its summed decoder weights. Claims are then qualified by conventional enrichment, held-out-project prediction, full-pipeline seed retraining, and atlas-derived broad composition-proxy sensitivity. FDR and the five-check evidence status are complementary summaries, not interchangeable discovery thresholds.",
        "",
        "## Package contents",
        "",
        "- `asgsr_2026_abstract.md`: 2026 ASGSR-formatted abstract.",
        "- `manuscript.md`: full research manuscript with references and figure captions.",
        "- `submission_checklist.md`: current ASGSR requirements and unresolved submission details.",
        "- `supplementary_methods.md`: exact runs, effect definitions, safeguards, and rebuild instructions.",
        "- `aligned_complementary_story.md`: pathway-by-pathway evidence audit and alternative biological narrative.",
        "- `expanded_pathway_family_review.md`: systematic top-decile and stable-extension review after Reactome-family consolidation.",
        "- `reviewer_robustness_analysis.md`: method, held-out-project, and composition-proxy results.",
        "- `reviewer_pathway_evidence.md`: pathway-level five-check interpretation.",
        "- `reviewer_response.md`: reviewer-concern response matrix, evidence locations, and unresolved limitations.",
        "- `tissue_selection_audit.md`: reconstruction of the eight-tissue screen and the implications for focused manuscript scope.",
        "- `requirements.txt`: pinned document-rendering dependencies.",
        "- `figures/`: publication figures in PNG and vector PDF formats.",
        "- `figures/source/`: source artwork and generation prompt for the conceptual biological-process figure.",
        "- `source_data/`: model, accession, pathway, gene-level, and source-verification tables.",
        "- `source_data/table_s9_systematic_pathway_screen.tsv`: all pathways ranked within tissue, including full-model four-project skin sensitivity fields, the top-20 screen, and main-figure review status.",
        "- `source_data/table_s10_expanded_pathway_review.tsv`: 153 pathway candidates with selection reason, sensitivity, gene support, nested-program conflicts, and family assignment.",
        "- `source_data/table_s11_nonredundant_pathway_families.tsv`: 37 reviewed process families and manuscript decisions.",
        "- `source_data/table_s24_pathway_robustness_evidence.tsv`: integrated robustness evidence for all 29 reviewed pathways.",
        "- `visual_audit.md`: rendered-page and standalone-figure quality-control record.",
        "- `source_data/source_verification.tsv`: manuscript DOI records checked against Crossref.",
        "",
        "## Rebuild",
        "",
        "```bash",
        "PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.build_asgsr_paper",
        "PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.reviewer_robustness_analysis",
        "PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.run_asgsr_seed_sensitivity",
        "PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.integrate_reviewer_robustness",
        "```",
        "",
        "## Model scope",
        "",
        "| tissue | ARCHS4 reference | OSDR query | accessions | genes | Reactome programs |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary.itertuples():
        lines.append(
            f"| {row.tissue} | {int(row.reference_samples):,} | {int(row.query_samples):,} | "
            f"{int(row.query_accessions)} | {int(row.genes_after_filter):,} | {int(row.reactome_programs):,} |"
        )
    lines += [
        "",
        "## Interpretation safeguards",
        "",
        "- Higher or lower refers to the expiMap latent pathway score, not uniform expression of every member gene.",
        "- Raw expiMap latent signs are arbitrary. Every reported direction is multiplied by `EXPIMAP.latent_directions(method=\"sum\")`, following the package's upregulation-orientation method.",
        "- This package supersedes older presentation annotations based on raw latent signs.",
        "- Magnitudes are ranked within each tissue because independently trained latent scales are not directly comparable across tissues.",
        "- OSD-289 thymus and OSD-714 soleus confound condition with strain and are excluded in a restricted sensitivity analysis.",
        "- The primary liver query excludes OSD-164 and OSD-168 because they overlap cohorts represented by OSD-47, OSD-48, and OSD-137. The full 12-accession run is retained as sensitivity evidence.",
        "- Skin protocol-context results separate MHU-2 microgravity from onboard artificial 1 g and stratify recovery, terminal collection, duration, and strain where sample names support those contrasts.",
        "- Literature roles and protocol sensitivity are separate annotations. The curated pathway table reports pooled skin effects beside MHU-2 true-microgravity and artificial-1-g effects, and de-duplicated liver effects beside the original 12-accession sensitivity.",
        "- Figure S1 displays the 20 largest absolute active-program effects per tissue without assigning unreviewed terms a biological evidence label.",
        "- Figure S2 applies the paired-site project-balance sensitivity to all 319 skin pathways, not only the eight manually reviewed programs.",
        "- Figure S3 consolidates 153 expanded-review pathways into 37 nonredundant process families; nested opposite-direction programs are treated as unresolved rather than independent mechanisms.",
        "- Figures S4-S7 benchmark conventional scoring, test held-out projects, adjust broad composition proxies, retrain three full pipelines, and integrate the results for every reviewed pathway.",
        "- `Triangulated` requires support from ssGSEA, preranked GSEA, held-out projects, all three seeds, and composition-proxy adjustment. It is a descriptive reproducibility label, not a significance level.",
        "- No reviewed soleus pathway passes all five checks; soleus results are retained as sensitivity-dependent follow-up candidates.",
        "- The 2026 ASGSR abstract limit is 300 words including acknowledgments; figures and tables are not allowed in the abstract submission.",
        "- The 2026 investigator abstract deadline was June 14, 2026, so the abstract is formatted for the meeting but requires confirmation that a submission record already exists or that ASGSR will accept a late/alternate submission.",
        "",
    ]
    (PAPER_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def run() -> None:
    import matplotlib

    matplotlib.use("Agg")
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    gmt = parse_gmt(ROOT / "data/pathways/reactome_current_mouse_ensembl.gmt")

    summaries = []
    covariates = []
    all_gene_results = []
    all_pathways = []
    all_selected = []
    all_accession_effects = []
    skin_context = None
    for config in CONFIGS:
        summaries.append(model_summary(config))
        covariates.append(covariate_audit(config))
        genes = gene_level_results(config)
        all_gene_results.append(genes)
        directions = latent_directions(config)
        if config.tissue == "skin":
            skin_context = skin_protocol_context_results(config, directions)
        pathways, selected = pathway_results(config, genes, gmt, directions)
        all_pathways.append(pathways)
        all_selected.append(selected)
        effects = pd.read_csv(
            config.run_dir / "accession_validation/per_accession_effects.tsv", sep="\t"
        )
        effects = effects.merge(directions, on="term", how="left")
        effects["raw_flight_minus_ground"] = effects["flight_minus_ground"]
        effects["flight_minus_ground"] = (
            effects["flight_minus_ground"] * effects["latent_orientation"]
        )
        effects.loc[
            effects["latent_orientation"].eq(0), "flight_minus_ground"
        ] = np.nan
        effects.insert(0, "tissue", config.tissue)
        effects["condition_strain_confounded"] = effects["id.accession"].astype(str).isin(
            config.confounded_accessions
        )
        all_accession_effects.append(effects)

    summary = pd.DataFrame(summaries)
    covariate = pd.concat(covariates, ignore_index=True)
    gene_results = pd.concat(all_gene_results, ignore_index=True)
    pathway = pd.concat(all_pathways, ignore_index=True)
    selected = pd.concat(all_selected, ignore_index=True)
    accession_effects = pd.concat(all_accession_effects, ignore_index=True)
    if skin_context is None:
        raise RuntimeError("Skin protocol-context results were not generated.")
    selected = add_context_aware_annotations(
        selected, skin_context, accession_effects
    )
    screen = systematic_pathway_screen(pathway, selected, accession_effects)

    summary.to_csv(SOURCE_DIR / "table_s1_model_summary.tsv", sep="\t", index=False)
    covariate.to_csv(SOURCE_DIR / "table_s2_accession_covariate_audit.tsv", sep="\t", index=False)
    pathway.to_csv(SOURCE_DIR / "table_s3_all_pathway_effects.tsv", sep="\t", index=False)
    selected.to_csv(SOURCE_DIR / "table_1_curated_pathway_results.tsv", sep="\t", index=False)
    gene_results.to_csv(SOURCE_DIR / "table_s4_gene_level_results.tsv.gz", sep="\t", index=False, compression="gzip")
    accession_effects.to_csv(SOURCE_DIR / "table_s5_accession_pathway_effects.tsv.gz", sep="\t", index=False, compression="gzip")
    skin_context.to_csv(
        SOURCE_DIR / "table_s8_skin_protocol_context_effects.tsv",
        sep="\t",
        index=False,
    )
    screen.to_csv(
        SOURCE_DIR / "table_s9_systematic_pathway_screen.tsv",
        sep="\t",
        index=False,
    )

    plot_workflow(summary)
    plot_pathway_shifts(selected, accession_effects)
    plot_evidence_map(selected)
    plot_sensitivity(selected)
    plot_conceptual_summary()
    plot_generated_process_illustration()
    plot_skin_protocol_context(skin_context)
    plot_broad_pathway_screen(screen)
    plot_skin_project_balance(screen)
    from expiMap_scarches.nasa_mouse_expimap.review_expanded_pathway_screen import (
        run as run_expanded_review,
    )

    run_expanded_review()
    write_readme(summary)
    print(
        json.dumps(
            {
                "paper_dir": str(PAPER_DIR.relative_to(ROOT)),
                "models": len(summary),
                "all_pathways": len(pathway),
                "curated_pathways": len(selected),
                "main_figures": 7,
                "supplementary_figures": 3,
            },
            indent=2,
        )
    )


def main() -> None:
    run()


if __name__ == "__main__":
    main()
