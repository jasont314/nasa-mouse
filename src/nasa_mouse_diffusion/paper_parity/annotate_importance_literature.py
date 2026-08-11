"""Annotate matched gene and grouped pathway importance results with literature."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

from .annotate_promoted_gene_literature import (
    ALLOWED_CLASSIFICATIONS,
    INTERPRETIVE_ROLE_BY_CLASSIFICATION,
    SOURCES as CONSENSUS_SOURCES,
    LiteratureSource,
)


ROOT = Path(__file__).resolve().parents[3]
SOURCE_DIR = ROOT / "paper" / "synthetic_guided_spaceflight" / "source_data"
GROUPED_DIR = (
    ROOT
    / "outputs"
    / "generative_benchmark"
    / "analyses"
    / "grouped_pathway_importance_osdr_disjoint_v1"
)
MATCHED_INPUT = SOURCE_DIR / "table_s19_matched_all_gene_candidates.tsv"
CONSENSUS_ANNOTATION_INPUT = (
    SOURCE_DIR / "table_s16_promoted_gene_literature_annotations.tsv"
)
GROUPED_INPUT = GROUPED_DIR / "eligible_synthetic_pathways.tsv.gz"
NONREDUNDANT_INPUT = GROUPED_DIR / "top_nonredundant_pathways.tsv"
GENE_OUTPUT = SOURCE_DIR / "table_s22_matched_gene_literature_annotations.tsv"
PATHWAY_OUTPUT = SOURCE_DIR / "table_s23_grouped_pathway_literature_annotations.tsv"
SOURCE_OUTPUT = SOURCE_DIR / "table_s24_importance_literature_sources.tsv"
SEARCH_DATE = "2026-08-05"

EXPECTED_GENE_COUNTS = {
    "aligning": 9,
    "complementary": 9,
    "ambiguous": 1,
    "unmatched": 2,
}
EXPECTED_PATHWAY_COUNTS = {
    "aligning": 6,
    "complementary": 2,
    "ambiguous": 2,
    "unmatched": 0,
}


@dataclass(frozen=True)
class ImportanceAnnotation:
    tissue: str
    feature: str
    literature_classification: str
    evidence_scope: str
    evidence_relationship: str
    source_ids: tuple[str, ...]
    literature_summary: str
    interpretation: str


NEW_SOURCES = (
    LiteratureSource(
        "beheshti_2019_liver",
        "Beheshti A, Chakravarty K, Fogle H, et al. Multi-omics analysis of multiple missions to space reveal a theme of lipid dysregulation in mouse liver. Scientific Reports. 2019;9:19195.",
        2019,
        "10.1038/s41598-019-55869-2",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC6915713/",
        "Multi-mission mouse liver evidence for lipid, insulin, and growth-factor signaling changes after spaceflight.",
        "Uses public mouse spaceflight cohorts; treated as related process evidence rather than independent validation.",
    ),
    LiteratureSource(
        "blaber_2017_liver_proteostasis",
        "Blaber EA, Pecaut MJ, Jonscher KR. Spaceflight activates autophagy programs and the proteasome in mouse liver. International Journal of Molecular Sciences. 2017;18:2062.",
        2017,
        "10.3390/ijms18102062",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC5666744/",
        "Mouse liver evidence for altered autophagy, proteasome activity, and protein homeostasis after spaceflight.",
        "Published mouse flight cohort; process context rather than an exact PPIC replication.",
    ),
    LiteratureSource(
        "vitry_2022_liver_muscle",
        "Vitry G, Finch R, McStay G, et al. Muscle atrophy phenotype gene expression during spaceflight is linked to a metabolic crosstalk in both the liver and the muscle in mice. iScience. 2022;25:105213.",
        2022,
        "10.1016/j.isci.2022.105213",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC9576569/",
        "Reanalysis of mouse flight liver and muscle reporting lower RNA-polymerase and protein-metabolism pathways.",
        "Uses public spaceflight data that may overlap the present aggregate; process-level evidence only.",
    ),
    LiteratureSource(
        "da_silveira_2020_mitochondrial",
        "da Silveira WA, Fazelinia H, Rosenthal SB, et al. Comprehensive multi-omics analysis reveals mitochondrial stress as a central biological hub for spaceflight impact. Cell. 2020;183:1185-1201.e20.",
        2020,
        "10.1016/j.cell.2020.11.002",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC7870178/",
        "Cross-tissue spaceflight multi-omics evidence that includes reduced interferon-related signatures in mouse liver.",
        "Uses public spaceflight cohorts; broad immune-process context rather than an exact H2-DMa replication.",
    ),
    LiteratureSource(
        "luo_2018_grb10",
        "Luo L, Jiang W, Liu H, et al. De-silencing Grb10 contributes to acute ER stress-induced steatosis in mouse liver. Journal of Molecular Endocrinology. 2018;60:285-297.",
        2018,
        "10.1530/JME-18-0018",
        "https://pubmed.ncbi.nlm.nih.gov/29555819/",
        "Mouse liver mechanism connecting GRB10 with insulin or IGF signaling, ER stress, and steatosis.",
        "Mechanistic liver context only; not a spaceflight replication.",
    ),
    LiteratureSource(
        "yang_2021_ppic",
        "Yang X, Shu B, Zhou Y, Li Z, He C. Ppic modulates CCl4-induced liver fibrosis and TGF-beta-caused mouse hepatic stellate cell activation and is regulated by miR-137-3p. Toxicology Letters. 2021;350:52-61.",
        2021,
        "10.1016/j.toxlet.2021.06.021",
        "https://doi.org/10.1016/j.toxlet.2021.06.021",
        "Mouse liver injury evidence connecting PPIC with stellate-cell activation and fibrotic remodeling.",
        "Mechanistic liver-injury context only; not a spaceflight replication.",
    ),
    LiteratureSource(
        "felix_2000_h2dma",
        "Felix NJ, Brickey WJ, Griffiths R, et al. H2-DMalpha(-/-) mice show the importance of major histocompatibility complex-bound peptide in cardiac allograft rejection. Journal of Experimental Medicine. 2000;192:31-40.",
        2000,
        "10.1084/jem.192.1.31",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC1887714/",
        "Primary mouse evidence that H2-DMalpha is required for effective MHC class II peptide presentation.",
        "Mechanistic antigen-presentation context only; not a liver or spaceflight replication.",
    ),
    LiteratureSource(
        "rusnac_2018_klhdc2",
        "Rusnac DV, Lin HC, Canzani D, et al. Recognition of the diglycine C-end degron by CRL2-KLHDC2 ubiquitin ligase. Molecular Cell. 2018;72:813-822.e4.",
        2018,
        "10.1016/j.molcel.2018.10.021",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC6294321/",
        "Structural evidence that KLHDC2 recognizes C-end degrons in ubiquitin-mediated protein quality control.",
        "Mechanistic protein-homeostasis context only; not a thymus or spaceflight replication.",
    ),
    LiteratureSource(
        "tokoro_2001_tspan3",
        "Tokoro Y, Shibuya K, Osawa M, et al. Molecular cloning and characterization of mouse Tspan-3, a novel member of the tetraspanin superfamily, expressed on resting dendritic cells. Biochemical and Biophysical Research Communications. 2001;288:178-183.",
        2001,
        "10.1006/bbrc.2001.5742",
        "https://pubmed.ncbi.nlm.nih.gov/11594770/",
        "Primary evidence that mouse TSPAN3 is expressed on resting dendritic cells.",
        "Immune-cell identity context only; not a thymus or spaceflight replication.",
    ),
    LiteratureSource(
        "knosp_2011_socs2",
        "Knosp CA, Carroll HP, Elliott J, et al. SOCS2 regulates T helper type 2 differentiation and the generation of type 2 allergic responses. Journal of Experimental Medicine. 2011;208:1523-1531.",
        2011,
        "10.1084/jem.20101167",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC3135359/",
        "Primary evidence that SOCS2 regulates cytokine signaling and CD4 T-cell differentiation.",
        "Mechanistic T-cell context only; not a thymus or spaceflight replication.",
    ),
    LiteratureSource(
        "kumari_2021_skin_necroptosis",
        "Kumari S, Van TM, Preukschat D, et al. NF-kappaB inhibition in keratinocytes causes RIPK1-mediated necroptosis and skin inflammation. Life Science Alliance. 2021;4:e202000956.",
        2021,
        "10.26508/lsa.202000956",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8091601/",
        "Primary skin evidence linking RIPK1-mediated necroptosis to keratinocyte death and inflammation.",
        "Mechanistic skin context only; not a spaceflight replication.",
    ),
    LiteratureSource(
        "gridley_2013_thymus_spleen",
        "Gridley DS, Mao XW, Stodieck LS, et al. Changes in mouse thymus and spleen after return from the space mission STS-135. PLOS ONE. 2013;8:e75097.",
        2013,
        "10.1371/journal.pone.0075097",
        "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075097",
        "Mouse flight evidence for lower thymic G1/S regulators, altered immune signaling, and nominal growth-factor receptor changes.",
        "Independent shuttle mission and targeted expression platform relative to most present RNA-seq cohorts.",
    ),
    LiteratureSource(
        "hughes_fulford_2015_tcell",
        "Hughes-Fulford M, Chang TT, Martinez EM, Li CF. Spaceflight alters expression of microRNA during T-cell activation. FASEB Journal. 2015;29:4893-4900.",
        2015,
        "10.1096/fj.15-277392",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC4653058/",
        "Primary spaceflight evidence for impaired early T-cell activation involving AP-1 and NF-kappaB regulation.",
        "Human T-cell mechanism and different assay context; not a mouse spleen replication.",
    ),
)


GENE_ANNOTATIONS = (
    ImportanceAnnotation(
        "liver",
        "Grb10",
        "complementary",
        "same_tissue_spaceflight_process_and_gene_mechanism",
        "Public mouse flight data provide process context; the GRB10 study is mechanistic and does not reproduce the observed flight direction.",
        ("beheshti_2019_liver", "luo_2018_grb10"),
        "Mouse spaceflight liver studies report lipid and insulin-signaling dysregulation. GRB10 regulates insulin or IGF signaling and can contribute to ER-stress steatosis in mouse liver, but no prior directional liver-flight Grb10 result was found.",
        "The flight-lower association adds a growth-factor-control candidate to hepatic metabolic remodeling without establishing GRB10 as a flight mechanism.",
    ),
    ImportanceAnnotation(
        "liver",
        "Ppic",
        "complementary",
        "same_tissue_spaceflight_process_and_gene_mechanism",
        "The flight study supports hepatic proteostasis changes; the PPIC study is a liver-injury mechanism rather than a flight replication.",
        ("blaber_2017_liver_proteostasis", "yang_2021_ppic"),
        "Mouse spaceflight liver studies report altered autophagy, proteasome activity, and protein homeostasis. PPIC promotes stellate-cell activation in mouse liver injury, but no prior directional flight-liver Ppic result was found.",
        "The flight-lower association links proteostasis with a candidate hepatic remodeling mechanism that requires direct validation.",
    ),
    ImportanceAnnotation(
        "liver",
        "H2-DMa",
        "complementary",
        "related_spaceflight_immune_process_and_gene_mechanism",
        "The spaceflight source provides liver immune-process context; H2-DMalpha function is mechanistic and not a directional flight match.",
        ("da_silveira_2020_mitochondrial", "felix_2000_h2dma"),
        "Spaceflight multi-omics reports reduced interferon-related signatures in mouse liver. H2-DMalpha is required for MHC class II peptide loading, but no direct liver-flight H2-DMa result was found.",
        "The flight-lower association adds an antigen-presentation or immune-composition candidate to the liver response.",
    ),
    ImportanceAnnotation(
        "liver",
        "Gtf2a2",
        "aligning",
        "same_tissue_process_same_direction",
        "Public mouse flight data may overlap the present aggregate; agreement is at the transcription-process level rather than the exact gene.",
        ("vitry_2022_liver_muscle",),
        "A reanalysis of mouse flight liver reported lower RNA-polymerase and protein-metabolism pathways. GTF2A2 is a TFIIA subunit used in RNA polymerase II transcription initiation.",
        "Lower Gtf2a2 agrees with previously reported suppression of a broad hepatic transcriptional process, not an exact published gene result.",
    ),
    ImportanceAnnotation(
        "thymus",
        "Klhdc2",
        "unmatched",
        "targeted_search_no_thymus_spaceflight_match",
        "The thymus source defines the searched flight context; KLHDC2 evidence is mechanistic only.",
        ("horie_2019_thymus", "rusnac_2018_klhdc2"),
        "KLHDC2 recognizes C-end degrons in ubiquitin-mediated protein quality control. The targeted search found no prior thymus-spaceflight association for Klhdc2.",
        "A literature-unmatched thymus candidate that should be kept separate from the supported mitotic program.",
    ),
    ImportanceAnnotation(
        "thymus",
        "E2f2",
        "aligning",
        "same_tissue_process_same_direction",
        "One source may overlap OSDR cohorts and one is an independent shuttle mission; neither reports the exact E2f2 result.",
        ("horie_2019_thymus", "gridley_2013_thymus_spleen"),
        "Mouse flight studies report lower thymic cell-cycle activity and lower E2F-related G1/S regulation. The current E2f2 direction agrees at the same-tissue process level.",
        "Flight-lower E2f2 extends the replicated thymic cell-cycle response to another E2F-family regulator.",
    ),
    ImportanceAnnotation(
        "thymus",
        "Plscr1",
        "complementary",
        "related_thymus_immune_process_and_gene_mechanism",
        "The thymus study may overlap OSDR cohorts; PLSCR1 evidence is mechanistic and not a directional thymus-flight replication.",
        ("horie_2019_thymus", "dong_2004_plscr1"),
        "Mouse flight thymus shows immune and composition changes. PLSCR1 is interferon inducible and can amplify interferon-stimulated transcription, but no prior thymus-flight Plscr1 direction was found.",
        "The flight-higher association adds an interferon-response candidate distinct from the flight-lower mitotic panel.",
    ),
    ImportanceAnnotation(
        "thymus",
        "Cdc20",
        "aligning",
        "same_tissue_process_same_direction",
        "One source may overlap OSDR cohorts and one is an independent shuttle mission; agreement is at the APC/C and cell-cycle process level.",
        ("horie_2019_thymus", "gridley_2013_thymus_spleen"),
        "Mouse flight studies report reduced thymic cell-cycle progression. CDC20 activates APC/C-mediated mitotic protein turnover, matching the direction of the broader flight-lower mitotic program.",
        "Flight-lower Cdc20 reinforces a same-tissue cell-cycle process reported in prior missions.",
    ),
    ImportanceAnnotation(
        "thymus",
        "Tspan3",
        "complementary",
        "related_thymus_composition_and_gene_mechanism",
        "The thymus study provides flight context; TSPAN3 evidence is immune-cell identity context only.",
        ("horie_2019_thymus", "tokoro_2001_tspan3"),
        "Spaceflight changes thymic cellularity and immune state. Mouse TSPAN3 is expressed on resting dendritic cells, but no directional thymus-flight Tspan3 result was found.",
        "The flight-higher association may reflect an immune-state or cell-composition shift and needs cell-resolved follow-up.",
    ),
    ImportanceAnnotation(
        "thymus",
        "Socs2",
        "complementary",
        "related_thymus_immune_process_and_gene_mechanism",
        "The thymus study provides flight context; SOCS2 evidence is mechanistic T-cell context only.",
        ("horie_2019_thymus", "knosp_2011_socs2"),
        "Spaceflight alters thymic immune state. SOCS2 regulates cytokine signaling and CD4 T-cell differentiation, but no prior directional thymus-flight Socs2 result was found.",
        "The flight-higher association adds a cytokine and T-cell-state hypothesis alongside the mitotic response.",
    ),
)


PATHWAY_ANNOTATIONS = (
    ImportanceAnnotation(
        "skin",
        "R-MMU-5213460_RIPK1_MEDIATED_REGULATED_NECROSIS",
        "complementary",
        "same_tissue_spaceflight_process_and_pathway_mechanism",
        "The skin study uses public OSDR cohorts; the necroptosis study is a skin mechanism rather than a flight replication.",
        ("cope_2024_skin", "kumari_2021_skin_necroptosis"),
        "Spaceflight skin studies report immune, barrier, mitochondrial, and damage responses. RIPK1-mediated necroptosis can drive keratinocyte death and skin inflammation, but direct spaceflight-necroptosis evidence was not found.",
        "The flight-higher group provides a testable regulated-cell-death mechanism for the broader skin response.",
    ),
    ImportanceAnnotation(
        "skin",
        "R-MMU-5675482_REGULATION_OF_NECROPTOTIC_CELL_DEATH",
        "complementary",
        "same_tissue_spaceflight_process_and_pathway_mechanism",
        "The skin study uses public OSDR cohorts; the necroptosis study is a skin mechanism rather than a flight replication.",
        ("cope_2024_skin", "kumari_2021_skin_necroptosis"),
        "Spaceflight skin studies report immune, barrier, mitochondrial, and damage responses. RIPK1-dependent necroptosis can produce inflammatory keratinocyte loss, but direct spaceflight-necroptosis evidence was not found.",
        "This overlapping flight-higher pathway adds regulated necroptotic death as a testable skin hypothesis.",
    ),
    ImportanceAnnotation(
        "spleen",
        "R-MMU-450341_ACTIVATION_OF_THE_AP_1_FAMILY_OF_TRANSCRIPTION_FACTORS",
        "ambiguous",
        "mixed_same_system_process_evidence",
        "Prior evidence comes from a mouse shuttle mission and human T-cell activation; neither is an exact whole-spleen RNA-seq replication.",
        ("gridley_2013_thymus_spleen", "hughes_fulford_2015_tcell"),
        "The current whole-spleen group is flight higher, whereas prior spaceflight work often reports impaired T-cell activation and lower AP-1-related signaling. Mission, cell composition, and activation context differ.",
        "The AP-1 result is biologically relevant but directionally context dependent and should not be described as straightforward replication.",
    ),
    ImportanceAnnotation(
        "thymus",
        "R-MMU-8847993_ERBB2_ACTIVATES_PTK6_SIGNALING",
        "ambiguous",
        "mixed_same_tissue_growth_signaling_evidence",
        "One source is an independent shuttle mission and one may overlap OSDR cohorts; evidence concerns related growth signaling rather than this exact Reactome term.",
        ("gridley_2013_thymus_spleen", "horie_2019_thymus"),
        "STS-135 thymus showed a nominal flight-higher Egfr signal, while longer-duration mouse studies emphasized reduced thymic proliferation. Evidence for flight-higher ERBB-family signaling is therefore limited and context dependent.",
        "The pathway is a growth-signaling hypothesis, not a settled counterpart to the flight-lower mitotic program.",
    ),
    *(
        ImportanceAnnotation(
            "thymus",
            term,
            "aligning",
            "same_tissue_process_same_direction",
            "One source may overlap public OSDR cohorts and one is an independent shuttle mission; agreement is at the pathway level.",
            ("horie_2019_thymus", "gridley_2013_thymus_spleen"),
            "Independent mouse flight studies report lower thymic cell-cycle progression and proliferation. The current APC/C, chromosome-condensation, or G2/M group is flight lower in the same tissue.",
            "This pathway reinforces the established flight-lower thymic mitotic program while identifying a more specific regulatory module.",
        )
        for term in (
            "R-MMU-174154_APC_C_CDC20_MEDIATED_DEGRADATION_OF_SECURIN",
            "R-MMU-174184_CDC20_PHOSPHO_APC_C_MEDIATED_DEGRADATION_OF_CYCLIN_A",
            "R-MMU-176408_REGULATION_OF_APC_C_ACTIVATORS_BETWEEN_G1_S_AND_EARLY_ANAPHASE",
            "R-MMU-179419_APC_CDC20_MEDIATED_DEGRADATION_OF_CELL_CYCLE_PROTEINS_PRIOR_TO_SATISFATION_OF_THE_CELL_CYCLE_CHECKPOINT",
            "R-MMU-2299718_CONDENSATION_OF_PROPHASE_CHROMOSOMES",
            "R-MMU-69478_G2_M_DNA_REPLICATION_CHECKPOINT",
        )
    ),
)


def _only(group: pd.DataFrame, column: str):
    values = group[column].drop_duplicates()
    if len(values) != 1:
        raise ValueError(
            f"Expected one {column} for {group[['tissue']].iloc[0].to_dict()}; "
            f"found {values.tolist()}"
        )
    return values.iloc[0]


def _arm_values(group: pd.DataFrame, value_column: str) -> str:
    return ";".join(
        f"{row.arm}:{getattr(row, value_column)}"
        for row in group.itertuples(index=False)
    )


def _collapse_matched() -> pd.DataFrame:
    matched = pd.read_csv(MATCHED_INPUT, sep="\t")
    keys = ["analysis_scope", "tissue", "gene", "symbol"]
    if matched.duplicated(keys + ["arm"]).any():
        raise ValueError("Matched input contains duplicate tissue-gene-arm rows")

    rows = []
    for order, (_, group) in enumerate(matched.groupby(keys, sort=False)):
        rows.append(
            {
                "analysis_scope": _only(group, "analysis_scope"),
                "scope": _only(group, "scope"),
                "tissue": _only(group, "tissue"),
                "gene": _only(group, "gene"),
                "symbol": _only(group, "symbol"),
                "flt_gc_direction": _only(group, "flt_gc_direction"),
                "n_accessions": _only(group, "n_accessions"),
                "meta_effect": _only(group, "meta_effect"),
                "meta_fdr": _only(group, "meta_fdr"),
                "supporting_arms": ";".join(group["arm"]),
                "importance_patterns": _arm_values(group, "pattern"),
                "matched_importance_interpretations": _arm_values(
                    group, "matched_importance_interpretation"
                ),
                "matched_statuses": _arm_values(group, "matched_status"),
                "real_only_auroc_permutation_loss_mean": group[
                    "real_only_permutation_roc_auc_mean"
                ].mean(),
                "synthetic_arm_auroc_permutation_loss_mean": group[
                    "arm_real_permutation_roc_auc_mean"
                ].mean(),
                "synthetic_arm_auroc_permutation_loss_max": group[
                    "arm_real_permutation_roc_auc_mean"
                ].max(),
                "synthetic_arm_shap_flt_minus_gc_mean": group[
                    "arm_real_linear_shap_flight_minus_ground"
                ].mean(),
                "joint_utility_nonworse_all_arms": bool(
                    group["joint_mean_all_metrics_nonworse"].all()
                ),
                "_source_order": order,
            }
        )
    collapsed = pd.DataFrame(rows)
    if len(collapsed) != 21:
        raise ValueError(f"Expected 21 matched associations, found {len(collapsed)}")
    return collapsed


def _annotation_frame(annotations: tuple[ImportanceAnnotation, ...]) -> pd.DataFrame:
    rows = []
    for annotation in annotations:
        row = asdict(annotation)
        row["source_ids"] = ";".join(annotation.source_ids)
        rows.append(row)
    return pd.DataFrame(rows)


def _build_gene_table() -> pd.DataFrame:
    matched = _collapse_matched()
    consensus = pd.read_csv(CONSENSUS_ANNOTATION_INPUT, sep="\t")
    consensus = consensus[
        [
            "tissue",
            "symbol",
            "flt_gc_direction",
            "literature_classification",
            "evidence_scope",
            "evidence_relationship",
            "source_ids",
            "literature_summary",
            "interpretation",
        ]
    ].copy()
    consensus = consensus.rename(
        columns={"flt_gc_direction": "annotation_flt_gc_direction"}
    )
    consensus["annotation_origin"] = "reused_consensus_annotation"
    consensus = consensus.rename(columns={"symbol": "feature"})

    curated = _annotation_frame(GENE_ANNOTATIONS)
    curated["annotation_origin"] = "matched_specific_review"
    if curated.duplicated(["tissue", "feature"]).any():
        raise ValueError("Matched-specific literature annotations are duplicated")

    available = pd.concat([consensus, curated], ignore_index=True)
    matched_keys = set(zip(matched["tissue"], matched["symbol"]))
    annotations = available[
        available.apply(lambda row: (row["tissue"], row["feature"]) in matched_keys, axis=1)
    ].copy()
    if annotations.duplicated(["tissue", "feature"]).any():
        duplicates = annotations.loc[
            annotations.duplicated(["tissue", "feature"], keep=False),
            ["tissue", "feature"],
        ]
        raise ValueError(f"Matched annotations overlap: {duplicates.to_dict('records')}")

    observed_keys = set(zip(annotations["tissue"], annotations["feature"]))
    if matched_keys != observed_keys:
        raise ValueError(
            "Matched annotation coverage mismatch; "
            f"missing={sorted(matched_keys - observed_keys)}, "
            f"extra={sorted(observed_keys - matched_keys)}"
        )

    reused = matched.merge(
        consensus[["tissue", "feature", "annotation_flt_gc_direction"]],
        left_on=["tissue", "symbol"],
        right_on=["tissue", "feature"],
        how="inner",
        validate="one_to_one",
    )
    direction_mismatch = reused[
        reused["flt_gc_direction"] != reused["annotation_flt_gc_direction"]
    ]
    if not direction_mismatch.empty:
        raise ValueError(
            "Reused consensus annotations have a different FLT/GC direction: "
            f"{direction_mismatch[['tissue', 'symbol']].to_dict('records')}"
        )

    annotations = annotations.rename(columns={"feature": "symbol"})
    annotations = annotations.drop(columns="annotation_flt_gc_direction")
    merged = matched.merge(
        annotations,
        on=["tissue", "symbol"],
        how="left",
        validate="one_to_one",
    )
    merged.insert(
        merged.columns.get_loc("literature_classification") + 1,
        "interpretive_role",
        merged["literature_classification"].map(
            INTERPRETIVE_ROLE_BY_CLASSIFICATION
        ),
    )
    merged["literature_search_date"] = SEARCH_DATE
    merged = merged.sort_values("_source_order").drop(columns="_source_order")
    _validate_class_counts(merged, EXPECTED_GENE_COUNTS, "matched gene")
    return merged


def _collapse_grouped() -> pd.DataFrame:
    grouped = pd.read_csv(GROUPED_INPUT, sep="\t")
    nonredundant = pd.read_csv(NONREDUNDANT_INPUT, sep="\t")
    keys = ["scope", "tissue", "term"]
    if grouped.duplicated(keys + ["arm"]).any():
        raise ValueError("Grouped input contains duplicate tissue-pathway-arm rows")
    nonredundant_keys = set(zip(nonredundant["tissue"], nonredundant["term"]))

    rows = []
    for order, (_, group) in enumerate(grouped.groupby(keys, sort=False)):
        tissue = _only(group, "tissue")
        term = _only(group, "term")
        rows.append(
            {
                "scope": _only(group, "scope"),
                "tissue": tissue,
                "term": term,
                "description": _only(group, "description"),
                "url": _only(group, "url"),
                "genes": _only(group, "genes"),
                "symbols": _only(group, "symbols"),
                "n_accessions": _only(group, "n_accessions"),
                "meta_effect": _only(group, "meta_effect"),
                "meta_fdr": _only(group, "meta_fdr"),
                "flt_gc_direction": _only(group, "flt_gc_direction"),
                "supporting_arms": ";".join(group["arm"]),
                "group_importance_patterns": _arm_values(
                    group, "group_importance_pattern"
                ),
                "real_only_group_auroc_permutation_loss_mean": group[
                    "real_only_permutation_roc_auc_mean"
                ].mean(),
                "synthetic_arm_group_auroc_permutation_loss_mean": group[
                    "arm_permutation_roc_auc_mean"
                ].mean(),
                "synthetic_arm_group_auroc_permutation_loss_max": group[
                    "arm_permutation_roc_auc_mean"
                ].max(),
                "synthetic_arm_group_shap_flt_minus_gc_mean": group[
                    "arm_group_shap_flight_minus_ground"
                ].mean(),
                "joint_utility_nonworse_all_arms": bool(
                    group["joint_mean_all_metrics_nonworse"].all()
                ),
                "group_shap_supports_flt_gc_all_arms": bool(
                    group["group_shap_supports_flt_gc"].all()
                ),
                "is_nonredundant_pathway": (tissue, term) in nonredundant_keys,
                "_source_order": order,
            }
        )
    collapsed = pd.DataFrame(rows)
    if len(collapsed) != 10:
        raise ValueError(f"Expected 10 grouped pathway associations, found {len(collapsed)}")
    if collapsed["is_nonredundant_pathway"].sum() != 9:
        raise ValueError("Expected nine nonredundant grouped pathway associations")
    return collapsed


def _build_pathway_table() -> pd.DataFrame:
    grouped = _collapse_grouped()
    annotations = _annotation_frame(PATHWAY_ANNOTATIONS).rename(
        columns={"feature": "term"}
    )
    annotations["annotation_origin"] = "grouped_pathway_review"
    keys = ["tissue", "term"]
    if annotations.duplicated(keys).any():
        raise ValueError("Grouped pathway literature annotations are duplicated")
    expected_keys = set(map(tuple, grouped[keys].itertuples(index=False, name=None)))
    observed_keys = set(
        map(tuple, annotations[keys].itertuples(index=False, name=None))
    )
    if expected_keys != observed_keys:
        raise ValueError(
            "Grouped annotation coverage mismatch; "
            f"missing={sorted(expected_keys - observed_keys)}, "
            f"extra={sorted(observed_keys - expected_keys)}"
        )

    merged = grouped.merge(annotations, on=keys, how="left", validate="one_to_one")
    merged.insert(
        merged.columns.get_loc("literature_classification") + 1,
        "interpretive_role",
        merged["literature_classification"].map(
            INTERPRETIVE_ROLE_BY_CLASSIFICATION
        ),
    )
    merged["literature_search_date"] = SEARCH_DATE
    merged = merged.sort_values("_source_order").drop(columns="_source_order")
    _validate_class_counts(merged, EXPECTED_PATHWAY_COUNTS, "grouped pathway")
    return merged


def _validate_class_counts(
    frame: pd.DataFrame, expected: dict[str, int], label: str
) -> None:
    observed_classes = set(frame["literature_classification"])
    if not observed_classes <= ALLOWED_CLASSIFICATIONS:
        raise ValueError(
            f"Unexpected {label} classifications: "
            f"{sorted(observed_classes - ALLOWED_CLASSIFICATIONS)}"
        )
    counts = (
        frame["literature_classification"]
        .value_counts()
        .reindex(sorted(ALLOWED_CLASSIFICATIONS), fill_value=0)
        .to_dict()
    )
    expected_sorted = {key: expected[key] for key in sorted(expected)}
    if counts != expected_sorted:
        raise ValueError(f"Unexpected {label} classification counts: {counts}")


def _build_source_table(
    genes: pd.DataFrame, pathways: pd.DataFrame
) -> pd.DataFrame:
    all_sources = (*CONSENSUS_SOURCES, *NEW_SOURCES)
    source_by_id: dict[str, LiteratureSource] = {}
    for source in all_sources:
        if source.source_id in source_by_id:
            raise ValueError(f"Duplicate literature source ID: {source.source_id}")
        source_by_id[source.source_id] = source

    usage: dict[str, set[str]] = {}
    for label, frame in (
        ("matched_gene", genes),
        ("grouped_pathway", pathways),
    ):
        for source_ids in frame["source_ids"]:
            for source_id in source_ids.split(";"):
                usage.setdefault(source_id, set()).add(label)
    missing = set(usage) - set(source_by_id)
    if missing:
        raise ValueError(f"Unresolved literature source IDs: {sorted(missing)}")

    rows = []
    for source_id in sorted(usage):
        row = asdict(source_by_id[source_id])
        row["annotation_sets"] = ";".join(sorted(usage[source_id]))
        row["literature_search_date"] = SEARCH_DATE
        rows.append(row)
    return pd.DataFrame(rows)


def build_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes = _build_gene_table()
    pathways = _build_pathway_table()
    sources = _build_source_table(genes, pathways)
    return genes, pathways, sources


def write_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes, pathways, sources = build_tables()
    genes.to_csv(GENE_OUTPUT, sep="\t", index=False)
    pathways.to_csv(PATHWAY_OUTPUT, sep="\t", index=False)
    sources.to_csv(SOURCE_OUTPUT, sep="\t", index=False)
    return genes, pathways, sources


def check_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expected = build_tables()
    observed = (
        pd.read_csv(GENE_OUTPUT, sep="\t"),
        pd.read_csv(PATHWAY_OUTPUT, sep="\t"),
        pd.read_csv(SOURCE_OUTPUT, sep="\t"),
    )
    for observed_frame, expected_frame in zip(observed, expected, strict=True):
        assert_frame_equal(observed_frame, expected_frame, check_dtype=False)
    return observed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify that the committed tables match the curated annotations.",
    )
    args = parser.parse_args()

    genes, pathways, sources = check_tables() if args.check else write_tables()
    print(
        f"Validated {len(genes)} matched genes and {len(pathways)} grouped "
        f"pathways against {len(sources)} sources; "
        f"gene classes={genes['literature_classification'].value_counts().to_dict()}; "
        f"pathway classes={pathways['literature_classification'].value_counts().to_dict()}"
    )


if __name__ == "__main__":
    main()
