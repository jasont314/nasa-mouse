"""Audit all aligned and complementary programs used to construct the paper story."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
PAPER_DIR = ROOT / "paper/asgsr_expimap_hvg"
SOURCE_DIR = PAPER_DIR / "source_data"

RUNS = {
    "thymus": ROOT
    / "outputs/expimap_archs4_reference_osdr_query_thymus/tutorial_hvg_2000/query_nb_250epoch_seed2020",
    "skin": ROOT
    / "outputs/expimap_archs4_reference_osdr_query_skin/tutorial_hvg_2000/query_nb_250epoch_seed2020",
    "liver": ROOT
    / "outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/query_nb_250epoch_seed2020_primary_deduplicated",
    "soleus": ROOT
    / "outputs/expimap_muscle_targeted_combined_min8/tutorial_hvg_soleus_2000/query_nb_250epoch_seed2020",
}

# Every program that was aligned or complementary before this deeper audit is retained here,
# including terms subsequently downgraded to context-sensitive.
REVIEW = {
    ("thymus", "Mitotic cell cycle"): (
        "direct_tissue_spaceflight",
        "core",
        "adaptive_proliferative_niche",
        "Lower proliferation is directly established in two flight cohorts and is the anchor phenotype.",
        "Horie 2019; Gridley 2013",
        "A lower bulk score can also reflect fewer proliferating thymocytes.",
    ),
    ("thymus", "DNA repair"): (
        "process_supported",
        "core",
        "adaptive_proliferative_niche",
        "Spaceflight thymus shows DNA fragmentation and stress, but lower repair capacity was not directly measured.",
        "Gridley 2013; Horie 2019",
        "Interpret as a reduced repair-associated state, not proof of defective repair kinetics.",
    ),
    ("thymus", "T-cell receptor signaling"): (
        "direct_process_spaceflight",
        "core",
        "adaptive_proliferative_niche",
        "Altered thymic T-cell signaling and reduced thymic output support the lower adaptive score.",
        "Lebsack 2010; Horie 2019",
        "Bulk abundance of developing T cells can drive this score.",
    ),
    ("thymus", "Innate TLR signaling"): (
        "indirect_spaceflight",
        "core",
        "innate_stromal_counterresponse",
        "The direction is exceptionally consistent, but direct thymus-specific TLR validation is absent.",
        "Shimizu 2023; Okamura 2024",
        "Treat as an innate counter-response hypothesis rather than established inflammation.",
    ),
    ("thymus", "Lymphoid-stromal interactions"): (
        "tissue_mechanism_only",
        "core",
        "adaptive_proliferative_niche",
        "Thymocyte development requires stromal contact and migration; direct spaceflight validation is lacking.",
        "Horie 2019; Savino 2000",
        "Cannot distinguish weaker interaction state from loss of interacting cell populations.",
    ),
    ("thymus", "RHOA cytoskeletal cycle"): (
        "tissue_mechanism_only",
        "supporting",
        "adaptive_proliferative_niche",
        "RHOA is biologically compatible with thymocyte movement and adhesion but was not directly tested in flight thymus.",
        "Savino 2003; Horie 2019",
        "Use as support for the niche-coordination axis, not as an independent mechanism.",
    ),
    ("thymus", "Extracellular matrix organization"): (
        "direct_process_spaceflight",
        "core",
        "innate_stromal_counterresponse",
        "Spaceflight tissue adaptation converges on ECM and developmental regulation, supporting a stromal response.",
        "Grandke 2026; Horie 2019",
        "Higher organization score does not establish fibrosis or increased matrix mass.",
    ),
    ("skin", "Cell-cell junction organization"): (
        "direct_tissue_spaceflight",
        "core",
        "communication_barrier_maintenance",
        "Human spatial skin data directly report lower junctional and barrier transcripts after flight.",
        "Park 2024; Cope 2024",
        "Mouse bulk composition and human post-flight sampling differ from the present contrasts.",
    ),
    ("skin", "Keratinization"): (
        "direct_tissue_spaceflight",
        "core",
        "communication_barrier_maintenance",
        "Lower filaggrin and CASP14 in most murine skin subsets directly support reduced cornified-barrier differentiation.",
        "Cope 2024; Neutelings 2015",
        "The RR-5 dorsal recovery group reverses direction, and individual keratin genes can differ by compartment.",
    ),
    ("skin", "Chromatin-modifying enzymes"): (
        "direct_process_spaceflight",
        "supporting",
        "communication_barrier_maintenance",
        "Mouse flight skin shows epigenetic plasticity, but not a validated global decrease in chromatin-modifier activity.",
        "Sarkar 2024; Cope 2024",
        "Direction and functional consequence are unresolved.",
    ),
    ("skin", "Gap-junction trafficking"): (
        "process_supported",
        "core",
        "communication_barrier_maintenance",
        "Flight and microgravity studies support disrupted cell adhesion and wound-cell communication, though not this exact Reactome term.",
        "Park 2024; Zhao 2025; Chen 2023",
        "Gap-junction trafficking is a specific hypothesis within a broader communication phenotype.",
    ),
    ("skin", "Phase II detoxification"): (
        "direction_ambiguous",
        "not_core",
        "communication_barrier_maintenance",
        "The broad phase-II node is lower, but its fully nested glutathione-conjugation child is higher, so a uniform detoxification direction is not identifiable.",
        "Mao 2014",
        "Retain as context-sensitive; neither latent node is a direct enzyme-activity, glutathione-abundance, or metabolite-flux measurement.",
    ),
    ("skin", "Hedgehog signaling"): (
        "tissue_mechanism_only",
        "supporting",
        "communication_barrier_maintenance",
        "Flight changes hair cycling, and Hedgehog is required for follicular regeneration, but flight-specific Hedgehog activity was not measured.",
        "Neutelings 2015; Lim 2018",
        "This is a regenerative-niche hypothesis, not direct pathway confirmation.",
    ),
    ("skin", "Sphingolipid metabolism"): (
        "tissue_mechanism_only",
        "supporting",
        "communication_barrier_maintenance",
        "Sphingolipids are required for epidermal permeability, but skin-specific flight direction has not been independently established.",
        "Holleran 1991; Lyu 2022; Cope 2024",
        "Do not equate this score with measured ceramide abundance.",
    ),
    ("liver", "Regulation of insulin secretion"): (
        "direct_tissue_spaceflight",
        "core",
        "metabolic_homeostasis",
        "Spaceflight liver studies directly report inhibition of insulin-related signaling and metabolic dysfunction.",
        "Mathyk 2024; Beheshti 2019",
        "The Reactome term does not imply that hepatocytes secrete insulin.",
    ),
    ("liver", "Glutathione conjugation"): (
        "direction_conflict",
        "not_core",
        "metabolic_homeostasis",
        "A compensatory increase is possible, but direct flight data show lower hepatic GSH pools and glutathione reducibility.",
        "Kurosawa 2021; Jonscher 2016",
        "Heterogeneous accession direction and weak gene support preclude a clear complementary claim.",
    ),
    ("liver", "MHC class II antigen presentation"): (
        "indirect_spaceflight",
        "core",
        "immune_mechanical_coordination",
        "Systemic spaceflight studies support diminished antigen-presentation capacity; liver-specific localization is untested.",
        "Gridley 2015; da Silveira 2020",
        "Likely sensitive to hepatic immune-cell abundance.",
    ),
    ("liver", "T-cell receptor signaling"): (
        "indirect_spaceflight",
        "core",
        "immune_mechanical_coordination",
        "The lower score coheres with immune dysfunction but probably reports resident or infiltrating lymphocytes rather than hepatocytes.",
        "da Silveira 2020; Gridley 2015",
        "Treat MHC II and TCR as one immune-composition axis, not two independent mechanisms.",
    ),
    ("liver", "Rho-family GTPase cycle"): (
        "process_supported",
        "supporting",
        "immune_mechanical_coordination",
        "Recent work supports hepatic mechanotransduction in spaceflight, but Rho-family causality was not tested.",
        "Li 2026; Jonscher 2016",
        "Use the broader term cytoskeletal-mechanical regulation.",
    ),
    ("liver", "Extracellular matrix organization"): (
        "direction_ambiguous",
        "not_core",
        "immune_mechanical_coordination",
        "Prior work supports hepatic ECM remodeling but does not establish a uniform lower organization direction.",
        "Jonscher 2016; Grandke 2026",
        "Seven of ten independent cohort sources agree and gene support is weak; retain as context-sensitive.",
    ),
    ("soleus", "Extracellular matrix degradation"): (
        "direct_process_spaceflight",
        "supporting",
        "matrix_damage_response",
        "Human soleus and muscle proteomics support ECM loss or remodeling after flight.",
        "Gambara 2017; Tascher 2017; Salanova 2023; Murgia 2024",
        "Pathway-member gene support is sparse and opposite-direction genes are present.",
    ),
    ("soleus", "Glycosaminoglycan metabolism"): (
        "tissue_mechanism_only",
        "exploratory",
        "matrix_damage_response",
        "Glycosaminoglycans are matrix components, but direct soleus spaceflight evidence for this pathway is absent.",
        "Murgia 2024",
        "Only three accessions and minimal gene support.",
    ),
    ("soleus", "DNA repair"): (
        "indirect_spaceflight",
        "exploratory",
        "matrix_damage_response",
        "Radiation and oxidative stress make DNA repair plausible, but soleus-specific validation is absent.",
        "Gambara 2017; Tascher 2017",
        "Only three accessions, modest effect after restriction, and weak gene support.",
    ),
}


def parse_gmt() -> dict[str, set[str]]:
    pathways = {}
    with (ROOT / "data/pathways/reactome_current_mouse_ensembl.gmt").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            pathways[fields[0]] = set(fields[2:])
    return pathways


def run() -> None:
    selected = pd.read_csv(SOURCE_DIR / "table_1_curated_pathway_results.tsv", sep="\t")
    gmt = parse_gmt()
    audit_rows = []
    pair_rows = []

    for tissue, run_dir in RUNS.items():
        keys = [key for key in REVIEW if key[0] == tissue]
        subset = selected.loc[
            selected["tissue"].eq(tissue)
            & selected["short_label"].isin([key[1] for key in keys])
        ].copy()
        terms = subset["term"].tolist()
        scores = pd.read_csv(
            run_dir / "query_pathway_scores.tsv",
            sep="\t",
            usecols=["id.accession"] + terms,
        )
        orientation = dict(zip(subset["term"], subset["latent_orientation"]))
        for term in terms:
            scores[term] = scores[term] * orientation[term]
        residual = scores[terms] - scores.groupby(scores["id.accession"])[terms].transform("mean")
        correlation = residual.corr()

        for first_index, first in enumerate(terms):
            for second in terms[first_index + 1 :]:
                first_genes = gmt[first]
                second_genes = gmt[second]
                shared = first_genes & second_genes
                union = first_genes | second_genes
                pair_rows.append(
                    {
                        "tissue": tissue,
                        "first_term": first,
                        "first_label": subset.loc[subset["term"].eq(first), "short_label"].iloc[0],
                        "second_term": second,
                        "second_label": subset.loc[subset["term"].eq(second), "short_label"].iloc[0],
                        "within_accession_residual_score_correlation": correlation.loc[first, second],
                        "shared_reactome_genes": len(shared),
                        "reactome_gene_jaccard": len(shared) / len(union) if union else np.nan,
                    }
                )

        for row in subset.itertuples():
            key = (tissue, row.short_label)
            evidence_tier, story_status, story_axis, assessment, sources, caution = REVIEW[key]
            other_terms = [term for term in terms if term != row.term]
            mean_abs_correlation = (
                correlation.loc[row.term, other_terms].abs().mean() if other_terms else np.nan
            )
            audit_rows.append(
                {
                    "tissue": tissue,
                    "term": row.term,
                    "short_label": row.short_label,
                    "reviewed_role": row.evidence_role,
                    "direction": "higher_in_flight" if row.mean_accession_effect > 0 else "lower_in_flight",
                    "study_balanced_effect": row.mean_accession_effect,
                    "accessions_same_direction": row.n_accession_same_direction,
                    "accessions_opposite_direction": row.n_accession_opposite_direction,
                    "meta_fdr": row.meta_fdr,
                    "maximum_leave_one_out_fdr": row.maximum_leave_one_out_fdr,
                    "gene_support_fraction": row.fraction_pathway_genes_pooled_fdr_lt_005,
                    "mean_absolute_within_accession_score_correlation": mean_abs_correlation,
                    "evidence_tier": evidence_tier,
                    "story_status": story_status,
                    "story_axis": story_axis,
                    "literature_assessment": assessment,
                    "supporting_sources": sources,
                    "interpretive_caution": caution,
                }
            )

    audit = pd.DataFrame(audit_rows).sort_values(["tissue", "story_status", "short_label"])
    pairs = pd.DataFrame(pair_rows).sort_values(
        ["tissue", "within_accession_residual_score_correlation"], ascending=[True, False]
    )
    audit.to_csv(SOURCE_DIR / "table_s6_program_story_audit.tsv", sep="\t", index=False)
    pairs.to_csv(SOURCE_DIR / "table_s7_program_pairwise_structure.tsv", sep="\t", index=False)
    print(
        f"Audited {len(audit)} programs and {len(pairs)} within-tissue pairs; "
        f"core={audit['story_status'].eq('core').sum()}"
    )


if __name__ == "__main__":
    run()
