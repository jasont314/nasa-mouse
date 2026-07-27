"""Build the expanded, family-level pathway review for the ASGSR paper."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
PAPER_DIR = ROOT / "paper/asgsr_expimap_hvg"
SOURCE_DIR = PAPER_DIR / "source_data"
FIGURE_DIR = PAPER_DIR / "figures"

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

ROLE_COLORS = {
    "aligned": "#24834b",
    "complementary": "#2869a8",
    "context_sensitive": "#b36a19",
    "not_interpretable": "#777777",
}


@dataclass(frozen=True)
class FamilyReview:
    name: str
    role: str
    decision: str
    representative: str
    synthesis: str
    references: str
    caution: str
    labels: tuple[str, ...]


# Every label is assigned explicitly. This is a pathway-by-pathway review map, not a
# keyword classifier. The grouping uses Reactome hierarchy, gene-set overlap, tissue
# context, and the literature assessment recorded below.
FAMILIES = {
    "liver": {
        "xenobiotic_redox": FamilyReview(
            "Xenobiotic and redox metabolism",
            "context_sensitive",
            "existing_context",
            "Cytochrome P450",
            "Known hepatic xenobiotic biology is recovered, but CYP, biological-oxidation, drug-specific ADME, and glutathione terms do not define one uniform direction.",
            "17;19;20",
            "CYP directions split by cohort, prednisone is exposure-specific, and glutathione score does not measure hepatic glutathione abundance.",
            (
                "Cytochrome P450",
                "Prednisone Adme",
                "Biological Oxidations",
                "Glutathione conjugation",
            ),
        ),
        "lipid_endocrine": FamilyReview(
            "Lipid and endocrine regulation",
            "context_sensitive",
            "existing_context",
            "Regulation of insulin secretion",
            "The family recovers established lipid and endocrine disruption, but its child programs separate into lower insulin and cholesterol scores and higher steroid, phospholipid, peroxisomal, and AKT-associated scores.",
            "17;18;21;23",
            "The Reactome insulin-secretion label is not evidence that hepatocytes secrete insulin, and nested signaling scores should not be treated as independent mechanisms.",
            (
                "Metabolism Of Steroid Hormones",
                "Peroxisomal Protein Import",
                "Pip3 Activates Akt Signaling",
                "Regulation of insulin secretion",
                "Cholesterol Biosynthesis",
                "Glycerophospholipid Biosynthesis",
                "Binding And Uptake Of Ligands By Scavenger Receptors",
            ),
        ),
        "immune_effector": FamilyReview(
            "Adaptive and innate effector programs",
            "complementary",
            "existing_core_axis",
            "MHC class II antigen presentation",
            "Lower MHC-I and MHC-II presentation, T-cell receptor, interleukin, complement, and neutrophil-degranulation scores extend the current lower hepatic immune-communication axis.",
            "10;22",
            "Bulk liver cannot distinguish regulation within resident cells from altered abundance of lymphocytes, Kupffer cells, neutrophils, or other non-parenchymal cells.",
            (
                "Immune System",
                "MHC class II antigen presentation",
                "Neutrophil Degranulation",
                "Complement Cascade",
                "Class I MHC Mediated Antigen Processing Presentation",
                "Signaling By Interleukins",
                "T-cell receptor signaling",
            ),
        ),
        "toll_like": FamilyReview(
            "Toll-like receptor branches",
            "context_sensitive",
            "supplement_only",
            "Toll Like Receptor 4 TLR4 Cascade",
            "The TLR4-specific node is higher while the broad Toll-like-receptor node is lower, so the expanded screen does not support a pathway-wide TLR direction.",
            "18;22",
            "The TLR4 genes are fully nested within the broad TLR set; opposite latent directions may reflect branch allocation or correlated-decoder competition rather than selective pathway activation.",
            (
                "Toll Like Receptor 4 TLR4 Cascade",
                "Toll Like Receptor Cascades",
            ),
        ),
        "matrix_mechanical_vascular": FamilyReview(
            "Matrix, mechanical, and vascular regulation",
            "complementary",
            "existing_support_axis",
            "Rho-family GTPase cycle",
            "Predominantly lower Rho-family, matrix, proteoglycan, VEGF, vascular-wall, and platelet-calcium scores reinforce a structural and non-parenchymal regulation layer beside metabolic dysfunction.",
            "18;23",
            "The family does not establish lower matrix abundance, lower fibrosis, impaired blood flow, or a causal mechanical stimulus.",
            (
                "Cell Surface Interactions At The Vascular Wall",
                "Rho-family GTPase cycle",
                "Extracellular matrix organization",
                "Response To Elevated Platelet Cytosolic Ca2",
                "Rho GTPase Effectors",
                "Signaling By Vegf",
                "ECM Proteoglycans",
            ),
        ),
        "proliferation_repair": FamilyReview(
            "Proliferation and homology-directed repair",
            "context_sensitive",
            "supplement_only",
            "Homology Directed Repair",
            "Mitotic, broad cell-cycle, and homology-directed-repair programs are lower, but this family is not a prespecified hepatic spaceflight phenotype.",
            "22",
            "The pattern may reflect cell composition or proliferative state and lacks direct liver-specific functional validation.",
            ("Cell Cycle Mitotic", "Cell Cycle", "Homology Directed Repair"),
        ),
        "trafficking_proteostasis": FamilyReview(
            "Organelle trafficking and proteostasis",
            "context_sensitive",
            "supplement_only",
            "Copi Dependent Golgi To ER Retrograde Traffic",
            "Lower COPI trafficking contrasts with higher mitochondrial-protein degradation and neddylation scores, providing no coherent family-wide direction.",
            "22",
            "These broad latent programs overlap many cell types and are not direct measurements of organelle flux or protein turnover.",
            (
                "Copi Dependent Golgi To ER Retrograde Traffic",
                "Mitochondrial Protein Degradation",
                "Neddylation",
            ),
        ),
        "broad_receptor_signaling": FamilyReview(
            "Broad receptor and stress signaling",
            "not_interpretable",
            "exclude",
            "Signaling By GPCR",
            "Umbrella receptor, second-messenger, and stimulus-response terms are too broad to support a distinct hepatic mechanism.",
            "21;22",
            "The constituent receptors and responding cell populations are unresolved.",
            (
                "Peptide Ligand Binding Receptors",
                "Intracellular Signaling By Second Messengers",
                "Signaling By GPCR",
                "Cellular Responses To Stimuli",
            ),
        ),
        "off_context": FamilyReview(
            "Tissue-incongruent labels",
            "not_interpretable",
            "exclude",
            "Cardiac Conduction",
            "Cardiac-conduction and neuronal-system labels do not provide credible liver-specific pathway interpretations.",
            "",
            "Shared ion-channel or receptor genes can generate these labels without cardiac or neuronal tissue activity.",
            ("Cardiac Conduction", "Neuronal System"),
        ),
    },
    "skin": {
        "epidermal_development": FamilyReview(
            "Epidermal differentiation and regeneration",
            "aligned",
            "existing_core_axis",
            "Keratinization",
            "Lower keratinization, developmental, Hedgehog, and neural-development programs support reduced epidermal differentiation and regenerative-niche activity.",
            "12;13;15",
            "The broad developmental term combines multiple compartments, and sphingolipid score is not a direct ceramide measurement.",
            (
                "Keratinization",
                "Developmental Biology",
                "Nervous System Development",
                "Hedgehog signaling",
                "Sphingolipid metabolism",
            ),
        ),
        "communication_adhesion": FamilyReview(
            "Cell communication, adhesion, and cytoskeleton",
            "complementary",
            "existing_core_axis",
            "Gap-junction trafficking",
            "Lower gap-junction, CDH1, Rho-family, second-messenger, and cell-junction scores form a coherent lower tissue-coordination family.",
            "13;16",
            "The result does not identify a specific skin layer or prove reduced functional junctional coupling.",
            (
                "Gap-junction trafficking",
                "Regulation Of CDH1 Expression And Function",
                "Negative Regulation Of CDH1 Gene Transcription",
                "Signaling By Rho GTPases",
                "Rho GTPase Cycle",
                "Intracellular Signaling By Second Messengers",
                "Cell-cell junction organization",
            ),
        ),
        "genome_regulation": FamilyReview(
            "Chromatin, proliferation, and genome maintenance",
            "context_sensitive",
            "existing_support_axis",
            "Chromatin-modifying enzymes",
            "Chromatin, nucleotide, cell-cycle, senescence, and DNA-repair scores are predominantly lower and reinforce a reduced regulatory and proliferative state.",
            "12;13",
            "Several terms divide three versus three across accessions, so the family is a protocol-qualified state rather than a universal repair defect.",
            (
                "Chromatin-modifying enzymes",
                "Metabolism Of Nucleotides",
                "Cell Cycle",
                "DNA Damage Telomere Stress Induced Senescence",
                "Chromatin Organization",
                "M Phase",
                "DNA repair",
                "Mitotic G2 G2 M Phases",
            ),
        ),
        "cutaneous_muscle": FamilyReview(
            "Cutaneous striated muscle",
            "aligned",
            "add_aligned_context",
            "Muscle Contraction",
            "Higher muscle and striated-muscle contraction scores recover the known panniculus-carnosus response in full-thickness mouse skin.",
            "12;15",
            "This is a skin-compartment signal, not evidence that epidermal cells acquired a muscle program; the direction is not universal across all projects.",
            ("Muscle Contraction", "Striated Muscle Contraction"),
        ),
        "innate_phagocytosis": FamilyReview(
            "Innate immune and phagocytic response",
            "aligned",
            "add_aligned_context",
            "Fcgamma Receptor Fcgr Dependent Phagocytosis",
            "Higher Fc-receptor phagocytosis and broad innate-immune scores recover a known inflammatory and immune component of flight-exposed skin.",
            "12;13",
            "Whole-skin scores can reflect immune-cell abundance and do not establish activation within keratinocytes.",
            (
                "Fcgamma Receptor Fcgr Dependent Phagocytosis",
                "Innate Immune System",
            ),
        ),
        "detox_endocrine": FamilyReview(
            "Detoxification and endocrine metabolism",
            "context_sensitive",
            "revise_existing_claim",
            "Phase II detoxification",
            "Lower broad phase-II conjugation contrasts with higher nested glutathione-conjugation, peptide-hormone, and steroid-metabolism scores; no uniform detoxification direction is supported.",
            "14",
            "Glutathione conjugation is fully nested in phase-II conjugation and the two latent scores are strongly anticorrelated, so branch-specific mechanistic claims are not identifiable here.",
            (
                "Peptide Hormone Metabolism",
                "Glutathione Conjugation",
                "Phase II detoxification",
                "Metabolism Of Steroids",
            ),
        ),
        "protein_processing": FamilyReview(
            "Secretory trafficking and protein processing",
            "context_sensitive",
            "supplement_only",
            "Copi Dependent Golgi To ER Retrograde Traffic",
            "Lower COPI and N-linked-glycosylation scores contrast with higher broad post-translational modification, preventing a coherent process-level direction.",
            "12",
            "These programs are broad and do not measure protein-processing flux.",
            (
                "Copi Dependent Golgi To ER Retrograde Traffic",
                "Asparagine N Linked Glycosylation",
                "Post Translational Protein Modification",
            ),
        ),
        "sensory_broad": FamilyReview(
            "Sensory and broad signaling labels",
            "not_interpretable",
            "exclude",
            "Visual Phototransduction",
            "Visual, neuronal, potassium-channel, NRAGE-JNK, and broad GPCR labels do not form a specific skin-spaceflight mechanism.",
            "15",
            "Cutaneous sensory cells exist, but these aggregate labels cannot localize them and include tissue-incongruent annotation.",
            (
                "Visual Phototransduction",
                "Neuronal System",
                "Nrage Signals Death Through Jnk",
                "GPCR Downstream Signalling",
                "Potassium Channels",
            ),
        ),
    },
    "soleus": {
        "immune_inflammatory": FamilyReview(
            "Immune and cytokine programs",
            "context_sensitive",
            "existing_context",
            "Immune-system signaling",
            "Broad immune and neutrophil scores are driven by the confounded accession, whereas interferon and cytokine scores are lower across the restricted comparison.",
            "24;27",
            "Nested immune programs point in opposite directions and cannot support a universal inflammatory claim.",
            (
                "Immune-system signaling",
                "Neutrophil Degranulation",
                "Interferon Signaling",
                "Cytokine signaling",
            ),
        ),
        "epithelial_incongruent": FamilyReview(
            "Cornified-envelope and keratin labels",
            "not_interpretable",
            "exclude",
            "Formation Of The Cornified Envelope",
            "Cornified-envelope and keratinization scores are tissue-incongruent for dissected soleus and are not interpreted as muscle biology.",
            "",
            "The signal could reflect dissection contamination, shared structural genes, or model annotation rather than a muscle mechanism.",
            ("Formation Of The Cornified Envelope", "Keratinization"),
        ),
        "matrix_turnover": FamilyReview(
            "Matrix disassembly and support",
            "complementary",
            "existing_core_axis",
            "Extracellular matrix degradation",
            "Higher collagen and extracellular-matrix degradation together with lower elastic-fibre and glycosaminoglycan scores support a matrix-disassembly and reduced-support state.",
            "24;26;28",
            "Only two unconfounded accessions remain and member-gene support is sparse; matrix abundance and protease activity were not measured.",
            (
                "Collagen Degradation",
                "Extracellular matrix degradation",
                "Molecules Associated With Elastic Fibres",
                "Glycosaminoglycan metabolism",
            ),
        ),
        "growth_damage": FamilyReview(
            "Growth, apoptosis, and DNA-damage response",
            "context_sensitive",
            "supplement_only",
            "DNA repair",
            "Higher apoptosis and DNA-repair scores persist after restriction, while broad cell-cycle and MLL3-MLL4 scores remain accession sensitive.",
            "24;25",
            "The family combines distinct processes and does not establish myofiber apoptosis or repair activity.",
            (
                "Cell Cycle",
                "Apoptosis",
                "DNA repair",
                "Epigenetic Regulation Of Gene Expression By Mll3 And Mll4 Complexes",
            ),
        ),
        "hemostatic": FamilyReview(
            "Hemostatic and platelet programs",
            "context_sensitive",
            "supplement_only",
            "Hemostasis",
            "Higher hemostasis and platelet-degranulation scores may reflect vascular or blood-cell composition rather than a myofiber program.",
            "",
            "Direct soleus-specific spaceflight validation is absent.",
            ("Hemostasis", "Platelet Degranulation"),
        ),
        "metabolic_transport": FamilyReview(
            "Metabolism and transmembrane transport",
            "context_sensitive",
            "existing_context",
            "Fatty-acid metabolism",
            "Broad metabolism and transport parents oppose several lipid, oxidation, vitamin, and SLC children, and the fatty-acid direction changes after accession restriction.",
            "24;25;26",
            "The family is not directionally identifiable as one metabolic response.",
            (
                "Fatty-acid metabolism",
                "Metabolism",
                "Biological Oxidations",
                "Transport Of Small Molecules",
                "Metabolism Of Lipids",
                "Slc Mediated Transmembrane Transport",
                "Metabolism Of Vitamins And Cofactors",
            ),
        ),
        "trophic_vascular": FamilyReview(
            "Trophic, adhesion, and vascular signaling",
            "aligned",
            "add_aligned_context",
            "Pip3 Activates Akt Signaling",
            "Lower PI3K-AKT, WNT, and VEGF programs after restriction align with reduced anabolic, focal-adhesion, and vascular-support signaling during soleus unloading.",
            "25;27",
            "PI3K-AKT attenuates strongly after excluding OSD-714, and broad GPCR and Rho parents are accession sensitive; only the lower trophic subfamily is retained.",
            (
                "Pip3 Activates Akt Signaling",
                "Signaling By Wnt",
                "Signaling By Vegf",
            ),
        ),
        "broad_gpcr_cytoskeleton": FamilyReview(
            "Broad GPCR and Rho-family signaling",
            "context_sensitive",
            "supplement_only",
            "G Alpha I Signalling Events",
            "Broad GPCR and Rho-family parents move in accession-sensitive or opposing directions and do not define the lower trophic subfamily.",
            "24;25",
            "Class-A receptors and G-alpha-i are nested within the GPCR parent, while the broad Rho term changes direction after restriction.",
            (
                "Signaling By GPCR",
                "Signaling By Rho GTPases",
                "G Alpha I Signalling Events",
                "Class A 1 Rhodopsin Like Receptors",
            ),
        ),
        "trafficking_cytoskeleton": FamilyReview(
            "Membrane trafficking and intracellular transport",
            "context_sensitive",
            "supplement_only",
            "Vesicle Mediated Transport",
            "GPI-anchor, membrane, kinesin, and vesicle programs are higher or preserved after restriction but lack direct soleus-spaceflight validation.",
            "26",
            "The terms are broad and should not be interpreted as measured trafficking rates.",
            (
                "Post Translational Modification Synthesis Of Gpi Anchored Proteins",
                "Membrane Trafficking",
                "Kinesins",
                "Vesicle Mediated Transport",
            ),
        ),
        "neuromuscular": FamilyReview(
            "Contractile and sensory programs",
            "context_sensitive",
            "existing_context",
            "Striated-muscle contraction",
            "The contraction direction changes after exclusion, while the broad sensory score remains higher without a specific neuromuscular interpretation.",
            "24;25;29;30",
            "Fiber-type transitions can move individual contractile genes in opposite directions.",
            ("Sensory Perception", "Striated-muscle contraction"),
        ),
    },
    "thymus": {
        "proliferation_genome": FamilyReview(
            "Proliferation and genome maintenance",
            "aligned",
            "existing_core_axis",
            "Mitotic cell cycle",
            "Lower mitotic, nucleotide, chromatin, transcription, telomere, rRNA, and DNA-repair programs form one dominant reduced proliferative and genome-maintenance state.",
            "6;9",
            "Many child terms may reflect fewer cycling thymocytes rather than independent defects in every molecular process.",
            (
                "Mitotic cell cycle",
                "DNA Damage Telomere Stress Induced Senescence",
                "Metabolism Of Nucleotides",
                "DNA repair",
                "Epigenetic Regulation Of Adipogenesis Genes By Mll3 And Mll4 Complexes",
                "Chromatin Modifying Enzymes",
                "Generic Transcription Pathway",
                "Cell Cycle",
                "Telomere Maintenance",
                "Major Pathway Of Rrna Processing In The Nucleolus And Cytosol",
            ),
        ),
        "death_stress": FamilyReview(
            "Cell death and stress response",
            "context_sensitive",
            "add_interpretive_context",
            "Apoptosis",
            "Lower apoptosis, programmed-cell-death, and stress-response scores support a proliferation-led long-duration phenotype rather than universal ongoing apoptosis.",
            "7;9",
            "Post-return STS-135 tissue showed DNA fragmentation, so timing and mission context prevent a general anti-apoptosis claim.",
            ("Programmed Cell Death", "Apoptosis", "Cellular Responses To Stress"),
        ),
        "trafficking_proteostasis": FamilyReview(
            "Secretory trafficking and proteostasis",
            "context_sensitive",
            "supplement_only",
            "Copi Dependent Golgi To ER Retrograde Traffic",
            "Lower retrograde trafficking and protein-metabolism scores contrast with higher N-linked glycosylation and broad post-translational modification.",
            "22",
            "Opposing nested terms and absent thymus-specific validation preclude a process-level direction.",
            (
                "Copi Dependent Golgi To ER Retrograde Traffic",
                "Golgi To ER Retrograde Transport",
                "Asparagine N Linked Glycosylation",
                "Post Translational Protein Modification",
                "Neddylation",
                "Metabolism Of Proteins",
            ),
        ),
        "niche_cytoskeleton": FamilyReview(
            "Thymocyte-niche interaction and cytoskeleton",
            "complementary",
            "existing_core_axis",
            "RHOA cytoskeletal cycle",
            "Lower RHOA, broader Rho-family, and lymphoid-stromal interaction scores reinforce reduced thymocyte migration, adhesion, or niche coordination.",
            "6;11",
            "Bulk tissue cannot separate weaker interaction state from loss of interacting thymocytes.",
            (
                "RHOA cytoskeletal cycle",
                "Rho GTPase Cycle",
                "Lymphoid-stromal interactions",
            ),
        ),
        "stromal_matrix": FamilyReview(
            "Stromal matrix and TGF-beta response",
            "complementary",
            "existing_core_axis",
            "Extracellular matrix organization",
            "Higher extracellular-matrix and TGF-beta-receptor scores reinforce a stromal response accompanying thymocyte loss.",
            "9;11",
            "Higher latent scores do not establish fibrosis, increased matrix mass, or the responding cell population.",
            ("Extracellular matrix organization", "Signaling By Tgf Beta Receptor Complex"),
        ),
        "metabolic_endocrine": FamilyReview(
            "Respiratory, metabolic, and endocrine state",
            "context_sensitive",
            "supplement_only",
            "Aerobic Respiration And Respiratory Electron Transport",
            "Higher respiratory and broad metabolic scores coexist with lower steroid and ESR-associated scores, suggesting composition or state redistribution rather than globally improved metabolism.",
            "8;22",
            "Latent respiratory score is not a direct measurement of mitochondrial function or ATP production.",
            (
                "Aerobic Respiration And Respiratory Electron Transport",
                "Metabolism",
                "Metabolism Of Steroids",
                "Esr Mediated Signaling",
                "Transport Of Small Molecules",
            ),
        ),
        "innate_hemostatic": FamilyReview(
            "Innate, hemostatic, and vascular-associated response",
            "complementary",
            "existing_support_axis",
            "Innate TLR signaling",
            "Higher broad immune, TLR, neutrophil, hemostasis, platelet, and heme-scavenging scores reinforce an innate and vascular-associated counter-response.",
            "9;10",
            "Blood contamination or altered immune-cell composition may contribute, and these terms should not be presented as six independent mechanisms.",
            (
                "Scavenging Of Heme From Plasma",
                "Hemostasis",
                "Immune System",
                "Neutrophil Degranulation",
                "Platelet Activation Signaling And Aggregation",
                "Innate TLR signaling",
            ),
        ),
        "adaptive_tcr": FamilyReview(
            "Adaptive T-cell signaling",
            "aligned",
            "existing_core_axis",
            "T-cell receptor signaling",
            "Lower T-cell receptor score anchors the established reduction in adaptive thymic output and signaling.",
            "8;10",
            "The score is sensitive to the abundance of developing T cells in bulk tissue.",
            ("T-cell receptor signaling",),
        ),
        "broad_signaling": FamilyReview(
            "Broad receptor signaling",
            "not_interpretable",
            "exclude",
            "Signal Transduction",
            "Signal-transduction, GPCR, and G-alpha-i umbrellas are too broad to define a thymus mechanism.",
            "",
            "The terms contain many receptors and cell types with unrelated functions.",
            ("Signal Transduction", "Signaling By GPCR", "G Alpha I Signalling Events"),
        ),
        "off_context": FamilyReview(
            "Tissue-incongruent labels",
            "not_interpretable",
            "exclude",
            "Neurotransmitter Receptors And Postsynaptic Signal Transmission",
            "Muscle-contraction, postsynaptic, and potassium-channel labels do not support a specific thymic mechanism.",
            "",
            "Shared receptors, ion channels, or stromal genes can produce these labels without muscle or neuronal tissue activity.",
            (
                "Muscle Contraction",
                "Neurotransmitter Receptors And Postsynaptic Signal Transmission",
                "Potassium Channels",
            ),
        ),
    },
}


def parse_gmt() -> dict[str, set[str]]:
    pathways: dict[str, set[str]] = {}
    with (ROOT / "data/pathways/reactome_current_mouse_ensembl.gmt").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            pathways[fields[0]] = set(fields[2:])
    return pathways


def candidate_screen(screen: pd.DataFrame) -> pd.DataFrame:
    candidates = []
    for tissue, group in screen.groupby("tissue", sort=False):
        active = group.loc[group["active_latent_program"]].copy()
        top_n = int(np.ceil(len(active) * 0.10))
        top_decile = set(
            active.nsmallest(top_n, "absolute_effect_rank_active")["term"]
        )
        stable_extension = set(
            active.loc[
                active["absolute_effect_rank_active"].le(40)
                & active["primary_direction_agreement"].ge(0.8),
                "term",
            ]
        )
        prespecified = set(group.loc[group["curated_for_main_figures"], "term"])
        selected = top_decile | stable_extension | prespecified
        subset = group.loc[group["term"].isin(selected)].copy()
        subset["candidate_reason"] = subset["term"].map(
            lambda term: ";".join(
                reason
                for reason, terms in (
                    ("top_decile_magnitude", top_decile),
                    ("stable_rank_40_extension", stable_extension),
                    ("prespecified_main_review", prespecified),
                )
                if term in terms
            )
        )
        subset["top_decile_cutoff"] = top_n
        candidates.append(subset)
    return pd.concat(candidates, ignore_index=True)


def family_lookup() -> dict[tuple[str, str], tuple[str, FamilyReview]]:
    lookup: dict[tuple[str, str], tuple[str, FamilyReview]] = {}
    for tissue, families in FAMILIES.items():
        for family_id, review in families.items():
            for label in review.labels:
                key = (tissue, label)
                if key in lookup:
                    raise RuntimeError(f"Duplicate family assignment: {key}")
                lookup[key] = (family_id, review)
    return lookup


def add_gene_support(candidates: pd.DataFrame, gmt: dict[str, set[str]]) -> pd.DataFrame:
    genes = pd.read_csv(SOURCE_DIR / "table_s4_gene_level_results.tsv.gz", sep="\t")
    rows = []
    for row in candidates.itertuples():
        measured = genes.loc[
            genes["tissue"].eq(row.tissue) & genes["gene_id"].isin(gmt[row.term])
        ]
        significant = measured.loc[measured["pooled_fdr"].lt(0.05)]
        pathway_sign = np.sign(row.mean_accession_effect)
        gene_sign = np.sign(significant["study_balanced_log2cpm_flight_minus_ground"])
        rows.append(
            {
                "tissue": row.tissue,
                "term": row.term,
                "pathway_genes_measured": len(measured),
                "pathway_genes_pooled_fdr_lt_005": len(significant),
                "pathway_gene_support_fraction": (
                    len(significant) / len(measured) if len(measured) else np.nan
                ),
                "significant_genes_same_direction": int((gene_sign == pathway_sign).sum()),
                "significant_genes_opposite_direction": int((gene_sign == -pathway_sign).sum()),
            }
        )
    return candidates.merge(pd.DataFrame(rows), on=["tissue", "term"], how="left")


def add_sensitivity(candidates: pd.DataFrame) -> pd.DataFrame:
    result = candidates.copy()
    result["sensitivity_analysis"] = ""
    result["sensitivity_effect"] = np.nan
    for tissue, label, column in (
        ("thymus", "exclude OSD-289", "restricted_study_balanced_effect"),
        ("skin", "four-project paired-site balance", "skin_project_balanced_effect"),
        ("liver", "original 12-accession input", "full_input_study_balanced_effect"),
        ("soleus", "exclude OSD-714", "restricted_study_balanced_effect"),
    ):
        mask = result["tissue"].eq(tissue)
        result.loc[mask, "sensitivity_analysis"] = label
        result.loc[mask, "sensitivity_effect"] = result.loc[mask, column]
    result["sensitivity_direction_preserved"] = (
        np.sign(result["mean_accession_effect"])
        == np.sign(result["sensitivity_effect"])
    )
    result["sensitivity_magnitude_ratio"] = (
        result["sensitivity_effect"].abs()
        / result["mean_accession_effect"].abs().replace(0, np.nan)
    )
    return result


def add_overlap_flags(candidates: pd.DataFrame, gmt: dict[str, set[str]]) -> pd.DataFrame:
    result = candidates.copy()
    result["opposite_nested_program"] = False
    result["opposite_nested_terms"] = ""
    result["maximum_opposite_nested_overlap_coefficient"] = np.nan
    for tissue, group in result.groupby("tissue"):
        terms = group["term"].tolist()
        effects = group.set_index("term")["mean_accession_effect"]
        opposite: dict[str, list[str]] = {term: [] for term in terms}
        maximum: dict[str, float] = {term: np.nan for term in terms}
        for first_index, first in enumerate(terms):
            for second in terms[first_index + 1 :]:
                first_genes = gmt[first]
                second_genes = gmt[second]
                shared = len(first_genes & second_genes)
                if not shared:
                    continue
                overlap = shared / min(len(first_genes), len(second_genes))
                if overlap < 0.8 or np.sign(effects[first]) == np.sign(effects[second]):
                    continue
                opposite[first].append(second)
                opposite[second].append(first)
                maximum[first] = np.nanmax([maximum[first], overlap])
                maximum[second] = np.nanmax([maximum[second], overlap])
        for term in terms:
            mask = result["tissue"].eq(tissue) & result["term"].eq(term)
            result.loc[mask, "opposite_nested_program"] = bool(opposite[term])
            result.loc[mask, "opposite_nested_terms"] = ";".join(opposite[term])
            result.loc[mask, "maximum_opposite_nested_overlap_coefficient"] = maximum[term]
    return result


def centered_score_correlations(
    candidates: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    correlations = {}
    for tissue, group in candidates.groupby("tissue"):
        terms = group["term"].tolist()
        scores = pd.read_csv(
            RUNS[tissue] / "query_pathway_scores.tsv",
            sep="\t",
            usecols=["id.accession"] + terms,
        )
        orientation = group.set_index("term")["latent_orientation"]
        scores[terms] = scores[terms].mul(orientation[terms], axis=1)
        centered = scores[terms] - scores.groupby(scores["id.accession"])[terms].transform(
            "mean"
        )
        correlations[tissue] = centered.corr()
    return correlations


def annotate_families(candidates: pd.DataFrame) -> pd.DataFrame:
    lookup = family_lookup()
    expected = set(zip(candidates["tissue"], candidates["screen_display_label"]))
    assigned = set(lookup)
    missing = sorted(expected - assigned)
    extra = sorted(assigned - expected)
    if missing or extra:
        raise RuntimeError(f"Family review mismatch; missing={missing}, extra={extra}")

    rows = []
    for row in candidates.itertuples():
        family_id, review = lookup[(row.tissue, row.screen_display_label)]
        values = row._asdict()
        values.update(
            {
                "family_id": family_id,
                "process_family": review.name,
                "family_literature_role": review.role,
                "family_narrative_decision": review.decision,
                "family_synthesis": review.synthesis,
                "family_supporting_references": review.references,
                "family_caution": review.caution,
                "family_representative": row.screen_display_label == review.representative,
                "pathway_review_disposition": (
                    "representative"
                    if row.screen_display_label == review.representative
                    else (
                        "excluded_tissue_incongruent_or_nonspecific"
                        if review.decision == "exclude"
                        else "family_support_or_redundant_child"
                    )
                ),
            }
        )
        rows.append(values)
    return pd.DataFrame(rows)


def summarize_families(
    pathways: pd.DataFrame,
    correlations: dict[str, pd.DataFrame],
    gmt: dict[str, set[str]],
) -> pd.DataFrame:
    rows = []
    for (tissue, family_id), group in pathways.groupby(["tissue", "family_id"]):
        review = FAMILIES[tissue][family_id]
        representative = group.loc[group["family_representative"]].iloc[0]
        terms = group["term"].tolist()
        pair_correlations = []
        pair_jaccards = []
        pair_overlap = []
        for first_index, first in enumerate(terms):
            for second in terms[first_index + 1 :]:
                pair_correlations.append(abs(correlations[tissue].loc[first, second]))
                shared = len(gmt[first] & gmt[second])
                pair_jaccards.append(shared / len(gmt[first] | gmt[second]))
                pair_overlap.append(shared / min(len(gmt[first]), len(gmt[second])))
        rows.append(
            {
                "tissue": tissue,
                "family_id": family_id,
                "process_family": review.name,
                "family_literature_role": review.role,
                "family_narrative_decision": review.decision,
                "n_candidate_pathways": len(group),
                "representative_term": representative["term"],
                "representative_label": representative["screen_display_label"],
                "representative_effect": representative["mean_accession_effect"],
                "representative_rank": representative["absolute_effect_rank_active"],
                "minimum_family_effect": group["mean_accession_effect"].min(),
                "maximum_family_effect": group["mean_accession_effect"].max(),
                "pathways_higher_in_flight": int(group["mean_accession_effect"].gt(0).sum()),
                "pathways_lower_in_flight": int(group["mean_accession_effect"].lt(0).sum()),
                "pathways_direction_preserved_in_sensitivity": int(
                    group["sensitivity_direction_preserved"].sum()
                ),
                "pathways_with_opposite_nested_program": int(
                    group["opposite_nested_program"].sum()
                ),
                "median_absolute_within_accession_score_correlation": (
                    float(np.median(pair_correlations)) if pair_correlations else np.nan
                ),
                "maximum_reactome_gene_jaccard": (
                    float(np.max(pair_jaccards)) if pair_jaccards else np.nan
                ),
                "maximum_reactome_gene_overlap_coefficient": (
                    float(np.max(pair_overlap)) if pair_overlap else np.nan
                ),
                "family_synthesis": review.synthesis,
                "supporting_references": review.references,
                "interpretive_caution": review.caution,
                "member_labels": ";".join(group["screen_display_label"]),
                "member_terms": ";".join(group["term"]),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["tissue", "family_narrative_decision", "representative_rank"]
    )


def plot_family_review(families: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 11.2))
    for ax, tissue in zip(axes.flat, ("thymus", "skin", "liver", "soleus")):
        subset = families.loc[families["tissue"].eq(tissue)].copy()
        priority = {
            "existing_core_axis": 0,
            "existing_support_axis": 1,
            "add_aligned_context": 2,
            "add_interpretive_context": 2,
            "revise_existing_claim": 3,
            "existing_context": 3,
            "supplement_only": 4,
            "exclude": 5,
        }
        subset["priority"] = subset["family_narrative_decision"].map(priority)
        subset = subset.sort_values(["priority", "representative_rank"], ascending=[False, False])
        y = np.arange(len(subset))
        for y_value, row in zip(y, subset.itertuples()):
            color = ROLE_COLORS[row.family_literature_role]
            ax.plot(
                [row.minimum_family_effect, row.maximum_family_effect],
                [y_value, y_value],
                color=color,
                alpha=0.55,
                linewidth=3,
                solid_capstyle="round",
            )
            marker = "x" if row.family_narrative_decision == "exclude" else "o"
            face = "none" if row.family_narrative_decision in {"supplement_only", "exclude"} else color
            scatter_kwargs = {
                "s": 45 + 16 * row.n_candidate_pathways,
                "marker": marker,
                "linewidths": 1.4,
                "zorder": 3,
            }
            if marker == "x":
                scatter_kwargs["color"] = color
            else:
                scatter_kwargs["facecolors"] = face
                scatter_kwargs["edgecolors"] = color
            ax.scatter(row.representative_effect, y_value, **scatter_kwargs)
        labels = ["\n".join(textwrap.wrap(name, 29)) for name in subset["process_family"]]
        ax.set_yticks(y, labels, fontsize=8.2)
        ax.axvline(0, color="#333333", linewidth=0.8)
        ax.grid(axis="x", color="#dddddd", linewidth=0.7)
        ax.set_axisbelow(True)
        ax.set_title(
            f"{tissue.title()} | {len(subset)} nonredundant families",
            fontsize=12,
            fontweight="bold",
        )
        ax.set_xlabel("Decoder-oriented FLT - GC pathway score", fontsize=9)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)

    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color,
               markeredgecolor=color, label=label, markersize=7)
        for label, color in (
            ("Literature-aligned", ROLE_COLORS["aligned"]),
            ("Complementary", ROLE_COLORS["complementary"]),
            ("Context-sensitive", ROLE_COLORS["context_sensitive"]),
            ("Not interpretable", ROLE_COLORS["not_interpretable"]),
        )
    ]
    fig.legend(
        handles=legend,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.943),
        ncol=4,
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(
        "Expanded pathway review after Reactome-family consolidation",
        fontsize=15,
        fontweight="bold",
        y=0.982,
    )
    fig.text(
        0.5,
        0.015,
        "Dots show a literature-reviewed representative; horizontal segments span candidate pathways in the same family. "
        "Dot area reflects family size. Hollow points are supplementary-only; x marks excluded broad or tissue-incongruent families. "
        "Scales are tissue-specific.",
        ha="center",
        fontsize=8.5,
    )
    fig.subplots_adjust(left=0.21, right=0.985, top=0.89, bottom=0.08, hspace=0.34, wspace=0.42)
    for suffix, kwargs in (("png", {"dpi": 300}), ("pdf", {})):
        fig.savefig(
            FIGURE_DIR / f"figure_s3_expanded_family_review.{suffix}",
            bbox_inches="tight",
            **kwargs,
        )
    plt.close(fig)


def write_review_summary(pathways: pd.DataFrame, families: pd.DataFrame) -> None:
    lines = [
        "# Expanded pathway-family review",
        "",
        "This audit starts from every active program, selects the top within-tissue decile by absolute study-balanced effect, adds directionally stable programs through rank 40, and retains all pathways prespecified for the main figures. It then consolidates overlapping Reactome terms into manually reviewed process families.",
        "",
        f"The expanded set contains {len(pathways)} pathway records and {len(families)} nonredundant tissue-family records. It is a review queue and evidence synthesis, not a new significance threshold.",
        "",
        "## Main conclusions",
        "",
        "- Thymus: the expanded terms reinforce one lower proliferation and genome-maintenance axis plus the existing lower niche and higher innate-stromal response. Lower apoptosis and stress scores favor a proliferation-led interpretation in the longer-duration cohorts but remain timing-sensitive.",
        "- Skin: higher cutaneous-muscle and innate-phagocytic programs recover known full-thickness-skin compartments. Lower communication and epidermal programs remain coherent. Opposite nested glutathione and phase-II scores mean that detoxification cannot be assigned one uniform direction.",
        "- Liver: additional MHC-I, complement, neutrophil, and interleukin terms strengthen the lower immune-effector family. Opposite nested TLR4 and broad-TLR scores are not interpreted as selective TLR4 activation.",
        "- Soleus: lower PI3K-AKT, WNT, and VEGF scores add literature-aligned trophic context to the matrix-disassembly hypothesis, but the PI3K-AKT magnitude attenuates after OSD-714 exclusion. Cornified-envelope and keratin terms are excluded as tissue-incongruent.",
        "",
        "No newly reviewed complementary family is sufficiently distinct and mission-stable to replace the manuscript's central tissue narratives. The review adds aligned context, strengthens existing families, and narrows claims where nested pathways conflict.",
        "",
        "## Family decisions",
        "",
        "| Tissue | Process family | Role | Decision | Representative | Paths | FLT-GC range |",
        "| --- | --- | --- | --- | --- | ---: | ---: |",
    ]
    for row in families.sort_values(["tissue", "representative_rank"]).itertuples():
        lines.append(
            f"| {row.tissue.title()} | {row.process_family} | {row.family_literature_role.replace('_', ' ')} | "
            f"{row.family_narrative_decision.replace('_', ' ')} | {row.representative_label} | "
            f"{row.n_candidate_pathways} | {row.minimum_family_effect:+.3f} to {row.maximum_family_effect:+.3f} |"
        )
    (PAPER_DIR / "expanded_pathway_family_review.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> None:
    screen = pd.read_csv(SOURCE_DIR / "table_s9_systematic_pathway_screen.tsv", sep="\t")
    gmt = parse_gmt()
    candidates = candidate_screen(screen)
    candidates = add_gene_support(candidates, gmt)
    candidates = add_sensitivity(candidates)
    candidates = add_overlap_flags(candidates, gmt)
    candidates = annotate_families(candidates)
    correlations = centered_score_correlations(candidates)
    families = summarize_families(candidates, correlations, gmt)

    pathway_columns = [
        "tissue",
        "term",
        "screen_display_label",
        "candidate_reason",
        "absolute_effect_rank_active",
        "within_tissue_magnitude_percentile",
        "mean_accession_effect",
        "primary_direction_accessions",
        "primary_direction_agreement",
        "sensitivity_analysis",
        "sensitivity_effect",
        "sensitivity_direction_preserved",
        "sensitivity_magnitude_ratio",
        "meta_fdr",
        "maximum_leave_one_out_fdr",
        "pathway_genes_measured",
        "pathway_genes_pooled_fdr_lt_005",
        "pathway_gene_support_fraction",
        "significant_genes_same_direction",
        "significant_genes_opposite_direction",
        "opposite_nested_program",
        "opposite_nested_terms",
        "maximum_opposite_nested_overlap_coefficient",
        "family_id",
        "process_family",
        "family_literature_role",
        "family_narrative_decision",
        "family_representative",
        "pathway_review_disposition",
        "family_synthesis",
        "family_supporting_references",
        "family_caution",
        "curated_for_main_figures",
        "evidence_role",
    ]
    candidates[pathway_columns].sort_values(
        ["tissue", "absolute_effect_rank_active"]
    ).to_csv(SOURCE_DIR / "table_s10_expanded_pathway_review.tsv", sep="\t", index=False)
    families.to_csv(
        SOURCE_DIR / "table_s11_nonredundant_pathway_families.tsv", sep="\t", index=False
    )
    plot_family_review(families)
    write_review_summary(candidates, families)
    print(
        f"Expanded review: pathways={len(candidates)}, families={len(families)}, "
        f"opposite_nested={int(candidates['opposite_nested_program'].sum())}"
    )


if __name__ == "__main__":
    run()
