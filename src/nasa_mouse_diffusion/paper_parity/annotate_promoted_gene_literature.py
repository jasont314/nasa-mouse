"""Build literature annotations for all synthetic-informed gene associations."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal


ROOT = Path(__file__).resolve().parents[3]
SOURCE_DIR = ROOT / "paper" / "synthetic_guided_spaceflight" / "source_data"
SYNTHETIC_INFORMED_INPUT = (
    SOURCE_DIR / "table_s10_synthetic_informed_bh_fdr_genes.tsv"
)
ANNOTATION_OUTPUT = SOURCE_DIR / "table_s16_promoted_gene_literature_annotations.tsv"
SOURCE_OUTPUT = SOURCE_DIR / "table_s17_promoted_gene_literature_sources.tsv"
SEARCH_DATE = "2026-08-03"

ALLOWED_CLASSIFICATIONS = {
    "aligning",
    "complementary",
    "ambiguous",
    "unmatched",
}
EXPECTED_COUNTS = {
    "aligning": 22,
    "complementary": 19,
    "ambiguous": 4,
    "unmatched": 4,
}
INTERPRETIVE_ROLE_BY_CLASSIFICATION = {
    "aligning": "recovery_of_prior_evidence",
    "complementary": "mechanistic_or_process_hypothesis_extension",
    "ambiguous": "context_dependent_or_mixed_prior_evidence",
    "unmatched": "literature_unmatched_candidate",
}
SELECTION_STATUS = {
    "synthetic_promoted": "promoted",
    "reinforced_real_and_synthetic": "reinforced",
}


@dataclass(frozen=True)
class LiteratureSource:
    source_id: str
    citation: str
    year: int
    doi: str
    url: str
    evidence_role: str
    data_relationship: str


@dataclass(frozen=True)
class LiteratureAnnotation:
    analysis_scope: str
    tissue: str
    symbol: str
    literature_classification: str
    evidence_scope: str
    evidence_relationship: str
    source_ids: tuple[str, ...]
    literature_summary: str
    interpretation: str


SOURCES = (
    LiteratureSource(
        "horie_2019_thymus",
        "Horie K, Kato T, Kudo T, et al. Impact of spaceflight on the murine thymus and mitigation by exposure to artificial gravity during spaceflight. Scientific Reports. 2019;9:19866.",
        2019,
        "10.1038/s41598-019-56432-9",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC6934594/",
        "Two-mission mouse thymus RNA-seq with exact cyclin results and process-level cell-cycle evidence.",
        "Published OSDR-related source cohorts; independence from the present aggregate is not assumed.",
    ),
    LiteratureSource(
        "gridley_2009_thymus_spleen",
        "Gridley DS, Slater JM, Luo-Owen X, et al. Spaceflight effects on T lymphocyte distribution, function and gene expression. Journal of Applied Physiology. 2009;106:194-202.",
        2009,
        "10.1152/japplphysiol.91126.2008",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC2636934/",
        "STS-118 thymus gene panel and spleen immune phenotyping; reported flight-higher Birc5.",
        "Independent shuttle mission and targeted expression platform relative to the present bulk RNA-seq aggregate.",
    ),
    LiteratureSource(
        "allen_2009_muscle",
        "Allen DL, Bandstra ER, Harrison BC, et al. Effects of spaceflight on murine skeletal muscle gene expression. Journal of Applied Physiology. 2009;106:582-595.",
        2009,
        "10.1152/japplphysiol.90780.2008",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC2644242/",
        "STS-108 gastrocnemius microarray and PCR study; reported higher Nfkbia and C/EBP-delta after flight.",
        "Independent shuttle mission and microarray/PCR platform relative to the present RNA-seq aggregate.",
    ),
    LiteratureSource(
        "li_2023_muscle",
        "Li K, Desai R, Scott RT, et al. Explainable machine learning identifies multi-omics signatures of muscle response to spaceflight in mice. npj Microgravity. 2023;9:95.",
        2023,
        "10.1038/s41526-023-00337-5",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC10719374/",
        "OSDR muscle multi-omics study linking spaceflight to calcium handling and listing Klhl21 among tibialis-anterior model features.",
        "Uses public OSDR muscle cohorts; treated as related evidence rather than independent validation.",
    ),
    LiteratureSource(
        "new_2003_mapkapk5",
        "New L, Jiang Y, Han J. Regulation of PRAK subcellular location by p38 MAP kinases. Molecular Biology of the Cell. 2003;14:2603-2616.",
        2003,
        "10.1091/mbc.E02-08-0538",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC194907/",
        "Primary mechanistic evidence that MAPKAPK5/PRAK functions downstream of stress-responsive p38 signaling.",
        "Mechanistic context only; not a spaceflight or skeletal-muscle replication study.",
    ),
    LiteratureSource(
        "yao_2018_reep5",
        "Yao L, Xie D, Geng L, et al. REEP5 acts as a sarcoplasmic reticulum membrane sculptor to modulate cardiac function. Journal of the American Heart Association. 2018;7:e007205.",
        2018,
        "10.1161/JAHA.117.007205",
        "https://pubmed.ncbi.nlm.nih.gov/29431104/",
        "Primary evidence connecting REEP5 to sarcoplasmic-reticulum architecture and calcium handling.",
        "Mechanistic cardiac-muscle context only; not a spaceflight skeletal-muscle replication.",
    ),
    LiteratureSource(
        "vitadello_2020_unloading",
        "Vitadello M, Sorge M, Percivalle E, et al. Loss of melusin is a novel, neuronal NO synthase/FoxO3-independent master switch of unloading-induced muscle atrophy. Journal of Cachexia, Sarcopenia and Muscle. 2020;11:802-819.",
        2020,
        "10.1002/jcsm.12546",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC7296270/",
        "Primary evidence that an integrin-associated mechanotransduction protein responds to unloading and modifies muscle atrophy.",
        "Ground unloading and mechanistic context only; ITGB5 itself was not tested.",
    ),
    LiteratureSource(
        "liu_2019_fhl2",
        "Liu Z, Han S, Wang Y, et al. The LIM-only protein FHL2 is involved in autophagy to regulate the development of skeletal muscle cell. International Journal of Biological Sciences. 2019;15:838-846.",
        2019,
        "10.7150/ijbs.31371",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC6429013/",
        "Primary evidence linking FHL2 to myogenic differentiation and autophagy.",
        "Mechanistic muscle context only; not a spaceflight replication study.",
    ),
    LiteratureSource(
        "cope_2024_skin",
        "Cope H, Elsborg J, Demharter S, et al. Transcriptomics analysis reveals molecular alterations underpinning spaceflight dermatology. Communications Medicine. 2024;4:106.",
        2024,
        "10.1038/s43856-024-00532-9",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC11166967/",
        "Multi-study OSDR skin analysis reporting immune, barrier, mitochondrial, collagen, and DNA-repair responses.",
        "Uses public OSDR skin cohorts; treated as related evidence rather than independent validation.",
    ),
    LiteratureSource(
        "dong_2004_plscr1",
        "Dong B, Zhou Q, Zhao J, et al. Phospholipid scramblase 1 potentiates the antiviral activity of interferon. Journal of Virology. 2004;78:8983-8993.",
        2004,
        "10.1128/JVI.78.17.8983-8993.2004",
        "https://pubmed.ncbi.nlm.nih.gov/15308695/",
        "Primary evidence that PLSCR1 is interferon inducible and amplifies a subset of interferon-stimulated genes.",
        "Mechanistic immune context only; not a spaceflight skin replication.",
    ),
    LiteratureSource(
        "hammond_2018_kidney",
        "Hammond TG, Allen PL, Birdsall HH. Effects of space flight on mouse liver versus kidney: gene pathway analyses. International Journal of Molecular Sciences. 2018;19:4106.",
        2018,
        "10.3390/ijms19124106",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC6321533/",
        "Mouse flight kidney analysis reporting altered protein-kinase regulation and hormone or peptide responses.",
        "Published spaceflight expression cohort; treated as process context rather than an independent INPP4B replication.",
    ),
    LiteratureSource(
        "kofuji_2015_inpp4b",
        "Kofuji S, Kimura H, Nakanishi H, et al. INPP4B is a PtdIns(3,4,5)P3 phosphatase that can act as a tumor suppressor. Cancer Discovery. 2015;5:730-739.",
        2015,
        "10.1158/2159-8290.CD-14-1329",
        "https://pubmed.ncbi.nlm.nih.gov/25883023/",
        "Primary biochemical evidence that INPP4B dephosphorylates PIP3 and restrains downstream AKT signaling.",
        "Mechanistic signaling context only; not a spaceflight kidney replication.",
    ),
    LiteratureSource(
        "jeong_2024_rai14",
        "Jeong W, Kwon H, Park SK, Lee IS, Jho EH. Retinoic acid-induced protein 14 links mechanical forces to Hippo signaling. EMBO Reports. 2024;25:4033-4061.",
        2024,
        "10.1038/s44319-024-00228-0",
        "https://pubmed.ncbi.nlm.nih.gov/39160347/",
        "Primary evidence that RAI14 senses mechanical forces through F-actin and regulates Hippo signaling.",
        "Mechanistic context only; not a spleen or spaceflight replication.",
    ),
    LiteratureSource(
        "hayashizaki_2016_myl9",
        "Hayashizaki K, Kimura MY, Tokoyoda K, et al. Myosin light chains 9 and 12 are functional ligands for CD69 that regulate airway inflammation. Science Immunology. 2016;1:eaaf9154.",
        2016,
        "10.1126/sciimmunol.aaf9154",
        "https://pubmed.ncbi.nlm.nih.gov/28783682/",
        "Primary evidence that platelet-derived MYL9 structures recruit CD69-positive inflammatory cells.",
        "Mechanistic inflammatory context only; not a spleen or spaceflight replication.",
    ),
    LiteratureSource(
        "fearnley_2019_ptprk",
        "Fearnley GW, Young KA, Edgar JR, et al. The homophilic receptor PTPRK selectively dephosphorylates multiple junctional regulators to promote cell-cell adhesion. eLife. 2019;8:e44597.",
        2019,
        "10.7554/eLife.44597",
        "https://pubmed.ncbi.nlm.nih.gov/30924770/",
        "Primary evidence linking PTPRK to junctional regulation and cell-cell adhesion.",
        "Mechanistic adhesion context only; not a spleen or spaceflight replication.",
    ),
    LiteratureSource(
        "mathyk_2024_multi_organ",
        "Mathyk BA, Tabetah M, Karim R, et al. Spaceflight induces changes in gene expression profiles linked to insulin and estrogen. Communications Biology. 2024;7:674.",
        2024,
        "10.1038/s42003-023-05213-2",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC11166981/",
        "Multi-organ OSDR analysis reporting 61 adrenal differentially expressed genes and broad endocrine or metabolic responses.",
        "Uses public OSDR cohorts; used to define the searched adrenal context, not as independent validation.",
    ),
    LiteratureSource(
        "keenan_2025_hsd17b11",
        "Keenan SN, Suriani ND, Fidelito G, et al. HSD17B11 regulates PLIN5-ATGL mediated lipolysis, but not hepatic lipid metabolism in mice. Journal of Lipid Research. 2025;66:100943.",
        2025,
        "10.1016/j.jlr.2025.100943",
        "https://pubmed.ncbi.nlm.nih.gov/41238190/",
        "Primary evidence connecting HSD17B11 to lipid droplets, steroid metabolism, and regulated lipolysis.",
        "Mechanistic metabolic context only; not a thymus or spaceflight replication study.",
    ),
    LiteratureSource(
        "shi_2026_etv1",
        "Shi Y, Wang S, Yan Y, et al. ETV1 drives CD4+ T cell-mediated intestinal inflammation in inflammatory bowel disease through amino acid transporter Slc7a5. Advanced Science. 2026;13:e11595.",
        2026,
        "10.1002/advs.202511595",
        "https://pubmed.ncbi.nlm.nih.gov/41347630/",
        "Primary evidence that ETV1 regulates CD4+ T-cell activation, proliferation, differentiation, and amino-acid uptake.",
        "Mechanistic T-cell context only; not a thymus or spaceflight replication study.",
    ),
    LiteratureSource(
        "kitamura_2011_psmb8",
        "Kitamura A, Maekawa Y, Uehara H, et al. A mutation in the immunoproteasome subunit PSMB8 causes autoinflammation and lipodystrophy in humans. Journal of Clinical Investigation. 2011;121:4150-4160.",
        2011,
        "10.1172/JCI58414",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC3195477/",
        "Primary evidence connecting interferon-inducible PSMB8 to immunoproteasome function, inflammation, and adipocyte homeostasis.",
        "Mechanistic immune and adipocyte context only; not an adrenal or spaceflight replication study.",
    ),
    LiteratureSource(
        "gambara_2017_longissimus",
        "Gambara G, Salanova M, Ciciliot S, et al. Microgravity-induced transcriptome adaptation in mouse paraspinal longissimus dorsi muscle highlights insulin resistance-linked genes. Frontiers in Physiology. 2017;8:279.",
        2017,
        "10.3389/fphys.2017.00279",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC5418220/",
        "Mouse flight-muscle transcriptomics with direct Cebpd and Arid5b results and qPCR support for Sesn1 and Prkcd.",
        "Published BION-M1 muscle cohort that may overlap the public aggregate; used as literature alignment, not independent validation.",
    ),
    LiteratureSource(
        "gambara_2017_soleus",
        "Gambara G, Salanova M, Ciciliot S, et al. Gene expression profiling in slow-type calf soleus muscle of 30 days space-flown mice. PLOS ONE. 2017;12:e0169314.",
        2017,
        "10.1371/journal.pone.0169314",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC5226702/",
        "Thirty-day mouse flight soleus transcriptomics reporting a slow-to-fast shift and broad oxidative-metabolism remodeling.",
        "Published BION-M1 soleus cohort that may overlap the public aggregate; used as process-level alignment, not independent validation.",
    ),
    LiteratureSource(
        "oommen_2024_muscle",
        "Oommen AM, Stafford P, Joshi L. Profiling muscle transcriptome in mice exposed to microgravity using gene set enrichment analysis. npj Microgravity. 2024;10:94.",
        2024,
        "10.1038/s41526-024-00434-z",
        "https://www.nature.com/articles/s41526-024-00434-z",
        "Public mouse flight-muscle synthesis reporting Sox4, Arid5b, Cdkn1a, and glycosylation-related responses.",
        "Reanalysis of public GeneLab muscle studies; overlap with the present aggregate is expected.",
    ),
    LiteratureSource(
        "ohira_2021_soleus",
        "Ohira T, Ino Y, Kimura Y, et al. Effects of microgravity exposure and fructo-oligosaccharide ingestion on the proteome of soleus and extensor digitorum longus muscles in developing mice. npj Microgravity. 2021;7:34.",
        2021,
        "10.1038/s41526-021-00164-6",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8448765/",
        "Mouse flight proteomics reporting lower ECH1 and DECR1 in soleus under microgravity and reduced oxidative-metabolism proteins.",
        "Independent assay-level support; the flight cohort may be represented in related public data.",
    ),
    LiteratureSource(
        "rosa_caldwell_2021_unloading",
        "Rosa-Caldwell ME, Lim S, Haynie WS, et al. Mitochondrial aberrations during the progression of disuse atrophy differentially affect male and female mice. Journal of Cachexia, Sarcopenia and Muscle. 2021;12:2056-2068.",
        2021,
        "10.1002/jcsm.12809",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8718086/",
        "Mouse hindlimb-unloading study showing sex-dependent Bnip3 responses in soleus.",
        "Independent ground analog; used to identify context dependence, not as spaceflight replication.",
    ),
    LiteratureSource(
        "chen_2021_retina",
        "Chen Z, Stanbouly S, Nishiyama NC, et al. Spaceflight decelerates the epigenetic clock orchestrated with a global alteration in DNA methylome and transcriptome in the mouse retina. Precision Clinical Medicine. 2021;4:93-108.",
        2021,
        "10.1093/pcmedi/pbab012",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8220224/",
        "Mouse flight retina multi-omics reporting altered proliferation, cell-cycle, and DNA-repair programs.",
        "Published public flight-eye cohort; tissue and sample overlap with the present eye aggregate may occur.",
    ),
    LiteratureSource(
        "maerki_2009_klhl21",
        "Maerki S, Olma MH, Staubli T, et al. The Cul3-KLHL21 E3 ubiquitin ligase targets Aurora B to midzone microtubules in anaphase and is required for cytokinesis. Journal of Cell Biology. 2009;187:791-800.",
        2009,
        "10.1083/jcb.200906117",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC2806313/",
        "Primary mechanistic evidence that KLHL21 is required for chromosome alignment and cytokinesis.",
        "Mechanistic context only; not a flight-eye replication.",
    ),
    LiteratureSource(
        "dacierno_2022_slc37a4",
        "D'Acierno M, Resaz R, Iervolino A, et al. Dapagliflozin prevents kidney glycogen accumulation and improves renal proximal tubule cell functions in a mouse model of glycogen storage disease type 1b. Journal of the American Society of Nephrology. 2022;33:1864-1875.",
        2022,
        "10.1681/ASN.2021070935",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC9528317/",
        "Mouse evidence that loss of the SLC37A4 glucose-6-phosphate transporter disrupts proximal-tubule metabolism and function.",
        "Mechanistic kidney context only; not a spaceflight replication.",
    ),
    LiteratureSource(
        "wiltshire_2002_sh3bp5",
        "Wiltshire C, Matsushita M, Tsukada S, Gillespie DAF, May GHW. A new c-Jun N-terminal kinase-interacting protein, Sab (SH3BP5), associates with mitochondria. Biochemical Journal. 2002;367:577-585.",
        2002,
        "10.1042/BJ20020553",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC1222945/",
        "Primary evidence that SH3BP5/SAB associates with mitochondria and binds stress-responsive JNK.",
        "Mechanistic context only; not a spaceflight-muscle replication.",
    ),
    LiteratureSource(
        "ramasamy_2016_tle1",
        "Ramasamy S, Saez B, Mukhopadhyay S, et al. Tle1 tumor suppressor negatively regulates inflammation in vivo and modulates NF-kappaB inflammatory pathway. Proceedings of the National Academy of Sciences USA. 2016;113:1871-1876.",
        2016,
        "10.1073/pnas.1511380113",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC4763742/",
        "Mouse knockout evidence that TLE1 counterregulates inflammatory NF-kappaB signaling.",
        "Mechanistic context only; not a spaceflight-muscle replication.",
    ),
    LiteratureSource(
        "han_2023_spleen_thymus",
        "Han Y, Shi S, Liu S, Gu X. Effects of spaceflight on the spleen and thymus of mice: gene pathway analysis and immune infiltration analysis. Mathematical Biosciences and Engineering. 2023;20:8531-8545.",
        2023,
        "10.3934/mbe.2023374",
        "https://pubmed.ncbi.nlm.nih.gov/37161210/",
        "Public mouse flight-data reanalysis reporting spleen and thymus immune, platelet, and infiltration responses.",
        "Reanalysis of public mouse flight studies; overlap with the present aggregate is expected.",
    ),
    LiteratureSource(
        "li_2021_loxl1",
        "Li Y, Wu B, Zou X. Mass cytometry and transcriptomic profiling reveal body-wide pathology induced by Loxl1 deficiency. Cell Proliferation. 2021;54:e13077.",
        2021,
        "10.1111/cpr.13077",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8249785/",
        "Mouse multi-organ evidence connecting LOXL1-dependent extracellular matrix homeostasis to splenic immune and cell-cycle state.",
        "Mechanistic spleen context only; not a spaceflight replication.",
    ),
    LiteratureSource(
        "yamashita_2003_st3gal5",
        "Yamashita T, Hashiramoto A, Haluzik M, et al. Enhanced insulin sensitivity in mice lacking ganglioside GM3. Proceedings of the National Academy of Sciences USA. 2003;100:3445-3449.",
        2003,
        "10.1073/pnas.0635898100",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC152312/",
        "Mouse evidence that loss of ST3GAL5-dependent GM3 synthesis enhances insulin-receptor phosphorylation in skeletal muscle.",
        "Mechanistic muscle context only; not a tibialis-anterior spaceflight replication.",
    ),
    LiteratureSource(
        "taylor_2002_rat_muscle",
        "Taylor WE, Bhasin S, Lalani R, Datta A, Gonzalez-Cadavid NF. Alteration of gene expression profiles in skeletal muscle of rats exposed to microgravity during a spaceflight. Journal of Gravitational Physiology. 2002;9:61-70.",
        2002,
        "not_assigned",
        "https://ntrs.nasa.gov/citations/20040087468",
        "Rat STS-90 tibialis-anterior and gastrocnemius profiling reporting lower p21/Cdkn1a-related proliferation signaling after flight.",
        "Independent species, mission, and expression platform; used as opposing directional context.",
    ),
)


def _annotation(
    tissue: str,
    symbol: str,
    classification: str,
    scope: str,
    relationship: str,
    source_ids: tuple[str, ...],
    summary: str,
    interpretation: str,
) -> LiteratureAnnotation:
    muscle_groups = {
        "extensor_digitorum_longus",
        "gastrocnemius",
        "quadriceps",
        "soleus",
        "tibialis_anterior",
    }
    return LiteratureAnnotation(
        "skeletal_muscle_group" if tissue in muscle_groups else "canonical_tissue",
        tissue,
        symbol,
        classification,
        scope,
        relationship,
        source_ids,
        summary,
        interpretation,
    )


THYMUS_PROCESS_SUMMARY = (
    "Two-mission mouse thymus RNA-seq reported flight-lower cell-cycle and "
    "chromosome-organization programs. This gene fits that mitotic or replication "
    "program, but the paper did not report it as an exact directional replication."
)
THYMUS_PROCESS_INTERPRETATION = (
    "Aligns at the same-tissue process level; it is not an exact prior gene replication."
)


ANNOTATIONS = (
    _annotation(
        "adrenal_gland",
        "Psmb8",
        "unmatched",
        "targeted_search_no_adrenal_specific_match",
        "Mechanistic support exists, but no adrenal spaceflight match was found in the targeted search.",
        ("mathyk_2024_multi_organ", "kitamura_2011_psmb8"),
        "PSMB8 is an interferon-inducible immunoproteasome subunit with roles in inflammation and adipocyte homeostasis. Published mouse adrenal analysis reported a small endocrine and metabolic response, but the targeted search found no prior adrenal spaceflight Psmb8 association.",
        "A biologically plausible immune, proteostasis, or composition candidate that remains unmatched in the searched adrenal spaceflight literature.",
    ),
    _annotation(
        "kidney",
        "Inpp4b",
        "complementary",
        "related_spaceflight_process_and_gene_mechanism",
        "Contextual support only; no same-gene kidney spaceflight replication was found.",
        ("hammond_2018_kidney", "kofuji_2015_inpp4b"),
        "Spaceflight kidney work implicated protein-kinase and hormone-response pathways, while INPP4B directly regulates PIP3-AKT signaling. The targeted search found no prior same-gene flight-kidney direction.",
        "Adds a focused phosphoinositide-signaling candidate to prior renal pathway findings.",
    ),
    _annotation(
        "skeletal_muscle",
        "Klhl21",
        "complementary",
        "same_gene_other_muscle_context_direction_unreported",
        "Public OSDR cohorts may overlap; this is related evidence, not independent validation.",
        ("li_2023_muscle",),
        "KLHL21 appeared among tibialis-anterior spaceflight model features, but its direction was not established and the present association is for pooled skeletal muscle.",
        "Extends a prior muscle-model feature to a directionally defined pooled-muscle association.",
    ),
    _annotation(
        "skeletal_muscle",
        "Mapkapk5",
        "complementary",
        "related_gene_mechanism",
        "Mechanistic context only; no direct spaceflight muscle replication was found.",
        ("new_2003_mapkapk5", "allen_2009_muscle"),
        "MAPKAPK5/PRAK functions downstream of stress-responsive p38 signaling, while prior flight muscle work reported broad stress-pathway remodeling. No prior same-gene flight direction was found.",
        "Nominates a stress-kinase component within established unloading-responsive muscle biology.",
    ),
    _annotation(
        "skeletal_muscle",
        "Reep5",
        "complementary",
        "related_gene_mechanism",
        "Mechanistic and OSDR context only; no direct spaceflight skeletal-muscle replication was found.",
        ("yao_2018_reep5", "li_2023_muscle"),
        "REEP5 shapes sarcoplasmic-reticulum membranes and calcium handling in muscle, while OSDR muscle studies implicate calcium and SERCA physiology. No prior flight-muscle REEP5 direction was found.",
        "Adds a membrane-architecture candidate to established calcium-handling responses.",
    ),
    _annotation(
        "skeletal_muscle",
        "Itgb5",
        "complementary",
        "related_unloading_mechanotransduction",
        "Ground-unloading mechanism only; ITGB5 itself was not replicated.",
        ("vitadello_2020_unloading", "allen_2009_muscle"),
        "An integrin-associated mechanotransduction protein responds early to unloading and modifies atrophy, but the targeted search found no direct ITGB5 spaceflight-muscle report.",
        "Adds a specific integrin subunit to established unloading-sensitive adhesion and mechanotransduction biology.",
    ),
    _annotation(
        "skin",
        "Plscr1",
        "complementary",
        "related_spaceflight_process_and_gene_mechanism",
        "OSDR skin cohorts may overlap; the gene mechanism is contextual rather than a replication.",
        ("cope_2024_skin", "dong_2004_plscr1"),
        "Mouse spaceflight skin studies reported immune and damage-response changes, and PLSCR1 is interferon inducible and amplifies interferon-stimulated transcription. No prior directional skin-flight match was found.",
        "Adds an interferon-linked candidate to the published skin immune-response context.",
    ),
    _annotation(
        "spleen",
        "Rai14",
        "complementary",
        "related_spaceflight_process_and_gene_mechanism",
        "Mechanistic context only; no same-gene spleen spaceflight replication was found.",
        ("gridley_2009_thymus_spleen", "jeong_2024_rai14"),
        "Spaceflight alters splenic immune composition and function, while RAI14 links F-actin-dependent force sensing to Hippo signaling. No prior spleen-flight RAI14 direction was found.",
        "Adds a mechanosensing candidate to the observed splenic structural and immune response.",
    ),
    _annotation(
        "spleen",
        "Myl9",
        "complementary",
        "related_spaceflight_process_and_gene_mechanism",
        "Mechanistic context only; no same-gene spleen spaceflight replication was found.",
        ("gridley_2009_thymus_spleen", "hayashizaki_2016_myl9"),
        "Spaceflight alters splenic immune-cell distributions, and platelet-derived MYL9 structures can recruit CD69-positive inflammatory cells. No prior spleen-flight MYL9 direction was found.",
        "Suggests a platelet-actomyosin link to the broader splenic immune response.",
    ),
    _annotation(
        "spleen",
        "Ptprk",
        "complementary",
        "related_spaceflight_process_and_gene_mechanism",
        "Mechanistic context only; no same-gene spleen spaceflight replication was found.",
        ("gridley_2009_thymus_spleen", "fearnley_2019_ptprk"),
        "Spaceflight alters splenic immune composition, while PTPRK regulates junctional proteins and cell-cell adhesion. No prior spleen-flight PTPRK direction was found.",
        "Adds a junctional-regulation candidate to the splenic structural hypothesis.",
    ),
    *(
        _annotation(
            "thymus",
            symbol,
            "aligning",
            "same_tissue_process_same_direction",
            "Published OSDR-related cohorts may overlap; this is process-level corroboration, not independent validation.",
            ("horie_2019_thymus",),
            THYMUS_PROCESS_SUMMARY,
            THYMUS_PROCESS_INTERPRETATION,
        )
        for symbol in ("Nusap1", "Stmn1", "Cdk1", "Top2a", "Aurka", "Kif20a", "Pcna", "Ccnf")
    ),
    _annotation(
        "thymus",
        "Birc5",
        "ambiguous",
        "mixed_exact_gene_and_process_evidence",
        "Includes an independent direct conflict and OSDR-related process-level agreement.",
        ("gridley_2009_thymus_spleen", "horie_2019_thymus"),
        "STS-118 targeted expression reported flight-higher Birc5, opposite the present flight-lower effect. Longer ISS RNA-seq later reported a broad flight-lower thymus cell-cycle program.",
        "Prior evidence is mixed across mission duration and assay, so the association is neither solely aligning nor solely contradictory.",
    ),
    *(
        _annotation(
            "thymus",
            symbol,
            "aligning",
            "direct_same_gene_same_tissue_same_direction",
            "Published OSDR-related source cohorts may overlap; this is an exact literature match but not independent validation.",
            ("horie_2019_thymus",),
            f"Two-mission mouse thymus RNA-seq reported flight-lower {symbol} expression, matching the present gene, tissue, and direction.",
            "Exact same-gene, same-tissue, same-direction literature alignment.",
        )
        for symbol in ("Ccnb2", "Ccne2")
    ),
    _annotation(
        "thymus",
        "Hsd17b11",
        "complementary",
        "related_lipid_mechanism_and_thymus_composition",
        "Mechanistic and tissue-composition context only; no direct thymus spaceflight match was found.",
        ("keenan_2025_hsd17b11", "horie_2019_thymus"),
        "HSD17B11 localizes to lipid droplets and supports regulated lipolysis in human cells, while mouse flight thymus work reported lower proliferation and expression shifts consistent with altered cellular composition. No prior directional thymus-flight Hsd17b11 result was found.",
        "A plausible lipid-handling or tissue-composition marker, not evidence that HSD17B11 drives the thymus response.",
    ),
    _annotation(
        "thymus",
        "Etv1",
        "complementary",
        "related_t_cell_mechanism_and_thymus_composition",
        "Mechanistic and tissue-composition context only; no direct thymus spaceflight match was found.",
        ("shi_2026_etv1", "horie_2019_thymus"),
        "ETV1 regulates CD4+ T-cell activation and proliferation, while mouse flight thymus work reported lower proliferation and expression shifts consistent with altered cellular composition. No prior directional thymus-flight Etv1 result was found.",
        "A plausible T-cell-state or tissue-composition marker, not a direct spaceflight mechanism.",
    ),
    _annotation(
        "gastrocnemius",
        "Nfkbia",
        "aligning",
        "direct_same_gene_same_tissue_same_direction",
        "Independent shuttle mission and expression platform.",
        ("allen_2009_muscle",),
        "STS-108 gastrocnemius profiling reported a 2.28-fold flight increase in Nfkbia, matching the present gene, tissue, and direction.",
        "Exact same-gene, same-tissue, same-direction alignment in an independent mission and platform.",
    ),
    _annotation(
        "gastrocnemius",
        "Fhl2",
        "complementary",
        "related_gene_mechanism",
        "Mechanistic context only; no direct gastrocnemius spaceflight direction was found.",
        ("liu_2019_fhl2", "allen_2009_muscle"),
        "FHL2 regulates myogenic differentiation and autophagy, processes relevant to unloading-sensitive muscle remodeling, but no prior directional gastrocnemius flight match was found.",
        "Adds an autophagy and myogenesis candidate to established flight-muscle remodeling.",
    ),
    _annotation(
        "tibialis_anterior",
        "Cebpd",
        "complementary",
        "same_gene_same_direction_other_muscle_group",
        "Independent shuttle study, but the prior result was in gastrocnemius rather than tibialis anterior.",
        ("allen_2009_muscle",),
        "STS-108 gastrocnemius profiling reported higher C/EBP-delta after flight, matching the gene and direction but not the present tibialis-anterior tissue subgroup.",
        "Extends an independent same-direction muscle result to a different anatomical muscle group.",
    ),
    _annotation(
        "adrenal_gland",
        "Tspan4",
        "unmatched",
        "targeted_search_no_adrenal_specific_match",
        "No direct or strong process-level adrenal spaceflight match was found in the targeted search.",
        ("mathyk_2024_multi_organ",),
        "Published mouse adrenal analyses describe endocrine and metabolic responses, but the targeted search found no prior adrenal spaceflight association for Tspan4.",
        "A literature-unmatched adrenal candidate; its statistical association does not by itself establish a mechanism.",
    ),
    _annotation(
        "eye",
        "Klhl21",
        "aligning",
        "same_tissue_process_same_direction",
        "Process-level agreement in a related eye tissue; no exact flight-eye KLHL21 result was found.",
        ("chen_2021_retina", "maerki_2009_klhl21"),
        "Mouse flight retina data reported lower proliferation-related programs, and KLHL21 is required for chromosome alignment and cytokinesis. The present flight-lower eye association fits that direction at the process level.",
        "Aligns with a flight-lower proliferative program in the eye, but is not an exact prior gene replication.",
    ),
    _annotation(
        "kidney",
        "Slc37a4",
        "complementary",
        "related_spaceflight_process_and_gene_mechanism",
        "Mechanistic kidney and spaceflight-metabolism context only; no directional flight-kidney replication was found.",
        ("hammond_2018_kidney", "dacierno_2022_slc37a4"),
        "Spaceflight kidney studies report metabolic remodeling, while mouse SLC37A4 loss disrupts proximal-tubule glucose handling and function. The targeted search found no prior flight-kidney Slc37a4 direction.",
        "Adds an endoplasmic-reticulum glucose-transport candidate to the renal metabolic response.",
    ),
    _annotation(
        "skeletal_muscle",
        "Sox4",
        "aligning",
        "direct_same_gene_same_tissue_same_direction",
        "Public GeneLab muscle studies may overlap; this is literature alignment, not independent validation.",
        ("oommen_2024_muscle",),
        "A public mouse flight-muscle synthesis reported higher Sox4, matching the present pooled-muscle gene and direction.",
        "Same-gene, broad-muscle, same-direction literature alignment.",
    ),
    _annotation(
        "skeletal_muscle",
        "Cebpd",
        "aligning",
        "same_gene_same_direction_other_muscle_context",
        "Published gastrocnemius and longissimus cohorts may partly overlap the public aggregate.",
        ("allen_2009_muscle", "gambara_2017_longissimus"),
        "Mouse flight studies reported higher Cebpd in gastrocnemius and longissimus dorsi, matching the present pooled-muscle direction.",
        "Same-gene and same-direction alignment across multiple anatomical muscle contexts.",
    ),
    _annotation(
        "skeletal_muscle",
        "Sh3bp5",
        "complementary",
        "related_gene_mechanism",
        "Mitochondrial stress mechanism only; no direct flight-muscle match was found.",
        ("wiltshire_2002_sh3bp5", "oommen_2024_muscle"),
        "SH3BP5/SAB localizes to mitochondria and binds stress-responsive JNK, while flight-muscle studies report mitochondrial and stress remodeling. No prior directional flight-muscle Sh3bp5 result was found.",
        "Adds a mitochondrial stress-signaling candidate to established muscle remodeling.",
    ),
    _annotation(
        "skeletal_muscle",
        "Bphl",
        "unmatched",
        "targeted_search_no_muscle_spaceflight_match",
        "No direct or strong process-level flight-muscle match was found in the targeted search.",
        ("allen_2009_muscle", "oommen_2024_muscle"),
        "The targeted search of mouse flight-muscle studies found no prior Bphl association or sufficiently specific process link.",
        "A literature-unmatched pooled-muscle candidate that requires independent functional context.",
    ),
    _annotation(
        "skeletal_muscle",
        "Prkcd",
        "aligning",
        "same_gene_same_direction_other_muscle_context",
        "Published BION-M1 muscle data may overlap the public aggregate.",
        ("gambara_2017_longissimus",),
        "Longissimus dorsi profiling and qPCR reported higher Prkcd after flight, matching the present pooled-muscle direction.",
        "Same-gene and same-direction alignment in another anatomical muscle context.",
    ),
    _annotation(
        "skeletal_muscle",
        "Arid5b",
        "aligning",
        "direct_same_gene_same_tissue_same_direction",
        "Both sources use public or published mouse flight-muscle cohorts that may overlap the present aggregate.",
        ("gambara_2017_longissimus", "oommen_2024_muscle"),
        "Mouse flight-muscle studies reported higher Arid5b, matching the present pooled-muscle gene and direction.",
        "Same-gene, broad-muscle, same-direction literature alignment.",
    ),
    _annotation(
        "skeletal_muscle",
        "Sesn1",
        "aligning",
        "same_gene_same_direction_other_muscle_context",
        "Published BION-M1 muscle data may overlap the public aggregate.",
        ("gambara_2017_longissimus",),
        "Longissimus dorsi profiling and qPCR reported higher Sesn1 after flight, matching the present pooled-muscle direction.",
        "Same-gene and same-direction alignment in another anatomical muscle context.",
    ),
    _annotation(
        "skeletal_muscle",
        "Tle1",
        "complementary",
        "related_gene_mechanism",
        "Inflammatory mechanism only; no direct flight-muscle Tle1 result was found.",
        ("ramasamy_2016_tle1", "oommen_2024_muscle"),
        "Mouse knockout work identifies TLE1 as a counterregulator of NF-kappaB inflammation, a process relevant to flight-muscle remodeling. No prior directional flight-muscle Tle1 result was found.",
        "Adds a transcriptional inflammatory-regulation candidate to pooled muscle.",
    ),
    _annotation(
        "spleen",
        "Loxl1",
        "complementary",
        "related_spaceflight_process_and_gene_mechanism",
        "Mechanistic spleen and public flight-process context only; no direct spleen-flight direction was found.",
        ("han_2023_spleen_thymus", "li_2021_loxl1"),
        "Mouse flight spleen studies implicate immune and platelet processes, while Loxl1 deficiency alters splenic extracellular-matrix, immune, and cell-cycle states. No prior directional spleen-flight Loxl1 result was found.",
        "Adds an extracellular-matrix and immune-organization candidate to the splenic response.",
    ),
    _annotation(
        "thymus",
        "Snx7",
        "unmatched",
        "targeted_search_no_thymus_spaceflight_match",
        "No direct or sufficiently specific process-level thymus spaceflight match was found.",
        ("horie_2019_thymus", "han_2023_spleen_thymus"),
        "The targeted search of mouse flight-thymus studies found no prior Snx7 association and no sufficiently specific mechanism to assign it to the established cell-cycle program.",
        "A literature-unmatched thymus candidate that should remain separate from the coherent mitotic panel.",
    ),
    _annotation(
        "thymus",
        "Ube2c",
        "aligning",
        "same_tissue_process_same_direction",
        "Published OSDR-related cohorts may overlap; this is process-level corroboration, not independent validation.",
        ("horie_2019_thymus",),
        THYMUS_PROCESS_SUMMARY,
        THYMUS_PROCESS_INTERPRETATION,
    ),
    _annotation(
        "thymus",
        "Gmnn",
        "aligning",
        "same_tissue_process_same_direction",
        "Published OSDR-related cohorts may overlap; this is process-level corroboration, not independent validation.",
        ("horie_2019_thymus",),
        THYMUS_PROCESS_SUMMARY,
        THYMUS_PROCESS_INTERPRETATION,
    ),
    _annotation(
        "soleus",
        "Bdh1",
        "aligning",
        "same_tissue_process_same_direction",
        "Published flight-soleus cohorts may overlap; this is process-level rather than exact-gene corroboration.",
        ("gambara_2017_soleus", "ohira_2021_soleus"),
        "Mouse flight soleus studies reported reduced oxidative and lipid-metabolism programs. Flight-lower Bdh1 fits that same-tissue metabolic direction, although no exact prior soleus Bdh1 result was found.",
        "Aligns with the established flight-lower oxidative-metabolism program in soleus.",
    ),
    _annotation(
        "soleus",
        "Ech1",
        "aligning",
        "same_gene_same_tissue_same_direction_cross_assay",
        "Same flight tissue and direction in proteomics; assay-level support rather than independent transcriptomic validation.",
        ("ohira_2021_soleus",),
        "Mouse flight proteomics reported lower ECH1 in soleus under microgravity, matching the present gene, tissue, and direction across assays.",
        "Exact gene/protein, same-tissue, same-direction cross-assay alignment.",
    ),
    _annotation(
        "soleus",
        "Decr1",
        "aligning",
        "same_gene_same_tissue_same_direction_cross_assay",
        "Same flight tissue and direction in proteomics; assay-level support rather than independent transcriptomic validation.",
        ("ohira_2021_soleus",),
        "Mouse flight proteomics reported lower DECR1 in soleus under microgravity, matching the present gene, tissue, and direction across assays.",
        "Exact gene/protein, same-tissue, same-direction cross-assay alignment.",
    ),
    _annotation(
        "soleus",
        "Bnip3",
        "ambiguous",
        "mixed_direction_ground_analog_evidence",
        "Independent unloading evidence is sex dependent and does not provide one stable directional precedent.",
        ("rosa_caldwell_2021_unloading", "gambara_2017_soleus"),
        "Hindlimb unloading lowered soleus Bnip3 in male mice but raised it in female mice. The present flight-lower association agrees with one context and conflicts with another.",
        "Prior direction depends on sex and unloading context, so the literature relationship is ambiguous.",
    ),
    _annotation(
        "soleus",
        "Tpm1",
        "ambiguous",
        "mixed_gene_direction_and_process_evidence",
        "The present direction fits contractile remodeling but differs from prior exact-gene results in other muscle contexts.",
        ("allen_2009_muscle", "gambara_2017_soleus"),
        "Flight-higher Tpm1 is compatible with contractile and fiber-type remodeling in soleus, but prior flight-muscle reports include lower Tpm1 in other anatomical contexts.",
        "Process-level agreement and exact-gene directional disagreement make the precedent ambiguous.",
    ),
    _annotation(
        "tibialis_anterior",
        "Cdkn1a",
        "ambiguous",
        "mixed_exact_gene_direction_across_species_and_studies",
        "Public mouse synthesis and an independent rat mission provide different directional precedents.",
        ("oommen_2024_muscle", "taylor_2002_rat_muscle"),
        "A public mouse flight-muscle synthesis reported Cdkn1a among altered or higher common genes, while rat STS-90 tibialis-anterior and gastrocnemius profiling reported lower p21/Cdkn1a signaling.",
        "Direction varies across species, studies, and muscle contexts, so the literature relationship is ambiguous.",
    ),
    _annotation(
        "tibialis_anterior",
        "St3gal5",
        "complementary",
        "related_gene_mechanism",
        "Metabolic muscle mechanism only; no direct tibialis-anterior flight match was found.",
        ("yamashita_2003_st3gal5", "oommen_2024_muscle"),
        "ST3GAL5 produces GM3 ganglioside, which restrains insulin-receptor signaling in mouse skeletal muscle, while flight-muscle studies implicate glycosylation and metabolic adaptation. No prior directional tibialis-anterior flight result was found.",
        "Adds a ganglioside and insulin-signaling candidate to tibialis-anterior remodeling.",
    ),
    _annotation(
        "tibialis_anterior",
        "Bnip3",
        "complementary",
        "related_unloading_mitophagy_mechanism",
        "Independent unloading evidence supports context dependence but not a direct tibialis-anterior flight replication.",
        ("rosa_caldwell_2021_unloading", "oommen_2024_muscle"),
        "BNIP3 participates in mitochondrial quality control during unloading, but prior direction is sex and muscle dependent and no exact tibialis-anterior flight match was found.",
        "Adds a mitophagy candidate to tibialis-anterior adaptation without claiming directional replication.",
    ),
)


def build_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    synthetic_informed = pd.read_csv(SYNTHETIC_INFORMED_INPUT, sep="\t")
    unexpected_selection = set(synthetic_informed["selection_interpretation"]) - set(
        SELECTION_STATUS
    )
    if unexpected_selection:
        raise ValueError(
            f"Unexpected synthetic selection labels: {sorted(unexpected_selection)}"
        )
    synthetic_informed.insert(
        synthetic_informed.columns.get_loc("selection_interpretation") + 1,
        "selection_status",
        synthetic_informed["selection_interpretation"].map(SELECTION_STATUS),
    )
    synthetic_informed["_source_order"] = range(len(synthetic_informed))

    annotation_rows = []
    for annotation in ANNOTATIONS:
        row = asdict(annotation)
        row["source_ids"] = ";".join(annotation.source_ids)
        annotation_rows.append(row)
    annotations = pd.DataFrame(annotation_rows)
    annotations.insert(
        3,
        "interpretive_role",
        annotations["literature_classification"].map(
            INTERPRETIVE_ROLE_BY_CLASSIFICATION
        ),
    )
    if annotations["interpretive_role"].isna().any():
        raise ValueError("Every literature classification needs an interpretive role")

    source_rows = [asdict(source) for source in SOURCES]
    sources = pd.DataFrame(source_rows)

    key_columns = ["analysis_scope", "tissue", "symbol"]
    if synthetic_informed.duplicated(key_columns).any():
        raise ValueError(
            "Synthetic-informed input contains duplicate analysis-tissue-gene keys"
        )
    if annotations.duplicated(key_columns).any():
        raise ValueError("Literature annotations contain duplicate tissue-gene keys")
    if sources["source_id"].duplicated().any():
        raise ValueError("Literature source IDs must be unique")

    synthetic_informed_keys = set(
        map(
            tuple,
            synthetic_informed[key_columns].itertuples(index=False, name=None),
        )
    )
    annotation_keys = set(
        map(tuple, annotations[key_columns].itertuples(index=False, name=None))
    )
    if synthetic_informed_keys != annotation_keys:
        missing = sorted(synthetic_informed_keys - annotation_keys)
        extra = sorted(annotation_keys - synthetic_informed_keys)
        raise ValueError(f"Annotation coverage mismatch; missing={missing}, extra={extra}")

    observed_classes = set(annotations["literature_classification"])
    if not observed_classes <= ALLOWED_CLASSIFICATIONS:
        raise ValueError(
            f"Unexpected literature classifications: {observed_classes - ALLOWED_CLASSIFICATIONS}"
        )
    counts = (
        annotations["literature_classification"]
        .value_counts()
        .reindex(sorted(ALLOWED_CLASSIFICATIONS), fill_value=0)
        .to_dict()
    )
    if counts != {key: EXPECTED_COUNTS[key] for key in sorted(EXPECTED_COUNTS)}:
        raise ValueError(f"Unexpected classification counts: {counts}")

    referenced_sources = {
        source_id
        for value in annotations["source_ids"]
        for source_id in value.split(";")
    }
    missing_sources = referenced_sources - set(sources["source_id"])
    if missing_sources:
        raise ValueError(f"Unresolved literature source IDs: {sorted(missing_sources)}")

    merged = synthetic_informed.merge(
        annotations,
        on=key_columns,
        how="left",
        validate="one_to_one",
    )
    merged["literature_search_date"] = SEARCH_DATE
    merged = merged.sort_values("_source_order").drop(columns="_source_order")
    sources["literature_search_date"] = SEARCH_DATE
    return merged, sources


def write_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    annotations, sources = build_tables()
    annotations.to_csv(ANNOTATION_OUTPUT, sep="\t", index=False)
    sources.to_csv(SOURCE_OUTPUT, sep="\t", index=False)
    return annotations, sources


def check_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    expected_annotations, expected_sources = build_tables()
    observed_annotations = pd.read_csv(ANNOTATION_OUTPUT, sep="\t")
    observed_sources = pd.read_csv(SOURCE_OUTPUT, sep="\t")
    assert_frame_equal(observed_annotations, expected_annotations, check_dtype=False)
    assert_frame_equal(observed_sources, expected_sources, check_dtype=False)
    return observed_annotations, observed_sources


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify that the committed tables match the curated annotations.",
    )
    args = parser.parse_args()

    annotations, sources = check_tables() if args.check else write_tables()
    counts = annotations["literature_classification"].value_counts().to_dict()
    direct = annotations["evidence_scope"].eq(
        "direct_same_gene_same_tissue_same_direction"
    ).sum()
    print(
        f"Validated {len(annotations)} synthetic-informed associations against "
        f"{len(sources)} sources; direct matches={direct}; classes={counts}"
    )


if __name__ == "__main__":
    main()
