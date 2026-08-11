# Distributed tissue-response hypotheses

## Purpose

This note connects synthetic-supported genes that do not form one enriched Reactome pathway to prior experimental literature. A distributed response is an organ-level hypothesis in which several cell types or functional compartments respond to the same stress in parallel. It does not imply that the selected genes interact directly.

The main matched all-gene analysis identified four flight-lower liver genes, one flight-higher skin gene, and one flight-higher spleen gene. The spleen consensus analysis added three genes, and the skin grouped analysis added RIPK1-regulated cell-death terms. The statistical and selection evidence remains in:

- `paper/synthetic_guided_spaceflight/source_data/table_s19_matched_all_gene_candidates.tsv`
- `paper/synthetic_guided_spaceflight/source_data/table_s22_matched_gene_literature_annotations.tsv`
- `paper/synthetic_guided_spaceflight/source_data/table_s23_grouped_pathway_literature_annotations.tsv`
- `paper/synthetic_guided_spaceflight/source_data/table_s16_promoted_gene_literature_annotations.tsv`

## Summary

| Tissue | Selected observations | Organ-level hypothesis | Main unresolved issue |
|---|---|---|---|
| Liver | Flight-lower *Grb10*, *Ppic*, *H2-DMa*, and *Gtf2a2* | Multicompartment adaptation involving metabolic feedback, stellate remodeling, antigen presentation, and transcriptional capacity | Regulation within cells versus altered liver cell composition |
| Spleen | Flight-higher *Loxl1*, *Rai14*, *Ptprk*, and *Myl9* | Remodeling of the stromal or vascular immune niche, with changes in matrix, mechanics, junctions, and immune-cell positioning | Whether the genes arise from one anatomical niche or several unrelated cell populations |
| Skin | Flight-higher *Plscr1* and RIPK1-regulated cell-death groups, with lower expiMap maintenance programs | Barrier stress accompanied by interferon signaling and altered control of damaged-cell fate | Regulated-death signaling versus active keratinocyte necroptosis |
| Thymus | Strong flight-lower mitotic program plus a smaller flight-higher immune-state set | Reduced proliferative renewal with a parallel stress, cytokine, or composition response | Loss of cycling thymocytes versus transcriptional repression within them |

## Liver: multicompartment adaptation to metabolic and cellular stress

### Evidence behind each gene

- The random-effects estimate for *Grb10* was lower in flight using 12 accessions. GRB10 regulates insulin, IGF1, and mTORC1 feedback. In adult mouse liver, ER stress can reactivate *Grb10*, and liver-specific loss of GRB10 reduced acute ER-stress lipogenesis and steatosis. Lower *Grb10* in the present analysis may therefore reflect altered growth-factor feedback or a compensatory response to metabolic stress. It should not be described simply as lower insulin signaling.
- The random-effects estimate for *Ppic* was lower in flight using 12 accessions. Experimental knockdown of *Ppic* reduced CCl4-induced fibrosis and TGF-beta-driven hepatic stellate-cell activation. The current direction raises the possibility of less PPIC-associated stellate activation, but it does not show reduced fibrosis in flight animals.
- The random-effects estimate for *H2-DMa* was lower in flight using 12 accessions. H2-DMalpha is needed for effective MHC class II peptide loading. The bulk signal is compatible with lower antigen-presentation activity or a lower abundance of antigen-presenting cells in liver.
- The random-effects estimate for *Gtf2a2* was lower in flight using 12 accessions. GTF2A2 is part of the TFIIA transcription-initiation complex. Prior mouse spaceflight reanalysis reported lower RNA-polymerase and protein-metabolism pathways, so this is process-level directional agreement rather than a published exact-gene replication.

Mouse spaceflight studies have independently reported hepatic lipid accumulation, altered insulin and growth-factor signaling, autophagy and proteasome responses, and changes in RNA-polymerase and protein-metabolism pathways. The four selected genes fit different parts of that background.

### Integrated hypothesis

The panel may describe parallel changes in several liver compartments:

1. Hepatocytes alter insulin, IGF, mTOR, and ER-stress feedback through *Grb10*.
2. Stellate cells reduce or reorganize PPIC-associated remodeling.
3. Resident or infiltrating antigen-presenting cells show lower MHC class II peptide-loading capacity, or become less abundant.
4. One or more compartments reduce basal transcriptional investment through *Gtf2a2*.

These changes could be part of an energy-allocation response to spaceflight-associated metabolic and cellular stress. They could also result from changes in cell fractions. The data do not support a direct four-gene pathway.

### Tests

- Use single-nucleus RNA-seq or spatial transcriptomics to assign each signal to hepatocytes, stellate cells, Kupffer cells, endothelial cells, or other immune cells.
- Measure GRB10 with insulin-AKT-mTOR activity and hepatic lipid burden.
- Measure PPIC with alpha-SMA, collagen, and stellate-cell abundance.
- Measure H2-DMalpha with MHC class II-positive antigen-presenting cells.
- Compare *Gtf2a2* with nascent transcription and RNA-polymerase pathway activity.

## Spleen: mechanical remodeling of an immune niche

### Evidence behind each gene

- *Loxl1* was higher in flight and passed the matched all-gene analysis. LOXL1 is required for elastic-fiber homeostasis. Loss of *Loxl1* changes splenic extracellular-matrix, immune, and cell-cycle expression.
- *Rai14* was higher in flight in the consensus analysis. RAI14 responds to F-actin and extracellular mechanical force and regulates Hippo-YAP signaling.
- *Ptprk* was higher in flight in the consensus analysis. PTPRK localizes to cell contacts and maintains junctional organization by dephosphorylating adhesion regulators.
- *Myl9* was higher in five of six spleen accessions. MYL9 contributes to actomyosin mechanics. Platelet-derived MYL9 structures can also recruit or retain CD69-positive inflammatory cells in vascular spaces.

Splenic fibroblastic stromal cells build white-pulp and red-pulp niches that support immune-cell localization and activation. This gives the four genes a plausible organ-level context even though they do not form one Reactome term.

### Integrated hypothesis

The panel may indicate remodeling of a splenic stromal or vascular niche. Higher LOXL1 could alter elastic-matrix organization, RAI14 could report altered force sensing, PTPRK could stabilize or reorganize junctions, and MYL9 could reflect contractile or platelet-associated changes that influence CD69-positive cell positioning. This provides a possible bridge to the lower T-cell receptor, C-type lectin, and neutrophil programs found by expiMap: altered tissue architecture could change where immune cells reside and how they receive signals.

This bridge is an inference. The genes may originate from stromal cells, endothelial cells, smooth-muscle-like cells, platelets, or immune cells in different splenic regions. Higher expression could also reflect altered cell abundance rather than remodeling within a fixed cell population.

### Tests

- Use spatial transcriptomics or multiplex imaging across white pulp, marginal zone, red pulp, and vasculature.
- Co-localize LOXL1 and elastic fibers with RAI14/YAP, PTPRK junctional proteins, MYL9, platelets, and CD69-positive immune cells.
- Quantify white-pulp area, stromal-network density, vascular structure, and immune-cell positioning.
- Repeat the expression analysis after estimating stromal, platelet, endothelial, and immune-cell fractions.

## Skin: barrier stress with interferon and cell-death control

### Evidence behind each result

- *Plscr1* was higher in flight in six skin accessions. PLSCR1 is interferon inducible and amplifies a subset of interferon-stimulated genes.
- Two overlapping RIPK1-regulated necrosis and necroptosis groups were higher in flight. Their selected members were *Cflar*, *Fas*, *Birc2*, and *Stub1*.
- expiMap found generally lower chromatin regulation, DNA repair, Hedgehog, sphingolipid, and junction programs.
- Prior multi-study spaceflight skin analysis reported changes in DNA repair, mitochondrial function, barrier genes, and collagen or extracellular matrix.

### Integrated hypothesis

Spaceflight skin may have reduced repair and barrier-maintenance capacity while increasing interferon-linked defense and cell-death checkpoint signaling. PLSCR1 could amplify the interferon arm. The RIPK1 groups could mark altered decisions between survival, apoptosis, and necroptosis in stressed keratinocytes. Because several members of these groups regulate or inhibit death as well as execute it, a higher pathway score is not evidence that necroptosis occurred.

### Tests

- Localize PLSCR1 and interferon-stimulated genes to epidermal and dermal compartments.
- Measure RIPK1, RIPK3, and MLKL activation together with keratinocyte death and inflammatory histology.
- Pair these measurements with collagen, junction, barrier, and DNA-damage assays.
- Stratify by dorsal versus femoral skin, mission, strain, diet, and recovery interval.

## Thymus

Thymus is not primarily a distributed-result case because its flight-lower genes form a coherent mitotic and DNA-replication program. A smaller flight-higher set, including *Plscr1*, *Socs2*, *Etv1*, and *Tspan3*, may represent a parallel interferon, cytokine, or immune-composition response. The most direct test is single-cell or spatial measurement of cycling thymocytes, dendritic cells, and stromal compartments.

## Exploratory secondary panels

The panels below came from the secondary consensus workflow and did not pass the final matched all-gene gate. They are useful for designing follow-up analyses, but they should not be presented at the same evidence level as liver, skin, spleen, or the central thymus result.

### Kidney

Flight-higher *Slc37a4* and *Inpp4b* suggest a renal metabolic-signaling response. SLC37A4 transports glucose-6-phosphate across the endoplasmic-reticulum membrane and is required for normal proximal-tubule glucose handling. INPP4B regulates PIP3-AKT signaling. Their shared direction could reflect coordinated adjustment of renal glucose handling and insulin or growth-factor signaling, consistent with prior spaceflight kidney reports of metabolic, kinase, and hormone-response changes. This is a two-gene hypothesis, and neither gene passed the matched all-gene gate.

Useful tests include proximal-tubule localization, glucose-6-phosphate transport, renal glycogen, AKT phosphorylation, and analysis of whether the effects persist after accounting for nephron-segment composition.

### Tibialis anterior

The flight-higher panel contained *Cdkn1a*, *St3gal5*, *Cebpd*, and *Bnip3*. These genes span stress-related cell-cycle arrest, GM3 ganglioside synthesis and insulin-receptor regulation, stress or inflammatory transcription, and mitochondrial quality control. A plausible distributed interpretation is that tibialis anterior enters a stress-adaptation state that combines growth arrest, altered insulin responsiveness, and mitophagy. The direction of *Cdkn1a* varies across flight and unloading studies, and the panel is consensus-only, so the combined mechanism remains exploratory.

Useful tests include fiber-type-resolved expression, p21 protein, insulin-AKT signaling, GM3 abundance, BNIP3-dependent mitophagy, and mitochondrial morphology.

### Gastrocnemius

Flight-higher *Nfkbia* and flight-lower *Fhl2* form a tentative two-axis response. Higher *Nfkbia* matches a prior gastrocnemius spaceflight result and may reflect negative feedback on NF-kappaB signaling. FHL2 participates in myogenic differentiation and autophagy, so lower *Fhl2* could accompany altered muscle maintenance. The pair suggests inflammatory feedback coupled to myogenic or autophagic remodeling, but it does not define one pathway and neither gene passed the matched all-gene gate.

### Adrenal gland

Flight-lower *Psmb8* and *Tspan4* do not yet support a coherent distributed hypothesis. PSMB8 provides an immune and proteostasis lead because it is an interferon-inducible immunoproteasome subunit. The targeted search found no sufficiently specific adrenal-flight mechanism for TSPAN4. Forcing these genes into one story would overstate the evidence.

Soleus is not analogous because its selected genes already form a coherent oxidative, mitochondrial, and contractile panel. Pooled skeletal muscle also produced a larger stress, metabolic, membrane, and adhesion panel with substantial prior alignment, but no individual gene passed the final all-gene importance gate. Eye had one selected gene, while retina, lung, EDL, and quadriceps did not yield an equivalent supported multi-gene panel.

## Primary literature

- Beheshti et al. 2019, mouse spaceflight liver lipid dysregulation: <https://doi.org/10.1038/s41598-019-55869-2>
- Blaber et al. 2017, mouse spaceflight liver proteostasis: <https://doi.org/10.3390/ijms18102062>
- Luo et al. 2018, GRB10 and ER-stress steatosis: <https://doi.org/10.1530/JME-18-0018>
- Yang et al. 2021, PPIC and hepatic stellate-cell activation: <https://doi.org/10.1016/j.toxlet.2021.06.021>
- Felix et al. 2000, H2-DMalpha and MHC class II peptide presentation: <https://doi.org/10.1084/jem.192.1.31>
- Vitry et al. 2022, liver and muscle metabolic response to spaceflight: <https://doi.org/10.1016/j.isci.2022.105213>
- Jeong et al. 2024, RAI14 mechanosensing and Hippo signaling: <https://doi.org/10.1038/s44319-024-00228-0>
- Fearnley et al. 2019, PTPRK and cell-cell adhesion: <https://doi.org/10.7554/eLife.44597>
- Hayashizaki et al. 2016, MYL9/12 and CD69-positive cell recruitment: <https://doi.org/10.1126/sciimmunol.aaf9154>
- Li et al. 2021, LOXL1-dependent matrix and spleen state: <https://doi.org/10.1111/cpr.13077>
- Alexandre et al. 2022, splenic fibroblastic stromal niches: <https://doi.org/10.1126/sciimmunol.abj0641>
- Cope et al. 2024, multi-study spaceflight skin analysis: <https://doi.org/10.1038/s43856-024-00532-9>
- Dong et al. 2004, PLSCR1 amplification of interferon signaling: <https://doi.org/10.1128/JVI.78.17.8983-8993.2004>
- Kumari et al. 2021, RIPK1-mediated keratinocyte necroptosis: <https://doi.org/10.26508/lsa.202000956>
- Hammond et al. 2018, mouse spaceflight liver and kidney pathway analysis: <https://doi.org/10.3390/ijms19124106>
- Kofuji et al. 2015, INPP4B regulation of PIP3-AKT signaling: <https://doi.org/10.1158/2159-8290.CD-14-1329>
- D'Acierno et al. 2022, SLC37A4 and proximal-tubule glucose handling: <https://doi.org/10.1681/ASN.2021070935>
- Allen et al. 2009, murine skeletal-muscle expression after spaceflight: <https://doi.org/10.1152/japplphysiol.90780.2008>
- Yamashita et al. 2003, GM3 and insulin sensitivity: <https://doi.org/10.1073/pnas.0635898100>
- Rosa-Caldwell et al. 2021, mitochondrial quality control during unloading: <https://doi.org/10.1002/jcsm.12809>
- Liu et al. 2019, FHL2 in myogenesis and autophagy: <https://doi.org/10.7150/ijbs.31371>
- Kitamura et al. 2011, PSMB8 in inflammation and proteostasis: <https://doi.org/10.1172/JCI58414>
