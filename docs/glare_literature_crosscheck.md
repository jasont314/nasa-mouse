# GLARE Literature Cross-Check

This note cross-checks the multi-tissue NASA mouse GLARE outputs against prior spaceflight and related rodent omics literature.

Local outputs reviewed:

- `outputs/glare_multi_tissue_api/validation_stack/candidate_modules.tsv`
- `outputs/glare_multi_tissue_api/validation_stack/candidate_module_score_meta.tsv`
- `outputs/glare_multi_tissue_api/validation_stack/candidate_module_panglao_enrichment.tsv`
- `outputs/glare_multi_tissue_api/validation_stack/metascape_gene_lists/candidate_module_gene_list_manifest.tsv`
- `outputs/glare_multi_tissue_api/*/dgea_comparison/GLARE_VS_DGEA_PER_STUDY.md`
- `outputs/glare_multi_tissue_api/*/dgea_comparison/recurring_dgea_glare_pathway_overlap.tsv`
- `outputs/glare_multi_tissue_api/*/dgea_comparison/significant_glare_reactome_terms_by_study.tsv`
- Cluster annotation TSVs under `outputs/glare_multi_tissue_api/*/.../plots/`

## Bottom Line

The most defensible GLARE results match known spaceflight biology: mitochondrial respiration, translation/protein metabolism, liver metabolic remodeling, and immune/lymphoid disruption. The strongest support is for modules that intersect DGEA and GLARE clustering.

Current validation summary:

| Module class | Significant modules | Median empirical p among all tested | Median empirical p among FDR-significant modules |
| --- | ---: | ---: | ---: |
| DGEA-intersection | 32 / 72 | 0.0050 | 0.0000 |
| GLARE-only | 17 / 72 | 0.1575 | 0.0500 |

Interpretation:

- DGEA-intersection modules are the main evidence.
- GLARE-only modules are useful hypotheses, but most need more validation.
- Liver olfactory/chemosensory terms should be retained as high-caution candidates, not automatically discarded.

## Tissue-Level Cross-Check

### Liver

Local GLARE/DGEA signal:

- Intersection modules: amino-acid metabolism up; immunoregulatory interactions, second-messenger signaling, peptide elongation, translational regulation, and SRP-dependent targeting down.
- GLARE-only modules: respiratory electron transport, platelet calcium response, Cyclin E/G1-S, Phase I functionalization, bile-acid/bile-salt synthesis, transcription termination.

Literature support: strong for liver metabolic remodeling, lipid handling, bile-acid/xenobiotic biology, oxidative stress, autophagy/proteostasis, and endocrine/metabolic stress.

Supporting sources:

- Beheshti et al. 2019 reported lipid dysregulation and fatty-acid processing in mouse liver spaceflight data. PMID:31844325, DOI: [10.1038/s41598-019-55869-2](https://doi.org/10.1038/s41598-019-55869-2)
- Jonscher et al. 2016 found spaceflight activation of lipotoxic pathways in mouse liver, including PPAR-related lipid and bile-acid changes. PMID:27097220, DOI: [10.1371/journal.pone.0152877](https://doi.org/10.1371/journal.pone.0152877)
- Blaber et al. 2017 reported liver autophagy/proteasome, aminoacyl-tRNA, and oxidative-defense shifts. PMID:28953266, DOI: [10.3390/ijms18102062](https://doi.org/10.3390/ijms18102062)
- Kurosawa et al. 2021 reported sulfur/glutathione metabolism changes. PMID:34750416, DOI: [10.1038/s41598-021-01129-1](https://doi.org/10.1038/s41598-021-01129-1)
- Uruno et al. 2021 reported Nrf2-linked lipid effects. PMID:34887485, DOI: [10.1038/s42003-021-02904-6](https://doi.org/10.1038/s42003-021-02904-6)

Assessment:

- Supported: amino-acid metabolism, lipid/bile/Phase I metabolism, translation/proteostasis, oxidative/metabolic stress.
- Partly supported: immune/second-messenger modules, because immune-metabolic coupling is well described but the module could reflect composition or study design.
- Ambiguous: liver olfactory/chemosensory receptor modules.

### Liver Olfactory/Chemosensory Terms

Local signal:

- Olfactory signaling appears repeatedly in liver clustering outputs.
- It also appears in many non-liver tissues, often as very large receptor-family clusters.
- Existing full validation outputs were generated before the code change that keeps liver olfactory terms as high-caution candidates, so rerunning validation is needed to score these modules directly.

Literature support: plausible but not yet validated in our data.

Supporting sources:

- OR10J5 has been linked to hepatocyte lipid handling. PMID:28842679, DOI: [10.1038/s41598-017-10379-x](https://doi.org/10.1038/s41598-017-10379-x)
- Olfr43 has been linked to hepatocyte lipid metabolism. PMID:30639733, DOI: [10.1016/j.bbalip.2019.01.004](https://doi.org/10.1016/j.bbalip.2019.01.004)

Assessment:

- Do not discard outright.
- Do not claim as a hidden GLARE discovery yet.
- Treat as high-caution until individual receptor genes show detectable expression, FLT-vs-GC shifts, acceptable mappability, liver specificity, and downstream signaling support.

### Skeletal Muscle

Local GLARE/DGEA signal:

- Aggregate skeletal muscle: respiratory electron transport, TCA/OXPHOS, protein metabolism, SRP/translation down; developmental and G1/S/cell-cycle modules up.
- Soleus: mitochondrial fatty-acid beta oxidation is one of the strongest validated modules.
- Gastrocnemius/quadriceps: TCA/OXPHOS and translation signals recur.
- EDL/tibialis anterior: cell-cycle, antigen/NF-kB/Wnt/apoptosis modules appear, but study counts are smaller.
- GLARE-only modules include respiratory/OXPHOS, Cyclin E, circadian clock, NGF/neuronal system, carbohydrate metabolism, and muscle contraction.

Literature support: strong for mitochondrial, contractile, calcium, oxidative, inflammatory, and subtype-specific muscle responses.

Supporting sources:

- Gambara et al. 2017 profiled slow-type calf soleus after 30-day spaceflight and supports soleus sensitivity, calcium/contractile, metabolic, inflammatory, and oxidative changes. PMID:28076365, DOI: [10.1371/journal.pone.0169314](https://doi.org/10.1371/journal.pone.0169314)
- Allen et al. 2009 analyzed murine skeletal muscle gene expression after spaceflight. DOI: [10.1152/japplphysiol.90780.2008](https://doi.org/10.1152/japplphysiol.90780.2008)
- Radugina et al. 2018 supports quadriceps atrophy and impaired regeneration after spaceflight. PMID:29475516, DOI: [10.1016/j.lssr.2017.08.005](https://doi.org/10.1016/j.lssr.2017.08.005)
- Vitry et al. 2022 supports liver-muscle crosstalk and muscle translation/autophagy/energy remodeling. PMID:36267920, DOI: [10.1016/j.isci.2022.105213](https://doi.org/10.1016/j.isci.2022.105213)

Assessment:

- Supported: mitochondrial/OXPHOS down, translation/protein metabolism shifts, muscle contraction/contractile remodeling, soleus-specific oxidative metabolism.
- Plausible but needs gene-level review: circadian, NGF/neuronal/NMJ-like modules.
- Ambiguous: amyloid and broad neuronal-system terms, especially when driven by generic Reactome membership rather than muscle-specific genes.

### Kidney

Local GLARE/DGEA signal:

- Intersection modules: translation/peptide elongation and mitochondrial respiration down.
- GLARE-only modules: membrane trafficking, mitochondrial fatty-acid beta oxidation, fatty-acid/TAG/ketone metabolism, Cyclin E, prolactin receptor signaling.

Literature support: moderate to strong for renal stress, lipid remodeling, ECM/tubular remodeling, and kidney risk in spaceflight.

Supporting sources:

- Siew et al. 2024 integrated pan-omics and reported renal transporter dephosphorylation, nephron remodeling, and renal dysfunction with simulated galactic cosmic radiation. PMID:38862484, DOI: [10.1038/s41467-024-49212-1](https://doi.org/10.1038/s41467-024-49212-1)
- Finch et al. 2025 reported mouse kidney transcriptomic lipid/ECM dysregulation, TGF-beta signaling, and strain-dependent effects. PMID:40133368, DOI: [10.1038/s41526-025-00465-0](https://doi.org/10.1038/s41526-025-00465-0)
- Suzuki et al. 2022 reported kidney effects after space travel including bone-mineralization, blood-pressure, lipid-metabolism, and Ugt1a-related changes. PMID:34767829, DOI: [10.1016/j.kint.2021.09.031](https://doi.org/10.1016/j.kint.2021.09.031)

Assessment:

- Supported: lipid/fatty-acid metabolism, renal stress/remodeling, membrane/tubular biology.
- Plausible: OXPHOS/translation down, but kidney literature more directly emphasizes lipid/ECM/tubule remodeling than translation alone.
- Weak/ambiguous: prolactin receptor signaling.

### Lung

Local GLARE/DGEA signal:

- Intersection modules: adaptive immune system, cytokine signaling, interferon signaling, innate immune system, GPVI activation, interleukin signaling are down.
- GLARE-only modules: respiratory electron transport, ERBB2/FGFR, platelet calcium response, circadian clock, smooth muscle contraction.

Literature support: partial. Lung-specific prior work emphasizes ECM, adhesion, profibrotic genes, remodeling, and injury; immune/interferon down is consistent with broader spaceflight immune dysregulation.

Supporting sources:

- Tian et al. 2010 reported ECM, adhesion, and profibrotic changes in mouse lung after spaceflight. PMID:19850731, DOI: [10.1152/japplphysiol.00730.2009](https://doi.org/10.1152/japplphysiol.00730.2009)
- Gridley et al. 2015 reported lung injury/remodeling after STS-135. PMID:26130787
- A 2025 lung transcriptome study reported few DEGs, ectopic olfactory/vomeronasal receptor upregulation, and protein-folding/circadian downregulation. PMID:40679517, DOI: [10.1007/s10517-025-06422-x](https://doi.org/10.1007/s10517-025-06422-x)

Assessment:

- Supported: immune/cytokine/interferon changes as part of systemic spaceflight immune disruption.
- Partly supported: smooth-muscle/contraction and circadian modules, but local score support is weak.
- Missing relative to prior lung literature: stronger ECM/profibrotic signal would be expected.

### Skin

Local GLARE/DGEA signal:

- Intersection modules: translation, peptide elongation, SRP-dependent targeting, immunoregulatory interaction, cell cycle, mitotic cell cycle.
- GLARE-only modules: respiratory electron transport, FGFR/FGFR disease, platelet calcium response, prolactin receptor signaling, circadian clock.

Literature support: moderate for oxidative stress, ECM/collagen/barrier, DNA repair, and mitochondrial dysregulation.

Supporting sources:

- Mao et al. 2014 reported ROS/antioxidant, ECM remodeling, and metabolic/cell-signaling changes in STS-135 mouse skin. PMID:24796731, DOI: [10.3109/10715762.2014.920086](https://doi.org/10.3109/10715762.2014.920086)
- Cope et al. 2024 reported spaceflight dermatology signals in skin barrier, collagen, DNA damage/repair, and mitochondria. PMID:38862781, DOI: [10.1038/s43856-024-00532-9](https://doi.org/10.1038/s43856-024-00532-9)

Assessment:

- Supported: cell-cycle/stress, immune, mitochondrial/oxidative, and barrier/ECM-related biology.
- Weak/ambiguous: prolactin receptor, FGFR, and circadian GLARE-only modules without stronger gene-level validation.

### Spleen And Thymus

Local GLARE/DGEA signal:

- Thymus: cell cycle, mitotic cell cycle, DNA replication, mRNA processing, capped pre-mRNA processing, mitotic M/M-G1 are strong and mostly down.
- Spleen: developmental biology, ECM organization, adaptive immune system, FGFR signaling, and TCR signaling.
- GLARE-only modules: thymus Golgi vesicle biogenesis, platelet calcium response, Cyclin E, transcription termination; spleen BCR signaling, xenobiotics, amyloids, biological oxidations.

Literature support: strong for immune and lymphoid disruption; thymus cell-cycle/DNA replication down is one of the clearest matches.

Supporting sources:

- Gridley et al. 2013 reported STS-135 thymus/spleen changes including thymic DNA fragmentation, spleen mass decrease, and altered thymus genes including Cdc25a, E2f1, and Myc down. PMID:24069384, DOI: [10.1371/journal.pone.0075097](https://doi.org/10.1371/journal.pone.0075097)
- Gridley et al. 2009 reported spaceflight effects on T lymphocyte distribution, function, and gene expression. PMID:18988762, DOI: [10.1152/japplphysiol.91126.2008](https://doi.org/10.1152/japplphysiol.91126.2008)
- Baqai et al. 2009 reported altered innate immune function and antioxidant gene expression after spaceflight. PMID:19342437, DOI: [10.1152/japplphysiol.91361.2008](https://doi.org/10.1152/japplphysiol.91361.2008)
- Pecaut et al. 2017 linked spaceflight immune dysfunction with systemic metabolic changes. PMID:28542236, DOI: [10.1371/journal.pone.0174174](https://doi.org/10.1371/journal.pone.0174174)

Assessment:

- Strongly supported: thymus cell-cycle/DNA replication down; immune/T-cell disruption.
- Plausible: spleen adaptive/TCR changes.
- Needs composition adjustment: spleen ECM/FGFR/development could reflect stromal/cell-population shifts rather than within-cell transcriptional regulation.

### Retina

Local GLARE status:

- Retina was skipped because the current TMS FACS pretraining source has no matching retina cells.
- No local retina GLARE module exists to validate.

Literature support: spaceflight retina biology is real, but not assessed by our GLARE run.

Supporting sources:

- Overbey et al. 2019 reported retinal DEGs, phototransduction/photoreceptor integrity, oxidative stress, and thinning after ISS flight. PMID:31527661, DOI: [10.1038/s41598-019-49453-x](https://doi.org/10.1038/s41598-019-49453-x)
- Chen et al. 2021 reported epigenetic/transcriptomic clock changes. PMID:34179686, DOI: [10.1093/pcmedi/pbab012](https://doi.org/10.1093/pcmedi/pbab012)
- Kremsky et al. 2024 reported artificial-gravity attenuation in optic nerve/retina. PMID:39596110, DOI: [10.3390/ijms252212041](https://doi.org/10.3390/ijms252212041)

Assessment:

- Retina should not be included in claims from this GLARE analysis.
- A retina-specific run needs a suitable retina single-cell pretraining source.

## GLARE-Only Versus DGEA-Intersection

DGEA-intersection modules are better supported by both validation statistics and prior literature.

Best-supported intersection families:

- Muscle and kidney OXPHOS/translation down.
- Soleus mitochondrial fatty-acid beta oxidation down.
- Liver amino-acid/lipid/bile/proteostasis/translation modules.
- Lung immune/interferon modules.
- Thymus DNA replication/cell-cycle down.
- Spleen adaptive/TCR modules.

GLARE-only modules with plausible follow-up value:

- Muscle and soleus OXPHOS/circadian/NMJ-like modules.
- Kidney fatty-acid/TAG/ketone metabolism and membrane trafficking.
- Liver Phase I and bile-acid modules.
- Thymus Golgi vesicle biogenesis and platelet calcium response.

GLARE-only modules that are currently weak or ambiguous:

- Olfactory/GPCR receptor-family clusters.
- Amyloids.
- Prolactin receptor signaling.
- Broad neuronal/NGF labels outside nervous tissue.
- Platelet cytosolic calcium response without matching platelet/hemostasis gene-level evidence.
- Xenobiotics/Phase I outside liver.
- Viral/influenza labels likely driven by host RNA/translation machinery.

## Criteria Before Claiming A Hidden GLARE Module

A GLARE-only hidden module should not be claimed as a discovery unless it passes most of these checks:

1. Consistent FLT-vs-GC module-score shift in held-out studies.
2. Stronger than matched random gene sets.
3. Direction consistency across studies or a clear design reason for study-specific direction.
4. Per-gene expression above detection, not only pathway membership.
5. Mappability/paralog checks for receptor families, ribosomal families, and other repetitive gene families.
6. Enrichment remains after removing broad OR/GPCR/ribosomal/viral superfamilies.
7. Cell-composition or tissue-composition adjustment.
8. Tissue specificity, not the same giant pathway cluster repeated across unrelated tissues.
9. Orthogonal validation by qPCR, protein, metabolite, histology, or public independent study.
10. For liver olfactory claims specifically: validate the actual receptor genes and downstream cAMP/PKA/CREB or hepatocyte-localized signaling, not just Reactome olfactory membership.

## Current Interpretation

The analysis appears biologically meaningful, but not primarily because of exotic hidden modules. The strongest finding is that GLARE organizes known spaceflight-responsive biology into coherent latent gene modules. The most publication-ready claims should focus on DGEA-supported GLARE modules, while GLARE-only modules should be framed as hypotheses for follow-up validation.
