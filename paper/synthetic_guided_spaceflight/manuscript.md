<div class="title-page">

<h1>Cross-study synthetic-guided transcriptomics identifies thymic proliferative suppression, soleus metabolic remodeling, and splenic <em>Igfbp3</em> elevation in spaceflown mice</h1>

<p class="subtitle">Organ-specific flight responses across NASA OSDR mouse transcriptomes</p>

<p class="authors">Jason Trinh</p>

<p class="affiliation">Space Life Sciences Training Program, NASA Ames Research Center, Moffett Field, California, USA</p>

<p>Correspondence: jasontrinh@berkeley.edu</p>

<p class="draft-note"><strong>Research manuscript draft for author review.</strong> Author list, acknowledgments, repository release URL, and archival DOI require final review before submission.</p>

</div>

## Abstract

**Background:** Mouse spaceflight studies provide access to tissues that cannot be sampled extensively from astronauts, but individual experiments are small and differ in design. We asked whether synthetic gene-expression models could help identify reproducible, organ-specific flight responses without treating generated profiles as new animals.

**Methods:** We assembled 1,610 mouse flight and ground-control bulk RNA-seq profiles through the NASA Open Science Data Repository API. For the leakage-corrected held-out test, a conditional diffusion model was pretrained on 17,244 tissue-diverse ARCHS4 mouse profiles after excluding every GEO series linked to OSDR, and complete GEO series were assigned to separate reference roles. It was adapted only with designated OSDR development studies. The broader all-tissue feature screen preceded this overlap audit and is interpreted as hypothesis-generating. Generated profiles guided gene selection, while flight effects were estimated exclusively from real samples within each study and then compared across studies. Public spleen stromal references were used to investigate the likely cellular source of a focused spleen signal.

**Results:** The strongest result was a leakage-corrected accession-held-out thymus test. Flight samples showed lower expression of eight mitotic genes, including *Cdk1*, *Ccnb1*, *Ccnb2*, *Birc5*, and *Ube2c*, together with lower APC/C, G2/M-checkpoint, and DNA-replication programs. Synthetic-guided balanced accuracy increased from 0.500 to 0.833 and AUROC from 0.840 to 0.979. The pattern was shared by wild-type and Nrf2-knockout mice and is consistent with reduced thymic proliferative renewal. These genes were already prominent in real-only ranking, and synthetic guidance organized them into a more effective coherent panel. Anatomical separation of skeletal muscle exposed a soleus response involving lower oxidative fuel handling, mitochondrial quality control, and slow-muscle identity, with higher *Tpm1* suggesting contractile remodeling. Synthetic guidance reinforced five soleus genes found by real-only selection and additionally promoted *Mef2c* and *Pxmp2*. In spleen, synthetic guidance promoted *Igfbp3*, which was not stable in real-only selection; *Igfbp3* was subsequently found to be higher in flight in all six eligible real studies and remained significant when each study was removed in turn. Healthy mouse references localized baseline *Igfbp3* predominantly to fibroblastic reticular and perivascular stromal populations, suggesting a testable splenic-niche remodeling hypothesis. A broader real-data BH-FDR screen identified panel-level associations in 10 canonical tissues and all five anatomical muscle groups, while the synthetic-guided biological interpretation remained narrower.

**Conclusions:** Synthetic-guided analysis added the most information when it prioritized genes for testing in real samples rather than increasing the apparent sample size. The findings support reduced thymic proliferative renewal, identify a soleus metabolic and contractile response, and nominate flight-associated splenic *Igfbp3* as a focused stromal hypothesis.

**Keywords:** spaceflight; bulk RNA-seq; synthetic data; thymus; soleus; spleen; *Igfbp3*; NASA OSDR; ARCHS4

## Introduction

Spaceflight affects immune, musculoskeletal, metabolic, and barrier tissues through a combination of microgravity, radiation, confinement, altered nutrition, stress, and mission-specific procedures. Mouse flight experiments provide tissue access that is unavailable in astronauts, but their transcriptomic interpretation is difficult. Individual studies are small, missions differ in strain and duration, and condition labels can be entangled with study, material, genotype, or collection protocol. Pooling samples without preserving these design variables can convert study effects into apparent flight biology.

The NASA Open Science Data Repository (OSDR) now exposes sample metadata and processed assay data through a queryable biological API [1]. This makes it possible to assemble a cross-study cohort without relying on a precombined raw HDF5 object and to retain accession-level provenance. Public reference resources offer a second opportunity. ARCHS4 uniformly processes a large fraction of public human and mouse RNA-seq data [2], providing tissue-diverse reference profiles for pretraining models that would be underdetermined on OSDR alone.

Deep generative models can learn high-dimensional expression distributions. Conditional WGAN-GP models have reproduced tissue and cancer properties in GTEx and TCGA [3]. More recently, Lacan and colleagues adapted denoising diffusion probabilistic and implicit models to bulk transcriptomics and reported strong gene-correlation, neighborhood, adversarial, and downstream classification metrics [4]. GeneJEPA instead learns masked-gene representations without reconstructing expression [5]. These approaches solve different problems: a generator can sample expression, whereas a representation learner needs an additional decoder or generative objective before it can do so.

Synthetic expression is commonly presented as a remedy for small sample size. Generated profiles, however, are not new biological replicates. We instead used synthetic expression as a model-derived view of the data that could prioritize genes and pathways for testing in real flight samples.

Our primary biological questions were whether this approach could clarify the tissue response to spaceflight, whether anatomical separation would expose muscle-specific responses hidden by pooling, and whether the resulting signals complemented pathway-level findings from expiMap.

The clearest findings arose in thymus, soleus, and spleen. Thymus supplied a leakage-corrected test in a study excluded from ARCHS4 pretraining, OSDR adaptation, and feature-policy development. Because that outcome had been examined before the leakage correction, it is a retrospective sensitivity analysis rather than a prospective untouched confirmation. Soleus supplied a cross-study metabolic signal that is biologically coherent but still needs confirmation in a newly collected or fully unseen study. Spleen supplied a narrower but highly consistent gene-level association that motivated investigation of the splenic stromal compartment. Results from the remaining tissues define the exploratory boundary of the approach.

![Study design and evidence ladder.](figures/figure_1_study_design.png)

<p class="caption"><strong>Figure 1. Synthetic-guided biological analysis.</strong> (A) A mouse tissue reference model was trained with ARCHS4 and adapted to API-derived OSDR flight and ground-control profiles. Generated expression prioritized genes, while biological effects were measured in real samples. (B) The analysis focused on organ-specific responses in thymus, anatomically separated skeletal muscle, spleen, and the remaining tissue cohort.</p>

## Materials and methods

### Data sources

The OSDR Biological Data API was used to identify *Mus musculus* bulk RNA-seq assays with flight or ground-control labels [1]. Tissue and material names were harmonized while preserving study provenance. The resulting cohort contained 1,610 biological profiles from 75 accessions: 835 flight and 775 ground control. Full-transcriptome expression was converted to transcripts per million before selecting a 974-gene mouse landmark panel. No raw integrated OSDR H5 file was used.

The local ARCHS4 mouse resource contained 997,515 public RNA-seq profiles [2]. For the corrected holdout experiment, all nine GEO series linked from the API-derived OSDR metadata were excluded before reference selection, removing 108 otherwise eligible profiles and preventing cross-resource sample overlap. A healthy-preferred, tissue-balanced subset of 17,244 profiles spanning 20 tissue classes was then used for model pretraining. Complete GEO series were assigned to one reference role, producing 10,150 training, 2,466 validation, and 4,628 test profiles with no series overlap. The earlier broad development screen preceded this audit; its synthetic-selection labels are retained as exploratory annotations, whereas its real-data BH-FDR estimates do not depend on the generator.

**Table 1. Data scope.**

| Source | Profiles used | Biological scope | Role |
|---|---:|---|---|
| ARCHS4 mouse v2.5 | 17,244 | 20 tissue classes | Mouse tissue pretraining |
| NASA OSDR API | 1,610 | 75 accessions; 835 flight and 775 ground control | Spaceflight analysis |

### Conditional expression model

We reproduced the bulk-expression diffusion architecture of Lacan et al. [4], based on denoising diffusion and implicit sampling [13,14]. The ARCHS4-pretrained model was adapted to OSDR while representing tissue, flight status, study, and material type as separate conditions. This allowed the same model to generate flight or ground-control profiles for represented biological and study contexts. The leakage-corrected backbone and adaptation were used for the fixed lung/thymus test; the earlier broad adaptation supplied the all-tissue development screens.

Generated profiles were checked for preservation of tissue identity, expression structure, diversity, and flight-related differences on withheld real profiles. The evaluation included correlation, neighborhood overlap, adversarial separability, and distributional distance [15,16]. The exact architecture, transformations, split construction, calibration, thresholds, and comparator-model results are reported in the supplementary methods.

### Synthetic-guided biological analysis

We first tested generated profiles as additional training rows and then as a guide for gene selection. Direct augmentation was not beneficial. The retained workflow instead used synthetic expression to rank candidate features and fitted the final predictive models using real profiles.

To distinguish reinforcement from synthetic-specific prioritization, we compared repeated feature selection with and without the synthetic ranking. Genes selected by both approaches were interpreted as reinforced real-data signals. Genes selected only with synthetic guidance were treated as model-nominated candidates. Throughout this report, "synthetic-promoted" is shorthand for promoted across the repeated feature-selection threshold; it does not mean that a gene was absent from the real data or biologically novel. Generated profiles were never counted as biological replicates.

Flight-minus-ground effects were estimated within each OSDR study and then summarized with a random-effects model [17]. Within each tissue, the 974 gene-level meta-analysis P values were adjusted with the Benjamini-Hochberg procedure, and BH FDR below 0.05 defined the primary statistical inclusion rule [6]. Benjamini-Hochberg correction controls the expected proportion of false discoveries among the reported genes while retaining more power than family-wise error correction. Synthetic selection status was then recorded separately as reinforced, synthetic-promoted, real-only, or not stably selected. Accession-direction agreement, between-study heterogeneity, and leave-one-accession-out results were retained as interpretation and sensitivity measures rather than inclusion requirements. The complete BH-FDR inventory is provided in Supplementary Table S17, its synthetic-guided subset in Table S16, and tissue-level counts in Table S18. Reactome was used to group selected genes into biological processes [7].

### Evidence hierarchy

The analyses answer three different questions and were interpreted in that order. The highest evidence level required an OSDR study to be excluded from ARCHS4 pretraining, OSDR adaptation, and feature-policy development. The thymus rerun satisfies those computational exclusions, but it is described as leakage-corrected rather than prospectively independent because its outcome was known before retraining. Direction agreement and LOO stability refine confidence within a tier but do not move a result between tiers.

**Table 2. Analysis and evidence levels.**

| Tier | Analysis | Question answered | Gene-level output | Permitted interpretation |
|---|---|---|---|---|
| 1 | Leakage-corrected accession holdout | Does the synthetic-guided policy survive complete reference-overlap removal and accession holdout? | Eight directional thymus core genes in OSD-457 | Strongest retrospective evidence of transferable feature selection; not a prospective untouched confirmation |
| 2 | Repeated development screen plus real-data BH FDR | Which BH-significant genes were repeatedly promoted or reinforced by synthetic guidance? | 52 tissue-gene results: 28 promoted and 24 reinforced | Cross-study developmental hypotheses; not independent synthetic generalization |
| 3 | Complete real-data random-effects screen | Which panel genes show a pooled FLT-GC association? | 459 tissue-gene results across canonical and muscle-group analyses | Conventional panel-level association; does not show that synthetic data added information |

Signals observed across represented studies but lacking a Tier 1 test are described as developmental or exploratory. Generative findings were also compared with the separate expiMap pathway analysis [12] to identify convergent and complementary tissue responses.

For the spleen follow-up, study-level *Igfbp3* effects were checked against full-transcriptome TPM values. Healthy sorted-cell and single-cell references, GSE156162 and E-MTAB-7703, were then examined to identify the normal spleen populations expressing *Igfbp3* [18,19]. These references were used for source localization after discovery and were not part of model training or gene selection.

## Results

### Synthetic expression preserved tissue structure and prioritized biological features

The leakage-corrected ARCHS4 model generated profiles that retained broad mouse tissue identity: classifiers trained on real and generated reference profiles both achieved tissue balanced accuracy 0.781 on the complete-series test split. Adversarial accuracy was 0.515, precision was 0.951, and recall was 0.890. However, gene-correlation-matrix agreement was 0.878 and did not pass the prespecified 0.952 gate. The reference is therefore adequate for the fixed feature-guidance sensitivity test, but not validated as a general substitute for real expression. Complete distributional results, denoising trajectories, and comparator-model screens are reported in Supplementary Figures S1-S4 and Supplementary Tables S1-S3.

Simply adding generated rows to the training data did not improve flight-versus-ground classification. The more useful strategy was to let the synthetic profiles influence which genes were considered, while fitting the final classifier and estimating biological effects only from real samples. This strategy transferred clearly to thymus and produced mixed results in lung. Detailed predictive comparisons are provided in Supplementary Figure S5 and Supplementary Table S4.

The Tier 3 real-data screen yielded 202 BH-FDR tissue-gene associations across 10 of 22 canonical tissue analyses and 257 across the five anatomical muscle groups. These are tissue-gene results rather than 459 unique or independent discoveries: a gene can be significant in more than one tissue, and pooled skeletal muscle overlaps its anatomical subgroups. Fifty-two results also entered Tier 2 stable synthetic-informed selection: 28 were synthetic-promoted and 24 were reinforced by both real-only and synthetic-guided selection. The complete set, including associations with mixed study directions, remains available in Supplementary Tables S16-S18. The biological sections below emphasize the Tier 1 result and smaller Tier 2 process-level patterns rather than treating every FDR result as a mechanistic discovery.

### Thymus shows lower proliferative renewal during spaceflight

OSD-457 provided the strongest result after a leakage audit found that its GEO series had been present in the original ARCHS4 reference. The reference was retrained from scratch after excluding every OSDR-linked GEO series and assigning complete GEO series to separate roles. OSD-457 was also excluded from OSDR adaptation and feature-policy development. In this corrected test, synthetic-guided gene selection improved separation of flight and ground-control thymus samples from balanced accuracy 0.500 to 0.833 and AUROC 0.840 to 0.979. The result held in both wild-type and Nrf2-knockout mice. Flight effects were closely aligned between the two genotypes.

The core genes *Birc5*, *Ccne2*, *Gmnn*, *Ube2c*, *Cdk1*, *Nusap1*, *Ccnb1*, and *Ccnb2* were lower in flight in both strata (Fig. 2A). These were not genes absent from the real-data analysis: all eight were among its 26 highest-ranked candidates, including four among the top ten. Synthetic guidance instead assembled a broader, coherent cell-cycle panel that transferred more effectively to the study-excluded test. Together, these genes regulate DNA replication, chromosome progression, mitotic entry, and completion of cell division. Reactome analysis connected them to APC/C-mediated cyclin degradation, G2/M checkpoints, DNA synthesis, and broader cell-cycle control (Fig. 2B).

The leakage-corrected Tier 1 result is therefore interpreted as panel-level reinforcement and organization, not de novo gene discovery. In the separate Tier 2 repeated screen, *Ccne2*, *Cdk1*, *Nusap1*, *Ccnb1*, and *Ccnb2* crossed the stable-selection threshold only with synthetic guidance; *Gmnn* was reinforced by both arms, *Ube2c* was stable in the real-only arm, and *Birc5* was BH-significant but did not cross either stability threshold. These labels describe repeated selection behavior and do not override the stronger held-out result. The complete cross-analysis mapping is provided in Supplementary Table S19.

This result is aligned with, but more specific than, prior thymus studies. STS-135 mouse thymus showed changes in cell-cycle and DNA-damage programs, including lower checkpoint-related expression [8]. A later ISS experiment reported marked thymus mass loss and partial artificial-gravity rescue of cell-cycle expression [9]. The current signature emphasizes mitotic completion and replication rather than acute apoptosis alone. Agreement between genotype strata suggests that the predictive signature is not confined to one Nrf2 background, but it does not prove Nrf2 independence.

Bulk thymus expression cannot distinguish lower transcription within proliferating thymocytes from loss or redistribution of proliferating cell populations. The defensible biological conclusion is lower abundance of a mitotic transcript program in flight, consistent with reduced proliferative renewal. Cell-resolved or histological confirmation is required to assign the effect to a cell-intrinsic mechanism.

![Thymus biology.](figures/figure_4_thymus_biology.png)

<p class="caption"><strong>Figure 2. Thymus response to spaceflight.</strong> (A) Flight-minus-ground effects in real OSD-457 profiles for eight cell-cycle genes that were lower in both wild-type and Nrf2-knockout mice. (B) The genes converge on a mitotic and DNA-replication process family.</p>

### Anatomical separation exposes a soleus-specific metabolic program

Aggregate skeletal muscle concealed substantial anatomical heterogeneity. We therefore examined extensor digitorum longus, gastrocnemius, quadriceps, soleus, and tibialis anterior separately. Soleus produced the clearest biological pattern: its selected genes showed consistent flight effects across three accessions and converged on related metabolic processes.

Seven synthetic-selected genes remained directionally consistent when each accession was examined in turn and agreed with the generated flight effect. *Bdh1*, *Bnip3*, *Mef2c*, *Ech1*, *Pxmp2*, and *Gmnn* were lower in flight, while *Tpm1* was higher (Fig. 3B). These genes were associated with mitochondrial protein turnover, fatty-acid oxidation, and lipid metabolism (Fig. 3C).

Five of these genes, *Gmnn*, *Tpm1*, *Bdh1*, *Ech1*, and *Bnip3*, were stable in both real-only and synthetic-guided selection. Synthetic guidance specifically added *Mef2c* and *Pxmp2* to the retained set; their flight effects were then supported using the real profiles. Soleus therefore illustrates both reinforcement of an existing signal and nomination of additional genes within the same biological program.

Using BH FDR without requiring leave-one-accession-out significance added reinforced *Decr1*, which was also lower in flight in all three soleus accessions. Its weaker study-removal margin makes it an extension of the metabolic set rather than one of the seven central genes.

The pattern links ketone or lipid utilization (*Bdh1*, *Ech1*), mitochondrial quality control (*Bnip3*), slow oxidative muscle identity (*Mef2c*), peroxisomal transport (*Pxmp2*), and contractile remodeling (*Tpm1*). Prior 30-day spaceflight profiling of mouse soleus reported a slow-to-fast shift and broad changes in oxidative metabolism, PPAR signaling, and contractile genes [10]. Unloading studies have also reported reduced soleus fatty-acid oxidation [11]. The pathway is therefore literature-aligned, while the compact gene prioritization and peroxisome-mitochondria emphasis are exploratory refinements rather than wholly de novo biology.

Unlike thymus, soleus was represented during model development. Its cross-study consistency makes it a focused biological hypothesis, but an entirely unseen soleus study is still needed for independent confirmation.

![Soleus biology.](figures/figure_5_soleus_biology.png)

<p class="caption"><strong>Figure 3. Skeletal-muscle and soleus response.</strong> (A) Number of synthetic-prioritized genes with consistent real effects across studies in each anatomical muscle group. (B) Seven soleus genes with consistent real flight effects. (C) Their strongest shared biological processes center on mitochondrial turnover and lipid metabolism.</p>

### Other muscle groups provide narrower hypotheses

The other muscle groups produced narrower hypotheses. Quadriceps contained four synthetic-informed BH-FDR genes that were higher in flight in all four accessions. Synthetic guidance promoted *Cebpd*, *Rbm6*, and *Sh3bp5*, while *Gpatch8* was reinforced by both selection approaches. *Rbm6* was retained in six of eight synthetic-guided analyses and was the only one of the four to remain significant after every individual study was removed. The gene functions suggest two possible components: *Rbm6* and *Gpatch8* implicate RNA processing, with RBM6 also linked experimentally to DNA double-strand-break repair, while SH3BP5 is a mitochondrial JNK-interacting protein [28-30]. Together with stress-responsive *Cebpd*, these genes nominate an RNA-processing and mitochondrial-stress response in quadriceps. G1/S and TP53-related pathways were suggestive but did not meet the corrected significance threshold, and *Rbm6* itself was not represented in the mouse Reactome pathway file. The proposed connection is therefore a functional hypothesis rather than a significant shared pathway.

EDL showed lower *Abcc5*, *Lsm6*, *Polr2i*, and *Tsc22d3* in flight across its two accessions, suggesting RNA-processing and nuclear-receptor responses. Tibialis anterior showed higher *Cdkn1a*, *St3gal5*, *Cebpd*, *Pdhx*, and *Bnip3* across both accessions, consistent with stress and metabolic remodeling, but synthetic guidance added little to an already strong real-data separation. Gastrocnemius produced one synthetic-promoted gene, *Cxcr4*, that was higher in flight in all three accessions but lacked a significant shared pathway.

### Flight-associated spleen *Igfbp3* points to stromal-niche remodeling

Spleen produced the strongest single-gene result outside the thymus and soleus programs. Synthetic-guided analysis repeatedly selected *Igfbp3*, whereas it was not a stable feature in the real-only selection arm. The biological association itself was then evaluated in real expression data. *Igfbp3* was higher in flight in all six eligible spleen studies, with a random-effects FDR of `1.76e-9`; the largest FDR after removing one study at a time was `0.00385`. A separate full-transcriptome TPM calculation showed 9% to 63% higher group means in flight across the same six studies.

The broader BH-FDR spleen set also contained synthetic-promoted *Rai14* and *Ptprk* and reinforced *Bace2* and *Loxl1*. Each was higher in flight in all six studies, but none remained as statistically stable as *Igfbp3* when studies were removed. Although these genes did not enrich one corrected pathway, three provide a tentative architectural extension of the *Igfbp3* hypothesis: RAI14 links actin-associated mechanosensing to Hippo signaling, PTPRK regulates cell-cell junctions, and LOXL1 supports elastic extracellular-matrix maintenance [25-27]. Their shared direction is compatible with altered stromal mechanics, adhesion, or matrix organization, but those functions were established in other biological systems and do not prove that the genes change in the same splenic cell population. *Bace2* remains outside this proposed mechanism. Synthetic-promoted *Snca* also passed BH FDR, but its pooled flight-lower effect agreed in only four of six studies. It remains in the primary FDR inventory but is not used to extend the stromal mechanism.

The result did not extend to a significant IGF/IGFBP Reactome program, indicating a focused gene association rather than a coordinated pathway shift. Healthy spleen references nevertheless provided a specific source hypothesis. In sorted populations from GSE156162, *Igfbp3* expression was highest in white-pulp mesenchymal cells and lower in red-pulp mesenchymal cells, with little expression in endothelial cells or red-pulp macrophages [18]. E-MTAB-7703 independently localized expression to fibroblastic reticular, collagen-producing, and perivascular stromal populations [19]. Existing flight splenocyte, blood, and marrow single-cell preparations contained too little *Igfbp3* to test this source, but those preparations also exclude or strongly deplete the relevant nonhematopoietic stroma.

Spaceflight proteomics has reported increased IGFBP3 in serum and femur under microgravity relative to onboard artificial gravity [20]. This provides systemic context but does not establish the spleen as the source. The current data support a testable hypothesis that spaceflight alters *Igfbp3* abundance within the splenic stromal niche through increased expression per stromal cell, altered stromal abundance or architecture, or both.

### Other organ responses were heterogeneous

Lung highlighted cell-cycle, senescence, and PI3K/AKT-related genes, but the pattern varied by genotype and study and no gene passed the primary BH-FDR screen. Skin selected G1/S and DNA-repair genes, including *Ccne2*, *Ccnd1*, *Chek1*, *Brca1*, and *Topbp1*. This reproduces the cell-cycle and DNA-damage themes reported in an independent analysis of the same class of OSDR skin studies, which also found substantial variation by mission, strain, recovery interval, and anatomical site [21]. In the real-data panel screen, flight was associated with higher *Nfkbie* and *Plscr1* and lower *Arid4b* at BH FDR below 0.05; only *Plscr1* had the same direction in all six studies, and none of the three was stably selected through synthetic guidance. The skin generated-only classifier did not improve the reserved profiles, so this is literature-aligned support rather than a new synthetic-guided discovery.

Kidney produced a more focused secondary lead. *Slc37a4* was selected in every repeated real and generated-guided analysis and was higher in flight in all six kidney accessions. It passed BH FDR in both the prespecified 974-gene panel and a full-transcriptome TPM sensitivity analysis, although its significance weakened when individual studies were removed. *Slc37a4* encodes the glucose-6-phosphate translocase used in renal glucose production and supports a metabolic-adaptation hypothesis. *Inpp4b* was an additional strong real-data signal but was not selected through synthetic guidance. The earlier *Hmox1*/*Alas1* porphyrin enrichment was not supported by significant real effects for either gene and is therefore retained only as an oxidative-stress hypothesis. The broader kidney pattern agrees with independent reports of strain-dependent lipid, ECM, TGF-beta, and Wnt responses and nephron remodeling after spaceflight [22-24].

Retina independently produced synthetic-promoted *Slc37a4*, which was higher in flight in all four accessions and passed BH FDR. Its enriched shear-stress terms were driven by different selected genes, so the gene and pathway observations do not yet form one mechanism. In additional API-derived tissues, adrenal gland showed flight-lower synthetic-promoted *Psmb8*, *Ticam1*, and *Pmaip1* together with reinforced *Tspan4* across all three accessions; eye showed reinforced flight-lower *Klhl21*. These are included as exploratory findings because they lack an unseen study and a sufficiently developed organ-level interpretation. Liver had 19 real-data BH-FDR genes, eight higher and 11 lower in flight, but none entered the stable synthetic-informed set and none had the same direction in all 12 studies. Liver is therefore negative for a coherent synthetic-guided result, not for real panel-level FLT-GC association.

The evidence distribution therefore differs from the separate expiMap analysis. Both approaches prioritize thymus. expiMap produced its broadest pathway evidence in thymus, skin, spleen, and kidney, whereas the generative workflow added its clearest complementary result in soleus. The methods examine different biological representations: expiMap tests predefined pathway activity, while synthetic guidance helps prioritize individual genes and their shared processes.

![Spleen Igfbp3 and tissue evidence hierarchy.](figures/figure_6_tissue_evidence.png)

<p class="caption"><strong>Figure 4. Spleen <em>Igfbp3</em> and the tissue evidence hierarchy.</strong> (A) Real flight-minus-ground <em>Igfbp3</em> effects were positive in all six spleen studies; error bars show 95% intervals and the diamond shows the random-effects estimate. Parentheses give flight/ground-control sample counts. (B) Healthy sorted-cell data localize baseline expression primarily to white- and red-pulp mesenchymal populations. (C) Thymus remains the strongest result, while soleus, quadriceps, and spleen provide complementary cross-study findings. Kidney and retina are secondary candidates; lung and skin are heterogeneous; liver lacks a coherent synthetic-guided pattern.</p>

**Table 3. Selected synthetic-guided biological interpretations by tissue. Complete real-data BH-FDR results are in Supplementary Tables S17-S18.**

| Tissue | Main signal | Interpretation |
|---|---|---|
| Thymus | Lower mitotic genes; APC/C, G2/M, and DNA replication | Strongest result; supports reduced proliferative renewal |
| Soleus | FLT-lower *Bdh1*, *Ech1*, *Bnip3*, *Gmnn*, *Mef2c*, *Pxmp2*, and *Decr1*; FLT-higher *Tpm1* | Coherent metabolic program; *Mef2c* and *Pxmp2* were synthetic-promoted |
| Quadriceps | FLT-higher *Cebpd*, *Rbm6*, *Sh3bp5*, and *Gpatch8* across four studies | RNA-processing and mitochondrial-stress hypothesis; no significant shared pathway |
| Spleen | FLT-higher *Igfbp3*, *Rai14*, *Ptprk*, *Bace2*, and *Loxl1* across six studies | *Igfbp3*-centered stromal-niche hypothesis with exploratory adhesion and ECM support |
| EDL | FLT-lower *Abcc5*, *Lsm6*, *Polr2i*, and *Tsc22d3* across two studies | Exploratory RNA-regulation response |
| Tibialis anterior | FLT-higher *Cdkn1a*, *St3gal5*, *Cebpd*, *Pdhx*, and *Bnip3* across two studies | Exploratory stress and metabolic response |
| Gastrocnemius | FLT-higher *Cxcr4* across three studies | Isolated synthetic-promoted gene |
| Kidney | Flight-higher *Slc37a4* and renal metabolic remodeling | Credible secondary hypothesis; requires an unseen kidney study |
| Retina | Flight-higher *Slc37a4* across four studies | Exploratory; enriched pathways involve different genes |
| Adrenal gland | FLT-lower *Psmb8*, *Ticam1*, *Pmaip1*, and *Tspan4* across three studies | Exploratory stress and translation-regulation response |
| Lung | Cell-cycle, senescence, and PI3K/AKT candidates | Mixed across genotype and study |
| Skin | Cell-cycle and DNA-repair candidates | Literature-aligned but heterogeneous; no new robust synthetic-guided gene |
| Liver | No coherent retained pattern | No synthetic-guided biological conclusion |

## Discussion

### What synthetic data added

The generated profiles did not create stronger evidence simply by increasing the number of training rows. Their useful contribution was either to reinforce coherent combinations already present in the real data or to promote candidates that were not stable under real-only selection. The synthetic model therefore acted as a feature-prior: it influenced what to examine, while the biological conclusions remained based on real flight and ground-control samples.

Thymus demonstrates panel-level reinforcement: its central genes already ranked strongly in real data, but the synthetic-guided panel transferred after the study was removed from both reference pretraining and OSDR adaptation and converged on one interpretable cell-cycle process. The Tier 2 gene labels describe which genes crossed repeated selection thresholds and are not claims of biological novelty. Soleus demonstrates a mixed contribution: synthetic guidance reinforced five real-selected genes and promoted *Mef2c* and *Pxmp2*, localizing a broad muscle response to oxidative soleus biology. Neither result depends on counting synthetic profiles as animals.

Spleen demonstrates a third use. Synthetic guidance elevated *Igfbp3*, which was not stable in real-only feature selection, and the resulting hypothesis was supported by a consistent association across six real studies. The broader set makes that hypothesis more biologically connected: synthetic-promoted *Rai14* and *Ptprk*, together with reinforced *Loxl1*, span mechanosensing, cell contact, and extracellular-matrix organization. This convergence complements the stromal localization of *Igfbp3*, but it remains a hypothesis assembled from gene functions rather than a significant pathway or cell-resolved mechanism. The synthetic profiles nominated the genes; they did not supply the biological replication.

Quadriceps provides a parallel but narrower example. Synthetic guidance elevated *Cebpd*, *Rbm6*, and *Sh3bp5*, and all three were higher in flight across four real studies; reinforced *Gpatch8* followed the same direction. The pairing of *Rbm6* with *Gpatch8* supplies an RNA-processing component, while *Sh3bp5* supplies a mitochondrial stress-signaling component and *Cebpd* is consistent with stress-responsive transcription. This yields a testable hypothesis of altered RNA handling or genome maintenance coupled to mitochondrial stress, but not a single demonstrated pathway. It is therefore retained as a focused secondary lead rather than interpreted at the same biological depth as the spleen stromal hypothesis.

### Thymus, soleus, and spleen define distinct organ responses

The thymus result refines established spaceflight immune biology. Prior work documents thymic involution and altered cell-cycle expression [8,9]. The current leakage-corrected signature concentrates on cyclins, CDK1, UBE2C, BIRC5, NUSAP1, geminin, APC/C-mediated protein turnover, and G2/M control. Together these support lower proliferative renewal or a lower proportion of cycling thymocytes. Because the data are bulk, composition and cell-intrinsic regulation remain inseparable.

The soleus result addresses a different physiological axis. Weight-bearing slow muscle is especially sensitive to unloading. The selected genes describe lower oxidative substrate handling, altered mitochondrial turnover, reduced slow-muscle transcriptional identity, and contractile remodeling. This is compatible with known soleus atrophy and slow-to-fast transition [10,11], while nominating a compact set of genes for targeted validation. It is complementary to the expiMap result because synthetic guidance localized a broad fatty-acid-oxidation response to a smaller set observed consistently across soleus studies.

The spleen result adds an immune-organ remodeling hypothesis distinct from the thymic proliferative program. Fibroblastic reticular and perivascular stromal cells organize lymphocyte niches, vascular interfaces, and white-pulp architecture [18,19]. Higher whole-spleen *Igfbp3* could therefore reflect altered signaling within those cells or a change in their abundance and organization. IGFBP3 can regulate IGF availability and also exert IGF-independent effects, but the present data do not establish which function, if any, changes in spleen. A stroma-preserving assay is required before connecting the transcript association to immune function.

### Why the other tissues remain useful

The other tissues constrain the method's scope. Lung shows that a model can separate conditions without yielding one reproducible pathway or BH-FDR gene. Skin shows that synthetic feature selection can recover a known biological theme without producing a transferable classifier or a new stable synthetic-guided gene claim. Kidney shows a useful but less stable complementarity: synthetic guidance reinforces *Slc37a4*, while expiMap and prior work place that candidate within broader renal metabolic and remodeling responses. Retina shows that a generated-promoted gene and an enriched pathway can coexist without defining the same mechanism. Liver shows that a tissue can contain BH-FDR real-data associations while a broadly realistic generator still fails to define a coherent synthetic-guided program.

These results still provide actionable hypotheses. Lung cell-cycle and PI3K/AKT candidates can be tested in a newly held-out study with genotype modeled prospectively. Skin requires mission- and anatomical-site-aware analysis rather than one pooled effect. Kidney *Slc37a4* should be tested prospectively together with renal glucose handling, lipid metabolism, and ECM remodeling in a new study. Quadriceps *Rbm6* should first receive full-transcriptome confirmation and targeted RNA-processing follow-up in a new study. For spleen, the direct next experiment is measurement of *Igfbp3* RNA and IGFBP3 protein in stroma-preserving sections or sorted CD45-negative fractions, together with fibroblastic reticular and perivascular markers.

### Limitations

No tissue currently has a completely prospective untouched biological test. OSD-457 was excluded from the leakage-corrected reference, adaptation, and feature-policy development, but its outcome had already been examined before the correction. It therefore provides strong retrospective sensitivity evidence rather than prospective confirmation. All other BH-FDR findings were identified within the broader model-development domain and still lack equivalent study-excluded testing. Agreement across represented real studies supports focused follow-up but does not establish transfer to a new mission or study.

The analysis used a 974-gene landmark panel, so relevant genes outside that panel could not be discovered. BH FDR was controlled separately within each declared tissue family, not once across every tissue-gene combination in the project. Direction disagreement and high heterogeneity do not remove a BH-significant association, but they reduce confidence that its pooled effect represents a common response across missions. Bulk tissue also cannot distinguish a transcriptional change within a cell type from a change in cell composition. The thymus result could reflect lower expression in proliferating thymocytes, fewer proliferating thymocytes, or both. Likewise, the spleen result cannot distinguish stromal induction of *Igfbp3* from altered stromal abundance or tissue architecture.

Finally, the generator represents tissues and study contexts available during training. The corrected model failed the strict gene-correlation-structure gate, and the broader all-tissue synthetic-selection screen was not rerun after the cross-resource overlap audit. Those tissue-level synthetic labels are therefore hypothesis-generating; the associated real-data BH-FDR effects remain valid as real-data analyses. The model should not be assumed to reproduce a new mission, strain, or sample-processing protocol without additional testing. Exact model limitations, sensitivity analyses, and statistical safeguards are documented in the supplementary methods.

## Conclusions

Synthetic expression was most informative as a guide to biological feature selection, not as a replacement for real animals. The leakage-corrected analysis supports a flight-lower thymus mitotic program, identifies a soleus response centered on oxidative metabolism and contractile remodeling, and nominates flight-higher splenic *Igfbp3* as a focused stromal-niche hypothesis.

Thymus remains the strongest result because it survived complete OSDR-series removal from reference pretraining and accession exclusion from adaptation. Soleus supplies the clearest complementary metabolic program, while spleen supplies the most developed secondary gene association. Quadriceps provides a narrower four-gene synthetic-guided set, and kidney, retina, adrenal gland, EDL, gastrocnemius, and tibialis anterior provide exploratory BH-FDR hypotheses. Skin remains heterogeneous, lung did not retain a synthetic-guided classifier after correction, and liver sets the clearest negative boundary for synthetic-guided interpretation despite containing real-data BH-FDR associations.

## Data and code availability

OSDR data were accessed through the public Biological Data API [1]. ARCHS4 and Reactome are public resources [2,7]. Repository scripts are under `src/nasa_mouse_rna_diffusion/`; frozen run outputs are under `outputs/generative_benchmark/`; this paper package is under `paper/synthetic_guided_spaceflight/`. The figure and table builder is `nasa_mouse_rna_diffusion.build_synthetic_guided_paper`. SHA-256 hashes of every frozen analysis input are included in `source_data/frozen_input_manifest.tsv`. Public repository URL and archival DOI will be added after author review.

## Ethics statement

This study is a secondary computational analysis of publicly available animal-study data. No new animals were used.

## Competing interests

The author declares no competing interests.

## Acknowledgments

The author acknowledges the NASA Open Science Data Repository, GeneLab data-processing teams, original flight-study investigators, ARCHS4, Reactome, and the developers of the evaluated generative-model implementations. Program and mentor acknowledgments require final author review.

## References

1. NASA Open Science Data Repository. OSDR Biological Data API. <https://visualization.osdr.nasa.gov/biodata/api/>. Accessed July 28, 2026.
2. Lachmann A, Torre D, Keenan AB, et al. Massive mining of publicly available RNA-seq data from human and mouse. *Nature Communications*. 2018;9:1366. <https://doi.org/10.1038/s41467-018-03751-6>.
3. Viñas R, Andrés-Terré H, Liò P, Bryson K. Adversarial generation of gene expression data. *Bioinformatics*. 2022;38:730-737. <https://doi.org/10.1093/bioinformatics/btab035>.
4. Lacan A, André R, Sebag M, Hanczar B. In silico generation of gene expression profiles using diffusion models. *BMC Bioinformatics*. 2026. <https://doi.org/10.1186/s12859-026-06470-8>.
5. Litman E, Myers T, Agarwal V, et al. GeneJEPA: A predictive world model of the transcriptome. *bioRxiv*. 2025. <https://doi.org/10.1101/2025.10.14.682378>.
6. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300. <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>.
7. Milacic M, Beavers D, Conley P, et al. The Reactome Pathway Knowledgebase 2024. *Nucleic Acids Research*. 2024;52:D672-D678. <https://doi.org/10.1093/nar/gkad1025>.
8. Gridley DS, Mao XW, Stodieck LS, et al. Changes in mouse thymus and spleen after return from the STS-135 mission in space. *PLoS ONE*. 2013;8:e75097. <https://doi.org/10.1371/journal.pone.0075097>.
9. Horie K, Kato T, Kudo T, et al. Impact of spaceflight on the murine thymus and mitigation by exposure to artificial gravity during spaceflight. *Scientific Reports*. 2019;9:19866. <https://doi.org/10.1038/s41598-019-56432-9>.
10. Gambara G, Salanova M, Ciciliot S, et al. Gene expression profiling in slow-type calf soleus muscle of 30 days space-flown mice. *PLoS ONE*. 2017;12:e0169314. <https://doi.org/10.1371/journal.pone.0169314>.
11. Stein TP, Schluter MD, Grindeland RE, Moran MM, Baer LA, Wade CE. Rate controlling steps in fatty acid oxidation by unloaded rodent soleus muscle. *Journal of Gravitational Physiology*. 2002;9:P165-P166. PMID: 15002531.
12. Lotfollahi M, Rybakov S, Hrovatin K, et al. Biologically informed deep learning to query gene programs in single-cell atlases. *Nature Cell Biology*. 2023;25:337-350. <https://doi.org/10.1038/s41556-022-01072-x>.
13. Ho J, Jain A, Abbeel P. Denoising diffusion probabilistic models. *Advances in Neural Information Processing Systems*. 2020;33:6840-6851.
14. Song J, Meng C, Ermon S. Denoising diffusion implicit models. *International Conference on Learning Representations*. 2021. <https://arxiv.org/abs/2010.02502>.
15. Heusel M, Ramsauer H, Unterthiner T, Nessler B, Hochreiter S. GANs trained by a two time-scale update rule converge to a local Nash equilibrium. *Advances in Neural Information Processing Systems*. 2017;30.
16. Kynkäänniemi T, Karras T, Laine S, Lehtinen J, Aila T. Improved precision and recall metric for assessing generative models. *Advances in Neural Information Processing Systems*. 2019;32.
17. DerSimonian R, Laird N. Meta-analysis in clinical trials. *Controlled Clinical Trials*. 1986;7:177-188. <https://doi.org/10.1016/0197-2456(86)90046-2>.
18. Pezoldt J, Wiechers C, Erhard F, et al. Single-cell transcriptional profiling of splenic fibroblasts reveals subset-specific innate immune signatures in homeostasis and during viral infection. *Communications Biology*. 2021;4:1355. <https://doi.org/10.1038/s42003-021-02882-9>.
19. Cheng HW, Onder L, Novkovic M, et al. Origin and differentiation trajectories of fibroblastic reticular cells in the splenic white pulp. *Nature Communications*. 2019;10:1739. <https://doi.org/10.1038/s41467-019-09728-3>.
20. Kimura Y, Nakai Y, Ino Y, Akiyama T, Ryo A, Hirano H. Identification of gravity-responsive serum proteins in spaceflight mice using a quantitative proteomic approach with data-independent acquisition mass spectrometry. *Proteomics*. 2024;24:2300214. <https://doi.org/10.1002/pmic.202300214>.
21. Cope H, Elsborg J, Demharter S, et al. Transcriptomics analysis reveals molecular alterations underpinning spaceflight dermatology. *Communications Medicine*. 2024;4:106. <https://doi.org/10.1038/s43856-024-00532-9>.
22. Finch RH, Vitry G, Siew K, et al. Spaceflight causes strain-dependent gene expression changes in the kidneys of mice. *npj Microgravity*. 2025;11:11. <https://doi.org/10.1038/s41526-025-00465-0>.
23. Siew K, et al. Cosmic kidney disease: an integrated pan-omic, physiological and morphological study into spaceflight-induced renal dysfunction. *Nature Communications*. 2024;15:4923. <https://doi.org/10.1038/s41467-024-49212-1>.
24. Suzuki N, Iwamura Y, Nakai T, et al. Gene expression changes related to bone mineralization, blood pressure and lipid metabolism in mouse kidneys after space travel. *Kidney International*. 2022;101:92-105. <https://doi.org/10.1016/j.kint.2021.09.031>.
25. Jeong W, Kwon H, Park SK, Lee IS, Jho EH. Retinoic acid-induced protein 14 links mechanical forces to Hippo signaling. *EMBO Reports*. 2024;25:4033-4061. <https://doi.org/10.1038/s44319-024-00228-0>.
26. Fearnley GW, Young KA, Edgar JR, et al. The homophilic receptor PTPRK selectively dephosphorylates multiple junctional regulators to promote cell-cell adhesion. *eLife*. 2019;8:e44597. <https://doi.org/10.7554/eLife.44597>.
27. Liu X, Zhao Y, Gao J, et al. Elastic fiber homeostasis requires lysyl oxidase-like 1 protein. *Nature Genetics*. 2004;36:178-182. <https://doi.org/10.1038/ng1297>.
28. Machour FE, Abu-Zhayia ER, Awwad SW, et al. RBM6 splicing factor promotes homologous recombination repair of double-strand breaks and modulates sensitivity to chemotherapeutic drugs. *Nucleic Acids Research*. 2021;49:11708-11727. <https://doi.org/10.1093/nar/gkab976>.
29. Benbarche S, Bello Pineda JM, Baquero Galvis L, et al. GPATCH8 modulates mutant SF3B1 mis-splicing and pathogenicity in hematologic malignancies. *Molecular Cell*. 2024;84:1886-1903.e10. <https://doi.org/10.1016/j.molcel.2024.04.006>.
30. Wiltshire C, Matsushita M, Tsukada S, Gillespie DAF, May GHW. A new c-Jun N-terminal kinase (JNK)-interacting protein, Sab (SH3BP5), associates with mitochondria. *Biochemical Journal*. 2002;367:577-585. <https://doi.org/10.1042/BJ20020553>.
