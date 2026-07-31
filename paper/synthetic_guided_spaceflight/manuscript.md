<div class="title-page">

<h1>Cross-study synthetic-guided transcriptomics identifies thymic proliferative suppression and soleus metabolic remodeling in spaceflown mice</h1>

<p class="subtitle">Organ-specific flight responses across NASA OSDR mouse transcriptomes</p>

<p class="authors">Jason Trinh</p>

<p class="affiliation">Space Life Sciences Training Program, NASA Ames Research Center, Moffett Field, California, USA</p>

<p>Correspondence: jasontrinh@berkeley.edu</p>

<p class="draft-note"><strong>Research manuscript draft for author review.</strong> Author list, acknowledgments, repository release URL, and archival DOI require final review before submission.</p>

</div>

## Abstract

**Background:** Mouse spaceflight studies provide access to tissues that cannot be sampled extensively from astronauts, but individual experiments are small and differ in design. We asked whether synthetic gene-expression models could help identify reproducible, organ-specific flight responses without treating generated profiles as new animals.

**Methods:** We assembled 1,610 mouse flight and ground-control bulk RNA-seq profiles through the NASA Open Science Data Repository API. A conditional diffusion model was pretrained on 17,244 tissue-diverse ARCHS4 mouse profiles after excluding every GEO series linked to OSDR, with complete GEO series assigned to separate reference roles. The corrected backbone was adapted to OSDR with tissue, flight status, study, and material-type conditioning. Generated profiles guided gene selection, while flight effects were estimated exclusively from real samples within each study and then compared across studies.

**Results:** The strongest result was a leakage-corrected accession-held-out thymus test. Flight samples showed lower expression of eight mitotic genes, including *Cdk1*, *Ccnb1*, *Ccnb2*, *Birc5*, and *Ube2c*, together with lower APC/C, G2/M-checkpoint, and DNA-replication programs. Synthetic-guided balanced accuracy increased from 0.500 to 0.833 and AUROC from 0.840 to 0.979. The pattern was shared by wild-type and Nrf2-knockout mice and is consistent with reduced thymic proliferative renewal. Anatomical separation of skeletal muscle exposed a soleus response involving lower *Bdh1*, *Ech1*, *Bnip3*, and *Decr1* with higher *Tpm1*. All five genes were reinforced by real-only and synthetic-guided selection and converged on mitochondrial protein turnover, fatty-acid oxidation, and lipid metabolism. Kidney provided a secondary result: synthetic guidance promoted flight-higher *Inpp4b*, which passed leave-one-accession-out sensitivity, and reinforced flight-higher *Slc37a4*. The corrected spleen screen instead promoted *Rai14*, *Ptprk*, and *Myl9* and reinforced *Loxl1*; none passed leave-one-accession-out FDR, and *Igfbp3* remained a strong real-data association but was not selected by the corrected synthetic workflow. Across tissues, 49 BH-FDR associations also entered stable synthetic-informed selection: 26 promoted and 23 reinforced.

**Conclusions:** Synthetic-guided analysis added the most information when it prioritized genes for testing in real samples rather than increasing the apparent sample size. The corrected findings support reduced thymic proliferative renewal, reinforce a soleus mitochondrial and contractile response, and nominate focused renal and spleen candidates for prospective testing.

**Keywords:** spaceflight; bulk RNA-seq; synthetic data; thymus; soleus; kidney; NASA OSDR; ARCHS4

## Introduction

Spaceflight affects immune, musculoskeletal, metabolic, and barrier tissues through a combination of microgravity, radiation, confinement, altered nutrition, stress, and mission-specific procedures. Mouse flight experiments provide tissue access that is unavailable in astronauts, but their transcriptomic interpretation is difficult. Individual studies are small, missions differ in strain and duration, and condition labels can be entangled with study, material, genotype, or collection protocol. Pooling samples without preserving these design variables can convert study effects into apparent flight biology.

The NASA Open Science Data Repository (OSDR) now exposes sample metadata and processed assay data through a queryable biological API [1]. This makes it possible to assemble a cross-study cohort without relying on a precombined raw HDF5 object and to retain accession-level provenance. Public reference resources offer a second opportunity. ARCHS4 uniformly processes a large fraction of public human and mouse RNA-seq data [2], providing tissue-diverse reference profiles for pretraining models that would be underdetermined on OSDR alone.

Deep generative models can learn high-dimensional expression distributions. Conditional WGAN-GP models have reproduced tissue and cancer properties in GTEx and TCGA [3]. More recently, Lacan and colleagues adapted denoising diffusion probabilistic and implicit models to bulk transcriptomics and reported strong gene-correlation, neighborhood, adversarial, and downstream classification metrics [4]. GeneJEPA instead learns masked-gene representations without reconstructing expression [5]. These approaches solve different problems: a generator can sample expression, whereas a representation learner needs an additional decoder or generative objective before it can do so.

Synthetic expression is commonly presented as a remedy for small sample size. Generated profiles, however, are not new biological replicates. We instead used synthetic expression as a model-derived view of the data that could prioritize genes and pathways for testing in real flight samples.

Our primary biological questions were whether this approach could clarify the tissue response to spaceflight, whether anatomical separation would expose muscle-specific responses hidden by pooling, and whether the resulting signals complemented pathway-level findings from expiMap.

The clearest findings arose in thymus and soleus. Thymus supplied a leakage-corrected test in a study excluded from ARCHS4 pretraining, OSDR adaptation, and feature-policy development. Because that outcome had been examined before the leakage correction, it is a retrospective sensitivity analysis rather than a prospective untouched confirmation. Soleus supplied a cross-study metabolic signal that is biologically coherent but still needs confirmation in a newly collected or fully unseen study. Kidney, pooled skeletal muscle, spleen, skin, and adrenal gland supplied secondary developmental candidates. Results from the remaining tissues define the exploratory boundary of the approach.

![Study design and evidence ladder.](figures/figure_1_study_design.png)

<p class="caption"><strong>Figure 1. Synthetic-guided biological analysis.</strong> (A) A mouse tissue reference model was trained with ARCHS4 and adapted to API-derived OSDR flight and ground-control profiles. Generated expression prioritized genes, while biological effects were measured in real samples. (B) The analysis focused on organ-specific responses in thymus, anatomically separated skeletal muscle, spleen, and the remaining tissue cohort.</p>

## Materials and methods

### Data sources

The OSDR Biological Data API was used to identify *Mus musculus* bulk RNA-seq assays with flight or ground-control labels [1]. Tissue and material names were harmonized while preserving study provenance. The resulting cohort contained 1,610 biological profiles from 75 accessions: 835 flight and 775 ground control. Full-transcriptome expression was converted to transcripts per million before selecting a 974-gene mouse landmark panel. No raw integrated OSDR H5 file was used.

The local ARCHS4 mouse resource contained 997,515 public RNA-seq profiles [2]. All nine GEO series linked from the API-derived OSDR metadata were excluded before reference selection, removing 108 otherwise eligible profiles and preventing cross-resource sample overlap. A healthy-preferred, tissue-balanced subset of 17,244 profiles spanning 20 tissue classes was then used for model pretraining. Complete GEO series were assigned to one reference role, producing 10,150 training, 2,466 validation, and 4,628 test profiles with no series overlap. The corrected backbone was used for both the held-out confirmation and the broad development screens.

**Table 1. Data scope.**

| Source | Profiles used | Biological scope | Role |
|---|---:|---|---|
| ARCHS4 mouse v2.5 | 17,244 | 20 tissue classes | Mouse tissue pretraining |
| NASA OSDR API | 1,610 | 75 accessions; 835 flight and 775 ground control | Spaceflight analysis |

### Conditional expression model

We reproduced the bulk-expression diffusion architecture of Lacan et al. [4], based on denoising diffusion and implicit sampling [13,14]. The corrected ARCHS4-pretrained model was adapted to OSDR while representing tissue, flight status, study, and material type as separate conditions. This allowed the same model lineage to generate flight or ground-control profiles for represented biological and study contexts. A designated accession-excluded adaptation was used for the fixed lung/thymus test, while a separate corrected broad adaptation supplied the repeated all-tissue development screens.

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
| 2 | Repeated development screen plus real-data BH FDR | Which BH-significant genes were repeatedly promoted or reinforced by synthetic guidance? | 49 tissue-gene results: 26 promoted and 23 reinforced | Cross-study developmental hypotheses; not independent synthetic generalization |
| 3 | Complete real-data random-effects screen | Which panel genes show a pooled FLT-GC association? | 459 tissue-gene results across canonical and muscle-group analyses | Conventional panel-level association; does not show that synthetic data added information |

Signals observed across represented studies but lacking a Tier 1 test are described as developmental or exploratory. Generative findings were also compared with the separate expiMap pathway analysis [12] to identify convergent and complementary tissue responses.

For a secondary spleen comparison, study-level *Igfbp3* effects were checked against full-transcriptome TPM values. Healthy sorted-cell and single-cell references, GSE156162 and E-MTAB-7703, were then examined to identify the normal spleen populations expressing *Igfbp3* [18,19]. The corrected model did not select *Igfbp3*, so this analysis is reported as real-data context rather than a synthetic-guided finding.

## Results

### Synthetic expression preserved tissue structure and prioritized biological features

The leakage-corrected ARCHS4 model generated profiles that retained broad mouse tissue identity: classifiers trained on real and generated reference profiles both achieved tissue balanced accuracy 0.781 on the complete-series test split. Adversarial accuracy was 0.515, precision was 0.951, and recall was 0.890. However, reference gene-correlation-matrix agreement was 0.878 and did not pass the prespecified 0.952 gate.

After broad OSDR adaptation, all four locked within-study repeats passed the finite-sample fidelity rule. Mean gene-correlation agreement was 0.974, precision 0.997, recall 0.996, F1 0.997, adversarial accuracy 0.475, and FD relative to the real-split P95 was 0.074. Pooled FLT-GC effect recovery passed three of four repeats and skeletal-muscle accession recovery passed all four. The preceding validation screen was less consistent: correlation narrowly missed its finite-sample floor and muscle accession-effect recovery failed. The generator is therefore suitable for the declared within-study feature-guidance analysis, not as a universal substitute for real expression or as evidence of unseen-study generation. Complete distributional results, denoising trajectories, and comparator-model screens are reported in Supplementary Figures S1-S4 and Supplementary Tables S1-S3.

Simply adding generated rows to the training data did not improve flight-versus-ground classification. On the locked test, real-only balanced accuracy was 0.754, compared with 0.695 for synthetic-only training and 0.737 for real-plus-synthetic training. The more useful strategy was to let the synthetic profiles influence which genes were considered, while fitting the final classifier and estimating biological effects only from real samples. This strategy transferred clearly to thymus and produced mixed results in lung. Detailed predictive comparisons are provided in Supplementary Figure S5 and Supplementary Table S4.

The Tier 3 real-data screen yielded 202 BH-FDR tissue-gene associations across 10 of 22 canonical tissue analyses and 257 across the five anatomical muscle groups. These are tissue-gene results rather than 459 unique or independent discoveries: a gene can be significant in more than one tissue, and pooled skeletal muscle overlaps its anatomical subgroups. Forty-nine results also entered Tier 2 stable synthetic-informed selection: 26 were synthetic-promoted and 23 were reinforced by both real-only and synthetic-guided selection. Synthetic labels were assigned only where the selected generated arm passed the prespecified metric gate; quadriceps, EDL, and liver therefore contribute real-data associations but no corrected synthetic-informed claims. The complete set, including associations with mixed study directions, remains available in Supplementary Tables S16-S18.

### Thymus shows lower proliferative renewal during spaceflight

OSD-457 provided the strongest result after a leakage audit found that its GEO series had been present in the original ARCHS4 reference. The reference was retrained from scratch after excluding every OSDR-linked GEO series and assigning complete GEO series to separate roles. OSD-457 was also excluded from OSDR adaptation and feature-policy development. In this corrected test, synthetic-guided gene selection improved separation of flight and ground-control thymus samples from balanced accuracy 0.500 to 0.833 and AUROC 0.840 to 0.979. The result held in both wild-type and Nrf2-knockout mice. Flight effects were closely aligned between the two genotypes.

The core genes *Birc5*, *Ccne2*, *Gmnn*, *Ube2c*, *Cdk1*, *Nusap1*, *Ccnb1*, and *Ccnb2* were lower in flight in both strata (Fig. 2A). These were not genes absent from the real-data analysis: all eight were among its 26 highest-ranked candidates, including four among the top ten. Synthetic guidance instead assembled a broader, coherent cell-cycle panel that transferred more effectively to the study-excluded test. Together, these genes regulate DNA replication, chromosome progression, mitotic entry, and completion of cell division. Reactome analysis connected them to APC/C-mediated cyclin degradation, G2/M checkpoints, DNA synthesis, and broader cell-cycle control (Fig. 2B).

The leakage-corrected Tier 1 result is therefore interpreted as panel-level reinforcement and organization, not de novo gene discovery. In the corrected Tier 2 repeated screen, *Birc5*, *Cdk1*, *Nusap1*, *Ccne2*, and *Ccnb2* crossed the stable-selection threshold only with synthetic guidance; *Gmnn* and *Ube2c* were reinforced by both arms; and *Ccnb1* remained BH-significant but did not cross either stability threshold. These labels describe repeated selection behavior and do not override the stronger held-out result. The complete cross-analysis mapping is provided in Supplementary Table S19.

This result is aligned with, but more specific than, prior thymus studies. STS-135 mouse thymus showed changes in cell-cycle and DNA-damage programs, including lower checkpoint-related expression [8]. A later ISS experiment reported marked thymus mass loss and partial artificial-gravity rescue of cell-cycle expression [9]. The current signature emphasizes mitotic completion and replication rather than acute apoptosis alone. Agreement between genotype strata suggests that the predictive signature is not confined to one Nrf2 background, but it does not prove Nrf2 independence.

Bulk thymus expression cannot distinguish lower transcription within proliferating thymocytes from loss or redistribution of proliferating cell populations. The defensible biological conclusion is lower abundance of a mitotic transcript program in flight, consistent with reduced proliferative renewal. Cell-resolved or histological confirmation is required to assign the effect to a cell-intrinsic mechanism.

![Thymus biology.](figures/figure_4_thymus_biology.png)

<p class="caption"><strong>Figure 2. Thymus response to spaceflight.</strong> (A) Flight-minus-ground effects in real OSD-457 profiles for eight cell-cycle genes that were lower in both wild-type and Nrf2-knockout mice. (B) The genes converge on a mitotic and DNA-replication process family.</p>

### Anatomical separation exposes a soleus-specific metabolic program

Aggregate skeletal muscle concealed substantial anatomical heterogeneity. We therefore examined extensor digitorum longus, gastrocnemius, quadriceps, soleus, and tibialis anterior separately. Soleus produced the clearest biological pattern: its selected genes showed consistent flight effects across three accessions and converged on related metabolic processes.

The corrected soleus screen selected real-plus-generated training. Balanced accuracy increased from 0.925 to 0.963, AUROC remained 0.980, and average precision increased from 0.980 to 0.986. Five BH-FDR genes were stable in both real-only and synthetic-guided selection: *Bdh1*, *Ech1*, *Bnip3*, and *Decr1* were lower in flight in all three accessions, while *Tpm1* was higher (Fig. 3B). Four genes passed leave-one-accession-out FDR; *Decr1* did not.

The corrected model did not promote a soleus gene absent from stable real-only selection. Its contribution was reinforcement of a coherent existing pattern. Reactome connected the retained genes to mitochondrial protein turnover, mitochondrial fatty-acid beta-oxidation, and lipid metabolism (Fig. 3C). *Bdh1* and *Ech1* support oxidative substrate handling, *Bnip3* supports mitochondrial quality control, *Decr1* contributes to unsaturated fatty-acid oxidation, and *Tpm1* suggests contractile remodeling.

Prior 30-day spaceflight profiling of mouse soleus reported a slow-to-fast shift and broad changes in oxidative metabolism, PPAR signaling, and contractile genes [10]. Unloading studies have also reported reduced soleus fatty-acid oxidation [11]. The corrected result is therefore literature-aligned panel reinforcement rather than de novo gene discovery.

Unlike thymus, soleus was represented during model development. Its cross-study consistency makes it a focused biological hypothesis, but an entirely unseen soleus study is still needed for independent confirmation.

![Soleus biology.](figures/figure_5_soleus_biology.png)

<p class="caption"><strong>Figure 3. Skeletal-muscle and soleus response.</strong> (A) Number of corrected synthetic-informed genes that also pass real leave-one-accession-out FDR in each anatomical muscle group. (B) Five reinforced soleus BH-FDR genes with consistent real flight effects. (C) Their strongest shared biological processes center on mitochondrial turnover and lipid metabolism.</p>

### Other muscle groups provide narrower hypotheses

The pooled skeletal-muscle screen selected synthetic-guided feature ranking and improved balanced accuracy, AUROC, and average precision by 0.071, 0.036, and 0.037. Twelve BH-FDR genes were synthetic-informed, nine of which passed leave-one-accession-out sensitivity. The associated gene sets included interferon signaling and sialic-acid metabolism. Because pooled muscle combines anatomical groups with distinct physiology, this result is interpreted alongside, rather than instead of, the soleus analysis.

Among the remaining anatomical groups, tibialis anterior retained a real-plus-generated arm and gastrocnemius retained a low-weight guided arm. Tibialis anterior reinforced *Cdkn1a*, *St3gal5*, and *Bnip3* and promoted *Cebpd*; gastrocnemius promoted *Fhl2* and *Nfkbia*. None of these genes passed leave-one-accession-out FDR. EDL and quadriceps selected real-only arms, so their BH-FDR genes, including quadriceps *Rbm6*, are retained in the complete real-data inventory but are not attributed to the corrected synthetic workflow.

### Corrected spleen screen nominates adhesion and cytoskeletal genes

The corrected spleen screen selected real-plus-generated training. Balanced accuracy increased from 0.517 to 0.648, AUROC from 0.540 to 0.704, and average precision from 0.589 to 0.749. Synthetic guidance promoted flight-higher *Rai14*, *Ptprk*, and *Myl9* and reinforced flight-higher *Loxl1*. *Rai14*, *Ptprk*, and *Loxl1* had the same direction in all six studies; *Myl9* agreed in five of six. None passed leave-one-accession-out FDR, and the set did not produce significant Reactome enrichment.

The four genes provide a tentative structural hypothesis rather than a pathway claim. RAI14 links actin-associated mechanosensing to Hippo signaling, PTPRK regulates cell-cell junctions, MYL9 contributes to actomyosin contractility, and LOXL1 supports elastic extracellular-matrix maintenance [25-27]. Their shared direction is compatible with altered adhesion, mechanics, or tissue architecture, but the present bulk data cannot identify a common source cell or establish one mechanism.

*Igfbp3* remains a notable but separate real-data association. It was higher in flight in all six eligible spleen studies, with a random-effects FDR of `1.76e-9`; the largest FDR after removing one study at a time was `0.00385`. A full-transcriptome TPM sensitivity calculation showed 9% to 63% higher group means in flight. Healthy references localized baseline *Igfbp3* mainly to fibroblastic reticular, collagen-producing, and perivascular stromal populations [18,19], and spaceflight proteomics has reported increased IGFBP3 in serum and femur under microgravity relative to onboard artificial gravity [20]. However, the corrected synthetic workflow did not select *Igfbp3*. It is therefore retained as a real-data stromal hypothesis, not evidence that synthetic guidance discovered or reinforced the gene.

### Other organ responses were heterogeneous

Lung produced large repeated development-screen gains and cell-cycle, senescence, and PI3K/AKT-related enrichment, but no selected gene passed the real-data BH-FDR screen. More importantly, the fixed OSD-900 test lost balanced accuracy and failed its validation gate. Lung therefore demonstrates why within-study separation and pathway enrichment cannot substitute for study-held-out transfer.

Skin selected real-plus-generated training and promoted flight-higher *Plscr1*. *Plscr1* passed BH FDR and had the same direction in all six studies, although it did not pass leave-one-accession-out FDR. Broader selected sets were enriched for G1/S and DNA-repair processes, consistent with prior OSDR skin analyses reporting strong mission, strain, recovery-interval, and anatomical-site dependence [21]. This is a literature-aligned developmental candidate rather than a confirmed transferable signal.

Kidney produced the clearest secondary gene-level result. *Slc37a4* was reinforced by both selection arms and was higher in flight in all six kidney accessions. Synthetic guidance additionally promoted flight-higher *Inpp4b*; it passed BH FDR and remained significant in leave-one-accession-out analysis, although one of six accession effects was negative. The selected pair did not yield significant Reactome enrichment. *Slc37a4* supports renal glucose handling, while *Inpp4b* nominates phosphoinositide signaling. Both require prospective study-level confirmation and should be interpreted alongside prior reports of strain-dependent lipid, ECM, TGF-beta, Wnt, and nephron-remodeling responses [22-24].

Adrenal gland produced flight-lower synthetic-promoted *Psmb8* and reinforced *Tspan4*, each with the same direction in all three studies but neither passing leave-one-accession-out FDR. Retina showed repeated predictive gains and selected-set enrichment but no corrected synthetic-informed BH-FDR gene. Liver selected the real-only arm: its 19 real-data BH-FDR genes remain in the complete inventory, but none support a corrected synthetic-guided claim. Quadriceps and EDL likewise selected real-only arms.

The evidence distribution therefore differs from the separate expiMap analysis. Both approaches prioritize thymus. expiMap produced its broadest pathway evidence in thymus, skin, spleen, and kidney, whereas the corrected generative workflow added its clearest complementary result in soleus and a focused kidney pair. The methods examine different biological representations: expiMap tests predefined pathway activity, while synthetic guidance helps prioritize individual genes and their shared processes.

![Corrected tissue evidence hierarchy.](figures/figure_6_tissue_evidence.png)

<p class="caption"><strong>Figure 4. Corrected-model tissue evidence.</strong> (A) Repeated development-screen changes relative to real-only models for selected tissues. (B) BH-FDR genes promoted or reinforced by the corrected synthetic workflow; synthetic labels are suppressed where the generated arm failed its metric gate. (C) Thymus remains the strongest result, soleus and pooled muscle provide the clearest process-level complement, kidney supplies a focused secondary pair, and spleen, skin, and adrenal gland remain developmental.</p>

**Table 3. Selected synthetic-guided biological interpretations by tissue. Complete real-data BH-FDR results are in Supplementary Tables S17-S18.**

| Tissue | Main signal | Interpretation |
|---|---|---|
| Thymus | Lower mitotic genes; APC/C, G2/M, and DNA replication | Strongest result; supports reduced proliferative renewal |
| Soleus | Reinforced FLT-lower *Bdh1*, *Ech1*, *Bnip3*, and *Decr1*; FLT-higher *Tpm1* | Coherent mitochondrial and lipid-metabolism program; no corrected promoted gene |
| Skeletal muscle, pooled | 12 synthetic-informed BH-FDR genes; interferon and sialic-acid terms | Developmental complement to the anatomically specific soleus result |
| Kidney | Promoted *Inpp4b* and reinforced *Slc37a4*, both FLT-higher | Focused renal metabolic-signaling hypothesis; *Inpp4b* passes LOO FDR |
| Spleen | Promoted *Rai14*, *Ptprk*, and *Myl9*; reinforced *Loxl1* | Adhesion/cytoskeletal hypothesis; no gene passes LOO FDR and no significant pathway |
| Skin | Promoted FLT-higher *Plscr1* plus cell-cycle/DNA-repair enrichment | Literature-aligned developmental candidate; not LOO-stable |
| Adrenal gland | Promoted *Psmb8* and reinforced *Tspan4*, both FLT-lower | Three-study developmental candidate; neither passes LOO FDR |
| Tibialis anterior | Reinforced *Cdkn1a*, *St3gal5*, and *Bnip3*; promoted *Cebpd* | Exploratory stress and metabolic response; no LOO-stable selected gene |
| Gastrocnemius | Promoted *Fhl2* and *Nfkbia* | Exploratory two-gene result without significant pathway or LOO stability |
| Lung | Within-study pathway candidates but failed OSD-900 transfer | No retained study-held-out synthetic-guided result |
| Retina | Predictive and pathway gains without a synthetic-informed BH-FDR gene | Exploratory gene/pathway mismatch |
| Quadriceps, EDL, liver | Corrected screen retained real-only arms | Real-data BH-FDR associations are not synthetic-guided findings |

## Discussion

### What synthetic data added

The generated profiles did not create stronger evidence simply by increasing the number of training rows. Their useful contribution was either to reinforce coherent combinations already present in the real data or to promote candidates that were not stable under real-only selection. The synthetic model therefore acted as a feature-prior: it influenced what to examine, while the biological conclusions remained based on real flight and ground-control samples.

Thymus demonstrates panel-level organization: its central genes already ranked strongly in real data, but the synthetic-guided panel transferred after the study was removed from both reference pretraining and OSDR adaptation and converged on one interpretable cell-cycle process. The Tier 2 labels describe which genes crossed repeated selection thresholds and are not claims of biological novelty.

Soleus demonstrates reinforcement rather than discovery. The corrected synthetic arm retained the same five BH-FDR genes as stable real-only selection and sharpened balanced accuracy without improving AUROC. The value is the convergence of independently selected genes on mitochondrial turnover and fatty-acid metabolism, not nomination of a new soleus gene.

Kidney provides the clearest corrected promotion result outside thymus. Synthetic guidance added *Inpp4b* to reinforced *Slc37a4*, and the real-data association for *Inpp4b* survived leave-one-accession-out analysis. Spleen provides a weaker multi-gene example: *Rai14*, *Ptprk*, and *Myl9* were promoted and *Loxl1* reinforced, but none passed LOO FDR and no pathway was significant. The previously emphasized *Igfbp3* association remains biologically interesting but does not demonstrate synthetic contribution. Quadriceps and EDL selected real-only arms and therefore no longer support synthetic-guided claims.

### Thymus and soleus define the strongest organ responses

The thymus result refines established spaceflight immune biology. Prior work documents thymic involution and altered cell-cycle expression [8,9]. The current leakage-corrected signature concentrates on cyclins, CDK1, UBE2C, BIRC5, NUSAP1, geminin, APC/C-mediated protein turnover, and G2/M control. Together these support lower proliferative renewal or a lower proportion of cycling thymocytes. Because the data are bulk, composition and cell-intrinsic regulation remain inseparable.

The soleus result addresses a different physiological axis. Weight-bearing slow muscle is especially sensitive to unloading. Lower *Bdh1*, *Ech1*, *Bnip3*, and *Decr1* together with higher *Tpm1* describe reduced oxidative substrate handling, altered mitochondrial quality control, and contractile remodeling. This is compatible with known soleus atrophy, altered lipid metabolism, and slow-to-fast transition [10,11]. It is complementary to expiMap because synthetic guidance reinforced a compact fatty-acid-oxidation panel observed consistently across soleus studies.

### Why the other tissues remain useful

The other tissues constrain the method's scope. Kidney combines one promoted and one reinforced gene but lacks pathway enrichment. Spleen combines strong within-study predictive gains with four BH-FDR genes, yet none survives the stricter LOO criterion. Skin promotes one directionally consistent gene while broader pathway enrichment remains heterogeneous. Adrenal gland has only three studies. Retina has predictive and pathway gains without a synthetic-informed BH-FDR gene. Lung fails the study-held-out test despite strong development-screen separation. Liver, quadriceps, and EDL select real-only arms.

These results define direct follow-up priorities. A new thymus study should test the fixed eight-gene panel prospectively. Soleus experiments should quantify the five-gene mitochondrial and contractile panel together with fatty-acid oxidation. Kidney should test *Inpp4b* and *Slc37a4* jointly in an unseen mission. Spleen should test *Rai14*, *Ptprk*, *Myl9*, and *Loxl1* in stroma-preserving assays, while *Igfbp3* should be evaluated separately as a real-data-derived stromal hypothesis. Lung requires a prospectively genotype-stratified study rather than further tuning on OSD-900.

### Limitations

No tissue currently has a completely prospective untouched biological test. OSD-457 was excluded from the leakage-corrected reference, adaptation, and feature-policy development, but its outcome had already been examined before the correction. It therefore provides strong retrospective sensitivity evidence rather than prospective confirmation. All other BH-FDR findings were identified within the broader model-development domain and still lack equivalent study-excluded testing. Agreement across represented real studies supports focused follow-up but does not establish transfer to a new mission or study.

The analysis used a 974-gene landmark panel, so relevant genes outside that panel could not be discovered. BH FDR was controlled separately within each declared tissue family, not once across every tissue-gene combination in the project. Direction disagreement and high heterogeneity do not remove a BH-significant association, but they reduce confidence that its pooled effect represents a common response across missions. Bulk tissue also cannot distinguish a transcriptional change within a cell type from a change in cell composition. The thymus result could reflect lower expression in proliferating thymocytes, fewer proliferating thymocytes, or both. The spleen structural candidates likewise cannot distinguish altered expression from altered stromal or contractile-cell abundance.

Finally, the generator represents tissues and study contexts available during training. The corrected ARCHS4 reference failed its strict correlation-structure gate. After OSDR adaptation, the locked within-study test passed its finite-sample fidelity rule, but the preceding validation screen missed its correlation floor and muscle accession-effect gate. The broad tissue screens measure interpolation within represented studies and remain hypothesis-generating. The model should not be assumed to reproduce a new mission, strain, or sample-processing protocol without additional testing. Exact model limitations, sensitivity analyses, and statistical safeguards are documented in the supplementary methods.

## Conclusions

Synthetic expression was most informative as a guide to biological feature selection, not as a replacement for real animals. The leakage-corrected analysis supports a flight-lower thymus mitotic program and reinforces a soleus response centered on mitochondrial lipid metabolism and contractile remodeling. Kidney supplies the clearest secondary promoted gene, *Inpp4b*, together with reinforced *Slc37a4*.

Thymus remains the strongest result because it survived complete OSDR-series removal from reference pretraining and accession exclusion from adaptation. Soleus supplies the clearest complementary process-level result. Pooled skeletal muscle, kidney, spleen, skin, adrenal gland, gastrocnemius, and tibialis anterior provide developmental hypotheses requiring unseen-study validation. Quadriceps, EDL, and liver no longer support corrected synthetic-guided claims because their real-only arms were retained. The lung result did not transfer to OSD-900. These revisions narrow the biological claims while preserving a reproducible framework for using synthetic expression as a feature prior.

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
