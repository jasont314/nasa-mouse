<div class="title-page">

<h1>Synthetic-guided feature discovery in mouse spaceflight transcriptomics prioritizes thymic cell-cycle suppression and soleus metabolic remodeling</h1>

<p class="subtitle">Thymic proliferative suppression and soleus metabolic remodeling revealed with synthetic-guided analysis</p>

<p class="authors">Jason Trinh</p>

<p class="affiliation">Space Life Sciences Training Program, NASA Ames Research Center, Moffett Field, California, USA</p>

<p>Correspondence: jasontrinh@berkeley.edu</p>

<p class="draft-note"><strong>Manuscript draft for author review.</strong> The analyses and numerical results are frozen. Author list, acknowledgments, repository release URL, and archival DOI require final review before submission.</p>

</div>

## Abstract

**Background:** Mouse spaceflight studies provide access to tissues that cannot be sampled extensively from astronauts, but individual experiments are small and differ in design. We asked whether synthetic gene-expression models could help reveal reproducible flight-associated biology without treating generated profiles as new animals.

**Methods:** We assembled 1,610 mouse flight and ground-control bulk RNA-seq profiles through the NASA Open Science Data Repository API. A conditional diffusion model was pretrained on 17,244 tissue-diverse ARCHS4 mouse profiles and adapted to the OSDR studies. Generated profiles were used to guide gene selection, while biological effects were estimated from real samples and checked for consistency across studies. Full model, statistical, and sensitivity procedures are provided in the supplementary methods.

**Results:** The strongest result was in an independently held-out thymus study. Flight samples showed lower expression of eight mitotic genes, including *Cdk1*, *Ccnb1*, *Ccnb2*, *Birc5*, and *Ube2c*, together with lower APC/C, G2/M-checkpoint, and DNA-replication programs. The pattern was shared by wild-type and Nrf2-knockout mice and is consistent with reduced thymic proliferative renewal. Separating skeletal muscle by anatomical source exposed a complementary soleus response involving lower oxidative fuel handling, mitochondrial quality control, and slow-muscle identity, with higher *Tpm1* suggesting contractile remodeling. Lung produced mixed evidence; spleen, skin, and kidney yielded narrower hypotheses; liver and retina did not produce a coherent synthetic-guided result.

**Conclusions:** Synthetic-guided analysis added the most information when it helped prioritize genes rather than increasing the apparent sample size. The thymus result provides independent support for reduced proliferative renewal during spaceflight, while soleus identifies a complementary metabolic and contractile hypothesis that requires confirmation in a new study.

**Keywords:** spaceflight; bulk RNA-seq; diffusion model; synthetic data; thymus; soleus; NASA OSDR; ARCHS4; feature selection

## Introduction

Spaceflight affects immune, musculoskeletal, metabolic, and barrier tissues through a combination of microgravity, radiation, confinement, altered nutrition, stress, and mission-specific procedures. Mouse flight experiments provide tissue access that is unavailable in astronauts, but their transcriptomic interpretation is difficult. Individual studies are small, missions differ in strain and duration, and condition labels can be entangled with study, material, genotype, or collection protocol. Pooling samples without preserving these design variables can convert study effects into apparent flight biology.

The NASA Open Science Data Repository (OSDR) now exposes sample metadata and processed assay data through a queryable biological API [1]. This makes it possible to assemble a cross-study cohort without relying on a precombined raw HDF5 object and to retain accession-level provenance. Public reference resources offer a second opportunity. ARCHS4 uniformly processes a large fraction of public human and mouse RNA-seq data [2], providing tissue-diverse reference profiles for pretraining models that would be underdetermined on OSDR alone.

Deep generative models can learn high-dimensional expression distributions. Conditional WGAN-GP models have reproduced tissue and cancer properties in GTEx and TCGA [3]. More recently, Lacan and colleagues adapted denoising diffusion probabilistic and implicit models to bulk transcriptomics and reported strong gene-correlation, neighborhood, adversarial, and downstream classification metrics [4]. GeneJEPA instead learns masked-gene representations without reconstructing expression [5]. These approaches solve different problems: a generator can sample expression, whereas a representation learner needs an additional decoder or generative objective before it can do so.

Synthetic expression is commonly presented as a remedy for small sample size. Generated profiles, however, are not new biological replicates. We instead used synthetic expression as a model-derived view of the data that could prioritize genes and pathways for testing in real flight samples.

Our primary biological questions were whether this approach could clarify the tissue response to spaceflight, whether anatomical separation would expose muscle-specific responses hidden by pooling, and whether the resulting signals complemented pathway-level findings from expiMap.

The clearest findings arose in thymus and soleus. Thymus supplied an independent test in a study excluded from model adaptation. Soleus supplied a cross-study metabolic signal that is biologically coherent but still needs confirmation in a newly collected or fully unseen study. Results from the remaining tissues define the exploratory boundary of the approach.

![Study design and evidence ladder.](figures/figure_1_study_design.png)

<p class="caption"><strong>Figure 1. Synthetic-guided biological analysis.</strong> (A) A mouse tissue reference model was trained with ARCHS4 and adapted to API-derived OSDR flight and ground-control profiles. Generated expression prioritized genes, while biological effects were measured in real samples. (B) The analysis focused on tissue-wide responses, thymus, anatomically separated skeletal muscle, and hypotheses from other tissues.</p>

## Materials and methods

### Data sources

The OSDR Biological Data API was used to identify *Mus musculus* bulk RNA-seq assays with flight or ground-control labels [1]. Tissue and material names were harmonized while preserving study provenance. The resulting cohort contained 1,610 biological profiles from 75 accessions: 835 flight and 775 ground control. Full-transcriptome expression was converted to transcripts per million before selecting a 974-gene mouse landmark panel. No raw integrated OSDR H5 file was used.

The local ARCHS4 mouse resource contained 997,515 public RNA-seq profiles [2]. A healthy-preferred, tissue-balanced subset of 17,244 profiles spanning 20 tissue classes was used for model pretraining. GEO studies were kept intact when constructing training and held-out sets.

**Table 1. Data scope.**

| Source | Profiles used | Biological scope | Role |
|---|---:|---|---|
| ARCHS4 mouse v2.5 | 17,244 | 20 tissue classes | Mouse tissue pretraining |
| NASA OSDR API | 1,610 | 75 accessions; 835 flight and 775 ground control | Spaceflight analysis |

### Conditional expression model

We reproduced the bulk-expression diffusion architecture of Lacan et al. [4], based on denoising diffusion and implicit sampling [13,14]. The ARCHS4-pretrained model was adapted to OSDR while representing tissue, flight status, study, and material type as separate conditions. This allowed the same model to generate flight or ground-control profiles for represented biological and study contexts.

Generated profiles were checked for preservation of tissue identity, expression structure, diversity, and flight-related differences on withheld real profiles. The evaluation included correlation, neighborhood overlap, adversarial separability, and distributional distance [15,16]. The exact architecture, transformations, split construction, calibration, thresholds, and comparator-model results are reported in the supplementary methods.

### Synthetic-guided biological analysis

We first tested generated profiles as additional training rows and then as a guide for gene selection. Direct augmentation was not beneficial. The retained workflow instead used synthetic expression to rank candidate features and fitted the final predictive models using real profiles.

Flight-minus-ground effects were estimated within each OSDR study and then summarized across studies. Candidate genes were retained when their direction was consistent across the available real studies, with random-effects and multiple-testing procedures used as safeguards [6,17]. Reactome was used to group selected genes into biological processes [7]. Exact statistical definitions and sensitivity analyses are provided in the supplement.

An OSDR study was considered an independent test only when it had been excluded from model adaptation and feature-policy development. Signals observed across represented studies but lacking such a test are described as developmental or exploratory. Generative findings were also compared with the separate expiMap pathway analysis [12] to identify convergent and complementary tissue responses.

## Results

### Synthetic expression preserved tissue structure but worked best as a guide

The ARCHS4-pretrained model generated profiles that retained broad mouse tissue identity: a classifier learned from synthetic profiles predicted held-out real tissues nearly as well as one learned from real profiles. After adaptation to OSDR, generated profiles also reproduced the overall structure and diversity of the represented studies. Complete distributional results, repeated generations, denoising trajectories, and comparator-model screens are reported in Supplementary Figures S1-S4 and Supplementary Tables S1-S3.

Simply adding generated rows to the training data did not improve flight-versus-ground classification. The more useful strategy was to let the synthetic profiles influence which genes were considered, while fitting the final classifier and estimating biological effects only from real samples. This strategy transferred clearly to thymus and produced mixed results in lung. Detailed predictive comparisons are provided in Supplementary Figure S5 and Supplementary Table S4.

### Thymus shows lower proliferative renewal during spaceflight

OSD-457 provided the strongest result because it was excluded from model adaptation before testing. Synthetic-guided gene selection improved separation of flight and ground-control thymus samples, and the result held in both wild-type and Nrf2-knockout mice. Flight effects were closely aligned between the two genotypes.

The core genes *Birc5*, *Ccne2*, *Gmnn*, *Ube2c*, *Cdk1*, *Nusap1*, *Ccnb1*, and *Ccnb2* were lower in flight in both strata (Fig. 2A). Together, these genes regulate DNA replication, chromosome progression, mitotic entry, and completion of cell division. Reactome analysis connected them to APC/C-mediated cyclin degradation, G2/M checkpoints, DNA synthesis, and broader cell-cycle control (Fig. 2B).

This result is aligned with, but more specific than, prior thymus studies. STS-135 mouse thymus showed changes in cell-cycle and DNA-damage programs, including lower checkpoint-related expression [8]. A later ISS experiment reported marked thymus mass loss and partial artificial-gravity rescue of cell-cycle expression [9]. The current signature emphasizes mitotic completion and replication rather than acute apoptosis alone. Agreement between genotype strata suggests that the predictive signature is not confined to one Nrf2 background, but it does not prove Nrf2 independence.

Bulk thymus expression cannot distinguish lower transcription within proliferating thymocytes from loss or redistribution of proliferating cell populations. The defensible biological conclusion is lower abundance of a mitotic transcript program in flight, consistent with reduced proliferative renewal. Cell-resolved or histological confirmation is required to assign the effect to a cell-intrinsic mechanism.

![Thymus biology.](figures/figure_4_thymus_biology.png)

<p class="caption"><strong>Figure 2. Thymus response to spaceflight.</strong> (A) Flight-minus-ground effects in real OSD-457 profiles for eight cell-cycle genes that were lower in both wild-type and Nrf2-knockout mice. (B) The genes converge on a mitotic and DNA-replication process family.</p>

### Anatomical separation exposes a soleus-specific metabolic program

Aggregate skeletal muscle concealed substantial anatomical heterogeneity. We therefore examined extensor digitorum longus, gastrocnemius, quadriceps, soleus, and tibialis anterior separately. Soleus produced the clearest biological pattern: its selected genes showed consistent flight effects across three accessions and converged on related metabolic processes.

Seven synthetic-selected genes remained directionally consistent when each accession was examined in turn and agreed with the generated flight effect. *Bdh1*, *Bnip3*, *Mef2c*, *Ech1*, *Pxmp2*, and *Gmnn* were lower in flight, while *Tpm1* was higher (Fig. 3B). These genes were associated with mitochondrial protein turnover, fatty-acid oxidation, and lipid metabolism (Fig. 3C).

The pattern links ketone or lipid utilization (*Bdh1*, *Ech1*), mitochondrial quality control (*Bnip3*), slow oxidative muscle identity (*Mef2c*), peroxisomal transport (*Pxmp2*), and contractile remodeling (*Tpm1*). Prior 30-day spaceflight profiling of mouse soleus reported a slow-to-fast shift and broad changes in oxidative metabolism, PPAR signaling, and contractile genes [10]. Unloading studies have also reported reduced soleus fatty-acid oxidation [11]. The pathway is therefore literature-aligned, while the compact gene prioritization and peroxisome-mitochondria emphasis are exploratory refinements rather than wholly de novo biology.

Unlike thymus, soleus was represented during model development. Its cross-study consistency makes it a focused biological hypothesis, but an entirely unseen soleus study is still needed for independent confirmation.

![Soleus biology.](figures/figure_5_soleus_biology.png)

<p class="caption"><strong>Figure 3. Skeletal-muscle and soleus response.</strong> (A) Number of synthetic-prioritized genes with consistent real effects across studies in each anatomical muscle group. (B) Seven soleus genes with consistent real flight effects. (C) Their strongest shared biological processes center on mitochondrial turnover and lipid metabolism.</p>

### Other muscle groups provide narrower hypotheses

The other muscle groups produced narrower hypotheses. Quadriceps retained *Rbm6* without a broader pathway pattern. EDL showed lower *Abcc5*, *Lsm6*, *Polr2i*, and *Tsc22d3* in flight across its two accessions, suggesting RNA-processing and nuclear-receptor responses. Tibialis anterior showed higher *Cdkn1a*, *St3gal5*, *Cebpd*, *Pdhx*, and *Bnip3*, consistent with stress and metabolic remodeling, but synthetic guidance added little to an already strong real-data separation. Gastrocnemius did not produce a coherent retained gene set.

### The remaining tissue screen is complementary, not a superset of expiMap

The remaining tissues yielded candidate signals rather than complete biological stories. Lung highlighted cell-cycle, senescence, and PI3K/AKT-related genes, but the pattern varied by genotype and study. Spleen repeatedly selected *Igfbp3* without a larger coherent pathway. Skin suggested cell-cycle and DNA-repair responses, while kidney nominated the porphyrin-related genes *Hmox1* and *Alas1*. These candidates were not consistently supported across the available real studies. Liver and retina did not yield a coherent synthetic-guided gene or pathway pattern.

The evidence distribution therefore differs from the separate expiMap analysis. Both approaches prioritize thymus. expiMap produced its broadest pathway evidence in thymus, skin, spleen, and kidney, whereas the generative workflow added its clearest complementary result in soleus. The methods examine different biological representations: expiMap tests predefined pathway activity, while synthetic guidance helps prioritize individual genes and their shared processes.

![Tissue evidence matrix.](figures/figure_6_tissue_evidence.png)

<p class="caption"><strong>Figure 4. Biological findings across tissues.</strong> Thymus provides the strongest result, centered on lower cell proliferation. Soleus provides a complementary metabolic and contractile response. Lung is mixed; spleen, skin, and kidney provide narrower hypotheses; liver and retina do not provide a coherent synthetic-guided result.</p>

**Table 2. Biological interpretation by tissue.**

| Tissue | Main signal | Interpretation |
|---|---|---|
| Thymus | Lower mitotic genes; APC/C, G2/M, and DNA replication | Strongest result; supports reduced proliferative renewal |
| Soleus | Lower oxidative fuel handling and mitochondrial quality control; contractile remodeling | Coherent multi-study hypothesis requiring a new independent study |
| Lung | Cell-cycle, senescence, and PI3K/AKT candidates | Mixed across genotype and study |
| Spleen | *Igfbp3* | Isolated candidate without a shared process |
| Skin | Cell-cycle and DNA-repair candidates | Exploratory and complementary to expiMap |
| Kidney | *Hmox1* and *Alas1* porphyrin response | Exploratory |
| Liver, retina | No coherent retained pattern | No synthetic-guided biological conclusion |

## Discussion

### What synthetic data added

The generated profiles did not create stronger evidence simply by increasing the number of training rows. Their useful contribution was to expose combinations of genes that were less apparent in the small real cohorts. The synthetic model therefore acted as a feature-prior: it influenced what to examine, while the biological conclusions remained based on real flight and ground-control samples.

Thymus demonstrates this role most clearly. The synthetic-guided genes transferred to a study the model had not seen and converged on one interpretable cell-cycle process. Soleus demonstrates a second use: separating anatomically distinct muscle groups and combining synthetic prioritization with consistency across existing studies localized a broad muscle response to oxidative soleus biology. Neither result depends on counting synthetic profiles as animals.

### Thymus and soleus offer complementary biological stories

The thymus result refines established spaceflight immune biology. Prior work documents thymic involution and altered cell-cycle expression [8,9]. The current independently tested signature concentrates on cyclins, CDK1, UBE2C, BIRC5, NUSAP1, geminin, APC/C-mediated protein turnover, and G2/M control. Together these support lower proliferative renewal or a lower proportion of cycling thymocytes. Because the data are bulk, composition and cell-intrinsic regulation remain inseparable.

The soleus result addresses a different physiological axis. Weight-bearing slow muscle is especially sensitive to unloading. The selected genes describe lower oxidative substrate handling, altered mitochondrial turnover, reduced slow-muscle transcriptional identity, and contractile remodeling. This is compatible with known soleus atrophy and slow-to-fast transition [10,11], while nominating a compact set of genes for targeted validation. It is complementary to the expiMap result because synthetic guidance localized a broad fatty-acid-oxidation response to a smaller set observed consistently across soleus studies.

### Why the other tissues remain useful

The other tissues constrain the method's scope. Lung shows that a model can separate conditions without yielding one reproducible pathway. Spleen shows that an isolated gene can remain after the broader process signal disappears. Skin and kidney show that a pathway model and a generative feature model can prioritize different aspects of the same tissue response. Retina and liver show that a broadly realistic generator does not force every tissue contrast to become biologically informative.

These results still provide actionable hypotheses. Lung cell-cycle and PI3K/AKT candidates can be tested in a newly held-out study with genotype modeled prospectively. Spleen *Igfbp3* can be examined alongside the expiMap immune-program result, but should not be presented as a synthetic discovery by itself. Skin and kidney should be revisited only after prespecifying a transfer study and cell-composition audit. The main manuscript includes these results to prevent a success-only narrative.

### Limitations

Only one thymus study provided a fully unseen biological test, and the soleus result still lacks an equivalent independent study. Both findings should therefore be treated as focused evidence for follow-up rather than a comprehensive map of mouse spaceflight biology.

The analysis used a 974-gene landmark panel, so relevant genes outside that panel could not be discovered. Bulk tissue also cannot distinguish a transcriptional change within a cell type from a change in cell composition. The thymus result, for example, could reflect reduced expression in proliferating thymocytes, fewer proliferating thymocytes, or both.

Finally, the generator represents tissues and study contexts available during training. It should not be assumed to reproduce a new mission, strain, or sample-processing protocol without additional testing. Exact model limitations, sensitivity analyses, and statistical safeguards are documented in the supplementary methods.

## Conclusions

Synthetic expression was most informative as a guide to biological feature selection, not as a replacement for real animals. This approach independently supports a flight-lower thymus mitotic program and identifies a soleus response centered on oxidative metabolism, mitochondrial quality control, and contractile remodeling.

Thymus is the strongest synthetic-guided result and agrees with the strongest cross-method immune signal. Soleus is the main complementary finding beyond expiMap and now has a focused gene set for testing in a new study. Lung, spleen, skin, and kidney provide narrower hypotheses, while liver and retina set a clear negative boundary.

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
