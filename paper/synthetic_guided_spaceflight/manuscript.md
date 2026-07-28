<div class="title-page">

<h1>Synthetic-guided feature discovery in mouse spaceflight transcriptomics prioritizes thymic cell-cycle suppression and soleus metabolic remodeling</h1>

<p class="subtitle">A conditional diffusion benchmark with independently held-out and cross-accession validation</p>

<p class="authors">Jason Trinh</p>

<p class="affiliation">Space Life Sciences Training Program, NASA Ames Research Center, Moffett Field, California, USA</p>

<p>Correspondence: jasontrinh@berkeley.edu</p>

<p class="draft-note"><strong>Manuscript draft for author review.</strong> The analyses and numerical results are frozen. Author list, acknowledgments, repository release URL, and archival DOI require final review before submission.</p>

</div>

## Abstract

**Background:** Bulk RNA-sequencing studies of mouse spaceflight contain strong tissue, study, and material effects, while individual flight-versus-ground-control comparisons are often small. Synthetic expression could regularize feature discovery, but generated profiles cannot be counted as new biological replicates and high marginal fidelity does not establish recovery of a biological contrast.

**Methods:** We queried the NASA Open Science Data Repository (OSDR) Biological Data API for *Mus musculus* bulk RNA-seq profiles labeled flight or ground control. A paper-parity denoising diffusion implicit model (DDIM) was pretrained on 17,244 healthy-preferred ARCHS4 mouse profiles spanning 20 tissue classes and 974 landmark genes. The model was adapted to OSDR with factorized tissue, flight status, accession, and material-type conditioning. Distribution fidelity was evaluated on a locked 293-profile test set over four prespecified generation seeds. Synthetic data were then tested as direct augmentation and as a feature-selection prior. Biological support always came from real within-accession flight-minus-ground effects, random-effects meta-analysis, false-discovery-rate control, and leave-one-accession-out sensitivity. A fixed feature-guidance policy was tested in OSD-900 lung and OSD-457 thymus after both accessions were excluded from generator adaptation.

**Results:** The ARCHS4 DDIM achieved synthetic-train-to-real-test tissue balanced accuracy 0.869 versus 0.895 for real training. On the locked OSDR test, mean gene-correlation agreement, precision, recall, F1, adversarial accuracy, and normalized Frechet distance were 0.977, 0.998, 0.997, 0.997, 0.458, and 0.075. Direct synthetic augmentation did not improve pooled flight classification. Synthetic-guided ranking improved held-out thymus balanced accuracy from 0.500 to 0.833 and AUROC from 0.840 to 0.972. Eight flight-lower mitotic genes were concordant across wild-type and Nrf2-knockout strata, with Reactome enrichment for APC/C regulation, G2/M checkpoints, and DNA replication. Lung transfer was mixed. In cross-accession development analysis, anatomical separation of skeletal muscle identified seven synthetic-selected, real-data leave-one-accession-out-stable soleus genes and enrichment for mitochondrial fatty-acid oxidation and protein turnover. Other tissues produced exploratory or negative results.

**Conclusions:** Conditional diffusion added biological information most defensibly by guiding feature selection while final inference remained anchored to real samples. Thymus provides independent support for lower proliferative renewal during spaceflight. Soleus supplies a complementary oxidative-metabolism hypothesis that requires an unseen-accession confirmation. Synthetic profiles should be treated as model-derived priors, not additional animals.

**Keywords:** spaceflight; bulk RNA-seq; diffusion model; synthetic data; thymus; soleus; NASA OSDR; ARCHS4; feature selection

## Introduction

Spaceflight affects immune, musculoskeletal, metabolic, and barrier tissues through a combination of microgravity, radiation, confinement, altered nutrition, stress, and mission-specific procedures. Mouse flight experiments provide tissue access that is unavailable in astronauts, but their transcriptomic interpretation is difficult. Individual studies are small, missions differ in strain and duration, and condition labels can be entangled with study, material, genotype, or collection protocol. Pooling samples without preserving these design variables can convert study effects into apparent flight biology.

The NASA Open Science Data Repository (OSDR) now exposes sample metadata and processed assay data through a queryable biological API [1]. This makes it possible to assemble a cross-study cohort without relying on a precombined raw HDF5 object and to retain accession-level provenance. Public reference resources offer a second opportunity. ARCHS4 uniformly processes a large fraction of public human and mouse RNA-seq data [2], providing tissue-diverse reference profiles for pretraining models that would be underdetermined on OSDR alone.

Deep generative models can learn high-dimensional expression distributions. Conditional WGAN-GP models have reproduced tissue and cancer properties in GTEx and TCGA [3]. More recently, Lacan and colleagues adapted denoising diffusion probabilistic and implicit models to bulk transcriptomics and reported strong gene-correlation, neighborhood, adversarial, and downstream classification metrics [4]. GeneJEPA instead learns masked-gene representations without reconstructing expression [5]. These approaches solve different problems: a generator can sample expression, whereas a representation learner needs an additional decoder or generative objective before it can do so.

Synthetic expression is commonly motivated as a remedy for small sample size. That framing is unsafe for biological inference. Multiple profiles sampled from one fitted model are not independent animals, and their apparent sample size cannot justify narrower confidence intervals or lower differential-expression P values. The more defensible question is whether a generator can expose a stable view of the training distribution that improves model selection or feature ranking, followed by validation using real, study-aware data.

We therefore separated three questions. First, can a paper-parity diffusion model learn broad mouse tissue expression from ARCHS4 and preserve tissue identity in held-out studies? Second, after OSDR adaptation, can it generate flight and ground-control profiles that pass distribution, diversity, memorization, and conditional-effect checks on a locked test? Third, can generated profiles improve flight-versus-ground feature discovery under a policy that does not count them as biological replicates?

The resulting evidence is deliberately tiered. OSD-457 thymus is an independently held-out test of a fixed generator and feature policy. A soleus analysis provides cross-accession real-data support but remains developmental because its domain contributed to generator adaptation. Lung transfer is mixed, and other tissues are reported as exploratory or negative. This structure distinguishes a publishable methods result with biological demonstrations from a broad claim that synthetic data discovered validated biology in every tissue.

![Study design and evidence ladder.](figures/figure_1_study_design.png)

<p class="caption"><strong>Figure 1. Training design and evidence ladder.</strong> (A) A paper-parity DDIM was pretrained on a healthy-preferred, tissue-balanced ARCHS4 mouse cohort and adapted to API-derived OSDR profiles with factorized tissue, flight status, study, and material conditioning. (B) Distribution validity, predictive utility, real-data biological support, and independent accession transfer were evaluated separately. Generated profiles were never treated as additional animals.</p>

## Materials and methods

### OSDR API cohort

The OSDR Biological Data API was used to discover assays, sample metadata, and processed expression [1]. Inclusion required organism *Mus musculus*, transcriptomic bulk RNA sequencing, a resolvable flight or ground-control label, and processed RSEM expected counts. Tissue and material labels were canonicalized from API metadata and audited aliases. All eligible OSDR data sources were considered. No raw integrated OSDR H5 file was used.

The API inventory contained 1,631 profile rows. Twenty-one technical replicate rows were aggregated to yield 1,610 biological profiles: 835 flight and 775 ground control, representing 75 accessions and 24 canonical material classes. The shared source matrix contained 48,694 mouse genes. Full-transcriptome transcripts per million (TPM) values were calculated with GENCODE M39 gene lengths before selecting the 974-gene mouse landmark panel.

The primary factorized model used 781 profiles for training, 536 for validation, and 293 for the locked test. Splits were made within accession, tissue, and condition, and every locked-test accession was represented in training. The locked test therefore measures within-study interpolation, not generalization to unseen studies. Study-level transfer was evaluated separately.

### ARCHS4 reference cohort

The local ARCHS4 mouse v2.5 HDF5 file contained 997,515 profiles and 53,511 genes. All profiles were audited, tissue labels were canonicalized, and tissue classes represented in eligible OSDR data were identified. Three reference cohorts were recorded: 23,614 control-only profiles from 3,213 GEO series, 62,299 healthy-preferred profiles from 5,307 series, and 134,250 broad profiles from 15,111 series.

The paper-parity DDIM used a deterministic balanced subset of 17,244 healthy-preferred profiles across 20 tissue classes. The split contained 9,796 training, 2,448 validation, and 5,000 held-out profiles, with GEO series assigned as whole units. Full-transcriptome TPM was computed before landmark selection. Training-set MaxAbs scaling was then applied to the 974-gene matrix. This matched the input dimensionality and preprocessing order of the reference diffusion implementation while adapting the landmark map from human to one-to-one or otherwise auditable mouse Ensembl genes.

**Table 1. Data scope.**

| Source | Available profiles | Analysis profiles | Split | Classes/accessions | Genes | Role |
|---|---:|---:|---|---:|---:|---|
| ARCHS4 mouse v2.5 | 997,515 | 17,244 | 9,796 / 2,448 / 5,000 | 20 tissues | 974 | Healthy-preferred tissue pretraining |
| NASA OSDR API | 1,631 rows | 1,610 biological profiles | 781 / 536 / 293 | 75 accessions | 974 | Conditional adaptation and FLT/GC analysis |

### Paper-parity diffusion pretraining

The DDIM implementation was taken from the official Lacan et al. source at commit `cde890154698fcea96c924804aaff04af3351b48` [4]. It follows the denoising-diffusion and implicit-sampling formulations of Ho et al. and Song et al. [13,14]. The model contained 227,109,786 parameters, two 8,192-unit residual hidden layers, 0.1 dropout, a scalar time embedding, a two-dimensional class embedding, and 1,000 quadratic diffusion steps from beta 0.0001 to 0.02. Training used the summed noise-prediction mean-squared error, antithetic timestep sampling, Adam with learning rate 0.0004783833151836702, OneCycle scheduling, automatic mixed precision, exponential moving average decay 0.999, and batch size 2,048. The model completed 15,000 epochs and 75,000 optimizer steps on an NVIDIA A100-SXM4-40GB.

### Factorized OSDR adaptation

The pretrained network was extended with embeddings for tissue, flight status, accession, and material type. A rank-512 low-rank domain adapter was fitted for 4,000 steps, followed by 1,000 condition-refinement steps. Sampling was balanced across tissue and accession during domain adaptation and across tissue, accession, and condition during refinement. Correlation regularization was applied to 256 genes at low diffusion timesteps. The adapted model retained the original DDIM backbone and 1,000-step noise schedule.

Generated expression was calibrated with train-only global and hierarchically shrunk accession/tissue means plus a positive missing-covariance residual. Flight status was not used to estimate calibrator means or covariance. Final accepted expression was clipped at zero before downstream biological use, and that policy was recorded. This calibration is part of the represented-study simulator and prevents the result from being interpreted as a generator for a new accession.

### Comparator models

We implemented the conditional WGAN-GP topology of Viñas et al. with a 64-dimensional noise vector, two 256-unit generator and critic layers, five critic updates per generator update, gradient-penalty weight 10, RMSProp learning rate 0.0005, and batch size 32 [3]. Twelve train-only calibration variants were evaluated. Because the strongest validation result remained distinguishable by an external adversary and failed accession-aware condition recovery, the WGAN locked test was not opened.

The exact released GeneJEPA architecture was screened with 4,096 train-selected highly variable genes and 43,744 ARCHS4 training exposures [5]. GeneJEPA has no expression decoder, so it was evaluated as a representation model rather than as a generator. Its held-out tissue embedding was compared with a classifier trained directly on processed expression.

### Generator evaluation

Distribution metrics followed the diffusion benchmark where transferable to mouse data [4]. They included gene-correlation-matrix agreement, manifold precision and recall, F1, nearest-neighbor adversarial accuracy, and Frechet distance (FD), using the distributional and manifold-evaluation principles introduced in references [15,16]. FD was computed in a train-fitted 50-dimensional PCA space and divided by the 95th percentile of real-versus-real split FD to account for finite-sample variability. Adversarial accuracy was required to fall between 0.40 and 0.60. Memorization was screened by comparing generated-to-training nearest-neighbor distances with the first percentile of training leave-one-out distances.

Four locked synthetic cohorts were sampled with prespecified seeds 5020-5023 against the same 293 real test profiles. Metrics were gated independently; no composite score was used. Gene-correlation agreement used the lower of an absolute 0.98 target and a real-data bootstrap fifth-percentile floor. Flight-minus-ground effects were compared in pooled data and, for skeletal muscle, after within-accession estimation.

### Downstream utility and feature guidance

Direct augmentation compared logistic classifiers trained on real profiles, generated profiles, or their union and evaluated all models on real held-out profiles. The feature-guidance workflow then separated feature discovery from classifier fitting. Five candidate arms were considered: real only; generated only; equal-weight real plus generated; a real-only classifier using generated-informed feature ranking; and a real classifier with generated profiles assigned 0.05 total training weight.

For development analyses, repeated outer splits were made within accession-by-condition strata. Inner splits selected the feature count, regularization, and ranking method. Balanced accuracy, AUROC, and average precision were compared independently. A generated-informed arm was eligible only when its mean was nonworse on all three metrics; repeat-level nonworse rates were retained because overlapping repeats are not independent.

Genes were considered selection-stable when they appeared in at least 50% of repeats and their classifier coefficient had at least 75% sign agreement. Stable generated-informed genes were not called biological findings until their effects were supported in real samples.

### Independent study confirmation

The fixed transfer policy selected one lung and one thymus accession for confirmation. OSD-900 lung contained 10 flight and 10 ground-control profiles; OSD-457 thymus contained 12 and 12. Neither accession had been used in earlier classifier tests. Both were removed from all OSDR generator-adaptation roles, the completed ARCHS4 checkpoint was reused, and the OSDR adaptation was rerun. OSD-464 lung and OSD-244 thymus served as fixed validation accessions. The feature policy and model-selection rule were declared before the two confirmation tests were opened.

Confirmation required balanced accuracy, AUROC, and average precision to be nonworse than the real-only baseline. Genotype was audited post hoc from sample names and GEO metadata because the initial API table did not expose it consistently. Each genotype stratum had equal flight and ground-control counts.

### Real-data effects, FDR, and pathway analysis

For each selected gene, flight-minus-ground effects were estimated separately within accession. Accession effects were combined with DerSimonian-Laird random-effects meta-analysis [17]. Benjamini-Hochberg false discovery rate (FDR) was applied within each declared gene or pathway family [6]. Leave-one-accession-out (LOO) sensitivity refitted the meta-analysis after removing each accession. A gene passed the strict stability rule only if its maximum LOO FDR remained below 0.05 and its effect direction did not reverse.

Reactome enrichment used the official current mouse Ensembl GMT generated from `ReactomePathways.txt` and `Ensembl2Reactome_All_Levels.txt` [7]. The tested background was the 974-gene panel. Enrichment P values were calculated by the hypergeometric test and adjusted within each tissue and selected-gene set. Reactome is hierarchical, so significant parent and child rows were interpreted as one process family rather than independent discoveries.

### Evidence tiers and expiMap triangulation

Independent confirmation required a test accession excluded from model adaptation and policy selection. Cross-accession development evidence required real random-effects and LOO support but did not establish unseen-study generator generalization. Exploratory evidence included predictive gains, ordinary FDR, or enrichment without a LOO-stable real gene set. Negative results lacked a coherent synthetic-guided gene or pathway result under the declared rules.

The generative results were compared qualitatively with a separate expiMap reference-query analysis of the same OSDR program [12]. Agreement was treated as triangulation across distinct model classes, not as independent replication when the same real samples contributed to both analyses.

## Results

### Broad ARCHS4 diffusion learned mouse tissue structure

The paper-parity ARCHS4 DDIM completed all 15,000 epochs in 5,987 seconds, with peak allocated GPU memory 5.93 GB. A classifier trained on generated profiles predicted held-out real tissue with balanced accuracy 0.869, compared with 0.895 when trained on real profiles (Fig. 2A). Gene mean and standard-deviation correlations were 0.997 and 0.944. Direct 974-gene precision and recall were 0.966 and 0.865; PCA-50 precision and recall were 0.986 and 0.943. Nearest-neighbor adversarial accuracy was 0.512 and PCA-50 FD was 0.0385.

Gene-correlation-matrix agreement was 0.879, below the stricter 0.98 target (Fig. 2B). The first two generated PCs had silhouette -0.271 even though held-out tissue classification remained strong, indicating that tissue information was not concentrated in the first two components. Raw model-scale samples contained 9.04% negative entries. These values are permissible in the unconstrained scaled space but are not physical TPM; downstream exports therefore require an explicit nonnegative policy.

The broad model was retained as a tissue-aware mouse reference generator but not claimed to reproduce every GTEx benchmark. Species, landmark mapping, cohort composition, split construction, and metric embedding differ from the human reference experiment.

### Adapted DDIM passed locked distribution gates but not exact-gene recovery

Across four locked generation seeds, mean gene-correlation agreement was 0.977 (range 0.974-0.979), precision 0.998, recall 0.997, F1 0.997, adversarial accuracy 0.458, and FD divided by real-split P95 was 0.075 (Fig. 2C,D; Table 2). All four generations passed finite-sample distribution, diversity, and memorization gates. Correlation remained slightly below the separate absolute 0.98 paper target.

Pooled flight-minus-ground effect recovery passed in three of four generations, with mean gene-effect correlation 0.598 and direction agreement 0.683. The skeletal-muscle accession-aware diagnostic passed in all four generations, with mean meta-effect correlation 0.608 across five accessions. Exact real-and-generated genes satisfying the strict LOO-FDR rule were sparse: zero, zero, zero, and one across the four seeds. The model therefore recovered the broad conditional contrast without reproducing a stable exact differential-gene list.

![Generator validation.](figures/figure_2_generator_validation.png)

<p class="caption"><strong>Figure 2. Generator validation.</strong> (A) Tissue balanced accuracy when a classifier was trained on held-out ARCHS4 real or synthetic profiles. (B) Broad-reference distribution metrics. The dashed line marks the strict correlation target. (C) Four OSDR locked-test generations; vertical marks show metric gates. (D) External adversarial accuracy and pooled or accession-aware flight-effect recovery. The shaded interval is the accepted adversarial-accuracy range. Full source values are in Tables S1-S2.</p>

**Table 2. Locked OSDR DDIM metrics over four prespecified generation seeds.**

| Metric | Mean | Range | Repeats passing |
|---|---:|---:|---:|
| Gene-correlation agreement | 0.977 | 0.974-0.979 | 4/4 finite-sample gate |
| Precision | 0.998 | 0.997-1.000 | 4/4 |
| Recall | 0.997 | 0.997-0.997 | 4/4 |
| F1 | 0.997 | 0.997-0.998 | 4/4 |
| Adversarial accuracy | 0.458 | 0.454-0.464 | 4/4 |
| FD / real-split P95 | 0.075 | 0.047-0.089 | 4/4 |
| Pooled FLT/GC effect recovery | r = 0.598 | 0.453-0.698 | 3/4 |
| Muscle accession-aware recovery | r = 0.608 | 0.498-0.664 | 4/4 |

### WGAN and GeneJEPA did not meet the generation objective

The strongest calibrated WGAN validation result achieved mean correlation 0.976, precision 0.976, recall 0.994, and F1 0.985, but external adversarial accuracy was 0.636, outside the prespecified 0.40-0.60 interval. Increasing residual variability moved adversarial accuracy toward chance while lowering correlation below its finite-sample floor. An earlier pooled flight-effect correlation of 0.805 disappeared after accession control (r = -0.022), and no evaluable tissue passed the accession-aware condition gate. The WGAN locked test was therefore left unopened.

The exact-architecture GeneJEPA duration screen processed 43,744 ARCHS4 exposures. Its held-out tissue embedding achieved balanced accuracy 0.703 and macro F1 0.701, below 0.839 and 0.840 from expression directly. Because GeneJEPA exposes no expression decoder, it could not generate flight or ground-control samples without adding and separately validating a decoder or guided generative model. It was not advanced for this objective.

### Direct augmentation did not improve flight classification

On the locked pooled OSDR test, real-only training achieved balanced accuracy 0.754 and AUROC 0.819. Synthetic-only training achieved 0.700 and 0.751, and equal real-plus-synthetic training achieved 0.734 and 0.801 (Fig. 3A). The generated data contained condition information, but naive augmentation diluted rather than improved the real classifier. Synthetic profiles were therefore evaluated as a low-weight or feature-ranking prior.

### Synthetic-guided features transferred strongly in thymus and inconsistently in lung

In independently held-out OSD-457 thymus, synthetic-guided feature ranking improved balanced accuracy from 0.500 to 0.833, AUROC from 0.840 to 0.972, and average precision from 0.876 to 0.976 (Fig. 3B). The final classifier used real rows only; generated profiles changed feature ranking but did not add training weight. In OSD-900 lung, low-weight guidance improved the same aggregate metrics from 0.400, 0.450, and 0.523 to 0.550, 0.550, and 0.635.

The post-hoc genotype audit separated each study into balanced strata (Fig. 3C). Thymus improved all three metrics in Nrf2-knockout mice and wild-type mice. Lung balanced accuracy and average precision improved in both strata, but knockout AUROC decreased from 0.560 to 0.520. Lung therefore met the pooled policy but did not provide uniformly validated subgroup performance.

Across the two studies, 11 profiles changed from incorrect to correct and none changed in the opposite direction. This profile-level comparison does not establish a study-level significance result because there was one independent accession per tissue and tissue selection followed an earlier transfer screen.

![Downstream utility.](figures/figure_3_downstream_utility.png)

<p class="caption"><strong>Figure 3. Downstream utility of generated expression.</strong> (A) Direct pooled augmentation on the locked real test. (B) Fixed synthetic-guided policies in independently held-out lung and thymus accessions. (C) Guided-minus-baseline metric changes after post-hoc genotype stratification. Thymus improved uniformly; lung knockout AUROC declined.</p>

### Held-out thymus effects support lower proliferative renewal

Selected-gene flight effects in OSD-457 correlated r = 0.975 between Nrf2-knockout and wild-type mice, and 86% had the same direction. The core genes *Birc5*, *Ccne2*, *Gmnn*, *Ube2c*, *Cdk1*, *Nusap1*, *Ccnb1*, and *Ccnb2* were lower in flight in both strata (Fig. 4A). Reactome rows for APC/C-CDC20-mediated cyclin degradation, G2/M checkpoints, DNA synthesis, and broader cell-cycle control passed FDR 0.05 (Fig. 4B).

This result is aligned with, but more specific than, prior thymus studies. STS-135 mouse thymus showed changes in cell-cycle and DNA-damage programs, including lower checkpoint-related expression [8]. A later ISS experiment reported marked thymus mass loss and partial artificial-gravity rescue of cell-cycle expression [9]. The current signature emphasizes mitotic completion and replication rather than acute apoptosis alone. Agreement between genotype strata suggests that the predictive signature is not confined to one Nrf2 background, but it does not prove Nrf2 independence.

Bulk thymus expression cannot distinguish lower transcription within proliferating thymocytes from loss or redistribution of proliferating cell populations. The defensible biological conclusion is lower abundance of a mitotic transcript program in flight, consistent with reduced proliferative renewal. Cell-resolved or histological confirmation is required to assign the effect to a cell-intrinsic mechanism.

![Thymus biology.](figures/figure_4_thymus_biology.png)

<p class="caption"><strong>Figure 4. Independently confirmed thymus biology.</strong> (A) Flight-minus-ground effects in real OSD-457 profiles for the eight core genes that were lower in both genotype strata. (B) Representative significant Reactome rows. Rows overlap hierarchically and describe one mitotic and DNA-replication process family rather than six independent pathways.</p>

### Anatomical separation exposes a soleus-specific metabolic program

Aggregate skeletal muscle concealed substantial anatomical heterogeneity. We therefore reused the fixed DDIM and three frozen synthetic views for five groups: extensor digitorum longus (EDL), gastrocnemius, quadriceps, soleus, and tibialis anterior. No neural network was retrained. Repeated nested development analysis included 24, 25, 35, 41, and 24 profiles, respectively.

All five groups selected a generated-informed arm, but only soleus combined nonworse predictive metrics with a multi-gene real-data LOO-stable set and coherent pathway enrichment (Fig. 5A). The soleus generated-only arm improved mean balanced accuracy, AUROC, and average precision by 0.025, 0.020, and 0.020 and was nonworse on all metrics in six of eight overlapping repeats.

Seven synthetic-selected genes passed real random-effects FDR 0.05, retained direction after every accession omission, and agreed with the generated effect: *Bdh1*, *Bnip3*, *Mef2c*, *Ech1*, *Pxmp2*, and *Gmnn* were lower in flight, while *Tpm1* was higher (Fig. 5B). *Arid5b* was an additional real-only LOO-stable gene. The selected sets were enriched for mitochondrial protein degradation, mitochondrial fatty-acid beta oxidation, and lipid metabolism (Fig. 5C).

The pattern links ketone or lipid utilization (*Bdh1*, *Ech1*), mitochondrial quality control (*Bnip3*), slow oxidative muscle identity (*Mef2c*), peroxisomal transport (*Pxmp2*), and contractile remodeling (*Tpm1*). Prior 30-day spaceflight profiling of mouse soleus reported a slow-to-fast shift and broad changes in oxidative metabolism, PPAR signaling, and contractile genes [10]. Unloading studies have also reported reduced soleus fatty-acid oxidation [11]. The pathway is therefore literature-aligned, while the compact gene prioritization and peroxisome-mitochondria emphasis are exploratory refinements rather than wholly de novo biology.

Soleus is not an independent generator transfer result. Its three accessions contributed to the development domain before nested feature analysis. LOO stability shows that the real gene effects are not driven by one accession, but it does not show that the DDIM generalizes to an unseen soleus study. A new accession excluded from adaptation and selection is required.

![Soleus biology.](figures/figure_5_soleus_biology.png)

<p class="caption"><strong>Figure 5. Skeletal-muscle group analysis.</strong> (A) Mean metric change for the selected generated-informed arm in five anatomical groups. (B) Seven synthetic-selected soleus genes that pass real-data random-effects FDR and the leave-one-accession-out stability rule. (C) Significant Reactome rows for the soleus core intersection. The rows are hierarchical and share genes.</p>

### Other muscle groups provide narrower hypotheses

Quadriceps retained one synthetic-selected, real LOO-stable gene, *Rbm6*, without a significant selected-set Reactome family. EDL showed flight-lower *Abcc5*, *Lsm6*, *Polr2i*, and *Tsc22d3* in both accessions, with small RNA-processing and nuclear-receptor enrichments; two accessions are insufficient for strict LOO confirmation. Tibialis anterior showed ordinary meta-analytic support for flight-higher *Cdkn1a*, *St3gal5*, *Cebpd*, *Pdhx*, and *Bnip3*, but its real-only classifier was already perfect, the selected arm tied it, and no pathway passed FDR. Gastrocnemius produced modest predictive gains without a real LOO-stable selected gene. These are reported as exploratory or negative.

### The remaining tissue screen is complementary, not a superset of expiMap

The broad developmental screen produced its largest mean classifier gains in spleen, thymus, retina, skin, and lung. These gains did not translate automatically into stable biology. Spleen selected one real LOO-stable gene, *Igfbp3*, without coherent Reactome enrichment. Skin, kidney, liver, retina, and lung had no selected gene passing the real LOO-FDR rule in this screen. Kidney showed a small porphyrin-related *Hmox1*/*Alas1* enrichment, but synthetic effect recovery failed and the genes were not LOO stable. Lung highlighted cell-cycle, senescence, and PI3K/AKT candidates, but no Reactome row passed FDR in the independent test and only 20%-28% of selected genes retained their training-study direction.

The evidence distribution therefore differs from the separate expiMap analysis. Both methods prioritize thymus. expiMap produced its broadest pathway evidence in thymus, skin, spleen, and kidney, whereas the generative workflow added its clearest complementary result in soleus. Skin, spleen, and kidney are not independently reproduced by synthetic guidance under the stricter real-gene rule. This is not a contradiction: expiMap tests constrained pathway scores, while synthetic guidance tests whether a learned expression distribution improves gene selection. Their null hypotheses, feature spaces, and regularization differ.

![Tissue evidence matrix.](figures/figure_6_tissue_evidence.png)

<p class="caption"><strong>Figure 6. Tissue evidence matrix.</strong> Thymus is the only tissue with a uniformly improved independent test and coherent held-out gene/pathway interpretation. Soleus has cross-accession developmental support but no unseen study. Lung is mixed; spleen, skin, and kidney remain exploratory; liver and retina do not provide coherent synthetic-guided biological results.</p>

**Table 3. Biological claim hierarchy.**

| Tissue | Evidence tier | Main signal | Defensible interpretation |
|---|---|---|---|
| Thymus | Independent held-out confirmation | Flight-lower mitotic genes; APC/C, G2/M, and DNA replication | Synthetic-guided ranking transferred to one unseen accession and both genotype strata |
| Soleus | Cross-accession development | Lower oxidative fuel handling and mitochondrial quality-control genes | Coherent real-data-stable hypothesis requiring unseen-accession confirmation |
| Lung | Mixed held-out exploratory | Predictive improvement; unstable directional genes | Potential multivariable predictor, not a validated pathway result |
| Spleen | Developmental exploratory | Large nested classifier gains; *Igfbp3* only | Prediction signal without a coherent stable pathway |
| Skin | Developmental exploratory | Cell-cycle/DNA-repair hypotheses | Does not independently reproduce the expiMap skin pathway result |
| Kidney | Developmental exploratory | Porphyrin candidates | Unstable and not independently reproduced |
| Liver, retina | Negative | No coherent retained gene/pathway set | No synthetic-guided biological claim |

## Discussion

### Synthetic generation is most useful here as a feature prior

The central finding is methodological. A conditional diffusion model can match broad OSDR expression distributions and recover an aggregate flight contrast, yet direct augmentation can still fail. Distribution fidelity, condition recovery, and downstream utility are distinct requirements. The accepted use was not to multiply the apparent number of animals. It was to provide a second, model-smoothed view of expression that changed feature ranking, after which a classifier and biological effect analysis remained anchored to real profiles.

This distinction explains why high precision, recall, and near-chance adversarial accuracy were necessary but insufficient. The locked DDIM reproduced neighborhoods and correlations but almost never reproduced the exact strict LOO-stable gene list. A generator can preserve a high-dimensional conditional distribution while uncertainty in small gene effects remains large. Conversely, the rejected WGAN showed high pooled flight-effect correlation but failed after accession control. Pooled agreement can be driven by study-condition structure rather than flight biology.

The independently held-out thymus result is the strongest evidence that synthetic guidance added information. OSD-457 was excluded from adaptation, the policy was fixed before testing, and performance improved in both genotype strata. The final classifier used only real profiles. This is narrower than claiming that generated thymus samples are biologically equivalent to new mice, but it is stronger and more reproducible than an in-sample augmentation gain.

### Thymus and soleus offer complementary biological stories

The thymus result refines established spaceflight immune biology. Prior work documents thymic involution and altered cell-cycle expression [8,9]. The current independently tested signature concentrates on cyclins, CDK1, UBE2C, BIRC5, NUSAP1, geminin, APC/C-mediated protein turnover, and G2/M control. Together these support lower proliferative renewal or a lower proportion of cycling thymocytes. Because the data are bulk, composition and cell-intrinsic regulation remain inseparable.

The soleus result addresses a different physiological axis. Weight-bearing slow muscle is especially sensitive to unloading. The selected genes describe lower oxidative substrate handling, altered mitochondrial turnover, reduced slow-muscle transcriptional identity, and contractile remodeling. This is compatible with known soleus atrophy and slow-to-fast transition [10,11], while nominating a compact set of genes for targeted validation. It is complementary to the expiMap result because synthetic guidance localized a broad, previously LOO-unstable fatty-acid-oxidation module to a smaller real-data-stable set.

### Why the other tissues remain useful

Negative and exploratory tissues constrain the method's scope. Lung demonstrates that predictive improvement can arise from covariance patterns even when univariate gene directions fail to replicate. Spleen shows that large within-study metric gains do not guarantee a coherent stable pathway. Skin and kidney show that a pathway model and a generative feature model can prioritize different aspects of the same tissue response. Retina and liver show that passing global generator metrics does not force every tissue contrast to become discoverable.

These results still provide actionable hypotheses. Lung cell-cycle and PI3K/AKT candidates can be tested in a newly held-out study with genotype modeled prospectively. Spleen *Igfbp3* can be examined alongside the expiMap immune-program result, but should not be presented as a synthetic discovery by itself. Skin and kidney should be revisited only after prespecifying a transfer study and cell-composition audit. The main manuscript includes these results to prevent a success-only narrative.

### Limitations

First, the OSDR locked test retains accessions on both sides of the split. It validates represented-study interpolation, not a new study. Study conditioning is necessary because study effects are substantial, but it also limits generation to observed accession profiles.

Second, the independent confirmation contains one lung and one thymus accession, and the tissues were selected after an earlier transfer screen. OSD-457 is a valid unseen-accession test of the frozen policy, but another prospectively selected thymus study is needed for study-level replication.

Third, the 974-gene landmark panel cannot discover effects outside that panel. It also compresses pathway coverage and may favor well-measured, broadly expressed genes. Full-transcriptome reconstruction was not used for biological claims.

Fourth, bulk tissue confounds cell abundance and cell-intrinsic regulation. Genotype was reconstructed post hoc for the confirmation studies because the first API-derived table lacked a reliable field. Future ingestion should capture genotype, sex, age, preservation method, and collection endpoint prospectively.

Fifth, the calibrator and factorized adapter use represented study domains. Calibrated within-domain fidelity correlations are internal diagnostics, not evidence of unseen-study generation. Soleus LOO analysis omits accessions from the real meta-analysis but does not retrain the generator for every omission.

Sixth, repeated nested splits overlap. Their nonworse rates describe stability, not independent P values. Reactome rows are hierarchical and cannot be counted as separate discoveries. The FDR rule controls the declared tested family, not all choices made during model development.

Finally, the broad ARCHS4 model does not reproduce every human GTEx metric, and its unconstrained output contains negative scaled values. The accepted OSDR export clips to nonnegative expression, which changes the generated distribution and must remain part of the model specification.

## Conclusions

A paper-parity DDIM pretrained on ARCHS4 and adapted with OSDR design covariates generated realistic represented-study mouse bulk expression. Synthetic cohorts did not improve flight classification by simple augmentation. Their useful role was narrower: guide feature selection, then require confirmation in real accession-aware data.

This workflow independently supports a flight-lower thymus mitotic program and develops a soleus hypothesis centered on oxidative metabolism, mitochondrial quality control, and contractile remodeling. Thymus is the strongest synthetic-data result and agrees with the strongest cross-method immune signal. Soleus is the main complementary result beyond expiMap. Lung, spleen, skin, kidney, liver, and retina define the method's exploratory and negative boundary. A manuscript is justified as a methods-plus-biological-demonstration study, provided these evidence tiers remain explicit and soleus is not described as independently validated.

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
