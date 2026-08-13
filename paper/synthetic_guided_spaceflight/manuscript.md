<div class="title-page">

<h1>A configurable generative transcriptomics framework identifies tissue-dependent synthetic utility in mouse spaceflight RNA-seq</h1>

<p class="subtitle">Benchmarking conditional WGAN and diffusion models across NASA OSDR mouse transcriptomes</p>

<p class="authors">Jason Trinh</p>

<p class="affiliation">Space Life Sciences Training Program, NASA Ames Research Center, Moffett Field, California, USA</p>

<p>Correspondence: jasontrinh@berkeley.edu</p>

<p class="draft-note"><strong>Research manuscript draft for mentor review.</strong> Confirm the author list and required NASA clearance before submission.</p>

</div>

## Abstract

**Background:** Mouse spaceflight studies provide access to tissues that cannot be sampled extensively from astronauts, but individual experiments are small and differ in design. We developed a configurable generative transcriptomics framework and asked whether synthetic expression could improve organ-specific flight analysis without treating generated profiles as new animals.

**Methods:** We assembled 1,610 mouse flight and ground-control bulk RNA-seq profiles through the NASA Open Science Data Repository API and audited all 997,515 profiles in ARCHS4 mouse. The reusable pipeline treated expression transformation, feature space, harmonization, training source, study scope, tissue structure, conditioning, and model family as configurable axes. Paper-based WGAN-GP and diffusion implementations were compared using fidelity, real-versus-synthetic separability, distributional-distance, and FLT/GC-effect metrics. The diffusion model used downstream was pretrained on 17,244 tissue-diverse ARCHS4 profiles. We compared real-only, synthetic-only, and real-plus-synthetic classifiers with the same 974 genes, preprocessing, splits, and regularization. Joint permutation and grouped SHAP tested correlated Reactome programs, and a secondary consensus analysis ranked compact gene panels. Flight associations were tested only in real profiles with accession-level random-effects models and Benjamini-Hochberg false-discovery control.

**Results:** None of nine matched liver harmonization arms provided adequate fidelity and conditional-effect recovery together. A calibrated study-conditioned WGAN-GP achieved correlation 0.976, precision 0.976, recall 0.994, F1 0.985, adversarial accuracy 0.636, and a Frechet-distance ratio of 0.144 on validation. Across four diffusion seeds, correlation was 0.974, precision 0.997, recall 0.996, F1 0.997, adversarial accuracy was 0.475, and the Frechet-distance ratio was 0.074. Diffusion was used downstream because it was less distinguishable from real data and had lower distributional distance while retaining high fidelity. Pooled augmentation reduced balanced accuracy from 0.754 to 0.737. In tissue-specific matched classifiers, real-plus-synthetic training was nonworse on pooled and accession-macro balanced accuracy, AUROC, and average precision in 18 of 27 analysis units; 16 improved at least one metric. Twenty-one real-data BH-FDR associations also had synthetic-supported importance: 15 in thymus, four in liver, and one each in skin and spleen. Seven thymus genes gained importance only after synthetic augmentation. The 15-gene thymus set enriched 26 Reactome terms, led by mitotic cell cycle (FDR 0.0047), and grouped importance retained six flight-lower thymus mitotic pathways. A secondary consensus analysis retained 49 associations, including coherent thymus and soleus panels; 11 associations were supported by both analyses.

**Conclusions:** Diffusion had the strongest joint validation metrics among the tested generators, but its downstream value depended on tissue and mode of use. Matched classifiers provided the direct test of synthetic contribution, while consensus ranking recovered correlated biological panels that can be diluted in an all-gene model. Thymus had the strongest support across both analyses. Liver, skin, and spleen produced narrower matched findings, and soleus produced a coherent secondary consensus result. Generated profiles did not add biological sample size, and all hypotheses require independent replication.

**Keywords:** spaceflight; bulk RNA-seq; synthetic data; diffusion model; skeletal muscle; thymus; NASA OSDR; ARCHS4

## Introduction

Spaceflight affects immune, musculoskeletal, metabolic, and barrier tissues through a combination of microgravity, radiation, confinement, altered nutrition, stress, and mission-specific procedures. Mouse flight experiments provide tissue access that is unavailable in astronauts, but their transcriptomic interpretation is difficult. Individual studies are small, missions differ in strain and duration, and condition labels can be entangled with study, material, genotype, or collection protocol. Pooling samples without preserving these design variables can convert study effects into apparent flight biology.

The NASA Open Science Data Repository (OSDR) exposes sample metadata and processed assay data through a queryable biological API [1]. We used it to assemble a cross-study cohort while retaining accession-level provenance. ARCHS4 uniformly processes a large collection of public human and mouse RNA-seq data [2] and supplied tissue-diverse reference profiles for models that would be underdetermined on OSDR alone.

Deep generative models can learn high-dimensional expression distributions. Conditional WGAN-GP models have reproduced tissue and cancer properties in GTEx and TCGA [3]. More recently, Lacan and colleagues adapted denoising diffusion probabilistic and implicit models to bulk transcriptomics and reported strong gene-correlation, neighborhood, adversarial, and downstream classification metrics [4]. We therefore compared these two generator families under a shared data and evaluation framework.

The model is only one part of the problem. Multi-study bulk RNA-seq can be represented as counts, CPM, TPM, or transformed and scaled expression; studies can be corrected, explicitly conditioned, modeled separately, or pooled. Published spaceflight workflows have used within-study standardization [23] and compared ComBat, ComBat-seq, and MBatch correction families [24]. MOBER offers a learned, inductive alternative based on an adversarial conditional variational autoencoder [25]. Any of these choices can improve one diagnostic while erasing flight-related structure or preserving study artifacts instead of biology.

We implemented WGAN-GP and DDIM in one configurable pipeline. Preprocessing, harmonization, training source, cohort structure, conditioning, and validation could be changed independently. We compared models using correlation, neighborhood, adversarial, distributional, diversity, memorization, and FLT/GC-effect metrics. The OSDR-adapted DDIM was less distinguishable from real profiles and had lower distributional distance while maintaining high fidelity, so we used it for downstream analysis.

Synthetic expression is commonly presented as a remedy for small sample size. Generated profiles, however, are not new biological replicates. After choosing diffusion for downstream analysis, we separated three questions: whether one pooled augmentation strategy helped at all, whether synthetic training changed tissue-specific classifiers when model choices were held fixed, and whether compact consensus panels offered additional biological interpretation.

Our primary biological question was whether tissue-specific synthetic training could add reproducible information about flight responses beyond real-only classifiers. Thymus had the strongest agreement between matched gene importance, real-data association, and pathway structure. Liver, skin, and spleen had narrower matched findings. A secondary consensus analysis recovered a coherent soleus metabolic panel that the all-gene marginal-importance test did not retain.

## Materials and methods

### Data sources

The OSDR Biological Data API was used to identify *Mus musculus* bulk RNA-seq assays with flight or ground-control labels [1]. Tissue and material names were harmonized while preserving study provenance. The resulting cohort contained 1,610 biological profiles from 75 accessions: 835 flight and 775 ground control. Full-transcriptome expression was converted to transcripts per million before selecting a 974-gene mouse landmark panel.

The local ARCHS4 mouse resource contained 997,515 public RNA-seq profiles [2]. All nine GEO series linked from the API-derived OSDR metadata were excluded before reference selection, removing 108 otherwise eligible profiles and preventing cross-resource sample overlap. A healthy-preferred, tissue-balanced subset of 17,244 profiles spanning 20 tissue classes was then used for model pretraining. Complete GEO series were assigned to one reference role, producing 10,150 training, 2,466 validation, and 4,628 test profiles with no series overlap. This backbone was used for both the held-out test and the broad development screens.

**Table 1. Data scope.**

| Source | Profiles used | Biological scope | Role |
|---|---:|---|---|
| ARCHS4 mouse v2.5 | 17,244 | 20 tissue classes | Mouse tissue pretraining |
| NASA OSDR API | 1,610 | 75 accessions; 835 flight and 775 ground control | Spaceflight analysis |

### Configurable generative benchmark

The framework treated expression representation, feature space, harmonization, training source, cohort structure, conditioning, model family, and validation design as separate axes (Table 2). A 463-row experiment planner combined these options and dispatched model-specific preprocessing and training through a common evaluation interface. It was not an exhaustive Cartesian benchmark; lower-cost screens identified configurations for full training and repeat evaluation.

**Table 2. Configurable pipeline and selected analysis branch.**

| Axis | Alternatives represented in the framework | Selected branch for downstream analysis |
|---|---|---|
| Expression | Raw, CPM, TPM; log1p/log2p1; z-score, robust, or MaxAbs scaling | Full-transcriptome TPM, then train-fitted MaxAbs scaling |
| Feature space | All shared genes, fold-selected HVGs, Reactome genes, mapped mouse L1000 | 974 mapped mouse L1000 landmarks |
| Harmonization | None, two study-wise z-score schemes, ComBat, ComBat-seq, three MBatch methods, MOBER | No global batch correction; accession represented as a condition |
| Training data | OSDR only, ARCHS4 only, ARCHS4 pretraining plus OSDR adaptation | ARCHS4 pretraining plus OSDR adaptation |
| Cohort structure | One or multiple studies; pooled or per-tissue fitting | All eligible OSDR accessions in a pooled tissue-conditioned generator |
| Conditions | FLT/GC, tissue, study, material, muscle group, and available design covariates | Tissue, FLT/GC, accession, and material |
| Model | WGAN-GP, DDIM | Factorized conditional DDIM |
| Validation | Grouped evaluation, repeat generation, unconditional controls | Four-seed 293-profile OSDR test |

### Preprocessing and harmonization

Each generator received its paper-based preprocessing and a shared benchmark representation. The WGAN-GP branch used log-transformed expression with training-gene z-scores [3]. The diffusion branch used TPM, a mapped mouse landmark panel, and training-fitted MaxAbs scaling [4]. All feature selection and transform statistics were fitted on training partitions.

Nine harmonization arms were compared in a matched liver experiment: no correction, two study-wise standardization schemes, ComBat, ComBat-seq, three MBatch methods, and MOBER [23-25]. Advancement required preservation of expression fidelity and FLT/GC effects, not simply reduced study separation. Methods without a natural frozen transform for a new batch were treated as transductive sensitivity analyses. Full transform definitions and results are provided in Supplementary Methods S7 and Supplementary Tables S13-S14.

### Model training and selection

We implemented paper-based WGAN-GP and DDIM generators [3-5,12]. The DDIM was first trained on the 17,244-profile ARCHS4 reference, then adapted to OSDR with tissue, FLT/GC, accession, and material conditions. Exact architectures, optimization schedules, adaptation stages, seeds, and hardware records are provided in Supplementary Methods S4-S7.

Complete GEO series and OSDR accessions were grouped when study-level separation was required. Candidate generators were evaluated for correlation structure, neighborhood precision and recall, external real-versus-synthetic separability, distributional distance, diversity, memorization, and recovery of FLT/GC effects [13,14]. Metrics were assessed independently and each model's evaluation split was reported. Exact thresholds are provided in Supplementary Methods S6. The factorized DDIM had lower adversarial accuracy and distributional distance than WGAN-GP while retaining high correlation and neighborhood fidelity. It generated FLT or GC expression for represented tissue and study contexts, and those profiles were used in the downstream screens.

### Evaluation funnel and synthetic-guided analysis

We first tested pooled augmentation across all tissues. Individual profiles were combined in one FLT/GC classification problem rather than averaged by tissue. The locked benchmark compared classifiers trained with real profiles, generated profiles, or both.

The primary tissue-specific analysis then held the classifier design fixed. We evaluated 22 canonical tissues and five anatomical muscle groups with ridge logistic classifiers that all used the same 974 genes. Eight outer splits were used per analysis unit. The scaler and regularization value were fitted or selected from real outer-training data and then reused without change for real-only, synthetic-only, and real-plus-synthetic arms. In the combined arm, the full synthetic set received the same total weight as the real set. Every arm was scored on the same held-out real profiles. Performance was summarized both across pooled test profiles and as a mean of accession-level scores. A synthetic arm passed the joint utility gate only when mean balanced accuracy, AUROC, and average precision were all no worse than real only in both summaries.

Accession-blocked permutation importance measured the held-out-real AUROC loss after shuffling one gene, with ten shuffles per gene and fit. Exact linear SHAP provided the direction of each gene's contribution. Candidate genes had to pass real-data BH FDR, come from a synthetic arm that passed the joint utility gate, show a mean permutation loss of at least 0.001 with a positive loss in at least half of outer splits, and have coefficient and SHAP directions consistent with the observed FLT/GC effect. This matched all-gene analysis is the primary test of synthetic contribution.

We also tested Reactome pathways as joint feature groups in the same classifiers. All member genes were permuted together within accession, and grouped SHAP was calculated as the sum of their linear SHAP values. Real-data pathway scores were calculated independently by standardizing genes within accession and averaging the genes in each pathway before random-effects analysis. A grouped result required joint classifier utility, repeat-consistent pathway permutation loss, SHAP direction consistent with the observed effect, and real-data pathway BH FDR below 0.05. Overlapping Reactome terms were retained in the complete table and marked by a separate Jaccard-based nonredundancy filter.

Marginal importance can understate a correlated biological program. Ridge regression can divide predictive weight among genes with similar expression patterns. If gene A is shuffled while correlated gene B remains intact, B can carry much of the same information and the measured loss for A will be small. We therefore used a secondary consensus analysis to recover compact panels. It compared five tissue-specific uses: real-only, generated-only, equal-weight real plus generated, consensus-guided ranking followed by real-only fitting, and consensus-guided ranking with generated profiles at 0.05 total weight. Genes stable with real-only and synthetic-informed ranking were called reinforced; genes stable only with an eligible synthetic-informed arm were called synthetic-promoted. These labels describe repeated selection, not biological novelty.

Flight-minus-ground effects were estimated within each OSDR study and then summarized with a random-effects model [15]. This kept mission-level contrasts separate before meta-analysis and prevented generated profiles from increasing the biological sample count.

Within each tissue, the 974 real-data gene-level meta-analysis P values were adjusted with the Benjamini-Hochberg procedure, and BH FDR below 0.05 defined the statistical inclusion rule [6]. Generated profiles were never entered as biological replicates. Accession-direction agreement, between-study heterogeneity, and leave-one-accession-out results were retained as interpretation and sensitivity measures rather than inclusion requirements. The complete BH-FDR inventory is provided in Supplementary Table S11. Matched classifier utility, candidates, Reactome results, and the matched-consensus crosswalk are in Tables S18-S21. Literature annotations for matched genes and grouped pathways are in Tables S22-S24. Reactome was used to group selected genes into biological processes [7].

We performed targeted literature review for the 49 associations from the secondary consensus analysis, the 21 associations from the primary matched analysis, and the ten eligible grouped pathways. Eleven matched genes already had a tissue-specific consensus annotation, which was reused without changing its label; the ten matched-only genes and all grouped pathways were reviewed separately. Selection behavior and literature interpretation were recorded independently. Prior evidence was classified as aligning, complementary, ambiguous, or unmatched. The annotation pipeline fixed the tissue, feature, and observed direction before searching spaceflight or microgravity literature and relevant process terms. Evidence scope distinguished direct and process-level agreement from mechanistic context, and source relationship recorded whether a study was independent, reused public OSDR cohorts, or was used only for context. Unmatched means that the targeted search found no sufficiently specific tissue or process match; it is not proof of novelty. Supplementary Tables S16-S17 and S22-S24 contain all decisions and sources.

**Table 3. Evaluation stages and permitted interpretation.**

| Stage | Analysis | Data separation | Question answered | Permitted interpretation |
|---|---|---|---|---|
| 1 | Pooled utility benchmark | Locked profiles from represented studies | Does one augmentation policy help across tissues? | Method-level result; no biological claim |
| 2 | Matched all-gene classifiers | Same 974 genes, splits, preprocessing, and regularization; held-out real profiles | Does the training source change predictive utility and gene importance? | Primary evidence of synthetic contribution within represented studies |
| 3 | Consensus panel analysis | Five tissue-specific arms with repeated ranking | Which compact gene panels remain stable across real and generated views? | Secondary module and pathway interpretation |
| 4 | Real-data random-effects BH FDR | FLT-GC effects estimated separately within each accession | Are nominated genes associated with flight in observed studies? | Real biological association; synthetic status reported separately |

## Results

### Generator metrics favored the OSDR-adapted diffusion model

The broad ARCHS4 DDIM retained tissue identity and passed most distributional tests on 4,628 held-out profiles, but its gene-correlation agreement missed the prespecified floor. It was therefore retained only as tissue-conditioned initialization.

No liver harmonization arm passed the joint fidelity and conditional-effect criteria (Supplementary Tables S13-S14). Some methods improved a single diagnostic while degrading neighborhood fidelity, real-versus-synthetic separability, or distributional distance. The selected path therefore retained TPM/MaxAbs expression without global correction and represented accession explicitly during OSDR adaptation.

WGAN-GP achieved correlation 0.976, F1 0.985, adversarial accuracy 0.636, and FD/real-P95 0.144. DDIM achieved correlation 0.974, F1 0.997, adversarial accuracy 0.475, and FD/real-P95 0.074. DDIM was therefore used downstream because it was less distinguishable from real data and had lower distributional distance while retaining high fidelity (Table 4).

![Generator metrics and model choice.](figures/figure_1_generator_validation.png)

<p class="caption"><strong>Figure 1. Generator metrics and model choice.</strong> (A) Real-trained and synthetic-trained classifiers retained the same broad ARCHS4 tissue balanced accuracy. (B) WGAN-GP validation metrics and DDIM test metrics on their stated evaluation splits; adversarial accuracy closer to 0.5 indicates lower real-versus-synthetic separability. (C) Fraction of four DDIM generation seeds passing each general fidelity or FLT/GC-effect criterion. DDIM was used downstream because it combined near-chance adversarial accuracy, lower distributional distance, and high fidelity.</p>

**Table 4. Generator metrics and model choice. Values are reported on each model's stated evaluation split and are not paired on one common split.**

| Model | Evaluation split | Corr. | Precision | Recall | F1 | AA | FD/real P95 | FLT/GC recovery | Decision |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| Broad-reference DDIM | 4,628 held-out ARCHS4 profiles | 0.878 | 0.951 | 0.890 | 0.919 | 0.515 | 0.866 | NA | Initialization only |
| Study-conditioned WGAN-GP | 536-profile validation; 6 seeds | 0.976 | 0.976 | 0.994 | 0.985 | 0.636 | 0.144 | 6/6 | Not used downstream |
| Factorized DDIM | 293-profile OSDR test; 4 seeds | 0.974 | 0.997 | 0.996 | 0.997 | 0.475 | 0.074 | 3/4 | Used downstream |

### Diffusion generation recovered tissue-conditioned expression structure

The reverse trajectory provides a direct view of the model transforming noise into tissue-conditioned expression. In the ARCHS4 reference model, profiles moved from an overlapping noise cloud at timestep 1,000 toward the tissue-structured real-data manifold at timestep 0 (Fig. 2, top). After OSDR adaptation, locked real and generated profiles occupied similar tissue-defined regions, and flight and ground-control profiles remained interspersed within that broader structure (Fig. 2, bottom). These PCA projections are descriptive: model selection used the full correlation, precision, recall, adversarial, Frechet-distance, and conditional-effect gates rather than visual similarity.

### Pooled augmentation did not improve the multi-tissue classifier

We first asked whether synthetic profiles could simply expand one multi-tissue FLT/GC training cohort. On the OSDR test, balanced accuracy was 0.754 with real-only training, 0.695 with generated-only training, and 0.737 with real-plus-generated training. The corresponding AUROCs were 0.820, 0.751, and 0.791. Pooled augmentation therefore did not improve the broad classifier.

Tissue identity dominates bulk expression, while flight effects differ among organs. One decision rule can therefore dilute tissue-specific condition signals even when the generator is conditioned on tissue. This negative pooled result motivated separate analysis within each tissue.

<div class="figure-block">
  <div class="figure-composite">
    <img class="trajectory-panel" src="figures/figure_2a_archs4_denoising_trajectory.png" alt="ARCHS4 DDIM denoising trajectory across mouse tissues">
    <img src="figures/figure_2b_locked_real_vs_synthetic_pca.png" alt="Locked OSDR real and synthetic profiles in PCA space">
  </div>
  <p class="caption"><strong>Figure 2. Diffusion generation across reference pretraining and OSDR adaptation.</strong> Top: ARCHS4 tissue-conditioned profiles at DDIM timesteps 1,000, 200, and 0 in a PCA space fitted to real reference expression; gray points are real ARCHS4 profiles and colors identify generated tissue conditions. Bottom: OSDR test profiles for generation seed 5020; circles denote real profiles and crosses denote generated profiles, colored by tissue on the left and flight condition on the right. PCA views are descriptive and do not replace quantitative validation metrics.</p>
</div>

### Matched all-gene classifiers identified direct synthetic contribution

The matched analysis fitted 648 classifiers across 27 analysis units. Real-plus-synthetic training was no worse than real-only training on all three pooled and all three accession-macro metrics in 18 units. Sixteen of those units improved at least one metric; bone marrow and brown adipose tissue were exact ties. Synthetic-only training passed the same gate in six units. These comparisons used the same 974 genes, splits, scaler, regularization, and held-out real profiles, so the differences can be attributed to training source rather than a different feature set or classifier setting.

Real-plus-synthetic training produced clear gains in several tissues (Table 5). Eye, retina, lung, skin, thymus, spleen, liver, and pooled skeletal muscle improved pooled balanced accuracy by 0.010 to 0.100. Their accession-macro metrics were also nonworse. Synthetic-only transfer was less reliable. For example, its pooled balanced accuracy in skeletal muscle was 0.670, compared with 0.952 for real-only training, despite high AUROC.

**Table 5. Matched real-plus-synthetic classifier results on held-out real profiles. All displayed units were nonworse on pooled and accession-macro balanced accuracy, AUROC, and average precision.**

| Tissue | Balanced-accuracy change | AUROC change | Average-precision change |
|---|---:|---:|---:|
| Eye | +0.094 | +0.156 | +0.115 |
| Retina | +0.100 | +0.075 | +0.063 |
| Lung | +0.062 | +0.086 | +0.108 |
| Skin | +0.074 | +0.087 | +0.076 |
| Thymus | +0.061 | +0.046 | +0.042 |
| Spleen | +0.057 | +0.075 | +0.072 |
| Liver | +0.043 | +0.037 | +0.029 |
| Skeletal muscle, pooled | +0.010 | +0.011 | +0.014 |

Twenty-one unique tissue-gene associations passed the complete matched gate. All BH-FDR statistics came from real profiles. Thymus accounted for 15 associations, liver for four, and skin and spleen for one each. Seven thymus genes, *Klhdc2*, *Snx7*, *Etv1*, *Plscr1*, *Tspan3*, *Socs2*, and *Kif20a*, crossed the marginal-importance threshold after real-plus-synthetic training but not in the real-only classifier. The other matched associations had measurable importance in both arms. No anatomical muscle group passed the complete matched utility and gene-importance gate.

Literature review of these 21 matched associations classified nine as aligning, nine as complementary, one as ambiguous, and two as unmatched. The aligning set was concentrated in the flight-lower thymus cell-cycle response. The four liver genes connected to hepatic transcription, metabolic signaling, proteostasis, or antigen presentation, but only *Gtf2a2* aligned with a previously reported same-tissue process. *Klhdc2* and *Snx7* had no sufficiently specific thymus-spaceflight match in the targeted search.

Joint pathway importance retained ten unique tissue-pathway associations, nine after removing one highly overlapping thymus term for display. Six flight-lower thymus pathways covered APC/C-Cdc20 control, chromosome condensation, and the G2/M checkpoint; all six aligned with prior reports of reduced thymic proliferation. Two flight-higher skin terms described RIPK1-regulated necroptosis and were classified as complementary because the mechanism is established in skin but has not been directly shown in spaceflight. Flight-higher spleen AP-1 activation and thymus ERBB2-PTK6 signaling were classified as ambiguous because prior direction varied with immune, mission, or growth context. The pathway FDR values came only from real OSDR profiles; grouped permutation and SHAP described the fitted classifiers.

### Consensus ranking provided secondary panel-level interpretation

The consensus analysis asked a different question: which genes repeatedly enter a compact predictive panel after comparing real and generated rankings? It retained 49 real-data BH-FDR associations, with 26 synthetic-promoted and 23 reinforced selections. Eleven also passed the matched all-gene analysis, 38 were consensus-only, and ten were matched-only.

This difference follows from correlation among genes. An all-gene classifier can spread weight across several genes in the same program. When one is permuted, correlated partners can preserve prediction, which reduces that gene's marginal importance. Consensus ranking can retain several members of the same program because it ranks their repeated positions before fitting a compact classifier. The consensus results are therefore useful for pathway interpretation, but the matched results provide the more direct test of synthetic contribution.

The consensus screen covered 22 canonical tissues and five muscle groups. Its 49 associations occurred in adrenal gland, eye, kidney, pooled skeletal muscle, skin, spleen, thymus, gastrocnemius, soleus, and tibialis anterior. Literature review classified 22 as aligning, 19 as complementary, four as ambiguous, and four as unmatched. Selection status and literature status were independent. Full evidence scopes and citations are provided in Supplementary Tables S16-S17.

### Convergent thymus analyses identify lower proliferative renewal

Thymus had the strongest matched gene-level result. Real-plus-synthetic training improved balanced accuracy by 0.061, AUROC by 0.046, and average precision by 0.042 on held-out real profiles while passing every pooled and accession-macro metric gate. Fifteen real-data BH-FDR associations also had synthetic-supported marginal importance. Six flight-higher genes were retained (*Klhdc2*, *Snx7*, *Etv1*, *Plscr1*, *Tspan3*, and *Socs2*) together with nine flight-lower genes (*Nusap1*, *Stmn1*, *Birc5*, *Ccnb2*, *E2f2*, *Ube2c*, *Cdc20*, *Gmnn*, and *Kif20a*). Seven of these genes gained measurable marginal importance only after synthetic augmentation.

Reactome enrichment of the 15 matched genes found 26 significant terms. The leading term was mitotic cell cycle (FDR 0.0047), supported by *Ube2c*, *Kif20a*, *Cdc20*, *Gmnn*, *E2f2*, and *Ccnb2*. A separate grouped-importance test also retained six flight-lower mitotic pathways, which showed that their genes contributed collectively even when correlated genes could reduce single-gene permutation loss. The secondary consensus analysis retained 16 thymus associations. Nine overlapped the matched result, while the compact panel added correlated cell-cycle genes such as *Cdk1*, *Top2a*, *Aurka*, *Ccne2*, *Pcna*, and *Ccnf*. Figure 3 shows the consensus cell-cycle panel and its real accession-level effects.

The two methods therefore support different parts of the same interpretation. Matched importance shows that synthetic augmentation changed held-out-real prediction and retained individual genes. Consensus ranking assembles a broader mitotic program whose members can substitute for one another in an all-gene classifier. All association statistics still come from observed OSDR profiles.

Prior thymus studies report related biology, although the present result is more specific. STS-135 mouse thymus showed changes in cell-cycle and DNA-damage programs, including lower checkpoint-related expression [8]. A later ISS experiment reported marked thymus mass loss and partial artificial-gravity rescue of cell-cycle expression [9]. The current signature emphasizes mitotic completion and replication rather than acute apoptosis alone. Agreement across accessions suggests a shared thymus response despite study heterogeneity, but it does not identify the responsible cell population.

Two flight-higher genes extend the cell-cycle interpretation. The consensus panel included *Hsd17b11*, which localizes to lipid droplets and can support regulated lipolysis [28]. Both analyses retained *Etv1*, a regulator of CD4+ T-cell activation and proliferation [29]. In a shrinking thymus with altered cell populations, these signals could reflect a shift in lipid handling, T-cell state, or relative cell abundance. They do not establish either gene as a driver of thymic involution.

Bulk thymus expression cannot distinguish lower transcription within proliferating thymocytes from loss or redistribution of proliferating cell populations. The defensible biological conclusion is lower abundance of a mitotic transcript program in flight, consistent with reduced proliferative renewal. Cell-resolved or histological confirmation is required to assign the effect to a cell-intrinsic mechanism.

![Thymus biology.](figures/figure_3_thymus_biology.png)

<p class="caption"><strong>Figure 3. Secondary consensus view of the thymus response.</strong> (A) Real-data random-effects flight-minus-ground estimates for ten genes across five thymus accessions; coral denotes genes promoted by consensus ranking and teal denotes genes retained by both real-only and consensus ranking. (B) The compact panel converges on mitotic and DNA-replication processes. The primary matched all-gene result is summarized in Figure 5 and Supplementary Table S19.</p>

### Consensus ranking identifies a soleus metabolic program

Aggregate skeletal muscle concealed substantial anatomical heterogeneity. We therefore examined extensor digitorum longus, gastrocnemius, quadriceps, soleus, and tibialis anterior separately. Soleus produced the clearest biological pattern: its selected genes showed consistent flight effects across three accessions and converged on related metabolic processes.

The secondary consensus screen selected real-plus-generated training. Balanced accuracy increased from 0.925 to 0.963, AUROC remained 0.980, and average precision increased from 0.980 to 0.986. Five BH-FDR genes were stable in both real-only and consensus-guided selection: *Bdh1*, *Ech1*, *Bnip3*, and *Decr1* were lower in flight in all three accessions, while *Tpm1* was higher (Fig. 4B). Four genes passed leave-one-accession-out FDR; *Decr1* did not.

The consensus method did not promote a soleus gene absent from stable real-only selection. Its contribution was reinforcement of a coherent existing pattern. Reactome connected the retained genes to mitochondrial protein turnover, mitochondrial fatty-acid beta-oxidation, and lipid metabolism (Fig. 4C). *Bdh1* and *Ech1* support oxidative substrate handling, *Bnip3* supports mitochondrial quality control, *Decr1* contributes to unsaturated fatty-acid oxidation, and *Tpm1* suggests contractile remodeling.

Soleus did not pass the primary matched all-gene gate. Real-plus-synthetic training increased balanced accuracy by 0.013 but reduced AUROC by 0.015 and average precision by 0.009, and no soleus gene passed the full marginal-importance rule. The soleus result is therefore a secondary panel-level finding, not direct evidence that synthetic augmentation improved the fixed all-gene classifier.

Prior 30-day spaceflight profiling of mouse soleus reported a slow-to-fast shift and broad changes in oxidative metabolism, PPAR signaling, and contractile genes [10]. Unloading studies have also reported reduced soleus fatty-acid oxidation [11]. This is literature-supported panel reinforcement rather than de novo gene discovery.

Unlike thymus, soleus was represented during model development. Its cross-study consistency makes it a focused biological hypothesis, but an entirely unseen soleus study is still needed for independent confirmation.

![Soleus biology.](figures/figure_4_soleus_biology.png)

<p class="caption"><strong>Figure 4. Secondary consensus analysis of skeletal muscle.</strong> (A) Number of consensus-selected genes that also pass real leave-one-accession-out FDR in each anatomical muscle group. (B) Five reinforced soleus BH-FDR genes with consistent real flight effects. (C) Their strongest shared biological processes center on mitochondrial turnover and lipid metabolism.</p>

### Additional muscle findings were limited to the consensus analysis

The pooled skeletal-muscle consensus screen improved balanced accuracy, AUROC, and average precision by 0.071, 0.036, and 0.037. Twelve BH-FDR genes entered the compact panel, nine of which passed leave-one-accession-out sensitivity. The associated sets included interferon signaling and sialic-acid metabolism. In the matched all-gene analysis, real-plus-synthetic training improved all six pooled and accession-macro metrics, but no individual muscle gene passed the complete importance gate. The broad classifier result supports utility at the tissue level, while the gene panel remains secondary.

Among the remaining anatomical groups, the consensus analysis retained *Cdkn1a*, *St3gal5*, *Bnip3*, and *Cebpd* in tibialis anterior and *Fhl2* and *Nfkbia* in gastrocnemius. None passed leave-one-accession-out FDR or the matched all-gene gate. EDL and quadriceps produced no synthetic-supported gene result.

### Spleen separates matched and panel-level findings

In the matched all-gene analysis, real-plus-synthetic training improved spleen balanced accuracy by 0.057, AUROC by 0.075, and average precision by 0.072. Flight-higher *Loxl1* was the only gene to pass the complete matched gate. The secondary consensus panel also retained *Loxl1* and added flight-higher *Rai14*, *Ptprk*, and *Myl9*. *Rai14*, *Ptprk*, and *Loxl1* had the same direction in all six studies; *Myl9* agreed in five of six. The four-gene panel did not produce significant Reactome enrichment.

The consensus genes provide a tentative structural hypothesis rather than a pathway claim. RAI14 links actin-associated mechanosensing to Hippo signaling, PTPRK regulates cell-cell junctions, MYL9 contributes to actomyosin contractility, and LOXL1 supports elastic extracellular-matrix maintenance [20-22]. Their shared direction is compatible with altered adhesion, mechanics, or tissue architecture, but the present bulk data cannot identify a common source cell or establish one mechanism.

### Other organ responses were heterogeneous

Lung gained 0.062 balanced accuracy, 0.086 AUROC, and 0.108 average precision with matched real-plus-synthetic training, but no gene passed the real-data BH-FDR screen. This is predictive utility without a supported gene-level biological claim.

Skin gained 0.074 balanced accuracy, 0.087 AUROC, and 0.076 average precision with matched real-plus-synthetic training. Flight-higher *Plscr1* passed both matched and consensus selection and had the same direction in all six studies, although it did not pass leave-one-accession-out FDR. Broader consensus sets were enriched for G1/S and DNA-repair processes, consistent with prior OSDR skin analyses reporting strong mission, strain, recovery-interval, and anatomical-site dependence [16].

Liver gained 0.043 balanced accuracy, 0.037 AUROC, and 0.029 average precision with matched real-plus-synthetic training. Four flight-lower genes passed the complete matched gate: *Grb10*, *Ppic*, *H2-DMa*, and *Gtf2a2*. All were also important in the real-only classifier, so this is shared rather than newly promoted importance. The four genes did not produce significant Reactome enrichment and are interpreted as individual candidates rather than a pathway.

Kidney contributed two consensus-only genes. *Slc37a4* was higher in flight in all six kidney accessions, while *Inpp4b* was higher in the random-effects model with one negative accession effect. The pair did not yield significant Reactome enrichment and neither passed the matched all-gene gate. They remain secondary renal metabolic-signaling hypotheses [17-19].

The consensus analysis produced flight-lower *Psmb8* and *Tspan4* in adrenal gland, each with the same direction in all three studies. Neither passed the matched all-gene gate. PSMB8 is an interferon-inducible immunoproteasome subunit linked to inflammation and adipocyte homeostasis [30], making lower adrenal *Psmb8* a plausible immune, proteostasis, or tissue-composition hypothesis. Retina had strong matched classifier gains but no synthetic-supported BH-FDR gene. Eye also had large matched utility gains without a retained gene-level candidate.

![Tissue evidence hierarchy.](figures/figure_5_tissue_evidence.png)

<p class="caption"><strong>Figure 5. Matched classifier evidence and its relation to consensus ranking.</strong> (A) Changes in held-out-real performance for matched real-plus-synthetic classifiers. (B) Real-data BH-FDR genes with synthetic-supported marginal importance. (C) Eleven tissue-gene associations were retained by both the matched and consensus analyses. (D) Correlated genes can substitute for one another in an all-gene classifier, reducing the measured permutation loss for either gene; consensus ranking can retain the broader panel.</p>

**Table 6. Tissue interpretation across the primary matched and secondary consensus analyses. Complete real-data BH-FDR results are in Supplementary Tables S11-S12.**

| Tissue | Matched all-gene result | Secondary consensus result | Interpretation |
|---|---|---|---|
| Thymus | 15 genes; seven promoted; 26 significant Reactome terms | 16-gene panel; nine genes overlap matched result | Strongest joint result; lower mitotic renewal in flight |
| Liver | Four shared-importance FLT-lower genes | No retained panel | Matched gene-level candidates without pathway enrichment |
| Skin | FLT-higher *Plscr1* | *Plscr1* plus broader selected pathways | Gene-level result supported by both analyses |
| Spleen | FLT-higher *Loxl1* | *Loxl1*, *Rai14*, *Ptprk*, and *Myl9* | Matched anchor within a tentative adhesion and cytoskeletal panel |
| Skeletal muscle, pooled | Classifier utility improved; no gene passed the full gate | 12-gene interferon and sialic-acid panel | Predictive utility with secondary panel-level interpretation |
| Soleus | Joint utility and gene gate not passed | Reinforced five-gene mitochondrial and lipid panel | Coherent secondary consensus hypothesis |
| Eye, retina, lung | Classifier utility improved; no retained BH-FDR gene | No principal gene-level claim | Predictive result without a supported biological candidate |
| Kidney, adrenal, gastrocnemius, tibialis anterior | No retained matched gene | Narrow consensus-only candidates | Exploratory panel results |

## Discussion

### Generator metrics favored diffusion for downstream analysis

Both WGAN-GP and DDIM reproduced expression structure, but DDIM was less distinguishable from real profiles and had lower distributional distance. OSDR adaptation was needed because the broad ARCHS4 model retained tissue identity but missed the correlation target. Harmonization also showed that reducing study structure could damage expression neighborhoods or FLT/GC effects, so accession was represented explicitly. We therefore used DDIM for downstream analysis. This comparison applies to the present benchmark and does not establish that diffusion is universally superior to adversarial training.

### What synthetic data added

Pooled augmentation did not improve FLT/GC classification, but matched tissue-specific augmentation passed all six pooled and accession-macro metric checks in 18 of 27 analysis units. The fixed all-gene comparison is the clearest evidence that generated profiles changed prediction on held-out real samples. It retained 21 real-data BH-FDR associations in four tissues, with most concentrated in thymus.

The consensus analysis added a second view. It selected compact panels after combining real and generated rankings and retained 49 real-data associations. Only 11 overlapped the matched result. This is not a contradiction. Correlated genes share predictive information. Ridge can divide their coefficients, and permuting one gene leaves its partners available, so each gene can have a small individual loss even when the program matters. Consensus ranking is better suited to describing such panels, while the matched all-gene analysis is better suited to claims about synthetic contribution.

Generated profiles did not add biological sample size. All associations were tested in observed profiles after estimating FLT/GC effects separately by accession. Literature annotations for the consensus, matched-gene, and grouped-pathway results separate prior alignment from model-selection behavior. These labels organize interpretation; they do not establish novelty.

### Thymus is the strongest joint result; soleus is a secondary panel

The thymus result refines established spaceflight immune biology. Prior work documents thymic involution and altered cell-cycle expression [8,9]. The current signature concentrates on cyclins, CDK1, UBE2C, BIRC5, NUSAP1, geminin, APC/C-mediated protein turnover, and G2/M control. The combined pattern is consistent with lower proliferative renewal or a lower proportion of cycling thymocytes. Because the data are bulk, composition and cell-intrinsic regulation remain inseparable.

The soleus result addresses a different physiological axis. Weight-bearing slow muscle is especially sensitive to unloading. Lower *Bdh1*, *Ech1*, *Bnip3*, and *Decr1* together with higher *Tpm1* describe reduced oxidative substrate handling, altered mitochondrial quality control, and contractile remodeling. This is compatible with known soleus atrophy, altered lipid metabolism, and slow-to-fast transition [10,11]. Consensus ranking reinforced this compact panel across soleus studies, but matched all-gene augmentation did not pass its joint gate. Soleus is therefore biologically coherent but methodologically secondary to thymus.

### Additional tissue results

Liver, skin, and spleen had narrower matched gene-level results. Liver contributed four FLT-lower genes without significant set enrichment. Skin *Plscr1* and spleen *Loxl1* were supported by both matched and consensus analyses. Pooled skeletal muscle improved at the classifier level but had no matched gene, while kidney, adrenal gland, gastrocnemius, and tibialis anterior were limited to consensus findings. Eye, retina, and lung improved prediction without a retained BH-FDR candidate. Follow-up should prioritize the thymus mitotic program, then test the smaller liver, skin, and spleen findings and the secondary soleus panel.

### Limitations

The tissue screens withheld profiles within represented accessions, so they measure utility within represented studies rather than transfer to a new mission, strain, or processing protocol. The configurable matrix was a gated search rather than a full factorial benchmark, and the selected path should not be interpreted as a universal model optimum.

The 974-gene landmark panel excludes potentially relevant genes. FDR was controlled within each tissue family, and heterogeneous study effects weaken a common-response interpretation even when the pooled association is significant. Bulk data also cannot separate cell-intrinsic regulation from changes in cell composition. Finally, the generator represents only tissue, condition, accession, and material contexts seen during adaptation. Independent studies and targeted experiments are required to confirm the proposed programs; detailed sensitivity analyses and safeguards are provided in the supplement.

## Conclusions

Synthetic expression did not improve one pooled multi-tissue classifier, but tissue-specific matched augmentation improved held-out-real performance in 16 analysis units without reducing any mean pooled or accession-macro metric. Twenty-one real-data BH-FDR associations also had synthetic-supported marginal importance. Thymus supplied 15 of them and a significant mitotic program, making it the strongest result. Liver, skin, and spleen supplied narrower matched findings.

The secondary consensus analysis retained broader correlated panels, including a coherent soleus mitochondrial and lipid-metabolism result, but those panels are not equivalent to matched synthetic contribution. The matched analysis shows where generated expression changed prediction; the consensus analysis organizes correlated biological hypotheses. Independent biological data are needed to test those hypotheses.

## Data and code availability

OSDR data were accessed through the public Biological Data API [1]. ARCHS4 and Reactome are public resources [2,7]. The code repository is <https://github.com/jasont314/nasa-mouse>. Analysis scripts are under `src/nasa_mouse_diffusion/paper_parity/`, frozen run outputs are under `outputs/generative/benchmark/`, and this paper package is under `paper/synthetic_guided_spaceflight/`. The figure and table builder is `nasa_mouse_diffusion.paper_parity.build_synthetic_guided_paper`. Annotation prompts, accepted rationales, and source links are recorded in `docs/annotation_prompts.md` and `docs/annotation_provenance.md`. No archival DOI has been assigned.

## Ethics statement

This study is a secondary computational analysis of publicly available animal-study data. No new animals were used.

## Competing interests

The author declares no competing interests.

## Acknowledgments

The author acknowledges the NASA Open Science Data Repository, GeneLab data-processing teams, original flight-study investigators, ARCHS4, Reactome, and the developers of the evaluated generative-model implementations. Program and mentor acknowledgments should be confirmed before submission.

## References

1. NASA Open Science Data Repository. OSDR Biological Data API. <https://visualization.osdr.nasa.gov/biodata/api/>. Accessed July 28, 2026.
2. Lachmann A, Torre D, Keenan AB, et al. Massive mining of publicly available RNA-seq data from human and mouse. *Nature Communications*. 2018;9:1366. <https://doi.org/10.1038/s41467-018-03751-6>.
3. Viñas R, Andrés-Terré H, Liò P, Bryson K. Adversarial generation of gene expression data. *Bioinformatics*. 2022;38:730-737. <https://doi.org/10.1093/bioinformatics/btab035>.
4. Lacan A, André R, Sebag M, Hanczar B. In silico generation of gene expression profiles using diffusion models. *BMC Bioinformatics*. 2026. <https://doi.org/10.1186/s12859-026-06470-8>.
5. Song J, Meng C, Ermon S. Denoising diffusion implicit models. *International Conference on Learning Representations*. 2021. <https://arxiv.org/abs/2010.02502>.
6. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57:289-300. <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>.
7. Milacic M, Beavers D, Conley P, et al. The Reactome Pathway Knowledgebase 2024. *Nucleic Acids Research*. 2024;52:D672-D678. <https://doi.org/10.1093/nar/gkad1025>.
8. Gridley DS, Mao XW, Stodieck LS, et al. Changes in mouse thymus and spleen after return from the STS-135 mission in space. *PLoS ONE*. 2013;8:e75097. <https://doi.org/10.1371/journal.pone.0075097>.
9. Horie K, Kato T, Kudo T, et al. Impact of spaceflight on the murine thymus and mitigation by exposure to artificial gravity during spaceflight. *Scientific Reports*. 2019;9:19866. <https://doi.org/10.1038/s41598-019-56432-9>.
10. Gambara G, Salanova M, Ciciliot S, et al. Gene expression profiling in slow-type calf soleus muscle of 30 days space-flown mice. *PLoS ONE*. 2017;12:e0169314. <https://doi.org/10.1371/journal.pone.0169314>.
11. Stein TP, Schluter MD, Grindeland RE, Moran MM, Baer LA, Wade CE. Rate controlling steps in fatty acid oxidation by unloaded rodent soleus muscle. *Journal of Gravitational Physiology*. 2002;9:P165-P166. PMID: 15002531.
12. Ho J, Jain A, Abbeel P. Denoising diffusion probabilistic models. *Advances in Neural Information Processing Systems*. 2020;33:6840-6851.
13. Heusel M, Ramsauer H, Unterthiner T, Nessler B, Hochreiter S. GANs trained by a two time-scale update rule converge to a local Nash equilibrium. *Advances in Neural Information Processing Systems*. 2017;30.
14. Kynkäänniemi T, Karras T, Laine S, Lehtinen J, Aila T. Improved precision and recall metric for assessing generative models. *Advances in Neural Information Processing Systems*. 2019;32.
15. DerSimonian R, Laird N. Meta-analysis in clinical trials. *Controlled Clinical Trials*. 1986;7:177-188. <https://doi.org/10.1016/0197-2456(86)90046-2>.
16. Cope H, Elsborg J, Demharter S, et al. Transcriptomics analysis reveals molecular alterations underpinning spaceflight dermatology. *Communications Medicine*. 2024;4:106. <https://doi.org/10.1038/s43856-024-00532-9>.
17. Finch RH, Vitry G, Siew K, et al. Spaceflight causes strain-dependent gene expression changes in the kidneys of mice. *npj Microgravity*. 2025;11:11. <https://doi.org/10.1038/s41526-025-00465-0>.
18. Siew K, et al. Cosmic kidney disease: an integrated pan-omic, physiological and morphological study into spaceflight-induced renal dysfunction. *Nature Communications*. 2024;15:4923. <https://doi.org/10.1038/s41467-024-49212-1>.
19. Suzuki N, Iwamura Y, Nakai T, et al. Gene expression changes related to bone mineralization, blood pressure and lipid metabolism in mouse kidneys after space travel. *Kidney International*. 2022;101:92-105. <https://doi.org/10.1016/j.kint.2021.09.031>.
20. Jeong W, Kwon H, Park SK, Lee IS, Jho EH. Retinoic acid-induced protein 14 links mechanical forces to Hippo signaling. *EMBO Reports*. 2024;25:4033-4061. <https://doi.org/10.1038/s44319-024-00228-0>.
21. Fearnley GW, Young KA, Edgar JR, et al. The homophilic receptor PTPRK selectively dephosphorylates multiple junctional regulators to promote cell-cell adhesion. *eLife*. 2019;8:e44597. <https://doi.org/10.7554/eLife.44597>.
22. Liu X, Zhao Y, Gao J, et al. Elastic fiber homeostasis requires lysyl oxidase-like 1 protein. *Nature Genetics*. 2004;36:178-182. <https://doi.org/10.1038/ng1297>.
23. Ilangovan H, Kothiyal P, Hoadley KA, et al. Harmonizing heterogeneous transcriptomics datasets for machine learning-based analysis to identify spaceflown murine liver-specific changes. *npj Microgravity*. 2024;10:61. <https://doi.org/10.1038/s41526-024-00379-3>.
24. Sanders LM, Chok H, Samson F, et al. Batch effect correction methods for NASA GeneLab transcriptomic datasets. *Frontiers in Astronomy and Space Sciences*. 2023;10:1200132. <https://doi.org/10.3389/fspas.2023.1200132>.
25. Dimitrieva S, Janssens R, Li G, et al. Biologically relevant integration of transcriptomics profiles from cancer cell lines, patient-derived xenografts, and clinical tumors using deep learning. *Science Advances*. 2025;11:eadn5596. <https://doi.org/10.1126/sciadv.adn5596>.
26. Allen DL, Bandstra ER, Harrison BC, et al. Effects of spaceflight on murine skeletal muscle gene expression. *Journal of Applied Physiology*. 2009;106:582-595. <https://doi.org/10.1152/japplphysiol.90780.2008>.
27. Gridley DS, Slater JM, Luo-Owen X, et al. Spaceflight effects on T lymphocyte distribution, function and gene expression. *Journal of Applied Physiology*. 2009;106:194-202. <https://doi.org/10.1152/japplphysiol.91126.2008>.
28. Keenan SN, Suriani ND, Fidelito G, et al. HSD17B11 regulates PLIN5-ATGL mediated lipolysis, but not hepatic lipid metabolism in mice. *Journal of Lipid Research*. 2025;66:100943. <https://doi.org/10.1016/j.jlr.2025.100943>.
29. Shi Y, Wang S, Yan Y, et al. ETV1 drives CD4+ T cell-mediated intestinal inflammation in inflammatory bowel disease through amino acid transporter Slc7a5. *Advanced Science*. 2026;13:e11595. <https://doi.org/10.1002/advs.202511595>.
30. Kitamura A, Maekawa Y, Uehara H, et al. A mutation in the immunoproteasome subunit PSMB8 causes autoinflammation and lipodystrophy in humans. *Journal of Clinical Investigation*. 2011;121:4150-4160. <https://doi.org/10.1172/JCI58414>.
