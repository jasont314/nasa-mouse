# Technical Report: Generative Models for NASA Mouse Bulk RNA-seq

- **Status:** consolidated through 2026-07-28
- **Primary data:** NASA OSDR API-derived mouse bulk RNA-seq and mouse ARCHS4
- **Models evaluated:** Lacan et al. DDIM, Vinas et al. conditional WGAN-GP, and
  GeneJEPA

## Executive Summary

The benchmark produced one defensible broad mouse generator and one defensible
spaceflight-conditioned generator:

1. The paper-architecture Lacan DDIM trained on mouse ARCHS4 learned broad
   tissue-conditioned expression. Synthetic-trained tissue balanced accuracy was
   0.869 on real held-out profiles, compared with 0.895 for real-trained
   classification. Its precision, recall, adversarial accuracy, and gene moments
   were strong, although gene-correlation agreement was 0.879 rather than the
   stricter 0.98 target.
2. The ARCHS4-pretrained, OSDR-adapted factorized DDIM passed every predeclared
   finite-sample fidelity gate on four locked-test generations. It can generate
   synthetic FLT or GC profiles for tissue, study, and material profiles represented
   during OSDR training. This is within-study interpolation, not generation for a
   novel accession.
3. The study-conditioned WGAN-GP did not pass joint fidelity. Its marginal support
   metrics were strong, but real and synthetic samples remained distinguishable
   with adversarial accuracy 0.636, outside the accepted 0.40 to 0.60 interval. Its
   locked test was therefore not opened.
4. GeneJEPA is a representation model with no expression decoder. The exact
   one-epoch mouse screen did not beat direct expression for tissue classification,
   so it was not advanced as a diffusion-guidance backbone. It cannot generate
   expression without adding a separate, non-paper GeneJEPA decoder.
5. Synthetic cohort augmentation did not improve FLT/GC prediction overall. On the
   locked OSDR test, real plus synthetic training reduced balanced accuracy from
   0.754 to 0.734. A complete six-fold skeletal-muscle confirmation was also
   negative.
6. Generated data were more useful as a low-weight feature-selection view than as
   replacement biological samples. Independent held-out-study confirmation was
   strong in thymus and mixed in lung. The thymus result recovered a coherent
   FLT-down cell-cycle module with Reactome FDR below 0.05.

The practical conclusion is that DDIM is the best model in this benchmark. It is
approved for controlled simulation within represented OSDR studies and for
carefully validated feature guidance. It is not approved as a source of independent
replicates, an unseen-study generator, or a default way to improve FLT/GC
classifiers.

## Scope And Questions

The work addressed three distinct questions that should not be conflated:

- **Broad expression generation:** can a model reproduce mouse bulk RNA-seq across
  tissues?
- **Conditional spaceflight generation:** can it generate realistic FLT and GC
  profiles for specified tissue and design covariates?
- **Downstream utility:** do generated profiles improve real-sample FLT/GC
  classification or help identify stable genes and pathways?

The benchmark included ARCHS4-only training, OSDR-only training, and ARCHS4
pretraining followed by OSDR adaptation. It tested pooled tissue-conditioned models,
study conditioning, material conditioning, per-tissue behavior, multiple
normalization and harmonization methods, direct synthetic augmentation, and
generated-informed feature selection.

## Evidence Policy

Results are interpreted in four evidence tiers:

| Tier | Definition | Appropriate claim |
|---|---|---|
| Locked | Model and gates fixed before a one-time held-out evaluation | Model-level acceptance or rejection |
| Independent | Entire test accession removed from generator and classifier fitting | Evidence of transfer to that accession |
| Development | Repeated nested splits with accessions represented on both sides | Within-study screening and hypothesis generation |
| Historical | Earlier proxies, post-hoc tables, or models that failed a gate | Engineering provenance only |

No scalar composite score is used. Correlation, precision, recall, F1, adversarial
accuracy, Frechet distance, diversity, memorization, and conditional-effect recovery
must pass independently. A high value on one metric cannot compensate for failure on
another.

The earlier per-tissue "best use" table is a development artifact. It compared
generated-only and real-plus-generated arms on validation data, including WGAN arms
that subsequently failed model-level fidelity. It is not the final model-selection
or biological-evidence table.

## Data Provenance

### NASA OSDR

OSDR expression and metadata were obtained through the repository's NASA OSDR API
workflow. The deprecated integrated raw OSDR H5 file was not used.

- 1,631 API profile rows were found.
- 21 technical replicate rows were aggregated.
- The resulting matrix contains 1,610 biological profiles:
  835 FLT and 775 ground control.
- Profiles span 75 OSDR accessions and 24 canonical material classes.
- The source matrix contains 48,694 shared mouse genes and RSEM expected counts.
- Tissue aliases such as `lvr` and `liver` were canonicalized before modeling.

The principal tissues were liver, skeletal muscle, skin, kidney, thymus, spleen,
lung, and retina. Smaller API-returned tissues remained eligible for pooled
conditioning and exploratory analyses. Skeletal-muscle material labels were retained
where available.

The accepted conditional DDIM used a within-accession, tissue, and condition split:
781 training, 536 validation, and 293 locked-test profiles. Every represented test
accession therefore had training examples. This design measures interpolation while
controlling study context; it does not measure unseen-accession transfer.

### ARCHS4

The local mouse ARCHS4 file contains 997,515 profiles and 53,511 genes. Metadata for
the full file were inspected rather than treating the older 24,428-profile
eight-tissue subset as the complete reference.

Three reusable cohorts were defined:

| ARCHS4 cohort | Profiles | GEO series | Intended use |
|---|---:|---:|---|
| Control-only | 23,614 | 3,213 | Strict healthy/control sensitivity |
| Healthy-preferred | 62,299 | 5,307 | Primary reference pool |
| Broad | 134,250 | 15,111 | Diversity sensitivity |

The paper-parity DDIM selected 17,244 unique healthy-preferred profiles across 20
tissue classes, split into 9,796 train, 2,448 validation, and 5,000 test profiles.
GeneJEPA's bounded exact-architecture screen used 43,744 training exposures.

## Model Implementations

| Model | Paper-faithful core | NASA adaptation | Capability |
|---|---|---|---|
| Lacan ModelDDIM | Official 227M-parameter residual DDIM, 1,000 steps, paper optimizer and 15,000 epochs | Mouse landmarks; factorized OSDR condition adapters and calibration | Native expression generation |
| Vinas WGAN-GP | Released conditional generator/critic topology and training defaults | Mouse data; ARCHS4 pretraining and OSDR fine-tuning | Native expression generation |
| GeneJEPA | 768-wide, 512-latent, 24-block predictive encoder | Mouse vocabulary and bulk profiles | Representation only |

### DDIM

The ARCHS4 model imports the official `ModelDDIM` implementation at inspected commit
`cde890154698fcea96c924804aaff04af3351b48`. Its fixed configuration includes:

- two 8,192-unit residual layers and 227,109,786 parameters;
- 1,000 quadratic diffusion steps with beta 0.0001 to 0.02;
- summed noise-prediction MSE with antithetic timesteps;
- Adam at learning rate 0.0004783833151836702;
- OneCycle scheduling, AMP, EMA 0.999, and deterministic DDIM sampling;
- batch size 2,048 and 15,000 epochs.

Full-transcriptome TPM is computed before selecting a deterministic 974-gene mouse
landmark panel. MaxAbs scaling is fitted on training profiles only.

The accepted OSDR extension adds factorized tissue, FLT/GC, study, and material
conditioning through a rank-512 domain LoRA adapter. It uses a 4,000-step domain
adaptation followed by a 1,000-step correlation-regularized refinement. Train-only
calibration aligns global and shrunk accession/tissue means and adds partial
positive missing-covariance residual noise. FLT/GC labels are not used to fit
calibration means or covariances.

### Conditional WGAN-GP

The WGAN uses the released Vinas et al. topology:

- 64-dimensional noise;
- two 256-unit ReLU layers in the generator and critic;
- learned embeddings for categorical covariates;
- RMSprop at `5e-4`;
- five critic updates per generator update;
- gradient penalty 10 and batch size 32.

The OSDR model conditions on tissue, FLT/GC, accession, material, sex, and muscle
group where available. It is a conditional sampler. It is not a paired
FLT-to-GC counterfactual model with subject-level identity preservation.

### GeneJEPA

GeneJEPA encodes masked expression into a self-supervised representation. The public
model exposes embeddings but no expression decoder. A GeneJEPA-guided generator
would therefore require a separate diffusion or decoder and would be a new method,
not paper-faithful GeneJEPA.

## Evaluation Metrics

The main metrics answer different failure modes:

- **Correlation agreement:** similarity between real and synthetic gene-gene
  correlation structure. Higher is better.
- **Precision:** fraction of generated support that lies near real support. Low
  precision indicates unrealistic or over-broad samples.
- **Recall:** fraction of real support covered by generated support. Low recall
  indicates mode loss.
- **F1:** harmonic mean of precision and recall.
- **Adversarial accuracy (AA):** held-out real-versus-synthetic discrimination.
  Approximately 0.5 is desirable; the accepted interval was 0.40 to 0.60.
- **Frechet distance (FD):** distance between real and synthetic embedding
  distributions. Lower is better. Conditional evaluations normalize FD to a
  real-split reference.
- **Reverse validation:** train on synthetic, test on real.
- **FLT/GC effect recovery:** correlation and direction agreement between real and
  synthetic FLT-minus-GC gene effects.
- **Accession-aware recovery:** FLT-minus-GC effects estimated within each study and
  then combined, preventing between-study differences from being mislabeled as
  spaceflight biology.
- **LOO-FDR stability:** the FDR requirement must survive omission of each accession
  in turn. This is stricter than ordinary pooled FDR.

Generated rows are never counted as new biological replicates in FDR testing.

## Broad Mouse Tissue Generation

The exact paper-architecture ARCHS4 DDIM completed 15,000 epochs and 75,000 optimizer
steps on an A100 in 5,987 seconds, with peak allocated GPU memory of 5.93 GB.

| Metric | ARCHS4 DDIM result | Interpretation |
|---|---:|---|
| Real-train to real-test tissue BA | 0.8954 | Real-data ceiling for the probe |
| Synthetic-train to real-test tissue BA | 0.8687 | Strong tissue utility |
| Gene mean correlation | 0.9965 | Pass |
| Gene standard-deviation correlation | 0.9437 | Strong |
| Gene-correlation agreement | 0.8791 | Below strict 0.98 target |
| Direct precision | 0.9655 | Pass |
| Direct recall | 0.8645 | Pass finite-sample rule |
| Nearest-neighbor AA | 0.512 | Real and synthetic are difficult to separate |
| PCA-50 precision / recall | 0.9855 / 0.9430 | Strong embedded support |
| PCA-50 FD | 0.0385 | Low |

The model learned tissue-specific generation despite overlap in the first two PCA
components. The final two-dimensional PCA silhouette was -0.271, which is not
inconsistent with the strong held-out tissue probe because a nonlinear or
higher-dimensional boundary can retain information that PC1 and PC2 do not display.

The raw model-scale output contains 9.04% negative entries. Those values are valid in
the unconstrained scaled model space but are not physical TPM. Any export for
downstream expression analysis must use an explicit nonnegative clipping or
nonnegative-generation policy and record it.

**Decision:** useful as a defensible broad 20-tissue mouse reference generator, but
not accepted as an all-metric finalist because correlation-structure agreement
missed the strict target. It is not by itself a spaceflight generator.

## OSDR Conditional Generation

### Accepted Factorized DDIM

Four independent synthetic cohorts were generated against the fixed 293-profile
locked test.

| Metric | Mean | Range | Repeats passing |
|---|---:|---:|---:|
| Gene-correlation agreement | 0.977 | 0.974 to 0.979 | 4/4 |
| Precision | 0.998 | 0.997 to 1.000 | 4/4 |
| Recall | 0.997 | 0.997 to 0.997 | 4/4 |
| F1 | 0.997 | 0.997 to 0.998 | 4/4 |
| Adversarial accuracy | 0.458 | 0.454 to 0.464 | 4/4 |
| FD / real-split P95 | 0.075 | 0.047 to 0.089 | 4/4 |

The test-size bootstrap correlation floor was 0.950, so all repeats pass the
finite-sample rule. They remain just below the separate 0.98 paper target.

Pooled FLT/GC effect recovery passed three of four generations. Mean effect
correlation was 0.598 and mean direction agreement was 0.683. The
skeletal-muscle accession-aware diagnostic passed all four generations, with mean
meta-effect correlation 0.608 and direction agreement 0.606 across five accessions.
Exact jointly LOO-FDR-stable real/synthetic genes remained sparse: 0, 0, 0, and 1
across the four generations.

This combination matters. The model reproduces the broad conditional distribution
and an aggregate FLT/GC contrast, but it does not reproduce a stable exact gene list.

### Per-Tissue Validation Screen

The table below reports DDIM FLT/GC effect correlation, direction agreement, and
repeat passes on shared validation data. It is useful for tissue prioritization but
is below the locked evidence tier.

| Tissue | Profiles / accessions | Tissue-level `r / direction / pass` | Accession-aware `r / direction / pass` |
|---|---:|---|---|
| Liver | 79 / 7 | 0.44 / 0.59 / 4/6 | 0.36 / 0.58 / 5/6 |
| Skeletal muscle | 61 / 11 | 0.57 / 0.77 / 6/6 | 0.24 / 0.58 / 3/6 |
| Skin | 50 / 6 | 0.72 / 0.71 / 6/6 | 0.62 / 0.69 / 6/6 |
| Kidney | 45 / 6 | -0.41 / 0.34 / 0/6 | -0.49 / 0.30 / 0/6 |
| Thymus | 37 / 4 | 0.68 / 0.72 / 6/6 | 0.76 / 0.63 / 6/6 |
| Spleen | 35 / 4 | -0.22 / 0.45 / 0/6 | -0.38 / 0.32 / 0/6 |
| Lung | 25 / 3 | 0.19 / 0.48 / 2/6 | 0.12 / 0.48 / 0/6 |
| Retina | 26 / 3 | 0.81 / 0.75 / 6/6 | 0.48 / 0.61 / 4/6 |

Skin and thymus were the strongest validation tissues after accession control.
Retina was strong without accession adjustment but less stable after it. Liver was
moderate. Kidney, spleen, and lung did not recover the contrast in this generation
screen. The locked-test muscle result was stronger than its earlier validation
screen, but exact gene overlap remained negligible.

### Classifier Utility

| Training data | Locked-test balanced accuracy | Locked-test ROC AUC |
|---|---:|---:|
| Real OSDR train | 0.754 | 0.819 |
| Synthetic train | 0.700 | 0.751 |
| Real plus synthetic | 0.734 | 0.801 |

Synthetic profiles contain condition information, but adding them did not improve
the broad real-only classifier.

**Decision:** accepted for conditional simulation inside represented studies.
Rejected for unseen-study generation, exact gene replication, and default
classification augmentation.

## WGAN-GP Findings

The strongest study-conditioned WGAN validation result used the same 974-gene
comparison space and six generation seeds.

| Metric | Mean | Gate | Result |
|---|---:|---:|---|
| Correlation agreement | 0.9759 | repeat-specific floor 0.9756 | 4/6 pass |
| Precision | 0.9764 | at least 0.95 | 6/6 pass |
| Recall | 0.9938 | at least 0.85 | 6/6 pass |
| F1 | 0.9850 | at least 0.90 | 6/6 pass |
| Adversarial accuracy | 0.6362 | 0.40 to 0.60 | 0/6 pass |
| FD / real-split P95 | 0.1439 | at most 1.0 | 6/6 pass |

Twelve train-only calibration variants were evaluated. Added residual variance moved
AA toward 0.5 but reduced correlation below its required floor. Per-gene quantile
calibration also failed, indicating a joint-distribution mismatch rather than only a
marginal-range problem.

Earlier WGAN transfer produced pooled FLT/GC effect correlation 0.805. That value
was misleading in isolation: accession-aware correlation was -0.022 and zero of nine
evaluable tissues passed. The WGAN could exploit pooled study-condition structure
without reproducing the within-study FLT/GC effect.

Adversarial training does not guarantee that an external classifier reaches chance.
The generator and training critic have finite capacity, optimize an approximate
Wasserstein objective, and can stop before the full joint distribution is matched.
External held-out discrimination is therefore still required.

**Decision:** rejected. The locked test remains unopened, and WGAN-generated
profiles should not be used for augmentation or biological claims.

## GeneJEPA Findings

The exact mouse screen used 4,096 train-selected highly variable genes and 43,744
ARCHS4 training exposures. It ran on an A100 with 31.68 GB peak allocated memory.
One epoch required 1,074 training seconds and made only 238 optimizer updates, below
the paper's 2,000-step EMA warmup.

| Held-out tissue probe | Balanced accuracy | Macro F1 |
|---|---:|---:|
| GeneJEPA embedding | 0.703 | 0.701 |
| Direct expression | 0.839 | 0.840 |

Embedding and UMAP silhouettes were -0.176 and -0.215. The screen was much shorter
than the paper's approximately 50 million profile exposures, so it is not evidence
that a full-scale GeneJEPA model must fail. A paper-scale single-A100 run was
estimated at roughly 20 GPU-days, however, and the bounded run did not justify that
cost for this benchmark.

**Decision:** retain GeneJEPA only as an optional representation experiment. Do not
describe it as a generator, and do not add a decoder under the GeneJEPA name without
labeling the resulting method as new.

## Harmonization Findings

A matched liver experiment held architecture, seed, 974 genes, training duration,
and data split fixed while changing preprocessing. Each arm used 119 train and 50
validation profiles; the 70-profile OSD-379 test was not opened.

| Method | Corr | Precision | Recall | F1 | AA | FD ratio |
|---|---:|---:|---:|---:|---:|---:|
| No harmonization, TPM | 0.283 | 0.200 | 0.960 | 0.331 | 0.850 | 2.052 |
| Ilangovan within-study z-score | 0.278 | 0.160 | 1.000 | 0.276 | 0.770 | 0.716 |
| Mentor two-stage z-score | 0.348 | 0.440 | 1.000 | 0.611 | 0.690 | 0.977 |
| ComBat | 0.004 | 0.040 | 1.000 | 0.077 | 0.810 | 63.185 |
| ComBat-seq | 0.067 | 0.020 | 1.000 | 0.039 | 0.850 | 156.647 |
| MBatch median polish | 0.009 | 0.020 | 1.000 | 0.039 | 0.930 | 205.470 |
| MBatch empirical Bayes | 0.001 | 0.000 | 1.000 | 0.000 | 0.850 | 60.625 |
| MBatch ANOVA | -0.003 | 0.020 | 1.000 | 0.039 | 0.870 | 44.716 |
| MOBER | 0.808 | 0.260 | 1.000 | 0.413 | 0.770 | 33.311 |

No method passed all fidelity criteria or either conditional-effect gate. MOBER's
0.808 correlation was the highest, but its low precision and F1, high AA, and
33-fold FD ratio show that the full distribution was not reproduced. The mentor
two-stage transform was the closest balanced alternative, but it also failed.

This does not prove that harmonization is generally harmful. It shows that none of
these adapters solved conditional generation on this fixed liver benchmark. ComBat
and MBatch also remain transductive sensitivity methods for unseen batches, whereas
MOBER is the principal inductive complex harmonizer tested here.

## Downstream FLT/GC Utility

### Direct Synthetic Augmentation

Three independent checks reached the same practical conclusion:

- Broad locked test: real plus synthetic BA 0.734 versus real-only 0.754.
- Adaptive per-tissue screen: heart improved on a five-profile, one-accession test;
  retina failed; skeletal muscle tied at a ceiling.
- Complete skeletal-muscle confirmation across 11 held-out accessions and 159
  profiles: accession-macro BA changed from 0.655 to 0.658, AUROC from 0.718 to
  0.690, and AP from 0.768 to 0.764. The BA bootstrap interval was -0.022 to 0.030.

Generated cohorts therefore do not provide a validated general augmentation gain.
They must not be treated as added sample size in DGEA or pathway FDR.

### Generated-Informed Feature Selection

The more successful use was to let one or more synthetic views influence feature
ranking while keeping the final classifier anchored to real profiles or assigning
synthetic rows very low total weight.

An independent confirmation removed OSD-900 lung and OSD-457 thymus from every OSDR
generator fine-tuning role before evaluating them:

| Tissue | Held-out study | Profiles | Baseline BA / AUC / AP | Guided BA / AUC / AP |
|---|---|---:|---|---|
| Lung | OSD-900 | 20 | 0.400 / 0.450 / 0.523 | 0.550 / 0.550 / 0.635 |
| Thymus | OSD-457 | 24 | 0.500 / 0.840 / 0.876 | 0.833 / 0.972 / 0.976 |

Across both studies, 11 predictions changed from incorrect to correct and none
changed from correct to incorrect. Only one accession was tested per tissue, so this
is not a study-level significance test.

The genotype audit strengthened thymus: metrics improved in both Nrf2KO and WT
subgroups. Lung was mixed because KO AUROC fell from 0.56 to 0.52 even though its BA
and AP improved.

The thymus signature contained `Birc5`, `Ccne2`, `Gmnn`, `Ube2c`, `Cdk1`, `Nusap1`,
`Ccnb1`, and `Ccnb2`, consistently FLT-down in both genotypes. Reactome terms for
G2/M checkpoints, mitotic protein degradation, DNA replication, and cell-cycle
control passed FDR 0.05. Lung highlighted `Cdkn1a`, `Ccne2`, `Pik3r3`, `Slc25a4`,
`Mapk9`, and `Igfbp3`, but no lung Reactome term passed FDR 0.05.

### Exploratory Within-Study Stability

Repeated nested development splits were run for 22 tissues. Because accessions were
represented on both sides and the generator had seen the original training role,
these are exploratory within-study results.

| Tissue | Selected use | Mean BA / AUC / AP delta | Repeats all metrics nonworse | Real LOO-FDR stable genes |
|---|---|---|---:|---:|
| Kidney | Guided, low weight | +0.029 / +0.093 / +0.097 | 6/8 | 0 |
| Liver | Guided, low weight | +0.024 / +0.054 / +0.044 | 5/8 | 0 |
| Lung | Generated only | +0.086 / +0.156 / +0.157 | 7/8 | 0 |
| Retina | Guided, low weight | +0.117 / +0.121 / +0.071 | 7/8 | 0 |
| Skeletal muscle | Real plus generated | +0.043 / 0.000 / +0.013 | 5/8 | 7 |
| Skin | Generated only | +0.086 / +0.109 / +0.065 | 4/8 | 0 |
| Spleen | Guided, real only | +0.170 / +0.208 / +0.204 | 7/8 | 1 |
| Thymus | Generated only | +0.121 / +0.092 / +0.075 | 7/8 | 8 |

The seven muscle LOO-stable genes were `Sox4`, `Sh3bp5`, `Cebpd`, `Cdkn1a`, `Bphl`,
`Prkcd`, and `Arid5b`. The eight thymus genes were `Cenpe`, `Ccnb1`, `Nusap1`,
`Stmn1`, `Cdk1`, `Top2a`, `Ccnb2`, and `Ccne2`. Spleen retained only `Igfbp3`.

The 292 significant Reactome rows produced across all tissues and overlapping gene
sets are not 292 independent pathway discoveries. Reactome hierarchy, pathway
overlap, and nested feature sets create extensive dependence.

### Skeletal-Muscle Group Follow-Up

The synthetic-guided workflow was subsequently rerun separately for EDL,
gastrocnemius, quadriceps, soleus, and tibialis anterior. Soleus was the strongest
group: seven synthetic-selected genes passed real-data LOO FDR 0.05, and the
selected set converged on lower flight expression of mitochondrial fatty-acid
oxidation and related muscle-regulatory genes. Quadriceps retained one LOO-stable
gene; the remaining groups were exploratory or limited to two accessions. This is
within-study evidence, not unseen-accession validation. See the
[muscle-group analysis](synthetic_skeletal_muscle_group_analysis.md).

## What The Generator Can And Cannot Do

### Supported

- Generate broad mouse tissue-conditioned expression from the ARCHS4 DDIM.
- Generate FLT or GC profiles for OSDR tissue, accession, and material combinations
  represented during factorized-DDIM adaptation.
- Reproduce the overall locked OSDR distribution with strong precision, recall,
  AA, and FD.
- Recover a moderate pooled FLT/GC effect and a locked skeletal-muscle
  accession-aware effect.
- Provide a useful secondary feature-ranking view, with the strongest independent
  evidence currently in thymus.

### Not Supported

- Generation for a new OSDR study with no training representation.
- Treating generated profiles as independent animals or biological replicates.
- Exact recovery of a stable FLT/GC gene list.
- General improvement from adding large synthetic cohorts to classifier training.
- Subject-paired FLT-to-GC counterfactuals with all non-condition attributes fixed.
- WGAN use for final generation or biological interpretation.
- GeneJEPA expression generation without a new decoder.

## Recommended Operating Protocol

1. Freeze the accepted factorized DDIM and its preprocessing for reported
   within-study simulation. Do not choose new calibration or hyperparameters using
   the locked test.
2. For downstream FLT/GC analysis, use generated profiles for feature ranking or
   weak regularization only. Fit the final model primarily on real profiles.
3. Evaluate every biological claim with within-accession real-data effects,
   random-effects meta-analysis, ordinary FDR, and leave-one-accession-out
   sensitivity. Synthetic rows must not enter the biological replicate count.
4. Keep OSD-900 and OSD-457 untouched for further tuning. Replicate the thymus
   feature policy in a newly obtained thymus accession and treat lung as
   exploratory until its direction is stable across genotype and study.
5. Extend API metadata ingestion to capture genotype, sex, material, and muscle
   group directly. Avoid reconstructing these fields from profile names where the
   source API or repository metadata can provide them.
6. Declare the nonnegative export policy for DDIM outputs. Store both the native
   scaled output and the clipped inverse-scaled expression with provenance.
7. For true unseen-study generation, design a new accession-held-out adaptation
   benchmark. Study conditioning should include an explicit unknown-study or
   hierarchical study representation rather than reusing a known accession label.
8. Do not spend full-scale GeneJEPA compute unless a representation task first
   demonstrates value over direct expression on a fresh benchmark.

## Reproducibility And Artifact Index

- [Pipeline design](generative_models_pipeline.md)
- [Paper and code audit](generative_model_code_audit.md)
- [Data audit](generative_data_audit.md)
- [Canonical benchmark results](generative_benchmark_results.md)
- [Machine-readable model scoreboard](../outputs/generative_benchmark/summary/model_scoreboard.tsv)
- [DDIM paper-parity implementation](rna_diffusion_paper_parity.md)
- [Accepted OSDR conditional DDIM](osdr_conditional_diffusion_finalist.md)
- [Rejected study-conditioned WGAN](osdr_conditional_wgan_study.md)
- [Generated-feature workflow](generated_feature_guidance_workflow.md)
- [ARCHS4 DDIM run](../outputs/generative_benchmark/runs/lacan_diffusion/archs4_mouse_paper_parity_seed1234/)
- [Conditional DDIM locked test](../outputs/generative_benchmark/runs/lacan_diffusion/osdr_factorized_study_lora512_correlation_refine_seed2020/evaluation/final_locked_test/)
- [DDIM versus WGAN tissue comparison](../outputs/generative_benchmark/comparisons/ddim_vs_wgan_study_conditioned_validation_seed3020/)
- [Independent feature confirmation](../outputs/generative_benchmark/analyses/generated_feature_guidance_confirmation_v1/)
- [Within-study feature stability](../outputs/generative_benchmark/analyses/within_study_generated_feature_stability_v1/)
- [Synthetic-guided skeletal-muscle groups](synthetic_skeletal_muscle_group_analysis.md)
- [Skeletal-muscle augmentation confirmation](../outputs/generative_benchmark/analyses/fresh_holdout_contrastive_ddim_augmentation_v1/)
- [Matched liver harmonization results](../outputs/generative_benchmark/summary/liver_harmonization/)

The older reduced diffusion proxy, scalar-composite rankings, pooled WGAN effect
correlation, and per-tissue "best use" table remain useful for provenance but are
superseded for final decisions by the paper-parity DDIM, independent metric gates,
locked conditional test, complete augmentation confirmation, and independent
feature-guidance confirmation summarized here.
