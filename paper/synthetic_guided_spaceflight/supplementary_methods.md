<div class="title-page">

<h1>Supplementary methods</h1>

<p class="subtitle">A configurable generative transcriptomics framework for spaceflown mice</p>

<p class="authors">Jason Trinh</p>

<p class="draft-note"><strong>Supplementary methods.</strong> This document records the data roles, model architecture, evaluation gates, statistical safeguards, and complete supporting results.</p>

</div>

## S1. OSDR discovery and expression ingestion

OSDR data were obtained through the Biological Data API documented at
<https://visualization.osdr.nasa.gov/biodata/api/>.

Eligibility required:

- `Mus musculus`;
- transcription profiling by bulk RNA sequencing;
- a resolvable spaceflight/flight or ground-control label;
- processed RSEM expected-count output;
- sample-level tissue or material metadata sufficient for canonicalization.

No raw combined OSDR HDF5 file was used. The API returned 1,631 profile rows. Twenty-one technical replicates were
aggregated, leaving 1,610 biological profiles, 835 flight and 775 ground control,
from 75 accessions and 24 canonical material classes.

## S2. ARCHS4 cohort audit

The source file was `assets/archs4/mouse_gene_v2.5.h5`, containing 997,515
profiles and 53,511 genes.

The three eligible reference cohorts contained:

| Cohort | Profiles | GEO series | Intended use |
|---|---:|---:|---|
| Control only | 23,614 | 3,213 | Conservative sensitivity cohort |
| Healthy preferred | 62,299 | 5,307 | Primary pretraining source |
| Broad | 134,250 | 15,111 | Diversity sensitivity cohort |

Before reference selection, nine GEO series linked to eligible OSDR accessions
were excluded by accession metadata. This removed 108 otherwise eligible ARCHS4
profiles, including all 53 selected profiles from GSE152382 (OSD-457). The
paper-parity run then selected 17,244 healthy-preferred profiles across 20 tissue
classes. Complete GEO-series grouping assigned 10,150 profiles to training,
2,466 to validation, and 4,628 to test, with no series shared across splits.

## S3. Landmark panel and normalization

Full-transcriptome TPM was calculated with GENCODE vM39 mouse gene lengths.
Landmark selection occurred after TPM calculation and used a deterministic
974-gene mouse mapping of the human L1000 panel. Training-set MaxAbs scaling was
applied after landmark selection.

## S4. ARCHS4 DDIM configuration

The full configuration is
`configs/rna_diffusion/archs4_mouse_paper_parity_osdr_disjoint.yaml`. The neural
architecture and optimizer settings match the paper-parity implementation; only
the reference exclusion and grouped split contract changed.

| Component | Value |
|---|---|
| Input genes | 974 |
| Hidden layers | 8,192; 8,192 |
| Parameters | 227,109,786 |
| Dropout | 0.1 |
| Diffusion steps | 1,000 |
| Beta schedule | Quadratic, 0.0001 to 0.02 |
| Objective | Summed noise MSE |
| Optimizer | Adam |
| Learning rate | 0.0004783833151836702 |
| Scheduler | OneCycle |
| Batch size | 2,048 |
| Epochs / optimizer steps | 15,000 / 75,000 |
| AMP | Enabled |
| EMA | 0.999 |
| Device | NVIDIA A100-SXM4-40GB |
| Runtime | 6,083 seconds |
| Peak allocated GPU memory | 5.93 GB |

## S5. Factorized OSDR adaptation

The accepted configuration is
`configs/rna_diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint.yaml`.
It uses the ARCHS4 backbone described in Section S4. A base adapter was trained for 12,000 domain and 4,000 condition steps,
followed by the 4,000-domain-step and 1,000-condition-step correlation-refinement stage
reported here. The all-tissue and muscle-group development screens were then
regenerated from this adapter.

| Component | Value |
|---|---|
| Backbone | Completed ARCHS4 paper-parity DDIM |
| Conditioning | Tissue, FLT/GC, accession, material type |
| Domain adapter | LoRA rank 512, alpha 512 |
| Domain stage | 4,000 steps, LR 0.00002 |
| Condition stage | 1,000 steps, LR 0.00002 |
| Batch size | 512 |
| Condition dropout | 0.15 |
| Correlation regularization | Weight 10, 256 genes, timestep <= 200 |
| Sampling | 100 DDIM steps for evaluation |

The data split contained 781 training, 536 validation, and 293 locked-test
profiles. Every test accession was represented in training. Results therefore
measure within-study interpolation.

The accepted calibrator used:

1. train-only global alignment;
2. hierarchically shrunk accession and tissue means;
3. positive missing-covariance residual noise;
4. no condition-specific fit of calibrator means or covariance;
5. explicit clipping at zero for accepted downstream expression.

## S6. Distribution and condition gates

Four generation seeds, 5020-5023, were declared before the locked test was
opened. Every metric was gated independently.

| Metric | Gate |
|---|---|
| Gene-correlation agreement | Paper target >= 0.98; finite-sample real-bootstrap floor also reported |
| Precision | >= 0.95 |
| Recall | >= 0.85 |
| F1 | >= 0.90 |
| Adversarial accuracy | 0.40 to 0.60 |
| FD / real-split P95 | <= 1.0 |
| Pooled effect recovery | Correlation >= 0.30 and direction >= 0.55 |
| Muscle accession recovery | Correlation >= 0.30 and direction >= 0.55 |
| Memorization | Generated fraction below train LOO P01 <= 0.05 |

The exact repeat rows are in `source_data/table_s2_locked_ddim_repeats.tsv`.
On the locked within-study test, mean gene-correlation agreement was
0.9744, precision 0.9974, recall 0.9957, F1 0.9966, adversarial accuracy
0.4753, and FD/real-split-P95 ratio 0.0740. All four repeats passed the
finite-sample correlation floor of 0.9497, although the mean remained below the
paper target of 0.98. Pooled FLT/GC effect recovery passed in three of four
repeats and muscle accession-effect recovery passed in all four. The preceding
validation screen missed its stricter sample-specific correlation floor and
muscle accession-effect gate; the broad biological screens therefore remain
developmental.

The term "adversarial accuracy" refers to an external nearest-neighbor
real-versus-synthetic classifier, not the WGAN training critic. A result near 0.5
indicates that this external discriminator cannot reliably separate the two
cohorts.

## S7. Configurable benchmark and comparator models

The resolved planner contains 463 gated experiment rows. This is a staged plan,
not evidence that all 463 rows received paper-duration training. Smoke tests,
paper-architecture feasibility runs, preprocessing screens, accession-scope
tests, harmonization comparisons, study-conditioning experiments, and tissue
expansion were advanced only when the preceding gates supported the next cost.
The framework separates input units, library normalization, transformation,
scaling, feature space, harmonization, training regime, accession scope, tissue
mode, FLT/GC conditioning, study policy, balancing, seed, and synthetic-to-real
ratio.

### Liver harmonization benchmark

Nine arms used the same 974-gene conditional DDIM architecture, seed, and split:
119 training profiles, 50 validation profiles, and a 70-profile locked test that
was not opened for this screen. Every arm ran on an NVIDIA A100. No arm passed
all independent fidelity, pooled-condition, and accession-effect gates.

| Method | Corr. | Precision | Recall | F1 | AA | FD/real P95 | Fidelity gate | Accession gate |
|---|---:|---:|---:|---:|---:|---:|---|---|
| No harmonization, TPM | 0.283 | 0.200 | 0.960 | 0.331 | 0.850 | 2.052 | Fail | Fail |
| Ilangovan study z-score | 0.278 | 0.160 | 1.000 | 0.276 | 0.770 | 0.716 | Fail | Fail |
| Mentor two-stage z-score | 0.348 | 0.440 | 1.000 | 0.611 | 0.690 | 0.977 | Fail | Fail |
| ComBat | 0.004 | 0.040 | 1.000 | 0.077 | 0.810 | 63.185 | Fail | Fail |
| ComBat-seq | 0.067 | 0.020 | 1.000 | 0.039 | 0.850 | 156.647 | Fail | Fail |
| MBatch Median Polish | 0.009 | 0.020 | 1.000 | 0.039 | 0.930 | 205.470 | Fail | Fail |
| MBatch Empirical Bayes | 0.001 | 0.000 | 1.000 | 0.000 | 0.850 | 60.625 | Fail | Fail |
| MBatch ANOVA | -0.003 | 0.020 | 1.000 | 0.039 | 0.870 | 44.716 | Fail | Fail |
| MOBER | 0.808 | 0.260 | 1.000 | 0.413 | 0.770 | 33.311 | Fail | Fail |

ComBat, ComBat-seq, and the three MBatch validation transformations used training
anchors but remain transductive sensitivity analyses for a new batch. MOBER
was the inductive complex harmonizer. Its high correlation did not compensate
for low precision and F1, external separability, or excessive distributional
distance. Exact values and preprocessing labels are in Tables S13-S14.

### WGAN-GP and GeneJEPA screens

The WGAN used the Viñas et al. topology: 64-dimensional noise, two 256-unit
generator and critic layers, five critic updates, gradient-penalty weight 10,
RMSProp learning rate 0.0005, and batch size 32. The strongest study-conditioned
validation result was evaluated across six sampling seeds. Mean correlation,
precision, recall, and F1 were 0.9759, 0.9764, 0.9938, and 0.9850. Mean external
adversarial accuracy was 0.6362 and FD/real-P95 ratio was 0.1439. Pooled FLT/GC
effect recovery passed in six of six repeats, but accession-aware skeletal-muscle
recovery passed in zero of six. Because no repeat passed the full fidelity gate
and accession-aware recovery was unstable, its locked test remained unopened.
Exact repeat rows are in Table S15.

The exact-architecture GeneJEPA duration screen used 4,096 genes and 43,744
replacement-sampled training exposures. It reached 0.703 held-out tissue balanced
accuracy versus 0.839 from expression. It is representation-only and has no
expression decoder.

## S8. Generated-feature workflow

The evaluation funnel had three stages: pooled utility, tissue-specific
development, and real-data association testing. The pooled utility benchmark compared real-only, generated-only, and
real-plus-generated training. Tissue-specific development then compared five
arms:

1. `real_only`;
2. `generated_only`;
3. `real_plus_generated`, equal total real and synthetic weight;
4. `guided_real_only`, real classifier with real/synthetic consensus ranking;
5. `guided_low_weight`, real plus recentered synthetic profiles at 0.05 total
   synthetic weight.

Splits were nested within accession-by-condition strata. The inner loop selected
feature count, regularization, and rank method. The outer loop measured balanced
accuracy, AUROC, and average precision separately. No composite score was used.
Because the same accessions could contribute profiles to training and testing,
this stage measured within-study development utility. It did not test transfer
to a new study.

Stable genes had selection frequency at least 0.50 and coefficient-sign agreement
at least 0.75. Generated-supported status did not constitute biological evidence.
Real random-effects and LOO tests were applied afterward.

For the BH-FDR synthetic-informed subset, `reinforced` denotes a
`core_intersection` gene that was stable in both the real-only and selected
synthetic-guided arms, had matching coefficient directions, and had a supporting
real meta-effect. `Synthetic-promoted` denotes a `generated_supported` gene that
was stable in the selected synthetic-guided arm but did not cross the real-only
stability threshold, with a supporting real meta-effect. Thus
`synthetic-promoted` is a feature-selection classification, not a claim that the
gene was absent from real expression, statistically nonsignificant in real data,
or biologically novel.

<!-- BEGIN GENERATED TISSUE UTILITY TABLES -->

All five arms were fitted for every analysis unit below. Values are means across eight repeated outer splits, and every outer evaluation used real profiles. An eligible arm was nonworse than real-only training in balanced accuracy, AUROC, and average precision. An eligible tie met that rule without improving a mean metric. These development results use profiles from represented accessions.

**Supplementary Table S9. Complete canonical-tissue utility screen.**

| Tissue | n (FLT/GC) | Selected arm | BA real/selected | AUROC real/selected | AP real/selected | Status |
|---|---|---|---|---|---|---|
| Adrenal gland | 31 (16/15) | Generated only | 0.781 / 0.922 | 0.906 / 0.984 | 0.918 / 0.988 | Eligible improvement |
| Bone | 30 (15/15) | Guided ranking; real fit | 0.656 / 0.734 | 0.797 / 0.828 | 0.858 / 0.865 | Eligible improvement |
| Bone marrow | 20 (10/10) | Generated only | 0.969 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | Eligible improvement |
| Brain | 45 (22/23) | Generated only | 0.688 / 0.740 | 0.712 / 0.802 | 0.731 / 0.788 | Eligible improvement |
| Brown adipose tissue | 20 (10/10) | Generated only | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | Eligible tie |
| Cecum | 16 (8/8) | Real only | 0.875 / 0.875 | 0.938 / 0.938 | 0.958 / 0.958 | Real-only retained |
| Cerebellum | 44 (23/21) | Generated only | 0.529 / 0.602 | 0.621 / 0.696 | 0.729 / 0.764 | Eligible improvement |
| Colon | 45 (23/22) | Real only | 0.656 / 0.656 | 0.677 / 0.677 | 0.730 / 0.730 | Real-only retained |
| Eye | 18 (9/9) | Generated only | 0.594 / 0.781 | 0.562 / 0.781 | 0.729 / 0.865 | Eligible improvement |
| Heart | 42 (21/21) | Generated only | 0.594 / 0.719 | 0.625 / 0.750 | 0.642 / 0.778 | Eligible improvement |
| Hippocampus | 30 (16/14) | Generated only | 0.698 / 0.703 | 0.740 / 0.781 | 0.849 / 0.864 | Eligible improvement |
| Kidney | 111 (56/55) | Guided ranking; 5% synthetic | 0.654 / 0.707 | 0.680 / 0.771 | 0.683 / 0.799 | Eligible improvement |
| Liver | 197 (101/96) | Real only | 0.721 / 0.721 | 0.764 / 0.764 | 0.773 / 0.773 | Real-only retained |
| Lung | 63 (32/31) | Generated only | 0.664 / 0.742 | 0.680 / 0.830 | 0.684 / 0.833 | Eligible improvement |
| Mammary gland | 20 (8/12) | Generated only | 0.771 / 0.885 | 0.854 / 0.958 | 0.875 / 0.958 | Eligible improvement |
| Optic nerve | 29 (19/10) | Guided ranking; 5% synthetic | 0.769 / 0.819 | 0.863 / 0.950 | 0.950 / 0.983 | Eligible improvement |
| Retina | 63 (37/26) | Guided ranking; 5% synthetic | 0.587 / 0.723 | 0.620 / 0.762 | 0.756 / 0.846 | Eligible improvement |
| Skeletal muscle, pooled | 149 (74/75) | Guided ranking; real fit | 0.861 / 0.932 | 0.924 / 0.961 | 0.908 / 0.945 | Eligible improvement |
| Skin | 125 (66/59) | Real + generated | 0.671 / 0.756 | 0.691 / 0.768 | 0.747 / 0.808 | Eligible improvement |
| Spleen | 86 (44/42) | Real + generated | 0.517 / 0.648 | 0.540 / 0.704 | 0.589 / 0.749 | Eligible improvement |
| Thymus | 94 (51/43) | Guided ranking; 5% synthetic | 0.733 / 0.844 | 0.818 / 0.887 | 0.862 / 0.918 | Eligible improvement |
| White adipose tissue | 20 (10/10) | Generated only | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | Eligible tie |

**Supplementary Table S6. Complete anatomical muscle-group utility screen.**

| Tissue | n (FLT/GC) | Selected arm | BA real/selected | AUROC real/selected | AP real/selected | Status |
|---|---|---|---|---|---|---|
| EDL | 24 (12/12) | Real only | 0.917 / 0.917 | 1.000 / 1.000 | 1.000 / 1.000 | Real-only retained |
| Gastrocnemius | 25 (10/15) | Guided ranking; 5% synthetic | 0.604 / 0.646 | 0.646 / 0.750 | 0.731 / 0.798 | Eligible improvement |
| Quadriceps | 35 (18/17) | Real only | 0.750 / 0.750 | 0.880 / 0.880 | 0.916 / 0.916 | Real-only retained |
| Soleus | 41 (22/19) | Real + generated | 0.925 / 0.963 | 0.980 / 0.980 | 0.980 / 0.986 | Eligible improvement |
| Tibialis anterior | 24 (12/12) | Real + generated | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | Eligible tie |

Sample counts are shown as total development profiles followed by flight/ground-control counts. Small cohorts and ceiling-level scores remain exploratory even when a synthetic arm is eligible.

<!-- END GENERATED TISSUE UTILITY TABLES -->

## S9. Random-effects reporting and LOO sensitivity

For gene \(g\) in accession \(a\), the real flight effect was:

```text
delta[g,a] = mean(real expression[g] | FLT,a)
           - mean(real expression[g] | GC,a)
```

Accession effects were combined with a random-effects model. Within each tissue,
Benjamini-Hochberg (BH) adjustment was applied to the 974 gene-level
meta-analysis P values. BH sorts the P values and controls the expected proportion
of false discoveries among the rejected hypotheses. The primary inclusion rule
was BH FDR < 0.05; synthetic selection, study-direction agreement, heterogeneity,
and LOO stability were not additional significance gates.

After BH correction, each association was annotated as:

- reinforced by stable real-only and synthetic-guided selection;
- synthetic-promoted;
- selected by the real-only arm;
- synthetic-selected without real-direction support; or
- not stably selected by either feature-selection arm.

The accession-direction fraction, random-effects variance (`tau2`), heterogeneity
(`I2`), and exact study counts remain in the complete inventory. A gene with mixed
study directions remains BH-significant when its pooled random-effects test passes,
but the disagreement qualifies interpretation of the pooled effect.

The LOO analysis removed each accession, repeated the random-effects fit, and
retained the maximum FDR and any sign reversal. LOO was treated as a sensitivity
label rather than an inclusion requirement. A LOO-stable gene additionally
required maximum LOO FDR < 0.05 and no LOO sign reversal. Failure of this label
does not remove a gene from the primary BH-FDR table.

Generated profiles were never included in the random-effects model.

## S10. Reactome analysis

The official mouse GMT was generated from `ReactomePathways.txt` and
`Ensembl2Reactome_All_Levels.txt`, restricted to *Mus musculus*, `R-MMU-*`
pathways, and `ENSMUSG*` genes.

Hypergeometric enrichment used the 974-gene landmark panel as background.
Benjamini-Hochberg FDR was applied separately by tissue and selected-gene set.
Reactome parent and child terms overlap. Counts of significant rows are therefore
not counts of independent biological discoveries.

## S11. Skeletal-muscle group analysis

The factorized DDIM and three frozen synthetic development views were
reused. No additional neural network was trained for the muscle-group screen.

| Group | Profiles | Accessions | FLT | GC |
|---|---:|---:|---:|---:|
| EDL | 24 | 2 | 12 | 12 |
| Gastrocnemius | 25 | 3 | 10 | 15 |
| Quadriceps | 35 | 4 | 18 | 17 |
| Soleus | 41 | 3 | 22 | 19 |
| Tibialis anterior | 24 | 2 | 12 | 12 |

The selected soleus arm was real plus generated expression. It improved mean
balanced accuracy from 0.9250 to 0.9625, left AUROC at 0.9800, and increased
average precision from 0.9804 to 0.9865. Five BH-FDR genes were reinforced by
both stable real-only and selected-arm feature ranking and had the same real
effect direction in all three accessions: FLT-lower `Bdh1`, `Ech1`, `Bnip3`,
and `Decr1`, and FLT-higher `Tpm1`. `Bdh1`, `Ech1`, `Bnip3`, and `Tpm1` also
passed the LOO sensitivity criterion; `Decr1` did not. The screen did
not promote `Mef2c` or `Pxmp2`.

Gastrocnemius selected a guided low-weight arm and promoted `Fhl2` and `Nfkbia`.
Tibialis anterior selected real plus generated expression and reinforced
`Cdkn1a`, `St3gal5`, and `Bnip3` while promoting `Cebpd`. None of these
gastrocnemius or tibialis genes passed the LOO FDR sensitivity rule. EDL and
quadriceps retained real-only arms, so their BH-FDR genes are not attributed to
synthetic guidance.

LOO here is a real-data meta-analysis sensitivity test. It does not remove the
accession from the already completed generator adaptation. Soleus remains
developmental until a new accession is excluded from adaptation and selection.

## S12. Random-effects BH-FDR gene inventories

Supplementary Table S11 is the primary statistical inventory. It contains every
real-data random-effects association with BH FDR < 0.05, without requiring
synthetic selection, unanimous study direction, or LOO stability. Supplementary
Table S12 reports the corresponding counts for every tested analysis unit,
including tissues with zero BH-FDR genes.

The 459 retained tissue-gene results comprise 202 associations across 10 of 22
canonical tissues and 257 across all five anatomical muscle groups. This is an
inventory size, not a count of independent discoveries: genes can recur across
tissues, and pooled skeletal muscle overlaps the anatomical subgroup analyses. Of
these results, 363 had the same effect direction in every represented accession
and 96 did not. Direction agreement is reported as an annotation, not an
exclusion rule.

| Analysis unit | BH-FDR genes | FLT higher | FLT lower | Unanimous direction | Synthetic-promoted | Reinforced |
|---|---:|---:|---:|---:|---:|---:|
| Adrenal gland | 22 | 1 | 21 | 21 | 1 | 1 |
| Eye | 8 | 4 | 4 | 8 | 0 | 1 |
| Heart | 4 | 2 | 2 | 4 | 0 | 0 |
| Kidney | 3 | 3 | 0 | 1 | 1 | 1 |
| Liver | 19 | 8 | 11 | 0 | 0 | 0 |
| Retina | 5 | 5 | 0 | 5 | 0 | 0 |
| Skeletal muscle, pooled | 26 | 12 | 14 | 5 | 4 | 8 |
| Skin | 3 | 2 | 1 | 1 | 1 | 0 |
| Spleen | 34 | 32 | 2 | 21 | 3 | 1 |
| Thymus | 78 | 37 | 41 | 57 | 13 | 3 |
| EDL | 136 | 14 | 122 | 130 | 0 | 0 |
| Gastrocnemius | 15 | 8 | 7 | 13 | 2 | 0 |
| Quadriceps | 29 | 25 | 4 | 22 | 0 | 0 |
| Soleus | 29 | 5 | 24 | 27 | 0 | 5 |
| Tibialis anterior | 48 | 42 | 6 | 48 | 1 | 3 |

Supplementary Table S10 is the 49-row synthetic-informed subset of the primary
inventory: 26 synthetic-promoted and 23 reinforced tissue-gene results with
supporting real effects. Synthetic attribution is suppressed when the candidate
generated arm failed the three-metric eligibility gate. The table is retained to
describe what the generator changed in feature ranking, but it is not the primary
significance table. The full inventory additionally contains 34 real-only
selected results, 370 BH-significant results not stably selected by either arm,
and six synthetic-selected results without real-direction support.

The complete synthetic-informed gene set is shown below. `FLT higher` and `FLT
lower` refer to the sign of the real random-effects meta-estimate. An asterisk
marks a BH-significant pooled effect whose direction was not identical in every
accession.

**Synthetic-promoted genes**

| Analysis unit | FLT higher | FLT lower |
|---|---|---|
| Adrenal gland | None | `Psmb8` |
| Kidney | `Inpp4b`\* | None |
| Skeletal muscle, pooled | None | `Klhl21`\*, `Mapkapk5`\*, `Reep5`\*, `Itgb5`\* |
| Skin | `Plscr1` | None |
| Spleen | `Rai14`, `Myl9`\*, `Ptprk` | None |
| Thymus | `Hsd17b11`\*, `Etv1` | `Nusap1`, `Stmn1`, `Birc5`, `Cdk1`, `Top2a`, `Ccnb2`, `Aurka`, `Ccne2`, `Kif20a`, `Pcna`, `Ccnf` |
| Gastrocnemius | `Nfkbia` | `Fhl2` |
| Tibialis anterior | `Cebpd` | None |

**Reinforced genes**

| Analysis unit | FLT higher | FLT lower |
|---|---|---|
| Adrenal gland | None | `Tspan4` |
| Eye | None | `Klhl21` |
| Kidney | `Slc37a4` | None |
| Skeletal muscle, pooled | `Sox4`, `Cebpd`\*, `Sh3bp5`, `Prkcd`, `Arid5b`\*, `Sesn1`\*, `Tle1`\* | `Bphl`\* |
| Spleen | `Loxl1` | None |
| Thymus | `Snx7` | `Ube2c`, `Gmnn` |
| Soleus | `Tpm1` | `Bdh1`, `Ech1`, `Bnip3`, `Decr1` |
| Tibialis anterior | `Cdkn1a`, `St3gal5`, `Bnip3` | None |

Heart, liver, retina, EDL, and quadriceps had real-data BH-FDR genes but no
synthetic-informed gene. Bone, bone marrow, brain, brown adipose
tissue, cecum, cerebellum, colon, hippocampus, lung, mammary gland, optic nerve,
and white adipose tissue had no BH-FDR gene in the 974-gene panel.

The pooled skeletal-muscle results remain auditable in Tables S11-S12 and are
interpreted separately from anatomical groups because those groups have
different responses. Liver has 19 real-data BH-FDR associations but selected a
real-only arm. Skin has three real-data BH-FDR associations and one
synthetic-promoted gene, `Plscr1`. Lung has no BH-FDR gene in the 974-gene
panel.

## S13. Supplementary figures

![Muscle arm heatmap.](figures/figure_s1_muscle_arm_heatmap.png)

<p class="caption"><strong>Figure S1. Repeated nested muscle-group balanced accuracy.</strong> Each row is a muscle group and each column is a downstream use of real or generated profiles. Arm selection also required nonworse AUROC and average precision.</p>

![Downstream utility.](figures/figure_s2_downstream_utility.png)

<p class="caption"><strong>Figure S2. Downstream utility of generated expression.</strong> (A) Direct pooled augmentation on the locked real-profile test. (B) Selected-arm changes in balanced accuracy, AUROC, and average precision across repeated tissue-specific development splits. All evaluations used real profiles.</p>

## S14. Supplementary data tables

- `table_s1_archs4_ddim_metrics.tsv`
- `table_s2_locked_ddim_repeats.tsv`
- `table_s3_naive_augmentation.tsv`
- `table_s4_thymus_core_genes.tsv`
- `table_s5_thymus_reactome.tsv`
- `table_s6_muscle_group_summary.tsv`
- `table_s7_soleus_genes.tsv`
- `table_s8_muscle_reactome.tsv`
- `table_s9_all_tissue_development_screen.tsv`
- `table_s10_synthetic_informed_bh_fdr_genes.tsv`
- `table_s11_all_random_effects_bh_fdr_genes.tsv`
- `table_s12_bh_fdr_tissue_summary.tsv`
- `table_s13_liver_harmonization_benchmark.tsv`
- `table_s14_liver_harmonization_full_metrics.tsv`
- `table_s15_wgan_validation_repeats.tsv`
