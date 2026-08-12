# Generative validation

## Scope

This is the current technical summary for model selection. It describes the
OSDR-disjoint DDIM branch used by the final manuscript and presentation. This
page omits earlier standalone runs and the pre-correction reference.

## Leakage and split controls

Nine GEO series containing 108 eligible OSDR-linked ARCHS4 profiles were removed
before reference training. ARCHS4 train, validation, and test partitions contain
complete, nonoverlapping GEO series. The final OSDR test is a within-study test:
held-out profiles come from represented accessions, but never enter training or
preprocessing fits. It measures interpolation within known studies rather than
transfer to a new mission.

## Generator comparison

| Model | Evaluation | Corr. | Precision | Recall | F1 | Adversarial accuracy | FD ratio | Decision |
|---|---|---:|---:|---:|---:|---:|---:|---|
| ARCHS4 DDIM | 4,628 held-out ARCHS4 profiles | 0.878 | 0.951 | 0.890 | 0.919 | 0.515 | 0.866 | Reference initialization; correlation agreement missed its target |
| Study-conditioned WGAN-GP | 536 OSDR validation profiles, six draws | 0.976 | 0.976 | 0.994 | 0.985 | 0.636 | 0.144 | Not used downstream because real and synthetic profiles remained separable and muscle effect recovery failed |
| Factorized DDIM | 293 OSDR test profiles, four draws | 0.974 | 0.997 | 0.996 | 0.997 | 0.475 | 0.074 | Selected for downstream analysis |

For the selected DDIM, all four draws passed the finite-sample fidelity,
diversity, and memorization checks. Three of four passed pooled FLT/GC effect
recovery, and all four passed the accession-aware skeletal-muscle effect check.
The WGAN locked test was not opened because its validation criteria were not met.

No scalar composite score was used. Correlation, precision, recall, F1,
real-versus-synthetic adversarial accuracy, distributional distance, diversity,
memorization, and condition-effect recovery were checked separately.

## Preprocessing decision

Nine matched liver DDIM runs compared no correction, within-study scaling,
two-stage scaling, ComBat, ComBat-seq, three MBatch methods, and MOBER. No
harmonization arm passed every fidelity and FLT/GC criterion. The selected model
therefore uses train-fitted TPM/MaxAbs preprocessing and explicit study
conditioning without global batch correction.

## Downstream use

The DDIM was used in tissue-specific classifier and feature analyses. Across 27
tissues and anatomical muscle groups, 18 real-plus-synthetic classifier arms were
nonworse on all six pooled and accession-macro metrics, and 16 improved at least
one metric. Gene and pathway associations were then tested on real OSDR profiles;
synthetic profiles affected training or ranking but were never counted as
additional animals.

The publication-facing values and interpretations are maintained in:

- [`table_4_generator_model_selection.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_4_generator_model_selection.tsv)
- [`table_s1_archs4_ddim_metrics.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_s1_archs4_ddim_metrics.tsv)
- [`table_s2_locked_ddim_repeats.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_s2_locked_ddim_repeats.tsv)
- [`table_s15_wgan_validation_repeats.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_s15_wgan_validation_repeats.tsv)
- [`paper/synthetic_guided_spaceflight/`](../paper/synthetic_guided_spaceflight/)

Exact configurations, run directories, and commands are linked from
[`generative_pipeline.md`](generative_pipeline.md) and
[`outputs/COMMANDS.md`](../outputs/COMMANDS.md).
