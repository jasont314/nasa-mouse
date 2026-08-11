# OSDR Conditional Diffusion Finalist

## Status

The fixed ARCHS4-pretrained, OSDR-adapted ModelDDIM is the first OSDR conditional
generator in this benchmark to pass the predeclared broad locked-test rule. This is
an accepted **within-study interpolation** model, not an unseen-study generator and
not a claim that synthetic augmentation improves FLT/GC classification.

The model and locked-test protocol were fixed in commit `7e8dec7` before the test
split was opened. The test was evaluated once with four predeclared generation
seeds. No composite score was used.

## Model And Data

- NASA input: API-derived OSDR mouse bulk RNA-seq counts; the old integrated raw
  OSDR H5 was not used.
- Reference model: the exact 227-million-parameter Lacan et al. ModelDDIM trained
  for 15,000 epochs on the 20-tissue mouse ARCHS4 reference.
- OSDR adaptation: factorized tissue, FLT/GC, study, and material conditioning with
  a rank-512 domain LoRA adapter, followed by a short correlation-regularized
  refinement.
- Representation: full-transcriptome TPM followed by the ARCHS4-train MaxAbs scale
  on the 974-gene mouse landmark panel.
- OSDR split: 781 train, 536 validation, and 293 locked-test profiles from 75
  accessions. Samples are split within accession/tissue/condition strata, so test
  studies are represented during training.
- Calibration: train-only global plus shrunk accession/tissue mean alignment and
  half-strength positive missing-covariance residual noise. FLT/GC labels fit no
  mean or covariance parameter; they are used only to balance paired residual noise.
  Final scaled expression is clipped to nonnegative values.

## Locked-Test Results

Four independent synthetic cohorts were generated for the same 293 real test
profiles. Every fidelity metric passed its finite-sample rule in all four cohorts.

| metric | mean | range | repeats passing |
|---|---:|---:|---:|
| Gene-correlation agreement | 0.977 | 0.974-0.979 | 4/4 |
| Precision | 0.998 | 0.997-1.000 | 4/4 |
| Recall | 0.997 | 0.997-0.997 | 4/4 |
| F1 | 0.997 | 0.997-0.998 | 4/4 |
| Adversarial accuracy | 0.458 | 0.454-0.464 | 4/4 |
| FD / real-split P95 | 0.075 | 0.047-0.089 | 4/4 |

The test-size real-bootstrap Corr floor was 0.950. Corr therefore passes the
finite-sample gate, but all four generations remain below the separate strict 0.98
paper benchmark. Precision, recall, F1, adversarial indistinguishability, FD,
diversity, memorization, and nonnegativity all pass independently.

Pooled FLT-minus-GC effect recovery passed three of four generations, exactly the
predeclared 75% stability threshold. Mean expression-effect correlation was 0.598
and mean direction agreement was 0.683. The failed seed still had effect correlation
0.453, but its 0.505 direction agreement missed the 0.55 requirement.

The skeletal-muscle accession-aware diagnostic passed all four locked-test
generations: mean meta-effect correlation was 0.608 and mean direction agreement was
0.606 across five accessions. This is stronger than the unstable validation result,
but it is not exact gene-level replication. The number of genes jointly passing the
real and synthetic LOO-FDR rule was 0, 0, 0, and 1 across the four generations.

## Downstream Utility

FLT/GC classifiers evaluated on the real locked test gave:

| training data | balanced accuracy | ROC AUC |
|---|---:|---:|
| Real OSDR train | 0.754 | 0.819 |
| Synthetic train | 0.700 | 0.751 |
| Real plus synthetic | 0.734 | 0.801 |

Synthetic samples contain useful FLT/GC information, but adding them did not improve
the real-only classifier. The current model is therefore accepted for conditional
generation and controlled simulation, not as validated training-data augmentation.

## Decision

- **Accepted:** broad tissue-conditioned FLT/GC generation for covariate profiles
  represented by the OSDR training studies.
- **Not accepted:** unseen-study transfer, strict paper-level Corr, exact
  LOO-stable gene recovery, or improved FLT/GC classifier augmentation.
- **Rejected alternatives:** the WGAN-GP arms, OSDR-only DDIM, unseen-accession
  transfer DDIM, and small skeletal-muscle specialist models all failed one or more
  independent fidelity or stability gates.

Machine-readable results are under
`outputs/generative/benchmark/runs/lacan_diffusion/osdr_factorized_study_lora512_correlation_refine_seed2020/evaluation/final_locked_test/`.
The fixed configuration is
`configs/generative/diffusion/osdr_factorized_study_lora512_correlation_refine.yaml`.
