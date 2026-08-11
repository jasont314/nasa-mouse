# Study-conditioned OSDR WGAN

## Scope

This experiment trains the Viñas et al. conditional WGAN-GP as the direct
counterpart to the accepted study-conditioned DDIM experiment. It uses one
pooled model across tissues and studies, with categorical conditioning on:

- tissue;
- flight or ground-control condition;
- OSDR study/accession;
- material type;
- sex; and
- skeletal-muscle group.

The model therefore accepts a selected covariate profile and generates a
synthetic FLT or GC expression profile. It is a conditional generator, not a
paired FLT-to-GC counterfactual model.

## Data and split

The run reuses the fixed DDIM comparison data and gene order:

- 974 genes;
- 9,796 ARCHS4 reference-training profiles;
- 2,448 ARCHS4 reference-validation profiles;
- 781 OSDR training profiles;
- 536 OSDR validation profiles; and
- 293 locked OSDR test profiles, which were not loaded for WGAN selection.

OSDR expression and metadata come from the repository's API-derived prepared
data, not the deprecated raw integrated OSDR H5. The split is within
accession, tissue, and condition. Results therefore measure interpolation for
known studies, not transfer to an unseen study.

The paper transform is full-transcriptome TPM, `log1p`, and gene-wise z-score
using ARCHS4 training statistics. Evaluation is in the common nonnegative
MaxAbs expression space used for cross-model comparison.

## Architecture and training

The baseline uses the released paper topology and defaults: 64-dimensional
noise, two 256-unit ReLU layers, learned categorical embeddings, RMSprop at
`5e-4`, five critic steps, gradient penalty 10, batch size 32, and released-code
gamma early stopping. ARCHS4 pretraining stopped at epoch 161. OSDR fine-tuning
stopped at epoch 311. Training ran on an NVIDIA A100-SXM4-40GB.

Three OSDR refinement arms started from the same archived ARCHS4 checkpoint:

| Arm | OSDR LR | Augmentation | Stop epoch | Best native gamma |
|---|---:|---|---:|---:|
| no augmentation | `1e-4` | none | 401 | 0.9901 |
| mild augmentation | `1e-4` | `p=0.25`, SD 0.1 | 1,056 | 0.9944 |
| stronger augmentation | `1e-4` | `p=0.50`, SD 0.25 | 1,006 | 0.9935 |

The augmentation is the Gaussian profile augmentation exposed by the released
WGAN implementation. Refinement lineage is recorded in each run's
`reference_initialization.json`.

## Validation results

The strongest common-space result remains the original baseline with
train-only residual calibration, not the higher-gamma refinement:

| Metric | Mean over 6 seeds | Required | Repeat pass |
|---|---:|---:|---:|
| correlation agreement | 0.9759 | per-repeat >= 0.9756 | 4/6 |
| precision | 0.9764 | >= 0.95 | 6/6 |
| recall | 0.9938 | >= 0.85 | 6/6 |
| F1 | 0.9850 | >= 0.90 | 6/6 |
| adversarial accuracy | 0.6362 | 0.40-0.60 | 0/6 |
| FD / real-split p95 | 0.1439 | <= 1.0 | 6/6 |

The mild refinement passed precision, recall, F1, FD, memorization, and pooled
FLT/GC recovery in every repeat. Its calibrated means were correlation 0.9749,
precision 0.9798, recall 0.9941, F1 0.9869, adversarial accuracy 0.6391, and FD
ratio 0.1333. It did not improve the baseline's joint-distribution result.

All 12 predeclared train-only calibration variants failed independent model
selection. Increasing residual variance brought adversarial accuracy into its
target range, but lowered correlation below its finite-sample requirement.
Per-gene empirical quantile calibration also failed, indicating a joint
distribution mismatch rather than a marginal-range problem.

Pooled FLT/GC effect recovery was stable, but the skeletal-muscle
accession-aware diagnostic was not. The WGAN is therefore not approved for
synthetic augmentation or biological claims. Its locked test remains unopened.

## Interpretation

Study conditioning is implemented and the WGAN learns broad support and pooled
condition differences. It still leaves real and synthetic profiles separable:
raw generated profiles are too locally concentrated, while residual calibration
over-disperses them before correlation and AA can pass together. The accepted
study-conditioned DDIM remains the stronger generator because it passes
correlation, precision, recall, F1, AA, FD, condition recovery, and the muscle
diagnostic together on its fixed evaluation protocol.

## Commands

```bash
PYTHONPATH=src python -m nasa_mouse_wgan.matched_study train \
  --config configs/generative/wgan/wgan_matched_study_conditioned.yaml

PYTHONPATH=src python -m nasa_mouse_wgan.matched_study train \
  --config configs/generative/wgan/wgan_matched_study_refine_aug025.yaml

PYTHONPATH=src python -m nasa_mouse_wgan.matched_study evaluate-validation \
  --config configs/generative/wgan/wgan_matched_study_refine_aug025.yaml

PYTHONPATH=src python -m nasa_mouse_wgan.matched_study screen-calibration \
  --config configs/generative/wgan/wgan_matched_study_refine_aug025.yaml
```

Do not run `evaluate-test --unlock-test` unless a future fixed validation
candidate passes every independent gate.
