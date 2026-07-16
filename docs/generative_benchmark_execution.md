# Adaptive Generative Benchmark Execution

This document defines how the three paper families are executed without conflating
paper reproduction, mouse data substitution, and NASA-specific extensions.

## Paper Contracts

All upstreams are checked out under `assets/model_sources/` and verified by Git
commit plus SHA-256 hashes of architecture, loss, and training files. The immutable
contracts are implemented in `src/nasa_mouse_generative/paper_contracts.py`.

### Viñas conditional WGAN-GP

- Source commit: `94fa44dd1bd52d924efd3af0fcd8eeb18bd141a8`.
- Architecture: latent 64; categorical embedding size
  `int(sqrt(cardinality)) + 1`; one released-code constant numeric placeholder;
  two 256-unit ReLU layers; linear expression and critic outputs.
- Training: batch 32, five critic steps, gradient penalty 10, TensorFlow-equivalent
  RMSprop (`lr=5e-4`, `rho=0.9`, `epsilon=1e-7`), at most 2,000 epochs.
- Input transform: `log1p` applied directly to the supplied normalized expression,
  then training-gene standardization. This is the `model_native` contract. CPM
  normalization for depth-dependent NASA/ARCHS4 counts is tested separately as
  `wgan_nasa_cpm_zscore`; it is a declared mouse data-interface adaptation and is
  not scored as a strict preprocessing reproduction.

There is a real paper/source early-stopping discrepancy. The released zero-based loop
evaluates gamma after epochs 1, 6, 11, and so on, and uses ten unsuccessful checks,
equivalent to up to 50 epochs without improvement. The paper text says 30 epochs of patience.
Profiles `paper_native` and `paper_native_paper_text` expose these variants explicitly.

### Lacan DDIM

- Source commit: `cde890154698fcea96c924804aaff04af3351b48`.
- Architecture: official `ModelDDIM`, 974 landmarks, two 8,192-unit residual layers,
  tissue embedding dimension two, and dropout 0.1.
- Training: 1,000 quadratic diffusion steps, summed noise MSE, antithetic timesteps,
  Adam at `0.0004783833151836702`, exact OneCycle peak behavior, batch 2,048,
  15,000 epochs, AMP, parameter EMA 0.999, and no effective gradient clipping.
- Input transform: full-transcriptome TPM, landmark selection, then train-fitted
  MaxAbs scaling.

The exact implementation is `src/nasa_mouse_rna_diffusion/`, not the smaller generic
conditional adapter. A generic run labeled `paper_native` is rejected. The completed
mouse ARCHS4 run uses the manuscript 9,796/2,448/5,000 split. The released code's
12,244-profile train-plus-validation behavior remains a separate declared variant.
The local A100 run completed in 5,987 seconds; the paper reports approximately 3 h
7 min for its GTEx training environment.

The OSDR extension uses that same pinned `ModelDDIM`, summed noise objective,
OneCycle schedule, 15,000-epoch duration, AMP, and EMA. It substitutes the data and
joint conditioning class only. API-derived raw counts are converted to TPM using all
48,303 OSDR genes with positive GENCODE M39 lengths before selecting the 974 mouse
landmarks; MaxAbs is fitted on training accessions only. The default joint classes are
`tissue` plus `condition`, yielding 48 observed classes without held-out unknowns.
Study conditioning is not silently enabled because a held-out accession is an unseen
study class for this upstream joint-class architecture.

An epoch is not a comparable unit across the reference and query cohorts. The
ARCHS4 paper-parity training partition has five batches per epoch and therefore
executes about 75,000 optimizer steps. The OSDR partition fits in one batch and
executes 15,000 optimizer steps at the same nominal 15,000-epoch duration. This is
the paper-duration baseline; any step-matched extension must be named and reported
as a NASA duration adaptation.

The ARCHS4-transfer arms use the exact completed tissue model and the ARCHS4-fitted
MaxAbs scale. The `reference_only` baseline maps tissue parameters into
`tissue=<name>||condition=reference` classes and leaves FLT/GC condition columns at
their seeded initialization. Because the upstream embeds a one-hot condition before
the residual condition layers, this mapping is not function preserving. The
`function_preserving_tissue` strategy instead copies the upstream-used embedding
rows and every shared tissue condition column, adjusts each condition-layer bias,
and zero-initializes unmatched tissue columns. Numerical tests verify identical
source and expanded denoiser outputs before fine-tuning for every shared tissue and
expanded condition slot.

All parameters are then fine-tuned on accession-grouped OSDR training data. OSDR
values under the reference scale have median absolute value 0.00856, 95th percentile
0.0726, and only 0.007% of entries above one, so the transfer input scale is not
grossly out of distribution. AMP overflow handling advances the optimizer, EMA,
scheduler, and global step only when `GradScaler` actually executes an optimizer
step; skipped steps are audited in the history and run summary. Optional DDIM eta is
seeded deterministically and reported as a named evaluation variant; the released
default remains eta zero.

```bash
PYTHONPATH=src python -m nasa_mouse_rna_diffusion prepare-osdr
PYTHONPATH=src python -m nasa_mouse_rna_diffusion train-osdr
PYTHONPATH=src python -m nasa_mouse_rna_diffusion evaluate-osdr
```

Validation uses 12 unseen accessions. The 12-accession test partition remains locked;
both the configuration and `--unlock-test` must opt in after model selection.

### GeneJEPA

- Source commit: `a2f4d7218b17f2f52cc5f1cc94420c8ef1ae3265`.
- The official Perceiver, Fourier/value tokenizer, random block masking,
  student/EMA-teacher, predictor, cosine objective, and anti-collapse losses are used.
- Architecture: `d=768`, 512 latents, 24 blocks, 12 heads, mask ratio 0.45.
- Training: AdamW `1e-4`, weight decay `2e-4`, 5% linear warmup plus cosine decay,
  batch 92/device, accumulation two, bfloat16, gradient clipping one, and 50 epochs.
- Released data duration: one million sampled training profiles per epoch, or about
  50 million profile exposures. The paper uses four H100 80 GB GPUs and does not
  report a directly portable wall-clock time.
- Input transform: `log1p` nonzero expression and one train-fitted scalar mean/SD;
  zero remains the absent-gene token.

GeneJEPA is representation-only. Any GeneJEPA-guided diffusion is a new experimental
method and must be compared against the unguided exact DDIM.

The exact released architecture fits one A100 40 GB. Batch 92 with accumulation two
peaked at 29.88 GB and processed 184 profile exposures in 6.57 seconds. Extrapolating
that measured throughput to the released 50-million-exposure duration is roughly 20
GPU-days on this single device. Therefore a shorter run can test representation
quality, but must be labeled as a duration adaptation and cannot be described as a
paper-duration reproduction.

## Execution

Generate the concrete plan:

```bash
PYTHONPATH=src python -m nasa_mouse_generative experiment-plan
```

Preview executable and unresolved rows:

```bash
PYTHONPATH=src python -m nasa_mouse_generative matrix-run --dry-run --max-runs 0
```

Run at most four rows by default, persisting status after every transition:

```bash
PYTHONPATH=src python -m nasa_mouse_generative matrix-run \
  --phase phase_0_smoke --smoke
```

`experiment_status.tsv` records `planned`, `running`, `complete`, `failed`,
`awaiting_selection`, and `awaiting_accession_selection` states. Completed rows are
not rerun. Use `--retry-failed` explicitly. Single/selected-accession rows do not run
until accession lists are supplied. Later phases containing `best_*` placeholders do
not run until a prior phase has selected concrete values.

`training.repeats` expands deterministic consecutive seeds. Configured generated
samples and synthetic-to-real ratios are evaluated independently. Real-plus-synthetic
augmentation is calculated only after fidelity, diversity, and memorization gates
pass.

## Storage And Checkpoints

New unified runs default to one atomic `latest.pt` checkpoint every 100 epochs.
Only the latest checkpoint is overwritten; no epoch archive is accumulated. After a
successful final model and validation evaluation, the training checkpoint is removed
unless `execution.retain_training_checkpoint=true`. Prepared matrices are regenerated
deterministically by default rather than duplicated in every run directory; set
`execution.save_prepared_data=true` for a selected archival finalist.

Every run enforces minimum-free-space and maximum-run-size guards and records disk
state, runtime by stage, device, CUDA peak allocation, source hashes, selected profile
and accession hashes, split identity, preprocessing fit scope, and inverse-transform
policy. The completed exact DDIM output is preserved unchanged.

## Adaptive Decisions

Build the scoreboard after each bounded batch:

```bash
PYTHONPATH=src python -m nasa_mouse_generative scoreboard
```

Promotion is based on held-out fidelity subject to fixed gates selected before the
conditional results are inspected. Fidelity requires composite >= 0.70, gene-mean
correlation >= 0.80, gene-SD correlation >= 0.70, precision/recall F1 >= 0.50, and
adversarial indistinguishability >= 0.50. Diversity requires recall >= 0.10 and a
synthetic-to-real global-SD ratio >= 0.10; memorization requires no more than 5% of
synthetic samples below the first percentile of real training nearest-neighbor
distances. Conditional-effect eligibility additionally requires FLT-minus-GC delta
correlation >= 0.30 and direction agreement >= 0.55. FLT/GC classifier augmentation
is not scored as eligible unless all four gates pass.

The locked OSDR test accessions remain unavailable until preprocessing, architecture,
and seed policy are fixed. A failing screen should redirect the next run; it should
not trigger the full Cartesian matrix.

The exact ARCHS4 DDIM is currently the only full generator that passes the broad
tissue benchmark. The one-epoch study-conditioned WGAN orchestration check correctly
fails diversity and is not a biological result.

Neither direct nor ARCHS4-transferred conditional DDIM currently passes all OSDR
gates. The best function-preserving validation checkpoint passes pooled FLT/GC effect,
diversity, and memorization checks, but fails the fidelity composite and local
real-versus-synthetic indistinguishability checks. Consequently no conditional model
has been promoted to three seeds, augmentation claims, or the locked test.
