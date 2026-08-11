# NASA Mouse Generative Model Benchmark

Completed benchmark decisions and final metrics are summarized in
[`generative_benchmark_results.md`](generative_benchmark_results.md).

## Objective

The primary objective is realistic synthetic mouse bulk RNA-seq expression. A
secondary objective is determining whether synthetic data improves FLT versus GC
classification on real, unseen OSDR accessions. Conditional FLT/GC generation is
in scope for version 1. Paired counterfactual generation is deferred until a model
passes fidelity, diversity, and memorization gates.

OSDR data must come from the NASA OSDR Biological Data API. The older integrated
OSDR H5 is not an allowed input.

## Data Cohorts

The OSDR inventory includes every Mus musculus bulk RNA-seq profile returned for
Space Flight or Ground Control. Material labels are canonicalized while retaining
the original API values and mapping audit.

Per-tissue eligibility has three tiers:

1. `confirmatory_per_tissue`: at least 60 samples, 20 per condition, and five
   accessions containing both FLT and GC. This supports separate validation and
   locked test accessions.
2. `exploratory_pretrained_per_tissue`: at least 30 samples, 10 per condition,
   and two paired accessions. These runs use LOO-style reporting and are not called
   independently confirmed.
3. `pooled_only`: all remaining materials participate in the tissue-conditioned
   pooled model but do not receive a standalone model.

The entire ARCHS4 mouse H5 metadata catalog is scanned. Primary pretraining uses a
`healthy_preferred` cohort: bulk-like RNA-seq profiles matching an OSDR tissue,
excluding explicit disease, tumor, genetic perturbation, intervention,
developmental-mismatch, and spaceflight-leakage terms. Metadata without a clear
health label is retained and explicitly marked unknown because requiring the word
"healthy" would discard most usable GEO profiles.

A `broad` sensitivity cohort retains disease and perturbation diversity while still
excluding single-cell-like and spaceflight-leakage profiles. It is not the primary
cohort because its condition structure can conflict with healthy OSDR animals.

Reference balancing caps each tissue and GEO series and writes hierarchical sampling
weights. Training samples tissue uniformly, then study uniformly within tissue, then
sample uniformly within study. Rare tissues are retained rather than downsampling all
tissues to the smallest class.

Unreadable compressed ARCHS4 columns are excluded rather than zero-imputed. Their H5
sample indices are stored in the content-addressed cache and run manifest. A run stops
if the count exceeds `data.archs4_max_corrupt_profiles` (100 by default), preventing a
damaged source file from being silently treated as a usable cohort.

The completed scan selected 62,299 healthy-preferred profiles from 5,307 GEO
series across 23 matchable classes. A stricter 23,614-profile control-only cohort
and a 134,250-profile broad cohort are retained as sensitivity arms. See
`docs/generative_data_audit.md` for selection details and the comparison with the
older 24,428-profile reference.

## Configurable Axes

- model: Vinas conditional WGAN-GP or Lacan landmark-space diffusion;
- preprocessing: raw counts, CPM, TPM when mouse gene lengths are available,
  `log1p`, `log2(x+1)`, gene z-score, robust scaling, or MaxAbs scaling;
- feature space: all 48,694 shared genes, fold-selected HVGs, Reactome genes, or
  mapped mouse L1000 landmarks plus reconstruction targets;
- cohort: one accession, selected accessions, or all eligible accessions;
- harmonization: none, within-study z-score, within-study then pooled z-score,
  ComBat, ComBat-seq, MBatch Median Polish, MBatch Empirical Bayes, MBatch
  ANOVA, or MOBER;
- training: OSDR only, ARCHS4 only, or ARCHS4 pretrain then OSDR fine-tune;
- tissue mode: pooled tissue-conditioned or standalone per tissue;
- condition mode: FLT/GC conditioned or unconditional negative control;
- study policy: no study input or explicit study conditioning;
- balancing, seed, executable repeats, generated sample count, and synthetic-to-real ratio;
- technical-replicate handling (`keep`, `sum`, or `mean`, with `sum` as default);
- optional conditioning on tissue, material type, muscle group, study, sex, assay,
  platform, and data source.

Normalization, transformation, scaling, and harmonization are separate parameters.
All fitted statistics use training folds only. Study-wise scaling for an unseen study
uses training-global fallback by default. Estimating statistics from the held-out
study is an explicitly labeled transductive sensitivity analysis.

HVG ranking can be bounded independently with
`features.selection_sample_limit`. The default deterministically round-robins 5,000
training profiles across tissues for ranking, but does not reduce the profiles used
to train the model. Set it to `0` for an exact all-training-profile variance scan.

## Harmonization Arms

Primary implementation references are the
[Scanpy ComBat documentation](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.combat.html),
[Bioconductor `sva` package](https://bioconductor.org/packages/release/bioc/html/sva.html),
[official MBatch repository](https://github.com/MD-Anderson-Bioinformatics/MBatch),
and [official MOBER repository](https://github.com/Novartis/MOBER).

- `within_study_zscore` implements the study-by-study log and gene scaling design
  used by Ilangovan et al. 2024.
- `within_study_then_global_zscore` adds the mentor-proposed second scaling after
  concatenation.
- `combat`, `combat_seq`, `mbatch_median_polish`,
  `mbatch_empirical_bayes`, and `mbatch_anova` reproduce the five correction
  families assessed by Sanders et al. 2023. Their finding that
  library-preparation/ComBat ranked best was specific to seven liver datasets
  and is not treated as universal.
- `mober` uses a batch-aware VAE and projection to a trained source. It must be
  checked for preservation of FLT/GC effects, not merely batch mixing.

ComBat, ComBat-seq, and the three MBatch methods have no natural inductive
transform for a new batch. Any use on a held-out accession is marked transductive.
MOBER can project new samples onto a trained source and is evaluated as an
inductive model-based harmonizer.

All six model/correction adapters are executable and serialized with each run:

| Method | Fit and held-out behavior | Main controls |
| --- | --- | --- |
| `combat` | Parametric empirical-Bayes fit; frozen correction for known batches and unlabeled transductive estimation for a new batch | `batch_key`, `max_batches`, `confounded_covariate_policy` |
| `combat_seq` | Bioconductor `sva::ComBat_seq`; each held-out partition is corrected with training anchors | `rscript`, `anchor_samples`, `noninteger_policy`, `singleton_batch_policy`, `confounded_covariate_policy` |
| `mbatch_median_polish` | Official MBatch batch-wise Median Polish; held-out partition corrected with training anchors | `rscript`, `source_root`, `anchor_samples` |
| `mbatch_empirical_bayes` | Official MBatch parametric Empirical Bayes; held-out partition corrected with training anchors | `rscript`, `source_root`, `anchor_samples` |
| `mbatch_anova` | Official variance-adjusted MBatch ANOVA; held-out partition corrected with training anchors | `rscript`, `source_root`, `anchor_samples`, `nonfinite_policy` |
| `mober` | Batch-aware adversarial VAE fit; frozen encoder mean decoded onto one trained target batch | `target_batch`, `encoding_dim`, `epochs`, `batch_size`, `learning_rate`, `adversary_weight`, `kl_weight` |

`batch_key=auto` uses `source` when ARCHS4 and OSDR are jointly fitted, and
`study` for OSDR-only runs. In the pretrain/fine-tune regime, complex harmonizers
fit on ARCHS4 plus the OSDR training partition only. Validation and test samples
never enter that fit. MOBER automatically targets the OSDR source in this regime.

ComBat and MBatch variants require
`validation.allow_transductive_preprocessing=true`; this is a deliberate opt-in,
not a leakage-free inductive benchmark. Their preservation covariates default to
`condition`, `tissue`, and `sex`. A preservation variable that is rank-confounded
with batch fails by default. `confounded_covariate_policy=drop` is available only
as an audited sensitivity arm. Preserving `condition` also makes preprocessing
outcome-informed, so FLT/GC classifier metrics from that arm are labeled unsuitable
as blind prediction estimates. Run a complementary outcome-blind arm with
`preprocessing.harmonization_covariates=[tissue,sex]`.

ComBat-seq is stricter because it operates on counts. NASA API-derived count-like
profiles include fractional values, so the default is to stop. An exploratory run
may set `noninteger_policy=round`; the fraction and maximum distance rounded are
recorded per partition. One-sample batches also stop by default. Set
`singleton_batch_policy=identity` to leave those samples uncorrected or `pool` for
an explicitly exploratory pooled-singleton arm. The R and `sva` versions, anchor
size, rounding, singleton handling, and confounded covariates are written to
`harmonizer.json` and evaluation summaries.

The MBatch adapters invoke the correction functions from the pinned official R
source checkout rather than reimplementing their equations. Median Polish,
Empirical Bayes, and ANOVA use balanced training anchors when correcting a held-out
partition. Variance-adjusted ANOVA can return non-finite values for genes with
degenerate within-batch variance. The optional `nonfinite_policy=identity_gene`
restores each affected gene in full and records the count; it never fills isolated
cells silently.

MOBER follows the published encoder/decoder and adversarial batch-classifier design
in a self-contained PyTorch adapter. Projection is deterministic: the encoder mean
is decoded onto the selected target batch. MOBER does not consume biological
preservation covariates itself, so FLT/GC preservation must be checked empirically.
The implementation is architecture-compatible with the official repository but is
not a vendored invocation of its training script.

## Shared And Native Preprocessing

Each generator receives both:

1. a shared preprocessing benchmark for controlled comparison; and
2. its paper-native preprocessing, including Vinas log/z-score and Lacan
   TPM/L1000/MaxAbs where gene lengths and ortholog mappings permit it.

Lacan's native arm uses TPM followed by training-fold MaxAbs scaling. Mouse gene
lengths are generated from the union of GENCODE M39 exon intervals for each
versionless Ensembl gene ID. The table in
`data/reference/gencode_vM39_mouse_gene_lengths.tsv` covers all 980 mapped mouse
L1000 genes; its source URL and GTF SHA-256 are recorded in the adjacent manifest.

## Validation And Model Selection

No sample-random split is allowed when multiple accessions are pooled. Entire OSDR
accessions and ARCHS4 GEO series are grouped during splitting.

Checkpoint and hyperparameter selection report the following held-out-validation
metrics separately:

- gene mean, variance, zero fraction, and quantile agreement;
- gene-gene correlation and pathway-correlation agreement;
- precision, recall, and density/coverage in a train-fitted embedding;
- adversarial real-versus-synthetic accuracy;
- nearest-neighbor memorization and duplicate checks;
- tissue and condition consistency.

There is no composite score or compensating rank. For the Lacan-paper metric
protocol, correlation-matrix agreement must be at least 0.98, precision at least
0.95, recall at least 0.85, F1 at least 0.90, symmetric nearest-neighbor adversarial
accuracy between 0.40 and 0.60, and mouse PCA Frechet distance no greater than the
95th-percentile real-versus-real split reference. Every requirement must pass on
both the paper-style training-distribution comparison and the held-out-accession
comparison. The mouse PCA Frechet reference replaces the paper's pretrained human
GTEx classifier because that classifier is not transferable to mouse genes.

The absolute paper gate is always reported. Corr is also evaluated against a
finite-sample gate whose minimum is the lower of 0.98 and the fifth percentile of
same-size real-data bootstrap Corr values. This is necessary for small OSDR tissue
partitions where even two real draws cannot reach 0.98. It does not average or
compensate metrics: precision, recall, F1, AA, FD, diversity, and memorization retain
their original independent requirements, and passing the calibrated gate is never
described as matching the paper's absolute Corr benchmark.

Eligibility additionally requires a synthetic-to-real global standard-deviation
ratio between 0.5 and 2.0 and no more than 5% of synthetic samples closer to a
training profile than the training leave-one-out first percentile. FLT/GC pooled and
accession-aware effect recovery are separate mandatory gates for a conditional
generator. These thresholds are benchmark gates, not biological significance tests.

The final test accession remains locked until model and preprocessing choices are
fixed. Final reporting includes real-train/real-test, synthetic-train/real-test,
real-train/synthetic-test, and real-plus-synthetic augmentation performance. FLT/GC
gene effects, pathway effects, accession-aware random effects, and LOO direction/FDR
stability are reported without using the final test set for tuning.

The resolved split plan currently contains 51 pooled training, 12 validation, and
12 locked-test accessions. Six tissues meet the confirmatory standalone threshold;
14 tissues support 72 total leave-one-accession-out folds.

## Staged Execution

The benchmark is gated to avoid an uncontrolled Cartesian search:

1. smoke-test finite generation and inverse transforms;
2. screen shared and native preprocessing on pooled data;
3. test study-conditioning and harmonization finalists;
4. expand finalists to every tissue tier and five seeds;
5. run unconditional negative controls;
6. consider paired counterfactual generation only after fidelity gates pass.

Configurations live under `configs/generative/`. Data audit outputs and experiment
plans live under `outputs/generative/benchmark/`. Model-specific training output will
be written under a run ID derived from the complete resolved configuration.

On the A100 40 GB, one 225.6-million-parameter liver ModelDDIM records about 5.06 GB
peak PyTorch allocation and about 6.2 GiB process residency. Six concurrent models
fit in memory, but the GPU reaches 100% compute utilization and each model slows
substantially. Parallelism is therefore bounded by measured aggregate throughput,
not free memory alone. Distinct screening arms may run concurrently; repeat seeds
are deferred until an arm passes every single-seed validation gate.

## Implementation Status

Implemented in `src/nasa_mouse_generative/`:

- two-model capability and paper/code provenance registry;
- validated configuration schema and gated experiment plan;
- current OSDR API inventory, technical-replicate-aware expression builder, and
  tissue eligibility tiers;
- full ARCHS4 metadata scan and balanced reference manifests;
- accession-grouped locked and LOO split plans;
- fold-aware shared preprocessing and harmonization constraints;
- reloadable ComBat, R `sva::ComBat_seq`, three official MBatch, and MOBER
  harmonization adapters with fold behavior and fallback audits;
- executable Vinas WGAN-GP and Lacan diffusion adapters;
- full-H5 ARCHS4 extraction with content-addressed caching and hierarchical sampling;
- direct OSDR, ARCHS4-only, and ARCHS4-pretrain/OSDR-fine-tune stage orchestration;
- epoch checkpoints, deterministic resume, final model serialization, and GPU records;
- conditioned generation, inverse preprocessing, held-out fidelity/diversity/
  memorization metrics, and real/synthetic FLT/GC classifier evaluation;
- locked-test enforcement and pooled, individual-tissue, or all-eligible-tissue CLI runs;
- resumable experiment-matrix execution, bounded rows per invocation, durable status,
  disk guards, sparse checkpoint retention, and ranked scoreboard generation;
- ARCHS4-only GEO-series train/validation/test partitions and held-out tissue plots;
- Lacan-style DDIM `t=1000,200,0` PCA trajectory.

The configurable Lacan adapter was checked against official code commit
`cde890154698fcea96c924804aaff04af3351b48`, but it remains an extension for
arbitrary OSDR covariates and pretrain/fine-tune stages. Exact architecture and
training-procedure comparisons use the separate pinned implementation under
`src/nasa_mouse_diffusion/paper_parity/`; see `docs/rna_diffusion_paper_parity.md`.

The `paper_native` hyperparameter profiles enforce the source configurations and
reject locked-value overrides. Exact Lacan runs with 8,192-unit hidden layers,
batch 2,048, and 15,000 epochs are expensive. `practical_screen` retains the model family,
native transformation, conditioning, objective, diffusion schedule, and held-out
evaluation while reducing architecture or epochs. Outputs must state which profile
was used; a practical bulk adaptation is not an exact paper reproduction.

## ARCHS4 Tissue Benchmark Results

The first ARCHS4-only benchmark used 62,297 readable healthy-preferred profiles
from 5,307 GEO series and 23 canonical tissue classes. GEO series were assigned as
whole groups: 43,744 profiles for training, 10,003 for validation, and 8,550 for the
locked test. Two unreadable source columns were excluded and recorded. These runs
did not open or fine-tune on OSDR expression. All training and figure inference ran
on an NVIDIA A100-SXM4-40GB.

The earlier reduced Lacan proxy and its failed 100/500-epoch output directories were
removed. They changed architecture, optimizer behavior, sampling weights, training
duration, input normalization, and split policy simultaneously, so they could not
attribute a result to mouse data. The replacement uses the unmodified upstream
227,109,786-parameter `ModelDDIM`, the paper's 9,796/2,448/5,000 split sizes, and
full-transcriptome TPM before selecting 974 mouse landmark genes. Its run directory
is `outputs/generative/benchmark/runs/lacan_diffusion/archs4_mouse_paper_parity_seed1234/`.
The completed 15,000-epoch run reached 0.869 synthetic-to-real held-out tissue
balanced accuracy versus 0.895 for real-to-real, direct L974 precision/recall of
0.966/0.865, gene mean/SD/correlation-matrix agreement of 0.997/0.944/0.879, and
nearest-neighbor adversarial accuracy of 0.512. Thus the exact model learns useful
mouse tissue-conditioned generation, unlike the deleted proxy. The first two PCA
components still overlap across tissues (silhouette -0.271), and 9.04% of generated
scaled entries are negative; see `docs/rna_diffusion_paper_parity.md` for the full
interpretation and output paths.

## Installation

```bash
conda activate nasa-mouse
python -m pip install -r requirements-nasa-mouse-generative.txt
conda install -c conda-forge -c bioconda r-base=4.5 bioconductor-sva=3.58
PYTHONPATH=src python -m nasa_mouse_generative prepare-upstreams
```

The upstream command checks out and verifies the pinned WGAN-GP and diffusion
sources under ignored `assets/model_sources/`.

## Training Commands

Run the bounded end-to-end WGAN check first:

```bash
PYTHONPATH=src python -m nasa_mouse_generative train \
  --config configs/generative/benchmark/default.yaml \
  --set training.regime=osdr_only \
  --set execution.device=cuda \
  --smoke
```

Run ARCHS4 pretraining followed by OSDR fine-tuning with explicit parameter
overrides:

```bash
PYTHONPATH=src python -m nasa_mouse_generative train \
  --set training.model=vinas_wgan_gp \
  --set training.regime=archs4_pretrain_osdr_finetune \
  --set features.space=hvg \
  --set features.hvg_genes=2000 \
  --set training.model_parameters.reference_epochs=100 \
  --set training.model_parameters.finetune_epochs=50 \
  --set execution.device=cuda
```

Enable the inductive MOBER adapter on the same regime:

```bash
PYTHONPATH=src python -m nasa_mouse_generative train \
  --set training.regime=archs4_pretrain_osdr_finetune \
  --set preprocessing.harmonization=mober \
  --set preprocessing.harmonization_parameters.epochs=300 \
  --set execution.device=cuda
```

Run ComBat as a transductive sensitivity arm:

```bash
PYTHONPATH=src python -m nasa_mouse_generative train \
  --set training.regime=archs4_pretrain_osdr_finetune \
  --set preprocessing.harmonization=combat \
  --set 'preprocessing.harmonization_covariates=[tissue,sex]' \
  --set validation.allow_transductive_preprocessing=true \
  --set execution.device=cuda
```

Run ComBat-seq with explicit policies required by the current API-derived matrix:

```bash
PYTHONPATH=src python -m nasa_mouse_generative train \
  --set training.regime=osdr_only \
  --set preprocessing.harmonization=combat_seq \
  --set 'preprocessing.harmonization_covariates=[condition,sex]' \
  --set preprocessing.harmonization_parameters.noninteger_policy=round \
  --set preprocessing.harmonization_parameters.singleton_batch_policy=identity \
  --set validation.allow_transductive_preprocessing=true \
  --set execution.device=cuda
```

Use `training.model=lacan_diffusion` for Lacan. Named preprocessing profiles are
selected with `--set preprocessing.profile=shared_log1p_cpm_zscore` or
`--set preprocessing.profile=model_native`. Lacan's native TPM profile requires an
aligned mouse gene-length table in `preprocessing.gene_lengths`.

Run one tissue or every confirmatory/exploratory tissue:

```bash
PYTHONPATH=src python -m nasa_mouse_generative train \
  --set training.tissue_mode=per_tissue --tissue liver

PYTHONPATH=src python -m nasa_mouse_generative train --all-tissues
```

An ARCHS4-only generator must set `training.condition_on_flight=false`: ARCHS4 has
no flight labels, so it is only a tissue/reference baseline and cannot identify a
spaceflight effect.

Prepare and train the paper-parity tissue-conditioned ARCHS4 diffusion baseline:

```bash
PYTHONPATH=src python -m nasa_mouse_diffusion.paper_parity prepare
PYTHONPATH=src python -m nasa_mouse_diffusion.paper_parity train
PYTHONPATH=src python -m nasa_mouse_diffusion.paper_parity evaluate
```

Selected archival runs can save `prepared_data.h5`; loading remains compatible with
the older `prepared_osdr.h5` filename. The storage-safe default records profile and
accession hashes and deterministically rebuilds matrices for later evaluation. In an
ARCHS4-only run, `reference` is the ARCHS4 training partition and no OSDR expression
file is opened.

Create paper-style held-out tissue figures after the configuration is fixed:

```bash
PYTHONPATH=src python -m nasa_mouse_diffusion.paper_parity evaluate
```

The paper-parity diffusion evaluation fits one PCA basis only on
real ARCHS4 training profiles, then projects real profiles and exact 1,000-step EMA
DDIM trajectory snapshots. Tissue classifier and silhouette metrics are written
beside the figures.

One-study runs use `data.osdr_accession_scope=single` with exactly one
`data.osdr_include_accessions` value. Because accession holdout is impossible in
that design, the run manifest labels its deterministic condition-stratified sample
split. Multi-study runs retain accession-grouped validation whenever possible.

## Evaluation And Generation

Training automatically evaluates held-out validation accessions. Representation
metrics neutralize the supplied condition token before FLT/GC classification so the
model cannot score by reading its conditioning label. The final test remains locked:

Generator promotion requires fidelity, diversity, and memorization gates plus two
condition checks: pooled FLT/GC delta recovery and accession-aware random-effects
recovery. The accession gate requires at least two eligible accessions, meta-effect
correlation of at least 0.30, and direction agreement of at least 0.55. Augmentation
is not evaluated before all gates pass. Per-tissue fidelity and FLT/GC recovery are
saved beside the global summary to expose composition-confounded pooled scores.

```bash
PYTHONPATH=src python -m nasa_mouse_generative evaluate \
  --run-dir outputs/generative/benchmark/runs/vinas_wgan_gp/RUN_ID \
  --split test --unlock-test
```

Generate FLT samples from an observed joint covariate profile, overriding fields as
needed:

```bash
PYTHONPATH=src python -m nasa_mouse_generative generate \
  --run-dir outputs/generative/benchmark/runs/vinas_wgan_gp/RUN_ID \
  --condition flight --set tissue=liver --n 100
```

Generation rejects fields that were not model-conditioning inputs. `study` is the
only additional selector: when study conditioning is disabled, it selects a coherent
observed profile and supplies the study identifier only to inverse preprocessing. The
generation manifest records model inputs and the full generation profile separately,
and flags requested combinations not observed during training.

Each run contains the resolved configuration, fitted preprocessing and categorical
vocabularies, data/split identities, epoch history, final model, validation metrics,
PCA/fidelity plots, embeddings, runtime/device/storage records, and a concise README.
The latest training checkpoint is overwritten sparsely and removed after successful
finalization unless retention is explicitly enabled.

ComBat, ComBat-seq, all three official MBatch methods, and MOBER have executable,
reloadable adapters. ComBat-seq uses R 4.5 and Bioconductor `sva`; MBatch invokes the
pinned official R functions; ComBat and MOBER run in the Python environment. The
one-epoch outputs validate mechanics only. The matched 15,000-epoch liver comparison
and its independent-metric plots are under
`outputs/generative/benchmark/summary/liver_harmonization/`.
