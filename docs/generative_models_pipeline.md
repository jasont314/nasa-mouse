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

- model: Vinas conditional WGAN-GP, Lacan landmark-space diffusion, or GeneJEPA;
- preprocessing: raw counts, CPM, TPM when mouse gene lengths are available,
  `log1p`, `log2(x+1)`, gene z-score, robust scaling, or MaxAbs scaling;
- feature space: all 48,694 shared genes, fold-selected HVGs, Reactome genes, or
  mapped mouse L1000 landmarks plus reconstruction targets;
- cohort: one accession, selected accessions, or all eligible accessions;
- harmonization: none, within-study z-score, within-study then pooled z-score,
  ComBat, ComBat-seq, or MOBER;
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
and [official MOBER repository](https://github.com/Novartis/MOBER).

- `within_study_zscore` implements the study-by-study log and gene scaling design
  used by Ilangovan et al. 2024.
- `within_study_then_global_zscore` adds the mentor-proposed second scaling after
  concatenation.
- `combat` and `combat_seq` follow the candidate methods assessed by Sanders et al.
  2023. Their finding that library-preparation/ComBat ranked best was specific to
  seven liver datasets and is not treated as universal.
- `mober` uses a batch-aware VAE and projection to a trained source. It must be
  checked for preservation of FLT/GC effects, not merely batch mixing.

ComBat and ComBat-seq have no natural inductive transform for a new batch. Any use
on a held-out accession is marked transductive. MOBER can project new samples onto a
trained source and is evaluated as an inductive model-based harmonizer.

All three adapters are executable and serialized with each run:

| Method | Fit and held-out behavior | Main controls |
| --- | --- | --- |
| `combat` | Parametric empirical-Bayes fit; frozen correction for known batches and unlabeled transductive estimation for a new batch | `batch_key`, `max_batches`, `confounded_covariate_policy` |
| `combat_seq` | Bioconductor `sva::ComBat_seq`; each held-out partition is corrected with training anchors | `rscript`, `anchor_samples`, `noninteger_policy`, `singleton_batch_policy`, `confounded_covariate_policy` |
| `mober` | Batch-aware adversarial VAE fit; frozen encoder mean decoded onto one trained target batch | `target_batch`, `encoding_dim`, `epochs`, `batch_size`, `learning_rate`, `adversary_weight`, `kl_weight` |

`batch_key=auto` uses `source` when ARCHS4 and OSDR are jointly fitted, and
`study` for OSDR-only runs. In the pretrain/fine-tune regime, complex harmonizers
fit on ARCHS4 plus the OSDR training partition only. Validation and test samples
never enter that fit. MOBER automatically targets the OSDR source in this regime.

ComBat variants require
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

MOBER follows the published encoder/decoder and adversarial batch-classifier design
in a self-contained PyTorch adapter. Projection is deterministic: the encoder mean
is decoded onto the selected target batch. MOBER does not consume biological
preservation covariates itself, so FLT/GC preservation must be checked empirically.
The implementation is architecture-compatible with the official repository but is
not a vendored invocation of its training script.

## Shared And Native Preprocessing

Each native generator receives both:

1. a shared preprocessing benchmark for controlled comparison; and
2. its paper-native preprocessing, including Vinas log/z-score and Lacan
   TPM/L1000/MaxAbs where gene lengths and ortholog mappings permit it.

GeneJEPA receives `log1p` followed by one scalar mean and standard deviation fitted
only over nonzero training expression values. Zero remains the absent-token sentinel,
matching the official sparse input contract. Its single-cell ragged-token input and
lack of a decoder are not disguised as a native bulk generator.

Lacan's native arm uses TPM followed by training-fold MaxAbs scaling. Mouse gene
lengths are generated from the union of GENCODE M39 exon intervals for each
versionless Ensembl gene ID. The table in
`data/reference/gencode_vM39_mouse_gene_lengths.tsv` covers all 980 mapped mouse
L1000 genes; its source URL and GTF SHA-256 are recorded in the adjacent manifest.

## Validation And Model Selection

No sample-random split is allowed when multiple accessions are pooled. Entire OSDR
accessions and ARCHS4 GEO series are grouped during splitting.

Checkpoint and hyperparameter selection use a held-out-validation fidelity composite:

- gene mean, variance, zero fraction, and quantile agreement;
- gene-gene correlation and pathway-correlation agreement;
- precision, recall, and density/coverage in a train-fitted embedding;
- adversarial real-versus-synthetic accuracy;
- nearest-neighbor memorization and duplicate checks;
- tissue and condition consistency.

The composite is subject to hard diversity and memorization gates. A model cannot win
by collapsing toward an average profile. FLT/GC effect recovery and classifier utility
are secondary validation metrics, not the sole optimization target.

The current executable composite averages gene-mean correlation, gene-standard-
deviation correlation, precision/recall F1, and adversarial indistinguishability.
Eligibility additionally requires manifold recall of at least 0.1, synthetic global
standard deviation of at least 10% of held-out real data, and no more than 5% of
synthetic samples closer to training profiles than the training leave-one-out first
percentile. These thresholds are benchmark gates, not biological significance tests.

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
plans live under `outputs/generative_benchmark/`. Model-specific training output will
be written under a run ID derived from the complete resolved configuration.

## Implementation Status

Implemented in `src/nasa_mouse_generative/`:

- three-model capability and paper/code provenance registry;
- validated configuration schema and 463-row gated experiment plan;
- current OSDR API inventory, technical-replicate-aware expression builder, and
  tissue eligibility tiers;
- full ARCHS4 metadata scan and balanced reference manifests;
- accession-grouped locked and LOO split plans;
- fold-aware shared preprocessing and harmonization constraints;
- reloadable ComBat, R `sva::ComBat_seq`, and MOBER harmonization adapters with
  fold behavior and fallback audits;
- executable Vinas WGAN-GP, Lacan diffusion, and GeneJEPA representation adapters;
- full-H5 ARCHS4 extraction with content-addressed caching and hierarchical sampling;
- direct OSDR, ARCHS4-only, and ARCHS4-pretrain/OSDR-fine-tune stage orchestration;
- epoch checkpoints, deterministic resume, final model serialization, and GPU records;
- conditioned generation, inverse preprocessing, held-out fidelity/diversity/
  memorization metrics, and real/synthetic FLT/GC classifier evaluation;
- locked-test enforcement and pooled, individual-tissue, or all-eligible-tissue CLI runs;
- resumable experiment-matrix execution, bounded rows per invocation, durable status,
  disk guards, sparse checkpoint retention, and ranked scoreboard generation;
- ARCHS4-only GEO-series train/validation/test partitions and held-out tissue plots;
- GeneJEPA held-out tissue UMAP and Lacan-style DDIM `t=1000,200,0` PCA trajectory.

GeneJEPA uses the pinned official model source at commit
`a2f4d7218b17f2f52cc5f1cc94420c8ef1ae3265`. It is evaluated only as a
representation model and never appears in synthetic-expression results.

The configurable Lacan adapter was checked against official code commit
`cde890154698fcea96c924804aaff04af3351b48`, but it remains an extension for
arbitrary OSDR covariates and pretrain/fine-tune stages. Exact architecture and
training-procedure comparisons use the separate pinned implementation under
`src/nasa_mouse_rna_diffusion/`; see `docs/rna_diffusion_paper_parity.md`.

The `paper_native` hyperparameter profiles enforce the source configurations and
reject locked-value overrides. Full
GeneJEPA (`d=768`, 24 blocks, 50 epochs) and Lacan (`8192` hidden units, batch 2048,
15,000 epochs) runs are expensive. `practical_screen` retains the model family,
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

The practical GeneJEPA bulk adaptation used 4,096 train-fold HVGs, the official
nonzero `log1p` global standardization, a 128-dimensional/4-block model, and 10
epochs. On 20 tissues represented in held-out validation series, a linear probe on
the learned embedding reached 0.453 balanced accuracy and 0.423 macro F1. The same
probe on preprocessed expression reached 0.839 balanced accuracy and 0.840 macro
F1. Embedding and UMAP silhouettes were -0.223 and -0.233. The model retained some
tissue information, but the UMAP did not reproduce the paper's clean cell-type
clusters. Results and coordinates are under
`outputs/generative_benchmark/runs/genejepa/archs4_genejepa_tissues_paper_preprocessing_v4/figures/archs4_tissues_validation/`.

The bounded exact-architecture follow-up retained the released 768-wide,
512-latent, 24-block, 12-head model and processed 43,744 replacement-sampled
profiles. It completed in 1,074 training seconds with 31.68 GB peak allocated A100
memory. On the same balanced 5,815-training/971-held-out figure cohort, tissue
balanced accuracy and macro F1 improved to 0.703 and 0.701, but remained below the
expression baseline of 0.839 and 0.840. Embedding and UMAP silhouettes remained
negative at -0.176 and -0.215. The run made 238 optimizer updates, below the
released 2,000-step EMA warmup and far below the paper's 50-million-profile-exposure
duration. It is therefore a successful architecture/runtime test but does not pass
the representation gate for GeneJEPA-guided diffusion. Its outputs are under
`outputs/generative_benchmark/runs/genejepa/matrix_phase_0_genejepa_exact_mouse_one_epoch_f2e01cf1f130d5cb/`.

The earlier reduced Lacan proxy and its failed 100/500-epoch output directories were
removed. They changed architecture, optimizer behavior, sampling weights, training
duration, input normalization, and split policy simultaneously, so they could not
attribute a result to mouse data. The replacement uses the unmodified upstream
227,109,786-parameter `ModelDDIM`, the paper's 9,796/2,448/5,000 split sizes, and
full-transcriptome TPM before selecting 974 mouse landmark genes. Its run directory
is `outputs/generative_benchmark/runs/lacan_diffusion/archs4_mouse_paper_parity_seed1234/`.
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

The upstream command checks out the pinned GeneJEPA source under ignored
`assets/model_sources/`. WGAN and diffusion reuse the repository's existing
`nasa_mouse_wgan` and `nasa_mouse_diffusion` PyTorch cores.

## Training Commands

Run the bounded end-to-end WGAN check first:

```bash
PYTHONPATH=src python -m nasa_mouse_generative train \
  --config configs/generative/default.yaml \
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

For GeneJEPA, set both `training.model=genejepa` and
`training.task=representation`. An ARCHS4-only generator must set
`training.condition_on_flight=false`: ARCHS4 has no flight labels, so it is only a
tissue/reference baseline and cannot identify a spaceflight effect.

Train GeneJEPA on ARCHS4 only with tissue retained as an evaluation label. The model
itself is not tissue-conditioned and cannot generate expression:

```bash
PYTHONPATH=src python -m nasa_mouse_generative train \
  --set training.model=genejepa \
  --set training.task=representation \
  --set training.regime=archs4_only \
  --set training.condition_on_flight=false \
  --set 'training.conditioning_covariates=[tissue]' \
  --set preprocessing.profile=model_native \
  --set features.space=hvg \
  --set features.hvg_genes=4096 \
  --set execution.device=cuda \
  --run-name archs4_genejepa_tissues
```

Prepare and train the paper-parity tissue-conditioned ARCHS4 diffusion baseline:

```bash
PYTHONPATH=src python -m nasa_mouse_rna_diffusion prepare
PYTHONPATH=src python -m nasa_mouse_rna_diffusion train
PYTHONPATH=src python -m nasa_mouse_rna_diffusion evaluate
```

Selected archival runs can save `prepared_data.h5`; loading remains compatible with
the older `prepared_osdr.h5` filename. The storage-safe default records profile and
accession hashes and deterministically rebuilds matrices for later evaluation. In an
ARCHS4-only run, `reference` is the ARCHS4 training partition and no OSDR expression
file is opened.

Create paper-style held-out tissue figures after the configuration is fixed:

```bash
PYTHONPATH=src python -m nasa_mouse_generative archs4-figures \
  --run-dir outputs/generative_benchmark/runs/genejepa/RUN_ID \
  --split validation --device cuda

PYTHONPATH=src python -m nasa_mouse_rna_diffusion evaluate
```

The GeneJEPA command fits UMAP only on held-out embeddings using the official test
callback defaults. The paper-parity diffusion evaluation fits one PCA basis only on
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
  --run-dir outputs/generative_benchmark/runs/vinas_wgan_gp/RUN_ID \
  --split test --unlock-test
```

Generate FLT samples from an observed joint covariate profile, overriding fields as
needed:

```bash
PYTHONPATH=src python -m nasa_mouse_generative generate \
  --run-dir outputs/generative_benchmark/runs/vinas_wgan_gp/RUN_ID \
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

ComBat, ComBat-seq, and MOBER now have executable, reloadable adapters. ComBat-seq
uses R 4.5 and Bioconductor `sva`; the other two run in the Python environment.
One-epoch outputs under `outputs/generative_benchmark/runs/` validate mechanics and
GPU/R integration only. They are not biological results or model-selection evidence.
