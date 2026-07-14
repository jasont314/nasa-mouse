# NASA Mouse Generative Model Benchmark

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
- balancing, seed, repeats, generated sample count, and synthetic-to-real ratio;
- technical-replicate handling (`keep`, `sum`, or `mean`, with `sum` as default);
- optional conditioning on tissue, material type, muscle group, study, sex, assay,
  platform, and data source.

Normalization, transformation, scaling, and harmonization are separate parameters.
All fitted statistics use training folds only. Study-wise scaling for an unseen study
uses training-global fallback by default. Estimating statistics from the held-out
study is an explicitly labeled transductive sensitivity analysis.

## Harmonization Arms

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

## Shared And Native Preprocessing

Each native generator receives both:

1. a shared preprocessing benchmark for controlled comparison; and
2. its paper-native preprocessing, including Vinas log/z-score and Lacan
   TPM/L1000/MaxAbs where gene lengths and ortholog mappings permit it.

GeneJEPA receives its global scalar `log1p` standardization in its own adapter. Its
single-cell ragged-token input and lack of a decoder are not disguised as a native
bulk generator.

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
- validated configuration schema and 458-row gated experiment plan;
- current OSDR API inventory, technical-replicate-aware expression builder, and
  tissue eligibility tiers;
- full ARCHS4 metadata scan and balanced reference manifests;
- accession-grouped locked and LOO split plans;
- fold-aware shared preprocessing and harmonization constraints;
- executable Vinas WGAN-GP, Lacan diffusion, and GeneJEPA representation adapters;
- full-H5 ARCHS4 extraction with content-addressed caching and hierarchical sampling;
- direct OSDR, ARCHS4-only, and ARCHS4-pretrain/OSDR-fine-tune stage orchestration;
- epoch checkpoints, deterministic resume, final model serialization, and GPU records;
- conditioned generation, inverse preprocessing, held-out fidelity/diversity/
  memorization metrics, and real/synthetic FLT/GC classifier evaluation;
- locked-test enforcement and pooled, individual-tissue, or all-eligible-tissue CLI runs.

GeneJEPA uses the pinned official model source at commit
`a2f4d7218b17f2f52cc5f1cc94420c8ef1ae3265`. It is evaluated only as a
representation model and never appears in synthetic-expression results.

## Installation

```bash
conda activate nasa-mouse
python -m pip install -r requirements-nasa-mouse-generative.txt
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

One-study runs use `data.osdr_accession_scope=single` with exactly one
`data.osdr_include_accessions` value. Because accession holdout is impossible in
that design, the run manifest labels its deterministic condition-stratified sample
split. Multi-study runs retain accession-grouped validation whenever possible.

## Evaluation And Generation

Training automatically evaluates held-out validation accessions. Representation
metrics neutralize the supplied condition token before FLT/GC classification so the
model cannot score by reading its conditioning label. The final test remains locked:

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
vocabularies, prepared OSDR partitions, epoch history, latest resumable checkpoint,
final model, validation metrics, PCA/fidelity plots, embeddings, device record, and
concise README.

ComBat, ComBat-seq, and MOBER remain registered experimental harmonization arms but
do not yet have executable adapters in this runner. Selecting one fails explicitly;
`none`, `within_study_zscore`, and `within_study_then_global_zscore` are executable.
The completed one-epoch outputs under `outputs/generative_benchmark/runs/` are
mechanics tests, not biological results.
