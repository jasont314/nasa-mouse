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
- fold-aware shared preprocessing and harmonization constraints.

The existing `src/nasa_mouse_wgan/` and `src/nasa_mouse_diffusion/` model cores will
be adapted behind this interface. Full benchmark training is deliberately after the
data/split checkpoint so the old fixed 24,428-sample inputs and sample-random
validation are not silently reused.
