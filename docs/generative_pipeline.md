# Generative pipeline

## Purpose

The generative framework compares conditional DDIM and WGAN-GP models for mouse
bulk RNA-seq. It supports broad mouse tissue generation from ARCHS4 and
flight-versus-ground generation after adaptation to NASA OSDR. Data ingestion,
preprocessing, split construction, model training, and evaluation are separate
steps so that changes to one axis can be compared without changing the others.

The shared framework is implemented in `src/nasa_mouse_generative/`. DDIM and
WGAN-specific training code lives in `src/nasa_mouse_diffusion/` and
`src/nasa_mouse_wgan/`, respectively.

## Configuration axes

| Axis | Supported choices | Final DDIM branch |
|---|---|---|
| Expression | counts, CPM, TPM, log transforms, z-score, robust scaling, MaxAbs | full-transcriptome TPM followed by train-fitted MaxAbs scaling |
| Features | shared genes, fold-selected HVGs, Reactome genes, mapped mouse L1000 landmarks | 974 mapped mouse landmarks |
| Harmonization | none, within-study z-score, two-stage z-score, ComBat, ComBat-seq, three MBatch methods, MOBER | no global correction |
| Training source | OSDR only, ARCHS4 only, ARCHS4 pretraining plus OSDR adaptation | ARCHS4 pretraining plus OSDR adaptation |
| Cohort | selected studies, all studies, pooled tissues, or one tissue | all eligible studies in one tissue-conditioned model |
| Conditioning | FLT/GC, tissue, study, material, muscle group, sex, assay, platform, source | tissue, FLT/GC, accession, and material |
| Generator | conditional WGAN-GP or ModelDDIM | ModelDDIM |

Preprocessing parameters are fitted on training profiles only. ARCHS4 partitions
are grouped by GEO series. OSDR metadata retain accession identity so study-aware
splits and within-study FLT/GC evaluation remain possible.

## Final data branch

- The full ARCHS4 mouse file contains 997,515 profiles. Metadata filters select
  bulk-like, healthy-preferred, tissue-matched profiles.
- The selected OSDR-disjoint reference contains 17,244 profiles from 20 tissue
  classes. Complete GEO series define the 10,150/2,466/4,628
  train/validation/test split.
- The API-derived OSDR inventory contains 1,610 biological profiles from 75
  accessions after technical replicates are collapsed. The final development
  split contains 781 training, 536 validation, and 293 test profiles.
- GEO series linked to eligible OSDR accessions are excluded from the selected
  ARCHS4 reference before pretraining.

The current data inventory is frozen in
[`table_1_data_inventory.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_1_data_inventory.tsv).

## Selected models

The selected ARCHS4 model follows the Lacan et al. residual-MLP ModelDDIM: two
8,192-unit hidden layers, 1,000 diffusion steps, quadratic beta schedule, EMA,
mixed precision, and 15,000 epochs. OSDR adaptation adds factorized domain and
condition modules with study and material conditioning.

The WGAN comparator follows the Viñas et al. conditional WGAN-GP: 64-dimensional
noise, two 256-unit ReLU layers, categorical embeddings, five critic updates per
generator update, and gradient penalty 10.

The configurations used for the final runs are:

- `configs/generative/diffusion/archs4_mouse_paper_parity_osdr_disjoint.yaml`
- `configs/generative/diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint.yaml`
- `configs/generative/wgan/wgan_matched_study_conditioned.yaml`

## Execution flow

1. Query and canonicalize OSDR metadata through the NASA API.
2. Build the OSDR expression matrix and ARCHS4 reference catalog.
3. Fit preprocessing on training data and create grouped splits.
4. Train the ARCHS4 reference generator.
5. Adapt the conditional model to OSDR where required.
6. Generate repeated synthetic cohorts and evaluate each metric independently.
7. Run tissue-specific real-only and real-plus-synthetic classifiers and feature
   analyses using held-out real profiles.

The canonical commands are in [`outputs/COMMANDS.md`](../outputs/COMMANDS.md).
Selected run directories and locally retained checkpoints are listed in
[`outputs/README.md`](../outputs/README.md). Generated profiles are model outputs,
not additional biological replicates; biological association tests use real OSDR
profiles.
