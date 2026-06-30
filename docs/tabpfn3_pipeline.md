# TabPFN3 OSDR Classification Pipeline

This project evaluates TabPFN3 as an OSDR-only classifier for NASA mouse bulk
RNA-seq flight-vs-ground-control prediction. It intentionally does not use
ARCHS4 and does not use the older raw integrated OSDR H5 files.

## Paper/API Review

TabPFN3 is a tabular foundation model for supervised classification and
regression. The paper describes a transformer trained by prior-data fitted
networks on 2.75 million synthetic datasets, then used at inference time for
new tabular tasks through in-context learning rather than task-specific gradient
training. The official package exposes `TabPFNClassifier` and
`TabPFNRegressor`, and the package README states that the latest installation
uses TabPFN-3 by default.

Practical constraints for this repository:

- TabPFN3 is discriminative, not generative. It is appropriate here for
  predicting `flight` vs `ground_control` and ranking features, not for
  synthesizing RNA-seq samples.
- The official Python input format is standard tabular supervised learning:
  numeric `X` with shape `n_samples x n_features` and a target vector `y`.
  The classifier exposes scikit-learn-like `fit`, `predict`, and
  `predict_proba` methods.
- TabPFN3 does not pretrain on ARCHS4 in this pipeline. The user request for
  this task was OSDR-only.
- The package can run on CUDA. The local A100 was detected and the backend
  check selected `cuda`.
- Local inference requires Prior Labs license acceptance and model-weight
  access. The production run used `TABPFN_TOKEN` from the local ignored `.env`
  file and ran on CUDA.
- The installed package defaults to model version `v3`. The pipeline now
  requests `ModelVersion.V3` explicitly through
  `TabPFNClassifier.create_default_for_version`.
- Installed package inspection shows V3 checkpoints are loaded from the gated
  `Prior-Labs/tabpfn_3` repository, with default checkpoint
  `tabpfn-v3-classifier-v3_default.ckpt`.
- Exact V3 pretraining limits are stored in the loaded model inference config,
  so they cannot be confirmed locally until authorized weights are available.
  The installed classifier documentation does state that earlier v2.5 defaults
  support 50,000 rows and 2,000 columns, with 500-feature subsampling per
  estimator above 500 columns. For RNA-seq, this means feature selection and a
  sufficiently large estimator ensemble are required for interpretable coverage.
- Probability output is available via `predict_proba`; the package also exposes
  softmax temperature, probability balancing, and tuning options. This pipeline
  records Brier score as a calibration-sensitive metric but does not fit a
  separate calibration model.
- No TabPFN3-specific gene attribution method was identified as validated for
  high-dimensional RNA-seq. The implemented feature importance is held-out
  permutation importance, restricted to top fold-local candidates for runtime.

Sources:

- TabPFN3 paper: <https://arxiv.org/pdf/2605.13986>
- Official package repository: <https://github.com/PriorLabs/TabPFN>
- Prior Labs documentation: <https://docs.priorlabs.ai/>

## Implementation

New source package:

```text
src/nasa_mouse_tabpfn3/
```

Main entry point:

```bash
PYTHONPATH=src python -m nasa_mouse_tabpfn3.run_osdr_classification
```

The runner reads only API-derived OSDR files:

- metadata: `data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv`
- count CSV cache: `data/osdr_api/counts/`

The loader keeps `Mus musculus` bulk RNA-seq samples with
`flight`/`ground_control` labels from the existing NASA OSDR API workflow. By
default, feature matrices are restricted to mouse Ensembl gene IDs with prefix
`ENSMUSG` so ERCC spike-ins do not enter feature importance.

## Planned Cohorts

Primary tissues:

- liver
- skeletal_muscle
- skin
- kidney
- thymus
- spleen
- lung
- retina

Skeletal-muscle split groups:

- soleus
- gastrocnemius
- quadriceps
- edl
- tibialis_anterior

## Modeling Design

For each tissue or split group:

1. Build raw-count and `log1p(CPM)` matrices from OSDR API count tables.
2. Fit feature selection inside each training fold.
3. Evaluate both feature modes:
   - `all_expressed`: expressed/variable `ENSMUSG` genes after fold-local filtering
   - `hvg`: fold-local top highly variable genes after the same filter
4. Run CV schemes:
   - random stratified CV
   - accession-aware grouped CV
   - leave-one-accession-out CV where feasible
5. Write predictions, metrics, permutation feature importance, and plots.

The production backend is `tabpfn`. A `sklearn_logreg` backend exists only for
smoke validation of data loading, folds, plotting, and output writing.

Leakage controls:

- CPM/log transform is unsupervised and computed from each sample's own library
  size.
- Expression prevalence filtering, HVG/variance selection, and univariate
  candidate ranking are fit only on each training fold.
- Held-out fold samples are never used to select genes or tune feature ranking.
- Accession-aware grouped CV and leave-one-accession-out are written separately
  from random stratified CV so study-structure leakage can be seen directly.
- Optional design covariates are one-hot encoded inside each training fold with `--include-design-covariates`. Target-condition fields and sample/file/profile identifiers are not used as covariates.

## Production Command

The completed OSDR-only production run used:

```bash
PYTHONPATH=src python -m nasa_mouse_tabpfn3.run_osdr_classification \
  --output-root outputs/tabpfn3_osdr \
  --backend tabpfn \
  --device cuda \
  --feature-modes all_expressed hvg \
  --cv-schemes random grouped loo_accession \
  --hvg-top-n 500 \
  --max-features 500 \
  --importance-candidates 5 \
  --permutation-repeats 1 \
  --n-estimators 3
```

The runner also reads `TABPFN_TOKEN` from a local `.env` file when the variable
is not already exported. The `.env` file is ignored by git and should not be
committed.

This run intentionally capped each fold at 500 selected genes for practical
runtime. With `--hvg-top-n 500` and `--max-features 500`, the `all_expressed`
and `hvg` tracks collapse to the same fold-local top-variance feature set.

Expected summary outputs:

- `outputs/tabpfn3_osdr/summary/tabpfn3_metrics.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_aggregate_metrics.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_predictions.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_feature_importance.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_run_manifest.tsv`
- `outputs/tabpfn3_osdr/summary/osdr_tissue_inventory.tsv`
- `outputs/tabpfn3_osdr/summary/osdr_tissue_accession_inventory.tsv`
- `outputs/tabpfn3_osdr/summary/osdr_tissue_design_inventory.tsv`
- `outputs/tabpfn3_osdr/summary/osdr_muscle_split_inventory.tsv`
- `outputs/tabpfn3_osdr/summary/osdr_muscle_design_inventory.tsv`
- `outputs/tabpfn3_osdr/summary/osdr_sample_inventory.tsv`

The design inventories include condition counts by tissue, accession, material
type, sex, platform, assay, data source accession, and project type where those
fields are present in the NASA API metadata. The sample inventory preserves one
row per selected profile with the same covariates.

## Validation Runs

The official TabPFN3 backend was smoke-tested locally on CUDA after loading the
token from `.env`; the smoke fit completed and returned probability outputs.

The code also supports a non-biological `sklearn_logreg` backend for mechanics
testing:

```bash
PYTHONPATH=src python -m nasa_mouse_tabpfn3.run_osdr_classification \
  --output-root outputs/tabpfn3_osdr_smoke_ensembl \
  --backend sklearn_logreg \
  --tissues liver \
  --no-muscle-splits \
  --feature-modes hvg \
  --cv-schemes random \
  --hvg-top-n 50 \
  --max-features 50 \
  --importance-candidates 5 \
  --permutation-repeats 1
```

That smoke run is not a TabPFN3 biological result. It validates that the OSDR
API data path, fold-local feature selection, metrics, feature importance, and
plots work.

Production results, including the covariate-augmented comparison, are documented in `docs/tabpfn3_results.md`.
