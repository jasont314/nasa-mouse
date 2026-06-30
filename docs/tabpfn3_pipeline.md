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
- TabPFN3 does not pretrain on ARCHS4 in this pipeline. The user request for
  this task was OSDR-only.
- The package can run on CUDA. The local A100 was detected and the backend
  check selected `cuda`.
- Local inference requires Prior Labs license acceptance and model-weight
  access. In this non-interactive environment the official package stops unless
  `TABPFN_TOKEN` is set or weights are already cached.

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

## Production Command

After accepting the Prior Labs license and exporting a token:

```bash
export TABPFN_TOKEN="..."
PYTHONPATH=src python -m nasa_mouse_tabpfn3.run_osdr_classification \
  --output-root outputs/tabpfn3_osdr \
  --backend tabpfn \
  --device cuda \
  --feature-modes all_expressed hvg \
  --cv-schemes random grouped loo_accession \
  --importance-candidates 100 \
  --permutation-repeats 3
```

Expected summary outputs:

- `outputs/tabpfn3_osdr/summary/tabpfn3_metrics.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_predictions.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_feature_importance.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_run_manifest.tsv`
- `outputs/tabpfn3_osdr/summary/osdr_tissue_inventory.tsv`
- `outputs/tabpfn3_osdr/summary/osdr_tissue_accession_inventory.tsv`
- `outputs/tabpfn3_osdr/summary/osdr_muscle_split_inventory.tsv`

## Validation Run

Because `TABPFN_TOKEN` is not available locally, the actual TabPFN3 backend is
blocked before fitting. The code was smoke-tested with:

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

