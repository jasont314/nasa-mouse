# NASA Mouse GLARE Workflow

This package contains the project-specific GLARE adaptation and the shared
NASA OSDR API utilities used by the other modeling workflows. The active GLARE
pipeline reads per-accession mouse bulk RNA-seq count tables from the NASA OSDR
Biological Data API and uses Tabula Muris Senis (TMS) FACS profiles as its
single-cell reference.

The upstream GLARE implementation is vendored under
`assets/model_sources/glare/`. Project-specific data preparation, training,
and validation remain in this package.

## Environment

Run commands from the repository root:

```bash
conda env create -f environment.yml
conda activate nasa-mouse
export PYTHONPATH=src
```

Use `conda env update -f environment.yml --prune` to refresh an existing
environment.

## Inputs

The active workflow uses:

- `data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv`: API-derived
  mouse bulk RNA-seq sample inventory.
- `data/osdr_api/counts/`: downloaded per-accession unnormalized count tables.
- `assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad`: TMS FACS reference.
- `data/pathways/reactome_current_mouse_ensembl.gmt`: current mouse Reactome
  pathways.

Large inputs and generated matrices are ignored by Git. See
[`data/README.md`](../../data/README.md) for their provenance and regeneration
commands.

## Discover OSDR Samples

Refresh the eligible Mus musculus bulk RNA-seq inventory and download count
tables:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics \
  --download-counts
```

The inventory includes all API-returned data sources with identifiable Flight
and Ground Control labels. Tissue and material labels are canonicalized by the
shared API loader.

## Audit and Prepare GLARE Inputs

Inspect tissue coverage against the available TMS reference:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare audit
```

Prepare aggregate and per-study scopes for all eligible tissues:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare prepare \
  --tissue all \
  --download-counts \
  --prepare-per-study
```

Preparation writes raw-count and `log2(CPM+1)` OSDR bundles, DESeq2 inputs,
the matching TMS subset, and aligned GLARE matrices under
`outputs/glare/multi_tissue_api/<tissue>/`. Retina is audited but cannot use
this GLARE workflow unless a matching single-cell reference is supplied.
Skeletal-muscle groups use official OSDR material labels and the TMS limb
muscle reference.

To prepare one tissue or selected accessions:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare prepare \
  --tissue liver \
  --accessions OSD-168 OSD-379 \
  --download-counts \
  --prepare-per-study
```

## Train GLARE

Train one prepared aggregate scope:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare run-glare-scope \
  --scope-dir outputs/glare/multi_tissue_api/liver/aggregate
```

Run every prepared accession separately:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare run-per-study-glare \
  --tissue-dir outputs/glare/multi_tissue_api/liver
```

Compare per-study GLARE modules with DESeq2 and Reactome enrichment:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare run-dgea-comparison \
  --tissue-dir outputs/glare/multi_tissue_api/liver
```

## MOBER Comparison

For a prepared multi-study scope, train MOBER and run GLARE on the corrected
expression matrix:

```bash
PYTHONPATH=src:assets/model_sources/MOBER \
  python -m nasa_mouse_glare.multi_tissue_api_glare run-mober-scope \
  --scope-dir outputs/glare/multi_tissue_api/liver/aggregate
```

MOBER is trained for the selected cohort; it is not a fixed preprocessing
transform that can be reused unchanged across tissues.

## Validation

Run the multi-tissue validation stack after training:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_validation \
  --include-per-study \
  --include-mober \
  --shap-aggregate
```

This produces representation and clustering checks, FLT-versus-GC module
comparisons, DESeq2 overlap, pathway and marker enrichment, and Metascape-ready
gene lists. Outputs remain under `outputs/glare/multi_tissue_api/`.

## Reactome Architecture

Regenerate the current mouse Reactome GMT from the official Reactome mapping
files:

```bash
PYTHONPATH=src python src/nasa_mouse_glare/build_reactome_mouse_gmt.py
```

The final GMT is tracked. Its downloaded source files under
`data/pathways/reactome_current/raw/` are a regenerable local cache.

## Historical OSD-379 Reproduction

`paper_finetune.py`, `paper_analysis.py`, `osd379_tissue_qc.py`, and
`deseq_glare_comparison.py` preserve the first OSD-379 paper-reproduction
analysis. They are not part of the active multi-tissue pipeline. Historical
NASA normalized-expression, differential-expression, or ISA files must be
passed explicitly where requested; these scripts no longer assume a deleted
`assets/osdr/` directory.

## Output and Command Ledgers

- [`outputs/README.md`](../../outputs/README.md) identifies the retained final
  runs and analyses.
- [`outputs/COMMANDS.md`](../../outputs/COMMANDS.md) records reproduction
  commands for GLARE, expiMap, DDIM, WGAN, and downstream analyses.
- [`docs/osdr_api.md`](../../docs/osdr_api.md) documents API endpoints and
  filtering rules.
