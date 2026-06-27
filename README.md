# NASA Mouse Spaceflight Models

This repository collects project-specific workflows for mouse spaceflight gene
expression modeling. It is no longer only a GLARE workspace: the current repo
also tracks OSDR/TMS preprocessing, expiMap/scArches setup, Reactome mouse
pathway architecture files, ARCHS4 mouse resources, and downstream analysis.

## Repository Layout

- `src/nasa_mouse_glare/`: project code for OSDR/TMS preprocessing, GLARE
  adaptation, Reactome GMT generation, and analysis utilities.
- `src/glare/`: vendored GLARE source with local runtime fixes.
- `src/expiMap_scarches/`: expiMap/scArches source and handoff notes.
- [`data/pathways/reactome_current_mouse_ensembl.gmt`](data/pathways/reactome_current_mouse_ensembl.gmt):
  generated Reactome mouse GMT file for the expiMap architecture mask.
- `assets/archs4/mouse_gene_v2.5.h5`: local ARCHS4 mouse H5 resource; ignored
  by git because it is large.
- `data/osdr_api/`: NASA OSDR Biological Data API metadata and small
  manifests; downloaded count CSVs under `data/osdr_api/counts/` are ignored.
- [`literature.md`](literature.md): links for GLARE, VEGA, expiMap, OntoVAE,
  and MOBER.
- [`docs/osdr_api.md`](docs/osdr_api.md): NASA OSDR Biological Data API notes
  and examples.

## Setup

Run workflow commands from the repository root:

```bash
cd path/to/nasa-mouse
conda activate nasa-mouse
export PYTHONPATH=src
```

To create or refresh the local environment:

```bash
conda create -y -n nasa-mouse python=3.11
conda run -n nasa-mouse python -m pip install -r requirements-nasa-mouse-glare.txt
```

## Current Inputs

- OSDR mouse bulk RNA-seq FLT/GC metadata and count tables are discovered from
  the NASA OSDR Biological Data API, not from the older local integrated OSDR
  HDF5.
- ARCHS4 mouse gene expression H5:
  `assets/archs4/mouse_gene_v2.5.h5`
- Reactome mouse expiMap architecture GMT:
  `data/pathways/reactome_current_mouse_ensembl.gmt`

The Reactome GMT is generated from official current Reactome files:

- `ReactomePathways.txt`
- `Ensembl2Reactome_All_Levels.txt`

Regenerate it with:

```bash
PYTHONPATH=src python src/nasa_mouse_glare/build_reactome_mouse_gmt.py
```

The output GMT uses one row per mouse Reactome pathway:

```text
R-MMU-73857_RNA_POLYMERASE_II_TRANSCRIPTION    https://reactome.org/PathwayBrowser/#/R-MMU-73857    ENSMUSG...
```

Discover OSDR mouse bulk RNA-seq Space Flight/Ground Control samples:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics
```

Audit and prepare API-native multi-tissue GLARE inputs:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare audit

PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare prepare \
  --tissue all \
  --download-counts \
  --prepare-per-study
```

Outputs are written under `outputs/glare_multi_tissue_api/`. Retina is audited
but skipped for GLARE unless a matching TMS FACS retina pretraining source is
added. Skeletal-muscle subtype runs use official OSDR material-type labels and
the available TMS FACS `limb muscle` pretraining tissue.

Run one prepared aggregate scope:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare run-glare-scope \
  --scope-dir outputs/glare_multi_tissue_api/liver/aggregate
```

Run MOBER-corrected aggregate GLARE for a multi-study scope:

```bash
PYTHONPATH=src:src/MOBER python -m nasa_mouse_glare.multi_tissue_api_glare run-mober-scope \
  --scope-dir outputs/glare_multi_tissue_api/liver/aggregate
```

Run all per-study GLARE scopes for one tissue and compare against per-study
DESeq2:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare run-per-study-glare \
  --tissue-dir outputs/glare_multi_tissue_api/liver

PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare run-dgea-comparison \
  --tissue-dir outputs/glare_multi_tissue_api/liver
```

Prepare tissue-specific expiMap inputs from API count tables:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.prepare_expimap_osdr_tissue --tissue liver
PYTHONPATH=src python -m nasa_mouse_glare.prepare_expimap_osdr_tissue --tissue kidney
```

## Workflows

The GLARE-compatible preprocessing and fine-tuning workflow is documented in
[`src/nasa_mouse_glare/README.md`](src/nasa_mouse_glare/README.md).

The expiMap/scArches handoff and architecture notes are documented in
[`src/expiMap_scarches/EXPIMAP_HANDOFF.md`](src/expiMap_scarches/EXPIMAP_HANDOFF.md).

Current expiMap run summaries and preprocessing comparisons are documented in
[`docs/expimap_results.md`](docs/expimap_results.md).

For NASA OSDR programmatic data access, see
[`docs/osdr_api.md`](docs/osdr_api.md).
