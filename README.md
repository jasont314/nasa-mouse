# NASA Mouse Spaceflight Models

This repository collects project-specific workflows for mouse spaceflight gene
expression modeling. It is no longer only a GLARE workspace: the current repo
also tracks OSDR/TMS preprocessing, expiMap/scArches setup, Reactome mouse
pathway architecture files, ARCHS4 mouse resources, and downstream analysis.

## Repository Layout

- `src/nasa_mouse_glare/`: project code for shared OSDR/TMS preprocessing,
  GLARE adaptation, Reactome GMT generation, and analysis utilities.
- `assets/model_sources/glare/`: vendored GLARE source with local runtime fixes.
- `src/nasa_mouse_expimap/`: NASA mouse expiMap preparation, training, query
  mapping, analysis, and publication workflows. The model uses the installed
  `scarches` package.
- `src/nasa_mouse_generative/`: configurable WGAN and diffusion benchmark with
  OSDR-only, ARCHS4-only, and pretrain/fine-tune regimes.
- `src/nasa_mouse_wgan/`: conditional WGAN-GP training, generation, and feature
  analysis.
- `src/nasa_mouse_diffusion/`: practical conditional diffusion backend and the
  exact Lacan DDIM workflow under `paper_parity/`.
- [`configs/generative/`](configs/generative/): shared benchmark settings plus
  separate DDIM and WGAN run configurations.
- `data/reference/expimap/`: expiMap paper metadata used by GLARE and expiMap
  validation.
- [`data/pathways/reactome_current_mouse_ensembl.gmt`](data/pathways/reactome_current_mouse_ensembl.gmt):
  generated Reactome mouse GMT file for the expiMap architecture mask.
- `assets/archs4/mouse_gene_v2.5.h5`: local ARCHS4 mouse H5 resource; ignored
  by git because it is large.
- `data/osdr_api/`: NASA OSDR Biological Data API metadata and small
  manifests; downloaded count CSVs under `data/osdr_api/counts/` are ignored.
- [`literature.md`](literature.md): links for GLARE, expiMap, MOBER, WGAN-GP,
  and DDIM.
- [`docs/osdr_api.md`](docs/osdr_api.md): NASA OSDR Biological Data API notes
  and examples.
- [`docs/method_sources.md`](docs/method_sources.md): locations and pinned
  revisions for upstream model and harmonization implementations.
- [`outputs/README.md`](outputs/README.md): organized output layout and the
  exact runs selected for analysis.
- [`outputs/COMMANDS.md`](outputs/COMMANDS.md): canonical commands used to
  reproduce the retained output families.

## Setup

Run workflow commands from the repository root:

```bash
cd path/to/nasa-mouse
conda activate nasa-mouse
export PYTHONPATH=src
```

To create or refresh the local environment:

```bash
conda env create -f environment.yml
# For an existing environment:
conda env update -f environment.yml --prune
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

Outputs are written under `outputs/glare/multi_tissue_api/`. Retina is audited
but skipped for GLARE unless a matching TMS FACS retina pretraining source is
added. Skeletal-muscle subtype runs use official OSDR material-type labels and
the available TMS FACS `limb muscle` pretraining tissue.

Run one prepared aggregate scope:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare run-glare-scope \
  --scope-dir outputs/glare/multi_tissue_api/liver/aggregate
```

Run MOBER-corrected aggregate GLARE for a multi-study scope:

```bash
PYTHONPATH=src:assets/model_sources/MOBER python -m nasa_mouse_glare.multi_tissue_api_glare run-mober-scope \
  --scope-dir outputs/glare/multi_tissue_api/liver/aggregate
```

Run all per-study GLARE scopes for one tissue and compare against per-study
DESeq2:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare run-per-study-glare \
  --tissue-dir outputs/glare/multi_tissue_api/liver

PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_api_glare run-dgea-comparison \
  --tissue-dir outputs/glare/multi_tissue_api/liver
```

Run the paper-style validation stack on the generated multi-tissue GLARE
outputs:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.multi_tissue_validation \
  --include-per-study \
  --include-mober \
  --shap-aggregate
```

This writes XGBoost verification, representation QC, clustering QC,
DEG-enrichment comparisons, intersection-vs-GLARE-only module-score validation,
Panglao marker enrichment, and Metascape-ready gene lists to
`outputs/glare/multi_tissue_api/validation_stack/`.

Prepare tissue-specific expiMap inputs from API count tables:

```bash
PYTHONPATH=src python -m nasa_mouse_expimap.prepare_expimap_osdr_tissue --tissue liver
PYTHONPATH=src python -m nasa_mouse_expimap.prepare_expimap_osdr_tissue --tissue kidney
```

## Workflows

The GLARE-compatible preprocessing and fine-tuning workflow is documented in
[`src/nasa_mouse_glare/README.md`](src/nasa_mouse_glare/README.md).

The expiMap/scArches handoff and architecture notes are documented in
[`docs/expimap_handoff.md`](docs/expimap_handoff.md).

Current expiMap run summaries and preprocessing comparisons are documented in
[`docs/expimap_results.md`](docs/expimap_results.md).

The liver query-extension de novo-program analysis is documented in
[`docs/expimap_de_novo_liver.md`](docs/expimap_de_novo_liver.md).

The tutorial-style liver expiMap run with HVG filtering, a deeper reference
model, and HSIC de novo query nodes is documented in
[`docs/expimap_tutorial_style_liver.md`](docs/expimap_tutorial_style_liver.md).

Accession-aware direct-model validation and the larger ARCHS4 reference
seed-stability result are documented in
[`docs/expimap_accession_validation.md`](docs/expimap_accession_validation.md)
and [`docs/expimap_reference_seed_stability.md`](docs/expimap_reference_seed_stability.md).
Condition-specific GC/FLT expiMap clustering is documented in
[`docs/expimap_condition_clustering.md`](docs/expimap_condition_clustering.md).

The skeletal-muscle pathway prior-work check is documented in
[`docs/expimap_skeletal_muscle_prior_work.md`](docs/expimap_skeletal_muscle_prior_work.md).

The configurable bulk generative-model benchmark, ARCHS4 cohort audit,
conditioning, harmonization, training regimes, and tissue-figure commands are
documented in
[`docs/generative_models_pipeline.md`](docs/generative_models_pipeline.md).
Model-specific WGAN-GP and diffusion workflows are documented in
[`docs/wgan_pipeline.md`](docs/wgan_pipeline.md) and
[`docs/diffusion_pipeline.md`](docs/diffusion_pipeline.md).
The exact upstream-architecture mouse diffusion comparison is documented in
[`docs/rna_diffusion_paper_parity.md`](docs/rna_diffusion_paper_parity.md).

For NASA OSDR programmatic data access, see
[`docs/osdr_api.md`](docs/osdr_api.md).
