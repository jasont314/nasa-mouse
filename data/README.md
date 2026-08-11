# Generated Data and Local Inputs

This directory contains generated preprocessing outputs, pathway architecture
files, and small manifests. Large matrix artifacts should stay out of git; the
generated matrices are ignored by `.gitignore`.

Original GLARE input notes were generated on 2026-06-12 in the `nasa` conda
environment.

Reactome mouse pathway files were added on 2026-06-26 for expiMap/scArches.

## Inputs

- TMS FACS `.h5ad`: `assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad`
- OSDR mouse bulk RNA-seq FLT/GC metadata:
  `data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv`
- ARCHS4 mouse H5: `assets/archs4/mouse_gene_v2.5.h5`
- expiMap Reactome mouse GMT:
  `data/pathways/reactome_current_mouse_ensembl.gmt`
- paper-parity mouse diffusion landmarks:
  `data/diffusion/l974_mouse_paper_parity.tsv`

The expiMap OSDR inputs should be built from the NASA OSDR Biological Data API.
Do not use the older local integrated OSDR HDF5 as the expiMap OSDR source.
Downloaded API count CSVs are cached under `data/osdr_api/counts/` and ignored
by git.

The 974-gene diffusion panel is generated deterministically from the paper's GTEx
landmarks, the Ensembl human-to-mouse map, and ARCHS4 gene annotations. Its adjacent
manifest records the ten human landmarks without a direct mouse dimension and the
ten additional paralogs used to preserve the published input width. Large prepared
TPM matrices remain under `outputs/generative/benchmark/data/lacan_diffusion/`.

## Pathways

`data/pathways/reactome_current_mouse_ensembl.gmt` is generated from official
current Reactome files and is the expiMap architecture source. Each row is one
mouse Reactome pathway, the second column is the Reactome browser URL, and the
remaining columns are mouse Ensembl gene IDs (`ENSMUSG...`).

Regenerate it from the repository root:

```bash
PYTHONPATH=src python src/nasa_mouse_glare/build_reactome_mouse_gmt.py
```

## Commands

Run these from the repository root after activating the current `nasa-mouse`
environment:

```bash
cd path/to/nasa-mouse
```

Discover NASA OSDR API mouse bulk RNA-seq Space Flight/Ground Control samples:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics
```

Prepare direct expiMap OSDR tissue inputs:

```bash
PYTHONPATH=src python -m nasa_mouse_expimap.prepare_expimap_osdr_tissue --tissue liver
PYTHONPATH=src python -m nasa_mouse_expimap.prepare_expimap_osdr_tissue --tissue kidney
```

Prepare GLARE expression directly from NASA OSDR API metadata and per-accession
unnormalized count tables:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics \
  --download-counts

PYTHONPATH=src python -m nasa_mouse_glare.osdr \
  --tissue liver \
  --output-dir outputs/glare/api_inputs/liver \
  --download-counts
```

The GLARE API loader writes raw-count and `log2(CPM+1)` matrix bundles plus
study-aware DESeq2 inputs. The multi-tissue training wrapper also prepares and
aligns the matching TMS reference automatically.

## Output Shapes

- `data/processed/tms_facs_3552_cells`: `21025 genes x 3552 cells`
- `data/processed/osdr_mouse_bulk`: `53511 genes x 3315 samples`
- `data/processed/tms_facs_osdr_aligned.pretrain`: `21010 genes x 3552 cells`
- `data/processed/tms_facs_osdr_aligned.target`: `21010 genes x 3315 samples`
- `data/glare_inputs/tms_facs_pretrain.mtx`: `21010 genes x 3552 cells`
- `data/glare_inputs/osdr_finetune.csv`: `21010 genes x 3315 samples`

The GLARE source is vendored directly in `assets/model_sources/glare` with the `hpt.py` runtime
fix applied: direct script execution works, MatrixMarket input uses SciPy's
`.toarray()` API, fine-tuning reuses the pretraining architecture, and final
representation extraction applies the fine-tuning adapter.
