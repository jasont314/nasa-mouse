# expiMap / scArches Handoff

Date: 2026-06-26

This file summarizes the current NASA mouse spaceflight transcriptomics context for handing the work to another agent. The current pivot is from GLARE to expiMap/scArches.

## Repository State

- Repo root: `/media/volume/mouse/nasa/nasa-mouse`
- Main branch: `main`
- Latest relevant commit: `8df9641 Add scarches expiMap dependencies`
- Conda env: `/home/exouser/miniforge3/envs/nasa-mouse`
- Project Python requirements: `requirements-nasa-mouse-glare.txt`

The dependency pins have been added to `requirements-nasa-mouse-glare.txt`.

## Important Package Caveat

There are two scArches code sources in this repo/environment:

1. Installed package in the `nasa-mouse` conda env:
   - `scarches==0.6.1`
   - Verified to expose `sca.models.EXPIMAP`

2. Local checkout:
   - `src/expiMap_scarches`
   - This is an older scArches checkout.
   - It does **not** contain modern `EXPIMAP`.
   - It mostly has older trVAE/scGEN/cell-decoder code.

Do not run expiMap scripts from inside `src/expiMap_scarches`, and do not add `src/expiMap_scarches` to `PYTHONPATH`, or it may shadow the installed package. Run from the repo root and import the installed `scarches`.

Quick verification command:

```bash
/home/exouser/miniforge3/envs/nasa-mouse/bin/python - <<'PY'
import scarches as sca
import anndata
import zarr
print("scarches", getattr(sca, "__version__", "no __version__"), sca.__file__)
print("anndata", anndata.__version__)
print("zarr", zarr.__version__)
print("has EXPIMAP", hasattr(sca.models, "EXPIMAP"))
PY
```

Known compatible pins now in `requirements-nasa-mouse-glare.txt`:

```text
anndata>=0.10,<0.12
zarr<3
scanpy>=1.10,<1.12
scarches==0.6.1
```

The `anndata` and `zarr` pins are required because `scarches 0.6.1` imports `anndata.read`, which is not present in `anndata 0.12.x`, and `anndata 0.11.x` requires `zarr<3`.

## Current User Goal

The user wants to test expiMap as a pathway-aware model for NASA mouse OSDR spaceflight transcriptomics, using API-derived OSDR bulk RNA-seq data rather than the older local integrated OSDR HDF5. Liver and kidney should both be prepared and analyzed at minimum.

Two planned approaches:

1. Direct OSDR expiMap
   - Train expiMap directly on API-derived OSDR tissue bulk RNA-seq samples.
   - Extract pathway/gene-program scores.
   - Compare FLT vs GC pathway scores and clusters.
   - This is closest to the usage in the Scientific Reports paper the user linked.

2. Reference-query expiMap
   - Train a general mouse bulk reference model, likely using ARCHS4 mouse RNA-seq.
   - Map OSDR tissue bulk RNA-seq as the query using scArches surgery/query mapping.
   - This is closer to the original expiMap/scArches reference-query workflow.

Recommended order:

1. Direct OSDR expiMap first.
2. ARCHS4 reference-query expiMap second.

Reason: OSDR is bulk RNA-seq. ARCHS4 is also bulk RNA-seq and is a better reference candidate than single-cell TMS for this method. Single-cell references are possible later, but bulk query to single-cell reference is a stronger modality mismatch unless pseudobulk is introduced.

## Key Links

- expiMap advanced surgery docs:
  - https://docs.scarches.org/en/latest/expimap_surgery_pipeline_advanced.html
- expiMap reproducibility repo:
  - https://github.com/theislab/expiMap_reproducibility
- Local expiMap reproducibility clone:
  - `src/expiMap_reproducibility`
- Scientific Reports paper mentioned by user:
  - https://www.nature.com/articles/s41598-025-08649-0
- ARCHS4 download page:
  - https://maayanlab.cloud/archs4/download.html
- Direct ARCHS4 mouse gene H5:
  - https://s3.dev.maayanlab.cloud/archs4/files/mouse_gene_v2.5.h5
- NASA OSDR Biological Data API:
  - https://visualization.osdr.nasa.gov/biodata/api/
- Repo OSDR API notes:
  - `docs/osdr_api.md`
- Current expiMap results summary:
  - `docs/expimap_results.md`

ARCHS4 mouse gene file details from the download page:

- Filename: `mouse_gene_v2.5.h5`
- Date: `8-24-2024`
- Size: about 36 GB
- HEAD request returned content length `38960132574`
- SHA1 on page: `22605c9b6c4e7502b0861d4d8591ce128907c39f`

## Current Data Sources

Do not use the older raw/integrated OSDR HDF5 as the expiMap OSDR source. The direct OSDR expiMap inputs are now built from the NASA OSDR Biological Data API and its current unnormalized count tables.

API discovery script:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics
```

Discovery filters:

- `study.characteristics.organism = Mus musculus`
- `study.factor value.spaceflight = Space Flight | Ground Control`
- `file.datatype = unnormalized counts`
- bulk RNA-seq inferred from assay/file text
- single-cell and microarray-like rows excluded
- all OSDR data sources returned by the API retained

Small API provenance files:

```text
data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv
data/osdr_api/osdr_api_mouse_bulk_rnaseq_tissue_counts.tsv
data/osdr_api/osdr_api_mouse_bulk_rnaseq_tissue_accession_counts.tsv
data/osdr_api/osdr_api_mouse_bulk_rnaseq_files.tsv
data/osdr_api/osdr_api_mouse_bulk_rnaseq_summary.json
```

Downloaded per-accession API count CSVs are cached under:

```text
data/osdr_api/counts/
```

That cache is ignored by git.

Current discovery counts:

- 1631 selected FLT/GC samples
- 75 OSD accessions
- 75 count files
- 24 tissues
- liver before the 12-profile outlier filter: 125 flight, 118 ground control
- kidney: 68 flight, 67 ground control

Prepared expiMap tissue inputs:

```text
outputs/expimap_direct_osdr_liver/input/
outputs/expimap_direct_osdr_kidney/input/
```

After the liver 12-profile muscle/outlier filter:

- liver: 231 samples, 118 flight, 113 ground control, 12 accessions, 9321 genes, 1140 Reactome terms
- kidney: 135 samples, 68 flight, 67 ground control, 6 accessions, 9321 genes, 1140 Reactome terms

Generated OSDR transformations:

- `raw_counts`: primary input, raw API unnormalized counts, recommended `recon_loss=nb`
- `cpm`: sensitivity input, recommended `recon_loss=mse`
- `log1p_cpm`: sensitivity input, recommended `recon_loss=mse`

TPM/log1pTPM were not generated because the selected API count tables are unnormalized counts and no transcript-length/TPM field is used. Z-scored inputs were not generated because installed expiMap applies log-style preprocessing and expects nonnegative expression values.

## Pathway Architecture

For expiMap, the "architecture" is the pathway membership mask:

- `adata.X`: expression matrix, sample/cell by gene
- `adata.layers["counts"]`: raw count matrix
- `adata.obs`: sample/cell metadata
- `adata.var`: gene metadata
- `adata.varm["I"]`: gene by pathway mask
- `adata.uns["terms"]`: pathway names in the same order as mask columns

Use the same Reactome architecture for both direct OSDR and ARCHS4 reference-query experiments, so pathway scores are comparable.

Preferred Reactome source:

```text
https://reactome.org/download/current/ReactomePathways.txt
https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt
```

Build a native/current mouse Reactome GMT by filtering `Ensembl2Reactome_All_Levels.txt` to:

```text
species == Mus musculus
pathway stable IDs starting with R-MMU-
genes starting with ENSMUSG
```

Then join pathway IDs to names from `ReactomePathways.txt`, also filtered to `Mus musculus`. This gives a direct Reactome Mus musculus Ensembl-ID pathway architecture matching OSDR gene IDs.

Suggested generated architecture path:

```text
data/pathways/reactome_current_mouse_ensembl.gmt
```

Historical fallback/reference only:

```text
src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt
```

That fallback came from the expiMap reproducibility repo and appears to be an MSigDB Reactome v4.0 symbol set converted to mouse Ensembl IDs. It is useful for reproducing older expiMap examples, but for this NASA mouse project prefer the direct/current Reactome Mus musculus mapping.

Do not use `c2.cp.reactome.v7.5.1.symbols.gmt` unless doing a symbol-based ARCHS4-only test. If ARCHS4 uses symbols, map ARCHS4 symbols to Ensembl or build a careful shared gene table. The preferred comparable workflow is:

```text
OSDR genes ∩ ARCHS4 genes ∩ Reactome mouse Ensembl genes
```

Then use this shared gene universe/order for both direct and reference-query expiMap.

## GLARE Background

The GLARE work before this pivot:

- GLARE was reproduced/adapted for mouse liver spaceflight transcriptomics.
- Pretraining used Tabula Muris Senis FACS liver single-cell data.
- Fine-tuning used OSDR liver bulk RNA-seq.
- Per-study GLARE plus DESeq2 comparisons were more defensible than the aggregate analysis.
- GLARE clusters genes, not samples.
- PCA/UMAP plots for GLARE show genes as points and GLARE consensus gene clusters as colors.

Important GLARE output directories:

```text
outputs/glare_per_study_liver_noercc_12filter/
outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/
outputs/glare_cluster_visual_review/
```

Important GLARE interpretation:

- Stronger recurring themes: translation/ribosome, immune/MHC/stress, mitochondrial/metabolism/lipid/xenobiotic, platelet/hemostasis, RNA/DNA/cell-cycle/repair.
- Apoptosis appears in Reactome enrichment tables, but it is spread across multiple GLARE clusters and should be treated as secondary, not a headline GLARE module.
- Muscle/cytoskeleton signatures were an issue in aggregate liver GLARE; a 12-sample muscle outlier filter was used.

The 12-filter file:

```text
data/filters/aggregate_liver_12_muscle_candidate_profiles.txt
```

For expiMap, start with a liver FLT/GC dataset that removes these 12 candidate muscle outlier profiles.

## Direct OSDR expiMap Workflow

Scripts now implemented:

```text
src/nasa_mouse_glare/fetch_osdr_mouse_transcriptomics.py
src/nasa_mouse_glare/prepare_expimap_osdr_tissue.py
src/nasa_mouse_glare/train_expimap_direct.py
src/nasa_mouse_glare/analyze_expimap_pathways.py
src/nasa_mouse_glare/compare_expimap_transformations.py
```

Prepare API-derived direct OSDR tissue inputs:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.prepare_expimap_osdr_tissue --tissue liver
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.prepare_expimap_osdr_tissue --tissue kidney
```

Primary direct expiMap runs:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.train_expimap_direct \
  --input outputs/expimap_direct_osdr_liver/input/osdr_liver_flt_gc_reactome_raw_counts.h5ad \
  --output-dir outputs/expimap_direct_osdr_liver/raw_counts_nb_50epoch \
  --recon-loss nb \
  --epochs 50 \
  --hidden-layer-sizes 64

PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.train_expimap_direct \
  --input outputs/expimap_direct_osdr_kidney/input/osdr_kidney_flt_gc_reactome_raw_counts.h5ad \
  --output-dir outputs/expimap_direct_osdr_kidney/raw_counts_nb_50epoch \
  --recon-loss nb \
  --epochs 50 \
  --hidden-layer-sizes 64
```

Sensitivity runs:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.train_expimap_direct \
  --input outputs/expimap_direct_osdr_liver/input/osdr_liver_flt_gc_reactome_cpm.h5ad \
  --output-dir outputs/expimap_direct_osdr_liver/cpm_mse_50epoch \
  --recon-loss mse \
  --epochs 50 \
  --hidden-layer-sizes 64

PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.train_expimap_direct \
  --input outputs/expimap_direct_osdr_liver/input/osdr_liver_flt_gc_reactome_log1p_cpm.h5ad \
  --output-dir outputs/expimap_direct_osdr_liver/log1p_cpm_mse_50epoch \
  --recon-loss mse \
  --epochs 50 \
  --hidden-layer-sizes 64

PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.train_expimap_direct \
  --input outputs/expimap_direct_osdr_kidney/input/osdr_kidney_flt_gc_reactome_cpm.h5ad \
  --output-dir outputs/expimap_direct_osdr_kidney/cpm_mse_50epoch \
  --recon-loss mse \
  --epochs 50 \
  --hidden-layer-sizes 64

PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.train_expimap_direct \
  --input outputs/expimap_direct_osdr_kidney/input/osdr_kidney_flt_gc_reactome_log1p_cpm.h5ad \
  --output-dir outputs/expimap_direct_osdr_kidney/log1p_cpm_mse_50epoch \
  --recon-loss mse \
  --epochs 50 \
  --hidden-layer-sizes 64
```

Analysis commands:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.analyze_expimap_pathways \
  --scores outputs/expimap_direct_osdr_liver/raw_counts_nb_50epoch/pathway_scores.tsv \
  --output-dir outputs/expimap_direct_osdr_liver/raw_counts_nb_50epoch/analysis

PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.compare_expimap_transformations \
  --tissue liver \
  --tissue-dir outputs/expimap_direct_osdr_liver \
  --output-dir outputs/expimap_direct_osdr_liver/preprocessing_comparison_50epoch
```

Run the same analysis/compare commands for kidney by replacing the tissue and output paths.

Current direct outputs:

- `outputs/expimap_direct_osdr_liver/raw_counts_nb_50epoch/`
- `outputs/expimap_direct_osdr_liver/cpm_mse_50epoch/`
- `outputs/expimap_direct_osdr_liver/log1p_cpm_mse_50epoch/`
- `outputs/expimap_direct_osdr_liver/preprocessing_comparison_50epoch/`
- `outputs/expimap_direct_osdr_kidney/raw_counts_nb_50epoch/`
- `outputs/expimap_direct_osdr_kidney/cpm_mse_50epoch/`
- `outputs/expimap_direct_osdr_kidney/log1p_cpm_mse_50epoch/`
- `outputs/expimap_direct_osdr_kidney/preprocessing_comparison_50epoch/`

Each analysis directory contains FLT-vs-GC pathway tests, study-aware accession effect summaries, PCA/UMAP plots by condition and accession, and a top-pathway heatmap.

Preprocessing comparison snapshot:

- Liver direct 50-epoch runs nominate `R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION` as lower in flight; it is FDR-significant and the top aggregate term in raw counts, CPM, and log1p-CPM.
- Kidney direct 50-epoch runs have no Welch or Mann-Whitney pathway FDR below 0.10 in raw counts, CPM, or log1p-CPM.
- Direct 50-epoch effect ranks are highly correlated across raw-count NB and CPM/log1p-CPM MSE sensitivity runs. Treat normalized MSE inputs as sensitivity only, not count-likelihood primary analyses.

## ARCHS4 Reference-Query Workflow

Reference file to download:

```text
https://s3.dev.maayanlab.cloud/archs4/files/mouse_gene_v2.5.h5
```

Suggested local path:

```text
assets/archs4/mouse_gene_v2.5.h5
```

Current local status:

- File exists at `assets/archs4/mouse_gene_v2.5.h5`
- Size: 38,960,132,574 bytes
- SHA1 verified: `22605c9b6c4e7502b0861d4d8591ce128907c39f`
- `assets/` is ignored by git

Scripts now implemented:

```text
src/nasa_mouse_glare/inspect_archs4_mouse.py
src/nasa_mouse_glare/prepare_expimap_archs4_reference.py
src/nasa_mouse_glare/train_expimap_archs4_reference.py
src/nasa_mouse_glare/map_expimap_osdr_query.py
```

ARCHS4 inspection output:

```text
data/archs4/archs4_mouse_tissue_summary.tsv
```

Usable nonleakage bulk-like tissue candidates from the inspection:

- liver: 8970 samples
- spleen: 6289 samples
- lung: 5674 samples
- skin: 2593 samples
- kidney: 2464 samples
- skeletal_muscle: 1412 samples

Current bounded reference-query outputs:

- `outputs/expimap_archs4_reference_osdr_query_liver/reference_input_1000/`
- `outputs/expimap_archs4_reference_osdr_query_liver/reference_nb_1000_50epoch/`
- `outputs/expimap_archs4_reference_osdr_query_liver/query_nb_1000ref_50epoch/`
- `outputs/expimap_archs4_reference_osdr_query_kidney/reference_input_1000/`
- `outputs/expimap_archs4_reference_osdr_query_kidney/reference_nb_1000_50epoch/`
- `outputs/expimap_archs4_reference_osdr_query_kidney/query_nb_1000ref_50epoch/`

The current bounded reference inputs use 1000 ARCHS4 samples per tissue, 50 reference-training epochs, and 50 query-mapping epochs. Query mapping preserved the Reactome term structure and dropped 2 OSDR genes absent from the reference, retaining 9319 shared genes and 1140 pathways. Neither liver nor kidney reference-query analysis has an aggregate Welch pathway FDR below 0.10 in these bounded runs. Older 100-sample/1-epoch smoke outputs remain under `reference_input_smoke/`, `reference_smoke_nb/`, and `query_smoke_nb/` as mechanics-validation artifacts.

Reference-query preprocessing is raw-count NB only in the current workflow. CPM/log1p-CPM comparisons are direct-workflow sensitivity analyses, not reference-query surgery runs.

Direct/reference agreement note:

- Direct liver raw-count NB nominates `R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION` as lower in flight with Welch FDR about 0.012.
- The same term is also lower in flight after bounded ARCHS4 reference-query mapping, but its reference-query Welch FDR is about 0.97.
- Treat the direct liver signal as preprocessing-stable but not reference-query-confirmed.
- Kidney has no aggregate FDR-significant direct or reference-query pathway signal in the current runs.

## Comparison to GLARE

GLARE output:

- Gene-level latent representations.
- Gene consensus modules.
- Interpretation via DGEA overlap and pathway enrichment.

expiMap output:

- Sample/cell-level pathway scores.
- Each latent dimension is tied to a pathway/gene program through the mask.
- More directly suited to asking: "Which pathways shift in flight?"

This is why expiMap may be more interpretable than GLARE for the current biological question.

## Decisions Already Made

- Use ARCHS4 as likely general mouse bulk reference, not TMS single-cell, for the first reference-query expiMap.
- Use the same Reactome architecture for direct OSDR and reference-query experiments.
- Use NASA OSDR Biological Data API metadata and count tables for OSDR expiMap inputs.
- Prepare both liver and kidney FLT vs GC direct OSDR inputs.
- Use the 12 muscle outlier filter for liver.
- Do not treat MOBER-corrected values as primary count input for expiMap.
- Use direct raw-count NB expiMap as the primary direct OSDR run; CPM/log1p-CPM MSE runs are sensitivity checks only.
- Use `id.accession` as the expiMap condition/batch key for direct OSDR runs unless a later diagnostic shows a better covariate.

## Open Questions / Next Full Runs

- Should direct OSDR expiMap use all API-discovered tissue FLT/GC studies or only stronger primary studies?
  - Previous stronger liver studies included `OSD-379`, `OSD-245`, and `OSD-463`.
  - Prior per-study GLARE/DGEA also included `OSD-168`, `OSD-48`, and `OSD-137`.
  - Earlier aggregate liver work also considered `OSD-242`, `OSD-173`, `OSD-47`, and `OSD-686`.
  - Current API metadata contains exact FLT/GC counts in `data/osdr_api/osdr_api_mouse_bulk_rnaseq_tissue_accession_counts.tsv`.

- Should ARCHS4 be filtered to liver/general tissue subsets?
  - Current bounded runs use tissue-specific ARCHS4 subsets.
  - Follow-up reference models should use larger tissue-specific subsets or all available nonleakage tissue samples, then compare to a broader nonleakage reference if compute allows.

- How many epochs and what architecture should full runs use?
  - Current direct runs use 50 epochs and hidden size 64.
  - Full follow-up runs should still monitor convergence and preserve the same API-derived input manifests.

## Minimal Next Commands

From repo root:

```bash
cd /media/volume/mouse/nasa/nasa-mouse

/home/exouser/miniforge3/envs/nasa-mouse/bin/python - <<'PY'
import scarches as sca
print(hasattr(sca.models, "EXPIMAP"))
PY
```

Re-run API discovery and tissue input preparation:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.prepare_expimap_osdr_tissue --tissue liver
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.prepare_expimap_osdr_tissue --tissue kidney
```

Before a full training run, inspect:

- input manifests under `outputs/expimap_direct_osdr_{tissue}/input/input_manifest.json`
- number of samples retained
- FLT/GC counts
- samples per OSD accession
- number of genes retained after Reactome filtering
- number of Reactome terms retained
- PCA by OSD accession after the direct and reference-query runs

Validation command:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m py_compile \
  src/nasa_mouse_glare/build_reactome_mouse_gmt.py \
  src/nasa_mouse_glare/fetch_osdr_mouse_transcriptomics.py \
  src/nasa_mouse_glare/prepare_expimap_osdr_tissue.py \
  src/nasa_mouse_glare/train_expimap_direct.py \
  src/nasa_mouse_glare/analyze_expimap_pathways.py \
  src/nasa_mouse_glare/inspect_archs4_mouse.py \
  src/nasa_mouse_glare/prepare_expimap_archs4_reference.py \
  src/nasa_mouse_glare/map_expimap_osdr_query.py \
  src/nasa_mouse_glare/compare_expimap_transformations.py
```
