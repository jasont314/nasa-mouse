# expiMap / scArches Handoff

Date: 2026-06-25

This file summarizes the current NASA mouse spaceflight transcriptomics context for handing the work to another agent. The current pivot is from GLARE to expiMap/scArches.

## Repository State

- Repo root: `/Users/jasontrinh/Desktop/Code/Berkeley/nasa/nasa-mouse`
- Main branch: `main`
- Latest relevant commit: `8df9641 Add scarches expiMap dependencies`
- Conda env: `/opt/anaconda3/envs/nasa`
- Environment file: `environment.yml`
- Project Python requirements: `requirements-nasa-mouse-glare.txt`

The dependency update has been committed and pushed. `pip check` passed after installation.

## Important Package Caveat

There are two scArches code sources in this repo/environment:

1. Installed package in the `nasa` conda env:
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
/opt/anaconda3/envs/nasa/bin/python - <<'PY'
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

## User Goal

The user wants to test expiMap as a new pathway-aware model for mouse spaceflight liver transcriptomics.

Two planned approaches:

1. Direct OSDR expiMap
   - Train expiMap directly on OSDR liver bulk RNA-seq samples.
   - Extract pathway/gene-program scores.
   - Compare FLT vs GC pathway scores and clusters.
   - This is closest to the usage in the Scientific Reports paper the user linked.

2. Reference-query expiMap
   - Train a general mouse bulk reference model, likely using ARCHS4 mouse RNA-seq.
   - Map OSDR liver bulk RNA-seq as the query using scArches surgery/query mapping.
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

ARCHS4 mouse gene file details from the download page:

- Filename: `mouse_gene_v2.5.h5`
- Date: `8-24-2024`
- Size: about 36 GB
- HEAD request returned content length `38960132574`
- SHA1 on page: `22605c9b6c4e7502b0861d4d8591ce128907c39f`

## Existing Data

Raw integrated OSDR mouse RNA-seq HDF5:

```text
assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5
```

Relevant structure observed:

```text
data/expression: shape (53511, 3315), dtype uint32
meta/genes/ensembl_gene
meta/genes/symbol
meta/info/id.accession
...
```

The raw expression matrix is integer count-like data and is suitable as the starting point for expiMap. Do not use MOBER-corrected values as the primary expiMap count input.

Processed OSDR liver files already exist:

```text
data/processed/osdr_mouse_bulk_liver.matrix.npz
data/processed/osdr_mouse_bulk_liver.genes.tsv
data/processed/osdr_mouse_bulk_liver.profiles.tsv
data/processed/osdr_mouse_bulk_liver.profile_metadata.tsv
data/processed/osdr_mouse_bulk_liver.manifest.json
```

Processed liver matrix:

```text
shape: (53511, 628)
orientation: genes_x_profiles
dtype in NPZ check: float32
```

Gene IDs in the processed OSDR liver files are mouse Ensembl IDs:

```text
ENSMUSG00000000001
ENSMUSG00000000003
ENSMUSG00000000028
...
```

## Pathway Architecture

For expiMap, the "architecture" is the pathway membership mask:

- `adata.X`: expression matrix, sample/cell by gene
- `adata.layers["counts"]`: raw count matrix
- `adata.obs`: sample/cell metadata
- `adata.var`: gene metadata
- `adata.varm["I"]`: gene by pathway mask
- `adata.uns["terms"]`: pathway names in the same order as mask columns

Use the same Reactome architecture for both direct OSDR and ARCHS4 reference-query experiments, so pathway scores are comparable.

Recommended Reactome GMT:

```text
src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt
```

Reason:

- It uses mouse Ensembl IDs (`ENSMUSG...`).
- That matches current OSDR gene IDs.
- It avoids symbol mapping ambiguity for OSDR.

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

## Direct OSDR expiMap Plan

Target output directory:

```text
outputs/expimap_osdr_liver_12filter/
```

First script to write, suggested location:

```text
src/nasa_mouse_glare/prepare_expimap_osdr_liver.py
```

Suggested outputs:

```text
outputs/expimap_osdr_liver_12filter/input/osdr_liver_flt_gc_reactome.h5ad
outputs/expimap_osdr_liver_12filter/input/reactome_terms.tsv
outputs/expimap_osdr_liver_12filter/input/gene_universe.tsv
outputs/expimap_osdr_liver_12filter/input/profile_metadata.tsv
outputs/expimap_osdr_liver_12filter/input/input_manifest.json
```

Prep steps:

1. Load OSDR liver counts from processed files or raw HDF5.
2. Convert to sample by gene AnnData.
3. Add sample metadata from `data/processed/osdr_mouse_bulk_liver.profile_metadata.tsv`.
4. Keep only:
   - `condition_inferred == flight`
   - `condition_inferred == ground_control`
5. Remove profiles listed in:
   - `data/filters/aggregate_liver_12_muscle_candidate_profiles.txt`
6. Add Reactome annotations using:
   - `src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt`
7. Filter genes to those in at least one retained pathway.
8. Filter pathways by minimum gene count, likely `min_genes=12`.
9. Store raw counts in both `adata.X` and `adata.layers["counts"]`.
10. Store study accession/batch in `adata.obs["id.accession"]`.

Training script, suggested location:

```text
src/nasa_mouse_glare/train_expimap_osdr_liver.py
```

Initial model settings to mirror docs/reproducibility:

- model: `sca.models.EXPIMAP`
- count likelihood: negative binomial
- condition/batch key: likely `id.accession`
- pathways/mask: `adata.varm["I"]`
- save pathway latent scores to `adata.obsm`

Primary analyses after training:

1. Pathway score table:
   - rows = samples
   - columns = Reactome pathways
   - metadata = FLT/GC, OSD accession, mission/study, sex/strain if available

2. FLT vs GC pathway comparison:
   - Do not run one naive test alone.
   - Include study-aware testing or at least per-study sensitivity.
   - Compare direct aggregate result against per-study effects.

3. Plots:
   - PCA/UMAP of samples using expiMap pathway scores.
   - Color by FLT/GC.
   - Color by OSD accession to inspect batch/study effect.
   - Heatmap of top FLT-vs-GC pathway scores.
   - Dotplot or forest plot of pathway effect by study.

## ARCHS4 Reference-Query Plan

Reference file to download:

```text
https://s3.dev.maayanlab.cloud/archs4/files/mouse_gene_v2.5.h5
```

Suggested local path:

```text
assets/archs4/mouse_gene_v2.5.h5
```

Because the file is about 36 GB, add/confirm `.gitignore` coverage before downloading. `assets/` has generally been ignored in this repo, but confirm before downloading.

Key tasks before training:

1. Inspect ARCHS4 HDF5 structure.
2. Extract mouse gene expression matrix and sample metadata.
3. Identify gene IDs/symbols and map to the same mouse Ensembl IDs used by OSDR.
4. Filter out OSDR/GeneLab/spaceflight-like samples to avoid leakage.
5. Decide whether to train on all ARCHS4 mouse samples or a biologically focused subset.
6. Build the same Reactome architecture/gene order as the OSDR input.
7. Train expiMap reference model.
8. Map OSDR liver FLT/GC samples as query.
9. Compare query pathway scores by FLT/GC and by study.

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
- Start with OSDR liver FLT vs GC only.
- Use the 12 muscle outlier filter.
- Do not treat MOBER-corrected values as primary count input for expiMap.
- Direct OSDR expiMap should be implemented first because it has fewer moving parts.

## Open Questions

- Should direct OSDR expiMap use all liver FLT/GC studies or only stronger primary studies?
  - Previous stronger liver studies included `OSD-379`, `OSD-245`, and `OSD-463`.
  - Prior per-study GLARE/DGEA also included `OSD-168`, `OSD-48`, and `OSD-137`.
  - Earlier aggregate liver work also considered `OSD-242`, `OSD-173`, `OSD-47`, and `OSD-686`.
  - Re-derive exact FLT/GC counts from metadata before training.

- Should `id.accession` be used as the expiMap condition/batch key, or should model training use a different batch covariate?
  - Start with `id.accession`.
  - Then inspect whether pathway-score structure is dominated by study.

- Should ARCHS4 be filtered to liver/general tissue subsets?
  - Start broad only if feasible.
  - A focused liver/general tissue subset may be more tractable and less noisy.

## Minimal Next Commands

From repo root:

```bash
cd /Users/jasontrinh/Desktop/Code/Berkeley/nasa/nasa-mouse

/opt/anaconda3/envs/nasa/bin/python - <<'PY'
import scarches as sca
print(hasattr(sca.models, "EXPIMAP"))
PY
```

Then implement the OSDR input prep script and write:

```text
outputs/expimap_osdr_liver_12filter/input/osdr_liver_flt_gc_reactome.h5ad
```

Before a full training run, inspect:

- number of samples retained
- FLT/GC counts
- samples per OSD accession
- number of genes retained after Reactome filtering
- number of Reactome terms retained
- top/largest Reactome terms
- whether any duplicate gene IDs need collapsing

