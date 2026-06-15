# Generated GLARE Inputs

This directory contains local preprocessing outputs and should not store large
matrix artifacts in git. The generated matrices are ignored by `.gitignore`.

Generated on 2026-06-12 in the `nasa` conda environment.

## Inputs

- TMS FACS `.h5ad`: `assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad`
- OSDR mouse RNA-seq HDF5: `assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5`

## Commands

Run these from the repository root:

```bash
cd path/to/nasa-mouse
```

```bash
PYTHONPATH=src /opt/anaconda3/envs/nasa/bin/python -m nasa_mouse_glare.tms \
  --input assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad \
  --output-prefix data/processed/tms_facs_3552_cells \
  --matrix-source X \
  --max-cells 3552

PYTHONPATH=src /opt/anaconda3/envs/nasa/bin/python -m nasa_mouse_glare.osdr prep-h5 \
  --input assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5 \
  --output-prefix data/processed/osdr_mouse_bulk \
  --matrix-key /data/expression \
  --gene-key /meta/genes/ensembl_gene \
  --sample-key "/meta/info/id.sample name"

PYTHONPATH=src /opt/anaconda3/envs/nasa/bin/python -m nasa_mouse_glare.align \
  --pretrain data/processed/tms_facs_3552_cells.manifest.json \
  --target data/processed/osdr_mouse_bulk.manifest.json \
  --output-prefix data/processed/tms_facs_osdr_aligned

PYTHONPATH=src /opt/anaconda3/envs/nasa/bin/python -m nasa_mouse_glare.export mtx \
  --bundle data/processed/tms_facs_osdr_aligned.pretrain.manifest.json \
  --output data/glare_inputs/tms_facs_pretrain.mtx

PYTHONPATH=src /opt/anaconda3/envs/nasa/bin/python -m nasa_mouse_glare.export csv \
  --bundle data/processed/tms_facs_osdr_aligned.target.manifest.json \
  --output data/glare_inputs/osdr_finetune.csv
```

## Output Shapes

- `data/processed/tms_facs_3552_cells`: `21025 genes x 3552 cells`
- `data/processed/osdr_mouse_bulk`: `53511 genes x 3315 samples`
- `data/processed/tms_facs_osdr_aligned.pretrain`: `21010 genes x 3552 cells`
- `data/processed/tms_facs_osdr_aligned.target`: `21010 genes x 3315 samples`
- `data/glare_inputs/tms_facs_pretrain.mtx`: `21010 genes x 3552 cells`
- `data/glare_inputs/osdr_finetune.csv`: `21010 genes x 3315 samples`

The GLARE source is vendored directly in `src/glare` with the `hpt.py` runtime
fix applied: direct script execution works, MatrixMarket input uses SciPy's
`.toarray()` API, fine-tuning reuses the pretraining architecture, and final
representation extraction applies the fine-tuning adapter.
