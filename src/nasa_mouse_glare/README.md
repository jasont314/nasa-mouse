# nasa_mouse_glare

Project-specific GLARE adaptation for mouse spaceflight data.

The upstream GLARE code is kept as a submodule in `src/glare`. This package is
where mouse-specific data preparation and model changes live.

## Data model

GLARE pretraining used an EBI MatrixMarket archive:

- `E-CURD-5.aggregated_filtered_normalised_counts.mtx`
- shape: `20,347 genes x 3,552 cells`
- format: `.mtx`, not `.h5` or `.h5ad`

The format does not matter biologically. GLARE used `.mtx`; TMS is `.h5ad`.
What matters is converting TMS `.h5ad` to the same orientation:
`genes x cells`.

Use direct cells and export a MatrixMarket `.mtx` for `hpt.py`.

For TMS, prefer processed `.h5ad` over raw FASTQ:

- FACS/Smart-seq2: higher gene coverage, closer to bulk RNA-seq.
- Droplet/10x: more cells, more sparse, useful as a second pretraining run.

## Setup

```bash
pip install -r requirements-nasa-mouse-glare.txt
```

## Download links

```bash
python -m nasa_mouse_glare.downloads links
python -m nasa_mouse_glare.downloads download --kind facs --output-dir assets/tms
python -m nasa_mouse_glare.downloads download --kind droplet --output-dir assets/tms
```

## Prepare TMS pretraining matrix

Strict GLARE-compatible direct-cell prep, sampled to the same cell count as the
paper's original `E-CURD-5` pretraining set:

```bash
python -m nasa_mouse_glare.tms \
  --input assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad \
  --output-prefix data/processed/tms_facs_3552_cells \
  --matrix-source X \
  --max-cells 3552
```

Full direct-cell TMS can be exported too, but upstream GLARE densifies the
matrix, so full FACS/droplet is likely too large without changing training.

## Prepare OSDR fine-tuning matrix

Your current OSDR HDF5 has `/data/expression` with shape `53511 x 3315`.

```bash
python -m nasa_mouse_glare.osdr inspect-h5 \
  --input assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5

python -m nasa_mouse_glare.osdr prep-h5 \
  --input assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5 \
  --output-prefix data/processed/osdr_mouse_bulk \
  --matrix-key /data/expression \
  --gene-key /meta/genes/ensembl_gene \
  --sample-key "/meta/info/id.sample name"
```

## Align genes

```bash
python -m nasa_mouse_glare.align \
  --pretrain data/processed/tms_facs_3552_cells.manifest.json \
  --target data/processed/osdr_mouse_bulk.manifest.json \
  --output-prefix data/processed/tms_facs_osdr_aligned
```

## Export for upstream GLARE `hpt.py`

```bash
python -m nasa_mouse_glare.export mtx \
  --bundle data/processed/tms_facs_osdr_aligned.pretrain.manifest.json \
  --output data/glare_inputs/tms_facs_pretrain.mtx

python -m nasa_mouse_glare.export csv \
  --bundle data/processed/tms_facs_osdr_aligned.target.manifest.json \
  --output data/glare_inputs/osdr_finetune.csv
```

That produces the same broad input contract as GLARE:

```bash
cd src/glare/Manuscript_Code/glare/codes
python hpt.py \
  --data1 ../../../../../data/glare_inputs/osdr_finetune.csv \
  --data2 ../../../../../data/glare_inputs/tms_facs_pretrain.mtx
```

## Train

```bash
python -m nasa_mouse_glare.model pretrain \
  --matrix data/processed/tms_facs_osdr_aligned.pretrain.manifest.json \
  --output models/tms_facs_pretrained.pth \
  --epochs 30

# Use the printed pretrain_input_dim value from pretrain.
python -m nasa_mouse_glare.model finetune \
  --target-matrix data/processed/tms_facs_osdr_aligned.target.manifest.json \
  --pretrained models/tms_facs_pretrained.pth \
  --pretrain-input-dim <printed_pretrain_input_dim> \
  --output models/tms_facs_osdr_finetuned.pth \
  --epochs 30
```
