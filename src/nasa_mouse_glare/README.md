# nasa_mouse_glare

Project-specific GLARE adaptation for mouse spaceflight data.

The GLARE code is vendored directly in `src/glare` so project-specific runtime
fixes can be edited and committed in this repository. This package is where
mouse-specific data preparation and model changes live.

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

All commands below assume they are run from the repository root:

```bash
cd path/to/nasa-mouse
```

```bash
conda env create -f environment.yml
conda activate nasa
export PYTHONPATH=src
```

For the existing local environment, dependencies were installed with:

```bash
conda create -y -n nasa python=3.11
conda run -n nasa python -m pip install -r requirements-nasa-mouse-glare.txt
conda activate nasa
export PYTHONPATH=src
```

## Download links

```bash
PYTHONPATH=src python -m nasa_mouse_glare.downloads links
PYTHONPATH=src python -m nasa_mouse_glare.downloads download --kind facs --output-dir assets/tms
PYTHONPATH=src python -m nasa_mouse_glare.downloads download --kind droplet --output-dir assets/tms
```

## Prepare TMS pretraining matrix

Strict GLARE-compatible direct-cell prep, sampled to the same cell count as the
paper's original `E-CURD-5` pretraining set:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.tms \
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
PYTHONPATH=src python -m nasa_mouse_glare.osdr inspect-h5 \
  --input assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5

PYTHONPATH=src python -m nasa_mouse_glare.osdr prep-h5 \
  --input assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5 \
  --output-prefix data/processed/osdr_mouse_bulk \
  --matrix-key /data/expression \
  --gene-key /meta/genes/ensembl_gene \
  --sample-key "/meta/info/id.sample name"
```

## Align genes

```bash
PYTHONPATH=src python -m nasa_mouse_glare.align \
  --pretrain data/processed/tms_facs_3552_cells.manifest.json \
  --target data/processed/osdr_mouse_bulk.manifest.json \
  --output-prefix data/processed/tms_facs_osdr_aligned
```

## Export for upstream GLARE `hpt.py`

```bash
PYTHONPATH=src python -m nasa_mouse_glare.export mtx \
  --bundle data/processed/tms_facs_osdr_aligned.pretrain.manifest.json \
  --output data/glare_inputs/tms_facs_pretrain.mtx

PYTHONPATH=src python -m nasa_mouse_glare.export csv \
  --bundle data/processed/tms_facs_osdr_aligned.target.manifest.json \
  --output data/glare_inputs/osdr_finetune.csv
```

That produces the same broad input contract as GLARE:

```bash
python src/glare/Manuscript_Code/glare/codes/hpt.py \
  --data1 data/glare_inputs/osdr_finetune.csv \
  --data2 data/glare_inputs/tms_facs_pretrain.mtx \
  --log-every-epochs 1 \
  --output-dir outputs/glare_hpt_tms_facs_osdr
```

The vendored GLARE copy includes the `hpt.py` runtime fix for direct script
execution, MatrixMarket loading, architecture-compatible fine-tuning, and
adapter-aware representation extraction. During HPT, each completed config is
saved to `hpt_config_results.csv` and `hpt_config_results.json` inside
`--output-dir`; final weights, the representation array, and `hpt_summary.json`
are saved there too.

## Reproduce Original GLARE Pretraining

Download the Arabidopsis single-cell normalized MatrixMarket file used by
GLARE:

```bash
mkdir -p assets/glare_original

curl -L --fail --show-error \
  --output assets/glare_original/E-CURD-5.aggregated_filtered_normalised_counts.mtx.gz \
  https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-CURD-5/E-CURD-5.aggregated_filtered_normalised_counts.mtx.gz

gunzip -kf assets/glare_original/E-CURD-5.aggregated_filtered_normalised_counts.mtx.gz
```

Run the released GLARE pretraining config:

```bash
python src/nasa_mouse_glare/reproduce_glare_pretrain.py \
  --input assets/glare_original/E-CURD-5.aggregated_filtered_normalised_counts.mtx \
  --output-dir outputs/glare_original_pretrain_config5 \
  --epochs 30
```

This uses the released GLARE SAE setup: `[128, 64, 32, 16]`, LayerNorm, ELU,
Adam `lr=1e-3`, weight decay `1e-4`, sparsity penalty `1e-5`, and batch size
`16`.

## Train

```bash
PYTHONPATH=src python -m nasa_mouse_glare.model pretrain \
  --matrix data/processed/tms_facs_osdr_aligned.pretrain.manifest.json \
  --output models/tms_facs_pretrained.pth \
  --epochs 30

# Use the printed pretrain_input_dim value from pretrain.
PYTHONPATH=src python -m nasa_mouse_glare.model finetune \
  --target-matrix data/processed/tms_facs_osdr_aligned.target.manifest.json \
  --pretrained models/tms_facs_pretrained.pth \
  --pretrain-input-dim <printed_pretrain_input_dim> \
  --output models/tms_facs_osdr_finetuned.pth \
  --epochs 30
```
