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

To resume after configs `1` through `249` completed in the single-cell HPT
stage, copy the prior `hpt_config_results.json` into the same `--output-dir`
and continue from config `250`:

```bash
python src/glare/Manuscript_Code/glare/codes/hpt.py \
  --data1 data/glare_inputs/osdr_finetune.csv \
  --data2 data/glare_inputs/tms_facs_pretrain.mtx \
  --log-every-epochs 1 \
  --output-dir outputs/glare_hpt_tms_facs_osdr \
  --resume \
  --start-config 250
```

The JSON log is required for resume because it restores the best hyperparameter
configuration from the completed configs. If the CSV is missing, it is rebuilt
from JSON automatically.

The single-cell HPT stage has `486` configs. The OSDR fine-tuning HPT stage has
`27` configs after the pretraining architecture is fixed; resume that stage with
`--osdr-start-config <N>`.

The vendored GLARE copy includes the `hpt.py` runtime fix for direct script
execution, MatrixMarket loading, architecture-compatible fine-tuning, and
adapter-aware representation extraction. During HPT, each completed config is
saved to `hpt_config_results.csv` and `hpt_config_results.json` inside
`--output-dir`; final weights, the representation array, and `hpt_summary.json`
are saved there too.

## Post Fine-Tuning Analysis

The current GLARE-compatible OSDR CSV is `genes x samples`, so
`FTSAE_representation.npy` is expected to contain gene-level latent vectors.
After `hpt.py` finishes, summarize and cluster those fine-tuned gene
representations:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.post_finetune \
  --representation outputs/glare_hpt_tms_facs_osdr/FTSAE_representation.npy \
  --target-manifest data/processed/tms_facs_osdr_aligned.target.manifest.json \
  --osdr assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5 \
  --output-dir outputs/glare_hpt_tms_facs_osdr/post_finetune
```

Key outputs:

- `gene_latent.tsv`: gene IDs, latent dimensions, and gene cluster labels.
- `gene_pca.tsv`: PCA coordinates for plotting the gene latent space.
- `gene_cluster_summary.tsv`: cluster sizes and PCA centroids.
- `gene_cluster_expression_by_*.tsv`: cluster-level mean expression grouped by
  OSDR metadata, including inferred flight/ground labels where available.
- `post_finetune_summary.json`: input paths, shapes, clustering settings, and
  output file paths.

Run GLARE-style ensemble clustering with GMM, HDBSCAN, Spectral clustering, and
a consensus step:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.ensemble_clustering \
  --representation outputs/glare_hpt_tms_facs_osdr/FTSAE_representation.npy \
  --gene-latent outputs/glare_hpt_tms_facs_osdr/post_finetune/gene_latent.tsv \
  --gene-pca outputs/glare_hpt_tms_facs_osdr/post_finetune/gene_pca.tsv \
  --output-dir outputs/glare_hpt_tms_facs_osdr/post_finetune/ensemble_clustering
```

This writes `ensemble_clusters.tsv`, a consensus `gene_clusters.tsv`,
cluster/algorithm summaries, metrics, and `ensemble_pca_by_consensus.png`.
The mouse run uses diagonal-covariance GMM for numerical stability. HDBSCAN
keeps GLARE's `min_cluster_size=60` but uses `min_samples=10`; `30` labeled
all genes as noise on the mouse FTSAE latent space.

Download official OSDR ISA metadata and extract sample/accession tissues:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.osdr_tissues \
  --download \
  --metadata-dir assets/osdr_metadata \
  --profile-metadata outputs/glare_hpt_tms_facs_osdr/post_finetune/profile_metadata.tsv \
  --output-dir outputs/glare_hpt_tms_facs_osdr/post_finetune/osdr_tissues
```

This uses OSDR `Characteristics[Material Type]` from the official sample
tables, including paginated studies. When Material Type is generic, it uses
the study's explicit `Characteristics[Tissue Type]` or `Factor Value[Tissue]`
field. It writes `osdr_sample_tissues.tsv`, `osdr_accession_tissues.tsv`,
`osdr_tissue_accessions.tsv`, an official-versus-inferred validation table,
and a JSON summary. The downloaded JSON cache remains ignored under `assets/`.

Recompute expression summaries and run study/tissue-stratified significance
tests for the ensemble consensus clusters:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.cluster_stratified_analysis \
  --gene-clusters outputs/glare_hpt_tms_facs_osdr/post_finetune/ensemble_clustering/gene_clusters.tsv \
  --target-manifest data/processed/tms_facs_osdr_aligned.target.manifest.json \
  --profile-metadata outputs/glare_hpt_tms_facs_osdr/post_finetune/profile_metadata.tsv \
  --official-tissues outputs/glare_hpt_tms_facs_osdr/post_finetune/osdr_tissues/osdr_sample_tissues.tsv \
  --output-dir outputs/glare_hpt_tms_facs_osdr/post_finetune/ensemble_analysis
```

The OSDR HDF5 does not include its ISA tissue fields, so this step joins the
official tissue table generated above. Sample-name inference is retained only
as a fallback. Flight-vs-ground effects are computed within OSDR accessions
using log2 expression ratios, tested across paired accessions with a two-sided
Wilcoxon signed-rank test, and corrected with Benjamini-Hochberg FDR.

Run enrichment for all 15 ensemble consensus clusters:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.cluster_enrichment \
  --post-dir outputs/glare_hpt_tms_facs_osdr/post_finetune/ensemble_analysis \
  --output-dir outputs/glare_hpt_tms_facs_osdr/post_finetune/ensemble_analysis/enrichment \
  --clusters 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
```

Run enrichment and driver summaries for the strongest flight-shifted clusters:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.cluster_enrichment \
  --post-dir outputs/glare_hpt_tms_facs_osdr/post_finetune \
  --output-dir outputs/glare_hpt_tms_facs_osdr/post_finetune/enrichment \
  --clusters 13 10 8 0 6
```

This writes cluster gene lists, Reactome/Panglao over-representation results,
and accession/condition/source driver summaries under
`outputs/glare_hpt_tms_facs_osdr/post_finetune/enrichment`.

Run GLARE-style representation evaluation and flight-vs-ground verification:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.glare_evaluation \
  --representation outputs/glare_hpt_tms_facs_osdr/FTSAE_representation.npy \
  --target-manifest data/processed/tms_facs_osdr_aligned.target.manifest.json \
  --post-dir outputs/glare_hpt_tms_facs_osdr/post_finetune \
  --output-dir outputs/glare_hpt_tms_facs_osdr/post_finetune/evaluation
```

This mirrors GLARE's `evaluation.py` metrics where they apply to the mouse
gene-level output: KMeans silhouette, KNN predictability of KMeans labels, and
trustworthiness against the original expression space. It also runs a
verification study for `condition_inferred=flight` versus
`condition_inferred=ground_control`, including both random CV and
accession-grouped CV. GLARE's SHAP post-pipeline used XGBoost/SHAP; those are
optional dependencies here, so the current evaluation records SHAP as not run
and writes linear classifier coefficients instead.

## Liver-Only Fixed-Config Run

Use all TMS FACS liver cells for pretraining and only officially annotated
OSDR liver profiles for fine-tuning. The training command reuses the previous
winning configurations from `outputs/glare_hpt_tms_facs_osdr/hpt_summary.json`
and does not run either hyperparameter sweep:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.tms \
  --input assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad \
  --output-prefix data/processed/tms_facs_liver \
  --matrix-source X \
  --filter tissue=liver

PYTHONPATH=src python -m nasa_mouse_glare.subset_bundle \
  --bundle data/processed/osdr_mouse_bulk.manifest.json \
  --metadata outputs/glare_hpt_tms_facs_osdr/post_finetune/osdr_tissues/osdr_sample_tissues.tsv \
  --filter-column tissue_final \
  --filter-value liver \
  --output-prefix data/processed/osdr_mouse_bulk_liver

PYTHONPATH=src python -m nasa_mouse_glare.align \
  --pretrain data/processed/tms_facs_liver.manifest.json \
  --target data/processed/osdr_mouse_bulk_liver.manifest.json \
  --output-prefix data/processed/tms_facs_liver_osdr_liver_aligned

PYTHONPATH=src python -m nasa_mouse_glare.export mtx \
  --bundle data/processed/tms_facs_liver_osdr_liver_aligned.pretrain.manifest.json \
  --output data/glare_inputs/tms_facs_liver_pretrain.mtx

PYTHONPATH=src python -m nasa_mouse_glare.export csv \
  --bundle data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json \
  --output data/glare_inputs/osdr_liver_finetune.csv

python src/glare/Manuscript_Code/glare/codes/hpt.py \
  --data1 data/glare_inputs/osdr_liver_finetune.csv \
  --data2 data/glare_inputs/tms_facs_liver_pretrain.mtx \
  --reuse-best-configs-from outputs/glare_hpt_tms_facs_osdr/hpt_summary.json \
  --log-every-epochs 1 \
  --num-workers 0 \
  --output-dir outputs/glare_fixed_tms_facs_liver_osdr_liver
```

All liver-only training and analysis outputs use
`outputs/glare_fixed_tms_facs_liver_osdr_liver`, leaving the cross-tissue run
unchanged.

Run the liver-only post-training analyses:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.post_finetune \
  --representation outputs/glare_fixed_tms_facs_liver_osdr_liver/FTSAE_representation.npy \
  --target-manifest data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json \
  --osdr assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5 \
  --output-dir outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune

PYTHONPATH=src python -m nasa_mouse_glare.ensemble_clustering \
  --representation outputs/glare_fixed_tms_facs_liver_osdr_liver/FTSAE_representation.npy \
  --gene-latent outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/gene_latent.tsv \
  --gene-pca outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/gene_pca.tsv \
  --hdbscan-min-cluster-size 100 \
  --hdbscan-min-samples 1 \
  --output-dir outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/ensemble_clustering

PYTHONPATH=src python -m nasa_mouse_glare.osdr_tissues \
  --metadata-dir assets/osdr_metadata \
  --profile-metadata outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/profile_metadata.tsv \
  --output-dir outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/osdr_tissues

PYTHONPATH=src python -m nasa_mouse_glare.cluster_stratified_analysis \
  --gene-clusters outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/ensemble_clustering/gene_clusters.tsv \
  --target-manifest data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json \
  --profile-metadata outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/profile_metadata.tsv \
  --official-tissues outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/osdr_tissues/osdr_sample_tissues.tsv \
  --output-dir outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/ensemble_analysis

PYTHONPATH=src python -m nasa_mouse_glare.cluster_enrichment \
  --post-dir outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/ensemble_analysis \
  --target-manifest data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json \
  --output-dir outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/ensemble_analysis/enrichment \
  --clusters 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14

PYTHONPATH=src python -m nasa_mouse_glare.glare_evaluation \
  --representation outputs/glare_fixed_tms_facs_liver_osdr_liver/FTSAE_representation.npy \
  --target-manifest data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json \
  --post-dir outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune \
  --output-dir outputs/glare_fixed_tms_facs_liver_osdr_liver/post_finetune/evaluation
```

The liver latent space required a lower-density HDBSCAN setting than the
cross-tissue run. With `min_cluster_size=60` and `min_samples=10`, every gene
was labeled noise. The documented liver setting produces two density clusters
and 3,695 noise genes, retaining HDBSCAN as an active ensemble member.

## Paper-Style Controlled OSD-379 Run

This run addresses the main differences from the published GLARE analysis:

- one controlled liver study (`OSD-379`, RR-8), not 17 pooled studies;
- official NASA normalized counts, not the integrated raw-count matrix;
- GLARE's released `[128, 64, 32, 16]` SAE;
- separate FLT and GC fine-tuning and representations;
- melted-data XGBoost verification and SHAP;
- GMM/HDBSCAN/Spectral partitions with true co-association EAC and
  average-linkage consensus;
- official NASA DESeq2 contrasts, DEG proportions, and Metascape upload files.

OSD-379 provides 35 flight and 35 ground-control liver profiles balanced over
young/old and ISS-terminal/live-animal-return collection strata. Download its
official normalized counts and differential-expression table:

```bash
curl -L --fail --show-error \
  --output assets/osdr/GLDS-379_rna_seq_Normalized_Counts_GLbulkRNAseq.csv \
  "https://osdr.nasa.gov/geode-py/ws/studies/OSD-379/download?file=GLDS-379_rna_seq_Normalized_Counts_GLbulkRNAseq.csv&version=1"

curl -L --fail --show-error \
  --output assets/osdr/GLDS-379_rna_seq_differential_expression_GLbulkRNAseq.csv \
  "https://osdr.nasa.gov/geode-py/ws/studies/OSD-379/download?file=GLDS-379_rna_seq_differential_expression_GLbulkRNAseq.csv&version=1"

curl -L --fail --show-error \
  --output assets/osdr/OSD-379_metadata_OSD-379-ISA.zip \
  "https://osdr.nasa.gov/geode-py/ws/studies/OSD-379/download?file=OSD-379_metadata_OSD-379-ISA.zip&version=6"

unzip -p assets/osdr/OSD-379_metadata_OSD-379-ISA.zip \
  'a_OSD-379_transcription-profiling_rna-sequencing-(rna-seq)_Illumina NovaSeq.txt' \
  > assets/osdr/OSD-379_assay_metadata.tsv
```

Pretrain the released 16-dimensional GLARE SAE on TMS FACS liver cells:

```bash
conda run -n nasa python src/nasa_mouse_glare/reproduce_glare_pretrain.py \
  --input data/glare_inputs/tms_facs_liver_pretrain.mtx \
  --output-dir outputs/glare_paper_tms_liver_osd379/pretraining \
  --epochs 30 \
  --batch-size 16 \
  --seed 1996 \
  --num-workers 0
```

Prepare matched OSD-379 inputs, then run verification before representation
learning as in GLARE:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.paper_finetune \
  --target-manifest data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json \
  --accession OSD-379 \
  --normalized-counts assets/osdr/GLDS-379_rna_seq_Normalized_Counts_GLbulkRNAseq.csv \
  --pretrained-weights outputs/glare_paper_tms_liver_osd379/pretraining/sc_shulse_pretrained_reproduced.pth \
  --output-dir outputs/glare_paper_tms_liver_osd379 \
  --prepare-only

PYTHONPATH=src MPLCONFIGDIR=/tmp/nasa-matplotlib \
  python -m nasa_mouse_glare.paper_analysis verification \
  --run-dir outputs/glare_paper_tms_liver_osd379
```

Fine-tune separate model copies on FLT and GC and cluster both latent spaces:

```bash
PYTHONPATH=src OMP_NUM_THREADS=1 LOKY_MAX_CPU_COUNT=1 \
  python -m nasa_mouse_glare.paper_finetune \
  --target-manifest data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json \
  --accession OSD-379 \
  --normalized-counts assets/osdr/GLDS-379_rna_seq_Normalized_Counts_GLbulkRNAseq.csv \
  --pretrained-weights outputs/glare_paper_tms_liver_osd379/pretraining/sc_shulse_pretrained_reproduced.pth \
  --output-dir outputs/glare_paper_tms_liver_osd379 \
  --epochs 30 \
  --seed 1996

PYTHONPATH=src OMP_NUM_THREADS=1 LOKY_MAX_CPU_COUNT=1 \
  MPLCONFIGDIR=/tmp/nasa-matplotlib \
  python -m nasa_mouse_glare.paper_clustering \
  --run-dir outputs/glare_paper_tms_liver_osd379
```

Use NASA's four age- and collection-matched DESeq2 contrasts for DEG
proportions and export Metascape multi-list files:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.paper_analysis post \
  --run-dir outputs/glare_paper_tms_liver_osd379 \
  --official-de assets/osdr/GLDS-379_rna_seq_differential_expression_GLbulkRNAseq.csv
```

Metascape does not provide a public API. Upload
`biological_analysis/metascape_gene_lists/metascape_multiple_gene_lists.csv`
from the run directory after enabling **Multiple Gene Lists** and confirm
that the first row is detected as the column header. Choose **Custom
Analysis**, select `Mus musculus` for input and analysis, and paste
`biological_analysis/metascape_gene_lists/metascape_background.txt` into
the custom enrichment-background dialog. The CSV contains one foreground
column per eligible cluster; the background file contains all tested genes.
Both use mouse Ensembl gene IDs to avoid ambiguous cross-species symbol
mapping during upload.

See `outputs/glare_paper_tms_liver_osd379/RUN_SUMMARY.md` for results and
method deviations that remain specific to the mouse adaptation.

Audit liver tissue composition when muscle-related clusters or enrichment
appear:

```bash
conda run -n nasa env PYTHONPATH=src MPLCONFIGDIR=/tmp/nasa-matplotlib \
  python -m nasa_mouse_glare.osd379_tissue_qc \
  --run-dir outputs/glare_paper_tms_liver_osd379 \
  --normalized-counts assets/osdr/GLDS-379_rna_seq_Normalized_Counts_GLbulkRNAseq.csv \
  --official-de assets/osdr/GLDS-379_rna_seq_differential_expression_GLbulkRNAseq.csv \
  --tms-h5ad assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad \
  --isa-assay assets/osdr/OSD-379_assay_metadata.tsv
```

This scores all OSD-379 profiles with a fast-skeletal-muscle marker panel,
compares the signal with TMS FACS liver cells, and tests whether the affected
old ISS-terminal cluster remains robust after severe composition outliers are
excluded. It is a sensitivity analysis on normalized counts, not a replacement
DESeq2 analysis from raw integer counts.

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

Download the CARA/OSD-120 normalized counts used for fine-tuning:

```bash
curl -L --fail --show-error \
  --output assets/glare_original/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv \
  "https://osdr.nasa.gov/geode-py/ws/studies/OSD-120/download?file=GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv&version=1"
```

Run the released GLARE fine-tuning config for FLT and GC:

```bash
python src/nasa_mouse_glare/reproduce_glare_finetune.py \
  --data assets/glare_original/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv \
  --pretrained-weights outputs/glare_original_pretrain_config5/sc_shulse_pretrained_reproduced.pth \
  --output-dir outputs/glare_original_finetune_osd120 \
  --epochs 30
```

This reproduces GLARE's CARA restructuring, PCA/k-means outlier step, fixed
adapter transform to the pretraining dimension, and SAE fine-tuning with Adam
`lr=1e-3`, no weight decay, sparsity penalty `1e-5`, and batch size `16`.

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
