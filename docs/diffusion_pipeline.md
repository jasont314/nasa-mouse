# Diffusion Bulk RNA-seq Pipeline

This document tracks the NASA mouse adaptation of Lacan et al.,
"In silico generation of gene expression profiles using diffusion models"
(BMC Bioinformatics, 2026; DOI: https://doi.org/10.1186/s12859-026-06470-8).

The downloaded accepted-manuscript PDF is stored at:

- `docs/papers/lacan_diffusion_2026_reference.pdf`

Code lives in `src/nasa_mouse_diffusion/`, separate from
`src/nasa_mouse_glare/` and `src/nasa_mouse_wgan/`.

## Paper And Code Review

The paper proposes a bulk RNA-seq synthetic generation pipeline:

- train DDPM/DDIM in reduced L1000 landmark-gene space;
- condition the denoiser on tissue type;
- use a residual MLP denoiser instead of an image U-Net;
- use 1,000 diffusion timesteps, quadratic beta schedule, AMP, large batches,
  and DDIM sampling for faster generation;
- reconstruct remaining target genes from generated landmarks using linear
  regression or MLP;
- evaluate synthetic data with correlation score, precision/recall, Frechet
  distance in classifier embedding space, adversarial accuracy, UMAP/PCA, and
  reverse validation.

The official code was inspected from:

- https://forge.ibisc.univ-evry.fr/alacan/rna-diffusion.git

Important implementation details from the code:

- the usable DDIM implementation lives under `src/generation/ddim/`;
- model class: `ModelDDIM`, a residual MLP over tabular gene vectors;
- loss: direct Gaussian-noise prediction MSE;
- sampler: generalized DDIM with `eta=0`, and DDPM-like sampling with `eta=1`;
- beta schedules: quadratic, linear, constant, JSD, and sigmoid;
- training uses antithetic timestep sampling, EMA, AMP, Adam, and optional
  learning-rate scheduling;
- paper configs use MaxAbs scaling, 974 GTEx or 978 TCGA landmark genes,
  hidden dimensions 8192/4096, dropout 0.1, quadratic beta schedule, and no
  sinusoidal time embedding.

Paper/code caveats:

- several official scripts contain hardcoded absolute paths under the author's
  home directory;
- top-level `main_ddim.py` imports `runners.diffusion`, while the actual runner
  is under `src/generation/ddim/runners/diffusion.py`;
- Frechet scoring depends on pretrained tissue classifiers in unavailable
  hardcoded paths;
- the NASA implementation therefore reimplements the method in local repo style
  instead of vendoring the paper scripts directly.

## NASA Mouse Adaptation

Inputs must remain API-derived:

- OSDR query: `outputs/expimap_direct_osdr_<tissue>/input/*.h5ad`
- ARCHS4 reference:
  `outputs/expimap_archs4_reference_osdr_query_<tissue>/reference_input_*/`

The workflow does not use the older raw integrated OSDR H5 files.

Training tracks:

- `osdr_only`: train conditional diffusion directly on OSDR FLT/GC samples.
- `archs4_pretrain_osdr_finetune`: pretrain on ARCHS4 mouse reference samples,
  score OSDR through the pretrained denoiser, then fine-tune on OSDR.
- `archs4_only`: train an ARCHS4 reference-distribution generator/control.

Conditioning covariates match the WGAN conditional-generation design:

- `wgan_condition`: `flight`, `ground_control`, or `archs4_reference`
- `wgan_tissue`
- `wgan_material_type`
- `wgan_muscle_group`
- `wgan_accession`
- `wgan_sex`
- `wgan_assay`
- `wgan_platform`
- `wgan_data_source`

The names reuse the already harmonized WGAN columns for consistency with the
existing repo data helpers.

## Landmark Genes

The OSDR/ARCHS4 AnnData files currently carry mouse Ensembl IDs
(`ENSMUSG...`) but not gene symbols. The paper's L1000 genes are human
landmarks, so this workflow supports:

- mapping the paper's human L1000 Ensembl genes to mouse ortholog Ensembl IDs
  using Ensembl BioMart:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_diffusion.map_l1000_mouse
```

- fallback to top-variance mouse genes from the aligned OSDR/ARCHS4 matrix if
  the L1000 map is unavailable or insufficient.

The selected landmark list is written per run as `landmark_genes.tsv`, and the
training summary records whether the source was `mouse_l1000_orthologs` or
`variance_hvg_fallback`.

## Outputs

Default output root:

- `outputs/diffusion_conditional_generation/`

Per model:

- `model.pt`
- optional `reference_pretrained_model.pt`
- `diffusion_feature_scores.tsv`
- optional `pretrained_query_diffusion_feature_scores.tsv`
- `generated_quality.json` / `generated_quality.tsv`
- `training_summary.json`
- `observed_conditioning_profiles.tsv`
- `genes.tsv`
- `landmark_genes.tsv`
- `normalization_stats.npz`

Feature analyses write:

- FLT vs GC feature comparison
- accession-aware random-effects meta-analysis
- leave-one-accession-out validation
- PCA/UMAP feature plots
- top feature-shift heatmap

Synthetic generation writes scaled, log1p-CPM, and CPM matrices. The scaled
matrix is in the diffusion model's MaxAbs-scaled space. Inverse-transformed
log1p-CPM values are clipped to the valid CPM range `[0, log1p(1e6)]` before
CPM export, and each condition writes a `*_clip_report.json` file.

Production summary tables are written under:

- `outputs/diffusion_conditional_generation/summary/diffusion_training_summary.tsv`
- `outputs/diffusion_conditional_generation/summary/diffusion_analysis_summary.tsv`
- `outputs/diffusion_conditional_generation/summary/diffusion_synthetic_examples_summary.tsv`
- `outputs/diffusion_conditional_generation/summary/diffusion_subgroup_analysis_summary.tsv`

## Current Deviations From The Paper

- The local first pass uses a smaller residual MLP by default than the paper's
  4096/8192 hidden dimensions so smoke and production runs fit the current
  iterative workflow.
- Frechet distance is computed in PCA space rather than a paper-trained GTEx/TCGA
  classifier embedding because the official pretrained classifier weights are not
  distributed in the cloned repo.
- The implemented reconstruction path is ridge linear regression. MLP
  reconstruction remains a follow-up unless needed by the results.
- The main biological FLT/GC analyses use denoiser feature embeddings and
  generated counterfactual deltas. These features are label-conditioned, so they
  should be interpreted as synthetic-generation diagnostics rather than unlabelled
  pathway modules.
- Frozen ARCHS4 reference projection must not use OSDR-only covariate embeddings
  that were unseen during pretraining. The workflow now keeps query expression
  and tissue but replaces condition/accession/source/platform/material/sex with
  trained ARCHS4 defaults by tissue before scoring `pretrained_query_*` outputs.

## Commands

Dry-run manifest:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_diffusion.run_conditional_generation --dry-run
```

Production pan-tissue run:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_diffusion.run_conditional_generation
```

Generate matched synthetic samples:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_diffusion.generate_synthetic \
  --model-dir outputs/diffusion_conditional_generation/archs4_pretrain_osdr_finetune \
  --output-dir outputs/diffusion_conditional_generation/synthetic_custom/liver \
  --n 128 \
  --counterfactual ground_control flight \
  --set wgan_tissue=liver \
  --set wgan_material_type=Liver \
  --set wgan_muscle_group=not_skeletal_muscle \
  --set wgan_accession=OSD-137 \
  --set wgan_sex=Female \
  --set "wgan_assay=RNA Sequencing (RNA-Seq)" \
  --set wgan_platform=illumina_hiseq_4000 \
  --set wgan_data_source=osdr
```

Generate the standard per-tissue and skeletal-muscle-split examples:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_diffusion.generate_synthetic_examples --overwrite
```

Correct frozen ARCHS4 projection scores from a saved pretrained checkpoint:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_diffusion.rescore_reference_projection
```

Refresh compact result tables:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_diffusion.summarize_results
```
