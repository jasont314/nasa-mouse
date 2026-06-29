# WGAN Bulk RNA-seq Pipeline

This document tracks the conditional WGAN-GP workflow for NASA mouse OSDR bulk
RNA-seq. The code lives in `src/nasa_mouse_wgan/`, separate from
`src/nasa_mouse_glare/`.

## Paper Basis

Reference:

- Viñas et al., "Adversarial generation of gene expression data",
  Bioinformatics, 2022.
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC8756177/
- DOI: https://doi.org/10.1093/bioinformatics/btab035
- Code link reported by the paper:
  https://github.com/rvinas/adversarial-gene-expression

The paper describes a conditional Wasserstein GAN with gradient penalty
(WGAN-GP), not a plain GAN. The generator and critic are multilayer perceptrons
conditioned on numerical and categorical covariates. Categorical covariates are
represented with learned embeddings. The critic loss is the WGAN critic
objective plus a gradient penalty along interpolations between real and
generated samples.

The paper's human RNA-seq case study used integrated GTEx + TCGA expression
data, condition/tissue/sex categorical covariates, age as a numerical covariate,
standardized numerical covariates, two hidden layers with 256 units, ReLU hidden
activations, and linear generator/critic outputs.

## NASA Mouse Adaptation

This project uses API-derived OSDR/ARCHS4 AnnData inputs already generated in
the repository. It must not use the older raw integrated OSDR H5.

Training modes:

- `direct_osdr`: train a conditional WGAN-GP directly on OSDR FLT/GC samples.
- `archs4_pretrain_osdr_finetune`: pretrain on tissue-matched ARCHS4 mouse
  reference samples, score/query OSDR samples through the pretrained critic,
  then fine-tune on OSDR.

Inputs:

- OSDR query: `outputs/expimap_direct_osdr_<tissue>/input/*.h5ad`
- ARCHS4 reference:
  `outputs/expimap_archs4_reference_osdr_query_<tissue>/reference_input_*/`
- Skeletal-muscle split OSDR groups:
  `outputs/expimap_muscle_targeted_combined_min8/group_inputs_exploratory_2acc/`

Expression transform:

- raw counts from AnnData
- `log1p(CPM)`
- gene-wise z-score using reference statistics for pretrained runs, query
  statistics for direct runs

Categorical covariates:

- OSDR condition: `condition_inferred`
- OSDR accession: `id.accession`
- ARCHS4 reference source: `series_id`/`archs4_condition`

## Planned Outputs

Per run:

- `critic_feature_scores.tsv`: OSDR sample-level critic feature embeddings and
  critic score after final training
- `pretrained_query_critic_feature_scores.tsv`: OSDR sample-level critic
  features before OSDR fine-tuning, for ARCHS4-pretrained runs
- `generated_quality.tsv` / `generated_quality.json`: real-vs-generated
  distribution checks
- `training_summary.json`: device, sample counts, genes, covariates, losses
- `model.pt` and optional `reference_pretrained_model.pt`

Per score set analysis:

- FLT vs GC feature comparison
- accession-aware random-effects validation
- leave-one-accession-out validation
- PCA/UMAP plots of critic features
- heatmap of top FLT/GC shifted critic features

## Validation Checkpoints

Smoke validation was run for liver with 256 genes, one epoch, and a small
32-unit critic/generator to verify the workflow rather than call biology. Both
paths used CUDA on the NVIDIA A100-SXM4-40GB:

- direct OSDR training wrote scores, FLT-vs-GC tables, random-effects and LOO
  summaries, PCA/UMAP plots, and a top-feature heatmap under
  `outputs/wgan_smoke_liver/direct_osdr/`;
- ARCHS4-pretrained training wrote pretrained-query scores, post-fine-tune
  scores, the same analysis outputs, and both pretrained/final model files under
  `outputs/wgan_smoke_liver/archs4_pretrain_osdr_finetune/`.

Smoke outputs are excluded from the aggregate production summarizer by default.

## Limitations

The WGAN paper is a generative expression simulator. It does not natively
produce Reactome pathway scores, de novo gene programs, or scArches-style query
mapping. For NASA mouse comparison, the most defensible signal is therefore:

- learned critic feature separation between FLT and GC
- stability across OSDR accessions
- whether ARCHS4 pretraining improves stable feature signal
- whether top shifted features correlate with generated-expression/gene-level
  patterns that can be followed up with DGEA or pathway enrichment

WGAN outputs should be treated as complementary representation/generative
evidence, not a direct replacement for expiMap pathway latents or OntoVAE
Reactome-constrained programs.
