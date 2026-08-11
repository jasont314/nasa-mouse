# Generative Model Paper And Code Audit

This audit fixes the scope of the two generative-model adapters before implementation.
Repository commits were inspected on 2026-07-14.

## Vinas Conditional WGAN-GP

- Paper: [Adversarial generation of gene expression data](https://doi.org/10.1093/bioinformatics/btab035)
- Code: [rvinas/adversarial-gene-expression](https://github.com/rvinas/adversarial-gene-expression)
- Inspected commit: `94fa44dd1bd52d924efd3af0fcd8eeb18bd141a8`
- Native task: conditional bulk-expression generation.
- Generator and critic concatenate quantitative covariates and learned embeddings
  for categorical covariates.
- Human RNA-seq code uses `log(1+x)`, a train-fitted gene-wise standard score,
  latent dimension 64, two 256-unit hidden layers, batch size 32, five critic
  updates per generator update, RMSprop at `5e-4`, and at most 2,000 epochs.
- The source's zero-based loop evaluates its gamma score after epochs 1, 6, 11, and
  so on, with default patience ten checks. This is up to 50 no-improvement epochs,
  whereas the paper text states 30 epochs. Both variants are named explicitly in
  the benchmark.
- The paper conditions GTEx/TCGA generation on tissue and dataset. Holding latent
  noise fixed while changing a category is its counterfactual construction.
- The official split is sample-random. This benchmark replaces it with accession
  grouping to avoid study leakage.
- ARCHS4 pretraining followed by OSDR fine-tuning is not native and must be
  implemented and labeled as a NASA extension.

## Lacan Landmark-Space Diffusion

- Paper: [In silico generation of gene expression profiles using diffusion models](https://doi.org/10.1186/s12859-026-06470-8)
- Code: `https://forge.ibisc.univ-evry.fr/alacan/rna-diffusion.git`
- Inspected commit: `cde890154698fcea96c924804aaff04af3351b48`
- Native task: tissue-conditioned bulk-expression generation.
- A residual MLP DDPM/DDIM predicts Gaussian noise for approximately 1,000 L1000
  landmark genes. Linear regression or an MLP reconstructs target genes.
- Paper data are TPM-normalized GTEx or processed TCGA expression. The selected
  models use MaxAbs scaling, 1,000 diffusion steps, a quadratic beta schedule,
  large batches, AMP, EMA, and DDIM sampling.
- The paper tunes models using unsupervised precision/recall F1 and evaluates
  correlation, Frechet distance, precision/recall, adversarial accuracy, PCA/UMAP,
  and train-on-synthetic/test-on-real reverse validation.
- ARCHS4 pretraining followed by OSDR fine-tuning is a NASA extension.

## Consequence For The Benchmark

The synthetic-expression comparison contains two native generators: Vinas WGAN-GP
and Lacan diffusion. Pinned source files are verified by SHA-256 in addition to
checking Git commits.
