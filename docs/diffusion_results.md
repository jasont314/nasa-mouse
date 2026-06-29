# Diffusion Results

This file will summarize completed NASA mouse diffusion runs after smoke and
production training finish.

## Smoke Validation

Smoke validation was run with two OSDR tissues, 256 genes, 64 landmarks, 20
diffusion timesteps, 5 DDIM sample steps, a 32-unit residual MLP, and one query
epoch. It used CUDA on the NVIDIA A100-SXM4-40GB.

Smoke outputs:

- `outputs/diffusion_smoke/osdr_only/`
- `outputs/diffusion_smoke/osdr_only/analysis/`
- `outputs/diffusion_smoke/synthetic_counterfactual/`

Smoke counts:

- query samples: 422
- genes: 256
- landmark genes: 64
- denoiser features: 32

The smoke run wrote a model checkpoint, denoiser feature scores, generated
quality metrics, LR reconstruction metrics, reverse-validation metrics,
PCA/UMAP plots, feature-shift heatmap, and matched synthetic
ground-control/flight matrices. Because the smoke subset only used the first
256 genes, the L1000 ortholog map had insufficient coverage and the run used
the intended variance-HVG fallback.

Planned result sections:

- OSDR-only conditional diffusion
- ARCHS4 pretrain plus OSDR fine-tune
- ARCHS4-only reference generator/control
- generated-expression quality
- FLT vs GC denoiser-feature analysis
- synthetic counterfactual examples
- comparison to WGAN, OntoVAE, and expiMap
- skeletal-muscle and split-muscle interpretation
