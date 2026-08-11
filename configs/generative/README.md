# Generative Model Configurations

All generative-model configurations live under this directory. They are split
by responsibility so shared benchmark settings are not mixed with model runs.

| Directory | Purpose |
|---|---|
| `benchmark/` | Shared data, preprocessing, model-profile, and experiment-matrix settings |
| `diffusion/` | DDIM pretraining, OSDR adaptation, validation, and feature-analysis runs |
| `wgan/` | Conditional WGAN training and refinement runs |

The main shared entry point is `benchmark/default.yaml`. The final DDIM and
WGAN configurations used by the current analyses are documented in
[`outputs/README.md`](../../outputs/README.md), with runnable commands in
[`outputs/COMMANDS.md`](../../outputs/COMMANDS.md).

Smoke mode remains available in the runners for quick execution checks, but
smoke-run outputs are temporary and are not retained in the repository.
