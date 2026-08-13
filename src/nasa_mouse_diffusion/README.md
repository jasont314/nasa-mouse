# Diffusion workflows

This package contains two generations of bulk RNA-seq diffusion code.

- [`paper_parity/`](paper_parity/) is the final Lacan et al. ModelDDIM branch
  used by the generative manuscript and presentation.
- The modules at this package root form the earlier standalone conditional
  diffusion pipeline. They remain for development-history reproduction and are
  not the source of the final synthetic-analysis claims.

## Final workflow

Read [`paper_parity/README.md`](paper_parity/README.md) for ARCHS4 pretraining,
OSDR adaptation, locked evaluation, downstream feature analysis, literature
annotation, and document builds.

The selected configurations and commands are recorded in:

- [`configs/generative/README.md`](../../configs/generative/README.md)
- [`outputs/README.md`](../../outputs/README.md)
- [`outputs/COMMANDS.md`](../../outputs/COMMANDS.md)

## Standalone module map

| Area | Modules |
|---|---|
| Data and paths | `data.py`, `paths.py`, `map_l1000_mouse.py` |
| Model and sampling | `model.py`, `diffusion.py`, `reconstruction.py` |
| Training and generation | `train_diffusion.py`, `run_conditional_generation.py`, `generate_synthetic.py`, `generate_synthetic_examples.py` |
| Evaluation and analysis | `evaluate.py`, `analyze_features.py`, `analyze_subgroups.py`, `rescore_reference_projection.py`, `refresh_reverse_validation.py`, `summarize_results.py` |

These standalone modules support the earlier OSDR-only, ARCHS4-only, and
ARCHS4-pretrain plus OSDR-fine-tune experiments. Their commands remain under
the "Standalone Generative Runs" section of `outputs/COMMANDS.md`. Their output
tree is not part of the final retained model set.

## Shared dependencies

The package uses `nasa_mouse_glare` for common matrix/import helpers and
`nasa_mouse_generative` for paper-aligned metrics and benchmark contracts. The
shared framework imports this package's model and sampler through its diffusion
adapter.
