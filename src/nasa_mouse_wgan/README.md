# Conditional WGAN-GP workflows

This package implements conditional WGAN-GP generation for mouse bulk RNA-seq.
The final paper uses the matched study-conditioned comparator in
`matched_study.py`. Earlier standalone tissue and pan-tissue runners remain for
development-history reproduction.

## Selected comparator

`matched_study.py` trains and evaluates the WGAN on the same OSDR split and
conditioning structure used for the selected DDIM comparison. It uses the
shared adapter, categorical encoder, metrics, and training partitions from
`nasa_mouse_generative`.

Final commands and the selected run are in:

- [`outputs/COMMANDS.md`](../../outputs/COMMANDS.md)
- [`outputs/README.md`](../../outputs/README.md)
- [`configs/generative/wgan/`](../../configs/generative/wgan/)

## Module map

| Area | Modules |
|---|---|
| Final matched comparator | `matched_study.py` |
| Model implementation | `model.py`, `training.py`, `data.py`, `paths.py` |
| Standalone training and generation | `train_wgan.py`, `run_conditional_generation.py`, `run_pipeline.py`, `generate_synthetic.py` |
| Supporting analysis | `analyze_features.py`, `summarize_results.py` |

## Status boundaries

The matched comparator contributes generator-validation metrics to the final
generative manuscript. The downstream synthetic biological analysis uses the
selected DDIM, not WGAN-generated profiles. Standalone WGAN outputs are not
part of the retained final model set.

The implementation imports shared preprocessing and metrics from
`nasa_mouse_generative`, common IO helpers from `nasa_mouse_glare`, and selected
evaluation helpers from `nasa_mouse_diffusion`.
