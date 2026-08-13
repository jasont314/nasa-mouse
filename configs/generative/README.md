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

## Configuration Families

| Pattern | Role |
|---|---|
| `benchmark/default.yaml` | Shared paths, cohort rules, and benchmark defaults |
| `benchmark/preprocessing_profiles.yaml` | Named expression, feature, scaling, and harmonization choices |
| `benchmark/model_profiles.yaml` | Paper-native model defaults and allowed overrides |
| `benchmark/experiment_matrix.yaml` | Staged development-screen combinations |
| `diffusion/archs4_mouse_paper_parity_osdr_disjoint.yaml` | Selected OSDR-disjoint ARCHS4 backbone |
| `diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint.yaml` | Selected OSDR factorized adapter |
| `diffusion/osdr_factorized_*` | Adapter capacity, conditioning, calibration, and split screens |
| `diffusion/osdr_liver_*` and `diffusion/liver_harmonization_benchmark.yaml` | Liver preprocessing and harmonization comparisons |
| `diffusion/matched_all_gene_*`, `classifier_importance_*`, `grouped_pathway_*`, and `within_study_*` | Downstream synthetic analyses |
| `diffusion/generated_feature_*` and `osdr_whole_study_*` | Confirmatory and whole-study transfer experiments |
| `wgan/wgan_matched_study_conditioned.yaml` | Selected WGAN comparator |
| `wgan/wgan_matched_study_refine_*` | WGAN refinement screens |

Only configurations named as selected in `outputs/README.md` define the final
model set. The other files preserve development screens, ablations, and
sensitivity analyses; they are not additional final models.

Smoke mode remains available in the runners for quick execution checks, but
smoke-run outputs are temporary and are not retained in the repository.
