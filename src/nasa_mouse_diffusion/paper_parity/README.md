# Paper-parity DDIM and synthetic analysis

This directory contains the final diffusion implementation. It reproduces the
Lacan et al. residual-MLP ModelDDIM architecture for mouse ARCHS4, adapts the
frozen model to API-derived OSDR profiles, evaluates the selected generator,
and runs the paper's tissue-specific synthetic analyses.

## Final execution path

1. `data.py`, `landmarks.py`, and `config.py` prepare the OSDR-disjoint ARCHS4
   reference in the mapped 974-gene space.
2. `train.py` trains the pinned upstream ModelDDIM exposed by `upstream.py`.
3. `conditional_data.py` prepares OSDR and `factorized_train.py` trains the
   study-, tissue-, condition-, and material-aware residual adapter.
4. Factorized evaluation and calibration modules select and test the frozen
   adapter without fitting on locked-test profiles.
5. Matched classifiers and feature analyses compare real-only and
   real-plus-synthetic training on held-out real profiles.
6. Annotation and document builders freeze the publication record.

Run `python -m nasa_mouse_diffusion.paper_parity --help` for the grouped CLI.
Use the exact commands in
[`outputs/COMMANDS.md`](../../../outputs/COMMANDS.md) rather than relying on CLI
defaults for a paper reproduction.

## Module map

| Area | Modules |
|---|---|
| CLI | `__main__.py` |
| ARCHS4 ModelDDIM | `config.py`, `data.py`, `landmarks.py`, `upstream.py`, `train.py`, `evaluate.py` |
| Direct conditional OSDR branch | `conditional_config.py`, `conditional_data.py`, `conditional_train.py`, `conditional_evaluate.py`, `real_effect_ceiling.py` |
| Selected factorized adapter | `factorized_adapter.py`, `factorized_config.py`, `factorized_train.py`, `factorized_evaluate.py`, `factorized_calibrate.py`, `factorized_mean_calibrate.py`, `factorized_distribution_calibrate.py`, `factorized_final_evaluate.py`, `factorized_subset.py`, `factorized_trajectory.py` |
| Guidance and transfer experiments | `contrastive_guidance.py`, `generated_feature_guidance.py`, `generated_feature_transfer.py`, `confirmatory_augmentation.py`, `whole_study_transfer.py` |
| Final downstream analyses | `matched_all_gene_classifiers.py`, `classifier_importance.py`, `grouped_pathway_importance.py`, `within_study_feature_stability.py`, `harmonization_summary.py`, `compare_material_conditioning_ablation.py` |
| Literature and deliverables | `annotate_promoted_gene_literature.py`, `annotate_importance_literature.py`, `build_synthetic_guided_paper.py`, `build_slstp_presentation.py` |

## Evidence boundaries

The selected generator is identified by the model-comparison metrics in the
generative manuscript. Generated data can alter classifier training and feature
ranking, but all reported FLT/GC biological associations are tested in real
OSDR profiles. The annotation tables and sources are indexed in
[`docs/annotation_provenance.md`](../../../docs/annotation_provenance.md).

## Outputs

Final runs and analyses are under `outputs/generative/benchmark/`. The paper
uses frozen source tables under `paper/synthetic_guided_spaceflight/source_data/`.
Large checkpoints and generated matrices remain local and are listed in
[`ARTIFACTS.md`](../../../ARTIFACTS.md).
