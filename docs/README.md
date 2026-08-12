# Documentation index

The final papers are the main scientific record. These notes document the code,
intermediate decisions, sensitivity analyses, and audits behind them.

## Start with these

- [NASA OSDR API](osdr_api.md): data discovery, endpoints, and local metadata.
- [Method sources](method_sources.md): upstream repositories and pinned commits.
- [expiMap results](expimap_results.md): current pathway-analysis summary.
- [Generative benchmark results](generative_benchmark_results.md): model and
  preprocessing comparison.
- [Generative pipeline](generative_models_pipeline.md): configurable benchmark
  contract and execution model.
- [Distributed response hypotheses](distributed_response_hypotheses.md):
  biological interpretation developed for the internship report.

## GLARE and batch correction

- [GLARE validation stack](glare_validation_stack.md)
- [15-term GLARE candidate screen](glare_validation_stack_terms15.md)
- [Signed module direction](glare_signed_module_direction.md)
- [Literature cross-check](glare_literature_crosscheck.md)
- [Metascape wrapper](metascape.md)

GLARE is an exploratory first stage in the final report. Its main conclusion is
that accession structure remained prominent after correction.

## expiMap

- [Accession-aware validation](expimap_accession_validation.md)
- [Condition-specific clustering](expimap_condition_clustering.md)
- [Liver de novo programs](expimap_de_novo_liver.md)
- [Literature comparison](expimap_literature_comparison.md)
- [Reference seed stability](expimap_reference_seed_stability.md)
- [Skeletal-muscle prior-work check](expimap_skeletal_muscle_prior_work.md)
- [Tutorial-style liver run](expimap_tutorial_style_liver.md)

The current publication-facing evidence and exact model scope are in
`paper/asgsr_expimap_hvg/`.

## Generative methods and validation

- [Data audit](generative_data_audit.md)
- [Paper and code audit](generative_model_code_audit.md)
- [Benchmark decisions](generative_benchmark_decisions.md)
- [Benchmark execution log](generative_benchmark_execution.md)
- [Diffusion pipeline](diffusion_pipeline.md)
- [Diffusion results](diffusion_results.md)
- [RNA diffusion paper-parity baseline](rna_diffusion_paper_parity.md)
- [Selected OSDR diffusion model](osdr_conditional_diffusion_finalist.md)
- [Leakage-control confirmation](diffusion_leak_free_confirmation.md)
- [WGAN pipeline](wgan_pipeline.md)
- [WGAN conditional generation](wgan_conditional_generation.md)
- [WGAN results](wgan_results.md)
- [Study-conditioned WGAN](osdr_conditional_wgan_study.md)

## Synthetic-informed analysis

- [Matched all-gene classifiers](matched_all_gene_classifier_analysis.md)
- [Classifier feature importance](classifier_importance.md)
- [Generated feature-guidance workflow](generated_feature_guidance_workflow.md)
- [Promoted-gene literature annotation](promoted_gene_literature_annotation.md)
- [Skeletal-muscle group analysis](synthetic_skeletal_muscle_group_analysis.md)
- [Whole-study transfer](whole_study_transfer.md)
- [Thymus platelet-calcium follow-up](thymus_platelet_calcium_validation.md)

The final tables used by the manuscript are frozen under
`paper/synthetic_guided_spaceflight/source_data/`.

## Archived work logs

Files under [`archive/`](archive/) preserve transfer notes and bounded side
analyses from development. They may contain old machine paths, run names, and
instructions that no longer describe the final repository. Use the current
README files and `outputs/COMMANDS.md` for reproducibility.
