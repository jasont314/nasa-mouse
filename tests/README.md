# Test suite

The tests cover project-owned data contracts, model wrappers, downstream
analyses, paper builders, and final handoff integrity. They use small synthetic
fixtures and tracked source tables; they do not retrain the full models.

Run the suite from the repository root:

```bash
python -m pytest
```

## Test groups

| Area | Files |
|---|---|
| Shared generative configuration and runtime | `test_generative_pipeline.py`, `test_generative_orchestration.py`, `test_generative_runtime.py` |
| OSDR and condition figures | `test_glare_osdr_api.py`, `test_condition_figures.py` |
| DDIM and WGAN parity | `test_rna_diffusion_parity.py`, `test_wgan_matched_study.py` |
| Classifier and feature analyses | `test_matched_all_gene_classifiers.py`, `test_classifier_importance.py`, `test_grouped_pathway_importance.py` |
| Guidance and transfer | `test_generated_feature_guidance.py`, `test_whole_study_transfer.py`, `test_within_study_feature_stability.py` |
| Papers and handoff | `test_paper_workflow.py`, `test_handoff_integrity.py` |

The parity tests require the pinned upstream source trees under
`assets/model_sources/`. Restore missing optional checkouts with:

```bash
python -m nasa_mouse_generative prepare-upstreams
```

`test_handoff_integrity.py` checks documentation links, required deliverables,
annotation source resolution, repository paths, licensing, and paper manifests.
See [`REPRODUCIBILITY.md`](../REPRODUCIBILITY.md) for the full fresh-clone and
external-artifact verification sequence.
