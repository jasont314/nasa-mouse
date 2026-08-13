# expiMap workflow

This package implements the pathway-constrained analysis used by the expiMap
manuscript and poster. It prepares tissue-matched ARCHS4 references, maps NASA
OSDR query profiles, estimates flight-minus-ground pathway shifts, checks those
shifts across studies and methods, and builds the publication artifacts.

## Final workflow

The paper-facing path is:

1. Prepare tissue-specific ARCHS4 and OSDR inputs.
2. Select HVGs and build the mouse Reactome mask.
3. Train an ARCHS4 reference model and map the OSDR query.
4. Calculate project-aware pathway effects and member-gene support.
5. Run seed, conventional-enrichment, held-out-project, and composition checks.
6. Freeze source tables, figures, manuscript files, and the poster.

The selected model directories and complete commands are in
[`supplementary_methods.md`](../../paper/asgsr_expimap_hvg/supplementary_methods.md)
and [`outputs/COMMANDS.md`](../../outputs/COMMANDS.md).

## Module map

| Area | Modules |
|---|---|
| Reference and query preparation | `inspect_archs4_mouse.py`, `prepare_expimap_archs4_reference.py`, `prepare_expimap_osdr_tissue.py`, `prepare_expimap_tutorial_hvg.py`, `split_expimap_muscle_groups.py` |
| Training, mapping, and scoring | `train_expimap_archs4_reference.py`, `train_expimap_direct.py`, `map_expimap_osdr_query.py`, `export_expimap_scores.py`, `analyze_expimap_pathways.py`, `validate_expimap_accession_effects.py` |
| Variant and exploratory analyses | `run_expimap_tissue_variant_matrix.py`, `compare_expimap_transformations.py`, `cluster_expimap_condition_scores.py`, `run_expimap_latent_enrich_condition.py`, `summarize_expimap_latent_enrich_bf.py`, `summarize_expimap_de_novo.py`, `summarize_expimap_reference_seeds.py`, `summarize_expimap_results.py` |
| Muscle and pathway follow-up | `analyze_muscle_targeted_modules.py`, `plot_expimap_pathway_followup.py`, `plot_expimap_accession_heatmaps.py`, `plot_hvg_interpretation_heatmaps.py`, `plot_hvg_literature_review_heatmaps.py` |
| Robustness and review | `audit_aligned_complementary_programs.py`, `review_expanded_pathway_screen.py`, `reviewer_robustness_analysis.py`, `run_asgsr_seed_sensitivity.py`, `integrate_reviewer_robustness.py` |
| Kidney and spleen reassessment | `run_kidney_spleen_seed_sensitivity.py`, `analyze_kidney_spleen_reassessment.py`, `curate_kidney_spleen_reassessment.py`, `integrate_reassessed_tissues_paper.py` |
| Publication builders | `build_asgsr_paper.py`, `build_publication_figures.py`, `build_asgsr_poster.py`, `render_asgsr_documents.py`, `verify_asgsr_sources.py` |

## Final and supporting code

The ARCHS4 reference-query HVG models and robustness builders support the final
paper. Direct OSDR training, all-gene, de novo, latent-enrichment, clustering,
and transformation comparisons are retained as development or sensitivity
analyses. Their presence does not make them part of the final evidence set.
`outputs/README.md` identifies the runs that were selected.

## Inputs and outputs

- OSDR samples come from the NASA API utilities in `nasa_mouse_glare`.
- The official mouse Reactome GMT is
  `data/pathways/reactome_current_mouse_ensembl.gmt`.
- Tissue-matched ARCHS4 reference data use the local mouse H5 listed in
  `ARTIFACTS.md`.
- Selected outputs are under `outputs/expimap/`.
- Publication artifacts are under `paper/asgsr_expimap_hvg/` and
  `presentation/poster/`.

The literature annotation record is documented in
[`docs/annotation_provenance.md`](../../docs/annotation_provenance.md).
