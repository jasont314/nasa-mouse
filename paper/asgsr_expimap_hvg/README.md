# ASGSR expiMap HVG paper package

This directory contains the abstract, manuscript, source tables, and reproducible figures for six tissue-specific HVG expiMap models. Thymus, skin, liver, and spleen form the main positive-result set. Kidney is a secondary exploratory result, and soleus is retained as a supplementary negative sensitivity analysis.

The primary descriptive quantity is the equally weighted project- or accession-specific flight-minus-ground pathway-score difference after decoder sign orientation. Claims are qualified by conventional enrichment, held-out-project prediction, complete seed retraining, broad composition-proxy sensitivity, member-gene review, tissue fit, and prior literature. FDR and the five-check directional status are reported separately and are not interchangeable discovery thresholds.

## Package contents

- `asgsr_2026_abstract.md`: ASGSR-formatted abstract, under the 300-word limit.
- `manuscript.md`: full manuscript with references and embedded figure captions.
- `supplementary_methods.md`: exact runs, effect definitions, safeguards, and rebuild commands.
- `tissue_selection_audit.md`: reconstruction of the eight-tissue screen and corrected kidney/spleen reassessment.
- `figures/`: publication figures in PNG and vector PDF formats.
- `source_data/table_2_retained_pathway_evidence.tsv`: the 16 retained main and secondary pathway records.
- `source_data/table_s25_revised_model_scope.tsv`: six-model scope and main, secondary, or supplementary role.
- `source_data/table_s26_kidney_spleen_training_manifest.tsv`: corrected-model training provenance.
- `source_data/table_s27_kidney_spleen_pathway_evidence.tsv`: complete corrected kidney/spleen evidence matrix.
- `source_data/table_s28_kidney_spleen_manual_review.tsv`: manual review of every corrected-model top-decile pathway.
- `source_data/table_s29_kidney_spleen_member_gene_support.tsv`: member-gene and decoder-weight audit.
- `source_data/table_s30_kidney_spleen_literature_sources.tsv`: source links used in the reassessment.
- `source_data/table_s3_all_pathway_effects.tsv` and `table_s9_systematic_pathway_screen.tsv`: complete original thymus, skin, liver, and soleus results.
- `source_data/table_s24_pathway_robustness_evidence.tsv`: original four-model five-check matrix.
- `visual_audit.md`: rendered-page and standalone-figure quality-control record.

## Rebuild

Run the original four-model build and robustness analyses first, then the corrected kidney/spleen analysis, paper integration, and document renderer:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.build_asgsr_paper
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.reviewer_robustness_analysis
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.run_asgsr_seed_sensitivity
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.integrate_reviewer_robustness
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.run_kidney_spleen_seed_sensitivity
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.analyze_kidney_spleen_reassessment
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.curate_kidney_spleen_reassessment
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.integrate_reassessed_tissues_paper
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.render_asgsr_documents
```

## Model scope

| Tissue | Role | ARCHS4 reference | Mapped OSDR | Primary samples | Primary projects | Genes | Reactome programs |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| thymus | main | 1,362 | 117 | 117 | 5 | 1,994 | 387 |
| skin | main | 2,593 | 151 | 151 | 4 | 1,997 | 319 |
| liver | main | 5,000 | 197 | 197 | 9 | 1,995 | 364 |
| spleen | main | 6,289 | 109 | 100 | 5 | 1,994 | 360 |
| kidney | secondary exploratory | 2,464 | 135 | 135 | 6 | 1,996 | 336 |
| soleus | supplementary sensitivity | 1,412 | 53 | 53 | 3 | 1,975 | 357 |

## Interpretation safeguards

- Higher or lower refers to the decoder-oriented expiMap latent pathway score, not uniform expression of every member gene.
- Independently trained latent scales are interpreted within tissue and are not compared quantitatively across tissues.
- OSD-288 is excluded from the primary spleen effect because recorded strain is disjoint by condition.
- OSD-289 thymus and OSD-714 soleus are excluded in restricted sensitivities for the same type of recorded confounding.
- The primary liver query excludes overlapping OSD-164 and OSD-168 cohort representations; the original 12-accession mapping remains sensitivity evidence.
- Paired skin sites are collapsed within mission project for project-level checks, and protocol depooling separates microgravity, artificial 1 g, recovery, terminal collection, strain, and duration.
- `Triangulated` means five supported directions; it is not a significance level. Kidney passes five directional checks but remains exploratory because all three GSEA FDR values exceed 0.05 and composition adjustment attenuates the effects.
- Spleen T-cell receptor, neutrophil degranulation, and C-type lectin receptor programs are lower across five unconfounded projects and three trainings and each has GSEA FDR below 0.05.
- The neutrophil degranulation score is a transcriptomic program, not a direct assay of neutrophil degranulation or suppressive function.
- No reviewed soleus pathway passes all five checks; this is a negative result for the current three-project bulk model, not evidence that soleus is unaffected by spaceflight.
