# ASGSR expiMap HVG paper package

This directory contains the abstract, manuscript, source tables, and reproducible figures from a six-model tissue screen. Thymus, skin, liver, and spleen form the main positive-result set, and kidney is a secondary exploratory result. Soleus did not meet the advancement criteria and is retained only in supplementary screening and sensitivity records.

The primary descriptive quantity is the equally weighted project- or accession-specific flight-minus-ground pathway-score difference after decoder sign orientation. Claims are qualified by conventional enrichment, held-out-project prediction, complete seed retraining, broad composition-proxy sensitivity, member-gene review, tissue fit, and prior literature. FDR and the five-check directional status are reported separately and are not interchangeable discovery thresholds.

## Package contents

- `asgsr_2026_abstract.md`: ASGSR-formatted abstract, under the 300-word limit.
- `manuscript.md`: full manuscript with references and embedded figure captions.
- `poster/`: editable 48 x 27 inch PowerPoint poster based on the approved NASA template, print PDF, high-resolution preview, reusable architecture graphic, and 700-dpi embedded figure assets.
- `supplementary_methods.md`: exact runs, effect definitions, safeguards, and rebuild commands.
- `tissue_selection_audit.md`: reconstruction of the eight-tissue screen and corrected kidney/spleen reassessment.
- `figures/`: publication figures in PNG and vector PDF formats.
- `source_data/table_2_retained_pathway_evidence.tsv`: the 16 retained main and secondary pathway records.
- `source_data/table_s25_revised_model_scope.tsv`: complete six-model screening scope and analysis role.
- `source_data/table_s26_kidney_spleen_training_manifest.tsv`: corrected-model training provenance.
- `source_data/table_s27_kidney_spleen_pathway_evidence.tsv`: complete corrected kidney/spleen evidence matrix.
- `source_data/table_s28_kidney_spleen_manual_review.tsv`: manual review of every corrected-model top-decile pathway.
- `source_data/table_s29_kidney_spleen_member_gene_support.tsv`: member-gene and decoder-weight audit.
- `source_data/table_s30_kidney_spleen_literature_sources.tsv`: source links used in the reassessment.
- `source_data/table_s31_latent_mapping_coordinates.tsv.gz` and `table_s32_latent_mapping_qc.tsv`: reference-query coordinates and mapping-coverage diagnostics.
- `source_data/table_s33_representative_program_sample_scores.tsv.gz`: project-centered sample scores behind Figure S9.
- `source_data/table_s34_retained_pathway_member_gene_support.tsv` and `table_s35_retained_pathway_member_gene_effects.tsv.gz`: pathway-level and gene-level directional support behind Figure 4.
- `source_data/figure_build_manifest.tsv`: dimensions, file sizes, vector-copy status, and border checks for every figure.
- `source_data/table_s3_all_pathway_effects.tsv` and `table_s9_systematic_pathway_screen.tsv`: complete original thymus, skin, liver, and soleus results.
- `source_data/table_s24_pathway_robustness_evidence.tsv`: original four-model five-check matrix.
- `visual_audit.md`: rendered-page and standalone-figure quality-control record.

## Rebuild

Run the original four-model build and robustness analyses first, then the corrected kidney/spleen analysis, paper integration, final publication-figure build, and document renderer:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.build_asgsr_paper
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.reviewer_robustness_analysis
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.run_asgsr_seed_sensitivity
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.integrate_reviewer_robustness
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.run_kidney_spleen_seed_sensitivity
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.analyze_kidney_spleen_reassessment
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.curate_kidney_spleen_reassessment
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.integrate_reassessed_tissues_paper
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.build_publication_figures
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.render_asgsr_documents
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.build_asgsr_poster
```

`build_publication_figures` authors the final main figures and dense supplementary summaries at a 7.2-inch publication width, writes 300-dpi PNG and vector PDF copies, and checks image borders and dimensions. It also produces the latent-mapping, sample-score, and retained-member-gene source tables. Figure 6 is a deterministic Discussion schematic that explicitly separates observed pathway-score directions from tissue-state hypotheses. The broader non-evidentiary process summary remains separate in `presentation/expimap/asgsr_process_summary.*`.

`build_asgsr_poster` creates an editable, single-slide 48 x 27 inch landscape PowerPoint using the approved NASA poster template's 16:9 structure and branding. The project objective and cross-study framing are adapted from the midpoint presentation, while all biological claims use the final robustness-filtered analysis. Native PowerPoint objects are used for poster text, architecture, and tables; the retained-pathway and tissue-state figures are rendered from vector PDFs at 700 dpi, giving at least 300 effective pixels per inch at final placement. When LibreOffice is installed, the command also exports the print PDF, 4,800 x 2,700 pixel preview, and a 300-dpi standalone architecture graphic.

## Model scope

| Tissue | Role | ARCHS4 reference | Mapped OSDR | Primary samples | Primary projects | Genes | Reactome programs |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| thymus | main | 1,362 | 117 | 117 | 5 | 1,994 | 387 |
| skin | main | 2,593 | 151 | 151 | 4 | 1,997 | 319 |
| liver | main | 5,000 | 197 | 197 | 9 | 1,995 | 364 |
| spleen | main | 6,289 | 109 | 100 | 5 | 1,994 | 360 |
| kidney | secondary exploratory | 2,464 | 135 | 135 | 6 | 1,996 | 336 |
| soleus | not advanced; supplementary record | 1,412 | 53 | 53 | 3 | 1,975 | 357 |

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
- No reviewed soleus pathway passes all five checks; the model was not advanced into the paper-facing tissue set, but its complete screening record is preserved.
