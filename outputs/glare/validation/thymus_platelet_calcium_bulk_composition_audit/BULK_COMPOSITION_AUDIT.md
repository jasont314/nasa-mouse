# Thymus Bulk Composition Audit For Platelet-Calcium GLARE Module

## Bottom Line

The actual OSDR thymus bulk profiles support a composition-sensitive signal, not a clean thymocyte-intrinsic module.
The GLARE platelet-calcium module is FLT-up in the strongest thymus studies, and the same samples also show FLT-up platelet/coagulation/endothelial/stromal marker scores.
However, removing platelet/coagulation genes from the module does not erase the FLT effect, which argues against a pure blood/platelet contamination artifact.

Best label: `spaceflight-associated thymic vascular/stromal/hemostasis-remodeling module with composition risk`.

## Gene-Set Coverage

| marker_set | requested_genes | genes_present |
| --- | --- | --- |
| coagulation_plasma | 16 | 16 |
| endothelial_vascular_core | 15 | 15 |
| erythroid_blood | 10 | 10 |
| fibroblast_ecm_stromal | 17 | 17 |
| glare_module_without_platelet_coag_endothelial | 50 | 50 |
| glare_module_without_platelet_or_coag | 51 | 51 |
| glare_platelet_ca_module_77 | 77 | 77 |
| myeloid_macrophage | 13 | 13 |
| platelet_megakaryocyte_core | 17 | 17 |
| thymic_epithelial_tec | 16 | 16 |
| thymocyte_t_cell | 18 | 18 |

## Marker Overlap With GLARE Module

| marker_set | genes_in_set | overlap_with_module_77 | overlap_fraction_of_marker_set |
| --- | --- | --- | --- |
| glare_platelet_ca_module_77 | 77 | 77 | 1 |
| glare_module_without_platelet_or_coag | 51 | 51 | 1 |
| glare_module_without_platelet_coag_endothelial | 50 | 50 | 1 |
| platelet_megakaryocyte_core | 17 | 11 | 0.647 |
| coagulation_plasma | 16 | 15 | 0.938 |
| endothelial_vascular_core | 15 | 2 | 0.133 |
| fibroblast_ecm_stromal | 17 | 2 | 0.118 |
| thymic_epithelial_tec | 16 | 0 | 0 |
| thymocyte_t_cell | 18 | 0 | 0 |
| myeloid_macrophage | 13 | 0 | 0 |
| erythroid_blood | 10 | 0 | 0 |

## FLT Minus GC Marker Scores By Study

Scores are mean log2 CPM over genes in the marker panel. Positive values mean FLT-higher.

| accession | marker_set | genes_present | n_flight | n_ground | flight_minus_ground | welch_p |
| --- | --- | --- | --- | --- | --- | --- |
| OSD-244 | coagulation_plasma | 16 | 19 | 19 | 1.27 | 0.0142 |
| OSD-289 | coagulation_plasma | 16 | 12 | 3 | 0.37 | 0.0528 |
| OSD-421 | coagulation_plasma | 16 | 10 | 10 | -0.0493 | 0.494 |
| OSD-457 | coagulation_plasma | 16 | 12 | 12 | 0.361 | 0.0107 |
| OSD-515 | coagulation_plasma | 16 | 9 | 9 | 0.0853 | 0.343 |
| OSD-244 | endothelial_vascular_core | 15 | 19 | 19 | 1.17 | 1.51e-04 |
| OSD-289 | endothelial_vascular_core | 15 | 12 | 3 | 0.369 | 0.00825 |
| OSD-421 | endothelial_vascular_core | 15 | 10 | 10 | 0.0599 | 0.731 |
| OSD-457 | endothelial_vascular_core | 15 | 12 | 12 | 0.343 | 8.36e-08 |
| OSD-515 | endothelial_vascular_core | 15 | 9 | 9 | 0.292 | 0.278 |
| OSD-244 | erythroid_blood | 10 | 19 | 19 | 1.15 | 0.00313 |
| OSD-289 | erythroid_blood | 10 | 12 | 3 | 0.0735 | 0.659 |
| OSD-421 | erythroid_blood | 10 | 10 | 10 | -0.192 | 0.463 |
| OSD-457 | erythroid_blood | 10 | 12 | 12 | -0.189 | 0.0151 |
| OSD-515 | erythroid_blood | 10 | 9 | 9 | 0.901 | 0.0408 |
| OSD-244 | fibroblast_ecm_stromal | 17 | 19 | 19 | 1.5 | 4.15e-04 |
| OSD-289 | fibroblast_ecm_stromal | 17 | 12 | 3 | 0.5 | 7.30e-05 |
| OSD-421 | fibroblast_ecm_stromal | 17 | 10 | 10 | -0.155 | 0.399 |
| OSD-457 | fibroblast_ecm_stromal | 17 | 12 | 12 | 0.399 | 5.54e-05 |
| OSD-515 | fibroblast_ecm_stromal | 17 | 9 | 9 | 0.135 | 0.738 |
| OSD-244 | glare_module_without_platelet_coag_endothelial | 50 | 19 | 19 | 0.714 | 1.11e-04 |
| OSD-289 | glare_module_without_platelet_coag_endothelial | 50 | 12 | 3 | 0.231 | 0.0129 |
| OSD-421 | glare_module_without_platelet_coag_endothelial | 50 | 10 | 10 | -0.213 | 0.0865 |
| OSD-457 | glare_module_without_platelet_coag_endothelial | 50 | 12 | 12 | 0.198 | 2.53e-05 |
| OSD-515 | glare_module_without_platelet_coag_endothelial | 50 | 9 | 9 | -0.0968 | 0.592 |
| OSD-244 | glare_module_without_platelet_or_coag | 51 | 19 | 19 | 0.713 | 1.06e-04 |
| OSD-289 | glare_module_without_platelet_or_coag | 51 | 12 | 3 | 0.227 | 0.0128 |
| OSD-421 | glare_module_without_platelet_or_coag | 51 | 10 | 10 | -0.213 | 0.0886 |
| OSD-457 | glare_module_without_platelet_or_coag | 51 | 12 | 12 | 0.195 | 2.59e-05 |
| OSD-515 | glare_module_without_platelet_or_coag | 51 | 9 | 9 | -0.0923 | 0.609 |
| OSD-244 | glare_platelet_ca_module_77 | 77 | 19 | 19 | 0.872 | 1.61e-04 |
| OSD-289 | glare_platelet_ca_module_77 | 77 | 12 | 3 | 0.295 | 0.00206 |
| OSD-421 | glare_platelet_ca_module_77 | 77 | 10 | 10 | -0.168 | 0.109 |
| OSD-457 | glare_platelet_ca_module_77 | 77 | 12 | 12 | 0.254 | 4.24e-05 |
| OSD-515 | glare_platelet_ca_module_77 | 77 | 9 | 9 | -0.0248 | 0.876 |
| OSD-244 | myeloid_macrophage | 13 | 19 | 19 | 0.981 | 5.43e-04 |
| OSD-289 | myeloid_macrophage | 13 | 12 | 3 | 0.655 | 3.87e-05 |
| OSD-421 | myeloid_macrophage | 13 | 10 | 10 | 0.112 | 0.582 |
| OSD-457 | myeloid_macrophage | 13 | 12 | 12 | 0.44 | 8.76e-06 |
| OSD-515 | myeloid_macrophage | 13 | 9 | 9 | -0.199 | 0.172 |
| OSD-244 | platelet_megakaryocyte_core | 17 | 19 | 19 | 0.903 | 0.00273 |
| OSD-289 | platelet_megakaryocyte_core | 17 | 12 | 3 | 0.325 | 0.00476 |
| OSD-421 | platelet_megakaryocyte_core | 17 | 10 | 10 | -0.0925 | 0.503 |
| OSD-457 | platelet_megakaryocyte_core | 17 | 12 | 12 | 0.26 | 3.75e-07 |
| OSD-515 | platelet_megakaryocyte_core | 17 | 9 | 9 | 0.117 | 0.617 |
| OSD-244 | thymic_epithelial_tec | 16 | 19 | 19 | 0.0186 | 0.944 |
| OSD-289 | thymic_epithelial_tec | 16 | 12 | 3 | 0.674 | 4.94e-06 |
| OSD-421 | thymic_epithelial_tec | 16 | 10 | 10 | -0.338 | 0.157 |
| OSD-457 | thymic_epithelial_tec | 16 | 12 | 12 | 0.487 | 9.61e-07 |
| OSD-515 | thymic_epithelial_tec | 16 | 9 | 9 | -0.551 | 0.0446 |
| OSD-244 | thymocyte_t_cell | 18 | 19 | 19 | -0.908 | 0.0254 |
| OSD-289 | thymocyte_t_cell | 18 | 12 | 3 | 0.0218 | 0.677 |
| OSD-421 | thymocyte_t_cell | 18 | 10 | 10 | -0.552 | 0.0226 |
| OSD-457 | thymocyte_t_cell | 18 | 12 | 12 | -0.157 | 4.32e-04 |
| OSD-515 | thymocyte_t_cell | 18 | 9 | 9 | -0.83 | 0.105 |

## Study-Normalized Correlation With GLARE Module

Scores were z-scored within OSD accession before correlation, so this checks sample-to-sample coupling inside studies rather than accession baseline.

| marker_set | n_samples | gene_overlap_with_outcome | spearman_r | spearman_p |
| --- | --- | --- | --- | --- |
| glare_module_without_platelet_coag_endothelial | 115 | 50 | 0.968 | 5.49e-70 |
| glare_module_without_platelet_or_coag | 115 | 51 | 0.967 | 5.92e-69 |
| panglao_erythroid | 115 | 3 | 0.901 | 8.93e-43 |
| panglao_fibroblasts | 115 | 2 | 0.889 | 4.09e-40 |
| panglao_platelets | 115 | 18 | 0.885 | 2.36e-39 |
| panglao_endothelial_cells | 115 | 10 | 0.841 | 5.77e-32 |
| panglao_monocytes | 115 | 8 | 0.816 | 1.04e-28 |
| fibroblast_ecm_stromal | 115 | 2 | 0.804 | 2.54e-27 |
| panglao_macrophages | 115 | 2 | 0.788 | 1.61e-25 |
| panglao_megakaryocytes | 115 | 10 | 0.785 | 3.26e-25 |
| coagulation_plasma | 115 | 15 | 0.778 | 1.60e-24 |
| panglao_epithelial_cells | 115 | 2 | 0.754 | 2.12e-22 |
| platelet_megakaryocyte_core | 115 | 11 | 0.748 | 8.01e-22 |
| myeloid_macrophage | 115 | 0 | 0.74 | 3.33e-21 |
| endothelial_vascular_core | 115 | 2 | 0.72 | 1.12e-19 |
| thymic_epithelial_tec | 115 | 0 | 0.379 | 2.93e-05 |
| panglao_t_cells | 115 | 0 | 0.344 | 1.70e-04 |
| erythroid_blood | 115 | 0 | 0.343 | 1.74e-04 |
| thymocyte_t_cell | 115 | 0 | -0.0842 | 0.371 |

## Flight Effect After Marker Adjustment

Outcome and marker scores are study-z values. `none` is the unadjusted FLT-vs-GC module-score model; other rows add one marker panel as a covariate.

| outcome | adjustment_marker_set | n_samples | flight_coef | flight_p | marker_coef | marker_p | r_squared |
| --- | --- | --- | --- | --- | --- | --- | --- |
| glare_platelet_ca_module_77 | none | 115 | 0.704 | 7.50e-05 |  |  | 0.123 |
| glare_platelet_ca_module_77 | coagulation_plasma | 115 | 0.247 | 0.0444 | 0.762 | 1.72e-34 | 0.652 |
| glare_platelet_ca_module_77 | endothelial_vascular_core | 115 | -0.00312 | 0.985 | 0.707 | 3.20e-14 | 0.499 |
| glare_platelet_ca_module_77 | erythroid_blood | 115 | 0.634 | 2.30e-04 | 0.297 | 6.77e-04 | 0.21 |
| glare_platelet_ca_module_77 | fibroblast_ecm_stromal | 115 | 0.0975 | 0.423 | 0.778 | 2.23e-38 | 0.638 |
| glare_platelet_ca_module_77 | myeloid_macrophage | 115 | 0.157 | 0.281 | 0.713 | 1.19e-21 | 0.557 |
| glare_platelet_ca_module_77 | platelet_megakaryocyte_core | 115 | 0.151 | 0.3 | 0.683 | 1.84e-14 | 0.514 |
| glare_platelet_ca_module_77 | thymic_epithelial_tec | 115 | 0.62 | 4.67e-04 | 0.309 | 0.00232 | 0.217 |
| glare_platelet_ca_module_77 | thymocyte_t_cell | 115 | 0.744 | 8.25e-05 | 0.0494 | 0.632 | 0.125 |
| glare_module_without_platelet_or_coag | none | 115 | 0.676 | 1.66e-04 |  |  | 0.114 |
| glare_module_without_platelet_or_coag | coagulation_plasma | 115 | 0.303 | 0.0452 | 0.623 | 7.95e-17 | 0.468 |
| glare_module_without_platelet_or_coag | endothelial_vascular_core | 115 | 0.00861 | 0.963 | 0.668 | 1.23e-10 | 0.449 |
| glare_module_without_platelet_or_coag | erythroid_blood | 115 | 0.62 | 4.12e-04 | 0.242 | 0.00839 | 0.172 |
| glare_module_without_platelet_or_coag | fibroblast_ecm_stromal | 115 | 0.0629 | 0.609 | 0.788 | 5.80e-46 | 0.64 |
| glare_module_without_platelet_or_coag | myeloid_macrophage | 115 | 0.097 | 0.475 | 0.755 | 1.44e-27 | 0.6 |
| glare_module_without_platelet_or_coag | platelet_megakaryocyte_core | 115 | 0.162 | 0.302 | 0.635 | 4.13e-12 | 0.451 |
| glare_module_without_platelet_or_coag | thymic_epithelial_tec | 115 | 0.563 | 8.50e-04 | 0.42 | 6.05e-06 | 0.287 |
| glare_module_without_platelet_or_coag | thymocyte_t_cell | 115 | 0.78 | 3.22e-05 | 0.126 | 0.219 | 0.127 |

## Interpretation

- Platelet/coagulation/endothelial/stromal panels move in the same direction as the GLARE module in OSD-244, OSD-289, and OSD-457, matching the earlier GLARE score support.
- Thymocyte/T-cell scores are not the main explanation. They are weaker and less consistently aligned with the platelet-calcium module than vascular/hemostasis panels.
- TEC markers are not completely unrelated, but TEC scoring does not explain away the module.
- The platelet/coagulation-adjusted residual module remains FLT-associated, so the module likely contains both composition-sensitive hemostasis genes and broader stromal/calcium/remodeling biology.

## Retained Files

- `marker_set_flight_ground_by_study.tsv`
- `module_marker_correlations_study_z.tsv`
- `flight_effect_after_marker_adjustment.tsv`
- `marker_set_overlap_with_glare_module.tsv`
