# TMS single-cell check: thymus platelet-calcium module

Input: local Tabula Muris Senis FACS matrix `data/processed/tms_facs_3552_cells.matrix.npz`. Scores are mean log1p(CPM/10k) over module gene sets.

Gene-set coverage:
- full_platelet_ca_module: 77/77 genes present in TMS
- platelet_coag_subset: 36/36 genes present in TMS
- endothelial_ecm_subset: 22/22 genes present in TMS
- residual_no_platelet_coag: 41/41 genes present in TMS

Main interpretation:
- The full module and residual module score higher in non-thymocyte classes than in T/thymocyte or thymic epithelial classes.
- The endothelial/ECM subset is highest in endothelial and stromal/fibroblast classes.
- Within TMS thymus cells, fibroblasts score highest, thymocytes are lower, and thymic epithelial cells are lowest; thymus fibroblast and epithelial counts are small.
- This supports a composition/stromal-vascular interpretation more than a thymocyte-intrinsic interpretation. It does not prove blood contamination, because TMS is a normal reference and lacks spaceflight condition labels.

Files:
- `module_set_coverage.tsv`
- `tms_cell_type_module_scores.tsv`
- `tms_broad_class_module_scores.tsv`
- `tms_thymus_cell_type_module_scores.tsv`
- `module_gene_tms_top_expression.tsv`
