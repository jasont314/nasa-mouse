# GSE240016 Thymic Endothelial scRNA Check

Source: GSE240016 `CD45neg_thymic_stroma_d0_EC+annotation.h5ad`, 1,661 day-0 thymic endothelial cells.

## Bottom Line

This independent mouse thymic endothelial dataset supports the endothelial/stromal side of the GLARE thymus platelet-calcium module.
The full GLARE module and the platelet/coagulation subsets are expressed in bona fide annotated thymic endothelial cells, especially the venous EC subset, so the signal is not automatically proof of accidental platelet contamination.
The module also correlates with total counts/n_genes in this h5ad, so subtype-specific interpretation still needs care.
It still remains composition-sensitive because the module also overlaps heavily with platelet/coagulation genes.

## Coverage

| module_set | symbols_requested | symbols_found | missing_symbols |
| --- | --- | --- | --- |
| coagulation_plasma | 19 | 8 | A2M,ALB,F13A1,FGA,FGB,HRG,KNG1,PLG,SERPINA1,SERPINA3,TF |
| endothelial_vascular_core | 15 | 15 |  |
| erythroid_blood | 10 | 6 | GATA1,GYPA,KLF1,SLC4A1 |
| fibroblast_ecm_stromal | 17 | 17 |  |
| glare_module_without_platelet_coag_endothelial | 54 | 42 | HGF,LEFTY2,LOC128192,LOC131691,LOC643997,LOC652612,LOC654188,LOC727848,POTEM,PPIAP22,PPIAP8,STX4 |
| glare_module_without_platelet_or_coag | 55 | 43 | HGF,LEFTY2,LOC128192,LOC131691,LOC643997,LOC652612,LOC654188,LOC727848,POTEM,PPIAP22,PPIAP8,STX4 |
| glare_platelet_ca_module_77 | 77 | 57 | A2M,ALB,FGA,FGB,HGF,HRG,LEFTY2,LOC128192,LOC131691,LOC643997,LOC652612,LOC654188,LOC727848,PLG,POTEM,PPIAP22,PPIAP8,SERPINA1,STX4,TF |
| myeloid_macrophage | 13 | 11 | CSF1R,ITGAM |
| platelet_megakaryocyte_core | 17 | 17 |  |
| thymic_epithelial_tec | 16 | 14 | FOXN1,PSMB11 |
| thymocyte_t_cell | 18 | 14 | CD8A,IL7R,RAG2,ZAP70 |

## Scores By Endothelial Subtype

| group | module_set | n_cells | mean_score | median_score | p90_score | positive_fraction |
| --- | --- | --- | --- | --- | --- | --- |
| 2:venEC | coagulation_plasma | 2.51e+02 | 0.183 | 0.162 | 0.415 | 0.625 |
| 0:arEC | coagulation_plasma | 3.25e+02 | 0.075 | 0 | 0.264 | 0.308 |
| 1:capEC | coagulation_plasma | 1.08e+03 | 0.0642 | 0 | 0.239 | 0.279 |
| 2:venEC | endothelial_vascular_core | 2.51e+02 | 1.75 | 1.81 | 2.1 | 0.996 |
| 1:capEC | endothelial_vascular_core | 1.08e+03 | 1.72 | 1.73 | 2.09 | 0.999 |
| 0:arEC | endothelial_vascular_core | 3.25e+02 | 1.67 | 1.69 | 2.02 | 1 |
| 0:arEC | fibroblast_ecm_stromal | 3.25e+02 | 0.208 | 0.191 | 0.361 | 0.88 |
| 2:venEC | fibroblast_ecm_stromal | 2.51e+02 | 0.206 | 0.187 | 0.366 | 0.892 |
| 1:capEC | fibroblast_ecm_stromal | 1.08e+03 | 0.191 | 0.189 | 0.286 | 0.936 |
| 2:venEC | glare_module_without_platelet_or_coag | 2.51e+02 | 0.792 | 0.819 | 0.983 | 1 |
| 0:arEC | glare_module_without_platelet_or_coag | 3.25e+02 | 0.771 | 0.79 | 0.969 | 1 |
| 1:capEC | glare_module_without_platelet_or_coag | 1.08e+03 | 0.727 | 0.726 | 0.899 | 1 |
| 2:venEC | glare_platelet_ca_module_77 | 2.51e+02 | 0.714 | 0.726 | 0.881 | 1 |
| 0:arEC | glare_platelet_ca_module_77 | 3.25e+02 | 0.639 | 0.642 | 0.801 | 1 |
| 1:capEC | glare_platelet_ca_module_77 | 1.08e+03 | 0.593 | 0.59 | 0.731 | 1 |
| 2:venEC | platelet_megakaryocyte_core | 2.51e+02 | 0.555 | 0.567 | 0.792 | 0.988 |
| 0:arEC | platelet_megakaryocyte_core | 3.25e+02 | 0.393 | 0.386 | 0.556 | 0.991 |
| 1:capEC | platelet_megakaryocyte_core | 1.08e+03 | 0.356 | 0.361 | 0.495 | 0.992 |
| 0:arEC | thymic_epithelial_tec | 3.25e+02 | 0.0812 | 0.0757 | 0.208 | 0.511 |
| 2:venEC | thymic_epithelial_tec | 2.51e+02 | 0.0693 | 0.0629 | 0.171 | 0.514 |
| 1:capEC | thymic_epithelial_tec | 1.08e+03 | 0.0645 | 0 | 0.176 | 0.412 |
| 2:venEC | thymocyte_t_cell | 2.51e+02 | 0.00652 | 0 | 0 | 0.0518 |
| 0:arEC | thymocyte_t_cell | 3.25e+02 | 0.005 | 0 | 0 | 0.0462 |
| 1:capEC | thymocyte_t_cell | 1.08e+03 | 0.00365 | 0 | 0 | 0.0341 |

## Correlation With EC Metadata Scores

| metadata_score | n_cells | spearman_r | spearman_p |
| --- | --- | --- | --- |
| total_counts | 1.66e+03 | 0.654 | 3.97e-203 |
| n_genes_by_counts | 1.66e+03 | 0.642 | 1.45e-193 |
| ribo_frac | 1.66e+03 | 0.295 | 8.96e-35 |
| venular | 1.66e+03 | 0.203 | 5.81e-17 |
| arterial | 1.66e+03 | 0.102 | 3.21e-05 |
| mito_frac | 1.66e+03 | 0.0265 | 0.28 |
| lymphatic | 1.66e+03 | -0.167 | 6.60e-12 |
| capilary | 1.66e+03 | -0.309 | 5.14e-38 |
| hb_frac | 1.66e+03 |  |  |

## Top GLARE Module Genes Expressed In EC Cells

| symbol | mean_expression | positive_fraction | p90_expression |
| --- | --- | --- | --- |
| TMSB4X | 4.81 | 0.996 | 5.54 |
| PPIA | 2.95 | 0.96 | 3.73 |
| CALM1 | 2.84 | 0.945 | 3.71 |
| PECAM1 | 2.48 | 0.902 | 3.5 |
| CFL1 | 1.97 | 0.827 | 3 |
| PSAP | 1.82 | 0.801 | 2.91 |
| PFN1 | 1.79 | 0.778 | 2.89 |
| SRGN | 1.76 | 0.763 | 3 |
| APP | 1.56 | 0.729 | 2.72 |
| CALM2 | 1.46 | 0.692 | 2.7 |
| ALDOA | 1.29 | 0.636 | 2.53 |
| ACTN4 | 1.06 | 0.546 | 2.42 |
| HSPA5 | 0.995 | 0.509 | 2.42 |
| TGFB1 | 0.706 | 0.396 | 2.1 |
| SOD1 | 0.699 | 0.397 | 2.04 |
| CLU | 0.68 | 0.279 | 2.79 |
| VCL | 0.637 | 0.328 | 2.22 |
| PDGFB | 0.601 | 0.341 | 1.96 |
| VWF | 0.531 | 0.217 | 2.42 |
| TLN1 | 0.518 | 0.31 | 1.84 |
