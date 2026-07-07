# External single-cell check: GSE103967 thymus stroma WT

Dataset: GSE103967, Large-scale single-cell mapping of the thymic stroma. Used filtered `thymus_stroma_WT` cells only.
- Filtered WT stromal cells scored: 2021
- Broad classes are inferred from FACS markers: CD31 endothelial-like, CD34 mesenchymal-like, EPCAM TEC-like, or low/mixed.

Coverage:
- full_platelet_ca_module: 63/77 symbols found
- platelet_coag_subset: 30/36 symbols found
- endothelial_ecm_subset: 21/22 symbols found
- residual_no_platelet_coag: 33/41 symbols found

Top FACS class by module set:
- full_platelet_ca_module: CD31_endothelial_like mean=0.042, n=115
- platelet_coag_subset: CD31_endothelial_like mean=0.041, n=115
- endothelial_ecm_subset: CD31_endothelial_like mean=0.038, n=115
- residual_no_platelet_coag: CD31_endothelial_like mean=0.043, n=115

Interpretation: this external thymic stroma reference supports a stromal/vascular interpretation if endothelial/CD34 classes score above EPCAM/low-marker classes, and it undermines a pure thymocyte-intrinsic explanation because the dataset is stromal and still captures module expression.
