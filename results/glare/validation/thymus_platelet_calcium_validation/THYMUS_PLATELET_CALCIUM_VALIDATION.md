# Thymus Platelet-Calcium Module Validation

## Bottom Line

The thymus `Response To Elevated Platelet Cytosolic Ca2` GLARE-only signal is real enough to inspect, but it is not clean evidence for thymocyte-intrinsic biology.
It is best labeled as a recurring FLT-up platelet/hemostasis/endothelial-remodeling program with substantial composition risk.

## Reproducibility Notes

- The per-study GLARE enrichment uses the Reactome v4 mouse Ensembl term with 85 mapped IDs.
- Across all cluster rows, the observed thymus term union has 82 IDs; across FDR-significant cluster rows, the union has 77 IDs.
- The DGEA gene summaries and sample-level Reactome score below use the FDR-significant cluster union.
- Symbols are mapped from the paired Reactome v4 symbol and mouse Ensembl GMT files. The OSDR count matrices remain Ensembl-keyed.
- `log2FoldChange > 0` means higher in spaceflight.

## Study-Level Recurrence

- Significant GLARE Reactome enrichment appears in 5 accessions: OSD-244, OSD-289, OSD-421, OSD-457, OSD-515.
- Strict module-score FLT-up support appears in 3 accessions: OSD-244, OSD-289, OSD-457.
- Strict module-score GC-up support appears in 0 accessions: none.
- OSD-421 trends GC-higher by module score but is not strict-significant; OSD-515 is effectively flat.

## Composition Checks

Panglao enrichment for the same module strongly favors platelet/megakaryocyte/endothelial categories:

| panglao_term | overlap | fdr_bh |
| --- | --- | --- |
| panglao_Platelets | 18 | 2.732167446173337e-23 |
| panglao_Megakaryocytes | 10 | 5.184364944980144e-15 |
| panglao_Endothelial_cells | 11 | 6.407698850117057e-10 |
| panglao_Hepatocytes | 10 | 2.408851226369631e-09 |
| panglao_Luminal_epithelial_cells | 7 | 1.5950918501330028e-08 |
| panglao_Monocytes | 8 | 4.52308973033116e-08 |

Manual marker-category counts among study-gene DGEA rows:

| category | genes | tested_study_gene_pairs | sig_pairs | flt_up_sig_pairs | gc_up_sig_pairs | median_log2fc_tested | top_symbols |
| --- | --- | --- | --- | --- | --- | --- | --- |
| panglao_platelet | 18 | 85 | 23 | 21 | 2 | 0.396580004698097 | ACTN1,ACTN4,CALM2,CFD,F8,FGA,HRG,ITGB3,LEFTY2,LOC643997,LOC652612,LOC727848,PECAM1,STXBP3,TGFB3,TMSB4X,TUBA4A,VWF |
| coagulation_plasma | 14 | 66 | 23 | 14 | 9 | 0.19583415647902802 | A2M,ALB,APOA1,F8,FGA,FGB,FGG,HRG,PLG,PROS1,SERPINA1,SERPINF2,SERPING1,TF |
| ecm_growth_factor_remodeling | 12 | 54 | 21 | 19 | 2 | 0.445566076694305 | EGF,FN1,HGF,PDGFA,PDGFB,TGFB1,TGFB2,TGFB3,THBS1,TIMP1,VEGFA,VEGFB |
| calcium_secretion_signaling | 11 | 47 | 11 | 11 | 0 | 0.43863699277221 | ABCC4,CALM1,CALM2,CALM3,PLEK,PRKCA,PRKCB,PRKCG,SRGN,STX4,STXBP3 |
| endothelial_vascular | 10 | 46 | 16 | 14 | 2 | 0.394787467318104 | EGF,FN1,HGF,PDGFA,PDGFB,PECAM1,THBS1,VEGFA,VEGFB,VWF |
| panglao_endothelial | 10 | 45 | 16 | 15 | 1 | 0.457541739572633 | ACTN4,BRPF3,CALM2,F8,FN1,HRG,ITGB3,LEFTY2,PECAM1,TGFB3 |
| panglao_megakaryocyte | 10 | 50 | 14 | 12 | 2 | 0.257956566721755 | ACTN4,APOA1,CALM2,CFD,FGA,LEFTY2,LOC643997,TGFB3,TMSB4X,VWF |
| panglao_monocyte | 8 | 32 | 9 | 8 | 1 | 0.5330584639898455 | ACTN4,F8,HABP4,PPIA,PRKCA,PRKCB,SERPINF2,SRGN |
| platelet_megakaryocyte_core | 8 | 30 | 8 | 8 | 0 | 0.38051191927622796 | ITGA2B,ITGB3,MMRN1,PLEK,PPBP,SELP,SRGN,VWF |
| panglao_pericyte | 3 | 15 | 5 | 4 | 1 | 0.374366323642704 | ACTN4,LOC643997,THBS1 |
| panglao_erythroid | 3 | 15 | 3 | 3 | 0 | 0.555846070772466 | F8,LOC727848,PDGFB |
| erythroid_blood | 0 | 0 | 0 | 0 | 0 |  |  |
| thymic_epithelial_stromal | 0 | 0 | 0 | 0 | 0 |  |  |
| thymocyte_t_cell | 0 | 0 | 0 | 0 | 0 |  |  |

The relevant sample-level marker scores are not uniformly FLT-up across every study, but the platelet/megakaryocyte/endothelial marker sets move FLT-up in the strongest score-supported accessions.
That points to either a real vascular/hemostasis response, a blood/platelet composition shift, or both.

Sample-level mean marker-score shifts, FLT minus GC:

| accession | marker_set | genes_scored | flight_minus_ground |
| --- | --- | --- | --- |
| OSD-244 | panglao_endothelial | 11 | 1.0678362846374512 |
| OSD-289 | panglao_endothelial | 11 | 0.4352550506591797 |
| OSD-421 | panglao_endothelial | 11 | -0.028009891510009766 |
| OSD-457 | panglao_endothelial | 11 | 0.3251938819885254 |
| OSD-515 | panglao_endothelial | 11 | 0.2477421760559082 |
| OSD-244 | panglao_megakaryocyte | 10 | 0.8520774841308594 |
| OSD-289 | panglao_megakaryocyte | 10 | 0.3517308235168457 |
| OSD-421 | panglao_megakaryocyte | 10 | -0.13929080963134766 |
| OSD-457 | panglao_megakaryocyte | 10 | 0.25762033462524414 |
| OSD-515 | panglao_megakaryocyte | 10 | -0.01998138427734375 |
| OSD-244 | panglao_platelet | 18 | 1.09151029586792 |
| OSD-289 | panglao_platelet | 18 | 0.41201186180114746 |
| OSD-421 | panglao_platelet | 18 | -0.08781051635742188 |
| OSD-457 | panglao_platelet | 18 | 0.3018608093261719 |
| OSD-515 | panglao_platelet | 18 | 0.16361713409423828 |
| OSD-244 | reactome_v4_module | 77 | 0.8719925880432129 |
| OSD-289 | reactome_v4_module | 77 | 0.2948281764984131 |
| OSD-421 | reactome_v4_module | 77 | -0.12824678421020508 |
| OSD-457 | reactome_v4_module | 77 | 0.25427675247192383 |
| OSD-515 | reactome_v4_module | 77 | -0.015921592712402344 |

## Top Recurrent Genes

| symbol | gene_id | categories | flt_up_sig_studies | gc_up_sig_studies | median_log2fc | min_padj | accessions_flt_up_sig | accessions_gc_up_sig |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PROS1 | ENSMUSG00000022037 | coagulation_plasma | 4 | 0 | 0.868330039410038 | 1.24355160303872e-08 | OSD-244,OSD-289,OSD-421,OSD-457 |  |
| TIMP1 | ENSMUSG00000025351 | ecm_growth_factor_remodeling | 4 | 0 | 0.778673109996907 | 9.11612533901085e-09 | OSD-244,OSD-289,OSD-421,OSD-457 |  |
| PRKCG | ENSMUSG00000022892 | calcium_secretion_signaling | 4 | 0 | 0.726129247630576 | 1.16044043396369e-08 | OSD-244,OSD-289,OSD-421,OSD-457 |  |
| PDGFA | ENSMUSG00000021253 | ecm_growth_factor_remodeling;endothelial_vascular | 3 | 0 | 0.782845200388497 | 0.0073813300524223 | OSD-244,OSD-289,OSD-457 |  |
| TGFB3 | ENSMUSG00000030342 | ecm_growth_factor_remodeling;panglao_endothelial;panglao_megakaryocyte;panglao_platelet | 3 | 0 | 0.781244934171195 | 7.04689593512989e-05 | OSD-244,OSD-289,OSD-457 |  |
| TLN1 | ENSMUSG00000032554 |  | 3 | 0 | 0.58667287188732 | 7.5186536382837e-10 | OSD-244,OSD-289,OSD-457 |  |
| MMRN1 | ENSMUSG00000022875 | platelet_megakaryocyte_core | 2 | 0 | 2.397936331756635 | 0.0015998838484932 | OSD-244,OSD-457 |  |
| TF | ENSMUSG00000030111 | coagulation_plasma | 2 | 0 | 2.133087851333435 | 0.0032499834884066 | OSD-244,OSD-515 |  |
| FN1 | ENSMUSG00000031520 | ecm_growth_factor_remodeling;endothelial_vascular;panglao_endothelial | 2 | 0 | 1.64907233621552 | 0.0044546142684789 | OSD-244,OSD-515 |  |
| LOC128192 | ENSMUSG00000078816 |  | 2 | 0 | 1.481980033032905 | 0.0166417554766475 | OSD-244,OSD-515 |  |
| CLU | ENSMUSG00000001131 |  | 2 | 0 | 0.806504511261159 | 0.0119442455532848 | OSD-244,OSD-457 |  |
| THBS1 | ENSMUSG00000023224 | ecm_growth_factor_remodeling;endothelial_vascular;panglao_pericyte | 2 | 0 | 0.684381267538106 | 1.24922832413484e-07 | OSD-244,OSD-457 |  |
| FGA | ENSMUSG00000020120 | coagulation_plasma;panglao_megakaryocyte;panglao_platelet | 2 | 0 | 0.447562036892264 | 4.36448240092942e-07 | OSD-244,OSD-457 |  |
| PECAM1 | ENSMUSG00000018593 | endothelial_vascular;panglao_endothelial;panglao_platelet | 2 | 0 | 0.434838156119478 | 0.0071834363454702 | OSD-244,OSD-289 |  |
| PLG | ENSMUSG00000022912 | coagulation_plasma | 2 | 0 | 0.397741199252678 | 3.25367995445662e-05 | OSD-244,OSD-457 |  |

## Interpretation

- Evidence for a one-accession artifact is weak: the pathway recurs across the thymus accessions, and three accessions have strict FLT-up module-score support.
- Evidence for a composition/annotation effect is strong: the module is enriched for platelet, megakaryocyte, endothelial, coagulation/plasma, and vascular-remodeling markers, while canonical thymic epithelial and thymocyte marker sets are absent from the module.
- Evidence for thymus-intrinsic biology is therefore indirect. The signal could reflect thymic vascular remodeling or stromal response, but bulk RNA-seq cannot separate that from platelet/blood contamination or vascular-cell proportion shifts.

## Recommended Follow-Up

1. Treat this as `FLT-up hemostasis/platelet-calcium/endothelial-remodeling`, not simply `thymus platelet activation`.
2. Re-run the module score after removing platelet-core and plasma/coagulation genes to see whether the FLT-up signal persists.
3. Check sample-level platelet/endothelial marker scores against sample metadata and outliers before using this as a main biological claim.
4. Prefer this module as a hypothesis-generating GLARE-only finding unless confirmed by cell-composition deconvolution or independent thymus histology/flow/cell-type data.

## Output Files

- `significant_cluster_marker_summary.tsv`
- `module_gene_dgea_by_study.tsv`
- `module_gene_dgea_summary.tsv`
- `marker_category_summary.tsv`
- `module_score_by_study.tsv`
- `sample_marker_scores.tsv`
- `panglao_marker_enrichment.tsv`
