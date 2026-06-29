# OntoVAE Follow-up Report

This report is generated from completed OntoVAE outputs. It focuses on
strict leave-one-accession-out stable FLT vs GC pathway/program shifts,
then links directly to the PCA/UMAP/heatmap visualizations for manual
inspection.

Primary machine-readable outputs:

- `outputs/ontovae_pipeline/followup/ontovae_stable_followup_terms.tsv`
- `outputs/ontovae_pipeline/followup/ontovae_stable_followup_top_genes.tsv`

Strict LOO-stable means `meta_fdr < 0.05`,
`maximum_leave_one_out_fdr < 0.05`, and all leave-one-accession-out
effects keep the same direction.

## Focus Runs

| tissue | group | mode | stable terms | top stable term | effect | FDR | max LOO FDR |
| --- | --- | --- | ---: | --- | ---: | ---: | ---: |
| liver |  | archs4_pretrain_osdr_finetune | 13 | `R-MMU-5676590_NIK_NONCANONICAL_NF_KB_SIGNALING` | 0.277 | 3.38e-04 | 0.0179 |
| liver |  | direct_osdr | 2 | `R-MMU-3000170_SYNDECAN_INTERACTIONS` | -0.142 | 1.80e-04 | 0.0227 |
| liver |  | hvg_archs4_pretrain_osdr_finetune | 3 | `R-MMU-1266695_INTERLEUKIN_7_SIGNALING` | 0.166 | 4.90e-05 | 0.00445 |
| skeletal_muscle |  | archs4_pretrain_osdr_finetune | 1 | `R-MMU-5576893_PHASE_2_PLATEAU_PHASE` | 0.289 | 2.93e-04 | 0.00669 |
| skeletal_muscle |  | hvg_archs4_pretrain_osdr_finetune | 2 | `R-MMU-73927_DEPURINATION` | 0.125 | 2.60e-06 | 2.09e-04 |
| skeletal_muscle | quadriceps | archs4_pretrain_osdr_finetune | 3 | `R-MMU-399956_CRMPS_IN_SEMA3A_SIGNALING` | 0.544 | 8.96e-05 | 0.0218 |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | 27 | `R-MMU-629597_HIGHLY_CALCIUM_PERMEABLE_NICOTINIC_ACETYLCHOLINE_RECEPTORS` | 0.935 | 7.41e-17 | 0.0101 |
| spleen |  | archs4_pretrain_osdr_finetune | 127 | `R-MMU-194138_SIGNALING_BY_VEGF` | 0.592 | 1.06e-10 | 3.92e-04 |
| spleen |  | direct_osdr | 1 | `R-MMU-8878171_TRANSCRIPTIONAL_REGULATION_BY_RUNX1` | -0.0813 | 4.53e-06 | 0.00343 |
| spleen |  | hvg_archs4_pretrain_osdr_finetune | 7 | `R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES` | 0.353 | 3.26e-05 | 0.0123 |
| thymus |  | archs4_pretrain_osdr_finetune | 730 | `R-MMU-9614085_FOXO_MEDIATED_TRANSCRIPTION` | 0.82 | 2.74e-44 | 1.19e-28 |
| thymus |  | direct_osdr | 114 | `R-MMU-174143_APC_C_MEDIATED_DEGRADATION_OF_CELL_CYCLE_PROTEINS` | -0.821 | 3.22e-42 | 3.81e-27 |
| thymus |  | hvg_archs4_pretrain_osdr_finetune | 146 | `R-MMU-3858494_BETA_CATENIN_INDEPENDENT_WNT_SIGNALING` | -0.831 | 2.28e-42 | 1.51e-18 |

## Priority Pathways

| tissue | group | mode | score set | rank | term | effect | FDR | max LOO FDR | top decoder genes |
| --- | --- | --- | --- | ---: | --- | ---: | ---: | ---: | --- |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | finetuned_or_direct | 1 | `R-MMU-629597_HIGHLY_CALCIUM_PERMEABLE_NICOTINIC_ACETYLCHOLINE_RECEPTORS` | 0.935 | 7.41e-17 | 0.0101 | ENSMUSG00000027950, ENSMUSG00000031492, ENSMUSG00000035200, ENSMUSG00000032303, ENSMUSG00000031491 |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | finetuned_or_direct | 2 | `R-MMU-211979_EICOSANOIDS` | 1.24 | 1.12e-11 | 7.90e-06 | ENSMUSG00000028712, ENSMUSG00000063929, ENSMUSG00000083138, ENSMUSG00000066071, ENSMUSG00000090700 |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | finetuned_or_direct | 3 | `R-MMU-975634_RETINOID_METABOLISM_AND_TRANSPORT` | 1.09 | 1.39e-09 | 0.00522 | ENSMUSG00000027070, ENSMUSG00000024391, ENSMUSG00000028613, ENSMUSG00000028003, ENSMUSG00000040249 |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | finetuned_or_direct | 4 | `R-MMU-5668599_RHO_GTPASES_ACTIVATE_NADPH_OXIDASES` | 1.24 | 1.12e-08 | 1.15e-04 | ENSMUSG00000015950, ENSMUSG00000071715, ENSMUSG00000031257, ENSMUSG00000052889, ENSMUSG00000015340 |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | finetuned_or_direct | 5 | `R-MMU-9955298_SLC_MEDIATED_TRANSPORT_OF_ORGANIC_ANIONS` | 1.7 | 4.44e-08 | 2.31e-04 | ENSMUSG00000020102, ENSMUSG00000018459, ENSMUSG00000063796, ENSMUSG00000021728, ENSMUSG00000024650 |
| liver |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 1 | `R-MMU-5676590_NIK_NONCANONICAL_NF_KB_SIGNALING` | 0.277 | 3.38e-04 | 0.0179 | ENSMUSG00000030061, ENSMUSG00000036309, ENSMUSG00000060073, ENSMUSG00000021832, ENSMUSG00000030603 |
| liver |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 2 | `R-MMU-444473_FORMYL_PEPTIDE_RECEPTORS_BIND_FORMYL_PEPTIDES_AND_MANY_OTHER_LIGANDS` | -0.511 | 3.38e-04 | 0.0457 | ENSMUSG00000019122, ENSMUSG00000042770, ENSMUSG00000018927, ENSMUSG00000052270, ENSMUSG00000045551 |
| liver |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 3 | `R-MMU-205043_NRIF_SIGNALS_CELL_DEATH_FROM_THE_NUCLEUS` | -0.275 | 6.01e-04 | 0.00525 | ENSMUSG00000015750, ENSMUSG00000000120, ENSMUSG00000032375, ENSMUSG00000019969, ENSMUSG00000028549 |
| liver |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 4 | `R-MMU-9706369_NEGATIVE_REGULATION_OF_FLT3` | -0.497 | 6.01e-04 | 0.00948 | ENSMUSG00000027636, ENSMUSG00000022372, ENSMUSG00000110206, ENSMUSG00000020027, ENSMUSG00000034342 |
| liver |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 5 | `R-MMU-5676594_TNF_RECEPTOR_SUPERFAMILY_TNFSF_MEMBERS_MEDIATING_NON_CANONICAL_NF_KB_PATHWAY` | -0.561 | 0.00102 | 0.00749 | ENSMUSG00000097328, ENSMUSG00000026321, ENSMUSG00000030339, ENSMUSG00000031497, ENSMUSG00000017652 |
| skeletal_muscle |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 1 | `R-MMU-5576893_PHASE_2_PLATEAU_PHASE` | 0.289 | 2.93e-04 | 0.00669 | ENSMUSG00000047330, ENSMUSG00000040407, ENSMUSG00000051331, ENSMUSG00000035165, ENSMUSG00000039672 |
| skeletal_muscle | quadriceps | archs4_pretrain_osdr_finetune | finetuned_or_direct | 1 | `R-MMU-399956_CRMPS_IN_SEMA3A_SIGNALING` | 0.544 | 8.96e-05 | 0.0218 | ENSMUSG00000024501, ENSMUSG00000026640, ENSMUSG00000048895, ENSMUSG00000031398, ENSMUSG00000025810 |
| skeletal_muscle | quadriceps | archs4_pretrain_osdr_finetune | finetuned_or_direct | 2 | `R-MMU-179812_GRB2_EVENTS_IN_EGFR_SIGNALING` | -0.33 | 1.79e-04 | 0.0149 | ENSMUSG00000029378, ENSMUSG00000029377, ENSMUSG00000024241, ENSMUSG00000082361, ENSMUSG00000035020 |
| skeletal_muscle | quadriceps | archs4_pretrain_osdr_finetune | finetuned_or_direct | 3 | `R-MMU-5651801_PCNA_DEPENDENT_LONG_PATCH_BASE_EXCISION_REPAIR` | -0.292 | 4.27e-04 | 0.00829 | ENSMUSG00000030042, ENSMUSG00000028394, ENSMUSG00000035960, ENSMUSG00000012483, ENSMUSG00000023104 |
| spleen |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 1 | `R-MMU-194138_SIGNALING_BY_VEGF` | 0.592 | 1.06e-10 | 3.92e-04 | ENSMUSG00000029648, ENSMUSG00000062960, ENSMUSG00000061878, ENSMUSG00000020357, ENSMUSG00000004951 |
| spleen |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 2 | `R-MMU-416700_OTHER_SEMAPHORIN_INTERACTIONS` | -0.722 | 1.30e-10 | 5.47e-05 | ENSMUSG00000028064, ENSMUSG00000021451, ENSMUSG00000074785, ENSMUSG00000026395, ENSMUSG00000031385 |
| spleen |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 3 | `R-MMU-5213460_RIPK1_MEDIATED_REGULATED_NECROSIS` | -0.361 | 1.47e-10 | 0.0164 | ENSMUSG00000000817, ENSMUSG00000028249, ENSMUSG00000039304, ENSMUSG00000020134, ENSMUSG00000026942 |
| spleen |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 4 | `R-MMU-1296071_POTASSIUM_CHANNELS` | 0.829 | 3.73e-10 | 1.36e-04 | ENSMUSG00000000794, ENSMUSG00000058248, ENSMUSG00000045246, ENSMUSG00000038077, ENSMUSG00000035580 |
| spleen |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 5 | `R-MMU-1474290_COLLAGEN_FORMATION` | 1.1 | 5.48e-10 | 4.09e-06 | ENSMUSG00000070436, ENSMUSG00000022371, ENSMUSG00000020674, ENSMUSG00000022098, ENSMUSG00000023191 |
| thymus |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 1 | `R-MMU-9614085_FOXO_MEDIATED_TRANSCRIPTION` | 0.82 | 2.74e-44 | 1.19e-28 | ENSMUSG00000048756, ENSMUSG00000024515, ENSMUSG00000020950, ENSMUSG00000042903, ENSMUSG00000032402 |
| thymus |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 2 | `R-MMU-392851_PROSTACYCLIN_SIGNALLING_THROUGH_PROSTACYCLIN_RECEPTOR` | -0.741 | 5.31e-44 | 1.75e-23 | ENSMUSG00000038607, ENSMUSG00000043004, ENSMUSG00000043017, ENSMUSG00000068523, ENSMUSG00000027523 |
| thymus |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 3 | `R-MMU-428930_THROMBOXANE_SIGNALLING_THROUGH_TP_RECEPTOR` | -1.27 | 7.17e-44 | 7.32e-26 | ENSMUSG00000034881, ENSMUSG00000071658, ENSMUSG00000038811, ENSMUSG00000032192, ENSMUSG00000020611 |
| thymus |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 4 | `R-MMU-5627123_RHO_GTPASES_ACTIVATE_PAKS` | 1.27 | 6.48e-43 | 2.86e-23 | ENSMUSG00000030739, ENSMUSG00000020900, ENSMUSG00000009073, ENSMUSG00000067818, ENSMUSG00000034868 |
| thymus |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 5 | `R-MMU-193639_P75NTR_SIGNALS_VIA_NF_KB` | 1.32 | 1.81e-42 | 4.71e-19 | ENSMUSG00000000120, ENSMUSG00000015837, ENSMUSG00000031392, ENSMUSG00000019505, ENSMUSG00000021025 |

## Plot Review Queue

### skeletal_muscle soleus archs4_pretrain_osdr_finetune finetuned_or_direct

Top stable term: `R-MMU-629597_HIGHLY_CALCIUM_PERMEABLE_NICOTINIC_ACETYLCHOLINE_RECEPTORS` (effect 0.935, FDR 7.41e-17, max LOO FDR 0.0101).

- [top pathway heatmap](../outputs/ontovae_skeletal_muscle_splits/soleus/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png) | [UMAP by accession](../outputs/ontovae_skeletal_muscle_splits/soleus/archs4_pretrain_osdr_finetune/analysis/pathway_score_umap_by_accession.png) | [PCA by accession](../outputs/ontovae_skeletal_muscle_splits/soleus/archs4_pretrain_osdr_finetune/analysis/pathway_score_pca_by_accession.png)

![Top pathway heatmap](../outputs/ontovae_skeletal_muscle_splits/soleus/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png)

### liver  archs4_pretrain_osdr_finetune finetuned_or_direct

Top stable term: `R-MMU-5676590_NIK_NONCANONICAL_NF_KB_SIGNALING` (effect 0.277, FDR 3.38e-04, max LOO FDR 0.0179).

- [top pathway heatmap](../outputs/ontovae_liver/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png) | [UMAP by accession](../outputs/ontovae_liver/archs4_pretrain_osdr_finetune/analysis/pathway_score_umap_by_accession.png) | [PCA by accession](../outputs/ontovae_liver/archs4_pretrain_osdr_finetune/analysis/pathway_score_pca_by_accession.png)

![Top pathway heatmap](../outputs/ontovae_liver/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png)

### skeletal_muscle  archs4_pretrain_osdr_finetune finetuned_or_direct

Top stable term: `R-MMU-5576893_PHASE_2_PLATEAU_PHASE` (effect 0.289, FDR 2.93e-04, max LOO FDR 0.00669).

- [top pathway heatmap](../outputs/ontovae_skeletal_muscle/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png) | [UMAP by accession](../outputs/ontovae_skeletal_muscle/archs4_pretrain_osdr_finetune/analysis/pathway_score_umap_by_accession.png) | [PCA by accession](../outputs/ontovae_skeletal_muscle/archs4_pretrain_osdr_finetune/analysis/pathway_score_pca_by_accession.png)

![Top pathway heatmap](../outputs/ontovae_skeletal_muscle/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png)

### skeletal_muscle quadriceps archs4_pretrain_osdr_finetune finetuned_or_direct

Top stable term: `R-MMU-399956_CRMPS_IN_SEMA3A_SIGNALING` (effect 0.544, FDR 8.96e-05, max LOO FDR 0.0218).

- [top pathway heatmap](../outputs/ontovae_skeletal_muscle_splits/quadriceps/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png) | [UMAP by accession](../outputs/ontovae_skeletal_muscle_splits/quadriceps/archs4_pretrain_osdr_finetune/analysis/pathway_score_umap_by_accession.png) | [PCA by accession](../outputs/ontovae_skeletal_muscle_splits/quadriceps/archs4_pretrain_osdr_finetune/analysis/pathway_score_pca_by_accession.png)

![Top pathway heatmap](../outputs/ontovae_skeletal_muscle_splits/quadriceps/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png)

### spleen  archs4_pretrain_osdr_finetune finetuned_or_direct

Top stable term: `R-MMU-194138_SIGNALING_BY_VEGF` (effect 0.592, FDR 1.06e-10, max LOO FDR 3.92e-04).

- [top pathway heatmap](../outputs/ontovae_spleen/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png) | [UMAP by accession](../outputs/ontovae_spleen/archs4_pretrain_osdr_finetune/analysis/pathway_score_umap_by_accession.png) | [PCA by accession](../outputs/ontovae_spleen/archs4_pretrain_osdr_finetune/analysis/pathway_score_pca_by_accession.png)

![Top pathway heatmap](../outputs/ontovae_spleen/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png)

### thymus  archs4_pretrain_osdr_finetune finetuned_or_direct

Top stable term: `R-MMU-9614085_FOXO_MEDIATED_TRANSCRIPTION` (effect 0.82, FDR 2.74e-44, max LOO FDR 1.19e-28).

- [top pathway heatmap](../outputs/ontovae_thymus/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png) | [UMAP by accession](../outputs/ontovae_thymus/archs4_pretrain_osdr_finetune/analysis/pathway_score_umap_by_accession.png) | [PCA by accession](../outputs/ontovae_thymus/archs4_pretrain_osdr_finetune/analysis/pathway_score_pca_by_accession.png)

![Top pathway heatmap](../outputs/ontovae_thymus/archs4_pretrain_osdr_finetune/analysis/top_pathway_shift_heatmap.png)

## Frozen Projection vs Fine-tuned Scores

Rows with `score_set = pre_finetune_projection` are OSDR samples scored
by the ARCHS4-pretrained OntoVAE before OSDR fine-tuning. Rows with
`finetuned_or_direct` are final fine-tuned or direct OSDR scores.

| tissue | group | mode | score set | stable terms | top stable term | effect | FDR |
| --- | --- | --- | --- | ---: | --- | ---: | ---: |
| liver |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 13 | `R-MMU-5676590_NIK_NONCANONICAL_NF_KB_SIGNALING` | 0.277 | 3.38e-04 |
| liver |  | archs4_pretrain_osdr_finetune | pre_finetune_projection | 48 | `R-MMU-190828_GAP_JUNCTION_TRAFFICKING` | -0.223 | 1.25e-06 |
| liver |  | direct_osdr | finetuned_or_direct | 2 | `R-MMU-3000170_SYNDECAN_INTERACTIONS` | -0.142 | 1.80e-04 |
| liver |  | hvg_archs4_pretrain_osdr_finetune | finetuned_or_direct | 3 | `R-MMU-1266695_INTERLEUKIN_7_SIGNALING` | 0.166 | 4.90e-05 |
| liver |  | hvg_archs4_pretrain_osdr_finetune | pre_finetune_projection | 3 | `R-MMU-1266695_INTERLEUKIN_7_SIGNALING` | 0.165 | 4.25e-05 |
| skeletal_muscle |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 1 | `R-MMU-5576893_PHASE_2_PLATEAU_PHASE` | 0.289 | 2.93e-04 |
| skeletal_muscle |  | archs4_pretrain_osdr_finetune | pre_finetune_projection | 2 | `R-MMU-8878171_TRANSCRIPTIONAL_REGULATION_BY_RUNX1` | 0.197 | 8.33e-06 |
| skeletal_muscle |  | hvg_archs4_pretrain_osdr_finetune | finetuned_or_direct | 2 | `R-MMU-73927_DEPURINATION` | 0.125 | 2.60e-06 |
| skeletal_muscle | quadriceps | archs4_pretrain_osdr_finetune | finetuned_or_direct | 3 | `R-MMU-399956_CRMPS_IN_SEMA3A_SIGNALING` | 0.544 | 8.96e-05 |
| skeletal_muscle | quadriceps | archs4_pretrain_osdr_finetune | pre_finetune_projection | 1 | `R-MMU-179812_GRB2_EVENTS_IN_EGFR_SIGNALING` | -0.134 | 0.00155 |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | finetuned_or_direct | 27 | `R-MMU-629597_HIGHLY_CALCIUM_PERMEABLE_NICOTINIC_ACETYLCHOLINE_RECEPTORS` | 0.935 | 7.41e-17 |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | pre_finetune_projection | 7 | `R-MMU-211981_XENOBIOTICS` | 0.57 | 7.30e-06 |
| spleen |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 127 | `R-MMU-194138_SIGNALING_BY_VEGF` | 0.592 | 1.06e-10 |
| spleen |  | archs4_pretrain_osdr_finetune | pre_finetune_projection | 129 | `R-MMU-1296071_POTASSIUM_CHANNELS` | 0.847 | 3.50e-13 |
| spleen |  | direct_osdr | finetuned_or_direct | 1 | `R-MMU-8878171_TRANSCRIPTIONAL_REGULATION_BY_RUNX1` | -0.0813 | 4.53e-06 |
| spleen |  | hvg_archs4_pretrain_osdr_finetune | finetuned_or_direct | 7 | `R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES` | 0.353 | 3.26e-05 |
| spleen |  | hvg_archs4_pretrain_osdr_finetune | pre_finetune_projection | 9 | `R-MMU-1566948_ELASTIC_FIBRE_FORMATION` | 0.397 | 6.10e-05 |
| thymus |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 730 | `R-MMU-9614085_FOXO_MEDIATED_TRANSCRIPTION` | 0.82 | 2.74e-44 |
| thymus |  | archs4_pretrain_osdr_finetune | pre_finetune_projection | 797 | `R-MMU-6781823_FORMATION_OF_TC_NER_PRE_INCISION_COMPLEX` | -0.531 | 7.06e-51 |
| thymus |  | direct_osdr | finetuned_or_direct | 114 | `R-MMU-174143_APC_C_MEDIATED_DEGRADATION_OF_CELL_CYCLE_PROTEINS` | -0.821 | 3.22e-42 |
| thymus |  | hvg_archs4_pretrain_osdr_finetune | finetuned_or_direct | 146 | `R-MMU-3858494_BETA_CATENIN_INDEPENDENT_WNT_SIGNALING` | -0.831 | 2.28e-42 |
| thymus |  | hvg_archs4_pretrain_osdr_finetune | pre_finetune_projection | 139 | `R-MMU-195258_RHO_GTPASE_EFFECTORS` | -0.471 | 2.77e-41 |

## Caveats

- Decoder genes are listed as Ensembl mouse gene IDs because the OntoVAE
  AnnData inputs do not carry a project-wide gene-symbol annotation.
- These are pathway/program score shifts, not direct gene-level DGEA.
- OntoVAE here uses ARCHS4 pretraining plus OSDR fine-tuning; it is not
  native scArches query mapping.
- Skin, kidney, EDL, gastrocnemius, and tibialis anterior had exploratory
  random-effects hits but no strict LOO-stable final hits.
