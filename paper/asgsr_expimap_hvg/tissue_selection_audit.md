# Retrospective tissue-selection audit

## Scope

The repository contains complete direct, reference-query, reference-query with de novo nodes, 2,000-HVG reference-query, and 2,000-HVG reference-query with de novo nodes for eight broad tissues: liver, aggregate skeletal muscle, skin, kidney, thymus, spleen, lung, and retina. Soleus was added later by splitting aggregate skeletal muscle. The manuscript's tissue set was therefore selected retrospectively from a broader screen.

The table below summarizes the original seed-2020 2,000-HVG reference-query outputs. These counts describe different statistical summaries and must not be treated as the final five-check evidence status.

| Tissue | Query samples | Accessions | Pooled Welch FDR < 0.05 | Random-effects FDR < 0.05 | Strict leave-one-accession-out FDR-stable |
| --- | ---: | ---: | ---: | ---: | ---: |
| Liver, original 12-accession input | 231 | 12 | 34 | 74 | 33 |
| Aggregate skeletal muscle | 191 | 13 | 0 | 5 | 1 |
| Skin | 151 | 6 | 94 | 2 | 0 |
| Kidney | 135 | 6 | 0 | 66 | 25 |
| Thymus | 117 | 5 | 174 | 157 | 74 |
| Spleen | 109 | 6 | 0 | 30 | 0 |
| Lung | 78 | 3 | 0 | 0 | 0 |
| Retina | 76 | 4 | 0 | 0 | 0 |
| Soleus split | 53 | 3 | not used here | 22 | 5 |

Strict leave-one-accession-out stability requires the same direction in every omitted-accession fit and maximum leave-one-out FDR below 0.05. It is retained here only to reconstruct the historical screen; the manuscript's current evidence framework uses held-out-project direction as a sensitivity check rather than maximum leave-one-out FDR as a discovery gate.

## Interpretation

- **Lung and retina:** no pooled, random-effects, or strict leave-one-out pathway signal in the 2,000-HVG runs. These are defensible negative screens.
- **Aggregate skeletal muscle:** weak pathway signal with very small effects. Muscle-type splitting was biologically reasonable, but the later soleus model has no pathway that passes the manuscript's five checks.
- **Spleen:** the historical result was neither null nor reliable enough to interpret. Thirty random-effects pathways passed FDR, but none passed the historical strict leave-one-out rule, and HVG selection had silently fallen back to an unbatched ranking. The corrected model described below supersedes that output.
- **Kidney:** the historical result was not safely classifiable as junk, but it used only 1,000 of 2,464 eligible ARCHS4 kidney profiles. Its broad metabolic, platelet, neuronal, keratinization, and drug-label terms should not be carried forward. The corrected model described below supersedes that output.
- **Thymus, skin, and liver:** each retains at least one pathway after the current five-check review and therefore fits a results-focused manuscript.
- **Soleus:** was advanced because muscle-type specificity and prior biology made aggregate muscle an inadequate test. After uniform robustness review, it no longer meets a results-based main-tissue criterion.

## Corrected kidney and spleen reassessment

Kidney was rebuilt with all 2,464 eligible ARCHS4 reference samples. Spleen retained all 6,289 reference samples while 19 singleton ARCHS4 series were excluded from HVG ranking only, allowing genuine batch-aware HVG selection without discarding those profiles from training. Both references used approximately 2,000 HVGs, three complete model trainings, conventional ssGSEA and preranked GSEA, held-out-project direction checks, broad composition-proxy adjustment, and member-gene/decoder-weight review. Spleen OSD-288 was excluded from the primary contrast because recorded strain was disjoint by condition.

### Spleen

The corrected spleen model has a clear lower immune-activation pattern across five unconfounded projects. T-cell receptor signaling, neutrophil degranulation, and C-type lectin receptor signaling are lower in all three complete trainings and all five projects, retain direction after composition-proxy adjustment, agree with ssGSEA, and pass preranked-GSEA FDR below 0.05. Their GSEA FDR values are 0.042, <0.001, and 0.017, respectively. Measured member genes and decoder-weighted gene directions support all three programs. Lower T-cell receptor signaling aligns with prior splenic T-cell impairment; lower C-type lectin-receptor and degranulation programs add a complementary innate pathogen-sensing and effector layer. The latter are transcriptomic programs, not direct functional assays. Lower CD28-linked activation and lymphoid-to-nonlymphoid interaction scores remain supporting context because preranked GSEA points in the opposite direction.

Spleen therefore meets the manuscript's results-based inclusion rule and is a stronger positive-tissue result than soleus.

### Kidney

The corrected kidney model has a narrower higher structural and growth-factor-response pattern. ECM proteoglycans, aggregate WNT signaling, and IGF transport are higher across all three trainings; WNT is higher in all six projects, while the other two are higher in five of six. Eighty-one to 86% of measured member genes share these directions, and 63% to 80% of absolute decoder weight predicts the observed gene directions across seeds. This is coherent with prior renal ECM, fibrosis-related, nephron-remodeling, and WNT observations and adds a joint pathway-level repair or maladaptive-growth hypothesis.

The kidney claim is weaker than the spleen claim. The three programs are attenuated by composition-proxy adjustment, and their preranked-GSEA FDR values are 0.156, 1.000, and 0.433. Lower amino-acid metabolism is rejected: only 27% of measured genes share the latent direction, decoder-weighted concordance is 42%, independent ssGSEA and GSEA point higher, and prior kidney GSEA also reported higher amino-acid metabolism. Biological oxidations is retained only as heterogeneous CYP/UGT/ACSM regulation, not a broad decrease. Off-tissue neuronal, keratinization, broad GPCR, platelet, and drug-label pathways are not promoted.

Kidney is suitable for a concise secondary or exploratory manuscript result, but not as an equal-strength statistically confirmed claim.

## Consequence for manuscript scope

For a manuscript explicitly focused on tissues with interpretable results, the defensible positive core is thymus, skin, liver, and spleen. Kidney can be included as a clearly labeled secondary result because its structural-growth axis is coherent but lacks the conventional FDR support of the spleen core. Soleus does not meet the advancement criterion and should remain only in the supplementary screening and sensitivity record rather than the biological narrative. Lung and retina remain the defensible negative screens.

The final tissue-selection rule should be stated independently of biological familiarity. A defensible rule is to advance a tissue when at least one manually reviewed pathway is either triangulated across all five checks or internally robust with explicitly incomplete conventional support, provided the member-gene and tissue-context review does not contradict the latent interpretation. Literature should interpret pathways after this rule, not determine which tissues enter the manuscript.
