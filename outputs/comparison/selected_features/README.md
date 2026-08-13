# Selected gene and pathway comparison

This directory is the comparison-ready result bundle for a reader who wants to
match another model's genes or pathways against this project. Start with
[`selected_feature_comparison.xlsx`](selected_feature_comparison.xlsx). The
workbook contains the same records as the TSV files, split into labeled sheets.

## Start here

| File | Rows | Use |
|---|---:|---|
| [`gene_crosswalk.tsv`](gene_crosswalk.tsv) | 2,048 | One combined index of expiMap pathway members, all-arm and selected-arm stable generative features, 21 primary matched genes, and 49 secondary consensus genes |
| [`pathway_crosswalk.tsv`](pathway_crosswalk.tsv) | 26 | The 16 retained expiMap pathways and ten grouped generative Reactome results in common columns |
| [`selected_feature_comparison.xlsx`](selected_feature_comparison.xlsx) | 12 sheets | Spreadsheet version of the complete bundle |

Join gene lists on `gene_id` when possible. It is the versionless mouse Ensembl
identifier used by the models. `gene_symbol` is supplied for readability and
for comparisons with symbol-based lists. Join pathways on `pathway_id`, which is
the canonical Reactome `R-MMU-<number>` identifier. `pathway_term` preserves the
longer annotation key used by the model.

## expiMap files

| File | Rows | Contents |
|---|---:|---|
| [`expimap_retained_pathways.tsv`](expimap_retained_pathways.tsv) | 16 | Every final retained pathway, including tissue, direction, three-seed effects, project agreement, GSEA FDR, evidence role, and robustness status |
| [`expimap_retained_pathway_gene_summary.tsv`](expimap_retained_pathway_gene_summary.tsv) | 743 | One row per tissue and gene represented in a retained pathway |
| [`expimap_retained_pathway_members.tsv.gz`](expimap_retained_pathway_members.tsv.gz) | 813 | Every retained pathway-to-gene membership with the observed gene effect, gene FDR, pathway-direction agreement, and decoder weight when available |

The expiMap gene files require careful wording. expiMap selected pathways, not a
standalone gene panel. A gene in these files is a measured member of a retained
pathway. `independently_selected_gene` is therefore always false. Use
`member_support_class`, `n_concordant_pathways`, `any_concordant_bh_fdr`, and
`minimum_gene_fdr` to identify member genes that support a pathway's direction.
`retained_pathway_ids` contains canonical accessions;
`retained_pathway_terms` preserves the labeled GMT terms.

## Generative feature files

| File | Rows | Contents |
|---|---:|---|
| [`generative_analysis_coverage.tsv`](generative_analysis_coverage.tsv) | 27 | Selected classifier arm and feature-result counts for every tissue or muscle analysis unit |
| [`generative_all_arm_stable_features.tsv.gz`](generative_all_arm_stable_features.tsv.gz) | 3,262 | Arm-specific stable rows for every tested synthetic arm across all 27 units, representing 1,307 unique tissue-gene pairs |
| [`generative_selected_arm_stable_features.tsv`](generative_selected_arm_stable_features.tsv) | 679 | Narrower stable union from the real-only classifier and selected synthetic-supported arm in the 22 units where a synthetic arm was retained |
| [`generative_full_selected_feature_comparison.tsv.gz`](generative_full_selected_feature_comparison.tsv.gz) | 4,475 | Complete selected-arm feature union before the stability filter for those 22 units |
| [`generative_matched_genes.tsv`](generative_matched_genes.tsv) | 21 | Primary genes with real-data BH-FDR association and synthetic-supported marginal importance |
| [`generative_consensus_genes.tsv`](generative_consensus_genes.tsv) | 49 | Secondary promoted or reinforced consensus genes |
| [`generative_grouped_pathways.tsv`](generative_grouped_pathways.tsv) | 10 | Reactome groups retained by grouped permutation importance and SHAP |

The feature-importance run completed all 27 units. The utility workflow retained
the real-only arm for cecum, colon, liver, EDL, and quadriceps, so those five do
not have a selected synthetic-arm comparison in the selected-arm full or stable
tables. They are present in the all-arm stable table. The coverage file makes
this distinction explicit. Candidate tables can still contain a result from one
of these units because their matched and consensus screens evaluated qualifying
arm-level evidence separately.

In either stable table, a gene is stable in the real-only classifier, a
synthetic-supported classifier, or both. Stability required selection in at
least 50% of repeated fits and at least 75% agreement in coefficient sign. The
thresholds and observed frequencies are included in each row. In
`gene_crosswalk.tsv`, the `generative_any_arm_stable_feature` and
`generative_selected_arm_stable_feature` flags distinguish the complete and
narrower definitions.

For the stable and full feature tables, `real_permutation_roc_auc` is the score
loss after shuffling a gene in the real-only classifier. The
`arm_real_permutation_roc_auc` column is the same test for the named synthetic
arm, evaluated on held-out real profiles. Larger positive values indicate
greater predictive importance. The selected-arm tables also include SHAP
direction from the fitted classifier. SHAP values are not
differential-expression effects. Real FLT/GC effect and FDR fields are populated
only for genes that also pass the real-data association screen.

## Rebuild

The GENCODE symbol map and all paper-facing candidate tables are tracked. The
full feature-importance inputs are local analysis outputs. Recreate them first
if they are absent:

```bash
python -m nasa_mouse_diffusion.paper_parity.classifier_importance \
  --config configs/generative/diffusion/classifier_importance_osdr_disjoint.yaml
python -m nasa_mouse_internship_report.build_comparison_exports
```

[`manifest.json`](manifest.json) records the source and output checksums. Effect
magnitudes should be compared within an analysis, not between expiMap latent
scores, expression effects, permutation losses, and SHAP values.
