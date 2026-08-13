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

This page is the result inventory for model-to-model comparison. The broader
[`results_guide.md`](../../../docs/results_guide.md) connects headline findings
to figures, manuscripts, source tables, and biological interpretation across
the whole project.

Join gene lists on `gene_id` when possible. It is the versionless mouse Ensembl
identifier used by the models. `gene_symbol` is supplied for readability and
for comparisons with symbol-based lists. Join pathways on `pathway_id`, which is
the canonical Reactome `R-MMU-<number>` identifier. `pathway_term` preserves the
longer annotation key used by the model.

## How to read the result labels

- `FLT higher` and `FLT lower` are the observed flight-minus-ground direction.
- In expiMap, `retained` means the pathway entered the final paper-facing set
  after directional robustness filtering and pathway-family, member-gene,
  redundancy, and tissue-fit review. Every retained pathway passed the
  held-out-project, all-three-training-seed, and composition-proxy direction
  checks. Fourteen also agreed with both ssGSEA and preranked GSEA and are
  classified as `triangulated`. Skin cell-junction organization and thymus
  lymphoid-stromal interaction passed the internal checks but had incomplete
  conventional-method support.
- Retained does not mean significant by GSEA FDR, unanimous in every accession
  or protocol subgroup, or confirmed by prior literature. GSEA FDR and project
  agreement remain separate columns. Kidney pathways are marked `secondary
  exploratory` despite passing the five directional checks because conventional
  FDR support was weak and composition adjustment attenuated their effects.
- A primary matched generative gene has a real-data BH-FDR association and
  synthetic-supported marginal permutation or SHAP importance.
- A consensus `promoted` gene became stable only with an eligible
  synthetic-informed ranking arm. A `reinforced` gene was stable with both real
  and synthetic-informed ranking. These labels describe feature selection, not
  biological novelty.
- The expiMap literature types are `aligned`, `complementary`, and `context
  sensitive`. The generative types are `aligning`, `complementary`, `ambiguous`,
  and `unmatched`. Aligned or aligning means prior evidence supports a compatible
  feature, tissue, and direction. Complementary means prior work supports a
  related process rather than an exact replication. Context sensitive or
  ambiguous records have material variation or mixed published context.
  Unmatched means the targeted search found no specific prior match.

The literature labels organize follow-up. They are not additional statistical
evidence. The cited rationale for each row is in the linked detailed tables and
the [annotation record](../../../docs/annotation_provenance.md).

The tables below show a literature type for every retained expiMap pathway,
primary matched gene, grouped generative pathway, and secondary consensus gene.
The 743 expiMap pathway-member genes and 1,307 stable generative tissue-gene
pairs were not each reviewed as paper-facing biological candidates, so they do
not have individual literature types.

## expiMap retained pathways

These are all 16 retained pathway records. expiMap selected the pathway score;
the member genes listed later are supporting measurements rather than a
separately selected gene panel.

| Tissue | Role | Direction | Pathway | Reactome | Literature relation |
|---|---|---|---|---|---|
| Liver | main | FLT lower | MHC class II antigen presentation | [R-MMU-2132295](https://reactome.org/PathwayBrowser/#/R-MMU-2132295) | complementary |
| Liver | main | FLT lower | T-cell receptor signaling | [R-MMU-202403](https://reactome.org/PathwayBrowser/#/R-MMU-202403) | complementary |
| Skin | main | FLT lower | Chromatin-modifying enzymes | [R-MMU-3247509](https://reactome.org/PathwayBrowser/#/R-MMU-3247509) | complementary |
| Skin | main | FLT lower | Hedgehog signaling | [R-MMU-5358351](https://reactome.org/PathwayBrowser/#/R-MMU-5358351) | complementary |
| Skin | main | FLT lower | DNA repair | [R-MMU-73894](https://reactome.org/PathwayBrowser/#/R-MMU-73894) | context sensitive |
| Skin | main | FLT lower | Sphingolipid metabolism | [R-MMU-428157](https://reactome.org/PathwayBrowser/#/R-MMU-428157) | complementary |
| Skin | main | FLT lower | Cell-cell junction organization | [R-MMU-421270](https://reactome.org/PathwayBrowser/#/R-MMU-421270) | complementary |
| Spleen | main | FLT lower | T-cell receptor signaling | [R-MMU-202403](https://reactome.org/PathwayBrowser/#/R-MMU-202403) | aligned |
| Spleen | main | FLT lower | Neutrophil degranulation program | [R-MMU-6798695](https://reactome.org/PathwayBrowser/#/R-MMU-6798695) | complementary |
| Spleen | main | FLT lower | C-type lectin receptor signaling | [R-MMU-5621481](https://reactome.org/PathwayBrowser/#/R-MMU-5621481) | complementary |
| Thymus | main | FLT lower | DNA repair | [R-MMU-73894](https://reactome.org/PathwayBrowser/#/R-MMU-73894) | complementary |
| Thymus | main | FLT lower | Lymphoid-stromal interactions | [R-MMU-198933](https://reactome.org/PathwayBrowser/#/R-MMU-198933) | complementary |
| Thymus | main | FLT lower | RHOA cytoskeletal cycle | [R-MMU-8980692](https://reactome.org/PathwayBrowser/#/R-MMU-8980692) | complementary |
| Kidney | secondary exploratory | FLT higher | WNT signaling | [R-MMU-195721](https://reactome.org/PathwayBrowser/#/R-MMU-195721) | complementary |
| Kidney | secondary exploratory | FLT higher | IGF transport and uptake | [R-MMU-381426](https://reactome.org/PathwayBrowser/#/R-MMU-381426) | complementary |
| Kidney | secondary exploratory | FLT higher | ECM proteoglycans | [R-MMU-3000178](https://reactome.org/PathwayBrowser/#/R-MMU-3000178) | aligned |

The complete pathway statistics and rationales are in
[`expimap_retained_pathways.tsv`](expimap_retained_pathways.tsv).

## Primary generative genes

These 21 genes combine a real OSDR BH-FDR association with synthetic-supported
marginal feature importance. Parentheses give the literature relationship.

| Tissue | Direction | Genes |
|---|---|---|
| Liver | FLT lower | `Grb10` (complementary), `Ppic` (complementary), `H2-DMa` (complementary), `Gtf2a2` (aligning) |
| Skin | FLT higher | `Plscr1` (complementary) |
| Spleen | FLT higher | `Loxl1` (complementary) |
| Thymus | FLT lower | `Nusap1` (aligning), `Stmn1` (aligning), `Birc5` (ambiguous), `Ccnb2` (aligning), `E2f2` (aligning), `Ube2c` (aligning), `Cdc20` (aligning), `Gmnn` (aligning), `Kif20a` (aligning) |
| Thymus | FLT higher | `Klhdc2` (unmatched), `Snx7` (unmatched), `Etv1` (complementary), `Plscr1` (complementary), `Tspan3` (complementary), `Socs2` (complementary) |

[`generative_matched_genes.tsv`](generative_matched_genes.tsv) gives Ensembl
IDs, effect sizes, FDR values, feature-importance patterns, annotation rationale,
and literature source IDs for every gene.

## Grouped generative pathways

These ten Reactome groups passed the grouped permutation and SHAP workflow.

| Tissue | Direction | Pathway | Reactome | Literature relation |
|---|---|---|---|---|
| Skin | FLT higher | RIPK1-mediated regulated necrosis | [R-MMU-5213460](https://reactome.org/PathwayBrowser/#/R-MMU-5213460) | complementary |
| Skin | FLT higher | Regulation of necroptotic cell death | [R-MMU-5675482](https://reactome.org/PathwayBrowser/#/R-MMU-5675482) | complementary |
| Spleen | FLT higher | Activation of the AP-1 family of transcription factors | [R-MMU-450341](https://reactome.org/PathwayBrowser/#/R-MMU-450341) | ambiguous |
| Thymus | FLT lower | Regulation of APC/C activators between G1/S and early anaphase | [R-MMU-176408](https://reactome.org/PathwayBrowser/#/R-MMU-176408) | aligning |
| Thymus | FLT lower | G2/M DNA replication checkpoint | [R-MMU-69478](https://reactome.org/PathwayBrowser/#/R-MMU-69478) | aligning |
| Thymus | FLT lower | Condensation of prophase chromosomes | [R-MMU-2299718](https://reactome.org/PathwayBrowser/#/R-MMU-2299718) | aligning |
| Thymus | FLT lower | Cdc20:Phospho-APC/C-mediated degradation of Cyclin A | [R-MMU-174184](https://reactome.org/PathwayBrowser/#/R-MMU-174184) | aligning |
| Thymus | FLT lower | APC:Cdc20-mediated degradation of cell-cycle proteins before checkpoint satisfaction | [R-MMU-179419](https://reactome.org/PathwayBrowser/#/R-MMU-179419) | aligning |
| Thymus | FLT higher | ERBB2 activates PTK6 signaling | [R-MMU-8847993](https://reactome.org/PathwayBrowser/#/R-MMU-8847993) | ambiguous |
| Thymus | FLT lower | APC/C:Cdc20-mediated degradation of Securin | [R-MMU-174154](https://reactome.org/PathwayBrowser/#/R-MMU-174154) | aligning |

Member genes, grouped importance, real-data effects, FDR values, and annotation
rationales are in
[`generative_grouped_pathways.tsv`](generative_grouped_pathways.tsv).

## Secondary consensus genes

The consensus analysis produced 49 tissue-gene records: 26 promoted and 23
reinforced. This table lists all of them. The complete statistics and literature
annotations are in
[`generative_consensus_genes.tsv`](generative_consensus_genes.tsv).

| Scope | Tissue | Status | Direction | Genes |
|---|---|---|---|---|
| Tissue | Adrenal gland | promoted | FLT lower | `Psmb8` (unmatched) |
| Tissue | Adrenal gland | reinforced | FLT lower | `Tspan4` (unmatched) |
| Tissue | Eye | reinforced | FLT lower | `Klhl21` (aligning) |
| Tissue | Kidney | promoted | FLT higher | `Inpp4b` (complementary) |
| Tissue | Kidney | reinforced | FLT higher | `Slc37a4` (complementary) |
| Tissue | Skeletal muscle | promoted | FLT lower | `Klhl21` (complementary), `Mapkapk5` (complementary), `Reep5` (complementary), `Itgb5` (complementary) |
| Tissue | Skeletal muscle | reinforced | FLT higher | `Sox4` (aligning), `Cebpd` (aligning), `Sh3bp5` (complementary), `Prkcd` (aligning), `Arid5b` (aligning), `Sesn1` (aligning), `Tle1` (complementary) |
| Tissue | Skeletal muscle | reinforced | FLT lower | `Bphl` (unmatched) |
| Tissue | Skin | promoted | FLT higher | `Plscr1` (complementary) |
| Tissue | Spleen | promoted | FLT higher | `Rai14` (complementary), `Myl9` (complementary), `Ptprk` (complementary) |
| Tissue | Spleen | reinforced | FLT higher | `Loxl1` (complementary) |
| Tissue | Thymus | promoted | FLT higher | `Hsd17b11` (complementary), `Etv1` (complementary) |
| Tissue | Thymus | promoted | FLT lower | `Nusap1` (aligning), `Stmn1` (aligning), `Birc5` (ambiguous), `Cdk1` (aligning), `Top2a` (aligning), `Ccnb2` (aligning), `Aurka` (aligning), `Ccne2` (aligning), `Kif20a` (aligning), `Pcna` (aligning), `Ccnf` (aligning) |
| Tissue | Thymus | reinforced | FLT higher | `Snx7` (unmatched) |
| Tissue | Thymus | reinforced | FLT lower | `Ube2c` (aligning), `Gmnn` (aligning) |
| Muscle group | Gastrocnemius | promoted | FLT higher | `Nfkbia` (aligning) |
| Muscle group | Gastrocnemius | promoted | FLT lower | `Fhl2` (complementary) |
| Muscle group | Soleus | reinforced | FLT higher | `Tpm1` (ambiguous) |
| Muscle group | Soleus | reinforced | FLT lower | `Bdh1` (aligning), `Ech1` (aligning), `Bnip3` (ambiguous), `Decr1` (aligning) |
| Muscle group | Tibialis anterior | promoted | FLT higher | `Cebpd` (complementary) |
| Muscle group | Tibialis anterior | reinforced | FLT higher | `Cdkn1a` (ambiguous), `St3gal5` (complementary), `Bnip3` (complementary) |

## Stable feature coverage

The full all-arm table has 3,262 arm-specific stable rows, representing 1,307
unique tissue-gene pairs. The chosen-arm table has 679 rows from the 22 units
where the utility workflow retained a synthetic-supported arm. Counts below are
not claims of biological association; they show the available feature lists.

| Scope | Tissue | All-arm stable genes | Chosen-arm stable genes | Primary genes | Consensus genes | Grouped pathways |
|---|---|---:|---:|---:|---:|---:|
| Tissue | Adrenal gland | 25 | 18 | 0 | 2 | 0 |
| Tissue | Bone | 36 | 14 | 0 | 0 | 0 |
| Tissue | Bone marrow | 22 | 19 | 0 | 0 | 0 |
| Tissue | Brain | 48 | 41 | 0 | 0 | 0 |
| Tissue | Brown adipose tissue | 28 | 19 | 0 | 0 | 0 |
| Tissue | Cecum | 32 | 0 | 0 | 0 | 0 |
| Tissue | Cerebellum | 38 | 29 | 0 | 0 | 0 |
| Tissue | Colon | 59 | 0 | 0 | 0 | 0 |
| Tissue | Eye | 26 | 24 | 0 | 1 | 0 |
| Tissue | Heart | 28 | 21 | 0 | 0 | 0 |
| Tissue | Hippocampus | 28 | 15 | 0 | 0 | 0 |
| Tissue | Kidney | 90 | 43 | 0 | 2 | 0 |
| Tissue | Liver | 92 | 0 | 4 | 0 | 0 |
| Tissue | Lung | 69 | 45 | 0 | 0 | 0 |
| Tissue | Mammary gland | 21 | 15 | 0 | 0 | 0 |
| Tissue | Optic nerve | 26 | 14 | 0 | 0 | 0 |
| Tissue | Retina | 70 | 50 | 0 | 0 | 0 |
| Tissue | Skeletal muscle | 136 | 98 | 0 | 12 | 0 |
| Tissue | Skin | 91 | 46 | 1 | 1 | 2 |
| Tissue | Spleen | 70 | 50 | 1 | 4 | 1 |
| Tissue | Thymus | 86 | 45 | 15 | 16 | 7 |
| Tissue | White adipose tissue | 25 | 16 | 0 | 0 | 0 |
| Muscle group | EDL | 24 | 0 | 0 | 0 | 0 |
| Muscle group | Gastrocnemius | 28 | 21 | 0 | 2 | 0 |
| Muscle group | Quadriceps | 63 | 0 | 0 | 0 | 0 |
| Muscle group | Soleus | 29 | 22 | 0 | 5 | 0 |
| Muscle group | Tibialis anterior | 17 | 14 | 0 | 4 | 0 |

Use [`gene_crosswalk.tsv`](gene_crosswalk.tsv) for the unique gene-level union,
[`generative_all_arm_stable_features.tsv.gz`](generative_all_arm_stable_features.tsv.gz)
for every arm-specific row, and
[`generative_selected_arm_stable_features.tsv`](generative_selected_arm_stable_features.tsv)
for the narrower chosen-arm result.

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
