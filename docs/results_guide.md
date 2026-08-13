# Results guide

This page is the quickest way to find the final numbers, plots, and biological
interpretations. The frozen tables in each paper's `source_data/` directory are
the numeric record. The figures visualize those tables, and the manuscripts
explain how the results were obtained and what can be concluded from them.

For a continuous account of the internship, read the
[internship report](../paper/slstp_internship_report/manuscript.pdf). For a
shorter visual account, use the
[final presentation](../presentation/final/SLSTP_2026_Generative_Transcriptomics.pdf).

The standalone final plots are grouped in three places:

- [combined internship-report figures](../paper/slstp_internship_report/figures/);
- [expiMap main and supplementary figures](../paper/asgsr_expimap_hvg/figures/);
- [generative main and supplementary figures](../paper/synthetic_guided_spaceflight/figures/).

## Gene and pathway comparison files

For comparison with another model, start with the
[selected-feature workbook](../outputs/comparison/selected_features/selected_feature_comparison.xlsx)
or the [bundle README](../outputs/comparison/selected_features/README.md). The
workbook has separate sheets for expiMap pathways, expiMap pathway-member genes,
generative feature importance, primary matched genes, secondary consensus
genes, grouped Reactome pathways, and coverage of all 27 generative analysis
units.

Two TSV files provide direct cross-method join points:

- [`gene_crosswalk.tsv`](../outputs/comparison/selected_features/gene_crosswalk.tsv)
  uses versionless mouse Ensembl IDs and gene symbols across both methods;
- [`pathway_crosswalk.tsv`](../outputs/comparison/selected_features/pathway_crosswalk.tsv)
  uses Reactome `R-MMU-*` pathway IDs.

The expiMap model selected pathways. Its gene rows are measured members of
those retained pathways, not an independently selected gene panel. The
generative bundle separately identifies 1,307 unique stable tissue-gene pairs
across all tested arms, the 679 selected-arm stable rows, the 21 primary matched
genes, and the 49 secondary consensus genes.
[`generative_analysis_coverage.tsv`](../outputs/comparison/selected_features/generative_analysis_coverage.tsv)
shows which units retained a synthetic-supported arm and the result count from
each feature analysis.

This guide is the narrative index for the whole project: findings, figures,
source tables, and interpretation. The comparison bundle is narrower. Its README
lists the selected genes and pathways themselves and its workbook provides
filterable cross-method tables.

## GLARE and MOBER

GLARE is an informative negative result in the final report. In pooled skeletal
muscle, MOBER reduced the accession silhouette from 0.272 to 0.075, but the
FLT/GC silhouette remained near zero, changing from 0.008 to -0.003. Study
structure weakened without producing a clear condition separation.

| Item | Location |
|---|---|
| Exact tissue-level silhouette values | [Aggregate versus MOBER table](../outputs/glare/study_effects/aggregate_vs_mober_study_effect_summary.tsv) |
| Main plot | [GLARE batch-effect figure](../paper/slstp_internship_report/figures/figure_2_glare_batch_effects.png) |
| Interpretation | [Internship report, Section 3](../paper/slstp_internship_report/manuscript.md#3-glare-and-batch-correction) |
| Exploratory module records | [Curated GLARE outputs](../outputs/glare/README.md) |

The exploratory GLARE module tables are preserved for audit and follow-up.
They are not part of the final biological claim set.

## expiMap pathway results

The final expiMap set contains 16 pathway records from four primary tissues and
one exploratory tissue. Exact effects for all three training seeds, project
directions, GSEA values, composition sensitivity, and evidence labels are in
[Table 2](../paper/asgsr_expimap_hvg/source_data/table_2_retained_pathway_evidence.tsv).

| Tissue | Reported pattern | Study coverage |
|---|---|---|
| Thymus | Lower DNA repair, RHOA cytoskeletal cycle, and lymphoid-stromal interaction scores | 5 projects |
| Skin | Lower chromatin regulation, DNA repair, Hedgehog, sphingolipid, and cell-junction scores | 4 project summaries, with protocol-dependent exceptions |
| Liver | Lower MHC class II antigen presentation and T-cell receptor signaling | Lower in 8 of 9 project summaries |
| Spleen | Lower T-cell receptor, neutrophil degranulation, and C-type lectin receptor scores | Lower in all 5 unconfounded projects; all three GSEA FDR values are below 0.05 |
| Kidney | Higher ECM proteoglycan, WNT, and IGF transport scores | Exploratory; conventional GSEA FDR values exceed 0.05 and composition adjustment attenuates the effects |

| Question | Plot | Exact data or explanation |
|---|---|---|
| Which pathways changed, and in which direction? | [Figure 3: pathway shifts](../paper/asgsr_expimap_hvg/figures/figure_3_tissue_pathway_shifts.png) | [Table 2: retained pathways](../paper/asgsr_expimap_hvg/source_data/table_2_retained_pathway_evidence.tsv) |
| Did the direction survive the supporting checks? | [Figure 4: evidence and member genes](../paper/asgsr_expimap_hvg/figures/figure_4_evidence_gene_support.png) | [Table S24: pathway evidence](../paper/asgsr_expimap_hvg/source_data/table_s24_pathway_robustness_evidence.tsv) and [Table S34: member-gene support](../paper/asgsr_expimap_hvg/source_data/table_s34_retained_pathway_member_gene_support.tsv) |
| How did skin protocol context affect the result? | [Figure 5: skin protocol contrasts](../paper/asgsr_expimap_hvg/figures/figure_5_skin_protocol_context.png) | [Table S8: protocol effects](../paper/asgsr_expimap_hvg/source_data/table_s8_skin_protocol_context_effects.tsv) |
| What biological hypotheses follow from the pathway results? | [Figure 6: tissue hypotheses](../paper/asgsr_expimap_hvg/figures/figure_6_tissue_state_hypotheses.png) | [Detailed expiMap manuscript](../paper/asgsr_expimap_hvg/manuscript.pdf) |
| What did each accession contribute? | [Internship-report heatmap](../paper/slstp_internship_report/figures/figure_4_expimap_pathway_heatmap.png) | [Table S5: accession effects](../paper/asgsr_expimap_hvg/source_data/table_s5_accession_pathway_effects.tsv.gz) |

The [expiMap package README](../paper/asgsr_expimap_hvg/README.md) lists the
complete model scope and explains why kidney is exploratory and soleus was not
advanced into the paper-facing result set.

## Generator validation

The generative inventory contains 1,610 OSDR profiles from 75 accessions: 835
flight and 775 ground-control profiles. The cohort total and model splits are in
[Table 1](../paper/synthetic_guided_spaceflight/source_data/table_1_data_inventory.tsv).
The [OSDR tissue inventory](../outputs/generative/benchmark/data_audit/osdr/osdr_tissue_inventory.tsv)
gives FLT, GC, accession, and training counts for every canonical tissue.

The OSDR-adapted DDIM was chosen for downstream work. It had near-chance
real-versus-synthetic adversarial accuracy and the lowest reported
Frechet-to-real-neighborhood ratio while preserving high correlation,
precision, recall, and F1.

| Model | Correlation | Precision | Recall | F1 | Adversarial accuracy | FD / real P95 |
|---|---:|---:|---:|---:|---:|---:|
| ARCHS4 DDIM initialization | 0.878 | 0.951 | 0.890 | 0.919 | 0.515 | 0.866 |
| Study-conditioned WGAN-GP | 0.976 | 0.976 | 0.994 | 0.985 | 0.636 | 0.144 |
| OSDR-adapted factorized DDIM | 0.974 | 0.997 | 0.996 | 0.997 | 0.475 | 0.074 |

Adversarial accuracy closer to 0.5 means that an external classifier had more
difficulty separating real from generated profiles. Lower FD / real P95 is
better. The models used their stated evaluation splits, so this table supports
model choice within this project rather than a universal model ranking.

| Item | Location |
|---|---|
| Exact model metrics and decisions | [Table 4](../paper/synthetic_guided_spaceflight/source_data/table_4_generator_model_selection.tsv) |
| Metric plot | [Figure 1](../paper/synthetic_guided_spaceflight/figures/figure_1_generator_validation.png) |
| ARCHS4 denoising trajectory | [Figure 2A](../paper/synthetic_guided_spaceflight/figures/figure_2a_archs4_denoising_trajectory.png) |
| Real and generated OSDR profiles | [Figure 2B](../paper/synthetic_guided_spaceflight/figures/figure_2b_locked_real_vs_synthetic_pca.png) |
| Full validation explanation | [Generative manuscript](../paper/synthetic_guided_spaceflight/manuscript.pdf) |

## Synthetic-data utility

The matched benchmark fitted 648 classifiers across 27 tissues and anatomical
muscle groups. Real-plus-synthetic training was no worse than real-only training
on all six pooled and accession-macro metrics in 18 units. Sixteen of those
improved at least one metric. The table below shows the eight units emphasized
in the manuscript.

| Tissue | Balanced-accuracy change | AUROC change | Average-precision change |
|---|---:|---:|---:|
| Eye | +0.094 | +0.156 | +0.115 |
| Retina | +0.100 | +0.075 | +0.063 |
| Lung | +0.062 | +0.086 | +0.108 |
| Skin | +0.074 | +0.087 | +0.076 |
| Thymus | +0.061 | +0.046 | +0.042 |
| Spleen | +0.057 | +0.075 | +0.072 |
| Liver | +0.043 | +0.037 | +0.029 |
| Skeletal muscle, pooled | +0.010 | +0.011 | +0.014 |

These are changes on held-out real profiles for real-plus-synthetic training
relative to real-only training. The complete 27-unit result, including negative
and tied cases, is in
[Table S18](../paper/synthetic_guided_spaceflight/source_data/table_s18_matched_all_gene_utility.tsv).
[Figure S2](../paper/synthetic_guided_spaceflight/figures/figure_s2_downstream_utility.png)
shows the full downstream utility summary.

## Synthetic-informed biology

The primary matched analysis identified 21 tissue-gene associations that had a
real-data BH-FDR association and synthetic-supported classifier importance.
Thymus contributed 15, liver four, skin one, and spleen one. The secondary
consensus analysis retained 49 associations: 26 promoted and 23 reinforced.
All association statistics came from real OSDR profiles.

| Tissue | Main result | Plot | Exact records |
|---|---|---|---|
| Thymus | 15 matched genes, including a flight-lower mitotic program; seven genes gained measurable marginal importance after adding synthetic training profiles | [Figure 3](../paper/synthetic_guided_spaceflight/figures/figure_3_thymus_biology.png) | [Matched genes](../paper/synthetic_guided_spaceflight/source_data/table_s19_matched_all_gene_candidates.tsv), [matched Reactome enrichment](../paper/synthetic_guided_spaceflight/source_data/table_s20_matched_candidate_reactome.tsv), and [secondary consensus panel](../paper/synthetic_guided_spaceflight/source_data/table_s4_thymus_core_genes.tsv) |
| Liver | Four flight-lower matched genes: *Grb10*, *Ppic*, *H2-DMa*, and *Gtf2a2*; no shared Reactome enrichment | [Figure 5](../paper/synthetic_guided_spaceflight/figures/figure_5_tissue_evidence.png) | [Matched genes](../paper/synthetic_guided_spaceflight/source_data/table_s19_matched_all_gene_candidates.tsv) |
| Skin | Flight-higher *Plscr1* in both matched and consensus analyses; grouped importance also identified two complementary necroptosis-related terms | [Figure 5](../paper/synthetic_guided_spaceflight/figures/figure_5_tissue_evidence.png) | [Matched genes](../paper/synthetic_guided_spaceflight/source_data/table_s19_matched_all_gene_candidates.tsv) and [grouped pathways](../paper/synthetic_guided_spaceflight/source_data/table_s23_grouped_pathway_literature_annotations.tsv) |
| Spleen | Flight-higher *Loxl1* passed the matched analysis; the consensus panel added *Rai14*, *Ptprk*, and *Myl9*, suggesting a tentative adhesion and tissue-mechanics response | [Figure 5](../paper/synthetic_guided_spaceflight/figures/figure_5_tissue_evidence.png) | [Matched genes](../paper/synthetic_guided_spaceflight/source_data/table_s19_matched_all_gene_candidates.tsv) and [consensus associations](../paper/synthetic_guided_spaceflight/source_data/table_s10_synthetic_informed_bh_fdr_genes.tsv) |
| Soleus | Secondary consensus reinforcement of lower *Bdh1*, *Ech1*, *Bnip3*, and *Decr1* with higher *Tpm1*, connecting oxidative metabolism, mitochondrial turnover, and contractile remodeling | [Figure 4](../paper/synthetic_guided_spaceflight/figures/figure_4_soleus_biology.png) | [Soleus genes](../paper/synthetic_guided_spaceflight/source_data/table_s7_soleus_genes.tsv) and [Reactome terms](../paper/synthetic_guided_spaceflight/source_data/table_s8_muscle_reactome.tsv) |

[Table 6](../paper/synthetic_guided_spaceflight/source_data/table_6_tissue_evidence.tsv)
is the compact tissue-by-tissue interpretation. The full real-data association
screen is in
[Table S11](../paper/synthetic_guided_spaceflight/source_data/table_s11_all_random_effects_bh_fdr_genes.tsv),
including genes that did not receive synthetic feature-selection support.

## Literature annotations

The annotation tables retain both the label and the reason for it:

| Analysis | Annotation table |
|---|---|
| expiMap pathway interpretation | [Tissue-specific records and sources](../paper/asgsr_expimap_hvg/source_data/literature_review/README.md) |
| Matched individual genes | [Table S22](../paper/synthetic_guided_spaceflight/source_data/table_s22_matched_gene_literature_annotations.tsv) |
| Grouped Reactome pathways | [Table S23](../paper/synthetic_guided_spaceflight/source_data/table_s23_grouped_pathway_literature_annotations.tsv) |
| Secondary consensus genes | [Table S16](../paper/synthetic_guided_spaceflight/source_data/table_s16_promoted_gene_literature_annotations.tsv) |
| Full generative source catalog | [Table S24](../paper/synthetic_guided_spaceflight/source_data/table_s24_importance_literature_sources.tsv) |

The [annotation provenance guide](annotation_provenance.md) explains every
field and rebuild path. The [prompt record](annotation_prompts.md) contains the
review instructions and label definitions.

## Which files are authoritative

Use paper `source_data/` tables for final numeric values and paper figures for
the corresponding plots. Use `outputs/` to trace a result back to its selected
run or analysis. Development runs and exploratory GLARE tables remain useful
for audit, but they do not replace the frozen paper records.
