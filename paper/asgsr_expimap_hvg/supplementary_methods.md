# Supplementary Methods and Reproducibility Notes

## Analysis scope

The analysis screened six retrospectively reassessed reference-query, 2,000-HVG models. Thymus, skin, liver, and spleen form the positive-result set, and kidney is secondary exploratory evidence. Soleus did not meet the advancement criteria and is retained only as a supplementary screening and sensitivity record. The broader repository screen also included lung, retina, and aggregate skeletal muscle, and the tissue set was not prespecified. Each model contains only constrained Reactome latent dimensions. The word *complementary* describes a Reactome pathway interpretation that adds tissue context beyond the principal prior-literature phenotype. It does not denote an unconstrained de novo latent program. De novo variants were evaluated separately during development but are outside this comparison because their leading programs were less stable and often overlapped existing pathway biology.

This package supersedes older presentation annotations that assigned biological direction directly from raw latent signs. Raw signs are arbitrary; only the decoder-oriented effects in this package should be used for higher or lower pathway statements.

The compact main review and broad screening reviews serve different purposes. The revised paper retains 16 pathways: 13 across the four main tissues and three secondary kidney pathways. The original expanded review starts with every active thymus, skin, liver, and soleus program, includes the top 10% by absolute effect, adds directionally stable programs through rank 40, and retains every original display pathway. Its union contains 153 records: 42 thymus, 36 skin, 39 liver, and 36 soleus pathways. The corrected kidney and spleen review instead evaluates every primary top-decile pathway, then applies complete-label tissue review, seed and project checks, conventional enrichment, composition sensitivity, member-gene direction, decoder-weight direction, and redundancy review.

Each of the original 153 pathways was assigned explicitly to a tissue-specific process family using the current Reactome hierarchy, pathway gene overlap, tissue context, and primary literature; no keyword classifier assigned the biological role. This produced 37 nonredundant tissue-family records. Opposite nested pairs were retained in source data but not interpreted as independent branches. The corrected kidney and spleen manual review likewise rejects tissue-incongruent neuronal, keratinization, broad GPCR, drug-label, and unstable parent terms rather than promoting them on magnitude alone.

## Exact model directories

| Tissue | Model directory |
| --- | --- |
| Thymus | `outputs/expimap_archs4_reference_osdr_query_thymus/tutorial_hvg_2000/query_nb_250epoch_seed2020` |
| Skin | `outputs/expimap_archs4_reference_osdr_query_skin/tutorial_hvg_2000/query_nb_250epoch_seed2020` |
| Liver primary, de-duplicated | `outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/query_nb_250epoch_seed2020_primary_deduplicated` |
| Liver full-input sensitivity | `outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/query_nb_250epoch_seed2020` |
| Spleen corrected | `outputs/expimap_archs4_reference_osdr_query_spleen/reassessment_hvg_2000/query_nb_250epoch_seed2020` |
| Kidney corrected, secondary | `outputs/expimap_archs4_reference_osdr_query_kidney/reassessment_hvg_2000/query_nb_250epoch_seed2020` |
| Soleus, non-advanced supplementary screen | `outputs/expimap_muscle_targeted_combined_min8/tutorial_hvg_soleus_2000/query_nb_250epoch_seed2020` |

## Direction and effect definitions

Raw variational latent signs are not intrinsically identifiable. For every pathway, scores were multiplied by the sign returned by `EXPIMAP.latent_directions(method="sum")`. This orients a positive latent direction toward the net positive decoder-weight direction. A reported positive flight-minus-ground effect therefore means a higher decoder-oriented pathway score in flight; it does not mean that every pathway member gene increased.

For accession *a* and pathway *p*, the within-accession effect was

`delta[a,p] = mean(score[p] | flight, a) - mean(score[p] | ground control, a)`.

The descriptive effect was the arithmetic mean of `delta[a,p]` over accessions, so a large accession did not dominate a small accession. For held-out and corrected-model project summaries, accessions sharing one mission project were first averaged. This collapses paired skin sites and overlapping project sources without claiming that projects are otherwise exchangeable. The all-sample mean difference, Welch test, Mann-Whitney test, Wilcoxon signed-rank test, random-effects meta-analysis, and leave-one-accession-out analysis were retained as sensitivity summaries. They were not used as hard discovery gates.

The initial liver query contained 231 samples from 12 accessions. OSD-168 was removed because its RR-1 and RR-3 cohorts are represented by OSD-48 and OSD-137 and because it includes with-ERCC and without-ERCC technical variants of the same RR-1 samples. OSD-164 was removed because its RR-1 animals overlap OSD-47. The resulting 197-sample, 10-cohort query was remapped for 250 epochs against the unchanged liver reference and is the primary liver run. The original 12-accession mapping is retained as a sensitivity analysis. In repository workflows where ERCC technical pairs must be collapsed, `prefer_noercc` retains the no-ERCC profile; the paper primary analysis instead excludes all of OSD-168 because its biological cohorts are duplicated elsewhere.

The corrected spleen model mapped all 109 samples, then excluded OSD-288 from the primary effect because its recorded strain was disjoint by condition. This leaves 100 samples and five projects. The full mapping remains available for sensitivity. The corrected kidney model uses all 135 mapped samples and six projects.

## Reviewer-directed robustness analyses

Conventional pathway benchmarks used the exact query samples, retained HVG genes, and Reactome memberships from each primary model. Rank-normalized ssGSEA used GSEApy 1.3.0, weight 0.25, and gene-set sizes of 5-500. Preranked GSEA ranked genes by the mean project-balanced `log2(CPM + 1)` flight-ground effect and used 1,000 permutations with weight 1. Project balancing first averaged accession effects that shared a mission project. In leave-one-project-out folds, effects and the top absolute decile were learned from the remaining projects only; the held-out project was used solely for directional evaluation.

Full-pipeline seed sensitivity used seeds 2020, 2021, and 2022. The two added runs retrained each 400-epoch ARCHS4 reference and each 250-epoch OSDR query with matched architecture and optimization settings. Seed-specific latent signs were oriented from the fitted decoder before calculating accession-balanced effects. Additional outputs follow the primary directory names with `seed2020` replaced by `seed2021` or `seed2022`.

Broad composition markers were derived from the Tabula Muris Senis Smart-seq2 CELLxGENE asset `be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad`. Cell types were grouped into the broad compartments listed in the manuscript, and the top 30 compartment-specific genes available in each query HVG matrix were retained. Within-accession pathway effects were adjusted for up to three principal components of the sample marker scores explaining at least 90% of marker-score variance. This is a proxy sensitivity, not cell-type deconvolution.

For Tables S24 and S27, composition support requires direction retention and an adjusted-to-unadjusted absolute effect ratio of at least 0.25; held-out support requires at least two-thirds project-direction concordance. The four robustness statuses are descriptive. They do not constitute statistical significance levels or retroactively prespecify tissue and pathway selection. GSEA FDR is therefore reported separately. This distinction is consequential for kidney: the three retained programs pass all five directional checks but have GSEA FDR of 0.156, 1.000, and 0.433.

## Covariate audit

Sex was fixed or balanced within each accession-condition comparison, but sex differed across missions and was not modeled as a single pooled effect. Age was not retained in the assembled analysis metadata and could not be adjusted. OSD-289 thymus and OSD-714 soleus had condition-disjoint strain labels. Restricted sensitivity estimates exclude these accessions. OSD-288 spleen also had condition-disjoint strain and was excluded from the corrected primary spleen effect. The broad flight label in OSD-289 includes distinct onboard gravity conditions and is another reason to keep its exclusion visible.

Skin protocol-context contrasts were derived from sample names and checked against official OSDR ISA metadata. OSD-238 and OSD-239 MHU-2 samples were separated into microgravity and onboard centrifuge-generated artificial 1 g, each against the same Earth 1 g controls. OSD-240 and OSD-241 RR-5 samples were interpreted as a 30-day flight followed by approximately 30 days of Earth recovery. Sample identifiers show that the MHU-2 dorsal and femoral accessions and the RR-5 dorsal and femoral accessions are paired anatomical samples from shared experimental cohorts. They are not duplicate expression profiles, but they are not independent mission effects. A project-balanced sensitivity therefore averaged paired-site effects within MHU-2 and RR-5 before averaging MHU-2, RR-5, RR-6, and RR-7. This sensitivity was applied to all 319 skin pathways; 291 were active latent programs. OSD-243 RR-6 samples were split into approximately 30-day live-return and approximately 60-day ISS-terminal groups. OSD-254 RR-7 samples were split by C3H/HeJ versus C57BL/6J and 25- versus 75-day terminal collection. These post hoc contrasts share controls in MHU-2 and are descriptive sensitivity analyses, not independent tests.

## Gene-level support

Unnormalized OSDR API counts were converted to `log2(counts per million + 1)` for a secondary gene-level check. For each gene, pooled flight-ground Welch tests were adjusted by Benjamini-Hochberg. Study-balanced gene effects were also computed as the mean of accession-specific expression differences. The manuscript treats these analyses as support for a latent pathway interpretation, not as proof that a pathway score is equivalent to differential expression of all member genes.

## Rebuild

From the repository root:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.build_asgsr_paper
```

The script writes all figures in PNG and vector PDF formats and rebuilds Tables S1-S5 and S8. It requires the `nasa-mouse` Conda environment and the existing trained models. The audited runs used an NVIDIA A100-SXM4-40GB GPU for training and query mapping.

The main build also invokes the expanded pathway-family review. To run that review alone after Table S9 exists:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.review_expanded_pathway_screen
```

To render the abstract and manuscript as HTML and PDF:

```bash
/home/exouser/miniforge3/envs/nasa-mouse/bin/python -m pip install \
  -r paper/asgsr_expimap_hvg/requirements.txt
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.render_asgsr_documents
```

To verify every manuscript DOI against Crossref and regenerate the source audit:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.verify_asgsr_sources
```

To rebuild the aligned and complementary program audit, including pairwise score correlations and Reactome gene overlap:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.audit_aligned_complementary_programs
```

To rebuild the conventional, held-out-project, and composition-proxy analyses:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.reviewer_robustness_analysis
```

To create missing full-pipeline seed runs and rebuild their summaries:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.run_asgsr_seed_sensitivity
```

To integrate all checks into the reviewed-pathway evidence matrix:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.integrate_reviewer_robustness
```

To rebuild the corrected kidney and spleen models, evidence audit, manual curation, and revised paper-facing figures and tables:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.run_kidney_spleen_seed_sensitivity
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.analyze_kidney_spleen_reassessment
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.curate_kidney_spleen_reassessment
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.integrate_reassessed_tissues_paper
```

The integration command must run after `build_asgsr_paper`, because it assembles the corrected main-tissue evidence and copies the original soleus material to supplementary figure names. The final publication figures and their added source tables are then generated with:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m expiMap_scarches.nasa_mouse_expimap.build_publication_figures
```

This final step authors the main figures and Figures S1, S7, and S9 at a 7.2-inch publication width, writes both 300-dpi PNG and vector PDF copies, and records dimensions and border checks in `figure_build_manifest.tsv`. It also removes superseded duplicate figures and retired generated artwork. A deterministic process synthesis is written to `presentation/expimap/asgsr_process_summary.*` for presentation use only and is not included in the scientific manuscript.

## Source-data index

- `table_s1_model_summary.tsv`: model samples, programs, training, and hardware.
- `table_s2_accession_covariate_audit.tsv`: accession sample counts and available covariates.
- `table_s3_all_pathway_effects.tsv`: all 1,427 tissue-pathway result records.
- `table_1_curated_pathway_results.tsv`: 29 literature-reviewed representative pathways, including separate evidence-role and protocol-sensitivity annotations.
- `table_s4_gene_level_results.tsv.gz`: gene-level expression sensitivity results.
- `table_s5_accession_pathway_effects.tsv.gz`: accession-specific pathway effects.
- `table_s6_program_story_audit.tsv`: evidence tier and story decision for the 23 focused aligned, complementary, or subsequently downgraded programs.
- `table_s7_program_pairwise_structure.tsv`: within-accession score correlations and Reactome gene overlap.
- `table_s8_skin_protocol_context_effects.tsv`: decoder-oriented skin effects separated by gravity, anatomical site, recovery or terminal collection, strain, and duration.
- `table_s9_systematic_pathway_screen.tsv`: all 1,427 pathways ranked within tissue, with active-program rank, accession-direction agreement, top-20 screen status, main-figure review status, and four-project skin sensitivity fields for all 319 skin pathways.
- `table_s10_expanded_pathway_review.tsv`: 153 candidate pathways with selection reason, decoder-oriented direction, tissue-specific sensitivity, gene support, nested-opposition flags, explicit family assignment, and pathway-level disposition.
- `table_s11_nonredundant_pathway_families.tsv`: 37 process-family summaries with representative term, effect range, overlap and correlation measures, literature role, narrative decision, sources, and caution.
- `table_s12_conventional_pathway_accession_effects.tsv`: accession-level expiMap and ssGSEA effects.
- `table_s13_method_benchmark.tsv`: pathway-level expiMap, ssGSEA, and preranked-GSEA comparison.
- `table_s14_method_benchmark_summary.tsv`: rank and direction agreement by tissue and pathway set.
- `table_s15_project_heldout_predictions.tsv.gz`: every leave-one-project-out pathway prediction.
- `table_s16_project_heldout_summary.tsv`: held-out direction summaries by tissue, method, and pathway set.
- `table_s17_tms_compartment_markers.tsv`: atlas-derived broad-compartment marker definitions.
- `table_s18_sample_composition_proxy_scores.tsv.gz`: sample marker scores and retained composition PCs.
- `table_s19_composition_proxy_adjusted_effects.tsv`: unadjusted and marker-PC-adjusted expiMap effects.
- `table_s20_training_seed_pathway_effects.tsv.gz`: accession-balanced pathway effects for all three full-pipeline seeds.
- `table_s21_training_seed_consensus.tsv`: pathway-level seed effect ranges and direction concordance.
- `table_s22_training_seed_summary.tsv`: seed stability by tissue and pathway set.
- `table_s23_training_seed_manifest.tsv`: inputs, outputs, epochs, timing, and GPU provenance for 12 reference-query runs.
- `table_s24_pathway_robustness_evidence.tsv`: five-check evidence matrix and descriptive robustness status for 29 reviewed pathways.
- `table_2_retained_pathway_evidence.tsv`: standardized evidence fields for the 13 main-tissue and 3 secondary kidney pathways retained in the revised manuscript.
- `table_s25_revised_model_scope.tsv`: complete six-model screening samples, roles, programs, training, and hardware provenance.
- `table_s26_kidney_spleen_training_manifest.tsv`: corrected kidney and spleen three-seed training provenance.
- `table_s27_kidney_spleen_pathway_evidence.tsv`: complete corrected-model pathway evidence matrix.
- `table_s28_kidney_spleen_manual_review.tsv`: all primary top-decile programs with complete-label tissue review and disposition.
- `table_s29_kidney_spleen_member_gene_support.tsv`: member-gene and decoder-weight support for reviewed corrected-model pathways.
- `table_s30_kidney_spleen_literature_sources.tsv`: primary literature links used for corrected-model interpretation.
- `table_s31_latent_mapping_coordinates.tsv.gz`: tissue-specific reference and query PCA coordinates, project labels, and 20-PC nearest-reference distances.
- `table_s32_latent_mapping_qc.tsv`: reference and query counts, explained variance, nearest-neighbor threshold, and query coverage for each main tissue.
- `table_s33_representative_program_sample_scores.tsv.gz`: decoder-oriented and project-centered sample scores for the representative programs in Figure S9.
- `table_s34_retained_pathway_member_gene_support.tsv`: retained-pathway member-gene direction fractions, measured-gene counts, GSEA support, and leading concordant genes.
- `table_s35_retained_pathway_member_gene_effects.tsv.gz`: complete project-balanced member-gene effects behind Table S34 and Figure 4.
- `figure_build_manifest.tsv`: final figure dimensions, file sizes, vector-copy status, and nonwhite-border check.
