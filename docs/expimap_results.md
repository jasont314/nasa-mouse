# expiMap Results Summary

Generated from current local manifests and ignored output artifacts.
This report intentionally distinguishes pipeline validation from final biological inference.

## OSDR API Discovery

OSDR inputs are built from the NASA OSDR Biological Data API, not from the older local integrated OSDR HDF5.

- API-selected samples: 1,631
- OSD accessions: 75
- count files: 75
- tissues: 24

| tissue | flight | ground_control | total |
| --- | --- | --- | --- |
| liver | 125 | 118 | 243 |
| skeletal_muscle | 95 | 96 | 191 |
| skin | 80 | 71 | 151 |
| kidney | 68 | 67 | 135 |
| thymus | 63 | 54 | 117 |
| spleen | 55 | 54 | 109 |
| lung | 40 | 38 | 78 |
| retina | 45 | 31 | 76 |
| brain | 28 | 29 | 57 |
| cerebellum | 29 | 27 | 56 |
| colon | 28 | 27 | 55 |
| heart | 26 | 25 | 51 |

Top additional tissues by available FLT/GC samples:

| tissue | flight | ground_control | total |
| --- | --- | --- | --- |
| skeletal_muscle | 95 | 96 | 191 |
| skin | 80 | 71 | 151 |
| thymus | 63 | 54 | 117 |
| spleen | 55 | 54 | 109 |
| lung | 40 | 38 | 78 |
| retina | 45 | 31 | 76 |
| brain | 28 | 29 | 57 |
| cerebellum | 29 | 27 | 56 |

## Direct OSDR expiMap

These runs train expiMap directly on each tissue-specific OSDR FLT/GC dataset.
The current direct summaries use matched 50-epoch runs. Raw counts with negative-binomial loss are the primary analysis; CPM and log1p(CPM) use MSE loss and are sensitivity checks.

| tissue | samples | flight | ground | accessions | genes | pathways | input |
| --- | --- | --- | --- | --- | --- | --- | --- |
| liver | 231 | 118 | 113 | 12 | 9,321 | 1,140 | `outputs/expimap_direct_osdr_liver/input` |
| kidney | 135 | 68 | 67 | 6 | 9,321 | 1,140 | `outputs/expimap_direct_osdr_kidney/input` |

| tissue | transform | role | loss | epochs | min Welch FDR | Welch FDR<0.10 | min study FDR | study FDR<0.10 | top aggregate term |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| liver | raw_counts | primary | nb | 50 | 0.012 | 1 | 0.85 | 0 | R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION |
| liver | cpm | sensitivity | mse | 50 | 0.0111 | 1 | 0.827 | 0 | R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION |
| liver | log1p_cpm | sensitivity | mse | 50 | 0.0113 | 1 | 0.85 | 0 | R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION |
| kidney | raw_counts | primary | nb | 50 | 0.612 | 0 | 0.959 | 0 | R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES |
| kidney | cpm | sensitivity | mse | 50 | 0.631 | 0 | 0.959 | 0 | R-MMU-9933947_FORMATION_OF_THE_NON_CANONICAL_BAF_NCBAF_COMPLEX |
| kidney | log1p_cpm | sensitivity | mse | 50 | 0.621 | 0 | 0.967 | 0 | R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES |

Current direct-run interpretation:

- Each current direct analysis directory contains PCA and UMAP coordinates plus condition/accession-colored plots.
- Liver has one aggregate FDR-significant pathway in the historical 50-epoch stochastic-score runs: RNA Polymerase II transcription elongation is lower in flight and is the top term across raw counts, CPM, and log1p(CPM).
- Kidney has no aggregate pathway FDR < 0.10 in the matched 50-epoch direct runs.
- 50-epoch effect rankings are highly correlated across raw-count NB and CPM/log1p(CPM) MSE sensitivity runs, unlike the earlier 3-epoch validation runs.
- Posterior-mean accession-aware validation does not support the Polymerase II result; see [accession validation](expimap_accession_validation.md).

## Preprocessing Comparison

### liver

| transform | loss | validity | min Welch FDR | Welch FDR<0.10 | top term |
| --- | --- | --- | --- | --- | --- |
| raw_counts | nb | primary_count_likelihood | 0.012 | 1 | R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION |
| cpm | mse | sensitivity_only_not_count_likelihood | 0.0111 | 1 | R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION |
| log1p_cpm | mse | sensitivity_only_not_count_likelihood | 0.0113 | 1 | R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION |

| comparison | Spearman effect rho | top50 overlap | top50 Jaccard |
| --- | --- | --- | --- |
| raw_counts vs cpm | 0.983 | 24 | 0.316 |
| raw_counts vs log1p_cpm | 0.996 | 37 | 0.587 |
| cpm vs log1p_cpm | 0.985 | 26 | 0.351 |

### kidney

| transform | loss | validity | min Welch FDR | Welch FDR<0.10 | top term |
| --- | --- | --- | --- | --- | --- |
| raw_counts | nb | primary_count_likelihood | 0.612 | 0 | R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES |
| cpm | mse | sensitivity_only_not_count_likelihood | 0.631 | 0 | R-MMU-9933947_FORMATION_OF_THE_NON_CANONICAL_BAF_NCBAF_COMPLEX |
| log1p_cpm | mse | sensitivity_only_not_count_likelihood | 0.621 | 0 | R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES |

| comparison | Spearman effect rho | top50 overlap | top50 Jaccard |
| --- | --- | --- | --- |
| raw_counts vs cpm | 0.986 | 25 | 0.333 |
| raw_counts vs log1p_cpm | 0.996 | 34 | 0.515 |
| cpm vs log1p_cpm | 0.985 | 24 | 0.316 |

## ARCHS4 Reference-Query expiMap

These runs train a tissue-filtered ARCHS4 mouse bulk reference and map the OSDR tissue dataset as query.
The current liver/kidney reference-query summaries use 1000 ARCHS4 tissue-filtered reference samples, 50 reference-training epochs, and 50 query-mapping epochs.

| tissue | ARCHS4 ref samples | ref epochs | ref genes | pathways | query samples | query genes | query epochs | min Welch FDR | Welch FDR<0.10 | top query term |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| liver | 1,000 | 50 | 9,319 | 1,140 | 231 | 9,319 | 50 | 0.894 | 0 | R-MMU-5654706_FRS_MEDIATED_FGFR3_SIGNALING |
| kidney | 1,000 | 50 | 9,319 | 1,140 | 135 | 9,319 | 50 | 0.777 | 0 | R-MMU-4641262_DISASSEMBLY_OF_THE_DESTRUCTION_COMPLEX_AND_RECRUITMENT_OF_AXIN_TO_THE_MEMBRANE |

Reference-query interpretation:

- These runs use tissue-filtered, leakage-excluded ARCHS4 references with the same Reactome architecture as the direct OSDR runs.
- Reference-query preprocessing is raw-count NB only in the current workflow; CPM/log1p(CPM) comparisons are direct-workflow sensitivity analyses.
- Each current query analysis directory contains PCA and UMAP coordinates plus condition/accession-colored plots.
- They are still bounded 1000-sample reference runs rather than all available ARCHS4 tissue samples, so compare them with direct OSDR results before treating a pathway as robust.

## Direct vs Reference-Query Agreement

| tissue | direct raw-count top term | direct effect | direct Welch FDR | reference-query effect | reference-query Welch FDR | same direction | reference FDR<0.10 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| liver | R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION | -0.538 | 0.012 | -0.179 | 0.972 | yes | no |
| kidney | R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES | -0.358 | 0.612 | 0.119 | 0.978 | no | no |

Current workflow-agreement interpretation:

- The direct liver top term has the same negative direction in the bounded ARCHS4 reference-query run, but it is not FDR-significant there.
- Kidney has no aggregate FDR-significant direct or reference-query pathway signal.
- Treat the direct liver signal as preprocessing-stable but not yet reference-query-confirmed.

## Posterior-Mean Accession Validation

The historical direct and bounded-reference tables above used stochastic latent
scores. Posterior means are now the default for future workflow outputs.

For the direct liver raw-count NB model, accession-aware random-effects
validation of the 12 OSDR accessions does not support Polymerase II transcription
elongation: FDR = 0.670, I2 = 0.904, and only seven of 12 accession effects
match the pooled meta-analysis direction. The 381 posterior-mean terms with
meta-analysis FDR < 0.05 remain model-score candidates pending direct-model
seed and count-level replication. Details and artifacts are in
[accession validation](expimap_accession_validation.md).

## 5,000-Sample ARCHS4 Liver Reference

A larger proportional series-stratified liver reference now uses 5,000 of 8,970
eligible nonleakage ARCHS4 samples across all 518 eligible series. Three GPU
reference seeds were trained for up to 200 epochs with early stopping, then each
mapped the 231-sample OSDR liver query for 50 epochs.

One term, `R-MMU-416700_OTHER_SEMAPHORIN_INTERACTIONS`, is FDR-significant in
all three seeds and every leave-one-accession-out run, with lower flight scores.
It is only a follow-up candidate: its accession heterogeneity remains high and
the all-term seed-to-seed effect ranks are nearly uncorrelated. Polymerase II
transcription elongation is neither significant nor direction-stable in these
three reference-query runs. See [reference seed stability](expimap_reference_seed_stability.md).

## Condition-Specific Clustering

As an expiMap analogue of separately clustering GC and FLT states, posterior-mean
liver pathway scores were clustered independently within each condition and then
compared by nearest pathway-score centroid.

The direct OSDR model selected two GC clusters and three FLT clusters. The two
large FLT clusters mostly mirror the same accession composition as the two GC
clusters, while the extra FLT cluster contains only five samples. The larger
ARCHS4 reference-query runs split FLT into one large cluster plus a small
mostly OSD-379 cluster across all three seeds, while GC cluster counts and
pathway shifts vary by seed. This makes the analysis useful QC, but not a new
confirmed pathway result. See [condition clustering](expimap_condition_clustering.md).

## Tutorial-Style Liver Run

A separate liver run now follows the scArches expiMap tutorial mechanics more
closely: 2,000 reference HVGs selected with `batch_key="archs4_condition"`,
strict post-HVG Reactome term filtering, a `300,300,300` reference network,
400-epoch reference training with early stopping, and 250-epoch query surgery
with three unconstrained de novo nodes plus HSIC one-vs-all regularization.

This retained 1,995 genes and 364 Reactome terms. The reference early-stopped
after 248 epochs, best epoch 196. Query surgery completed all 250 epochs after
a repo runtime patch stabilized the installed scArches HSIC gamma-ratio
calculation.

The run produced 70 accession-aware random-effects meta-analysis terms with
FDR < 0.05, and 32 also passed every leave-one-accession-out FDR < 0.05 in the
same meta-analysis direction. However, the study-aware Wilcoxon test across
accession effects did not reach FDR 0.05, expiMap `latent_enrich` condition
Bayes-factor-like scores were weak, and the three de novo programs were not
FLT/GC-significant. See [tutorial-style liver run](expimap_tutorial_style_liver.md).

## ARCHS4 Tissue Candidates

| tissue | usable nonleakage bulk-like samples |
| --- | --- |
| liver | 8,970 |
| spleen | 6,289 |
| lung | 5,674 |
| skin | 2,593 |
| kidney | 2,464 |
| skeletal_muscle | 1,412 |

## Current Bottom Line

- API-derived OSDR liver and kidney direct expiMap analyses are implemented and runnable.
- The historical pooled direct liver score nominated Polymerase II transcription elongation, but posterior-mean accession-aware validation does not support it.
- The matched 50-epoch direct kidney runs do not show aggregate FDR-significant FLT-vs-GC pathway shifts.
- CPM/log1p(CPM) sensitivity results are now rank-stable against raw-count NB for the matched 50-epoch direct runs.
- A 5,000-sample, three-seed ARCHS4 liver reference-query experiment identifies one stringent semaphorin candidate, but it is not yet a biological conclusion.
- GC/FLT condition-specific expiMap clustering mostly recovers accession/sample
  substructure rather than a clean flight-only pathway program.
- The tutorial-style HVG/deep/HSIC run is now complete; it generates candidate
  lipid/chromatin/ECM/Rho pathways but does not yet pass the full robustness bar.

## Next Full-Run Gate

Before making stronger scientific claims, replicate posterior-mean results across
direct-model seeds and fit a count-level, accession-aware model. Validate the
semaphorin candidate with an independent cohort or a materially different
reference sampling strategy.
