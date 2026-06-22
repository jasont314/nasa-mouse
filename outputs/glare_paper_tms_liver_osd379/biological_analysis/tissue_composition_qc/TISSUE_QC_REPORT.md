# OSD-379 Liver Tissue-Composition QC

## Conclusion

OSD-379 is a liver dataset: the study is titled "Transcriptional profiling of
livers", and the samples are annotated as left lobe of the liver. The
muscle-enriched GLARE result is nevertheless not robust evidence of a
spaceflight-induced liver program.

The signal is concentrated in 2 old ISS-terminal flight
samples: RR8_LVR_FLT_ISS-T_OLD_FI17, RR8_LVR_FLT_ISS-T_OLD_FI16. These profiles retain strong liver-marker
expression, so they are not simply mislabeled muscle samples. They instead look
like liver RNA mixed with a substantial amount of skeletal-muscle RNA.

## Evidence

- FLT cluster 12 contains 115 genes.
- NASA's official old ISS-terminal contrast calls
  93/115 of them
  significant.
- Across all four matched contrasts, the cluster significance counts are
  93 (ISS-terminal old), 15 (live-animal-return young), 1 (live-animal-return old), 0 (ISS-terminal young).
- Before excluding the severe composition outliers,
  82
  cluster genes have an absolute mean log2 difference of at least 1.
- After excluding them, that falls to
  11;
  0
  cluster genes pass FDR < 0.05 and absolute effect >= 1 in the normalized-count
  sensitivity test.
- Across all tested genes, 42 non-cluster-specific
  candidates still pass that sensitivity threshold. The outlier finding
  invalidates the muscle-cluster interpretation, not every possible
  spaceflight response in the old ISS-terminal stratum.
- The implicated samples still express liver markers strongly:
  RR8_LVR_FLT_ISS-T_OLD_FI17: muscle score 13.92, liver score 16.36, RR8_LVR_FLT_ISS-T_OLD_FI16: muscle score 13.16, liver score 16.29.
- In 2,859 TMS FACS liver cells, the corresponding fast-muscle
  markers are sparse rather than a coherent liver-cell program. The median
  detection fraction for the marker panel is
  0.0017.
- Similar high muscle-marker scores occur in ground, vivarium, and baseline
  samples, including RR8_LVR_GC_LAR_YNG_GL5. That broader distribution argues
  against a flight-specific hepatic program and suggests a recurrent
  tissue-composition problem.
- RR8_LVR_FLT_ISS-T_OLD_FI16 had RIN 6.4, library prep 24-Feb-21, 66,531,392 reads, and 3.45% rRNA contamination. RR8_LVR_FLT_ISS-T_OLD_FI17 had RIN 6.4, library prep 29-Oct-20, 59,533,317 reads, and 0.64% rRNA contamination. The implicated samples were prepared in different library
  batches, making one shared library-preparation spillover event less likely.

## Recommended Filter

For the filtered GLARE and DESeq2 rerun, exclude profiles meeting both
predeclared criteria:

- At least 10/20 skeletal-muscle markers above 100
  normalized counts.
- At least 0.005% of total normalized abundance
  assigned to the 20-marker panel.

This directly flags 14 FLT/GC profiles. The filtered
analysis removes those profiles independently rather than discarding clean
animals from the opposite condition. The matched-slot table remains available
for methods that specifically require equal feature dimensions.

## Likely Cause

The most likely explanation is physical tissue admixture during dissection or
sample trimming. ISS-terminal animals were frozen as whole carcasses, later
thawed for approximately 60-90 minutes, and then dissected. Liver material was
trimmed to approximately 25 mg for RNA extraction. A small fragment from the
adjacent diaphragm or body-wall muscle can contribute a large muscle RNA signal
while leaving abundant liver transcripts intact.

This analysis cannot identify the exact contamination stage. Cross-sample
carryover during tissue handling, homogenization, or library preparation
remains possible. Normal read depth and mapping metrics do not rule out either
form of tissue-composition contamination.

## Impact On GLARE

The SAE and clustering captured a real expression module present in the input;
the muscle enrichment is not a software error. The biological interpretation
is the problem: the official DEG union is dominated by one contrast whose
effect is strongly leveraged by a few composition outliers.

Keep the current run as an unfiltered audit trail. For biological conclusions,
rerun the OSD-379 differential expression and GLARE fine-tuning after applying
a documented tissue-composition QC rule, and compare the filtered result with
the original.

## Files

- `sample_muscle_marker_scores.tsv`: all 141 liver profiles and marker scores.
- `sample_muscle_marker_scores.png`: cohort-wide score distribution.
- `old_iss_muscle_marker_heatmap.png`: marker expression in the affected stratum.
- `flt_cluster_12_significance_by_contrast.tsv`: localization by NASA contrast.
- `flt_cluster_12_old_iss_sensitivity.tsv`: before/after exclusion effects.
- `old_iss_cleaned_candidate_genes.tsv`: candidates retained after exclusion.
- `tms_liver_marker_detection.tsv`: single-cell liver reference detection.
- `implicated_sample_assay_qc.tsv`: official ISA assay metadata for FI16/FI17.
- `recommended_sample_exclusions.tsv`: profiles directly meeting the rule.
- `recommended_matched_slot_exclusions.tsv`: balanced slots removed downstream.

## Sources

- OSD-379: https://osdr.nasa.gov/bio/repo/data/studies/OSD-379
- RRRM-1/RR-8 payload: https://osdr.nasa.gov/bio/repo/data/payloads/RRRM-1%20%28RR-8%29
- RR-8 tissue handling details: https://www.nature.com/articles/s41467-026-68737-1
