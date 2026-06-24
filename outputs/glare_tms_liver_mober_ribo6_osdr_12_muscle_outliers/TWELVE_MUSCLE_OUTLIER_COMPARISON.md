# Twelve Muscle-Outlier Filter Sensitivity

This run removes the 12 candidate muscle-composition outliers identified in
the full six-study aggregate liver QC:

- the original 8 severe OSD-379 outliers;
- 2 OSD-245 FLT profiles;
- 1 OSD-137 FLT profile;
- 1 OSD-242 GC profile.

The run reuses the same TMS liver pretrained GLARE SAE, the same six OSDR
liver accessions, and the same `OSD-379` MOBER projection target. MOBER was
trained with batch size `24` instead of `32` because the 144-sample filtered
split produced a singleton final mini-batch at batch size `32`, which is
invalid for batch normalization.

## Sample Counts

```tsv
run	samples	FLT	GC	OSD-137_FLT	OSD-137_GC	OSD-242_FLT	OSD-242_GC	OSD-245_FLT	OSD-245_GC	OSD-379_FLT	OSD-379_GC	OSD-463_FLT	OSD-463_GC
baseline	156	80	76	6	6	5	4	20	19	35	35	12	10
8_removed	148	76	72	6	6	5	4	20	19	31	31	12	10
12_removed	144	73	71	5	6	5	3	18	19	31	31	12	10
```

## Residual Muscle QC

```tsv
run	high_muscle_abundance	candidate_muscle_outliers	broad_review_muscle_outliers
baseline	19	12	6
8_removed	11	6	5
12_removed	7	2	5
```

After removing all 12 candidates, two lower-abundance OSD-379 profiles remain
candidate outliers under the same robust-z rule:

- `RR8_LVR_GC_LAR_OLD_GL7`
- `RR8_LVR_FLT_ISS-T_OLD_FI12`

The five broad-review profiles are all OSD-379 LAR FLT profiles. They are
high relative to broad OSD-379 groups, but not high relative to their
collection/age strata.

## MOBER QC

Lower `data_source` silhouette after projection is consistent with reduced
accession/batch separation.

```tsv
run	space	data_source_silhouette	condition_silhouette	sex_silhouette	strain_silhouette
baseline	mober_projected_onto_OSD-379	-0.1392	-0.0085	0.5773	0.0036
8_removed	mober_projected_onto_OSD-379	0.1413	-0.0069	0.2742	0.1198
12_removed	mober_projected_onto_OSD-379	-0.0172	0.0011	0.2263	0.0245
```

## GLARE Fine-Tuning And Clustering

```tsv
run	FLT_profiles	GC_profiles	FLT_best_loss	GC_best_loss	FLT_clusters	FLT_silhouette	GC_clusters	GC_silhouette	FLT_GC_ARI	FLT_GC_NMI	median_latent_shift
baseline	80	76	0.08153502	0.08426146	16	0.2792	15	0.4334	0.5770	0.7757	1.4535
8_removed	76	72	0.08221147	0.08075465	16	0.2937	15	0.2732	0.7170	0.7445	1.1598
12_removed	73	71	0.08380858	0.08114952	16	0.3244	15	0.3032	0.6935	0.6997	1.2480
```

The 12-filter run keeps the same consensus cluster counts as the baseline and
8-filter runs. FLT cluster silhouette improves monotonically across the filters;
GC silhouette drops after filtering but is slightly better than the 8-filter
run. FLT/GC cluster agreement remains higher than baseline but lower than the
8-filter run.

## Batch/Study Driver Check

```tsv
run	location	h5_accession_median_eta	h5_accession_max_eta	strong_accession_driver_clusters
baseline	FLT	0.8179	0.9981	16
baseline	GC	0.7847	0.9988	15
8_removed	FLT	0.4819	0.8323	14
8_removed	GC	0.5253	0.9325	11
12_removed	FLT	0.3545	0.7385	14
12_removed	GC	0.4452	0.8780	14
```

The 12-filter run reduces median accession effect size further than the
8-filter run, but many clusters are still accession-driven. That means batch
and study structure remain a major interpretation caveat.

## Raw-Count DESeq2 Meta-Analysis

These numbers use raw-count TSV inputs exported from the MOBER prep directories,
not the MOBER-projected expression used for GLARE fine-tuning.

```tsv
run	genes_tested	genes_eligible_min_2_studies	significant_fdr05_abs_log2fc1	up	down
baseline	16521	15683	38	28	10
8_removed	16276	15631	57	45	12
12_removed	16255	15597	47	36	11
```

Top 12-filter significant genes still include muscle/contractile genes:
`Myl1`, `Tnnt3`, `Acta1`, `Myh4`, `Myh2`, `Mybpc1`, and `Myh1`. The muscle
signature is weaker than the 8-filter run but still present, so removing the
12 candidate samples does not eliminate the muscle/contractile DGEA pattern.

## GLARE/DESeq2 Overlap

Corrected raw-count DESeq2 hits were joined to GLARE consensus clusters:

- Eligible genes with GLARE clusters: 15,597
- Significant DESeq2 meta genes: 47
- Strongest FLT enrichment: cluster 13, 11 significant genes, Fisher FDR
  `2.70e-08`, top genes include `Myl1`, `Tnnt3`, `Acta1`, `Myh4`, and `Myh2`.
- Strongest GC enrichment: cluster 1, 29 significant genes, Fisher FDR
  `7.20e-20`, top genes include `Myl1`, `Tnnt3`, `Acta1`, `Tg`, `Myh4`, and
  `Myh2`.

## Interpretation

The 12-filter sensitivity run is useful and should be kept, but it should not
be treated as a clean final answer by itself. It reduces accession effects and
removes the non-OSD-379 candidate muscle outliers, but residual muscle/contractile
signal remains in both DGEA and GLARE cluster overlap.

Pragmatic next comparison:

1. keep baseline, 8-filter, and 12-filter as sensitivity tiers;
2. do not automatically remove the five broad-review OSD-379 LAR FLT profiles;
3. separately test a 14-filter run only if we want to evaluate the two residual
   lower-abundance candidate profiles.
