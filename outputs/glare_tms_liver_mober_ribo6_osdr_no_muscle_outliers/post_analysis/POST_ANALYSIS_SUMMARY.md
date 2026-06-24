# Aggregated Liver Post-Analysis

## Study/Batch QC

Cluster mean expression was tested against accession, mission, library,
sex, strain, and sequencing instrument within FLT and GC separately.
`mission` is the OSDR project identifier field.

```tsv
location	variable	clusters_tested	median_eta_squared	max_eta_squared	strong_driver_clusters
FLT	h5_accession	16	0.48185603257631504	0.832277240568975	14
FLT	library_layout	16			0
FLT	library_selection	16	0.07884673046141108	0.6315259285670453	4
FLT	mission	16	0.48185603257631504	0.832277240568975	14
FLT	sequencing_instrument	16	0.016247003272984766	0.07569606549052883	0
FLT	sex	16	0.12151917997465067	0.6014103293273758	6
FLT	strain	16	0.400799845933852	0.748876757203388	12
GC	h5_accession	15	0.5252744187054522	0.9325261600156073	11
GC	library_layout	15			0
GC	library_selection	15	0.04169805711108807	0.5274259956551978	2
GC	mission	15	0.5252744187054522	0.9325261600156073	11
GC	sequencing_instrument	15	0.011763702311501496	0.04196789995171593	0
GC	sex	15	0.17834491072766373	0.6736235052739846	4
GC	strain	15	0.4390549525565892	0.6819286642135689	10
```

## FLT vs GC Cluster Structure

- Genes compared: 21,010
- Adjusted Rand index: 0.7170
- Normalized mutual information: 0.7445
- Median Procrustes latent shift: 1.1598

## Exploratory Meta-DGEA

This is not DESeq2 from raw counts. It is an exploratory fixed-effect
meta-analysis over per-study Welch tests on log2(normalized expression + 1)
from the integrated HDF5/aligned target matrix.

- Accessions tested: OSD-137, OSD-173, OSD-242, OSD-245, OSD-379, OSD-463
- Genes tested: 21,010
- Significant genes at FDR < 0.05 and abs(log2FC) >= 1: 0
- Up: 0
- Down: 0

Top meta-analysis genes are in `top_meta_dgea_genes.tsv`.

## Raw-Count DESeq2 Inputs

The HDF5 expression matrix is integer count-like data. Per-study DESeq2
inputs were exported for a study-aware FLT-vs-GC analysis:

- Counts: `outputs/glare_tms_liver_mober_ribo6_osdr_no_muscle_outliers/post_analysis/deseq2_inputs/counts.tsv`
- Metadata: `outputs/glare_tms_liver_mober_ribo6_osdr_no_muscle_outliers/post_analysis/deseq2_inputs/sample_metadata.tsv`
- Gene symbols: `outputs/glare_tms_liver_mober_ribo6_osdr_no_muscle_outliers/post_analysis/deseq2_inputs/gene_symbols.tsv`
