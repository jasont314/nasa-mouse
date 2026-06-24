# Aggregated Liver Post-Analysis

## Study/Batch QC

Cluster mean expression was tested against accession, mission, library,
sex, strain, and sequencing instrument within FLT and GC separately.
`mission` is the OSDR project identifier field.

```tsv
location	variable	clusters_tested	median_eta_squared	max_eta_squared	strong_driver_clusters
FLT	h5_accession	16	0.35453496512984173	0.7384904604421718	14
FLT	library_layout	16			0
FLT	library_selection	16	0.11750152611594641	0.5100042871834795	1
FLT	mission	16	0.35453496512984173	0.7384904604421718	14
FLT	sequencing_instrument	16	0.0025393272227721498	0.025973046300187052	0
FLT	sex	16	0.16597475882480486	0.4446718217297849	5
FLT	strain	16	0.33876779304690663	0.7187373137390758	12
GC	h5_accession	15	0.4452274994088763	0.8779941371158937	14
GC	library_layout	15			0
GC	library_selection	15	0.11437719201421703	0.4790901022578162	4
GC	mission	15	0.4452274994088763	0.8779941371158937	14
GC	sequencing_instrument	15	0.014322278817618037	0.05727445981294702	0
GC	sex	15	0.11186290275024778	0.34516606046176	3
GC	strain	15	0.3861326708644622	0.8658051831793452	13
```

## FLT vs GC Cluster Structure

- Genes compared: 21,010
- Adjusted Rand index: 0.6935
- Normalized mutual information: 0.6997
- Median Procrustes latent shift: 1.2480

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

- Counts: `outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/deseq2_inputs/counts.tsv`
- Metadata: `outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/deseq2_inputs/sample_metadata.tsv`
- Gene symbols: `outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/deseq2_inputs/gene_symbols.tsv`
