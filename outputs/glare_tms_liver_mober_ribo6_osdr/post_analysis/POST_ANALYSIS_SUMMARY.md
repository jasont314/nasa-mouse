# Aggregated Liver Post-Analysis

## Study/Batch QC

Cluster mean expression was tested against accession, mission, library,
sex, strain, and sequencing instrument within FLT and GC separately.
`mission` is the OSDR project identifier field.

```tsv
location	variable	clusters_tested	median_eta_squared	max_eta_squared	strong_driver_clusters
FLT	h5_accession	16	0.8179324888170041	0.9981291893038707	16
FLT	library_layout	16			0
FLT	library_selection	16	0.0354808002188126	0.09855834292479994	0
FLT	mission	16	0.8179324888170041	0.9981291893038707	16
FLT	sequencing_instrument	16	0.03091760742332091	0.38786703581715737	1
FLT	sex	16	0.560228988487533	0.9976567727236778	8
FLT	strain	16	0.8133176021837443	0.9981291893038707	16
GC	h5_accession	15	0.7846828360724234	0.998766551188354	15
GC	library_layout	15			0
GC	library_selection	15	0.057441433863495844	0.14547864210083483	0
GC	mission	15	0.7846828360724234	0.998766551188354	15
GC	sequencing_instrument	15	0.027776051106069132	0.10644193778190508	0
GC	sex	15	0.48046099733301356	0.9982993792252786	8
GC	strain	15	0.7846472308568385	0.998766551188354	15
```

## FLT vs GC Cluster Structure

- Genes compared: 21,010
- Adjusted Rand index: 0.5770
- Normalized mutual information: 0.7757
- Median Procrustes latent shift: 1.4535

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

- Counts: `outputs/glare_tms_liver_mober_ribo6_osdr/post_analysis/deseq2_inputs/counts.tsv`
- Metadata: `outputs/glare_tms_liver_mober_ribo6_osdr/post_analysis/deseq2_inputs/sample_metadata.tsv`
- Gene symbols: `outputs/glare_tms_liver_mober_ribo6_osdr/post_analysis/deseq2_inputs/gene_symbols.tsv`
