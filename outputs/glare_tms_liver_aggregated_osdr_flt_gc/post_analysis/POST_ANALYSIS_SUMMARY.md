# Aggregated Liver Post-Analysis

## Study/Batch QC

Cluster mean expression was tested against accession, mission, library,
sex, strain, and sequencing instrument within FLT and GC separately.
`mission` is the OSDR project identifier field.

```tsv
location	variable	clusters_tested	median_eta_squared	max_eta_squared	strong_driver_clusters
FLT	h5_accession	16	0.6040957831022402	0.9411800252979861	14
FLT	library_layout	16	0.050545364271487286	0.08933089525746188	0
FLT	library_selection	16	0.2678224392664441	0.588310549415426	10
FLT	mission	16	0.6040957831022402	0.9411800252979861	14
FLT	sequencing_instrument	16	0.060649868010412986	0.13443131908421727	0
FLT	sex	16	0.003835155196789554	0.8153395916286307	2
FLT	strain	16	0.10496459875240635	0.8909768549847394	2
GC	h5_accession	15	0.5740413315553256	0.9654684118780771	12
GC	library_layout	15	0.006924995333579297	0.06418716559954386	0
GC	library_selection	15	0.1866234679013586	0.4634718419655892	5
GC	mission	15	0.5740413315553256	0.9654684118780771	12
GC	sequencing_instrument	15	0.03230167413728423	0.06990306550268909	0
GC	sex	15	0.018204084490171184	0.8426523671371801	1
GC	strain	15	0.1701109925074266	0.9131987194471801	3
```

## FLT vs GC Cluster Structure

- Genes compared: 21,010
- Adjusted Rand index: 0.4202
- Normalized mutual information: 0.6098
- Median Procrustes latent shift: 0.2205

## Exploratory Meta-DGEA

This is not DESeq2 from raw counts. It is an exploratory fixed-effect
meta-analysis over per-study Welch tests on log2(normalized expression + 1)
from the integrated HDF5/aligned target matrix.

- Accessions tested: OSD-137, OSD-173, OSD-242, OSD-245, OSD-379, OSD-463, OSD-47, OSD-686
- Genes tested: 21,010
- Significant genes at FDR < 0.05 and abs(log2FC) >= 1: 30
- Up: 8
- Down: 22

Top meta-analysis genes are in `top_meta_dgea_genes.tsv`.

## Raw-Count DESeq2 Inputs

The HDF5 expression matrix is integer count-like data. Per-study DESeq2
inputs were exported for a study-aware FLT-vs-GC analysis:

- Counts: `outputs/glare_tms_liver_aggregated_osdr_flt_gc/post_analysis/deseq2_inputs/counts.tsv`
- Metadata: `outputs/glare_tms_liver_aggregated_osdr_flt_gc/post_analysis/deseq2_inputs/sample_metadata.tsv`
- Gene symbols: `outputs/glare_tms_liver_aggregated_osdr_flt_gc/post_analysis/deseq2_inputs/gene_symbols.tsv`
