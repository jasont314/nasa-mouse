# MOBER Liver Ribo-Depletion Six-Dataset Input

Prepared aggregate OSDR liver FLT/GC expression for MOBER.

- Input h5ad: `outputs/mober_liver_ribo6_osdr_no_muscle_outliers/mober_liver_ribo6_input.h5ad`
- Shape: 148 samples x 21,010 genes
- Normalization: log2(CPM+1) from OSDR count-like HDF5 expression
- Batch/data_source column: `h5_accession`
- Data sources: OSD-137, OSD-173, OSD-242, OSD-245, OSD-379, OSD-463
- Excluded profiles matched: 8 of 8 requested

```tsv
h5_accession	FLT	GC
OSD-137	6	6
OSD-173	2	2
OSD-242	5	4
OSD-245	20	19
OSD-379	31	31
OSD-463	12	10
```
