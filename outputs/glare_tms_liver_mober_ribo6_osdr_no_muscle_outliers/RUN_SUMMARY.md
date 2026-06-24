# GLARE on MOBER-Corrected Aggregate Liver Data

- MOBER input: `outputs/mober_liver_ribo6_osdr_no_muscle_outliers/projection/mober_projected_onto_OSD-379.h5ad`
- Projection target: `OSD-379`
- Shape: 148 samples x 21,010 genes

## Condition Counts

```tsv
h5_accession	FLT	GC
OSD-137	6	6
OSD-173	2	2
OSD-242	5	4
OSD-245	20	19
OSD-379	31	31
OSD-463	12	10
```

## Fine-Tuning

```tsv
location	genes	profiles	best_loss	best_epoch	epochs
FLT	21010	76	0.08221147	30	30
GC	21010	72	0.08075465	30	30
```
