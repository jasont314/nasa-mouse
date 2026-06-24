# GLARE on MOBER-Corrected Aggregate Liver Data

- MOBER input: `outputs/mober_liver_ribo6_osdr_12_muscle_outliers/projection/mober_projected_onto_OSD-379.h5ad`
- Projection target: `OSD-379`
- Shape: 144 samples x 21,010 genes

## Condition Counts

```tsv
h5_accession	FLT	GC
OSD-137	5	6
OSD-173	2	2
OSD-242	5	3
OSD-245	18	19
OSD-379	31	31
OSD-463	12	10
```

## Fine-Tuning

```tsv
location	genes	profiles	best_loss	best_epoch	epochs
FLT	21010	73	0.08380858	25	30
GC	21010	71	0.08114952	30	30
```
