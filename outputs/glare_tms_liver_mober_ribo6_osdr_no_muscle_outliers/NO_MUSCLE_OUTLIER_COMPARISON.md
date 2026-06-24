# MOBER-GLARE With OSD-379 Muscle Outliers Removed

This reruns the six-study ribo-depletion liver MOBER and GLARE workflow after
removing the 8 recommended OSD-379 FLT/GC skeletal-muscle composition outliers
from `data/filters/osd379_severe_muscle_outlier_profiles.txt`.

## Sample Set

```tsv
h5_accession	FLT	GC	total
OSD-379	31	31	62
OSD-245	20	19	39
OSD-463	12	10	22
OSD-242	5	4	9
OSD-137	6	6	12
OSD-173	2	2	4
```

Total samples: 148, down from 156 in the unfiltered six-study run.

## GLARE Fine-Tuning

```tsv
run	location	profiles	best_loss	best_epoch
unfiltered	FLT	80	0.08153502	29
unfiltered	GC	76	0.08426146	30
no_muscle_outliers	FLT	76	0.08221147	30
no_muscle_outliers	GC	72	0.08075465	30
```

## Consensus Clustering

```tsv
run	location	consensus_clusters	silhouette	hdbscan_hard_clusters	hdbscan_noise_points
unfiltered	FLT	16	0.2792308032512665	9	15437
unfiltered	GC	15	0.43341612815856934	15	19346
no_muscle_outliers	FLT	16	0.2937171161174774	11	8883
no_muscle_outliers	GC	15	0.2731858491897583	15	18299
```

FLT clustering improves slightly and HDBSCAN leaves fewer genes as noise. GC
silhouette drops, although HDBSCAN still finds 15 hard clusters.

## FLT/GC Cluster Agreement

```tsv
run	adjusted_rand_index	normalized_mutual_information	median_latent_shift	mean_latent_shift
unfiltered	0.5770211637870556	0.7757338989488549	1.4535023935240612	1.5832095348933521
no_muscle_outliers	0.7169817618492853	0.7445459020540262	1.1597717828690959	1.2386517704824067
```

Removing the outliers makes FLT and GC partitions more similar by ARI and
reduces the latent shift.

## Raw-Count DESeq2 Meta-Analysis

Raw-count DESeq2 was run from the pre-MOBER count matrix, not the
MOBER-projected expression.

```tsv
run	genes_tested	eligible_min_2_studies	significant_fdr05_abs_log2fc1	up	down
unfiltered	16521	15683	38	28	10
no_muscle_outliers	16276	15631	57	45	12
```

The muscle/contractile genes do not disappear. They become stronger in the
filtered raw-count meta-analysis:

```tsv
gene	unfiltered_log2fc	unfiltered_fdr	filtered_log2fc	filtered_fdr
Acta1	2.3207264006391	0.000233416924610688	4.9231242313581	4.29860682032181e-22
Tnnt3	2.36685551913223	0.000285963858008652	4.07415845656012	3.52983264539031e-13
Myh4	2.04388981998208	0.00480277112563862	3.07096510753825	3.49438867003617e-07
Myh1	1.42215129145568	0.00458413762321109	2.22111794346504	5.26719450947532e-07
Myot	3.14100092971512	0.00145129604289235	4.63214959468764	4.12524381295157e-07
Xirp2	2.32449638505732	0.000389694612286542	3.26178417834645	6.38249718003887e-08
```

## DESeq2/GLARE Overlap

Top no-outlier DE-enriched GLARE clusters:

```tsv
location	cluster	eligible_genes	significant_genes	fisher_fdr_bh	top_significant_genes
FLT	2	2015	25	1.2271293549514216e-07	Tnnc2,Mb,Tmem182,Myot,Myl1,Tnnt3,Tcap,Xirp2,Myh4,Trdn
FLT	1	1197	16	2.9018593136370334e-05	Mucl1,Asb11,Ucp3,Trpm5,4931429L15Rik,Cdk3,1700027J07Rik,Il6,Fgf10,Prtn3
GC	2	2671	29	3.77716101791875e-08	Tnnc2,Mb,Acta1,Myot,Myl1,Tnnt3,Tcap,Tg,Xirp2,Myh4
GC	3	1054	20	5.218402027468359e-09	Tmem182,Mucl1,Asb11,Ucp3,Trpm5,Tbx15,4931429L15Rik,Cdk3,1700027J07Rik,Il6
```

Conclusion: removing the recommended OSD-379 muscle-composition outlier samples
does not remove the muscle/contractile biology from the six-study aggregate. It
reduces some sample-level leverage and makes FLT/GC latent partitions more
similar, but the muscle module persists and is stronger by study-aware
raw-count DESeq2.
