# Study/Batch Signal Audit

This audit asks whether the aggregate liver FLT-vs-GC signal looks like a
portable spaceflight transcriptomic signal or a study/mission/batch-specific
signal.

Input run:
`outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers`

Samples:

```tsv
accession	FLT	GC
OSD-137	5	6
OSD-173	2	2
OSD-242	5	3
OSD-245	18	19
OSD-379	31	31
OSD-463	12	10
```

## Main Result

The aggregate run contains a real expression signal, but it is not yet a
portable cross-study spaceflight signature.

- Random sample-level FLT-vs-GC expression CV is weak/moderate:
  ROC AUC `0.615`.
- Accession-centered random CV remains weak/moderate:
  ROC AUC `0.652`.
- Leave-one-accession-out FLT-vs-GC expression prediction fails:
  pooled ROC AUC `0.445`.
- Leave-one-accession-out after accession-centering also fails:
  pooled ROC AUC `0.465`.
- Expression predicts OSD accession strongly:
  4-fold accession classification accuracy `0.778`,
  balanced accuracy `0.686`.

Interpretation: within-study FLT/GC differences exist in some studies, but the
direction/signature does not transfer reliably to unseen OSD studies. Study or
mission-specific expression structure is stronger than the shared aggregate
spaceflight signature.

## FLT-vs-GC Sample-Level Tests

```tsv
test	n	accuracy	balanced_accuracy	f1	roc_auc
expression_random_stratified_5fold	144	0.5833333333333334	0.5840246961219371	0.5652173913043478	0.6154736639012155
expression_accession_centered_random_5fold	144	0.5763888888888888	0.5773683195060776	0.5481481481481482	0.6515531545437006
expression_leave_one_accession_out_pooled	144	0.4791666666666667	0.4805132162840054	0.42748091603053434	0.445109010225738
expression_accession_centered_leave_one_accession_out_pooled	144	0.4861111111111111	0.48543314682616245	0.5131578947368421	0.46517460929963345
metadata_only_random_stratified_5fold	144	0.4930555555555556	0.49035307736831946	0.5780346820809249	0.48196025467875747
metadata_only_leave_one_accession_out_pooled	144	0.5138888888888888	0.508392822689562	0.6534653465346535	0.5408064827320085
accession_only_random_stratified_5fold	144	0.4444444444444444	0.44375844105730267	0.47368421052631576	0.418869380667567
dummy_prior_random_5fold	144	0.5069444444444444	0.5	0.6728110599078341	0.48620490063669686
```

Metadata-only models used accession, project identifier, sex, strain,
genotype, age at launch, duration, library selection/layout, and sequencing
instrument. Sample/source names and all explicit condition labels were
excluded.

## Per-Study Generalization

```tsv
test	heldout_accession	n	accuracy	balanced_accuracy	f1	roc_auc	n_splits
expression_leave_one_accession_out	OSD-137	11	0.45454545454545453	0.4166666666666667	0.0	0.23333333333333334	
expression_leave_one_accession_out	OSD-173	4	0.25	0.25	0.4	0.25	
expression_leave_one_accession_out	OSD-242	8	0.625	0.6333333333333333	0.6666666666666666	0.6666666666666667	
expression_leave_one_accession_out	OSD-245	37	0.5135135135135135	0.5014619883040935	0.1	0.3654970760233918	
expression_leave_one_accession_out	OSD-379	62	0.41935483870967744	0.4193548387096774	0.4	0.3881373569198751	
expression_leave_one_accession_out	OSD-463	22	0.5909090909090909	0.5583333333333333	0.7096774193548387	0.525	
expression_within_accession_cv	OSD-137	11	0.45454545454545453	0.45	0.4	0.2	5.0
expression_within_accession_cv	OSD-173	4	0.5	0.5	0.5	0.25	2.0
expression_within_accession_cv	OSD-242	8	0.5	0.4666666666666667	0.6	0.7333333333333334	3.0
expression_within_accession_cv	OSD-245	37	0.6486486486486487	0.6491228070175439	0.6486486486486487	0.7865497076023391	5.0
expression_within_accession_cv	OSD-379	62	0.6451612903225806	0.6451612903225806	0.6451612903225806	0.6826222684703434	5.0
expression_within_accession_cv	OSD-463	22	0.5454545454545454	0.5416666666666667	0.5833333333333334	0.6916666666666667	5.0
```

The larger studies have some within-study FLT/GC signal, especially OSD-245,
OSD-379, and OSD-463. But models trained on the other studies do not predict
held-out OSD-245 or OSD-379 well. That is the key failure mode for a shared
aggregate spaceflight signature.

## Accession Predictability

Expression predicts OSD accession better than it predicts FLT-vs-GC:

```tsv
test	n	classes	accuracy	balanced_accuracy	macro_f1
expression_predict_accession_random_4fold	144	6	0.7777777777777778	0.6855168951943145	0.6691384983420381
```

Confusion matrix:

```tsv
	OSD-137	OSD-173	OSD-242	OSD-245	OSD-379	OSD-463
OSD-137	7	0	0	0	4	0
OSD-173	0	2	0	2	0	0
OSD-242	0	0	4	1	0	3
OSD-245	0	1	1	35	0	0
OSD-379	15	0	0	0	47	0
OSD-463	0	0	5	0	0	17
```

This agrees with the earlier cluster-level batch QC: accession/mission and
strain explain large fractions of cluster mean-expression variation in both
FLT and GC partitions.

## Relationship To The GLARE Verifier

The GLARE-original-style XGBoost verifier remains highly separable
(`ROC AUC ~0.976`), but it is a gene-row classifier: each row is a gene under
FLT or GC, and feature columns are matched sample slots. That verifies that the
FLT and GC expression matrices have separable gene-pattern structure.

The sample-level audit here asks a different question: can a model trained on
some biological samples or studies identify FLT-vs-GC in new samples or unseen
studies? For that question, leave-one-accession-out performance is poor.

## Conclusion

Current best interpretation:

- There is within-study transcriptomic FLT-vs-GC signal.
- There is strong study/mission/batch structure.
- The shared aggregate liver spaceflight signature is not robust enough yet to
  generalize across held-out OSD studies.
- Treat the current aggregate GLARE result as exploratory and batch-sensitive,
  not definitive evidence of a universal liver spaceflight response.

For a stronger biological claim, analyze controlled studies separately, then
meta-analyze pathway-level effects across studies rather than pooling all
samples into one aggregate classifier.
