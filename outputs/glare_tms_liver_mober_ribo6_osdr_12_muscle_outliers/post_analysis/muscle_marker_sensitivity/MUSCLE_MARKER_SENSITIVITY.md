# Muscle/Cytoskeleton Marker Sensitivity

Same GLARE-style balanced FLT/GC XGBoost verification, after removing muscle/cytoskeleton-associated genes.

```tsv
filter_name	genes_kept	genes_removed	elapsed_seconds	scheme	accuracy_mean	accuracy_std	f1_mean	f1_std	roc_auc_mean	roc_auc_std
baseline_all_genes	21010	0	40.9758620262146	gene_grouped_kfold_audit	0.8912660637791529	0.0021058373922540705	0.8786688714133462	0.002658293258499477	0.9762452026991975	0.0009069227968313274
baseline_all_genes	21010	0	40.97666907310486	random_kfold_glare	0.8879581151832461	0.0038101268511721016	0.8843507714648758	0.012993851470745329	0.9761855771581726	0.00105187847015215
remove_curated_20_skeletal_muscle_markers	20990	20	40.87403988838196	gene_grouped_kfold_audit	0.8914244878513579	0.002901113662056643	0.8787331679397354	0.003696080984888809	0.976289470405041	0.0012185067834538938
remove_curated_20_skeletal_muscle_markers	20990	20	40.87414884567261	random_kfold_glare	0.8898523106241066	0.003286071129977774	0.8818333669355839	0.010213351764430788	0.9758455719045385	0.001538433125145722
remove_metascape_muscle_cytoskeleton_hits	19929	1081	39.68536901473999	gene_grouped_kfold_audit	0.8858445147728468	0.0035555000818097973	0.8717199846926451	0.004343849170069933	0.9737951144888773	0.0016318210781477698
remove_metascape_muscle_cytoskeleton_hits	19929	1081	39.68546915054321	random_kfold_glare	0.883160261769285	0.00375074732274765	0.8744511837867723	0.009169987393651763	0.9732280537354214	0.001050443249723133
```

- Curated muscle marker genes removed: 20
- Broad Metascape muscle/cytoskeleton hit genes removed: 1081

Interpretation: compare the grouped-by-gene ROC AUC across filters; this is the conservative audit because held-out genes are not seen during classifier training.
