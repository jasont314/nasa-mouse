# GLARE Original-Style Analysis

This run intentionally follows the original GLARE post-analysis style:
direct FLT and GC representations are analyzed separately, direct
consensus clusters are interpreted, and XGBoost/SHAP is used as a
verification study. It does not use the later paired-cluster split
analysis such as `FLT14_not_GC8`.

Because the 12-filter aggregate liver target has unequal profile counts
after filtering (`73` FLT and `71` GC), the verification feature matrix
was balanced within each OSD accession before running the GLARE-style
melted-data classifier. This keeps the method usable without mixing
studies as feature positions.

## Verification

- Genes: 21,010
- Balanced feature columns: 69
- Rows in melted classifier table: 42,020

```tsv
scheme	accuracy_mean	accuracy_std	f1_mean	f1_std	roc_auc_mean	roc_auc_std
gene_grouped_kfold_audit	0.8912660637791529	0.002105837392254	0.8786688714133462	0.0026582932584994	0.9762452026991976	0.0009069227968313
random_kfold_glare	0.8879581151832461	0.0038101268511721	0.8843507714648758	0.0129938514707453	0.9761855771581726	0.0010518784701521
```

## Direct Cluster Interpretation

- Direct FLT/GC clusters summarized: 31
- Metascape-eligible direct cluster lists: 29

```tsv
location	cluster	gene_count	eligible_dgea_genes	significant_dgea_genes	significant_up_genes	significant_down_genes	mean_meta_log2fc_sig	included_for_metascape	top_significant_genes
GC	1	2185	1205	29	23	6	1.48495545914316	True	Myl1,Tnnt3,Acta1,Tg,Myh4,Myh2,Mybpc1,Cox8b,Mucl1,Myh1,Xirp2,Casq2,Trpm5,Sh3bgr,Ptgs2,Cdk3,Trim72,Il6,Fgf10,Ckmt1
FLT	3	1968	1052	13	7	6	0.3967970373293361	True	Tg,Mucl1,Trpm5,Ptgs2,Cdk3,Wnt11,Il6,Fgf10,Prtn3,Cxcl2,Gm13490,Mypn,3830403N18Rik
FLT	13	326	314	11	10	1	2.4801709163962355	True	Myl1,Tnnt3,Acta1,Myh4,Myh2,Mybpc1,Myh1,Xirp2,Adipoq,Cyp2c53-ps,Ckmt1
FLT	15	278	199	8	8	0	1.528459834224285	True	Cox8b,Casq2,Sh3bgr,Trim72,Alox15,Myo18b,Hrc,Chrdl2
GC	7	1340	1326	7	5	2	0.5730478948302444	True	Treh,Wnt11,Pitx3,Adipoq,Cyp2c53-ps,Trim80,Tuba8
GC	8	1173	1173	4	2	2	-0.0920660124374049	True	Mup15,Cyp4a14,Mup17,Cdkn1a
FLT	9	913	913	4	3	1	0.6017808348822775	True	Treh,Pitx3,Trim80,Tuba8
FLT	14	303	303	4	2	2	0.001323952775635	True	Apoa4,Mup15,Mup17,Cdkn1a
GC	0	4424	11	3	2	1	0.4958350661358632	False	Barx2,4931429L15Rik,1700027J07Rik
FLT	0	4412	9	3	2	1	0.4958350661358632	False	Barx2,4931429L15Rik,1700027J07Rik
GC	2	2104	2104	1	1	0	1.05546561291316	True	Adcy1
FLT	2	2085	2085	1	1	0	1.3788833687061	True	Sult1e1
```

## Outputs

- `verification/xgboost_verification_summary.tsv`
- `verification/shap_feature_importance.tsv`
- `verification/shap_gene_condition.tsv`
- `verification/shap_beeswarm.png`
- `direct_clusters/direct_cluster_summary.tsv`
- `direct_clusters/metascape_gene_lists/metascape_direct_cluster_gene_lists.csv`
