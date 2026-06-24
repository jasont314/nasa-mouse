# Paired FLT/GC GLARE Cluster Report

This report compares flight and ground-control liver GLARE gene
representations after MOBER correction and removal of the 12 candidate
muscle-outlier profiles. FLT and GC were fine-tuned as separate SAE
models, so latent-space distances use the existing Procrustes alignment
rather than direct unaligned coordinate subtraction.

## Overall Structure

- Genes compared: 21,010
- FLT clusters: 16
- GC clusters: 15
- Adjusted Rand index: 0.6935
- Normalized mutual information: 0.6997
- Median Procrustes latent shift: 1.2480

A high match fraction means most genes in a FLT cluster remain together
in one GC cluster. A high latent shift means those genes move farther
after aligning the two SAE spaces. The strongest candidates are therefore
clusters with both split structure and high latent shift.

## Top FLT Reorganization Candidates

```tsv
flt_cluster	matched_gc_cluster	gene_count	matched_fraction	gc_clusters_spanned	mean_latent_shift	significant_dgea_genes	cluster_status	top_shift_genes
14	8	303	0.538	8	3.26	4	split_reorganized	Myh11,Snord17,Scarna6,Neb,Slc25a4,Snora23,Snora74a,Fhl1,Sult2a1,Snord22,Cyp3a44,Hao2
9	2	913	0.45	6	2.72	4	split_reorganized	Ttn,Tpm2,Des,Eno3,Synm,Rnu12,Pygm,Cryab,Xirp1,Eps8l1,Acta2,Snora44
6	2	1422	0.543	6	2.35	1	split_reorganized	Ankrd23,Atp1a2,Synpo2,Hk2,Tmem38a,Snhg20,Carmn,Speg,Lynx1,4933411K16Rik,Cyp2g1,Sptb
15	1	278	0.601	4	3.17	8	split_reorganized	Actn3,Myoz1,ENSMUSG00000074199,Myo18b,Casq1,Cacna1s,Asprv1,Apobec2,Myot,Hrc,Trdn,Tgm3
13	7	326	0.773	6	3.33	11	partially_split_shifted	Ckm,Acta1,Myh1,Tnnc2,Tnnt3,Tnni2,Myh8,Cmya5,Mb,Obscn,Atp2a1,Ryr1
2	4	2085	0.538	5	1.73	1	split_low_shift	Firre,Dpt,Lpl,Ppp1r12b,Flvcr2,Ctla2a,Fads3,Gt(ROSA)26Sor,Gm12250,4632427E13Rik,Tagln,Cebpb
5	11	1431	0.517	6	1.36	0	split_low_shift	Tpm1,Shfl,Rps16,Snrnp70,Depp1,Rps21,Celsr1,Eps8l2,Kank2,Amotl1,Tpmt,Rps15a
1	5	2137	0.589	6	1.58	0	split_low_shift	Rbis,Gm14391,Taf1d,Cyp4a32,Tomm5,Zeb2,Mafb,Dnajc15,Fus,Immp1l,Pet100,Ccnl2
4	10	1619	0.557	5	1.25	1	split_low_shift	9330159F19Rik,Atcayos,Mex3d,Cdc25b,Ly6m,Ccdc166,4930554C24Rik,Tmprss11d,Ripply3,Arhgef4,Trpv2,Entpd2
8	3	934	0.681	4	1.39	0	partially_split	Fam193b,Leng8,Josd2,1810019D21Rik,Mapk8ip3,Dynll1,Gsn,Tjp3,Dct,Snrpg,Id3,Mid1ip1
```

## Displaced FLT Modules

`FLT##_not_GC##` means genes from a FLT cluster that do not map to that
cluster's primary GC counterpart.

```tsv
module	gene_count	fraction_of_flt_cluster	gc_clusters_spanned	mean_latent_shift	significant_dgea_genes	top_shift_genes
FLT13_not_GC7	74	0.227	5	5.26	9	Acta1,Myh1,Tnnc2,Tnnt3,Tnni2,Mb,Tcap,Amy2a5,Xirp2,ENSMUSG00000044041,Mybpc2,Myh2,Myl1,Sln,Mybpc1,ENSMUSG00000041984,Myh4,Eef1a2,Myom2,Ldb3
FLT14_not_GC8	140	0.462	7	3.07	1	Myh11,Snord17,Neb,Snora74a,Fhl1,Sult2a1,Hao2,Sult2a2,Eno1b,Rps2,Dio1,Cyp2a22,Rpl22l1,Slc22a27,Ndufb4,Gpat3,Ankrd55,Snord15b,Gbp11,Vldlr
FLT9_not_GC2	502	0.55	5	2.91	4	Des,Synm,Rnu12,Pygm,Xirp1,Eps8l1,Snora44,Tuba8,Snora78,Snora7a,Gstp2,Snord118,Snora3,Scand1,Gvin1,Adgrf1,AI506816,Ctse,Acpp,H2ac18
FLT15_not_GC1	111	0.399	3	2.49	0	ENSMUSG00000074199,Asprv1,Tgm3,Dmkn,ENSMUSG00000059956,ENSMUSG00000044594,ENSMUSG00000031757,Aldh3a1,Smtnl1,Dsg1a,ENSMUSG00000042306,ENSMUSG00000032807,ENSMUSG00000043430,ENSMUSG00000017204,Actc1,ENSMUSG00000045545,Cbr2,Neurl1a,Apol6,Aloxe3
FLT6_not_GC2	650	0.457	5	2.4	0	Ankrd23,Carmn,Speg,4933411K16Rik,Cyp2g1,Vaultrc5,Limch1,Zfp651,Lmnb1,Med31,Eci3,A730063M14Rik,Slc22a29,Prss23,Cap2,Gbp2b,Gngt2,Slc1a4,Phf13,S100a9
FLT2_not_GC4	964	0.462	4	2.11	1	Dpt,Flvcr2,Ctla2a,Fads3,Gt(ROSA)26Sor,Gm12250,Tagln,Cebpb,Tgtp1,Slc25a48,Car1,Arrdc2,Mef2c,Mir99ahg,Zbp1,Tspyl4,Zbtb9,Stmn1,Cerkl,Acot2
FLT12_not_GC6	24	0.0403	3	1.81	0	Synpo,1110038B12Rik,Tuba1a,H1f3,Tmem141,Tsr3,9030025P20Rik,Slc5a6,Endog,Dedd2,Gstm5,Prex1,Haus7,Uba7,Tmed1,Limd2,Vasp,Dqx1,Ogfod2,Psmb9
FLT7_not_GC2	272	0.255	2	1.73	0	Slc4a1,Spta1,Gm5111,Dnah11,Gzmb,Nek2,Ccnb2,Ttll13,Noxo1,A630076J17Rik,D7Ertd443e,Zfp879,Ak7,Ube2c,Ncaph,H2-Bl,Hs3st3a1,BC049715,Chaf1b,Sectm1a
FLT8_not_GC3	298	0.319	3	1.61	0	Josd2,1810019D21Rik,Mapk8ip3,Dynll1,Tjp3,Dct,Snrpg,AI661453,1810037I17Rik,Plin4,Uqcc2,Tnk2,Pfkm,Atoh8,Capn15,Pnpla6,Zfp710,Adck5,Camta2,Cep250
FLT11_not_GC12	135	0.212	3	1.51	1	Rps27l,Rps12,Rpl41,Atp5md,Rps23,Cox7a2,Atp5l,Fabp2,Rpl31,Rps25,Cox7c,Rps27a,Ndufa4,Rps8,Cox7b,Rpl38,Rpl22,Lpin1,Rpl9,Rps27
```

## Highest-Shift FLT-to-GC Edges

```tsv
flt_cluster	gc_cluster	edge_role	gene_count	fraction_of_flt_cluster	fraction_of_gc_cluster	mean_latent_shift	significant_dgea_genes	top_shift_genes
13	0	secondary	2	0.00613	0.000452	8.1	0	ENSMUSG00000044041,ENSMUSG00000041984
13	1	secondary	32	0.0982	0.0146	6.69	9	Acta1,Myh1,Tnnc2,Tnnt3,Tnni2,Mb,Tcap,Xirp2,Mybpc2,Myh2,Myl1,Sln
15	0	secondary	18	0.0647	0.00407	6.16	0	ENSMUSG00000074199,Asprv1,Tgm3,Dmkn,ENSMUSG00000059956,ENSMUSG00000044594,ENSMUSG00000031757,Aldh3a1,Smtnl1,Dsg1a,ENSMUSG00000042306,ENSMUSG00000032807
13	10	secondary	17	0.0521	0.0163	4.35	0	Ldb3,Perm1,Itgb1bp2,Synpo2l,Myh7,Nos1,Abca12,Trim63,Aspn,Il36g,Sgcd,Igdcc4
13	14	secondary	2	0.00613	0.0238	4.21	0	Pnlip,Gm13889
14	4	secondary	4	0.0132	0.00249	3.74	0	Fhl1,Serpinb6a,Fkbp11,F830016B08Rik
14	6	secondary	36	0.119	0.0268	3.67	0	Myh11,Sult2a1,Hao2,Sult2a2,Eno1b,Cyp2a22,Slc22a27,Ankrd55,Vldlr,Sult3a2,Abhd1,Syt1
13	2	secondary	21	0.0644	0.00998	3.63	0	Amy2a5,Myom1,Actn2,Klhl31,Dtna,Mfap3l,Akap6,Klf4,Scn7a,Nol3,Snora31,Meg3
15	1	flt_primary_only	167	0.601	0.0764	3.62	8	Actn3,Myoz1,Myo18b,Casq1,Cacna1s,Apobec2,Myot,Hrc,Trdn,Cox8b,Cox6a2,Smpx
12	8	secondary	3	0.00504	0.00256	3.47	0	Synpo,1110038B12Rik,H1f3
14	8	flt_primary_only	163	0.538	0.139	3.42	3	Scarna6,Slc25a4,Snora23,Snord22,Cyp3a44,Cyp3a16,Snord104,Cyp2a4,Cyp2b13,Cyp2b9,Ifi27l2a,Sult3a1
9	14	secondary	6	0.00657	0.0714	3.38	0	Scand1,H2bc12,H4c6,H2bc13,H4c1,Guca1a
```

## Interpretation

- FLT cluster 14 is the strongest non-obvious reorganization candidate:
  it has high latent shift and only about half of its genes map to its
  top GC counterpart. It still needs gene-level review because some
  high-shift genes are contractile or smooth-muscle associated.
- FLT clusters 13 and 15 are high-shift and DGEA-enriched, but their
  top genes are dominated by contractile/muscle markers, so they should
  be treated as residual composition or contamination-sensitive modules.
- FLT cluster 3 is mostly conserved with GC cluster 1, despite having many
  DGEA-overlapping genes. That pattern is more consistent with a stable
  module whose expression changes than with a strongly reorganized module.
- Cluster-level reorganization is not the same thing as flight-up
  expression. Direction still needs to be read from the DGEA columns and
  per-study consistency.

## Outputs

- `flt_to_gc_paired_cluster_summary.tsv`
- `gc_to_flt_paired_cluster_summary.tsv`
- `flt_gc_cluster_edges.tsv`
- `flt_displaced_modules.tsv`
- `paired_gene_level_table.tsv`
- `gene_lists/paired_reorganized_gene_lists.csv`
