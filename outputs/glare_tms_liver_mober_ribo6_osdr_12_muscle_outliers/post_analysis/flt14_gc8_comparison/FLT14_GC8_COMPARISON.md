# FLT14 vs GC8 Focused Comparison

This compares `FLT14_all`, `FLT14_not_GC8`, and `GC8_all` from the
12-filter MOBER + GLARE liver run. `FLT14_and_GC8` and
`GC8_not_FLT14` are included as controls so the split is explicit.

## Main Answer

- `GC8` is the primary GC match for `FLT14`, but only 163/303 FLT14 genes map there (53.8%).
- `FLT14_not_GC8` contains 140 genes (46.2% of FLT14).
- `GC8_all` contains 1,173 genes, but only 13.9% of GC8 comes from FLT14.
- Therefore GC8 is the closest ground counterpart, not a clean
  one-to-one liver equivalent.

## Module Summary

```tsv
module	gene_count	overlap_with_FLT14_genes	fraction_of_FLT14	overlap_with_GC8_genes	fraction_of_GC8	mean_latent_shift	significant_dgea_genes	significant_up_genes	significant_down_genes	top_significant_genes
FLT14_all	303	303	1	163	0.139	3.259	4	2	2	Apoa4,Mup15,Mup17,Cdkn1a
FLT14_and_GC8	163	163	0.538	163	0.139	3.419	3	1	2	Mup15,Mup17,Cdkn1a
FLT14_not_GC8	140	140	0.462	0	0	3.072	1	1	0	Apoa4
GC8_all	1173	163	0.538	1173	1	2.505	4	2	2	Mup15,Cyp4a14,Mup17,Cdkn1a
GC8_not_FLT14	1010	0	0	1010	0.861	2.357	1	1	0	Cyp4a14
```

## Where FLT14 Genes Go In GC

```tsv
gc_cluster	gene_count	fraction	mean_latent_shift	significant_dgea_genes	top_shift_genes	top_significant_genes
8	163	0.538	3.419	3	Scarna6,Slc25a4,Snora23,Snord22,Cyp3a44,Cyp3a16,Snord104,Cyp2a4,Cyp2b13,Cyp2b9,Ifi27l2a,Sult3a1	Mup15,Mup17,Cdkn1a
3	39	0.1287	2.649	0	Gpat3,Snord15b,Gbp11,Gbp10,Lsm7,Inhba,Slc13a3,Fmo4,Cd163,Sco2,Morf4l1-ps1,Hmgn5	
6	36	0.1188	3.673	0	Myh11,Sult2a1,Hao2,Sult2a2,Eno1b,Cyp2a22,Slc22a27,Ankrd55,Vldlr,Sult3a2,Abhd1,Syt1	
5	33	0.1089	3.346	0	Snord17,Neb,Snora74a,Rps2,Dio1,Rpl22l1,Prg4,Cops9,Lrtm1,Slc35d2,Tsku,Vnn1	
11	19	0.06271	2.273	0	Smlr1,Rpl27,Rarres1,Rplp2,Gsta4,Rpl35,Rps29,Tenm3,Atox1,Rpl28,Serpina11,G0s2	
9	6	0.0198	3.2	0	Ndufb4,Timm8b,Elovl6,Rpl39,Ndufs4,Rpl30	
4	4	0.0132	3.738	0	Fhl1,Serpinb6a,Fkbp11,F830016B08Rik	
12	3	0.009901	2.258	1	Serpina6,Apoa4,Prlr	Apoa4
```

## What GC8 Is Made Of

```tsv
flt_cluster	gene_count	fraction	mean_latent_shift	significant_dgea_genes	top_shift_genes	top_significant_genes
2	653	0.5567	2.176	0	Dpt,Ctla2a,Gt(ROSA)26Sor,Tagln,Cebpb,Tgtp1,Arrdc2,Mef2c,Mir99ahg,Tspyl4,Zbtb9,Stmn1	
9	239	0.2038	2.841	0	Rnu12,Snora44,Snora78,Snora7a,Gstp2,Snord118,Gvin1,Adgrf1,AI506816,Mcm6,Trim30d,Trim12a	
14	163	0.139	3.419	3	Scarna6,Slc25a4,Snora23,Snord22,Cyp3a44,Cyp3a16,Snord104,Cyp2a4,Cyp2b13,Cyp2b9,Ifi27l2a,Sult3a1	Mup15,Mup17,Cdkn1a
6	61	0.052	2.687	0	Vaultrc5,Med31,A730063M14Rik,Prss23,Gngt2,Arhgap10,Acot11,Cd52,Il1rn,Pirb,Thbd,Spata22	
1	23	0.01961	2.018	0	Calcrl,Chka,Lyz2,Ptprc,Fam135a,Per2,Slc16a5,F2r,Slc41a2,Ss18l2,Marcks,Eif3m	
8	13	0.01108	2.058	0	Snrpg,Plin4,Uqcc2,Avpi1,Timm10,Cox19,Sdf2l1,Pltp,Tspo,Rcan1,Mrps24,Atp6v0e	
5	11	0.009378	1.912	0	Eps8l2,Osgin1,Use1,Syvn1,H2-T22,Nlrp12,Ttc39c,C8b,H1f2,B4galt5,Crcp	
11	5	0.004263	1.943	1	Lpin1,Rpl13a,Cyp3a11,Cyp4a14,Egfr	Cyp4a14
12	3	0.002558	3.471	0	Synpo,1110038B12Rik,H1f3	
10	2	0.001705	1.365	0	Msmo1,Cyp51
```

## Set Overlaps

```tsv
left	right	left_genes	right_genes	intersection_genes	jaccard
GC8_all	GC8_not_FLT14	1173	1010	1010	0.861
FLT14_all	FLT14_and_GC8	303	163	163	0.538
FLT14_all	GC8_all	303	1173	163	0.1241
FLT14_and_GC8	GC8_all	163	1173	163	0.139
FLT14_all	FLT14_not_GC8	303	140	140	0.462
FLT14_all	GC8_not_FLT14	303	1010	0	0
FLT14_and_GC8	FLT14_not_GC8	163	140	0	0
FLT14_and_GC8	GC8_not_FLT14	163	1010	0	0
FLT14_not_GC8	GC8_all	140	1173	0	0
FLT14_not_GC8	GC8_not_FLT14	140	1010	0	0
```

## Interpretation

- `FLT14_and_GC8` is the matched portion of the module. It carries
  `Mup15`, `Mup17`, and `Cdkn1a` among the significant DGEA-overlap
  genes.
- `FLT14_not_GC8` is the split/reorganized portion. Its only significant
  DGEA-overlap gene in this analysis is `Apoa4`, but its high-shift gene
  list includes liver metabolism genes mixed with contractile-associated
  genes.
- `GC8_all` is broad. Most GC8 genes come from FLT2 and FLT9 rather than
  FLT14, so enrichment on all GC8 should be interpreted as the GC
  neighborhood around FLT14, not as FLT14's direct equivalent.

## Outputs

- `flt14_gc8_module_summary.tsv`
- `flt14_destination_gc_clusters.tsv`
- `gc8_source_flt_clusters.tsv`
- `flt14_gc8_set_overlaps.tsv`
- `flt14_gc8_gene_membership.tsv`
- `gene_lists/flt14_gc8_gene_lists.csv`
