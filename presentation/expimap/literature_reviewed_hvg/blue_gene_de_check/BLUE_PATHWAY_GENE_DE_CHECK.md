# Blue Pathway Gene-Level DE Check

This report checks whether the genes inside the literature-reviewed blue expiMap HVG pathways are themselves differentially expressed in the OSDR FLT vs GC query data.

Scope and method:

- Models checked: `liver_hvg`, `skin_hvg`, `thymus_hvg`, and `soleus_hvg`.
- Gene universe: only genes retained in each HVG expiMap query `.h5ad`; this is not an all-gene DGE pass.
- Gene test: log2(CPM + 1) FLT vs GC Welch test with BH FDR, plus an accession-level FLT-GC effect test where paired accessions exist.
- `n_any_fdr05` counts genes with pooled FDR < 0.05 or accession-level FDR < 0.05.
- Direction columns compare gene-level FLT-GC sign to the pathway-level expiMap FLT-GC direction.

## Dataset Summary

| model | n_samples | n_genes_hvg | n_flight | n_ground | n_accessions_with_both | accessions_with_both | n_pooled_gene_fdr05 | n_accession_gene_fdr05 | n_any_gene_fdr05 | n_nominal_gene_p05 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| liver_hvg | 231 | 1995 | 118 | 113 | 12 | OSD-137,OSD-164,OSD-168,OSD-173,OSD-242,OSD-245,OSD-379,OSD-457,OSD-463,OSD-47,OSD-48,OSD-686 | 109 | 1 | 110 | 587 |
| skin_hvg | 151 | 1997 | 80 | 71 | 6 | OSD-238,OSD-239,OSD-240,OSD-241,OSD-243,OSD-254 | 197 | 0 | 197 | 659 |
| thymus_hvg | 117 | 1994 | 63 | 54 | 5 | OSD-244,OSD-289,OSD-421,OSD-457,OSD-515 | 497 | 6 | 500 | 927 |
| soleus_hvg | 53 | 1975 | 28 | 25 | 3 | OSD-104,OSD-714,OSD-770 | 147 | 0 | 147 | 407 |

## Blue Pathway Gene DE Summary

| model | blue_pathways | with_any_hvg_gene_fdr05 | with_2plus_hvg_gene_fdr05 | with_any_nominal_gene_p05 | median_fraction_any_fdr05 | median_fraction_nominal_p05 | total_unique_sig_genes_in_blue_pathways |
| --- | --- | --- | --- | --- | --- | --- | --- |
| liver_hvg | 46 | 35 | 26 | 46 | 0.0435 | 0.306 | 61 |
| skin_hvg | 102 | 84 | 68 | 102 | 0.0688 | 0.3 | 151 |
| thymus_hvg | 86 | 85 | 81 | 86 | 0.215 | 0.41 | 340 |
| soleus_hvg | 44 | 44 | 36 | 44 | 0.0396 | 0.18 | 92 |

## Direction Summary

| model | blue_pathways | any_fdr_gene | two_plus_fdr_genes | strong_gene_support | moderate_gene_support | weak_or_no_gene_support | zero_fdr_gene | sig_gene_mostly_same_direction | sig_gene_mostly_opposite_direction | sig_gene_mixed_tie | median_hvg_same_direction | total_unique_sig_genes_in_blue_pathways |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| liver_hvg | 46 | 35 | 26 | 12 | 14 | 20 | 11 | 16 | 18 | 1 | 23 | 61 |
| skin_hvg | 102 | 84 | 68 | 38 | 30 | 34 | 18 | 43 | 35 | 6 | 58 | 151 |
| soleus_hvg | 44 | 44 | 36 | 14 | 22 | 8 | 0 | 20 | 17 | 7 | 28 | 92 |
| thymus_hvg | 86 | 85 | 81 | 70 | 11 | 5 | 1 | 41 | 41 | 3 | 46 | 340 |

## Interpretation

- Liver: most blue pathways contain at least one DE HVG gene, but the gene-level direction is mixed. Stronger gene-level support is in transport, MHC/TCR/interferon, Rho/cytoskeleton, fatty-acid, and vitamin/cofactor metabolism. Several ECM, PDGF, MET, heat-stress, porphyrin, and insulin terms have weak or no FDR-level member genes and should remain hypothesis-level.
- Skin: many blue pathways have DE member genes, especially vascular/GPCR/RTK, Rho/cytoskeleton, immune, lipid, carbohydrate, and SLC terms. However, accession-level gene FDR did not pass for individual genes, so these are pooled-sample signals and should be interpreted as screen-level support.
- Thymus: this has the strongest gene-level backing. Almost every blue pathway has FDR-significant member genes, and the signal is dense in immune system, neutrophil/degranulation, ECM, RTK/EGFR/TGFA, cell-cycle/protein-modification, and trafficking pathways. Direction is split, which means some blue pathways are capturing immune remodeling rather than simple whole-pathway activation/suppression.
- Soleus: every blue pathway has at least one FDR-significant HVG gene. The most convincing gene-level support is in metabolism, biological oxidations, small-molecule transport, calcium/conduction-like muscle genes, and developmental/neuromuscular terms. Rho-family terms have many DE genes but often in the opposite sign from the pathway-level shift, so those should be treated cautiously.

## Strongest Blue Pathway Examples

### liver_hvg
| term | observed_direction | gene_level_support_strength | fdr_gene_direction_call | n_pathway_genes_measured_hvg | n_any_fdr05 | n_nominal_p05 | fraction_any_fdr05 | top_gene_examples |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R-MMU-382551_TRANSPORT_OF_SMALL_MOLECULES | FLT_up | strong | mostly_opposite_direction | 176 | 13 | 47 | 0.0739 | Apoa4:+1.04,q=0.000497; Slc34a2:+0.72,q=0.000934; Slc5a6:-0.44,q=0.00123; Slc38a3:-0.39,q=0.0134; ENSMUSG00000026614:+0.44,q=0.0134; Slc43a1:-0.22,q=0.0177; Aqp9:-0.22,q=0.023; Nedd4l:+0.32,q=0.023 |
| R-MMU-2132295_MHC_CLASS_II_ANTIGEN_PRESENTATION | FLT_down | strong | mostly_same_direction | 33 | 11 | 16 | 0.333 | Tubb2a:+0.92,q=6.29e-05; H2-Eb1:-0.51,q=0.00263; Tubb4b:+0.44,q=0.00727; H2-DMa:-0.35,q=0.0134; Cd74:-0.41,q=0.0144; Tuba1c:+0.38,q=0.0227; Tuba8:+0.52,q=0.0253; H2-Ab1:-0.40,q=0.0289 |
| R-MMU-194315_SIGNALING_BY_RHO_GTPASES | FLT_down | strong | mostly_opposite_direction | 173 | 10 | 41 | 0.0578 | Tubb2a:+0.92,q=6.29e-05; Tubb4b:+0.44,q=0.00727; Net1:+0.27,q=0.0177; Tuba1c:+0.38,q=0.0227; Cybb:-0.35,q=0.0235; Ncf1:-0.26,q=0.0235; Tuba8:+0.52,q=0.0253; Fgd2:-0.31,q=0.0398 |
| R-MMU-9716542_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3 | FLT_up | strong | mostly_same_direction | 174 | 10 | 42 | 0.0575 | Tubb2a:+0.92,q=6.29e-05; Tubb4b:+0.44,q=0.00727; Net1:+0.27,q=0.0177; Tuba1c:+0.38,q=0.0227; Cybb:-0.35,q=0.0235; Ncf1:-0.26,q=0.0235; Tuba8:+0.52,q=0.0253; Fgd2:-0.31,q=0.0398 |
| R-MMU-913531_INTERFERON_SIGNALING | FLT_down | strong | mostly_opposite_direction | 37 | 8 | 12 | 0.216 | Tubb2a:+0.92,q=6.29e-05; Tubb4b:+0.44,q=0.00727; Gbp5:-0.28,q=0.0188; Tuba1c:+0.38,q=0.0227; Tuba8:+0.52,q=0.0253; Gbp2:-0.37,q=0.0335; Tuba4a:+0.35,q=0.0398; Ifngr1:-0.17,q=0.0463 |
| R-MMU-195258_RHO_GTPASE_EFFECTORS | FLT_down | strong | mostly_opposite_direction | 70 | 8 | 18 | 0.114 | Tubb2a:+0.92,q=6.29e-05; Tubb4b:+0.44,q=0.00727; Tuba1c:+0.38,q=0.0227; Cybb:-0.35,q=0.0235; Ncf1:-0.26,q=0.0235; Tuba8:+0.52,q=0.0253; Tuba4a:+0.35,q=0.0398; Nckap1l:-0.28,q=0.0433 |
| R-MMU-425407_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT | FLT_up | strong | mostly_opposite_direction | 69 | 7 | 22 | 0.101 | Slc34a2:+0.72,q=0.000934; Slc5a6:-0.44,q=0.00123; Slc38a3:-0.39,q=0.0134; ENSMUSG00000026614:+0.44,q=0.0134; Slc43a1:-0.22,q=0.0177; Slc45a3:+0.46,q=0.028; Slc22a5:-0.33,q=0.0428; Slc2a2:-0.22,q=0.0509 |
| R-MMU-5358351_SIGNALING_BY_HEDGEHOG | FLT_up | strong | mostly_same_direction | 24 | 6 | 8 | 0.25 | Tubb2a:+0.92,q=6.29e-05; Tubb4b:+0.44,q=0.00727; Adcy7:-0.34,q=0.013; Tuba1c:+0.38,q=0.0227; Adcy1:+0.41,q=0.0253; Tuba4a:+0.35,q=0.0398; Ptch1:-0.13,q=0.233; Tubb2b:+0.06,q=0.281 |
| R-MMU-8978868_FATTY_ACID_METABOLISM | FLT_down | strong | mostly_same_direction | 71 | 5 | 21 | 0.0704 | Cyp4a14:+1.10,q=0.00123; Thrsp:-0.65,q=0.0376; Cpt1a:-0.29,q=0.0386; Slc22a5:-0.33,q=0.0428; ENSMUSG00000066072:+0.49,q=0.0439; Cyp1a2:-0.34,q=0.0541; Fads2:-0.26,q=0.0599; Cyp8b1:+0.41,q=0.0654 |
| R-MMU-9012999_RHO_GTPASE_CYCLE | FLT_up | strong | mostly_opposite_direction | 126 | 5 | 33 | 0.0397 | Net1:+0.27,q=0.0177; Cybb:-0.35,q=0.0235; Ncf1:-0.26,q=0.0235; Fgd2:-0.31,q=0.0398; Nckap1l:-0.28,q=0.0433; Ncf2:-0.27,q=0.0548; ENSMUSG00000048865:-0.23,q=0.0565; Vim:-0.23,q=0.0625 |
### skin_hvg
| term | observed_direction | gene_level_support_strength | fdr_gene_direction_call | n_pathway_genes_measured_hvg | n_any_fdr05 | n_nominal_p05 | fraction_any_fdr05 | top_gene_examples |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R-MMU-168256_IMMUNE_SYSTEM | FLT_down | strong | mostly_same_direction | 445 | 28 | 111 | 0.0629 | Itpr1:+0.30,q=0.00232; Dsc1:-0.43,q=0.0112; Ptprb:+0.35,q=0.0113; Nos2:+0.48,q=0.0123; Il22ra2:+0.65,q=0.0123; Il6st:+0.20,q=0.0123; Pecam1:+0.34,q=0.016; Dtx4:+0.25,q=0.0213 |
| R-MMU-9716542_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3 | FLT_down | strong | mostly_same_direction | 150 | 20 | 66 | 0.133 | Dlc1:+0.29,q=0.0123; Kntc1:-0.46,q=0.0272; Fgd5:+0.36,q=0.0283; Arhgap31:+0.15,q=0.0324; Akap12:+0.37,q=0.0345; Cenpe:-0.59,q=0.0377; Cenpf:-0.51,q=0.0401; Sgo1:-0.51,q=0.0435 |
| R-MMU-194315_SIGNALING_BY_RHO_GTPASES | FLT_down | strong | mostly_same_direction | 150 | 20 | 66 | 0.133 | Dlc1:+0.29,q=0.0123; Kntc1:-0.46,q=0.0272; Fgd5:+0.36,q=0.0283; Arhgap31:+0.15,q=0.0324; Akap12:+0.37,q=0.0345; Cenpe:-0.59,q=0.0377; Cenpf:-0.51,q=0.0401; Sgo1:-0.51,q=0.0435 |
| R-MMU-109582_HEMOSTASIS | FLT_down | strong | mixed_tie | 198 | 18 | 55 | 0.0909 | Itpr1:+0.30,q=0.00232; Vegfa:+0.34,q=0.0123; Nos2:+0.48,q=0.0123; Arrb1:+0.25,q=0.0123; Pecam1:+0.34,q=0.016; Kif4:-0.41,q=0.024; Slc7a11:-1.18,q=0.0252; Cenpe:-0.59,q=0.0377 |
| R-MMU-9012999_RHO_GTPASE_CYCLE | FLT_up | strong | mostly_opposite_direction | 108 | 13 | 40 | 0.12 | Dlc1:+0.29,q=0.0123; Fgd5:+0.36,q=0.0283; Arhgap31:+0.15,q=0.0324; Akap12:+0.37,q=0.0345; Arhgef2:+0.19,q=0.0441; Dsp:-0.19,q=0.0441; Arhgap11a:-0.26,q=0.0482; Gja1:-0.56,q=0.0482 |
| R-MMU-372790_SIGNALING_BY_GPCR | FLT_up | strong | mostly_same_direction | 155 | 13 | 44 | 0.0839 | Itpr1:+0.30,q=0.00232; Arrb1:+0.25,q=0.0123; Ccl5:+0.50,q=0.0126; Gabbr1:+0.20,q=0.0177; Pde4d:+0.34,q=0.0198; Mc5r:-0.43,q=0.0216; ENSMUSG00000049112:+0.48,q=0.0239; S1pr1:+0.35,q=0.0328 |
| R-MMU-425407_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT | FLT_up | strong | mostly_opposite_direction | 67 | 12 | 34 | 0.179 | Slc34a2:-0.78,q=0.0132; ENSMUSG00000021565:-0.68,q=0.0236; Slc7a11:-1.18,q=0.0252; Slc39a8:-0.61,q=0.0278; Slc30a1:-0.56,q=0.0392; Slc24a5:-1.04,q=0.0441; Slc16a7:-0.48,q=0.0441; Slc24a4:-1.01,q=0.047 |
| R-MMU-71387_METABOLISM_OF_CARBOHYDRATES_AND_CARBOHYDRATE_DERIVATIVES | FLT_down | strong | mixed_tie | 76 | 12 | 26 | 0.158 | Prelp:+0.45,q=0.0168; Slc37a2:-0.30,q=0.0198; Ppp1r3c:+0.64,q=0.0236; Fbp2:+0.74,q=0.0363; St3gal1:+0.26,q=0.04; Hmmr:-0.52,q=0.0402; Chst1:+0.44,q=0.0435; Cd44:-0.29,q=0.0435 |
| R-MMU-388396_GPCR_DOWNSTREAM_SIGNALLING | FLT_up | strong | mostly_same_direction | 142 | 12 | 42 | 0.0845 | Itpr1:+0.30,q=0.00232; Arrb1:+0.25,q=0.0123; Ccl5:+0.50,q=0.0126; Gabbr1:+0.20,q=0.0177; Pde4d:+0.34,q=0.0198; Mc5r:-0.43,q=0.0216; ENSMUSG00000049112:+0.48,q=0.0239; Adcy5:+0.45,q=0.0345 |
| R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES | FLT_up | strong | mostly_same_direction | 128 | 11 | 40 | 0.0859 | Flt1:+0.36,q=0.00492; Cdh5:+0.40,q=0.00791; Vegfa:+0.34,q=0.0123; Kdr:+0.40,q=0.0142; ENSMUSG00000062991:-0.63,q=0.0208; Tgfbr3:+0.25,q=0.0265; Flrt3:-0.41,q=0.038; Nrp1:+0.25,q=0.04 |
### soleus_hvg
| term | observed_direction | gene_level_support_strength | fdr_gene_direction_call | n_pathway_genes_measured_hvg | n_any_fdr05 | n_nominal_p05 | fraction_any_fdr05 | top_gene_examples |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R-MMU-1430728_METABOLISM | FLT_down | strong | mostly_same_direction | 396 | 47 | 95 | 0.119 | Fabp3:-1.46,q=4.79e-09; Ckmt2:-0.96,q=6.55e-09; Cidea:-1.36,q=4.07e-08; Ldhb:-1.32,q=5.92e-08; Idh2:-1.64,q=7.36e-08; Ephx1:-1.18,q=2.98e-07; Nqo1:-1.50,q=2.98e-07; ENSMUSG00000109372:-1.17,q=4.98e-07 |
| R-MMU-382551_TRANSPORT_OF_SMALL_MOLECULES | FLT_down | strong | mostly_same_direction | 173 | 19 | 52 | 0.11 | Pln:-2.14,q=6.55e-09; Slc38a4:+0.71,q=2.78e-06; Stom:-0.59,q=1.22e-05; Fxyd6:-1.20,q=1.48e-05; Lpl:-1.20,q=2.66e-05; Slc40a1:-0.92,q=2.84e-05; Car2:-0.75,q=3.81e-05; Mb:-0.82,q=0.000753 |
| R-MMU-211859_BIOLOGICAL_OXIDATIONS | FLT_down | strong | mostly_same_direction | 55 | 9 | 15 | 0.164 | Ephx1:-1.18,q=2.98e-07; Acss2:-0.95,q=1.44e-05; Acss1:-0.59,q=0.000477; Gsta2:-0.73,q=0.000979; Arnt2:-1.09,q=0.00335; Aldh1a1:+0.44,q=0.0133; Gsto1:-0.58,q=0.0133; Oplah:-1.00,q=0.0461 |
| R-MMU-1266738_DEVELOPMENTAL_BIOLOGY | FLT_up | strong | mostly_same_direction | 142 | 7 | 29 | 0.0493 | Myog:+1.48,q=2.99e-06; Lgi1:-1.30,q=8.64e-05; Enah:+1.28,q=0.000191; Gdnf:+0.87,q=0.000196; Mef2c:-0.34,q=0.000281; Myf6:+0.66,q=0.00323; Ank1:-0.35,q=0.01; Rap1gap:-0.63,q=0.0645 |
| R-MMU-194315_SIGNALING_BY_RHO_GTPASES | FLT_up | strong | mostly_opposite_direction | 162 | 7 | 26 | 0.0432 | Sgo2a:-0.99,q=2.98e-07; Cyfip2:-1.96,q=1.03e-06; Rhobtb1:-1.19,q=3.96e-06; Stom:-0.59,q=1.22e-05; Fgd4:-0.48,q=0.00133; Fgd3:-1.18,q=0.01; Rhou:-0.42,q=0.0128; Anln:+0.69,q=0.0557 |
| R-MMU-9716542_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3 | FLT_up | strong | mostly_opposite_direction | 164 | 7 | 27 | 0.0427 | Sgo2a:-0.99,q=2.98e-07; Cyfip2:-1.96,q=1.03e-06; Rhobtb1:-1.19,q=3.96e-06; Stom:-0.59,q=1.22e-05; Fgd4:-0.48,q=0.00133; Fgd3:-1.18,q=0.01; Rhou:-0.42,q=0.0128; Anln:+0.69,q=0.0557 |
| R-MMU-5576891_CARDIAC_CONDUCTION | FLT_up | strong | mostly_opposite_direction | 41 | 6 | 11 | 0.146 | Pln:-2.14,q=6.55e-09; Fxyd6:-1.20,q=1.48e-05; Kcnj2:-0.50,q=0.00189; Casq2:-0.69,q=0.00744; Atp2a2:-0.88,q=0.0199; Slc8a1:+0.54,q=0.0426; ENSMUSG00000057378:-0.64,q=0.107; Camk2a:+0.51,q=0.115 |
| R-MMU-112316_NEURONAL_SYSTEM | FLT_up | strong | mostly_opposite_direction | 67 | 6 | 18 | 0.0896 | Kcnq5:+0.52,q=0.000616; Kcnj2:-0.50,q=0.00189; Chrna4:-0.31,q=0.0318; ENSMUSG00000029205:+0.91,q=0.0379; Homer2:-0.61,q=0.0379; Abat:-0.77,q=0.0429; Chrna1:-0.64,q=0.0564; Slc1a1:+0.46,q=0.0589 |
| R-MMU-425407_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT | FLT_down | strong | mixed_tie | 68 | 6 | 19 | 0.0882 | Slc38a4:+0.71,q=2.78e-06; Slc40a1:-0.92,q=2.84e-05; Cp:-0.43,q=0.00391; Slc4a4:-0.60,q=0.00464; Slc8a1:+0.54,q=0.0426; ENSMUSG00000021565:+1.36,q=0.0431; Slc1a1:+0.46,q=0.0589; Slc16a3:+1.12,q=0.0608 |
| R-MMU-9012999_RHO_GTPASE_CYCLE | FLT_up | strong | mostly_opposite_direction | 118 | 6 | 18 | 0.0508 | Cyfip2:-1.96,q=1.03e-06; Rhobtb1:-1.19,q=3.96e-06; Stom:-0.59,q=1.22e-05; Fgd4:-0.48,q=0.00133; Fgd3:-1.18,q=0.01; Rhou:-0.42,q=0.0128; Anln:+0.69,q=0.0557; Fam13a:-0.52,q=0.09 |
### thymus_hvg
| term | observed_direction | gene_level_support_strength | fdr_gene_direction_call | n_pathway_genes_measured_hvg | n_any_fdr05 | n_nominal_p05 | fraction_any_fdr05 | top_gene_examples |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R-MMU-168256_IMMUNE_SYSTEM | FLT_up | strong | mostly_same_direction | 565 | 123 | 237 | 0.218 | Ube2c:-1.35,q=5.17e-05; App:+0.74,q=6.61e-05; Cdk1:-1.49,q=7.94e-05; Pglyrp1:+1.04,q=0.000133; Il6ra:+0.74,q=0.0002; Tnfrsf11b:+0.74,q=0.000212; Serpina3n:+1.23,q=0.000213; Pigr:+1.50,q=0.000218 |
| R-MMU-597592_POST_TRANSLATIONAL_PROTEIN_MODIFICATION | FLT_down | strong | mostly_opposite_direction | 254 | 75 | 129 | 0.295 | Smad7:+0.63,q=4.22e-05; Ube2c:-1.35,q=5.17e-05; App:+0.74,q=6.61e-05; Cdk1:-1.49,q=7.94e-05; Igfbp5:+0.88,q=0.000157; Top2a:-1.41,q=0.000281; Cdh2:+0.74,q=0.00032; Apoe:+0.90,q=0.00032 |
| R-MMU-9716542_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3 | FLT_down | strong | mostly_same_direction | 177 | 57 | 88 | 0.322 | Spc25:-1.47,q=1.91e-06; Jag1:+0.66,q=8.1e-06; Spc24:-1.47,q=8.33e-06; Sgo1:-1.50,q=3.91e-05; Prc1:-1.39,q=0.000151; Cenpa:-1.02,q=0.000157; Mad2l1:-1.04,q=0.000183; Bub1:-1.35,q=0.000191 |
| R-MMU-194315_SIGNALING_BY_RHO_GTPASES | FLT_down | strong | mostly_same_direction | 174 | 56 | 87 | 0.322 | Spc25:-1.47,q=1.91e-06; Jag1:+0.66,q=8.1e-06; Spc24:-1.47,q=8.33e-06; Sgo1:-1.50,q=3.91e-05; Prc1:-1.39,q=0.000151; Cenpa:-1.02,q=0.000157; Mad2l1:-1.04,q=0.000183; Bub1:-1.35,q=0.000191 |
| R-MMU-6798695_NEUTROPHIL_DEGRANULATION | FLT_up | strong | mostly_same_direction | 161 | 38 | 71 | 0.236 | Pglyrp1:+1.04,q=0.000133; Serpina3n:+1.23,q=0.000213; Pigr:+1.50,q=0.000218; ENSMUSG00000025701:+1.05,q=0.000968; Ctsz:+0.58,q=0.00104; Arsb:+0.57,q=0.00135; Ifi204:+0.61,q=0.0022; Plaur:+0.82,q=0.00248 |
| R-MMU-5653656_VESICLE_MEDIATED_TRANSPORT | FLT_down | strong | mostly_opposite_direction | 142 | 35 | 57 | 0.246 | App:+0.74,q=6.61e-05; Ap3s1:-0.70,q=0.000122; Cenpe:-1.27,q=0.000308; Apoe:+0.90,q=0.00032; ENSMUSG00000012443:-1.37,q=0.000958; ENSMUSG00000074899:+0.92,q=0.00098; Ctsz:+0.58,q=0.00104; Trf:+0.80,q=0.00114 |
| R-MMU-1474244_EXTRACELLULAR_MATRIX_ORGANIZATION | FLT_up | strong | mostly_same_direction | 115 | 34 | 59 | 0.296 | App:+0.74,q=6.61e-05; F11r:+0.48,q=0.00102; Itgb8:+0.82,q=0.00102; Ltbp3:+0.93,q=0.00155; Htra1:+0.68,q=0.0016; Sdc4:+0.33,q=0.00237; Mmp14:+0.62,q=0.00247; Adam15:+0.89,q=0.00332 |
| R-MMU-9012999_RHO_GTPASE_CYCLE | FLT_down | strong | mostly_opposite_direction | 128 | 33 | 59 | 0.258 | Jag1:+0.66,q=8.1e-06; Ect2:-1.15,q=0.00105; Diaph3:-1.13,q=0.00159; Pcdh7:+0.80,q=0.00165; Dsg2:+0.58,q=0.00282; Cyfip1:+0.41,q=0.00364; Arhgef12:+0.64,q=0.00415; Arhgef28:+0.82,q=0.00575 |
| R-MMU-195258_RHO_GTPASE_EFFECTORS | FLT_up | strong | mostly_opposite_direction | 63 | 31 | 37 | 0.492 | Spc25:-1.47,q=1.91e-06; Spc24:-1.47,q=8.33e-06; Sgo1:-1.50,q=3.91e-05; Prc1:-1.39,q=0.000151; Cenpa:-1.02,q=0.000157; Mad2l1:-1.04,q=0.000183; Bub1:-1.35,q=0.000191; Cenpe:-1.27,q=0.000308 |
| R-MMU-199991_MEMBRANE_TRAFFICKING | FLT_down | strong | mostly_opposite_direction | 91 | 30 | 46 | 0.33 | App:+0.74,q=6.61e-05; Ap3s1:-0.70,q=0.000122; Cenpe:-1.27,q=0.000308; ENSMUSG00000012443:-1.37,q=0.000958; ENSMUSG00000074899:+0.92,q=0.00098; Ctsz:+0.58,q=0.00104; Trf:+0.80,q=0.00114; Kif15:-1.24,q=0.00154 |

## Weakest / No-FDR Blue Pathway Examples

### liver_hvg
| term | observed_direction | gene_level_support_strength | n_pathway_genes_measured_hvg | n_any_fdr05 | n_nominal_p05 | median_hvg_direction_vs_pathway | top_gene_examples |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R-MMU-422356_REGULATION_OF_INSULIN_SECRETION | FLT_down | none | 15 | 0 | 2 | same | Slc2a2:-0.22,q=0.0509; Marcks:-0.21,q=0.22; Adra2a:+0.06,q=0.247; Gng2:-0.14,q=0.339; Itpr3:-0.15,q=0.348; Cacnb3:-0.08,q=0.446; Adcy6:-0.13,q=0.45; Rapgef4:-0.11,q=0.519 |
| R-MMU-189445_METABOLISM_OF_PORPHYRINS | FLT_up | none | 14 | 0 | 2 | opposite | Slco1b2:-0.24,q=0.128; Ugt1a5:-0.29,q=0.196; Ugt1a1:-0.15,q=0.235; Abcc2:-0.16,q=0.297; Abcc1:-0.17,q=0.312; Fabp1:-0.17,q=0.439; Cpox:+0.08,q=0.617; Alb:-0.04,q=0.642 |
| R-MMU-3371556_CELLULAR_RESPONSE_TO_HEAT_STRESS | FLT_up | none | 14 | 0 | 3 | same | Nup210:-0.18,q=0.297; Hspa1b:-0.35,q=0.324; Hspa1a:-0.23,q=0.366; Hspa5:+0.13,q=0.498; Bag3:+0.13,q=0.603; Cryab:+0.09,q=0.694; Hsph1:-0.11,q=0.715; Nupl1:+0.01,q=0.757 |
| R-MMU-9955298_SLC_MEDIATED_TRANSPORT_OF_ORGANIC_ANIONS | FLT_down | none | 16 | 0 | 5 | same | Slco1b2:-0.24,q=0.128; Slc16a1:+0.21,q=0.179; ENSMUSG00000020805:-0.19,q=0.25; Slco3a1:-0.10,q=0.329; Slc16a7:-0.14,q=0.418; Slc16a2:-0.13,q=0.431; Slc44a3:-0.07,q=0.441; Emb:-0.10,q=0.445 |
| R-MMU-373755_SEMAPHORIN_INTERACTIONS | FLT_up | none | 24 | 0 | 7 | opposite | Ptprc:-0.40,q=0.0548; Fyn:-0.16,q=0.13; Nrp1:-0.22,q=0.133; Sema6d:-0.16,q=0.167; Cd72:-0.19,q=0.177; Tyrobp:-0.28,q=0.227; Pak1:-0.14,q=0.239; Sema4d:-0.17,q=0.28 |
| R-MMU-6806834_SIGNALING_BY_MET | FLT_up | none | 19 | 0 | 7 | opposite | Col3a1:-0.37,q=0.0532; Col5a2:-0.26,q=0.0906; Col5a1:-0.21,q=0.227; Col5a3:+0.35,q=0.244; Sh3kbp1:-0.11,q=0.285; Src:-0.14,q=0.309; Met:-0.15,q=0.32; Tns3:-0.15,q=0.334 |
| R-MMU-382556_ABC_FAMILY_PROTEIN_MEDIATED_TRANSPORT | FLT_down | none | 21 | 0 | 8 | same | Abcb1a:-0.39,q=0.0933; Abcc4:-0.25,q=0.13; Abcg1:-0.22,q=0.144; Abca9:-0.18,q=0.193; Abcg5:-0.19,q=0.196; Abcg8:-0.21,q=0.22; Abca6:-0.18,q=0.293; Abcc2:-0.16,q=0.297 |
| R-MMU-2173789_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS | FLT_up | none | 21 | 0 | 9 | opposite | Tgfbr3:-0.22,q=0.164; Itga8:-0.32,q=0.167; Ltbp3:-0.18,q=0.238; Ltbp1:-0.17,q=0.274; Tgfb3:+0.16,q=0.282; Tgfb2:+0.08,q=0.36; Ltbp4:-0.22,q=0.36; Ltbp2:-0.13,q=0.394 |
### skin_hvg
| term | observed_direction | gene_level_support_strength | n_pathway_genes_measured_hvg | n_any_fdr05 | n_nominal_p05 | median_hvg_direction_vs_pathway | top_gene_examples |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R-MMU-191273_CHOLESTEROL_BIOSYNTHESIS | FLT_down | none | 13 | 0 | 1 | same | Sqle:-0.20,q=0.09; Fdps:-0.11,q=0.4; Cyp51:-0.10,q=0.448; Msmo1:-0.11,q=0.453; Fdft1:-0.06,q=0.5; ENSMUSG00000021670:-0.08,q=0.621; Hmgcs1:-0.05,q=0.677; Acat2:-0.05,q=0.698 |
| R-MMU-166658_COMPLEMENT_CASCADE | FLT_down | none | 19 | 0 | 2 | opposite | C3:+0.30,q=0.129; Hc:+0.51,q=0.155; C4b:+0.30,q=0.199; Fcna:+0.33,q=0.207; Clu:-0.25,q=0.264; C5ar1:+0.21,q=0.286; Serping1:+0.20,q=0.314; C1ra:+0.16,q=0.335 |
| R-MMU-1236975_ANTIGEN_PROCESSING_CROSS_PRESENTATION | FLT_up | none | 16 | 0 | 2 | same | Cd36:+0.35,q=0.0974; H2-Q10:+0.23,q=0.174; Mrc2:+0.27,q=0.244; H2-T22:+0.13,q=0.28; Psmb8:+0.17,q=0.343; H2-Q7:+0.52,q=0.348; Fcgr1:+0.19,q=0.356; H2-M5:+0.16,q=0.463 |
| R-MMU-977606_REGULATION_OF_COMPLEMENT_CASCADE | FLT_down | none | 16 | 0 | 2 | opposite | C3:+0.30,q=0.129; Hc:+0.51,q=0.155; C4b:+0.30,q=0.199; Clu:-0.25,q=0.264; C5ar1:+0.21,q=0.286; Serping1:+0.20,q=0.314; C1ra:+0.16,q=0.335; C1qb:+0.21,q=0.356 |
| R-MMU-8878171_TRANSCRIPTIONAL_REGULATION_BY_RUNX1 | FLT_up | none | 44 | 0 | 3 | opposite | ENSMUSG00000002028:+0.14,q=0.0622; Ccnd1:-0.23,q=0.1; Esr1:+0.24,q=0.146; ENSMUSG00000029673:+0.13,q=0.167; Ccnd2:-0.11,q=0.209; Ctsl:-0.13,q=0.236; H3c2:-0.59,q=0.265; H2bc12:-0.60,q=0.292 |
| R-MMU-186797_SIGNALING_BY_PDGF | FLT_up | none | 26 | 0 | 3 | same | Pdgfc:+0.21,q=0.101; Pdgfrb:+0.28,q=0.106; Col2a1:-0.16,q=0.141; Pik3r1:+0.12,q=0.246; Col9a2:+0.20,q=0.259; Col6a3:-0.13,q=0.378; Thbs4:-0.16,q=0.429; Thbs3:+0.13,q=0.434 |
| R-MMU-9909648_REGULATION_OF_PD_L1_CD274_EXPRESSION | FLT_up | none | 35 | 0 | 4 | opposite | Nek2:-0.45,q=0.0566; Ccnd1:-0.23,q=0.1; ENSMUSG00000028518:+0.36,q=0.127; H3c2:-0.59,q=0.265; H2bc12:-0.60,q=0.292; H2bc13:-0.61,q=0.301; Cd274:+0.14,q=0.301; ENSMUSG00000094777:-0.56,q=0.304 |
| R-MMU-9769739_REGULATION_OF_CLOTTING_CASCADE | FLT_up | none | 23 | 0 | 4 | same | Sdc4:+0.12,q=0.0553; Serpine2:+0.22,q=0.0655; Procr:-0.37,q=0.151; Vwf:+0.29,q=0.155; Ano5:+0.47,q=0.231; F5:-0.28,q=0.248; Sdc1:-0.21,q=0.256; Serpine1:-0.41,q=0.301 |
### soleus_hvg
| term | observed_direction | gene_level_support_strength | n_pathway_genes_measured_hvg | n_any_fdr05 | n_nominal_p05 | median_hvg_direction_vs_pathway | top_gene_examples |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R-MMU-5668541_TNFR2_NON_CANONICAL_NF_KB_PATHWAY | FLT_up | weak | 16 | 1 | 2 | same | Eda2r:+1.17,q=7.6e-05; Tnfrsf12a:+0.95,q=0.127; Tnfrsf13b:+0.36,q=0.508; Tnfrsf11b:-0.28,q=0.664; Map3k14:+0.37,q=0.672; Relb:+0.21,q=0.902; Tnfrsf13c:-0.03,q=0.921; Ltb:+0.04,q=0.93 |
| R-MMU-216083_INTEGRIN_CELL_SURFACE_INTERACTIONS | FLT_down | weak | 49 | 1 | 6 | opposite | Col4a5:-0.79,q=0.0126; Col7a1:+0.84,q=0.172; Col9a1:+0.22,q=0.194; Itga5:+0.70,q=0.254; Comp:-1.14,q=0.278; Jam2:-0.42,q=0.285; Col8a2:-1.00,q=0.317; Itga11:-0.60,q=0.359 |
| R-MMU-449147_SIGNALING_BY_INTERLEUKINS | FLT_down | weak | 67 | 1 | 7 | opposite | Map2k6:-0.91,q=0.000911; Jak2:+0.27,q=0.0506; H3c8:+0.85,q=0.0605; Lifr:-0.53,q=0.0708; Il10:+0.15,q=0.11; H3c10:+1.05,q=0.198; Fos:+0.63,q=0.283; Jun:+0.45,q=0.314 |
| R-MMU-9013149_RAC1_GTPASE_CYCLE | FLT_up | weak | 63 | 1 | 7 | same | Cyfip2:-1.96,q=1.03e-06; Fam13a:-0.52,q=0.09; Pak1:+0.49,q=0.0983; Abr:-0.62,q=0.185; Ect2:+0.52,q=0.219; Tiam2:+0.60,q=0.221; Tfrc:-0.70,q=0.241; Iqgap3:+0.53,q=0.266 |
| R-MMU-375276_PEPTIDE_LIGAND_BINDING_RECEPTORS | FLT_up | weak | 63 | 1 | 8 | same | Mchr1:-0.80,q=0.0261; F2:-0.53,q=0.0692; Ccr3:+0.43,q=0.0816; Ackr3:-0.46,q=0.0886; Fpr2:+0.45,q=0.125; Cxcl11:-0.75,q=0.14; Agtr2:-0.10,q=0.145; Avpr2:+0.04,q=0.146 |
| R-MMU-198933_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL | FLT_up | weak | 41 | 1 | 8 | same | Cd1d1:-0.89,q=0.00828; Treml4:+0.21,q=0.146; Cd200:-0.61,q=0.16; ENSMUSG00000069609:+0.66,q=0.178; Cd300lg:-0.60,q=0.227; Cd300ld3:+0.65,q=0.239; ENSMUSG00000071068:+0.36,q=0.25; Cd300lf:+0.54,q=0.251 |
| R-MMU-8856828_CLATHRIN_MEDIATED_ENDOCYTOSIS | FLT_up | weak | 31 | 1 | 8 | same | Hbegf:+1.58,q=0.0247; Areg:+0.29,q=0.0528; Ldlr:+1.05,q=0.0602; Synj2:+0.36,q=0.11; Adrb2:+0.39,q=0.141; Avpr2:+0.04,q=0.146; Fcho1:+0.33,q=0.209; Tfrc:-0.70,q=0.241 |
| R-MMU-1650814_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES | FLT_up | weak | 42 | 1 | 9 | opposite | Col4a5:-0.79,q=0.0126; ENSMUSG00000028197:-0.47,q=0.113; Col7a1:+0.84,q=0.172; Col11a1:-1.06,q=0.187; Col9a1:+0.22,q=0.194; Col22a1:-0.94,q=0.221; Col11a2:-1.27,q=0.229; Col8a2:-1.00,q=0.317 |
### thymus_hvg
| term | observed_direction | gene_level_support_strength | n_pathway_genes_measured_hvg | n_any_fdr05 | n_nominal_p05 | median_hvg_direction_vs_pathway | top_gene_examples |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R-MMU-9007101_RAB_REGULATION_OF_TRAFFICKING | FLT_down | none | 14 | 0 | 3 | opposite | Sbf2:+0.48,q=0.0598; Dennd5a:+0.41,q=0.0963; Rin3:+0.41,q=0.176; Rab7b:+0.31,q=0.198; Dennd2d:+0.28,q=0.237; Dennd3:+0.40,q=0.32; Dennd4a:+0.26,q=0.449; Dennd4c:-0.19,q=0.519 |
| R-MMU-5621481_C_TYPE_LECTIN_RECEPTORS_CLRS | FLT_up | weak | 19 | 1 | 3 | same | Relb:+0.51,q=0.0473; Lyn:+0.46,q=0.0673; Nfkb2:+0.51,q=0.097; Clec4e:+0.44,q=0.148; Src:+0.66,q=0.166; Plcg2:+0.54,q=0.166; Fcer1g:+0.25,q=0.174; Nfatc2:+0.49,q=0.217 |
| R-MMU-1428517_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT | FLT_down | weak | 24 | 1 | 4 | opposite | Rmnd5a:-0.63,q=0.0206; ENSMUSG00000064368:+0.34,q=0.23; Ubc:+0.23,q=0.479; Me1:+0.29,q=0.487; Ldhb:-0.12,q=0.571; Idh2:+0.13,q=0.591; ENSMUSG00000064357:-0.29,q=0.594; ENSMUSG00000064356:+0.17,q=0.656 |
| R-MMU-9851695_EPIGENETIC_REGULATION_OF_ADIPOGENESIS_GENES_BY_MLL3_AND_MLL4_COMPLEXES | FLT_up | weak | 40 | 1 | 13 | opposite | H2ac23:-1.43,q=0.238; H3c14:-0.95,q=0.0525; H2ax:-0.90,q=0.0606; H2ac11:-1.27,q=0.0677; H2bu2:-0.86,q=0.116; H2ac4:-1.06,q=0.126; H2bc3:-1.17,q=0.134; ENSMUSG00000069303:-1.02,q=0.144 |
| R-MMU-9940951_INTERACTION_OF_NURD_COMPLEXES_WITH_TRANSCRIPTION_FACTORS | FLT_up | weak | 45 | 1 | 15 | opposite | H2ac23:-1.43,q=0.238; H3c14:-0.95,q=0.0525; H2ax:-0.90,q=0.0606; Zfp827:+0.57,q=0.0606; H2ac11:-1.27,q=0.0677; Cdk2ap1:-0.32,q=0.0761; H2bu2:-0.86,q=0.116; H2ac4:-1.06,q=0.126 |
| R-MMU-2454202_FC_EPSILON_RECEPTOR_FCERI_SIGNALING | FLT_up | moderate | 52 | 2 | 5 | same | Jun:+0.48,q=0.00592; Lat2:+0.55,q=0.0226; ENSMUSG00000087642:+0.68,q=0.0606; ENSMUSG00000096108:+0.90,q=0.0662; Lyn:+0.46,q=0.0673; Plcg2:+0.54,q=0.166; Fcer1g:+0.25,q=0.174; Lat:-0.44,q=0.189 |
| R-MMU-2172127_DAP12_INTERACTIONS | FLT_down | moderate | 16 | 2 | 6 | opposite | Cd300lb:+0.70,q=0.0079; Klrc1:+0.58,q=0.0339; Tyrobp:+0.49,q=0.0525; B2m:+0.23,q=0.0722; Sirpa:+0.46,q=0.0867; Clec5a:+0.48,q=0.113; Lck:-0.50,q=0.134; Klrc2:+0.44,q=0.166 |
| R-MMU-1296071_POTASSIUM_CHANNELS | FLT_up | moderate | 19 | 2 | 7 | same | Kcnq4:+0.83,q=0.00597; Kcnc1:+0.70,q=0.0201; Gnb3:+0.49,q=0.0513; Kcnj15:+0.47,q=0.0595; Kcnj8:+0.57,q=0.0671; Gng2:+0.35,q=0.0857; ENSMUSG00000040136:+0.54,q=0.0955; Kcnf1:+0.32,q=0.165 |

## Output Files

- `blue_pathway_gene_de_overlap.tsv`: full blue pathway by gene-DE overlap table.
- `blue_pathway_significant_genes.tsv`: FDR-significant member genes per blue pathway.
- `blue_pathway_gene_de_manual_calls.tsv`: compact manual-review calls per blue pathway.
- `blue_pathway_gene_de_direction_summary.tsv`: tissue-level direction/support summary.
