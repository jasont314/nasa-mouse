# Per-Study GLARE vs DGEA

Each study was analyzed separately. DESeq2 inputs were exported from the
same per-study GLARE target matrices, so this comparison uses the 12
muscle-outlier filter and the same ERCC/noERCC selections as the GLARE runs.

## Inputs

- Counts: `outputs/glare_per_study_liver_noercc_12filter/dgea_comparison/deseq2_inputs/counts.tsv`
- Metadata: `outputs/glare_per_study_liver_noercc_12filter/dgea_comparison/deseq2_inputs/sample_metadata.tsv`
- Gene symbols: `outputs/glare_per_study_liver_noercc_12filter/dgea_comparison/deseq2_inputs/gene_symbols.tsv`
- DESeq2 alpha: 0.05
- Cluster DEG strict flag: adjusted p < 0.05 and abs(log2FC) >= 1.0

## Study Summary

| accession | n_flight | n_ground | deseq2_design | genes_tested_dgea | significant_padj05 | significant_padj05_abs_lfc | deg_enriched_glare_clusters_fdr05 | abs_stat_vs_latent_shift_spearman_rho | latent_shift_deg_roc_auc | rank_pathway_glare_overlaps |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| OSD-379 | 31 | 31 | ~ stratum + condition | 15381 | 1622 | 121 | 10 | 0.0782 | 0.562 | 148 |
| OSD-245 | 18 | 19 | ~ stratum + condition | 15798 | 700 | 54 | 8 | 0.023 | 0.598 | 61 |
| OSD-463 | 12 | 10 | ~ condition | 15362 | 2704 | 472 | 11 | 0.208 | 0.644 | 81 |
| OSD-168 | 9 | 9 | ~ condition | 16339 | 6 | 3 | 0 | -0.0307 | 0.86 | 227 |
| OSD-48 | 7 | 7 | ~ condition | 13691 | 82 | 42 | 1 | 0.0358 | 0.713 | 62 |
| OSD-137 | 5 | 6 | ~ condition | 14625 | 1 | 0 | 0 | -0.0668 | 0.412 | 101 |

## DEG-Enriched GLARE Clusters

| accession | location | cluster | tested_dgea_genes | significant_padj05_genes | significant_padj05_abs_lfc_genes | fisher_odds_ratio | fisher_fdr_bh | top_significant_genes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| OSD-463 | FLT | 6 | 806 | 240 | 35 | 2.08 | 2.93e-17 | Cyp4a14,Cyp8b1,Apoa4,Aco2,Cdo1,Mat1a,Saa4,Phb2,Abcb11,Gpt |
| OSD-463 | GC | 4 | 2087 | 510 | 72 | 1.63 | 1.5e-16 | Cyp8b1,Saa2,Saa1,Abca8a,Apoa4,Slc17a2,Bhmt,Aco2,Hectd2os,St6gal1 |
| OSD-245 | GC | 6 | 1311 | 119 | 4 | 2.39 | 2.41e-13 | Nedd4l,Cldn1,Ss18l2,Dnaaf9,Vamp3,Dclre1a,Ncoa1,Mup12,Bptf,Cgrrf1 |
| OSD-463 | FLT | 1 | 3043 | 669 | 43 | 1.42 | 1.78e-11 | Rogdi,Cntnap1,Cirbp,Cyp4a31,Inca1,Get4,Bcap29,Scara5,Elac1,Iqsec2 |
| OSD-463 | GC | 9 | 781 | 207 | 38 | 1.75 | 8.07e-10 | Slco1a4,Car14,Abcd1,Dpy19l3,Aqp4,Htatip2,H2-Eb1,Slc22a28,Pogk,Tstd3 |
| OSD-245 | FLT | 3 | 2301 | 162 | 6 | 1.82 | 5.94e-09 | Nedd4l,Chka,Sdr9c7,Trim2,Net1,Syne1,Cldn1,Dnaaf9,Rnf144b,Clock |
| OSD-245 | FLT | 1 | 2958 | 191 | 11 | 1.67 | 6.35e-08 | Fam47e,Dtx4,Pdk4,Tefm,Ss18l2,Rcor1,Vamp3,Rai14,Dclre1a,Ncoa1 |
| OSD-379 | FLT | 3 | 2810 | 381 | 22 | 1.43 | 1.75e-07 | Slc22a5,Inhba,Leo1,Cars,Tmie,Gm11437,Arhgap24,Pde9a,Cd9,Arntl |
| OSD-463 | GC | 2 | 2794 | 582 | 26 | 1.3 | 2.72e-06 | Rogdi,Cyp4a32,Inca1,Get4,Bcap29,Elac1,Acot3,Iqsec2,Col13a1,Taf1b |
| OSD-379 | FLT | 6 | 1333 | 196 | 3 | 1.53 | 3.01e-06 | Svil,Pter,Tars,Gpcpd1,Tubb2a,Ppp1r3c,St3gal1,Acaa1a,Tab2,Ern1 |
| OSD-245 | GC | 7 | 1282 | 94 | 9 | 1.82 | 4.86e-06 | Gclc,Insig2,Cyp2a5,Mup10,Wfdc17,Gch1,Cyp2a4,Rnase4,Tef,Gys2 |
| OSD-379 | GC | 5 | 1984 | 273 | 2 | 1.43 | 1.2e-05 | Svil,Pter,Tars,Gpcpd1,1810008I18Rik,Tubb2a,Ppp1r3c,Dsp,4932438A13Rik,Fam126b |
| OSD-463 | GC | 6 | 26 | 14 | 14 | 5.48 | 9.92e-05 | Sult3a1,Atp13a5,Krt79,Scg2,Cyp3a57,Shox2,Patl2,A530016L24Rik,Sdcbp2,Fbxw15 |
| OSD-463 | FLT | 3 | 1837 | 388 | 61 | 1.3 | 9.98e-05 | Smyd4,Flvcr2,Kdf1,Trim80,1810010H24Rik,Ppp3cc,E330011O21Rik,A530088E08Rik,Srrm4os,Mettl15 |
| OSD-245 | GC | 2 | 2429 | 146 | 3 | 1.48 | 0.000216 | Chka,Sdr9c7,Trim2,Tefm,Rcor1,Rai14,Gmnn,Slc5a6,Cebpzos,Tfcp2l1 |
| OSD-463 | FLT | 5 | 998 | 222 | 25 | 1.37 | 0.000235 | AI463170,Treh,Bcl6,Cyp2c39,Serpinh1,Dpy19l3,Oxld1,Taf1b,Cpt1b,H2-Eb1 |
| OSD-379 | GC | 4 | 2182 | 282 | 3 | 1.31 | 0.000403 | Inhba,Cars,Tmie,Pde9a,Trim2,Cd9,Erbb4,Chka,Pdcl3,Ston2 |
| OSD-379 | GC | 2 | 225 | 43 | 41 | 2.03 | 0.000403 | Tnnt3,Tnni2,Sln,Mypn,Mybpc1,Tnnc2,Xirp2,Mybpc2,Mb,Tbx15 |
| OSD-245 | FLT | 8 | 1032 | 70 | 8 | 1.63 | 0.00108 | Clpx,Gclc,Insig2,Cyp2a5,Mup10,Cyp2a4,Rnase4,Tef,Wfdc21,Sdhd |
| OSD-245 | GC | 9 | 905 | 60 | 4 | 1.58 | 0.00396 | Clpx,Clock,Zfp652,Spcs1,Nhlrc2,Pdcd6,Tomm22,Ubqln1,Fgd6,Aox1 |
| OSD-379 | FLT | 5 | 1578 | 203 | 9 | 1.29 | 0.0044 | Avpr1a,Bach1,Asb13,Usp6nl,Irf6,Ola1,Trim2,Fam126b,Rdh11,Pmvk |
| OSD-379 | FLT | 8 | 874 | 120 | 2 | 1.38 | 0.0044 | 1810008I18Rik,Upp2,Dsp,Sc5d,4932438A13Rik,Slc25a47,Slc7a2,Cyp2c70,Cpt1a,Abcb11 |
| OSD-463 | FLT | 4 | 1633 | 330 | 28 | 1.21 | 0.00654 | Arsg,Rbm3,Atf5,Slco1a4,Abca8a,Slc17a2,Tmem120a,Glrx,Hs3st3b1,Steap4 |
| OSD-379 | FLT | 11 | 100 | 20 | 14 | 2.13 | 0.00924 | Npas2,Ankrd23,Hspb1,Cox8b,Nrap,Rfx4,Pgam2,Srl,Wee1,Ppard |

## Recurring DGEA/GLARE Pathway Concordance

| direction | clean_term | study_count | accessions | best_dgea_fdr_bh | best_glare_fdr_bh | mean_wald_stat_shift_mean | example_glare_clusters |
| --- | --- | --- | --- | --- | --- | --- | --- |
| down_in_flight | Peptide Chain Elongation | 4 | OSD-168,OSD-245,OSD-379,OSD-48 | 5.28e-22 | 5.44e-59 | -1.19 | GC7;FLT8;FLT10;GC3;GC9;FLT8;GC9;GC5;FLT6;GC8;FLT7;FLT4;GC5 |
| down_in_flight | Influenza Viral Rna Transcription And Replication | 4 | OSD-168,OSD-245,OSD-379,OSD-48 | 5.28e-22 | 4.35e-51 | -1.02 | FLT10;GC3;GC9;GC7;FLT8;FLT8;GC9;GC5;FLT6;GC8;FLT7;FLT4;GC5 |
| down_in_flight | 3 Utr Mediated Translational Regulation | 4 | OSD-168,OSD-245,OSD-379,OSD-48 | 4.46e-17 | 5.44e-59 | -0.903 | GC7;FLT8;FLT10;GC3;FLT7;FLT8;GC9;GC5;FLT6;GC8;FLT4;FLT7;GC5 |
| down_in_flight | Platelet Activation Signaling And Aggregation | 4 | OSD-137,OSD-168,OSD-245,OSD-463 | 6.37e-11 | 2.31e-13 | -0.437 | GC4;FLT6;FLT4;GC10;FLT3;GC9;FLT8;GC8;FLT7;GC8;GC5;FLT4;GC6;FLT8;GC7;GC5;FLT10;FLT4 |
| down_in_flight | Response To Elevated Platelet Cytosolic Ca2  | 4 | OSD-137,OSD-168,OSD-245,OSD-463 | 1.68e-08 | 4.11e-18 | -0.622 | GC4;FLT6;FLT4;GC10;FLT3;FLT8;GC6;GC8;FLT7;GC8;FLT4;GC5;FLT8;GC7;GC5 |
| down_in_flight | Formation Of Fibrin Clot Clotting Cascade | 4 | OSD-137,OSD-168,OSD-245,OSD-463 | 1.49e-06 | 3.8e-12 | -1.06 | GC10;FLT3;FLT6;GC4;FLT7;GC8;FLT8;GC7 |
| up_in_flight | Scfskp2 Mediated Degradation Of P27 P21 | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 1.62e-06 | 1.94e-18 | 0.703 | GC3;FLT7;FLT5;FLT10;GC5;FLT6;GC4;FLT4;FLT6;FLT8;GC8;GC6;FLT3 |
| up_in_flight | Orc1 Removal From Chromatin | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 1.9e-06 | 6.61e-15 | 0.663 | FLT7;GC3;FLT5;GC5;FLT6;GC4;FLT4;FLT6;FLT8;GC8;GC6;FLT3 |
| up_in_flight | Vif Mediated Degradation Of Apobec3G | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 2.22e-06 | 5.81e-22 | 0.813 | GC3;FLT7;FLT5;FLT10;GC5;FLT6;FLT5;GC6;GC4;FLT4;FLT6;GC7;FLT8;GC6;FLT3;GC8;GC7 |
| up_in_flight | Scf Beta Trcp Mediated Degradation Of Emi1 | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 2.22e-06 | 9.95e-22 | 0.812 | GC3;FLT7;FLT5;FLT10;GC5;FLT6;FLT5;GC4;FLT4;FLT6;FLT8;GC6;GC8;FLT3 |
| up_in_flight | Regulation Of Apoptosis | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 2.22e-06 | 8.06e-20 | 0.758 | GC3;FLT7;FLT5;GC5;FLT6;FLT5;FLT4;GC4;FLT6;GC7;FLT8;GC6;GC8;FLT3;GC7 |
| up_in_flight | Cross Presentation Of Soluble Exogenous Antigens Endosomes | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 2.22e-06 | 1.33e-19 | 0.739 | GC5;FLT6;FLT5;FLT7;GC3;FLT5;FLT4;GC4;FLT8;GC8;FLT3;GC6;GC7 |
| up_in_flight | Assembly Of The Pre Replicative Complex | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 2.22e-06 | 1.82e-15 | 0.633 | FLT7;GC3;FLT5;GC5;FLT6;GC4;FLT4;FLT6;FLT8;GC8;GC6;FLT3 |
| up_in_flight | Synthesis Of Dna | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 2.22e-06 | 9.51e-12 | 0.614 | FLT7;GC3;FLT1;GC6;GC5;FLT6;GC3;FLT4;GC4;FLT8;GC8;GC6;FLT3;GC2 |
| up_in_flight | P53 Dependent G1 Dna Damage Response | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 3.94e-06 | 4.04e-18 | 0.739 | GC3;FLT7;FLT5;GC5;FLT6;GC4;FLT4;FLT6;FLT8;GC8;GC6;FLT3 |
| up_in_flight | Autodegradation Of The E3 Ubiquitin Ligase Cop1 | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 6.57e-06 | 9.42e-21 | 0.802 | GC3;FLT7;FLT5;FLT10;GC5;FLT6;FLT5;GC4;FLT4;FLT6;FLT8;GC6;GC8;FLT3;GC7 |
| up_in_flight | Activation Of Nf Kappab In B Cells | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 1e-05 | 3.55e-21 | 0.627 | GC3;FLT7;FLT5;FLT10;GC5;FLT6;FLT5;GC6;GC4;FLT4;GC7;FLT6;FLT8;GC6;FLT3;GC8;GC7 |
| up_in_flight | Cdk Mediated Phosphorylation And Removal Of Cdc6 | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 1.03e-05 | 9.42e-21 | 0.76 | GC3;FLT7;FLT5;FLT10;GC5;FLT6;FLT5;GC4;FLT4;FLT6;FLT8;GC6;GC8;FLT3;GC7 |
| up_in_flight | P53 Independent G1 S Dna Damage Checkpoint | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 1.89e-05 | 8.74e-20 | 0.741 | GC3;FLT7;FLT5;GC5;FLT6;FLT5;GC4;FLT4;FLT6;FLT8;GC6;GC8;FLT3 |
| up_in_flight | Cell Cycle Checkpoints | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 1.95e-05 | 2.41e-10 | 0.423 | FLT7;GC3;GC5;FLT6;FLT3;FLT4;GC4;GC8;FLT8;GC6;FLT3;GC2;GC4;FLT4 |
| up_in_flight | Autodegradation Of Cdh1 By Cdh1 Apc C | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 2.51e-05 | 2.1e-17 | 0.744 | GC5;FLT6;GC6;FLT3;FLT5;GC3;FLT7;FLT5;FLT4;GC4;GC7;FLT6;FLT2;FLT8;GC8;GC6;FLT3 |
| up_in_flight | Cdt1 Association With The Cdc6 Orc Origin Complex | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 3.61e-05 | 2.27e-17 | 0.689 | GC3;FLT7;FLT5;GC5;FLT6;GC4;FLT4;FLT6;FLT8;GC8;GC6;FLT3 |
| up_in_flight | Apc C Cdh1 Mediated Degradation Of Cdc20 And Other Apc C Cdh1 Targeted Proteins In Late Mitosis Early G1 | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 4.89e-05 | 5.19e-17 | 0.686 | GC5;FLT6;GC6;FLT7;GC3;FLT5;FLT4;GC4;GC7;FLT2;FLT6;FLT8;GC8;GC6;FLT3 |
| up_in_flight | Regulation Of Ornithine Decarboxylase Odc | 4 | OSD-137,OSD-379,OSD-463,OSD-48 | 8.05e-05 | 8.74e-20 | 0.722 | GC3;FLT7;FLT5;GC5;FLT6;FLT5;GC4;FLT4;FLT6;FLT8;GC6;GC8;FLT3;GC7 |

## Per-Study DGEA/GLARE Pathway Overlaps

| accession | direction | clean_term | fdr_bh | mean_wald_stat_shift | glare_best_fdr_bh | glare_locations | glare_example_clusters |
| --- | --- | --- | --- | --- | --- | --- | --- |
| OSD-137 | down_in_flight | Transmembrane Transport Of Small Molecules | 8.29e-12 | -0.34 | 0.000712 | FLT,GC | FLT3;GC10;GC6;FLT8 |
| OSD-137 | up_in_flight | Olfactory Signaling Pathway | 1.09e-08 | 0.486 | 1.52e-101 | FLT,GC | GC0;FLT1;FLT5;FLT0;FLT10 |
| OSD-137 | down_in_flight | Metabolism Of Lipids And Lipoproteins | 1.1e-07 | -0.217 | 1.26e-22 | FLT,GC | FLT3;GC10;GC6;FLT8;GC9;FLT13;GC8;GC5 |
| OSD-137 | down_in_flight | Post Translational Protein Modification | 8.12e-07 | -0.355 | 0.000174 | FLT,GC | GC6;FLT8;GC1;FLT3;FLT7 |
| OSD-137 | down_in_flight | Integrin Cell Surface Interactions | 8.12e-07 | -0.58 | 0.0034 | FLT,GC | FLT3;GC8;GC3 |
| OSD-137 | down_in_flight | Response To Elevated Platelet Cytosolic Ca2  | 9.06e-07 | -0.535 | 8.36e-15 | FLT,GC | GC10;FLT3;FLT8;GC6;GC8 |
| OSD-137 | down_in_flight | Platelet Activation Signaling And Aggregation | 1.01e-06 | -0.337 | 1.95e-09 | FLT,GC | GC10;FLT3;GC9;FLT8;GC8 |
| OSD-137 | down_in_flight | Slc Mediated Transmembrane Transport | 2.27e-06 | -0.32 | 0.028 | GC | GC9;GC10 |
| OSD-137 | down_in_flight | Hemostasis | 5.11e-06 | -0.212 | 4.22e-08 | FLT,GC | FLT3;GC10;GC9;FLT8;GC7 |
| OSD-137 | down_in_flight | Phospholipid Metabolism | 1.96e-05 | -0.284 | 4.69e-05 | FLT,GC | FLT3;GC5;GC9;GC6;FLT9;GC7;FLT4 |
| OSD-137 | down_in_flight | Formation Of Fibrin Clot Clotting Cascade | 4.85e-05 | -0.698 | 3.8e-12 | FLT,GC | GC10;FLT3 |
| OSD-137 | down_in_flight | Signaling By Pdgf | 0.000114 | -0.336 | 0.000406 | FLT,GC | GC8;FLT4;FLT3;GC9;GC5;GC4 |
| OSD-137 | down_in_flight | Sphingolipid Metabolism | 0.000159 | -0.505 | 0.000142 | FLT,GC | GC5;FLT9;FLT3;GC8 |
| OSD-137 | down_in_flight | Metabolism Of Carbohydrates | 0.000287 | -0.26 | 9.8e-05 | FLT,GC | FLT3;GC10;GC6;GC5;GC9;FLT8;FLT9 |
| OSD-137 | up_in_flight | Synthesis Of Dna | 0.000297 | 0.385 | 2.34e-06 | FLT,GC | FLT8;GC8;GC6;FLT3;GC2 |
| OSD-137 | down_in_flight | Axon Guidance | 0.000351 | -0.255 | 0.00406 | FLT,GC | FLT3;FLT4;GC4;GC6;GC8 |
| OSD-137 | down_in_flight | Glycosphingolipid Metabolism | 0.00038 | -0.753 | 0.00172 | FLT,GC | GC5;FLT9;GC1;FLT2 |
| OSD-137 | down_in_flight | Lipoprotein Metabolism | 0.000767 | -0.658 | 1.78e-08 | FLT,GC | GC10;GC6;FLT3;FLT8 |
| OSD-137 | down_in_flight | Developmental Biology | 0.00102 | -0.181 | 0.000202 | FLT,GC | FLT3;FLT4;GC4;GC6;GC8;GC5 |
| OSD-137 | down_in_flight | Biological Oxidations | 0.00114 | -0.303 | 2.41e-16 | FLT,GC | GC10;FLT8;FLT3;GC6;GC9;GC8 |
| OSD-137 | down_in_flight | Post Translational Modification Synthesis Of Gpi Anchored Proteins | 0.0012 | -0.567 | 0.00103 | FLT,GC | FLT2;GC1;GC4 |
| OSD-137 | down_in_flight | Abca Transporters In Lipid Homeostasis | 0.00123 | -0.701 | 0.00914 | FLT | FLT3 |
| OSD-137 | down_in_flight | Lipid Digestion Mobilization And Transport | 0.00159 | -0.481 | 8.05e-07 | FLT,GC | GC10;GC6;FLT8;FLT3 |
| OSD-137 | down_in_flight | Initial Triggering Of Complement | 0.00159 | -0.762 | 0.00481 | FLT | FLT3 |
| OSD-137 | down_in_flight | Cell Surface Interactions At The Vascular Wall | 0.00159 | -0.384 | 0.0161 | GC | GC7 |
| OSD-137 | down_in_flight | Asparagine N Linked Glycosylation | 0.00216 | -0.336 | 1.66e-06 | FLT,GC | GC6;FLT8;FLT3;FLT7 |
| OSD-137 | down_in_flight | Intrinsic Pathway | 0.00242 | -0.784 | 1.04e-05 | FLT,GC | GC10;FLT3 |
| OSD-137 | down_in_flight | Hyaluronan Uptake And Degradation | 0.003 | -1.29 | 0.00964 | FLT | FLT4 |
| OSD-137 | down_in_flight | Glycerophospholipid Biosynthesis | 0.00323 | -0.31 | 1.56e-05 | FLT,GC | GC6;FLT8;FLT3;GC9 |
| OSD-137 | down_in_flight | Hdl Mediated Lipid Transport | 0.00487 | -0.788 | 0.0147 | FLT | FLT3 |
| OSD-137 | down_in_flight | Abc Family Proteins Mediated Transport | 0.005 | -0.458 | 0.000176 | FLT,GC | GC10;FLT3;FLT8 |
| OSD-137 | down_in_flight | Amyloids | 0.00649 | -0.773 | 8.79e-09 | FLT,GC | FLT1;GC10;GC0 |
| OSD-137 | up_in_flight | Mrna Splicing | 0.00671 | 0.309 | 1.09e-05 | FLT,GC | FLT4;GC5;FLT3;GC1;FLT9;GC4;GC8;GC6 |
| OSD-137 | up_in_flight | Dna Strand Elongation | 0.00679 | 0.58 | 5.19e-06 | FLT,GC | FLT2;GC2;FLT7 |
| OSD-137 | up_in_flight | M G1 Transition | 0.00701 | 0.315 | 3.87e-07 | FLT,GC | FLT8;GC8;GC6;FLT3 |
| OSD-137 | up_in_flight | Assembly Of The Pre Replicative Complex | 0.00715 | 0.346 | 2.52e-08 | FLT,GC | FLT8;GC8;GC6;FLT3 |
| OSD-137 | down_in_flight | Iron Uptake And Transport | 0.00717 | -0.407 | 0.0143 | FLT,GC | GC5;FLT3;FLT4 |
| OSD-137 | up_in_flight | Activation Of Nf Kappab In B Cells | 0.00726 | 0.349 | 2.7e-09 | FLT,GC | FLT8;GC6;FLT3;GC8;GC7 |
| OSD-137 | up_in_flight | Metabolism Of Amino Acids And Derivatives | 0.00914 | 0.204 | 7.48e-23 | FLT,GC | GC10;FLT3;GC6;FLT8;GC8;GC7 |
| OSD-137 | up_in_flight | Scf Beta Trcp Mediated Degradation Of Emi1 | 0.0103 | 0.361 | 8.46e-10 | FLT,GC | FLT8;GC6;GC8;FLT3 |

## Interpretation

DESeq2 and GLARE are being used for different evidence layers.
DESeq2 is the per-study sample-level FLT-vs-GC statistical test.
GLARE is treated as module discovery: it is most useful where its
clusters are enriched for DESeq2 genes or recover Reactome pathways
that are also shifted in per-study DESeq2 rankings.

Small studies can show module structure but have limited DESeq2 power;
interpret `OSD-48`, `OSD-137`, and `OSD-168` as support/sensitivity
unless their signals recur in the larger studies.
