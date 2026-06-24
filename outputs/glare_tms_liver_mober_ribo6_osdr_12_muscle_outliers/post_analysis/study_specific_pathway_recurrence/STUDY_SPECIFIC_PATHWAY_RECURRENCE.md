# Study-Specific FLT-vs-GC Pathway Recurrence

This analysis treats OSD-379, OSD-245, and OSD-463 as separate primary studies.
Aggregate GLARE is used only as exploratory support after study-specific raw-count DESeq2 and Reactome pathway testing.

## Inputs

- Raw expression HDF5: `assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5`
- Liver metadata: `data/processed/osdr_mouse_bulk_liver.profile_metadata.tsv`
- Excluded muscle-outlier profiles: `data/filters/aggregate_liver_12_muscle_candidate_profiles.txt`
- Reactome GMT: `src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt`
- DESeq2/pathway alpha: 0.05
- Rank-sum pathway size range: 10-500 tested genes

## Sample Counts After 12-Outlier Filter

| accession | flight | ground |
| --- | --- | --- |
| OSD-245 | 18 | 19 |
| OSD-379 | 31 | 31 |
| OSD-463 | 12 | 10 |

## Study-Internal Strata

| accession | stratum | flight | ground |
| --- | --- | --- | --- |
| OSD-245 | ISS-T | 8 | 10 |
| OSD-245 | LAR | 10 | 9 |
| OSD-379 | ISS-T_OLD | 8 | 10 |
| OSD-379 | ISS-T_YNG | 8 | 9 |
| OSD-379 | LAR_OLD | 7 | 6 |
| OSD-379 | LAR_YNG | 8 | 6 |
| OSD-463 | all | 12 | 10 |

## DESeq2 Summary

| accession | n_flight | n_ground | genes_tested | significant_padj05 | significant_up | significant_down | design | dispersion_fit |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| OSD-245 | 18 | 19 | 24291 | 735 | 288 | 447 | ~ stratum + condition | default |
| OSD-379 | 31 | 31 | 22370 | 1632 | 1087 | 545 | ~ stratum + condition | default |
| OSD-463 | 12 | 10 | 22420 | 2908 | 1603 | 1305 | ~ condition | default |

## Recurring Reactome Pathways

| direction | clean_term | study_count | accessions | best_fdr_bh | mean_wald_stat_shift_mean | aggregate_glare_support_categories |
| --- | --- | --- | --- | --- | --- | --- |
| down_in_flight | Immunoregulatory Interactions Between A Lymphoid And A Non Lymphoid Cell | 3 | OSD-245,OSD-379,OSD-463 | 4.39e-11 | -1.21 | immune_complement_coagulation |
| down_in_flight | Inflammasomes | 3 | OSD-245,OSD-379,OSD-463 | 0.00104 | -1.4 | immune_complement_coagulation |
| down_in_flight | Generation Of Second Messenger Molecules | 3 | OSD-245,OSD-379,OSD-463 | 0.00301 | -0.979 |  |
| up_in_flight | Mrna Processing | 2 | OSD-245,OSD-379 | 4.81e-15 | 0.714 |  |
| up_in_flight | Processing Of Capped Intron Containing Pre Mrna | 2 | OSD-245,OSD-379 | 8.06e-14 | 0.738 |  |
| up_in_flight | Metabolism Of Amino Acids And Derivatives | 2 | OSD-379,OSD-463 | 6.93e-12 | 0.733 |  |
| up_in_flight | Cell Cycle | 2 | OSD-379,OSD-463 | 1.97e-08 | 0.385 |  |
| up_in_flight | Cell Cycle Mitotic | 2 | OSD-379,OSD-463 | 7.24e-08 | 0.383 |  |
| down_in_flight | Formation Of Fibrin Clot Clotting Cascade | 2 | OSD-245,OSD-463 | 2.9e-07 | -1.23 | immune_complement_coagulation |
| up_in_flight | Transport Of Mature Transcript To Cytoplasm | 2 | OSD-245,OSD-379 | 7.59e-07 | 0.906 |  |
| up_in_flight | Synthesis Of Dna | 2 | OSD-379,OSD-463 | 2.3e-06 | 0.676 |  |
| up_in_flight | Transport Of Mature Mrna Derived From An Intronless Transcript | 2 | OSD-245,OSD-379 | 6.48e-06 | 1.06 |  |
| down_in_flight | Interferon Gamma Signaling | 2 | OSD-245,OSD-463 | 1.05e-05 | -1.09 | immune_complement_coagulation |
| down_in_flight | Innate Immune System | 2 | OSD-245,OSD-463 | 2.06e-05 | -0.463 | immune_complement_coagulation |
| down_in_flight | Cholesterol Biosynthesis | 2 | OSD-245,OSD-463 | 2.28e-05 | -1.64 | lipid_steroid_xenobiotic_metabolism |
| up_in_flight | Apc C Cdc20 Mediated Degradation Of Mitotic Proteins | 2 | OSD-379,OSD-463 | 3.38e-05 | 0.827 |  |
| up_in_flight | Regulation Of Mitotic Cell Cycle | 2 | OSD-379,OSD-463 | 3.38e-05 | 0.755 |  |
| down_in_flight | Complement Cascade | 2 | OSD-245,OSD-463 | 4.15e-05 | -1.3 | immune_complement_coagulation |
| down_in_flight | Chemokine Receptors Bind Chemokines | 2 | OSD-245,OSD-463 | 6.32e-05 | -1.2 |  |
| up_in_flight | Signaling By Wnt | 2 | OSD-379,OSD-463 | 6.81e-05 | 0.799 |  |

## Top Per-Study Rank-Sum Pathways

Reactome rank-sum tests ask whether all genes in a pathway are shifted up or down in the DESeq2 Wald-statistic ranking.

| accession | direction | clean_term | pathway_genes_tested | mean_wald_stat_shift | p_value | fdr_bh |
| --- | --- | --- | --- | --- | --- | --- |
| OSD-245 | down_in_flight | Peptide Chain Elongation | 56 | -1.95 | 3.73e-25 | 2.26e-22 |
| OSD-245 | down_in_flight | Influenza Viral Rna Transcription And Replication | 70 | -1.7 | 2.42e-24 | 6.93e-22 |
| OSD-245 | down_in_flight | Srp Dependent Cotranslational Protein Targeting To Membrane | 79 | -1.64 | 3.43e-24 | 6.93e-22 |
| OSD-245 | down_in_flight | Metabolism Of Proteins | 361 | -0.716 | 3.81e-23 | 5.77e-21 |
| OSD-245 | down_in_flight | Translation | 117 | -1.22 | 1.17e-20 | 1.41e-18 |
| OSD-245 | up_in_flight | Olfactory Signaling Pathway | 195 | 0.676 | 2.97e-19 | 1.8e-16 |
| OSD-245 | up_in_flight | Circadian Repression Of Expression By Rev Erba | 22 | 2.01 | 2.98e-06 | 0.0009 |
| OSD-245 | up_in_flight | Transport Of Mature Transcript To Cytoplasm | 51 | 0.694 | 8.17e-06 | 0.00165 |
| OSD-245 | up_in_flight | Processing Of Capped Intron Containing Pre Mrna | 127 | 0.377 | 4.57e-05 | 0.00673 |
| OSD-245 | up_in_flight | Mrna 3 End Processing | 34 | 0.785 | 7.46e-05 | 0.00673 |
| OSD-379 | down_in_flight | Gpcr Downstream Signaling | 452 | -0.416 | 4.57e-10 | 2.73e-07 |
| OSD-379 | down_in_flight | Immunoregulatory Interactions Between A Lymphoid And A Non Lymphoid Cell | 54 | -1.17 | 7.71e-09 | 2.3e-06 |
| OSD-379 | down_in_flight | Olfactory Signaling Pathway | 158 | -0.46 | 2.61e-05 | 0.00518 |
| OSD-379 | down_in_flight | G Alpha Q Signalling Events | 108 | -0.489 | 0.000221 | 0.033 |
| OSD-379 | down_in_flight | Class A1 Rhodopsin Like Receptors | 141 | -0.434 | 0.000289 | 0.0344 |
| OSD-379 | up_in_flight | Mrna Processing | 147 | 1.07 | 8.07e-18 | 4.81e-15 |
| OSD-379 | up_in_flight | Processing Of Capped Intron Containing Pre Mrna | 128 | 1.1 | 2.7e-16 | 8.06e-14 |
| OSD-379 | up_in_flight | Muscle Contraction | 44 | 2.16 | 1.68e-14 | 3.35e-12 |
| OSD-379 | up_in_flight | Mrna Splicing | 101 | 1.06 | 5.54e-12 | 8.26e-10 |
| OSD-379 | up_in_flight | Antigen Processing Ubiquitination Proteasome Degradation | 182 | 0.747 | 8.39e-11 | 1e-08 |
| OSD-463 | down_in_flight | Cytokine Signaling In Immune System | 234 | -0.79 | 3.98e-10 | 2.37e-07 |
| OSD-463 | down_in_flight | Interferon Gamma Signaling | 51 | -1.56 | 3.52e-08 | 1.05e-05 |
| OSD-463 | down_in_flight | Innate Immune System | 196 | -0.753 | 1.04e-07 | 2.06e-05 |
| OSD-463 | down_in_flight | Cholesterol Biosynthesis | 21 | -2.16 | 1.53e-07 | 2.28e-05 |
| OSD-463 | down_in_flight | Signaling By Ils | 102 | -0.936 | 2.15e-07 | 2.56e-05 |
| OSD-463 | up_in_flight | Respiratory Electron Transport | 65 | 2.37 | 2.05e-26 | 1.09e-23 |
| OSD-463 | up_in_flight | Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins  | 66 | 2.33 | 3.64e-26 | 1.09e-23 |
| OSD-463 | up_in_flight | Tca Cycle And Respiratory Electron Transport | 102 | 1.84 | 3.08e-24 | 6.11e-22 |
| OSD-463 | up_in_flight | Metabolism Of Amino Acids And Derivatives | 174 | 1.13 | 4.65e-14 | 6.93e-12 |
| OSD-463 | up_in_flight | Mitochondrial Protein Import | 44 | 1.92 | 1.4e-10 | 1.66e-08 |

## ORA On Significant DESeq2 Genes

| accession | direction | clean_term | query_genes | overlap | p_value | fdr_bh |
| --- | --- | --- | --- | --- | --- | --- |
| OSD-245 | down_in_flight | Respiratory Electron Transport | 447 | 10 | 3.39e-06 | 0.000137 |
| OSD-245 | down_in_flight | Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins  | 447 | 10 | 3.92e-06 | 0.000137 |
| OSD-245 | down_in_flight | Tca Cycle And Respiratory Electron Transport | 447 | 11 | 3.84e-05 | 0.000896 |
| OSD-245 | down_in_flight | Mitochondrial Protein Import | 447 | 7 | 9.08e-05 | 0.00159 |
| OSD-245 | down_in_flight | Srp Dependent Cotranslational Protein Targeting To Membrane | 447 | 9 | 0.000124 | 0.00173 |
| OSD-245 | up_in_flight | Circadian Repression Of Expression By Rev Erba | 288 | 6 | 9.08e-07 | 6.45e-05 |
| OSD-245 | up_in_flight | Circadian Clock | 288 | 7 | 1.14e-05 | 0.000403 |
| OSD-245 | up_in_flight | Rora Activates Circadian Expression | 288 | 5 | 2.57e-05 | 0.000609 |
| OSD-245 | up_in_flight | Bmal1 Clock Npas2 Activates Circadian Expression | 288 | 5 | 0.000159 | 0.00226 |
| OSD-245 | up_in_flight | Nuclear Signaling By Erbb4 | 288 | 5 | 0.000159 | 0.00226 |
| OSD-379 | down_in_flight | Metabolism Of Lipids And Lipoproteins | 545 | 31 | 1.34e-06 | 0.000107 |
| OSD-379 | down_in_flight | Fatty Acid Triacylglycerol And Ketone Body Metabolism | 545 | 16 | 3.09e-05 | 0.00124 |
| OSD-379 | down_in_flight | Glycerophospholipid Biosynthesis | 545 | 8 | 0.000522 | 0.0124 |
| OSD-379 | down_in_flight | Xenobiotics | 545 | 5 | 0.000709 | 0.0124 |
| OSD-379 | down_in_flight | Cytochrome P450 Arranged By Substrate Type | 545 | 7 | 0.000772 | 0.0124 |
| OSD-379 | up_in_flight | Striated Muscle Contraction | 1087 | 14 | 1.51e-12 | 3.5e-10 |
| OSD-379 | up_in_flight | Muscle Contraction | 1087 | 15 | 5.63e-09 | 6.53e-07 |
| OSD-379 | up_in_flight | Post Chaperonin Tubulin Folding Pathway | 1087 | 7 | 4.26e-06 | 0.00033 |
| OSD-379 | up_in_flight | Cholesterol Biosynthesis | 1087 | 8 | 6.76e-06 | 0.000392 |
| OSD-379 | up_in_flight | Loss Of Nlp From Mitotic Centrosomes | 1087 | 12 | 2.98e-05 | 0.00138 |
| OSD-463 | down_in_flight | Immune System | 1305 | 100 | 3.05e-11 | 6.29e-09 |
| OSD-463 | down_in_flight | Interferon Gamma Signaling | 1305 | 17 | 1.29e-08 | 1.32e-06 |
| OSD-463 | down_in_flight | Cytokine Signaling In Immune System | 1305 | 40 | 1.93e-08 | 1.32e-06 |
| OSD-463 | down_in_flight | Cholesterol Biosynthesis | 1305 | 9 | 1.96e-06 | 9.92e-05 |
| OSD-463 | down_in_flight | Chemokine Receptors Bind Chemokines | 1305 | 11 | 2.41e-06 | 9.92e-05 |
| OSD-463 | up_in_flight | Respiratory Electron Transport | 1603 | 29 | 2.4e-15 | 4.95e-13 |
| OSD-463 | up_in_flight | Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins  | 1603 | 29 | 3.95e-15 | 4.95e-13 |
| OSD-463 | up_in_flight | Tca Cycle And Respiratory Electron Transport | 1603 | 35 | 5.07e-14 | 4.25e-12 |
| OSD-463 | up_in_flight | Metabolism Of Lipids And Lipoproteins | 1603 | 66 | 1.33e-07 | 8.35e-06 |
| OSD-463 | up_in_flight | Metabolism Of Amino Acids And Derivatives | 1603 | 34 | 7.87e-07 | 3.95e-05 |

## Aggregate GLARE Support

Aggregate GLARE is not used as the main evidence here because prior audits showed strong study/mission/batch structure.
It is only used to flag whether study-specific Reactome terms fall into broad categories also seen in aggregate GLARE/Metascape clusters.
Aggregate GLARE support categories detected: immune_complement_coagulation, lipid_steroid_xenobiotic_metabolism, muscle_cytoskeleton, protein_processing_localization, translation_ribosome

## Interpretation Rule

Use a pathway as stronger evidence only when it appears in the same direction in multiple individual studies.
Use one-study pathways as study-specific findings, and use aggregate GLARE clusters only as exploratory context for gene modules or follow-up hypotheses.
