# Moderate Candidate Deep Check

## Bottom Line

- **Soleus circadian** stays moderate, not strong. It has significant combined score support and two of three studies are FLT-lower, but OSD-714 reverses direction and gene-level DGEA is mixed rather than a clean coordinated circadian suppression.
- **Kidney fatty-acid/TAG/ketone metabolism** also stays moderate. It is FLT-higher in five of six studies and has positive random-set support, but only OSD-457/OSD-771 approach per-study score significance and DGEA direction is mixed across lipid genes.

## Summary Metrics

| candidate | module_genes_exported | studies_tested | mean_flight_minus_ground | direction_consistency | combined_welch_fdr_bh | median_empirical_abs_p | mean_random_effect_z | positive_shift_studies | negative_shift_studies | strict_supported_studies | direction_call |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| soleus_circadian | 49 | 3 | -0.0882 | 0.667 | 2.33e-04 | 0.07 | -1.13 | 1 | 2 | 1 | mixed/reorganized DEG |
| kidney_fatty_acid_tag_ketone | 1.65e+02 | 6 | 0.0557 | 0.833 | 0.0338 | 0.065 | 3.42 | 5 | 1 | 0 | mixed/reorganized DEG |

## Per-Study Scores

| candidate | accession | module_genes | n_flight | n_ground | flight_minus_ground | welch_p_value | empirical_abs_p | welch_fdr_bh | direction_label | flt_up_sig_genes | gc_up_sig_genes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| soleus_circadian | OSD-104 | 47 | 6 | 6 | -0.124 | 8.97e-04 | 0 | 0.0254 | mixed/reorganized DEG | 12 | 11 |
| soleus_circadian | OSD-714 | 45 | 12 | 9 | 0.0898 | 0.0193 | 0.08 | 0.121 | mixed/reorganized DEG | 7 | 4 |
| soleus_circadian | OSD-770 | 46 | 9 | 9 | -0.23 | 0.0103 | 0.07 | 0.087 | FLT-up DEG-supported | 15 | 5 |
| kidney_fatty_acid_tag_ketone | OSD-102 | 1.65e+02 | 6 | 6 | -0.0457 | 0.106 | 0.13 | 0.32 | mixed/reorganized DEG | 9 | 6 |
| kidney_fatty_acid_tag_ketone | OSD-163 | 1.65e+02 | 6 | 6 | 0.00309 | 0.961 | 0.74 | 0.982 | no clear DGEA direction | 0 | 0 |
| kidney_fatty_acid_tag_ketone | OSD-253 | 1.65e+02 | 20 | 19 | 0.144 | 0.322 | 0 | 0.611 | no clear DGEA direction | 1 | 0 |
| kidney_fatty_acid_tag_ketone | OSD-457 | 1.65e+02 | 6 | 6 | 0.0655 | 0.00911 | 0 | 0.0802 | mixed/reorganized DEG | 4 | 4 |
| kidney_fatty_acid_tag_ketone | OSD-513 | 1.65e+02 | 9 | 9 | 0.0236 | 0.852 | 0.24 | 0.946 | mixed/reorganized DEG | 21 | 23 |
| kidney_fatty_acid_tag_ketone | OSD-771 | 1.65e+02 | 20 | 20 | 0.144 | 0.00678 | 0 | 0.0689 | FLT-up DEG-supported | 4 | 0 |

## Top Recurrent Gene-Level Signals

### soleus_circadian

| symbol | gene_id | studies_tested | sig_studies | flt_up_sig_studies | gc_up_sig_studies | median_log2fc_tested | min_padj | accessions_flt_up_sig | accessions_gc_up_sig |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Npas2 | ENSMUSG00000026077 | 3 | 2 | 2 | 0 | 1.01 | 1.96e-08 | OSD-104,OSD-770 |  |
| Serpine1 | ENSMUSG00000037411 | 3 | 2 | 1 | 1 | 0.867 | 1.49e-15 | OSD-714 | OSD-104 |
| Dbp | ENSMUSG00000059824 | 3 | 2 | 0 | 2 | -1.47 | 8.28e-09 |  | OSD-104,OSD-770 |
| Bhlhe41 | ENSMUSG00000030256 | 3 | 1 | 1 | 0 | 0.0662 | 9.12e-16 | OSD-770 |  |
| Bmal1 | ENSMUSG00000055116 | 3 | 1 | 1 | 0 | 0.489 | 8.02e-09 | OSD-770 |  |
| Nr1d1 | ENSMUSG00000020889 | 3 | 1 | 1 | 0 | 0.0635 | 4.42e-06 | OSD-770 |  |
| Bhlhe40 | ENSMUSG00000030103 | 3 | 1 | 1 | 0 | -0.257 | 4.19e-05 | OSD-714 |  |
| Ppara | ENSMUSG00000022383 | 3 | 1 | 0 | 1 | -0.757 | 1.05e-11 |  | OSD-770 |
| Per2 | ENSMUSG00000055866 | 3 | 1 | 0 | 1 | -0.682 | 1.88e-06 |  | OSD-104 |
| Ncoa1 | ENSMUSG00000020647 | 3 | 0 | 0 | 0 | 0.0874 | 7.03e-16 |  |  |
| Smarcd3 | ENSMUSG00000028949 | 3 | 0 | 0 | 0 | 0.401 | 1.66e-14 |  |  |
| Cry2 | ENSMUSG00000068742 | 3 | 0 | 0 | 0 | -0.33 | 8.27e-12 |  |  |
| Cry1 | ENSMUSG00000020038 | 3 | 0 | 0 | 0 | 0.134 | 4.72e-09 |  |  |
| Clock | ENSMUSG00000029238 | 3 | 0 | 0 | 0 | 0.304 | 6.36e-09 |  |  |
| Srebf1 | ENSMUSG00000020538 | 3 | 0 | 0 | 0 | -0.251 | 8.73e-08 |  |  |
| Hif1a | ENSMUSG00000021109 | 3 | 0 | 0 | 0 | 0.168 | 8.92e-08 |  |  |
| Csnk1e | ENSMUSG00000022433 | 3 | 0 | 0 | 0 | 0.319 | 7.91e-07 |  |  |
| Crebbp | ENSMUSG00000022521 | 3 | 0 | 0 | 0 | 0.207 | 1.87e-06 |  |  |
| Rora | ENSMUSG00000032238 | 3 | 0 | 0 | 0 | 0.309 | 1.97e-06 |  |  |
| Cul1 | ENSMUSG00000029686 | 3 | 0 | 0 | 0 | 0.0779 | 7.47e-06 |  |  |

### kidney_fatty_acid_tag_ketone

| symbol | gene_id | studies_tested | sig_studies | flt_up_sig_studies | gc_up_sig_studies | median_log2fc_tested | min_padj | accessions_flt_up_sig | accessions_gc_up_sig |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Npas2 | ENSMUSG00000026077 | 6 | 2 | 2 | 0 | 0.59 | 2.71e-09 | OSD-102,OSD-771 |  |
| Bmal1 | ENSMUSG00000055116 | 6 | 2 | 2 | 0 | 0.451 | 1.99e-07 | OSD-102,OSD-771 |  |
| Ugt1a10 | ENSMUSG00000090165 | 6 | 2 | 2 | 0 | 0.00634 | 6.73e-07 | OSD-457,OSD-513 |  |
| Hmgcs2 | ENSMUSG00000027875 | 6 | 2 | 1 | 1 | 0.161 | 3.10e-06 | OSD-513 | OSD-102 |
| Fabp1 | ENSMUSG00000054422 | 6 | 1 | 1 | 0 | -0.439 | 1.85e-07 | OSD-513 |  |
| Apoa1 | ENSMUSG00000032083 | 6 | 1 | 1 | 0 | 0.339 | 6.00e-05 | OSD-513 |  |
| Apoa5 | ENSMUSG00000032079 | 2 | 1 | 1 | 0 | 1.12 | 0.0383 | OSD-771 |  |
| Ep300 | ENSMUSG00000055024 | 6 | 0 | 0 | 0 | -0.0473 | 1.18e-10 |  |  |
| Slc25a20 | ENSMUSG00000032602 | 6 | 0 | 0 | 0 | -0.123 | 1.77e-08 |  |  |
| Decr1 | ENSMUSG00000028223 | 6 | 0 | 0 | 0 | -0.0785 | 1.02e-07 |  |  |
| Slc27a1 | ENSMUSG00000031808 | 6 | 0 | 0 | 0 | 0.0568 | 1.77e-06 |  |  |
| Nr1d1 | ENSMUSG00000020889 | 6 | 0 | 0 | 0 | -0.0337 | 2.32e-06 |  |  |
| Gpd1l | ENSMUSG00000050627 | 6 | 0 | 0 | 0 | 7.19e-04 | 2.57e-06 |  |  |
| Oxct1 | ENSMUSG00000022186 | 6 | 0 | 0 | 0 | -0.14 | 2.49e-05 |  |  |
| Fhl2 | ENSMUSG00000008136 | 6 | 0 | 0 | 0 | 0.0759 | 2.61e-05 |  |  |
| Tead1 | ENSMUSG00000055320 | 6 | 0 | 0 | 0 | 0.0284 | 3.61e-05 |  |  |
| Ncor1 | ENSMUSG00000018501 | 6 | 0 | 0 | 0 | 6.90e-04 | 7.60e-05 |  |  |
| Bdh1 | ENSMUSG00000046598 | 6 | 0 | 0 | 0 | 0.0238 | 8.43e-05 |  |  |
| Gpd2 | ENSMUSG00000026827 | 6 | 0 | 0 | 0 | 0.237 | 8.76e-05 |  |  |
| Elovl2 | ENSMUSG00000021364 | 6 | 0 | 0 | 0 | -0.0196 | 1.07e-04 |  |  |

## Interpretation

### Soleus circadian

This is plausible but not clean enough to claim as a strong hidden module. The module-level score is FLT-lower in OSD-104 and OSD-770, but FLT-higher in OSD-714. Gene-level DGEA is mixed: across studies there are more FLT-up significant genes than GC-up genes, even though the module score is mostly FLT-lower. That pattern suggests circadian reorganization or phase/design effects rather than uniform suppression.

### Kidney fatty-acid/TAG/ketone

This is more coherent at the module-score level than at the gene-level DGEA level. Five of six studies have positive FLT-GC module shifts, but only the random-set comparison is consistently favorable; per-study Welch FDR is borderline. Gene-level DGEA includes both FLT-up and GC-up lipid genes, so the safest claim is renal lipid-handling reorganization, not simple activation of fatty-acid/ketone metabolism.

## Retained Files

- `moderate_candidate_deep_check_summary.tsv`
- `moderate_candidate_per_study_scores.tsv`
