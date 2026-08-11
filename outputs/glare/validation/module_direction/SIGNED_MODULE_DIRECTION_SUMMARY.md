# Signed Direction Summary For GLARE Candidate Modules

This report adds within-module DGEA direction checks to the validated GLARE candidate modules.
Positive log2FC means higher in spaceflight; negative log2FC means higher in ground control.

The GLARE paper used DEG proportion, enrichment, heatmaps, and follow-up feature analyses to validate modules.
This signed table is a stricter extension for the multi-study mouse analysis, because enrichment alone does not say whether a pathway is activated, suppressed, or mixed.

## Direction Call Counts

| module_class | direction_call | modules |
| --- | --- | --- |
| glare_only | mixed/reorganized DEG | 101 |
| glare_only | ambiguous/no validated direction | 35 |
| glare_only | FLT-up DEG-supported | 23 |
| glare_only | GC-up DEG-supported | 21 |
| intersection | GC-up DEG-supported | 84 |
| intersection | FLT-up DEG-supported | 53 |
| intersection | mixed/reorganized DEG | 38 |
| intersection | ambiguous/no validated direction | 5 |

## Strict GLARE-Only Hidden Direction Candidates

No rows.

## Mixed/Reorganized DEG Modules

| tissue | module_class | clean_term | flt_up_sig_genes_sum | gc_up_sig_genes_sum | combined_welch_fdr_bh | mean_flight_minus_ground |
| --- | --- | --- | --- | --- | --- | --- |
| thymus | glare_only | Cell Death Signalling Via Nrage Nrif And Nade | 32 | 23 | 7.47e-10 | 0.04345 |
| thymus | glare_only | Downstream Signaling Of Activated Fgfr | 43 | 52 | 7.47e-10 | 0.03871 |
| thymus | glare_only | Circadian Clock | 24 | 24 | 1.293e-09 | 0.06845 |
| thymus | glare_only | Golgi Associated Vesicle Biogenesis | 26 | 19 | 2.16e-08 | 0.05332 |
| thymus | glare_only | Biosynthesis Of The N Glycan Precursor Dolichol Lipid Linked Oligosaccharide Llo And Transfer To A Nascent Protein | 9 | 5 | 2.289e-07 | -0.01246 |
| thymus | glare_only | Synthesis Of Pips At The Plasma Membrane | 10 | 14 | 1.849e-06 | -0.01827 |
| skeletal_muscle | glare_only | Bmal1 Clock Npas2 Activates Circadian Expression | 58 | 60 | 2.873e-06 | -0.07905 |
| thymus | glare_only | Iron Uptake And Transport | 24 | 15 | 4.282e-06 | 0.1242 |
| skeletal_muscle | glare_only | Cyclin E Associated Events During G1 S Transition | 131 | 57 | 2.767e-05 | 0.03854 |
| thymus | glare_only | Synthesis Of Substrates In N Glycan Biosythesis | 6 | 3 | 0.0001489 | 0.02134 |
| skeletal_muscle_soleus | glare_only | Circadian Clock | 34 | 20 | 0.0002328 | -0.0882 |
| skeletal_muscle_soleus | glare_only | Cyclin E Associated Events During G1 S Transition | 47 | 31 | 0.0007156 | -0.003754 |
| skeletal_muscle | glare_only | L1Cam Interactions | 94 | 77 | 0.001035 | 0.008231 |
| skeletal_muscle | glare_only | Formation Of Rna Pol Ii Elongation Complex | 45 | 40 | 0.001079 | 0.001188 |
| skeletal_muscle | glare_only | Creb Phosphorylation Through The Activation Of Ras | 25 | 27 | 0.001362 | -0.01133 |
| skeletal_muscle | glare_only | Response To Elevated Platelet Cytosolic Ca2 | 87 | 83 | 0.001382 | 0.00276 |
| skeletal_muscle | glare_only | Trif Mediated Tlr3 Signaling | 85 | 64 | 0.008553 | 0.01612 |
| skeletal_muscle | glare_only | Mitochondrial Trna Aminoacylation | 12 | 18 | 0.01295 | -0.02508 |
| kidney | glare_only | Membrane Trafficking | 13 | 15 | 0.01798 | 0.03382 |
| kidney | glare_only | Signaling By Insulin Receptor | 8 | 18 | 0.02038 | 0.02823 |
| skeletal_muscle_soleus | glare_only | Signalling By Ngf | 111 | 100 | 0.02038 | -0.07781 |
| skeletal_muscle | glare_only | Prolactin Receptor Signaling | 12 | 11 | 0.02038 | 0.003814 |
| liver | glare_only | Cleavage Of Growing Transcript In The Termination Region | 7 | 4 | 0.02462 | -0.01242 |
| kidney | glare_only | Fatty Acid Triacylglycerol And Ketone Body Metabolism | 39 | 33 | 0.03383 | 0.05573 |
| skeletal_muscle_quadriceps | glare_only | Cyclin E Associated Events During G1 S Transition | 20 | 14 | 0.03585 | 0.01906 |
| skeletal_muscle_soleus | glare_only | Response To Elevated Platelet Cytosolic Ca2 | 55 | 25 | 0.05594 | 0.01534 |
| skeletal_muscle_soleus | glare_only | Ngf Signalling Via Trka From The Plasma Membrane | 84 | 63 | 0.05721 | -0.071 |
| kidney | glare_only | Response To Elevated Platelet Cytosolic Ca2 | 22 | 11 | 0.05764 | 0.08573 |
| kidney | glare_only | Rna Pol Iii Transcription | 2 | 3 | 0.05764 | 0.02895 |
| skeletal_muscle_soleus | glare_only | Prolactin Receptor Signaling | 5 | 4 | 0.06582 | -0.02285 |

## Output Files

- `candidate_module_signed_dgea_by_study.tsv`: per-study signed counts for each module.
- `candidate_module_signed_dgea_meta.tsv`: cross-study signed summary joined to module-score validation.
