# Multi-Tissue GLARE Validation Stack

This report validates existing API-derived multi-tissue GLARE outputs using a paper-style stack.

## What Was Run

- XGBoost melted FLT-vs-GC verification with both GLARE-like random folds and gene-grouped audit folds.
- Representation QC for raw PCA, FT-SAE latent, and FT-SAE PCA representations.
- Consensus/base-clustering QC, including sampled average-linkage EAC agreement.
- DEG-enrichment proportion by GMM, HDBSCAN, Spectral, and consensus labels.
- Candidate module-score validation for DGEA-intersection and GLARE-only recurring modules.
- Random-gene-set controls and Metascape-ready gene-list export.
- PanglaoDB marker enrichment as a cell-type-marker proxy for the paper's cell-type follow-up.
- TF/stress-network validation is documented as unavailable here because the repo does not contain a curated mouse spaceflight stress/TF network.

## Key Counts

- Verification summaries: 178 scheme rows.
- Representation QC rows: 534.
- Clustering QC rows: 178.
- Candidate modules tested: 360.
- Module-score meta rows: 360.

## Verification Snapshot

| tissue | scope | scheme | n_folds | accuracy_mean | accuracy_std | f1_mean | f1_std | roc_auc_mean | roc_auc_std | n_genes | balanced_features | shap_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kidney | aggregate | gene_grouped_kfold_audit | 5 | 0.9487 | 0.0007197 | 0.9498 | 0.0006905 | 0.9906 | 0.0004806 | 21024 | 66 | ok |
| kidney | aggregate | random_kfold_glare | 5 | 0.9475 | 0.004022 | 0.9486 | 0.003982 | 0.9905 | 0.001179 | 21024 | 66 | ok |
| liver | aggregate | gene_grouped_kfold_audit | 5 | 0.9424 | 0.003274 | 0.9397 | 0.003506 | 0.9907 | 0.001008 | 20984 | 112 | ok |
| liver | aggregate | random_kfold_glare | 5 | 0.9431 | 0.001983 | 0.9405 | 0.001986 | 0.9909 | 0.0007524 | 20984 | 112 | ok |
| lung | aggregate | gene_grouped_kfold_audit | 5 | 0.9166 | 0.004012 | 0.9135 | 0.004492 | 0.9802 | 0.001584 | 20984 | 37 | ok |
| lung | aggregate | random_kfold_glare | 5 | 0.9177 | 0.001516 | 0.9147 | 0.001952 | 0.9802 | 0.0003812 | 20984 | 37 | ok |
| skeletal_muscle | aggregate | gene_grouped_kfold_audit | 5 | 0.9136 | 0.001735 | 0.9084 | 0.002109 | 0.9793 | 0.001463 | 21024 | 89 | ok |
| skeletal_muscle | aggregate | random_kfold_glare | 5 | 0.9129 | 0.005102 | 0.9077 | 0.005889 | 0.9789 | 0.002154 | 21024 | 89 | ok |
| skeletal_muscle_edl | aggregate | gene_grouped_kfold_audit | 5 | 0.7908 | 0.003832 | 0.7714 | 0.005012 | 0.8805 | 0.004192 | 21024 | 15 | ok |
| skeletal_muscle_edl | aggregate | random_kfold_glare | 5 | 0.7902 | 0.003092 | 0.7709 | 0.00462 | 0.8806 | 0.002952 | 21024 | 15 | ok |
| skeletal_muscle_gastrocnemius | aggregate | gene_grouped_kfold_audit | 5 | 0.8088 | 0.004694 | 0.7881 | 0.006136 | 0.9113 | 0.003287 | 21024 | 13 | ok |
| skeletal_muscle_gastrocnemius | aggregate | random_kfold_glare | 5 | 0.8082 | 0.004939 | 0.788 | 0.007027 | 0.9095 | 0.002483 | 21024 | 13 | ok |
| skeletal_muscle_quadriceps | aggregate | gene_grouped_kfold_audit | 5 | 0.8319 | 0.007219 | 0.8133 | 0.009237 | 0.9283 | 0.004363 | 21024 | 22 | ok |
| skeletal_muscle_quadriceps | aggregate | random_kfold_glare | 5 | 0.8315 | 0.003752 | 0.8131 | 0.003404 | 0.9271 | 0.001707 | 21024 | 22 | ok |
| skeletal_muscle_soleus | aggregate | gene_grouped_kfold_audit | 5 | 0.8901 | 0.00449 | 0.8814 | 0.005294 | 0.9656 | 0.001998 | 21024 | 24 | ok |
| skeletal_muscle_soleus | aggregate | random_kfold_glare | 5 | 0.8886 | 0.004282 | 0.8797 | 0.004657 | 0.9641 | 0.00324 | 21024 | 24 | ok |
| skeletal_muscle_tibialis_anterior | aggregate | gene_grouped_kfold_audit | 5 | 0.7974 | 0.00419 | 0.7831 | 0.005406 | 0.8997 | 0.002679 | 21024 | 15 | ok |
| skeletal_muscle_tibialis_anterior | aggregate | random_kfold_glare | 5 | 0.7935 | 0.003603 | 0.7792 | 0.004353 | 0.8962 | 0.004259 | 21024 | 15 | ok |
| skin | aggregate | gene_grouped_kfold_audit | 5 | 0.9213 | 0.004215 | 0.9248 | 0.003958 | 0.9828 | 0.001457 | 20984 | 69 | ok |
| skin | aggregate | random_kfold_glare | 5 | 0.921 | 0.002978 | 0.9245 | 0.003036 | 0.9824 | 0.0006923 | 20984 | 69 | ok |
| spleen | aggregate | gene_grouped_kfold_audit | 5 | 0.9551 | 0.001539 | 0.954 | 0.001553 | 0.9944 | 0.0002725 | 21024 | 52 | ok |
| spleen | aggregate | random_kfold_glare | 5 | 0.9543 | 0.001853 | 0.9532 | 0.002181 | 0.9942 | 0.0006457 | 21024 | 52 | ok |
| thymus | aggregate | gene_grouped_kfold_audit | 5 | 0.9701 | 0.001941 | 0.9694 | 0.002025 | 0.9976 | 0.0003623 | 21024 | 53 | ok |
| thymus | aggregate | random_kfold_glare | 5 | 0.9701 | 0.001562 | 0.9693 | 0.001633 | 0.9975 | 0.0003611 | 21024 | 53 | ok |

## Strongest Module-Score Meta Results

| tissue | module_class | clean_term | studies_tested | combined_welch_fdr_bh | mean_flight_minus_ground | direction_consistency | median_empirical_abs_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| thymus | glare_only | Cell Death Signalling Via Nrage Nrif And Nade | 5 | 7.47e-10 | 0.04345 | 0.6 | 0.05 |
| thymus | glare_only | Downstream Signaling Of Activated Fgfr | 5 | 7.47e-10 | 0.03871 | 0.6 | 0.07 |
| thymus | glare_only | Circadian Clock | 5 | 1.293e-09 | 0.06845 | 0.6 | 0.08 |
| thymus | glare_only | Golgi Associated Vesicle Biogenesis | 5 | 2.16e-08 | 0.05332 | 0.6 | 0.08 |
| thymus | glare_only | Biosynthesis Of The N Glycan Precursor Dolichol Lipid Linked Oligosaccharide Llo And Transfer To A Nascent Protein | 5 | 2.289e-07 | -0.01246 | 0.4 | 0.13 |
| thymus | glare_only | Response To Elevated Platelet Cytosolic Ca2 | 5 | 8.398e-07 | 0.2642 | 0.8 | 0 |
| thymus | glare_only | Synthesis Of Pips At The Plasma Membrane | 5 | 1.849e-06 | -0.01827 | 0.4 | 0.39 |
| skeletal_muscle | glare_only | Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins | 13 | 2.125e-06 | -0.1009 | 0.6154 | 0.01 |
| skeletal_muscle | glare_only | Bmal1 Clock Npas2 Activates Circadian Expression | 13 | 2.873e-06 | -0.07905 | 0.8462 | 0.01 |
| thymus | glare_only | Iron Uptake And Transport | 5 | 4.282e-06 | 0.1242 | 0.8 | 0.17 |
| thymus | glare_only | Pkb Mediated Events | 5 | 6.928e-06 | -0.01812 | 0.4 | 0.19 |
| skeletal_muscle_soleus | glare_only | Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins | 3 | 1.057e-05 | -0.3304 | 1 | 0 |
| skeletal_muscle_soleus | intersection | Mitochondrial Fatty Acid Beta Oxidation | 3 | 4.944e-08 | -0.5805 | 1 | 0 |
| skeletal_muscle_soleus | intersection | N Glycan Trimming In The Er And Calnexin Calreticulin Cycle | 3 | 4.944e-08 | -0.1561 | 0.6667 | 0 |
| thymus | intersection | Mitotic Prometaphase | 5 | 6.554e-08 | -0.5876 | 1 | 0 |
| skeletal_muscle | intersection | Orc1 Removal From Chromatin | 13 | 6.849e-07 | 0.05531 | 0.7692 | 0 |
| skeletal_muscle | intersection | Respiratory Electron Transport | 13 | 6.849e-07 | -0.1017 | 0.6154 | 0 |
| skeletal_muscle | intersection | Metabolism Of Proteins | 13 | 1.609e-06 | -0.01179 | 0.5385 | 0 |
| skeletal_muscle | intersection | G1 S Transition | 13 | 1.849e-06 | 0.03675 | 0.6923 | 0.01 |
| skeletal_muscle | intersection | Post Translational Protein Modification | 13 | 1.849e-06 | -0.02016 | 0.5385 | 0.03 |
| skeletal_muscle | intersection | Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins  | 13 | 2.125e-06 | -0.1009 | 0.6154 | 0.01 |
| skeletal_muscle | intersection | Downstream Signaling Events Of B Cell Receptor Bcr | 13 | 5.118e-06 | 0.0322 | 0.7692 | 0.01 |
| skeletal_muscle_soleus | intersection | Asparagine N Linked Glycosylation | 3 | 5.622e-06 | -0.1253 | 0.6667 | 0 |
| skeletal_muscle | intersection | Mitotic G1 G1 S Phases | 13 | 6.861e-06 | 0.02945 | 0.6923 | 0.02 |

## Interpretation Rules

- Intersection modules have the strongest support because they recur in both per-study DGEA and GLARE cluster enrichment.
- GLARE-only modules are candidate hidden modules only when module-score tests are consistent across studies and stronger than random gene sets.
- Generic GPCR, defensin, viral, and similarly broad labels are excluded from candidate hidden-module selection.
- Liver olfactory/chemosensory labels are retained as high-caution candidates because liver expression can be biologically relevant, but large receptor gene families can also dominate enrichment.
- Sampled EAC is a scalability audit of the GLARE paper's average-linkage EAC idea, not a full dense 20k-gene co-association matrix.
