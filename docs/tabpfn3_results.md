# TabPFN3 Results

## Status

The TabPFN3 OSDR-only production run completed on the local A100 using the
official `tabpfn==8.0.8` package with model version `v3`.

Output root:

```text
outputs/tabpfn3_osdr
```

Key summary files:

- `outputs/tabpfn3_osdr/summary/tabpfn3_aggregate_metrics.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_metrics.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_predictions.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_feature_importance.tsv`
- `outputs/tabpfn3_osdr/summary/tabpfn3_backend_status.json`
- `outputs/tabpfn3_osdr/summary/osdr_sample_inventory.tsv`

The manifest also contains older `blocked_backend_unavailable` rows from
pre-token attempts. The completed production pass is represented by 26
`completed` dataset/mode rows: 13 datasets times 2 feature modes.

## Production Settings

The completed run used:

```bash
PYTHONPATH=src python -m nasa_mouse_tabpfn3.run_osdr_classification \
  --output-root outputs/tabpfn3_osdr \
  --backend tabpfn \
  --device cuda \
  --feature-modes all_expressed hvg \
  --cv-schemes random grouped loo_accession \
  --hvg-top-n 500 \
  --max-features 500 \
  --importance-candidates 5 \
  --permutation-repeats 1 \
  --n-estimators 3
```

Because `--max-features 500` and `--hvg-top-n 500` were both used, the
`all_expressed` and `hvg` tracks selected the same fold-local top-variance
500-gene feature sets. The current production result is therefore a practical
500-feature TabPFN3 run, not an uncapped all-gene run.

## Main Tissue Results

Main interpretation should use accession-aware grouped CV and
leave-one-accession-out CV. Random CV is useful as an optimism/leakage check
only because samples from the same accession can appear in both train and test.

The table below reports the `hvg` track; `all_expressed` matched under the
500-feature cap.

| dataset | cv | n | n_flight | n_gc | bal_acc | auroc | auprc | f1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| liver | grouped | 243 | 125 | 118 | 0.673 | 0.703 | 0.699 | 0.697 |
| liver | loo_accession | 243 | 125 | 118 | 0.737 | 0.822 | 0.831 | 0.695 |
| skeletal_muscle | grouped | 191 | 95 | 96 | 0.848 | 0.929 | 0.940 | 0.828 |
| skeletal_muscle | loo_accession | 191 | 95 | 96 | 0.853 | 0.950 | 0.955 | 0.841 |
| skin | grouped | 151 | 80 | 71 | 0.565 | 0.557 | 0.572 | 0.652 |
| skin | loo_accession | 151 | 80 | 71 | 0.585 | 0.617 | 0.606 | 0.667 |
| kidney | grouped | 135 | 68 | 67 | 0.534 | 0.489 | 0.502 | 0.479 |
| kidney | loo_accession | 135 | 68 | 67 | 0.518 | 0.521 | 0.532 | 0.564 |
| thymus | grouped | 117 | 63 | 54 | 0.553 | 0.706 | 0.735 | 0.671 |
| thymus | loo_accession | 117 | 63 | 54 | 0.556 | 0.678 | 0.732 | 0.662 |
| spleen | grouped | 109 | 55 | 54 | 0.457 | 0.492 | 0.484 | 0.528 |
| spleen | loo_accession | 109 | 55 | 54 | 0.458 | 0.524 | 0.517 | 0.512 |
| lung | grouped | 78 | 40 | 38 | 0.500 | 0.459 | 0.515 | 0.506 |
| lung | loo_accession | 78 | 40 | 38 | 0.482 | 0.448 | 0.501 | 0.583 |
| retina | grouped | 76 | 45 | 31 | 0.570 | 0.552 | 0.673 | 0.765 |
| retina | loo_accession | 76 | 45 | 31 | 0.516 | 0.566 | 0.619 | 0.750 |

Interpretation:

- Strong accession-aware signal: `skeletal_muscle`.
- Moderate accession-aware signal: `liver`.
- Weak or threshold-sensitive signal: `soleus` split, and possibly quadriceps
  by AUROC but not by 0.5-threshold balanced accuracy.
- No robust accession-aware classifier signal in this run: skin, kidney,
  thymus, spleen, lung, retina, EDL, gastrocnemius, tibialis anterior.

TabPFN3 does not produce FDR values here. These are predictive-validation
results, not differential-expression or pathway-enrichment significance tests.

## Skeletal-Muscle Splits

| dataset | cv | n | n_flight | n_gc | bal_acc | auroc | auprc | f1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| skeletal_muscle__soleus | grouped | 53 | 28 | 25 | 0.726 | 0.774 | 0.847 | 0.681 |
| skeletal_muscle__soleus | loo_accession | 53 | 28 | 25 | 0.775 | 0.811 | 0.877 | 0.778 |
| skeletal_muscle__gastrocnemius | grouped | 30 | 13 | 17 | 0.448 | 0.552 | 0.593 | 0.333 |
| skeletal_muscle__gastrocnemius | loo_accession | 30 | 13 | 17 | 0.536 | 0.561 | 0.610 | 0.381 |
| skeletal_muscle__quadriceps | grouped | 46 | 23 | 23 | 0.609 | 0.847 | 0.843 | 0.357 |
| skeletal_muscle__quadriceps | loo_accession | 46 | 23 | 23 | 0.587 | 0.735 | 0.759 | 0.457 |
| skeletal_muscle__edl | grouped | 32 | 16 | 16 | 0.531 | 0.406 | 0.539 | 0.118 |
| skeletal_muscle__edl | loo_accession | 32 | 16 | 16 | 0.562 | 0.516 | 0.659 | 0.222 |
| skeletal_muscle__tibialis_anterior | grouped | 30 | 15 | 15 | 0.433 | 0.511 | 0.491 | 0.320 |
| skeletal_muscle__tibialis_anterior | loo_accession | 30 | 15 | 15 | 0.367 | 0.373 | 0.448 | 0.174 |

The combined skeletal-muscle classifier is much stronger than most individual
muscle splits. Soleus is the clearest split-level signal. EDL and tibialis
anterior look perfect under random CV but fail accession-aware CV, which is a
strong warning that the random result is study/split leakage or very small-n
overfit.

## Random CV Check

Random CV balanced accuracy is much higher than accession-aware CV for several
small datasets:

| dataset | n | bal_acc | auroc | auprc |
| --- | ---: | ---: | ---: | ---: |
| liver | 243 | 0.856 | 0.923 | 0.932 |
| skeletal_muscle | 191 | 0.942 | 0.993 | 0.994 |
| skin | 151 | 0.607 | 0.685 | 0.714 |
| kidney | 135 | 0.667 | 0.751 | 0.751 |
| thymus | 117 | 0.862 | 0.922 | 0.941 |
| spleen | 109 | 0.688 | 0.798 | 0.810 |
| lung | 78 | 0.678 | 0.670 | 0.659 |
| retina | 76 | 0.647 | 0.763 | 0.831 |
| skeletal_muscle__soleus | 53 | 0.982 | 1.000 | 1.000 |
| skeletal_muscle__gastrocnemius | 30 | 0.767 | 0.851 | 0.851 |
| skeletal_muscle__quadriceps | 46 | 0.978 | 1.000 | 1.000 |
| skeletal_muscle__edl | 32 | 1.000 | 1.000 | 1.000 |
| skeletal_muscle__tibialis_anterior | 30 | 1.000 | 1.000 | 1.000 |

This confirms that random CV is too optimistic for biological claims. The
accession-aware results should be treated as primary.

## Feature Importance

Feature importance uses held-out permutation of a small candidate set:
`--importance-candidates 5` and `--permutation-repeats 1`. It is useful for
triage only, not as a stable gene-ranking analysis.

Top leave-one-accession-out genes from the `hvg` track:

| dataset | top genes as ENSMUSG:mean_drop |
|---|---|
| skeletal_muscle | ENSMUSG00000031762:0.017; ENSMUSG00000004939:0.009; ENSMUSG00000035963:0.008 |
| liver | ENSMUSG00000058672:0.057; ENSMUSG00000028715:0.039; ENSMUSG00000032080:0.022; ENSMUSG00000059824:0.017; ENSMUSG00000023067:0.016 |
| skeletal_muscle__soleus | ENSMUSG00000030541:0.050; ENSMUSG00000052276:0.050; ENSMUSG00000053025:0.042; ENSMUSG00000112071:0.042 |
| skeletal_muscle__quadriceps | ENSMUSG00000026822:0.042 |
| thymus | ENSMUSG00000031004:0.042; ENSMUSG00000073418:0.025; ENSMUSG00000027715:0.023 |
| skin | ENSMUSG00000104423:0.013; ENSMUSG00000020183:0.013; ENSMUSG00000027737:0.004 |

The repo does not currently contain a general mouse Ensembl-to-symbol mapping
for these OSDR features, so the production summaries report Ensembl IDs.

## Covariate-Augmented Check

A follow-up run added fold-local one-hot encoded OSDR design covariates to the same 500-gene `hvg` feature set. It used accession, tissue, muscle group, material type, sex, strain, genotype, platform, assay, data source, project identifier, and project type. It still excluded `condition_inferred`, `study.factor value.spaceflight`, sample names, profile IDs, and file names because those fields either are the target or can encode `FLT`/`GC` in the identifier text.

Output root:

```text
outputs/tabpfn3_osdr_covariates
```

Leave-one-accession-out comparison:

| dataset | expr LOO acc | +cov LOO acc | delta acc | expr LOO AUC | +cov LOO AUC | delta AUC |
|---|---:|---:|---:|---:|---:|---:|
| liver | 0.733 | 0.683 | -0.049 | 0.822 | 0.762 | -0.060 |
| skeletal_muscle | 0.853 | 0.906 | +0.052 | 0.950 | 0.950 | -0.001 |
| skin | 0.596 | 0.570 | -0.026 | 0.617 | 0.551 | -0.066 |
| kidney | 0.519 | 0.504 | -0.015 | 0.521 | 0.517 | -0.004 |
| thymus | 0.573 | 0.632 | +0.060 | 0.678 | 0.666 | -0.012 |
| spleen | 0.459 | 0.505 | +0.046 | 0.524 | 0.557 | +0.033 |
| lung | 0.487 | 0.462 | -0.026 | 0.448 | 0.443 | -0.005 |
| retina | 0.605 | 0.605 | +0.000 | 0.566 | 0.537 | -0.029 |
| skeletal_muscle__soleus | 0.774 | 0.774 | +0.000 | 0.811 | 0.823 | +0.011 |
| skeletal_muscle__gastrocnemius | 0.567 | 0.433 | -0.133 | 0.561 | 0.489 | -0.072 |
| skeletal_muscle__quadriceps | 0.587 | 0.717 | +0.130 | 0.735 | 0.767 | +0.032 |
| skeletal_muscle__edl | 0.562 | 0.562 | +0.000 | 0.516 | 0.516 | +0.000 |
| skeletal_muscle__tibialis_anterior | 0.367 | 0.367 | +0.000 | 0.373 | 0.373 | +0.000 |

Interpretation: metadata covariates do not broadly rescue weak tissues. They improve LOO accuracy for skeletal muscle and quadriceps and slightly improve spleen AUC, but liver, skin, lung, gastrocnemius, and retina get worse by AUC. This suggests most low grouped/LOO accuracy is not simply because the model lacked visible design covariates.

Feature importance was also run for the covariate model. The aggregated top-candidate table contains 501 rows, only 5 of which are covariate features. The top LOO features are still genes, not metadata, under the quick `--importance-candidates 5 --permutation-repeats 1` screen.

## Comparison To Earlier Methods

Compared with the expiMap/OntoVAE work, TabPFN3 gives the strongest practical
classification signal in skeletal muscle and recovers the expected broad
spaceflight sensitivity there more clearly than the pathway-level FDR analyses.
The tradeoff is interpretability: TabPFN3 is a discriminative classifier with
permutation gene importance, not a pathway-latent model with FDR-tested modules.

For biology follow-up, use TabPFN3 to prioritize tissues and candidate genes,
then validate with differential expression, pathway enrichment, and the existing
expiMap/OntoVAE pathway outputs.
