# Classifier feature importance

## Purpose

The synthetic-guided analysis originally used consensus ranks to choose genes and
logistic-regression coefficients to show classifier direction. Those quantities
are not post-selection feature importance. Consensus measures how repeatedly a
gene is selected; a coefficient records its fitted weight.

This follow-up asks a different question: once a classifier and its gene set have
been selected, how much does each gene contribute to FLT/GC prediction?

## Analysis

The analysis reconstructs the five classifier arms from their saved nested-run
hyperparameters and selected genes. It covers 22 canonical tissues and five
skeletal-muscle groups. No DDIM is retrained and no classifier hyperparameter is
reselected.

For every tissue, repeat, and arm:

1. The classifier is refit on the original outer-training profiles.
2. Its saved feature coefficients are restored.
3. Held-out metrics are checked against the original nested-run metrics.
4. Each selected gene is permuted 20 times on held-out real profiles and matched
   generated profiles. Synthetic permutations remain separate within each DDIM
   draw.
5. The loss in balanced accuracy, AUROC, and average precision is recorded.
6. Exact linear SHAP contributions are calculated on the log-odds scale as
   `coefficient * (expression - training background)`.

The 1,080 reconstructed classifiers reproduced their saved held-out metrics to a
maximum absolute error of `1.11e-16`. The final output contains 82,600
repeat-level gene/domain rows.

## Main Result

The manuscript contains 49 tissue-gene associations that both pass real-data BH
FDR and were promoted or reinforced by the selected synthetic-informed arm.
Thirty-three of those 49 have a positive mean held-out-real AUROC permutation
loss in at least half of the repeats where they were selected. All 49 have a
positive mean linear SHAP FLT-minus-GC contribution.

| Analysis unit | Promoted | Reinforced | Total | Positive permutation |
|---|---:|---:|---:|---:|
| Adrenal gland | 1 | 1 | 2 | 1 |
| Eye | 0 | 1 | 1 | 0 |
| Kidney | 1 | 1 | 2 | 2 |
| Skeletal muscle, pooled | 4 | 8 | 12 | 10 |
| Skin | 1 | 0 | 1 | 1 |
| Spleen | 3 | 1 | 4 | 2 |
| Thymus | 13 | 3 | 16 | 14 |
| Gastrocnemius | 2 | 0 | 2 | 1 |
| Soleus | 0 | 5 | 5 | 2 |
| Tibialis anterior | 1 | 3 | 4 | 0 |
| **Total** | **26** | **23** | **49** | **33** |

Among promoted genes, 19 of 26 have positive permutation support. Among
reinforced genes, 14 of 23 do. The clearest coherent panel remains thymus: 14 of
16 genes have repeat-consistent marginal importance, including `Stmn1`, `Ube2c`,
`Ccne2`, `Nusap1`, `Birc5`, `Gmnn`, `Ccnb2`, and `Cdk1`. Pooled skeletal muscle
retains support for 10 of 12 genes, and both kidney genes (`Slc37a4` and `Inpp4b`)
are supported.

Soleus illustrates why permutation and consensus answer different questions.
`Tpm1` and `Bdh1` have clear marginal importance, while `Ech1`, `Bnip3`, and
`Decr1` have positive SHAP direction but little single-gene permutation loss.
Those genes belong to a correlated metabolic panel, so the remaining genes can
partly replace one permuted feature. Lack of marginal permutation loss therefore
does not negate their real-data BH-FDR association or repeated selection.

The same distinction applies to spleen, where `Loxl1` and `Myl9` have positive
permutation support but `Rai14` and `Ptprk` do not meet the repeat-consistency
rule. Tibialis-anterior genes remain directionally aligned in SHAP but do not
show repeat-consistent marginal permutation importance.

## Interpretation

- **Consensus rank** determines which genes enter a classifier repeatedly.
- **Coefficient** gives the direction and scale of the fitted linear weight.
- **Permutation importance** measures marginal predictive dependence on a gene.
- **Linear SHAP** assigns each gene a directional contribution to sample
  log-odds relative to the training background.

Permutation magnitudes should not be compared directly between tissues because
sample counts, selected feature counts, and baseline classifier performance
differ. Correlated genes can divide or mask permutation importance. SHAP values
also explain the fitted classifier; they are not an independent biological test.
The real-data random-effects BH-FDR result remains the biological association
test.

Synthetic-domain and held-out-real importance rankings were only weakly
correlated across arms (median Spearman correlations about `0.13-0.20`). This is
why synthetic-only importance is not treated as biological evidence. The
held-out-real permutation columns are the relevant post-selection check.
"Held out" here refers to the classifier's outer split. The fixed DDIM predates
the nested classifier splits and is not independently validated by this test.

## Outputs

The complete analysis is under
`outputs/generative_benchmark/analyses/classifier_importance_osdr_disjoint_v1/`.
Key files are:

- `synthetic_informed_bh_fdr_gene_importance.tsv`: 49-gene manuscript crosswalk.
- `synthetic_informed_gene_importance.png`: permutation and SHAP comparison.
- `selected_arm_importance.tsv.gz`: retained-arm importance for every unit.
- `importance_summary.tsv.gz`: all arms, genes, tissues, and domains.
- `importance_by_repeat.tsv.gz`: repeat-level values.
- `<scope>/<unit>/classifier_importance_heatmaps.png`: per-unit arm comparison.

Reproduce the analysis with:

```bash
PYTHONPATH=src conda run -n nasa-mouse python \
  -m nasa_mouse_rna_diffusion.classifier_importance \
  --config configs/rna_diffusion/classifier_importance_osdr_disjoint.yaml
```
