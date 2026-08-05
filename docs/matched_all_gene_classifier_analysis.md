# Matched all-gene classifier analysis

## Question

The consensus workflow ranks genes before classifier fitting. It is useful for
finding small, repeatable panels, but it does not directly isolate the effect of
the training source. A synthetic arm can use a different gene set, feature
count, or regularization value from the real-only arm.

This analysis holds those choices fixed. It asks whether a classifier trained
on synthetic profiles, or on real and synthetic profiles together, predicts
held-out real FLT and GC profiles better than the same classifier trained on
real profiles alone. It also compares gene importance across the three arms.

## Design

The analysis covers 22 canonical tissues and five skeletal-muscle groups. Each
classifier uses all 974 expression landmarks.

| Arm | Training profiles | Role in the analysis |
|---|---|---|
| Real only | Real OSDR profiles | Reference classifier |
| Synthetic only | Three matched DDIM draws | Test whether the generated distribution transfers to real profiles |
| Real + synthetic | Real profiles plus three DDIM draws | Test synthetic augmentation |

Eight nested splits are used per analysis unit. The scaler is fitted on the real
outer-training profiles. Ridge logistic regularization is chosen from inner
real-only data, then the same value is used for every arm in that split. The
real-plus-synthetic arm gives the complete synthetic set the same total weight
as the real set.

Every arm is scored on the same outer real profiles. Performance is reported
both after pooling those profiles and after calculating each accession
separately and averaging the accession scores. An arm passes the joint utility
gate when its mean balanced accuracy, AUROC, and average precision are all no
worse than real only in both summaries.

Gene importance is calculated in two ways:

1. Accession-blocked permutation importance measures the loss in classifier
   performance when one gene is shuffled. Ten permutations are used for every
   gene and fit. Synthetic permutations also remain separate by DDIM draw.
2. Exact linear SHAP assigns directional log-odds contributions relative to the
   shared real-training background.

A gene has repeat-consistent marginal importance when its mean held-out-real
AUROC permutation loss is at least `0.001` and that loss is positive in at least
half of the outer splits. The practical `0.001` cutoff excludes numerical
changes near zero.

The biological candidate table starts from real-data random-effects BH-FDR
associations. Synthetic profiles do not enter the BH-FDR calculation. A row is
retained only when the synthetic arm passes the joint utility gate, has
repeat-consistent importance on held-out real profiles, has a coefficient that
matches the real FLT versus GC effect direction, and has positive FLT versus GC
SHAP separation.

## Classifier utility

The run fitted 648 classifiers: 27 units, eight outer splits, and three arms.
Real plus synthetic was no worse on all six mean metrics in 18 units. Sixteen of
those units improved at least one metric; bone marrow and brown adipose tissue
were exact ties. Synthetic only passed the same gate in six units. Five improved
at least one metric, while brown adipose tissue tied.

| Analysis unit | Arm | Pooled BA change | Pooled AUROC change | Pooled AP change | Accession-macro BA change | Accession-macro AUROC change | Accession-macro AP change |
|---|---|---:|---:|---:|---:|---:|---:|
| Eye | Real + synthetic | +0.094 | +0.156 | +0.115 | +0.094 | +0.062 | +0.031 |
| Retina | Real + synthetic | +0.100 | +0.075 | +0.063 | +0.105 | +0.117 | +0.086 |
| Lung | Real + synthetic | +0.062 | +0.086 | +0.108 | +0.068 | +0.062 | +0.054 |
| Skin | Real + synthetic | +0.074 | +0.087 | +0.076 | +0.096 | +0.105 | +0.059 |
| Thymus | Real + synthetic | +0.061 | +0.046 | +0.042 | +0.066 | +0.028 | +0.021 |
| Spleen | Real + synthetic | +0.057 | +0.075 | +0.072 | +0.042 | +0.028 | +0.014 |
| Liver | Real + synthetic | +0.043 | +0.037 | +0.029 | +0.044 | +0.045 | +0.025 |
| Skeletal muscle, pooled | Real + synthetic | +0.010 | +0.011 | +0.014 | +0.012 | +0.026 | +0.020 |

Synthetic-only performance was most convincing in eye, lung, skin, and thymus.
It was not uniformly reliable. In pooled skeletal muscle, for example,
synthetic-only balanced accuracy was `0.670`, compared with `0.952` for real
only, even though its AUROC was high. This makes synthetic only a useful test of
generator transfer, but not a default replacement for real training data.

## BH-FDR gene results

Twenty-one unique tissue-gene associations pass the complete matched gate. Two
associations are supported by both synthetic arms, producing 23 arm-level rows.
The results occur in four tissues.

| Tissue | FLT higher | FLT lower | Unique associations |
|---|---|---|---:|
| Thymus | `Klhdc2`, `Snx7`, `Etv1`, `Plscr1`, `Tspan3`, `Socs2` | `Nusap1`, `Stmn1`, `Birc5`, `Ccnb2`, `E2f2`, `Ube2c`, `Cdc20`, `Gmnn`, `Kif20a` | 15 |
| Liver | None | `Grb10`, `Ppic`, `H2-DMa`, `Gtf2a2` | 4 |
| Skin | `Plscr1` | None | 1 |
| Spleen | `Loxl1` | None | 1 |

Seven thymus genes are promoted under the matched definition: `Klhdc2`, `Snx7`,
`Etv1`, `Plscr1`, `Tspan3`, `Socs2`, and `Kif20a`. They cross the marginal
importance gate in the real-plus-synthetic classifier but not in the real-only
classifier. The other associations retain real-only importance after synthetic
training, with some increasing and some decreasing in magnitude.

`Ccnb2` in thymus and `Plscr1` in skin are supported by both the synthetic-only
and real-plus-synthetic classifiers. No skeletal-muscle group passes the full
matched BH-FDR and utility gate. That does not remove the real-data muscle
associations; it means this all-gene analysis does not show added gene-level
support from synthetic training for those associations.

## Reactome results

Reactome enrichment used the 974 landmark genes as background. The 15 thymus
candidates tested 36 terms, and 26 had FDR below 0.05. The leading term was
mitotic cell cycle (`R-MMU-69278`, FDR `0.004744`), supported by `Ube2c`,
`Kif20a`, `Cdc20`, `Gmnn`, `E2f2`, and `Ccnb2`. The flight-lower thymus subset
had 33 significant terms. The four liver candidates tested two terms, neither
significant. Single skin and spleen candidates were not treated as gene sets.

These enrichment rows are overlapping Reactome parent and child terms, not 26
independent discoveries. They show that the matched thymus genes form a
cell-cycle program; the narrower liver, skin, and spleen results do not yet form
comparable pathway-level claims.

The crosswalk joins annotations by Ensembl gene ID. It keeps both symbols when
the BH inventory and expression annotation use different aliases. This occurs
for `ENSMUSG00000061559`, labeled `Skic8` in the BH inventory and `Wdr61` in the
expression annotation.

## Relation to consensus ranking

The earlier consensus analysis reported 49 synthetic-informed BH-FDR
associations. Eleven also pass this matched all-gene analysis. Thirty-eight are
specific to the consensus workflow, while ten associations are specific to the
matched workflow.

The difference is expected because the methods estimate different quantities.
Consensus ranking asks which genes repeatedly enter a small predictive panel.
The matched workflow asks whether one gene has measurable marginal importance
when all 974 genes are available. Correlated genes can replace one another in
an all-gene ridge classifier, reducing the permutation loss for any single gene.
The matched workflow also requires the complete synthetic arm to preserve all
pooled and accession-macro metrics.

For claims about whether synthetic training changes a classifier or transfers
gene importance to real profiles, the matched all-gene workflow should be the
primary analysis. The consensus workflow remains useful as a secondary analysis
of compact, correlated gene panels. Real-data BH-FDR remains the association
test in both cases.

## Limits

The outer profiles are held out from classifier fitting and regularization
selection, but this is not an independent test of the generator. The fixed DDIM
was trained before these classifier splits and saw the broader development
pool. The outer split also retains represented accessions; it is not a held-out
study test.

Repeated outer splits overlap. Some small units have only four profiles in an
outer test split. Marginal permutation importance can understate the value of
correlated genes. SHAP explains the fitted classifier and is not an independent
biological test. Generated profiles are model draws, not biological replicates.

## Outputs

Results are under
`outputs/generative_benchmark/analyses/matched_all_gene_classifiers_osdr_disjoint_v1/`.

- `eligible_bh_fdr_candidates.tsv` contains the 23 retained arm-level rows.
- `bh_fdr_matched_importance.tsv.gz` contains all 459 BH-FDR associations for
  both synthetic arms.
- `arm_utility.tsv` contains paired performance changes by unit and arm.
- `arm_gene_comparison.tsv.gz` contains matched importance for all 974 genes.
- `importance_summary.tsv.gz` contains arm and domain summaries.
- `importance_by_repeat.tsv.gz` contains repeat-level importance values.
- `matched_classifier_metric_deltas.png` summarizes tissue-level performance.
- `<scope>/<unit>/matched_classifier_importance.png` shows each unit's
  importance comparison and top-gene heatmaps.

Run the workflow with:

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
NUMEXPR_NUM_THREADS=1 PYTHONPATH=src \
conda run --no-capture-output -n nasa-mouse python \
  -m nasa_mouse_rna_diffusion.matched_all_gene_classifiers \
  --config configs/rna_diffusion/matched_all_gene_classifiers_osdr_disjoint.yaml
```
