# DDIM-guided FLT/GC feature workflow

## Purpose

This workflow asks whether a fixed conditional DDIM can improve FLT-versus-GC
feature discovery or classification without treating generated profiles as new
biological replicates. It complements the accession-held-out confirmation in
`outputs/generative/benchmark/analyses/generated_feature_guidance_confirmation_v1`.

The earlier per-tissue "best use" table is retained as an exploratory screen and
provenance record. It does not choose the final arm in this workflow.

## Data and fixed model

- Prepared OSDR panel: 974 mouse landmark genes.
- Development roles: 781 training and 536 validation profiles.
- Historical test role: 293 profiles, used only for descriptive reporting.
- Generator: ARCHS4-pretrained DDIM, then OSDR-fine-tuned with tissue, study,
  material-type, and FLT/GC conditioning.
- Synthetic views: three calibrated draws with seeds 3020, 3022, and 3023.
- Neural-network retraining in this workflow: none.

OSD-900 and OSD-457 are not used to tune this workflow. Their independent lung
and thymus confirmation remains the stronger result for those two accessions.

## Analysis design

For every eligible tissue, each repeat splits profiles inside
`accession x condition` strata. An inner split selects feature count,
regularization, and rank method; the matched outer split compares five arms:

1. `real_only`: real profiles and real-data feature ranking.
2. `generated_only`: generated profiles and generated-data feature ranking.
3. `real_plus_generated`: real plus generated profiles at 1:1 total weight.
4. `guided_real_only`: real classifier with generated-informed feature ranking.
5. `guided_low_weight`: real profiles plus recentered generated profiles at 0.05
   total weight.

The selected generated-informed arm must match or improve real-only mean balanced
accuracy, AUROC, and average precision. No composite score is used. Repeat-level
nonworse and strict-win rates are reported separately because the eight outer
splits overlap and are not independent observations.

Accessions remain on both sides of the primary classifier splits. This measures
within-study interpolation, as requested. Accession-aware random-effects and LOO
analyses are then applied to real data only to test whether selected gene
directions are supported across studies.

## Stable gene sets

Genes must be selected in at least 50% of repeats and have classifier-coefficient
sign agreement of at least 75%.

- `core_intersection`: stable in real and generated arms, matching coefficient
  directions, and supported by the real accession-level direction.
- `generated_supported`: stable in the generated-informed arm and supported by
  the real accession-level direction, but not in the core intersection.
- `exploratory_union`: all other genes stable in either arm.

The union is for hypothesis generation. A gene is not a biological finding merely
because it was generated-supported; real random-effects FDR and LOO stability are
separate requirements.

## Current results

Twenty-two of 24 tissue categories completed eight repeats. `cells` lacked valid
nested strata and `cultured_cells` had only four development profiles.

| Tissue | Selected arm | Mean delta BA | Mean delta AUROC | Mean delta AP | All metrics nonworse | Stable genes | Stable genes LOO FDR |
|---|---|---:|---:|---:|---:|---:|---:|
| kidney | guided low weight | +0.029 | +0.093 | +0.097 | 6/8 | 39 | 0 |
| liver | guided low weight | +0.024 | +0.054 | +0.044 | 5/8 | 54 | 0 |
| lung | generated only | +0.086 | +0.156 | +0.157 | 7/8 | 50 | 0 |
| retina | guided low weight | +0.117 | +0.121 | +0.071 | 7/8 | 47 | 0 |
| skeletal muscle | real plus generated | +0.043 | +0.000 | +0.013 | 5/8 | 97 | 7 |
| skin | generated only | +0.086 | +0.109 | +0.065 | 4/8 | 60 | 0 |
| spleen | guided real only | +0.170 | +0.208 | +0.204 | 7/8 | 57 | 1 |
| thymus | generated only | +0.121 | +0.092 | +0.075 | 7/8 | 62 | 8 |

The generated-informed arm is an exploratory per-tissue choice. In particular,
the independently confirmed thymus policy remains a real-only classifier with
generated-guided features, and the independently tested lung policy remains a
low-weight guided model. The nested screen must not be used to retune either
policy against OSD-457 or OSD-900.

### Evidence interpretation

- **Thymus:** strongest cross-method support. Eight selected genes pass the real
  LOO-FDR rule: `Cenpe`, `Ccnb1`, `Nusap1`, `Stmn1`, `Cdk1`, `Top2a`, `Ccnb2`,
  and `Ccne2`. Twenty-five generated-supported Reactome terms pass FDR 0.05 and
  are dominated by mitosis, G2/M checkpoints, and nuclear-envelope reformation.
  This agrees with the independent OSD-457 confirmation.
- **Skeletal muscle:** seven selected genes pass real LOO FDR: `Sox4`, `Sh3bp5`,
  `Cebpd`, `Cdkn1a`, `Bphl`, `Prkcd`, and `Arid5b`. Core pathways include TP53
  metabolic regulation and ERBB2 signaling. This is a strong follow-up signal,
  although the selected classifier arm was nonworse on all metrics in only 5/8
  overlapping repeats.
- **Spleen:** large nested classifier gains and one LOO-stable selected gene,
  `Igfbp3`, but no stable-set Reactome enrichment at FDR 0.05.
- **Lung:** nested prediction improves, and a core heat-shock/chaperone pathway
  enrichment passes FDR. No selected gene passes real random-effects FDR or LOO,
  consistent with the prior mixed genotype-stratified OSD-900 result.
- **Kidney, liver, retina, and skin:** useful exploratory feature/pathway screens,
  but no selected gene passes the real LOO-FDR rule. Retina and skin pathway
  enrichments therefore remain hypotheses, not validated pathway claims.

The 292 significant Reactome rows across all tissues and nested gene sets are not
292 independent discoveries. Reactome is hierarchical, terms overlap, and the
three reported gene sets are nested for enrichment.

## Independence limitations

The fixed DDIM was fine-tuned before these nested splits and saw the original
training role. Consequently, repeated-split performance can be optimistic and is
not an independent estimate of generator generalization. The historical test role
was also examined in earlier work, so `descriptive_original_test.tsv` is a sanity
check only. It is not a second confirmation set.

Generated profiles are never entered into random-effects meta-analysis as if they
were animals. Biological support comes from real within-accession FLT-minus-GC
effects. Entire-accession holdout remains necessary for claims about unseen-study
generalization, even though it is not the primary split used here.

## Run

```bash
PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python \
  -m nasa_mouse_diffusion.paper_parity.within_study_feature_stability \
  --config configs/generative/diffusion/within_study_generated_feature_stability.yaml
```

Outputs are under
`outputs/generative/benchmark/analyses/within_study_generated_feature_stability_v1`.
The main entry points are `tissue_arm_choices.tsv`, `paired_repeat_support.tsv`,
`stable_gene_sets.tsv.gz`, `biological_support_summary.tsv`,
`real_random_effects.tsv.gz`, and `reactome_enrichment.tsv.gz`.
