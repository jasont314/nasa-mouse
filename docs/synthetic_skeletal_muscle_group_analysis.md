# Synthetic-Guided Skeletal-Muscle Group Analysis

Date: 2026-07-28

## Question

Does the accepted ARCHS4-pretrained, OSDR-adapted DDIM provide useful information
when skeletal muscle is separated into anatomical groups rather than pooled?

## Design

The analysis reuses the fixed accepted DDIM and three independently sampled,
calibrated synthetic development views. No neural network was retrained. For each
muscle group, five downstream uses were compared over eight repeated nested
within-accession splits:

- real only;
- generated only;
- real plus generated;
- synthetic-guided feature selection with a real-only classifier;
- synthetic-guided feature selection with 5% synthetic training weight.

The original test role was not used for arm or feature selection. Biological
support was assessed afterward using only real FLT and ground-control profiles,
accession-level random-effects estimates, and leave-one-accession-out (LOO)
sensitivity checks. Reactome tests are enrichment tests on selected genes, not
independent measurements of pathway activity.

| Muscle group | Development profiles | Accessions | FLT | GC |
|---|---:|---:|---:|---:|
| EDL | 24 | 2 | 12 | 12 |
| Gastrocnemius | 25 | 3 | 10 | 15 |
| Quadriceps | 35 | 4 | 18 | 17 |
| Soleus | 41 | 3 | 22 | 19 |
| Tibialis anterior | 24 | 2 | 12 | 12 |

## Results

| Group | Selected synthetic use | Mean BA / AUC / AP change | Repeats all metrics nonworse | Mean pooled / accession effect correlation | Synthetic-selected genes with real LOO FDR < 0.05 |
|---|---|---|---:|---|---:|
| EDL | Generated only | +0.021 / 0.000 / 0.000 | 6/8 | 0.970 / 0.959 | 0 |
| Gastrocnemius | Guided, 5% synthetic weight | +0.063 / +0.042 / +0.009 | 5/8 | 0.978 / 0.762 | 0 |
| Quadriceps | Guided, 5% synthetic weight | +0.050 / +0.040 / +0.026 | 5/8 | 0.285 / 0.708 | 1 |
| Soleus | Generated only | +0.025 / +0.020 / +0.020 | 6/8 | 0.955 / 0.893 | 7 |
| Tibialis anterior | Real plus generated | 0.000 / 0.000 / 0.000 | 8/8 | 0.726 / 0.860 | 0 |

The fidelity correlations are internal diagnostics. Calibration and DDIM
adaptation used the development domain, so they do not establish unseen-study
generalization.

### Soleus

Soleus is the strongest grouped result. Seven synthetic-selected genes were also
significant and directionally stable in every real-data LOO refit:

| Gene | Real FLT minus GC effect | Direction |
|---|---:|---|
| `Bdh1` | -0.0645 | Lower in flight |
| `Bnip3` | -0.0358 | Lower in flight |
| `Mef2c` | -0.0335 | Lower in flight |
| `Ech1` | -0.0298 | Lower in flight |
| `Pxmp2` | -0.0191 | Lower in flight |
| `Gmnn` | -0.0011 | Lower in flight |
| `Tpm1` | +0.0260 | Higher in flight |

`Arid5b` was an additional real-only stable gene. The selected sets were enriched
for mitochondrial fatty-acid beta oxidation, lipid metabolism, and mitochondrial
protein-turnover terms. The pathway rows are hierarchical duplicates rather than
independent discoveries.

The clearest interpretation is reduced soleus oxidative fuel handling with
concurrent changes in mitochondrial quality control, muscle transcription, and
contractile structure. This agrees with the earlier count-level result in
`docs/expimap_skeletal_muscle_prior_work.md`, where soleus fatty-acid oxidation
was lower in flight but the broad module did not pass the stricter LOO FDR rule.
Synthetic guidance therefore helped localize that broad signal to a smaller set
of real-data-stable genes.

### Other Muscle Groups

- **Quadriceps:** `Rbm6` was the only synthetic-selected, real LOO-stable gene.
  Pooled synthetic effect recovery varied substantially by draw, and no selected
  Reactome family passed FDR 0.05. This is not yet a coherent biological result.
- **EDL:** flight-lower `Abcc5`, `Lsm6`, `Polr2i`, and `Tsc22d3` were consistent
  across both accessions. RNA processing and nuclear-receptor terms were enriched,
  but the terms were driven by small overlaps and cannot be LOO-confirmed with
  only two accessions.
- **Tibialis anterior:** `Cdkn1a`, `St3gal5`, `Cebpd`, `Pdhx`, and `Bnip3` were
  higher in flight in both accessions and passed ordinary meta FDR 0.05. The
  synthetic arm tied an already-perfect real-only classifier, no Reactome family
  passed FDR, and two accessions are insufficient for a strong LOO claim.
- **Gastrocnemius:** downstream metrics improved modestly, but no selected gene
  passed real LOO FDR and no Reactome family passed FDR. Treat this as negative or
  exploratory.

## Conclusion

Separating skeletal muscle materially improves the synthetic-data interpretation.
The pooled muscle result concealed a coherent soleus-specific signal. The
defensible claim is that synthetic-guided feature selection prioritizes a
real-data-stable soleus program involving lower mitochondrial lipid oxidation and
associated muscle-regulatory changes in flight.

This remains within-study evidence. A new soleus accession that was excluded from
DDIM adaptation and all feature selection is required for independent
confirmation. Synthetic profiles must not be counted as additional animals.

## Reproduction

```bash
PYTHONPATH=src python -m nasa_mouse_diffusion.paper_parity.within_study_feature_stability \
  --config configs/generative/diffusion/within_study_generated_feature_stability_muscle_groups.yaml
```

Outputs are written to:

`outputs/generative/benchmark/analyses/within_study_generated_feature_stability_muscle_groups_v1/`

Key files:

- `tissue_arm_choices.tsv`
- `paired_repeat_support.tsv`
- `synthetic_group_validation.tsv`
- `stable_gene_sets.tsv.gz`
- `real_random_effects.tsv.gz`
- `reactome_enrichment.tsv.gz`
- `biological_support_summary.tsv`
