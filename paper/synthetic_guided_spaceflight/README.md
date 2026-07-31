# Synthetic-guided spaceflight paper package

This package contains a research manuscript on organ-specific mouse spaceflight
responses identified with synthetic-guided transcriptomics. The central result is
a leakage-corrected retrospective thymus test. The OSD-457 accession was excluded
from OSDR adaptation, and all OSDR-linked GEO series were excluded from ARCHS4
pretraining before the model was retrained from scratch. Soleus provides
cross-accession metabolic evidence, while kidney provides a focused promoted
`Inpp4b` and reinforced `Slc37a4` result.

The corrected OSDR-disjoint backbone was used for the lung/thymus held-out
experiment and for regenerated all-tissue and muscle-group development screens.
All BH-FDR effects use real OSDR samples only. Synthetic attribution is retained
only where the selected generated arm passed the balanced-accuracy, AUROC, and
average-precision eligibility gate.

Generated profiles are not treated as additional animals. Biological support is
calculated from real NASA OSDR profiles using within-accession effects,
random-effects meta-analysis, and Benjamini-Hochberg FDR. BH FDR below 0.05 is
the primary gene-level inclusion rule. Synthetic-selection status, study-direction
agreement, heterogeneity, and leave-one-accession-out results are reported as
separate interpretation or sensitivity annotations rather than inclusion gates.

## Contents

- `manuscript.md`, `manuscript.html`, and `manuscript.pdf`: full paper draft.
- `supplementary_methods.md`, `.html`, and `.pdf`: exact implementation,
  evaluation gates, output provenance, and limitations.
- `figures/`: four biology-focused main figures and five supplementary figures
  in PNG/PDF. Model-performance panels are retained as Figures S4-S5.
- `source_data/`: manuscript tables, figure source tables, and SHA-256 manifests
  for every frozen analysis input.

## Rebuild

The build only reads completed outputs. It does not retrain any model or rerun
feature selection.

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m nasa_mouse_rna_diffusion.build_synthetic_guided_paper
```

The builder checks key frozen values before writing figures, then renders the
Markdown manuscripts with WeasyPrint.

## Analysis tiers

- **Tier 1, leakage-corrected held-out study test:** a fixed generator and feature
  policy were tested in an OSDR accession excluded from adaptation and policy
  development after removing every OSDR-linked GEO series from ARCHS4. The
  eight-gene OSD-457 thymus panel is the central result. Because outcomes were
  seen in the earlier overlapping analysis, a new untouched study is still
  required for prospective confirmation.
- **Tier 2, synthetic-informed development:** real-data BH-FDR genes were also
  repeatedly synthetic-promoted or reinforced. This contains 49 tissue-gene
  results and supports developmental hypotheses rather than independent transfer.
- **Tier 3, complete real-data screen:** all 459 random-effects BH-FDR
  tissue-gene results, regardless of synthetic feature-selection status.

## Tissue evidence labels

- **Cross-accession development:** selected BH-FDR genes form a coherent pattern
  across represented studies, but the generator was adapted in the same
  development domain. This applies most clearly to soleus, pooled skeletal
  muscle, and kidney.
- **Exploratory:** predictive or enrichment signals without a stable real
  accession-aware gene set or without unseen-study transfer. This applies to
  spleen, skin, adrenal gland, gastrocnemius, and tibialis anterior.
- **Negative:** no coherent synthetic-guided biological result under the declared
  rules. Quadriceps, EDL, and liver retained real-only arms; lung failed its
  accession-held-out test.
