# Synthetic-guided spaceflight paper package

This package reports a frozen benchmark of conditional bulk RNA-seq generation
and its use for mouse spaceflight feature discovery. The central result is an
independently held-out thymus confirmation. The soleus result is cross-accession
developmental evidence and is deliberately labeled as requiring a new unseen
accession.

Generated profiles are not treated as additional animals. Biological support is
calculated from real NASA OSDR profiles using within-accession effects,
random-effects meta-analysis, FDR, and leave-one-accession-out sensitivity.

## Contents

- `manuscript.md`, `manuscript.html`, and `manuscript.pdf`: full paper draft.
- `supplementary_methods.md`, `.html`, and `.pdf`: exact implementation,
  evaluation gates, output provenance, and limitations.
- `figures/`: six main figures and three supplementary figures in PNG/PDF.
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

## Evidence labels

- **Independent confirmation:** the generator and feature policy did not use the
  test accession. This applies to OSD-457 thymus; lung OSD-900 is mixed.
- **Cross-accession development:** real random-effects and leave-one-accession-out
  evidence are available, but the generator was adapted in the same development
  domain. This applies to soleus.
- **Exploratory:** predictive or enrichment signals without a stable real
  accession-aware gene set.
- **Negative:** no coherent synthetic-guided biological result under the declared
  rules.

