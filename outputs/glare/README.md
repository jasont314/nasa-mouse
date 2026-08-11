# Curated GLARE Results

This directory is the tracked GLARE result bundle for the NASA mouse analysis. It keeps the current interpretable summaries and validation tables, while large generated artifacts remain ignored under `outputs/`.

## What Is Included

- `summary/`: multi-tissue run summaries, local QC, and OSDR/TMS input audits.
- `validation/`: GLARE validation-stack summaries, candidate module tables, representation/clustering QC, and module direction tables.
- `validation/novel_candidate_validation/`: validated and near-validated complementary GLARE candidates.
- `validation/module_direction/`: signed FLT-vs-GC direction summaries for candidate modules.
- `validation/thymus_platelet_calcium_*`: focused thymus platelet/calcium/endothelial/stromal checks, including bulk marker adjustment and single-cell support.
- `study_effects/`: sample-level GLARE module-score PCA/UMAP audit summaries for accession/study separation.

Presentation plots are tracked separately in `presentation/glare/`.

## Current Interpretation

Use these results as a curated evidence layer, not as raw training output. The strongest GLARE-supported leads are:

- skeletal muscle circadian and Cyclin E/G1-S regulatory modules;
- soleus NGF signaling;
- kidney membrane trafficking and insulin receptor signaling;
- thymus platelet/calcium/endothelial/stromal remodeling, with composition sensitivity noted.

The full working tree under `outputs/glare/multi_tissue_api/` contains raw intermediate files, model outputs, and larger diagnostic tables. Those files are regenerable and intentionally not tracked.

## Cleanup Policy

Tracked `outputs/glare/*` directories and old liver MOBER intermediates were older exploratory runs or stale liver-focused artifacts. New GLARE outputs should stay under ignored `outputs/` until promoted into this curated results bundle.
