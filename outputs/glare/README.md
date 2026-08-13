# Curated GLARE Results

This directory is the tracked GLARE result bundle for the NASA mouse analysis. It keeps the current interpretable summaries and validation tables, while large generated artifacts remain ignored under `outputs/`.

## What Is Included

- `summary/`: multi-tissue run summaries, local QC, and OSDR/TMS input audits.
- `validation/`: GLARE validation-stack summaries, candidate module tables, representation/clustering QC, and module direction tables.
- `validation/novel_candidate_validation/`: validated and near-validated complementary GLARE candidates.
- `validation/module_direction/`: signed FLT-vs-GC direction summaries for candidate modules.
- `validation/thymus_platelet_calcium_*`: focused thymus platelet/calcium/endothelial/stromal checks, including bulk marker adjustment and single-cell support.
- `study_effects/`: sample-level GLARE module-score PCA/UMAP audit summaries for accession/study separation.

The GLARE panel used by the internship report is stored with that report's
source data.

## Final handoff interpretation

The internship report uses GLARE as a batch-effect and representation-learning
result. MOBER reduced accession separation in several tissues, but it did not
produce clear FLT/GC separation. The final biological claims therefore come
from the expiMap and synthetic-informed analyses.

The following GLARE-supported leads are preserved as exploratory records:

- skeletal muscle circadian and Cyclin E/G1-S regulatory modules;
- soleus NGF signaling;
- kidney membrane trafficking and insulin receptor signaling;
- thymus platelet/calcium/endothelial/stromal remodeling, with composition sensitivity noted.

They should not be read as part of the final biological claim set. The full
working tree under `outputs/glare/multi_tissue_api/` contains regenerable raw
intermediates, model outputs, and larger diagnostic tables that are not tracked.

## Cleanup Policy

The retained files document the final batch-effect summary and selected
exploratory follow-up. New GLARE outputs should remain ignored until they are
reviewed and deliberately added to this curated bundle.
