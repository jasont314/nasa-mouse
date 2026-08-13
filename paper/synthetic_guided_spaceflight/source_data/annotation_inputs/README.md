# Frozen annotation inputs

This directory contains compact candidate sets needed to verify the literature
annotation tables without restoring ignored analysis outputs.

`grouped_pathway_candidates.tsv` contains the ten grouped Reactome candidates
that entered literature review. It stores statistical and feature-importance
results only; the literature labels and rationale are added by
`annotate_importance_literature.py`.

Refresh this snapshot only after rerunning the grouped pathway analysis:

```bash
python \
  -m nasa_mouse_diffusion.paper_parity.annotate_importance_literature \
  --refresh-frozen-input
```

The annotation prompt and label definitions are recorded in
[`docs/annotation_prompts.md`](../../../../docs/annotation_prompts.md).
