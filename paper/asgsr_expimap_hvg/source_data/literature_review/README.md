# expiMap literature annotations

This directory contains the source-checked pathway annotations used by the
expiMap paper. OpenAI Codex (GPT-5) assisted with organizing literature and
drafting concise labels. A source-level review was required before a label was
kept. The annotations help interpret model results; they did not calculate
pathway effects or determine statistical significance.

The repository stores the auditable rationale rather than a model chat
transcript. Each reviewed row records the observed direction, literature
relationship, direction assessment, confidence, rationale, and citation keys.

- [`manual/`](manual/) contains the review decisions entered for each tissue.
- [`sources/`](sources/) maps citation keys to titles, DOI or PMID records, and
  URLs.
- [`final/`](final/) merges the review fields with observed expiMap pathway
  effects. The four files contain 1,427 reviewed pathways.
- [`../table_2_retained_pathway_evidence.tsv`](../table_2_retained_pathway_evidence.tsv)
  contains the 16 pathways retained for the main and secondary paper results.
- [`../table_s30_kidney_spleen_literature_sources.tsv`](../table_s30_kidney_spleen_literature_sources.tsv)
  contains the source catalog for the corrected kidney and spleen review.

The project-wide annotation record, including the synthetic feature-importance
review, is in
[`docs/annotation_provenance.md`](../../../../docs/annotation_provenance.md).

Run the following command after refreshing the selected expiMap analyses:

```bash
python -m nasa_mouse_expimap.plot_hvg_literature_review_heatmaps
```

It updates `final/` and writes review heatmaps and diagnostics to
`outputs/expimap/analyses/literature_review/`.
