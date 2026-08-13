# Literature annotation record

The project used LLM assistance to organize literature for selected expiMap
pathways and synthetic-analysis genes and pathways. Annotation happened after
the quantitative candidate sets were fixed. It did not calculate expression
effects, FDR values, classifier scores, or feature importance.

The tracked record contains concise, reviewable rationale rather than a model's
private chain of thought. Each annotation can be traced from its result row to a
source identifier and then to a citation, DOI, and URL. Labels organize the
interpretation; they are not independent biological evidence.

## expiMap pathway annotations

OpenAI Codex (GPT-5) assisted with the expiMap literature review. Source-level
review was required before an annotation was used in a figure, table, or paper.
The review uses three paper-facing roles: literature aligned, complementary,
and context sensitive.

| Record | Contents |
|---|---|
| [`paper/asgsr_expimap_hvg/source_data/literature_review/manual/`](../paper/asgsr_expimap_hvg/source_data/literature_review/manual/) | Per-pathway classification, direction assessment, confidence, rationale, and citation keys |
| [`paper/asgsr_expimap_hvg/source_data/literature_review/final/`](../paper/asgsr_expimap_hvg/source_data/literature_review/final/) | The same annotations merged with 1,427 observed pathway-effect rows |
| [`paper/asgsr_expimap_hvg/source_data/literature_review/sources/`](../paper/asgsr_expimap_hvg/source_data/literature_review/sources/) | Full source records for thymus, skin, liver, and soleus |
| [`paper/asgsr_expimap_hvg/source_data/table_2_retained_pathway_evidence.tsv`](../paper/asgsr_expimap_hvg/source_data/table_2_retained_pathway_evidence.tsv) | Rationale, literature keys, and quantitative evidence for 16 retained pathways |
| [`paper/asgsr_expimap_hvg/source_data/table_s30_kidney_spleen_literature_sources.tsv`](../paper/asgsr_expimap_hvg/source_data/table_s30_kidney_spleen_literature_sources.tsv) | Full source records for the corrected kidney and spleen review |
| [`paper/asgsr_expimap_hvg/source_data/table_s34_retained_pathway_member_gene_support.tsv`](../paper/asgsr_expimap_hvg/source_data/table_s34_retained_pathway_member_gene_support.tsv) | Member-gene support used to check pathway-level interpretations |

The `manual` and `final` tables expose the editorial reasoning through
`literature_alignment`, `direction_assessment`, `confidence`, `rationale` or
`review_rationale`, and `citations`. The retained-pathway table uses
`manual_rationale` and `literature_keys`.

## Synthetic feature-importance annotations

The generative analysis keeps selection behavior separate from literature
interpretation. `Promoted` and `reinforced` describe what changed after adding
synthetic training data. `Aligning`, `complementary`, `ambiguous`, and
`unmatched` describe the relationship to prior literature.

| Record | Contents |
|---|---|
| [`paper/synthetic_guided_spaceflight/source_data/table_s16_promoted_gene_literature_annotations.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_s16_promoted_gene_literature_annotations.tsv) | 49 consensus-ranking gene annotations |
| [`paper/synthetic_guided_spaceflight/source_data/table_s17_promoted_gene_literature_sources.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_s17_promoted_gene_literature_sources.tsv) | Source catalog for the consensus annotations |
| [`paper/synthetic_guided_spaceflight/source_data/table_s22_matched_gene_literature_annotations.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_s22_matched_gene_literature_annotations.tsv) | 21 individual-gene permutation and SHAP annotations |
| [`paper/synthetic_guided_spaceflight/source_data/table_s23_grouped_pathway_literature_annotations.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_s23_grouped_pathway_literature_annotations.tsv) | 10 grouped Reactome permutation and SHAP annotations |
| [`paper/synthetic_guided_spaceflight/source_data/table_s24_importance_literature_sources.tsv`](../paper/synthetic_guided_spaceflight/source_data/table_s24_importance_literature_sources.tsv) | Source catalog for matched genes and grouped pathways |

These tables record `evidence_scope`, `evidence_relationship`, `source_ids`,
`literature_summary`, `interpretation`, `annotation_origin`, and the search date.
The source tables also state whether a paper is independent, reuses public OSDR
cohorts, or supplies mechanistic context only. The exact point-model name was
not stored per generative annotation row, so the repository does not assign one
retroactively.

The structured annotation definitions and table checks are implemented in:

- [`annotate_promoted_gene_literature.py`](../src/nasa_mouse_diffusion/paper_parity/annotate_promoted_gene_literature.py)
- [`annotate_importance_literature.py`](../src/nasa_mouse_diffusion/paper_parity/annotate_importance_literature.py)

## Rebuild and checks

```bash
python -m nasa_mouse_expimap.plot_hvg_literature_review_heatmaps
python -m nasa_mouse_diffusion.paper_parity.annotate_promoted_gene_literature
python -m nasa_mouse_diffusion.paper_parity.annotate_importance_literature
python -m pytest -q tests/test_handoff_integrity.py
```

The integrity tests require every annotation to have explanatory text and every
source identifier to resolve to a tracked source record.
