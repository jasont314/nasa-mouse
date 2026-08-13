# Literature annotation prompts

This file records how LLM-assisted literature annotations were produced and how
to repeat the review. It should be read with
[`annotation_provenance.md`](annotation_provenance.md), which points to the
completed annotations and source catalogs.

## Historical record

The original review was conducted interactively rather than through one saved
API request. The repository therefore cannot provide a verbatim system prompt,
model snapshot, or decoding parameters for every row. The instruction sequence
was preserved in the project conversation and had three final requirements:

1. Review the genes and pathways selected by the quantitative analyses against
   prior spaceflight and relevant mechanistic literature.
2. Keep synthetic-selection status (`promoted` or `reinforced`) independent of
   literature relationship (`aligning`, `complementary`, `ambiguous`, or
   `unmatched`).
3. Record the reason and sources for every label, cover all eligible tissues,
   and reserve the paper narrative for the strongest results.

The prompts below consolidate those instructions into rerunnable review
protocols. They reproduce the decision rules and output fields used in the
committed tables, but they are not presented as verbatim chat transcripts.
Literature searches and model outputs can change over time. The committed
tables, source records, search dates, and validation scripts are the frozen
historical result.

## Shared review rules

Both annotation protocols used these rules:

- Start only after the statistical or feature-selection candidate set is
  fixed. Do not add or remove a candidate because its biology seems familiar.
- Do not calculate expression effects, FDR values, classifier scores, or
  feature importance with the LLM.
- Prefer primary papers and authoritative database records. Use reviews only
  to locate primary evidence or provide broad context.
- Check whether a source is independent of the analyzed OSDR data, reuses a
  public OSDR cohort, or supplies mechanistic context only.
- Distinguish an exact same-tissue, same-feature, same-direction match from a
  process-level match and from evidence in another tissue, species, or model.
- Do not call a result novel merely because a targeted search found no exact
  match. Use `unmatched` and describe the search boundary.
- Return a concise rationale that a reader can audit. Do not return hidden
  reasoning or unsupported causal claims.
- Preserve the observed FLT-minus-GC direction supplied in the input.

## Prompt `expimap-pathway-review-v1`

This protocol applies to the expiMap pathway review under
[`paper/asgsr_expimap_hvg/source_data/literature_review/`](../paper/asgsr_expimap_hvg/source_data/literature_review/)
and the corrected kidney and spleen source records.

```text
You are reviewing a fixed set of tissue-specific mouse spaceflight pathway
results. The quantitative analysis is complete. Do not change pathway scores,
effects, FDR values, or the candidate set.

For each input row, use:
- tissue
- Reactome pathway identifier and name
- decoder-oriented FLT-minus-GC direction
- project-level and sensitivity evidence supplied with the row

Search primary literature for the same tissue and biological process in
spaceflight, microgravity, unloading, radiation, or closely related mechanisms.
For every row:
1. classify the literature relationship as literature aligned,
   complementary, or context sensitive;
2. state whether the reported direction agrees, conflicts, is mixed, or cannot
   be compared directly;
3. assign a confidence level based on evidence scope;
4. write a short rationale that distinguishes direct tissue evidence from
   process-level or mechanistic context;
5. list stable citation keys that resolve to title, DOI or PMID, and URL in the
   source catalog.

Do not infer that a pathway score proves pathway activity, a cell-type change,
or causality. Do not use literature familiarity as a statistical validation.
Return one structured row per supplied pathway and do not omit negative or
context-sensitive rows.
```

The manually reviewed tables store the resulting fields as
`literature_alignment`, `direction_assessment`, `confidence`, `rationale` or
`review_rationale`, and `citations`.

## Prompt `synthetic-feature-review-v2`

This protocol applies to consensus-ranked genes, matched permutation and SHAP
genes, and grouped Reactome permutation and SHAP pathways in the synthetic
analysis.

```text
You are reviewing a fixed set of mouse spaceflight transcriptomic candidates
selected after classifier analysis. Do not rerank candidates, alter selection
status, calculate statistics, or treat generated profiles as biological
replicates.

For each input row, use:
- analysis scope and tissue or anatomical muscle group
- gene symbol or Reactome pathway
- FLT-minus-GC direction measured in real OSDR samples
- selection status from the analysis, kept as a separate field
- classifier or grouped-feature evidence supplied with the row

Search primary spaceflight literature first, followed by unloading,
microgravity, same-tissue disease or injury, and direct gene or pathway
mechanism papers when needed. Classify the literature relationship using
exactly one label:
- aligning: prior evidence supports the same gene or process in a compatible
  tissue and direction;
- complementary: prior work supports a related process or mechanism, while the
  exact tissue-gene-direction result is not a replication;
- ambiguous: prior directions or contexts are mixed or materially
  incompatible;
- unmatched: the targeted search found no relevant prior match beyond general
  plausibility.

For each row, report:
- evidence_scope
- evidence_relationship
- source_ids
- literature_summary
- interpretation

State whether each source is independent, reuses public OSDR cohorts, or is
mechanistic context only. Keep promoted or reinforced status independent from
the literature label. Do not describe unmatched candidates as proven novelty.
Do not make a causal or functional claim from bulk RNA-seq alone. Return one
structured row for every supplied candidate.
```

The final four labels replaced an earlier draft category named
`contradictory`. Mixed or opposing evidence is recorded as `ambiguous`, with the
specific disagreement explained in `evidence_relationship` and
`literature_summary`.

## Repeating the review

The scripts do not call an LLM. They freeze the reviewed decisions and verify
that candidate coverage, labels, and source identifiers remain intact:

```bash
python -m nasa_mouse_diffusion.paper_parity.annotate_promoted_gene_literature --check
python -m nasa_mouse_diffusion.paper_parity.annotate_importance_literature --check
python -m pytest -q tests/test_handoff_integrity.py
```

For a new review, export the fixed candidate rows, run the appropriate prompt,
check every cited source manually, and encode the accepted rows in the matching
annotation script or expiMap manual table. Keep rejected LLM suggestions out of
the frozen record rather than silently changing the quantitative candidate set.
