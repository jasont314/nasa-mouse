# Editorial review record

Project-authored reader-facing prose was reviewed with
[`blader/humanizer`](https://github.com/blader/humanizer), version 2.9.1. The
latest pass was completed on 2026-08-13 in file mode under its no-fabrication
rule.

## Scope

The review covered:

- the [combined internship report](../paper/slstp_internship_report/manuscript.md);
- the [expiMap manuscript](../paper/asgsr_expimap_hvg/manuscript.md) and
  [synthetic-guided manuscript](../paper/synthetic_guided_spaceflight/manuscript.md),
  including their supplementary methods;
- the root README, mentor handoff, artifact inventory, reproducibility record,
  and maintained technical guides under `docs/`;
- paper, source-package, data, configuration, output, test, and presentation
  READMEs; and
- the final presentation's speaker notes and project-authored slide text during
  the earlier presentation audit.

Generated HTML and PDF files were rebuilt from the reviewed sources. Frozen TSV
and JSON results, configuration values, code, command examples, citations,
vendored documentation, and automatically generated output reports were not
rewritten as prose.

## Editing rules

The pass removed inflated claims, formulaic transitions, repetitive summaries,
vague references, unnecessary section previews, and wording that sounded like
an agent work log. It favored direct subjects and verbs and split long sentences
when that improved readability. Em dashes were not introduced.

The following material was held fixed:

- sample counts, metrics, effect directions, FDR values, and model settings;
- gene, pathway, tissue, accession, and software names;
- literature labels, source relationships, and citation links;
- distinctions among measured results, model outputs, and biological
  hypotheses; and
- limitations and uncertainty statements.

Humanizer was used only as an editorial checklist. It did not select results,
annotate literature, calculate statistics, or supply scientific evidence.

## Final checks

The reviewed Markdown is checked for local-link integrity by
`tests/test_handoff_integrity.py`. The document builders then render the final
PDFs from those Markdown sources. Annotation prompts and evidence provenance are
recorded separately in
[`annotation_prompts.md`](annotation_prompts.md) and
[`annotation_provenance.md`](annotation_provenance.md).
