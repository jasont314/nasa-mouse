# NASA mouse internship report

This directory contains the combined GLARE, expiMap, and conditional-generative-model internship report.

The expiMap section includes a constrained latent-space architecture, accession-level pathway heatmaps for thymus, skin, spleen, and kidney, and a schematic of the evidence and literature-annotation workflow.

The literature rationale and proposed validation tests for distributed liver, spleen, skin, and exploratory secondary responses are in `docs/distributed_response_hypotheses.md`.

The structured expiMap and feature-importance literature annotations, including
their rationale fields and source catalogs, are indexed in
[`docs/annotation_provenance.md`](../../docs/annotation_provenance.md).
The consolidated review instructions and final label definitions are in
[`docs/annotation_prompts.md`](../../docs/annotation_prompts.md).

The compact tables required to rebuild the report are tracked under
`source_data/`; the builder does not depend on ignored training outputs.

Build the figures, source manifest, HTML, and PDF from the repository root:

```bash
python \
  -m nasa_mouse_internship_report.build_report
```

The main-text page limit excludes the references and appendices. The build audit records the page on which references begin and is stored in `PAGE_AUDIT.md` after final rendering.
