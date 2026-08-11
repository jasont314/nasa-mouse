# NASA mouse internship report

This directory contains the combined GLARE, expiMap, and conditional-generative-model internship report.

The expiMap section includes a constrained latent-space architecture, accession-level pathway heatmaps for thymus, skin, spleen, and kidney, and a schematic of the evidence and literature-annotation workflow.

The literature rationale and proposed validation tests for distributed liver, spleen, skin, and exploratory secondary responses are in `docs/distributed_response_hypotheses.md`.

Build the figures, source manifest, HTML, and PDF from the repository root:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m nasa_mouse_internship_report.build_report
```

The main-text page limit excludes the references and appendices. The build audit records the page on which references begin and is stored in `PAGE_AUDIT.md` after final rendering.
