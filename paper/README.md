# Paper and report packages

The three directories serve different readers. The internship report is the
recommended first document; the other two packages contain the detailed
project-specific analyses.

| Package | Status | Read this for |
|---|---|---|
| [`slstp_internship_report/`](slstp_internship_report/) | Final internship handoff | Combined GLARE, expiMap, and generative-model narrative |
| [`asgsr_expimap_hvg/`](asgsr_expimap_hvg/) | Research manuscript and poster package | Cross-mission Reactome pathway analysis |
| [`synthetic_guided_spaceflight/`](synthetic_guided_spaceflight/) | Research manuscript draft | Generator validation and synthetic-informed FLT/GC analysis |

Each package contains Markdown source and a rendered PDF. The project-specific
packages also contain publication figures and source tables. Read each package's
README before rebuilding it.

## Rebuild

Activate the `nasa-mouse` environment and install the repository in editable
mode first. Then run:

```bash
python -m nasa_mouse_internship_report.build_report
python -m nasa_mouse_expimap.render_asgsr_documents
python -m nasa_mouse_diffusion.paper_parity.build_synthetic_guided_paper --render-only
```

The project packages include more specialized figure and audit commands in
their own README files. Rendering does not retrain models. Refreshing the
generative paper's source tables and figures requires its ignored final analysis
outputs.

## Submission status

The files are suitable for mentor review and project preservation. Before any
external submission, review the author list, acknowledgments, repository URL,
data-access wording, and target-journal format.
