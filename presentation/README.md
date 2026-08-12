# Final internship presentation

The final presentation covers the full project, beginning with GLARE and
expiMap and then presenting the generative-model benchmark and downstream
analysis.

## Deliverables

- [`SLSTP_2026_Generative_Transcriptomics.pptx`](SLSTP_2026_Generative_Transcriptomics.pptx):
  editable 29-slide deck with embedded notes.
- [`SLSTP_2026_Generative_Transcriptomics.pdf`](SLSTP_2026_Generative_Transcriptomics.pdf):
  static 29-page export.
- [`generative_slstp_2026/speaker_notes.md`](generative_slstp_2026/speaker_notes.md):
  plain-text copy of the notes with a 13:35 pacing plan.
- [`SLSTP_template_2026.pptx`](SLSTP_template_2026.pptx): NASA presentation
  template used by the builder.

The images and compact source tables required by the builder are tracked under
`generative_slstp_2026/`. Their provenance is recorded in
`generative_slstp_2026/assets/SOURCES.md`.

## Rebuild

```bash
python -m nasa_mouse_diffusion.paper_parity.build_slstp_presentation
libreoffice --headless --convert-to pdf --outdir presentation \
  presentation/SLSTP_2026_Generative_Transcriptomics.pptx
```

The builder does not need the original midpoint PPTX. It uses tracked extracts
from that deck so a fresh clone can reproduce the current slides.

## Verification

After a rebuild, confirm that the PPTX and PDF both contain 29 pages, that every
slide has notes where expected, and that a rendered contact sheet has no clipped
or overlapping elements. The final files in this directory passed that visual
check before handoff.
