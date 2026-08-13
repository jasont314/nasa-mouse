# Posters and presentations

This directory contains the three visual deliverables from the internship.

| Directory | Deliverable |
|---|---|
| [`poster/`](poster/) | Editable ASGSR expiMap poster, print PDF, preview, and embedded source assets |
| [`midpoint/`](midpoint/) | Midpoint presentation in editable PowerPoint and PDF formats |
| [`final/`](final/) | Final 29-slide presentation, PDF, speaker notes, template, and reproducible source files |

## Rebuild

Build the poster:

```bash
python -m nasa_mouse_expimap.build_asgsr_poster
```

Build and export the final presentation:

```bash
python -m nasa_mouse_diffusion.paper_parity.build_slstp_presentation
libreoffice --headless --convert-to pdf --outdir presentation/final \
  presentation/final/SLSTP_2026_Generative_Transcriptomics.pptx
```

The midpoint presentation is the original delivered deck. The final builder
uses selected tracked extracts from it under `final/source/assets/`; it does not
modify the midpoint deck. The midpoint PPTX is its editable source and can be
exported again with LibreOffice:

```bash
libreoffice --headless --convert-to pdf --outdir presentation/midpoint \
  presentation/midpoint/SLSTP_2026_Midpoint_Presentation.pptx
```

The poster and final presentation were visually checked after rendering. The
PowerPoint files are the editable sources; the PDFs are the review and print
copies.

See [`docs/figure_reproduction.md`](../docs/figure_reproduction.md) for the
figure-level source boundary and full-analysis commands.
