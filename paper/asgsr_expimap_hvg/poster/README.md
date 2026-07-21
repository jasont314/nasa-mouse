# ASGSR expiMap poster

`asgsr_expimap_poster.pptx` is an editable, single-slide scientific poster sized at 48 x 36 inches in landscape orientation. The PDF is the print-ready export, and `asgsr_expimap_poster_preview.png` is a 4,800 x 3,600 pixel visual-check copy.

The title, workflow, tables, evidence summaries, interpretation text, and footer remain native PowerPoint objects. The two data-dense figures are rendered at 700 dpi from the paper's vector PDFs and retain more than 300 effective pixels per inch at their placed poster dimensions.

Rebuild from the repository root:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m expiMap_scarches.nasa_mouse_expimap.build_asgsr_poster
```

The builder requires `python-pptx` and uses `pdftocairo` for high-resolution figure rendering. When LibreOffice is available, it also regenerates the PDF and PNG preview.

Interpretation safeguards are embedded in the poster: lower scores are decoder-oriented latent changes rather than proof of biochemical inhibition, bulk RNA-seq cannot separate cell abundance from cell state, and the displayed tissue-state interpretations require independent cell-resolved and functional validation.
