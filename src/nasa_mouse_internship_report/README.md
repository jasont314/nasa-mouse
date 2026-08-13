# Internship report builder

This package contains the combined report builder and the cross-method selected
feature exporter. `build_report.py` rebuilds the report figures, source
manifest, HTML, and PDF from tracked source tables and frozen figure inputs.

```bash
python -m nasa_mouse_internship_report.build_report
```

`build_comparison_exports.py` creates the mentor-facing expiMap and generative
gene/pathway workbook and TSV bundle:

```bash
python -m nasa_mouse_internship_report.build_comparison_exports
```

The comparison exporter needs the completed local generative classifier
importance tables. The tracked comparison bundle can be inspected without
those local inputs.

The builder does not train GLARE, expiMap, DDIM, or WGAN models. Its inputs are
under `paper/slstp_internship_report/source_data/`, and its outputs remain in
that paper package. See
[`paper/slstp_internship_report/README.md`](../../paper/slstp_internship_report/README.md)
for the document scope and page audit.
