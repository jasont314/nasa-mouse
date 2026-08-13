# Internship report builder

This package contains one entry point, `build_report.py`. It rebuilds the
combined internship report figures, source manifest, HTML, and PDF from tracked
source tables and frozen figure inputs.

```bash
python -m nasa_mouse_internship_report.build_report
```

The builder does not train GLARE, expiMap, DDIM, or WGAN models. Its inputs are
under `paper/slstp_internship_report/source_data/`, and its outputs remain in
that paper package. See
[`paper/slstp_internship_report/README.md`](../../paper/slstp_internship_report/README.md)
for the document scope and page audit.
