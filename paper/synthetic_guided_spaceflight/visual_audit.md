# Visual rendering audit

Audit date: August 13, 2026 UTC.

The 18-page manuscript PDF was rebuilt from the reviewed Markdown and tracked
figures. Every page was rendered to PNG and inspected in a numbered contact
sheet. Figure pages and the title page were also checked at their rendered page
size.

The audit found:

- no clipped or overlapping text;
- no missing or blank figures;
- no split tables that obscured a row or heading;
- readable axes, labels, legends, captions, and references;
- no malformed replacement characters in extracted PDF text; and
- no unresolved bracketed drafting placeholders.

The paper contains five main figures and two supplementary figures. PNG and
vector PDF copies of each figure are tracked under `figures/`, and their hashes
are stored in `source_data/figure_build_manifest.tsv`. Figures 1, 3, 4, and 5
and Figure S2 can be redrawn from tracked source tables. Figure 2 and Figure S1
are preserved model-output graphics; their full regeneration path is documented
in [`docs/figure_reproduction.md`](../../docs/figure_reproduction.md).

Rebuild command:

```bash
python \
  -m nasa_mouse_diffusion.paper_parity.build_synthetic_guided_paper \
  --figures-from-frozen-source
```
