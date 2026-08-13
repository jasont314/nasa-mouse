# Literature review inputs

This directory contains the compact, hand-reviewed pathway annotations used by
the expiMap paper builder.

- `manual/`: review decisions entered for each tissue.
- `sources/`: citation-key lists supporting those decisions.
- `final/`: review labels merged with the observed expiMap pathway effects.

Run `python -m nasa_mouse_expimap.plot_hvg_literature_review_heatmaps` after
refreshing the selected expiMap analyses. It updates `final/` and writes review
heatmaps and diagnostics to
`outputs/expimap/analyses/literature_review/`. Those generated diagnostics are
not part of the presentation source tree.
