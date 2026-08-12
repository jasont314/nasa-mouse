# GLARE source snapshot

This directory contains the subset of the upstream GLARE manuscript code used
by the NASA mouse workflow. Project-owned data preparation, training wrappers,
and analyses are under `src/nasa_mouse_glare/`.

## Local changes

The retained source includes compatibility fixes needed by the current runtime:

- direct script execution resolves local modules correctly;
- sparse MatrixMarket input uses SciPy's current dense-conversion API;
- fine-tuning reuses the pretraining architecture; and
- final representation extraction applies the fine-tuning adapter.

Upstream demonstrations, example datasets, historical weights, and manuscript
media are omitted. The upstream license is preserved in this directory and in
`Manuscript_Code/glare/LICENSE`.

See `docs/method_sources.md` for method provenance and
`src/nasa_mouse_glare/README.md` for the active NASA OSDR workflow.
