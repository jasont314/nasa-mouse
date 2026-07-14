# Generative Benchmark Outputs

This directory contains current OSDR/ARCHS4 data audits, resolved experiment plans,
and future model runs for `src/nasa_mouse_generative/`.

Large compressed per-profile ARCHS4 manifests are generated locally and ignored by
Git. Small summaries and manifests retain the selection rules and provenance.

Current audited inputs:

- OSDR API: 1,631 profile rows, 1,610 biological profiles after technical-replicate
  summation, 75 accessions, 24 canonical classes, and 48,694 ARCHS4-shared genes;
- ARCHS4 healthy-preferred: 62,299 profiles from 5,307 GEO series across 23
  matchable classes;
- ARCHS4 control-only sensitivity: 23,614 profiles;
- ARCHS4 broad-diversity sensitivity: 134,250 profiles.

`splits/` contains accession-grouped validation plans. No model should use the
locked-test assignments for preprocessing, checkpoint selection, or hyperparameter
tuning.
