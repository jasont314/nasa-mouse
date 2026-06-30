# TabPFN3 Results Status

## Current Status

The TabPFN3 project is implemented, but production TabPFN3 fitting has not run
because the official `tabpfn` package requires Prior Labs license acceptance and
a `TABPFN_TOKEN` in this non-interactive environment.

Observed backend status:

- package installed: `tabpfn==8.0.8`
- selected device: `cuda`
- local GPU: A100 detected
- local token/cache: no `TABPFN_TOKEN`; a TabPFN cache directory exists after
  failed initialization, but no usable authorized model weights are available
- production output root: `outputs/tabpfn3_osdr`
- status file: `outputs/tabpfn3_osdr/summary/tabpfn3_backend_status.json`
- run manifest: `outputs/tabpfn3_osdr/summary/tabpfn3_run_manifest.tsv`

The manifest contains blocked rows for every planned tissue/split and both
feature modes, rather than substituting a different model.

## OSDR API Inventory

Primary tissue counts from the generated inventory:

| tissue | flight | ground_control | total |
|---|---:|---:|---:|
| liver | 125 | 118 | 243 |
| skeletal_muscle | 95 | 96 | 191 |
| skin | 80 | 71 | 151 |
| kidney | 68 | 67 | 135 |
| thymus | 63 | 54 | 117 |
| spleen | 55 | 54 | 109 |
| lung | 40 | 38 | 78 |
| retina | 45 | 31 | 76 |

Skeletal-muscle split counts:

| muscle_group | flight | ground_control | total |
|---|---:|---:|---:|
| soleus | 28 | 25 | 53 |
| quadriceps | 23 | 23 | 46 |
| edl | 16 | 16 | 32 |
| gastrocnemius | 13 | 17 | 30 |
| tibialis_anterior | 15 | 15 | 30 |

## Smoke Validation

A small `sklearn_logreg` smoke run completed under
`outputs/tabpfn3_osdr_smoke_ensembl`. It used liver only, 50 HVG features, random
CV, and 5 permutation-importance candidates. This confirms code mechanics only;
it is not a TabPFN3 result and should not be compared biologically with expiMap,
OntoVAE, WGAN, or diffusion outputs.

The Ensembl-only filter was validated in the smoke feature-importance table:
all reported genes started with `ENSMUSG`.

## Next Run

To produce real TabPFN3 metrics:

1. Accept the Prior Labs license.
2. Export `TABPFN_TOKEN`.
3. Rerun `PYTHONPATH=src python -m nasa_mouse_tabpfn3.run_osdr_classification`.
4. Compare accession-aware metrics before interpreting tissue biology.

No significant TabPFN3 biological signals can be claimed until that run
completes.

