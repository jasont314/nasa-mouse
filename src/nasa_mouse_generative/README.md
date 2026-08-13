# Shared generative framework

This package provides the common data, configuration, preprocessing,
conditioning, harmonization, training, and evaluation contract used by the
DDIM and WGAN workflows. It does not define the final neural architectures;
those live in `nasa_mouse_diffusion` and `nasa_mouse_wgan`.

## Command interface

Run `python -m nasa_mouse_generative <command> --help`. The main commands are:

| Command | Module | Purpose |
|---|---|---|
| `osdr-inventory` | `osdr_inventory.py` | Query and canonicalize eligible OSDR samples |
| `osdr-expression` | `osdr_expression.py` | Build the ARCHS4-aligned OSDR count matrix |
| `archs4-catalog` | `archs4_catalog.py` | Audit the full ARCHS4 mouse file and select references |
| `split-plan` | `split_plan.py` | Create grouped train, validation, and test assignments |
| `experiment-plan` | `experiment_plan.py` | Expand benchmark choices into concrete runs |
| `matrix-run` | `matrix_runner.py` | Execute or resume experiment-plan rows |
| `train`, `evaluate`, `generate` | `runner.py`, `evaluate.py`, `generate.py` | Run one resolved model configuration |
| `scoreboard` | `scoreboard.py` | Summarize completed screens using hard metric gates |
| `prepare-references` | `public_references.py` | Resume and checksum-verify the public ARCHS4 and TMS inputs |
| `prepare-upstreams` | `upstreams.py` | Restore pinned optional model-source checkouts |
| `archs4-figures`, `condition-figures` | figure modules | Build evaluation figures from completed runs |
| `gene-lengths` | `gene_lengths.py` | Build the mouse Ensembl gene-length table |

The exact final commands are in
[`outputs/COMMANDS.md`](../../outputs/COMMANDS.md). Selected configurations and
runs are listed in [`outputs/README.md`](../../outputs/README.md).

## Module map

| Area | Modules |
|---|---|
| CLI and data inventory | `__main__.py`, `public_references.py`, `osdr_inventory.py`, `osdr_expression.py`, `archs4_catalog.py`, `gene_lengths.py`, `tissues.py` |
| Configuration and planning | `config.py`, `profiles.py`, `models.py`, `paper_contracts.py`, `experiment_plan.py`, `split_plan.py`, `matrix_runner.py`, `scoreboard.py`, `upstreams.py` |
| Preprocessing and harmonization | `preprocessing.py`, `conditioning.py`, `training_data.py`, `harmonization.py`, `harmonizers.py`, `mober_harmonizer.py`, `mbatch_harmonizer.py` |
| Model execution | `runner.py`, `evaluate.py`, `generate.py`, `adapters/base.py`, `adapters/diffusion.py`, `adapters/wgan.py` |
| Metrics and figures | `metrics.py`, `paper_metrics.py`, `effect_validation.py`, `archs4_figures.py`, `condition_figures.py` |

## Data and leakage rules

- OSDR discovery uses the NASA Biological Data API.
- ARCHS4 and OSDR overlap exclusions are applied before reference selection.
- GEO series define ARCHS4 groups; accessions define OSDR groups.
- Preprocessing and harmonization fit only on training profiles.
- Locked-test assignments cannot be used for model selection.
- Generated profiles are never counted as biological replicates in association
  tests.

## Configuration status

The selected paper runs use the OSDR-disjoint ARCHS4 DDIM backbone, the
factorized OSDR DDIM adapter, and the matched study-conditioned WGAN comparator.
The remaining experiment-matrix combinations and harmonization profiles are
development or sensitivity records. See
[`configs/generative/README.md`](../../configs/generative/README.md) for the
configuration families.
