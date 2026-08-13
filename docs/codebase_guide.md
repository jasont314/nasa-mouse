# Codebase guide

This guide explains how the retained repository fits together. Read it after
the root README and `MENTOR_HANDOFF.md`. Scientific results belong in the paper
packages; this document is for navigating code, configurations, data, outputs,
and rebuild tools.

## End-to-end flow

```text
NASA OSDR API -----------------------> API metadata and count-table cache
ARCHS4 / TMS / Reactome ------------> local references and tracked manifests
                                           |
                 +-------------------------+-------------------------+
                 |                         |                         |
               GLARE                    expiMap              DDIM and WGAN
                 |                         |                         |
                 +-------------------------+-------------------------+
                                           |
                              selected outputs and source tables
                                           |
                         papers, poster, and final presentation
```

The NASA API remains the source of OSDR sample discovery. Large downloaded
matrices and checkpoints stay outside Git. Paper-facing source tables, figures,
manuscripts, and visual deliverables are tracked.

## Repository areas

| Path | Role | Main documentation |
|---|---|---|
| [`src/`](../src/README.md) | Project-owned Python and R code | Package READMEs under each `nasa_mouse_*` directory |
| [`configs/generative/`](../configs/generative/README.md) | Shared benchmark, DDIM, and WGAN run configurations | Configuration README and selected-run list |
| [`data/`](../data/README.md) | Tracked manifests, pathway files, gene maps, and local API cache | Data README |
| [`assets/`](../assets/README.md) | Pinned upstream methods and large local references | Assets README and [`ARTIFACTS.md`](../ARTIFACTS.md) |
| [`outputs/`](../outputs/README.md) | Selected runs, analyses, compact summaries, and command history | Output README and [`COMMANDS.md`](../outputs/COMMANDS.md) |
| [`paper/`](../paper/README.md) | Internship report and project manuscripts | README in each paper package |
| [`presentation/`](../presentation/README.md) | Poster, midpoint deck, final deck, and rebuild sources | Presentation README |
| [`tests/`](../tests/README.md) | Unit, integration, parity, and handoff checks | Test-suite README |

## Source packages

### `nasa_mouse_glare`

This package owns NASA OSDR API discovery, canonical tissue labels, common
matrix-loading helpers, TMS preparation, and the GLARE/MOBER workflow. Several
other packages import its API and IO utilities. Its README separates the active
multi-tissue API workflow from the original OSD-379 and aggregate-liver work.

Read [`src/nasa_mouse_glare/README.md`](../src/nasa_mouse_glare/README.md).

### `nasa_mouse_expimap`

This package prepares tissue-specific ARCHS4 references and OSDR queries,
trains or maps expiMap models, calculates pathway effects, runs robustness
checks, and builds the expiMap paper and poster. It uses the GLARE package for
OSDR access and a few shared statistical and import helpers.

Read [`src/nasa_mouse_expimap/README.md`](../src/nasa_mouse_expimap/README.md).

### `nasa_mouse_generative`

This is the common generative framework. It defines tissue labels,
conditioning, preprocessing, harmonization, grouped splits, experiment plans,
model adapters, metrics, and run orchestration. It dispatches model work to the
diffusion and WGAN packages.

Read [`src/nasa_mouse_generative/README.md`](../src/nasa_mouse_generative/README.md).

### `nasa_mouse_diffusion`

The package root contains the earlier standalone conditional diffusion
pipeline. The final paper-matched ModelDDIM implementation and OSDR adapter are
under `paper_parity/`. Downstream classifier, feature-importance, grouped
Reactome, annotation, manuscript, and presentation code also lives there.

Read [`src/nasa_mouse_diffusion/README.md`](../src/nasa_mouse_diffusion/README.md)
and
[`src/nasa_mouse_diffusion/paper_parity/README.md`](../src/nasa_mouse_diffusion/paper_parity/README.md).

### `nasa_mouse_wgan`

This package contains the conditional WGAN-GP implementation. The selected
paper comparator uses `matched_study.py`; the other runners preserve earlier
per-tissue and pan-tissue experiments.

Read [`src/nasa_mouse_wgan/README.md`](../src/nasa_mouse_wgan/README.md).

### `nasa_mouse_internship_report`

This small package builds the combined internship report from tracked source
tables and figures. It does not train models.

Read
[`src/nasa_mouse_internship_report/README.md`](../src/nasa_mouse_internship_report/README.md).

## Dependency shape

The packages are organized by workflow, not as isolated libraries:

- GLARE provides shared OSDR and matrix utilities.
- The shared generative framework imports DDIM and WGAN adapters.
- The diffusion and WGAN analysis modules reuse shared metrics and effect
  validation.
- expiMap reuses OSDR loading, FDR, and accession-effect helpers but maintains
  its own model pipeline.
- Paper and presentation builders read frozen source tables instead of
  retraining models.

Because these packages share utilities, install the repository as one editable
project rather than installing package directories independently.

## Active and retained code

The repository keeps three kinds of project-owned code:

1. Final workflow code produces the analyses cited by the current papers.
2. Supporting code prepares data, evaluates sensitivity, creates figures, or
   checks publication artifacts.
3. Development code preserves earlier GLARE, direct expiMap, standalone DDIM,
   and standalone WGAN experiments recorded in the command ledger.

Development code is not evidence for the final conclusions unless a paper or
`outputs/README.md` points to its result. Package READMEs label these groups.

## Configuration and output lifecycle

Generative runs begin with YAML under `configs/generative/`. The shared runner
resolves each configuration, fits preprocessing only on its training split,
writes checkpoints and metrics under `outputs/generative/benchmark/`, and adds
compact results to the retained source tables. The selected configurations and
output directories are named in `outputs/README.md`; other YAML files preserve
development screens and sensitivity analyses.

expiMap uses command-line arguments rather than one shared YAML contract. Its
selected model directories and exact commands are recorded in the expiMap
supplementary methods and `outputs/COMMANDS.md`.

## How to trace a result

Use this order when checking a claim:

1. Find the statement and table or figure number in the relevant manuscript.
2. Open the adjacent `source_data/` table.
3. Use `outputs/README.md` to identify the selected analysis or model.
4. Use `outputs/COMMANDS.md` to find the generating command and configuration.
5. Open the package README to identify the entry point and supporting modules.
6. Consult `ARTIFACTS.md` if the command requires an ignored matrix or model.

Literature labels and their cited sources follow the separate map in
[`annotation_provenance.md`](annotation_provenance.md).

## What a fresh clone can do

A fresh clone can run tests, inspect all final source tables, and rebuild the
documents and presentation. ARCHS4 and TMS can be downloaded with
`python -m nasa_mouse_generative prepare-references`, and OSDR inputs can be
recreated through the NASA API. Selected checkpoints are needed only to rerun
inference without retraining. The exact environment and verification commands
are in
[`REPRODUCIBILITY.md`](../REPRODUCIBILITY.md).
