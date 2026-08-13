# NASA mouse spaceflight transcriptomics

This repository contains the code, frozen analysis tables, manuscripts, and
presentation from a 2026 NASA Space Life Sciences Training Program internship.
The project tested three approaches to mouse bulk RNA-seq from NASA OSDR:
GLARE representation learning, Reactome-constrained expiMap models, and
conditional WGAN-GP and diffusion models.

For the completed-project handoff, start with
[`MENTOR_HANDOFF.md`](MENTOR_HANDOFF.md). It explains what is final, what is
exploratory, and which large files are stored only on the original workstation.

## Start here

| Deliverable | Purpose |
|---|---|
| [Internship report](paper/slstp_internship_report/manuscript.pdf) | Shortest complete account of the full project |
| [Poster](presentation/poster/asgsr_expimap_poster.pdf) | Print-ready expiMap research poster |
| [Midpoint presentation](presentation/midpoint/SLSTP_2026_Midpoint_Presentation.pdf) | Project status and methods at the internship midpoint |
| [Final presentation](presentation/final/SLSTP_2026_Generative_Transcriptomics.pdf) | 29-slide overview for a biology audience |
| [Editable final presentation](presentation/final/SLSTP_2026_Generative_Transcriptomics.pptx) | Slides with embedded speaker notes |
| [expiMap manuscript](paper/asgsr_expimap_hvg/manuscript.pdf) | Pathway-constrained cross-mission analysis |
| [Generative manuscript](paper/synthetic_guided_spaceflight/manuscript.pdf) | Generator validation and synthetic-informed analysis |
| [Output selections](outputs/README.md) | Final runs and analyses used by the papers |
| [Command ledger](outputs/COMMANDS.md) | Commands used to produce the retained output families |

## Main findings

- GLARE was useful as an initial representation-learning test, but study identity
  remained a major source of structure after batch correction.
- expiMap produced tissue-specific pathway hypotheses in thymus, skin, liver,
  spleen, and kidney. The paper reports the supporting study-level and member-gene
  checks rather than treating every pathway score as a biological finding.
- The selected conditional DDIM matched the real expression distribution more
  closely than the WGAN comparator on the reported validation metrics.
- Synthetic training changed FLT/GC prediction and feature ranking differently
  by tissue. Thymus had the clearest agreement between pathway analysis,
  synthetic-supported genes, and prior biology.

These are transcriptomic hypotheses. Bulk tissue composition, study design, and
limited independent replication constrain the biological interpretation.

## Repository map

| Path | Contents |
|---|---|
| [`src/`](src/README.md) | Project-owned Python packages and entry points |
| [`configs/generative/`](configs/generative/README.md) | Shared benchmark, DDIM, and WGAN configurations |
| [`data/`](data/README.md) | Tracked manifests, Reactome files, gene maps, and small reference tables |
| [`assets/`](ARTIFACTS.md) | Local large references and pinned upstream method sources |
| [`outputs/`](outputs/README.md) | Curated results plus ignored local model and matrix artifacts |
| [`docs/`](docs/README.md) | OSDR access, method provenance, and current technical references |
| [`paper/`](paper/README.md) | Internship report and two project-specific manuscripts |
| [`presentation/`](presentation/README.md) | Final deck, PDF, speaker notes, and rebuild inputs |
| [`tests/`](tests) | Unit and integration tests for data, configuration, and analysis utilities |

## Fresh-clone boundary

A clone contains the final papers, presentation, publication source tables,
small data manifests, and project code. It does not contain the 39 GB ARCHS4
file, the 2.5 GB Tabula Muris Senis reference, downloaded OSDR count tables, or
the selected model checkpoints. See [`ARTIFACTS.md`](ARTIFACTS.md) for sources,
checksums, and the task-to-artifact matrix.

## Setup

Run commands from the repository root:

```bash
conda env create -f environment.yml
conda activate nasa-mouse
python -m pip install -e .
```

For an existing environment:

```bash
conda env update -f environment.yml --prune
python -m pip install -e .
```

[`environment-lock.yml`](environment-lock.yml) records the complete Linux
environment used for the final handoff. It is a snapshot, not a portable
cross-platform specification. Use `environment.yml` for normal installation and
see [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for hardware and verification
details.

## Common tasks

Discover eligible mouse bulk RNA-seq samples through the NASA OSDR Biological
Data API:

```bash
python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics
```

Download missing OSDR count tables and prepare the shared expression cache:

```bash
python -m nasa_mouse_generative osdr-expression --download-missing
```

Run the test suite:

```bash
python -m pytest
```

The project-owned tests run from a fresh clone. Restore the pinned upstream
method checkouts before running the source-parity checks:

```bash
python -m nasa_mouse_generative prepare-upstreams
python -m pytest
```

Rebuild the internship report and presentation:

```bash
python -m nasa_mouse_internship_report.build_report
python -m nasa_mouse_diffusion.paper_parity.build_slstp_presentation
```

The complete training and analysis sequence is in
[`outputs/COMMANDS.md`](outputs/COMMANDS.md). Model-specific background is
indexed in [`docs/README.md`](docs/README.md).

## Data and release policy

OSDR expression is discovered through the NASA API. Large downloaded matrices,
checkpoints, and generated samples remain outside ordinary Git history. Curated
tables and figures needed to inspect the final claims are tracked.

The final handoff is tagged `internship-final-2026-08-13-r2`. Project-authored
content is released under the [MIT License](LICENSE). Vendored code and public
data retain their original terms; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
