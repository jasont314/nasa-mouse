# Reproducibility record

This file records the final handoff environment and the minimum checks for a
fresh clone. It complements the data and checkpoint inventory in
`ARTIFACTS.md`.

## Environment

The supported setup starts from `environment.yml`. The exact Linux environment
present on 2026-08-12 is frozen in `environment-lock.yml`.

The handoff machine reported:

| Component | Value |
|---|---|
| Operating system | Linux 7.0.0-28-generic, x86_64, glibc 2.39 |
| Python | 3.11.15 |
| PyTorch | 2.12.1+cu130 |
| PyTorch CUDA build | 13.0 |
| GPU | NVIDIA A100-SXM4-40GB |
| DESeq2 | 1.50.2 |
| sva | 3.58.0 |

The final expiMap and generative training records also identify an NVIDIA
A100-SXM4-40GB. CPU-only inspection, document builds, and most tests are
supported, but model training requires a CUDA-capable environment and was not
benchmarked on CPU.

## Install

```bash
conda env create -f environment.yml
conda activate nasa-mouse
python -m pip install -e .
```

Use `environment-lock.yml` only when recreating the final Linux package set is
more important than portability. CUDA packages and system libraries make that
snapshot unsuitable for macOS and most CPU-only machines.

## Verification

Run these checks from the repository root:

```bash
python -m pytest
python -m compileall -q src
sha256sum --check assets/EXTERNAL_ARTIFACTS.sha256
sha256sum --check outputs/MODEL_ARTIFACTS.sha256
```

The checksum commands require the ignored local assets. A fresh clone should
skip the corresponding check until those files have been transferred.
Tests that compare directly against the pinned upstream WGAN and DDIM sources
skip until `python -m nasa_mouse_generative prepare-upstreams` restores those
checkouts.

Document-only verification does not require the large data or models:

```bash
python -m nasa_mouse_internship_report.build_report
python -m nasa_mouse_diffusion.paper_parity.build_slstp_presentation
python -m nasa_mouse_expimap.render_asgsr_documents
python -m nasa_mouse_diffusion.paper_parity.build_synthetic_guided_paper --render-only
```

These commands read tracked source tables, figures, and manuscript files. They
do not retrain a model. A full refresh of the detailed generative paper requires
the ignored final analysis outputs listed in that paper's source manifest.

## Provenance conventions

- YAML configs are the intended inputs for config-driven model runs.
- Each completed run stores a resolved config or training summary in its output
  directory.
- `outputs/COMMANDS.md` records the final command sequence.
- `outputs/README.md` identifies which runs were used by the final papers.
- Publication-facing source tables are frozen under each paper directory.
- Historical absolute paths in archived logs and provenance manifests describe
  the original workstation and are not required runtime paths.
