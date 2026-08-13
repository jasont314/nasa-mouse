# Data and model artifacts

Git contains the code, configs, compact result tables, figures, manuscripts, and
presentation. Large matrices and model checkpoints are ignored to keep clones
manageable. Paths below are relative to the repository root.

## Public reference inputs

| Local path | Size on handoff machine | Source | SHA-256 |
|---|---:|---|---|
| `assets/archs4/mouse_gene_v2.5.h5` | 38,960,132,574 bytes | [ARCHS4 download page](https://maayanlab.cloud/archs4/download.html), [direct mouse file](https://s3.dev.maayanlab.cloud/archs4/files/mouse_gene_v2.5.h5) | `74b509f82623bced395119244becf30df601a24fcaaf905691e2716bf83118b8` |
| `assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad` | 2,548,190,251 bytes | [CELLxGENE-hosted Tabula Muris Senis FACS file](https://datasets.cellxgene.cziscience.com/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad) | `1d7fd90acb33269c3337dc5031b4a89d9aa4f72806a45b9c12e768fedc8acf8f` |

The same checksums are stored in `assets/EXTERNAL_ARTIFACTS.sha256` so local
copies can be checked with:

```bash
python -m nasa_mouse_generative prepare-references --check
```

Both files are public and can be restored directly from their versioned source
URLs. The downloader resumes partial transfers and checks the expected byte
count and SHA-256 digest before installing either file:

```bash
python -m nasa_mouse_generative prepare-references --list
python -m nasa_mouse_generative prepare-references
```

The two downloads require about 41.5 GB of disk space before filesystem and
temporary-file overhead. They do not need to be copied from the handoff machine.
The equivalent manifest-only check remains available as
`sha256sum --check assets/EXTERNAL_ARTIFACTS.sha256`.

The supplied NASA poster template is tracked at
`assets/poster_template/00 Poster Session Student Template_Approved by Legal.pptx`,
so the ASGSR poster can be rebuilt from a fresh clone. The generated poster PPTX
and PDF are tracked under `presentation/poster/`.

## NASA OSDR cache

Small OSDR metadata and inventory files are tracked under `data/osdr_api/` and
`outputs/osdr/inventory/`. Per-accession count CSVs under
`data/osdr_api/counts/` are ignored and can be downloaded again:

```bash
python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics --download-counts
python -m nasa_mouse_generative osdr-expression --download-missing
```

OSDR is a live service. A later download may include revised metadata or newly
released studies, so retained manifests should be used to audit differences.

## Model checkpoints

The final local allowlist is recorded in
[`outputs/MODEL_ARTIFACTS.sha256`](outputs/MODEL_ARTIFACTS.sha256). It covers:

- the selected ARCHS4 DDIM backbone;
- the selected OSDR-conditioned DDIM adapter;
- the WGAN comparison checkpoints;
- the expiMap reference and query models used for thymus, skin, liver, and
  soleus, including their gene-order metadata.

Checkpoints are ignored by Git and currently have no public download URL. A
fresh clone can inspect the frozen analyses and rebuild the documents, but it
cannot rerun model inference until these files are copied from the handoff
machine or placed in managed storage.

Refreshing every table and figure in the detailed generative manuscript also
requires the compact files listed in
`paper/synthetic_guided_spaceflight/source_data/frozen_input_manifest.tsv`.
The final rendered paper, figures, and derived source tables are tracked and do
not require those files for review or render-only builds.

## Upstream method code

GLARE and MOBER source snapshots required by the active workflows are tracked
under `assets/model_sources/`. Larger upstream development checkouts for the
WGAN, RNA diffusion, MBatch, and TRRAC comparisons are ignored. Their repository
URLs and exact commits are recorded in [docs/method_sources.md](docs/method_sources.md).
They can be restored with:

```bash
python -m nasa_mouse_generative prepare-upstreams
```

## What each task needs

| Task | ARCHS4 | TMS | OSDR count cache | Checkpoints | Final analysis outputs |
|---|---:|---:|---:|---:|---:|
| Read final papers and presentation | No | No | No | No | No |
| Redraw frozen-source figures; render papers, poster, report, and final presentation | No | No | No | No | No |
| Refresh OSDR inventory | No | No | Downloaded as needed | No | No |
| Refresh detailed generative-paper tables and figures | No | No | No | No | Yes |
| Run selected model inference | No | No | Depends on analysis | Yes | No |
| Retrain GLARE | No | Public download | API download | No | No |
| Retrain expiMap reference-query models | Public download | No | API download | No | No |
| Retrain ARCHS4 DDIM and OSDR adapter | Public download | No | API download | No | No |

## Transfer recommendation

ARCHS4, TMS, OSDR, and the optional upstream method repositories can all be
restored from public sources. They may be mirrored for convenience, but they are
not unique handoff artifacts.

Before deleting the original workstation, preserve the files named in
`outputs/MODEL_ARTIFACTS.sha256` if avoiding full model retraining matters.
Preserving the compact files in the generative paper's
`frozen_input_manifest.tsv` is also recommended for an exact source-table
refresh. Record the storage URL, access policy, and retrieval date here. The
ignored root `.env` is machine-specific and may contain credentials; do not
include it in a handoff archive.

The boundary between figures rebuilt from tracked tables and graphics that
require model outputs is documented in
[`docs/figure_reproduction.md`](docs/figure_reproduction.md).
