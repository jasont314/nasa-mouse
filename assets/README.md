# Assets

This directory separates upstream method code from large local reference files.

| Path | Git status | Purpose |
|---|---|---|
| `model_sources/glare/` | Tracked | GLARE runtime snapshot with local compatibility fixes |
| `model_sources/MOBER/` | Tracked | MOBER source used by the batch-correction workflow |
| `model_sources/MBatch/` | Ignored | Pinned development checkout for R harmonization adapters |
| `model_sources/adversarial-gene-expression/` | Ignored | WGAN paper implementation used for architecture review |
| `model_sources/rna-diffusion/` | Ignored | DDIM paper implementation used for parity checks |
| `model_sources/trrac/` | Ignored | Spaceflight batch-correction reference implementation |
| `archs4/` | Ignored | Full ARCHS4 mouse expression H5 |
| `tms/` | Ignored | Tabula Muris Senis reference H5AD |
| `poster_template/` | Ignored | Supplied NASA poster source templates |

See [`../ARTIFACTS.md`](../ARTIFACTS.md) for download locations, checksums, and
the tasks that require each large file. Upstream repository commits are listed
in [`../docs/method_sources.md`](../docs/method_sources.md).
