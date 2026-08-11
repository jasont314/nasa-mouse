# Method Repositories

Upstream method code and reference material used by the NASA mouse workflows.
Project-owned packages remain under `src/`; upstream source is kept under
`assets/model_sources/` or installed as a pinned dependency.

## GLARE

- `assets/model_sources/glare`
  - purpose: original GLARE manuscript code and model assets.
  - local changes: runtime compatibility fixes used by
    `src/nasa_mouse_glare/`.

## expiMap

Paper: "Biologically informed deep learning to query gene programs in
single-cell atlases"

Article: https://www.nature.com/articles/s41556-022-01072-x

- `scarches==0.6.1`
  - upstream: https://github.com/theislab/scarches
  - purpose: installed expiMap implementation used by
    `src/nasa_mouse_expimap/`.
- `data/reference/expimap/paper_metadata`
  - upstream: https://github.com/theislab/expiMap_reproducibility
  - source commit: `295ac3c0fff29b8c9e33bc412c8e8282201b0be2`
  - purpose: Reactome and PanglaoDB reference files used by validation code.

## MOBER

Method: "Multi Origin Batch Effect Remover"

- `assets/model_sources/MOBER`
  - upstream: https://github.com/Novartis/MOBER
  - branch: `main`
  - commit: `81a628322044eb53e95bb606a65b4359fa645085`
  - purpose: batch/source correction for aggregate OSDR liver RNA-seq.
  - local changes: Scanpy and MLflow are treated as optional dependencies so
    the package can run in the existing `nasa` environment with AnnData and
    local file logging.

## Generative model references

These pinned upstream checkouts are local development assets and are ignored by
Git. Project-owned adapters, training code, and evaluation workflows remain in
`src/nasa_mouse_generative/`, `src/nasa_mouse_wgan/`, and
`src/nasa_mouse_diffusion/`.

- `assets/model_sources/adversarial-gene-expression`
  - upstream: https://github.com/rvinas/adversarial-gene-expression
  - commit: `94fa44dd1bd52d924efd3af0fcd8eeb18bd141a8`
  - purpose: paper implementation used to verify the WGAN-GP architecture and
    training contract.
- `assets/model_sources/rna-diffusion`
  - upstream: https://forge.ibisc.univ-evry.fr/alacan/rna-diffusion.git
  - commit: `cde890154698fcea96c924804aaff04af3351b48`
  - purpose: paper implementation and L1000 metadata used by the DDIM
    paper-parity workflow.
- `assets/model_sources/MBatch`
  - upstream: https://github.com/MD-Anderson-Bioinformatics/MBatch
  - commit: `93cddd2ba18ed8781b9865ba0259fafa057bcc17`
  - purpose: official R implementations used by the MBatch harmonization
    adapters.
- `assets/model_sources/trrac`
  - upstream: https://github.com/nasa/trrac
  - commit: `abcf12d57d68a36a4628f83dec191e2b2a6b778e`
  - purpose: reference implementation used by the spaceflight batch-correction
    adapter.
