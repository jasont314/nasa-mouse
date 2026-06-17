# Method Repositories

Vendored method code and reproducibility material used as references for
biological interpretation of mouse spaceflight expression data.

These directories are copied as editable source trees, not git submodules.

## expiMap

Paper: "Biologically informed deep learning to query gene programs in
single-cell atlases"

Article: https://www.nature.com/articles/s41556-022-01072-x

- `expiMap_scarches`
  - upstream: https://github.com/theislab/scarches
  - branch: `soft_new_mask`
  - commit: `80358fcd0b303e93d5e437bed04c455210f3bfed`
  - purpose: expiMap/scArches method implementation referenced by the
    reproducibility repo.
- `expiMap_reproducibility`
  - upstream: https://github.com/theislab/expiMap_reproducibility
  - branch: `main`
  - commit: `295ac3c0fff29b8c9e33bc412c8e8282201b0be2`
  - copied content: `README.md`, `LICENSE`, `scripts/`, and `metadata/`.
  - note: the full upstream repo contains about 600 MB of executed notebooks,
    so those large notebook directories were not vendored here.

## OntoVAE

Paper: "OntoVAE: ontology-guided variational autoencoders for
interpretable representation learning and perturbation prediction"

Article: https://academic.oup.com/bioinformatics/article/39/6/btad387/7199588

- `onto-vae`
  - upstream: https://github.com/hdsu-bioquant/onto-vae
  - branch: `main`
  - commit: `a0007555d3c92b45288f453c14f418b339f17c79`
  - purpose: Python package for ontology preprocessing and OntoVAE training.
- `OntoVAE_manuscript`
  - upstream: https://github.com/hdsu-bioquant/OntoVAE_manuscript
  - branch: `main`
  - commit: `2e0cb284288d46b86ea83341bfc00e879fb05fa8`
  - purpose: manuscript analysis/reproduction code.

## VEGA

Paper: "VEGA is an interpretable generative model for inferring biological
network activity in single-cell transcriptomics"

Article: https://www.nature.com/articles/s41467-021-26017-0

- `vega`
  - upstream: https://github.com/LucasESBS/vega
  - branch: `main`
  - commit: `146c6c80aa2904138d2cfac7459861d2cb812cb1`
  - purpose: Python package/API for VEGA.
  - note: the paper's code availability section cites
    `https://github.com/LucasESBS/vega/tree/vega_dev`, but that branch is no
    longer available upstream; the current package source was vendored from
    `main`.
- `vega-reproducibility`
  - upstream: https://github.com/LucasESBS/vega-reproducibility
  - branch: `main`
  - commit: `2474b51dc50654f0c6a5509cdebd179747c9e8f5`
  - purpose: manuscript reproduction code, small data assets, trained-model
    placeholders, and plotting scripts.
