# RNA Diffusion Paper-Parity Mouse Baseline

This project tests whether the GTEx DDIM from Lacan et al. transfers to mouse bulk
RNA-seq when the model and training procedure are held fixed and only the data
interface is adapted. It replaces the earlier reduced proxy experiment, whose model
and output directories were removed.

## Source Contract

The implementation imports `ModelDDIM` directly from the vendored official source
at commit `cde890154698fcea96c924804aaff04af3351b48`. Before import it verifies the
commit plus SHA-256 hashes for the upstream model, loss, and denoising files. Unit
tests compare the local loss and final DDIM state numerically with those upstream
functions.

The fixed model/training settings are:

- two 8,192-unit residual layers and 227,109,786 trainable parameters for 20 tissues;
- the upstream full-one-hot tissue embedding behavior with embedding dimension 2;
- 1,000 quadratic diffusion steps, beta range 0.0001 to 0.02, and DDIM `eta=0`;
- direct summed noise-prediction MSE, antithetic timesteps, AMP, and EMA 0.999;
- Adam at the retained search rate 0.0004783833151836702;
- OneCycle with its maximum at global optimizer step 1,000;
- batch size 2,048, 15,000 epochs, unweighted shuffled sampling, and no gradient
  clipping.

## Mouse Data Substitution

The source data are healthy-preferred, bulk-like mouse ARCHS4 profiles selected from
the full 997,515-profile file and its audited metadata. Tissues with fewer than 100
eligible profiles are excluded, leaving 20 classes. Seeded within-tissue shuffling
and round-robin selection provide 17,244 unique profiles, followed by the paper's
9,796 train, 2,448 validation, and 5,000 test sample counts.

The accepted manuscript explicitly reports those three split sizes. The current
upstream preprocessing script first writes a 12,244-profile train-plus-validation
file, and its final DDIM loader fits on that whole file. This run follows the reported
manuscript split and keeps validation separate. The distinction is recorded because
it prevents describing this as a bit-for-bit execution of the authors' private final
run, even though the published architecture, objective, optimizer, scheduler, and
training duration are retained.

TPM is computed using the denominator across all 52,848 ARCHS4 genes with positive
GENCODE M39 union-exon lengths. Landmark genes are selected only after TPM, and
MaxAbs scaling is fitted only on training profiles. Held-out values are not clipped;
they may exceed one, as with the paper's train-fitted scaler.

The deterministic mouse panel has 974 unique Ensembl genes. It covers 954 human
landmarks through the existing Ensembl orthology map and ten more through exact
case-insensitive mouse symbol matches. Ten additional mapped paralogs preserve the
paper's 974 input dimensions for human landmarks lacking a recoverable direct mouse
counterpart. This mapping is a declared data-interface limitation, not an exact
one-to-one orthology claim. See `data/diffusion/l974_mouse_paper_parity.tsv`.

## Commands

```bash
PYTHONPATH=src python -m nasa_mouse_rna_diffusion prepare \
  --config configs/rna_diffusion/archs4_mouse_paper_parity.yaml

PYTHONPATH=src python -m nasa_mouse_rna_diffusion train \
  --config configs/rna_diffusion/archs4_mouse_paper_parity.yaml

PYTHONPATH=src python -m nasa_mouse_rna_diffusion evaluate \
  --config configs/rna_diffusion/archs4_mouse_paper_parity.yaml
```

Training resumes from `checkpoints/latest.pt`. Evaluation loads EMA weights, runs all
1,000 DDIM steps, and writes the paper-style 1,024-profile `t=1000`, `t=200`, and
`t=0` PCA trajectory. It separately generates 9,796 profiles following the real
training-label distribution for tissue utility and quantitative metrics.

## Outputs

Run outputs are under:

`outputs/generative_benchmark/runs/lacan_diffusion/archs4_mouse_paper_parity_seed1234/`

Prepared data are under:

`outputs/generative_benchmark/data/lacan_diffusion/`

The human GTEx tissue-classifier weights are not transferable to mouse genes.
Therefore Frechet distance is reported in a PCA-50 embedding fitted only on mouse
training data, precision/recall are reported in both scaled L974 and PCA-50 spaces,
and reverse tissue utility uses a fixed logistic probe. These metrics are labeled as
mouse equivalents rather than exact human classifier-embedding reproductions.

The paper's unsupervised precision, recall, and gene-correlation agreement are
computed directly in scaled L974 space over the configured full quality cohort.
Nearest-neighbor adversarial accuracy is separately limited to 2,048 profiles,
matching the released evaluation code. Any older run evaluated with the former
2,000-profile common cap is identified by `metric_samples` in its summary and must
be regenerated before a strict paper-protocol comparison.

## Completed Run

The seed-1234 A100 run completed all 15,000 epochs and 75,000 optimizer steps in
5,987 seconds. Final summed noise MSE was 2.5354, noise MAE was 0.02255, and peak
allocated CUDA memory was 5.93 GB.

The 9,796-profile synthetic cohort produced the following locked-test and
distribution results:

- real-train to real-test tissue balanced accuracy: 0.8954;
- synthetic-train to real-test tissue balanced accuracy: 0.8687;
- gene mean, gene standard-deviation, and gene-correlation-matrix agreement:
  0.9965, 0.9437, and 0.8791;
- direct scaled-L974 precision and recall: 0.9655 and 0.8645;
- nearest-neighbor adversarial accuracy: 0.512, close to the desired 0.5;
- PCA-50 precision, recall, and Frechet distance: 0.9855, 0.9430, and 0.0385.

Synthetic reverse-validation recall matched the real-data recall for skeletal muscle
(0.9721). It was also high for cerebellum (0.9771), brown adipose (0.9588), retina
(0.9543), white adipose (0.9542), liver (0.9466), and adrenal gland (0.9455). Bone
(0.6357) and brain (0.7017) were weakest, but their real-data probe ceilings were
also the two lowest. See `evaluation/per_tissue_reverse_validation.tsv` and its
PNG/PDF plot.

This result reverses the conclusion from the deleted reduced proxy: the full
published architecture and training duration do learn tissue-conditioned mouse
ARCHS4 generation. It does not reproduce every GTEx headline metric. Direct L974
recall is 0.8645 versus the paper's reported 0.8923 for human GTEx, and correlation-
matrix agreement is lower than the paper's reported correlation score. Those values
are not strict cross-species comparisons because the mouse panel, classes, cohort,
split handling, and metric sample limit differ.

Two limitations remain important. The two-dimensional `t=0` PCA silhouette is
-0.271, so tissues overlap in the first two PCs even though the held-out classifier
recovers them. Also, 9.04% of generated scaled entries are negative. The output keeps
these values and inverse-scaled TPM explicitly unclipped; downstream use must apply a
declared nonnegativity/count-generation policy rather than silently treating the
unclipped matrix as physical TPM.
