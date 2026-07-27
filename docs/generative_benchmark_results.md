# Generative Benchmark Results

This report summarizes the adaptive mouse bulk RNA-seq benchmark completed through
2026-07-17. OSDR expression came from the NASA OSDR API-derived count matrix; the
old integrated raw OSDR H5 was not used. OSDR partitions were grouped by accession,
ARCHS4 partitions by GEO series, and preprocessing was fitted on training data only.

## Decision Summary

Earlier development used a scalar fidelity composite. That score is retired and is
not used for selection. Current decisions require every paper-aligned fidelity
metric, diversity, memorization, and conditional-effect gate to pass independently.

| model and arm | regime | independent result | decision |
|---|---|---|---|
| Exact Lacan ModelDDIM | ARCHS4 only | Precision/recall/F1 0.966/0.865/0.912 and AA 0.512 pass; correlation-matrix agreement 0.879 fails 0.98 | Useful tissue baseline; not an all-metric finalist |
| Factorized ModelDDIM plus train-only calibration | ARCHS4 then OSDR | Locked-test mean Corr/P/R/F1/AA 0.977/0.998/0.997/0.997/0.458; FD ratio 0.075; pooled FLT/GC recovery passes 3/4 generations | Accepted for within-study conditional generation; not unseen-study transfer |
| Paper-native Vinas WGAN-GP | ARCHS4 only | Mean/SD correlations 0.323/0.563; AA 0.950 | Rejected |
| Exact ModelDDIM condition extension | OSDR only | FLT/GC delta correlation 0.125; SD correlation 0.153 | Rejected |
| Function-preserving ModelDDIM transfer | ARCHS4 then OSDR | Training Corr/P/R/F1/AA 0.949/0.986/0.966/0.976/0.646; held-out 0.767/0.454/0.788/0.576/0.985 | Rejected |
| Liver ModelDDIM with study conditioning | ARCHS4 then OSDR | Held-out Corr/P/R/F1/AA 0.635/0.974/0.974/0.974/0.618; FLT/GC delta correlation 0.270 | Rejected |
| Liver ModelDDIM with within-study plus global z-score | OSDR only | Held-out Corr/P/R/F1/AA 0.188/0.708/1.000/0.829/0.708; FLT/GC delta correlation -0.097 | Rejected |
| Vinas WGAN-GP, CPM adaptation | OSDR only | PR F1 0.238 and AA 0.998 | Rejected |
| Vinas WGAN-GP transfer | ARCHS4 then OSDR | Pooled delta correlation 0.805, accession-aware correlation -0.022, zero of nine tissues pass | Rejected |

The exact paper-architecture ARCHS4 ModelDDIM remains the strongest broad-tissue
baseline, but it does not pass the stricter all-metric rule because its direct
correlation-matrix agreement is 0.879. Its direct scaled-L974 precision/recall are
0.966/0.865, mean/SD correlations are 0.997/0.944, and nearest-neighbor two-sample
accuracy is 0.512. It is a defensible 20-tissue ARCHS4 generator, not an OSDR
spaceflight generator. Generated scaled values include 9.04% negatives, so export
must use the documented nonnegative clipping policy before treating inverse-scaled
values as TPM.

The factorized DDIM is the first OSDR FLT/GC generator to pass the fixed broad
locked-test rule. It passed every finite-sample fidelity metric in four of four test
generations and pooled FLT/GC recovery in three of four. Skeletal-muscle
accession-aware recovery also passed in four of four test generations, although
joint real/synthetic LOO-stable gene hits remained 0, 0, 0, and 1. The mean Corr of
0.977 passes the test-size real-bootstrap floor of 0.950 but remains below the
separate strict 0.98 paper benchmark.

This result is limited to within-study interpolation: the 781/536/293-profile split
places samples from each eligible accession in train, validation, and test strata.
The train-only calibrator does not fit FLT/GC effects, but it does use known
accession/tissue means. It therefore does not establish generation for a novel
study. Real-plus-synthetic classifier training also failed to improve over real-only
training (balanced accuracy 0.734 versus 0.754), so augmentation is not promoted.
See `osdr_conditional_diffusion_finalist.md` for the fixed protocol and full result.

## Matched Liver Harmonization Benchmark

Nine preprocessing arms were trained for 15,000 epochs with the same 225.6-million-
parameter conditional ModelDDIM, random seed, 974 genes, and API-derived liver
partition. All prepared matrices have identical source-row and gene hashes: 119
training profiles, 50 validation profiles from OSD-137/457/48, and 70 locked test
profiles from OSD-379. The test partition was not evaluated.

| preprocessing / harmonization | Corr | precision | recall | F1 | AA | FD / real reference | FLT/GC gate | accession gate |
|---|---:|---:|---:|---:|---:|---:|---|---|
| No harmonization, TPM | 0.283 | 0.200 | 0.960 | 0.331 | 0.850 | 2.052 | fail | fail |
| Ilangovan per-study z-score | 0.278 | 0.160 | 1.000 | 0.276 | 0.770 | 0.716 | fail | fail |
| Mentor two-stage z-score | 0.348 | 0.440 | 1.000 | 0.611 | 0.690 | 0.977 | fail | fail |
| ComBat by study | 0.004 | 0.040 | 1.000 | 0.077 | 0.810 | 63.185 | fail | fail |
| ComBat-seq by study | 0.067 | 0.020 | 1.000 | 0.039 | 0.850 | 156.647 | fail | fail |
| MBatch Median Polish by study | 0.009 | 0.020 | 1.000 | 0.039 | 0.930 | 205.470 | fail | fail |
| MBatch Empirical Bayes by study | 0.001 | 0.000 | 1.000 | 0.000 | 0.850 | 60.625 | fail | fail |
| MBatch ANOVA by study | -0.003 | 0.020 | 1.000 | 0.039 | 0.870 | 44.716 | fail | fail |
| MOBER by study | 0.808 | 0.260 | 1.000 | 0.413 | 0.770 | 33.311 | fail | fail |

No arm passes all six independent fidelity criteria, and no arm passes either
condition-effect requirement. The mentor two-stage transform is the closest balanced
fidelity result, but it still fails Corr, precision, F1, and AA. MOBER's high Corr
does not compensate for its low precision/F1, separability, or 33.3-fold FD ratio.
The uniformly high recall with very low precision in the batch-correction arms is
consistent with synthetic support that is much broader than the real validation
distribution, not successful generation.

ComBat, ComBat-seq, and the three MBatch held-out transforms use training anchors
but remain transductive sensitivity analyses. MOBER is the only complex inductive
harmonizer in this table. Machine-readable results and independent-metric plots are
under `outputs/generative_benchmark/summary/liver_harmonization/`.

## Adaptive Comparisons

- **ARCHS4 pretraining:** the exact tissue model plus factorized OSDR adaptation and
  train-only calibration is now accepted for within-study conditional generation.
  It does not solve unseen-accession transfer or improve classifier augmentation.
- **Preprocessing:** CPM plus `log1p` and gene z-score improved direct WGAN over the
  strict paper transform (0.468 versus 0.424), but both failed fidelity. Paper-native
  full-transcriptome TPM plus train-fitted MaxAbs worked for broad ModelDDIM.
- **Study conditioning and harmonization:** all requested matched liver arms now
  have full 15,000-epoch evaluations: no correction, Ilangovan normalization,
  mentor two-stage scaling, ComBat, ComBat-seq, all three official MBatch methods,
  and MOBER. None passes all independent fidelity and effect-recovery gates. Exact
  liver study conditioning separately reduced held-out AA from 0.985 to 0.618 but
  also failed Corr and FLT/GC recovery.
- **Per-tissue behavior:** the WGAN transfer failed FLT/GC recovery in all nine
  evaluable validation tissues. Small skeletal-muscle DDIM specialists were less
  stable than the pooled factorized model and were rejected. The pooled model's
  muscle accession diagnostic passed on locked test, but exact LOO-stable gene
  overlap remained negligible.

## GeneJEPA

GeneJEPA remains representation-only and has no expression decoder. The exact
768-wide, 512-latent, 24-block screen completed 43,744 ARCHS4 exposures in 1,074
training seconds with 31.68 GB peak allocated A100 memory. On the balanced
held-out-series figure cohort, its tissue balanced accuracy/macro F1 were
0.703/0.701 versus 0.839/0.840 from expression. Embedding and UMAP silhouettes were
-0.176/-0.215. This improves over the practical model but fails the predeclared
guidance gate, so no GeneJEPA-guided diffusion experiment was started. The bounded
run made 238 optimizer updates, below the paper's 2,000-step EMA warmup, and is not
presented as a fully trained paper reproduction.

## Reproducibility

The ranked table is at `outputs/generative_benchmark/scoreboard.tsv`. Each model run
stores a resolved configuration, source and data identities, split hashes, fitted
preprocessing, device/runtime records, model checkpoint, validation summary, and
plots under `outputs/generative_benchmark/runs/`. The experiment matrix is resumable
and content-addressed. The finalist protocol implementation passed all 111
repository tests.

See `generative_models_pipeline.md`, `generative_benchmark_decisions.md`,
`rna_diffusion_paper_parity.md`, and `osdr_conditional_diffusion_finalist.md` for
implementation details and experiment history.
