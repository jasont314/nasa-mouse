# Generative Benchmark Results

This report summarizes the adaptive mouse bulk RNA-seq benchmark completed on
2026-07-16. OSDR expression came from the NASA OSDR API-derived count matrix; the
old integrated raw OSDR H5 was not used. OSDR partitions were grouped by accession,
ARCHS4 partitions by GEO series, and preprocessing was fitted on training data only.

## Decision Summary

Earlier development used a scalar fidelity composite. That score is retired and is
not used for selection. Current decisions require every paper-aligned fidelity
metric, diversity, memorization, and conditional-effect gate to pass independently.

| model and arm | regime | independent result | decision |
|---|---|---|---|
| Exact Lacan ModelDDIM | ARCHS4 only | Precision/recall/F1 0.966/0.865/0.912 and AA 0.512 pass; correlation-matrix agreement 0.879 fails 0.98 | Useful tissue baseline; not an all-metric finalist |
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

No OSDR FLT/GC generator passes all independent fidelity gates. The OSDR locked
test therefore remains unopened and augmentation is blocked. Pooled FLT/GC recovery
alone is not accepted:
promotion now also requires accession-aware meta-effect correlation at least 0.30
and direction agreement at least 0.55.

## Adaptive Comparisons

- **ARCHS4 pretraining:** improved several exact-DDIM distribution and pooled
  condition metrics, but did not solve real-versus-synthetic separation on unseen
  accessions. There is no accepted conditional benefit yet.
- **Preprocessing:** CPM plus `log1p` and gene z-score improved direct WGAN over the
  strict paper transform (0.468 versus 0.424), but both failed fidelity. Paper-native
  full-transcriptome TPM plus train-fitted MaxAbs worked for broad ModelDDIM.
- **Study conditioning and harmonization:** ComBat, ComBat-seq, and MOBER adapters
  passed integration smoke tests only. Exact liver study conditioning reduced
  held-out AA from 0.985 to 0.618 but failed Corr and FLT/GC recovery. A full
  mentor-style `log1p` TPM, within-study z-score, then pooled z-score run reached
  held-out AA 0.708 but failed Corr, precision, F1, and condition-effect recovery.
  There is no decision-quality evidence that batch harmonization has produced an
  acceptable generator.
- **Per-tissue behavior:** the WGAN transfer failed FLT/GC recovery in all nine
  evaluable validation tissues. Standalone per-tissue expansion was stopped because
  the pooled candidate failed the fixed gate.

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
and content-addressed. The final validation run passed all 82 repository tests.

See `generative_models_pipeline.md`, `generative_benchmark_decisions.md`, and
`rna_diffusion_paper_parity.md` for implementation details and experiment history.
