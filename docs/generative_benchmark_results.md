# Generative Benchmark Results

This report summarizes the adaptive mouse bulk RNA-seq benchmark completed on
2026-07-16. OSDR expression came from the NASA OSDR API-derived count matrix; the
old integrated raw OSDR H5 was not used. OSDR partitions were grouped by accession,
ARCHS4 partitions by GEO series, and preprocessing was fitted on training data only.

## Decision Summary

| model and arm | regime | held-out composite | main result | decision |
|---|---|---:|---|---|
| Exact Lacan ModelDDIM | ARCHS4 only | 0.957 | Tissue BA 0.869 from synthetic training versus 0.895 from real training | Broad-tissue generator accepted |
| Paper-native Vinas WGAN-GP | ARCHS4 only | 0.380 | Mean/SD correlations 0.323/0.563; two-sample accuracy 0.950 | Rejected |
| Exact ModelDDIM condition extension | OSDR only | 0.367 | FLT/GC delta correlation 0.125; SD correlation 0.153 | Rejected |
| Exact ModelDDIM transfer | ARCHS4 then OSDR | 0.606 | Mean/SD correlations 0.937/0.841 but two-sample accuracy 0.973 | Rejected |
| Function-preserving ModelDDIM transfer | ARCHS4 then OSDR | 0.610 | Pooled FLT/GC delta correlation 0.542 but two-sample accuracy 0.967 | Rejected |
| Vinas WGAN-GP, CPM adaptation | OSDR only | 0.468 | Better than strict preprocessing, but PR F1 0.238 and two-sample accuracy 0.998 | Rejected |
| Vinas WGAN-GP transfer | ARCHS4 then OSDR | 0.437 | Pooled delta correlation 0.805, accession-aware correlation -0.022, zero of nine tissues pass | Rejected |

Only the exact paper-architecture ARCHS4 ModelDDIM passes the fixed fidelity,
diversity, and memorization gates. Its direct scaled-L974 precision/recall are
0.966/0.865, mean/SD correlations are 0.997/0.944, and nearest-neighbor two-sample
accuracy is 0.512. It is a defensible 20-tissue ARCHS4 generator, not an OSDR
spaceflight generator. Generated scaled values include 9.04% negatives, so export
must use the documented nonnegative clipping policy before treating inverse-scaled
values as TPM.

No OSDR FLT/GC generator passes the minimum 0.70 fidelity composite. The OSDR
locked test therefore remains unopened, augmentation is blocked, and no candidate
is promoted to additional seeds. Pooled FLT/GC recovery alone is not accepted:
promotion now also requires accession-aware meta-effect correlation at least 0.30
and direction agreement at least 0.55.

## Adaptive Comparisons

- **ARCHS4 pretraining:** improved exact DDIM OSDR composite from 0.367 to about
  0.61 and improved pooled condition recovery, but did not solve real-versus-
  synthetic separation. It reduced WGAN composite from 0.468 to 0.437. There is no
  accepted conditional benefit.
- **Preprocessing:** CPM plus `log1p` and gene z-score improved direct WGAN over the
  strict paper transform (0.468 versus 0.424), but both failed fidelity. Paper-native
  full-transcriptome TPM plus train-fitted MaxAbs worked for broad ModelDDIM.
- **Study conditioning and harmonization:** ComBat, ComBat-seq, MOBER, and explicit
  study-conditioning paths passed integration smoke tests. They were not promoted
  to full-duration comparisons because the prerequisite pooled generators failed;
  no decision-quality evidence says that any helped generation.
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
