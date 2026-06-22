# Controlled OSD-379 GLARE Run

## Data

- Pretraining: 21,010 shared genes by 2,859 TMS FACS liver cells from 9 mice.
- Fine-tuning: OSD-379/RR-8 official normalized counts.
- Target intersection: 16,553 genes.
- Conditions: 35 flight and 35 ground-control liver profiles.
- Matched strata: ISS-terminal young/old and live-animal-return young/old.

## Model

- SAE architecture: 2,859 -> 128 -> 64 -> 32 -> 16.
- TMS pretraining best loss: 0.861717 at epoch 28/30.
- FLT fine-tuning best loss: 0.154424 at epoch 28/30.
- GC fine-tuning best loss: 0.154984 at epoch 30/30.
- Seed: 1996, reset identically for FLT and GC.

## Verification

- GLARE shuffled five-fold XGBoost ROC AUC: 0.9797 +/- 0.0026.
- Gene-grouped five-fold XGBoost ROC AUC: 0.9802 +/- 0.0018.
- SHAP values, feature rankings, gene-condition rankings, and a beeswarm plot
  are under `verification/`.

The grouped split keeps both condition rows for each gene in one fold. This
rules out paired-row leakage, but the task is still GLARE's gene-level
condition verification rather than held-out-mouse prediction.

## Clustering

- FLT: GMM 20, Spectral 25, EAC consensus 16; silhouette 0.201.
- GC: GMM 25, Spectral 20, EAC consensus 15; silhouette 0.129.
- FLT-vs-GC consensus agreement: ARI 0.489, NMI 0.625.

EAC is the paper's co-association distance followed by average linkage. It is
computed with lossless partition-signature compression instead of allocating a
full gene-by-gene matrix.

At the paper's HDBSCAN settings, every normalized mouse gene was labeled
noise in both conditions. The implementation records a neutral constant
partition; because that adds the same co-association value to every pair, it
does not alter average-linkage merge ordering. No HDBSCAN thresholds were
tuned after seeing the result.

## Biology

NASA's official DESeq2 results were filtered to four matched Space Flight vs
Ground Control contrasts. Their strict DEG counts are:

- ISS-terminal, young: 10.
- ISS-terminal, old: 471.
- Live-animal-return, young: 102.
- Live-animal-return, old: 93.
- Union: 636 genes; 441 flight-up, 174 flight-down, and 21 with mixed
  direction across significant strata.

Clusters with the largest DEG proportions include FLT cluster 12 (100/115,
86.96%), GC cluster 12 (121/164, 73.78%), and GC cluster 11 (114/242,
47.11%). FLT cluster 15 contains five genes and is 100% DEG, but is excluded
from enrichment because it is below the 10-gene threshold.

Metascape Custom Analysis session `tdv0higf0` completed with 26 eligible
cluster lists and the 16,319-gene converted custom background. The final
tables and interpretation are under
`biological_analysis/metascape_results/tdv0higf0/`. Clusters over 3,000 genes
and under 10 genes are excluded from whole-cluster enrichment.

The strongest DEG-enriched clusters are dominated by muscle contraction,
cytoskeleton, and muscle-development annotations. Because the input is liver
and the muscle-gene signal is highly heterogeneous across samples, this result
requires sample-level QC before it is interpreted as a spaceflight-induced
liver program.

## Remaining Deviations

- TMS mouse liver replaces Arabidopsis root-tip single-cell data.
- OSD-379 replaces OSD-120 and has age/collection strata rather than plant
  genotype/light strata.
- GLARE's three manually selected Arabidopsis outlier genes are not removed;
  mouse PCA/k-means audit tables are exported instead.
- HDBSCAN is neutral in this run because the paper's thresholds produce only
  noise on the mouse latent spaces.
- Metascape web execution is pending; all required upload files are present.
