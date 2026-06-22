# Filtered OSD-379 GLARE and DESeq2 Run

## Tissue-Composition Filter

- Input study: OSD-379/RR-8 left-lobe liver RNA-seq.
- Direct exclusion rule: at least 10/20 canonical skeletal-muscle markers above
  100 normalized counts and at least 0.005% total marker-panel abundance.
- Excluded: 14 profiles.
- Retained: 26 flight and 30 ground-control profiles.
- Clean animals in the opposite condition were retained; filtering was
  independent rather than matched-slot deletion.

The final DESeq2 model also includes the within-stratum standardized muscle
score to adjust for residual composition below the hard-filter threshold.

## DESeq2

- Input: official RSEM unnormalized integer counts.
- Global design: `~ stratum + muscle_score_z + condition`.
- Strata: ISS-terminal young/old and live-animal-return young/old.
- Prefilter: count >= 10 in at least three retained samples.
- Tested globally: 17,112 genes.
- Significant globally: 96 genes at adjusted p < 0.05 and absolute log2 fold
  change >= 1; 60 flight-up and 36 flight-down.
- None of the 20 canonical muscle markers is significant.

The global adjusted result is primary. Stratum-specific results are secondary,
especially live-animal-return old, which has four flight and five ground
profiles after tissue QC.

Top biological signals include circadian regulation (`Npas2`, `Bmal1`), lipid
and energy metabolism, and genes including `Dhrs9`, `Tubb2a`, `Upp2`, `Chka`,
`Acmsd`, `Cyp2a5`, and `Fabp5`.

## GLARE

- Pretraining weights: existing TMS FACS liver 16-dimensional SAE.
- Gene universe: the same 16,553 genes as the original controlled run.
- Fine-tuning: separate FLT and GC adapters, 30 epochs, seed 1996.
- FLT: 26 profiles, best loss 0.15614395 at epoch 30.
- GC: 30 profiles, best loss 0.15580788 at epoch 28.
- FLT consensus: 16 clusters, silhouette 0.117.
- GC consensus: 15 clusters, silhouette 0.353.
- Filtered FLT-vs-GC partitions: ARI 0.423, NMI 0.671.

The XGBoost verification stage was not run because independent QC leaves
unequal FLT/GC feature dimensions. Removing clean counterpart animals only to
equalize dimensions would be methodologically unjustified.

## DESeq2 vs GLARE

- Of 96 global DEGs, 48 occur in the controlled GLARE gene universe.
- No complete GLARE consensus cluster is DEG-enriched at BH FDR < 0.05.
- DEGs have larger Procrustes-aligned FLT-vs-GC latent shifts than non-DEGs:
  Mann-Whitney p = 1.60e-8 and ROC AUC = 0.731.
- Absolute DESeq2 effect and latent shift are weakly correlated:
  Spearman rho = 0.074.
- The original muscle-heavy FLT cluster 12 splits across 11 filtered clusters,
  and none of its 115 genes is a composition-adjusted DEG.
- Filtered-vs-original cluster agreement is moderate: FLT ARI 0.386 and GC ARI
  0.413, showing that tissue composition materially affected the old result.

Reactome enrichment gives complementary views:

- DESeq2 emphasizes circadian clock, lipid metabolism, and energy metabolism.
- Top GLARE latent shifts emphasize platelet/hemostasis structure, glucose
  metabolism, and lipid metabolism.
- Six pathways are significant in both, including lipid/lipoprotein
  metabolism, fatty-acyl-CoA synthesis, triglyceride biosynthesis, integrated
  energy metabolism, and small-molecule transport.

The methods therefore agree most clearly at the metabolic-pathway level.
DESeq2 identifies sparse mean-expression changes, while GLARE captures broader
distributed representation changes rather than one DEG-dense cluster.
