# Filtered DESeq2 vs GLARE

- DESeq2 tested 12,962 genes represented
  in GLARE and called 48 significant at
  adjusted p < 0.05 and absolute log2 fold change >= 1
  (32 up, 16 down).
- 0 FLT/GC consensus clusters are
  enriched for DEGs by one-sided Fisher test with BH FDR < 0.05.
- Absolute DESeq2 effect and aligned GLARE latent shift have Spearman rho
  0.074
  (p=3.28e-17).
- GLARE latent shift separates DEGs from non-DEGs with ROC AUC
  0.731; the one-sided Mann-Whitney p-value
  is 1.6e-08.
- FLT and GC consensus partitions have ARI 0.423 and NMI
  0.671.
- Reactome ORA found 27 significant terms for
  mapped DESeq2 DEGs and 30 for an equal-sized set
  of top GLARE-shift genes; 6
  terms are shared.
- The original muscle-heavy FLT cluster 12 now spans
  11 filtered FLT clusters and contains
  0
  composition-adjusted DEGs.
- Filtered-vs-original clustering ARI is
  0.386 for FLT
  and 0.413 for GC.

DESeq2 and GLARE answer different questions. DESeq2 tests sample-level mean
expression changes after adjusting for age/collection stratum. GLARE groups
genes by nonlinear expression representations and can capture coordinated
structure even when individual genes do not pass a DEG threshold. Agreement is
therefore assessed through DEG enrichment within GLARE clusters and through
the relationship between DESeq2 effect size and aligned latent-space shift.
