# Metascape t1oeisjsa Summary

Source report: https://metascape.org/gp/index.html#/reportfinal/t1oeisjsa

This is the 12-muscle-outlier-filtered aggregate liver GLARE/MOBER run.
Metascape was run as a multi-list custom analysis with M. musculus as both
input and analysis species. The report used a user-supplied background of
20,979 genes, p < 0.01, minimum overlap of 3, and enrichment factor > 1.5.

## Input Lists

| List | Submitted genes | Unique genes |
| --- | ---: | ---: |
| DESeq2_meta_sig_12filter | 47 | 47 |
| FLT_cluster_13_muscle_contractile | 326 | 326 |
| FLT_cluster_14_liver_secreted | 303 | 303 |
| FLT_cluster_15_contractile_mito | 278 | 278 |
| FLT_cluster_3_mixed_stress | 1,968 | 1,965 |
| GC_cluster_1_DE_enriched | 2,185 | 2,182 |

## Main Interpretation

The strongest combined signal is still skeletal/striated-muscle biology, even
after removing the 12 highest muscle-signature profiles. This signal appears in
the DESeq2 significant genes and in the FLT cluster 13/15 gene modules.

The cleanest liver-specific GLARE result is FLT cluster 14. It is enriched for
xenobiotic metabolism, steroid metabolism, sulfur-compound metabolism, aflatoxin
activation/detoxification, ribosomal/translation terms, and PaGenBase
tissue-specific liver genes.

FLT cluster 3 and GC cluster 1 share large cilium/flagellum motility,
neuroactive ligand-receptor, GPCR, and ion-transport modules. These are
statistically enriched, but because the clusters are very large they should be
treated as secondary until checked against study, sex, strain, mission, and
library effects.

## Top Signals By List

| List | Strongest terms |
| --- | --- |
| DESeq2_meta_sig_12filter | striated muscle contraction; cytoskeleton in muscle cells; motor proteins; actomyosin organization; myofibril assembly |
| FLT_cluster_13_muscle_contractile | cytoskeleton in muscle cells; striated muscle contraction; muscle contraction; myofibril assembly; sarcomere organization |
| FLT_cluster_14_liver_secreted | xenobiotic metabolism; cellular response to xenobiotic stimulus; steroid metabolism; ribosome/translation; sulfur metabolism; aflatoxin detoxification |
| FLT_cluster_15_contractile_mito | muscle tissue/structure development; muscle cell differentiation; skeletal muscle adaptation; cytoskeleton in muscle cells |
| FLT_cluster_3_mixed_stress | cilium movement; flagellated sperm motility terms; metal ion transport; neuroactive ligand-receptor interaction; GPCR terms |
| GC_cluster_1_DE_enriched | neuroactive ligand-receptor interaction; metal ion transport; cilium movement; GPCR ligand binding; striated muscle contraction |

## QC Notes

PaGenBase confirms two important patterns:

| Term | Main matching lists |
| --- | --- |
| Tissue-specific: liver | FLT_cluster_14_liver_secreted |
| Cell-specific: c2c12 / tissue-specific: skeletal muscle | DESeq2_meta_sig_12filter, FLT_cluster_13_muscle_contractile, FLT_cluster_15_contractile_mito, GC_cluster_1_DE_enriched |

TRRUST highlights expected muscle regulators for the muscle-like lists, including
Myog, Mef2c, and Myod1, and liver/metabolic regulators for FLT cluster 14,
including Srebf1, Nr1d1/Clock, and Hnf4a.

## Recommended Use

Use FLT cluster 14 as the main biologically interpretable GLARE liver module
from this Metascape run.

Treat DESeq2_meta_sig_12filter, FLT cluster 13, and FLT cluster 15 as residual
muscle/contractile modules rather than primary spaceflight-liver biology.

Use FLT cluster 3 and GC cluster 1 only after study/batch/tissue checks, because
their enriched modules are broad and their gene lists are much larger than the
other inputs.

## Key Files

| File | Purpose |
| --- | --- |
| `AnalysisReport.html` | Local HTML version of the Metascape report |
| `metascape_result.xlsx` | Metascape workbook export |
| `Enrichment_GO/_FINAL_GO.csv` | Clustered representative enriched terms |
| `Enrichment_GO/GO_AllLists.csv` | Per-list enrichment results |
| `Enrichment_heatmap/HeatmapSelectedGO.csv` | Selected heatmap terms across lists |
| `Enrichment_PPI/GO_MCODE.csv` | Enrichment of PPI/MCODE modules |
| `Enrichment_PPI/_FINAL_MCODE.csv` | Final PPI MCODE module genes |
| `Enrichment_QC/GO_PaGenBase.csv` | Tissue/cell-specific enrichment QC |
| `Enrichment_QC/GO_TRRUST.csv` | Transcription-factor enrichment QC |
