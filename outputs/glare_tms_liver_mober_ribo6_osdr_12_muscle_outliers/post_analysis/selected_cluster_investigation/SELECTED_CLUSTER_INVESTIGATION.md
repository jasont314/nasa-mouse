# Selected Cluster Investigation

Run: `glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers`

This checks whether selected GLARE modules are higher in flight and whether
FLT cluster 3 / GC cluster 1 look biologically interpretable after batch checks.
Module scores are mean processed training expression over the genes in each
module, evaluated in both FLT and GC sample matrices.

## Directionality

FLT cluster 14 is liver-enriched, but it is not higher in flight as a module.

| Module | Genes | FLT mean | GC mean | FLT-GC | Accession-adjusted FLT-GC | Adjusted p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| FLT_cluster_14_liver_secreted | 303 | 5.0064 | 5.0453 | -0.0389 | -0.0254 | 0.4093 |
| FLT_cluster_3_mixed_stress | 1,968 | 0.1402 | 0.1424 | -0.0021 | 0.0011 | 0.8728 |
| GC_cluster_1_DE_enriched | 2,185 | 0.1526 | 0.1441 | 0.0086 | 0.0117 | 0.1418 |
| FLT3_GC1_overlap | 1,874 | 0.1364 | 0.1378 | -0.0014 | 0.0017 | 0.8026 |

Interpretation: the selected modules are not convincingly shifted up or down as
whole modules after accounting for OSD accession. FLT cluster 14 is a
flight-trained latent cluster with liver/metabolism enrichment, not a
flight-up expression module.

## FLT Cluster 14

Metascape enrichment is liver-like:

- xenobiotic metabolic process
- cellular response to xenobiotic stimulus
- steroid metabolic process
- sulfur compound metabolic process
- aflatoxin activation and detoxification
- ribosome/translation terms
- PaGenBase tissue-specific liver

Raw-count DESeq2/meta-analysis finds only 4 significant genes in this 303-gene
module:

| Gene | Meta log2FC | FDR | Direction across studies |
| --- | ---: | ---: | --- |
| Apoa4 | 1.5906 | 8.84e-39 | mixed, 5 up / 1 down |
| Cdkn1a | 1.1292 | 3.76e-08 | mixed, 4 up / 2 down |
| Mup17 | -1.1798 | 7.73e-08 | mixed, 2 up / 4 down |
| Mup15 | -1.5347 | 3.23e-06 | all_down, 0 up / 5 down |

So the biological content is plausible liver metabolism, but it is not a
uniformly flight-induced module.

## FLT Cluster 3 And GC Cluster 1

These two modules are largely the same gene set:

| Query | Genes | Top matching cluster | Overlap | Overlap fraction |
| --- | ---: | ---: | ---: | ---: |
| FLT cluster 3 to GC | 1,968 | GC cluster 1 | 1,874 | 0.9522 |
| GC cluster 1 to FLT | 2,185 | FLT cluster 3 | 1,874 | 0.8577 |

Metascape enrichment is dominated by:

- cilium movement / flagellated sperm motility
- neuroactive ligand-receptor interaction
- GPCR signaling
- metal ion transport / membrane potential
- for GC cluster 1, additional striated muscle contraction terms

This does not look like a clean liver-specific spaceflight module.

## Batch And Study Checks

The module is strongly study/batch-associated.

| Module | Strongest batch variable | Eta-squared | FDR | Strong driver |
| --- | --- | ---: | ---: | --- |
| FLT cluster 3 | OSD accession / mission | 0.5650 | 2.10e-10 | yes |
| FLT cluster 3 | strain | 0.5482 | 1.43e-10 | yes |
| FLT cluster 3 | sex | 0.3869 | 2.28e-08 | yes |
| GC cluster 1 | OSD accession / mission | 0.3135 | 1.61e-04 | yes |
| GC cluster 1 | strain | 0.2570 | 6.06e-04 | yes |
| FLT cluster 14 | OSD accession / mission | 0.2755 | 6.27e-04 | yes |
| FLT cluster 14 | strain | 0.2656 | 3.69e-04 | yes |

The FLT3/GC1 module score also varies more by study than by condition. For
example, FLT cluster 3 is highest in OSD-245 and OSD-379, lower in OSD-173 and
OSD-463, and has mixed FLT-GC directions by study.

## DESeq2 Context

FLT cluster 3 has 13 significant DESeq2/meta genes. They include a mixed set:
`Tg`, `Mucl1`, `Trpm5`, `Wnt11`, `Fgf10` are flight-up in available studies,
while `Il6`, `Ptgs2`, `Cxcl2`, and `Cdk3` are flight-down.

GC cluster 1 has 29 significant DESeq2/meta genes, but the strongest ones include
the same residual muscle/contractile genes already flagged elsewhere:
`Myl1`, `Tnnt3`, `Acta1`, `Myh4`, `Myh2`, `Mybpc1`, `Myh1`, plus `Tg`,
`Cox8b`, and `Mucl1`.

Important label note: `GC_cluster_1_DE_enriched` means this cluster was learned
from the GC representation and is enriched for DESeq2-significant genes. It does
not mean the genes are higher in GC.

## Conclusion

Use FLT cluster 14 cautiously as a liver/metabolism module, not as evidence that
the whole module is flight-up.

Treat FLT cluster 3 and GC cluster 1 as a shared broad/batch-sensitive module.
It may contain real biology in individual genes, but the full module is too
study-, strain-, sex-, and residual-muscle-associated for a primary
spaceflight-liver interpretation.

## Generated Tables

| File | Purpose |
| --- | --- |
| `selected_module_score_directionality.tsv` | Overall FLT vs GC module score tests |
| `selected_module_score_by_study.tsv` | Per-accession module score means |
| `selected_module_scores_by_sample.tsv` | Sample-level module scores |
| `selected_cluster_overlap.tsv` | FLT/GC cluster overlap checks |
| `selected_cluster_batch_qc.tsv` | Batch-driver rows for selected clusters |
| `selected_cluster_dgea_genes.tsv` | DESeq2/meta rows for all selected-cluster genes |
| `selected_cluster_significant_dgea_genes.tsv` | Significant DESeq2/meta genes in selected clusters |
| `selected_cluster_metascape_top_terms.tsv` | Top Metascape terms for selected clusters |
