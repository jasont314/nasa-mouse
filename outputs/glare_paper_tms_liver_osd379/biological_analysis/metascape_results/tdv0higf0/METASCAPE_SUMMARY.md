# Metascape Custom-Background Summary

## Configuration

- Session: `tdv0higf0`
- Analysis: Custom Analysis, multiple gene lists
- Species: *Mus musculus*
- Foreground: 26 eligible FLT/GC consensus clusters
- Submitted background: 16,553 mouse Ensembl gene IDs
- Converted background: 16,319 unique mouse Entrez genes
- Foreground genes added to background: 0
- Thresholds: p-value < 0.01, enrichment factor > 1.5, overlap >= 3

The custom background is the set of genes shared by the TMS pretraining and
OSD-379 fine-tuning inputs. This is the final enrichment run; session
`tjb92xfo7` used the whole mouse genome as background.

## Main Results

The custom-background analysis returned 8,583 significant gene-list/term
associations. The whole-genome run returned 14,448. Of the custom-background
associations, 8,304 were also found in the whole-genome run.

The strongest functional modules are:

- `FLT_cluster_12`: muscle-cell cytoskeleton, muscle contraction, striated
  muscle contraction, and myofibril assembly.
- `GC_cluster_12`: the same muscle contraction/cytoskeleton module.
- `GC_cluster_11` and `GC_cluster_10`: muscle development, differentiation,
  and actin organization.
- `FLT_cluster_14`: muscle cytoskeleton plus extracellular-matrix and elastic
  fiber assembly.
- `FLT_cluster_13`: electron transport chain and oxidative phosphorylation.
- `FLT_cluster_11`, `FLT_cluster_09`, `GC_cluster_06`, and `GC_cluster_13`:
  small-molecule, fatty-acid, monocarboxylic-acid, and xenobiotic metabolism.
- `FLT_cluster_02`: mitochondrial gene expression, RNA metabolism, splicing,
  ribosome biogenesis, and endosomal transport.
- `FLT_cluster_06`, `FLT_cluster_10`, and `GC_cluster_09`: translation,
  ribosome-associated quality control, and nonsense-mediated decay.

## Integration With Differential Expression

Functional enrichment describes what a cluster contains; it does not establish
that the function changed during spaceflight. The strongest DEG-enriched
clusters are:

| Cluster | DEGs / genes | DEG proportion | Main Metascape annotation |
| --- | ---: | ---: | --- |
| FLT 12 | 100 / 115 | 86.96% | Muscle contraction and cytoskeleton |
| GC 12 | 121 / 164 | 73.78% | Muscle contraction and cytoskeleton |
| GC 11 | 114 / 242 | 47.11% | Muscle development and actin organization |
| FLT 14 | 24 / 58 | 41.38% | Muscle cytoskeleton and extracellular matrix |
| GC 10 | 83 / 277 | 29.96% | Skeletal-muscle development |
| FLT 8 | 36 / 324 | 11.11% | Skeletal/muscle development |

FLT cluster 12 and GC cluster 12 are largely the same stable module: 105 of
the 115 FLT-cluster-12 genes map to GC cluster 12. Their separate enrichment
does not by itself demonstrate different activity between FLT and GC.

## Interpretation Caveat

The strongest DEG-enriched modules contain skeletal-muscle genes despite the
samples being labeled liver. Their expression is highly heterogeneous across
samples and is driven partly by extreme profiles. Sample-level QC and
age/collection-stratified checks are required before interpreting this module
as a liver response to spaceflight rather than tissue contamination,
dissection heterogeneity, or another technical/biological mixture.

