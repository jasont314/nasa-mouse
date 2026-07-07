# GLARE Study-Effect Visualizations

These plots score each OSDR sample by mean z-scored expression of genes in each
GLARE consensus cluster, then reduce the sample-by-module score matrix with PCA
or UMAP. Points are samples, colors are OSDR accessions, and marker shape is FLT/GC.

Tracked presentation panels are in `presentation/glare/study_effects/`.
Full per-scope coordinates and module scores are generated under ignored
`outputs/glare_multi_tissue_api/study_effects/`.

Lower accession silhouette after MOBER is consistent with reduced study/batch
separation. Positive condition silhouette indicates stronger FLT/GC separation
in the GLARE-module score space.
