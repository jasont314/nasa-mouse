# Reviewer robustness analysis

This analysis was specified before running the reviewer-directed robustness outputs. It does not retroactively make the original exploratory pathway review prespecified.

## Frozen decisions

- Same query samples, 2,000-HVG gene universe, Reactome memberships, and decoder-oriented primary expiMap scores as the manuscript.
- Conventional benchmarks: rank-normalized ssGSEA per sample and project-balanced log2-CPM preranked GSEA.
- Internal validation: leave-one-project-out direction prediction, with the top decile selected using training projects only.
- Composition sensitivity: broad compartment markers derived from the independent Tabula Muris Senis Smart-seq2 atlas, followed by within-accession regression on marker-score principal components.
- These are triangulation and sensitivity analyses. They are not an external replication cohort, causal adjustment, or cell-type deconvolution.

## Results

### Thymus

Across active programs, expiMap and ssGSEA had Spearman r=0.25 and 58% directional agreement. Among the curated programs, expiMap agreed in direction with ssGSEA for 86% and with preranked GSEA for 86%.
Training-only top-decile pathways predicted the held-out project direction in 99% of fold-pathway comparisons. After adjustment for 3 atlas-derived broad compartment scores, 86% of curated expiMap directions were retained.

### Skin

Across active programs, expiMap and ssGSEA had Spearman r=0.24 and 60% directional agreement. Among the curated programs, expiMap agreed in direction with ssGSEA for 88% and with preranked GSEA for 100%.
Training-only top-decile pathways predicted the held-out project direction in 75% of fold-pathway comparisons. After adjustment for 3 atlas-derived broad compartment scores, 100% of curated expiMap directions were retained.

### Liver

Across active programs, expiMap and ssGSEA had Spearman r=0.22 and 57% directional agreement. Among the curated programs, expiMap agreed in direction with ssGSEA for 86% and with preranked GSEA for 100%.
Training-only top-decile pathways predicted the held-out project direction in 75% of fold-pathway comparisons. After adjustment for 4 atlas-derived broad compartment scores, 86% of curated expiMap directions were retained.

### Soleus

Across active programs, expiMap and ssGSEA had Spearman r=0.12 and 54% directional agreement. Among the curated programs, expiMap agreed in direction with ssGSEA for 67% and with preranked GSEA for 67%.
Training-only top-decile pathways predicted the held-out project direction in 28% of fold-pathway comparisons. After adjustment for 5 atlas-derived broad compartment scores, 57% of curated expiMap directions were retained.

## Interpretation limits

Project-wise cross-validation reuses the same repository and cannot replace validation in newly generated missions. The atlas marker analysis is a proxy sensitivity: composition can be a biological mediator of spaceflight, the atlas lacks some mature cell states (notably mature myofibers), and removing marker-associated variation can remove real tissue response as well as composition bias.
