# Accession Split Plans

All splits group complete OSDR accessions and were assigned without reading
expression values. `pooled_accession_split.tsv` is the pooled model plan;
`per_tissue_locked_accession_splits.tsv` covers confirmatory tissues; and
`per_tissue_loo_accession_folds.tsv` defines confirmatory and exploratory LOO folds.

Locked-test accessions must not be used for preprocessing, feature selection,
checkpoint selection, or hyperparameter tuning.
