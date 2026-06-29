# GLARE Validation Stack

This note records the validation pass run for the multi-tissue NASA mouse GLARE analysis.

## Command

```bash
PYTHONUNBUFFERED=1 PYTHONPATH=src /opt/anaconda3/envs/nasa/bin/python -m nasa_mouse_glare.multi_tissue_validation \
  --include-per-study \
  --include-mober \
  --shap-aggregate \
  --verification-estimators 80 \
  --max-eval-genes 4000 \
  --max-eac-genes 1200 \
  --terms-per-class 6 \
  --random-sets 100
```

## Output Location

The validation outputs are under:

```text
outputs/glare_multi_tissue_api/validation_stack/
```

Important files:

- `VALIDATION_STACK_SUMMARY.md`
- `xgboost_verification_summary.tsv`
- `xgboost_verification_folds.tsv`
- `representation_qc.tsv`
- `clustering_qc.tsv`
- `base_cluster_dgea_enrichment.tsv`
- `base_cluster_dgea_summary.tsv`
- `candidate_modules.tsv`
- `candidate_module_score_validation.tsv`
- `candidate_module_random_set_effects.tsv`
- `candidate_module_score_meta.tsv`
- `candidate_module_panglao_enrichment.tsv`
- `metascape_gene_lists/candidate_module_gene_lists.csv`
- `metascape_gene_lists/candidate_module_gene_list_manifest.tsv`
- `metascape_gene_lists/background_all_candidate_tissue_genes.txt`

The `outputs/` directory is ignored by git, so this note tracks the run and interpretation while the full result tables remain local/generated artifacts.

## Scope

The run completed successfully for 89 scopes:

- 12 aggregate GLARE scopes
- 12 aggregate plus MOBER scopes
- 65 per-study GLARE scopes

The validation summary includes:

- 178 XGBoost verification summary rows
- 534 representation QC rows
- 178 clustering QC rows
- 144 candidate modules tested
- 144 module-score meta-analysis rows

## Validation Components

The stack follows the GLARE paper's validation logic as closely as possible for this dataset:

1. Verification classifier: XGBoost FLT-vs-GC prediction using melted gene expression.
2. Representation QC: FT-SAE latent space compared with raw PCA using silhouette, KNN pseudo-label accuracy, and trustworthiness.
3. Consensus clustering QC: consensus clusters compared with base methods and sampled average-linkage EAC agreement.
4. DEG enrichment QC: cluster DEG enrichment compared across GMM, HDBSCAN, spectral clustering, and consensus clustering.
5. Biological enrichment setup: candidate gene lists exported for Metascape.
6. Hidden-module follow-up: GLARE-only and DGEA-intersection modules scored across studies and compared with random gene-set controls.
7. Cell-type proxy: Panglao marker enrichment run as a weak marker-based proxy.
8. SHAP: aggregate scopes produced XGBoost feature importance and SHAP outputs.

## Main Results

The verification classifier found strong learnable FLT-vs-GC signal in aggregate tissue analyses. Grouped gene-fold AUCs for aggregate scopes were:

| Tissue | Grouped AUC | Accuracy |
| --- | ---: | ---: |
| thymus | 0.997610 | 0.970082 |
| spleen | 0.994418 | 0.955099 |
| liver | 0.990744 | 0.942385 |
| kidney | 0.990576 | 0.948678 |
| skin | 0.982778 | 0.921321 |
| lung | 0.980160 | 0.916579 |
| skeletal_muscle | 0.979250 | 0.913599 |
| skeletal_muscle_soleus | 0.965613 | 0.890149 |
| skeletal_muscle_quadriceps | 0.928345 | 0.831883 |
| skeletal_muscle_gastrocnemius | 0.911346 | 0.808838 |
| skeletal_muscle_tibialis_anterior | 0.899698 | 0.797351 |
| skeletal_muscle_edl | 0.880470 | 0.790834 |

Across all 89 scopes, grouped-fold verification had mean AUC 0.893430 and median AUC 0.898143. The weakest AUCs were mostly small per-study scopes.

Representation QC was mixed. FT-SAE latent representations usually had high KNN pseudo-label accuracy and high trustworthiness, but cluster silhouettes were often weak. Aggregate MOBER scopes improved median silhouette relative to aggregate-only scopes.

Consensus clustering was not uniformly superior to base methods. GMM produced the most DEG-enriched clusters in most tissues, while consensus won only in gastrocnemius and spleen. This differs from the original GLARE paper, where EAC had the highest DEG proportion.

## GLARE-Only Versus DGEA-Intersection Signal

The strongest signal is in DGEA-intersection modules, not GLARE-only hidden modules:

| Module class | Terms tested | FDR < 0.05 | Median FDR | Median empirical p |
| --- | ---: | ---: | ---: | ---: |
| DGEA intersection | 72 | 32 | 0.081459 | 0.0050 |
| GLARE only | 72 | 17 | 0.193479 | 0.1575 |

Interpretation:

- DGEA-intersection modules are better validated and should be the main biological evidence.
- GLARE-only modules include some plausible candidates, but they are weaker overall and should be treated as exploratory.

Strong DGEA-intersection examples:

- `skeletal_muscle_soleus`: mitochondrial fatty acid beta oxidation
- `skeletal_muscle`: respiratory electron transport, protein metabolism, mitotic/cell-cycle modules
- `thymus`: DNA replication and mitotic cell-cycle modules
- `liver`: immunoregulatory and second-messenger modules

Promising GLARE-only candidates:

- `thymus`: Golgi vesicle biogenesis and platelet calcium response
- `skeletal_muscle` / `soleus`: respiratory and Cyclin E/circadian modules
- `kidney`: membrane trafficking and fatty acid/ketone metabolism

## Caveats

This run is close to the GLARE validation stack, but not identical:

- Full dense EAC over all genes was not run because it is too large for these gene sets; validation uses sampled average-linkage EAC agreement.
- Metascape gene lists were exported, but Metascape reports still require submission or automated client execution.
- No curated mouse spaceflight stress-network or transcription-factor validation database is included in this repo.
- Panglao enrichment is only a marker proxy and should not be interpreted as definitive cell-type deconvolution.
- Liver olfactory/chemosensory labels should not be automatically discarded. Olfactory receptor expression can be biologically relevant in liver, but these terms remain high-caution because large receptor gene families can dominate enrichment.
- Since GMM often outperformed consensus by DEG-enriched cluster count, we should not claim that the consensus clustering is superior in the same way as the original GLARE paper.

## Bottom Line

The validation stack supports real learnable FLT-vs-GC structure across several tissues. The most defensible biological conclusions come from modules that are supported by both DGEA and GLARE clustering. GLARE-only modules are useful for hypothesis generation, but they need Metascape, stress/TF-network, and external validation before being treated as hidden discoveries.
