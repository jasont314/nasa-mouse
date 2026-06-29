# OntoVAE Pipeline

This note documents the OntoVAE adaptation for NASA mouse OSDR bulk RNA-seq.
It is parallel to the existing expiMap runs and does not overwrite expiMap
outputs.

## Sources Inspected

- Local OntoVAE package: `src/onto-vae/`
- Local OntoVAE model API: `src/onto-vae/onto_vae/vae_model.py`
- Local method inventory: `src/METHOD_REPOS.md`
- Upstream package: https://github.com/hdsu-bioquant/onto-vae
- Paper: https://academic.oup.com/bioinformatics/article/39/6/btad387/7199588

The vendored package expects an `Ontobj`-like object containing a sample-by-gene
expression matrix, ontology masks, retained genes, and annotation tables. The
core `OntoVAE` model uses an unconstrained encoder and an ontology-masked
decoder, trains with MSE reconstruction plus KL loss, masks invalid decoder
gradients, and clamps ontology decoder weights to be nonnegative.

## Adaptation

The implemented runner builds a flat Reactome mouse decoder mask from the
existing AnnData `varm["I"]` architecture. Each retained Reactome pathway is a
latent program connected to its member genes. This uses the official mouse
Reactome GMT already generated for the project:

`data/pathways/reactome_current_mouse_ensembl.gmt`

OntoVAE does not provide a native scArches-style query-mapping mode, query
extension nodes, HSIC regularization, or de novo query programs. The closest
defensible reference-query equivalent implemented here is:

1. Pretrain OntoVAE on tissue-matched ARCHS4 mouse reference samples.
2. Project OSDR samples through the pretrained model and save those scores.
3. Fine-tune the pretrained model on OSDR FLT/GC samples and save final scores.

A direct OSDR-only baseline is kept separate. HVG variants use the existing
HVG-filtered AnnData inputs from the expiMap tutorial-style workflows.

Expression input is `log1p(CPM)` followed by z-scoring. Reference-pretrain runs
use reference-derived gene means/standard deviations for both reference and
query. Direct OSDR runs use query-derived gene scaling. This is suitable for
bulk expression matrices, but it is not a raw-count negative-binomial model.

OSDR inputs are the API-derived tissue AnnData files produced by this repo's
NASA OSDR API workflow. The raw integrated OSDR H5 file is not used.

## Scripts

- `src/nasa_mouse_glare/train_ontovae.py`
  - trains one OntoVAE run
  - writes latent/program scores, top decoder gene weights, model checkpoints,
    normalization statistics, and a JSON training summary
- `src/nasa_mouse_glare/run_ontovae_pipeline.py`
  - runs the per-tissue batch and calls the existing FLT/GC analysis and
    accession-validation helpers
- `src/nasa_mouse_glare/summarize_ontovae_results.py`
  - aggregates completed OntoVAE outputs and compares against available
    expiMap random-effects validation outputs

Run command used for the completed batch:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m nasa_mouse_glare.run_ontovae_pipeline \
  --include-reference \
  --include-direct \
  --include-hvg \
  --include-muscle-splits \
  --reference-epochs 60 \
  --query-epochs 60 \
  --batch-size 256 \
  --min-term-genes 5
```

Summary command:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m nasa_mouse_glare.summarize_ontovae_results \
  --search-root outputs \
  --output-dir outputs/ontovae_pipeline/summary
```

## Outputs

Main summary outputs:

- `outputs/ontovae_pipeline/summary/ontovae_run_manifest.tsv`
- `outputs/ontovae_pipeline/summary/ontovae_run_summary.tsv`
- `outputs/ontovae_pipeline/summary/ontovae_top_terms.tsv`
- `outputs/ontovae_pipeline/summary/ontovae_top_gene_weights_for_top_terms.tsv`
- `outputs/ontovae_pipeline/summary/ontovae_vs_expimap_random_effects.tsv`

Per-run outputs are under:

- `outputs/ontovae_<tissue>/direct_osdr/`
- `outputs/ontovae_<tissue>/archs4_pretrain_osdr_finetune/`
- `outputs/ontovae_<tissue>/hvg_archs4_pretrain_osdr_finetune/`
- `outputs/ontovae_skeletal_muscle_splits/<group>/direct_osdr/`
- `outputs/ontovae_skeletal_muscle_splits/<group>/archs4_pretrain_osdr_finetune/`

Each analyzed score set includes:

- `pathway_scores.tsv` or `pretrained_query_pathway_scores.tsv`
- `analysis/flt_vs_gc_pathway_comparison.tsv`
- `analysis/pathway_score_pca.png`
- `analysis/pathway_score_pca_by_accession.png`
- `analysis/pathway_score_umap.png`
- `analysis/pathway_score_umap_by_accession.png`
- `analysis/top_pathway_shift_heatmap.png`
- `accession_validation/random_effects_meta_analysis.tsv`
- `accession_validation/leave_one_out_summary.tsv`

The completed batch wrote 34 model runs, 55 analyzed score sets, and 280 plot
PNGs. All training summaries report CUDA on `NVIDIA A100-SXM4-40GB`.

## Completed Tissues

Main tissues:

- liver
- skeletal_muscle
- skin
- kidney
- thymus
- spleen
- lung
- retina

Skeletal-muscle split groups:

- soleus
- gastrocnemius
- quadriceps
- edl
- tibialis_anterior

## Result Summary

The strongest strict leave-one-accession-out stable OntoVAE signals are:

| tissue | group | mode | score set | LOO-stable hits | top stable term | effect |
| --- | --- | --- | --- | ---: | --- | ---: |
| liver |  | ARCHS4 pretrain + OSDR fine-tune | final | 13 | `R-MMU-5676590_NIK_NONCANONICAL_NF_KB_SIGNALING` | 0.277 |
| liver |  | direct OSDR | final | 2 | `R-MMU-3000170_SYNDECAN_INTERACTIONS` | -0.142 |
| liver |  | HVG ARCHS4 pretrain + fine-tune | final | 3 | `R-MMU-1266695_INTERLEUKIN_7_SIGNALING` | 0.166 |
| skeletal_muscle |  | ARCHS4 pretrain + OSDR fine-tune | final | 1 | `R-MMU-5576893_PHASE_2_PLATEAU_PHASE` | 0.289 |
| skeletal_muscle |  | HVG ARCHS4 pretrain + fine-tune | final | 2 | `R-MMU-73927_DEPURINATION` | 0.125 |
| skeletal_muscle | quadriceps | ARCHS4 pretrain + OSDR fine-tune | final | 3 | `R-MMU-399956_CRMPS_IN_SEMA3A_SIGNALING` | 0.544 |
| skeletal_muscle | soleus | ARCHS4 pretrain + OSDR fine-tune | final | 27 | `R-MMU-629597_HIGHLY_CALCIUM_PERMEABLE_NICOTINIC_ACETYLCHOLINE_RECEPTORS` | 0.935 |
| spleen |  | ARCHS4 pretrain + OSDR fine-tune | final | 127 | `R-MMU-194138_SIGNALING_BY_VEGF` | 0.592 |
| spleen |  | direct OSDR | final | 1 | `R-MMU-8878171_TRANSCRIPTIONAL_REGULATION_BY_RUNX1` | -0.081 |
| spleen |  | HVG ARCHS4 pretrain + fine-tune | final | 7 | `R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES` | 0.353 |
| thymus |  | ARCHS4 pretrain + OSDR fine-tune | final | 730 | `R-MMU-9614085_FOXO_MEDIATED_TRANSCRIPTION` | 0.820 |
| thymus |  | direct OSDR | final | 114 | `R-MMU-174143_APC_C_MEDIATED_DEGRADATION_OF_CELL_CYCLE_PROTEINS` | -0.821 |
| thymus |  | HVG ARCHS4 pretrain + fine-tune | final | 146 | `R-MMU-3858494_BETA_CATENIN_INDEPENDENT_WNT_SIGNALING` | -0.831 |

Exploratory random-effects FDR hits that did not pass strict LOO stability:

- skin: 420 random-effects hits in the full ARCHS4-pretrained final model and
  41 in the HVG-pretrained final model, but 0 strict LOO-stable hits.
- kidney: 14 random-effects hits in the HVG-pretrained final model, but 0
  strict LOO-stable hits.
- EDL: 65 random-effects hits in the ARCHS4-pretrained final model and 25 in
  the direct model, but 0 strict LOO-stable hits.
- gastrocnemius: 54 random-effects hits in the ARCHS4-pretrained final model,
  but 0 strict LOO-stable hits.
- tibialis_anterior: 8 random-effects hits in the ARCHS4-pretrained final
  model, but 0 strict LOO-stable hits.

Lung and retina did not show OntoVAE FLT/GC signal by ordinary FDR,
random-effects FDR, or strict LOO stability in the completed runs.

## ARCHS4 Pretraining

ARCHS4 pretraining generally helped by exposing FLT/GC pathway shifts that were
weak or absent in direct OSDR-only runs:

- aggregate skeletal muscle: direct had 0 random-effects hits; ARCHS4-pretrained
  final had 6 random-effects hits and 1 strict LOO-stable hit.
- soleus: direct had 0 random-effects hits; ARCHS4-pretrained final had 60
  random-effects hits and 27 strict LOO-stable hits.
- liver: ARCHS4-pretrained final had 82 random-effects hits and 13 strict
  LOO-stable hits versus 10 and 2 for direct OSDR.

The exceptions are broad immune tissues where direct OSDR already had strong
signals, especially thymus and spleen. In those cases pretraining changed the
ranking and often increased stable-hit counts, but it was not required to see a
condition effect.

## Comparison With expiMap

`outputs/ontovae_pipeline/summary/ontovae_vs_expimap_random_effects.tsv`
contains the side-by-side random-effects/LOO counts against available expiMap
runs.

At a high level:

- OntoVAE recovered strict stable signal in thymus, spleen, liver, aggregate
  skeletal muscle, soleus, and quadriceps.
- expiMap had strong stable thymus and spleen results and strong liver results,
  but aggregate skeletal muscle remained weak.
- For muscle split groups, both methods suggest that splitting muscle is
  important. OntoVAE's strongest stable muscle result is soleus; expiMap also
  showed its strongest split-level stable result in soleus.
- OntoVAE did not clearly outperform expiMap for broad immune tissues. Its
  value here is as a parallel pathway-scoring model with different assumptions,
  not as a drop-in replacement.

## Limitations

- Reactome is used as a flat pathway mask, not as OntoVAE's original GO/HPO DAG
  hierarchy.
- There is no native OntoVAE query-mapping API. Fine-tuning is not equivalent
  to scArches surgery/query mapping.
- There are no de novo OntoVAE programs in this implementation. Reported
  top-gene tables are decoder weights inside existing pathway masks, not new
  genes outside Reactome modules.
- The model uses MSE on transformed expression rather than a count likelihood.
- Broad random-effects signals, especially in thymus, spleen, and skin, should
  still be inspected with PCA/UMAP, accession-colored plots, heatmaps, and
  LOO validation before biological interpretation.
