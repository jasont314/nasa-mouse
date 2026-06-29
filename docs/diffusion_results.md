# Diffusion Results

Generated from `outputs/diffusion_conditional_generation/` on 2026-06-29.

The production run is a pan-tissue conditional DDPM/DDIM model, not eight
independent tissue models. Tissue, material type, muscle group, accession, sex,
platform, source, and condition are conditioning covariates. The run covers all
requested OSDR tissues: liver, skeletal muscle, skin, kidney, thymus, spleen,
lung, and retina. Synthetic examples were generated separately for each tissue
and for skeletal-muscle splits: soleus, gastrocnemius, quadriceps, EDL, and
tibialis anterior.

Machine-readable summaries:

- `outputs/diffusion_conditional_generation/summary/diffusion_training_summary.tsv`
- `outputs/diffusion_conditional_generation/summary/diffusion_analysis_summary.tsv`
- `outputs/diffusion_conditional_generation/summary/diffusion_synthetic_examples_summary.tsv`
- `outputs/diffusion_conditional_generation/summary/diffusion_subgroup_analysis_summary.tsv`
- `outputs/diffusion_conditional_generation/summary/diffusion_reverse_validation_refresh.tsv`

## Completed Tracks

| track | query samples | ARCHS4 ref | genes | landmarks | device | notes |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| OSDR-only | 1,088 | 0 | 9,321 | 512 | A100 CUDA | direct conditional model |
| ARCHS4 pretrain + OSDR fine-tune | 1,088 | 24,428 | 9,319 | 512 | A100 CUDA | reference pretrain, corrected frozen projection, then fine-tune |
| ARCHS4-only | 24,428 | 0 | 9,319 | 512 | A100 CUDA | reference generator/control, no FLT/GC query labels |

Landmarks came from the mouse ortholog mapping of human L1000 genes:
`data/diffusion/l1000_human_to_mouse_ensembl.tsv`.

## FLT vs GC Feature Signals

These are diffusion denoiser-feature signals, not named pathways. Because the
model is conditional, these should be treated as generation/representation
diagnostics, not independent biological modules.

| score set | ordinary Welch FDR < 0.05 | random-effects FDR < 0.05 | LOO-stable FDR < 0.05 |
| --- | ---: | ---: | ---: |
| OSDR-only post-training | 277 | 300 | 271 |
| ARCHS4-pretrained, OSDR fine-tuned | 168 | 130 | 112 |
| Frozen ARCHS4 projection, corrected reference covariates | 27 | 67 | 52 |

Important correction: the frozen projection now replaces OSDR-only
condition/accession/source/platform/material/sex covariates with trained ARCHS4
reference defaults by tissue. The previous uncorrected frozen projection was
not interpretable because it used embeddings never seen during ARCHS4 pretraining.

## Tissue And Muscle-Split Signals

Per-tissue and split-muscle analyses are under:

- `outputs/diffusion_conditional_generation/<track>/analysis_by_subgroup/<score_set>/<subgroup>/`
- `outputs/diffusion_conditional_generation/summary/diffusion_subgroup_analysis_summary.tsv`

All 39 subgroup analyses completed. Rows below show LOO-stable feature counts;
again, these are learned-feature diagnostics rather than named pathway calls.

| subgroup | OSDR-only | ARCHS4 + fine-tune | frozen ARCHS4 projection |
| --- | ---: | ---: | ---: |
| liver | 95 | 39 | 10 |
| skeletal_muscle | 89 | 41 | 49 |
| skin | 85 | 58 | 2 |
| kidney | 53 | 30 | 1 |
| thymus | 41 | 42 | 50 |
| spleen | 49 | 41 | 2 |
| lung | 11 | 18 | 0 |
| retina | 28 | 13 | 0 |
| soleus | 21 | 35 | 42 |
| gastrocnemius | 2 | 3 | 3 |
| quadriceps | 11 | 6 | 3 |
| EDL | 0 | 0 | 0 |
| tibialis anterior | 0 | 0 | 0 |

The most relevant muscle result is that soleus retains many LOO-stable learned
features in the ARCHS4-fine-tuned and frozen-reference score sets. That is
consistent with prior postural-muscle expectations, but diffusion still does
not identify named mitochondrial/contractile/calcium pathways without a separate
attribution step.

## Generated Quality

Generated-expression quality is currently weak. ARCHS4 pretraining helped some
global metrics, but not enough to call the synthetic samples high fidelity.

| track | gene mean corr | gene std corr | Frechet PCA | adversarial accuracy | synthetic-train real-test acc | real-train synthetic-test acc |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OSDR-only | -0.017 | 0.493 | 1,135,921 | 0.800 | 0.492 | 0.498 |
| ARCHS4 + fine-tune | 0.032 | 0.558 | 199,792 | 0.994 | 0.520 | 0.508 |
| Frozen ARCHS4 projection | 0.135 | 0.456 | 138,038 | 0.985 | not run | not run |
| ARCHS4-only | 0.072 | 0.661 | 237,549 | 0.972 | not applicable | not applicable |

The synthetic examples had substantial clipping during inverse transformation
to valid CPM/log1p-CPM space. OSDR-only examples clipped about 91% of exported
values; ARCHS4 + fine-tune examples clipped about 85%. That means the current
samples are useful as pipeline artifacts and counterfactual diagnostics, but
not as high-confidence replacement biological samples.

Synthetic example root:

- `outputs/diffusion_conditional_generation/synthetic_examples/`

All 26 requested model/profile examples completed for OSDR-only and ARCHS4
fine-tuned models. Each directory contains scaled, log1p-CPM, CPM, profile,
clip-report, and mean FLT-minus-GC delta files.

## Comparison To Existing Methods

Diffusion is currently complementary, not stronger.

- Versus WGAN: diffusion has more post-fine-tune LOO-stable pan-tissue feature
  hits than the WGAN conditional post-fine-tune model, but WGAN has better
  established result documentation and stronger split-muscle evidence. Both
  methods produce unnamed learned features requiring attribution.
- Versus OntoVAE: OntoVAE remains more biology-facing because it reports
  pathway/program scores and decoder gene associations. OntoVAE also recovered
  stronger thymus/spleen/liver and soleus signals in prior summaries.
- Versus expiMap: expiMap remains the most interpretable Reactome-module method.
  Diffusion does not yet recover named skeletal-muscle pathways or de novo
  gene modules better than expiMap/OntoVAE.

For prior-literature alignment, diffusion currently supports that FLT/GC states
are separable in a conditional generator, including muscle split profiles, but
it does not by itself recover interpretable mitochondrial, contractile, calcium,
immune, or ECM pathways. Use GLARE/OntoVAE/expiMap/DGEA for named biology and
use diffusion mainly for synthetic-generation experiments after improving
calibration.

## Outputs And Plots

Main analysis plots:

- `outputs/diffusion_conditional_generation/osdr_only/analysis/diffusion_feature_pca.png`
- `outputs/diffusion_conditional_generation/osdr_only/analysis/diffusion_feature_umap.png`
- `outputs/diffusion_conditional_generation/osdr_only/analysis/top_diffusion_feature_shift_heatmap.png`
- `outputs/diffusion_conditional_generation/archs4_pretrain_osdr_finetune/analysis/diffusion_feature_pca.png`
- `outputs/diffusion_conditional_generation/archs4_pretrain_osdr_finetune/analysis/diffusion_feature_umap.png`
- `outputs/diffusion_conditional_generation/archs4_pretrain_osdr_finetune/analysis/top_diffusion_feature_shift_heatmap.png`
- `outputs/diffusion_conditional_generation/archs4_pretrain_osdr_finetune/pretrained_query_analysis/diffusion_feature_pca.png`
- `outputs/diffusion_conditional_generation/archs4_pretrain_osdr_finetune/pretrained_query_analysis/diffusion_feature_umap.png`
- `outputs/diffusion_conditional_generation/archs4_pretrain_osdr_finetune/pretrained_query_analysis/top_diffusion_feature_shift_heatmap.png`

## Limitations

- Learned diffusion features are not pathway modules.
- Condition is an input covariate, so feature-level FLT/GC tests are
  model-diagnostic rather than unsupervised biology.
- The denoiser is much smaller than the paper's largest hidden layers.
- Frechet distance uses PCA rather than the paper's unavailable pretrained
  classifier embedding.
- Reconstruction is ridge LR, not MLP.
- Synthetic expression calibration is poor; clipping reports should be checked
  before using generated CPM matrices downstream.

Best next step: add condition-held-fixed scoring and feature/gene attribution,
then retrain with stronger calibration before treating synthetic samples as
biologically realistic.
