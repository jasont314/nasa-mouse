# Tutorial-Style expiMap Liver Run

This run follows the scArches expiMap tutorial workflow more closely than the
full Reactome liver runs.

Reference workflow:

- https://docs.scarches.org/en/latest/expimap_surgery_pipeline_basic.html
- https://docs.scarches.org/en/stable/expimap_surgery_pipeline_advanced.html

## What Changed

- HVG selection was applied to the ARCHS4 reference after normalization/log1p,
  using `scanpy.pp.highly_variable_genes(n_top_genes=2000,
  batch_key="archs4_condition")`.
- Reactome terms were re-filtered after HVG selection with the tutorial-style
  strict `>12` genes rule.
- The reference model used `hidden_layer_sizes=[300, 300, 300]`, raw-count NB
  loss, `alpha=0.7`, `alpha_kl=0.5`, `alpha_epoch_anneal=100`, and a 400-epoch
  early-stopping schedule.
- Query surgery used three unconstrained de novo nodes, `gamma_ext=0.7`,
  `gamma_epoch_anneal=50`, `alpha_kl=0.22`, 250 epochs, and HSIC
  one-vs-all regularization with `beta=3`.

The advanced tutorial also adds one constrained B-cell extension node from a
withheld B-cell receptor mask. There is no direct NASA liver analogue for that
withheld mask, so this run uses only unconstrained de novo extension nodes.

## Inputs And Training

Input artifacts:

- `outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_5000/input/`
- `archs4_mouse_liver_reference_tutorial_hvg_raw_counts.h5ad`
- `osdr_liver_query_tutorial_hvg_raw_counts.h5ad`

HVG filtering retained 1,995 genes and 364 Reactome terms from the 5,000-sample
ARCHS4 liver reference and the 231-sample OSDR liver query.

Reference output:

- `outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_5000/reference_nb_400epoch_seed2020/`

The reference run early-stopped after 248 epochs, with best epoch 196.

Query output:

- `outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_5000/query_denovo3_hsic_250epoch_seed2020/`

The query run completed all 250 epochs with 367 score dimensions: 364 Reactome
terms plus `unconstrained_0`, `unconstrained_1`, and `unconstrained_2`.

## HSIC Patch

The installed scArches HSIC implementation computes a gamma-function ratio
directly. At this latent dimensionality, that ratio overflows and produces NaN
values. The mapper now applies a runtime patch when HSIC is enabled, replacing
the direct gamma ratio with the equivalent stable log-gamma form and clamping
the final HSIC quantity before square root.

The query summary records this as `"stable_hsic_patch": true`.

## Results

The pooled FLT-vs-GC Welch test found several FDR-significant pathway score
differences. The top pooled terms include:

- `R-MMU-141444_AMPLIFICATION_OF_SIGNAL_FROM_UNATTACHED_KINETOCHORES_VIA_A_MAD2_INHIBITORY_SIGNAL`
- `R-MMU-8848021_SIGNALING_BY_PTK6`
- `R-MMU-70326_GLUCOSE_METABOLISM`
- `R-MMU-422356_REGULATION_OF_INSULIN_SECRETION`
- `R-MMU-9013148_CDC42_GTPASE_CYCLE`

Accession-aware random-effects meta-analysis found 70 Reactome terms with
meta-analysis FDR < 0.05. Thirty-two terms also had every leave-one-accession-out
run below FDR 0.05 and in the same meta-analysis direction. The strongest such
terms include:

- `R-MMU-212300_PRC2_METHYLATES_HISTONES_AND_DNA`
- `R-MMU-8963899_PLASMA_LIPOPROTEIN_REMODELING`
- `R-MMU-1483257_PHOSPHOLIPID_METABOLISM`
- `R-MMU-1630316_GLYCOSAMINOGLYCAN_METABOLISM`
- `R-MMU-194315_SIGNALING_BY_RHO_GTPASES`

The study-aware Wilcoxon test across accession effects did not reach FDR 0.05;
its minimum FDR was 0.114.

The expiMap `latent_enrich` condition test ran with `ground_control` as the
comparison group. The largest absolute Bayes-factor-like score was small
(`abs(bf) = 0.442` for `R-MMU-9012999_RHO_GTPASE_CYCLE`), so this does not
provide strong condition-level enrichment evidence. Applying the local expiMap
paper-style cutoff of `abs(bf) >= 2.3` gives zero passing terms; the ranked
summary is in
`outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_5000/query_denovo3_hsic_250epoch_seed2020/latent_enrich_condition/paper_style_bf_summary/`.

The three de novo programs were not FLT/GC-significant:

| program | Welch FDR | study-aware FDR | note |
| --- | ---: | ---: | --- |
| `unconstrained_0` | 0.531 | 0.205 | diffuse, protein/peroxisomal correlation |
| `unconstrained_1` | 0.440 | 0.437 | diffuse, platelet-development correlation |
| `unconstrained_2` | 0.448 | 0.940 | diffuse, drug-ADME correlation |

The previous Polymerase II elongation and semaphorin candidates are absent from
this HVG-reduced architecture after term filtering, so this run cannot validate
or refute those exact terms.

## Interpretation

This run now matches the expiMap tutorial mechanics much more closely:
HVG-reduced architecture, deeper reference model, longer training, de novo query
nodes, HSIC one-vs-all regularization, and expiMap latent enrichment.

It still should not be treated as final biology. The strongest robust-looking
signals are model-score candidates from one reference seed. They need seed
replication and count-level accession-aware validation before promotion.
