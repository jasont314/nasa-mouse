# Archived material-type conditioning ablation

Last updated: 2026-08-04

This file preserves a bounded development analysis. Current final-model and
publication references are listed in `outputs/README.md` and
`paper/synthetic_guided_spaceflight/README.md`.

## Purpose

This document hands off the side-conversation analysis comparing the primary
material-conditioned DDIM against a matched DDIM trained without material-type
conditioning. It records the exact model and output locations, aggregate
generator metrics, downstream classifier and gene-selection changes, targeted
literature review of the no-material-only genes, limitations, and the recommended
interpretation.

The immediate question was whether conditioning on `material_type` is better
than conditioning on tissue without material type.

## Bottom line

- The two models are effectively tied on aggregate generator fidelity. Neither
  model clearly wins correlation, precision/recall/F1, adversarial accuracy, and
  Frechet ratio as a group.
- The material-conditioned model is slightly better at preserving FLT-GC effects:
  condition-effect correlation is 0.751 versus 0.743, and condition-direction
  agreement is 0.789 versus 0.777.
- The no-material model is slightly closer to ideal on development-validation
  adversarial accuracy and Frechet ratio, and better on the reported
  muscle-accession metrics. These differences do not amount to a global fidelity
  win.
- Explicit material conditioning is preferable as the primary scientific model
  because material labels encode real anatomical distinctions, especially among
  skeletal-muscle groups. It preserves a coherent soleus mitochondrial/lipid
  panel that disappears without material conditioning.
- The no-material model has useful local results. It recovers strong,
  literature-aligned thymus `Ccna2` and `Ccnb1` signals and pooled-muscle
  `Cdkn1a`.
- The thymus cell-cycle interpretation is robust to the conditioning choice. The
  soleus synthetic-informed interpretation is conditioning-sensitive.
- Use the material-conditioned model as primary and the no-material model as a
  sensitivity analysis. Do not claim that material conditioning is superior on
  aggregate metrics.

## Critical terminology

### Tissue versus material type

`tissue` is the broad organ label, such as thymus, spleen, kidney, or pooled
skeletal muscle. `material_type` records a more specific sampled anatomical
source when available. For skeletal muscle, this distinction is important
because soleus, gastrocnemius, EDL, quadriceps, and tibialis anterior differ in
fiber composition, loading, metabolism, and response to unloading.

### What the no-material model actually conditions on

The ablation is not an unconditional or strictly tissue-only model. It retains:

- tissue
- FLT versus ground-control condition
- study/accession

It removes only:

- material type

This distinction matters because accession/study identifies many material
labels. Study conditioning can therefore act as a proxy for material in this
dataset. The experiment measures the incremental value of explicit material
conditioning after tissue and study are already known. It does not prove that
material information is generally unnecessary.

## Matched training design

Both models use the same ARCHS4 mouse pretrained DDIM checkpoint, prepared OSDR
data, train/validation split, seed, architecture, optimizer settings, and training
schedule. The only intended conditioning difference is `material_type`.

| Property | Material | No material |
|---|---:|---:|
| Tissue conditioning | yes | yes |
| FLT/GC conditioning | yes | yes |
| Study conditioning | yes | yes |
| Material conditioning | yes | no |
| Training profiles | 781 | 781 |
| Validation profiles | 536 | 536 |
| Domain-refinement steps | 4,000 | 4,000 |
| Condition-refinement steps | 1,000 | 1,000 |
| Learning rate | 2e-5 | 2e-5 |
| Batch size | 512 | 512 |
| GPU | NVIDIA A100-SXM4-40GB | NVIDIA A100-SXM4-40GB |
| Base parameters | 227,109,786 | 227,109,786 |
| Adapter parameters | 69,987,328 | 65,563,648 |
| Locked test opened for this run | no | no |

Shared ARCHS4 checkpoint:

`outputs/generative/benchmark/runs/lacan_diffusion/archs4_mouse_paper_parity_osdr_disjoint_seed1234/model.pt`

The resolved metadata records epoch 15,000, global step 75,000, EMA weights, and
checkpoint SHA-256
`24eedb834fd756781e2dc404868f28ff37f4c4aaac954b11edfc307c57fffc3d`.

## Files and outputs

### Material-conditioned model

- Refine config:
  `configs/generative/diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint.yaml`
- Model:
  `outputs/generative/benchmark/runs/lacan_diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint_seed2020/model.pt`
- Canonical-tissue analysis:
  `outputs/generative/benchmark/analyses/within_study_generated_feature_stability_osdr_disjoint_v1/`
- Muscle-group analysis:
  `outputs/generative/benchmark/analyses/within_study_generated_feature_stability_muscle_groups_osdr_disjoint_v1/`

### No-material model

- Initial matched-training config:
  `configs/generative/diffusion/osdr_factorized_study_lora512_replicated_validation_osdr_disjoint_no_material.yaml`
- Refine config:
  `configs/generative/diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint_no_material.yaml`
- Canonical-tissue feature config:
  `configs/generative/diffusion/within_study_generated_feature_stability_no_material_osdr_disjoint.yaml`
- Muscle-group feature config:
  `configs/generative/diffusion/within_study_generated_feature_stability_muscle_groups_no_material_osdr_disjoint.yaml`
- Model run:
  `outputs/generative/benchmark/runs/lacan_diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint_no_material_seed2020/`
  The nonselected checkpoint was pruned after its metrics and analysis tables
  were frozen; rerun the refine config to recreate it.
- Canonical-tissue analysis:
  `outputs/generative/benchmark/analyses/within_study_generated_feature_stability_no_material_osdr_disjoint_v1/`
- Muscle-group analysis:
  `outputs/generative/benchmark/analyses/within_study_generated_feature_stability_muscle_groups_no_material_osdr_disjoint_v1/`

### Direct comparison

- Comparison script:
  `src/nasa_mouse_diffusion/paper_parity/compare_material_conditioning_ablation.py`
- Comparison root:
  `outputs/generative/benchmark/analyses/material_type_conditioning_ablation_v1/`
- Main saved summary:
  `outputs/generative/benchmark/analyses/material_type_conditioning_ablation_v1/README.md`
- Machine-readable summary:
  `outputs/generative/benchmark/analyses/material_type_conditioning_ablation_v1/summary.json`
- Generator metrics:
  `outputs/generative/benchmark/analyses/material_type_conditioning_ablation_v1/calibrated_generator_metrics.tsv`
- Arm comparison:
  `outputs/generative/benchmark/analyses/material_type_conditioning_ablation_v1/arm_choice_comparison.tsv`
- Gene comparison:
  `outputs/generative/benchmark/analyses/material_type_conditioning_ablation_v1/synthetic_informed_gene_comparison.tsv`
- Full no-material BH-FDR inventory:
  `outputs/generative/benchmark/analyses/material_type_conditioning_ablation_v1/no_material_all_bh_fdr_genes.tsv.gz`
- Pathway comparison:
  `outputs/generative/benchmark/analyses/material_type_conditioning_ablation_v1/significant_pathway_comparison.tsv.gz`

To regenerate the comparison from existing outputs:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m nasa_mouse_diffusion.paper_parity.compare_material_conditioning_ablation
```

## Generator validation metrics

These are matched repeated development-validation results.

| Metric | Desired direction | Material | No material | Interpretation |
|---|---|---:|---:|---|
| Correlation | higher | 0.972556 | 0.972529 | material higher by 0.000027; tied in practice |
| Precision | higher | 0.995802 | 0.996269 | no material higher by 0.000466 |
| Recall | higher | 0.995802 | 0.995336 | material higher by 0.000466 |
| F1 | higher | 0.995797 | 0.995796 | effectively identical |
| Adversarial accuracy | closer to 0.5 | 0.586754 | 0.586521 | no material closer by 0.000233; effectively identical |
| Frechet ratio | lower | 0.110703 | 0.108474 | no material lower by 0.002229 |
| FLT-GC delta correlation | higher | 0.751358 | 0.742988 | material better by 0.008369 |
| FLT-GC direction agreement | higher | 0.789014 | 0.776694 | material better by 0.012320 |
| Muscle-accession correlation | higher | 0.143768 | 0.217025 | no material better by 0.073257 |
| Muscle-accession direction | higher | 0.573409 | 0.591376 | no material better by 0.017967 |

Interpretation:

- Global synthetic fidelity is a tie.
- Material conditioning has a small advantage for the FLT-GC contrast, which is
  the principal biological target.
- No-material has a local advantage for the reported muscle-accession metrics.
- The metric differences are mixed. Selecting the material model solely because
  its genes look more biologically familiar would be circular; the scientific
  reason to retain it is that material is a known anatomical design variable and
  the primary model was specified with that variable.

### Adversarial-accuracy caveat

For adversarial accuracy, 0.5 is ideal because a classifier cannot distinguish
real from generated profiles. The approximately 0.587 values above are
development-validation values. The no-material model was not evaluated on a new
confirmatory locked test. Its 0.587 must not replace or be compared as if it were
the primary model's previously reported locked-test AA of 0.475.

## Downstream feature-selection workflow

The existing workflow was rerun without material conditioning. It used the same:

- three synthetic draws
- eight repeated nested development splits
- 27 analysis units: 22 canonical tissues plus five muscle groups
- five synthetic-use arms
- balanced accuracy, AUROC, and average precision eligibility rule
- stable-feature thresholds
- real-data accession-specific effects and random-effects BH FDR
- Reactome analysis

The five arms are:

1. `real_only`: classifier and feature ranking from observed profiles only.
2. `generated_only`: classifier and ranking from generated profiles.
3. `real_plus_generated`: equal-use real and generated training.
4. `guided_real_only`: generated data influence feature ranking, but the final
   classifier is fit only on real profiles.
5. `guided_low_weight`: generated data influence ranking and enter training at
   0.05 total weight after recentering.

Stable genes require selection frequency at least 0.50 and coefficient-sign
agreement at least 0.75.

- `synthetic_promoted`: stable only with an eligible synthetic-informed arm and
  supported by a real-data meta-effect.
- `reinforced_real_and_synthetic`: stable in real-only and synthetic-informed
  selection with matching direction and a real-data meta-effect.

These labels describe selection behavior. They do not mean that generated
profiles constitute biological replicates, and `promoted` does not mean novel.
All reported association directions and BH-FDR values come from observed OSDR
profiles, not synthetic profiles.

## Downstream summary

| Result | Value |
|---|---:|
| Analysis units | 27 |
| Selected arm changed | 13 |
| Material synthetic-informed associations | 49 |
| No-material synthetic-informed associations | 45 |
| Shared with same promoted/reinforced status | 32 |
| Material-only | 17 |
| No-material-only | 13 |
| Union | 62 |
| Jaccard overlap | 0.516 |
| Material promoted/reinforced | 26 / 23 |
| No-material promoted/reinforced | 26 / 19 |
| Real-data statistics changed | no |

The 32 retained genes are 65.3% of the 49 material-conditioned associations.
Thresholded biological attribution is therefore only moderately robust even
though aggregate generator metrics are almost identical.

### Analysis units whose selected arm changed

The balanced-accuracy difference below is no-material selected-arm BA minus
material selected-arm BA. It is not a new held-out-study estimate.

| Unit | Material arm | No-material arm | BA difference |
|---|---|---|---:|
| Bone | guided real only | generated only | +0.016 |
| Cerebellum | generated only | guided low weight | +0.069 |
| Hippocampus | generated only | real only | -0.005 |
| Kidney | guided low weight | guided real only | -0.010 |
| Liver | real only | guided real only | +0.034 |
| Mammary gland | generated only | real plus generated | +0.042 |
| Retina | guided low weight | generated only | +0.034 |
| Skeletal muscle, pooled | guided real only | guided low weight | -0.014 |
| Skin | real plus generated | guided real only | -0.031 |
| Thymus | guided low weight | generated only | +0.027 |
| EDL | real only | guided low weight | +0.042 |
| Gastrocnemius | guided low weight | generated only | +0.052 |
| Soleus | real plus generated | real only | -0.037 |

The arm changes show that downstream choice is more sensitive than aggregate
distribution metrics. Some BA changes also trade against AUROC or average
precision, so BA alone should not determine a winner.

## Gene-level comparison

### Shared 32 associations

These retained the same promoted or reinforced status under both models:

- Adrenal gland: `Tspan4`
- Eye: `Klhl21`
- Kidney: `Slc37a4`
- Skeletal muscle, pooled: `Sox4`, `Cebpd`, `Sh3bp5`, `Bphl`, `Klhl21`,
  `Prkcd`, `Mapkapk5`, `Arid5b`, `Reep5`, `Sesn1`, `Tle1`, `Itgb5`
- Skin: `Plscr1`
- Spleen: `Myl9`, `Loxl1`
- Thymus: `Nusap1`, `Stmn1`, `Birc5`, `Cdk1`, `Top2a`, `Ccnb2`, `Aurka`,
  `Ccne2`, `Gmnn`, `Ccnf`
- Tibialis anterior: `Cdkn1a`, `St3gal5`, `Cebpd`, `Bnip3`

### Material-only 17 associations

- Adrenal gland: `Psmb8`
- Kidney: `Inpp4b`
- Spleen: `Rai14`, `Ptprk`
- Thymus: `Snx7`, `Hsd17b11`, `Ube2c`, `Etv1`, `Kif20a`, `Pcna`
- Gastrocnemius: `Nfkbia`, `Fhl2`
- Soleus: `Bdh1`, `Ech1`, `Bnip3`, `Decr1`, `Tpm1`

### No-material-only 13 associations

| Tissue | Gene | Selection | Real FLT direction | BH FDR | Accessions agreeing | LOO status | Literature class |
|---|---|---|---|---:|---:|---|---|
| Adrenal | `Herpud1` | promoted | lower | 0.00379 | 3/3 | not stable | unmatched |
| Retina | `Inpp1` | promoted | higher | 0.0334 | 4/4 | not stable | complementary |
| Skeletal muscle | `Cdkn1a` | reinforced | higher | 0.000660 | 12/13 | stable | aligning |
| Skeletal muscle | `Dusp3` | promoted | lower | 0.00523 | 11/13 | not stable | unmatched |
| Skeletal muscle | `Gpatch8` | promoted | higher | 0.0377 | 11/13 | not stable | unmatched |
| Skeletal muscle | `Pex11a` | promoted | lower | 0.0467 | 10/13 | not stable | complementary |
| Thymus | `Ccnb1` | promoted | lower | 3.40e-10 | 5/5 | stable | aligning |
| Thymus | `Ccna2` | promoted | lower | 7.95e-7 | 5/5 | stable | aligning |
| Thymus | `Socs2` | promoted | higher | 0.0443 | 4/5 | not stable | complementary |
| EDL | `Abcc5` | reinforced | lower | 1.43e-24 | 2/2 | not assessable robustly | unmatched |
| EDL | `Rsu1` | reinforced | lower | 3.99e-8 | 2/2 | not assessable robustly | complementary |
| EDL | `Zfp131` | promoted | lower | 1.68e-5 | 2/2 | not assessable robustly | unmatched |
| EDL | `Bnip3` | promoted | higher | 0.00113 | 2/2 | not assessable robustly | ambiguous |

For EDL, two accessions are insufficient for a meaningful leave-one-accession-out
stability claim. The stored boolean is false, but the appropriate interpretation
is limited identifiability rather than a conventional LOO failure.

## Targeted review of the no-material-only genes

The 13 no-material-only associations were reviewed against targeted primary
literature. These classifications are captured here but have not yet been added
to the deterministic manuscript annotation table.

### Aligning

- Thymus `Ccnb1`, FLT lower: directly coherent with prior mouse-flight thymus
  suppression of proliferative/cell-cycle genes across missions.
- Thymus `Ccna2`, FLT lower: same interpretation and unusually strong real-data
  support here: 5/5 accessions, BH FDR 7.95e-7, LOO stable.
- Pooled-muscle `Cdkn1a`, FLT higher: consistent with published flight and
  unloading stress/senescence responses. It is supported in 12/13 accessions and
  passes LOO sensitivity.

### Complementary

- Retina `Inpp1`, FLT higher: a plausible inositol-signaling extension to prior
  spaceflight retina signaling changes, but not an exact same-gene replication.
- Pooled-muscle `Pex11a`, FLT lower: mechanistically compatible with reduced
  peroxisomal/oxidative capacity, but no exact pooled-muscle spaceflight match was
  found.
- Thymus `Socs2`, FLT higher: connects cytokine/JAK-STAT regulation and T-cell
  biology to the thymus response, without a direct same-gene thymus-flight match.
- EDL `Rsu1`, FLT lower: plausible focal-adhesion/costamere signaling extension;
  no direct EDL-flight replication and no significant EDL Reactome module.

### Ambiguous

- EDL `Bnip3`, FLT higher: unloading and flight evidence for BNIP3/autophagy is
  sex-, duration-, assay-, and muscle-dependent. The no-material EDL direction
  also contrasts with the material-conditioned soleus `Bnip3` FLT-lower result.

### Unmatched after targeted search

- Adrenal `Herpud1`
- Pooled-muscle `Dusp3`
- Pooled-muscle `Gpatch8`
- EDL `Abcc5`
- EDL `Zfp131`

`Unmatched` means that the targeted search did not find a sufficiently specific
tissue-spaceflight or process-level match. It does not establish novelty and does
not mean the gene is biologically implausible.

### Literature-class totals

The 32 shared genes retain their prior material-model literature labels. Combining
those labels with the targeted review above gives this provisional no-material
breakdown:

| Literature class | Material model, persisted | No-material model, provisional |
|---|---:|---:|
| Aligning | 22 | 18 |
| Complementary | 19 | 17 |
| Ambiguous | 4 | 3 |
| Unmatched | 4 | 7 |
| Total | 49 | 45 |

The material-conditioned result has a lower unmatched fraction: 4/49 versus 7/45.
That supports its interpretability, but this should not be treated as an unbiased
performance metric because literature concordance was assessed after gene
selection.

### Key literature used in the targeted review

- Horie et al., mouse spaceflight thymus and cell-cycle suppression, including
  `Ccna2` and `Ccnb1`: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6934594/>
- Allen et al., mouse flight skeletal muscle, including higher `Cdkn1a` and
  `Socs2` in gastrocnemius:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC2644242/>
- Oommen et al., public mouse spaceflight muscle synthesis and recurrent
  `Cdkn1a`: <https://www.nature.com/articles/s41526-024-00434-z>
- `Pex11a` and peroxisomal lipid/oxidative biology:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6378206/>
- Spaceflight skeletal-muscle oxidative capacity:
  <https://pubmed.ncbi.nlm.nih.gov/12882990/>
- SOCS2 regulation in T-cell biology:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC3135359/>
- STS-135 retina signaling context:
  <https://reference-global.com/download/article/10.2478/gsr-2014-0001.pdf>
- RSU1 and focal-adhesion signaling:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC7909951/>
- BNIP3 response to unloading:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC8718086/>
- Long-duration spaceflight, EDL, and autophagy context:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC3314659/>

The primary material-conditioned 49-gene annotation is already persisted at:

`paper/synthetic_guided_spaceflight/source_data/table_s16_promoted_gene_literature_annotations.tsv`

Its source inventory is:

`paper/synthetic_guided_spaceflight/source_data/table_s17_promoted_gene_literature_sources.tsv`

## Pathway changes

### Thymus

- Material model: 45 significant eligible synthetic-supported Reactome terms.
- No-material model: 71.
- Shared: 29.
- Both preserve a coherent cell-cycle signal, including mitotic cell cycle, DNA
  replication, APC/C control, G2/M transition, and checkpoint regulation.
- No material adds `Ccna2` and `Ccnb1`, which strengthen the same core biological
  interpretation rather than creating a contradictory story.

Conclusion: thymus is robust at the process level and is the strongest result
under either conditioning choice.

### Soleus

- Material model: five significant pathways, including mitochondrial fatty-acid
  beta oxidation, lipid metabolism, and mitochondrial protein degradation.
- No-material model: zero.
- The material model's five reinforced genes, `Bdh1`, `Ech1`, `Bnip3`, `Decr1`,
  and `Tpm1`, lose synthetic attribution because soleus changes to a real-only
  arm.

Conclusion: the soleus biological panel is coherent under explicit material
conditioning but not robust to removal of that label. This is expected to some
degree because soleus is itself an anatomical material distinction.

### Pooled skeletal muscle

- Material model: six significant terms.
- No-material model: four.
- Only beta-catenin complex deactivation is shared.
- The material model carries interferon and sialic-acid terms.
- The no-material model instead carries ERBB2/SHC1 and TCF/WNT-related terms.

Conclusion: pooled-muscle pathway interpretation is conditioning-sensitive even
though 12 gene associations are shared.

### EDL

The no-material model adds `Abcc5`, `Bnip3`, `Rsu1`, and `Zfp131`, but these do not
form a significant Reactome pathway. With only two accessions, this is exploratory
and less coherent than the material-conditioned soleus result.

### Kidney and spleen

- No-material kidney gains four pathways centered on cholesterol, steroid, and
  lipid metabolism.
- No-material spleen gains two collagen-formation pathways.
- These enrichments are conditioning-sensitive and should not be presented as
  robust discoveries without further validation.

## Interpretation for manuscript or presentation

Use language equivalent to:

> Removing material-type conditioning left aggregate DDIM fidelity essentially
> unchanged but altered downstream arm choice and synthetic-informed feature
> attribution. The thymus mitotic and DNA-replication program was preserved and
> expanded to include Ccna2 and Ccnb1. In contrast, the soleus mitochondrial and
> lipid-metabolism panel required explicit material conditioning. We therefore
> retained material conditioning for the primary analysis and treated the
> no-material model as a sensitivity analysis.

Avoid these claims:

- "Material conditioning clearly improves all model metrics."
- "The no-material model proves material is irrelevant."
- "Synthetic samples made the genes statistically significant."
- "No-material-only genes are novel."
- "EDL is confirmed," given only two accessions and no coherent pathway.

## Recommended next actions for the main session

1. Keep the material-conditioned DDIM as the primary generator.
2. Add the no-material result as a compact sensitivity analysis if the manuscript
   discusses conditioning design. The key message is robust thymus and
   conditioning-sensitive soleus.
3. If no-material genes will appear in a paper or slide, persist their 13-row
   literature annotations and source relationships in deterministic source-data
   tables rather than relying only on this handoff.
4. Do not open or report the no-material locked test as though it were a fresh
   confirmatory comparison without first specifying the decision rule. Prefer a
   genuinely independent future accession for model-choice confirmation.
5. For a clean causal assessment of material conditioning, run a factorial
   design with study on/off and material on/off, or evaluate held-out
   study-material combinations. The current ablation only tests material's
   incremental value when study is retained.
6. Do not overwrite the primary 49-gene tables with the no-material 45-gene set.
   Keep the ablation separate and label it as sensitivity analysis.
