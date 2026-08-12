# Archived skin and kidney generative-analysis handoff

This file preserves a bounded follow-up analysis. The final synthetic-informed
results and selected runs are documented in the generative manuscript and
`outputs/README.md`.

Last updated: 2026-07-29

## Bottom line

Both tissues contain biologically relevant spaceflight signals, but they should
not be presented at the same evidence level as thymus or soleus.

- **Skin:** the generated-feature analysis recovers a published
  cell-cycle/DNA-repair theme, but it does not provide a robust new
  synthetic-guided gene claim. The strongest real-data candidate is `Plscr1`,
  which is flight-up in all six studies, but it was not selected by the
  generated-feature model and does not survive the strict leave-one-study
  sensitivity rule.
- **Kidney:** the best synthetic-guided candidate is `Slc37a4`, not the
  previously emphasized `Hmox1`/`Alas1` pair. `Slc37a4` is selected repeatedly
  from both real and generated data and is flight-up in all six kidney studies.
  It passes ordinary FDR in both the 974-gene model panel and a
  full-transcriptome TPM sensitivity analysis, but not the strict
  leave-one-study FDR rule. `Inpp4b` is an additional strong real-data result,
  but it was not discovered through synthetic guidance.

The recommended manuscript role is therefore:

1. Kidney as a meaningful secondary, hypothesis-generating result centered on
   renal glucose/metabolic regulation and supported by broader lipid/ECM
   pathway evidence.
2. Skin as a literature-aligned replication of a heterogeneous
   cell-cycle/DNA-repair response, not a new headline result.

## Cohort audit

The API-derived cohort is balanced within every included accession. Flight and
ground-control samples from a given accession have the same recorded material,
sex, and strain.

| Tissue | Profiles | Flight | Ground control | Accessions |
|---|---:|---:|---:|---:|
| Skin | 151 | 80 | 71 | 6 |
| Kidney | 133 | 67 | 66 | 6 |

Skin contains 115 dorsal-skin profiles from four accessions and 36 femoral or
femoral-lateral-skin profiles from two accessions. Kidney contains left, right,
and unspecified kidney samples, but material type is matched between conditions
within each study.

## Skin

### Generated-feature result

The nested development screen selected the `generated_only` arm:

| Metric | Real only | Selected arm | Mean change |
|---|---:|---:|---:|
| Balanced accuracy | 0.671 | 0.757 | +0.086 |
| AUROC | 0.691 | 0.800 | +0.109 |
| Average precision | 0.747 | 0.812 | +0.065 |

These are within-study development results. On the previously reserved
26-profile split, the selected arm was worse than real only:

| Metric | Real only | Generated only |
|---|---:|---:|
| Balanced accuracy | 0.619 | 0.577 |
| AUROC | 0.649 | 0.542 |
| Average precision | 0.738 | 0.601 |

No synthetic-selected skin gene passes the real-data FDR threshold. The selected
genes nevertheless form coherent groups:

- G1/S and DNA-damage control: `Ccne2`, `Ccnd1`, and `Chek1`.
- DNA repair: `Brca1`, `Chek1`, and `Topbp1`.
- Hyaluronan-related signaling: `Cd44` and `Abcc5`.

Reactome over-representation passes FDR for these small selected-gene groups,
but that statistic describes enrichment of a model-selected list. It does not
establish a significant flight effect for the pathway in real samples.

### Real-data candidate

`Plscr1` is the strongest interpretable real-data skin candidate:

- Flight-up in all six accessions.
- Approximately 18% to 34% higher mean TPM in flight in each accession.
- FDR `0.00367` within the prespecified 974-gene panel.
- FDR `0.0335` in an unfiltered 48,303-gene TPM sensitivity analysis.
- Does not pass the strict leave-one-accession-out FDR rule.
- Was not selected as a stable generated-guided feature.

`Plscr1` is interferon inducible and contributes to type-I-interferon signaling,
making it relevant to the immune and damage-response biology reported in
spaceflight skin. It is a real-data follow-up candidate, not evidence that
synthetic generation discovered a stable skin mechanism.

### Relationship to expiMap and prior work

The ARCHS4-reference expiMap result is concentrated in dorsal skin:

- 327 dorsal-skin pathways pass ordinary FDR; 121 have the same direction in
  all four dorsal accessions.
- No dorsal pathway passes the leave-one-accession-out FDR rule.
- Femoral skin has no FDR-significant pathway.
- Direct OSDR expiMap has no FDR-significant pathway in either material group.

The dorsal candidates include mitochondrial translation, mitochondrial
quality control, homologous recombination, and other DNA-repair terms.

This is highly relevant to the published OSDR skin analysis by Cope et al.,
which found 189 cross-mission genes dominated by cell-cycle processes and
reported DNA-repair, mitochondrial, barrier, and collagen/ECM responses.
That paper also found opposing DNA-repair and ECM directions among missions,
strains, recovery times, and anatomical sites. The heterogeneity in the current
analysis is therefore biologically expected rather than evidence that skin is
unaffected.

Recommended interpretation:

> Skin reproduces a known, study-dependent cell-cycle and DNA-damage response.
> It supports the biological relevance of the workflow, but it does not provide
> a robust new synthetic-guided skin discovery.

## Kidney

### Generated-feature result

The nested development screen selected `guided_low_weight`:

| Metric | Real only | Selected arm | Mean change |
|---|---:|---:|---:|
| Balanced accuracy | 0.654 | 0.683 | +0.029 |
| AUROC | 0.680 | 0.773 | +0.093 |
| Average precision | 0.683 | 0.781 | +0.097 |

The same direction of improvement appears on the previously reserved
22-profile split, although that split contains represented studies rather than
an unseen accession.

The primary synthetic-guided gene is `Slc37a4`:

- Selected in all eight real-only and all eight generated-guided repeats.
- Positive classifier coefficient in both arms.
- Flight-up in all six kidney accessions.
- Mean TPM is approximately 3% to 32% higher in flight across the six studies.
- FDR `0.00118` within the 974-gene panel.
- FDR `0.00769` in an unfiltered 48,303-gene TPM sensitivity analysis.
- Does not pass the strict leave-one-accession-out FDR rule.

`Slc37a4` encodes the glucose-6-phosphate translocase used in renal and hepatic
glucose production. Its consistent increase supports a renal metabolic
adaptation hypothesis and is complementary to prior reports of altered renal
lipid and energy metabolism after spaceflight.

Two additional ordinary-FDR genes were found in the real data:

- `Inpp4b`: FDR `0.000556` in the full-transcriptome sensitivity analysis,
  flight-up in five of six studies. It passes the strict rule when correction
  is limited to the 974-gene panel, but not when the full transcriptome is
  reconsidered after each study omission. It was not a stable
  generated-guided feature.
- `Mamld1`: FDR `0.00476` in the full-transcriptome sensitivity analysis,
  flight-up in five of six studies. Its renal interpretation is less clear and
  it was not generated-supported.

The prior porphyrin interpretation should be downgraded. `Hmox1` and `Alas1`
produce a selected-list Reactome enrichment at FDR `0.0469`, but neither gene
has a significant real cross-study effect. `Hmox1` reverses strongly in
OSD-513 and is approximately unchanged in OSD-771. This is an exploratory
oxidative-stress hypothesis, not the main kidney result.

### Relationship to expiMap and prior work

The corrected kidney expiMap reassessment supplies broader pathway context:

- Higher ECM proteoglycan, IGF/IGFBP transport, Wnt, and ECM-degradation
  program scores.
- Lower fatty-acid-metabolism scores.
- Directional support across repeated training seeds and most held-out
  projects.

These are triangulated pathway rankings rather than strict pathway-level FDR
discoveries. They closely match prior kidney studies:

- Finch et al. reported strain-dependent lipid, cholesterol, ECM, TGF-beta,
  Wnt, and inflammatory responses in RR-1 and RR-3 kidney.
- Siew et al. reported renal transporter changes, nephron remodeling, and
  convergent lipid/ECM and oxidative-stress biology across spaceflight
  datasets.
- Suzuki et al. reported altered renal lipid metabolism, blood-pressure, and
  bone-mineralization programs after flight.

The generated core also contains `Egr1`, which Finch et al. independently
reported as flight-up in BALB/c kidney, although `Egr1` is not significant in
the current six-study meta-analysis.

Recommended interpretation:

> Kidney provides a credible synthetic-guided `Slc37a4` metabolic hypothesis
> and independent real-data `Inpp4b` evidence. The broader lipid/ECM/IGF/Wnt
> story is supported by expiMap and prior literature, but an entirely unseen
> kidney study is needed before presenting `Slc37a4` as a confirmed discovery.

## Evidence ranking

| Tissue | Relevance | Statistical interpretation | Recommended paper role |
|---|---|---|---|
| Skin | High for known cell-cycle/DNA-repair biology | Ordinary real-data signal, no strict or transfer confirmation | Literature-aligned supporting result |
| Kidney | High for renal metabolic and remodeling biology | `Slc37a4` passes ordinary FDR in all-study analyses; no strict full-transcriptome LOO or unseen-study confirmation | Stronger secondary hypothesis |

Neither tissue should replace thymus or soleus as the main generative-model
story. Kidney should be promoted above the porphyrin-only wording, while skin
should be retained as evidence that biological heterogeneity across mission,
strain, recovery interval, and anatomical site matters.

## Repository evidence

- `outputs/generative/benchmark/analyses/within_study_generated_feature_stability_v1/`
- `outputs/expimap/analyses/skin_material_split/README.md`
- `outputs/expimap/analyses/kidney_spleen_reassessment/README.md`
- `docs/glare_literature_crosscheck.md`
- `docs/expimap_literature_comparison.md`

The full-transcriptome sensitivity calculation used:

- `outputs/generative/benchmark/data/osdr/osdr_api_raw_counts.h5ad`
- `data/reference/gencode_vM39_mouse_gene_lengths.tsv`

TPM was calculated from full-transcriptome counts and GENCODE gene lengths.
Flight-minus-ground effects were estimated inside each accession and combined
with the same random-effects implementation used by the generated-feature
workflow. The 48,303-gene calculation was an unfiltered sensitivity analysis,
so the 974-gene prespecified-panel and full-transcriptome FDR results are both
reported rather than treating either as uniquely definitive.

## Primary literature

- Cope et al. 2024, *Transcriptomics analysis reveals molecular alterations
  underpinning spaceflight dermatology*:
  https://doi.org/10.1038/s43856-024-00532-9
- Finch et al. 2025, *Spaceflight causes strain-dependent gene expression
  changes in the kidneys of mice*:
  https://doi.org/10.1038/s41526-025-00465-0
- Siew et al. 2024, *Cosmic kidney disease: an integrated pan-omic,
  physiological and morphological study into spaceflight-induced renal
  dysfunction*:
  https://doi.org/10.1038/s41467-024-49212-1
- Suzuki et al. 2022, *Gene expression changes related to bone
  mineralization, blood pressure and lipid metabolism in mouse kidneys after
  space travel*:
  https://doi.org/10.1016/j.kint.2021.09.031
- Talukder et al. 2012, *Phospholipid Scramblase 1 regulates Toll-like
  receptor 9-mediated type I interferon production in plasmacytoid dendritic
  cells*:
  https://doi.org/10.1038/cr.2012.45
