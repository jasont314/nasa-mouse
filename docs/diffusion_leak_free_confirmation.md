# Leakage-Corrected DDIM Confirmation

## Purpose

This rerun tests whether the lung and thymus generated-feature results persist
after removing OSDR overlap from ARCHS4 pretraining. It preserves the original
DDIM architecture, 15,000-epoch reference duration, 5,000-step OSDR adaptation,
contrastive generation seeds, and validation-selected feature-policy grid.

This is a leakage-corrected retest, not a new prospective confirmation, because
the OSD-457 and OSD-900 outcomes had already been inspected before this rerun.

## Leakage controls

The ARCHS4 reference excludes every GEO series linked from the API-derived OSDR
metadata:

```text
data/diffusion/osdr_archs4_overlap_exclusions.tsv
```

Nine series and 108 eligible ARCHS4 profiles were removed, including OSD-457's
`GSE152382`. No excluded series remains in the prepared reference. ARCHS4 roles
were assigned by complete GEO series and have zero train/validation/test series
overlap.

| Reference role | Profiles |
|---|---:|
| Train | 10,150 |
| Validation | 2,466 |
| Test | 4,628 |

OSDR adaptation used 92 profiles from OSD-248, OSD-289, OSD-421, and OSD-515.
Validation used OSD-244 and OSD-464. All 24 OSD-457 and 20 OSD-900 profiles
occurred only in the test role.

## Reference training

The 227,109,786-parameter ARCHS4 DDIM completed 15,000 epochs and 75,000
optimizer steps on an NVIDIA A100-SXM4-40GB in 6,083 seconds.

| Metric | Original reference | OSDR-disjoint reference |
|---|---:|---:|
| Final loss | 2.535 | 2.586 |
| Final noise MAE | 0.02255 | 0.02275 |
| Real-trained tissue BA | 0.895 | 0.781 |
| Synthetic-trained tissue BA | 0.869 | 0.781 |
| Real/synthetic adversarial accuracy | 0.512 | 0.515 |
| Precision | 0.966 | 0.951 |
| Recall | 0.865 | 0.890 |
| Gene-correlation-matrix agreement | 0.879 | 0.878 |
| Memorization fraction below train LOO P01 | 0.0035 | 0.0035 |

The lower real-trained tissue BA reflects the harder complete-series test split.
Synthetic and real training produced essentially identical tissue BA on that
split. The strict correlation-matrix gate still fails (`0.878` versus `0.952`
required), but the original model also scored `0.879`; this is an existing model
limitation rather than an effect of leakage removal.

## Confirmation results

| Tissue | Test accession | Baseline BA | Guided BA | Baseline AUROC | Guided AUROC | Baseline AP | Guided AP | Result |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Lung | OSD-900 | 0.400 | 0.350 | 0.450 | 0.470 | 0.523 | 0.578 | Rejected by validation gate |
| Thymus | OSD-457 | 0.500 | 0.833 | 0.840 | 0.979 | 0.876 | 0.983 | Retained |

The lung candidate lost balanced accuracy. Its validation gate did not pass, so
the deployed policy correctly remained the real-only baseline.

The thymus classifier used real training profiles only; generated profiles
changed feature ranking. Performance improved in both genotype strata:

| Genotype | Baseline BA | Guided BA | Baseline AUROC | Guided AUROC |
|---|---:|---:|---:|---:|
| Nrf2KO | 0.500 | 0.917 | 0.944 | 1.000 |
| WT | 0.500 | 0.750 | 0.778 | 0.972 |

All eight previously reported thymus genes, `Birc5`, `Cdk1`, `Ccnb2`, `Nusap1`,
`Ccnb1`, `Gmnn`, `Ccne2`, and `Ube2c`, remained in the selected 100-gene panel.
Their real and generated training effects agreed in direction, and their real
OSD-457 effects remained flight-lower in both genotypes. Eighty-nine Reactome
terms passed pathway FDR 0.05, led by mitotic cell cycle, cell-cycle checkpoints,
APC/C-mediated degradation, and DNA replication.

## Interpretation

The original lung improvement does not survive the stricter retraining and
should not be presented as confirmed. The thymus signal does survive removal of
the exact OSD-457 profiles and all other OSDR-linked GEO series from ARCHS4
pretraining. This substantially strengthens the evidence that synthetic-guided
feature ranking organizes a transferable thymus cell-cycle panel.

It does not restore prospective independence: OSD-457 had already been examined
when this correction was designed. A newly reserved thymus accession or external
dataset remains necessary for a genuinely untouched confirmation.

This rerun covers the fixed lung/thymus confirmation experiment. The broader
all-tissue synthetic-selection screens were not regenerated with this corrected
backbone and remain hypothesis-generating. Their real-only random-effects BH-FDR
results are unaffected by the generator overlap.

## Outputs

- Reference run:
  `outputs/generative_benchmark/runs/lacan_diffusion/archs4_mouse_paper_parity_osdr_disjoint_seed1234/`
- OSDR adaptation:
  `outputs/generative_benchmark/runs/lacan_diffusion/osdr_generated_feature_confirmation_disjoint_5000_seed3036/`
- Confirmation:
  `outputs/generative_benchmark/analyses/generated_feature_guidance_confirmation_disjoint_v1/`
- Reference config:
  `configs/rna_diffusion/archs4_mouse_paper_parity_osdr_disjoint.yaml`
- OSDR config:
  `configs/rna_diffusion/osdr_generated_feature_confirmation_disjoint_5000.yaml`
- Feature-policy config:
  `configs/rna_diffusion/generated_feature_guidance_confirmation_disjoint.yaml`
