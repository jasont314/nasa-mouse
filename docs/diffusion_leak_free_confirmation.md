# Leakage-Corrected DDIM Confirmation

## Purpose

This work tests whether generated-feature results persist after removing OSDR
overlap from ARCHS4 pretraining. The first phase repeated the frozen lung and
thymus held-out experiment. A subsequent phase retrained the broad factorized
OSDR adapter and regenerated the all-tissue and muscle-group development screens
from the same OSDR-disjoint ARCHS4 backbone. Both phases preserve the original
DDIM architecture and declared feature-policy rules.

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

## Corrected broad OSDR rerun

The broad base adapter completed 12,000 domain and 4,000 condition steps. Its
correlation-refinement continuation completed 4,000 domain and 1,000 condition
steps. Both stages ran on an NVIDIA A100-SXM4-40GB. The corrected locked
within-study test produced:

| Metric | Mean | Repeat pass count |
|---|---:|---:|
| Gene-correlation agreement | 0.9744 | 4/4 against finite-sample floor |
| Precision | 0.9974 | 4/4 |
| Recall | 0.9957 | 4/4 |
| F1 | 0.9966 | 4/4 |
| Adversarial accuracy | 0.4753 | 4/4 |
| FD / real-split P95 | 0.0740 | 4/4 |
| Pooled FLT/GC effect correlation | 0.5408 | 3/4 |
| Muscle accession-effect correlation | 0.5930 | 4/4 |

The finite-sample correlation floor was 0.9497; the stricter paper target was
0.98. A preceding validation calibration missed its sample-specific correlation
floor and muscle accession-effect gate. The all-tissue findings are therefore
developmental, even though the corrected locked output passed its implemented
broad-finalist rule.

Direct augmentation remained unhelpful: real-only balanced accuracy/AUROC was
0.7544/0.8196, compared with 0.6947/0.7514 for synthetic-only and
0.7372/0.7914 for real plus synthetic. The useful mode remained synthetic-guided
feature ranking rather than treating generated profiles as extra animals.

The corrected all-tissue and muscle-group screens changed several prior
attributions. Across 459 real-data BH-FDR tissue-gene associations, 49 met the
corrected synthetic-informed definition: 26 synthetic-promoted and 23
reinforced. Soleus retained five reinforced genes (`Bdh1`, `Ech1`, `Bnip3`,
`Decr1`, and `Tpm1`) but no promoted genes. Kidney promoted `Inpp4b` and
reinforced `Slc37a4`. Spleen promoted `Rai14`, `Ptprk`, and `Myl9` and reinforced
`Loxl1`; `Igfbp3` remained a strong real-data-only association. Quadriceps, EDL,
and liver selected real-only arms and no longer support synthetic-guided claims.

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

The fixed lung/thymus experiment and broader development screens now both use
the OSDR-disjoint reference backbone. Only the lung/thymus experiment excludes
an entire OSDR test accession from adaptation and policy development. The
all-tissue screens interpolate within represented studies and remain
hypothesis-generating. Their random-effects BH-FDR values are based exclusively
on real OSDR profiles.

## Outputs

- Reference run:
  `outputs/generative/benchmark/runs/lacan_diffusion/archs4_mouse_paper_parity_osdr_disjoint_seed1234/`
- OSDR adaptation:
  `outputs/generative/benchmark/runs/lacan_diffusion/osdr_generated_feature_confirmation_disjoint_5000_seed3036/`
- Confirmation:
  `outputs/generative/benchmark/analyses/generated_feature_guidance_confirmation_disjoint_v1/`
- Corrected broad OSDR adaptation:
  `outputs/generative/benchmark/runs/lacan_diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint_seed2020/`
- Corrected broad tissue screen:
  `outputs/generative/benchmark/analyses/within_study_generated_feature_stability_osdr_disjoint_v1/`
- Corrected muscle-group screen:
  `outputs/generative/benchmark/analyses/within_study_generated_feature_stability_muscle_groups_osdr_disjoint_v1/`
- Reference config:
  `configs/generative/diffusion/archs4_mouse_paper_parity_osdr_disjoint.yaml`
- OSDR config:
  `configs/generative/diffusion/osdr_generated_feature_confirmation_disjoint_5000.yaml`
- Feature-policy config:
  `configs/generative/diffusion/generated_feature_guidance_confirmation_disjoint.yaml`
