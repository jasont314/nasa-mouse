# GLARE Validation Stack: 15-Term Candidate Screen

This note records the expanded candidate-module validation run using 15 candidate Reactome terms per tissue per module class. It was run after updating the validation code to retain liver olfactory/chemosensory terms as high-caution candidates.

## Command

```bash
PYTHONUNBUFFERED=1 PYTHONPATH=src /opt/anaconda3/envs/nasa/bin/python -m nasa_mouse_glare.multi_tissue_validation \
  --include-per-study \
  --include-mober \
  --shap-aggregate \
  --verification-estimators 80 \
  --max-eval-genes 4000 \
  --max-eac-genes 1200 \
  --terms-per-class 15 \
  --random-sets 100 \
  --output-dir outputs/glare_multi_tissue_api/validation_stack_terms15
```

## Output Location

Generated outputs are under:

```text
outputs/glare_multi_tissue_api/validation_stack_terms15/
```

The `outputs/` directory is ignored by git; this document tracks the run and interpretation.

## Scope

The run completed successfully for the same 89 scopes:

- 12 aggregate GLARE scopes
- 12 aggregate plus MOBER scopes
- 65 per-study GLARE scopes

The expanded run tested 360 candidate modules:

- 180 DGEA-intersection modules
- 180 GLARE-only modules

The earlier validation run tested 144 modules total, using 6 terms per class.

## Overall Results

| Module class | FDR < 0.05 | Strict: FDR < 0.05 and empirical p <= 0.05 | Median empirical p among all tested | Median empirical p among FDR-significant |
| --- | ---: | ---: | ---: | ---: |
| DGEA-intersection | 76 / 180 | 59 / 180 | 0.0100 | 0.0000 |
| GLARE-only | 34 / 180 | 13 / 180 | 0.2525 | 0.0700 |

Interpretation:

- Expanding from 6 to 15 terms increased the raw GLARE-only FDR-significant count from 17 to 34.
- The broader GLARE-only set is weaker against random gene-set controls: median empirical p among all GLARE-only tested modules worsened from 0.1575 to 0.2525.
- The strongest evidence still comes from DGEA-intersection modules.
- A conservative hidden-module shortlist should use both FDR and empirical random-set support.

## GLARE-Only Results By Tissue

| Tissue | FDR-significant GLARE-only modules | Strict GLARE-only modules |
| --- | ---: | ---: |
| kidney | 3 / 15 | 2 / 15 |
| liver | 1 / 15 | 1 / 15 |
| lung | 0 / 15 | 0 / 15 |
| skeletal_muscle | 10 / 15 | 4 / 15 |
| skeletal_muscle_edl | 0 / 15 | 0 / 15 |
| skeletal_muscle_gastrocnemius | 1 / 15 | 0 / 15 |
| skeletal_muscle_quadriceps | 1 / 15 | 0 / 15 |
| skeletal_muscle_soleus | 5 / 15 | 3 / 15 |
| skeletal_muscle_tibialis_anterior | 2 / 15 | 1 / 15 |
| skin | 0 / 15 | 0 / 15 |
| spleen | 0 / 15 | 0 / 15 |
| thymus | 11 / 15 | 2 / 15 |

Strict means `combined_welch_fdr_bh < 0.05` and `median_empirical_abs_p <= 0.05`.

## Conservative GLARE-Only Shortlist

These are the GLARE-only modules that passed both FDR and empirical random-set support:

| Tissue | Module | FDR | Empirical p | Direction consistency | Mean FLT-GC |
| --- | --- | ---: | ---: | ---: | ---: |
| thymus | Cell Death Signalling Via NRAGE, NRIF and NADE | 7.47e-10 | 0.05 | 0.60 | 0.043 |
| thymus | Response To Elevated Platelet Cytosolic Ca2 | 8.40e-07 | 0.00 | 0.80 | 0.264 |
| skeletal_muscle | Respiratory Electron Transport / ATP Synthesis | 2.12e-06 | 0.01 | 0.62 | -0.101 |
| skeletal_muscle | BMAL1:CLOCK,NPAS2 Activates Circadian Expression | 2.87e-06 | 0.01 | 0.85 | -0.079 |
| skeletal_muscle_soleus | Respiratory Electron Transport / ATP Synthesis | 1.06e-05 | 0.00 | 1.00 | -0.330 |
| skeletal_muscle | Cyclin E Associated Events During G1/S Transition | 2.77e-05 | 0.00 | 0.69 | 0.039 |
| skeletal_muscle_soleus | Cyclin E Associated Events During G1/S Transition | 7.16e-04 | 0.00 | 0.67 | -0.004 |
| skeletal_muscle | L1CAM Interactions | 1.04e-03 | 0.02 | 0.54 | 0.008 |
| kidney | Membrane Trafficking | 1.80e-02 | 0.025 | 0.67 | 0.034 |
| skeletal_muscle_tibialis_anterior | Cyclin E Associated Events During G1/S Transition | 1.93e-02 | 0.00 | 1.00 | 0.113 |
| kidney | Signaling By Insulin Receptor | 2.04e-02 | 0.02 | 0.67 | 0.028 |
| skeletal_muscle_soleus | Signalling By NGF | 2.04e-02 | 0.04 | 0.67 | -0.078 |
| liver | Cleavage Of Growing Transcript In The Termination Region | 2.46e-02 | 0.025 | 0.50 | -0.012 |

Best-supported hidden-module candidates are therefore concentrated in skeletal muscle/soleus, thymus, and kidney. Liver has only one strict GLARE-only module, and its effect is small with only 50% direction consistency.

## Liver GLARE-Only Modules

The expanded liver GLARE-only screen produced one FDR-significant module:

| Module | FDR | Empirical p | Direction consistency | Mean FLT-GC |
| --- | ---: | ---: | ---: | ---: |
| Cleavage Of Growing Transcript In The Termination Region | 0.0246 | 0.025 | 0.50 | -0.012 |

Liver biologically interesting modules did not pass FDR:

| Module | FDR | Empirical p | Direction consistency |
| --- | ---: | ---: | ---: |
| Phase I Functionalization Of Compounds | 0.169 | 0.030 | 0.42 |
| Synthesis Of Bile Acids And Bile Salts | 0.215 | 0.065 | 0.50 |
| Respiratory Electron Transport / ATP Synthesis | 0.246 | 0.020 | 0.42 |
| Olfactory Signaling Pathway | 0.701 | 1.000 | 0.25 |

The liver olfactory term is now included in the candidate table and is clearly not supported by module-score validation in this run.

## Comparison To Original GLARE Cluster Scale

Original GLARE's released OSD-120 cluster outputs contain:

- FLT: 16 consensus clusters
- GC: 15 consensus clusters

The NASA mouse aggregate runs use the same final-cluster scale:

- FLT: 16 clusters per aggregate tissue
- GC: 15 clusters per aggregate tissue

The `--terms-per-class` setting is not the number of GLARE clusters. It controls how many recurring Reactome terms per tissue/class are selected for downstream module-score validation.

## Interpretation

Using 15 terms per class gives a more complete screen and makes the analysis closer to the original GLARE cluster scale, but it does not change the main conclusion:

- DGEA-intersection modules remain much stronger than GLARE-only modules.
- GLARE-only hidden-module support exists, but it is tissue-specific and narrow.
- The strongest hidden candidates are muscle/soleus circadian/respiratory/cell-cycle modules, thymus signaling/cell-death/platelet-calcium modules, and kidney membrane/insulin/fatty-acid modules.
- Liver does not currently have strong hidden-module support. Its strongest biological story remains DGEA-supported metabolic remodeling.
