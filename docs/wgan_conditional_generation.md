# WGAN Conditional Generation

Generated on 2026-06-29 from API-derived OSDR inputs and existing ARCHS4
mouse reference inputs. Code lives under `src/nasa_mouse_wgan/`; trained
outputs are under `outputs/wgan_conditional_generation/`.

## Conditioning Design

The pan-tissue generator is conditioned on:

- `wgan_condition`: `flight`, `ground_control`, or `archs4_reference`
- `wgan_tissue`
- `wgan_material_type`
- `wgan_muscle_group`
- `wgan_accession`
- `wgan_sex`
- `wgan_assay`
- `wgan_platform`
- `wgan_data_source`

This supports counterfactual-style generation by holding the nuisance/design
covariates fixed and flipping `wgan_condition` between `ground_control` and
`flight`.

## Trained Tracks

Command:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_wgan.run_conditional_generation
```

| track | samples | genes | epochs | device | purpose |
| --- | ---: | ---: | ---: | --- | --- |
| `osdr_only` | 1,088 OSDR | 9,321 | 100 OSDR | A100 CUDA | Direct OSDR conditional generator |
| `archs4_pretrain_osdr_finetune` | 24,428 ARCHS4 + 1,088 OSDR | 9,319 | 100 ARCHS4 + 100 OSDR | A100 CUDA | Reference-pretrained generator fine-tuned to OSDR |
| `archs4_only` | 24,428 ARCHS4 | 9,319 | 100 ARCHS4 | A100 CUDA | ARCHS4 reference-distribution generator/control |

The ARCHS4-only model has only the `archs4_reference` condition. It cannot
generate flight or ground-control samples unless an OSDR fine-tune stage adds
those condition labels.

## Fit And FLT/GC Signal

| score set | gene mean r | gene std r | RE FDR<0.05 | LOO-stable FDR<0.05 |
| --- | ---: | ---: | ---: | ---: |
| OSDR-only final | -0.0025 | 0.4004 | 2 | 0 |
| ARCHS4-pretrained frozen OSDR projection | -0.2618 | -0.1027 | 214 | 209 |
| ARCHS4-pretrained + OSDR fine-tuned final | 0.9991 | 0.9691 | 67 | 50 |
| ARCHS4-only | -0.0675 | 0.7018 | n/a | n/a |

The fine-tuned ARCHS4 model fits the OSDR expression distribution much better
than the OSDR-only model. The frozen pretrained projection shows the strongest
condition separation, while the fine-tuned model is the better generator for
synthetic OSDR-like FLT/GC samples.

## Synthetic Examples

Small matched counterfactual examples were generated for all eight main tissues
from both OSDR-capable models:

- `outputs/wgan_conditional_generation/synthetic_examples/osdr_only/<tissue>/`
- `outputs/wgan_conditional_generation/synthetic_examples/archs4_pretrain_osdr_finetune/<tissue>/`

Each tissue directory contains:

- `ground_control_zscore.tsv.gz`, `ground_control_log1p_cpm.tsv.gz`,
  `ground_control_cpm.tsv.gz`
- `flight_zscore.tsv.gz`, `flight_log1p_cpm.tsv.gz`, `flight_cpm.tsv.gz`
- `flight_minus_ground_control_mean_log1p_cpm_delta.tsv`
- JSON files recording the conditioning profiles used

The example generation used `n=8` per condition, real observed OSDR covariate
profiles, and matched random noise for the ground-control and flight pair.

To generate more samples for a specific observed profile:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_wgan.generate_synthetic \
  --model-dir outputs/wgan_conditional_generation/archs4_pretrain_osdr_finetune \
  --output-dir outputs/wgan_conditional_generation/synthetic_custom/liver \
  --n 128 \
  --counterfactual ground_control flight \
  --set wgan_tissue=liver \
  --set wgan_material_type=Liver \
  --set wgan_muscle_group=not_skeletal_muscle \
  --set wgan_accession=OSD-137 \
  --set wgan_sex=Female \
  --set "wgan_assay=RNA Sequencing (RNA-Seq)" \
  --set wgan_platform=illumina_hiseq_4000 \
  --set wgan_data_source=osdr
```

## Interpretation Limits

These models are conditional generators, not causal estimators. A
ground-control-to-flight counterfactual holds observed metadata fixed and flips
the condition label, but unmeasured study design effects can still be embedded
in accession/platform/material covariates.

The generated matrices are continuous expression values in z-score,
log1p-CPM, and CPM spaces. They are not integer read counts. WGAN critic
features are learned representation features, not named pathways or Reactome
programs; gene/pathway interpretation still needs attribution or downstream
enrichment.
