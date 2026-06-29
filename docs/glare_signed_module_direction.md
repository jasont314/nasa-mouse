# GLARE Signed Module Direction

This note records the signed within-module DGEA direction check for the expanded
15-term validation stack.

## Command

```bash
PYTHONPATH=src /opt/anaconda3/envs/nasa/bin/python -m nasa_mouse_glare.module_direction_summary \
  --validation-dir outputs/glare_multi_tissue_api/validation_stack_terms15
```

## Outputs

Generated files are under:

```text
outputs/glare_multi_tissue_api/validation_stack_terms15/module_direction/
```

Key files:

- `candidate_module_signed_dgea_by_study.tsv`
- `candidate_module_signed_dgea_meta.tsv`
- `SIGNED_MODULE_DIRECTION_SUMMARY.md`

`outputs/` is ignored by git, so this document tracks the interpretation.

## Why This Was Added

Reactome/Metascape enrichment tells us that a module contains genes from a
pathway, but it does not tell us whether that pathway is activated, suppressed,
or internally reorganized. A single module can contain both FLT-up and GC-up
genes from the same pathway.

The GLARE paper used related validation layers: DEG proportion in clusters,
GO/Metascape enrichment, expression heatmaps, biological network/TF follow-up,
cell-type prediction, and SHAP-style feature follow-up. It did not present this
exact signed FLT-up/GC-up/mixed table. This table is a stricter extension for
our multi-study mouse setting.

## Direction Calls

Positive `log2FoldChange` means higher in spaceflight. Negative
`log2FoldChange` means higher in ground control.

| Module class | Direction call | Modules |
| --- | --- | ---: |
| GLARE-only | mixed/reorganized DEG | 101 |
| GLARE-only | ambiguous/no validated direction | 35 |
| GLARE-only | FLT-up DEG-supported | 23 |
| GLARE-only | GC-up DEG-supported | 21 |
| DGEA-intersection | GC-up DEG-supported | 84 |
| DGEA-intersection | FLT-up DEG-supported | 53 |
| DGEA-intersection | mixed/reorganized DEG | 38 |
| DGEA-intersection | ambiguous/no validated direction | 5 |

Important interpretation:

- `DGEA-intersection` means the pathway recurred in both DGEA pathway analysis
  and GLARE cluster enrichment.
- `GLARE-only` means the pathway was selected from recurring GLARE enrichment
  but not from recurring DGEA pathway overlap. It can still contain individual
  DE genes.
- `mixed/reorganized DEG` means the module contains both FLT-up and GC-up
  significant genes, so enrichment should not be interpreted as simple pathway
  activation or suppression.

## Main Interpretation

The signed check makes the current result more cautious:

- Many validated modules are directional, especially DGEA-intersection modules.
- Many GLARE-only modules are not clean hidden activation/suppression signals;
  they are mixed/reorganized pathway groupings.
- Skeletal muscle respiratory modules are directionally coherent: mostly
  GC-higher / FLT-suppressed.
- Thymus platelet-calcium signaling is FLT-up DEG-supported.
- Several thymus and skeletal muscle GLARE-only modules remain mixed, including
  cell-death, circadian, and some signaling modules.
- Liver hidden-module support remains weak. The olfactory module is still
  ambiguous/no validated direction.

## Examples

| Tissue | Module class | Term | Signed result |
| --- | --- | --- | --- |
| skeletal_muscle | GLARE-only | Respiratory electron transport / ATP synthesis | GC-up DEG-supported |
| skeletal_muscle_soleus | GLARE-only | Respiratory electron transport / ATP synthesis | GC-up DEG-supported |
| thymus | GLARE-only | Response to elevated platelet cytosolic Ca2 | FLT-up DEG-supported |
| thymus | GLARE-only | Cell death signalling via NRAGE/NRIF/NADE | mixed/reorganized DEG |
| skeletal_muscle | GLARE-only | BMAL1:CLOCK/NPAS2 circadian expression | mixed/reorganized DEG |
| kidney | GLARE-only | Signaling by insulin receptor | mixed/reorganized DEG |
| liver | GLARE-only | Cleavage of growing transcript in the termination region | mixed/reorganized DEG |
| liver | GLARE-only | Olfactory signaling pathway | ambiguous/no validated direction |

## Practical Rule

For biological interpretation, use:

1. Pathway enrichment for module identity.
2. Within-module signed DGEA counts for direction.
3. Module-score FLT-GC effect and direction consistency for hidden/non-DEG
   signal.
4. Random gene-set empirical p-values to reject generic large-gene-set effects.

Do not infer activation or suppression from enrichment alone.
