# Thymus Platelet-Calcium Validation

This note records the follow-up validation for the thymus GLARE-only Reactome
module:

`REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_`

## Command

```bash
PYTHONPATH=src /opt/anaconda3/envs/nasa/bin/python -m nasa_mouse_glare.validate_thymus_platelet_calcium
```

Generated files are under:

```text
outputs/glare_multi_tissue_api/validation_stack_terms15/thymus_platelet_calcium_validation/
```

`outputs/` is ignored by Git, so this document tracks the conclusion.

## Bottom Line

The signal is real enough to keep as a GLARE finding, but it is not clean
evidence for thymocyte-intrinsic biology.

Best label:

```text
FLT-up platelet/hemostasis/endothelial-remodeling module with substantial
composition risk
```

## Key Checks

- Significant Reactome enrichment appears in all five thymus accessions:
  `OSD-244`, `OSD-289`, `OSD-421`, `OSD-457`, and `OSD-515`.
- Strict module-score FLT-up support appears in three accessions:
  `OSD-244`, `OSD-289`, and `OSD-457`.
- No accession has strict GC-up support.
- `OSD-421` trends GC-higher but is not strict-significant.
- `OSD-515` is effectively flat.

The Reactome v4 mouse Ensembl term has 85 mapped IDs. In the thymus GLARE
outputs, 82 IDs appear across all cluster rows and 77 IDs appear in
FDR-significant cluster rows.

## Marker Evidence

Panglao marker enrichment strongly supports a composition/vascular component:

| Marker set | Overlap | FDR |
| --- | ---: | ---: |
| Platelets | 18 | `2.73e-23` |
| Megakaryocytes | 10 | `5.18e-15` |
| Endothelial cells | 11 | `6.41e-10` |

Manual marker checks show:

| Category | Genes | Significant pairs | FLT-up pairs | GC-up pairs |
| --- | ---: | ---: | ---: | ---: |
| Panglao platelet | 18 | 23 | 21 | 2 |
| Coagulation/plasma | 14 | 23 | 14 | 9 |
| ECM/growth-factor remodeling | 12 | 21 | 19 | 2 |
| Calcium/secretion/signaling | 11 | 11 | 11 | 0 |
| Endothelial/vascular | 10 | 16 | 14 | 2 |
| Thymic epithelial/stromal | 0 | 0 | 0 | 0 |
| Thymocyte/T cell | 0 | 0 | 0 | 0 |

Top recurrent FLT-up genes include `PROS1`, `TIMP1`, `PRKCG`, `PDGFA`,
`TGFB3`, `TLN1`, `MMRN1`, `TF`, `FN1`, `THBS1`, `FGA`, `PECAM1`, and `PLG`.

## Interpretation

Evidence against a one-accession artifact is reasonably good: the term recurs
across studies and three accessions have strict FLT-up module-score support.

Evidence for composition or vascular biology is stronger than evidence for
direct thymocyte biology. The module is enriched for platelet,
megakaryocyte, endothelial, coagulation/plasma, and remodeling genes, while
canonical thymocyte and thymic epithelial markers are absent from the module.

Practical interpretation:

```text
This is a recurring spaceflight-associated hemostasis/platelet-calcium/
endothelial-remodeling signature in thymus bulk RNA-seq. It could reflect
vascular remodeling, platelet/blood content, cell-composition shifts, or a
mixture of those effects.
```

Do not call it simply "thymus platelet activation" without additional
cell-composition or histology/flow validation.
