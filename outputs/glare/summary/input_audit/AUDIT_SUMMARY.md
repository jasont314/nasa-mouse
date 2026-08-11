# Multi-Tissue API GLARE Input Audit

This audit uses NASA OSDR Biological Data API metadata and TMS FACS h5ad metadata.

- API metadata: `data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv`
- TMS h5ad: `assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad`
- Minimum TMS cells: 100

## Requested Tissue Status

| tissue_slug | label | tissue_final | material_terms | tms_tissue | tms_cells | space_flight | ground_control | accessions | mober_eligible_ge2_studies | pretraining_status | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| liver | liver | liver |  | liver | 2859 | 125 | 118 | 12 | True | ok |  |
| skeletal_muscle | skeletal muscle | skeletal_muscle |  | limb muscle | 3855 | 95 | 96 | 13 | True | ok |  |
| skin | skin | skin |  | skin of body | 4860 | 80 | 71 | 6 | True | ok |  |
| kidney | kidney | kidney |  | kidney | 1833 | 68 | 67 | 6 | True | ok |  |
| thymus | thymus | thymus |  | thymus | 4047 | 63 | 54 | 5 | True | ok |  |
| spleen | spleen | spleen |  | spleen | 3834 | 55 | 54 | 6 | True | ok |  |
| lung | lung | lung |  | lung | 5218 | 40 | 38 | 3 | True | ok |  |
| retina | retina | retina |  |  | 0 | 45 | 31 | 4 | True | skip_no_matching_or_too_few_tms_cells | No matching TMS FACS tissue in current h5ad |
| skeletal_muscle_soleus | skeletal muscle: soleus | skeletal_muscle | Right soleus; Soleus; Soleus-both sides | limb muscle | 3855 | 28 | 25 | 3 | True | ok | OSDR sub-tissue run; TMS pretraining uses combined limb muscle FACS |
| skeletal_muscle_quadriceps | skeletal muscle: quadriceps | skeletal_muscle | Left quadriceps femoris; Quadriceps femoris; Right quadriceps femoris | limb muscle | 3855 | 23 | 23 | 4 | True | ok | OSDR sub-tissue run; TMS pretraining uses combined limb muscle FACS |
| skeletal_muscle_gastrocnemius | skeletal muscle: gastrocnemius | skeletal_muscle | Gastrocnemius; Left gastrocnemius; Right gastrocnemius | limb muscle | 3855 | 13 | 17 | 3 | True | ok | OSDR sub-tissue run; TMS pretraining uses combined limb muscle FACS |
| skeletal_muscle_edl | skeletal muscle: edl | skeletal_muscle | Extensor digitorum longus- both sides; Right extensor digitorum longus | limb muscle | 3855 | 16 | 16 | 2 | True | ok | OSDR sub-tissue run; TMS pretraining uses combined limb muscle FACS |
| skeletal_muscle_tibialis_anterior | skeletal muscle: tibialis anterior | skeletal_muscle | Left tibialis anterior; Right tibialis anterior | limb muscle | 3855 | 15 | 15 | 2 | True | ok | OSDR sub-tissue run; TMS pretraining uses combined limb muscle FACS |

## Notes

- MOBER is marked eligible when at least two OSDR studies are available.
- Skeletal muscle sub-tissue runs use exact OSDR material-type labels.
- Retina has OSDR bulk FLT/GC data but no matching TMS FACS tissue in the current h5ad, so it is skipped for GLARE pretraining.
