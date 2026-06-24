# OSDR Mouse Bulk RNA-seq FLT/GC Counts By Tissue

Source metadata: `/Users/jasontrinh/Desktop/Code/Berkeley/nasa/codex-analysis/osdr_mus_musculus_transcription_profiling/ml_ready/bulk_rna_seq_sample_level_spaceflight_vs_ground_control/bulk_rna_seq_sample_level_spaceflight_vs_ground_control_sample_metadata.tsv`

Scope: Mus musculus sample-level bulk RNA-seq Spaceflight Study profiles included in the local `bulk_rna_seq_sample_level_spaceflight_vs_ground_control` matrix. Counts are expression-profile/sample rows, not guaranteed unique animals.

Total included profiles: 1,382; spaceflight: 699; ground control: 683; included studies: 68; normalized tissues: 24.

## Ranked Tissues

| Rank | Tissue | FLT | GC | Total | Studies | Accessions |
| ---: | --- | ---: | ---: | ---: | ---: | --- |
| 1 | liver | 102 | 98 | 200 | 9 | OSD-137,OSD-168,OSD-173,OSD-245,OSD-379,OSD-463,OSD-47,OSD-48,OSD-686 |
| 2 | kidney | 71 | 105 | 176 | 6 | OSD-102,OSD-163,OSD-253,OSD-462,OSD-513,OSD-771 |
| 3 | dorsal_skin | 60 | 55 | 115 | 4 | OSD-238,OSD-240,OSD-243,OSD-254 |
| 4 | thymus | 46 | 40 | 86 | 4 | OSD-244,OSD-289,OSD-421,OSD-515 |
| 5 | retina | 45 | 31 | 76 | 4 | OSD-194,OSD-255,OSD-397,OSD-758 |
| 6 | lung | 39 | 37 | 76 | 3 | OSD-248,OSD-464,OSD-900 |
| 7 | spleen | 39 | 37 | 76 | 4 | OSD-246,OSD-288,OSD-420,OSD-506 |
| 8 | colon | 28 | 27 | 55 | 2 | OSD-247,OSD-667 |
| 9 | cerebellum | 28 | 26 | 54 | 3 | OSD-525,OSD-561,OSD-563 |
| 10 | soleus_muscle | 27 | 24 | 51 | 3 | OSD-104,OSD-714,OSD-770 |
| 11 | heart | 26 | 25 | 51 | 3 | OSD-270,OSD-580,OSD-599 |
| 12 | adrenal_gland | 20 | 19 | 39 | 3 | OSD-161,OSD-512,OSD-98 |
| 13 | quadriceps_muscle | 19 | 19 | 38 | 3 | OSD-103,OSD-326,OSD-666 |
| 14 | hippocampus | 20 | 17 | 37 | 2 | OSD-562,OSD-564 |
| 15 | femoral_skin | 20 | 16 | 36 | 2 | OSD-239,OSD-241 |
| 16 | optic_nerve | 23 | 12 | 35 | 1 | OSD-759 |
| 17 | extensor_digitorum_longus_muscle | 15 | 15 | 30 | 2 | OSD-665,OSD-99 |
| 18 | tibialis_anterior_muscle | 15 | 15 | 30 | 2 | OSD-105,OSD-576 |
| 19 | mammary_gland | 10 | 15 | 25 | 1 | OSD-511 |
| 20 | bone_marrow | 12 | 12 | 24 | 1 | OSD-690 |
| 21 | gastrocnemius_muscle | 10 | 14 | 24 | 2 | OSD-101,OSD-419 |
| 22 | eye | 11 | 11 | 22 | 2 | OSD-100,OSD-162 |
| 23 | cecum | 10 | 10 | 20 | 1 | OSD-899 |
| 24 | spermatogonia | 3 | 3 | 6 | 1 | OSD-901 |

## Tissue-System View

| Rank | Tissue system | FLT | GC | Total | Tissues | Studies |
| ---: | --- | ---: | ---: | ---: | --- | ---: |
| 1 | metabolic_digestive | 102 | 98 | 200 | liver | 9 |
| 2 | renal | 71 | 105 | 176 | kidney | 6 |
| 3 | skeletal_muscle | 86 | 87 | 173 | extensor_digitorum_longus_muscle,gastrocnemius_muscle,quadriceps_muscle,soleus_muscle,tibialis_anterior_muscle | 12 |
| 4 | immune | 85 | 77 | 162 | spleen,thymus | 8 |
| 5 | integumentary | 80 | 71 | 151 | dorsal_skin,femoral_skin | 6 |
| 6 | visual | 79 | 54 | 133 | eye,optic_nerve,retina | 7 |
| 7 | nervous_system | 48 | 43 | 91 | cerebellum,hippocampus | 5 |
| 8 | respiratory | 39 | 37 | 76 | lung | 3 |
| 9 | gastrointestinal | 38 | 37 | 75 | cecum,colon | 3 |
| 10 | cardiovascular | 26 | 25 | 51 | heart | 3 |
| 11 | endocrine | 20 | 19 | 39 | adrenal_gland | 3 |
| 12 | reproductive | 10 | 15 | 25 | mammary_gland | 1 |
| 13 | hematopoietic_immune | 12 | 12 | 24 | bone_marrow | 1 |
| 14 | reproductive_cells | 3 | 3 | 6 | spermatogonia | 1 |

## Caveats

- `tissue` is the normalized tissue label from the local codex-analysis matrix. Raw aliases such as left/right lobes are collapsed when the upstream analysis normalized them.
- Basal/vivarium/other controls are excluded here; only rows classified as `spaceflight` or `ground_control` are counted.
- Count matrices include VST/normalized GLbulkRNAseq and normalized processed tables; 21 rows are from studies flagged with technical replicates, and many studies use ERCC-processed matrices. Treat counts as available expression profiles, not unique animals.
- Skeletal muscle appears split by muscle type at the tissue level. If aggregated to `skeletal_muscle`, it has 86 FLT and 87 GC profiles, 173 total.
