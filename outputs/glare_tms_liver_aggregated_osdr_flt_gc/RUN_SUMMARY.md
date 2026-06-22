# Aggregated OSDR Liver FLT/GC GLARE Run

- Pretrained model: `outputs/glare_paper_tms_liver_osd379/pretraining/sc_shulse_pretrained_reproduced.pth`
- Target expression: `data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json`
- Tissue field: `/meta/samples/characteristics/study.characteristics.material type`
- Condition field: `/meta/samples/factors/study.factor value.spaceflight`
- Included conditions: `Space Flight` and `Ground Control`
- Included accessions: `OSD-379`, `OSD-245`, `OSD-463`, `OSD-242`, `OSD-137`, `OSD-47`, `OSD-686`, `OSD-173`
- Gene universe: 21,010 genes shared by TMS liver pretraining and OSDR liver fine-tuning
- Profiles: 91 FLT and 82 GC

## Per-Accession Counts

| Accession | FLT | GC | Total |
| --- | ---: | ---: | ---: |
| OSD-379 | 35 | 35 | 70 |
| OSD-245 | 20 | 19 | 39 |
| OSD-463 | 12 | 10 | 22 |
| OSD-242 | 5 | 4 | 9 |
| OSD-137 | 6 | 6 | 12 |
| OSD-47 | 5 | 3 | 8 |
| OSD-686 | 6 | 3 | 9 |
| OSD-173 | 2 | 2 | 4 |

## Fine-Tuning

| Location | Profiles | Best loss | Best epoch | Elapsed |
| --- | ---: | ---: | ---: | ---: |
| FLT | 91 | 0.15166330 | 27 | 1m 15s |
| GC | 82 | 0.14730988 | 30 | 1m 20s |

## Clustering

GLARE paper-style GMM/HDBSCAN/Spectral clustering followed by exact
evidence-accumulation consensus clustering was run for both FLT and GC.

| Location | Consensus clusters | Silhouette | Notes |
| --- | ---: | ---: | --- |
| FLT | 16 | 0.1168 | HDBSCAN found 13 hard clusters with 19,245 initial noise points |
| GC | 15 | -0.2033 | HDBSCAN found 9 hard clusters with 20,132 initial noise points |

The negative GC silhouette means the GC consensus partition should be treated
as exploratory and checked for study/batch structure before biological
interpretation.
