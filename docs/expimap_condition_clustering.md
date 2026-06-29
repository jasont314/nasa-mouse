# expiMap Condition-Specific Clustering

This is the expiMap analogue of clustering ground-control and flight states
separately, then comparing the resulting structures.

Important difference from the GLARE-style gene-cluster workflow: these expiMap
tables cluster samples by Reactome pathway scores, not genes by latent
representation. Numeric cluster IDs are condition-specific, so comparisons are
made through nearest pathway-score centroids and accession composition.

## Direct Liver Model

Artifact:
`outputs/expimap_direct_osdr_liver/raw_counts_nb_50epoch/condition_cluster_comparison`

The direct posterior-mean liver pathway scores cluster into:

| condition | selected k | silhouette |
| --- | ---: | ---: |
| ground_control | 2 | 0.483 |
| flight | 3 | 0.440 |

The two main clusters mostly preserve the same accession structure in both
conditions:

| condition | cluster | samples | top accessions |
| --- | ---: | ---: | --- |
| ground_control | 0 | 85 | OSD-379, OSD-245, OSD-168, OSD-463, OSD-137 |
| flight | 0 | 83 | OSD-379, OSD-245, OSD-168, OSD-463, OSD-137 |
| ground_control | 1 | 28 | OSD-457, OSD-48, OSD-164, OSD-47, OSD-686 |
| flight | 1 | 30 | OSD-457, OSD-48, OSD-686, OSD-164, OSD-47 |

The extra flight cluster has only five samples, from OSD-168, OSD-379, and
OSD-48. Its large pathway deltas are therefore not a reliable flight-specific
program without accession-level replication.

## 5,000-Sample ARCHS4 Reference Query

Artifacts:

- `outputs/expimap_archs4_reference_osdr_query_liver/query_nb_5000stratified_seed2020_50epoch/condition_cluster_comparison`
- `outputs/expimap_archs4_reference_osdr_query_liver/query_nb_5000stratified_seed2021_50epoch/condition_cluster_comparison`
- `outputs/expimap_archs4_reference_osdr_query_liver/query_nb_5000stratified_seed2022_50epoch/condition_cluster_comparison`

Across the three reference-query seeds, flight repeatedly splits into one large
cluster and one small cluster:

| seed | GC k | FLT k | small FLT cluster |
| --- | ---: | ---: | --- |
| 2020 | 5 | 2 | 8 samples, mostly OSD-379 |
| 2021 | 4 | 2 | 7 samples, mostly OSD-379 |
| 2022 | 4 | 2 | 7 samples, mostly OSD-379 |

Ground-control k and pathway shifts vary across seeds. Several GC clusters are
accession-heavy or accession-specific, including OSD-245-only and
OSD-457/OSD-686-heavy clusters. This is consistent with the broader
reference-query seed-stability result: cluster structure is useful as QC and
hypothesis triage, but is not yet a stable biological conclusion.

## Interpretation

Condition-specific clustering is worth doing, but the current result is mostly
a study/accession-structure diagnostic. It does not rescue the unsupported
Polymerase II result, and it does not promote a new pathway to confirmed status.

Use this analysis to:

- identify accession-driven substructure before interpreting pooled pathway
  shifts;
- flag small condition-specific sample groups for QC;
- generate pathway candidates that still require accession-aware validation.

For a stronger version of this paper-style comparison, cluster within each
accession or require each condition-specific cluster signature to recur across
multiple accessions before ranking pathways biologically.
