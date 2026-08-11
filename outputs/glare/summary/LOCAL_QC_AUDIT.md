# Local QC Audit

Review-agent audit could not complete because the subagent hit a usage limit. This local audit checks the same output classes from disk.

| tissue | status | FLT/GC | studies | per-study | DGEA | MOBER | plots | sil direct FLT/GC | sil MOBER FLT/GC | caveats |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kidney | PASS | 67/66 | 6 | 6/6 | yes | yes | 20 legend; 4 UMAP | 0.0915/0.0514 | 0.323/0.3 | none |
| liver | PASS | 117/112 | 12 | 12/12 | yes | yes | 32 legend; 4 UMAP | 0.0825/0.142 | 0.342/0.369 | none |
| lung | PASS | 39/37 | 3 | 3/3 | yes | yes | 14 legend; 4 UMAP | 0.238/0.161 | 0.302/0.304 | none |
| skeletal muscle | PASS | 92/93 | 13 | 13/13 | yes | yes | 34 legend; 4 UMAP | -0.00706/-0.0991 | 0.42/0.368 | FLT aggregate silhouette negative (-0.00706); GC aggregate silhouette negative (-0.0991) |
| skeletal muscle: edl | PASS | 15/15 | 2 | 2/2 | yes | yes | 12 legend; 4 UMAP | -0.0318/0.0996 | -0.0942/-0.098 | FLT aggregate silhouette negative (-0.0318); MOBER FLT silhouette negative (-0.0942); MOBER GC silhouette negative (-0.098); small muscle-subtype MOBER fit showed high VAE loss/early checkpoint in run log; treat MOBER branch as sensitivity only |
| skeletal muscle: gastrocnemius | PASS | 13/17 | 3 | 3/3 | yes | yes | 14 legend; 4 UMAP | 0.203/0.0998 | -0.0651/-0.101 | MOBER FLT silhouette negative (-0.0651); MOBER GC silhouette negative (-0.101); small muscle-subtype MOBER fit showed high VAE loss/early checkpoint in run log; treat MOBER branch as sensitivity only |
| skeletal muscle: quadriceps | PASS | 22/22 | 4 | 4/4 | yes | yes | 16 legend; 4 UMAP | 0.19/0.139 | 0.0344/0.14 | small muscle-subtype MOBER fit showed high VAE loss/early checkpoint in run log; treat MOBER branch as sensitivity only |
| skeletal muscle: soleus | PASS | 27/24 | 3 | 3/3 | yes | yes | 14 legend; 4 UMAP | 0.0385/-0.117 | 0.286/0.106 | GC aggregate silhouette negative (-0.117); small muscle-subtype MOBER fit showed high VAE loss/early checkpoint in run log; treat MOBER branch as sensitivity only |
| skeletal muscle: tibialis anterior | PASS | 15/15 | 2 | 2/2 | yes | yes | 12 legend; 4 UMAP | 0.174/0.162 | -0.0295/-0.0247 | MOBER FLT silhouette negative (-0.0295); MOBER GC silhouette negative (-0.0247); small muscle-subtype MOBER fit showed high VAE loss/early checkpoint in run log; treat MOBER branch as sensitivity only |
| skin | PASS | 80/71 | 6 | 6/6 | yes | yes | 20 legend; 4 UMAP | 0.266/-0.0343 | 0.325/0.334 | GC aggregate silhouette negative (-0.0343) |
| spleen | PASS | 55/53 | 6 | 6/6 | yes | yes | 20 legend; 4 UMAP | 0.159/0.1 | 0.377/0.38 | none |
| thymus | PASS | 62/53 | 5 | 5/5 | yes | yes | 18 legend; 4 UMAP | -0.0221/0.0797 | 0.335/0.319 | FLT aggregate silhouette negative (-0.0221) |

## Notes

- `status=PASS` means required aggregate/per-study GLARE, raw-count availability, cluster files, DGEA comparison outputs, aggregate MOBER outputs, and legend plots are present and internally consistent.
- Negative silhouette is a QC caveat for cluster geometry, not an automatic failure.
- Small skeletal-muscle subtype MOBER runs completed, but several had high MOBER VAE loss or epoch-0/early checkpoints; use direct aggregate and per-study results as primary evidence for those subtypes.
- PCA legend plots were generated for aggregate, MOBER, and per-study scopes. UMAP legend plots were generated for aggregate and MOBER scopes.
