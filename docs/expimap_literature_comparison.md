# expiMap Literature Comparison

Date: 2026-07-01

This note compares the current expiMap FLT-vs-GC OSDR mouse results with
published spaceflight mouse/OSDR literature through 2026-07-01. It is limited
to expiMap outputs and does not use GLARE results.

## Result Summary

Counts below are accession-aware random-effects pathway tests at FDR < 0.05.
LOO means the pathway also passed leave-one-accession-out FDR < 0.05 in the
same direction.

| tissue | direct meta-FDR terms | direct LOO-pass terms | ARCHS4 ref-query meta-FDR terms | ARCHS4 ref-query LOO-pass terms | interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| liver | 381 | 15 | 6 | 0 | Direct has signal; ref-query is weaker. The more stable 5k liver reference nominated semaphorin interaction as a follow-up candidate, not the canonical lipid/metabolic theme. |
| kidney | 0 | 0 | 0 | 0 | Literature predicts renal lipid/ECM/TGF-beta/oxidative biology, but expiMap did not recover a stable pathway-level FLT-vs-GC signal. |
| skeletal_muscle | 2 | 0 | 0 | 0 | Aggregate muscle is weak. Split-muscle and targeted-module runs align better with prior muscle biology but mostly fail strict LOO stability. |
| skin | 0 | 0 | 208 | 0 | Reference-query has broad candidates, but none are LOO-stable. Treat as follow-up, not final biology. |
| thymus | 953 | 221 | 725 | 466 | Strongest expiMap result. This is the clearest match to prior immune/lymphoid and cell-cycle/DNA-replication literature. |
| spleen | 194 | 5 | 104 | 4 | Moderate robust immune signal, with fewer LOO-stable terms than thymus. |
| lung | 0 | 0 | 0 | 0 | Literature supports lung remodeling/injury signals, but current expiMap did not recover them. |
| retina | 0 | 0 | 0 | 0 | Literature supports retina/phototransduction/oxidative effects, but current expiMap did not recover stable pathway shifts. |

## Literature Alignment

| tissue | prior literature expectation | expiMap relationship |
| --- | --- | --- |
| liver | Lipid accumulation, fatty-acid/lipid processing, bile/xenobiotic/metabolic and proteostasis stress are repeatedly reported in mouse liver spaceflight datasets. | Partial/complementary. Direct OSDR expiMap detects many accession-aware terms, but the stronger liver reference-query follow-up does not cleanly recover the canonical lipid theme. The semaphorin candidate is best treated as a new hypothesis. |
| skeletal_muscle | Spaceflight muscle work emphasizes atrophy, calcium/SERCA dysregulation, contractile remodeling, mitochondrial metabolism, proteostasis/autophagy, ECM, and IGF/AKT-related signaling, often muscle-type-specific. | Candidate-level support. Aggregate skeletal muscle is weak, but split-muscle targeted expiMap finds plausible directions: gastrocnemius atrophy/slow-contractile changes, soleus ECM/contractile/metabolic changes, quadriceps inflammation/oxidative stress, and exploratory tibialis IGF/AKT/mTOR. None pass the stricter targeted LOO rule. |
| skin | Skin literature points to DNA damage/repair, mitochondrial dysregulation, barrier function, collagen/ECM and oxidative stress. | Candidate-level support only. ARCHS4 ref-query skin has many FDR candidates but no LOO-stable pathway, so it should be plotted and followed up before calling biology. |
| thymus | Prior mouse studies report strong thymus/spleen immune disruption, lymphoid effects, DNA fragmentation and cell-cycle/transcriptional changes. | Strong support. Both direct and ref-query expiMap show many accession-aware and LOO-stable thymus terms. This is currently the most defensible expiMap tissue. |
| spleen | Prior work supports adaptive/innate immune disruption and lymphoid remodeling. | Moderate support. Spleen has many meta-FDR pathways and a small LOO-stable core in both direct and ref-query runs. |
| kidney | Recent kidney work highlights strain-dependent transcriptomic shifts in lipid metabolism, ECM remodeling, TGF-beta signaling, oxidative stress and renal-risk pathways. | Not recovered. Current expiMap is negative in both direct and ref-query modes, so we should not claim kidney pathway biology from expiMap yet. |
| lung | Mouse lung studies report ECM/adhesion/profibrotic, injury/remodeling, immune, circadian/protein-folding and receptor-family changes. | Not recovered. No stable expiMap pathway signal in the current lung runs. |
| retina | Retina work reports phototransduction/visual-perception, retinal disease genes, oxidative stress, photoreceptor integrity and retinal-layer thinning after flight. | Not recovered. No stable expiMap pathway signal in the current retina runs. |

## Best-Supported vs Complementary Calls

Best-supported by both prior literature and expiMap:

- Thymus: immune/lymphoid and cell-cycle/DNA-replication biology.
- Spleen: immune and lymphoid pathway remodeling, weaker than thymus.

Plausible but not final:

- Skeletal muscle split analyses: expected muscle biology appears when muscle
  type is separated, but targeted modules do not pass strict LOO FDR.
- Skin: reference-query candidates match skin stress themes, but fail LOO.

Complementary/new hypotheses:

- Liver semaphorin interaction from the larger ARCHS4 liver reference-query
  stability follow-up. This is not the main published liver lipid/metabolism
  theme and needs independent validation.
- Muscle-type-specific expiMap candidates, especially soleus and gastrocnemius,
  are biologically plausible but need stronger accession stability.

Not recovered by current expiMap:

- Kidney, lung and retina. Existing literature supports real spaceflight
  biology in these tissues, but this method/run did not produce stable
  pathway-level FLT-vs-GC calls.

## Figure Inventory

Every listed analysis directory contains:

- `pathway_score_pca.png`: pathway/latent-score PCA colored by FLT vs GC.
- `pathway_score_umap.png`: pathway/latent-score UMAP colored by FLT vs GC.
- `pathway_score_pca_by_accession.png`: PCA colored by accession.
- `pathway_score_umap_by_accession.png`: UMAP colored by accession.
- `top_pathway_shift_heatmap.png`: heatmap of top pathway/program shifts.

Primary direct OSDR expiMap plot directories:

- `outputs/expimap/runs/direct/liver/raw_counts_nb_50epoch/analysis/`
- `outputs/expimap/runs/direct/kidney/raw_counts_nb_50epoch/analysis/`
- `outputs/expimap/runs/direct/skeletal_muscle/raw_counts_nb_50epoch/analysis/`
- `outputs/expimap/runs/direct/skin/raw_counts_nb_50epoch/analysis/`
- `outputs/expimap/runs/direct/thymus/raw_counts_nb_50epoch/analysis/`
- `outputs/expimap/runs/direct/spleen/raw_counts_nb_50epoch/analysis/`
- `outputs/expimap/runs/direct/lung/raw_counts_nb_50epoch/analysis/`
- `outputs/expimap/runs/direct/retina/raw_counts_nb_50epoch/analysis/`

Primary ARCHS4 reference-query expiMap plot directories:

- `outputs/expimap/runs/reference_query/liver/query_nb_5000stratified_seed2020_50epoch/analysis/`
- `outputs/expimap/runs/reference_query/kidney/query_nb_1000ref_50epoch/analysis/`
- `outputs/expimap/runs/reference_query/skeletal_muscle/query_nb_allref_50epoch/analysis/`
- `outputs/expimap/runs/reference_query/skin/query_nb_allref_50epoch/analysis/`
- `outputs/expimap/runs/reference_query/thymus/query_nb_allref_50epoch/analysis/`
- `outputs/expimap/runs/reference_query/spleen/query_nb_allref_50epoch/analysis/`
- `outputs/expimap/runs/reference_query/lung/query_nb_allref_50epoch/analysis/`
- `outputs/expimap/runs/reference_query/retina/query_nb_allref_50epoch/analysis/`

Split-muscle targeted expiMap plot directories:

- `outputs/expimap/runs/muscle_groups/combined_min8/direct_gastrocnemius_nb_100epoch/analysis/`
- `outputs/expimap/runs/muscle_groups/combined_min8/direct_quadriceps_nb_100epoch/analysis/`
- `outputs/expimap/runs/muscle_groups/combined_min8/direct_soleus_nb_100epoch/analysis/`
- `outputs/expimap/runs/muscle_groups/combined_min8/query_gastrocnemius_nb_allref_50epoch/analysis/`
- `outputs/expimap/runs/muscle_groups/combined_min8/query_quadriceps_nb_allref_50epoch/analysis/`
- `outputs/expimap/runs/muscle_groups/combined_min8/query_soleus_nb_allref_50epoch/analysis/`
- `outputs/expimap/runs/muscle_groups/combined_min8/direct_edl_nb_100epoch/analysis/`
- `outputs/expimap/runs/muscle_groups/combined_min8/direct_tibialis_anterior_nb_100epoch/analysis/`
- `outputs/expimap/runs/muscle_groups/combined_min8/query_edl_nb_allref_50epoch/analysis/`
- `outputs/expimap/runs/muscle_groups/combined_min8/query_tibialis_anterior_nb_allref_50epoch/analysis/`

## Representative Literature Links

- Liver lipid/metabolic remodeling: Beheshti et al. 2019, Scientific Reports,
  https://doi.org/10.1038/s41598-019-55869-2
- Liver lipotoxic/PPAR response after STS-135: Jonscher et al. 2016,
  https://doi.org/10.1371/journal.pone.0152877
- Muscle calcium/SERCA and multi-omics signatures: Li et al. 2023,
  https://www.nature.com/articles/s41526-023-00337-5
- Muscle transcriptome/GSEA and atrophy systems biology: Oommen et al. 2024,
  https://www.nature.com/articles/s41526-024-00434-z
- Skin dermatology, barrier/collagen/DNA repair/mitochondria: Cope et al. 2024,
  https://www.nature.com/articles/s43856-024-00532-9
- Kidney lipid/ECM/TGF-beta strain-dependent response: Finch et al. 2025,
  https://doi.org/10.1038/s41526-025-00465-0
- Thymus/spleen STS-135 immune/cell-cycle effects: Gridley et al. 2013,
  https://doi.org/10.1371/journal.pone.0075097
- Lung ECM/profibrotic response: Tian et al. 2010,
  https://doi.org/10.1152/japplphysiol.00730.2009
- Retina phototransduction/oxidative stress/photoreceptor changes: Overbey et
  al. 2019, https://www.nature.com/articles/s41598-019-49453-x
