# Output Directory

Generated artifacts are grouped by workflow and purpose. New runs should use
this layout instead of adding another prefixed directory directly under
`outputs/`.

For final numbers and plots, start with the
[results guide](../docs/results_guide.md). This directory records run selection
and provenance; the paper `source_data/` tables are the final numeric record.

## Layout

| Path | Contents |
|---|---|
| `osdr/inventory/` | OSDR sample inventories and tissue/material summaries |
| `glare/` | GLARE and MOBER runs when regenerated |
| `expimap/runs/direct/` | OSDR-only expiMap runs, organized by tissue |
| `expimap/runs/reference_query/` | ARCHS4 reference to OSDR query runs |
| `expimap/runs/muscle_groups/` | Muscle-group expiMap runs |
| `expimap/analyses/` | Cross-run summaries, follow-up analyses, and figures |
| `expimap/summary/` | Compact multi-tissue result tables |
| `generative/benchmark/` | Unified DDIM/WGAN benchmark, data audit, runs, and analyses |
| `generative/standalone/` | Regeneration target for the earlier direct diffusion and WGAN runners |

Large data, checkpoints, and generated samples remain local and are ignored by
Git. Curated inventories and summary tables already versioned by the project
remain trackable.

## Final Analysis Selections

The paths below are the runs consumed by the current paper, report, or
presentation builders. Other experiments remain documented in the command
ledger, but their local run directories are regenerable and are not retained.

### expiMap

Core HVG reference-query results:

- Thymus:
  `expimap/runs/reference_query/thymus/tutorial_hvg_2000/query_nb_250epoch_seed2020/`
- Skin:
  `expimap/runs/reference_query/skin/tutorial_hvg_2000/query_nb_250epoch_seed2020/`
- Liver primary deduplicated analysis:
  `expimap/runs/reference_query/liver/tutorial_hvg_2000/query_nb_250epoch_seed2020_primary_deduplicated/`
- Soleus:
  `expimap/runs/muscle_groups/combined_min8/tutorial_hvg_soleus_2000/query_nb_250epoch_seed2020/`
- Kidney and spleen reassessment:
  `expimap/analyses/kidney_spleen_reassessment/`

The command ledger preserves the direct OSDR, all-gene, de novo, alternate-seed,
and preprocessing experiments. Their local run trees were removed after the
final selections were frozen into paper source tables.

### Generative Models

Selected ARCHS4 DDIM backbone:

`generative/benchmark/runs/lacan_diffusion/archs4_mouse_paper_parity_osdr_disjoint_seed1234/`

Selected OSDR-conditioned DDIM adapter:

`generative/benchmark/runs/lacan_diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint_seed2020/`

WGAN comparison run:

`generative/benchmark/runs/vinas_wgan_gp/osdr_matched_study_conditioned_seed2020/`

The DDIM was selected for downstream synthetic analysis. The WGAN remains a
model comparison, not the source of the final synthetic profiles.

Selected downstream analyses:

- Tissue and muscle consensus screens:
  `generative/benchmark/analyses/within_study_generated_feature_stability_osdr_disjoint_v1/`
  and
  `generative/benchmark/analyses/within_study_generated_feature_stability_muscle_groups_osdr_disjoint_v1/`
- Real versus real-plus-synthetic classifier analysis:
  `generative/benchmark/analyses/matched_all_gene_classifiers_osdr_disjoint_v1/`
- Individual-gene permutation and SHAP analysis:
  `generative/benchmark/analyses/classifier_importance_osdr_disjoint_v1/`
- Reactome-grouped importance analysis:
  `generative/benchmark/analyses/grouped_pathway_importance_osdr_disjoint_v1/`

The manuscript uses frozen tables under
`paper/synthetic_guided_spaceflight/source_data/`. Those tables are publication
artifacts derived from the selected output directories above.

## Local Artifact Retention

Only model artifacts required by the selections above are kept locally:

- the selected ARCHS4 DDIM backbone;
- the selected final OSDR DDIM adapter;
- the selected WGAN model, reference model, and best checkpoints;
- the thymus, skin, liver, and soleus expiMap query models;
- the reference models required by those four expiMap queries.

`MODEL_ARTIFACTS.sha256` records hashes for this local allowlist, including the
expiMap gene-order metadata needed to reload each model. Run
`sha256sum --check outputs/MODEL_ARTIFACTS.sha256` to verify it. The models are
ignored by Git, so the check is intended for the machine that holds the local
artifacts.

The selected run inputs and scored expiMap objects needed by the paper builders
remain on the handoff machine. Other prepared data, caches, rejected runs, and
smoke runs can be recreated from the configs and commands in `COMMANDS.md`.
The current OSDR inventory is generated from the NASA OSDR API under
`generative/benchmark/data_audit/osdr/`.

The standalone output tree is not part of the final analysis and is not kept
locally. Its older runners remain available for reproducing the development
experiments recorded in the command ledger.

See [COMMANDS.md](COMMANDS.md) for the command ledger and reproduction order.
