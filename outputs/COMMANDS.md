# Output Command Ledger

This is the canonical rerun ledger for the retained output families. Agent tool
calls were not recorded as a literal terminal transcript, so commands were
reconstructed from the current CLI definitions, YAML configs, run summaries,
and handoff documents. For a config-driven run, the YAML file and the resolved
config saved in its run directory are authoritative.

Commands are run from the repository root in the `nasa-mouse` environment:

```bash
conda activate nasa-mouse
export PYTHONPATH=src
```

## Final Selected Runs

These are the model and analysis commands that reproduce the outputs selected
for the current generative analysis.

### ARCHS4 DDIM backbone

```bash
CFG=configs/generative/diffusion/archs4_mouse_paper_parity_osdr_disjoint.yaml
python -m nasa_mouse_diffusion.paper_parity prepare --config "$CFG"
python -m nasa_mouse_diffusion.paper_parity train --config "$CFG"
python -m nasa_mouse_diffusion.paper_parity evaluate --config "$CFG"
```

### OSDR preparation and DDIM adaptation

```bash
DATA_CFG=configs/generative/diffusion/osdr_factorized_within_study_replicated_validation_data.yaml
BASE_CFG=configs/generative/diffusion/osdr_factorized_study_lora512_replicated_validation_osdr_disjoint.yaml
FINAL_CFG=configs/generative/diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint.yaml

python -m nasa_mouse_diffusion.paper_parity prepare-osdr --config "$DATA_CFG"
python -m nasa_mouse_diffusion.paper_parity train-adapter --config "$BASE_CFG"
python -m nasa_mouse_diffusion.paper_parity train-adapter --config "$FINAL_CFG"

python -m nasa_mouse_diffusion.paper_parity evaluate-adapter \
  --config "$FINAL_CFG" --guidance-scales 1 \
  --validation-sampling-seed 3020 --train-sampling-seed 4020
python -m nasa_mouse_diffusion.paper_parity evaluate-adapter \
  --config "$FINAL_CFG" --guidance-scales 1 \
  --validation-sampling-seed 3022 --train-sampling-seed 4022 \
  --evaluation-variant seed3022
python -m nasa_mouse_diffusion.paper_parity evaluate-adapter \
  --config "$FINAL_CFG" --guidance-scales 1 \
  --validation-sampling-seed 3023 --train-sampling-seed 4023 \
  --evaluation-variant seed3023

python -m nasa_mouse_diffusion.paper_parity calibrate-distribution-adapter \
  --config "$FINAL_CFG" --guidance-scale 1 \
  --prior-strength 5 --residual-scale 0.5 --residual-seed 9100
python -m nasa_mouse_diffusion.paper_parity plot-adapter-trajectory \
  --config "$FINAL_CFG" --snapshot-timesteps 1000 200 0
```

The one-time locked-test output already present under `final_locked_test/` was
produced with the selected calibrator:

```bash
python -m nasa_mouse_diffusion.paper_parity evaluate-finalist-test \
  --config "$FINAL_CFG" \
  --calibrator-dir outputs/generative/benchmark/runs/lacan_diffusion/osdr_factorized_study_lora512_correlation_refine_osdr_disjoint_seed2020/evaluation/repeated_distribution_calibration/prior_5_residual_0.5 \
  --sampling-seeds 5020 5021 5022 5023 --residual-seed 15020 \
  --unlock-test
```

### WGAN comparison

```bash
WGAN_CFG=configs/generative/wgan/wgan_matched_study_conditioned.yaml
python -m nasa_mouse_wgan.matched_study train --config "$WGAN_CFG"
python -m nasa_mouse_wgan.matched_study evaluate-validation --config "$WGAN_CFG"
python -m nasa_mouse_wgan.matched_study screen-calibration --config "$WGAN_CFG"
```

The WGAN run is retained as the comparator. Its locked test was not selected
for the final synthetic analysis.

### Synthetic analyses

```bash
python -m nasa_mouse_diffusion.paper_parity.within_study_feature_stability \
  --config configs/generative/diffusion/within_study_generated_feature_stability_osdr_disjoint.yaml
python -m nasa_mouse_diffusion.paper_parity.within_study_feature_stability \
  --config configs/generative/diffusion/within_study_generated_feature_stability_muscle_groups_osdr_disjoint.yaml

OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
NUMEXPR_NUM_THREADS=1 python \
  -m nasa_mouse_diffusion.paper_parity.matched_all_gene_classifiers \
  --config configs/generative/diffusion/matched_all_gene_classifiers_osdr_disjoint.yaml

python -m nasa_mouse_diffusion.paper_parity.classifier_importance \
  --config configs/generative/diffusion/classifier_importance_osdr_disjoint.yaml
python -m nasa_mouse_diffusion.paper_parity.grouped_pathway_importance \
  --config configs/generative/diffusion/grouped_pathway_importance_osdr_disjoint.yaml
python -m nasa_mouse_diffusion.paper_parity.annotate_importance_literature --check
```

## OSDR API and ARCHS4 Inventories

```bash
python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics
python -m nasa_mouse_generative osdr-inventory --refresh
python -m nasa_mouse_generative osdr-expression --download-missing
python -m nasa_mouse_generative archs4-catalog
python -m nasa_mouse_generative split-plan
python -m nasa_mouse_generative experiment-plan
```

Current API-derived inventory outputs are under
`outputs/generative/benchmark/data_audit/osdr/`. The separate
`outputs/osdr/inventory/legacy_bulk_tissue_counts/` tables are historical and
are not regenerated by the current pipeline.

## GLARE and MOBER

The retained GLARE stage was exploratory and has no output directory designated
as final. These commands reproduce its main API, batch-correction, and
per-study workflow; detailed optional analysis commands remain in
`src/nasa_mouse_glare/README.md`.

```bash
python -m nasa_mouse_glare.multi_tissue_api_glare audit
python -m nasa_mouse_glare.multi_tissue_api_glare prepare \
  --tissues liver skeletal_muscle skin kidney thymus spleen lung retina
python -m nasa_mouse_glare.multi_tissue_api_glare run-glare-scope \
  --tissues liver skeletal_muscle skin kidney thymus spleen lung retina
PYTHONPATH=src:assets/model_sources/MOBER \
python -m nasa_mouse_glare.multi_tissue_api_glare run-mober-scope \
  --tissues liver skeletal_muscle skin kidney thymus spleen lung retina
python -m nasa_mouse_glare.multi_tissue_api_glare run-per-study-glare
python -m nasa_mouse_glare.multi_tissue_api_glare run-dgea-comparison
python -m nasa_mouse_glare.multi_tissue_validation
```

These commands write future results below `outputs/glare/`.

## expiMap Runs

### API-derived OSDR-only baselines

```bash
tissues=(liver skeletal_muscle skin kidney thymus spleen lung retina)
for tissue in "${tissues[@]}"; do
  python -m nasa_mouse_expimap.prepare_expimap_osdr_tissue --tissue "$tissue"
  python -m nasa_mouse_expimap.train_expimap_direct \
    --input "outputs/expimap/runs/direct/$tissue/input/osdr_${tissue}_flt_gc_reactome_raw_counts.h5ad" \
    --output-dir "outputs/expimap/runs/direct/$tissue/raw_counts_nb_50epoch" \
    --recon-loss nb --epochs 50 --hidden-layer-sizes 64
  python -m nasa_mouse_expimap.analyze_expimap_pathways \
    --scores "outputs/expimap/runs/direct/$tissue/raw_counts_nb_50epoch/pathway_scores.tsv" \
    --output-dir "outputs/expimap/runs/direct/$tissue/raw_counts_nb_50epoch/analysis"
done
```

CPM and log1p-CPM MSE runs under the same tissue directories are preprocessing
sensitivities. Their exact input, loss, epoch, and output values are recorded in
each `training_summary.json`.

### ARCHS4 reference and OSDR query

The reference-query runs follow this sequence:

```bash
python -m nasa_mouse_expimap.prepare_expimap_archs4_reference \
  --query-h5ad <osdr_query.h5ad> --tissue <tissue> \
  --output-dir <tissue_run>/reference_input_all --max-samples 0
python -m nasa_mouse_expimap.prepare_expimap_tutorial_hvg \
  --reference-h5ad <reference_input.h5ad> \
  --query-h5ad <osdr_query.h5ad> \
  --output-dir <tissue_run>/tutorial_hvg_2000/input \
  --n-top-genes 2000
python -m nasa_mouse_expimap.train_expimap_archs4_reference \
  --input <tutorial_reference.h5ad> \
  --output-dir <tissue_run>/tutorial_hvg_2000/reference_nb_400epoch_seed2020 \
  --recon-loss nb --epochs 400 --seed 2020
python -m nasa_mouse_expimap.map_expimap_osdr_query \
  --reference-model <reference_model> --query-h5ad <tutorial_query.h5ad> \
  --output-dir <tissue_run>/tutorial_hvg_2000/query_nb_250epoch_seed2020 \
  --epochs 250 --seed 2020
```

This sequence was run for liver, kidney, lung, retina, skeletal muscle, skin,
spleen, and thymus, with seeds 2020, 2021, and 2022 where seed sensitivity was
required. Final liver analysis used the primary-deduplicated query. Kidney and
spleen used the corrected `reassessment_hvg_2000` runs. Exact file names and
arguments are in `tutorial_hvg_input_manifest.json`, `training_summary.json`,
and `query_mapping_summary.json` inside each run.

### Query de novo programs

```bash
python -m nasa_mouse_expimap.map_expimap_osdr_query \
  --reference-model <reference_model> --query-h5ad <tutorial_query.h5ad> \
  --output-dir <query_denovo_output> --epochs 250 --seed 2020 \
  --n-de-novo-programs 3 --gamma-ext 1.0 --use-hsic --hsic-one-vs-all
python -m nasa_mouse_expimap.analyze_expimap_pathways \
  --scores <query_denovo_output>/query_pathway_scores.tsv \
  --output-dir <query_denovo_output>/analysis --include-de-novo
python -m nasa_mouse_expimap.summarize_expimap_de_novo \
  --mapped-h5ad <query_denovo_output>/mapped_query_with_scores.h5ad \
  --scores <query_denovo_output>/query_pathway_scores.tsv \
  --comparison <query_denovo_output>/analysis/flt_vs_gc_pathway_comparison.tsv \
  --study-tests <query_denovo_output>/analysis/flight_ground_study_aware_tests.tsv \
  --programs <query_denovo_output>/de_novo_programs.tsv \
  --gene-loadings <query_denovo_output>/de_novo_program_gene_loadings.tsv \
  --output-dir <query_denovo_output>/de_novo_analysis
```

These outputs are sensitivity and hypothesis-generation runs, not the primary
pathway models.

### Muscle groups and cross-run analyses

```bash
python -m nasa_mouse_expimap.split_expimap_muscle_groups
python -m nasa_mouse_expimap.analyze_muscle_targeted_modules
python -m nasa_mouse_expimap.run_expimap_tissue_variant_matrix --execute
python -m nasa_mouse_expimap.plot_expimap_pathway_followup \
  --tissue thymus --tissue skin --tissue spleen
python -m nasa_mouse_expimap.run_kidney_spleen_seed_sensitivity
python -m nasa_mouse_expimap.analyze_kidney_spleen_reassessment
python -m nasa_mouse_expimap.curate_kidney_spleen_reassessment
python -m nasa_mouse_expimap.summarize_expimap_results
```

## Generative Benchmark Screens

The configurable benchmark and its parameter matrix were launched with:

```bash
python -m nasa_mouse_generative prepare-upstreams
python -m nasa_mouse_generative experiment-plan
python -m nasa_mouse_generative matrix-run --dry-run --max-runs 0
python -m nasa_mouse_generative matrix-run --max-runs 0
python -m nasa_mouse_generative scoreboard
```

Every matrix row is recorded in
`outputs/generative/benchmark/summary/experiment_plan.tsv`. DDIM development,
harmonization, tissue conditioning, study conditioning, material-type ablation,
and whole-study transfer runs use the YAML files under
`configs/generative/diffusion/`. The command form is determined by each config's
`contract`:

```bash
# ARCHS4 reference configs
python -m nasa_mouse_diffusion.paper_parity prepare --config <config.yaml>
python -m nasa_mouse_diffusion.paper_parity train --config <config.yaml>
python -m nasa_mouse_diffusion.paper_parity evaluate --config <config.yaml>

# OSDR extension configs
python -m nasa_mouse_diffusion.paper_parity prepare-osdr --config <config.yaml>
python -m nasa_mouse_diffusion.paper_parity train-osdr --config <config.yaml>
python -m nasa_mouse_diffusion.paper_parity evaluate-osdr --config <config.yaml>

# Factorized adapter configs
python -m nasa_mouse_diffusion.paper_parity train-adapter --config <config.yaml>
python -m nasa_mouse_diffusion.paper_parity evaluate-adapter --config <config.yaml>
```

Run directories contain the resolved config and summaries needed to identify
the exact command variant. Smoke runs are temporary and are not retained.
Obsolete harmonization screens, earlier LoRA ranks, no-material ablations, and
alternate seeds are not final selections.

## Standalone Generative Runs

These commands produced the earlier model-specific trees now grouped under
`outputs/generative/standalone/`:

```bash
python -m nasa_mouse_diffusion.run_conditional_generation
python -m nasa_mouse_diffusion.generate_synthetic_examples --overwrite
python -m nasa_mouse_diffusion.summarize_results

python -m nasa_mouse_wgan.run_conditional_generation
python -m nasa_mouse_wgan.generate_synthetic \
  --model-dir <trained_wgan_model_dir> --output-dir <generation_output_dir> \
  --n 1000 --condition flight
python -m nasa_mouse_wgan.run_pipeline \
  --include-direct --include-reference --include-muscle-splits
python -m nasa_mouse_wgan.summarize_results
```

These runs were superseded by the unified benchmark and are not used in the
final analysis.

## Manuscripts, Figures, and Presentation

```bash
python -m nasa_mouse_expimap.build_asgsr_paper
python -m nasa_mouse_expimap.build_publication_figures
python -m nasa_mouse_diffusion.paper_parity.build_synthetic_guided_paper
python -m nasa_mouse_diffusion.paper_parity.build_slstp_presentation
python -m nasa_mouse_internship_report.build_report
```

Publication builders read only the final selections listed in
`outputs/README.md`; they do not choose a model by scanning every run directory.
