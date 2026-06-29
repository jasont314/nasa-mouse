# Larger ARCHS4 Liver Reference: Seed Stability

## Design

The bounded 1,000-sample reference was expanded to a 5,000-sample mouse liver
reference. The input is a proportional series-stratified sample of 8,970
eligible nonleakage ARCHS4 samples and retains all 518 contributing ARCHS4
series. It has 9,319 shared Ensembl genes and 1,140 current Reactome mouse
terms.

Three raw-count NB expiMap reference models were trained on GPU with a 64-unit
hidden layer, 90/10 train-validation split, batch size 128, 200 maximum epochs,
and early stopping. They stopped after 166, 138, and 93 epochs (best epochs
134, 106, and 61) for seeds 2020, 2021, and 2022. Each was then used to map the
same 231-sample OSDR liver query for 50 epochs. Results use posterior-mean
query scores and the same accession-aware random-effects plus
leave-one-accession-out validation as the direct run.

## Result

Only three terms were FDR-significant with the same effect direction in all
three seeds. One also passed every leave-one-accession-out test:
`R-MMU-416700_OTHER_SEMAPHORIN_INTERACTIONS`. Its FLT-minus-ground model-score
effects were -0.160, -0.171, and -0.152 with FDRs 0.00716, 0.000564, and
0.00728. Eleven of 12 accession effects had that direction in every seed.

This remains a candidate, not a validated pathway result. Its between-accession
heterogeneity is substantial (I2 0.73, 0.68, and 0.63 across the three seeds),
and the all-term effect-rank Spearman correlations between seeds are near zero
(-0.009 to 0.010). The model is therefore seed-sensitive overall.

`R-MMU-2172127_DAP12_INTERACTIONS` and
`R-MMU-9943965_CHD3_CHD4_CHD5_SUBFAMILY` were significant in every seed but
failed the strict leave-one-accession-out gate. Polymerase II transcription
elongation was not significant in any seed and its effect direction changed
between seeds, so the larger reference does not support it.

## Output Artifacts

- `outputs/expimap_archs4_reference_osdr_query_liver/reference_input_5000_stratified/reference_input_manifest.json`
- `outputs/expimap_archs4_reference_osdr_query_liver/reference_seed_stability_5000/README.md`
- `outputs/expimap_archs4_reference_osdr_query_liver/reference_seed_stability_5000/reference_seed_pathway_stability.tsv`
- `outputs/expimap_archs4_reference_osdr_query_liver/query_nb_5000stratified_seed{2020,2021,2022}_50epoch/posterior_mean_accession_validation/`

The preparation and cross-seed summary code is
`src/nasa_mouse_glare/prepare_expimap_archs4_reference.py` and
`src/nasa_mouse_glare/summarize_expimap_reference_seeds.py`.

## Inference Gate

Use the stable semaphorin result to prioritize a count-level, study-aware
follow-up rather than as a final annotation. Replicate it with another reference
sampling strategy or an independent mouse liver cohort before treating it as a
spaceflight-associated pathway.
