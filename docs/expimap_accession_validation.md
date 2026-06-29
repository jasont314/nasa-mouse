# Accession-Aware expiMap Validation

## Purpose

The pooled FLT-versus-ground comparison can be driven by differences between
OSDR accessions. This validation computes a FLT-minus-ground effect within each
accession, combines the 12 liver-accession effects with a random-effects
meta-analysis, and repeats the analysis while leaving each accession out once.

Scores are posterior means exported from the saved direct raw-count NB model,
not a stochastic latent draw. The input contains 231 liver samples (118 FLT,
113 ground control), 9,321 genes, and 1,140 Reactome terms.

## Current Result

The earlier pooled result for
`R-MMU-75955_RNA_POLYMERASE_II_TRANSCRIPTION_ELONGATION` does not pass this
validation gate. Its posterior-mean random-effects estimate is +0.000286,
with P = 0.574, FDR = 0.670, and I2 = 0.904; seven of 12 accession effects
have the meta-analysis direction and five oppose it. No leave-one-accession-out
analysis makes it FDR-significant.

The old direct summary used a stochastic latent sample, which can make a small
pooled difference look more definite than it is. It is retained as a historical
workflow result, not evidence that Polymerase II elongation changes in flight.

The posterior-mean direct meta-analysis has 381 terms at FDR < 0.05. This is a
model-score result on one trained model, not a list of validated liver biology:
the effects are small on the latent-score scale and have not yet been replicated
across direct-model seeds or a count-level differential-expression model.

## Output Artifacts

- `outputs/expimap_direct_osdr_liver/raw_counts_nb_50epoch/posterior_mean_pathway_scores.tsv`
- `outputs/expimap_direct_osdr_liver/raw_counts_nb_50epoch/accession_validation/random_effects_meta_analysis.tsv`
- `outputs/expimap_direct_osdr_liver/raw_counts_nb_50epoch/accession_validation/per_accession_effects.tsv`
- `outputs/expimap_direct_osdr_liver/raw_counts_nb_50epoch/accession_validation/leave_one_accession_out.tsv`

The implementation is `src/nasa_mouse_glare/export_expimap_scores.py` and
`src/nasa_mouse_glare/validate_expimap_accession_effects.py`.

## Inference Gate

Treat a pathway as a candidate only when its direction is stable across direct
model seeds, accession-aware analyses, and a count-level model that accounts
for accession or study. Test it in an independent OSDR subset before making a
biological claim.
