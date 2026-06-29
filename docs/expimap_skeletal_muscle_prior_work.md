# Skeletal Muscle expiMap Prior-Work Check

Date: 2026-06-26

## Bottom Line

Yes. Skeletal muscle is expected to be one of the strongest spaceflight-affected systems, especially unloaded/postural muscle. Prior work points to muscle atrophy, contractile remodeling, calcium handling, mitochondrial/metabolic stress, autophagy/proteostasis, extracellular matrix remodeling, and insulin/IGF/AKT-related metabolism.

The current aggregate skeletal-muscle expiMap result does not yet give a clean Reactome-level FDR-significant muscle pathway call in the preferred all-eligible ARCHS4 reference -> OSDR query run. The expected biology is present in the ranked pathway list, but mostly below FDR significance and with accession/muscle-type heterogeneity.

## Local expiMap Result

Primary local comparison:

- OSDR query samples: 191 skeletal-muscle FLT/GC samples.
- Accessions: 13.
- ARCHS4 all-eligible skeletal-muscle reference: 1,412 samples from 88 series.
- Architecture: 1,140 mouse Reactome pathways, 9,319 matched genes.
- Preferred run: `outputs/expimap_archs4_reference_osdr_query_skeletal_muscle/query_nb_allref_50epoch`.

Summary across skeletal-muscle runs:

| run | Reactome terms | FDR < 0.05 | min FDR | interpretation |
| --- | ---: | ---: | ---: | --- |
| direct OSDR NB 50 epoch | 1,140 | 2 | 0.00538 | Two hits, but not a canonical muscle-atrophy result. |
| ARCHS4 1k reference -> OSDR query | 1,140 | 39 | 0.00000782 | Reference-size-sensitive; not carried forward as final. |
| ARCHS4 all-eligible reference -> OSDR query | 1,140 | 0 | 0.10076 | Preferred run; no FDR-significant Reactome pathway. |

Top preferred-run terms:

| term | meta effect | raw p | FDR | same direction |
| --- | ---: | ---: | ---: | ---: |
| `R-MMU-4086398_CA2_PATHWAY` | 0.00578 | 0.000099 | 0.10076 | 11/13 |
| `R-MMU-9013026_RHOB_GTPASE_CYCLE` | 0.01192 | 0.000177 | 0.10076 | 12/13 |
| `R-MMU-6804115_TP53_REGULATES_TRANSCRIPTION_OF_ADDITIONAL_CELL_CYCLE_GENES_WHOSE_EXACT_ROLE_IN_THE_P53_PATHWAY_REMAIN_UNCERTAIN` | 0.00276 | 0.000537 | 0.10288 | 12/13 |
| `R-MMU-5654228_PHOSPHOLIPASE_C_MEDIATED_CASCADE_FGFR4` | -0.00347 | 0.000590 | 0.10288 | 11/13 |
| `R-MMU-936440_NEGATIVE_REGULATORS_OF_DDX58_IFIH1_SIGNALING` | -0.00982 | 0.000655 | 0.10288 | 11/13 |

Expected muscle themes in the preferred-run ranking:

| theme | representative terms | FDR range | note |
| --- | --- | ---: | --- |
| Calcium/contractile handling | `CA2_PATHWAY`, `PLATELET_CALCIUM_HOMEOSTASIS`, `MITOCHONDRIAL_CALCIUM_ION_TRANSPORT`, `STRIATED_MUSCLE_CONTRACTION` | 0.1008-0.3626 | Calcium signal is near the top; striated contraction is heterogeneous. |
| Proteostasis/autophagy | `AUTOPHAGY`, `MACROAUTOPHAGY`, `FOXO_MEDIATED_TRANSCRIPTION`, ubiquitin degradation terms, lysosome terms | 0.1347-0.1903 for the better-ranked terms | Present but not FDR-significant. |
| Mitochondria/metabolism | mitochondrial fatty-acid beta oxidation, TCA cycle, mitochondrial quality control | 0.1347-0.1978 | Direction varies by accession/muscle. |
| IGF/insulin/AKT axis | IGF1R/IRS and IGFBP terms | 0.1347-0.1803 | IGF terms are high-ranked; mTOR terms are weak in this Reactome run. |
| ECM/fibrosis | collagen fibril assembly, collagen biosynthesis, collagen formation | 0.1347-0.2063 | Similar to known remodeling biology, not significant here. |

The follow-up plots and leave-one-accession-out checks for the preferred muscle run are in:

- `outputs/expimap_pathway_followup_muscle_allref/`
- `outputs/expimap_pathway_followup_muscle_allref/skeletal_muscle/plots/`

None of the selected preferred-run muscle pathways pass the stricter follow-up rule: FDR < 0.05, same direction across all accessions, same direction under every leave-one-accession-out refit, and max leave-one-out FDR < 0.05.

## Why The Aggregate expiMap Result Can Miss Expected Muscle Biology

The current `skeletal_muscle` label aggregates different muscles and missions:

| material type | flight | ground control |
| --- | ---: | ---: |
| Extensor digitorum longus, both sides | 6 | 6 |
| Gastrocnemius | 3 | 3 |
| Left gastrocnemius | 6 | 6 |
| Left quadriceps femoris | 6 | 6 |
| Left tibialis anterior | 6 | 6 |
| Quadriceps femoris | 7 | 7 |
| Right extensor digitorum longus | 10 | 10 |
| Right gastrocnemius | 4 | 8 |
| Right quadriceps femoris | 10 | 10 |
| Right soleus | 10 | 10 |
| Right tibialis anterior | 9 | 9 |
| Soleus | 12 | 9 |
| Soleus, both sides | 6 | 6 |

This matters because soleus, EDL, tibialis anterior, gastrocnemius, and quadriceps do not respond identically. Some are slow/oxidative and unloading-sensitive; others are fast/glycolytic or more resilient. A random-effects accession-aware test will intentionally downweight pathways whose direction changes by accession or muscle.

Reactome also does not give a dedicated "spaceflight muscle atrophy" program. The most relevant biology is split across broad terms such as calcium, autophagy, ubiquitin, FOXO, collagen, TCA, mitochondrial transport, and IGF/insulin signaling. A targeted muscle gene-set analysis may be more sensitive than a broad Reactome mask.

## Prior Work Alignment

The literature supports skeletal muscle as a major spaceflight target:

- Vitry et al. 2022, iScience, analyzed NASA RR1 liver and quadriceps after 37 days of spaceflight and reported that muscle and liver are among the most affected tissues; they linked impaired liver lipid metabolism with muscle atrophy gene expression, including muscle autophagy/translation/DNA-repair processes. Link: https://researchonline.ljmu.ac.uk/id/eprint/18458/
- Oommen et al. 2024, npj Microgravity, profiled mouse muscle transcriptomes using GeneLab data and highlighted muscle atrophy biology including growth-factor signaling, proteolysis, oxidative stress, catabolism, mitochondrial dysfunction, ubiquitination, autophagy, ion transport, ECM interactions, and muscle-specific markers. Link: https://www.nature.com/articles/s41526-024-00434-z
- NASA/data.gov OSD-326 describes RR quadriceps spaceflight data where flight quadriceps wet weight was reduced, myosin/troponin genes were suppressed, and networks linked to sarcomeric integrity, immune fitness, oxidative stress, metabolism, and ATP synthesis/hydrolysis were inhibited. Link: https://catalog.data.gov/dataset/gene-metabolite-network-linked-to-inhibited-bioenergetics-in-association-with-spaceflight-
- Li et al. 2023, npj Microgravity, used NASA OSDR multi-omics data from soleus and tibialis anterior during 30+ day ISS missions and focused on muscle atrophy, calcium dysregulation, SERCA pump function, and calcium reuptake biomarkers. Link: https://www.nature.com/articles/s41526-023-00337-5
- da Silveira et al. 2020, Cell, found mitochondrial stress as a central spaceflight phenotype across GeneLab and astronaut data, with mitochondrial activity, innate immunity, chronic inflammation, cell cycle, and circadian rhythm among enriched pathways. Link: https://doi.org/10.1016/j.cell.2020.11.002
- Beheshti et al. 2023, Communications Biology, found the liver had the largest DEG count in the RR1 cross-tissue analysis, but soleus was second and EDL was also highly affected; the paper explicitly notes muscle-specific insulin/estrogen/metabolic pathway variation between soleus and EDL. Link: https://www.nature.com/articles/s42003-023-05213-2

## Interpretation

Muscle is still a strong biological candidate. The current result says the aggregate Reactome expiMap model is not yet capturing a stable, cross-accession skeletal-muscle pathway at FDR < 0.05 in the preferred all-reference run.

The most plausible explanation is not "muscle has no spaceflight effect." It is that the effect is muscle-type-specific, gene/module-specific, and heterogeneous across accessions. That is exactly where the prior literature points: soleus and EDL can show very different responses, and contractile/calcium/mitochondrial/proteostasis biology can appear as marker-gene or module-level changes rather than one broad Reactome pathway.

## Targeted Module Follow-up

Completed follow-up output:

- `outputs/osdr_skeletal_muscle_targeted_modules/`
- `outputs/osdr_skeletal_muscle_targeted_modules/README.md`
- `outputs/osdr_skeletal_muscle_targeted_modules/plots/targeted_module_effect_heatmap.png`
- `outputs/osdr_skeletal_muscle_targeted_modules/plots/top_targeted_module_boxplots.png`
- `data/pathways/mouse_muscle_targeted_modules.gmt`
- `data/pathways/reactome_current_mouse_ensembl_plus_muscle_targeted.gmt`

This analysis reconstructs the full downloaded OSDR unnormalized count matrix for the 191 skeletal-muscle FLT/GC samples, maps 140 targeted muscle gene symbols to mouse Ensembl IDs, scores 12 targeted modules as mean z-scored log1p(CPM), and runs accession-aware random-effects tests globally and by material-type-derived muscle group.

It also writes a 12-term targeted muscle GMT and a combined 1,848-term Reactome+targeted GMT for follow-up expiMap retraining.

Per-muscle sample sizes are small, so these are follow-up candidates rather than final calls: EDL and tibialis anterior have two accessions, soleus and gastrocnemius have three, and quadriceps has four. The `strict_candidate` flag requires FDR < 0.05, at least three accessions, and all accession effects in the same direction. The stronger `strict_loo_candidate` flag additionally requires every leave-one-accession-out refit to keep FDR < 0.05 and the same direction.

Targeted module candidates that pass FDR/direction before the LOO filter:

| analysis group | module | effect | FDR | accessions |
| --- | --- | ---: | ---: | ---: |
| gastrocnemius | `contractile_slow_postural` | -0.2255 | 4.99e-14 | 3/3 |
| soleus | `contractile_fast_sarcomere` | 0.4182 | 3.05e-4 | 3/3 |
| soleus | `fatty_acid_oxidation` | -0.7325 | 0.0270 | 3/3 |
| soleus | `myogenesis_regeneration` | 0.4900 | 0.0270 | 3/3 |
| soleus | `contractile_slow_postural` | -0.1328 | 0.0441 | 3/3 |
| quadriceps | `contractile_fast_sarcomere` | -0.0826 | 0.0453 | 4/4 |

None of these targeted modules pass the stricter LOO FDR stability rule. Their leave-one-out effects keep the same direction, but the maximum leave-one-out FDR ranges from 0.125 to 0.582. This supports muscle-type-specific follow-up, not a final robust module claim.

Top strict gene-level candidates align with known spaceflight muscle biology:

- Gastrocnemius: lower `Myl3`, `Atp2a2`, `Myl2`, `Tnni1`, and `Tnc` in flight.
- Soleus: higher `Myh1`, `Cacna1s`, `Ckm`, `Camk2d`, `Musk`, `Myod1`, and `Myf5`; lower `Tnnt1`, `Myl2`, `Myl3`, `Myh7`, `Atp2a2`, `Pln`, `Slc25a20`, `Acadm`, `Hadha`, and `Hadhb`.
- Quadriceps: targeted module effects include lower fast sarcomere and mitochondrial/OXPHOS/TCA trends, plus IGF/AKT/mTOR and regeneration signals, but some are not strict because not all four accessions match direction.

The stratified expiMap Reactome score analysis also finds strict candidates, mainly in soleus and gastrocnemius, but it repeats 1,140 pathway tests inside small muscle-type subsets. Those Reactome hits should get leave-one-accession-out checks before being called.

## Targeted expiMap Retraining: Direct vs Reference-Query

The earlier targeted-module result above was not a full expiMap retraining in both modes. It was a count-level targeted module score analysis, plus a stratified look at existing aggregate expiMap Reactome scores. The follow-up below reruns expiMap with the combined Reactome+targeted muscle architecture and separates the adequate muscle groups before testing FLT vs GC.

New follow-up output:

- `outputs/expimap_muscle_targeted_combined_min8/`
- `outputs/expimap_muscle_targeted_combined_min8/summary/README.md`
- `outputs/expimap_muscle_targeted_combined_min8/summary/run_validation_summary.tsv`
- `outputs/expimap_muscle_targeted_combined_min8/summary/targeted_muscle_module_validation_summary.tsv`

Setup:

- Architecture: 1,359 terms = official mouse Reactome terms plus 12 targeted muscle modules.
- Direct OSDR split runs: gastrocnemius, quadriceps, and soleus, each trained for 100 epochs.
- ARCHS4 reference-query runs: 1,412 ARCHS4 skeletal-muscle reference samples trained for 100 epochs, then each OSDR muscle group mapped for 50 epochs.
- Hardware: all new direct, reference, and query runs used `NVIDIA A100-SXM4-40GB` with trainer device `cuda`.
- Query mapping dropped two query-only genes because the ARCHS4 reference contained 9,342 of the 9,344 query genes.

Run-level validation summary:

| mode | group | samples | all-term FDR < 0.05 | all-term LOO pass | targeted FDR < 0.05 | targeted LOO pass |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| direct OSDR split | gastrocnemius | 30 | 1 | 0 | 0 | 0 |
| direct OSDR split | quadriceps | 46 | 23 | 0 | 1 | 0 |
| direct OSDR split | soleus | 53 | 6 | 1 | 0 | 0 |
| ARCHS4 reference -> OSDR query split | gastrocnemius | 30 | 190 | 0 | 2 | 0 |
| ARCHS4 reference -> OSDR query split | quadriceps | 46 | 4 | 0 | 0 | 0 |
| ARCHS4 reference -> OSDR query split | soleus | 53 | 245 | 89 | 1 | 0 |

Targeted module hits with meta FDR < 0.05:

| mode | group | module | effect | FDR | same-direction accessions | max LOO FDR |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| ARCHS4 reference -> OSDR query split | gastrocnemius | `atrophy_ubiquitin_foxo` | 0.2086 | 0.0225 | 3/3 | 0.779 |
| ARCHS4 reference -> OSDR query split | gastrocnemius | `contractile_slow_postural` | -0.2439 | 0.0476 | 3/3 | 0.779 |
| ARCHS4 reference -> OSDR query split | soleus | `ecm_fibrosis_remodeling` | -0.3605 | 4.33e-5 | 3/3 | 0.0569 |
| direct OSDR split | quadriceps | `inflammation_oxidative_stress` | 0.000130 | 0.0366 | 4/4 | 0.775 |

None of the targeted modules pass the stricter leave-one-accession-out FDR rule in either direct or ARCHS4-reference query mode. The directions are stable, but at least one held-out-accession refit loses FDR < 0.05.

The failures are accession-leverage failures, not sign flips:

| mode | group | module | weakest held-out accession | LOO effect | LOO FDR |
| --- | --- | --- | --- | ---: | ---: |
| ARCHS4 reference -> OSDR query split | gastrocnemius | `atrophy_ubiquitin_foxo` | OSD-401 | 0.1238 | 0.779 |
| ARCHS4 reference -> OSDR query split | gastrocnemius | `contractile_slow_postural` | OSD-401 | -0.1374 | 0.779 |
| ARCHS4 reference -> OSDR query split | soleus | `ecm_fibrosis_remodeling` | OSD-104 | -0.2619 | 0.0569 |
| direct OSDR split | quadriceps | `inflammation_oxidative_stress` | OSD-666 | 0.000116 | 0.775 |

This explains the apparent mismatch with prior literature. The expected muscle biology is not absent: atrophy/ubiquitin/FOXO, contractile slow/postural, ECM remodeling, inflammation/oxidative-stress, IGF, calcium, mitochondrial, and sarcomere signals appear in the targeted count-level results or the split expiMap rankings. What fails is the stricter claim that the same module remains FDR-significant after removing any one accession from a per-muscle analysis with only three or four accessions. Most prior studies did not require cross-accession leave-one-study-out FDR stability; they often analyzed one mission/muscle at a time, or combined transcriptomics with wet weight, metabolomics, calcium physiology, or other endpoints.

The ARCHS4-reference query split is more sensitive than direct OSDR training for some muscles, especially soleus, but it also produces many broad Reactome hits. Those all-term hits should be interpreted separately from the targeted modules because the reference model can amplify stable latent shifts that are not necessarily canonical spaceflight muscle-atrophy modules.

### latent_enrich Bayes-Factor Check

The other expiMap paper-style method is `latent_enrich`, which reports a Bayes-factor-style score `log(p_h0 / p_h1)` for condition-level latent-score separation. I ran it on the same split-muscle direct and ARCHS4-reference query models with `condition_inferred` as the grouping variable, `ground_control` as comparison, direction-corrected latent scores, exact probabilities, and 5,000 samples.

| mode | group | max abs(BF) | terms abs(BF) >= 0.5 | terms abs(BF) >= 1.0 | terms abs(BF) >= 2.3 |
| --- | --- | ---: | ---: | ---: | ---: |
| direct OSDR split | gastrocnemius | 0.0351 | 0 | 0 | 0 |
| direct OSDR split | quadriceps | 0.000577 | 0 | 0 | 0 |
| direct OSDR split | soleus | 0.00196 | 0 | 0 | 0 |
| ARCHS4 reference -> OSDR query split | gastrocnemius | 0.684 | 27 | 0 | 0 |
| ARCHS4 reference -> OSDR query split | quadriceps | 0.205 | 0 | 0 | 0 |
| ARCHS4 reference -> OSDR query split | soleus | 1.371 | 62 | 4 | 0 |

No split-muscle run passes the local paper-style cutoff of `abs(BF) >= 2.3`. Direct OSDR split models are essentially flat by this method. The ARCHS4-reference query soleus run has moderate all-term BF movement, with top broad terms including amino-acid metabolism, immune system, olfactory signaling, and post-translational protein modification, but still no paper-style BF call. Targeted muscle modules are weaker: the largest targeted values are 0.525 for gastrocnemius `MUSCLE_CONTRACTILE_SLOW_POSTURAL`, 0.536 for soleus `MUSCLE_CONTRACTILE_SLOW_POSTURAL`, and 0.134 for quadriceps `MUSCLE_CONTRACTILE_FAST_SARCOMERE` in the ARCHS4-reference query runs.

This reinforces the same interpretation as the accession-aware validation: the split-muscle results contain follow-up signal, especially in reference-query soleus and gastrocnemius, but the paper-style BF method does not support a final strong condition-level module call.

### Exploratory EDL and Tibialis Anterior Runs

EDL and tibialis anterior were originally skipped in the split-muscle expiMap
retraining because each has only two paired accessions. I trained them anyway
as exploratory follow-up, keeping the same combined Reactome+targeted
architecture, NB loss, direct 100-epoch runs, ARCHS4-reference query 50-epoch
mapping, accession-aware validation, plots, and `latent_enrich` BF checks.

New exploratory output:

- `outputs/expimap_muscle_targeted_combined_min8/group_inputs_exploratory_2acc/`
- `outputs/expimap_muscle_targeted_combined_min8/direct_edl_nb_100epoch/`
- `outputs/expimap_muscle_targeted_combined_min8/direct_tibialis_anterior_nb_100epoch/`
- `outputs/expimap_muscle_targeted_combined_min8/query_edl_nb_allref_50epoch/`
- `outputs/expimap_muscle_targeted_combined_min8/query_tibialis_anterior_nb_allref_50epoch/`
- `outputs/expimap_muscle_targeted_combined_min8/summary/exploratory_two_accession_muscle_split_summary.tsv`

Run-level exploratory summary:

| mode | group | samples | accessions | all-term FDR < 0.05 | all-term LOO pass | targeted FDR < 0.05 | targeted LOO pass | max abs(BF) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| direct OSDR split | EDL | 32 | 2 | 9 | 0 | 0 | 0 | 0.00174 |
| direct OSDR split | tibialis anterior | 30 | 2 | 56 | 0 | 1 | 0 | 0.000568 |
| ARCHS4 reference -> OSDR query split | EDL | 32 | 2 | 1 | 0 | 0 | 0 | 0.329 |
| ARCHS4 reference -> OSDR query split | tibialis anterior | 30 | 2 | 0 | 0 | 0 | 0 | 0.265 |

The direct tibialis anterior run has one targeted module with ordinary
accession-aware meta FDR < 0.05: `MUSCLE_IGF_AKT_MTOR` is lower in flight
(effect -0.000201, FDR 0.00478, 2/2 accessions same direction). It does not
pass leave-one-accession-out stability; each single-accession refit has FDR
1.0. The broad all-term hits in direct EDL and direct tibialis also fail LOO.
No exploratory EDL/tibialis run passes the local paper-style `abs(BF) >= 2.3`
criterion.

Interpretation: these runs are useful sanity checks and point to tibialis
IGF/AKT/mTOR as a candidate worth inspecting at gene/sample level, but they do
not change the main conclusion. The robust split-muscle evidence still rests
more on the three- and four-accession groups, especially soleus,
gastrocnemius, and quadriceps.

I then ran the missing tutorial-style variants for EDL and tibialis anterior:
all-gene reference-query with 3 de novo extension programs, HVG reference-query,
and HVG reference-query with 3 de novo extension programs. These used the same
settings as the full tissue variant matrix: 250 query epochs, `alpha_kl=0.22`,
no group-lasso alpha in query mapping, 3 de novo programs, `gamma_ext=0.7`,
HSIC one-vs-all regularization, and 2,000 reference-selected HVGs. HVG filtering
retained 1,975 genes and 357 annotated terms.

Variant summary:

| variant | group | all-term FDR < 0.05 | targeted FDR < 0.05 | LOO pass | de novo condition FDR < 0.05 | max abs(BF) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| reference-query + de novo | EDL | 1 | 0 | 0 | 0 | 0.331 |
| reference-query + de novo | tibialis anterior | 0 | 0 | 0 | 0 | 0.264 |
| HVG reference-query | EDL | 8 | 0 | 0 | 0 | 0.333 |
| HVG reference-query | tibialis anterior | 0 | 0 | 0 | 0 | 0.650 |
| HVG reference-query + de novo | EDL | 8 | 0 | 0 | 0 | 0.333 |
| HVG reference-query + de novo | tibialis anterior | 0 | 0 | 0 | 0 | 0.652 |

Output tables:

- `outputs/expimap_muscle_targeted_combined_min8/summary/exploratory_two_accession_variant_matrix_summary.tsv`
- `outputs/expimap_muscle_targeted_combined_min8/summary/exploratory_two_accession_de_novo_program_summary.tsv`

These runs explain why the ARCHS4 reference did not provide a stronger result
for EDL/tibialis. The reference-query model maps two-accession query data into
an aggregate ARCHS4 skeletal-muscle latent space. That can stabilize broad
program scores, but it also limits query-specific condition movement unless the
FLT-vs-GC shift is strong and consistent across both accessions. The de novo
extension nodes did learn weighted gene programs, but none separated FLT vs GC
at aggregate Welch FDR < 0.05 or study-aware FDR < 0.05. HVG selection made EDL
somewhat more sensitive for broad annotated terms, but it still did not produce
a targeted muscle-module hit, a LOO-stable term, or a paper-style BF call. So
the reference-query negative result is not just because we forgot HVG/de novo;
the signal remains weak under those tutorial-style variants.

## Recommended Next Checks

1. For the ARCHS4-reference query soleus hits that pass all-term LOO, inspect pathway loadings and sample-level plots before calling biology.
2. Run gene-level accession-aware meta-analysis for the broader muscle transcriptome, not only the targeted module genes.
3. Compare FLT and GC clusters within each muscle type, not only across the merged `skeletal_muscle` label.
4. Treat targeted modules that fail LOO as follow-up candidates and check which held-out accession drives each failure.
