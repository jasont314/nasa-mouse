# WGAN Results Summary

Generated from local WGAN outputs on 2026-06-29.

Primary machine-readable outputs:

- `outputs/wgan_pipeline/summary/wgan_run_summary.tsv`
- `outputs/wgan_pipeline/summary/wgan_top_features.tsv`

Each WGAN run writes PCA, UMAP, accession-colored plots, and a top-feature heatmap under its `analysis/` directory. ARCHS4-pretrained runs also write the same plots for the frozen pretrained OSDR projection under `pretrained_query_analysis/`.

## Interpretation Rules

WGAN uses learned critic features, not Reactome pathways. A WGAN feature can show FLT-vs-GC separation, but it is not itself a named pathway or gene module. Treat WGAN as complementary representation evidence unless followed by feature attribution, DGEA, or pathway enrichment.

Strict LOO-stable means random-effects FDR < 0.05, every leave-one-accession-out FDR < 0.05, and all LOO effects keep the same direction.

## Main Tissue Results

| tissue | mode | score set | query | ARCHS4 ref | RE FDR hits | LOO-stable hits | min RE FDR | top feature | effect | top max LOO FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| spleen | archs4_pretrain_osdr_finetune | pre_finetune_projection | 109 | 6289 | 178 | 137 | 5.79e-35 | WGAN_FEATURE_192 | 2.934 | 4.73e-17 |
| thymus | archs4_pretrain_osdr_finetune | finetuned_or_direct | 117 | 1362 | 177 | 132 | 1.17e-17 | WGAN_FEATURE_046 | 1.774 | 7.71e-05 |
| skeletal_muscle | archs4_pretrain_osdr_finetune | pre_finetune_projection | 191 | 1412 | 141 | 120 | 2.01e-13 | WGAN_FEATURE_231 | -1.851 | 6.01e-08 |
| thymus | archs4_pretrain_osdr_finetune | pre_finetune_projection | 117 | 1362 | 182 | 117 | 3.96e-16 | WGAN_FEATURE_008 | 1.292 | 1.07e-05 |
| spleen | archs4_pretrain_osdr_finetune | finetuned_or_direct | 109 | 6289 | 168 | 107 | 6.10e-31 | WGAN_FEATURE_160 | 2.509 | 2.58e-14 |
| liver | archs4_pretrain_osdr_finetune | pre_finetune_projection | 231 | 5000 | 116 | 77 | 3.62e-22 | WGAN_FEATURE_249 | -1.213 | 1.62e-13 |
| skeletal_muscle | direct_osdr | finetuned_or_direct | 191 | 0 | 105 | 44 | 6.58e-04 | WGAN_FEATURE_091 | 5.352 | 0.00614 |
| skin | archs4_pretrain_osdr_finetune | pre_finetune_projection | 151 | 2593 | 103 | 23 | 5.75e-24 | WGAN_FEATURE_254 | 1.334 | 1.75e-15 |
| liver | archs4_pretrain_osdr_finetune | finetuned_or_direct | 231 | 5000 | 28 | 18 | 4.88e-31 | WGAN_FEATURE_111 | -1.074 | 2.73e-21 |
| liver | direct_osdr | finetuned_or_direct | 231 | 0 | 50 | 4 | 5.69e-05 | WGAN_FEATURE_248 | 3.671 | 0.03993 |
| lung | archs4_pretrain_osdr_finetune | pre_finetune_projection | 78 | 5674 | 19 | 4 | 3.00e-25 | WGAN_FEATURE_134 | 5.576 | 0.001005 |
| skin | archs4_pretrain_osdr_finetune | finetuned_or_direct | 151 | 2593 | 35 | 3 | 3.43e-04 | WGAN_FEATURE_084 | -1.346 | 0.01099 |
| lung | archs4_pretrain_osdr_finetune | finetuned_or_direct | 78 | 5674 | 18 | 3 | 7.15e-18 | WGAN_FEATURE_042 | 2.491 | 0.003922 |
| spleen | direct_osdr | finetuned_or_direct | 109 | 0 | 52 | 0 | 4.67e-09 | WGAN_FEATURE_028 | 2.974 | 0.9654 |
| kidney | archs4_pretrain_osdr_finetune | pre_finetune_projection | 135 | 1000 | 29 | 0 | 0.001872 | WGAN_FEATURE_068 | 0.364 | 0.1276 |
| thymus | direct_osdr | finetuned_or_direct | 117 | 0 | 20 | 0 | 4.13e-12 | WGAN_FEATURE_159 | 3.541 | 0.09074 |
| retina | archs4_pretrain_osdr_finetune | pre_finetune_projection | 76 | 1187 | 9 | 0 | 2.31e-04 | WGAN_FEATURE_041 | 1.017 | 0.3763 |
| lung | direct_osdr | finetuned_or_direct | 78 | 0 | 1 | 0 | 0.001327 | WGAN_FEATURE_088 | -2.473 | 1 |
| skin | direct_osdr | finetuned_or_direct | 151 | 0 | 1 | 0 | 6.51e-09 | WGAN_FEATURE_030 | -9.095 | 1 |
| kidney | archs4_pretrain_osdr_finetune | finetuned_or_direct | 135 | 1000 | 0 | 0 | 0.1342 | WGAN_FEATURE_040 | -0.423 | 0.917 |
| kidney | direct_osdr | finetuned_or_direct | 135 | 0 | 0 | 0 | 0.2432 | WGAN_FEATURE_055 | -0.8608 | 0.7914 |
| retina | archs4_pretrain_osdr_finetune | finetuned_or_direct | 76 | 1187 | 0 | 0 | 0.5477 | WGAN_FEATURE_130 | 1.276 | 1 |
| retina | direct_osdr | finetuned_or_direct | 76 | 0 | 0 | 0 | 0.4958 | WGAN_FEATURE_028 | -2.151 | 1 |
| skeletal_muscle | archs4_pretrain_osdr_finetune | finetuned_or_direct | 191 | 1412 | 0 | 0 | 0.1614 | WGAN_FEATURE_049 | 1.494 | 0.6045 |

Strict WGAN signals in the main tissues:

- spleen `archs4_pretrain_osdr_finetune` `pre_finetune_projection`: 137 LOO-stable critic features, top WGAN_FEATURE_192.
- thymus `archs4_pretrain_osdr_finetune` `finetuned_or_direct`: 132 LOO-stable critic features, top WGAN_FEATURE_046.
- skeletal_muscle `archs4_pretrain_osdr_finetune` `pre_finetune_projection`: 120 LOO-stable critic features, top WGAN_FEATURE_231.
- thymus `archs4_pretrain_osdr_finetune` `pre_finetune_projection`: 117 LOO-stable critic features, top WGAN_FEATURE_008.
- spleen `archs4_pretrain_osdr_finetune` `finetuned_or_direct`: 107 LOO-stable critic features, top WGAN_FEATURE_160.
- liver `archs4_pretrain_osdr_finetune` `pre_finetune_projection`: 77 LOO-stable critic features, top WGAN_FEATURE_249.
- skeletal_muscle `direct_osdr` `finetuned_or_direct`: 44 LOO-stable critic features, top WGAN_FEATURE_091.
- skin `archs4_pretrain_osdr_finetune` `pre_finetune_projection`: 23 LOO-stable critic features, top WGAN_FEATURE_254.
- liver `archs4_pretrain_osdr_finetune` `finetuned_or_direct`: 18 LOO-stable critic features, top WGAN_FEATURE_111.
- liver `direct_osdr` `finetuned_or_direct`: 4 LOO-stable critic features, top WGAN_FEATURE_248.
- lung `archs4_pretrain_osdr_finetune` `pre_finetune_projection`: 4 LOO-stable critic features, top WGAN_FEATURE_134.
- skin `archs4_pretrain_osdr_finetune` `finetuned_or_direct`: 3 LOO-stable critic features, top WGAN_FEATURE_084.
- lung `archs4_pretrain_osdr_finetune` `finetuned_or_direct`: 3 LOO-stable critic features, top WGAN_FEATURE_042.

No strict WGAN signal was observed for kidney or retina. Kidney and retina had exploratory pretrained-projection random-effects hits, but no LOO-stable features.

## Skeletal Muscle Split Results

| muscle group | mode | score set | query | ARCHS4 ref | RE FDR hits | LOO-stable hits | min RE FDR | top feature | effect | top max LOO FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| edl | archs4_pretrain_osdr_finetune | finetuned_or_direct | 32 | 1412 | 22 | 0 | 5.92e-08 | WGAN_FEATURE_040 | 1.264 | 1 |
| edl | archs4_pretrain_osdr_finetune | pre_finetune_projection | 32 | 1412 | 14 | 0 | 1.95e-10 | WGAN_FEATURE_111 | 1.544 | 1 |
| edl | direct_osdr | finetuned_or_direct | 32 | 0 | 0 | 0 | 1 | WGAN_FEATURE_114 | -7.533 | 1 |
| gastrocnemius | archs4_pretrain_osdr_finetune | finetuned_or_direct | 30 | 1412 | 5 | 0 | 0.001081 | WGAN_FEATURE_194 | 1.721 | 1 |
| gastrocnemius | archs4_pretrain_osdr_finetune | pre_finetune_projection | 30 | 1412 | 1 | 1 | 2.93e-06 | WGAN_FEATURE_063 | -1.299 | 0.004069 |
| gastrocnemius | direct_osdr | finetuned_or_direct | 30 | 0 | 0 | 0 | 1 | WGAN_FEATURE_250 | -0.3855 | 1 |
| quadriceps | archs4_pretrain_osdr_finetune | finetuned_or_direct | 46 | 1412 | 1 | 0 | 0.004613 | WGAN_FEATURE_086 | 0.7265 | 1 |
| quadriceps | archs4_pretrain_osdr_finetune | pre_finetune_projection | 46 | 1412 | 22 | 1 | 6.03e-09 | WGAN_FEATURE_217 | 1.34 | 0.05714 |
| quadriceps | direct_osdr | finetuned_or_direct | 46 | 0 | 3 | 0 | 4.64e-16 | WGAN_FEATURE_086 | 2.497 | 0.7435 |
| soleus | archs4_pretrain_osdr_finetune | finetuned_or_direct | 53 | 1412 | 107 | 12 | 1.39e-07 | WGAN_FEATURE_253 | 4.724 | 7.18e-05 |
| soleus | archs4_pretrain_osdr_finetune | pre_finetune_projection | 53 | 1412 | 11 | 2 | 9.40e-11 | WGAN_FEATURE_231 | -3.404 | 1.17e-05 |
| soleus | direct_osdr | finetuned_or_direct | 53 | 0 | 33 | 0 | 1.22e-09 | WGAN_FEATURE_145 | 5.636 | 1 |
| tibialis_anterior | archs4_pretrain_osdr_finetune | finetuned_or_direct | 30 | 1412 | 0 | 0 | 0.1136 | WGAN_FEATURE_119 | -0.4267 | 1 |
| tibialis_anterior | archs4_pretrain_osdr_finetune | pre_finetune_projection | 30 | 1412 | 0 | 0 | 1 | WGAN_FEATURE_099 | 0.3317 | 1 |
| tibialis_anterior | direct_osdr | finetuned_or_direct | 30 | 0 | 0 | 0 | 0.9387 | WGAN_FEATURE_178 | 9.421 | 1 |

Strict split-muscle WGAN signals:

- gastrocnemius `archs4_pretrain_osdr_finetune` `pre_finetune_projection`: 1 LOO-stable critic features, top WGAN_FEATURE_063.
- quadriceps `archs4_pretrain_osdr_finetune` `pre_finetune_projection`: 1 LOO-stable critic features, top WGAN_FEATURE_217.
- soleus `archs4_pretrain_osdr_finetune` `finetuned_or_direct`: 12 LOO-stable critic features, top WGAN_FEATURE_253.
- soleus `archs4_pretrain_osdr_finetune` `pre_finetune_projection`: 2 LOO-stable critic features, top WGAN_FEATURE_231.

Soleus is the clearest split-muscle WGAN result. Gastrocnemius and quadriceps have only one strict frozen-projection feature each. EDL and tibialis anterior do not pass strict LOO stability by WGAN, although EDL has exploratory random-effects hits.

## Did ARCHS4 Pretraining Help

Yes, but the most useful score set varies by tissue. ARCHS4 pretraining strongly increased strict WGAN signal in spleen, thymus, liver, skin, lung, and the frozen skeletal-muscle projection. It also enabled the soleus split signal. Fine-tuning improved generated-expression fit to OSDR but sometimes reduced condition separation, most notably whole skeletal muscle where the frozen projection had 120 strict features and the fine-tuned score set had zero.

This means both WGAN score sets should be retained: `pre_finetune_projection` asks how OSDR separates in the ARCHS4-trained reference critic, while `finetuned_or_direct` asks how separation remains after adapting the critic to OSDR.

## Comparison To expiMap

Current expiMap multi-tissue FDR summary:

| tissue | direct_osdr_nb_50epoch | archs4_1000ref_nb_ref50_query50 | archs4_allref_nb_ref100_query50 |
| --- | --- | --- | --- |
| lung | 0 | 0 | 0 |
| retina | 0 | 0 | 0 |
| skeletal_muscle | 2 | 39 | 0 |
| skin | 0 | 561 | 208 |
| spleen | 194 | 671 | 104 |
| thymus | 953 | 797 | 725 |

WGAN is broader than expiMap for representation-level signal, especially in skeletal muscle and liver, but expiMap remains more interpretable because its latent dimensions are Reactome or targeted modules. The main disagreement is skeletal muscle: expiMap all-reference Reactome did not produce a final strict aggregate muscle pathway call, while WGAN shows strong whole-muscle direct and frozen-reference representation separation plus a strict soleus split signal. That supports muscle-type-specific follow-up, not a named pathway claim from WGAN alone.

## Comparison To OntoVAE

Top strict OntoVAE rows from the current local summary:

| tissue | group | mode | score set | RE FDR hits | LOO-stable hits | min RE FDR |
| --- | --- | --- | --- | --- | --- | --- |
| thymus |  | archs4_pretrain_osdr_finetune | pre_finetune_projection | 928 | 797 | 7.06e-51 |
| thymus |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 902 | 730 | 2.74e-44 |
| thymus |  | hvg_archs4_pretrain_osdr_finetune | finetuned_or_direct | 260 | 146 | 2.28e-42 |
| thymus |  | hvg_archs4_pretrain_osdr_finetune | pre_finetune_projection | 261 | 139 | 2.77e-41 |
| spleen |  | archs4_pretrain_osdr_finetune | pre_finetune_projection | 339 | 129 | 3.50e-13 |
| spleen |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 350 | 127 | 1.06e-10 |
| thymus |  | direct_osdr | finetuned_or_direct | 731 | 114 | 3.22e-42 |
| liver |  | archs4_pretrain_osdr_finetune | pre_finetune_projection | 120 | 48 | 1.25e-06 |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | finetuned_or_direct | 60 | 27 | 7.41e-17 |
| liver |  | archs4_pretrain_osdr_finetune | finetuned_or_direct | 82 | 13 | 3.38e-04 |
| spleen |  | hvg_archs4_pretrain_osdr_finetune | pre_finetune_projection | 48 | 9 | 1.05e-05 |
| spleen |  | hvg_archs4_pretrain_osdr_finetune | finetuned_or_direct | 54 | 7 | 8.12e-07 |
| skeletal_muscle | soleus | archs4_pretrain_osdr_finetune | pre_finetune_projection | 24 | 7 | 7.30e-06 |
| skeletal_muscle | quadriceps | archs4_pretrain_osdr_finetune | finetuned_or_direct | 63 | 3 | 8.96e-05 |
| liver |  | hvg_archs4_pretrain_osdr_finetune | pre_finetune_projection | 18 | 3 | 4.25e-05 |
| liver |  | hvg_archs4_pretrain_osdr_finetune | finetuned_or_direct | 17 | 3 | 4.90e-05 |
| liver |  | direct_osdr | finetuned_or_direct | 10 | 2 | 1.80e-04 |
| skeletal_muscle |  | archs4_pretrain_osdr_finetune | pre_finetune_projection | 6 | 2 | 8.33e-06 |
| skeletal_muscle |  | hvg_archs4_pretrain_osdr_finetune | finetuned_or_direct | 2 | 2 | 2.60e-06 |
| spleen |  | direct_osdr | finetuned_or_direct | 534 | 1 | 9.06e-08 |

OntoVAE is currently the stronger biology-facing companion to expiMap because it reports named Reactome pathway/program scores and decoder gene associations. WGAN is complementary: it often finds stronger or broader FLT-vs-GC separation, but the features need attribution before they become biological modules.

## Prior-Literature Alignment

The WGAN muscle results align with the prior-work expectation that spaceflight effects are muscle-type-specific and strongest in unloaded/postural muscle. The clearest split result is soleus, matching the prior emphasis on postural muscle, calcium handling, contractile remodeling, mitochondrial/metabolic stress, proteostasis, and ECM remodeling. Because WGAN features are unnamed, this is evidence of separability rather than direct recovery of those pathways.

Thymus and spleen remain broad immune-organ signals across WGAN, expiMap, and OntoVAE. Liver is supported by WGAN and OntoVAE, with expiMap yielding more limited pathway evidence. Lung and skin have smaller WGAN signals; kidney and retina remain weak or non-stable across methods.

## Plot Locations

Representative WGAN visualizations:

- whole skeletal muscle direct: `outputs/wgan_skeletal_muscle/direct_osdr/analysis/top_wgan_feature_shift_heatmap.png`
- whole skeletal muscle frozen reference projection: `outputs/wgan_skeletal_muscle/archs4_pretrain_osdr_finetune/pretrained_query_analysis/top_wgan_feature_shift_heatmap.png`
- soleus reference fine-tuned: `outputs/wgan_skeletal_muscle_splits/soleus/archs4_pretrain_osdr_finetune/analysis/top_wgan_feature_shift_heatmap.png`
- spleen frozen reference projection: `outputs/wgan_spleen/archs4_pretrain_osdr_finetune/pretrained_query_analysis/top_wgan_feature_shift_heatmap.png`
- thymus reference fine-tuned: `outputs/wgan_thymus/archs4_pretrain_osdr_finetune/analysis/top_wgan_feature_shift_heatmap.png`

Each analysis directory also contains `wgan_feature_pca.png`, `wgan_feature_pca_by_accession.png`, `wgan_feature_umap.png`, and `wgan_feature_umap_by_accession.png`.

## Remaining Limitations

- WGAN features are not pathway modules. Add critic-gradient or permutation attribution before calling genes or pathways.
- These are single-seed WGAN runs. Stable claims should be replicated across seeds.
- The frozen projection can show stronger condition separation than the fine-tuned model. That is useful but should be interpreted as reference-space separation, not query-adapted biology.
- Small split groups, especially gastrocnemius, EDL, and tibialis anterior, are accession-limited.
