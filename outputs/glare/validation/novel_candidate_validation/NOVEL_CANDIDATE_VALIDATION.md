# Novel/Complementary GLARE Candidate Validation

## Bottom Line

Validated or near-validated candidates:

- **Aggregate skeletal muscle regulatory modules** are the best non-thymus lead: circadian/BMAL1 and Cyclin E/G1-S validate across 13 studies. L1CAM is a weaker follow-up candidate; TLR3 and CREB/RAS are weaker despite nominal FDR support because empirical/random-set and direction consistency are not as strong.
- **Soleus hidden modules** are mixed: NGF signaling validates, circadian is moderate, Cyclin E/G1-S is statistically significant but has almost zero mean effect, and DARPP-32 is weak. OXPHOS/fatty-acid oxidation remains the validation anchor, not the novel result.
- **Kidney hidden modules** are plausible but weaker than muscle: membrane trafficking and insulin receptor signaling pass the current criteria; fatty-acid/TAG/ketone is directionally consistent but empirical support is borderline; organic ion transport and mitochondrial beta oxidation do not validate yet.
- **Muscle subtype comparison** supports subtype specificity: soleus has the clearest OXPHOS/fatty-acid suppression and NGF-like signal; aggregate muscle has the clearest circadian and L1CAM signals; TA has an underpowered Cyclin E/G1-S signal; EDL is mostly not validated.

## Candidate Summary

| candidate_group | tissue | module_class | clean_term | theme | studies_tested | mean_flight_minus_ground | direction_consistency | median_empirical_abs_p | combined_welch_fdr_bh | strict_supported_studies | direction_call | validation_call |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kidney_hidden_modules | kidney | glare_only | Fatty Acid Triacylglycerol And Ketone Body Metabolism | lipid/TAG/ketone metabolism | 6 | 0.0557 | 0.833 | 0.065 | 0.0338 | 0 | mixed/reorganized DEG | moderate validated candidate |
| kidney_hidden_modules | kidney | glare_only | Organic Cation Anion Zwitterion Transport | organic ion transport | 6 | 0.00711 | 0.5 | 0.025 | 0.276 | 0 | mixed/reorganized DEG | not validated by current stack |
| kidney_hidden_modules | kidney | glare_only | Mitochondrial Fatty Acid Beta Oxidation | mitochondrial fatty-acid beta oxidation | 6 | 0.0244 | 0.667 | 0.2 | 0.284 | 0 | ambiguous/no validated direction | not validated by current stack |
| kidney_hidden_modules | kidney | glare_only | Membrane Trafficking | membrane trafficking/tubular handling | 6 | 0.0338 | 0.667 | 0.025 | 0.018 | 1 | mixed/reorganized DEG | strong validated candidate |
| kidney_hidden_modules | kidney | glare_only | Signaling By Insulin Receptor | insulin receptor/metabolic signaling | 6 | 0.0282 | 0.667 | 0.02 | 0.0204 | 0 | mixed/reorganized DEG | strong validated candidate |
| kidney_hidden_modules | kidney | glare_only | Trans Golgi Network Vesicle Budding | vesicle budding/trafficking | 6 | 0.0301 | 0.667 | 0.045 | 0.0727 | 0 | mixed/reorganized DEG | weak/needs follow-up |
| kidney_validation_anchor | kidney | intersection | Peptide Chain Elongation | known translation validation anchor | 6 | -0.0255 | 0.333 | 0.285 | 0.114 | 1 | GC-up DEG-supported | not validated by current stack |
| kidney_validation_anchor | kidney | intersection | Tca Cycle And Respiratory Electron Transport | known mitochondrial validation anchor | 6 | 0.0206 | 0.5 | 0.125 | 0.0511 | 1 | mixed/reorganized DEG | weak/needs follow-up |
| skeletal_muscle_regulatory | skeletal_muscle | glare_only | Creb Phosphorylation Through The Activation Of Ras | CREB/RAS signaling | 13 | -0.0113 | 0.538 | 0.38 | 0.00136 | 0 | mixed/reorganized DEG | not validated by current stack |
| skeletal_muscle_regulatory | skeletal_muscle | glare_only | Bmal1 Clock Npas2 Activates Circadian Expression | circadian clock/regulatory | 13 | -0.079 | 0.846 | 0.01 | 2.87e-06 | 1 | mixed/reorganized DEG | strong validated candidate |
| skeletal_muscle_regulatory | skeletal_muscle | glare_only | Cyclin E Associated Events During G1 S Transition | cell-cycle/regeneration | 13 | 0.0385 | 0.692 | 0 | 2.77e-05 | 2 | mixed/reorganized DEG | strong validated candidate |
| skeletal_muscle_regulatory | skeletal_muscle | glare_only | L1Cam Interactions | adhesion/NMJ-like signaling | 13 | 0.00823 | 0.538 | 0.02 | 0.00104 | 1 | mixed/reorganized DEG | weak/needs follow-up |
| skeletal_muscle_regulatory | skeletal_muscle | glare_only | Trif Mediated Tlr3 Signaling | innate antiviral/stress signaling | 13 | 0.0161 | 0.538 | 0.12 | 0.00855 | 0 | mixed/reorganized DEG | weak/needs follow-up |
| soleus_hidden_modules | skeletal_muscle_soleus | glare_only | Circadian Clock | circadian clock | 3 | -0.0882 | 0.667 | 0.07 | 2.33e-04 | 1 | mixed/reorganized DEG | moderate validated candidate |
| soleus_hidden_modules | skeletal_muscle_soleus | glare_only | Darpp 32 Events | DARPP-32 signaling | 3 | 0.0215 | 0.667 | 0.39 | 0.0472 | 0 | FLT-up DEG-supported | not validated by current stack |
| soleus_hidden_modules | skeletal_muscle_soleus | glare_only | Cyclin E Associated Events During G1 S Transition | cell-cycle/regeneration | 3 | -0.00375 | 0.667 | 0 | 7.16e-04 | 1 | mixed/reorganized DEG | strong validated candidate |
| soleus_hidden_modules | skeletal_muscle_soleus | glare_only | Signalling By Ngf | NGF signaling | 3 | -0.0778 | 0.667 | 0.04 | 0.0204 | 0 | mixed/reorganized DEG | strong validated candidate |
| soleus_hidden_modules | skeletal_muscle_soleus | glare_only | Ngf Signalling Via Trka From The Plasma Membrane | NGF/TRKA signaling | 3 | -0.071 | 0.667 | 0.16 | 0.0572 | 0 | mixed/reorganized DEG | weak/needs follow-up |
| soleus_validation_anchor | skeletal_muscle_soleus | intersection | Mitochondrial Fatty Acid Beta Oxidation | known OXPHOS/fatty-acid validation anchor | 3 | -0.581 | 1 | 0 | 4.94e-08 | 2 | GC-up DEG-supported | strong validated candidate |
| soleus_validation_anchor | skeletal_muscle_soleus | glare_only | Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins | known OXPHOS validation anchor | 3 | -0.33 | 1 | 0 | 1.06e-05 | 1 | GC-up DEG-supported | strong validated candidate |

## Muscle Subtype Matrix

| clean_term | tissue | module_class | studies_tested | mean_flight_minus_ground | direction_consistency | median_empirical_abs_p | combined_welch_fdr_bh |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Bmal1 Clock Npas2 Activates Circadian Expression | skeletal_muscle | glare_only | 13 | -0.079 | 0.846 | 0.01 | 2.87e-06 |
| Circadian Clock | skeletal_muscle_soleus | glare_only | 3 | -0.0882 | 0.667 | 0.07 | 2.33e-04 |
| Circadian Clock | skeletal_muscle_tibialis_anterior | glare_only | 2 | -0.0452 | 1 | 0.275 | 0.372 |
| Cyclin E Associated Events During G1 S Transition | skeletal_muscle | glare_only | 13 | 0.0385 | 0.692 | 0 | 2.77e-05 |
| Cyclin E Associated Events During G1 S Transition | skeletal_muscle_soleus | glare_only | 3 | -0.00375 | 0.667 | 0 | 7.16e-04 |
| Cyclin E Associated Events During G1 S Transition | skeletal_muscle_tibialis_anterior | glare_only | 2 | 0.113 | 1 | 0 | 0.0193 |
| Cyclin E Associated Events During G1 S Transition | skeletal_muscle_quadriceps | glare_only | 4 | 0.0191 | 0.5 | 0.205 | 0.0359 |
| Cyclin E Associated Events During G1 S Transition | skeletal_muscle_gastrocnemius | glare_only | 3 | -0.00217 | 0.333 | 0.49 | 0.274 |
| Mitochondrial Fatty Acid Beta Oxidation | skeletal_muscle_soleus | intersection | 3 | -0.581 | 1 | 0 | 4.94e-08 |
| Ngf Signalling Via Trka From The Plasma Membrane | skeletal_muscle_soleus | glare_only | 3 | -0.071 | 0.667 | 0.16 | 0.0572 |
| Ngf Signalling Via Trka From The Plasma Membrane | skeletal_muscle_quadriceps | glare_only | 4 | 0.027 | 0.75 | 0.565 | 0.274 |
| Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins | skeletal_muscle | glare_only | 13 | -0.101 | 0.615 | 0.01 | 2.12e-06 |
| Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins | skeletal_muscle_soleus | glare_only | 3 | -0.33 | 1 | 0 | 1.06e-05 |
| Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins | skeletal_muscle_gastrocnemius | glare_only | 3 | -0.053 | 0.667 | 0 | 0.0988 |
| Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins | skeletal_muscle_quadriceps | glare_only | 4 | -0.0829 | 0.5 | 0.155 | 0.0993 |
| Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins | skeletal_muscle_tibialis_anterior | glare_only | 2 | 0.00281 | 0.5 | 0.03 | 0.11 |
| Respiratory Electron Transport Atp Synthesis By Chemiosmotic Coupling And Heat Production By Uncoupling Proteins | skeletal_muscle_edl | glare_only | 2 | -0.0299 | 0.5 | 0.1 | 0.372 |
| Signalling By Ngf | skeletal_muscle_soleus | glare_only | 3 | -0.0778 | 0.667 | 0.04 | 0.0204 |
| Signalling By Ngf | skeletal_muscle_quadriceps | glare_only | 4 | 0.0326 | 0.75 | 0.43 | 0.176 |

## Per-Candidate Interpretation

### 1. Aggregate Skeletal Muscle Regulatory Modules

**Validated:** BMAL1/CLOCK/NPAS2 circadian and Cyclin E/G1-S. These are the best complementary muscle candidates because they recur across 13 studies and pass random-set checks.
**Weaker:** L1CAM, TLR3, and CREB/RAS have nominal combined FDR but weaker direction consistency or empirical/random-set support, so treat them as follow-up rather than claims.

### 2. Soleus Hidden Modules

**Validation anchor:** mitochondrial fatty-acid beta oxidation and OXPHOS are strongly GC-up / FLT-suppressed.
**Novel candidate:** NGF signaling is the most defensible soleus hidden module. Circadian is moderate; Cyclin E/G1-S has small mean effect despite statistical support; DARPP-32 is currently weak.

### 3. Kidney GLARE-Only Modules

**Validated candidates:** membrane trafficking and insulin receptor signaling.
**Borderline:** fatty-acid/TAG/ketone metabolism and trans-Golgi vesicle budding.
**Not validated yet:** organic cation/anion/zwitterion transport and mitochondrial fatty-acid beta oxidation. These may still be biologically plausible but do not pass the current stack.

### 4. Muscle Subtype Comparison

Subtype comparison is useful as an organizing analysis, not just another candidate module. The current validation suggests:
- Soleus: strongest OXPHOS/fatty-acid suppression and NGF-like candidate.
- Aggregate muscle: strongest circadian/L1CAM/Cyclin E evidence.
- TA: Cyclin E/G1-S appears but only two studies, so underpowered.
- Quadriceps/gastrocnemius/EDL: useful context but weaker candidate validation.

## Retained Files

- `candidate_validation_summary.tsv`
- `candidate_per_study_validation.tsv`
- `muscle_subtype_candidate_matrix.tsv`
