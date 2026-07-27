# Expanded pathway-family review

> **Scope update (July 14, 2026):** This family review describes the original thymus, skin, liver, and soleus seed-2020 screen and is retained as a transparent review queue. Corrected spleen and kidney pathways were reviewed separately in Tables S27-S29. The revised paper-facing evidence set is `source_data/table_2_retained_pathway_evidence.tsv`.

This audit starts from every active program, selects the top within-tissue decile by absolute study-balanced effect, adds directionally stable programs through rank 40, and retains all pathways prespecified for the main figures. It then consolidates overlapping Reactome terms into manually reviewed process families.

The expanded set contains 153 pathway records and 37 nonredundant tissue-family records. It is a review queue and evidence synthesis, not a new significance threshold.

## Main conclusions

- Thymus: the expanded primary-model terms suggested one lower proliferation and genome-maintenance axis plus lower niche and higher innate-stromal responses. Later seed and composition checks retain only a narrower lower DNA-repair and RHOA cytoskeletal core, with lower lymphoid-stromal interaction as an internally robust but conventionally incomplete hypothesis.
- Skin: higher cutaneous-muscle and innate-phagocytic programs recover known full-thickness-skin compartments. Lower communication and epidermal programs remain coherent. Opposite nested glutathione and phase-II scores mean that detoxification cannot be assigned one uniform direction.
- Liver: additional MHC-I, complement, neutrophil, and interleukin terms provide context for a lower immune-effector family. The integrated analysis retains lower MHC class II antigen presentation and T-cell receptor signaling; matrix and Rho-family terms remain follow-up hypotheses. Opposite nested TLR4 and broad-TLR scores are not interpreted as selective TLR4 activation.
- Soleus: lower PI3K-AKT, WNT, and VEGF scores add literature-aligned trophic context to the primary matrix-disassembly hypothesis, but no reviewed soleus pathway passes all five later robustness checks. Cornified-envelope and keratin terms are excluded as tissue-incongruent.

No newly reviewed complementary family is sufficiently distinct and mission-stable to replace the manuscript's central tissue narratives. The review adds aligned context, strengthens existing families, and narrows claims where nested pathways conflict.

## Family decisions

| Tissue | Process family | Role | Decision | Representative | Paths | FLT-GC range |
| --- | --- | --- | --- | --- | ---: | ---: |
| Liver | Xenobiotic and redox metabolism | context sensitive | existing context | Cytochrome P450 | 4 | -0.096 to +0.239 |
| Liver | Toll-like receptor branches | context sensitive | supplement only | Toll Like Receptor 4 TLR4 Cascade | 2 | -0.091 to +0.195 |
| Liver | Organelle trafficking and proteostasis | context sensitive | supplement only | Copi Dependent Golgi To ER Retrograde Traffic | 3 | -0.166 to +0.100 |
| Liver | Adaptive and innate effector programs | complementary | existing core axis | MHC class II antigen presentation | 7 | -0.201 to -0.097 |
| Liver | Matrix, mechanical, and vascular regulation | complementary | existing support axis | Rho-family GTPase cycle | 7 | -0.148 to -0.096 |
| Liver | Proliferation and homology-directed repair | context sensitive | supplement only | Homology Directed Repair | 3 | -0.128 to -0.125 |
| Liver | Tissue-incongruent labels | not interpretable | exclude | Cardiac Conduction | 2 | -0.114 to -0.085 |
| Liver | Lipid and endocrine regulation | context sensitive | existing context | Regulation of insulin secretion | 7 | -0.101 to +0.151 |
| Liver | Broad receptor and stress signaling | not interpretable | exclude | Signaling By GPCR | 4 | -0.096 to +0.090 |
| Skin | Epidermal differentiation and regeneration | aligned | existing core axis | Keratinization | 5 | -1.221 to -0.194 |
| Skin | Chromatin, proliferation, and genome maintenance | context sensitive | existing support axis | Chromatin-modifying enzymes | 8 | -0.663 to -0.223 |
| Skin | Cell communication, adhesion, and cytoskeleton | complementary | existing core axis | Gap-junction trafficking | 7 | -0.408 to -0.197 |
| Skin | Cutaneous striated muscle | aligned | add aligned context | Muscle Contraction | 2 | +0.297 to +0.339 |
| Skin | Sensory and broad signaling labels | not interpretable | exclude | Visual Phototransduction | 5 | -0.325 to +0.215 |
| Skin | Innate immune and phagocytic response | aligned | add aligned context | Fcgamma Receptor Fcgr Dependent Phagocytosis | 2 | +0.243 to +0.284 |
| Skin | Secretory trafficking and protein processing | context sensitive | supplement only | Copi Dependent Golgi To ER Retrograde Traffic | 3 | -0.273 to +0.239 |
| Skin | Detoxification and endocrine metabolism | context sensitive | revise existing claim | Phase II detoxification | 4 | -0.255 to +0.306 |
| Soleus | Immune and cytokine programs | context sensitive | existing context | Immune-system signaling | 4 | -0.117 to +0.435 |
| Soleus | Cornified-envelope and keratin labels | not interpretable | exclude | Formation Of The Cornified Envelope | 2 | +0.260 to +0.403 |
| Soleus | Hemostatic and platelet programs | context sensitive | supplement only | Hemostasis | 2 | +0.182 to +0.230 |
| Soleus | Metabolism and transmembrane transport | context sensitive | existing context | Fatty-acid metabolism | 7 | -0.198 to +0.181 |
| Soleus | Broad GPCR and Rho-family signaling | context sensitive | supplement only | G Alpha I Signalling Events | 4 | -0.174 to +0.347 |
| Soleus | Matrix disassembly and support | complementary | existing core axis | Extracellular matrix degradation | 4 | -0.122 to +0.207 |
| Soleus | Growth, apoptosis, and DNA-damage response | context sensitive | supplement only | DNA repair | 4 | +0.140 to +0.338 |
| Soleus | Trophic, adhesion, and vascular signaling | aligned | add aligned context | Pip3 Activates Akt Signaling | 3 | -0.122 to -0.100 |
| Soleus | Membrane trafficking and intracellular transport | context sensitive | supplement only | Vesicle Mediated Transport | 4 | +0.098 to +0.163 |
| Soleus | Contractile and sensory programs | context sensitive | existing context | Striated-muscle contraction | 2 | +0.090 to +0.162 |
| Thymus | Proliferation and genome maintenance | aligned | existing core axis | Mitotic cell cycle | 10 | -1.438 to -0.437 |
| Thymus | Secretory trafficking and proteostasis | context sensitive | supplement only | Copi Dependent Golgi To ER Retrograde Traffic | 6 | -0.929 to +0.439 |
| Thymus | Thymocyte-niche interaction and cytoskeleton | complementary | existing core axis | RHOA cytoskeletal cycle | 3 | -0.812 to -0.501 |
| Thymus | Respiratory, metabolic, and endocrine state | context sensitive | supplement only | Aerobic Respiration And Respiratory Electron Transport | 5 | -0.580 to +0.716 |
| Thymus | Broad receptor signaling | not interpretable | exclude | Signal Transduction | 3 | +0.331 to +0.452 |
| Thymus | Stromal matrix and TGF-beta response | complementary | existing core axis | Extracellular matrix organization | 2 | +0.353 to +0.432 |
| Thymus | Cell death and stress response | context sensitive | add interpretive context | Apoptosis | 3 | -0.448 to -0.381 |
| Thymus | Tissue-incongruent labels | not interpretable | exclude | Neurotransmitter Receptors And Postsynaptic Signal Transmission | 3 | -0.358 to +0.385 |
| Thymus | Innate, hemostatic, and vascular-associated response | complementary | existing support axis | Innate TLR signaling | 6 | +0.312 to +0.506 |
| Thymus | Adaptive T-cell signaling | aligned | existing core axis | T-cell receptor signaling | 1 | -0.255 to -0.255 |
