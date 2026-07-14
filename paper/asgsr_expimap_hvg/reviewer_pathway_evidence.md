# Pathway-level robustness interpretation

The labels below integrate five directional checks: ssGSEA, preranked GSEA, leave-one-project-out prediction, three full reference-query training seeds, and adjustment for atlas-derived broad composition proxies. They are descriptive evidence categories, not hypothesis-test significance levels.

> **Scope update (July 14, 2026):** The sections below document the original thymus, skin, liver, and soleus matrix in Table S24. Corrected spleen and kidney evidence is in Table S27 and the revised 16-pathway paper-facing set is in `table_2_retained_pathway_evidence.tsv`. Spleen is now a main positive tissue, kidney is secondary, and soleus is supplementary sensitivity evidence.

- **Triangulated:** all five checks support the primary direction.
- **Internally robust, incomplete conventional support:** held-out, seed, and composition checks support the direction, but one or both conventional methods do not. These are the clearest expiMap-specific complementary hypotheses.
- **Method-supported, model-sensitive:** conventional methods and held-out projects support the direction, but full-pipeline seed or composition sensitivity does not.
- **Sensitivity-dependent:** the pathway does not meet either reproducibility pattern.

## Thymus

- **triangulated:** DNA repair, RHOA cytoskeletal cycle
- **method-supported, model-sensitive:** Extracellular matrix organization, Innate TLR signaling, Mitotic cell cycle, T-cell receptor signaling
- **internally robust, incomplete conventional support:** Lymphoid-stromal interactions

## Skin

- **triangulated:** Chromatin-modifying enzymes, DNA repair, Hedgehog signaling, Sphingolipid metabolism
- **internally robust, incomplete conventional support:** Cell-cell junction organization
- **method-supported, model-sensitive:** Gap-junction trafficking, Keratinization, Phase II detoxification

## Liver

- **triangulated:** MHC class II antigen presentation, T-cell receptor signaling
- **sensitivity-dependent:** Cytochrome P450, Glutathione conjugation, Regulation of insulin secretion
- **method-supported, model-sensitive:** Extracellular matrix organization, Rho-family GTPase cycle

## Soleus

- **method-supported, model-sensitive:** DNA repair, Glycosaminoglycan metabolism
- **sensitivity-dependent:** Fatty-acid metabolism, Cytokine signaling, Striated-muscle contraction, Extracellular matrix degradation, Immune-system signaling
