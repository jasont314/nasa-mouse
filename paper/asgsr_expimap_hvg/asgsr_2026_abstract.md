# ASGSR 2026 Abstract

## Title

Reference-guided expiMap prioritizes reproducible tissue-specific pathway shifts in mouse spaceflight transcriptomes

## Authors and affiliations

**Jason Trinh** (presenting author)<sup>1</sup>

<sup>1</sup> NASA Space Life Sciences Training Program (NASA/SLSTP), NASA Ames Research Center, Moffett Field, California, USA

## Abstract

Spaceflight affects multiple organs, but differences among missions can obscure biological responses that recur across studies. We asked whether expiMap, a pathway-constrained variational autoencoder, could identify reproducible biological programs in mouse spaceflight transcriptomes. NASA Open Science Data Repository bulk RNA-sequencing samples were mapped into tissue-matched expiMap models trained on non-spaceflight ARCHS4 data, with approximately 2,000 highly variable genes linked to current mouse Reactome pathways. Each accession or mission project received equal weight in the summary of flight-ground pathway shifts, which were evaluated using conventional gene-set enrichment, held-out-project prediction, three complete reference-query trainings, cell-composition sensitivity analyses, and member-gene review. Four tissues produced reproducible patterns. In skin, flight was associated with lower chromatin regulation, DNA repair, Hedgehog signaling, sphingolipid metabolism, and cell-junction programs, consistent with reduced tissue maintenance and barrier coordination. Thymus showed lower DNA-repair and cytoskeletal programs, together with a model-specific lower lymphoid-stromal interaction score. Liver showed lower MHC class II antigen-presentation and T-cell receptor scores despite heterogeneous metabolic responses. Spleen showed coordinated lower T-cell receptor signaling, neutrophil degranulation, and C-type lectin receptor signaling across five unconfounded projects and three trainings; all three programs had preranked GSEA false discovery rates below 0.05. These results recover established tissue responses and add testable hypotheses about skin and thymic tissue maintenance, hepatic immune communication, and splenic innate effector programs. Independent cell-resolved and functional studies are needed to identify which cells drive these pathway shifts and whether they affect biological function. Acknowledgment: NASA/SLSTP.
