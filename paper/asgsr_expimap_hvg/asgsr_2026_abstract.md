# Cross-Mission expiMap Analysis Recovers Established Tissue Responses and Identifies Complementary Pathway Shifts in Mouse Spaceflight Transcriptomes

## Authors and affiliations

**Jason Trinh**<sup>1,2</sup> (presenting author); **James A. Casaletto**<sup>3</sup>; **Walter Alvarado**<sup>4</sup>

<sup>1</sup> Space Life Sciences Training Program, NASA Ames Research Center, Moffett Field, California, USA

<sup>2</sup> University of California, Berkeley, Berkeley, California, USA

<sup>3</sup> Blue Marble Space Institute of Science, NASA Ames Research Center, Mountain View, California, USA

<sup>4</sup> NASA Ames Research Center, Moffett Field, California, USA

## Abstract

Mission differences can obscure recurring organ responses to spaceflight. We asked whether expiMap, a pathway-constrained variational autoencoder, could identify reproducible biological programs in mouse spaceflight transcriptomes. NASA Open Science Data Repository bulk RNA-sequencing samples were mapped into tissue-matched expiMap models trained on non-spaceflight ARCHS4 data, with approximately 2,000 highly variable genes linked to current mouse Reactome pathways. Each accession or mission project received equal weight in the summary of flight-ground pathway shifts, which were evaluated using conventional gene-set enrichment, held-out-project prediction, three complete reference-query trainings, cell-composition sensitivity analyses, and member-gene review. Four tissues produced reproducible patterns. In skin, flight was associated with lower chromatin regulation, DNA repair, Hedgehog signaling, sphingolipid metabolism, and cell-junction programs, consistent with reduced tissue maintenance and barrier coordination. Thymus showed lower DNA-repair and cytoskeletal programs, together with a model-specific lower lymphoid-stromal interaction score. Liver showed lower MHC class II antigen-presentation and T-cell receptor scores despite heterogeneous metabolic responses. Spleen showed coordinated lower T-cell receptor signaling, neutrophil degranulation, and C-type lectin receptor signaling across five unconfounded projects and three trainings; all three programs had preranked GSEA false discovery rates below 0.05. These results recover established tissue responses and add testable hypotheses about skin and thymic tissue maintenance, hepatic immune communication, and splenic innate effector programs. Independent cell-resolved and functional studies are needed to identify which cells drive these pathway shifts and whether they affect biological function.

## Acknowledgment

This work was supported by [funding agency] (grant #[grant number] to [mentor initials]). Funding for the Space Life Sciences Training Program was provided by the NASA Science Mission Directorate's Biological and Physical Sciences Division and NASA Ames Research Center.

## AI usage disclosure

Microsoft Copilot ([model/version]) assisted with code generation, data processing and analysis, interpretation, visualization generation, and text editing. Jason Trinh reviewed and verified all included AI-assisted output.
