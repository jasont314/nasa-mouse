# ASGSR 2026 Abstract

## Title

Reference-guided expiMap identifies tissue-specific pathway responses to spaceflight in mouse transcriptomes

## Authors and affiliations

**Jason Trinh** (presenting author)<sup>1</sup>

<sup>1</sup> NASA Space Life Sciences Training Program (NASA/SLSTP), NASA Ames Research Center, Moffett Field, California, USA

## Abstract

Spaceflight alters many organs, but study-specific designs can obscure biological processes shared across missions. We adapted expiMap, a pathway-constrained variational autoencoder, to NASA Open Science Data Repository mouse bulk RNA sequencing. Six tissue-matched non-spaceflight ARCHS4 references were trained on approximately 2,000 highly variable genes connected to current mouse Reactome programs. Thymus, skin, liver, and spleen formed the positive-result set; kidney was secondary and soleus was a supplementary sensitivity analysis. Decoder-oriented flight-minus-ground scores were balanced across accessions or mission projects. We tested directions using single-sample and preranked gene-set enrichment, held-out-project prediction, three complete reference-query trainings, broad cell-composition proxies, and member-gene review. Flight skin had lower chromatin-regulatory, DNA-repair, Hedgehog, sphingolipid, and cell-junction scores, supporting a coordinated lower tissue-maintenance state. Flight thymus had lower DNA-repair and RHOA cytoskeletal scores; lower lymphoid-stromal interaction remained an internally reproducible expiMap-specific hypothesis. Flight liver had lower MHC class II antigen-presentation and T-cell receptor scores beside heterogeneous metabolic responses. Flight spleen had lower T-cell receptor, neutrophil degranulation, and C-type lectin receptor programs across five unconfounded projects and three trainings; all three had preranked-GSEA FDR below 0.05. Kidney showed a reproducible higher ECM-proteoglycan, WNT, and IGF-transport axis, but these effects were composition-attenuated and lacked conventional FDR support. No soleus pathway passed all checks. Reference-guided pathway modeling therefore recovers known tissue biology and adds regulatory, structural, and innate-immune hypotheses while exposing model- and mission-sensitive findings. Independent cell-resolved and functional validation is required. Acknowledgment: NASA/SLSTP.
