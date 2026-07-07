# Multi-Tissue API GLARE Summary

All scopes use NASA OSDR API-derived expression/count inputs. GLARE uses log2(CPM+1) expression aligned to the matching Tabula Muris Senis FACS tissue where available; DESeq2 uses the matched raw-count inputs.

## Run Status

| tissue | FLT/GC profiles | studies | per-study GLARE | DGEA | MOBER | aggregate silhouette FLT/GC | top recurring DGEA/GLARE pathways |
| --- | --- | --- | --- | --- | --- | --- | --- |
| kidney | 67/66 | 6 | 6/6 | yes | yes | 0.0915/0.0514 | Peptide Chain Elongation; 3 Utr Mediated Translational Regulation; Nonsense Mediated Decay Enhanced By The Exon Junctio... |
| liver | 117/112 | 12 | 12/12 | yes | yes | 0.0825/0.142 | Metabolism Of Amino Acids And Derivatives; Immunoregulatory Interactions Between A Lymphoid And A Non Lymphoid Cell; Gp... |
| lung | 39/37 | 3 | 3/3 | yes | yes | 0.238/0.161 | Adaptive Immune System; Cytokine Signaling In Immune System; Interferon Signaling; Innate Immune System; Gpvi Mediated... |
| skeletal muscle | 92/93 | 13 | 13/13 | yes | yes | -0.00706/-0.0991 | Metabolism Of Proteins; Influenza Life Cycle; Developmental Biology; Respiratory Electron Transport; Respiratory Electr... |
| skeletal muscle: edl | 15/15 | 2 | 2/2 | yes | yes | -0.0318/0.0996 | Orc1 Removal From Chromatin; Scfskp2 Mediated Degradation Of P27 P21; Activation Of Nf Kappab In B Cells; Assembly Of T... |
| skeletal muscle: gastrocnemius | 13/17 | 3 | 3/3 | yes | yes | 0.203/0.0998 | Axon Guidance; Developmental Biology; Tca Cycle And Respiratory Electron Transport; Metabolism Of Proteins; Respiratory... |
| skeletal muscle: quadriceps | 22/22 | 4 | 4/4 | yes | yes | 0.19/0.139 | Tca Cycle And Respiratory Electron Transport; Respiratory Electron Transport; Respiratory Electron Transport Atp Synthe... |
| skeletal muscle: soleus | 27/24 | 3 | 3/3 | yes | yes | 0.0385/-0.117 | Mitochondrial Fatty Acid Beta Oxidation; Translation; Metabolism Of Mrna; Srp Dependent Cotranslational Protein Targeti... |
| skeletal muscle: tibialis anterior | 15/15 | 2 | 2/2 | yes | yes | 0.174/0.162 | Destabilization Of Mrna By Auf1 Hnrnp D0; Regulation Of Apoptosis; Downstream Signaling Events Of B Cell Receptor Bcr;... |
| skin | 80/71 | 6 | 6/6 | yes | yes | 0.266/-0.0343 | Translation; Influenza Viral Rna Transcription And Replication; Peptide Chain Elongation; Srp Dependent Cotranslational... |
| spleen | 55/53 | 6 | 6/6 | yes | yes | 0.159/0.1 | Developmental Biology; Extracellular Matrix Organization; Adaptive Immune System; Signaling By Fgfr; Downstream Signali... |
| thymus | 62/53 | 5 | 5/5 | yes | yes | -0.0221/0.0797 | Cell Cycle; Cell Cycle Mitotic; Dna Replication; Mrna Processing; Processing Of Capped Intron Containing Pre Mrna |

## Interpretation Notes

- The PCA/UMAP points are genes, not samples. The plot subtitles list the number of source FLT or GC profiles used to learn that gene representation.
- Cluster descriptions are evidence labels: DGEA-enriched clusters have significant same-study DEG over-representation; Reactome-enriched clusters have significant pathway over-representation; ambiguous clusters lack either label at FDR < 0.05.
- Per-study DESeq2 recurrence is the primary FLT-vs-GC evidence. GLARE is used as module discovery that supports the DGEA signal when the same pathways recur in both layers.
- Retina is skipped unless a matching TMS FACS pretraining tissue becomes available.
