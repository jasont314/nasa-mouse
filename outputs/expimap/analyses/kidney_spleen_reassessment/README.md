# Kidney and spleen expiMap reassessment

This reassessment corrects the two main comparability issues in the historical screen: kidney now uses all eligible ARCHS4 reference samples, and spleen uses batch-aware HVG selection after excluding singleton ARCHS4 series from HVG ranking only. The primary spleen contrast also excludes OSD-288 because recorded strain is disjoint between flight and ground groups; the full query remains available as a sensitivity. It applies three full training seeds, ssGSEA, preranked GSEA, held-out-project direction checks, and broad atlas-derived composition-proxy adjustment.

The pathway ranking below uses relative effect magnitude rather than a hard FDR gate. A top-decile label means the pathway is among the largest absolute seed-2020 expiMap shifts within that tissue.

## Kidney

The corrected model uses 2,464 ARCHS4 samples, 336 Reactome programs, and 135 mapped OSDR samples. The primary effect summary uses 6 unconfounded projects. Across the primary top-decile pathways, 83% retained one direction across all three complete training runs.
Compared with the historical HVG model, the corrected model has Spearman r=0.49, 64% directional agreement, and top-decile Jaccard overlap 0.30.

### Largest pathways after triangulation

- **Platelet Degranulation:** median shift +0.143; 5/5 directional checks; triangulated; held-out projects 83% concordant.
- **Metabolism:** median shift +0.128; 5/5 directional checks; triangulated; held-out projects 67% concordant.
- **Ecm Proteoglycans:** median shift +0.096; 5/5 directional checks; triangulated; held-out projects 83% concordant.
- **Regulation Of Insulin Like Growth Factor Igf Transport And Uptake By Insulin Like Growth Factor Binding Proteins Igfbps:** median shift +0.077; 5/5 directional checks; triangulated; held-out projects 83% concordant.
- **Signaling By Wnt:** median shift +0.067; 5/5 directional checks; triangulated; held-out projects 100% concordant.
- **Metabolism Of Proteins:** median shift +0.136; 4/5 directional checks; method-supported, model-sensitive; held-out projects 67% concordant.
- **Response To Elevated Platelet Cytosolic Ca2:** median shift +0.135; 4/5 directional checks; method-supported, model-sensitive; held-out projects 83% concordant.
- **Degradation Of The Extracellular Matrix:** median shift +0.088; 4/5 directional checks; method-supported, model-sensitive; held-out projects 83% concordant.
- **Smooth Muscle Contraction:** median shift +0.084; 4/5 directional checks; method-supported, model-sensitive; held-out projects 67% concordant.
- **Developmental Biology:** median shift +0.048; 4/5 directional checks; method-supported, model-sensitive; held-out projects 100% concordant.
- **G Alpha I Signalling Events:** median shift -0.093; 4/5 directional checks; internally robust, incomplete conventional support; held-out projects 83% concordant.
- **Fatty Acid Metabolism:** median shift -0.111; 4/5 directional checks; method-supported, model-sensitive; held-out projects 67% concordant.
- **Regulation Of Expression And Function Of Type I Classical Cadherins:** median shift +0.104; 3/5 directional checks; sensitivity-dependent; held-out projects 67% concordant.
- **Plasma Lipoprotein Assembly Remodeling And Clearance:** median shift +0.099; 3/5 directional checks; sensitivity-dependent; held-out projects 83% concordant.
- **Metabolism Of Steroid Hormones:** median shift +0.081; 3/5 directional checks; sensitivity-dependent; held-out projects 83% concordant.
- **Keratinization:** median shift +0.066; 3/5 directional checks; sensitivity-dependent; held-out projects 100% concordant.
- **Cell Junction Organization:** median shift +0.057; 3/5 directional checks; sensitivity-dependent; held-out projects 83% concordant.
- **Prednisone Adme:** median shift +0.016; 3/5 directional checks; method-supported, model-sensitive; held-out projects 100% concordant.
- **Golgi To Er Retrograde Transport:** median shift -0.037; 3/5 directional checks; sensitivity-dependent; held-out projects 83% concordant.
- **Cellular Response To Chemical Stress:** median shift -0.068; 3/5 directional checks; sensitivity-dependent; held-out projects 83% concordant.

## Spleen

The corrected model uses 6,289 ARCHS4 samples, 360 Reactome programs, and 109 mapped OSDR samples. The primary effect summary uses 5 unconfounded projects. Across the primary top-decile pathways, 68% retained one direction across all three complete training runs.
Compared with the historical HVG model, the corrected model has Spearman r=0.20, 56% directional agreement, and top-decile Jaccard overlap 0.16.

### Largest pathways after triangulation

- **Axon Guidance:** median shift +0.118; 5/5 directional checks; triangulated; held-out projects 80% concordant.
- **Slc Mediated Transport Of Amino Acids:** median shift +0.071; 5/5 directional checks; triangulated; held-out projects 100% concordant.
- **Collagen Formation:** median shift +0.036; 5/5 directional checks; triangulated; held-out projects 80% concordant.
- **Dap12 Interactions:** median shift -0.003; 5/5 directional checks; triangulated; held-out projects 100% concordant.
- **Post Translational Modification Synthesis Of Gpi Anchored Proteins:** median shift -0.072; 5/5 directional checks; triangulated; held-out projects 80% concordant.
- **Class A 1 Rhodopsin Like Receptors:** median shift -0.080; 5/5 directional checks; triangulated; held-out projects 100% concordant.
- **C Type Lectin Receptors Clrs:** median shift -0.105; 5/5 directional checks; triangulated; held-out projects 100% concordant.
- **Neutrophil Degranulation:** median shift -0.182; 5/5 directional checks; triangulated; held-out projects 100% concordant.
- **Tcr Signaling:** median shift -0.216; 5/5 directional checks; triangulated; held-out projects 100% concordant.
- **Signaling By Rho Gtpases Miro Gtpases And Rhobtb3:** median shift +0.135; 4/5 directional checks; sensitivity-dependent; held-out projects 40% concordant.
- **Intra Golgi And Retrograde Golgi To Er Traffic:** median shift +0.090; 4/5 directional checks; sensitivity-dependent; held-out projects 60% concordant.
- **Fc Epsilon Receptor Fceri Signaling:** median shift +0.035; 4/5 directional checks; method-supported, model-sensitive; held-out projects 80% concordant.
- **Signaling By The B Cell Receptor Bcr:** median shift +0.004; 4/5 directional checks; method-supported, model-sensitive; held-out projects 100% concordant.
- **Platelet Activation Signaling And Aggregation:** median shift -0.038; 4/5 directional checks; internally robust, incomplete conventional support; held-out projects 100% concordant.
- **Metabolism Of Lipids:** median shift -0.061; 4/5 directional checks; internally robust, incomplete conventional support; held-out projects 100% concordant.
- **Interleukin 1 Signaling:** median shift -0.113; 4/5 directional checks; sensitivity-dependent; held-out projects 60% concordant.
- **Tnfs Bind Their Physiological Receptors:** median shift -0.127; 4/5 directional checks; method-supported, model-sensitive; held-out projects 100% concordant.
- **Peptide Hormone Metabolism:** median shift -0.132; 4/5 directional checks; sensitivity-dependent; held-out projects 60% concordant.
- **Regulation Of T Cell Activation By Cd28 Family:** median shift -0.132; 4/5 directional checks; internally robust, incomplete conventional support; held-out projects 100% concordant.
- **Mhc Class Ii Antigen Presentation:** median shift -0.147; 4/5 directional checks; sensitivity-dependent; held-out projects 60% concordant.

## Interpretation rule

A pathway is suitable for biological follow-up when it combines a large relative effect with repeatable direction across seeds and projects and support from at least one conventional enrichment method. Off-tissue labels, tiny effects, and pathways that reverse across complete training runs are not promoted even if an isolated statistic is favorable.
