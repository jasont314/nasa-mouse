# SLSTP 2026 mouse spaceflight transcriptomics speaker notes

Target length: 12-15 minutes. Suggested pacing totals 13:35.

## 1. Interpretable and generative models for mouse spaceflight (0:10)

This project uses machine learning to study mouse bulk RNA-seq from NASA spaceflight experiments. One model asks which pathways change. The other asks whether realistic synthetic samples can improve a tissue-specific comparison of flight and ground control.

## 2. Learn how spaceflight changes living systems (0:30)

The project has two connected goals. First, identify which genes, pathways, and biological systems differ between flight and ground-control mice. Second, learn meaningful expression patterns and generate realistic profiles that can support that comparison.

## 3. Autoencoders compress thousands of genes into a few features (0:35)

An autoencoder compresses thousands of gene measurements into a smaller set of features and then reconstructs the original profile. The output bars are close to, but not exactly the same as, the input bars. Samples with similar compressed profiles sit near one another. MOBER uses this compressed space to reduce study identity, while expiMap connects its features to known Reactome gene programs.

## 4. Study identity dominates the expression structure (0:20)

This is EDL skeletal muscle from two OSDR studies. Color marks the study and shape marks flight or ground control. The two colors form separate UMAP clusters, while triangles and circles overlap within each cluster. Study identity is the strongest visible structure.

## 5. MOBER tries to remove study identity (0:30)

Study differences can be larger than the biological effect. MOBER combines an autoencoder with a source discriminator that tries to identify the originating study. The encoder makes study labels harder to predict while still reconstructing the expression profiles. The animation starts with separated studies and moves them into a shared distribution.

## 6. MOBER reduces study separation, but FLT and GC still overlap (0:25)

This UMAP compares the same two muscle studies before and after MOBER. Before correction, study identity dominates the layout. After MOBER, samples from the studies mix more closely, but flight and ground control still overlap. Removing visible study structure does not automatically reveal a strong spaceflight axis.

## 7. expiMap assigns each latent feature to a known pathway (0:30)

GLARE and MOBER learn their latent structure from expression. expiMap instead connects each latent node to a known mouse Reactome pathway through a masked decoder. That lets us compare flight and ground control one named biological program at a time. The next slide shows how those program scores are summarized.

## 8. Program scores summarize pathway changes within one tissue (0:35)

The table is an illustrative example for one tissue, not a project result. Each sample receives one score per Reactome program. We compare the average flight and ground-control scores and subtract GC from FLT. The rows are ordered from higher red shifts to lower blue shifts. The heatmap at right shows the complete program-by-sample output that these comparisons summarize.

## 9. Each program was compared with prior literature (0:25)

The top heatmap shows the pathway scores alone. The lower heatmap shows the same scores after each pathway name was reviewed against spaceflight literature. Green agrees with prior work, blue adds a related or complementary interpretation, orange is uncertain, red conflicts, and gray has little effect. The score values do not change; the colors describe the evidence attached to each pathway.

## 10. Five tissues showed recurring pathway patterns (0:45)

Thymus showed lower repair, cytoskeletal, and stromal-interaction programs. Skin showed lower chromatin regulation, repair, Hedgehog, sphingolipid, and cell-junction programs. Liver showed lower MHC class II and T-cell receptor scores. Spleen combined lower T-cell receptor, neutrophil-degranulation, and C-type lectin programs. Kidney showed higher ECM proteoglycan, WNT, and IGF-transport programs, suggesting structural and growth-factor remodeling.

## 11. What is synthetic transcriptomics? (0:20)

A generator learns patterns from measured RNA-seq and creates new numeric expression profiles for a chosen tissue and FLT or GC condition. These profiles are useful model outputs, but they are not new animals or independent biological measurements.

## 12. Diffusion: denoise from noise (0:25)

The left animation shows diffusion for images: start with random noise and remove it step by step until an image appears. The right side applies the same idea to gene expression. The phrase data manifold on the slide means the region occupied by real expression profiles.

## 13. Conditional WGAN-GP: generator versus critic (0:25)

A GAN trains two networks against each other. The generator makes expression profiles. The critic, which is similar to a discriminator, scores how much each profile resembles measured data. Tissue, FLT or GC, study, and sample material tell the generator what type of profile to make.

## 14. Building the RNA-seq generator (0:35)

We compared different ways to prepare the data, handle study differences, train the generator, and specify what it should make. Feature-space options included all shared genes, highly variable genes, Reactome genes, and a mouse mapping of the human L1000 landmark panel. The selected diffusion model uses the 974-gene L1000 mapping with TPM and ARCHS4-training MaxAbs scaling. It first learns broad mouse tissue patterns from ARCHS4, then adapts to OSDR.

## 15. Diffusion best reproduced the measured expression distribution (0:35)

We compared WGAN-GP with diffusion. Correlation measures gene-expression agreement. Coverage F1 asks whether generated samples cover the same regions as measured samples. Real-versus-synthetic accuracy asks whether a prediction model can tell the two apart, where 0.5 is ideal. Distribution distance measures overall separation, so lower is better. Diffusion had better coverage, near-chance discrimination, and the smaller distance, so we used it for the biological analysis.

## 16. Diffusion learns tissue structure from noise (0:20)

The same 1,024 generated samples begin as noise, develop tissue structure by step 200, and reach their final expression profiles at step zero. PCA compresses the 974 gene values to two axes so the movement can be plotted.

## 17. FLT and GC overlap in the global PCA view (0:20)

We first ask whether flight and ground-control profiles separate when all tissues and studies are viewed together. They overlap across the first two principal components, so there is no clear global condition axis. This motivates a closer look at the tissue and study structure.

## 18. Tissue and study structure dominate the PCA space (0:25)

Recoloring the same profiles by tissue and study reveals much stronger structure. The left panel shows tissue and the right shows study. Circles are real samples and crosses are generated samples. We therefore compare flight and ground control within each tissue and account for study in the downstream analysis.

## 19. Does adding synthetic data improve FLT vs GC prediction? (0:25)

We train the same tissue-specific classifier twice: once with real profiles and once with real plus synthetic profiles. Both versions are tested on the same real samples. This directly asks whether adding synthetic profiles improves flight-versus-ground-control prediction.

## 20. Real + synthetic vs real-only balanced accuracy across 27 tissues (0:35)

Each line compares real-only training with real-plus-synthetic training for flight-versus-ground-control classification in one tissue analysis. All scores come from held-out real OSDR samples. Balanced accuracy averages flight sensitivity and ground-control specificity. Teal points identify the 18 analyses where balanced accuracy, AUROC, and average precision all held or improved both overall and after giving each study equal weight. Coral marks the nine mixed results where at least one measure declined.

## 21. Compare what each classifier finds important (0:25)

Both classifiers see the same 974 genes and are evaluated on the same real samples. Features important in both classifiers are reinforced. Features that become important only after adding synthetic samples are promoted.

## 22. We compared feature importance at three levels (0:30)

Each analysis asks which features support flight-versus-ground-control prediction. Individual-gene permutation and SHAP score one gene at a time. Grouped permutation and grouped SHAP score a Reactome pathway together. Consensus ranking compares real and generated gene rankings and tests compact panels.

## 23. The three analyses retained different types of results (0:30)

The individual-gene analysis retained 21 associations. The grouped analysis retained 10 Reactome pathways, and consensus ranking retained 49 gene associations in compact panels. Promoted means the feature became important after synthetic samples were added. Reinforced means it was important with and without synthetic samples.

## 24. All three analyses found candidates in thymus, skin and spleen (0:35)

Thymus, skin, and spleen have retained candidates in every column. This agreement is why the biological interpretation focuses on these three tissues. The full analysis retains the single-method candidates for follow-up.

## 25. Feature selection and literature review answer separate questions (0:30)

The model comparison labels a candidate as promoted or reinforced. The literature review then fixes the candidate, tissue, and flight direction and asks whether prior work agrees, supports a related mechanism, is ambiguous, or offers no close match.

## 26. Thymus: all three analyses point to lower cell-cycle activity in flight (0:40)

Thymus has the clearest agreement. Fifteen individual genes, seven Reactome groups, and a 16-gene compact panel point toward lower mitotic activity in flight. This could reflect reduced proliferative renewal or a smaller fraction of cycling thymocytes. The expression data cannot distinguish those explanations.

## 27. Skin and spleen show different flight-associated patterns (0:40)

In skin, flight-higher PLSCR1 and two flight-higher necroptosis pathways suggest an interferon-linked regulated cell-death response. This is a complementary hypothesis because direct spaceflight evidence for skin necroptosis is limited. In spleen, flight-higher LOXL1, RAI14, PTPRK, and MYL9 suggest extracellular-matrix, adhesion, and cytoskeletal remodeling. This is also complementary because prior work supports the component mechanisms but not the same spleen-flight gene directions. The four genes did not form a significant shared Reactome pathway.

## 28. What the models can and cannot tell us (0:25)

The models summarize patterns in the RNA-seq data and narrow the list of tissues, pathways, and genes to examine. They do not establish mechanism. Wet-lab experiments must confirm the biological changes and distinguish gene regulation from shifts in cell composition.

## 29. Thank you (0:10)

Acknowledge James Casaletto, SLSTP, NASA OSDR, ARCHS4, Reactome, ChatGPT, and Claude, then invite questions.
