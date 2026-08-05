# SLSTP 2026 mouse spaceflight transcriptomics speaker notes

Target length: 12-15 minutes. Planned speaking time: about 13 minutes.

## 1. Interpretable and generative models for mouse spaceflight (0:15)

This project uses machine learning to study mouse bulk RNA-seq from NASA spaceflight experiments. One model asks which pathways change. The other asks whether realistic synthetic samples can improve a tissue-specific comparison of flight and ground control.

## 2. One dataset, two machine-learning questions (0:35)

Both parts start with the same OSDR flight and ground-control data. expiMap summarizes expression as Reactome pathways. The generative pipeline learns the broader expression patterns and tests whether generated samples improve FLT-GC prediction and gene ranking.

## 3. Autoencoders compress thousands of genes into a few features (0:40)

An autoencoder compresses thousands of gene measurements into a smaller set of features and then reconstructs the original profile. This compressed representation is often called latent space. Samples with similar expression sit near one another. expiMap gives each compressed feature a biological meaning by connecting it to a known Reactome gene program.

## 4. Program scores become pathway-level changes (0:45)

Each OSDR sample receives one score per Reactome program. The table shows average flight and ground-control scores, followed by flight minus ground control. We first compare FLT and GC within each study, where the samples share the same experimental setting, and then combine those study-level changes. A lower program score in flight does not by itself prove that the pathway is biochemically inhibited.

## 5. Four tissues produced the clearest pathway patterns (0:50)

Thymus showed lower repair, cytoskeletal, and stromal-interaction programs. Skin showed lower chromatin regulation, repair, Hedgehog, sphingolipid, and cell-junction programs. Liver showed lower MHC class II and T-cell receptor scores. Spleen combined lower T-cell receptor, neutrophil-degranulation, and C-type lectin programs.

## 6. What is synthetic transcriptomics? (0:25)

A generator learns patterns from measured RNA-seq and creates new numeric expression profiles for a chosen tissue and FLT or GC condition. These profiles are useful model outputs, but they are not new animals or independent biological measurements.

## 7. Small studies and study effects complicate tissue comparisons (0:30)

OSDR contains 1,610 profiles from 75 studies. NASA identifiers for those studies are sometimes called accessions. ARCHS4 supplies a much larger mouse reference. The generator must reproduce tissue and condition patterns without mistaking differences between studies for spaceflight biology.

## 8. Diffusion: denoise from noise (0:30)

The left animation shows diffusion for images: start with random noise and remove it step by step until an image appears. The right side applies the same idea to gene expression. The phrase data manifold on the slide means the region occupied by real expression profiles.

## 9. Conditional WGAN-GP: generator versus critic (0:30)

A GAN trains two networks against each other. The generator makes expression profiles. The critic, which is similar to a discriminator, scores how much each profile resembles measured data. Tissue, FLT or GC, study, and sample material tell the generator what type of profile to make.

## 10. We built a configurable bulk RNA-seq generation pipeline (0:45)

We tested different data sources, processing methods, study corrections, models, and conditioning variables. The selected branch converts expression to TPM, rescales each gene to a common range, and uses 974 informative genes. It first learns general mouse tissue patterns from ARCHS4, then adapts to OSDR and receives tissue, FLT or GC, study, and sample-material labels.

## 11. Diffusion best reproduced the measured expression distribution (0:45)

We compared WGAN-GP with diffusion. Correlation measures gene-expression agreement. Coverage F1 asks whether generated samples cover the same regions as measured samples. Real-versus-synthetic accuracy asks whether a prediction model can tell the two apart, where 0.5 is ideal. Distribution distance measures overall separation, so lower is better. Diffusion had better coverage, near-chance discrimination, and the smaller distance, so we used it for the biological analysis.

## 12. Diffusion learns tissue structure from noise (0:25)

The same 1,024 generated samples begin as noise, develop tissue structure by step 200, and reach their final expression profiles at step zero. PCA compresses the 974 gene values to two axes so the movement can be plotted.

## 13. Tissue and study structure dominate the PCA space (0:25)

The left panel colors each sample by tissue and the right colors it by study. Circles are measured samples and crosses are generated samples. Their overlap shows that the generator reproduced both tissue differences and study-related variation.

## 14. Flight condition is subtler than study structure (0:25)

The left panel colors the same samples by FLT or GC. The right repeats the study view. FLT and GC overlap much more than the study clusters, so we estimate the FLT-GC difference separately within each study before combining the evidence.

## 15. A fair test changes only the training data (0:45)

This is a controlled comparison, not a claim that individual animals were matched. Every prediction model uses the same 974 genes, processing, train-test splits, and settings. We train it on real samples, synthetic samples, or both, then test all three versions on the same unseen real samples. Balanced accuracy gives equal weight to FLT and GC accuracy. AUROC measures ranking, and average precision focuses on correct FLT predictions.

## 16. Combine real and synthetic rankings to choose a small gene panel (0:30)

We also rank genes separately in the real and synthetic data. A gene that ranks well in both moves upward in the combined ranking. We then test small panels containing the top 10, 25, 50, or 100 genes. This analysis is useful for finding related gene groups, but the controlled all-gene test is the clearer test of whether synthetic training helps.

## 17. Correlated genes can hide one another in one-at-a-time tests (0:35)

Several genes in one pathway can carry nearly the same information. If we shuffle Gene A, Genes B and C may still support the prediction, so performance changes little. The one-gene score can therefore be small even when the group matters. The combined ranking helps keep these correlated genes together as a panel.

## 18. Adding synthetic samples helped many tissues, but not all (0:40)

A single prediction model across all tissues became slightly worse after synthetic samples were added. Separate tissue models worked better: 18 of 27 analyses improved or held steady on all six checks. The checks cover overall prediction and the average result across studies. Teal marks tissues that passed every check; coral marks mixed results.

## 19. All-gene tests and compact panels answer different questions (0:35)

The controlled all-gene test retained 21 tissue-gene associations. The combined ranking retained 49, and 11 appeared in both. Every association had FDR below 0.05 in measured OSDR data. FDR is the expected fraction of false positives among results called significant after testing many genes.

## 20. Selection and literature evidence are separate (0:35)

Promoted means a gene entered the selected panel only after the synthetic ranking was added. Reinforced means the real-only and combined rankings both selected it. The literature label asks a different question: does prior work directly agree, support a related mechanism, give mixed evidence, or provide no close match?

## 21. The controlled test covered all 27 tissue analyses (0:30)

Four tissue analyses improved prediction and retained at least one significant gene that helped classification. Fourteen more improved prediction but did not retain a gene under both criteria. Nine had a mixed prediction result. Better classification does not automatically produce a biological finding.

## 22. The combined ranking spans ten tissue analyses (0:25)

This slide lists all 49 combined-ranking associations. Rows separate genes that were higher or lower in flight and whether they were promoted or reinforced. Gene color shows the separate literature interpretation.

## 23. Thymus is strongest across both analyses (0:45)

The controlled test retained 15 thymus genes, including seven added after synthetic training. Reactome analysis found 26 significant terms led by mitotic cell cycle. Nine of those genes also appear in the 16-gene combined panel, which adds related cell-cycle genes. Together they suggest less proliferative renewal or fewer cycling thymocytes in flight.

## 24. Soleus has a coherent secondary gene panel (0:35)

The combined ranking found lower Bdh1, Ech1, Bnip3, and Decr1 with higher Tpm1. That pattern suggests changes in mitochondrial turnover and fatty-acid metabolism. However, only one of the three prediction measures improved when synthetic samples were added, so soleus remains a biological hypothesis rather than confirmed synthetic benefit.

## 25. Liver, skin, spleen, and pooled muscle add narrower results (0:30)

Liver contributes four flight-lower genes without pathway enrichment. Skin Plscr1 and spleen Loxl1 appear in both analyses. Pooled muscle prediction improved after adding synthetic samples, but the individual genes come from the secondary combined ranking.

## 26. Additional panels remain exploratory (0:25)

Kidney, adrenal gland, gastrocnemius, and tibialis anterior produced smaller tissue-specific panels. These panels did not meet the controlled all-gene criteria, so we treat them as hypotheses for follow-up rather than firm evidence of synthetic benefit.

## 27. Takeaways (0:30)

Diffusion generated realistic tissue-specific expression profiles. Adding those profiles improved all six prediction checks in 18 of 27 tissue analyses. The controlled and combined-ranking analyses both point most strongly to a lower cell-cycle signal in thymus. All statistical claims about FLT-GC biology still come from measured OSDR samples.

## 28. Thank you (0:10)

Acknowledge James Casaletto, SLSTP, NASA OSDR, ARCHS4, Reactome, and NASA Ames compute, then invite questions.
