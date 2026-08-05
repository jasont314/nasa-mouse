# SLSTP 2026 mouse spaceflight transcriptomics speaker notes

Target length: 12-15 minutes. Planned speaking time: about 13 minutes.

## 1. Interpretable and generative models for mouse spaceflight (0:15)

This project uses machine learning to study mouse bulk RNA-seq from NASA spaceflight experiments. The first model asks which pathways change. The second asks whether realistic synthetic profiles can improve a tissue-specific FLT versus GC analysis.

## 2. One dataset, two machine-learning questions (0:35)

Both parts begin with the same OSDR flight and ground-control data. expiMap organizes expression into interpretable Reactome programs. The generative pipeline models the expression distribution and then tests whether generated profiles improve classification and gene ranking.

## 3. Autoencoders compress expression into a latent space (0:40)

An autoencoder takes a high-dimensional gene-expression profile, compresses it into a small set of latent variables, and reconstructs the original profile. Nearby points in latent space have similar expression. In a standard autoencoder those axes may have no clear biological meaning. expiMap constrains them with known gene programs.

## 4. Program scores become pathway-level changes (0:45)

These are observed OSDR examples from the expiMap analysis. Each sample receives a latent score for each Reactome program. The table shows project-centered mean scores for flight and ground control, followed by flight minus ground. The formal analysis estimates the change inside each project and then combines project effects. All four examples are lower in flight; a lower latent score does not by itself prove biochemical inhibition.

## 5. Four tissues produced the clearest pathway patterns (0:50)

Thymus showed lower repair, cytoskeletal, and stromal-interaction programs. Skin showed lower chromatin regulation, repair, Hedgehog, sphingolipid, and cell-junction programs. Liver showed lower MHC class II and T-cell receptor scores. Spleen combined lower T-cell receptor, neutrophil-degranulation, and C-type lectin programs.

## 6. What is synthetic transcriptomics? (0:25)

A generator learns the distribution of measured expression and samples new numeric profiles for a chosen tissue and FLT or GC context. We use those profiles to test classifiers and rank genes. They are model output, not new animals or independent biological measurements.

## 7. Small studies and study effects complicate tissue comparisons (0:30)

OSDR covers many tissues, but its 1,610 profiles are spread across 75 accessions. ARCHS4 supplies a much larger mouse reference. The model must preserve tissue and condition structure without simply learning study identity.

## 8. A generator learns by competing with a critic (0:30)

A WGAN-GP alternates between two updates. The generator turns noise and biological conditions into an expression profile. The critic compares measured and generated profiles. Their competition teaches the generator to match the observed distribution while retaining the requested tissue, flight or ground-control, accession, and material context.

## 9. We built a configurable bulk RNA-seq generation pipeline (0:45)

The pipeline can change data scope, transformation, harmonization, model, training source, and conditioning. The branch used here applies TPM, MaxAbs scaling, 974 landmarks, ARCHS4 pretraining, OSDR adaptation, and conditioning on tissue, FLT or GC, accession, and material type.

## 10. Diffusion turns noise into a conditioned sample (0:25)

This teaching animation starts from random points. During reverse diffusion, the model repeatedly predicts and removes noise while receiving the requested condition. In the RNA-seq model, the output is a vector of 974 gene values conditioned on tissue, flight or ground control, accession, and material type.

## 11. DDIM matched expression and reduced separability (0:40)

We compared WGAN-GP with the conditional diffusion model. Both matched expression well. DDIM had higher F1, adversarial accuracy close to chance, and lower distributional distance, so the remaining analysis uses DDIM.

## 12. Diffusion learns tissue structure from noise (0:25)

The same generated profiles begin as noise, develop structure by timestep 200, and reach their tissue-conditioned regions at timestep zero. All three panels share PCA axes.

## 13. Tissue and study structure dominate the PCA space (0:25)

These panels use the same locked coordinates. The left panel colors each profile by tissue, while the right colors it by OSDR accession. Circles are observed profiles and crosses are matched DDIM profiles. Both tissue and study structure are reproduced in the generated data.

## 14. Flight condition is subtler than study structure (0:25)

The left panel now colors the same profiles by flight or ground-control condition, while the right repeats the accession view. FLT and GC overlap much more than the study clusters. This is why the downstream analysis estimates FLT-GC effects within accession rather than treating the pooled separation as biology.

## 15. The primary comparison changes only the training source (0:40)

Every matched classifier uses all 974 genes, the same real-fitted scaler, the same outer split, and one regularization value selected from real training data. Real-only, synthetic-only, and real-plus-synthetic models are evaluated on the same held-out real profiles. This isolates training source within the classifier analysis.

## 16. Consensus ranking is a secondary panel analysis (0:25)

Real and generated profiles rank the same 974 genes. Combining those rankings can move a gene into or out of a compact top-k panel. This is useful for pathway interpretation, but it does not isolate training source as directly as the matched all-gene comparison.

## 17. Correlated genes dilute marginal importance (0:35)

Suppose several genes carry the same pathway signal. Ridge can divide weight among them. If I shuffle Gene A, Genes B and C remain, so held-out performance changes little and Gene A receives low permutation importance. Consensus ranking can still retain the group. Low individual importance means replaceable in this classifier, not biologically irrelevant.

## 18. Matched augmentation helped many tissues, but not all (0:35)

The pooled multi-tissue classifier declined with augmentation. In separate tissue models, real plus synthetic passed all pooled and accession-macro balanced-accuracy, AUROC, and average-precision checks in 18 of 27 analyses. The coral examples show why a balanced-accuracy gain alone is not enough when another metric declines.

## 19. Matched and consensus results overlap only partly (0:30)

The matched analysis retained 21 BH-FDR tissue-gene associations, and consensus ranking retained 49. Eleven appear in both. Matched results are primary evidence that synthetic training changes classifier behavior. Consensus-only results are secondary panel evidence.

## 20. Consensus selection and literature are separate dimensions (0:30)

For the 49 consensus associations, promoted or reinforced describes feature selection. Aligning, complementary, ambiguous, or unmatched describes the literature review. These labels answer different questions and can occur in any combination.

## 21. The matched screen covered all 27 tissue analyses (0:25)

Four tissues had both matched utility and retained BH-FDR genes. Fourteen more passed the utility gate without a retained gene. Nine failed at least one mean metric. Predictive improvement does not automatically produce a biological candidate.

## 22. The secondary consensus inventory spans ten tissue analyses (0:20)

This slide lists all 49 consensus associations. Rows separate FLT direction and selection status. Gene color gives the independent literature classification.

## 23. Thymus is strongest across both analyses (0:45)

The matched analysis retained 15 thymus genes, seven promoted after augmentation, and 26 significant Reactome terms led by mitotic cell cycle. Nine genes overlap the 16-gene consensus panel, which adds correlated cell-cycle members. Together they support lower proliferative renewal or fewer cycling thymocytes in flight.

## 24. Soleus remains a secondary consensus result (0:35)

The consensus analysis reinforces lower Bdh1, Ech1, Bnip3, and Decr1 with higher Tpm1, a coherent mitochondrial and lipid-metabolism panel. In the fixed all-gene comparison, balanced accuracy rose slightly but AUROC and average precision fell, so the matched gate did not pass.

## 25. Three additional tissues have matched genes (0:25)

Liver contributes four flight-lower shared-importance genes without pathway enrichment. Skin Plscr1 and spleen Loxl1 are supported by both analyses. Pooled muscle improves prediction in the matched classifier, but its individual gene interpretation remains consensus-level.

## 26. Additional panels remain consensus-only (0:20)

Kidney, adrenal gland, gastrocnemius, and tibialis anterior add smaller tissue-specific consensus candidates. They did not pass the matched gene gate, so they remain exploratory panel-level hypotheses.

## 27. Use matched tests for contribution and consensus for programs (0:30)

Conditional DDIM produced realistic profiles. Matched classifiers show where synthetic training changes held-out-real prediction and gene importance. Consensus ranking organizes correlated biological panels. Thymus is strongest across both; all association statistics still come from observed OSDR samples.

## 28. Thank you (0:10)

Acknowledge James Casaletto, SLSTP, NASA OSDR, ARCHS4, Reactome, and NASA Ames compute, then invite questions.
