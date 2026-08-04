# SLSTP 2026 mouse spaceflight transcriptomics speaker notes

Target length: 12-15 minutes. Planned speaking time: about 13 minutes.

## 1. Interpretable and generative models for mouse spaceflight (0:15)

This project uses machine learning to study mouse bulk RNA-seq from NASA spaceflight experiments. The first model asks which pathways change. The second asks whether realistic synthetic profiles can improve a tissue-specific FLT versus GC analysis.

## 2. One dataset, two machine-learning questions (0:35)

Both parts begin with the same OSDR flight and ground-control data. expiMap organizes expression into interpretable Reactome programs. The generative pipeline models the expression distribution and then tests whether generated profiles improve classification and gene ranking.

## 3. Autoencoders compress expression into a latent space (0:40)

An autoencoder takes a high-dimensional gene-expression profile, compresses it into a small set of latent variables, and reconstructs the original profile. Nearby points in latent space have similar expression. In a standard autoencoder those axes may have no clear biological meaning. expiMap constrains them with known gene programs.

## 4. Reactome pathways make the latent space interpretable (0:45)

For each tissue, ARCHS4 supplies a non-spaceflight reference and Reactome supplies the gene-program mask. expiMap learns the reference, then maps OSDR flight and ground samples as an accession-conditioned query. The analysis used about two thousand highly variable genes and several hundred Reactome programs per tissue. The implementation is under src/expiMap_scarches/nasa_mouse_expimap.

## 5. Four tissues produced the clearest pathway patterns (0:50)

Thymus showed lower repair, cytoskeletal, and stromal-interaction programs. Skin showed lower chromatin regulation, repair, Hedgehog, sphingolipid, and cell-junction programs. Liver showed lower MHC class II and T-cell receptor scores. Spleen combined lower T-cell receptor, neutrophil-degranulation, and C-type lectin programs.

## 6. What is synthetic transcriptomics? (0:25)

A generator learns the distribution of measured expression and samples new numeric profiles for a chosen tissue and FLT or GC context. We use those profiles to test classifiers and rank genes. They are model output, not new animals or independent biological measurements.

## 7. Small studies and study effects complicate tissue comparisons (0:30)

OSDR covers many tissues, but its 1,610 profiles are spread across 75 accessions. ARCHS4 supplies a much larger mouse reference. The model must preserve tissue and condition structure without simply learning study identity.

## 8. We built a configurable bulk RNA-seq generation pipeline (0:45)

The pipeline can change data scope, transformation, harmonization, model, training source, and conditioning. The branch used here applies TPM, MaxAbs scaling, 974 landmarks, ARCHS4 pretraining, OSDR adaptation, and conditioning on tissue, FLT or GC, accession, and material type.

## 9. DDIM matched expression and reduced separability (0:40)

We compared WGAN-GP with the conditional diffusion model. Both matched expression well. DDIM had higher F1, adversarial accuracy close to chance, and lower distributional distance, so the remaining analysis uses DDIM.

## 10. Diffusion learns tissue structure from noise (0:25)

The same generated profiles begin as noise, develop structure by timestep 200, and reach their tissue-conditioned regions at timestep zero. All three panels share PCA axes.

## 11. Tissue and study structure dominate the PCA space (0:25)

These panels use the same locked coordinates. The left panel colors each profile by tissue, while the right colors it by OSDR accession. Circles are observed profiles and crosses are matched DDIM profiles. Both tissue and study structure are reproduced in the generated data.

## 12. Flight condition is subtler than study structure (0:25)

The left panel now colors the same profiles by flight or ground-control condition, while the right repeats the accession view. FLT and GC overlap much more than the study clusters. This is why the downstream analysis estimates FLT-GC effects within accession rather than treating the pooled separation as biology.

## 13. Five arms separate gene ranking from classifier fitting (0:40)

The five arms vary which profiles rank genes and which profiles fit the classifier. The guided arms use real and generated rankings together. One then fits on real profiles only; the other gives generated profiles five percent of the training weight. Association tests still use observed OSDR data.

## 14. Consensus ranking chooses the classifier input genes (0:25)

Real and generated profiles rank the same 974 genes. Combining those rankings can move a gene into or out of the selected top set. The classifier is then trained using only the selected expression columns.

## 15. Synthetic guidance can shift the FLT/GC boundary (0:15)

Opaque points are training profiles and transparent points are held-out real profiles. The panels use the same held-out samples. A useful synthetic-guided feature set changes the fitted boundary so that more held-out labels fall on the correct side.

## 16. Synthetic use helped some tissues and hurt others (0:35)

The pooled classifier declined with generated profiles. Within tissues, the best synthetic-informed candidate improved balanced accuracy in examples such as spleen, thymus, and skin, but declined in cecum, colon, and slightly in liver. The selection rule retained a synthetic arm only when balanced accuracy, AUROC, and average precision were all non-worse than real-only training.

## 17. Synthetic guidance changed ranking, not statistical evidence (0:30)

Twenty-three associations were selected with and without guidance, so we call them reinforced. Twenty-six crossed the repeated selection threshold only with guidance, so we call them promoted. All 49 passed BH FDR in observed OSDR profiles.

## 18. Selection and literature are separate dimensions (0:30)

Promoted or reinforced describes feature selection. Aligning, complementary, ambiguous, or unmatched describes the literature review. These labels answer different questions and can occur in any combination.

## 19. The screen covered all 27 completed tissue analyses (0:20)

The full screen includes 22 canonical tissues and five anatomical muscle groups. Ten analyses had a synthetic-informed BH-FDR association, five had real-data BH-FDR genes without synthetic support, and 12 had no BH-FDR gene in the landmark panel.

## 20. Ten tissue analyses contained synthetic-informed genes (0:20)

This slide lists all 49 synthetic-informed associations. Rows separate FLT direction and selection status. Gene color gives the independent literature classification.

## 21. Thymus points to lower proliferative renewal (0:40)

Thymus produced the clearest promoted panel. Lower mitotic and DNA-replication genes fit prior reports of thymic involution. Hsd17b11 and Etv1 add lipid-handling and T-cell-state hypotheses. Bulk RNA-seq cannot separate transcriptional regulation from cell-composition change.

## 22. Soleus reinforces a mitochondrial and lipid program (0:35)

Soleus improved with real plus generated training. Lower Bdh1, Ech1, Bnip3, and Decr1 with higher Tpm1 support altered oxidative metabolism and contractile remodeling. The real association remains, but its synthetic reinforcement was sensitive to material conditioning.

## 23. Additional tissues produced distinct hypotheses (0:20)

Pooled muscle, kidney, spleen, and skin each produced a separate result. Promoted and reinforced genes are kept on separate rows so the selection claim stays clear.

## 24. Eye, adrenal and muscle-group results remain tissue-specific (0:20)

Eye, adrenal gland, gastrocnemius, and tibialis anterior add smaller tissue-specific candidates. These are follow-up hypotheses rather than one shared systemic signature.

## 25. Synthetic data worked best as a tissue-specific prior (0:25)

Conditional DDIM produced realistic profiles. Tissue-specific ranking and light synthetic training improved held-out prediction in selected tissues. The final biological associations and FDR still come from observed OSDR samples.

## 26. Thank you (0:10)

Acknowledge James Casaletto, SLSTP, NASA OSDR, ARCHS4, Reactome, and NASA Ames compute, then invite questions.
