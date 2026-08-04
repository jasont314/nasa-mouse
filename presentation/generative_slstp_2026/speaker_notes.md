# SLSTP 2026 generative transcriptomics speaker notes

Target length: 12-15 minutes. Planned speaking time: about 15 minutes.

## 1. Synthetic transcriptomics for mouse spaceflight (0:20)

Introduce the central question: can generated expression help us find tissue-specific FLT versus GC biology? Synthetic profiles support the analysis, but they do not count as additional animals.

## 2. What is synthetic transcriptomics? (0:40)

Synthetic transcriptomes are numeric gene-expression vectors sampled from a model trained on measured RNA-seq data. Conditioning lets us request a tissue and FLT or GC context. These profiles can stress-test a classifier and guide gene ranking. They do not add biological replicates or independent evidence.

## 3. Small studies and study effects complicate tissue comparisons (0:45)

OSDR gives broad tissue coverage, but the profiles are spread across 75 accessions with different mission and assay contexts. ARCHS4 supplies a much larger mouse reference. The challenge is to use that reference without confusing study structure with spaceflight biology.

## 4. Match tissue distributions, then test FLT versus GC biology (0:45)

The generator has two jobs. First, synthetic bulk RNA-seq should reproduce the tissue-defined distributions seen in real data. Second, it should preserve the smaller FLT versus GC signal within each tissue and improve prediction on held-out real samples. Gene effects and BH FDR are always calculated from observed OSDR profiles.

## 5. We built a configurable bulk RNA-seq generation pipeline (1:00)

Each column shows the alternatives available at one pipeline stage; outlines identify the downstream branch. We used ARCHS4 and NASA OSDR across multiple studies and all tissues. The selected path used TPM, training-fitted MaxAbs scaling, 974 mouse landmarks, no global correction, ARCHS4 pretraining plus OSDR adaptation, and a DDIM conditioned on tissue, FLT/GC, accession and material. WGAN-GP and the other preprocessing and harmonization choices remained benchmark alternatives.

## 6. DDIM matched expression and reduced separability (0:50)

The left image is Figure 1C from Lacan and colleagues. Their residual dense denoiser predicts the noise added to a sample using diffusion timestep and tissue. Our implementation adds FLT/GC, accession and material context during OSDR LoRA adaptation. Both generators had high correlation and F1. DDIM had adversarial accuracy near 0.5 and a lower Frechet-distance ratio, so it was harder to separate from real profiles and closer in distribution.

## 7. Diffusion learns tissue structure from noise (0:40)

Read the panels from left to right. The same generated profiles begin as noise, develop structure by timestep 200 and approach tissue-conditioned regions at timestep zero. The axes are shared, so the visual change does not come from rescaling each panel.

## 8. Generated profiles track the real OSDR PCA manifold (0:40)

Circles are locked real OSDR profiles and crosses are matched DDIM profiles in the same PCA space. Generated samples follow the tissue-defined branches. FLT and GC overlap more because condition effects are smaller than tissue effects. The numerical validation on the previous slide tests fidelity directly.

## 9. Five arms separate gene ranking from classifier fitting (0:50)

Each arm makes two decisions: which profiles rank the genes and which profiles fit the classifier. In both guided arms, real and synthetic evidence jointly rank genes. Guided real fit then trains only on observed profiles. Guided 5% also uses condition-recentered synthetic profiles, but they contribute only 5% of total classifier weight. Held-out real profiles determine eligibility, and FLT/GC effects and BH FDR come from observed OSDR profiles only.

## 10. Consensus ranking chooses the classifier input genes (0:35)

Each ranking orders the same 974 genes. Real-only, generated-only and consensus ranking can therefore produce different candidate gene sets. For each ranking, we test the top 10, 25, 50 and 100 genes, and held-out validation chooses the feature count and regularization. Logistic regression is then fitted using only those selected gene-expression columns. Ranking chooses which genes are available to the classifier; classifier training separately learns their coefficients and decision boundary.

## 11. Held-out real samples decide whether guidance helps (0:25)

Opaque profiles belong to the training subset. Transparent profiles are held-out real OSDR samples that never enter gene ranking, classifier fitting or top-k selection. The same held-out profiles are scored by the real-only and synthetic-guided candidates. This schematic shows the intended outcome: a guided boundary that predicts more held-out labels correctly. Balanced accuracy, AUROC and average precision determine whether a synthetic arm is eligible for that tissue; the observed tissue-specific results follow on the next slide.

## 12. Pooling tissues hid useful signal (0:50)

The pooled augmentation test was negative: balanced accuracy fell from 0.754 to 0.737 with real plus synthetic training. Tissue-specific analysis changed the result. Different tissues benefited from different synthetic uses, which argues against one global augmentation policy.

## 13. Synthetic guidance changed ranking, not statistical evidence (0:45)

The blue set contains genes selected stably by real-only ranking, and the teal set contains genes selected stably by the eligible synthetic-guided arm. Thirty-four were real-only, 23 were selected by both arms and classified as reinforced, and 26 were selected only with guidance and classified as promoted. Promoted does not mean biologically novel. All 49 synthetic-informed tissue-gene associations passed BH FDR in observed OSDR profiles.

## 14. Selection and literature are separate dimensions (0:50)

Every association has two labels. Promoted or reinforced describes repeated feature selection. Aligning, complementary, ambiguous or unmatched describes prior literature. Across all 49 associations, 22 aligned, 19 were complementary, four were ambiguous and four were unmatched. Table S16 records the gene-level rationale and source IDs; Table S17 records the citations and evidence relationship.

## 15. The screen covered all 27 completed tissue analyses (0:35)

This is the full analysis coverage: 22 canonical tissues and five anatomical muscle groups. Ten had a synthetic-informed BH-FDR association, five had real BH-FDR genes without synthetic-informed selection, and 12 had no BH-FDR gene in the landmark panel. Every completed tissue result remains visible here.

## 16. Ten tissue analyses contained synthetic-informed genes (0:40)

This is the complete 49-association inventory. Separate rows show FLT-higher or FLT-lower direction and promoted or reinforced selection status. Gene color independently shows aligning, complementary, ambiguous or unmatched literature. FLT directions come from real-data meta-analysis.

## 17. Thymus points to lower proliferative renewal (1:00)

Thymus produced the clearest promoted panel. The lower mitotic and DNA-replication genes agree with prior reports of thymic involution and altered cell-cycle expression after flight. Higher Hsd17b11 and Etv1 add lipid-handling and T-cell-state hypotheses. A matched sensitivity model without material-type conditioning preserved the cell-cycle interpretation. Because this is bulk RNA-seq, the pattern may reflect transcription, cell composition or both.

## 18. Soleus reinforces a mitochondrial and lipid program (0:55)

Soleus improved with real plus generated training. The selected genes were already stable in real-only analysis, so synthetic data reinforced rather than introduced the panel. Lower Bdh1, Ech1, Bnip3 and Decr1, with higher Tpm1, support altered oxidative metabolism and contractile remodeling. In the no-material sensitivity model, the synthetic attribution disappeared. The real OSDR association remains, but the generated-profile contribution is conditioning-sensitive.

## 19. Additional tissues produced distinct hypotheses (0:35)

Promoted and reinforced genes are shown on separate subrows for each tissue. Pooled muscle, kidney, spleen and skin each produced a distinct synthetic-informed result. The rows share a slide for presentation space; each remains a separate hypothesis.

## 20. Eye, adrenal and muscle-group results remain tissue-specific (0:35)

Promoted and reinforced genes remain separated here as well. Eye reinforces lower cytokinesis, adrenal contributes two unmatched candidates, gastrocnemius combines an NF-kappa-B stress signal with an autophagy or myogenesis candidate, and tibialis anterior spans stress, cell-cycle, ganglioside and mitophagy hypotheses.

## 21. Synthetic data worked best as a tissue-specific prior (0:35)

Conditional DDIM generated realistic expression profiles. Tissue-specific consensus ranking and light synthetic regularization improved held-out prediction in selected tissues. Synthetic-informed selection prioritized promoted and reinforced genes, while literature annotation separated prior alignment from complementary hypotheses. Biological evidence and FDR remained grounded in observed OSDR profiles.

## 22. Thank you (0:10)

Acknowledge the mentor, SLSTP, NASA OSDR, ARCHS4, Reactome and NASA Ames compute. Invite questions.
