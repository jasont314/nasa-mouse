# SLSTP 2026 generative transcriptomics speaker notes

Target length: 12-15 minutes. Planned speaking time: about 14 minutes 30 seconds.

## 1. Synthetic transcriptomics for mouse spaceflight (0:25)

Introduce the question. This talk asks whether generated expression can help analyze tissue-specific FLT versus GC biology without counting synthetic profiles as additional animals.

## 2. Small studies and study effects complicate tissue comparisons (0:50)

OSDR gives broad tissue coverage, but the data are spread across 75 accessions with different mission and assay contexts. ARCHS4 supplies a much larger mouse reference. The challenge is to use that reference without letting study structure masquerade as spaceflight biology.

## 3. We built a configurable bulk RNA-seq generation pipeline (1:10)

Each column shows the alternatives available at one pipeline stage; outlines identify the downstream branch. We used both ARCHS4 and NASA OSDR across multiple studies and all tissues. The selected path used TPM, training-fitted MaxAbs scaling, 974 mouse landmarks, no global correction, ARCHS4 pretraining plus OSDR adaptation, and a DDIM conditioned on tissue, FLT/GC, accession and material. WGAN-GP and the other preprocessing and harmonization choices remained benchmark alternatives.

## 4. DDIM matched expression and reduced separability (0:55)

The left image is Figure 1C from Lacan and colleagues. It shows their residual dense denoiser, conditioned on diffusion timestep and tissue, which predicts the noise added to a sample. Our implementation follows that architecture and adds FLT/GC, accession and material context during OSDR LoRA adaptation. Both generators had high correlation and F1. DDIM had adversarial accuracy near 0.5 and a lower Frechet-distance ratio, so it was harder to separate from real profiles and closer in distribution.

## 5. Diffusion learns tissue structure from noise (0:50)

Read the panels from left to right. The same generated profiles begin as noise, develop structure by timestep 200 and approach tissue-conditioned regions at timestep zero. The axes are shared, so the visual change is not caused by rescaling each panel.

## 6. Generated profiles track the real OSDR PCA manifold (0:50)

Circles are locked real OSDR profiles and crosses are matched DDIM profiles in the same PCA space. Generated samples track the tissue-defined branches. FLT and GC remain more interspersed because condition effects are smaller than tissue effects. Visual overlap complements the quantitative validation.

## 7. Synthetic profiles entered the analysis in five different ways (0:55)

Synthetic data can be used for direct training, mixed training or feature guidance. Each tissue could choose among five arms. The eligibility check used held-out real profiles. Once features were nominated, FLT/GC effects and BH FDR were computed from observed OSDR samples only.

## 8. Pooling tissues hid useful signal (0:55)

The pooled augmentation test was negative: balanced accuracy fell from 0.754 to 0.737 with real plus synthetic training. The bars use a true zero-to-one balanced-accuracy scale. Tissue-specific analysis changed the result. Different tissues benefited from different synthetic uses, which argues against one global augmentation policy.

## 9. Synthetic guidance changed ranking, not statistical evidence (0:50)

Reinforced genes were selected with and without synthetic guidance. Promoted genes crossed the stable-selection rule only with synthetic guidance. Promoted does not mean biologically novel. All 49 synthetic-informed tissue-gene associations also had BH FDR below 0.05 in real data.

## 10. Selection and literature are separate dimensions (0:55)

Every association has two labels. Promoted or reinforced describes repeated feature selection. Aligning, complementary, ambiguous or unmatched describes prior literature. Across all 49 associations, 22 aligned, 19 were complementary, four were ambiguous and four were unmatched. Table S16 records the gene-level rationale, evidence scope and source IDs; Table S17 records citations, DOI or URL and whether the evidence is independent, overlapping or mechanistic context.

## 11. The screen covered all 27 completed tissue analyses (0:45)

This is the full analysis coverage: 22 canonical tissues and five anatomical muscle groups. Ten had a synthetic-informed BH-FDR association, five had real BH-FDR genes without synthetic-informed selection, and 12 had no BH-FDR gene in the landmark panel. The discussion focuses later, but every completed tissue result remains visible here.

## 12. Ten tissue analyses contained synthetic-informed genes (0:45)

This is the complete 49-association inventory. Separate rows show FLT-higher or FLT-lower direction and promoted or reinforced selection status. Gene color independently shows aligning, complementary, ambiguous or unmatched literature. FLT directions come from real-data meta-analysis.

## 13. Thymus points to lower proliferative renewal (1:10)

Thymus produced the clearest promoted panel. The lower mitotic and DNA-replication genes agree with prior reports of thymic involution and altered cell-cycle expression after flight. Higher Hsd17b11 and Etv1 add lipid-handling and T-cell-state hypotheses. Because this is bulk RNA-seq, the pattern may reflect transcription, cell composition or both.

## 14. Soleus reinforces a mitochondrial and lipid program (1:00)

Soleus improved with real plus generated training. The selected genes were already stable in real-only analysis, so synthetic data reinforced rather than introduced the panel. Lower Bdh1, Ech1, Bnip3 and Decr1, with higher Tpm1, support altered oxidative metabolism and contractile remodeling. The literature is mixed for Bnip3 and Tpm1, which is recorded explicitly.

## 15. Additional tissues produced distinct hypotheses (0:40)

Promoted and reinforced genes are shown on separate subrows for each tissue. Pooled muscle, kidney, spleen and skin each produced a distinct synthetic-informed result. The rows share a slide for presentation space; each remains a separate hypothesis. Pooled muscle is heterogeneous, kidney suggests phosphoinositide and glucose handling, spleen suggests adhesion and extracellular-matrix or immune organization, and skin contributes a single interferon-linked candidate.

## 16. Eye, adrenal and muscle-group results remain tissue-specific (0:40)

Promoted and reinforced genes remain separated here as well. Eye reinforces lower cytokinesis, adrenal contributes two unmatched candidates, gastrocnemius combines an NF-kappa-B stress signal with an autophagy or myogenesis candidate, and tibialis anterior spans stress, cell-cycle, ganglioside and mitophagy hypotheses.

## 17. Synthetic data worked best as a tissue-specific prior (0:40)

Synthetic data was useful for tissue-specific feature ranking and limited regularization. It did not increase biological sample size. Literature annotation separated exact recovery, process-level agreement and complementary hypotheses. Independent samples and cell-resolved experiments are the next tests.

## 18. Thank you (0:10)

Acknowledge the mentor, SLSTP, NASA OSDR, ARCHS4, Reactome and NASA Ames compute. Invite questions.
