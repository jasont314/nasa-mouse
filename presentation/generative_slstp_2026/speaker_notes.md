# SLSTP 2026 generative transcriptomics speaker notes

Target length: 12-15 minutes. Planned speaking time: about 13 minutes 35 seconds.

## 1. Synthetic transcriptomics for mouse spaceflight (0:25)

Introduce the question. This talk asks whether generated expression can help analyze tissue-specific FLT versus GC biology without counting synthetic profiles as additional animals.

## 2. Small studies and study effects complicate tissue comparisons (0:55)

OSDR gives broad tissue coverage, but the data are spread across 75 accessions with different mission and assay contexts. ARCHS4 supplies a much larger mouse reference. The challenge is to use that reference without letting study structure masquerade as spaceflight biology.

## 3. We built a configurable bulk RNA-seq generation pipeline (1:15)

We built one pipeline that can change data source, transformation, feature set, harmonization, model, training scope and conditions without changing the evaluation contract. It supports OSDR-only, ARCHS4-only and ARCHS4-pretrained plus OSDR-adapted runs, with pooled or per-tissue cohorts. We completed WGAN-GP and DDIM generator branches. The selected path used TPM, 974 mouse landmarks, no global correction and a DDIM conditioned on tissue, FLT/GC, accession and material.

## 4. DDIM matched expression and reduced separability (1:00)

Both WGAN-GP and DDIM had high correlation and F1. DDIM had adversarial accuracy near 0.5 and a lower Frechet-distance ratio, so it was harder to separate from real profiles and closer in distribution. The metrics use each model's stated evaluation split, so this is a model-choice summary rather than a paired significance test. These results are why the remaining analyses use DDIM.

## 5. Diffusion learns tissue structure from noise (0:55)

Having selected DDIM, read the panels from left to right. The same generated profiles begin as noise, develop structure by timestep 200 and approach tissue-conditioned regions at timestep zero. The axes are shared, so the visual change is not caused by rescaling each panel.

## 6. Generated profiles track the real OSDR PCA manifold (0:55)

This is PCA, not UMAP. Circles are locked real OSDR profiles and crosses are matched DDIM profiles in the same PCA space. On the left, generated samples track the tissue-defined branches. On the right, FLT and GC remain interspersed because condition effects are smaller than tissue effects. Visual overlap is descriptive and complements the quantitative validation shown earlier.

## 7. Synthetic profiles entered the analysis in five different ways (1:00)

Synthetic data can be used for direct training, mixed training or feature guidance. Each tissue could choose among five arms. The eligibility check used held-out real profiles. Once features were nominated, FLT/GC effects and BH FDR were computed from observed OSDR samples only.

## 8. Pooling tissues hid useful signal (1:00)

The simplest pooled augmentation test was negative: balanced accuracy fell from 0.754 to 0.737 with real plus synthetic training. Tissue-specific analysis changed the result. Different tissues benefited from different synthetic uses, which argues against one global augmentation policy.

## 9. Synthetic guidance changed ranking, not statistical evidence (0:55)

Reinforced genes were selected with and without synthetic guidance. Promoted genes crossed the stable-selection rule only with synthetic guidance. Promoted does not mean biologically novel. All 49 synthetic-informed tissue-gene associations also had a supporting effect and BH FDR below 0.05 in real data.

## 10. Annotation separates recovery from hypothesis extension (1:05)

This slide separates four questions. Promoted tells us that synthetic guidance changed stable feature ranking. BH FDR tells us the association is present in real OSDR profiles. The literature label describes whether prior work is exact, related, mixed or unmatched. The biological interpretation remains a hypothesis. Eleven associations aligned, 13 were complementary, one was ambiguous and one was literature unmatched. Only Ccnb2, Ccne2 and Nfkbia were exact gene-tissue-direction matches. Psmb8 was unmatched in adrenal spaceflight literature but remains mechanistically plausible.

## 11. Thymus points to lower proliferative renewal (1:15)

Thymus produced the clearest promoted panel. The lower mitotic and DNA-replication genes agree with prior reports of thymic involution and altered cell-cycle expression after flight. Higher Hsd17b11 and Etv1 add lipid-handling and T-cell-state hypotheses, but neither is a direct prior flight replication or an established driver. Because this is bulk RNA-seq, the pattern may reflect transcription, cell composition or both.

## 12. Soleus reinforces a mitochondrial and lipid program (1:05)

Soleus improved with real plus generated training. The selected genes were already stable in real-only analysis, so synthetic data reinforced rather than introduced the panel. Lower Bdh1, Ech1, Bnip3 and Decr1, with higher Tpm1, support altered oxidative metabolism and contractile remodeling.

## 13. Kidney adds a focused signal; other tissues are exploratory (0:55)

Kidney supplied the strongest secondary pair: promoted Inpp4b and reinforced Slc37a4. Spleen and skin results are narrower. Adrenal Psmb8 has plausible immunoproteasome biology but no prior adrenal flight match. Lung and retina improved predictively without a synthetic-informed BH-FDR gene, while liver, EDL and quadriceps retained real-only models.

## 14. Synthetic data helped most as a tissue-specific prior (0:45)

The useful role was tissue-specific feature ranking or limited regularization, not increasing biological sample size. Literature annotation then separated exact recovery, process-level agreement and complementary hypotheses. Independent samples and cell-resolved experiments are the next tests.

## 15. Thank you (0:10)

Acknowledge the mentor, SLSTP, NASA OSDR, ARCHS4, Reactome and NASA Ames compute. Invite questions.
