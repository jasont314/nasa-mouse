# SLSTP 2026 generative transcriptomics speaker notes

Target length: 12-15 minutes. Planned speaking time: about 13 minutes 30 seconds.

## 1. Synthetic transcriptomics for mouse spaceflight (0:25)

Introduce the question. This talk asks whether generated expression can help analyze tissue-specific FLT versus GC biology without counting synthetic profiles as additional animals.

## 2. Small studies and study effects complicate tissue comparisons (1:00)

OSDR gives broad tissue coverage, but the data are spread across 75 accessions with different mission and assay contexts. ARCHS4 supplies a much larger mouse reference. The challenge is to use that reference without letting study structure masquerade as spaceflight biology.

## 3. We compared three model families and used conditional diffusion (1:10)

The framework kept the published model architectures and varied the surrounding choices. GeneJEPA was useful for representation but had no expression decoder. WGAN-GP and DDIM could generate profiles. The chosen DDIM used TPM, 974 mouse landmarks, ARCHS4 pretraining and OSDR adaptation conditioned on tissue, FLT/GC, accession and material.

## 4. Diffusion learns tissue structure from noise (0:55)

Read the panels from left to right. The same generated profiles begin as noise, develop structure by timestep 200 and approach tissue-conditioned regions at timestep zero. The axes are shared, so the visual change is not caused by rescaling each panel.

## 5. DDIM matched expression and reduced separability (1:15)

Both WGAN-GP and DDIM had high correlation and F1. DDIM had adversarial accuracy near 0.5 and a lower Frechet-distance ratio, so it was harder to separate from real profiles and closer in distribution. The metrics use each model's stated evaluation split, so this is a model-choice summary rather than a paired significance test.

## 6. Synthetic profiles entered the analysis in five different ways (1:15)

Synthetic data can be used for direct training, mixed training or feature guidance. Each tissue could choose among five arms. The eligibility check used held-out real profiles. Once features were nominated, FLT/GC effects and BH FDR were computed from observed OSDR samples only.

## 7. Pooling tissues hid useful signal (1:10)

The simplest pooled augmentation test was negative: balanced accuracy fell from 0.754 to 0.737 with real plus synthetic training. Tissue-specific analysis changed the result. Different tissues benefited from different synthetic uses, which argues against one global augmentation policy.

## 8. Synthetic guidance changed ranking, not statistical evidence (1:05)

Reinforced genes were selected with and without synthetic guidance. Promoted genes crossed the stable-selection rule only with synthetic guidance. Promoted does not mean biologically novel. All 49 synthetic-informed tissue-gene associations also had a supporting effect and BH FDR below 0.05 in real data.

## 9. Thymus points to lower proliferative renewal (1:25)

Thymus produced the clearest promoted panel. The lower mitotic and DNA-replication genes agree with prior reports of thymic involution and altered cell-cycle expression after flight. Because this is bulk RNA-seq, the signal could reflect fewer cycling thymocytes, lower transcription within those cells or both.

## 10. Soleus reinforces a mitochondrial and lipid program (1:15)

Soleus improved with real plus generated training. The selected genes were already stable in real-only analysis, so synthetic data reinforced rather than introduced the panel. Lower Bdh1, Ech1, Bnip3 and Decr1, with higher Tpm1, support altered oxidative metabolism and contractile remodeling.

## 11. Kidney adds a focused signal; other tissues are exploratory (1:00)

Kidney supplied the strongest secondary pair: promoted Inpp4b and reinforced Slc37a4. Spleen, skin and adrenal results are narrower. Lung and retina improved predictively without a synthetic-informed BH-FDR gene, while liver, EDL and quadriceps retained real-only models.

## 12. Synthetic data helped most as a tissue-specific prior (0:55)

The main result is not that synthetic samples increase biological n. The useful role was tissue-specific feature ranking or limited regularization. The next decisive test is complete-accession holdout followed by validation of the thymus, soleus and kidney hypotheses in independent data.

## 13. Thank you (0:10)

Acknowledge the mentor, SLSTP, NASA OSDR, ARCHS4, Reactome and NASA Ames compute. Invite questions.
