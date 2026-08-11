# Biological review of corrected kidney and spleen models

Every pathway in the primary top decile was manually reviewed as a complete Reactome program. Magnitude, direction across projects and seeds, conventional enrichment, composition sensitivity, tissue fit, redundancy, and prior literature were considered separately. A favorable numerical score alone was not enough to promote a pathway.

## Kidney

The coherent kidney result is a higher structural and growth-factor-response axis: ECM proteoglycans, aggregate WNT signaling, and IGF transport are higher in flight. These programs agree with prior reports of renal ECM dysregulation, fibrosis-related signaling, nephron remodeling, and WNT involvement, while their joint expression adds a pathway-level repair or maladaptive-growth hypothesis. All three have strong member-gene directional support, but none passes conventional preranked-GSEA FDR below 0.05 and each is attenuated by broad composition-proxy adjustment, so the axis remains a corroborating or complementary hypothesis rather than a statistically confirmed pathway discovery. Lower amino-acid-metabolism and biological-oxidation latent scores are not promoted as broad decreases. Amino-acid metabolism is rejected because raw genes, ssGSEA, preranked GSEA, and prior kidney GSEA point higher. Biological oxidation is retained only as a heterogeneous enzyme subset because several influential CYP, UGT, and ACSM genes are lower while the full gene set points higher. Platelet degranulation is numerically robust but remains a vascular or blood-composition signal rather than a kidney-cell claim.

## Spleen

The coherent spleen result is lower adaptive activation plus lower innate pathogen-sensing and degranulation programs. T-cell receptor signaling, neutrophil degranulation, and C-type lectin receptor signaling are lower in all five unconfounded projects, all three seeds, both conventional methods, and composition-adjusted analyses; each also passes conventional preranked-GSEA FDR below 0.05 and has majority member-gene support. This recovers prior T-cell activation deficits and adds a coordinated innate-effector layer. The transcriptomic degranulation score does not directly measure neutrophil function. CD28-linked activation and lymphoid-to-nonlymphoid interaction scores are lower but disagree with preranked GSEA, while MHC class II antigen presentation is directionally mixed across projects; these remain supporting context rather than central claims.

## Scope recommendation

Spleen should no longer be described as a junk or null tissue. Its corrected HVG reference-query model has a coherent immune-suppression pattern across five unconfounded projects and is a stronger manuscript candidate than soleus. Kidney also has a coherent three-program structural and growth-factor axis, but it is composition-sensitive and lacks conventional pathway FDR support; it is suitable as a secondary corroborating and pathway-integration result with conservative language, not as an equal-strength discovery claim. Lung and retina remain the actual null tissues from the original screen.

## Main-candidate member-gene direction

### Kidney

- **Signaling By Wnt:** 84% of measured member genes move with the pathway score. Across the three trained decoders, a median 72% of absolute decoder weight predicts the observed member-gene direction (minimum 67%). Largest concordant seed-2020 effects: Wnt7a (raw +0.04, weight +0.93); Wnt7b (raw +0.22, weight +0.56); Fzd2 (raw +0.13, weight +0.55); Wnt4 (raw +0.25, weight +0.38); Sox4 (raw +0.08, weight +0.36); Rspo3 (raw +0.17, weight +0.34); Wnt5a (raw +0.06, weight +0.29); Fzd4 (raw +0.02, weight +0.23).
- **Ecm Proteoglycans:** 86% of measured member genes move with the pathway score. Across the three trained decoders, a median 80% of absolute decoder weight predicts the observed member-gene direction (minimum 77%). Largest concordant seed-2020 effects: Tnxb (raw +0.14, weight +0.58); Bgn (raw +0.29, weight +0.44); Col2a1 (raw -0.02, weight -0.42); Tnn (raw +0.38, weight +0.41); Col5a3 (raw +0.33, weight +0.35); Col6a3 (raw +0.18, weight +0.34); Dcn (raw +0.03, weight +0.34); Col6a1 (raw +0.24, weight +0.32).
- **Regulation Of Insulin Like Growth Factor Igf Transport And Uptake By Insulin Like Growth Factor Binding Proteins Igfbps:** 81% of measured member genes move with the pathway score. Across the three trained decoders, a median 63% of absolute decoder weight predicts the observed member-gene direction (minimum 60%). Largest concordant seed-2020 effects: Chgb (raw +0.06, weight +0.67); Igfbp2 (raw +0.31, weight +0.66); Spp2 (raw -0.03, weight -0.56); Fstl1 (raw +0.13, weight +0.34); Sparcl1 (raw +0.23, weight +0.34); Vgf (raw +0.00, weight +0.33); Lgals1 (raw +0.12, weight +0.32); Rcn1 (raw +0.02, weight +0.31).

### Spleen

- **Tcr Signaling:** 62% of measured member genes move with the pathway score. Across the three trained decoders, a median 69% of absolute decoder weight predicts the observed member-gene direction (minimum 62%). Largest concordant seed-2020 effects: Trbv16 (raw -0.01, weight +0.33); Cd4 (raw -0.00, weight +0.32); Trat1 (raw -0.18, weight +0.27); Cd247 (raw -0.08, weight +0.25); Cd101 (raw -0.09, weight +0.25); Cd3e (raw -0.00, weight +0.24); Trac (raw -0.17, weight +0.23); Trbc1 (raw -0.20, weight +0.23).
- **C Type Lectin Receptors Clrs:** 62% of measured member genes move with the pathway score. Across the three trained decoders, a median 76% of absolute decoder weight predicts the observed member-gene direction (minimum 61%). Largest concordant seed-2020 effects: Clec4e (raw -0.36, weight +0.67); Clec4d (raw -0.34, weight +0.49); Il1b (raw -0.08, weight +0.39); Cd209a (raw -0.08, weight +0.36); Fcer1g (raw -0.09, weight +0.29); Clec4n (raw -0.07, weight +0.18); Malt1 (raw -0.02, weight +0.08); Map3k14 (raw -0.05, weight +0.06).
- **Neutrophil Degranulation:** 63% of measured member genes move with the pathway score. Across the three trained decoders, a median 78% of absolute decoder weight predicts the observed member-gene direction (minimum 69%). Largest concordant seed-2020 effects: Lyz2 (raw -0.19, weight +0.59); S100a8 (raw -0.74, weight +0.58); S100a9 (raw -0.69, weight +0.57); Camp (raw -0.78, weight +0.55); Ltf (raw -0.81, weight +0.55); Cd177 (raw -0.64, weight +0.54); Hp (raw -0.32, weight +0.51); Elane (raw -0.80, weight +0.50).

## Kidney metabolic-label audit

- **Metabolism Of Amino Acids And Derivatives:** the median latent shift is -0.140, but only 27% of measured genes move in that same direction. The median decoder-weighted gene-direction match is 42% across seeds, and the median decoder-predicted versus observed gene-effect correlation is 0.22.
- **Biological Oxidations:** the median latent shift is -0.110, but only 31% of measured genes move in that same direction. The median decoder-weighted gene-direction match is 55% across seeds, and the median decoder-predicted versus observed gene-effect correlation is 0.31.

## Literature

- **Finch2025:** [Finch et al. Spaceflight causes strain-dependent gene expression changes in the kidneys of mice.](https://doi.org/10.1038/s41526-025-00465-0)
- **Siew2024:** [Siew et al. Cosmic kidney disease: an integrated pan-omic, physiological and morphological study into spaceflight-induced renal dysfunction.](https://doi.org/10.1038/s41467-024-49212-1)
- **Hammond2018:** [Hammond et al. Effects of space flight on mouse liver versus kidney: gene pathway analyses.](https://doi.org/10.3390/ijms19124106)
- **Horie2019:** [Horie et al. Down-regulation of GATA1-dependent erythrocyte-related genes in the spleens of mice exposed to space travel.](https://doi.org/10.1038/s41598-019-44067-9)
- **Gridley2009:** [Gridley et al. Spaceflight effects on T lymphocyte distribution, function and gene expression.](https://doi.org/10.1152/japplphysiol.91126.2008)
- **Martinez2015:** [Martinez et al. Spaceflight and simulated microgravity cause a significant reduction of key gene expression in early T-cell activation.](https://doi.org/10.1152/ajpregu.00449.2014)
- **Hwang2015:** [Hwang et al. Post-spaceflight mouse splenocytes demonstrate altered activation properties and surface molecule expression.](https://doi.org/10.1371/journal.pone.0124380)
- **Wu2024:** [Wu et al. Single-cell analysis identifies conserved features of immune dysfunction in simulated microgravity and spaceflight.](https://doi.org/10.1038/s41467-023-42013-y)
- **Buchheim2026:** [Buchheim et al. Spaceflight alters the immune regulatory functions of neutrophil granulocytes on T lymphocytes.](https://doi.org/10.1016/j.isci.2025.114380)
