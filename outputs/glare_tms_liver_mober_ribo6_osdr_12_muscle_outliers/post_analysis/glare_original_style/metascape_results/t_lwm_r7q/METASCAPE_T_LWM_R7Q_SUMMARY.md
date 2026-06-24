# GLARE Direct-Cluster Metascape Summary

- Session: https://metascape.org/gp/index.html#/reportfinal/t_lwm_r7q
- Input: direct FLT and GC consensus cluster gene lists from the GLARE original-style analysis.
- Lists submitted: 29
- Converted input genes: 16,572
- Custom background: 20,803 genes (20,808 recognized from 21,010 submitted).
- Metascape settings: Mus musculus input and analysis species, p < 0.01, minimum overlap 3, enrichment factor > 1.5, custom background enabled.
- Note: Metascape did not provide `Enrichment_PPI/_FINAL_MCODE.csv` for this session; GO, heatmap, report HTML, and `GO_MCODE.csv` downloaded successfully.

## XGBoost Verification

| scheme | accuracy mean | F1 mean | ROC AUC mean |
|---|---:|---:|---:|
| gene_grouped_kfold_audit | 0.891 | 0.879 | 0.976 |
| random_kfold_glare | 0.888 | 0.884 | 0.976 |

## Top Shared Enrichment Terms

| term | description | best log10(P) | significant lists | strongest lists |
|---|---|---:|---:|---|
| GO:0044282 | small molecule catabolic process | -62.08 | 6 | GC_cluster_12:-62.08; FLT_cluster_11:-60.77; GC_cluster_11:-9.39; FLT_cluster_05:-8.04; GC_cluster_09:-4.10; FLT_cluster_10:-2.04 |
| GO:0032787 | monocarboxylic acid metabolic process | -56.66 | 8 | FLT_cluster_11:-56.66; GC_cluster_12:-53.40; FLT_cluster_14:-13.44; GC_cluster_09:-9.02; FLT_cluster_05:-6.93; GC_cluster_11:-6.33 |
| GO:0044283 | small molecule biosynthetic process | -42.69 | 6 | FLT_cluster_11:-42.69; GC_cluster_12:-39.27; FLT_cluster_05:-5.38; FLT_cluster_14:-5.26; GC_cluster_11:-4.16; GC_cluster_09:-3.13 |
| GO:0008202 | steroid metabolic process | -35.34 | 7 | GC_cluster_12:-35.34; FLT_cluster_11:-32.59; FLT_cluster_14:-14.36; FLT_cluster_10:-3.55; GC_cluster_09:-2.83; GC_cluster_11:-2.41 |
| mmu04610 | Complement and coagulation cascades - Mus musculus (house mouse) | -31.94 | 3 | GC_cluster_12:-31.94; FLT_cluster_11:-31.43; FLT_cluster_09:-2.33 |
| mmu05012 | Parkinson disease - Mus musculus (house mouse) | -31.25 | 6 | FLT_cluster_05:-31.25; GC_cluster_11:-26.26; GC_cluster_09:-19.53; FLT_cluster_11:-9.28; FLT_cluster_10:-4.43; GC_cluster_12:-3.16 |
| mmu04820 | Cytoskeleton in muscle cells - Mus musculus (house mouse) | -30.51 | 4 | FLT_cluster_13:-30.51; FLT_cluster_15:-5.82; GC_cluster_07:-4.46; GC_cluster_01:-3.07 |
| GO:0009410 | response to xenobiotic stimulus | -29.70 | 4 | FLT_cluster_11:-29.70; GC_cluster_12:-25.90; FLT_cluster_14:-10.62; GC_cluster_08:-3.29 |
| GO:0006790 | sulfur compound metabolic process | -29.53 | 7 | GC_cluster_12:-29.53; FLT_cluster_11:-23.59; FLT_cluster_14:-9.87; FLT_cluster_05:-5.98; GC_cluster_09:-3.61; FLT_cluster_10:-2.52 |
| GO:0072521 | purine-containing compound metabolic process | -27.35 | 7 | FLT_cluster_11:-27.35; GC_cluster_12:-25.26; FLT_cluster_05:-11.41; GC_cluster_09:-8.87; GC_cluster_11:-8.19; FLT_cluster_14:-3.45 |
| R-MMU-72766 | Translation | -26.62 | 11 | FLT_cluster_05:-26.62; GC_cluster_11:-20.13; GC_cluster_09:-18.17; FLT_cluster_11:-14.22; FLT_cluster_14:-6.78; GC_cluster_05:-5.48 |
| mmu05171 | Coronavirus disease - Mus musculus (house mouse) | -25.61 | 6 | FLT_cluster_11:-25.61; FLT_cluster_14:-11.39; GC_cluster_11:-11.18; GC_cluster_12:-7.70; GC_cluster_09:-6.90; FLT_cluster_05:-5.69 |

## DGEA-Overlap Direct Clusters

| cluster | genes | significant DGEA genes | top significant genes | top Metascape terms |
|---|---:|---:|---|---|
| GC_cluster_01 | 2185 | 29 | Myl1,Tnnt3,Acta1,Tg,Myh4,Myh2,Mybpc1,Cox8b,Mucl1,Myh1,Xirp2,Casq2,Trpm5,Sh3bgr,Ptgs2,Cdk3,Trim72,Il6,Fgf10,... | Neuroactive ligand-receptor interaction - Mus musculus (house mouse) (-12.0); metal ion transport (-12.0); monoatomic cation transmembrane transport (-10.0) |
| FLT_cluster_03 | 1968 | 13 | Tg,Mucl1,Trpm5,Ptgs2,Cdk3,Wnt11,Il6,Fgf10,Prtn3,Cxcl2,Gm13490,Mypn,3830403N18Rik | cilium movement (-12.0); metal ion transport (-11.0); cilium movement involved in cell motility (-11.0) |
| FLT_cluster_13 | 326 | 11 | Myl1,Tnnt3,Acta1,Myh4,Myh2,Mybpc1,Myh1,Xirp2,Adipoq,Cyp2c53-ps,Ckmt1 | Cytoskeleton in muscle cells - Mus musculus (house mouse) (-31.0); Striated muscle contraction (-22.0); muscle contraction (-22.0) |
| FLT_cluster_15 | 278 | 8 | Cox8b,Casq2,Sh3bgr,Trim72,Alox15,Myo18b,Hrc,Chrdl2 | muscle structure development (-8.6); muscle cell differentiation (-8.3); muscle tissue development (-7.6) |
| GC_cluster_07 | 1340 | 7 | Treh,Wnt11,Pitx3,Adipoq,Cyp2c53-ps,Trim80,Tuba8 | Regulation of PD-L1(CD274) transcription (-12.0); Condensation of Prophase Chromosomes (-12.0); MLL4 and MLL3 complexes regulate expression of PPARG target genes in adipogenesis and hepatic steatosis (-12.0) |
| GC_cluster_08 | 1173 | 4 | Mup15,Cyp4a14,Mup17,Cdkn1a | Retinol metabolism - Mus musculus (house mouse) (-8.4); Miscellaneous substrates (-7.6); oxidative demethylation (-6.4) |
| FLT_cluster_09 | 913 | 4 | Treh,Pitx3,Trim80,Tuba8 | Regulation of endogenous retroelements by KRAB-ZFP proteins (-19.0); Regulation of endogenous retroelements (-19.0); Condensation of Prophase Chromosomes (-17.0) |
| FLT_cluster_14 | 303 | 4 | Apoa4,Mup15,Mup17,Cdkn1a | xenobiotic metabolic process (-18.0); Cytoplasmic ribosomal proteins (-17.0); cellular response to xenobiotic stimulus (-15.0) |
| GC_cluster_02 | 2104 | 1 | Adcy1 | leukocyte activation (-9.8); Hematopoietic cell lineage - Mus musculus (house mouse) (-8.2); lymphocyte activation (-7.8) |
| FLT_cluster_02 | 2085 | 1 | Sult1e1 | DNA metabolic process (-17.0); tRNA metabolic process (-16.0); DNA repair (-13.0) |
| FLT_cluster_04 | 1619 | 1 | Olfr827 | cell-cell signaling (-8.5); anterograde trans-synaptic signaling (-7.5); chemical synaptic transmission (-7.5) |
| FLT_cluster_06 | 1422 | 1 | Adcy1 | cell projection assembly (-8.9); plasma membrane bounded cell projection assembly (-8.4); cilium organization (-6.1) |

## Interpretation

- The strongest direct-cluster enrichments are liver metabolic programs: small-molecule, monocarboxylic-acid, fatty-acid, steroid, amino-acid, sulfur-compound, xenobiotic, and biological-oxidation terms.
- Muscle/cytoskeleton terms remain visible in the direct clusters even after the 12-profile muscle-outlier filter, especially in direct clusters with `Myl1`, `Tnnt3`, `Acta1`, and myosin-family genes.
- This summary follows the GLARE-style direct FLT/GC cluster interpretation. The earlier `FLT14_not_GC8` paired-cluster split is a separate extension and is not used here.

## Files

- `Enrichment_GO/GO_AllLists.csv`
- `Enrichment_GO/_FINAL_GO.csv`
- `Enrichment_heatmap/HeatmapSelectedGO.csv`
- `top_heatmap_terms.tsv`
- `top_terms_by_list.tsv`
- `AnalysisReport.html`
