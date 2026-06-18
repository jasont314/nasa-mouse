# Liver-Only GLARE Run

- Pretraining: 2,859 Tabula Muris Senis FACS liver cells from 9 mice.
- Fine-tuning: 628 official OSDR liver profiles from 17 studies.
- Shared genes: 21,010.
- Hyperparameters: reused from the previous cross-tissue HPT winners.
- Pretraining best loss: 0.811362 over 50 epochs.
- Fine-tuning best loss: 0.153114 over 50 epochs.
- Total training time: 22m 30s on CPU.
- Final representation: 21,010 genes x 64 latent dimensions.
- KMeans latent silhouette: 0.540910.
- Flight versus ground profiles: 107 versus 103 across 10 paired studies.
- Nominal Wilcoxon cluster tests: 3 of 15 at p < 0.05.
- FDR-significant cluster tests: 0 of 15 at q < 0.05.
- Grouped flight/ground classifier ROC AUC: 0.580981.
- Random-CV flight/ground classifier ROC AUC: 0.663250.

The strongest nominal negative flight shifts are consensus clusters 6, 3, and
9. These results are exploratory because none remains significant after
Benjamini-Hochberg correction.
