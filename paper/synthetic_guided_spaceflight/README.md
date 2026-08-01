# Synthetic-guided spaceflight paper package

This package contains a research manuscript on a configurable bulk RNA-seq
generation framework and its use in organ-specific mouse spaceflight analysis.
The paper first describes the expression, feature, harmonization, cohort,
conditioning, training, and model choices represented by the framework. It then
reports the WGAN-GP, DDIM, GeneJEPA representation, and nine-arm liver
harmonization screens before introducing the selected diffusion model and the
downstream biology.

The factorized DDIM was used downstream because it combined high fidelity,
near-chance real-versus-synthetic adversarial accuracy, and lower distributional
distance. WGAN-GP retained strong correlation and neighborhood metrics but
remained more distinguishable from real profiles. Metrics are reported on each
model's stated evaluation split rather than as a paired comparison.

The downstream analysis compares pooled and tissue-specific uses of generated
expression. Thymus provides a real-supported cell-cycle panel, soleus provides a
cross-accession metabolic program, and kidney provides a focused promoted
`Inpp4b` and reinforced `Slc37a4` result.

Using the DDIM generator, the paper follows a three-part downstream funnel:

1. a pooled benchmark tests whether one synthetic-data policy works across
   tissues;
2. tissue-specific development compares five synthetic-data uses;
3. real-only random-effects testing evaluates prioritized genes.

The pooled augmentation benchmark was negative, while several tissues improved
during tissue-specific development. All BH-FDR effects use real OSDR samples
only. Synthetic attribution is retained only where the selected generated arm
passed the balanced-accuracy, AUROC, and average-precision eligibility gate.

Generated profiles are not treated as additional animals. Biological support is
calculated from real NASA OSDR profiles using within-accession effects,
random-effects meta-analysis, and Benjamini-Hochberg FDR. BH FDR below 0.05 is
the primary gene-level inclusion rule. Synthetic-selection status, study-direction
agreement, heterogeneity, and leave-one-accession-out results are reported as
separate interpretation or sensitivity annotations rather than inclusion gates.

## Contents

- `manuscript.md`, `manuscript.html`, and `manuscript.pdf`: full paper draft.
- `supplementary_methods.md`, `.html`, and `.pdf`: implementation details,
  evaluation metrics, complete supporting results, and limitations.
- `figures/`: five main figures and two supplementary figures in PNG/PDF.
  Figures 1-2 report generator metrics and diffusion trajectories; Figures 3-5
  report the tissue biology.
- `source_data/`: manuscript and supplementary data tables. SHA-256 manifests
  are retained for repository auditing but are not part of the formal supplement.

## Rebuild

The build only reads completed outputs. It does not retrain any model or rerun
feature selection.

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m nasa_mouse_rna_diffusion.build_synthetic_guided_paper
```

The builder checks key frozen values before writing figures, then renders the
Markdown manuscripts with WeasyPrint.

## Editorial provenance

The manuscript prose was audited with `blader/humanizer` version 2.9.1 under its
no-fabrication rule. The pass removed formulaic transitions and repetitive
sentence structures while preserving the scientific register. It did not alter
frozen analysis outputs, source-table values, gene names, citations, or
statistical interpretations. Humanizer was an editorial check, not part of the
scientific method.

## Evidence interpretation

- **Synthetic-informed development:** real-data BH-FDR genes were also repeatedly
  synthetic-promoted or reinforced. This contains 49 tissue-gene results and
  supports developmental hypotheses rather than independent biological evidence.
- **Complete real-data screen:** all 459 random-effects BH-FDR tissue-gene
  results, regardless of synthetic feature-selection status.

## Tissue evidence labels

- **Cross-accession development:** selected BH-FDR genes form a coherent pattern
  across represented studies, but the generator was adapted in the same
  development domain. This applies most clearly to soleus, pooled skeletal
  muscle, and kidney.
- **Exploratory:** predictive or enrichment signals without a stable real
  accession-aware gene set. This applies to
  spleen, skin, adrenal gland, gastrocnemius, and tibialis anterior.
- **Negative:** no coherent synthetic-guided biological result under the declared
  rules. Quadriceps, EDL, and liver retained real-only arms. Lung improved during
  development but had no BH-FDR gene in the 974-gene panel.
