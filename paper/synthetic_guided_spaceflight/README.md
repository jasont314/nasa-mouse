# Synthetic-guided spaceflight paper package

This package contains a research manuscript on a configurable bulk RNA-seq
generation framework and its use in organ-specific mouse spaceflight analysis.
The paper first describes the expression, feature, harmonization, cohort,
conditioning, training, and model choices represented by the framework. It then
reports the WGAN-GP, DDIM, and nine-arm liver harmonization screens before
introducing the selected diffusion model and the
downstream biology.

The factorized DDIM was used downstream because it combined high fidelity,
near-chance real-versus-synthetic adversarial accuracy, and lower distributional
distance. WGAN-GP retained strong correlation and neighborhood metrics but
remained more distinguishable from real profiles. Metrics are reported on each
model's stated evaluation split rather than as a paired comparison.

The downstream analysis first compares real-only, synthetic-only, and
real-plus-synthetic classifiers using the same 974 genes and classifier
settings. It then uses consensus ranking as a secondary analysis of compact,
correlated gene panels. Thymus has the strongest support across both analyses.
Liver, skin, and spleen provide narrower matched findings, while soleus provides
a coherent secondary consensus panel.

Targeted literature review covers all 49 associations from the secondary
consensus analysis, all 21 matched gene associations, and all ten eligible
grouped Reactome associations. Selection behavior is kept separate from the
aligning, complementary, ambiguous, or unmatched literature label. Evidence
scope distinguishes exact matches, process-level agreement, and mechanistic
context.

Using the DDIM generator, the paper follows a four-part downstream funnel:

1. a pooled benchmark tests whether one synthetic-data policy works across
   tissues;
2. matched all-gene classifiers isolate the effect of training source;
3. tissue-specific consensus ranking identifies compact gene panels;
4. real-only random-effects testing evaluates prioritized genes.

The pooled augmentation benchmark was negative. In matched tissue-specific
analysis, real-plus-synthetic training passed all pooled and accession-macro
metric checks in 18 of 27 analysis units and improved at least one metric in 16.
Twenty-one BH-FDR associations also had synthetic-supported marginal
importance. All BH-FDR effects use real OSDR samples only.

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
  Tables S16-S17 contain the consensus-gene literature annotations and source
  inventory. Tables S18-S21 contain the matched classifier results. Tables
  S22-S24 contain matched-gene and grouped-pathway annotations and their source
  inventory.

## Rebuild

The build only reads completed outputs. It does not retrain any model or rerun
feature selection.

```bash
python \
  -m nasa_mouse_diffusion.paper_parity.build_synthetic_guided_paper
```

The builder checks key frozen values before writing figures, then renders the
Markdown manuscripts with WeasyPrint.

## Editorial provenance

The current manuscript, presentation text, and speaker notes were audited on
2026-08-03 with
[`blader/humanizer`](https://github.com/blader/humanizer) version 2.9.1 under
its no-fabrication rule. The file-mode pass simplified formulaic transitions and
repetitive sentence structures while preserving the scientific register. It did
not alter frozen analysis outputs, source-table values, gene names, citations,
or statistical interpretations. Humanizer was an editorial check, not part of
the scientific method. This audit is recorded in the commit named
`Reaudit manuscript and slides with Humanizer`; it supersedes an earlier
provenance note whose cited commit is not present in this repository.

## Evidence interpretation

- **Matched all-gene evidence:** 21 real-data BH-FDR associations had measurable
  synthetic-supported importance in a fixed 974-gene classifier. This is the
  primary analysis of synthetic contribution.
- **Consensus panel evidence:** 49 real-data BH-FDR associations were repeatedly
  promoted or reinforced in compact panels. This is the secondary analysis for
  correlated genes and pathways.
- **Complete real-data screen:** all 459 random-effects BH-FDR tissue-gene
  results, regardless of synthetic feature-selection status.
- **Literature interpretation:** all 49 consensus associations, 21 matched gene
  associations, and ten eligible grouped pathways receive one mutually
  exclusive literature label. The tables distinguish direct matches from
  process-level agreement, mechanistic context, and unmatched candidates.

## Tissue evidence labels

- **Strongest joint result:** thymus combines matched gene importance, a coherent
  cell-cycle program, and consensus support.
- **Narrower matched findings:** liver, skin, and spleen have retained genes but
  less pathway coherence.
- **Secondary consensus findings:** soleus provides a coherent metabolic panel;
  pooled muscle, kidney, adrenal gland, gastrocnemius, and tibialis anterior have
  narrower panels that did not pass the matched gene gate.
