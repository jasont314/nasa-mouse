# Synthetic-guided spaceflight paper package

This package contains a research manuscript on a configurable bulk RNA-seq
generation framework and its use in organ-specific mouse spaceflight analysis.
The paper first describes the expression, feature, harmonization, cohort,
conditioning, training, and model choices represented by the framework. It then
reports the WGAN-GP, DDIM, GeneJEPA representation, and nine-arm liver
harmonization screens before introducing the selected diffusion model and the
downstream biology.

The factorized DDIM was selected because it was the only generator to pass the
final joint locked fidelity and effect-recovery gates. The WGAN-GP retained
strong validation correlation and neighborhood metrics but failed external
indistinguishability and accession-aware effect recovery, so its locked test was
not opened. These rows are consecutive selection stages, not a paired test-set
comparison.

The central biological result is the study-excluded thymus test. OSD-457 was
excluded from OSDR adaptation, and all OSDR-linked GEO series were excluded from
ARCHS4 pretraining. Soleus provides cross-accession metabolic evidence, while
kidney provides a focused promoted `Inpp4b` and reinforced `Slc37a4` result.

After generator selection, the paper follows a four-stage downstream funnel:

1. a pooled benchmark tests whether one synthetic-data policy works across
   tissues;
2. tissue-specific development compares five synthetic-data uses;
3. real-only random-effects testing evaluates prioritized genes;
4. complete-accession tests evaluate transfer outside the development studies.

The pooled augmentation benchmark was negative. Several tissues improved during
tissue-specific development, but most gains did not establish whole-study
transfer. Thymus supplied the strongest retained transfer result. An initial
pooled-muscle augmentation gain did not generalize when the same frozen recipe
was extended to 11 held-out accessions.

The same pretrained backbone was used for the lung/thymus held-out experiment
and for regenerated all-tissue and muscle-group development screens. All BH-FDR
effects use real OSDR samples only. Synthetic attribution is retained only where
the selected generated arm passed the balanced-accuracy, AUROC, and
average-precision eligibility gate.

Generated profiles are not treated as additional animals. Biological support is
calculated from real NASA OSDR profiles using within-accession effects,
random-effects meta-analysis, and Benjamini-Hochberg FDR. BH FDR below 0.05 is
the primary gene-level inclusion rule. Synthetic-selection status, study-direction
agreement, heterogeneity, and leave-one-accession-out results are reported as
separate interpretation or sensitivity annotations rather than inclusion gates.

## Contents

- `manuscript.md`, `manuscript.html`, and `manuscript.pdf`: full paper draft.
- `supplementary_methods.md`, `.html`, and `.pdf`: exact implementation,
  evaluation gates, output provenance, and limitations.
- `figures/`: six main figures and four supplementary figures in PNG/PDF.
  Figures 1-3 describe the configurable pipeline, generator selection, and
  diffusion trajectories; Figures 4-6 report the tissue biology.
- `source_data/`: manuscript tables, figure source tables, and SHA-256 manifests
  for every frozen analysis input.

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

The final manuscript prose was audited with `blader/humanizer` version 2.9.1 at
repository commit `523374dee72d67c7b2b5f858ea0094ffda49c3ac`. The skill was
applied in file mode under its no-fabrication rule. The pass removed formulaic
transitions and repetitive sentence structures while preserving the scientific
register. It did not alter frozen analysis outputs, source-table values, gene
names, citations, or statistical interpretations. Humanizer was an editorial
check, not part of the scientific method.

## Evidence interpretation

- **Whole-study transfer:** a fixed generator and feature policy were tested in
  an OSDR accession excluded from adaptation and policy development after
  removing every OSDR-linked GEO series from ARCHS4. The eight-gene OSD-457
  thymus panel is the central result. A new study is still required for
  prospective replication.
- **Synthetic-informed development:** real-data BH-FDR genes were also repeatedly
  synthetic-promoted or reinforced. This contains 49 tissue-gene results and
  supports developmental hypotheses rather than independent transfer.
- **Complete real-data screen:** all 459 random-effects BH-FDR tissue-gene
  results, regardless of synthetic feature-selection status.

## Tissue evidence labels

- **Cross-accession development:** selected BH-FDR genes form a coherent pattern
  across represented studies, but the generator was adapted in the same
  development domain. This applies most clearly to soleus, pooled skeletal
  muscle, and kidney.
- **Exploratory:** predictive or enrichment signals without a stable real
  accession-aware gene set or without unseen-study transfer. This applies to
  spleen, skin, adrenal gland, gastrocnemius, and tibialis anterior.
- **Negative:** no coherent synthetic-guided biological result under the declared
  rules. Quadriceps, EDL, and liver retained real-only arms; lung failed its
  accession-held-out test.
