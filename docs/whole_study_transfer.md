# Twelve-tissue whole-study transfer

## Purpose

This analysis evaluates whether DDIM-derived feature guidance transfers to real
OSDR studies that were absent from generator adaptation, feature selection, and
classifier fitting. It uses all canonical tissues with at least three eligible
FLT/GC accessions:

- adrenal gland
- brain
- cerebellum
- heart
- kidney
- liver
- lung
- retina
- skeletal muscle
- skin
- spleen
- thymus

The OSDR input remains the NASA OSDR API-derived expression object. No raw
integrated OSDR H5 file is used.

## Split contract

Eligible accessions are assigned to three global outer folds. Every accession is
in the test role exactly once. Within a fold, validation and test accessions are
disjoint, and every tissue retains FLT and GC profiles in train, validation, and
test. An accession has one role globally even when it contributes multiple
tissues.

Each fold starts from the 15,000-epoch ARCHS4 mouse checkpoint from which all
GEO series linked to eligible OSDR accessions were removed. The OSDR adapter is
then trained for 5,000 epochs using only that fold's training accessions.

The three model configurations are:

```text
configs/rna_diffusion/osdr_whole_study_transfer_12t_fold0.yaml
configs/rna_diffusion/osdr_whole_study_transfer_12t_fold1.yaml
configs/rna_diffusion/osdr_whole_study_transfer_12t_fold2.yaml
```

The common downstream configuration is:

```text
configs/rna_diffusion/generated_feature_whole_study_transfer_12t.yaml
```

## Evaluation levels

Pooled, tissue-specific, and tissue-by-accession effect recovery are not
interchangeable.

**Pooled FLT/GC recovery** combines every test tissue in one fold, computes one
974-gene real FLT-minus-GC vector and one synthetic vector, and compares them.
It is a broad check that the conditional label is not ignored or inverted. It
can be dominated by abundant tissues and can hide opposing tissue responses.

**Per-tissue recovery** estimates real and synthetic FLT-minus-GC effects from
the held-out accessions separately for each tissue. It asks whether the generator
preserves a tissue's average response without allowing other tissues to dominate.

**Per-tissue, per-accession recovery** compares the 974-gene effect vectors
inside each held-out study. It asks whether the generated condition shift agrees
with the local study effect and reveals heterogeneity hidden by either pooled or
tissue-level averaging.

The downstream classifier analysis independently compares real-only and
synthetic-guided policies using balanced accuracy, AUROC, and average precision
on held-out real profiles. Generated profiles are not biological replicates and
are never included in random-effects P values or BH-FDR calculations.

## Results

The three folds evaluated 68 tissue-accession pairs from 63 unique accessions
and 1,284 profiles. Five accessions contributed more than one tissue. Every pair
appeared in the test role once.

Across accessions, the real-only baseline had balanced accuracy 0.581, AUROC
0.620, and average precision 0.694. Validation-gated deployment changed these to
0.583, 0.631, and 0.708. The gain was therefore small overall. Skin, pooled
skeletal muscle, and lung met the predeclared tissue rule: the synthetic-informed
candidate passed the inner validation gate in at least two folds, aggregate
balanced accuracy increased, and neither aggregate AUROC nor average precision
decreased.

| Tissue | Accessions | BA change | AUROC change | AP change | Rule passed |
|---|---:|---:|---:|---:|---|
| Adrenal gland | 3 | 0.000 | 0.000 | 0.000 | No |
| Brain | 4 | 0.024 | -0.104 | -0.092 | No |
| Cerebellum | 3 | 0.000 | -0.041 | -0.011 | No |
| Heart | 3 | -0.033 | -0.049 | -0.036 | No |
| Kidney | 6 | 0.000 | 0.030 | 0.027 | No |
| Liver | 12 | -0.004 | -0.074 | -0.037 | No |
| Lung | 3 | 0.056 | 0.164 | 0.080 | Yes |
| Retina | 4 | -0.014 | 0.062 | 0.032 | No |
| Skeletal muscle | 13 | 0.015 | 0.055 | 0.060 | Yes |
| Skin | 6 | 0.077 | 0.110 | 0.079 | Yes |
| Spleen | 6 | -0.040 | 0.019 | 0.013 | No |
| Thymus | 5 | -0.060 | -0.017 | 0.009 | No |

The three effect summaries gave different answers:

- Pooled fold-level effect correlation averaged 0.232 and direction agreement
  averaged 0.575. The three fold correlations were 0.380, 0.158, and 0.160.
- Per-tissue random-effects correlation averaged 0.489 and direction agreement
  averaged 0.682. Thymus (0.776), skin (0.733), and pooled skeletal muscle
  (0.598) had the highest tissue-level correlations.
- Direct per-tissue, per-accession correlation averaged 0.033 and direction
  agreement averaged 0.468 across 68 pairs.

The tissue-level result shows that the generator captures some average
tissue-specific response structure after effects are combined across studies.
The near-zero accession-level result shows that it does not reproduce the local
FLT/GC effect of an arbitrary unseen study. The pooled number is not a substitute
for either result: tissue abundance and opposing organ responses can change it.
The field name `delta_correlation` in the pooled output means correlation between
the real and synthetic FLT-minus-GC vectors; it is not a performance difference
between two classifiers.

## Execution

For each fold, prepare, train, evaluate the held-out test role, and generate the
fixed scale-2 training draws. Fold-specific seeds are declared in the common
configuration.

```bash
python -m nasa_mouse_rna_diffusion prepare-osdr --config <fold-config>
python -m nasa_mouse_rna_diffusion train-osdr --config <fold-config>
python -m nasa_mouse_rna_diffusion evaluate-osdr --config <fold-config> --unlock-test --evaluation-variant whole_study
python -m nasa_mouse_rna_diffusion generate-contrastive-osdr --config <fold-config> --guidance-scales 2 --seeds <fold-seeds>
```

After all three folds:

```bash
python -m nasa_mouse_rna_diffusion whole-study-transfer \
  --config configs/rna_diffusion/generated_feature_whole_study_transfer_12t.yaml
```

Outputs are written under:

```text
outputs/generative_benchmark/analyses/whole_study_transfer_12_tissue_v1/
```

This is a uniform retrospective cross-study evaluation. It is stronger than
within-accession profile splitting but is not a prospective validation on a new
mission collected after the analysis was specified.
