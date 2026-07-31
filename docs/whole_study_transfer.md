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
