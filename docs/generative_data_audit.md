# Generative Benchmark Data Audit

Audit date: 2026-07-14.

## OSDR API Cohort

The cohort was refreshed from the NASA OSDR Biological Data API with these
filters:

- organism: `Mus musculus`;
- assay: bulk RNA-seq inferred from assay and file metadata;
- data type: unnormalized counts;
- condition: Space Flight or Ground Control;
- source: every OSDR source returned by the query.

The API returned 1,631 profile rows from 75 accessions and 24 canonical material
classes. Twenty-one rows are technical replicates. The default training policy
sums those replicates, leaving 1,610 biological profiles: 835 FLT and 775 GC.

All 75 count tables are cached under `data/osdr_api/counts/`. The query-data
endpoint returned HTTP 500 for `OSD-759` on two attempts. Its count file was
downloaded from the direct NASA URL returned by the official
`/dataset/OSD-759/files/` REST record. No count columns were missing afterward.

The generated API-derived matrix has 1,610 samples and 48,694 Ensembl mouse genes
shared across all 75 count tables and ARCHS4. Values are RSEM unnormalized expected
counts and can be fractional. The older integrated OSDR H5 was not read.

## Tissue Tiers

| Tier | Requirement after technical-replicate collapse | Tissues |
|---|---|---|
| Confirmatory per tissue | at least 60 total, 20 per condition, and 5 paired accessions | liver, skeletal muscle, skin, kidney, thymus, spleen |
| Exploratory pretrained per tissue | at least 30 total, 10 per condition, and 2 paired accessions | lung, retina, brain, cerebellum, colon, heart, adrenal gland, hippocampus |
| Pooled only | below the exploratory threshold | bone, optic nerve, mammary gland, bone marrow, brown adipose, white adipose, eye, cecum, cells, cultured cells |

The confirmatory threshold permits one validation accession, one locked test
accession, and at least three training accessions. The exploratory threshold only
supports leave-one-accession-out reporting. Smaller tissues remain useful in the
pooled tissue-conditioned model but cannot support a defensible standalone model.

The full counts and accession lists are in
`outputs/generative_benchmark/data_audit/osdr/osdr_tissue_inventory.tsv`.

## ARCHS4 Scan

The local `mouse_gene_v2.5.h5` is the complete 997,515-profile ARCHS4 mouse file,
with 53,511 genes. Every profile's metadata was inspected. A profile is eligible
for a reference only when it:

- maps to one of the OSDR canonical tissues;
- is mouse RNA-seq with transcriptomic library source;
- has `singlecellprobability < 0.5` and no explicit scRNA-seq, snRNA-seq, 10x,
  Drop-seq, Chromium, or Smart-seq technology label;
- has no NASA, GeneLab, OSDR, spaceflight, microgravity, ISS, Rodent Research, or
  hindlimb-unloading leakage term.

Three reproducible cohorts were generated:

| Cohort | Samples | GEO series | Purpose |
|---|---:|---:|---|
| `control_only` | 23,614 | 3,213 | high-precision sensitivity arm requiring explicit normal/control-like metadata |
| `healthy_preferred` | 62,299 | 5,307 | primary arm; excludes explicit disease, intervention, genetic perturbation, developmental, and cell-technology terms while retaining unknown health metadata |
| `broad` | 134,250 | 15,111 | diversity sensitivity arm retaining disease and perturbation studies |

All three cover 23 matchable OSDR classes. The generic OSDR `cells` label has no
defensible ARCHS4 tissue match. Each cohort is capped at 10,000 samples per tissue
and 100 per GEO series. Stored sampling weights select tissue uniformly, then GEO
series uniformly within tissue, then samples uniformly within series.

The health labels are metadata heuristics, not clinical verification. That is why
`healthy_preferred` is the primary compromise, `control_only` tests sensitivity to
stricter curation, and `broad` tests whether additional perturbation diversity helps
or teaches the generator unrelated condition effects.

## Why 24,428 Was Not The Full Reference

The earlier WGAN and diffusion runs concatenated already-extracted ARCHS4 references
for eight tissues: liver, skeletal muscle, skin, kidney, thymus, spleen, lung, and
retina. Their combined 24,428 samples were a selected working subset, not the full
997,515-profile ARCHS4 catalog. They omitted all other OSDR tissues and used older
tissue/single-cell filtering.

The new primary reference is 62,299 profiles across every matchable OSDR class.
The strict arm is similar in size to the old reference but has broader tissue
coverage and stronger metadata exclusions.

## Validation Splits

Splits use OSDR accession, never individual samples, and were assigned without
looking at expression values.

- pooled model: 51 training, 12 validation, and 12 locked-test accessions;
- standalone confirmatory models: six tissues with training/validation/locked-test
  accession assignments;
- LOO validation: 72 folds across all 14 confirmatory or exploratory tissues.

The pooled splitter guarantees that every observed tissue-condition pair remains in
training. The final test accessions stay locked until preprocessing and model choices
are fixed.

## Reproduction

```bash
PYTHONPATH=src python -m nasa_mouse_generative osdr-inventory --refresh
PYTHONPATH=src python -m nasa_mouse_generative osdr-expression --download-missing
PYTHONPATH=src python -m nasa_mouse_generative archs4-catalog
PYTHONPATH=src python -m nasa_mouse_generative split-plan
PYTHONPATH=src python -m nasa_mouse_generative experiment-plan
```

Use the `nasa-mouse` Conda environment when `python` is not on the base shell path.
