# Metascape Internal API Wrapper

This repo includes a small Metascape web-client wrapper at
`src/nasa_mouse_glare/metascape_client.py`. It automates the same workflow that
we were previously doing manually in the browser:

1. create a Metascape session
2. upload a multi-gene-list CSV
3. select mouse input and analysis species
4. upload/convert a custom enrichment background
5. start enrichment
6. poll progress for GO/PPI jobs
7. generate the report
8. download compact result tables

Important: Metascape does not document these web-app endpoints as a stable
public API. This wrapper mirrors the current browser workflow against
`https://metascape.org/gp_server`. If Metascape changes its frontend contract,
the wrapper may need to be updated.

## Quick Start

Run commands from the repo root, `nasa-mouse`.

For the current 12-filter aggregate liver GLARE run:

```bash
conda run -n nasa env PYTHONPATH=src \
  python -m nasa_mouse_glare.metascape_client submit \
  --gene-lists outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/metascape_gene_lists/metascape_12filter_priority_gene_lists.csv \
  --background outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/metascape_gene_lists/metascape_background_all_glare_genes.txt \
  --output-dir 'outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/metascape_results/{session_id}' \
  --input-species 10090 \
  --analysis-species 10090
```

`10090` is the NCBI taxonomy ID for `Mus musculus`.

The `{session_id}` placeholder is replaced after the wrapper creates a
Metascape session, so each run writes to a separate output directory.

## Validate Without Starting Enrichment

Use `--prepare-only` when you want to test upload, species conversion, custom
background conversion, and ontology/source selection without launching the
enrichment job:

```bash
conda run -n nasa env PYTHONPATH=src \
  python -m nasa_mouse_glare.metascape_client submit \
  --gene-lists outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/metascape_gene_lists/metascape_12filter_priority_gene_lists.csv \
  --background outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/metascape_gene_lists/metascape_background_all_glare_genes.txt \
  --output-dir '/tmp/metascape_client_prepare_test/{session_id}' \
  --input-species 10090 \
  --analysis-species 10090 \
  --prepare-only
```

This writes:

- `upload_response_summary.json`
- `background_summary.json`
- `metascape_run_summary.json`

## Check Or Download Existing Sessions

Print report/session metadata without embedding the large HTML report:

```bash
conda run -n nasa env PYTHONPATH=src \
  python -m nasa_mouse_glare.metascape_client status \
  --session-id <metascape_session_id>
```

Download compact result files for a browser-created or previous CLI-created
session:

```bash
conda run -n nasa env PYTHONPATH=src \
  python -m nasa_mouse_glare.metascape_client download \
  --session-id <metascape_session_id> \
  --output-dir outputs/metascape_runs/<metascape_session_id>
```

Use `--include-zip` if you also want Metascape's full `all.zip` archive.

## Inputs

`--gene-lists`

A CSV where each column is one foreground gene list. This is the same format as
the Metascape "Multiple Gene Lists" browser upload. For GLARE outputs, these
files are produced under:

```text
<run-dir>/post_analysis/metascape_gene_lists/
```

`--background`

A plain text, CSV, TSV, or semicolon-separated list of background genes. The
wrapper sends this through Metascape's background conversion endpoint and then
adds any converted foreground genes that were not already present in the
background. This matches the manual custom-background workflow and avoids
dropping uploaded foreground genes from the analysis universe.

## Outputs

For a full `submit` run, the wrapper writes:

```text
AnalysisReport.html
report_information.json
metascape_run_summary.json
upload_response_summary.json
background_summary.json
Enrichment_GO/_FINAL_GO.csv
Enrichment_GO/GO_AllLists.csv
Enrichment_heatmap/HeatmapSelectedGO.csv
Enrichment_PPI/GO_MCODE.csv
Enrichment_PPI/_FINAL_MCODE.csv
Enrichment_PPI/MCODE.csv
Enrichment_QC/GO_PaGenBase.csv
Enrichment_QC/GO_TRRUST.csv
```

The browser report is available at:

```text
https://metascape.org/gp/index.html#/reportfinal/<session_id>
```

## Progress And Failure Behavior

The wrapper prints timestamped progress logs. During enrichment, it polls
Metascape job status every `--poll-interval` seconds.

Example from a successful run:

```text
Metascape GO status: 0/47 GO
Metascape GO status: 17/47 EVIDENCE
Metascape GO status: 27/47 GPEC_GO
Requesting Metascape report generation
Metascape report: https://metascape.org/gp/index.html#/reportfinal/t0efo3prt
```

Useful runtime flags:

```text
--timeout-minutes 180
--poll-interval 30
--disable-ppi
--no-download
--include-zip
```

If Metascape returns an HTTP error, empty conversion result, no checked
ontology terms, or a job timeout, the wrapper raises `MetascapeError` and exits
non-zero.

## Verified Run

The wrapper was smoke-tested and then used for a full 12-filter aggregate liver
GLARE enrichment run.

Session:

```text
t0efo3prt
```

Report:

```text
https://metascape.org/gp/index.html#/reportfinal/t0efo3prt
```

Local output directory:

```text
outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers/post_analysis/metascape_results/t0efo3prt
```

Run summary:

```text
converted input genes: 2,995
final custom background IDs: 20,979
selected ontology/source terms: 8
runtime: about 2.5 minutes
```

## Endpoint Map

The wrapper currently mirrors these Metascape web-app calls:

```text
GET  /get_session_id
POST /upload_excel
POST /apply_species
POST /convert_background_list_2_gene_id
POST /termMembershipCount
POST /enrichmentanalysismultiplelist
GET  /get_job_status
GET  /make_analysis_report
GET  /get_report_information
GET  /get_file
```

Keep this implementation isolated. If Metascape changes one endpoint, update
`src/nasa_mouse_glare/metascape_client.py` rather than spreading Metascape
request logic through the analysis scripts.

## Manual Fallback

If the wrapper breaks, use the browser workflow:

1. open `https://metascape.org/gp/index.html#/main/step1`
2. upload the multi-list CSV
3. enable "Multiple Gene Lists"
4. choose "Custom Analysis"
5. set input species and analysis species to `M. musculus`
6. paste the background gene list into custom background
7. run enrichment and download the result archive

