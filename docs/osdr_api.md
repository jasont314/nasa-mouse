# NASA OSDR API Notes

These notes cover the NASA Open Science Data Repository Biological Data API
used for discovering OSDR metadata, assay/sample structure, and public data
files relevant to this project.

Sources checked on 2026-06-26:

- [Official API landing page](https://visualization.osdr.nasa.gov/biodata/api/)
- [NASA OSDR API announcement](https://www.nasa.gov/osdr-api/)
- [OSDR tutorial guide](https://osdr-tutorials.readthedocs.io/en/latest/pages/guides/biological_data_api.html)
- [OSDR repository](https://osdr.nasa.gov/bio/repo)

## What the API Provides

The Biological Data API gives programmatic access to GeneLab and ALSDA data in
NASA OSDR. It has two useful access patterns:

- REST traversal: structured JSON-style traversal of datasets, assays,
  samples, and metadata.
- Query endpoints: table-oriented metadata and data queries with filters,
  intended for analysis and file discovery.

The query interface is usually the better fit for this repo because we often
need to filter by organism, tissue, assay type, spaceflight condition, or file
datatype.

## Base URL

```text
https://visualization.osdr.nasa.gov/biodata/api/v2/
```

Useful endpoints:

```text
/query/metadata/
/query/data/
/query/assays/
/dataset/{OSD_ACCESSION}/
/dataset/{OSD_ACCESSION}/assay/*/sample/*/
```

## Common Fields

OSDR metadata follows the ISA model. The API exposes ISA fields as
dot-delimited query keys:

```text
id.accession
id.assay name
id.sample name
study.characteristics.organism
study.characteristics.material type
study.factor value.spaceflight
file.datatype
file.filename
```

Spaces in field names or values must be URL-encoded as `%20`.

## Example Queries

All Mus musculus metadata rows with unnormalized count files:

```text
https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?study.characteristics.organism=mus%20musculus&file.datatype=unnormalized%20counts
```

Mouse left-kidney spaceflight metadata with unnormalized count files:

```text
https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?study.characteristics.organism=mus%20musculus&study.characteristics.material%20type=left%20kidney&study.factor%20value.spaceflight&file.datatype=unnormalized%20counts
```

Matching data query for the same filter:

```text
https://visualization.osdr.nasa.gov/biodata/api/v2/query/data/?study.characteristics.organism=mus%20musculus&study.characteristics.material%20type=left%20kidney&study.factor%20value.spaceflight&file.datatype=unnormalized%20counts
```

Inspect all assays and samples for one dataset:

```text
https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-48/assay/*/sample/*/
```

Request a tabular view explicitly:

```text
https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?id.accession=OSD-48&study.characteristics&format=csv
```

## Project Use

Use the OSDR API for provenance and discovery:

- find public mouse OSDR studies and assay names;
- inspect ISA sample metadata such as tissue, organism, and flight condition;
- locate file records and data types before downloading large artifacts;
- reproduce metadata joins used by local preprocessing.

For expiMap OSDR inputs in this repo, use the API path rather than the older
local integrated OSDR HDF5:

```bash
PYTHONPATH=src python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics
PYTHONPATH=src python -m nasa_mouse_glare.prepare_expimap_osdr_tissue --tissue liver
PYTHONPATH=src python -m nasa_mouse_glare.prepare_expimap_osdr_tissue --tissue kidney
```

Keep downloaded API count CSVs under `data/osdr_api/counts/`; that directory is
ignored by git. Store small provenance files, manifests, and scripts in git so
the OSDR inputs can be rediscovered or regenerated.
