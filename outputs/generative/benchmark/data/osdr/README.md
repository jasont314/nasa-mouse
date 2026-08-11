# API-Derived OSDR Matrix

`osdr_api_raw_counts.h5ad` is generated from the 75 NASA OSDR API count tables;
it is not the older integrated OSDR H5. Technical replicates are summed by default.
The matrix is ignored by Git; `osdr_api_expression_summary.json` records its shape,
count semantics, provenance, and the `OSD-759` REST-download fallback.
