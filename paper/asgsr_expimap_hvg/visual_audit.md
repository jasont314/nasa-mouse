# Visual rendering audit

Audit date: July 14, 2026 UTC.

## Scope

- Rebuilt the manuscript HTML and PDF after replacing the legacy main figures with a five-figure quantitative set and adding Figure S10.
- Rendered and inspected all 34 manuscript pages at 120 dpi as a contact sheet and then inspected pages 20-34 individually at full rendered resolution.
- Inspected the one-page 2026 conference abstract at 160 dpi.
- Checked that all 15 manuscript image references resolve, every final figure has a vector PDF copy, all PDF fonts are embedded, and the manuscript PDF contains no suspect objects or encryption.
- Re-ran the DOI audit; all 42 manuscript DOI records resolved through Crossref with exact title matches.

## Result

No text block, table cell, pathway label, legend, plot, caption, page margin, or page boundary overlaps or clips in the final render. Main figures are authored at a 7.2-inch final width and exported as 300-dpi PNG plus vector PDF. Dense legacy supplementary panels use a Letter-landscape page with a 7.55-inch plot column and a separate caption column; labels remain legible at rendered size.

The main set now exposes the workflow and model scope, tissue-specific reference-query support, project and complete-training effects, orthogonal directional checks, retained-pathway member-gene support, and skin protocol context. Figure S10 exposes individual sample distributions without treating samples as independent mission replicates. Figure S8 is a deterministic vector synthesis; no generated biological artwork remains in the final package.

Color is not the sole encoding of a reported distinction. Flight and ground use triangle and circle markers, literature roles use distinct marker shapes, and directional-support cells contain plus or minus symbols. Axes, legends, panel labels, and captions remain readable at manuscript scale.

| Figure | PNG dimensions | PDF copy | Rendering result |
| --- | ---: | --- | --- |
| Figure 1, workflow and architecture | 2160 x 1395 | Yes | Pass |
| Figure 2, latent mapping | 2160 x 1755 | Yes | Pass |
| Figure 3, project and training shifts | 2160 x 1875 | Yes | Pass |
| Figure 4, directional and gene support | 2160 x 2025 | Yes | Pass |
| Figure 5, skin protocol context | 2160 x 1905 | Yes | Pass |
| Figure S1, broad pathway screen | 2160 x 2640 | Yes | Pass; long terms shortened only in the figure |
| Figure S2, skin project balance | 4954 x 2722 | Yes | Pass; landscape |
| Figure S3, expanded family review | 4699 x 3317 | Yes | Pass; landscape |
| Figure S4, conventional and held-out checks | 5400 x 2700 | Yes | Pass; landscape |
| Figure S5, composition-proxy sensitivity | 4800 x 3300 | Yes | Pass; landscape |
| Figure S6, complete-training sensitivity | 4200 x 3600 | Yes | Pass; landscape |
| Figure S7, pathway robustness matrix | 2160 x 2235 | Yes | Pass; landscape |
| Figure S8, vector process synthesis | 2160 x 1395 | Yes | Pass |
| Figure S9, confound and overlap sensitivity | 4894 x 2074 | Yes | Pass; landscape |
| Figure S10, sample-score distributions | 2160 x 1575 | Yes | Pass |

The abstract remains one page and renders without clipping. The manuscript still contains visible placeholders for the repository DOI, funding, acknowledgments, contribution review, and final clearance confirmation. These are submission and editorial tasks, not figure-rendering defects.
