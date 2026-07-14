# Visual rendering audit

Audit date: July 14, 2026 UTC.

## Scope

- Rendered and inspected all 33 manuscript PDF pages after integrating the corrected spleen and kidney models, retaining soleus as supplementary sensitivity evidence, and revising the evidence hierarchy.
- Inspected a contact sheet spanning every page, then inspected the quantitative main figures, conceptual synthesis, skin protocol-context figure, and final supplementary pages at higher resolution.
- Inspected the one-page 2026 conference abstract PDF at 160 dpi.
- Checked that all 15 manuscript image references exist and that each HTML image container is balanced.
- Confirmed the final-page Figure S9 caption at 200 dpi after its smaller contact-sheet rendering made the figure number appear truncated.

## Result

No text block, table cell, pathway label, legend, plot, caption, page margin, or page boundary overlaps or clips in the final render. Figure 3 was regenerated after its footer approached the bottom row labels; the revised version is clear at standalone and manuscript scale. Figure S9 contains its complete figure number and caption. The abstract title, author, NASA Ames affiliation, correspondence information, and body remain within the page boundary.

Figure 5 is explicitly conceptual and is not presented as quantitative evidence. Figures 2-4 provide the current quantitative main-tissue, evidence-matrix, and corrected kidney/spleen views. Figures S4-S7 retain the original four-model reviewer checks, Figure S8 retains the generated biological-process artwork as an explicitly illustrative supplement, and Figure S9 records the original-tissue confound and cohort-overlap sensitivity analysis.

| Figure | PNG dimensions | PDF/vector copy | Rendering result |
| --- | ---: | --- | --- |
| Figure 1, workflow | 3129 x 1263 | Yes | Pass |
| Figure 2, retained pathway shifts | 4174 x 3089 | Yes | Pass; landscape in manuscript |
| Figure 3, directional evidence map | 4089 x 2613 | Yes | Pass; landscape in manuscript; footer spacing verified |
| Figure 4, kidney/spleen reassessment | 4235 x 1886 | Yes | Pass; landscape in manuscript |
| Figure 5, complementary process model | 4534 x 2796 | Yes | Pass; landscape in manuscript |
| Figure 6, skin protocol context | 4824 x 2553 | Yes | Pass; landscape in manuscript |
| Figure S1, broad pathway screen | 5194 x 4665 | Yes | Pass; landscape in manuscript |
| Figure S2, skin project balance | 4954 x 2722 | Yes | Pass; landscape in manuscript |
| Figure S3, expanded family review | 4699 x 3317 | Yes | Pass; landscape in manuscript |
| Figure S4, conventional and held-out checks | 5400 x 2700 | Yes | Pass |
| Figure S5, composition-proxy sensitivity | 4800 x 3300 | Yes | Pass |
| Figure S6, full-training-seed sensitivity | 4200 x 3600 | Yes | Pass |
| Figure S7, pathway robustness matrix | 5700 x 2550 | Yes | Pass; footer collision fixed |
| Figure S8, generated biological processes | 4494 x 2635 | Yes | Pass; illustrative status explicit |
| Figure S9, original-tissue sensitivity | 4894 x 2074 | Yes | Pass; caption verified at 200 dpi |

The one-page abstract contains 237 words including its acknowledgment and renders without clipping. All 42 manuscript DOI records resolved through Crossref with exact title matches. The manuscript still contains visible placeholders for the repository DOI, funding, acknowledgments, competing-interest confirmation, and clearance details; these are editorial placeholders, not rendering defects.
