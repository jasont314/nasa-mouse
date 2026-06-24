# FLT14/GC8 Focused Metascape Summary

Session: `tbvuggxmr`

Report URL:

```text
https://metascape.org/gp/index.html#/reportfinal/tbvuggxmr
```

Inputs:

- `FLT14_all`
- `FLT14_and_GC8`
- `FLT14_not_GC8`
- `GC8_all`
- `GC8_not_FLT14`

Run settings:

- Input species: `M. musculus` (`10090`)
- Analysis species: `M. musculus` (`10090`)
- Converted foreground genes: `1,313`
- Final custom background IDs: `20,979`
- p cutoff: `0.01`
- minimum overlap: `3`
- minimum enrichment: `1.5`

## Main Result

`GC8` is the primary ground-control cluster counterpart for `FLT14`, but the
focused enrichment shows it is not a clean one-to-one equivalent.

- `FLT14_and_GC8` is strongly liver/xenobiotic/CYP enriched.
- `FLT14_not_GC8` is dominated by ribosomal/translation terms rather than the
  CYP/retinol signature.
- `GC8_all` mixes the CYP/retinol signal with a large non-FLT14 component.
- `GC8_not_FLT14` is enriched for DNA metabolism, DNA replication, DNA damage,
  and cilium/basal-body terms.

This supports treating `GC8` as the nearest GC neighborhood for `FLT14`, not as
the direct biological equivalent of the whole FLT14 module.

## Top Terms By List

```tsv
GeneList	Top terms
FLT14_all	xenobiotic metabolic process; cytoplasmic ribosomal proteins; cellular response to xenobiotic stimulus; steroid metabolic process; ribosome; monocarboxylic acid metabolic process
FLT14_and_GC8	xenobiotic metabolic process; cellular response to xenobiotic stimulus; retinol metabolism; xenobiotic catabolic process; steroid hormone biosynthesis; chemical carcinogenesis/DNA adducts
FLT14_not_GC8	cytoplasmic ribosomal proteins; negative regulation of myoblast fusion; ribosome rescue/NMD/SRP-dependent cotranslational targeting; cytoplasmic translation
GC8_all	retinol metabolism; cytochrome P450 substrates; oxidative demethylation; xenobiotic metabolic process; basal body/cilium terms; DNA metabolic process
GC8_not_FLT14	DNA metabolic process; DNA-templated DNA replication; primary cilium assembly; cilium assembly; basal-body anchoring; DNA damage response
```

## Files

- `Enrichment_GO/GO_AllLists.csv`
- `Enrichment_GO/_FINAL_GO.csv`
- `metascape_top_terms_by_list.tsv`
- `AnalysisReport.html`
- `metascape_run_summary.json`

