# Synthetic-promoted gene literature annotation

## Purpose

The synthetic-guided screen produced 26 promoted tissue-gene associations. A
promoted gene crossed the repeated feature-selection rule only with an eligible
synthetic-informed arm. Promotion does not mean that the gene was absent from
the real data, and it does not establish biological novelty.

This review separates prior-literature agreement from new hypothesis generation.
The search was completed on 2026-08-03.

## Decision rules

Each tissue-gene association receives one mutually exclusive classification:

- `aligning`: prior evidence agrees with the observed direction or the same-tissue
  biological process.
- `contradictory`: the relevant prior result is opposite and there is no comparable
  supporting evidence.
- `complementary`: prior work supports a related spaceflight process or gene
  mechanism without reproducing the same gene, tissue, and direction.
- `ambiguous`: relevant evidence is mixed or depends on mission or assay context.
- `unsupported/potentially_novel`: the targeted search found no relevant direct
  match. This is not proof that the association has never been reported.

The `evidence_scope` column records whether agreement is exact or process-level.
The `evidence_relationship` column states whether the source is independent,
potentially overlaps the OSDR aggregate, or supplies mechanistic context only.
The `interpretive_role` column states whether the result recovers prior evidence,
extends a mechanism or process, is context dependent, or remains unmatched in the
searched tissue literature. Literature status and biological plausibility are
separate questions.

## Results

| Classification | Associations |
|---|---:|
| Aligning | 11 |
| Complementary | 13 |
| Ambiguous | 1 |
| Unsupported/potentially novel | 1 |
| Contradictory | 0 |

Only three associations have direct same-gene, same-tissue, same-direction
support:

- thymus `Ccnb2`, flight lower;
- thymus `Ccne2`, flight lower;
- gastrocnemius `Nfkbia`, flight higher.

The two thymus matches come from published two-mission RNA-seq that may overlap
the current OSDR aggregate, so they are not independent confirmation. The
`Nfkbia` match comes from an independent shuttle mission and microarray/PCR
platform. Eight other flight-lower thymus genes align at the cell-cycle process
level rather than as exact prior gene results.

Thymus `Birc5` is ambiguous. An earlier shuttle study reported higher expression
after flight, whereas later ISS work reported a broad lower cell-cycle program.
Thymus `Hsd17b11` and `Etv1` are complementary. Primary studies connect
HSD17B11 to lipid-droplet lipolysis and ETV1 to CD4+ T-cell activation and
proliferation, but neither supplies a directional thymus spaceflight match.

The 13 complementary associations connect kidney phosphoinositide signaling,
muscle stress and mechanotransduction, skin interferon response, and spleen
adhesion or inflammatory mechanics to earlier findings. These links support
hypothesis development, not replication.

Adrenal `Psmb8` is the only literature-unmatched candidate. PSMB8 is an
interferon-inducible immunoproteasome subunit linked to inflammation and
adipocyte homeostasis, so the association is biologically plausible. The
targeted search found no prior adrenal spaceflight result for the gene.

## Reproducibility

The curated decisions and source metadata live in
`src/nasa_mouse_diffusion/paper_parity/annotate_promoted_gene_literature.py`. The command
below rebuilds the two paper tables and validates complete coverage of the frozen
promoted-gene inventory:

```bash
PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python \
  -m nasa_mouse_diffusion.paper_parity.annotate_promoted_gene_literature
```

Use `--check` to verify that committed tables match the curated definitions.

- `paper/synthetic_guided_spaceflight/source_data/table_s16_promoted_gene_literature_annotations.tsv`
- `paper/synthetic_guided_spaceflight/source_data/table_s17_promoted_gene_literature_sources.tsv`
