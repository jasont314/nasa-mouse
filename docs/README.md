# Documentation

This directory contains the small set of technical references that remain useful
after the internship handoff. Scientific results, figures, and detailed evidence
tables live with the papers rather than being duplicated here.

## Maintained references

| File | Purpose |
|---|---|
| [Codebase guide](codebase_guide.md) | Package dependencies, data flow, status boundaries, and result tracing |
| [NASA OSDR API](osdr_api.md) | Dataset discovery, metadata fields, endpoints, and local API use |
| [Method sources](method_sources.md) | Upstream repositories, pinned commits, and method provenance |
| [Generative pipeline](generative_pipeline.md) | Configurable DDIM/WGAN workflow and selected implementation branch |
| [Generative validation](generative_validation.md) | Final OSDR-disjoint model comparison and validation summary |
| [Figure reproduction](figure_reproduction.md) | Fresh-clone figure builds, model-dependent exceptions, and visual checks |
| [Biological hypotheses](distributed_response_hypotheses.md) | Literature-supported interpretation used by the internship report |
| [Literature annotations](annotation_provenance.md) | LLM-assisted annotation tables, rationale fields, source catalogs, and rebuild scripts |
| [Annotation prompts](annotation_prompts.md) | Original instruction sequence, final label definitions, canonical rerun prompts, and historical limits |

## Scientific record

Use the repository's paper packages for conclusions and publication-facing
methods:

- [`paper/slstp_internship_report/`](../paper/slstp_internship_report/) combines
  the GLARE, expiMap, and generative work.
- [`paper/asgsr_expimap_hvg/`](../paper/asgsr_expimap_hvg/) contains the detailed
  Reactome pathway analysis.
- [`paper/synthetic_guided_spaceflight/`](../paper/synthetic_guided_spaceflight/)
  contains generator validation and synthetic-informed gene and pathway results.

Run selection and exact commands are documented in
[`outputs/README.md`](../outputs/README.md) and
[`outputs/COMMANDS.md`](../outputs/COMMANDS.md). Model configuration files are
under [`configs/generative/`](../configs/generative/).

Superseded agent handoffs, presentation work logs, intermediate result summaries,
and rejected-model narratives were removed from the working tree. They remain
available in Git history and are not current operating instructions.
