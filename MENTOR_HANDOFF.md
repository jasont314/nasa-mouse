# Mentor handoff

Last reviewed: 2026-08-13

This is the completed handoff for Jason Trinh's 2026 NASA Space Life Sciences
Training Program internship. The repository preserves the final reports and
figures in Git while keeping large public reference matrices and trained model
files outside normal Git history.

## Recommended reading order

1. Read the [internship report](paper/slstp_internship_report/manuscript.pdf) for
   the project question, methods, main figures, and biological interpretation.
2. Use the [results guide](docs/results_guide.md) to move directly from each
   headline finding to its exact source table, plot, and interpretation.
   The [selected gene and pathway bundle](outputs/comparison/selected_features/README.md)
   lists the comparison results by tissue. Its
   [workbook](outputs/comparison/selected_features/selected_feature_comparison.xlsx)
   contains all 16 retained expiMap pathways, their 743 unique tissue-gene
   members, 1,307 stable tissue-gene pairs across all generative classifier
   arms, the narrower 679-row selected-arm set, 21 primary matched genes, 49
   secondary consensus genes, ten grouped Reactome results, and an
   analysis-coverage sheet.
3. Review the [final presentation](presentation/final/SLSTP_2026_Generative_Transcriptomics.pdf)
   for the concise visual narrative. The editable deck contains speaker notes.
4. Use the [expiMap manuscript](paper/asgsr_expimap_hvg/manuscript.pdf) for the
   pathway analysis and its study-level checks.
5. Use the [generative manuscript](paper/synthetic_guided_spaceflight/manuscript.pdf)
   for DDIM/WGAN validation, tissue-specific classifier tests, and
   synthetic-informed gene and pathway analyses.
6. Read the [codebase guide](docs/codebase_guide.md) before navigating scripts,
   configurations, or retained development experiments.
7. Use the [literature annotation record](docs/annotation_provenance.md) and
   [annotation prompts](docs/annotation_prompts.md) to trace expiMap pathways
   and synthetic feature-importance genes to their rationale, cited sources,
   and review instructions.
8. Consult [outputs/README.md](outputs/README.md) and
   [outputs/COMMANDS.md](outputs/COMMANDS.md) only when tracing a result back to
   a run or rebuilding an analysis.
9. Use the [figure reproduction guide](docs/figure_reproduction.md) for the
   exact fresh-clone and full-analysis figure paths.
10. See the [editorial review record](docs/editorial_review.md) for the files
   reviewed with Humanizer and the content that the prose pass did not alter.

## What was completed

### OSDR data layer

The repository discovers *Mus musculus* bulk transcriptomic RNA-seq profiles
through the NASA OSDR Biological Data API. It canonicalizes tissue and material
labels while retaining accession, condition, assay, sex when available, and
muscle-group information. The final generative inventory contains 1,610 FLT and
ground-control profiles from 75 accessions.

### GLARE

GLARE was the first modeling approach. It learned expression representations
after Tabula Muris Senis pretraining, with MOBER tested as a batch-correction
step. Study effects weakened in several views, but FLT and GC remained mixed.
This stage is reported as an informative negative result rather than a final
biological discovery workflow.

### expiMap

Tissue-matched ARCHS4 references were trained with current mouse Reactome masks,
and OSDR profiles were mapped as tissue-specific queries. The final paper focuses
on thymus, skin, liver, and spleen, with kidney retained as exploratory evidence.
The analysis includes accession-level effects, repeated training, conventional
gene-set scoring, member-gene review, and composition sensitivity checks.

### Generative models

The configurable benchmark tested expression transformations, harmonization,
feature sets, cohort scope, study and tissue conditioning, and training source.
WGAN-GP and DDIM implementations followed their paper architectures. The final
branch pretrained DDIM on ARCHS4, adapted it to OSDR with tissue, condition,
study, and material conditioning, and evaluated generated profiles against real
validation samples. DDIM was used for downstream analysis because it had the
stronger set of reported fidelity and real-versus-synthetic metrics.

Synthetic data were then tested as training support for tissue-specific FLT/GC
classifiers and as a source of gene and pathway ranking. Thymus produced the
strongest combined interpretation. Liver, skin, and spleen had narrower matched
gene findings, while soleus appeared in the secondary compact-panel analysis.
All biological association tests were calculated from real OSDR profiles.

## Reproducibility levels

| Task | Recoverable from a clone | Extra requirements |
|---|---:|---|
| Read PDFs and inspect source tables | Yes | None |
| Edit the PPTX or manuscripts | Yes | Office/PDF tools as needed |
| Redraw frozen-source figures; render manuscripts, poster, and final presentation | Yes | `nasa-mouse` environment; LibreOffice for PPTX-to-PDF export |
| Rebuild the historical generative tables exactly | No, not from Git alone | Restore the 3.9 MB compact input snapshot named in `frozen_input_manifest.tsv` |
| Recompute the generative analyses | Yes | Restore checkpoints or retrain, then run `outputs/COMMANDS.md` |
| Rerun OSDR ingestion | Yes | Network access to NASA OSDR |
| Rerun final model inference without training | No | Selected local checkpoints |
| Retrain GLARE or run composition checks | Yes | Public TMS download and OSDR API data |
| Retrain expiMap or DDIM on ARCHS4 | Yes | Public ARCHS4 download, OSDR API data, and GPU |

The exact files and checksums are listed in [ARTIFACTS.md](ARTIFACTS.md).

## Repository status

- `main` contains the final report, presentation, manuscripts, source tables,
  code, configs, and selected compact outputs.
- `internship-final-2026-08-13-r6` marks the completed handoff state. It
  supersedes the earlier handoff tags after licensing, repository cleanup, and
  final verification.
- Superseded development handoffs and intermediate analyses remain available in
  Git history; they are not current operating instructions.
- A clean checkout contains about 108 MiB of tracked files, but the full Git
  object pack is about 884 MB because older commits contain deleted datasets
  and intermediate outputs. Use a shallow clone for review-only access. History
  was not rewritten because doing so would replace every shared commit
  identifier.
- Project-authored content is released under the MIT License. Vendored code,
  reference files, and public datasets retain their upstream terms as recorded
  in `THIRD_PARTY_NOTICES.md`.
- Manuscripts remain research drafts until the author list, acknowledgments,
  repository release URL, and submission requirements receive final review.

## Suggested preservation step

The Git tag preserves the code and compact results. ARCHS4, TMS, and GENCODE are
public and can be restored with
`python -m nasa_mouse_generative prepare-references`; they do not require a
private handoff copy. Long-term preservation should focus on the files listed
in `outputs/MODEL_ARTIFACTS.sha256`, which avoid roughly 2.3 GB of checkpoint
regeneration, and the compact ignored inputs named by the generative paper's
`frozen_input_manifest.tsv`. The repository has no public checkpoint download
location yet. `ARTIFACTS.md` includes an archive and restore recipe for these
local files.

Contact: Jason Trinh, `jasontrinh@berkeley.edu`.
