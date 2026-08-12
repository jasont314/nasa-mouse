# Mentor handoff

Last reviewed: 2026-08-12

This is the completed handoff for Jason Trinh's 2026 NASA Space Life Sciences
Training Program internship. The repository preserves the final reports and
figures in Git while keeping large public reference matrices and trained model
files outside normal Git history.

## Recommended reading order

1. Read the [internship report](paper/slstp_internship_report/manuscript.pdf) for
   the project question, methods, main figures, and biological interpretation.
2. Review the [final presentation](presentation/SLSTP_2026_Generative_Transcriptomics.pdf)
   for the concise visual narrative. The editable deck contains speaker notes.
3. Use the [expiMap manuscript](paper/asgsr_expimap_hvg/manuscript.pdf) for the
   pathway analysis and its study-level checks.
4. Use the [generative manuscript](paper/synthetic_guided_spaceflight/manuscript.pdf)
   for DDIM/WGAN validation, tissue-specific classifier tests, and
   synthetic-informed gene and pathway analyses.
5. Consult [outputs/README.md](outputs/README.md) and
   [outputs/COMMANDS.md](outputs/COMMANDS.md) only when tracing a result back to
   a run or rebuilding an analysis.

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

| Task | Fresh clone | Extra requirements |
|---|---:|---|
| Read PDFs and inspect source tables | Yes | None |
| Edit the PPTX or manuscripts | Yes | Office/PDF tools as needed |
| Rebuild the main report and final presentation | Yes | `nasa-mouse` environment |
| Rebuild publication figures from frozen tables | Mostly | See each paper README |
| Rerun OSDR ingestion | Yes | Network access to NASA OSDR |
| Rerun final model inference | No | Selected local checkpoints |
| Retrain GLARE or run composition checks | No | Tabula Muris Senis H5AD |
| Retrain expiMap or DDIM on ARCHS4 | No | Full ARCHS4 mouse H5 and GPU |

The exact files and checksums are listed in [ARTIFACTS.md](ARTIFACTS.md).

## Repository status

- `main` contains the final report, presentation, manuscripts, source tables,
  code, configs, and selected compact outputs.
- `internship-final-2026-08-12` marks the completed handoff state.
- Historical transfer notes are under `docs/archive/`; they are not current
  operating instructions.
- No license has been selected. The mentor and author should choose one before
  distributing the code for reuse.
- Manuscripts remain research drafts until the author list, acknowledgments,
  repository release URL, and submission requirements receive final review.

## Suggested preservation step

The Git tag preserves the code and compact results. Long-term reproducibility
also requires copying the files listed in `outputs/MODEL_ARTIFACTS.sha256` and
`assets/EXTERNAL_ARTIFACTS.sha256` to managed storage, then adding that storage
location to `ARTIFACTS.md`. At present, the checkpoint manifest identifies the
files but the repository has no public checkpoint download location.

Contact: Jason Trinh, `jasontrinh@berkeley.edu`.
