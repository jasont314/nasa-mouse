# SLSTP presentation and project handoff

Last updated: 2026-08-04

## Resume prompt

In a new Codex conversation, use:

> Read `/media/volume/mouse/nasa/nasa-mouse/docs/handoffs/presentation.md`, inspect the referenced files, and continue from the current repository state. Do not reset or stage unrelated worktree changes. Push after each new commit.

This file transfers the working context and decisions. It does not contain the complete chat transcript, and a new conversation will not inherit hidden conversational state.

## Git state

- Branch: `main`
- Current presentation commit: use `git log -1 --oneline`; the latest deck and handoff are pushed to `origin/main`.
- Inspect `git status` before editing and do not reset unrelated worktree changes.
- User preference: create checkpoint commits for meaningful progress and push after every commit.

## Current deliverables

- Presentation: `presentation/SLSTP_2026_Generative_Transcriptomics.pptx`
- PDF: `presentation/SLSTP_2026_Generative_Transcriptomics.pdf`
- Speaker notes: `presentation/generative_slstp_2026/speaker_notes.md`
- Rebuild script: `src/nasa_mouse_diffusion/paper_parity/build_slstp_presentation.py`
- NASA template: `presentation/SLSTP_template_2026.pptx`
- Midpoint reference: `presentation/Biomedical_Foundation_Models_SLSTP_Midpoint_portable (2).pptx`

The current deck has 26 slides, 26 embedded note pages, and a planned runtime of 12:45. All slides were rendered and visually checked after the most recent rebuild.

## Working on another device

A fresh clone contains everything needed to open and continue editing the current presentation. It also contains the inputs used by the presentation builder:

- the NASA PowerPoint template;
- the generated PPTX and PDF;
- the presentation builder and speaker notes;
- the expiMap manuscript figures and retained-pathway source tables;
- the generative manuscript figures and source tables;
- the DDIM architecture excerpt used in the validation slide; and
- `environment.yml` and the repository requirements file.

The following local files are not part of the Git repository:

- `presentation/Biomedical_Foundation_Models_SLSTP_Midpoint_portable (2).pptx` is an ignored local visual reference. The builder does not require it.
- Raw H5/H5AD expression files and large model artifacts are ignored by
  `.gitignore`; curated output summaries and publication inputs are tracked.
- Multi-gigabyte DDIM and WGAN checkpoints under `outputs/` are local only.

Therefore, another device can edit the deck and rebuild its current figures and slides from the checked-in summary data. It cannot retrain models or rerun analyses that require raw ARCHS4/OSDR matrices or saved checkpoints until those large files are transferred through separate storage. Do not add multi-gigabyte checkpoints to ordinary Git history.

## Presentation flow

1. Title: interpretable and generative models for mouse spaceflight.
2. One OSDR dataset, two machine-learning questions.
3. Neural-network and autoencoder introduction.
4. expiMap architecture using tissue-matched ARCHS4 references and mouse Reactome programs.
5. expiMap pathway findings in thymus, skin, liver, and spleen.
6. Synthetic transcriptomics introduction.
7. OSDR and ARCHS4 data scope.
8. Configurable preprocessing, harmonization, model, training, and conditioning pipeline.
9. WGAN-GP versus DDIM validation metrics and DDIM architecture.
10. DDIM reverse-diffusion PCA trajectory.
11. Locked PCA colored by tissue beside the same PCA colored by OSDR accession.
12. Locked PCA colored by FLT/GC condition beside the same PCA colored by OSDR accession.
13. Five tissue-specific analysis arms.
14. Real, generated, and consensus gene ranking with top-k classifier inputs.
15. Schematic real-only versus synthetic-guided FLT/GC classifier boundaries.
16. Pooled and tissue-specific classifier utility, including positive and negative examples.
17. Reinforced, promoted, and real-only feature-selection sets.
18. Independent literature-interpretation categories.
19. Coverage across 27 completed tissue analyses.
20. Complete 49-association synthetic-informed gene inventory.
21. Thymus biological interpretation.
22. Soleus biological interpretation.
23. Additional findings in pooled muscle, kidney, spleen, and skin.
24. Additional findings in eye, adrenal gland, gastrocnemius, and tibialis anterior.
25. Takeaways.
26. Acknowledgments and questions.

## Presentation decisions

- The presentation covers the entire project, not only generative modeling.
- Roughly the first quarter introduces neural networks, autoencoders, and expiMap. The remaining slides cover the generative pipeline and analysis.
- GLARE is omitted from the short project review.
- The title slide follows the supplied NASA template.
- Slide 6 keeps only the enlarged three-stage synthetic-transcriptomics diagram. The explanation that generated profiles are not new biological evidence is in the speaker notes.
- Slides 11 and 12 use the same locked coordinates for all 293 observed and 293 matched DDIM profiles. Slide 11 pairs tissue and accession coloring; slide 12 pairs FLT/GC condition and accession coloring. Circles are observed profiles and crosses are DDIM profiles. The model was conditioned on accession, so the accession panels check conditional fidelity rather than batch removal.
- Slide 15 contains only the two classifier diagrams and their legend. The removed comparison panel and detailed leakage-control sentence remain explained verbally.
- The same transparent held-out profiles appear in both slide 15 panels. They are scored, not fitted.
- Slide 16 shows both gains and declines. Cecum, colon, and liver are the displayed non-improving synthetic-informed candidates; non-worse performance across BA, AUROC, and AP remains the retention rule.
- Slide 25 has only Generate, Use, and Interpret. It has no next-test section and does not mention pooled augmentation.
- Wording was reviewed using the local humanizer guidance. Avoid formulaic phrasing, promotional claims, em dashes, and unnecessary technical caveats on the slides.

## expiMap summary used in the deck

The expiMap manuscript is under `paper/asgsr_expimap_hvg/`.

- Reference scope: 17,708 tissue-matched ARCHS4 profiles.
- Primary effect scope: 700 OSDR profiles.
- Models retained about 2,000 highly variable genes and 319 to 387 Reactome programs per primary tissue.
- Thymus: lower DNA repair, RHOA cytoskeletal-cycle, and lymphoid-stromal interaction programs.
- Skin: lower chromatin regulation, DNA repair, Hedgehog signaling, sphingolipid metabolism, and cell-junction organization.
- Liver: lower MHC class II antigen presentation and T-cell receptor signaling, with heterogeneous metabolic responses.
- Spleen: lower T-cell receptor signaling, neutrophil degranulation, and C-type lectin receptor signaling. All three passed conventional GSEA FDR below 0.05 and all five directional checks.
- The supporting evidence tables remain checked in, but the dedicated expiMap evidence-matrix slide was removed from the presentation.
- expiMap implementation code is under `src/nasa_mouse_expimap/`.

## Generative summary used in the deck

- The benchmark compares WGAN-GP and conditional DDIM.
- The selected branch uses TPM, training-fitted MaxAbs scaling, 974 mouse landmark genes, ARCHS4 pretraining, and OSDR adaptation.
- The OSDR model is conditioned on tissue, FLT/GC, accession, and material type.
- DDIM was selected because it retained high expression fidelity while giving higher F1, adversarial accuracy near chance, and lower distributional distance than WGAN-GP in the reported benchmark.
- Synthetic profiles are used through generated-only, real-plus-generated, consensus-guided, and lightly weighted synthetic-training arms.
- Final FLT/GC effects and Benjamini-Hochberg FDR are calculated from observed OSDR profiles.
- The feature inventory contains 23 reinforced and 26 promoted tissue-gene associations, for 49 synthetic-informed associations in total.
- Literature status is independent of selection status: aligning, complementary, ambiguous, or unmatched.

## Build and verification

Rebuild the deck:

```bash
/home/exouser/miniforge3/condabin/conda run -n nasa-mouse \
  python src/nasa_mouse_diffusion/paper_parity/build_slstp_presentation.py
```

Regenerate the PDF:

```bash
libreoffice --headless --convert-to pdf --outdir presentation \
  presentation/SLSTP_2026_Generative_Transcriptomics.pptx
```

Minimum checks before committing:

```bash
/home/exouser/miniforge3/condabin/conda run -n nasa-mouse \
  python -m py_compile src/nasa_mouse_diffusion/paper_parity/build_slstp_presentation.py
git diff --check
```

Render every PDF page with `pdftoppm`, inspect a contact sheet, and inspect changed slides at full resolution. The most recent verified state had 26 PPTX slides, 26 embedded note pages, and 26 PDF pages.

Only stage these files for presentation-only changes unless the task explicitly requires more:

```text
presentation/SLSTP_2026_Generative_Transcriptomics.pdf
presentation/SLSTP_2026_Generative_Transcriptomics.pptx
presentation/generative_slstp_2026/assets/SOURCES.md
presentation/generative_slstp_2026/speaker_notes.md
presentation/generative_slstp_2026/source_data/*.tsv
src/nasa_mouse_diffusion/paper_parity/build_slstp_presentation.py
docs/handoffs/presentation.md
```

## Recent presentation commits

```text
9a70a76 Add presentation transfer handoff
d4006f1 Expand presentation to full project scope
2b5f257 Clarify held-out evaluation visualization
a825023 Show five-gene consensus selection example
a5ab1d5 Explain top-k classifier feature selection
a885332 Clarify generated draws in presentation
691dc31 Split consensus ranking and classifier visuals
0268237 Clarify consensus ranking visualization
8d64e63 Explain consensus ranking in presentation
221b9d3 Visualize synthetic guidance mechanism
6b59457 Fix slide 2 title spacing
e42052a Clean up slide 4 outcome label
68a7088 Add synthetic transcriptomics introduction slides
19fa571 Polish PCA caption and utility labels
68a1a7e Refine feature selection slide spacing
01a8c4c Clarify synthetic analysis workflow slides
e082915 Enlarge diffusion architecture figure
ab076be Rebalance diffusion validation slide
e542028 Use paper DDIM architecture and template title
efe17c5 Redesign presentation for scientific clarity
e9e357d Clarify selected pipeline configuration
```

Use `git log --oneline --decorate -20` to refresh this list after new presentation commits.
