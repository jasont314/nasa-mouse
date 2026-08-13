# Figure reproduction

The repository keeps final figures beside each paper and records the tables or
model outputs used to make them. Two kinds of rebuild are available:

- A **frozen-source rebuild** redraws a figure from committed tables. It checks
  plotting code and document assembly without rerunning a model.
- A **full analysis rebuild** recreates the source tables from downloaded data
  and trained models, then redraws the figures.

These are different claims. A frozen-source rebuild should be the first check
after cloning the repository. Use the full path only when auditing the complete
analysis.

## Fresh-clone rebuild

Run these commands from the repository root after installing the environment:

```bash
python -m nasa_mouse_expimap.build_publication_figures --from-frozen-source
python -m nasa_mouse_expimap.render_asgsr_documents
python -m nasa_mouse_expimap.build_asgsr_poster

python \
  -m nasa_mouse_diffusion.paper_parity.build_synthetic_guided_paper \
  --figures-from-frozen-source

python -m nasa_mouse_internship_report.build_report
python -m nasa_mouse_diffusion.paper_parity.build_slstp_presentation
```

The first command redraws the final expiMap architecture, latent-map, pathway,
member-gene, protocol, hypothesis, and selected supplementary figures from
`paper/asgsr_expimap_hvg/source_data/`. The generative command redraws Figures
1, 3, 4, and 5 and Figure S2 from that paper's tracked tables. It preserves the
three model-output graphics described below and then renders the manuscript and
supplement.

The internship report builder redraws all report figures from committed source
tables or committed paper figures. The final presentation builder reads only
tracked source tables, paper figures, templates, and assets under
`presentation/final/source/`.

## Model-output graphics

Three generative-paper graphics display outputs that cannot be reconstructed
from summary tables alone:

| Figure | Tracked copy | Full regeneration |
|---|---|---|
| Figure 2A, ARCHS4 denoising trajectory | `paper/synthetic_guided_spaceflight/figures/figure_2a_archs4_denoising_trajectory.*` | Run the selected ARCHS4 DDIM `prepare`, `train`, and `evaluate` commands in `outputs/COMMANDS.md` |
| Figure 2B, real and synthetic OSDR PCA | `paper/synthetic_guided_spaceflight/figures/figure_2b_locked_real_vs_synthetic_pca.*` | Run the selected adapter training, calibration, and one-time locked-test command |
| Figure S1, muscle-arm heatmap | `paper/synthetic_guided_spaceflight/figures/figure_s1_muscle_arm_heatmap.*` | Run `within_study_feature_stability` with the muscle-group configuration |

Their final PNG and vector PDF copies are committed, so papers and slides can
still be rendered without checkpoints. Exact model commands and paths are in
[`outputs/COMMANDS.md`](../outputs/COMMANDS.md); checkpoint requirements are in
[`ARTIFACTS.md`](../ARTIFACTS.md).

The expiMap supplementary sensitivity figures S2 through S6 and S8 are also
committed. A full regeneration uses `build_asgsr_paper`,
`review_expanded_pathway_screen`, `reviewer_robustness_analysis`,
`run_asgsr_seed_sensitivity`, and `integrate_reassessed_tissues_paper` after the
selected model outputs have been recreated. The complete order is in the
expiMap paper README and the command ledger.

## Poster and presentations

The poster builder uses the tracked NASA template and expiMap publication
figures. The final presentation builder uses its tracked template, source
tables, and assets. Both are reproducible from a clone.

The midpoint presentation is an archival, manually authored PowerPoint rather
than a code-generated deck. Its PPTX is the editable source. A new review PDF
can be exported with:

```bash
libreoffice --headless --convert-to pdf --outdir presentation/midpoint \
  presentation/midpoint/SLSTP_2026_Midpoint_Presentation.pptx
```

## Visual checks

The builders check dimensions and required files. After rebuilding a PDF or
PPTX, render its pages or slides to images and inspect text wrapping, clipping,
legends, and image placement. The recorded manuscript checks are in the paper
packages' visual-audit files. A successful build verifies reproducibility of
the artifact, not the biological conclusion; the source tables and analysis
commands remain the evidence trail.
