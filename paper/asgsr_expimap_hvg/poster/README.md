# ASGSR expiMap poster

`asgsr_expimap_poster.pptx` is an editable, single-slide scientific poster sized at 48 x 27 inches in landscape orientation. It preserves the 16:9 proportions, header treatment, NASA branding, three-column organization, and acknowledgement footer of `assets/poster_template/00 Poster Session Student Template_Approved by Legal.pptx`. The objective, cross-study confounding problem, and pathway-wiring explanation are adapted from `assets/poster_template/Biomedical_Foundation_Models_SLSTP_Midpoint_portable.pptx`; the preliminary midpoint findings are replaced by the final robustness-filtered thymus, skin, liver, and spleen results.

The PDF is the print-ready export, and `asgsr_expimap_poster_preview.png` is a 4,800 x 2,700 pixel visual-check copy. `assets/expimap_architecture_visualization_300dpi.png` is a separate high-resolution export of the compact reference-training and query-mapping panel for reuse in talks or documents.

The title, architecture, tables, evidence summaries, interpretation text, and footer remain native PowerPoint objects. The pathway plot is regenerated from the paper's source data at a wider poster-specific aspect ratio and 400 dpi. It shows only the 13 primary retained programs from manuscript Table 2. The program-score workflow also contains a 400-dpi heatmap of the actual seed-2020 spleen project effects for the three retained pathways; its values come from `outputs/expimap/analyses/kidney_spleen_reassessment/seed_accession_effects.tsv.gz`, with OSD-288 excluded. The tissue-state figure is rendered at 700 dpi from the paper's vector PDF. All embedded data figures retain more than 300 effective pixels per inch at their placed poster dimensions.

The poster abstract is condensed for viewing distance but retains the manuscript abstract's methods, tissue directions, and strength-of-evidence language. The central results figure displays the retained literature-aligned, complementary, and context-sensitive pathways that support the paper's main interpretation across thymus, skin, liver, and spleen. Literature-supported but non-retained pathways, including skin keratinization, are not plotted as findings. Kidney is absent from the primary scope and results displays and remains only an exploratory result in the conclusion.

The annotation workflow identifies OpenAI Codex (GPT-5) as the LLM used to organize pathway effects, mission context, and primary literature. It also shows the required source-verification step before assigning the manuscript's literature-aligned, complementary, or context-sensitive roles; the LLM output is not presented as independent biological evidence.

Rebuild from the repository root:

```bash
python -m nasa_mouse_expimap.build_asgsr_poster
```

The builder requires `python-pptx` and uses `pdftocairo` for high-resolution figure rendering. When LibreOffice is available, it also regenerates the PDF and PNG preview.

Interpretation safeguards are embedded in the poster: lower scores are decoder-oriented latent changes rather than proof of biochemical inhibition, bulk RNA-seq cannot separate cell abundance from cell state, and the displayed tissue-state interpretations require independent cell-resolved and functional validation.
