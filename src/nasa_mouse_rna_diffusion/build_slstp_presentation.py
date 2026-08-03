"""Build the 2026 SLSTP presentation for the generative transcriptomics study."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE, PP_PLACEHOLDER
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "presentation/SLSTP_template_2026.pptx"
OUTPUT = ROOT / "presentation/SLSTP_2026_Generative_Transcriptomics.pptx"
PACKAGE_DIR = ROOT / "presentation/generative_slstp_2026"
ASSET_DIR = PACKAGE_DIR / "assets"
NOTES_PATH = PACKAGE_DIR / "speaker_notes.md"
PAPER_DIR = ROOT / "paper/synthetic_guided_spaceflight"

SLIDE_W = 13.333333
SLIDE_H = 7.5
FONT = "Arial"

NAVY = "082B55"
BLUE = "2D6496"
TEAL = "178681"
GREEN = "478A64"
ORANGE = "D37B00"
GOLD = "9A6B20"
CORAL = "D96552"
PURPLE = "8063A6"
DARK = "263746"
GRAY = "59697A"
MID_GRAY = "8A969E"
LIGHT = "F2F5F7"
PALE_BLUE = "E9F1F7"
PALE_TEAL = "E8F3F1"
PALE_GOLD = "F8F1E2"
PALE_CORAL = "F8ECE9"
WHITE = "FFFFFF"


@dataclass(frozen=True)
class SlideNote:
    number: int
    title: str
    time: str
    text: str


def _rgb(value: str) -> RGBColor:
    return RGBColor.from_string(value)


def _set_fill(shape, color: str, transparency: int | None = None) -> None:
    shape.fill.solid()
    shape.fill.fore_color.rgb = _rgb(color)
    if transparency is not None:
        shape.fill.transparency = transparency


def _set_line(shape, color: str, width: float = 1.0) -> None:
    shape.line.color.rgb = _rgb(color)
    shape.line.width = Pt(width)


def _add_text(
    slide,
    text: str,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    size: float = 18,
    color: str = DARK,
    bold: bool = False,
    italic: bool = False,
    align: PP_ALIGN = PP_ALIGN.LEFT,
    valign: MSO_ANCHOR = MSO_ANCHOR.TOP,
    margin: float = 0.02,
    name: str | None = None,
):
    shape = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        shape.name = name
    frame = shape.text_frame
    frame.clear()
    frame.word_wrap = True
    frame.margin_left = Inches(margin)
    frame.margin_right = Inches(margin)
    frame.margin_top = Inches(margin)
    frame.margin_bottom = Inches(margin)
    frame.vertical_anchor = valign
    paragraph = frame.paragraphs[0]
    paragraph.alignment = align
    paragraph.space_after = Pt(0)
    run = paragraph.add_run()
    run.text = text
    run.font.name = FONT
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = _rgb(color)
    return shape


def _add_rich_text(
    slide,
    segments: Iterable[tuple[str, dict]],
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    size: float = 18,
    color: str = DARK,
    valign: MSO_ANCHOR = MSO_ANCHOR.TOP,
    margin: float = 0.02,
):
    shape = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    frame = shape.text_frame
    frame.clear()
    frame.word_wrap = True
    frame.margin_left = Inches(margin)
    frame.margin_right = Inches(margin)
    frame.margin_top = Inches(margin)
    frame.margin_bottom = Inches(margin)
    frame.vertical_anchor = valign
    paragraph = frame.paragraphs[0]
    paragraph.space_after = Pt(0)
    for text, style in segments:
        run = paragraph.add_run()
        run.text = text
        run.font.name = FONT
        run.font.size = Pt(style.get("size", size))
        run.font.bold = style.get("bold", False)
        run.font.italic = style.get("italic", False)
        run.font.color.rgb = _rgb(style.get("color", color))
    return shape


def _add_panel(
    slide,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    fill: str = LIGHT,
    line: str = "DDE4E8",
    radius: bool = True,
):
    shape_type = MSO_SHAPE.ROUNDED_RECTANGLE if radius else MSO_SHAPE.RECTANGLE
    shape = slide.shapes.add_shape(
        shape_type, Inches(x), Inches(y), Inches(w), Inches(h)
    )
    _set_fill(shape, fill)
    _set_line(shape, line, 0.8)
    return shape


def _add_rule(slide, x: float, y: float, w: float, color: str, height: float = 0.04):
    shape = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(height)
    )
    _set_fill(shape, color)
    shape.line.fill.background()
    return shape


def _add_circle(slide, x: float, y: float, diameter: float, color: str):
    shape = slide.shapes.add_shape(
        MSO_SHAPE.OVAL,
        Inches(x),
        Inches(y),
        Inches(diameter),
        Inches(diameter),
    )
    _set_fill(shape, color)
    shape.line.fill.background()
    return shape


def _add_arrow(slide, x: float, y: float, w: float, h: float, color: str = MID_GRAY):
    shape = slide.shapes.add_shape(
        MSO_SHAPE.RIGHT_ARROW, Inches(x), Inches(y), Inches(w), Inches(h)
    )
    _set_fill(shape, color)
    shape.line.fill.background()
    return shape


def _add_bullet_rows(
    slide,
    rows: list[str],
    x: float,
    y: float,
    w: float,
    *,
    size: float = 16,
    color: str = DARK,
    bullet_color: str = TEAL,
    row_h: float = 0.48,
):
    for index, row in enumerate(rows):
        top = y + index * row_h
        _add_circle(slide, x, top + 0.12, 0.10, bullet_color)
        _add_text(
            slide,
            row,
            x + 0.18,
            top,
            w - 0.18,
            row_h,
            size=size,
            color=color,
            valign=MSO_ANCHOR.MIDDLE,
        )


def _add_picture_contain(
    slide,
    path: Path,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    alt: str,
):
    with Image.open(path) as image:
        aspect = image.width / image.height
    box_aspect = w / h
    if aspect >= box_aspect:
        width = w
        height = w / aspect
        left = x
        top = y + (h - height) / 2
    else:
        height = h
        width = h * aspect
        left = x + (w - width) / 2
        top = y
    shape = slide.shapes.add_picture(
        str(path), Inches(left), Inches(top), Inches(width), Inches(height)
    )
    shape.name = alt
    shape._element.nvPicPr.cNvPr.set("descr", alt)
    return shape


def _add_slide_title(slide, eyebrow: str, title: str, subtitle: str | None = None):
    _add_text(
        slide,
        eyebrow.upper(),
        0.32,
        0.58,
        5.5,
        0.25,
        size=10.5,
        color=GOLD,
        bold=True,
    )
    _add_text(
        slide,
        title,
        0.32,
        0.86,
        12.5,
        0.64,
        size=27,
        color=NAVY,
        bold=True,
        valign=MSO_ANCHOR.MIDDLE,
    )
    if subtitle:
        _add_text(
            slide,
            subtitle,
            0.34,
            1.48,
            12.3,
            0.34,
            size=14.5,
            color=GRAY,
            valign=MSO_ANCHOR.MIDDLE,
        )


def _prepare_content_slide(slide, number: int):
    for placeholder in slide.placeholders:
        kind = placeholder.placeholder_format.type
        if kind in (PP_PLACEHOLDER.TITLE, PP_PLACEHOLDER.OBJECT):
            placeholder.text = ""
        elif kind == PP_PLACEHOLDER.FOOTER:
            placeholder.text = "Biological & Physical Sciences"
        elif kind == PP_PLACEHOLDER.SLIDE_NUMBER:
            placeholder.text = str(number)
    return slide


def _add_source(slide, text: str):
    _add_text(
        slide,
        text,
        0.36,
        7.12,
        12.0,
        0.19,
        size=7.2,
        color=MID_GRAY,
        valign=MSO_ANCHOR.MIDDLE,
    )


def _build_tissue_utility_chart() -> Path:
    table = pd.read_csv(PAPER_DIR / "source_data/table_s9_all_tissue_development_screen.tsv", sep="\t")
    order = [
        "thymus",
        "skeletal_muscle",
        "soleus",
        "kidney",
        "spleen",
        "skin",
        "lung",
        "adrenal_gland",
    ]
    muscle = pd.read_csv(PAPER_DIR / "source_data/table_s6_muscle_group_summary.tsv", sep="\t")
    rows = []
    for tissue in order:
        source = muscle if tissue == "soleus" else table
        row = source.loc[source["tissue"].eq(tissue)].iloc[0]
        rows.append(row)
    data = pd.DataFrame(rows).reset_index(drop=True)
    display = {
        "thymus": "Thymus",
        "skeletal_muscle": "Skeletal muscle",
        "soleus": "Soleus",
        "kidney": "Kidney",
        "spleen": "Spleen",
        "skin": "Skin",
        "lung": "Lung",
        "adrenal_gland": "Adrenal gland",
    }
    real = data["real_mean_balanced_accuracy"].to_numpy(float)
    selected = data["selected_mean_balanced_accuracy"].to_numpy(float)
    y = np.arange(len(data))
    fig, ax = plt.subplots(figsize=(8.1, 4.6))
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")
    for yi, start, end in zip(y, real, selected):
        ax.plot([start, end], [yi, yi], color="#AEB9BF", lw=2.0, zorder=1)
    ax.scatter(real, y, s=62, color="#7E8A92", label="Real only", zorder=3)
    ax.scatter(selected, y, s=72, color="#178681", label="Selected synthetic use", zorder=4)
    for yi, start, end in zip(y, real, selected):
        ax.text(end + 0.009, yi, f"+{end-start:.3f}", va="center", fontsize=9, color="#263746")
    ax.set_yticks(y, [display[value] for value in data["tissue"]])
    ax.invert_yaxis()
    ax.set_xlim(0.45, 1.01)
    ax.set_xlabel("Balanced accuracy on held-out real profiles", fontsize=10)
    ax.grid(axis="x", color="#DDE4E8", lw=0.8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=10)
    ax.tick_params(axis="x", labelsize=9)
    ax.legend(
        frameon=False,
        ncol=2,
        fontsize=9,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
    )
    fig.tight_layout()
    path = ASSET_DIR / "tissue_balanced_accuracy_dumbbell.png"
    fig.savefig(path, dpi=220, transparent=True, bbox_inches="tight")
    plt.close(fig)
    return path


def _build_trajectory_crop() -> Path:
    source = PAPER_DIR / "figures/figure_2a_archs4_denoising_trajectory.png"
    with Image.open(source) as image:
        top = int(image.height * 0.07)
        bottom = int(image.height * 0.84)
        cropped = image.crop((0, top, image.width, bottom))
        path = ASSET_DIR / "ddim_trajectory_panels.png"
        cropped.save(path)
    return path


def _set_title_slide(slide):
    main = next(
        placeholder
        for placeholder in slide.placeholders
        if placeholder.placeholder_format.type == PP_PLACEHOLDER.BODY
    )
    main.text_frame.clear()
    paragraph = main.text_frame.paragraphs[0]
    paragraph.space_after = Pt(0)
    run = paragraph.add_run()
    run.text = "Synthetic transcriptomics\nfor mouse spaceflight"
    run.font.name = FONT
    run.font.size = Pt(30)
    run.font.bold = True
    run.font.color.rgb = _rgb(WHITE)

    subtitle = next(
        placeholder
        for placeholder in slide.placeholders
        if placeholder.placeholder_format.type == PP_PLACEHOLDER.SUBTITLE
    )
    subtitle.text_frame.clear()
    subtitle.text_frame.word_wrap = True
    entries = [
        ("Can generated RNA-seq improve tissue-specific FLT vs GC analysis?", 17, False, "D8E1EB"),
        ("", 5, False, WHITE),
        ("Jason Trinh", 18, True, WHITE),
        ("EECS | UC Berkeley", 12, False, "B9C7D5"),
        ("Mentor: James Casaletto | August 2026", 12, False, "B9C7D5"),
    ]
    for index, (text, size, bold, color) in enumerate(entries):
        paragraph = subtitle.text_frame.paragraphs[0] if index == 0 else subtitle.text_frame.add_paragraph()
        paragraph.space_after = Pt(2)
        run = paragraph.add_run()
        run.text = text
        run.font.name = FONT
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.color.rgb = _rgb(color)
    _add_text(slide, "SLSTP FINAL PRESENTATION", 0.54, 1.08, 3.8, 0.25, size=10.5, color="D5A64A", bold=True)


def _slide_2(slide):
    _add_slide_title(
        slide,
        "Question",
        "Small studies and study effects complicate tissue comparisons",
        "Can a generator help rank spaceflight signal without pretending it creates new animals?",
    )
    _add_panel(slide, 0.35, 2.02, 6.05, 2.35, fill=PALE_BLUE, line="C9DCE9")
    _add_rule(slide, 0.35, 2.02, 6.05, BLUE)
    _add_text(slide, "NASA OSDR", 0.65, 2.31, 2.2, 0.32, size=18, color=BLUE, bold=True)
    _add_text(slide, "1,610", 0.62, 2.68, 2.25, 0.72, size=36, color=NAVY, bold=True)
    _add_text(slide, "biological profiles", 0.66, 3.32, 2.3, 0.28, size=13, color=GRAY)
    _add_text(slide, "75", 3.08, 2.70, 1.3, 0.55, size=30, color=NAVY, bold=True)
    _add_text(slide, "accessions", 3.10, 3.26, 1.5, 0.27, size=13, color=GRAY)
    _add_text(slide, "835 FLT", 4.62, 2.73, 1.35, 0.38, size=20, color=CORAL, bold=True)
    _add_text(slide, "775 GC", 4.62, 3.22, 1.35, 0.38, size=20, color=BLUE, bold=True)

    _add_panel(slide, 6.72, 2.02, 6.25, 2.35, fill=PALE_TEAL, line="C7DDD8")
    _add_rule(slide, 6.72, 2.02, 6.25, TEAL)
    _add_text(slide, "ARCHS4 mouse", 7.02, 2.31, 2.5, 0.32, size=18, color=TEAL, bold=True)
    _add_text(slide, "997,515", 7.00, 2.68, 2.45, 0.72, size=34, color=NAVY, bold=True)
    _add_text(slide, "profiles audited", 7.03, 3.32, 2.2, 0.28, size=13, color=GRAY)
    _add_text(slide, "17,244", 9.50, 2.70, 1.95, 0.55, size=29, color=NAVY, bold=True)
    _add_text(slide, "selected", 9.53, 3.26, 1.25, 0.27, size=13, color=GRAY)
    _add_text(slide, "20 tissues", 11.25, 2.75, 1.45, 0.38, size=19, color=TEAL, bold=True)

    _add_panel(slide, 0.35, 4.67, 12.62, 1.75, fill=WHITE, line="DDE4E8")
    _add_bullet_rows(
        slide,
        [
            "Mission, strain, material and assay differences can resemble a flight effect.",
            "The generator must preserve tissue and FLT/GC structure without memorizing samples.",
            "Statistical evidence must still come from observed OSDR profiles.",
        ],
        0.75,
        4.94,
        11.8,
        size=16,
        bullet_color=ORANGE,
        row_h=0.45,
    )
    _add_source(slide, "Sources: NASA OSDR Biological Data API; ARCHS4 mouse v2.5.")


def _slide_3(slide):
    _add_slide_title(
        slide,
        "Method",
        "We compared three model families and used conditional diffusion",
        "The model architecture followed each paper; preprocessing, harmonization and conditioning remained configurable.",
    )
    stages = [
        (0.40, 2.20, 2.65, "Expression", "TPM\n974 mouse landmarks\ntrain-fitted MaxAbs", BLUE, PALE_BLUE),
        (3.34, 2.20, 2.65, "Models", "WGAN-GP\nDDIM\nGeneJEPA screen", ORANGE, PALE_GOLD),
        (6.28, 2.20, 2.65, "Training", "ARCHS4 pretrain\nOSDR adaptation\npooled tissue model", TEAL, PALE_TEAL),
        (9.22, 2.20, 3.70, "Conditions", "tissue + FLT/GC\naccession + material\nmatched profile generation", GREEN, "ECF4ED"),
    ]
    for x, y, w, title, body, color, fill in stages:
        _add_panel(slide, x, y, w, 2.25, fill=fill, line=color)
        _add_rule(slide, x, y, w, color)
        _add_text(slide, title, x + 0.22, y + 0.28, w - 0.44, 0.33, size=17, color=color, bold=True)
        _add_text(slide, body, x + 0.22, y + 0.79, w - 0.44, 1.18, size=16, color=DARK)
    for x in [3.07, 6.01, 8.95]:
        _add_arrow(slide, x, 3.05, 0.22, 0.32)

    _add_panel(slide, 0.40, 4.78, 12.52, 1.45, fill="F7F9FA", line="DDE4E8")
    _add_text(slide, "Selected branch", 0.70, 5.04, 1.7, 0.28, size=13, color=GOLD, bold=True)
    _add_rich_text(
        slide,
        [
            ("No global batch correction. ", {"bold": True, "color": NAVY}),
            ("Accession was supplied as a model condition so study structure was represented instead of erased.", {"color": DARK}),
        ],
        2.25,
        4.99,
        9.95,
        0.55,
        size=15.5,
    )
    _add_text(
        slide,
        "Nine liver harmonization arms and a 463-row experiment plan were screened before full training.",
        2.25,
        5.58,
        9.95,
        0.33,
        size=12.5,
        color=GRAY,
    )
    _add_source(slide, "Model sources: Vinas et al., Bioinformatics (2022); Lacan et al., BMC Bioinformatics (2026); GeneJEPA (2025).")


def _slide_4(slide, trajectory: Path):
    _add_slide_title(
        slide,
        "Diffusion",
        "Diffusion learns tissue structure from noise",
        "The same PCA axes follow 1,024 generated profiles through reverse diffusion.",
    )
    _add_picture_contain(
        slide,
        trajectory,
        0.34,
        1.78,
        12.64,
        4.72,
        alt="DDIM reverse trajectory at timesteps 1000, 200 and 0",
    )
    _add_text(slide, "t = 1000: noise", 0.78, 6.50, 2.5, 0.28, size=13, color=GRAY, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "t = 200: partial structure", 5.23, 6.50, 2.9, 0.28, size=13, color=GRAY, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "t = 0: tissue-conditioned profiles", 9.15, 6.50, 3.4, 0.28, size=13, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
    _add_source(slide, "Gray points are real ARCHS4 profiles; colors identify generated tissue conditions. PC1 and PC2 limits are shared across panels.")


def _slide_5(slide):
    _add_slide_title(
        slide,
        "Validation",
        "DDIM matched expression and was harder to distinguish from real data",
        "AA closer to 0.5 means a classifier has less success separating real and generated profiles; lower FD is better.",
    )
    _add_panel(slide, 0.38, 2.02, 4.05, 3.72, fill=PALE_BLUE, line="C9DCE9")
    _add_text(slide, "ARCHS4 tissue probe", 0.67, 2.30, 3.45, 0.35, size=17, color=BLUE, bold=True)
    _add_text(slide, "Balanced accuracy on held-out real profiles", 0.68, 2.72, 3.35, 0.30, size=12.5, color=GRAY)
    base_y = 5.18
    chart_top = 3.17
    chart_h = 1.82
    for index, (label, value, color) in enumerate([
        ("Real-trained", 0.781, MID_GRAY),
        ("Synthetic-trained", 0.781, TEAL),
    ]):
        x = 1.04 + index * 1.72
        bar_h = chart_h * value
        shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(base_y - bar_h), Inches(0.86), Inches(bar_h))
        _set_fill(shape, color)
        shape.line.fill.background()
        _add_text(slide, f"{value:.3f}", x - 0.03, base_y - bar_h - 0.36, 0.95, 0.30, size=16, color=NAVY, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, label, x - 0.30, base_y + 0.08, 1.45, 0.40, size=12.5, color=DARK, align=PP_ALIGN.CENTER)
    _add_rule(slide, 0.82, base_y, 3.15, "A9B5BC", 0.015)

    _add_panel(slide, 4.72, 2.02, 8.23, 3.72, fill=WHITE, line="DDE4E8")
    _add_text(slide, "Metric", 5.03, 2.27, 1.55, 0.28, size=12.5, color=GRAY, bold=True)
    _add_text(slide, "WGAN-GP", 7.15, 2.27, 1.35, 0.28, size=12.5, color=CORAL, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "DDIM", 9.03, 2.27, 1.35, 0.28, size=12.5, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "Reading", 10.72, 2.27, 1.78, 0.28, size=12.5, color=GRAY, bold=True)
    metrics = [
        ("Correlation", "0.976", "0.974", "similar"),
        ("F1", "0.985", "0.997", "higher"),
        ("Adversarial accuracy", "0.636", "0.475", "closer to 0.5"),
        ("FD / real P95", "0.144", "0.074", "lower"),
    ]
    for index, (label, wgan, ddim, reading) in enumerate(metrics):
        y = 2.78 + index * 0.63
        if index % 2 == 0:
            shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(4.94), Inches(y - 0.07), Inches(7.72), Inches(0.53))
            _set_fill(shape, "F5F7F8")
            shape.line.fill.background()
        _add_text(slide, label, 5.03, y, 1.95, 0.34, size=14, color=DARK, bold=index >= 2)
        _add_text(slide, wgan, 7.13, y, 1.38, 0.34, size=16, color=CORAL, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, ddim, 9.01, y, 1.38, 0.34, size=16, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, reading, 10.74, y, 1.62, 0.34, size=13.5, color=GRAY)

    _add_panel(slide, 0.38, 6.04, 12.57, 0.73, fill=NAVY, line=NAVY)
    _add_text(slide, "Use DDIM downstream", 0.72, 6.19, 2.55, 0.34, size=18, color=WHITE, bold=True, valign=MSO_ANCHOR.MIDDLE)
    _add_text(slide, "Lower separability and distributional distance, with comparable correlation and higher F1.", 3.12, 6.20, 9.18, 0.32, size=15, color="DCE7F2", valign=MSO_ANCHOR.MIDDLE)
    _add_source(slide, "WGAN values use validation data and DDIM values use the stated OSDR test; they are not paired on one common split.")


def _slide_6(slide):
    _add_slide_title(
        slide,
        "Analysis",
        "Synthetic profiles entered the analysis in five different ways",
        "Every arm was judged on held-out real profiles before synthetic attribution was allowed.",
    )
    _add_panel(slide, 0.38, 2.05, 2.42, 1.28, fill=PALE_BLUE, line=BLUE)
    _add_text(slide, "Real OSDR", 0.64, 2.30, 1.88, 0.30, size=17, color=BLUE, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "observed FLT / GC", 0.64, 2.74, 1.88, 0.25, size=12.5, color=GRAY, align=PP_ALIGN.CENTER)
    _add_panel(slide, 0.38, 3.63, 2.42, 1.28, fill=PALE_CORAL, line=CORAL)
    _add_text(slide, "DDIM generated", 0.57, 3.88, 2.05, 0.30, size=17, color=CORAL, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "matched conditions", 0.64, 4.31, 1.88, 0.25, size=12.5, color=GRAY, align=PP_ALIGN.CENTER)
    _add_arrow(slide, 2.93, 3.12, 0.52, 0.38)

    cards = [
        (3.62, 2.10, "Real only", "rank: real\nfit: real", BLUE, PALE_BLUE),
        (5.47, 2.10, "Generated only", "rank: generated\nfit: generated", CORAL, PALE_CORAL),
        (7.32, 2.10, "Real + generated", "rank: consensus\nfit: equal weight", TEAL, PALE_TEAL),
        (4.55, 3.73, "Guided, real fit", "rank: consensus\nfit: real", PURPLE, "F0ECF6"),
        (6.40, 3.73, "Guided, 5% synth.", "rank: consensus\nfit: real + 5%", ORANGE, PALE_GOLD),
    ]
    for x, y, title, body, color, fill in cards:
        _add_panel(slide, x, y, 1.65, 1.28, fill=fill, line=color)
        _add_text(slide, title, x + 0.10, y + 0.18, 1.45, 0.35, size=13.5, color=color, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, body, x + 0.13, y + 0.66, 1.39, 0.50, size=11.5, color=DARK, align=PP_ALIGN.CENTER)

    _add_arrow(slide, 9.32, 3.12, 0.52, 0.38)
    _add_panel(slide, 10.00, 2.38, 2.88, 2.25, fill="ECF4ED", line=GREEN)
    _add_text(slide, "Held-out real profiles", 10.23, 2.69, 2.42, 0.50, size=17, color=GREEN, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "BA | AUROC | AP", 10.28, 3.35, 2.30, 0.36, size=16, color=NAVY, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "select one eligible arm\nfor each tissue", 10.28, 3.83, 2.30, 0.50, size=12.5, color=GRAY, align=PP_ALIGN.CENTER)

    _add_panel(slide, 0.38, 5.38, 12.50, 1.15, fill=NAVY, line=NAVY)
    _add_text(slide, "Biology check", 0.70, 5.65, 1.75, 0.34, size=17, color=WHITE, bold=True)
    _add_text(slide, "FLT vs GC effects, random-effects meta-analysis and BH FDR use real OSDR profiles only.", 2.52, 5.66, 7.40, 0.32, size=15.5, color="DCE7F2")
    _add_text(slide, "Synthetic profiles never increase animal n.", 10.05, 5.54, 2.46, 0.56, size=14, color="FFD69A", bold=True, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE)


def _slide_7(slide, utility_chart: Path):
    _add_slide_title(
        slide,
        "Utility",
        "Pooling tissues hid useful signal",
        "A single FLT/GC classifier did not improve, but tissue-specific synthetic use often did.",
    )
    _add_panel(slide, 0.38, 2.02, 3.50, 4.62, fill=PALE_GOLD, line="E8D6AF")
    _add_text(slide, "Pooled across tissues", 0.67, 2.30, 2.90, 0.35, size=17, color=GOLD, bold=True)
    labels = ["Real only", "Generated only", "Real + synth."]
    values = [0.754, 0.695, 0.737]
    colors = [BLUE, CORAL, TEAL]
    baseline = 6.00
    top = 3.12
    scale = 3.45
    for index, (label, value, color) in enumerate(zip(labels, values, colors)):
        x = 0.78 + index * 0.96
        bar_h = max(0.05, (value - 0.64) * scale)
        shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(baseline - bar_h), Inches(0.62), Inches(bar_h))
        _set_fill(shape, color)
        shape.line.fill.background()
        _add_text(slide, f"{value:.3f}", x - 0.08, baseline - bar_h - 0.35, 0.78, 0.28, size=14, color=NAVY, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, label, x - 0.20, baseline + 0.08, 1.05, 0.55, size=10.5, color=DARK, align=PP_ALIGN.CENTER)
    _add_rule(slide, 0.63, baseline, 2.95, "A9B5BC", 0.015)
    _add_text(slide, "Balanced accuracy", 0.72, 2.78, 2.78, 0.25, size=12.5, color=GRAY)

    _add_panel(slide, 4.12, 2.02, 8.84, 4.62, fill=WHITE, line="DDE4E8")
    _add_text(slide, "Selected use within each tissue", 4.42, 2.29, 4.2, 0.34, size=17, color=TEAL, bold=True)
    _add_picture_contain(
        slide,
        utility_chart,
        4.32,
        2.65,
        8.32,
        3.70,
        alt="Balanced accuracy before and after tissue-specific synthetic use",
    )
    _add_text(slide, "The useful arm differed by tissue.", 4.44, 6.25, 4.5, 0.27, size=12.5, color=GRAY, italic=True)
    _add_source(slide, "Development results with held-out profiles from represented accessions; they do not measure transfer to a new mission.")


def _slide_8(slide):
    _add_slide_title(
        slide,
        "Feature analysis",
        "Synthetic guidance changed ranking, not statistical evidence",
        "Promoted and reinforced describe repeated feature selection; both require a supporting effect in real OSDR data.",
    )
    _add_text(slide, "Repeated stable selection", 0.55, 2.02, 4.0, 0.34, size=16, color=GRAY, bold=True, align=PP_ALIGN.CENTER)
    left = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.62), Inches(2.52), Inches(3.00), Inches(2.35))
    left.fill.background()
    _set_line(left, BLUE, 2.4)
    right = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(2.48), Inches(2.52), Inches(3.00), Inches(2.35))
    right.fill.background()
    _set_line(right, TEAL, 2.4)
    _add_text(slide, "Real-only\nranking", 0.82, 3.03, 1.35, 0.72, size=17, color=BLUE, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "Synthetic-informed\nranking", 3.90, 3.03, 1.35, 0.72, size=16, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
    _add_panel(slide, 2.18, 3.15, 1.72, 0.92, fill=PALE_GOLD, line=ORANGE)
    _add_text(slide, "23\nreinforced", 2.38, 3.28, 1.32, 0.60, size=17, color=ORANGE, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "26 promoted", 3.72, 4.38, 1.46, 0.30, size=15, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "selected only with\nsynthetic guidance", 3.62, 4.74, 1.66, 0.50, size=11.5, color=GRAY, align=PP_ALIGN.CENTER)

    _add_arrow(slide, 5.78, 3.30, 0.65, 0.48, TEAL)
    _add_panel(slide, 6.68, 2.23, 5.98, 3.67, fill=PALE_BLUE, line=BLUE)
    _add_text(slide, "Real OSDR association test", 7.02, 2.52, 5.30, 0.40, size=20, color=BLUE, bold=True, align=PP_ALIGN.CENTER)
    _add_bullet_rows(
        slide,
        [
            "Estimate FLT vs GC inside each accession",
            "Combine accession effects with a random-effects model",
            "Apply BH FDR within each tissue",
        ],
        7.22,
        3.18,
        4.95,
        size=15,
        bullet_color=BLUE,
        row_h=0.58,
    )
    _add_panel(slide, 7.16, 5.12, 5.00, 0.52, fill=NAVY, line=NAVY)
    _add_text(slide, "49 synthetic-informed associations", 7.34, 5.21, 4.65, 0.29, size=16, color=WHITE, bold=True, align=PP_ALIGN.CENTER)

    _add_panel(slide, 0.52, 6.18, 12.15, 0.58, fill="F7F9FA", line="DDE4E8")
    _add_text(slide, "These are 49 tissue-gene associations among 459 BH-FDR rows. A gene can appear in more than one tissue.", 0.72, 6.31, 11.75, 0.30, size=13.5, color=GRAY, align=PP_ALIGN.CENTER)


def _slide_9(slide):
    _add_slide_title(
        slide,
        "Thymus",
        "Thymus points to lower proliferative renewal",
        "Synthetic guidance organized real-data associations into a coherent mitotic and DNA-replication panel.",
    )
    figure = PAPER_DIR / "figures/figure_3_thymus_biology.png"
    _add_picture_contain(slide, figure, 0.35, 1.93, 8.55, 4.78, alt="Thymus gene effects and Reactome processes")
    _add_panel(slide, 9.10, 2.02, 3.82, 4.57, fill=PALE_CORAL, line="E8C9C2")
    _add_text(slide, "16", 9.43, 2.34, 1.10, 0.65, size=35, color=CORAL, bold=True)
    _add_text(slide, "synthetic-informed\nBH-FDR associations", 10.38, 2.38, 2.12, 0.62, size=14, color=DARK, bold=True)
    _add_text(slide, "13 promoted | 3 reinforced", 9.44, 3.13, 3.05, 0.30, size=14, color=TEAL, bold=True)
    _add_bullet_rows(
        slide,
        [
            "Lower Nusap1, Stmn1, Birc5, Cdk1, Top2a and related genes",
            "Cell cycle, DNA replication, APC/C and G2/M processes",
            "Bulk RNA-seq cannot separate cell loss from lower transcription",
        ],
        9.44,
        3.62,
        3.03,
        size=13.2,
        bullet_color=CORAL,
        row_h=0.77,
    )
    _add_source(slide, "Prior thymus studies report involution and altered cell-cycle expression after spaceflight (Gridley et al. 2013; Horie et al. 2019).")


def _slide_10(slide):
    _add_slide_title(
        slide,
        "Soleus",
        "Soleus reinforces a mitochondrial and lipid program",
        "This result strengthens an existing real-data pattern rather than promoting a new soleus gene.",
    )
    figure = PAPER_DIR / "figures/figure_4_soleus_biology.png"
    _add_picture_contain(slide, figure, 0.34, 1.96, 9.15, 4.70, alt="Soleus gene effects and Reactome processes")
    _add_panel(slide, 9.63, 2.04, 3.25, 4.52, fill=PALE_TEAL, line="C7DDD8")
    _add_text(slide, "0.925 -> 0.963", 9.89, 2.34, 2.72, 0.42, size=20, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "balanced accuracy", 9.98, 2.78, 2.55, 0.28, size=12.5, color=GRAY, align=PP_ALIGN.CENTER)
    _add_text(slide, "5 reinforced | 0 promoted", 9.87, 3.27, 2.77, 0.34, size=15, color=NAVY, bold=True, align=PP_ALIGN.CENTER)
    _add_bullet_rows(
        slide,
        [
            "Lower Bdh1, Ech1, Bnip3 and Decr1",
            "Higher Tpm1",
            "Mitochondrial turnover and fatty-acid metabolism",
            "Compatible with unloading and contractile remodeling",
        ],
        9.94,
        3.78,
        2.54,
        size=12.8,
        bullet_color=TEAL,
        row_h=0.58,
    )
    _add_source(slide, "Literature context: Gambara et al. (2017) and Stein et al. (2002).")


def _slide_11(slide):
    _add_slide_title(
        slide,
        "Other tissues",
        "Kidney adds a focused signal; other tissues remain exploratory",
        "The strongest process-level stories were thymus and soleus, but several tissues produced narrower candidates.",
    )
    rows = [
        ("Kidney", "Inpp4b promoted, Slc37a4 reinforced; both higher in flight", "Focused secondary", BLUE, PALE_BLUE),
        ("Spleen", "Rai14, Ptprk and Myl9 promoted; Loxl1 reinforced; no Reactome enrichment", "Exploratory", PURPLE, "F0ECF6"),
        ("Skin + adrenal", "Plscr1 higher in skin; Psmb8 and Tspan4 lower in adrenal", "Exploratory", ORANGE, PALE_GOLD),
        ("Lung + retina", "Predictive gains without a synthetic-informed BH-FDR gene", "Prediction only", TEAL, PALE_TEAL),
        ("Liver + EDL + quadriceps", "Real-only arm retained", "No synthetic claim", MID_GRAY, "F2F4F5"),
    ]
    for index, (tissue, finding, label, color, fill) in enumerate(rows):
        y = 2.02 + index * 0.89
        _add_panel(slide, 0.43, y, 12.44, 0.73, fill=fill, line="DDE4E8", radius=False)
        _add_rule(slide, 0.43, y, 0.08, color, 0.73)
        _add_text(slide, tissue, 0.70, y + 0.17, 2.12, 0.34, size=16, color=color, bold=True)
        _add_text(slide, finding, 2.79, y + 0.14, 7.70, 0.40, size=14, color=DARK, valign=MSO_ANCHOR.MIDDLE)
        _add_text(slide, label, 10.67, y + 0.16, 1.84, 0.34, size=12.5, color=color, bold=True, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE)
    _add_text(slide, "Independent missions are still needed before any candidate becomes a transferable spaceflight biomarker.", 0.75, 6.60, 11.80, 0.33, size=14, color=GRAY, italic=True, align=PP_ALIGN.CENTER)


def _slide_12(slide):
    _add_slide_title(
        slide,
        "Takeaways",
        "Synthetic data helped most as a tissue-specific prior",
        "The generator created useful structure, not new biological replication.",
    )
    columns = [
        (0.42, "1", "Generate", "Conditional DDIM produced high-fidelity profiles with near-chance real-versus-synthetic separation.", BLUE, PALE_BLUE),
        (4.44, "2", "Use", "Pooled augmentation failed. Tissue-specific ranking and low-weight training were more useful.", TEAL, PALE_TEAL),
        (8.46, "3", "Interpret", "Thymus and soleus formed the clearest programs; kidney supplied a focused secondary pair.", ORANGE, PALE_GOLD),
    ]
    for x, number, heading, body, color, fill in columns:
        _add_panel(slide, x, 2.10, 3.74, 2.72, fill=fill, line=color)
        _add_circle(slide, x + 0.26, 2.37, 0.52, color)
        _add_text(slide, number, x + 0.26, 2.43, 0.52, 0.30, size=16, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, heading, x + 0.91, 2.38, 2.45, 0.36, size=19, color=color, bold=True)
        _add_text(slide, body, x + 0.28, 3.08, 3.15, 1.33, size=15, color=DARK)
    _add_panel(slide, 0.42, 5.15, 12.00, 1.30, fill=NAVY, line=NAVY)
    _add_text(slide, "Next test", 0.75, 5.48, 1.38, 0.34, size=17, color="FFD69A", bold=True)
    _add_text(slide, "Hold out complete accessions, then validate the thymus panel, soleus metabolism and kidney Inpp4b / Slc37a4 in independent data.", 2.05, 5.40, 9.92, 0.58, size=16, color=WHITE, valign=MSO_ANCHOR.MIDDLE)
    _add_text(slide, "Generated profiles never enter the biological sample count or the BH-FDR test.", 2.05, 6.05, 9.92, 0.28, size=12.5, color="BFD0E1")


def _slide_13(slide):
    body = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(0.38), Inches(SLIDE_W), Inches(SLIDE_H - 0.38))
    _set_fill(body, NAVY)
    body.line.fill.background()
    _add_text(slide, "Thank you", 0.65, 1.22, 5.40, 0.80, size=34, color=WHITE, bold=True)
    _add_text(slide, "Questions?", 0.67, 2.10, 4.0, 0.45, size=21, color="BFD0E1")
    _add_rule(slide, 0.68, 2.84, 3.0, ORANGE, 0.05)
    _add_text(slide, "Jason Trinh", 0.68, 3.21, 3.0, 0.35, size=18, color=WHITE, bold=True)
    _add_text(slide, "jasontrinh@berkeley.edu", 0.68, 3.63, 3.5, 0.30, size=14, color="BFD0E1")
    _add_text(slide, "Mentor", 6.02, 1.48, 1.30, 0.27, size=11, color="7FA0BC", bold=True)
    _add_text(slide, "James Casaletto", 7.43, 1.45, 4.55, 0.34, size=17, color=WHITE, bold=True)
    _add_text(slide, "Program", 6.02, 2.16, 1.30, 0.27, size=11, color="7FA0BC", bold=True)
    _add_text(slide, "NASA Space Life Sciences Training Program", 7.43, 2.13, 4.82, 0.52, size=15, color=WHITE)
    _add_text(slide, "Data", 6.02, 3.04, 1.30, 0.27, size=11, color="7FA0BC", bold=True)
    _add_text(slide, "NASA OSDR | ARCHS4 | Reactome", 7.43, 3.01, 4.55, 0.34, size=15, color=WHITE)
    _add_text(slide, "Compute", 6.02, 3.73, 1.30, 0.27, size=11, color="7FA0BC", bold=True)
    _add_text(slide, "NASA Ames Research Center", 7.43, 3.70, 4.55, 0.34, size=15, color=WHITE)
    _add_text(slide, "Code and manuscript", 6.02, 4.42, 1.30, 0.40, size=11, color="7FA0BC", bold=True)
    _add_text(slide, "github.com/jasont314/nasa-mouse", 7.43, 4.40, 4.55, 0.34, size=15, color=WHITE)
    _add_text(slide, "All biological associations were tested in observed OSDR profiles.", 6.02, 5.50, 6.00, 0.40, size=13.5, color="BFD0E1", italic=True)


def _add_notes(slide, note: SlideNote) -> None:
    frame = slide.notes_slide.notes_text_frame
    frame.text = f"Target time: {note.time}\n\n{note.text}"


def _write_notes(notes: list[SlideNote]) -> None:
    lines = [
        "# SLSTP 2026 generative transcriptomics speaker notes",
        "",
        "Target length: 12-15 minutes. Planned speaking time: about 13 minutes 30 seconds.",
        "",
    ]
    for note in notes:
        lines.extend([
            f"## {note.number}. {note.title} ({note.time})",
            "",
            note.text,
            "",
        ])
    NOTES_PATH.write_text("\n".join(lines), encoding="utf-8")


def build() -> Path:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    utility_chart = _build_tissue_utility_chart()
    trajectory = _build_trajectory_crop()
    presentation = Presentation(TEMPLATE)
    _set_title_slide(presentation.slides[0])
    _prepare_content_slide(presentation.slides[1], 2)
    while len(presentation.slides) < 13:
        number = len(presentation.slides) + 1
        slide = presentation.slides.add_slide(presentation.slide_layouts[3])
        _prepare_content_slide(slide, number)

    builders = [
        None,
        _slide_2,
        _slide_3,
        lambda slide: _slide_4(slide, trajectory),
        _slide_5,
        _slide_6,
        lambda slide: _slide_7(slide, utility_chart),
        _slide_8,
        _slide_9,
        _slide_10,
        _slide_11,
        _slide_12,
        _slide_13,
    ]
    for index, builder in enumerate(builders):
        if builder is not None:
            builder(presentation.slides[index])

    notes = [
        SlideNote(1, "Synthetic transcriptomics for mouse spaceflight", "0:25", "Introduce the question. This talk asks whether generated expression can help analyze tissue-specific FLT versus GC biology without counting synthetic profiles as additional animals."),
        SlideNote(2, "Small studies and study effects complicate tissue comparisons", "1:00", "OSDR gives broad tissue coverage, but the data are spread across 75 accessions with different mission and assay contexts. ARCHS4 supplies a much larger mouse reference. The challenge is to use that reference without letting study structure masquerade as spaceflight biology."),
        SlideNote(3, "We compared three model families and used conditional diffusion", "1:10", "The framework kept the published model architectures and varied the surrounding choices. GeneJEPA was useful for representation but had no expression decoder. WGAN-GP and DDIM could generate profiles. The chosen DDIM used TPM, 974 mouse landmarks, ARCHS4 pretraining and OSDR adaptation conditioned on tissue, FLT/GC, accession and material."),
        SlideNote(4, "Diffusion learns tissue structure from noise", "0:55", "Read the panels from left to right. The same generated profiles begin as noise, develop structure by timestep 200 and approach tissue-conditioned regions at timestep zero. The axes are shared, so the visual change is not caused by rescaling each panel."),
        SlideNote(5, "DDIM matched expression and was harder to distinguish from real data", "1:15", "Both WGAN-GP and DDIM had high correlation and F1. DDIM had adversarial accuracy near 0.5 and a lower Frechet-distance ratio, so it was harder to separate from real profiles and closer in distribution. The metrics use each model's stated evaluation split, so this is a model-choice summary rather than a paired significance test."),
        SlideNote(6, "Synthetic profiles entered the analysis in five different ways", "1:15", "Synthetic data can be used for direct training, mixed training or feature guidance. Each tissue could choose among five arms. The eligibility check used held-out real profiles. Once features were nominated, FLT/GC effects and BH FDR were computed from observed OSDR samples only."),
        SlideNote(7, "Pooling tissues hid useful signal", "1:10", "The simplest pooled augmentation test was negative: balanced accuracy fell from 0.754 to 0.737 with real plus synthetic training. Tissue-specific analysis changed the result. Different tissues benefited from different synthetic uses, which argues against one global augmentation policy."),
        SlideNote(8, "Synthetic guidance changed ranking, not statistical evidence", "1:05", "Reinforced genes were selected with and without synthetic guidance. Promoted genes crossed the stable-selection rule only with synthetic guidance. Promoted does not mean biologically novel. All 49 synthetic-informed tissue-gene associations also had a supporting effect and BH FDR below 0.05 in real data."),
        SlideNote(9, "Thymus points to lower proliferative renewal", "1:25", "Thymus produced the clearest promoted panel. The lower mitotic and DNA-replication genes agree with prior reports of thymic involution and altered cell-cycle expression after flight. Because this is bulk RNA-seq, the signal could reflect fewer cycling thymocytes, lower transcription within those cells or both."),
        SlideNote(10, "Soleus reinforces a mitochondrial and lipid program", "1:15", "Soleus improved with real plus generated training. The selected genes were already stable in real-only analysis, so synthetic data reinforced rather than introduced the panel. Lower Bdh1, Ech1, Bnip3 and Decr1, with higher Tpm1, support altered oxidative metabolism and contractile remodeling."),
        SlideNote(11, "Kidney adds a focused signal; other tissues remain exploratory", "1:00", "Kidney supplied the strongest secondary pair: promoted Inpp4b and reinforced Slc37a4. Spleen, skin and adrenal results are narrower. Lung and retina improved predictively without a synthetic-informed BH-FDR gene, while liver, EDL and quadriceps retained real-only models."),
        SlideNote(12, "Synthetic data helped most as a tissue-specific prior", "0:55", "The main result is not that synthetic samples increase biological n. The useful role was tissue-specific feature ranking or limited regularization. The next decisive test is complete-accession holdout followed by validation of the thymus, soleus and kidney hypotheses in independent data."),
        SlideNote(13, "Thank you", "0:10", "Acknowledge the mentor, SLSTP, NASA OSDR, ARCHS4, Reactome and NASA Ames compute. Invite questions."),
    ]
    for note, slide in zip(notes, presentation.slides):
        _add_notes(slide, note)
    _write_notes(notes)
    presentation.save(OUTPUT)
    return OUTPUT


if __name__ == "__main__":
    print(build())
