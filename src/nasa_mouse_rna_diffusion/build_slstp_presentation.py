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
    cleanup = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(0),
        Inches(5.42),
        Inches(4.30),
        Inches(SLIDE_H - 5.42),
    )
    cleanup.name = "Title background artifact mask"
    _set_fill(cleanup, "000000")
    cleanup.line.fill.background()
    shape_tree = slide.shapes._spTree
    shape_tree.remove(cleanup._element)
    shape_tree.insert(2, cleanup._element)

    main = next(
        placeholder
        for placeholder in slide.placeholders
        if placeholder.placeholder_format.type == PP_PLACEHOLDER.BODY
    )
    main.left = Inches(0.62)
    main.top = Inches(2.36)
    main.width = Inches(4.70)
    main.height = Inches(1.55)
    main.text_frame.clear()
    main.text_frame.word_wrap = True
    main.text_frame.margin_left = Inches(0)
    main.text_frame.margin_right = Inches(0)
    main.text_frame.margin_top = Inches(0)
    main.text_frame.margin_bottom = Inches(0)
    main.text_frame.vertical_anchor = MSO_ANCHOR.TOP
    paragraph = main.text_frame.paragraphs[0]
    paragraph.space_before = Pt(0)
    paragraph.space_after = Pt(0)
    paragraph.line_spacing = 1.0
    run = paragraph.add_run()
    run.text = "Synthetic transcriptomics for mouse spaceflight"
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
        ("Can generated RNA-seq improve FLT vs GC analysis within each tissue?", 17, False, "D8E1EB"),
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
    _add_text(slide, "Biological & Physical Sciences", 0.62, 7.12, 3.30, 0.22, size=10.5, color="B9C7D5")


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
        "We built a configurable bulk RNA-seq generation pipeline",
        "Each stage exposes alternatives; a strong outline marks the configuration used downstream.",
    )

    def option_tile(label, x, y, w, color, selected=False, *, size=9.1):
        tile = slide.shapes.add_shape(
            MSO_SHAPE.ROUNDED_RECTANGLE,
            Inches(x),
            Inches(y),
            Inches(w),
            Inches(0.36),
        )
        _set_fill(tile, WHITE if selected else "F7F9FA")
        _set_line(tile, color if selected else "C9D2D8", 2.0 if selected else 0.65)
        text_x = x + (0.15 if selected else 0.06)
        text_w = w - (0.18 if selected else 0.12)
        if selected:
            _add_circle(slide, x + 0.055, y + 0.152, 0.055, color)
        _add_text(
            slide,
            label,
            text_x,
            y + 0.055,
            text_w,
            0.23,
            size=size,
            color=DARK if selected else GRAY,
            bold=selected,
            align=PP_ALIGN.CENTER,
            valign=MSO_ANCHOR.MIDDLE,
            margin=0,
        )

    def segmented_row(x, y, label, options, color):
        _add_text(slide, label, x + 0.16, y, 2.02, 0.18, size=8.4, color=GRAY, bold=True, margin=0)
        gap = 0.05
        available = 2.03
        tile_width = (available - gap * (len(options) - 1)) / len(options)
        for index, (option, selected) in enumerate(options):
            option_tile(
                option,
                x + 0.16 + index * (tile_width + gap),
                y + 0.22,
                tile_width,
                color,
                selected,
                size=8.2 if len(options) == 3 else 8.8,
            )

    def stage_panel(x, number, title, color, fill):
        _add_panel(slide, x, 2.00, 2.35, 3.35, fill=fill, line="DDE4E8")
        _add_circle(slide, x + 0.16, 2.18, 0.32, color)
        _add_text(slide, number, x + 0.16, 2.225, 0.32, 0.20, size=10.5, color=WHITE, bold=True, align=PP_ALIGN.CENTER, margin=0)
        _add_text(slide, title, x + 0.58, 2.19, 1.57, 0.30, size=14.0, color=color, bold=True, valign=MSO_ANCHOR.MIDDLE)
        _add_rule(slide, x + 0.16, 2.58, 2.03, color, 0.025)

    xs = [0.35, 2.90, 5.45, 8.00, 10.55]

    stage_panel(xs[0], "1", "Data scope", BLUE, PALE_BLUE)
    segmented_row(xs[0], 2.72, "Sources", [("OSDR API", True), ("ARCHS4", True)], BLUE)
    segmented_row(xs[0], 3.48, "Studies", [("Single", False), ("Multiple", True)], BLUE)
    segmented_row(xs[0], 4.24, "Tissues", [("Per tissue", False), ("All tissues", True)], BLUE)

    stage_panel(xs[1], "2", "Transform", TEAL, PALE_TEAL)
    segmented_row(xs[1], 2.72, "Expression", [("Raw", False), ("CPM", False), ("TPM", True)], TEAL)
    segmented_row(xs[1], 3.48, "Scaling", [("None", False), ("Z-score", False), ("MaxAbs", True)], TEAL)
    segmented_row(xs[1], 4.24, "Features", [("All", False), ("HVG", False), ("L1000", True)], TEAL)

    stage_panel(xs[2], "3", "Harmonization", PURPLE, "F0ECF6")
    _add_text(slide, "Global or study correction", xs[2] + 0.16, 2.72, 2.02, 0.18, size=8.4, color=GRAY, bold=True, margin=0)
    for index, (label, selected) in enumerate([
        ("None", True),
        ("Within-study z-score", False),
        ("ComBat / MBatch", False),
        ("MOBER", False),
    ]):
        option_tile(label, xs[2] + 0.16, 2.94 + index * 0.50, 2.03, PURPLE, selected, size=9.2)

    stage_panel(xs[3], "4", "Model training", ORANGE, PALE_GOLD)
    segmented_row(xs[3], 2.72, "Generator", [("WGAN-GP", False), ("DDIM", True)], ORANGE)
    _add_text(slide, "Training source", xs[3] + 0.16, 3.48, 2.02, 0.18, size=8.4, color=GRAY, bold=True, margin=0)
    for index, (label, selected) in enumerate([
        ("OSDR only", False),
        ("ARCHS4 only", False),
        ("Pretrain + adapt", True),
    ]):
        option_tile(label, xs[3] + 0.16, 3.70 + index * 0.50, 2.03, ORANGE, selected, size=9.2)

    stage_panel(xs[4], "5", "Conditioning", GREEN, "ECF4ED")
    _add_text(slide, "Model inputs", xs[4] + 0.16, 2.72, 2.02, 0.18, size=8.4, color=GRAY, bold=True, margin=0)
    for index, (label, selected) in enumerate([
        ("Tissue", True),
        ("FLT / GC", True),
        ("Accession", True),
        ("Material type", True),
        ("Sex / age", False),
    ]):
        option_tile(label, xs[4] + 0.16, 2.92 + index * 0.42, 2.03, GREEN, selected, size=8.9)

    _add_panel(slide, 0.35, 5.60, 12.55, 1.20, fill="F7F9FA", line="DDE4E8")
    _add_text(slide, "Selected branch", 0.64, 5.81, 1.42, 0.25, size=10.8, color=GOLD, bold=True)
    _add_text(slide, "Used downstream", 0.64, 6.19, 1.35, 0.24, size=9.6, color=GRAY)

    def selected_step(x, width, heading, value, color):
        _add_panel(slide, x, 5.73, width, 0.49, fill=WHITE, line=color)
        _add_text(slide, heading, x + 0.10, 5.78, width - 0.20, 0.12, size=7.1, color=color, bold=True, margin=0)
        _add_text(slide, value, x + 0.10, 5.94, width - 0.20, 0.19, size=9.1, color=DARK, bold=True, align=PP_ALIGN.CENTER, margin=0)

    selected_path = [
        (2.10, 2.15, "PREPROCESS", "TPM / MaxAbs / 974 landmarks", TEAL),
        (4.65, 1.85, "PRETRAIN", "ARCHS4", TEAL),
        (6.90, 1.85, "ADAPT", "OSDR", BLUE),
        (9.15, 2.45, "GENERATE", "Conditional DDIM", ORANGE),
    ]
    for index, (x, width, heading, value, color) in enumerate(selected_path):
        selected_step(x, width, heading, value, color)
        if index < len(selected_path) - 1:
            _add_arrow(slide, x + width + 0.07, 5.86, 0.23, 0.20, MID_GRAY)

    _add_text(slide, "CONDITION", 2.10, 6.37, 0.72, 0.16, size=7.3, color=GREEN, bold=True, margin=0)
    condition_tiles = [
        (2.88, 0.90, "Tissue"),
        (3.85, 1.00, "FLT / GC"),
        (4.92, 1.14, "Accession"),
        (6.13, 1.38, "Material type"),
    ]
    for x, width, label in condition_tiles:
        option_tile(label, x, 6.27, width, GREEN, True, size=7.8)

    _add_text(slide, "HARMONIZE", 7.76, 6.37, 0.82, 0.16, size=7.3, color=PURPLE, bold=True, margin=0)
    option_tile("None", 8.65, 6.27, 0.78, PURPLE, True, size=7.8)
    _add_text(slide, "SCOPE", 9.70, 6.37, 0.46, 0.16, size=7.3, color=BLUE, bold=True, margin=0)
    option_tile("All tissues", 10.24, 6.27, 1.20, BLUE, True, size=7.8)
    _add_source(slide, "Screen: 463 planned configurations and nine matched liver harmonization arms. Models: Vinas et al. (2022); Lacan et al. (2026).")


def _slide_5(slide, trajectory: Path):
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


def _slide_6(slide):
    _add_slide_title(
        slide,
        "Diffusion output",
        "Generated profiles track the real OSDR PCA manifold",
        "Circles are locked real OSDR profiles; crosses are matched DDIM profiles from generation seed 5020.",
    )
    figure = PAPER_DIR / "figures/figure_2b_locked_real_vs_synthetic_pca.png"
    _add_picture_contain(
        slide,
        figure,
        0.34,
        1.86,
        12.64,
        4.56,
        alt="PCA of locked real and matched DDIM OSDR profiles by tissue and condition",
    )
    _add_panel(slide, 0.48, 6.46, 12.37, 0.48, fill=PALE_BLUE, line="C9DCE9")
    _add_text(
        slide,
        "Interpretation: tissue structure dominates; FLT and GC are subtler. PCA is descriptive, so quantitative metrics determine validation.",
        0.73,
        6.55,
        11.86,
        0.27,
        size=13.3,
        color=NAVY,
        bold=True,
        align=PP_ALIGN.CENTER,
    )
    _add_source(slide, "PCA was fitted in the common locked OSDR expression space. This is PCA, not UMAP.")


def _slide_4(slide):
    _add_slide_title(
        slide,
        "Validation",
        "DDIM matched expression and reduced separability",
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


def _slide_7(slide):
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
        (3.62, 2.10, "Real only", "Real", "Real", BLUE, PALE_BLUE),
        (5.47, 2.10, "Generated only", "Generated", "Generated", CORAL, PALE_CORAL),
        (7.32, 2.10, "Real + generated", "Consensus", "Equal weight", TEAL, PALE_TEAL),
        (4.55, 3.73, "Guided: real fit", "Consensus", "Real", PURPLE, "F0ECF6"),
        (6.40, 3.73, "Guided: 5% synthetic", "Consensus", "Real + 5% synthetic", ORANGE, PALE_GOLD),
    ]
    for x, y, title, rank_value, fit_value, color, fill in cards:
        _add_panel(slide, x, y, 1.65, 1.28, fill=fill, line=color)
        _add_text(slide, title, x + 0.09, y + 0.13, 1.47, 0.34, size=12.4, color=color, bold=True, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE)
        _add_rule(slide, x + 0.13, y + 0.53, 1.39, "D7DEE2", 0.014)
        _add_text(slide, "Rank", x + 0.14, y + 0.66, 0.40, 0.20, size=8.8, color=GRAY, bold=True, margin=0)
        _add_text(slide, rank_value, x + 0.56, y + 0.62, 0.91, 0.26, size=10.6, color=DARK, bold=True, align=PP_ALIGN.RIGHT, valign=MSO_ANCHOR.MIDDLE, margin=0)
        _add_text(slide, "Fit", x + 0.14, y + 0.96, 0.40, 0.20, size=8.8, color=GRAY, bold=True, margin=0)
        _add_text(slide, fit_value, x + 0.56, y + 0.91, 0.91, 0.31, size=9.7, color=DARK, bold=True, align=PP_ALIGN.RIGHT, valign=MSO_ANCHOR.MIDDLE, margin=0)

    _add_arrow(slide, 9.32, 3.12, 0.52, 0.38)
    _add_panel(slide, 10.00, 2.38, 2.88, 2.25, fill="ECF4ED", line=GREEN)
    _add_text(slide, "Held-out real profiles", 10.23, 2.64, 2.42, 0.38, size=16.5, color=GREEN, bold=True, align=PP_ALIGN.CENTER)
    _add_rule(slide, 10.28, 3.11, 2.30, "C9D9CF", 0.018)
    _add_text(slide, "Metrics", 10.30, 3.30, 0.58, 0.20, size=9.3, color=GRAY, bold=True, margin=0)
    _add_text(slide, "BA | AUROC | AP", 10.92, 3.25, 1.60, 0.30, size=12.0, color=NAVY, bold=True, align=PP_ALIGN.RIGHT, valign=MSO_ANCHOR.MIDDLE, margin=0)
    _add_text(slide, "Decision", 10.30, 3.79, 0.58, 0.20, size=9.3, color=GRAY, bold=True, margin=0)
    _add_text(slide, "Eligible arm\nfor each tissue", 10.92, 3.70, 1.60, 0.52, size=11.2, color=DARK, bold=True, align=PP_ALIGN.RIGHT, valign=MSO_ANCHOR.MIDDLE, margin=0)

    _add_panel(slide, 0.38, 5.38, 12.50, 1.15, fill=NAVY, line=NAVY)
    _add_text(slide, "Biology check", 0.70, 5.65, 1.75, 0.34, size=17, color=WHITE, bold=True)
    _add_text(slide, "FLT vs GC effects, random-effects meta-analysis and BH FDR use real OSDR profiles only.", 2.52, 5.66, 7.40, 0.32, size=15.5, color="DCE7F2")
    _add_text(slide, "Synthetic profiles never increase animal n.", 10.05, 5.54, 2.46, 0.56, size=14, color="FFD69A", bold=True, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE)


def _slide_8(slide, utility_chart: Path):
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
    chart_height = 3.00
    chart_left = 0.70
    chart_right = 3.58
    for tick, label in [(0.0, "0"), (0.5, "0.5"), (1.0, "1.0")]:
        tick_y = baseline - tick * chart_height
        _add_rule(slide, chart_left, tick_y, chart_right - chart_left, "D8DEE2", 0.012)
        _add_text(
            slide,
            label,
            0.43,
            tick_y - 0.10,
            0.22,
            0.20,
            size=8.5,
            color=MID_GRAY,
            align=PP_ALIGN.RIGHT,
            margin=0,
        )
    for index, (label, value, color) in enumerate(zip(labels, values, colors)):
        x = 0.78 + index * 0.96
        bar_h = value * chart_height
        shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(baseline - bar_h), Inches(0.62), Inches(bar_h))
        _set_fill(shape, color)
        shape.line.fill.background()
        _add_text(slide, f"{value:.3f}", x - 0.08, baseline - bar_h - 0.35, 0.78, 0.28, size=14, color=NAVY, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, label, x - 0.20, baseline + 0.08, 1.05, 0.55, size=10.5, color=DARK, align=PP_ALIGN.CENTER)
    _add_rule(slide, chart_left, baseline, chart_right - chart_left, "86949D", 0.02)
    _add_text(slide, "Balanced accuracy (0-1)", 0.72, 2.78, 2.78, 0.25, size=12.5, color=GRAY)

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


def _slide_9(slide):
    _add_slide_title(
        slide,
        "Feature analysis",
        "Synthetic guidance changed ranking, not statistical evidence",
        "Promoted and reinforced describe repeated feature selection; both require a supporting effect in real OSDR data.",
    )
    _add_text(slide, "Repeated stable selection", 0.62, 2.02, 4.86, 0.34, size=16, color=GRAY, bold=True, align=PP_ALIGN.CENTER)
    left = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.62), Inches(2.52), Inches(3.00), Inches(2.35))
    left.fill.background()
    _set_line(left, BLUE, 2.4)
    right = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(2.48), Inches(2.52), Inches(3.00), Inches(2.35))
    right.fill.background()
    _set_line(right, TEAL, 2.4)
    # Separate line boxes prevent PowerPoint from wrapping one Venn label
    # differently from the other when fonts are substituted during rendering.
    _add_text(slide, "Real-only", 0.71, 3.08, 1.76, 0.29, size=16, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "ranking", 0.71, 3.40, 1.76, 0.29, size=16, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "Synthetic-guided", 3.61, 3.08, 1.88, 0.29, size=13.5, color=TEAL, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "ranking", 3.61, 3.40, 1.88, 0.29, size=16, color=TEAL, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "23", 2.52, 3.28, 1.06, 0.24, size=12, color=ORANGE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "reinforced", 2.52, 3.54, 1.06, 0.24, size=12, color=ORANGE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_panel(slide, 3.52, 4.96, 1.86, 0.78, fill=PALE_TEAL, line=TEAL)
    _add_text(slide, "26 promoted", 3.66, 5.06, 1.58, 0.24, size=13.5, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "selected only with\nsynthetic guidance", 3.68, 5.34, 1.54, 0.32, size=9.5, color=GRAY, align=PP_ALIGN.CENTER)

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


def _slide_10(slide):
    _add_slide_title(
        slide,
        "Literature interpretation",
        "Selection and literature are separate dimensions",
        "All 49 synthetic-informed associations received both a selection status and a literature interpretation.",
    )

    steps = [
        ("1", "Candidate", "Fix", "Gene, tissue, and FLT direction", BLUE, PALE_BLUE),
        ("2", "Search", "Query", "Spaceflight and mechanism literature", TEAL, PALE_TEAL),
        ("3", "Classify", "Match", "Result with published evidence", ORANGE, PALE_GOLD),
        ("4", "Interpret", "Record", "Category, rationale, and source IDs", PURPLE, "F0ECF6"),
    ]
    for index, (number, heading, action, detail, color, fill) in enumerate(steps):
        x = 0.45 + index * 3.12
        _add_panel(slide, x, 2.00, 2.55, 0.98, fill=fill, line=color)
        _add_circle(slide, x + 0.15, 2.14, 0.31, color)
        _add_text(slide, number, x + 0.15, 2.185, 0.31, 0.19, size=10.4, color=WHITE, bold=True, align=PP_ALIGN.CENTER, margin=0)
        _add_text(slide, heading, x + 0.58, 2.13, 1.73, 0.28, size=14.0, color=color, bold=True, valign=MSO_ANCHOR.MIDDLE)
        _add_rule(slide, x + 0.15, 2.51, 2.25, "D7DEE2", 0.014)
        _add_text(slide, action, x + 0.16, 2.65, 0.52, 0.19, size=8.7, color=GRAY, bold=True, margin=0)
        _add_text(slide, detail, x + 0.74, 2.59, 1.64, 0.31, size=9.8, color=DARK, bold=True, align=PP_ALIGN.RIGHT, valign=MSO_ANCHOR.MIDDLE, margin=0)
        if index < len(steps) - 1:
            _add_arrow(slide, x + 2.67, 2.31, 0.31, 0.28, MID_GRAY)

    annotations = pd.read_csv(
        PAPER_DIR / "source_data/table_s16_promoted_gene_literature_annotations.tsv",
        sep="\t",
    )
    counts = annotations["literature_classification"].value_counts().to_dict()
    selection_counts = annotations["selection_status"].value_counts().to_dict()
    direct = int(
        annotations["evidence_scope"]
        .eq("direct_same_gene_same_tissue_same_direction")
        .sum()
    )
    if len(annotations) != 49 or direct != 5:
        raise ValueError("Unexpected synthetic-informed literature inventory")

    _add_panel(slide, 0.45, 3.19, 4.14, 3.26, fill=PALE_BLUE, line="C9D9EA")
    _add_text(slide, "Selection status", 0.73, 3.43, 3.58, 0.34, size=17, color=BLUE, bold=True)
    selection_rows = [
        ("Promoted", selection_counts.get("promoted", 0), "Stable only with synthetic guidance", CORAL),
        ("Reinforced", selection_counts.get("reinforced", 0), "Stable with and without guidance", TEAL),
    ]
    for index, (label, count, detail, color) in enumerate(selection_rows):
        y = 4.08 + index * 0.93
        _add_circle(slide, 0.75, y, 0.50, color)
        _add_text(slide, str(count), 0.75, y + 0.08, 0.50, 0.25, size=14, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, label, 1.45, y - 0.02, 2.65, 0.30, size=15, color=color, bold=True)
        _add_text(slide, detail, 1.45, y + 0.31, 2.65, 0.38, size=11.5, color=DARK)

    _add_panel(slide, 4.80, 3.19, 8.07, 3.26, fill=WHITE, line="DDE4E8")
    _add_text(slide, "Literature interpretation", 5.08, 3.43, 7.52, 0.34, size=17, color=NAVY, bold=True)
    category_rows = [
        ("Aligning", counts.get("aligning", 0), "Direct or same-tissue process agreement", GREEN),
        ("Complementary", counts.get("complementary", 0), "Related process or mechanism", BLUE),
        ("Ambiguous", counts.get("ambiguous", 0), "Mixed or context-dependent evidence", ORANGE),
        ("Unmatched", counts.get("unmatched", 0), "No sufficiently specific match found", PURPLE),
    ]
    for index, (label, count, detail, color) in enumerate(category_rows):
        y = 3.91 + index * 0.57
        _add_circle(slide, 5.09, y + 0.02, 0.37, color)
        _add_text(slide, str(count), 5.09, y + 0.07, 0.37, 0.22, size=12, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, label, 5.62, y, 1.74, 0.34, size=14.2, color=color, bold=True, valign=MSO_ANCHOR.MIDDLE)
        _add_text(slide, detail, 7.42, y, 4.93, 0.34, size=12.8, color=DARK, valign=MSO_ANCHOR.MIDDLE)

    _add_panel(slide, 0.45, 6.57, 12.42, 0.43, fill=NAVY, line=NAVY)
    _add_text(
        slide,
        "A reinforced gene can be complementary, ambiguous, or unmatched; an aligning gene can be promoted.",
        0.71,
        6.64,
        11.90,
        0.25,
        size=14,
        color=WHITE,
        bold=True,
        align=PP_ALIGN.CENTER,
    )
    _add_source(
        slide,
        "Table S16 records gene-level rationale, evidence scope and source IDs; Table S17 provides 33 citations with DOI/URL and data relationship.",
    )


def _inventory_gene_segments(
    rows: pd.DataFrame,
    annotation_lookup: dict[tuple[str, str, str], str],
) -> list[tuple[str, dict]]:
    class_colors = {
        "aligning": GREEN,
        "complementary": BLUE,
        "ambiguous": ORANGE,
        "unmatched": PURPLE,
    }
    ordered = rows.sort_values(["real_meta_fdr", "symbol"])
    segments: list[tuple[str, dict]] = []
    for index, row in enumerate(ordered.itertuples(index=False)):
        if index:
            segments.append((", ", {"size": 9.0, "color": MID_GRAY}))
        key = (row.analysis_scope, row.tissue, row.gene)
        literature_class = annotation_lookup.get(key)
        if literature_class not in class_colors:
            raise ValueError(f"Missing synthetic-informed annotation for {key}")
        color = class_colors[literature_class]
        segments.append(
            (
                row.symbol,
                {
                    "size": 9.0,
                    "color": color,
                    "bold": True,
                    "italic": True,
                },
            )
        )
    return segments


def _add_gene_inventory_block(
    slide,
    inventory: pd.DataFrame,
    annotation_lookup: dict[tuple[str, str, str], str],
    *,
    scope: str,
    tissue: str,
    label: str,
    x: float,
    y: float,
    w: float,
    h: float,
) -> None:
    rows = inventory[
        inventory["analysis_scope"].eq(scope) & inventory["tissue"].eq(tissue)
    ]
    if rows.empty:
        raise ValueError(f"No synthetic-informed genes for {scope}/{tissue}")

    _add_rule(slide, x, y, w, "DDE4E8", 0.018)
    count_label = "gene" if len(rows) == 1 else "genes"
    label_size = 10.8 if len(label.replace("\n", "")) > 12 else 11.8
    _add_rich_text(
        slide,
        [
            (f"{label}\n", {"size": label_size, "color": NAVY, "bold": True}),
            (f"{len(rows)} {count_label}", {"size": 9.0, "color": MID_GRAY}),
        ],
        x + 0.06,
        y + 0.07,
        1.36,
        h - 0.10,
        margin=0,
    )

    groups = []
    for direction, direction_label, direction_color in [
        ("FLT_higher", "FLT higher", CORAL),
        ("FLT_lower", "FLT lower", TEAL),
    ]:
        for selection, selection_label, selection_color in [
            ("synthetic_promoted", "Promoted", CORAL),
            ("reinforced_real_and_synthetic", "Reinforced", TEAL),
        ]:
            group_rows = rows[
                rows["flt_gc_direction"].eq(direction)
                & rows["selection_interpretation"].eq(selection)
            ]
            if group_rows.empty:
                continue
            symbol_characters = int(group_rows["symbol"].str.len().sum()) + 2 * max(0, len(group_rows) - 1)
            estimated_lines = max(1, int(np.ceil(symbol_characters / 38)))
            groups.append(
                (
                    direction_label,
                    direction_color,
                    selection_label,
                    selection_color,
                    group_rows,
                    1.0 + 0.75 * (estimated_lines - 1),
                )
            )

    available_height = h - 0.11
    weight_total = sum(group[-1] for group in groups)
    row_y = y + 0.06
    for index, (direction_label, direction_color, selection_label, selection_color, group_rows, weight) in enumerate(groups):
        row_h = available_height * weight / weight_total
        if index:
            _add_rule(slide, x + 1.44, row_y, w - 1.50, "E2E7EA", 0.010)
        _add_text(
            slide,
            direction_label,
            x + 1.47,
            row_y + 0.01,
            0.72,
            row_h - 0.02,
            size=8.0,
            color=direction_color,
            bold=True,
            valign=MSO_ANCHOR.MIDDLE,
            margin=0,
        )
        _add_text(
            slide,
            selection_label,
            x + 2.22,
            row_y + 0.01,
            0.78,
            row_h - 0.02,
            size=7.9,
            color=selection_color,
            bold=True,
            valign=MSO_ANCHOR.MIDDLE,
            margin=0,
        )
        _add_rich_text(
            slide,
            _inventory_gene_segments(group_rows, annotation_lookup),
            x + 3.04,
            row_y + 0.01,
            w - 3.12,
            row_h - 0.02,
            size=9.0,
            valign=MSO_ANCHOR.MIDDLE,
            margin=0,
        )
        row_y += row_h


def _slide_all_tissue_coverage(slide):
    _add_slide_title(
        slide,
        "Analysis coverage",
        "All 27 completed tissue analyses were retained",
        "The biological narrative narrows only after reporting synthetic-informed, real-only, and null outcomes.",
    )
    summary = pd.read_csv(
        PAPER_DIR / "source_data/table_s12_bh_fdr_tissue_summary.tsv",
        sep="\t",
    )
    if len(summary) != 27:
        raise ValueError("Expected 27 completed tissue analysis units")

    synthetic = summary.loc[
        summary["n_reinforced_real_and_synthetic"]
        .add(summary["n_synthetic_promoted"])
        .gt(0)
    ]
    real_only = summary.loc[
        summary["n_bh_fdr_genes"].gt(0)
        & ~summary.index.to_series().isin(synthetic.index)
    ]
    null = summary.loc[summary["n_bh_fdr_genes"].eq(0)]
    if (len(synthetic), len(real_only), len(null)) != (10, 5, 12):
        raise ValueError("Unexpected all-tissue coverage partition")

    display = {
        "skeletal_muscle": "Skeletal muscle (pooled)",
        "edl": "EDL",
    }

    def names(frame: pd.DataFrame) -> list[str]:
        return [
            display.get(value, value.replace("_", " ").title())
            for value in frame["tissue"]
        ]

    columns = [
        (
            0.43,
            "10",
            "Synthetic-informed",
            names(synthetic),
            TEAL,
            PALE_TEAL,
            "Promoted or reinforced BH-FDR gene",
        ),
        (
            4.49,
            "5",
            "Real BH-FDR only",
            names(real_only),
            BLUE,
            PALE_BLUE,
            "Real association, no synthetic-informed selection",
        ),
        (
            8.55,
            "12",
            "No BH-FDR gene",
            names(null),
            ORANGE,
            PALE_GOLD,
            "No BH-FDR gene in the 974-gene panel",
        ),
    ]
    for x, count, heading, tissue_names, color, fill, detail in columns:
        _add_panel(slide, x, 2.02, 3.83, 4.70, fill=fill, line=color)
        _add_circle(slide, x + 0.25, 2.25, 0.52, color)
        _add_text(slide, count, x + 0.25, 2.34, 0.52, 0.25, size=14.5, color=WHITE, bold=True, align=PP_ALIGN.CENTER, margin=0)
        _add_text(slide, heading, x + 0.92, 2.27, 2.58, 0.38, size=16.3, color=color, bold=True, valign=MSO_ANCHOR.MIDDLE)
        _add_rule(slide, x + 0.28, 2.83, 3.27, "D7DEE2", 0.014)
        _add_text(slide, "Criterion", x + 0.30, 3.00, 0.76, 0.18, size=8.7, color=GRAY, bold=True, margin=0)
        _add_text(slide, detail, x + 1.08, 2.94, 2.42, 0.39, size=10.2, color=DARK, bold=True, align=PP_ALIGN.RIGHT, valign=MSO_ANCHOR.MIDDLE, margin=0)
        _add_rule(slide, x + 0.28, 3.48, 3.27, "D7DEE2", 0.014)
        _add_text(slide, "Analysis units", x + 0.30, 3.66, 1.20, 0.19, size=8.8, color=GRAY, bold=True, margin=0)

        list_top = 3.98
        list_height = 2.42
        row_height = min(0.36, list_height / len(tissue_names))
        list_size = 9.8 if len(tissue_names) >= 12 else 10.4 if len(tissue_names) >= 10 else 11.2
        for index, tissue_name in enumerate(tissue_names):
            row_y = list_top + index * row_height
            _add_circle(slide, x + 0.31, row_y + (row_height - 0.06) / 2, 0.06, color)
            _add_text(
                slide,
                tissue_name,
                x + 0.48,
                row_y,
                3.02,
                row_height,
                size=list_size,
                color=DARK,
                valign=MSO_ANCHOR.MIDDLE,
                margin=0,
            )
    _add_source(slide, "Twenty-two canonical tissues plus EDL, gastrocnemius, quadriceps, soleus and tibialis anterior; complete counts in Table S12.")


def _slide_11(slide):
    _add_slide_title(
        slide,
        "Gene inventory",
        "Synthetic-informed genes spanned 10 tissue analyses",
        "All 49 associations passed BH FDR in real OSDR data; synthetic profiles affected prioritization, not the statistical test.",
    )
    inventory = pd.read_csv(
        PAPER_DIR / "source_data/table_s10_synthetic_informed_bh_fdr_genes.tsv",
        sep="\t",
    )
    annotations = pd.read_csv(
        PAPER_DIR / "source_data/table_s16_promoted_gene_literature_annotations.tsv",
        sep="\t",
    )
    if len(inventory) != 49 or len(annotations) != 49:
        raise ValueError("Unexpected synthetic-informed gene inventory")
    if inventory[["analysis_scope", "tissue"]].drop_duplicates().shape[0] != 10:
        raise ValueError("Expected 10 tissue analyses in synthetic-informed inventory")
    annotation_lookup = {
        (row.analysis_scope, row.tissue, row.gene): row.literature_classification
        for row in annotations.itertuples(index=False)
    }

    _add_text(slide, "Rows: FLT direction + status", 0.45, 1.86, 2.42, 0.25, size=9.8, color=NAVY, bold=True)
    legend = [
        ("aligning", GREEN, 3.02, 0.84),
        ("complementary", BLUE, 4.02, 1.15),
        ("ambiguous", ORANGE, 5.37, 0.96),
        ("unmatched", PURPLE, 6.53, 0.96),
    ]
    for label, color, x, width in legend:
        _add_rule(slide, x, 1.91, 0.14, color, 0.05)
        _add_text(slide, label, x + 0.19, 1.85, width, 0.25, size=9.8, color=color, bold=True)
    _add_text(slide, "Gene color = literature | Non-empty rows shown", 8.64, 1.85, 4.20, 0.25, size=9.8, color=MID_GRAY, italic=True, align=PP_ALIGN.RIGHT)

    left_blocks = [
        ("canonical_tissue", "thymus", "Thymus", 2.24, 1.45),
        ("canonical_tissue", "spleen", "Spleen", 3.78, 0.67),
        ("canonical_tissue", "kidney", "Kidney", 4.54, 0.61),
        ("skeletal_muscle_group", "gastrocnemius", "Gastrocnemius", 5.24, 0.61),
        ("canonical_tissue", "eye", "Eye", 5.94, 0.51),
    ]
    right_blocks = [
        ("canonical_tissue", "skeletal_muscle", "Skeletal muscle\n(pooled)", 2.24, 1.25),
        ("skeletal_muscle_group", "soleus", "Soleus", 3.58, 0.68),
        ("skeletal_muscle_group", "tibialis_anterior", "Tibialis anterior", 4.35, 0.69),
        ("canonical_tissue", "adrenal_gland", "Adrenal gland", 5.13, 0.61),
        ("canonical_tissue", "skin", "Skin", 5.83, 0.55),
    ]
    for scope, tissue, label, y, h in left_blocks:
        _add_gene_inventory_block(
            slide,
            inventory,
            annotation_lookup,
            scope=scope,
            tissue=tissue,
            label=label,
            x=0.43,
            y=y,
            w=6.15,
            h=h,
        )
    for scope, tissue, label, y, h in right_blocks:
        _add_gene_inventory_block(
            slide,
            inventory,
            annotation_lookup,
            scope=scope,
            tissue=tissue,
            label=label,
            x=6.76,
            y=y,
            w=6.15,
            h=h,
        )
    _add_source(slide, "Separate rows report FLT direction and selection status; gene color independently reports literature interpretation for all 49 associations.")


def _slide_12(slide):
    _add_slide_title(
        slide,
        "Thymus",
        "Thymus points to lower proliferative renewal",
        "The core panel recovers cell-cycle biology; two flight-higher genes extend the hypothesis.",
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
            "Higher Hsd17b11 and Etv1 suggest lipid or T-cell-state shifts",
            "Bulk RNA-seq cannot separate cell loss from lower transcription",
        ],
        9.44,
        3.53,
        3.03,
        size=12.2,
        bullet_color=CORAL,
        row_h=0.64,
    )
    _add_source(slide, "Context: Gridley et al. (2013), Horie et al. (2019), Keenan et al. (2025), and Shi et al. (2026).")


def _slide_13(slide):
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


def _add_additional_finding_rows(slide, rows):
    _add_text(slide, "Analysis unit", 0.72, 2.00, 1.70, 0.28, size=11.5, color=GRAY, bold=True)
    _add_text(slide, "Selection", 2.48, 2.00, 1.15, 0.28, size=11.5, color=GRAY, bold=True)
    _add_text(slide, "Genes and FLT direction", 3.78, 2.00, 3.15, 0.28, size=11.5, color=GRAY, bold=True)
    _add_text(slide, "Interpretation", 7.20, 2.00, 5.10, 0.28, size=11.5, color=GRAY, bold=True)
    for index, (tissue, promoted, reinforced, interpretation, color, fill) in enumerate(rows):
        y = 2.35 + index * 1.06
        _add_panel(slide, 0.43, y, 12.44, 0.88, fill=fill, line="DDE4E8", radius=False)
        _add_rule(slide, 0.43, y, 0.08, color, 0.88)
        _add_text(slide, tissue, 0.72, y + 0.20, 1.62, 0.44, size=15.0, color=color, bold=True, valign=MSO_ANCHOR.MIDDLE)
        _add_rule(slide, 2.42, y + 0.44, 4.42, "DDE4E8", 0.012)
        _add_text(slide, "Promoted", 2.48, y + 0.08, 1.12, 0.25, size=10.5, color=CORAL, bold=True, valign=MSO_ANCHOR.MIDDLE)
        _add_text(slide, promoted, 3.78, y + 0.05, 3.08, 0.35, size=10.5, color=DARK, valign=MSO_ANCHOR.MIDDLE)
        _add_text(slide, "Reinforced", 2.48, y + 0.52, 1.12, 0.25, size=10.5, color=TEAL, bold=True, valign=MSO_ANCHOR.MIDDLE)
        _add_text(slide, reinforced, 3.78, y + 0.48, 3.08, 0.35, size=10.5, color=DARK, valign=MSO_ANCHOR.MIDDLE)
        _add_text(slide, interpretation, 7.20, y + 0.10, 5.20, 0.66, size=12.2, color=DARK, valign=MSO_ANCHOR.MIDDLE)


def _slide_14(slide):
    _add_slide_title(
        slide,
        "Additional findings I",
        "Each additional tissue defines a separate hypothesis",
        "Promoted and reinforced selections are separated within every tissue.",
    )
    rows = [
        (
            "Pooled muscle",
            "Lower: Klhl21, Mapkapk5, Reep5, Itgb5",
            "Higher: Sox4, Cebpd, Sh3bp5, Prkcd, Arid5b, Sesn1, Tle1; lower: Bphl",
            "Heterogeneous stress, differentiation, interferon and sialic-acid response; interpret with the anatomical muscle groups.",
            BLUE,
            PALE_BLUE,
        ),
        (
            "Kidney",
            "Inpp4b higher",
            "Slc37a4 higher",
            "Renal phosphoinositide signaling and glucose-handling hypothesis; the pair had no shared Reactome enrichment.",
            PURPLE,
            "F0ECF6",
        ),
        (
            "Spleen",
            "Rai14, Myl9, Ptprk higher",
            "Loxl1 higher",
            "Adhesion, actomyosin and extracellular-matrix or immune-organization hypothesis; no coherent pathway enrichment.",
            TEAL,
            PALE_TEAL,
        ),
        (
            "Skin",
            "Plscr1 higher",
            "None",
            "Interferon-linked skin candidate within broader cell-cycle and DNA-repair responses; a single-gene result.",
            ORANGE,
            PALE_GOLD,
        ),
    ]
    _add_additional_finding_rows(slide, rows)
    _add_source(slide, "All listed effects passed BH FDR in observed OSDR profiles; 'None' means no gene in that selection category.")


def _slide_additional_2(slide):
    _add_slide_title(
        slide,
        "Additional findings II",
        "Eye, adrenal and muscle-group results remain tissue-specific",
        "Promoted and reinforced selections are separated within every tissue.",
    )
    rows = [
        (
            "Eye",
            "None",
            "Klhl21 lower",
            "Process-level alignment with lower proliferation and cytokinesis in flight eye tissue; not an exact prior gene replication.",
            GREEN,
            "ECF4ED",
        ),
        (
            "Adrenal gland",
            "Psmb8 lower",
            "Tspan4 lower",
            "Literature-unmatched immune, proteostasis or tissue-composition candidates; no adrenal mechanism is established.",
            MID_GRAY,
            "F2F4F5",
        ),
        (
            "Gastrocnemius",
            "Nfkbia higher; Fhl2 lower",
            "None",
            "NF-kappaB stress alignment plus an autophagy or myogenesis candidate; the two genes do not form a coherent pathway.",
            BLUE,
            PALE_BLUE,
        ),
        (
            "Tibialis anterior",
            "Cebpd higher",
            "Cdkn1a, St3gal5, Bnip3 higher",
            "Stress, cell-cycle, ganglioside-signaling and mitophagy candidates with mixed or complementary prior evidence.",
            TEAL,
            PALE_TEAL,
        ),
    ]
    _add_additional_finding_rows(slide, rows)
    _add_source(slide, "Interpretations follow the independent literature annotations in Table S16; 'None' means no gene in that selection category.")


def _slide_15(slide):
    _add_slide_title(
        slide,
        "Takeaways",
        "Synthetic data helped most as a tissue-specific prior",
        "The generator created useful structure, not new biological replication.",
    )
    columns = [
        (0.42, "1", "Generate", "Conditional DDIM produced high-fidelity profiles with near-chance real-versus-synthetic separation.", BLUE, PALE_BLUE),
        (4.44, "2", "Use", "Pooled augmentation failed. Tissue-specific ranking and low-weight training were more useful.", TEAL, PALE_TEAL),
        (8.46, "3", "Interpret", "Thymus and soleus were clearest. Literature review separated recovery from complementary hypotheses.", ORANGE, PALE_GOLD),
    ]
    for x, number, heading, body, color, fill in columns:
        _add_panel(slide, x, 2.10, 3.74, 2.72, fill=fill, line=color)
        _add_circle(slide, x + 0.26, 2.37, 0.52, color)
        _add_text(slide, number, x + 0.26, 2.43, 0.52, 0.30, size=16, color=WHITE, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, heading, x + 0.91, 2.38, 2.45, 0.36, size=19, color=color, bold=True)
        _add_text(slide, body, x + 0.28, 3.08, 3.15, 1.33, size=15, color=DARK)
    _add_panel(slide, 0.42, 5.15, 12.00, 1.30, fill=NAVY, line=NAVY)
    _add_text(slide, "Next test", 0.75, 5.48, 1.38, 0.34, size=17, color="FFD69A", bold=True)
    _add_text(slide, "Use independent samples and cell-resolved assays to test thymus proliferation, soleus metabolism and prioritized candidates from additional tissues.", 2.05, 5.40, 9.92, 0.58, size=16, color=WHITE, valign=MSO_ANCHOR.MIDDLE)
    _add_text(slide, "Generated profiles never enter the biological sample count or the BH-FDR test.", 2.05, 6.05, 9.92, 0.28, size=12.5, color="BFD0E1")


def _slide_16(slide):
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
        "Target length: 12-15 minutes. Planned speaking time: about 14 minutes 30 seconds.",
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
    while len(presentation.slides) < 18:
        number = len(presentation.slides) + 1
        slide = presentation.slides.add_slide(presentation.slide_layouts[3])
        _prepare_content_slide(slide, number)

    builders = [
        None,
        _slide_2,
        _slide_3,
        _slide_4,
        lambda slide: _slide_5(slide, trajectory),
        _slide_6,
        _slide_7,
        lambda slide: _slide_8(slide, utility_chart),
        _slide_9,
        _slide_10,
        _slide_all_tissue_coverage,
        _slide_11,
        _slide_12,
        _slide_13,
        _slide_14,
        _slide_additional_2,
        _slide_15,
        _slide_16,
    ]
    for index, builder in enumerate(builders):
        if builder is not None:
            builder(presentation.slides[index])

    notes = [
        SlideNote(1, "Synthetic transcriptomics for mouse spaceflight", "0:25", "Introduce the question. This talk asks whether generated expression can help analyze tissue-specific FLT versus GC biology without counting synthetic profiles as additional animals."),
        SlideNote(2, "Small studies and study effects complicate tissue comparisons", "0:50", "OSDR gives broad tissue coverage, but the data are spread across 75 accessions with different mission and assay contexts. ARCHS4 supplies a much larger mouse reference. The challenge is to use that reference without letting study structure masquerade as spaceflight biology."),
        SlideNote(3, "We built a configurable bulk RNA-seq generation pipeline", "1:10", "Each card shows alternatives available at one pipeline stage, and the heavy outlines identify the downstream branch. We used both ARCHS4 and NASA OSDR across multiple studies and all tissues. The selected path used TPM, training-fitted MaxAbs scaling, 974 mouse landmarks, no global correction, ARCHS4 pretraining plus OSDR adaptation, and a DDIM conditioned on tissue, FLT/GC, accession and material. WGAN-GP and the other preprocessing and harmonization choices remained benchmark alternatives."),
        SlideNote(4, "DDIM matched expression and reduced separability", "0:55", "Both WGAN-GP and DDIM had high correlation and F1. DDIM had adversarial accuracy near 0.5 and a lower Frechet-distance ratio, so it was harder to separate from real profiles and closer in distribution. The metrics use each model's stated evaluation split. DDIM was used for the remaining analyses."),
        SlideNote(5, "Diffusion learns tissue structure from noise", "0:50", "Read the panels from left to right. The same generated profiles begin as noise, develop structure by timestep 200 and approach tissue-conditioned regions at timestep zero. The axes are shared, so the visual change is not caused by rescaling each panel."),
        SlideNote(6, "Generated profiles track the real OSDR PCA manifold", "0:50", "Circles are locked real OSDR profiles and crosses are matched DDIM profiles in the same PCA space. Generated samples track the tissue-defined branches. FLT and GC remain more interspersed because condition effects are smaller than tissue effects. Visual overlap complements the quantitative validation."),
        SlideNote(7, "Synthetic profiles entered the analysis in five different ways", "0:55", "Synthetic data can be used for direct training, mixed training or feature guidance. Each tissue could choose among five arms. The eligibility check used held-out real profiles. Once features were nominated, FLT/GC effects and BH FDR were computed from observed OSDR samples only."),
        SlideNote(8, "Pooling tissues hid useful signal", "0:55", "The pooled augmentation test was negative: balanced accuracy fell from 0.754 to 0.737 with real plus synthetic training. The bars use a true zero-to-one balanced-accuracy scale. Tissue-specific analysis changed the result. Different tissues benefited from different synthetic uses, which argues against one global augmentation policy."),
        SlideNote(9, "Synthetic guidance changed ranking, not statistical evidence", "0:50", "Reinforced genes were selected with and without synthetic guidance. Promoted genes crossed the stable-selection rule only with synthetic guidance. Promoted does not mean biologically novel. All 49 synthetic-informed tissue-gene associations also had BH FDR below 0.05 in real data."),
        SlideNote(10, "Selection and literature are separate dimensions", "0:55", "Every association has two labels. Promoted or reinforced describes repeated feature selection. Aligning, complementary, ambiguous or unmatched describes prior literature. Across all 49 associations, 22 aligned, 19 were complementary, four were ambiguous and four were unmatched. Table S16 records the gene-level rationale, evidence scope and source IDs; Table S17 records citations, DOI or URL and whether the evidence is independent, overlapping or mechanistic context."),
        SlideNote(11, "All 27 completed tissue analyses were retained", "0:45", "This is the full analysis coverage: 22 canonical tissues and five anatomical muscle groups. Ten had a synthetic-informed BH-FDR association, five had real BH-FDR genes without synthetic-informed selection, and 12 had no BH-FDR gene in the landmark panel. The narrative focuses later, but no completed tissue result is hidden."),
        SlideNote(12, "Synthetic-informed genes spanned 10 tissue analyses", "0:45", "This is the complete 49-association inventory. Separate rows show FLT-higher or FLT-lower direction and promoted or reinforced selection status. Gene color independently shows aligning, complementary, ambiguous or unmatched literature. FLT directions come from real-data meta-analysis."),
        SlideNote(13, "Thymus points to lower proliferative renewal", "1:10", "Thymus produced the clearest promoted panel. The lower mitotic and DNA-replication genes agree with prior reports of thymic involution and altered cell-cycle expression after flight. Higher Hsd17b11 and Etv1 add lipid-handling and T-cell-state hypotheses. Because this is bulk RNA-seq, the pattern may reflect transcription, cell composition or both."),
        SlideNote(14, "Soleus reinforces a mitochondrial and lipid program", "1:00", "Soleus improved with real plus generated training. The selected genes were already stable in real-only analysis, so synthetic data reinforced rather than introduced the panel. Lower Bdh1, Ech1, Bnip3 and Decr1, with higher Tpm1, support altered oxidative metabolism and contractile remodeling. The literature is mixed for Bnip3 and Tpm1, which is recorded explicitly."),
        SlideNote(15, "Each additional tissue defines a separate hypothesis", "0:40", "Promoted and reinforced genes are shown on separate subrows for each tissue. Pooled muscle, kidney, spleen and skin each produced a distinct synthetic-informed result. The rows share a slide for presentation space; each remains a separate hypothesis. Pooled muscle is heterogeneous, kidney suggests phosphoinositide and glucose handling, spleen suggests adhesion and extracellular-matrix or immune organization, and skin contributes a single interferon-linked candidate."),
        SlideNote(16, "Eye, adrenal and muscle-group results remain tissue-specific", "0:40", "Promoted and reinforced genes remain separated here as well. Eye reinforces lower cytokinesis, adrenal contributes two unmatched candidates, gastrocnemius combines an NF-kappa-B stress signal with an autophagy or myogenesis candidate, and tibialis anterior spans stress, cell-cycle, ganglioside and mitophagy hypotheses."),
        SlideNote(17, "Synthetic data helped most as a tissue-specific prior", "0:40", "Synthetic data was useful for tissue-specific feature ranking and limited regularization. It did not increase biological sample size. Literature annotation separated exact recovery, process-level agreement and complementary hypotheses. Independent samples and cell-resolved experiments are the next tests."),
        SlideNote(18, "Thank you", "0:10", "Acknowledge the mentor, SLSTP, NASA OSDR, ARCHS4, Reactome and NASA Ames compute. Invite questions."),
    ]
    for note, slide in zip(notes, presentation.slides):
        _add_notes(slide, note)
    _write_notes(notes)
    presentation.save(OUTPUT)
    return OUTPUT


if __name__ == "__main__":
    print(build())
