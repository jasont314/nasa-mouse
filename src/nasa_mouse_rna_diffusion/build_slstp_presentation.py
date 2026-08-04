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
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE, PP_PLACEHOLDER
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.xmlchemy import OxmlElement
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
        color_nodes = shape._element.spPr.xpath("./a:solidFill/a:srgbClr")
        if color_nodes:
            alpha = OxmlElement("a:alpha")
            alpha.set("val", str(int((100 - transparency) * 1000)))
            color_nodes[0].append(alpha)


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


def _add_data_badge(slide, label: str, x: float, y: float, color: str, diameter: float = 0.24):
    _add_circle(slide, x, y, diameter, color)
    _add_text(
        slide,
        label,
        x,
        y + 0.005,
        diameter,
        diameter - 0.01,
        size=8.2,
        color=WHITE,
        bold=True,
        align=PP_ALIGN.CENTER,
        valign=MSO_ANCHOR.MIDDLE,
        margin=0,
    )


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
        _add_circle(slide, x, top + 0.10, 0.10, bullet_color)
        _add_text(
            slide,
            row,
            x + 0.20,
            top,
            w - 0.20,
            row_h,
            size=size,
            color=color,
            valign=MSO_ANCHOR.TOP,
            margin=0,
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
            placeholder.text = ""
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
    main.text_frame.word_wrap = True
    main.text_frame.vertical_anchor = MSO_ANCHOR.BOTTOM
    paragraph = main.text_frame.paragraphs[0]
    paragraph.space_before = Pt(0)
    paragraph.space_after = Pt(0)
    paragraph.line_spacing = 1.0
    run = paragraph.add_run()
    run.text = "Synthetic transcriptomics\nfor mouse spaceflight"
    run.font.name = FONT
    run.font.size = Pt(29)
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
        ("Jason Trinh", 21, True, WHITE),
        ("EECS | UC Berkeley", 14, False, WHITE),
        ("", 6, False, WHITE),
        ("Mentor: James Casaletto | August 2026", 13, False, WHITE),
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


def _slide_why_synthetic(slide):
    _add_slide_title(
        slide,
        "Background",
        "What is synthetic transcriptomics?",
        "A generator learns patterns in measured RNA-seq data, then samples new expression vectors under chosen conditions.",
    )

    columns = [
        (
            0.50,
            "01",
            "Observed profiles",
            "Measured gene expression from real mouse tissue",
            BLUE,
            PALE_BLUE,
        ),
        (
            4.62,
            "02",
            "Learn the distribution",
            "Capture gene relationships and variation by tissue and FLT/GC condition",
            ORANGE,
            PALE_GOLD,
        ),
        (
            8.74,
            "03",
            "Synthetic profiles",
            "Sample new numeric expression vectors for a selected context",
            TEAL,
            PALE_TEAL,
        ),
    ]
    for index, (x, number, heading, body, color, fill) in enumerate(columns):
        if index:
            _add_arrow(slide, x - 0.47, 3.07, 0.34, 0.24, MID_GRAY)
        _add_text(slide, number, x, 2.10, 0.42, 0.25, size=11, color=color, bold=True, margin=0)
        _add_text(slide, heading, x + 0.50, 2.04, 3.08, 0.38, size=17.2, color=NAVY, bold=True, margin=0)
        _add_panel(slide, x, 2.58, 3.58, 1.26, fill=fill, line=fill, radius=False)

        if index in (0, 2):
            cell_colors = [color, "AFC3D4", color, "D8E2E8", color]
            for row in range(4):
                for col in range(5):
                    cell_color = cell_colors[(row * 2 + col + index) % len(cell_colors)]
                    square = slide.shapes.add_shape(
                        MSO_SHAPE.RECTANGLE,
                        Inches(x + 0.22 + col * 0.29),
                        Inches(2.82 + row * 0.20),
                        Inches(0.20),
                        Inches(0.12),
                    )
                    _set_fill(square, cell_color)
                    square.line.fill.background()
            _add_text(slide, "genes", x + 0.31, 3.61, 1.12, 0.16, size=8.2, color=GRAY, align=PP_ALIGN.CENTER, margin=0)
            _add_text(slide, "profiles", x + 1.73, 3.00, 0.62, 0.18, size=8.2, color=GRAY, align=PP_ALIGN.CENTER, margin=0)
        else:
            node_positions = [
                (x + 0.38, 2.85),
                (x + 1.03, 2.75),
                (x + 1.03, 3.25),
                (x + 1.72, 3.00),
                (x + 2.40, 2.75),
                (x + 2.40, 3.25),
                (x + 3.02, 3.00),
            ]
            connections = [(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (3, 5), (4, 6), (5, 6)]
            for start, end in connections:
                x1, y1 = node_positions[start]
                x2, y2 = node_positions[end]
                line = slide.shapes.add_connector(
                    MSO_CONNECTOR.STRAIGHT,
                    Inches(x1 + 0.08),
                    Inches(y1 + 0.08),
                    Inches(x2 + 0.08),
                    Inches(y2 + 0.08),
                )
                _set_line(line, "C6B07A", 1.1)
            for node_index, (node_x, node_y) in enumerate(node_positions):
                _add_circle(slide, node_x, node_y, 0.16, ORANGE if node_index in (0, 3, 6) else GOLD)

        _add_text(slide, body, x + 0.02, 4.03, 3.50, 0.66, size=13.2, color=DARK, margin=0)

    _add_rule(slide, 0.53, 4.95, 12.20, "D6DEE3", 0.018)
    _add_text(slide, "It can help", 0.55, 5.20, 1.60, 0.30, size=16, color=TEAL, bold=True, margin=0)
    _add_text(
        slide,
        "test classifiers, rank candidate genes and compare matched FLT and GC profiles",
        2.03,
        5.16,
        4.17,
        0.55,
        size=14,
        color=DARK,
        margin=0,
    )
    _add_rule(slide, 6.49, 5.18, 0.012, "D6DEE3", 0.70)
    _add_text(slide, "It does not add", 6.82, 5.20, 1.85, 0.30, size=16, color=CORAL, bold=True, margin=0)
    _add_text(
        slide,
        "new animals, new measurements or independent biological evidence",
        8.51,
        5.16,
        3.93,
        0.55,
        size=14,
        color=DARK,
        margin=0,
    )


def _slide_scientific_objective(slide):
    _add_slide_title(
        slide,
        "Study goals",
        "Match tissue distributions, then test FLT versus GC biology",
        "The first goal validates the generator. The second asks whether synthetic data improves a real-data analysis.",
    )

    _add_text(slide, "01", 0.58, 2.08, 0.42, 0.26, size=11, color=TEAL, bold=True, margin=0)
    _add_text(slide, "Reproduce tissue structure", 1.08, 2.02, 4.50, 0.38, size=18, color=NAVY, bold=True, margin=0)
    _add_text(
        slide,
        "Synthetic bulk RNA-seq should occupy the same tissue-defined expression space as real profiles.",
        0.58,
        2.49,
        5.58,
        0.56,
        size=13.4,
        color=DARK,
        margin=0,
    )

    _add_rule(slide, 0.83, 4.58, 4.92, MID_GRAY, 0.018)
    _add_rule(slide, 0.83, 3.14, 0.018, MID_GRAY, 1.46)

    def add_cross(x, y, color):
        first = _add_rule(slide, x, y + 0.06, 0.18, color, 0.025)
        first.rotation = 45
        second = _add_rule(slide, x, y + 0.06, 0.18, color, 0.025)
        second.rotation = -45

    clusters = [
        (1.52, 3.98, TEAL, "Liver"),
        (3.12, 3.43, ORANGE, "Muscle"),
        (4.77, 4.03, PURPLE, "Thymus"),
    ]
    offsets = [(-0.23, -0.10), (0.05, -0.19), (0.22, 0.02), (-0.08, 0.18)]
    for center_x, center_y, color, label in clusters:
        for point_index, (dx, dy) in enumerate(offsets):
            _add_circle(slide, center_x + dx, center_y + dy, 0.13, color)
            add_cross(center_x + dx + 0.10, center_y + dy + (0.08 if point_index % 2 else -0.03), color)
        _add_text(slide, label, center_x - 0.37, center_y + 0.40, 0.98, 0.22, size=9.2, color=GRAY, align=PP_ALIGN.CENTER, margin=0)
    _add_circle(slide, 0.93, 3.05, 0.12, NAVY)
    _add_text(slide, "Real", 1.10, 3.02, 0.55, 0.20, size=9.2, color=GRAY, margin=0)
    add_cross(1.70, 3.03, NAVY)
    _add_text(slide, "Generated", 1.95, 3.02, 0.90, 0.20, size=9.2, color=GRAY, margin=0)

    _add_rule(slide, 6.50, 2.05, 0.015, "D6DEE3", 2.72)

    _add_text(slide, "02", 6.86, 2.08, 0.42, 0.26, size=11, color=CORAL, bold=True, margin=0)
    _add_text(slide, "Preserve the FLT/GC difference", 7.36, 2.02, 4.96, 0.38, size=18, color=NAVY, bold=True, margin=0)
    _add_text(
        slide,
        "Within each tissue, generated profiles should retain the smaller condition signal.",
        6.86,
        2.49,
        5.60,
        0.56,
        size=13.4,
        color=DARK,
        margin=0,
    )
    _add_text(slide, "FLT", 7.18, 3.24, 0.72, 0.28, size=14, color=CORAL, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "GC", 9.08, 3.24, 0.72, 0.28, size=14, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    flt_points = [(7.18, 3.73), (7.52, 3.58), (7.79, 3.81), (7.37, 4.04)]
    gc_points = [(9.09, 3.73), (9.43, 3.58), (9.70, 3.81), (9.28, 4.04)]
    for x, y in flt_points:
        _add_circle(slide, x, y, 0.17, CORAL)
    for x, y in gc_points:
        _add_circle(slide, x, y, 0.17, BLUE)
    _add_text(slide, "vs", 8.38, 3.75, 0.42, 0.24, size=11.5, color=MID_GRAY, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_arrow(slide, 9.98, 3.67, 0.31, 0.28, MID_GRAY)
    _add_panel(slide, 10.49, 3.30, 2.00, 1.14, fill=PALE_BLUE, line="D4E2EC", radius=False)
    _add_text(slide, "OUTPUT", 10.68, 3.47, 1.63, 0.18, size=9.2, color=BLUE, bold=True, margin=0)
    _add_text(
        slide,
        "Tissue-specific\ngenes and pathways",
        10.68,
        3.75,
        1.63,
        0.54,
        size=12.2,
        color=NAVY,
        bold=True,
        margin=0,
    )

    _add_rule(slide, 0.60, 5.12, 12.08, TEAL, 0.035)
    _add_text(slide, "Model choice", 0.62, 5.43, 1.55, 0.30, size=16.5, color=TEAL, bold=True, margin=0)
    _add_text(
        slide,
        "Use the pipeline that matches tissue structure and performs best on held-out real FLT/GC samples.",
        2.18,
        5.34,
        7.40,
        0.58,
        size=15,
        color=NAVY,
        valign=MSO_ANCHOR.MIDDLE,
        margin=0,
    )
    _add_text(slide, "Gene effects and BH FDR use observed OSDR data.", 9.88, 5.42, 2.57, 0.46, size=11.5, color=GRAY, margin=0)


def _slide_2(slide):
    _add_slide_title(
        slide,
        "Question",
        "Small studies and study effects complicate tissue comparisons",
        "Can a generator help rank spaceflight signal without pretending it creates new animals?",
    )
    _add_text(slide, "NASA OSDR", 0.62, 2.04, 2.20, 0.30, size=17, color=BLUE, bold=True)
    _add_text(slide, "Observed spaceflight cohort", 0.62, 2.37, 2.80, 0.24, size=11.5, color=GRAY)
    _add_text(slide, "1,610", 0.58, 2.64, 2.25, 0.69, size=42, color=NAVY, bold=True)
    _add_text(slide, "profiles", 0.64, 3.29, 1.15, 0.25, size=12.5, color=GRAY)
    _add_text(slide, "75", 3.01, 2.72, 1.20, 0.54, size=31, color=NAVY, bold=True)
    _add_text(slide, "accessions", 3.04, 3.26, 1.40, 0.25, size=12.5, color=GRAY)

    bar_x = 0.64
    bar_y = 4.03
    bar_w = 5.24
    flt_w = bar_w * 835 / 1610
    _add_text(slide, "835 FLT", bar_x, 3.65, 1.10, 0.26, size=13.5, color=CORAL, bold=True)
    _add_text(slide, "775 GC", bar_x + bar_w - 1.10, 3.65, 1.10, 0.26, size=13.5, color=BLUE, bold=True, align=PP_ALIGN.RIGHT)
    _add_rule(slide, bar_x, bar_y, flt_w, CORAL, 0.12)
    _add_rule(slide, bar_x + flt_w, bar_y, bar_w - flt_w, BLUE, 0.12)

    _add_rule(slide, 6.49, 2.02, 0.018, "CBD4DA", 2.21)
    _add_text(slide, "ARCHS4 mouse", 6.88, 2.04, 2.40, 0.30, size=17, color=TEAL, bold=True)
    _add_text(slide, "Reference pretraining cohort", 6.88, 2.37, 2.80, 0.24, size=11.5, color=GRAY)
    _add_text(slide, "997,515", 6.84, 2.64, 2.55, 0.69, size=39, color=NAVY, bold=True)
    _add_text(slide, "profiles screened", 6.90, 3.29, 1.75, 0.25, size=12.5, color=GRAY)
    _add_arrow(slide, 9.26, 2.88, 0.42, 0.30, MID_GRAY)
    _add_text(slide, "17,244", 9.86, 2.66, 1.70, 0.54, size=29, color=NAVY, bold=True)
    _add_text(slide, "selected", 9.89, 3.24, 1.10, 0.25, size=12.5, color=GRAY)
    _add_text(slide, "20 tissues", 11.42, 2.73, 1.25, 0.36, size=17, color=TEAL, bold=True)

    _add_rule(slide, 0.62, 4.57, 12.08, "D6DEE3", 0.018)
    _add_text(slide, "Study effects can look like spaceflight biology.", 0.62, 4.82, 11.80, 0.42, size=22, color=NAVY, bold=True)
    challenge_rows = [
        ("01", "Preserve tissue and FLT/GC structure."),
        ("02", "Avoid memorizing individual profiles."),
        ("03", "Test biological effects in observed OSDR samples."),
    ]
    for index, (number, text) in enumerate(challenge_rows):
        y = 5.42 + index * 0.42
        _add_text(slide, number, 0.64, y, 0.42, 0.26, size=10.5, color=ORANGE, bold=True, margin=0)
        _add_text(slide, text, 1.17, y - 0.03, 10.95, 0.31, size=14.5, color=DARK, margin=0)
    _add_source(slide, "Sources: NASA OSDR Biological Data API; ARCHS4 mouse v2.5.")


def _slide_3(slide):
    _add_slide_title(
        slide,
        "Method",
        "We built a configurable bulk RNA-seq generation pipeline",
        "The downstream analysis used the outlined options; the remaining choices stay configurable.",
    )

    def option_tile(label, x, y, w, color, selected=False, *, size=9.1):
        if selected:
            tile = slide.shapes.add_shape(
                MSO_SHAPE.RECTANGLE,
                Inches(x),
                Inches(y),
                Inches(w),
                Inches(0.36),
            )
            _set_fill(tile, WHITE)
            _set_line(tile, color, 1.8)
            _add_rule(slide, x, y, 0.035, color, 0.36)
        _add_text(
            slide,
            label,
            x + 0.05,
            y + 0.055,
            w - 0.10,
            0.23,
            size=size,
            color=NAVY if selected else MID_GRAY,
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
        _add_text(slide, number, x + 0.14, 2.15, 0.34, 0.24, size=11.5, color=MID_GRAY, bold=True, margin=0)
        _add_text(slide, title, x + 0.56, 2.11, 1.64, 0.32, size=14.2, color=NAVY, bold=True, valign=MSO_ANCHOR.MIDDLE)
        _add_rule(slide, x + 0.14, 2.55, 2.05, color, 0.022)

    xs = [0.35, 2.90, 5.45, 8.00, 10.55]
    for separator_x in [2.79, 5.34, 7.89, 10.44]:
        _add_rule(slide, separator_x, 2.03, 0.012, "D9E0E4", 3.24)

    stage_panel(xs[0], "01", "Data scope", BLUE, WHITE)
    segmented_row(xs[0], 2.72, "Sources", [("OSDR API", True), ("ARCHS4", True)], BLUE)
    segmented_row(xs[0], 3.48, "Studies", [("Single", False), ("Multiple", True)], BLUE)
    segmented_row(xs[0], 4.24, "Tissues", [("Per tissue", False), ("All tissues", True)], BLUE)

    stage_panel(xs[1], "02", "Transform", TEAL, WHITE)
    segmented_row(xs[1], 2.72, "Expression", [("Raw", False), ("CPM", False), ("TPM", True)], TEAL)
    segmented_row(xs[1], 3.48, "Scaling", [("None", False), ("Z-score", False), ("MaxAbs", True)], TEAL)
    segmented_row(xs[1], 4.24, "Features", [("All", False), ("HVG", False), ("L1000", True)], TEAL)

    stage_panel(xs[2], "03", "Harmonization", TEAL, WHITE)
    _add_text(slide, "Global or study correction", xs[2] + 0.16, 2.72, 2.02, 0.18, size=8.4, color=GRAY, bold=True, margin=0)
    for index, (label, selected) in enumerate([
        ("None", True),
        ("Within-study z-score", False),
        ("ComBat / MBatch", False),
        ("MOBER", False),
    ]):
        option_tile(label, xs[2] + 0.16, 2.94 + index * 0.50, 2.03, TEAL, selected, size=9.2)

    stage_panel(xs[3], "04", "Model training", ORANGE, WHITE)
    segmented_row(xs[3], 2.72, "Generator", [("WGAN-GP", False), ("DDIM", True)], ORANGE)
    _add_text(slide, "Training source", xs[3] + 0.16, 3.48, 2.02, 0.18, size=8.4, color=GRAY, bold=True, margin=0)
    for index, (label, selected) in enumerate([
        ("OSDR only", False),
        ("ARCHS4 only", False),
        ("Pretrain + adapt", True),
    ]):
        option_tile(label, xs[3] + 0.16, 3.70 + index * 0.50, 2.03, ORANGE, selected, size=9.2)

    stage_panel(xs[4], "05", "Conditioning", TEAL, WHITE)
    _add_text(slide, "Model inputs", xs[4] + 0.16, 2.72, 2.02, 0.18, size=8.4, color=GRAY, bold=True, margin=0)
    for index, (label, selected) in enumerate([
        ("Tissue", True),
        ("FLT / GC", True),
        ("Accession", True),
        ("Material type", True),
        ("Sex / age", False),
    ]):
        option_tile(label, xs[4] + 0.16, 2.92 + index * 0.42, 2.03, TEAL, selected, size=8.9)

    _add_panel(slide, 0.35, 5.57, 12.55, 1.23, fill="F7F9FA", line="F7F9FA", radius=False)
    _add_rule(slide, 0.35, 5.57, 12.55, "C9D3D9", 0.018)
    _add_text(slide, "Downstream", 0.64, 5.80, 1.42, 0.25, size=10.8, color=NAVY, bold=True)
    _add_text(slide, "configuration", 0.64, 6.17, 1.35, 0.24, size=9.6, color=GRAY)

    def selected_step(x, width, heading, value, color):
        _add_panel(slide, x, 5.72, width, 0.49, fill=WHITE, line=color, radius=False)
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

    _add_text(slide, "CONDITION", 2.10, 6.37, 0.72, 0.16, size=7.3, color=TEAL, bold=True, margin=0)
    condition_tiles = [
        (2.88, 0.90, "Tissue"),
        (3.85, 1.00, "FLT / GC"),
        (4.92, 1.14, "Accession"),
        (6.13, 1.38, "Material type"),
    ]
    for x, width, label in condition_tiles:
        option_tile(label, x, 6.27, width, TEAL, True, size=7.8)

    _add_text(slide, "HARMONIZE", 7.76, 6.37, 0.82, 0.16, size=7.3, color=TEAL, bold=True, margin=0)
    option_tile("None", 8.65, 6.27, 0.78, TEAL, True, size=7.8)
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
    _add_panel(slide, 0.48, 6.46, 12.37, 0.48, fill=PALE_BLUE, line="C9DCE9", radius=False)
    _add_text(
        slide,
        "Tissues separate clearly here. Flight and ground-control samples overlap much more, so we test that difference statistically.",
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


def _slide_4(slide, architecture_figure: Path):
    _add_slide_title(
        slide,
        "Validation",
        "DDIM had lower separability and distributional distance",
        "Both generators matched expression. DDIM was harder to distinguish from real profiles and had higher F1.",
    )
    _add_text(slide, "Architecture from Lacan et al.", 0.52, 2.08, 7.12, 0.32, size=16.5, color=NAVY, bold=True)
    _add_picture_contain(
        slide,
        architecture_figure,
        0.20,
        2.39,
        7.48,
        2.73,
        alt="Lacan et al. Figure 1C residual diffusion generator architecture",
    )
    _add_text(slide, "Residual dense denoiser conditioned on timestep and tissue.", 0.56, 5.17, 6.90, 0.28, size=11.3, color=DARK, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "OSDR adaptation adds FLT/GC, accession and material context through LoRA.", 0.56, 5.51, 6.90, 0.34, size=10.6, color=GRAY, align=PP_ALIGN.CENTER)

    _add_rule(slide, 7.83, 2.05, 0.015, "D5DDE2", 4.32)
    _add_text(slide, "Metric", 8.04, 2.17, 1.42, 0.28, size=10.0, color=GRAY, bold=True)
    _add_text(slide, "WGAN-GP", 9.57, 2.17, 0.85, 0.28, size=10.0, color=CORAL, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "DDIM", 10.52, 2.17, 0.76, 0.28, size=10.0, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
    _add_text(slide, "Read", 11.43, 2.17, 1.10, 0.28, size=10.0, color=GRAY, bold=True)
    _add_rule(slide, 8.04, 2.56, 4.58, NAVY, 0.022)
    metrics = [
        ("Correlation", "0.976", "0.974", "similar"),
        ("F1", "0.985", "0.997", "higher"),
        ("Adversarial acc.", "0.636", "0.475", "near 0.5"),
        ("FD / real P95", "0.144", "0.074", "lower"),
    ]
    for index, (label, wgan, ddim, reading) in enumerate(metrics):
        y = 2.80 + index * 0.59
        if index % 2 == 0:
            shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(8.02), Inches(y - 0.07), Inches(4.61), Inches(0.50))
            _set_fill(shape, "F5F7F8")
            shape.line.fill.background()
        _add_text(slide, label, 8.10, y, 1.39, 0.32, size=10.3, color=DARK, bold=index >= 2)
        _add_text(slide, wgan, 9.57, y, 0.85, 0.32, size=12.6, color=CORAL, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, ddim, 10.52, y, 0.76, 0.32, size=12.6, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, reading, 11.43, y, 1.10, 0.32, size=10.2, color=GRAY)

    _add_rule(slide, 0.55, 5.91, 7.03, "D5DDE2", 0.015)
    _add_text(slide, "ARCHS4 tissue probe", 0.55, 6.08, 1.75, 0.24, size=10.4, color=BLUE, bold=True, margin=0)
    _add_text(slide, "Real 0.781  |  Synthetic 0.781", 4.30, 6.06, 3.12, 0.27, size=10.4, color=DARK, align=PP_ALIGN.RIGHT, margin=0)
    _add_text(slide, "Balanced accuracy is preserved in generated tissue labels.", 0.55, 6.39, 6.87, 0.25, size=9.7, color=GRAY, align=PP_ALIGN.CENTER, margin=0)
    _add_rule(slide, 8.04, 5.36, 0.055, TEAL, 0.78)
    _add_text(slide, "Use DDIM", 8.29, 5.40, 1.32, 0.31, size=16, color=NAVY, bold=True, margin=0)
    _add_text(slide, "Higher F1 with lower separability and FD.", 9.73, 5.38, 2.78, 0.47, size=11.3, color=DARK, valign=MSO_ANCHOR.MIDDLE, margin=0)
    _add_text(slide, "AA near 0.5 means real and generated profiles are difficult to separate.", 8.29, 5.97, 4.22, 0.40, size=9.5, color=GRAY, margin=0)
    _add_source(slide, "Architecture excerpt: Lacan et al. (2026), Fig. 1C, doi:10.1186/s12859-026-06470-8. Metrics use each model's stated split and are not paired.")


def _slide_7(slide):
    _add_slide_title(
        slide,
        "Analysis",
        "Five arms separate gene ranking from classifier fitting",
        "Guided arms use synthetic profiles to help choose genes; they do not treat generated profiles as new animals.",
    )
    _add_text(slide, "INPUTS", 0.47, 2.00, 0.66, 0.22, size=8.8, color=GRAY, bold=True, margin=0)
    _add_data_badge(slide, "R", 1.20, 1.98, BLUE)
    _add_text(slide, "real OSDR", 1.52, 1.98, 1.18, 0.24, size=10.8, color=BLUE, bold=True, margin=0)
    _add_data_badge(slide, "S", 2.91, 1.98, CORAL)
    _add_text(slide, "matched DDIM", 3.23, 1.98, 1.42, 0.24, size=10.8, color=CORAL, bold=True, margin=0)
    _add_text(slide, "GUIDED", 6.10, 1.99, 0.70, 0.22, size=8.8, color=TEAL, bold=True, margin=0)
    _add_text(slide, "R + S help rank genes", 6.85, 1.97, 2.05, 0.25, size=10.5, color=DARK, margin=0)
    _add_text(slide, "5%", 9.41, 1.99, 0.38, 0.22, size=8.8, color=ORANGE, bold=True, margin=0)
    _add_text(slide, "S carries 5% of total fit weight", 9.84, 1.97, 2.78, 0.25, size=10.5, color=DARK, margin=0)

    headers = [
        ("ANALYSIS ARM", 0.66, 2.18),
        ("GENE RANKING", 3.07, 2.90),
        ("CLASSIFIER FIT", 6.43, 2.95),
        ("WHAT IT TESTS", 9.82, 2.72),
    ]
    for heading, x, width in headers:
        _add_text(slide, heading, x, 2.38, width, 0.22, size=9.0, color=GRAY, bold=True, margin=0)
    _add_rule(slide, 0.47, 2.64, 12.18, NAVY, 0.022)

    arms = [
        ("Real only", BLUE, ("R",), "Rank from real", ("R",), "Fit real", "Baseline"),
        ("Generated only", CORAL, ("S",), "Rank from generated", ("S",), "Fit generated", "Synthetic-only stress test"),
        ("Real + generated", TEAL, ("R", "S"), "Consensus rank", ("R", "S"), "Equal total weight", "Direct augmentation"),
        ("Guided: real fit", TEAL, ("R", "S"), "Consensus rank", ("R",), "Fit real only", "Feature guidance only"),
        ("Guided: 5% synthetic", ORANGE, ("R", "S"), "Consensus rank", ("R", "S"), "S = 5% weight", "Guidance + light regularization"),
    ]

    def add_sources(tokens: tuple[str, ...], x: float, y: float) -> float:
        if len(tokens) == 1:
            color = BLUE if tokens[0] == "R" else CORAL
            _add_data_badge(slide, tokens[0], x, y, color)
            return x + 0.34
        _add_data_badge(slide, "R", x, y, BLUE)
        _add_text(slide, "+", x + 0.27, y + 0.01, 0.18, 0.20, size=9.5, color=GRAY, bold=True, align=PP_ALIGN.CENTER, margin=0)
        _add_data_badge(slide, "S", x + 0.47, y, CORAL)
        return x + 0.82

    for index, (arm, color, rank_tokens, rank_label, fit_tokens, fit_label, purpose) in enumerate(arms):
        top = 2.73 + index * 0.57
        if index % 2 == 0:
            shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.47), Inches(top), Inches(12.18), Inches(0.50))
            _set_fill(shape, "F4F6F7")
            shape.line.fill.background()
        _add_rule(slide, 0.47, top, 0.055, color, 0.50)
        _add_text(slide, arm, 0.67, top + 0.08, 2.17, 0.33, size=11.2, color=color, bold=True, valign=MSO_ANCHOR.MIDDLE, margin=0)

        rank_end = add_sources(rank_tokens, 3.08, top + 0.13)
        _add_arrow(slide, rank_end, top + 0.17, 0.26, 0.14, MID_GRAY)
        _add_text(slide, rank_label, rank_end + 0.37, top + 0.08, 1.78, 0.33, size=10.4, color=DARK, bold=True, valign=MSO_ANCHOR.MIDDLE, margin=0)

        fit_end = add_sources(fit_tokens, 6.45, top + 0.13)
        _add_arrow(slide, fit_end, top + 0.17, 0.26, 0.14, MID_GRAY)
        _add_text(slide, fit_label, fit_end + 0.37, top + 0.08, 1.71, 0.33, size=10.4, color=DARK, bold=True, valign=MSO_ANCHOR.MIDDLE, margin=0)
        _add_text(slide, purpose, 9.83, top + 0.08, 2.70, 0.33, size=10.5, color=DARK, valign=MSO_ANCHOR.MIDDLE, margin=0)

    _add_rule(slide, 0.47, 5.64, 12.18, "CDD6DC", 0.018)
    _add_text(slide, "All five arms", 0.62, 5.83, 1.36, 0.27, size=12.5, color=NAVY, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_arrow(slide, 2.13, 5.87, 0.34, 0.19, MID_GRAY)
    _add_text(slide, "Held-out real profiles", 2.66, 5.82, 2.00, 0.30, size=12.5, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_arrow(slide, 4.83, 5.87, 0.34, 0.19, MID_GRAY)
    _add_text(slide, "BA  |  AUROC  |  AP", 5.37, 5.83, 1.88, 0.27, size=12.2, color=DARK, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_arrow(slide, 7.41, 5.87, 0.34, 0.19, MID_GRAY)
    _add_text(slide, "Choose an eligible arm within each tissue", 7.92, 5.80, 4.40, 0.34, size=12.8, color=TEAL, bold=True, align=PP_ALIGN.CENTER, margin=0)

    _add_panel(slide, 0.46, 6.25, 12.22, 0.58, fill=NAVY, line=NAVY, radius=False)
    _add_text(slide, "Association test", 0.74, 6.40, 1.56, 0.25, size=13.5, color=WHITE, bold=True, margin=0)
    _add_text(slide, "FLT vs GC effects and BH FDR use real OSDR profiles only.", 2.47, 6.39, 6.40, 0.27, size=12.9, color="DCE7F2", margin=0)
    _add_text(slide, "Animal n stays unchanged.", 9.19, 6.39, 3.13, 0.27, size=11.5, color="FFD69A", bold=True, align=PP_ALIGN.RIGHT, margin=0)


def _slide_guidance_mechanism(slide):
    _add_slide_title(
        slide,
        "Analysis",
        "Consensus ranking chooses the classifier input genes",
        "Real and generated rankings combine, then only the top k genes are used for FLT-GC classification.",
    )

    def add_rank_track(x: float, y: float, label: str, position: float, color: str):
        track_x = x + 0.90
        track_w = 1.92
        _add_text(slide, label, x, y - 0.055, 0.78, 0.23, size=10.4, color=DARK, bold=True, margin=0)
        _add_rule(slide, track_x, y + 0.035, track_w, "C9D2D8", 0.025)
        _add_circle(slide, track_x + track_w * position - 0.075, y - 0.025, 0.15, color)

    panel_y = 2.05
    panel_h = 2.60
    _add_panel(slide, 0.48, panel_y, 3.61, panel_h, fill=PALE_BLUE, line="C9DCE9", radius=False)
    _add_text(slide, "1  REAL RANKING", 0.72, 2.27, 1.86, 0.25, size=14.0, color=BLUE, bold=True, margin=0)
    _add_text(slide, "Rank all 974 genes using real OSDR", 0.72, 2.61, 3.02, 0.23, size=10.4, color=DARK, margin=0)
    _add_text(slide, "lower-ranked", 1.59, 2.96, 0.82, 0.18, size=8.5, color=GRAY, margin=0)
    _add_text(slide, "top-ranked", 3.02, 2.96, 0.72, 0.18, size=8.5, color=GRAY, align=PP_ALIGN.RIGHT, margin=0)
    example_genes = ["Gene 1", "Gene 2", "Gene 3", "Gene 4", "Gene 5"]
    real_positions = [0.90, 0.82, 0.72, 0.56, 0.35]
    for index, (gene, position) in enumerate(zip(example_genes, real_positions)):
        add_rank_track(0.72, 3.18 + index * 0.25, gene, position, BLUE)

    _add_arrow(slide, 4.16, 2.93, 0.27, 0.20, MID_GRAY)
    _add_panel(slide, 4.50, panel_y, 3.61, panel_h, fill=PALE_CORAL, line="E6CEC8", radius=False)
    _add_text(slide, "2  GENERATED RANKING", 4.74, 2.27, 2.48, 0.25, size=14.0, color=CORAL, bold=True, margin=0)
    _add_text(slide, "Rank the same 974 genes using generated profiles", 4.74, 2.61, 3.07, 0.23, size=10.4, color=DARK, margin=0)
    _add_text(slide, "lower-ranked", 5.61, 2.96, 0.82, 0.18, size=8.5, color=GRAY, margin=0)
    _add_text(slide, "top-ranked", 7.04, 2.96, 0.72, 0.18, size=8.5, color=GRAY, align=PP_ALIGN.RIGHT, margin=0)
    generated_positions = [0.88, 0.45, 0.86, 0.68, 0.32]
    for index, (gene, position) in enumerate(zip(example_genes, generated_positions)):
        add_rank_track(4.74, 3.18 + index * 0.25, gene, position, CORAL)
    _add_arrow(slide, 8.18, 2.93, 0.27, 0.20, MID_GRAY)
    _add_panel(slide, 8.52, panel_y, 4.13, panel_h, fill=PALE_TEAL, line="C9DFDB", radius=False)
    _add_text(slide, "3  CONSENSUS RANKING", 8.76, 2.27, 2.55, 0.25, size=14.0, color=TEAL, bold=True, margin=0)
    _add_text(slide, "Genes supported by both move upward", 8.76, 2.61, 2.82, 0.23, size=10.4, color=DARK, margin=0)
    _add_text(slide, "lower-ranked", 9.66, 2.96, 0.82, 0.18, size=8.5, color=GRAY, margin=0)
    _add_text(slide, "top-ranked", 11.61, 2.96, 0.72, 0.18, size=8.5, color=GRAY, align=PP_ALIGN.RIGHT, margin=0)
    shared_rows = [
        ("Gene 1", 0.89, TEAL, "selected"),
        ("Gene 2", 0.61, MID_GRAY, "not selected"),
        ("Gene 3", 0.79, TEAL, "selected"),
        ("Gene 4", 0.62, MID_GRAY, "not selected"),
        ("Gene 5", 0.33, MID_GRAY, "not selected"),
    ]
    for index, (gene, position, color, outcome) in enumerate(shared_rows):
        y = 3.18 + index * 0.25
        add_rank_track(8.76, y, gene, position, color)
        _add_text(slide, outcome, 11.72, y - 0.055, 0.80, 0.21, size=7.8, color=color, bold=True, align=PP_ALIGN.RIGHT, margin=0)
    shared_cutoff_x = 8.76 + 0.90 + 1.92 * 0.75
    _add_rule(slide, shared_cutoff_x, 3.10, 0.022, TEAL, 1.27)
    _add_text(slide, "selection cutoff", 11.15, 4.42, 1.20, 0.18, size=8.5, color=TEAL, italic=True, align=PP_ALIGN.RIGHT, margin=0)

    _add_text(slide, "Top-ranked genes become classifier inputs", 0.52, 4.84, 4.65, 0.28, size=16.0, color=NAVY, bold=True, margin=0)
    _add_panel(slide, 0.51, 5.18, 3.37, 1.28, fill="F4F8FA", line="D5E1E7", radius=False)
    _add_text(slide, "Rank 974 genes", 0.77, 5.45, 2.84, 0.25, size=13.0, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "real, generated, or consensus", 0.77, 5.88, 2.84, 0.22, size=10.0, color=GRAY, align=PP_ALIGN.CENTER, margin=0)
    _add_arrow(slide, 3.99, 5.68, 0.31, 0.20, MID_GRAY)

    _add_panel(slide, 4.41, 5.18, 3.37, 1.28, fill=PALE_TEAL, line="C9DFDB", radius=False)
    _add_text(slide, "Keep the top k", 4.67, 5.39, 2.84, 0.25, size=13.0, color=TEAL, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "k = 10, 25, 50, or 100", 4.67, 5.76, 2.84, 0.22, size=10.2, color=DARK, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "validation chooses k", 4.67, 6.06, 2.84, 0.18, size=8.8, color=GRAY, italic=True, align=PP_ALIGN.CENTER, margin=0)
    _add_arrow(slide, 7.89, 5.68, 0.31, 0.20, MID_GRAY)

    _add_panel(slide, 8.31, 5.18, 4.34, 1.28, fill="F6F7F8", line="D9DFE3", radius=False)
    _add_text(slide, "Fit FLT-GC classifier", 8.57, 5.45, 3.82, 0.25, size=13.0, color=NAVY, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "using only the selected gene columns", 8.57, 5.88, 3.82, 0.22, size=10.0, color=DARK, align=PP_ALIGN.CENTER, margin=0)

    _add_source(slide, "Five genes are schematic examples from the 974-gene ranking. Real and generated normalized ranks are combined.")


def _slide_guidance_boundary(slide):
    _add_slide_title(
        slide,
        "Analysis",
        "Held-out real samples decide whether guidance helps",
        "Opaque profiles train the classifier; transparent profiles are used only to score generalization.",
    )

    def add_point(x: float, y: float, color: str, *, held_out: bool):
        if not held_out:
            _add_circle(slide, x, y, 0.17, color)
            return
        shape = slide.shapes.add_shape(
            MSO_SHAPE.OVAL,
            Inches(x),
            Inches(y),
            Inches(0.20),
            Inches(0.20),
        )
        _set_fill(shape, color, transparency=58)
        _set_line(shape, color, 1.2)

    def add_scatter(x: float, y: float, *, guided: bool):
        plot_w = 4.35
        plot_h = 1.95
        _add_rule(slide, x, y + plot_h, plot_w, MID_GRAY, 0.022)
        _add_rule(slide, x, y, 0.022, MID_GRAY, plot_h + 0.02)
        train_flight = [(0.42, 0.34), (0.75, 1.27), (1.07, 0.78), (1.37, 1.66), (1.58, 0.43)]
        train_ground = [(2.58, 0.36), (2.88, 1.28), (3.19, 0.76), (3.55, 1.65), (3.88, 0.46)]
        held_flight = [(1.84, 0.86), (2.03, 1.42)]
        held_ground = [(2.29, 0.59), (2.48, 1.55)]
        for px, py in train_flight:
            add_point(x + px, y + py, CORAL, held_out=False)
        for px, py in train_ground:
            add_point(x + px, y + py, BLUE, held_out=False)
        for px, py in held_flight:
            add_point(x + px, y + py, CORAL, held_out=True)
        for px, py in held_ground:
            add_point(x + px, y + py, BLUE, held_out=True)
        if guided:
            boundary = _add_rule(slide, x + 2.15, y + 0.12, 0.045, TEAL, 1.66)
        else:
            boundary = _add_rule(slide, x + 1.70, y + 0.12, 0.035, GRAY, 1.66)
        _add_text(slide, "predicted FLT", x + 0.10, y + 2.04, 1.02, 0.18, size=8.4, color=CORAL, margin=0)
        _add_text(slide, "predicted GC", x + 3.19, y + 2.04, 1.02, 0.18, size=8.4, color=BLUE, align=PP_ALIGN.RIGHT, margin=0)

    _add_circle(slide, 8.79, 1.96, 0.13, CORAL)
    _add_text(slide, "FLT", 8.98, 1.92, 0.36, 0.21, size=9.2, color=GRAY, margin=0)
    _add_circle(slide, 9.50, 1.96, 0.13, BLUE)
    _add_text(slide, "GC", 9.69, 1.92, 0.34, 0.21, size=9.2, color=GRAY, margin=0)
    _add_circle(slide, 10.30, 1.96, 0.13, MID_GRAY)
    _add_text(slide, "training", 10.49, 1.92, 0.62, 0.21, size=9.2, color=GRAY, margin=0)
    held_legend = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(11.28), Inches(1.94), Inches(0.16), Inches(0.16))
    _set_fill(held_legend, MID_GRAY, transparency=58)
    _set_line(held_legend, MID_GRAY, 1.0)
    _add_text(slide, "held-out", 11.52, 1.92, 0.72, 0.21, size=9.2, color=GRAY, margin=0)

    _add_text(slide, "REAL-ONLY CLASSIFIER", 0.78, 2.22, 4.35, 0.27, size=13.0, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "SYNTHETIC-GUIDED CANDIDATE", 7.10, 2.22, 4.35, 0.27, size=13.0, color=TEAL, bold=True, align=PP_ALIGN.CENTER, margin=0)
    add_scatter(0.78, 2.56, guided=False)
    add_scatter(7.10, 2.56, guided=True)
    _add_arrow(slide, 5.53, 3.20, 0.72, 0.40, MID_GRAY)
    _add_text(slide, "same held-out profiles", 5.20, 3.73, 1.38, 0.37, size=9.4, color=GRAY, bold=True, align=PP_ALIGN.CENTER, margin=0)

    _add_panel(slide, 0.57, 5.10, 12.08, 1.25, fill="F5F7F8", line="D5DDE2", radius=False)
    _add_text(slide, "Compare on held-out real profiles", 0.84, 5.31, 2.42, 0.28, size=13.0, color=NAVY, bold=True, margin=0)
    _add_text(slide, "Balanced accuracy", 3.47, 5.30, 1.32, 0.24, size=10.3, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_rule(slide, 4.91, 5.28, 0.015, "D1D9DE", 0.35)
    _add_text(slide, "AUROC", 5.12, 5.30, 0.80, 0.24, size=10.3, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_rule(slide, 6.10, 5.28, 0.015, "D1D9DE", 0.35)
    _add_text(slide, "Average precision", 6.31, 5.30, 1.34, 0.24, size=10.3, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_arrow(slide, 7.94, 5.32, 0.34, 0.20, MID_GRAY)
    _add_text(slide, "Retain the better eligible arm for that tissue", 8.55, 5.25, 3.68, 0.36, size=12.0, color=TEAL, bold=True, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE, margin=0)
    _add_rule(slide, 0.84, 5.77, 11.38, "D9E0E4", 0.015)
    _add_text(slide, "Held-out profiles never enter gene ranking, classifier fitting, or k selection.", 0.84, 5.92, 11.38, 0.24, size=10.3, color=DARK, bold=True, align=PP_ALIGN.CENTER, margin=0)

    _add_source(slide, "Schematic nested evaluation. Observed tissue-specific held-out results follow on the next slide.")


def _slide_8(slide, utility_chart: Path):
    _add_slide_title(
        slide,
        "Utility",
        "Pooled training missed improvements seen within tissues",
        "The pooled FLT/GC classifier declined, while several tissue-specific classifiers improved.",
    )
    _add_rule(slide, 0.48, 2.05, 3.30, GOLD, 0.035)
    _add_text(slide, "Pooled across tissues", 0.67, 2.25, 2.90, 0.35, size=17, color=NAVY, bold=True)
    labels = ["Real", "Generated", "Real + synth."]
    values = [0.754, 0.695, 0.737]
    colors = [BLUE, CORAL, TEAL]
    baseline = 6.00
    chart_height = 3.00
    chart_left = 0.96
    chart_right = 3.62
    axis_title = _add_text(
        slide,
        "Balanced accuracy",
        -0.28,
        4.35,
        1.55,
        0.28,
        size=10.5,
        color=GRAY,
        bold=True,
        align=PP_ALIGN.CENTER,
        margin=0,
    )
    axis_title.rotation = 270
    for tick, label in [(0.0, "0"), (0.5, "0.5"), (1.0, "1.0")]:
        tick_y = baseline - tick * chart_height
        _add_rule(slide, chart_left, tick_y, chart_right - chart_left, "D8DEE2", 0.012)
        _add_text(
            slide,
            label,
            0.69,
            tick_y - 0.10,
            0.22,
            0.20,
            size=8.5,
            color=MID_GRAY,
            align=PP_ALIGN.RIGHT,
            margin=0,
        )
    for index, (label, value, color) in enumerate(zip(labels, values, colors)):
        x = 0.99 + index * 0.89
        bar_h = value * chart_height
        shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(baseline - bar_h), Inches(0.62), Inches(bar_h))
        _set_fill(shape, color)
        shape.line.fill.background()
        _add_text(slide, f"{value:.3f}", x - 0.08, baseline - bar_h - 0.35, 0.78, 0.28, size=14, color=NAVY, bold=True, align=PP_ALIGN.CENTER)
        _add_text(slide, label, x - 0.13, baseline + 0.08, 0.88, 0.32, size=9.0, color=DARK, align=PP_ALIGN.CENTER, margin=0)
    _add_rule(slide, chart_left, baseline, chart_right - chart_left, "86949D", 0.02)

    _add_rule(slide, 3.98, 2.05, 0.015, "D5DDE2", 4.55)
    _add_rule(slide, 4.25, 2.05, 8.48, TEAL, 0.035)
    _add_text(slide, "Selected use within each tissue", 4.42, 2.25, 4.2, 0.34, size=17, color=NAVY, bold=True)
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
        "Blue marks real-only selection; teal marks synthetic-guided selection. Every displayed association has real OSDR support.",
    )
    _add_text(slide, "Stable selection after real-data support", 0.62, 2.02, 4.30, 0.30, size=15.5, color=NAVY, bold=True, margin=0)
    _add_circle(slide, 0.64, 2.42, 0.13, BLUE)
    _add_text(slide, "Blue: stable in real-only ranking", 0.86, 2.36, 2.42, 0.25, size=10.5, color=BLUE, bold=True, margin=0)
    _add_circle(slide, 3.43, 2.42, 0.13, TEAL)
    _add_text(slide, "Teal: stable in synthetic-guided ranking", 3.65, 2.36, 2.98, 0.25, size=10.5, color=TEAL, bold=True, margin=0)

    left = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.62), Inches(2.78), Inches(3.86), Inches(2.78))
    _set_fill(left, BLUE, transparency=82)
    _set_line(left, BLUE, 2.6)
    right = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(3.00), Inches(2.78), Inches(3.86), Inches(2.78))
    _set_fill(right, TEAL, transparency=80)
    _set_line(right, TEAL, 2.6)

    _add_text(slide, "34", 1.16, 3.42, 1.42, 0.48, size=29, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "real-only", 1.12, 3.93, 1.50, 0.30, size=14.5, color=NAVY, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "WITHOUT GUIDANCE", 1.03, 4.31, 1.68, 0.24, size=9.2, color=BLUE, bold=True, align=PP_ALIGN.CENTER, margin=0)

    _add_text(slide, "23", 3.16, 3.42, 1.17, 0.48, size=29, color=ORANGE, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "reinforced", 3.10, 3.93, 1.29, 0.30, size=13.2, color=NAVY, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "BOTH RANKINGS", 3.04, 4.31, 1.40, 0.24, size=9.2, color=ORANGE, bold=True, align=PP_ALIGN.CENTER, margin=0)

    _add_text(slide, "26", 4.89, 3.42, 1.50, 0.48, size=29, color=TEAL, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "promoted", 4.89, 3.93, 1.50, 0.30, size=14.5, color=NAVY, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "GUIDANCE ONLY", 4.80, 4.31, 1.68, 0.24, size=9.2, color=TEAL, bold=True, align=PP_ALIGN.CENTER, margin=0)

    _add_rule(slide, 7.18, 2.10, 0.018, "D5DDE2", 3.98)
    _add_text(slide, "Real OSDR evidence gate", 7.53, 2.43, 4.92, 0.40, size=20, color=NAVY, bold=True, align=PP_ALIGN.CENTER)
    evidence_rows = [
        (3.18, 0.38, "Estimate FLT vs GC inside each accession"),
        (3.82, 0.58, "Combine accession effects with a\nrandom-effects model"),
        (4.63, 0.38, "Apply BH FDR within each tissue"),
    ]
    for y, height, label in evidence_rows:
        _add_circle(slide, 7.79, y + 0.10, 0.10, BLUE)
        _add_text(
            slide,
            label,
            8.05,
            y,
            4.15,
            height,
            size=14.2,
            color=DARK,
            valign=MSO_ANCHOR.TOP,
            margin=0,
        )
    _add_panel(slide, 7.68, 5.14, 4.46, 0.78, fill=NAVY, line=NAVY, radius=False)
    _add_text(slide, "23 reinforced + 26 promoted", 7.90, 5.26, 4.02, 0.24, size=11.5, color="BFD0E1", bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "49 synthetic-informed associations", 7.88, 5.52, 4.06, 0.28, size=15.5, color=WHITE, bold=True, align=PP_ALIGN.CENTER, margin=0)

    _add_rule(slide, 0.52, 6.26, 12.15, "D6DEE3", 0.018)
    _add_text(slide, "The 34 real-only associations remain real-data findings; synthetic guidance did not reinforce them.", 0.72, 6.39, 11.75, 0.26, size=12.6, color=GRAY, align=PP_ALIGN.CENTER)
    _add_source(slide, "Counts are tissue-gene associations after BH FDR and compatible real-effect direction; genes may recur across tissues.")


def _slide_10(slide):
    _add_slide_title(
        slide,
        "Literature interpretation",
        "Selection and literature are separate dimensions",
        "All 49 synthetic-informed associations received both a selection status and a literature interpretation.",
    )

    steps = [
        ("01", "Fix the candidate", "Gene, tissue, and FLT direction", BLUE),
        ("02", "Search", "Spaceflight and mechanism literature", TEAL),
        ("03", "Compare", "Observed result with published evidence", ORANGE),
        ("04", "Record", "Category, rationale, and source IDs", PURPLE),
    ]
    for index, (number, heading, detail, color) in enumerate(steps):
        x = 0.48 + index * 3.16
        if index:
            _add_rule(slide, x - 0.28, 2.02, 0.012, "D8DFE3", 0.95)
        _add_text(slide, number, x, 2.07, 0.40, 0.24, size=10.8, color=color, bold=True, margin=0)
        _add_text(slide, heading, x + 0.48, 2.02, 2.30, 0.30, size=14.5, color=NAVY, bold=True, margin=0)
        _add_text(slide, detail, x + 0.48, 2.40, 2.30, 0.42, size=10.8, color=GRAY, margin=0)
        _add_rule(slide, x, 2.95, 2.72, color, 0.025)

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

    _add_text(slide, "Selection status", 0.55, 3.35, 3.58, 0.34, size=17, color=NAVY, bold=True)
    selection_rows = [
        ("Promoted", selection_counts.get("promoted", 0), "Stable only with synthetic guidance", CORAL),
        ("Reinforced", selection_counts.get("reinforced", 0), "Stable with and without guidance", TEAL),
    ]
    for index, (label, count, detail, color) in enumerate(selection_rows):
        y = 4.02 + index * 1.05
        if index:
            _add_rule(slide, 0.55, y - 0.18, 3.58, "E0E5E8", 0.012)
        _add_text(slide, str(count), 0.55, y - 0.05, 0.68, 0.52, size=29, color=color, bold=True, margin=0)
        _add_text(slide, label, 1.38, y - 0.02, 2.65, 0.30, size=15, color=color, bold=True, margin=0)
        _add_text(slide, detail, 1.38, y + 0.31, 2.70, 0.38, size=11.5, color=DARK, margin=0)

    _add_rule(slide, 4.45, 3.33, 0.015, "D5DDE2", 2.77)
    _add_text(slide, "Literature interpretation", 4.82, 3.35, 7.52, 0.34, size=17, color=NAVY, bold=True)
    _add_text(slide, "Class", 4.84, 3.86, 1.55, 0.22, size=9.5, color=GRAY, bold=True, margin=0)
    _add_text(slide, "Count", 6.86, 3.86, 0.70, 0.22, size=9.5, color=GRAY, bold=True, align=PP_ALIGN.CENTER, margin=0)
    _add_text(slide, "Meaning", 7.69, 3.86, 4.31, 0.22, size=9.5, color=GRAY, bold=True, margin=0)
    _add_rule(slide, 4.82, 4.13, 7.54, NAVY, 0.020)
    category_rows = [
        ("Aligning", counts.get("aligning", 0), "Direct or same-tissue process agreement", GREEN),
        ("Complementary", counts.get("complementary", 0), "Related process or mechanism", BLUE),
        ("Ambiguous", counts.get("ambiguous", 0), "Mixed or context-dependent evidence", ORANGE),
        ("Unmatched", counts.get("unmatched", 0), "No sufficiently specific match found", PURPLE),
    ]
    for index, (label, count, detail, color) in enumerate(category_rows):
        y = 4.27 + index * 0.48
        if index % 2 == 0:
            shade = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(4.82), Inches(y - 0.04), Inches(7.54), Inches(0.42))
            _set_fill(shade, "F5F7F8")
            shade.line.fill.background()
        _add_rule(slide, 4.86, y + 0.11, 0.12, color, 0.12)
        _add_text(slide, label, 5.10, y, 1.68, 0.34, size=12.6, color=color, bold=True, valign=MSO_ANCHOR.MIDDLE, margin=0)
        _add_text(slide, str(count), 6.87, y, 0.66, 0.34, size=14.0, color=DARK, bold=True, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE, margin=0)
        _add_text(slide, detail, 7.69, y, 4.41, 0.34, size=11.8, color=DARK, valign=MSO_ANCHOR.MIDDLE, margin=0)

    _add_rule(slide, 0.55, 6.35, 12.02, "CED7DD", 0.018)
    _add_text(
        slide,
        "A reinforced gene can be complementary, ambiguous, or unmatched; an aligning gene can be promoted.",
        0.65,
        6.49,
        11.82,
        0.28,
        size=13.5,
        color=NAVY,
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
        "The screen covered all 27 completed tissue analyses",
        "The biological discussion narrows only after reporting synthetic-informed, real-only, and null outcomes.",
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
    for index, (x, count, heading, tissue_names, color, fill, detail) in enumerate(columns):
        if index:
            _add_rule(slide, x - 0.18, 2.05, 0.012, "D8DFE3", 4.58)
        _add_rule(slide, x, 2.04, 3.65, color, 0.035)
        _add_text(slide, count, x + 0.02, 2.22, 0.70, 0.55, size=31, color=color, bold=True, margin=0)
        _add_text(slide, heading, x + 0.84, 2.27, 2.72, 0.38, size=16.0, color=NAVY, bold=True, valign=MSO_ANCHOR.MIDDLE, margin=0)
        _add_text(slide, detail, x + 0.02, 2.92, 3.48, 0.42, size=10.5, color=GRAY, valign=MSO_ANCHOR.MIDDLE, margin=0)
        _add_rule(slide, x, 3.50, 3.58, "D7DEE2", 0.014)
        _add_text(slide, "ANALYSIS UNITS", x + 0.02, 3.66, 1.20, 0.19, size=8.8, color=GRAY, bold=True, margin=0)

        list_top = 3.96
        list_height = 2.42
        row_height = min(0.36, list_height / len(tissue_names))
        list_size = 9.8 if len(tissue_names) >= 12 else 10.4 if len(tissue_names) >= 10 else 11.2
        for index, tissue_name in enumerate(tissue_names):
            row_y = list_top + index * row_height
            _add_rule(slide, x + 0.03, row_y + (row_height - 0.045) / 2, 0.13, color, 0.045)
            _add_text(
                slide,
                tissue_name,
                x + 0.25,
                row_y,
                3.25,
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
        "Ten tissue analyses contained synthetic-informed genes",
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
        "The cell-cycle program persisted without material conditioning; two FLT-higher genes extend the hypothesis.",
    )
    figure = PAPER_DIR / "figures/figure_3_thymus_biology.png"
    _add_picture_contain(slide, figure, 0.35, 1.93, 8.55, 4.78, alt="Thymus gene effects and Reactome processes")
    _add_panel(slide, 9.10, 2.02, 3.82, 4.57, fill=PALE_CORAL, line="E8C9C2", radius=False)
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
        "The real-data program is coherent, but its synthetic reinforcement required explicit material conditioning.",
    )
    figure = PAPER_DIR / "figures/figure_4_soleus_biology.png"
    _add_picture_contain(slide, figure, 0.34, 1.96, 9.15, 4.70, alt="Soleus gene effects and Reactome processes")
    _add_panel(slide, 9.63, 2.04, 3.25, 4.52, fill=PALE_TEAL, line="C7DDD8", radius=False)
    _add_text(slide, "0.925 to 0.963", 9.89, 2.34, 2.72, 0.42, size=20, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
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
        "Additional tissues produced distinct hypotheses",
        "Promoted and reinforced selections are shown separately for each tissue.",
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
        "Synthetic data worked best as a tissue-specific prior",
        "The generator helped prioritize observed signal. Biological claims still come from OSDR profiles.",
    )
    columns = [
        (0.42, "1", "Generate", "Conditional DDIM produced high-fidelity profiles with near-chance real-versus-synthetic separation.", BLUE, PALE_BLUE),
        (4.44, "2", "Use", "Tissue-specific consensus ranking and light synthetic regularization improved held-out prediction in selected tissues.", TEAL, PALE_TEAL),
        (8.46, "3", "Interpret", "Promoted and reinforced genes supported tissue-specific hypotheses. Literature review separated prior alignment from complementary findings.", ORANGE, PALE_GOLD),
    ]
    for index, (x, number, heading, body, color, _fill) in enumerate(columns):
        if index:
            _add_rule(slide, x - 0.20, 2.10, 0.012, "D8DFE3", 3.53)
        _add_rule(slide, x, 2.10, 3.55, color, 0.035)
        _add_text(slide, f"0{number}", x + 0.02, 2.33, 0.68, 0.40, size=21, color=color, bold=True, margin=0)
        _add_text(slide, heading, x + 0.86, 2.35, 2.46, 0.36, size=19, color=NAVY, bold=True, margin=0)
        _add_text(slide, body, x + 0.03, 3.12, 3.28, 2.10, size=16, color=DARK, margin=0)


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
        "Target length: 12-15 minutes. Planned speaking time: about 15 minutes.",
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
    architecture_figure = ASSET_DIR / "lacan_figure1c_generator_architecture.png"
    if not architecture_figure.exists():
        raise FileNotFoundError(
            "Missing Lacan et al. architecture excerpt; see presentation/generative_slstp_2026/assets/SOURCES.md"
        )
    presentation = Presentation(TEMPLATE)
    _set_title_slide(presentation.slides[0])
    _prepare_content_slide(presentation.slides[1], 2)
    while len(presentation.slides) < 22:
        number = len(presentation.slides) + 1
        slide = presentation.slides.add_slide(presentation.slide_layouts[3])
        _prepare_content_slide(slide, number)

    builders = [
        None,
        _slide_why_synthetic,
        _slide_2,
        _slide_scientific_objective,
        _slide_3,
        lambda slide: _slide_4(slide, architecture_figure),
        lambda slide: _slide_5(slide, trajectory),
        _slide_6,
        _slide_7,
        _slide_guidance_mechanism,
        _slide_guidance_boundary,
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
        SlideNote(1, "Synthetic transcriptomics for mouse spaceflight", "0:20", "Introduce the central question: can generated expression help us find tissue-specific FLT versus GC biology? Synthetic profiles support the analysis, but they do not count as additional animals."),
        SlideNote(2, "What is synthetic transcriptomics?", "0:40", "Synthetic transcriptomes are numeric gene-expression vectors sampled from a model trained on measured RNA-seq data. Conditioning lets us request a tissue and FLT or GC context. These profiles can stress-test a classifier and guide gene ranking. They do not add biological replicates or independent evidence."),
        SlideNote(3, "Small studies and study effects complicate tissue comparisons", "0:45", "OSDR gives broad tissue coverage, but the profiles are spread across 75 accessions with different mission and assay contexts. ARCHS4 supplies a much larger mouse reference. The challenge is to use that reference without confusing study structure with spaceflight biology."),
        SlideNote(4, "Match tissue distributions, then test FLT versus GC biology", "0:45", "The generator has two jobs. First, synthetic bulk RNA-seq should reproduce the tissue-defined distributions seen in real data. Second, it should preserve the smaller FLT versus GC signal within each tissue and improve prediction on held-out real samples. Gene effects and BH FDR are always calculated from observed OSDR profiles."),
        SlideNote(5, "We built a configurable bulk RNA-seq generation pipeline", "1:00", "Each column shows the alternatives available at one pipeline stage; outlines identify the downstream branch. We used ARCHS4 and NASA OSDR across multiple studies and all tissues. The selected path used TPM, training-fitted MaxAbs scaling, 974 mouse landmarks, no global correction, ARCHS4 pretraining plus OSDR adaptation, and a DDIM conditioned on tissue, FLT/GC, accession and material. WGAN-GP and the other preprocessing and harmonization choices remained benchmark alternatives."),
        SlideNote(6, "DDIM matched expression and reduced separability", "0:50", "The left image is Figure 1C from Lacan and colleagues. Their residual dense denoiser predicts the noise added to a sample using diffusion timestep and tissue. Our implementation adds FLT/GC, accession and material context during OSDR LoRA adaptation. Both generators had high correlation and F1. DDIM had adversarial accuracy near 0.5 and a lower Frechet-distance ratio, so it was harder to separate from real profiles and closer in distribution."),
        SlideNote(7, "Diffusion learns tissue structure from noise", "0:40", "Read the panels from left to right. The same generated profiles begin as noise, develop structure by timestep 200 and approach tissue-conditioned regions at timestep zero. The axes are shared, so the visual change does not come from rescaling each panel."),
        SlideNote(8, "Generated profiles track the real OSDR PCA manifold", "0:40", "Circles are locked real OSDR profiles and crosses are matched DDIM profiles in the same PCA space. Generated samples follow the tissue-defined branches. FLT and GC overlap more because condition effects are smaller than tissue effects. The numerical validation on the previous slide tests fidelity directly."),
        SlideNote(9, "Five arms separate gene ranking from classifier fitting", "0:50", "Each arm makes two decisions: which profiles rank the genes and which profiles fit the classifier. In both guided arms, real and synthetic evidence jointly rank genes. Guided real fit then trains only on observed profiles. Guided 5% also uses condition-recentered synthetic profiles, but they contribute only 5% of total classifier weight. Held-out real profiles determine eligibility, and FLT/GC effects and BH FDR come from observed OSDR profiles only."),
        SlideNote(10, "Consensus ranking chooses the classifier input genes", "0:35", "Each ranking orders the same 974 genes. Real-only, generated-only and consensus ranking can therefore produce different candidate gene sets. For each ranking, we test the top 10, 25, 50 and 100 genes, and held-out validation chooses the feature count and regularization. Logistic regression is then fitted using only those selected gene-expression columns. Ranking chooses which genes are available to the classifier; classifier training separately learns their coefficients and decision boundary."),
        SlideNote(11, "Held-out real samples decide whether guidance helps", "0:25", "Opaque profiles belong to the training subset. Transparent profiles are held-out real OSDR samples that never enter gene ranking, classifier fitting or top-k selection. The same held-out profiles are scored by the real-only and synthetic-guided candidates. This schematic shows the intended outcome: a guided boundary that predicts more held-out labels correctly. Balanced accuracy, AUROC and average precision determine whether a synthetic arm is eligible for that tissue; the observed tissue-specific results follow on the next slide."),
        SlideNote(12, "Pooling tissues hid useful signal", "0:50", "The pooled augmentation test was negative: balanced accuracy fell from 0.754 to 0.737 with real plus synthetic training. Tissue-specific analysis changed the result. Different tissues benefited from different synthetic uses, which argues against one global augmentation policy."),
        SlideNote(13, "Synthetic guidance changed ranking, not statistical evidence", "0:45", "The blue set contains genes selected stably by real-only ranking, and the teal set contains genes selected stably by the eligible synthetic-guided arm. Thirty-four were real-only, 23 were selected by both arms and classified as reinforced, and 26 were selected only with guidance and classified as promoted. Promoted does not mean biologically novel. All 49 synthetic-informed tissue-gene associations passed BH FDR in observed OSDR profiles."),
        SlideNote(14, "Selection and literature are separate dimensions", "0:50", "Every association has two labels. Promoted or reinforced describes repeated feature selection. Aligning, complementary, ambiguous or unmatched describes prior literature. Across all 49 associations, 22 aligned, 19 were complementary, four were ambiguous and four were unmatched. Table S16 records the gene-level rationale and source IDs; Table S17 records the citations and evidence relationship."),
        SlideNote(15, "The screen covered all 27 completed tissue analyses", "0:35", "This is the full analysis coverage: 22 canonical tissues and five anatomical muscle groups. Ten had a synthetic-informed BH-FDR association, five had real BH-FDR genes without synthetic-informed selection, and 12 had no BH-FDR gene in the landmark panel. Every completed tissue result remains visible here."),
        SlideNote(16, "Ten tissue analyses contained synthetic-informed genes", "0:40", "This is the complete 49-association inventory. Separate rows show FLT-higher or FLT-lower direction and promoted or reinforced selection status. Gene color independently shows aligning, complementary, ambiguous or unmatched literature. FLT directions come from real-data meta-analysis."),
        SlideNote(17, "Thymus points to lower proliferative renewal", "1:00", "Thymus produced the clearest promoted panel. The lower mitotic and DNA-replication genes agree with prior reports of thymic involution and altered cell-cycle expression after flight. Higher Hsd17b11 and Etv1 add lipid-handling and T-cell-state hypotheses. A matched sensitivity model without material-type conditioning preserved the cell-cycle interpretation. Because this is bulk RNA-seq, the pattern may reflect transcription, cell composition or both."),
        SlideNote(18, "Soleus reinforces a mitochondrial and lipid program", "0:55", "Soleus improved with real plus generated training. The selected genes were already stable in real-only analysis, so synthetic data reinforced rather than introduced the panel. Lower Bdh1, Ech1, Bnip3 and Decr1, with higher Tpm1, support altered oxidative metabolism and contractile remodeling. In the no-material sensitivity model, the synthetic attribution disappeared. The real OSDR association remains, but the generated-profile contribution is conditioning-sensitive."),
        SlideNote(19, "Additional tissues produced distinct hypotheses", "0:35", "Promoted and reinforced genes are shown on separate subrows for each tissue. Pooled muscle, kidney, spleen and skin each produced a distinct synthetic-informed result. The rows share a slide for presentation space; each remains a separate hypothesis."),
        SlideNote(20, "Eye, adrenal and muscle-group results remain tissue-specific", "0:35", "Promoted and reinforced genes remain separated here as well. Eye reinforces lower cytokinesis, adrenal contributes two unmatched candidates, gastrocnemius combines an NF-kappa-B stress signal with an autophagy or myogenesis candidate, and tibialis anterior spans stress, cell-cycle, ganglioside and mitophagy hypotheses."),
        SlideNote(21, "Synthetic data worked best as a tissue-specific prior", "0:35", "Conditional DDIM generated realistic expression profiles. Tissue-specific consensus ranking and light synthetic regularization improved held-out prediction in selected tissues. Synthetic-informed selection prioritized promoted and reinforced genes, while literature annotation separated prior alignment from complementary hypotheses. Biological evidence and FDR remained grounded in observed OSDR profiles."),
        SlideNote(22, "Thank you", "0:10", "Acknowledge the mentor, SLSTP, NASA OSDR, ARCHS4, Reactome and NASA Ames compute. Invite questions."),
    ]
    for note, slide in zip(notes, presentation.slides):
        _add_notes(slide, note)
    _write_notes(notes)
    presentation.save(OUTPUT)
    return OUTPUT


if __name__ == "__main__":
    print(build())
