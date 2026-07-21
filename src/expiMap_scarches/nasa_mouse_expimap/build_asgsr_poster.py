"""Build and render a 48 x 36 inch ASGSR expiMap scientific poster."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt


ROOT = Path(__file__).resolve().parents[3]
PAPER_DIR = ROOT / "paper/asgsr_expimap_hvg"
FIGURE_DIR = PAPER_DIR / "figures"
POSTER_DIR = PAPER_DIR / "poster"
ASSET_DIR = POSTER_DIR / "assets"

TITLE = (
    "Cross-mission expiMap analysis recovers established tissue responses "
    "and identifies complementary pathway shifts in mouse spaceflight "
    "transcriptomes"
)
SLIDE_W = 48.0
SLIDE_H = 36.0
FONT = "Arial"

WHITE = "FFFFFF"
INK = "182126"
MUTED = "59676D"
LIGHT = "F1F4F5"
RULE = "CDD5D8"
BLUE = "1D6FA5"
BLUE_DARK = "174B6B"
BLUE_PALE = "E7F1F6"
GREEN = "087F6C"
ORANGE = "D96725"
THYMUS = "6A4C93"
SKIN = "C55A35"
LIVER = "0E8177"
SPLEEN = "A64D67"


def rgb(value: str) -> RGBColor:
    return RGBColor.from_string(value)


def fill(shape, color: str) -> None:
    shape.fill.solid()
    shape.fill.fore_color.rgb = rgb(color)


def outline(shape, color: str, width: float = 1.0) -> None:
    shape.line.color.rgb = rgb(color)
    shape.line.width = Pt(width)


def text(
    slide,
    value: str,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    size: float,
    color: str = INK,
    bold: bool = False,
    italic: bool = False,
    align=PP_ALIGN.LEFT,
    valign=MSO_ANCHOR.TOP,
    margin: float = 0.0,
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
    paragraph.alignment = align
    paragraph.space_before = Pt(0)
    paragraph.space_after = Pt(0)
    paragraph.line_spacing = 1.0
    run = paragraph.add_run()
    run.text = value
    run.font.name = FONT
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = rgb(color)
    return shape


def paragraphs(
    slide,
    values: list[str],
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    size: float,
    color: str = INK,
    gap: float = 5.0,
):
    shape = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    frame = shape.text_frame
    frame.clear()
    frame.word_wrap = True
    frame.margin_left = 0
    frame.margin_right = 0
    frame.margin_top = 0
    frame.margin_bottom = 0
    for index, value in enumerate(values):
        paragraph = frame.paragraphs[0] if index == 0 else frame.add_paragraph()
        paragraph.space_before = Pt(0)
        paragraph.space_after = Pt(gap)
        paragraph.line_spacing = 1.0
        run = paragraph.add_run()
        run.text = value
        run.font.name = FONT
        run.font.size = Pt(size)
        run.font.color.rgb = rgb(color)
    return shape


def rule(slide, x: float, y: float, w: float, color: str = RULE) -> None:
    shape = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(0.035)
    )
    fill(shape, color)
    shape.line.fill.background()


def section(slide, label: str, x: float, y: float, w: float, color: str) -> None:
    marker = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(0.16), Inches(0.54)
    )
    fill(marker, color)
    marker.line.fill.background()
    text(
        slide,
        label,
        x + 0.30,
        y - 0.02,
        w - 0.30,
        0.58,
        size=25,
        bold=True,
        valign=MSO_ANCHOR.MIDDLE,
    )
    rule(slide, x, y + 0.67, w)


def labeled_box(
    slide,
    title: str,
    subtitle: str,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    box_fill: str,
    box_line: str,
) -> None:
    shape = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h)
    )
    fill(shape, box_fill)
    outline(shape, box_line, 1.5)
    frame = shape.text_frame
    frame.clear()
    frame.word_wrap = True
    frame.margin_left = Inches(0.14)
    frame.margin_right = Inches(0.14)
    frame.margin_top = Inches(0.08)
    frame.margin_bottom = Inches(0.06)
    frame.vertical_anchor = MSO_ANCHOR.MIDDLE
    first = frame.paragraphs[0]
    first.alignment = PP_ALIGN.CENTER
    first.space_after = Pt(3)
    run = first.add_run()
    run.text = title
    run.font.name = FONT
    run.font.size = Pt(21)
    run.font.bold = True
    run.font.color.rgb = rgb(INK)
    second = frame.add_paragraph()
    second.alignment = PP_ALIGN.CENTER
    second.space_before = Pt(0)
    second.space_after = Pt(0)
    run = second.add_run()
    run.text = subtitle
    run.font.name = FONT
    run.font.size = Pt(16.5)
    run.font.color.rgb = rgb(MUTED)


def down_arrow(slide, center_x: float, y: float) -> None:
    shape = slide.shapes.add_shape(
        MSO_SHAPE.DOWN_ARROW,
        Inches(center_x - 0.25),
        Inches(y),
        Inches(0.50),
        Inches(0.43),
    )
    fill(shape, BLUE)
    shape.line.fill.background()


def badge(
    slide, label: str, x: float, y: float, w: float, color: str = BLUE
) -> None:
    shape = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE,
        Inches(x),
        Inches(y),
        Inches(w),
        Inches(0.77),
    )
    fill(shape, BLUE_PALE if color == BLUE else "EDF5F2")
    outline(shape, color, 1.0)
    frame = shape.text_frame
    frame.clear()
    frame.word_wrap = True
    frame.margin_left = Inches(0.06)
    frame.margin_right = Inches(0.06)
    frame.margin_top = 0
    frame.margin_bottom = 0
    frame.vertical_anchor = MSO_ANCHOR.MIDDLE
    paragraph = frame.paragraphs[0]
    paragraph.alignment = PP_ALIGN.CENTER
    run = paragraph.add_run()
    run.text = label
    run.font.name = FONT
    run.font.size = Pt(15.5)
    run.font.bold = True
    run.font.color.rgb = rgb(BLUE_DARK)


def set_cell(cell, value: str, *, size: float, color: str, bold: bool, center: bool):
    cell.text = ""
    cell.margin_left = Inches(0.08)
    cell.margin_right = Inches(0.06)
    cell.margin_top = Inches(0.04)
    cell.margin_bottom = Inches(0.04)
    cell.text_frame.vertical_anchor = MSO_ANCHOR.MIDDLE
    paragraph = cell.text_frame.paragraphs[0]
    paragraph.alignment = PP_ALIGN.CENTER if center else PP_ALIGN.LEFT
    run = paragraph.add_run()
    run.text = value
    run.font.name = FONT
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = rgb(color)


def scope_table(slide, x: float, y: float, w: float, h: float) -> None:
    rows = [
        ("Thymus", "1,362", "117", "5", "387", THYMUS),
        ("Skin", "2,593", "151", "4", "319", SKIN),
        ("Liver", "5,000", "197", "9", "364", LIVER),
        ("Spleen", "6,289", "100", "5", "360", SPLEEN),
    ]
    shape = slide.shapes.add_table(5, 5, Inches(x), Inches(y), Inches(w), Inches(h))
    table = shape.table
    widths = [3.15, 3.00, 2.55, 2.45, 2.70]
    for index, value in enumerate(widths):
        table.columns[index].width = Inches(value)
    for column, value in enumerate(("Tissue", "ARCHS4", "OSDR", "Projects", "Programs")):
        cell = table.cell(0, column)
        cell.fill.solid()
        cell.fill.fore_color.rgb = rgb(BLUE_DARK)
        set_cell(cell, value, size=16, color=WHITE, bold=True, center=column > 0)
    for row_index, row in enumerate(rows, 1):
        for column, value in enumerate(row[:5]):
            cell = table.cell(row_index, column)
            cell.fill.solid()
            cell.fill.fore_color.rgb = rgb(WHITE if row_index % 2 else LIGHT)
            set_cell(
                cell,
                value,
                size=16.5,
                color=row[5] if column == 0 else INK,
                bold=column == 0,
                center=column > 0,
            )
    for row in table.rows:
        row.height = Inches(h / 5)


def render_pdf(source: Path, output: Path, dpi: int = 700) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "pdftocairo",
            "-singlefile",
            "-png",
            "-r",
            str(dpi),
            str(source),
            str(output.with_suffix("")),
        ],
        check=True,
    )


def picture_contain(
    slide,
    path: Path,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    name: str,
    description: str,
) -> float:
    with Image.open(path) as image:
        px_w, px_h = image.size
    image_ratio = px_w / px_h
    region_ratio = w / h
    if image_ratio >= region_ratio:
        placed_w = w
        placed_h = w / image_ratio
        placed_x = x
        placed_y = y + (h - placed_h) / 2
    else:
        placed_h = h
        placed_w = h * image_ratio
        placed_x = x + (w - placed_w) / 2
        placed_y = y
    picture = slide.shapes.add_picture(
        str(path),
        Inches(placed_x),
        Inches(placed_y),
        Inches(placed_w),
        Inches(placed_h),
    )
    properties = picture._element.xpath(".//p:cNvPr")[0]
    properties.set("name", name)
    properties.set("descr", description)
    ppi = min(px_w / placed_w, px_h / placed_h)
    if ppi < 300:
        raise ValueError(f"{path.name} has only {ppi:.0f} effective ppi")
    return ppi


def evidence_row(
    slide,
    tissue: str,
    color: str,
    finding: str,
    evidence: str,
    x: float,
    y: float,
    w: float,
) -> None:
    h = 2.25
    rule(slide, x, y + h - 0.03, w)
    marker = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(x),
        Inches(y + 0.12),
        Inches(0.16),
        Inches(h - 0.30),
    )
    fill(marker, color)
    marker.line.fill.background()
    text(slide, tissue.upper(), x + 0.34, y + 0.05, 2.15, 0.45, size=18, color=color, bold=True)
    text(slide, finding, x + 0.34, y + 0.52, w - 4.0, 1.55, size=17.5)
    text(
        slide,
        evidence,
        x + w - 3.40,
        y + 0.10,
        3.22,
        1.95,
        size=15.5,
        color=MUTED,
        bold=True,
        align=PP_ALIGN.RIGHT,
        valign=MSO_ANCHOR.MIDDLE,
    )


def story_row(
    slide, tissue: str, color: str, finding: str, x: float, y: float, w: float
) -> None:
    dot = slide.shapes.add_shape(
        MSO_SHAPE.OVAL, Inches(x), Inches(y + 0.12), Inches(0.34), Inches(0.34)
    )
    fill(dot, color)
    dot.line.fill.background()
    text(slide, tissue + ":", x + 0.52, y, 1.80, 1.42, size=18, color=color, bold=True)
    text(slide, finding, x + 2.05, y, w - 2.05, 1.42, size=17.5)


def validate(prs: Presentation) -> None:
    if len(prs.slides) != 1:
        raise ValueError("Poster must contain exactly one slide")
    for shape in prs.slides[0].shapes:
        if shape.left < 0 or shape.top < 0:
            raise ValueError(f"{shape.name!r} starts outside the slide")
        if shape.left + shape.width > prs.slide_width + 1:
            raise ValueError(f"{shape.name!r} exceeds the slide width")
        if shape.top + shape.height > prs.slide_height + 1:
            raise ValueError(f"{shape.name!r} exceeds the slide height")


def render_poster(pptx_path: Path) -> tuple[Path | None, Path | None]:
    libreoffice = shutil.which("libreoffice")
    pdftocairo = shutil.which("pdftocairo")
    if not libreoffice:
        return None, None
    pdf_path = pptx_path.with_suffix(".pdf")
    pdf_path.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix="asgsr_poster_lo_") as profile:
        subprocess.run(
            [
                libreoffice,
                "--headless",
                f"-env:UserInstallation={Path(profile).resolve().as_uri()}",
                "--convert-to",
                "pdf",
                "--outdir",
                str(pptx_path.parent),
                str(pptx_path),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if not pdf_path.exists():
        raise FileNotFoundError(f"LibreOffice did not create {pdf_path}")
    if not pdftocairo:
        return pdf_path, None
    preview = pptx_path.with_name(pptx_path.stem + "_preview.png")
    subprocess.run(
        [
            pdftocairo,
            "-singlefile",
            "-png",
            "-r",
            "100",
            str(pdf_path),
            str(preview.with_suffix("")),
        ],
        check=True,
    )
    return pdf_path, preview


def build() -> tuple[Path, Path | None, Path | None, list[float]]:
    POSTER_DIR.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    result_asset = ASSET_DIR / "figure_3_tissue_pathway_shifts_700dpi.png"
    hypothesis_asset = ASSET_DIR / "figure_6_tissue_state_hypotheses_700dpi.png"
    render_pdf(FIGURE_DIR / "figure_3_tissue_pathway_shifts.pdf", result_asset)
    render_pdf(FIGURE_DIR / "figure_6_tissue_state_hypotheses.pdf", hypothesis_asset)

    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)
    prs.core_properties.title = TITLE
    prs.core_properties.author = "Jason Trinh"
    prs.core_properties.subject = (
        "Cross-mission expiMap analysis of NASA OSDR mouse transcriptomes"
    )
    prs.core_properties.keywords = (
        "spaceflight, expiMap, OSDR, ARCHS4, Reactome, transcriptomics"
    )
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = rgb(WHITE)

    # Header
    bar = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, 0, 0, Inches(SLIDE_W), Inches(0.22)
    )
    fill(bar, BLUE_DARK)
    bar.line.fill.background()
    for index, color in enumerate((THYMUS, SKIN, LIVER, SPLEEN)):
        segment = slide.shapes.add_shape(
            MSO_SHAPE.RECTANGLE,
            Inches(index * 12.0),
            Inches(0.22),
            Inches(12.0),
            Inches(0.08),
        )
        fill(segment, color)
        segment.line.fill.background()

    text(
        slide,
        "ASGSR 2026  |  NASA OPEN SCIENCE DATA REPOSITORY",
        0.78,
        0.48,
        30.0,
        0.40,
        size=16,
        color=BLUE,
        bold=True,
        valign=MSO_ANCHOR.MIDDLE,
    )
    text(
        slide,
        TITLE.replace(" and identifies", "\nand identifies"),
        0.76,
        0.95,
        46.4,
        2.15,
        size=46,
        bold=True,
        valign=MSO_ANCHOR.MIDDLE,
    )
    text(
        slide,
        (
            "Jason Trinh  |  NASA Space Life Sciences Training Program, "
            "NASA Ames Research Center"
        ),
        0.78,
        3.13,
        32.4,
        0.50,
        size=20,
        bold=True,
        valign=MSO_ANCHOR.MIDDLE,
    )
    text(
        slide,
        "jasontrinh@berkeley.edu\ngithub.com/jasont314/nasa-mouse",
        33.4,
        3.04,
        13.8,
        0.72,
        size=17,
        color=BLUE_DARK,
        bold=True,
        align=PP_ALIGN.RIGHT,
        valign=MSO_ANCHOR.MIDDLE,
    )

    take_home = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(0.72),
        Inches(3.92),
        Inches(46.55),
        Inches(1.26),
    )
    fill(take_home, BLUE_PALE)
    take_home.line.fill.background()
    accent = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(0.72),
        Inches(3.92),
        Inches(0.24),
        Inches(1.26),
    )
    fill(accent, BLUE)
    accent.line.fill.background()
    text(
        slide,
        "TAKE-HOME",
        1.22,
        4.10,
        3.05,
        0.44,
        size=19,
        color=BLUE_DARK,
        bold=True,
        valign=MSO_ANCHOR.MIDDLE,
    )
    text(
        slide,
        (
            "Cross-mission pathway mapping recovered established tissue responses "
            "and prioritized complementary maintenance and immune hypotheses; "
            "spleen showed the strongest multi-method evidence."
        ),
        4.38,
        4.04,
        42.35,
        0.86,
        size=22,
        bold=True,
        valign=MSO_ANCHOR.MIDDLE,
    )

    # Main three-column grid
    col1_x, col1_w = 0.72, 14.45
    col2_x, col2_w = 15.82, 15.28
    col3_x, col3_w = 31.77, 15.50
    body_top = 5.62
    body_bottom = 32.38
    for x in (15.49, 31.43):
        divider = slide.shapes.add_shape(
            MSO_SHAPE.RECTANGLE,
            Inches(x),
            Inches(body_top),
            Inches(0.025),
            Inches(body_bottom - body_top),
        )
        fill(divider, RULE)
        divider.line.fill.background()

    # Column 1: question, workflow, scope, and validation
    section(slide, "QUESTION & APPROACH", col1_x, body_top, col1_w, BLUE)
    text(
        slide,
        (
            "Can a tissue-matched, pathway-constrained reference distinguish "
            "recurring flight-associated programs from mission-specific variation?"
        ),
        col1_x,
        6.48,
        col1_w,
        1.32,
        size=20,
        bold=True,
    )
    text(
        slide,
        (
            "Public mouse RNA-seq studies differ in mission, strain, sex, duration, "
            "hardware, collection endpoint, and sequencing protocol."
        ),
        col1_x,
        7.82,
        col1_w,
        1.05,
        size=17.5,
        color=MUTED,
    )
    text(
        slide,
        "REFERENCE-MAPPING WORKFLOW",
        col1_x,
        8.82,
        col1_w,
        0.42,
        size=18,
        color=BLUE_DARK,
        bold=True,
    )
    labeled_box(
        slide,
        "ARCHS4",
        "tissue-matched non-spaceflight RNA-seq",
        col1_x,
        9.35,
        6.90,
        1.42,
        box_fill="E9F2F6",
        box_line=BLUE,
    )
    labeled_box(
        slide,
        "Reactome",
        "current mouse gene-program mask",
        col1_x + 7.25,
        9.35,
        7.20,
        1.42,
        box_fill="EAF4EF",
        box_line=GREEN,
    )
    down_arrow(slide, col1_x + col1_w / 2, 10.86)
    labeled_box(
        slide,
        "expiMap reference",
        "negative binomial model | approximately 2,000 HVGs",
        col1_x + 1.20,
        11.34,
        col1_w - 2.40,
        1.50,
        box_fill="EEEAF5",
        box_line=THYMUS,
    )
    down_arrow(slide, col1_x + col1_w / 2, 12.93)
    labeled_box(
        slide,
        "NASA OSDR query",
        "flight and ground samples | accession conditioned",
        col1_x + 1.20,
        13.42,
        col1_w - 2.40,
        1.50,
        box_fill="F8EEE8",
        box_line=ORANGE,
    )
    down_arrow(slide, col1_x + col1_w / 2, 15.01)
    labeled_box(
        slide,
        "Project-balanced pathway shifts",
        "decoder-oriented mean flight minus ground score",
        col1_x + 1.20,
        15.50,
        col1_w - 2.40,
        1.50,
        box_fill=BLUE_PALE,
        box_line=BLUE,
    )

    text(
        slide,
        "PRIMARY MODEL SCOPE",
        col1_x,
        17.40,
        col1_w,
        0.44,
        size=18,
        color=BLUE_DARK,
        bold=True,
    )
    scope_table(slide, col1_x, 17.95, col1_w, 4.25)
    text(
        slide,
        (
            "OSDR counts are primary analysis samples; spleen excludes a "
            "condition-strain-confounded project."
        ),
        col1_x,
        22.30,
        col1_w,
        0.65,
        size=14.5,
        color=MUTED,
        italic=True,
    )

    text(
        slide,
        "ROBUSTNESS GATE",
        col1_x,
        23.10,
        col1_w,
        0.44,
        size=18,
        color=BLUE_DARK,
        bold=True,
    )
    badge_w = 4.53
    badge_gap = 0.40
    for index, label in enumerate(("ssGSEA", "Preranked GSEA", "Held-out projects")):
        badge(
            slide,
            label,
            col1_x + index * (badge_w + badge_gap),
            23.68,
            badge_w,
        )
    for index, label in enumerate(
        ("3 full trainings", "Composition proxies", "Member genes")
    ):
        badge(
            slide,
            label,
            col1_x + index * (badge_w + badge_gap),
            24.66,
            badge_w,
            THYMUS if index == 0 else GREEN,
        )
    text(
        slide,
        (
            "Pathways were interpreted after project direction, complete retraining, "
            "conventional enrichment, composition sensitivity, and member-gene review."
        ),
        col1_x,
        25.73,
        col1_w,
        1.33,
        size=17,
    )
    training = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(col1_x),
        Inches(27.20),
        Inches(col1_w),
        Inches(1.18),
    )
    fill(training, LIGHT)
    training.line.fill.background()
    text(
        slide,
        "TRAINING",
        col1_x + 0.24,
        27.38,
        1.80,
        0.38,
        size=15.5,
        color=BLUE_DARK,
        bold=True,
    )
    text(
        slide,
        (
            "3 x 300 hidden units  |  reference <=400 epochs  |  "
            "query 250 epochs  |  NVIDIA A100"
        ),
        col1_x + 2.00,
        27.30,
        col1_w - 2.20,
        0.60,
        size=15.5,
        valign=MSO_ANCHOR.MIDDLE,
    )
    text(
        slide,
        (
            "Complementary programs are annotated Reactome nodes that extend prior "
            "literature; de novo nodes were not retained in the final models."
        ),
        col1_x,
        28.63,
        col1_w,
        1.40,
        size=16.5,
        color=MUTED,
    )
    scope_note = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(col1_x),
        Inches(30.20),
        Inches(col1_w),
        Inches(1.70),
    )
    fill(scope_note, "FFF5E8")
    scope_note.line.fill.background()
    text(
        slide,
        "ANALYSIS SCOPE",
        col1_x + 0.24,
        30.39,
        2.65,
        0.40,
        size=15.5,
        color=ORANGE,
        bold=True,
    )
    text(
        slide,
        (
            "Kidney remained exploratory. Soleus, lung, and retina were not "
            "advanced to the primary biological narrative."
        ),
        col1_x + 0.24,
        30.84,
        col1_w - 0.48,
        0.78,
        size=16,
    )

    # Column 2: pathway shifts and evidence
    section(slide, "PATHWAY SHIFTS", col2_x, body_top, col2_w, BLUE)
    text(
        slide,
        "Retained decoder-oriented flight minus ground effects",
        col2_x,
        6.48,
        col2_w,
        0.48,
        size=19,
        bold=True,
        align=PP_ALIGN.CENTER,
    )
    result_ppi = picture_contain(
        slide,
        result_asset,
        col2_x,
        7.10,
        col2_w,
        12.84,
        name="Retained pathway shifts",
        description=(
            "Project-balanced expiMap pathway shifts for thymus, skin, liver, and "
            "spleen, with project effects and three-training ranges."
        ),
    )
    text(
        slide,
        (
            "Open circles: OSDR projects. Colored ranges: three complete "
            "reference-query trainings. Axes are tissue specific."
        ),
        col2_x + 0.25,
        19.90,
        col2_w - 0.50,
        0.66,
        size=14.5,
        color=MUTED,
        italic=True,
        align=PP_ALIGN.CENTER,
    )
    text(
        slide,
        "EVIDENCE AT A GLANCE",
        col2_x,
        20.75,
        col2_w,
        0.44,
        size=18,
        color=BLUE_DARK,
        bold=True,
    )
    evidence = [
        (
            "Thymus",
            THYMUS,
            (
                "Lower DNA repair and cytoskeletal scores; expiMap also suggests "
                "lower niche interaction, but GSEA does not."
            ),
            "5/5 projects\n3/3 trainings\nFDR <0.001 / 0.154 / 1.000",
        ),
        (
            "Skin",
            SKIN,
            "Lower regulation, repair, Hedgehog, sphingolipid, and junction scores.",
            "3/4 projects\n3/3 trainings\n3 of 5 FDR <0.05",
        ),
        (
            "Liver",
            LIVER,
            (
                "Lower MHC class II antigen presentation and T-cell receptor "
                "scores amid mixed metabolism."
            ),
            "8/9 projects\n3/3 trainings\nFDR 0.121 / 0.051",
        ),
        (
            "Spleen",
            SPLEEN,
            (
                "Lower T-cell receptor, degranulation, and C-type lectin "
                "receptor programs."
            ),
            "5/5 projects\n3/3 trainings\nall FDR <0.05",
        ),
    ]
    y = 21.30
    for tissue, color, finding, summary in evidence:
        evidence_row(slide, tissue, color, finding, summary, col2_x, y, col2_w)
        y += 2.25
    text(
        slide,
        (
            "FDR values are from preranked GSEA. Directional robustness and FDR "
            "are complementary evidence, not interchangeable thresholds."
        ),
        col2_x,
        30.55,
        col2_w,
        1.18,
        size=15.5,
        color=MUTED,
        italic=True,
        align=PP_ALIGN.CENTER,
    )

    # Column 3: tissue-state interpretation and next tests
    section(
        slide,
        "BIOLOGICAL INTERPRETATION",
        col3_x,
        body_top,
        col3_w,
        GREEN,
    )
    hypothesis_ppi = picture_contain(
        slide,
        hypothesis_asset,
        col3_x + 0.10,
        6.45,
        col3_w - 0.20,
        12.70,
        name="Tissue-state hypotheses",
        description=(
            "Observed lower expiMap program scores and qualified tissue-state "
            "hypotheses for thymus, skin, liver, and spleen."
        ),
    )
    text(
        slide,
        "WHAT THIS ADDS TO PRIOR LITERATURE",
        col3_x,
        19.55,
        col3_w,
        0.44,
        size=18,
        color=BLUE_DARK,
        bold=True,
    )
    stories = [
        (
            "Thymus",
            THYMUS,
            (
                "Repair and cytoskeletal reduction adds a possible "
                "niche-coordination layer to known involution."
            ),
        ),
        (
            "Skin",
            SKIN,
            (
                "Multiple programs converge on lower tissue maintenance and "
                "barrier coordination."
            ),
        ),
        (
            "Liver",
            LIVER,
            (
                "A lower immune-communication axis complements the better-known "
                "metabolic response."
            ),
        ),
        (
            "Spleen",
            SPLEEN,
            (
                "Strongest evidence: adaptive activation and innate sensing and "
                "effector transcription decrease together."
            ),
        ),
    ]
    y = 20.20
    for tissue, color, finding in stories:
        story_row(slide, tissue, color, finding, col3_x, y, col3_w)
        y += 1.62

    text(
        slide,
        "INTERPRETATION AND NEXT TESTS",
        col3_x,
        26.88,
        col3_w,
        0.44,
        size=18,
        color=BLUE_DARK,
        bold=True,
    )
    paragraphs(
        slide,
        [
            "\u2022 Lower means a lower decoder-oriented latent score, not biochemical inhibition.",
            "\u2022 Bulk tissue cannot distinguish cell abundance from altered cell state.",
            (
                "\u2022 Validate with independent single-cell, spatial, biochemical, "
                "and functional studies."
            ),
        ],
        col3_x,
        27.45,
        col3_w,
        2.60,
        size=16.5,
    )
    conclusion = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(col3_x),
        Inches(30.28),
        Inches(col3_w),
        Inches(1.62),
    )
    fill(conclusion, "EAF4EF")
    conclusion.line.fill.background()
    text(
        slide,
        "CONCLUSION",
        col3_x + 0.24,
        30.46,
        2.15,
        0.40,
        size=15.5,
        color=GREEN,
        bold=True,
    )
    text(
        slide,
        (
            "The cross-mission design narrows broad tissue phenotypes into "
            "testable pathway-level hypotheses without claiming causality or "
            "cell-type resolution."
        ),
        col3_x + 0.24,
        30.87,
        col3_w - 0.48,
        0.78,
        size=16,
    )

    # Footer
    footer_y = 32.55
    rule(slide, 0.72, footer_y, 46.55, BLUE_DARK)
    text(
        slide,
        "SELECTED REFERENCES",
        0.72,
        32.82,
        4.0,
        0.36,
        size=14.5,
        color=BLUE_DARK,
        bold=True,
    )
    references = (
        "1 Gebre et al., Nucleic Acids Res 2025 (OSDR).  "
        "2 Lotfollahi et al., Nat Cell Biol 2023 (expiMap).  "
        "3 Lotfollahi et al., Nat Biotechnol 2022 (scArches).  "
        "4 Lachmann et al., Nat Commun 2018 (ARCHS4).\n"
        "5 Milacic et al., Nucleic Acids Res 2024 (Reactome).  "
        "6 Horie et al., Sci Rep 2019 (thymus).  "
        "7 Cope et al., Commun Med 2024 (skin).  "
        "8 Beheshti et al., Sci Rep 2019 (liver).  "
        "9 Gridley et al., J Appl Physiol 2009 (spleen)."
    )
    text(
        slide,
        references,
        0.72,
        33.22,
        31.0,
        1.68,
        size=13.5,
        color=MUTED,
    )
    text(
        slide,
        "DATA, CODE, AND CONTACT",
        32.20,
        32.82,
        7.0,
        0.36,
        size=14.5,
        color=BLUE_DARK,
        bold=True,
    )
    text(
        slide,
        (
            "NASA OSDR API | ARCHS4 | Reactome\n"
            "github.com/jasont314/nasa-mouse\n"
            "jasontrinh@berkeley.edu"
        ),
        32.20,
        33.22,
        9.15,
        1.62,
        size=13.5,
        bold=True,
    )
    acknowledgment = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE,
        Inches(41.55),
        Inches(32.82),
        Inches(5.72),
        Inches(2.02),
    )
    fill(acknowledgment, BLUE_DARK)
    acknowledgment.line.fill.background()
    text(
        slide,
        "NASA/SLSTP",
        41.80,
        33.06,
        5.20,
        0.55,
        size=22,
        color=WHITE,
        bold=True,
        align=PP_ALIGN.CENTER,
        valign=MSO_ANCHOR.MIDDLE,
    )
    text(
        slide,
        "NASA Ames Research Center\nMoffett Field, California",
        41.80,
        33.65,
        5.20,
        0.82,
        size=13.5,
        color=WHITE,
        align=PP_ALIGN.CENTER,
        valign=MSO_ANCHOR.MIDDLE,
    )

    validate(prs)
    output = POSTER_DIR / "asgsr_expimap_poster.pptx"
    prs.save(output)
    pdf_path, preview_path = render_poster(output)
    return output, pdf_path, preview_path, [result_ppi, hypothesis_ppi]


def run() -> None:
    output, pdf_path, preview_path, ppi = build()
    print(output)
    if pdf_path:
        print(pdf_path)
    if preview_path:
        print(preview_path)
    print(
        "Embedded data-figure effective resolution: "
        + ", ".join(f"{value:.0f} ppi" for value in ppi)
    )


if __name__ == "__main__":
    run()
