"""Build the approved-template ASGSR expiMap scientific poster."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from io import BytesIO
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE, MSO_SHAPE_TYPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt


ROOT = Path(__file__).resolve().parents[3]
PAPER_DIR = ROOT / "paper/asgsr_expimap_hvg"
FIGURE_DIR = PAPER_DIR / "figures"
POSTER_DIR = PAPER_DIR / "poster"
ASSET_DIR = POSTER_DIR / "assets"
TEMPLATE_PATH = (
    ROOT
    / "assets/poster_template/00 Poster Session Student Template_Approved by Legal.pptx"
)

TITLE = (
    "Cross-mission expiMap analysis recovers established tissue responses "
    "and identifies complementary pathway shifts in mouse spaceflight "
    "transcriptomes"
)
DISPLAY_TITLE = (
    "Cross-mission expiMap analysis recovers established tissue responses\n"
    "and identifies complementary pathway shifts in mouse spaceflight transcriptomes"
)

SLIDE_W = 48.0
SLIDE_H = 27.0
FONT = "Arial"

WHITE = "FFFFFF"
INK = "071C3F"
BODY = "33445D"
MUTED = "697889"
HEADER = "DDECF7"
PANEL = "DCEAF4"
PANEL_LIGHT = "EEF4F8"
RULE = "AEBCC5"
NAVY = "082552"
BLUE = "236DA2"
GOLD = "9A6D1D"
GREEN = "178C6A"
ORANGE = "D97800"
THYMUS = "6A4C93"
SKIN = "C65D37"
LIVER = "13877C"
SPLEEN = "A64D67"


def rgb(value: str) -> RGBColor:
    return RGBColor.from_string(value)


def set_fill(shape, color: str) -> None:
    shape.fill.solid()
    shape.fill.fore_color.rgb = rgb(color)


def set_line(shape, color: str, width: float = 1.0) -> None:
    shape.line.color.rgb = rgb(color)
    shape.line.width = Pt(width)


def add_text(
    slide,
    value: str,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    size: float,
    color: str = BODY,
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


def add_paragraphs(
    slide,
    values: list[str],
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    size: float,
    color: str = BODY,
    gap: float = 4.0,
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


def add_rect(
    slide,
    x: float,
    y: float,
    w: float,
    h: float,
    fill: str,
    *,
    line: str | None = None,
    line_width: float = 1.0,
    rounded: bool = False,
):
    kind = MSO_SHAPE.ROUNDED_RECTANGLE if rounded else MSO_SHAPE.RECTANGLE
    shape = slide.shapes.add_shape(
        kind, Inches(x), Inches(y), Inches(w), Inches(h)
    )
    set_fill(shape, fill)
    if line:
        set_line(shape, line, line_width)
    else:
        shape.line.fill.background()
    return shape


def add_panel(slide, x: float, y: float, w: float, h: float) -> None:
    add_rect(slide, x, y, w, h, PANEL)


def add_section(
    slide, label: str, x: float, y: float, w: float, *, size: float = 22
) -> None:
    add_rect(slide, x, y, w, 0.72, WHITE, line=RULE, line_width=1.1)
    add_text(
        slide,
        label,
        x + 0.15,
        y + 0.02,
        w - 0.30,
        0.66,
        size=size,
        color="111111",
        bold=True,
        align=PP_ALIGN.CENTER,
        valign=MSO_ANCHOR.MIDDLE,
    )


def add_badge(
    slide,
    label: str,
    x: float,
    y: float,
    w: float,
    *,
    color: str = BLUE,
) -> None:
    shape = add_rect(
        slide,
        x,
        y,
        w,
        0.72,
        WHITE,
        line=color,
        line_width=1.2,
        rounded=True,
    )
    shape.adjustments[0] = 0.08
    add_text(
        slide,
        label,
        x + 0.08,
        y + 0.03,
        w - 0.16,
        0.63,
        size=15,
        color=color,
        bold=True,
        align=PP_ALIGN.CENTER,
        valign=MSO_ANCHOR.MIDDLE,
    )


def add_box(
    slide,
    title: str,
    subtitle: str,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    fill: str,
    line: str,
) -> None:
    shape = add_rect(
        slide, x, y, w, h, fill, line=line, line_width=1.4, rounded=True
    )
    shape.adjustments[0] = 0.06
    frame = shape.text_frame
    frame.clear()
    frame.word_wrap = True
    frame.margin_left = Inches(0.10)
    frame.margin_right = Inches(0.10)
    frame.margin_top = Inches(0.05)
    frame.margin_bottom = Inches(0.04)
    frame.vertical_anchor = MSO_ANCHOR.MIDDLE
    first = frame.paragraphs[0]
    first.alignment = PP_ALIGN.CENTER
    first.space_after = Pt(2)
    run = first.add_run()
    run.text = title
    run.font.name = FONT
    run.font.size = Pt(16)
    run.font.bold = True
    run.font.color.rgb = rgb(INK)
    second = frame.add_paragraph()
    second.alignment = PP_ALIGN.CENTER
    second.space_before = Pt(0)
    second.space_after = Pt(0)
    run = second.add_run()
    run.text = subtitle
    run.font.name = FONT
    run.font.size = Pt(12.5)
    run.font.color.rgb = rgb(MUTED)


def add_down_arrow(slide, center_x: float, y: float, color: str = BLUE) -> None:
    shape = slide.shapes.add_shape(
        MSO_SHAPE.DOWN_ARROW,
        Inches(center_x - 0.18),
        Inches(y),
        Inches(0.36),
        Inches(0.34),
    )
    set_fill(shape, color)
    shape.line.fill.background()


def add_right_arrow(
    slide, x: float, y: float, w: float, h: float, color: str = BLUE
) -> None:
    shape = slide.shapes.add_shape(
        MSO_SHAPE.RIGHT_ARROW, Inches(x), Inches(y), Inches(w), Inches(h)
    )
    set_fill(shape, color)
    shape.line.fill.background()


def add_connector(
    slide,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    color: str,
    width: float = 1.1,
    dashed: bool = False,
):
    connector = slide.shapes.add_connector(
        MSO_CONNECTOR.STRAIGHT,
        Inches(x1),
        Inches(y1),
        Inches(x2),
        Inches(y2),
    )
    set_line(connector, color, width)
    if dashed:
        connector.line.dash_style = 2
    return connector


def set_table_cell(
    cell,
    value: str,
    *,
    size: float,
    color: str,
    bold: bool,
    center: bool,
) -> None:
    cell.text = ""
    cell.margin_left = Inches(0.07)
    cell.margin_right = Inches(0.06)
    cell.margin_top = Inches(0.03)
    cell.margin_bottom = Inches(0.03)
    cell.text_frame.vertical_anchor = MSO_ANCHOR.MIDDLE
    paragraph = cell.text_frame.paragraphs[0]
    paragraph.alignment = PP_ALIGN.CENTER if center else PP_ALIGN.LEFT
    run = paragraph.add_run()
    run.text = value
    run.font.name = FONT
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = rgb(color)


def add_scope_table(slide, x: float, y: float, w: float, h: float) -> None:
    rows = [
        ("Thymus", "1,362", "117", "5", "387", THYMUS),
        ("Skin", "2,593", "151", "4", "319", SKIN),
        ("Liver", "5,000", "197", "9", "364", LIVER),
        ("Spleen", "6,289", "100", "5", "360", SPLEEN),
    ]
    shape = slide.shapes.add_table(5, 5, Inches(x), Inches(y), Inches(w), Inches(h))
    table = shape.table
    widths = [3.10, 2.75, 2.35, 2.25, 2.65]
    for index, value in enumerate(widths):
        table.columns[index].width = Inches(value)
    for column, value in enumerate(("Tissue", "ARCHS4", "OSDR", "Projects", "Programs")):
        cell = table.cell(0, column)
        cell.fill.solid()
        cell.fill.fore_color.rgb = rgb(NAVY)
        set_table_cell(
            cell,
            value,
            size=14,
            color=WHITE,
            bold=True,
            center=column > 0,
        )
    for row_index, row in enumerate(rows, 1):
        for column, value in enumerate(row[:5]):
            cell = table.cell(row_index, column)
            cell.fill.solid()
            cell.fill.fore_color.rgb = rgb(WHITE if row_index % 2 else PANEL_LIGHT)
            set_table_cell(
                cell,
                value,
                size=14.5,
                color=row[5] if column == 0 else BODY,
                bold=column == 0,
                center=column > 0,
            )
    for row in table.rows:
        row.height = Inches(h / 5)


def add_nasa_logo(slide, template: Presentation) -> None:
    template_slide = template.slides[0]
    candidates = [
        shape
        for shape in template_slide.shapes
        if shape.shape_type == MSO_SHAPE_TYPE.PICTURE
        and shape.left / template.slide_width > 0.85
        and shape.top / template.slide_height < 0.25
    ]
    if not candidates:
        raise RuntimeError("NASA logo was not found in the approved template")
    logo = max(candidates, key=lambda shape: shape.width * shape.height)
    slide.shapes.add_picture(BytesIO(logo.image.blob), Inches(44.20), Inches(0.78), width=Inches(3.05))


def add_challenge_equation(slide, x: float, y: float, w: float) -> None:
    add_text(
        slide,
        "THE CENTRAL CHALLENGE",
        x,
        y,
        w,
        0.40,
        size=17,
        color=GOLD,
        bold=True,
    )
    add_text(
        slide,
        "Observed expression",
        x,
        y + 0.48,
        2.70,
        0.60,
        size=16,
        color=INK,
        bold=True,
        valign=MSO_ANCHOR.MIDDLE,
    )
    parts = [
        ("Spaceflight biology", ORANGE, 2.75),
        ("Tissue biology", NAVY, 2.20),
        ("Animal variation", MUTED, 2.35),
        ("Study / protocol", GOLD, 2.55),
    ]
    cursor = x + 2.85
    for index, (label, color, width) in enumerate(parts):
        if index:
            add_text(
                slide,
                "+",
                cursor,
                y + 0.48,
                0.38,
                0.60,
                size=21,
                color=MUTED,
                bold=True,
                align=PP_ALIGN.CENTER,
                valign=MSO_ANCHOR.MIDDLE,
            )
            cursor += 0.43
        add_rect(slide, cursor, y + 0.48, width, 0.60, color, rounded=True)
        add_text(
            slide,
            label,
            cursor + 0.07,
            y + 0.50,
            width - 0.14,
            0.54,
            size=14,
            color=WHITE,
            bold=True,
            align=PP_ALIGN.CENTER,
            valign=MSO_ANCHOR.MIDDLE,
        )
        cursor += width


def add_architecture(slide, x: float, y: float, w: float, h: float) -> None:
    add_text(
        slide,
        "Pathway structure is wired into the decoder",
        x,
        y,
        w,
        0.44,
        size=18.5,
        color=INK,
        bold=True,
        align=PP_ALIGN.CENTER,
    )
    add_text(
        slide,
        "Train a tissue reference, then map flight and ground samples into the same annotated program space.",
        x + 0.35,
        y + 0.47,
        w - 0.70,
        0.45,
        size=14.5,
        color=MUTED,
        align=PP_ALIGN.CENTER,
    )

    top_y = y + 1.02
    add_box(
        slide,
        "ARCHS4 reference",
        "tissue-matched non-spaceflight counts",
        x + 0.25,
        top_y,
        4.35,
        0.90,
        fill="E9F3F8",
        line=BLUE,
    )
    add_box(
        slide,
        "Reactome mask",
        "mouse gene-to-program memberships",
        x + 5.15,
        top_y,
        4.35,
        0.90,
        fill="E9F4EE",
        line=GREEN,
    )
    add_box(
        slide,
        "NASA OSDR query",
        "FLT and GC | accession conditioned",
        x + 10.05,
        top_y,
        4.35,
        0.90,
        fill="FAEFE7",
        line=ORANGE,
    )

    network_y = y + 2.25
    input_x = x + 0.85
    encoder_x = x + 2.20
    latent_x = x + 5.45
    output_x = x + 8.20
    score_x = x + 10.30

    add_down_arrow(slide, x + 2.43, top_y + 0.94, BLUE)
    add_down_arrow(slide, x + 7.33, top_y + 0.94, GREEN)
    add_down_arrow(slide, x + 12.23, top_y + 0.94, ORANGE)

    gene_ys = [network_y + 0.28 + 0.55 * index for index in range(5)]
    for gene_y in gene_ys:
        node = slide.shapes.add_shape(
            MSO_SHAPE.OVAL,
            Inches(input_x),
            Inches(gene_y),
            Inches(0.28),
            Inches(0.28),
        )
        set_fill(node, "A9BBD0")
        node.line.fill.background()

    encoder = slide.shapes.add_shape(
        MSO_SHAPE.TRAPEZOID,
        Inches(encoder_x),
        Inches(network_y + 0.10),
        Inches(2.10),
        Inches(2.85),
    )
    set_fill(encoder, "DCE6F2")
    set_line(encoder, BLUE, 1.4)
    add_text(
        slide,
        "dense\nencoder",
        encoder_x + 0.26,
        network_y + 0.90,
        1.55,
        0.92,
        size=15.5,
        color=INK,
        bold=True,
        align=PP_ALIGN.CENTER,
        valign=MSO_ANCHOR.MIDDLE,
    )
    for gene_y in gene_ys:
        add_connector(
            slide,
            input_x + 0.28,
            gene_y + 0.14,
            encoder_x + 0.25,
            network_y + 1.52,
            color="A9BBD0",
            width=0.8,
        )

    latent_rows = [
        ("DNA repair", THYMUS, network_y + 0.35),
        ("T-cell signaling", SPLEEN, network_y + 1.25),
        ("Cell junctions", SKIN, network_y + 2.15),
    ]
    for label, color, node_y in latent_rows:
        add_connector(
            slide,
            encoder_x + 1.85,
            network_y + 1.52,
            latent_x,
            node_y + 0.18,
            color="A9BBD0",
            width=0.9,
        )
        node = slide.shapes.add_shape(
            MSO_SHAPE.OVAL,
            Inches(latent_x),
            Inches(node_y),
            Inches(0.36),
            Inches(0.36),
        )
        set_fill(node, color)
        node.line.fill.background()
        add_text(
            slide,
            label,
            latent_x + 0.48,
            node_y - 0.03,
            1.95,
            0.42,
            size=13.5,
            color=color,
            bold=True,
            valign=MSO_ANCHOR.MIDDLE,
        )

    output_ys = [network_y + 0.20 + 0.65 * index for index in range(5)]
    for out_y in output_ys:
        node = slide.shapes.add_shape(
            MSO_SHAPE.OVAL,
            Inches(output_x),
            Inches(out_y),
            Inches(0.28),
            Inches(0.28),
        )
        set_fill(node, "A9BBD0")
        node.line.fill.background()

    connection_map = [(0, 0), (0, 2), (1, 1), (1, 3), (2, 2), (2, 4)]
    for latent_index, output_index in connection_map:
        _, color, node_y = latent_rows[latent_index]
        add_connector(
            slide,
            latent_x + 0.36,
            node_y + 0.18,
            output_x,
            output_ys[output_index] + 0.14,
            color=color,
            width=1.2,
        )

    add_text(
        slide,
        "genes",
        input_x - 0.20,
        network_y + 3.00,
        0.72,
        0.35,
        size=12.5,
        color=MUTED,
        align=PP_ALIGN.CENTER,
    )
    add_text(
        slide,
        "annotated latent programs",
        latent_x - 0.20,
        network_y + 3.00,
        2.80,
        0.35,
        size=12.5,
        color=MUTED,
        align=PP_ALIGN.CENTER,
    )
    add_text(
        slide,
        "masked decoder",
        output_x - 0.55,
        network_y + 3.00,
        1.50,
        0.35,
        size=12.5,
        color=MUTED,
        align=PP_ALIGN.CENTER,
    )

    add_right_arrow(slide, output_x + 0.55, network_y + 1.25, 1.08, 0.42, BLUE)
    add_box(
        slide,
        "Program scores",
        "posterior mean by sample",
        score_x,
        network_y + 0.55,
        3.75,
        1.02,
        fill=WHITE,
        line=NAVY,
    )
    add_down_arrow(slide, score_x + 1.88, network_y + 1.65, BLUE)
    add_box(
        slide,
        "Project-balanced shift",
        "mean FLT minus GC score",
        score_x,
        network_y + 2.02,
        3.75,
        1.02,
        fill=WHITE,
        line=BLUE,
    )

    add_text(
        slide,
        (
            "approximately 2,000 HVGs | 319-387 retained Reactome programs | "
            "negative-binomial reference | 250-epoch query map"
        ),
        x + 0.25,
        y + h - 0.52,
        w - 0.50,
        0.38,
        size=12.5,
        color=MUTED,
        align=PP_ALIGN.CENTER,
    )


def render_pdf_asset(source: Path, output: Path, dpi: int = 700) -> None:
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


def add_picture_contain(
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


def render_architecture_crop(pdf_path: Path) -> Path:
    output = ASSET_DIR / "expimap_architecture_visualization_300dpi.png"
    dpi = 300
    x, y, w, h = 16.18, 6.28, 15.05, 7.08
    subprocess.run(
        [
            "pdftocairo",
            "-singlefile",
            "-png",
            "-r",
            str(dpi),
            "-x",
            str(round(x * dpi)),
            "-y",
            str(round(y * dpi)),
            "-W",
            str(round(w * dpi)),
            "-H",
            str(round(h * dpi)),
            str(pdf_path),
            str(output.with_suffix("")),
        ],
        check=True,
    )
    return output


def build() -> tuple[Path, Path | None, Path | None, Path | None, list[float]]:
    POSTER_DIR.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    result_asset = ASSET_DIR / "figure_3_tissue_pathway_shifts_700dpi.png"
    hypothesis_asset = ASSET_DIR / "figure_6_tissue_state_hypotheses_700dpi.png"
    render_pdf_asset(FIGURE_DIR / "figure_3_tissue_pathway_shifts.pdf", result_asset)
    render_pdf_asset(
        FIGURE_DIR / "figure_6_tissue_state_hypotheses.pdf", hypothesis_asset
    )

    template = Presentation(TEMPLATE_PATH)
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)
    prs.core_properties.title = TITLE
    prs.core_properties.author = "Jason Trinh"
    prs.core_properties.subject = (
        "Template-based cross-mission expiMap poster using NASA OSDR mouse RNA-seq"
    )
    prs.core_properties.keywords = (
        "spaceflight, expiMap, OSDR, ARCHS4, Reactome, transcriptomics, NASA"
    )
    prs.core_properties.comments = (
        "Layout and NASA branding reference 00 Poster Session Student "
        "Template_Approved by Legal.pptx."
    )
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = rgb(WHITE)

    # Approved-template header
    add_rect(slide, 0, 0, SLIDE_W, 5.05, HEADER)
    add_text(
        slide,
        "National Aeronautics and Space Administration",
        0.80,
        0.28,
        18.0,
        0.48,
        size=15,
        color=MUTED,
    )
    add_text(
        slide,
        DISPLAY_TITLE,
        0.80,
        0.92,
        42.40,
        2.40,
        size=50,
        color="111111",
        bold=True,
        align=PP_ALIGN.CENTER,
        valign=MSO_ANCHOR.MIDDLE,
    )
    add_text(
        slide,
        (
            "Jason Trinh / University of California, Berkeley | "
            "NASA Space Life Sciences Training Program, NASA Ames Research Center"
        ),
        1.50,
        3.70,
        41.0,
        0.72,
        size=24,
        color=MUTED,
        bold=True,
        align=PP_ALIGN.CENTER,
        valign=MSO_ANCHOR.MIDDLE,
    )
    add_nasa_logo(slide, template)

    # Template-style panels
    left_x, left_w = 0.35, 15.15
    center_x, center_w = 15.85, 15.75
    right_x, right_w = 31.95, 15.70

    add_panel(slide, left_x, 5.35, left_w, 3.78)
    add_section(slide, "ABSTRACT", left_x + 0.28, 5.55, left_w - 0.56)
    abstract = (
        "Spaceflight affects multiple organs, but mission differences can obscure "
        "responses that recur across studies. We mapped NASA OSDR mouse bulk "
        "RNA-seq samples into tissue-matched expiMap references trained on ARCHS4 "
        "and constrained by current mouse Reactome programs. Project-balanced "
        "flight-ground shifts were checked with conventional enrichment, held-out "
        "projects, three complete trainings, composition proxies, and member-gene "
        "review. Thymus, skin, liver, and spleen produced reproducible tissue-specific "
        "patterns; spleen had the strongest multi-pathway evidence."
    )
    add_text(
        slide,
        abstract,
        left_x + 0.50,
        6.48,
        left_w - 1.00,
        2.35,
        size=18.5,
        color=BODY,
    )

    add_panel(slide, left_x, 9.43, left_w, 15.25)
    add_section(slide, "INTRODUCTION", left_x + 0.28, 9.63, left_w - 0.56)
    add_text(
        slide,
        "PROJECT OBJECTIVE",
        left_x + 0.55,
        10.62,
        left_w - 1.10,
        0.42,
        size=17,
        color=GOLD,
        bold=True,
    )
    add_text(
        slide,
        (
            "Learn how spaceflight changes living systems by asking which "
            "gene programs shift consistently across missions."
        ),
        left_x + 0.55,
        11.05,
        left_w - 1.10,
        1.10,
        size=23,
        color=INK,
        bold=True,
    )
    add_challenge_equation(slide, left_x + 0.55, 12.30, left_w - 1.10)
    add_text(
        slide,
        (
            "Study identity and protocol can dominate an unconstrained expression "
            "representation. Tissue-matched reference mapping, accession conditioning, "
            "and equal project weighting reduce this bias without claiming to remove "
            "all mission confounding."
        ),
        left_x + 0.55,
        13.62,
        left_w - 1.10,
        1.55,
        size=17,
        color=BODY,
    )
    add_text(
        slide,
        "FINAL ANALYSIS SCOPE",
        left_x + 0.55,
        15.32,
        left_w - 1.10,
        0.40,
        size=17,
        color=GOLD,
        bold=True,
    )
    add_scope_table(
        slide, left_x + 0.55, 15.83, left_w - 1.10, 4.05
    )
    add_text(
        slide,
        (
            "OSDR counts are primary analysis samples. Spleen excludes one "
            "condition-strain-confounded project."
        ),
        left_x + 0.55,
        19.97,
        left_w - 1.10,
        0.55,
        size=13.5,
        color=MUTED,
        italic=True,
    )
    add_text(
        slide,
        "ROBUSTNESS GATE",
        left_x + 0.55,
        20.65,
        left_w - 1.10,
        0.40,
        size=17,
        color=GOLD,
        bold=True,
    )
    badge_y = 21.18
    badge_w = 4.24
    badge_gap = 0.42
    for index, label in enumerate(("ssGSEA + GSEA", "Held-out projects", "3 full trainings")):
        add_badge(
            slide,
            label,
            left_x + 0.55 + index * (badge_w + badge_gap),
            badge_y,
            badge_w,
            color=BLUE if index < 2 else THYMUS,
        )
    for index, label in enumerate(("Composition proxies", "Member genes", "Literature review")):
        add_badge(
            slide,
            label,
            left_x + 0.55 + index * (badge_w + badge_gap),
            badge_y + 0.90,
            badge_w,
            color=GREEN,
        )
    add_text(
        slide,
        (
            "Complementary means an annotated Reactome program that adds a plausible "
            "perspective beyond the dominant phenotype in prior literature. De novo "
            "nodes were not retained in the final models."
        ),
        left_x + 0.55,
        23.08,
        left_w - 1.10,
        1.15,
        size=15.5,
        color=BODY,
    )

    add_panel(slide, center_x, 5.35, center_w, 8.05)
    add_section(slide, "METHODS", center_x + 0.28, 5.55, center_w - 0.56)
    add_architecture(
        slide,
        center_x + 0.40,
        6.42,
        center_w - 0.80,
        6.65,
    )

    add_panel(slide, center_x, 13.70, center_w, 10.98)
    add_section(
        slide,
        "RESULTS: CROSS-MISSION PATHWAY SHIFTS",
        center_x + 0.28,
        13.90,
        center_w - 0.56,
        size=20,
    )
    result_ppi = add_picture_contain(
        slide,
        result_asset,
        center_x + 0.35,
        14.83,
        center_w - 0.70,
        8.82,
        name="Retained pathway shifts",
        description=(
            "Project-balanced expiMap pathway shifts for thymus, skin, liver, and "
            "spleen, with project effects and three-training ranges."
        ),
    )
    add_text(
        slide,
        (
            "Open circles are OSDR projects; colored ranges span three complete "
            "reference-query trainings. All displayed shifts are lower in flight; "
            "axes are tissue specific."
        ),
        center_x + 0.60,
        23.78,
        center_w - 1.20,
        0.55,
        size=13.5,
        color=MUTED,
        italic=True,
        align=PP_ALIGN.CENTER,
    )

    add_panel(slide, right_x, 5.35, right_w, 13.28)
    add_section(
        slide,
        "RESULTS: BIOLOGICAL INTERPRETATION",
        right_x + 0.28,
        5.55,
        right_w - 0.56,
        size=20,
    )
    add_text(
        slide,
        "Final robustness-filtered tissue stories",
        right_x + 0.50,
        6.42,
        right_w - 1.00,
        0.42,
        size=17,
        color=BLUE,
        bold=True,
        align=PP_ALIGN.CENTER,
    )
    hypothesis_ppi = add_picture_contain(
        slide,
        hypothesis_asset,
        right_x + 0.35,
        6.87,
        right_w - 0.70,
        11.22,
        name="Tissue-state hypotheses",
        description=(
            "Observed lower expiMap program scores and qualified tissue-state "
            "hypotheses for thymus, skin, liver, and spleen."
        ),
    )

    add_panel(slide, right_x, 18.93, right_w, 3.55)
    add_section(slide, "CONCLUSIONS", right_x + 0.28, 19.13, right_w - 0.56)
    add_paragraphs(
        slide,
        [
            (
                "\u2022 Spleen was strongest: T-cell receptor, neutrophil "
                "degranulation, and C-type lectin receptor programs were lower "
                "across five projects and three trainings; all had GSEA FDR <0.05."
            ),
            (
                "\u2022 Skin and thymus emphasized lower tissue-maintenance programs; "
                "liver added a lower adaptive-immune axis to heterogeneous metabolism."
            ),
            (
                "\u2022 Bulk latent scores are testable hypotheses, not causal or "
                "cell-type-resolved pathway activity."
            ),
        ],
        right_x + 0.55,
        20.05,
        right_w - 1.10,
        2.25,
        size=15.5,
        color=BODY,
        gap=4,
    )

    add_panel(slide, right_x, 22.75, right_w, 1.93)
    add_section(
        slide,
        "REFERENCES",
        right_x + 0.28,
        22.95,
        right_w - 0.56,
        size=18,
    )
    references = (
        "1 Gebre et al., NAR 2025 (OSDR). 2 Lotfollahi et al., Nat Cell Biol "
        "2023 (expiMap). 3 Lachmann et al., Nat Commun 2018 (ARCHS4). "
        "4 Milacic et al., NAR 2024 (Reactome). 5 Horie et al., Sci Rep 2019. "
        "6 Cope et al., Commun Med 2024. 7 Beheshti et al., Sci Rep 2019. "
        "8 Gridley et al., J Appl Physiol 2009."
    )
    add_text(
        slide,
        references,
        right_x + 0.52,
        23.78,
        right_w - 1.04,
        0.78,
        size=13.5,
        color=BODY,
    )

    # Approved-template footer and acknowledgements
    add_rect(slide, center_x, 24.92, center_w + right_w + 0.05, 2.08, HEADER)
    add_text(
        slide,
        "Acknowledgements",
        center_x + 0.35,
        25.02,
        5.0,
        0.45,
        size=20,
        color=MUTED,
        bold=True,
    )
    add_text(
        slide,
        (
            "Mentor: James Casaletto | NASA Space Life Sciences Training Program | "
            "NASA OSDR, ARCHS4, and Reactome investigators and curation teams"
        ),
        center_x + 0.35,
        25.52,
        20.60,
        0.88,
        size=14.5,
        color=BODY,
    )
    add_text(
        slide,
        "jasontrinh@berkeley.edu\ngithub.com/jasont314/nasa-mouse",
        39.25,
        25.32,
        7.92,
        0.90,
        size=14.5,
        color=NAVY,
        bold=True,
        align=PP_ALIGN.RIGHT,
    )
    add_text(
        slide,
        "www.nasa.gov",
        0.70,
        26.10,
        3.60,
        0.40,
        size=14,
        color=MUTED,
    )

    validate(prs)
    output = POSTER_DIR / "asgsr_expimap_poster.pptx"
    prs.save(output)
    pdf_path, preview_path = render_poster(output)
    architecture_path = render_architecture_crop(pdf_path) if pdf_path else None
    return (
        output,
        pdf_path,
        preview_path,
        architecture_path,
        [result_ppi, hypothesis_ppi],
    )


def run() -> None:
    output, pdf_path, preview_path, architecture_path, ppi = build()
    print(output)
    if pdf_path:
        print(pdf_path)
    if preview_path:
        print(preview_path)
    if architecture_path:
        print(architecture_path)
    print(
        "Embedded data-figure effective resolution: "
        + ", ".join(f"{value:.0f} ppi" for value in ppi)
    )


if __name__ == "__main__":
    run()
