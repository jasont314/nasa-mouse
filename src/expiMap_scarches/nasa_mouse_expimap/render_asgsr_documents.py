"""Render the ASGSR Markdown documents as print-ready HTML and PDF."""

from __future__ import annotations

from pathlib import Path

import markdown


ROOT = Path(__file__).resolve().parents[3]
PAPER_DIR = ROOT / "paper/asgsr_expimap_hvg"

STYLE = """
@page { size: Letter; margin: 0.72in 0.76in 0.75in; }
@page figure-wide { size: Letter landscape; margin: 0.52in 0.60in 0.55in; }
* { box-sizing: border-box; }
body {
  color: #171b1d;
  font-family: "DejaVu Serif", Georgia, serif;
  font-size: 10.2pt;
  line-height: 1.46;
  margin: 0 auto;
  max-width: 8.1in;
}
h1, h2, h3 { color: #111517; font-family: "DejaVu Sans", Arial, sans-serif; }
h1 { font-size: 20pt; line-height: 1.16; margin: 0 0 14pt; }
h2 { border-bottom: 0.7pt solid #aeb5b8; font-size: 13pt; margin: 20pt 0 7pt; padding-bottom: 2pt; }
h3 { font-size: 10.7pt; margin: 14pt 0 4pt; }
p { margin: 0 0 7pt; orphans: 3; widows: 3; }
a { color: #205f89; text-decoration: none; }
table { border-collapse: collapse; font-size: 8.5pt; margin: 8pt 0 12pt; width: 100%; }
th { background: #e8eef0; font-family: "DejaVu Sans", Arial, sans-serif; }
th, td { border: 0.6pt solid #aeb5b8; padding: 4pt 5pt; text-align: left; vertical-align: top; }
code { background: #f1f3f3; font-family: "DejaVu Sans Mono", monospace; font-size: 8.3pt; padding: 0 2pt; }
pre { background: #f1f3f3; border-left: 2pt solid #607d86; overflow-wrap: anywhere; padding: 7pt; white-space: pre-wrap; }
.figure { break-before: page; break-inside: avoid; margin: 0 0 7pt; text-align: center; }
h2 + .figure { break-before: avoid; }
.figure-wide {
  align-items: start;
  column-gap: 0.22in;
  display: grid;
  grid-template-columns: 7.55in 2.03in;
  margin-left: -0.74in;
  page: figure-wide;
  width: 9.80in;
}
.figure img { height: auto; max-height: 8.2in; max-width: 100%; object-fit: contain; }
.figure p { text-align: left; }
.figure-wide p:first-child { grid-column: 1; margin: 0; text-align: center; }
.figure-wide p:last-child { font-size: 8.1pt; grid-column: 2; line-height: 1.30; margin: 0; }
.figure-wide img { max-height: 6.75in; width: 100%; }
ul { margin-top: 4pt; }
"""


def render(source_name: str, title: str) -> Path:
    source = PAPER_DIR / source_name
    body = markdown.markdown(
        source.read_text(encoding="utf-8"),
        extensions=["extra", "md_in_html", "sane_lists", "smarty"],
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{STYLE}</style>
</head>
<body>{body}</body>
</html>
"""
    output = source.with_suffix(".html")
    output.write_text(document, encoding="utf-8")
    return output


def run() -> None:
    outputs = [
        render(
            "asgsr_2026_abstract.md",
            "Cross-Mission expiMap Analysis Recovers Established Tissue Responses "
            "and Identifies Complementary Pathway Shifts in Mouse Spaceflight "
            "Transcriptomes",
        ),
        render(
            "manuscript.md",
            "Cross-Mission expiMap Analysis Recovers Established Tissue Responses "
            "and Identifies Complementary Pathway Shifts in Mouse Spaceflight "
            "Transcriptomes",
        ),
    ]
    try:
        from weasyprint import HTML
    except ImportError:
        for path in outputs:
            print(path)
        print("WeasyPrint unavailable; HTML rendered, PDF skipped.")
        return

    for path in outputs:
        pdf_path = path.with_suffix(".pdf")
        HTML(filename=path).write_pdf(pdf_path)
        print(pdf_path)


if __name__ == "__main__":
    run()
