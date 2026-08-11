"""Verify manuscript DOI metadata against Crossref and write an audit table."""

from __future__ import annotations

from datetime import datetime, timezone
import difflib
import json
from pathlib import Path
import re
import urllib.parse
import urllib.request

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PAPER_DIR = ROOT / "paper/asgsr_expimap_hvg"
MANUSCRIPT = PAPER_DIR / "manuscript.md"
OUTPUT = PAPER_DIR / "source_data/source_verification.tsv"


def normalize(value: str) -> str:
    value = re.sub(r"<[^>]+>", " ", value)
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def clean_registry_text(value: str) -> str:
    value = re.sub(r"<[^>]+>", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def references() -> list[tuple[int, str, str]]:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    block = text.split("## References", 1)[1].split("## Figure captions", 1)[0]
    parsed = []
    for number, entry in re.findall(r"(?ms)^([0-9]+)\. (.*?)(?=^[0-9]+\. |\Z)", block):
        doi_match = re.search(r"doi:(10\.\S+?)\.$", entry.strip())
        if not doi_match:
            raise RuntimeError(f"Reference {number} has no terminal DOI")
        prefix = entry.split(". *", 1)[0]
        manuscript_title = prefix.split(". ", 1)[1]
        parsed.append((int(number), manuscript_title, doi_match.group(1)))
    return parsed


def crossref(doi: str) -> dict:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "nasa-mouse-source-audit/1.0 (mailto:noreply@example.com)"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)["message"]


def run() -> None:
    rows = []
    failures = []
    verified_at = datetime.now(timezone.utc).isoformat()
    for number, manuscript_title, doi in references():
        try:
            record = crossref(doi)
            registry_title = clean_registry_text((record.get("title") or [""])[0])
            similarity = difflib.SequenceMatcher(
                None, normalize(manuscript_title), normalize(registry_title)
            ).ratio()
            status = "verified" if similarity >= 0.95 else "manual_review"
            journal = clean_registry_text((record.get("container-title") or [""])[0])
            date_parts = (
                record.get("published-print", record.get("published-online", record.get("published", {})))
                .get("date-parts", [[None]])[0]
            )
            year = date_parts[0]
        except Exception as error:  # Network and registry errors must remain visible.
            registry_title = ""
            similarity = 0.0
            status = "error"
            journal = ""
            year = None
            failures.append(f"Reference {number}: {error}")
        rows.append(
            {
                "reference_number": number,
                "status": status,
                "title_similarity": similarity,
                "doi": doi,
                "doi_url": f"https://doi.org/{doi}",
                "manuscript_title": manuscript_title,
                "registry_title": registry_title,
                "registry_journal": journal,
                "registry_year": year,
                "verified_at_utc": verified_at,
            }
        )
    frame = pd.DataFrame(rows)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUTPUT, sep="\t", index=False)
    if failures or not frame["status"].eq("verified").all():
        raise RuntimeError("; ".join(failures) or "One or more titles require manual review")
    print(f"Verified {len(frame)} DOI records: {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    run()
