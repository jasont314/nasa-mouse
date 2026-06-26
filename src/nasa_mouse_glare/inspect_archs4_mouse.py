"""Inspect ARCHS4 mouse metadata for tissue-specific reference subsets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

from .io import require_import


DEFAULT_ARCHS4 = "assets/archs4/mouse_gene_v2.5.h5"
DEFAULT_OUTPUT_DIR = "data/archs4"
LEAKAGE_TERMS = [
    "NASA",
    "GeneLab",
    "OSD",
    "GLDS",
    "spaceflight",
    "microgravity",
    "ISS",
    "Rodent Research",
    "RR-",
    "hindlimb unloading",
]
TISSUE_KEYWORDS = {
    "liver": ["liver", "hepatic", "hepatocyte"],
    "kidney": ["kidney", "renal", "nephron"],
    "skin": ["skin", "dermal", "epiderm"],
    "skeletal_muscle": ["skeletal muscle", "soleus", "gastrocnemius", "quadriceps"],
    "spleen": ["spleen", "splenic"],
    "lung": ["lung", "pulmonary"],
}


def decode_array(values) -> list[str]:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8", "replace"))
        else:
            decoded.append(str(value))
    return decoded


def clean_token(value: str) -> str:
    return " ".join(str(value).lower().split())


def tissue_match_text(row_text: str, tissue: str) -> bool:
    keywords = TISSUE_KEYWORDS.get(tissue, [tissue.replace("_", " ")])
    text = clean_token(row_text)
    return any(keyword in text for keyword in keywords)


def has_leakage(row_text: str) -> bool:
    text = clean_token(row_text)
    return any(term.lower() in text for term in LEAKAGE_TERMS)


def load_sample_metadata(path: str | Path):
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    columns = [
        "geo_accession",
        "series_id",
        "title",
        "source_name_ch1",
        "characteristics_ch1",
        "library_strategy",
        "library_source",
        "organism_ch1",
        "singlecellprobability",
    ]
    rows = {}
    with h5py.File(path, "r") as handle:
        samples = handle["/meta/samples"]
        for column in columns:
            if column not in samples:
                continue
            values = samples[column][:]
            if values.dtype.kind in {"S", "O", "U"}:
                rows[column] = decode_array(values)
            else:
                rows[column] = values
    metadata = pd.DataFrame(rows)
    metadata.insert(0, "archs4_sample_index", range(len(metadata)))
    return metadata


def classify(metadata):
    metadata = metadata.copy()
    text_cols = [
        column
        for column in [
            "geo_accession",
            "series_id",
            "title",
            "source_name_ch1",
            "characteristics_ch1",
            "library_strategy",
            "library_source",
            "organism_ch1",
        ]
        if column in metadata
    ]
    metadata["filter_text"] = metadata[text_cols].astype(str).agg(" | ".join, axis=1)
    metadata["leakage_excluded"] = metadata["filter_text"].map(has_leakage)
    if "singlecellprobability" in metadata:
        metadata["singlecellprobability"] = metadata["singlecellprobability"].astype(float)
    else:
        metadata["singlecellprobability"] = 0.0
    for tissue in TISSUE_KEYWORDS:
        metadata[f"matches_{tissue}"] = metadata["filter_text"].map(
            lambda text, tissue=tissue: tissue_match_text(text, tissue)
        )
    return metadata


def write_outputs(metadata, output_dir: Path, max_rows: int) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    for tissue in TISSUE_KEYWORDS:
        matched = metadata[metadata[f"matches_{tissue}"]]
        usable = matched[
            ~matched["leakage_excluded"]
            & matched["singlecellprobability"].fillna(0).lt(0.5)
        ]
        summary_rows.append(
            {
                "tissue": tissue,
                "matched_samples": int(len(matched)),
                "usable_nonleakage_bulk_like_samples": int(len(usable)),
                "unique_series": int(usable["series_id"].astype(str).nunique())
                if "series_id" in usable
                else 0,
            }
        )
        sample_path = output_dir / f"archs4_mouse_{tissue}_candidate_samples.tsv"
        usable.head(max_rows).drop(columns=["filter_text"]).to_csv(
            sample_path,
            sep="\t",
            index=False,
        )
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    summary = pd.DataFrame(summary_rows).sort_values(
        "usable_nonleakage_bulk_like_samples",
        ascending=False,
    )
    summary_path = output_dir / "archs4_mouse_tissue_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    manifest = {
        "leakage_terms": LEAKAGE_TERMS,
        "tissue_keywords": TISSUE_KEYWORDS,
        "outputs": {"summary": str(summary_path)},
    }
    manifest_path = output_dir / "archs4_mouse_inspection_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def run(args) -> Path:
    metadata = classify(load_sample_metadata(args.input))
    manifest = write_outputs(metadata, Path(args.output_dir), args.max_rows)
    print(json.dumps(manifest, indent=2))
    return Path(manifest["outputs"]["manifest"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect ARCHS4 mouse metadata for tissue reference subsets."
    )
    parser.add_argument("--input", default=DEFAULT_ARCHS4)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--max-rows", type=int, default=5000)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
