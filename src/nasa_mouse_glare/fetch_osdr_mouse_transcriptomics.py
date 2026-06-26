"""Discover and download mouse bulk OSDR transcriptomics via NASA OSDR API."""

from __future__ import annotations

import argparse
from io import BytesIO
import json
from pathlib import Path
import re
from urllib.parse import quote
from urllib.request import Request, urlopen

from .cluster_stratified_analysis import infer_tissue, tissue_inference_rule
from .io import require_import
from .osdr_tissues import canonical_tissue


API_ROOT = "https://visualization.osdr.nasa.gov/biodata/api/v2"
DEFAULT_OUTPUT_DIR = "data/osdr_api"
USER_AGENT = "nasa-mouse/1.0"

METADATA_FIELDS = [
    "id.accession",
    "id.assay name",
    "id.sample name",
    "investigation.study assays.study assay technology type",
    "investigation.study.comment.data source accession",
    "investigation.study.comment.project identifier",
    "investigation.study.comment.project type",
    "study.characteristics.material type",
    "study.characteristics.organism",
    "study.characteristics.sex",
    "study.characteristics.strain",
    "study.characteristics.genotype",
    "study.factor value.spaceflight",
    "file.datatype",
    "file.filename",
]


def api_url(endpoint: str, params: list[str]) -> str:
    return f"{API_ROOT}{endpoint}/?{'&'.join(params)}"


def read_url_bytes(url: str, timeout: int) -> bytes:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=timeout) as response:
        return response.read()


def read_api_csv(url: str, timeout: int):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    payload = read_url_bytes(url, timeout)
    if payload.lstrip().startswith(b"<"):
        raise RuntimeError(f"API returned HTML instead of CSV for {url}")
    return pd.read_csv(BytesIO(payload), keep_default_na=False)


def metadata_query_url() -> str:
    params = [
        "study.characteristics.organism=mus%20musculus",
        "study.factor%20value.spaceflight=Space%20Flight|Ground%20Control",
        "file.datatype=unnormalized%20counts",
    ]
    for field in METADATA_FIELDS:
        params.append(quote(field, safe="."))
    params.append("format=csv")
    return api_url("/query/metadata", params)


def normalize_condition(value: str) -> str:
    normalized = " ".join(str(value).strip().lower().split())
    if normalized == "space flight":
        return "flight"
    if normalized == "ground control":
        return "ground_control"
    return "excluded_control_or_unknown"


def normalize_tissue(material_type: str, sample_name: str) -> tuple[str, str, str]:
    official = canonical_tissue(material_type)
    if official not in {"unknown", "unspecified", "cells", "cultured_cells"}:
        return official, "osdr_api_material_type", f"OSDR material type: {material_type}"
    inferred = infer_tissue(sample_name)
    if inferred != "unknown":
        return inferred, "sample_name_inference", tissue_inference_rule(sample_name)
    return official, "unassigned", ""


def is_bulk_rna_seq(row) -> bool:
    assay = str(row.get("id.assay name", "")).lower()
    technology = str(
        row.get("investigation.study assays.study assay technology type", "")
    ).lower()
    filename = str(row.get("file.filename", "")).lower()
    text = " ".join([assay, technology, filename])
    if "single-cell" in text or "single cell" in text or "scrna" in text:
        return False
    if "microarray" in text:
        return False
    return "rna" in text and "unnormalized" in filename


def discover_metadata(timeout: int):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    url = metadata_query_url()
    metadata = read_api_csv(url, timeout)
    metadata["api_metadata_query_url"] = url
    metadata["profile"] = metadata["id.sample name"].astype(str)
    metadata["condition_inferred"] = metadata[
        "study.factor value.spaceflight"
    ].map(normalize_condition)
    tissue_records = [
        normalize_tissue(
            row.get("study.characteristics.material type", ""),
            row.get("id.sample name", ""),
        )
        for _, row in metadata.iterrows()
    ]
    metadata["tissue_final"] = [record[0] for record in tissue_records]
    metadata["tissue_source"] = [record[1] for record in tissue_records]
    metadata["tissue_inference_rule"] = [record[2] for record in tissue_records]
    metadata["bulk_rna_seq_inferred"] = metadata.apply(is_bulk_rna_seq, axis=1)
    metadata = metadata.loc[
        metadata["condition_inferred"].isin(["flight", "ground_control"])
        & metadata["bulk_rna_seq_inferred"]
    ].copy()
    metadata = metadata.drop_duplicates(
        subset=["id.accession", "id.assay name", "id.sample name"],
        keep="first",
    ).reset_index(drop=True)
    metadata["profile_id"] = (
        metadata["id.accession"].astype(str)
        + "/"
        + metadata["id.sample name"].astype(str)
    )
    return metadata


def write_metadata_outputs(metadata, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = output_dir / "osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)

    tissue_counts = (
        metadata.groupby(["tissue_final", "condition_inferred"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    for condition in ["flight", "ground_control"]:
        if condition not in tissue_counts:
            tissue_counts[condition] = 0
    tissue_counts["total_flt_gc"] = (
        tissue_counts["flight"] + tissue_counts["ground_control"]
    )
    tissue_counts = tissue_counts.sort_values(
        ["total_flt_gc", "flight", "ground_control"],
        ascending=False,
    )
    tissue_counts_path = output_dir / "osdr_api_mouse_bulk_rnaseq_tissue_counts.tsv"
    tissue_counts.to_csv(tissue_counts_path, sep="\t")

    accession_counts = (
        metadata.groupby(["tissue_final", "id.accession", "condition_inferred"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    for condition in ["flight", "ground_control"]:
        if condition not in accession_counts:
            accession_counts[condition] = 0
    accession_counts["total_flt_gc"] = (
        accession_counts["flight"] + accession_counts["ground_control"]
    )
    accession_counts = accession_counts.sort_values(
        ["tissue_final", "total_flt_gc"],
        ascending=[True, False],
    )
    accession_counts_path = (
        output_dir / "osdr_api_mouse_bulk_rnaseq_tissue_accession_counts.tsv"
    )
    accession_counts.to_csv(accession_counts_path, sep="\t")

    file_counts = (
        metadata.groupby(
            ["id.accession", "id.assay name", "file.filename"],
            dropna=False,
        )
        .size()
        .rename("n_selected_samples")
        .reset_index()
        .sort_values(["id.accession", "file.filename"])
    )
    file_counts_path = output_dir / "osdr_api_mouse_bulk_rnaseq_files.tsv"
    file_counts.to_csv(file_counts_path, sep="\t", index=False)

    summary = {
        "api_root": API_ROOT,
        "metadata_query_url": metadata["api_metadata_query_url"].iloc[0]
        if len(metadata)
        else metadata_query_url(),
        "filters": {
            "organism": "Mus musculus",
            "spaceflight": ["Space Flight", "Ground Control"],
            "file.datatype": "unnormalized counts",
            "bulk_rna_seq_inferred": (
                "include RNA-seq/count files, exclude single-cell and microarray "
                "by assay/file text"
            ),
            "data_sources": "all OSDR sources returned by the API",
        },
        "counts": {
            "samples": int(len(metadata)),
            "accessions": int(metadata["id.accession"].nunique()),
            "files": int(file_counts.shape[0]),
            "tissues": int(metadata["tissue_final"].nunique()),
        },
        "outputs": {
            "metadata": str(metadata_path),
            "tissue_counts": str(tissue_counts_path),
            "tissue_accession_counts": str(accession_counts_path),
            "files": str(file_counts_path),
        },
    }
    summary_path = output_dir / "osdr_api_mouse_bulk_rnaseq_summary.json"
    summary["outputs"]["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def data_query_url(accession: str) -> str:
    params = [
        f"id.accession={quote(accession, safe='-')}",
        "file.datatype=unnormalized%20counts",
        "column.*",
        "format=csv",
    ]
    return api_url("/query/data", params)


def download_count_tables(
    metadata,
    output_dir: Path,
    timeout: int,
    accessions: list[str] | None = None,
    overwrite: bool = False,
) -> dict[str, str]:
    counts_dir = output_dir / "counts"
    counts_dir.mkdir(parents=True, exist_ok=True)
    selected_accessions = sorted(metadata["id.accession"].dropna().astype(str).unique())
    if accessions:
        wanted = set(accessions)
        selected_accessions = [a for a in selected_accessions if a in wanted]
    outputs = {}
    for index, accession in enumerate(selected_accessions, start=1):
        path = counts_dir / f"{accession}_unnormalized_counts.csv"
        if path.exists() and not overwrite:
            outputs[accession] = str(path)
            continue
        url = data_query_url(accession)
        print(f"download {index}/{len(selected_accessions)} {accession}", flush=True)
        payload = read_url_bytes(url, timeout)
        if payload.lstrip().startswith(b"<"):
            raise RuntimeError(f"API returned HTML instead of CSV for {accession}: {url}")
        path.write_bytes(payload)
        outputs[accession] = str(path)
    manifest = {
        "api_root": API_ROOT,
        "accessions": selected_accessions,
        "outputs": outputs,
        "query_url_template": data_query_url("{ACCESSION}"),
    }
    manifest_path = counts_dir / "download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    outputs["_manifest"] = str(manifest_path)
    return outputs


def parse_accessions(values: list[str] | None) -> list[str] | None:
    if not values:
        return None
    accessions = []
    for value in values:
        accessions.extend(token for token in re.split(r"[,\s]+", value) if token)
    return sorted(set(accessions))


def run(args) -> Path:
    output_dir = Path(args.output_dir)
    metadata = discover_metadata(args.timeout)
    summary = write_metadata_outputs(metadata, output_dir)
    if args.download_counts:
        outputs = download_count_tables(
            metadata,
            output_dir,
            args.timeout,
            accessions=parse_accessions(args.accession),
            overwrite=args.overwrite,
        )
        summary["outputs"]["count_tables"] = outputs
        Path(summary["outputs"]["summary"]).write_text(
            json.dumps(summary, indent=2) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(summary, indent=2))
    return Path(summary["outputs"]["metadata"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Discover Mus musculus bulk RNA-seq Space Flight/Ground Control "
            "samples and optionally download unnormalized count tables through "
            "the NASA OSDR Biological Data API."
        )
    )
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--download-counts", action="store_true")
    parser.add_argument(
        "--accession",
        action="append",
        help="Restrict count downloads to accession(s). Repeat or comma-separate.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
