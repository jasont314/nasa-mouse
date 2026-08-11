"""Load NASA OSDR mouse bulk RNA-seq expression through the Biological Data API."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from .fetch_osdr_mouse_transcriptomics import (
    discover_metadata,
    download_count_tables,
)
from .io import write_matrix_bundle


DEFAULT_API_METADATA = (
    "data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv"
)
DEFAULT_COUNTS_DIR = "data/osdr_api/counts"
DEFAULT_TIMEOUT = 180
DEFAULT_GLARE_LIVER_TARGET_MANIFEST = (
    "outputs/glare/multi_tissue_api/liver/aggregate/inputs/"
    "aligned_tms_api.target.manifest.json"
)


def clean_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return token or "profile"


def biological_sample_name(sample: str) -> str:
    """Remove the API suffix used to identify technical sequencing replicates."""
    return re.sub(r"_techrep\d+$", "", str(sample))


def read_api_metadata(
    path: str | Path = DEFAULT_API_METADATA,
    *,
    refresh: bool = False,
    timeout: int = DEFAULT_TIMEOUT,
) -> pd.DataFrame:
    """Read the cached API inventory, refreshing it from NASA when requested."""
    path = Path(path)
    if refresh or not path.exists():
        metadata = discover_metadata(timeout)
        path.parent.mkdir(parents=True, exist_ok=True)
        metadata.to_csv(path, sep="\t", index=False)
    else:
        metadata = pd.read_csv(path, sep="\t", keep_default_na=False)
    return normalize_api_metadata(metadata)


def normalize_api_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Add stable aliases used by the older GLARE analysis scripts."""
    metadata = metadata.copy()
    required = {
        "id.accession",
        "id.sample name",
        "condition_inferred",
        "tissue_final",
    }
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(f"OSDR API metadata is missing columns: {missing}")

    metadata["id.accession"] = metadata["id.accession"].astype(str)
    metadata["id.sample name"] = metadata["id.sample name"].astype(str)
    metadata["profile"] = metadata["id.sample name"].astype(str)
    metadata["accession"] = metadata["id.accession"].astype(str)
    # Keep these aliases while historical GLARE result readers still use them.
    metadata["h5_accession"] = metadata["id.accession"].astype(str)
    metadata["h5_sample_name"] = metadata["id.sample name"].astype(str)
    metadata["h5_accession_sample_name"] = (
        metadata["id.accession"].astype(str)
        + "/"
        + metadata["id.sample name"].astype(str)
    )

    aliases = {
        "project_identifier": "investigation.study.comment.project identifier",
        "project_type": "investigation.study.comment.project type",
        "assay_technology": (
            "investigation.study assays.study assay technology type"
        ),
        "material_type": "study.characteristics.material type",
        "spaceflight_factor": "study.factor value.spaceflight",
        "sex": "study.characteristics.sex",
        "strain": "study.characteristics.strain",
        "genotype": "study.characteristics.genotype",
    }
    for alias, source in aliases.items():
        metadata[alias] = (
            metadata[source].astype(str) if source in metadata.columns else ""
        )
    for column in [
        "source_name",
        "tissue_type",
        "age_at_launch",
        "age",
        "duration",
        "sample_preservation_method",
        "library_selection",
        "library_layout",
        "sequencing_instrument",
    ]:
        if column not in metadata.columns:
            metadata[column] = ""

    metadata["condition_label"] = metadata["condition_inferred"].map(
        {"flight": "FLT", "ground_control": "GC"}
    )
    metadata["condition"] = metadata["condition_inferred"].map(
        {"flight": "flight", "ground_control": "ground"}
    )
    metadata["profile_id"] = (
        metadata["id.accession"].astype(str)
        + "/"
        + metadata["id.sample name"].astype(str)
    )
    return metadata


def select_api_metadata(
    metadata: pd.DataFrame,
    *,
    tissue: str | None = None,
    accessions: list[str] | None = None,
    material_terms: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Select eligible FLT/GC profiles from the API inventory."""
    selected = normalize_api_metadata(metadata)
    mask = selected["condition_inferred"].isin(["flight", "ground_control"])
    if tissue:
        mask &= selected["tissue_final"].astype(str).eq(str(tissue))
    if accessions:
        mask &= selected["id.accession"].isin(set(map(str, accessions)))
    if material_terms:
        material_column = "study.characteristics.material type"
        terms = {str(term).strip().lower() for term in material_terms}
        mask &= selected[material_column].astype(str).str.strip().str.lower().isin(terms)
    selected = selected.loc[mask].copy()
    selected = selected.drop_duplicates(
        subset=["id.accession", "id.sample name"], keep="first"
    )
    return selected.reset_index(drop=True)


def ensure_count_tables(
    metadata: pd.DataFrame,
    counts_dir: str | Path = DEFAULT_COUNTS_DIR,
    *,
    timeout: int = DEFAULT_TIMEOUT,
    download_missing: bool = False,
) -> dict[str, Path]:
    """Resolve one API unnormalized-count table per selected accession."""
    counts_dir = Path(counts_dir)
    paths = {
        accession: counts_dir / f"{accession}_unnormalized_counts.csv"
        for accession in sorted(metadata["id.accession"].astype(str).unique())
    }
    missing = [accession for accession, path in paths.items() if not path.exists()]
    if missing and download_missing:
        if counts_dir.name != "counts":
            raise ValueError(
                "API downloads require --counts-dir to end in 'counts'; "
                f"received {counts_dir}"
            )
        download_count_tables(
            metadata,
            counts_dir.parent,
            timeout,
            accessions=missing,
            overwrite=False,
        )
    missing = [accession for accession, path in paths.items() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing NASA API count CSVs. Re-run with --download-counts: "
            + ", ".join(missing[:20])
        )
    return paths


def read_count_csv(path: str | Path) -> pd.DataFrame:
    """Read one NASA API count response as a genes-by-profiles table."""
    frame = pd.read_csv(path)
    if frame.shape[1] < 2:
        raise ValueError(f"Count table has no sample columns: {path}")
    gene_column = frame.columns[0]
    frame = frame.rename(columns={gene_column: "gene_id"})
    frame["gene_id"] = frame["gene_id"].astype(str)
    value_columns = [column for column in frame.columns if column != "gene_id"]
    frame[value_columns] = frame[value_columns].apply(
        pd.to_numeric, errors="coerce"
    ).fillna(0)
    if frame["gene_id"].duplicated().any():
        frame = frame.groupby("gene_id", as_index=False)[value_columns].sum()
    return frame.set_index("gene_id")


def sample_column_map(columns: list[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for column in columns:
        suffix = str(column).split("/")[-1]
        mapping.setdefault(suffix, str(column))
    return mapping


def log2_cpm(matrix: np.ndarray) -> np.ndarray:
    matrix = matrix.astype(np.float32, copy=False)
    library_sizes = matrix.sum(axis=0, keepdims=True)
    library_sizes[library_sizes <= 0] = 1.0
    cpm = matrix / library_sizes * 1_000_000.0
    return np.log2(cpm + 1.0).astype(np.float32, copy=False)


def load_api_expression(
    selected: pd.DataFrame,
    *,
    counts_dir: str | Path = DEFAULT_COUNTS_DIR,
    timeout: int = DEFAULT_TIMEOUT,
    download_missing: bool = False,
) -> dict[str, object]:
    """Load selected API profiles and collapse only explicit technical replicates."""
    selected = normalize_api_metadata(selected)
    if selected.empty:
        raise ValueError("No OSDR API profiles were selected")
    count_paths = ensure_count_tables(
        selected,
        counts_dir,
        timeout=timeout,
        download_missing=download_missing,
    )
    selected = selected.sort_values(
        ["id.accession", "condition_label", "profile"]
    ).copy()

    blocks: list[pd.DataFrame] = []
    retained_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, str]] = []
    for accession in selected["id.accession"].drop_duplicates().astype(str):
        accession_rows = selected.loc[selected["id.accession"].eq(accession)].copy()
        accession_rows["biological_sample"] = accession_rows["profile"].map(
            biological_sample_name
        )
        table = read_count_csv(count_paths[accession])
        columns = sample_column_map([str(column) for column in table.columns])
        block_columns: list[pd.Series] = []
        for biological_sample, group in accession_rows.groupby(
            "biological_sample", sort=False
        ):
            count_columns = []
            for sample in group["profile"].astype(str):
                column = columns.get(sample)
                if column is None:
                    missing_rows.append(
                        {
                            "id.accession": accession,
                            "id.sample name": sample,
                            "count_csv": str(count_paths[accession]),
                        }
                    )
                else:
                    count_columns.append(column)
            if len(count_columns) != len(group):
                continue
            feature = f"{accession}__{clean_token(biological_sample)}"
            summed = table.loc[:, count_columns].sum(axis=1)
            summed.name = feature
            block_columns.append(summed)

            retained = group.iloc[0].to_dict()
            retained["feature"] = feature
            retained["profile"] = biological_sample
            retained["id.sample name"] = biological_sample
            retained["count_column"] = ";".join(count_columns)
            retained["technical_replicate_count"] = int(len(count_columns))
            retained["technical_replicate_samples"] = ";".join(
                group["profile"].astype(str)
            )
            retained_rows.append(retained)
        if block_columns:
            blocks.append(pd.concat(block_columns, axis=1))

    if not blocks:
        raise ValueError("No API count columns matched the selected metadata")
    common = set(blocks[0].index)
    for block in blocks[1:]:
        common.intersection_update(block.index)
    common_genes = [gene for gene in blocks[0].index if gene in common]
    if not common_genes:
        raise ValueError("No shared genes across selected API count tables")

    raw = pd.concat([block.loc[common_genes] for block in blocks], axis=1)
    if raw.columns.duplicated().any():
        duplicates = raw.columns[raw.columns.duplicated()].astype(str).tolist()
        raise ValueError(f"Duplicate API feature IDs after merging: {duplicates[:10]}")
    raw = raw.astype(np.float32)
    metadata = pd.DataFrame(retained_rows).set_index("feature").loc[raw.columns]
    metadata.index.name = "feature"
    metadata = metadata.reset_index()
    metadata["sample"] = metadata["feature"].astype(str)
    metadata["technical_replicate_group"] = metadata["feature"].astype(str)
    metadata = normalize_api_metadata(metadata)
    metadata["feature"] = raw.columns.astype(str)
    metadata["sample"] = raw.columns.astype(str)
    metadata["stratum"] = "all"

    normalized = pd.DataFrame(
        log2_cpm(raw.to_numpy(dtype=np.float32, copy=False)),
        index=raw.index,
        columns=raw.columns,
    )
    return {
        "raw_counts": raw,
        "log2_cpm": normalized,
        "metadata": metadata,
        "count_paths": count_paths,
        "missing_count_columns": pd.DataFrame(missing_rows),
    }


def write_api_expression_bundles(
    selected: pd.DataFrame,
    output_dir: str | Path,
    *,
    counts_dir: str | Path = DEFAULT_COUNTS_DIR,
    timeout: int = DEFAULT_TIMEOUT,
    download_missing: bool = False,
) -> dict[str, object]:
    """Write GLARE matrix bundles plus raw DESeq2 inputs for an API selection."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    loaded = load_api_expression(
        selected,
        counts_dir=counts_dir,
        timeout=timeout,
        download_missing=download_missing,
    )
    raw = loaded["raw_counts"]
    normalized = loaded["log2_cpm"]
    metadata = loaded["metadata"]
    missing = loaded["missing_count_columns"]
    assert isinstance(raw, pd.DataFrame)
    assert isinstance(normalized, pd.DataFrame)
    assert isinstance(metadata, pd.DataFrame)
    assert isinstance(missing, pd.DataFrame)

    metadata.to_csv(output_dir / "retained_profile_features.tsv", sep="\t", index=False)
    if not missing.empty:
        missing.to_csv(output_dir / "missing_count_columns.tsv", sep="\t", index=False)

    raw_manifest = write_matrix_bundle(
        output_dir / "api_raw_counts",
        raw.to_numpy(dtype=np.float32, copy=False),
        genes=raw.index.astype(str).tolist(),
        profiles=raw.columns.astype(str).tolist(),
        profile_metadata=metadata,
        description="NASA OSDR Biological Data API unnormalized counts",
    )
    normalized_manifest = write_matrix_bundle(
        output_dir / "api_log2_cpm",
        normalized.to_numpy(dtype=np.float32, copy=False),
        genes=normalized.index.astype(str).tolist(),
        profiles=normalized.columns.astype(str).tolist(),
        profile_metadata=metadata,
        description="NASA OSDR Biological Data API log2(CPM+1) expression",
    )

    raw_dir = output_dir / "raw_deseq2_inputs"
    raw_dir.mkdir(parents=True, exist_ok=True)
    counts = raw.copy()
    counts.index.name = "gene_id"
    counts.round().astype(np.int64).to_csv(raw_dir / "counts.tsv", sep="\t")
    metadata.to_csv(raw_dir / "sample_metadata.tsv", sep="\t", index=False)
    pd.DataFrame(
        {"gene_id": raw.index.astype(str), "gene_symbol": raw.index.astype(str)}
    ).to_csv(raw_dir / "gene_symbols.tsv", sep="\t", index=False)

    summary = {
        "raw_manifest": str(raw_manifest),
        "log2_cpm_manifest": str(normalized_manifest),
        "raw_deseq2_inputs": {
            "counts": str(raw_dir / "counts.tsv"),
            "metadata": str(raw_dir / "sample_metadata.tsv"),
            "gene_symbols": str(raw_dir / "gene_symbols.tsv"),
        },
        "genes": int(raw.shape[0]),
        "samples": int(raw.shape[1]),
        "accessions": metadata["id.accession"].drop_duplicates().astype(str).tolist(),
        "missing_count_columns": int(len(missing)),
        "input_kind": "nasa_osdr_api_counts_and_log2_cpm",
    }
    (output_dir / "api_expression_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare NASA OSDR API mouse bulk RNA-seq expression for GLARE."
    )
    parser.add_argument("--metadata", default=DEFAULT_API_METADATA)
    parser.add_argument("--counts-dir", default=DEFAULT_COUNTS_DIR)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--tissue")
    parser.add_argument("--accessions", nargs="+")
    parser.add_argument("--material-terms", nargs="+")
    parser.add_argument("--refresh-metadata", action="store_true")
    parser.add_argument("--download-counts", action="store_true")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metadata = read_api_metadata(
        args.metadata,
        refresh=args.refresh_metadata,
        timeout=args.timeout,
    )
    selected = select_api_metadata(
        metadata,
        tissue=args.tissue,
        accessions=args.accessions,
        material_terms=args.material_terms,
    )
    summary = write_api_expression_bundles(
        selected,
        args.output_dir,
        counts_dir=args.counts_dir,
        timeout=args.timeout,
        download_missing=args.download_counts,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
