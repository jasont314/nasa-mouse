"""Build an API-derived, ARCHS4-aligned OSDR raw-count matrix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from urllib.parse import quote
from urllib.request import Request, urlopen

import h5py
import numpy as np
import pandas as pd

from nasa_mouse_glare.fetch_osdr_mouse_transcriptomics import download_count_tables


ENSEMBL_MOUSE = re.compile(r"^ENSMUSG\d+$")
OSDR_API_ROOT = "https://visualization.osdr.nasa.gov/biodata/api/v2"
USER_AGENT = "nasa-mouse/1.0"


def biological_sample_name(value: object) -> str:
    return re.sub(r"_techrep\d+$", "", str(value), flags=re.IGNORECASE)


def normalize_gene_ids(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace(r"\.\d+$", "", regex=True)


def load_archs4_genes(path: str | Path) -> list[str]:
    with h5py.File(path, "r") as handle:
        values = handle["meta/genes/ensembl_gene"][:]
    genes = [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    ]
    return list(dict.fromkeys(genes))


def count_table_path(root: str | Path, accession: str) -> Path:
    return Path(root) / "counts" / f"{accession}_unnormalized_counts.csv"


def ensure_count_tables(
    metadata: pd.DataFrame,
    *,
    api_dir: str | Path,
    download_missing: bool,
    timeout: int,
) -> tuple[dict[str, Path], list[dict]]:
    accessions = sorted(metadata["id.accession"].astype(str).unique())
    paths = {accession: count_table_path(api_dir, accession) for accession in accessions}
    missing = [accession for accession, path in paths.items() if not path.exists()]
    fallback_manifest = Path(api_dir) / "counts" / "rest_download_fallbacks.json"
    if fallback_manifest.exists():
        fallbacks = json.loads(fallback_manifest.read_text(encoding="utf-8"))
        if not isinstance(fallbacks, list):
            raise ValueError(f"Invalid fallback manifest: {fallback_manifest}")
    else:
        fallbacks = []
    if missing and download_missing:
        for accession in missing:
            subset = metadata.loc[
                metadata["id.accession"].astype(str).eq(accession)
            ].copy()
            try:
                download_count_tables(
                    subset,
                    Path(api_dir),
                    timeout,
                    accessions=[accession],
                    overwrite=False,
                )
            except Exception as query_error:
                filenames = subset.get("file.filename", pd.Series(dtype=str))
                filenames = filenames.dropna().astype(str).unique().tolist()
                if len(filenames) != 1:
                    raise RuntimeError(
                        f"Could not identify one count file for {accession}: {filenames}"
                    ) from query_error
                fallback = download_via_rest_file_record(
                    accession,
                    filenames[0],
                    paths[accession],
                    timeout=timeout,
                )
                fallback["query_error"] = repr(query_error)
                fallbacks.append(fallback)
        fallback_manifest.write_text(
            json.dumps(fallbacks, indent=2) + "\n", encoding="utf-8"
        )
    missing = [accession for accession, path in paths.items() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing API count tables; use --download-missing: " + ", ".join(missing)
        )
    return paths, fallbacks


def download_via_rest_file_record(
    accession: str,
    filename: str,
    output: Path,
    *,
    timeout: int,
) -> dict:
    listing_url = f"{OSDR_API_ROOT}/dataset/{quote(accession, safe='-')}/files/"
    request = Request(listing_url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=timeout) as response:
        listing = json.load(response)
    try:
        file_record = listing[accession]["files"][filename]
        download_url = str(file_record["URL"])
    except (KeyError, TypeError) as error:
        raise RuntimeError(
            f"OSDR REST file listing did not contain {accession}/{filename}"
        ) from error
    request = Request(download_url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=timeout) as response:
        payload = response.read()
    if not payload or payload.lstrip().startswith(b"<"):
        raise RuntimeError(f"OSDR REST download did not return a count CSV: {download_url}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)
    return {
        "id.accession": accession,
        "filename": filename,
        "rest_listing_url": listing_url,
        "download_url": download_url,
        "output": str(output),
    }


def _count_column_map(columns: list[object]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for column in columns:
        value = str(column)
        mapping.setdefault(value.split("/")[-1], value)
    return mapping


def _read_accession_block(
    path: Path,
    rows: pd.DataFrame,
    *,
    technical_replicate_policy: str,
) -> tuple[pd.DataFrame, list[dict], list[dict]]:
    table = pd.read_csv(path, low_memory=False)
    gene_column = table.columns[0]
    table = table.rename(columns={gene_column: "gene_id"})
    table["gene_id"] = normalize_gene_ids(table["gene_id"])
    table = table.loc[table["gene_id"].map(lambda value: bool(ENSEMBL_MOUSE.match(value)))]
    value_columns = [column for column in table.columns if column != "gene_id"]
    table[value_columns] = table[value_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    if table["gene_id"].duplicated().any():
        table = table.groupby("gene_id", as_index=False, sort=False)[value_columns].sum()
    table = table.set_index("gene_id")
    column_map = _count_column_map(list(table.columns))

    selected_columns: list[pd.Series] = []
    retained: list[dict] = []
    missing: list[dict] = []
    rows = rows.copy()
    if technical_replicate_policy == "keep":
        rows["biological_sample_name"] = rows["id.sample name"].astype(str)
    else:
        rows["biological_sample_name"] = rows["id.sample name"].map(
            biological_sample_name
        )

    for biological_sample, group in rows.groupby("biological_sample_name", sort=False):
        for column in ("condition_inferred", "tissue_canonical"):
            if column in group and group[column].astype(str).nunique() != 1:
                raise ValueError(
                    f"Technical-replicate group {biological_sample!r} spans multiple {column} values"
                )
        sample_names = group["id.sample name"].astype(str).tolist()
        columns = [column_map.get(sample) for sample in sample_names]
        absent = [sample for sample, column in zip(sample_names, columns) if column is None]
        if absent:
            missing.extend(
                {
                    "id.accession": str(group.iloc[0]["id.accession"]),
                    "id.sample name": sample,
                    "count_table": str(path),
                }
                for sample in absent
            )
            continue
        concrete_columns = [str(column) for column in columns]
        values = table[concrete_columns]
        if technical_replicate_policy == "mean":
            expression = values.mean(axis=1)
        else:
            expression = values.sum(axis=1)
        accession = str(group.iloc[0]["id.accession"])
        profile_id = f"{accession}/{biological_sample}"
        expression.name = profile_id
        selected_columns.append(expression)

        record = group.iloc[0].to_dict()
        record["profile_id"] = profile_id
        record["biological_sample_name"] = str(biological_sample)
        record["api_profile_count"] = int(len(group))
        record["api_sample_names"] = ";".join(sample_names)
        record["count_columns"] = ";".join(concrete_columns)
        retained.append(record)

    if not selected_columns:
        return pd.DataFrame(), retained, missing
    return pd.concat(selected_columns, axis=1), retained, missing


def build_expression(
    metadata: pd.DataFrame,
    count_paths: dict[str, Path],
    *,
    archs4_genes: list[str],
    technical_replicate_policy: str,
) -> tuple[np.ndarray, pd.DataFrame, list[str], pd.DataFrame]:
    blocks: list[pd.DataFrame] = []
    records: list[dict] = []
    missing: list[dict] = []
    for accession in sorted(count_paths):
        rows = metadata.loc[metadata["id.accession"].astype(str).eq(accession)].copy()
        block, retained, absent = _read_accession_block(
            count_paths[accession],
            rows,
            technical_replicate_policy=technical_replicate_policy,
        )
        if not block.empty:
            blocks.append(block)
            records.extend(retained)
        missing.extend(absent)
    if not blocks:
        raise ValueError("No OSDR API count columns matched the selected metadata")

    common = set(blocks[0].index.astype(str))
    for block in blocks[1:]:
        common.intersection_update(block.index.astype(str))
    genes = [gene for gene in archs4_genes if gene in common]
    if not genes:
        raise ValueError("No shared mouse Ensembl genes across OSDR and ARCHS4")
    aligned = pd.concat([block.loc[genes] for block in blocks], axis=1)
    obs = pd.DataFrame(records).set_index("profile_id", drop=False)
    obs = obs.loc[aligned.columns.astype(str)].copy()
    expression = aligned.T.to_numpy(dtype=np.float32, copy=True)
    return expression, obs, genes, pd.DataFrame(missing)


def run(args: argparse.Namespace) -> Path:
    ad = __import__("anndata")
    metadata = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    if args.tissue:
        metadata = metadata.loc[
            metadata["tissue_canonical"].astype(str).isin(set(args.tissue))
        ].copy()
    if args.accession:
        metadata = metadata.loc[
            metadata["id.accession"].astype(str).isin(set(args.accession))
        ].copy()
    if args.exclude_accession:
        metadata = metadata.loc[
            ~metadata["id.accession"].astype(str).isin(set(args.exclude_accession))
        ].copy()
    if metadata.empty:
        raise ValueError("No OSDR profiles remain after filtering")

    count_paths, download_fallbacks = ensure_count_tables(
        metadata,
        api_dir=args.api_dir,
        download_missing=args.download_missing,
        timeout=args.timeout,
    )
    expression, obs, genes, missing = build_expression(
        metadata,
        count_paths,
        archs4_genes=load_archs4_genes(args.archs4_h5),
        technical_replicate_policy=args.technical_replicate_policy,
    )
    var = pd.DataFrame({"gene_id": genes}, index=pd.Index(genes, name="gene_id"))
    for column in obs.select_dtypes(include=["object"]).columns:
        obs[column] = obs[column].fillna("").astype(str)
    adata = ad.AnnData(X=expression, obs=obs, var=var)
    adata.uns["source"] = "NASA OSDR Biological Data API unnormalized count tables"
    adata.uns["technical_replicate_policy"] = args.technical_replicate_policy
    adata.uns["count_semantics"] = (
        "RSEM unnormalized expected counts; non-negative values may be fractional"
    )
    adata.uns["gene_universe"] = "intersection across selected OSDR accessions and ARCHS4 mouse Ensembl genes"

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output, compression="gzip")
    missing_path = output.with_name("missing_api_count_columns.tsv")
    missing.to_csv(missing_path, sep="\t", index=False)
    summary = {
        "source": "NASA OSDR Biological Data API",
        "raw_integrated_osdr_h5_used": False,
        "metadata": str(args.metadata),
        "api_count_directory": str(Path(args.api_dir) / "counts"),
        "technical_replicate_policy": args.technical_replicate_policy,
        "count_semantics": (
            "RSEM unnormalized expected counts; non-negative values may be fractional"
        ),
        "api_download_fallbacks": download_fallbacks,
        "counts": {
            "api_metadata_profiles_selected": int(len(metadata)),
            "training_profiles": int(adata.n_obs),
            "accessions": int(adata.obs["id.accession"].astype(str).nunique()),
            "tissues": int(adata.obs["tissue_canonical"].astype(str).nunique()),
            "genes": int(adata.n_vars),
            "missing_count_columns": int(len(missing)),
        },
        "outputs": {
            "h5ad": str(output),
            "missing_count_columns": str(missing_path),
        },
    }
    summary_path = output.with_name("osdr_api_expression_summary.json")
    summary["outputs"]["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata",
        default="outputs/generative/benchmark/data_audit/osdr/osdr_canonical_metadata.tsv",
    )
    parser.add_argument("--api-dir", default="data/osdr_api")
    parser.add_argument("--archs4-h5", default="assets/archs4/mouse_gene_v2.5.h5")
    parser.add_argument(
        "--output",
        default="outputs/generative/benchmark/data/osdr/osdr_api_raw_counts.h5ad",
    )
    parser.add_argument(
        "--technical-replicate-policy", choices=["keep", "sum", "mean"], default="sum"
    )
    parser.add_argument("--tissue", action="append")
    parser.add_argument("--accession", action="append")
    parser.add_argument("--exclude-accession", action="append")
    parser.add_argument("--download-missing", action="store_true")
    parser.add_argument("--timeout", type=int, default=180)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
