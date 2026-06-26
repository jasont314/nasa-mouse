"""Prepare tissue-specific NASA OSDR API AnnData inputs for expiMap."""

from __future__ import annotations

import argparse
from io import BytesIO
import json
from pathlib import Path

from .fetch_osdr_mouse_transcriptomics import (
    DEFAULT_OUTPUT_DIR as DEFAULT_OSDR_API_DIR,
    data_query_url,
    discover_metadata,
    download_count_tables,
    read_url_bytes,
    write_metadata_outputs,
)
from .io import require_import


DEFAULT_METADATA = (
    "data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv"
)
DEFAULT_GMT = "data/pathways/reactome_current_mouse_ensembl.gmt"
DEFAULT_LIVER_EXCLUDE = "data/filters/aggregate_liver_12_muscle_candidate_profiles.txt"
PRIMARY_TRANSFORMS = ("raw_counts",)
SENSITIVITY_TRANSFORMS = ("cpm", "log1p_cpm")
UNAVAILABLE_TRANSFORMS = ("tpm", "log1p_tpm")


def read_gmt(path: str | Path) -> list[dict]:
    rows = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term, description, *genes = parts
            rows.append(
                {
                    "term": term,
                    "description": description,
                    "genes": list(dict.fromkeys(g for g in genes if g)),
                }
            )
    return rows


def safe_tissue_name(value: str) -> str:
    return str(value).strip().lower().replace(" ", "_").replace("/", "_")


def load_excluded_profiles(path: str | Path | None) -> set[str]:
    if not path:
        return set()
    path = Path(path)
    if not path.exists():
        return set()
    excluded = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        token = line.strip()
        if token and not token.startswith("#"):
            excluded.add(token)
    return excluded


def load_or_discover_metadata(path: str | Path, timeout: int):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    path = Path(path)
    if path.exists():
        return pd.read_csv(path, sep="\t", keep_default_na=False)
    metadata = discover_metadata(timeout)
    write_metadata_outputs(metadata, path.parent)
    return metadata


def select_profiles(metadata, tissue: str, excluded: set[str], min_per_condition: int):
    tissue = safe_tissue_name(tissue)
    selected = metadata.loc[metadata["tissue_final"].astype(str).eq(tissue)].copy()
    selected = selected.loc[
        selected["condition_inferred"].isin(["flight", "ground_control"])
    ].copy()
    selected = selected.loc[
        selected["bulk_rna_seq_inferred"].map(
            lambda value: value is True or str(value).strip().lower() == "true"
        )
    ].copy()
    if "profile_id" not in selected:
        selected["profile_id"] = (
            selected["id.accession"].astype(str)
            + "/"
            + selected["id.sample name"].astype(str)
        )
    selected = selected.drop_duplicates(
        subset=["id.accession", "id.assay name", "id.sample name"],
        keep="first",
    ).copy()

    if excluded:
        columns = [
            column
            for column in [
                "profile",
                "id.sample name",
                "id.accession_sample name",
            ]
            if column in selected
        ]
        keep = []
        for _, row in selected.iterrows():
            tokens = {str(row.get(column, "")) for column in columns}
            keep.append(not bool(tokens & excluded))
        selected = selected.loc[keep].copy()

    counts = selected["condition_inferred"].value_counts()
    for condition in ["flight", "ground_control"]:
        if int(counts.get(condition, 0)) < min_per_condition:
            raise SystemExit(
                f"{tissue} has {counts.get(condition, 0)} {condition} samples; "
                f"minimum is {min_per_condition}."
            )
    return selected.sort_values(
        ["id.accession", "condition_inferred", "profile"],
        kind="stable",
    ).reset_index(drop=True)


def count_table_path(api_dir: Path, accession: str) -> Path:
    return api_dir / "counts" / f"{accession}_unnormalized_counts.csv"


def ensure_count_tables(
    metadata,
    api_dir: Path,
    timeout: int,
    overwrite: bool,
) -> dict[str, Path]:
    accessions = sorted(metadata["id.accession"].astype(str).unique())
    missing = [
        accession
        for accession in accessions
        if overwrite or not count_table_path(api_dir, accession).exists()
    ]
    if missing:
        download_count_tables(
            metadata,
            api_dir,
            timeout,
            accessions=missing,
            overwrite=overwrite,
        )
    return {accession: count_table_path(api_dir, accession) for accession in accessions}


def sample_name_from_count_column(column: str) -> str:
    return str(column).split("/")[-1]


def load_counts_from_api_tables(selected, count_paths: dict[str, Path]):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    frames = []
    sample_to_profile = {
        (str(row["id.accession"]), str(row["id.sample name"])): str(row["profile_id"])
        for _, row in selected.iterrows()
    }
    for accession, path in count_paths.items():
        if not path.exists():
            raise FileNotFoundError(path)
        table = pd.read_csv(path, keep_default_na=False)
        gene_col = table.columns[0]
        rename = {}
        keep_cols = [gene_col]
        for column in table.columns[1:]:
            sample = sample_name_from_count_column(column)
            key = (str(accession), sample)
            if key in sample_to_profile:
                keep_cols.append(column)
                rename[column] = sample_to_profile[key]
        if len(keep_cols) == 1:
            continue
        sub = table[keep_cols].rename(columns={gene_col: "gene_id", **rename})
        sub = sub.set_index("gene_id")
        frames.append(sub)
    if not frames:
        raise SystemExit("No selected samples were found in downloaded OSDR API counts.")

    matrix = pd.concat(frames, axis=1, join="outer")
    matrix = matrix.loc[:, selected["profile_id"].astype(str).tolist()]
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    duplicated = matrix.index.duplicated(keep="first")
    if duplicated.any():
        matrix = matrix.loc[~duplicated]
    return matrix


def build_architecture(count_genes: list[str], gmt_rows: list[dict], min_genes: int, max_terms: int | None):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    gene_to_index = {gene: index for index, gene in enumerate(count_genes)}
    retained_terms = []
    for row in gmt_rows:
        indices = sorted({gene_to_index[gene] for gene in row["genes"] if gene in gene_to_index})
        if len(indices) >= min_genes:
            retained_terms.append({**row, "indices": indices, "n_genes": len(indices)})
    retained_terms.sort(key=lambda row: (-row["n_genes"], row["term"]))
    if max_terms is not None:
        retained_terms = retained_terms[:max_terms]
    retained_terms.sort(key=lambda row: row["term"])

    retained_gene_indices = sorted(
        {index for row in retained_terms for index in row["indices"]}
    )
    old_to_new = {old: new for new, old in enumerate(retained_gene_indices)}
    mask = np.zeros((len(retained_gene_indices), len(retained_terms)), dtype="float32")
    for term_index, row in enumerate(retained_terms):
        for old_index in row["indices"]:
            new_index = old_to_new.get(old_index)
            if new_index is not None:
                mask[new_index, term_index] = 1.0

    var = pd.DataFrame(
        {
            "gene_id": [count_genes[index] for index in retained_gene_indices],
            "reactome_pathway_count": mask.sum(axis=1).astype(int),
            "source_gene_index": retained_gene_indices,
        }
    )
    terms = pd.DataFrame(
        {
            "term": [row["term"] for row in retained_terms],
            "description": [row["description"] for row in retained_terms],
            "n_genes": [int(mask[:, i].sum()) for i in range(mask.shape[1])],
        }
    )
    return retained_gene_indices, var, terms, mask


def transformed_matrix(counts, transform: str):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    counts = np.asarray(counts, dtype="float32")
    if transform == "raw_counts":
        return counts, {
            "transformation": transform,
            "normalization_method": "none",
            "library_size_handling": "raw OSDR API count table library sizes",
            "recommended_recon_loss": "nb",
            "validity": "primary_count_likelihood",
        }
    library_size = counts.sum(axis=1, keepdims=True)
    library_size[library_size <= 0] = 1.0
    cpm = counts / library_size * 1_000_000.0
    if transform == "cpm":
        return cpm.astype("float32", copy=False), {
            "transformation": transform,
            "normalization_method": "counts_per_million",
            "library_size_handling": "per-sample total-count normalization to 1e6",
            "recommended_recon_loss": "mse",
            "validity": "sensitivity_only_not_count_likelihood",
        }
    if transform == "log1p_cpm":
        return np.log1p(cpm).astype("float32", copy=False), {
            "transformation": transform,
            "normalization_method": "log1p(counts_per_million)",
            "library_size_handling": "per-sample total-count normalization to 1e6",
            "recommended_recon_loss": "mse",
            "validity": "sensitivity_only_not_count_likelihood",
        }
    raise ValueError(f"Unsupported generated transform: {transform}")


def write_counts_tables(selected, output_dir: Path) -> dict[str, str]:
    paths = {}
    counts = (
        selected.groupby(["id.accession", "condition_inferred"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    for column in ["flight", "ground_control"]:
        if column not in counts:
            counts[column] = 0
    counts["total"] = counts["flight"] + counts["ground_control"]
    path = output_dir / "sample_counts_by_accession_condition.tsv"
    counts.to_csv(path, sep="\t")
    paths["sample_counts_by_accession_condition"] = str(path)

    strata_cols = [
        column
        for column in [
            "tissue_final",
            "id.accession",
            "condition_inferred",
            "study.characteristics.sex",
            "study.characteristics.strain",
            "study.characteristics.genotype",
            "investigation.study assays.study assay technology type",
            "investigation.study.comment.data source accession",
            "investigation.study.comment.project type",
        ]
        if column in selected
    ]
    strata = (
        selected.groupby(strata_cols, dropna=False)
        .size()
        .rename("n_samples")
        .reset_index()
        .sort_values("n_samples", ascending=False)
    )
    path = output_dir / "sample_counts_by_metadata.tsv"
    strata.to_csv(path, sep="\t", index=False)
    paths["sample_counts_by_metadata"] = str(path)
    return paths


def run(args) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")

    tissue = safe_tissue_name(args.tissue)
    api_dir = Path(args.api_dir)
    output_dir = Path(args.output_dir or f"outputs/expimap_direct_osdr_{tissue}/input")
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = load_or_discover_metadata(args.metadata, args.timeout)
    exclude_path = args.exclude_profiles_file
    if exclude_path is None and tissue == "liver" and Path(DEFAULT_LIVER_EXCLUDE).exists():
        exclude_path = DEFAULT_LIVER_EXCLUDE
    selected = select_profiles(
        metadata,
        tissue,
        load_excluded_profiles(exclude_path),
        args.min_per_condition,
    )
    count_paths = ensure_count_tables(selected, api_dir, args.timeout, args.overwrite_counts)
    count_matrix = load_counts_from_api_tables(selected, count_paths)

    gmt_rows = read_gmt(args.gmt)
    gene_indices, var, terms, mask = build_architecture(
        count_matrix.index.astype(str).tolist(),
        gmt_rows,
        args.min_genes,
        args.max_terms,
    )
    if not gene_indices or terms.empty:
        raise SystemExit("No genes/pathways retained after Reactome filtering.")

    retained_counts = count_matrix.iloc[gene_indices, :].T.to_numpy(dtype="float32")
    obs = selected.copy()
    obs.index = obs["profile_id"].astype(str)
    var = var.copy()
    var.index = var["gene_id"].astype(str)

    generated = []
    skipped = []
    transforms = args.transform or list(PRIMARY_TRANSFORMS + SENSITIVITY_TRANSFORMS)
    for transform in transforms:
        if transform in UNAVAILABLE_TRANSFORMS:
            skipped.append(
                {
                    "transformation": transform,
                    "reason": (
                        "TPM unavailable from the selected OSDR API unnormalized "
                        "count tables; no transcript-length or TPM field was used."
                    ),
                }
            )
            continue
        x, preprocessing = transformed_matrix(retained_counts, transform)
        adata = ad.AnnData(X=x, obs=obs.copy(), var=var.copy())
        adata.layers["counts"] = retained_counts.copy()
        adata.varm["I"] = mask.copy()
        adata.uns["terms"] = terms["term"].astype(str).tolist()
        adata.uns["term_descriptions"] = terms["description"].astype(str).tolist()
        adata.uns["expimap_preprocessing"] = {
            **preprocessing,
            "gene_filtering": (
                f"genes retained if present in >=1 Reactome pathway after "
                f"min_genes={args.min_genes} pathway filtering"
            ),
            "pathway_architecture": str(args.gmt),
            "min_genes": args.min_genes,
            "n_input_samples": int(adata.n_obs),
            "n_retained_genes": int(adata.n_vars),
            "n_retained_pathways": int(len(adata.uns["terms"])),
            "counts_layer": "NASA OSDR API unnormalized counts",
            "osdr_api_data_query_template": data_query_url("{ACCESSION}"),
        }
        output_h5ad = output_dir / f"osdr_{tissue}_flt_gc_reactome_{transform}.h5ad"
        adata.write_h5ad(output_h5ad)
        generated.append(
            {
                "transformation": transform,
                "h5ad": str(output_h5ad),
                "recommended_recon_loss": preprocessing["recommended_recon_loss"],
                "validity": preprocessing["validity"],
            }
        )

    metadata_path = output_dir / "profile_metadata.tsv"
    selected.to_csv(metadata_path, sep="\t", index=False)
    gene_path = output_dir / "gene_universe.tsv"
    var.reset_index(drop=True).to_csv(gene_path, sep="\t", index=False)
    terms_path = output_dir / "reactome_terms.tsv"
    terms.to_csv(terms_path, sep="\t", index=False)
    count_paths_summary = write_counts_tables(selected, output_dir)

    manifest = {
        "tissue": tissue,
        "source": "NASA OSDR Biological Data API",
        "metadata": str(args.metadata),
        "gmt": str(args.gmt),
        "exclude_profiles_file": str(exclude_path or ""),
        "filters": {
            "organism": "Mus musculus",
            "file.datatype": "unnormalized counts",
            "assay": "bulk RNA-seq inferred from API assay/file text",
            "conditions": ["Space Flight", "Ground Control"],
            "controls_excluded": ["Basal Control", "Vivarium Control", "unknown"],
            "all_data_sources": True,
            "min_genes_per_pathway": args.min_genes,
            "min_per_condition": args.min_per_condition,
            "max_terms": args.max_terms,
        },
        "counts": {
            "samples": int(len(selected)),
            "flight": int((selected["condition_inferred"] == "flight").sum()),
            "ground_control": int(
                (selected["condition_inferred"] == "ground_control").sum()
            ),
            "accessions": int(selected["id.accession"].nunique()),
            "genes": int(len(var)),
            "pathways": int(len(terms)),
        },
        "outputs": {
            "profile_metadata": str(metadata_path),
            "gene_universe": str(gene_path),
            "reactome_terms": str(terms_path),
            **count_paths_summary,
            "h5ad": generated,
            "skipped_transformations": skipped,
        },
    }
    manifest_path = output_dir / "input_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare tissue-specific OSDR API FLT/GC expiMap AnnData inputs."
    )
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--api-dir", default=DEFAULT_OSDR_API_DIR)
    parser.add_argument("--gmt", default=DEFAULT_GMT)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--min-genes", type=int, default=12)
    parser.add_argument("--min-per-condition", type=int, default=3)
    parser.add_argument("--max-terms", type=int, default=None)
    parser.add_argument(
        "--transform",
        action="append",
        choices=list(PRIMARY_TRANSFORMS + SENSITIVITY_TRANSFORMS + UNAVAILABLE_TRANSFORMS),
        help="Transformation to generate. Repeat to generate multiple.",
    )
    parser.add_argument("--exclude-profiles-file", default=None)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--overwrite-counts", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
