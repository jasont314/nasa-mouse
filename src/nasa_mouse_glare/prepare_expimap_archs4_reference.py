"""Prepare tissue-filtered ARCHS4 reference AnnData for expiMap."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .inspect_archs4_mouse import (
    DEFAULT_ARCHS4,
    classify,
    load_sample_metadata,
    TISSUE_KEYWORDS,
)
from .io import require_import


def decode_array(values) -> list[str]:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8", "replace"))
        else:
            decoded.append(str(value))
    return decoded


def stratified_series_sample(selected, max_samples: int, seed: int):
    """Sample proportional to ARCHS4 series while retaining every series when possible."""
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    groups = list(selected.groupby("series_id", dropna=False, sort=True))
    sizes = np.asarray([len(group) for _, group in groups], dtype=int)
    if max_samples >= int(sizes.sum()):
        return selected.copy(), "all_eligible"
    if max_samples < len(groups):
        raise SystemExit(
            f"max_samples={max_samples} is smaller than the {len(groups)} ARCHS4 series."
        )

    quotas = np.ones(len(groups), dtype=int)
    remaining = max_samples - int(quotas.sum())
    capacity = sizes - quotas
    allocation = capacity / capacity.sum() * remaining
    quotas += np.floor(allocation).astype(int)
    remainder = max_samples - int(quotas.sum())
    order = np.argsort(-(allocation - np.floor(allocation)), kind="stable")
    for index in order:
        if remainder == 0:
            break
        if quotas[index] < sizes[index]:
            quotas[index] += 1
            remainder -= 1

    sampled = [
        group.sample(n=int(quota), random_state=seed + index)
        for index, (_, group), quota in zip(range(len(groups)), groups, quotas)
    ]
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    return pd.concat(sampled, ignore_index=False), "proportional_series_stratified"


def select_archs4_samples(
    metadata,
    tissue: str,
    max_samples: int,
    sample_seed: int,
):
    match_col = f"matches_{tissue}"
    if match_col not in metadata:
        raise SystemExit(f"Unsupported tissue {tissue!r}; known: {sorted(TISSUE_KEYWORDS)}")
    selected = metadata[
        metadata[match_col]
        & ~metadata["leakage_excluded"]
        & metadata["singlecellprobability"].fillna(0).lt(0.5)
    ].copy()
    selected = selected.sort_values(["series_id", "geo_accession"], kind="stable")
    eligible_samples = int(len(selected))
    eligible_series = int(selected["series_id"].astype(str).nunique())
    if max_samples:
        selected, sampling_method = stratified_series_sample(
            selected,
            max_samples=max_samples,
            seed=sample_seed,
        )
    else:
        sampling_method = "all_eligible"
    if selected.empty:
        raise SystemExit(f"No ARCHS4 samples selected for tissue {tissue}.")
    selected = selected.sort_values(["series_id", "geo_accession"], kind="stable")
    selected["archs4_condition"] = selected["series_id"].astype(str)
    selected["profile_id"] = selected["geo_accession"].astype(str)
    selected.index = selected["profile_id"]
    selection = {
        "method": sampling_method,
        "seed": sample_seed,
        "eligible_samples": eligible_samples,
        "eligible_series": eligible_series,
        "selected_samples": int(len(selected)),
        "selected_series": int(selected["series_id"].astype(str).nunique()),
    }
    return selected, selection


def load_reference_counts(archs4_h5: str | Path, selected, query_adata):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")

    query_genes = query_adata.var_names.astype(str).tolist()
    with h5py.File(archs4_h5, "r") as handle:
        archs4_genes = decode_array(handle["/meta/genes/ensembl_gene"][:])
        gene_to_index = {gene: index for index, gene in enumerate(archs4_genes)}
        gene_indices = [gene_to_index[gene] for gene in query_genes if gene in gene_to_index]
        retained_genes = [gene for gene in query_genes if gene in gene_to_index]
        if retained_genes != query_genes:
            query_positions = [query_genes.index(gene) for gene in retained_genes]
        else:
            query_positions = list(range(len(query_genes)))
        sample_indices = selected["archs4_sample_index"].astype(int).tolist()
        sorted_pairs = sorted(enumerate(sample_indices), key=lambda item: item[1])
        expression = handle["/data/expression"]
        # ARCHS4 is chunked by individual sample columns. Streaming columns avoids
        # h5py's memory-heavy multi-column fancy-index implementation.
        counts = np.empty((len(sample_indices), len(gene_indices)), dtype="float32")
        for original_pos, sample_index in sorted_pairs:
            column = np.asarray(expression[:, sample_index], dtype="float32")
            counts[original_pos] = column[gene_indices]
    mask = query_adata.varm["I"][query_positions, :].copy()
    var = query_adata.var.iloc[query_positions, :].copy()
    return counts, var, mask, retained_genes


def run(args) -> Path:
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")

    tissue = args.tissue.strip().lower().replace(" ", "_")
    query_adata = ad.read_h5ad(args.query_h5ad)
    metadata = classify(load_sample_metadata(args.input))
    selected, selection = select_archs4_samples(
        metadata,
        tissue,
        args.max_samples,
        args.sample_seed,
    )
    counts, var, mask, retained_genes = load_reference_counts(
        args.input,
        selected,
        query_adata,
    )
    adata = ad.AnnData(X=counts, obs=selected.drop(columns=["filter_text"]), var=var)
    adata.layers["counts"] = counts.copy()
    adata.varm["I"] = mask
    adata.uns["terms"] = list(map(str, query_adata.uns["terms"]))
    adata.uns["term_descriptions"] = list(
        map(str, query_adata.uns.get("term_descriptions", query_adata.uns["terms"]))
    )
    adata.uns["expimap_preprocessing"] = {
        "source": "ARCHS4 mouse_gene_v2.5.h5",
        "tissue": tissue,
        "leakage_exclusion": "NASA/GeneLab/OSD/GLDS/spaceflight/microgravity/ISS/RR/hindlimb unloading terms",
        "reference_query_gene_order": str(args.query_h5ad),
        "transformation": "raw_counts",
        "recommended_recon_loss": "nb",
        "validity": "reference_primary_count_likelihood",
    }

    output_dir = Path(args.output_dir or f"outputs/expimap_archs4_reference_osdr_query_{tissue}/reference_input")
    output_dir.mkdir(parents=True, exist_ok=True)
    h5ad_path = output_dir / f"archs4_mouse_{tissue}_reference_reactome_raw_counts.h5ad"
    adata.write_h5ad(h5ad_path)
    selected_path = output_dir / "archs4_selected_samples.tsv"
    selected.drop(columns=["filter_text"]).to_csv(selected_path, sep="\t", index=False)
    manifest = {
        "tissue": tissue,
        "archs4_h5": str(args.input),
        "query_h5ad": str(args.query_h5ad),
        "max_samples": args.max_samples,
        "selection": selection,
        "counts": {
            "samples": int(adata.n_obs),
            "genes": int(adata.n_vars),
            "pathways": int(len(adata.uns["terms"])),
        },
        "outputs": {
            "h5ad": str(h5ad_path),
            "selected_samples": str(selected_path),
        },
    }
    manifest_path = output_dir / "reference_input_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare ARCHS4 tissue reference AnnData matching an OSDR expiMap input."
    )
    parser.add_argument("--input", default=DEFAULT_ARCHS4)
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--max-samples", type=int, default=500)
    parser.add_argument("--sample-seed", type=int, default=2020)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
