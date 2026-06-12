"""Prepare Tabula Muris Senis data for GLARE-style pretraining."""

from __future__ import annotations

import argparse
from pathlib import Path

from .io import require_import, write_matrix_bundle


def parse_filter(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("filters must be KEY=VALUE")
    key, filter_value = value.split("=", 1)
    return key, filter_value


def _matrix_from_adata(adata, matrix_source: str):
    if matrix_source == "X":
        return adata.X, adata.var.copy()
    if matrix_source == "raw.X":
        if adata.raw is None:
            raise ValueError("Requested raw.X, but this h5ad has no raw matrix")
        return adata.raw.X, adata.raw.var.copy()
    if matrix_source.startswith("layers/"):
        layer = matrix_source.split("/", 1)[1]
        return adata.layers[layer], adata.var.copy()
    raise ValueError("matrix_source must be X, raw.X, or layers/<name>")


def _gene_ids(var, gene_field: str) -> list[str]:
    if gene_field == "index":
        return [str(v) for v in var.index]
    if gene_field not in var.columns:
        raise ValueError(
            f"Gene field '{gene_field}' not in var columns. Available: "
            f"{', '.join(map(str, var.columns[:20]))}"
        )
    return [str(v) for v in var[gene_field].tolist()]


def prepare_tms(
    input_h5ad: str | Path,
    output_prefix: str | Path,
    matrix_source: str = "X",
    gene_field: str = "index",
    max_cells: int | None = None,
    obs_filters: list[tuple[str, str]] | None = None,
    random_seed: int = 2026,
) -> Path:
    """Write a direct-cell TMS matrix bundle for GLARE-style pretraining."""
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_sparse = require_import(
        "scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt"
    )
    anndata = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")

    adata_backed = anndata.read_h5ad(input_h5ad, backed="r")
    obs = adata_backed.obs.copy()

    mask = pd.Series(True, index=obs.index)
    for key, value in obs_filters or []:
        if key not in obs.columns:
            raise ValueError(f"Filter field '{key}' not found in obs")
        mask &= obs[key].astype(str).eq(value)

    selected_obs = obs.index[mask].to_numpy()
    if max_cells is not None and len(selected_obs) > max_cells:
        rng = np.random.default_rng(random_seed)
        selected_obs = rng.choice(selected_obs, size=max_cells, replace=False)

    adata = adata_backed[selected_obs, :].to_memory()
    try:
        adata_backed.file.close()
    except AttributeError:
        pass
    obs = adata.obs.copy()
    X, var = _matrix_from_adata(adata, matrix_source)
    genes = _gene_ids(var, gene_field)

    if max_cells is None and X.shape[0] > 10_000:
        raise SystemExit(
            "Direct-cell pretraining matrix is large. Provide --max-cells for "
            "a GLARE-size sampled run, or intentionally raise this guard in code "
            "after confirming memory/GPU capacity."
        )

    matrix = X.T.tocsr() if scipy_sparse.issparse(X) else np.asarray(X).T
    profile_names = [str(v) for v in obs.index.tolist()]
    index_name = "cell_id" if "cell" in obs.columns else "cell"
    metadata = obs.reset_index(names=index_name)
    description = f"TMS {matrix_source} direct cell profiles"

    return write_matrix_bundle(
        output_prefix,
        matrix,
        genes=genes,
        profiles=profile_names,
        profile_metadata=metadata,
        description=description,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare TMS h5ad for GLARE pretraining.")
    parser.add_argument("--input", required=True, help="TMS .h5ad file")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument(
        "--matrix-source",
        default="X",
        help="X, raw.X, or layers/<name>. GLARE used normalized counts; inspect the h5ad before choosing.",
    )
    parser.add_argument("--gene-field", default="index")
    parser.add_argument("--max-cells", type=int)
    parser.add_argument("--filter", action="append", type=parse_filter, default=[])
    parser.add_argument("--random-seed", type=int, default=2026)
    args = parser.parse_args()

    manifest = prepare_tms(
        input_h5ad=args.input,
        output_prefix=args.output_prefix,
        matrix_source=args.matrix_source,
        gene_field=args.gene_field,
        max_cells=args.max_cells,
        obs_filters=args.filter,
        random_seed=args.random_seed,
    )
    print(manifest)


if __name__ == "__main__":
    main()
