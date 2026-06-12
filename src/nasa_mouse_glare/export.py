"""Export prepared bundles to upstream GLARE-compatible files."""

from __future__ import annotations

import argparse
from pathlib import Path

from .io import dense_matrix, load_matrix_bundle, require_import


def export_mtx(bundle_manifest: str | Path, output_mtx: str | Path) -> None:
    """Export a bundle matrix as MatrixMarket for upstream GLARE hpt.py."""
    scipy_io = require_import("scipy.io", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_sparse = require_import(
        "scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt"
    )
    bundle = load_matrix_bundle(bundle_manifest)
    matrix = bundle.matrix
    if not scipy_sparse.issparse(matrix):
        matrix = scipy_sparse.coo_matrix(matrix)

    output_mtx = Path(output_mtx)
    output_mtx.parent.mkdir(parents=True, exist_ok=True)
    scipy_io.mmwrite(output_mtx, matrix)
    output_mtx.with_suffix(output_mtx.suffix + ".genes.tsv").write_text(
        "\n".join(bundle.genes) + "\n", encoding="utf-8"
    )
    output_mtx.with_suffix(output_mtx.suffix + ".profiles.tsv").write_text(
        "\n".join(bundle.profiles) + "\n", encoding="utf-8"
    )


def export_csv(
    bundle_manifest: str | Path,
    output_csv: str | Path,
    include_gene_id: bool = False,
    max_dense_gb: float = 8.0,
) -> None:
    """Export a bundle matrix as CSV for upstream GLARE fine-tuning scripts."""
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    bundle = load_matrix_bundle(bundle_manifest)
    matrix = dense_matrix(bundle.matrix, max_dense_gb=max_dense_gb)
    df = pd.DataFrame(matrix, columns=bundle.profiles)
    if include_gene_id:
        df.insert(0, "gene_id", bundle.genes)
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Export nasa_mouse_glare bundles.")
    sub = parser.add_subparsers(dest="command", required=True)

    mtx = sub.add_parser("mtx")
    mtx.add_argument("--bundle", required=True)
    mtx.add_argument("--output", required=True)

    csv_parser = sub.add_parser("csv")
    csv_parser.add_argument("--bundle", required=True)
    csv_parser.add_argument("--output", required=True)
    csv_parser.add_argument("--include-gene-id", action="store_true")
    csv_parser.add_argument("--max-dense-gb", type=float, default=8.0)

    args = parser.parse_args()
    if args.command == "mtx":
        export_mtx(args.bundle, args.output)
    elif args.command == "csv":
        export_csv(args.bundle, args.output, args.include_gene_id, args.max_dense_gb)


if __name__ == "__main__":
    main()
