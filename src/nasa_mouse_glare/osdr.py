"""Prepare OSDR mouse expression matrices for GLARE fine-tuning."""

from __future__ import annotations

import argparse
from pathlib import Path

from .io import require_import, write_matrix_bundle


DEFAULT_OSDR_MATRIX_KEY = "/data/expression"
DEFAULT_OSDR_GENE_KEY = "/meta/genes/ensembl_gene"
DEFAULT_OSDR_SAMPLE_KEY = "/meta/info/id.sample name"


def _decode_array(values) -> list[str]:
    out = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8", "replace"))
        else:
            out.append(str(value))
    return out


def inspect_h5(input_h5: str | Path) -> None:
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")
    with h5py.File(input_h5, "r") as handle:
        def visit(name, obj):
            if hasattr(obj, "shape"):
                print(f"/{name}\t{obj.shape}\t{obj.dtype}")
            else:
                print(f"/{name}/")

        handle.visititems(visit)


def prepare_osdr_h5(
    input_h5: str | Path,
    output_prefix: str | Path,
    matrix_key: str = DEFAULT_OSDR_MATRIX_KEY,
    gene_key: str = DEFAULT_OSDR_GENE_KEY,
    sample_key: str = DEFAULT_OSDR_SAMPLE_KEY,
    log1p: bool = False,
) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    with h5py.File(input_h5, "r") as handle:
        matrix = np.asarray(handle[matrix_key], dtype="float32")
        genes = _decode_array(handle[gene_key][:])
        profiles = _decode_array(handle[sample_key][:])

        metadata = pd.DataFrame({"profile": profiles})
        for key in [
            "/meta/info/id.accession",
            "/meta/info/investigation.study.comment.project type",
            "/meta/info/investigation.study assays.study assay technology type",
        ]:
            if key in handle:
                metadata[Path(key).name] = _decode_array(handle[key][:])

    if log1p:
        matrix = np.log1p(matrix)

    return write_matrix_bundle(
        output_prefix,
        matrix,
        genes=genes,
        profiles=profiles,
        profile_metadata=metadata,
        description=f"OSDR expression from {input_h5}",
    )


def prepare_osdr_csv(
    input_csv: str | Path,
    output_prefix: str | Path,
    gene_col: str = "gene_id",
    log1p: bool = False,
) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    df = pd.read_csv(input_csv)
    if gene_col not in df.columns:
        raise ValueError(f"gene_col '{gene_col}' not present in CSV")
    genes = df[gene_col].astype(str).tolist()
    expression = df.drop(columns=[gene_col])
    profiles = expression.columns.astype(str).tolist()
    matrix = expression.to_numpy(dtype="float32")
    if log1p:
        matrix = np.log1p(matrix)

    return write_matrix_bundle(
        output_prefix,
        matrix,
        genes=genes,
        profiles=profiles,
        description=f"OSDR CSV expression from {input_csv}",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare OSDR expression for GLARE.")
    sub = parser.add_subparsers(dest="command", required=True)

    inspect = sub.add_parser("inspect-h5")
    inspect.add_argument("--input", required=True)

    h5 = sub.add_parser("prep-h5")
    h5.add_argument("--input", required=True)
    h5.add_argument("--output-prefix", required=True)
    h5.add_argument("--matrix-key", default=DEFAULT_OSDR_MATRIX_KEY)
    h5.add_argument("--gene-key", default=DEFAULT_OSDR_GENE_KEY)
    h5.add_argument("--sample-key", default=DEFAULT_OSDR_SAMPLE_KEY)
    h5.add_argument("--log1p", action="store_true")

    csv_parser = sub.add_parser("prep-csv")
    csv_parser.add_argument("--input", required=True)
    csv_parser.add_argument("--output-prefix", required=True)
    csv_parser.add_argument("--gene-col", default="gene_id")
    csv_parser.add_argument("--log1p", action="store_true")

    args = parser.parse_args()
    if args.command == "inspect-h5":
        inspect_h5(args.input)
    elif args.command == "prep-h5":
        print(
            prepare_osdr_h5(
                args.input,
                args.output_prefix,
                args.matrix_key,
                args.gene_key,
                args.sample_key,
                args.log1p,
            )
        )
    elif args.command == "prep-csv":
        print(prepare_osdr_csv(args.input, args.output_prefix, args.gene_col, args.log1p))


if __name__ == "__main__":
    main()
