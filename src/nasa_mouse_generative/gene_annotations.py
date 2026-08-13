"""Build a mouse Ensembl gene-symbol map from an official GENCODE GTF."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re

import pandas as pd


GENE_ID = re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')
GENE_NAME = re.compile(r'(?:^|;\s*)gene_name\s+"([^"]+)"')
GENE_TYPE = re.compile(r'(?:^|;\s*)gene_type\s+"([^"]+)"')


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_gene_annotations(gtf: Path) -> pd.DataFrame:
    """Return one versionless Ensembl ID, symbol, and type per mouse gene."""

    records: dict[str, tuple[str, str]] = {}
    with _open_text(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene_match = GENE_ID.search(fields[8])
            if gene_match is None:
                continue
            gene_id = gene_match.group(1).split(".", 1)[0]
            if not gene_id.startswith("ENSMUSG"):
                continue
            name_match = GENE_NAME.search(fields[8])
            type_match = GENE_TYPE.search(fields[8])
            symbol = name_match.group(1) if name_match else gene_id
            gene_type = type_match.group(1) if type_match else ""
            previous = records.get(gene_id)
            current = (symbol, gene_type)
            if previous is not None and previous != current:
                raise ValueError(
                    f"Conflicting GENCODE annotations for {gene_id}: "
                    f"{previous!r} versus {current!r}"
                )
            records[gene_id] = current

    rows = [
        {"gene_id": gene_id, "gene_symbol": values[0], "gene_type": values[1]}
        for gene_id, values in records.items()
    ]
    if not rows:
        raise ValueError(f"No mouse Ensembl gene annotations found in {gtf}")
    return pd.DataFrame(rows).sort_values("gene_id", kind="stable").reset_index(
        drop=True
    )


def _manifest_path(output: Path) -> Path:
    stem = output.with_suffix("") if output.suffix == ".gz" else output
    return stem.with_suffix(".manifest.json")


def run(args: argparse.Namespace) -> Path:
    gtf = Path(args.gtf)
    if not gtf.is_file():
        raise FileNotFoundError(gtf)
    table = build_gene_annotations(gtf)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    compression = {"method": "gzip", "mtime": 0} if output.suffix == ".gz" else None
    table.to_csv(
        output,
        sep="\t",
        index=False,
        compression=compression,
        lineterminator="\n",
    )
    manifest = {
        "source_url": args.source_url,
        "source_file": str(gtf),
        "source_sha256": _sha256(gtf),
        "genes": int(len(table)),
        "output": str(output),
    }
    _manifest_path(output).write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-url", default="")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
