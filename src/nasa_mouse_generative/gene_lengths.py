"""Build versionless mouse Ensembl gene lengths from an official GENCODE GTF."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re

import pandas as pd


GENE_ID = re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')


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


def build_gene_lengths(gtf: Path) -> pd.DataFrame:
    intervals: dict[str, list[tuple[str, int, int]]] = {}
    with _open_text(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            match = GENE_ID.search(fields[8])
            if match is None:
                continue
            gene = match.group(1).split(".", 1)[0]
            if not gene.startswith("ENSMUSG"):
                continue
            intervals.setdefault(gene, []).append(
                (fields[0], int(fields[3]), int(fields[4]))
            )

    rows = []
    for gene, gene_intervals in intervals.items():
        merged_length = 0
        interval_count = 0
        chromosomes = set()
        by_chromosome: dict[str, list[tuple[int, int]]] = {}
        for chromosome, start, end in gene_intervals:
            chromosomes.add(chromosome)
            by_chromosome.setdefault(chromosome, []).append((start, end))
        for chromosome_intervals in by_chromosome.values():
            current_start = -1
            current_end = -1
            for start, end in sorted(chromosome_intervals):
                if start > current_end + 1:
                    if current_start >= 0:
                        merged_length += current_end - current_start + 1
                        interval_count += 1
                    current_start, current_end = start, end
                else:
                    current_end = max(current_end, end)
            if current_start >= 0:
                merged_length += current_end - current_start + 1
                interval_count += 1
        rows.append(
            {
                "gene_id": gene,
                "length_bp": merged_length,
                "merged_exon_intervals": interval_count,
                "chromosomes": ";".join(sorted(chromosomes)),
            }
        )
    return pd.DataFrame(rows).sort_values("gene_id", kind="stable").reset_index(
        drop=True
    )


def run(args: argparse.Namespace) -> Path:
    gtf = Path(args.gtf)
    if not gtf.exists():
        raise FileNotFoundError(gtf)
    table = build_gene_lengths(gtf)
    if table.empty:
        raise ValueError(f"No mouse Ensembl exon annotations found in {gtf}")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index=False)
    manifest = {
        "source_url": args.source_url,
        "source_file": str(gtf),
        "source_sha256": _sha256(gtf),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "genes": int(len(table)),
        "length_definition": "union of all annotated exon intervals per versionless gene_id",
        "output": str(output),
    }
    manifest_path = output.with_suffix(".manifest.json")
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
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
