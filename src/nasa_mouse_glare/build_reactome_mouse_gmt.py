"""Build a current mouse Reactome GMT using Ensembl gene IDs."""

from __future__ import annotations

import argparse
import json
import re
import urllib.request
from collections import defaultdict
from pathlib import Path


REACTOME_PATHWAYS_URL = "https://reactome.org/download/current/ReactomePathways.txt"
ENSEMBL2REACTOME_URL = (
    "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt"
)
PATHWAY_BROWSER_URL = "https://reactome.org/PathwayBrowser/#/"


def sanitize_name(name: str) -> str:
    label = re.sub(r"[^A-Za-z0-9]+", "_", name.upper()).strip("_")
    return re.sub(r"_+", "_", label)


def download_if_missing(url: str, path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response, path.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def read_mouse_pathways(path: Path) -> dict[str, str]:
    pathway_names: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            reactome_id, name, species = parts[:3]
            if species == "Mus musculus" and reactome_id.startswith("R-MMU-"):
                pathway_names[reactome_id] = name
    return pathway_names


def read_mouse_ensembl_mapping(path: Path, pathway_names: dict[str, str]) -> dict[str, set[str]]:
    pathway_to_genes: dict[str, set[str]] = defaultdict(set)
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            gene_id, pathway_id, _url, pathway_name, _evidence, species = parts[:6]
            if species != "Mus musculus":
                continue
            if not pathway_id.startswith("R-MMU-"):
                continue
            if not gene_id.startswith("ENSMUSG"):
                continue
            pathway_names.setdefault(pathway_id, pathway_name)
            pathway_to_genes[pathway_id].add(gene_id)
    return pathway_to_genes


def build_gmt(raw_dir: Path, output_dir: Path, download: bool) -> None:
    pathways_file = raw_dir / "ReactomePathways.txt"
    mapping_file = raw_dir / "Ensembl2Reactome_All_Levels.txt"

    if download:
        download_if_missing(REACTOME_PATHWAYS_URL, pathways_file)
        download_if_missing(ENSEMBL2REACTOME_URL, mapping_file)

    if not pathways_file.exists():
        raise FileNotFoundError(pathways_file)
    if not mapping_file.exists():
        raise FileNotFoundError(mapping_file)

    output_dir.mkdir(parents=True, exist_ok=True)
    pathway_names = read_mouse_pathways(pathways_file)
    pathway_to_genes = read_mouse_ensembl_mapping(mapping_file, pathway_names)

    rows = []
    for reactome_id, genes in pathway_to_genes.items():
        if not genes:
            continue
        name = pathway_names.get(reactome_id, reactome_id)
        term = f"{reactome_id}_{sanitize_name(name)}"
        url = f"{PATHWAY_BROWSER_URL}{reactome_id}"
        rows.append((reactome_id, term, name, url, sorted(genes)))
    rows.sort(key=lambda row: row[1])

    gmt_path = output_dir / "reactome_current_mouse_ensembl.gmt"
    terms_path = output_dir / "reactome_current_mouse_ensembl_terms.tsv"
    genes_path = output_dir / "reactome_current_mouse_ensembl_genes.tsv"
    manifest_path = output_dir / "reactome_current_mouse_ensembl_manifest.json"

    with gmt_path.open("w") as out:
        for _reactome_id, term, _name, url, genes in rows:
            out.write("\t".join([term, url, *genes]) + "\n")

    with terms_path.open("w") as out:
        out.write("term\treactome_id\tname\turl\tn_genes\n")
        for reactome_id, term, name, url, genes in rows:
            out.write(f"{term}\t{reactome_id}\t{name}\t{url}\t{len(genes)}\n")

    all_genes = sorted({gene for *_prefix, genes in rows for gene in genes})
    with genes_path.open("w") as out:
        out.write("ensembl_gene\n")
        for gene in all_genes:
            out.write(gene + "\n")

    manifest = {
        "source": "Reactome current download",
        "source_urls": {
            pathways_file.name: REACTOME_PATHWAYS_URL,
            mapping_file.name: ENSEMBL2REACTOME_URL,
        },
        "filters": {
            "species": "Mus musculus",
            "pathway_id_prefix": "R-MMU-",
            "gene_id_prefix": "ENSMUSG",
        },
        "format": {
            "gmt_columns": [
                "R-MMU pathway ID plus sanitized pathway name",
                "Reactome Pathway Browser URL",
                "mouse Ensembl gene IDs",
            ],
            "example_term": "R-MMU-73857_RNA_POLYMERASE_II_TRANSCRIPTION",
        },
        "outputs": {
            "gmt": str(gmt_path),
            "terms": str(terms_path),
            "genes": str(genes_path),
        },
        "counts": {
            "mouse_pathways_in_reactome_pathways": len(
                [pid for pid in pathway_names if pid.startswith("R-MMU-")]
            ),
            "pathways_with_ensembl_genes": len(rows),
            "unique_ensembl_genes": len(all_genes),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    print(json.dumps(manifest["counts"], indent=2))
    print(f"wrote {gmt_path}")
    print(f"wrote {terms_path}")
    print(f"wrote {genes_path}")
    print(f"wrote {manifest_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a current Mus musculus Reactome GMT with ENSMUSG gene IDs."
    )
    parser.add_argument(
        "--raw-dir",
        default="data/pathways/reactome_current/raw",
        help="Directory containing ReactomePathways.txt and Ensembl2Reactome_All_Levels.txt.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/pathways",
        help="Directory for generated GMT, term table, gene table, and manifest.",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Require raw Reactome files to already exist instead of downloading missing files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_gmt(
        raw_dir=Path(args.raw_dir),
        output_dir=Path(args.output_dir),
        download=not args.no_download,
    )


if __name__ == "__main__":
    main()
