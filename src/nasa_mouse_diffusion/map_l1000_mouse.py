"""Map the paper's human L1000 landmark genes to mouse Ensembl IDs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import textwrap

from nasa_mouse_glare.io import require_import


def load_human_landmarks(path: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    table = pd.read_csv(path)
    if {"Description", "ensembl_id", "Type"}.issubset(table.columns):
        table = table.loc[table["Type"].astype(str).str.lower().eq("landmark")]
        return table[["Description", "ensembl_id"]].rename(
            columns={"Description": "human_symbol", "ensembl_id": "human_ensembl_gene"}
        )
    if {"Symbol", "Type"}.issubset(table.columns):
        table = table.loc[table["Type"].astype(str).str.lower().eq("landmark")]
        table = table.rename(columns={"Symbol": "human_symbol"})
        table["human_ensembl_gene"] = ""
        return table[["human_symbol", "human_ensembl_gene"]]
    raise SystemExit(f"Unsupported landmark table columns: {table.columns.tolist()}")


def biomart_query(human_ensembl_ids: list[str], *, timeout: int):
    requests = require_import("requests", "pip install requests")
    query = f"""
    <!DOCTYPE Query>
    <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
      <Dataset name="hsapiens_gene_ensembl" interface="default">
        <Filter name="ensembl_gene_id" value="{','.join(human_ensembl_ids)}"/>
        <Attribute name="ensembl_gene_id"/>
        <Attribute name="external_gene_name"/>
        <Attribute name="mmusculus_homolog_ensembl_gene"/>
        <Attribute name="mmusculus_homolog_associated_gene_name"/>
        <Attribute name="mmusculus_homolog_orthology_type"/>
        <Attribute name="mmusculus_homolog_perc_id"/>
      </Dataset>
    </Query>
    """
    response = requests.get(
        "https://www.ensembl.org/biomart/martservice",
        params={"query": textwrap.dedent(query).strip()},
        timeout=timeout,
    )
    response.raise_for_status()
    return response.text


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    landmarks = load_human_landmarks(Path(args.landmark_table))
    ids = [value for value in landmarks["human_ensembl_gene"].dropna().astype(str).tolist() if value.startswith("ENSG")]
    if not ids:
        raise SystemExit("Landmark table does not contain human Ensembl IDs; use GTEx description table from the paper repo.")
    rows = []
    for start in range(0, len(ids), args.chunk_size):
        text = biomart_query(ids[start : start + args.chunk_size], timeout=args.timeout)
        chunk = pd.read_csv(__import__("io").StringIO(text), sep="\t")
        rows.append(chunk)
    mapped = pd.concat(rows, ignore_index=True)
    mapped.columns = [
        "human_ensembl_gene",
        "human_symbol_biomart",
        "mouse_ensembl_gene",
        "mouse_symbol",
        "orthology_type",
        "mouse_perc_id",
    ]
    merged = landmarks.merge(mapped, on="human_ensembl_gene", how="left")
    merged = merged.loc[merged["mouse_ensembl_gene"].fillna("").astype(str).str.startswith("ENSMUSG")]
    merged = merged.drop_duplicates(["human_ensembl_gene", "mouse_ensembl_gene"])
    merged.to_csv(output, sep="\t", index=False)
    manifest = {
        "landmark_table": str(args.landmark_table),
        "output": str(output),
        "n_human_landmarks": int(len(landmarks)),
        "n_mapped_rows": int(len(merged)),
        "n_unique_mouse_genes": int(merged["mouse_ensembl_gene"].nunique()) if not merged.empty else 0,
    }
    manifest_path = output.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--landmark-table", default="/tmp/rna-diffusion-paper/data/gtex_description.csv")
    parser.add_argument("--output", default="data/diffusion/l1000_human_to_mouse_ensembl.tsv")
    parser.add_argument("--chunk-size", type=int, default=200)
    parser.add_argument("--timeout", type=int, default=60)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
