"""Create a documented expiMap query subset by excluding accessions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nasa_mouse_glare.io import require_import


def run(args: argparse.Namespace) -> Path:
    ad = require_import("anndata", "pip install -r requirements.txt")

    input_path = Path(args.input)
    output_path = Path(args.output)
    if input_path.resolve() == output_path.resolve():
        raise SystemExit("Input and output paths must differ.")

    query = ad.read_h5ad(input_path)
    if args.accession_column not in query.obs:
        raise SystemExit(
            f"Input is missing accession column {args.accession_column!r}."
        )

    excluded = tuple(dict.fromkeys(map(str, args.exclude_accession)))
    if not excluded:
        raise SystemExit("At least one --exclude-accession value is required.")

    accessions = query.obs[args.accession_column].astype(str)
    present = set(accessions)
    missing = sorted(set(excluded) - present)
    if missing:
        raise SystemExit(
            "Requested exclusions are absent from the input: " + ", ".join(missing)
        )

    exclusion_counts = {
        accession: int(accessions.eq(accession).sum()) for accession in excluded
    }
    keep = ~accessions.isin(excluded)
    subset = query[keep.to_numpy()].copy()
    if subset.n_obs == 0:
        raise SystemExit("The requested exclusions removed every query sample.")
    if args.expected_samples is not None and subset.n_obs != args.expected_samples:
        raise SystemExit(
            f"Expected {args.expected_samples} retained samples, found {subset.n_obs}."
        )

    provenance = {
        "source": str(input_path),
        "accession_column": args.accession_column,
        "excluded_accessions": list(excluded),
        "excluded_sample_counts": exclusion_counts,
        "samples_before": int(query.n_obs),
        "samples_after": int(subset.n_obs),
        "genes": int(subset.n_vars),
    }
    subset.uns["expimap_query_subset"] = provenance

    output_path.parent.mkdir(parents=True, exist_ok=True)
    subset.write_h5ad(output_path)

    manifest_path = (
        Path(args.manifest)
        if args.manifest
        else output_path.with_name(f"{output_path.stem}_manifest.json")
    )
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {**provenance, "output": str(output_path), "manifest": str(manifest_path)}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Prepared expiMap query H5AD.")
    parser.add_argument("--output", required=True, help="Subset query H5AD to write.")
    parser.add_argument(
        "--exclude-accession",
        action="append",
        default=[],
        help="Accession to remove. Repeat for multiple accessions.",
    )
    parser.add_argument("--accession-column", default="id.accession")
    parser.add_argument("--expected-samples", type=int)
    parser.add_argument("--manifest")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
