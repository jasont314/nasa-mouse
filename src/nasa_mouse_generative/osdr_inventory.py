"""Build current, canonical OSDR tissue inventories and training tiers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import pandas as pd

from nasa_mouse_glare.fetch_osdr_mouse_transcriptomics import discover_metadata

from .config import DataConfig
from .tissues import canonicalize_material


REQUIRED_COLUMNS = {
    "id.accession",
    "id.sample name",
    "condition_inferred",
    "study.characteristics.material type",
}


def load_metadata(path: str | Path, *, refresh: bool, timeout: int) -> pd.DataFrame:
    if refresh:
        return discover_metadata(timeout)
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise ValueError(f"OSDR metadata is missing columns: {sorted(missing)}")
    return frame


def canonicalize(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["tissue_api_original"] = result.get("tissue_final", "unknown").astype(str)
    result["material_type_original"] = result["study.characteristics.material type"].fillna("").astype(str)
    result["tissue_canonical"] = result["material_type_original"].map(canonicalize_material)
    unresolved = result["tissue_canonical"].isin({"unknown", "unspecified"})
    result.loc[unresolved, "tissue_canonical"] = result.loc[unresolved, "tissue_api_original"]
    result["tissue_alias_changed"] = result["tissue_canonical"] != result["tissue_api_original"]
    result["biological_sample_name"] = result["id.sample name"].astype(str).map(
        lambda value: re.sub(r"_techrep\d+$", "", value, flags=re.IGNORECASE)
    )
    result["biological_profile_id"] = (
        result["id.accession"].astype(str)
        + "/"
        + result["biological_sample_name"].astype(str)
    )
    return result


def _tier(row: pd.Series, thresholds: DataConfig) -> str:
    confirmatory = (
        row["training_samples"] >= thresholds.min_confirmatory_total
        and row["training_min_condition_samples"]
        >= thresholds.min_confirmatory_per_condition
        and row["accessions_with_both_conditions"] >= thresholds.min_confirmatory_accessions
    )
    if confirmatory:
        return "confirmatory_per_tissue"
    exploratory = (
        row["training_samples"] >= thresholds.min_exploratory_total
        and row["training_min_condition_samples"]
        >= thresholds.min_exploratory_per_condition
        and row["accessions_with_both_conditions"] >= thresholds.min_exploratory_accessions
    )
    if exploratory:
        return "exploratory_pretrained_per_tissue"
    return "pooled_only"


def build_inventory(frame: pd.DataFrame, thresholds: DataConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    canonical = canonicalize(frame)
    counts = (
        canonical.groupby(["tissue_canonical", "condition_inferred"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    for condition in ("flight", "ground_control"):
        if condition not in counts:
            counts[condition] = 0
    training_counts = (
        canonical.groupby(["tissue_canonical", "condition_inferred"], dropna=False)[
            "biological_profile_id"
        ]
        .nunique()
        .unstack(fill_value=0)
    )
    for condition in ("flight", "ground_control"):
        if condition not in training_counts:
            training_counts[condition] = 0
    training_counts = training_counts.rename(
        columns={
            "flight": "training_flight",
            "ground_control": "training_ground_control",
        }
    )
    accessions = (
        canonical.groupby(["tissue_canonical", "id.accession", "condition_inferred"])
        .size()
        .unstack(fill_value=0)
    )
    for condition in ("flight", "ground_control"):
        if condition not in accessions:
            accessions[condition] = 0
    accessions["has_both"] = (accessions["flight"] > 0) & (accessions["ground_control"] > 0)
    accession_summary = accessions.reset_index().groupby("tissue_canonical").agg(
        accessions=("id.accession", "nunique"),
        accessions_with_both_conditions=("has_both", "sum"),
    )
    inventory = counts.join(training_counts).join(accession_summary).reset_index()
    inventory["total_samples"] = inventory["flight"] + inventory["ground_control"]
    inventory["min_condition_samples"] = inventory[["flight", "ground_control"]].min(axis=1)
    inventory["training_samples"] = (
        inventory["training_flight"] + inventory["training_ground_control"]
    )
    inventory["training_min_condition_samples"] = inventory[
        ["training_flight", "training_ground_control"]
    ].min(axis=1)
    accession_lists = canonical.groupby("tissue_canonical")["id.accession"].agg(
        lambda values: ";".join(sorted(set(map(str, values))))
    )
    inventory["accession_list"] = inventory["tissue_canonical"].map(accession_lists)
    inventory["training_tier"] = inventory.apply(_tier, axis=1, thresholds=thresholds)
    inventory["pooled_conditioned_included"] = True
    inventory = inventory.sort_values(
        ["total_samples", "tissue_canonical"], ascending=[False, True]
    ).reset_index(drop=True)
    return canonical, inventory


def run(args: argparse.Namespace) -> Path:
    thresholds = DataConfig(
        osdr_metadata=args.metadata,
        min_confirmatory_total=args.min_confirmatory_total,
        min_confirmatory_per_condition=args.min_confirmatory_per_condition,
        min_confirmatory_accessions=args.min_confirmatory_accessions,
        min_exploratory_total=args.min_exploratory_total,
        min_exploratory_per_condition=args.min_exploratory_per_condition,
        min_exploratory_accessions=args.min_exploratory_accessions,
    )
    frame = load_metadata(args.metadata, refresh=args.refresh, timeout=args.timeout)
    canonical, inventory = build_inventory(frame, thresholds)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = output_dir / "osdr_canonical_metadata.tsv"
    inventory_path = output_dir / "osdr_tissue_inventory.tsv"
    alias_path = output_dir / "osdr_tissue_alias_audit.tsv"
    canonical.to_csv(metadata_path, sep="\t", index=False)
    inventory.to_csv(inventory_path, sep="\t", index=False)
    (
        canonical.groupby(
            ["material_type_original", "tissue_api_original", "tissue_canonical"],
            dropna=False,
        )
        .size()
        .rename("samples")
        .reset_index()
        .sort_values(["tissue_canonical", "samples"], ascending=[True, False])
        .to_csv(alias_path, sep="\t", index=False)
    )
    summary = {
        "source": "NASA OSDR Biological Data API" if args.refresh else str(args.metadata),
        "samples": int(len(canonical)),
        "biological_samples_after_technical_replicate_collapse": int(
            canonical["biological_profile_id"].nunique()
        ),
        "accessions": int(canonical["id.accession"].nunique()),
        "canonical_tissues": int(inventory.shape[0]),
        "training_tiers": {
            str(key): int(value)
            for key, value in inventory["training_tier"].value_counts().sort_index().items()
        },
        "outputs": {
            "metadata": str(metadata_path),
            "inventory": str(inventory_path),
            "alias_audit": str(alias_path),
        },
    }
    summary_path = output_dir / "osdr_inventory_summary.json"
    summary["outputs"]["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return inventory_path


def parse_args() -> argparse.Namespace:
    defaults = DataConfig()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default=defaults.osdr_metadata)
    parser.add_argument("--output-dir", default="outputs/generative_benchmark/data_audit/osdr")
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--min-confirmatory-total", type=int, default=defaults.min_confirmatory_total)
    parser.add_argument(
        "--min-confirmatory-per-condition",
        type=int,
        default=defaults.min_confirmatory_per_condition,
    )
    parser.add_argument(
        "--min-confirmatory-accessions", type=int, default=defaults.min_confirmatory_accessions
    )
    parser.add_argument("--min-exploratory-total", type=int, default=defaults.min_exploratory_total)
    parser.add_argument(
        "--min-exploratory-per-condition",
        type=int,
        default=defaults.min_exploratory_per_condition,
    )
    parser.add_argument(
        "--min-exploratory-accessions", type=int, default=defaults.min_exploratory_accessions
    )
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
