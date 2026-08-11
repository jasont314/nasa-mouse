"""Create accession-grouped validation plans without inspecting expression values."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import pandas as pd


def _stable_key(seed: int, *parts: object) -> int:
    payload = ":".join([str(seed), *(str(part) for part in parts)]).encode(
        "utf-8", "replace"
    )
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big", signed=False)


def _accession_condition_table(metadata: pd.DataFrame) -> pd.DataFrame:
    required = {"id.accession", "tissue_canonical", "condition_inferred"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Canonical OSDR metadata is missing {sorted(missing)}")
    sample_column = (
        "biological_profile_id"
        if "biological_profile_id" in metadata
        else "profile_id"
    )
    table = metadata.groupby(
        ["tissue_canonical", "id.accession", "condition_inferred"],
        dropna=False,
    )[sample_column].nunique().unstack(fill_value=0)
    for condition in ("flight", "ground_control"):
        if condition not in table:
            table[condition] = 0
    table = table.reset_index()
    table["has_both_conditions"] = (table["flight"] > 0) & (
        table["ground_control"] > 0
    )
    table["samples"] = table["flight"] + table["ground_control"]
    return table


def build_per_tissue_plans(
    metadata: pd.DataFrame, inventory: pd.DataFrame, *, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    accession_counts = _accession_condition_table(metadata)
    tier_by_tissue = inventory.set_index("tissue_canonical")["training_tier"].to_dict()
    locked_rows: list[dict] = []
    loo_rows: list[dict] = []

    for tissue, frame in accession_counts.groupby("tissue_canonical", sort=True):
        frame = frame.rename(columns={"id.accession": "accession"})
        tier = str(tier_by_tissue.get(tissue, "pooled_only"))
        paired = frame.loc[frame["has_both_conditions"]].copy()
        paired["stable_key"] = paired["accession"].map(
            lambda accession: _stable_key(seed, tissue, accession)
        )
        paired = paired.sort_values("stable_key", kind="stable")
        accessions = paired["accession"].astype(str).tolist()

        if tier == "confirmatory_per_tissue" and len(accessions) >= 5:
            role_by_accession = {
                accession: (
                    "locked_test"
                    if index == 0
                    else "validation"
                    if index == 1
                    else "training"
                )
                for index, accession in enumerate(accessions)
            }
            for row in frame.itertuples(index=False):
                accession = str(row.accession)
                locked_rows.append(
                    {
                        "tissue_canonical": tissue,
                        "training_tier": tier,
                        "id.accession": accession,
                        "role": role_by_accession.get(accession, "training_unpaired"),
                        "flight": int(row.flight),
                        "ground_control": int(row.ground_control),
                        "samples": int(row.samples),
                    }
                )

        if len(accessions) >= 2 and tier != "pooled_only":
            for held_out in accessions:
                fold_id = f"{tissue}__holdout_{held_out}"
                for row in frame.itertuples(index=False):
                    accession = str(row.accession)
                    loo_rows.append(
                        {
                            "fold_id": fold_id,
                            "tissue_canonical": tissue,
                            "training_tier": tier,
                            "id.accession": accession,
                            "role": "held_out" if accession == held_out else "training",
                            "flight": int(row.flight),
                            "ground_control": int(row.ground_control),
                            "samples": int(row.samples),
                        }
                    )

    return pd.DataFrame(locked_rows), pd.DataFrame(loo_rows)


def build_pooled_plan(
    metadata: pd.DataFrame,
    *,
    seed: int,
    validation_fraction: float,
    test_fraction: float,
) -> pd.DataFrame:
    membership = (
        metadata.groupby(
            ["id.accession", "tissue_canonical", "condition_inferred"],
            dropna=False,
        )
        .size()
        .rename("samples")
        .reset_index()
    )
    accessions = sorted(membership["id.accession"].astype(str).unique())
    roles = {accession: "training" for accession in accessions}
    pairs_by_accession = {
        str(accession): set(
            zip(frame["tissue_canonical"].astype(str), frame["condition_inferred"].astype(str))
        )
        for accession, frame in membership.groupby("id.accession")
    }

    def can_remove_from_training(accession: str) -> bool:
        for pair in pairs_by_accession[accession]:
            remaining = sum(
                roles[other] == "training" and pair in pairs_by_accession[other]
                for other in accessions
                if other != accession
            )
            if remaining < 1:
                return False
        return True

    def assign(role: str, fraction: float, salt: str) -> None:
        target = int(math.ceil(len(accessions) * fraction))
        candidates = sorted(
            accessions,
            key=lambda accession: _stable_key(seed, salt, accession),
        )
        assigned = 0
        for accession in candidates:
            if assigned >= target:
                break
            if roles[accession] != "training" or not can_remove_from_training(accession):
                continue
            roles[accession] = role
            assigned += 1

    assign("locked_test", test_fraction, "test")
    assign("validation", validation_fraction, "validation")

    sample_column = (
        "biological_profile_id"
        if "biological_profile_id" in metadata
        else "profile_id"
    )
    accession_summary = metadata.groupby("id.accession").agg(
        samples=(sample_column, "nunique"),
        tissues=("tissue_canonical", lambda values: ";".join(sorted(set(map(str, values))))),
        conditions=(
            "condition_inferred",
            lambda values: ";".join(sorted(set(map(str, values)))),
        ),
    )
    result = accession_summary.reset_index()
    result["role"] = result["id.accession"].astype(str).map(roles)
    result["seed"] = int(seed)
    return result.sort_values(["role", "id.accession"], kind="stable").reset_index(
        drop=True
    )


def run(args: argparse.Namespace) -> Path:
    metadata = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    inventory = pd.read_csv(args.inventory, sep="\t")
    locked, loo = build_per_tissue_plans(metadata, inventory, seed=args.seed)
    pooled = build_pooled_plan(
        metadata,
        seed=args.seed,
        validation_fraction=args.validation_fraction,
        test_fraction=args.test_fraction,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    locked_path = output_dir / "per_tissue_locked_accession_splits.tsv"
    loo_path = output_dir / "per_tissue_loo_accession_folds.tsv"
    pooled_path = output_dir / "pooled_accession_split.tsv"
    locked.to_csv(locked_path, sep="\t", index=False)
    loo.to_csv(loo_path, sep="\t", index=False)
    pooled.to_csv(pooled_path, sep="\t", index=False)

    summary = {
        "split_unit": "OSDR accession",
        "seed": int(args.seed),
        "expression_values_used_for_split": False,
        "pooled_role_counts": {
            str(key): int(value) for key, value in pooled["role"].value_counts().items()
        },
        "per_tissue_locked_tissues": int(locked["tissue_canonical"].nunique())
        if not locked.empty
        else 0,
        "per_tissue_loo_tissues": int(loo["tissue_canonical"].nunique())
        if not loo.empty
        else 0,
        "per_tissue_loo_folds": int(loo["fold_id"].nunique()) if not loo.empty else 0,
        "outputs": {
            "pooled": str(pooled_path),
            "per_tissue_locked": str(locked_path),
            "per_tissue_loo": str(loo_path),
        },
    }
    summary_path = output_dir / "split_plan_summary.json"
    summary["outputs"]["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata",
        default="outputs/generative/benchmark/data_audit/osdr/osdr_canonical_metadata.tsv",
    )
    parser.add_argument(
        "--inventory",
        default="outputs/generative/benchmark/data_audit/osdr/osdr_tissue_inventory.tsv",
    )
    parser.add_argument(
        "--output-dir", default="outputs/generative/benchmark/splits"
    )
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--validation-fraction", type=float, default=0.15)
    parser.add_argument("--test-fraction", type=float, default=0.15)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
