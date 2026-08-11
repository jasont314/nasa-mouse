"""Fine-tune GLARE on aggregated OSDR liver Space Flight vs Ground Control."""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from .osdr import (
    DEFAULT_API_METADATA,
    DEFAULT_COUNTS_DIR,
    load_api_expression,
    read_api_metadata,
    select_api_metadata,
)
from .paper_finetune import (
    finetune_location,
    format_elapsed,
    infer_pretrain_input_dim,
    log,
    write_outlier_audit,
)


DEFAULT_OUTPUT_DIR = "outputs/glare/tms_liver_aggregated_osdr_flt_gc"
DEFAULT_ACCESSIONS = [
    "OSD-379",
    "OSD-245",
    "OSD-463",
    "OSD-242",
    "OSD-137",
    "OSD-47",
    "OSD-686",
    "OSD-173",
]


def clean_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return token or "profile"


def load_excluded_profiles(
    exclude_profiles_file: str | Path | None = None,
    exclude_profile: list[str] | None = None,
) -> set[str]:
    excluded: set[str] = set()
    if exclude_profiles_file:
        path = Path(exclude_profiles_file)
        for line in path.read_text(encoding="utf-8").splitlines():
            token = line.strip()
            if token and not token.startswith("#"):
                excluded.add(token)
    for token in exclude_profile or []:
        token = str(token).strip()
        if token:
            excluded.add(token)
    return excluded


def is_liver_material(value: str) -> bool:
    return "liver" in str(value).strip().lower()


def ercc_status(profile: str) -> str:
    profile = str(profile)
    if "_noERCC_" in profile or profile.endswith("_noERCC"):
        return "noERCC"
    if "_wERCC_" in profile or profile.endswith("_wERCC"):
        return "wERCC"
    return "not_annotated"


def ercc_biological_key(profile: str) -> str:
    return re.sub(r"_(?:wERCC|noERCC)(?=_|$)", "", str(profile))


def apply_ercc_policy(
    selected: pd.DataFrame,
    output_dir: Path,
    ercc_policy: str,
) -> tuple[pd.DataFrame, dict]:
    if ercc_policy == "keep_all":
        return selected, {
            "policy": ercc_policy,
            "profiles_before": int(len(selected)),
            "profiles_after": int(len(selected)),
            "profiles_dropped": 0,
            "duplicate_groups": 0,
            "unique_wERCC_profiles_retained": int(
                selected["profile"].astype(str).map(ercc_status).eq("wERCC").sum()
            ),
        }
    if ercc_policy != "prefer_noercc":
        raise ValueError(f"Unsupported ERCC policy: {ercc_policy}")

    output_dir.mkdir(parents=True, exist_ok=True)
    working = selected.copy()
    working["ercc_status"] = working["profile"].astype(str).map(ercc_status)
    working["ercc_biological_key"] = working["profile"].astype(str).map(
        ercc_biological_key
    )
    group_columns = ["h5_accession", "condition_label", "ercc_biological_key"]
    audit_rows = []
    retained_indices = []
    for group_key, group in working.groupby(group_columns, dropna=False, sort=False):
        statuses = set(group["ercc_status"])
        duplicate_group = len(group) > 1
        if "noERCC" in statuses:
            keep = group.loc[group["ercc_status"].eq("noERCC")].sort_values(
                "profile"
            ).head(1)
            decision = "kept_noERCC"
        else:
            keep = group.sort_values("profile").head(1)
            decision = (
                "kept_unique_wERCC"
                if "wERCC" in statuses
                else "kept_unannotated"
            )
        retained_indices.extend(keep.index.tolist())
        kept_profiles = set(keep["profile"].astype(str))
        for row in group.itertuples():
            audit_rows.append(
                {
                    "h5_accession": group_key[0],
                    "condition_label": group_key[1],
                    "ercc_biological_key": group_key[2],
                    "profile": row.profile,
                    "ercc_status": row.ercc_status,
                    "duplicate_group": duplicate_group,
                    "retained": str(row.profile) in kept_profiles,
                    "decision": decision,
                }
            )

    audit = pd.DataFrame(audit_rows)
    audit.to_csv(output_dir / "ercc_profile_policy.tsv", sep="\t", index=False)
    retained = working.loc[retained_indices].copy()
    retained = retained.drop(columns=["ercc_status", "ercc_biological_key"])
    dropped = audit.loc[~audit["retained"]]
    unique_w_ercc = audit.loc[
        audit["retained"]
        & audit["ercc_status"].eq("wERCC")
        & ~audit["duplicate_group"]
    ]
    summary = {
        "policy": ercc_policy,
        "profiles_before": int(len(selected)),
        "profiles_after": int(len(retained)),
        "profiles_dropped": int(len(dropped)),
        "duplicate_groups": int(audit.loc[audit["duplicate_group"], group_columns].drop_duplicates().shape[0]),
        "unique_wERCC_profiles_retained": int(len(unique_w_ercc)),
        "audit_path": str(output_dir / "ercc_profile_policy.tsv"),
    }
    return retained, summary


def select_aggregate_profiles(
    api_metadata: str | Path,
    counts_dir: str | Path,
    accessions: list[str],
    output_dir: Path,
    exclude_profiles: set[str] | None = None,
    ercc_policy: str = "keep_all",
    *,
    refresh_metadata: bool = False,
    download_counts: bool = False,
    timeout: int = 180,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_api_metadata(
        api_metadata,
        refresh=refresh_metadata,
        timeout=timeout,
    )
    selected = select_api_metadata(
        metadata,
        tissue="liver",
        accessions=accessions,
    )
    accession_set = set(accessions)
    selected = selected.loc[selected["material_type"].map(is_liver_material)].copy()
    if selected.empty:
        raise ValueError("No aggregate FLT/GC liver profiles matched the selection")

    excluded_requested = set(exclude_profiles or set())
    excluded_selected = pd.DataFrame(columns=selected.columns)
    if excluded_requested:
        exclude_columns = [
            column
            for column in [
                "profile",
                "h5_sample_name",
                "h5_accession_sample_name",
                "official_sample_name",
            ]
            if column in selected.columns
        ]
        exclude_mask = pd.Series(False, index=selected.index)
        for column in exclude_columns:
            exclude_mask |= selected[column].astype(str).isin(excluded_requested)
        excluded_selected = selected.loc[exclude_mask].copy()
        selected = selected.loc[~exclude_mask].copy()
        output_dir.mkdir(parents=True, exist_ok=True)
        excluded_selected.to_csv(
            output_dir / "excluded_profile_features.tsv", sep="\t", index=False
        )
        matched = set()
        for column in exclude_columns:
            matched.update(excluded_selected[column].dropna().astype(str).tolist())
        unmatched = sorted(excluded_requested - matched)
        if unmatched:
            (output_dir / "unmatched_excluded_profiles.txt").write_text(
                "\n".join(unmatched) + "\n", encoding="utf-8"
            )
        if selected.empty:
            raise ValueError("All selected aggregate FLT/GC liver profiles were excluded")

    missing_accessions = sorted(accession_set - set(selected["id.accession"]))
    if missing_accessions:
        raise ValueError(
            "No selected FLT/GC liver profiles for accessions: "
            f"{missing_accessions}"
        )

    selected, ercc_summary = apply_ercc_policy(selected, output_dir, ercc_policy)
    missing_accessions = sorted(accession_set - set(selected["id.accession"]))
    if missing_accessions:
        raise ValueError(
            "No selected FLT/GC liver profiles remain after ERCC filtering for "
            f"accessions: {missing_accessions}"
        )

    selected["accession_order"] = selected["id.accession"].map(
        {accession: index for index, accession in enumerate(accessions)}
    )
    selected = selected.sort_values(
        ["condition_label", "accession_order", "profile"]
    ).reset_index(drop=True)

    loaded = load_api_expression(
        selected,
        counts_dir=counts_dir,
        timeout=timeout,
        download_missing=download_counts,
    )
    raw = loaded["raw_counts"]
    retained_metadata = loaded["metadata"]
    missing_count_columns = loaded["missing_count_columns"]
    assert isinstance(raw, pd.DataFrame)
    assert isinstance(retained_metadata, pd.DataFrame)
    assert isinstance(missing_count_columns, pd.DataFrame)
    if not missing_count_columns.empty:
        missing_count_columns.to_csv(
            output_dir / "missing_count_columns.tsv", sep="\t", index=False
        )

    genes = raw.index.astype(str).tolist()
    matrices = {}
    features = {}
    retained_rows = []
    for location in ("FLT", "GC"):
        rows = retained_metadata.loc[
            retained_metadata["condition_label"].eq(location)
        ].copy()
        features[location] = rows["feature"].astype(str).tolist()
        matrices[location] = raw.loc[:, features[location]].to_numpy(
            dtype=np.float32, copy=False
        )
        retained = rows.copy()
        retained.insert(0, "location", location)
        retained_rows.append(retained)

    output_dir.mkdir(parents=True, exist_ok=True)
    target_path = output_dir / "controlled_target.npz"
    np.savez_compressed(
        target_path,
        flt=matrices["FLT"],
        gc=matrices["GC"],
        genes=np.asarray(genes, dtype=str),
        flt_features=np.asarray(features["FLT"], dtype=str),
        gc_features=np.asarray(features["GC"], dtype=str),
        input_kind=np.asarray("nasa_osdr_api_unnormalized_counts"),
        input_path=np.asarray(str(api_metadata)),
    )

    retained_profile_features = pd.concat(retained_rows, ignore_index=True)
    retained_profile_features.to_csv(
        output_dir / "retained_profile_features.tsv", sep="\t", index=False
    )
    retained_metadata.to_csv(
        output_dir / "study_profile_metadata.tsv", sep="\t", index=False
    )
    counts = (
        retained_metadata.groupby(["id.accession", "condition_label"])
        .size()
        .unstack(fill_value=0)
        .reindex(accessions)
        .fillna(0)
        .astype(int)
    )
    counts["total"] = counts.sum(axis=1)
    counts.to_csv(output_dir / "aggregate_condition_counts.tsv", sep="\t")

    return {
        "path": target_path,
        "genes": genes,
        "matrices": matrices,
        "features": features,
        "metadata": retained_metadata,
        "counts": counts,
        "input_path": str(api_metadata),
        "api_metadata": str(api_metadata),
        "counts_dir": str(counts_dir),
        "excluded_profiles_requested": sorted(excluded_requested),
        "excluded_profiles_matched": int(len(excluded_selected)),
        "ercc_policy": ercc_summary,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fine-tune GLARE separately on aggregated OSDR liver Space Flight "
            "and Ground Control profiles."
        )
    )
    parser.add_argument("--api-metadata", default=DEFAULT_API_METADATA)
    parser.add_argument("--counts-dir", default=DEFAULT_COUNTS_DIR)
    parser.add_argument("--refresh-metadata", action="store_true")
    parser.add_argument("--download-counts", action="store_true")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument(
        "--accessions",
        nargs="+",
        default=DEFAULT_ACCESSIONS,
        help="OSD accessions to include, in output/report order.",
    )
    parser.add_argument("--pretrained-weights", required=True)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1996)
    parser.add_argument(
        "--exclude-profiles-file",
        help="Text file of profile/sample IDs to exclude, one per line.",
    )
    parser.add_argument(
        "--exclude-profile",
        action="append",
        default=[],
        help="Profile/sample ID to exclude. Can be supplied multiple times.",
    )
    parser.add_argument(
        "--ercc-policy",
        choices=["keep_all", "prefer_noercc"],
        default="keep_all",
        help=(
            "prefer_noercc collapses wERCC/noERCC duplicate profiles by keeping "
            "noERCC when both are present."
        ),
    )
    parser.add_argument("--prepare-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_start = time.perf_counter()
    output_dir = Path(args.output_dir)
    prepared = select_aggregate_profiles(
        args.api_metadata,
        args.counts_dir,
        args.accessions,
        output_dir,
        load_excluded_profiles(args.exclude_profiles_file, args.exclude_profile),
        args.ercc_policy,
        refresh_metadata=args.refresh_metadata,
        download_counts=args.download_counts,
        timeout=args.timeout,
    )
    log(
        "Prepared aggregate liver FLT/GC target: "
        f"{len(prepared['genes'])} genes, "
        f"{prepared['matrices']['FLT'].shape[1]} FLT and "
        f"{prepared['matrices']['GC'].shape[1]} GC profiles"
    )
    if args.prepare_only:
        return

    pretrained_weights = Path(args.pretrained_weights)
    input_dim = infer_pretrain_input_dim(pretrained_weights)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    for location, matrix in prepared["matrices"].items():
        write_outlier_audit(matrix, prepared["genes"], location, output_dir)

    locations = [
        finetune_location(
            prepared["matrices"][location],
            prepared["genes"],
            location,
            pretrained_weights,
            output_dir,
            device,
            input_dim,
            args.epochs,
            args.batch_size,
            args.seed,
        )
        for location in ("FLT", "GC")
    ]
    counts = prepared["counts"].reset_index().rename(columns={"index": "accession"})
    summary = {
        "method": "GLARE released 16-dimensional SAE with aggregated liver FLT/GC fine-tuning",
        "accessions": args.accessions,
        "selection": {
            "material_field": "study.characteristics.material type",
            "condition_field": "condition_inferred",
            "included_conditions": ["flight", "ground_control"],
            "excluded_profiles_requested": prepared["excluded_profiles_requested"],
            "excluded_profiles_matched": prepared["excluded_profiles_matched"],
            "ercc_policy": prepared["ercc_policy"],
        },
        "condition_counts": counts.to_dict(orient="records"),
        "api_metadata": prepared["api_metadata"],
        "counts_dir": prepared["counts_dir"],
        "target_expression_input": prepared["input_path"],
        "target_expression_kind": "nasa_osdr_api_unnormalized_counts",
        "pretrained_weights": str(pretrained_weights),
        "pretrained_input_dim": input_dim,
        "device": str(device),
        "seed_reused_for_each_location": args.seed,
        "architecture": [128, 64, 32, 16],
        "learning_rate": 1e-3,
        "weight_decay": 0,
        "sparsity_penalty": 1e-5,
        "batch_size": args.batch_size,
        "outlier_policy": (
            "PCA/k-means audit exported; no genes removed because GLARE's three "
            "fixed Arabidopsis outlier IDs do not transfer to mouse"
        ),
        "locations": locations,
        "elapsed_seconds": round(time.perf_counter() - run_start, 3),
        "elapsed": format_elapsed(time.perf_counter() - run_start),
    }
    (output_dir / "finetune_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    log(f"Saved summary: {output_dir / 'finetune_summary.json'}")


if __name__ == "__main__":
    main()
