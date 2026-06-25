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

from .io import dense_matrix, load_matrix_bundle, require_import
from .paper_finetune import (
    finetune_location,
    format_elapsed,
    infer_pretrain_input_dim,
    log,
    write_outlier_audit,
)


DEFAULT_TARGET_MANIFEST = (
    "data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json"
)
DEFAULT_OSDR_H5 = "assets/osdr/OSDR_mouse_RNAseq_Feb2026.h5"
DEFAULT_OUTPUT_DIR = "outputs/glare_tms_liver_aggregated_osdr_flt_gc"
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
CONDITION_MAP = {
    "space flight": ("FLT", "flight"),
    "ground control": ("GC", "ground_control"),
    "ground control group": ("GC", "ground_control"),
}
H5_FIELDS = {
    "h5_accession": "/meta/info/id.accession",
    "h5_sample_name": "/meta/info/id.sample name",
    "h5_accession_sample_name": "/meta/info/id.accession_sample name",
    "project_identifier": "/meta/info/investigation.study.comment.project identifier",
    "project_type": "/meta/info/investigation.study.comment.project type",
    "assay_technology": (
        "/meta/info/investigation.study assays.study assay technology type"
    ),
    "source_name": "/meta/info/study.source name",
    "material_type": (
        "/meta/samples/characteristics/study.characteristics.material type"
    ),
    "tissue_type": "/meta/samples/characteristics/study.characteristics.tissue type",
    "spaceflight_factor": "/meta/samples/factors/study.factor value.spaceflight",
    "sex": "/meta/samples/characteristics/study.characteristics.sex",
    "strain": "/meta/samples/characteristics/study.characteristics.strain",
    "genotype": "/meta/samples/characteristics/study.characteristics.genotype",
    "age_at_launch": (
        "/meta/samples/characteristics/study.characteristics.age at launch"
    ),
    "age": "/meta/samples/characteristics/study.characteristics.age",
    "duration": "/meta/samples/parameters/study.parameter value.duration",
    "sample_preservation_method": (
        "/meta/samples/parameters/study.parameter value.sample preservation method"
    ),
    "library_selection": (
        "/meta/samples/assay_parameters/assay.parameter value.library selection"
    ),
    "library_layout": (
        "/meta/samples/assay_parameters/assay.parameter value.library layout"
    ),
    "sequencing_instrument": (
        "/meta/samples/assay_parameters/assay.parameter value.sequencing instrument"
    ),
}


def clean_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return token or "profile"


def decode_array(values) -> list[str]:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8", "replace"))
        else:
            decoded.append(str(value))
    return decoded


def normalize_spaceflight(value: str) -> tuple[str, str] | None:
    cleaned = re.sub(r"\s+", " ", str(value).strip().lower())
    return CONDITION_MAP.get(cleaned)


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


def load_h5_profile_metadata(osdr_h5: str | Path, metadata: pd.DataFrame) -> pd.DataFrame:
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")

    with h5py.File(osdr_h5, "r") as handle:
        source_profile_indices = None
        if "source_profile_index" in metadata.columns:
            source_profile_indices = metadata["source_profile_index"].astype(int).tolist()
        else:
            source_profiles = decode_array(handle["/meta/info/id.sample name"][:])
            profile_to_index = {
                profile: index for index, profile in enumerate(source_profiles)
            }
            source_profile_indices = [
                profile_to_index.get(str(profile), -1)
                for profile in metadata["profile"].astype(str)
            ]
            if any(index < 0 for index in source_profile_indices):
                missing = metadata.loc[
                    [index < 0 for index in source_profile_indices], "profile"
                ].head(5)
                raise ValueError(
                    "Could not map bundle profiles to HDF5 profile indices; "
                    f"first missing: {missing.tolist()}"
                )

        rows = {"profile": metadata["profile"].astype(str).tolist()}
        for column, h5_path in H5_FIELDS.items():
            if h5_path not in handle:
                rows[column] = [""] * len(source_profile_indices)
                continue
            values = np.asarray(decode_array(handle[h5_path][:]), dtype=object)
            rows[column] = values[source_profile_indices].astype(str).tolist()

    h5_metadata = pd.DataFrame(rows)
    merged = metadata.reset_index(drop=True).copy()
    for column in h5_metadata.columns:
        if column == "profile":
            continue
        if column in merged.columns:
            merged[column] = merged[column].fillna(h5_metadata[column])
        else:
            merged[column] = h5_metadata[column]
    return merged


def select_aggregate_profiles(
    target_manifest: str | Path,
    osdr_h5: str | Path,
    accessions: list[str],
    output_dir: Path,
    exclude_profiles: set[str] | None = None,
    ercc_policy: str = "keep_all",
) -> dict:
    bundle = load_matrix_bundle(target_manifest)
    if bundle.profile_metadata is None:
        metadata = pd.DataFrame({"profile": bundle.profiles})
    else:
        metadata = bundle.profile_metadata.copy()
    if "profile" not in metadata.columns:
        metadata.insert(0, "profile", bundle.profiles)
    metadata = metadata.reset_index(drop=True)
    metadata["profile"] = [str(profile) for profile in bundle.profiles]

    metadata = load_h5_profile_metadata(osdr_h5, metadata)
    metadata["condition_label"] = ""
    metadata["condition_inferred"] = ""
    for index, value in metadata["spaceflight_factor"].items():
        normalized = normalize_spaceflight(value)
        if normalized is None:
            continue
        metadata.at[index, "condition_label"] = normalized[0]
        metadata.at[index, "condition_inferred"] = normalized[1]

    accession_set = set(accessions)
    selected_mask = (
        metadata["h5_accession"].isin(accession_set)
        & metadata["material_type"].map(is_liver_material)
        & metadata["condition_label"].isin(["FLT", "GC"])
    )
    selected = metadata.loc[selected_mask].copy()
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

    missing_accessions = sorted(accession_set - set(selected["h5_accession"]))
    if missing_accessions:
        raise ValueError(
            "No selected FLT/GC liver profiles for accessions: "
            f"{missing_accessions}"
        )

    selected, ercc_summary = apply_ercc_policy(selected, output_dir, ercc_policy)
    missing_accessions = sorted(accession_set - set(selected["h5_accession"]))
    if missing_accessions:
        raise ValueError(
            "No selected FLT/GC liver profiles remain after ERCC filtering for "
            f"accessions: {missing_accessions}"
        )

    selected["accession_order"] = selected["h5_accession"].map(
        {accession: index for index, accession in enumerate(accessions)}
    )
    selected = selected.sort_values(
        ["condition_label", "accession_order", "profile"]
    ).reset_index(drop=True)

    matrix = dense_matrix(bundle.matrix)
    profile_to_index = {str(profile): index for index, profile in enumerate(bundle.profiles)}
    genes = [str(gene) for gene in bundle.genes]
    matrices = {}
    features = {}
    retained_rows = []
    for location in ("FLT", "GC"):
        rows = selected.loc[selected["condition_label"].eq(location)].copy()
        profiles = rows["profile"].astype(str).tolist()
        indices = [profile_to_index[profile] for profile in profiles]
        matrices[location] = matrix[:, indices].astype(np.float32, copy=False)
        features[location] = [
            f"{row.h5_accession}_{clean_token(row.profile)}"
            for row in rows.itertuples()
        ]
        retained = rows.copy()
        retained.insert(0, "location", location)
        retained.insert(1, "feature", features[location])
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
        input_kind=np.asarray("aligned_osdr_liver_hdf5_expression"),
        input_path=np.asarray(str(target_manifest)),
    )

    retained_profile_features = pd.concat(retained_rows, ignore_index=True)
    retained_profile_features.to_csv(
        output_dir / "retained_profile_features.tsv", sep="\t", index=False
    )
    selected.to_csv(output_dir / "study_profile_metadata.tsv", sep="\t", index=False)
    counts = (
        selected.groupby(["h5_accession", "condition_label"])
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
        "metadata": selected,
        "counts": counts,
        "input_path": str(target_manifest),
        "target_manifest": str(target_manifest),
        "osdr_h5": str(osdr_h5),
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
    parser.add_argument("--target-manifest", default=DEFAULT_TARGET_MANIFEST)
    parser.add_argument("--osdr-h5", default=DEFAULT_OSDR_H5)
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
        args.target_manifest,
        args.osdr_h5,
        args.accessions,
        output_dir,
        load_excluded_profiles(args.exclude_profiles_file, args.exclude_profile),
        args.ercc_policy,
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
            "material_field": (
                "/meta/samples/characteristics/study.characteristics.material type"
            ),
            "condition_field": "/meta/samples/factors/study.factor value.spaceflight",
            "included_conditions": ["Space Flight", "Ground Control"],
            "excluded_profiles_requested": prepared["excluded_profiles_requested"],
            "excluded_profiles_matched": prepared["excluded_profiles_matched"],
            "ercc_policy": prepared["ercc_policy"],
        },
        "condition_counts": counts.to_dict(orient="records"),
        "target_manifest": prepared["target_manifest"],
        "target_expression_input": prepared["input_path"],
        "target_expression_kind": "aligned_osdr_liver_hdf5_expression",
        "osdr_h5": prepared["osdr_h5"],
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
