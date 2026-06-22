"""GLARE-style separate-condition fine-tuning for one controlled OSDR study."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader

from .io import dense_matrix, load_matrix_bundle
from .reproduce_glare_finetune import Adapter, SparseAutoEncoder, load_state_dict


CONDITIONS = {"FLT": "flight", "GC": "ground_control"}
COHORT_ORDER = {
    ("ISS-T", "OLD"): 0,
    ("ISS-T", "YNG"): 1,
    ("LAR", "OLD"): 2,
    ("LAR", "YNG"): 3,
}


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def format_elapsed(seconds: float) -> str:
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes}m {seconds}s"
    if minutes:
        return f"{minutes}m {seconds}s"
    return f"{seconds}s"


def profile_cohort(profile: str) -> tuple[str, str, int]:
    """Return preservation cohort, age, and animal number from RR-8 names."""
    platform = "ISS-T" if "ISS-T" in profile else "LAR" if "LAR" in profile else ""
    age = "OLD" if "_OLD_" in profile else "YNG" if "_YNG_" in profile else ""
    match = re.search(r"(\d+)$", profile)
    if not platform or not age or not match:
        raise ValueError(f"Cannot derive OSD-379 cohort slot from profile: {profile}")
    return platform, age, int(match.group(1))


def _ordered_condition_rows(metadata: pd.DataFrame, condition: str) -> pd.DataFrame:
    selected = metadata.loc[metadata["condition_inferred"].eq(condition)].copy()
    if selected.empty:
        raise ValueError(f"No profiles found for condition_inferred={condition!r}")
    cohorts = selected["profile"].astype(str).map(profile_cohort)
    selected["cohort"] = cohorts.map(lambda value: value[0])
    selected["age"] = cohorts.map(lambda value: value[1])
    selected["animal_number"] = cohorts.map(lambda value: value[2])
    selected["cohort_order"] = [
        COHORT_ORDER[(cohort, age)]
        for cohort, age in zip(selected["cohort"], selected["age"])
    ]
    return selected.sort_values(
        ["cohort_order", "animal_number", "profile"]
    ).reset_index(drop=True)


def prepare_controlled_target(
    target_manifest: str | Path,
    accession: str,
    output_dir: Path,
    normalized_counts: str | Path | None = None,
    excluded_profiles: set[str] | None = None,
    filter_mode: str = "matched",
) -> dict:
    """Select balanced FLT/GC profiles and assign comparable cohort slots."""
    bundle = load_matrix_bundle(target_manifest)
    if bundle.profile_metadata is None:
        raise ValueError("Target bundle must include profile metadata")
    metadata = bundle.profile_metadata.copy()
    required = {"profile", "id.accession", "condition_inferred"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Target metadata is missing columns: {sorted(missing)}")

    study_metadata = metadata.loc[metadata["id.accession"].eq(accession)].copy()
    if study_metadata.empty:
        raise ValueError(f"No profiles found for accession {accession}")
    if normalized_counts:
        normalized_counts = Path(normalized_counts)
        counts = pd.read_csv(normalized_counts)
        gene_column = counts.columns[0]
        counts[gene_column] = counts[gene_column].astype(str)
        if counts[gene_column].duplicated().any():
            raise ValueError(f"Duplicate gene IDs in normalized counts: {normalized_counts}")
        counts = counts.set_index(gene_column)
        selected_genes = [gene for gene in bundle.genes if gene in counts.index]
        if not selected_genes:
            raise ValueError("Normalized counts have no genes in common with the bundle")
        target = counts.loc[selected_genes].to_numpy(dtype=np.float32)
        source_profiles = counts.columns.astype(str).tolist()
        input_kind = "official_osdr_normalized_counts"
        input_path = str(normalized_counts)
    else:
        target = dense_matrix(bundle.matrix)
        selected_genes = bundle.genes
        source_profiles = [str(profile) for profile in bundle.profiles]
        input_kind = "matrix_bundle_expression"
        input_path = str(target_manifest)
    profile_to_index = {profile: index for index, profile in enumerate(source_profiles)}

    ordered = {
        label: _ordered_condition_rows(study_metadata, condition)
        for label, condition in CONDITIONS.items()
    }
    cohort_counts = {
        label: rows.groupby(["cohort", "age"]).size().to_dict()
        for label, rows in ordered.items()
    }
    if cohort_counts["FLT"] != cohort_counts["GC"]:
        raise ValueError(
            "FLT and GC cohort counts differ; matched melted features cannot be built: "
            f"{cohort_counts}"
        )

    feature_rows = []
    for (cohort, age), count in sorted(
        cohort_counts["FLT"].items(), key=lambda item: COHORT_ORDER[item[0]]
    ):
        flt_rows = ordered["FLT"].loc[
            ordered["FLT"]["cohort"].eq(cohort) & ordered["FLT"]["age"].eq(age)
        ]
        gc_rows = ordered["GC"].loc[
            ordered["GC"]["cohort"].eq(cohort) & ordered["GC"]["age"].eq(age)
        ]
        for rank, (flt_row, gc_row) in enumerate(
            zip(flt_rows.itertuples(), gc_rows.itertuples()), start=1
        ):
            feature_rows.append(
                {
                    "feature_index": len(feature_rows),
                    "feature": f"{cohort}_{age}_rep{rank:02d}",
                    "cohort": cohort,
                    "age": age,
                    "flt_profile": flt_row.profile,
                    "gc_profile": gc_row.profile,
                    "flt_animal_number": flt_row.animal_number,
                    "gc_animal_number": gc_row.animal_number,
                }
            )

    feature_map = pd.DataFrame(feature_rows)
    excluded_profiles = excluded_profiles or set()
    selected_profiles = set(feature_map["flt_profile"]) | set(feature_map["gc_profile"])
    unknown_exclusions = excluded_profiles - selected_profiles
    if unknown_exclusions:
        raise ValueError(
            "Excluded profiles are not part of the selected FLT/GC study cohort: "
            f"{sorted(unknown_exclusions)}"
        )
    affected_slot_mask = (
        feature_map["flt_profile"].isin(excluded_profiles)
        | feature_map["gc_profile"].isin(excluded_profiles)
    )
    excluded_slots = feature_map.loc[affected_slot_mask].copy()
    feature_map["flt_included"] = ~feature_map["flt_profile"].isin(
        excluded_profiles
    )
    feature_map["gc_included"] = ~feature_map["gc_profile"].isin(excluded_profiles)
    if filter_mode == "matched":
        retained_feature_map = feature_map.loc[~affected_slot_mask].copy()
        profile_maps = {
            label: retained_feature_map[
                [
                    "feature",
                    "cohort",
                    "age",
                    f"{label.lower()}_profile",
                    f"{label.lower()}_animal_number",
                ]
            ].rename(
                columns={
                    f"{label.lower()}_profile": "profile",
                    f"{label.lower()}_animal_number": "animal_number",
                }
            )
            for label in CONDITIONS
        }
    elif filter_mode == "independent":
        retained_feature_map = feature_map.loc[
            feature_map["flt_included"] & feature_map["gc_included"]
        ].copy()
        profile_maps = {}
        for label in CONDITIONS:
            included_column = f"{label.lower()}_included"
            profile_column = f"{label.lower()}_profile"
            animal_column = f"{label.lower()}_animal_number"
            profile_maps[label] = feature_map.loc[
                feature_map[included_column],
                ["feature", "cohort", "age", profile_column, animal_column],
            ].rename(
                columns={
                    profile_column: "profile",
                    animal_column: "animal_number",
                }
            )
    else:
        raise ValueError(f"Unsupported filter_mode: {filter_mode}")

    retained_feature_map["original_feature_index"] = retained_feature_map[
        "feature_index"
    ]
    retained_feature_map["feature_index"] = np.arange(len(retained_feature_map))
    retained_feature_map = retained_feature_map.reset_index(drop=True)

    matrices = {}
    for label in CONDITIONS:
        profiles = profile_maps[label]["profile"].astype(str).tolist()
        missing_profiles = [
            profile for profile in profiles if profile not in profile_to_index
        ]
        if missing_profiles:
            raise ValueError(
                f"Expression input is missing {len(missing_profiles)} {label} profiles; "
                f"first missing: {missing_profiles[0]}"
            )
        indices = [profile_to_index[profile] for profile in profiles]
        matrices[label] = target[:, indices].astype(np.float32, copy=False)

    filtered_cohort_counts = {
        label: {
            f"{cohort}_{age}": int(count)
            for (cohort, age), count in profile_maps[label].groupby(
                ["cohort", "age"]
            ).size().items()
        }
        for label in CONDITIONS
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    target_path = output_dir / "controlled_target.npz"
    np.savez_compressed(
        target_path,
        flt=matrices["FLT"],
        gc=matrices["GC"],
        genes=np.asarray(selected_genes, dtype=str),
        features=retained_feature_map["feature"].to_numpy(dtype=str),
        flt_features=profile_maps["FLT"]["feature"].to_numpy(dtype=str),
        gc_features=profile_maps["GC"]["feature"].to_numpy(dtype=str),
        input_kind=np.asarray(input_kind),
        input_path=np.asarray(input_path),
    )
    retained_feature_map.to_csv(
        output_dir / "matched_feature_slots.tsv", sep="\t", index=False
    )
    excluded_slots.to_csv(
        output_dir / "excluded_matched_feature_slots.tsv", sep="\t", index=False
    )
    feature_map.to_csv(
        output_dir / "all_feature_slots_with_filter.tsv", sep="\t", index=False
    )
    pd.concat(
        [
            profile_maps[label].assign(location=label)
            for label in CONDITIONS
        ],
        ignore_index=True,
    ).to_csv(output_dir / "retained_profile_features.tsv", sep="\t", index=False)
    study_metadata.to_csv(output_dir / "study_profile_metadata.tsv", sep="\t", index=False)
    return {
        "path": target_path,
        "genes": selected_genes,
        "genes_dropped_from_broad_alignment": len(bundle.genes) - len(selected_genes),
        "features": {
            label: profile_maps[label]["feature"].tolist()
            for label in CONDITIONS
        },
        "matrices": matrices,
        "cohort_counts": filtered_cohort_counts,
        "study_profile_count": int(len(study_metadata)),
        "input_kind": input_kind,
        "input_path": input_path,
        "excluded_profiles_requested": sorted(excluded_profiles),
        "affected_matched_slots": int(len(excluded_slots)),
        "excluded_matched_slots": (
            int(len(excluded_slots)) if filter_mode == "matched" else 0
        ),
        "excluded_profiles_total": (
            int(2 * len(excluded_slots))
            if filter_mode == "matched"
            else int(len(excluded_profiles))
        ),
        "filter_mode": filter_mode,
    }


def load_excluded_profiles(path: str | Path | None) -> set[str]:
    if not path:
        return set()
    path = Path(path)
    table = pd.read_csv(path, sep=None, engine="python")
    for column in ("sample", "directly_flagged_profile"):
        if column in table.columns:
            return set(table[column].dropna().astype(str))
    if table.shape[1] == 1:
        return set(table.iloc[:, 0].dropna().astype(str))
    raise ValueError(
        f"Exclusion table must contain sample or directly_flagged_profile: {path}"
    )


def write_epoch_logs(output_dir: Path, location: str, records: list[dict]) -> None:
    fields = [
        "location",
        "epoch",
        "loss",
        "best_loss",
        "elapsed_seconds",
        "elapsed",
        "epoch_seconds",
    ]
    with (output_dir / f"{location}_epoch_losses.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    (output_dir / f"{location}_epoch_losses.json").write_text(
        json.dumps(records, indent=2) + "\n", encoding="utf-8"
    )


def write_outlier_audit(
    matrix: np.ndarray,
    genes: list[str],
    location: str,
    output_dir: Path,
) -> None:
    """Export GLARE's PCA/k-means inspection without arbitrary gene removal."""
    scaled = StandardScaler().fit_transform(matrix)
    pca = PCA(n_components=3, random_state=1).fit_transform(scaled)
    n_clusters = 5 if location == "FLT" else 4
    labels = KMeans(
        n_clusters=n_clusters,
        init="k-means++",
        n_init=10,
        random_state=1,
    ).fit_predict(pca)
    pd.DataFrame(
        {
            "gene_id": genes,
            "pc1": pca[:, 0],
            "pc2": pca[:, 1],
            "pc3": pca[:, 2],
            "inspection_cluster": labels,
        }
    ).to_csv(output_dir / f"{location}_outlier_audit.tsv", sep="\t", index=False)


def infer_pretrain_input_dim(weights: Path) -> int:
    state = load_state_dict(weights)
    key = "encoder.0.weight"
    if key not in state:
        raise ValueError(f"Cannot infer pretraining input dimension from {weights}")
    return int(state[key].shape[1])


def finetune_location(
    matrix: np.ndarray,
    genes: list[str],
    location: str,
    pretrained_weights: Path,
    output_dir: Path,
    device: torch.device,
    input_dim: int,
    epochs: int,
    batch_size: int,
    seed: int,
) -> dict:
    start = time.perf_counter()
    torch.manual_seed(seed)
    np.random.seed(seed)

    scaler = StandardScaler()
    scaled = scaler.fit_transform(matrix).astype(np.float32, copy=False)
    np.savez_compressed(
        output_dir / f"{location}_standard_scaler.npz",
        mean=scaler.mean_,
        scale=scaler.scale_,
    )
    source = torch.tensor(scaled, dtype=torch.float32)
    adapter = Adapter(source.shape[1], input_dim)
    adapted = adapter(source).clone().detach()
    loader = DataLoader(adapted, batch_size=batch_size, shuffle=True, num_workers=0)

    model = SparseAutoEncoder(input_dim).to(device)
    model.load_state_dict(load_state_dict(pretrained_weights))
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    records = []
    best_loss = float("inf")
    best_epoch = 0
    best_state = None
    log(
        f"{location}: fine-tuning {matrix.shape[0]} genes x {matrix.shape[1]} "
        f"profiles through adapter -> {input_dim}"
    )
    for epoch in range(epochs):
        epoch_start = time.perf_counter()
        total_loss = 0.0
        model.train()
        for batch in loader:
            batch = batch.to(device)
            optimizer.zero_grad()
            outputs = model(batch)
            loss = criterion(outputs, batch)
            # This intentionally matches GLARE's released sparsity calculation.
            encoded_for_l1 = model.encoder[-1](batch)
            loss += 1e-5 * torch.mean(torch.abs(encoded_for_l1))
            loss.backward()
            optimizer.step()
            total_loss += float(loss.item())

        average_loss = total_loss / len(loader)
        if average_loss < best_loss:
            best_loss = average_loss
            best_epoch = epoch + 1
            best_state = copy.deepcopy(model.state_dict())
        elapsed = time.perf_counter() - start
        records.append(
            {
                "location": location,
                "epoch": epoch + 1,
                "loss": round(average_loss, 8),
                "best_loss": round(best_loss, 8),
                "elapsed_seconds": round(elapsed, 3),
                "elapsed": format_elapsed(elapsed),
                "epoch_seconds": round(time.perf_counter() - epoch_start, 3),
            }
        )
        write_epoch_logs(output_dir, location, records)
        log(
            f"{location}: epoch {epoch + 1}/{epochs} loss={average_loss:.8f} "
            f"best={best_loss:.8f}"
        )

    if best_state is None:
        raise AssertionError("Fine-tuning completed without a best model")
    model.load_state_dict(best_state)
    weights_path = output_dir / f"{location}_finetuned_sae.pth"
    adapter_path = output_dir / f"{location}_adapter.pth"
    torch.save(best_state, weights_path)
    torch.save(adapter.state_dict(), adapter_path)
    model.eval()
    with torch.no_grad():
        representation = model.encoder(adapted.to(device)).cpu().numpy()
    representation_path = output_dir / f"{location}_FTSAE_representation.npy"
    np.save(representation_path, representation)
    pd.DataFrame(
        representation,
        index=genes,
        columns=[f"latent_{index + 1}" for index in range(representation.shape[1])],
    ).rename_axis("gene_id").to_csv(
        output_dir / f"{location}_gene_latent.tsv", sep="\t"
    )
    return {
        "location": location,
        "genes": int(matrix.shape[0]),
        "profiles": int(matrix.shape[1]),
        "best_loss": round(best_loss, 8),
        "best_epoch": best_epoch,
        "epochs": epochs,
        "weights": str(weights_path),
        "adapter": str(adapter_path),
        "representation": str(representation_path),
        "elapsed_seconds": round(time.perf_counter() - start, 3),
        "elapsed": format_elapsed(time.perf_counter() - start),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fine-tune GLARE separately on FLT and GC from one OSDR study."
    )
    parser.add_argument(
        "--target-manifest",
        default="data/processed/tms_facs_liver_osdr_liver_aligned.target.manifest.json",
    )
    parser.add_argument("--accession", default="OSD-379")
    parser.add_argument(
        "--normalized-counts",
        help="Official OSDR normalized-counts CSV; preferred over raw bundle expression.",
    )
    parser.add_argument("--pretrained-weights", required=True)
    parser.add_argument(
        "--output-dir",
        default="outputs/glare_paper_tms_liver_osd379",
    )
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1996)
    parser.add_argument(
        "--exclude-samples",
        help=(
            "CSV/TSV containing sample or directly_flagged_profile. Removal "
            "behavior is controlled by --filter-mode."
        ),
    )
    parser.add_argument(
        "--filter-mode",
        choices=["matched", "independent"],
        default="matched",
        help=(
            "matched removes both profiles in affected slots; independent removes "
            "only directly flagged profiles and permits unequal FLT/GC dimensions."
        ),
    )
    parser.add_argument("--prepare-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_start = time.perf_counter()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    excluded_profiles = load_excluded_profiles(args.exclude_samples)
    prepared = prepare_controlled_target(
        args.target_manifest,
        args.accession,
        output_dir,
        args.normalized_counts,
        excluded_profiles,
        args.filter_mode,
    )
    log(
        f"Prepared {args.accession}: {len(prepared['genes'])} genes, "
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
        for location in CONDITIONS
    ]
    summary = {
        "method": "GLARE released 16-dimensional SAE with separate FLT/GC fine-tuning",
        "accession": args.accession,
        "target_manifest": args.target_manifest,
        "target_expression_input": prepared["input_path"],
        "target_expression_kind": prepared["input_kind"],
        "pretrained_weights": str(pretrained_weights),
        "pretrained_input_dim": input_dim,
        "device": str(device),
        "seed_reused_for_each_location": args.seed,
        "architecture": [128, 64, 32, 16],
        "learning_rate": 1e-3,
        "weight_decay": 0,
        "sparsity_penalty": 1e-5,
        "batch_size": args.batch_size,
        "cohort_counts": prepared["cohort_counts"],
        "sample_filter": {
            "exclusion_file": args.exclude_samples or "",
            "directly_flagged_profiles": prepared["excluded_profiles_requested"],
            "affected_matched_slots": prepared["affected_matched_slots"],
            "matched_slots_excluded": prepared["excluded_matched_slots"],
            "total_profiles_excluded": prepared["excluded_profiles_total"],
            "filter_mode": prepared["filter_mode"],
        },
        "study_profile_count_all_conditions": prepared["study_profile_count"],
        "genes_dropped_from_broad_alignment": prepared[
            "genes_dropped_from_broad_alignment"
        ],
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
