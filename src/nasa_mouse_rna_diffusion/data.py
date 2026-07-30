"""Prepare the ARCHS4 mouse counterpart to the paper's GTEx DDIM data."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time
from typing import Any

import h5py
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, train_test_split

from .config import load_config
from .landmarks import build_mouse_landmark_panel


def _decode(values: np.ndarray) -> list[str]:
    return [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _balanced_candidates(
    metadata: pd.DataFrame, *, minimum: int, seed: int
) -> tuple[pd.DataFrame, dict[str, int]]:
    counts = metadata.groupby("canonical_tissue").size()
    eligible = sorted(counts[counts >= int(minimum)].index.astype(str))
    metadata = metadata.loc[
        metadata["canonical_tissue"].astype(str).isin(eligible)
    ].copy()
    rng = np.random.default_rng(int(seed))
    groups: dict[str, pd.DataFrame] = {}
    for tissue in eligible:
        group = metadata.loc[metadata["canonical_tissue"].astype(str).eq(tissue)]
        order = rng.permutation(len(group))
        groups[tissue] = group.iloc[order].reset_index(drop=True)
    rows: list[pd.Series] = []
    round_index = 0
    while True:
        added = False
        for tissue in eligible:
            group = groups[tissue]
            if round_index < len(group):
                rows.append(group.iloc[round_index])
                added = True
        if not added:
            break
        round_index += 1
    result = pd.DataFrame(rows).reset_index(drop=True)
    return result, {tissue: int(len(groups[tissue])) for tissue in eligible}


def _targeted_candidates(
    metadata: pd.DataFrame, *, quotas: dict[str, int], seed: int
) -> tuple[pd.DataFrame, dict[str, int]]:
    available = {
        tissue: int(count)
        for tissue, count in metadata.groupby("canonical_tissue").size().items()
    }
    selected: list[pd.DataFrame] = []
    rng = np.random.default_rng(int(seed))
    for tissue, quota in sorted(quotas.items()):
        group = metadata.loc[
            metadata["canonical_tissue"].astype(str).eq(str(tissue))
        ].copy()
        if len(group) < int(quota):
            raise ValueError(
                f"ARCHS4 has {len(group)} selected {tissue} profiles; need {quota}"
            )
        if "selection_rank_within_tissue" in group:
            group = group.sort_values(
                ["selection_rank_within_tissue", "archs4_sample_index"],
                kind="stable",
            )
        else:
            group = group.iloc[rng.permutation(len(group))]
        selected.append(group.iloc[: int(quota)])
    result = pd.concat(selected, ignore_index=True)
    result = result.iloc[rng.permutation(len(result))].reset_index(drop=True)
    return result, available


def _length_vector(archs4_genes: list[str], path: str | Path) -> np.ndarray:
    table = pd.read_csv(path, sep="\t")
    lengths = table.set_index("gene_id")["length_bp"]
    return np.asarray([lengths.get(gene, np.nan) for gene in archs4_genes], dtype=np.float64)


def _extract_full_transcriptome_tpm(
    *,
    archs4_h5: str | Path,
    candidates: pd.DataFrame,
    landmark_genes: list[str],
    gene_lengths: str | Path,
    profiles: int,
) -> tuple[np.ndarray, pd.DataFrame, dict[str, Any]]:
    started = time.time()
    retained_rows: list[int] = []
    skipped_indices: list[int] = []
    sample_qc: list[tuple[float, float]] = []
    output = np.empty((int(profiles), len(landmark_genes)), dtype=np.float32)
    with h5py.File(archs4_h5, "r") as handle:
        expression = handle["data/expression"]
        archs4_genes = _decode(handle["meta/genes/ensembl_gene"][:])
        gene_map = {gene: index for index, gene in enumerate(archs4_genes)}
        missing = [gene for gene in landmark_genes if gene not in gene_map]
        if missing:
            raise ValueError(f"ARCHS4 lacks {len(missing)} paper-parity landmarks")
        lengths = _length_vector(archs4_genes, gene_lengths)
        valid_lengths = np.isfinite(lengths) & (lengths > 0)
        landmark_indices = np.asarray([gene_map[gene] for gene in landmark_genes])
        landmark_lengths_kb = lengths[landmark_indices] / 1000.0
        if not np.isfinite(landmark_lengths_kb).all():
            raise ValueError("One or more mouse landmarks lack a positive gene length")
        for candidate_row, row in candidates.iterrows():
            if len(retained_rows) >= int(profiles):
                break
            sample_index = int(row["archs4_sample_index"])
            try:
                counts = np.asarray(expression[:, sample_index], dtype=np.float64)
            except OSError:
                skipped_indices.append(sample_index)
                continue
            rpk_sum = float(np.sum(counts[valid_lengths] / (lengths[valid_lengths] / 1000.0)))
            if not np.isfinite(rpk_sum) or rpk_sum <= 0:
                skipped_indices.append(sample_index)
                continue
            output[len(retained_rows)] = (
                counts[landmark_indices] / landmark_lengths_kb / rpk_sum * 1_000_000.0
            ).astype(np.float32)
            retained_rows.append(int(candidate_row))
            sample_qc.append((float(counts.sum()), rpk_sum))
            if len(retained_rows) % 500 == 0:
                print(
                    f"[rna-diffusion:data] TPM {len(retained_rows)}/{profiles}",
                    flush=True,
                )
    if len(retained_rows) != int(profiles):
        raise RuntimeError(
            f"Only {len(retained_rows)} readable ARCHS4 profiles; need {profiles}"
        )
    retained = candidates.iloc[retained_rows].reset_index(drop=True).copy()
    retained["raw_count_sum"] = [value[0] for value in sample_qc]
    retained["full_transcriptome_rpk_sum"] = [value[1] for value in sample_qc]
    return output, retained, {
        "archs4_genes": len(archs4_genes),
        "genes_with_lengths": int(valid_lengths.sum()),
        "genes_without_lengths": int((~valid_lengths).sum()),
        "skipped_sample_indices": skipped_indices,
        "seconds": float(time.time() - started),
    }


def _split_indices(total: int, split: dict[str, int], seed: int) -> dict[str, np.ndarray]:
    indices = np.arange(int(total), dtype=np.int64)
    train_validation, test = train_test_split(
        indices,
        test_size=int(split["test"]),
        random_state=int(seed),
        shuffle=True,
    )
    train, validation = train_test_split(
        train_validation,
        test_size=int(split["validation"]),
        random_state=int(seed),
        shuffle=True,
    )
    result = {
        "train": np.sort(train),
        "validation": np.sort(validation),
        "test": np.sort(test),
    }
    for role, expected in split.items():
        if len(result[role]) != int(expected):
            raise RuntimeError(f"Split {role} has {len(result[role])}, expected {expected}")
    return result


def _series_group(value: object, fallback: str) -> str:
    if pd.isna(value) or not str(value).strip():
        return fallback
    return str(value).split(",", maxsplit=1)[0].strip()


def _series_ids(value: object) -> set[str]:
    if pd.isna(value) or not str(value).strip():
        return set()
    return {
        token.strip()
        for token in str(value).replace(";", ",").split(",")
        if token.strip()
    }


def _load_excluded_series_ids(data: dict[str, Any]) -> set[str]:
    excluded = {
        str(value).strip()
        for value in data.get("exclude_series_ids", [])
        if str(value).strip()
    }
    source = data.get("exclude_series_ids_file")
    if source:
        table = pd.read_csv(source, sep="\t", dtype=str)
        column = str(data.get("exclude_series_ids_column", "series_id"))
        if column not in table:
            raise ValueError(f"Series exclusion file lacks column {column!r}")
        excluded.update(
            value.strip()
            for value in table[column].dropna().astype(str)
            if value.strip()
        )
    return excluded


def _filter_excluded_series(
    metadata: pd.DataFrame,
    *,
    excluded: set[str],
    group_column: str = "series_id",
) -> tuple[pd.DataFrame, int]:
    if not excluded:
        return metadata.copy(), 0
    if group_column not in metadata:
        raise ValueError(
            f"ARCHS4 metadata lacks series exclusion column {group_column!r}"
        )
    mask = metadata[group_column].map(
        lambda value: bool(_series_ids(value) & excluded)
    )
    return metadata.loc[~mask].copy(), int(mask.sum())


def _group_split_indices(
    metadata: pd.DataFrame,
    *,
    fractions: dict[str, float],
    seed: int,
    group_column: str,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    if group_column not in metadata:
        raise ValueError(f"ARCHS4 metadata lacks split group column {group_column!r}")
    groups = np.asarray(
        [
            _series_group(value, f"profile:{index}")
            for index, value in enumerate(metadata[group_column])
        ],
        dtype=str,
    )
    classes = metadata["canonical_tissue"].astype(str).to_numpy()
    indices = np.arange(len(metadata), dtype=np.int64)
    expected_classes = set(classes)
    validation_relative = float(fractions["validation"]) / (
        float(fractions["train"]) + float(fractions["validation"])
    )
    for attempt in range(100):
        first = GroupShuffleSplit(
            n_splits=1,
            test_size=float(fractions["test"]),
            random_state=int(seed) + attempt,
        )
        train_validation, test = next(first.split(indices, groups=groups))
        second = GroupShuffleSplit(
            n_splits=1,
            test_size=validation_relative,
            random_state=int(seed) + 10_000 + attempt,
        )
        train_relative, validation_relative_indices = next(
            second.split(train_validation, groups=groups[train_validation])
        )
        result = {
            "train": np.sort(train_validation[train_relative]),
            "validation": np.sort(train_validation[validation_relative_indices]),
            "test": np.sort(test),
        }
        if all(set(classes[part]) == expected_classes for part in result.values()):
            role_groups = {role: set(groups[part]) for role, part in result.items()}
            if (
                role_groups["train"].isdisjoint(role_groups["validation"])
                and role_groups["train"].isdisjoint(role_groups["test"])
                and role_groups["validation"].isdisjoint(role_groups["test"])
            ):
                return result, groups
    raise RuntimeError(
        "Could not construct a series-held-out split containing every tissue class"
    )


def prepare(config_path: str | Path, *, force: bool = False) -> Path:
    config = load_config(config_path)
    data = config["data"]
    output = Path(data["prepared_h5"])
    manifest_path = output.with_suffix(".manifest.json")
    if output.exists() and manifest_path.exists() and not force:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        expected_split = (
            {key: int(value) for key, value in data["split"].items()}
            if "split" in data
            else None
        )
        if (
            int(manifest.get("profiles", -1)) != int(data["profiles"])
            or int(manifest.get("genes", -1)) != 974
            or (expected_split is not None and manifest.get("split") != expected_split)
            or manifest.get("split_strategy", "random_profile")
            != data.get("split_strategy", "random_profile")
        ):
            raise ValueError(
                f"Existing prepared data do not match {config_path}; rerun with --force"
            )
        print(f"[rna-diffusion:data] using existing {output}", flush=True)
        return output

    panel_path, panel_manifest = build_mouse_landmark_panel(
        human_landmarks=data["human_landmarks"],
        source_map=data["source_landmark_map"],
        archs4_h5=data["archs4_h5"],
        output=data["mouse_landmark_panel"],
        dimensions=974,
    )
    panel = pd.read_csv(panel_path, sep="\t")
    landmark_genes = panel["mouse_ensembl_gene"].astype(str).tolist()
    metadata = pd.read_csv(data["cohort_metadata"], sep="\t", low_memory=False)
    excluded_series_ids = _load_excluded_series_ids(data)
    metadata, excluded_series_profile_count = _filter_excluded_series(
        metadata,
        excluded=excluded_series_ids,
        group_column=str(data.get("exclude_series_group_column", "series_id")),
    )
    excluded_sample_indices = {
        int(value) for value in data.get("exclude_archs4_sample_indices", [])
    }
    if excluded_sample_indices:
        metadata = metadata.loc[
            ~metadata["archs4_sample_index"].astype(int).isin(excluded_sample_indices)
        ].copy()
    if data.get("profiles_per_tissue"):
        candidates, available_by_tissue = _targeted_candidates(
            metadata,
            quotas={
                str(tissue): int(count)
                for tissue, count in data["profiles_per_tissue"].items()
            },
            seed=int(data["selection_seed"]),
        )
    else:
        candidates, available_by_tissue = _balanced_candidates(
            metadata,
            minimum=int(data["minimum_tissue_profiles"]),
            seed=int(data["selection_seed"]),
        )
    tpm, retained, extraction = _extract_full_transcriptome_tpm(
        archs4_h5=data["archs4_h5"],
        candidates=candidates,
        landmark_genes=landmark_genes,
        gene_lengths=data["gene_lengths"],
        profiles=int(data["profiles"]),
    )
    classes = sorted(retained["canonical_tissue"].astype(str).unique())
    class_map = {name: index for index, name in enumerate(classes)}
    retained["class_index"] = (
        retained["canonical_tissue"].astype(str).map(class_map).astype(int)
    )
    if data.get("split_strategy") == "series_holdout":
        partitions, split_groups = _group_split_indices(
            retained,
            fractions={
                role: float(value)
                for role, value in data["split_fractions"].items()
            },
            seed=int(data["split_seed"]),
            group_column=str(data.get("split_group_column", "series_id")),
        )
        retained["split_group"] = split_groups
    else:
        partitions = _split_indices(
            len(retained), data["split"], seed=int(data["split_seed"])
        )
        retained["split_group"] = retained["geo_accession"].astype(str)
    retained["role"] = ""
    for role, indices in partitions.items():
        retained.loc[indices, "role"] = role
    coverage = pd.crosstab(retained["role"], retained["canonical_tissue"]).reindex(
        index=("train", "validation", "test"), columns=classes, fill_value=0
    )
    if (coverage == 0).any().any():
        raise RuntimeError("A paper-parity split omitted a tissue class")

    scale = np.max(np.abs(tpm[partitions["train"]]), axis=0).astype(np.float32)
    scale[scale == 0] = 1.0
    scaled = (tpm / scale.reshape(1, -1)).astype(np.float32)
    if not np.isfinite(scaled).all():
        raise FloatingPointError("Prepared paper-parity expression contains non-finite values")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp.h5")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["format"] = "nasa_mouse_lacan_paper_parity_v1"
        handle.attrs["normalization"] = str(data["normalization"])
        handle.create_dataset(
            "genes",
            data=np.asarray(landmark_genes, dtype=h5py.string_dtype("utf-8")),
        )
        handle.create_dataset(
            "classes", data=np.asarray(classes, dtype=h5py.string_dtype("utf-8"))
        )
        handle.create_dataset("maxabs_scale", data=scale)
        for role, indices in partitions.items():
            group = handle.create_group(role)
            group.create_dataset(
                "expression",
                data=scaled[indices],
                chunks=(min(256, len(indices)), len(landmark_genes)),
                compression="lzf",
            )
            group.create_dataset(
                "class_index", data=retained.loc[indices, "class_index"].to_numpy()
            )
            group.create_dataset("source_row", data=indices)
    temporary.replace(output)
    samples_path = output.with_suffix(".samples.tsv.gz")
    retained.to_csv(samples_path, sep="\t", index=False, compression="gzip")

    counts = retained.groupby(["role", "canonical_tissue"]).size().unstack(fill_value=0)
    manifest = {
        "prepared_h5": str(output),
        "prepared_h5_sha256": _sha256(output),
        "samples": str(samples_path),
        "profiles": len(retained),
        "genes": len(landmark_genes),
        "classes": classes,
        "class_count": len(classes),
        "split": {role: int(len(indices)) for role, indices in partitions.items()},
        "split_strategy": str(data.get("split_strategy", "random_profile")),
        "split_group_column": str(data.get("split_group_column", "geo_accession")),
        "split_group_overlap": {
            "train_validation": int(
                len(
                    set(retained.loc[partitions["train"], "split_group"])
                    & set(retained.loc[partitions["validation"], "split_group"])
                )
            ),
            "train_test": int(
                len(
                    set(retained.loc[partitions["train"], "split_group"])
                    & set(retained.loc[partitions["test"], "split_group"])
                )
            ),
            "validation_test": int(
                len(
                    set(retained.loc[partitions["validation"], "split_group"])
                    & set(retained.loc[partitions["test"], "split_group"])
                )
            ),
        },
        "profiles_by_role_and_tissue": {
            str(role): {str(tissue): int(value) for tissue, value in row.items()}
            for role, row in counts.iterrows()
        },
        "eligible_profiles_by_tissue": available_by_tissue,
        "selection": (
            "Configured tissue quotas from the full-catalog healthy-preferred "
            "selection, followed by a series-held-out split."
            if data.get("profiles_per_tissue")
            else "Seeded within-tissue shuffle followed by round-robin unique-profile "
            "selection across tissues with at least 100 eligible profiles."
        ),
        "normalization": (
            "TPM denominator over all 52,848 ARCHS4 genes with GENCODE M39 union-exon "
            f"lengths, landmark subset second, MaxAbs fitted on the "
            f"{len(partitions['train']):,} training profiles."
        ),
        "landmark_panel": panel_manifest,
        "extraction": extraction,
        "source": {
            "archs4_h5": str(data["archs4_h5"]),
            "cohort_metadata": str(data["cohort_metadata"]),
            "gene_lengths": str(data["gene_lengths"]),
            "excluded_series_ids_file": data.get("exclude_series_ids_file"),
            "excluded_series_ids": sorted(excluded_series_ids),
            "excluded_series_profile_count": excluded_series_profile_count,
            "excluded_unreadable_archs4_sample_indices": sorted(
                excluded_sample_indices
            ),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return output


def load_prepared(path: str | Path) -> dict[str, Any]:
    result: dict[str, Any] = {}
    with h5py.File(path, "r") as handle:
        if handle.attrs.get("format") != "nasa_mouse_lacan_paper_parity_v1":
            raise ValueError(f"Unsupported prepared-data format: {path}")
        result["genes"] = _decode(handle["genes"][:])
        result["classes"] = _decode(handle["classes"][:])
        result["maxabs_scale"] = np.asarray(handle["maxabs_scale"][:], dtype=np.float32)
        for role in ("train", "validation", "test"):
            result[role] = {
                "expression": np.asarray(handle[f"{role}/expression"][:], dtype=np.float32),
                "class_index": np.asarray(handle[f"{role}/class_index"][:], dtype=np.int64),
                "source_row": np.asarray(handle[f"{role}/source_row"][:], dtype=np.int64),
            }
            if len(result[role]["expression"]) != len(result[role]["class_index"]):
                raise RuntimeError(f"Prepared {role} expression/label lengths differ")
    if len(result["genes"]) != 974:
        raise ValueError(f"Prepared data have {len(result['genes'])} genes, expected 974")
    return result
