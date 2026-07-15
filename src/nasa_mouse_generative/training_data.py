"""Leakage-safe OSDR/ARCHS4 preparation for executable model adapters."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Iterable

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .conditioning import (
    CategoryEncoder,
    archs4_conditioning_frame,
    osdr_conditioning_frame,
)
from .config import BenchmarkConfig
from .preprocessing import FittedPreprocessor
from .split_plan import build_pooled_plan


PARTITION_NAMES = ("train", "validation", "test")


@dataclass
class DataPartition:
    name: str
    matrix: np.ndarray
    obs: pd.DataFrame
    categories: np.ndarray
    weights: np.ndarray

    def __len__(self) -> int:
        return int(self.matrix.shape[0])


@dataclass
class PreparedTrainingData:
    genes: list[str]
    covariates: tuple[str, ...]
    encoder: CategoryEncoder
    preprocessor: FittedPreprocessor
    partitions: dict[str, DataPartition]
    reference: DataPartition | None = None
    metadata: dict[str, object] = field(default_factory=dict)

    @property
    def train(self) -> DataPartition:
        return self.partitions["train"]


def _dense(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values, dtype=np.float32)


def _stable_hash(*parts: object) -> int:
    value = "|".join(map(str, parts)).encode("utf-8", "replace")
    return int.from_bytes(hashlib.sha256(value).digest()[:8], "big", signed=False)


def effective_covariates(config: BenchmarkConfig) -> tuple[str, ...]:
    covariates = list(config.training.conditioning_covariates)
    if not config.training.condition_on_flight:
        covariates = [name for name in covariates if name != "condition"]
    return tuple(covariates)


def _load_split_roles(
    config: BenchmarkConfig, tissue: str | None
) -> tuple[dict[str, str], dict[str, object]]:
    split_dir = Path(config.data.split_dir)
    if config.training.tissue_mode == "pooled_conditioned":
        table = pd.read_csv(split_dir / "pooled_accession_split.tsv", sep="\t")
        role_map = {
            "training": "train",
            "validation": "validation",
            "locked_test": "test",
        }
        roles = {
            str(row["id.accession"]): role_map[str(row["role"])]
            for _, row in table.iterrows()
        }
        return roles, {"kind": "pooled_accession", "fold_id": ""}

    if not tissue:
        raise ValueError("per_tissue training requires exactly one tissue")
    locked_path = split_dir / "per_tissue_locked_accession_splits.tsv"
    locked = pd.read_csv(locked_path, sep="\t")
    locked = locked.loc[locked["tissue_canonical"].astype(str).eq(tissue)]
    if not locked.empty:
        role_map = {
            "training": "train",
            "training_unpaired": "train",
            "validation": "validation",
            "locked_test": "test",
        }
        roles = {
            str(row["id.accession"]): role_map[str(row["role"])]
            for _, row in locked.iterrows()
        }
        return roles, {"kind": "per_tissue_locked", "fold_id": ""}

    loo = pd.read_csv(
        split_dir / "per_tissue_loo_accession_folds.tsv", sep="\t"
    )
    loo = loo.loc[loo["tissue_canonical"].astype(str).eq(tissue)]
    if loo.empty:
        raise ValueError(
            f"{tissue!r} is pooled-only and has no defensible standalone split"
        )
    fold_id = config.validation.fold_id or sorted(loo["fold_id"].astype(str).unique())[0]
    fold = loo.loc[loo["fold_id"].astype(str).eq(fold_id)]
    if fold.empty:
        raise ValueError(f"Unknown LOO fold {fold_id!r} for tissue {tissue}")
    roles = {
        str(row["id.accession"]): (
            "validation" if str(row["role"]) == "held_out" else "train"
        )
        for _, row in fold.iterrows()
    }
    return roles, {"kind": "per_tissue_loo", "fold_id": fold_id}


def _stratified_limit(frame: pd.DataFrame, limit: int, seed: int) -> pd.DataFrame:
    if limit <= 0 or len(frame) <= limit:
        return frame.copy()
    result = frame.copy()
    strata = [column for column in ("tissue", "condition", "accession") if column in result]
    result["_stable_key"] = [
        _stable_hash(seed, index, *[result.iloc[index][column] for column in strata])
        for index in range(len(result))
    ]
    result["_within_stratum"] = result.groupby(strata, dropna=False)[
        "_stable_key"
    ].rank(method="first")
    result = result.sort_values(
        ["_within_stratum", "_stable_key"], kind="stable"
    ).head(limit)
    return result.drop(columns=["_stable_key", "_within_stratum"])


def _osdr_rows(
    config: BenchmarkConfig, tissue: str | None
) -> tuple[ad.AnnData, pd.DataFrame, dict[str, object]]:
    path = Path(config.data.osdr_h5ad)
    if not path.exists():
        raise FileNotFoundError(
            f"API-derived OSDR matrix not found: {path}. Run osdr-expression first."
        )
    adata = ad.read_h5ad(path)
    normalized = osdr_conditioning_frame(adata.obs)
    normalized["_row_index"] = np.arange(adata.n_obs, dtype=np.int64)
    tissues = tuple(config.data.osdr_tissues)
    if tissue:
        tissues = (tissue,)
    if tissues:
        normalized = normalized.loc[normalized["tissue"].isin(tissues)].copy()
    include = set(config.data.osdr_include_accessions)
    exclude = set(config.data.osdr_exclude_accessions)
    if include:
        normalized = normalized.loc[normalized["accession"].isin(include)].copy()
    if exclude:
        normalized = normalized.loc[~normalized["accession"].isin(exclude)].copy()
    scope = config.data.osdr_accession_scope
    if scope == "all_eligible":
        roles, split_metadata = _load_split_roles(config, tissue)
        normalized["role"] = normalized["accession"].map(roles)
    elif normalized["accession"].nunique() == 1:
        normalized["role"] = _single_accession_roles(
            normalized,
            seed=config.training.seed,
            validation_fraction=config.validation.pooled_validation_fraction,
            test_fraction=config.validation.pooled_test_fraction,
        )
        split_metadata = {
            "kind": "single_accession_stratified_sample_fallback",
            "fold_id": "",
            "limitation": "accession-held-out validation is impossible with one accession",
        }
    else:
        selected_metadata = normalized.rename(
            columns={
                "accession": "id.accession",
                "tissue": "tissue_canonical",
                "condition": "condition_inferred",
            }
        )
        plan = build_pooled_plan(
            selected_metadata,
            seed=config.training.seed,
            validation_fraction=config.validation.pooled_validation_fraction,
            test_fraction=config.validation.pooled_test_fraction,
        )
        role_map = {
            "training": "train",
            "validation": "validation",
            "locked_test": "test",
        }
        roles = {
            str(row["id.accession"]): role_map[str(row["role"])]
            for _, row in plan.iterrows()
        }
        normalized["role"] = normalized["accession"].map(roles)
        if not normalized["role"].eq("validation").any():
            normalized["role"] = _single_accession_roles(
                normalized,
                seed=config.training.seed,
                validation_fraction=config.validation.pooled_validation_fraction,
                test_fraction=config.validation.pooled_test_fraction,
            )
            split_metadata = {
                "kind": "selected_accessions_stratified_sample_fallback",
                "fold_id": "",
                "accessions": int(normalized["accession"].nunique()),
                "limitation": "selected accessions could not preserve all strata under accession holdout",
            }
        else:
            split_metadata = {
                "kind": "selected_accessions_grouped",
                "fold_id": "",
                "accessions": int(normalized["accession"].nunique()),
            }
    normalized = normalized.dropna(subset=["role"]).copy()
    if normalized.empty:
        raise ValueError("No OSDR profiles remain after tissue/accession/split filters")
    return adata, normalized, split_metadata


def _single_accession_roles(
    frame: pd.DataFrame,
    *,
    seed: int,
    validation_fraction: float,
    test_fraction: float,
) -> pd.Series:
    """Stratify one accession when study-level holdout is impossible."""

    roles = pd.Series("train", index=frame.index, dtype="object")
    strata = ["tissue", "condition"]
    for keys, group in frame.groupby(strata, dropna=False, sort=True):
        ordered = sorted(
            group.index,
            key=lambda index: _stable_hash(seed, "single", keys, frame.loc[index, "profile_id"]),
        )
        count = len(ordered)
        if count < 3:
            continue
        test_count = min(max(1, round(count * test_fraction)), count - 2)
        validation_count = min(
            max(1, round(count * validation_fraction)), count - test_count - 1
        )
        for index in ordered[:test_count]:
            roles.loc[index] = "test"
        for index in ordered[test_count : test_count + validation_count]:
            roles.loc[index] = "validation"
    return roles


def _read_gmt_genes(path: str | Path) -> set[str]:
    genes: set[str] = set()
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            genes.update(gene for gene in fields[2:] if gene.startswith("ENSMUSG"))
    return genes


def _read_l1000_genes(path: str | Path) -> list[str]:
    table = pd.read_csv(path, sep="\t")
    if "mouse_ensembl_gene" not in table:
        raise ValueError(f"L1000 map lacks mouse_ensembl_gene: {path}")
    genes: list[str] = []
    seen: set[str] = set()
    for value in table["mouse_ensembl_gene"].dropna().astype(str):
        for gene in value.split(";"):
            gene = gene.strip()
            if gene.startswith("ENSMUSG") and gene not in seen:
                seen.add(gene)
                genes.append(gene)
    return genes


def _log1p_cpm(values: np.ndarray) -> np.ndarray:
    matrix = np.asarray(values, dtype=np.float64)
    library = np.maximum(matrix.sum(axis=1, keepdims=True), 1e-8)
    return np.log1p(matrix / library * 1_000_000.0).astype(np.float32)


def _variance_rank(matrix: np.ndarray) -> np.ndarray:
    variances = np.var(_log1p_cpm(matrix), axis=0, dtype=np.float64)
    return np.argsort(-variances, kind="stable")


def _archs4_selection_path(config: BenchmarkConfig) -> Path:
    return (
        Path(config.data.archs4_catalog_dir)
        / f"archs4_{config.data.archs4_cohort}_balanced.tsv.gz"
    )


def load_archs4_selection(
    config: BenchmarkConfig, tissues: Iterable[str]
) -> pd.DataFrame:
    path = _archs4_selection_path(config)
    if not path.exists():
        raise FileNotFoundError(f"ARCHS4 balanced cohort not found: {path}")
    table = pd.read_csv(path, sep="\t", low_memory=False)
    tissue_set = set(map(str, tissues))
    if tissue_set:
        table = table.loc[table["canonical_tissue"].astype(str).isin(tissue_set)]
    if table.empty:
        raise ValueError("No ARCHS4 profiles match the requested OSDR tissues")
    limit = config.data.archs4_sample_limit
    if limit > 0 and len(table) > limit:
        table = table.copy()
        table["_stable_key"] = table["geo_accession"].map(
            lambda value: _stable_hash(config.training.seed, value)
        )
        table["_round"] = table.groupby("canonical_tissue", dropna=False).cumcount()
        table = (
            table.sort_values(["_round", "_stable_key"], kind="stable")
            .head(limit)
            .drop(columns=["_stable_key", "_round"])
        )
    return table.reset_index(drop=True)


def _archs4_gene_map(path: str | Path) -> tuple[list[str], dict[str, int]]:
    with h5py.File(path, "r") as handle:
        raw = handle["meta/genes/ensembl_gene"][:]
    genes = [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in raw
    ]
    return genes, {gene: index for index, gene in enumerate(genes)}


def _stream_archs4_variance(
    path: str | Path, sample_indices: np.ndarray, gene_indices: np.ndarray
) -> np.ndarray:
    """Compute log1p-CPM variance without materializing the full reference."""

    order = np.argsort(gene_indices)
    sorted_indices = gene_indices[order]
    count = 0
    mean = np.zeros(len(gene_indices), dtype=np.float64)
    m2 = np.zeros(len(gene_indices), dtype=np.float64)
    with h5py.File(path, "r") as handle:
        expression = handle["data/expression"]
        for sample_index in sample_indices:
            sorted_values = np.asarray(
                expression[sorted_indices, int(sample_index)], dtype=np.float64
            )
            values = np.empty_like(sorted_values)
            values[order] = sorted_values
            library = max(float(values.sum()), 1e-8)
            values = np.log1p(values / library * 1_000_000.0)
            count += 1
            delta = values - mean
            mean += delta / count
            m2 += delta * (values - mean)
    return m2 / max(count, 1)


def select_features(
    config: BenchmarkConfig,
    adata: ad.AnnData,
    osdr_rows: pd.DataFrame,
    reference_metadata: pd.DataFrame | None,
) -> tuple[list[str], dict[str, object]]:
    osdr_genes = [str(value) for value in adata.var_names]
    available = set(osdr_genes)
    space = config.features.space
    if space == "reactome_shared":
        allowed = _read_gmt_genes(config.features.reactome_gmt)
        candidates = [gene for gene in osdr_genes if gene in allowed]
    elif space == "l1000_landmarks":
        candidates = [
            gene for gene in _read_l1000_genes(config.features.l1000_map) if gene in available
        ]
    else:
        candidates = osdr_genes
    if not candidates:
        raise ValueError(f"Feature space {space!r} selected no shared genes")

    target = 0
    if space == "hvg":
        target = config.features.hvg_genes
    elif config.features.max_genes > 0:
        target = config.features.max_genes
    target = min(target, len(candidates)) if target else 0

    source = config.features.selection_source
    if source == "auto":
        source = (
            "archs4_reference"
            if config.training.regime == "archs4_only"
            else "osdr_train"
        )
    if target and target < len(candidates):
        if source == "osdr_train":
            train_rows = osdr_rows.loc[osdr_rows["role"].eq("train")]
            row_indices = train_rows["_row_index"].to_numpy(dtype=int)
            osdr_gene_map = {gene: index for index, gene in enumerate(osdr_genes)}
            positions = np.asarray([osdr_gene_map[gene] for gene in candidates])
            matrix = _dense(adata.X[row_indices][:, positions])
            order = _variance_rank(matrix)
        else:
            if reference_metadata is None:
                raise ValueError("ARCHS4 feature selection requires reference metadata")
            _, archs4_map = _archs4_gene_map(config.data.archs4_h5)
            gene_indices = np.asarray([archs4_map[gene] for gene in candidates], dtype=int)
            variances = _stream_archs4_variance(
                config.data.archs4_h5,
                reference_metadata["archs4_sample_index"].to_numpy(dtype=int),
                gene_indices,
            )
            order = np.argsort(-variances, kind="stable")
        candidates = [candidates[int(index)] for index in order[:target]]
    return candidates, {
        "space": space,
        "selection_source": source,
        "selected_genes": len(candidates),
        "selection_fit_roles": (
            ["train"] if source == "osdr_train" else ["archs4_reference"]
        ),
    }


def extract_archs4_matrix(
    config: BenchmarkConfig, metadata: pd.DataFrame, genes: list[str]
) -> tuple[np.ndarray, dict[str, object]]:
    source = Path(config.data.archs4_h5)
    _, gene_map = _archs4_gene_map(source)
    missing = [gene for gene in genes if gene not in gene_map]
    if missing:
        raise ValueError(f"{len(missing)} selected genes are absent from ARCHS4")
    sample_indices = metadata["archs4_sample_index"].to_numpy(dtype=np.int64)
    gene_indices = np.asarray([gene_map[gene] for gene in genes], dtype=np.int64)
    digest = hashlib.sha256()
    digest.update(str(source.resolve()).encode())
    source_stat = source.stat()
    digest.update(f"{source_stat.st_size}:{source_stat.st_mtime_ns}".encode())
    digest.update(np.asarray(sample_indices, dtype="<i8").tobytes())
    digest.update("\n".join(genes).encode())
    cache_key = digest.hexdigest()[:20]
    cache_dir = Path(config.output_root) / "cache" / "archs4"
    cache_path = cache_dir / f"{config.data.archs4_cohort}_{cache_key}.h5"
    if config.execution.cache_archs4 and cache_path.exists():
        with h5py.File(cache_path, "r") as handle:
            matrix = np.asarray(handle["expression"][:], dtype=np.float32)
        if matrix.shape == (len(metadata), len(genes)):
            return matrix, {"cache": str(cache_path), "cache_hit": True}

    matrix = np.empty((len(metadata), len(genes)), dtype=np.float32)
    order = np.argsort(gene_indices)
    sorted_gene_indices = gene_indices[order]
    with h5py.File(source, "r") as handle:
        expression = handle["data/expression"]
        for row, sample_index in enumerate(sample_indices):
            values = np.asarray(
                expression[sorted_gene_indices, int(sample_index)], dtype=np.float32
            )
            matrix[row, order] = values
    if config.execution.cache_archs4:
        cache_dir.mkdir(parents=True, exist_ok=True)
        temporary = cache_path.with_suffix(".tmp.h5")
        with h5py.File(temporary, "w") as handle:
            handle.create_dataset(
                "expression",
                data=matrix,
                chunks=(min(256, len(matrix)), min(1024, len(genes))),
                compression="lzf",
            )
            handle.attrs["genes_sha256"] = hashlib.sha256(
                "\n".join(genes).encode()
            ).hexdigest()
        temporary.replace(cache_path)
    return matrix, {"cache": str(cache_path), "cache_hit": False}


def _gene_lengths(config: BenchmarkConfig, genes: list[str]) -> np.ndarray | None:
    if config.preprocessing.library_normalization != "tpm":
        return None
    path = Path(config.preprocessing.gene_lengths)
    table = pd.read_csv(path, sep=None, engine="python")
    gene_column = next(
        (column for column in ("gene_id", "gene", "ensembl_gene") if column in table),
        None,
    )
    length_column = next(
        (column for column in ("length", "gene_length", "length_bp") if column in table),
        None,
    )
    if not gene_column or not length_column:
        raise ValueError("Gene-length table needs gene_id and length columns")
    lengths = table.set_index(gene_column)[length_column]
    missing = [gene for gene in genes if gene not in lengths.index]
    if missing:
        raise ValueError(f"Gene lengths missing for {len(missing)} selected genes")
    return lengths.loc[genes].to_numpy(dtype=np.float32)


def _sampling_weights(obs: pd.DataFrame) -> np.ndarray:
    if len(obs) == 0:
        return np.empty(0, dtype=np.float32)
    tissue_counts = obs.groupby("tissue")["profile_id"].transform("count")
    accession_counts = obs.groupby(["tissue", "accession"])["profile_id"].transform(
        "count"
    )
    accession_n = obs.groupby("tissue")["accession"].transform("nunique")
    weights = 1.0 / accession_n.clip(lower=1) / accession_counts.clip(lower=1)
    if "tissue" in obs and obs["tissue"].nunique() > 1:
        weights = weights / obs["tissue"].nunique()
    del tissue_counts
    values = weights.to_numpy(dtype=np.float64)
    return (values / values.sum()).astype(np.float32)


def prepare_training_data(
    config: BenchmarkConfig, *, tissue: str | None = None
) -> PreparedTrainingData:
    if config.training.tissue_mode == "per_tissue":
        configured = tuple(config.data.osdr_tissues)
        if tissue is None and len(configured) == 1:
            tissue = configured[0]
        if not tissue:
            raise ValueError("Pass --tissue for a per_tissue run")
    adata, rows, split_metadata = _osdr_rows(config, tissue)
    tissues = sorted(rows["tissue"].astype(str).unique())
    needs_reference = config.training.regime in {
        "archs4_only",
        "archs4_pretrain_osdr_finetune",
    }
    reference_metadata = load_archs4_selection(config, tissues) if needs_reference else None
    genes, feature_metadata = select_features(
        config, adata, rows, reference_metadata
    )
    gene_position = {str(gene): index for index, gene in enumerate(adata.var_names)}
    positions = np.asarray([gene_position[gene] for gene in genes], dtype=int)
    covariates = effective_covariates(config)

    raw_partitions: dict[str, tuple[np.ndarray, pd.DataFrame]] = {}
    for role in PARTITION_NAMES:
        subset = rows.loc[rows["role"].eq(role)].copy()
        subset = _stratified_limit(
            subset, config.data.osdr_sample_limit, config.training.seed
        )
        indices = subset["_row_index"].to_numpy(dtype=int)
        matrix = (
            _dense(adata.X[indices][:, positions])
            if len(indices)
            else np.empty((0, len(genes)), dtype=np.float32)
        )
        raw_partitions[role] = (
            matrix,
            subset.drop(columns=["_row_index"]).reset_index(drop=True),
        )
    if len(raw_partitions["train"][0]) == 0:
        raise ValueError("The selected split has no OSDR training profiles")

    reference_matrix = None
    reference_obs = None
    cache_metadata: dict[str, object] = {}
    if reference_metadata is not None:
        reference_matrix, cache_metadata = extract_archs4_matrix(
            config, reference_metadata, genes
        )
        reference_obs = archs4_conditioning_frame(reference_metadata).reset_index(
            drop=True
        )

    train_obs = raw_partitions["train"][1]
    encoder_frames = [train_obs]
    if reference_obs is not None:
        encoder_frames.insert(0, reference_obs)
    encoder = CategoryEncoder.fit(encoder_frames, covariates)
    processor = FittedPreprocessor(
        config.preprocessing,
        device_spec=config.execution.device,
        seed=config.training.seed,
    )
    lengths = _gene_lengths(config, genes)
    dedicated_harmonizer = config.preprocessing.harmonization in {
        "combat",
        "combat_seq",
        "mober",
    }
    joint_harmonizer_fit = dedicated_harmonizer and config.training.regime == (
        "archs4_pretrain_osdr_finetune"
    )

    if reference_matrix is not None:
        if joint_harmonizer_fit:
            combined_matrix = np.concatenate(
                [reference_matrix, raw_partitions["train"][0]]
            )
            combined_obs = pd.concat(
                [reference_obs, train_obs], ignore_index=True, sort=False
            ).fillna("__missing__")
            combined_transformed = processor.fit_transform(
                combined_matrix,
                combined_obs["study"],
                gene_lengths=lengths,
                metadata=combined_obs,
            )
            reference_transformed = combined_transformed[: len(reference_matrix)]
            train_transformed = combined_transformed[len(reference_matrix) :]
        else:
            reference_transformed = processor.fit_transform(
                reference_matrix,
                reference_obs["study"],
                gene_lengths=lengths,
                metadata=reference_obs,
            )
            processor.fit_additional_study_stats(
                raw_partitions["train"][0],
                train_obs["study"],
                gene_lengths=lengths,
            )
    else:
        reference_transformed = None
        train_transformed = processor.fit_transform(
            raw_partitions["train"][0],
            train_obs["study"],
            gene_lengths=lengths,
            metadata=train_obs,
        )

    partitions: dict[str, DataPartition] = {}
    for role, (matrix, obs) in raw_partitions.items():
        if role == "train" and (reference_matrix is None or joint_harmonizer_fit):
            transformed = train_transformed
        elif len(matrix) == 0:
            transformed = np.empty((0, len(genes)), dtype=np.float32)
        else:
            transformed = processor.transform(
                matrix,
                obs["study"],
                gene_lengths=lengths,
                allow_transductive=(
                    config.validation.allow_transductive_preprocessing
                    and (
                        role in {"validation", "test"}
                        or config.training.regime == "archs4_only"
                    )
                ),
                metadata=obs,
            )
        partitions[role] = DataPartition(
            name=role,
            matrix=transformed,
            obs=obs,
            categories=encoder.transform(obs),
            weights=_sampling_weights(obs),
        )

    reference = None
    if reference_obs is not None and reference_transformed is not None:
        tissue_studies = reference_metadata.groupby("canonical_tissue")[
            "series_id"
        ].transform("nunique")
        study_sizes = reference_metadata.groupby(
            ["canonical_tissue", "series_id"]
        )["geo_accession"].transform("size")
        tissue_count = max(reference_metadata["canonical_tissue"].nunique(), 1)
        weights = (
            1.0
            / tissue_count
            / tissue_studies.clip(lower=1)
            / study_sizes.clip(lower=1)
        ).to_numpy(dtype=np.float64)
        weights = (weights / weights.sum()).astype(np.float32)
        reference = DataPartition(
            name="reference",
            matrix=reference_transformed,
            obs=reference_obs,
            categories=encoder.transform(reference_obs),
            weights=weights,
        )

    conditions = set(partitions["train"].obs["condition"].astype(str))
    if config.training.condition_on_flight and not {
        "flight",
        "ground_control",
    }.issubset(conditions):
        raise ValueError("OSDR training data must retain both flight and ground_control")
    return PreparedTrainingData(
        genes=genes,
        covariates=covariates,
        encoder=encoder,
        preprocessor=processor,
        partitions=partitions,
        reference=reference,
        metadata={
            "source": "NASA OSDR API plus full-catalog ARCHS4 selection",
            "raw_integrated_osdr_h5_used": False,
            "tissues": tissues,
            "split": split_metadata,
            "features": feature_metadata,
            "archs4": cache_metadata,
            "harmonization": processor.audit(),
            "transductive_preprocessing_enabled": (
                config.validation.allow_transductive_preprocessing
            ),
            "partition_samples": {
                name: len(partition) for name, partition in partitions.items()
            },
            "reference_samples": len(reference) if reference is not None else 0,
            "preprocessing_fit_source": (
                "archs4_reference_plus_osdr_train"
                if reference is not None and joint_harmonizer_fit
                else "archs4_reference"
                if reference is not None
                else "osdr_train"
            ),
        },
    )


def save_prepared_osdr(data: PreparedTrainingData, directory: str | Path) -> Path:
    output_dir = Path(directory)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "prepared_osdr.h5"
    with h5py.File(path, "w") as handle:
        handle.create_dataset("genes", data=np.asarray(data.genes, dtype="S"))
        for name, partition in data.partitions.items():
            group = handle.create_group(name)
            group.create_dataset("matrix", data=partition.matrix, compression="lzf")
            group.create_dataset("categories", data=partition.categories)
            group.create_dataset("weights", data=partition.weights)
            partition.obs.to_csv(
                output_dir / f"{name}_obs.tsv.gz",
                sep="\t",
                index=False,
                compression="gzip",
            )
    (output_dir / "prepared_data_manifest.json").write_text(
        json.dumps(data.metadata, indent=2) + "\n", encoding="utf-8"
    )
    return path


def load_prepared_osdr(directory: str | Path) -> tuple[list[str], dict[str, DataPartition]]:
    input_dir = Path(directory)
    partitions: dict[str, DataPartition] = {}
    with h5py.File(input_dir / "prepared_osdr.h5", "r") as handle:
        genes = [value.decode("utf-8") for value in handle["genes"][:]]
        for name in PARTITION_NAMES:
            group = handle[name]
            obs = pd.read_csv(input_dir / f"{name}_obs.tsv.gz", sep="\t")
            partitions[name] = DataPartition(
                name=name,
                matrix=np.asarray(group["matrix"][:], dtype=np.float32),
                obs=obs,
                categories=np.asarray(group["categories"][:], dtype=np.int64),
                weights=np.asarray(group["weights"][:], dtype=np.float32),
            )
    return genes, partitions
