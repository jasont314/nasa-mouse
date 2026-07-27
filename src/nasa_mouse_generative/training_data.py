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


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_identity(path: str | Path, *, hash_limit_gb: float = 1.0) -> dict[str, object]:
    source = Path(path)
    if not source.exists():
        return {"path": str(source), "exists": False}
    stat = source.stat()
    result: dict[str, object] = {
        "path": str(source.resolve()),
        "exists": True,
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }
    if stat.st_size <= hash_limit_gb * 1024**3:
        result["sha256"] = _sha256_file(source)
    else:
        result["sha256"] = None
        result["hash_policy"] = "omitted_for_source_larger_than_hash_limit"
    return result


def _values_sha256(values: Iterable[object]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8", "replace"))
        digest.update(b"\n")
    return digest.hexdigest()


def _prepared_identity(
    config: BenchmarkConfig,
    genes: list[str],
    partitions: dict[str, "DataPartition"],
    reference: "DataPartition | None",
) -> dict[str, object]:
    selection = _archs4_selection_path(config)
    osdr_used = config.training.regime != "archs4_only"
    return {
        "genes_sha256": _values_sha256(genes),
        "partition_profile_ids_sha256": {
            name: _values_sha256(partition.obs.get("profile_id", []))
            for name, partition in partitions.items()
        },
        "partition_accessions_sha256": {
            name: _values_sha256(partition.obs.get("accession", []))
            for name, partition in partitions.items()
        },
        "reference_profile_ids_sha256": (
            _values_sha256(reference.obs.get("profile_id", []))
            if reference is not None
            else ""
        ),
        "sources": {
            "osdr_api_expression": (
                _file_identity(config.data.osdr_h5ad)
                if osdr_used
                else {"path": config.data.osdr_h5ad, "used": False}
            ),
            "osdr_api_metadata": (
                _file_identity(config.data.osdr_metadata)
                if osdr_used
                else {"path": config.data.osdr_metadata, "used": False}
            ),
            "archs4_h5": _file_identity(config.data.archs4_h5),
            "archs4_selection": _file_identity(selection),
        },
    }


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


def _read_archs4_column(
    expression: h5py.Dataset,
    sorted_gene_indices: np.ndarray,
    sample_index: int,
) -> np.ndarray:
    """Read one column using the source H5 chunk layout efficiently."""

    indices = np.asarray(sorted_gene_indices, dtype=np.int64)
    if len(indices) == 0:
        return np.empty(0, dtype=expression.dtype)
    contiguous = bool(indices[-1] - indices[0] + 1 == len(indices))
    if contiguous:
        return np.asarray(
            expression[int(indices[0]) : int(indices[-1]) + 1, int(sample_index)]
        )

    try:
        return np.asarray(expression[indices, int(sample_index)])
    except OSError as indexed_error:
        try:
            full_column = np.asarray(expression[:, int(sample_index)])
        except OSError as full_error:
            raise OSError(
                "ARCHS4 expression decompression failed for sample column "
                f"{sample_index} using both indexed and contiguous reads"
            ) from full_error
        if full_column.shape[0] != expression.shape[0]:
            raise OSError(
                f"ARCHS4 fallback read for sample {sample_index} was truncated"
            ) from indexed_error
        return full_column[indices]


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
                _read_archs4_column(
                    expression, sorted_indices, int(sample_index)
                ),
                dtype=np.float64,
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
            skipped = (
                np.asarray(handle["skipped_sample_indices"][:], dtype=np.int64).tolist()
                if "skipped_sample_indices" in handle
                else []
            )
        if matrix.shape == (len(metadata) - len(skipped), len(genes)):
            return matrix, {
                "cache": str(cache_path),
                "cache_hit": True,
                "requested_profiles": int(len(metadata)),
                "retained_profiles": int(len(matrix)),
                "skipped_corrupt_sample_indices": skipped,
            }

    matrix = np.empty((len(metadata), len(genes)), dtype=np.float32)
    retained = np.ones(len(metadata), dtype=bool)
    skipped: list[int] = []
    order = np.argsort(gene_indices)
    sorted_gene_indices = gene_indices[order]
    with h5py.File(source, "r") as handle:
        expression = handle["data/expression"]
        for row in np.argsort(sample_indices, kind="stable"):
            sample_index = sample_indices[row]
            try:
                values = np.asarray(
                    _read_archs4_column(
                        expression, sorted_gene_indices, int(sample_index)
                    ),
                    dtype=np.float32,
                )
            except OSError as error:
                retained[row] = False
                skipped.append(int(sample_index))
                if len(skipped) > config.data.archs4_max_corrupt_profiles:
                    raise OSError(
                        "ARCHS4 unreadable-profile count exceeded "
                        f"data.archs4_max_corrupt_profiles="
                        f"{config.data.archs4_max_corrupt_profiles}"
                    ) from error
                continue
            matrix[row, order] = values
    matrix = matrix[retained]
    if not len(matrix):
        raise OSError("Every selected ARCHS4 profile was unreadable")
    if skipped:
        print(
            "[archs4] excluded unreadable sample columns: "
            + ",".join(map(str, skipped)),
            flush=True,
        )
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
            handle.create_dataset(
                "skipped_sample_indices",
                data=np.asarray(skipped, dtype=np.int64),
            )
            handle.attrs["genes_sha256"] = hashlib.sha256(
                "\n".join(genes).encode()
            ).hexdigest()
        temporary.replace(cache_path)
    return matrix, {
        "cache": str(cache_path),
        "cache_hit": False,
        "requested_profiles": int(len(metadata)),
        "retained_profiles": int(len(matrix)),
        "skipped_corrupt_sample_indices": skipped,
    }


def _retain_readable_archs4_metadata(
    metadata: pd.DataFrame,
    matrix: np.ndarray,
    cache_metadata: dict[str, object],
) -> pd.DataFrame:
    skipped = set(
        map(int, cache_metadata.get("skipped_corrupt_sample_indices", []))
    )
    if skipped:
        retained = metadata.loc[
            ~metadata["archs4_sample_index"].astype(int).isin(skipped)
        ].copy()
    else:
        retained = metadata.copy()
    retained = retained.reset_index(drop=True)
    if len(retained) != len(matrix):
        raise RuntimeError(
            "ARCHS4 cache metadata is not aligned with the extracted matrix"
        )
    return retained


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


def _select_archs4_only_features(
    config: BenchmarkConfig, training_metadata: pd.DataFrame
) -> tuple[list[str], dict[str, object]]:
    archs4_genes, gene_map = _archs4_gene_map(config.data.archs4_h5)
    space = config.features.space
    if space == "reactome_shared":
        allowed = _read_gmt_genes(config.features.reactome_gmt)
        candidates = [gene for gene in archs4_genes if gene in allowed]
    elif space == "l1000_landmarks":
        candidates = [
            gene
            for gene in _read_l1000_genes(config.features.l1000_map)
            if gene in gene_map
        ]
    else:
        candidates = archs4_genes
    if not candidates:
        raise ValueError(f"Feature space {space!r} selected no ARCHS4 genes")

    target = 0
    if space == "hvg":
        target = config.features.hvg_genes
    elif config.features.max_genes > 0:
        target = config.features.max_genes
    target = min(target, len(candidates)) if target else 0
    ranking_required = bool(target and target < len(candidates))
    ranking_metadata = training_metadata
    selection_limit = int(config.features.selection_sample_limit)
    if (
        ranking_required
        and selection_limit > 0
        and len(ranking_metadata) > selection_limit
    ):
        ranking_metadata = ranking_metadata.copy()
        ranking_metadata["_stable_key"] = ranking_metadata["geo_accession"].map(
            lambda value: _stable_hash(config.training.seed, "feature_rank", value)
        )
        ranking_metadata = ranking_metadata.sort_values(
            ["canonical_tissue", "_stable_key"], kind="stable"
        )
        ranking_metadata["_within_tissue"] = ranking_metadata.groupby(
            "canonical_tissue", sort=True
        ).cumcount()
        ranking_metadata = (
            ranking_metadata.sort_values(
                ["_within_tissue", "_stable_key"], kind="stable"
            )
            .head(selection_limit)
            .drop(columns=["_stable_key", "_within_tissue"])
        )
    selection_cache = ""
    if ranking_required and config.execution.cache_archs4:
        source = Path(config.data.archs4_h5)
        source_stat = source.stat()
        cache_digest = hashlib.sha256()
        cache_digest.update(
            f"{source.resolve()}:{source_stat.st_size}:{source_stat.st_mtime_ns}".encode()
        )
        cache_digest.update(
            np.asarray(
                ranking_metadata["archs4_sample_index"], dtype="<i8"
            ).tobytes()
        )
        cache_digest.update("\n".join(candidates).encode())
        cache_digest.update(
            f"{space}:{target}:{selection_limit}:{config.training.seed}".encode()
        )
        cache_dir = Path(config.output_root) / "cache" / "feature_selection"
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / f"archs4_{cache_digest.hexdigest()[:20]}.json"
        selection_cache = str(cache_path)
        if cache_path.exists():
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
            selected = list(map(str, payload.get("genes", [])))
            if len(selected) == target and set(selected).issubset(candidates):
                return selected, {
                    "space": space,
                    "selection_source": "archs4_training_series",
                    "selected_genes": len(selected),
                    "selection_fit_roles": ["archs4_train"],
                    "selection_profiles_available": int(len(training_metadata)),
                    "selection_profiles_used": int(len(ranking_metadata)),
                    "selection_sample_limit": selection_limit,
                    "selection_cache": selection_cache,
                    "selection_cache_hit": True,
                }
    if ranking_required:
        candidate_indices = np.asarray(
            [gene_map[gene] for gene in candidates], dtype=np.int64
        )
        variances = _stream_archs4_variance(
            config.data.archs4_h5,
            ranking_metadata["archs4_sample_index"].to_numpy(dtype=np.int64),
            candidate_indices,
        )
        order = np.argsort(-variances, kind="stable")[:target]
        candidates = [candidates[int(index)] for index in order]
        if selection_cache:
            Path(selection_cache).write_text(
                json.dumps(
                    {
                        "genes": candidates,
                        "profiles_sha256": _values_sha256(
                            ranking_metadata["geo_accession"]
                        ),
                        "fit_role": "archs4_train",
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
    return candidates, {
        "space": space,
        "selection_source": "archs4_training_series",
        "selected_genes": len(candidates),
        "selection_fit_roles": ["archs4_train"],
        "selection_profiles_available": int(len(training_metadata)),
        "selection_profiles_used": (
            int(len(ranking_metadata)) if ranking_required else 0
        ),
        "selection_sample_limit": selection_limit,
        "selection_cache": selection_cache,
        "selection_cache_hit": False,
    }


def _split_archs4_selection(
    metadata: pd.DataFrame, config: BenchmarkConfig
) -> tuple[pd.DataFrame, dict[str, object]]:
    split_input = metadata.rename(
        columns={
            "series_id": "id.accession",
            "canonical_tissue": "tissue_canonical",
            "geo_accession": "profile_id",
        }
    ).copy()
    split_input["condition_inferred"] = "archs4_reference"
    plan = build_pooled_plan(
        split_input,
        seed=config.training.seed,
        validation_fraction=config.validation.pooled_validation_fraction,
        test_fraction=config.validation.pooled_test_fraction,
    )
    role_names = {
        "training": "train",
        "validation": "validation",
        "locked_test": "test",
    }
    role_by_series = {
        str(row["id.accession"]): role_names[str(row["role"])]
        for _, row in plan.iterrows()
    }
    result = metadata.copy()
    result["role"] = result["series_id"].astype(str).map(role_by_series)
    if result["role"].isna().any():
        raise RuntimeError("ARCHS4 split omitted one or more GEO series")
    if not result["role"].eq("train").any():
        raise ValueError("ARCHS4 split contains no training profiles")
    return result, {
        "kind": "archs4_geo_series_grouped",
        "split_unit": "GEO series_id",
        "expression_values_used_for_split": False,
        "series_by_role": {
            role: int(
                result.loc[result["role"].eq(role), "series_id"].astype(str).nunique()
            )
            for role in PARTITION_NAMES
        },
        "profiles_by_role": {
            role: int(result["role"].eq(role).sum()) for role in PARTITION_NAMES
        },
    }


def _prepare_archs4_only(
    config: BenchmarkConfig, *, tissue: str | None
) -> PreparedTrainingData:
    requested_tissues = tuple(config.data.osdr_tissues)
    if tissue:
        requested_tissues = (tissue,)
    metadata = load_archs4_selection(config, requested_tissues)
    metadata, split_metadata = _split_archs4_selection(metadata, config)
    training_metadata = metadata.loc[metadata["role"].eq("train")].copy()
    genes, feature_metadata = _select_archs4_only_features(
        config, training_metadata
    )
    raw_matrix, cache_metadata = extract_archs4_matrix(config, metadata, genes)
    metadata = _retain_readable_archs4_metadata(
        metadata, raw_matrix, cache_metadata
    )
    if cache_metadata.get("skipped_corrupt_sample_indices"):
        split_metadata["profiles_by_role_before_readability_filter"] = dict(
            split_metadata["profiles_by_role"]
        )
        split_metadata["profiles_by_role"] = {
            role: int(metadata["role"].eq(role).sum()) for role in PARTITION_NAMES
        }
    lengths = _gene_lengths(config, genes)
    obs = archs4_conditioning_frame(metadata).reset_index(drop=True)
    obs["role"] = metadata["role"].to_numpy()
    covariates = effective_covariates(config)
    train_mask = obs["role"].eq("train").to_numpy()
    train_obs = obs.loc[train_mask].reset_index(drop=True)
    encoder = CategoryEncoder.fit([train_obs], covariates)
    processor = FittedPreprocessor(
        config.preprocessing,
        device_spec=config.execution.device,
        seed=config.training.seed,
    )
    train_transformed = processor.fit_transform(
        raw_matrix[train_mask],
        train_obs["study"],
        gene_lengths=lengths,
        metadata=train_obs,
    )

    partitions: dict[str, DataPartition] = {}
    for role in PARTITION_NAMES:
        mask = obs["role"].eq(role).to_numpy()
        role_obs = obs.loc[mask].reset_index(drop=True)
        if role == "train":
            transformed = train_transformed
        elif mask.any():
            transformed = processor.transform(
                raw_matrix[mask],
                role_obs["study"],
                gene_lengths=lengths,
                allow_transductive=config.validation.allow_transductive_preprocessing,
                metadata=role_obs,
            )
        else:
            transformed = np.empty((0, len(genes)), dtype=np.float32)
        partitions[role] = DataPartition(
            name=role,
            matrix=transformed,
            obs=role_obs,
            categories=encoder.transform(role_obs),
            weights=_sampling_weights(role_obs),
        )

    tissues = sorted(metadata["canonical_tissue"].astype(str).unique())
    return PreparedTrainingData(
        genes=genes,
        covariates=covariates,
        encoder=encoder,
        preprocessor=processor,
        partitions=partitions,
        reference=partitions["train"],
        metadata={
            "source": "full-catalog ARCHS4 selection",
            "raw_integrated_osdr_h5_used": False,
            "osdr_expression_used": False,
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
            "reference_samples": len(partitions["train"]),
            "preprocessing_fit_source": "archs4_training_series",
            "data_identity": _prepared_identity(
                config, genes, partitions, partitions["train"]
            ),
        },
    )


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
    if config.training.regime == "archs4_only":
        return _prepare_archs4_only(config, tissue=tissue)
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
        reference_metadata = _retain_readable_archs4_metadata(
            reference_metadata, reference_matrix, cache_metadata
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
            "data_identity": _prepared_identity(
                config, genes, partitions, reference
            ),
        },
    )


def save_prepared_osdr(
    data: PreparedTrainingData,
    directory: str | Path,
    *,
    include_matrix: bool = True,
) -> Path:
    output_dir = Path(directory)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "prepared_data.h5"
    for name, partition in data.partitions.items():
        partition.obs.to_csv(
            output_dir / f"{name}_obs.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
    if include_matrix:
        with h5py.File(path, "w") as handle:
            handle.create_dataset("genes", data=np.asarray(data.genes, dtype="S"))
            for name, partition in data.partitions.items():
                group = handle.create_group(name)
                group.create_dataset(
                    "matrix", data=partition.matrix, compression="lzf"
                )
                group.create_dataset("categories", data=partition.categories)
                group.create_dataset("weights", data=partition.weights)
        data.metadata.setdefault("data_identity", {})["prepared_data_h5"] = {
            "path": str(path.resolve()),
            "size_bytes": int(path.stat().st_size),
            "sha256": _sha256_file(path),
        }
    else:
        data.metadata.setdefault("data_identity", {})["prepared_data_h5"] = {
            "saved": False,
            "reason": "execution.save_prepared_data=false",
            "reconstruction": "deterministic_reprepare_from_resolved_config",
        }
    manifest_path = output_dir / "prepared_data_manifest.json"
    manifest_path.write_text(
        json.dumps(data.metadata, indent=2) + "\n", encoding="utf-8"
    )
    return path if include_matrix else manifest_path


def load_prepared_osdr(directory: str | Path) -> tuple[list[str], dict[str, DataPartition]]:
    input_dir = Path(directory)
    partitions: dict[str, DataPartition] = {}
    path = input_dir / "prepared_data.h5"
    if not path.exists():
        path = input_dir / "prepared_osdr.h5"
    with h5py.File(path, "r") as handle:
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
