"""Data loading and preprocessing for conditional WGAN-GP runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from nasa_mouse_glare.io import require_import


DEFAULT_CATEGORICAL_COVARIATES = ("wgan_condition", "wgan_accession")


@dataclass
class PreparedData:
    """Aligned expression arrays and conditional covariates."""

    query_x: object
    query_obs: object
    genes: list[str]
    categorical_covariates: tuple[str, ...]
    vocabularies: dict[str, list[str]]
    query_categories: object
    reference_x: object | None = None
    reference_obs: object | None = None
    reference_categories: object | None = None
    mean: object | None = None
    std: object | None = None


def counts_matrix(adata):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_sparse = require_import("scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt")
    x = adata.X
    if scipy_sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def log1p_cpm(counts):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    counts = np.asarray(counts, dtype=np.float32)
    library = counts.sum(axis=1, keepdims=True)
    library = np.maximum(library, 1.0)
    return np.log1p(counts / library * 1_000_000.0).astype(np.float32)


def standardize(*arrays, clip: float = 10.0):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    reference = np.asarray(arrays[0], dtype=np.float32)
    mean = reference.mean(axis=0, dtype=np.float64).astype(np.float32)
    std = reference.std(axis=0, dtype=np.float64).astype(np.float32)
    std = np.where(std < 1e-6, 1.0, std).astype(np.float32)
    transformed = []
    for array in arrays:
        values = ((np.asarray(array, dtype=np.float32) - mean) / std).astype(np.float32)
        if clip and clip > 0:
            values = np.clip(values, -clip, clip).astype(np.float32)
        transformed.append(values)
    return tuple(transformed), mean, std


def gene_ids(adata) -> list[str]:
    if "gene_id" in adata.var:
        return [str(value) for value in adata.var["gene_id"].tolist()]
    return [str(value) for value in adata.var_names.tolist()]


def subset_genes(adata, genes: list[str]):
    current = gene_ids(adata)
    index = {gene: idx for idx, gene in enumerate(current)}
    positions = [index[gene] for gene in genes]
    return adata[:, positions].copy()


def align_anndata(query, reference=None, *, max_genes: int = 0):
    query_genes = gene_ids(query)
    if reference is None:
        genes = query_genes
    else:
        reference_genes = set(gene_ids(reference))
        genes = [gene for gene in query_genes if gene in reference_genes]
    if max_genes and max_genes > 0:
        genes = genes[:max_genes]
    query = subset_genes(query, genes)
    reference = subset_genes(reference, genes) if reference is not None else None
    return query, reference, genes


def _first_existing_column(frame, columns: tuple[str, ...]) -> str | None:
    for column in columns:
        if column in frame:
            return column
    return None


def harmonize_obs(obs, *, source: str):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    frame = obs.copy()
    if "profile_id" not in frame:
        frame["profile_id"] = frame.index.astype(str)

    if source == "osdr":
        condition_col = _first_existing_column(frame, ("condition_inferred",))
        accession_col = _first_existing_column(frame, ("id.accession",))
        frame["wgan_condition"] = (
            frame[condition_col].astype(str) if condition_col else "unknown_condition"
        )
        frame["wgan_accession"] = (
            frame[accession_col].astype(str) if accession_col else "unknown_accession"
        )
    else:
        source_col = _first_existing_column(frame, ("series_id", "archs4_condition", "geo_accession"))
        frame["wgan_condition"] = "archs4_reference"
        frame["wgan_accession"] = (
            frame[source_col].astype(str) if source_col else "unknown_archs4_source"
        )
        if "condition_inferred" not in frame:
            frame["condition_inferred"] = "archs4_reference"
        if "id.accession" not in frame:
            frame["id.accession"] = frame["wgan_accession"]

    frame["wgan_source"] = source
    for column in DEFAULT_CATEGORICAL_COVARIATES:
        frame[column] = frame[column].fillna("missing").astype(str)
    return pd.DataFrame(frame)


def build_vocabularies(frames, covariates: tuple[str, ...]) -> dict[str, list[str]]:
    vocabularies = {}
    for covariate in covariates:
        values = []
        for frame in frames:
            values.extend(frame[covariate].fillna("missing").astype(str).tolist())
        vocabularies[covariate] = sorted(set(values))
    return vocabularies


def encode_categories(frame, vocabularies: dict[str, list[str]], covariates: tuple[str, ...]):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    columns = []
    for covariate in covariates:
        mapping = {value: idx for idx, value in enumerate(vocabularies[covariate])}
        columns.append(
            frame[covariate].fillna("missing").astype(str).map(mapping).fillna(0).to_numpy()
        )
    return np.stack(columns, axis=1).astype("int64")


def load_prepared_data(
    *,
    query_h5ad: str | Path,
    reference_h5ad: str | Path | None = None,
    categorical_covariates: tuple[str, ...] = DEFAULT_CATEGORICAL_COVARIATES,
    clip: float = 10.0,
    max_genes: int = 0,
) -> PreparedData:
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    query = ad.read_h5ad(query_h5ad)
    reference = ad.read_h5ad(reference_h5ad) if reference_h5ad else None
    query, reference, genes = align_anndata(query, reference, max_genes=max_genes)

    query_obs = harmonize_obs(query.obs, source="osdr")
    reference_obs = harmonize_obs(reference.obs, source="archs4") if reference is not None else None

    query_log = log1p_cpm(counts_matrix(query))
    if reference is not None:
        reference_log = log1p_cpm(counts_matrix(reference))
        (reference_x, query_x), mean, std = standardize(reference_log, query_log, clip=clip)
        frames = [reference_obs, query_obs]
    else:
        (query_x,), mean, std = standardize(query_log, clip=clip)
        reference_x = None
        frames = [query_obs]

    vocabularies = build_vocabularies(frames, categorical_covariates)
    query_categories = encode_categories(query_obs, vocabularies, categorical_covariates)
    reference_categories = (
        encode_categories(reference_obs, vocabularies, categorical_covariates)
        if reference_obs is not None
        else None
    )
    return PreparedData(
        query_x=query_x,
        query_obs=query_obs,
        reference_x=reference_x,
        reference_obs=reference_obs,
        genes=genes,
        categorical_covariates=categorical_covariates,
        vocabularies=vocabularies,
        query_categories=query_categories,
        reference_categories=reference_categories,
        mean=mean,
        std=std,
    )
