"""Data loading and landmark/reconstruction preparation for diffusion runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from nasa_mouse_glare.io import require_import
from nasa_mouse_wgan.data import (
    CONDITIONAL_GENERATION_COVARIATES,
    align_anndata,
    build_vocabularies,
    counts_matrix,
    encode_categories,
    gene_ids,
    harmonize_obs,
    load_h5ads,
    log1p_cpm,
)


DIFFUSION_COVARIATES = CONDITIONAL_GENERATION_COVARIATES


@dataclass
class PreparedDiffusionData:
    """Aligned full-gene and landmark matrices plus encoded covariates."""

    query_full: object
    query_landmark: object
    query_obs: object
    query_categories: object
    genes: list[str]
    landmark_genes: list[str]
    target_genes: list[str]
    landmark_indices: object
    target_indices: object
    categorical_covariates: tuple[str, ...]
    vocabularies: dict[str, list[str]]
    full_center: object
    full_scale: object
    reference_full: object | None = None
    reference_landmark: object | None = None
    reference_obs: object | None = None
    reference_categories: object | None = None
    landmark_strategy: str = ""
    landmark_source: str = ""


def joined_paths(value) -> list[Path]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [Path(item) for item in value if str(item)]
    return [Path(item) for item in str(value).split(",") if str(item)]


def maxabs_scale(reference, *arrays):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    reference = np.asarray(reference, dtype=np.float32)
    center = np.zeros(reference.shape[1], dtype=np.float32)
    scale = np.max(np.abs(reference), axis=0).astype(np.float32)
    scale = np.where(scale < 1e-6, 1.0, scale).astype(np.float32)
    scaled = []
    for array in arrays:
        value = (np.asarray(array, dtype=np.float32) - center.reshape(1, -1)) / scale.reshape(1, -1)
        scaled.append(np.clip(value, -1.5, 1.5).astype(np.float32))
    return tuple(scaled), center, scale


def load_l1000_mouse_genes(path: str | Path | None, available: set[str]) -> list[str]:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    if not path:
        return []
    path = Path(path)
    if not path.exists():
        return []
    table = pd.read_csv(path, sep="\t")
    column = "mouse_ensembl_gene"
    if column not in table:
        return []
    genes = []
    seen = set()
    for value in table[column].dropna().astype(str):
        for gene in value.split(";"):
            gene = gene.strip()
            if gene.startswith("ENSMUSG") and gene in available and gene not in seen:
                genes.append(gene)
                seen.add(gene)
    return genes


def select_landmarks(
    matrix,
    genes: list[str],
    *,
    n_landmarks: int,
    strategy: str,
    l1000_map: str | Path | None = None,
    min_l1000: int = 300,
):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    available = set(genes)
    source = strategy
    selected: list[str] = []
    if strategy in {"l1000", "l1000_or_hvg"}:
        selected = load_l1000_mouse_genes(l1000_map, available)
        source = "mouse_l1000_orthologs"
        if strategy == "l1000" and len(selected) < min_l1000:
            raise SystemExit(
                f"Only {len(selected)} mapped mouse L1000 genes available; "
                f"need at least {min_l1000}."
            )
    if not selected or (strategy == "l1000_or_hvg" and len(selected) < min_l1000):
        variances = np.var(np.asarray(matrix, dtype=np.float32), axis=0)
        order = np.argsort(-variances)
        selected = [genes[int(idx)] for idx in order[: int(n_landmarks)]]
        source = "variance_hvg_fallback"
    else:
        selected = selected[: int(n_landmarks)]
    index = {gene: idx for idx, gene in enumerate(genes)}
    landmark_indices = np.asarray([index[gene] for gene in selected], dtype=np.int64)
    target_indices = np.asarray([idx for idx, gene in enumerate(genes) if gene not in set(selected)], dtype=np.int64)
    target_genes = [genes[int(idx)] for idx in target_indices]
    return selected, target_genes, landmark_indices, target_indices, source


def prepare_diffusion_data(
    *,
    query_h5ad: str,
    query_source: str,
    reference_h5ad: str = "",
    reference_source: str = "archs4",
    categorical_covariates: tuple[str, ...] = DIFFUSION_COVARIATES,
    max_genes: int = 0,
    n_landmarks: int = 512,
    landmark_strategy: str = "l1000_or_hvg",
    l1000_map: str | Path | None = "data/diffusion/l1000_human_to_mouse_ensembl.tsv",
):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    query = load_h5ads(query_h5ad, source=query_source)
    reference = load_h5ads(reference_h5ad, source=reference_source) if reference_h5ad else None
    if query is None:
        raise SystemExit("No query h5ad input supplied.")
    query, reference, genes = align_anndata(query, reference, max_genes=max_genes)

    query_obs = harmonize_obs(query.obs, source=query_source)
    reference_obs = harmonize_obs(reference.obs, source=reference_source) if reference is not None else None
    frames = [query_obs] + ([reference_obs] if reference_obs is not None else [])
    vocabularies = build_vocabularies(frames, categorical_covariates)
    query_categories = encode_categories(query_obs, vocabularies, categorical_covariates)
    reference_categories = (
        encode_categories(reference_obs, vocabularies, categorical_covariates)
        if reference_obs is not None
        else None
    )

    query_log = log1p_cpm(counts_matrix(query))
    reference_log = log1p_cpm(counts_matrix(reference)) if reference is not None else None
    scale_reference = reference_log if reference_log is not None else query_log
    if reference_log is not None:
        (reference_full, query_full), center, scale = maxabs_scale(scale_reference, reference_log, query_log)
    else:
        (query_full,), center, scale = maxabs_scale(scale_reference, query_log)
        reference_full = None

    landmark_source_matrix = reference_full if reference_full is not None else query_full
    landmark_genes, target_genes, landmark_indices, target_indices, landmark_source = select_landmarks(
        landmark_source_matrix,
        genes,
        n_landmarks=n_landmarks,
        strategy=landmark_strategy,
        l1000_map=l1000_map,
    )
    query_landmark = np.asarray(query_full[:, landmark_indices], dtype=np.float32)
    reference_landmark = (
        np.asarray(reference_full[:, landmark_indices], dtype=np.float32)
        if reference_full is not None
        else None
    )
    return PreparedDiffusionData(
        query_full=np.asarray(query_full, dtype=np.float32),
        query_landmark=query_landmark,
        query_obs=query_obs,
        query_categories=query_categories,
        genes=genes,
        landmark_genes=landmark_genes,
        target_genes=target_genes,
        landmark_indices=landmark_indices,
        target_indices=target_indices,
        categorical_covariates=categorical_covariates,
        vocabularies=vocabularies,
        full_center=center,
        full_scale=scale,
        reference_full=np.asarray(reference_full, dtype=np.float32) if reference_full is not None else None,
        reference_landmark=reference_landmark,
        reference_obs=reference_obs,
        reference_categories=reference_categories,
        landmark_strategy=landmark_strategy,
        landmark_source=landmark_source,
    )


def _mode_or_default(series, default: str) -> str:
    values = series.dropna().astype(str)
    if values.empty:
        return default
    counts = values.value_counts()
    if counts.empty:
        return default
    return str(counts.index[0])


def reference_projection_obs(prepared: PreparedDiffusionData):
    """Replace query-specific covariates with trained reference defaults.

    ARCHS4 pretraining never observes OSDR-only condition/accession/source
    categories. Frozen reference projection therefore must not use those
    untrained embeddings. The projection keeps query expression and tissue, but
    uses per-tissue ARCHS4 defaults for all other categorical covariates.
    """

    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    if prepared.reference_obs is None:
        return prepared.query_obs.copy()
    ref = prepared.reference_obs
    query = prepared.query_obs.copy()
    global_defaults = {
        covariate: _mode_or_default(ref[covariate], "unknown")
        for covariate in prepared.categorical_covariates
    }
    by_tissue: dict[str, dict[str, str]] = {}
    if "wgan_tissue" in ref:
        for tissue, frame in ref.groupby("wgan_tissue", dropna=False):
            by_tissue[str(tissue)] = {
                covariate: _mode_or_default(frame[covariate], global_defaults[covariate])
                for covariate in prepared.categorical_covariates
            }
    rows = []
    observed_tissues = set(ref["wgan_tissue"].dropna().astype(str).tolist()) if "wgan_tissue" in ref else set()
    for _, row in query.iterrows():
        tissue = str(row.get("wgan_tissue", global_defaults.get("wgan_tissue", "unknown")))
        defaults = by_tissue.get(tissue, global_defaults)
        values = row.copy()
        for covariate in prepared.categorical_covariates:
            if covariate == "wgan_tissue" and tissue in observed_tissues:
                values[covariate] = tissue
            else:
                values[covariate] = defaults[covariate]
        rows.append(values)
    return pd.DataFrame(rows, index=query.index)


def reference_projection_categories(prepared: PreparedDiffusionData):
    projected = reference_projection_obs(prepared)
    return encode_categories(projected, prepared.vocabularies, prepared.categorical_covariates)


def write_observed_profiles(path: Path, prepared: PreparedDiffusionData) -> None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    frames = []
    query = prepared.query_obs[list(prepared.categorical_covariates)].drop_duplicates().copy()
    query.insert(0, "training_source", "query")
    frames.append(query)
    if prepared.reference_obs is not None:
        ref = prepared.reference_obs[list(prepared.categorical_covariates)].drop_duplicates().copy()
        ref.insert(0, "training_source", "reference")
        frames.append(ref)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(frames, ignore_index=True).to_csv(path, sep="\t", index=False)
