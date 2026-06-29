"""Data loading and preprocessing for conditional WGAN-GP runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from nasa_mouse_glare.io import require_import


DEFAULT_CATEGORICAL_COVARIATES = ("wgan_condition", "wgan_accession")
CONDITIONAL_GENERATION_COVARIATES = (
    "wgan_condition",
    "wgan_tissue",
    "wgan_material_type",
    "wgan_muscle_group",
    "wgan_accession",
    "wgan_sex",
    "wgan_assay",
    "wgan_platform",
    "wgan_data_source",
)


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


def split_h5ad_paths(value) -> list[Path]:
    """Parse one or more h5ad paths.

    Existing CLIs pass a single string. The conditional-generation workflow uses
    a comma-separated list to train a pan-tissue model without first materializing
    a combined AnnData file.
    """

    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [Path(item) for item in value if str(item)]
    return [Path(item) for item in str(value).split(",") if item]


def infer_tissue_from_path(path: str | Path) -> str:
    name = Path(path).name
    patterns = (
        r"osdr_(.+?)_flt_gc",
        r"archs4_mouse_(.+?)_reference",
    )
    for pattern in patterns:
        match = re.search(pattern, name)
        if match:
            return match.group(1)
    return "unknown_tissue"


def _common_genes(adatas) -> list[str]:
    genes = gene_ids(adatas[0])
    common = set(genes)
    for adata in adatas[1:]:
        common &= set(gene_ids(adata))
    return [gene for gene in genes if gene in common]


def load_h5ads(value, *, source: str):
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    paths = split_h5ad_paths(value)
    if not paths:
        return None
    adatas = []
    seen_profiles = set()
    for path in paths:
        adata = ad.read_h5ad(path)
        obs = adata.obs.copy()
        if "profile_id" not in obs:
            obs["profile_id"] = obs.index.astype(str)
        obs["input_h5ad"] = str(path)
        obs["wgan_input_tissue"] = infer_tissue_from_path(path)
        if "tissue_final" not in obs:
            obs["tissue_final"] = obs["wgan_input_tissue"]
        keep = ~obs["profile_id"].astype(str).isin(seen_profiles)
        if not bool(keep.all()):
            adata = adata[keep.to_numpy()].copy()
            obs = obs.loc[keep].copy()
        seen_profiles.update(obs["profile_id"].astype(str).tolist())
        adata.obs = obs
        adatas.append(adata)
    if len(adatas) == 1:
        return adatas[0]

    genes = _common_genes(adatas)
    aligned = []
    for adata in adatas:
        item = subset_genes(adata, genes)
        item.var_names = genes
        aligned.append(item)
    return ad.concat(aligned, join="inner", merge="first", index_unique=None)


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


def _series_from_columns(frame, columns: tuple[str, ...], default: str):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    column = _first_existing_column(frame, columns)
    if column is None:
        return pd.Series(default, index=frame.index, dtype="object")
    return frame[column].fillna(default).astype(str)


def infer_muscle_group(material, tissue) -> str:
    text = str(material).lower().replace("-", " ").replace("_", " ")
    tissue_text = str(tissue).lower()
    if "soleus" in text:
        return "soleus"
    if "gastrocnemius" in text:
        return "gastrocnemius"
    if "quadriceps" in text:
        return "quadriceps"
    if "tibialis" in text:
        return "tibialis_anterior"
    if "extensor digitorum" in text or "edl" in text:
        return "edl"
    if "skeletal_muscle" in tissue_text or "skeletal muscle" in tissue_text:
        return "skeletal_muscle_other"
    return "not_skeletal_muscle"


def simplify_platform(value: str) -> str:
    text = str(value).lower()
    if "hiseq 4000" in text:
        return "illumina_hiseq_4000"
    if "hiseq 2500" in text:
        return "illumina_hiseq_2500"
    if "hiseq" in text:
        return "illumina_hiseq"
    if "nextseq" in text:
        return "illumina_nextseq"
    if "novaseq" in text:
        return "illumina_novaseq"
    if "illumina" in text:
        return "illumina"
    if "rna-seq" in text or "rna seq" in text:
        return "rna_seq"
    if not text or text == "nan":
        return "unknown_platform"
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")[:80] or "unknown_platform"


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
        frame["wgan_tissue"] = _series_from_columns(
            frame, ("tissue_final", "wgan_input_tissue"), "unknown_tissue"
        )
        frame["wgan_material_type"] = _series_from_columns(
            frame, ("study.characteristics.material type",), "unknown_material_type"
        )
        frame["wgan_sex"] = _series_from_columns(
            frame, ("study.characteristics.sex",), "unknown_sex"
        )
        frame["wgan_assay"] = _series_from_columns(
            frame,
            ("investigation.study assays.study assay technology type", "file.datatype"),
            "unknown_assay",
        )
        raw_platform = _series_from_columns(frame, ("id.assay name",), "unknown_platform")
        frame["wgan_platform"] = raw_platform.map(simplify_platform)
    else:
        source_col = _first_existing_column(frame, ("series_id", "archs4_condition", "geo_accession"))
        frame["wgan_condition"] = "archs4_reference"
        frame["wgan_accession"] = (
            frame[source_col].astype(str) if source_col else "unknown_archs4_source"
        )
        frame["wgan_tissue"] = _series_from_columns(
            frame, ("tissue_final", "wgan_input_tissue"), "unknown_tissue"
        )
        frame["wgan_material_type"] = frame["wgan_tissue"]
        frame["wgan_sex"] = "unknown_sex"
        frame["wgan_assay"] = _series_from_columns(
            frame, ("library_strategy", "library_source"), "unknown_assay"
        )
        frame["wgan_platform"] = _series_from_columns(
            frame, ("library_source",), "unknown_platform"
        ).map(simplify_platform)
        if "condition_inferred" not in frame:
            frame["condition_inferred"] = "archs4_reference"
        if "id.accession" not in frame:
            frame["id.accession"] = frame["wgan_accession"]

    frame["wgan_source"] = source
    frame["wgan_data_source"] = source
    if "muscle_group" in frame:
        frame["wgan_muscle_group"] = frame["muscle_group"].fillna("").astype(str)
    else:
        frame["wgan_muscle_group"] = [
            infer_muscle_group(material, tissue)
            for material, tissue in zip(frame["wgan_material_type"], frame["wgan_tissue"])
        ]
    frame.loc[frame["wgan_muscle_group"].eq(""), "wgan_muscle_group"] = [
        infer_muscle_group(material, tissue)
        for material, tissue in zip(frame.loc[frame["wgan_muscle_group"].eq(""), "wgan_material_type"],
                                    frame.loc[frame["wgan_muscle_group"].eq(""), "wgan_tissue"])
    ]

    for column in set(DEFAULT_CATEGORICAL_COVARIATES) | set(CONDITIONAL_GENERATION_COVARIATES):
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
    query_source: str = "osdr",
    reference_source: str = "archs4",
    categorical_covariates: tuple[str, ...] = DEFAULT_CATEGORICAL_COVARIATES,
    clip: float = 10.0,
    max_genes: int = 0,
) -> PreparedData:
    query = load_h5ads(query_h5ad, source=query_source)
    reference = load_h5ads(reference_h5ad, source=reference_source) if reference_h5ad else None
    query, reference, genes = align_anndata(query, reference, max_genes=max_genes)

    query_obs = harmonize_obs(query.obs, source=query_source)
    reference_obs = harmonize_obs(reference.obs, source=reference_source) if reference is not None else None

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
