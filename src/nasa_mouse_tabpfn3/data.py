"""NASA OSDR API-derived data loading for TabPFN3 classifiers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from nasa_mouse_glare.io import require_import
from nasa_mouse_glare.prepare_expimap_osdr_tissue import (
    count_table_path,
    ensure_count_tables,
    load_counts_from_api_tables,
)

from .paths import CONDITIONS, DEFAULT_METADATA, DEFAULT_OSDR_API_DIR, MUSCLE_GROUPS


@dataclass
class OsdrExpressionDataset:
    """Expression matrix and metadata for one tissue or tissue subgroup."""

    dataset_id: str
    tissue: str
    split_group: str
    counts: object
    obs: object
    genes: list[str]

    @property
    def n_samples(self) -> int:
        return int(self.counts.shape[0])

    @property
    def n_genes(self) -> int:
        return int(self.counts.shape[1])


def safe_name(value: str) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_") or "unknown"


def infer_muscle_group(material: str, tissue: str = "") -> str:
    text = str(material).lower().replace("-", " ").replace("_", " ")
    tissue_text = str(tissue).lower().replace(" ", "_")
    if "soleus" in text:
        return "soleus"
    if "gastrocnemius" in text:
        return "gastrocnemius"
    if "quadriceps" in text:
        return "quadriceps"
    if "tibialis" in text:
        return "tibialis_anterior"
    if "extensor digitorum" in text or re.search(r"\bedl\b", text):
        return "edl"
    if tissue_text == "skeletal_muscle":
        return "skeletal_muscle_other"
    return "not_skeletal_muscle"


def load_metadata(path: str | Path = DEFAULT_METADATA):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} does not exist. Run "
            "`PYTHONPATH=src python -m nasa_mouse_glare.fetch_osdr_mouse_transcriptomics`."
        )
    metadata = pd.read_csv(path, sep="\t", keep_default_na=False)
    if "profile_id" not in metadata:
        metadata["profile_id"] = (
            metadata["id.accession"].astype(str)
            + "/"
            + metadata["id.sample name"].astype(str)
        )
    metadata["tabpfn_condition"] = metadata["condition_inferred"].astype(str)
    metadata["tabpfn_label"] = metadata["tabpfn_condition"].map(
        {"ground_control": 0, "flight": 1}
    )
    metadata["tabpfn_tissue"] = metadata["tissue_final"].map(safe_name)
    metadata["tabpfn_material_type"] = metadata[
        "study.characteristics.material type"
    ].astype(str)
    metadata["tabpfn_muscle_group"] = [
        infer_muscle_group(material, tissue)
        for material, tissue in zip(
            metadata["tabpfn_material_type"], metadata["tabpfn_tissue"]
        )
    ]
    return metadata


def selected_metadata(metadata, tissue: str, split_group: str = ""):
    selected = metadata.loc[
        metadata["tabpfn_tissue"].eq(safe_name(tissue))
        & metadata["tabpfn_condition"].isin(CONDITIONS)
    ].copy()
    if split_group:
        selected = selected.loc[
            selected["tabpfn_muscle_group"].eq(safe_name(split_group))
        ].copy()
    if "bulk_rna_seq_inferred" in selected:
        selected = selected.loc[
            selected["bulk_rna_seq_inferred"].map(
                lambda value: value is True or str(value).lower() == "true"
            )
        ].copy()
    selected = selected.drop_duplicates(
        subset=["id.accession", "id.assay name", "id.sample name"],
        keep="first",
    )
    return selected.sort_values(
        ["id.accession", "tabpfn_condition", "profile_id"],
        kind="stable",
    ).reset_index(drop=True)


def _validate_conditions(selected, min_per_condition: int) -> tuple[bool, str]:
    counts = selected["tabpfn_condition"].value_counts()
    missing = [
        condition
        for condition in CONDITIONS
        if int(counts.get(condition, 0)) < min_per_condition
    ]
    if missing:
        detail = ", ".join(
            f"{condition}={int(counts.get(condition, 0))}" for condition in CONDITIONS
        )
        return False, f"insufficient samples per condition ({detail})"
    return True, ""


def load_dataset(
    *,
    tissue: str,
    split_group: str = "",
    metadata_path: str | Path = DEFAULT_METADATA,
    api_dir: str | Path = DEFAULT_OSDR_API_DIR,
    min_per_condition: int = 3,
    timeout: int = 120,
    overwrite_counts: bool = False,
    gene_prefix: str = "ENSMUSG",
) -> OsdrExpressionDataset:
    metadata = load_metadata(metadata_path)
    selected = selected_metadata(metadata, tissue, split_group)
    valid, reason = _validate_conditions(selected, min_per_condition)
    if not valid:
        raise ValueError(f"{tissue}/{split_group or 'all'}: {reason}")

    api_dir = Path(api_dir)
    accessions = selected["id.accession"].astype(str).unique().tolist()
    missing = [
        accession
        for accession in accessions
        if overwrite_counts or not count_table_path(api_dir, accession).exists()
    ]
    if missing:
        ensure_count_tables(selected, api_dir, timeout, overwrite_counts)
    count_paths = {accession: count_table_path(api_dir, accession) for accession in accessions}
    count_matrix = load_counts_from_api_tables(selected, count_paths)
    count_matrix = count_matrix.loc[:, selected["profile_id"].astype(str).tolist()]
    if gene_prefix:
        keep_genes = count_matrix.index.astype(str).str.startswith(str(gene_prefix))
        count_matrix = count_matrix.loc[keep_genes].copy()
        if count_matrix.empty:
            raise ValueError(f"No genes matched prefix {gene_prefix!r}.")
    counts = count_matrix.transpose()

    obs = selected.set_index("profile_id", drop=False).loc[counts.index].copy()
    obs["tabpfn_dataset_id"] = dataset_id(tissue, split_group)
    genes = [str(gene) for gene in counts.columns.tolist()]
    return OsdrExpressionDataset(
        dataset_id=dataset_id(tissue, split_group),
        tissue=safe_name(tissue),
        split_group=safe_name(split_group) if split_group else "",
        counts=counts,
        obs=obs,
        genes=genes,
    )


def dataset_id(tissue: str, split_group: str = "") -> str:
    tissue_name = safe_name(tissue)
    if split_group:
        return f"{tissue_name}__{safe_name(split_group)}"
    return tissue_name


def log1p_cpm(counts):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    values = counts.to_numpy(dtype="float32", copy=True)
    library = values.sum(axis=1, keepdims=True)
    library[library <= 0] = 1.0
    return np.log1p(values / library * 1_000_000.0).astype("float32")


def planned_datasets(
    tissues: tuple[str, ...],
    *,
    include_muscle_splits: bool = True,
    muscle_groups: tuple[str, ...] = MUSCLE_GROUPS,
) -> list[tuple[str, str]]:
    plan = [(safe_name(tissue), "") for tissue in tissues]
    if include_muscle_splits and "skeletal_muscle" in {safe_name(t) for t in tissues}:
        plan.extend(("skeletal_muscle", group) for group in muscle_groups)
    return plan


def inventory_tables(
    metadata_path: str | Path = DEFAULT_METADATA,
    *,
    tissues: tuple[str, ...] | None = None,
):
    metadata = load_metadata(metadata_path)
    if tissues is not None:
        wanted = {safe_name(tissue) for tissue in tissues}
        metadata = metadata.loc[metadata["tabpfn_tissue"].isin(wanted)].copy()

    tissue_cols = ["tabpfn_tissue", "tabpfn_condition"]
    tissue_counts = metadata.groupby(tissue_cols, dropna=False).size().unstack(fill_value=0)
    for condition in CONDITIONS:
        if condition not in tissue_counts:
            tissue_counts[condition] = 0
    tissue_counts["total"] = tissue_counts[list(CONDITIONS)].sum(axis=1)
    tissue_counts = tissue_counts.reset_index().rename(columns={"tabpfn_tissue": "tissue"})

    accession_cols = ["tabpfn_tissue", "id.accession", "tabpfn_condition"]
    accession_counts = (
        metadata.groupby(accession_cols, dropna=False).size().unstack(fill_value=0)
    )
    for condition in CONDITIONS:
        if condition not in accession_counts:
            accession_counts[condition] = 0
    accession_counts["total"] = accession_counts[list(CONDITIONS)].sum(axis=1)
    accession_counts = accession_counts.reset_index().rename(
        columns={"tabpfn_tissue": "tissue"}
    )

    muscle = metadata.loc[metadata["tabpfn_tissue"].eq("skeletal_muscle")].copy()
    muscle_cols = ["tabpfn_muscle_group", "tabpfn_condition"]
    muscle_counts = muscle.groupby(muscle_cols, dropna=False).size().unstack(fill_value=0)
    for condition in CONDITIONS:
        if condition not in muscle_counts:
            muscle_counts[condition] = 0
    muscle_counts["total"] = muscle_counts[list(CONDITIONS)].sum(axis=1)
    muscle_counts = muscle_counts.reset_index().rename(
        columns={"tabpfn_muscle_group": "muscle_group"}
    )
    return tissue_counts, accession_counts, muscle_counts

