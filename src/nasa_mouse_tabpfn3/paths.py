"""Path and cohort constants for the TabPFN3 OSDR project."""

from __future__ import annotations

from pathlib import Path


DEFAULT_METADATA = Path("data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv")
DEFAULT_OSDR_API_DIR = Path("data/osdr_api")
DEFAULT_OUTPUT_ROOT = Path("outputs/tabpfn3_osdr")

TARGET_TISSUES = (
    "liver",
    "skeletal_muscle",
    "skin",
    "kidney",
    "thymus",
    "spleen",
    "lung",
    "retina",
)

MUSCLE_GROUPS = (
    "soleus",
    "gastrocnemius",
    "quadriceps",
    "edl",
    "tibialis_anterior",
)

CONDITIONS = ("ground_control", "flight")

DEFAULT_DESIGN_COVARIATES = (
    "id.accession",
    "tabpfn_tissue",
    "tabpfn_muscle_group",
    "tabpfn_material_type",
    "tabpfn_sex",
    "tabpfn_strain",
    "tabpfn_genotype",
    "tabpfn_platform",
    "tabpfn_assay",
    "tabpfn_data_source",
    "tabpfn_project_identifier",
    "tabpfn_project_type",
)

