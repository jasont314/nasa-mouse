"""Prepare multi-tissue GLARE inputs from the NASA OSDR Biological Data API."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from types import SimpleNamespace

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import torch

from .aggregate_liver_mober import ensure_mober_importable, sanitize_anndata_strings
from .aggregate_liver_mober_glare import prepare_mober_target
from .aggregate_liver_finetune import apply_ercc_policy, clean_token
from .align import align_bundles
from .export import export_mtx
from .fetch_osdr_mouse_transcriptomics import (
    DEFAULT_OUTPUT_DIR as DEFAULT_API_DIR,
    discover_metadata,
    download_count_tables,
)
from .io import dense_matrix, load_matrix_bundle, require_import, write_matrix_bundle
from .paper_finetune import (
    finetune_location,
    format_elapsed,
    infer_pretrain_input_dim,
    log,
    write_outlier_audit,
)
from .tms import prepare_tms


DEFAULT_OUTPUT_DIR = "outputs/glare_multi_tissue_api"
DEFAULT_TMS_H5AD = "assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad"
DEFAULT_API_METADATA = (
    "data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv"
)
DEFAULT_COUNTS_DIR = "data/osdr_api/counts"
DEFAULT_REACTOME_GMT = (
    "src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt"
)
DEFAULT_MIN_TMS_CELLS = 100
DEFAULT_MIN_PROFILES_PER_CONDITION = 2
DEFAULT_PRETRAIN_EPOCHS = 30
DEFAULT_FINETUNE_EPOCHS = 30
DEFAULT_BATCH_SIZE = 16
DEFAULT_SEED = 1996


@dataclass(frozen=True)
class TissueSpec:
    slug: str
    label: str
    tissue_final: str
    tms_tissue: str | None
    material_terms: tuple[str, ...] = ()
    note: str = ""


BASE_SPECS = [
    TissueSpec("liver", "liver", "liver", "liver"),
    TissueSpec("skeletal_muscle", "skeletal muscle", "skeletal_muscle", "limb muscle"),
    TissueSpec("skin", "skin", "skin", "skin of body"),
    TissueSpec("kidney", "kidney", "kidney", "kidney"),
    TissueSpec("thymus", "thymus", "thymus", "thymus"),
    TissueSpec("spleen", "spleen", "spleen", "spleen"),
    TissueSpec("lung", "lung", "lung", "lung"),
    TissueSpec("retina", "retina", "retina", None, note="No matching TMS FACS tissue in current h5ad"),
]

MUSCLE_SUBTYPE_PATTERNS = {
    "soleus": ("soleus",),
    "quadriceps": ("quadriceps",),
    "gastrocnemius": ("gastrocnemius",),
    "edl": ("extensor digitorum longus", "edl"),
    "tibialis_anterior": ("tibialis anterior",),
}


def slugify(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_") or "value"


def read_api_metadata(path: str | Path, refresh: bool = False, timeout: int = 180) -> pd.DataFrame:
    path = Path(path)
    if refresh or not path.exists():
        output_dir = path.parent if path.name.endswith(".tsv") else Path(DEFAULT_API_DIR)
        metadata = discover_metadata(timeout)
        output_dir.mkdir(parents=True, exist_ok=True)
        metadata.to_csv(path, sep="\t", index=False)
        return metadata
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def tms_tissue_counts(tms_h5ad: str | Path) -> pd.DataFrame:
    anndata = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    adata = anndata.read_h5ad(tms_h5ad, backed="r")
    counts = (
        adata.obs["tissue"]
        .astype(str)
        .value_counts()
        .rename_axis("tms_tissue")
        .rename("cells")
        .reset_index()
    )
    try:
        adata.file.close()
    except AttributeError:
        pass
    return counts


def discover_muscle_subtype_specs(metadata: pd.DataFrame) -> list[TissueSpec]:
    muscle = metadata.loc[metadata["tissue_final"].eq("skeletal_muscle")].copy()
    material_values = sorted(
        value
        for value in muscle["study.characteristics.material type"].dropna().astype(str).unique()
        if value
    )
    specs: list[TissueSpec] = []
    for slug, patterns in MUSCLE_SUBTYPE_PATTERNS.items():
        matched_terms = tuple(
            value
            for value in material_values
            if any(pattern in value.lower() for pattern in patterns)
        )
        if matched_terms:
            specs.append(
                TissueSpec(
                    slug=f"skeletal_muscle_{slug}",
                    label=f"skeletal muscle: {slug.replace('_', ' ')}",
                    tissue_final="skeletal_muscle",
                    tms_tissue="limb muscle",
                    material_terms=matched_terms,
                    note="OSDR sub-tissue run; TMS pretraining uses combined limb muscle FACS",
                )
            )
    return specs


def all_specs(metadata: pd.DataFrame) -> list[TissueSpec]:
    return [*BASE_SPECS, *discover_muscle_subtype_specs(metadata)]


def spec_map(metadata: pd.DataFrame) -> dict[str, TissueSpec]:
    return {spec.slug: spec for spec in all_specs(metadata)}


def select_spec_metadata(metadata: pd.DataFrame, spec: TissueSpec) -> pd.DataFrame:
    selected = metadata.loc[
        metadata["condition_inferred"].isin(["flight", "ground_control"])
        & metadata["tissue_final"].eq(spec.tissue_final)
    ].copy()
    if spec.material_terms:
        selected = selected.loc[
            selected["study.characteristics.material type"].astype(str).isin(spec.material_terms)
        ].copy()
    selected["condition_label"] = selected["condition_inferred"].map(
        {"flight": "FLT", "ground_control": "GC"}
    )
    selected["h5_accession"] = selected["id.accession"].astype(str)
    selected["profile"] = selected["id.sample name"].astype(str)
    selected["profile_id"] = (
        selected["id.accession"].astype(str) + "/" + selected["id.sample name"].astype(str)
    )
    return selected.reset_index(drop=True)


def condition_counts(selected: pd.DataFrame) -> pd.DataFrame:
    if selected.empty:
        return pd.DataFrame(columns=["id.accession", "FLT", "GC", "total"])
    counts = (
        selected.groupby(["id.accession", "condition_label"])
        .size()
        .unstack(fill_value=0)
    )
    for column in ["FLT", "GC"]:
        if column not in counts:
            counts[column] = 0
    counts["total"] = counts["FLT"] + counts["GC"]
    return counts.reset_index().sort_values(["total", "FLT", "GC"], ascending=False)


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "No rows."
    display = frame.copy()
    columns = [str(column) for column in display.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in display.iterrows():
        values = []
        for column in display.columns:
            value = str(row[column]).replace("\n", " ")
            values.append(value.replace("|", "\\|"))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_audit(args: argparse.Namespace) -> Path:
    metadata = read_api_metadata(args.metadata, refresh=args.refresh_metadata, timeout=args.timeout)
    output_dir = Path(args.output_dir) / "input_audit"
    output_dir.mkdir(parents=True, exist_ok=True)

    specs = all_specs(metadata)
    tms_counts = tms_tissue_counts(args.tms_h5ad)
    tms_count_map = dict(zip(tms_counts["tms_tissue"], tms_counts["cells"]))
    tms_counts.to_csv(output_dir / "tms_facs_tissue_counts.tsv", sep="\t", index=False)

    tissue_rows = []
    accession_tables = []
    for spec in specs:
        selected = select_spec_metadata(metadata, spec)
        counts = condition_counts(selected)
        if not counts.empty:
            counts.insert(0, "tissue_slug", spec.slug)
            accession_tables.append(counts)
        tms_cells = int(tms_count_map.get(spec.tms_tissue, 0)) if spec.tms_tissue else 0
        tissue_rows.append(
            {
                "tissue_slug": spec.slug,
                "label": spec.label,
                "tissue_final": spec.tissue_final,
                "material_terms": "; ".join(spec.material_terms),
                "tms_tissue": spec.tms_tissue or "",
                "tms_cells": tms_cells,
                "space_flight": int(selected["condition_label"].eq("FLT").sum()),
                "ground_control": int(selected["condition_label"].eq("GC").sum()),
                "accessions": int(selected["id.accession"].nunique()),
                "mober_eligible_ge2_studies": bool(selected["id.accession"].nunique() >= 2),
                "pretraining_status": (
                    "ok"
                    if tms_cells >= args.min_tms_cells
                    else "skip_no_matching_or_too_few_tms_cells"
                ),
                "note": spec.note,
            }
        )

    tissue_table = pd.DataFrame(tissue_rows)
    tissue_table.to_csv(output_dir / "requested_tissue_input_status.tsv", sep="\t", index=False)
    if accession_tables:
        pd.concat(accession_tables, ignore_index=True).to_csv(
            output_dir / "requested_tissue_accession_counts.tsv",
            sep="\t",
            index=False,
        )

    muscle_material = (
        metadata.loc[metadata["tissue_final"].eq("skeletal_muscle")]
        .groupby(["study.characteristics.material type", "condition_inferred"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    for column in ["flight", "ground_control"]:
        if column not in muscle_material:
            muscle_material[column] = 0
    muscle_material["total"] = muscle_material["flight"] + muscle_material["ground_control"]
    muscle_material.sort_values("total", ascending=False).to_csv(
        output_dir / "skeletal_muscle_material_counts.tsv",
        sep="\t",
    )

    summary_lines = [
        "# Multi-Tissue API GLARE Input Audit",
        "",
        "This audit uses NASA OSDR Biological Data API metadata and TMS FACS h5ad metadata.",
        "",
        f"- API metadata: `{args.metadata}`",
        f"- TMS h5ad: `{args.tms_h5ad}`",
        f"- Minimum TMS cells: {args.min_tms_cells}",
        "",
        "## Requested Tissue Status",
        "",
        markdown_table(tissue_table),
        "",
        "## Notes",
        "",
        "- MOBER is marked eligible when at least two OSDR studies are available.",
        "- Skeletal muscle sub-tissue runs use exact OSDR material-type labels.",
        "- Retina has OSDR bulk FLT/GC data but no matching TMS FACS tissue in the current h5ad, so it is skipped for GLARE pretraining.",
    ]
    (output_dir / "AUDIT_SUMMARY.md").write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )
    return output_dir


def ensure_count_tables(
    metadata: pd.DataFrame,
    output_dir: Path,
    timeout: int,
    download_missing: bool,
) -> dict[str, Path]:
    counts_dir = output_dir / "counts"
    paths = {
        accession: counts_dir / f"{accession}_unnormalized_counts.csv"
        for accession in sorted(metadata["id.accession"].astype(str).unique())
    }
    missing = [accession for accession, path in paths.items() if not path.exists()]
    if missing and download_missing:
        download_count_tables(
            metadata,
            output_dir,
            timeout,
            accessions=missing,
            overwrite=False,
        )
    missing = [accession for accession, path in paths.items() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing API count CSVs. Re-run with --download-counts or fetch them first: "
            + ", ".join(missing[:20])
        )
    return paths


def read_count_csv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    gene_column = frame.columns[0]
    frame = frame.rename(columns={gene_column: "gene_id"})
    frame["gene_id"] = frame["gene_id"].astype(str)
    value_columns = [column for column in frame.columns if column != "gene_id"]
    frame[value_columns] = frame[value_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    if frame["gene_id"].duplicated().any():
        frame = frame.groupby("gene_id", as_index=False)[value_columns].sum()
    return frame.set_index("gene_id")


def sample_column_map(columns: list[str]) -> dict[str, str]:
    mapping = {}
    for column in columns:
        suffix = str(column).split("/")[-1]
        if suffix not in mapping:
            mapping[suffix] = column
    return mapping


def biological_sample_name(sample: str) -> str:
    return re.sub(r"_techrep\d+$", "", str(sample))


def log2_cpm(matrix: np.ndarray) -> np.ndarray:
    matrix = matrix.astype(np.float32, copy=False)
    library_sizes = matrix.sum(axis=0, keepdims=True)
    library_sizes[library_sizes <= 0] = 1.0
    cpm = matrix / library_sizes * 1_000_000.0
    return np.log2(cpm + 1.0).astype(np.float32, copy=False)


def build_api_expression_bundles(
    selected: pd.DataFrame,
    output_dir: Path,
    api_output_dir: Path,
    timeout: int,
    download_counts: bool,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    count_paths = ensure_count_tables(selected, api_output_dir, timeout, download_counts)
    selected = selected.sort_values(["id.accession", "condition_label", "profile"]).copy()
    accession_order = selected["id.accession"].drop_duplicates().astype(str).tolist()

    blocks = []
    retained_rows = []
    missing_rows = []
    for accession in accession_order:
        accession_rows = selected.loc[selected["id.accession"].eq(accession)].copy()
        accession_rows["biological_sample"] = accession_rows["profile"].map(
            biological_sample_name
        )
        table = read_count_csv(count_paths[accession])
        colmap = sample_column_map([str(column) for column in table.columns])
        block_columns = []
        for biological_sample, group in accession_rows.groupby(
            "biological_sample", sort=False
        ):
            keep_columns = []
            missing_samples = []
            for sample in group["profile"].astype(str):
                column = colmap.get(sample)
                if column is None:
                    missing_samples.append(sample)
                else:
                    keep_columns.append(column)
            if missing_samples:
                for sample in missing_samples:
                    missing_rows.append(
                        {
                            "id.accession": accession,
                            "id.sample name": sample,
                            "count_csv": str(count_paths[accession]),
                        }
                    )
                continue
            if not keep_columns:
                continue
            feature = f"{accession}__{clean_token(biological_sample)}"
            summed = table.loc[:, keep_columns].sum(axis=1)
            summed.name = feature
            block_columns.append(summed)
            retained = group.iloc[0].to_dict()
            retained["feature"] = feature
            retained["profile"] = biological_sample
            retained["id.sample name"] = biological_sample
            retained["count_column"] = ";".join(keep_columns)
            retained["technical_replicate_count"] = int(len(keep_columns))
            retained["technical_replicate_samples"] = ";".join(
                group["profile"].astype(str).tolist()
            )
            retained_rows.append(retained)
        if block_columns:
            blocks.append(pd.concat(block_columns, axis=1))

    if missing_rows:
        pd.DataFrame(missing_rows).to_csv(
            output_dir / "missing_count_columns.tsv",
            sep="\t",
            index=False,
        )
    if not blocks:
        raise ValueError("No API count columns matched selected metadata")

    common_genes = list(blocks[0].index)
    common_gene_set = set(common_genes)
    for block in blocks[1:]:
        common_gene_set &= set(block.index)
    common_genes = [gene for gene in common_genes if gene in common_gene_set]
    if not common_genes:
        raise ValueError("No shared genes across selected API count tables")
    raw = pd.concat([block.loc[common_genes] for block in blocks], axis=1)
    raw = raw.astype(np.float32)
    log_expr = pd.DataFrame(
        log2_cpm(raw.to_numpy(dtype=np.float32, copy=False)),
        index=raw.index,
        columns=raw.columns,
    )

    retained = pd.DataFrame(retained_rows)
    retained = retained.set_index("feature").loc[raw.columns]
    retained.index.name = "feature"
    retained = retained.reset_index()
    if "feature" not in retained.columns and "index" in retained.columns:
        retained = retained.rename(columns={"index": "feature"})
    retained["profile"] = retained["id.sample name"].astype(str)
    retained["condition_label"] = retained["condition_inferred"].map(
        {"flight": "FLT", "ground_control": "GC"}
    )
    retained["condition"] = retained["condition_inferred"].map(
        {"flight": "flight", "ground_control": "ground"}
    )
    retained["sample"] = retained["feature"]
    retained["technical_replicate_group"] = retained["feature"]
    retained.to_csv(output_dir / "retained_profile_features.tsv", sep="\t", index=False)

    raw_manifest = write_matrix_bundle(
        output_dir / "api_raw_counts",
        raw.to_numpy(dtype=np.float32, copy=False),
        genes=raw.index.astype(str).tolist(),
        profiles=raw.columns.astype(str).tolist(),
        profile_metadata=retained,
        description="NASA OSDR API unnormalized counts",
    )
    log_manifest = write_matrix_bundle(
        output_dir / "api_log2_cpm",
        log_expr.to_numpy(dtype=np.float32, copy=False),
        genes=log_expr.index.astype(str).tolist(),
        profiles=log_expr.columns.astype(str).tolist(),
        profile_metadata=retained,
        description="NASA OSDR API log2(CPM+1) expression",
    )

    raw_counts_dir = output_dir / "raw_deseq2_inputs"
    raw_counts_dir.mkdir(parents=True, exist_ok=True)
    raw_for_r = raw.copy()
    raw_for_r.index.name = "gene_id"
    raw_for_r.round().astype(np.int64).to_csv(raw_counts_dir / "counts.tsv", sep="\t")
    metadata_for_r = retained[
        [
            "sample",
            "id.accession",
            "condition",
            "condition_label",
            "profile",
            "study.characteristics.material type",
            "study.characteristics.sex",
            "study.characteristics.strain",
            "study.characteristics.genotype",
            "tissue_final",
            "tissue_source",
        ]
    ].rename(columns={"id.accession": "accession"})
    metadata_for_r["stratum"] = "all"
    metadata_for_r.to_csv(raw_counts_dir / "sample_metadata.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": raw.index.astype(str), "gene_symbol": ""}).to_csv(
        raw_counts_dir / "gene_symbols.tsv",
        sep="\t",
        index=False,
    )
    condition_table = condition_counts(retained)
    condition_table.to_csv(output_dir / "condition_counts.tsv", sep="\t", index=False)

    summary = {
        "raw_manifest": str(raw_manifest),
        "log2_cpm_manifest": str(log_manifest),
        "raw_deseq2_inputs": {
            "counts": str(raw_counts_dir / "counts.tsv"),
            "metadata": str(raw_counts_dir / "sample_metadata.tsv"),
            "gene_symbols": str(raw_counts_dir / "gene_symbols.tsv"),
        },
        "genes": int(raw.shape[0]),
        "samples": int(raw.shape[1]),
        "accessions": accession_order,
        "missing_count_columns": len(missing_rows),
        "input_kind": "api_log2_cpm_expression_for_glare_and_raw_counts_for_deseq2",
    }
    (output_dir / "api_expression_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def write_controlled_target_from_bundle(target_manifest: str | Path, output_dir: Path) -> dict:
    bundle = load_matrix_bundle(target_manifest)
    if bundle.profile_metadata is None:
        raise ValueError("Aligned target bundle is missing profile metadata")
    metadata = bundle.profile_metadata.copy()
    if "condition_label" not in metadata.columns:
        raise ValueError("Target metadata is missing condition_label")
    matrix = dense_matrix(bundle.matrix)
    genes = [str(gene) for gene in bundle.genes]
    matrices = {}
    features = {}
    retained_rows = []
    profile_to_index = {str(profile): index for index, profile in enumerate(bundle.profiles)}
    for location in ("FLT", "GC"):
        rows = metadata.loc[metadata["condition_label"].eq(location)].copy()
        profiles = rows["feature"].astype(str).tolist() if "feature" in rows else rows["sample"].astype(str).tolist()
        indices = [profile_to_index[profile] for profile in profiles]
        matrices[location] = matrix[:, indices].astype(np.float32, copy=False)
        features[location] = profiles
        rows.insert(0, "location", location)
        retained_rows.append(rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / "controlled_target.npz",
        flt=matrices["FLT"],
        gc=matrices["GC"],
        genes=np.asarray(genes, dtype=str),
        flt_features=np.asarray(features["FLT"], dtype=str),
        gc_features=np.asarray(features["GC"], dtype=str),
        input_kind=np.asarray("api_log2_cpm_expression"),
        input_path=np.asarray(str(target_manifest)),
    )
    retained = pd.concat(retained_rows, ignore_index=True)
    retained.to_csv(output_dir / "retained_profile_features.tsv", sep="\t", index=False)
    return {
        "genes": genes,
        "matrices": matrices,
        "features": features,
        "target_manifest": str(target_manifest),
        "retained_profile_features": str(output_dir / "retained_profile_features.tsv"),
    }


def prepare_tms_bundle(spec: TissueSpec, output_dir: Path, tms_h5ad: str | Path) -> Path | None:
    if not spec.tms_tissue:
        return None
    pretrain_dir = output_dir / "pretrain"
    manifest = pretrain_dir / "tms_facs.manifest.json"
    if manifest.exists():
        return manifest
    return prepare_tms(
        input_h5ad=tms_h5ad,
        output_prefix=pretrain_dir / "tms_facs",
        matrix_source="X",
        obs_filters=[("tissue", spec.tms_tissue)],
    )


def prepare_scope(
    spec: TissueSpec,
    selected: pd.DataFrame,
    pretrain_manifest: Path,
    scope_dir: Path,
    api_output_dir: Path,
    timeout: int,
    download_counts: bool,
) -> dict:
    inputs_dir = scope_dir / "inputs"
    expr_summary = build_api_expression_bundles(
        selected,
        inputs_dir,
        api_output_dir,
        timeout,
        download_counts,
    )
    aligned_prefix = inputs_dir / "aligned_tms_api"
    aligned_pretrain, aligned_target = align_bundles(
        pretrain_manifest,
        expr_summary["log2_cpm_manifest"],
        aligned_prefix,
        prefer="pretrain",
    )
    mtx_path = inputs_dir / "tms_facs_pretrain.mtx"
    export_mtx(aligned_pretrain, mtx_path)
    controlled = write_controlled_target_from_bundle(aligned_target, scope_dir)
    counts = pd.read_csv(inputs_dir / "condition_counts.tsv", sep="\t")
    scope_summary = {
        "tissue_slug": spec.slug,
        "tissue_label": spec.label,
        "scope_dir": str(scope_dir),
        "pretrain_manifest": str(aligned_pretrain),
        "pretrain_mtx": str(mtx_path),
        "target_manifest": str(aligned_target),
        "controlled_target": str(scope_dir / "controlled_target.npz"),
        "expression_summary": expr_summary,
        "condition_counts": counts.to_dict(orient="records"),
        "genes": len(controlled["genes"]),
        "flt_profiles": int(controlled["matrices"]["FLT"].shape[1]),
        "gc_profiles": int(controlled["matrices"]["GC"].shape[1]),
    }
    (scope_dir / "input_summary.json").write_text(
        json.dumps(scope_summary, indent=2) + "\n",
        encoding="utf-8",
    )
    return scope_summary


def write_skip(output_dir: Path, spec: TissueSpec, reason: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "tissue_slug": spec.slug,
        "tissue_label": spec.label,
        "status": "skipped",
        "reason": reason,
        "spec": spec.__dict__,
    }
    (output_dir / "SKIPPED.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    (output_dir / "SKIPPED.md").write_text(
        f"# {spec.label}\n\nSkipped: {reason}\n",
        encoding="utf-8",
    )


def prepare_spec(args: argparse.Namespace, metadata: pd.DataFrame, spec: TissueSpec) -> dict:
    root = Path(args.output_dir) / spec.slug
    tms_counts = tms_tissue_counts(args.tms_h5ad)
    tms_count_map = dict(zip(tms_counts["tms_tissue"], tms_counts["cells"]))
    tms_cells = int(tms_count_map.get(spec.tms_tissue, 0)) if spec.tms_tissue else 0
    if tms_cells < args.min_tms_cells:
        reason = (
            f"TMS FACS tissue {spec.tms_tissue!r} has {tms_cells} cells; "
            f"minimum is {args.min_tms_cells}"
        )
        write_skip(root, spec, reason)
        return {"tissue_slug": spec.slug, "status": "skipped", "reason": reason}

    selected = select_spec_metadata(metadata, spec)
    if selected.empty:
        reason = "No OSDR API FLT/GC metadata rows matched this tissue spec"
        write_skip(root, spec, reason)
        return {"tissue_slug": spec.slug, "status": "skipped", "reason": reason}

    selected, ercc_summary = apply_ercc_policy(selected, root / "ercc_audit", args.ercc_policy)
    if args.accessions:
        selected = selected.loc[
            selected["id.accession"].astype(str).isin(set(args.accessions))
        ].copy()
        if selected.empty:
            reason = f"No selected rows remained after --accessions filter: {args.accessions}"
            write_skip(root, spec, reason)
            return {"tissue_slug": spec.slug, "status": "skipped", "reason": reason}
    pretrain_manifest = prepare_tms_bundle(spec, root, args.tms_h5ad)
    if pretrain_manifest is None:
        reason = "No matching TMS FACS tissue was configured"
        write_skip(root, spec, reason)
        return {"tissue_slug": spec.slug, "status": "skipped", "reason": reason}

    count_table = condition_counts(selected)
    count_table.to_csv(root / "selected_accession_counts.tsv", sep="\t", index=False)
    selected.to_csv(root / "selected_api_metadata.tsv", sep="\t", index=False)

    accessions_with_both = count_table.loc[
        count_table["FLT"].ge(args.min_profiles_per_condition)
        & count_table["GC"].ge(args.min_profiles_per_condition),
        "id.accession",
    ].astype(str).tolist()
    if not accessions_with_both:
        reason = (
            "No accessions had at least "
            f"{args.min_profiles_per_condition} FLT and GC profiles"
        )
        write_skip(root, spec, reason)
        return {"tissue_slug": spec.slug, "status": "skipped", "reason": reason}

    aggregate_selected = selected.loc[selected["id.accession"].isin(accessions_with_both)].copy()
    aggregate_summary = prepare_scope(
        spec,
        aggregate_selected,
        pretrain_manifest,
        root / "aggregate",
        Path(args.api_output_dir),
        args.timeout,
        args.download_counts,
    )

    per_study_summaries = []
    if args.prepare_per_study:
        for accession in accessions_with_both:
            study_selected = selected.loc[selected["id.accession"].eq(accession)].copy()
            per_study_summaries.append(
                prepare_scope(
                    spec,
                    study_selected,
                    pretrain_manifest,
                    root / "per_study" / accession,
                    Path(args.api_output_dir),
                    args.timeout,
                    args.download_counts,
                )
            )

    summary = {
        "tissue_slug": spec.slug,
        "tissue_label": spec.label,
        "status": "prepared",
        "spec": spec.__dict__,
        "tms_cells": tms_cells,
        "ercc_policy": ercc_summary,
        "accessions_with_both_conditions": accessions_with_both,
        "mober_eligible_ge2_studies": len(accessions_with_both) >= 2,
        "aggregate": aggregate_summary,
        "per_study": per_study_summaries,
    }
    (root / "PREP_SUMMARY.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    write_prepare_report(root, summary)
    return summary


def write_prepare_report(root: Path, summary: dict) -> None:
    lines = [
        f"# {summary['tissue_label']} API GLARE Preparation",
        "",
        f"- Status: `{summary['status']}`",
        f"- TMS FACS cells: {summary['tms_cells']:,}",
        f"- MOBER eligible: {summary['mober_eligible_ge2_studies']}",
        f"- Accessions with usable FLT/GC counts: {', '.join(summary['accessions_with_both_conditions'])}",
        "",
        "## Aggregate",
        "",
        f"- Genes: {summary['aggregate']['genes']:,}",
        f"- FLT profiles: {summary['aggregate']['flt_profiles']:,}",
        f"- GC profiles: {summary['aggregate']['gc_profiles']:,}",
        f"- Controlled target: `{summary['aggregate']['controlled_target']}`",
        f"- Pretrain MTX: `{summary['aggregate']['pretrain_mtx']}`",
        "",
        "## Per-Study Inputs",
        "",
    ]
    if summary["per_study"]:
        table = pd.DataFrame(
            {
                "scope": item["scope_dir"],
                "genes": item["genes"],
                "FLT": item["flt_profiles"],
                "GC": item["gc_profiles"],
            }
            for item in summary["per_study"]
        )
        lines.append(markdown_table(table))
    else:
        lines.append("Per-study inputs were not requested.")
    (root / "PREP_SUMMARY.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def parse_spec_selection(values: list[str] | None, specs: dict[str, TissueSpec]) -> list[TissueSpec]:
    if not values or values == ["all"]:
        return list(specs.values())
    selected = []
    for value in values:
        for token in re.split(r"[,\s]+", value):
            if not token:
                continue
            if token not in specs:
                raise SystemExit(f"Unknown tissue spec {token!r}. Available: {', '.join(sorted(specs))}")
            selected.append(specs[token])
    seen = set()
    unique = []
    for spec in selected:
        if spec.slug not in seen:
            unique.append(spec)
            seen.add(spec.slug)
    return unique


def run_prepare(args: argparse.Namespace) -> Path:
    start = time.perf_counter()
    metadata = read_api_metadata(args.metadata, refresh=args.refresh_metadata, timeout=args.timeout)
    specs = spec_map(metadata)
    selected_specs = parse_spec_selection(args.tissue, specs)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for spec in selected_specs:
        log(f"Preparing {spec.slug}")
        summaries.append(prepare_spec(args, metadata, spec))
    manifest = {
        "status": "prepared",
        "output_dir": str(output_dir),
        "metadata": args.metadata,
        "api_output_dir": args.api_output_dir,
        "tissues": summaries,
        "elapsed_seconds": round(time.perf_counter() - start, 3),
        "elapsed": format_elapsed(time.perf_counter() - start),
    }
    (output_dir / "PREP_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return output_dir


def run_glare_scope(
    scope_dir: Path,
    pretrained_weights: Path,
    epochs: int,
    batch_size: int,
    seed: int,
    skip_clustering: bool,
) -> dict:
    target = np.load(scope_dir / "controlled_target.npz")
    genes = target["genes"].astype(str).tolist()
    matrices = {"FLT": target["flt"], "GC": target["gc"]}
    input_dim = infer_pretrain_input_dim(pretrained_weights)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    for location, matrix in matrices.items():
        write_outlier_audit(matrix, genes, location, scope_dir)
    locations = [
        finetune_location(
            matrices[location],
            genes,
            location,
            pretrained_weights,
            scope_dir,
            device,
            input_dim,
            epochs,
            batch_size,
            seed,
        )
        for location in ("FLT", "GC")
    ]
    if not skip_clustering:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "nasa_mouse_glare.paper_clustering",
                "--run-dir",
                str(scope_dir),
                "--skip-tsne",
            ],
            check=True,
        )
    return {
        "scope_dir": str(scope_dir),
        "pretrained_weights": str(pretrained_weights),
        "pretrained_input_dim": input_dim,
        "device": str(device),
        "locations": locations,
    }


def infer_tissue_root(scope_dir: Path) -> Path | None:
    for candidate in [scope_dir, *scope_dir.parents]:
        if (candidate / "pretrain" / "tms_facs.manifest.json").exists():
            return candidate
    return None


def default_pretraining_dir(scope_dir: Path) -> Path:
    tissue_root = infer_tissue_root(scope_dir)
    if tissue_root is not None:
        return tissue_root / "glare_pretraining"
    return scope_dir / "pretraining"


def run_glare(args: argparse.Namespace) -> Path:
    root = Path(args.scope_dir)
    input_summary = json.loads((root / "input_summary.json").read_text(encoding="utf-8"))
    pretraining_dir = (
        Path(args.pretraining_dir)
        if args.pretraining_dir
        else default_pretraining_dir(root)
    )
    pretrained_weights = (
        Path(args.pretrained_weights) if args.pretrained_weights else pretraining_dir / "sc_shulse_pretrained_reproduced.pth"
    )
    if not pretrained_weights.exists():
        pretraining_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                sys.executable,
                "-m",
                "nasa_mouse_glare.reproduce_glare_pretrain",
                "--input",
                input_summary["pretrain_mtx"],
                "--output-dir",
                str(pretraining_dir),
                "--epochs",
                str(args.pretrain_epochs),
                "--batch-size",
                str(args.batch_size),
                "--seed",
                str(args.seed),
                "--num-workers",
                "0",
            ],
            check=True,
        )
    summary = run_glare_scope(
        root,
        pretrained_weights,
        args.finetune_epochs,
        args.batch_size,
        args.seed,
        args.skip_clustering,
    )
    (root / "finetune_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    return root


def bundle_to_mober_h5ad(
    target_manifest: str | Path,
    output_path: Path,
    batch_column: str,
) -> dict:
    anndata = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    bundle = load_matrix_bundle(target_manifest)
    if bundle.profile_metadata is None:
        raise ValueError("Target bundle is missing profile metadata")
    metadata = bundle.profile_metadata.copy()
    metadata.index = pd.Index([str(profile) for profile in bundle.profiles], dtype=object)
    metadata.index.name = "sample"
    if "feature" not in metadata.columns:
        metadata.insert(0, "feature", metadata.index.astype(str))
    if "profile" not in metadata.columns:
        metadata.insert(1, "profile", metadata.index.astype(str))
    if batch_column not in metadata.columns:
        raise ValueError(f"MOBER batch column not found: {batch_column}")
    if "condition_label" not in metadata.columns:
        raise ValueError("Target bundle metadata is missing condition_label")
    metadata["data_source"] = metadata[batch_column].astype(str)
    metadata["location"] = metadata["condition_label"].astype(str)
    metadata["condition"] = metadata["location"]
    for column in metadata.columns:
        metadata[column] = metadata[column].map(
            lambda value: "" if pd.isna(value) else str(value)
        ).astype(object)

    x = dense_matrix(bundle.matrix).T.astype(np.float32, copy=False)
    var = pd.DataFrame(
        {"gene_id": [str(gene) for gene in bundle.genes]},
        index=pd.Index([str(gene) for gene in bundle.genes], dtype=object),
    )
    adata = anndata.AnnData(X=x, obs=metadata, var=var)
    adata.uns["normalization"] = "log2(CPM+1) from NASA OSDR API raw counts"
    adata.uns["source_target_manifest"] = str(target_manifest)
    adata.uns["batch_column"] = batch_column
    sanitize_anndata_strings(adata)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    return {
        "input_h5ad": str(output_path),
        "target_manifest": str(target_manifest),
        "shape_samples_x_genes": list(adata.shape),
        "batch_column": batch_column,
        "data_sources": sorted(metadata["data_source"].unique().tolist()),
    }


def choose_mober_onto(h5ad_path: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    anndata = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    adata = anndata.read_h5ad(h5ad_path, backed="r")
    counts = adata.obs["data_source"].astype(str).value_counts()
    try:
        adata.file.close()
    except AttributeError:
        pass
    if counts.empty:
        raise ValueError("Cannot choose MOBER projection target without data_source values")
    return str(counts.index[0])


def train_mober(
    train_file: Path,
    output_dir: Path,
    args: argparse.Namespace,
) -> Path:
    ensure_mober_importable()
    from mober.core import train as mober_train

    train_output = output_dir / "mober_train"
    train_args = SimpleNamespace(
        train_file=str(train_file),
        use_sparse_mat=False,
        src_adv_weight=args.src_adv_weight,
        src_adv_lr=args.src_adv_lr,
        batch_ae_lr=args.batch_ae_lr,
        val_set_size=args.val_set_size,
        encoding_dim=args.mober_encoding_dim,
        balanced_sources_ae=False,
        balanced_sources_src_adv=args.balanced_sources_src_adv,
        batch_size=args.mober_batch_size,
        epochs=args.mober_epochs,
        random_seed=args.seed,
        kl_weight=args.kl_weight,
        patience=args.mober_patience,
        output_dir=str(train_output),
        use_mlflow=False,
        mlflow_storage_path="",
        experiment_name="mober",
        run_name="run",
        tmp_dir="tmp",
    )
    log(f"Training MOBER: {train_file} -> {train_output}")
    mober_train.main(train_args)
    return train_output / "models"


def project_mober(
    train_file: Path,
    model_dir: Path,
    output_dir: Path,
    onto: str,
    decimals: int,
    projection_batch_size: int,
) -> dict:
    ensure_mober_importable()
    import anndata as ad
    import torch as torch_module
    from mober.core.projection import do_projection, load_model

    projection_dir = output_dir / "projection"
    projection_dir.mkdir(parents=True, exist_ok=True)
    device = torch_module.device("cuda" if torch_module.cuda.is_available() else "cpu")
    model, features, label_encode = load_model(model_dir, device)
    adata = ad.read_h5ad(train_file)
    adata = adata[:, features].copy()
    proj_adata, z_adata = do_projection(
        model,
        adata,
        onto,
        label_encode,
        device,
        decimals=decimals,
        batch_size=projection_batch_size,
    )
    sanitize_anndata_strings(proj_adata)
    sanitize_anndata_strings(z_adata)
    projected_path = projection_dir / f"mober_projected_onto_{onto}.h5ad"
    latent_path = projection_dir / f"mober_latent_onto_{onto}.h5ad"
    proj_adata.write_h5ad(projected_path)
    z_adata.write_h5ad(latent_path)
    latent = pd.DataFrame(
        z_adata.X,
        index=z_adata.obs_names,
        columns=z_adata.var_names,
    )
    latent.insert(0, "sample", latent.index)
    for column in ["data_source", "location", "condition", "tissue_final", "profile"]:
        if column in z_adata.obs:
            latent.insert(1, column, z_adata.obs[column].astype(str).to_numpy())
    latent_tsv = projection_dir / f"mober_latent_onto_{onto}.tsv"
    latent.to_csv(latent_tsv, sep="\t", index=False)
    summary = {
        "projected_h5ad": str(projected_path),
        "latent_h5ad": str(latent_path),
        "latent_tsv": str(latent_tsv),
        "model_dir": str(model_dir),
        "projection_file": str(train_file),
        "onto": onto,
        "shape_samples_x_genes": list(proj_adata.shape),
        "latent_shape": list(z_adata.shape),
    }
    (projection_dir / "mober_projection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def run_mober_scope(args: argparse.Namespace) -> Path:
    source_scope = Path(args.scope_dir)
    source_summary = json.loads(
        (source_scope / "input_summary.json").read_text(encoding="utf-8")
    )
    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else source_scope.parent / f"{source_scope.name}_mober"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    mober_dir = output_dir / "mober"
    h5ad_path = mober_dir / "mober_input.h5ad"
    prep = bundle_to_mober_h5ad(
        source_summary["target_manifest"],
        h5ad_path,
        args.batch_column,
    )
    onto = choose_mober_onto(h5ad_path, args.onto)
    model_dir = (
        Path(args.mober_model_dir)
        if args.mober_model_dir
        else train_mober(h5ad_path, mober_dir, args)
    )
    projection = project_mober(
        h5ad_path,
        model_dir,
        mober_dir,
        onto,
        args.decimals,
        args.projection_batch_size,
    )
    prepared = prepare_mober_target(projection["projected_h5ad"], output_dir)
    counts = prepared["counts"]
    counts.to_csv(output_dir / "condition_counts.tsv", sep="\t", index=False)
    input_summary = {
        "source_scope_dir": str(source_scope),
        "mober_dir": str(mober_dir),
        "mober_prepare": prep,
        "mober_projection": projection,
        "pretrain_mtx": source_summary["pretrain_mtx"],
        "pretrain_manifest": source_summary["pretrain_manifest"],
        "target_manifest": projection["projected_h5ad"],
        "controlled_target": str(output_dir / "controlled_target.npz"),
        "genes": len(prepared["genes"]),
        "flt_profiles": int(prepared["matrices"]["FLT"].shape[1]),
        "gc_profiles": int(prepared["matrices"]["GC"].shape[1]),
        "input_kind": "mober_projected_log2_cpm_expression",
        "mober_projection_target": onto,
    }
    (output_dir / "input_summary.json").write_text(
        json.dumps(input_summary, indent=2) + "\n",
        encoding="utf-8",
    )
    if not args.prepare_only:
        run_args = SimpleNamespace(
            scope_dir=str(output_dir),
            pretrained_weights=args.pretrained_weights,
            pretraining_dir=(
                args.pretraining_dir
                if args.pretraining_dir
                else str(default_pretraining_dir(source_scope))
            ),
            pretrain_epochs=args.pretrain_epochs,
            finetune_epochs=args.finetune_epochs,
            batch_size=args.batch_size,
            seed=args.seed,
            skip_clustering=args.skip_clustering,
        )
        run_glare(run_args)
    return output_dir


def run_per_study_glare(args: argparse.Namespace) -> Path:
    root = Path(args.tissue_dir)
    per_study_root = root / "per_study"
    if not per_study_root.exists():
        raise FileNotFoundError(per_study_root)
    summaries = []
    for input_summary_path in sorted(per_study_root.glob("*/input_summary.json")):
        scope_dir = input_summary_path.parent
        if args.studies and scope_dir.name not in args.studies:
            continue
        log(f"Running per-study GLARE scope: {scope_dir}")
        run_args = SimpleNamespace(
            scope_dir=str(scope_dir),
            pretrained_weights=args.pretrained_weights,
            pretraining_dir=args.pretraining_dir,
            pretrain_epochs=args.pretrain_epochs,
            finetune_epochs=args.finetune_epochs,
            batch_size=args.batch_size,
            seed=args.seed,
            skip_clustering=args.skip_clustering,
        )
        run_glare(run_args)
        summaries.append(str(scope_dir / "finetune_summary.json"))
    summary = {
        "tissue_dir": str(root),
        "studies": [Path(path).parent.name for path in summaries],
        "summaries": summaries,
    }
    (root / "per_study_glare_run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    return root


def run_dgea_comparison(args: argparse.Namespace) -> Path:
    from .per_study_glare import cluster_enrichment_for_run, summarize_recurring_terms
    from .per_study_glare_dgea_comparison import (
        compare_pathways,
        default_rscript,
        run_deseq2,
        study_gene_comparison,
        write_report as write_dgea_glare_report,
    )

    tissue_dir = Path(args.tissue_dir)
    per_study_dir = tissue_dir / "per_study"
    aggregate_inputs = tissue_dir / "aggregate" / "inputs" / "raw_deseq2_inputs"
    output_dir = Path(args.output_dir) if args.output_dir else tissue_dir / "dgea_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)
    input_paths = {
        "counts": str(aggregate_inputs / "counts.tsv"),
        "metadata": str(aggregate_inputs / "sample_metadata.tsv"),
        "gene_symbols": str(aggregate_inputs / "gene_symbols.tsv"),
    }
    for path in input_paths.values():
        if not Path(path).exists():
            raise FileNotFoundError(path)
    studies = args.studies or sorted(
        path.name for path in per_study_dir.iterdir() if (path / "controlled_target.npz").exists()
    )
    if not args.skip_deseq2:
        run_deseq2(
            args.rscript or default_rscript(),
            input_paths,
            output_dir,
            args.alpha,
            args.min_count,
            args.min_samples,
        )
    deseq_path = output_dir / "deseq2" / "per_study_deseq2.tsv"
    if not deseq_path.exists():
        raise FileNotFoundError(deseq_path)
    deseq = pd.read_csv(deseq_path, sep="\t")
    deseq["gene_id"] = deseq["gene_id"].astype(str)

    enrichment_tables = []
    for accession in studies:
        run_dir = per_study_dir / accession
        if not (run_dir / "clustering" / "FLT_gene_clusters.tsv").exists():
            continue
        table = cluster_enrichment_for_run(run_dir, Path(args.reactome_gmt), args.min_overlap)
        if not table.empty:
            enrichment_tables.append(table)
    if enrichment_tables:
        enrichment = pd.concat(enrichment_tables, ignore_index=True)
    else:
        enrichment = pd.DataFrame()
    recurrence = summarize_recurring_terms(enrichment, args.alpha)
    enrichment.to_csv(
        per_study_dir / "glare_cluster_reactome_enrichment.tsv",
        sep="\t",
        index=False,
    )
    recurrence.to_csv(
        per_study_dir / "recurring_glare_reactome_terms.tsv",
        sep="\t",
        index=False,
    )

    gene_tables = []
    cluster_tables = []
    summaries = []
    for accession in studies:
        run_dir = per_study_dir / accession
        if not (run_dir / "clustering" / "FLT_gene_clusters.tsv").exists():
            continue
        gene_table, cluster_summary, study_summary = study_gene_comparison(
            per_study_dir,
            accession,
            deseq,
            args.alpha,
            args.lfc_cutoff,
        )
        gene_table.to_csv(
            output_dir / f"{accession}_gene_level_glare_dgea.tsv",
            sep="\t",
            index=False,
        )
        cluster_summary.to_csv(
            output_dir / f"{accession}_cluster_dgea_enrichment.tsv",
            sep="\t",
            index=False,
        )
        gene_tables.append(gene_table.assign(accession_for_table=accession))
        cluster_tables.append(cluster_summary)
        summaries.append(study_summary)
    all_clusters = pd.concat(cluster_tables, ignore_index=True) if cluster_tables else pd.DataFrame()
    all_clusters.to_csv(
        output_dir / "cluster_dgea_enrichment_all_studies.tsv",
        sep="\t",
        index=False,
    )
    rank, ora, rank_overlap, recurring_overlap = compare_pathways(
        deseq,
        per_study_dir,
        output_dir,
        studies,
        Path(args.reactome_gmt),
        args.alpha,
        args.min_pathway_size,
        args.max_pathway_size,
        args.min_overlap,
    )
    summary = pd.DataFrame(summaries)
    deseq_summary_path = output_dir / "deseq2" / "study_deseq2_summary.tsv"
    if deseq_summary_path.exists() and not summary.empty:
        deseq_summary = pd.read_csv(deseq_summary_path, sep="\t")
        summary = summary.merge(
            deseq_summary[
                [
                    "accession",
                    "n_flight",
                    "n_ground",
                    "genes_tested",
                    "design",
                    "dispersion_fit",
                ]
            ].rename(
                columns={
                    "design": "deseq2_design",
                    "dispersion_fit": "deseq2_dispersion_fit",
                }
            ),
            on="accession",
            how="left",
            validate="one_to_one",
        )
    if not summary.empty:
        overlap_counts = (
            rank_overlap.groupby("accession").size() if not rank_overlap.empty else pd.Series(dtype=int)
        )
        summary["rank_pathway_glare_overlaps"] = (
            summary["accession"].map(overlap_counts).fillna(0).astype(int)
        )
        summary.to_csv(output_dir / "per_study_glare_dgea_summary.tsv", sep="\t", index=False)
        write_dgea_glare_report(
            output_dir,
            summary,
            all_clusters,
            recurring_overlap,
            rank_overlap,
            input_paths,
            args.alpha,
            args.lfc_cutoff,
        )
    return output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default=DEFAULT_API_METADATA)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--tms-h5ad", default=DEFAULT_TMS_H5AD)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--refresh-metadata", action="store_true")
    parser.add_argument("--min-tms-cells", type=int, default=DEFAULT_MIN_TMS_CELLS)

    subparsers = parser.add_subparsers(dest="command", required=True)
    audit = subparsers.add_parser("audit")

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--tissue", nargs="+", default=["all"])
    prepare.add_argument("--api-output-dir", default=DEFAULT_API_DIR)
    prepare.add_argument("--download-counts", action="store_true")
    prepare.add_argument("--accessions", nargs="+")
    prepare.add_argument("--prepare-per-study", action="store_true")
    prepare.add_argument(
        "--ercc-policy",
        choices=["keep_all", "prefer_noercc"],
        default="prefer_noercc",
    )
    prepare.add_argument(
        "--min-profiles-per-condition",
        type=int,
        default=DEFAULT_MIN_PROFILES_PER_CONDITION,
    )

    glare = subparsers.add_parser("run-glare-scope")
    glare.add_argument("--scope-dir", required=True)
    glare.add_argument("--pretrained-weights")
    glare.add_argument("--pretraining-dir")
    glare.add_argument("--pretrain-epochs", type=int, default=DEFAULT_PRETRAIN_EPOCHS)
    glare.add_argument("--finetune-epochs", type=int, default=DEFAULT_FINETUNE_EPOCHS)
    glare.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    glare.add_argument("--seed", type=int, default=DEFAULT_SEED)
    glare.add_argument("--skip-clustering", action="store_true")

    mober = subparsers.add_parser("run-mober-scope")
    mober.add_argument("--scope-dir", required=True)
    mober.add_argument("--output-dir")
    mober.add_argument("--batch-column", default="id.accession")
    mober.add_argument("--onto", default="auto")
    mober.add_argument("--mober-model-dir")
    mober.add_argument("--mober-epochs", type=int, default=300)
    mober.add_argument("--mober-batch-size", type=int, default=32)
    mober.add_argument("--mober-patience", type=int, default=50)
    mober.add_argument("--mober-encoding-dim", type=int, default=64)
    mober.add_argument("--src-adv-weight", type=float, default=0.01)
    mober.add_argument("--src-adv-lr", type=float, default=1e-3)
    mober.add_argument("--batch-ae-lr", type=float, default=1e-3)
    mober.add_argument("--kl-weight", type=float, default=1e-5)
    mober.add_argument("--val-set-size", type=float, default=0.1)
    mober.add_argument("--balanced-sources-src-adv", action="store_true")
    mober.add_argument("--projection-batch-size", type=int, default=64)
    mober.add_argument("--decimals", type=int, default=4)
    mober.add_argument("--pretrained-weights")
    mober.add_argument("--pretraining-dir")
    mober.add_argument("--pretrain-epochs", type=int, default=DEFAULT_PRETRAIN_EPOCHS)
    mober.add_argument("--finetune-epochs", type=int, default=DEFAULT_FINETUNE_EPOCHS)
    mober.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    mober.add_argument("--seed", type=int, default=DEFAULT_SEED)
    mober.add_argument("--skip-clustering", action="store_true")
    mober.add_argument("--prepare-only", action="store_true")

    per_study = subparsers.add_parser("run-per-study-glare")
    per_study.add_argument("--tissue-dir", required=True)
    per_study.add_argument("--studies", nargs="+")
    per_study.add_argument("--pretrained-weights")
    per_study.add_argument("--pretraining-dir")
    per_study.add_argument("--pretrain-epochs", type=int, default=DEFAULT_PRETRAIN_EPOCHS)
    per_study.add_argument("--finetune-epochs", type=int, default=DEFAULT_FINETUNE_EPOCHS)
    per_study.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    per_study.add_argument("--seed", type=int, default=DEFAULT_SEED)
    per_study.add_argument("--skip-clustering", action="store_true")

    compare = subparsers.add_parser("run-dgea-comparison")
    compare.add_argument("--tissue-dir", required=True)
    compare.add_argument("--output-dir")
    compare.add_argument("--studies", nargs="+")
    compare.add_argument("--rscript")
    compare.add_argument("--skip-deseq2", action="store_true")
    compare.add_argument("--reactome-gmt", default=DEFAULT_REACTOME_GMT)
    compare.add_argument("--alpha", type=float, default=0.05)
    compare.add_argument("--lfc-cutoff", type=float, default=1.0)
    compare.add_argument("--min-count", type=int, default=10)
    compare.add_argument("--min-samples", type=int, default=3)
    compare.add_argument("--min-overlap", type=int, default=3)
    compare.add_argument("--min-pathway-size", type=int, default=10)
    compare.add_argument("--max-pathway-size", type=int, default=500)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "audit":
        output = write_audit(args)
    elif args.command == "prepare":
        output = run_prepare(args)
    elif args.command == "run-glare-scope":
        output = run_glare(args)
    elif args.command == "run-mober-scope":
        output = run_mober_scope(args)
    elif args.command == "run-per-study-glare":
        output = run_per_study_glare(args)
    elif args.command == "run-dgea-comparison":
        output = run_dgea_comparison(args)
    else:
        raise ValueError(args.command)
    print(output)


if __name__ == "__main__":
    main()
