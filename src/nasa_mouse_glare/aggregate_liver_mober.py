"""Prepare and run MOBER on the ribo-depletion aggregate liver OSDR set."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from .aggregate_liver_finetune import (
    DEFAULT_OSDR_H5,
    DEFAULT_TARGET_MANIFEST,
    select_aggregate_profiles,
)
from .io import require_import


DEFAULT_ACCESSIONS = [
    "OSD-379",
    "OSD-245",
    "OSD-463",
    "OSD-242",
    "OSD-137",
    "OSD-173",
]
DEFAULT_OUTPUT_DIR = "outputs/mober_liver_ribo6_osdr"
DEFAULT_BATCH_COLUMN = "h5_accession"
DEFAULT_ONTO = "OSD-379"


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def mober_source_path() -> Path:
    return Path(__file__).resolve().parents[1] / "MOBER"


def ensure_mober_importable() -> None:
    source_path = str(mober_source_path())
    if source_path not in sys.path:
        sys.path.insert(0, source_path)


def decode_array(values) -> list[str]:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8", "replace"))
        else:
            decoded.append(str(value))
    return decoded


def load_gene_symbols(osdr_h5: str | Path) -> dict[str, str]:
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")
    with h5py.File(osdr_h5, "r") as handle:
        genes = decode_array(handle["/meta/genes/ensembl_gene"][:])
        symbols = decode_array(handle["/meta/genes/symbol"][:])
    return dict(zip(genes, symbols))


def log2_cpm(matrix: np.ndarray) -> np.ndarray:
    """Convert genes x samples count-like matrix to samples x genes log2(CPM+1)."""
    matrix = matrix.astype(np.float32, copy=False)
    library_sizes = matrix.sum(axis=0, keepdims=True)
    library_sizes[library_sizes <= 0] = 1.0
    cpm = matrix / library_sizes * 1_000_000.0
    return np.log2(cpm + 1.0).T.astype(np.float32, copy=False)


def prepare(args: argparse.Namespace) -> dict:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prepared = select_aggregate_profiles(
        args.target_manifest,
        args.osdr_h5,
        args.accessions,
        output_dir,
    )
    target = np.load(output_dir / "controlled_target.npz")
    genes = target["genes"].astype(str)
    flt = target["flt"]
    gc = target["gc"]
    flt_features = target["flt_features"].astype(str)
    gc_features = target["gc_features"].astype(str)
    metadata = pd.read_csv(output_dir / "retained_profile_features.tsv", sep="\t")
    if len(metadata) != len(flt_features) + len(gc_features):
        raise ValueError(
            "MOBER metadata rows do not match FLT+GC feature count: "
            f"{len(metadata)} vs {len(flt_features) + len(gc_features)}"
        )
    if args.batch_column not in metadata.columns:
        raise ValueError(f"Batch column not found in metadata: {args.batch_column}")

    matrix = np.concatenate([flt, gc], axis=1)
    x = log2_cpm(matrix)
    metadata = metadata.copy()
    metadata.index = pd.Index(
        [str(feature) for feature in list(flt_features) + list(gc_features)],
        dtype=object,
    )
    metadata.index.name = "sample"
    metadata["data_source"] = metadata[args.batch_column].astype(str)
    metadata["condition"] = metadata["location"].astype(str)
    for column in metadata.columns:
        metadata[column] = metadata[column].map(
            lambda value: "" if pd.isna(value) else str(value)
        ).astype(object)

    symbols = load_gene_symbols(args.osdr_h5)
    var = pd.DataFrame(
        {
            "gene_id": [str(gene) for gene in genes],
            "gene_symbol": [symbols.get(gene, "") for gene in genes],
        },
        index=pd.Index([str(gene) for gene in genes], dtype=object),
    )
    for column in var.columns:
        var[column] = var[column].map(
            lambda value: "" if pd.isna(value) else str(value)
        ).astype(object)
    anndata = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    adata = anndata.AnnData(X=x, obs=metadata, var=var)
    adata.uns["normalization"] = "log2(CPM+1) from OSDR count-like HDF5 expression"
    adata.uns["batch_column"] = args.batch_column
    adata.uns["accessions"] = args.accessions
    h5ad_path = output_dir / "mober_liver_ribo6_input.h5ad"
    adata.write_h5ad(h5ad_path)

    counts = metadata.groupby(["h5_accession", "location"]).size().unstack(fill_value=0)
    summary = {
        "input_h5ad": str(h5ad_path),
        "output_dir": str(output_dir),
        "target_manifest": args.target_manifest,
        "osdr_h5": args.osdr_h5,
        "accessions": args.accessions,
        "shape_samples_x_genes": list(adata.shape),
        "normalization": adata.uns["normalization"],
        "batch_column": args.batch_column,
        "data_sources": sorted(metadata["data_source"].unique().tolist()),
        "conditions": sorted(metadata["condition"].unique().tolist()),
        "condition_counts": counts.reset_index().to_dict(orient="records"),
    }
    (output_dir / "mober_prep_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_prepare_report(output_dir, summary)
    log(f"Wrote MOBER input: {h5ad_path} ({adata.n_obs} samples x {adata.n_vars} genes)")
    return summary


def write_prepare_report(output_dir: Path, summary: dict) -> None:
    counts = pd.DataFrame(summary["condition_counts"])
    count_lines = counts.to_csv(sep="\t", index=False).strip().splitlines()
    text = [
        "# MOBER Liver Ribo-Depletion Six-Dataset Input",
        "",
        "Prepared aggregate OSDR liver FLT/GC expression for MOBER.",
        "",
        f"- Input h5ad: `{summary['input_h5ad']}`",
        f"- Shape: {summary['shape_samples_x_genes'][0]:,} samples x "
        f"{summary['shape_samples_x_genes'][1]:,} genes",
        f"- Normalization: {summary['normalization']}",
        f"- Batch/data_source column: `{summary['batch_column']}`",
        f"- Data sources: {', '.join(summary['data_sources'])}",
        "",
        "```tsv",
        *count_lines,
        "```",
    ]
    (output_dir / "MOBER_PREP_SUMMARY.md").write_text(
        "\n".join(text) + "\n", encoding="utf-8"
    )


def sanitize_anndata_strings(adata):
    adata.obs.index = pd.Index([str(value) for value in adata.obs.index], dtype=object)
    adata.var.index = pd.Index([str(value) for value in adata.var.index], dtype=object)
    for frame in (adata.obs, adata.var):
        for column in frame.columns:
            if (
                pd.api.types.is_object_dtype(frame[column])
                or pd.api.types.is_string_dtype(frame[column])
                or pd.api.types.is_categorical_dtype(frame[column])
            ):
                frame[column] = frame[column].map(
                    lambda value: "" if pd.isna(value) else str(value)
                ).astype(object)


def train(args: argparse.Namespace) -> None:
    ensure_mober_importable()
    from mober.core import train as mober_train

    output_dir = Path(args.output_dir)
    train_file = Path(args.train_file) if args.train_file else output_dir / "mober_liver_ribo6_input.h5ad"
    train_output = Path(args.mober_output_dir) if args.mober_output_dir else output_dir / "mober_train"
    train_args = SimpleNamespace(
        train_file=str(train_file),
        use_sparse_mat=False,
        src_adv_weight=args.src_adv_weight,
        src_adv_lr=args.src_adv_lr,
        batch_ae_lr=args.batch_ae_lr,
        val_set_size=args.val_set_size,
        encoding_dim=args.encoding_dim,
        balanced_sources_ae=False,
        balanced_sources_src_adv=args.balanced_sources_src_adv,
        batch_size=args.batch_size,
        epochs=args.epochs,
        random_seed=args.seed,
        kl_weight=args.kl_weight,
        patience=args.patience,
        output_dir=str(train_output),
        use_mlflow=False,
        mlflow_storage_path="",
        experiment_name="mober",
        run_name="run",
        tmp_dir="tmp",
    )
    log(f"Training MOBER: {train_file} -> {train_output}")
    mober_train.main(train_args)
    log(f"Saved MOBER model under {train_output / 'models'}")


def project(args: argparse.Namespace) -> None:
    ensure_mober_importable()
    import torch
    import anndata as ad
    from mober.core.projection import do_projection, load_model

    output_dir = Path(args.output_dir)
    train_file = Path(args.train_file) if args.train_file else output_dir / "mober_liver_ribo6_input.h5ad"
    model_dir = Path(args.model_dir) if args.model_dir else output_dir / "mober_train" / "models"
    projection_dir = output_dir / "projection"
    projection_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, features, label_encode = load_model(model_dir, device)
    adata = ad.read_h5ad(train_file)
    adata = adata[:, features].copy()
    proj_adata, z_adata = do_projection(
        model,
        adata,
        args.onto,
        label_encode,
        device,
        decimals=args.decimals,
        batch_size=args.projection_batch_size,
    )
    sanitize_anndata_strings(proj_adata)
    sanitize_anndata_strings(z_adata)
    projected_path = projection_dir / f"mober_projected_onto_{args.onto}.h5ad"
    latent_path = projection_dir / f"mober_latent_onto_{args.onto}.h5ad"
    proj_adata.write_h5ad(projected_path)
    z_adata.write_h5ad(latent_path)

    latent = pd.DataFrame(
        z_adata.X,
        index=z_adata.obs_names,
        columns=z_adata.var_names,
    )
    latent.insert(0, "sample", latent.index)
    for column in ["h5_accession", "project_identifier", "location", "condition", "sex", "strain"]:
        if column in z_adata.obs:
            latent.insert(1, column, z_adata.obs[column].astype(str).to_numpy())
    latent.to_csv(projection_dir / f"mober_latent_onto_{args.onto}.tsv", sep="\t", index=False)

    summary = {
        "projected_h5ad": str(projected_path),
        "latent_h5ad": str(latent_path),
        "latent_tsv": str(projection_dir / f"mober_latent_onto_{args.onto}.tsv"),
        "model_dir": str(model_dir),
        "projection_file": str(train_file),
        "onto": args.onto,
        "shape_samples_x_genes": list(proj_adata.shape),
        "latent_shape": list(z_adata.shape),
    }
    (projection_dir / "mober_projection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_projection_report(projection_dir, summary)
    log(f"Wrote MOBER projection: {projected_path}")
    log(f"Wrote MOBER latent embeddings: {latent_path}")


def write_projection_report(output_dir: Path, summary: dict) -> None:
    text = [
        "# MOBER Projection Summary",
        "",
        f"- Projected onto: `{summary['onto']}`",
        f"- Projected h5ad: `{summary['projected_h5ad']}`",
        f"- Latent h5ad: `{summary['latent_h5ad']}`",
        f"- Latent TSV: `{summary['latent_tsv']}`",
        f"- Projected shape: {summary['shape_samples_x_genes'][0]:,} samples x "
        f"{summary['shape_samples_x_genes'][1]:,} genes",
        f"- Latent shape: {summary['latent_shape'][0]:,} samples x "
        f"{summary['latent_shape'][1]:,} dimensions",
    ]
    (output_dir / "MOBER_PROJECTION_SUMMARY.md").write_text(
        "\n".join(text) + "\n", encoding="utf-8"
    )


def run(args: argparse.Namespace) -> None:
    prepare(args)
    train(args)
    project(args)


def add_shared_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--target-manifest", default=DEFAULT_TARGET_MANIFEST)
    parser.add_argument("--osdr-h5", default=DEFAULT_OSDR_H5)
    parser.add_argument("--accessions", nargs="+", default=DEFAULT_ACCESSIONS)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--batch-column", default=DEFAULT_BATCH_COLUMN)


def add_train_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--train-file")
    parser.add_argument("--mober-output-dir")
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--val-set-size", type=float, default=0.1)
    parser.add_argument("--patience", type=int, default=50)
    parser.add_argument("--encoding-dim", type=int, default=64)
    parser.add_argument("--src-adv-weight", type=float, default=0.01)
    parser.add_argument("--src-adv-lr", type=float, default=1e-3)
    parser.add_argument("--batch-ae-lr", type=float, default=1e-3)
    parser.add_argument("--kl-weight", type=float, default=1e-5)
    parser.add_argument("--seed", type=int, default=1996)
    parser.add_argument("--balanced-sources-src-adv", action="store_true")


def add_project_arguments(
    parser: argparse.ArgumentParser, include_train_file: bool = True
) -> None:
    if include_train_file:
        parser.add_argument("--train-file")
    parser.add_argument("--model-dir")
    parser.add_argument("--onto", default=DEFAULT_ONTO)
    parser.add_argument("--projection-batch-size", type=int, default=64)
    parser.add_argument("--decimals", type=int, default=4)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare/train/project MOBER for six ribo-depletion OSDR liver datasets."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prep_parser = subparsers.add_parser("prepare")
    add_shared_arguments(prep_parser)

    train_parser = subparsers.add_parser("train")
    train_parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    add_train_arguments(train_parser)

    project_parser = subparsers.add_parser("project")
    project_parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    add_project_arguments(project_parser)

    run_parser = subparsers.add_parser("run")
    add_shared_arguments(run_parser)
    add_train_arguments(run_parser)
    add_project_arguments(run_parser, include_train_file=False)

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "prepare":
        prepare(args)
    elif args.command == "train":
        train(args)
    elif args.command == "project":
        project(args)
    elif args.command == "run":
        run(args)
    else:
        raise ValueError(args.command)


if __name__ == "__main__":
    main()
