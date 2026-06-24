"""Run GLARE fine-tuning on MOBER-corrected aggregate liver expression."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from .paper_finetune import (
    finetune_location,
    format_elapsed,
    infer_pretrain_input_dim,
    log,
    write_outlier_audit,
)


DEFAULT_MOBER_H5AD = (
    "outputs/mober_liver_ribo6_osdr/projection/mober_projected_onto_OSD-379.h5ad"
)
DEFAULT_PRETRAINED_WEIGHTS = (
    "outputs/glare_paper_tms_liver_osd379/pretraining/sc_shulse_pretrained_reproduced.pth"
)
DEFAULT_OUTPUT_DIR = "outputs/glare_tms_liver_mober_ribo6_osdr"


def as_array(matrix) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return np.asarray(matrix)


def prepare_mober_target(
    mober_h5ad: str | Path,
    output_dir: Path,
) -> dict:
    import anndata as ad

    adata = ad.read_h5ad(mober_h5ad)
    required = {"location", "data_source"}
    missing = required - set(adata.obs.columns)
    if missing:
        raise ValueError(f"MOBER h5ad obs is missing columns: {sorted(missing)}")

    x = as_array(adata.X).astype(np.float32, copy=False)
    genes = adata.var_names.astype(str).tolist()
    metadata = adata.obs.copy()
    metadata.index = metadata.index.astype(str)
    if "feature" not in metadata.columns:
        metadata.insert(0, "feature", metadata.index)
    if "profile" not in metadata.columns:
        metadata.insert(1, "profile", metadata.index)

    matrices = {}
    features = {}
    retained_rows = []
    for location in ("FLT", "GC"):
        mask = metadata["location"].astype(str).eq(location).to_numpy()
        if not mask.any():
            raise ValueError(f"No {location} samples found in MOBER h5ad")
        matrices[location] = x[mask, :].T
        features[location] = metadata.index[mask].astype(str).tolist()
        rows = metadata.loc[mask].copy()
        rows.insert(0, "glare_location", location)
        retained_rows.append(rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / "controlled_target.npz",
        flt=matrices["FLT"],
        gc=matrices["GC"],
        genes=np.asarray(genes, dtype=str),
        flt_features=np.asarray(features["FLT"], dtype=str),
        gc_features=np.asarray(features["GC"], dtype=str),
        input_kind=np.asarray("mober_projected_log2_cpm_expression"),
        input_path=np.asarray(str(mober_h5ad)),
    )
    retained = pd.concat(retained_rows, ignore_index=True)
    retained.to_csv(output_dir / "retained_profile_features.tsv", sep="\t", index=False)
    counts = (
        retained.groupby(["data_source", "glare_location"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .rename(columns={"data_source": "h5_accession"})
    )
    counts.to_csv(output_dir / "mober_glare_condition_counts.tsv", sep="\t", index=False)
    return {
        "mober_h5ad": str(mober_h5ad),
        "genes": genes,
        "matrices": matrices,
        "features": features,
        "counts": counts,
        "shape_samples_x_genes": list(adata.shape),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fine-tune GLARE on MOBER-projected aggregate liver expression."
    )
    parser.add_argument("--mober-h5ad", default=DEFAULT_MOBER_H5AD)
    parser.add_argument("--pretrained-weights", default=DEFAULT_PRETRAINED_WEIGHTS)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1996)
    parser.add_argument("--prepare-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_start = time.perf_counter()
    output_dir = Path(args.output_dir)
    prepared = prepare_mober_target(args.mober_h5ad, output_dir)
    log(
        "Prepared MOBER-corrected GLARE target: "
        f"{len(prepared['genes'])} genes, "
        f"{prepared['matrices']['FLT'].shape[1]} FLT and "
        f"{prepared['matrices']['GC'].shape[1]} GC samples"
    )
    if args.prepare_only:
        return

    pretrained_weights = Path(args.pretrained_weights)
    input_dim = infer_pretrain_input_dim(pretrained_weights)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    for location, matrix in prepared["matrices"].items():
        write_outlier_audit(matrix, prepared["genes"], location, output_dir)

    locations = [
        finetune_location(
            prepared["matrices"][location],
            prepared["genes"],
            location,
            pretrained_weights,
            output_dir,
            device,
            input_dim,
            args.epochs,
            args.batch_size,
            args.seed,
        )
        for location in ("FLT", "GC")
    ]
    summary = {
        "method": "GLARE released 16-dimensional SAE fine-tuned on MOBER-projected aggregate liver FLT/GC expression",
        "mober_h5ad": prepared["mober_h5ad"],
        "mober_projection_target": "OSD-379",
        "target_expression_kind": "mober_projected_log2_cpm_expression",
        "shape_samples_x_genes": prepared["shape_samples_x_genes"],
        "condition_counts": prepared["counts"].to_dict(orient="records"),
        "pretrained_weights": str(pretrained_weights),
        "pretrained_input_dim": input_dim,
        "device": str(device),
        "seed_reused_for_each_location": args.seed,
        "architecture": [128, 64, 32, 16],
        "learning_rate": 1e-3,
        "weight_decay": 0,
        "sparsity_penalty": 1e-5,
        "batch_size": args.batch_size,
        "locations": locations,
        "elapsed_seconds": round(time.perf_counter() - run_start, 3),
        "elapsed": format_elapsed(time.perf_counter() - run_start),
    }
    (output_dir / "finetune_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_report(output_dir, summary)
    log(f"Saved summary: {output_dir / 'finetune_summary.json'}")


def write_report(output_dir: Path, summary: dict) -> None:
    counts = pd.DataFrame(summary["condition_counts"])
    count_lines = counts.to_csv(sep="\t", index=False).strip().splitlines()
    location_rows = []
    for location in summary["locations"]:
        location_rows.append(
            {
                "location": location["location"],
                "genes": location["genes"],
                "profiles": location["profiles"],
                "best_loss": location["best_loss"],
                "best_epoch": location["best_epoch"],
                "epochs": location["epochs"],
            }
        )
    location_lines = pd.DataFrame(location_rows).to_csv(
        sep="\t", index=False
    ).strip().splitlines()
    text = [
        "# GLARE on MOBER-Corrected Aggregate Liver Data",
        "",
        f"- MOBER input: `{summary['mober_h5ad']}`",
        f"- Projection target: `{summary['mober_projection_target']}`",
        f"- Shape: {summary['shape_samples_x_genes'][0]:,} samples x "
        f"{summary['shape_samples_x_genes'][1]:,} genes",
        "",
        "## Condition Counts",
        "",
        "```tsv",
        *count_lines,
        "```",
        "",
        "## Fine-Tuning",
        "",
        "```tsv",
        *location_lines,
        "```",
    ]
    (output_dir / "RUN_SUMMARY.md").write_text(
        "\n".join(text) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
