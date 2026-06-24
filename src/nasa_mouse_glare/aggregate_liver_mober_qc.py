"""QC metrics for MOBER-corrected aggregate liver data."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


DEFAULT_RUN_DIR = "outputs/mober_liver_ribo6_osdr"
DEFAULT_ONTO = "OSD-379"


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def as_array(matrix) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return np.asarray(matrix)


def pca_table(name: str, x: np.ndarray, obs: pd.DataFrame, output_dir: Path) -> dict:
    n_components = min(20, x.shape[0] - 1, x.shape[1])
    scaled = StandardScaler().fit_transform(x)
    pca = PCA(n_components=n_components, random_state=1996)
    coords = pca.fit_transform(scaled)
    columns = [f"PC{i + 1}" for i in range(coords.shape[1])]
    table = pd.DataFrame(coords, index=obs.index, columns=columns)
    table.insert(0, "sample", table.index)
    for column in ["data_source", "h5_accession", "project_identifier", "location", "condition", "sex", "strain"]:
        if column in obs.columns:
            table.insert(1, column, obs[column].astype(str).to_numpy())
    table.to_csv(output_dir / f"{name}_pca.tsv", sep="\t", index=False)
    return {
        "name": name,
        "n_components": int(n_components),
        "pc1_variance_ratio": float(pca.explained_variance_ratio_[0]),
        "pc2_variance_ratio": float(pca.explained_variance_ratio_[1]) if n_components > 1 else np.nan,
        "coords": coords,
    }


def silhouette(coords: np.ndarray, labels: pd.Series) -> float:
    labels = labels.astype(str)
    if labels.nunique() < 2:
        return np.nan
    if labels.value_counts().min() < 2:
        return np.nan
    return float(silhouette_score(coords, labels))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="QC MOBER aggregate liver outputs.")
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR)
    parser.add_argument("--onto", default=DEFAULT_ONTO)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    projection_dir = run_dir / "projection"
    output_dir = run_dir / "mober_qc"
    output_dir.mkdir(parents=True, exist_ok=True)

    import anndata as ad

    log("Loading MOBER input and projected data")
    raw = ad.read_h5ad(run_dir / "mober_liver_ribo6_input.h5ad")
    projected = ad.read_h5ad(projection_dir / f"mober_projected_onto_{args.onto}.h5ad")
    latent = ad.read_h5ad(projection_dir / f"mober_latent_onto_{args.onto}.h5ad")

    raw_info = pca_table("input_log2_cpm", as_array(raw.X), raw.obs, output_dir)
    projected_info = pca_table(
        f"mober_projected_onto_{args.onto}",
        as_array(projected.X),
        projected.obs,
        output_dir,
    )
    latent_x = as_array(latent.X)
    latent_table = pd.DataFrame(
        latent_x,
        index=latent.obs_names,
        columns=[f"z_{i}" for i in range(latent_x.shape[1])],
    )
    latent_table.insert(0, "sample", latent_table.index)
    for column in ["data_source", "h5_accession", "project_identifier", "location", "condition", "sex", "strain"]:
        if column in latent.obs.columns:
            latent_table.insert(1, column, latent.obs[column].astype(str).to_numpy())
    latent_table.to_csv(output_dir / f"mober_latent_onto_{args.onto}.tsv", sep="\t", index=False)

    metric_rows = []
    for info, obs in [(raw_info, raw.obs), (projected_info, projected.obs)]:
        coords = info["coords"][:, : min(10, info["coords"].shape[1])]
        metric_rows.append(
            {
                "space": info["name"],
                "pc1_variance_ratio": info["pc1_variance_ratio"],
                "pc2_variance_ratio": info["pc2_variance_ratio"],
                "data_source_silhouette_pc1_10": silhouette(coords, obs["data_source"]),
                "condition_silhouette_pc1_10": silhouette(coords, obs["condition"]),
                "sex_silhouette_pc1_10": silhouette(coords, obs["sex"]),
                "strain_silhouette_pc1_10": silhouette(coords, obs["strain"]),
            }
        )
    metric_rows.append(
        {
            "space": f"mober_latent_onto_{args.onto}",
            "pc1_variance_ratio": np.nan,
            "pc2_variance_ratio": np.nan,
            "data_source_silhouette_pc1_10": silhouette(latent_x, latent.obs["data_source"]),
            "condition_silhouette_pc1_10": silhouette(latent_x, latent.obs["condition"]),
            "sex_silhouette_pc1_10": silhouette(latent_x, latent.obs["sex"]),
            "strain_silhouette_pc1_10": silhouette(latent_x, latent.obs["strain"]),
        }
    )
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(output_dir / "mober_qc_metrics.tsv", sep="\t", index=False)
    summary = {
        "run_dir": str(run_dir),
        "onto": args.onto,
        "metrics": metrics.replace({np.nan: None}).to_dict(orient="records"),
        "notes": [
            "Lower data_source silhouette after projection is consistent with reduced batch/source separation.",
            "Condition, sex, and strain are partly confounded with accession in this design.",
        ],
    }
    (output_dir / "mober_qc_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_report(output_dir, metrics, args.onto)
    log(f"Saved MOBER QC to {output_dir}")


def write_report(output_dir: Path, metrics: pd.DataFrame, onto: str) -> None:
    metric_lines = metrics.to_csv(sep="\t", index=False).strip().splitlines()
    text = [
        "# MOBER QC Summary",
        "",
        f"MOBER projection target: `{onto}`.",
        "",
        "Silhouette scores use PCA dimensions 1-10 for expression spaces and all",
        "64 dimensions for the MOBER latent space. Lower `data_source` silhouette",
        "after projection is consistent with less accession/batch separation.",
        "",
        "```tsv",
        *metric_lines,
        "```",
    ]
    (output_dir / "MOBER_QC_SUMMARY.md").write_text(
        "\n".join(text) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
