"""Shared-split FLT/GC PCA figures for calibrated DDIM and WGAN samples."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from nasa_mouse_rna_diffusion.factorized_adapter import load_factorized_role
from nasa_mouse_rna_diffusion.factorized_calibrate import _aligned_expression


COLORS = {"flight": "#C14924", "ground_control": "#176B87"}
MODEL_LABELS = {
    "real": "Real OSDR",
    "diffusion": "Calibrated DDIM",
    "wgan": "Calibrated WGAN-GP",
}


def _center_within_accession(
    expression: np.ndarray, accessions: pd.Series
) -> np.ndarray:
    centered = np.asarray(expression, dtype=np.float32).copy()
    values = accessions.astype(str).to_numpy()
    for accession in np.unique(values):
        mask = values == accession
        centered[mask] -= centered[mask].mean(axis=0, keepdims=True)
    return centered


def _eligible_accessions(
    samples: pd.DataFrame, minimum_per_condition: int = 2
) -> list[str]:
    counts = (
        samples.groupby(["accession", "condition"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    for condition in COLORS:
        if condition not in counts:
            counts[condition] = 0
    eligible = counts.index[
        (counts["flight"] >= int(minimum_per_condition))
        & (counts["ground_control"] >= int(minimum_per_condition))
    ]
    return sorted(map(str, eligible))


def _condition_scatter(
    axis: plt.Axes,
    coordinates: np.ndarray,
    conditions: pd.Series,
    *,
    first: int = 0,
    second: int = 1,
    show_centroids: bool = True,
) -> None:
    condition_values = conditions.astype(str).to_numpy()
    for condition, color in COLORS.items():
        mask = condition_values == condition
        axis.scatter(
            coordinates[mask, first],
            coordinates[mask, second],
            s=16,
            alpha=0.58,
            color=color,
            edgecolors="none",
            rasterized=True,
            label=condition.replace("_", " "),
        )
        if show_centroids and mask.any():
            centroid = coordinates[mask][:, [first, second]].mean(axis=0)
            axis.scatter(
                centroid[0],
                centroid[1],
                s=78,
                marker="P",
                color=color,
                edgecolors="white",
                linewidths=0.8,
                zorder=4,
            )


def _axis_limits(values: list[np.ndarray], column: int) -> tuple[float, float]:
    observed = np.concatenate([value[:, column] for value in values])
    lower, upper = np.quantile(observed, [0.005, 0.995])
    padding = max(float(upper - lower) * 0.08, 1e-3)
    return float(lower - padding), float(upper + padding)


def _plot_global_three_pc(
    matrices: Mapping[str, np.ndarray], samples: pd.DataFrame, output: Path
) -> dict[str, object]:
    pca = PCA(n_components=3, random_state=0).fit(matrices["real"])
    coordinates = {name: pca.transform(values) for name, values in matrices.items()}
    limits = [_axis_limits(list(coordinates.values()), index) for index in range(3)]

    figure = plt.figure(figsize=(16.2, 5.5))
    for panel, name in enumerate(("real", "diffusion", "wgan"), start=1):
        axis = figure.add_subplot(1, 3, panel, projection="3d")
        condition_values = samples["condition"].astype(str).to_numpy()
        for condition, color in COLORS.items():
            mask = condition_values == condition
            points = coordinates[name][mask]
            axis.scatter(
                points[:, 0],
                points[:, 1],
                points[:, 2],
                s=12,
                alpha=0.55,
                color=color,
                edgecolors="none",
                rasterized=True,
                label=condition.replace("_", " "),
            )
        axis.set_title(MODEL_LABELS[name], fontweight="bold")
        axis.set_xlabel("PC1")
        axis.set_ylabel("PC2")
        axis.set_zlabel("PC3")
        axis.set_xlim(*limits[0])
        axis.set_ylim(*limits[1])
        axis.set_zlim(*limits[2])
        axis.view_init(elev=22, azim=38)
    handles, labels = figure.axes[-1].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    figure.suptitle(
        "FLT and GC in a shared real-fitted three-PC space",
        fontsize=15,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0.07, 1, 0.95))
    path = output / "global_condition_pca_3d.png"
    figure.savefig(path, dpi=240, bbox_inches="tight")
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)

    pairs = ((0, 1), (0, 2), (1, 2))
    figure, axes = plt.subplots(3, 3, figsize=(14.5, 13.2))
    for row, name in enumerate(("real", "diffusion", "wgan")):
        for column, (first, second) in enumerate(pairs):
            axis = axes[row, column]
            _condition_scatter(
                axis,
                coordinates[name],
                samples["condition"],
                first=first,
                second=second,
            )
            axis.set_xlim(*limits[first])
            axis.set_ylim(*limits[second])
            axis.set_xlabel(f"PC{first + 1}")
            axis.set_ylabel(f"PC{second + 1}")
            axis.grid(alpha=0.15, linewidth=0.7)
            if column == 0:
                axis.set_title(MODEL_LABELS[name], loc="left", fontweight="bold")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    figure.suptitle(
        "Pairwise views of the first three real-fitted PCs",
        fontsize=15,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0.04, 1, 0.97))
    pairwise_path = output / "global_condition_pca_pairwise.png"
    figure.savefig(pairwise_path, dpi=220, bbox_inches="tight")
    figure.savefig(pairwise_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return {
        "three_dimensional": str(path),
        "pairwise": str(pairwise_path),
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
    }


def _metric_text(table: pd.DataFrame, tissue: str, *, accession: bool) -> str:
    row = table.loc[table["tissue"].astype(str).eq(tissue)]
    if row.empty:
        return "unavailable"
    row = row.iloc[0]
    correlation_column = (
        "accession_meta_effect_correlation" if accession else "delta_correlation"
    )
    direction_column = (
        "accession_meta_direction_agreement"
        if accession
        else "direction_agreement"
    )
    correlation = row.get(correlation_column, np.nan)
    direction = row.get(direction_column, np.nan)
    if not np.isfinite(correlation) or not np.isfinite(direction):
        return "insufficient data"
    return f"r={correlation:.2f}, direction={direction:.2f}"


def _plot_tissue(
    tissue: str,
    matrices: Mapping[str, np.ndarray],
    samples: pd.DataFrame,
    diffusion_metrics: pd.DataFrame,
    wgan_metrics: pd.DataFrame,
    output: Path,
    *,
    accession_centered: bool,
) -> Path | None:
    tissue_mask = samples["tissue"].astype(str).eq(tissue).to_numpy()
    tissue_samples = samples.loc[tissue_mask].reset_index(drop=True)
    if accession_centered:
        eligible = _eligible_accessions(tissue_samples)
        if len(eligible) < 2:
            return None
        eligible_mask = tissue_samples["accession"].astype(str).isin(eligible).to_numpy()
        original_positions = np.flatnonzero(tissue_mask)[eligible_mask]
        tissue_samples = tissue_samples.loc[eligible_mask].reset_index(drop=True)
    else:
        eligible = []
        original_positions = np.flatnonzero(tissue_mask)
    if len(tissue_samples) < 3:
        return None

    selected = {name: values[original_positions] for name, values in matrices.items()}
    if accession_centered:
        selected = {
            name: _center_within_accession(values, tissue_samples["accession"])
            for name, values in selected.items()
        }
    pca = PCA(n_components=2, random_state=0).fit(selected["real"])
    coordinates = {name: pca.transform(values) for name, values in selected.items()}
    x_limits = _axis_limits(list(coordinates.values()), 0)
    y_limits = _axis_limits(list(coordinates.values()), 1)
    figure, axes = plt.subplots(1, 3, figsize=(13.8, 4.4), sharex=True, sharey=True)
    for axis, name in zip(axes, ("real", "diffusion", "wgan")):
        _condition_scatter(axis, coordinates[name], tissue_samples["condition"])
        if accession_centered:
            accession_values = tissue_samples["accession"].astype(str).to_numpy()
            condition_values = tissue_samples["condition"].astype(str).to_numpy()
            for accession in eligible:
                group = accession_values == accession
                ground = coordinates[name][
                    group & (condition_values == "ground_control")
                ].mean(axis=0)
                flight = coordinates[name][group & (condition_values == "flight")].mean(
                    axis=0
                )
                axis.annotate(
                    "",
                    xy=flight,
                    xytext=ground,
                    arrowprops={"arrowstyle": "->", "color": "#222222", "alpha": 0.5},
                    zorder=3,
                )
        axis.set_title(MODEL_LABELS[name], fontweight="bold")
        axis.set_xlim(*x_limits)
        axis.set_ylim(*y_limits)
        axis.set_xlabel("PC1")
        axis.grid(alpha=0.15, linewidth=0.7)
    axes[0].set_ylabel("PC2")
    handles, labels = axes[-1].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    mode = "accession-centered" if accession_centered else "raw within-tissue"
    metric_mode = accession_centered
    figure.suptitle(
        f"{tissue.replace('_', ' ')}: {mode} FLT/GC PCA\n"
        f"DDIM {_metric_text(diffusion_metrics, tissue, accession=metric_mode)}; "
        f"WGAN {_metric_text(wgan_metrics, tissue, accession=metric_mode)}; "
        f"n={len(tissue_samples)}"
        + (f", eligible accessions={len(eligible)}" if accession_centered else ""),
        fontsize=12,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0.08, 1, 0.87))
    path = output / f"{tissue}.png"
    figure.savefig(path, dpi=220, bbox_inches="tight")
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return path


def create_condition_figures(
    *,
    prepared_h5: str | Path,
    samples_tsv: str | Path,
    diffusion_expression: str | Path,
    wgan_expression: str | Path,
    diffusion_tissue_metrics: str | Path,
    wgan_tissue_metrics: str | Path,
    output_dir: str | Path,
) -> Path:
    role = load_factorized_role(prepared_h5, samples_tsv, "validation")
    diffusion_payload = np.load(diffusion_expression)
    diffusion_rows = np.asarray(diffusion_payload["source_row"], dtype=np.int64)
    real, samples = _aligned_expression(role, diffusion_rows)
    diffusion = np.asarray(diffusion_payload["scaled_expression"], dtype=np.float32)
    wgan_payload = np.load(wgan_expression)
    wgan_role_order = np.asarray(wgan_payload["scaled_expression"], dtype=np.float32)
    if len(wgan_role_order) != len(role["expression"]):
        raise ValueError("WGAN samples do not match the validation role length")
    role_lookup = {
        int(source): index for index, source in enumerate(role["source_row"])
    }
    wgan = wgan_role_order[[role_lookup[int(source)] for source in diffusion_rows]]
    if real.shape != diffusion.shape or real.shape != wgan.shape:
        raise ValueError("Real, DDIM, and WGAN matrices must have identical shapes")
    genes = np.asarray(role["genes"]).astype(str)
    if not np.array_equal(np.asarray(diffusion_payload["genes"]).astype(str), genes):
        raise ValueError("Diffusion genes do not match prepared validation genes")

    matrices = {"real": real, "diffusion": diffusion, "wgan": wgan}
    diffusion_metrics = pd.read_csv(diffusion_tissue_metrics, sep="\t")
    wgan_metrics = pd.read_csv(wgan_tissue_metrics, sep="\t")
    output = Path(output_dir)
    raw_output = output / "per_tissue" / "raw"
    centered_output = output / "per_tissue" / "accession_centered"
    raw_output.mkdir(parents=True, exist_ok=True)
    centered_output.mkdir(parents=True, exist_ok=True)

    global_outputs = _plot_global_three_pc(matrices, samples, output)
    eligibility_rows: list[dict[str, object]] = []
    raw_paths: dict[str, str] = {}
    centered_paths: dict[str, str] = {}
    for tissue in sorted(samples["tissue"].astype(str).unique()):
        tissue_samples = samples.loc[samples["tissue"].astype(str).eq(tissue)]
        eligible = _eligible_accessions(tissue_samples)
        eligibility_rows.append(
            {
                "tissue": tissue,
                "profiles": int(len(tissue_samples)),
                "accessions": int(tissue_samples["accession"].nunique()),
                "eligible_accessions": int(len(eligible)),
                "eligible_accession_ids": ",".join(eligible),
            }
        )
        raw_path = _plot_tissue(
            tissue,
            matrices,
            samples,
            diffusion_metrics,
            wgan_metrics,
            raw_output,
            accession_centered=False,
        )
        if raw_path is not None:
            raw_paths[tissue] = str(raw_path)
        centered_path = _plot_tissue(
            tissue,
            matrices,
            samples,
            diffusion_metrics,
            wgan_metrics,
            centered_output,
            accession_centered=True,
        )
        if centered_path is not None:
            centered_paths[tissue] = str(centered_path)

    eligibility = pd.DataFrame(eligibility_rows)
    eligibility.to_csv(output / "accession_eligibility.tsv", sep="\t", index=False)
    summary = {
        "status": "complete",
        "split": "validation",
        "profiles": int(len(samples)),
        "genes": int(len(genes)),
        "conditions": samples["condition"].value_counts().to_dict(),
        "global": global_outputs,
        "raw_per_tissue_figures": raw_paths,
        "accession_centered_per_tissue_figures": centered_paths,
        "accession_centering": (
            "Each real or synthetic matrix is centered independently within accession. "
            "Only accessions with at least two FLT and two GC profiles are shown; "
            "arrows connect the GC centroid to the FLT centroid within each accession."
        ),
        "interpretation": (
            "PCA is descriptive. FLT/GC effect correlations and direction agreement "
            "remain the inferential recovery diagnostics."
        ),
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-h5", required=True)
    parser.add_argument("--samples-tsv", required=True)
    parser.add_argument("--diffusion-expression", required=True)
    parser.add_argument("--wgan-expression", required=True)
    parser.add_argument("--diffusion-tissue-metrics", required=True)
    parser.add_argument("--wgan-tissue-metrics", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    create_condition_figures(
        prepared_h5=args.prepared_h5,
        samples_tsv=args.samples_tsv,
        diffusion_expression=args.diffusion_expression,
        wgan_expression=args.wgan_expression,
        diffusion_tissue_metrics=args.diffusion_tissue_metrics,
        wgan_tissue_metrics=args.wgan_tissue_metrics,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
