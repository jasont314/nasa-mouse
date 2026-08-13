"""Plot sample-level study separation from GLARE cluster/module scores."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcache")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

try:
    from umap import UMAP
except Exception:  # pragma: no cover - optional plotting dependency
    UMAP = None


SCOPES = ("aggregate", "aggregate_mober")
CONDITION_MARKERS = {"FLT": "^", "GC": "o", "flight": "^", "ground": "o"}


def clean_label(value: object) -> str:
    text = str(value) if value is not None else ""
    text = text.strip()
    return text if text and text.lower() != "nan" else "unknown"


def tissue_label(tissue_slug: str) -> str:
    return tissue_slug.replace("_", " ")


def read_vector(path: Path, column_name: str) -> pd.Series:
    return pd.read_csv(path, sep="\t", header=None, names=[column_name])[column_name].astype(str)


def load_scope(scope_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    controlled_path = scope_dir / "controlled_target.npz"
    metadata_path = scope_dir / "inputs" / "aligned_tms_api.target.profile_metadata.tsv"
    if not metadata_path.exists() and scope_dir.name == "aggregate_mober":
        metadata_path = (
            scope_dir.parent
            / "aggregate"
            / "inputs"
            / "aligned_tms_api.target.profile_metadata.tsv"
        )
    if not controlled_path.exists():
        raise FileNotFoundError(controlled_path)
    if not metadata_path.exists():
        raise FileNotFoundError(metadata_path)

    controlled = np.load(controlled_path)
    genes = controlled["genes"].astype(str)
    flt = controlled["flt"].astype(float)
    gc = controlled["gc"].astype(float)
    flt_features = controlled["flt_features"].astype(str)
    gc_features = controlled["gc_features"].astype(str)

    expression = np.concatenate([flt, gc], axis=1)
    sample_ids = np.concatenate([flt_features, gc_features])
    expression_df = pd.DataFrame(expression, index=genes, columns=sample_ids)

    metadata = pd.read_csv(metadata_path, sep="\t")
    metadata["feature"] = metadata["feature"].astype(str)
    metadata = metadata.set_index("feature").reindex(sample_ids).reset_index()
    if metadata["id.accession"].isna().any():
        missing = metadata.loc[metadata["id.accession"].isna(), "feature"].head(5).tolist()
        raise ValueError(f"Missing metadata for {len(missing)} samples, examples: {missing}")

    metadata["condition_label"] = metadata.get("condition_label", "").map(clean_label)
    if metadata["condition_label"].eq("unknown").any() and "condition" in metadata:
        metadata["condition_label"] = metadata["condition"].map(clean_label)
    metadata["id.accession"] = metadata["id.accession"].map(clean_label)
    metadata["sample_id"] = sample_ids
    return expression_df, metadata


def cluster_modules(scope_dir: Path) -> pd.DataFrame:
    frames = []
    for condition in ("FLT", "GC"):
        path = scope_dir / "clustering" / f"{condition}_gene_clusters.tsv"
        if not path.exists():
            raise FileNotFoundError(path)
        clusters = pd.read_csv(path, sep="\t")
        if "gene_id" not in clusters or "consensus" not in clusters:
            raise ValueError(f"Expected gene_id and consensus columns in {path}")
        clusters = clusters[["gene_id", "consensus"]].copy()
        clusters["condition_module"] = condition
        clusters["module"] = [
            f"{condition}_C{int(cluster):02d}" for cluster in clusters["consensus"]
        ]
        frames.append(clusters)
    return pd.concat(frames, ignore_index=True)


def compute_module_scores(expression: pd.DataFrame, modules: pd.DataFrame) -> pd.DataFrame:
    values = expression.to_numpy(dtype=float)
    means = np.nanmean(values, axis=1, keepdims=True)
    stds = np.nanstd(values, axis=1, keepdims=True)
    stds[stds == 0] = 1.0
    z = (values - means) / stds
    z_expression = pd.DataFrame(z, index=expression.index, columns=expression.columns)

    score_columns: dict[str, np.ndarray] = {}
    for module, group in modules.groupby("module", sort=True):
        genes = [gene for gene in group["gene_id"].astype(str) if gene in z_expression.index]
        if len(genes) < 2:
            continue
        score_columns[module] = z_expression.loc[genes].mean(axis=0).to_numpy()
    if not score_columns:
        raise ValueError("No module scores could be computed")
    scores = pd.DataFrame(score_columns, index=expression.columns)
    scores.index.name = "sample_id"
    return scores


def safe_silhouette(matrix: np.ndarray, labels: pd.Series) -> float:
    labels = labels.astype(str).to_numpy()
    unique, counts = np.unique(labels, return_counts=True)
    if len(unique) < 2 or len(unique) >= len(labels) or np.any(counts < 2):
        return float("nan")
    try:
        return float(silhouette_score(matrix, labels))
    except Exception:
        return float("nan")


def reduce_scores(scores: pd.DataFrame, seed: int) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
    scaled = StandardScaler().fit_transform(scores.to_numpy(dtype=float))
    pca_model = PCA(n_components=2, random_state=seed)
    pca_coords = pca_model.fit_transform(scaled)
    umap_coords = None
    if UMAP is not None and len(scores) >= 5:
        neighbors = max(2, min(15, len(scores) - 1))
        umap_coords = UMAP(
            n_components=2,
            n_neighbors=neighbors,
            min_dist=0.3,
            metric="euclidean",
            random_state=seed,
        ).fit_transform(scaled)
    return pca_coords, umap_coords, pca_model.explained_variance_ratio_


def color_map(labels: list[str]) -> dict[str, tuple[float, float, float, float]]:
    palette = list(plt.get_cmap("tab20").colors)
    if len(labels) > len(palette):
        palette = [plt.get_cmap("hsv")(i / len(labels)) for i in range(len(labels))]
    return {label: palette[index % len(palette)] for index, label in enumerate(labels)}


def plot_single(
    coords: np.ndarray,
    metadata: pd.DataFrame,
    colors: dict[str, object],
    title: str,
    xlabel: str,
    ylabel: str,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 6.2), dpi=170)
    plot_axis(ax, coords, metadata, colors, title, xlabel, ylabel, show_legends=True)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_axis(
    ax: plt.Axes,
    coords: np.ndarray,
    metadata: pd.DataFrame,
    colors: dict[str, object],
    title: str,
    xlabel: str,
    ylabel: str,
    show_legends: bool,
) -> None:
    accessions = metadata["id.accession"].astype(str).tolist()
    conditions = metadata["condition_label"].astype(str).tolist()
    for accession in sorted(set(accessions)):
        for condition in sorted(set(conditions)):
            mask = (metadata["id.accession"].astype(str).eq(accession)) & (
                metadata["condition_label"].astype(str).eq(condition)
            )
            if not mask.any():
                continue
            marker = CONDITION_MARKERS.get(condition, "s")
            ax.scatter(
                coords[mask.to_numpy(), 0],
                coords[mask.to_numpy(), 1],
                s=42,
                marker=marker,
                color=colors[accession],
                edgecolor="white",
                linewidth=0.45,
                alpha=0.86,
                label=accession if condition == sorted(set(conditions))[0] else None,
            )
    ax.axhline(0, color="#d0d0d0", linewidth=0.6, zorder=0)
    ax.axvline(0, color="#d0d0d0", linewidth=0.6, zorder=0)
    ax.set_title(title, fontsize=10, pad=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, color="#eeeeee", linewidth=0.5)
    if show_legends:
        accession_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=colors[accession],
                markeredgecolor="white",
                label=accession,
                markersize=7,
            )
            for accession in sorted(set(accessions))
        ]
        condition_handles = [
            plt.Line2D(
                [0],
                [0],
                marker=CONDITION_MARKERS.get(condition, "s"),
                color="#4a4a4a",
                linestyle="none",
                label=condition,
                markersize=7,
            )
            for condition in sorted(set(conditions))
        ]
        leg1 = ax.legend(
            handles=accession_handles,
            title="OSDR accession",
            loc="center left",
            bbox_to_anchor=(1.02, 0.52),
            frameon=False,
            fontsize=7,
            title_fontsize=8,
        )
        ax.add_artist(leg1)
        ax.legend(
            handles=condition_handles,
            title="Condition",
            loc="center left",
            bbox_to_anchor=(1.02, 0.08),
            frameon=False,
            fontsize=7,
            title_fontsize=8,
        )


def plot_pair(
    direct: dict[str, object],
    mober: dict[str, object],
    coord_key: str,
    title: str,
    output_path: Path,
) -> None:
    metadata_all = pd.concat([direct["metadata"], mober["metadata"]], ignore_index=True)
    colors = color_map(sorted(metadata_all["id.accession"].astype(str).unique()))
    fig, axes = plt.subplots(1, 2, figsize=(14.8, 5.9), dpi=170, sharex=False, sharey=False)
    for ax, item, label in [
        (axes[0], direct, "Aggregate"),
        (axes[1], mober, "Aggregate + MOBER"),
    ]:
        coords = item[coord_key]
        if coords is None:
            ax.axis("off")
            ax.set_title(f"{label}: unavailable")
            continue
        metrics = item["metrics"]
        subtitle = (
            f"{label} | accession sil {metrics['accession_silhouette']:.3g}; "
            f"condition sil {metrics['condition_silhouette']:.3g}; n={metrics['samples']}"
        )
        if coord_key == "pca":
            xlab = f"PC1 ({item['pca_var'][0] * 100:.1f}%)"
            ylab = f"PC2 ({item['pca_var'][1] * 100:.1f}%)"
        else:
            xlab = "UMAP1"
            ylab = "UMAP2"
        plot_axis(ax, coords, item["metadata"], colors, subtitle, xlab, ylab, show_legends=False)

    accessions = sorted(metadata_all["id.accession"].astype(str).unique())
    accession_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=colors[accession],
            markeredgecolor="white",
            label=accession,
            markersize=7,
        )
        for accession in accessions
    ]
    conditions = sorted(metadata_all["condition_label"].astype(str).unique())
    condition_handles = [
        plt.Line2D(
            [0],
            [0],
            marker=CONDITION_MARKERS.get(condition, "s"),
            color="#4a4a4a",
            linestyle="none",
            label=condition,
            markersize=7,
        )
        for condition in conditions
    ]
    fig.suptitle(title, fontsize=12, y=0.98)
    fig.legend(
        handles=accession_handles,
        title="OSDR accession",
        loc="center left",
        bbox_to_anchor=(0.995, 0.56),
        frameon=False,
        fontsize=7,
        title_fontsize=8,
    )
    fig.legend(
        handles=condition_handles,
        title="Condition",
        loc="center left",
        bbox_to_anchor=(0.995, 0.18),
        frameon=False,
        fontsize=7,
        title_fontsize=8,
    )
    fig.tight_layout(rect=(0, 0, 0.86, 0.94))
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def process_scope(tissue: str, scope: str, scope_dir: Path, output_dir: Path, seed: int) -> dict[str, object]:
    expression, metadata = load_scope(scope_dir)
    modules = cluster_modules(scope_dir)
    scores = compute_module_scores(expression, modules)
    metadata = metadata.set_index("sample_id").reindex(scores.index).reset_index()
    pca_coords, umap_coords, pca_var = reduce_scores(scores, seed=seed)
    scaled = StandardScaler().fit_transform(scores.to_numpy(dtype=float))

    metrics = {
        "tissue": tissue,
        "scope": scope,
        "samples": int(scores.shape[0]),
        "modules": int(scores.shape[1]),
        "accessions": int(metadata["id.accession"].nunique()),
        "conditions": int(metadata["condition_label"].nunique()),
        "accession_silhouette": safe_silhouette(scaled, metadata["id.accession"]),
        "condition_silhouette": safe_silhouette(scaled, metadata["condition_label"]),
        "pca_accession_silhouette": safe_silhouette(pca_coords, metadata["id.accession"]),
        "pca_condition_silhouette": safe_silhouette(pca_coords, metadata["condition_label"]),
        "pc1_var": float(pca_var[0]),
        "pc2_var": float(pca_var[1]),
    }

    accession_colors = color_map(sorted(metadata["id.accession"].astype(str).unique()))
    prefix = f"{tissue}__{scope}"
    score_out = scores.reset_index().merge(
        metadata[["sample_id", "id.accession", "condition_label"]],
        on="sample_id",
        how="left",
    )
    score_out.to_csv(output_dir / f"{prefix}__module_scores.tsv", sep="\t", index=False)
    coord_df = metadata[["sample_id", "id.accession", "condition_label"]].copy()
    coord_df["pca1"] = pca_coords[:, 0]
    coord_df["pca2"] = pca_coords[:, 1]
    if umap_coords is not None:
        coord_df["umap1"] = umap_coords[:, 0]
        coord_df["umap2"] = umap_coords[:, 1]
    coord_df.to_csv(output_dir / f"{prefix}__sample_coordinates.tsv", sep="\t", index=False)

    plot_single(
        pca_coords,
        metadata,
        accession_colors,
        (
            f"{tissue_label(tissue)} {scope}: sample GLARE-module PCA by study\n"
            f"accession sil {metrics['accession_silhouette']:.3g}; "
            f"condition sil {metrics['condition_silhouette']:.3g}; n={metrics['samples']}"
        ),
        f"PC1 ({pca_var[0] * 100:.1f}%)",
        f"PC2 ({pca_var[1] * 100:.1f}%)",
        output_dir / f"{prefix}__pca_by_accession.png",
    )
    if umap_coords is not None:
        plot_single(
            umap_coords,
            metadata,
            accession_colors,
            (
                f"{tissue_label(tissue)} {scope}: sample GLARE-module UMAP by study\n"
                f"accession sil {metrics['accession_silhouette']:.3g}; "
                f"condition sil {metrics['condition_silhouette']:.3g}; n={metrics['samples']}"
            ),
            "UMAP1",
            "UMAP2",
            output_dir / f"{prefix}__umap_by_accession.png",
        )

    return {
        "metadata": metadata,
        "scores": scores,
        "pca": pca_coords,
        "umap": umap_coords,
        "pca_var": pca_var,
        "metrics": metrics,
    }


def discover_tissues(root: Path) -> list[str]:
    tissues = []
    for path in sorted(root.iterdir()):
        if not path.is_dir():
            continue
        if (path / "aggregate" / "controlled_target.npz").exists():
            tissues.append(path.name)
    return tissues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("outputs/glare/multi_tissue_api"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/glare/multi_tissue_api/study_effects"))
    parser.add_argument(
        "--report-source",
        type=Path,
        default=Path(
            "paper/slstp_internship_report/source_data/glare/"
            "skeletal_muscle_aggregate_vs_mober_umap_by_accession.png"
        ),
    )
    parser.add_argument("--results-dir", type=Path, default=Path("outputs/glare/study_effects"))
    parser.add_argument("--seed", type=int, default=1996)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_source.parent.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    pair_rows = []
    for tissue in discover_tissues(args.root):
        tissue_dir = args.root / tissue
        processed: dict[str, dict[str, object]] = {}
        for scope in SCOPES:
            scope_dir = tissue_dir / scope
            if not (scope_dir / "controlled_target.npz").exists():
                continue
            print(f"Processing {tissue} {scope}", flush=True)
            item = process_scope(tissue, scope, scope_dir, args.output_dir, args.seed)
            processed[scope] = item
            rows.append(item["metrics"])

        if all(scope in processed for scope in SCOPES):
            direct = processed["aggregate"]
            mober = processed["aggregate_mober"]
            direct_metrics = direct["metrics"]
            mober_metrics = mober["metrics"]
            pair_rows.append(
                {
                    "tissue": tissue,
                    "samples": direct_metrics["samples"],
                    "accessions": direct_metrics["accessions"],
                    "aggregate_accession_silhouette": direct_metrics["accession_silhouette"],
                    "mober_accession_silhouette": mober_metrics["accession_silhouette"],
                    "mober_minus_aggregate_accession_silhouette": (
                        mober_metrics["accession_silhouette"]
                        - direct_metrics["accession_silhouette"]
                    ),
                    "aggregate_condition_silhouette": direct_metrics["condition_silhouette"],
                    "mober_condition_silhouette": mober_metrics["condition_silhouette"],
                    "mober_minus_aggregate_condition_silhouette": (
                        mober_metrics["condition_silhouette"]
                        - direct_metrics["condition_silhouette"]
                    ),
                }
            )
            for coord_key in ("pca", "umap"):
                if direct[coord_key] is None or mober[coord_key] is None:
                    continue
                out_name = f"{tissue}__aggregate_vs_mober__{coord_key}_by_accession.png"
                plot_pair(
                    direct,
                    mober,
                    coord_key,
                    f"{tissue_label(tissue)}: sample GLARE-module {coord_key.upper()} colored by study",
                    args.output_dir / out_name,
                )
                if tissue == "skeletal_muscle" and coord_key == "umap":
                    plot_pair(
                        direct,
                        mober,
                        coord_key,
                        f"{tissue_label(tissue)}: sample GLARE-module {coord_key.upper()} colored by study",
                        args.report_source,
                    )

    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary = summary.sort_values(["tissue", "scope"])
    summary.to_csv(args.output_dir / "study_effect_summary.tsv", sep="\t", index=False)
    summary.to_csv(args.results_dir / "study_effect_summary.tsv", sep="\t", index=False)

    pair_summary = pd.DataFrame(pair_rows)
    if not pair_summary.empty:
        pair_summary = pair_summary.sort_values("tissue")
    pair_summary.to_csv(args.output_dir / "aggregate_vs_mober_study_effect_summary.tsv", sep="\t", index=False)
    pair_summary.to_csv(args.results_dir / "aggregate_vs_mober_study_effect_summary.tsv", sep="\t", index=False)

    readme = args.results_dir / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# GLARE Study-Effect Visualizations",
                "",
                "These plots score each OSDR sample by mean z-scored expression of genes in each",
                "GLARE consensus cluster, then reduce the sample-by-module score matrix with PCA",
                "or UMAP. Points are samples, colors are OSDR accessions, and marker shape is FLT/GC.",
                "",
                "The skeletal-muscle UMAP used by the internship report is frozen under",
                "`paper/slstp_internship_report/source_data/glare/`.",
                "Full per-scope coordinates and module scores are generated under ignored",
                "`outputs/glare/multi_tissue_api/study_effects/`.",
                "",
                "Lower accession silhouette after MOBER is consistent with reduced study/batch",
                "separation. Positive condition silhouette indicates stronger FLT/GC separation",
                "in the GLARE-module score space.",
                "",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
