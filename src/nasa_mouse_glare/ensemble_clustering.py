"""GLARE-style ensemble clustering for fine-tuned mouse GLARE representations."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from .io import require_import

os.environ.setdefault("LOKY_MAX_CPU_COUNT", str(os.cpu_count() or 1))


DEFAULT_REPRESENTATION = "outputs/glare_hpt_tms_facs_osdr/FTSAE_representation.npy"
DEFAULT_GENE_LATENT = "outputs/glare_hpt_tms_facs_osdr/post_finetune/gene_latent.tsv"
DEFAULT_GENE_PCA = "outputs/glare_hpt_tms_facs_osdr/post_finetune/gene_pca.tsv"
DEFAULT_OUTPUT_DIR = "outputs/glare_hpt_tms_facs_osdr/post_finetune/ensemble_clustering"


def scale_latent(latent, no_scale: bool):
    if no_scale:
        return latent
    StandardScaler = require_import(
        "sklearn.preprocessing", "pip install -r requirements-nasa-mouse-glare.txt"
    ).StandardScaler
    return StandardScaler().fit_transform(latent)


def run_gmm(
    latent,
    n_clusters: int,
    seed: int,
    max_iter: int,
    covariance_type: str,
    reg_covar: float,
    fallback_covariance_type: str,
):
    GaussianMixture = require_import(
        "sklearn.mixture", "pip install -r requirements-nasa-mouse-glare.txt"
    ).GaussianMixture
    model = GaussianMixture(
        n_components=n_clusters,
        covariance_type=covariance_type,
        reg_covar=reg_covar,
        random_state=seed,
        n_init=min(20, n_clusters),
        max_iter=max_iter,
    )
    try:
        return model.fit_predict(latent), model, covariance_type
    except ValueError:
        if fallback_covariance_type == covariance_type:
            raise
        fallback = GaussianMixture(
            n_components=n_clusters,
            covariance_type=fallback_covariance_type,
            reg_covar=reg_covar,
            random_state=seed,
            n_init=min(20, n_clusters),
            max_iter=max_iter,
        )
        return fallback.fit_predict(latent), fallback, fallback_covariance_type


def run_hdbscan(latent, min_cluster_size: int, min_samples: int, n_jobs: int):
    cluster_mod = require_import(
        "sklearn.cluster", "pip install 'scikit-learn>=1.3'"
    )
    HDBSCAN = getattr(cluster_mod, "HDBSCAN", None)
    if HDBSCAN is None:
        raise SystemExit(
            "sklearn.cluster.HDBSCAN is unavailable. Install scikit-learn>=1.3 "
            "or run this in the nasa conda env."
        )
    model = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_method="leaf",
        n_jobs=n_jobs,
        copy=False,
    )
    return model.fit_predict(latent), model


def run_spectral(latent, n_clusters: int, n_neighbors: int, seed: int, n_jobs: int):
    SpectralClustering = require_import(
        "sklearn.cluster", "pip install -r requirements-nasa-mouse-glare.txt"
    ).SpectralClustering
    model = SpectralClustering(
        n_clusters=n_clusters,
        affinity="nearest_neighbors",
        n_neighbors=n_neighbors,
        assign_labels="kmeans",
        random_state=seed,
        n_jobs=n_jobs,
    )
    return model.fit_predict(latent), model


def one_hot_membership(label_arrays: dict[str, object]):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    sparse = require_import("scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt")

    names = list(label_arrays)
    n_rows = len(next(iter(label_arrays.values())))
    row_parts = []
    col_parts = []
    data_parts = []
    feature_names = []
    offset = 0

    for name in names:
        labels = np.asarray(label_arrays[name])
        valid = labels >= 0
        unique_labels = sorted(set(labels[valid].tolist()))
        label_to_col = {label: idx + offset for idx, label in enumerate(unique_labels)}
        rows = np.flatnonzero(valid)
        cols = np.array([label_to_col[label] for label in labels[valid]], dtype="int64")
        row_parts.append(rows)
        col_parts.append(cols)
        data_parts.append(np.ones(len(rows), dtype="float32"))
        feature_names.extend([f"{name}:{label}" for label in unique_labels])
        offset += len(unique_labels)

    if offset == 0:
        raise SystemExit("No non-noise base cluster labels were available for consensus")

    matrix = sparse.csr_matrix(
        (
            np.concatenate(data_parts),
            (np.concatenate(row_parts), np.concatenate(col_parts)),
        ),
        shape=(n_rows, offset),
        dtype="float32",
    )
    return matrix, feature_names


def consensus_from_base_labels(label_arrays: dict[str, object], n_clusters: int, seed: int):
    MiniBatchKMeans = require_import(
        "sklearn.cluster", "pip install -r requirements-nasa-mouse-glare.txt"
    ).MiniBatchKMeans
    membership, feature_names = one_hot_membership(label_arrays)
    model = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=seed,
        n_init=20,
        batch_size=4096,
    )
    labels = model.fit_predict(membership)
    return labels, model, feature_names


def silhouette(latent, labels, sample_size: int, seed: int):
    silhouette_score = require_import(
        "sklearn.metrics", "pip install -r requirements-nasa-mouse-glare.txt"
    ).silhouette_score
    unique = sorted(set(labels.tolist()))
    if len(unique) < 2:
        return None
    kwargs = {"metric": "euclidean", "random_state": seed}
    if sample_size and len(labels) > sample_size:
        kwargs["sample_size"] = sample_size
    return float(silhouette_score(latent, labels, **kwargs))


def cluster_count(labels) -> int:
    return len({int(label) for label in labels.tolist() if int(label) >= 0})


def noise_count(labels) -> int:
    return int((labels < 0).sum())


def write_outputs(
    genes: list[str],
    latent,
    pca_df,
    labels: dict[str, object],
    consensus_labels,
    metrics_rows: list[dict[str, object]],
    feature_names: list[str],
    args,
    output_dir: Path,
):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    output_dir.mkdir(parents=True, exist_ok=True)

    cluster_df = pd.DataFrame({"gene": genes})
    for name, values in labels.items():
        cluster_df[f"{name}_cluster"] = values.astype(int)
    cluster_df["consensus_cluster"] = consensus_labels.astype(int)
    clusters_path = output_dir / "ensemble_clusters.tsv"
    cluster_df.to_csv(clusters_path, sep="\t", index=False)

    gene_clusters = cluster_df[["gene", "consensus_cluster"]].rename(
        columns={"consensus_cluster": "gene_cluster"}
    )
    gene_clusters_path = output_dir / "gene_clusters.tsv"
    gene_clusters.to_csv(gene_clusters_path, sep="\t", index=False)

    summary_rows = []
    for cluster in sorted(set(consensus_labels.tolist())):
        mask = consensus_labels == cluster
        row = {
            "consensus_cluster": int(cluster),
            "n_genes": int(mask.sum()),
            "latent_centroid_norm": float(np.linalg.norm(latent[mask].mean(axis=0))),
        }
        if pca_df is not None:
            pca_masked = pca_df.loc[mask]
            for col in ["pc1", "pc2"]:
                if col in pca_masked:
                    row[f"{col}_mean"] = float(pca_masked[col].mean())
        summary_rows.append(row)
    summary_path = output_dir / "ensemble_cluster_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

    base_rows = []
    for name, values in labels.items():
        for cluster, count in pd.Series(values).value_counts().sort_index().items():
            base_rows.append(
                {
                    "algorithm": name,
                    "cluster": int(cluster),
                    "n_genes": int(count),
                    "is_noise": bool(int(cluster) < 0),
                }
            )
    base_summary_path = output_dir / "ensemble_base_cluster_summary.tsv"
    pd.DataFrame(base_rows).to_csv(base_summary_path, sep="\t", index=False)

    metrics_path = output_dir / "ensemble_metrics.tsv"
    pd.DataFrame(metrics_rows).to_csv(metrics_path, sep="\t", index=False)

    feature_path = output_dir / "ensemble_consensus_features.txt"
    feature_path.write_text("\n".join(feature_names) + "\n", encoding="utf-8")

    plot_path = None
    if pca_df is not None and not args.no_plots:
        plot_path = output_dir / "ensemble_pca_by_consensus.png"
        write_pca_plot(pca_df, consensus_labels, plot_path)

    summary = {
        "representation": str(args.representation),
        "gene_latent": str(args.gene_latent),
        "gene_pca": str(args.gene_pca),
        "n_genes": len(genes),
        "n_clusters": args.n_clusters,
        "seed": args.seed,
        "scaled_latent": not args.no_scale,
        "base_algorithms": {
            "gmm": {
                "n_components": args.n_clusters,
                "max_iter": args.gmm_max_iter,
                "covariance_type_requested": args.gmm_covariance_type,
                "reg_covar": args.gmm_reg_covar,
                "fallback_covariance_type": args.gmm_fallback_covariance_type,
                "covariance_type_used": args.gmm_covariance_type_used,
            },
            "hdbscan": {
                "min_cluster_size": args.hdbscan_min_cluster_size,
                "min_samples": args.hdbscan_min_samples,
                "cluster_selection_method": "leaf",
                "implementation": "sklearn.cluster.HDBSCAN",
            },
            "spectral": {
                "n_clusters": args.n_clusters,
                "affinity": "nearest_neighbors",
                "n_neighbors": args.spectral_neighbors,
            },
        },
        "consensus": {
            "method": (
                "MiniBatchKMeans on sparse one-hot base-cluster membership; "
                "equivalent feature embedding for EAC co-association labels"
            ),
            "n_clusters": args.n_clusters,
            "hdbscan_noise_handling": (
                "noise points are omitted from the HDBSCAN one-hot block so they "
                "do not form one artificial noise cluster"
            ),
        },
        "outputs": {
            "ensemble_clusters": str(clusters_path),
            "gene_clusters": str(gene_clusters_path),
            "ensemble_cluster_summary": str(summary_path),
            "ensemble_base_cluster_summary": str(base_summary_path),
            "ensemble_metrics": str(metrics_path),
            "ensemble_consensus_features": str(feature_path),
            "ensemble_pca_by_consensus": str(plot_path) if plot_path else None,
        },
    }
    summary_path_json = output_dir / "ensemble_summary.json"
    summary_path_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary_path_json


def write_pca_plot(pca_df, consensus_labels, path: Path):
    mpl_config_dir = path.parent / ".matplotlib"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(mpl_config_dir))
    matplotlib = require_import("matplotlib", "pip install matplotlib")
    matplotlib.use("Agg")
    plt = require_import("matplotlib.pyplot", "pip install matplotlib")

    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    cmap = plt.get_cmap("tab20")
    for idx, cluster in enumerate(sorted(set(consensus_labels.tolist()))):
        mask = consensus_labels == cluster
        ax.scatter(
            pca_df.loc[mask, "pc1"],
            pca_df.loc[mask, "pc2"],
            s=5,
            alpha=0.65,
            color=cmap(idx % 20),
            label=str(cluster),
        )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("GLARE-style ensemble consensus clusters")
    if len(set(consensus_labels.tolist())) <= 30:
        ax.legend(title="consensus", markerscale=2, fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def run(args) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    adjusted_mutual_info_score = require_import(
        "sklearn.metrics", "pip install -r requirements-nasa-mouse-glare.txt"
    ).adjusted_mutual_info_score

    latent = np.load(args.representation)
    latent = scale_latent(latent, args.no_scale)

    gene_latent = pd.read_csv(args.gene_latent, sep="\t")
    if len(gene_latent) != latent.shape[0]:
        raise SystemExit(
            f"gene_latent rows do not match representation rows: "
            f"{len(gene_latent)} vs {latent.shape[0]}"
        )
    genes = gene_latent["gene"].astype(str).tolist()

    pca_df = None
    gene_pca_path = Path(args.gene_pca)
    if gene_pca_path.exists():
        pca_df = pd.read_csv(gene_pca_path, sep="\t")
        if len(pca_df) != latent.shape[0]:
            pca_df = None

    print("running_gmm", flush=True)
    gmm_labels, _, gmm_covariance_type_used = run_gmm(
        latent,
        args.n_clusters,
        args.seed,
        args.gmm_max_iter,
        args.gmm_covariance_type,
        args.gmm_reg_covar,
        args.gmm_fallback_covariance_type,
    )
    args.gmm_covariance_type_used = gmm_covariance_type_used
    print("running_hdbscan", flush=True)
    hdbscan_labels, _ = run_hdbscan(
        latent,
        args.hdbscan_min_cluster_size,
        args.hdbscan_min_samples,
        args.n_jobs,
    )
    print("running_spectral", flush=True)
    spectral_labels, _ = run_spectral(
        latent,
        args.n_clusters,
        args.spectral_neighbors,
        args.seed,
        args.n_jobs,
    )

    base_labels = {
        "gmm": np.asarray(gmm_labels, dtype=int),
        "hdbscan": np.asarray(hdbscan_labels, dtype=int),
        "spectral": np.asarray(spectral_labels, dtype=int),
    }

    print("running_consensus", flush=True)
    consensus_labels, _, feature_names = consensus_from_base_labels(
        base_labels,
        args.n_clusters,
        args.seed,
    )
    consensus_labels = np.asarray(consensus_labels, dtype=int)

    metrics_rows = []
    for name, values in [*base_labels.items(), ("consensus", consensus_labels)]:
        values = np.asarray(values, dtype=int)
        non_noise = values >= 0
        sil = None
        if non_noise.sum() > 1 and len(set(values[non_noise].tolist())) > 1:
            sil = silhouette(latent[non_noise], values[non_noise], args.silhouette_sample_size, args.seed)
        metrics_rows.append(
            {
                "algorithm": name,
                "n_clusters": cluster_count(values),
                "noise_genes": noise_count(values),
                "silhouette": sil,
                "adjusted_mutual_info_vs_consensus": (
                    1.0
                    if name == "consensus"
                    else float(adjusted_mutual_info_score(consensus_labels, values))
                ),
            }
        )

    output_dir = Path(args.output_dir)
    summary_path = write_outputs(
        genes,
        latent,
        pca_df,
        base_labels,
        consensus_labels,
        metrics_rows,
        feature_names,
        args,
        output_dir,
    )

    print(f"summary={summary_path}", flush=True)
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run GLARE-style GMM/HDBSCAN/Spectral ensemble clustering."
    )
    parser.add_argument("--representation", default=DEFAULT_REPRESENTATION)
    parser.add_argument("--gene-latent", default=DEFAULT_GENE_LATENT)
    parser.add_argument("--gene-pca", default=DEFAULT_GENE_PCA)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-clusters", type=int, default=15)
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--gmm-max-iter", type=int, default=100)
    parser.add_argument("--gmm-covariance-type", default="diag")
    parser.add_argument("--gmm-reg-covar", type=float, default=1e-4)
    parser.add_argument("--gmm-fallback-covariance-type", default="diag")
    parser.add_argument("--hdbscan-min-cluster-size", type=int, default=60)
    parser.add_argument("--hdbscan-min-samples", type=int, default=10)
    parser.add_argument("--spectral-neighbors", type=int, default=15)
    parser.add_argument("--silhouette-sample-size", type=int, default=5000)
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--no-scale", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
