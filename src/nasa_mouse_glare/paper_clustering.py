"""GLARE base clustering and exact evidence-accumulation consensus clustering."""

from __future__ import annotations

import argparse
import heapq
import json
import time
from pathlib import Path

import hdbscan
import matplotlib
import numpy as np
import pandas as pd
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture

matplotlib.use("Agg")
import matplotlib.pyplot as plt


LOCATION_CONFIG = {
    "FLT": {
        "gmm_clusters": 20,
        "hdbscan_min_cluster_size": 60,
        "spectral_clusters": 25,
        "consensus_clusters": 16,
    },
    "GC": {
        "gmm_clusters": 25,
        "hdbscan_min_cluster_size": 50,
        "spectral_clusters": 20,
        "consensus_clusters": 15,
    },
}


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def compressed_eac(base_labels: np.ndarray, n_clusters: int) -> np.ndarray:
    """Run exact average-linkage EAC after lossless signature compression.

    Genes with identical labels in every base partition have zero EAC distance
    and identical distance to every other gene. Collapsing those duplicates and
    retaining their counts therefore gives the same average-linkage updates as
    the full gene-by-gene co-association matrix without quadratic gene memory.
    """
    labels = np.asarray(base_labels)
    if labels.ndim != 2 or labels.shape[1] < 2:
        raise ValueError("base_labels must have shape (objects, clusterings)")
    signatures, inverse, counts = np.unique(
        labels, axis=0, return_inverse=True, return_counts=True
    )
    n_signatures = len(signatures)
    if not 1 <= n_clusters <= n_signatures:
        raise ValueError(
            f"n_clusters must be in [1, {n_signatures}], got {n_clusters}"
        )

    active = set(range(n_signatures))
    weights = {index: int(counts[index]) for index in active}
    members = {index: [index] for index in active}
    distances: dict[tuple[int, int], float] = {}
    heap: list[tuple[float, int, int]] = []

    for left in range(n_signatures):
        matches = np.mean(signatures[left + 1 :] == signatures[left], axis=1)
        for offset, similarity in enumerate(matches, start=left + 1):
            key = (left, offset)
            distance = float(1.0 - similarity)
            distances[key] = distance
            heapq.heappush(heap, (distance, left, offset))

    next_id = n_signatures
    while len(active) > n_clusters:
        while heap:
            distance, left, right = heapq.heappop(heap)
            if left in active and right in active:
                break
        else:
            raise RuntimeError("EAC priority queue emptied before clustering finished")

        left_weight = weights[left]
        right_weight = weights[right]
        others = list(active - {left, right})
        active.remove(left)
        active.remove(right)
        merged = next_id
        next_id += 1
        active.add(merged)
        weights[merged] = left_weight + right_weight
        members[merged] = members[left] + members[right]

        for other in others:
            left_key = (min(left, other), max(left, other))
            right_key = (min(right, other), max(right, other))
            merged_distance = (
                left_weight * distances[left_key]
                + right_weight * distances[right_key]
            ) / (left_weight + right_weight)
            merged_key = (min(merged, other), max(merged, other))
            distances[merged_key] = merged_distance
            heapq.heappush(
                heap,
                (merged_distance, merged_key[0], merged_key[1]),
            )

    signature_labels = np.empty(n_signatures, dtype=np.int32)
    final_groups = sorted(
        (members[cluster] for cluster in active),
        key=lambda group: -sum(int(counts[index]) for index in group),
    )
    for label, group in enumerate(final_groups):
        signature_labels[group] = label
    return signature_labels[inverse]


def run_base_clusterings(
    representation: np.ndarray, location: str
) -> tuple[np.ndarray, dict]:
    config = LOCATION_CONFIG[location]
    log(f"{location}: fitting {config['gmm_clusters']}-component GMM")
    gmm = GaussianMixture(
        n_components=config["gmm_clusters"],
        random_state=2024,
        n_init=config["gmm_clusters"],
        max_iter=100,
    ).fit_predict(representation)

    min_cluster_size = config["hdbscan_min_cluster_size"]
    log(f"{location}: fitting HDBSCAN min_cluster_size={min_cluster_size}")
    density_model = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_cluster_size // 2,
        cluster_selection_method="leaf",
        prediction_data=True,
    ).fit(representation)
    membership = hdbscan.all_points_membership_vectors(density_model)
    if membership.ndim != 2 or membership.shape[1] == 0:
        raise RuntimeError(f"{location}: HDBSCAN produced no soft clusters")
    density = np.argmax(membership, axis=1)

    log(f"{location}: fitting {config['spectral_clusters']}-cluster spectral model")
    spectral = SpectralClustering(
        n_clusters=config["spectral_clusters"],
        affinity="nearest_neighbors",
        n_jobs=-1,
        random_state=2024,
    ).fit_predict(representation)
    labels = np.column_stack([gmm, density, spectral])
    diagnostics = {
        "gmm_clusters_observed": int(np.unique(gmm).size),
        "hdbscan_hard_clusters_observed": int(
            np.unique(density_model.labels_[density_model.labels_ >= 0]).size
        ),
        "hdbscan_soft_clusters_observed": int(np.unique(density).size),
        "hdbscan_noise_points_before_soft_assignment": int(
            np.sum(density_model.labels_ < 0)
        ),
        "spectral_clusters_observed": int(np.unique(spectral).size),
        "unique_partition_signatures": int(np.unique(labels, axis=0).shape[0]),
    }
    return labels, diagnostics


def save_visualization(
    representation: np.ndarray,
    genes: list[str],
    consensus: np.ndarray,
    location: str,
    output_dir: Path,
    skip_tsne: bool,
) -> dict:
    pca = PCA(n_components=2, random_state=2024).fit_transform(representation)
    pca_path = output_dir / f"{location}_pca.tsv"
    pd.DataFrame(
        {
            "gene_id": genes,
            "x": pca[:, 0],
            "y": pca[:, 1],
            "consensus": consensus,
        }
    ).to_csv(pca_path, sep="\t", index=False)
    _scatter(pca, consensus, f"{location} GLARE consensus (PCA)", output_dir / f"{location}_pca.png")
    outputs = {"pca": str(pca_path), "pca_plot": str(output_dir / f"{location}_pca.png")}

    if not skip_tsne:
        log(f"{location}: computing GLARE-style 2D t-SNE visualization")
        tsne = TSNE(
            n_components=2,
            random_state=1996,
            n_jobs=-1,
            learning_rate="auto",
            init="pca",
        ).fit_transform(representation)
        tsne_path = output_dir / f"{location}_tsne.tsv"
        pd.DataFrame(
            {
                "gene_id": genes,
                "x": tsne[:, 0],
                "y": tsne[:, 1],
                "consensus": consensus,
            }
        ).to_csv(tsne_path, sep="\t", index=False)
        _scatter(
            tsne,
            consensus,
            f"{location} GLARE consensus (t-SNE)",
            output_dir / f"{location}_tsne.png",
        )
        outputs.update(
            {"tsne": str(tsne_path), "tsne_plot": str(output_dir / f"{location}_tsne.png")}
        )
    return outputs


def _scatter(
    coordinates: np.ndarray,
    labels: np.ndarray,
    title: str,
    output_path: Path,
) -> None:
    figure, axis = plt.subplots(figsize=(8, 6), constrained_layout=True)
    points = axis.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        c=labels,
        cmap="tab20",
        s=3,
        alpha=0.65,
        linewidths=0,
    )
    axis.set_title(title)
    axis.set_xlabel("Dimension 1")
    axis.set_ylabel("Dimension 2")
    figure.colorbar(points, ax=axis, label="Consensus cluster")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def cluster_location(
    representation_path: Path,
    genes: list[str],
    location: str,
    output_dir: Path,
    skip_tsne: bool,
) -> dict:
    start = time.perf_counter()
    representation = np.load(representation_path)
    if representation.shape[0] != len(genes):
        raise ValueError(
            f"{location} representation has {representation.shape[0]} rows but "
            f"there are {len(genes)} genes"
        )
    base_labels, diagnostics = run_base_clusterings(representation, location)
    requested = LOCATION_CONFIG[location]["consensus_clusters"]
    log(
        f"{location}: exact EAC average-linkage consensus -> {requested} clusters"
    )
    consensus = compressed_eac(base_labels, requested)
    clusters = pd.DataFrame(
        {
            "gene_id": genes,
            "gmm": base_labels[:, 0],
            "hdbscan": base_labels[:, 1],
            "spectral": base_labels[:, 2],
            "consensus": consensus,
        }
    )
    clusters_path = output_dir / f"{location}_gene_clusters.tsv"
    clusters.to_csv(clusters_path, sep="\t", index=False)
    summary_table = (
        clusters.groupby("consensus", as_index=False)
        .size()
        .rename(columns={"size": "gene_count"})
    )
    summary_table.to_csv(
        output_dir / f"{location}_cluster_summary.tsv", sep="\t", index=False
    )

    sample_size = min(5000, len(representation))
    sample = np.random.default_rng(2024).choice(
        len(representation), sample_size, replace=False
    )
    silhouette = silhouette_score(representation[sample], consensus[sample])
    visualization = save_visualization(
        representation, genes, consensus, location, output_dir, skip_tsne
    )
    return {
        "location": location,
        "representation": str(representation_path),
        "shape": list(representation.shape),
        "base_configuration": LOCATION_CONFIG[location],
        "diagnostics": diagnostics,
        "consensus_clusters_observed": int(np.unique(consensus).size),
        "consensus_silhouette_sample_size": sample_size,
        "consensus_silhouette": float(silhouette),
        "clusters": str(clusters_path),
        "visualization": visualization,
        "elapsed_seconds": round(time.perf_counter() - start, 3),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run paper-faithful GLARE EAC clustering.")
    parser.add_argument(
        "--run-dir", default="outputs/glare_paper_tms_liver_osd379"
    )
    parser.add_argument("--skip-tsne", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    target = np.load(run_dir / "controlled_target.npz")
    genes = target["genes"].astype(str).tolist()
    output_dir = run_dir / "clustering"
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = [
        cluster_location(
            run_dir / f"{location}_FTSAE_representation.npy",
            genes,
            location,
            output_dir,
            args.skip_tsne,
        )
        for location in ("FLT", "GC")
    ]
    summary = {
        "method": (
            "GLARE GMM/HDBSCAN/Spectral partitions followed by co-association "
            "distance and average-linkage EAC"
        ),
        "implementation": (
            "Exact lossless partition-signature compression avoids the full "
            "gene-by-gene co-association matrix"
        ),
        "locations": summaries,
    }
    (output_dir / "clustering_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    log(f"Saved clustering summary: {output_dir / 'clustering_summary.json'}")


if __name__ == "__main__":
    main()
