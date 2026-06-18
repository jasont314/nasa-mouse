"""Post-fine-tuning summaries for the mouse GLARE workflow."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re

from .io import load_matrix_bundle, require_import

os.environ.setdefault("LOKY_MAX_CPU_COUNT", str(os.cpu_count() or 1))


DEFAULT_TARGET_MANIFEST = "data/processed/tms_facs_osdr_aligned.target.manifest.json"
DEFAULT_OSDR_SAMPLE_KEY = "/meta/info/id.sample name"
DEFAULT_EXPRESSION_GROUP_COLS = [
    "condition_inferred",
    "flight_status_inferred",
    "id.accession",
    "investigation.study assays.study assay technology type",
]


def _decode_array(values) -> list[str]:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8", "replace"))
        else:
            decoded.append(str(value))
    return decoded


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return cleaned[:80] or "value"


def infer_condition(profile: str) -> str:
    tokens = {
        token
        for token in re.split(r"[^A-Z0-9]+", str(profile).upper())
        if token
    }
    if tokens & {"FLT", "FLIGHT", "SPACEFLIGHT"}:
        return "flight"
    if tokens & {"GC", "GROUND", "GRD"}:
        return "ground_control"
    if tokens & {"BSL", "BASELINE"}:
        return "baseline_control"
    if tokens & {"AEM"}:
        return "aem_control"
    if tokens & {"VIV", "VIVARIUM"}:
        return "vivarium_control"
    return "unknown"


def infer_flight_status(condition: str) -> str:
    if condition == "flight":
        return "flight"
    if condition.endswith("_control"):
        return "ground_or_control"
    return "unknown"


def load_osdr_metadata(
    osdr_h5: str | Path,
    profiles: list[str],
    sample_key: str = DEFAULT_OSDR_SAMPLE_KEY,
    source_profile_indices=None,
):
    """Load 1D /meta/info arrays that align with profile order."""
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = {}
    with h5py.File(osdr_h5, "r") as handle:
        info = handle.get("/meta/info")
        if info is None:
            return pd.DataFrame({"profile": profiles})

        sample_col = Path(sample_key).name
        if sample_col not in info:
            return pd.DataFrame({"profile": profiles})
        source_length = len(info[sample_col])
        if source_profile_indices is not None:
            source_profile_indices = list(map(int, source_profile_indices))
            if len(source_profile_indices) != len(profiles):
                raise ValueError(
                    "source_profile_indices must match the requested profile count"
                )
            if source_profile_indices and max(source_profile_indices) >= source_length:
                raise ValueError(
                    "source_profile_index exceeds the OSDR HDF5 metadata length"
                )

        for key, dataset in info.items():
            shape = getattr(dataset, "shape", ())
            if len(shape) != 1:
                continue
            if source_profile_indices is not None and shape[0] == source_length:
                values = dataset[:]
                rows[key] = _decode_array(values[source_profile_indices])
            elif source_profile_indices is None and shape[0] == len(profiles):
                rows[key] = _decode_array(dataset[:])

    metadata = pd.DataFrame(rows)
    if sample_col in metadata:
        metadata.insert(0, "profile", metadata[sample_col].astype(str))
    else:
        metadata.insert(0, "profile", profiles)
    return metadata


def merge_profile_metadata(bundle, osdr_h5: str | Path | None, sample_key: str):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    if bundle.profile_metadata is not None:
        metadata = bundle.profile_metadata.copy()
    else:
        metadata = pd.DataFrame({"profile": bundle.profiles})
    if "profile" not in metadata:
        metadata.insert(0, "profile", bundle.profiles)
    metadata = metadata.reset_index(drop=True)

    if osdr_h5:
        source_profile_indices = (
            metadata["source_profile_index"].tolist()
            if "source_profile_index" in metadata
            else None
        )
        h5_metadata = load_osdr_metadata(
            osdr_h5,
            bundle.profiles,
            sample_key,
            source_profile_indices=source_profile_indices,
        )
        h5_profiles = h5_metadata["profile"].astype(str).tolist()
        bundle_profiles = [str(profile) for profile in bundle.profiles]
        if len(h5_metadata) == len(metadata) and h5_profiles == bundle_profiles:
            h5_metadata = h5_metadata.reset_index(drop=True)
            for column in h5_metadata.columns:
                if column == "profile":
                    continue
                if column in metadata:
                    metadata[column] = metadata[column].fillna(h5_metadata[column])
                else:
                    metadata[column] = h5_metadata[column]
        else:
            h5_metadata = h5_metadata.drop_duplicates(subset=["profile"], keep="first")
            metadata = metadata.merge(
                h5_metadata,
                on="profile",
                how="left",
                suffixes=("", "_h5"),
            )
            for column in list(metadata.columns):
                if not column.endswith("_h5"):
                    continue
                base = column[:-3]
                if base in metadata:
                    metadata[base] = metadata[base].fillna(metadata[column])
                    metadata = metadata.drop(columns=[column])
                else:
                    metadata = metadata.rename(columns={column: base})

    if len(metadata) != len(bundle.profiles):
        raise SystemExit(
            "Profile metadata row count changed during merge: "
            f"{len(metadata)} rows for {len(bundle.profiles)} profiles"
        )
    metadata["profile"] = bundle.profiles

    metadata["condition_inferred"] = metadata["profile"].map(infer_condition)
    metadata["flight_status_inferred"] = metadata["condition_inferred"].map(
        infer_flight_status
    )
    return metadata


def resolve_entity_axis(representation, bundle, requested_axis: str) -> str:
    if requested_axis != "auto":
        expected = len(bundle.genes) if requested_axis == "genes" else len(bundle.profiles)
        if representation.shape[0] != expected:
            raise SystemExit(
                f"--entity-axis {requested_axis} expects {expected} rows, but "
                f"representation has {representation.shape[0]} rows"
            )
        return requested_axis

    matches_genes = representation.shape[0] == len(bundle.genes)
    matches_profiles = representation.shape[0] == len(bundle.profiles)
    if matches_genes and not matches_profiles:
        return "genes"
    if matches_profiles and not matches_genes:
        return "profiles"
    if matches_genes and matches_profiles:
        return "genes"
    raise SystemExit(
        "Representation row count does not match target genes or profiles: "
        f"{representation.shape[0]} rows, {len(bundle.genes)} genes, "
        f"{len(bundle.profiles)} profiles"
    )


def scaled_latent(representation, scale: bool):
    if not scale:
        return representation
    StandardScaler = require_import(
        "sklearn.preprocessing", "pip install -r requirements-nasa-mouse-glare.txt"
    ).StandardScaler
    return StandardScaler().fit_transform(representation)


def cluster_latent(latent, n_clusters: int, seed: int):
    KMeans = require_import(
        "sklearn.cluster", "pip install -r requirements-nasa-mouse-glare.txt"
    ).KMeans
    n_clusters = min(n_clusters, latent.shape[0])
    if n_clusters < 2:
        raise SystemExit("Need at least two rows to cluster latent representations")
    model = KMeans(n_clusters=n_clusters, random_state=seed, n_init=20)
    return model.fit_predict(latent), model


def compute_pca(latent, n_components: int, seed: int):
    PCA = require_import(
        "sklearn.decomposition", "pip install -r requirements-nasa-mouse-glare.txt"
    ).PCA
    n_components = min(n_components, latent.shape[0], latent.shape[1])
    if n_components < 1:
        raise SystemExit("Cannot compute PCA for an empty latent matrix")
    model = PCA(n_components=n_components, random_state=seed)
    return model.fit_transform(latent), model


def compute_silhouette(latent, clusters, sample_size: int, seed: int) -> float | None:
    metrics = require_import(
        "sklearn.metrics", "pip install -r requirements-nasa-mouse-glare.txt"
    )
    n_clusters = len(set(clusters))
    if n_clusters < 2 or latent.shape[0] <= n_clusters:
        return None
    kwargs = {"metric": "euclidean", "random_state": seed}
    if sample_size and latent.shape[0] > sample_size:
        kwargs["sample_size"] = sample_size
    return float(metrics.silhouette_score(latent, clusters, **kwargs))


def write_entity_tables(
    representation,
    pca_values,
    pca_model,
    clusters,
    entity_axis: str,
    bundle,
    profile_metadata,
    output_dir: Path,
):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    if entity_axis == "genes":
        ids = bundle.genes
        id_col = "gene"
    else:
        ids = bundle.profiles
        id_col = "profile"

    cluster_col = f"{entity_axis[:-1]}_cluster"
    latent_cols = [f"z{i + 1}" for i in range(representation.shape[1])]
    pca_cols = [f"pc{i + 1}" for i in range(pca_values.shape[1])]

    latent_df = pd.DataFrame(representation, columns=latent_cols)
    latent_df.insert(0, id_col, ids)
    latent_df.insert(1, cluster_col, clusters)

    pca_df = pd.DataFrame(pca_values, columns=pca_cols)
    pca_df.insert(0, id_col, ids)
    pca_df.insert(1, cluster_col, clusters)

    if entity_axis == "profiles":
        latent_df = latent_df.merge(profile_metadata, on="profile", how="left")
        pca_df = pca_df.merge(profile_metadata, on="profile", how="left")

    cluster_df = latent_df[[id_col, cluster_col]].copy()
    summary = (
        pca_df.groupby(cluster_col)
        .size()
        .rename("n_entities")
        .reset_index()
        .sort_values(cluster_col)
    )
    for col in pca_cols[:2]:
        summary[f"{col}_mean"] = (
            pca_df.groupby(cluster_col)[col].mean().reindex(summary[cluster_col]).values
        )

    variance_df = pd.DataFrame(
        {
            "component": np.arange(1, len(pca_model.explained_variance_ratio_) + 1),
            "explained_variance_ratio": pca_model.explained_variance_ratio_,
        }
    )

    paths = {
        "latent_table": output_dir / f"{entity_axis[:-1]}_latent.tsv",
        "pca_table": output_dir / f"{entity_axis[:-1]}_pca.tsv",
        "clusters": output_dir / f"{entity_axis[:-1]}_clusters.tsv",
        "cluster_summary": output_dir / f"{entity_axis[:-1]}_cluster_summary.tsv",
        "pca_variance": output_dir / "pca_variance.tsv",
    }
    latent_df.to_csv(paths["latent_table"], sep="\t", index=False)
    pca_df.to_csv(paths["pca_table"], sep="\t", index=False)
    cluster_df.to_csv(paths["clusters"], sep="\t", index=False)
    summary.to_csv(paths["cluster_summary"], sep="\t", index=False)
    variance_df.to_csv(paths["pca_variance"], sep="\t", index=False)
    return paths, pca_df, cluster_col


def write_profile_cluster_crosstabs(
    pca_df,
    cluster_col: str,
    output_dir: Path,
    max_groups: int,
):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    paths = {}
    metadata_cols = [
        col
        for col in pca_df.columns
        if col not in {"profile", cluster_col} and not col.startswith(("z", "pc"))
    ]
    for col in metadata_cols:
        values = pca_df[col].fillna("NA").astype(str)
        if values.nunique() < 2 or values.nunique() > max_groups:
            continue
        table = pd.crosstab(pca_df[cluster_col], values)
        path = output_dir / f"profile_cluster_by_{safe_name(col)}.tsv"
        table.to_csv(path, sep="\t")
        paths[f"profile_cluster_by_{col}"] = path
    return paths


def write_gene_expression_summaries(
    bundle,
    profile_metadata,
    gene_clusters,
    output_dir: Path,
    group_cols: list[str],
    max_groups: int,
):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_sparse = require_import(
        "scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt"
    )

    matrix = bundle.matrix
    clusters = np.asarray(gene_clusters)
    unique_clusters = sorted(set(clusters.tolist()))
    paths = {}

    for col in group_cols:
        if col not in profile_metadata:
            continue
        groups = profile_metadata[col].fillna("NA").astype(str)
        group_values = sorted(groups.unique().tolist())
        if len(group_values) < 2 or len(group_values) > max_groups:
            continue

        rows = []
        for cluster in unique_clusters:
            gene_mask = clusters == cluster
            cluster_matrix = matrix[gene_mask, :]
            if scipy_sparse.issparse(cluster_matrix):
                sample_means = np.asarray(cluster_matrix.mean(axis=0)).ravel()
            else:
                sample_means = np.asarray(cluster_matrix.mean(axis=0)).ravel()
            row = {"gene_cluster": cluster, "n_genes": int(gene_mask.sum())}
            for value in group_values:
                sample_mask = groups.eq(value).to_numpy()
                row[value] = float(sample_means[sample_mask].mean())
            rows.append(row)

        df = pd.DataFrame(rows)
        path = output_dir / f"gene_cluster_expression_by_{safe_name(col)}.tsv"
        df.to_csv(path, sep="\t", index=False)
        paths[f"gene_cluster_expression_by_{col}"] = path

    if "flight_status_inferred" in profile_metadata:
        status = profile_metadata["flight_status_inferred"].fillna("unknown").astype(str)
        flight_mask = status.eq("flight").to_numpy()
        ground_mask = status.eq("ground_or_control").to_numpy()
        if flight_mask.any() and ground_mask.any():
            rows = []
            for cluster in unique_clusters:
                gene_mask = clusters == cluster
                cluster_matrix = matrix[gene_mask, :]
                sample_means = np.asarray(cluster_matrix.mean(axis=0)).ravel()
                flight_mean = float(sample_means[flight_mask].mean())
                ground_mean = float(sample_means[ground_mask].mean())
                rows.append(
                    {
                        "gene_cluster": cluster,
                        "n_genes": int(gene_mask.sum()),
                        "flight_mean": flight_mean,
                        "ground_or_control_mean": ground_mean,
                        "flight_minus_ground_or_control": flight_mean - ground_mean,
                    }
                )
            path = output_dir / "gene_cluster_flight_ground_summary.tsv"
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
            paths["gene_cluster_flight_ground_summary"] = path

    return paths


def maybe_write_pca_plot(pca_df, cluster_col: str, output_dir: Path):
    mpl_config_dir = output_dir / ".matplotlib"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(mpl_config_dir))
    try:
        matplotlib = require_import("matplotlib", "pip install matplotlib")
        matplotlib.use("Agg")
        plt = require_import("matplotlib.pyplot", "pip install matplotlib")
    except SystemExit:
        return {}

    if "pc1" not in pca_df or "pc2" not in pca_df:
        return {}

    path = output_dir / "pca_by_cluster.png"
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    clusters = sorted(pca_df[cluster_col].unique().tolist())
    cmap = plt.get_cmap("tab20")
    for idx, cluster in enumerate(clusters):
        sub = pca_df[pca_df[cluster_col] == cluster]
        ax.scatter(
            sub["pc1"],
            sub["pc2"],
            s=5,
            alpha=0.65,
            color=cmap(idx % 20),
            label=str(cluster),
        )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Fine-tuned GLARE latent PCA")
    if len(clusters) <= 30:
        ax.legend(title=cluster_col, markerscale=2, fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return {"pca_by_cluster_plot": path}


def maybe_write_tsne(
    latent,
    ids,
    clusters,
    entity_axis: str,
    output_dir: Path,
    perplexity: float,
    seed: int,
):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    TSNE = require_import(
        "sklearn.manifold", "pip install -r requirements-nasa-mouse-glare.txt"
    ).TSNE

    if latent.shape[0] <= 3:
        return {}
    perplexity = min(perplexity, max(1.0, (latent.shape[0] - 1) / 3))
    values = TSNE(
        n_components=2,
        perplexity=perplexity,
        init="pca",
        learning_rate="auto",
        random_state=seed,
    ).fit_transform(latent)
    id_col = "gene" if entity_axis == "genes" else "profile"
    cluster_col = f"{entity_axis[:-1]}_cluster"
    df = pd.DataFrame(
        {
            id_col: ids,
            cluster_col: clusters,
            "tsne1": values[:, 0],
            "tsne2": values[:, 1],
        }
    )
    path = output_dir / f"{entity_axis[:-1]}_tsne.tsv"
    df.to_csv(path, sep="\t", index=False)
    return {"tsne_table": path}


def path_strings(paths: dict[str, Path]) -> dict[str, str]:
    return {key: str(value) for key, value in paths.items()}


def run(args) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bundle = load_matrix_bundle(args.target_manifest)
    representation = np.load(args.representation)
    if representation.ndim != 2:
        raise SystemExit(
            f"Expected 2D representation array, got shape {representation.shape}"
        )

    entity_axis = resolve_entity_axis(representation, bundle, args.entity_axis)
    profile_metadata = merge_profile_metadata(bundle, args.osdr, args.sample_key)
    profile_metadata_path = output_dir / "profile_metadata.tsv"
    profile_metadata.to_csv(profile_metadata_path, sep="\t", index=False)

    latent_for_analysis = scaled_latent(representation, scale=not args.no_scale)
    clusters, _ = cluster_latent(latent_for_analysis, args.n_clusters, args.seed)
    pca_values, pca_model = compute_pca(
        latent_for_analysis,
        args.pca_components,
        args.seed,
    )
    silhouette = compute_silhouette(
        latent_for_analysis,
        clusters,
        args.silhouette_sample_size,
        args.seed,
    )

    paths, pca_df, cluster_col = write_entity_tables(
        representation,
        pca_values,
        pca_model,
        clusters,
        entity_axis,
        bundle,
        profile_metadata,
        output_dir,
    )
    paths["profile_metadata"] = profile_metadata_path

    if entity_axis == "profiles":
        paths.update(
            write_profile_cluster_crosstabs(
                pca_df,
                cluster_col,
                output_dir,
                args.max_expression_groups,
            )
        )
    elif not args.skip_expression_summary:
        group_cols = args.expression_group_cols or DEFAULT_EXPRESSION_GROUP_COLS
        paths.update(
            write_gene_expression_summaries(
                bundle,
                profile_metadata,
                clusters,
                output_dir,
                group_cols,
                args.max_expression_groups,
            )
        )

    if not args.no_plots:
        paths.update(maybe_write_pca_plot(pca_df, cluster_col, output_dir))

    if args.run_tsne:
        ids = bundle.genes if entity_axis == "genes" else bundle.profiles
        paths.update(
            maybe_write_tsne(
                latent_for_analysis,
                ids,
                clusters,
                entity_axis,
                output_dir,
                args.tsne_perplexity,
                args.seed,
            )
        )

    summary = {
        "representation": str(args.representation),
        "target_manifest": str(args.target_manifest),
        "osdr": str(args.osdr) if args.osdr else "",
        "entity_axis": entity_axis,
        "representation_shape": list(representation.shape),
        "target_shape": list(bundle.matrix.shape),
        "n_clusters": int(len(set(clusters.tolist()))),
        "seed": args.seed,
        "scaled_latent": not args.no_scale,
        "silhouette": silhouette,
        "outputs": path_strings(paths),
    }
    summary_path = output_dir / "post_finetune_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"entity_axis={entity_axis}")
    print(f"representation_shape={representation.shape}")
    print(f"n_clusters={summary['n_clusters']}")
    if silhouette is not None:
        print(f"silhouette={silhouette:.4f}")
    print(f"summary={summary_path}")
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize GLARE fine-tuned latent representations."
    )
    parser.add_argument("--representation", required=True, help="Path to .npy latent array.")
    parser.add_argument(
        "--target-manifest",
        default=DEFAULT_TARGET_MANIFEST,
        help="Aligned OSDR target manifest used to export the fine-tuning CSV.",
    )
    parser.add_argument(
        "--osdr",
        default="",
        help="Optional OSDR HDF5 file for additional /meta/info sample metadata.",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--entity-axis",
        choices=["auto", "genes", "profiles"],
        default="auto",
        help="Whether latent rows represent genes or profiles. Auto detects from manifest.",
    )
    parser.add_argument("--n-clusters", type=int, default=15)
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--pca-components", type=int, default=10)
    parser.add_argument("--silhouette-sample-size", type=int, default=5000)
    parser.add_argument(
        "--expression-group-cols",
        nargs="*",
        default=None,
        help="Profile metadata columns for gene-cluster expression summaries.",
    )
    parser.add_argument("--max-expression-groups", type=int, default=250)
    parser.add_argument("--skip-expression-summary", action="store_true")
    parser.add_argument("--no-scale", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--run-tsne", action="store_true")
    parser.add_argument("--tsne-perplexity", type=float, default=30.0)
    parser.add_argument("--sample-key", default=DEFAULT_OSDR_SAMPLE_KEY)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
