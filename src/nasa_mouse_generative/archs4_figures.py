"""Create held-out ARCHS4 tissue figures for the Lacan diffusion model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    silhouette_score,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .adapters import load_adapter
from .adapters.diffusion import DiffusionAdapter
from .config import load_config
from .training_data import DataPartition, load_prepared_osdr, prepare_training_data


def _subset(partition: DataPartition, indices: np.ndarray) -> DataPartition:
    indices = np.asarray(indices, dtype=np.int64)
    return DataPartition(
        name=partition.name,
        matrix=partition.matrix[indices],
        obs=partition.obs.iloc[indices].reset_index(drop=True),
        categories=partition.categories[indices],
        weights=partition.weights[indices],
    )


def _balanced_indices(
    obs: pd.DataFrame, *, per_tissue: int, seed: int
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    for _, group in obs.groupby("tissue", sort=True):
        positions = group.index.to_numpy(dtype=np.int64)
        take = min(len(positions), int(per_tissue))
        if take:
            selected.extend(rng.choice(positions, take, replace=False).tolist())
    return np.asarray(sorted(selected), dtype=np.int64)


def _classification_metrics(
    train_values: np.ndarray,
    train_labels: np.ndarray,
    evaluation_values: np.ndarray,
    evaluation_labels: np.ndarray,
) -> dict[str, float | int]:
    train_labels = np.asarray(train_labels, dtype=str)
    evaluation_labels = np.asarray(evaluation_labels, dtype=str)
    classifier = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            max_iter=4000,
            class_weight="balanced",
            random_state=0,
        ),
    )
    classifier.fit(train_values, train_labels)
    predictions = classifier.predict(evaluation_values)
    return {
        "train_samples": int(len(train_values)),
        "evaluation_samples": int(len(evaluation_values)),
        "classes": int(len(np.unique(evaluation_labels))),
        "accuracy": float(accuracy_score(evaluation_labels, predictions)),
        "balanced_accuracy": float(
            balanced_accuracy_score(evaluation_labels, predictions)
        ),
        "macro_f1": float(
            f1_score(evaluation_labels, predictions, average="macro")
        ),
    }


def _silhouette(values: np.ndarray, labels: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=str)
    if len(values) < 3 or len(np.unique(labels)) < 2:
        return float("nan")
    return float(silhouette_score(values, labels))


def _matrix_stats(values: np.ndarray) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    return {
        "minimum": float(array.min()),
        "maximum": float(array.max()),
        "mean": float(array.mean()),
        "standard_deviation": float(array.std()),
    }


def _correlation(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=np.float64)
    second = np.asarray(second, dtype=np.float64)
    if first.size < 2 or first.std() == 0 or second.std() == 0:
        return float("nan")
    return float(np.corrcoef(first, second)[0, 1])


def _palette(labels: np.ndarray) -> dict[str, object]:
    import matplotlib.pyplot as plt

    names = sorted(set(map(str, labels)))
    color_map = plt.get_cmap("tab20", max(len(names), 1))
    return {name: color_map(index) for index, name in enumerate(names)}


def _write_coordinates(
    path: Path,
    obs: pd.DataFrame,
    coordinates: np.ndarray,
    *,
    first: str,
    second: str,
    source: str,
) -> Path:
    table = obs.reset_index(drop=True).copy()
    table[first] = coordinates[:, 0]
    table[second] = coordinates[:, 1]
    table["plot_source"] = source
    table.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def _diffusion_figure(
    adapter: DiffusionAdapter,
    train: DataPartition,
    evaluation: DataPartition,
    output: Path,
    *,
    seed: int,
    background_samples: int,
) -> dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(seed)
    background_indices = np.arange(len(train), dtype=np.int64)
    if len(background_indices) > background_samples:
        background_indices = rng.choice(
            background_indices, background_samples, replace=False
        )
    background = train.matrix[background_indices]
    pca = PCA(n_components=2, random_state=seed).fit(background)
    background_coordinates = pca.transform(background)
    total_timesteps = int(adapter.model_config["num_timesteps"])
    middle_timestep = 200 if total_timesteps >= 200 else max(1, total_timesteps // 5)
    requested = (total_timesteps, middle_timestep, 0)
    trajectory = adapter.generate_trajectory(
        evaluation.categories,
        seed=seed + 31,
        snapshot_timesteps=requested,
        sample_steps=adapter.model_config["num_timesteps"],
    )
    labels = evaluation.obs["tissue"].astype(str).to_numpy()
    colors = _palette(labels)
    figure, axes = plt.subplots(1, 3, figsize=(15.6, 5.0), sharex=True, sharey=True)
    coordinate_paths: dict[str, str] = {}
    for axis, timestep in zip(axes, requested):
        coordinates = pca.transform(trajectory[timestep])
        axis.scatter(
            background_coordinates[:, 0],
            background_coordinates[:, 1],
            s=5,
            alpha=0.10,
            color="#606060",
            edgecolors="none",
            label="real ARCHS4 train",
            rasterized=True,
        )
        for tissue in sorted(colors):
            mask = labels == tissue
            axis.scatter(
                coordinates[mask, 0],
                coordinates[mask, 1],
                s=11,
                alpha=0.82,
                color=colors[tissue],
                edgecolors="none",
                label=tissue,
                rasterized=True,
            )
        axis.set_title(f"t = {timestep}", fontweight="bold")
        axis.set_xlabel("PC1", fontweight="bold")
        path = _write_coordinates(
            output / f"diffusion_archs4_tissue_pca_t{timestep}.tsv.gz",
            evaluation.obs,
            coordinates,
            first="PC1",
            second="PC2",
            source=f"synthetic_t{timestep}",
        )
        coordinate_paths[str(timestep)] = str(path)
    axes[0].set_ylabel("PC2", fontweight="bold")
    handles, legend_labels = axes[-1].get_legend_handles_labels()
    figure.legend(
        handles,
        legend_labels,
        frameon=False,
        fontsize=8,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.03),
    )
    figure.suptitle(
        "Lacan-style DDIM trajectory - ARCHS4 mouse tissues",
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0.10, 1, 0.95))
    figure_path = output / "diffusion_archs4_tissue_trajectory_pca.png"
    figure.savefig(figure_path, dpi=220, bbox_inches="tight")
    plt.close(figure)

    background_table = train.obs.iloc[background_indices].reset_index(drop=True)
    background_path = _write_coordinates(
        output / "diffusion_archs4_tissue_pca_real_background.tsv.gz",
        background_table,
        background_coordinates,
        first="PC1",
        second="PC2",
        source="real_train",
    )
    real_baseline = _classification_metrics(
        background,
        background_table["tissue"].astype(str).to_numpy(),
        evaluation.matrix,
        labels,
    )
    synthetic_tstr = _classification_metrics(
        trajectory[0], labels, evaluation.matrix, labels
    )
    synthetic_final_coordinates = pca.transform(trajectory[0])
    trajectory_stats = {
        str(timestep): _matrix_stats(values)
        for timestep, values in trajectory.items()
    }
    moment_fidelity = {
        "gene_mean_correlation": _correlation(
            evaluation.matrix.mean(axis=0), trajectory[0].mean(axis=0)
        ),
        "gene_std_correlation": _correlation(
            evaluation.matrix.std(axis=0), trajectory[0].std(axis=0)
        ),
    }
    return {
        "method": (
            "DDIM x_t snapshots projected through one PCA fitted only to real "
            "ARCHS4 training profiles"
        ),
        "figure": str(figure_path),
        "coordinates": coordinate_paths,
        "real_background_coordinates": str(background_path),
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "real_train_to_heldout_tissue_metrics": real_baseline,
        "synthetic_train_to_real_heldout_tissue_metrics": synthetic_tstr,
        "synthetic_t0_pca_silhouette": _silhouette(
            synthetic_final_coordinates, labels
        ),
        "real_heldout_expression_stats": _matrix_stats(evaluation.matrix),
        "trajectory_expression_stats": trajectory_stats,
        "synthetic_t0_moment_fidelity": moment_fidelity,
    }


def run(args: argparse.Namespace) -> Path:
    run_dir = Path(args.run_dir)
    config = load_config(run_dir / "resolved_config.yaml")
    if config.training.regime != "archs4_only":
        raise ValueError("ARCHS4 tissue figures require training.regime=archs4_only")
    if (
        args.split == "test"
        and config.validation.final_test_locked
        and not args.unlock_test
    ):
        raise SystemExit(
            "The final ARCHS4 test split is locked. Pass --unlock-test only after "
            "the model and plot settings are fixed."
        )
    if (run_dir / "prepared_data.h5").exists() or (
        run_dir / "prepared_osdr.h5"
    ).exists():
        _, partitions = load_prepared_osdr(run_dir)
    else:
        prepared = prepare_training_data(config)
        partitions = prepared.partitions
    adapter = load_adapter(
        run_dir, device_spec=args.device or config.execution.device
    )
    evaluation_indices = _balanced_indices(
        partitions[args.split].obs,
        per_tissue=args.samples_per_tissue,
        seed=config.training.seed + 2,
    )
    evaluation = _subset(partitions[args.split], evaluation_indices)
    if evaluation.obs["tissue"].nunique() < 2:
        raise ValueError("At least two held-out tissues are required for this figure")
    evaluation_tissues = set(evaluation.obs["tissue"].astype(str))
    eligible_train = np.flatnonzero(
        partitions["train"].obs["tissue"].astype(str).isin(evaluation_tissues)
    )
    train_pool = _subset(partitions["train"], eligible_train)
    train_indices = _balanced_indices(
        train_pool.obs,
        per_tissue=args.train_per_tissue,
        seed=config.training.seed + 1,
    )
    train = _subset(train_pool, train_indices)

    output = run_dir / "figures" / f"archs4_tissues_{args.split}"
    output.mkdir(parents=True, exist_ok=True)
    if isinstance(adapter, DiffusionAdapter):
        result = _diffusion_figure(
            adapter,
            train,
            evaluation,
            output,
            seed=config.training.seed,
            background_samples=args.background_samples,
        )
    else:
        raise ValueError("ARCHS4 tissue figures require lacan_diffusion")
    summary = {
        "run_dir": str(run_dir),
        "adapter_id": adapter.adapter_id,
        "regime": config.training.regime,
        "split": args.split,
        "split_unit": "held-out GEO series",
        "tissues": sorted(evaluation.obs["tissue"].astype(str).unique()),
        "balanced_training_profiles": len(train),
        "balanced_evaluation_profiles": len(evaluation),
        "device": adapter.device_summary(),
        "result": result,
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "README.md").write_text(
        "# ARCHS4 Tissue Figure\n\n"
        f"Model: `{adapter.adapter_id}`  \n"
        f"Evaluation split: `{args.split}` with GEO series held out  \n"
        f"Tissues: {', '.join(summary['tissues'])}  \n"
        f"Training profiles plotted/scored: {len(train)}  \n"
        f"Held-out profiles plotted/scored: {len(evaluation)}\n\n"
        "See `summary.json` for quantitative tissue-separation metrics.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--split", choices=["validation", "test"], default="validation")
    parser.add_argument("--unlock-test", action="store_true")
    parser.add_argument("--samples-per-tissue", type=int, default=50)
    parser.add_argument("--train-per-tissue", type=int, default=300)
    parser.add_argument("--background-samples", type=int, default=5000)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
