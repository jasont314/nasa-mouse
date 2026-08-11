"""Evaluate a frozen DDIM augmentation recipe on held-out OSDR studies."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    roc_auc_score,
)
from sklearn.preprocessing import StandardScaler


def _decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values]
    )


def _load_role(
    handle: h5py.File,
    samples: pd.DataFrame,
    role: str,
    tissue: str,
) -> tuple[np.ndarray, pd.DataFrame]:
    source_rows = np.asarray(handle[f"{role}/source_row"][:], dtype=np.int64)
    lookup = samples.set_index("_row_index", drop=False)
    missing = sorted(set(map(int, source_rows)) - set(lookup.index.astype(int)))
    if missing:
        raise ValueError(f"Metadata is missing source rows: {missing[:5]}")
    metadata = lookup.loc[source_rows].reset_index(drop=True)
    expression = np.asarray(handle[f"{role}/expression"][:], dtype=np.float32)
    keep = metadata["tissue"].astype(str).eq(tissue).to_numpy()
    return expression[keep], metadata.loc[keep].reset_index(drop=True)


def _labels(metadata: pd.DataFrame) -> np.ndarray:
    conditions = metadata["condition"].astype(str)
    unexpected = sorted(set(conditions) - {"flight", "ground_control"})
    if unexpected:
        raise ValueError(f"Unexpected conditions in evaluation data: {unexpected}")
    return conditions.eq("flight").to_numpy(dtype=np.int64)


def _metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    return {
        "balanced_accuracy": float(
            balanced_accuracy_score(labels, probabilities >= 0.5)
        ),
        "roc_auc": float(roc_auc_score(labels, probabilities)),
        "average_precision": float(average_precision_score(labels, probabilities)),
    }


def _accession_metrics(
    metadata: pd.DataFrame,
    labels: np.ndarray,
    probabilities: np.ndarray,
    arm: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for accession, indices in metadata.groupby("accession", observed=True).groups.items():
        positions = np.asarray(list(indices), dtype=np.int64)
        accession_labels = labels[positions]
        if np.unique(accession_labels).size != 2:
            continue
        rows.append(
            {
                "accession": str(accession),
                "arm": arm,
                "profiles": int(len(positions)),
                **_metrics(accession_labels, probabilities[positions]),
            }
        )
    if not rows:
        raise ValueError("No held-out accession contains both FLT and GC profiles")
    return pd.DataFrame(rows)


def _macro_metrics(table: pd.DataFrame) -> dict[str, float]:
    return {
        metric: float(table[metric].mean())
        for metric in ("balanced_accuracy", "roc_auc", "average_precision")
    }


def _generated_path(directory: Path, guidance_scale: float, seed: int) -> Path:
    scale = f"{guidance_scale:g}".replace(".", "p")
    return directory / f"scale_{scale}_seed_{seed}.npz"


def _load_recentered_draw(
    path: Path,
    source_rows: np.ndarray,
    real_expression: np.ndarray,
    labels: np.ndarray,
    residual_scale: float,
) -> np.ndarray:
    with np.load(path) as generated:
        generated_rows = np.asarray(generated["source_row"], dtype=np.int64)
        lookup = {int(row): index for index, row in enumerate(generated_rows)}
        missing = sorted(set(map(int, source_rows)) - set(lookup))
        if missing:
            raise ValueError(f"Generated draw {path} is missing rows: {missing[:5]}")
        indices = np.asarray([lookup[int(row)] for row in source_rows], dtype=np.int64)
        synthetic = np.asarray(generated["scaled_expression"][indices], dtype=np.float32)
    recentered = np.empty_like(synthetic)
    for condition in (0, 1):
        mask = labels == condition
        if not mask.any():
            raise ValueError("Both FLT and GC training profiles are required")
        real_mean = real_expression[mask].mean(axis=0)
        synthetic_mean = synthetic[mask].mean(axis=0)
        recentered[mask] = real_mean + residual_scale * (
            synthetic[mask] - synthetic_mean
        )
    if not np.isfinite(recentered).all():
        raise ValueError(f"Non-finite values in generated draw {path}")
    return recentered


def evaluate(
    *,
    fold: str,
    prepared_h5: Path,
    samples_tsv: Path,
    generated_dir: Path,
    output_dir: Path,
    seeds: list[int],
    tissue: str,
    guidance_scale: float,
    feature_count: int,
    regularization_c: float,
    synthetic_weight: float,
    residual_scale: float,
    classifier_seed: int,
) -> dict[str, object]:
    samples = pd.read_csv(samples_tsv, sep="\t")
    with h5py.File(prepared_h5, "r") as handle:
        genes = _decode(handle["genes"][:])
        train_expression, train_metadata = _load_role(
            handle, samples, "train", tissue
        )
        test_expression, test_metadata = _load_role(handle, samples, "test", tissue)

    train_labels = _labels(train_metadata)
    test_labels = _labels(test_metadata)
    scaler = StandardScaler().fit(train_expression)
    scaled_train = scaler.transform(train_expression)
    scaled_test = scaler.transform(test_expression)
    selector = SelectKBest(f_classif, k=min(feature_count, train_expression.shape[1]))
    selected_train = selector.fit_transform(scaled_train, train_labels)
    selected_test = selector.transform(scaled_test)

    classifier_options = {
        "C": regularization_c,
        "class_weight": "balanced",
        "max_iter": 5000,
        "random_state": classifier_seed,
        "solver": "lbfgs",
    }
    real_classifier = LogisticRegression(**classifier_options).fit(
        selected_train, train_labels
    )
    real_probability = real_classifier.predict_proba(selected_test)[:, 1]

    synthetic_blocks = [
        _load_recentered_draw(
            _generated_path(generated_dir, guidance_scale, seed),
            train_metadata["_row_index"].to_numpy(dtype=np.int64),
            train_expression,
            train_labels,
            residual_scale,
        )
        for seed in seeds
    ]
    synthetic_expression = np.concatenate(synthetic_blocks)
    synthetic_labels = np.tile(train_labels, len(seeds))
    selected_synthetic = selector.transform(scaler.transform(synthetic_expression))
    augmented_expression = np.concatenate((selected_train, selected_synthetic))
    augmented_labels = np.concatenate((train_labels, synthetic_labels))
    sample_weight = np.concatenate(
        (
            np.ones(len(train_labels), dtype=np.float64),
            np.full(
                len(synthetic_labels),
                synthetic_weight / len(seeds),
                dtype=np.float64,
            ),
        )
    )
    augmented_classifier = LogisticRegression(**classifier_options).fit(
        augmented_expression,
        augmented_labels,
        sample_weight=sample_weight,
    )
    augmented_probability = augmented_classifier.predict_proba(selected_test)[:, 1]

    real_accessions = _accession_metrics(
        test_metadata, test_labels, real_probability, "real_only"
    )
    augmented_accessions = _accession_metrics(
        test_metadata, test_labels, augmented_probability, "real_plus_synthetic"
    )
    real_macro = _macro_metrics(real_accessions)
    augmented_macro = _macro_metrics(augmented_accessions)
    deltas = {
        metric: augmented_macro[metric] - real_macro[metric]
        for metric in real_macro
    }
    success = bool(
        deltas["balanced_accuracy"] >= 0.03 and deltas["roc_auc"] >= -1e-12
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    per_accession = pd.concat((real_accessions, augmented_accessions), ignore_index=True)
    per_accession.to_csv(output_dir / "per_accession.tsv", sep="\t", index=False)
    predictions = test_metadata[
        ["_row_index", "profile_id", "accession", "condition", "tissue"]
    ].copy()
    predictions["label"] = test_labels
    predictions["real_probability"] = real_probability
    predictions["augmented_probability"] = augmented_probability
    predictions.to_csv(output_dir / "predictions.tsv.gz", sep="\t", index=False)
    pd.DataFrame(
        {"gene": genes[selector.get_support()], "f_score": selector.scores_[selector.get_support()]}
    ).sort_values("f_score", ascending=False).to_csv(
        output_dir / "selected_features.tsv", sep="\t", index=False
    )

    summary: dict[str, object] = {
        "fold": fold,
        "tissue": tissue,
        "protocol": "pooled_training_with_accession_heldout_test",
        "train_profiles": int(len(train_labels)),
        "test_profiles": int(len(test_labels)),
        "test_accessions": sorted(test_metadata["accession"].astype(str).unique()),
        "recipe": {
            "guidance_scale": guidance_scale,
            "seeds": seeds,
            "variant": "condition_recentered_residual",
            "residual_scale": residual_scale,
            "total_synthetic_weight": synthetic_weight,
            "feature_count": int(selector.get_support().sum()),
            "regularization_c": regularization_c,
            "classifier": "standard_scaler_f_classif_l2_logistic_regression",
        },
        "real_only": {
            "pooled": _metrics(test_labels, real_probability),
            "accession_macro": real_macro,
        },
        "real_plus_synthetic": {
            "pooled": _metrics(test_labels, augmented_probability),
            "accession_macro": augmented_macro,
        },
        "accession_macro_delta": deltas,
        "success_rule": "balanced_accuracy_delta>=0.03 and roc_auc_delta>=0",
        "success": success,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fold", required=True)
    parser.add_argument("--prepared-h5", required=True, type=Path)
    parser.add_argument("--samples-tsv", required=True, type=Path)
    parser.add_argument("--generated-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seeds", required=True, nargs="+", type=int)
    parser.add_argument("--tissue", default="skeletal_muscle")
    parser.add_argument("--guidance-scale", default=2.0, type=float)
    parser.add_argument("--feature-count", default=100, type=int)
    parser.add_argument("--regularization-c", default=1.0, type=float)
    parser.add_argument("--synthetic-weight", default=0.05, type=float)
    parser.add_argument("--residual-scale", default=1.0, type=float)
    parser.add_argument("--classifier-seed", default=1234, type=int)
    return parser


def main() -> None:
    arguments = vars(build_parser().parse_args())
    summary = evaluate(**arguments)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
