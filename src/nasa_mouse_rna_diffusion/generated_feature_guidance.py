"""Cross-fit DDIM-guided feature selection on held-out OSDR accessions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.stats import hypergeom, rankdata
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    roc_auc_score,
)
from sklearn.preprocessing import StandardScaler


METRICS = ("balanced_accuracy", "roc_auc", "average_precision")
RANK_METHODS = (
    "synthetic_f",
    "f_rank_consensus",
    "effect_consensus",
    "coefficient_consensus",
)


def _decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values]
    )


def _labels(metadata: pd.DataFrame) -> np.ndarray:
    condition = metadata["condition"].astype(str)
    unexpected = sorted(set(condition) - {"flight", "ground_control"})
    if unexpected:
        raise ValueError(f"Unexpected conditions: {unexpected}")
    return condition.eq("flight").to_numpy(dtype=np.int64)


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


def _metric_set(labels: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    if np.unique(labels).size != 2:
        raise ValueError("FLT/GC metrics require both classes")
    return {
        "balanced_accuracy": float(
            balanced_accuracy_score(labels, probability >= 0.5)
        ),
        "roc_auc": float(roc_auc_score(labels, probability)),
        "average_precision": float(average_precision_score(labels, probability)),
    }


def _accession_metrics(
    metadata: pd.DataFrame,
    labels: np.ndarray,
    probability: np.ndarray,
    arm: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for accession, indices in metadata.groupby("accession", observed=True).groups.items():
        positions = np.asarray(list(indices), dtype=np.int64)
        subset_labels = labels[positions]
        if np.unique(subset_labels).size != 2:
            continue
        rows.append(
            {
                "accession": str(accession),
                "arm": arm,
                "profiles": int(len(positions)),
                **_metric_set(subset_labels, probability[positions]),
            }
        )
    if not rows:
        raise ValueError("No accession contains both FLT and GC profiles")
    return pd.DataFrame(rows)


def _macro_metrics(table: pd.DataFrame) -> dict[str, float]:
    return {metric: float(table[metric].mean()) for metric in METRICS}


def _safe_f_score(expression: np.ndarray, labels: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        score = f_classif(expression, labels)[0]
    return np.nan_to_num(score, nan=0.0, posinf=0.0, neginf=0.0)


def _percentile_rank(values: np.ndarray) -> np.ndarray:
    values = np.nan_to_num(np.asarray(values, dtype=float), nan=0.0)
    return rankdata(values, method="average") / len(values)


def _fit_classifier(
    expression: np.ndarray,
    labels: np.ndarray,
    *,
    regularization_c: float,
    seed: int,
    sample_weight: np.ndarray | None = None,
) -> LogisticRegression:
    classifier = LogisticRegression(
        C=float(regularization_c),
        class_weight="balanced",
        max_iter=5000,
        random_state=int(seed),
        solver="lbfgs",
    )
    classifier.fit(expression, labels, sample_weight=sample_weight)
    return classifier


def _accession_effects(
    expression: np.ndarray, labels: np.ndarray, metadata: pd.DataFrame
) -> np.ndarray:
    effects: list[np.ndarray] = []
    for _, indices in metadata.groupby("accession", observed=True).groups.items():
        positions = np.asarray(list(indices), dtype=np.int64)
        subset_labels = labels[positions]
        if np.unique(subset_labels).size != 2:
            continue
        subset = expression[positions]
        effects.append(subset[subset_labels == 1].mean(0) - subset[subset_labels == 0].mean(0))
    if not effects:
        raise ValueError("No training accession has both FLT and GC")
    return np.stack(effects)


def _generated_path(directory: Path, seed: int) -> Path:
    return directory / f"scale_2_seed_{seed}.npz"


def _load_generated_draw(
    path: Path,
    source_rows: np.ndarray,
) -> np.ndarray:
    with np.load(path) as generated:
        generated_rows = np.asarray(generated["source_row"], dtype=np.int64)
        lookup = {int(row): index for index, row in enumerate(generated_rows)}
        missing = sorted(set(map(int, source_rows)) - set(lookup))
        if missing:
            raise ValueError(f"Generated draw {path} is missing rows: {missing[:5]}")
        indices = np.asarray([lookup[int(row)] for row in source_rows], dtype=np.int64)
        expression = np.asarray(
            generated["scaled_expression"][indices], dtype=np.float32
        )
    if not np.isfinite(expression).all():
        raise ValueError(f"Generated draw {path} contains non-finite values")
    return expression


def _recenter_draw(
    synthetic: np.ndarray,
    real: np.ndarray,
    labels: np.ndarray,
) -> np.ndarray:
    recentered = np.empty_like(synthetic)
    for condition in (0, 1):
        mask = labels == condition
        recentered[mask] = real[mask].mean(0) + (
            synthetic[mask] - synthetic[mask].mean(0)
        )
    return recentered


def _build_rankings(
    real_train: np.ndarray,
    labels: np.ndarray,
    metadata: pd.DataFrame,
    synthetic_draws: list[np.ndarray],
    *,
    seed: int,
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    real_f = _safe_f_score(real_train, labels)
    synthetic_f_by_seed = np.stack(
        [_safe_f_score(draw, labels) for draw in synthetic_draws]
    )
    synthetic_f = np.median(synthetic_f_by_seed, axis=0)

    real_effects = _accession_effects(real_train, labels, metadata)
    real_effect = np.median(real_effects, axis=0)
    real_stability = np.abs(np.mean(np.sign(real_effects), axis=0))
    synthetic_effects = np.stack(
        [draw[labels == 1].mean(0) - draw[labels == 0].mean(0) for draw in synthetic_draws]
    )
    synthetic_effect = np.median(synthetic_effects, axis=0)
    synthetic_stability = np.abs(np.mean(np.sign(synthetic_effects), axis=0))
    effect_direction_match = np.sign(real_effect) == np.sign(synthetic_effect)

    real_coefficient = _fit_classifier(
        real_train, labels, regularization_c=0.1, seed=seed
    ).coef_[0]
    synthetic_coefficients = np.stack(
        [
            _fit_classifier(draw, labels, regularization_c=0.1, seed=seed + offset + 1).coef_[0]
            for offset, draw in enumerate(synthetic_draws)
        ]
    )
    synthetic_coefficient = np.median(synthetic_coefficients, axis=0)
    coefficient_stability = np.abs(np.mean(np.sign(synthetic_coefficients), axis=0))
    coefficient_direction_match = (
        np.sign(real_coefficient) == np.sign(synthetic_coefficient)
    )

    rankings = {
        "real_f": real_f,
        "synthetic_f": synthetic_f,
        "f_rank_consensus": np.sqrt(
            _percentile_rank(real_f) * _percentile_rank(synthetic_f)
        ),
        "effect_consensus": np.sqrt(
            _percentile_rank(np.abs(real_effect))
            * _percentile_rank(np.abs(synthetic_effect))
        )
        * real_stability
        * synthetic_stability
        * effect_direction_match,
        "coefficient_consensus": np.sqrt(
            _percentile_rank(np.abs(real_coefficient))
            * _percentile_rank(np.abs(synthetic_coefficient))
        )
        * coefficient_stability
        * coefficient_direction_match,
    }
    diagnostics = pd.DataFrame(
        {
            "real_f": real_f,
            "synthetic_f": synthetic_f,
            "real_effect": real_effect,
            "synthetic_effect": synthetic_effect,
            "real_effect_sign_stability": real_stability,
            "synthetic_effect_sign_stability": synthetic_stability,
            "effect_direction_match": effect_direction_match,
            "real_full_coefficient": real_coefficient,
            "synthetic_full_coefficient": synthetic_coefficient,
            "synthetic_coefficient_sign_stability": coefficient_stability,
            "coefficient_direction_match": coefficient_direction_match,
        }
    )
    return rankings, diagnostics


def _selected_indices(score: np.ndarray, count: int) -> np.ndarray:
    count = min(int(count), len(score))
    return np.argsort(-np.asarray(score), kind="stable")[:count]


def _fit_candidate(
    candidate: pd.Series | dict[str, object],
    real_train: np.ndarray,
    train_labels: np.ndarray,
    synthetic_raw: np.ndarray,
    synthetic_recentered: np.ndarray,
    rankings: dict[str, np.ndarray],
    *,
    seed: int,
) -> tuple[LogisticRegression, np.ndarray]:
    rank_method = str(candidate["rank_method"])
    selected = _selected_indices(
        rankings[rank_method], int(candidate["feature_count"])
    )
    arm = str(candidate["training_arm"])
    if arm == "real_only":
        classifier = _fit_classifier(
            real_train[:, selected],
            train_labels,
            regularization_c=float(candidate["regularization_c"]),
            seed=seed,
        )
        return classifier, selected

    synthetic = synthetic_raw if arm == "raw_synthetic" else synthetic_recentered
    synthetic_labels = np.tile(train_labels, synthetic.shape[0] // len(train_labels))
    expression = np.concatenate((real_train[:, selected], synthetic[:, selected]))
    labels = np.concatenate((train_labels, synthetic_labels))
    number_of_draws = synthetic.shape[0] // len(train_labels)
    weights = np.concatenate(
        (
            np.ones(len(train_labels), dtype=float),
            np.full(
                len(synthetic_labels),
                float(candidate["synthetic_weight"]) / number_of_draws,
                dtype=float,
            ),
        )
    )
    classifier = _fit_classifier(
        expression,
        labels,
        regularization_c=float(candidate["regularization_c"]),
        seed=seed,
        sample_weight=weights,
    )
    return classifier, selected


def _candidate_metrics(
    candidate: dict[str, object],
    real_train: np.ndarray,
    train_labels: np.ndarray,
    validation: np.ndarray,
    validation_labels: np.ndarray,
    validation_metadata: pd.DataFrame,
    synthetic_raw: np.ndarray,
    synthetic_recentered: np.ndarray,
    rankings: dict[str, np.ndarray],
    *,
    seed: int,
) -> dict[str, object]:
    classifier, selected = _fit_candidate(
        candidate,
        real_train,
        train_labels,
        synthetic_raw,
        synthetic_recentered,
        rankings,
        seed=seed,
    )
    probability = classifier.predict_proba(validation[:, selected])[:, 1]
    macro = _macro_metrics(
        _accession_metrics(
            validation_metadata,
            validation_labels,
            probability,
            str(candidate["training_arm"]),
        )
    )
    return {**candidate, **{f"validation_{key}": value for key, value in macro.items()}}


def _select_best(table: pd.DataFrame) -> pd.Series:
    ordered = table.sort_values(
        [
            "validation_balanced_accuracy",
            "validation_roc_auc",
            "validation_average_precision",
            "feature_count",
            "synthetic_weight",
            "rank_method",
            "training_arm",
        ],
        ascending=[False, False, False, True, True, True, True],
        kind="stable",
    )
    return ordered.iloc[0]


def _symbol_mapping(path: Path) -> dict[str, str]:
    with h5py.File(path, "r") as handle:
        genes = _decode(handle["meta/genes/ensembl_gene"][:])
        symbols = _decode(handle["meta/genes/symbol"][:])
    return dict(zip(genes, symbols))


def _adjust_bh(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def _reactome_enrichment(
    selected_genes: list[str],
    background_genes: list[str],
    gmt_path: Path,
    symbols: dict[str, str] | None = None,
) -> pd.DataFrame:
    selected = set(selected_genes)
    background = set(background_genes)
    rows: list[dict[str, object]] = []
    with gmt_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            term, description, *genes = fields
            pathway = set(genes) & background
            overlap = pathway & selected
            if len(pathway) < 3 or len(overlap) < 2:
                continue
            rows.append(
                {
                    "term": term,
                    "description": description,
                    "pathway_genes_in_background": len(pathway),
                    "selected_genes": len(selected),
                    "overlap": len(overlap),
                    "overlap_genes": ",".join(sorted(overlap)),
                    "overlap_symbols": ",".join(
                        symbols.get(gene, gene) if symbols else gene
                        for gene in sorted(overlap)
                    ),
                    "p_value": float(
                        hypergeom.sf(
                            len(overlap) - 1,
                            len(background),
                            len(pathway),
                            len(selected),
                        )
                    ),
                }
            )
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    table["fdr"] = _adjust_bh(table["p_value"])
    return table.sort_values(["fdr", "p_value", "term"])


def _plot_performance(folds: pd.DataFrame, output: Path) -> None:
    positions = np.arange(len(folds))
    width = 0.25
    figure, axis = plt.subplots(figsize=(10, 5.5))
    axis.bar(
        positions - width,
        folds["baseline_balanced_accuracy"],
        width,
        label="Real-only baseline",
        color="#4C78A8",
    )
    axis.bar(
        positions,
        folds["generated_balanced_accuracy"],
        width,
        label="Best generated candidate",
        color="#F58518",
    )
    axis.bar(
        positions + width,
        folds["deployed_balanced_accuracy"],
        width,
        label="Inner-gated deployment",
        color="#54A24B",
    )
    axis.axhline(0.5, color="#333333", linestyle="--", linewidth=1)
    axis.set_xticks(positions, folds["fold"])
    axis.set_ylim(0, 1.03)
    axis.set_ylabel("Outer accession-macro balanced accuracy")
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output / "outer_fold_balanced_accuracy.png", dpi=220)
    figure.savefig(output / "outer_fold_balanced_accuracy.pdf")
    plt.close(figure)


def _plot_features(features: pd.DataFrame, output: Path, *, tissue: str) -> None:
    if features.empty:
        return
    top = features.head(25).copy()
    labels = top["symbol"].where(top["symbol"].ne(""), top["gene"])
    figure, axis = plt.subplots(figsize=(9, 8))
    if int(features["selected_folds"].max()) == 1:
        top = top.sort_values("mean_absolute_classifier_coefficient")
        labels = top["symbol"].where(top["symbol"].ne(""), top["gene"])
        values = top["mean_classifier_coefficient"]
        axis.barh(
            np.arange(len(top)),
            values,
            color=np.where(values >= 0, "#D64E4E", "#3977B8"),
        )
        axis.axvline(0, color="#333333", linewidth=0.8)
        axis.set_xlabel("Standardized logistic coefficient (positive = FLT)")
    else:
        top = top.sort_values(
            ["selected_folds", "mean_absolute_classifier_coefficient"]
        )
        labels = top["symbol"].where(top["symbol"].ne(""), top["gene"])
        axis.barh(
            np.arange(len(top)),
            top["selected_folds"],
            color=np.where(
                top["mean_classifier_coefficient"] >= 0, "#D64E4E", "#3977B8"
            ),
        )
        axis.set_xlabel("Selected outer folds")
        axis.set_xlim(0, int(top["selected_folds"].max()) + 1)
    axis.set_yticks(np.arange(len(top)), labels)
    axis.set_title(f"Top DDIM-guided {tissue.replace('_', ' ')} features")
    figure.tight_layout()
    figure.savefig(output / "stable_feature_importance.png", dpi=220)
    figure.savefig(output / "stable_feature_importance.pdf")
    plt.close(figure)


def _evaluate_fold(
    fold: dict[str, object],
    config: dict[str, object],
    symbols: dict[str, str],
    *,
    seed: int,
) -> dict[str, object]:
    label = str(fold["label"])
    output = Path(str(config["output_dir"])) / "folds" / label
    output.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(str(fold["samples_tsv"]), sep="\t")
    with h5py.File(str(fold["prepared_h5"]), "r") as handle:
        genes = _decode(handle["genes"][:])
        train_expression, train_metadata = _load_role(
            handle, samples, "train", str(config["tissue"])
        )
        validation_expression, validation_metadata = _load_role(
            handle, samples, "validation", str(config["tissue"])
        )
        test_expression, test_metadata = _load_role(
            handle, samples, "test", str(config["tissue"])
        )
    train_labels = _labels(train_metadata)
    validation_labels = _labels(validation_metadata)
    test_labels = _labels(test_metadata)
    scaler = StandardScaler().fit(train_expression)
    train = scaler.transform(train_expression)
    validation = scaler.transform(validation_expression)
    test = scaler.transform(test_expression)

    source_rows = train_metadata["_row_index"].to_numpy(dtype=np.int64)
    synthetic_unscaled = [
        _load_generated_draw(
            _generated_path(Path(str(fold["generated_dir"])), int(draw_seed)),
            source_rows,
        )
        for draw_seed in fold["seeds"]
    ]
    synthetic_draws = [scaler.transform(draw) for draw in synthetic_unscaled]
    synthetic_recentered_draws = [
        scaler.transform(_recenter_draw(draw, train_expression, train_labels))
        for draw in synthetic_unscaled
    ]
    synthetic_raw = np.concatenate(synthetic_draws)
    synthetic_recentered = np.concatenate(synthetic_recentered_draws)
    rankings, diagnostics = _build_rankings(
        train,
        train_labels,
        train_metadata,
        synthetic_draws,
        seed=seed,
    )
    diagnostics.insert(0, "symbol", [symbols.get(gene, "") for gene in genes])
    diagnostics.insert(0, "gene", genes)
    for name, score in rankings.items():
        diagnostics[f"ranking_{name}"] = score
    diagnostics.to_csv(output / "train_feature_diagnostics.tsv.gz", sep="\t", index=False)

    grid = config["grid"]
    baseline_rows: list[dict[str, object]] = []
    for feature_count in grid["feature_counts"]:
        for regularization_c in grid["regularization_c"]:
            candidate = {
                "rank_method": "real_f",
                "training_arm": "real_only",
                "synthetic_weight": 0.0,
                "feature_count": int(feature_count),
                "regularization_c": float(regularization_c),
            }
            baseline_rows.append(
                _candidate_metrics(
                    candidate,
                    train,
                    train_labels,
                    validation,
                    validation_labels,
                    validation_metadata,
                    synthetic_raw,
                    synthetic_recentered,
                    rankings,
                    seed=seed,
                )
            )
    baseline_table = pd.DataFrame(baseline_rows)
    best_baseline = _select_best(baseline_table)

    generated_rows: list[dict[str, object]] = []
    for rank_method in RANK_METHODS:
        for feature_count in grid["feature_counts"]:
            for regularization_c in grid["regularization_c"]:
                arms = [("real_only", 0.0)] + [
                    (arm, float(weight))
                    for arm in ("raw_synthetic", "recentered_synthetic")
                    for weight in grid["synthetic_weights"]
                ]
                for training_arm, synthetic_weight in arms:
                    candidate = {
                        "rank_method": rank_method,
                        "training_arm": training_arm,
                        "synthetic_weight": synthetic_weight,
                        "feature_count": int(feature_count),
                        "regularization_c": float(regularization_c),
                    }
                    generated_rows.append(
                        _candidate_metrics(
                            candidate,
                            train,
                            train_labels,
                            validation,
                            validation_labels,
                            validation_metadata,
                            synthetic_raw,
                            synthetic_recentered,
                            rankings,
                            seed=seed,
                        )
                    )
    generated_table = pd.DataFrame(generated_rows)
    best_generated = _select_best(generated_table)
    gate = bool(
        best_generated["validation_balanced_accuracy"]
        >= best_baseline["validation_balanced_accuracy"]
        + float(grid["minimum_validation_ba_gain"])
        and best_generated["validation_roc_auc"]
        >= best_baseline["validation_roc_auc"]
        - float(grid["maximum_validation_auc_loss"])
    )
    baseline_table.assign(role="baseline").to_csv(
        output / "inner_baseline_candidates.tsv", sep="\t", index=False
    )
    generated_table.assign(role="generated").to_csv(
        output / "inner_generated_candidates.tsv.gz", sep="\t", index=False
    )

    predictions: dict[str, np.ndarray] = {}
    models: dict[str, tuple[LogisticRegression, np.ndarray]] = {}
    for arm, candidate in (
        ("baseline", best_baseline),
        ("generated", best_generated),
    ):
        classifier, selected = _fit_candidate(
            candidate,
            train,
            train_labels,
            synthetic_raw,
            synthetic_recentered,
            rankings,
            seed=seed,
        )
        models[arm] = classifier, selected
        predictions[arm] = classifier.predict_proba(test[:, selected])[:, 1]
    predictions["deployed"] = predictions["generated"] if gate else predictions["baseline"]

    accession_tables = []
    metrics: dict[str, dict[str, float]] = {}
    for arm, probability in predictions.items():
        table = _accession_metrics(
            test_metadata, test_labels, probability, arm
        )
        accession_tables.append(table)
        metrics[arm] = _macro_metrics(table)
    pd.concat(accession_tables, ignore_index=True).to_csv(
        output / "outer_accession_metrics.tsv", sep="\t", index=False
    )
    prediction_table = test_metadata[
        ["_row_index", "profile_id", "accession", "condition", "tissue"]
    ].copy()
    prediction_table["label"] = test_labels
    for arm, probability in predictions.items():
        prediction_table[f"{arm}_probability"] = probability
    prediction_table.to_csv(output / "outer_predictions.tsv.gz", sep="\t", index=False)

    generated_classifier, generated_indices = models["generated"]
    feature_table = diagnostics.iloc[generated_indices].copy()
    feature_table.insert(0, "feature_rank", np.arange(1, len(feature_table) + 1))
    feature_table.insert(0, "fold", label)
    feature_table["classifier_coefficient"] = generated_classifier.coef_[0]
    feature_table["inner_gate_passed"] = gate
    feature_table["rank_method"] = str(best_generated["rank_method"])
    feature_table["training_arm"] = str(best_generated["training_arm"])
    feature_table.to_csv(output / "selected_generated_features.tsv", sep="\t", index=False)

    selected = {
        "fold": label,
        "test_accessions": sorted(test_metadata["accession"].astype(str).unique()),
        "train_profiles": int(len(train)),
        "validation_profiles": int(len(validation)),
        "test_profiles": int(len(test)),
        "baseline": best_baseline.to_dict(),
        "generated": best_generated.to_dict(),
        "inner_gate_passed": gate,
        "outer_accession_macro": metrics,
    }
    (output / "summary.json").write_text(
        json.dumps(selected, indent=2) + "\n", encoding="utf-8"
    )
    return selected


def _aggregate(config: dict[str, object], symbols: dict[str, str]) -> dict[str, object]:
    output = Path(str(config["output_dir"]))
    fold_rows: list[dict[str, object]] = []
    accession_tables: list[pd.DataFrame] = []
    prediction_tables: list[pd.DataFrame] = []
    feature_tables: list[pd.DataFrame] = []
    for fold in config["folds"]:
        label = str(fold["label"])
        fold_output = output / "folds" / label
        summary = json.loads((fold_output / "summary.json").read_text())
        row: dict[str, object] = {
            "fold": label,
            "test_accessions": ",".join(summary["test_accessions"]),
            "test_profiles": summary["test_profiles"],
            "inner_gate_passed": summary["inner_gate_passed"],
            "generated_rank_method": summary["generated"]["rank_method"],
            "generated_training_arm": summary["generated"]["training_arm"],
            "generated_feature_count": summary["generated"]["feature_count"],
            "generated_synthetic_weight": summary["generated"]["synthetic_weight"],
        }
        for arm in ("baseline", "generated", "deployed"):
            for metric, value in summary["outer_accession_macro"][arm].items():
                row[f"{arm}_{metric}"] = value
        fold_rows.append(row)
        accession = pd.read_csv(fold_output / "outer_accession_metrics.tsv", sep="\t")
        accession.insert(0, "fold", label)
        accession_tables.append(accession)
        prediction = pd.read_csv(fold_output / "outer_predictions.tsv.gz", sep="\t")
        prediction.insert(0, "fold", label)
        prediction_tables.append(prediction)
        feature_tables.append(
            pd.read_csv(fold_output / "selected_generated_features.tsv", sep="\t")
        )
    folds = pd.DataFrame(fold_rows)
    accessions = pd.concat(accession_tables, ignore_index=True)
    predictions = pd.concat(prediction_tables, ignore_index=True)
    selected_features = pd.concat(feature_tables, ignore_index=True)
    folds.to_csv(output / "outer_fold_results.tsv", sep="\t", index=False)
    accessions.to_csv(output / "outer_accession_results.tsv", sep="\t", index=False)
    predictions.to_csv(output / "outer_predictions.tsv.gz", sep="\t", index=False)
    selected_features.to_csv(output / "selected_features_all_folds.tsv.gz", sep="\t", index=False)

    feature_stability = (
        selected_features.groupby(["gene", "symbol"], observed=True)
        .agg(
            selected_folds=("fold", "nunique"),
            gated_folds=("inner_gate_passed", "sum"),
            mean_classifier_coefficient=("classifier_coefficient", "mean"),
            mean_absolute_classifier_coefficient=(
                "classifier_coefficient",
                lambda values: float(np.mean(np.abs(values))),
            ),
            coefficient_sign_agreement=(
                "classifier_coefficient",
                lambda values: float(abs(np.mean(np.sign(values)))),
            ),
            mean_real_effect=("real_effect", "mean"),
            mean_synthetic_effect=("synthetic_effect", "mean"),
            effect_direction_match_rate=("effect_direction_match", "mean"),
            mean_real_effect_sign_stability=("real_effect_sign_stability", "mean"),
            mean_synthetic_effect_sign_stability=(
                "synthetic_effect_sign_stability",
                "mean",
            ),
        )
        .reset_index()
        .sort_values(
            [
                "gated_folds",
                "selected_folds",
                "coefficient_sign_agreement",
                "mean_absolute_classifier_coefficient",
            ],
            ascending=[False, False, False, False],
        )
    )
    feature_stability.to_csv(output / "feature_stability.tsv", sep="\t", index=False)

    with h5py.File(str(config["folds"][0]["prepared_h5"]), "r") as handle:
        background = _decode(handle["genes"][:]).tolist()
    top_genes = feature_stability.head(min(100, len(feature_stability)))["gene"].tolist()
    enrichment = _reactome_enrichment(
        top_genes,
        background,
        Path(str(config["annotations"]["reactome_gmt"])),
        symbols,
    )
    enrichment.to_csv(output / "reactome_enrichment.tsv", sep="\t", index=False)

    macro: dict[str, dict[str, float]] = {}
    for arm in ("baseline", "generated", "deployed"):
        subset = accessions.loc[accessions["arm"].eq(arm)]
        macro[arm] = _macro_metrics(subset)
    generated_delta = {
        metric: macro["generated"][metric] - macro["baseline"][metric]
        for metric in METRICS
    }
    deployed_delta = {
        metric: macro["deployed"][metric] - macro["baseline"][metric]
        for metric in METRICS
    }
    promising = bool(
        folds["inner_gate_passed"].sum() >= 2
        and deployed_delta["balanced_accuracy"] > 0
        and deployed_delta["roc_auc"] >= 0
        and deployed_delta["average_precision"] >= 0
    )
    summary = {
        "status": "complete",
        "design": "post_hoc_cross_fitted_outer_accession_evaluation",
        "tissue": str(config["tissue"]),
        "outer_folds": int(len(folds)),
        "outer_accessions": int(
            accessions.loc[accessions["arm"].eq("baseline"), "accession"].nunique()
        ),
        "outer_profiles": int(len(predictions)),
        "inner_gate_passed_folds": int(folds["inner_gate_passed"].sum()),
        "accession_macro": macro,
        "generated_minus_baseline": generated_delta,
        "deployed_minus_baseline": deployed_delta,
        "promising_by_predeclared_rule": promising,
        "top_feature_genes": feature_stability.head(20)[
            ["gene", "symbol", "selected_folds", "gated_folds"]
        ].to_dict(orient="records"),
        "significant_reactome_terms_fdr_0_05": (
            int((enrichment["fdr"] < 0.05).sum()) if not enrichment.empty else 0
        ),
        "interpretation": (
            "Exploratory only because the analyst had seen outer-fold outcomes before "
            "this cross-fitted protocol was specified."
        ),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    _plot_performance(folds, output)
    _plot_features(feature_stability, output, tissue=str(config["tissue"]))
    return summary


def run(config_path: Path, *, seed: int = 4040) -> dict[str, object]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    output = Path(str(config["output_dir"]))
    output.mkdir(parents=True, exist_ok=True)
    symbols = _symbol_mapping(Path(str(config["annotations"]["archs4_h5"])))
    for offset, fold in enumerate(config["folds"]):
        print(f"[generated-features] evaluating {fold['label']}", flush=True)
        _evaluate_fold(fold, config, symbols, seed=seed + offset)
    summary = _aggregate(config, symbols)
    print(json.dumps(summary, indent=2))
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--seed", default=4040, type=int)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    run(arguments.config, seed=arguments.seed)


if __name__ == "__main__":
    main()
