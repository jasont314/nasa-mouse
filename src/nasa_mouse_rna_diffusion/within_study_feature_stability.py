"""Repeated within-study DDIM-guided FLT/GC feature stability workflow."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import textwrap
from typing import Any, Iterable
import warnings

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    roc_auc_score,
)
from sklearn.preprocessing import StandardScaler

from nasa_mouse_generative.effect_validation import (
    accession_effects,
    leave_one_accession_out,
    random_effects_table,
)

from .factorized_adapter import load_factorized_role
from .factorized_calibrate import _aligned_expression
from .factorized_distribution_calibrate import PositiveResidualCalibrator
from .generated_feature_guidance import (
    _build_rankings,
    _fit_classifier,
    _reactome_enrichment,
    _recenter_draw,
    _selected_indices,
    _symbol_mapping,
)


METRICS = ("balanced_accuracy", "roc_auc", "average_precision")
GENERATED_ARMS = (
    "generated_only",
    "real_plus_generated",
    "guided_real_only",
    "guided_low_weight",
)


@dataclass(frozen=True)
class WorkflowData:
    genes: list[str]
    symbols: dict[str, str]
    development_expression: np.ndarray
    development_samples: pd.DataFrame
    test_expression: np.ndarray
    test_samples: pd.DataFrame
    all_expression: np.ndarray
    all_samples: pd.DataFrame
    synthetic_draws: dict[str, np.ndarray]


def _build_rankings_quiet(*args: Any, **kwargs: Any) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    """Build feature rankings without sklearn's expected constant-feature noise."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Features .* are constant\.",
            category=UserWarning,
        )
        return _build_rankings(*args, **kwargs)


def _labels(samples: pd.DataFrame) -> np.ndarray:
    condition = samples["condition"].astype(str)
    unexpected = sorted(set(condition) - {"flight", "ground_control"})
    if unexpected:
        raise ValueError(f"Unexpected conditions: {unexpected}")
    return condition.eq("flight").to_numpy(dtype=np.int64)


def _metric_set(labels: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    if np.unique(labels).size != 2:
        return {name: float("nan") for name in METRICS}
    return {
        "balanced_accuracy": float(
            balanced_accuracy_score(labels, probability >= 0.5)
        ),
        "roc_auc": float(roc_auc_score(labels, probability)),
        "average_precision": float(average_precision_score(labels, probability)),
    }


def _within_stratum_split(
    samples: pd.DataFrame,
    *,
    fraction: float,
    seed: int,
    strata: Iterable[str] = ("accession", "condition"),
) -> tuple[np.ndarray, np.ndarray]:
    """Split samples inside accession/condition strata while retaining training rows."""

    if not 0.0 < float(fraction) < 0.5:
        raise ValueError("Split fraction must be in (0, 0.5)")
    columns = tuple(map(str, strata))
    missing = set(columns).difference(samples.columns)
    if missing:
        raise ValueError(f"Split metadata lacks columns: {sorted(missing)}")
    rng = np.random.default_rng(int(seed))
    held_out: list[int] = []
    grouped = samples.groupby(list(columns), sort=True, observed=True, dropna=False)
    for _, frame in grouped:
        positions = frame.index.to_numpy(dtype=np.int64)
        if len(positions) < 2:
            continue
        shuffled = rng.permutation(positions)
        count = max(1, int(round(len(positions) * float(fraction))))
        count = min(count, len(positions) - 1)
        held_out.extend(map(int, shuffled[:count]))
    held_out_array = np.asarray(sorted(set(held_out)), dtype=np.int64)
    retained = np.setdiff1d(
        samples.index.to_numpy(dtype=np.int64), held_out_array, assume_unique=False
    )
    return retained, held_out_array


def _align_generated(
    role: dict[str, object], path: Path
) -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    arrays = np.load(path)
    source_rows = np.asarray(arrays["source_row"], dtype=np.int64)
    real, samples = _aligned_expression(role, source_rows)
    synthetic = np.asarray(arrays["scaled_expression"], dtype=np.float32)
    if synthetic.shape != real.shape:
        raise ValueError(f"Generated and real matrices differ for {path}")
    return synthetic, samples, source_rows


def _ordered_by_source_row(
    matrix: np.ndarray,
    source_rows: np.ndarray,
    expected_rows: np.ndarray,
) -> np.ndarray:
    lookup = {int(row): index for index, row in enumerate(source_rows)}
    missing = sorted(set(map(int, expected_rows)) - set(lookup))
    if missing:
        raise ValueError(f"Generated matrix is missing source rows: {missing[:5]}")
    indices = np.asarray([lookup[int(row)] for row in expected_rows], dtype=np.int64)
    return np.asarray(matrix[indices], dtype=np.float32)


def _load_data(config: dict[str, Any]) -> WorkflowData:
    options = config["data"]
    prepared_h5 = Path(options["prepared_h5"])
    samples_tsv = Path(options["samples_tsv"])
    roles = {
        role: load_factorized_role(prepared_h5, samples_tsv, role)
        for role in ("train", "validation", "test")
    }
    genes = list(map(str, roles["train"]["genes"]))
    if any(list(map(str, role["genes"])) != genes for role in roles.values()):
        raise ValueError("Prepared role gene orders differ")

    development_expression = np.concatenate(
        [roles["train"]["expression"], roles["validation"]["expression"]]
    ).astype(np.float32)
    development_samples = pd.concat(
        [roles["train"]["samples"], roles["validation"]["samples"]],
        ignore_index=True,
    )
    test_expression = np.asarray(roles["test"]["expression"], dtype=np.float32)
    test_samples = roles["test"]["samples"].reset_index(drop=True)
    all_expression = np.concatenate((development_expression, test_expression))
    all_samples = pd.concat((development_samples, test_samples), ignore_index=True)
    expected_rows = development_samples["_row_index"].to_numpy(dtype=np.int64)

    calibrator = PositiveResidualCalibrator.load(options["calibrator"])
    synthetic_draws: dict[str, np.ndarray] = {}
    for draw in options["synthetic_draws"]:
        raw_parts: list[np.ndarray] = []
        sample_parts: list[pd.DataFrame] = []
        source_parts: list[np.ndarray] = []
        for role_name in ("train", "validation"):
            synthetic, samples, source_rows = _align_generated(
                roles[role_name], Path(draw[f"{role_name}_path"])
            )
            raw_parts.append(synthetic)
            sample_parts.append(samples)
            source_parts.append(source_rows)
        raw = np.concatenate(raw_parts)
        metadata = pd.concat(sample_parts, ignore_index=True)
        source_rows = np.concatenate(source_parts)
        calibrated = calibrator.apply(
            raw, metadata, seed=int(draw["residual_seed"])
        )
        synthetic_draws[str(draw["name"])] = _ordered_by_source_row(
            calibrated, source_rows, expected_rows
        )

    symbols = _symbol_mapping(Path(config["annotations"]["archs4_h5"]))
    return WorkflowData(
        genes=genes,
        symbols=symbols,
        development_expression=development_expression,
        development_samples=development_samples,
        test_expression=test_expression,
        test_samples=test_samples,
        all_expression=all_expression,
        all_samples=all_samples,
        synthetic_draws=synthetic_draws,
    )


def _arm_specs(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    defaults: dict[str, dict[str, Any]] = {
        "real_only": {
            "rank_methods": ["real_f"],
            "training": "real",
            "synthetic_weight": 0.0,
            "recenter": False,
        },
        "generated_only": {
            "rank_methods": ["synthetic_f"],
            "training": "synthetic",
            "synthetic_weight": 1.0,
            "recenter": False,
        },
        "real_plus_generated": {
            "rank_methods": ["f_rank_consensus"],
            "training": "combined",
            "synthetic_weight": 1.0,
            "recenter": False,
        },
        "guided_real_only": {
            "rank_methods": ["f_rank_consensus", "effect_consensus"],
            "training": "real",
            "synthetic_weight": 0.0,
            "recenter": False,
        },
        "guided_low_weight": {
            "rank_methods": ["f_rank_consensus", "effect_consensus"],
            "training": "combined",
            "synthetic_weight": 0.05,
            "recenter": True,
        },
    }
    configured = config.get("arms", {})
    for name, values in configured.items():
        if name not in defaults:
            raise ValueError(f"Unknown workflow arm: {name}")
        defaults[name].update(values or {})
    return defaults


def _fit_arm(
    arm: dict[str, Any],
    *,
    rank_method: str,
    feature_count: int,
    regularization_c: float,
    real: np.ndarray,
    labels: np.ndarray,
    synthetic_draws: list[np.ndarray],
    recentered_draws: list[np.ndarray],
    rankings: dict[str, np.ndarray],
    seed: int,
) -> tuple[LogisticRegression, np.ndarray]:
    selected = _selected_indices(rankings[rank_method], feature_count)
    training = str(arm["training"])
    if training == "real":
        return (
            _fit_classifier(
                real[:, selected],
                labels,
                regularization_c=regularization_c,
                seed=seed,
            ),
            selected,
        )

    selected_draws = recentered_draws if bool(arm["recenter"]) else synthetic_draws
    synthetic = np.concatenate([draw[:, selected] for draw in selected_draws])
    synthetic_labels = np.tile(labels, len(selected_draws))
    if training == "synthetic":
        return (
            _fit_classifier(
                synthetic,
                synthetic_labels,
                regularization_c=regularization_c,
                seed=seed,
            ),
            selected,
        )
    if training != "combined":
        raise ValueError(f"Unsupported arm training mode: {training}")
    weight = float(arm["synthetic_weight"]) / len(selected_draws)
    expression = np.concatenate((real[:, selected], synthetic))
    combined_labels = np.concatenate((labels, synthetic_labels))
    sample_weight = np.concatenate(
        (
            np.ones(len(real), dtype=float),
            np.full(len(synthetic), weight, dtype=float),
        )
    )
    return (
        _fit_classifier(
            expression,
            combined_labels,
            regularization_c=regularization_c,
            seed=seed,
            sample_weight=sample_weight,
        ),
        selected,
    )


def _scaled_views(
    real: np.ndarray,
    labels: np.ndarray,
    synthetic: list[np.ndarray],
    scaler: StandardScaler,
) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray]]:
    scaled_real = scaler.transform(real)
    scaled_synthetic = [scaler.transform(draw) for draw in synthetic]
    recentered = [
        scaler.transform(_recenter_draw(draw, real, labels)) for draw in synthetic
    ]
    return scaled_real, scaled_synthetic, recentered


def _select_candidate(table: pd.DataFrame) -> pd.Series:
    return table.sort_values(
        [
            "balanced_accuracy",
            "roc_auc",
            "average_precision",
            "feature_count",
            "regularization_c",
            "rank_method",
        ],
        ascending=[False, False, False, True, True, True],
        kind="stable",
    ).iloc[0]


def _candidate_grid(
    arm_name: str,
    arm: dict[str, Any],
    *,
    real_train: np.ndarray,
    train_labels: np.ndarray,
    train_metadata: pd.DataFrame,
    synthetic_train: list[np.ndarray],
    real_validation: np.ndarray,
    validation_labels: np.ndarray,
    grid: dict[str, Any],
    seed: int,
) -> tuple[pd.Series, pd.DataFrame]:
    scaler = StandardScaler().fit(real_train)
    real_scaled, synthetic_scaled, recentered_scaled = _scaled_views(
        real_train, train_labels, synthetic_train, scaler
    )
    rankings, _ = _build_rankings_quiet(
        real_scaled,
        train_labels,
        train_metadata.reset_index(drop=True),
        synthetic_scaled,
        seed=seed,
    )
    validation_scaled = scaler.transform(real_validation)
    rows: list[dict[str, object]] = []
    for rank_method in arm["rank_methods"]:
        for feature_count in grid["feature_counts"]:
            for regularization_c in grid["regularization_c"]:
                classifier, selected = _fit_arm(
                    arm,
                    rank_method=str(rank_method),
                    feature_count=int(feature_count),
                    regularization_c=float(regularization_c),
                    real=real_scaled,
                    labels=train_labels,
                    synthetic_draws=synthetic_scaled,
                    recentered_draws=recentered_scaled,
                    rankings=rankings,
                    seed=seed,
                )
                probability = classifier.predict_proba(
                    validation_scaled[:, selected]
                )[:, 1]
                rows.append(
                    {
                        "arm": arm_name,
                        "rank_method": str(rank_method),
                        "feature_count": int(feature_count),
                        "regularization_c": float(regularization_c),
                        **_metric_set(validation_labels, probability),
                    }
                )
    table = pd.DataFrame(rows)
    return _select_candidate(table), table


def _fit_selected(
    arm_name: str,
    arm: dict[str, Any],
    candidate: pd.Series | dict[str, object],
    *,
    real_train: np.ndarray,
    train_labels: np.ndarray,
    train_metadata: pd.DataFrame,
    synthetic_train: list[np.ndarray],
    real_evaluation: np.ndarray,
    evaluation_labels: np.ndarray,
    genes: list[str],
    symbols: dict[str, str],
    seed: int,
) -> tuple[dict[str, float], pd.DataFrame, np.ndarray]:
    scaler = StandardScaler().fit(real_train)
    real_scaled, synthetic_scaled, recentered_scaled = _scaled_views(
        real_train, train_labels, synthetic_train, scaler
    )
    rankings, diagnostics = _build_rankings_quiet(
        real_scaled,
        train_labels,
        train_metadata.reset_index(drop=True),
        synthetic_scaled,
        seed=seed,
    )
    classifier, selected = _fit_arm(
        arm,
        rank_method=str(candidate["rank_method"]),
        feature_count=int(candidate["feature_count"]),
        regularization_c=float(candidate["regularization_c"]),
        real=real_scaled,
        labels=train_labels,
        synthetic_draws=synthetic_scaled,
        recentered_draws=recentered_scaled,
        rankings=rankings,
        seed=seed,
    )
    probability = classifier.predict_proba(
        scaler.transform(real_evaluation)[:, selected]
    )[:, 1]
    features = diagnostics.iloc[selected].copy()
    features.insert(0, "symbol", [symbols.get(genes[index], "") for index in selected])
    features.insert(0, "gene", [genes[index] for index in selected])
    features["classifier_coefficient"] = classifier.coef_[0]
    features["arm"] = arm_name
    features["rank_method"] = str(candidate["rank_method"])
    features["feature_count"] = int(candidate["feature_count"])
    features["regularization_c"] = float(candidate["regularization_c"])
    return _metric_set(evaluation_labels, probability), features, probability


def _valid_nested_split(
    samples: pd.DataFrame,
    *,
    outer_fraction: float,
    inner_fraction: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    outer_train, outer_test = _within_stratum_split(
        samples, fraction=outer_fraction, seed=seed
    )
    if len(outer_test) < 4 or np.unique(_labels(samples.loc[outer_test])).size != 2:
        return None
    inner_source = samples.loc[outer_train].copy()
    inner_source.index = outer_train
    inner_train, inner_validation = _within_stratum_split(
        inner_source, fraction=inner_fraction, seed=seed + 10_000
    )
    if (
        len(inner_validation) < 4
        or np.unique(_labels(samples.loc[inner_train])).size != 2
        or np.unique(_labels(samples.loc[inner_validation])).size != 2
    ):
        return None
    return inner_train, inner_validation, outer_test


def _run_tissue(
    tissue: str,
    data: WorkflowData,
    config: dict[str, Any],
    output: Path,
    *,
    repeats: int,
    seed: int,
) -> dict[str, Any]:
    development_mask = data.development_samples["tissue"].astype(str).eq(tissue)
    development_indices = np.flatnonzero(development_mask.to_numpy())
    real = data.development_expression[development_indices]
    samples = data.development_samples.loc[development_mask].reset_index(drop=True)
    synthetic = [
        draw[development_indices] for draw in data.synthetic_draws.values()
    ]
    labels = _labels(samples)
    minimum = int(config["analysis"].get("minimum_development_profiles", 12))
    minimum_class = int(config["analysis"].get("minimum_profiles_per_condition", 5))
    if len(samples) < minimum or min(np.bincount(labels, minlength=2)) < minimum_class:
        return {
            "status": "insufficient_development_profiles",
            "tissue": tissue,
            "development_profiles": int(len(samples)),
            "flight": int(labels.sum()),
            "ground_control": int((labels == 0).sum()),
        }

    arms = _arm_specs(config)
    metric_rows: list[dict[str, object]] = []
    feature_tables: list[pd.DataFrame] = []
    candidate_tables: list[pd.DataFrame] = []
    completed = 0
    for repeat in range(int(repeats)):
        split = _valid_nested_split(
            samples,
            outer_fraction=float(config["analysis"]["outer_fraction"]),
            inner_fraction=float(config["analysis"]["inner_fraction"]),
            seed=seed + repeat,
        )
        if split is None:
            continue
        inner_train, inner_validation, outer_test = split
        outer_train = np.setdiff1d(np.arange(len(samples)), outer_test)
        for arm_offset, (arm_name, arm) in enumerate(arms.items()):
            selected, candidates = _candidate_grid(
                arm_name,
                arm,
                real_train=real[inner_train],
                train_labels=labels[inner_train],
                train_metadata=samples.loc[inner_train],
                synthetic_train=[draw[inner_train] for draw in synthetic],
                real_validation=real[inner_validation],
                validation_labels=labels[inner_validation],
                grid=config["grid"],
                seed=seed + repeat * 100 + arm_offset,
            )
            candidates.insert(0, "repeat", repeat)
            candidates.insert(0, "tissue", tissue)
            candidate_tables.append(candidates)
            metrics, features, _ = _fit_selected(
                arm_name,
                arm,
                selected,
                real_train=real[outer_train],
                train_labels=labels[outer_train],
                train_metadata=samples.loc[outer_train],
                synthetic_train=[draw[outer_train] for draw in synthetic],
                real_evaluation=real[outer_test],
                evaluation_labels=labels[outer_test],
                genes=data.genes,
                symbols=data.symbols,
                seed=seed + repeat * 100 + arm_offset,
            )
            metric_rows.append(
                {
                    "tissue": tissue,
                    "repeat": repeat,
                    "arm": arm_name,
                    "outer_train_profiles": int(len(outer_train)),
                    "outer_test_profiles": int(len(outer_test)),
                    "rank_method": str(selected["rank_method"]),
                    "feature_count": int(selected["feature_count"]),
                    "regularization_c": float(selected["regularization_c"]),
                    **metrics,
                }
            )
            features.insert(0, "repeat", repeat)
            features.insert(0, "tissue", tissue)
            feature_tables.append(features)
        completed += 1

    tissue_output = output / tissue
    tissue_output.mkdir(parents=True, exist_ok=True)
    if not metric_rows:
        return {
            "status": "insufficient_nested_splits",
            "tissue": tissue,
            "development_profiles": int(len(samples)),
        }
    metrics = pd.DataFrame(metric_rows)
    selected_features = pd.concat(feature_tables, ignore_index=True)
    candidates = pd.concat(candidate_tables, ignore_index=True)
    metrics.to_csv(tissue_output / "nested_arm_metrics.tsv", sep="\t", index=False)
    selected_features.to_csv(
        tissue_output / "nested_selected_features.tsv.gz",
        sep="\t",
        index=False,
    )
    candidates.to_csv(
        tissue_output / "inner_candidate_metrics.tsv.gz", sep="\t", index=False
    )
    return {
        "status": "complete",
        "tissue": tissue,
        "development_profiles": int(len(samples)),
        "flight": int(labels.sum()),
        "ground_control": int((labels == 0).sum()),
        "completed_repeats": int(completed),
        "metrics": metrics,
        "selected_features": selected_features,
    }


def _arm_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    return (
        metrics.groupby(["tissue", "arm"], observed=True)
        .agg(
            repeats=("repeat", "nunique"),
            mean_balanced_accuracy=("balanced_accuracy", "mean"),
            sd_balanced_accuracy=("balanced_accuracy", "std"),
            mean_roc_auc=("roc_auc", "mean"),
            sd_roc_auc=("roc_auc", "std"),
            mean_average_precision=("average_precision", "mean"),
            sd_average_precision=("average_precision", "std"),
        )
        .reset_index()
    )


def _choose_arms(
    arm_summary: pd.DataFrame, prior: pd.DataFrame
) -> pd.DataFrame:
    prior_lookup = prior.set_index("tissue").to_dict(orient="index") if not prior.empty else {}
    rows: list[dict[str, object]] = []
    for tissue, frame in arm_summary.groupby("tissue", sort=True):
        baseline = frame.loc[frame["arm"].eq("real_only")].iloc[0]
        generated = frame.loc[frame["arm"].isin(GENERATED_ARMS)].copy()
        eligible = generated.loc[
            generated["mean_balanced_accuracy"].ge(baseline["mean_balanced_accuracy"])
            & generated["mean_roc_auc"].ge(baseline["mean_roc_auc"])
            & generated["mean_average_precision"].ge(
                baseline["mean_average_precision"]
            )
        ]
        selected = (
            eligible.sort_values(
                [
                    "mean_balanced_accuracy",
                    "mean_roc_auc",
                    "mean_average_precision",
                    "arm",
                ],
                ascending=[False, False, False, True],
                kind="stable",
            ).iloc[0]
            if not eligible.empty
            else baseline
        )
        best_generated = generated.sort_values(
            [
                "mean_balanced_accuracy",
                "mean_roc_auc",
                "mean_average_precision",
                "arm",
            ],
            ascending=[False, False, False, True],
            kind="stable",
        ).iloc[0]
        prior_row = prior_lookup.get(str(tissue), {})
        prior_mode = str(prior_row.get("best_mode", ""))
        prior_arm = {
            "generated": "generated_only",
            "real+generated": "real_plus_generated",
        }.get(prior_mode, prior_mode)
        rows.append(
            {
                "tissue": str(tissue),
                "selected_arm": str(selected["arm"]),
                "best_generated_arm": str(best_generated["arm"]),
                "prior_best_use_arm": prior_arm,
                "prior_comparison_status": prior_row.get("comparison_status", ""),
                "prior_test_n": prior_row.get("test_n", np.nan),
                "real_mean_balanced_accuracy": float(
                    baseline["mean_balanced_accuracy"]
                ),
                "selected_mean_balanced_accuracy": float(
                    selected["mean_balanced_accuracy"]
                ),
                "real_mean_roc_auc": float(baseline["mean_roc_auc"]),
                "selected_mean_roc_auc": float(selected["mean_roc_auc"]),
                "real_mean_average_precision": float(
                    baseline["mean_average_precision"]
                ),
                "selected_mean_average_precision": float(
                    selected["mean_average_precision"]
                ),
                "generated_arm_eligible_all_metrics": bool(
                    str(selected["arm"]) != "real_only"
                ),
            }
        )
    return pd.DataFrame(rows)


def _paired_repeat_support(
    metrics: pd.DataFrame, choices: pd.DataFrame
) -> pd.DataFrame:
    """Summarize selected-versus-real differences on matched outer splits."""

    rows: list[dict[str, object]] = []
    for choice in choices.itertuples(index=False):
        tissue = str(choice.tissue)
        arm = str(choice.selected_arm)
        baseline = metrics.loc[
            metrics["tissue"].eq(tissue) & metrics["arm"].eq("real_only"),
            ["repeat", *METRICS],
        ]
        selected = metrics.loc[
            metrics["tissue"].eq(tissue) & metrics["arm"].eq(arm),
            ["repeat", *METRICS],
        ]
        paired = baseline.merge(
            selected,
            on="repeat",
            how="inner",
            suffixes=("_real", "_selected"),
            validate="one_to_one",
        )
        record: dict[str, object] = {
            "tissue": tissue,
            "selected_arm": arm,
            "paired_repeats": int(len(paired)),
        }
        deltas: list[np.ndarray] = []
        for metric in METRICS:
            delta = (
                paired[f"{metric}_selected"] - paired[f"{metric}_real"]
            ).to_numpy(dtype=float)
            deltas.append(delta)
            record[f"mean_delta_{metric}"] = float(np.mean(delta))
            record[f"nonworse_rate_{metric}"] = float(np.mean(delta >= 0.0))
            record[f"strict_win_rate_{metric}"] = float(np.mean(delta > 0.0))
        matrix = np.column_stack(deltas)
        record["all_metrics_nonworse_rate"] = float(
            np.mean(np.all(matrix >= 0.0, axis=1))
        )
        record["all_metrics_strict_win_rate"] = float(
            np.mean(np.all(matrix > 0.0, axis=1))
        )
        rows.append(record)
    return pd.DataFrame(rows)


def _feature_stability(
    selected_features: pd.DataFrame, repeats_by_tissue: dict[str, int]
) -> pd.DataFrame:
    table = (
        selected_features.groupby(["tissue", "arm", "gene", "symbol"], observed=True)
        .agg(
            selected_runs=("repeat", "nunique"),
            median_classifier_coefficient=("classifier_coefficient", "median"),
            mean_absolute_classifier_coefficient=(
                "classifier_coefficient",
                lambda values: float(np.mean(np.abs(values))),
            ),
            coefficient_sign_agreement=(
                "classifier_coefficient",
                lambda values: float(abs(np.mean(np.sign(values)))),
            ),
            median_real_effect=("real_effect", "median"),
            median_synthetic_effect=("synthetic_effect", "median"),
            effect_direction_match_rate=("effect_direction_match", "mean"),
        )
        .reset_index()
    )
    table["completed_repeats"] = table["tissue"].map(repeats_by_tissue)
    table["selection_frequency"] = (
        table["selected_runs"] / table["completed_repeats"]
    )
    return table.sort_values(
        [
            "tissue",
            "arm",
            "selection_frequency",
            "coefficient_sign_agreement",
            "mean_absolute_classifier_coefficient",
        ],
        ascending=[True, True, False, False, False],
        kind="stable",
    )


def _real_effect_summary(
    tissue: str, data: WorkflowData
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = data.all_samples["tissue"].astype(str).eq(tissue)
    samples = data.all_samples.loc[mask].reset_index(drop=True)
    expression = data.all_expression[np.flatnonzero(mask.to_numpy())]
    effects = accession_effects(expression, samples, data.genes)
    meta = random_effects_table(effects)
    if meta.empty:
        return effects, pd.DataFrame({"gene": data.genes})
    _, loo = leave_one_accession_out(effects, meta)
    summary = meta.rename(columns={"feature": "gene"})
    if not loo.empty:
        summary = summary.merge(
            loo.rename(columns={"feature": "gene"}), on="gene", how="left"
        )
    summary["accession_direction_fraction"] = (
        summary["n_accession_same_direction"] / summary["n_accessions"]
    )
    summary["loo_direction_stable"] = summary.get(
        "n_same_direction", pd.Series(0, index=summary.index)
    ).eq(summary.get("n_leave_one_out", pd.Series(-1, index=summary.index)))
    summary["loo_fdr_stable_0_05"] = (
        summary["meta_fdr"].lt(0.05)
        & summary.get(
            "maximum_leave_one_out_fdr", pd.Series(np.inf, index=summary.index)
        ).lt(0.05)
        & summary["loo_direction_stable"]
    )
    return effects, summary


def _stable_gene_sets(
    tissue: str,
    stability: pd.DataFrame,
    choices: pd.DataFrame,
    effects: pd.DataFrame,
    *,
    minimum_frequency: float,
    minimum_sign_agreement: float,
    minimum_accession_direction: float,
) -> pd.DataFrame:
    choice = choices.loc[choices["tissue"].eq(tissue)].iloc[0]
    best_generated = str(choice["best_generated_arm"])
    subset = stability.loc[stability["tissue"].eq(tissue)].copy()
    stable = subset.loc[
        subset["selection_frequency"].ge(minimum_frequency)
        & subset["coefficient_sign_agreement"].ge(minimum_sign_agreement)
    ]
    real = stable.loc[stable["arm"].eq("real_only")].set_index("gene")
    generated = stable.loc[stable["arm"].eq(best_generated)].set_index("gene")
    effect_lookup = effects.set_index("gene") if not effects.empty else pd.DataFrame()
    genes = sorted(set(real.index) | set(generated.index))
    rows: list[dict[str, object]] = []
    for gene in genes:
        real_row = real.loc[gene] if gene in real.index else None
        generated_row = generated.loc[gene] if gene in generated.index else None
        effect = (
            effect_lookup.loc[gene]
            if gene in effect_lookup.index and "meta_effect" in effect_lookup.columns
            else None
        )
        real_coefficient = (
            float(real_row["median_classifier_coefficient"])
            if real_row is not None
            else float("nan")
        )
        generated_coefficient = (
            float(generated_row["median_classifier_coefficient"])
            if generated_row is not None
            else float("nan")
        )
        meta_effect = float(effect["meta_effect"]) if effect is not None else float("nan")
        direction_fraction = (
            float(effect["accession_direction_fraction"])
            if effect is not None
            else float("nan")
        )
        coefficient_match = bool(
            real_row is not None
            and generated_row is not None
            and np.sign(real_coefficient) == np.sign(generated_coefficient)
        )
        real_effect_supports_generated = bool(
            generated_row is not None
            and np.isfinite(meta_effect)
            and np.sign(generated_coefficient) == np.sign(meta_effect)
            and direction_fraction >= minimum_accession_direction
        )
        if coefficient_match and real_effect_supports_generated:
            gene_set = "core_intersection"
        elif generated_row is not None and real_effect_supports_generated:
            gene_set = "generated_supported"
        else:
            gene_set = "exploratory_union"
        symbol = ""
        if real_row is not None:
            symbol = str(real_row["symbol"])
        elif generated_row is not None:
            symbol = str(generated_row["symbol"])
        rows.append(
            {
                "tissue": tissue,
                "gene": gene,
                "symbol": symbol,
                "gene_set": gene_set,
                "best_generated_arm": best_generated,
                "stable_real": real_row is not None,
                "stable_generated": generated_row is not None,
                "real_selection_frequency": (
                    float(real_row["selection_frequency"])
                    if real_row is not None
                    else 0.0
                ),
                "generated_selection_frequency": (
                    float(generated_row["selection_frequency"])
                    if generated_row is not None
                    else 0.0
                ),
                "real_coefficient": real_coefficient,
                "generated_coefficient": generated_coefficient,
                "coefficient_direction_match": coefficient_match,
                "real_meta_effect": meta_effect,
                "real_meta_fdr": (
                    float(effect["meta_fdr"]) if effect is not None else float("nan")
                ),
                "real_accession_direction_fraction": direction_fraction,
                "real_loo_fdr_stable_0_05": (
                    bool(effect["loo_fdr_stable_0_05"])
                    if effect is not None
                    else False
                ),
                "real_effect_supports_generated": real_effect_supports_generated,
            }
        )
    return pd.DataFrame(rows)


def _enrichment_for_sets(
    tissue: str,
    genes: pd.DataFrame,
    background: list[str],
    config: dict[str, Any],
    symbols: dict[str, str],
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    set_members = {
        "core_intersection": genes.loc[
            genes["gene_set"].eq("core_intersection"), "gene"
        ].tolist(),
        "generated_supported": genes.loc[
            genes["gene_set"].isin(["core_intersection", "generated_supported"]),
            "gene",
        ].tolist(),
        "exploratory_union": genes["gene"].tolist(),
    }
    for name, selected in set_members.items():
        table = _reactome_enrichment(
            selected,
            background,
            Path(config["annotations"]["reactome_gmt"]),
            symbols,
        )
        if table.empty:
            continue
        table.insert(0, "gene_set", name)
        table.insert(0, "tissue", tissue)
        rows.append(table)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _biological_support_summary(
    choices: pd.DataFrame,
    gene_sets: pd.DataFrame,
    real_effects: pd.DataFrame,
    enrichment: pd.DataFrame,
) -> pd.DataFrame:
    """Count independent real-data and pathway support for each tissue screen."""

    rows: list[dict[str, object]] = []
    for choice in choices.itertuples(index=False):
        tissue = str(choice.tissue)
        genes = gene_sets.loc[gene_sets["tissue"].eq(tissue)]
        effects = real_effects.loc[real_effects["tissue"].eq(tissue)]
        pathways = enrichment.loc[enrichment["tissue"].eq(tissue)]
        record: dict[str, object] = {
            "tissue": tissue,
            "selected_arm": str(choice.selected_arm),
            "stable_genes": int(len(genes)),
            "core_intersection_genes": int(
                genes["gene_set"].eq("core_intersection").sum()
            ),
            "generated_supported_genes": int(
                genes["gene_set"].eq("generated_supported").sum()
            ),
            "exploratory_union_genes": int(
                genes["gene_set"].eq("exploratory_union").sum()
            ),
            "stable_genes_real_fdr_0_05": int(genes["real_meta_fdr"].lt(0.05).sum()),
            "stable_genes_real_loo_fdr_0_05": int(
                genes["real_loo_fdr_stable_0_05"].fillna(False).sum()
            ),
            "all_panel_genes_real_fdr_0_05": int(effects["meta_fdr"].lt(0.05).sum()),
            "all_panel_genes_real_loo_fdr_0_05": int(
                effects["loo_fdr_stable_0_05"].fillna(False).sum()
            ),
        }
        for gene_set in (
            "core_intersection",
            "generated_supported",
            "exploratory_union",
        ):
            record[f"reactome_fdr_0_05_{gene_set}"] = int(
                (
                    pathways["gene_set"].eq(gene_set)
                    & pathways["fdr"].lt(0.05)
                ).sum()
            )
        rows.append(record)
    return pd.DataFrame(rows)


def _modal_candidate(metrics: pd.DataFrame, tissue: str, arm: str) -> dict[str, object]:
    subset = metrics.loc[metrics["tissue"].eq(tissue) & metrics["arm"].eq(arm)]
    counts = (
        subset.groupby(["rank_method", "feature_count", "regularization_c"])
        .size()
        .rename("count")
        .reset_index()
        .sort_values(
            ["count", "feature_count", "regularization_c", "rank_method"],
            ascending=[False, True, True, True],
            kind="stable",
        )
    )
    return counts.iloc[0].to_dict()


def _descriptive_test(
    data: WorkflowData,
    config: dict[str, Any],
    metrics: pd.DataFrame,
    choices: pd.DataFrame,
    tissues: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    arms = _arm_specs(config)
    rows: list[dict[str, object]] = []
    features: list[pd.DataFrame] = []
    for offset, tissue in enumerate(tissues):
        development_mask = data.development_samples["tissue"].astype(str).eq(tissue)
        test_mask = data.test_samples["tissue"].astype(str).eq(tissue)
        if test_mask.sum() < 4:
            continue
        train_labels = _labels(data.development_samples.loc[development_mask])
        test_labels = _labels(data.test_samples.loc[test_mask])
        if np.unique(test_labels).size != 2:
            continue
        train_indices = np.flatnonzero(development_mask.to_numpy())
        choice = choices.loc[choices["tissue"].eq(tissue)].iloc[0]
        for role, arm_name in (
            ("baseline", "real_only"),
            ("selected", str(choice["selected_arm"])),
        ):
            candidate = _modal_candidate(metrics, tissue, arm_name)
            result, selected, probability = _fit_selected(
                arm_name,
                arms[arm_name],
                candidate,
                real_train=data.development_expression[train_indices],
                train_labels=train_labels,
                train_metadata=data.development_samples.loc[
                    development_mask
                ].reset_index(drop=True),
                synthetic_train=[
                    draw[train_indices] for draw in data.synthetic_draws.values()
                ],
                real_evaluation=data.test_expression[np.flatnonzero(test_mask.to_numpy())],
                evaluation_labels=test_labels,
                genes=data.genes,
                symbols=data.symbols,
                seed=int(config["run"]["seed"]) + offset,
            )
            rows.append(
                {
                    "tissue": tissue,
                    "role": role,
                    "arm": arm_name,
                    "test_profiles": int(test_mask.sum()),
                    **result,
                }
            )
            selected.insert(0, "role", role)
            selected.insert(0, "tissue", tissue)
            features.append(selected)
    return pd.DataFrame(rows), pd.concat(features, ignore_index=True)


def _plot_arm_summary(summary: pd.DataFrame, output: Path) -> None:
    tissues = sorted(summary["tissue"].unique())
    arms = ["real_only", *GENERATED_ARMS]
    pivot = summary.pivot(index="tissue", columns="arm", values="mean_balanced_accuracy")
    matrix = pivot.reindex(index=tissues, columns=arms).to_numpy(dtype=float)
    figure, axis = plt.subplots(figsize=(10, max(6, 0.35 * len(tissues))))
    image = axis.imshow(matrix, aspect="auto", vmin=0.3, vmax=1.0, cmap="viridis")
    axis.set_xticks(np.arange(len(arms)), [name.replace("_", " ") for name in arms], rotation=30, ha="right")
    axis.set_yticks(np.arange(len(tissues)), [name.replace("_", " ") for name in tissues])
    axis.set_title("Nested within-study FLT/GC balanced accuracy")
    figure.colorbar(image, ax=axis, label="Balanced accuracy")
    figure.tight_layout()
    figure.savefig(output / "arm_balanced_accuracy_heatmap.png", dpi=220)
    figure.savefig(output / "arm_balanced_accuracy_heatmap.pdf")
    plt.close(figure)


def _plot_gene_set_counts(genes: pd.DataFrame, output: Path) -> None:
    if genes.empty:
        return
    counts = genes.groupby(["tissue", "gene_set"]).size().unstack(fill_value=0)
    order = ["core_intersection", "generated_supported", "exploratory_union"]
    counts = counts.reindex(columns=order, fill_value=0)
    figure, axis = plt.subplots(figsize=(11, max(6, 0.35 * len(counts))))
    left = np.zeros(len(counts))
    colors = ["#2A9D8F", "#E9C46A", "#A8A8A8"]
    for column, color in zip(order, colors):
        axis.barh(np.arange(len(counts)), counts[column], left=left, label=column.replace("_", " "), color=color)
        left += counts[column].to_numpy()
    axis.set_yticks(np.arange(len(counts)), [value.replace("_", " ") for value in counts.index])
    axis.set_xlabel("Stable genes")
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output / "stable_gene_set_counts.png", dpi=220)
    figure.savefig(output / "stable_gene_set_counts.pdf")
    plt.close(figure)


def _plot_tissue_gene_support(genes: pd.DataFrame, output: Path) -> None:
    colors = {
        "core_intersection": "#2A9D8F",
        "generated_supported": "#E9C46A",
        "exploratory_union": "#A8A8A8",
    }
    for tissue, frame in genes.groupby("tissue", sort=True, observed=True):
        table = frame.copy()
        table["display_coefficient"] = table["generated_coefficient"].where(
            table["generated_coefficient"].notna(), table["real_coefficient"]
        )
        table = table.loc[table["display_coefficient"].notna()].copy()
        table["strength"] = table["display_coefficient"].abs()
        loo_stable = table.loc[table["real_loo_fdr_stable_0_05"].fillna(False)]
        remaining = table.drop(index=loo_stable.index).nlargest(
            max(0, 25 - len(loo_stable)), "strength"
        )
        table = pd.concat((loo_stable, remaining)).sort_values(
            "display_coefficient"
        )
        if table.empty:
            continue
        labels = []
        for row in table.itertuples(index=False):
            label = str(row.symbol) if str(row.symbol) else str(row.gene)
            if bool(row.real_loo_fdr_stable_0_05):
                label += " *"
            labels.append(label)
        figure, axis = plt.subplots(figsize=(9, max(5.5, 0.3 * len(table))))
        positions = np.arange(len(table))
        axis.barh(
            positions,
            table["display_coefficient"],
            color=[colors[str(value)] for value in table["gene_set"]],
        )
        axis.axvline(0.0, color="#333333", linewidth=0.8)
        axis.set_yticks(positions, labels)
        axis.set_xlabel("Median standardized classifier coefficient (positive = FLT)")
        axis.set_title(f"{str(tissue).replace('_', ' ')} stable FLT/GC genes")
        handles = [
            plt.Rectangle((0, 0), 1, 1, color=color, label=name.replace("_", " "))
            for name, color in colors.items()
        ]
        axis.legend(handles=handles, frameon=False, fontsize=8, loc="best")
        axis.text(
            0.99,
            0.01,
            "* real-data LOO FDR < 0.05",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
        )
        figure.tight_layout()
        tissue_output = output / str(tissue)
        figure.savefig(tissue_output / "stable_gene_coefficients.png", dpi=220)
        figure.savefig(tissue_output / "stable_gene_coefficients.pdf")
        plt.close(figure)


def _plot_tissue_pathways(enrichment: pd.DataFrame, output: Path) -> None:
    for tissue, frame in enrichment.groupby("tissue", sort=True, observed=True):
        preferred = frame.loc[frame["gene_set"].eq("generated_supported")]
        if preferred.empty:
            preferred = frame.loc[frame["gene_set"].eq("core_intersection")]
        if preferred.empty:
            preferred = frame.loc[frame["gene_set"].eq("exploratory_union")]
        table = (
            preferred.sort_values(["fdr", "p_value", "term"], kind="stable")
            .drop_duplicates("term")
            .head(15)
            .copy()
        )
        if table.empty:
            continue
        table["score"] = -np.log10(table["fdr"].clip(lower=1e-300))
        table = table.sort_values("score")
        labels = []
        for term in table["term"].astype(str):
            name = term.split("_", 1)[1] if "_" in term else term
            labels.append(textwrap.fill(name.replace("_", " "), width=48))
        figure, axis = plt.subplots(figsize=(11, max(5.5, 0.48 * len(table))))
        positions = np.arange(len(table))
        axis.barh(
            positions,
            table["score"],
            color=np.where(table["fdr"].lt(0.05), "#2A9D8F", "#A8A8A8"),
        )
        axis.axvline(-np.log10(0.05), color="#333333", linestyle="--", linewidth=1)
        axis.set_yticks(positions, labels)
        axis.set_xlabel("-log10 Reactome FDR")
        gene_set = str(table["gene_set"].iloc[0]).replace("_", " ")
        axis.set_title(
            f"{str(tissue).replace('_', ' ')} Reactome enrichment ({gene_set})"
        )
        figure.tight_layout()
        tissue_output = output / str(tissue)
        figure.savefig(tissue_output / "reactome_enrichment.png", dpi=220)
        figure.savefig(tissue_output / "reactome_enrichment.pdf")
        plt.close(figure)


def _write_readme(output: Path, summary: dict[str, Any]) -> None:
    improved = summary["selected_generated_arms"]
    text = f"""# Within-study generated-feature stability

This analysis reuses the fixed ARCHS4-pretrained, OSDR-fine-tuned DDIM and does
not perform additional neural-network training. The original train and validation
roles are pooled for repeated nested development splits; the original test role is
never used to select an arm or a feature. Because accessions are represented on
both sides of each split, performance measures within-study interpolation.

The earlier best-use table is retained as a prior/sanity-check column. Arm choices
are re-estimated from nested development folds using balanced accuracy, AUROC, and
average precision independently. No composite score is used.

Completed tissues: {summary['completed_tissues']}

Tissues selecting a generated-informed arm: {', '.join(improved) if improved else 'none'}

## Outputs

- `arm_summary.tsv`: repeated nested metrics by tissue and arm.
- `tissue_arm_choices.tsv`: selected arm and earlier best-use result.
- `paired_repeat_support.tsv`: selected-versus-real differences on matched splits.
- `feature_stability.tsv.gz`: selection frequency and coefficient stability.
- `stable_gene_sets.tsv.gz`: core, generated-supported, and exploratory genes.
- `real_accession_effects.tsv.gz`: all-data within-accession FLT-GC effects.
- `real_random_effects.tsv.gz`: random-effects FDR and LOO summaries.
- `biological_support_summary.tsv`: real-gene and pathway support counts by tissue.
- `reactome_enrichment.tsv.gz`: enrichment using the 974-gene panel as background.
- `<tissue>/stable_gene_coefficients.png`: top stable signed gene coefficients.
- `<tissue>/reactome_enrichment.png`: top pathway enrichments and FDR threshold.
- `descriptive_original_test.tsv`: post-selection descriptive test metrics; not a
  fresh confirmatory result because this test role had been examined previously.

Generated profiles are model outputs, not independent biological replicates.
Gene claims therefore require support from real accession-level effects. LOO
statistics are sensitivity checks and do not imply the DDIM generalized to an
unseen accession. The DDIM was fine-tuned before these nested splits and saw the
original training role, so repeated-split gains are exploratory and can be
optimistic. Repeated splits overlap and are not independent observations.
"""
    (output / "README.md").write_text(text, encoding="utf-8")


def run(
    config_path: Path,
    *,
    tissues_override: list[str] | None = None,
    repeats_override: int | None = None,
    output_override: Path | None = None,
) -> Path:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    output = output_override or Path(config["run"]["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    data = _load_data(config)
    tissues = tissues_override or config["analysis"].get("tissues") or sorted(
        data.development_samples["tissue"].astype(str).unique()
    )
    tissues = list(map(str, tissues))
    repeats = int(repeats_override or config["analysis"]["repeats"])
    seed = int(config["run"]["seed"])

    tissue_results: list[dict[str, Any]] = []
    metric_tables: list[pd.DataFrame] = []
    feature_tables: list[pd.DataFrame] = []
    for offset, tissue in enumerate(tissues):
        print(f"[within-study-features] {tissue}", flush=True)
        result = _run_tissue(
            tissue,
            data,
            config,
            output,
            repeats=repeats,
            seed=seed + offset * 100_000,
        )
        tissue_results.append(
            {key: value for key, value in result.items() if key not in {"metrics", "selected_features"}}
        )
        if result["status"] == "complete":
            metric_tables.append(result["metrics"])
            feature_tables.append(result["selected_features"])

    inventory = pd.DataFrame(tissue_results)
    inventory.to_csv(output / "tissue_inventory.tsv", sep="\t", index=False)
    if not metric_tables:
        raise RuntimeError("No tissue completed the nested workflow")
    metrics = pd.concat(metric_tables, ignore_index=True)
    selected_features = pd.concat(feature_tables, ignore_index=True)
    metrics.to_csv(output / "nested_arm_metrics.tsv.gz", sep="\t", index=False)
    selected_features.to_csv(
        output / "nested_selected_features.tsv.gz", sep="\t", index=False
    )

    arm_summary = _arm_summary(metrics)
    prior_path = Path(config["prior"]["best_use_table"])
    prior = pd.read_csv(prior_path, sep="\t") if prior_path.exists() else pd.DataFrame()
    choices = _choose_arms(arm_summary, prior)
    paired_support = _paired_repeat_support(metrics, choices)
    completed = sorted(choices["tissue"].astype(str).unique())
    repeats_by_tissue = (
        metrics.groupby("tissue")["repeat"].nunique().astype(int).to_dict()
    )
    stability = _feature_stability(selected_features, repeats_by_tissue)
    arm_summary.to_csv(output / "arm_summary.tsv", sep="\t", index=False)
    choices.to_csv(output / "tissue_arm_choices.tsv", sep="\t", index=False)
    paired_support.to_csv(
        output / "paired_repeat_support.tsv", sep="\t", index=False
    )
    stability.to_csv(output / "feature_stability.tsv.gz", sep="\t", index=False)

    effect_tables: list[pd.DataFrame] = []
    effect_summaries: list[pd.DataFrame] = []
    gene_sets: list[pd.DataFrame] = []
    enrichments: list[pd.DataFrame] = []
    thresholds = config["stability"]
    for tissue in completed:
        effects, effect_summary = _real_effect_summary(tissue, data)
        if not effects.empty:
            effects.insert(0, "analysis_tissue", tissue)
            effect_tables.append(effects)
        if not effect_summary.empty:
            effect_summary.insert(0, "tissue", tissue)
            effect_summaries.append(effect_summary)
        sets = _stable_gene_sets(
            tissue,
            stability,
            choices,
            effect_summary,
            minimum_frequency=float(thresholds["minimum_selection_frequency"]),
            minimum_sign_agreement=float(thresholds["minimum_coefficient_sign_agreement"]),
            minimum_accession_direction=float(thresholds["minimum_accession_direction_fraction"]),
        )
        if not sets.empty:
            gene_sets.append(sets)
            enrichment = _enrichment_for_sets(
                tissue, sets, data.genes, config, data.symbols
            )
            if not enrichment.empty:
                enrichments.append(enrichment)

    all_effects = pd.concat(effect_tables, ignore_index=True) if effect_tables else pd.DataFrame()
    all_effect_summaries = pd.concat(effect_summaries, ignore_index=True) if effect_summaries else pd.DataFrame()
    all_gene_sets = pd.concat(gene_sets, ignore_index=True) if gene_sets else pd.DataFrame()
    all_enrichment = pd.concat(enrichments, ignore_index=True) if enrichments else pd.DataFrame()
    all_effects.to_csv(output / "real_accession_effects.tsv.gz", sep="\t", index=False)
    all_effect_summaries.to_csv(output / "real_random_effects.tsv.gz", sep="\t", index=False)
    all_gene_sets.to_csv(output / "stable_gene_sets.tsv.gz", sep="\t", index=False)
    all_enrichment.to_csv(output / "reactome_enrichment.tsv.gz", sep="\t", index=False)
    biological_support = _biological_support_summary(
        choices, all_gene_sets, all_effect_summaries, all_enrichment
    )
    biological_support.to_csv(
        output / "biological_support_summary.tsv", sep="\t", index=False
    )

    descriptive_test, final_features = _descriptive_test(
        data, config, metrics, choices, completed
    )
    descriptive_test.to_csv(
        output / "descriptive_original_test.tsv", sep="\t", index=False
    )
    final_features.to_csv(
        output / "final_development_selected_features.tsv.gz", sep="\t", index=False
    )
    _plot_arm_summary(arm_summary, output)
    _plot_gene_set_counts(all_gene_sets, output)
    _plot_tissue_gene_support(all_gene_sets, output)
    _plot_tissue_pathways(all_enrichment, output)

    selected_generated = choices.loc[
        choices["selected_arm"].ne("real_only"), "tissue"
    ].astype(str).tolist()
    summary = {
        "status": "complete",
        "design": "repeated_nested_within_accession_feature_stability",
        "config": str(config_path.resolve()),
        "fixed_ddim": config["data"].get("model", ""),
        "neural_network_retrained": False,
        "development_roles": ["train", "validation"],
        "arm_selection_uses_original_test": False,
        "descriptive_test_previously_opened": True,
        "synthetic_draws": list(data.synthetic_draws),
        "requested_tissues": len(tissues),
        "completed_tissues": len(completed),
        "selected_generated_arms": selected_generated,
        "selected_arm_counts": choices["selected_arm"].value_counts().to_dict(),
        "selected_arms_nonworse_all_metrics_in_at_least_75pct_repeats": (
            paired_support.loc[
                paired_support["selected_arm"].ne("real_only")
                & paired_support["all_metrics_nonworse_rate"].ge(0.75),
                "tissue",
            ].astype(str).tolist()
        ),
        "gene_set_counts": (
            all_gene_sets["gene_set"].value_counts().to_dict()
            if not all_gene_sets.empty
            else {}
        ),
        "significant_reactome_terms_fdr_0_05": (
            int(all_enrichment["fdr"].lt(0.05).sum())
            if not all_enrichment.empty
            else 0
        ),
        "limitations": [
            "Splits retain accessions on both sides and measure within-study interpolation.",
            "The fixed DDIM saw the original training role before the nested splits.",
            "Repeated outer splits overlap and are not independent observations.",
            "The original test role was previously examined; its metrics are descriptive.",
            "Generated profiles are not independent biological replicates.",
            "Random-effects and LOO summaries use real data for biological support.",
        ],
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    _write_readme(output, summary)
    print(json.dumps(summary, indent=2), flush=True)
    return output / "summary.json"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--tissues", nargs="+")
    parser.add_argument("--repeats", type=int)
    parser.add_argument("--output", type=Path)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    run(
        arguments.config,
        tissues_override=arguments.tissues,
        repeats_override=arguments.repeats,
        output_override=arguments.output,
    )


if __name__ == "__main__":
    main()
