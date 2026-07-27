"""Accession-aware real-versus-synthetic condition-effect validation."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


def _bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    clean = np.where(np.isfinite(values), values, 1.0)
    order = np.argsort(clean, kind="stable")
    ranked = clean[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def _random_effects(effect: np.ndarray, variance: np.ndarray) -> dict[str, float | int]:
    valid = np.isfinite(effect) & np.isfinite(variance) & (variance > 0)
    values = np.asarray(effect, dtype=float)[valid]
    variances = np.asarray(variance, dtype=float)[valid]
    if len(values) < 2:
        return {
            "n_accessions": int(len(values)),
            "meta_effect": float("nan"),
            "meta_se": float("nan"),
            "meta_p": float("nan"),
            "tau2": float("nan"),
            "i2": float("nan"),
        }
    fixed_weights = 1.0 / variances
    fixed_effect = float(np.sum(fixed_weights * values) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * np.square(values - fixed_effect)))
    degrees_freedom = len(values) - 1
    c_value = float(
        np.sum(fixed_weights)
        - np.sum(np.square(fixed_weights)) / np.sum(fixed_weights)
    )
    tau2 = max((q - degrees_freedom) / c_value, 0.0) if c_value > 0 else 0.0
    weights = 1.0 / (variances + tau2)
    meta_effect = float(np.sum(weights * values) / np.sum(weights))
    meta_se = float(np.sqrt(1.0 / np.sum(weights)))
    z_score = meta_effect / meta_se if meta_se > 0 else float("nan")
    meta_p = (
        float(2.0 * stats.norm.sf(abs(z_score)))
        if np.isfinite(z_score)
        else float("nan")
    )
    i2 = max((q - degrees_freedom) / q, 0.0) if q > 0 else 0.0
    return {
        "n_accessions": int(len(values)),
        "meta_effect": meta_effect,
        "meta_se": meta_se,
        "meta_p": meta_p,
        "tau2": float(tau2),
        "i2": float(i2),
    }


def accession_effects(
    matrix: np.ndarray,
    samples: pd.DataFrame,
    feature_names: list[str],
) -> pd.DataFrame:
    """Compute FLT-minus-GC effects and sampling variances within accessions."""

    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (len(samples), len(feature_names)):
        raise ValueError("Expression, sample metadata, and feature names do not align")
    required = {"accession", "tissue", "condition"}
    missing = required.difference(samples.columns)
    if missing:
        raise ValueError(f"Samples lack accession-effect columns: {sorted(missing)}")
    rows: list[pd.DataFrame] = []
    groups = samples.groupby(["accession", "tissue"], sort=True, dropna=False)
    for (accession, tissue), frame in groups:
        positions = frame.index.to_numpy(dtype=int)
        conditions = frame["condition"].astype(str).to_numpy()
        flight = positions[conditions == "flight"]
        ground = positions[conditions == "ground_control"]
        if len(flight) < 2 or len(ground) < 2:
            continue
        flight_values = matrix[flight]
        ground_values = matrix[ground]
        effect = flight_values.mean(axis=0) - ground_values.mean(axis=0)
        variance = (
            flight_values.var(axis=0, ddof=1) / len(flight)
            + ground_values.var(axis=0, ddof=1) / len(ground)
        )
        rows.append(
            pd.DataFrame(
                {
                    "accession": str(accession),
                    "tissue": str(tissue),
                    "feature": feature_names,
                    "n_flight": len(flight),
                    "n_ground_control": len(ground),
                    "flight_minus_ground": effect,
                    "effect_variance": variance,
                }
            )
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "accession",
                "tissue",
                "feature",
                "n_flight",
                "n_ground_control",
                "flight_minus_ground",
                "effect_variance",
            ]
        )
    return pd.concat(rows, ignore_index=True)


def random_effects_table(effects: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for feature, frame in effects.groupby("feature", sort=False):
        result = _random_effects(
            frame["flight_minus_ground"].to_numpy(),
            frame["effect_variance"].to_numpy(),
        )
        direction = np.sign(float(result["meta_effect"]))
        observed = np.sign(frame["flight_minus_ground"].to_numpy(dtype=float))
        rows.append(
            {
                "feature": feature,
                **result,
                "n_accession_same_direction": int(np.sum(observed == direction)),
                "n_accession_opposite_direction": int(
                    np.sum(observed == -direction)
                ),
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["meta_fdr"] = _bh_fdr(result["meta_p"].to_numpy())
    return result.sort_values(["meta_fdr", "meta_p", "feature"], kind="stable")


def leave_one_accession_out(
    effects: pd.DataFrame,
    primary: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    accessions = sorted(effects["accession"].astype(str).unique())
    for held_out in accessions:
        retained = effects.loc[effects["accession"].astype(str).ne(held_out)]
        for feature, frame in retained.groupby("feature", sort=False):
            rows.append(
                {
                    "held_out_accession": held_out,
                    "feature": feature,
                    **_random_effects(
                        frame["flight_minus_ground"].to_numpy(),
                        frame["effect_variance"].to_numpy(),
                    ),
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result, pd.DataFrame()
    result["meta_fdr"] = result.groupby("held_out_accession", sort=False)[
        "meta_p"
    ].transform(lambda values: _bh_fdr(values.to_numpy()))
    primary_effect = primary.set_index("feature")["meta_effect"]
    summary_rows: list[dict[str, Any]] = []
    for feature, frame in result.groupby("feature", sort=False):
        direction = np.sign(float(primary_effect.get(feature, 0.0)))
        loo_direction = np.sign(frame["meta_effect"].to_numpy(dtype=float))
        summary_rows.append(
            {
                "feature": feature,
                "n_leave_one_out": int(len(frame)),
                "n_same_direction": int(np.sum(loo_direction == direction)),
                "minimum_leave_one_out_fdr": float(frame["meta_fdr"].min()),
                "maximum_leave_one_out_fdr": float(frame["meta_fdr"].max()),
            }
        )
    return result, pd.DataFrame(summary_rows)


def compare_real_synthetic_effects(
    real: np.ndarray,
    synthetic: np.ndarray,
    samples: pd.DataFrame,
    feature_names: list[str],
) -> tuple[dict[str, pd.DataFrame], dict[str, Any]]:
    """Run accession-aware meta-analysis independently and compare effect recovery."""

    local_samples = samples.reset_index(drop=True)
    real_effects = accession_effects(real, local_samples, feature_names)
    synthetic_effects = accession_effects(synthetic, local_samples, feature_names)
    real_meta = random_effects_table(real_effects)
    synthetic_meta = random_effects_table(synthetic_effects)
    real_loo, real_loo_summary = leave_one_accession_out(real_effects, real_meta)
    synthetic_loo, synthetic_loo_summary = leave_one_accession_out(
        synthetic_effects, synthetic_meta
    )
    if real_meta.empty or synthetic_meta.empty:
        comparison = pd.DataFrame({"feature": feature_names})
        return {
            "real_per_accession": real_effects,
            "synthetic_per_accession": synthetic_effects,
            "real_meta": real_meta,
            "synthetic_meta": synthetic_meta,
            "real_leave_one_out": real_loo,
            "synthetic_leave_one_out": synthetic_loo,
            "comparison": comparison,
        }, {
            "features": int(len(feature_names)),
            "accessions": int(real_effects["accession"].nunique()),
            "status": "insufficient_accessions_with_both_conditions",
            "meta_effect_correlation": float("nan"),
            "meta_direction_agreement": float("nan"),
            "real_random_effects_fdr_lt_005": 0,
            "synthetic_random_effects_fdr_lt_005": 0,
            "real_loo_stable_fdr_lt_005": 0,
            "synthetic_loo_stable_fdr_lt_005": 0,
            "concordant_loo_stable_fdr_lt_005": 0,
        }
    comparison = real_meta.add_prefix("real_").rename(
        columns={"real_feature": "feature"}
    ).merge(
        synthetic_meta.add_prefix("synthetic_").rename(
            columns={"synthetic_feature": "feature"}
        ),
        on="feature",
        how="outer",
    )
    comparison["direction_agrees"] = np.sign(
        comparison["real_meta_effect"]
    ).eq(np.sign(comparison["synthetic_meta_effect"]))
    stability_columns = [
        "feature",
        "n_leave_one_out",
        "n_same_direction",
        "minimum_leave_one_out_fdr",
        "maximum_leave_one_out_fdr",
    ]

    def prefixed_stability(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
        normalized = frame.reindex(columns=stability_columns)
        return normalized.add_prefix(f"{prefix}_").rename(
            columns={f"{prefix}_feature": "feature"}
        )

    real_stability = prefixed_stability(real_loo_summary, "real")
    synthetic_stability = prefixed_stability(
        synthetic_loo_summary, "synthetic"
    )
    stability = real_stability.merge(
        synthetic_stability, on="feature", how="outer"
    )
    comparison = comparison.merge(stability, on="feature", how="left")
    real_effect = comparison["real_meta_effect"].to_numpy(dtype=float)
    synthetic_effect = comparison["synthetic_meta_effect"].to_numpy(dtype=float)
    valid = np.isfinite(real_effect) & np.isfinite(synthetic_effect)
    correlation = (
        float(np.corrcoef(real_effect[valid], synthetic_effect[valid])[0, 1])
        if valid.sum() >= 2
        and np.std(real_effect[valid]) > 0
        and np.std(synthetic_effect[valid]) > 0
        else float("nan")
    )
    real_stable = (
        comparison["real_meta_fdr"].lt(0.05)
        & comparison["real_maximum_leave_one_out_fdr"].lt(0.05)
        & comparison["real_n_same_direction"].eq(
            comparison["real_n_leave_one_out"]
        )
    )
    synthetic_stable = (
        comparison["synthetic_meta_fdr"].lt(0.05)
        & comparison["synthetic_maximum_leave_one_out_fdr"].lt(0.05)
        & comparison["synthetic_n_same_direction"].eq(
            comparison["synthetic_n_leave_one_out"]
        )
    )
    summary = {
        "features": int(len(comparison)),
        "accessions": int(real_effects["accession"].nunique()),
        "meta_effect_correlation": correlation,
        "meta_direction_agreement": float(comparison["direction_agrees"].mean()),
        "real_random_effects_fdr_lt_005": int(
            comparison["real_meta_fdr"].lt(0.05).sum()
        ),
        "synthetic_random_effects_fdr_lt_005": int(
            comparison["synthetic_meta_fdr"].lt(0.05).sum()
        ),
        "real_loo_stable_fdr_lt_005": int(real_stable.sum()),
        "synthetic_loo_stable_fdr_lt_005": int(synthetic_stable.sum()),
        "concordant_loo_stable_fdr_lt_005": int(
            (real_stable & synthetic_stable & comparison["direction_agrees"]).sum()
        ),
    }
    return {
        "real_per_accession": real_effects,
        "synthetic_per_accession": synthetic_effects,
        "real_meta": real_meta,
        "synthetic_meta": synthetic_meta,
        "real_leave_one_out": real_loo,
        "synthetic_leave_one_out": synthetic_loo,
        "comparison": comparison,
    }, summary
