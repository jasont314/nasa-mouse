"""Estimate the cross-study FLT/GC signal ceiling in real OSDR profiles."""

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
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


CONDITIONS = ("ground_control", "flight")


def _correlation(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    finite = np.isfinite(left) & np.isfinite(right)
    if finite.sum() < 3 or np.std(left[finite]) == 0 or np.std(right[finite]) == 0:
        return float("nan")
    return float(spearmanr(left[finite], right[finite]).statistic)


def _direction_agreement(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    finite = np.isfinite(left) & np.isfinite(right)
    if not finite.any():
        return float("nan")
    return float(np.mean(np.sign(left[finite]) == np.sign(right[finite])))


def _condition_effect(expression: np.ndarray, conditions: np.ndarray) -> np.ndarray:
    conditions = np.asarray(conditions, dtype=str)
    flight = conditions == "flight"
    ground = conditions == "ground_control"
    if not flight.any() or not ground.any():
        raise ValueError("Both flight and ground_control profiles are required")
    return expression[flight].mean(axis=0) - expression[ground].mean(axis=0)


def _eligible_accessions(
    samples: pd.DataFrame, *, minimum_per_condition: int
) -> list[str]:
    counts = (
        samples.groupby(["accession", "condition"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    for condition in CONDITIONS:
        if condition not in counts:
            counts[condition] = 0
    return sorted(
        counts.index[
            counts.loc[:, list(CONDITIONS)].min(axis=1) >= minimum_per_condition
        ].astype(str)
    )


def _accession_effects(
    expression: np.ndarray,
    samples: pd.DataFrame,
    *,
    minimum_per_condition: int,
) -> dict[str, np.ndarray]:
    effects: dict[str, np.ndarray] = {}
    for accession in _eligible_accessions(
        samples, minimum_per_condition=minimum_per_condition
    ):
        mask = samples["accession"].astype(str).eq(accession).to_numpy()
        effects[accession] = _condition_effect(
            expression[mask], samples.loc[mask, "condition"].astype(str).to_numpy()
        )
    return effects


def _loo_effect_rows(effects: dict[str, np.ndarray]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for accession, heldout in effects.items():
        others = [effect for key, effect in effects.items() if key != accession]
        if not others:
            continue
        reference = np.mean(np.stack(others), axis=0)
        rows.append(
            {
                "accession": accession,
                "effect_correlation": _correlation(reference, heldout),
                "effect_direction_agreement": _direction_agreement(
                    reference, heldout
                ),
            }
        )
    return rows


def _loo_classifier_rows(
    expression: np.ndarray,
    samples: pd.DataFrame,
    accessions: Iterable[str],
    *,
    seed: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    labels = samples["condition"].astype(str).eq("flight").to_numpy(dtype=int)
    accession_values = samples["accession"].astype(str).to_numpy()
    for offset, accession in enumerate(accessions):
        test = accession_values == accession
        train = np.isin(accession_values, list(accessions)) & ~test
        if min(np.unique(labels[train]).size, np.unique(labels[test]).size) < 2:
            continue
        classifier = make_pipeline(
            StandardScaler(),
            LogisticRegression(
                C=0.1,
                class_weight="balanced",
                max_iter=4000,
                random_state=seed + offset,
                solver="liblinear",
            ),
        )
        classifier.fit(expression[train], labels[train])
        probability = classifier.predict_proba(expression[test])[:, 1]
        prediction = probability >= 0.5
        rows.append(
            {
                "accession": accession,
                "balanced_accuracy": float(
                    balanced_accuracy_score(labels[test], prediction)
                ),
                "roc_auc": float(roc_auc_score(labels[test], probability)),
                "test_profiles": int(test.sum()),
            }
        )
    return rows


def _shuffle_within_accession(
    samples: pd.DataFrame, rng: np.random.Generator
) -> pd.DataFrame:
    shuffled = samples.copy()
    for _, indices in shuffled.groupby("accession", observed=True).groups.items():
        positions = np.asarray(list(indices), dtype=int)
        shuffled.loc[positions, "condition"] = rng.permutation(
            shuffled.loc[positions, "condition"].to_numpy()
        )
    return shuffled


def _finite_median(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    return float(np.median(array)) if len(array) else float("nan")


def _empirical_upper_p(observed: float, null: Iterable[float]) -> float:
    values = np.asarray(list(null), dtype=float)
    values = values[np.isfinite(values)]
    if not np.isfinite(observed) or not len(values):
        return float("nan")
    return float((1 + np.sum(values >= observed)) / (len(values) + 1))


def analyze_tissue(
    expression: np.ndarray,
    samples: pd.DataFrame,
    *,
    tissue: str,
    minimum_per_condition: int = 2,
    permutation_repeats: int = 100,
    seed: int = 1234,
) -> tuple[dict[str, object], pd.DataFrame]:
    """Evaluate cross-accession condition reproducibility for one tissue."""

    mask = samples["tissue"].astype(str).eq(tissue).to_numpy()
    tissue_expression = np.asarray(expression[mask], dtype=np.float64)
    tissue_samples = samples.loc[mask].reset_index(drop=True)
    effects = _accession_effects(
        tissue_expression,
        tissue_samples,
        minimum_per_condition=minimum_per_condition,
    )
    effect_rows = _loo_effect_rows(effects)
    classifier_rows = _loo_classifier_rows(
        tissue_expression,
        tissue_samples,
        effects,
        seed=seed,
    )
    effect_table = pd.DataFrame(
        effect_rows,
        columns=(
            "accession",
            "effect_correlation",
            "effect_direction_agreement",
        ),
    )
    classifier_table = pd.DataFrame(
        classifier_rows,
        columns=("accession", "balanced_accuracy", "roc_auc", "test_profiles"),
    )
    detail = effect_table.merge(classifier_table, on="accession", how="outer")
    detail.insert(0, "tissue", tissue)

    observed_correlation = _finite_median(
        row["effect_correlation"] for row in effect_rows
    )
    observed_direction = _finite_median(
        row["effect_direction_agreement"] for row in effect_rows
    )
    null_correlation: list[float] = []
    null_direction: list[float] = []
    rng = np.random.default_rng(seed)
    for _ in range(int(permutation_repeats)):
        shuffled = _shuffle_within_accession(tissue_samples, rng)
        shuffled_effects = _accession_effects(
            tissue_expression,
            shuffled,
            minimum_per_condition=minimum_per_condition,
        )
        shuffled_rows = _loo_effect_rows(shuffled_effects)
        null_correlation.append(
            _finite_median(row["effect_correlation"] for row in shuffled_rows)
        )
        null_direction.append(
            _finite_median(
                row["effect_direction_agreement"] for row in shuffled_rows
            )
        )

    summary = {
        "tissue": tissue,
        "profiles": int(len(tissue_samples)),
        "accessions": int(tissue_samples["accession"].nunique()),
        "eligible_accessions": int(len(effects)),
        "loo_effect_correlation_median": observed_correlation,
        "loo_effect_correlation_minimum": (
            float(detail["effect_correlation"].min()) if len(detail) else float("nan")
        ),
        "loo_direction_agreement_median": observed_direction,
        "loo_balanced_accuracy_median": _finite_median(
            row["balanced_accuracy"] for row in classifier_rows
        ),
        "loo_balanced_accuracy_minimum": (
            float(detail["balanced_accuracy"].min()) if len(detail) else float("nan")
        ),
        "loo_roc_auc_median": _finite_median(
            row["roc_auc"] for row in classifier_rows
        ),
        "permutation_repeats": int(permutation_repeats),
        "effect_correlation_null_median": _finite_median(null_correlation),
        "effect_correlation_empirical_p": _empirical_upper_p(
            observed_correlation, null_correlation
        ),
        "direction_agreement_null_median": _finite_median(null_direction),
        "direction_agreement_empirical_p": _empirical_upper_p(
            observed_direction, null_direction
        ),
    }
    return summary, detail


def load_development_expression(
    prepared_h5: str | Path,
    samples_tsv: str | Path,
    *,
    roles: Iterable[str] = ("train",),
    transform: str = "log1p",
) -> tuple[np.ndarray, pd.DataFrame]:
    """Load declared non-test roles and align them by immutable source row."""

    requested = tuple(map(str, roles))
    if "test" in requested:
        raise ValueError(
            "The real-data ceiling command cannot open the locked test role"
        )
    samples = pd.read_csv(samples_tsv, sep="\t")
    expression_blocks: list[np.ndarray] = []
    metadata_blocks: list[pd.DataFrame] = []
    with h5py.File(prepared_h5, "r") as handle:
        for role in requested:
            if role not in handle:
                raise ValueError(f"Prepared data has no role named {role!r}")
            group = handle[role]
            key = "analysis_expression" if "analysis_expression" in group else "tpm"
            values = np.asarray(group[key][:], dtype=np.float64)
            source_rows = np.asarray(group["source_row"][:], dtype=np.int64)
            role_samples = samples.loc[samples["role"].astype(str).eq(role)].copy()
            lookup = role_samples.set_index("_row_index", drop=False)
            missing = sorted(set(source_rows) - set(lookup.index.astype(int)))
            if missing:
                raise ValueError(f"Metadata is missing source rows: {missing[:5]}")
            aligned = lookup.loc[source_rows].reset_index(drop=True)
            if len(aligned) != len(values):
                raise ValueError("Expression and metadata lengths differ")
            expression_blocks.append(values)
            metadata_blocks.append(aligned)
    expression = np.concatenate(expression_blocks)
    metadata = pd.concat(metadata_blocks, ignore_index=True)
    if transform == "log1p":
        if (expression < 0).any():
            raise ValueError("log1p ceiling analysis requires nonnegative expression")
        expression = np.log1p(expression)
    elif transform != "none":
        raise ValueError(f"Unsupported ceiling transform: {transform}")
    return expression.astype(np.float32), metadata


def _plot_summary(summary: pd.DataFrame, output: Path) -> None:
    if summary.empty:
        return
    ordered = summary.sort_values("loo_effect_correlation_median")
    positions = np.arange(len(ordered))
    figure, axes = plt.subplots(1, 2, figsize=(13.0, max(5.0, len(ordered) * 0.34)))
    axes[0].barh(positions, ordered["loo_effect_correlation_median"], color="#2878B5")
    axes[0].axvline(0.0, color="#333333", linewidth=1)
    axes[0].set_yticks(positions, ordered["tissue"])
    axes[0].set_xlabel("Median LOO effect Spearman correlation")
    axes[1].barh(positions, ordered["loo_balanced_accuracy_median"], color="#C14924")
    axes[1].axvline(0.5, color="#333333", linewidth=1, linestyle="--")
    axes[1].set_yticks(positions, [""] * len(positions))
    axes[1].set_xlim(0.0, 1.0)
    axes[1].set_xlabel("Median LOO balanced accuracy")
    figure.suptitle(
        "Real OSDR cross-accession FLT/GC signal ceiling", fontweight="bold"
    )
    figure.tight_layout()
    figure.savefig(output / "real_flt_gc_ceiling.png", dpi=220, bbox_inches="tight")
    figure.savefig(output / "real_flt_gc_ceiling.pdf", bbox_inches="tight")
    plt.close(figure)


def run_ceiling(
    prepared_h5: str | Path,
    samples_tsv: str | Path,
    output_dir: str | Path,
    *,
    roles: Iterable[str] = ("train",),
    transform: str = "log1p",
    minimum_per_condition: int = 2,
    minimum_accessions: int = 2,
    permutation_repeats: int = 100,
    seed: int = 1234,
) -> Path:
    expression, samples = load_development_expression(
        prepared_h5, samples_tsv, roles=roles, transform=transform
    )
    summaries: list[dict[str, object]] = []
    details: list[pd.DataFrame] = []
    for tissue in sorted(samples["tissue"].dropna().astype(str).unique()):
        summary, detail = analyze_tissue(
            expression,
            samples,
            tissue=tissue,
            minimum_per_condition=minimum_per_condition,
            permutation_repeats=permutation_repeats,
            seed=seed,
        )
        if int(summary["eligible_accessions"]) >= int(minimum_accessions):
            summaries.append(summary)
            details.append(detail)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    summary_table = pd.DataFrame(summaries).sort_values(
        ["loo_effect_correlation_median", "tissue"], ascending=[False, True]
    )
    detail_table = pd.concat(details, ignore_index=True) if details else pd.DataFrame()
    summary_table.to_csv(output / "tissue_summary.tsv", sep="\t", index=False)
    detail_table.to_csv(output / "accession_loo.tsv", sep="\t", index=False)
    _plot_summary(summary_table, output)
    audit = {
        "prepared_h5": str(Path(prepared_h5).resolve()),
        "samples_tsv": str(Path(samples_tsv).resolve()),
        "roles": list(map(str, roles)),
        "locked_test_opened": False,
        "transform": transform,
        "minimum_profiles_per_condition_per_accession": minimum_per_condition,
        "minimum_accessions": minimum_accessions,
        "permutation_repeats": permutation_repeats,
        "seed": seed,
        "profiles": int(len(samples)),
        "tissues_reported": int(len(summary_table)),
    }
    (output / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-h5", required=True)
    parser.add_argument("--samples-tsv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--roles", nargs="+", default=["train"])
    parser.add_argument("--transform", choices=("log1p", "none"), default="log1p")
    parser.add_argument("--minimum-per-condition", type=int, default=2)
    parser.add_argument("--minimum-accessions", type=int, default=2)
    parser.add_argument("--permutation-repeats", type=int, default=100)
    parser.add_argument("--seed", type=int, default=1234)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_ceiling(
        args.prepared_h5,
        args.samples_tsv,
        args.output_dir,
        roles=args.roles,
        transform=args.transform,
        minimum_per_condition=args.minimum_per_condition,
        minimum_accessions=args.minimum_accessions,
        permutation_repeats=args.permutation_repeats,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
