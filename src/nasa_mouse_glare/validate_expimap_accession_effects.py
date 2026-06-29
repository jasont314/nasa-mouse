"""Validate expiMap condition shifts with accession-aware meta-analysis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .cluster_enrichment import bh_fdr
from .io import require_import


def score_columns(frame) -> list[str]:
    return [
        column
        for column in frame
        if str(column).startswith(("R-MMU-", "MUSCLE_"))
    ]


def random_effects(effect, variance, scipy_stats):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    valid = np.isfinite(effect) & np.isfinite(variance) & (variance > 0)
    y = np.asarray(effect)[valid]
    v = np.asarray(variance)[valid]
    if len(y) < 2:
        return {"n_accessions": int(len(y)), "meta_effect": float("nan"), "meta_se": float("nan"), "meta_p": float("nan"), "tau2": float("nan"), "i2": float("nan")}
    fixed_weights = 1.0 / v
    fixed_effect = float((fixed_weights * y).sum() / fixed_weights.sum())
    q = float((fixed_weights * (y - fixed_effect) ** 2).sum())
    degrees_freedom = len(y) - 1
    c_value = float(fixed_weights.sum() - (fixed_weights**2).sum() / fixed_weights.sum())
    tau2 = max((q - degrees_freedom) / c_value, 0.0) if c_value > 0 else 0.0
    weights = 1.0 / (v + tau2)
    meta_effect = float((weights * y).sum() / weights.sum())
    meta_se = float((1.0 / weights.sum()) ** 0.5)
    z_score = meta_effect / meta_se if meta_se > 0 else float("nan")
    meta_p = float(2 * scipy_stats.norm.sf(abs(z_score))) if np.isfinite(z_score) else float("nan")
    i2 = max((q - degrees_freedom) / q, 0.0) if q > 0 else 0.0
    return {
        "n_accessions": int(len(y)),
        "meta_effect": meta_effect,
        "meta_se": meta_se,
        "meta_p": meta_p,
        "tau2": float(tau2),
        "i2": float(i2),
    }


def accession_effects(scores, terms):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    rows = []
    for accession, frame in scores.groupby("id.accession", dropna=False):
        flight = frame.loc[frame["condition_inferred"].eq("flight")]
        ground = frame.loc[frame["condition_inferred"].eq("ground_control")]
        if flight.empty or ground.empty:
            continue
        for term in terms:
            x = flight[term].astype(float).to_numpy()
            y = ground[term].astype(float).to_numpy()
            variance = float(np.var(x, ddof=1) / len(x) + np.var(y, ddof=1) / len(y))
            rows.append(
                {
                    "id.accession": accession,
                    "term": term,
                    "n_flight": int(len(x)),
                    "n_ground_control": int(len(y)),
                    "flight_minus_ground": float(np.mean(x) - np.mean(y)),
                    "effect_variance": variance,
                }
            )
    return pd.DataFrame(rows)


def meta_analysis(effects):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")
    rows = []
    for term, frame in effects.groupby("term", sort=False):
        result = random_effects(
            frame["flight_minus_ground"].to_numpy(),
            frame["effect_variance"].to_numpy(),
            scipy_stats,
        )
        rows.append({"term": term, **result})
    result = pd.DataFrame(rows)
    direction_rows = []
    for term, frame in effects.groupby("term", sort=False):
        direction = float(result.loc[result["term"].eq(term), "meta_effect"].iloc[0])
        individual = np.sign(frame["flight_minus_ground"].astype(float).to_numpy())
        direction_rows.append(
            {
                "term": term,
                "n_accession_same_direction": int((individual == np.sign(direction)).sum()),
                "n_accession_opposite_direction": int((individual == -np.sign(direction)).sum()),
            }
        )
    result = result.merge(pd.DataFrame(direction_rows), on="term", how="left")
    result["meta_fdr"] = bh_fdr(result["meta_p"].fillna(1.0).to_numpy())
    return result.sort_values(["meta_fdr", "meta_p"], kind="stable")


def leave_one_accession_out(effects):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")
    rows = []
    accessions = sorted(effects["id.accession"].astype(str).unique())
    for held_out in accessions:
        retained = effects.loc[effects["id.accession"].astype(str).ne(held_out)]
        for term, frame in retained.groupby("term", sort=False):
            result = random_effects(
                frame["flight_minus_ground"].to_numpy(),
                frame["effect_variance"].to_numpy(),
                scipy_stats,
            )
            rows.append({"held_out_accession": held_out, "term": term, **result})
    result = pd.DataFrame(rows)
    result["meta_fdr"] = result.groupby("held_out_accession", sort=False)["meta_p"].transform(
        lambda values: bh_fdr(values.fillna(1.0).to_numpy())
    )
    primary = meta_analysis(effects).set_index("term")["meta_effect"]
    summary_rows = []
    for term, frame in result.groupby("term", sort=False):
        direction = float(np.sign(primary.get(term, 0.0)))
        loo_direction = np.sign(frame["meta_effect"].astype(float).to_numpy())
        summary_rows.append(
            {
                "term": term,
                "n_leave_one_out": int(len(frame)),
                "n_same_direction": int((loo_direction == direction).sum()),
                "minimum_leave_one_out_fdr": float(frame["meta_fdr"].min()),
                "maximum_leave_one_out_fdr": float(frame["meta_fdr"].max()),
            }
        )
    return result, pd.DataFrame(summary_rows)


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(args.scores, sep="\t")
    required = {"id.accession", "condition_inferred"}
    missing = required.difference(scores.columns)
    if missing:
        raise SystemExit(f"Scores missing required columns: {sorted(missing)}")
    terms = score_columns(scores)
    if not terms:
        raise SystemExit("No Reactome score columns found.")
    effects = accession_effects(scores, terms)
    meta = meta_analysis(effects)
    loo, loo_summary = leave_one_accession_out(effects)
    paths = {
        "per_accession_effects": output_dir / "per_accession_effects.tsv",
        "random_effects_meta_analysis": output_dir / "random_effects_meta_analysis.tsv",
        "leave_one_accession_out": output_dir / "leave_one_accession_out.tsv",
        "leave_one_out_summary": output_dir / "leave_one_out_summary.tsv",
    }
    effects.to_csv(paths["per_accession_effects"], sep="\t", index=False)
    meta.to_csv(paths["random_effects_meta_analysis"], sep="\t", index=False)
    loo.to_csv(paths["leave_one_accession_out"], sep="\t", index=False)
    loo_summary.to_csv(paths["leave_one_out_summary"], sep="\t", index=False)
    merged = meta.merge(loo_summary, on="term", how="left")
    top = merged.head(10)
    lines = [
        "# Accession-Aware expiMap Validation",
        "",
        f"- Scores: `{args.scores}`",
        f"- Reactome terms: {len(terms)}",
        f"- Accessions with both conditions: {effects['id.accession'].nunique()}",
        f"- Random-effects meta-analysis FDR < 0.05: {int((meta['meta_fdr'] < 0.05).sum())}",
        "",
        "| term | meta effect | meta FDR | I2 | same-direction accessions | same-direction LOO | max LOO FDR |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in top.iterrows():
        lines.append(
            f"| {row.term} | {row.meta_effect:.4g} | {row.meta_fdr:.4g} | {row.i2:.3g} | "
            f"{int(row.n_accession_same_direction)}/{int(row.n_accessions)} | "
            f"{int(row.n_same_direction)}/{int(row.n_leave_one_out)} | {row.maximum_leave_one_out_fdr:.4g} |"
        )
    readme_path = output_dir / "README.md"
    readme_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    summary = {
        "scores": str(args.scores),
        "terms": int(len(terms)),
        "accessions": int(effects["id.accession"].nunique()),
        "meta_fdr_lt_005": int((meta["meta_fdr"] < 0.05).sum()),
        "outputs": {key: str(path) for key, path in paths.items()} | {"readme": str(readme_path)},
    }
    summary_path = output_dir / "validation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run accession-aware random-effects and leave-one-out expiMap validation."
    )
    parser.add_argument("--scores", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
