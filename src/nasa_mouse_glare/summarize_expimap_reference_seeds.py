"""Summarize seed stability for ARCHS4-reference expiMap query results."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

from .io import require_import


def seed_label(run_dir: Path) -> str:
    match = re.search(r"seed(\d+)", run_dir.name)
    return match.group(1) if match else run_dir.name


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    merged = None
    run_rows = []
    seeds = []
    for raw_run_dir in args.run_dir:
        run_dir = Path(raw_run_dir)
        seed = seed_label(run_dir)
        validation_dir = run_dir / "posterior_mean_accession_validation"
        meta = pd.read_csv(validation_dir / "random_effects_meta_analysis.tsv", sep="\t")
        loo = pd.read_csv(validation_dir / "leave_one_out_summary.tsv", sep="\t")
        table = meta.merge(loo, on="term", how="left")
        columns = {column: f"{column}_seed{seed}" for column in table if column != "term"}
        table = table.rename(columns=columns)
        merged = table if merged is None else merged.merge(table, on="term", how="inner")
        training_summary = run_dir.parent / f"reference_nb_5000_stratified_200epoch_seed{seed}" / "training_summary.json"
        run_rows.append(
            {
                "seed": seed,
                "run_dir": str(run_dir),
                "meta_fdr_lt_threshold": int((meta["meta_fdr"] < args.fdr_threshold).sum()),
                "training_summary": str(training_summary),
            }
        )
        seeds.append(seed)

    fdr_columns = [f"meta_fdr_seed{seed}" for seed in seeds]
    effect_columns = [f"meta_effect_seed{seed}" for seed in seeds]
    loo_columns = [f"maximum_leave_one_out_fdr_seed{seed}" for seed in seeds]
    merged["same_effect_direction_all_seeds"] = (
        merged[effect_columns].gt(0).nunique(axis=1).eq(1)
    )
    merged["fdr_significant_all_seeds"] = merged[fdr_columns].lt(args.fdr_threshold).all(axis=1)
    merged["loo_significant_all_seeds"] = merged[loo_columns].lt(args.fdr_threshold).all(axis=1)
    merged["stable_all_seeds"] = (
        merged["same_effect_direction_all_seeds"]
        & merged["fdr_significant_all_seeds"]
        & merged["loo_significant_all_seeds"]
    )
    merged["maximum_seed_fdr"] = merged[fdr_columns].max(axis=1)
    merged["minimum_seed_abs_effect"] = merged[effect_columns].abs().min(axis=1)
    merged = merged.sort_values(
        ["stable_all_seeds", "maximum_seed_fdr", "minimum_seed_abs_effect"],
        ascending=[False, True, False],
        kind="stable",
    )

    correlations = merged[effect_columns].corr(method="spearman")
    correlations.index = [index.replace("meta_effect_seed", "seed") for index in correlations.index]
    correlations.columns = [column.replace("meta_effect_seed", "seed") for column in correlations.columns]
    stability_path = output_dir / "reference_seed_pathway_stability.tsv"
    correlations_path = output_dir / "reference_seed_effect_spearman.tsv"
    runs_path = output_dir / "reference_seed_runs.tsv"
    merged.to_csv(stability_path, sep="\t", index=False)
    correlations.to_csv(correlations_path, sep="\t")
    pd.DataFrame(run_rows).to_csv(runs_path, sep="\t", index=False)

    stable = merged.loc[merged["stable_all_seeds"]]
    lines = [
        "# ARCHS4 Reference Seed Stability",
        "",
        f"- Query runs compared: {len(seeds)}",
        f"- Meta-analysis FDR threshold: {args.fdr_threshold}",
        f"- Terms FDR-significant in every seed with one shared direction: {int((merged['fdr_significant_all_seeds'] & merged['same_effect_direction_all_seeds']).sum())}",
        f"- Terms additionally significant in every leave-one-accession-out run: {len(stable)}",
        "",
        "## Strictly Stable Terms",
        "",
        "| term | worst seed FDR | smallest absolute effect |",
        "| --- | --- | --- |",
    ]
    if stable.empty:
        lines.append("| None | NA | NA |")
    else:
        for _, row in stable.iterrows():
            lines.append(
                f"| {row.term} | {row.maximum_seed_fdr:.4g} | {row.minimum_seed_abs_effect:.4g} |"
            )
    lines.extend(["", "## Effect-Rank Spearman Correlations", ""])
    headers = ["seed", *correlations.columns.tolist()]
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for seed, row in correlations.iterrows():
        lines.append(
            "| " + " | ".join([str(seed), *[f"{float(value):.4g}" for value in row]]) + " |"
        )
    readme_path = output_dir / "README.md"
    readme_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(readme_path)
    return readme_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize posterior-mean ARCHS4 reference query seed stability."
    )
    parser.add_argument("--run-dir", action="append", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
