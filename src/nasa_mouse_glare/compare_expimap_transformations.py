"""Compare expiMap preprocessing/transformation sensitivity runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import require_import


DEFAULT_RUNS = {
    "raw_counts": "raw_counts_nb_50epoch",
    "cpm": "cpm_mse_50epoch",
    "log1p_cpm": "log1p_cpm_mse_50epoch",
}


def load_run(tissue_dir: Path, label: str, run_name: str):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    run_dir = tissue_dir / run_name
    training_path = run_dir / "training_summary.json"
    analysis_path = run_dir / "analysis" / "analysis_summary.json"
    comparison_path = run_dir / "analysis" / "flt_vs_gc_pathway_comparison.tsv"
    if not training_path.exists():
        raise FileNotFoundError(training_path)
    if not comparison_path.exists():
        raise FileNotFoundError(comparison_path)

    training = json.loads(training_path.read_text(encoding="utf-8"))
    analysis = (
        json.loads(analysis_path.read_text(encoding="utf-8"))
        if analysis_path.exists()
        else {}
    )
    comparison = pd.read_csv(comparison_path, sep="\t")
    comparison["transformation"] = label
    return {
        "label": label,
        "run_name": run_name,
        "run_dir": str(run_dir),
        "training": training,
        "analysis": analysis,
        "comparison": comparison,
    }


def fdr_count(frame, column: str = "welch_fdr", threshold: float = 0.1) -> int:
    if column not in frame:
        return 0
    return int(frame[column].fillna(1.0).le(threshold).sum())


def summarize_runs(runs):
    rows = []
    for run in runs:
        training = run["training"]
        input_path = Path(training["input"])
        preprocessing = {}
        try:
            ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
            adata = ad.read_h5ad(input_path, backed="r")
            preprocessing = dict(adata.uns.get("expimap_preprocessing", {}))
            adata.file.close()
        except Exception:
            preprocessing = {}
        comparison = run["comparison"]
        top = comparison.sort_values("welch_fdr").head(1)
        rows.append(
            {
                "transformation": run["label"],
                "run_dir": run["run_dir"],
                "recon_loss": training.get("recon_loss", ""),
                "epochs": training.get("epochs", ""),
                "validity": preprocessing.get("validity", ""),
                "normalization_method": preprocessing.get("normalization_method", ""),
                "library_size_handling": preprocessing.get("library_size_handling", ""),
                "n_samples": training.get("n_samples", ""),
                "n_genes": training.get("n_genes", ""),
                "n_terms": training.get("n_terms", ""),
                "min_welch_fdr": float(comparison["welch_fdr"].min())
                if len(comparison)
                else None,
                "n_welch_fdr_0_10": fdr_count(comparison, "welch_fdr", 0.10),
                "n_mannwhitney_fdr_0_10": fdr_count(
                    comparison,
                    "mannwhitney_fdr",
                    0.10,
                ),
                "top_welch_term": top["term"].iloc[0] if len(top) else "",
                "top_welch_effect": float(top["flight_minus_ground"].iloc[0])
                if len(top)
                else None,
            }
        )
    return rows


def compare_effects(runs):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    for i, left in enumerate(runs):
        left_effect = left["comparison"].set_index("term")["flight_minus_ground"]
        left_top = set(
            left["comparison"]
            .assign(abs_effect=left["comparison"]["flight_minus_ground"].abs())
            .sort_values(["abs_effect", "welch_fdr"], ascending=[False, True])
            .head(50)["term"]
        )
        for right in runs[i + 1 :]:
            right_effect = right["comparison"].set_index("term")["flight_minus_ground"]
            shared = left_effect.index.intersection(right_effect.index)
            rho = scipy_stats.spearmanr(
                left_effect.loc[shared],
                right_effect.loc[shared],
                nan_policy="omit",
            )
            right_top = set(
                right["comparison"]
                .assign(abs_effect=right["comparison"]["flight_minus_ground"].abs())
                .sort_values(["abs_effect", "welch_fdr"], ascending=[False, True])
                .head(50)["term"]
            )
            rows.append(
                {
                    "left_transformation": left["label"],
                    "right_transformation": right["label"],
                    "shared_terms": int(len(shared)),
                    "spearman_effect_rho": float(rho.statistic),
                    "spearman_effect_p": float(rho.pvalue),
                    "top50_overlap": int(len(left_top & right_top)),
                    "top50_jaccard": (
                        float(len(left_top & right_top) / len(left_top | right_top))
                        if left_top | right_top
                        else 0.0
                    ),
                }
            )
    return pd.DataFrame(rows)


def write_markdown(output_dir: Path, tissue: str, summary, correlations) -> Path:
    lines = [
        f"# expiMap Preprocessing Comparison: {tissue}",
        "",
        "These runs compare direct OSDR expiMap inputs built from the NASA OSDR API.",
        "`raw_counts` with NB likelihood is the primary analysis. `cpm` and",
        "`log1p_cpm` use MSE and are sensitivity analyses only.",
        "",
        "TPM/log1p(TPM) were not run because the selected OSDR API files are",
        "unnormalized count tables and no transcript length or TPM layer was used.",
        "Z-scored inputs were not run because installed expiMap applies log-style",
        "handling in the MSE path and z-scores can be negative; that would violate",
        "the model's expected nonnegative expression input.",
        "",
        "## Run Summary",
        "",
    ]
    summary_cols = [
        "transformation",
        "recon_loss",
        "validity",
        "min_welch_fdr",
        "n_welch_fdr_0_10",
        "top_welch_term",
    ]
    lines.append("| " + " | ".join(summary_cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(summary_cols)) + " |")
    for row in summary:
        values = []
        for column in summary_cols:
            value = row.get(column, "")
            if isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")

    lines += ["", "## Effect Correlations", ""]
    if correlations.empty:
        lines.append("No pairwise correlations available.")
    else:
        corr_cols = [
            "left_transformation",
            "right_transformation",
            "spearman_effect_rho",
            "top50_overlap",
            "top50_jaccard",
        ]
        lines.append("| " + " | ".join(corr_cols) + " |")
        lines.append("| " + " | ".join(["---"] * len(corr_cols)) + " |")
        for _, row in correlations[corr_cols].iterrows():
            values = []
            for column in corr_cols:
                value = row[column]
                if isinstance(value, float):
                    values.append(f"{value:.6g}")
                else:
                    values.append(str(value))
            lines.append("| " + " | ".join(values) + " |")
    path = output_dir / "README.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    tissue_dir = Path(args.tissue_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_map = dict(DEFAULT_RUNS)
    for item in args.run or []:
        label, run_name = item.split("=", 1)
        run_map[label] = run_name

    runs = [load_run(tissue_dir, label, run_name) for label, run_name in run_map.items()]
    summary = summarize_runs(runs)
    summary_frame = pd.DataFrame(summary)
    summary_path = output_dir / "preprocessing_run_summary.tsv"
    summary_frame.to_csv(summary_path, sep="\t", index=False)

    correlations = compare_effects(runs)
    correlation_path = output_dir / "preprocessing_effect_correlations.tsv"
    correlations.to_csv(correlation_path, sep="\t", index=False)

    combined = pd.concat([run["comparison"] for run in runs], ignore_index=True)
    combined_path = output_dir / "preprocessing_all_pathway_comparisons.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)

    readme_path = write_markdown(output_dir, args.tissue, summary, correlations)
    manifest = {
        "tissue": args.tissue,
        "tissue_dir": str(tissue_dir),
        "runs": run_map,
        "outputs": {
            "summary": str(summary_path),
            "effect_correlations": str(correlation_path),
            "all_pathway_comparisons": str(combined_path),
            "readme": str(readme_path),
        },
    }
    manifest_path = output_dir / "preprocessing_comparison_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare direct expiMap preprocessing transformation runs."
    )
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--tissue-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--run",
        action="append",
        help="Optional label=run_dir_name override/addition.",
    )
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
