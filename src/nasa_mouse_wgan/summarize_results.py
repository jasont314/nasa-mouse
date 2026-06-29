"""Summarize completed WGAN-GP tissue outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nasa_mouse_glare.io import require_import


DEFAULT_SEARCH_ROOT = Path("outputs")
DEFAULT_OUTPUT_DIR = Path("outputs/wgan_pipeline/summary")
EXCLUDED_OUTPUT_ROOTS = {"wgan_pipeline"}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_output_dir(path: Path) -> dict:
    parts = path.parts
    tissue = ""
    group = ""
    if len(parts) >= 2 and parts[0] == "outputs":
        root = parts[1]
        if root == "wgan_skeletal_muscle_splits" and len(parts) >= 3:
            tissue = "skeletal_muscle"
            group = parts[2]
        elif root.startswith("wgan_"):
            tissue = root.removeprefix("wgan_")
    return {"tissue": tissue, "group": group}


def discover_run_dirs(root: Path) -> list[Path]:
    run_dirs = []
    for path in root.glob("wgan_*/**/training_summary.json"):
        relative = path.relative_to(root)
        output_root = relative.parts[0] if relative.parts else ""
        if output_root in EXCLUDED_OUTPUT_ROOTS or output_root.startswith("wgan_smoke"):
            continue
        run_dirs.append(path.parent)
    return sorted(run_dirs)


def stable_count(meta, loo):
    if meta.empty or loo.empty:
        return 0
    merged = meta.merge(loo, on="feature", how="inner")
    return int(
        (
            (merged["meta_fdr"] < 0.05)
            & (merged["maximum_leave_one_out_fdr"] < 0.05)
            & (merged["n_same_direction"] == merged["n_leave_one_out"])
        ).sum()
    )


def summarize_score_set(run_dir: Path, score_set: str, analysis_dir: Path, training: dict):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    comparison_path = analysis_dir / "flt_vs_gc_wgan_feature_comparison.tsv"
    meta_path = analysis_dir / "random_effects_meta_analysis.tsv"
    loo_path = analysis_dir / "leave_one_out_summary.tsv"
    comparison = pd.read_csv(comparison_path, sep="\t") if comparison_path.exists() else pd.DataFrame()
    meta = pd.read_csv(meta_path, sep="\t") if meta_path.exists() else pd.DataFrame()
    loo = pd.read_csv(loo_path, sep="\t") if loo_path.exists() else pd.DataFrame()

    run_info = parse_output_dir(run_dir)
    counts = training.get("counts", {})
    torch = training.get("torch", {})
    row = {
        **run_info,
        "mode": training.get("mode", run_dir.name),
        "score_set": score_set,
        "output_dir": str(run_dir),
        "query_samples": counts.get("query_samples", ""),
        "reference_samples": counts.get("reference_samples", ""),
        "genes": counts.get("genes", ""),
        "critic_features": counts.get("critic_features", ""),
        "model_device": torch.get("model_device", ""),
        "cuda_device_name": torch.get("cuda_device_name", ""),
        "ordinary_welch_fdr_lt_005": 0,
        "random_effects_fdr_lt_005": 0,
        "loo_stable_fdr_lt_005": 0,
        "min_welch_fdr": "",
        "top_welch_feature": "",
        "top_welch_effect": "",
        "min_meta_fdr": "",
        "top_meta_feature": "",
        "top_meta_effect": "",
        "top_meta_max_loo_fdr": "",
    }
    top_rows = []
    if not comparison.empty:
        comparison = comparison.sort_values(["welch_fdr", "welch_p"], kind="stable")
        row["ordinary_welch_fdr_lt_005"] = int((comparison["welch_fdr"] < 0.05).sum())
        top = comparison.iloc[0]
        row["min_welch_fdr"] = float(top["welch_fdr"])
        row["top_welch_feature"] = str(top["feature"])
        row["top_welch_effect"] = float(top["flight_minus_ground"])
        for _, item in comparison.head(10).iterrows():
            top_rows.append(
                {
                    **run_info,
                    "mode": row["mode"],
                    "score_set": score_set,
                    "rank_source": "ordinary_welch",
                    "feature": item["feature"],
                    "effect": item["flight_minus_ground"],
                    "fdr": item["welch_fdr"],
                    "output_dir": str(run_dir),
                }
            )
    if not meta.empty:
        meta = meta.sort_values(["meta_fdr", "meta_p"], kind="stable")
        row["random_effects_fdr_lt_005"] = int((meta["meta_fdr"] < 0.05).sum())
        row["loo_stable_fdr_lt_005"] = stable_count(meta, loo)
        top = meta.iloc[0]
        row["min_meta_fdr"] = float(top["meta_fdr"])
        row["top_meta_feature"] = str(top["feature"])
        row["top_meta_effect"] = float(top["meta_effect"])
        if not loo.empty:
            match = loo.loc[loo["feature"].eq(top["feature"])]
            if not match.empty:
                row["top_meta_max_loo_fdr"] = float(match.iloc[0]["maximum_leave_one_out_fdr"])
        for _, item in meta.head(10).iterrows():
            top_rows.append(
                {
                    **run_info,
                    "mode": row["mode"],
                    "score_set": score_set,
                    "rank_source": "random_effects",
                    "feature": item["feature"],
                    "effect": item["meta_effect"],
                    "fdr": item["meta_fdr"],
                    "output_dir": str(run_dir),
                }
            )
    return row, top_rows


def summarize_run(run_dir: Path):
    training = load_json(run_dir / "training_summary.json")
    rows = []
    top_rows = []
    row, top = summarize_score_set(
        run_dir, "finetuned_or_direct", run_dir / "analysis", training
    )
    rows.append(row)
    top_rows.extend(top)
    if (run_dir / "pretrained_query_critic_feature_scores.tsv").exists():
        row, top = summarize_score_set(
            run_dir,
            "pre_finetune_projection",
            run_dir / "pretrained_query_analysis",
            training,
        )
        rows.append(row)
        top_rows.extend(top)
    return rows, top_rows


def write_readme(output_dir: Path, summary, top_terms) -> None:
    lines = [
        "# WGAN-GP Summary",
        "",
        f"- Score sets summarized: {len(summary)}",
        f"- Runs with random-effects FDR < 0.05: {int((summary['random_effects_fdr_lt_005'] > 0).sum()) if not summary.empty else 0}",
        f"- Runs with strict LOO-stable FDR < 0.05: {int((summary['loo_stable_fdr_lt_005'] > 0).sum()) if not summary.empty else 0}",
        "",
        "Primary files:",
        "",
        "- `wgan_run_summary.tsv`",
        "- `wgan_top_features.tsv`",
    ]
    if not summary.empty:
        lines.extend(
            [
                "",
                "| tissue | group | mode | score set | random FDR hits | LOO-stable hits | top meta feature |",
                "| --- | --- | --- | --- | ---: | ---: | --- |",
            ]
        )
        for _, row in summary.sort_values(
            ["loo_stable_fdr_lt_005", "random_effects_fdr_lt_005"],
            ascending=False,
        ).head(15).iterrows():
            lines.append(
                f"| {row['tissue']} | {row['group']} | {row['mode']} | {row['score_set']} | "
                f"{row['random_effects_fdr_lt_005']} | {row['loo_stable_fdr_lt_005']} | "
                f"{row['top_meta_feature']} |"
            )
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    top_rows = []
    for run_dir in discover_run_dirs(Path(args.search_root)):
        run_rows, run_top = summarize_run(run_dir)
        rows.extend(run_rows)
        top_rows.extend(run_top)
    summary = pd.DataFrame(rows)
    top = pd.DataFrame(top_rows)
    summary_path = output_dir / "wgan_run_summary.tsv"
    top_path = output_dir / "wgan_top_features.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    top.to_csv(top_path, sep="\t", index=False)
    write_readme(output_dir, summary, top)
    print(
        json.dumps(
            {
                "score_sets": int(len(summary)),
                "top_rows": int(len(top)),
                "summary": str(summary_path),
                "top_features": str(top_path),
            },
            indent=2,
        )
    )
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search-root", type=Path, default=DEFAULT_SEARCH_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
