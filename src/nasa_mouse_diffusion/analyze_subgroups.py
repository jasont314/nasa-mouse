"""Run diffusion FLT/GC feature analysis within tissues and muscle splits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from nasa_mouse_glare.io import require_import

from . import analyze_features
from .paths import DEFAULT_TISSUES, MUSCLE_GROUPS, OUTPUT_ROOT, SUMMARY_DIR


SCORE_SETS = (
    ("osdr_only", "post_training", "diffusion_feature_scores.tsv"),
    ("archs4_pretrain_osdr_finetune", "post_finetune", "diffusion_feature_scores.tsv"),
    ("archs4_pretrain_osdr_finetune", "frozen_reference_projection", "pretrained_query_diffusion_feature_scores.tsv"),
)


def subgroup_frames(scores, tissues: tuple[str, ...], muscle_groups: tuple[str, ...]):
    for tissue in tissues:
        frame = scores.loc[scores["wgan_tissue"].astype(str).eq(tissue)].copy()
        yield tissue, frame
    muscle = scores.loc[scores["wgan_tissue"].astype(str).eq("skeletal_muscle")]
    for group in muscle_groups:
        frame = muscle.loc[muscle["wgan_muscle_group"].astype(str).eq(group)].copy()
        yield f"skeletal_muscle_{group}", frame


def has_two_conditions(frame) -> bool:
    counts = frame["condition_inferred"].astype(str).value_counts()
    return int(counts.get("flight", 0)) >= 2 and int(counts.get("ground_control", 0)) >= 2


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_root = Path(args.output_root)
    rows = []
    for track, score_set, filename in SCORE_SETS:
        scores_path = output_root / track / filename
        if not scores_path.exists():
            rows.append({"track": track, "score_set": score_set, "status": "missing_scores", "scores": str(scores_path)})
            continue
        scores = pd.read_csv(scores_path, sep="\t")
        for label, frame in subgroup_frames(scores, tuple(args.tissues), tuple(args.muscle_groups)):
            output_dir = output_root / track / "analysis_by_subgroup" / score_set / label
            row = {
                "track": track,
                "score_set": score_set,
                "subgroup": label,
                "scores": str(scores_path),
                "output_dir": str(output_dir),
                "n_samples": int(len(frame)),
                "n_flight": int(frame["condition_inferred"].astype(str).eq("flight").sum()) if not frame.empty else 0,
                "n_ground_control": int(frame["condition_inferred"].astype(str).eq("ground_control").sum()) if not frame.empty else 0,
            }
            if len(frame) < args.min_samples or not has_two_conditions(frame):
                row["status"] = "skipped_small_or_one_condition"
                rows.append(row)
                continue
            output_dir.mkdir(parents=True, exist_ok=True)
            filtered_scores = output_dir / "subgroup_diffusion_feature_scores.tsv"
            frame.to_csv(filtered_scores, sep="\t", index=False)
            analyze_features.run(
                SimpleNamespace(
                    scores=str(filtered_scores),
                    output_dir=str(output_dir),
                    top_features=args.top_features,
                )
            )
            summary_path = output_dir / "analysis_summary.json"
            summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
            row.update(
                {
                    "status": "completed",
                    "ordinary_welch_fdr_lt_005": summary.get("ordinary_welch_fdr_lt_005", 0),
                    "ordinary_mannwhitney_fdr_lt_005": summary.get("ordinary_mannwhitney_fdr_lt_005", 0),
                    "random_effects_fdr_lt_005": summary.get("random_effects_fdr_lt_005", 0),
                    "loo_stable_fdr_lt_005": summary.get("loo_stable_fdr_lt_005", 0),
                }
            )
            rows.append(row)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = SUMMARY_DIR / "diffusion_subgroup_analysis_summary.tsv"
    pd.DataFrame(rows).to_csv(manifest_path, sep="\t", index=False)
    (SUMMARY_DIR / "diffusion_subgroup_analysis_summary.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": str(manifest_path), "rows": len(rows)}, indent=2))
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default=str(OUTPUT_ROOT))
    parser.add_argument("--tissues", nargs="+", default=list(DEFAULT_TISSUES))
    parser.add_argument("--muscle-groups", nargs="+", default=list(MUSCLE_GROUPS))
    parser.add_argument("--min-samples", type=int, default=10)
    parser.add_argument("--top-features", type=int, default=30)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
