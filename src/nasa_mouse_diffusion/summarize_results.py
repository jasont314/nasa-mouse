"""Write compact summary tables for diffusion production outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nasa_mouse_glare.io import require_import

from .paths import OUTPUT_ROOT, SUMMARY_DIR


TRACKS = ("osdr_only", "archs4_pretrain_osdr_finetune", "archs4_only")
ANALYSES = (
    ("osdr_only", "post_training", "analysis"),
    ("archs4_pretrain_osdr_finetune", "post_finetune", "analysis"),
    ("archs4_pretrain_osdr_finetune", "frozen_reference_projection", "pretrained_query_analysis"),
)


def flatten(prefix: str, values: dict) -> dict:
    row = {}
    for key, value in values.items():
        if isinstance(value, dict):
            row.update(flatten(f"{prefix}{key}_", value))
        else:
            row[f"{prefix}{key}"] = value
    return row


def training_rows(root: Path) -> list[dict]:
    rows = []
    for track in TRACKS:
        path = root / track / "training_summary.json"
        if not path.exists():
            continue
        summary = json.loads(path.read_text())
        row = {
            "track": track,
            "mode": summary.get("mode", ""),
            "output_dir": summary.get("output_dir", ""),
            "landmark_source": summary.get("landmarks", {}).get("source", ""),
            "device": summary.get("torch", {}).get("model_device", ""),
            "cuda_device": summary.get("torch", {}).get("cuda_device_name", ""),
        }
        row.update(flatten("count_", summary.get("counts", {})))
        row.update(flatten("train_", summary.get("training_design", {})))
        row.update(flatten("quality_", summary.get("quality", {})))
        row.update(flatten("pretrained_query_quality_", summary.get("pretrained_query_quality") or {}))
        row.update(flatten("reconstruction_", summary.get("reconstruction", {})))
        row.update(flatten("reverse_validation_", summary.get("reverse_validation", {})))
        rows.append(row)
    return rows


def analysis_rows(root: Path) -> list[dict]:
    rows = []
    for track, score_set, analysis_dir in ANALYSES:
        path = root / track / analysis_dir / "analysis_summary.json"
        if not path.exists():
            continue
        summary = json.loads(path.read_text())
        row = {
            "track": track,
            "score_set": score_set,
            "analysis_dir": str(root / track / analysis_dir),
            "scores": summary.get("scores", ""),
            "n_samples": summary.get("n_samples", 0),
            "n_features": summary.get("n_features", 0),
            "ordinary_welch_fdr_lt_005": summary.get("ordinary_welch_fdr_lt_005", 0),
            "ordinary_mannwhitney_fdr_lt_005": summary.get("ordinary_mannwhitney_fdr_lt_005", 0),
            "random_effects_fdr_lt_005": summary.get("random_effects_fdr_lt_005", 0),
            "loo_stable_fdr_lt_005": summary.get("loo_stable_fdr_lt_005", 0),
        }
        rows.append(row)
    return rows


def synthetic_rows(root: Path) -> list[dict]:
    manifest = root / "synthetic_examples" / "synthetic_examples_manifest.tsv"
    if not manifest.exists():
        return []
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    rows = []
    table = pd.read_csv(manifest, sep="\t")
    for _, item in table.iterrows():
        row = item.to_dict()
        output_dir = Path(str(row.get("output_dir", "")))
        reports = list(output_dir.glob("*_clip_report.json"))
        lows = []
        highs = []
        values = []
        for report in reports:
            payload = json.loads(report.read_text())
            lows.append(payload.get("n_clipped_low", 0))
            highs.append(payload.get("n_clipped_high", 0))
            values.append(payload.get("n_values", 0))
        total = int(np.sum(values)) if values else 0
        row["clip_report_count"] = len(reports)
        row["n_clipped_low"] = int(np.sum(lows)) if lows else 0
        row["n_clipped_high"] = int(np.sum(highs)) if highs else 0
        row["clip_fraction"] = float((np.sum(lows) + np.sum(highs)) / total) if total else 0.0
        rows.append(row)
    return rows


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    root = Path(args.root)
    summary_dir = Path(args.summary_dir)
    summary_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "diffusion_training_summary.tsv": training_rows(root),
        "diffusion_analysis_summary.tsv": analysis_rows(root),
        "diffusion_synthetic_examples_summary.tsv": synthetic_rows(root),
    }
    manifest_rows = []
    for name, rows in outputs.items():
        path = summary_dir / name
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        manifest_rows.append({"file": str(path), "rows": len(rows)})
    subgroup_path = summary_dir / "diffusion_subgroup_analysis_summary.tsv"
    if subgroup_path.exists():
        subgroup_rows = len(pd.read_csv(subgroup_path, sep="\t"))
        manifest_rows.append({"file": str(subgroup_path), "rows": subgroup_rows})
    reverse_path = summary_dir / "diffusion_reverse_validation_refresh.tsv"
    if reverse_path.exists():
        reverse_rows = len(pd.read_csv(reverse_path, sep="\t"))
        manifest_rows.append({"file": str(reverse_path), "rows": reverse_rows})
    manifest_path = summary_dir / "diffusion_summary_manifest.tsv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)
    print(json.dumps({"manifest": str(manifest_path), "outputs": manifest_rows}, indent=2))
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(OUTPUT_ROOT))
    parser.add_argument("--summary-dir", default=str(SUMMARY_DIR))
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
