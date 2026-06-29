"""Generate per-tissue synthetic FLT/GC examples from trained diffusion models."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from types import SimpleNamespace

from nasa_mouse_glare.io import require_import

from . import generate_synthetic
from .paths import DEFAULT_TISSUES, MUSCLE_GROUPS, OUTPUT_ROOT


DEFAULT_MODELS = ("osdr_only", "archs4_pretrain_osdr_finetune")


def slugify(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return text.strip("_").lower() or "unknown"


def profile_sets(model_dir: Path, tissues: tuple[str, ...], muscle_groups: tuple[str, ...]) -> list[dict]:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    profiles_path = model_dir / "observed_conditioning_profiles.tsv"
    if not profiles_path.exists():
        return []
    profiles = pd.read_csv(profiles_path, sep="\t")
    if "training_source" in profiles.columns:
        query = profiles.loc[profiles["training_source"].astype(str).eq("query")].copy()
        if not query.empty:
            profiles = query
    rows = []
    for tissue in tissues:
        tissue_frame = profiles.loc[profiles["wgan_tissue"].astype(str).eq(tissue)]
        if tissue_frame.empty:
            rows.append({"label": tissue, "status": "missing_profile", "reason": "no observed query profile"})
            continue
        rows.append({"label": tissue, "status": "planned", "profile": tissue_frame.iloc[0].to_dict()})
    muscle_frame = profiles.loc[profiles["wgan_tissue"].astype(str).eq("skeletal_muscle")]
    for group in muscle_groups:
        group_frame = muscle_frame.loc[muscle_frame["wgan_muscle_group"].astype(str).eq(group)]
        label = f"skeletal_muscle_{group}"
        if group_frame.empty:
            rows.append({"label": label, "status": "missing_profile", "reason": "no observed query profile"})
            continue
        rows.append({"label": label, "status": "planned", "profile": group_frame.iloc[0].to_dict()})
    return rows


def profile_to_set_args(profile: dict) -> list[str]:
    values = []
    for key, value in profile.items():
        if key == "training_source" or key == "wgan_condition":
            continue
        if key.startswith("wgan_"):
            values.append(f"{key}={value}")
    return values


def write_readme(output_dir: Path, manifest_path: Path, rows: list[dict]) -> None:
    completed = sum(1 for row in rows if row.get("status") == "completed")
    lines = [
        "# Diffusion Synthetic Examples",
        "",
        "Matched counterfactual examples generated from trained conditional diffusion models.",
        "Each profile holds tissue/material/accession/sex/assay/platform/source fixed and flips `ground_control` to `flight`.",
        "",
        f"- Completed profile/model runs: {completed}",
        f"- Manifest: `{manifest_path.name}`",
        "",
        "Each completed directory contains scaled, log1p CPM, CPM, condition profile JSON files, and a mean FLT-GC log1p CPM delta table.",
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for model_name in args.models:
        model_dir = Path(args.models_root) / model_name
        for profile_row in profile_sets(model_dir, tuple(args.tissues), tuple(args.muscle_groups)):
            label = profile_row["label"]
            row = {"model": model_name, "profile_label": label, "model_dir": str(model_dir)}
            if profile_row.get("status") != "planned":
                row.update(profile_row)
                rows.append(row)
                continue
            profile = profile_row["profile"]
            output_dir = output_root / model_name / slugify(label)
            summary_path = output_dir / "synthetic_generation_summary.json"
            if summary_path.exists() and not args.overwrite:
                row.update({"status": "completed", "output_dir": str(output_dir), "summary": str(summary_path)})
                rows.append(row)
                continue
            summary = generate_synthetic.run(
                SimpleNamespace(
                    model_dir=str(model_dir),
                    output_dir=str(output_dir),
                    n=args.n,
                    condition="",
                    counterfactual=("ground_control", "flight"),
                    set=profile_to_set_args(profile),
                    sample_steps=args.sample_steps,
                    eta=args.eta,
                    batch_size=args.batch_size,
                    seed=args.seed,
                    cpu=args.cpu,
                )
            )
            row.update(
                {
                    "status": "completed",
                    "output_dir": str(output_dir),
                    "summary": str(summary),
                    "wgan_tissue": profile.get("wgan_tissue", ""),
                    "wgan_material_type": profile.get("wgan_material_type", ""),
                    "wgan_muscle_group": profile.get("wgan_muscle_group", ""),
                    "wgan_accession": profile.get("wgan_accession", ""),
                    "wgan_sex": profile.get("wgan_sex", ""),
                    "wgan_platform": profile.get("wgan_platform", ""),
                }
            )
            rows.append(row)
    manifest_path = output_root / "synthetic_examples_manifest.tsv"
    pd.DataFrame(rows).to_csv(manifest_path, sep="\t", index=False)
    (output_root / "synthetic_examples_manifest.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    write_readme(output_root, manifest_path, rows)
    print(json.dumps({"manifest": str(manifest_path), "rows": rows}, indent=2))
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models-root", default=str(OUTPUT_ROOT))
    parser.add_argument("--output-root", default=str(OUTPUT_ROOT / "synthetic_examples"))
    parser.add_argument("--models", nargs="+", default=list(DEFAULT_MODELS))
    parser.add_argument("--tissues", nargs="+", default=list(DEFAULT_TISSUES))
    parser.add_argument("--muscle-groups", nargs="+", default=list(MUSCLE_GROUPS))
    parser.add_argument("--n", type=int, default=16)
    parser.add_argument("--sample-steps", type=int, default=50)
    parser.add_argument("--eta", type=float, default=0.0)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--cpu", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
