"""Run WGAN-GP direct and ARCHS4-pretrained workflows per tissue."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from . import analyze_features
from . import train_wgan
from .paths import DEFAULT_TISSUES, MUSCLE_GROUPS, build_run_specs


SUMMARY_ROOT = Path("outputs/wgan_pipeline/summary")


def run_complete(output_dir: Path) -> bool:
    return (
        (output_dir / "training_summary.json").exists()
        and (output_dir / "analysis" / "analysis_summary.json").exists()
    )


def write_analysis(output_dir: Path, *, top_features: int) -> None:
    scores = output_dir / "critic_feature_scores.tsv"
    analyze_features.run(
        SimpleNamespace(
            scores=str(scores),
            output_dir=str(output_dir / "analysis"),
            top_features=top_features,
        )
    )
    pretrained = output_dir / "pretrained_query_critic_feature_scores.tsv"
    if pretrained.exists():
        analyze_features.run(
            SimpleNamespace(
                scores=str(pretrained),
                output_dir=str(output_dir / "pretrained_query_analysis"),
                top_features=top_features,
            )
        )


def write_manifest(specs, output_dir: Path) -> Path:
    pd = __import__("pandas")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in specs:
        rows.append(
            {
                "label": spec.label,
                "tissue": spec.tissue,
                "group": spec.group,
                "mode": spec.mode,
                "query_h5ad": str(spec.query_h5ad),
                "reference_h5ad": str(spec.reference_h5ad or ""),
                "output_dir": str(spec.output_dir),
                "status": "completed" if run_complete(spec.output_dir) else "pending",
            }
        )
    path = output_dir / "wgan_run_manifest.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    (output_dir / "wgan_run_manifest.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )
    return path


def run(args) -> Path:
    output_dir = Path(args.summary_dir)
    specs = build_run_specs(
        tissues=tuple(args.tissues),
        muscle_groups=tuple(args.muscle_groups),
        include_direct=args.include_direct,
        include_reference=args.include_reference,
        include_muscle_splits=args.include_muscle_splits,
    )
    manifest = write_manifest(specs, output_dir)
    if args.dry_run:
        print(f"Wrote dry-run manifest: {manifest}")
        return manifest

    for index, spec in enumerate(specs, start=1):
        if run_complete(spec.output_dir) and not args.overwrite:
            print(f"[{index}/{len(specs)}] skip completed {spec.label}")
            continue
        print(f"[{index}/{len(specs)}] {spec.label} -> {spec.output_dir}")
        train_wgan.run(
            SimpleNamespace(
                query_h5ad=str(spec.query_h5ad),
                reference_h5ad=str(spec.reference_h5ad or ""),
                output_dir=str(spec.output_dir),
                run_mode=spec.mode,
                reference_epochs=args.reference_epochs,
                query_epochs=args.query_epochs,
                batch_size=args.batch_size,
                learning_rate=args.learning_rate,
                finetune_learning_rate=args.finetune_learning_rate,
                critic_steps=args.critic_steps,
                gradient_penalty=args.gradient_penalty,
                noise_dim=args.noise_dim,
                hidden_dims=args.hidden_dims,
                clip=args.clip,
                max_genes=args.max_genes,
                seed=args.seed,
                cpu=args.cpu,
            )
        )
        write_analysis(spec.output_dir, top_features=args.top_features)
        write_manifest(specs, output_dir)
    return write_manifest(specs, output_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tissues", nargs="+", default=list(DEFAULT_TISSUES))
    parser.add_argument("--muscle-groups", nargs="+", default=list(MUSCLE_GROUPS))
    parser.add_argument("--include-direct", action="store_true")
    parser.add_argument("--include-reference", action="store_true")
    parser.add_argument("--include-muscle-splits", action="store_true")
    parser.add_argument("--reference-epochs", type=int, default=100)
    parser.add_argument("--query-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--finetune-learning-rate", type=float, default=0.0)
    parser.add_argument("--critic-steps", type=int, default=5)
    parser.add_argument("--gradient-penalty", type=float, default=10.0)
    parser.add_argument("--noise-dim", type=int, default=128)
    parser.add_argument("--hidden-dims", default="256,256")
    parser.add_argument("--clip", type=float, default=10.0)
    parser.add_argument("--max-genes", type=int, default=0)
    parser.add_argument("--top-features", type=int, default=30)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--summary-dir", default=str(SUMMARY_ROOT))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--cpu", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
