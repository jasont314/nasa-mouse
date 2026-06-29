"""Train pan-tissue conditional diffusion models for synthetic generation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from . import analyze_features
from . import train_diffusion
from .paths import DEFAULT_TISSUES, OUTPUT_ROOT, SUMMARY_DIR, existing_pan_tissue_inputs


def joined(paths: list[Path]) -> str:
    return ",".join(str(path) for path in paths)


def run_complete(output_dir: Path, *, expect_analysis: bool) -> bool:
    if not (output_dir / "training_summary.json").exists():
        return False
    if expect_analysis and not (output_dir / "analysis" / "analysis_summary.json").exists():
        return False
    return True


def analyze_if_osdr(output_dir: Path, *, top_features: int) -> None:
    analyze_features.run(
        SimpleNamespace(
            scores=str(output_dir / "diffusion_feature_scores.tsv"),
            output_dir=str(output_dir / "analysis"),
            top_features=top_features,
        )
    )
    pretrained = output_dir / "pretrained_query_diffusion_feature_scores.tsv"
    if pretrained.exists():
        analyze_features.run(
            SimpleNamespace(
                scores=str(pretrained),
                output_dir=str(output_dir / "pretrained_query_analysis"),
                top_features=top_features,
            )
        )


def train_one(args, *, label: str, query_h5ad: str, query_source: str, output_dir: Path, reference_h5ad: str = "", expect_analysis: bool = True):
    if run_complete(output_dir, expect_analysis=expect_analysis) and not args.overwrite:
        return {"label": label, "status": "completed", "output_dir": str(output_dir)}
    output_dir.mkdir(parents=True, exist_ok=True)
    train_diffusion.run(
        SimpleNamespace(
            query_h5ad=query_h5ad,
            reference_h5ad=reference_h5ad,
            query_source=query_source,
            reference_source="archs4",
            output_dir=str(output_dir),
            run_mode=label,
            reference_epochs=args.reference_epochs,
            query_epochs=args.query_epochs,
            batch_size=args.batch_size,
            learning_rate=args.learning_rate,
            finetune_learning_rate=args.finetune_learning_rate,
            hidden_dims=args.hidden_dims,
            dropout=args.dropout,
            time_embedding_dim=args.time_embedding_dim,
            categorical_embedding_dim=args.categorical_embedding_dim,
            sinusoidal_time=args.sinusoidal_time,
            diffusion_timesteps=args.diffusion_timesteps,
            sample_steps=args.sample_steps,
            eta=args.eta,
            beta_schedule=args.beta_schedule,
            beta_start=args.beta_start,
            beta_end=args.beta_end,
            n_landmarks=args.n_landmarks,
            landmark_strategy=args.landmark_strategy,
            l1000_map=args.l1000_map,
            reconstruction_alpha=args.reconstruction_alpha,
            max_genes=args.max_genes,
            max_metric_samples=args.max_metric_samples,
            seed=args.seed,
            cpu=args.cpu,
            no_amp=args.no_amp,
        )
    )
    if expect_analysis:
        analyze_if_osdr(output_dir, top_features=args.top_features)
    return {"label": label, "status": "completed", "output_dir": str(output_dir)}


def write_manifest(rows: list[dict]) -> Path:
    pd = __import__("pandas")
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    path = SUMMARY_DIR / "conditional_generation_manifest.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    (SUMMARY_DIR / "conditional_generation_manifest.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    return path


def run(args) -> Path:
    query_paths, reference_paths = existing_pan_tissue_inputs(tuple(args.tissues))
    if not query_paths:
        raise SystemExit("No OSDR query h5ad inputs found.")
    if not reference_paths:
        raise SystemExit("No ARCHS4 reference h5ad inputs found.")
    rows = [
        {
            "label": "conditional_osdr_only",
            "query_h5ad": joined(query_paths),
            "reference_h5ad": "",
            "output_dir": str(OUTPUT_ROOT / "osdr_only"),
            "status": "planned",
        },
        {
            "label": "conditional_archs4_pretrain_osdr_finetune",
            "query_h5ad": joined(query_paths),
            "reference_h5ad": joined(reference_paths),
            "output_dir": str(OUTPUT_ROOT / "archs4_pretrain_osdr_finetune"),
            "status": "planned",
        },
        {
            "label": "conditional_archs4_only",
            "query_h5ad": joined(reference_paths),
            "reference_h5ad": "",
            "output_dir": str(OUTPUT_ROOT / "archs4_only"),
            "status": "planned",
        },
    ]
    if args.dry_run:
        return write_manifest(rows)
    results = []
    if args.include_osdr_only:
        results.append(
            train_one(
                args,
                label="conditional_osdr_only",
                query_h5ad=joined(query_paths),
                query_source="osdr",
                output_dir=OUTPUT_ROOT / "osdr_only",
                expect_analysis=True,
            )
        )
    if args.include_archs4_pretrain:
        results.append(
            train_one(
                args,
                label="conditional_archs4_pretrain_osdr_finetune",
                query_h5ad=joined(query_paths),
                query_source="osdr",
                reference_h5ad=joined(reference_paths),
                output_dir=OUTPUT_ROOT / "archs4_pretrain_osdr_finetune",
                expect_analysis=True,
            )
        )
    if args.include_archs4_only:
        results.append(
            train_one(
                args,
                label="conditional_archs4_only",
                query_h5ad=joined(reference_paths),
                query_source="archs4",
                output_dir=OUTPUT_ROOT / "archs4_only",
                expect_analysis=False,
            )
        )
    return write_manifest(results)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tissues", nargs="+", default=list(DEFAULT_TISSUES))
    parser.add_argument("--include-osdr-only", action="store_true", default=True)
    parser.add_argument("--include-archs4-pretrain", action="store_true", default=True)
    parser.add_argument("--include-archs4-only", action="store_true", default=True)
    parser.add_argument("--reference-epochs", type=int, default=100)
    parser.add_argument("--query-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--learning-rate", type=float, default=5e-4)
    parser.add_argument("--finetune-learning-rate", type=float, default=0.0)
    parser.add_argument("--hidden-dims", default="512,512")
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--time-embedding-dim", type=int, default=64)
    parser.add_argument("--categorical-embedding-dim", type=int, default=8)
    parser.add_argument("--sinusoidal-time", action="store_true")
    parser.add_argument("--diffusion-timesteps", type=int, default=1000)
    parser.add_argument("--sample-steps", type=int, default=50)
    parser.add_argument("--eta", type=float, default=0.0)
    parser.add_argument("--beta-schedule", default="quad")
    parser.add_argument("--beta-start", type=float, default=0.0001)
    parser.add_argument("--beta-end", type=float, default=0.02)
    parser.add_argument("--n-landmarks", type=int, default=512)
    parser.add_argument("--landmark-strategy", choices=["l1000", "l1000_or_hvg", "hvg"], default="l1000_or_hvg")
    parser.add_argument("--l1000-map", default="data/diffusion/l1000_human_to_mouse_ensembl.tsv")
    parser.add_argument("--reconstruction-alpha", type=float, default=1.0)
    parser.add_argument("--max-genes", type=int, default=0)
    parser.add_argument("--max-metric-samples", type=int, default=2000)
    parser.add_argument("--top-features", type=int, default=30)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--cpu", action="store_true")
    parser.add_argument("--no-amp", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
