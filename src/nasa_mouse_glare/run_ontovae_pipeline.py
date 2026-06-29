"""Run OntoVAE tissue workflows on API-derived OSDR and ARCHS4 inputs."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from . import analyze_expimap_pathways
from . import train_ontovae
from . import validate_expimap_accession_effects
from .io import require_import


DEFAULT_TISSUES = (
    "liver",
    "skeletal_muscle",
    "skin",
    "kidney",
    "thymus",
    "spleen",
    "lung",
    "retina",
)
MUSCLE_GROUPS = (
    "soleus",
    "gastrocnemius",
    "quadriceps",
    "edl",
    "tibialis_anterior",
)
SUMMARY_ROOT = Path("outputs/ontovae_pipeline/summary")


@dataclass(frozen=True)
class OntoVAERun:
    label: str
    tissue: str
    group: str
    mode: str
    query_h5ad: Path
    reference_h5ad: Path | None
    output_dir: Path

    @property
    def complete(self) -> bool:
        return (
            (self.output_dir / "training_summary.json").exists()
            and (self.output_dir / "analysis" / "analysis_summary.json").exists()
            and (self.output_dir / "accession_validation" / "validation_summary.json").exists()
        )


def direct_query_h5ad(tissue: str) -> Path:
    return (
        Path(f"outputs/expimap_direct_osdr_{tissue}")
        / "input"
        / f"osdr_{tissue}_flt_gc_reactome_raw_counts.h5ad"
    )


def reference_candidates(tissue: str) -> list[Path]:
    root = Path(f"outputs/expimap_archs4_reference_osdr_query_{tissue}")
    return [
        root / "reference_input_all" / f"archs4_mouse_{tissue}_reference_reactome_raw_counts.h5ad",
        root / "reference_input_5000_stratified" / f"archs4_mouse_{tissue}_reference_reactome_raw_counts.h5ad",
        root / "reference_input_1000" / f"archs4_mouse_{tissue}_reference_reactome_raw_counts.h5ad",
    ]


def first_existing(paths: list[Path]) -> Path | None:
    return next((path for path in paths if path.exists()), None)


def hvg_query_h5ad(tissue: str) -> Path:
    label = "5000" if tissue == "liver" else "2000"
    return (
        Path(f"outputs/expimap_archs4_reference_osdr_query_{tissue}")
        / f"tutorial_hvg_{label}"
        / "input"
        / f"osdr_{tissue}_query_tutorial_hvg_raw_counts.h5ad"
    )


def hvg_reference_h5ad(tissue: str) -> Path:
    label = "5000" if tissue == "liver" else "2000"
    return (
        Path(f"outputs/expimap_archs4_reference_osdr_query_{tissue}")
        / f"tutorial_hvg_{label}"
        / "input"
        / f"archs4_mouse_{tissue}_reference_tutorial_hvg_raw_counts.h5ad"
    )


def muscle_group_query_h5ad(group: str) -> Path:
    return (
        Path("outputs/expimap_muscle_targeted_combined_min8")
        / "group_inputs_exploratory_2acc"
        / f"osdr_skeletal_muscle_{group}_flt_gc_reactome_plus_muscle_raw_counts.h5ad"
    )


def muscle_group_reference_h5ad() -> Path:
    return (
        Path("outputs/expimap_muscle_targeted_combined_min8")
        / "reference_input_all"
        / "archs4_mouse_skeletal_muscle_reference_reactome_raw_counts.h5ad"
    )


def build_runs(args) -> list[OntoVAERun]:
    runs: list[OntoVAERun] = []
    for tissue in args.tissues:
        query = direct_query_h5ad(tissue)
        output_root = Path(f"outputs/ontovae_{tissue}")
        if args.include_direct and query.exists():
            runs.append(
                OntoVAERun(
                    label=f"{tissue}:direct",
                    tissue=tissue,
                    group="",
                    mode="direct_osdr",
                    query_h5ad=query,
                    reference_h5ad=None,
                    output_dir=output_root / "direct_osdr",
                )
            )
        reference = first_existing(reference_candidates(tissue))
        if args.include_reference and query.exists() and reference is not None:
            runs.append(
                OntoVAERun(
                    label=f"{tissue}:reference",
                    tissue=tissue,
                    group="",
                    mode="archs4_pretrain_osdr_finetune",
                    query_h5ad=query,
                    reference_h5ad=reference,
                    output_dir=output_root / "archs4_pretrain_osdr_finetune",
                )
            )
        hvg_query = hvg_query_h5ad(tissue)
        hvg_reference = hvg_reference_h5ad(tissue)
        if args.include_hvg and hvg_query.exists() and hvg_reference.exists():
            runs.append(
                OntoVAERun(
                    label=f"{tissue}:hvg_reference",
                    tissue=tissue,
                    group="",
                    mode="hvg_archs4_pretrain_osdr_finetune",
                    query_h5ad=hvg_query,
                    reference_h5ad=hvg_reference,
                    output_dir=output_root / "hvg_archs4_pretrain_osdr_finetune",
                )
            )

    if args.include_muscle_splits:
        reference = muscle_group_reference_h5ad()
        for group in args.muscle_groups:
            query = muscle_group_query_h5ad(group)
            output_root = Path("outputs/ontovae_skeletal_muscle_splits") / group
            if args.include_direct and query.exists():
                runs.append(
                    OntoVAERun(
                        label=f"skeletal_muscle:{group}:direct",
                        tissue="skeletal_muscle",
                        group=group,
                        mode="direct_osdr",
                        query_h5ad=query,
                        reference_h5ad=None,
                        output_dir=output_root / "direct_osdr",
                    )
                )
            if args.include_reference and query.exists() and reference.exists():
                runs.append(
                    OntoVAERun(
                        label=f"skeletal_muscle:{group}:reference",
                        tissue="skeletal_muscle",
                        group=group,
                        mode="archs4_pretrain_osdr_finetune",
                        query_h5ad=query,
                        reference_h5ad=reference,
                        output_dir=output_root / "archs4_pretrain_osdr_finetune",
                    )
                )
    return runs


def write_score_analyses(scores: Path, output_dir: Path) -> None:
    analyze_expimap_pathways.run(
        SimpleNamespace(
            scores=str(scores),
            output_dir=str(output_dir / "analysis"),
            include_de_novo=False,
        )
    )
    validate_expimap_accession_effects.run(
        SimpleNamespace(
            scores=str(scores),
            output_dir=str(output_dir / "accession_validation"),
        )
    )


def write_pretrained_score_analyses(scores: Path, output_dir: Path) -> None:
    if not scores.exists():
        return
    analyze_expimap_pathways.run(
        SimpleNamespace(
            scores=str(scores),
            output_dir=str(output_dir / "pretrained_query_analysis"),
            include_de_novo=False,
        )
    )
    validate_expimap_accession_effects.run(
        SimpleNamespace(
            scores=str(scores),
            output_dir=str(output_dir / "pretrained_query_accession_validation"),
        )
    )


def run_one(run: OntoVAERun, args) -> dict:
    if run.complete and not args.overwrite:
        return {
            "label": run.label,
            "tissue": run.tissue,
            "group": run.group,
            "mode": run.mode,
            "status": "skipped_complete",
            "output_dir": str(run.output_dir),
        }

    train_summary = train_ontovae.run(
        SimpleNamespace(
            query_h5ad=str(run.query_h5ad),
            reference_h5ad=str(run.reference_h5ad or ""),
            output_dir=str(run.output_dir),
            reference_epochs=args.reference_epochs,
            query_epochs=args.query_epochs,
            batch_size=args.batch_size,
            learning_rate=args.learning_rate,
            finetune_learning_rate=args.finetune_learning_rate,
            kl_coeff=args.kl_coeff,
            train_frac=args.train_frac,
            seed=args.seed,
            neuronnum=args.neuronnum,
            drop=args.drop,
            z_drop=args.z_drop,
            clip=args.clip,
            min_term_genes=args.min_term_genes,
            top_genes=args.top_genes,
            run_mode=run.mode,
        )
    )
    scores = run.output_dir / "pathway_scores.tsv"
    write_score_analyses(scores, run.output_dir)
    write_pretrained_score_analyses(
        run.output_dir / "pretrained_query_pathway_scores.tsv",
        run.output_dir,
    )
    try:
        torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    except SystemExit:
        pass

    return {
        "label": run.label,
        "tissue": run.tissue,
        "group": run.group,
        "mode": run.mode,
        "status": "completed",
        "query_h5ad": str(run.query_h5ad),
        "reference_h5ad": str(run.reference_h5ad or ""),
        "output_dir": str(run.output_dir),
        "training_summary": str(train_summary),
    }


def write_manifest(rows: list[dict], output_dir: Path) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "ontovae_run_manifest.tsv"
    pd.DataFrame(rows).to_csv(manifest_path, sep="\t", index=False)
    json_path = output_dir / "ontovae_run_manifest.json"
    json_path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def write_audit(runs: list[OntoVAERun], output_dir: Path) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "label": run.label,
            "tissue": run.tissue,
            "group": run.group,
            "mode": run.mode,
            "query_h5ad": str(run.query_h5ad),
            "query_exists": run.query_h5ad.exists(),
            "reference_h5ad": str(run.reference_h5ad or ""),
            "reference_exists": bool(run.reference_h5ad and run.reference_h5ad.exists()),
            "output_dir": str(run.output_dir),
            "complete": run.complete,
        }
        for run in runs
    ]
    path = output_dir / "ontovae_run_audit.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def run(args) -> Path:
    runs = build_runs(args)
    audit_path = write_audit(runs, Path(args.summary_dir))
    if args.dry_run:
        print(f"Wrote audit: {audit_path}")
        return audit_path

    rows = []
    for index, run in enumerate(runs, start=1):
        print(f"[{index}/{len(runs)}] {run.label} -> {run.output_dir}")
        rows.append(run_one(run, args))
        write_manifest(rows, Path(args.summary_dir))
    manifest_path = write_manifest(rows, Path(args.summary_dir))
    print(f"Wrote manifest: {manifest_path}")
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run OntoVAE per-tissue workflows.")
    parser.add_argument("--tissues", nargs="+", default=list(DEFAULT_TISSUES))
    parser.add_argument("--muscle-groups", nargs="+", default=list(MUSCLE_GROUPS))
    parser.add_argument("--summary-dir", default=str(SUMMARY_ROOT))
    parser.add_argument("--reference-epochs", type=int, default=60)
    parser.add_argument("--query-epochs", type=int, default=60)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--finetune-learning-rate", type=float, default=None)
    parser.add_argument("--kl-coeff", type=float, default=1e-4)
    parser.add_argument("--train-frac", type=float, default=0.9)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--neuronnum", type=int, default=1)
    parser.add_argument("--drop", type=float, default=0.1)
    parser.add_argument("--z-drop", type=float, default=0.1)
    parser.add_argument("--clip", type=float, default=10.0)
    parser.add_argument("--min-term-genes", type=int, default=5)
    parser.add_argument("--top-genes", type=int, default=30)
    parser.add_argument("--include-direct", action="store_true")
    parser.add_argument("--include-reference", action="store_true")
    parser.add_argument("--include-hvg", action="store_true")
    parser.add_argument("--include-muscle-splits", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not args.include_direct and not args.include_reference and not args.include_hvg:
        args.include_reference = True
    return args


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
