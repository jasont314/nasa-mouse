"""Audit and run the multi-tissue expiMap variant matrix.

The matrix tracks five experiment variants per tissue:

1. direct OSDR training
2. ARCHS4 reference -> OSDR query
3. ARCHS4 reference -> OSDR query with de novo query nodes
4. HVG-filtered ARCHS4 reference -> OSDR query
5. HVG-filtered ARCHS4 reference -> OSDR query with de novo query nodes

By default this script writes an audit table and a shell script containing
commands for missing outputs. Use ``--execute`` to run missing steps
sequentially.
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path
import shlex
import subprocess
import sys


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
DEFAULT_OUTPUT_ROOT = Path("outputs/expimap_tissue_variant_matrix")
DEFAULT_PYTHON = Path("/home/exouser/miniforge3/envs/nasa-mouse/bin/python")
VARIANT_NAMES = (
    "direct",
    "reference_query",
    "reference_query_de_novo",
    "hvg_reference_query",
    "hvg_reference_query_de_novo",
)


@dataclass(frozen=True)
class Step:
    tissue: str
    variant: str
    name: str
    command: tuple[str, ...]
    outputs: tuple[Path, ...]
    kind: str = "core"

    @property
    def complete(self) -> bool:
        return all(path.exists() for path in self.outputs)

    @property
    def shell(self) -> str:
        return shlex.join(self.command)


@dataclass(frozen=True)
class Variant:
    tissue: str
    name: str
    output_dir: Path
    steps: tuple[Step, ...]

    @property
    def complete(self) -> bool:
        return all(step.complete for step in self.steps)

    @property
    def core_complete(self) -> bool:
        return all(step.complete for step in self.steps if step.kind == "core")

    @property
    def missing_steps(self) -> tuple[Step, ...]:
        return tuple(step for step in self.steps if not step.complete)

    @property
    def missing_core_steps(self) -> tuple[Step, ...]:
        return tuple(step for step in self.steps if step.kind == "core" and not step.complete)


def module_cmd(python: Path, module: str, *args: str | Path) -> tuple[str, ...]:
    return (str(python), "-m", f"nasa_mouse_glare.{module}", *(str(arg) for arg in args))


def direct_input(tissue: str) -> Path:
    return (
        Path(f"outputs/expimap_direct_osdr_{tissue}")
        / "input"
        / f"osdr_{tissue}_flt_gc_reactome_raw_counts.h5ad"
    )


def direct_run_dir(tissue: str) -> Path:
    return Path(f"outputs/expimap_direct_osdr_{tissue}") / "raw_counts_nb_50epoch"


def reference_root(tissue: str) -> Path:
    return Path(f"outputs/expimap_archs4_reference_osdr_query_{tissue}")


def canonical_reference_label(tissue: str) -> str:
    root = reference_root(tissue)
    if (root / "reference_nb_all_100epoch" / "training_summary.json").exists():
        return "all"
    if tissue == "liver" and (
        root / "reference_nb_5000_stratified_200epoch_seed2020" / "training_summary.json"
    ).exists():
        return "5000_stratified"
    return "1000"


def reference_input_dir(tissue: str, label: str) -> Path:
    root = reference_root(tissue)
    if label == "all":
        return root / "reference_input_all"
    if label == "5000_stratified":
        return root / "reference_input_5000_stratified"
    return root / "reference_input_1000"


def reference_h5ad(tissue: str, label: str) -> Path:
    return reference_input_dir(tissue, label) / f"archs4_mouse_{tissue}_reference_reactome_raw_counts.h5ad"


def reference_run_dir(tissue: str, label: str) -> Path:
    root = reference_root(tissue)
    if label == "all":
        return root / "reference_nb_all_100epoch"
    if label == "5000_stratified":
        return root / "reference_nb_5000_stratified_200epoch_seed2020"
    return root / "reference_nb_1000_50epoch"


def reference_query_run_dir(tissue: str, label: str) -> Path:
    root = reference_root(tissue)
    if label == "all":
        return root / "query_nb_allref_50epoch"
    if label == "5000_stratified":
        return root / "query_nb_5000stratified_seed2020_50epoch"
    return root / "query_nb_1000ref_50epoch"


def hvg_root(tissue: str) -> Path:
    if tissue == "liver":
        return reference_root(tissue) / "tutorial_hvg_5000"
    return reference_root(tissue) / "tutorial_hvg_2000"


def hvg_reference_h5ad(tissue: str) -> Path:
    return hvg_root(tissue) / "input" / f"archs4_mouse_{tissue}_reference_tutorial_hvg_raw_counts.h5ad"


def hvg_query_h5ad(tissue: str) -> Path:
    return hvg_root(tissue) / "input" / f"osdr_{tissue}_query_tutorial_hvg_raw_counts.h5ad"


def hvg_reference_run_dir(tissue: str) -> Path:
    return hvg_root(tissue) / "reference_nb_400epoch_seed2020"


def query_analysis_steps(
    *,
    python: Path,
    tissue: str,
    variant: str,
    query_dir: Path,
    include_de_novo: bool,
) -> list[Step]:
    scores = query_dir / "query_pathway_scores.tsv"
    analysis = query_dir / "analysis"
    steps = [
        Step(
            tissue,
            variant,
            "analyze_query_scores",
            module_cmd(
                python,
                "analyze_expimap_pathways",
                "--scores",
                scores,
                "--output-dir",
                analysis,
                *(("--include-de-novo",) if include_de_novo else ()),
            ),
            (analysis / "analysis_summary.json",),
        ),
        Step(
            tissue,
            variant,
            "validate_accession_effects",
            module_cmd(
                python,
                "validate_expimap_accession_effects",
                "--scores",
                scores,
                "--output-dir",
                query_dir / "accession_validation",
            ),
            (query_dir / "accession_validation" / "validation_summary.json",),
        ),
        Step(
            tissue,
            variant,
            "latent_enrich_condition",
            module_cmd(
                python,
                "run_expimap_latent_enrich_condition",
                "--model-dir",
                query_dir / "query_model",
                "--h5ad",
                query_dir / "mapped_query_with_scores.h5ad",
                "--output-dir",
                query_dir / "latent_enrich_condition",
            ),
            (query_dir / "latent_enrich_condition" / "latent_enrich_summary.json",),
            "optional_bf",
        ),
    ]
    if include_de_novo:
        steps.append(
            Step(
                tissue,
                variant,
                "summarize_de_novo_programs",
                module_cmd(
                    python,
                    "summarize_expimap_de_novo",
                    "--mapped-h5ad",
                    query_dir / "mapped_query_with_scores.h5ad",
                    "--scores",
                    scores,
                    "--comparison",
                    analysis / "flt_vs_gc_pathway_comparison.tsv",
                    "--study-tests",
                    analysis / "flight_ground_study_aware_tests.tsv",
                    "--programs",
                    query_dir / "de_novo_programs.tsv",
                    "--gene-loadings",
                    query_dir / "de_novo_program_gene_loadings.tsv",
                    "--output-dir",
                    query_dir / "de_novo_analysis",
                ),
                (query_dir / "de_novo_analysis" / "de_novo_program_summary.tsv",),
            )
        )
    return steps


def direct_variant(python: Path, tissue: str) -> Variant:
    output_dir = direct_run_dir(tissue)
    steps = [
        Step(
            tissue,
            "direct",
            "prepare_osdr_input",
            module_cmd(
                python,
                "prepare_expimap_osdr_tissue",
                "--tissue",
                tissue,
                "--transform",
                "raw_counts",
                "--output-dir",
                output_dir.parent / "input",
            ),
            (direct_input(tissue),),
        ),
        Step(
            tissue,
            "direct",
            "train_direct",
            module_cmd(
                python,
                "train_expimap_direct",
                "--input",
                direct_input(tissue),
                "--output-dir",
                output_dir,
                "--recon-loss",
                "nb",
                "--epochs",
                "50",
            ),
            (output_dir / "training_summary.json", output_dir / "pathway_scores.tsv"),
        ),
        Step(
            tissue,
            "direct",
            "analyze_direct_scores",
            module_cmd(
                python,
                "analyze_expimap_pathways",
                "--scores",
                output_dir / "pathway_scores.tsv",
                "--output-dir",
                output_dir / "analysis",
            ),
            (output_dir / "analysis" / "analysis_summary.json",),
        ),
        Step(
            tissue,
            "direct",
            "validate_accession_effects",
            module_cmd(
                python,
                "validate_expimap_accession_effects",
                "--scores",
                output_dir / "pathway_scores.tsv",
                "--output-dir",
                output_dir / "accession_validation",
            ),
            (output_dir / "accession_validation" / "validation_summary.json",),
        ),
    ]
    return Variant(tissue, "direct", output_dir, tuple(steps))


def reference_query_variant(python: Path, tissue: str) -> Variant:
    label = canonical_reference_label(tissue)
    query_dir = reference_query_run_dir(tissue, label)
    steps = query_analysis_steps(
        python=python,
        tissue=tissue,
        variant="reference_query",
        query_dir=query_dir,
        include_de_novo=False,
    )
    return Variant(tissue, "reference_query", query_dir, tuple(steps))


def reference_query_denovo_variant(python: Path, tissue: str) -> Variant:
    label = canonical_reference_label(tissue)
    query_dir = reference_root(tissue) / "query_denovo3_hsic_250epoch_seed2020"
    steps = [
        Step(
            tissue,
            "reference_query_de_novo",
            "map_query_de_novo",
            module_cmd(
                python,
                "map_expimap_osdr_query",
                "--reference-model",
                reference_run_dir(tissue, label) / "model",
                "--query-h5ad",
                direct_input(tissue),
                "--output-dir",
                query_dir,
                "--epochs",
                "250",
                "--alpha-kl",
                "0.22",
                "--n-de-novo-programs",
                "3",
                "--gamma-ext",
                "0.7",
                "--gamma-epoch-anneal",
                "50",
                "--use-hsic",
                "--hsic-one-vs-all",
                "--no-alpha",
            ),
            (query_dir / "query_mapping_summary.json", query_dir / "de_novo_programs.tsv"),
        )
    ]
    steps.extend(
        query_analysis_steps(
            python=python,
            tissue=tissue,
            variant="reference_query_de_novo",
            query_dir=query_dir,
            include_de_novo=True,
        )
    )
    return Variant(tissue, "reference_query_de_novo", query_dir, tuple(steps))


def hvg_input_steps(python: Path, tissue: str) -> Step:
    label = canonical_reference_label(tissue)
    return Step(
        tissue,
        "hvg_shared",
        "prepare_hvg_inputs",
        module_cmd(
            python,
            "prepare_expimap_tutorial_hvg",
            "--reference-h5ad",
            reference_h5ad(tissue, label),
            "--query-h5ad",
            direct_input(tissue),
            "--output-dir",
            hvg_root(tissue) / "input",
            "--n-top-genes",
            "2000",
            "--allow-no-batch-fallback",
        ),
        (hvg_reference_h5ad(tissue), hvg_query_h5ad(tissue)),
    )


def hvg_reference_step(python: Path, tissue: str, variant: str) -> Step:
    output_dir = hvg_reference_run_dir(tissue)
    return Step(
        tissue,
        variant,
        "train_hvg_reference",
        module_cmd(
            python,
            "train_expimap_archs4_reference",
            "--input",
            hvg_reference_h5ad(tissue),
            "--output-dir",
            output_dir,
            "--recon-loss",
            "nb",
            "--epochs",
            "400",
            "--hidden-layers",
            "300,300,300",
            "--alpha-kl",
            "0.5",
            "--alpha-epoch-anneal",
            "100",
            "--early-stopping",
            "--early-stopping-patience",
            "50",
        ),
        (output_dir / "training_summary.json", output_dir / "model" / "model_params.pt"),
    )


def hvg_query_variant(python: Path, tissue: str) -> Variant:
    query_dir = hvg_root(tissue) / "query_nb_250epoch_seed2020"
    steps = [
        hvg_input_steps(python, tissue),
        hvg_reference_step(python, tissue, "hvg_reference_query"),
        Step(
            tissue,
            "hvg_reference_query",
            "map_hvg_query",
            module_cmd(
                python,
                "map_expimap_osdr_query",
                "--reference-model",
                hvg_reference_run_dir(tissue) / "model",
                "--query-h5ad",
                hvg_query_h5ad(tissue),
                "--output-dir",
                query_dir,
                "--epochs",
                "250",
                "--alpha-kl",
                "0.22",
                "--no-alpha",
            ),
            (query_dir / "query_mapping_summary.json",),
        ),
    ]
    steps.extend(
        query_analysis_steps(
            python=python,
            tissue=tissue,
            variant="hvg_reference_query",
            query_dir=query_dir,
            include_de_novo=False,
        )
    )
    return Variant(tissue, "hvg_reference_query", query_dir, tuple(steps))


def hvg_denovo_variant(python: Path, tissue: str) -> Variant:
    query_dir = hvg_root(tissue) / "query_denovo3_hsic_250epoch_seed2020"
    steps = [
        hvg_input_steps(python, tissue),
        hvg_reference_step(python, tissue, "hvg_reference_query_de_novo"),
        Step(
            tissue,
            "hvg_reference_query_de_novo",
            "map_hvg_query_de_novo",
            module_cmd(
                python,
                "map_expimap_osdr_query",
                "--reference-model",
                hvg_reference_run_dir(tissue) / "model",
                "--query-h5ad",
                hvg_query_h5ad(tissue),
                "--output-dir",
                query_dir,
                "--epochs",
                "250",
                "--alpha-kl",
                "0.22",
                "--n-de-novo-programs",
                "3",
                "--gamma-ext",
                "0.7",
                "--gamma-epoch-anneal",
                "50",
                "--use-hsic",
                "--hsic-one-vs-all",
                "--no-alpha",
            ),
            (query_dir / "query_mapping_summary.json", query_dir / "de_novo_programs.tsv"),
        ),
    ]
    steps.extend(
        query_analysis_steps(
            python=python,
            tissue=tissue,
            variant="hvg_reference_query_de_novo",
            query_dir=query_dir,
            include_de_novo=True,
        )
    )
    return Variant(tissue, "hvg_reference_query_de_novo", query_dir, tuple(steps))


def build_variants(python: Path, tissues: list[str], variant_filter: set[str] | None) -> list[Variant]:
    variants = []
    for tissue in tissues:
        tissue_variants = [
            direct_variant(python, tissue),
            reference_query_variant(python, tissue),
            reference_query_denovo_variant(python, tissue),
            hvg_query_variant(python, tissue),
            hvg_denovo_variant(python, tissue),
        ]
        if variant_filter is not None:
            tissue_variants = [
                variant for variant in tissue_variants if variant.name in variant_filter
            ]
        variants.extend(tissue_variants)
    return variants


def write_status(variants: list[Variant], output_root: Path) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    status_path = output_root / "variant_matrix_status.tsv"
    with status_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tissue",
                "variant",
                "complete",
                "core_complete",
                "output_dir",
                "n_steps",
                "n_missing_steps",
                "n_missing_core_steps",
                "missing_steps",
                "missing_core_steps",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for variant in variants:
            writer.writerow(
                {
                    "tissue": variant.tissue,
                    "variant": variant.name,
                    "complete": int(variant.complete),
                    "core_complete": int(variant.core_complete),
                    "output_dir": str(variant.output_dir),
                    "n_steps": len(variant.steps),
                    "n_missing_steps": len(variant.missing_steps),
                    "n_missing_core_steps": len(variant.missing_core_steps),
                    "missing_steps": ";".join(step.name for step in variant.missing_steps),
                    "missing_core_steps": ";".join(
                        step.name for step in variant.missing_core_steps
                    ),
                }
            )
    return status_path


def missing_steps(variants: list[Variant], core_only: bool) -> list[Step]:
    seen_outputs: set[tuple[str, ...]] = set()
    steps: list[Step] = []
    for variant in variants:
        for step in variant.missing_steps:
            if core_only and step.kind != "core":
                continue
            key = tuple(str(path) for path in step.outputs)
            if key in seen_outputs:
                continue
            seen_outputs.add(key)
            steps.append(step)
    return steps


def write_missing_shell(steps: list[Step], output_root: Path) -> Path:
    script_path = output_root / "run_missing_matrix.sh"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"',
        "",
    ]
    for step in steps:
        lines.append(f"echo '[{step.tissue} / {step.variant}] {step.name}'")
        lines.append(step.shell)
        lines.append("")
    script_path.write_text("\n".join(lines), encoding="utf-8")
    script_path.chmod(0o755)
    return script_path


def write_readme(variants: list[Variant], steps: list[Step], output_root: Path) -> Path:
    complete = sum(1 for variant in variants if variant.complete)
    core_complete = sum(1 for variant in variants if variant.core_complete)
    by_variant: dict[str, tuple[int, int]] = {}
    by_variant_core: dict[str, tuple[int, int]] = {}
    for name in sorted({variant.name for variant in variants}):
        subset = [variant for variant in variants if variant.name == name]
        by_variant[name] = (sum(1 for variant in subset if variant.complete), len(subset))
        by_variant_core[name] = (
            sum(1 for variant in subset if variant.core_complete),
            len(subset),
        )

    lines = [
        "# expiMap Tissue Variant Matrix",
        "",
        "Canonical variants:",
        "",
        "- `direct`: OSDR FLT/GC direct expiMap on retained Reactome genes.",
        "- `reference_query`: ARCHS4 tissue reference mapped to OSDR query.",
        "- `reference_query_de_novo`: ARCHS4 reference-query plus three HSIC-regularized de novo query nodes.",
        "- `hvg_reference_query`: tutorial-style 2,000-HVG reference-query.",
        "- `hvg_reference_query_de_novo`: tutorial-style 2,000-HVG reference-query plus three HSIC-regularized de novo nodes.",
        "",
        f"Core-complete variants: {core_complete}/{len(variants)}",
        f"Full-complete variants, including optional latent enrichment: {complete}/{len(variants)}",
        f"Missing executable steps: {len(steps)}",
        "",
        "Core status by variant:",
        "",
    ]
    for name, (done, total) in by_variant_core.items():
        lines.append(f"- `{name}`: {done}/{total}")
    lines.extend(["", "Full status by variant, including optional latent enrichment:", ""])
    for name, (done, total) in by_variant.items():
        lines.append(f"- `{name}`: {done}/{total}")
    lines.extend(
        [
            "",
            "Files:",
            "",
            "- `variant_matrix_status.tsv`: one row per tissue/variant.",
            "- `run_missing_matrix.sh`: sequential commands for incomplete outputs.",
            "",
            "Run the audit only:",
            "",
            "```bash",
            "python3 src/nasa_mouse_glare/run_expimap_tissue_variant_matrix.py",
            "```",
            "",
            "Run missing steps deliberately:",
            "",
            "```bash",
            "PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.run_expimap_tissue_variant_matrix --execute",
            "```",
            "",
            "Use `--max-steps N` to submit a bounded chunk of work.",
            "",
            "Run only one core variant block, for example all-gene de novo query mapping:",
            "",
            "```bash",
            "PYTHONPATH=src /home/exouser/miniforge3/envs/nasa-mouse/bin/python -m nasa_mouse_glare.run_expimap_tissue_variant_matrix --variants reference_query_de_novo --core-only --execute",
            "```",
            "",
        ]
    )
    readme_path = output_root / "README.md"
    readme_path.write_text("\n".join(lines), encoding="utf-8")
    return readme_path


def execute_steps(steps: list[Step], max_steps: int | None) -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = "src" + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    selected = steps if max_steps is None else steps[:max_steps]
    for index, step in enumerate(selected, start=1):
        print(f"[{index}/{len(selected)}] {step.tissue} / {step.variant}: {step.name}", flush=True)
        subprocess.run(step.command, check=True, env=env)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tissues", nargs="+", default=list(DEFAULT_TISSUES))
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--variants", nargs="+", choices=VARIANT_NAMES, default=None)
    parser.add_argument(
        "--core-only",
        action="store_true",
        help="Generate or execute only core model, analysis, validation, and de novo-summary steps.",
    )
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--max-steps", type=int, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    python = args.python if args.python.exists() else Path(sys.executable)
    tissues = [tissue.strip().lower().replace(" ", "_") for tissue in args.tissues]
    variant_filter = set(args.variants) if args.variants else None
    variants = build_variants(python, tissues, variant_filter)
    steps = missing_steps(variants, args.core_only)
    status_path = write_status(variants, args.output_root)
    script_path = write_missing_shell(steps, args.output_root)
    readme_path = write_readme(variants, steps, args.output_root)
    print(f"Wrote {status_path}")
    print(f"Wrote {script_path}")
    print(f"Wrote {readme_path}")
    print(
        "Core-complete variants: "
        f"{sum(1 for variant in variants if variant.core_complete)}/{len(variants)}"
    )
    print(f"Full-complete variants: {sum(1 for variant in variants if variant.complete)}/{len(variants)}")
    print(f"Missing steps: {len(steps)}")
    if args.execute:
        execute_steps(steps, args.max_steps)


if __name__ == "__main__":
    main()
