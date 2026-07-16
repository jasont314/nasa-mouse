"""Build a ranked scoreboard from completed unified and exact baseline runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _nested(payload: dict[str, Any], *keys: str, default: Any = np.nan) -> Any:
    value: Any = payload
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def _unified_row(summary_path: Path) -> dict[str, Any]:
    run = json.loads(summary_path.read_text(encoding="utf-8"))
    validation_path = str(_nested(run, "outputs", "validation", default=""))
    validation_file = Path(validation_path) if validation_path else None
    validation = (
        json.loads(validation_file.read_text(encoding="utf-8"))
        if validation_file is not None and validation_file.exists()
        else {}
    )
    generation = validation.get("generation", {})
    selection = generation.get("model_selection", {})
    fidelity = generation.get("fidelity_transformed", {})
    effect = generation.get("flt_gc_effect_recovery", {})
    utilities = generation.get("expression_flt_gc_utility_by_ratio", {})
    ratio_one = utilities.get("ratio_1", {})
    augmented = ratio_one.get("real_plus_synthetic_train_real_evaluation", {})
    real = ratio_one.get("real_train_real_evaluation", {})
    return {
        "run_id": run.get("run_id", summary_path.parent.name),
        "model": run.get("model", ""),
        "implementation": _nested(
            run, "model_provenance", "display_name", default=""
        ),
        "model_profile": run.get("model_profile", ""),
        "regime": run.get("regime", ""),
        "tissue_mode": run.get("tissue_mode", ""),
        "tissues": ";".join(map(str, _nested(run, "data", "tissues", default=[]))),
        "genes": run.get("genes", np.nan),
        "reference_samples": _nested(run, "data", "reference_samples"),
        "train_samples": _nested(run, "data", "partition_samples", "train"),
        "training_seconds": run.get("training_seconds", np.nan),
        "cuda_peak_memory_gb": run.get("cuda_peak_memory_gb", np.nan),
        "generation_status": generation.get("status", "missing"),
        "heldout_fidelity_composite": selection.get(
            "heldout_fidelity_composite", np.nan
        ),
        "eligible_for_model_selection": selection.get(
            "eligible_for_model_selection", False
        ),
        "diversity_gate": _nested(selection, "diversity_gate", "passed", default=False),
        "memorization_gate": _nested(
            selection, "memorization_gate", "passed", default=False
        ),
        "gene_mean_correlation": fidelity.get("gene_mean_correlation", np.nan),
        "gene_std_correlation": fidelity.get("gene_std_correlation", np.nan),
        "precision": fidelity.get("precision", np.nan),
        "recall": fidelity.get("recall", np.nan),
        "adversarial_accuracy": fidelity.get("adversarial_accuracy", np.nan),
        "flt_gc_delta_correlation": effect.get("delta_correlation", np.nan),
        "real_validation_balanced_accuracy": real.get("balanced_accuracy", np.nan),
        "augmented_validation_balanced_accuracy": augmented.get(
            "balanced_accuracy", np.nan
        ),
        "augmentation_status": ratio_one.get("augmentation_status", ""),
        "summary_path": str(summary_path),
    }


def _exact_diffusion_row(summary_path: Path) -> dict[str, Any]:
    run = json.loads(summary_path.read_text(encoding="utf-8"))
    evaluation_path = summary_path.parent / "evaluation/summary.json"
    evaluation = (
        json.loads(evaluation_path.read_text(encoding="utf-8"))
        if evaluation_path.exists()
        else {}
    )
    quality = evaluation.get("quality", evaluation)
    direct_pr = quality.get("precision_recall_in_scaled_l974", {})
    precision = float(direct_pr.get("precision", np.nan))
    recall = float(direct_pr.get("recall", np.nan))
    f1 = (
        2.0 * precision * recall / (precision + recall)
        if np.isfinite(precision + recall) and precision + recall > 0
        else np.nan
    )
    adversarial = float(
        quality.get(
            "nearest_neighbor_adversarial_accuracy_in_scaled_l974", np.nan
        )
    )
    adversarial_score = (
        float(np.clip(1.0 - 2.0 * abs(adversarial - 0.5), 0.0, 1.0))
        if np.isfinite(adversarial)
        else np.nan
    )
    components = [
        float(quality.get("gene_mean_correlation", np.nan)),
        float(quality.get("gene_standard_deviation_correlation", np.nan)),
        f1,
        adversarial_score,
    ]
    legacy_composite = (
        float(np.mean(components)) if np.isfinite(components).all() else np.nan
    )
    selection = quality.get("model_selection", {})
    composite = selection.get("heldout_fidelity_composite", legacy_composite)
    eligible = bool(selection.get("eligible_for_model_selection", False))
    return {
        "run_id": summary_path.parent.name,
        "model": "lacan_diffusion",
        "implementation": "upstream ModelDDIM paper-parity",
        "model_profile": "paper_native",
        "regime": "archs4_only",
        "tissue_mode": "pooled_conditioned",
        "tissues": ";".join(map(str, run.get("classes", []))),
        "genes": 974,
        "reference_samples": 0,
        "train_samples": _nested(run, "profiles", "train"),
        "training_seconds": run.get("training_seconds_this_invocation", np.nan),
        "cuda_peak_memory_gb": run.get("cuda_peak_memory_gb", np.nan),
        "generation_status": "complete" if evaluation else "missing",
        "heldout_fidelity_composite": composite,
        "eligible_for_model_selection": eligible,
        "diversity_gate": bool(
            _nested(selection, "diversity_gate", "passed", default=False)
        ),
        "memorization_gate": bool(
            _nested(selection, "memorization_gate", "passed", default=False)
        ),
        "gene_mean_correlation": quality.get("gene_mean_correlation", np.nan),
        "gene_std_correlation": quality.get(
            "gene_standard_deviation_correlation", np.nan
        ),
        "precision": precision,
        "recall": recall,
        "adversarial_accuracy": adversarial,
        "flt_gc_delta_correlation": np.nan,
        "real_validation_balanced_accuracy": np.nan,
        "augmented_validation_balanced_accuracy": np.nan,
        "augmentation_status": "not_applicable_archs4_tissue_baseline",
        "summary_path": str(summary_path),
    }


def _conditional_diffusion_row(summary_path: Path) -> dict[str, Any]:
    run = json.loads(summary_path.read_text(encoding="utf-8"))
    evaluation_path = summary_path.parent / "evaluation/validation/summary.json"
    evaluation = (
        json.loads(evaluation_path.read_text(encoding="utf-8"))
        if evaluation_path.exists()
        else {}
    )
    fidelity = evaluation.get("fidelity_transformed", {})
    selection = evaluation.get("model_selection", {})
    utility = evaluation.get("flt_gc_classifier_utility", {})
    real_utility = utility.get("real_train_real_evaluation", {})
    augmented_utility = utility.get(
        "real_plus_synthetic_train_real_evaluation", {}
    )
    classes = list(map(str, run.get("classes", [])))
    tissues = sorted(
        {
            field.split("=", 1)[1]
            for label in classes
            for field in label.split("||")
            if field.startswith("tissue=")
        }
    )
    regime = str(run.get("regime", "osdr_only"))
    return {
        "run_id": summary_path.parent.name,
        "model": "lacan_diffusion",
        "implementation": "upstream ModelDDIM NASA condition extension",
        "model_profile": "paper_architecture_osdr_extension",
        "regime": regime,
        "tissue_mode": "pooled_conditioned",
        "tissues": ";".join(tissues),
        "genes": 974,
        "reference_samples": 9796 if regime.startswith("archs4_pretrain") else 0,
        "train_samples": _nested(run, "profiles", "train"),
        "training_seconds": run.get("training_seconds_this_invocation", np.nan),
        "cuda_peak_memory_gb": run.get("cuda_peak_memory_gb", np.nan),
        "generation_status": "complete" if evaluation else "missing",
        "heldout_fidelity_composite": selection.get(
            "heldout_fidelity_composite", np.nan
        ),
        "eligible_for_model_selection": selection.get(
            "eligible_for_model_selection", False
        ),
        "diversity_gate": _nested(
            selection, "diversity_gate", "passed", default=False
        ),
        "memorization_gate": _nested(
            selection, "memorization_gate", "passed", default=False
        ),
        "gene_mean_correlation": fidelity.get("gene_mean_correlation", np.nan),
        "gene_std_correlation": fidelity.get("gene_std_correlation", np.nan),
        "precision": fidelity.get("precision", np.nan),
        "recall": fidelity.get("recall", np.nan),
        "adversarial_accuracy": fidelity.get("adversarial_accuracy", np.nan),
        "flt_gc_delta_correlation": _nested(
            evaluation, "flt_gc_effect_recovery", "delta_correlation"
        ),
        "real_validation_balanced_accuracy": real_utility.get(
            "balanced_accuracy", np.nan
        ),
        "augmented_validation_balanced_accuracy": augmented_utility.get(
            "balanced_accuracy", np.nan
        ),
        "augmentation_status": utility.get("augmentation_status", ""),
        "summary_path": str(summary_path),
    }


def build_scoreboard(output_root: str | Path) -> pd.DataFrame:
    root = Path(output_root)
    rows: list[dict[str, Any]] = []
    for path in sorted((root / "runs").glob("*/*/run_summary.json")):
        if path.parent.name == "archs4_mouse_paper_parity_seed1234":
            rows.append(_exact_diffusion_row(path))
        elif "ModelDDIM OSDR condition extension" in json.loads(
            path.read_text(encoding="utf-8")
        ).get("model", ""):
            rows.append(_conditional_diffusion_row(path))
        else:
            rows.append(_unified_row(path))
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    table["rank"] = table["heldout_fidelity_composite"].where(
        table["eligible_for_model_selection"].astype(bool)
    ).rank(method="min", ascending=False)
    return table.sort_values(
        ["eligible_for_model_selection", "heldout_fidelity_composite"],
        ascending=[False, False],
        na_position="last",
    ).reset_index(drop=True)


def run(args: argparse.Namespace) -> Path:
    table = build_scoreboard(args.output_root)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index=False)
    summary = {
        "runs": int(len(table)),
        "generation_complete": int(
            table.get("generation_status", pd.Series(dtype=str)).eq("complete").sum()
        ),
        "eligible": int(
            table.get("eligible_for_model_selection", pd.Series(dtype=bool))
            .astype(bool)
            .sum()
        ),
        "scoreboard": str(output),
    }
    output.with_suffix(".json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="outputs/generative_benchmark")
    parser.add_argument(
        "--output",
        default="outputs/generative_benchmark/summary/model_scoreboard.tsv",
    )
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
