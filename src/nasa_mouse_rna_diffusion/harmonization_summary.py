"""Summarize matched OSDR DDIM harmonization runs without composite scoring."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


FIDELITY_METRICS = (
    ("heldout_corr", "Corr", 0.98, "minimum"),
    ("heldout_precision", "Precision", 0.95, "minimum"),
    ("heldout_recall", "Recall", 0.85, "minimum"),
    ("heldout_f1", "F1", 0.90, "minimum"),
    ("heldout_aa", "Adversarial accuracy", (0.40, 0.60), "range"),
    ("heldout_fd_ratio", "FD / real-split p95", 1.0, "maximum"),
)


def _get(payload: dict[str, Any], *keys: str, default: Any = np.nan) -> Any:
    value: Any = payload
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _record(label: str, config_path: Path, split: str) -> dict[str, Any]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    run_dir = Path(config["run"]["output_dir"])
    summary_path = run_dir / "evaluation" / split / "summary.json"
    record: dict[str, Any] = {
        "label": label,
        "config": str(config_path),
        "run_dir": str(run_dir),
        "summary": str(summary_path),
        "normalization": config["data"].get("normalization", ""),
        "expression_representation": config["data"].get(
            "expression_representation", "full_transcriptome_tpm"
        ),
        "harmonization": _get(
            config, "data", "preprocessing", "harmonization", default="none"
        ),
        "transductive_preprocessing": (
            _get(
                config,
                "data",
                "preprocessing",
                "unseen_study_policy",
                default="global_train_fallback",
            )
            == "transductive_unlabeled"
        ),
        "status": "missing_evaluation",
    }
    if not summary_path.exists():
        return record
    summary = _load_json(summary_path)
    heldout = summary.get("heldout_fidelity_transformed", {})
    train = summary.get("paper_train_fidelity_transformed", {})
    heldout_selection = summary.get("heldout_model_selection", {})
    train_selection = summary.get("paper_train_model_selection", {})
    effect = summary.get("flt_gc_effect_recovery", {})
    accession = summary.get("accession_effect_validation", {})
    gene_accession = accession.get("gene", {})
    pathway_accession = accession.get("pathway", {})
    utility = summary.get("flt_gc_classifier_utility", {})
    record.update(
        {
            "status": str(summary.get("status", "complete")),
            "device": str(summary.get("device", "")),
            "validation_samples": int(_get(summary, "samples", default=0)),
            "heldout_corr": heldout.get("correlation_matrix_agreement", np.nan),
            "heldout_precision": heldout.get("precision", np.nan),
            "heldout_recall": heldout.get("recall", np.nan),
            "heldout_f1": heldout.get("f1", np.nan),
            "heldout_aa": heldout.get("adversarial_accuracy", np.nan),
            "heldout_fd_ratio": heldout.get(
                "frechet_ratio_to_real_split_p95", np.nan
            ),
            "heldout_absolute_gate": _get(
                heldout_selection,
                "absolute_paper_fidelity_gate",
                "passed",
                default=False,
            ),
            "heldout_diversity_gate": _get(
                heldout_selection, "diversity_gate", "passed", default=False
            ),
            "heldout_memorization_gate": _get(
                heldout_selection, "memorization_gate", "passed", default=False
            ),
            "train_corr": train.get("correlation_matrix_agreement", np.nan),
            "train_precision": train.get("precision", np.nan),
            "train_recall": train.get("recall", np.nan),
            "train_f1": train.get("f1", np.nan),
            "train_aa": train.get("adversarial_accuracy", np.nan),
            "train_fd_ratio": train.get(
                "frechet_ratio_to_real_split_p95", np.nan
            ),
            "train_absolute_gate": _get(
                train_selection,
                "absolute_paper_fidelity_gate",
                "passed",
                default=False,
            ),
            "delta_correlation": effect.get("delta_correlation", np.nan),
            "direction_agreement": effect.get("direction_agreement", np.nan),
            "conditional_effect_gate": _get(
                summary, "conditional_effect_gate", "passed", default=False
            ),
            "accession_meta_correlation": gene_accession.get(
                "meta_effect_correlation", np.nan
            ),
            "accession_meta_direction": gene_accession.get(
                "meta_direction_agreement", np.nan
            ),
            "accession_effect_gate": _get(
                summary, "accession_effect_gate", "passed", default=False
            ),
            "real_gene_fdr_lt_005": gene_accession.get(
                "real_random_effects_fdr_lt_005", 0
            ),
            "real_gene_loo_stable_fdr_lt_005": gene_accession.get(
                "real_loo_stable_fdr_lt_005", 0
            ),
            "synthetic_gene_fdr_lt_005": gene_accession.get(
                "synthetic_random_effects_fdr_lt_005", 0
            ),
            "real_pathway_fdr_lt_005": pathway_accession.get(
                "real_random_effects_fdr_lt_005", 0
            ),
            "real_pathway_loo_stable_fdr_lt_005": pathway_accession.get(
                "real_loo_stable_fdr_lt_005", 0
            ),
            "synthetic_pathway_fdr_lt_005": pathway_accession.get(
                "synthetic_random_effects_fdr_lt_005", 0
            ),
            "real_train_real_balanced_accuracy": _get(
                utility,
                "real_train_real_evaluation",
                "balanced_accuracy",
            ),
            "synthetic_train_real_balanced_accuracy": _get(
                utility,
                "synthetic_train_real_evaluation",
                "balanced_accuracy",
            ),
            "condition_consistency_balanced_accuracy": _get(
                summary, "condition_consistency", "balanced_accuracy"
            ),
        }
    )
    checks = _get(
        heldout_selection,
        "absolute_paper_fidelity_gate",
        "checks",
        default={},
    )
    for field, _, _, _ in FIDELITY_METRICS:
        metric_key = {
            "heldout_corr": "correlation_matrix_agreement",
            "heldout_precision": "precision",
            "heldout_recall": "recall",
            "heldout_f1": "f1",
            "heldout_aa": "adversarial_accuracy",
            "heldout_fd_ratio": "frechet_ratio_to_real_split_p95",
        }[field]
        record[f"{field}_passed"] = bool(checks.get(metric_key, False))
    return record


def _plot_fidelity(frame: pd.DataFrame, output: Path) -> None:
    complete = frame.loc[frame["status"].eq("complete")].copy()
    if complete.empty:
        return
    labels = complete["label"].astype(str).tolist()
    positions = np.arange(len(labels))
    figure, axes = plt.subplots(3, 2, figsize=(15, 14), constrained_layout=True)
    for axis, (field, title, threshold, policy) in zip(axes.flat, FIDELITY_METRICS):
        values = pd.to_numeric(complete[field], errors="coerce").to_numpy()
        passed = complete[f"{field}_passed"].astype(bool).to_numpy()
        colors = np.where(passed, "#1b7837", "#b2182b")
        axis.scatter(values, positions, c=colors, s=62, zorder=3)
        if policy == "minimum":
            axis.axvspan(float(threshold), 1.05, color="#d9f0d3", alpha=0.45)
            axis.axvline(float(threshold), color="#1b7837", linestyle="--")
            axis.set_xlim(0, 1.05)
        elif policy == "maximum":
            maximum = max(1.2, float(np.nanmax(values)) * 1.08)
            axis.axvspan(0, float(threshold), color="#d9f0d3", alpha=0.45)
            axis.axvline(float(threshold), color="#1b7837", linestyle="--")
            axis.set_xlim(0, maximum)
        else:
            lower, upper = threshold
            axis.axvspan(lower, upper, color="#d9f0d3", alpha=0.45)
            axis.axvline(lower, color="#1b7837", linestyle="--")
            axis.axvline(upper, color="#1b7837", linestyle="--")
            axis.set_xlim(0, 1)
        axis.set_title(title)
        axis.set_yticks(positions, labels)
        axis.invert_yaxis()
        axis.grid(axis="x", alpha=0.2)
    figure.suptitle("Held-out liver fidelity: independent acceptance metrics", fontsize=16)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def _plot_effects(frame: pd.DataFrame, output: Path) -> None:
    complete = frame.loc[frame["status"].eq("complete")].copy()
    if complete.empty:
        return
    labels = complete["label"].astype(str).tolist()
    positions = np.arange(len(labels))
    panels = (
        ("delta_correlation", "FLT-GC delta correlation", 0.30),
        ("direction_agreement", "FLT-GC direction agreement", 0.55),
        ("accession_meta_correlation", "Accession meta-effect correlation", 0.30),
        (
            "synthetic_train_real_balanced_accuracy",
            "Synthetic-train / real-validation BA",
            0.50,
        ),
    )
    figure, axes = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
    for axis, (field, title, threshold) in zip(axes.flat, panels):
        values = pd.to_numeric(complete[field], errors="coerce").to_numpy()
        axis.scatter(values, positions, color="#2166ac", s=62, zorder=3)
        axis.axvline(threshold, color="#1b7837", linestyle="--")
        lower = min(-0.1, float(np.nanmin(values)) - 0.05)
        axis.set_xlim(lower, 1.05)
        axis.set_title(title)
        axis.set_yticks(positions, labels)
        axis.invert_yaxis()
        axis.grid(axis="x", alpha=0.2)
    figure.suptitle("Held-out liver conditional-effect diagnostics", fontsize=16)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def summarize(manifest_path: str | Path) -> Path:
    manifest_file = Path(manifest_path)
    manifest = yaml.safe_load(manifest_file.read_text(encoding="utf-8"))
    split = str(manifest.get("evaluation_split", "validation"))
    records = [
        _record(str(item["label"]), Path(item["config"]), split)
        for item in manifest["runs"]
    ]
    frame = pd.DataFrame(records)
    output = Path(manifest["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    table_path = output / "independent_metrics.tsv"
    frame.to_csv(table_path, sep="\t", index=False, na_rep="NA")
    json_records = json.loads(frame.to_json(orient="records"))
    (output / "independent_metrics.json").write_text(
        json.dumps(json_records, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _plot_fidelity(frame, output / "heldout_fidelity_independent_metrics.png")
    _plot_effects(frame, output / "conditional_effect_diagnostics.png")
    completed = int(frame["status"].eq("complete").sum())
    readme = f"""# Liver harmonization benchmark

Matched conditional DDIM runs on the same OSDR API-derived liver cohort and
accession split. Completed evaluations: {completed}/{len(frame)}.

- Validation accessions: OSD-137, OSD-457, OSD-48.
- Locked test accession: OSD-379; it was not evaluated.
- Every fidelity metric is reported and gated independently. There is no
  composite score.
- ComBat, ComBat-seq, and MBatch held-out transforms are explicitly
  transductive sensitivity analyses.
- `independent_metrics.tsv` is the machine-readable comparison table.
- `heldout_fidelity_independent_metrics.png` shows the six independent paper
  fidelity criteria.
- `conditional_effect_diagnostics.png` shows FLT/GC and accession diagnostics.
"""
    (output / "README.md").write_text(readme, encoding="utf-8")
    return table_path
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="configs/rna_diffusion/liver_harmonization_benchmark.yaml",
    )
    args = parser.parse_args()
    print(summarize(args.manifest))


if __name__ == "__main__":
    main()
