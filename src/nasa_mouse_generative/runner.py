"""Configure, train, checkpoint, and validate one generative benchmark run."""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import re
import time

import yaml
import pandas as pd

from .adapters import create_adapter
from .config import BenchmarkConfig, load_config_with_overrides
from .metrics import evaluate_model
from .profiles import (
    epochs_for_stage,
    learning_rate_for_stage,
    load_model_parameters,
    resolve_preprocessing_profile,
)
from .training_data import prepare_training_data, save_prepared_osdr


def _smoke_config(config: BenchmarkConfig) -> BenchmarkConfig:
    parameters = dict(config.training.model_parameters)
    common = {"batch_size": 16, "epochs": 1, "reference_epochs": 1, "finetune_epochs": 1}
    if config.training.model == "vinas_wgan_gp":
        specific = {
            "hidden_dims": [32],
            "noise_dim": 16,
            "critic_steps": 1,
            "gradient_penalty": 1.0,
        }
    elif config.training.model == "lacan_diffusion":
        specific = {
            "hidden_dim": 32,
            "n_blocks": 1,
            "diffusion_timesteps": 20,
            "sample_steps": 5,
            "n_landmarks": 16,
            "landmark_strategy": "hvg",
            "reconstruction_samples": 64,
            "use_amp": False,
        }
    else:
        specific = {
            "d": 32,
            "latents_L": 16,
            "blocks_D": 1,
            "heads_h": 4,
            "min_context_genes": 16,
            "min_target_genes_per_block": 4,
            "max_tokens": 64,
            "fourier_num_frequencies": 8,
            "predictor_depth": 2,
            "predictor_expansion_factor": 2,
            "ema_warmup_steps": 0,
        }
    parameters.update(common)
    parameters.update(specific)
    return replace(
        config,
        data=replace(
            config.data,
            archs4_sample_limit=min(config.data.archs4_sample_limit or 64, 64),
            osdr_sample_limit=min(config.data.osdr_sample_limit or 96, 96),
        ),
        features=replace(
            config.features,
            max_genes=min(config.features.max_genes or 64, 64),
            hvg_genes=min(config.features.hvg_genes, 64),
        ),
        training=replace(config.training, model_parameters=parameters),
        validation=replace(config.validation, max_metric_samples=64),
        execution=replace(config.execution, checkpoint_every_epochs=1),
    )


def _run_identity(
    config: BenchmarkConfig, parameters: dict, tissue: str | None, run_name: str
) -> tuple[str, str]:
    payload = {
        "config": config.to_dict(),
        "resolved_model_parameters": parameters,
        "tissue_override": tissue or "",
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(serialized.encode()).hexdigest()
    if run_name:
        identifier = re.sub(r"[^A-Za-z0-9_.-]+", "_", run_name).strip("_")
        if not identifier:
            raise ValueError("--run-name contains no usable characters")
    else:
        tissue_label = tissue or config.training.tissue_mode
        identifier = f"{tissue_label}_{config.training.regime}_{digest[:12]}"
    return identifier, digest


def _stages(config: BenchmarkConfig, data) -> list[tuple[str, object]]:
    if config.training.regime == "osdr_only":
        return [("osdr", data.train)]
    if config.training.regime == "archs4_only":
        if data.reference is None:
            raise ValueError("archs4_only requires a prepared reference")
        return [("reference", data.reference)]
    if data.reference is None:
        raise ValueError("ARCHS4 pretraining requires a prepared reference")
    return [("reference", data.reference), ("osdr_finetune", data.train)]


def _write_readme(run_dir: Path, summary: dict) -> None:
    lines = [
        "# Generative Benchmark Run",
        "",
        f"- Model: `{summary['model']}`",
        f"- Task: `{summary['task']}`",
        f"- Regime: `{summary['regime']}`",
        f"- Tissue mode: `{summary['tissue_mode']}`",
        f"- Tissues: `{', '.join(summary['data']['tissues'])}`",
        f"- Device: `{summary['device']['device']}` {summary['device']['cuda_device_name']}",
        f"- Reference profiles: {summary['data']['reference_samples']}",
        f"- OSDR training profiles: {summary['data']['partition_samples']['train']}",
        f"- Genes: {summary['genes']}",
        "",
        "The validation split is accession-held-out. The locked test split is not",
        "evaluated automatically and requires an explicit `evaluate --unlock-test` call.",
    ]
    (run_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _claim_run_identity(
    run_dir: Path, *, identifier: str, digest: str, model: str
) -> None:
    """Validate every durable run marker before writing a new identity file."""

    recorded: dict[str, str] = {}
    identity_path = run_dir / "run_identity.json"
    if identity_path.exists():
        payload = json.loads(identity_path.read_text(encoding="utf-8"))
        recorded["run_identity.json"] = str(payload.get("run_sha256", ""))
    summary_path = run_dir / "run_summary.json"
    if summary_path.exists():
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        recorded["run_summary.json"] = str(payload.get("run_sha256", ""))
    resolved_path = run_dir / "resolved_config.yaml"
    if resolved_path.exists():
        payload = yaml.safe_load(resolved_path.read_text(encoding="utf-8")) or {}
        recorded["resolved_config.yaml"] = str(
            payload.get("run", {}).get("sha256", "")
        )

    missing = sorted(name for name, value in recorded.items() if not value)
    if missing:
        raise ValueError(f"Run markers lack a configuration digest: {missing}")
    existing_digests = set(recorded.values())
    if len(existing_digests) > 1:
        raise ValueError(
            f"Run directory {run_dir} has inconsistent configuration markers: "
            f"{recorded}"
        )
    if existing_digests and digest not in existing_digests:
        existing_digest = next(iter(existing_digests))
        raise ValueError(
            f"Run name {identifier!r} already belongs to configuration "
            f"{existing_digest[:12]}; choose another --run-name"
        )
    if not identity_path.exists():
        identity_path.write_text(
            json.dumps(
                {
                    "run_id": identifier,
                    "run_sha256": digest,
                    "model": model,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )


def train_one(
    config: BenchmarkConfig,
    *,
    tissue: str | None = None,
    run_name: str = "",
) -> Path:
    config.validate()
    parameters = load_model_parameters(config)
    identifier, digest = _run_identity(config, parameters, tissue, run_name)
    run_dir = Path(config.output_root) / "runs" / config.training.model / identifier
    run_dir.mkdir(parents=True, exist_ok=True)
    _claim_run_identity(
        run_dir,
        identifier=identifier,
        digest=digest,
        model=config.training.model,
    )
    resolved = config.to_dict()
    resolved["resolved_model_parameters"] = parameters
    resolved["run"] = {
        "id": identifier,
        "sha256": digest,
        "tissue_override": tissue or "",
    }
    (run_dir / "resolved_config.yaml").write_text(
        yaml.safe_dump(resolved, sort_keys=False), encoding="utf-8"
    )
    started = time.time()
    data = prepare_training_data(config, tissue=tissue)
    data.encoder.save(run_dir / "categorical_encoder.json")
    data.preprocessor.save(run_dir)
    save_prepared_osdr(data, run_dir)
    (run_dir / "genes.tsv").write_text(
        "gene_id\n" + "\n".join(data.genes) + "\n", encoding="utf-8"
    )
    adapter = create_adapter(config, data, parameters, run_dir)
    for stage, partition in _stages(config, data):
        adapter.fit_stage(
            partition,
            stage=stage,
            epochs=epochs_for_stage(parameters, stage),
            learning_rate=learning_rate_for_stage(parameters, stage),
        )
    model_path = adapter.save_final()
    validation_path = ""
    if config.execution.evaluate_after_training and len(data.partitions["validation"]):
        validation_path = str(
            evaluate_model(
                adapter,
                data.partitions,
                data.preprocessor,
                split="validation",
                output_dir=run_dir / "evaluation",
                seed=config.training.seed,
                max_samples=config.validation.max_metric_samples,
                save_generated_matrix=config.execution.save_generated_matrix,
            )
        )
    summary = {
        "run_id": identifier,
        "run_sha256": digest,
        "model": config.training.model,
        "model_profile": config.training.model_profile,
        "task": config.training.task,
        "regime": config.training.regime,
        "tissue_mode": config.training.tissue_mode,
        "parameters": parameters,
        "genes": len(data.genes),
        "covariates": list(data.covariates),
        "device": adapter.device_summary(),
        "data": data.metadata,
        "completed_epochs": adapter.state.completed_epochs,
        "training_seconds": float(time.time() - started),
        "outputs": {
            "run_dir": str(run_dir),
            "model": str(model_path),
            "validation": validation_path,
            "resolved_config": str(run_dir / "resolved_config.yaml"),
        },
    }
    summary_path = run_dir / "run_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_readme(run_dir, summary)
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/generative/default.yaml")
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Dotted YAML override; may be repeated.",
    )
    parser.add_argument("--tissue", default="")
    parser.add_argument(
        "--all-tissues",
        action="store_true",
        help="Run every confirmatory/exploratory standalone tissue sequentially.",
    )
    parser.add_argument("--run-name", default="")
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Bound samples/features and train one epoch per stage.",
    )
    return parser.parse_args()


def run(args: argparse.Namespace) -> Path:
    config = load_config_with_overrides(args.config, args.set)
    config = resolve_preprocessing_profile(config)
    if args.smoke:
        config = _smoke_config(config)
        config.validate()
    if args.all_tissues:
        if args.tissue:
            raise ValueError("Use either --tissue or --all-tissues, not both")
        config = replace(
            config, training=replace(config.training, tissue_mode="per_tissue")
        )
        inventory = pd.read_csv(
            Path(config.output_root)
            / "data_audit"
            / "osdr"
            / "osdr_tissue_inventory.tsv",
            sep="\t",
        )
        tissues = inventory.loc[
            inventory["training_tier"].ne("pooled_only"), "tissue_canonical"
        ].astype(str)
        outputs = []
        for tissue in tissues:
            name = f"{args.run_name}_{tissue}" if args.run_name else ""
            outputs.append(
                str(train_one(config, tissue=tissue, run_name=name))
            )
        summary_path = (
            Path(config.output_root)
            / "runs"
            / config.training.model
            / "per_tissue_batch_summary.json"
        )
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(
            json.dumps({"tissues": list(tissues), "run_summaries": outputs}, indent=2)
            + "\n",
            encoding="utf-8",
        )
        return summary_path
    return train_one(
        config,
        tissue=args.tissue or None,
        run_name=args.run_name,
    )


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
