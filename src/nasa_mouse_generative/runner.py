"""Configure, train, checkpoint, and validate one generative benchmark run."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import shutil
import time

import yaml
import pandas as pd
import torch

from .adapters import create_adapter
from .config import BenchmarkConfig, load_config_with_overrides
from .metrics import evaluate_model
from .models import MODEL_REGISTRY
from .profiles import (
    epochs_for_stage,
    learning_rate_for_stage,
    load_model_parameters,
    resolve_preprocessing_profile,
)
from .training_data import prepare_training_data, save_prepared_osdr


GIB = 1024**3


@lru_cache(maxsize=None)
def _pipeline_source_manifest(model: str) -> dict[str, object]:
    source_root = Path(__file__).resolve().parents[1]
    package_names = ["nasa_mouse_generative"]
    package_names.extend(
        {
            "vinas_wgan_gp": ["nasa_mouse_wgan"],
            "lacan_diffusion": ["nasa_mouse_diffusion"],
            "genejepa": [],
        }[model]
    )
    files = sorted(
        path
        for package in package_names
        for path in (source_root / package).rglob("*.py")
    )
    hashes: dict[str, str] = {}
    combined = hashlib.sha256()
    for path in files:
        relative = str(path.relative_to(source_root))
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        hashes[relative] = digest
        combined.update(relative.encode())
        combined.update(b"\0")
        combined.update(digest.encode())
        combined.update(b"\n")
    return {"sha256": combined.hexdigest(), "files": hashes}


def _directory_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _storage_snapshot(path: Path) -> dict[str, float]:
    probe = path
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    usage = shutil.disk_usage(probe)
    return {
        "total_gb": float(usage.total / GIB),
        "used_gb": float(usage.used / GIB),
        "free_gb": float(usage.free / GIB),
    }


def _enforce_storage(config: BenchmarkConfig, run_dir: Path, *, stage: str) -> dict:
    snapshot = _storage_snapshot(run_dir)
    size_gb = _directory_bytes(run_dir) / GIB
    if snapshot["free_gb"] < config.execution.min_free_space_gb:
        raise RuntimeError(
            f"Storage guard stopped {stage}: {snapshot['free_gb']:.2f} GiB free, "
            f"minimum is {config.execution.min_free_space_gb:.2f} GiB"
        )
    if size_gb > config.execution.max_run_output_gb:
        raise RuntimeError(
            f"Storage guard stopped {stage}: run uses {size_gb:.2f} GiB, "
            f"maximum is {config.execution.max_run_output_gb:.2f} GiB"
        )
    return {**snapshot, "run_size_gb": float(size_gb), "stage": stage}


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
            "samples_per_epoch": 64,
            "num_workers": 0,
        }
    parameters.update(common)
    parameters.update(specific)
    harmonization_parameters = dict(
        config.preprocessing.harmonization_parameters
    )
    if config.preprocessing.harmonization == "mober":
        harmonization_parameters.update(
            {
                "epochs": 1,
                "batch_size": 16,
                "encoding_dim": 8,
                "projection_batch_size": 64,
            }
        )
    elif config.preprocessing.harmonization == "combat":
        harmonization_parameters.setdefault("confounded_covariate_policy", "drop")
    elif config.preprocessing.harmonization == "combat_seq":
        harmonization_parameters["anchor_samples"] = min(
            int(harmonization_parameters.get("anchor_samples", 32)), 32
        )
        harmonization_parameters.setdefault("singleton_batch_policy", "pool")
        harmonization_parameters.setdefault("confounded_covariate_policy", "drop")
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
        preprocessing=replace(
            config.preprocessing,
            harmonization_parameters=harmonization_parameters,
        ),
        training=replace(
            config.training,
            model_profile="practical_screen",
            model_parameters=parameters,
        ),
        validation=replace(config.validation, max_metric_samples=64),
        execution=replace(config.execution, checkpoint_every_epochs=1),
    )


def _run_identity(
    config: BenchmarkConfig, parameters: dict, tissue: str | None, run_name: str
) -> tuple[str, str]:
    pipeline_source = _pipeline_source_manifest(config.training.model)
    payload = {
        "config": config.to_dict(),
        "resolved_model_parameters": parameters,
        "tissue_override": tissue or "",
        "pipeline_source_sha256": pipeline_source["sha256"],
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
        f"- Training profiles: {summary['data']['partition_samples']['train']}",
        f"- Genes: {summary['genes']}",
        "",
        f"The split unit is {summary['data']['split'].get('split_unit', 'OSDR accession')}.",
        "The locked test split is not",
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
    storage_before = _enforce_storage(config, run_dir, stage="before_data_preparation")
    _claim_run_identity(
        run_dir,
        identifier=identifier,
        digest=digest,
        model=config.training.model,
    )
    resolved = config.to_dict()
    resolved["resolved_model_parameters"] = parameters
    resolved["pipeline_source"] = _pipeline_source_manifest(
        config.training.model
    )
    resolved["run"] = {
        "id": identifier,
        "sha256": digest,
        "tissue_override": tissue or "",
    }
    (run_dir / "resolved_config.yaml").write_text(
        yaml.safe_dump(resolved, sort_keys=False), encoding="utf-8"
    )
    started = time.time()
    print(
        f"[{identifier}] preparing {config.training.regime} data",
        flush=True,
    )
    data = prepare_training_data(config, tissue=tissue)
    print(
        f"[{identifier}] prepared train={len(data.train)} "
        f"validation={len(data.partitions['validation'])} "
        f"test={len(data.partitions['test'])} genes={len(data.genes)}",
        flush=True,
    )
    data.encoder.save(run_dir / "categorical_encoder.json")
    data.preprocessor.save(run_dir)
    save_prepared_osdr(
        data,
        run_dir,
        include_matrix=config.execution.save_prepared_data,
    )
    (run_dir / "genes.tsv").write_text(
        "gene_id\n" + "\n".join(data.genes) + "\n", encoding="utf-8"
    )
    adapter = create_adapter(config, data, parameters, run_dir)
    if adapter.device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(adapter.device)
    stage_runtime: dict[str, float] = {}
    for stage, partition in _stages(config, data):
        print(
            f"[{identifier}] training stage={stage} profiles={len(partition)} "
            f"epochs={epochs_for_stage(parameters, stage)} device={adapter.device}",
            flush=True,
        )
        stage_started = time.time()
        adapter.fit_stage(
            partition,
            stage=stage,
            epochs=epochs_for_stage(parameters, stage),
            learning_rate=learning_rate_for_stage(parameters, stage),
        )
        stage_runtime[stage] = float(time.time() - stage_started)
        _enforce_storage(config, run_dir, stage=f"after_{stage}")
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
                samples_per_covariate_profile=(
                    config.generation.samples_per_covariate_profile
                ),
                synthetic_to_real_ratios=(
                    config.generation.synthetic_to_real_ratios
                ),
            )
        )
    peak_memory_gb = (
        float(torch.cuda.max_memory_allocated(adapter.device) / GIB)
        if adapter.device.type == "cuda"
        else 0.0
    )
    storage_after = _enforce_storage(config, run_dir, stage="after_evaluation")
    checkpoint_retained = bool(config.execution.retain_training_checkpoint)
    if not checkpoint_retained and adapter.checkpoint_dir.exists():
        shutil.rmtree(adapter.checkpoint_dir)
        storage_after = _enforce_storage(
            config, run_dir, stage="after_checkpoint_cleanup"
        )
    summary = {
        "run_id": identifier,
        "run_sha256": digest,
        "model": config.training.model,
        "model_provenance": asdict(MODEL_REGISTRY[config.training.model]),
        "pipeline_source": _pipeline_source_manifest(config.training.model),
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
        "stage_training_seconds": stage_runtime,
        "cuda_peak_memory_gb": peak_memory_gb,
        "storage": {
            "before": storage_before,
            "after": storage_after,
            "checkpoint_every_epochs": config.execution.checkpoint_every_epochs,
            "training_checkpoint_retained": checkpoint_retained,
        },
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
    repeat_count = int(config.training.repeats)
    repeat_configs = [
        replace(
            config,
            training=replace(
                config.training,
                seed=config.training.seed + index,
                repeats=1,
            ),
        )
        for index in range(repeat_count)
    ]
    if args.all_tissues:
        if args.tissue:
            raise ValueError("Use either --tissue or --all-tissues, not both")
        repeat_configs = [
            replace(
                item, training=replace(item.training, tissue_mode="per_tissue")
            )
            for item in repeat_configs
        ]
        inventory = pd.read_csv(
            Path(config.output_root)
            / "data_audit"
            / "osdr"
            / "osdr_tissue_inventory.tsv",
            sep="\t",
        )
        tissues = inventory.loc[
            inventory["training_tier"].ne("pooled_only"), "tissue_canonical"
        ].astype(str).drop_duplicates()
        outputs = []
        for tissue in tissues:
            for repeat_config in repeat_configs:
                name = f"{args.run_name}_{tissue}" if args.run_name else ""
                if repeat_count > 1 and name:
                    name = f"{name}_seed{repeat_config.training.seed}"
                outputs.append(
                    str(train_one(repeat_config, tissue=tissue, run_name=name))
                )
        summary_path = (
            Path(config.output_root)
            / "runs"
            / config.training.model
            / "per_tissue_batch_summary.json"
        )
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(
            json.dumps(
                {
                    "tissues": list(tissues),
                    "requested_repeats": repeat_count,
                    "seeds": [item.training.seed for item in repeat_configs],
                    "run_summaries": outputs,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return summary_path
    if repeat_count == 1:
        return train_one(
            repeat_configs[0],
            tissue=args.tissue or None,
            run_name=args.run_name,
        )
    outputs = []
    for repeat_config in repeat_configs:
        name = args.run_name
        if name:
            name = f"{name}_seed{repeat_config.training.seed}"
        outputs.append(
            str(
                train_one(
                    repeat_config,
                    tissue=args.tissue or None,
                    run_name=name,
                )
            )
        )
    summary_path = (
        Path(config.output_root)
        / "runs"
        / config.training.model
        / "repeat_batch_summary.json"
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(
            {
                "requested_repeats": repeat_count,
                "seeds": [item.training.seed for item in repeat_configs],
                "run_summaries": outputs,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return summary_path


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
