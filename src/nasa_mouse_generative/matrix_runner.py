"""Execute staged experiment-matrix rows with durable resumable status."""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import traceback
from typing import Any

import pandas as pd
import yaml

from .config import BenchmarkConfig, load_config
from .experiment_plan import expand_matrix
from .profiles import load_model_parameters, resolve_preprocessing_profile
from .runner import _run_identity, _smoke_config, train_one


SELECTION_PLACEHOLDERS = {
    "best_shared",
    "selected_from_phase_1",
    "best_validated",
    "best_from_phase_2",
}

CONDITIONING_PROFILES = {
    "condition_tissue": ("condition", "tissue"),
    "condition_tissue_sex": ("condition", "tissue", "sex"),
}


def _row_id(row: dict[str, Any]) -> str:
    payload = {
        key: value
        for key, value in row.items()
        if key not in {"status", "purpose", "native_expression_generator", "row_id"}
    }
    encoded = json.dumps(payload, sort_keys=True, default=str).encode()
    return hashlib.sha256(encoded).hexdigest()[:16]


def _expand_tissues(
    table: pd.DataFrame, *, inventory_path: Path, tissue_filter: str
) -> pd.DataFrame:
    if not table["tissue_mode"].astype(str).eq("per_tissue").any():
        result = table.copy()
        result["tissue"] = ""
        result["row_id"] = [
            _row_id(row) for row in result.to_dict(orient="records")
        ]
        return result
    inventory = pd.read_csv(inventory_path, sep="\t")
    tissues = (
        inventory.loc[
            inventory["training_tier"].astype(str).ne("pooled_only"),
            "tissue_canonical",
        ]
        .astype(str)
        .drop_duplicates()
        .tolist()
    )
    if tissue_filter:
        if tissue_filter not in tissues:
            raise ValueError(f"Tissue {tissue_filter!r} is not standalone-eligible")
        tissues = [tissue_filter]
    rows: list[dict[str, Any]] = []
    for row in table.to_dict(orient="records"):
        row_tissues = tissues if str(row.get("tissue_mode")) == "per_tissue" else [""]
        for tissue in row_tissues:
            expanded = dict(row)
            expanded["tissue"] = tissue
            expanded["row_id"] = _row_id(expanded)
            rows.append(expanded)
    return pd.DataFrame(rows)


def _unresolved_reason(row: dict[str, Any], base: BenchmarkConfig) -> str:
    for key in ("preprocessing_profile", "feature_space", "harmonization", "study_policy"):
        if str(row.get(key, "")) in SELECTION_PLACEHOLDERS:
            return f"awaiting_selection:{key}"
    scope = str(row.get("accession_scope", "all_eligible"))
    if scope == "single" and len(base.data.osdr_include_accessions) != 1:
        return "awaiting_accession_selection:single"
    if scope == "selected" and not base.data.osdr_include_accessions:
        return "awaiting_accession_selection:selected"
    return ""


def config_for_row(base: BenchmarkConfig, row: dict[str, Any]) -> BenchmarkConfig:
    """Resolve one executable matrix row into the validated benchmark schema."""

    reason = _unresolved_reason(row, base)
    if reason:
        raise ValueError(reason)
    preprocessing_name = str(row.get("preprocessing_profile", "custom"))
    if preprocessing_name == "genejepa_native":
        preprocessing_name = "model_native"
    feature_name = str(row.get("feature_space", base.features.space))
    hvg_genes = base.features.hvg_genes
    if feature_name == "hvg_2000":
        feature_name = "hvg"
        hvg_genes = 2000
    elif feature_name == "hvg_4096":
        feature_name = "hvg"
        hvg_genes = 4096
    elif feature_name == "all_shared_tokens":
        feature_name = "all_shared"

    study_policy = str(row.get("study_policy", base.training.study_policy))
    conditioning_profile = str(
        row.get("conditioning_profile", "all_configured")
    )
    if conditioning_profile == "all_configured":
        covariates = list(base.training.conditioning_covariates)
    elif conditioning_profile in CONDITIONING_PROFILES:
        covariates = list(CONDITIONING_PROFILES[conditioning_profile])
    else:
        raise ValueError(
            f"Unknown conditioning_profile {conditioning_profile!r}; choose from "
            f"{['all_configured', *sorted(CONDITIONING_PROFILES)]}"
        )
    if study_policy == "conditioned" and "study" not in covariates:
        covariates.append("study")
    if study_policy == "not_conditioned":
        covariates = [value for value in covariates if value != "study"]
    condition_on_flight = bool(row.get("condition_on_flight", True))
    task = str(row.get("task", base.training.task))
    model = str(row.get("model", base.training.model))
    model_profile = str(row.get("model_profile", base.training.model_profile))
    regime = str(row.get("training_regime", base.training.regime))
    harmonization = str(row.get("harmonization", base.preprocessing.harmonization))
    transductive = bool(
        base.validation.allow_transductive_preprocessing
        or harmonization in {"combat", "combat_seq"}
    )
    config = replace(
        base,
        preprocessing=replace(
            base.preprocessing,
            profile=preprocessing_name,
            harmonization=harmonization,
        ),
        data=replace(
            base.data,
            osdr_accession_scope=str(
                row.get("accession_scope", base.data.osdr_accession_scope)
            ),
        ),
        features=replace(
            base.features,
            space=feature_name,
            hvg_genes=hvg_genes,
        ),
        training=replace(
            base.training,
            model=model,
            model_profile=model_profile,
            task=task,
            regime=regime,
            tissue_mode=str(row.get("tissue_mode", base.training.tissue_mode)),
            condition_on_flight=condition_on_flight,
            study_policy=study_policy,
            conditioning_covariates=tuple(covariates),
            seed=int(row.get("seed", base.training.seed)),
            repeats=1,
        ),
        validation=replace(
            base.validation,
            allow_transductive_preprocessing=transductive,
        ),
    )
    config = resolve_preprocessing_profile(config)
    config.validate()
    return config


def _resolve_row_ids(table: pd.DataFrame, base: BenchmarkConfig) -> pd.DataFrame:
    result = table.copy()
    identifiers: list[str] = []
    for row in result.to_dict(orient="records"):
        if str(row.get("status", "")) == "capability_blocked" or _unresolved_reason(
            row, base
        ):
            identifiers.append(_row_id(row))
            continue
        config = config_for_row(base, row)
        parameters = load_model_parameters(config)
        _, digest = _run_identity(
            config,
            parameters,
            str(row.get("tissue", "")) or None,
            "",
        )
        backend_config = str(row.get("backend_config", ""))
        if backend_config:
            backend_path = Path(backend_config)
            if backend_path.exists():
                digest = hashlib.sha256(
                    (
                        digest
                        + ":"
                        + hashlib.sha256(backend_path.read_bytes()).hexdigest()
                    ).encode()
                ).hexdigest()
        identifiers.append(digest[:16])
    result["row_id"] = identifiers
    return result


def _write_status(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    table.to_csv(temporary, sep="\t", index=False)
    temporary.replace(path)


def _initial_status(plan: pd.DataFrame, existing: pd.DataFrame | None) -> pd.DataFrame:
    result = plan.copy()
    result["in_current_plan"] = True
    result["status"] = result.get("status", "planned")
    for column, default in {
        "started_utc": "",
        "finished_utc": "",
        "run_summary": "",
        "error": "",
    }.items():
        result[column] = default
    if existing is None or existing.empty or "row_id" not in existing:
        return result
    by_id = existing.set_index("row_id")
    for index, row in result.iterrows():
        row_id = row["row_id"]
        if row_id not in by_id.index:
            continue
        previous = by_id.loc[row_id]
        if isinstance(previous, pd.DataFrame):
            previous = previous.iloc[-1]
        for column in ("status", "started_utc", "finished_utc", "run_summary", "error"):
            result.at[index, column] = previous.get(column, result.at[index, column])
    retired = existing.loc[~existing["row_id"].isin(result["row_id"])].copy()
    if not retired.empty:
        never_started = retired["status"].astype(str).isin(
            {
                "planned",
                "ready",
                "awaiting_selection",
                "awaiting_accession_selection",
                "capability_blocked",
            }
        )
        has_execution_record = pd.Series(False, index=retired.index)
        for column in ("started_utc", "run_summary", "error"):
            if column in retired:
                has_execution_record |= retired[column].astype(str).str.len().gt(0)
        retired = retired.loc[~never_started | has_execution_record].copy()
        retired["in_current_plan"] = False
        result = pd.concat([result, retired], ignore_index=True, sort=False)
    return result


def _execution_mask(
    status: pd.DataFrame, *, phases: list[str], tissue_filter: str
) -> pd.Series:
    mask = status["in_current_plan"].fillna(False).astype(bool)
    if phases:
        mask &= status["phase"].astype(str).isin(phases)
    if tissue_filter:
        mask &= status["tissue_mode"].astype(str).ne("per_tissue") | status[
            "tissue"
        ].astype(str).eq(tissue_filter)
    return mask


def run(args: argparse.Namespace) -> Path:
    payload = yaml.safe_load(Path(args.matrix).read_text(encoding="utf-8")) or {}
    plan = expand_matrix(payload)
    base = load_config(args.config)
    plan = _expand_tissues(
        plan,
        inventory_path=(
            Path(base.output_root)
            / "data_audit/osdr/osdr_tissue_inventory.tsv"
        ),
        tissue_filter="",
    )
    plan = _resolve_row_ids(plan, base)
    status_path = Path(args.status)
    existing = (
        pd.read_csv(status_path, sep="\t", keep_default_na=False)
        if status_path.exists()
        else None
    )
    status = _initial_status(plan, existing)
    selected = _execution_mask(
        status, phases=list(args.phase), tissue_filter=args.tissue
    )
    executed = 0
    for index, row_series in status.loc[selected].iterrows():
        row = row_series.to_dict()
        current = str(row.get("status", "planned"))
        if current == "complete":
            continue
        if current == "capability_blocked":
            continue
        if current == "failed" and not args.retry_failed:
            continue
        reason = _unresolved_reason(row, base)
        if reason:
            status.at[index, "status"] = reason.split(":", 1)[0]
            status.at[index, "error"] = reason
            continue
        if args.max_runs > 0 and executed >= args.max_runs:
            break
        try:
            config = config_for_row(base, row)
            if args.smoke:
                config = _smoke_config(config)
            if args.dry_run:
                status.at[index, "status"] = "ready"
                continue
            status.at[index, "status"] = "running"
            status.at[index, "started_utc"] = pd.Timestamp.utcnow().isoformat()
            status.at[index, "error"] = ""
            _write_status(status, status_path)
            name = f"matrix_{row['phase']}_{row['row_id']}"
            if str(row.get("execution_backend", "")) == "nasa_mouse_rna_diffusion":
                baseline_config = str(
                    row.get(
                        "backend_config",
                        "configs/rna_diffusion/archs4_mouse_paper_parity.yaml",
                    )
                )
                baseline_summary = Path(
                    "outputs/generative_benchmark/runs/lacan_diffusion/"
                    "archs4_mouse_paper_parity_seed1234/run_summary.json"
                )
                if baseline_summary.exists() and json.loads(
                    baseline_summary.read_text(encoding="utf-8")
                ).get("status") == "complete":
                    summary = baseline_summary
                else:
                    from nasa_mouse_rna_diffusion.train import train

                    train(baseline_config)
                    summary = baseline_summary
            else:
                summary = train_one(
                    config,
                    tissue=str(row.get("tissue", "")) or None,
                    run_name=name,
                )
            status.at[index, "status"] = "complete"
            status.at[index, "run_summary"] = str(summary)
            executed += 1
        except KeyboardInterrupt:
            status.at[index, "status"] = "interrupted"
            status.at[index, "error"] = "KeyboardInterrupt"
            raise
        except Exception as error:
            status.at[index, "status"] = "failed"
            status.at[index, "error"] = "".join(
                traceback.format_exception_only(type(error), error)
            ).strip()
            executed += 1
            if args.stop_on_failure:
                raise
        finally:
            status.at[index, "finished_utc"] = pd.Timestamp.utcnow().isoformat()
            _write_status(status, status_path)
    _write_status(status, status_path)
    print(status["status"].value_counts(dropna=False).to_string())
    return status_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", default="configs/generative/experiment_matrix.yaml")
    parser.add_argument("--config", default="configs/generative/default.yaml")
    parser.add_argument(
        "--status",
        default="outputs/generative_benchmark/summary/experiment_status.tsv",
    )
    parser.add_argument("--phase", action="append", default=[])
    parser.add_argument("--tissue", default="")
    parser.add_argument(
        "--max-runs",
        type=int,
        default=4,
        help="Maximum rows per invocation; pass 0 only for an intentional unlimited run.",
    )
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--stop-on-failure", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
