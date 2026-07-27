"""Serializable fold-aware harmonization adapters for bulk expression."""

from __future__ import annotations

from abc import ABC, abstractmethod
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Iterable

import numpy as np
import pandas as pd


EPSILON = 1e-8


def _metadata_frame(
    studies: Iterable[object], metadata: pd.DataFrame | None, samples: int
) -> pd.DataFrame:
    studies_array = np.asarray(list(studies), dtype=str)
    if len(studies_array) != samples:
        raise ValueError("studies must contain one value per sample")
    if metadata is None:
        frame = pd.DataFrame(index=np.arange(samples))
    else:
        frame = metadata.reset_index(drop=True).copy()
        if len(frame) != samples:
            raise ValueError("harmonization metadata must align to expression rows")
    frame["study"] = studies_array
    for column in frame.columns:
        frame[column] = frame[column].fillna("__missing__").astype(str)
    return frame


def _resolve_batch_key(frame: pd.DataFrame, parameters: dict[str, Any]) -> str:
    requested = str(parameters.get("batch_key", "auto"))
    if requested != "auto":
        if requested not in frame:
            raise ValueError(f"Harmonization batch key is absent: {requested}")
        return requested
    if "source" in frame and frame["source"].nunique() > 1:
        return "source"
    return "study"


def _batch_values(
    frame: pd.DataFrame, batch_key: str, *, known_key: str | None = None
) -> np.ndarray:
    key = known_key or batch_key
    if key not in frame:
        if key == "study":
            return frame["study"].astype(str).to_numpy()
        raise ValueError(f"Harmonization metadata lacks batch column {key!r}")
    return frame[key].astype(str).to_numpy()


def _apply_covariate_schema(
    frame: pd.DataFrame, schema: list[dict[str, object]]
) -> np.ndarray:
    columns: list[np.ndarray] = []
    for item in schema:
        name = str(item["name"])
        values = (
            frame[name].fillna("__missing__").astype(str)
            if name in frame
            else pd.Series("__missing__", index=frame.index)
        )
        for level in list(item["levels"])[1:]:
            columns.append(values.eq(str(level)).to_numpy(dtype=np.float64))
    return (
        np.column_stack(columns)
        if columns
        else np.empty((len(frame), 0), dtype=np.float64)
    )


def _rank_safe_covariates(
    frame: pd.DataFrame,
    covariates: tuple[str, ...],
    batch_design: np.ndarray,
    *,
    policy: str,
) -> tuple[list[dict[str, object]], np.ndarray, list[str]]:
    design = np.asarray(batch_design, dtype=np.float64)
    schema: list[dict[str, object]] = []
    retained_columns: list[np.ndarray] = []
    dropped: list[str] = []
    for covariate in covariates:
        if covariate not in frame or frame[covariate].nunique() < 2:
            continue
        values = frame[covariate].fillna("__missing__").astype(str)
        levels = sorted(values.unique())
        encoded = np.column_stack(
            [values.eq(level).to_numpy(dtype=np.float64) for level in levels[1:]]
        )
        previous_rank = np.linalg.matrix_rank(design)
        candidate = np.concatenate([design, encoded], axis=1)
        if np.linalg.matrix_rank(candidate) - previous_rank != encoded.shape[1]:
            dropped.append(covariate)
            continue
        design = candidate
        schema.append({"name": covariate, "levels": levels})
        retained_columns.extend(encoded[:, index] for index in range(encoded.shape[1]))
    if dropped and policy != "drop":
        raise ValueError(
            "Harmonization covariates are confounded with batch: "
            f"{dropped}. Set confounded_covariate_policy=drop only as an "
            "explicitly audited sensitivity analysis."
        )
    matrix = (
        np.column_stack(retained_columns)
        if retained_columns
        else np.empty((len(frame), 0), dtype=np.float64)
    )
    return schema, matrix, dropped


class Harmonizer(ABC):
    method_id: str
    requires_raw_counts = False
    is_transductive = False

    def __init__(
        self,
        *,
        covariates: tuple[str, ...],
        parameters: dict[str, Any],
        device_spec: str,
        seed: int,
    ) -> None:
        self.covariates = tuple(covariates)
        self.parameters = dict(parameters)
        self.device_spec = str(device_spec)
        self.seed = int(seed)

    @abstractmethod
    def fit_transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
    ) -> np.ndarray:
        raise NotImplementedError

    @abstractmethod
    def transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
        *,
        allow_transductive: bool,
    ) -> np.ndarray:
        raise NotImplementedError

    @abstractmethod
    def save(self, directory: Path) -> dict[str, object]:
        raise NotImplementedError

    @abstractmethod
    def audit(self) -> dict[str, object]:
        raise NotImplementedError


def _combat_posterior(
    standardized: np.ndarray,
    gamma_hat: np.ndarray,
    delta_hat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    gamma_bar = float(np.mean(gamma_hat))
    t2 = max(float(np.var(gamma_hat)), EPSILON)
    delta = np.nan_to_num(delta_hat, nan=1.0, posinf=1.0, neginf=1.0)
    delta = np.maximum(delta, EPSILON)
    mean_delta = max(float(np.mean(delta)), EPSILON)
    variance_delta = max(float(np.var(delta)), EPSILON)
    a_prior = (2.0 * variance_delta + mean_delta**2) / variance_delta
    b_prior = (mean_delta * variance_delta + mean_delta**3) / variance_delta
    sample_count = standardized.shape[0]
    gamma = gamma_hat.astype(np.float64, copy=True)
    updated_delta = delta.astype(np.float64, copy=True)
    for _ in range(100):
        next_gamma = (
            t2 * sample_count * gamma_hat + updated_delta * gamma_bar
        ) / (t2 * sample_count + updated_delta)
        residual = standardized - next_gamma.reshape(1, -1)
        next_delta = (0.5 * np.square(residual).sum(axis=0) + b_prior) / (
            sample_count / 2.0 + a_prior - 1.0
        )
        gamma_change = np.max(
            np.abs(next_gamma - gamma) / np.maximum(np.abs(gamma), EPSILON)
        )
        delta_change = np.max(
            np.abs(next_delta - updated_delta)
            / np.maximum(np.abs(updated_delta), EPSILON)
        )
        gamma, updated_delta = next_gamma, np.maximum(next_delta, EPSILON)
        if max(float(gamma_change), float(delta_change)) < 1e-4:
            break
    return gamma.astype(np.float32), updated_delta.astype(np.float32)


class CombatHarmonizer(Harmonizer):
    method_id = "combat"
    is_transductive = True

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.batch_key = ""
        self.batch_levels: list[str] = []
        self.covariate_design: list[dict[str, object]] = []
        self.grand_mean: np.ndarray | None = None
        self.pooled_variance: np.ndarray | None = None
        self.covariate_coefficients: np.ndarray | None = None
        self.gamma: np.ndarray | None = None
        self.delta: np.ndarray | None = None
        self.dropped_covariates: list[str] = []

    def fit_transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
    ) -> np.ndarray:
        values = np.asarray(matrix, dtype=np.float64)
        frame = _metadata_frame(studies, metadata, len(values))
        self.batch_key = _resolve_batch_key(frame, self.parameters)
        batches = _batch_values(frame, self.batch_key)
        self.batch_levels = sorted(set(batches))
        max_batches = int(self.parameters.get("max_batches", 512))
        if len(self.batch_levels) > max_batches:
            raise ValueError(
                f"ComBat resolved {len(self.batch_levels)} batches, above max_batches="
                f"{max_batches}. Use batch_key=source, restrict the cohort, or raise "
                "the limit deliberately."
            )
        batch_matrix = np.column_stack(
            [batches == level for level in self.batch_levels]
        ).astype(np.float64)
        self.covariate_design, covariates, self.dropped_covariates = (
            _rank_safe_covariates(
                frame,
                self.covariates,
                batch_matrix,
                policy=str(
                    self.parameters.get("confounded_covariate_policy", "error")
                ),
            )
        )
        design = np.concatenate([batch_matrix, covariates], axis=1)
        coefficients = np.linalg.pinv(design) @ values
        batch_coefficients = coefficients[: len(self.batch_levels)]
        batch_counts = batch_matrix.sum(axis=0)
        self.grand_mean = (
            batch_counts / max(float(len(values)), 1.0)
        ) @ batch_coefficients
        self.covariate_coefficients = coefficients[len(self.batch_levels) :].astype(
            np.float32
        )
        stand_mean = self._stand_mean(covariates)
        residual = values - design @ coefficients
        pooled = np.mean(np.square(residual), axis=0)
        self.pooled_variance = np.maximum(pooled, EPSILON).astype(np.float32)
        standardized = (values - stand_mean) / np.sqrt(self.pooled_variance)
        gamma_rows: list[np.ndarray] = []
        delta_rows: list[np.ndarray] = []
        for level in self.batch_levels:
            block = standardized[batches == level]
            gamma_hat = block.mean(axis=0)
            delta_hat = (
                block.var(axis=0, ddof=1)
                if len(block) > 1
                else np.ones(values.shape[1], dtype=np.float64)
            )
            gamma, delta = _combat_posterior(block, gamma_hat, delta_hat)
            gamma_rows.append(gamma)
            delta_rows.append(delta)
        self.gamma = np.stack(gamma_rows)
        self.delta = np.stack(delta_rows)
        return self._adjust(values, batches, covariates, allow_transductive=False)

    def _stand_mean(self, covariates: np.ndarray) -> np.ndarray:
        if self.grand_mean is None or self.covariate_coefficients is None:
            raise RuntimeError("ComBat is not fitted")
        result = np.broadcast_to(
            self.grand_mean, (len(covariates), len(self.grand_mean))
        ).copy()
        if covariates.shape[1]:
            result += covariates @ self.covariate_coefficients
        return result

    def _adjust(
        self,
        values: np.ndarray,
        batches: np.ndarray,
        covariates: np.ndarray,
        *,
        allow_transductive: bool,
    ) -> np.ndarray:
        if self.pooled_variance is None or self.gamma is None or self.delta is None:
            raise RuntimeError("ComBat is not fitted")
        stand_mean = self._stand_mean(covariates)
        standardized = (values - stand_mean) / np.sqrt(self.pooled_variance)
        adjusted = np.empty_like(standardized)
        level_to_index = {level: index for index, level in enumerate(self.batch_levels)}
        for level in sorted(set(batches)):
            mask = batches == level
            if level in level_to_index:
                index = level_to_index[level]
                gamma, delta = self.gamma[index], self.delta[index]
            else:
                if not allow_transductive:
                    raise ValueError(
                        f"ComBat encountered unseen batch {level!r}; enable explicit "
                        "transductive preprocessing to estimate its adjustment."
                    )
                block = standardized[mask]
                gamma_hat = block.mean(axis=0)
                delta_hat = (
                    block.var(axis=0, ddof=1)
                    if len(block) > 1
                    else np.ones(values.shape[1], dtype=np.float64)
                )
                gamma, delta = _combat_posterior(block, gamma_hat, delta_hat)
            adjusted[mask] = (standardized[mask] - gamma) / np.sqrt(delta)
        result = adjusted * np.sqrt(self.pooled_variance) + stand_mean
        return np.nan_to_num(result, copy=False).astype(np.float32)

    def transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
        *,
        allow_transductive: bool,
    ) -> np.ndarray:
        values = np.asarray(matrix, dtype=np.float64)
        frame = _metadata_frame(studies, metadata, len(values))
        batches = _batch_values(frame, self.batch_key)
        covariates = _apply_covariate_schema(frame, self.covariate_design)
        return self._adjust(
            values, batches, covariates, allow_transductive=allow_transductive
        )

    def save(self, directory: Path) -> dict[str, object]:
        if any(
            value is None
            for value in (
                self.grand_mean,
                self.pooled_variance,
                self.covariate_coefficients,
                self.gamma,
                self.delta,
            )
        ):
            raise RuntimeError("Cannot save an unfitted ComBat harmonizer")
        path = directory / "combat_harmonizer.npz"
        np.savez_compressed(
            path,
            grand_mean=self.grand_mean,
            pooled_variance=self.pooled_variance,
            covariate_coefficients=self.covariate_coefficients,
            gamma=self.gamma,
            delta=self.delta,
        )
        return {**self.audit(), "artifact": path.name}

    @classmethod
    def load(cls, directory: Path, payload: dict[str, object]) -> "CombatHarmonizer":
        harmonizer = cls(
            covariates=tuple(payload.get("covariates", ())),
            parameters=dict(payload.get("parameters", {})),
            device_spec="cpu",
            seed=int(payload.get("seed", 0)),
        )
        harmonizer.batch_key = str(payload["batch_key"])
        harmonizer.batch_levels = list(map(str, payload["batch_levels"]))
        harmonizer.covariate_design = list(payload.get("covariate_design", []))
        harmonizer.dropped_covariates = list(
            map(str, payload.get("dropped_covariates", []))
        )
        arrays = np.load(directory / str(payload["artifact"]))
        harmonizer.grand_mean = np.asarray(arrays["grand_mean"], dtype=np.float32)
        harmonizer.pooled_variance = np.asarray(
            arrays["pooled_variance"], dtype=np.float32
        )
        harmonizer.covariate_coefficients = np.asarray(
            arrays["covariate_coefficients"], dtype=np.float32
        )
        harmonizer.gamma = np.asarray(arrays["gamma"], dtype=np.float32)
        harmonizer.delta = np.asarray(arrays["delta"], dtype=np.float32)
        return harmonizer

    def audit(self) -> dict[str, object]:
        return {
            "method": self.method_id,
            "fold_behavior": "frozen_known_batches_transductive_unseen_batches",
            "batch_key": self.batch_key,
            "batch_levels": self.batch_levels,
            "covariates": list(self.covariates),
            "covariate_design": self.covariate_design,
            "retained_covariates": [
                str(item["name"]) for item in self.covariate_design
            ],
            "dropped_covariates": self.dropped_covariates,
            "outcome_informed": any(
                str(item["name"]) == "condition"
                for item in self.covariate_design
            ),
            "outcome_covariate_requested": "condition" in self.covariates,
            "parameters": self.parameters,
            "seed": self.seed,
        }


class CombatSeqHarmonizer(Harmonizer):
    method_id = "combat_seq"
    requires_raw_counts = True
    is_transductive = True

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.batch_key = ""
        self.batch_levels: list[str] = []
        self.anchor_matrix: np.ndarray | None = None
        self.anchor_metadata: pd.DataFrame | None = None
        self.rounding_audit: list[dict[str, object]] = []
        self.singleton_batch_audit: list[dict[str, object]] = []
        self.confounded_covariate_audit: list[dict[str, object]] = []
        self.r_runtime: dict[str, str] = {}

    def _integer_counts(self, matrix: np.ndarray) -> np.ndarray:
        values = np.asarray(matrix, dtype=np.float64)
        if (values < 0).any() or not np.isfinite(values).all():
            raise ValueError("ComBat-seq requires finite non-negative raw counts")
        rounded = np.rint(values)
        noninteger = np.abs(values - rounded) > 1e-6
        if noninteger.any():
            policy = str(self.parameters.get("noninteger_policy", "error"))
            if policy != "round":
                raise ValueError(
                    "ComBat-seq requires integer raw counts, but the selected API "
                    "matrix contains fractional count-like values. Set "
                    "preprocessing.harmonization_parameters.noninteger_policy=round "
                    "only as an explicitly audited sensitivity analysis."
                )
            self.rounding_audit.append(
                {
                    "samples": int(values.shape[0]),
                    "genes": int(values.shape[1]),
                    "entries_rounded": int(noninteger.sum()),
                    "fraction_rounded": float(noninteger.mean()),
                    "maximum_fractional_distance": float(
                        np.abs(values - rounded).max()
                    ),
                }
            )
        return rounded.astype(np.int64)

    def _run(self, matrix: np.ndarray, frame: pd.DataFrame) -> np.ndarray:
        rscript = str(self.parameters.get("rscript", "Rscript"))
        executable = shutil.which(rscript) if "/" not in rscript else rscript
        if not executable and rscript == "Rscript":
            environment_rscript = Path(sys.executable).with_name("Rscript")
            executable = (
                str(environment_rscript) if environment_rscript.exists() else None
            )
        if not executable or not Path(executable).exists():
            raise RuntimeError(
                "ComBat-seq requires Rscript and the Bioconductor sva package. "
                "Install the generative environment requirements first."
            )
        batches = _batch_values(frame, self.batch_key)
        batch_counts = pd.Series(batches).value_counts()
        singleton_levels = set(
            batch_counts.loc[batch_counts.eq(1)].index.astype(str)
        )
        passthrough: np.ndarray | None = None
        passthrough_mask: np.ndarray | None = None
        if singleton_levels:
            policy = str(self.parameters.get("singleton_batch_policy", "error"))
            if policy not in {"error", "identity", "pool"}:
                raise ValueError(
                    "singleton_batch_policy must be error, identity, or pool"
                )
            if policy == "error":
                raise ValueError(
                    "ComBat-seq does not support one-sample batches. Set "
                    "singleton_batch_policy=identity to leave them uncorrected or "
                    "singleton_batch_policy=pool for an exploratory pooled arm."
                )
            singleton_mask = np.isin(batches, list(singleton_levels))
            effective_policy = policy
            if policy == "pool" and int(singleton_mask.sum()) > 1:
                batches = batches.copy()
                batches[singleton_mask] = "__pooled_singleton_batches__"
            else:
                effective_policy = "identity"
                passthrough = matrix.astype(np.float32, copy=True)
                passthrough_mask = ~singleton_mask
                matrix = matrix[passthrough_mask]
                frame = frame.loc[passthrough_mask].reset_index(drop=True)
                batches = batches[passthrough_mask]
            self.singleton_batch_audit.append(
                {
                    "policy": effective_policy,
                    "singleton_batches": len(singleton_levels),
                    "singleton_samples": int(singleton_mask.sum()),
                }
            )
        if len(set(batches)) < 2:
            return (
                passthrough
                if passthrough is not None
                else matrix.astype(np.float32)
            )
        metadata = pd.DataFrame({"batch": batches})
        batch_design = pd.get_dummies(
            pd.Series(batches, dtype="object"), drop_first=False, dtype=float
        ).to_numpy()
        design = batch_design
        retained: list[str] = []
        confounded: list[str] = []
        covariate_columns: list[tuple[str, np.ndarray]] = []
        for covariate in self.covariates:
            if covariate not in frame or frame[covariate].nunique() < 2:
                continue
            encoded_frame = pd.get_dummies(
                frame[covariate].astype(str), drop_first=True, dtype=float
            )
            encoded = encoded_frame.to_numpy()
            previous_rank = np.linalg.matrix_rank(design)
            candidate = np.concatenate([design, encoded], axis=1)
            if np.linalg.matrix_rank(candidate) - previous_rank != encoded.shape[1]:
                confounded.append(covariate)
                continue
            design = candidate
            retained.append(covariate)
            if covariate == "condition":
                metadata["group"] = frame[covariate].astype(str).to_numpy()
            else:
                for level_index in range(encoded.shape[1]):
                    covariate_columns.append(
                        (f"cov_{len(covariate_columns):03d}", encoded[:, level_index])
                    )
        if confounded:
            policy = str(
                self.parameters.get("confounded_covariate_policy", "error")
            )
            if policy != "drop":
                raise ValueError(
                    "ComBat-seq preservation covariates are confounded with batch: "
                    f"{confounded}. Set confounded_covariate_policy=drop only as an "
                    "explicitly audited sensitivity analysis."
                )
        for name, values in covariate_columns:
            metadata[name] = values
        self.confounded_covariate_audit.append(
            {"retained": retained, "dropped_as_confounded": confounded}
        )
        script = r'''args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(sva))
counts_samples <- as.matrix(read.table(args[1], header=FALSE, sep="\t", check.names=FALSE))
meta <- read.delim(args[2], stringsAsFactors=FALSE, check.names=FALSE)
batch <- factor(meta$batch)
group <- if ("group" %in% names(meta)) factor(meta$group) else NULL
cov_names <- grep("^cov_", names(meta), value=TRUE)
covar_mod <- NULL
if (length(cov_names) > 0) covar_mod <- as.matrix(meta[, cov_names, drop=FALSE])
adjusted <- ComBat_seq(
  counts=t(counts_samples), batch=batch, group=group, covar_mod=covar_mod,
  full_mod=TRUE, shrink=FALSE, shrink.disp=FALSE
)
write.table(t(adjusted), file=args[3], sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
writeLines(c(R.version.string, as.character(packageVersion("sva"))), con=args[4])
'''
        with tempfile.TemporaryDirectory(prefix="nasa_mouse_combat_seq_") as tmp:
            root = Path(tmp)
            counts_path = root / "counts.tsv"
            metadata_path = root / "metadata.tsv"
            output_path = root / "adjusted.tsv"
            runtime_path = root / "runtime.txt"
            script_path = root / "combat_seq.R"
            np.savetxt(counts_path, matrix, delimiter="\t", fmt="%d")
            metadata.to_csv(metadata_path, sep="\t", index=False)
            script_path.write_text(script, encoding="utf-8")
            completed = subprocess.run(
                [
                    str(executable),
                    str(script_path),
                    str(counts_path),
                    str(metadata_path),
                    str(output_path),
                    str(runtime_path),
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=int(self.parameters.get("timeout_seconds", 7200)),
            )
            if completed.returncode != 0:
                raise RuntimeError(
                    "sva::ComBat_seq failed:\n"
                    + (completed.stderr or completed.stdout)[-4000:]
                )
            adjusted = np.loadtxt(output_path, delimiter="\t", ndmin=2)
            runtime_lines = runtime_path.read_text(encoding="utf-8").splitlines()
            self.r_runtime = {
                "r": runtime_lines[0] if runtime_lines else "unknown",
                "sva": runtime_lines[1] if len(runtime_lines) > 1 else "unknown",
            }
        if adjusted.shape != matrix.shape:
            raise RuntimeError(
                f"ComBat-seq returned {adjusted.shape}, expected {matrix.shape}"
            )
        if passthrough is not None and passthrough_mask is not None:
            passthrough[passthrough_mask] = adjusted
            return passthrough
        return adjusted.astype(np.float32)

    def _select_anchors(
        self, matrix: np.ndarray, frame: pd.DataFrame
    ) -> tuple[np.ndarray, pd.DataFrame]:
        limit = int(self.parameters.get("anchor_samples", 512))
        if limit <= 0 or len(frame) <= limit:
            return matrix.copy(), frame.copy()
        rng = np.random.default_rng(self.seed)
        batches = _batch_values(frame, self.batch_key)
        selected: list[int] = []
        groups = {
            level: np.flatnonzero(batches == level)
            for level in sorted(set(batches))
        }
        while len(selected) < limit and any(len(indices) for indices in groups.values()):
            for level in list(groups):
                indices = groups[level]
                if not len(indices) or len(selected) >= limit:
                    continue
                choice = int(rng.choice(indices))
                selected.append(choice)
                groups[level] = indices[indices != choice]
        selected = sorted(selected)
        return matrix[selected].copy(), frame.iloc[selected].reset_index(drop=True)

    def fit_transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
    ) -> np.ndarray:
        counts = self._integer_counts(matrix)
        frame = _metadata_frame(studies, metadata, len(counts))
        self.batch_key = _resolve_batch_key(frame, self.parameters)
        self.batch_levels = sorted(set(_batch_values(frame, self.batch_key)))
        self.anchor_matrix, self.anchor_metadata = self._select_anchors(counts, frame)
        return self._run(counts, frame)

    def transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
        *,
        allow_transductive: bool,
    ) -> np.ndarray:
        if not allow_transductive:
            raise ValueError(
                "ComBat-seq transformation is transductive and requires explicit permission"
            )
        if self.anchor_matrix is None or self.anchor_metadata is None:
            raise RuntimeError("ComBat-seq is not fitted")
        counts = self._integer_counts(matrix)
        frame = _metadata_frame(studies, metadata, len(counts))
        combined = np.concatenate([self.anchor_matrix, counts])
        combined_frame = pd.concat(
            [self.anchor_metadata, frame], ignore_index=True, sort=False
        ).fillna("__missing__")
        adjusted = self._run(combined, combined_frame)
        return adjusted[len(self.anchor_matrix) :]

    def save(self, directory: Path) -> dict[str, object]:
        if self.anchor_matrix is None or self.anchor_metadata is None:
            raise RuntimeError("Cannot save an unfitted ComBat-seq harmonizer")
        arrays_path = directory / "combat_seq_harmonizer.npz"
        metadata_path = directory / "combat_seq_anchor_metadata.tsv.gz"
        np.savez_compressed(arrays_path, anchor_matrix=self.anchor_matrix)
        self.anchor_metadata.to_csv(
            metadata_path, sep="\t", index=False, compression="gzip"
        )
        return {
            **self.audit(),
            "artifact": arrays_path.name,
            "anchor_metadata": metadata_path.name,
        }

    @classmethod
    def load(cls, directory: Path, payload: dict[str, object]) -> "CombatSeqHarmonizer":
        harmonizer = cls(
            covariates=tuple(payload.get("covariates", ())),
            parameters=dict(payload.get("parameters", {})),
            device_spec="cpu",
            seed=int(payload.get("seed", 0)),
        )
        harmonizer.batch_key = str(payload["batch_key"])
        harmonizer.batch_levels = list(map(str, payload["batch_levels"]))
        harmonizer.rounding_audit = list(payload.get("rounding_audit", []))
        harmonizer.singleton_batch_audit = list(
            payload.get("singleton_batch_audit", [])
        )
        harmonizer.confounded_covariate_audit = list(
            payload.get("confounded_covariate_audit", [])
        )
        harmonizer.r_runtime = dict(payload.get("r_runtime", {}))
        arrays = np.load(directory / str(payload["artifact"]))
        harmonizer.anchor_matrix = np.asarray(arrays["anchor_matrix"], dtype=np.int64)
        harmonizer.anchor_metadata = pd.read_csv(
            directory / str(payload["anchor_metadata"]), sep="\t"
        ).fillna("__missing__")
        return harmonizer

    def audit(self) -> dict[str, object]:
        return {
            "method": self.method_id,
            "fold_behavior": "transductive_with_training_anchors",
            "batch_key": self.batch_key,
            "batch_levels": self.batch_levels,
            "covariates": list(self.covariates),
            "outcome_informed": any(
                "condition" in item.get("retained", [])
                for item in self.confounded_covariate_audit
            ),
            "outcome_covariate_requested": "condition" in self.covariates,
            "parameters": self.parameters,
            "seed": self.seed,
            "anchor_samples": (
                0 if self.anchor_matrix is None else len(self.anchor_matrix)
            ),
            "rounding_audit": self.rounding_audit,
            "singleton_batch_audit": self.singleton_batch_audit,
            "confounded_covariate_audit": self.confounded_covariate_audit,
            "r_runtime": self.r_runtime,
        }


def create_harmonizer(
    method: str,
    *,
    covariates: tuple[str, ...],
    parameters: dict[str, Any],
    device_spec: str,
    seed: int,
) -> Harmonizer:
    common = {
        "covariates": covariates,
        "parameters": parameters,
        "device_spec": device_spec,
        "seed": seed,
    }
    if method == "combat":
        return CombatHarmonizer(**common)
    if method == "combat_seq":
        return CombatSeqHarmonizer(**common)
    if method in {
        "mbatch_median_polish",
        "mbatch_empirical_bayes",
        "mbatch_anova",
    }:
        from .mbatch_harmonizer import MBatchHarmonizer

        return MBatchHarmonizer(method_id=method, **common)
    if method == "mober":
        from .mober_harmonizer import MoberHarmonizer

        return MoberHarmonizer(**common)
    raise ValueError(f"No dedicated harmonizer adapter for {method!r}")


def save_harmonizer(harmonizer: Harmonizer, directory: Path) -> Path:
    payload = harmonizer.save(directory)
    path = directory / "harmonizer.json"
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def load_harmonizer(directory: Path) -> Harmonizer | None:
    path = directory / "harmonizer.json"
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    method = str(payload["method"])
    if method == "combat":
        return CombatHarmonizer.load(directory, payload)
    if method == "combat_seq":
        return CombatSeqHarmonizer.load(directory, payload)
    if method in {
        "mbatch_median_polish",
        "mbatch_empirical_bayes",
        "mbatch_anova",
    }:
        from .mbatch_harmonizer import MBatchHarmonizer

        return MBatchHarmonizer.load(directory, payload)
    if method == "mober":
        from .mober_harmonizer import MoberHarmonizer

        return MoberHarmonizer.load(directory, payload)
    raise ValueError(f"Unknown serialized harmonizer: {method}")
