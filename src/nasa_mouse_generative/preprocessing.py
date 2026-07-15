"""Fold-aware preprocessing primitives shared by model adapters."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .config import PreprocessingConfig
from .harmonizers import (
    Harmonizer,
    create_harmonizer,
    load_harmonizer,
    save_harmonizer,
)


EPSILON = 1e-8


@dataclass
class ScaleStats:
    center: np.ndarray
    scale: np.ndarray


@dataclass
class FittedPreprocessor:
    spec: PreprocessingConfig
    device_spec: str = "auto"
    seed: int = 0
    global_stats: ScaleStats | None = None
    study_stats: dict[str, ScaleStats] = field(default_factory=dict)
    post_harmonization_stats: ScaleStats | None = None
    final_stats: ScaleStats | None = None
    harmonizer: Harmonizer | None = field(default=None, repr=False)

    def fit_transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        *,
        gene_lengths: np.ndarray | None = None,
        metadata: pd.DataFrame | None = None,
    ) -> np.ndarray:
        studies_array = np.asarray(list(studies), dtype=str)
        raw_values = np.asarray(matrix)
        if raw_values.shape[0] != studies_array.shape[0]:
            raise ValueError("studies must contain one value per sample")

        if self.spec.harmonization in {"combat", "combat_seq", "mober"}:
            self.harmonizer = create_harmonizer(
                self.spec.harmonization,
                covariates=self.spec.harmonization_covariates,
                parameters=self.spec.harmonization_parameters,
                device_spec=self.device_spec,
                seed=self.seed,
            )
            if self.harmonizer.requires_raw_counts:
                harmonized_raw = self.harmonizer.fit_transform(
                    raw_values, studies_array, metadata
                )
                values = self._base_transform(
                    harmonized_raw, gene_lengths=gene_lengths
                )
            else:
                values = self._base_transform(
                    raw_values, gene_lengths=gene_lengths
                )
                values = self.harmonizer.fit_transform(
                    values, studies_array, metadata
                )
            return self._fit_final_scaler(values)

        values = self._base_transform(raw_values, gene_lengths=gene_lengths)

        if self.spec.harmonization in {
            "within_study_zscore",
            "within_study_then_global_zscore",
        }:
            self.global_stats = fit_stats(values, "zscore")
            harmonized = np.empty_like(values, dtype=np.float32)
            for study in sorted(set(studies_array)):
                mask = studies_array == study
                stats = fit_stats(values[mask], "zscore")
                self.study_stats[study] = stats
                harmonized[mask] = apply_stats(values[mask], stats)
            if self.spec.harmonization == "within_study_then_global_zscore":
                self.post_harmonization_stats = fit_stats(harmonized, "zscore")
                harmonized = apply_stats(harmonized, self.post_harmonization_stats)
            return self._fit_final_scaler(harmonized)

        if self.spec.harmonization != "none":
            raise ValueError(f"Unsupported harmonization: {self.spec.harmonization}")
        return self._fit_final_scaler(values)

    def transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        *,
        gene_lengths: np.ndarray | None = None,
        allow_transductive: bool = False,
        metadata: pd.DataFrame | None = None,
    ) -> np.ndarray:
        studies_array = np.asarray(list(studies), dtype=str)
        raw_values = np.asarray(matrix)
        if raw_values.shape[0] != studies_array.shape[0]:
            raise ValueError("studies must contain one value per sample")
        if self.spec.harmonization in {"combat", "combat_seq", "mober"}:
            if self.harmonizer is None:
                raise RuntimeError("Dedicated harmonizer is not fitted")
            if self.harmonizer.requires_raw_counts:
                harmonized_raw = self.harmonizer.transform(
                    raw_values,
                    studies_array,
                    metadata,
                    allow_transductive=allow_transductive,
                )
                values = self._base_transform(
                    harmonized_raw, gene_lengths=gene_lengths
                )
            else:
                values = self._base_transform(
                    raw_values, gene_lengths=gene_lengths
                )
                values = self.harmonizer.transform(
                    values,
                    studies_array,
                    metadata,
                    allow_transductive=allow_transductive,
                )
            return self._apply_final_scaler(values)

        values = self._base_transform(raw_values, gene_lengths=gene_lengths)
        if self.spec.harmonization in {
            "within_study_zscore",
            "within_study_then_global_zscore",
        }:
            if self.global_stats is None:
                raise RuntimeError("Preprocessor is not fitted")
            harmonized = np.empty_like(values, dtype=np.float32)
            for study in sorted(set(studies_array)):
                mask = studies_array == study
                stats = self.study_stats.get(study)
                if stats is None:
                    if self.spec.unseen_study_policy == "transductive_unlabeled":
                        if not allow_transductive:
                            raise ValueError(
                                "Unseen-study transductive scaling was requested but not allowed."
                            )
                        stats = fit_stats(values[mask], "zscore")
                    else:
                        stats = self.global_stats
                harmonized[mask] = apply_stats(values[mask], stats)
            if self.post_harmonization_stats is not None:
                harmonized = apply_stats(harmonized, self.post_harmonization_stats)
            return self._apply_final_scaler(harmonized)
        return self._apply_final_scaler(values)

    def fit_additional_study_stats(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        *,
        gene_lengths: np.ndarray | None = None,
    ) -> None:
        """Fit harmonization statistics for additional training-only studies.

        This is used after reference fitting so OSDR fine-tuning accessions can
        receive training-fitted study transforms while unseen validation and test
        accessions still use the reference-global fallback.
        """

        if self.spec.harmonization not in {
            "within_study_zscore",
            "within_study_then_global_zscore",
        }:
            return
        studies_array = np.asarray(list(studies), dtype=str)
        values = self._base_transform(matrix, gene_lengths=gene_lengths)
        if len(studies_array) != len(values):
            raise ValueError("studies must contain one value per sample")
        for study in sorted(set(studies_array)):
            mask = studies_array == study
            self.study_stats[study] = fit_stats(values[mask], "zscore")

    def inverse_transform(
        self, values: np.ndarray, studies: Iterable[object]
    ) -> np.ndarray:
        """Return generated values in the configured normalized input units.

        Library-size normalization is not invertible without a target library
        size. Consequently CPM and TPM runs return CPM/TPM, while ``none`` returns
        the original input units after reversing transform/scaling.
        """

        result = np.asarray(values, dtype=np.float32)
        studies_array = np.asarray(list(studies), dtype=str)
        if result.ndim != 2 or len(studies_array) != len(result):
            raise ValueError("values and studies must be aligned samples x genes")
        if self.final_stats is not None:
            result = invert_stats(result, self.final_stats)
        if self.post_harmonization_stats is not None:
            result = invert_stats(result, self.post_harmonization_stats)
        if self.spec.harmonization in {
            "within_study_zscore",
            "within_study_then_global_zscore",
        }:
            if self.global_stats is None:
                raise RuntimeError("Preprocessor is not fitted")
            restored = np.empty_like(result, dtype=np.float32)
            for study in sorted(set(studies_array)):
                mask = studies_array == study
                stats = self.study_stats.get(study, self.global_stats)
                restored[mask] = invert_stats(result[mask], stats)
            result = restored
        if self.spec.transform == "log1p":
            result = np.expm1(np.clip(result, -30.0, 30.0))
        elif self.spec.transform == "log2p1":
            result = np.exp2(np.clip(result, -30.0, 30.0)) - 1.0
        result = np.nan_to_num(result, nan=0.0, posinf=np.finfo(np.float32).max)
        if self.spec.input_units == "raw_counts" or self.spec.library_normalization != "none":
            result = np.maximum(result, 0.0)
        if self.spec.library_normalization in {"cpm", "tpm"}:
            totals = result.sum(axis=1, keepdims=True, dtype=np.float64)
            result = result / np.maximum(totals, EPSILON) * 1_000_000.0
        return result.astype(np.float32)

    @property
    def output_units(self) -> str:
        return (
            self.spec.library_normalization
            if self.spec.library_normalization != "none"
            else self.spec.input_units
        )

    def save(self, directory: str | Path) -> tuple[Path, Path]:
        output_dir = Path(directory)
        output_dir.mkdir(parents=True, exist_ok=True)
        arrays: dict[str, np.ndarray] = {}
        metadata: dict[str, object] = {
            "spec": asdict(self.spec),
            "output_units": self.output_units,
            "runtime": {"device_spec": self.device_spec, "seed": self.seed},
            "stats": {},
            "study_stats": {},
        }

        def add_stats(label: str, stats: ScaleStats | None) -> None:
            if stats is None:
                return
            arrays[f"{label}_center"] = stats.center
            arrays[f"{label}_scale"] = stats.scale
            metadata["stats"][label] = True

        add_stats("global", self.global_stats)
        add_stats("post_harmonization", self.post_harmonization_stats)
        add_stats("final", self.final_stats)
        for index, (study, stats) in enumerate(sorted(self.study_stats.items())):
            label = f"study_{index:05d}"
            arrays[f"{label}_center"] = stats.center
            arrays[f"{label}_scale"] = stats.scale
            metadata["study_stats"][study] = label
        arrays_path = output_dir / "preprocessing_stats.npz"
        metadata_path = output_dir / "preprocessing.json"
        if self.harmonizer is not None:
            save_harmonizer(self.harmonizer, output_dir)
            metadata["harmonizer"] = self.harmonizer.audit()
        np.savez_compressed(arrays_path, **arrays)
        metadata_path.write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )
        return metadata_path, arrays_path

    @classmethod
    def load(cls, directory: str | Path) -> "FittedPreprocessor":
        input_dir = Path(directory)
        metadata = json.loads(
            (input_dir / "preprocessing.json").read_text(encoding="utf-8")
        )
        arrays = np.load(input_dir / "preprocessing_stats.npz")
        spec_options = dict(metadata["spec"])
        if "harmonization_covariates" in spec_options:
            spec_options["harmonization_covariates"] = tuple(
                spec_options["harmonization_covariates"]
            )
        runtime = metadata.get("runtime", {})
        processor = cls(
            PreprocessingConfig(**spec_options),
            device_spec=str(runtime.get("device_spec", "auto")),
            seed=int(runtime.get("seed", 0)),
        )

        def read_stats(label: str) -> ScaleStats | None:
            center_key = f"{label}_center"
            if center_key not in arrays:
                return None
            return ScaleStats(
                center=np.asarray(arrays[center_key], dtype=np.float32),
                scale=np.asarray(arrays[f"{label}_scale"], dtype=np.float32),
            )

        processor.global_stats = read_stats("global")
        processor.post_harmonization_stats = read_stats("post_harmonization")
        processor.final_stats = read_stats("final")
        processor.study_stats = {
            study: read_stats(label)
            for study, label in metadata.get("study_stats", {}).items()
        }
        if any(stats is None for stats in processor.study_stats.values()):
            raise ValueError("Serialized preprocessing study statistics are incomplete")
        processor.harmonizer = load_harmonizer(input_dir)
        return processor

    def audit(self) -> dict[str, object]:
        return {
            "method": self.spec.harmonization,
            "covariates": list(self.spec.harmonization_covariates),
            "outcome_informed": (
                self.harmonizer.audit().get("outcome_informed", False)
                if self.harmonizer is not None
                else False
            ),
            "adapter": self.harmonizer.audit() if self.harmonizer is not None else None,
        }

    def _base_transform(
        self, matrix: np.ndarray, *, gene_lengths: np.ndarray | None
    ) -> np.ndarray:
        values = np.asarray(matrix, dtype=np.float64)
        if values.ndim != 2:
            raise ValueError("expression matrix must be samples x genes")
        if not np.isfinite(values).all():
            raise ValueError("expression matrix contains non-finite values")
        if (values < 0).any() and self.spec.input_units == "raw_counts":
            raise ValueError("raw counts cannot be negative")

        normalization = self.spec.library_normalization
        if normalization == "none":
            pass
        elif normalization == "cpm":
            library = values.sum(axis=1, keepdims=True)
            values = values / np.maximum(library, EPSILON) * 1_000_000.0
        elif normalization == "tpm":
            if gene_lengths is None:
                raise ValueError("TPM normalization requires gene lengths")
            lengths = np.asarray(gene_lengths, dtype=np.float64).reshape(1, -1)
            if lengths.shape[1] != values.shape[1] or (lengths <= 0).any():
                raise ValueError("gene lengths must be positive and aligned to matrix columns")
            rates = values / (lengths / 1000.0)
            values = rates / np.maximum(rates.sum(axis=1, keepdims=True), EPSILON) * 1_000_000.0
        else:
            raise ValueError(f"Unsupported library normalization: {normalization}")

        if self.spec.transform == "log1p":
            values = np.log1p(values)
        elif self.spec.transform == "log2p1":
            values = np.log2(values + 1.0)
        elif self.spec.transform != "none":
            raise ValueError(f"Unsupported transform: {self.spec.transform}")
        return values.astype(np.float32)

    def _fit_final_scaler(self, values: np.ndarray) -> np.ndarray:
        if self.spec.scaler == "none":
            return values.astype(np.float32, copy=False)
        stats = fit_stats(values, self.spec.scaler)
        # Keep final scaling separate from the mentor method's second z-score.
        self.final_stats = stats
        return apply_stats(values, stats)

    def _apply_final_scaler(self, values: np.ndarray) -> np.ndarray:
        if self.spec.scaler == "none":
            return values.astype(np.float32, copy=False)
        if self.final_stats is None:
            raise RuntimeError("Preprocessor is not fitted")
        return apply_stats(values, self.final_stats)


def fit_stats(values: np.ndarray, method: str) -> ScaleStats:
    array = np.asarray(values, dtype=np.float64)
    if method == "zscore":
        center = array.mean(axis=0)
        scale = array.std(axis=0)
    elif method == "global_zscore":
        center = np.asarray([array.mean()], dtype=np.float64)
        scale = np.asarray([array.std()], dtype=np.float64)
    elif method == "robust":
        center = np.median(array, axis=0)
        q25, q75 = np.percentile(array, [25, 75], axis=0)
        scale = q75 - q25
    elif method == "maxabs":
        center = np.zeros(array.shape[1], dtype=np.float64)
        scale = np.max(np.abs(array), axis=0)
    else:
        raise ValueError(f"Unsupported scaler: {method}")
    scale = np.where(np.abs(scale) < EPSILON, 1.0, scale)
    return ScaleStats(center=center.astype(np.float32), scale=scale.astype(np.float32))


def apply_stats(values: np.ndarray, stats: ScaleStats) -> np.ndarray:
    return ((np.asarray(values, dtype=np.float32) - stats.center) / stats.scale).astype(
        np.float32
    )


def invert_stats(values: np.ndarray, stats: ScaleStats) -> np.ndarray:
    return (
        np.asarray(values, dtype=np.float32) * stats.scale + stats.center
    ).astype(np.float32)
