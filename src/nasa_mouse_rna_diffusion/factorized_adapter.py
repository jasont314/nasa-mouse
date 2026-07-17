"""Factorized residual conditioning for the pinned Lacan ModelDDIM."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np
import pandas as pd
import torch
from torch import nn


TISSUE_FALLBACKS = {
    "cecum": "colon",
    "cells": "cultured_cells",
    "eye": "retina",
    "optic_nerve": "retina",
}


def _decode(values: Iterable[object]) -> list[str]:
    return [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    ]


@dataclass(frozen=True)
class FactorizedSchema:
    """Names and slices of pretrained and residual conditioning features."""

    base_classes: tuple[str, ...]
    groups: dict[str, tuple[str, ...]]
    tissue_to_base: dict[str, str]

    @property
    def base_width(self) -> int:
        return len(self.base_classes)

    @property
    def adapter_width(self) -> int:
        return sum(len(values) for values in self.groups.values())

    @property
    def total_width(self) -> int:
        return self.base_width + self.adapter_width

    def group_slices(self) -> dict[str, slice]:
        offset = 0
        result: dict[str, slice] = {}
        for group, names in self.groups.items():
            result[group] = slice(offset, offset + len(names))
            offset += len(names)
        return result

    def as_dict(self) -> dict[str, object]:
        return {
            "base_classes": list(self.base_classes),
            "groups": {key: list(value) for key, value in self.groups.items()},
            "tissue_to_base": dict(self.tissue_to_base),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "FactorizedSchema":
        return cls(
            base_classes=tuple(map(str, payload["base_classes"])),
            groups={
                str(key): tuple(map(str, value))
                for key, value in payload["groups"].items()
            },
            tissue_to_base={
                str(key): str(value)
                for key, value in payload["tissue_to_base"].items()
            },
        )


def build_factorized_schema(
    train_samples: pd.DataFrame,
    base_classes: Iterable[str],
    *,
    include_study: bool = False,
    include_material_type: bool = False,
) -> FactorizedSchema:
    base = tuple(map(str, base_classes))
    base_set = set(base)
    tissues = tuple(sorted(train_samples["tissue"].dropna().astype(str).unique()))
    tissue_to_base: dict[str, str] = {}
    for tissue in tissues:
        mapped = tissue if tissue in base_set else TISSUE_FALLBACKS.get(tissue, "")
        if mapped not in base_set:
            raise ValueError(
                f"No pretrained tissue or declared fallback exists for {tissue!r}"
            )
        tissue_to_base[tissue] = mapped

    sexes = tuple(
        value
        for value in ("female", "male")
        if value in set(train_samples.get("sex", pd.Series(dtype=str)).astype(str))
    )
    muscle_groups = tuple(
        sorted(
            value
            for value in train_samples.get(
                "muscle_group", pd.Series(dtype=str)
            ).dropna().astype(str).unique()
            if value not in {"", "unknown", "nan", "none"}
        )
    )
    studies = (
        tuple(sorted(train_samples["accession"].dropna().astype(str).unique()))
        if include_study
        else ()
    )
    material_types = (
        tuple(
            sorted(
                value
                for value in train_samples.get(
                    "material_type", pd.Series(dtype=str)
                ).dropna().astype(str).unique()
                if value not in {"", "unknown", "nan", "none"}
            )
        )
        if include_material_type
        else ()
    )
    condition_values = ("ground_control", "flight")
    domain = (
        ("domain=osdr",)
        + tuple(f"tissue_residual={value}" for value in tissues)
        + tuple(f"sex={value}" for value in sexes)
        + tuple(f"muscle_group={value}" for value in muscle_groups)
        + tuple(f"study={value}" for value in studies)
        + tuple(f"material_type={value}" for value in material_types)
    )
    condition = (
        tuple(f"condition={value}" for value in condition_values)
        + tuple(
            f"tissue_condition={tissue}::{condition_value}"
            for tissue in tissues
            for condition_value in condition_values
        )
        + tuple(
            f"muscle_condition={muscle_group}::{condition_value}"
            for muscle_group in muscle_groups
            for condition_value in condition_values
        )
        + tuple(
            f"study_condition={study}::{condition_value}"
            for study in studies
            for condition_value in condition_values
        )
    )
    return FactorizedSchema(
        base_classes=base,
        groups={"domain": domain, "condition": condition},
        tissue_to_base=tissue_to_base,
    )


def encode_factorized_labels(
    samples: pd.DataFrame, schema: FactorizedSchema
) -> np.ndarray:
    base_map = {name: index for index, name in enumerate(schema.base_classes)}
    labels = np.zeros((len(samples), schema.total_width), dtype=np.int64)
    tissues = samples["tissue"].astype(str).to_numpy()
    for row, tissue in enumerate(tissues):
        mapped = schema.tissue_to_base.get(tissue)
        if mapped is None:
            raise ValueError(f"Tissue {tissue!r} was absent from the training schema")
        labels[row, base_map[mapped]] = 1

    adapter_offset = schema.base_width
    for group, names in schema.groups.items():
        for local_index, name in enumerate(names):
            column = adapter_offset + local_index
            if name == "domain=osdr":
                labels[:, column] = 1
            elif name.startswith("tissue_residual="):
                value = name.split("=", 1)[1]
                labels[:, column] = tissues == value
            elif name.startswith("sex="):
                value = name.split("=", 1)[1]
                labels[:, column] = (
                    samples.get("sex", pd.Series("unknown", index=samples.index))
                    .astype(str)
                    .to_numpy()
                    == value
                )
            elif name.startswith("muscle_group="):
                value = name.split("=", 1)[1]
                labels[:, column] = (
                    samples.get(
                        "muscle_group", pd.Series("unknown", index=samples.index)
                    )
                    .astype(str)
                    .to_numpy()
                    == value
                )
            elif name.startswith("study="):
                value = name.split("=", 1)[1]
                labels[:, column] = (
                    samples["accession"].astype(str).to_numpy() == value
                )
            elif name.startswith("material_type="):
                value = name.split("=", 1)[1]
                labels[:, column] = (
                    samples.get(
                        "material_type", pd.Series("unknown", index=samples.index)
                    )
                    .astype(str)
                    .to_numpy()
                    == value
                )
            elif name.startswith("condition="):
                value = name.split("=", 1)[1]
                labels[:, column] = samples["condition"].astype(str).to_numpy() == value
            elif name.startswith("tissue_condition="):
                tissue, condition = name.split("=", 1)[1].split("::", 1)
                labels[:, column] = (tissues == tissue) & (
                    samples["condition"].astype(str).to_numpy() == condition
                )
            elif name.startswith("muscle_condition="):
                muscle_group, condition = name.split("=", 1)[1].split("::", 1)
                observed_group = samples.get(
                    "muscle_group", pd.Series("unknown", index=samples.index)
                ).astype(str).to_numpy()
                labels[:, column] = (observed_group == muscle_group) & (
                    samples["condition"].astype(str).to_numpy() == condition
                )
            elif name.startswith("study_condition="):
                study, condition = name.split("=", 1)[1].split("::", 1)
                labels[:, column] = (
                    samples["accession"].astype(str).to_numpy() == study
                ) & (samples["condition"].astype(str).to_numpy() == condition)
            else:
                raise ValueError(f"Unsupported factorized feature: {name}")
        adapter_offset += len(names)
    return labels


def load_factorized_role(
    prepared_h5: str | Path,
    samples_tsv: str | Path,
    role: str,
) -> dict[str, object]:
    """Read exactly one declared role; test data are never loaded implicitly."""

    samples = pd.read_csv(samples_tsv, sep="\t")
    role_samples = samples.loc[samples["role"].astype(str).eq(role)].copy()
    with h5py.File(prepared_h5, "r") as handle:
        if role not in handle:
            raise ValueError(f"Prepared data has no role named {role!r}")
        group = handle[role]
        expression = np.asarray(group["expression"][:], dtype=np.float32)
        analysis_key = (
            "analysis_expression" if "analysis_expression" in group else "tpm"
        )
        analysis_expression = np.asarray(group[analysis_key][:], dtype=np.float32)
        source_rows = np.asarray(group["source_row"][:], dtype=np.int64)
        lookup = role_samples.set_index("_row_index", drop=False)
        aligned = lookup.loc[source_rows].reset_index(drop=True)
        genes = _decode(handle["genes"][:])
        scale = np.asarray(handle["maxabs_scale"][:], dtype=np.float32)
        analysis_units = str(handle.attrs.get("analysis_units", "tpm"))
    if len(aligned) != len(expression):
        raise ValueError("Expression and metadata do not align")
    return {
        "expression": expression,
        "analysis_expression": analysis_expression,
        "source_row": source_rows,
        "samples": aligned,
        "genes": genes,
        "maxabs_scale": scale,
        "analysis_units": analysis_units,
    }


def _swish(values: torch.Tensor) -> torch.Tensor:
    return values * torch.sigmoid(values)


class _GroupedProjection(nn.Module):
    def __init__(
        self,
        group_widths: dict[str, int],
        embedding_dim: int,
        output_dim: int,
    ) -> None:
        super().__init__()
        self.layers = nn.ModuleDict(
            {
                group: nn.Linear(width * embedding_dim, output_dim, bias=False)
                for group, width in group_widths.items()
                if width
            }
        )
        for layer in self.layers.values():
            nn.init.zeros_(layer.weight)

    def forward(self, inputs: dict[str, torch.Tensor]) -> torch.Tensor:
        result = None
        for group, layer in self.layers.items():
            projected = layer(inputs[group])
            result = projected if result is None else result + projected
        if result is None:
            raise RuntimeError("Factorized adapter has no feature groups")
        return result


class FactorizedAdapterDDIM(nn.Module):
    """Frozen upstream DDIM with residual projections for factorized metadata."""

    def __init__(self, base_model: nn.Module, schema: FactorizedSchema) -> None:
        super().__init__()
        self.base_model = base_model
        self.schema = schema
        for parameter in self.base_model.parameters():
            parameter.requires_grad_(False)
        embedding_dim = int(base_model.dim_y_emb)
        widths = {key: len(value) for key, value in schema.groups.items()}
        self.w1_adapters = nn.ModuleList(
            [
                _GroupedProjection(widths, embedding_dim, block.out_channels)
                for block in base_model.mid
            ]
        )
        self.w2_adapters = nn.ModuleList(
            [
                _GroupedProjection(widths, embedding_dim, block.out_channels)
                for block in base_model.mid
            ]
        )
        self.x_adapters = nn.ModuleDict(
            {
                str(index): _GroupedProjection(
                    widths, embedding_dim, block.out_channels
                )
                for index, block in enumerate(base_model.mid)
                if block.in_channels != block.out_channels
            }
        )
        self.set_trainable_groups(())

    def set_trainable_groups(self, groups: Iterable[str]) -> None:
        selected = set(map(str, groups))
        unknown = selected.difference(self.schema.groups)
        if unknown:
            raise ValueError(f"Unknown adapter groups: {sorted(unknown)}")
        for name, parameter in self.named_parameters():
            if name.startswith("base_model."):
                parameter.requires_grad_(False)
                continue
            parameter.requires_grad_(
                any(f"layers.{group}." in name for group in selected)
            )

    def trainable_parameter_count(self) -> int:
        return sum(
            parameter.numel()
            for parameter in self.parameters()
            if parameter.requires_grad
        )

    def adapter_state_dict(self) -> dict[str, torch.Tensor]:
        return {
            name: value.detach().cpu()
            for name, value in self.state_dict().items()
            if not name.startswith("base_model.")
        }

    def load_adapter_state_dict(self, state: dict[str, torch.Tensor]) -> None:
        observed = self.state_dict()
        observed.update(state)
        self.load_state_dict(observed)

    def _embedded_labels(
        self, labels: torch.Tensor
    ) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
        if labels.ndim != 2 or labels.shape[1] != self.schema.total_width:
            raise ValueError(
                f"Expected labels with width {self.schema.total_width}, "
                f"observed {tuple(labels.shape)}"
            )
        base_labels = labels[:, : self.schema.base_width].long()
        adapter_labels = labels[:, self.schema.base_width :].long()
        base_embedding = _swish(self.base_model.y_emb(base_labels)).reshape(
            len(labels), -1
        )
        adapter_embedding = _swish(
            self.base_model.y_emb(adapter_labels)
        )
        inputs: dict[str, torch.Tensor] = {}
        for group, bounds in self.schema.group_slices().items():
            inputs[group] = adapter_embedding[:, bounds, :].reshape(len(labels), -1)
        return base_embedding, inputs

    def forward(
        self, x: torch.Tensor, timesteps: torch.Tensor, labels: torch.Tensor
    ) -> torch.Tensor:
        base = self.base_model
        if base.is_time_embed:
            raise ValueError("The pinned paper model uses scalar timestep conditioning")
        temporal = (timesteps / base.timestep_max).unsqueeze(1)
        base_labels, adapter_inputs = self._embedded_labels(labels)
        hidden = x
        for index, block in enumerate(base.mid):
            projected = torch.cat((hidden, temporal), dim=1)
            projected = block.batch_norm_temb1(projected)
            projected = block.w1(torch.cat((projected, base_labels), dim=1))
            projected = projected + self.w1_adapters[index](adapter_inputs)
            projected = block.batch_norm1(projected)
            projected = block.dropout(block.relu(projected))
            projected = torch.cat((projected, temporal), dim=1)
            projected = block.batch_norm_temb2(projected)
            projected = block.w2(torch.cat((projected, base_labels), dim=1))
            projected = projected + self.w2_adapters[index](adapter_inputs)
            projected = block.batch_norm2(projected)
            projected = block.dropout(block.relu(projected))
            residual = hidden
            if block.in_channels != block.out_channels:
                residual = block.x_proj(torch.cat((hidden, base_labels), dim=1))
                residual = residual + self.x_adapters[str(index)](adapter_inputs)
            hidden = residual + projected
        hidden = _swish(base.norm_out(hidden))
        return base.lin_out(hidden)


def neutralize_group(
    labels: torch.Tensor, schema: FactorizedSchema, group: str
) -> torch.Tensor:
    result = labels.clone()
    bounds = schema.group_slices()[group]
    start = schema.base_width + int(bounds.start)
    stop = schema.base_width + int(bounds.stop)
    result[:, start:stop] = 0
    return result
