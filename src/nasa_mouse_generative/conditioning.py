"""Canonical categorical conditioning for OSDR and ARCHS4 profiles."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd


MISSING = "__missing__"
UNKNOWN = "__unknown__"


def _column(frame: pd.DataFrame, candidates: Iterable[str], default: str) -> pd.Series:
    for candidate in candidates:
        if candidate in frame:
            return frame[candidate].fillna(default).astype(str)
    return pd.Series(default, index=frame.index, dtype="object")


def _clean(values: pd.Series, default: str = MISSING) -> pd.Series:
    result = values.fillna(default).astype(str).str.strip()
    return result.mask(result.eq(""), default)


def infer_muscle_group(material: object, tissue: object) -> str:
    text = re.sub(r"[^a-z0-9]+", " ", str(material).lower()).strip()
    if "soleus" in text:
        return "soleus"
    if "gastrocnemius" in text:
        return "gastrocnemius"
    if "quadriceps" in text:
        return "quadriceps"
    if "tibialis" in text:
        return "tibialis_anterior"
    if "extensor digitorum" in text or re.search(r"\bedl\b", text):
        return "edl"
    if str(tissue) == "skeletal_muscle":
        return "unspecified_skeletal_muscle"
    return "not_applicable"


def osdr_conditioning_frame(obs: pd.DataFrame) -> pd.DataFrame:
    result = pd.DataFrame(index=obs.index.copy())
    result["profile_id"] = _column(
        obs, ("biological_profile_id", "profile_id"), "unknown_profile"
    )
    result["accession"] = _column(obs, ("id.accession",), "unknown_accession")
    result["condition"] = _column(
        obs, ("condition_inferred",), "unknown_condition"
    ).str.lower()
    result["tissue"] = _column(
        obs, ("tissue_canonical", "tissue_final"), "unknown_tissue"
    ).str.lower()
    result["material_type"] = _column(
        obs,
        ("study.characteristics.material type", "material_type_original"),
        "unknown_material",
    )
    result["study"] = result["accession"]
    result["sex"] = _column(
        obs, ("study.characteristics.sex",), "unknown_sex"
    ).str.lower()
    result["age"] = _column(
        obs,
        (
            "study.characteristics.age",
            "study.characteristics.age at launch",
            "age",
        ),
        "unknown_age",
    ).str.lower()
    result["assay"] = _column(
        obs,
        (
            "investigation.study assays.study assay technology type",
            "id.assay name",
            "file.datatype",
        ),
        "bulk_rna_seq",
    ).str.lower()
    result["platform"] = _column(
        obs, ("platform", "instrument_model"), "unknown_platform"
    ).str.lower()
    result["data_source"] = _column(
        obs,
        (
            "investigation.study.comment.project type",
            "investigation.study.comment.data source accession",
        ),
        "nasa_osdr_api",
    ).str.lower()
    if "muscle_group" in obs:
        result["muscle_group"] = _clean(obs["muscle_group"])
    else:
        result["muscle_group"] = [
            infer_muscle_group(material, tissue)
            for material, tissue in zip(result["material_type"], result["tissue"])
        ]
    result["source"] = "osdr"
    for column in result.columns:
        result[column] = _clean(result[column])
    return result


def archs4_conditioning_frame(metadata: pd.DataFrame) -> pd.DataFrame:
    result = pd.DataFrame(index=metadata.index.copy())
    result["profile_id"] = _column(
        metadata, ("geo_accession",), "unknown_archs4_profile"
    )
    result["accession"] = _column(
        metadata, ("series_id",), "unknown_archs4_series"
    )
    result["condition"] = "archs4_reference"
    result["tissue"] = _column(
        metadata, ("canonical_tissue",), "unknown_tissue"
    ).str.lower()
    result["material_type"] = _column(
        metadata, ("source_name_ch1",), "unknown_material"
    )
    result["study"] = result["accession"]
    characteristics = _column(metadata, ("characteristics_ch1",), "")
    sex = characteristics.str.extract(
        r"(?:sex|gender)\s*[:=]\s*([^,;|]+)", flags=re.IGNORECASE, expand=False
    )
    result["sex"] = sex.fillna("unknown_sex").str.lower()
    age = characteristics.str.extract(
        r"(?:age|developmental stage)\s*[:=]\s*([^,;|]+)",
        flags=re.IGNORECASE,
        expand=False,
    )
    result["age"] = age.fillna("unknown_age").str.lower()
    result["assay"] = _column(
        metadata, ("library_strategy",), "rna_seq"
    ).str.lower()
    result["platform"] = _column(
        metadata, ("instrument_model",), "unknown_platform"
    ).str.lower()
    result["data_source"] = "archs4"
    result["muscle_group"] = [
        infer_muscle_group(material, tissue)
        for material, tissue in zip(result["material_type"], result["tissue"])
    ]
    result["source"] = "archs4"
    for column in result.columns:
        result[column] = _clean(result[column])
    return result


@dataclass
class CategoryEncoder:
    covariates: tuple[str, ...]
    vocabularies: dict[str, list[str]]

    @classmethod
    def fit(
        cls, frames: Iterable[pd.DataFrame], covariates: Iterable[str]
    ) -> "CategoryEncoder":
        frame_list = list(frames)
        names = tuple(map(str, covariates))
        vocabularies: dict[str, list[str]] = {}
        for covariate in names:
            values = {MISSING, UNKNOWN}
            for frame in frame_list:
                if covariate not in frame:
                    raise ValueError(f"Conditioning frame lacks {covariate!r}")
                values.update(_clean(frame[covariate]).tolist())
            vocabularies[covariate] = [MISSING, UNKNOWN] + sorted(
                values.difference({MISSING, UNKNOWN})
            )
        return cls(covariates=names, vocabularies=vocabularies)

    @property
    def cardinalities(self) -> list[int]:
        return [len(self.vocabularies[name]) for name in self.covariates]

    def transform(self, frame: pd.DataFrame) -> np.ndarray:
        if not self.covariates:
            return np.empty((len(frame), 0), dtype=np.int64)
        columns = []
        for covariate in self.covariates:
            if covariate not in frame:
                values = pd.Series(MISSING, index=frame.index)
            else:
                values = _clean(frame[covariate])
            vocabulary = self.vocabularies[covariate]
            mapping = {value: index for index, value in enumerate(vocabulary)}
            columns.append(
                values.map(mapping).fillna(mapping[UNKNOWN]).to_numpy(dtype=np.int64)
            )
        return np.stack(columns, axis=1)

    def encode_profiles(self, profiles: list[dict[str, str]]) -> np.ndarray:
        return self.transform(pd.DataFrame(profiles))

    def save(self, path: str | Path) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "covariates": list(self.covariates),
            "vocabularies": self.vocabularies,
            "unknown_policy": "map unseen values to __unknown__",
        }
        output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        return output

    @classmethod
    def load(cls, path: str | Path) -> "CategoryEncoder":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(
            covariates=tuple(payload["covariates"]),
            vocabularies={
                str(key): list(map(str, values))
                for key, values in payload["vocabularies"].items()
            },
        )
