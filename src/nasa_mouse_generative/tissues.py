"""Shared tissue ontology for OSDR and ARCHS4 metadata."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Iterable


@dataclass(frozen=True)
class TissueRule:
    canonical: str
    aliases: tuple[str, ...]
    material_class: str = "tissue"

    @property
    def pattern(self) -> re.Pattern[str]:
        parts = [r"(?<![a-z0-9])" + re.escape(alias) + r"(?![a-z0-9])" for alias in self.aliases]
        return re.compile("|".join(parts), flags=re.IGNORECASE)


# Specific anatomical structures precede their parent organs.
TISSUE_RULES = (
    TissueRule("optic_nerve", ("optic nerve",)),
    TissueRule("retina", ("retina", "retinal")),
    TissueRule("cerebellum", ("cerebellum", "cerebellar")),
    TissueRule("hippocampus", ("hippocampus", "hippocampal")),
    TissueRule("bone_marrow", ("bone marrow",)),
    TissueRule("brown_adipose_tissue", ("brown adipose", "brown fat")),
    TissueRule("white_adipose_tissue", ("white adipose", "white fat")),
    TissueRule(
        "skeletal_muscle",
        (
            "skeletal muscle",
            "soleus",
            "gastrocnemius",
            "quadriceps",
            "tibialis anterior",
            "extensor digitorum longus",
            "edl muscle",
        ),
    ),
    TissueRule("adrenal_gland", ("adrenal gland", "adrenal")),
    TissueRule("mammary_gland", ("mammary gland", "mammary")),
    TissueRule("cecum", ("cecum", "caecum")),
    TissueRule("colon", ("colon", "colonic")),
    TissueRule("liver", ("liver", "hepatic", "hepatocyte")),
    TissueRule("kidney", ("kidney", "renal", "nephron")),
    TissueRule("spleen", ("spleen", "splenic")),
    TissueRule("thymus", ("thymus", "thymic")),
    TissueRule("lung", ("lung", "pulmonary")),
    TissueRule("skin", ("skin", "dermal", "epidermal", "cutaneous")),
    TissueRule("heart", ("heart", "cardiac", "myocardium", "ventricle")),
    TissueRule("brain", ("brain", "cerebrum", "cerebral hemisphere")),
    TissueRule("eye", ("eye", "ocular")),
    TissueRule("bone", ("bone", "mandible")),
    TissueRule(
        "cultured_cells",
        ("cultured cells", "cells cultured", "cell culture", "cell line"),
        "cell_model",
    ),
)

RULE_BY_NAME = {rule.canonical: rule for rule in TISSUE_RULES}
NONMATCHABLE_OSDR_CLASSES = {"cells"}


def normalize_text(value: object) -> str:
    text = re.sub(r"[^a-z0-9]+", " ", str(value or "").lower())
    return " ".join(text.split())


def canonicalize_material(value: object) -> str:
    text = normalize_text(value)
    if not text:
        return "unknown"
    if text == "cells":
        return "cells"
    for rule in TISSUE_RULES:
        if rule.pattern.search(text):
            return rule.canonical
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_") or "unknown"


def rules_for_tissues(tissues: Iterable[str]) -> tuple[TissueRule, ...]:
    wanted = set(map(str, tissues))
    return tuple(rule for rule in TISSUE_RULES if rule.canonical in wanted)
