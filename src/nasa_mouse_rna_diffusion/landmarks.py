"""Build a deterministic 974-gene mouse analogue of the GTEx L1000 panel."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


ORTHOLOGY_RANK = {
    "ortholog_one2one": 0,
    "ortholog_one2many": 1,
    "ortholog_many2many": 2,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _decode(values: np.ndarray) -> list[str]:
    return [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def _archs4_genes(path: str | Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        genes = _decode(handle["meta/genes/ensembl_gene"][:])
        symbols = _decode(handle["meta/genes/symbol"][:])
    return pd.DataFrame({"mouse_ensembl_gene": genes, "mouse_symbol": symbols})


def build_mouse_landmark_panel(
    *,
    human_landmarks: str | Path,
    source_map: str | Path,
    archs4_h5: str | Path,
    output: str | Path,
    dimensions: int = 974,
) -> tuple[Path, dict[str, Any]]:
    human_path = Path(human_landmarks)
    map_path = Path(source_map)
    archs4_path = Path(archs4_h5)
    output_path = Path(output)
    human = pd.read_csv(human_path)
    required = {"Description", "ensembl_id", "Type"}
    if not required.issubset(human.columns):
        raise ValueError(f"Human landmark table lacks {sorted(required)}")
    human = human.loc[human["Type"].astype(str).str.lower().eq("landmark")].copy()
    human = human[["Description", "ensembl_id"]].rename(
        columns={"Description": "human_symbol", "ensembl_id": "human_ensembl_gene"}
    )
    if len(human) != dimensions or human["human_ensembl_gene"].nunique() != dimensions:
        raise ValueError(
            f"Expected {dimensions} unique GTEx landmarks, observed {len(human)} rows"
        )

    mapped = pd.read_csv(map_path, sep="\t")
    required_map = {
        "human_symbol",
        "human_ensembl_gene",
        "mouse_ensembl_gene",
        "mouse_symbol",
        "orthology_type",
        "mouse_perc_id",
    }
    if not required_map.issubset(mapped.columns):
        raise ValueError(f"Landmark map lacks {sorted(required_map)}")
    mapped = mapped.copy()
    mapped["_orthology_rank"] = (
        mapped["orthology_type"].map(ORTHOLOGY_RANK).fillna(99).astype(int)
    )
    mapped["_identity"] = pd.to_numeric(mapped["mouse_perc_id"], errors="coerce").fillna(-1)
    mapped = mapped.sort_values(
        ["human_ensembl_gene", "_orthology_rank", "_identity", "mouse_ensembl_gene"],
        ascending=[True, True, False, True],
        kind="stable",
    )

    archs4 = _archs4_genes(archs4_path)
    available = set(archs4["mouse_ensembl_gene"])
    symbol_groups = {
        symbol: group.sort_values("mouse_ensembl_gene", kind="stable")
        for symbol, group in archs4.assign(
            _symbol_key=archs4["mouse_symbol"].str.casefold()
        ).groupby("_symbol_key", sort=False)
    }
    mapped = mapped.loc[mapped["mouse_ensembl_gene"].isin(available)].copy()
    mapped_groups = {
        gene: group for gene, group in mapped.groupby("human_ensembl_gene", sort=False)
    }

    selected: list[dict[str, Any]] = []
    used_mouse: set[str] = set()
    covered_human: set[str] = set()
    for row in human.itertuples(index=False):
        candidates = mapped_groups.get(row.human_ensembl_gene)
        chosen = None
        if candidates is not None:
            for candidate in candidates.itertuples(index=False):
                if candidate.mouse_ensembl_gene not in used_mouse:
                    chosen = {
                        "human_symbol": row.human_symbol,
                        "human_ensembl_gene": row.human_ensembl_gene,
                        "mouse_ensembl_gene": candidate.mouse_ensembl_gene,
                        "mouse_symbol": candidate.mouse_symbol,
                        "orthology_type": candidate.orthology_type,
                        "mouse_perc_id": candidate.mouse_perc_id,
                        "selection_role": "best_ensembl_ortholog",
                    }
                    break
        if chosen is None:
            exact = symbol_groups.get(str(row.human_symbol).casefold())
            if exact is not None:
                exact = exact.loc[~exact["mouse_ensembl_gene"].isin(used_mouse)]
                if not exact.empty:
                    candidate = exact.iloc[0]
                    chosen = {
                        "human_symbol": row.human_symbol,
                        "human_ensembl_gene": row.human_ensembl_gene,
                        "mouse_ensembl_gene": candidate["mouse_ensembl_gene"],
                        "mouse_symbol": candidate["mouse_symbol"],
                        "orthology_type": "case_insensitive_symbol_fallback",
                        "mouse_perc_id": np.nan,
                        "selection_role": "exact_symbol_fallback",
                    }
        if chosen is not None:
            selected.append(chosen)
            used_mouse.add(str(chosen["mouse_ensembl_gene"]))
            covered_human.add(str(row.human_ensembl_gene))

    extras = mapped.loc[~mapped["mouse_ensembl_gene"].isin(used_mouse)].sort_values(
        ["_orthology_rank", "_identity", "human_ensembl_gene", "mouse_ensembl_gene"],
        ascending=[True, False, True, True],
        kind="stable",
    )
    for candidate in extras.itertuples(index=False):
        if len(selected) >= dimensions:
            break
        if candidate.mouse_ensembl_gene in used_mouse:
            continue
        selected.append(
            {
                "human_symbol": candidate.human_symbol,
                "human_ensembl_gene": candidate.human_ensembl_gene,
                "mouse_ensembl_gene": candidate.mouse_ensembl_gene,
                "mouse_symbol": candidate.mouse_symbol,
                "orthology_type": candidate.orthology_type,
                "mouse_perc_id": candidate.mouse_perc_id,
                "selection_role": "additional_paralog_for_fixed_dimension",
            }
        )
        used_mouse.add(str(candidate.mouse_ensembl_gene))

    panel = pd.DataFrame(selected)
    if len(panel) != dimensions:
        raise ValueError(
            f"Could only construct {len(panel)} unique mouse landmarks; need {dimensions}"
        )
    if panel["mouse_ensembl_gene"].nunique() != dimensions:
        raise RuntimeError("Mouse landmark panel contains duplicate genes")
    panel.insert(0, "paper_dimension", np.arange(dimensions, dtype=int))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(output_path, sep="\t", index=False)

    missing = human.loc[~human["human_ensembl_gene"].isin(covered_human)]
    manifest = {
        "output": str(output_path),
        "dimensions": dimensions,
        "unique_mouse_genes": int(panel["mouse_ensembl_gene"].nunique()),
        "human_landmarks_covered": int(len(covered_human)),
        "human_landmarks_without_direct_dimension": int(len(missing)),
        "human_landmarks_without_direct_dimension_symbols": missing[
            "human_symbol"
        ].tolist(),
        "selection_roles": {
            str(key): int(value)
            for key, value in panel["selection_role"].value_counts().items()
        },
        "adaptation": (
            "The fixed 974-dimensional input contains one best mouse mapping per "
            "recoverable human landmark, then additional mapped paralogs for human "
            "landmarks with no direct mouse counterpart."
        ),
        "sources": {
            "human_landmarks": str(human_path),
            "human_landmarks_sha256": _sha256(human_path),
            "source_map": str(map_path),
            "source_map_sha256": _sha256(map_path),
            "archs4_h5": str(archs4_path),
        },
    }
    manifest_path = output_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return output_path, manifest
