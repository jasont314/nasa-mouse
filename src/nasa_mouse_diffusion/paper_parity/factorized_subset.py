"""Create development-only tissue subsets of prepared factorized DDIM data."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd


def subset_factorized_data(
    source_h5: str | Path,
    samples_tsv: str | Path,
    output_h5: str | Path,
    *,
    tissues: Iterable[str],
) -> Path:
    selected_tissues = tuple(sorted(set(map(str, tissues))))
    if not selected_tissues:
        raise ValueError("At least one tissue is required")
    samples = pd.read_csv(samples_tsv, sep="\t")
    if "_row_index" not in samples or "tissue" not in samples:
        raise ValueError("Samples table lacks _row_index or tissue")
    tissue_by_row = samples.set_index("_row_index")["tissue"].astype(str)
    output = Path(output_h5)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    role_counts: dict[str, int] = {}
    with h5py.File(source_h5, "r") as source, h5py.File(temporary, "w") as target:
        for key, value in source.attrs.items():
            target.attrs[key] = value
        target.attrs["subset_tissues"] = json.dumps(selected_tissues)
        target.attrs["locked_test_copied"] = False
        for name, dataset in source.items():
            if isinstance(dataset, h5py.Dataset):
                source.copy(dataset, target, name=name)
        for role in ("train", "validation"):
            if role not in source:
                continue
            source_rows = np.asarray(source[role]["source_row"][:], dtype=np.int64)
            aligned_tissues = tissue_by_row.loc[source_rows].to_numpy(dtype=str)
            keep = np.isin(aligned_tissues, selected_tissues)
            group = target.create_group(role)
            for name, dataset in source[role].items():
                group.create_dataset(
                    name,
                    data=dataset[:][keep],
                    compression="gzip",
                    compression_opts=4,
                    shuffle=True,
                )
            role_counts[role] = int(keep.sum())
        target.attrs["subset_role_counts"] = json.dumps(role_counts, sort_keys=True)
    temporary.replace(output)
    manifest = {
        "source_h5": str(Path(source_h5).resolve()),
        "samples_tsv": str(Path(samples_tsv).resolve()),
        "output_h5": str(output.resolve()),
        "tissues": list(selected_tissues),
        "role_counts": role_counts,
        "locked_test_copied": False,
    }
    output.with_suffix(output.suffix + ".manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return output
