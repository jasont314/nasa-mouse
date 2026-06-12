"""Gene-intersection utilities for TMS pretraining and OSDR fine-tuning."""

from __future__ import annotations

import argparse
from pathlib import Path

from .io import load_matrix_bundle, require_import, write_matrix_bundle


def _row_slice(matrix, indices):
    scipy_sparse = require_import(
        "scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt"
    )
    if scipy_sparse.issparse(matrix):
        return matrix[indices, :]
    return matrix[indices, :]


def align_bundles(
    pretrain_manifest: str | Path,
    target_manifest: str | Path,
    output_prefix: str | Path,
    prefer: str = "pretrain",
) -> tuple[Path, Path]:
    """Write aligned pretrain and target bundles over the shared gene set."""
    pretrain = load_matrix_bundle(pretrain_manifest)
    target = load_matrix_bundle(target_manifest)

    pretrain_gene_to_idx = {gene: idx for idx, gene in enumerate(pretrain.genes)}
    target_gene_to_idx = {gene: idx for idx, gene in enumerate(target.genes)}
    shared = set(pretrain_gene_to_idx) & set(target_gene_to_idx)
    if not shared:
        raise SystemExit("No overlapping genes between pretraining and target bundles")

    if prefer == "target":
        ordered_genes = [gene for gene in target.genes if gene in shared]
    else:
        ordered_genes = [gene for gene in pretrain.genes if gene in shared]

    pre_idx = [pretrain_gene_to_idx[gene] for gene in ordered_genes]
    target_idx = [target_gene_to_idx[gene] for gene in ordered_genes]

    output_prefix = Path(output_prefix)
    pre_manifest = write_matrix_bundle(
        f"{output_prefix}.pretrain",
        _row_slice(pretrain.matrix, pre_idx),
        genes=ordered_genes,
        profiles=pretrain.profiles,
        profile_metadata=pretrain.profile_metadata,
        description=f"Pretraining matrix aligned to {len(ordered_genes)} shared genes",
    )
    target_manifest_out = write_matrix_bundle(
        f"{output_prefix}.target",
        _row_slice(target.matrix, target_idx),
        genes=ordered_genes,
        profiles=target.profiles,
        profile_metadata=target.profile_metadata,
        description=f"Target matrix aligned to {len(ordered_genes)} shared genes",
    )
    return pre_manifest, target_manifest_out


def main() -> None:
    parser = argparse.ArgumentParser(description="Align pretraining and OSDR genes.")
    parser.add_argument("--pretrain", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--prefer", choices=["pretrain", "target"], default="pretrain")
    args = parser.parse_args()

    pre, target = align_bundles(args.pretrain, args.target, args.output_prefix, args.prefer)
    print(pre)
    print(target)


if __name__ == "__main__":
    main()
