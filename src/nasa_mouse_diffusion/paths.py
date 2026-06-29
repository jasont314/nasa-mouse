"""Input/output path helpers for diffusion bulk RNA-seq runs."""

from __future__ import annotations

from pathlib import Path

from nasa_mouse_wgan.paths import (
    DEFAULT_TISSUES,
    MUSCLE_GROUPS,
    direct_query_h5ad,
    first_existing,
    muscle_group_query_h5ad,
    muscle_group_reference_h5ad,
    reference_candidates,
)


OUTPUT_ROOT = Path("outputs/diffusion_conditional_generation")
SUMMARY_DIR = OUTPUT_ROOT / "summary"
PAPER_PDF = Path("docs/papers/lacan_diffusion_2026_reference.pdf")


def existing_pan_tissue_inputs(tissues: tuple[str, ...]):
    query_paths = []
    reference_paths = []
    for tissue in tissues:
        query = direct_query_h5ad(tissue)
        reference = first_existing(reference_candidates(tissue))
        if query.exists():
            query_paths.append(query)
        if reference is not None:
            reference_paths.append(reference)
    return query_paths, reference_paths
