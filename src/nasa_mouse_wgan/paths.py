"""Input/output path helpers for the WGAN bulk RNA-seq project."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


DEFAULT_TISSUES = (
    "liver",
    "skeletal_muscle",
    "skin",
    "kidney",
    "thymus",
    "spleen",
    "lung",
    "retina",
)

MUSCLE_GROUPS = (
    "soleus",
    "gastrocnemius",
    "quadriceps",
    "edl",
    "tibialis_anterior",
)


@dataclass(frozen=True)
class WGANRunSpec:
    """One direct or ARCHS4-pretrained WGAN run."""

    label: str
    tissue: str
    group: str
    mode: str
    query_h5ad: Path
    reference_h5ad: Path | None
    output_dir: Path

    @property
    def is_reference_run(self) -> bool:
        return self.reference_h5ad is not None


def direct_query_h5ad(tissue: str) -> Path:
    return (
        Path(f"outputs/expimap_direct_osdr_{tissue}")
        / "input"
        / f"osdr_{tissue}_flt_gc_reactome_raw_counts.h5ad"
    )


def reference_candidates(tissue: str) -> list[Path]:
    root = Path(f"outputs/expimap_archs4_reference_osdr_query_{tissue}")
    return [
        root / "reference_input_all" / f"archs4_mouse_{tissue}_reference_reactome_raw_counts.h5ad",
        root
        / "reference_input_5000_stratified"
        / f"archs4_mouse_{tissue}_reference_reactome_raw_counts.h5ad",
        root / "reference_input_1000" / f"archs4_mouse_{tissue}_reference_reactome_raw_counts.h5ad",
    ]


def first_existing(paths: list[Path]) -> Path | None:
    return next((path for path in paths if path.exists()), None)


def muscle_group_query_h5ad(group: str) -> Path:
    return (
        Path("outputs/expimap_muscle_targeted_combined_min8")
        / "group_inputs_exploratory_2acc"
        / f"osdr_skeletal_muscle_{group}_flt_gc_reactome_plus_muscle_raw_counts.h5ad"
    )


def muscle_group_reference_h5ad() -> Path:
    return (
        Path("outputs/expimap_muscle_targeted_combined_min8")
        / "reference_input_all"
        / "archs4_mouse_skeletal_muscle_reference_reactome_raw_counts.h5ad"
    )


def build_run_specs(
    *,
    tissues: tuple[str, ...] = DEFAULT_TISSUES,
    muscle_groups: tuple[str, ...] = MUSCLE_GROUPS,
    include_direct: bool = True,
    include_reference: bool = True,
    include_muscle_splits: bool = True,
) -> list[WGANRunSpec]:
    """Build runnable specs from existing API-derived input files."""

    specs: list[WGANRunSpec] = []
    for tissue in tissues:
        query = direct_query_h5ad(tissue)
        output_root = Path(f"outputs/wgan_{tissue}")
        if include_direct and query.exists():
            specs.append(
                WGANRunSpec(
                    label=f"{tissue}:direct",
                    tissue=tissue,
                    group="",
                    mode="direct_osdr",
                    query_h5ad=query,
                    reference_h5ad=None,
                    output_dir=output_root / "direct_osdr",
                )
            )
        reference = first_existing(reference_candidates(tissue))
        if include_reference and query.exists() and reference is not None:
            specs.append(
                WGANRunSpec(
                    label=f"{tissue}:reference",
                    tissue=tissue,
                    group="",
                    mode="archs4_pretrain_osdr_finetune",
                    query_h5ad=query,
                    reference_h5ad=reference,
                    output_dir=output_root / "archs4_pretrain_osdr_finetune",
                )
            )

    if include_muscle_splits:
        reference = muscle_group_reference_h5ad()
        for group in muscle_groups:
            query = muscle_group_query_h5ad(group)
            output_root = Path("outputs/wgan_skeletal_muscle_splits") / group
            if include_direct and query.exists():
                specs.append(
                    WGANRunSpec(
                        label=f"skeletal_muscle:{group}:direct",
                        tissue="skeletal_muscle",
                        group=group,
                        mode="direct_osdr",
                        query_h5ad=query,
                        reference_h5ad=None,
                        output_dir=output_root / "direct_osdr",
                    )
                )
            if include_reference and query.exists() and reference.exists():
                specs.append(
                    WGANRunSpec(
                        label=f"skeletal_muscle:{group}:reference",
                        tissue="skeletal_muscle",
                        group=group,
                        mode="archs4_pretrain_osdr_finetune",
                        query_h5ad=query,
                        reference_h5ad=reference,
                        output_dir=output_root / "archs4_pretrain_osdr_finetune",
                    )
                )
    return specs
