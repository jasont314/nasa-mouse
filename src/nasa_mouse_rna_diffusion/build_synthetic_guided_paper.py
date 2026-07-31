"""Build the synthetic-guided spaceflight manuscript package.

The builder intentionally consumes frozen analysis outputs. It does not train a
model, resample a cohort, or rerun feature selection.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PAPER_DIR = ROOT / "paper" / "synthetic_guided_spaceflight"
FIGURE_DIR = PAPER_DIR / "figures"
SOURCE_DIR = PAPER_DIR / "source_data"
UTILITY_TABLES_BEGIN = "<!-- BEGIN GENERATED TISSUE UTILITY TABLES -->"
UTILITY_TABLES_END = "<!-- END GENERATED TISSUE UTILITY TABLES -->"
LANDMARK_PANEL = ROOT / "data/diffusion/l974_mouse_paper_parity.tsv"

ARCHS4_RUN = (
    ROOT
    / "outputs/generative_benchmark/runs/lacan_diffusion/"
    "archs4_mouse_paper_parity_osdr_disjoint_seed1234"
)
OSDR_RUN = (
    ROOT
    / "outputs/generative_benchmark/runs/lacan_diffusion/"
    "osdr_factorized_study_lora512_correlation_refine_osdr_disjoint_seed2020"
)
LOCKED_DIR = OSDR_RUN / "evaluation/final_locked_test"
CONFIRM_DIR = (
    ROOT
    / "outputs/generative_benchmark/analyses/"
    "generated_feature_guidance_confirmation_disjoint_v1"
)
TRANSFER_DIR = (
    ROOT
    / "outputs/generative_benchmark/analyses/"
    "generated_feature_guidance_transfer_v1"
)
ADAPTIVE_HOLDOUT_DIR = (
    ROOT
    / "outputs/generative_benchmark/analyses/"
    "adaptive_per_tissue_ddim_augmentation_v1"
)
FRESH_HOLDOUT_DIR = (
    ROOT
    / "outputs/generative_benchmark/analyses/"
    "fresh_holdout_contrastive_ddim_augmentation_v1"
)
MUSCLE_DIR = (
    ROOT
    / "outputs/generative_benchmark/analyses/"
    "within_study_generated_feature_stability_muscle_groups_osdr_disjoint_v1"
)
TISSUE_DIR = (
    ROOT
    / "outputs/generative_benchmark/analyses/"
    "within_study_generated_feature_stability_osdr_disjoint_v1"
)
WGAN_DIR = (
    ROOT
    / "outputs/generative_benchmark/runs/vinas_wgan_gp/"
    "osdr_matched_study_conditioned_seed2020/evaluation/matched_validation"
)
HARMONIZATION_DIR = (
    ROOT / "outputs/generative_benchmark/summary/liver_harmonization"
)
GENEJEPA_DIR = (
    ROOT
    / "outputs/generative_benchmark/runs/genejepa/"
    "matrix_phase_0_genejepa_exact_mouse_one_epoch_f2e01cf1f130d5cb"
)

COLORS = {
    "navy": "#23445D",
    "teal": "#17807E",
    "blue": "#3E78A8",
    "coral": "#D96552",
    "gold": "#D69A2D",
    "green": "#4F845C",
    "purple": "#8064A2",
    "gray": "#7B858C",
    "light": "#EDF2F3",
    "dark": "#23313A",
}


def _required(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Required frozen input is missing: {path}")
    return path


def _read_json(path: Path) -> dict[str, Any]:
    with _required(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(_required(path), sep="\t")


def _write_tsv(frame: pd.DataFrame, name: str) -> Path:
    path = SOURCE_DIR / name
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def _save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def _clean_term(term: str) -> str:
    text = re.sub(r"^R-MMU-\d+_", "", term).replace("_", " ").lower()
    replacements = {
        "apc c": "APC/C",
        "cdc20": "CDC20",
        "g2 m": "G2/M",
        "dna": "DNA",
        "rna": "RNA",
        "tp53": "TP53",
        "ub ": "UB ",
    }
    for source, target in replacements.items():
        text = text.replace(source, target)
    return text[:1].upper() + text[1:]


def _assert_close(value: float, expected: float, label: str, tolerance: float = 1e-5) -> None:
    if not math.isclose(float(value), expected, abs_tol=tolerance, rel_tol=tolerance):
        raise ValueError(f"{label} changed: observed {value}, expected {expected}")


def _gate_synthetic_selection(
    stable_genes: pd.DataFrame,
    arm_choices: pd.DataFrame,
) -> pd.DataFrame:
    """Remove synthetic attribution where the generated arm failed its metric gate."""
    eligible_tissues = set(
        arm_choices.loc[
            arm_choices["generated_arm_eligible_all_metrics"].astype(bool),
            "tissue",
        ].astype(str)
    )
    gated = stable_genes.copy()
    ineligible = ~gated["tissue"].astype(str).isin(eligible_tissues)
    gated.loc[ineligible, "stable_generated"] = False
    gated.loc[ineligible, "generated_selection_frequency"] = 0.0
    gated.loc[ineligible, "real_effect_supports_generated"] = False
    gated.loc[
        ineligible
        & gated["gene_set"].isin(["core_intersection", "generated_supported"]),
        "gene_set",
    ] = "exploratory_union"
    return gated


def _ordinary_fdr_gene_inventory(
    stable_genes: pd.DataFrame,
    random_effects: pd.DataFrame,
    analysis_scope: str,
) -> pd.DataFrame:
    effect_supported = (
        stable_genes["real_effect_supports_generated"]
        .astype(str)
        .str.lower()
        .eq("true")
    )
    retained = stable_genes.loc[
        stable_genes["gene_set"].isin(["core_intersection", "generated_supported"])
        & stable_genes["real_meta_fdr"].lt(0.05)
        & effect_supported
    ].copy()
    retained["selection_interpretation"] = retained["gene_set"].map(
        {
            "core_intersection": "reinforced_real_and_synthetic",
            "generated_supported": "synthetic_promoted",
        }
    )
    retained["flt_gc_direction"] = np.where(
        retained["real_meta_effect"].gt(0),
        "FLT_higher",
        "FLT_lower",
    )
    retained["all_accessions_same_direction"] = np.isclose(
        retained["real_accession_direction_fraction"],
        1.0,
    )
    retained["loo_fdr_stable_0_05"] = (
        retained["real_loo_fdr_stable_0_05"]
        .astype(str)
        .str.lower()
        .eq("true")
    )

    random_effect_columns = random_effects[
        [
            "tissue",
            "gene",
            "n_accessions",
            "minimum_leave_one_out_fdr",
            "maximum_leave_one_out_fdr",
        ]
    ].copy()
    retained = retained.merge(
        random_effect_columns,
        on=["tissue", "gene"],
        how="left",
        validate="one_to_one",
    )
    retained.insert(0, "analysis_scope", analysis_scope)
    return retained[
        [
            "analysis_scope",
            "tissue",
            "gene",
            "symbol",
            "selection_interpretation",
            "n_accessions",
            "flt_gc_direction",
            "real_selection_frequency",
            "generated_selection_frequency",
            "real_meta_effect",
            "real_meta_fdr",
            "real_accession_direction_fraction",
            "all_accessions_same_direction",
            "loo_fdr_stable_0_05",
            "minimum_leave_one_out_fdr",
            "maximum_leave_one_out_fdr",
        ]
    ]


def _all_bh_fdr_gene_inventory(
    stable_genes: pd.DataFrame,
    random_effects: pd.DataFrame,
    landmark_panel: pd.DataFrame,
    analysis_scope: str,
) -> pd.DataFrame:
    family_sizes = random_effects.groupby("tissue", observed=True).size()
    invalid_families = family_sizes.loc[family_sizes.ne(974)]
    if not invalid_families.empty:
        raise ValueError(
            "Expected BH correction over 974 genes within every tissue; observed "
            f"{invalid_families.to_dict()}"
        )

    symbols = landmark_panel[
        ["mouse_ensembl_gene", "mouse_symbol"]
    ].drop_duplicates("mouse_ensembl_gene")
    if len(symbols) != 974:
        raise ValueError(
            f"Expected 974 unique mouse landmark symbols, found {len(symbols)}"
        )

    selection_columns = stable_genes[
        [
            "tissue",
            "gene",
            "gene_set",
            "stable_real",
            "stable_generated",
            "real_selection_frequency",
            "generated_selection_frequency",
            "real_effect_supports_generated",
        ]
    ].copy()
    retained = random_effects.loc[random_effects["meta_fdr"].lt(0.05)].copy()
    retained = retained.merge(
        symbols,
        left_on="gene",
        right_on="mouse_ensembl_gene",
        how="left",
        validate="many_to_one",
    ).drop(columns="mouse_ensembl_gene")
    retained = retained.rename(columns={"mouse_symbol": "symbol"})
    if retained["symbol"].isna().any():
        missing = retained.loc[retained["symbol"].isna(), "gene"].unique().tolist()
        raise ValueError(f"Missing landmark symbols for BH-FDR genes: {missing[:5]}")

    retained = retained.merge(
        selection_columns,
        on=["tissue", "gene"],
        how="left",
        validate="one_to_one",
    )
    retained["gene_set"] = retained["gene_set"].fillna(
        "not_in_stable_selection_union"
    )
    for column in ["stable_real", "stable_generated"]:
        retained[column] = retained[column].fillna(False).astype(bool)
    for column in ["real_selection_frequency", "generated_selection_frequency"]:
        retained[column] = retained[column].fillna(0.0)
    retained["real_effect_supports_generated"] = (
        retained["real_effect_supports_generated"].fillna(False).astype(bool)
    )

    retained["selection_interpretation"] = "not_stably_selected"
    retained.loc[
        retained["stable_real"] & ~retained["stable_generated"],
        "selection_interpretation",
    ] = "real_only_selected"
    retained.loc[
        ~retained["stable_real"]
        & retained["stable_generated"]
        & ~retained["real_effect_supports_generated"],
        "selection_interpretation",
    ] = "synthetic_selected_without_real_direction_support"
    retained.loc[
        retained["stable_real"]
        & retained["stable_generated"]
        & ~retained["real_effect_supports_generated"],
        "selection_interpretation",
    ] = "selected_both_direction_discordant"
    retained.loc[
        retained["gene_set"].eq("generated_supported"),
        "selection_interpretation",
    ] = "synthetic_promoted"
    retained.loc[
        retained["gene_set"].eq("core_intersection"),
        "selection_interpretation",
    ] = "reinforced_real_and_synthetic"

    retained["flt_gc_direction"] = np.where(
        retained["meta_effect"].gt(0),
        "FLT_higher",
        "FLT_lower",
    )
    retained["all_accessions_same_direction"] = np.isclose(
        retained["accession_direction_fraction"],
        1.0,
    )
    retained["directionally_heterogeneous"] = (
        ~retained["all_accessions_same_direction"]
    )
    retained["direction_consistency"] = np.where(
        retained["all_accessions_same_direction"],
        "unanimous",
        "non_unanimous",
    )
    retained.insert(0, "analysis_scope", analysis_scope)
    retained.insert(4, "bh_family_size", 974)
    retained.insert(5, "bh_scope", "within_tissue_974_gene_panel")
    retained.insert(6, "bh_fdr_threshold", 0.05)
    return retained[
        [
            "analysis_scope",
            "tissue",
            "gene",
            "symbol",
            "bh_family_size",
            "bh_scope",
            "bh_fdr_threshold",
            "n_accessions",
            "flt_gc_direction",
            "meta_effect",
            "meta_se",
            "meta_p",
            "meta_fdr",
            "tau2",
            "i2",
            "n_accession_same_direction",
            "n_accession_opposite_direction",
            "accession_direction_fraction",
            "direction_consistency",
            "all_accessions_same_direction",
            "directionally_heterogeneous",
            "selection_interpretation",
            "gene_set",
            "stable_real",
            "stable_generated",
            "real_selection_frequency",
            "generated_selection_frequency",
            "real_effect_supports_generated",
            "n_leave_one_out",
            "n_same_direction",
            "minimum_leave_one_out_fdr",
            "maximum_leave_one_out_fdr",
            "loo_direction_stable",
            "loo_fdr_stable_0_05",
        ]
    ]


def _bh_fdr_tissue_summary(
    inventory: pd.DataFrame,
    tested_families: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for family in tested_families.itertuples(index=False):
        analysis_scope = str(family.analysis_scope)
        tissue = str(family.tissue)
        frame = inventory.loc[
            inventory["analysis_scope"].eq(analysis_scope)
            & inventory["tissue"].eq(tissue)
        ]
        interpretation = frame["selection_interpretation"]
        rows.append(
            {
                "analysis_scope": analysis_scope,
                "tissue": tissue,
                "bh_family_size": int(family.bh_family_size),
                "n_bh_fdr_genes": len(frame),
                "n_flt_higher": int(frame["flt_gc_direction"].eq("FLT_higher").sum()),
                "n_flt_lower": int(frame["flt_gc_direction"].eq("FLT_lower").sum()),
                "n_unanimous_direction": int(
                    frame["all_accessions_same_direction"].sum()
                ),
                "n_directionally_heterogeneous": int(
                    frame["directionally_heterogeneous"].sum()
                ),
                "n_reinforced_real_and_synthetic": int(
                    interpretation.eq("reinforced_real_and_synthetic").sum()
                ),
                "n_synthetic_promoted": int(
                    interpretation.eq("synthetic_promoted").sum()
                ),
                "n_real_only_selected": int(
                    interpretation.eq("real_only_selected").sum()
                ),
                "n_synthetic_selected_direction_discordant": int(
                    interpretation.isin(
                        [
                            "synthetic_selected_without_real_direction_support",
                            "selected_both_direction_discordant",
                        ]
                    ).sum()
                ),
                "n_not_stably_selected": int(
                    interpretation.eq("not_stably_selected").sum()
                ),
                "n_loo_fdr_stable_0_05": int(
                    frame["loo_fdr_stable_0_05"].astype(bool).sum()
                ),
                "minimum_bh_fdr": (
                    float(frame["meta_fdr"].min()) if len(frame) else float("nan")
                ),
                "median_i2": (
                    float(frame["i2"].median()) if len(frame) else float("nan")
                ),
            }
        )
    return pd.DataFrame(rows)


def build_source_tables() -> dict[str, pd.DataFrame]:
    arch_eval = _read_json(ARCHS4_RUN / "evaluation/summary.json")
    arch_run = _read_json(ARCHS4_RUN / "run_summary.json")
    locked = _read_tsv(LOCKED_DIR / "repeat_metrics.tsv")
    locked_run = _read_json(LOCKED_DIR / "summary.json")
    wgan_run = _read_json(WGAN_DIR / "summary.json")
    wgan_repeats = _read_tsv(WGAN_DIR / "calibrated_repeat_metrics.tsv")
    harmonization = _read_tsv(HARMONIZATION_DIR / "independent_metrics.tsv")
    confirmation = _read_tsv(CONFIRM_DIR / "tissue_results.tsv")
    transfer_screen = _read_tsv(TRANSFER_DIR / "tissue_results.tsv")
    adaptive_holdout = _read_json(ADAPTIVE_HOLDOUT_DIR / "final_summary.json")
    fresh_holdout = _read_json(FRESH_HOLDOUT_DIR / "final_summary.json")
    muscle_holdout_extension = _read_json(
        FRESH_HOLDOUT_DIR / "confirmatory_folds/aggregate_summary.json"
    )
    genotype = _read_tsv(CONFIRM_DIR / "genotype_subgroup_results.tsv")
    thymus_features = _read_tsv(CONFIRM_DIR / "thymus/feature_stability.tsv")
    thymus_reactome = _read_tsv(CONFIRM_DIR / "thymus/reactome_enrichment.tsv")
    muscle_choices = _read_tsv(MUSCLE_DIR / "tissue_arm_choices.tsv")
    muscle_repeats = _read_tsv(MUSCLE_DIR / "paired_repeat_support.tsv")
    muscle_genes = _gate_synthetic_selection(
        _read_tsv(MUSCLE_DIR / "stable_gene_sets.tsv.gz"),
        muscle_choices,
    )
    muscle_reactome = _read_tsv(MUSCLE_DIR / "reactome_enrichment.tsv.gz")
    muscle_inventory = _read_tsv(MUSCLE_DIR / "tissue_inventory.tsv")
    muscle_accession_effects = _read_tsv(MUSCLE_DIR / "real_accession_effects.tsv.gz")
    muscle_random_effects = _read_tsv(MUSCLE_DIR / "real_random_effects.tsv.gz")
    tissue_choices = _read_tsv(TISSUE_DIR / "tissue_arm_choices.tsv")
    tissue_repeats = _read_tsv(TISSUE_DIR / "paired_repeat_support.tsv")
    tissue_biology = _read_tsv(TISSUE_DIR / "biological_support_summary.tsv")
    tissue_inventory = _read_tsv(TISSUE_DIR / "tissue_inventory.tsv")
    tissue_genes = _gate_synthetic_selection(
        _read_tsv(TISSUE_DIR / "stable_gene_sets.tsv.gz"),
        tissue_choices,
    )
    tissue_accession_effects = _read_tsv(TISSUE_DIR / "real_accession_effects.tsv.gz")
    tissue_random_effects = _read_tsv(TISSUE_DIR / "real_random_effects.tsv.gz")
    landmark_panel = _read_tsv(LANDMARK_PANEL)

    _assert_close(
        arch_eval["synthetic_to_real_test_tissue_classifier"]["balanced_accuracy"],
        0.7810085910974481,
        "ARCHS4 reverse-validation balanced accuracy",
    )
    _assert_close(
        confirmation.loc[confirmation["tissue"] == "thymus", "generated_roc_auc"].iloc[0],
        0.9791666666666666,
        "held-out thymus guided AUROC",
    )
    _assert_close(
        wgan_run["calibrated"]["metric_repeat_stability"]
        ["adversarial_accuracy"]["mean"],
        0.6361940298507462,
        "calibrated WGAN adversarial accuracy",
    )
    if len(harmonization) != 9:
        raise ValueError(
            "Expected nine liver harmonization arms, "
            f"found {len(harmonization)}"
        )
    if harmonization["heldout_absolute_gate"].astype(bool).any():
        raise ValueError("No liver harmonization arm should pass the absolute gate")

    inventory = pd.DataFrame(
        [
            {
                "source": "ARCHS4 mouse v2.5",
                "available_profiles": 997_515,
                "analysis_profiles": sum(arch_run["profiles"].values()),
                "analysis_split": (
                    f"{arch_run['profiles']['train']} train / "
                    f"{arch_run['profiles']['validation']} validation / "
                    f"{arch_run['profiles']['test']} test"
                ),
                "classes_or_accessions": 20,
                "genes": 974,
                "role": "healthy-preferred tissue pretraining",
            },
            {
                "source": "NASA OSDR API",
                "available_profiles": 1_631,
                "analysis_profiles": 1_610,
                "analysis_split": "781 train / 536 validation / 293 locked test",
                "classes_or_accessions": 75,
                "genes": 974,
                "role": "conditional adaptation and FLT/GC analysis",
            },
        ]
    )

    pipeline_design = pd.DataFrame(
        [
            {
                "axis": "Expression input",
                "configurable_options": (
                    "raw counts; CPM; TPM; log1p/log2p1; gene z-score; "
                    "robust or MaxAbs scaling"
                ),
                "evaluated_scope": (
                    "shared preprocessing screens plus each model's native contract"
                ),
                "selected_branch": "full-transcriptome TPM, 974 landmarks, train MaxAbs",
            },
            {
                "axis": "Feature space",
                "configurable_options": (
                    "all shared genes; fold-selected HVGs; Reactome genes; "
                    "mapped mouse L1000 landmarks"
                ),
                "evaluated_scope": "gated screens; selection statistics fitted on training folds",
                "selected_branch": "974 mapped mouse L1000 landmarks",
            },
            {
                "axis": "Harmonization",
                "configurable_options": (
                    "none; within-study z-score; within-study then global z-score; "
                    "ComBat; ComBat-seq; three MBatch methods; MOBER"
                ),
                "evaluated_scope": "nine matched 15,000-epoch liver DDIM arms",
                "selected_branch": "no global correction; explicit study conditioning",
            },
            {
                "axis": "Training data",
                "configurable_options": (
                    "OSDR only; ARCHS4 only; ARCHS4 pretraining then OSDR adaptation"
                ),
                "evaluated_scope": "paper-native and practical gated phases",
                "selected_branch": "ARCHS4 pretraining then OSDR adaptation",
            },
            {
                "axis": "Cohort scope",
                "configurable_options": (
                    "single accession; selected accessions; all eligible accessions"
                ),
                "evaluated_scope": "accession-scope screen followed by all eligible OSDR",
                "selected_branch": "all eligible OSDR accessions",
            },
            {
                "axis": "Tissue structure",
                "configurable_options": "pooled tissue-conditioned; standalone per tissue",
                "evaluated_scope": (
                    "pooled generator plus tissue-specific downstream policies"
                ),
                "selected_branch": "pooled tissue-conditioned generator",
            },
            {
                "axis": "Conditioning",
                "configurable_options": (
                    "FLT/GC; tissue; material; muscle group; study; sex; assay; "
                    "platform; source"
                ),
                "evaluated_scope": (
                    "condition and study policies screened; available covariates audited"
                ),
                "selected_branch": "tissue, FLT/GC, accession, and material",
            },
            {
                "axis": "Generator",
                "configurable_options": "WGAN-GP; DDIM; GeneJEPA representation screen",
                "evaluated_scope": (
                    "paper-reproduced architectures, then staged independent gates"
                ),
                "selected_branch": "ARCHS4-pretrained, OSDR-adapted DDIM",
            },
            {
                "axis": "Validation",
                "configurable_options": (
                    "GEO-series or accession-grouped validation; locked test; "
                    "unconditional controls; multiple generation seeds"
                ),
                "evaluated_scope": "no sample-random model-selection split",
                "selected_branch": "four-seed 293-profile locked OSDR test",
            },
        ]
    )

    harmonization_summary = harmonization[
        [
            "label",
            "normalization",
            "harmonization",
            "transductive_preprocessing",
            "heldout_corr",
            "heldout_precision",
            "heldout_recall",
            "heldout_f1",
            "heldout_aa",
            "heldout_fd_ratio",
            "heldout_absolute_gate",
            "delta_correlation",
            "direction_agreement",
            "conditional_effect_gate",
            "accession_meta_correlation",
            "accession_meta_direction",
            "accession_effect_gate",
        ]
    ].copy()
    harmonization_summary = harmonization_summary.rename(
        columns={
            "label": "method",
            "heldout_corr": "correlation",
            "heldout_precision": "precision",
            "heldout_recall": "recall",
            "heldout_f1": "f1",
            "heldout_aa": "adversarial_accuracy",
            "heldout_fd_ratio": "frechet_ratio",
            "heldout_absolute_gate": "fidelity_gate",
            "delta_correlation": "condition_effect_correlation",
            "direction_agreement": "condition_direction_agreement",
            "accession_meta_correlation": "accession_effect_correlation",
            "accession_meta_direction": "accession_direction_agreement",
        }
    )
    harmonization_summary["all_required_gates"] = (
        harmonization_summary["fidelity_gate"].astype(bool)
        & harmonization_summary["conditional_effect_gate"].astype(bool)
        & harmonization_summary["accession_effect_gate"].astype(bool)
    )

    locked_summary = pd.DataFrame(
        [
            {
                "metric": "gene_correlation_agreement",
                "mean": locked["correlation"].mean(),
                "minimum": locked["correlation"].min(),
                "maximum": locked["correlation"].max(),
                "target": "finite-sample floor >= 0.950; paper target >= 0.980",
                "repeats_passing": int(locked["fidelity_pass"].sum()),
                "repeats": len(locked),
            },
            {
                "metric": "precision",
                "mean": locked["precision"].mean(),
                "minimum": locked["precision"].min(),
                "maximum": locked["precision"].max(),
                "target": ">= 0.950",
                "repeats_passing": int((locked["precision"] >= 0.95).sum()),
                "repeats": len(locked),
            },
            {
                "metric": "recall",
                "mean": locked["recall"].mean(),
                "minimum": locked["recall"].min(),
                "maximum": locked["recall"].max(),
                "target": ">= 0.850",
                "repeats_passing": int((locked["recall"] >= 0.85).sum()),
                "repeats": len(locked),
            },
            {
                "metric": "f1",
                "mean": locked["f1"].mean(),
                "minimum": locked["f1"].min(),
                "maximum": locked["f1"].max(),
                "target": ">= 0.900",
                "repeats_passing": int((locked["f1"] >= 0.90).sum()),
                "repeats": len(locked),
            },
            {
                "metric": "adversarial_accuracy",
                "mean": locked["adversarial_accuracy"].mean(),
                "minimum": locked["adversarial_accuracy"].min(),
                "maximum": locked["adversarial_accuracy"].max(),
                "target": "0.400 to 0.600",
                "repeats_passing": int(
                    locked["adversarial_accuracy"].between(0.40, 0.60).sum()
                ),
                "repeats": len(locked),
            },
            {
                "metric": "frechet_ratio_to_real_split_p95",
                "mean": locked["frechet_ratio"].mean(),
                "minimum": locked["frechet_ratio"].min(),
                "maximum": locked["frechet_ratio"].max(),
                "target": "<= 1.000",
                "repeats_passing": int((locked["frechet_ratio"] <= 1.0).sum()),
                "repeats": len(locked),
            },
            {
                "metric": "pooled_flt_gc_effect_correlation",
                "mean": locked["condition_delta_correlation"].mean(),
                "minimum": locked["condition_delta_correlation"].min(),
                "maximum": locked["condition_delta_correlation"].max(),
                "target": "r >= 0.300 and direction >= 0.550",
                "repeats_passing": int(locked["condition_effect_pass"].sum()),
                "repeats": len(locked),
            },
            {
                "metric": "muscle_accession_effect_correlation",
                "mean": locked["muscle_accession_correlation"].mean(),
                "minimum": locked["muscle_accession_correlation"].min(),
                "maximum": locked["muscle_accession_correlation"].max(),
                "target": "r >= 0.300 and direction >= 0.550",
                "repeats_passing": int(locked["muscle_accession_pass"].sum()),
                "repeats": len(locked),
            },
        ]
    )

    arch_summary = pd.DataFrame(
        [
            ("Real train to real test tissue BA", arch_eval["real_train_to_test_tissue_classifier"]["balanced_accuracy"]),
            ("Synthetic train to real test tissue BA", arch_eval["synthetic_to_real_test_tissue_classifier"]["balanced_accuracy"]),
            ("Gene mean correlation", arch_eval["gene_mean_correlation"]),
            ("Gene SD correlation", arch_eval["gene_standard_deviation_correlation"]),
            ("Gene correlation agreement", arch_eval["gene_correlation_matrix_agreement"]),
            ("Precision, scaled L974", arch_eval["precision_recall_in_scaled_l974"]["precision"]),
            ("Recall, scaled L974", arch_eval["precision_recall_in_scaled_l974"]["recall"]),
            ("Adversarial accuracy", arch_eval["nearest_neighbor_adversarial_accuracy_in_scaled_l974"]),
            ("Precision, PCA50", arch_eval["precision_recall_in_train_pca50"]["precision"]),
            ("Recall, PCA50", arch_eval["precision_recall_in_train_pca50"]["recall"]),
            ("Frechet distance, PCA50", arch_eval["frechet_distance_in_train_pca50"]),
        ],
        columns=["metric", "value"],
    )

    wgan_metrics = wgan_run["calibrated"]["metric_repeat_stability"]
    arch_precision = arch_eval["precision_recall_in_scaled_l974"]["precision"]
    arch_recall = arch_eval["precision_recall_in_scaled_l974"]["recall"]
    arch_f1 = 2 * arch_precision * arch_recall / (arch_precision + arch_recall)
    model_screen = pd.DataFrame(
        [
            {
                "model": "Broad-reference DDIM",
                "training_regime": "ARCHS4 only",
                "evaluation_split": "4,628 held-out ARCHS4 profiles; complete GEO series",
                "generation_repeats": 1,
                "correlation": arch_eval["gene_correlation_matrix_agreement"],
                "precision": arch_precision,
                "recall": arch_recall,
                "f1": arch_f1,
                "adversarial_accuracy": arch_eval[
                    "nearest_neighbor_adversarial_accuracy_in_scaled_l974"
                ],
                "frechet_ratio": arch_eval["frechet_ratio_to_real_split_p95"],
                "fidelity_repeats_passing": "0/1",
                "condition_repeats_passing": "not applicable",
                "accession_repeats_passing": "not applicable",
                "locked_test_opened": "not applicable",
                "decision": (
                    "retained as tissue-conditioned initialization despite failed "
                    "correlation-structure gate"
                ),
            },
            {
                "model": "Study-conditioned WGAN-GP",
                "training_regime": "OSDR matched study-conditioned",
                "evaluation_split": "536-profile validation; test unopened",
                "generation_repeats": len(wgan_repeats),
                "correlation": wgan_metrics["correlation"]["mean"],
                "precision": wgan_metrics["precision"]["mean"],
                "recall": wgan_metrics["recall"]["mean"],
                "f1": wgan_metrics["f1"]["mean"],
                "adversarial_accuracy": wgan_metrics["adversarial_accuracy"]["mean"],
                "frechet_ratio": wgan_metrics["frechet_ratio"]["mean"],
                "fidelity_repeats_passing": (
                    f"{int(wgan_repeats['fidelity_pass'].sum())}/{len(wgan_repeats)}"
                ),
                "condition_repeats_passing": (
                    f"{int(wgan_repeats['condition_effect_pass'].sum())}/{len(wgan_repeats)}"
                ),
                "accession_repeats_passing": (
                    f"{int(wgan_repeats['muscle_accession_pass'].sum())}/{len(wgan_repeats)}"
                ),
                "locked_test_opened": bool(wgan_run["locked_test_opened"]),
                "decision": (
                    "rejected on validation: external separability and unstable "
                    "accession-aware effect recovery"
                ),
            },
            {
                "model": "Factorized DDIM",
                "training_regime": "ARCHS4 pretraining then OSDR adaptation",
                "evaluation_split": "293-profile locked within-study test",
                "generation_repeats": len(locked),
                "correlation": locked["correlation"].mean(),
                "precision": locked["precision"].mean(),
                "recall": locked["recall"].mean(),
                "f1": locked["f1"].mean(),
                "adversarial_accuracy": locked["adversarial_accuracy"].mean(),
                "frechet_ratio": locked["frechet_ratio"].mean(),
                "fidelity_repeats_passing": (
                    f"{int(locked['fidelity_pass'].sum())}/{len(locked)}"
                ),
                "condition_repeats_passing": (
                    f"{int(locked['condition_effect_pass'].sum())}/{len(locked)}"
                ),
                "accession_repeats_passing": (
                    f"{int(locked['muscle_accession_pass'].sum())}/{len(locked)}"
                ),
                "locked_test_opened": bool(locked_run["locked_test_opened"]),
                "decision": (
                    "selected: only generator to pass the final joint locked gates"
                ),
            },
            {
                "model": "GeneJEPA",
                "training_regime": "ARCHS4 representation screen",
                "evaluation_split": "held-out ARCHS4 series",
                "generation_repeats": 0,
                "correlation": np.nan,
                "precision": np.nan,
                "recall": np.nan,
                "f1": np.nan,
                "adversarial_accuracy": np.nan,
                "frechet_ratio": np.nan,
                "fidelity_repeats_passing": "not applicable",
                "condition_repeats_passing": "not applicable",
                "accession_repeats_passing": "not applicable",
                "locked_test_opened": False,
                "decision": "not a generator: released architecture has no expression decoder",
            },
        ]
    )

    utility = locked_run["classifier_utility"]
    naive_utility = pd.DataFrame(
        [
            (
                "Real OSDR",
                utility["real_train_real_evaluation"]["balanced_accuracy"],
                utility["real_train_real_evaluation"]["roc_auc"],
            ),
            (
                "Synthetic",
                utility["synthetic_train_real_evaluation"]["balanced_accuracy"],
                utility["synthetic_train_real_evaluation"]["roc_auc"],
            ),
            (
                "Real + synthetic",
                utility["real_plus_synthetic_train_real_evaluation"][
                    "balanced_accuracy"
                ],
                utility["real_plus_synthetic_train_real_evaluation"]["roc_auc"],
            ),
        ],
        columns=["training_data", "balanced_accuracy", "roc_auc"],
    )

    core_symbols = ["Birc5", "Ccne2", "Gmnn", "Ube2c", "Cdk1", "Nusap1", "Ccnb1", "Ccnb2"]
    thymus_core = (
        thymus_features.loc[thymus_features["symbol"].isin(core_symbols)]
        .copy()
        .sort_values("mean_real_effect")
    )
    if set(thymus_core["symbol"]) != set(core_symbols):
        missing = sorted(set(core_symbols) - set(thymus_core["symbol"]))
        raise ValueError(f"Missing held-out thymus core genes: {missing}")

    soleus_genes = muscle_genes.loc[
        (muscle_genes["tissue"] == "soleus")
        & muscle_genes["gene_set"].eq("core_intersection")
        & muscle_genes["real_meta_fdr"].lt(0.05)
        & muscle_genes["real_effect_supports_generated"].astype(bool)
    ].copy()
    soleus_genes = soleus_genes.sort_values("real_meta_effect")
    if len(soleus_genes) != 5:
        raise ValueError(
            "Expected five reinforced soleus genes, "
            f"found {len(soleus_genes)}"
        )

    rbm6_gene = "ENSMUSG00000032582"
    quadriceps_rbm6_effects = muscle_accession_effects.loc[
        muscle_accession_effects["analysis_tissue"].eq("quadriceps")
        & muscle_accession_effects["feature"].eq(rbm6_gene)
    ].copy()
    quadriceps_rbm6_effects = quadriceps_rbm6_effects.sort_values("accession")
    if len(quadriceps_rbm6_effects) != 4:
        raise ValueError(
            "Expected four quadriceps Rbm6 accession effects, "
            f"found {len(quadriceps_rbm6_effects)}"
        )
    if not quadriceps_rbm6_effects["flight_minus_ground"].gt(0).all():
        raise ValueError("Expected every quadriceps Rbm6 accession effect to be positive")

    quadriceps_rbm6_meta = muscle_random_effects.loc[
        muscle_random_effects["tissue"].eq("quadriceps")
        & muscle_random_effects["gene"].eq(rbm6_gene)
    ].copy()
    if len(quadriceps_rbm6_meta) != 1:
        raise ValueError(
            "Expected one quadriceps Rbm6 random-effects row, "
            f"found {len(quadriceps_rbm6_meta)}"
        )
    _assert_close(
        quadriceps_rbm6_meta["meta_fdr"].iloc[0],
        0.0007503147625256085,
        "quadriceps Rbm6 random-effects FDR",
    )

    supported_gene_counts = (
        muscle_genes.loc[
            muscle_genes["real_loo_fdr_stable_0_05"].astype(bool)
            & muscle_genes["real_effect_supports_generated"].astype(bool)
        ]
        .groupby("tissue")
        .size()
        .rename("cross_study_supported_genes")
        .reset_index()
    )
    muscle_summary = muscle_choices.merge(
        muscle_repeats,
        on=["tissue", "selected_arm"],
        how="left",
        validate="one_to_one",
    ).merge(
        muscle_inventory,
        on="tissue",
        how="left",
        validate="one_to_one",
    ).merge(
        supported_gene_counts,
        on="tissue",
        how="left",
        validate="one_to_one",
    )
    muscle_summary["cross_study_supported_genes"] = (
        muscle_summary["cross_study_supported_genes"].fillna(0).astype(int)
    )

    ordinary_fdr_genes = pd.concat(
        [
            _ordinary_fdr_gene_inventory(
                tissue_genes,
                tissue_random_effects,
                "canonical_tissue",
            ),
            _ordinary_fdr_gene_inventory(
                muscle_genes,
                muscle_random_effects,
                "skeletal_muscle_group",
            ),
        ],
        ignore_index=True,
    ).sort_values(
        ["analysis_scope", "tissue", "real_meta_fdr", "symbol"],
        ignore_index=True,
    )

    all_bh_fdr_genes = pd.concat(
        [
            _all_bh_fdr_gene_inventory(
                tissue_genes,
                tissue_random_effects,
                landmark_panel,
                "canonical_tissue",
            ),
            _all_bh_fdr_gene_inventory(
                muscle_genes,
                muscle_random_effects,
                landmark_panel,
                "skeletal_muscle_group",
            ),
        ],
        ignore_index=True,
    ).sort_values(
        ["analysis_scope", "tissue", "meta_fdr", "symbol"],
        ignore_index=True,
    )
    tested_bh_families = pd.concat(
        [
            tissue_random_effects.groupby("tissue", observed=True)
            .size()
            .rename("bh_family_size")
            .reset_index()
            .assign(analysis_scope="canonical_tissue"),
            muscle_random_effects.groupby("tissue", observed=True)
            .size()
            .rename("bh_family_size")
            .reset_index()
            .assign(analysis_scope="skeletal_muscle_group"),
        ],
        ignore_index=True,
    )[["analysis_scope", "tissue", "bh_family_size"]].sort_values(
        ["analysis_scope", "tissue"],
        ignore_index=True,
    )
    bh_fdr_tissue_summary = _bh_fdr_tissue_summary(
        all_bh_fdr_genes,
        tested_bh_families,
    )
    scope_counts = all_bh_fdr_genes.groupby("analysis_scope").size().to_dict()
    if scope_counts != {
        "canonical_tissue": 202,
        "skeletal_muscle_group": 257,
    }:
        raise ValueError(f"Unexpected BH-FDR inventory counts: {scope_counts}")
    selection_counts = all_bh_fdr_genes["selection_interpretation"].value_counts()
    if (
        int(selection_counts.get("synthetic_promoted", 0)) != 26
        or int(selection_counts.get("reinforced_real_and_synthetic", 0)) != 23
    ):
        raise ValueError(
            "Unexpected synthetic-informed BH-FDR counts: "
            f"{selection_counts.to_dict()}"
        )
    if len(bh_fdr_tissue_summary) != 27:
        raise ValueError(
            "Expected BH-FDR summaries for 22 canonical tissues and five "
            f"muscle groups, found {len(bh_fdr_tissue_summary)}"
        )

    thymus_development = all_bh_fdr_genes.loc[
        all_bh_fdr_genes["analysis_scope"].eq("canonical_tissue")
        & all_bh_fdr_genes["tissue"].eq("thymus")
        & all_bh_fdr_genes["symbol"].isin(core_symbols),
        [
            "gene",
            "symbol",
            "selection_interpretation",
            "stable_real",
            "stable_generated",
            "real_selection_frequency",
            "generated_selection_frequency",
            "meta_effect",
            "meta_fdr",
            "accession_direction_fraction",
            "all_accessions_same_direction",
        ],
    ].copy()
    thymus_evidence_mapping = thymus_core[
        [
            "gene",
            "symbol",
            "mean_classifier_coefficient",
            "mean_real_effect",
            "mean_synthetic_effect",
        ]
    ].merge(
        thymus_development,
        on=["gene", "symbol"],
        how="left",
        validate="one_to_one",
    )
    if len(thymus_evidence_mapping) != 8 or thymus_evidence_mapping[
        "selection_interpretation"
    ].isna().any():
        raise ValueError(
            "Expected all eight held-out thymus core genes to map to the "
            "cross-study BH-FDR inventory"
        )
    thymus_evidence_mapping.insert(0, "evidence_tier", 1)
    thymus_evidence_mapping.insert(1, "heldout_accession", "OSD-457")
    thymus_evidence_mapping.insert(
        4,
        "heldout_direction_in_wt_and_nrf2ko",
        "FLT_lower",
    )
    thymus_evidence_mapping = thymus_evidence_mapping.rename(
        columns={
            "mean_classifier_coefficient": "heldout_classifier_coefficient",
            "mean_real_effect": "heldout_mean_real_effect",
            "mean_synthetic_effect": "development_mean_synthetic_effect",
            "selection_interpretation": "tier_2_selection_interpretation",
            "stable_real": "tier_2_stable_real",
            "stable_generated": "tier_2_stable_generated",
            "real_selection_frequency": "tier_2_real_selection_frequency",
            "generated_selection_frequency": (
                "tier_2_generated_selection_frequency"
            ),
            "meta_effect": "real_cross_study_meta_effect",
            "meta_fdr": "real_cross_study_bh_fdr",
            "accession_direction_fraction": (
                "real_cross_study_direction_fraction"
            ),
            "all_accessions_same_direction": (
                "real_cross_study_unanimous_direction"
            ),
        }
    )

    tissue_summary = tissue_choices.merge(
        tissue_repeats,
        on=["tissue", "selected_arm"],
        how="left",
        validate="one_to_one",
    ).merge(
        tissue_biology,
        on=["tissue", "selected_arm"],
        how="left",
        validate="one_to_one",
    ).merge(
        tissue_inventory,
        on="tissue",
        how="left",
        validate="one_to_one",
    )

    development_tissues = [
        "thymus",
        "skeletal_muscle",
        "kidney",
        "spleen",
        "skin",
        "lung",
        "retina",
        "adrenal_gland",
    ]
    development_highlights = tissue_summary.loc[
        tissue_summary["tissue"].isin(development_tissues),
        [
            "tissue",
            "selected_arm",
            "real_mean_balanced_accuracy",
            "selected_mean_balanced_accuracy",
            "real_mean_roc_auc",
            "selected_mean_roc_auc",
            "real_mean_average_precision",
            "selected_mean_average_precision",
            "generated_arm_eligible_all_metrics",
        ],
    ].copy()
    development_highlights.insert(0, "analysis_scope", "canonical_tissue")
    muscle_highlights = muscle_summary.loc[
        muscle_summary["tissue"].isin(["soleus", "gastrocnemius"]),
        development_highlights.columns.drop("analysis_scope").tolist(),
    ].copy()
    muscle_highlights.insert(0, "analysis_scope", "skeletal_muscle_group")
    development_highlights = pd.concat(
        [development_highlights, muscle_highlights],
        ignore_index=True,
    )
    development_order = development_tissues + ["soleus", "gastrocnemius"]
    development_highlights["_order"] = development_highlights["tissue"].map(
        {tissue: index for index, tissue in enumerate(development_order)}
    )
    development_highlights = (
        development_highlights.sort_values("_order")
        .drop(columns="_order")
        .reset_index(drop=True)
    )
    if len(development_highlights) != 10:
        raise ValueError(
            "Expected ten tissue-specific utility highlights, "
            f"found {len(development_highlights)}"
        )
    if not development_highlights["generated_arm_eligible_all_metrics"].all():
        raise ValueError(
            "Every tissue-specific utility highlight must pass the arm gate"
        )

    transfer_rows: list[dict[str, Any]] = []
    for row in transfer_screen.itertuples(index=False):
        transfer_rows.append(
            {
                "experiment": "cross_tissue_feature_guidance_screen",
                "tissue": row.tissue,
                "test_accessions": int(row.test_accessions),
                "test_profiles": int(row.test_profiles),
                "synthetic_use": "feature_guidance",
                "real_balanced_accuracy": row.baseline_balanced_accuracy,
                "synthetic_informed_balanced_accuracy": (
                    row.generated_balanced_accuracy
                ),
                "real_roc_auc": row.baseline_roc_auc,
                "synthetic_informed_roc_auc": row.generated_roc_auc,
                "real_average_precision": row.baseline_average_precision,
                "synthetic_informed_average_precision": (
                    row.generated_average_precision
                ),
                "passed_declared_rule": bool(row.improved_ba_without_auc_loss),
                "interpretation": (
                    "advanced to fixed confirmation"
                    if row.improved_ba_without_auc_loss
                    else "did not advance"
                ),
            }
        )
    for row in confirmation.itertuples(index=False):
        passed = bool(row.improved_ba_without_auc_loss)
        transfer_rows.append(
            {
                "experiment": "fixed_lung_thymus_confirmation",
                "tissue": row.tissue,
                "test_accessions": int(row.test_accessions),
                "test_profiles": int(row.test_profiles),
                "synthetic_use": (
                    "guided_real_only"
                    if row.tissue == "thymus"
                    else "guided_low_weight"
                ),
                "real_balanced_accuracy": row.baseline_balanced_accuracy,
                "synthetic_informed_balanced_accuracy": (
                    row.generated_balanced_accuracy
                ),
                "real_roc_auc": row.baseline_roc_auc,
                "synthetic_informed_roc_auc": row.generated_roc_auc,
                "real_average_precision": row.baseline_average_precision,
                "synthetic_informed_average_precision": (
                    row.generated_average_precision
                ),
                "passed_declared_rule": passed,
                "interpretation": (
                    "retained study-held-out result"
                    if passed
                    else "failed fixed confirmation"
                ),
            }
        )
    for tissue, result in adaptive_holdout["tissues"].items():
        real = result["selected_real_only"]
        augmented = result["real_plus_synthetic"]
        transfer_rows.append(
            {
                "experiment": "adaptive_augmentation_screen",
                "tissue": tissue,
                "test_accessions": result["test_accessions"],
                "test_profiles": result["test_profiles"],
                "synthetic_use": "real_plus_generated",
                "real_balanced_accuracy": real["balanced_accuracy"],
                "synthetic_informed_balanced_accuracy": augmented[
                    "balanced_accuracy"
                ],
                "real_roc_auc": real["roc_auc"],
                "synthetic_informed_roc_auc": augmented["roc_auc"],
                "real_average_precision": real["average_precision"],
                "synthetic_informed_average_precision": augmented[
                    "average_precision"
                ],
                "passed_declared_rule": bool(result["success"]),
                "interpretation": (
                    "exploratory five-profile gain"
                    if tissue == "heart" and result["success"]
                    else "no improvement"
                ),
            }
        )
    for tissue, result in fresh_holdout["tissues"].items():
        real = result["selected_real_only"]
        augmented = result["real_plus_synthetic"]
        transfer_rows.append(
            {
                "experiment": "frozen_augmentation_initial_test",
                "tissue": tissue,
                "test_accessions": result["test_accessions"],
                "test_profiles": result["test_profiles"],
                "synthetic_use": "real_plus_generated",
                "real_balanced_accuracy": real["balanced_accuracy"],
                "synthetic_informed_balanced_accuracy": augmented[
                    "balanced_accuracy"
                ],
                "real_roc_auc": real["roc_auc"],
                "synthetic_informed_roc_auc": augmented["roc_auc"],
                "real_average_precision": real["average_precision"],
                "synthetic_informed_average_precision": augmented[
                    "average_precision"
                ],
                "passed_declared_rule": bool(result["success"]),
                "interpretation": (
                    "initial gain; expanded below"
                    if result["success"]
                    else "failed initial test"
                ),
            }
        )
    extension = muscle_holdout_extension["accession_macro"]
    transfer_rows.append(
        {
            "experiment": "frozen_augmentation_full_extension",
            "tissue": muscle_holdout_extension["tissue"],
            "test_accessions": muscle_holdout_extension["test_accessions"],
            "test_profiles": muscle_holdout_extension["test_profiles"],
            "synthetic_use": "real_plus_generated",
            "real_balanced_accuracy": extension["real_only"]["balanced_accuracy"],
            "synthetic_informed_balanced_accuracy": extension[
                "real_plus_synthetic"
            ]["balanced_accuracy"],
            "real_roc_auc": extension["real_only"]["roc_auc"],
            "synthetic_informed_roc_auc": extension["real_plus_synthetic"][
                "roc_auc"
            ],
            "real_average_precision": extension["real_only"][
                "average_precision"
            ],
            "synthetic_informed_average_precision": extension[
                "real_plus_synthetic"
            ]["average_precision"],
            "passed_declared_rule": False,
            "interpretation": "initial gain did not generalize across all accessions",
        }
    )
    study_holdout_context = pd.DataFrame(transfer_rows)
    _assert_close(
        study_holdout_context.loc[
            study_holdout_context["experiment"].eq(
                "frozen_augmentation_full_extension"
            ),
            "synthetic_informed_roc_auc",
        ].iloc[0],
        0.6897095959595959,
        "extended held-out muscle augmented AUROC",
    )
    if len(study_holdout_context) != 14:
        raise ValueError(
            "Expected fourteen whole-study transfer context rows, "
            f"found {len(study_holdout_context)}"
        )

    igfbp3_gene = "ENSMUSG00000020427"
    spleen_igfbp3_effects = tissue_accession_effects.loc[
        tissue_accession_effects["analysis_tissue"].eq("spleen")
        & tissue_accession_effects["feature"].eq(igfbp3_gene)
    ].copy()
    spleen_igfbp3_effects["standard_error"] = np.sqrt(
        spleen_igfbp3_effects["effect_variance"]
    )
    spleen_igfbp3_effects["ci_low"] = (
        spleen_igfbp3_effects["flight_minus_ground"]
        - 1.96 * spleen_igfbp3_effects["standard_error"]
    )
    spleen_igfbp3_effects["ci_high"] = (
        spleen_igfbp3_effects["flight_minus_ground"]
        + 1.96 * spleen_igfbp3_effects["standard_error"]
    )
    spleen_igfbp3_effects = spleen_igfbp3_effects[
        [
            "accession",
            "n_flight",
            "n_ground_control",
            "flight_minus_ground",
            "standard_error",
            "ci_low",
            "ci_high",
        ]
    ].sort_values("accession")
    if len(spleen_igfbp3_effects) != 6:
        raise ValueError(
            "Expected six spleen Igfbp3 accession effects, "
            f"found {len(spleen_igfbp3_effects)}"
        )
    if not spleen_igfbp3_effects["flight_minus_ground"].gt(0).all():
        raise ValueError("Expected every spleen Igfbp3 accession effect to be positive")

    spleen_igfbp3_meta = tissue_random_effects.loc[
        tissue_random_effects["tissue"].eq("spleen")
        & tissue_random_effects["gene"].eq(igfbp3_gene)
    ].copy()
    if len(spleen_igfbp3_meta) != 1:
        raise ValueError(
            "Expected one spleen Igfbp3 random-effects row, "
            f"found {len(spleen_igfbp3_meta)}"
        )
    _assert_close(
        spleen_igfbp3_meta["meta_fdr"].iloc[0],
        1.7592185680837808e-09,
        "spleen Igfbp3 random-effects FDR",
    )

    spleen_reference_expression = pd.DataFrame(
        [
            ("White-pulp mesenchymal", 1476.95),
            ("Red-pulp mesenchymal", 414.99),
            ("Endothelial", 8.01),
            ("Red-pulp macrophage", 2.37),
        ],
        columns=["population", "mean_igfbp3_rpkm"],
    )
    spleen_reference_expression.insert(0, "dataset", "GSE156162")

    evidence = pd.DataFrame(
        [
            {
                "tissue": "thymus",
                "tier": "held-out study validation",
                "tier_score": 3,
                "predictive_result": "BA 0.500 to 0.833; AUROC 0.840 to 0.979",
                "real_gene_support": "8 core FLT-down genes; genotype effect r=0.975",
                "pathway_support": "G2/M, APC/C, DNA replication; Reactome FDR < 0.05",
                "interpretation": (
                    "study-held-out feature-guidance result; prospective replication required"
                ),
            },
            {
                "tissue": "soleus",
                "tier": "cross-accession development",
                "tier_score": 2,
                "predictive_result": "real-plus-generated delta BA/AUROC/AP +0.038/+0.000/+0.006",
                "real_gene_support": "5 reinforced BH-FDR genes; 4 pass LOO FDR",
                "pathway_support": "mitochondrial lipid oxidation/protein turnover; FDR < 0.05",
                "interpretation": "coherent hypothesis; requires unseen-accession confirmation",
            },
            {
                "tissue": "kidney",
                "tier": "cross-accession development",
                "tier_score": 2,
                "predictive_result": "guided delta BA/AUROC/AP +0.053/+0.091/+0.115",
                "real_gene_support": "Inpp4b promoted and LOO-stable; Slc37a4 reinforced",
                "pathway_support": "no stable-set Reactome term at FDR < 0.05",
                "interpretation": "focused renal metabolic-signaling hypothesis",
            },
            {
                "tissue": "skeletal_muscle",
                "tier": "cross-accession development",
                "tier_score": 2,
                "predictive_result": "guided delta BA/AUROC/AP +0.071/+0.036/+0.037",
                "real_gene_support": "12 synthetic-informed BH-FDR genes; 9 pass LOO FDR",
                "pathway_support": "interferon signaling and sialic-acid metabolism; FDR < 0.05",
                "interpretation": "pooled signal complements anatomical soleus result",
            },
            {
                "tissue": "lung",
                "tier": "mixed held-out exploratory",
                "tier_score": 1,
                "predictive_result": (
                    "BA decreased; modest AUROC/AP gains; validation gate rejected"
                ),
                "real_gene_support": "no directionally confirmed core gene set",
                "pathway_support": "no Reactome term at FDR < 0.05",
                "interpretation": "no retained synthetic-guided classifier result",
            },
            {
                "tissue": "spleen",
                "tier": "developmental exploratory",
                "tier_score": 1,
                "predictive_result": "real-plus-generated delta BA/AUROC/AP +0.131/+0.163/+0.160",
                "real_gene_support": "Rai14, Ptprk, Myl9 promoted; Loxl1 reinforced; none pass LOO",
                "pathway_support": "no coherent stable-set Reactome enrichment",
                "interpretation": "adhesion/cytoskeletal hypothesis; Igfbp3 is real-data-only",
            },
            {
                "tissue": "skin",
                "tier": "developmental exploratory",
                "tier_score": 1,
                "predictive_result": "real-plus-generated delta BA/AUROC/AP +0.085/+0.076/+0.062",
                "real_gene_support": "Plscr1 is promoted and FLT-up in 6/6 studies; not LOO-stable",
                "pathway_support": "cell-cycle/DNA-repair theme matches published skin analyses",
                "interpretation": "literature-aligned developmental candidate",
            },
            {
                "tissue": "adrenal_gland",
                "tier": "developmental exploratory",
                "tier_score": 1,
                "predictive_result": "generated-only delta BA/AUROC/AP +0.141/+0.078/+0.070",
                "real_gene_support": "Psmb8 promoted and Tspan4 reinforced; FLT-lower in 3/3 studies",
                "pathway_support": "heat-shock/RNA-regulation enrichment; FDR < 0.05",
                "interpretation": "small-study developmental candidate; neither gene passes LOO",
            },
            {
                "tissue": "liver",
                "tier": "negative",
                "tier_score": 0,
                "predictive_result": "real-only arm retained by screen",
                "real_gene_support": "no synthetic-informed BH-FDR genes",
                "pathway_support": "no retained coherent synthetic-guided pathway",
                "interpretation": "no convincing synthetic-guided biological result",
            },
        ]
    )
    evidence_order = [
        "thymus",
        "soleus",
        "kidney",
        "skeletal_muscle",
        "spleen",
        "skin",
        "adrenal_gland",
        "lung",
        "liver",
    ]
    evidence["_order"] = evidence["tissue"].map(
        {tissue: index for index, tissue in enumerate(evidence_order)}
    )
    evidence = evidence.sort_values("_order").drop(columns="_order").reset_index(drop=True)

    tables = {
        "inventory": inventory,
        "pipeline_design": pipeline_design,
        "harmonization_summary": harmonization_summary,
        "harmonization_full": harmonization,
        "arch_summary": arch_summary,
        "locked_repeats": locked,
        "locked_summary": locked_summary,
        "model_screen": model_screen,
        "wgan_repeats": wgan_repeats,
        "naive_utility": naive_utility,
        "confirmation": confirmation,
        "genotype": genotype,
        "thymus_core": thymus_core,
        "thymus_reactome": thymus_reactome,
        "muscle_summary": muscle_summary,
        "soleus_genes": soleus_genes,
        "quadriceps_rbm6_effects": quadriceps_rbm6_effects,
        "quadriceps_rbm6_meta": quadriceps_rbm6_meta,
        "muscle_reactome": muscle_reactome,
        "tissue_summary": tissue_summary,
        "development_highlights": development_highlights,
        "study_holdout_context": study_holdout_context,
        "ordinary_fdr_genes": ordinary_fdr_genes,
        "all_bh_fdr_genes": all_bh_fdr_genes,
        "bh_fdr_tissue_summary": bh_fdr_tissue_summary,
        "thymus_evidence_mapping": thymus_evidence_mapping,
        "spleen_igfbp3_effects": spleen_igfbp3_effects,
        "spleen_igfbp3_meta": spleen_igfbp3_meta,
        "spleen_reference_expression": spleen_reference_expression,
        "evidence": evidence,
    }

    names = {
        "inventory": "table_1_data_inventory.tsv",
        "pipeline_design": "table_2_pipeline_design_space.tsv",
        "model_screen": "table_4_generator_model_selection.tsv",
        "locked_summary": "table_s24_locked_ddim_metric_summary.tsv",
        "confirmation": "table_s25_heldout_study_confirmation.tsv",
        "evidence": "table_7_tissue_evidence.tsv",
        "study_holdout_context": "table_6_whole_study_transfer_context.tsv",
        "arch_summary": "table_s1_archs4_ddim_metrics.tsv",
        "locked_repeats": "table_s2_locked_ddim_repeats.tsv",
        "naive_utility": "table_s3_naive_augmentation.tsv",
        "genotype": "table_s4_confirmation_genotypes.tsv",
        "thymus_core": "table_s5_thymus_core_genes.tsv",
        "thymus_reactome": "table_s6_thymus_reactome.tsv",
        "muscle_summary": "table_s7_muscle_group_summary.tsv",
        "soleus_genes": "table_s8_soleus_genes.tsv",
        "muscle_reactome": "table_s9_muscle_reactome.tsv",
        "tissue_summary": "table_s10_all_tissue_development_screen.tsv",
        "spleen_igfbp3_effects": "table_s11_spleen_igfbp3_accession_effects.tsv",
        "spleen_igfbp3_meta": "table_s12_spleen_igfbp3_random_effects.tsv",
        "spleen_reference_expression": "table_s13_spleen_reference_expression.tsv",
        "quadriceps_rbm6_effects": "table_s14_quadriceps_rbm6_accession_effects.tsv",
        "quadriceps_rbm6_meta": "table_s15_quadriceps_rbm6_random_effects.tsv",
        "ordinary_fdr_genes": "table_s16_ordinary_fdr_directional_genes.tsv",
        "all_bh_fdr_genes": "table_s17_all_random_effects_bh_fdr_genes.tsv",
        "bh_fdr_tissue_summary": "table_s18_bh_fdr_tissue_summary.tsv",
        "thymus_evidence_mapping": "table_s19_thymus_evidence_level_mapping.tsv",
        "development_highlights": "table_s20_tissue_utility_highlights.tsv",
        "harmonization_summary": "table_s21_liver_harmonization_benchmark.tsv",
        "harmonization_full": "table_s22_liver_harmonization_full_metrics.tsv",
        "wgan_repeats": "table_s23_wgan_validation_repeats.tsv",
    }
    for key, name in names.items():
        _write_tsv(tables[key], name)
    return tables


def _utility_display_name(value: str) -> str:
    labels = {
        "edl": "EDL",
        "skeletal_muscle": "Skeletal muscle, pooled",
        "tibialis_anterior": "Tibialis anterior",
    }
    return labels.get(value, value.replace("_", " ").capitalize())


def _utility_arm_name(value: str) -> str:
    labels = {
        "real_only": "Real only",
        "generated_only": "Generated only",
        "real_plus_generated": "Real + generated",
        "guided_real_only": "Guided ranking; real fit",
        "guided_low_weight": "Guided ranking; 5% synthetic",
    }
    if value not in labels:
        raise ValueError(f"Unexpected tissue utility arm: {value}")
    return labels[value]


def _utility_metric_pair(real: float, selected: float) -> str:
    return f"{float(real):.3f} / {float(selected):.3f}"


def _compact_utility_rows(frame: pd.DataFrame) -> list[list[str]]:
    metric_pairs = [
        ("real_mean_balanced_accuracy", "selected_mean_balanced_accuracy"),
        ("real_mean_roc_auc", "selected_mean_roc_auc"),
        ("real_mean_average_precision", "selected_mean_average_precision"),
    ]
    rows: list[list[str]] = []
    for row in frame.sort_values("tissue").itertuples(index=False):
        eligible = bool(row.generated_arm_eligible_all_metrics)
        deltas = [
            float(getattr(row, selected)) - float(getattr(row, real))
            for real, selected in metric_pairs
        ]
        if not eligible:
            status = "Real-only retained"
        elif all(math.isclose(delta, 0.0, abs_tol=1e-12) for delta in deltas):
            status = "Eligible tie"
        else:
            status = "Eligible improvement"
        rows.append(
            [
                _utility_display_name(str(row.tissue)),
                (
                    f"{int(row.development_profiles)} "
                    f"({int(row.flight)}/{int(row.ground_control)})"
                ),
                _utility_arm_name(str(row.selected_arm)),
                _utility_metric_pair(
                    row.real_mean_balanced_accuracy,
                    row.selected_mean_balanced_accuracy,
                ),
                _utility_metric_pair(
                    row.real_mean_roc_auc,
                    row.selected_mean_roc_auc,
                ),
                _utility_metric_pair(
                    row.real_mean_average_precision,
                    row.selected_mean_average_precision,
                ),
                status,
            ]
        )
    return rows


def _markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    def clean(value: str) -> str:
        return str(value).replace("|", "\\|").replace("\n", " ")

    rendered = [
        "| " + " | ".join(map(clean, headers)) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
    ]
    rendered.extend(
        "| " + " | ".join(clean(value) for value in row) + " |" for row in rows
    )
    return "\n".join(rendered)


def update_supplementary_utility_tables(tables: dict[str, pd.DataFrame]) -> None:
    canonical = tables["tissue_summary"]
    muscle = tables["muscle_summary"]
    if len(canonical) != 22 or len(muscle) != 5:
        raise ValueError(
            "Expected utility summaries for 22 canonical tissues and five "
            f"muscle groups, found {len(canonical)} and {len(muscle)}"
        )

    headers = [
        "Tissue",
        "n (FLT/GC)",
        "Selected arm",
        "BA real/selected",
        "AUROC real/selected",
        "AP real/selected",
        "Status",
    ]
    block = "\n\n".join(
        [
            (
                "All five arms were fitted for every analysis unit below. Values "
                "are means across eight repeated outer splits, and every outer "
                "evaluation used real profiles. An eligible arm was nonworse than "
                "real-only training in balanced accuracy, AUROC, and average "
                "precision. An eligible tie met that rule without improving a mean "
                "metric. These development results are not complete-study transfer "
                "tests."
            ),
            "**Supplementary Table S10. Complete canonical-tissue utility screen.**",
            _markdown_table(headers, _compact_utility_rows(canonical)),
            (
                "**Supplementary Table S7. Complete anatomical muscle-group "
                "utility screen.**"
            ),
            _markdown_table(headers, _compact_utility_rows(muscle)),
            (
                "Sample counts are shown as total development profiles followed by "
                "flight/ground-control counts. Small cohorts and ceiling-level "
                "scores remain exploratory even when a synthetic arm is eligible."
            ),
        ]
    )

    path = _required(PAPER_DIR / "supplementary_methods.md")
    text = path.read_text(encoding="utf-8")
    if text.count(UTILITY_TABLES_BEGIN) != 1 or text.count(UTILITY_TABLES_END) != 1:
        raise ValueError(
            "Supplementary utility-table markers are missing or duplicated"
        )
    prefix, remainder = text.split(UTILITY_TABLES_BEGIN, maxsplit=1)
    _, suffix = remainder.split(UTILITY_TABLES_END, maxsplit=1)
    path.write_text(
        (
            f"{prefix}{UTILITY_TABLES_BEGIN}\n\n{block}\n\n"
            f"{UTILITY_TABLES_END}{suffix}"
        ),
        encoding="utf-8",
    )


def figure_1_workflow() -> None:
    fig, ax = plt.subplots(figsize=(8.1, 6.5))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10.5)
    ax.axis("off")

    def box(
        x: float,
        y: float,
        width: float,
        height: float,
        title: str,
        body: str,
        color: str,
    ) -> None:
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.025,rounding_size=0.08",
            linewidth=1.2,
            edgecolor=color,
            facecolor="white",
        )
        ax.add_patch(patch)
        ax.text(
            x + 0.18,
            y + height - 0.23,
            title,
            color=color,
            weight="bold",
            va="top",
            fontsize=8.4,
        )
        ax.text(
            x + 0.18,
            y + height - (0.82 if "\n" in title else 0.62),
            body,
            color=COLORS["dark"],
            va="top",
            fontsize=6.8,
            linespacing=1.2,
        )

    def arrow(x1: float, y1: float, x2: float, y2: float) -> None:
        ax.annotate(
            "",
            xy=(x2, y2),
            xytext=(x1, y1),
            arrowprops={"arrowstyle": "-|>", "lw": 1.25, "color": COLORS["gray"]},
        )

    ax.text(0, 10.1, "A", weight="bold", fontsize=13)
    ax.text(0.45, 10.1, "Data sources", weight="bold", fontsize=11)
    box(
        0.2,
        7.8,
        3.35,
        1.65,
        "ARCHS4 mouse",
        "997,515 profiles audited\n17,244 selected across 20 tissues\nComplete GEO-series splits",
        COLORS["blue"],
    )
    box(
        4.0,
        7.8,
        3.35,
        1.65,
        "NASA OSDR API",
        "1,610 biological profiles\n75 accessions; FLT/GC labels\nStudy and material retained",
        COLORS["coral"],
    )
    box(
        7.8,
        7.8,
        4.0,
        1.65,
        "Biological scope",
        "Pooled multi-tissue generation\nTissue-specific analysis\nSkeletal-muscle groups retained",
        COLORS["green"],
    )

    ax.text(0, 7.25, "B", weight="bold", fontsize=13)
    ax.text(0.45, 7.25, "Configurable generative benchmark", weight="bold", fontsize=11)
    axes = [
        (
            "Expression",
            "Raw, CPM, TPM\nlog transforms\nz-score, robust, MaxAbs",
            COLORS["blue"],
        ),
        (
            "Harmonization",
            "None; two study z-scores\nComBat/ComBat-seq\nMBatch; MOBER",
            COLORS["teal"],
        ),
        (
            "Generator",
            "WGAN-GP; DDIM\nGeneJEPA screened as\nrepresentation only",
            COLORS["gold"],
        ),
        (
            "Training scope",
            "OSDR only; ARCHS4 only\nARCHS4 pretrain plus\nOSDR adaptation",
            COLORS["coral"],
        ),
        (
            "Cohort structure",
            "One or many studies\npooled or per tissue\naccession balancing",
            COLORS["purple"],
        ),
        (
            "Conditioning",
            "FLT/GC; tissue; study\nmaterial; muscle group\nand available covariates",
            COLORS["green"],
        ),
    ]
    positions = [(0.2, 5.25), (4.0, 5.25), (7.8, 5.25), (0.2, 3.15), (4.0, 3.15), (7.8, 3.15)]
    for (x, y), (title, body, color) in zip(positions, axes):
        box(x, y, 3.35, 1.65, title, body, color)

    ax.text(0, 2.45, "C", weight="bold", fontsize=13)
    ax.text(0.45, 2.45, "Staged selection and biological use", weight="bold", fontsize=11)
    stages = [
        (
            "Grouped splits",
            "GEO series or OSDR\naccessions kept intact",
            COLORS["blue"],
        ),
        (
            "Independent gates",
            "Fidelity, diversity,\nmemorization, effects",
            COLORS["gold"],
        ),
        (
            "Selected DDIM",
            "Only candidate passing\nfinal joint locked gates",
            COLORS["teal"],
        ),
        (
            "Synthetic-guided\nanalysis",
            "Per-tissue use selected;\nreal OSDR defines effects",
            COLORS["coral"],
        ),
    ]
    x_positions = [0.2, 3.15, 6.1, 9.05]
    for x, (title, body, color) in zip(x_positions, stages):
        box(x, 0.55, 2.55, 1.45, title, body, color)
    for left, right in zip(x_positions[:-1], x_positions[1:]):
        arrow(left + 2.55, 1.28, right, 1.28)
    ax.text(
        0.2,
        0.12,
        "The matrix defined a gated search, not an exhaustive Cartesian sweep. Synthetic profiles were never counted as additional animals.",
        fontsize=7.6,
        color=COLORS["gray"],
    )
    _save_figure(fig, "figure_1_study_design")


def figure_2_validation(tables: dict[str, pd.DataFrame]) -> None:
    arch = tables["arch_summary"].set_index("metric")["value"]
    locked = tables["locked_repeats"]
    model_screen = tables["model_screen"].set_index("model")
    harmonization = tables["harmonization_summary"]
    fig = plt.figure(figsize=(7.6, 6.7))
    grid = fig.add_gridspec(2, 2, width_ratios=[0.86, 1.14], hspace=0.45, wspace=0.34)

    ax = fig.add_subplot(grid[0, 0])
    values = [
        arch["Real train to real test tissue BA"],
        arch["Synthetic train to real test tissue BA"],
    ]
    bars = ax.bar(
        ["Real", "Synthetic"],
        values,
        color=[COLORS["gray"], COLORS["teal"]],
        width=0.62,
    )
    ax.set_ylim(0, 1)
    ax.set_ylabel("Balanced accuracy")
    ax.set_title("A  ARCHS4 tissue reverse validation", loc="left")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.025, f"{value:.3f}", ha="center")

    ax = fig.add_subplot(grid[0, 1])
    metric_columns = [
        ("Corr", "correlation"),
        ("Precision", "precision"),
        ("Recall", "recall"),
        ("F1", "f1"),
        ("AA", "adversarial_accuracy"),
        ("FD ratio", "frechet_ratio"),
    ]
    labels = [item[0] for item in metric_columns]
    x = np.arange(len(labels))
    width = 0.36
    wgan_values = [
        model_screen.loc["Study-conditioned WGAN-GP", column]
        for _, column in metric_columns
    ]
    ddim_values = [
        model_screen.loc["Factorized DDIM", column]
        for _, column in metric_columns
    ]
    ax.bar(
        x - width / 2,
        wgan_values,
        width,
        label="WGAN, validation",
        color=COLORS["coral"],
    )
    ax.bar(
        x + width / 2,
        ddim_values,
        width,
        label="DDIM, locked test",
        color=COLORS["teal"],
    )
    ax.axhspan(0.4, 0.6, xmin=0.66, xmax=0.84, color=COLORS["light"], zorder=0)
    ax.set_xticks(x, labels, rotation=25, ha="right", fontsize=7)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Metric value")
    ax.set_title("B  Generator candidates", loc="left")
    ax.legend(frameon=False, fontsize=7, loc="lower left")

    ax = fig.add_subplot(grid[1, 0])
    gate_labels = ["Corr", "Precision", "Recall", "F1", "AA", "FD", "FLT/GC", "Muscle"]
    gate_passes = [
        (locked["correlation"] >= locked["correlation_minimum"]).mean(),
        (locked["precision"] >= 0.95).mean(),
        (locked["recall"] >= 0.85).mean(),
        (locked["f1"] >= 0.90).mean(),
        locked["adversarial_accuracy"].between(0.40, 0.60).mean(),
        (locked["frechet_ratio"] <= 1.0).mean(),
        locked["condition_effect_pass"].astype(bool).mean(),
        locked["muscle_accession_pass"].astype(bool).mean(),
    ]
    y = np.arange(len(gate_labels))
    bars = ax.barh(
        y,
        gate_passes,
        color=[COLORS["teal"]] * 6 + [COLORS["gold"], COLORS["gold"]],
    )
    ax.set_yticks(y, gate_labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Fraction of four seeds passing")
    ax.set_title("C  Selected DDIM locked gates", loc="left")
    for bar, value in zip(bars, gate_passes):
        ax.text(value + 0.025, bar.get_y() + bar.get_height() / 2, f"{int(value * 4)}/4", va="center", fontsize=7)

    ax = fig.add_subplot(grid[1, 1])
    method_colors = [COLORS["gray"]] * len(harmonization)
    method_labels = harmonization["method"].astype(str).tolist()
    for index, label in enumerate(method_labels):
        if label == "Mentor two-stage z-score":
            method_colors[index] = COLORS["teal"]
        elif label == "MOBER (study)":
            method_colors[index] = COLORS["coral"]
        elif label == "No harmonization (TPM)":
            method_colors[index] = COLORS["blue"]
    ax.scatter(
        harmonization["correlation"],
        harmonization["f1"],
        c=method_colors,
        s=42,
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
    )
    for label in ["No harmonization (TPM)", "Mentor two-stage z-score", "MOBER (study)"]:
        row = harmonization.loc[harmonization["method"].eq(label)].iloc[0]
        short = {
            "No harmonization (TPM)": "None",
            "Mentor two-stage z-score": "Two-stage",
            "MOBER (study)": "MOBER",
        }[label]
        ax.annotate(
            short,
            (row["correlation"], row["f1"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )
    ax.axvline(0.98, color=COLORS["coral"], lw=1, ls="--")
    ax.axhline(0.90, color=COLORS["coral"], lw=1, ls="--")
    ax.set_xlim(-0.06, 1.03)
    ax.set_ylim(-0.04, 1.02)
    ax.set_xlabel("Correlation agreement")
    ax.set_ylabel("F1")
    ax.set_title("D  Matched liver harmonization", loc="left")

    fig.suptitle(
        "Staged benchmarking selected diffusion after WGAN and harmonization screens",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
    )
    fig.text(
        0.5,
        0.005,
        "WGAN values are validation results; DDIM values are locked-test results after staged selection, not a paired test-set comparison.",
        ha="center",
        fontsize=6.8,
        color=COLORS["gray"],
    )
    _save_figure(fig, "figure_2_generator_validation")


def figure_3_utility(tables: dict[str, pd.DataFrame]) -> None:
    naive = tables["naive_utility"]
    confirm = tables["confirmation"]
    genotype = tables["genotype"]
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.75), gridspec_kw={"width_ratios": [0.9, 1.25, 1.15]})

    ax = axes[0]
    x = np.arange(len(naive))
    width = 0.36
    ax.bar(x - width / 2, naive["balanced_accuracy"], width, label="BA", color=COLORS["teal"])
    ax.bar(x + width / 2, naive["roc_auc"], width, label="AUROC", color=COLORS["blue"])
    ax.set_xticks(x, ["Real", "Synth.", "Real +\nsynth."], fontsize=7)
    ax.set_ylim(0.55, 0.88)
    ax.set_title("A  Naive augmentation", loc="left")
    ax.legend(frameon=False, ncol=2, fontsize=7, loc="upper center")

    ax = axes[1]
    rows = []
    for _, row in confirm.iterrows():
        for metric, baseline, guided in [
            ("BA", "baseline_balanced_accuracy", "generated_balanced_accuracy"),
            ("AUROC", "baseline_roc_auc", "generated_roc_auc"),
            ("AP", "baseline_average_precision", "generated_average_precision"),
        ]:
            rows.append((row["tissue"], metric, row[baseline], row[guided]))
    plot = pd.DataFrame(rows, columns=["tissue", "metric", "baseline", "guided"])
    labels = [f"{row.tissue.title()} {row.metric}" for row in plot.itertuples()]
    y = np.arange(len(plot))
    ax.hlines(y, plot["baseline"], plot["guided"], color=COLORS["gray"], lw=1.2)
    ax.scatter(plot["baseline"], y, color=COLORS["gray"], label="Real-only baseline", s=28)
    ax.scatter(plot["guided"], y, color=COLORS["coral"], label="Synthetic-guided", s=30)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0.3, 1.01)
    ax.set_title("B  Held-out study tests", loc="left")
    ax.legend(frameon=False, fontsize=7, loc="lower right")

    ax = axes[2]
    metric_columns = {
        "BA": "delta_balanced_accuracy",
        "AUROC": "delta_roc_auc",
        "AP": "delta_average_precision",
    }
    ypos = np.arange(len(genotype))
    offsets = [-0.18, 0.0, 0.18]
    colors = [COLORS["teal"], COLORS["blue"], COLORS["gold"]]
    for offset, (metric, column), color in zip(offsets, metric_columns.items(), colors):
        ax.scatter(genotype[column], ypos + offset, label=metric, color=color, s=28)
    ax.axvline(0, color=COLORS["gray"], lw=1)
    ax.set_yticks(ypos, [f"{r.tissue.title()} {r.genotype}" for r in genotype.itertuples()])
    ax.invert_yaxis()
    ax.set_xlim(-0.08, 0.54)
    ax.set_xlabel("Guided minus baseline")
    ax.set_title("C  Genotype strata", loc="left")
    ax.legend(frameon=False, ncol=3, fontsize=7, loc="upper right")

    fig.suptitle(
        "Synthetic profiles did not improve naive augmentation, but fixed feature guidance transferred in thymus",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    _save_figure(fig, "figure_3_downstream_utility")


def figure_4_thymus(tables: dict[str, pd.DataFrame]) -> None:
    genes = tables["thymus_core"].sort_values("mean_real_effect")
    pathways = tables["thymus_reactome"].copy()
    selected_ids = [
        "R-MMU-174048_APC_C_CDC20_MEDIATED_DEGRADATION_OF_CYCLIN_B",
        "R-MMU-69478_G2_M_DNA_REPLICATION_CHECKPOINT",
        "R-MMU-69620_CELL_CYCLE_CHECKPOINTS",
        "R-MMU-69239_SYNTHESIS_OF_DNA",
        "R-MMU-69481_G2_M_CHECKPOINTS",
        "R-MMU-73894_DNA_REPAIR",
    ]
    pathways = pathways.loc[pathways["term"].isin(selected_ids)].sort_values("fdr", ascending=False)
    if len(pathways) != len(selected_ids):
        raise ValueError("A fixed thymus Reactome row is missing")

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.2), gridspec_kw={"width_ratios": [0.88, 1.32]})
    ax = axes[0]
    y = np.arange(len(genes))
    colors = [COLORS["coral"] if effect < 0 else COLORS["teal"] for effect in genes["mean_real_effect"]]
    ax.barh(y, genes["mean_real_effect"], color=colors)
    ax.axvline(0, color=COLORS["gray"], lw=0.9)
    ax.set_yticks(y, genes["symbol"])
    ax.set_xlabel("Held-out real FLT - GC effect")
    ax.set_title("A  Core genes in OSD-457", loc="left")
    for yi, effect in zip(y, genes["mean_real_effect"]):
        ax.text(
            effect + 0.035 if effect < 0 else effect + 0.025,
            yi,
            f"{effect:.2f}",
            va="center",
            ha="left",
            color="white" if effect < -0.15 else COLORS["dark"],
            fontsize=7,
        )

    ax = axes[1]
    y = np.arange(len(pathways))
    scores = pathways["overlap"] / pathways["pathway_genes_in_background"]
    ax.barh(y, scores, color=COLORS["gold"])
    labels = [_clean_term(term) for term in pathways["term"]]
    ax.set_yticks(y, labels)
    ax.set_xlim(0, max(scores) * 1.28)
    ax.set_xlabel("Fraction of pathway genes represented")
    ax.set_title("B  Shared cell-cycle processes", loc="left")
    for yi, score, overlap, total in zip(
        y,
        scores,
        pathways["overlap"],
        pathways["pathway_genes_in_background"],
    ):
        ax.text(score + 0.015, yi, f"{overlap}/{total} genes", va="center", fontsize=7)

    fig.suptitle(
        "Held-out thymus test prioritizes a flight-lower mitotic program in both genotypes",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    _save_figure(fig, "figure_4_thymus_biology")


def figure_5_soleus(tables: dict[str, pd.DataFrame]) -> None:
    genes = tables["soleus_genes"].sort_values("real_meta_effect")
    summary = tables["muscle_summary"].copy()
    pathways = tables["muscle_reactome"]
    pathways = pathways.loc[
        (pathways["tissue"] == "soleus")
        & (pathways["gene_set"] == "core_intersection")
        & (pathways["fdr"] < 0.05)
    ].nsmallest(4, "fdr").sort_values("fdr", ascending=False)

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(9.0, 4.15),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [0.9, 0.9, 1.3]},
    )
    ax = axes[0]
    y = np.arange(len(summary))
    counts = summary["cross_study_supported_genes"]
    bar_colors = [
        COLORS["gold"] if tissue == "soleus" else COLORS["blue"]
        for tissue in summary["tissue"]
    ]
    ax.barh(y, counts, color=bar_colors)
    ax.set_yticks(y, [name.replace("_", " ").title() for name in summary["tissue"]])
    ax.invert_yaxis()
    ax.set_xlim(0, max(counts.max() + 1.5, 2))
    ax.set_xlabel("Consistent genes")
    ax.set_title("A  Muscle groups", loc="left", fontsize=9)
    for yi, value in zip(y, counts):
        ax.text(value + 0.18, yi, f"{int(value)}", va="center", fontsize=7)

    ax = axes[1]
    y = np.arange(len(genes))
    colors = [COLORS["coral"] if effect < 0 else COLORS["teal"] for effect in genes["real_meta_effect"]]
    ax.barh(y, genes["real_meta_effect"], color=colors)
    ax.axvline(0, color=COLORS["gray"], lw=0.9)
    ax.set_yticks(y, genes["symbol"])
    ax.set_xlabel("Flight - ground")
    ax.set_title("B  Soleus genes", loc="left", fontsize=9)
    for yi, effect in zip(y, genes["real_meta_effect"]):
        if effect < -0.01:
            ax.text(
                effect / 2,
                yi,
                f"{effect:.3f}",
                va="center",
                ha="center",
                color="white",
                fontsize=6.8,
            )
        elif effect < 0:
            ax.text(
                effect - 0.002,
                yi,
                f"{effect:.3f}",
                va="center",
                ha="right",
                color=COLORS["dark"],
                fontsize=6.8,
            )
        else:
            ax.text(
                effect + 0.0015,
                yi,
                f"{effect:.3f}",
                va="center",
                ha="left",
                color=COLORS["dark"],
                fontsize=6.8,
            )

    ax = axes[2]
    y = np.arange(len(pathways))
    scores = pathways["overlap"] / pathways["pathway_genes_in_background"]
    ax.barh(y, scores, color=COLORS["purple"])
    pathway_labels = {
        "R-MMU-9837999_MITOCHONDRIAL_PROTEIN_DEGRADATION": "Mitochondrial protein turnover",
        "R-MMU-556833_METABOLISM_OF_LIPIDS": "Lipid metabolism",
        "R-MMU-77289_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION": "Mitochondrial FA beta oxidation",
        "R-MMU-8978868_FATTY_ACID_METABOLISM": "Fatty acid metabolism",
    }
    ax.set_yticks(y, [pathway_labels.get(term, _clean_term(term)) for term in pathways["term"]])
    ax.set_xlim(0, max(scores) * 1.35)
    ax.set_xlabel("Pathway genes represented")
    ax.set_title("C  Biological processes", loc="left", fontsize=9)
    for yi, score, overlap, total in zip(
        y,
        scores,
        pathways["overlap"],
        pathways["pathway_genes_in_background"],
    ):
        ax.text(score + 0.015, yi, f"{overlap}/{total} genes", va="center", fontsize=7)

    fig.suptitle(
        "Anatomical separation reveals a soleus-specific oxidative-metabolism hypothesis",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
    )
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.91))
    _save_figure(fig, "figure_5_soleus_biology")


def figure_6_evidence(tables: dict[str, pd.DataFrame]) -> None:
    evidence = tables["evidence"].copy()
    scores = evidence["tier_score"].to_numpy()
    palette = {
        0: "#B8C0C5",
        1: COLORS["gold"],
        2: COLORS["teal"],
        3: COLORS["coral"],
    }

    fig = plt.figure(figsize=(7.4, 7.6))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.92, 1.2],
        width_ratios=[1.12, 0.88],
        hspace=0.58,
        wspace=0.42,
    )

    ax = fig.add_subplot(grid[0, 0])
    performance = pd.concat(
        [tables["tissue_summary"], tables["muscle_summary"]],
        ignore_index=True,
        sort=False,
    ).drop_duplicates("tissue", keep="last")
    performance_order = [
        "soleus",
        "kidney",
        "skeletal_muscle",
        "spleen",
        "skin",
        "adrenal_gland",
    ]
    performance = (
        performance.set_index("tissue").loc[performance_order].reset_index()
    )
    y_performance = np.arange(len(performance))
    metric_columns = [
        ("mean_delta_balanced_accuracy", "Balanced accuracy", COLORS["teal"]),
        ("mean_delta_roc_auc", "AUROC", COLORS["blue"]),
        ("mean_delta_average_precision", "Average precision", COLORS["gold"]),
    ]
    offsets = [-0.22, 0.0, 0.22]
    for (column, label, color), offset in zip(metric_columns, offsets):
        ax.scatter(
            performance[column],
            y_performance + offset,
            color=color,
            s=28,
            label=label,
            zorder=3,
        )
    ax.axvline(0, color=COLORS["gray"], linewidth=0.9, linestyle="--")
    ax.grid(axis="x", color="#E5E9EB", linewidth=0.8)
    ax.set_yticks(
        y_performance,
        [name.replace("_", " ").title() for name in performance["tissue"]],
    )
    ax.invert_yaxis()
    ax.set_xlabel("Selected arm - real-only")
    ax.set_title("A  Repeated development-screen gains", loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=6.8, ncol=1, loc="lower right")

    ax = fig.add_subplot(grid[0, 1])
    genes = tables["ordinary_fdr_genes"]
    gene_order = [
        "thymus",
        "skeletal_muscle",
        "soleus",
        "kidney",
        "spleen",
        "skin",
        "adrenal_gland",
    ]
    gene_counts = (
        genes.loc[genes["tissue"].isin(gene_order)]
        .groupby(["tissue", "selection_interpretation"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(gene_order, fill_value=0)
    )
    promoted = gene_counts.get(
        "synthetic_promoted",
        pd.Series(0, index=gene_counts.index),
    )
    reinforced = gene_counts.get(
        "reinforced_real_and_synthetic",
        pd.Series(0, index=gene_counts.index),
    )
    y_genes = np.arange(len(gene_counts))
    ax.barh(y_genes, reinforced, color=COLORS["teal"], label="Reinforced")
    ax.barh(
        y_genes,
        promoted,
        left=reinforced,
        color=COLORS["coral"],
        label="Promoted",
    )
    ax.set_yticks(
        y_genes,
        [name.replace("_", " ").title() for name in gene_counts.index],
    )
    ax.invert_yaxis()
    ax.set_xlabel("BH-FDR genes")
    ax.set_title("B  Synthetic-informed genes", loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.grid(axis="x", color="#E5E9EB", linewidth=0.8)

    ax = fig.add_subplot(grid[1, :])
    y = np.arange(len(evidence))
    ax.scatter(scores, y, s=135, color=[palette[int(score)] for score in scores], zorder=3)
    for yi, score, theme in zip(y, scores, evidence["interpretation"]):
        ax.plot([score, 3.25], [yi, yi], color="#D9DEE1", lw=0.8, zorder=1)
        ax.text(3.35, yi, theme, va="center", fontsize=7.4, color=COLORS["dark"])
    ax.set_yticks(y, [name.title() for name in evidence["tissue"]])
    ax.invert_yaxis()
    ax.set_xlim(-0.25, 6.6)
    ax.set_xticks(
        [0, 1, 2, 3],
        ["No clear\nsignal", "Candidate", "Promising", "Strongest"],
    )
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.grid(axis="x", color="#E5E9EB", linewidth=0.8)
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.spines["left"].set_color("#D9DEE1")
    ax.set_title(
        "C  Tissue-level biological interpretation",
        loc="left",
        pad=14,
        fontsize=9,
    )
    fig.suptitle(
        "Evidence is strongest in thymus and soleus, with secondary tissue candidates",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
    )
    fig.subplots_adjust(top=0.92, bottom=0.06, left=0.16, right=0.98)
    _save_figure(fig, "figure_6_tissue_evidence")


def copy_supplementary_figures() -> None:
    copies = {
        ARCHS4_RUN / "evaluation/archs4_mouse_ddim_trajectory_pca.png": "figure_s1_archs4_denoising_trajectory.png",
        ARCHS4_RUN / "evaluation/archs4_mouse_ddim_trajectory_pca.pdf": "figure_s1_archs4_denoising_trajectory.pdf",
        LOCKED_DIR / "seed5020/real_vs_synthetic_pca.png": "figure_s2_locked_real_vs_synthetic_pca.png",
        LOCKED_DIR / "seed5020/real_vs_synthetic_pca.pdf": "figure_s2_locked_real_vs_synthetic_pca.pdf",
        MUSCLE_DIR / "arm_balanced_accuracy_heatmap.png": "figure_s3_muscle_arm_heatmap.png",
        MUSCLE_DIR / "arm_balanced_accuracy_heatmap.pdf": "figure_s3_muscle_arm_heatmap.pdf",
    }
    for source, target in copies.items():
        shutil.copy2(_required(source), FIGURE_DIR / target)


def build_manifest() -> None:
    tracked_inputs = [
        ARCHS4_RUN / "evaluation/summary.json",
        ARCHS4_RUN / "run_summary.json",
        LOCKED_DIR / "repeat_metrics.tsv",
        LOCKED_DIR / "summary.json",
        CONFIRM_DIR / "tissue_results.tsv",
        CONFIRM_DIR / "genotype_subgroup_results.tsv",
        CONFIRM_DIR / "thymus/feature_stability.tsv",
        CONFIRM_DIR / "thymus/reactome_enrichment.tsv",
        TRANSFER_DIR / "tissue_results.tsv",
        ADAPTIVE_HOLDOUT_DIR / "final_summary.json",
        FRESH_HOLDOUT_DIR / "final_summary.json",
        FRESH_HOLDOUT_DIR / "confirmatory_folds/aggregate_summary.json",
        MUSCLE_DIR / "tissue_arm_choices.tsv",
        MUSCLE_DIR / "paired_repeat_support.tsv",
        MUSCLE_DIR / "stable_gene_sets.tsv.gz",
        MUSCLE_DIR / "reactome_enrichment.tsv.gz",
        MUSCLE_DIR / "real_accession_effects.tsv.gz",
        MUSCLE_DIR / "real_random_effects.tsv.gz",
        TISSUE_DIR / "tissue_arm_choices.tsv",
        TISSUE_DIR / "paired_repeat_support.tsv",
        TISSUE_DIR / "biological_support_summary.tsv",
        TISSUE_DIR / "stable_gene_sets.tsv.gz",
        TISSUE_DIR / "real_accession_effects.tsv.gz",
        TISSUE_DIR / "real_random_effects.tsv.gz",
        LANDMARK_PANEL,
        WGAN_DIR / "summary.json",
        WGAN_DIR / "calibrated_repeat_metrics.tsv",
        HARMONIZATION_DIR / "independent_metrics.tsv",
        ROOT / "configs/generative/preprocessing_profiles.yaml",
        ROOT / "configs/generative/model_profiles.yaml",
        ROOT / "configs/generative/experiment_matrix.yaml",
        GENEJEPA_DIR / "figures/archs4_tissues_validation/summary.json",
    ]
    rows = []
    for path in tracked_inputs:
        _required(path)
        rows.append(
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )
    _write_tsv(pd.DataFrame(rows), "frozen_input_manifest.tsv")

    figure_rows = []
    for path in sorted(FIGURE_DIR.glob("*")):
        figure_rows.append(
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )
    _write_tsv(pd.DataFrame(figure_rows), "figure_build_manifest.tsv")


def render_document(markdown_path: Path, title: str) -> tuple[Path, Path]:
    try:
        import markdown
        from weasyprint import HTML
    except ImportError as error:
        raise RuntimeError(
            "Rendering requires the markdown and weasyprint packages in the nasa-mouse environment"
        ) from error

    body = markdown.markdown(
        markdown_path.read_text(encoding="utf-8"),
        extensions=["tables", "fenced_code", "sane_lists", "md_in_html"],
    )
    css = """
    @page {
      size: A4;
      margin: 19mm 17mm 19mm 17mm;
      @bottom-center {
        content: counter(page);
        color: #68737a;
        font-size: 8pt;
      }
    }
    @page:first { @bottom-center { content: none; } }
    html { font-family: "DejaVu Sans", Arial, sans-serif; color: #23313a; }
    body { font-size: 9.2pt; line-height: 1.43; }
    h1 { color: #23445d; font-size: 23pt; line-height: 1.15; margin: 0 0 8mm 0; }
    h2 { color: #23445d; font-size: 14pt; margin: 7mm 0 2.5mm 0; border-bottom: 0.5pt solid #b9c5ca; padding-bottom: 1.2mm; }
    h3 { color: #17807e; font-size: 10.8pt; margin: 5mm 0 1.5mm 0; }
    h4 { color: #23313a; font-size: 9.6pt; margin: 4mm 0 1mm 0; }
    p { margin: 0 0 2.7mm 0; text-align: justify; hyphens: auto; }
    ul, ol { margin: 1.5mm 0 3mm 5mm; padding-left: 4mm; }
    li { margin-bottom: 1mm; }
    a { color: #236a8e; text-decoration: none; }
    code { font-family: "DejaVu Sans Mono", monospace; font-size: 8pt; background: #f0f3f4; padding: 0.2mm 0.5mm; }
    pre { white-space: pre-wrap; background: #f0f3f4; padding: 3mm; border-left: 2pt solid #17807e; font-size: 7.6pt; }
    table { width: 100%; border-collapse: collapse; margin: 3mm 0 5mm 0; font-size: 7.4pt; line-height: 1.25; break-before: avoid; break-inside: avoid; }
    th { background: #e8eef0; color: #23445d; font-weight: bold; text-align: left; padding: 1.5mm; border-bottom: 1pt solid #7b858c; }
    td { padding: 1.35mm; border-bottom: 0.45pt solid #ced6d9; vertical-align: top; }
    img { max-width: 100%; max-height: 215mm; display: block; margin: 3mm auto 2mm auto; }
    .figure-block { break-inside: avoid; }
    .figure-composite { break-inside: avoid; margin: 3mm 0 2mm 0; }
    .figure-composite img { width: 88%; max-height: none; margin: 1mm auto; }
    blockquote { margin: 3mm 0; padding: 2.5mm 4mm; border-left: 2.5pt solid #d69a2d; background: #f7f4ea; }
    .title-page { min-height: 238mm; display: flex; flex-direction: column; justify-content: center; page-break-after: always; }
    .title-page h1 { font-size: 25pt; }
    .title-page p { text-align: left; }
    .subtitle { font-size: 12pt; color: #4d5f69; }
    .authors { margin-top: 9mm; font-size: 11pt; font-weight: bold; }
    .affiliation { color: #4d5f69; }
    .draft-note { margin-top: 15mm; padding: 3mm; border-left: 3pt solid #d96552; background: #f7ece9; }
    .caption { font-size: 8pt; line-height: 1.35; color: #3c4b53; text-align: left; margin-bottom: 5mm; }
    .page-break { page-break-before: always; }
    .small { font-size: 8pt; color: #4d5f69; }
    """
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{title}</title>
<style>{css}</style></head><body>{body}</body></html>"""
    html_path = markdown_path.with_suffix(".html")
    pdf_path = markdown_path.with_suffix(".pdf")
    html_path.write_text(html_text, encoding="utf-8")
    HTML(string=html_text, base_url=str(PAPER_DIR)).write_pdf(str(pdf_path))
    return html_path, pdf_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-render",
        action="store_true",
        help="Build source tables and figures without rendering HTML/PDF.",
    )
    args = parser.parse_args()

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    _style()
    tables = build_source_tables()
    update_supplementary_utility_tables(tables)
    figure_1_workflow()
    figure_2_validation(tables)
    figure_3_utility(tables)
    figure_4_thymus(tables)
    figure_5_soleus(tables)
    figure_6_evidence(tables)
    copy_supplementary_figures()
    build_manifest()

    if not args.skip_render:
        render_document(
            _required(PAPER_DIR / "manuscript.md"),
            (
                "A configurable generative transcriptomics framework reveals "
                "thymic proliferative suppression and soleus metabolic "
                "remodeling in spaceflown mice"
            ),
        )
        render_document(
            _required(PAPER_DIR / "supplementary_methods.md"),
            (
                "Supplementary methods: configurable generative transcriptomics "
                "in spaceflown mice"
            ),
        )

    print(f"Built paper package: {PAPER_DIR}")


if __name__ == "__main__":
    main()
