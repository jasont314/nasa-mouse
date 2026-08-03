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

from .annotate_promoted_gene_literature import write_tables as write_literature_tables


ROOT = Path(__file__).resolve().parents[2]
PAPER_DIR = ROOT / "paper" / "synthetic_guided_spaceflight"
FIGURE_DIR = PAPER_DIR / "figures"
SOURCE_DIR = PAPER_DIR / "source_data"
UTILITY_TABLES_BEGIN = "<!-- BEGIN GENERATED TISSUE UTILITY TABLES -->"
UTILITY_TABLES_END = "<!-- END GENERATED TISSUE UTILITY TABLES -->"
LANDMARK_PANEL = ROOT / "data/diffusion/l974_mouse_paper_parity.tsv"
OBSOLETE_PAPER_ARTIFACTS = (
    "figures/figure_1_study_design.png",
    "figures/figure_1_study_design.pdf",
    "figures/figure_2_generator_validation.png",
    "figures/figure_2_generator_validation.pdf",
    "figures/figure_3a_archs4_denoising_trajectory.png",
    "figures/figure_3a_archs4_denoising_trajectory.pdf",
    "figures/figure_3b_locked_real_vs_synthetic_pca.png",
    "figures/figure_3b_locked_real_vs_synthetic_pca.pdf",
    "figures/figure_4_thymus_biology.png",
    "figures/figure_4_thymus_biology.pdf",
    "figures/figure_5_soleus_biology.png",
    "figures/figure_5_soleus_biology.pdf",
    "figures/figure_6_tissue_evidence.png",
    "figures/figure_6_tissue_evidence.pdf",
    "figures/figure_3_downstream_utility.png",
    "figures/figure_3_downstream_utility.pdf",
    "figures/figure_s1_archs4_denoising_trajectory.png",
    "figures/figure_s1_archs4_denoising_trajectory.pdf",
    "figures/figure_s2_locked_real_vs_synthetic_pca.png",
    "figures/figure_s2_locked_real_vs_synthetic_pca.pdf",
    "figures/figure_s3_muscle_arm_heatmap.png",
    "figures/figure_s3_muscle_arm_heatmap.pdf",
    "figures/figure_s5_whole_study_transfer.png",
    "figures/figure_s5_whole_study_transfer.pdf",
    "figures/figure_s6_effect_recovery_levels.png",
    "figures/figure_s6_effect_recovery_levels.pdf",
    "source_data/table_6_whole_study_transfer_context.tsv",
    "source_data/table_7_tissue_evidence.tsv",
    "source_data/table_s4_confirmation_genotypes.tsv",
    "source_data/table_s5_thymus_core_genes.tsv",
    "source_data/table_s6_thymus_reactome.tsv",
    "source_data/table_s7_muscle_group_summary.tsv",
    "source_data/table_s8_soleus_genes.tsv",
    "source_data/table_s9_muscle_reactome.tsv",
    "source_data/table_s10_all_tissue_development_screen.tsv",
    "source_data/table_s11_spleen_igfbp3_accession_effects.tsv",
    "source_data/table_s12_spleen_igfbp3_random_effects.tsv",
    "source_data/table_s13_spleen_reference_expression.tsv",
    "source_data/table_s14_quadriceps_rbm6_accession_effects.tsv",
    "source_data/table_s15_quadriceps_rbm6_random_effects.tsv",
    "source_data/table_s16_ordinary_fdr_directional_genes.tsv",
    "source_data/table_s17_all_random_effects_bh_fdr_genes.tsv",
    "source_data/table_s18_bh_fdr_tissue_summary.tsv",
    "source_data/table_s19_thymus_evidence_level_mapping.tsv",
    "source_data/table_s20_tissue_utility_highlights.tsv",
    "source_data/table_s21_liver_harmonization_benchmark.tsv",
    "source_data/table_s22_liver_harmonization_full_metrics.tsv",
    "source_data/table_s23_wgan_validation_repeats.tsv",
    "source_data/table_s24_locked_ddim_metric_summary.tsv",
    "source_data/table_s25_heldout_study_confirmation.tsv",
    "source_data/table_s26_prior_transfer_experiments.tsv",
    "source_data/table_s27_whole_study_tissue_effect_recovery.tsv",
    "source_data/table_s28_whole_study_accession_effect_recovery.tsv",
    "source_data/table_s29_whole_study_pooled_effect_recovery.tsv",
)

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
    muscle_choices = _read_tsv(MUSCLE_DIR / "tissue_arm_choices.tsv")
    muscle_repeats = _read_tsv(MUSCLE_DIR / "paired_repeat_support.tsv")
    muscle_genes = _gate_synthetic_selection(
        _read_tsv(MUSCLE_DIR / "stable_gene_sets.tsv.gz"),
        muscle_choices,
    )
    muscle_reactome = _read_tsv(MUSCLE_DIR / "reactome_enrichment.tsv.gz")
    muscle_inventory = _read_tsv(MUSCLE_DIR / "tissue_inventory.tsv")
    muscle_random_effects = _read_tsv(MUSCLE_DIR / "real_random_effects.tsv.gz")
    tissue_choices = _read_tsv(TISSUE_DIR / "tissue_arm_choices.tsv")
    tissue_repeats = _read_tsv(TISSUE_DIR / "paired_repeat_support.tsv")
    tissue_biology = _read_tsv(TISSUE_DIR / "biological_support_summary.tsv")
    tissue_inventory = _read_tsv(TISSUE_DIR / "tissue_inventory.tsv")
    tissue_genes = _gate_synthetic_selection(
        _read_tsv(TISSUE_DIR / "stable_gene_sets.tsv.gz"),
        tissue_choices,
    )
    tissue_reactome = _read_tsv(TISSUE_DIR / "reactome_enrichment.tsv.gz")
    tissue_random_effects = _read_tsv(TISSUE_DIR / "real_random_effects.tsv.gz")
    landmark_panel = _read_tsv(LANDMARK_PANEL)

    _assert_close(
        arch_eval["synthetic_to_real_test_tissue_classifier"]["balanced_accuracy"],
        0.7810085910974481,
        "ARCHS4 reverse-validation balanced accuracy",
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
                "analysis_split": "781 train / 536 validation / 293 test",
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
                "configurable_options": "WGAN-GP; DDIM",
                "evaluated_scope": "paper-reproduced architectures and independent metrics",
                "selected_branch": "ARCHS4-pretrained, OSDR-adapted DDIM",
            },
            {
                "axis": "Validation",
                "configurable_options": (
                    "GEO-series or accession-grouped evaluation; test set; "
                    "unconditional controls; multiple generation seeds"
                ),
                "evaluated_scope": "no sample-random model-selection split",
                "selected_branch": "four-seed 293-profile OSDR test",
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
                    "retained as tissue-conditioned initialization; correlation "
                    "agreement was below the target"
                ),
            },
            {
                "model": "Study-conditioned WGAN-GP",
                "training_regime": "OSDR matched study-conditioned",
                "evaluation_split": "536-profile validation",
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
                    "not used downstream: higher external separability and no "
                    "accession-aware muscle-effect recovery"
                ),
            },
            {
                "model": "Factorized DDIM",
                "training_regime": "ARCHS4 pretraining then OSDR adaptation",
                "evaluation_split": "293-profile within-study OSDR test",
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
                    "used downstream: lower adversarial accuracy and Frechet ratio "
                    "with muscle-effect recovery in 4/4 repeats"
                ),
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

    core_symbols = [
        "Nusap1",
        "Stmn1",
        "Birc5",
        "Cdk1",
        "Top2a",
        "Ccnb2",
        "Aurka",
        "Ccne2",
        "Ube2c",
        "Gmnn",
    ]
    thymus_core = tissue_genes.loc[
        tissue_genes["tissue"].eq("thymus")
        & tissue_genes["symbol"].isin(core_symbols)
        & tissue_genes["stable_generated"].astype(bool)
        & tissue_genes["real_meta_fdr"].lt(0.05)
        & tissue_genes["real_effect_supports_generated"].astype(bool)
    ].copy()
    thymus_core = thymus_core.sort_values("real_meta_effect")
    if set(thymus_core["symbol"]) != set(core_symbols):
        missing = sorted(set(core_symbols) - set(thymus_core["symbol"]))
        raise ValueError(f"Missing development-screen thymus core genes: {missing}")
    thymus_reactome = tissue_reactome.loc[
        tissue_reactome["tissue"].eq("thymus")
        & tissue_reactome["gene_set"].eq("generated_supported")
        & tissue_reactome["fdr"].lt(0.05)
    ].copy()

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

    evidence = pd.DataFrame(
        [
            {
                "tissue": "thymus",
                "tier": "coherent development panel",
                "tier_score": 2,
                "predictive_result": "guided delta BA/AUROC/AP +0.111/+0.069/+0.056",
                "real_gene_support": "13 promoted and 3 reinforced BH-FDR genes",
                "pathway_support": "mitosis, G2/M, APC/C, DNA replication; FDR < 0.05",
                "interpretation": "coherent flight-lower proliferative-renewal hypothesis",
            },
            {
                "tissue": "soleus",
                "tier": "cross-accession development",
                "tier_score": 2,
                "predictive_result": "real-plus-generated delta BA/AUROC/AP +0.038/+0.000/+0.006",
                "real_gene_support": "5 reinforced BH-FDR genes; 4 pass LOO FDR",
                "pathway_support": "mitochondrial lipid oxidation/protein turnover; FDR < 0.05",
                "interpretation": "coherent metabolic hypothesis; requires independent confirmation",
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
                "interpretation": "pooled-muscle development complements anatomical soleus result",
            },
            {
                "tissue": "lung",
                "tier": "predictive development only",
                "tier_score": 1,
                "predictive_result": "generated-only delta BA/AUROC/AP +0.078/+0.150/+0.148",
                "real_gene_support": "no BH-FDR gene in the 974-gene panel",
                "pathway_support": "no Reactome term at FDR < 0.05",
                "interpretation": "development gain without a supported biological panel",
            },
            {
                "tissue": "spleen",
                "tier": "developmental exploratory",
                "tier_score": 1,
                "predictive_result": "real-plus-generated delta BA/AUROC/AP +0.131/+0.163/+0.160",
                "real_gene_support": "Rai14, Ptprk, Myl9 promoted; Loxl1 reinforced; none pass LOO",
                "pathway_support": "no coherent stable-set Reactome enrichment",
                "interpretation": "adhesion/cytoskeletal hypothesis",
            },
            {
                "tissue": "skin",
                "tier": "developmental exploratory",
                "tier_score": 1,
                "predictive_result": "real-plus-generated delta BA/AUROC/AP +0.085/+0.077/+0.061",
                "real_gene_support": "Plscr1 is promoted and FLT-up in 6/6 studies; not LOO-stable",
                "pathway_support": "cell-cycle/DNA-repair theme matches published skin analyses",
                "interpretation": "literature-aligned exploratory biological support",
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
        "thymus_core": thymus_core,
        "thymus_reactome": thymus_reactome,
        "muscle_summary": muscle_summary,
        "soleus_genes": soleus_genes,
        "muscle_reactome": muscle_reactome,
        "tissue_summary": tissue_summary,
        "development_highlights": development_highlights,
        "ordinary_fdr_genes": ordinary_fdr_genes,
        "all_bh_fdr_genes": all_bh_fdr_genes,
        "bh_fdr_tissue_summary": bh_fdr_tissue_summary,
        "evidence": evidence,
    }

    names = {
        "inventory": "table_1_data_inventory.tsv",
        "pipeline_design": "table_2_pipeline_design_space.tsv",
        "model_screen": "table_4_generator_model_selection.tsv",
        "evidence": "table_6_tissue_evidence.tsv",
        "arch_summary": "table_s1_archs4_ddim_metrics.tsv",
        "locked_repeats": "table_s2_locked_ddim_repeats.tsv",
        "naive_utility": "table_s3_naive_augmentation.tsv",
        "thymus_core": "table_s4_thymus_core_genes.tsv",
        "thymus_reactome": "table_s5_thymus_reactome.tsv",
        "muscle_summary": "table_s6_muscle_group_summary.tsv",
        "soleus_genes": "table_s7_soleus_genes.tsv",
        "muscle_reactome": "table_s8_muscle_reactome.tsv",
        "tissue_summary": "table_s9_all_tissue_development_screen.tsv",
        "ordinary_fdr_genes": "table_s10_synthetic_informed_bh_fdr_genes.tsv",
        "all_bh_fdr_genes": "table_s11_all_random_effects_bh_fdr_genes.tsv",
        "bh_fdr_tissue_summary": "table_s12_bh_fdr_tissue_summary.tsv",
        "harmonization_summary": "table_s13_liver_harmonization_benchmark.tsv",
        "harmonization_full": "table_s14_liver_harmonization_full_metrics.tsv",
        "wgan_repeats": "table_s15_wgan_validation_repeats.tsv",
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
            "**Supplementary Table S9. Complete canonical-tissue utility screen.**",
            _markdown_table(headers, _compact_utility_rows(canonical)),
            '<div class="page-break"></div>',
            (
                "**Supplementary Table S6. Complete anatomical muscle-group "
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
    fig, ax = plt.subplots(figsize=(8.1, 9.1))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 14.5)
    ax.axis("off")
    top = 3.9

    def box(
        x: float,
        y: float,
        width: float,
        height: float,
        title: str,
        body: str,
        color: str,
        title_fontsize: float = 8.4,
        body_fontsize: float = 6.8,
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
            fontsize=title_fontsize,
        )
        ax.text(
            x + 0.18,
            y + height - (0.82 if "\n" in title else 0.62),
            body,
            color=COLORS["dark"],
            va="top",
            fontsize=body_fontsize,
            linespacing=1.2,
        )

    def arrow(x1: float, y1: float, x2: float, y2: float) -> None:
        ax.annotate(
            "",
            xy=(x2, y2),
            xytext=(x1, y1),
            arrowprops={"arrowstyle": "-|>", "lw": 1.25, "color": COLORS["gray"]},
        )

    ax.text(0, 10.1 + top, "A", weight="bold", fontsize=13)
    ax.text(0.45, 10.1 + top, "Data sources", weight="bold", fontsize=11)
    box(
        0.2,
        7.8 + top,
        3.35,
        1.65,
        "ARCHS4 mouse",
        "997,515 profiles audited\n17,244 selected across 20 tissues\nComplete GEO-series splits",
        COLORS["blue"],
    )
    box(
        4.0,
        7.8 + top,
        3.35,
        1.65,
        "NASA OSDR API",
        "1,610 biological profiles\n75 accessions; FLT/GC labels\nStudy and material retained",
        COLORS["coral"],
    )
    box(
        7.8,
        7.8 + top,
        4.0,
        1.65,
        "Biological scope",
        "Pooled multi-tissue generation\nTissue-specific analysis\nSkeletal-muscle groups retained",
        COLORS["green"],
    )

    ax.text(0, 7.25 + top, "B", weight="bold", fontsize=13)
    ax.text(
        0.45,
        7.25 + top,
        "Configurable generative benchmark",
        weight="bold",
        fontsize=11,
    )
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
            "WGAN-GP; DDIM\npaper-based architectures\ncommon evaluation contract",
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
    positions = [
        (0.2, 5.25 + top),
        (4.0, 5.25 + top),
        (7.8, 5.25 + top),
        (0.2, 3.15 + top),
        (4.0, 3.15 + top),
        (7.8, 3.15 + top),
    ]
    for (x, y), (title, body, color) in zip(positions, axes):
        box(x, y, 3.35, 1.65, title, body, color)

    ax.text(0, 2.45 + top, "C", weight="bold", fontsize=13)
    ax.text(
        0.45,
        2.45 + top,
        "Generator metrics and model choice",
        weight="bold",
        fontsize=11,
    )
    comparisons = [
        (
            0.2,
            3.35,
            "WGAN-GP",
            "Corr. 0.976; F1 0.985\nAA 0.636; muscle recovery 0/6",
            COLORS["coral"],
        ),
        (
            4.0,
            3.35,
            "DDIM",
            "Corr. 0.974; F1 0.997\nAA 0.475; muscle recovery 4/4",
            COLORS["teal"],
        ),
        (
            7.8,
            4.0,
            "Downstream model",
            "Lower AA and FD; muscle recovery 4/4\nDDIM used for downstream analysis",
            COLORS["green"],
        ),
    ]
    for x, width, title, body, color in comparisons:
        box(x, 0.55 + top, width, 1.45, title, body, color)
    arrow(7.35, 1.28 + top, 7.8, 1.28 + top)
    ax.text(
        0.2,
        0.12 + top,
        "AA near 0.5 indicates lower separability; generated profiles were never counted as additional animals.",
        fontsize=7.6,
        color=COLORS["gray"],
    )

    ax.text(0, 3.52, "D", weight="bold", fontsize=13)
    ax.text(
        0.45,
        3.52,
        "Five tissue-specific uses of generated expression",
        weight="bold",
        fontsize=11,
    )
    ax.text(0.2, 3.15, "Training views", weight="bold", fontsize=7.2)
    box(
        0.2,
        2.05,
        2.0,
        0.9,
        "Real OSDR",
        "Observed FLT/GC profiles",
        COLORS["blue"],
        title_fontsize=7.2,
        body_fontsize=5.7,
    )
    box(
        0.2,
        0.95,
        2.0,
        0.9,
        "DDIM generated",
        "Matched conditional profiles",
        COLORS["coral"],
        title_fontsize=7.2,
        body_fontsize=5.7,
    )

    ax.text(2.75, 3.15, "Candidate arm fitted within each tissue", weight="bold", fontsize=7.2)
    arms = [
        (
            2.75,
            2.0,
            "Real only",
            "Rank: real\nFit: real",
            COLORS["blue"],
        ),
        (
            4.62,
            2.0,
            "Generated only",
            "Rank: generated\nFit: generated",
            COLORS["coral"],
        ),
        (
            6.49,
            2.0,
            "Real + synth.",
            "Rank: consensus\nFit: equal real/synth.",
            COLORS["teal"],
        ),
        (
            3.69,
            0.83,
            "Guided; real fit",
            "Rank: consensus\nFit: real",
            COLORS["purple"],
        ),
        (
            5.56,
            0.83,
            "Guided; 5% synth.",
            "Rank: consensus\nFit: real + 5% synth.",
            COLORS["gold"],
        ),
    ]
    for x, y, title, body, color in arms:
        box(
            x,
            y,
            1.68,
            0.95,
            title,
            body,
            color,
            title_fontsize=6.4,
            body_fontsize=5.2,
        )

    arrow(2.2, 2.5, 2.56, 1.97)
    arrow(2.2, 1.4, 2.56, 1.83)
    arrow(2.58, 1.9, 2.68, 1.9)
    arrow(8.2, 1.9, 8.75, 1.9)
    box(
        8.8,
        0.85,
        3.0,
        2.15,
        "Within-study real-profile\nevaluation",
        (
            "BA, AUROC, AP on real holdouts\n"
            "Profiles split within accession\n"
            "Eligible arm selected per tissue"
        ),
        COLORS["green"],
        title_fontsize=7.4,
        body_fontsize=5.8,
    )
    arrow(10.3, 0.83, 10.3, 0.61)
    endpoint = FancyBboxPatch(
        (2.75, 0.06),
        9.05,
        0.5,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        linewidth=1.1,
        edgecolor=COLORS["coral"],
        facecolor="white",
    )
    ax.add_patch(endpoint)
    ax.text(
        7.275,
        0.31,
        (
            "Selected arm -> repeated stable genes; FLT/GC effects and BH-FDR "
            "are calculated from real OSDR only"
        ),
        ha="center",
        va="center",
        fontsize=5.9,
        color=COLORS["dark"],
        weight="bold",
    )
    _save_figure(fig, "figure_1_study_design")


def figure_1_validation(tables: dict[str, pd.DataFrame]) -> None:
    arch = tables["arch_summary"].set_index("metric")["value"]
    locked = tables["locked_repeats"]
    model_screen = tables["model_screen"].set_index("model")
    fig = plt.figure(figsize=(7.6, 5.2))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.05, 0.82],
        width_ratios=[0.86, 1.14],
        hspace=0.52,
        wspace=0.34,
    )

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
        label="DDIM, OSDR test",
        color=COLORS["teal"],
    )
    ax.axhspan(0.4, 0.6, xmin=0.66, xmax=0.84, color=COLORS["light"], zorder=0)
    ax.set_xticks(x, labels, rotation=25, ha="right", fontsize=7)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Metric value")
    ax.set_title("B  Generator metrics", loc="left")
    ax.legend(frameon=False, fontsize=7, loc="lower left")

    ax = fig.add_subplot(grid[1, :])
    gate_labels = ["Corr", "Precision", "Recall", "F1", "AA", "FD", "FLT/GC"]
    gate_passes = [
        (locked["correlation"] >= locked["correlation_minimum"]).mean(),
        (locked["precision"] >= 0.95).mean(),
        (locked["recall"] >= 0.85).mean(),
        (locked["f1"] >= 0.90).mean(),
        locked["adversarial_accuracy"].between(0.40, 0.60).mean(),
        (locked["frechet_ratio"] <= 1.0).mean(),
        locked["condition_effect_pass"].astype(bool).mean(),
    ]
    y = np.arange(len(gate_labels))
    bars = ax.barh(
        y,
        gate_passes,
        color=[COLORS["teal"]] * 6 + [COLORS["gold"]],
    )
    ax.set_yticks(y, gate_labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Fraction of four seeds passing")
    ax.set_title("C  DDIM repeat performance", loc="left")
    for bar, value in zip(bars, gate_passes):
        ax.text(value + 0.025, bar.get_y() + bar.get_height() / 2, f"{int(value * 4)}/4", va="center", fontsize=7)

    fig.suptitle(
        "Generator metrics support DDIM use in downstream analysis",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
    )
    fig.text(
        0.5,
        0.005,
        "WGAN values use validation data and DDIM values use the stated test data; metrics are not paired on one common split.",
        ha="center",
        fontsize=6.8,
        color=COLORS["gray"],
    )
    _save_figure(fig, "figure_1_generator_validation")


def figure_s2_utility(tables: dict[str, pd.DataFrame]) -> None:
    naive = tables["naive_utility"]
    development = tables["development_highlights"].copy()
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.4, 4.6),
        gridspec_kw={"width_ratios": [0.85, 1.65]},
    )

    ax = axes[0]
    x = np.arange(len(naive))
    width = 0.36
    ax.bar(x - width / 2, naive["balanced_accuracy"], width, label="BA", color=COLORS["teal"])
    ax.bar(x + width / 2, naive["roc_auc"], width, label="AUROC", color=COLORS["blue"])
    ax.set_xticks(x, ["Real", "Synth.", "Real +\nsynth."], fontsize=7)
    ax.set_ylim(0.55, 0.88)
    ax.set_title("A  Naive augmentation", loc="left")
    ax.legend(frameon=False, ncol=2, fontsize=7, loc="upper center")

    metric_columns = {
        "BA": ("real_mean_balanced_accuracy", "selected_mean_balanced_accuracy"),
        "AUROC": ("real_mean_roc_auc", "selected_mean_roc_auc"),
        "AP": ("real_mean_average_precision", "selected_mean_average_precision"),
    }
    ypos = np.arange(len(development))
    offsets = [-0.18, 0.0, 0.18]
    colors = [COLORS["teal"], COLORS["blue"], COLORS["gold"]]
    ax = axes[1]
    for offset, (metric, columns), color in zip(offsets, metric_columns.items(), colors):
        real_column, selected_column = columns
        delta = development[selected_column] - development[real_column]
        ax.scatter(delta, ypos + offset, label=metric, color=color, s=28)
    ax.axvline(0, color=COLORS["gray"], lw=1)
    ax.set_yticks(
        ypos,
        [_utility_display_name(tissue) for tissue in development["tissue"]],
    )
    ax.invert_yaxis()
    ax.set_xlabel("Selected arm minus real-only")
    ax.set_title("B  Tissue-specific development", loc="left")
    ax.legend(frameon=False, ncol=3, fontsize=7, loc="upper right")

    fig.suptitle(
        "Pooled augmentation failed, while tissue-specific synthetic use varied",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    _save_figure(fig, "figure_s2_downstream_utility")


def figure_3_thymus(tables: dict[str, pd.DataFrame]) -> None:
    genes = tables["thymus_core"].sort_values("real_meta_effect")
    pathways = tables["thymus_reactome"].copy()
    selected_ids = [
        "R-MMU-69278_CELL_CYCLE_MITOTIC",
        "R-MMU-1640170_CELL_CYCLE",
        "R-MMU-69239_SYNTHESIS_OF_DNA",
        "R-MMU-69306_DNA_REPLICATION",
        "R-MMU-68886_M_PHASE",
        "R-MMU-174048_APC_C_CDC20_MEDIATED_DEGRADATION_OF_CYCLIN_B",
    ]
    pathways = pathways.loc[pathways["term"].isin(selected_ids)].sort_values("fdr", ascending=False)
    if len(pathways) != len(selected_ids):
        raise ValueError("A fixed thymus Reactome row is missing")

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.2), gridspec_kw={"width_ratios": [0.88, 1.32]})
    ax = axes[0]
    y = np.arange(len(genes))
    colors = [
        COLORS["coral"]
        if label == "generated_supported"
        else COLORS["teal"]
        for label in genes["gene_set"]
    ]
    ax.barh(y, genes["real_meta_effect"], color=colors)
    ax.axvline(0, color=COLORS["gray"], lw=0.9)
    ax.set_yticks(y, genes["symbol"])
    ax.set_xlabel("Random-effects FLT - GC estimate")
    ax.set_title("A  Cross-study thymus effects", loc="left")
    for yi, effect in zip(y, genes["real_meta_effect"]):
        ax.text(
            effect + 0.002,
            yi,
            f"{effect:.3f}",
            va="center",
            ha="left",
            color=COLORS["dark"],
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
        "Synthetic-informed thymus genes converge on a flight-lower mitotic program",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    _save_figure(fig, "figure_3_thymus_biology")


def figure_4_soleus(tables: dict[str, pd.DataFrame]) -> None:
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
    _save_figure(fig, "figure_4_soleus_biology")


def figure_5_evidence(tables: dict[str, pd.DataFrame]) -> None:
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
    ax.set_yticks(
        y,
        [name.replace("_", " ").title() for name in evidence["tissue"]],
    )
    ax.invert_yaxis()
    ax.set_xlim(-0.25, 6.6)
    ax.set_xticks(
        [0, 1, 2],
        ["No clear\nsignal", "Exploratory", "Coherent\ndevelopment"],
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
        "Synthetic-guided development and real-data support identify tissue priorities",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
    )
    fig.subplots_adjust(top=0.92, bottom=0.06, left=0.16, right=0.98)
    _save_figure(fig, "figure_5_tissue_evidence")


def copy_publication_figures() -> None:
    copies = {
        ARCHS4_RUN / "evaluation/archs4_mouse_ddim_trajectory_pca.png": "figure_2a_archs4_denoising_trajectory.png",
        ARCHS4_RUN / "evaluation/archs4_mouse_ddim_trajectory_pca.pdf": "figure_2a_archs4_denoising_trajectory.pdf",
        LOCKED_DIR / "seed5020/real_vs_synthetic_pca.png": "figure_2b_locked_real_vs_synthetic_pca.png",
        LOCKED_DIR / "seed5020/real_vs_synthetic_pca.pdf": "figure_2b_locked_real_vs_synthetic_pca.pdf",
        MUSCLE_DIR / "arm_balanced_accuracy_heatmap.png": "figure_s1_muscle_arm_heatmap.png",
        MUSCLE_DIR / "arm_balanced_accuracy_heatmap.pdf": "figure_s1_muscle_arm_heatmap.pdf",
    }
    for source, target in copies.items():
        shutil.copy2(_required(source), FIGURE_DIR / target)


def build_manifest() -> None:
    tracked_inputs = [
        ARCHS4_RUN / "evaluation/summary.json",
        ARCHS4_RUN / "run_summary.json",
        LOCKED_DIR / "repeat_metrics.tsv",
        LOCKED_DIR / "summary.json",
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
        TISSUE_DIR / "reactome_enrichment.tsv.gz",
        TISSUE_DIR / "real_accession_effects.tsv.gz",
        TISSUE_DIR / "real_random_effects.tsv.gz",
        LANDMARK_PANEL,
        WGAN_DIR / "summary.json",
        WGAN_DIR / "calibrated_repeat_metrics.tsv",
        HARMONIZATION_DIR / "independent_metrics.tsv",
        ROOT / "configs/generative/preprocessing_profiles.yaml",
        ROOT / "configs/generative/model_profiles.yaml",
        ROOT / "configs/generative/experiment_matrix.yaml",
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
    .figure-composite img.trajectory-panel { width: 112%; max-width: none; margin: 1mm 0 1mm -6%; }
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
    for relative_path in OBSOLETE_PAPER_ARTIFACTS:
        (PAPER_DIR / relative_path).unlink(missing_ok=True)
    _style()
    tables = build_source_tables()
    write_literature_tables()
    update_supplementary_utility_tables(tables)
    figure_1_validation(tables)
    figure_s2_utility(tables)
    figure_3_thymus(tables)
    figure_4_soleus(tables)
    figure_5_evidence(tables)
    copy_publication_figures()
    build_manifest()

    if not args.skip_render:
        render_document(
            _required(PAPER_DIR / "manuscript.md"),
            (
                "A configurable generative transcriptomics framework identifies "
                "tissue-dependent synthetic utility in mouse spaceflight RNA-seq"
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
