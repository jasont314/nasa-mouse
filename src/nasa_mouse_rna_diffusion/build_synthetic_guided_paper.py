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
LANDMARK_PANEL = ROOT / "data/diffusion/l974_mouse_paper_parity.tsv"

ARCHS4_RUN = (
    ROOT
    / "outputs/generative_benchmark/runs/lacan_diffusion/"
    "archs4_mouse_paper_parity_osdr_disjoint_seed1234"
)
OSDR_RUN = (
    ROOT
    / "outputs/generative_benchmark/runs/lacan_diffusion/"
    "osdr_factorized_study_lora512_correlation_refine_seed2020"
)
LOCKED_DIR = OSDR_RUN / "evaluation/final_locked_test"
CONFIRM_DIR = (
    ROOT
    / "outputs/generative_benchmark/analyses/"
    "generated_feature_guidance_confirmation_disjoint_v1"
)
MUSCLE_DIR = (
    ROOT
    / "outputs/generative_benchmark/analyses/"
    "within_study_generated_feature_stability_muscle_groups_v1"
)
TISSUE_DIR = (
    ROOT
    / "outputs/generative_benchmark/analyses/"
    "within_study_generated_feature_stability_v1"
)
WGAN_DIR = (
    ROOT
    / "outputs/generative_benchmark/runs/vinas_wgan_gp/"
    "osdr_matched_study_conditioned_seed2020/evaluation/matched_validation"
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
    confirmation = _read_tsv(CONFIRM_DIR / "tissue_results.tsv")
    genotype = _read_tsv(CONFIRM_DIR / "genotype_subgroup_results.tsv")
    thymus_features = _read_tsv(CONFIRM_DIR / "thymus/feature_stability.tsv")
    thymus_reactome = _read_tsv(CONFIRM_DIR / "thymus/reactome_enrichment.tsv")
    muscle_choices = _read_tsv(MUSCLE_DIR / "tissue_arm_choices.tsv")
    muscle_repeats = _read_tsv(MUSCLE_DIR / "paired_repeat_support.tsv")
    muscle_genes = _read_tsv(MUSCLE_DIR / "stable_gene_sets.tsv.gz")
    muscle_reactome = _read_tsv(MUSCLE_DIR / "reactome_enrichment.tsv.gz")
    muscle_inventory = _read_tsv(MUSCLE_DIR / "tissue_inventory.tsv")
    muscle_accession_effects = _read_tsv(MUSCLE_DIR / "real_accession_effects.tsv.gz")
    muscle_random_effects = _read_tsv(MUSCLE_DIR / "real_random_effects.tsv.gz")
    tissue_choices = _read_tsv(TISSUE_DIR / "tissue_arm_choices.tsv")
    tissue_repeats = _read_tsv(TISSUE_DIR / "paired_repeat_support.tsv")
    tissue_biology = _read_tsv(TISSUE_DIR / "biological_support_summary.tsv")
    tissue_genes = _read_tsv(TISSUE_DIR / "stable_gene_sets.tsv.gz")
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

    model_screen = pd.DataFrame(
        [
            {
                "model": "ARCHS4-pretrained DDIM",
                "evaluation": "held-out ARCHS4 series",
                "correlation": arch_eval["gene_correlation_matrix_agreement"],
                "precision": arch_eval["precision_recall_in_scaled_l974"]["precision"],
                "recall": arch_eval["precision_recall_in_scaled_l974"]["recall"],
                "adversarial_accuracy": arch_eval[
                    "nearest_neighbor_adversarial_accuracy_in_scaled_l974"
                ],
                "decision": (
                    "retained for tissue-conditioned initialization; strict "
                    "correlation-matrix gate failed"
                ),
            },
            {
                "model": "OSDR-adapted factorized DDIM",
                "evaluation": "293-profile locked within-study test",
                "correlation": locked["correlation"].mean(),
                "precision": locked["precision"].mean(),
                "recall": locked["recall"].mean(),
                "adversarial_accuracy": locked["adversarial_accuracy"].mean(),
                "decision": "accepted for represented-study conditional simulation",
            },
            {
                "model": "OSDR study-conditioned WGAN-GP",
                "evaluation": "validation only; locked test unopened",
                "correlation": 0.9759,
                "precision": 0.9764,
                "recall": 0.9938,
                "adversarial_accuracy": 0.6362,
                "decision": "rejected; externally separable and no accession-aware effect recovery",
            },
            {
                "model": "GeneJEPA exact-architecture duration screen",
                "evaluation": "held-out ARCHS4 series",
                "correlation": np.nan,
                "precision": np.nan,
                "recall": np.nan,
                "adversarial_accuracy": np.nan,
                "decision": "representation only; no expression decoder; not advanced",
            },
        ]
    )

    naive_utility = pd.DataFrame(
        [
            ("Real OSDR", 0.754, 0.819),
            ("Synthetic", 0.700, 0.751),
            ("Real + synthetic", 0.734, 0.801),
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
        & muscle_genes["real_loo_fdr_stable_0_05"].astype(bool)
        & muscle_genes["real_effect_supports_generated"].astype(bool)
    ].copy()
    soleus_genes = soleus_genes.sort_values("real_meta_effect")
    if len(soleus_genes) != 7:
        raise ValueError(f"Expected 7 generated-supported soleus LOO genes, found {len(soleus_genes)}")

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
        int(selection_counts.get("synthetic_promoted", 0)) != 28
        or int(selection_counts.get("reinforced_real_and_synthetic", 0)) != 24
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
                "tier": "leakage-corrected held-out confirmation",
                "tier_score": 3,
                "predictive_result": "BA 0.500 to 0.833; AUROC 0.840 to 0.979",
                "real_gene_support": "8 core FLT-down genes; genotype effect r=0.975",
                "pathway_support": "G2/M, APC/C, DNA replication; Reactome FDR < 0.05",
                "interpretation": (
                    "leakage-corrected feature-guidance result; outcomes were "
                    "known before retraining"
                ),
            },
            {
                "tissue": "soleus",
                "tier": "cross-accession development",
                "tier_score": 2,
                "predictive_result": "generated-only delta BA/AUROC/AP +0.025/+0.020/+0.020",
                "real_gene_support": "8 unanimous ordinary-FDR genes; 7 pass LOO FDR",
                "pathway_support": "mitochondrial lipid oxidation/protein turnover; FDR < 0.05",
                "interpretation": "coherent hypothesis; requires unseen-accession confirmation",
            },
            {
                "tissue": "quadriceps",
                "tier": "cross-accession gene-level evidence",
                "tier_score": 2,
                "predictive_result": "guided delta BA/AUROC/AP +0.050/+0.040/+0.026",
                "real_gene_support": "4 unanimous ordinary-FDR genes; Rbm6 passes LOO FDR",
                "pathway_support": "G1/S and TP53 terms were suggestive; minimum FDR 0.0598",
                "interpretation": "secondary four-gene association; mechanism unconfirmed",
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
                "tier": "cross-accession gene-level evidence",
                "tier_score": 2,
                "predictive_result": "nested delta BA/AUROC/AP +0.170/+0.208/+0.204",
                "real_gene_support": "5 unanimous ordinary-FDR genes; Igfbp3 passes LOO FDR",
                "pathway_support": "no coherent stable-set Reactome enrichment",
                "interpretation": "strong single-gene association; stromal mechanism unconfirmed",
            },
            {
                "tissue": "skin",
                "tier": "developmental exploratory",
                "tier_score": 1,
                "predictive_result": "nested gains did not reproduce on the reserved profiles",
                "real_gene_support": "Plscr1 is FLT-up in 6/6 studies but is not synthetic-selected",
                "pathway_support": "cell-cycle/DNA-repair theme matches published skin analyses",
                "interpretation": "literature-aligned heterogeneous response; not a new synthetic claim",
            },
            {
                "tissue": "kidney",
                "tier": "developmental exploratory",
                "tier_score": 1,
                "predictive_result": "nested delta BA/AUROC/AP +0.029/+0.093/+0.097",
                "real_gene_support": "Slc37a4 is reinforced and FLT-up in 6/6 studies",
                "pathway_support": "renal glucose metabolism; lipid/ECM context from expiMap",
                "interpretation": "credible secondary hypothesis; requires unseen-study confirmation",
            },
            {
                "tissue": "liver",
                "tier": "negative",
                "tier_score": 0,
                "predictive_result": "small nested gains; 5/8 repeats nonworse",
                "real_gene_support": "0 selected genes pass real LOO FDR",
                "pathway_support": "no retained coherent synthetic-guided pathway",
                "interpretation": "no convincing synthetic-guided biological result",
            },
            {
                "tissue": "retina",
                "tier": "developmental exploratory",
                "tier_score": 1,
                "predictive_result": "nested gains but no real LOO-stable selected genes",
                "real_gene_support": "Slc37a4 is FLT-up in 4/4 studies; ordinary FDR 0.0238",
                "pathway_support": "shear-stress enrichment involves a different gene set",
                "interpretation": "exploratory gene/pathway mismatch; no integrated claim",
            },
        ]
    )
    evidence_order = [
        "thymus",
        "soleus",
        "quadriceps",
        "spleen",
        "kidney",
        "lung",
        "skin",
        "liver",
        "retina",
    ]
    evidence["_order"] = evidence["tissue"].map(
        {tissue: index for index, tissue in enumerate(evidence_order)}
    )
    evidence = evidence.sort_values("_order").drop(columns="_order").reset_index(drop=True)

    tables = {
        "inventory": inventory,
        "arch_summary": arch_summary,
        "locked_repeats": locked,
        "locked_summary": locked_summary,
        "model_screen": model_screen,
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
        "model_screen": "table_2_model_screen.tsv",
        "locked_summary": "table_3_locked_ddim_metrics.tsv",
        "confirmation": "table_4_leakage_corrected_confirmation.tsv",
        "evidence": "table_5_tissue_evidence.tsv",
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
    }
    for key, name in names.items():
        _write_tsv(tables[key], name)
    return tables


def figure_1_workflow() -> None:
    fig, ax = plt.subplots(figsize=(8.1, 5.2))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 9)
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
        ax.text(x + 0.18, y + height - 0.28, title, color=color, weight="bold", va="top")
        ax.text(
            x + 0.18,
            y + height - (0.95 if "\n" in title else 0.72),
            body,
            color=COLORS["dark"],
            va="top",
            fontsize=7.4,
            linespacing=1.25,
        )

    def arrow(x1: float, y1: float, x2: float, y2: float) -> None:
        ax.annotate(
            "",
            xy=(x2, y2),
            xytext=(x1, y1),
            arrowprops={"arrowstyle": "-|>", "lw": 1.25, "color": COLORS["gray"]},
        )

    ax.text(0, 8.55, "A", weight="bold", fontsize=13)
    ax.text(0.45, 8.55, "Data and synthetic-guided discovery", weight="bold", fontsize=11)
    box(
        0.2,
        5.25,
        2.55,
        2.55,
        "ARCHS4 mouse",
        "17,244 reference profiles\n20 tissue classes\nLearns baseline tissue\nexpression",
        COLORS["blue"],
    )
    box(
        3.15,
        5.25,
        2.55,
        2.55,
        "Conditional\ndiffusion",
        "Generates tissue-aware\nmouse expression\nAdapts to spaceflight\ncohorts",
        COLORS["teal"],
    )
    box(
        6.1,
        5.25,
        2.55,
        2.55,
        "NASA OSDR API",
        "1,610 flight/control\nprofiles\n75 studies\n24 material classes",
        COLORS["coral"],
    )
    box(
        9.05,
        5.25,
        2.55,
        2.55,
        "Biological\nanalysis",
        "Synthetic profiles\nprioritize genes\nEffects measured in\nreal samples",
        COLORS["gold"],
    )
    arrow(2.75, 6.5, 3.15, 6.5)
    arrow(5.7, 6.5, 6.1, 6.5)
    arrow(8.65, 6.5, 9.05, 6.5)

    ax.text(0, 4.45, "B", weight="bold", fontsize=13)
    ax.text(0.45, 4.45, "Biological questions and findings", weight="bold", fontsize=11)
    levels = [
        (
            "1. Tissue\nresponse",
            "Which organs show\nflight-associated\nexpression changes?",
            COLORS["blue"],
        ),
        (
            "2. Thymus",
            "Cell division and\nproliferative renewal",
            COLORS["teal"],
        ),
        (
            "3. Skeletal\nmuscle",
            "Soleus metabolism,\nmitochondrial turnover,\nand contractile identity",
            COLORS["gold"],
        ),
        (
            "4. Other\ntissues",
            "Splenic IGFBP3,\nrenal metabolism,\nand heterogeneous\norgan responses",
            COLORS["coral"],
        ),
    ]
    x_positions = [0.2, 3.15, 6.1, 9.05]
    for x, (title, body, color) in zip(x_positions, levels):
        box(x, 1.45, 2.55, 2.25, title, body, color)
    for left, right in zip(x_positions[:-1], x_positions[1:]):
        arrow(left + 2.55, 2.55, right, 2.55)
    ax.text(
        0.2,
        0.55,
        "Generated expression prioritizes hypotheses; real flight and ground-control samples determine the biology.",
        fontsize=8.5,
        weight="bold",
        color=COLORS["dark"],
    )
    _save_figure(fig, "figure_1_study_design")


def figure_2_validation(tables: dict[str, pd.DataFrame]) -> None:
    arch = tables["arch_summary"].set_index("metric")["value"]
    locked = tables["locked_repeats"]
    fig = plt.figure(figsize=(7.4, 6.2))
    grid = fig.add_gridspec(2, 2, width_ratios=[0.92, 1.08], hspace=0.38, wspace=0.32)

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
    metrics = ["Gene mean", "Gene SD", "Correlation", "Precision", "Recall"]
    values = [
        arch["Gene mean correlation"],
        arch["Gene SD correlation"],
        arch["Gene correlation agreement"],
        arch["Precision, scaled L974"],
        arch["Recall, scaled L974"],
    ]
    y = np.arange(len(metrics))
    ax.barh(y, values, color=[COLORS["blue"], COLORS["blue"], COLORS["gold"], COLORS["teal"], COLORS["teal"]])
    ax.set_yticks(y, metrics)
    ax.invert_yaxis()
    ax.set_xlim(0.75, 1.01)
    ax.axvline(0.98, color=COLORS["coral"], lw=1, ls="--")
    ax.set_xlabel("Agreement or neighborhood fraction")
    ax.set_title("B  Broad-reference fidelity", loc="left")
    for yi, value in zip(y, values):
        ax.text(value - 0.004, yi, f"{value:.3f}", ha="right", va="center", color="white", fontsize=8)

    selected = [
        ("Correlation", "correlation", (0.94, 1.005), 0.949716),
        ("Precision", "precision", (0.94, 1.005), 0.95),
        ("Recall", "recall", (0.94, 1.005), 0.85),
        ("F1", "f1", (0.94, 1.005), 0.90),
    ]
    ax = fig.add_subplot(grid[1, 0])
    for idx, (label, column, _, threshold) in enumerate(selected):
        vals = locked[column].to_numpy()
        ax.scatter(vals, np.full_like(vals, idx), color=COLORS["teal"], s=28, zorder=3)
        ax.plot([vals.min(), vals.max()], [idx, idx], color=COLORS["teal"], lw=2)
        ax.plot(threshold, idx, marker="|", color=COLORS["coral"], markersize=12, mew=1.5)
    ax.set_yticks(range(len(selected)), [item[0] for item in selected])
    ax.invert_yaxis()
    ax.set_xlim(0.84, 1.005)
    ax.set_xlabel("Four locked generation seeds")
    ax.set_title("C  OSDR locked fidelity", loc="left")

    ax = fig.add_subplot(grid[1, 1])
    effect_metrics = [
        ("AA", "adversarial_accuracy", 0.5),
        ("Pooled FLT/GC r", "condition_delta_correlation", 0.3),
        ("Muscle accession r", "muscle_accession_correlation", 0.3),
    ]
    for idx, (label, column, target) in enumerate(effect_metrics):
        vals = locked[column].to_numpy()
        ax.scatter(vals, np.full_like(vals, idx), color=COLORS["coral"] if idx == 0 else COLORS["gold"], s=28, zorder=3)
        ax.plot([vals.min(), vals.max()], [idx, idx], color=COLORS["gray"], lw=2)
        ax.plot(target, idx, marker="|", color=COLORS["navy"], markersize=12, mew=1.5)
    ax.axvspan(0.4, 0.6, ymin=0.72, ymax=1.0, color=COLORS["light"], zorder=0)
    ax.set_yticks(range(len(effect_metrics)), [item[0] for item in effect_metrics])
    ax.invert_yaxis()
    ax.set_xlim(0.25, 0.82)
    ax.set_xlabel("Accuracy or effect correlation")
    ax.set_title("D  Indistinguishability and condition recovery", loc="left")

    fig.suptitle(
        "The ARCHS4 DDIM preserves tissue information and the adapted DDIM passes locked distribution gates",
        x=0.02,
        ha="left",
        fontsize=11,
        weight="bold",
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
    ax.set_title("B  Leakage-corrected held-out studies", loc="left")
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
        "Independent thymus confirmation prioritizes a flight-lower mitotic program in both genotypes",
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
    ].sort_values("fdr", ascending=False)

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
    themes = [
        "Lower cell division and proliferative renewal",
        "Oxidative metabolism and contractile remodeling",
        "FLT-higher four-gene set; G1/S and TP53 suggestive",
        "FLT-higher IGFBP3; stromal-niche hypothesis",
        "Renal glucose metabolism; lipid/ECM context",
        "Cell cycle, senescence, and PI3K/AKT (mixed)",
        "Cell cycle and DNA repair candidates",
        "No coherent retained pattern",
        "FLT-higher SLC37A4; pathway context differs",
    ]
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
    effects = tables["spleen_igfbp3_effects"].copy()
    meta = tables["spleen_igfbp3_meta"].iloc[0]
    scale = 1_000.0
    y_effects = np.arange(len(effects))
    ax.errorbar(
        effects["flight_minus_ground"] * scale,
        y_effects,
        xerr=1.96 * effects["standard_error"] * scale,
        fmt="o",
        color=COLORS["teal"],
        ecolor="#86AAA8",
        elinewidth=1.1,
        capsize=2.5,
        markersize=5,
        zorder=3,
    )
    meta_y = len(effects) + 0.35
    ax.errorbar(
        float(meta["meta_effect"]) * scale,
        meta_y,
        xerr=1.96 * float(meta["meta_se"]) * scale,
        fmt="D",
        color=COLORS["coral"],
        ecolor=COLORS["coral"],
        elinewidth=1.4,
        capsize=3,
        markersize=5.5,
        zorder=4,
    )
    effect_labels = [
        f"{row.accession}  ({int(row.n_flight)}/{int(row.n_ground_control)})"
        for row in effects.itertuples()
    ]
    ax.set_yticks(
        [*y_effects, meta_y],
        [*effect_labels, "Random-effects estimate"],
    )
    ax.invert_yaxis()
    ax.axvline(0, color=COLORS["gray"], linewidth=0.9, linestyle="--")
    ax.grid(axis="x", color="#E5E9EB", linewidth=0.8)
    ax.set_xlabel(r"FLT - GC ($\times 10^{-3}$, model scale)")
    ax.set_title("A  Cross-study spleen $Igfbp3$ effect", loc="left", fontsize=9)

    ax = fig.add_subplot(grid[0, 1])
    reference = tables["spleen_reference_expression"].iloc[::-1].copy()
    y_reference = np.arange(len(reference))
    ax.barh(
        y_reference,
        reference["mean_igfbp3_rpkm"],
        color=[COLORS["coral"], COLORS["gold"], COLORS["teal"], COLORS["blue"]],
    )
    ax.set_xscale("log")
    ax.set_xlim(1, 2_500)
    ax.set_yticks(y_reference, reference["population"])
    ax.set_xlabel("Mean $Igfbp3$ (RPKM, log scale)")
    ax.set_title(
        "B  Healthy spleen source (GSE156162)",
        loc="left",
        fontsize=9,
    )
    ax.grid(axis="x", color="#E5E9EB", linewidth=0.8)
    for yi, value in zip(y_reference, reference["mean_igfbp3_rpkm"]):
        ax.text(
            float(value) * 1.12,
            yi,
            f"{float(value):,.1f}",
            va="center",
            fontsize=6.8,
            color=COLORS["dark"],
        )

    ax = fig.add_subplot(grid[1, :])
    y = np.arange(len(evidence))
    ax.scatter(scores, y, s=135, color=[palette[int(score)] for score in scores], zorder=3)
    for yi, score, theme in zip(y, scores, themes):
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
        "Spleen $Igfbp3$ adds a focused cross-study result to the tissue response hierarchy",
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
        CONFIRM_DIR / "tissue_results.tsv",
        CONFIRM_DIR / "genotype_subgroup_results.tsv",
        CONFIRM_DIR / "thymus/feature_stability.tsv",
        CONFIRM_DIR / "thymus/reactome_enrichment.tsv",
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
    table { width: 100%; border-collapse: collapse; margin: 3mm 0 5mm 0; font-size: 7.4pt; line-height: 1.25; break-inside: avoid; }
    th { background: #e8eef0; color: #23445d; font-weight: bold; text-align: left; padding: 1.5mm; border-bottom: 1pt solid #7b858c; }
    td { padding: 1.35mm; border-bottom: 0.45pt solid #ced6d9; vertical-align: top; }
    img { max-width: 100%; max-height: 215mm; display: block; margin: 3mm auto 2mm auto; }
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
            "Cross-study synthetic-guided transcriptomics of spaceflown mice",
        )
        render_document(
            _required(PAPER_DIR / "supplementary_methods.md"),
            "Supplementary methods: synthetic-guided spaceflight transcriptomics",
        )

    print(f"Built paper package: {PAPER_DIR}")


if __name__ == "__main__":
    main()
