"""Integrate reviewer-directed sensitivity analyses into an evidence matrix."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .build_asgsr_paper import CONFIGS, FIGURE_DIR, SOURCE_DIR


COMPOSITION_MINIMUM_EFFECT_RATIO = 0.25
HELDOUT_MINIMUM_CONCORDANCE = 2 / 3


def build_evidence_matrix() -> pd.DataFrame:
    benchmark = pd.read_csv(SOURCE_DIR / "table_s13_method_benchmark.tsv", sep="\t")
    heldout = pd.read_csv(
        SOURCE_DIR / "table_s15_project_heldout_predictions.tsv.gz", sep="\t"
    )
    seeds = pd.read_csv(SOURCE_DIR / "table_s21_training_seed_consensus.tsv", sep="\t")
    composition = pd.read_csv(
        SOURCE_DIR / "table_s19_composition_proxy_adjusted_effects.tsv", sep="\t"
    )

    heldout = (
        heldout.loc[heldout["method"].eq("expimap")]
        .groupby(["tissue", "term"])["direction_concordant"]
        .agg(
            heldout_project_direction_concordance="mean",
            heldout_projects="size",
        )
        .reset_index()
    )
    seed_columns = [
        "tissue",
        "term",
        "effect_seed2020",
        "effect_seed2021",
        "effect_seed2022",
        "all_three_seeds_available",
        "all_three_direction_concordant",
    ]
    composition_columns = [
        "tissue",
        "term",
        "within_accession_unadjusted_effect",
        "composition_proxy_adjusted_effect",
        "adjusted_direction_matches_unadjusted",
        "absolute_effect_ratio_adjusted_to_unadjusted",
    ]
    matrix = (
        benchmark.loc[benchmark["curated_pathway"]]
        .merge(heldout, on=["tissue", "term"], how="left")
        .merge(seeds[seed_columns], on=["tissue", "term"], how="left")
        .merge(composition[composition_columns], on=["tissue", "term"], how="left")
    )
    matrix["ssgsea_direction_support"] = matrix[
        "expimap_ssgsea_direction_match"
    ].fillna(False)
    matrix["preranked_gsea_direction_support"] = matrix[
        "expimap_gsea_direction_match"
    ].fillna(False)
    matrix["heldout_direction_support"] = matrix[
        "heldout_project_direction_concordance"
    ].ge(HELDOUT_MINIMUM_CONCORDANCE)
    matrix["seed_direction_support"] = (
        matrix["all_three_seeds_available"].fillna(False)
        & matrix["all_three_direction_concordant"].fillna(False)
    )
    matrix["composition_proxy_support"] = (
        matrix["adjusted_direction_matches_unadjusted"].fillna(False)
        & matrix["absolute_effect_ratio_adjusted_to_unadjusted"].ge(
            COMPOSITION_MINIMUM_EFFECT_RATIO
        )
    )
    conventional_count = matrix[
        ["ssgsea_direction_support", "preranked_gsea_direction_support"]
    ].sum(axis=1)
    internal = (
        matrix["heldout_direction_support"]
        & matrix["seed_direction_support"]
        & matrix["composition_proxy_support"]
    )
    matrix["robustness_status"] = "sensitivity-dependent"
    matrix.loc[
        internal & conventional_count.lt(2), "robustness_status"
    ] = "internally robust, incomplete conventional support"
    matrix.loc[
        internal & conventional_count.eq(2), "robustness_status"
    ] = "triangulated"
    method_supported_sensitive = (
        conventional_count.eq(2)
        & matrix["heldout_direction_support"]
        & ~internal
    )
    matrix.loc[
        method_supported_sensitive, "robustness_status"
    ] = "method-supported, model-sensitive"
    matrix["robustness_support_count"] = matrix[
        [
            "ssgsea_direction_support",
            "preranked_gsea_direction_support",
            "heldout_direction_support",
            "seed_direction_support",
            "composition_proxy_support",
        ]
    ].sum(axis=1)
    return matrix.sort_values(
        ["tissue", "robustness_support_count", "short_label"],
        ascending=[True, False, True],
    )


def plot_evidence_matrix(matrix: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(19, 8.5), constrained_layout=True)
    columns = [
        "ssgsea_direction_support",
        "preranked_gsea_direction_support",
        "heldout_direction_support",
        "seed_direction_support",
        "composition_proxy_support",
    ]
    labels = ["ssGSEA", "GSEA", "Held-out\nproject", "3 seeds", "Composition\nproxy"]
    for ax, config in zip(axes, CONFIGS):
        frame = matrix.loc[matrix["tissue"].eq(config.tissue)].copy()
        status_order = {
            "triangulated": 0,
            "internally robust, incomplete conventional support": 1,
            "method-supported, model-sensitive": 2,
            "sensitivity-dependent": 3,
        }
        frame["status_order"] = frame["robustness_status"].map(status_order)
        frame = frame.sort_values(
            ["status_order", "robustness_support_count", "short_label"],
            ascending=[True, False, True],
        )
        values = frame[columns].astype(float).to_numpy()
        color = np.zeros((*values.shape, 3), dtype=float)
        color[values == 1] = np.array([0.18, 0.55, 0.34])
        color[values == 0] = np.array([0.78, 0.32, 0.28])
        ax.imshow(color, aspect="auto", interpolation="nearest")
        ax.set_xticks(np.arange(len(labels)), labels, fontsize=8)
        ax.set_yticks(np.arange(len(frame)), frame["short_label"], fontsize=8)
        ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False)
        for row, record in enumerate(frame.itertuples(index=False)):
            heldout_value = record.heldout_project_direction_concordance
            ratio = record.absolute_effect_ratio_adjusted_to_unadjusted
            for column, value in enumerate(values[row]):
                text = "yes" if value else "no"
                if column == 2:
                    text = f"{heldout_value:.0%}"
                elif column == 4 and np.isfinite(ratio):
                    text = f"{ratio:.0%}"
                ax.text(
                    column,
                    row,
                    text,
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=7,
                    fontweight="bold",
                )
        ax.set_title(config.display_name, fontsize=12, fontweight="bold", pad=45)
        ax.set_xticks(np.arange(-0.5, len(labels), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(frame), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.2)
        ax.tick_params(which="minor", bottom=False, left=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle(
        "Pathway-level robustness matrix for the 29 literature-reviewed programs",
        fontsize=16,
        fontweight="bold",
    )
    fig.savefig(FIGURE_DIR / "figure_s7_pathway_robustness_matrix.png", dpi=300)
    fig.savefig(FIGURE_DIR / "figure_s7_pathway_robustness_matrix.pdf")
    plt.close(fig)


def write_summary(matrix: pd.DataFrame) -> None:
    lines = [
        "# Pathway-level robustness interpretation",
        "",
        "The labels below integrate five directional checks: ssGSEA, preranked GSEA, leave-one-project-out prediction, three full reference-query training seeds, and adjustment for atlas-derived broad composition proxies. They are descriptive evidence categories, not hypothesis-test significance levels.",
        "",
        "- **Triangulated:** all five checks support the primary direction.",
        "- **Internally robust, incomplete conventional support:** held-out, seed, and composition checks support the direction, but one or both conventional methods do not. These are the clearest expiMap-specific complementary hypotheses.",
        "- **Method-supported, model-sensitive:** conventional methods and held-out projects support the direction, but full-pipeline seed or composition sensitivity does not.",
        "- **Sensitivity-dependent:** the pathway does not meet either reproducibility pattern.",
        "",
    ]
    for config in CONFIGS:
        lines.extend([f"## {config.display_name}", ""])
        frame = matrix.loc[matrix["tissue"].eq(config.tissue)]
        for status, group in frame.groupby("robustness_status", sort=False):
            labels = ", ".join(group["short_label"])
            lines.append(f"- **{status}:** {labels}")
        lines.append("")
    (Path(SOURCE_DIR).parent / "reviewer_pathway_evidence.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    matrix = build_evidence_matrix()
    matrix.to_csv(
        SOURCE_DIR / "table_s24_pathway_robustness_evidence.tsv", sep="\t", index=False
    )
    plot_evidence_matrix(matrix)
    write_summary(matrix)
    print(
        matrix.groupby(["tissue", "robustness_status"], observed=True)
        .size()
        .rename("pathways")
        .to_string()
    )


if __name__ == "__main__":
    main()
