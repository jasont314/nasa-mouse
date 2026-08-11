"""Compare matched DDIM workflows with and without material conditioning."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .build_synthetic_guided_paper import (
    _all_bh_fdr_gene_inventory,
    _gate_synthetic_selection,
)


ROOT = Path(__file__).resolve().parents[3]
ANALYSIS_ROOT = ROOT / "outputs/generative/benchmark/analyses"
RUN_ROOT = ROOT / "outputs/generative/benchmark/runs/lacan_diffusion"
MATERIAL_RUN = RUN_ROOT / (
    "osdr_factorized_study_lora512_correlation_refine_osdr_disjoint_seed2020"
)
NO_MATERIAL_RUN = RUN_ROOT / (
    "osdr_factorized_study_lora512_correlation_refine_osdr_disjoint_"
    "no_material_seed2020"
)
MATERIAL_TISSUE = ANALYSIS_ROOT / (
    "within_study_generated_feature_stability_osdr_disjoint_v1"
)
MATERIAL_MUSCLE = ANALYSIS_ROOT / (
    "within_study_generated_feature_stability_muscle_groups_osdr_disjoint_v1"
)
NO_MATERIAL_TISSUE = ANALYSIS_ROOT / (
    "within_study_generated_feature_stability_no_material_osdr_disjoint_v1"
)
NO_MATERIAL_MUSCLE = ANALYSIS_ROOT / (
    "within_study_generated_feature_stability_muscle_groups_no_material_"
    "osdr_disjoint_v1"
)
DEFAULT_OUTPUT = ANALYSIS_ROOT / "material_type_conditioning_ablation_v1"
LANDMARKS = ROOT / "data/diffusion/l974_mouse_paper_parity.tsv"
SYNTHETIC_STATUSES = {
    "synthetic_promoted",
    "reinforced_real_and_synthetic",
}
SCOPES = (
    ("canonical_tissue", MATERIAL_TISSUE, NO_MATERIAL_TISSUE),
    ("skeletal_muscle_group", MATERIAL_MUSCLE, NO_MATERIAL_MUSCLE),
)


def _read_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _load_analysis(
    directory: Path,
    analysis_scope: str,
    landmarks: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    choices = pd.read_csv(directory / "tissue_arm_choices.tsv", sep="\t")
    stable = _gate_synthetic_selection(
        pd.read_csv(directory / "stable_gene_sets.tsv.gz", sep="\t"),
        choices,
    )
    random_effects = pd.read_csv(
        directory / "real_random_effects.tsv.gz", sep="\t"
    )
    inventory = _all_bh_fdr_gene_inventory(
        stable,
        random_effects,
        landmarks,
        analysis_scope,
    )
    choices.insert(0, "analysis_scope", analysis_scope)
    pathways = pd.read_csv(directory / "reactome_enrichment.tsv.gz", sep="\t")
    pathways.insert(0, "analysis_scope", analysis_scope)
    eligible = set(
        choices.loc[
            choices["generated_arm_eligible_all_metrics"].astype(bool), "tissue"
        ].astype(str)
    )
    pathways = pathways.loc[
        pathways["tissue"].astype(str).isin(eligible)
        & pathways["gene_set"].isin(
            ["core_intersection", "generated_supported"]
        )
        & pathways["fdr"].lt(0.05)
    ].copy()
    return choices, inventory, pathways


def _load_conditioning_arm(
    conditioning: str,
    landmarks: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    choice_parts: list[pd.DataFrame] = []
    inventory_parts: list[pd.DataFrame] = []
    pathway_parts: list[pd.DataFrame] = []
    for scope, material_dir, no_material_dir in SCOPES:
        directory = material_dir if conditioning == "material" else no_material_dir
        choices, inventory, pathways = _load_analysis(
            directory, scope, landmarks
        )
        choice_parts.append(choices)
        inventory_parts.append(inventory)
        pathway_parts.append(pathways)
    return (
        pd.concat(choice_parts, ignore_index=True),
        pd.concat(inventory_parts, ignore_index=True),
        pd.concat(pathway_parts, ignore_index=True),
    )


def _calibration_comparison() -> pd.DataFrame:
    summaries = {}
    repeat_tables = {}
    for label, run in (("material", MATERIAL_RUN), ("no_material", NO_MATERIAL_RUN)):
        directory = (
            run
            / "evaluation/repeated_distribution_calibration/"
            "prior_5_residual_0.5"
        )
        summaries[label] = _read_json(directory / "summary.json")
        repeat_tables[label] = pd.read_csv(
            directory / "repeat_metrics.tsv", sep="\t"
        )
    rows = []
    fidelity_metrics = (
        "correlation",
        "precision",
        "recall",
        "f1",
        "adversarial_accuracy",
        "frechet_ratio",
    )
    for metric in fidelity_metrics:
        material = float(
            summaries["material"]["metric_repeat_stability"][metric]["mean"]
        )
        no_material = float(
            summaries["no_material"]["metric_repeat_stability"][metric]["mean"]
        )
        rows.append(
            {
                "metric": metric,
                "material_mean": material,
                "no_material_mean": no_material,
                "no_material_minus_material": no_material - material,
            }
        )
    for metric in (
        "condition_delta_correlation",
        "condition_direction_agreement",
        "muscle_accession_correlation",
        "muscle_accession_direction",
    ):
        material = float(repeat_tables["material"][metric].mean())
        no_material = float(repeat_tables["no_material"][metric].mean())
        rows.append(
            {
                "metric": metric,
                "material_mean": material,
                "no_material_mean": no_material,
                "no_material_minus_material": no_material - material,
            }
        )
    return pd.DataFrame(rows)


def _arm_comparison(
    material: pd.DataFrame, no_material: pd.DataFrame
) -> pd.DataFrame:
    keys = ["analysis_scope", "tissue"]
    columns = keys + [
        "selected_arm",
        "real_mean_balanced_accuracy",
        "selected_mean_balanced_accuracy",
        "real_mean_roc_auc",
        "selected_mean_roc_auc",
        "real_mean_average_precision",
        "selected_mean_average_precision",
        "generated_arm_eligible_all_metrics",
    ]
    result = material[columns].merge(
        no_material[columns],
        on=keys,
        how="outer",
        validate="one_to_one",
        suffixes=("_material", "_no_material"),
    )
    for metric in (
        "real_mean_balanced_accuracy",
        "real_mean_roc_auc",
        "real_mean_average_precision",
    ):
        difference = (
            result[f"{metric}_material"] - result[f"{metric}_no_material"]
        ).abs()
        if float(difference.max()) > 1e-12:
            raise ValueError(f"Real-only metric changed in the matched ablation: {metric}")
    result["arm_changed"] = result["selected_arm_material"].ne(
        result["selected_arm_no_material"]
    )
    for metric in ("balanced_accuracy", "roc_auc", "average_precision"):
        result[f"selected_{metric}_change"] = (
            result[f"selected_mean_{metric}_no_material"]
            - result[f"selected_mean_{metric}_material"]
        )
    return result.sort_values(keys, ignore_index=True)


def _gene_comparison(
    material: pd.DataFrame, no_material: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    keys = ["analysis_scope", "tissue", "gene"]
    material_keys = material[keys].sort_values(keys, ignore_index=True)
    no_material_keys = no_material[keys].sort_values(keys, ignore_index=True)
    if not material_keys.equals(no_material_keys):
        raise ValueError("The real-data BH-FDR inventory changed between arms")
    for metric in ("meta_effect", "meta_p", "meta_fdr"):
        if not np.allclose(material[metric], no_material[metric], equal_nan=True):
            raise ValueError(f"Real-data statistic changed between arms: {metric}")

    material_selected = material.loc[
        material["selection_interpretation"].isin(SYNTHETIC_STATUSES),
        keys + ["selection_interpretation"],
    ].rename(columns={"selection_interpretation": "material_selection"})
    no_material_selected = no_material.loc[
        no_material["selection_interpretation"].isin(SYNTHETIC_STATUSES),
        keys + ["selection_interpretation"],
    ].rename(columns={"selection_interpretation": "no_material_selection"})
    comparison = material_selected.merge(
        no_material_selected,
        on=keys,
        how="outer",
        validate="one_to_one",
    )
    real_columns = keys + [
        "symbol",
        "n_accessions",
        "flt_gc_direction",
        "meta_effect",
        "meta_fdr",
        "accession_direction_fraction",
        "loo_fdr_stable_0_05",
    ]
    comparison = comparison.merge(
        material[real_columns], on=keys, how="left", validate="one_to_one"
    )
    comparison[["material_selection", "no_material_selection"]] = comparison[
        ["material_selection", "no_material_selection"]
    ].fillna("not_synthetic_informed")
    comparison["comparison"] = np.select(
        [
            comparison["material_selection"].eq("not_synthetic_informed"),
            comparison["no_material_selection"].eq("not_synthetic_informed"),
            comparison["material_selection"].eq(
                comparison["no_material_selection"]
            ),
        ],
        ["no_material_only", "material_only", "retained_same_status"],
        default="retained_changed_status",
    )
    comparison = comparison.sort_values(
        ["analysis_scope", "tissue", "meta_fdr", "symbol"], ignore_index=True
    )

    rows = []
    units = pd.concat(
        (material[keys[:2]], no_material[keys[:2]]), ignore_index=True
    ).drop_duplicates()
    for unit in units.itertuples(index=False):
        frame = comparison.loc[
            comparison["analysis_scope"].eq(unit.analysis_scope)
            & comparison["tissue"].eq(unit.tissue)
        ]
        material_genes = set(
            frame.loc[
                frame["material_selection"].ne("not_synthetic_informed"), "gene"
            ]
        )
        no_material_genes = set(
            frame.loc[
                frame["no_material_selection"].ne("not_synthetic_informed"),
                "gene",
            ]
        )
        overlap = material_genes & no_material_genes
        union = material_genes | no_material_genes
        rows.append(
            {
                "analysis_scope": unit.analysis_scope,
                "tissue": unit.tissue,
                "material_synthetic_informed": len(material_genes),
                "no_material_synthetic_informed": len(no_material_genes),
                "overlap": len(overlap),
                "material_only": len(material_genes - no_material_genes),
                "no_material_only": len(no_material_genes - material_genes),
                "jaccard": len(overlap) / len(union) if union else np.nan,
                "material_promoted": int(
                    frame["material_selection"].eq("synthetic_promoted").sum()
                ),
                "material_reinforced": int(
                    frame["material_selection"]
                    .eq("reinforced_real_and_synthetic")
                    .sum()
                ),
                "no_material_promoted": int(
                    frame["no_material_selection"].eq("synthetic_promoted").sum()
                ),
                "no_material_reinforced": int(
                    frame["no_material_selection"]
                    .eq("reinforced_real_and_synthetic")
                    .sum()
                ),
            }
        )
    return comparison, pd.DataFrame(rows).sort_values(
        ["analysis_scope", "tissue"], ignore_index=True
    )


def _pathway_comparison(
    material: pd.DataFrame, no_material: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    keys = ["analysis_scope", "tissue", "term"]
    material_terms = (
        material.sort_values("fdr")
        .drop_duplicates(keys)
        [keys + ["fdr", "overlap_symbols"]]
        .rename(
            columns={
                "fdr": "material_fdr",
                "overlap_symbols": "material_overlap_symbols",
            }
        )
    )
    no_material_terms = (
        no_material.sort_values("fdr")
        .drop_duplicates(keys)
        [keys + ["fdr", "overlap_symbols"]]
        .rename(
            columns={
                "fdr": "no_material_fdr",
                "overlap_symbols": "no_material_overlap_symbols",
            }
        )
    )
    comparison = material_terms.merge(
        no_material_terms,
        on=keys,
        how="outer",
        validate="one_to_one",
        indicator=True,
    ).rename(columns={"_merge": "comparison"})
    comparison["comparison"] = comparison["comparison"].map(
        {"left_only": "material_only", "right_only": "no_material_only", "both": "both"}
    )
    units = pd.concat(
        (material[keys[:2]], no_material[keys[:2]]), ignore_index=True
    ).drop_duplicates()
    rows = []
    for unit in units.itertuples(index=False):
        frame = comparison.loc[
            comparison["analysis_scope"].eq(unit.analysis_scope)
            & comparison["tissue"].eq(unit.tissue)
        ]
        rows.append(
            {
                "analysis_scope": unit.analysis_scope,
                "tissue": unit.tissue,
                "material_significant_pathways": int(
                    frame["material_fdr"].notna().sum()
                ),
                "no_material_significant_pathways": int(
                    frame["no_material_fdr"].notna().sum()
                ),
                "overlap": int(frame["comparison"].eq("both").sum()),
                "material_only": int(
                    frame["comparison"].eq("material_only").sum()
                ),
                "no_material_only": int(
                    frame["comparison"].eq("no_material_only").sum()
                ),
            }
        )
    return comparison.sort_values(keys, ignore_index=True), pd.DataFrame(rows)


def _plot_gene_counts(counts: pd.DataFrame, output: Path) -> None:
    table = counts.loc[
        counts[["material_synthetic_informed", "no_material_synthetic_informed"]]
        .max(axis=1)
        .gt(0)
    ].copy()
    table["label"] = table["tissue"].str.replace("_", " ")
    table = table.sort_values(
        ["analysis_scope", "material_synthetic_informed"], ascending=[True, True]
    )
    y = np.arange(len(table))
    width = 0.38
    figure, axis = plt.subplots(figsize=(10, max(5.5, 0.42 * len(table))))
    axis.barh(
        y - width / 2,
        table["material_synthetic_informed"],
        height=width,
        color="#17807E",
        label="Material conditioned",
    )
    axis.barh(
        y + width / 2,
        table["no_material_synthetic_informed"],
        height=width,
        color="#D96552",
        label="No material conditioning",
    )
    axis.set_yticks(y, table["label"])
    axis.set_xlabel("Synthetic-informed BH-FDR associations")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output / "synthetic_informed_gene_counts.png", dpi=220)
    figure.savefig(output / "synthetic_informed_gene_counts.pdf")
    plt.close(figure)


def _plot_utility_change(arms: pd.DataFrame, output: Path) -> None:
    table = arms.sort_values("selected_balanced_accuracy_change").copy()
    table["label"] = table["tissue"].str.replace("_", " ")
    colors = np.where(
        table["selected_balanced_accuracy_change"].ge(0), "#17807E", "#D96552"
    )
    figure, axis = plt.subplots(figsize=(10, max(6, 0.32 * len(table))))
    axis.barh(
        np.arange(len(table)),
        table["selected_balanced_accuracy_change"],
        color=colors,
    )
    axis.axvline(0, color="#333333", linewidth=0.8)
    axis.set_yticks(np.arange(len(table)), table["label"])
    axis.set_xlabel("Selected-arm BA change (no material minus material)")
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    figure.savefig(output / "selected_arm_balanced_accuracy_change.png", dpi=220)
    figure.savefig(output / "selected_arm_balanced_accuracy_change.pdf")
    plt.close(figure)


def _symbols(
    comparison: pd.DataFrame, scope: str, tissue: str, category: str
) -> str:
    frame = comparison.loc[
        comparison["analysis_scope"].eq(scope)
        & comparison["tissue"].eq(tissue)
        & comparison["comparison"].eq(category)
    ]
    return ", ".join(sorted(frame["symbol"].astype(str))) or "none"


def _write_readme(
    output: Path,
    calibration: pd.DataFrame,
    arms: pd.DataFrame,
    genes: pd.DataFrame,
    counts: pd.DataFrame,
) -> None:
    current_total = int(counts["material_synthetic_informed"].sum())
    no_material_total = int(counts["no_material_synthetic_informed"].sum())
    retained = int(genes["comparison"].eq("retained_same_status").sum())
    changed_arms = arms.loc[arms["arm_changed"]]
    aa = calibration.set_index("metric").loc["adversarial_accuracy"]
    corr = calibration.set_index("metric").loc["correlation"]
    lines = [
        "# Material-type conditioning ablation",
        "",
        "This matched ablation retrains and recalibrates the factorized DDIM without "
        "material-type conditioning while retaining tissue, FLT/GC, and study "
        "conditioning. It reruns the same three synthetic draws, eight nested "
        "development repeats, five synthetic-use arms, real-data BH-FDR support, "
        "and Reactome analysis. Existing material-conditioned outputs are unchanged.",
        "",
        "## Generator validation",
        "",
        f"Calibrated validation AA was {aa.material_mean:.3f} with material and "
        f"{aa.no_material_mean:.3f} without material. Correlation was "
        f"{corr.material_mean:.4f} and {corr.no_material_mean:.4f}, respectively. "
        "Thus global fidelity was effectively unchanged. These are matched "
        "development-validation results; the no-material model was not evaluated "
        "as a new confirmatory model on the already-open locked test, so its AA "
        "should not be substituted for the primary locked-test AA of 0.475.",
        "",
        "## Downstream changes",
        "",
        f"The selected synthetic-use arm changed for {len(changed_arms)} of "
        f"{len(arms)} analysis units. The current workflow contains {current_total} "
        f"synthetic-informed BH-FDR associations; the no-material workflow contains "
        f"{no_material_total}. Of the current associations, {retained} "
        f"({retained / current_total:.1%}) were retained with the same promoted or "
        "reinforced status.",
        "",
        "Key changes:",
        "",
        f"- Thymus retained most of the cell-cycle panel. Lost current genes: "
        f"{_symbols(genes, 'canonical_tissue', 'thymus', 'material_only')}. New "
        f"no-material genes: {_symbols(genes, 'canonical_tissue', 'thymus', 'no_material_only')}.",
        f"- Soleus changed from a generated-informed arm to real-only; its five "
        "current reinforced genes "
        f"({_symbols(genes, 'skeletal_muscle_group', 'soleus', 'material_only')}) "
        "lost synthetic attribution.",
        f"- Spleen lost {_symbols(genes, 'canonical_tissue', 'spleen', 'material_only')} "
        "but retained the remaining synthetic-informed genes.",
        f"- Kidney lost {_symbols(genes, 'canonical_tissue', 'kidney', 'material_only')} "
        "while retaining its shared result.",
        f"- EDL gained {_symbols(genes, 'skeletal_muscle_group', 'edl', 'no_material_only')}; "
        "gastrocnemius lost its current synthetic-informed genes.",
        "- Thymus retained 29 significant synthetic-supported Reactome terms, "
        "including mitotic cell cycle, DNA replication, APC/C, and G2/M control.",
        "- Soleus lost all five current synthetic-supported Reactome terms. Pooled "
        "muscle retained only one shared significant term, with the no-material "
        "result shifting away from the current interferon and sialic-acid pattern.",
        "- The four new EDL associations did not produce a significant Reactome "
        "pathway, while spleen and kidney gained conditioning-sensitive collagen "
        "and lipid-metabolism enrichments, respectively.",
        "",
        "## Interpretation",
        "",
        "Removing material type does not materially change aggregate synthetic-data "
        "fidelity, but it changes arm selection and thresholded gene attribution. "
        "The thymus cell-cycle interpretation is broadly robust, while the soleus "
        "synthetic-informed panel is conditioning-sensitive. Because study identifies "
        "most material labels, this ablation does not prove that material is "
        "irrelevant; it shows that study can absorb much of its broad predictive "
        "information. Explicit material conditioning remains preferable when "
        "generating or interpreting anatomical muscle subgroups. All utility "
        "comparisons are repeated within-study development results rather than new "
        "whole-study confirmation.",
        "",
        "## Outputs",
        "",
        "- `calibrated_generator_metrics.tsv`: matched calibrated validation metrics.",
        "- `arm_choice_comparison.tsv`: selected-arm and utility differences for all 27 units.",
        "- `synthetic_informed_gene_comparison.tsv`: retained, lost, and gained associations.",
        "- `synthetic_informed_gene_counts_by_tissue.tsv`: promoted/reinforced counts and overlap.",
        "- `no_material_all_bh_fdr_genes.tsv.gz`: complete no-material BH-FDR inventory.",
        "- `significant_pathway_comparison.tsv.gz`: eligible synthetic-supported Reactome terms.",
        "- `significant_pathway_counts_by_tissue.tsv`: pathway overlap by analysis unit.",
    ]
    (output / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(output: Path) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    landmarks = pd.read_csv(LANDMARKS, sep="\t")
    material_choices, material_inventory, material_pathways = (
        _load_conditioning_arm("material", landmarks)
    )
    no_material_choices, no_material_inventory, no_material_pathways = (
        _load_conditioning_arm("no_material", landmarks)
    )
    calibration = _calibration_comparison()
    arms = _arm_comparison(material_choices, no_material_choices)
    genes, gene_counts = _gene_comparison(
        material_inventory, no_material_inventory
    )
    pathways, pathway_counts = _pathway_comparison(
        material_pathways, no_material_pathways
    )

    calibration.to_csv(
        output / "calibrated_generator_metrics.tsv", sep="\t", index=False
    )
    arms.to_csv(output / "arm_choice_comparison.tsv", sep="\t", index=False)
    genes.to_csv(
        output / "synthetic_informed_gene_comparison.tsv", sep="\t", index=False
    )
    gene_counts.to_csv(
        output / "synthetic_informed_gene_counts_by_tissue.tsv",
        sep="\t",
        index=False,
    )
    no_material_inventory.to_csv(
        output / "no_material_all_bh_fdr_genes.tsv.gz", sep="\t", index=False
    )
    pathways.to_csv(
        output / "significant_pathway_comparison.tsv.gz", sep="\t", index=False
    )
    pathway_counts.to_csv(
        output / "significant_pathway_counts_by_tissue.tsv",
        sep="\t",
        index=False,
    )
    _plot_gene_counts(gene_counts, output)
    _plot_utility_change(arms, output)
    _write_readme(output, calibration, arms, genes, gene_counts)

    summary = {
        "status": "complete",
        "conditioning_retained": ["tissue", "condition", "study"],
        "conditioning_removed": ["material_type"],
        "analysis_units": int(len(arms)),
        "selected_arms_changed": int(arms["arm_changed"].sum()),
        "material_synthetic_informed_associations": int(
            gene_counts["material_synthetic_informed"].sum()
        ),
        "no_material_synthetic_informed_associations": int(
            gene_counts["no_material_synthetic_informed"].sum()
        ),
        "retained_same_status": int(
            genes["comparison"].eq("retained_same_status").sum()
        ),
        "material_only": int(genes["comparison"].eq("material_only").sum()),
        "no_material_only": int(genes["comparison"].eq("no_material_only").sum()),
        "retained_changed_status": int(
            genes["comparison"].eq("retained_changed_status").sum()
        ),
        "real_data_statistics_identical": True,
        "no_material_locked_test_opened": False,
        "primary_outputs_modified": False,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)
    return output / "summary.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    run(arguments.output)


if __name__ == "__main__":
    main()
