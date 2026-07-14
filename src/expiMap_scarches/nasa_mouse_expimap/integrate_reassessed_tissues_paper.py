"""Integrate corrected kidney and spleen models into the ASGSR paper package."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import textwrap

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

from .build_asgsr_paper import ROOT


PAPER_DIR = ROOT / "paper/asgsr_expimap_hvg"
SOURCE_DIR = PAPER_DIR / "source_data"
FIGURE_DIR = PAPER_DIR / "figures"
REASSESSMENT_DIR = ROOT / "outputs/expimap_kidney_spleen_reassessment"

MAIN_TISSUES = ("thymus", "skin", "liver", "spleen")
TISSUE_ORDER = (*MAIN_TISSUES, "kidney")
ROLE_COLORS = {
    "aligned": "#24834b",
    "complementary": "#2869a8",
    "context_sensitive": "#b36a19",
}
TISSUE_COLORS = {
    "thymus": "#7855a6",
    "skin": "#bd5b2a",
    "liver": "#1a7a73",
    "spleen": "#9b4f65",
    "kidney": "#287b8e",
}

RETAINED_TERMS = {
    "thymus": (
        "R-MMU-73894_DNA_REPAIR",
        "R-MMU-8980692_RHOA_GTPASE_CYCLE",
        "R-MMU-198933_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
    ),
    "skin": (
        "R-MMU-3247509_CHROMATIN_MODIFYING_ENZYMES",
        "R-MMU-73894_DNA_REPAIR",
        "R-MMU-5358351_SIGNALING_BY_HEDGEHOG",
        "R-MMU-428157_SPHINGOLIPID_METABOLISM",
        "R-MMU-421270_CELL_CELL_JUNCTION_ORGANIZATION",
    ),
    "liver": (
        "R-MMU-2132295_MHC_CLASS_II_ANTIGEN_PRESENTATION",
        "R-MMU-202403_TCR_SIGNALING",
    ),
    "spleen": (
        "R-MMU-202403_TCR_SIGNALING",
        "R-MMU-6798695_NEUTROPHIL_DEGRANULATION",
        "R-MMU-5621481_C_TYPE_LECTIN_RECEPTORS_CLRS",
    ),
    "kidney": (
        "R-MMU-3000178_ECM_PROTEOGLYCANS",
        "R-MMU-195721_SIGNALING_BY_WNT",
        "R-MMU-381426_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS",
    ),
}

DISPLAY_LABELS = {
    "R-MMU-73894_DNA_REPAIR": "DNA repair",
    "R-MMU-8980692_RHOA_GTPASE_CYCLE": "RHOA cytoskeletal cycle",
    "R-MMU-198933_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL": "Lymphoid-stromal interactions",
    "R-MMU-3247509_CHROMATIN_MODIFYING_ENZYMES": "Chromatin-modifying enzymes",
    "R-MMU-5358351_SIGNALING_BY_HEDGEHOG": "Hedgehog signaling",
    "R-MMU-428157_SPHINGOLIPID_METABOLISM": "Sphingolipid metabolism",
    "R-MMU-421270_CELL_CELL_JUNCTION_ORGANIZATION": "Cell-cell junction organization",
    "R-MMU-2132295_MHC_CLASS_II_ANTIGEN_PRESENTATION": "MHC class II antigen presentation",
    "R-MMU-202403_TCR_SIGNALING": "T-cell receptor signaling",
    "R-MMU-6798695_NEUTROPHIL_DEGRANULATION": "Neutrophil degranulation program",
    "R-MMU-5621481_C_TYPE_LECTIN_RECEPTORS_CLRS": "C-type lectin receptor signaling",
    "R-MMU-3000178_ECM_PROTEOGLYCANS": "ECM proteoglycans",
    "R-MMU-195721_SIGNALING_BY_WNT": "WNT signaling",
    "R-MMU-381426_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS": "IGF transport and uptake",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_figure(fig: plt.Figure, name: str) -> None:
    fig.savefig(FIGURE_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{name}.pdf", bbox_inches="tight")


def build_model_scope() -> pd.DataFrame:
    original = pd.read_csv(SOURCE_DIR / "table_s1_model_summary.tsv", sep="\t")
    original["mapped_query_samples"] = original["query_samples"]
    original["primary_effect_samples"] = original["query_samples"]
    original["primary_effect_flight"] = original["query_flight"]
    original["primary_effect_ground_control"] = original["query_ground_control"]
    original["primary_effect_projects"] = original["query_accessions"]
    original["analysis_role"] = original["tissue"].map(
        {
            "thymus": "main",
            "skin": "main",
            "liver": "main",
            "soleus": "not_advanced_supplementary_record",
        }
    )

    new_rows = []
    seed_manifest = pd.read_csv(
        REASSESSMENT_DIR / "seed_training_manifest.tsv", sep="\t"
    )
    for tissue in ("kidney", "spleen"):
        base = ROOT / f"outputs/expimap_archs4_reference_osdr_query_{tissue}/reassessment_hvg_2000"
        manifest = load_json(base / "input/tutorial_hvg_input_manifest.json")
        query_path = base / "input" / f"osdr_{tissue}_query_tutorial_hvg_raw_counts.h5ad"
        obs = ad.read_h5ad(query_path, backed="r").obs.copy()
        retained = pd.Series(True, index=obs.index)
        excluded = ""
        if tissue == "spleen":
            retained = ~obs["id.accession"].astype(str).eq("OSD-288")
            excluded = "OSD-288"
        primary = obs.loc[retained]
        training = seed_manifest.loc[
            seed_manifest["tissue"].eq(tissue) & seed_manifest["seed"].eq(2020)
        ].iloc[0]
        new_rows.append(
            {
                "tissue": tissue,
                "reference_samples": manifest["n_reference_samples"],
                "reference_series": int(training["reference_series"]),
                "query_samples": manifest["n_query_samples"],
                "query_flight": int(
                    obs["condition_inferred"].astype(str).eq("flight").sum()
                ),
                "query_ground_control": int(
                    obs["condition_inferred"].astype(str).eq("ground_control").sum()
                ),
                "query_accessions": int(obs["id.accession"].nunique()),
                "source_query_samples_before_primary_filter": manifest[
                    "n_query_samples"
                ],
                "source_query_accessions_before_primary_filter": int(
                    obs["id.accession"].nunique()
                ),
                "primary_excluded_accessions": excluded,
                "hvg_requested": manifest["n_top_genes_requested"],
                "genes_after_filter": manifest["n_genes_after_term_filter"],
                "reactome_programs": manifest["n_terms_after_hvg_filter"],
                "reference_epochs_requested": 400,
                "reference_epochs_completed": int(
                    training["reference_epochs_completed"]
                ),
                "query_epochs": int(training["query_epochs"]),
                "reconstruction_loss": "nb",
                "hidden_layers": "300x300x300",
                "gpu": training["gpu"],
                "posterior_mean_scores": True,
                "mapped_query_samples": manifest["n_query_samples"],
                "primary_effect_samples": int(len(primary)),
                "primary_effect_flight": int(
                    primary["condition_inferred"].astype(str).eq("flight").sum()
                ),
                "primary_effect_ground_control": int(
                    primary["condition_inferred"]
                    .astype(str)
                    .eq("ground_control")
                    .sum()
                ),
                "primary_effect_projects": int(primary["id.accession"].nunique()),
                "analysis_role": "main" if tissue == "spleen" else "secondary_exploratory",
            }
        )
    combined = pd.concat([original, pd.DataFrame(new_rows)], ignore_index=True)
    order = ["thymus", "skin", "liver", "spleen", "kidney", "soleus"]
    combined["display_order"] = combined["tissue"].map(
        {tissue: index for index, tissue in enumerate(order)}
    )
    return combined.sort_values("display_order").drop(columns="display_order")


def load_retained_evidence() -> pd.DataFrame:
    original = pd.read_csv(
        SOURCE_DIR / "table_s24_pathway_robustness_evidence.tsv", sep="\t"
    )
    curated = pd.read_csv(
        SOURCE_DIR / "table_1_curated_pathway_results.tsv", sep="\t"
    )[
        [
            "tissue",
            "term",
            "paper_interpretation",
            "paper_citations",
            "within_tissue_magnitude_percentile",
        ]
    ]
    original = original.merge(curated, on=["tissue", "term"], how="left")
    original = original.loc[
        [
            row.term in RETAINED_TERMS.get(row.tissue, ())
            for row in original.itertuples(index=False)
        ]
    ].copy()
    original["analysis_role"] = "main"
    original["display_label"] = original["term"].map(DISPLAY_LABELS)
    original.loc[
        original["term"].eq("R-MMU-421270_CELL_CELL_JUNCTION_ORGANIZATION"),
        "evidence_role",
    ] = "complementary"
    original["seed_effect_median"] = original[
        ["effect_seed2020", "effect_seed2021", "effect_seed2022"]
    ].median(axis=1)
    original["seed_effect_minimum"] = original[
        ["effect_seed2020", "effect_seed2021", "effect_seed2022"]
    ].min(axis=1)
    original["seed_effect_maximum"] = original[
        ["effect_seed2020", "effect_seed2021", "effect_seed2022"]
    ].max(axis=1)
    original["manual_rationale"] = original["paper_interpretation"]
    original["literature_keys"] = original["paper_citations"]

    reassessed = pd.read_csv(
        REASSESSMENT_DIR / "top_decile_manual_review.tsv", sep="\t"
    )
    reassessed = reassessed.loc[
        [
            row.term in RETAINED_TERMS.get(row.tissue, ())
            for row in reassessed.itertuples(index=False)
        ]
    ].copy()
    reassessed["analysis_role"] = np.where(
        reassessed["tissue"].eq("spleen"), "main", "secondary_exploratory"
    )
    reassessed["display_label"] = reassessed["term"].map(DISPLAY_LABELS)
    reassessed["short_label"] = reassessed["display_label"]
    reassessed["evidence_role"] = np.where(
        reassessed["manual_category"].eq("coherent_aligned"),
        "aligned",
        "complementary",
    )
    reassessed["within_tissue_magnitude_percentile"] = reassessed[
        "primary_absolute_percentile"
    ]

    columns = sorted(set(original.columns).union(reassessed.columns))
    aligned = [original.reindex(columns=columns), reassessed.reindex(columns=columns)]
    boolean_columns = {
        column
        for frame in (original, reassessed)
        for column in frame.columns
        if pd.api.types.is_bool_dtype(frame[column])
    }
    for frame in aligned:
        for column in boolean_columns:
            frame[column] = frame[column].astype("boolean")
    combined = pd.concat(aligned, ignore_index=True)
    tissue_rank = {tissue: index for index, tissue in enumerate(TISSUE_ORDER)}
    term_rank = {
        (tissue, term): index
        for tissue, terms in RETAINED_TERMS.items()
        for index, term in enumerate(terms)
    }
    combined["tissue_order"] = combined["tissue"].map(tissue_rank)
    combined["term_order"] = [
        term_rank[(tissue, term)]
        for tissue, term in zip(combined["tissue"], combined["term"])
    ]
    return combined.sort_values(["tissue_order", "term_order"]).drop(
        columns=["tissue_order", "term_order"]
    )


def write_retained_table(evidence: pd.DataFrame) -> None:
    columns = [
        "tissue",
        "analysis_role",
        "term",
        "display_label",
        "evidence_role",
        "effect_seed2020",
        "effect_seed2021",
        "effect_seed2022",
        "seed_effect_median",
        "seed_effect_minimum",
        "seed_effect_maximum",
        "expimap_n_projects",
        "expimap_projects_positive",
        "expimap_projects_negative",
        "ssgsea_project_balanced_effect",
        "gsea_nes",
        "gsea_fdr",
        "heldout_project_direction_concordance",
        "absolute_effect_ratio_adjusted_to_unadjusted",
        "robustness_support_count",
        "robustness_status",
        "manual_rationale",
        "literature_keys",
    ]
    evidence.reindex(columns=columns).to_csv(
        SOURCE_DIR / "table_2_retained_pathway_evidence.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )


def copy_reassessment_tables() -> None:
    sources = {
        "seed_training_manifest.tsv": "table_s26_kidney_spleen_training_manifest.tsv",
        "pathway_evidence_matrix.tsv": "table_s27_kidney_spleen_pathway_evidence.tsv",
        "top_decile_manual_review.tsv": "table_s28_kidney_spleen_manual_review.tsv",
        "top_decile_member_gene_support.tsv": "table_s29_kidney_spleen_member_gene_support.tsv",
        "literature_sources.tsv": "table_s30_kidney_spleen_literature_sources.tsv",
    }
    for source, destination in sources.items():
        shutil.copy2(REASSESSMENT_DIR / source, SOURCE_DIR / destination)


def copy_reassigned_figures() -> None:
    pairs = {
        "figure_7_skin_protocol_context": "figure_6_skin_protocol_context",
        "figure_6_generated_biological_processes": "figure_s8_generated_biological_processes",
        "figure_4_primary_analysis_sensitivity": "figure_s9_original_tissue_sensitivity",
    }
    for source, destination in pairs.items():
        for suffix in ("png", "pdf"):
            shutil.copy2(
                FIGURE_DIR / f"{source}.{suffix}",
                FIGURE_DIR / f"{destination}.{suffix}",
            )
    for suffix in ("png", "pdf"):
        shutil.copy2(
            REASSESSMENT_DIR
            / f"curated_main_pathway_project_seed_evidence.{suffix}",
            FIGURE_DIR / f"figure_4_kidney_spleen_reassessment.{suffix}",
        )


def plot_workflow(scope: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(13.2, 4.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = (
        (0.02, 0.20, 0.18, 0.62, "ARCHS4 references", "Six tissue-matched mouse\nbulk RNA-seq references\nSpaceflight records excluded"),
        (0.23, 0.20, 0.18, 0.62, "Reactome architecture", "Current mouse Ensembl\npathways\n2,000 reference-selected HVGs"),
        (0.44, 0.20, 0.18, 0.62, "expiMap references", "Negative-binomial models\n3 x 300 hidden units\nA100 GPU"),
        (0.65, 0.20, 0.16, 0.62, "OSDR query mapping", "Flight and ground control\nAccession-conditioned mapping\nDecoder-oriented scores"),
        (0.84, 0.20, 0.14, 0.62, "Evidence review", "Enrichment benchmarks\nHeld-out projects\nThree full seeds\nComposition proxies"),
    )
    fills = ("#e7f0f5", "#edf3e8", "#e9e7f3", "#f6ece4", "#eef0f1")
    for index, (x, y, width, height, title, body) in enumerate(boxes):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle="round,pad=0.012,rounding_size=0.015",
                linewidth=1.1,
                edgecolor="#37434a",
                facecolor=fills[index],
            )
        )
        ax.text(
            x + width / 2,
            y + height - 0.12,
            title,
            ha="center",
            va="center",
            fontsize=10.6,
            weight="bold",
        )
        ax.text(
            x + width / 2,
            y + height / 2 - 0.03,
            body,
            ha="center",
            va="center",
            fontsize=7.6,
            linespacing=1.4,
        )
        if index < len(boxes) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x + width + 0.005, 0.51),
                    (boxes[index + 1][0] - 0.005, 0.51),
                    arrowstyle="-|>",
                    mutation_scale=14,
                    linewidth=1.2,
                    color="#4b565c",
                )
            )
    scope_lookup = scope.set_index("tissue")
    main_text = ", ".join(
        f"{tissue.title()} {int(scope_lookup.loc[tissue, 'primary_effect_samples'])}"
        for tissue in MAIN_TISSUES
    )
    ax.text(
        0.5,
        0.09,
        f"Main positive tissues: {main_text}  |  Secondary: Kidney 135  |  Non-advanced screens: supplement",
        ha="center",
        va="center",
        fontsize=8.3,
        color="#30383c",
    )
    ax.set_title(
        "Reference-guided pathway analysis of NASA OSDR mouse transcriptomes",
        fontsize=14,
        weight="bold",
        pad=12,
    )
    save_figure(fig, "figure_1_workflow")
    plt.close(fig)


def project_effects(tissue: str, terms: tuple[str, ...]) -> pd.DataFrame:
    if tissue == "spleen":
        effects = pd.read_csv(
            REASSESSMENT_DIR / "seed_accession_effects.tsv.gz", sep="\t"
        )
        effects = effects.loc[
            effects["tissue"].eq(tissue)
            & effects["seed"].eq(2020)
            & ~effects["accession"].astype(str).eq("OSD-288")
            & effects["term"].isin(terms)
        ]
        return (
            effects.groupby(["project", "term"], as_index=False)["effect"]
            .mean()
            .rename(columns={"project": "heldout_project", "effect": "project_effect"})
        )
    heldout = pd.read_csv(
        SOURCE_DIR / "table_s15_project_heldout_predictions.tsv.gz", sep="\t"
    )
    heldout = heldout.loc[
        heldout["tissue"].eq(tissue)
        & heldout["method"].eq("expimap")
        & heldout["term"].isin(terms)
    ]
    return heldout[
        ["heldout_project", "term", "heldout_project_effect"]
    ].rename(columns={"heldout_project_effect": "project_effect"})


def plot_main_pathway_shifts(evidence: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 10.2), constrained_layout=True)
    for ax, tissue in zip(axes.flat, MAIN_TISSUES):
        terms = RETAINED_TERMS[tissue]
        subset = evidence.loc[evidence["tissue"].eq(tissue)].set_index("term").loc[
            list(terms)
        ].reset_index()
        points = project_effects(tissue, terms)
        positions = np.arange(len(subset))[::-1]
        for position, row in zip(positions, subset.itertuples(index=False)):
            local = points.loc[points["term"].eq(row.term), "project_effect"]
            ax.scatter(
                local,
                np.full(len(local), position),
                s=42,
                color="#7d878c",
                edgecolor="white",
                linewidth=0.6,
                alpha=0.85,
                zorder=2,
            )
            ax.hlines(
                position,
                row.seed_effect_minimum,
                row.seed_effect_maximum,
                color=ROLE_COLORS[row.evidence_role],
                linewidth=3.2,
                zorder=3,
            )
            ax.scatter(
                row.seed_effect_median,
                position,
                marker="D",
                s=86,
                color=ROLE_COLORS[row.evidence_role],
                edgecolor="white",
                linewidth=0.8,
                zorder=4,
            )
        ax.axvline(0, color="#30383c", linewidth=0.9)
        ax.set_yticks(positions)
        ax.set_yticklabels(subset["display_label"], fontsize=9.2)
        for tick, role in zip(ax.get_yticklabels(), subset["evidence_role"]):
            tick.set_color(ROLE_COLORS[role])
        project_count = int(subset["expimap_n_projects"].max())
        ax.set_title(
            f"{tissue.title()}  |  {project_count} primary projects",
            loc="left",
            fontsize=12,
            weight="bold",
        )
        ax.set_xlabel("Flight minus ground expiMap pathway shift", fontsize=9)
        ax.grid(axis="x", color="#e1e4e5", linewidth=0.7)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#7d878c",
            markeredgecolor="white",
            markersize=7,
            label="OSDR project, primary seed",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="#4b555a",
            linewidth=3,
            markerfacecolor="#4b555a",
            markersize=7,
            label="Median and range across three complete trainings",
        ),
        *[
            Line2D([0], [0], color=color, linewidth=4, label=label)
            for label, color in (
                ("Prior-literature aligned", ROLE_COLORS["aligned"]),
                ("Complementary process", ROLE_COLORS["complementary"]),
                ("Context-sensitive", ROLE_COLORS["context_sensitive"]),
            )
        ],
    ]
    fig.legend(
        handles=handles,
        frameon=False,
        ncol=3,
        loc="outside lower center",
        fontsize=8.5,
    )
    fig.suptitle(
        "Retained pathway directions across projects and full training runs",
        fontsize=15,
        weight="bold",
    )
    save_figure(fig, "figure_2_tissue_pathway_shifts")
    plt.close(fig)


def plot_evidence_matrix(evidence: pd.DataFrame) -> None:
    columns = (
        ("ssgsea_direction_support", "ssGSEA"),
        ("preranked_gsea_direction_support", "GSEA\ndirection"),
        ("heldout_direction_support", "Held-out\nprojects"),
        ("seed_direction_support", "Three\ntrainings"),
        ("composition_proxy_support", "Composition\nproxy"),
    )
    frame = evidence.copy()
    labels = []
    for row in frame.itertuples(index=False):
        suffix = " (secondary)" if row.tissue == "kidney" else ""
        labels.append(f"{row.tissue.title()}{suffix}: {row.display_label}")
    matrix = frame[[name for name, _ in columns]].astype(bool).to_numpy()

    fig, ax = plt.subplots(figsize=(13.8, 9.2), constrained_layout=False)
    fig.subplots_adjust(left=0.31, right=0.97, top=0.83, bottom=0.12)
    image = ax.imshow(
        matrix,
        aspect="auto",
        cmap=plt.matplotlib.colors.ListedColormap(["#eceff0", "#3b8b61"]),
        vmin=0,
        vmax=1,
    )
    del image
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels([label for _, label in columns], fontsize=10)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(frame)))
    ax.set_yticklabels(labels, fontsize=9.1)
    for tick, role in zip(ax.get_yticklabels(), frame["evidence_role"]):
        tick.set_color(ROLE_COLORS[role])
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(frame), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    for row in range(len(frame)):
        for column in range(len(columns)):
            ax.text(
                column,
                row,
                "yes" if matrix[row, column] else "no",
                ha="center",
                va="center",
                fontsize=8.2,
                color="white" if matrix[row, column] else "#60696d",
                weight="bold" if matrix[row, column] else "normal",
            )
    for index in range(1, len(frame)):
        if frame.iloc[index]["tissue"] != frame.iloc[index - 1]["tissue"]:
            ax.axhline(index - 0.5, color="#586268", linewidth=1.4)
    ax.set_xlim(-0.5, len(columns) + 2.6)
    ax.text(len(columns) + 0.05, -0.8, "GSEA FDR", fontsize=10, weight="bold")
    ax.text(len(columns) + 1.25, -0.8, "Project direction", fontsize=10, weight="bold")
    for index, row in enumerate(frame.itertuples(index=False)):
        q_text = "<0.001" if float(row.gsea_fdr) == 0 else f"{float(row.gsea_fdr):.3f}"
        ax.text(len(columns) + 0.05, index, q_text, va="center", fontsize=8.8)
        project_text = (
            f"{int(row.expimap_projects_positive)}/{int(row.expimap_n_projects)} higher"
            if row.seed_effect_median >= 0
            else f"{int(row.expimap_projects_negative)}/{int(row.expimap_n_projects)} lower"
        )
        ax.text(len(columns) + 1.25, index, project_text, va="center", fontsize=8.8)
    ax.set_title(
        "Directional evidence for retained main and secondary pathways",
        fontsize=15,
        weight="bold",
        pad=34,
    )
    fig.text(
        0.5,
        0.025,
        "Green cells indicate directional support. GSEA FDR is reported separately because the five-check framework does not use FDR as a binary gate. "
        "Kidney directions pass the five directional checks but remain exploratory because GSEA FDR is above 0.05 and composition adjustment attenuates the effects.",
        ha="center",
        fontsize=8.8,
        wrap=True,
    )
    save_figure(fig, "figure_3_evidence_map")
    plt.close(fig)


def plot_conceptual_summary() -> None:
    panels = (
        (
            "Thymus",
            "Established context",
            "Thymic involution, reduced proliferation, and adaptive immune suppression.",
            "Retained pathway layer",
            "Lower DNA repair and cytoskeletal regulation, with lower lymphoid-stromal interaction as a model-specific niche hypothesis.",
            "Cell-resolved validation is needed to separate thymocyte state from stromal composition.",
        ),
        (
            "Skin",
            "Established context",
            "Barrier injury, inflammation, oxidative stress, and impaired regeneration.",
            "Retained pathway layer",
            "Lower chromatin regulation, DNA repair, Hedgehog, sphingolipid metabolism, and cell-junction organization describe a coordinated maintenance-state decrease.",
            "Protocol depooling shows that gravity, site, recovery, strain, and duration modify the aggregate response.",
        ),
        (
            "Liver",
            "Established context",
            "Lipid, xenobiotic, insulin-related, and mitochondrial dysregulation.",
            "Retained pathway layer",
            "Lower MHC class II antigen presentation and T-cell receptor signaling add a reproducible adaptive-immune communication axis.",
            "Bulk data cannot distinguish immune-cell abundance from signaling state.",
        ),
        (
            "Spleen",
            "Established context",
            "Reduced T-cell abundance and activation after spaceflight.",
            "Retained pathway layer",
            "Lower T-cell receptor, neutrophil degranulation, and C-type lectin receptor programs connect adaptive suppression with lower innate effector and pathogen-sensing states.",
            "The scores are transcriptomic programs, not direct functional measurements of neutrophils or pathogen response.",
        ),
    )
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 9.0), constrained_layout=True)
    for ax, panel in zip(axes.flat, panels):
        tissue, known_heading, known, added_heading, added, caution = panel
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.add_patch(
            FancyBboxPatch(
                (0.015, 0.02),
                0.97,
                0.95,
                boxstyle="round,pad=0.012,rounding_size=0.015",
                facecolor="#fbfcfc",
                edgecolor="#879195",
                linewidth=1.0,
            )
        )
        ax.add_patch(
            FancyBboxPatch(
                (0.04, 0.76),
                0.92,
                0.16,
                boxstyle="round,pad=0.008,rounding_size=0.01",
                facecolor="#edf1f2",
                edgecolor="none",
            )
        )
        ax.text(0.055, 0.88, tissue, fontsize=15, weight="bold", va="top")
        ax.text(0.055, 0.70, known_heading, fontsize=9.2, weight="bold", color="#4c575c")
        ax.text(
            0.055,
            0.63,
            textwrap.fill(known, width=70),
            fontsize=9.0,
            va="top",
            linespacing=1.35,
        )
        ax.text(
            0.055,
            0.43,
            added_heading,
            fontsize=9.3,
            weight="bold",
            color=TISSUE_COLORS[tissue.lower()],
        )
        ax.text(
            0.055,
            0.36,
            textwrap.fill(added, width=70),
            fontsize=9.1,
            va="top",
            linespacing=1.35,
        )
        ax.text(
            0.055,
            0.08,
            textwrap.fill("Interpretation boundary: " + caution, width=78),
            fontsize=8.0,
            color="#6a4a22",
            va="bottom",
            linespacing=1.25,
        )
    fig.suptitle(
        "Tissue-specific interpretation after project, seed, method, and composition review",
        fontsize=16,
        weight="bold",
    )
    fig.text(
        0.5,
        -0.015,
        "Conceptual synthesis only. Lower refers to the decoder-oriented flight-minus-ground pathway score and does not establish a causal mechanism.",
        ha="center",
        fontsize=8.8,
    )
    save_figure(fig, "figure_5_complementary_process_model")
    plt.close(fig)


def run() -> None:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    scope = build_model_scope()
    evidence = load_retained_evidence()
    scope.to_csv(
        SOURCE_DIR / "table_s25_revised_model_scope.tsv", sep="\t", index=False
    )
    write_retained_table(evidence)
    copy_reassessment_tables()
    copy_reassigned_figures()
    plot_workflow(scope)
    plot_main_pathway_shifts(evidence)
    plot_evidence_matrix(evidence)
    plot_conceptual_summary()
    print(
        json.dumps(
            {
                "main_tissues": list(MAIN_TISSUES),
                "secondary_tissues": ["kidney"],
                "non_advanced_screening_tissues": ["soleus"],
                "retained_pathways": int(len(evidence)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    run()
