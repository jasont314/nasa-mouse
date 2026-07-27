"""Create interpretation-colored copies of HVG expiMap heatmaps.

The output intentionally keeps the original all-label heatmap layout and only
changes y-axis pathway label colors by a manually curated interpretation
category.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]


CATEGORY_COLORS = {
    "prior_literature_match": "#16803c",
    "plausible_complementary": "#1f63b5",
    "broad_or_uncertain_signal": "#a6611a",
    "likely_artifact": "#c51b29",
    "low_or_no_effect": "#8f8f8f",
}

CATEGORY_LABELS = {
    "prior_literature_match": "prior literature match",
    "plausible_complementary": "plausible/complementary",
    "broad_or_uncertain_signal": "broad/uncertain signal",
    "likely_artifact": "likely artifact/junk",
    "low_or_no_effect": "low/no relative effect",
}

BASE_ARTIFACT_PATTERNS = (
    "olfactory",
    "sensory perception",
    "taste",
    "visual phototransduction",
    "phototransduction",
    "gaba",
    "nmda",
    "opioid",
    "fertilization",
    "reproduction",
)


EXACT_OVERRIDES = {
    "liver_hvg": {
        "R-MMU-211897_CYTOCHROME_P450_ARRANGED_BY_SUBSTRATE_TYPE": (
            "prior_literature_match",
            "literature override: liver CYP450/xenobiotic metabolism",
        ),
        "R-MMU-168256_IMMUNE_SYSTEM": (
            "prior_literature_match",
            "literature override: liver immune/inflammatory remodeling",
        ),
        "R-MMU-191273_CHOLESTEROL_BIOSYNTHESIS": (
            "prior_literature_match",
            "literature override: liver lipid/cholesterol metabolism",
        ),
        "R-MMU-156590_GLUTATHIONE_CONJUGATION": (
            "prior_literature_match",
            "literature override: liver detoxification/oxidative defense",
        ),
        "R-MMU-156588_GLUCURONIDATION": (
            "prior_literature_match",
            "literature override: liver phase II conjugation",
        ),
        "R-MMU-211945_PHASE_I_FUNCTIONALIZATION_OF_COMPOUNDS": (
            "prior_literature_match",
            "literature override: liver phase I xenobiotic metabolism",
        ),
        "R-MMU-156580_PHASE_II_CONJUGATION_OF_COMPOUNDS": (
            "prior_literature_match",
            "literature override: liver phase II xenobiotic metabolism",
        ),
    },
    "skin_hvg": {
        "R-MMU-6805567_KERATINIZATION": (
            "prior_literature_match",
            "literature override: skin barrier/keratinization",
        ),
        "R-MMU-6809371_FORMATION_OF_THE_CORNIFIED_ENVELOPE": (
            "prior_literature_match",
            "literature override: skin cornified envelope/barrier",
        ),
        "R-MMU-2187338_VISUAL_PHOTOTRANSDUCTION": (
            "likely_artifact",
            "literature override: off-context visual pathway in skin",
        ),
        "R-MMU-397014_MUSCLE_CONTRACTION": (
            "likely_artifact",
            "literature override: off-context muscle pathway in skin",
        ),
        "R-MMU-390522_STRIATED_MUSCLE_CONTRACTION": (
            "likely_artifact",
            "literature override: off-context striated muscle pathway in skin",
        ),
        "R-MMU-5576891_CARDIAC_CONDUCTION": (
            "likely_artifact",
            "literature override: off-context cardiac conduction pathway in skin",
        ),
    },
    "thymus_hvg": {
        "R-MMU-1640170_CELL_CYCLE": (
            "prior_literature_match",
            "literature override: thymic cell-cycle remodeling",
        ),
        "R-MMU-69278_CELL_CYCLE_MITOTIC": (
            "prior_literature_match",
            "literature override: thymic mitotic/cell-cycle remodeling",
        ),
        "R-MMU-73894_DNA_REPAIR": (
            "prior_literature_match",
            "literature override: thymic DNA repair/stress response",
        ),
        "R-MMU-202403_TCR_SIGNALING": (
            "prior_literature_match",
            "literature override: thymic T-cell receptor biology",
        ),
        "R-MMU-202427_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS": (
            "prior_literature_match",
            "literature override: thymic TCR/CD3 signaling",
        ),
        "R-MMU-168256_IMMUNE_SYSTEM": (
            "prior_literature_match",
            "literature override: thymic immune system remodeling",
        ),
        "R-MMU-9769740_COAGULATION_PATHWAY": (
            "plausible_complementary",
            "literature override: platelet/coagulation adjunct, not core thymus axis",
        ),
    },
    "soleus_hvg": {
        "R-MMU-6809371_FORMATION_OF_THE_CORNIFIED_ENVELOPE": (
            "likely_artifact",
            "literature override: off-context skin barrier pathway in soleus",
        ),
        "R-MMU-6805567_KERATINIZATION": (
            "likely_artifact",
            "literature override: off-context skin keratinization pathway in soleus",
        ),
        "R-MMU-168256_IMMUNE_SYSTEM": (
            "prior_literature_match",
            "literature override: muscle immune/inflammatory response",
        ),
        "R-MMU-166663_INITIAL_TRIGGERING_OF_COMPLEMENT": (
            "prior_literature_match",
            "literature override: muscle complement/inflammatory response",
        ),
        "R-MMU-8978868_FATTY_ACID_METABOLISM": (
            "prior_literature_match",
            "literature override: muscle fatty-acid metabolism",
        ),
        "R-MMU-1442490_COLLAGEN_DEGRADATION": (
            "prior_literature_match",
            "literature override: muscle ECM/collagen remodeling",
        ),
        "R-MMU-397014_MUSCLE_CONTRACTION": (
            "prior_literature_match",
            "literature override: muscle contraction",
        ),
        "R-MMU-390522_STRIATED_MUSCLE_CONTRACTION": (
            "prior_literature_match",
            "literature override: striated muscle contraction",
        ),
        "R-MMU-211897_CYTOCHROME_P450_ARRANGED_BY_SUBSTRATE_TYPE": (
            "low_or_no_effect",
            "literature override: not a prominent soleus heatmap effect",
        ),
        "R-MMU-202427_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS": (
            "low_or_no_effect",
            "literature override: low-effect immune/TCR row in soleus heatmap",
        ),
    },
}


@dataclass(frozen=True)
class ModelConfig:
    name: str
    title: str
    run_dir: Path
    prior_patterns: tuple[str, ...]
    plausible_patterns: tuple[str, ...]
    artifact_patterns: tuple[str, ...] = ()


MODELS = (
    ModelConfig(
        name="liver_hvg",
        title="Liver HVG expiMap FLT-GC pathway shifts",
        run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/query_nb_250epoch_seed2020",
        prior_patterns=(
            "cytochrome p450",
            "phase i",
            "phase ii",
            "drug adme",
            "xenobiotic",
            "lipid",
            "lipoprotein",
            "phospholipid",
            "bile acid",
            "bile salt",
            "steroid",
            "glucose",
            "insulin",
            "insulin secretion",
            "glutathione",
            "glucuronidation",
            "cholesterol",
            "nr1h",
            "ppar",
            "fatty",
            "peroxisom",
            "mitochond",
            "respiratory",
            "oxidative",
            "autophagy",
            "proteasome",
            "amino acids",
            "vitamins and cofactors",
            "detoxification",
            "immune system",
            "toll like receptor",
            "tlr",
            "interleukin",
            "interferon",
            "complement",
            "mhc",
            "antigen",
            "neutrophil",
        ),
        plausible_patterns=(
            "tcr",
            "pd 1",
            "co inhibition",
            "cytokine",
            "immune",
            "mapk",
            "akt",
            "rho",
            "rac",
            "gtpase",
            "gpcr",
            "vascular wall",
            "ecm",
            "extracellular matrix",
            "collagen",
            "integrin",
            "cadherin",
            "proteoglycans",
            "dna repair",
            "telomere",
            "senescence",
            "apoptosis",
            "atm",
            "base excision",
            "platelet",
            "coagulation",
            "hemostasis",
            "cell cycle",
            "ptk6",
        ),
        artifact_patterns=(
            "keratinization",
            "cornified",
            "keratan",
            "cardiac conduction",
            "muscle contraction",
            "neuronal",
            "synapse",
            "postsynaptic",
            "primary cilium",
            "cilium",
        ),
    ),
    ModelConfig(
        name="skin_hvg",
        title="Skin HVG expiMap FLT-GC pathway shifts",
        run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_skin/tutorial_hvg_2000/query_nb_250epoch_seed2020",
        prior_patterns=(
            "keratinization",
            "cornified envelope",
            "barrier",
            "epiderm",
            "collagen",
            "extracellular matrix",
            "ecm",
            "integrin",
            "cell cell junction",
            "cadherin",
            "gap junction",
            "hyaluronan",
            "dna repair",
            "dna damage",
            "dna double strand",
            "telomere",
            "senescence",
            "cell cycle",
            "mitotic",
            "apoptosis",
            "programmed cell death",
            "interleukin 1",
            "interleukin 17",
            "interleukin",
            "interferon",
            "cytokine",
            "mhc",
            "antigen",
            "innate immune",
            "adaptive immune",
            "immune",
            "glutathione",
            "phase ii",
            "melanocyte",
            "mitochond",
            "respiratory",
            "reactive oxygen",
            "oxidative",
            "cholesterol",
            "lipid",
            "sphingolipid",
            "steroid",
            "glucose",
            "insulin",
            "estrogen",
        ),
        plausible_patterns=(
            "retinoic acid",
            "rho",
            "rac",
            "gtpase",
            "actin",
            "gpi anchored",
            "sumoylation",
            "protein modification",
            "membrane trafficking",
            "gpcr",
            "g alpha",
            "hedgehog",
            "developmental biology",
            "cilium",
            "vitamins and cofactors",
            "chromatin",
            "neddylation",
            "glycosylation",
            "mucin",
            "vesicle",
            "clathrin",
            "endocytosis",
            "potassium channels",
        ),
        artifact_patterns=(
            "rhodopsin",
            "muscle contraction",
            "striated muscle contraction",
            "smooth muscle contraction",
            "cardiac conduction",
            "neuronal",
            "synapse",
            "postsynaptic",
        ),
    ),
    ModelConfig(
        name="thymus_hvg",
        title="Thymus HVG expiMap FLT-GC pathway shifts",
        run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_thymus/tutorial_hvg_2000/query_nb_250epoch_seed2020",
        prior_patterns=(
            "immune",
            "cytokine",
            "interleukin",
            "tcr",
            "cd3",
            "zap",
            "bcr",
            "antigen",
            "mhc",
            "interferon",
            "tlr",
            "toll",
            "myd88",
            "adaptive",
            "innate",
            "cell cycle",
            "mitotic",
            "proliferation",
            "checkpoint",
            "g1",
            "g2",
            "m phase",
            "dna replication",
            "dna repair",
            "dna damage",
            "telomere",
            "apoptosis",
            "programmed cell death",
            "death receptor",
            "nucleotide",
            "chromatin",
            "transcription",
            "cellular response to stress",
        ),
        plausible_patterns=(
            "coagulation",
            "platelet",
            "hemostasis",
            "neutrophil",
            "complement",
            "neddylation",
            "rho",
            "rac",
            "mapk",
            "wnt",
            "tgf beta",
            "gpcr",
            "g alpha",
            "gtpase",
            "actin",
            "phagocytic",
            "golgi",
            "er",
            "copi",
            "er retrograde",
            "vesicle",
            "membrane trafficking",
            "glycosaminoglycan",
            "steroid",
            "adipogenesis",
            "esr mediated",
            "metabolism of proteins",
            "protein modification",
            "ubiquitin",
            "sumoylation",
            "transport of small molecules",
        ),
        artifact_patterns=(
            "keratinization",
            "cornified envelope",
            "muscle contraction",
            "cardiac conduction",
            "neuronal",
            "synapse",
            "postsynaptic",
        ),
    ),
    ModelConfig(
        name="soleus_hvg",
        title="Soleus HVG expiMap FLT-GC pathway shifts",
        run_dir=ROOT
        / "outputs/expimap_muscle_targeted_combined_min8/tutorial_hvg_soleus_2000/query_nb_250epoch_seed2020",
        prior_patterns=(
            "muscle contraction",
            "striated muscle contraction",
            "myofibril",
            "calcium",
            "potassium channel",
            "ion channel",
            "igf",
            "insulin like growth factor",
            "akt",
            "pip3",
            "pi 3k",
            "pi3k",
            "fatty acid",
            "lipid",
            "lipoprotein",
            "carbohydrate",
            "glucose",
            "insulin",
            "estrogen",
            "mitochond",
            "respiratory",
            "oxidative",
            "ros",
            "reactive oxygen",
            "detoxification",
            "apoptosis",
            "programmed cell death",
            "ecm",
            "extracellular matrix",
            "collagen",
            "integrin",
            "pdgf",
            "tgf beta",
            "muscle",
            "myogenesis",
            "dystrophin",
            "dag1",
            "matriglycan",
            "nitric oxide",
            "autophagy",
            "ubiquitin",
            "proteasome",
            "amino acids",
            "immune system",
            "complement",
            "tlr",
            "interferon",
            "interleukin",
            "cytokine",
            "neutrophil",
        ),
        plausible_patterns=(
            "wnt",
            "arachidonate",
            "gap junction",
            "vesicle",
            "chromatin",
            "dna repair",
            "double strand break",
            "hemostasis",
            "coagulation",
            "scavenger receptors",
            "tcr",
            "cd3",
            "g alpha",
            "gpcr",
            "rho",
            "rac",
            "mapk",
            "gtpase",
            "neuronal",
            "synapse",
            "postsynaptic",
            "cardiac conduction",
            "membrane",
            "transport",
            "cell cycle",
            "m phase",
            "g2",
            "g1",
            "calcium",
        ),
        artifact_patterns=(
            "keratinization",
            "cornified envelope",
            "keratan",
            "rhodopsin",
            "sperm",
            "meiosis",
        ),
    ),
)


def normalize_term(term: str) -> str:
    return term.replace("_", " ").lower()


def first_match(text: str, patterns: tuple[str, ...]) -> str:
    for pattern in patterns:
        if pattern in text:
            return pattern
    return ""


def categorize_terms(
    config: ModelConfig,
    comparison: pd.DataFrame,
    matrix: pd.DataFrame,
) -> tuple[pd.DataFrame, float]:
    comparison = comparison.copy()
    comparison["abs_effect"] = comparison["flight_minus_ground"].abs()
    matrix_values = matrix.apply(pd.to_numeric, errors="coerce")
    matrix_stats = pd.DataFrame(
        {
            "term": matrix.index.astype(str),
            "mean_abs_accession_effect": matrix_values.abs().mean(axis=1).to_numpy(),
            "abs_mean_accession_effect": matrix_values.mean(axis=1).abs().to_numpy(),
        }
    )
    comparison = comparison.merge(matrix_stats, on="term", how="left")
    low_threshold = 0.10 * float(comparison["mean_abs_accession_effect"].quantile(0.95))
    if not np.isfinite(low_threshold) or low_threshold <= 0:
        low_threshold = 0.0

    records = []
    artifact_patterns = BASE_ARTIFACT_PATTERNS + config.artifact_patterns
    exact_overrides = EXACT_OVERRIDES.get(config.name, {})
    for row in comparison.itertuples(index=False):
        text = normalize_term(row.term)
        row_strength = row.mean_abs_accession_effect
        if not np.isfinite(row_strength):
            row_strength = abs(row.flight_minus_ground)
        is_low = bool(row_strength < low_threshold)

        category = "low_or_no_effect"
        rationale = "below model-relative accession-level heatmap effect threshold"
        matched = ""

        override = exact_overrides.get(row.term)
        if override and override[0] == "low_or_no_effect":
            category, rationale = override
            matched = "exact override"
        elif override and (not is_low or override[0] == "likely_artifact"):
            category, rationale = override
            matched = "exact override"
        elif not is_low:
            matched = first_match(text, artifact_patterns)
            if matched:
                category = "likely_artifact"
                rationale = f"context-incongruent Reactome term: {matched}"
            else:
                matched = first_match(text, config.prior_patterns)
                if matched:
                    category = "prior_literature_match"
                    rationale = f"matches expected tissue/spaceflight biology: {matched}"
                else:
                    matched = first_match(text, config.plausible_patterns)
                    if matched:
                        category = "plausible_complementary"
                        rationale = f"plausible complementary mechanism: {matched}"
                    else:
                        category = "broad_or_uncertain_signal"
                        rationale = "moderate/significant effect but not manually assigned to a literature family"

        records.append(
            {
                "term": row.term,
                "flight_minus_ground": row.flight_minus_ground,
                "welch_fdr": row.welch_fdr,
                "abs_effect": abs(row.flight_minus_ground),
                "mean_abs_accession_effect": row.mean_abs_accession_effect,
                "abs_mean_accession_effect": row.abs_mean_accession_effect,
                "low_effect_threshold": low_threshold,
                "category": category,
                "category_label": CATEGORY_LABELS[category],
                "matched_keyword": matched,
                "rationale": rationale,
            }
        )
    return pd.DataFrame(records), low_threshold


def plot_heatmap(config: ModelConfig) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    analysis_dir = config.run_dir / "analysis"
    matrix = pd.read_csv(
        analysis_dir / "all_program_accession_flt_minus_gc_matrix_signed_order.tsv",
        sep="\t",
    ).set_index("term")
    comparison = pd.read_csv(analysis_dir / "flt_vs_gc_pathway_comparison.tsv", sep="\t")
    labels, low_threshold = categorize_terms(config, comparison, matrix)
    label_lookup = labels.set_index("term")
    labels.to_csv(analysis_dir / "pathway_interpretation_labels.tsv", sep="\t", index=False)

    ordered_labels = matrix.index.astype(str).tolist()
    row_categories = [
        label_lookup.loc[term, "category"] if term in label_lookup.index else "low_or_no_effect"
        for term in ordered_labels
    ]
    row_colors = [CATEGORY_COLORS[category] for category in row_categories]

    data = matrix.to_numpy(dtype=float)
    vmax = 2.5

    n_rows, n_cols = matrix.shape
    fig_width = max(18.0, min(32.0, 15.0 + 1.0 * n_cols))
    fig_height = min(max(8.0, n_rows * 0.08 + 2.0), 140.0)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#f2f2f2")
    image = ax.imshow(
        np.ma.masked_invalid(data),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=-vmax,
        vmax=vmax,
    )

    accessions = matrix.columns.astype(str).tolist()
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(ordered_labels, fontsize=4.0)
    for tick, color in zip(ax.get_yticklabels(), row_colors):
        tick.set_color(color)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.set_xticks(np.arange(len(accessions)))
    ax.set_xticklabels(
        accessions,
        rotation=90,
        ha="center",
        va="top",
        fontsize=7.0,
    )
    ax.set_xlabel("OSD accession / study", fontsize=8)
    ax.set_title("All expiMap FLT-GC pathway shifts, signed-effect row order", fontsize=10, pad=10)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.015)
    cbar.set_label("mean FLT - mean GC pathway score", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.text(
        0.63,
        0.012,
        "Rows ordered by signed study-mean FLT-GC effect (FLT-up to FLT-down)",
        ha="center",
        va="bottom",
        fontsize=6,
    )
    fig.subplots_adjust(left=0.56, right=0.89, top=0.985, bottom=0.055)

    output = analysis_dir / "all_program_accession_flt_minus_gc_heatmap_signed_order_all_labels_interpretation.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)

    presentation_dir = ROOT / "presentation/expimap/annotated_hvg"
    presentation_dir.mkdir(parents=True, exist_ok=True)
    presentation_output = presentation_dir / f"{config.name}_interpretation_heatmap.png"
    shutil.copy2(output, presentation_output)
    shutil.copy2(
        analysis_dir / "pathway_interpretation_labels.tsv",
        presentation_dir / f"{config.name}_pathway_interpretation_labels.tsv",
    )
    return output


def main() -> None:
    for config in MODELS:
        output = plot_heatmap(config)
        print(output.relative_to(ROOT))


if __name__ == "__main__":
    main()
