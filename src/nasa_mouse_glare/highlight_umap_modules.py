"""Create UMAP copies with prior-work and hidden-module legend highlights."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

from .multi_tissue_reports import clean_reactome_term, cluster_colors


DEFAULT_ROOT = Path("outputs/glare/multi_tissue_api")
DEFAULT_VALIDATION_DIR = DEFAULT_ROOT / "validation_stack_terms15"

PRIOR_COLOR = "#2563eb"
HIDDEN_COLOR = "#f97316"
OTHER_COLOR = "#6b7280"
UNCLEAR_COLOR = "#9ca3af"

HIDDEN_PRIORITY_TERMS = {
    "skeletal_muscle": {
        "Bmal1 Clock Npas2 Activates Circadian Expression",
        "Cyclin E Associated Events During G1 S Transition",
    },
    "skeletal_muscle_soleus": {
        "Signalling By Ngf",
    },
    "kidney": {
        "Membrane Trafficking",
        "Signaling By Insulin Receptor",
    },
    "thymus": {
        "Response To Elevated Platelet Cytosolic Ca2",
    },
}

HIDDEN_TERM_THEMES = {
    "Bmal1 Clock Npas2 Activates Circadian Expression": "circadian clock/regulatory",
    "Cyclin E Associated Events During G1 S Transition": "G1/S cell-cycle/DNA replication",
    "Signalling By Ngf": "NGF/neurotrophic signaling",
    "Membrane Trafficking": "membrane/vesicle trafficking",
    "Signaling By Insulin Receptor": "insulin receptor/metabolic signaling",
    "Response To Elevated Platelet Cytosolic Ca2": "platelet-calcium/vascular remodeling",
}

THEME_RULES = [
    (
        "circadian clock/regulatory",
        ("circadian", "bmal1", "clock", "npas2", "rora", "nr1d1"),
    ),
    (
        "G1/S cell-cycle/DNA replication",
        ("cell cycle", "g1 s", "cyclin", "mitotic", "dna", "cdc6", "cdt1", "orc", "p53"),
    ),
    (
        "translation/RNA processing",
        ("translation", "peptide chain", "mrna", "rna", "splicing", "utr", "eif", "srp", "nonsense mediated"),
    ),
    (
        "mitochondrial/OXPHOS metabolism",
        ("tca cycle", "respiratory electron", "atp synthesis", "mitochondrial", "oxidative phosphorylation"),
    ),
    (
        "immune/antigen/cytokine signaling",
        ("immune", "cytokine", "interferon", "antigen", "mhc", "tcr", "b cell", "ils"),
    ),
    (
        "platelet/hemostasis/vascular",
        ("platelet", "hemostasis", "clot", "gpvi", "vascular", "endothelial"),
    ),
    (
        "membrane/vesicle trafficking",
        ("membrane trafficking", "vesicle", "golgi", "endosome", "lysosome", "transport"),
    ),
    (
        "lipid/xenobiotic/detox metabolism",
        ("lipid", "fatty acid", "ketone", "xenobiotic", "phase", "glucuronidation", "bile"),
    ),
    (
        "muscle/contractile/cytoskeleton",
        ("muscle contraction", "striated", "actin", "tubulin", "cytoskeleton", "sarcomere"),
    ),
    (
        "growth-factor/developmental signaling",
        ("ngf", "insulin receptor", "fgfr", "erbb", "axon", "developmental", "wnt", "tgf"),
    ),
    (
        "protein modification/glycosylation",
        ("post translational", "glycosylation", "n linked", "protein modification"),
    ),
    (
        "neuronal/sensory",
        ("olfactory", "neuronal", "rhodopsin", "sensory"),
    ),
]

BROAD_TOP_TERMS = {
    "Adaptive Immune System",
    "Immune System",
    "Metabolism",
    "Metabolism Of Proteins",
    "Transcription",
    "Translation",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Copy aggregate GLARE UMAPs with legend entries highlighted as "
            "prior-work aligned or hidden/novel."
        )
    )
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--validation-dir", type=Path, default=DEFAULT_VALIDATION_DIR)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument(
        "--scopes",
        nargs="+",
        default=["aggregate", "aggregate_mober"],
        help="Run scopes to scan under each tissue directory.",
    )
    return parser.parse_args()


def load_prior_terms(validation_dir: Path) -> dict[str, set[str]]:
    prior_terms: dict[str, set[str]] = {}
    candidate_modules = validation_dir / "candidate_modules.tsv"
    if candidate_modules.exists():
        table = pd.read_csv(candidate_modules, sep="\t")
        table = table.loc[table["module_class"].astype(str).eq("intersection")]
        for tissue, group in table.groupby("tissue"):
            prior_terms.setdefault(str(tissue), set()).update(group["clean_term"].dropna())

    candidate_summary = validation_dir / "novel_candidate_validation" / "candidate_validation_summary.tsv"
    if candidate_summary.exists():
        table = pd.read_csv(candidate_summary, sep="\t")
        known = table.loc[
            table["candidate_group"].astype(str).str.contains("validation_anchor", na=False)
            | table["theme"].astype(str).str.contains("known", case=False, na=False)
        ]
        for tissue, group in known.groupby("tissue"):
            prior_terms.setdefault(str(tissue), set()).update(group["clean_term"].dropna())
    return prior_terms


def load_direction_meta(validation_dir: Path) -> dict[tuple[str, str], dict]:
    path = validation_dir / "candidate_module_score_meta.tsv"
    if not path.exists():
        return {}
    table = pd.read_csv(path, sep="\t")
    directions: dict[tuple[str, str], dict] = {}
    for row in table.itertuples(index=False):
        data = row._asdict()
        tissue = str(data.get("tissue", ""))
        clean_term = str(data.get("clean_term", ""))
        if tissue and clean_term:
            directions[(tissue, clean_term)] = data
    return directions


def tissue_from_run(run_dir: Path) -> str:
    if run_dir.name in {"aggregate", "aggregate_mober"}:
        return run_dir.parent.name
    if run_dir.parent.name == "per_study":
        return run_dir.parent.parent.name
    return run_dir.parent.name


def load_enrichment_terms(run_dir: Path, location: str, alpha: float) -> pd.DataFrame:
    path = run_dir / "clustering" / "reactome_enrichment.tsv"
    if not path.exists():
        return pd.DataFrame(columns=["cluster", "clean_term", "fdr_bh"])
    table = pd.read_csv(path, sep="\t")
    if table.empty:
        return pd.DataFrame(columns=["cluster", "clean_term", "fdr_bh"])
    table = table.loc[table["location"].astype(str).eq(location)].copy()
    table = table.loc[pd.to_numeric(table["fdr_bh"], errors="coerce").le(alpha)].copy()
    if table.empty:
        return pd.DataFrame(columns=["cluster", "clean_term", "fdr_bh"])
    table["cluster"] = table["cluster"].astype(int)
    table["clean_term"] = table["term"].map(clean_reactome_term)
    keep = [
        column
        for column in [
            "cluster",
            "clean_term",
            "fdr_bh",
            "p_value",
            "overlap",
            "cluster_genes",
            "term_genes_in_universe",
        ]
        if column in table.columns
    ]
    return table[keep]


def rollup_themes(terms: list[str]) -> list[str]:
    themes = []
    for term in terms:
        lower_term = str(term).lower()
        for theme, needles in THEME_RULES:
            if theme in themes:
                continue
            if any(needle in lower_term for needle in needles):
                themes.append(theme)
    return themes


def format_term_list(terms: list[str], max_terms: int) -> str:
    return "; ".join(str(term) for term in terms[:max_terms])


def direction_summary(
    tissue: str,
    terms: list[str],
    direction_meta: dict[tuple[str, str], dict],
) -> str:
    parts = []
    for term in terms:
        row = direction_meta.get((tissue, term))
        if not row:
            continue
        try:
            delta = float(row.get("mean_flight_minus_ground"))
            fdr = float(row.get("combined_welch_fdr_bh"))
        except (TypeError, ValueError):
            continue
        direction = "FLT higher" if delta > 0 else "FLT lower" if delta < 0 else "no shift"
        studies = row.get("studies_tested", "")
        parts.append(f"{term}: {direction} ({delta:+.3g}, FDR={fdr:.2g}, studies={studies})")
    return "; ".join(parts)


def build_caveats(
    tissue: str,
    row: pd.Series,
    all_themes: list[str],
    hidden_matches: list[str],
    hidden_fdrs: list[float],
) -> str:
    caveats = []
    gene_count = int(row.get("gene_count", 0))
    top_term = str(row.get("best_reactome_term", ""))
    if gene_count >= 3000:
        caveats.append("large mixed cluster")
    if len(all_themes) >= 4:
        caveats.append("multiple enriched themes; top label alone is insufficient")
    if top_term in BROAD_TOP_TERMS:
        caveats.append("broad top Reactome label")
    if hidden_matches and top_term and top_term not in hidden_matches:
        caveats.append("hidden term differs from top Reactome label")
    if any(fdr > 0.01 for fdr in hidden_fdrs):
        caveats.append("one or more hidden-term enrichments are secondary/weak")
    if (
        tissue not in {"lung", "spleen", "thymus"}
        and "immune/antigen/cytokine signaling" in all_themes
    ):
        caveats.append("immune/composition-sensitive signal")
    if not bool(row.get("reactome_significant", False)):
        caveats.append("no significant Reactome label at threshold")
    return "; ".join(dict.fromkeys(caveats))


def build_legend_description(row: pd.Series) -> str:
    primary = str(row.get("primary_theme", ""))
    secondary = str(row.get("secondary_themes", ""))
    top = str(row.get("top_reactome_term", ""))
    caveat = str(row.get("caveat", ""))
    parts = [primary]
    if secondary and secondary != "nan":
        secondary_terms = [term.strip() for term in secondary.split(";") if term.strip()]
        parts.append(f"other: {'; '.join(secondary_terms[:3])}")
    if top and top != "nan" and top not in primary:
        parts.append(f"top: {top}")
    short_flags = []
    if "large mixed cluster" in caveat:
        short_flags.append("large/mixed")
    if "secondary/weak" in caveat:
        short_flags.append("weak secondary hidden term")
    if "top label alone is insufficient" in caveat:
        short_flags.append("multi-theme")
    if short_flags:
        parts.append(", ".join(short_flags))
    return "; ".join(part for part in parts if part)


def classify_clusters(
    annotations: pd.DataFrame,
    enrichment: pd.DataFrame,
    tissue: str,
    prior_terms: dict[str, set[str]],
    direction_meta: dict[tuple[str, str], dict],
) -> pd.DataFrame:
    hidden = HIDDEN_PRIORITY_TERMS.get(tissue, set())
    prior = prior_terms.get(tissue, set())
    rows = []
    for _, row in annotations.iterrows():
        cluster = int(row["cluster"])
        cluster_terms = enrichment.loc[enrichment["cluster"].eq(cluster)].copy()
        if not cluster_terms.empty:
            sort_columns = [
                column for column in ("fdr_bh", "p_value") if column in cluster_terms.columns
            ]
            if sort_columns:
                cluster_terms = cluster_terms.sort_values(sort_columns, na_position="last")
        terms = cluster_terms["clean_term"].dropna().astype(str).tolist()
        hidden_matches = list(dict.fromkeys(term for term in terms if term in hidden))
        prior_matches = list(dict.fromkeys(term for term in terms if term in prior))
        all_themes = rollup_themes(terms)
        hidden_theme_labels = [
            HIDDEN_TERM_THEMES.get(term, term) for term in hidden_matches
        ]
        primary_theme_terms = hidden_theme_labels if hidden_theme_labels else all_themes[:1]
        if not primary_theme_terms and terms:
            primary_theme_terms = terms[:1]
        if not primary_theme_terms:
            primary_theme_terms = [str(row.get("cluster_description", "ambiguous"))]
        primary_theme = "; ".join(primary_theme_terms)
        primary_theme_set = set(primary_theme_terms)
        secondary_themes = [
            theme for theme in all_themes if theme not in primary_theme_set
        ][:5]
        matched_terms = hidden_matches if hidden_matches else prior_matches
        hidden_fdrs = (
            cluster_terms.loc[cluster_terms["clean_term"].isin(hidden_matches), "fdr_bh"]
            .dropna()
            .astype(float)
            .tolist()
        )
        if hidden_matches:
            category = "hidden_novel"
            category_label = "Hidden/novel"
            category_color = HIDDEN_COLOR
        elif prior_matches:
            category = "prior_work"
            category_label = "Prior-work aligned"
            category_color = PRIOR_COLOR
        elif str(row.get("annotation_status", "")) == "ambiguous":
            category = "unclear"
            category_label = "Unclear/ambiguous"
            category_color = UNCLEAR_COLOR
            matched_terms = []
        else:
            category = "other"
            category_label = "Other"
            category_color = OTHER_COLOR
            matched_terms = []
        data = row.to_dict()
        data.update(
            {
                "primary_theme": primary_theme,
                "secondary_themes": "; ".join(secondary_themes),
                "top_reactome_term": row.get("best_reactome_term", ""),
                "top_reactome_fdr_bh": row.get("best_reactome_fdr_bh", ""),
                "top_reactome_overlap": row.get("best_reactome_overlap", ""),
                "top_enriched_terms": format_term_list(terms, 10),
                "highlight_category": category,
                "highlight_label": category_label,
                "highlight_color": category_color,
                "matched_hidden_terms": "; ".join(hidden_matches),
                "matched_prior_terms": "; ".join(prior_matches),
                "matched_highlight_terms": "; ".join(matched_terms[:4]),
                "matched_direction_summary": direction_summary(
                    tissue, matched_terms, direction_meta
                ),
                "caveat": build_caveats(
                    tissue,
                    row,
                    all_themes,
                    hidden_matches,
                    hidden_fdrs,
                ),
            }
        )
        data["legend_description"] = build_legend_description(pd.Series(data))
        rows.append(data)
    return pd.DataFrame(rows)


def plot_highlighted_umap(
    coords: pd.DataFrame,
    labels: np.ndarray,
    annotations: pd.DataFrame,
    title: str,
    output_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    clusters = annotations["cluster"].astype(int).tolist()
    colors = cluster_colors(clusters)
    grouped = [
        ("hidden_novel", "Hidden/novel modules"),
        ("prior_work", "Prior-work aligned modules"),
        ("other", "Other modules"),
        ("unclear", "Unclear/ambiguous modules"),
    ]

    wrapped_groups = []
    total_sections = 0
    total_label_lines = 0
    for category, legend_title in grouped:
        subset = annotations.loc[annotations["highlight_category"].eq(category)].copy()
        if subset.empty:
            continue
        entries = []
        total_sections += 1
        for row in subset.itertuples(index=False):
            cluster = int(row.cluster)
            description = str(getattr(row, "legend_description", ""))
            if not description or description == "nan":
                description = str(row.cluster_description)
            label = (
                f"C{cluster} ({int(row.gene_count)} genes): "
                f"{description}"
            )
            wrapped = textwrap.wrap(label, width=62)
            entries.append((row, wrapped))
            total_label_lines += max(1, len(wrapped))
        wrapped_groups.append((category, legend_title, entries))

    legend_height = (
        1.0
        + 0.4 * total_sections
        + 0.19 * total_label_lines
        + 0.09 * sum(len(entries) for _, _, entries in wrapped_groups)
    )
    figure_height = max(8.6, legend_height)

    figure = plt.figure(figsize=(15.5, figure_height), dpi=170)
    grid = figure.add_gridspec(1, 2, width_ratios=[2.0, 1.35], wspace=0.04)
    axis = figure.add_subplot(grid[0, 0])
    legend_axis = figure.add_subplot(grid[0, 1])
    legend_axis.set_xlim(0, 1)
    legend_axis.set_ylim(0, 1)
    legend_axis.axis("off")

    for cluster in clusters:
        mask = labels == cluster
        axis.scatter(
            coords.loc[mask, "x"],
            coords.loc[mask, "y"],
            s=3,
            alpha=0.58,
            linewidths=0,
            color=colors[cluster],
        )

    axis.set_title(title, fontsize=11, pad=14)
    axis.set_xlabel("UMAP 1")
    axis.set_ylabel("UMAP 2")

    legend_axis.text(
        0.0,
        0.995,
        "Cluster Legend",
        fontsize=11,
        fontweight="bold",
        va="top",
        transform=legend_axis.transAxes,
    )
    legend_axis.text(
        0.0,
        0.958,
        "Orange text = hidden/novel\nBlue text = prior-work aligned\nDots keep cluster colors",
        fontsize=7.5,
        color="#374151",
        va="top",
        transform=legend_axis.transAxes,
    )

    line_step = 0.14 / figure_height
    entry_gap = 0.07 / figure_height
    section_step = 0.28 / figure_height
    section_gap = 0.12 / figure_height
    y = 0.89
    for category, legend_title, entries in wrapped_groups:
        title_color = {
            "hidden_novel": HIDDEN_COLOR,
            "prior_work": PRIOR_COLOR,
            "other": "#374151",
            "unclear": UNCLEAR_COLOR,
        }[category]
        legend_axis.text(
            0.0,
            max(y, 0.035),
            legend_title,
            fontsize=8.2,
            fontweight="bold",
            color=title_color,
            va="top",
            transform=legend_axis.transAxes,
        )
        y -= section_step
        for row, wrapped in entries:
            cluster = int(row.cluster)
            text_color = (
                HIDDEN_COLOR
                if category == "hidden_novel"
                else PRIOR_COLOR
                if category == "prior_work"
                else UNCLEAR_COLOR
                if category == "unclear"
                else "#374151"
            )
            legend_axis.scatter(
                [0.018],
                [y - 0.004],
                s=32,
                color=[colors[cluster]],
                edgecolors=[colors[cluster]],
                linewidths=0.7,
                transform=legend_axis.transAxes,
                clip_on=True,
            )
            legend_axis.text(
                0.052,
                y,
                "\n".join(wrapped),
                fontsize=6.7,
                color=text_color,
                va="top",
                transform=legend_axis.transAxes,
                clip_on=True,
            )
            y -= line_step * max(1, len(wrapped)) + entry_gap
        y -= section_gap

    figure.savefig(output_path)
    plt.close(figure)


def process_run(
    run_dir: Path,
    prior_terms: dict[str, set[str]],
    direction_meta: dict[tuple[str, str], dict],
    alpha: float,
) -> list[dict]:
    tissue = tissue_from_run(run_dir)
    output_rows: list[dict] = []
    for location in ("FLT", "GC"):
        coords_path = run_dir / "plots" / f"{location}_umap_coordinates.tsv"
        clusters_path = run_dir / "clustering" / f"{location}_gene_clusters.tsv"
        annotations_path = run_dir / "plots" / f"{location}_cluster_annotations.tsv"
        if not coords_path.exists() or not clusters_path.exists() or not annotations_path.exists():
            continue
        coords = pd.read_csv(coords_path, sep="\t")
        clusters = pd.read_csv(clusters_path, sep="\t")
        labels = clusters["consensus"].astype(int).to_numpy()
        annotations = pd.read_csv(annotations_path, sep="\t")
        enrichment = load_enrichment_terms(run_dir, location, alpha)
        highlighted = classify_clusters(
            annotations,
            enrichment,
            tissue,
            prior_terms,
            direction_meta,
        )

        output_annotations = run_dir / "plots" / f"{location}_cluster_highlight_annotations.tsv"
        highlighted.to_csv(output_annotations, sep="\t", index=False)
        output_interpretation = run_dir / "plots" / f"{location}_cluster_interpretation.tsv"
        highlighted.to_csv(output_interpretation, sep="\t", index=False)

        output_plot = run_dir / "plots" / f"{location}_umap_prior_hidden_legend.png"
        plot_highlighted_umap(
            coords,
            labels,
            highlighted,
            f"{tissue} {run_dir.name} {location} GLARE clusters",
            output_plot,
        )

        counts = highlighted["highlight_category"].value_counts().to_dict()
        output_rows.append(
            {
                "tissue": tissue,
                "scope": run_dir.name,
                "location": location,
                "plot": str(output_plot),
                "annotations": str(output_annotations),
                "interpretation": str(output_interpretation),
                "hidden_novel_clusters": int(counts.get("hidden_novel", 0)),
                "prior_work_clusters": int(counts.get("prior_work", 0)),
                "other_clusters": int(counts.get("other", 0)),
                "unclear_clusters": int(counts.get("unclear", 0)),
            }
        )
    return output_rows


def main() -> None:
    args = parse_args()
    prior_terms = load_prior_terms(args.validation_dir)
    direction_meta = load_direction_meta(args.validation_dir)
    manifest_rows: list[dict] = []
    for tissue_dir in sorted(args.root.iterdir()):
        if not tissue_dir.is_dir():
            continue
        for scope in args.scopes:
            run_dir = tissue_dir / scope
            if run_dir.exists():
                manifest_rows.extend(
                    process_run(run_dir, prior_terms, direction_meta, args.alpha)
                )

    manifest = pd.DataFrame(manifest_rows)
    output_path = args.root / "priority_umap_highlight_manifest.tsv"
    manifest.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {len(manifest)} highlighted UMAP rows to {output_path}")


if __name__ == "__main__":
    main()
