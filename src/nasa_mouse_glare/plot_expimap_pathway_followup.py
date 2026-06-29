"""Create pathway-level expiMap follow-up plots for accession-aware hits."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

from .io import require_import


def safe_name(value: str, max_len: int = 120) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")
    return cleaned[:max_len]


def read_tables(base_dir: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scores = pd.read_csv(base_dir / "query_pathway_scores.tsv", sep="\t")
    validation = base_dir / "accession_validation"
    meta = pd.read_csv(validation / "random_effects_meta_analysis.tsv", sep="\t")
    effects = pd.read_csv(validation / "per_accession_effects.tsv", sep="\t")
    loo = pd.read_csv(validation / "leave_one_accession_out.tsv", sep="\t")
    loo_summary = pd.read_csv(validation / "leave_one_out_summary.tsv", sep="\t")
    return scores, meta, effects, loo, loo_summary


def select_terms(meta, loo_summary, top_n: int):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    merged = meta.merge(loo_summary, on="term", how="left")
    merged["all_accessions_same_direction"] = (
        merged["n_accession_same_direction"] == merged["n_accessions"]
    )
    merged["all_loo_same_direction"] = (
        merged["n_same_direction"] == merged["n_leave_one_out"]
    )
    merged["loo_fdr_stable"] = merged["maximum_leave_one_out_fdr"] < 0.05
    merged["robust_followup"] = (
        (merged["meta_fdr"] < 0.05)
        & merged["all_accessions_same_direction"]
        & merged["all_loo_same_direction"]
        & merged["loo_fdr_stable"]
    )
    robust = merged.loc[merged["robust_followup"]].sort_values(
        ["meta_fdr", "maximum_leave_one_out_fdr", "i2"],
        kind="stable",
    )
    if len(robust) >= top_n:
        selected = robust.head(top_n).copy()
    else:
        fillers = merged.loc[~merged["term"].isin(robust["term"])].sort_values(
            ["meta_fdr", "maximum_leave_one_out_fdr", "i2"],
            kind="stable",
        )
        selected = pd.concat([robust, fillers.head(top_n - len(robust))], ignore_index=True)
    return merged, selected


def plot_scores(scores, term: str, output_path: Path):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")

    accessions = sorted(scores["id.accession"].astype(str).unique())
    colors = {"flight": "#c43c39", "ground_control": "#2878b5"}
    offsets = {"ground_control": -0.18, "flight": 0.18}
    fig_width = max(7, len(accessions) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 4.8))
    rng = np.random.default_rng(2020)
    for index, accession in enumerate(accessions):
        accession_frame = scores.loc[scores["id.accession"].astype(str).eq(accession)]
        for condition in ["ground_control", "flight"]:
            values = accession_frame.loc[
                accession_frame["condition_inferred"].eq(condition),
                term,
            ].astype(float).to_numpy()
            if len(values) == 0:
                continue
            x = index + offsets[condition] + rng.normal(0.0, 0.035, size=len(values))
            ax.scatter(
                x,
                values,
                s=20,
                alpha=0.7,
                color=colors[condition],
                label=condition if index == 0 else None,
                edgecolors="none",
            )
            mean = float(np.mean(values))
            ax.plot(
                [index + offsets[condition] - 0.12, index + offsets[condition] + 0.12],
                [mean, mean],
                color=colors[condition],
                linewidth=2,
            )
    ax.axhline(0.0, color="#999999", linewidth=0.8, alpha=0.6)
    ax.set_xticks(range(len(accessions)))
    ax.set_xticklabels(accessions, rotation=45, ha="right")
    ax.set_ylabel("expiMap pathway score")
    ax.set_title(term.replace("_", " "))
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_forest(effects, meta_row, term: str, output_path: Path):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")

    frame = effects.loc[effects["term"].eq(term)].copy()
    frame = frame.sort_values("flight_minus_ground", kind="stable")
    y = np.arange(len(frame))
    effect = frame["flight_minus_ground"].astype(float).to_numpy()
    se = np.sqrt(frame["effect_variance"].astype(float).to_numpy())
    fig, ax = plt.subplots(figsize=(7, max(3.8, len(frame) * 0.42)))
    ax.errorbar(
        effect,
        y,
        xerr=1.96 * se,
        fmt="o",
        color="#333333",
        ecolor="#888888",
        capsize=2,
        markersize=4,
    )
    meta_effect = float(meta_row["meta_effect"])
    ax.axvline(0.0, color="#999999", linewidth=0.9)
    ax.axvline(meta_effect, color="#c43c39", linewidth=1.6)
    labels = [
        f"{row['id.accession']} ({int(row['n_flight'])}/{int(row['n_ground_control'])})"
        for _, row in frame.iterrows()
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Flight minus ground expiMap score")
    ax.set_title(
        f"{term.replace('_', ' ')}\n"
        f"meta effect={meta_effect:.4g}, FDR={float(meta_row['meta_fdr']):.3g}, I2={float(meta_row['i2']):.2g}"
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_loo(loo, meta_row, term: str, output_path: Path):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")

    frame = loo.loc[loo["term"].eq(term)].copy()
    frame = frame.sort_values("held_out_accession", kind="stable")
    x = np.arange(len(frame))
    colors = ["#2878b5" if value < 0.05 else "#c43c39" for value in frame["meta_fdr"]]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.scatter(x, frame["meta_effect"].astype(float), c=colors, s=45)
    ax.axhline(0.0, color="#999999", linewidth=0.9)
    ax.axhline(float(meta_row["meta_effect"]), color="#222222", linewidth=1.2)
    for xpos, (_, row) in zip(x, frame.iterrows()):
        ax.annotate(
            f"{float(row['meta_fdr']):.2g}",
            (xpos, float(row["meta_effect"])),
            textcoords="offset points",
            xytext=(0, 7),
            ha="center",
            fontsize=8,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(frame["held_out_accession"].astype(str), rotation=45, ha="right")
    ax.set_ylabel("LOO random-effects meta effect")
    ax.set_title(f"{term.replace('_', ' ')}\nlabels are leave-one-out FDR values")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_readme(tissue: str, selected, output_dir: Path, base_dir: Path) -> Path:
    lines = [
        f"# {tissue} expiMap Pathway Follow-up",
        "",
        f"- Source scores: `{base_dir / 'query_pathway_scores.tsv'}`",
        f"- Selected pathways: {len(selected)}",
        "- Robust flag: meta FDR < 0.05, all accessions same direction, all leave-one-out fits same direction, and maximum leave-one-out FDR < 0.05.",
        "",
        "| term | meta effect | meta FDR | I2 | same accessions | max LOO FDR | robust |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in selected.iterrows():
        lines.append(
            f"| {row.term} | {row.meta_effect:.4g} | {row.meta_fdr:.4g} | {row.i2:.3g} | "
            f"{int(row.n_accession_same_direction)}/{int(row.n_accessions)} | "
            f"{row.maximum_leave_one_out_fdr:.4g} | {bool(row.robust_followup)} |"
        )
    path = output_dir / "README.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run(args) -> None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_root = Path(args.output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for tissue in args.tissue:
        base_dir = Path(
            f"outputs/expimap_archs4_reference_osdr_query_{tissue}/query_nb_allref_50epoch"
        )
        tissue_dir = output_root / tissue
        plot_dir = tissue_dir / "plots"
        plot_dir.mkdir(parents=True, exist_ok=True)
        scores, meta, effects, loo, loo_summary = read_tables(base_dir)
        merged, selected = select_terms(meta, loo_summary, args.top_n)
        merged.to_csv(tissue_dir / "all_pathway_loo_merged.tsv", sep="\t", index=False)
        selected.to_csv(tissue_dir / "selected_pathway_followup.tsv", sep="\t", index=False)
        for _, row in selected.iterrows():
            term = str(row["term"])
            prefix = safe_name(term)
            plot_scores(scores, term, plot_dir / f"{prefix}_scores_by_accession.png")
            plot_forest(effects, row, term, plot_dir / f"{prefix}_accession_forest.png")
            plot_loo(loo, row, term, plot_dir / f"{prefix}_loo_effects.png")
        write_readme(tissue, selected, tissue_dir, base_dir)
        summary = selected.copy()
        summary.insert(0, "tissue", tissue)
        all_rows.append(summary)
    combined = pd.concat(all_rows, ignore_index=True)
    combined.to_csv(output_root / "selected_pathway_followup_all_tissues.tsv", sep="\t", index=False)
    write_readme("combined", combined, output_root, Path("outputs"))
    print(output_root / "selected_pathway_followup_all_tissues.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot selected expiMap pathway-level accession and LOO diagnostics."
    )
    parser.add_argument(
        "--tissue",
        action="append",
        required=True,
        help="Tissue name, e.g. thymus. Repeat for multiple tissues.",
    )
    parser.add_argument("--top-n", type=int, default=8)
    parser.add_argument("--output-dir", default="outputs/expimap_pathway_followup_allref")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
