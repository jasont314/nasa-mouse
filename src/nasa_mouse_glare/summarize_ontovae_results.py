"""Summarize OntoVAE tissue runs and compare them with available expiMap runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import require_import


DEFAULT_OUTPUT_DIR = Path("outputs/ontovae_pipeline/summary")
DEFAULT_SEARCH_ROOT = Path("outputs")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_output_dir(path: Path) -> dict:
    parts = path.parts
    tissue = ""
    group = ""
    if len(parts) >= 2 and parts[0] == "outputs":
        root = parts[1]
        if root == "ontovae_skeletal_muscle_splits" and len(parts) >= 3:
            tissue = "skeletal_muscle"
            group = parts[2]
        elif root.startswith("ontovae_"):
            tissue = root.removeprefix("ontovae_")
    return {"tissue": tissue, "group": group}


def significant_loo(meta, loo):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    if meta.empty or loo.empty:
        return pd.DataFrame()
    merged = meta.merge(loo, on="term", how="inner")
    if merged.empty:
        return merged
    return merged.loc[
        (merged["meta_fdr"] < 0.05)
        & (merged["maximum_leave_one_out_fdr"] < 0.05)
        & (merged["n_same_direction"] == merged["n_leave_one_out"])
    ].copy()


def summarize_score_set(
    *,
    run_dir: Path,
    score_set: str,
    analysis_dir: Path,
    validation_dir: Path,
    training: dict,
) -> tuple[dict, list[dict]]:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    comparison_path = analysis_dir / "flt_vs_gc_pathway_comparison.tsv"
    meta_path = validation_dir / "random_effects_meta_analysis.tsv"
    loo_path = validation_dir / "leave_one_out_summary.tsv"
    terms_path = run_dir / "terms.tsv"

    comparison = (
        pd.read_csv(comparison_path, sep="\t") if comparison_path.exists() else pd.DataFrame()
    )
    meta = pd.read_csv(meta_path, sep="\t") if meta_path.exists() else pd.DataFrame()
    loo = pd.read_csv(loo_path, sep="\t") if loo_path.exists() else pd.DataFrame()
    terms = pd.read_csv(terms_path, sep="\t") if terms_path.exists() else pd.DataFrame()
    desc = dict(zip(terms.get("term", []), terms.get("description", [])))

    run_info = parse_output_dir(run_dir)
    counts = training.get("counts", {})
    torch = training.get("torch", {})
    row = {
        **run_info,
        "mode": training.get("mode", run_dir.name),
        "score_set": score_set,
        "output_dir": str(run_dir),
        "query_samples": counts.get("query_samples", ""),
        "reference_samples": counts.get("reference_samples", ""),
        "genes": counts.get("genes", ""),
        "terms": counts.get("terms", ""),
        "model_device": torch.get("model_device", ""),
        "cuda_device_name": torch.get("cuda_device_name", ""),
        "ordinary_welch_fdr_lt_005": 0,
        "ordinary_mannwhitney_fdr_lt_005": 0,
        "random_effects_fdr_lt_005": 0,
        "loo_stable_fdr_lt_005": 0,
        "min_welch_fdr": "",
        "top_welch_term": "",
        "top_welch_description": "",
        "top_welch_flight_minus_ground": "",
        "min_meta_fdr": "",
        "top_meta_term": "",
        "top_meta_description": "",
        "top_meta_effect": "",
        "top_meta_max_loo_fdr": "",
    }

    top_rows: list[dict] = []
    if not comparison.empty:
        comparison = comparison.copy()
        comparison["description"] = comparison["term"].map(desc).fillna("")
        row["ordinary_welch_fdr_lt_005"] = int((comparison["welch_fdr"] < 0.05).sum())
        row["ordinary_mannwhitney_fdr_lt_005"] = int(
            (comparison["mannwhitney_fdr"] < 0.05).sum()
        )
        top = comparison.sort_values("welch_fdr").head(1).iloc[0]
        row["min_welch_fdr"] = float(top["welch_fdr"])
        row["top_welch_term"] = str(top["term"])
        row["top_welch_description"] = str(top.get("description", ""))
        row["top_welch_flight_minus_ground"] = float(top["flight_minus_ground"])
        for _, item in comparison.sort_values("welch_fdr").head(10).iterrows():
            top_rows.append(
                {
                    **run_info,
                    "mode": row["mode"],
                    "score_set": score_set,
                    "rank_source": "ordinary_welch",
                    "term": item["term"],
                    "description": item.get("description", ""),
                    "effect": item["flight_minus_ground"],
                    "fdr": item["welch_fdr"],
                    "output_dir": str(run_dir),
                }
            )

    if not meta.empty:
        meta = meta.copy()
        meta["description"] = meta["term"].map(desc).fillna("")
        row["random_effects_fdr_lt_005"] = int((meta["meta_fdr"] < 0.05).sum())
        top = meta.sort_values("meta_fdr").head(1).iloc[0]
        row["min_meta_fdr"] = float(top["meta_fdr"])
        row["top_meta_term"] = str(top["term"])
        row["top_meta_description"] = str(top.get("description", ""))
        row["top_meta_effect"] = float(top["meta_effect"])
        if not loo.empty:
            match = loo.loc[loo["term"].eq(top["term"])]
            if not match.empty:
                row["top_meta_max_loo_fdr"] = float(
                    match.iloc[0]["maximum_leave_one_out_fdr"]
                )
        stable = significant_loo(meta, loo)
        row["loo_stable_fdr_lt_005"] = int(len(stable))
        for _, item in meta.sort_values("meta_fdr").head(10).iterrows():
            top_rows.append(
                {
                    **run_info,
                    "mode": row["mode"],
                    "score_set": score_set,
                    "rank_source": "random_effects",
                    "term": item["term"],
                    "description": item.get("description", ""),
                    "effect": item["meta_effect"],
                    "fdr": item["meta_fdr"],
                    "output_dir": str(run_dir),
                }
            )
    return row, top_rows


def summarize_run(run_dir: Path) -> tuple[list[dict], list[dict]]:
    summary_path = run_dir / "training_summary.json"
    if not summary_path.exists():
        return [], []
    training = load_json(summary_path)
    rows = []
    top_rows = []
    row, top = summarize_score_set(
        run_dir=run_dir,
        score_set="finetuned_or_direct",
        analysis_dir=run_dir / "analysis",
        validation_dir=run_dir / "accession_validation",
        training=training,
    )
    rows.append(row)
    top_rows.extend(top)
    if (run_dir / "pretrained_query_pathway_scores.tsv").exists():
        row, top = summarize_score_set(
            run_dir=run_dir,
            score_set="pre_finetune_projection",
            analysis_dir=run_dir / "pretrained_query_analysis",
            validation_dir=run_dir / "pretrained_query_accession_validation",
            training=training,
        )
        rows.append(row)
        top_rows.extend(top)
    return rows, top_rows


def discover_ontovae_dirs(root: Path) -> list[Path]:
    return sorted(
        path.parent
        for path in root.glob("ontovae_*/**/training_summary.json")
        if "ontovae_smoke" not in str(path) and "ontovae_pipeline" not in str(path)
    )


def existing_expimap_validation_dirs(tissue: str, group: str, mode: str) -> list[Path]:
    if group:
        root = Path("outputs/expimap_muscle_targeted_combined_min8")
        if mode == "direct_osdr":
            return [root / f"direct_{group}_nb_100epoch" / "accession_validation"]
        return [root / f"query_{group}_nb_allref_50epoch" / "accession_validation"]
    if mode == "direct_osdr":
        return [
            Path(f"outputs/expimap_direct_osdr_{tissue}")
            / "raw_counts_nb_50epoch"
            / "accession_validation"
        ]
    root = Path(f"outputs/expimap_archs4_reference_osdr_query_{tissue}")
    return [
        root / "query_nb_allref_50epoch" / "posterior_mean_accession_validation",
        root / "query_nb_allref_50epoch" / "accession_validation",
        root / "query_nb_5000stratified_seed2020_50epoch" / "posterior_mean_accession_validation",
        root / "query_nb_5000stratified_seed2020_50epoch" / "accession_validation",
        root / "query_nb_1000ref_50epoch" / "accession_validation",
    ]


def summarize_expimap_counterparts(onto_rows: list[dict]) -> list[dict]:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    rows = []
    seen: set[tuple[str, str, str]] = set()
    for onto in onto_rows:
        if onto["score_set"] != "finetuned_or_direct":
            continue
        key = (onto["tissue"], onto["group"], onto["mode"])
        if key in seen:
            continue
        seen.add(key)
        for validation_dir in existing_expimap_validation_dirs(*key):
            meta_path = validation_dir / "random_effects_meta_analysis.tsv"
            loo_path = validation_dir / "leave_one_out_summary.tsv"
            if not meta_path.exists():
                continue
            meta = pd.read_csv(meta_path, sep="\t")
            loo = pd.read_csv(loo_path, sep="\t") if loo_path.exists() else pd.DataFrame()
            stable = significant_loo(meta, loo)
            top = meta.sort_values("meta_fdr").head(1)
            rows.append(
                {
                    "tissue": key[0],
                    "group": key[1],
                    "mode": key[2],
                    "method": "expiMap",
                    "validation_dir": str(validation_dir),
                    "random_effects_fdr_lt_005": int((meta["meta_fdr"] < 0.05).sum()),
                    "loo_stable_fdr_lt_005": int(len(stable)),
                    "top_meta_term": "" if top.empty else str(top.iloc[0]["term"]),
                    "min_meta_fdr": "" if top.empty else float(top.iloc[0]["meta_fdr"]),
                    "top_meta_effect": ""
                    if top.empty
                    else float(top.iloc[0]["meta_effect"]),
                }
            )
            break
    return rows


def write_top_gene_table(top_rows: list[dict], output_dir: Path) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    selected = {
        (row["output_dir"], row["term"])
        for row in top_rows
        if row["rank_source"] in {"ordinary_welch", "random_effects"}
    }
    rows = []
    for run_output_dir, term in sorted(selected):
        weights_path = Path(run_output_dir) / "term_gene_weights_top.tsv"
        if not weights_path.exists():
            continue
        weights = pd.read_csv(weights_path, sep="\t")
        subset = weights.loc[weights["term"].eq(term)].head(10)
        for _, item in subset.iterrows():
            rows.append(
                {
                    "output_dir": run_output_dir,
                    "term": term,
                    "rank": int(item["rank"]),
                    "gene": item["gene"],
                    "decoder_weight": float(item["decoder_weight"]),
                }
            )
    path = output_dir / "ontovae_top_gene_weights_for_top_terms.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def write_readme(output_dir: Path, summary, top_terms, expimap) -> Path:
    if summary.empty:
        readme = output_dir / "README.md"
        readme.write_text(
            "\n".join(
                [
                    "# OntoVAE Summary",
                    "",
                    "No completed non-smoke OntoVAE runs were found yet.",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        return readme
    significant = summary.loc[
        (summary["score_set"].eq("finetuned_or_direct"))
        & (
            (summary["ordinary_welch_fdr_lt_005"] > 0)
            | (summary["random_effects_fdr_lt_005"] > 0)
            | (summary["loo_stable_fdr_lt_005"] > 0)
        )
    ]
    lines = [
        "# OntoVAE Summary",
        "",
        f"- OntoVAE score sets summarized: {len(summary)}",
        f"- Fine-tuned/direct runs with any ordinary or accession-aware FDR < 0.05: {len(significant)}",
        f"- Runs with strict LOO-stable FDR < 0.05: {int((summary['loo_stable_fdr_lt_005'] > 0).sum())}",
        "",
        "Primary files:",
        "",
        "- `ontovae_run_summary.tsv`",
        "- `ontovae_top_terms.tsv`",
        "- `ontovae_top_gene_weights_for_top_terms.tsv`",
        "- `ontovae_vs_expimap_random_effects.tsv`",
    ]
    if not significant.empty:
        lines.extend(["", "Top OntoVAE rows:", "", "| tissue | group | mode | score set | ordinary FDR hits | random-effects FDR hits | LOO-stable hits | top meta term |"])
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for _, row in significant.head(12).iterrows():
            lines.append(
                f"| {row['tissue']} | {row['group']} | {row['mode']} | {row['score_set']} | "
                f"{row['ordinary_welch_fdr_lt_005']} | {row['random_effects_fdr_lt_005']} | "
                f"{row['loo_stable_fdr_lt_005']} | {row['top_meta_term']} |"
            )
    readme = output_dir / "README.md"
    readme.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return readme


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_dirs = discover_ontovae_dirs(Path(args.search_root))
    summary_rows: list[dict] = []
    top_rows: list[dict] = []
    for run_dir in run_dirs:
        rows, top = summarize_run(run_dir)
        summary_rows.extend(rows)
        top_rows.extend(top)

    summary = pd.DataFrame(summary_rows)
    top_terms = pd.DataFrame(top_rows)
    expimap = pd.DataFrame(summarize_expimap_counterparts(summary_rows))

    summary_path = output_dir / "ontovae_run_summary.tsv"
    top_path = output_dir / "ontovae_top_terms.tsv"
    expimap_path = output_dir / "ontovae_vs_expimap_random_effects.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    top_terms.to_csv(top_path, sep="\t", index=False)
    expimap.to_csv(expimap_path, sep="\t", index=False)
    write_top_gene_table(top_rows, output_dir)
    write_readme(output_dir, summary, top_terms, expimap)
    print(
        json.dumps(
            {
                "runs": int(len(summary)),
                "top_terms": int(len(top_terms)),
                "summary": str(summary_path),
                "top_terms_path": str(top_path),
                "expimap_comparison": str(expimap_path),
            },
            indent=2,
        )
    )
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize OntoVAE tissue outputs.")
    parser.add_argument("--search-root", default=str(DEFAULT_SEARCH_ROOT))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
