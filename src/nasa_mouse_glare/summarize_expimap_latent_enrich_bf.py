"""Summarize expiMap latent_enrich Bayes-factor-style pathway calls."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import require_import


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(args.bayes_factors, sep="\t")
    required = {"term", "group", "p_h0", "p_h1", "bf"}
    missing = required.difference(table.columns)
    if missing:
        raise SystemExit(
            f"Bayes factor table is missing required columns: {sorted(missing)}"
        )

    table["bf"] = table["bf"].astype(float)
    table["abs_bf"] = table["bf"].abs()
    table["passes_bf_threshold"] = table["abs_bf"].ge(args.threshold)
    table = table.sort_values("abs_bf", ascending=False, kind="stable")

    ranked_path = output_dir / "latent_enrich_bf_ranked.tsv"
    calls_path = output_dir / "latent_enrich_bf_pass_threshold.tsv"
    summary_path = output_dir / "latent_enrich_bf_summary.json"
    readme_path = output_dir / "README.md"

    passing = table.loc[table["passes_bf_threshold"]].copy()
    table.to_csv(ranked_path, sep="\t", index=False)
    passing.to_csv(calls_path, sep="\t", index=False)

    counts_by_cutoff = {
        str(cutoff): int(table["abs_bf"].ge(cutoff).sum())
        for cutoff in args.report_cutoff
    }
    top = table.head(args.top_n).copy()
    summary = {
        "bayes_factors": str(args.bayes_factors),
        "threshold": args.threshold,
        "n_terms": int(len(table)),
        "n_passing_threshold": int(len(passing)),
        "max_abs_bf": float(table["abs_bf"].max()) if len(table) else None,
        "counts_by_abs_bf_cutoff": counts_by_cutoff,
        "outputs": {
            "ranked": str(ranked_path),
            "passing": str(calls_path),
            "readme": str(readme_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# expiMap latent_enrich Bayes Factors",
        "",
        "This summarizes the expiMap paper-style `latent_enrich` output. The decision rule here is based on absolute log Bayes-factor-style score, not FDR.",
        "",
        f"- Terms tested: {len(table)}",
        f"- Threshold: abs(bf) >= {args.threshold:g}",
        f"- Terms passing threshold: {len(passing)}",
        f"- Maximum abs(bf): {float(table['abs_bf'].max()) if len(table) else float('nan'):.4g}",
        "",
        "## Counts By Cutoff",
        "",
        "| abs(bf) cutoff | terms |",
        "| --- | ---: |",
    ]
    for cutoff, count in counts_by_cutoff.items():
        lines.append(f"| {cutoff} | {count} |")

    lines.extend(
        [
            "",
            f"## Top {len(top)} Terms By abs(bf)",
            "",
            "| term | group | p_h0 | p_h1 | bf | abs(bf) | passes threshold |",
            "| --- | --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in top.itertuples(index=False):
        lines.append(
            f"| {row.term} | {row.group} | {float(row.p_h0):.4g} | {float(row.p_h1):.4g} | {float(row.bf):.4g} | {float(row.abs_bf):.4g} | {bool(row.passes_bf_threshold)} |"
        )
    readme_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(readme_path)
    return readme_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize expiMap latent_enrich Bayes-factor-style results."
    )
    parser.add_argument("--bayes-factors", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threshold", type=float, default=2.3)
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument(
        "--report-cutoff",
        type=float,
        action="append",
        default=[0.5, 1.0, 2.3],
        help="Additional abs(bf) cutoffs to count. Can be passed multiple times.",
    )
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
