#!/usr/bin/env python
"""Replot NASA OSDR expiMap accession heatmaps with readable row labels."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild all_program_accession_flt_minus_gc heatmaps from saved "
            "matrix TSVs using larger per-row strips and wider label margins."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("outputs"),
        help="Root directory to search for saved matrix TSVs.",
    )
    parser.add_argument(
        "--matrix-name",
        default="all_program_accession_flt_minus_gc_matrix_weighted_order.tsv",
        help="Matrix filename to plot inside each analysis directory.",
    )
    parser.add_argument(
        "--output-name",
        default="all_program_accession_flt_minus_gc_heatmap_weighted_order_all_labels.png",
        help="PNG filename to write next to each matrix.",
    )
    parser.add_argument(
        "--order",
        choices=("input", "abs", "signed"),
        default="input",
        help=(
            "Row order: keep the input TSV order, sort by absolute study-mean "
            "effect, or sort by signed study-mean effect from FLT-up to FLT-down."
        ),
    )
    parser.add_argument(
        "--ordered-matrix-name",
        default=None,
        help="Optional TSV filename for the row-ordered matrix.",
    )
    parser.add_argument(
        "--row-height",
        type=float,
        default=0.08,
        help="Figure inches allocated per pathway row.",
    )
    parser.add_argument(
        "--min-height",
        type=float,
        default=8.0,
        help="Minimum figure height in inches.",
    )
    parser.add_argument(
        "--max-height",
        type=float,
        default=140.0,
        help="Maximum figure height in inches.",
    )
    parser.add_argument(
        "--label-font-size",
        type=float,
        default=4.0,
        help="Y-axis pathway label font size.",
    )
    parser.add_argument(
        "--x-label-font-size",
        type=float,
        default=7.0,
        help="X-axis accession label font size.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=180,
        help="Output PNG DPI.",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=2.5,
        help="Symmetric color scale limit.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process the first N matrices; useful for smoke tests.",
    )
    return parser.parse_args()


def order_label(order: str) -> str:
    if order == "signed":
        return "signed study-mean FLT-GC effect (FLT-up to FLT-down)"
    if order == "abs":
        return "abs(study-mean FLT-GC effect)"
    return "input matrix order"


def order_title(order: str) -> str:
    if order == "signed":
        return "signed-effect row order"
    if order == "abs":
        return "absolute-effect row order"
    return "input row order"


def order_matrix(matrix, values, order: str):
    if order == "input":
        return matrix
    effect = values.mean(axis=1, skipna=True)
    ordered = matrix.assign(_study_mean_effect=effect)
    if order == "signed":
        ordered = ordered.sort_values(
            "_study_mean_effect",
            ascending=False,
            na_position="last",
            kind="mergesort",
        )
    elif order == "abs":
        ordered = (
            ordered.assign(_abs_study_mean_effect=ordered["_study_mean_effect"].abs())
            .sort_values(
                ["_abs_study_mean_effect", "_study_mean_effect"],
                ascending=[False, False],
                na_position="last",
                kind="mergesort",
            )
            .drop(columns=["_abs_study_mean_effect"])
        )
    return ordered.drop(columns=["_study_mean_effect"])


def figure_size(
    n_rows: int,
    n_cols: int,
    row_height: float,
    min_height: float,
    max_height: float,
) -> tuple[float, float]:
    height = min(max(min_height, n_rows * row_height + 2.0), max_height)
    width = max(18.0, min(32.0, 15.0 + 1.0 * n_cols))
    return width, height


def replot_matrix(matrix_path: Path, args: argparse.Namespace) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    matrix = pd.read_csv(matrix_path, sep="\t")
    if matrix.empty or "term" not in matrix.columns:
        raise ValueError(f"{matrix_path} is empty or lacks a 'term' column")

    values = matrix.drop(columns=["term"]).apply(pd.to_numeric, errors="coerce")
    matrix = order_matrix(matrix, values, args.order)
    if args.ordered_matrix_name:
        matrix.to_csv(matrix_path.with_name(args.ordered_matrix_name), sep="\t", index=False)

    labels = matrix["term"].astype(str).tolist()
    values = matrix.drop(columns=["term"]).apply(pd.to_numeric, errors="coerce")
    accessions = values.columns.astype(str).tolist()
    data = values.to_numpy(dtype=float)

    fig_width, fig_height = figure_size(
        n_rows=len(labels),
        n_cols=len(accessions),
        row_height=args.row_height,
        min_height=args.min_height,
        max_height=args.max_height,
    )

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#f2f2f2")
    image = ax.imshow(
        np.ma.masked_invalid(data),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=-args.vmax,
        vmax=args.vmax,
    )

    ax.set_title(
        f"All expiMap FLT-GC pathway shifts, {order_title(args.order)}",
        fontsize=10,
        pad=10,
    )
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=args.label_font_size)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.set_xticks(np.arange(len(accessions)))
    ax.set_xticklabels(
        accessions,
        rotation=90,
        ha="center",
        va="top",
        fontsize=args.x_label_font_size,
    )
    ax.set_xlabel("OSD accession / study", fontsize=8)

    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.015)
    cbar.set_label("mean FLT - mean GC pathway score", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.text(
        0.63,
        0.012,
        f"Rows ordered by {order_label(args.order)}",
        ha="center",
        va="bottom",
        fontsize=6,
    )
    fig.subplots_adjust(left=0.56, right=0.89, top=0.985, bottom=0.055)

    output_path = matrix_path.with_name(args.output_name)
    fig.savefig(output_path, dpi=args.dpi)
    plt.close(fig)
    return output_path


def main() -> None:
    args = parse_args()
    matrices = sorted(args.root.rglob(args.matrix_name))
    if args.limit is not None:
        matrices = matrices[: args.limit]
    if not matrices:
        raise SystemExit(f"No {args.matrix_name} files found under {args.root}")

    for index, matrix_path in enumerate(matrices, start=1):
        output_path = replot_matrix(matrix_path, args)
        print(f"[{index}/{len(matrices)}] wrote {output_path}", flush=True)


if __name__ == "__main__":
    main()
