"""Prepare expiMap inputs following the scArches tutorial-style HVG workflow."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from .io import require_import


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def matrix_sum(values: Any, axis: int) -> np.ndarray:
    summed = values.sum(axis=axis)
    if hasattr(summed, "A1"):
        return np.asarray(summed.A1)
    return np.asarray(summed).ravel()


def subset_terms(adata, selected_terms: np.ndarray):
    adata.uns["terms"] = np.asarray(adata.uns["terms"], dtype=object)[selected_terms].tolist()
    descriptions = np.asarray(
        adata.uns.get("term_descriptions", adata.uns["terms"]),
        dtype=object,
    )
    if len(descriptions) == len(selected_terms):
        adata.uns["term_descriptions"] = descriptions[selected_terms].tolist()
    else:
        adata.uns["term_descriptions"] = list(map(str, adata.uns["terms"]))
    adata.varm["I"] = np.asarray(adata.varm["I"])[:, selected_terms]


def infer_tissue(reference_path: Path, query_path: Path) -> str:
    reference_name = reference_path.name
    query_name = query_path.name
    reference_prefix = "archs4_mouse_"
    reference_suffix = "_reference_reactome_raw_counts.h5ad"
    query_prefix = "osdr_"
    query_suffix = "_flt_gc_reactome_raw_counts.h5ad"

    if reference_name.startswith(reference_prefix) and reference_name.endswith(reference_suffix):
        return reference_name[len(reference_prefix) : -len(reference_suffix)]
    if query_name.startswith(query_prefix) and query_name.endswith(query_suffix):
        return query_name[len(query_prefix) : -len(query_suffix)]
    raise SystemExit(
        "Could not infer tissue from input filenames. Expected names like "
        "archs4_mouse_<tissue>_reference_reactome_raw_counts.h5ad or "
        "osdr_<tissue>_flt_gc_reactome_raw_counts.h5ad."
    )


def prepare(args: argparse.Namespace) -> Path:
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    sc = require_import("scanpy", "pip install -r requirements-nasa-mouse-glare.txt")

    reference_path = Path(args.reference_h5ad)
    query_path = Path(args.query_h5ad)
    tissue = args.label or infer_tissue(reference_path, query_path)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log(f"Loading reference {reference_path}")
    reference = ad.read_h5ad(reference_path)
    log(f"Loading query {query_path}")
    query = ad.read_h5ad(query_path)

    if "counts" not in reference.layers or "counts" not in query.layers:
        raise SystemExit("Both reference and query inputs must contain layers['counts'].")
    if "I" not in reference.varm or "I" not in query.varm:
        raise SystemExit("Both reference and query inputs must contain varm['I'].")
    if "terms" not in reference.uns or "terms" not in query.uns:
        raise SystemExit("Both reference and query inputs must contain uns['terms'].")

    common_genes = reference.var_names.intersection(query.var_names)
    common_genes = reference.var_names[reference.var_names.isin(common_genes)]
    if len(common_genes) < args.n_top_genes:
        raise SystemExit(
            f"Only {len(common_genes)} common genes are available, fewer than "
            f"requested n_top_genes={args.n_top_genes}."
        )
    reference = reference[:, common_genes].copy()
    query = query[:, common_genes].copy()

    hvg_adata = reference.copy()
    hvg_adata.X = hvg_adata.layers["counts"].copy()
    sc.pp.normalize_total(hvg_adata)
    sc.pp.log1p(hvg_adata)

    batch_key = args.batch_key if args.batch_key in hvg_adata.obs else None
    hvg_error = None
    log(
        "Selecting HVGs with scanpy.pp.highly_variable_genes "
        f"n_top_genes={args.n_top_genes} batch_key={batch_key!r}"
    )
    try:
        sc.pp.highly_variable_genes(
            hvg_adata,
            n_top_genes=args.n_top_genes,
            batch_key=batch_key,
            subset=False,
        )
        hvg_method = "scanpy_highly_variable_genes_batch"
    except Exception as exc:  # pragma: no cover - only used for runtime robustness.
        if not args.allow_no_batch_fallback or batch_key is None:
            raise
        hvg_error = repr(exc)
        log(f"Batch-aware HVG failed; falling back to unbatched HVG: {hvg_error}")
        sc.pp.highly_variable_genes(
            hvg_adata,
            n_top_genes=args.n_top_genes,
            batch_key=None,
            subset=False,
        )
        hvg_method = "scanpy_highly_variable_genes_unbatched_fallback"

    hvg = hvg_adata.var["highly_variable"].to_numpy(dtype=bool)
    if hvg.sum() == 0:
        raise SystemExit("HVG selection returned zero genes.")
    hvg_genes = hvg_adata.var_names[hvg]
    reference = reference[:, hvg_genes].copy()
    query = query[:, hvg_genes].copy()

    selected_terms = matrix_sum(reference.varm["I"], axis=0) > args.min_genes_per_term
    if selected_terms.sum() == 0:
        raise SystemExit("Term filtering after HVG selection returned zero terms.")
    subset_terms(reference, selected_terms)
    subset_terms(query, selected_terms)

    gene_keep = matrix_sum(reference.varm["I"], axis=1) > 0
    reference = reference[:, gene_keep].copy()
    query = query[:, gene_keep].copy()

    reference.X = reference.layers["counts"].copy()
    query.X = query.layers["counts"].copy()
    reference.uns.setdefault("expimap_preprocessing", {})
    query.uns.setdefault("expimap_preprocessing", {})
    for adata in (reference, query):
        adata.uns["expimap_preprocessing"].update(
            {
                "tutorial_style_hvg": True,
                "hvg_method": hvg_method,
                "hvg_batch_key": batch_key,
                "n_top_genes_requested": int(args.n_top_genes),
                "min_genes_per_term_strict_gt": int(args.min_genes_per_term),
                "recommended_recon_loss": "nb",
            }
        )

    reference_out = output_dir / f"archs4_mouse_{tissue}_reference_tutorial_hvg_raw_counts.h5ad"
    query_out = output_dir / f"osdr_{tissue}_query_tutorial_hvg_raw_counts.h5ad"
    reference.write_h5ad(reference_out)
    query.write_h5ad(query_out)

    hvg_table = pd.DataFrame(
        {
            "gene_id": hvg_adata.var_names.astype(str),
            "highly_variable": hvg_adata.var["highly_variable"].to_numpy(dtype=bool),
        }
    )
    for column in [
        "means",
        "dispersions",
        "dispersions_norm",
        "highly_variable_nbatches",
        "highly_variable_intersection",
    ]:
        if column in hvg_adata.var:
            hvg_table[column] = hvg_adata.var[column].to_numpy()
    hvg_table.to_csv(output_dir / "hvg_selection.tsv", sep="\t", index=False)

    terms = list(map(str, reference.uns["terms"]))
    pd.DataFrame(
        {
            "term": terms,
            "description": list(map(str, reference.uns.get("term_descriptions", terms))),
            "n_genes": matrix_sum(reference.varm["I"], axis=0).astype(int),
        }
    ).to_csv(output_dir / "tutorial_hvg_terms.tsv", sep="\t", index=False)

    summary = {
        "reference_input": str(reference_path),
        "query_input": str(query_path),
        "reference_output": str(reference_out),
        "query_output": str(query_out),
        "n_reference_samples": int(reference.n_obs),
        "n_query_samples": int(query.n_obs),
        "n_common_genes_before_hvg": int(len(common_genes)),
        "n_top_genes_requested": int(args.n_top_genes),
        "n_hvg_selected": int(len(hvg_genes)),
        "n_genes_after_term_filter": int(reference.n_vars),
        "n_terms_after_hvg_filter": int(len(reference.uns["terms"])),
        "min_genes_per_term_strict_gt": int(args.min_genes_per_term),
        "hvg_method": hvg_method,
        "hvg_batch_key": batch_key,
        "hvg_error": hvg_error,
        "outputs": {
            "reference_h5ad": str(reference_out),
            "query_h5ad": str(query_out),
            "hvg_selection": str(output_dir / "hvg_selection.tsv"),
            "terms": str(output_dir / "tutorial_hvg_terms.tsv"),
            "summary": str(output_dir / "tutorial_hvg_input_manifest.json"),
        },
    }
    (output_dir / "tutorial_hvg_input_manifest.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    log(json.dumps(summary, indent=2))
    return output_dir / "tutorial_hvg_input_manifest.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create matched ARCHS4 reference and OSDR query h5ads using the "
            "scArches expiMap tutorial-style HVG and term filtering workflow."
        )
    )
    parser.add_argument("--reference-h5ad", required=True)
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--label",
        help=(
            "Output label for the generated reference/query filenames. "
            "Defaults to inferring the broad tissue from input filenames."
        ),
    )
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--batch-key", default="archs4_condition")
    parser.add_argument("--min-genes-per-term", type=int, default=12)
    parser.add_argument(
        "--allow-no-batch-fallback",
        action="store_true",
        help="Fall back to unbatched HVG selection if batch-aware HVG fails.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    prepare(parse_args())
