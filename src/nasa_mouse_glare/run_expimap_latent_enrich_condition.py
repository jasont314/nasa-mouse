"""Run expiMap latent_enrich condition Bayes-factor-style tests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import require_import


def bf_table(adata, key: str, terms: list[str]):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    scores = adata.uns[key]
    for group, values in scores.items():
        for term, p_h0, p_h1, bf in zip(
            terms,
            values["p_h0"],
            values["p_h1"],
            values["bf"],
            strict=False,
        ):
            rows.append(
                {
                    "term": str(term),
                    "group": str(group),
                    "p_h0": float(p_h0),
                    "p_h1": float(p_h1),
                    "bf": float(bf),
                }
            )
    return pd.DataFrame(rows)


def run(args) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    sca = require_import("scarches", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    np.random.seed(args.seed)

    adata = ad.read_h5ad(args.h5ad)
    if args.groups not in adata.obs:
        raise SystemExit(f"groups column {args.groups!r} not found in AnnData obs.")
    if args.comparison not in set(adata.obs[args.groups].astype(str)):
        raise SystemExit(
            f"comparison {args.comparison!r} not found in {args.groups!r} values."
        )

    model = sca.models.EXPIMAP.load(args.model_dir, adata=adata)
    model.latent_directions(method=args.direction_method, adata=model.adata)
    model.latent_enrich(
        groups=args.groups,
        comparison=args.comparison,
        n_sample=args.n_sample,
        use_directions=not args.no_directions,
        directions_key="directions",
        adata=model.adata,
        exact=not args.approximate,
        key_added="condition_bf_scores",
    )

    terms = list(map(str, model.adata.uns["terms"]))
    table = bf_table(model.adata, "condition_bf_scores", terms)
    bayes_path = output_dir / "condition_latent_enrich_bayes_factors.tsv"
    table.to_csv(bayes_path, sep="\t", index=False)

    h5ad_path = output_dir / "model_adata_with_condition_bf.h5ad"
    model.adata.write_h5ad(h5ad_path)

    summary = {
        "model_dir": str(args.model_dir),
        "h5ad": str(args.h5ad),
        "groups": args.groups,
        "comparison": args.comparison,
        "use_directions": bool(not args.no_directions),
        "direction_method": args.direction_method,
        "n_sample": int(args.n_sample),
        "exact": bool(not args.approximate),
        "seed": int(args.seed),
        "n_terms": int(len(terms)),
        "table_shape": [int(table.shape[0]), int(table.shape[1])],
        "max_abs_bf": float(table["bf"].abs().max()) if not table.empty else None,
        "torch": {
            "version": str(torch.__version__),
            "cuda_available": bool(torch.cuda.is_available()),
            "cuda_device_count": int(torch.cuda.device_count()),
            "cuda_device_name": (
                torch.cuda.get_device_name(0) if torch.cuda.is_available() else ""
            ),
            "model_parameter_device": str(next(model.model.parameters()).device),
        },
        "outputs": {
            "bayes_factors": str(bayes_path),
            "h5ad": str(h5ad_path),
        },
    }
    summary_path = output_dir / "latent_enrich_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run expiMap latent_enrich condition Bayes-factor-style tests."
    )
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--groups", default="condition_inferred")
    parser.add_argument("--comparison", default="ground_control")
    parser.add_argument("--n-sample", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--direction-method", choices=["sum", "counts"], default="sum")
    parser.add_argument("--no-directions", action="store_true")
    parser.add_argument("--approximate", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
