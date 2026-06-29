"""Export deterministic posterior-mean expiMap scores from a saved model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import require_import


def run(args) -> Path:
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    sca = require_import("scarches", "pip install -r requirements-nasa-mouse-glare.txt")

    adata = ad.read_h5ad(args.input)
    model = sca.models.EXPIMAP.load(args.model, adata=adata)
    latent = model.get_latent(mean=not args.sample_latent)
    if isinstance(latent, tuple):
        latent = latent[0]
    model.update_terms("terms")
    terms = list(map(str, model.adata.uns["terms"]))
    if latent.shape[1] != len(terms):
        raise SystemExit(
            f"Latent dimension ({latent.shape[1]}) does not match terms ({len(terms)})."
        )
    scores = pd.concat(
        [
            model.adata.obs.reset_index(names="obs_name"),
            pd.DataFrame(latent, columns=terms),
        ],
        axis=1,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(output, sep="\t", index=False)
    summary = {
        "model": str(args.model),
        "input": str(args.input),
        "output": str(output),
        "posterior_mean": not args.sample_latent,
        "samples": int(len(scores)),
        "terms": int(len(terms)),
    }
    summary_path = output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export expiMap model scores.")
    parser.add_argument("--model", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-latent", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
