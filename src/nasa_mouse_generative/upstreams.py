"""Install pinned external model source trees used by optional adapters."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

from .paper_contracts import PAPER_SOURCES, verify_pinned_source


GENEJEPA_URL = PAPER_SOURCES["genejepa"]["url"]


def prepare_source(model: str, path: str | Path) -> Path:
    contract = PAPER_SOURCES[model]
    target = Path(path)
    if not (target / ".git").exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "clone", contract["url"], str(target)], check=True)
    subprocess.run(["git", "-C", str(target), "fetch", "origin"], check=True)
    subprocess.run(
        ["git", "-C", str(target), "checkout", "--detach", contract["commit"]],
        check=True,
    )
    verify_pinned_source(model, target)
    return target


def prepare_genejepa(path: str | Path) -> Path:
    return prepare_source("genejepa", path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genejepa", default="assets/model_sources/GeneJEPA")
    parser.add_argument(
        "--wgan", default="assets/model_sources/adversarial-gene-expression"
    )
    parser.add_argument(
        "--diffusion", default="assets/model_sources/rna-diffusion"
    )
    parser.add_argument("--mbatch", default="assets/model_sources/MBatch")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = {
        "vinas_wgan_gp": prepare_source("vinas_wgan_gp", args.wgan),
        "lacan_diffusion": prepare_source("lacan_diffusion", args.diffusion),
        "genejepa": prepare_source("genejepa", args.genejepa),
        "mbatch": prepare_source("mbatch", args.mbatch),
    }
    for model, path in paths.items():
        print(f"{model}\t{path}")


if __name__ == "__main__":
    main()
