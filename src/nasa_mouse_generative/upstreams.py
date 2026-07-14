"""Install pinned external model source trees used by optional adapters."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

from .adapters.genejepa import PINNED_COMMIT


GENEJEPA_URL = "https://github.com/BiostateAI/GeneJEPA.git"


def prepare_genejepa(path: str | Path) -> Path:
    target = Path(path)
    if not (target / ".git").exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["git", "clone", GENEJEPA_URL, str(target)], check=True
        )
    subprocess.run(["git", "-C", str(target), "fetch", "origin"], check=True)
    subprocess.run(
        ["git", "-C", str(target), "checkout", "--detach", PINNED_COMMIT],
        check=True,
    )
    observed = subprocess.check_output(
        ["git", "-C", str(target), "rev-parse", "HEAD"], text=True
    ).strip()
    if observed != PINNED_COMMIT:
        raise RuntimeError(f"Expected GeneJEPA {PINNED_COMMIT}, observed {observed}")
    return target


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genejepa", default="assets/model_sources/GeneJEPA")
    return parser.parse_args()


def main() -> None:
    print(prepare_genejepa(parse_args().genejepa))


if __name__ == "__main__":
    main()
