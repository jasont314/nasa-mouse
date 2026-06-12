"""Public data links used by the mouse GLARE workflow."""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import urlretrieve


TMS_DATASETS = {
    "facs": {
        "label": "Tabula Muris Senis Smart-seq2/FACS",
        "cells": 110_824,
        "features": 21_025,
        "size_bytes": 2_548_190_251,
        "url": "https://datasets.cellxgene.cziscience.com/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad",
    },
    "droplet": {
        "label": "Tabula Muris Senis 10x/droplet",
        "cells": 245_389,
        "features": 17_943,
        "size_bytes": 3_668_818_742,
        "url": "https://datasets.cellxgene.cziscience.com/084058cc-4f17-43ce-b14e-1278df074013.h5ad",
    },
}


def print_links() -> None:
    for key, item in TMS_DATASETS.items():
        size_gb = item["size_bytes"] / (1024**3)
        print(
            f"{key}\t{item['cells']} cells\t{item['features']} genes\t"
            f"{size_gb:.2f} GiB\t{item['url']}"
        )


def download_dataset(kind: str, output_dir: str | Path) -> Path:
    item = TMS_DATASETS[kind]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / Path(item["url"]).name
    if output.exists() and output.stat().st_size > 0:
        print(f"exists\t{output}")
        return output
    print(f"downloading\t{item['url']}\n -> {output}")
    urlretrieve(item["url"], output)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Download or list TMS assets.")
    parser.add_argument("command", choices=["links", "download"])
    parser.add_argument("--kind", choices=sorted(TMS_DATASETS), default="facs")
    parser.add_argument("--output-dir", default="assets/tms")
    args = parser.parse_args()

    if args.command == "links":
        print_links()
    else:
        download_dataset(args.kind, args.output_dir)


if __name__ == "__main__":
    main()
