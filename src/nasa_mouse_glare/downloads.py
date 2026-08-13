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
    if kind == "facs":
        from nasa_mouse_generative.public_references import (
            PUBLIC_REFERENCES,
            download_reference,
        )

        return download_reference(PUBLIC_REFERENCES["tms"], destination=output)
    if output.exists() and output.stat().st_size == item["size_bytes"]:
        print(f"exists\t{output}")
        return output
    if output.exists():
        raise RuntimeError(
            f"Existing file has the wrong size: {output}. Remove it before retrying."
        )
    partial = output.with_name(f"{output.name}.part")
    print(f"downloading\t{item['url']}\n -> {partial}")
    urlretrieve(item["url"], partial)
    if partial.stat().st_size != item["size_bytes"]:
        raise RuntimeError(f"Downloaded file has the wrong size: {partial}")
    partial.replace(output)
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
