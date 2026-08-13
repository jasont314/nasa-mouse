"""Download and verify the public reference files used by the workflows."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import sys
from urllib.request import Request, urlopen


@dataclass(frozen=True)
class PublicReference:
    """Immutable identity and local destination for a public data file."""

    key: str
    label: str
    path: Path
    url: str
    source_page: str
    size_bytes: int
    sha256: str
    upstream_sha1: str | None = None


PUBLIC_REFERENCES = {
    "archs4": PublicReference(
        key="archs4",
        label="ARCHS4 mouse gene expression v2.5",
        path=Path("assets/archs4/mouse_gene_v2.5.h5"),
        url=(
            "https://s3.dev.maayanlab.cloud/archs4/files/"
            "mouse_gene_v2.5.h5"
        ),
        source_page="https://maayanlab.cloud/archs4/download.html",
        size_bytes=38_960_132_574,
        sha256="74b509f82623bced395119244becf30df601a24fcaaf905691e2716bf83118b8",
        upstream_sha1="22605c9b6c4e7502b0861d4d8591ce128907c39f",
    ),
    "tms": PublicReference(
        key="tms",
        label="Tabula Muris Senis Smart-seq2/FACS",
        path=Path("assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad"),
        url=(
            "https://datasets.cellxgene.cziscience.com/"
            "be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad"
        ),
        source_page=(
            "https://cellxgene.cziscience.com/collections/"
            "0b9d8a04-bb9d-44da-aa27-705bb65b54eb"
        ),
        size_bytes=2_548_190_251,
        sha256="1d7fd90acb33269c3337dc5031b4a89d9aa4f72806a45b9c12e768fedc8acf8f",
    ),
    "gencode": PublicReference(
        key="gencode",
        label="GENCODE mouse primary-assembly annotation vM39",
        path=Path(
            "assets/reference/gencode.vM39.primary_assembly.annotation.gtf.gz"
        ),
        url=(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
            "release_M39/gencode.vM39.primary_assembly.annotation.gtf.gz"
        ),
        source_page="https://www.gencodegenes.org/mouse/release_M39.html",
        size_bytes=91_741_340,
        sha256="d6da97913ce30f99883fc1216b111569f9947cf203886a5afb607b59228574d4",
    ),
}


def _digest(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_reference(reference: PublicReference, path: Path | None = None) -> None:
    """Raise when a local reference differs from its recorded identity."""

    target = path or reference.path
    if not target.is_file():
        raise FileNotFoundError(f"Missing {reference.label}: {target}")
    observed_size = target.stat().st_size
    if observed_size != reference.size_bytes:
        raise RuntimeError(
            f"Wrong size for {target}: expected {reference.size_bytes}, "
            f"observed {observed_size}"
        )
    observed_sha256 = _digest(target)
    if observed_sha256 != reference.sha256:
        raise RuntimeError(
            f"SHA-256 mismatch for {target}: expected {reference.sha256}, "
            f"observed {observed_sha256}"
        )
    if reference.upstream_sha1 is not None:
        observed_sha1 = _digest(target, "sha1")
        if observed_sha1 != reference.upstream_sha1:
            raise RuntimeError(
                f"SHA-1 mismatch for {target}: expected "
                f"{reference.upstream_sha1}, observed {observed_sha1}"
            )


def _prepare_partial(
    reference: PublicReference, target: Path, partial: Path, force: bool
) -> int:
    if force:
        target.unlink(missing_ok=True)
        partial.unlink(missing_ok=True)
        return 0

    if target.exists():
        try:
            verify_reference(reference, target)
        except RuntimeError:
            target_size = target.stat().st_size
            if target_size < reference.size_bytes and not partial.exists():
                target.replace(partial)
            else:
                raise RuntimeError(
                    f"Existing file failed verification: {target}. "
                    "Use --force to replace it."
                )
        else:
            return reference.size_bytes

    if partial.exists() and partial.stat().st_size > reference.size_bytes:
        raise RuntimeError(
            f"Partial file is larger than expected: {partial}. "
            "Use --force to restart it."
        )
    return partial.stat().st_size if partial.exists() else 0


def download_reference(
    reference: PublicReference,
    *,
    destination: Path | None = None,
    force: bool = False,
    timeout: int = 120,
) -> Path:
    """Resume, verify, and atomically install one public reference."""

    target = destination or reference.path
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(f"{target.name}.part")
    offset = _prepare_partial(reference, target, partial, force)
    if offset == reference.size_bytes:
        if target.exists():
            print(f"verified\t{target}")
            return target
        verify_reference(reference, partial)
        os.replace(partial, target)
        print(f"verified\t{target}")
        return target

    headers = {"User-Agent": "nasa-mouse-public-reference-downloader/1.0"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    request = Request(reference.url, headers=headers)
    print(
        f"downloading\t{reference.label}\n"
        f"source\t{reference.url}\n"
        f"target\t{target}\n"
        f"resume_bytes\t{offset}",
        flush=True,
    )
    with urlopen(request, timeout=timeout) as response:
        status = getattr(response, "status", response.getcode())
        append = bool(offset and status == 206)
        if append:
            content_range = response.headers.get("Content-Range", "")
            if not content_range.startswith(f"bytes {offset}-"):
                raise RuntimeError(
                    f"Server returned an unexpected byte range: {content_range!r}"
                )
        mode = "ab" if append else "wb"
        downloaded = offset if append else 0
        next_report = downloaded + 1024**3
        with partial.open(mode) as handle:
            while True:
                block = response.read(8 * 1024 * 1024)
                if not block:
                    break
                handle.write(block)
                downloaded += len(block)
                if downloaded >= next_report:
                    print(
                        f"progress_gib\t{downloaded / 1024**3:.1f}",
                        file=sys.stderr,
                        flush=True,
                    )
                    next_report = downloaded + 1024**3

    verify_reference(reference, partial)
    os.replace(partial, target)
    print(f"verified\t{target}")
    return target


def _selected(keys: list[str] | None) -> list[PublicReference]:
    requested = keys or list(PUBLIC_REFERENCES)
    return [PUBLIC_REFERENCES[key] for key in requested]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference",
        action="append",
        choices=sorted(PUBLIC_REFERENCES),
        help="Reference to process; repeat to select several (default: all)",
    )
    parser.add_argument(
        "--list", action="store_true", help="Print source and destination metadata"
    )
    parser.add_argument(
        "--check", action="store_true", help="Verify local files without downloading"
    )
    parser.add_argument(
        "--force", action="store_true", help="Discard an invalid or partial download"
    )
    parser.add_argument("--timeout", type=int, default=120)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    references = _selected(args.reference)
    if args.list:
        for reference in references:
            print(
                "\t".join(
                    [
                        reference.key,
                        str(reference.path),
                        str(reference.size_bytes),
                        reference.sha256,
                        reference.url,
                        reference.source_page,
                    ]
                )
            )
        return

    for reference in references:
        if args.check:
            verify_reference(reference)
            print(f"verified\t{reference.path}")
        else:
            download_reference(reference, force=args.force, timeout=args.timeout)


if __name__ == "__main__":
    main()
