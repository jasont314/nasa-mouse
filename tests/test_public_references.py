from dataclasses import replace
import hashlib
from io import BytesIO
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from nasa_mouse_generative.gene_annotations import build_gene_annotations
from nasa_mouse_generative.public_references import (
    PUBLIC_REFERENCES,
    PublicReference,
    download_reference,
    verify_reference,
)
from nasa_mouse_glare import downloads as glare_downloads


class _Response(BytesIO):
    def __init__(self, payload: bytes, status: int = 200, headers=None):
        super().__init__(payload)
        self.status = status
        self.headers = headers or {}

    def getcode(self):
        return self.status

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.close()


def _reference(path: Path, payload: bytes) -> PublicReference:
    return PublicReference(
        key="fixture",
        label="fixture",
        path=path,
        url="https://example.test/fixture.bin",
        source_page="https://example.test",
        size_bytes=len(payload),
        sha256=hashlib.sha256(payload).hexdigest(),
        upstream_sha1=hashlib.sha1(payload).hexdigest(),
    )


class PublicReferenceTests(unittest.TestCase):
    def test_registry_matches_the_tracked_checksum_manifest(self):
        manifest = {}
        for line in Path("assets/EXTERNAL_ARTIFACTS.sha256").read_text().splitlines():
            digest, path = line.split(maxsplit=1)
            manifest[path] = digest

        self.assertEqual(
            {str(item.path): item.sha256 for item in PUBLIC_REFERENCES.values()},
            manifest,
        )
        self.assertEqual(PUBLIC_REFERENCES["archs4"].size_bytes, 38_960_132_574)
        self.assertEqual(PUBLIC_REFERENCES["tms"].size_bytes, 2_548_190_251)
        self.assertEqual(PUBLIC_REFERENCES["gencode"].size_bytes, 91_741_340)

    def test_gencode_gene_annotations_are_versionless(self):
        fixture = (
            "##description: fixture\n"
            'chr1\tHAVANA\tgene\t1\t20\t.\t+\t.\tgene_id '
            '"ENSMUSG00000000001.7"; gene_type "protein_coding"; '
            'gene_name "Gnai3";\n'
            'chr2\tHAVANA\tgene\t30\t50\t.\t-\t.\tgene_id '
            '"ENSMUSG00000000028.15"; gene_type "protein_coding"; '
            'gene_name "Cdc45";\n'
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.gtf"
            path.write_text(fixture, encoding="utf-8")
            observed = build_gene_annotations(path)
        self.assertEqual(
            observed.to_dict(orient="records"),
            [
                {
                    "gene_id": "ENSMUSG00000000001",
                    "gene_symbol": "Gnai3",
                    "gene_type": "protein_coding",
                },
                {
                    "gene_id": "ENSMUSG00000000028",
                    "gene_symbol": "Cdc45",
                    "gene_type": "protein_coding",
                },
            ],
        )

    def test_download_verifies_and_atomically_installs(self):
        payload = b"public-reference-fixture"
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "reference.bin"
            reference = _reference(target, payload)
            with mock.patch(
                "nasa_mouse_generative.public_references.urlopen",
                return_value=_Response(payload),
            ):
                observed = download_reference(reference)
            self.assertEqual(observed.read_bytes(), payload)
            self.assertFalse(target.with_name("reference.bin.part").exists())
            verify_reference(reference)

    def test_download_resumes_a_partial_file(self):
        payload = b"abcdefghijklmnopqrstuvwxyz"
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "reference.bin"
            partial = target.with_name("reference.bin.part")
            partial.write_bytes(payload[:10])
            reference = _reference(target, payload)
            response = _Response(
                payload[10:], status=206, headers={"Content-Range": "bytes 10-25/26"}
            )
            with mock.patch(
                "nasa_mouse_generative.public_references.urlopen",
                return_value=response,
            ) as open_url:
                download_reference(reference)
            request = open_url.call_args.args[0]
            self.assertEqual(request.headers["Range"], "bytes=10-")
            self.assertEqual(target.read_bytes(), payload)

    def test_existing_invalid_complete_file_requires_force(self):
        payload = b"expected"
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "reference.bin"
            target.write_bytes(b"bad-data")
            reference = replace(_reference(target, payload), size_bytes=8)
            with self.assertRaisesRegex(RuntimeError, "Use --force"):
                download_reference(reference)

    @mock.patch("nasa_mouse_generative.public_references.download_reference")
    def test_glare_facs_download_uses_the_verified_reference(self, download):
        with tempfile.TemporaryDirectory() as directory:
            expected = Path(directory) / PUBLIC_REFERENCES["tms"].path.name
            download.return_value = expected
            observed = glare_downloads.download_dataset("facs", directory)
        self.assertEqual(observed, expected)
        download.assert_called_once_with(
            PUBLIC_REFERENCES["tms"], destination=expected
        )


if __name__ == "__main__":
    unittest.main()
