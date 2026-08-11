import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from nasa_mouse_glare.osdr import (
    load_api_expression,
    normalize_api_metadata,
    select_api_metadata,
    write_api_expression_bundles,
)


def _metadata() -> pd.DataFrame:
    rows = [
        ("OSD-1", "FLT1_techrep1", "flight", "Liver"),
        ("OSD-1", "FLT1_techrep2", "flight", "Liver"),
        ("OSD-1", "GC1", "ground_control", "Liver"),
        ("OSD-2", "FLT2", "flight", "Liver"),
        ("OSD-2", "GC2", "ground_control", "Liver"),
        ("OSD-2", "KID1", "flight", "Kidney"),
    ]
    return pd.DataFrame(
        {
            "id.accession": [row[0] for row in rows],
            "id.sample name": [row[1] for row in rows],
            "condition_inferred": [row[2] for row in rows],
            "tissue_final": [row[3].lower() for row in rows],
            "study.characteristics.material type": [row[3] for row in rows],
            "study.characteristics.sex": ["female"] * len(rows),
            "study.characteristics.strain": ["C57BL/6"] * len(rows),
            "study.characteristics.genotype": ["wild type"] * len(rows),
        }
    )


def _write_counts(path: Path, accession: str, values: dict[str, list[int]]) -> None:
    frame = pd.DataFrame({"gene_id": ["G1", "G2", "G3"], **values})
    frame = frame.rename(
        columns={
            sample: f"{accession}/assay/{sample}"
            for sample in values
        }
    )
    frame.to_csv(path, index=False)


class GlareOsdrApiTests(unittest.TestCase):
    def test_selection_adds_legacy_compatible_aliases(self):
        selected = select_api_metadata(
            normalize_api_metadata(_metadata()),
            tissue="liver",
            accessions=["OSD-1"],
        )
        self.assertEqual(len(selected), 3)
        self.assertEqual(set(selected["condition_label"]), {"FLT", "GC"})
        self.assertEqual(set(selected["h5_accession"]), {"OSD-1"})
        self.assertTrue(selected["material_type"].eq("Liver").all())

    def test_expression_loader_collapses_technical_replicates(self):
        with tempfile.TemporaryDirectory() as tmp:
            counts_dir = Path(tmp) / "counts"
            counts_dir.mkdir()
            _write_counts(
                counts_dir / "OSD-1_unnormalized_counts.csv",
                "OSD-1",
                {
                    "FLT1_techrep1": [1, 2, 3],
                    "FLT1_techrep2": [4, 5, 6],
                    "GC1": [7, 8, 9],
                },
            )
            _write_counts(
                counts_dir / "OSD-2_unnormalized_counts.csv",
                "OSD-2",
                {"FLT2": [2, 4, 6], "GC2": [3, 6, 9], "KID1": [1, 1, 1]},
            )
            selected = select_api_metadata(_metadata(), tissue="liver")
            loaded = load_api_expression(selected, counts_dir=counts_dir)

            raw = loaded["raw_counts"]
            metadata = loaded["metadata"]
            self.assertEqual(raw.shape, (3, 4))
            np.testing.assert_array_equal(
                raw["OSD-1__FLT1"].to_numpy(), np.asarray([5, 7, 9])
            )
            self.assertEqual(
                int(
                    metadata.set_index("feature").loc[
                        "OSD-1__FLT1", "technical_replicate_count"
                    ]
                ),
                2,
            )
            self.assertEqual(loaded["missing_count_columns"].shape[0], 0)
            self.assertTrue(np.isfinite(loaded["log2_cpm"].to_numpy()).all())

    def test_bundle_writer_exports_api_provenance_and_deseq_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            counts_dir = root / "counts"
            counts_dir.mkdir()
            _write_counts(
                counts_dir / "OSD-1_unnormalized_counts.csv",
                "OSD-1",
                {
                    "FLT1_techrep1": [1, 2, 3],
                    "FLT1_techrep2": [4, 5, 6],
                    "GC1": [7, 8, 9],
                },
            )
            selected = select_api_metadata(
                _metadata(), tissue="liver", accessions=["OSD-1"]
            )
            output_dir = root / "bundle"
            summary = write_api_expression_bundles(
                selected, output_dir, counts_dir=counts_dir
            )

            self.assertEqual(summary["input_kind"], "nasa_osdr_api_counts_and_log2_cpm")
            self.assertEqual(summary["samples"], 2)
            self.assertTrue(Path(summary["raw_manifest"]).exists())
            self.assertTrue(Path(summary["log2_cpm_manifest"]).exists())
            self.assertTrue(
                Path(summary["raw_deseq2_inputs"]["metadata"]).exists()
            )


if __name__ == "__main__":
    unittest.main()
