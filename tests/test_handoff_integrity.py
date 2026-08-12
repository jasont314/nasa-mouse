from pathlib import Path
import sys
import unittest
from unittest import mock

import pandas as pd

from nasa_mouse_diffusion.paper_parity import build_synthetic_guided_paper


ROOT = Path(__file__).resolve().parents[1]


class HandoffIntegrityTests(unittest.TestCase):
    def test_internship_report_manifest_inputs_exist(self):
        manifest = pd.read_csv(
            ROOT / "paper/slstp_internship_report/source_data/source_manifest.tsv",
            sep="\t",
        )
        missing = [path for path in manifest["path"] if not (ROOT / path).exists()]
        self.assertEqual(missing, [])

    def test_generative_render_only_skips_analysis_refresh(self):
        arguments = ["build_synthetic_guided_paper", "--render-only"]
        with (
            mock.patch.object(sys, "argv", arguments),
            mock.patch.object(
                build_synthetic_guided_paper, "render_paper_documents"
            ) as render,
            mock.patch.object(
                build_synthetic_guided_paper, "build_source_tables"
            ) as refresh,
        ):
            build_synthetic_guided_paper.main()

        render.assert_called_once_with()
        refresh.assert_not_called()


if __name__ == "__main__":
    unittest.main()
