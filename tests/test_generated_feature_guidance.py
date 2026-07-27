import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from nasa_mouse_rna_diffusion.generated_feature_guidance import (
    _build_rankings,
    _reactome_enrichment,
    _recenter_draw,
    _selected_indices,
)


class GeneratedFeatureGuidanceTest(unittest.TestCase):
    def setUp(self) -> None:
        self.labels = np.asarray([0, 0, 1, 1, 0, 0, 1, 1])
        self.metadata = pd.DataFrame(
            {
                "accession": ["A"] * 4 + ["B"] * 4,
                "condition": ["ground_control"] * 2
                + ["flight"] * 2
                + ["ground_control"] * 2
                + ["flight"] * 2,
            }
        )
        signal = np.asarray([-1.1, -0.9, 1.0, 1.2, -1.0, -1.2, 0.9, 1.1])
        opposite = -signal
        noise = np.asarray([0.1, -0.1, 0.0, 0.1, -0.1, 0.0, 0.1, -0.1])
        self.real = np.column_stack((signal, opposite, noise))

    def test_recenter_preserves_real_condition_means(self) -> None:
        synthetic = self.real * 2.0 + 3.0
        recentered = _recenter_draw(synthetic, self.real, self.labels)
        for condition in (0, 1):
            mask = self.labels == condition
            np.testing.assert_allclose(
                recentered[mask].mean(axis=0), self.real[mask].mean(axis=0)
            )

    def test_consensus_rankings_retain_matching_signal(self) -> None:
        draws = [self.real + offset for offset in (0.01, -0.02, 0.03)]
        rankings, diagnostics = _build_rankings(
            self.real,
            self.labels,
            self.metadata,
            draws,
            seed=11,
        )
        selected = _selected_indices(rankings["effect_consensus"], 2)
        self.assertEqual(set(selected), {0, 1})
        self.assertTrue(diagnostics.loc[0, "effect_direction_match"])
        self.assertTrue(diagnostics.loc[1, "effect_direction_match"])

    def test_selected_indices_are_deterministic_for_ties(self) -> None:
        selected = _selected_indices(np.asarray([1.0, 1.0, 0.5]), 2)
        np.testing.assert_array_equal(selected, np.asarray([0, 1]))

    def test_reactome_enrichment_reports_symbols(self) -> None:
        with TemporaryDirectory() as directory:
            gmt = Path(directory) / "test.gmt"
            gmt.write_text("term\tdescription\tg1\tg2\tg3\n", encoding="utf-8")
            table = _reactome_enrichment(
                ["g1", "g2"],
                ["g1", "g2", "g3", "g4"],
                gmt,
                {"g1": "Gene1", "g2": "Gene2"},
            )
        self.assertEqual(table.loc[0, "overlap_symbols"], "Gene1,Gene2")


if __name__ == "__main__":
    unittest.main()
