from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from nasa_mouse_generative.condition_figures import (
    _center_within_accession,
    _eligible_accessions,
)


class ConditionFigureTests(unittest.TestCase):
    def test_center_within_accession_removes_group_means(self):
        expression = np.asarray([[1, 3], [3, 5], [10, 20], [14, 24]], dtype=np.float32)
        accessions = pd.Series(["a", "a", "b", "b"])
        centered = _center_within_accession(expression, accessions)
        np.testing.assert_allclose(centered[:2].mean(axis=0), 0.0)
        np.testing.assert_allclose(centered[2:].mean(axis=0), 0.0)

    def test_eligible_accessions_require_both_conditions(self):
        samples = pd.DataFrame(
            {
                "accession": ["a"] * 4 + ["b"] * 3 + ["c"] * 4,
                "condition": (
                    ["flight"] * 2
                    + ["ground_control"] * 2
                    + ["flight"] * 2
                    + ["ground_control"]
                    + ["flight"]
                    + ["ground_control"] * 3
                ),
            }
        )
        self.assertEqual(_eligible_accessions(samples), ["a"])


if __name__ == "__main__":
    unittest.main()
