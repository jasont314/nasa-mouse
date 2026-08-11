import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from nasa_mouse_diffusion.paper_parity.whole_study_transfer import (
    _effect_recovery,
    _markdown_table,
    _safe_correlation,
)


class WholeStudyTransferTests(unittest.TestCase):
    def test_markdown_table_has_no_optional_dependency(self):
        frame = pd.DataFrame({"label": ["a|b"], "score": [0.1236], "missing": [np.nan]})
        rendered = _markdown_table(frame)
        self.assertIn("| label | score | missing |", rendered)
        self.assertIn(r"| a\|b | 0.124 | NA |", rendered)

    def test_safe_correlation(self):
        self.assertAlmostEqual(_safe_correlation([1, 2, 3], [2, 4, 6]), 1.0)
        self.assertTrue(np.isnan(_safe_correlation([1, 1], [2, 3])))

    def test_effect_recovery_is_tissue_and_accession_specific(self):
        rows = []
        for tissue in ("liver", "thymus"):
            for accession_index, accession in enumerate(("A", "B", "C")):
                for feature_index, feature in enumerate(("g1", "g2", "g3")):
                    effect = float((feature_index + 1) * (accession_index + 1))
                    rows.append(
                        {
                            "fold": f"fold{accession_index}",
                            "accession": f"{tissue}_{accession}",
                            "tissue": tissue,
                            "feature": feature,
                            "n_flight": 5,
                            "n_ground_control": 5,
                            "flight_minus_ground": effect,
                            "effect_variance": 0.1,
                        }
                    )
        real = pd.DataFrame(rows)
        synthetic = real.copy()
        comparison, accessions, tissues, meta = _effect_recovery(real, synthetic)
        self.assertEqual(len(comparison), 18)
        self.assertEqual(len(accessions), 6)
        self.assertEqual(len(tissues), 2)
        self.assertEqual(len(meta), 6)
        np.testing.assert_allclose(accessions["effect_correlation"], 1.0)
        np.testing.assert_allclose(tissues["meta_effect_correlation"], 1.0)
        np.testing.assert_allclose(tissues["meta_direction_agreement"], 1.0)

    def test_duplicate_outer_test_accession_is_rejected(self):
        real = pd.DataFrame(
            {
                "fold": ["fold0", "fold1"],
                "accession": ["A", "A"],
                "tissue": ["liver", "liver"],
                "feature": ["g1", "g1"],
                "n_flight": [5, 5],
                "n_ground_control": [5, 5],
                "flight_minus_ground": [1.0, 1.0],
                "effect_variance": [0.1, 0.1],
            }
        )
        with self.assertRaisesRegex(ValueError, "more than one outer fold"):
            _effect_recovery(real, real.copy())

    def test_declared_folds_cover_each_eligible_accession_once(self):
        tissues = {
            "adrenal_gland",
            "brain",
            "cerebellum",
            "heart",
            "kidney",
            "liver",
            "lung",
            "retina",
            "skeletal_muscle",
            "skin",
            "spleen",
            "thymus",
        }
        inventory = pd.read_csv(
            "data/osdr_api/osdr_api_mouse_bulk_rnaseq_tissue_accession_counts.tsv",
            sep="\t",
        )
        inventory = inventory.loc[inventory["tissue_final"].isin(tissues)]
        expected = set(zip(inventory["tissue_final"], inventory["id.accession"]))
        observed = set()
        for fold in range(3):
            path = Path(
                f"configs/generative/diffusion/osdr_whole_study_transfer_12t_fold{fold}.yaml"
            )
            config = yaml.safe_load(path.read_text(encoding="utf-8"))
            validation = set(config["data"]["validation_accessions"])
            test = set(config["data"]["test_accessions"])
            self.assertFalse(validation & test)
            for tissue, frame in inventory.groupby("tissue_final"):
                available = set(frame["id.accession"])
                roles = (
                    available - validation - test,
                    available & validation,
                    available & test,
                )
                self.assertTrue(all(roles), (fold, tissue, roles))
                for accession in roles[2]:
                    key = (tissue, accession)
                    self.assertNotIn(key, observed)
                    observed.add(key)
        self.assertEqual(observed, expected)


if __name__ == "__main__":
    unittest.main()
