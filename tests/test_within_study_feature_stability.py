import unittest

import numpy as np
import pandas as pd

from nasa_mouse_rna_diffusion.within_study_feature_stability import (
    _choose_arms,
    _metric_set,
    _paired_repeat_support,
    _within_stratum_split,
)


class WithinStudyFeatureStabilityTest(unittest.TestCase):
    def test_split_retains_each_accession_condition_stratum(self) -> None:
        samples = pd.DataFrame(
            {
                "accession": ["A"] * 8 + ["B"] * 8,
                "condition": (["flight"] * 4 + ["ground_control"] * 4) * 2,
            }
        )
        retained, held_out = _within_stratum_split(
            samples, fraction=0.25, seed=7
        )
        self.assertEqual(len(held_out), 4)
        retained_samples = samples.loc[retained]
        self.assertTrue(
            (retained_samples.groupby(["accession", "condition"]).size() >= 1).all()
        )

    def test_metric_set_reports_all_three_endpoints(self) -> None:
        observed = _metric_set(
            np.asarray([0, 0, 1, 1]),
            np.asarray([0.1, 0.4, 0.6, 0.9]),
        )
        self.assertEqual(set(observed), {
            "balanced_accuracy",
            "roc_auc",
            "average_precision",
        })
        self.assertEqual(observed["balanced_accuracy"], 1.0)

    def test_arm_choice_requires_every_metric_to_tie_or_improve(self) -> None:
        table = pd.DataFrame(
            [
                {
                    "tissue": "liver",
                    "arm": "real_only",
                    "mean_balanced_accuracy": 0.7,
                    "mean_roc_auc": 0.8,
                    "mean_average_precision": 0.75,
                },
                {
                    "tissue": "liver",
                    "arm": "real_plus_generated",
                    "mean_balanced_accuracy": 0.8,
                    "mean_roc_auc": 0.79,
                    "mean_average_precision": 0.9,
                },
                {
                    "tissue": "liver",
                    "arm": "guided_real_only",
                    "mean_balanced_accuracy": 0.72,
                    "mean_roc_auc": 0.81,
                    "mean_average_precision": 0.8,
                },
            ]
        )
        selected = _choose_arms(table, pd.DataFrame()).iloc[0]
        self.assertEqual(selected["selected_arm"], "guided_real_only")

    def test_paired_repeat_support_uses_matched_splits(self) -> None:
        metrics = pd.DataFrame(
            [
                {
                    "tissue": "liver",
                    "repeat": 0,
                    "arm": "real_only",
                    "balanced_accuracy": 0.6,
                    "roc_auc": 0.7,
                    "average_precision": 0.8,
                },
                {
                    "tissue": "liver",
                    "repeat": 1,
                    "arm": "real_only",
                    "balanced_accuracy": 0.7,
                    "roc_auc": 0.8,
                    "average_precision": 0.8,
                },
                {
                    "tissue": "liver",
                    "repeat": 0,
                    "arm": "guided_real_only",
                    "balanced_accuracy": 0.7,
                    "roc_auc": 0.8,
                    "average_precision": 0.8,
                },
                {
                    "tissue": "liver",
                    "repeat": 1,
                    "arm": "guided_real_only",
                    "balanced_accuracy": 0.6,
                    "roc_auc": 0.9,
                    "average_precision": 0.9,
                },
            ]
        )
        choices = pd.DataFrame(
            [{"tissue": "liver", "selected_arm": "guided_real_only"}]
        )
        support = _paired_repeat_support(metrics, choices).iloc[0]
        self.assertEqual(support["paired_repeats"], 2)
        self.assertAlmostEqual(support["mean_delta_balanced_accuracy"], 0.0)
        self.assertAlmostEqual(support["nonworse_rate_roc_auc"], 1.0)
        self.assertAlmostEqual(support["all_metrics_nonworse_rate"], 0.5)


if __name__ == "__main__":
    unittest.main()
