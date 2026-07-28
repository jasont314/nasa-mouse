import unittest

import numpy as np
import pandas as pd

from nasa_mouse_rna_diffusion.within_study_feature_stability import (
    WorkflowData,
    _choose_arms,
    _metric_set,
    _muscle_group_analysis_data,
    _paired_repeat_support,
    _safe_correlation,
    _within_stratum_split,
)


class WithinStudyFeatureStabilityTest(unittest.TestCase):
    def test_muscle_group_mode_relabels_only_selected_skeletal_samples(self) -> None:
        samples = pd.DataFrame(
            {
                "tissue": [
                    "skeletal_muscle",
                    "skeletal_muscle",
                    "skeletal_muscle",
                    "liver",
                ],
                "muscle_group": [
                    "soleus",
                    "soleus",
                    "quadriceps",
                    "not_applicable",
                ],
            }
        )
        data = WorkflowData(
            genes=["gene"],
            symbols={},
            development_expression=np.zeros((4, 1)),
            development_samples=samples.copy(),
            test_expression=np.zeros((4, 1)),
            test_samples=samples.copy(),
            all_expression=np.zeros((4, 1)),
            all_samples=samples.copy(),
            synthetic_draws={"draw": np.zeros((4, 1))},
        )

        grouped, groups = _muscle_group_analysis_data(data, ["soleus"])

        self.assertEqual(groups, ["soleus"])
        self.assertEqual(
            grouped.development_samples["tissue"].tolist(),
            ["soleus", "soleus", "skeletal_muscle", "liver"],
        )
        self.assertEqual(
            data.development_samples["tissue"].tolist(),
            ["skeletal_muscle", "skeletal_muscle", "skeletal_muscle", "liver"],
        )

    def test_muscle_group_mode_rejects_unknown_group(self) -> None:
        samples = pd.DataFrame(
            {"tissue": ["skeletal_muscle"], "muscle_group": ["soleus"]}
        )
        data = WorkflowData(
            genes=["gene"],
            symbols={},
            development_expression=np.zeros((1, 1)),
            development_samples=samples,
            test_expression=np.zeros((1, 1)),
            test_samples=samples,
            all_expression=np.zeros((1, 1)),
            all_samples=samples,
            synthetic_draws={},
        )
        with self.assertRaises(ValueError):
            _muscle_group_analysis_data(data, ["quadriceps"])

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

    def test_safe_correlation_handles_constant_vectors(self) -> None:
        self.assertTrue(
            np.isnan(_safe_correlation(np.ones(4), np.arange(4, dtype=float)))
        )
        self.assertAlmostEqual(
            _safe_correlation(np.arange(4, dtype=float), np.arange(4, dtype=float)),
            1.0,
        )

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
