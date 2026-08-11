import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score

from nasa_mouse_diffusion.paper_parity.classifier_importance import (
    _aggregate_importance,
    _comparison_pattern,
    _linear_shap_rows,
    _metric_matrix_from_logits,
    _permutation_rows,
)


def _classifier_fixture():
    rng = np.random.default_rng(17)
    labels = np.repeat([0, 1], 80)
    expression = rng.normal(size=(160, 3))
    expression[:, 0] += 3.0 * labels
    classifier = LogisticRegression(max_iter=2000).fit(expression, labels)
    return classifier, expression, labels


def test_permutation_importance_identifies_informative_feature():
    classifier, expression, labels = _classifier_fixture()
    table, baseline = _permutation_rows(
        classifier,
        expression,
        labels,
        ["gene_1", "gene_2", "gene_3"],
        {},
        permutation_repeats=20,
        seed=18,
    )

    importance = table.set_index("gene")["permutation_roc_auc_mean"]
    assert baseline["roc_auc"] > 0.95
    assert importance["gene_1"] > 0.25
    assert importance["gene_1"] > importance[["gene_2", "gene_3"]].max()


def test_vectorized_metrics_match_sklearn_with_score_ties():
    labels = np.array([0, 1, 0, 1, 0, 1])
    logits = np.array(
        [
            [-1.0, 0.5, -1.0, 0.5, 0.5, 2.0],
            [0.2, -0.1, 0.2, 1.0, -0.1, 1.0],
        ]
    )
    observed = _metric_matrix_from_logits(labels, logits)
    for row, scores in enumerate(logits):
        assert np.isclose(
            observed["balanced_accuracy"][row],
            balanced_accuracy_score(labels, scores >= 0.0),
        )
        assert np.isclose(observed["roc_auc"][row], roc_auc_score(labels, scores))
        assert np.isclose(
            observed["average_precision"][row],
            average_precision_score(labels, scores),
        )


def test_linear_shap_reconstructs_linear_log_odds():
    classifier, expression, labels = _classifier_fixture()
    background = expression.mean(axis=0)
    table, error = _linear_shap_rows(
        classifier,
        expression,
        labels,
        background,
        ["gene_1", "gene_2", "gene_3"],
        {},
    )

    informative = table.set_index("gene").loc["gene_1"]
    assert error < 1e-10
    assert informative["linear_shap_mean_absolute"] > 1.0
    assert informative["linear_shap_flight_minus_ground"] > 0.0


def test_comparison_pattern_requires_real_transfer_for_emergent_gene():
    common = {
        "real_stable": False,
        "arm_stable": True,
        "coefficient_match": False,
        "real_importance": 0.0,
        "arm_synthetic_importance": 0.2,
        "arm_synthetic_positive_fraction": 1.0,
    }
    assert (
        _comparison_pattern(
            **common,
            arm_real_importance=0.1,
            arm_real_positive_fraction=0.75,
        )
        == "synthetic_emergent_real_transfer"
    )
    assert (
        _comparison_pattern(
            **common,
            arm_real_importance=-0.1,
            arm_real_positive_fraction=0.0,
        )
        == "synthetic_domain_only"
    )


def test_shared_pattern_requires_positive_real_importance():
    assert (
        _comparison_pattern(
            real_stable=True,
            arm_stable=True,
            coefficient_match=True,
            real_importance=-0.02,
            arm_real_importance=-0.01,
            arm_real_positive_fraction=0.25,
            arm_synthetic_importance=0.05,
            arm_synthetic_positive_fraction=0.75,
        )
        == "shared_synthetic_domain_only"
    )
    assert (
        _comparison_pattern(
            real_stable=True,
            arm_stable=True,
            coefficient_match=True,
            real_importance=0.01,
            arm_real_importance=0.02,
            arm_real_positive_fraction=0.75,
            arm_synthetic_importance=0.05,
            arm_synthetic_positive_fraction=0.75,
        )
        == "shared_reinforced"
    )


def test_importance_aggregation_keeps_positive_fraction_keyed():
    rows = []
    for gene, values in {"gene_b": [-0.1, -0.2], "gene_a": [0.1, -0.1]}.items():
        for repeat, value in enumerate(values):
            row = {
                "scope": "tissue",
                "tissue": "thymus",
                "arm": "real_only",
                "domain": "real",
                "gene": gene,
                "symbol": gene,
                "repeat": repeat,
                "classifier_coefficient": 1.0,
                "linear_shap_mean_absolute": 0.1,
                "linear_shap_mean_flight": 0.1,
                "linear_shap_mean_ground_control": -0.1,
                "linear_shap_flight_minus_ground": 0.2,
                "linear_shap_reconstruction_max_error": 0.0,
            }
            for metric in ("balanced_accuracy", "roc_auc", "average_precision"):
                row[f"baseline_{metric}"] = 0.8
                row[f"permutation_{metric}_mean"] = value
                row[f"permutation_{metric}_sd"] = 0.0
                row[f"permutation_{metric}_positive_fraction"] = float(value > 0)
            rows.append(row)
    import pandas as pd

    table = _aggregate_importance(
        pd.DataFrame(rows),
        pd.DataFrame(
            {"scope": ["tissue"], "tissue": ["thymus"], "completed_repeats": [2]}
        ),
    ).set_index("gene")
    assert table.loc["gene_a", "permutation_roc_auc_positive_repeat_fraction"] == 0.5
    assert table.loc["gene_b", "permutation_roc_auc_positive_repeat_fraction"] == 0.0
