import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from nasa_mouse_diffusion.paper_parity.classifier_importance import _linear_shap_rows
from nasa_mouse_diffusion.paper_parity.grouped_pathway_importance import (
    PathwayGroup,
    _group_linear_shap_rows,
    _group_permutation_rows,
    _load_pathway_groups,
    _membership_matrix,
    _within_accession_pathway_scores,
)


def _correlated_example():
    rng = np.random.default_rng(44)
    labels = np.repeat([0, 1], 60)
    signal = labels + rng.normal(scale=0.25, size=len(labels))
    expression = np.column_stack(
        (
            signal + rng.normal(scale=0.02, size=len(labels)),
            signal + rng.normal(scale=0.02, size=len(labels)),
            rng.normal(size=len(labels)),
        )
    )
    classifier = LogisticRegression(C=0.1, max_iter=5000).fit(expression, labels)
    groups = [
        PathwayGroup(
            term="PAIR",
            description="correlated pair",
            url="https://example.test/pair",
            indices=(0, 1),
            genes=("G1", "G2"),
            symbols=("One", "Two"),
        )
    ]
    return expression, labels, classifier, groups


def test_group_permutation_jointly_removes_correlated_signal():
    expression, labels, classifier, groups = _correlated_example()
    membership = _membership_matrix(expression.shape[1], groups)
    grouped, baseline = _group_permutation_rows(
        classifier,
        expression,
        labels,
        groups,
        membership,
        permutation_repeats=20,
        seed=45,
        blocks=[np.arange(len(labels))],
    )
    assert baseline["roc_auc"] > 0.95
    assert grouped.loc[0, "permutation_roc_auc_mean"] > 0.2
    assert grouped.loc[0, "permutation_roc_auc_positive_fraction"] > 0.9


def test_grouped_linear_shap_equals_sum_of_member_means():
    expression, labels, classifier, groups = _correlated_example()
    membership = _membership_matrix(expression.shape[1], groups)
    background = expression.mean(axis=0)
    grouped = _group_linear_shap_rows(
        classifier,
        expression,
        labels,
        background,
        groups,
        membership,
    )
    genes, _ = _linear_shap_rows(
        classifier,
        expression,
        labels,
        background,
        ["G1", "G2", "G3"],
        {"G1": "One", "G2": "Two", "G3": "Three"},
    )
    expected = genes.loc[
        genes["gene"].isin(["G1", "G2"]), "linear_shap_flight_minus_ground"
    ].sum()
    assert np.isclose(grouped.loc[0, "group_shap_flight_minus_ground"], expected)


def test_pathway_loader_intersects_and_orders_classifier_genes(tmp_path):
    path = tmp_path / "mouse.gmt"
    path.write_text(
        "TERM\tdescription\tG3\tMISSING\tG1\tG2\n",
        encoding="utf-8",
    )
    groups = _load_pathway_groups(
        path,
        ["G1", "G2", "G3", "G4"],
        {"G1": "One", "G2": "Two", "G3": "Three"},
        minimum_genes=3,
        maximum_genes=10,
    )
    assert groups[0].genes == ("G1", "G2", "G3")
    assert groups[0].indices == (0, 1, 2)


def test_pathway_scores_are_centered_within_accession():
    expression = np.asarray(
        [
            [1.0, 2.0],
            [3.0, 4.0],
            [10.0, 20.0],
            [14.0, 24.0],
        ]
    )
    samples = pd.DataFrame({"accession": ["A", "A", "B", "B"]})
    membership = np.ones((2, 1), dtype=float)
    scores = _within_accession_pathway_scores(expression, samples, membership)
    assert np.allclose(scores[:2].mean(axis=0), 0.0)
    assert np.allclose(scores[2:].mean(axis=0), 0.0)
