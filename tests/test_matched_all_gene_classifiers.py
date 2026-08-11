import numpy as np
import pandas as pd

from nasa_mouse_diffusion.paper_parity.matched_all_gene_classifiers import (
    _accession_blocks,
    _bh_fdr_crosswalk,
    _candidate_reactome_enrichment,
    _fit_matched_arm,
    _importance_pattern,
    _metric_summary,
    _positive_importance,
    _select_shared_regularization,
    _synthetic_accession_blocks,
)


def _toy_expression():
    rng = np.random.default_rng(31)
    labels = np.repeat([0, 1], 20)
    expression = rng.normal(size=(40, 6))
    expression[:, 0] += 3.0 * labels
    samples = pd.DataFrame(
        {
            "accession": np.tile(np.repeat(["A", "B"], 10), 2),
            "condition": np.where(labels == 1, "flight", "ground_control"),
        }
    )
    return expression, labels, samples


def test_shared_regularization_uses_real_validation_and_returns_full_grid():
    expression, labels, samples = _toy_expression()
    selected, candidates = _select_shared_regularization(
        expression,
        labels,
        samples,
        inner_train=np.r_[0:8, 20:28],
        inner_validation=np.r_[8:20, 28:40],
        regularization_grid=[0.001, 0.01, 0.1],
        seed=32,
    )
    assert selected in {0.001, 0.01, 0.1}
    assert candidates["regularization_c"].tolist() == [0.001, 0.01, 0.1]
    assert candidates["roc_auc"].max() > 0.8


def test_all_matched_arms_retain_every_gene():
    expression, labels, _ = _toy_expression()
    synthetic = [expression + 0.01, expression - 0.01, expression + 0.02]
    classifiers = [
        _fit_matched_arm(
            arm,
            real_scaled=expression,
            labels=labels,
            synthetic_scaled=synthetic,
            regularization_c=0.01,
            synthetic_weight=1.0,
            seed=40 + index,
        )
        for index, arm in enumerate(
            ("real_only", "generated_only", "real_plus_generated")
        )
    ]
    assert all(classifier.coef_.shape == (1, expression.shape[1]) for classifier in classifiers)


def test_accession_blocks_partition_real_and_each_synthetic_draw():
    samples = pd.DataFrame({"accession": ["A", "A", "B", "B", "B"]})
    real_blocks = _accession_blocks(samples)
    synthetic_blocks = _synthetic_accession_blocks(samples, 3)
    assert sorted(np.concatenate(real_blocks).tolist()) == list(range(5))
    assert len(synthetic_blocks) == 6
    assert sorted(np.concatenate(synthetic_blocks).tolist()) == list(range(15))


def test_importance_pattern_separates_transfer_from_synthetic_domain_only():
    frame = pd.DataFrame(
        {
            "real_only_positive_importance": [False, False, True, True],
            "arm_real_positive_importance": [True, False, True, False],
            "arm_synthetic_positive_importance": [True, True, True, False],
            "coefficient_direction_match": [True, True, True, True],
            "real_only_permutation_roc_auc_mean": [0.0, 0.0, 0.02, 0.03],
            "arm_real_permutation_roc_auc_mean": [0.03, 0.0, 0.04, 0.0],
        }
    )
    assert _importance_pattern(frame).tolist() == [
        "synthetic_promoted_real_transfer",
        "synthetic_domain_only",
        "shared_reinforced",
        "real_only",
    ]


def test_positive_importance_requires_practical_loss_and_repeat_consistency():
    frame = pd.DataFrame(
        {
            "arm_permutation_roc_auc_mean": [0.0009, 0.0010, 0.02, 0.02],
            "arm_permutation_roc_auc_positive_repeat_fraction": [1.0, 0.5, 0.49, 0.5],
        }
    )
    assert _positive_importance(
        frame,
        "arm",
        minimum_importance=0.001,
        minimum_fraction=0.5,
    ).tolist() == [False, True, False, True]


def test_bh_crosswalk_joins_by_stable_gene_id_when_symbols_differ(tmp_path):
    inventory_path = tmp_path / "bh.tsv"
    pd.DataFrame(
        {
            "analysis_scope": ["canonical_tissue"],
            "tissue": ["thymus"],
            "gene": ["ENSMUSG1"],
            "symbol": ["OldName"],
            "meta_effect": [0.2],
        }
    ).to_csv(inventory_path, sep="\t", index=False)
    comparison = pd.DataFrame(
        {
            "scope": ["tissue"],
            "tissue": ["thymus"],
            "arm": ["real_plus_generated"],
            "gene": ["ENSMUSG1"],
            "symbol": ["NewName"],
            "arm_real_median_classifier_coefficient": [0.4],
            "arm_real_linear_shap_flight_minus_ground": [0.1],
            "pattern": ["synthetic_promoted_real_transfer"],
        }
    )
    utility = pd.DataFrame(
        {
            "scope": ["tissue"],
            "tissue": ["thymus"],
            "arm": ["real_plus_generated"],
            "pooled_mean_all_metrics_nonworse": [True],
            "macro_mean_all_metrics_nonworse": [True],
            "joint_mean_all_metrics_nonworse": [True],
        }
    )
    result = _bh_fdr_crosswalk(inventory_path, comparison, utility)
    assert result.loc[0, "symbol"] == "OldName"
    assert result.loc[0, "importance_symbol"] == "NewName"
    assert not bool(result.loc[0, "symbol_matches_importance_annotation"])
    assert bool(result.loc[0, "eligible_synthetic_biological_candidate"])


def test_candidate_reactome_enrichment_reports_all_and_direction_sets(tmp_path):
    gmt_path = tmp_path / "reactome.gmt"
    gmt_path.write_text(
        "TERM\tdescription\tG1\tG2\tG3\n",
        encoding="utf-8",
    )
    eligible = pd.DataFrame(
        {
            "scope": ["tissue", "tissue"],
            "tissue": ["thymus", "thymus"],
            "gene": ["G1", "G2"],
            "flt_gc_direction": ["FLT_lower", "FLT_lower"],
        }
    )
    result = _candidate_reactome_enrichment(
        eligible,
        background=["G1", "G2", "G3", "G4", "G5"],
        gmt_path=gmt_path,
        symbols={"G1": "One", "G2": "Two"},
    )
    assert set(result["gene_set"]) == {"all", "flt_lower"}
    assert set(result["overlap_symbols"]) == {"One,Two"}


def test_metric_summary_requires_pooled_and_macro_nonworse():
    rows = []
    for repeat in range(2):
        for arm, shift, macro_shift in (
            ("real_only", 0.0, 0.0),
            ("generated_only", 0.05, -0.01),
            ("real_plus_generated", 0.02, 0.01),
        ):
            row = {
                "scope": "tissue",
                "tissue": "thymus",
                "repeat": repeat,
                "arm": arm,
            }
            for metric in ("balanced_accuracy", "roc_auc", "average_precision"):
                row[metric] = 0.7 + shift
                row[f"accession_macro_{metric}"] = 0.7 + macro_shift
            rows.append(row)
    _, utility = _metric_summary(pd.DataFrame(rows))
    lookup = utility.set_index("arm")
    assert not bool(lookup.loc["generated_only", "joint_mean_all_metrics_nonworse"])
    assert bool(
        lookup.loc["real_plus_generated", "joint_mean_all_metrics_nonworse"]
    )
