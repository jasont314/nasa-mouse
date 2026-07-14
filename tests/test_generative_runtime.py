import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from nasa_mouse_generative.conditioning import CategoryEncoder, UNKNOWN
from nasa_mouse_generative.config import (
    BenchmarkConfig,
    DataConfig,
    ExecutionConfig,
    PreprocessingConfig,
    load_config_with_overrides,
)
from nasa_mouse_generative.preprocessing import FittedPreprocessor
from nasa_mouse_generative.training_data import extract_archs4_matrix


class RuntimeConfigTests(unittest.TestCase):
    def test_dotted_overrides_resolve_model_parameters(self):
        config = load_config_with_overrides(
            "configs/generative/default.yaml",
            [
                "training.regime=osdr_only",
                "training.model_parameters.epochs=3",
                "features.max_genes=64",
            ],
        )
        self.assertEqual(config.training.regime, "osdr_only")
        self.assertEqual(config.training.model_parameters["epochs"], 3)
        self.assertEqual(config.features.max_genes, 64)


class ConditioningTests(unittest.TestCase):
    def test_unseen_categories_use_explicit_unknown_code(self):
        encoder = CategoryEncoder.fit(
            [pd.DataFrame({"condition": ["flight", "ground_control"]})],
            ["condition"],
        )
        observed = encoder.transform(pd.DataFrame({"condition": ["new_condition"]}))
        self.assertEqual(
            int(observed[0, 0]), encoder.vocabularies["condition"].index(UNKNOWN)
        )


class PreprocessorSerializationTests(unittest.TestCase):
    def test_round_trip_and_inverse_transform(self):
        processor = FittedPreprocessor(
            PreprocessingConfig(
                input_units="raw_counts",
                library_normalization="none",
                transform="log1p",
                scaler="zscore",
                harmonization="none",
            )
        )
        matrix = np.asarray([[1.0, 2.0], [3.0, 8.0], [5.0, 4.0]])
        transformed = processor.fit_transform(matrix, ["A", "A", "B"])
        with tempfile.TemporaryDirectory() as directory:
            processor.save(directory)
            restored = FittedPreprocessor.load(directory)
            second = restored.transform(matrix, ["A", "A", "B"])
            inverted = restored.inverse_transform(second, ["A", "A", "B"])
        np.testing.assert_allclose(second, transformed)
        np.testing.assert_allclose(inverted, matrix, atol=1e-5)


class Archs4ExtractionTests(unittest.TestCase):
    def test_selected_columns_and_gene_order_are_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "archs4.h5"
            with h5py.File(source, "w") as handle:
                handle.create_dataset(
                    "meta/genes/ensembl_gene",
                    data=np.asarray([b"g1", b"g2", b"g3"]),
                )
                handle.create_dataset(
                    "data/expression",
                    data=np.asarray([[1, 2], [3, 4], [5, 6]], dtype=np.uint32),
                )
            config = BenchmarkConfig(
                output_root=str(root / "out"),
                data=replace(DataConfig(), archs4_h5=str(source)),
                execution=replace(ExecutionConfig(), cache_archs4=True),
            )
            metadata = pd.DataFrame(
                {
                    "archs4_sample_index": [1, 0],
                    "hierarchical_sampling_weight": [0.5, 0.5],
                }
            )
            matrix, first = extract_archs4_matrix(config, metadata, ["g3", "g1"])
            cached, second = extract_archs4_matrix(config, metadata, ["g3", "g1"])
        np.testing.assert_array_equal(matrix, [[6, 2], [5, 1]])
        np.testing.assert_array_equal(cached, matrix)
        self.assertFalse(first["cache_hit"])
        self.assertTrue(second["cache_hit"])


if __name__ == "__main__":
    unittest.main()
