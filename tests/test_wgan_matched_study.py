from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
import torch

from nasa_mouse_generative.adapters.wgan import WGANAdapter
from nasa_mouse_generative.conditioning import CategoryEncoder
from nasa_mouse_wgan.matched_study import (
    NativeLogZScaler,
    _require_test_unlock,
    initialize_query_embeddings,
    load_config,
)
from nasa_mouse_wgan.training import augment_expression


class MatchedStudyWGANTests(unittest.TestCase):
    def test_paper_augmentation_can_be_disabled_or_applied(self):
        expression = torch.zeros((8, 4), dtype=torch.float32)
        self.assertIs(
            augment_expression(expression, probability=0.0, noise_scale=0.5),
            expression,
        )
        torch.manual_seed(3)
        augmented = augment_expression(
            expression, probability=1.0, noise_scale=0.25
        )
        self.assertFalse(torch.equal(augmented, expression))
        with self.assertRaises(ValueError):
            augment_expression(expression, probability=1.1, noise_scale=0.25)

    def test_native_scaler_round_trip_to_common_scale(self):
        tpm = np.asarray(
            [[0.0, 2.0, 8.0], [1.0, 5.0, 4.0], [4.0, 9.0, 2.0]],
            dtype=np.float32,
        )
        maxabs = np.asarray([4.0, 9.0, 8.0], dtype=np.float32)
        scaler = NativeLogZScaler.fit(tpm)
        observed = scaler.inverse_to_scaled(scaler.transform(tpm), maxabs)
        np.testing.assert_allclose(observed, tpm / maxabs, rtol=1e-5, atol=1e-6)

    def test_query_embeddings_copy_reference_profiles_and_tissue_fallbacks(self):
        covariates = (
            "tissue",
            "condition",
            "study",
            "material_type",
            "sex",
            "muscle_group",
        )
        reference = pd.DataFrame(
            {
                "tissue": ["colon", "retina", "cultured_cells"],
                "condition": ["archs4_reference"] * 3,
                "study": ["archs4_reference"] * 3,
                "material_type": ["archs4_reference"] * 3,
                "sex": ["unknown_sex"] * 3,
                "muscle_group": ["not_applicable"] * 3,
            }
        )
        query = pd.DataFrame(
            {
                "tissue": ["cecum", "eye", "cells"],
                "condition": ["flight", "ground_control", "flight"],
                "study": ["OSD-1", "OSD-2", "OSD-1"],
                "material_type": ["Cecum", "Left eye", "Cells"],
                "sex": ["female", "male", "female"],
                "muscle_group": ["not_applicable"] * 3,
            }
        )
        encoder = CategoryEncoder.fit([reference, query], covariates)
        with tempfile.TemporaryDirectory() as directory:
            adapter = WGANAdapter(
                genes=["g1", "g2"],
                cardinalities=encoder.cardinalities,
                covariates=encoder.covariates,
                parameters={"noise_dim": 2, "hidden_dims": [4], "batch_size": 2},
                device_spec="cpu",
                output_dir=Path(directory),
                checkpoint_every=1,
                resume=False,
                seed=3,
                num_workers=0,
                source_path="",
                validation_partition=None,
            )
            initialize_query_embeddings(adapter, encoder)
            for owner in (
                adapter.model.generator.covariates,
                adapter.model.critic.covariates,
            ):
                tissue_vocab = encoder.vocabularies["tissue"]
                tissue_embedding = owner.embeddings[0].weight.detach()
                for target, source in {
                    "cecum": "colon",
                    "eye": "retina",
                    "cells": "cultured_cells",
                }.items():
                    torch.testing.assert_close(
                        tissue_embedding[tissue_vocab.index(target)],
                        tissue_embedding[tissue_vocab.index(source)],
                    )
                condition_vocab = encoder.vocabularies["condition"]
                condition_embedding = owner.embeddings[1].weight.detach()
                for target in ("flight", "ground_control"):
                    torch.testing.assert_close(
                        condition_embedding[condition_vocab.index(target)],
                        condition_embedding[condition_vocab.index("archs4_reference")],
                    )

    def test_test_evaluation_requires_explicit_unlock(self):
        with self.assertRaises(PermissionError):
            _require_test_unlock(False)
        _require_test_unlock(True)

    def test_checked_in_config_uses_paper_architecture_and_study(self):
        config = load_config(
            "configs/generative/wgan/wgan_matched_study_conditioned.yaml"
        )
        self.assertEqual(config["model"]["hidden_dims"], [256, 256])
        self.assertEqual(config["model"]["critic_steps"], 5)
        self.assertIn("study", config["conditioning"]["covariates"])
        self.assertEqual(config["evaluation"]["minimum_repeat_pass_fraction"], 0.75)
        self.assertEqual(
            config["evaluation"]["calibration_screen_prior_strengths"],
            [1.0, 5.0, 20.0],
        )


if __name__ == "__main__":
    unittest.main()
