import unittest

import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score

from nasa_mouse_glare.paper_clustering import compressed_eac
from nasa_mouse_glare.paper_finetune import profile_cohort


class PaperWorkflowTests(unittest.TestCase):
    def test_profile_cohort(self):
        self.assertEqual(
            profile_cohort("RR8_LVR_FLT_ISS-T_OLD_FI11"), ("ISS-T", "OLD", 11)
        )
        self.assertEqual(
            profile_cohort("RR8_LVR_GC_LAR_YNG_GL6"), ("LAR", "YNG", 6)
        )

    def test_compressed_eac_matches_full_average_linkage(self):
        base_labels = np.asarray(
            [
                [0, 0, 0],
                [0, 0, 0],
                [0, 1, 0],
                [1, 1, 1],
                [1, 1, 1],
                [1, 0, 1],
            ]
        )
        coassociation = np.mean(
            base_labels[:, None, :] == base_labels[None, :, :], axis=2
        )
        expected = AgglomerativeClustering(
            n_clusters=2,
            metric="precomputed",
            linkage="average",
        ).fit_predict(1.0 - coassociation)
        observed = compressed_eac(base_labels, n_clusters=2)
        self.assertEqual(adjusted_rand_score(expected, observed), 1.0)


if __name__ == "__main__":
    unittest.main()
