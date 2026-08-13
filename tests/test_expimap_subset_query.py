from argparse import Namespace
import json

import anndata as ad
import numpy as np
import pandas as pd

from nasa_mouse_expimap.subset_expimap_query import run


def test_subset_expimap_query_excludes_accessions_and_preserves_order(tmp_path):
    input_path = tmp_path / "query.h5ad"
    output_path = tmp_path / "query_subset.h5ad"
    query = ad.AnnData(
        X=np.arange(15, dtype=np.float32).reshape(5, 3),
        obs=pd.DataFrame(
            {"id.accession": ["OSD-1", "OSD-164", "OSD-2", "OSD-168", "OSD-2"]},
            index=[f"sample-{index}" for index in range(5)],
        ),
        var=pd.DataFrame(index=[f"gene-{index}" for index in range(3)]),
    )
    query.layers["counts"] = query.X.copy()
    query.write_h5ad(input_path)

    manifest_path = run(
        Namespace(
            input=str(input_path),
            output=str(output_path),
            exclude_accession=["OSD-164", "OSD-168"],
            accession_column="id.accession",
            expected_samples=3,
            manifest=None,
        )
    )

    subset = ad.read_h5ad(output_path)
    assert subset.obs_names.tolist() == ["sample-0", "sample-2", "sample-4"]
    assert subset.obs["id.accession"].tolist() == ["OSD-1", "OSD-2", "OSD-2"]
    np.testing.assert_array_equal(subset.layers["counts"], query.layers["counts"][[0, 2, 4]])
    assert list(subset.uns["expimap_query_subset"]["excluded_accessions"]) == [
        "OSD-164",
        "OSD-168",
    ]

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["samples_before"] == 5
    assert manifest["samples_after"] == 3
    assert manifest["excluded_sample_counts"] == {"OSD-164": 1, "OSD-168": 1}
