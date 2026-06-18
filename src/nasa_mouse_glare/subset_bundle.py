"""Subset a matrix bundle's profiles using an external metadata table."""

from __future__ import annotations

import argparse
from pathlib import Path

from .io import load_matrix_bundle, require_import, write_matrix_bundle


def subset_bundle(
    bundle_manifest: str | Path,
    metadata_path: str | Path,
    output_prefix: str | Path,
    filter_column: str,
    filter_value: str,
    profile_column: str = "profile",
) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    bundle = load_matrix_bundle(bundle_manifest)
    selector = pd.read_csv(metadata_path, sep="\t", keep_default_na=False)
    if profile_column not in selector:
        raise ValueError(
            f"Profile column '{profile_column}' not found in {metadata_path}"
        )
    if filter_column not in selector:
        raise ValueError(
            f"Filter column '{filter_column}' not found in {metadata_path}"
        )
    selector_profiles = selector[profile_column].astype(str).tolist()
    bundle_profiles = [str(profile) for profile in bundle.profiles]
    if selector_profiles == bundle_profiles:
        aligned = selector.reset_index(drop=True)
    else:
        if selector[profile_column].duplicated().any():
            duplicates = selector.loc[
                selector[profile_column].duplicated(keep=False),
                profile_column,
            ].astype(str)
            raise ValueError(
                "Selector metadata is not row-aligned and contains duplicate "
                "profiles: " + ", ".join(duplicates.head(10))
            )
        selector = selector.set_index(profile_column)
        missing = [
            profile for profile in bundle.profiles if profile not in selector.index
        ]
        if missing:
            raise ValueError(
                f"{len(missing)} bundle profiles are missing from selector metadata; "
                f"first missing profile: {missing[0]}"
            )
        aligned = selector.loc[bundle.profiles].reset_index()
    selected = aligned[filter_column].astype(str).eq(filter_value).to_numpy()
    if not selected.any():
        raise ValueError(
            f"No profiles matched {filter_column}={filter_value!r}"
        )

    matrix = bundle.matrix[:, selected]
    profiles = [
        profile
        for profile, keep in zip(bundle.profiles, selected)
        if keep
    ]
    if bundle.profile_metadata is None:
        profile_metadata = aligned.loc[selected].reset_index(drop=True)
    else:
        profile_metadata = bundle.profile_metadata.loc[selected].reset_index(drop=True)
        if "profile" not in profile_metadata:
            profile_metadata.insert(0, "profile", profiles)
        if profile_metadata["profile"].astype(str).tolist() != profiles:
            raise ValueError(
                "Bundle profile metadata is not aligned with the bundle profile order"
            )
        external = aligned.loc[selected].reset_index(drop=True)
        for column in external:
            if column not in profile_metadata:
                profile_metadata[column] = external[column]

    description = (
        f"{bundle_manifest} subset where {filter_column}={filter_value}"
    )
    return write_matrix_bundle(
        output_prefix,
        matrix,
        genes=bundle.genes,
        profiles=profiles,
        profile_metadata=profile_metadata,
        description=description,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Subset matrix-bundle profiles using a metadata field."
    )
    parser.add_argument("--bundle", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--filter-column", required=True)
    parser.add_argument("--filter-value", required=True)
    parser.add_argument("--profile-column", default="profile")
    args = parser.parse_args()
    manifest = subset_bundle(
        args.bundle,
        args.metadata,
        args.output_prefix,
        args.filter_column,
        args.filter_value,
        args.profile_column,
    )
    print(manifest)


if __name__ == "__main__":
    main()
