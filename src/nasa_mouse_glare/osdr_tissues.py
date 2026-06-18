"""Extract authoritative tissue metadata from cached OSDR ISA study records."""

from __future__ import annotations

import argparse
from collections import Counter
import html
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen

from .cluster_stratified_analysis import infer_tissue
from .io import require_import


DEFAULT_METADATA_DIR = "assets/osdr_metadata"
DEFAULT_PROFILE_METADATA = (
    "outputs/glare_hpt_tms_facs_osdr/post_finetune/profile_metadata.tsv"
)
DEFAULT_OUTPUT_DIR = (
    "outputs/glare_hpt_tms_facs_osdr/post_finetune/osdr_tissues"
)


def download_json(url: str, output: Path, timeout: int) -> None:
    request = Request(url, headers={"User-Agent": "nasa-mouse/1.0"})
    with urlopen(request, timeout=timeout) as response:
        payload = response.read()
    json.loads(payload)
    output.write_bytes(payload)


def download_metadata(
    profile_metadata: str | Path,
    metadata_dir: Path,
    timeout: int,
    refresh: bool,
) -> None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    profiles = pd.read_csv(profile_metadata, sep="\t", keep_default_na=False)
    accessions = sorted(
        {
            accession
            for _, row in profiles.iterrows()
            if (accession := resolved_accession(row))
        }
    )
    metadata_dir.mkdir(parents=True, exist_ok=True)
    for index, accession in enumerate(accessions, start=1):
        study_path = metadata_dir / f"{accession}.json"
        samples_path = metadata_dir / f"{accession}.samples.json"
        study_url = (
            f"https://osdr.nasa.gov/geode-py/ws/repo/studies/{accession}"
        )
        samples_url = (
            f"{study_url}/table?page=1&page_size=1000&table_name=sample"
        )
        print(f"download {index}/{len(accessions)} {accession}", flush=True)
        if refresh or not study_path.exists():
            download_json(study_url, study_path, timeout)
        if refresh or not samples_path.exists():
            download_json(samples_url, samples_path, timeout)


def clean_value(value) -> str:
    value = "" if value is None else str(value)
    value = re.sub(r"<[^>]+>", "", value)
    return " ".join(html.unescape(value).split())


def normalized_title(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def normalized_identifier(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", clean_value(value).lower())


def resolved_accession(row) -> str:
    accession = clean_value(row.get("id.accession", ""))
    if accession:
        return accession
    accession_sample = clean_value(row.get("id.accession_sample name", ""))
    match = re.match(r"^(OSD-\d+)(?:_|$)", accession_sample, flags=re.IGNORECASE)
    return match.group(1).upper() if match else ""


def field_for_title(header: list[dict], *titles: str) -> str:
    wanted = {normalized_title(title) for title in titles}
    for column in header:
        if normalized_title(column.get("title", "")) in wanted:
            return str(column.get("field", ""))
    return ""


def canonical_tissue(material_type: str) -> str:
    value = clean_value(material_type).lower()
    if not value:
        return "unknown"
    if "liver" in value:
        return "liver"
    if "kidney" in value:
        return "kidney"
    if "spleen" in value:
        return "spleen"
    if "thymus" in value:
        return "thymus"
    if "lung" in value:
        return "lung"
    if "retina" in value:
        return "retina"
    if "optic nerve" in value:
        return "optic_nerve"
    if value == "eye" or value.endswith(" eye"):
        return "eye"
    if "cerebell" in value:
        return "cerebellum"
    if "hippocamp" in value:
        return "hippocampus"
    if "cerebral hemisphere" in value or value in {"brain", "cerebrum"}:
        return "brain"
    if "ventricle" in value or value == "heart":
        return "heart"
    if "adrenal" in value:
        return "adrenal_gland"
    if "colon" in value:
        return "colon"
    if value == "cecum":
        return "cecum"
    if "skin" in value:
        return "skin"
    if "mammary tumor" in value:
        return "mammary_tumor"
    if any(
        token in value
        for token in [
            "soleus",
            "gastrocnemius",
            "quadriceps",
            "tibialis anterior",
            "extensor digitorum longus",
        ]
    ):
        return "skeletal_muscle"
    if "bone marrow" in value:
        return "bone_marrow"
    if "bone" in value or value == "mandible":
        return "bone"
    if value == "blood":
        return "blood"
    if value == "placenta":
        return "placenta"
    if value == "zygote":
        return "zygote"
    if "cultured" in value:
        return "cultured_cells"
    if value == "cells":
        return "cells"
    if value in {"not applicable", "tissue"}:
        return "unspecified"
    return re.sub(r"[^a-z0-9]+", "_", value).strip("_") or "unknown"


def unique_mapping(rows: list[dict], key: str) -> dict[str, dict]:
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        value = clean_value(row.get(key, ""))
        if value:
            grouped.setdefault(value, []).append(row)
    return {
        value: matches[0]
        for value, matches in grouped.items()
        if len(matches) == 1
    }


def unique_normalized_mapping(rows: list[dict], key: str) -> dict[str, dict]:
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        value = normalized_identifier(row.get(key, ""))
        if value:
            grouped.setdefault(value, []).append(row)
    return {
        value: matches[0]
        for value, matches in grouped.items()
        if len(matches) == 1
    }


def load_studies(metadata_dir: Path):
    studies = {}
    sample_rows = {}
    for study_path in sorted(metadata_dir.glob("OSD-*.json")):
        if study_path.name.endswith(".samples.json"):
            continue
        accession = study_path.stem
        samples_path = metadata_dir / f"{accession}.samples.json"
        if not samples_path.exists():
            continue
        study = json.loads(study_path.read_text(encoding="utf-8"))
        page = json.loads(samples_path.read_text(encoding="utf-8"))
        header = study.get("samples", {}).get("header", [])
        sample_field = field_for_title(header, "Sample Name")
        source_field = field_for_title(header, "Source Name")
        material_field = field_for_title(
            header,
            "Characteristics: Material Type",
            "Characteristics[Material Type]",
        )
        tissue_field = field_for_title(
            header,
            "Characteristics: Tissue Type",
            "Characteristics[Tissue Type]",
            "Factor Value: Tissue",
            "Factor Value[Tissue]",
        )
        rows = []
        for raw_row in page.get("tableData", {}).get("current", []):
            material = clean_value(raw_row.get(material_field, "")) if material_field else ""
            tissue_type = (
                clean_value(raw_row.get(tissue_field, "")) if tissue_field else ""
            )
            material_tissue = canonical_tissue(material)
            secondary_tissue = canonical_tissue(tissue_type)
            if material_tissue in {
                "unknown",
                "unspecified",
                "cells",
                "cultured_cells",
            } and secondary_tissue not in {"unknown", "unspecified"}:
                official_tissue = secondary_tissue
                official_tissue_source = "official_osdr_tissue_type"
            else:
                official_tissue = material_tissue
                official_tissue_source = "official_osdr_material_type"
            rows.append(
                {
                    "accession": accession,
                    "official_sample_name": clean_value(raw_row.get(sample_field, "")),
                    "official_source_name": clean_value(raw_row.get(source_field, "")),
                    "official_material_type": material,
                    "official_tissue_type": tissue_type,
                    "official_tissue": official_tissue,
                    "official_tissue_source": official_tissue_source,
                }
            )
        studies[accession] = {
            "accession": accession,
            "title": clean_value(study.get("title", "")),
            "description": clean_value(study.get("description", "")),
            "project_identifier": clean_value(study.get("projectIdentifier", "")),
            "url": f"https://osdr.nasa.gov/bio/repo/data/studies/{accession}",
            "api_url": f"https://osdr.nasa.gov/geode-py/ws/repo/studies/{accession}",
            "material_field_present": bool(material_field),
            "official_sample_count": int(page.get("totalRecords", len(rows))),
        }
        sample_rows[accession] = rows
    return studies, sample_rows


def match_official_sample(local_row, official_rows: list[dict]):
    by_sample = unique_mapping(official_rows, "official_sample_name")
    by_source = unique_mapping(official_rows, "official_source_name")
    by_sample_normalized = unique_normalized_mapping(
        official_rows,
        "official_sample_name",
    )
    by_source_normalized = unique_normalized_mapping(
        official_rows,
        "official_source_name",
    )
    candidates = [
        (
            "profile_to_sample",
            clean_value(local_row.get("profile", "")),
            by_sample,
            by_sample_normalized,
        ),
        (
            "id_sample_to_sample",
            clean_value(local_row.get("id.sample name", "")),
            by_sample,
            by_sample_normalized,
        ),
        (
            "profile_to_source",
            clean_value(local_row.get("profile", "")),
            by_source,
            by_source_normalized,
        ),
        (
            "study_source_to_source",
            clean_value(local_row.get("study.source name", "")),
            by_source,
            by_source_normalized,
        ),
    ]
    accession_sample = clean_value(local_row.get("id.accession_sample name", ""))
    accession = resolved_accession(local_row)
    prefix = f"{accession}_"
    if accession_sample.startswith(prefix):
        candidates.insert(
            2,
            (
                "accession_sample_to_sample",
                accession_sample[len(prefix) :],
                by_sample,
                by_sample_normalized,
            ),
        )
    for method, value, mapping, _ in candidates:
        if value and value in mapping:
            return mapping[value], method
    for method, value, _, normalized_mapping in candidates:
        normalized = normalized_identifier(value)
        if normalized and normalized in normalized_mapping:
            return normalized_mapping[normalized], f"{method}_normalized"
    return None, "unmatched"


def counter_text(values) -> str:
    counts = Counter(value for value in values if value)
    return "; ".join(f"{key} ({counts[key]})" for key in sorted(counts))


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    metadata_dir = Path(args.metadata_dir)
    if args.download:
        download_metadata(
            args.profile_metadata,
            metadata_dir,
            args.download_timeout,
            args.refresh,
        )
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles = pd.read_csv(args.profile_metadata, sep="\t", keep_default_na=False)
    studies, official_samples = load_studies(metadata_dir)

    rows = []
    for _, local in profiles.iterrows():
        accession = resolved_accession(local)
        inferred = infer_tissue(local.get("profile", ""))
        official, match_method = match_official_sample(
            local,
            official_samples.get(accession, []),
        )
        official_material = official["official_material_type"] if official else ""
        official_tissue_type = official["official_tissue_type"] if official else ""
        official_tissue = official["official_tissue"] if official else "unknown"
        if official_tissue not in {"unknown", "unspecified"}:
            final_tissue = official_tissue
            tissue_source = official["official_tissue_source"]
        elif inferred != "unknown":
            final_tissue = inferred
            tissue_source = "sample_name_inference"
        else:
            final_tissue = official_tissue
            tissue_source = (
                "official_osdr_material_type"
                if official_material
                else "unassigned"
            )
        rows.append(
            {
                "profile": clean_value(local.get("profile", "")),
                "id.accession": accession,
                "condition_inferred": clean_value(
                    local.get("condition_inferred", "")
                ),
                "official_sample_name": (
                    official["official_sample_name"] if official else ""
                ),
                "official_source_name": (
                    official["official_source_name"] if official else ""
                ),
                "official_material_type": official_material,
                "official_tissue_type": official_tissue_type,
                "official_tissue": official_tissue,
                "sample_name_tissue_inferred": inferred,
                "official_vs_inferred_match": (
                    official_tissue == inferred
                    if official_tissue not in {"unknown", "unspecified"}
                    and inferred != "unknown"
                    else ""
                ),
                "tissue_final": final_tissue,
                "tissue_source": tissue_source,
                "official_sample_match_method": match_method,
                "official_study_url": studies.get(accession, {}).get("url", ""),
                "official_api_url": studies.get(accession, {}).get("api_url", ""),
            }
        )
    sample_df = pd.DataFrame(rows)
    sample_path = output_dir / "osdr_sample_tissues.tsv"
    sample_df.to_csv(sample_path, sep="\t", index=False)

    accession_rows = []
    for accession, local_group in sample_df.groupby("id.accession", dropna=False):
        if not accession:
            continue
        study = studies.get(accession, {})
        matched = local_group["official_sample_match_method"].ne("unmatched")
        official_material = local_group.loc[
            local_group["official_material_type"].ne(""),
            "official_material_type",
        ]
        official_tissues = local_group.loc[
            ~local_group["official_tissue"].isin(["unknown", "unspecified"]),
            "official_tissue",
        ]
        final_tissues = local_group.loc[
            ~local_group["tissue_final"].isin(["unknown", "unspecified"]),
            "tissue_final",
        ]
        accession_rows.append(
            {
                "id.accession": accession,
                "study_title": study.get("title", ""),
                "project_identifier": study.get("project_identifier", ""),
                "n_local_profiles": int(len(local_group)),
                "n_official_sample_matches": int(matched.sum()),
                "official_match_fraction": float(matched.mean()),
                "official_material_types": counter_text(official_material),
                "official_tissues": counter_text(official_tissues),
                "final_tissues": counter_text(final_tissues),
                "material_field_present": study.get(
                    "material_field_present", False
                ),
                "official_sample_count": study.get("official_sample_count", 0),
                "official_study_url": study.get("url", ""),
                "official_api_url": study.get("api_url", ""),
            }
        )
    accession_df = pd.DataFrame(accession_rows).sort_values("id.accession")
    accession_path = output_dir / "osdr_accession_tissues.tsv"
    accession_df.to_csv(accession_path, sep="\t", index=False)

    tissue_accession_rows = []
    assigned = sample_df[
        ~sample_df["tissue_final"].isin(["unknown", "unspecified"])
    ]
    for tissue, group in assigned.groupby("tissue_final"):
        accessions = sorted(
            accession
            for accession in group["id.accession"].astype(str).unique()
            if accession
        )
        tissue_accession_rows.append(
            {
                "tissue": tissue,
                "n_accessions": len(accessions),
                "n_profiles": int(len(group)),
                "accessions": ",".join(accessions),
            }
        )
    tissue_accession_df = pd.DataFrame(tissue_accession_rows).sort_values(
        ["n_accessions", "n_profiles", "tissue"],
        ascending=[False, False, True],
    )
    tissue_accession_path = output_dir / "osdr_tissue_accessions.tsv"
    tissue_accession_df.to_csv(tissue_accession_path, sep="\t", index=False)

    validation = sample_df[
        sample_df["official_vs_inferred_match"].isin([True, False])
    ].copy()
    confusion = (
        validation.groupby(
            ["official_tissue", "sample_name_tissue_inferred"],
            dropna=False,
        )
        .size()
        .rename("n_profiles")
        .reset_index()
        .sort_values("n_profiles", ascending=False)
    )
    confusion_path = output_dir / "official_vs_inferred_tissue.tsv"
    confusion.to_csv(confusion_path, sep="\t", index=False)

    official_assigned = sample_df["official_tissue"].isin(
        [
            value
            for value in sample_df["official_tissue"].unique()
            if value not in {"unknown", "unspecified"}
        ]
    )
    matched = sample_df["official_sample_match_method"].ne("unmatched")
    summary = {
        "metadata_dir": str(metadata_dir),
        "profile_metadata": str(args.profile_metadata),
        "official_studies_loaded": len(studies),
        "local_profiles": len(sample_df),
        "official_sample_matches": int(matched.sum()),
        "official_sample_match_fraction": float(matched.mean()),
        "official_tissue_assignments": int(official_assigned.sum()),
        "official_tissue_assignment_fraction": float(official_assigned.mean()),
        "final_tissue_assignments": int(
            (~sample_df["tissue_final"].isin(["unknown", "unspecified"])).sum()
        ),
        "validation_profiles": len(validation),
        "sample_name_inference_accuracy_on_officially_labeled_profiles": (
            float(validation["official_vs_inferred_match"].mean())
            if len(validation)
            else None
        ),
        "outputs": {
            "sample_tissues": str(sample_path),
            "accession_tissues": str(accession_path),
            "tissue_accessions": str(tissue_accession_path),
            "official_vs_inferred": str(confusion_path),
        },
    }
    summary_path = output_dir / "osdr_tissue_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(
        f"official_sample_matches={summary['official_sample_matches']}/"
        f"{summary['local_profiles']}"
    )
    print(
        "official_tissue_assignments="
        f"{summary['official_tissue_assignments']}/{summary['local_profiles']}"
    )
    print(f"summary={summary_path}")
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract tissues from official OSDR ISA material and tissue metadata."
        )
    )
    parser.add_argument("--metadata-dir", default=DEFAULT_METADATA_DIR)
    parser.add_argument("--profile-metadata", default=DEFAULT_PROFILE_METADATA)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--download-timeout", type=int, default=120)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
