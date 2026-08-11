"""Inspect all ARCHS4 mouse profiles and build balanced tissue references."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import h5py
import numpy as np
import pandas as pd

from .config import DataConfig
from .tissues import NONMATCHABLE_OSDR_CLASSES, rules_for_tissues


METADATA_FIELDS = (
    "geo_accession",
    "series_id",
    "title",
    "source_name_ch1",
    "characteristics_ch1",
    "library_strategy",
    "library_source",
    "instrument_model",
    "organism_ch1",
    "singlecellprobability",
)
LEAKAGE_PATTERN = re.compile(
    r"\b(?:nasa|genelab|glds|osd[- _]?\d+|space\s*flight|microgravity|"
    r"international space station|iss|rodent research|rr[- _]?\d+|hindlimb unloading)\b",
    flags=re.IGNORECASE,
)
HEALTH_PATTERNS = {
    "disease_or_tumor": re.compile(
        r"\b(?:tumou?r|cancer|carcinoma|adenoma|sarcoma|leukemia|lymphoma|"
        r"metasta\w*|fibrosis|cirrhosis|disease|diabet\w*|obes\w*|infect\w*|"
        r"inflamm\w*|atherosclero\w*)\b",
        flags=re.IGNORECASE,
    ),
    "genetic_perturbation": re.compile(
        r"\b(?:knock out|knockout|knock down|knockdown|crispr|mutant|mutation|"
        r"mut\s*\d*|transgenic|deficien\w*|null|overexpress\w*|conditional ko|"
        r"cko|dko|tko|qko|ko|k o|[a-z0-9]*ko\d*|"
        r"[a-z0-9]*(?:flox|flx)[a-z0-9]*|cre|[a-z0-9]*creert[a-z0-9]*|"
        r"e?knockout|heterozyg\w*|homozyg\w*|transduced)\b",
        flags=re.IGNORECASE,
    ),
    "experimental_intervention": re.compile(
        r"\b(?:treated|treatment|drug|compound|dose|dosed|irradiat\w*|radiation|"
        r"injur\w*|surgery|resection|exposure|exposed|toxin|lipopolysaccharide|"
        r"lps|oxldl|high fat|hfd|fasted|fasting|exercise|stressed|stress model|"
        r"[a-z0-9]*sirna[a-z0-9]*|[a-z0-9]*shrna[a-z0-9]*|dox|injected|"
        r"injection|stimulated|activated|organoids?|in vitro|"
        r"transfect\w*|tamoxifen|vehicle|saline|agent)\b",
        flags=re.IGNORECASE,
    ),
    "developmental_mismatch": re.compile(
        r"\b(?:embryo\w*|embryonic|fetal|foetal|neonatal|newborn|postnatal day|"
        r"postnatal|post natal|embryoid|embroid|placenta|zygote|oocyte|"
        r"e\s*\d+(?:\s*\d+)?|p\s*\d+(?:\s*\d+)?)\b",
        flags=re.IGNORECASE,
    ),
}
CONTROL_PATTERN = re.compile(
    r"\b(?:healthy|normal|wild type|wildtype|wt|untreated|control|vehicle|sham|naive|baseline)\b",
    flags=re.IGNORECASE,
)
SINGLE_CELL_PATTERN = re.compile(
    r"\b(?:single cell|singlecell|scrna(?: seq)?|snrna(?: seq)?|single nucleus|"
    r"single nuclei|10x(?: genomics)?|chromium|drop seq|dropseq|"
    r"smart seq\w*|smartseq\w*)\b",
    flags=re.IGNORECASE,
)


def _decode(values) -> list[str]:
    return [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def load_archs4_metadata(path: str | Path) -> tuple[pd.DataFrame, dict]:
    rows: dict[str, object] = {}
    with h5py.File(path, "r") as handle:
        samples = handle["meta/samples"]
        for field in METADATA_FIELDS:
            values = samples[field][:]
            rows[field] = _decode(values) if values.dtype.kind in {"S", "O", "U"} else values
        expression_shape = tuple(map(int, handle["data/expression"].shape))
    frame = pd.DataFrame(rows)
    frame.insert(0, "archs4_sample_index", np.arange(len(frame), dtype=np.int64))
    metadata = {
        "profiles_inspected": int(len(frame)),
        "expression_shape_genes_by_samples": list(expression_shape),
        "metadata_fields": list(METADATA_FIELDS),
    }
    return frame, metadata


def _normalize_series(values: pd.Series) -> pd.Series:
    return (
        values.fillna("")
        .astype(str)
        .str.lower()
        .str.replace(r"[^a-z0-9]+", " ", regex=True)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )


def _combined_text(frame: pd.DataFrame) -> pd.Series:
    raw = (
        frame["title"].fillna("").astype(str)
        + " | "
        + frame["source_name_ch1"].fillna("").astype(str)
        + " | "
        + frame["characteristics_ch1"].fillna("").astype(str)
    )
    return _normalize_series(raw)


def classify_tissues(frame: pd.DataFrame, target_tissues: list[str]) -> pd.DataFrame:
    result = frame.copy()
    result["canonical_tissue"] = "unmatched"
    result["tissue_match_source"] = ""
    rules = rules_for_tissues(target_tissues)
    for source in ("source_name_ch1", "title", "characteristics_ch1"):
        values = _normalize_series(result[source])
        for rule in rules:
            available = result["canonical_tissue"].eq("unmatched")
            matched = available & values.str.contains(rule.pattern, na=False)
            if matched.any():
                result.loc[matched, "canonical_tissue"] = rule.canonical
                result.loc[matched, "tissue_match_source"] = source
    return result


def classify_eligibility(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    text = _combined_text(result)
    result["leakage_excluded"] = text.str.contains(LEAKAGE_PATTERN, na=False)
    health_reason = pd.Series("", index=result.index, dtype="object")
    adverse = pd.Series(False, index=result.index)
    for reason, pattern in HEALTH_PATTERNS.items():
        matched = text.str.contains(pattern, na=False)
        adverse |= matched
        empty = health_reason.eq("") & matched
        health_reason.loc[empty] = reason
    explicit_control = text.str.contains(CONTROL_PATTERN, na=False)
    result["health_status"] = np.where(
        adverse,
        "explicit_nonhealthy_or_perturbed",
        np.where(explicit_control, "explicit_control_like", "metadata_unknown"),
    )
    result["health_exclusion_reason"] = health_reason
    result["mouse_organism"] = result["organism_ch1"].str.contains(
        r"mus musculus", case=False, regex=True, na=False
    )
    result["rna_seq"] = result["library_strategy"].str.contains(
        r"rna[- ]?seq", case=False, regex=True, na=False
    )
    result["transcriptomic_library_source"] = _normalize_series(
        result["library_source"]
    ).eq("transcriptomic")
    result["explicit_single_cell"] = text.str.contains(SINGLE_CELL_PATTERN, na=False)
    result["bulk_like"] = (
        pd.to_numeric(result["singlecellprobability"], errors="coerce")
        .fillna(0.0)
        .lt(0.5)
        & result["transcriptomic_library_source"]
        & ~result["explicit_single_cell"]
    )
    common = (
        result["mouse_organism"]
        & result["rna_seq"]
        & result["bulk_like"]
        & ~result["leakage_excluded"]
        & result["canonical_tissue"].ne("unmatched")
    )
    result["eligible_broad"] = common
    result["eligible_healthy_preferred"] = common & result["health_status"].ne(
        "explicit_nonhealthy_or_perturbed"
    )
    result["eligible_control_only"] = common & result["health_status"].eq(
        "explicit_control_like"
    )
    return result


def _stable_key(seed: int, value: object) -> int:
    payload = f"{seed}:{value}".encode("utf-8", "replace")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big", signed=False)


def select_balanced(
    frame: pd.DataFrame,
    *,
    eligibility_column: str,
    max_per_tissue: int,
    max_per_series: int,
    seed: int,
) -> pd.DataFrame:
    selected = frame.loc[frame[eligibility_column]].copy()
    selected = selected.drop_duplicates("geo_accession", keep="first")
    selected["health_priority"] = selected["health_status"].map(
        {"explicit_control_like": 0, "metadata_unknown": 1, "explicit_nonhealthy_or_perturbed": 2}
    ).fillna(3)
    selected["stable_random_key"] = selected["geo_accession"].map(
        lambda value: _stable_key(seed, value)
    )
    selected = selected.sort_values(
        ["canonical_tissue", "series_id", "health_priority", "stable_random_key"],
        kind="stable",
    )
    if max_per_series > 0:
        selected = selected.loc[
            selected.groupby(["canonical_tissue", "series_id"], dropna=False).cumcount()
            < max_per_series
        ].copy()

    selected["within_series_rank"] = selected.groupby(
        ["canonical_tissue", "series_id"], dropna=False
    ).cumcount()
    selected["series_random_key"] = selected["series_id"].map(
        lambda value: _stable_key(seed + 1, value)
    )
    selected = selected.sort_values(
        [
            "canonical_tissue",
            "within_series_rank",
            "series_random_key",
            "health_priority",
            "stable_random_key",
        ],
        kind="stable",
    )
    if max_per_tissue > 0:
        selected = selected.loc[
            selected.groupby("canonical_tissue", dropna=False).cumcount() < max_per_tissue
        ].copy()
    selected["selection_rank_within_tissue"] = selected.groupby(
        "canonical_tissue", dropna=False
    ).cumcount()

    tissue_studies = selected.groupby("canonical_tissue")["series_id"].transform("nunique")
    study_sizes = selected.groupby(["canonical_tissue", "series_id"])["geo_accession"].transform(
        "size"
    )
    n_tissues = max(int(selected["canonical_tissue"].nunique()), 1)
    selected["hierarchical_sampling_weight"] = (
        1.0 / n_tissues / tissue_studies.clip(lower=1) / study_sizes.clip(lower=1)
    )
    return selected.sort_values(
        ["canonical_tissue", "series_id", "geo_accession"], kind="stable"
    ).reset_index(drop=True)


def _cohort_summary(
    classified: pd.DataFrame, selected: pd.DataFrame, target_tissues: list[str], eligibility: str
) -> pd.DataFrame:
    candidate = classified.loc[classified[eligibility]]
    candidate_counts = candidate.groupby("canonical_tissue").agg(
        candidate_samples=("geo_accession", "nunique"),
        candidate_series=("series_id", "nunique"),
        explicit_control_like=(
            "health_status",
            lambda values: int(values.eq("explicit_control_like").sum()),
        ),
        metadata_unknown=("health_status", lambda values: int(values.eq("metadata_unknown").sum())),
    )
    selected_counts = selected.groupby("canonical_tissue").agg(
        selected_samples=("geo_accession", "nunique"),
        selected_series=("series_id", "nunique"),
        sampling_weight_sum=("hierarchical_sampling_weight", "sum"),
    )
    summary = pd.DataFrame(index=pd.Index(target_tissues, name="canonical_tissue"))
    summary = summary.join(candidate_counts).join(selected_counts).fillna(0).reset_index()
    for column in (
        "candidate_samples",
        "candidate_series",
        "explicit_control_like",
        "metadata_unknown",
        "selected_samples",
        "selected_series",
    ):
        summary[column] = summary[column].astype(int)
    summary["archs4_reference_available"] = summary["selected_samples"] > 0
    summary["nonmatchable_reason"] = np.where(
        summary["canonical_tissue"].isin(NONMATCHABLE_OSDR_CLASSES),
        "generic OSDR material label has no defensible ARCHS4 tissue match",
        "none",
    )
    return summary


def run(args: argparse.Namespace) -> Path:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    inventory = pd.read_csv(args.osdr_inventory, sep="\t")
    target_tissues = inventory["tissue_canonical"].astype(str).tolist()
    frame, source_metadata = load_archs4_metadata(args.archs4_h5)
    classified = classify_eligibility(classify_tissues(frame, target_tissues))

    audit_columns = [
        "archs4_sample_index",
        "geo_accession",
        "series_id",
        "title",
        "source_name_ch1",
        "characteristics_ch1",
        "library_strategy",
        "library_source",
        "instrument_model",
        "organism_ch1",
        "singlecellprobability",
        "canonical_tissue",
        "tissue_match_source",
        "health_status",
        "health_exclusion_reason",
        "leakage_excluded",
        "mouse_organism",
        "rna_seq",
        "transcriptomic_library_source",
        "explicit_single_cell",
        "bulk_like",
        "eligible_broad",
        "eligible_healthy_preferred",
        "eligible_control_only",
    ]
    if not args.skip_full_audit:
        classified[audit_columns].to_csv(
            output_dir / "archs4_full_profile_audit.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    outputs: dict[str, str] = {}
    cohort_summaries = []
    for cohort, eligibility in (
        ("control_only", "eligible_control_only"),
        ("healthy_preferred", "eligible_healthy_preferred"),
        ("broad", "eligible_broad"),
    ):
        selected = select_balanced(
            classified,
            eligibility_column=eligibility,
            max_per_tissue=args.max_per_tissue,
            max_per_series=args.max_per_series,
            seed=args.seed,
        )
        selected_path = output_dir / f"archs4_{cohort}_balanced.tsv.gz"
        selected[audit_columns + [
            "selection_rank_within_tissue",
            "hierarchical_sampling_weight",
        ]].to_csv(selected_path, sep="\t", index=False, compression="gzip")
        summary = _cohort_summary(classified, selected, target_tissues, eligibility)
        summary.insert(0, "cohort", cohort)
        summary_path = output_dir / f"archs4_{cohort}_summary.tsv"
        summary.to_csv(summary_path, sep="\t", index=False)
        outputs[f"{cohort}_selected"] = str(selected_path)
        outputs[f"{cohort}_summary"] = str(summary_path)
        cohort_summaries.append(summary)

    combined_summary = pd.concat(cohort_summaries, ignore_index=True)
    combined_path = output_dir / "archs4_reference_summary.tsv"
    combined_summary.to_csv(combined_path, sep="\t", index=False)
    exclusion_counts = {
        "unmatched_tissue": int(classified["canonical_tissue"].eq("unmatched").sum()),
        "not_mouse": int((~classified["mouse_organism"]).sum()),
        "not_rna_seq": int((~classified["rna_seq"]).sum()),
        "single_cell_like": int((~classified["bulk_like"]).sum()),
        "explicit_single_cell": int(classified["explicit_single_cell"].sum()),
        "spaceflight_leakage": int(classified["leakage_excluded"].sum()),
        "explicit_nonhealthy_or_perturbed": int(
            classified["health_status"].eq("explicit_nonhealthy_or_perturbed").sum()
        ),
    }
    manifest = {
        "source": str(args.archs4_h5),
        **source_metadata,
        "target_osdr_tissues": target_tissues,
        "selection": {
            "bulk_like": (
                "singlecellprobability < 0.5, RNA-Seq strategy, transcriptomic "
                "library source, and no explicit scRNA/snRNA/10x/Drop-seq/"
                "Chromium/Smart-seq terms"
            ),
            "leakage_excluded": "NASA/GeneLab/spaceflight/microgravity/ISS/RR/hindlimb terms",
            "healthy_preferred": (
                "exclude explicit disease, tumor, genetic perturbation, intervention, "
                "and developmental-mismatch terms; retain control-like and unknown metadata"
            ),
            "control_only": (
                "healthy-preferred eligibility restricted to metadata explicitly labeled "
                "healthy, normal, wild type, untreated, control, sham, naive, or baseline"
            ),
            "max_per_tissue": args.max_per_tissue,
            "max_per_series": args.max_per_series,
            "seed": args.seed,
            "training_sampler": "uniform tissue, then uniform GEO series, then uniform sample",
        },
        "exclusion_flag_counts_nonexclusive": exclusion_counts,
        "outputs": {
            **outputs,
            "combined_summary": str(combined_path),
            "full_profile_audit": (
                "" if args.skip_full_audit else str(output_dir / "archs4_full_profile_audit.tsv.gz")
            ),
        },
    }
    manifest_path = output_dir / "archs4_catalog_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return combined_path


def parse_args() -> argparse.Namespace:
    defaults = DataConfig()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archs4-h5", default=defaults.archs4_h5)
    parser.add_argument(
        "--osdr-inventory",
        default="outputs/generative/benchmark/data_audit/osdr/osdr_tissue_inventory.tsv",
    )
    parser.add_argument("--output-dir", default="outputs/generative/benchmark/data_audit/archs4")
    parser.add_argument("--max-per-tissue", type=int, default=10000)
    parser.add_argument("--max-per-series", type=int, default=100)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--skip-full-audit", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
