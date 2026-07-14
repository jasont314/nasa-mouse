"""Retrain the kidney and spleen HVG reference-query models across three seeds."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from .build_asgsr_paper import ROOT
from .run_asgsr_seed_sensitivity import _run_logged


SEEDS = (2020, 2021, 2022)
OUTPUT_DIR = ROOT / "outputs/expimap_kidney_spleen_reassessment"
PROJECT_COLUMN = "investigation.study.comment.project identifier"


@dataclass(frozen=True)
class ReassessmentConfig:
    tissue: str
    base_dir: Path

    @property
    def reference_input(self) -> Path:
        return (
            self.base_dir
            / "reassessment_hvg_2000/input"
            / f"archs4_mouse_{self.tissue}_reference_tutorial_hvg_raw_counts.h5ad"
        )

    @property
    def query_input(self) -> Path:
        return (
            self.base_dir
            / "reassessment_hvg_2000/input"
            / f"osdr_{self.tissue}_query_tutorial_hvg_raw_counts.h5ad"
        )

    def reference_dir(self, seed: int) -> Path:
        return self.base_dir / "reassessment_hvg_2000" / f"reference_nb_400epoch_seed{seed}"

    def query_dir(self, seed: int) -> Path:
        return self.base_dir / "reassessment_hvg_2000" / f"query_nb_250epoch_seed{seed}"


CONFIGS = tuple(
    ReassessmentConfig(
        tissue=tissue,
        base_dir=ROOT / f"outputs/expimap_archs4_reference_osdr_query_{tissue}",
    )
    for tissue in ("kidney", "spleen")
)


def ensure_seed_run(config: ReassessmentConfig, seed: int) -> tuple[Path, Path]:
    reference_dir = config.reference_dir(seed)
    query_dir = config.query_dir(seed)
    if not (reference_dir / "training_summary.json").exists():
        print(f"training reference {config.tissue} seed {seed}", flush=True)
        _run_logged(
            [
                sys.executable,
                "-m",
                "expiMap_scarches.nasa_mouse_expimap.train_expimap_archs4_reference",
                "--input",
                str(config.reference_input),
                "--output-dir",
                str(reference_dir),
                "--condition-key",
                "archs4_condition",
                "--recon-loss",
                "nb",
                "--epochs",
                "400",
                "--hidden-layers",
                "300,300,300",
                "--seed",
                str(seed),
                "--early-stopping",
                "--early-stopping-patience",
                "50",
                "--batch-size",
                "128",
                "--train-frac",
                "0.9",
            ],
            reference_dir / "seed_sensitivity_training.log",
        )

    if not (query_dir / "query_mapping_summary.json").exists():
        print(f"mapping query {config.tissue} seed {seed}", flush=True)
        _run_logged(
            [
                sys.executable,
                "-m",
                "expiMap_scarches.nasa_mouse_expimap.map_expimap_osdr_query",
                "--reference-model",
                str(reference_dir / "model"),
                "--query-h5ad",
                str(config.query_input),
                "--output-dir",
                str(query_dir),
                "--reference-condition-key",
                "archs4_condition",
                "--query-condition-source",
                "id.accession",
                "--recon-loss",
                "nb",
                "--epochs",
                "250",
                "--seed",
                str(seed),
                "--mean-latent",
            ],
            query_dir / "seed_sensitivity_query.log",
        )
    return reference_dir, query_dir


def latent_directions(query_dir: Path) -> pd.Series:
    import scarches as sca

    mapped = ad.read_h5ad(query_dir / "mapped_query_with_scores.h5ad")
    model = sca.models.EXPIMAP.load(query_dir / "query_model", adata=mapped)
    model.latent_directions(method="sum", adata=model.adata)
    terms = list(map(str, model.adata.uns["terms"]))
    directions = np.asarray(model.adata.uns["directions"], dtype=float)
    return pd.Series(directions, index=terms, name="latent_orientation")


def seed_effects(
    config: ReassessmentConfig, seed: int, query_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.read_csv(query_dir / "query_pathway_scores.tsv", sep="\t")
    directions = latent_directions(query_dir)
    terms = [term for term in directions.index if term in scores and directions[term] != 0]
    oriented = scores[terms].astype(float).mul(directions.loc[terms], axis=1)
    condition = scores["condition_inferred"].astype(str)
    rows = []
    for accession, indexes in scores.groupby("id.accession", observed=True).indices.items():
        indexes = np.asarray(indexes)
        local = condition.iloc[indexes]
        flight = local.eq("flight").to_numpy()
        ground = local.eq("ground_control").to_numpy()
        if not flight.any() or not ground.any():
            continue
        effect = (
            oriented.iloc[indexes].iloc[flight].mean(axis=0)
            - oriented.iloc[indexes].iloc[ground].mean(axis=0)
        )
        project_values = scores.iloc[indexes][PROJECT_COLUMN].dropna().astype(str)
        project = project_values.mode().iloc[0] if not project_values.empty else str(accession)
        rows.extend(
            {
                "tissue": config.tissue,
                "seed": seed,
                "accession": str(accession),
                "project": project,
                "term": str(term),
                "effect": float(value),
            }
            for term, value in effect.items()
        )
    accession = pd.DataFrame(rows)
    project = (
        accession.groupby(["tissue", "seed", "project", "term"], as_index=False)[
            "effect"
        ]
        .mean()
        .sort_values(["tissue", "seed", "term", "project"])
    )
    summary = (
        accession.groupby(["tissue", "seed", "term"])["effect"]
        .agg(
            accession_balanced_effect="mean",
            n_accessions="size",
            accessions_positive=lambda values: int((values > 0).sum()),
            accessions_negative=lambda values: int((values < 0).sum()),
        )
        .reset_index()
    )
    project_summary = (
        project.groupby(["tissue", "seed", "term"])["effect"]
        .agg(
            project_balanced_effect="mean",
            n_projects="size",
            projects_positive=lambda values: int((values > 0).sum()),
            projects_negative=lambda values: int((values < 0).sum()),
        )
        .reset_index()
    )
    return accession, summary.merge(
        project_summary, on=["tissue", "seed", "term"], how="left"
    )


def consensus_table(effects: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide = effects.pivot_table(
        index=["tissue", "term"],
        columns="seed",
        values="accession_balanced_effect",
    ).reset_index()
    wide = wide.rename(columns={seed: f"effect_seed{seed}" for seed in SEEDS})
    seed_columns = [f"effect_seed{seed}" for seed in SEEDS]
    wide["all_three_seeds_available"] = wide[seed_columns].notna().all(axis=1)
    signs = np.sign(wide[seed_columns])
    wide["all_three_direction_concordant"] = (
        wide["all_three_seeds_available"] & signs.nunique(axis=1, dropna=True).eq(1)
    )
    wide["seed_effect_median"] = wide[seed_columns].median(axis=1)
    wide["seed_effect_minimum"] = wide[seed_columns].min(axis=1)
    wide["seed_effect_maximum"] = wide[seed_columns].max(axis=1)
    wide["primary_absolute_percentile"] = wide.groupby("tissue", observed=True)[
        "effect_seed2020"
    ].transform(lambda values: values.abs().rank(pct=True))

    rows = []
    for tissue, frame in wide.groupby("tissue", observed=True):
        complete = frame.loc[frame["all_three_seeds_available"]]
        for pathway_set, subset in (
            ("all_active", complete),
            (
                "primary_top_decile",
                complete.loc[complete["primary_absolute_percentile"].ge(0.9)],
            ),
        ):
            correlations = [
                stats.spearmanr(
                    subset["effect_seed2020"], subset[f"effect_seed{seed}"]
                ).statistic
                for seed in (2021, 2022)
            ]
            directions = [
                (
                    np.sign(subset["effect_seed2020"])
                    == np.sign(subset[f"effect_seed{seed}"])
                ).mean()
                for seed in (2021, 2022)
            ]
            rows.append(
                {
                    "tissue": tissue,
                    "pathway_set": pathway_set,
                    "n_pathways": int(len(subset)),
                    "minimum_seed_vs_primary_spearman_rho": float(min(correlations)),
                    "mean_seed_vs_primary_spearman_rho": float(np.mean(correlations)),
                    "minimum_seed_vs_primary_direction_agreement": float(min(directions)),
                    "mean_seed_vs_primary_direction_agreement": float(np.mean(directions)),
                    "all_three_seed_direction_agreement": float(
                        subset["all_three_direction_concordant"].mean()
                    ),
                }
            )
    return wide, pd.DataFrame(rows)


def training_manifest() -> pd.DataFrame:
    rows = []
    for config in CONFIGS:
        input_manifest = json.loads(
            (
                config.base_dir
                / "reassessment_hvg_2000/input/tutorial_hvg_input_manifest.json"
            ).read_text(encoding="utf-8")
        )
        for seed in SEEDS:
            reference_dir = config.reference_dir(seed)
            query_dir = config.query_dir(seed)
            reference = json.loads(
                (reference_dir / "training_summary.json").read_text(encoding="utf-8")
            )
            query = json.loads(
                (query_dir / "query_mapping_summary.json").read_text(encoding="utf-8")
            )
            rows.append(
                {
                    "tissue": config.tissue,
                    "seed": seed,
                    "reference_samples": input_manifest["n_reference_samples"],
                    "reference_series": reference["n_conditions"],
                    "reference_epochs_completed": reference["training"]["epochs_completed"],
                    "reference_best_epoch": reference["training"]["best_epoch"],
                    "reference_training_seconds": reference["training"]["training_seconds"],
                    "query_samples": query["n_query_samples"],
                    "query_epochs": query["epochs"],
                    "gpu": query["torch"]["cuda_device_name"],
                    "hvg_method": input_manifest["hvg_method"],
                    "hvg_excluded_samples": input_manifest["hvg_excluded_samples"],
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    accession_frames = []
    effect_frames = []
    for config in CONFIGS:
        for seed in SEEDS:
            _, query_dir = ensure_seed_run(config, seed)
            print(f"summarizing {config.tissue} seed {seed}", flush=True)
            accession, effects = seed_effects(config, seed, query_dir)
            accession_frames.append(accession)
            effect_frames.append(effects)

    accessions = pd.concat(accession_frames, ignore_index=True)
    effects = pd.concat(effect_frames, ignore_index=True)
    consensus, summary = consensus_table(effects)
    manifest = training_manifest()
    accessions.to_csv(
        OUTPUT_DIR / "seed_accession_effects.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    effects.to_csv(OUTPUT_DIR / "seed_pathway_effects.tsv", sep="\t", index=False)
    consensus.to_csv(OUTPUT_DIR / "seed_consensus.tsv", sep="\t", index=False)
    summary.to_csv(OUTPUT_DIR / "seed_summary.tsv", sep="\t", index=False)
    manifest.to_csv(OUTPUT_DIR / "seed_training_manifest.tsv", sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
