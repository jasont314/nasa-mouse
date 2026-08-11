"""Retrain the four ASGSR expiMap pipelines under matched random seeds."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .build_asgsr_paper import (
    CONFIGS,
    CURATED_PATHWAYS,
    FIGURE_DIR,
    ROLE_COLORS,
    SOURCE_DIR,
    latent_directions,
)


SEEDS = (2020, 2021, 2022)


def _curated() -> pd.DataFrame:
    frames = []
    for config in CONFIGS:
        frame = pd.DataFrame(
            CURATED_PATHWAYS[config.tissue],
            columns=[
                "term",
                "short_label",
                "evidence_role",
                "paper_interpretation",
                "paper_citations",
            ],
        )
        frame.insert(0, "tissue", config.tissue)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _family_terms() -> pd.DataFrame:
    frame = pd.read_csv(
        SOURCE_DIR / "table_s11_nonredundant_pathway_families.tsv", sep="\t"
    )
    return frame[["tissue", "representative_term"]].rename(
        columns={"representative_term": "term"}
    )


def _seed_path(path: Path, seed: int) -> Path:
    return Path(str(path).replace("seed2020", f"seed{seed}"))


def _run_logged(command: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    source_root = str(Path(__file__).resolve().parents[1])
    environment["PYTHONPATH"] = ":".join(
        value
        for value in (source_root, environment.get("PYTHONPATH", ""))
        if value
    )
    with log_path.open("w", encoding="utf-8") as handle:
        process = subprocess.run(
            command,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            env=environment,
            check=False,
        )
    if process.returncode != 0:
        tail = "\n".join(log_path.read_text(encoding="utf-8").splitlines()[-40:])
        raise RuntimeError(
            f"Command failed with exit code {process.returncode}: {' '.join(command)}\n{tail}"
        )


def ensure_seed_run(config, seed: int) -> tuple[Path, Path]:
    reference_dir = _seed_path(config.reference_summary.parent, seed)
    query_dir = _seed_path(config.run_dir, seed)
    if seed == 2020:
        return reference_dir, query_dir

    reference_summary = reference_dir / "training_summary.json"
    if not reference_summary.exists():
        primary_summary = json.loads(config.reference_summary.read_text(encoding="utf-8"))
        reference_input = str(primary_summary["input"])
        print(f"training reference {config.tissue} seed {seed}", flush=True)
        _run_logged(
            [
                sys.executable,
                "-m",
                "nasa_mouse_expimap.train_expimap_archs4_reference",
                "--input",
                reference_input,
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

    query_summary = query_dir / "query_mapping_summary.json"
    if not query_summary.exists():
        print(f"mapping query {config.tissue} seed {seed}", flush=True)
        _run_logged(
            [
                sys.executable,
                "-m",
                "nasa_mouse_expimap.map_expimap_osdr_query",
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


def _directions_for_run(config, seed: int, query_dir: Path) -> pd.Series:
    if seed == 2020:
        return latent_directions(config).set_index("term")["latent_orientation"]

    import scarches as sca

    mapped = ad.read_h5ad(query_dir / "mapped_query_with_scores.h5ad")
    model = sca.models.EXPIMAP.load(query_dir / "query_model", adata=mapped)
    model.latent_directions(method="sum", adata=model.adata)
    terms = list(map(str, model.adata.uns["terms"]))
    directions = np.asarray(model.adata.uns["directions"], dtype=float)
    return pd.Series(directions, index=terms, name="latent_orientation")


def seed_effects(config, seed: int, query_dir: Path) -> pd.DataFrame:
    scores = pd.read_csv(query_dir / "query_pathway_scores.tsv", sep="\t")
    directions = _directions_for_run(config, seed, query_dir)
    terms = [term for term in directions.index if term in scores and directions[term] != 0]
    oriented = scores[terms].astype(float).mul(directions.loc[terms], axis=1)
    condition = scores["condition_inferred"].astype(str)
    accession_rows = []
    for accession, indexes in scores.groupby("id.accession", observed=True).indices.items():
        indexes = np.asarray(indexes)
        local = condition.iloc[indexes]
        flight = local.eq("flight").to_numpy()
        ground = local.eq("ground_control").to_numpy()
        if not flight.any() or not ground.any():
            continue
        matrix = oriented.iloc[indexes]
        effect = matrix.iloc[flight].mean(axis=0) - matrix.iloc[ground].mean(axis=0)
        accession_rows.extend(
            {
                "tissue": config.tissue,
                "seed": seed,
                "accession": str(accession),
                "term": str(term),
                "effect": float(value),
            }
            for term, value in effect.items()
        )
    accession = pd.DataFrame(accession_rows)
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
    return summary


def add_flags(effects: pd.DataFrame) -> pd.DataFrame:
    curated = _curated()
    family = _family_terms()
    effects = effects.merge(
        curated[["tissue", "term", "short_label", "evidence_role"]],
        on=["tissue", "term"],
        how="left",
    )
    effects = effects.merge(
        family.assign(family_representative=True),
        on=["tissue", "term"],
        how="left",
    )
    effects["curated_pathway"] = effects["evidence_role"].notna()
    effects["family_representative"] = effects["family_representative"].fillna(False)
    primary = effects.loc[effects["seed"].eq(2020), [
        "tissue", "term", "accession_balanced_effect"
    ]].copy()
    primary["primary_absolute_percentile"] = primary.groupby(
        "tissue", observed=True
    )["accession_balanced_effect"].transform(lambda values: values.abs().rank(pct=True))
    return effects.merge(
        primary[["tissue", "term", "primary_absolute_percentile"]],
        on=["tissue", "term"],
        how="left",
    )


def seed_summary(effects: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide = effects.pivot_table(
        index=["tissue", "term"],
        columns="seed",
        values="accession_balanced_effect",
    ).reset_index()
    wide = wide.rename(columns={seed: f"effect_seed{seed}" for seed in SEEDS})
    flags = effects.loc[
        effects["seed"].eq(2020),
        [
            "tissue",
            "term",
            "short_label",
            "evidence_role",
            "curated_pathway",
            "family_representative",
            "primary_absolute_percentile",
        ],
    ]
    wide = wide.merge(flags, on=["tissue", "term"], how="left")
    wide["curated_pathway"] = wide["curated_pathway"].fillna(False).astype(bool)
    wide["family_representative"] = (
        wide["family_representative"].fillna(False).astype(bool)
    )
    seed_columns = [f"effect_seed{seed}" for seed in SEEDS]
    wide["all_three_seeds_available"] = wide[seed_columns].notna().all(axis=1)
    signs = np.sign(wide[seed_columns])
    wide["all_three_direction_concordant"] = signs.nunique(axis=1, dropna=True).eq(1)
    wide["seed_effect_median"] = wide[seed_columns].median(axis=1)
    wide["seed_effect_minimum"] = wide[seed_columns].min(axis=1)
    wide["seed_effect_maximum"] = wide[seed_columns].max(axis=1)

    rows = []
    masks = {
        "all_active": pd.Series(True, index=wide.index),
        "primary_top_decile": wide["primary_absolute_percentile"].ge(0.9).fillna(False),
        "curated": wide["curated_pathway"],
        "family_representatives": wide["family_representative"],
    }
    for pathway_set, mask in masks.items():
        for tissue, frame in wide.loc[mask].groupby("tissue", observed=True):
            complete = frame.loc[frame["all_three_seeds_available"]]
            correlations = []
            directions = []
            for seed in (2021, 2022):
                comparison = f"effect_seed{seed}"
                correlations.append(
                    stats.spearmanr(complete["effect_seed2020"], complete[comparison]).statistic
                )
                directions.append(
                    (
                        np.sign(complete["effect_seed2020"])
                        == np.sign(complete[comparison])
                    ).mean()
                )
            rows.append(
                {
                    "tissue": tissue,
                    "pathway_set": pathway_set,
                    "n_pathways": int(len(complete)),
                    "minimum_seed_vs_primary_spearman_rho": float(min(correlations)),
                    "mean_seed_vs_primary_spearman_rho": float(np.mean(correlations)),
                    "minimum_seed_vs_primary_direction_agreement": float(min(directions)),
                    "mean_seed_vs_primary_direction_agreement": float(np.mean(directions)),
                    "all_three_seed_direction_agreement": float(
                        complete["all_three_direction_concordant"].mean()
                    ),
                }
            )
    return wide, pd.DataFrame(rows)


def plot_seed_stability(consensus: pd.DataFrame, summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
    for ax, config in zip(axes.flat, CONFIGS):
        frame = consensus.loc[
            consensus["tissue"].eq(config.tissue)
            & consensus["all_three_seeds_available"]
        ]
        ax.vlines(
            frame["effect_seed2020"],
            frame["seed_effect_minimum"],
            frame["seed_effect_maximum"],
            color="#c4c8ca",
            alpha=0.45,
            lw=0.7,
        )
        ax.scatter(
            frame["effect_seed2020"],
            frame["seed_effect_median"],
            s=17,
            color="#a3aaad",
            alpha=0.55,
            linewidth=0,
        )
        selected = frame.loc[frame["curated_pathway"]]
        ax.scatter(
            selected["effect_seed2020"],
            selected["seed_effect_median"],
            s=55,
            c=selected["evidence_role"].map(ROLE_COLORS),
            edgecolor="white",
            linewidth=0.7,
        )
        bounds = np.nanmax(
            np.abs(
                frame[[
                    "effect_seed2020",
                    "seed_effect_median",
                    "seed_effect_minimum",
                    "seed_effect_maximum",
                ]].to_numpy()
            )
        )
        ax.plot([-bounds, bounds], [-bounds, bounds], color="#596064", lw=0.8)
        ax.axhline(0, color="#7b8083", lw=0.6)
        ax.axvline(0, color="#7b8083", lw=0.6)
        row = summary.loc[
            summary["tissue"].eq(config.tissue)
            & summary["pathway_set"].eq("curated")
        ].iloc[0]
        ax.set_title(
            f"{config.display_name} | selected directions stable: "
            f"{row.all_three_seed_direction_agreement:.0%}"
        )
        ax.set_xlabel("Seed 2020 study-balanced shift")
        ax.set_ylabel("Median shift across seeds")
    fig.suptitle(
        "Full reference and query training seed sensitivity",
        fontsize=16,
        fontweight="bold",
    )
    fig.savefig(FIGURE_DIR / "figure_s6_training_seed_sensitivity.png", dpi=300)
    fig.savefig(FIGURE_DIR / "figure_s6_training_seed_sensitivity.pdf")
    plt.close(fig)


def training_manifest() -> pd.DataFrame:
    rows = []
    for config in CONFIGS:
        for seed in SEEDS:
            reference_dir = _seed_path(config.reference_summary.parent, seed)
            query_dir = _seed_path(config.run_dir, seed)
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
                    "reference_input": reference["input"],
                    "reference_output": str(reference_dir),
                    "reference_epochs_completed": reference["training"]["epochs_completed"],
                    "reference_best_epoch": reference["training"]["best_epoch"],
                    "reference_training_seconds": reference["training"]["training_seconds"],
                    "reference_reconstruction_loss": reference["recon_loss"],
                    "query_input": query["query_h5ad"],
                    "query_output": str(query_dir),
                    "query_epochs": query["epochs"],
                    "gpu": query["torch"]["cuda_device_name"],
                    "posterior_mean_latent": query["posterior_mean_latent"],
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    frames = []
    for config in CONFIGS:
        for seed in SEEDS:
            _, query_dir = ensure_seed_run(config, seed)
            print(f"summarizing {config.tissue} seed {seed}", flush=True)
            frames.append(seed_effects(config, seed, query_dir))
    effects = add_flags(pd.concat(frames, ignore_index=True))
    consensus, summary = seed_summary(effects)
    manifest = training_manifest()
    effects.to_csv(
        SOURCE_DIR / "table_s20_training_seed_pathway_effects.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    consensus.to_csv(
        SOURCE_DIR / "table_s21_training_seed_consensus.tsv", sep="\t", index=False
    )
    summary.to_csv(
        SOURCE_DIR / "table_s22_training_seed_summary.tsv", sep="\t", index=False
    )
    manifest.to_csv(
        SOURCE_DIR / "table_s23_training_seed_manifest.tsv", sep="\t", index=False
    )
    plot_seed_stability(consensus, summary)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
