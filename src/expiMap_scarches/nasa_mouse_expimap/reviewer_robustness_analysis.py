"""Run method, held-out-project, and bulk-composition sensitivity analyses.

The analyses in this module are deliberately matched to the frozen four-tissue
ASGSR inputs. They are intended to test whether the reported expiMap directions
depend on the latent model, reuse of the same projects for pathway ranking, or
broad shifts in cell-composition proxies.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable

import anndata as ad
import gseapy as gp
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .build_asgsr_paper import (
    CONFIGS,
    CURATED_PATHWAYS,
    FIGURE_DIR,
    PAPER_DIR,
    ROLE_COLORS,
    SOURCE_DIR,
    TISSUE_COLORS,
    latent_directions,
)


TMS_PATH = (
    Path(__file__).resolve().parents[3]
    / "assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad"
)

ATLAS_TISSUE = {
    "thymus": "thymus",
    "skin": "skin of body",
    "liver": "liver",
    "soleus": "limb muscle",
    "kidney": "kidney",
    "spleen": "spleen",
}

BROAD_COMPARTMENTS = {
    "thymus": {
        "thymocyte": ("DN4 thymocyte", "thymocyte"),
        "myeloid": ("macrophage",),
        "stromal": (
            "epithelial cell of thymus",
            "fibroblast",
            "endothelial cell",
        ),
    },
    "skin": {
        "epidermal": (
            "keratinocyte stem cell",
            "basal cell of epidermis",
            "epidermal cell",
        ),
        "immune": ("T cell", "macrophage"),
        "stromal": ("fibroblast",),
    },
    "liver": {
        "hepatocyte": ("hepatocyte",),
        "endothelial": ("endothelial cell of hepatic sinusoid",),
        "myeloid": ("Kupffer cell", "myeloid leukocyte", "neutrophil"),
        "lymphoid": (
            "B cell",
            "mature NK T cell",
            "CD8-positive, alpha-beta T cell",
            "T cell",
            "CD4-positive, alpha-beta T cell",
            "natural killer cell",
        ),
    },
    "soleus": {
        "satellite": ("skeletal muscle satellite cell",),
        "stromal": ("mesenchymal stem cell",),
        "endothelial": ("endothelial cell",),
        "myeloid": ("macrophage",),
        "lymphoid": ("B cell", "T cell"),
    },
    "kidney": {
        "epithelial": (
            "kidney collecting duct epithelial cell",
            "epithelial cell of proximal tubule",
            "kidney collecting duct principal cell",
            "kidney loop of Henle ascending limb epithelial cell",
        ),
        "endothelial": ("fenestrated endothelial cell",),
        "stromal": ("mesangial cell", "kidney interstitial fibroblast"),
        "myeloid": ("macrophage",),
        "lymphoid": ("B cell", "T cell"),
    },
    "spleen": {
        "b_cell": ("B cell",),
        "t_cell": (
            "CD4-positive, alpha-beta T cell",
            "CD8-positive, alpha-beta T cell",
        ),
        "natural_killer": ("natural killer cell",),
        "erythroid": ("proerythroblast",),
        "granulocyte": ("granulocyte",),
    },
}

PROJECT_COLUMN = "investigation.study.comment.project identifier"
MIN_GENE_SET_SIZE = 5
MAX_GENE_SET_SIZE = 500
N_MARKERS_PER_COMPARTMENT = 30
RANDOM_SEED = 2026


@dataclass
class TissueData:
    tissue: str
    obs: pd.DataFrame
    genes: np.ndarray
    log2cpm: np.ndarray
    gene_sets: dict[str, list[str]]
    expimap_scores: pd.DataFrame


def _dense(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=float)


def _log2cpm(counts: np.ndarray) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals <= 0] = 1.0
    return np.log2(counts / totals * 1_000_000.0 + 1.0)


def _curated_frame(tissue: str) -> pd.DataFrame:
    return pd.DataFrame(
        CURATED_PATHWAYS[tissue],
        columns=[
            "term",
            "short_label",
            "evidence_role",
            "paper_interpretation",
            "paper_citations",
        ],
    )


def _family_terms() -> pd.DataFrame:
    path = SOURCE_DIR / "table_s11_nonredundant_pathway_families.tsv"
    frame = pd.read_csv(path, sep="\t")
    return frame[["tissue", "representative_term"]].rename(
        columns={"representative_term": "term"}
    )


def load_tissue_data(config) -> TissueData:
    adata = ad.read_h5ad(config.query_input)
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    counts = _dense(counts)
    genes = adata.var_names.astype(str).to_numpy()
    terms = list(map(str, adata.uns["terms"]))
    mask = np.asarray(adata.varm["I"])
    if mask.shape != (len(genes), len(terms)):
        raise RuntimeError(f"Unexpected pathway mask shape for {config.tissue}")
    gene_sets = {
        term: genes[np.flatnonzero(mask[:, index] != 0)].tolist()
        for index, term in enumerate(terms)
    }

    score_frame = pd.read_csv(config.run_dir / "query_pathway_scores.tsv", sep="\t")
    direction = latent_directions(config).set_index("term")["latent_orientation"]
    active_terms = [term for term in terms if float(direction.get(term, 0.0)) != 0.0]
    scores = score_frame[active_terms].astype(float).copy()
    scores = scores.mul(direction.loc[active_terms].astype(float), axis=1)
    scores.index = score_frame["obs_name"].astype(str)

    obs = adata.obs.copy()
    obs.index = obs.index.astype(str)
    if list(scores.index) != list(obs.index):
        scores = scores.reindex(obs.index)
    if scores.isna().any().any():
        raise RuntimeError(f"Could not align expiMap scores for {config.tissue}")
    return TissueData(
        tissue=config.tissue,
        obs=obs,
        genes=genes,
        log2cpm=_log2cpm(counts),
        gene_sets=gene_sets,
        expimap_scores=scores,
    )


def _project_lookup(obs: pd.DataFrame) -> dict[str, str]:
    lookup = {}
    for accession, frame in obs.groupby("id.accession", observed=True):
        values = frame[PROJECT_COLUMN].astype(str)
        values = values.loc[~values.isin({"", "nan", "None"})]
        project = values.mode().iloc[0] if not values.empty else str(accession)
        lookup[str(accession)] = str(project)
    return lookup


def accession_and_project_effects(
    scores: pd.DataFrame,
    obs: pd.DataFrame,
    tissue: str,
    method: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scores = scores.reindex(obs.index)
    condition = obs["condition_inferred"].astype(str)
    accession_rows = []
    for accession, indexes in obs.groupby("id.accession", observed=True).indices.items():
        indexes = np.asarray(indexes)
        local_condition = condition.iloc[indexes]
        flight = local_condition.eq("flight").to_numpy()
        ground = local_condition.eq("ground_control").to_numpy()
        if not flight.any() or not ground.any():
            continue
        matrix = scores.iloc[indexes]
        effects = matrix.iloc[flight].mean(axis=0) - matrix.iloc[ground].mean(axis=0)
        accession_rows.extend(
            {
                "tissue": tissue,
                "method": method,
                "accession": str(accession),
                "term": str(term),
                "effect": float(value),
            }
            for term, value in effects.items()
        )
    accession = pd.DataFrame(accession_rows)
    lookup = _project_lookup(obs)
    accession["project"] = accession["accession"].map(lookup)
    project = (
        accession.groupby(["tissue", "method", "project", "term"], as_index=False)[
            "effect"
        ]
        .mean()
        .sort_values(["term", "project"])
    )
    summary = (
        accession.groupby(["tissue", "method", "term"])["effect"]
        .agg(
            accession_balanced_effect="mean",
            n_accessions="size",
            accessions_positive=lambda values: int((values > 0).sum()),
            accessions_negative=lambda values: int((values < 0).sum()),
        )
        .reset_index()
    )
    project_summary = (
        project.groupby(["tissue", "method", "term"])["effect"]
        .agg(
            project_balanced_effect="mean",
            n_projects="size",
            projects_positive=lambda values: int((values > 0).sum()),
            projects_negative=lambda values: int((values < 0).sum()),
        )
        .reset_index()
    )
    summary = summary.merge(
        project_summary, on=["tissue", "method", "term"], how="left"
    )
    return accession, project, summary


def run_ssgsea(data: TissueData) -> pd.DataFrame:
    gene_sets = {
        term: genes
        for term, genes in data.gene_sets.items()
        if MIN_GENE_SET_SIZE <= len(genes) <= MAX_GENE_SET_SIZE
    }
    expression = pd.DataFrame(
        data.log2cpm.T,
        index=data.genes,
        columns=data.obs.index.astype(str),
    )
    result = gp.ssgsea(
        data=expression,
        gene_sets=gene_sets,
        outdir=None,
        sample_norm_method="rank",
        correl_norm_type="rank",
        min_size=MIN_GENE_SET_SIZE,
        max_size=MAX_GENE_SET_SIZE,
        weight=0.25,
        permutation_num=None,
        threads=4,
        no_plot=True,
        seed=RANDOM_SEED,
        verbose=False,
    ).res2d
    score_column = "NES" if "NES" in result else "ES"
    scores = result.pivot(index="Name", columns="Term", values=score_column)
    scores.index = scores.index.astype(str)
    scores = scores.reindex(data.obs.index.astype(str))
    return scores.astype(float)


def _project_gene_effects(data: TissueData) -> pd.DataFrame:
    expression = pd.DataFrame(data.log2cpm, index=data.obs.index, columns=data.genes)
    condition = data.obs["condition_inferred"].astype(str)
    rows = []
    for accession, indexes in data.obs.groupby("id.accession", observed=True).indices.items():
        indexes = np.asarray(indexes)
        local = condition.iloc[indexes]
        flight = local.eq("flight").to_numpy()
        ground = local.eq("ground_control").to_numpy()
        if not flight.any() or not ground.any():
            continue
        matrix = expression.iloc[indexes]
        effect = matrix.iloc[flight].mean(axis=0) - matrix.iloc[ground].mean(axis=0)
        row = effect.to_frame().T
        row.index = [str(accession)]
        rows.append(row)
    accession = pd.concat(rows)
    projects = pd.Series(_project_lookup(data.obs), name="project")
    accession = accession.join(projects, how="left")
    return accession.groupby("project", observed=True).mean(numeric_only=True)


def run_preranked_gsea(data: TissueData) -> pd.DataFrame:
    project_effects = _project_gene_effects(data)
    ranking = project_effects.mean(axis=0).sort_values(ascending=False)
    # Deterministic sub-machine-precision offsets make tied zero effects orderable.
    ranking = ranking + np.linspace(1e-12, -1e-12, len(ranking))
    gene_sets = {
        term: genes
        for term, genes in data.gene_sets.items()
        if MIN_GENE_SET_SIZE <= len(genes) <= MAX_GENE_SET_SIZE
    }
    result = gp.prerank(
        rnk=ranking,
        gene_sets=gene_sets,
        outdir=None,
        min_size=MIN_GENE_SET_SIZE,
        max_size=MAX_GENE_SET_SIZE,
        permutation_num=1000,
        weight=1.0,
        ascending=False,
        threads=4,
        seed=RANDOM_SEED,
        no_plot=True,
        verbose=False,
    ).res2d.copy()
    result = result.rename(
        columns={
            "Term": "term",
            "NES": "gsea_nes",
            "NOM p-val": "gsea_nominal_p",
            "FDR q-val": "gsea_fdr",
            "Lead_genes": "gsea_leading_edge_genes",
        }
    )
    columns = [
        "term",
        "gsea_nes",
        "gsea_nominal_p",
        "gsea_fdr",
        "gsea_leading_edge_genes",
    ]
    result = result[[column for column in columns if column in result]].copy()
    result.insert(0, "tissue", data.tissue)
    for column in ("gsea_nes", "gsea_nominal_p", "gsea_fdr"):
        result[column] = pd.to_numeric(result[column], errors="coerce")
    return result


def project_heldout_folds(
    project_effects: pd.DataFrame,
    curated: pd.DataFrame,
    family: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for (tissue, method), frame in project_effects.groupby(
        ["tissue", "method"], observed=True
    ):
        curated_terms = set(curated.loc[curated["tissue"].eq(tissue), "term"])
        family_terms = set(family.loc[family["tissue"].eq(tissue), "term"])
        wide = frame.pivot(index="project", columns="term", values="effect")
        for heldout_project in wide.index:
            training = wide.drop(index=heldout_project).mean(axis=0)
            heldout = wide.loc[heldout_project]
            percentile = training.abs().rank(pct=True, method="average")
            for term in wide.columns:
                train_value = float(training[term])
                test_value = float(heldout[term])
                rows.append(
                    {
                        "tissue": tissue,
                        "method": method,
                        "heldout_project": str(heldout_project),
                        "n_training_projects": int(len(wide.index) - 1),
                        "term": str(term),
                        "training_project_balanced_effect": train_value,
                        "heldout_project_effect": test_value,
                        "training_absolute_effect_percentile": float(percentile[term]),
                        "selected_in_training_top_decile": bool(percentile[term] >= 0.9),
                        "direction_concordant": bool(
                            np.sign(train_value) == np.sign(test_value)
                            and np.sign(train_value) != 0
                        ),
                        "curated_pathway": term in curated_terms,
                        "family_representative": term in family_terms,
                    }
                )
    return pd.DataFrame(rows)


def summarize_heldout(folds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    masks = {
        "all_active": pd.Series(True, index=folds.index),
        "training_top_decile": folds["selected_in_training_top_decile"],
        "curated": folds["curated_pathway"],
        "family_representatives": folds["family_representative"],
    }
    for name, mask in masks.items():
        subset = folds.loc[mask]
        for keys, frame in subset.groupby(["tissue", "method"], observed=True):
            tissue, method = keys
            per_project = frame.groupby("heldout_project")["direction_concordant"].mean()
            rows.append(
                {
                    "tissue": tissue,
                    "method": method,
                    "pathway_set": name,
                    "n_fold_pathway_predictions": int(len(frame)),
                    "n_heldout_projects": int(frame["heldout_project"].nunique()),
                    "direction_concordance": float(frame["direction_concordant"].mean()),
                    "minimum_project_concordance": float(per_project.min()),
                    "maximum_project_concordance": float(per_project.max()),
                }
            )
    return pd.DataFrame(rows)


def _read_csr_rows(path: Path, group: str, rows: np.ndarray) -> sparse.csr_matrix:
    rows = np.asarray(rows, dtype=int)
    if len(rows) == 0:
        raise ValueError("At least one row is required")
    rows = np.sort(rows)
    with h5py.File(path, "r") as handle:
        node = handle[group]
        shape = tuple(map(int, node.attrs["shape"]))
        indptr = node["indptr"][:]
        boundaries = np.flatnonzero(np.diff(rows) != 1) + 1
        runs = np.split(rows, boundaries)
        chunks = []
        for run in runs:
            start = int(run[0])
            stop = int(run[-1]) + 1
            left = int(indptr[start])
            right = int(indptr[stop])
            local_indptr = indptr[start : stop + 1] - left
            chunks.append(
                sparse.csr_matrix(
                    (
                        node["data"][left:right],
                        node["indices"][left:right],
                        local_indptr,
                    ),
                    shape=(stop - start, shape[1]),
                )
            )
    return sparse.vstack(chunks, format="csr")


def atlas_marker_signatures(
    data: TissueData,
    atlas_obs: pd.DataFrame,
    atlas_genes: pd.Index,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tissue_mask = atlas_obs["tissue"].astype(str).eq(ATLAS_TISSUE[data.tissue])
    atlas_gene_lookup = pd.Series(np.arange(len(atlas_genes)), index=atlas_genes)
    shared = [gene for gene in data.genes if gene in atlas_gene_lookup.index]
    atlas_columns = atlas_gene_lookup.loc[shared].astype(int).to_numpy()
    rows = []
    cell_counts = {}
    for compartment, cell_types in BROAD_COMPARTMENTS[data.tissue].items():
        mask = tissue_mask & atlas_obs["cell_type"].astype(str).isin(cell_types)
        indexes = np.flatnonzero(mask.to_numpy())
        matrix = _read_csr_rows(TMS_PATH, "raw/X", indexes)
        pseudobulk = np.asarray(matrix[:, atlas_columns].sum(axis=0)).ravel()
        cpm = pseudobulk / max(float(pseudobulk.sum()), 1.0) * 1_000_000.0
        rows.append(pd.Series(np.log2(cpm + 1.0), index=shared, name=compartment))
        cell_counts[compartment] = int(len(indexes))
    signatures = pd.DataFrame(rows)

    marker_rows = []
    for compartment in signatures.index:
        others = signatures.drop(index=compartment)
        specificity = signatures.loc[compartment] - others.max(axis=0)
        candidates = specificity.loc[
            signatures.loc[compartment].gt(np.log2(2.0)) & specificity.gt(0)
        ].sort_values(ascending=False)
        for rank, (gene, value) in enumerate(
            candidates.head(N_MARKERS_PER_COMPARTMENT).items(), start=1
        ):
            marker_rows.append(
                {
                    "tissue": data.tissue,
                    "atlas_tissue": ATLAS_TISSUE[data.tissue],
                    "compartment": compartment,
                    "atlas_cells": cell_counts[compartment],
                    "rank": rank,
                    "gene_id": gene,
                    "atlas_log2cpm": float(signatures.loc[compartment, gene]),
                    "specificity_log2_fold": float(value),
                }
            )
    markers = pd.DataFrame(marker_rows)
    return signatures, markers


def composition_marker_scores(data: TissueData, markers: pd.DataFrame) -> pd.DataFrame:
    expression = pd.DataFrame(data.log2cpm, index=data.obs.index, columns=data.genes)
    standardized = pd.DataFrame(
        StandardScaler().fit_transform(expression),
        index=expression.index,
        columns=expression.columns,
    )
    result = pd.DataFrame(index=expression.index)
    for compartment, frame in markers.groupby("compartment", observed=True):
        genes = [gene for gene in frame["gene_id"] if gene in standardized.columns]
        result[str(compartment)] = standardized[genes].mean(axis=1)
    return result


def _within_accession_center(frame: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    centered = frame.copy()
    groups = obs["id.accession"].astype(str)
    return centered - centered.groupby(groups, observed=True).transform("mean")


def composition_adjusted_effects(
    data: TissueData,
    marker_scores: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    condition = data.obs["condition_inferred"].astype(str).eq("flight").astype(float)
    condition = _within_accession_center(condition.to_frame("flight"), data.obs)[
        "flight"
    ]
    centered_markers = _within_accession_center(marker_scores, data.obs)
    scaled = StandardScaler().fit_transform(centered_markers)
    n_components = min(3, scaled.shape[1], scaled.shape[0] - 1)
    pca = PCA(n_components=n_components, random_state=RANDOM_SEED).fit(scaled)
    cumulative = np.cumsum(pca.explained_variance_ratio_)
    retained = int(min(n_components, np.searchsorted(cumulative, 0.9) + 1))
    pcs = pca.transform(scaled)[:, :retained]
    design_unadjusted = condition.to_numpy()[:, None]
    design_adjusted = np.column_stack([condition.to_numpy(), pcs])

    centered_scores = _within_accession_center(data.expimap_scores, data.obs)
    unadjusted = np.linalg.lstsq(
        design_unadjusted, centered_scores.to_numpy(), rcond=None
    )[0][0]
    adjusted = np.linalg.lstsq(
        design_adjusted, centered_scores.to_numpy(), rcond=None
    )[0][0]
    result = pd.DataFrame(
        {
            "tissue": data.tissue,
            "term": centered_scores.columns,
            "within_accession_unadjusted_effect": unadjusted,
            "composition_proxy_adjusted_effect": adjusted,
        }
    )
    result["adjusted_direction_matches_unadjusted"] = (
        np.sign(result["within_accession_unadjusted_effect"])
        == np.sign(result["composition_proxy_adjusted_effect"])
    )
    denominator = result["within_accession_unadjusted_effect"].abs().replace(0, np.nan)
    result["absolute_effect_ratio_adjusted_to_unadjusted"] = (
        result["composition_proxy_adjusted_effect"].abs() / denominator
    )
    result["absolute_attenuation_fraction"] = 1.0 - result[
        "absolute_effect_ratio_adjusted_to_unadjusted"
    ]

    sample = data.obs[
        ["id.accession", "condition_inferred", PROJECT_COLUMN]
    ].copy()
    sample.insert(0, "sample", sample.index.astype(str))
    sample.insert(0, "tissue", data.tissue)
    for column in marker_scores:
        sample[f"marker_score_{column}"] = marker_scores[column].to_numpy()
    for index in range(retained):
        sample[f"within_accession_composition_pc{index + 1}"] = pcs[:, index]
    sample["composition_pcs_retained"] = retained
    sample["composition_pc_variance_explained"] = float(cumulative[retained - 1])
    return result, sample


def add_selection_flags(frame: pd.DataFrame) -> pd.DataFrame:
    curated = pd.concat(
        [_curated_frame(config.tissue).assign(tissue=config.tissue) for config in CONFIGS],
        ignore_index=True,
    )
    families = _family_terms()
    frame = frame.merge(
        curated[["tissue", "term", "short_label", "evidence_role"]],
        on=["tissue", "term"],
        how="left",
    )
    frame = frame.merge(
        families.assign(family_representative=True),
        on=["tissue", "term"],
        how="left",
    )
    frame["curated_pathway"] = frame["evidence_role"].notna()
    frame["family_representative"] = frame["family_representative"].fillna(False)
    return frame


def benchmark_summary(benchmark: pd.DataFrame) -> pd.DataFrame:
    rows = []
    masks = {
        "all_active": pd.Series(True, index=benchmark.index),
        "expimap_top_decile": benchmark.groupby("tissue", observed=True)[
            "expimap_accession_balanced_effect"
        ].transform(lambda values: values.abs().rank(pct=True) >= 0.9),
        "curated": benchmark["curated_pathway"],
        "family_representatives": benchmark["family_representative"],
    }
    for pathway_set, mask in masks.items():
        for tissue, frame in benchmark.loc[mask].groupby("tissue", observed=True):
            complete_ss = frame.dropna(
                subset=[
                    "expimap_accession_balanced_effect",
                    "ssgsea_accession_balanced_effect",
                ]
            )
            complete_gsea = frame.dropna(
                subset=["expimap_accession_balanced_effect", "gsea_nes"]
            )
            rho_ss = stats.spearmanr(
                complete_ss["expimap_accession_balanced_effect"],
                complete_ss["ssgsea_accession_balanced_effect"],
            ).statistic
            rho_gsea = stats.spearmanr(
                complete_gsea["expimap_accession_balanced_effect"],
                complete_gsea["gsea_nes"],
            ).statistic
            rows.append(
                {
                    "tissue": tissue,
                    "pathway_set": pathway_set,
                    "n_pathways": int(len(frame)),
                    "n_ssgsea_comparable": int(len(complete_ss)),
                    "expimap_ssgsea_spearman_rho": float(rho_ss),
                    "expimap_ssgsea_direction_agreement": float(
                        (
                            np.sign(complete_ss["expimap_accession_balanced_effect"])
                            == np.sign(complete_ss["ssgsea_accession_balanced_effect"])
                        ).mean()
                    ),
                    "n_gsea_comparable": int(len(complete_gsea)),
                    "expimap_gsea_spearman_rho": float(rho_gsea),
                    "expimap_gsea_direction_agreement": float(
                        (
                            np.sign(complete_gsea["expimap_accession_balanced_effect"])
                            == np.sign(complete_gsea["gsea_nes"])
                        ).mean()
                    ),
                }
            )
    return pd.DataFrame(rows)


def plot_method_benchmark(
    benchmark: pd.DataFrame,
    summary: pd.DataFrame,
    heldout_summary: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(18, 9), constrained_layout=True)
    for column, config in enumerate(CONFIGS):
        tissue = config.tissue
        frame = benchmark.loc[benchmark["tissue"].eq(tissue)].dropna(
            subset=[
                "expimap_accession_balanced_effect",
                "ssgsea_accession_balanced_effect",
            ]
        )
        ax = axes[0, column]
        ax.scatter(
            frame["expimap_accession_balanced_effect"],
            frame["ssgsea_accession_balanced_effect"],
            s=18,
            color="#b9bec1",
            alpha=0.55,
            linewidth=0,
        )
        selected = frame.loc[frame["curated_pathway"]]
        ax.scatter(
            selected["expimap_accession_balanced_effect"],
            selected["ssgsea_accession_balanced_effect"],
            s=48,
            c=selected["evidence_role"].map(ROLE_COLORS),
            edgecolor="white",
            linewidth=0.7,
        )
        ax.axhline(0, color="#676d70", lw=0.7)
        ax.axvline(0, color="#676d70", lw=0.7)
        row = summary.loc[
            summary["tissue"].eq(tissue) & summary["pathway_set"].eq("all_active")
        ].iloc[0]
        ax.set_title(f"{config.display_name}\nSpearman r = {row.expimap_ssgsea_spearman_rho:.2f}")
        ax.set_xlabel("expiMap study-balanced shift")
        if column == 0:
            ax.set_ylabel("ssGSEA study-balanced shift")

        ax = axes[1, column]
        method_row = summary.loc[
            summary["tissue"].eq(tissue) & summary["pathway_set"].eq("curated")
        ].iloc[0]
        held = heldout_summary.loc[
            heldout_summary["tissue"].eq(tissue)
            & heldout_summary["method"].eq("expimap")
            & heldout_summary["pathway_set"].eq("training_top_decile")
        ].iloc[0]
        values = [
            method_row.expimap_ssgsea_direction_agreement,
            method_row.expimap_gsea_direction_agreement,
            held.direction_concordance,
        ]
        labels = ["ssGSEA\ncurated", "GSEA\ncurated", "Held-out\ntop decile"]
        ax.bar(labels, values, color=["#3c7f8f", "#7a6a9d", TISSUE_COLORS[tissue]])
        ax.axhline(0.5, color="#676d70", linestyle="--", lw=0.8)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Directional agreement" if column == 0 else "")
        ax.tick_params(axis="x", labelsize=8)
        for index, value in enumerate(values):
            ax.text(index, min(value + 0.035, 0.97), f"{value:.0%}", ha="center", fontsize=9)
    fig.suptitle(
        "Method-independent and held-out-project robustness",
        fontsize=16,
        fontweight="bold",
    )
    fig.savefig(FIGURE_DIR / "figure_s4_method_and_heldout_robustness.png", dpi=300)
    fig.savefig(FIGURE_DIR / "figure_s4_method_and_heldout_robustness.pdf")
    plt.close(fig)


def plot_composition_sensitivity(composition: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16, 11), constrained_layout=True)
    for ax, config in zip(axes.flat, CONFIGS):
        frame = composition.loc[
            composition["tissue"].eq(config.tissue) & composition["curated_pathway"]
        ].sort_values("within_accession_unadjusted_effect")
        y = np.arange(len(frame))
        ax.hlines(
            y,
            frame["within_accession_unadjusted_effect"],
            frame["composition_proxy_adjusted_effect"],
            color="#aab0b3",
            lw=1.5,
        )
        ax.scatter(
            frame["within_accession_unadjusted_effect"],
            y,
            marker="o",
            s=45,
            color="#90979a",
            label="Within-accession",
        )
        ax.scatter(
            frame["composition_proxy_adjusted_effect"],
            y,
            marker="D",
            s=48,
            c=frame["evidence_role"].map(ROLE_COLORS),
            edgecolor="white",
            linewidth=0.6,
            label="Composition-proxy adjusted",
        )
        ax.axvline(0, color="#565c5f", lw=0.8)
        ax.set_yticks(y, frame["short_label"], fontsize=9)
        ax.set_xlabel("Flight minus ground coefficient")
        retained = frame["adjusted_direction_matches_unadjusted"].mean()
        ax.set_title(f"{config.display_name} | {retained:.0%} selected directions retained")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    fig.suptitle(
        "Sensitivity to atlas-derived broad cell-composition proxies",
        fontsize=16,
        fontweight="bold",
    )
    fig.savefig(FIGURE_DIR / "figure_s5_composition_proxy_sensitivity.png", dpi=300)
    fig.savefig(FIGURE_DIR / "figure_s5_composition_proxy_sensitivity.pdf")
    plt.close(fig)


def write_report(
    benchmark_summary_frame: pd.DataFrame,
    heldout_summary: pd.DataFrame,
    composition: pd.DataFrame,
    marker_definitions: pd.DataFrame,
) -> None:
    lines = [
        "# Reviewer robustness analysis",
        "",
        "This analysis was specified before running the reviewer-directed robustness outputs. It does not retroactively make the original exploratory pathway review prespecified.",
        "",
        "## Frozen decisions",
        "",
        "- Same query samples, 2,000-HVG gene universe, Reactome memberships, and decoder-oriented primary expiMap scores as the manuscript.",
        "- Conventional benchmarks: rank-normalized ssGSEA per sample and project-balanced log2-CPM preranked GSEA.",
        "- Internal validation: leave-one-project-out direction prediction, with the top decile selected using training projects only.",
        "- Composition sensitivity: broad compartment markers derived from the independent Tabula Muris Senis Smart-seq2 atlas, followed by within-accession regression on marker-score principal components.",
        "- These are triangulation and sensitivity analyses. They are not an external replication cohort, causal adjustment, or cell-type deconvolution.",
        "",
        "## Results",
        "",
    ]
    for config in CONFIGS:
        tissue = config.tissue
        all_row = benchmark_summary_frame.loc[
            benchmark_summary_frame["tissue"].eq(tissue)
            & benchmark_summary_frame["pathway_set"].eq("all_active")
        ].iloc[0]
        curated_row = benchmark_summary_frame.loc[
            benchmark_summary_frame["tissue"].eq(tissue)
            & benchmark_summary_frame["pathway_set"].eq("curated")
        ].iloc[0]
        held_row = heldout_summary.loc[
            heldout_summary["tissue"].eq(tissue)
            & heldout_summary["method"].eq("expimap")
            & heldout_summary["pathway_set"].eq("training_top_decile")
        ].iloc[0]
        comp = composition.loc[
            composition["tissue"].eq(tissue) & composition["curated_pathway"]
        ]
        markers = marker_definitions.loc[marker_definitions["tissue"].eq(tissue)]
        lines.extend(
            [
                f"### {config.display_name}",
                "",
                f"Across active programs, expiMap and ssGSEA had Spearman r={all_row.expimap_ssgsea_spearman_rho:.2f} and {all_row.expimap_ssgsea_direction_agreement:.0%} directional agreement. Among the curated programs, expiMap agreed in direction with ssGSEA for {curated_row.expimap_ssgsea_direction_agreement:.0%} and with preranked GSEA for {curated_row.expimap_gsea_direction_agreement:.0%}.",
                f"Training-only top-decile pathways predicted the held-out project direction in {held_row.direction_concordance:.0%} of fold-pathway comparisons. After adjustment for {markers['compartment'].nunique()} atlas-derived broad compartment scores, {comp['adjusted_direction_matches_unadjusted'].mean():.0%} of curated expiMap directions were retained.",
                "",
            ]
        )
    lines.extend(
        [
            "## Interpretation limits",
            "",
            "Project-wise cross-validation reuses the same repository and cannot replace validation in newly generated missions. The atlas marker analysis is a proxy sensitivity: composition can be a biological mediator of spaceflight, the atlas lacks some mature cell states (notably mature myofibers), and removing marker-associated variation can remove real tissue response as well as composition bias.",
            "",
        ]
    )
    (PAPER_DIR / "reviewer_robustness_analysis.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    if not TMS_PATH.exists():
        raise SystemExit(f"Tabula Muris Senis asset not found: {TMS_PATH}")

    atlas = ad.read_h5ad(TMS_PATH, backed="r")
    atlas_obs = atlas.obs.copy()
    atlas_genes = atlas.var_names.astype(str)
    try:
        atlas.file.close()
    except AttributeError:
        pass

    all_accession = []
    all_project = []
    all_method_summary = []
    all_gsea = []
    all_markers = []
    all_marker_scores = []
    all_composition = []

    family = _family_terms()
    for config in CONFIGS:
        print(f"loading {config.tissue}", flush=True)
        data = load_tissue_data(config)
        curated = _curated_frame(config.tissue)
        family_tissue = family.loc[family["tissue"].eq(config.tissue)]

        expi_accession, expi_project, expi_summary = accession_and_project_effects(
            data.expimap_scores, data.obs, data.tissue, "expimap"
        )
        print(f"ssGSEA {config.tissue}", flush=True)
        ss_scores = run_ssgsea(data)
        ss_accession, ss_project, ss_summary = accession_and_project_effects(
            ss_scores, data.obs, data.tissue, "ssgsea"
        )
        ss_long = ss_scores.copy()
        ss_long.insert(0, "sample", ss_long.index)
        ss_long.insert(0, "tissue", data.tissue)
        ss_long.to_csv(
            SOURCE_DIR / f"ssgsea_sample_scores_{data.tissue}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        print(f"preranked GSEA {config.tissue}", flush=True)
        all_gsea.append(run_preranked_gsea(data))

        print(f"composition proxies {config.tissue}", flush=True)
        _, markers = atlas_marker_signatures(data, atlas_obs, atlas_genes)
        marker_scores = composition_marker_scores(data, markers)
        composition, sample_markers = composition_adjusted_effects(data, marker_scores)
        all_markers.append(markers)
        all_marker_scores.append(sample_markers)
        all_composition.append(composition)

        all_accession.extend([expi_accession, ss_accession])
        all_project.extend([expi_project, ss_project])
        all_method_summary.extend([expi_summary, ss_summary])

    accession = pd.concat(all_accession, ignore_index=True)
    project = pd.concat(all_project, ignore_index=True)
    method_summary = pd.concat(all_method_summary, ignore_index=True)
    gsea = pd.concat(all_gsea, ignore_index=True)
    markers = pd.concat(all_markers, ignore_index=True)
    sample_markers = pd.concat(all_marker_scores, ignore_index=True)
    composition = add_selection_flags(pd.concat(all_composition, ignore_index=True))

    expi = method_summary.loc[method_summary["method"].eq("expimap")].drop(
        columns="method"
    )
    expi = expi.rename(
        columns={
            column: f"expimap_{column}"
            for column in expi.columns
            if column not in {"tissue", "term"}
        }
    )
    ss = method_summary.loc[method_summary["method"].eq("ssgsea")].drop(
        columns="method"
    )
    ss = ss.rename(
        columns={
            column: f"ssgsea_{column}"
            for column in ss.columns
            if column not in {"tissue", "term"}
        }
    )
    benchmark = expi.merge(ss, on=["tissue", "term"], how="left").merge(
        gsea, on=["tissue", "term"], how="left"
    )
    benchmark = add_selection_flags(benchmark)
    benchmark["expimap_ssgsea_direction_match"] = (
        np.sign(benchmark["expimap_accession_balanced_effect"])
        == np.sign(benchmark["ssgsea_accession_balanced_effect"])
    )
    benchmark["expimap_gsea_direction_match"] = (
        np.sign(benchmark["expimap_accession_balanced_effect"])
        == np.sign(benchmark["gsea_nes"])
    )
    benchmark_summary_frame = benchmark_summary(benchmark)

    curated_all = pd.concat(
        [_curated_frame(config.tissue).assign(tissue=config.tissue) for config in CONFIGS],
        ignore_index=True,
    )
    folds = project_heldout_folds(project, curated_all, family)
    heldout_summary = summarize_heldout(folds)

    accession.to_csv(
        SOURCE_DIR / "table_s12_conventional_pathway_accession_effects.tsv",
        sep="\t",
        index=False,
    )
    benchmark.to_csv(
        SOURCE_DIR / "table_s13_method_benchmark.tsv", sep="\t", index=False
    )
    benchmark_summary_frame.to_csv(
        SOURCE_DIR / "table_s14_method_benchmark_summary.tsv", sep="\t", index=False
    )
    folds.to_csv(
        SOURCE_DIR / "table_s15_project_heldout_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    heldout_summary.to_csv(
        SOURCE_DIR / "table_s16_project_heldout_summary.tsv", sep="\t", index=False
    )
    markers.to_csv(
        SOURCE_DIR / "table_s17_tms_compartment_markers.tsv", sep="\t", index=False
    )
    sample_markers.to_csv(
        SOURCE_DIR / "table_s18_sample_composition_proxy_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    composition.to_csv(
        SOURCE_DIR / "table_s19_composition_proxy_adjusted_effects.tsv",
        sep="\t",
        index=False,
    )

    plot_method_benchmark(benchmark, benchmark_summary_frame, heldout_summary)
    plot_composition_sensitivity(composition)
    write_report(benchmark_summary_frame, heldout_summary, composition, markers)
    manifest = {
        "random_seed": RANDOM_SEED,
        "ssgsea": {
            "package": "gseapy",
            "version": gp.__version__,
            "sample_norm_method": "rank",
            "correl_norm_type": "rank",
            "weight": 0.25,
            "minimum_gene_set_size": MIN_GENE_SET_SIZE,
            "maximum_gene_set_size": MAX_GENE_SET_SIZE,
        },
        "preranked_gsea": {
            "ranking": "mean project-balanced log2-CPM flight-minus-ground effect",
            "permutations": 1000,
            "weight": 1.0,
        },
        "heldout": {
            "unit": "project identifier after averaging accession effects within project",
            "selection": "top decile of absolute training-project effect",
        },
        "composition_proxy": {
            "atlas": str(TMS_PATH),
            "atlas_citation": str(atlas.uns.get("citation", "")),
            "markers_per_compartment": N_MARKERS_PER_COMPARTMENT,
            "adjustment": "within-accession condition coefficient plus up to three marker-score PCs explaining at least 90% of marker-score variance",
        },
    }
    (SOURCE_DIR / "reviewer_robustness_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
