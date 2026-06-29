"""Targeted skeletal-muscle module analysis for OSDR FLT/GC profiles."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request

from .cluster_enrichment import bh_fdr
from .io import require_import
from .prepare_expimap_osdr_tissue import load_counts_from_api_tables
from .validate_expimap_accession_effects import random_effects


DEFAULT_METADATA = "outputs/expimap_direct_osdr_skeletal_muscle/input/profile_metadata.tsv"
DEFAULT_COUNTS_DIR = "data/osdr_api/counts"
DEFAULT_PATHWAY_SCORES = (
    "outputs/expimap_archs4_reference_osdr_query_skeletal_muscle/"
    "query_nb_allref_50epoch/query_pathway_scores.tsv"
)
DEFAULT_OUTPUT_DIR = "outputs/osdr_skeletal_muscle_targeted_modules"
DEFAULT_REACTOME_GMT = "data/pathways/reactome_current_mouse_ensembl.gmt"
DEFAULT_TARGETED_GMT = "data/pathways/mouse_muscle_targeted_modules.gmt"
DEFAULT_COMBINED_GMT = (
    "data/pathways/reactome_current_mouse_ensembl_plus_muscle_targeted.gmt"
)


MODULES: dict[str, list[str]] = {
    "atrophy_ubiquitin_foxo": [
        "Fbxo32",
        "Trim63",
        "Foxo1",
        "Foxo3",
        "Mstn",
        "Ctsl",
        "Capn1",
        "Capn2",
        "Ubb",
        "Ubc",
        "Ube2b",
        "Psma1",
        "Psma2",
        "Psma3",
        "Psmb5",
    ],
    "autophagy_lysosome": [
        "Map1lc3b",
        "Becn1",
        "Atg5",
        "Atg7",
        "Sqstm1",
        "Lamp1",
        "Lamp2",
        "Ctsb",
        "Ctsd",
        "Tfeb",
        "Ulk1",
    ],
    "igf_akt_mtor": [
        "Igf1",
        "Igf1r",
        "Irs1",
        "Akt1",
        "Akt2",
        "Pik3r1",
        "Pten",
        "Tsc1",
        "Tsc2",
        "Mtor",
        "Rptor",
        "Rictor",
        "Rps6kb1",
        "Eif4ebp1",
    ],
    "calcium_excitation_contraction": [
        "Atp2a1",
        "Atp2a2",
        "Ryr1",
        "Casq1",
        "Casq2",
        "Cacna1s",
        "Cacng1",
        "Sln",
        "Pln",
        "Camk2d",
        "Pvalb",
    ],
    "contractile_fast_sarcomere": [
        "Acta1",
        "Ckm",
        "Myh1",
        "Myh2",
        "Myh4",
        "Mybpc1",
        "Mybpc2",
        "Tnni2",
        "Tnnt3",
        "Tpm1",
    ],
    "contractile_slow_postural": [
        "Myh7",
        "Myh7b",
        "Mybpc1",
        "Tnni1",
        "Tnnt1",
        "Tpm3",
        "Actn2",
        "Des",
        "Myl2",
        "Myl3",
    ],
    "mitochondrial_oxphos_tca": [
        "Ppargc1a",
        "Tfam",
        "Cs",
        "Aco2",
        "Idh3a",
        "Sdha",
        "Sdhb",
        "Uqcrc1",
        "Cox4i1",
        "Cox5a",
        "Ndufs1",
        "Ndufb8",
        "Atp5f1a",
        "Atp5f1b",
        "Cycs",
    ],
    "fatty_acid_oxidation": [
        "Cpt1b",
        "Slc25a20",
        "Acadl",
        "Acadm",
        "Acads",
        "Hadha",
        "Hadhb",
        "Echs1",
        "Etfa",
        "Ppara",
        "Ppard",
    ],
    "ecm_fibrosis_remodeling": [
        "Col1a1",
        "Col1a2",
        "Col3a1",
        "Col4a1",
        "Col5a1",
        "Fn1",
        "Tnc",
        "Postn",
        "Mmp2",
        "Mmp9",
        "Tgfb1",
        "Ctgf",
    ],
    "inflammation_oxidative_stress": [
        "Il6",
        "Tnf",
        "Nfkb1",
        "Rela",
        "Ccl2",
        "Cxcl10",
        "Stat1",
        "Irf7",
        "Ddit3",
        "Hspa1a",
        "Hsp90aa1",
        "Nfe2l2",
        "Sod2",
        "Gpx1",
        "Cat",
    ],
    "myogenesis_regeneration": [
        "Pax7",
        "Myod1",
        "Myog",
        "Myf5",
        "Myf6",
        "Myh3",
        "Myh8",
        "Cdkn1a",
    ],
    "neuromuscular_junction": [
        "Chrna1",
        "Chrnb1",
        "Chrnd",
        "Chrne",
        "Chrng",
        "Musk",
        "Lrp4",
        "Agrn",
        "Dok7",
    ],
}


def muscle_group(value: object) -> str:
    text = str(value or "").strip().lower()
    if "soleus" in text:
        return "soleus"
    if "extensor digitorum" in text or text == "edl":
        return "edl"
    if "tibialis" in text:
        return "tibialis_anterior"
    if "gastrocnemius" in text:
        return "gastrocnemius"
    if "quadriceps" in text:
        return "quadriceps"
    cleaned = (
        text.replace("-", "_")
        .replace("/", "_")
        .replace(" ", "_")
        .replace("__", "_")
        .strip("_")
    )
    return cleaned or "unknown"


def read_metadata(path: str | Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    metadata = pd.read_csv(path, sep="\t", keep_default_na=False)
    if "profile_id" not in metadata:
        metadata["profile_id"] = (
            metadata["id.accession"].astype(str)
            + "/"
            + metadata["id.sample name"].astype(str)
        )
    material_col = "study.characteristics.material type"
    if material_col not in metadata:
        raise SystemExit(f"Metadata missing `{material_col}`.")
    metadata["muscle_group"] = metadata[material_col].map(muscle_group)
    return metadata


def count_paths_for_metadata(metadata, counts_dir: Path) -> dict[str, Path]:
    return {
        accession: counts_dir / f"{accession}_unnormalized_counts.csv"
        for accession in sorted(metadata["id.accession"].astype(str).unique())
    }


def query_ensembl_symbol(symbol: str, timeout: int) -> dict[str, object]:
    quoted = urllib.parse.quote(symbol)
    url = (
        "https://rest.ensembl.org/xrefs/symbol/mus_musculus/"
        f"{quoted}?content-type=application/json"
    )
    with urllib.request.urlopen(url, timeout=timeout) as response:
        data = json.loads(response.read().decode("utf-8"))
    candidates = [
        str(row.get("id", ""))
        for row in data
        if row.get("type") == "gene" and str(row.get("id", "")).startswith("ENSMUSG")
    ]
    return {
        "gene_symbol": symbol,
        "candidate_ensembl_genes": ";".join(dict.fromkeys(candidates)),
        "n_candidates": len(dict.fromkeys(candidates)),
    }


def load_or_query_symbol_map(
    symbols: list[str],
    count_genes: set[str],
    cache_path: Path,
    timeout: int,
):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        cache = pd.read_csv(cache_path, sep="\t", keep_default_na=False)
    else:
        cache = pd.DataFrame(
            columns=[
                "gene_symbol",
                "candidate_ensembl_genes",
                "n_candidates",
                "selected_ensembl_gene",
                "present_in_counts",
            ]
        )

    seen = set(cache["gene_symbol"].astype(str)) if not cache.empty else set()
    rows = []
    for symbol in symbols:
        if symbol in seen:
            continue
        try:
            result = query_ensembl_symbol(symbol, timeout)
        except Exception as exc:  # pragma: no cover - network-dependent path
            result = {
                "gene_symbol": symbol,
                "candidate_ensembl_genes": "",
                "n_candidates": 0,
                "mapping_error": str(exc),
            }
        candidates = [
            candidate
            for candidate in str(result.get("candidate_ensembl_genes", "")).split(";")
            if candidate
        ]
        present = [candidate for candidate in candidates if candidate in count_genes]
        selected = present[0] if present else (candidates[0] if candidates else "")
        result["selected_ensembl_gene"] = selected
        result["present_in_counts"] = bool(selected in count_genes)
        rows.append(result)
        time.sleep(0.05)

    if rows:
        cache = pd.concat([cache, pd.DataFrame(rows)], ignore_index=True)
        cache = cache.drop_duplicates("gene_symbol", keep="last").sort_values("gene_symbol")
        cache.to_csv(cache_path, sep="\t", index=False)

    return cache


def module_member_table(symbol_map, count_genes: set[str]):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    lookup = symbol_map.set_index("gene_symbol").to_dict(orient="index")
    rows = []
    for module, symbols in MODULES.items():
        for symbol in symbols:
            row = lookup.get(symbol, {})
            ensembl = str(row.get("selected_ensembl_gene", ""))
            rows.append(
                {
                    "module": module,
                    "gene_symbol": symbol,
                    "ensembl_gene": ensembl,
                    "present_in_counts": bool(ensembl in count_genes),
                    "candidate_ensembl_genes": row.get("candidate_ensembl_genes", ""),
                }
            )
    return pd.DataFrame(rows)


def write_targeted_gmts(members, reactome_gmt: Path, targeted_gmt: Path, combined_gmt: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    targeted_gmt.parent.mkdir(parents=True, exist_ok=True)
    combined_gmt.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    term_rows = []
    for module, group in members.groupby("module", sort=False):
        genes = list(
            dict.fromkeys(
                gene
                for gene in group.loc[group["present_in_counts"], "ensembl_gene"].astype(str)
                if gene
            )
        )
        term = f"MUSCLE_{module.upper()}"
        description = f"targeted_muscle_module:{module}"
        if genes:
            rows.append("\t".join([term, description, *genes]))
        term_rows.append(
            {
                "term": term,
                "description": description,
                "module": module,
                "n_genes": len(genes),
                "genes": ";".join(genes),
            }
        )

    targeted_gmt.write_text("\n".join(rows) + "\n", encoding="utf-8")
    terms = pd.DataFrame(term_rows)
    terms.to_csv(targeted_gmt.with_suffix(".terms.tsv"), sep="\t", index=False)

    reactome_text = reactome_gmt.read_text(encoding="utf-8").rstrip("\n")
    combined_text = reactome_text + "\n" + "\n".join(rows) + "\n"
    combined_gmt.write_text(combined_text, encoding="utf-8")
    return {
        "targeted_gmt": str(targeted_gmt),
        "targeted_terms": str(targeted_gmt.with_suffix(".terms.tsv")),
        "combined_gmt": str(combined_gmt),
        "n_targeted_terms": len(rows),
        "n_combined_terms": int(
            sum(1 for line in combined_text.splitlines() if line.strip())
        ),
    }


def log1p_cpm(counts):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    library = counts.sum(axis=0).astype("float64")
    library[library <= 0] = 1.0
    cpm = counts.div(library, axis=1) * 1_000_000.0
    return np.log1p(cpm)


def zscore_rows(frame):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    means = frame.mean(axis=1)
    stds = frame.std(axis=1, ddof=0).replace(0, np.nan)
    return frame.sub(means, axis=0).div(stds, axis=0).fillna(0.0)


def score_modules(metadata, counts, members):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    transformed = log1p_cpm(counts)
    z = zscore_rows(transformed)

    score_frame = metadata.copy()
    score_frame = score_frame.set_index("profile_id", drop=False)
    z = z.loc[:, score_frame.index.astype(str).tolist()]
    for module, group in members.groupby("module", sort=False):
        genes = [
            gene
            for gene in group.loc[group["present_in_counts"], "ensembl_gene"].astype(str)
            if gene in z.index
        ]
        if genes:
            score_frame[module] = z.loc[genes].mean(axis=0).to_numpy()
        else:
            score_frame[module] = float("nan")
    return score_frame.reset_index(drop=True), z


def variance_of_difference(flight, ground) -> float:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    x = np.asarray(flight, dtype=float)
    y = np.asarray(ground, dtype=float)
    x_var = float(np.var(x, ddof=1) / len(x)) if len(x) > 1 else 0.0
    y_var = float(np.var(y, ddof=1) / len(y)) if len(y) > 1 else 0.0
    variance = x_var + y_var
    return variance if variance > 0 else float("nan")


def effect_table(frame, terms: list[str], group_label: str):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    rows = []
    strata = [("all_skeletal_muscle", frame)]
    strata.extend((name, group) for name, group in frame.groupby("muscle_group", sort=True))
    for stratum, stratum_frame in strata:
        for accession, accession_frame in stratum_frame.groupby("id.accession", sort=True):
            flight_mask = accession_frame["condition_inferred"].eq("flight")
            ground_mask = accession_frame["condition_inferred"].eq("ground_control")
            if not flight_mask.any() or not ground_mask.any():
                continue
            for term in terms:
                flight = accession_frame.loc[flight_mask, term].astype(float)
                ground = accession_frame.loc[ground_mask, term].astype(float)
                rows.append(
                    {
                        "analysis_group": stratum,
                        "id.accession": accession,
                        "term": term,
                        "score_type": group_label,
                        "n_flight": int(len(flight)),
                        "n_ground_control": int(len(ground)),
                        "mean_flight": float(flight.mean()),
                        "mean_ground_control": float(ground.mean()),
                        "flight_minus_ground": float(flight.mean() - ground.mean()),
                        "effect_variance": variance_of_difference(flight, ground),
                    }
                )
    return pd.DataFrame(rows)


def meta_table(effects):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    for (analysis_group, term, score_type), group in effects.groupby(
        ["analysis_group", "term", "score_type"], sort=False
    ):
        result = random_effects(
            group["flight_minus_ground"].to_numpy(),
            group["effect_variance"].to_numpy(),
            scipy_stats,
        )
        direction = np.sign(result["meta_effect"])
        individual = np.sign(group["flight_minus_ground"].astype(float).to_numpy())
        rows.append(
            {
                "analysis_group": analysis_group,
                "term": term,
                "score_type": score_type,
                **result,
                "n_accession_same_direction": int((individual == direction).sum()),
                "n_accession_opposite_direction": int((individual == -direction).sum()),
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["meta_fdr"] = result.groupby(["analysis_group", "score_type"], sort=False)[
        "meta_p"
    ].transform(lambda values: bh_fdr(values.fillna(1.0).to_numpy()))
    result["direction_consistency"] = (
        result["n_accession_same_direction"] / result["n_accessions"]
    )
    result["strict_candidate"] = (
        result["meta_fdr"].lt(0.05)
        & result["n_accessions"].ge(3)
        & result["n_accession_same_direction"].eq(result["n_accessions"])
    )
    return result.sort_values(
        ["score_type", "analysis_group", "meta_fdr", "meta_p"],
        kind="stable",
    )


def leave_one_out_tables(effects, primary_meta):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    for (analysis_group, score_type), scope in effects.groupby(
        ["analysis_group", "score_type"], sort=False
    ):
        accessions = sorted(scope["id.accession"].astype(str).unique())
        if len(accessions) < 3:
            continue
        for held_out in accessions:
            retained = scope.loc[scope["id.accession"].astype(str).ne(held_out)]
            for term, group in retained.groupby("term", sort=False):
                result = random_effects(
                    group["flight_minus_ground"].to_numpy(),
                    group["effect_variance"].to_numpy(),
                    scipy_stats,
                )
                rows.append(
                    {
                        "analysis_group": analysis_group,
                        "score_type": score_type,
                        "held_out_accession": held_out,
                        "term": term,
                        **result,
                    }
                )
    loo = pd.DataFrame(rows)
    if loo.empty:
        return loo, pd.DataFrame()
    loo["meta_fdr"] = loo.groupby(
        ["analysis_group", "score_type", "held_out_accession"], sort=False
    )["meta_p"].transform(lambda values: bh_fdr(values.fillna(1.0).to_numpy()))

    primary = primary_meta.set_index(["analysis_group", "score_type", "term"])[
        "meta_effect"
    ].to_dict()
    summary_rows = []
    for (analysis_group, score_type, term), group in loo.groupby(
        ["analysis_group", "score_type", "term"], sort=False
    ):
        primary_direction = np.sign(primary.get((analysis_group, score_type, term), 0.0))
        loo_direction = np.sign(group["meta_effect"].astype(float).to_numpy())
        summary_rows.append(
            {
                "analysis_group": analysis_group,
                "score_type": score_type,
                "term": term,
                "n_leave_one_out": int(len(group)),
                "n_loo_same_direction": int((loo_direction == primary_direction).sum()),
                "minimum_leave_one_out_fdr": float(group["meta_fdr"].min()),
                "maximum_leave_one_out_fdr": float(group["meta_fdr"].max()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary["loo_all_same_direction"] = summary["n_loo_same_direction"].eq(
        summary["n_leave_one_out"]
    )
    summary["loo_fdr_stable"] = summary["maximum_leave_one_out_fdr"].lt(0.05)
    return loo, summary


def gene_score_frame(metadata, z, members):
    present = members.loc[members["present_in_counts"]].copy()
    present = present.drop_duplicates("ensembl_gene", keep="first")
    genes = [gene for gene in present["ensembl_gene"].astype(str) if gene in z.index]
    frame = metadata.copy().set_index("profile_id", drop=False)
    values = z.loc[genes, frame.index.astype(str).tolist()].T
    values.columns = genes
    return frame.join(values).reset_index(drop=True), present


def pathway_terms(frame) -> list[str]:
    return [column for column in frame.columns if str(column).startswith("R-MMU-")]


def add_muscle_group_to_scores(scores):
    material_col = "study.characteristics.material type"
    if material_col not in scores:
        raise SystemExit(f"Pathway scores missing `{material_col}`.")
    scores = scores.copy()
    scores["muscle_group"] = scores[material_col].map(muscle_group)
    return scores


def sample_counts(metadata):
    counts = (
        metadata.groupby(
            [
                "muscle_group",
                "id.accession",
                "condition_inferred",
                "study.characteristics.material type",
            ],
            dropna=False,
        )
        .size()
        .rename("n_samples")
        .reset_index()
        .sort_values(["muscle_group", "id.accession", "condition_inferred"])
    )
    return counts


def plot_module_heatmap(meta, output_dir: Path):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")

    frame = meta.loc[meta["score_type"].eq("targeted_module")].copy()
    if frame.empty:
        return None
    pivot = frame.pivot(index="term", columns="analysis_group", values="meta_effect")
    preferred_cols = [
        "all_skeletal_muscle",
        "soleus",
        "edl",
        "tibialis_anterior",
        "gastrocnemius",
        "quadriceps",
    ]
    columns = [column for column in preferred_cols if column in pivot.columns]
    columns += [column for column in pivot.columns if column not in columns]
    pivot = pivot[columns]
    values = pivot.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(values)) if np.isfinite(values).any() else 1.0
    vmax = max(float(vmax), 0.01)

    fig, ax = plt.subplots(figsize=(max(7, len(columns) * 1.2), 6))
    image = ax.imshow(values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels(columns, rotation=35, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Targeted module FLT-GC random-effects meta effect")
    ax.set_xlabel("muscle group")
    ax.set_ylabel("module")
    fig.colorbar(image, ax=ax, shrink=0.82, label="FLT-GC effect (z-score units)")
    fig.tight_layout()
    path = output_dir / "targeted_module_effect_heatmap.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_module_boxplots(scores, meta, output_dir: Path, top_n: int = 6):
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")

    frame = meta.loc[meta["score_type"].eq("targeted_module")].copy()
    frame = frame.loc[frame["analysis_group"].ne("all_skeletal_muscle")]
    if frame.empty:
        return None
    selected = (
        frame.sort_values(["meta_fdr", "meta_p"])
        .drop_duplicates(["analysis_group", "term"], keep="first")
        .head(top_n)
    )
    if selected.empty:
        return None
    fig, axes = plt.subplots(
        len(selected),
        1,
        figsize=(8, max(3.0, len(selected) * 2.2)),
        squeeze=False,
    )
    colors = {"flight": "#c43c39", "ground_control": "#2878b5"}
    for ax, (_, row) in zip(axes[:, 0], selected.iterrows()):
        subset = scores.loc[scores["muscle_group"].eq(row["analysis_group"])]
        positions = []
        labels = []
        data = []
        color_values = []
        for index, condition in enumerate(["ground_control", "flight"], start=1):
            values = subset.loc[subset["condition_inferred"].eq(condition), row["term"]]
            if values.empty:
                continue
            data.append(values.astype(float).to_numpy())
            positions.append(index)
            labels.append(condition)
            color_values.append(colors[condition])
        ax.boxplot(data, positions=positions, widths=0.5, patch_artist=False)
        for pos, values, color in zip(positions, data, color_values):
            jitter = [pos] * len(values)
            ax.scatter(jitter, values, s=18, alpha=0.75, color=color)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
        ax.set_ylabel(row["term"])
        ax.set_title(
            f"{row['analysis_group']} | FDR={row['meta_fdr']:.3g}, "
            f"effect={row['meta_effect']:.3g}"
        )
    fig.tight_layout()
    path = output_dir / "top_targeted_module_boxplots.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def markdown_table(frame, columns: list[str], max_rows: int = 10) -> list[str]:
    if frame.empty:
        return ["No rows."]
    rows = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for _, row in frame.head(max_rows).iterrows():
        values = []
        for column in columns:
            value = row.get(column, "")
            if isinstance(value, float):
                values.append(f"{value:.4g}")
            else:
                values.append(str(value))
        rows.append("| " + " | ".join(values) + " |")
    return rows


def write_readme(output_dir: Path, summary: dict, module_meta, gene_meta, reactome_meta):
    if reactome_meta is None:
        reactome_lines = ["Reactome pathway scores were not supplied."]
    else:
        reactome_display = reactome_meta.loc[reactome_meta["strict_candidate"]].copy()
        if reactome_display.empty:
            reactome_display = reactome_meta.copy()
        reactome_lines = markdown_table(
            reactome_display.sort_values(["meta_fdr", "meta_p"])
            if not reactome_display.empty
            else reactome_display,
            [
                "analysis_group",
                "term",
                "n_accessions",
                "meta_effect",
                "meta_p",
                "meta_fdr",
                "i2",
                "n_accession_same_direction",
                "direction_consistency",
                "strict_candidate",
            ],
            max_rows=20,
        )

    lines = [
        "# OSDR Skeletal-Muscle Targeted Module Analysis",
        "",
        "This analysis splits the OSDR skeletal-muscle FLT/GC profiles by muscle "
        "material type and scores targeted modules from the full downloaded OSDR "
        "unnormalized count tables. Scores are mean z-scored log1p(CPM) expression "
        "over mapped module genes.",
        "",
        f"- Samples: {summary['n_samples']}",
        f"- Accessions: {summary['n_accessions']}",
        f"- Count genes before module filtering: {summary['n_count_genes']}",
        f"- Targeted modules: {summary['n_modules']}",
        f"- Module gene entries present in counts: {summary['n_module_gene_entries_present']}",
        f"- Unique module genes present in counts: {summary['n_unique_module_genes_present']}",
        f"- Targeted GMT: `{summary['gmt_outputs']['targeted_gmt']}`",
        f"- Combined Reactome+targeted GMT: `{summary['gmt_outputs']['combined_gmt']}`",
        "",
        "## Sample Counts",
        "",
        *markdown_table(
            summary["sample_counts_preview"],
            [
                "muscle_group",
                "id.accession",
                "condition_inferred",
                "study.characteristics.material type",
                "n_samples",
            ],
            max_rows=20,
        ),
        "",
        "## Caution",
        "",
        "The all-skeletal-muscle result has 13 accessions. Per-muscle groups are "
        "smaller: EDL and tibialis anterior have two accessions, soleus and "
        "gastrocnemius have three, and quadriceps has four. Use the per-muscle "
        "results as targeted follow-up candidates, not final biology calls. "
        "`strict_candidate` means FDR < 0.05, at least three accessions, and all "
        "accession effects in the same direction. `strict_loo_candidate` also "
        "requires every leave-one-accession-out fit to keep the same direction "
        "and FDR < 0.05.",
        "",
        "## Targeted Module Meta-Analysis",
        "",
        *markdown_table(
            module_meta.sort_values(["meta_fdr", "meta_p"]),
            [
                "analysis_group",
                "term",
                "n_accessions",
                "meta_effect",
                "meta_p",
                "meta_fdr",
                "i2",
                "n_accession_same_direction",
                "direction_consistency",
                "strict_candidate",
                "maximum_leave_one_out_fdr",
                "loo_all_same_direction",
                "strict_loo_candidate",
            ],
            max_rows=20,
        ),
        "",
        "## Targeted Gene Meta-Analysis",
        "",
        *markdown_table(
            gene_meta.sort_values(["meta_fdr", "meta_p"]),
            [
                "analysis_group",
                "gene_symbol",
                "term",
                "n_accessions",
                "meta_effect",
                "meta_p",
                "meta_fdr",
                "i2",
                "direction_consistency",
                "strict_candidate",
            ],
            max_rows=20,
        ),
        "",
        "## Reactome expiMap By Muscle Group",
        "",
        "This table shows strict candidates first when any are present. Reactome "
        "stratification is exploratory because it repeats 1,140 pathway tests "
        "inside small muscle-type subsets.",
        "",
        *reactome_lines,
        "",
        "## Outputs",
        "",
        "- `module_scores.tsv`: sample-level targeted module scores.",
        "- `module_random_effects_meta.tsv`: accession-aware module meta-analysis.",
        "- `module_leave_one_out_summary.tsv`: module leave-one-accession-out stability.",
        "- `gene_random_effects_meta.tsv`: accession-aware module-gene meta-analysis.",
        "- `reactome_by_muscle_group_random_effects_meta.tsv`: stratified expiMap Reactome pathway meta-analysis, if pathway scores were supplied.",
        "- `plots/targeted_module_effect_heatmap.png`: targeted module effects by muscle group.",
        "- `plots/top_targeted_module_boxplots.png`: sample-level scores for top module/group candidates.",
        f"- `{summary['gmt_outputs']['targeted_gmt']}`: targeted muscle module GMT.",
        f"- `{summary['gmt_outputs']['combined_gmt']}`: Reactome GMT with targeted muscle modules appended.",
    ]
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    output_dir = Path(args.output_dir)
    plot_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    metadata = read_metadata(args.metadata)
    counts = load_counts_from_api_tables(
        metadata,
        count_paths_for_metadata(metadata, Path(args.counts_dir)),
    )
    count_genes = set(counts.index.astype(str))
    symbols = sorted({symbol for genes in MODULES.values() for symbol in genes})
    symbol_map = load_or_query_symbol_map(
        symbols,
        count_genes,
        output_dir / "gene_symbol_to_ensembl.tsv",
        args.timeout,
    )
    members = module_member_table(symbol_map, count_genes)
    members.to_csv(output_dir / "module_definitions.tsv", sep="\t", index=False)
    gmt_outputs = write_targeted_gmts(
        members,
        Path(args.reactome_gmt),
        Path(args.targeted_gmt_output),
        Path(args.combined_gmt_output),
    )

    counts_by_group = sample_counts(metadata)
    counts_by_group.to_csv(output_dir / "sample_counts_by_muscle_group.tsv", sep="\t", index=False)

    module_scores, z = score_modules(metadata, counts, members)
    module_terms = list(MODULES)
    module_scores.to_csv(output_dir / "module_scores.tsv", sep="\t", index=False)

    module_effects = effect_table(module_scores, module_terms, "targeted_module")
    module_meta = meta_table(module_effects)
    module_loo, module_loo_summary = leave_one_out_tables(module_effects, module_meta)
    if not module_loo_summary.empty:
        module_meta = module_meta.merge(
            module_loo_summary,
            on=["analysis_group", "score_type", "term"],
            how="left",
        )
    module_meta["strict_loo_candidate"] = (
        module_meta["strict_candidate"]
        & module_meta["loo_all_same_direction"].fillna(False)
        & module_meta["loo_fdr_stable"].fillna(False)
    )
    module_effects.to_csv(output_dir / "module_per_accession_effects.tsv", sep="\t", index=False)
    module_meta.to_csv(output_dir / "module_random_effects_meta.tsv", sep="\t", index=False)
    module_loo.to_csv(output_dir / "module_leave_one_accession_out.tsv", sep="\t", index=False)
    module_loo_summary.to_csv(
        output_dir / "module_leave_one_out_summary.tsv",
        sep="\t",
        index=False,
    )

    gene_scores, present_genes = gene_score_frame(metadata, z, members)
    gene_terms = present_genes["ensembl_gene"].astype(str).tolist()
    gene_effects = effect_table(gene_scores, gene_terms, "targeted_gene")
    gene_meta = meta_table(gene_effects)
    gene_annotations = present_genes[
        ["module", "gene_symbol", "ensembl_gene"]
    ].rename(columns={"ensembl_gene": "term"})
    gene_meta = gene_meta.merge(gene_annotations, on="term", how="left")
    gene_effects = gene_effects.merge(gene_annotations, on="term", how="left")
    gene_effects.to_csv(output_dir / "gene_per_accession_effects.tsv", sep="\t", index=False)
    gene_meta.to_csv(output_dir / "gene_random_effects_meta.tsv", sep="\t", index=False)

    reactome_meta = None
    pathway_scores_path = Path(args.pathway_scores) if args.pathway_scores else None
    if pathway_scores_path and pathway_scores_path.exists():
        pathway_scores = add_muscle_group_to_scores(
            pd.read_csv(pathway_scores_path, sep="\t", keep_default_na=False)
        )
        terms = pathway_terms(pathway_scores)
        reactome_effects = effect_table(pathway_scores, terms, "expimap_reactome")
        reactome_meta = meta_table(reactome_effects)
        reactome_effects.to_csv(
            output_dir / "reactome_by_muscle_group_per_accession_effects.tsv",
            sep="\t",
            index=False,
        )
        reactome_meta.to_csv(
            output_dir / "reactome_by_muscle_group_random_effects_meta.tsv",
            sep="\t",
            index=False,
        )

    heatmap = plot_module_heatmap(module_meta, plot_dir)
    boxplots = plot_module_boxplots(module_scores, module_meta, plot_dir)

    summary = {
        "metadata": str(args.metadata),
        "counts_dir": str(args.counts_dir),
        "pathway_scores": str(args.pathway_scores) if args.pathway_scores else None,
        "n_samples": int(metadata.shape[0]),
        "n_accessions": int(metadata["id.accession"].nunique()),
        "n_count_genes": int(counts.shape[0]),
        "n_modules": int(len(MODULES)),
        "n_module_gene_symbols": int(len(symbols)),
        "n_module_gene_entries_present": int(members["present_in_counts"].sum()),
        "n_unique_module_genes_present": int(
            members.loc[members["present_in_counts"], "ensembl_gene"].nunique()
        ),
        "module_meta_fdr_lt_005": int((module_meta["meta_fdr"] < 0.05).sum()),
        "module_strict_candidates": int(module_meta["strict_candidate"].sum()),
        "module_strict_loo_candidates": int(module_meta["strict_loo_candidate"].sum()),
        "gene_meta_fdr_lt_005": int((gene_meta["meta_fdr"] < 0.05).sum()),
        "gene_strict_candidates": int(gene_meta["strict_candidate"].sum()),
        "reactome_meta_fdr_lt_005": int((reactome_meta["meta_fdr"] < 0.05).sum())
        if reactome_meta is not None and not reactome_meta.empty
        else None,
        "reactome_strict_candidates": int(reactome_meta["strict_candidate"].sum())
        if reactome_meta is not None and not reactome_meta.empty
        else None,
        "gmt_outputs": gmt_outputs,
        "sample_counts_preview": counts_by_group.head(40),
        "plots": {
            "targeted_module_effect_heatmap": str(heatmap) if heatmap else None,
            "top_targeted_module_boxplots": str(boxplots) if boxplots else None,
        },
    }
    serializable_summary = {
        key: value
        for key, value in summary.items()
        if key != "sample_counts_preview"
    }
    (output_dir / "analysis_summary.json").write_text(
        json.dumps(serializable_summary, indent=2) + "\n",
        encoding="utf-8",
    )
    write_readme(output_dir, summary, module_meta, gene_meta, reactome_meta)
    print(json.dumps(serializable_summary, indent=2))
    return output_dir / "analysis_summary.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Score targeted skeletal-muscle modules and run accession-aware "
            "FLT-vs-GC tests globally and by muscle material type."
        )
    )
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--counts-dir", default=DEFAULT_COUNTS_DIR)
    parser.add_argument("--pathway-scores", default=DEFAULT_PATHWAY_SCORES)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--reactome-gmt", default=DEFAULT_REACTOME_GMT)
    parser.add_argument("--targeted-gmt-output", default=DEFAULT_TARGETED_GMT)
    parser.add_argument("--combined-gmt-output", default=DEFAULT_COMBINED_GMT)
    parser.add_argument("--timeout", type=int, default=30)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
