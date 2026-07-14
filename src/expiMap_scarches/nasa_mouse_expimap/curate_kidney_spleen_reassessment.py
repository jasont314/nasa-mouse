"""Manually curate the largest kidney and spleen reassessment pathways."""

from __future__ import annotations

import anndata as ad
import h5py
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats

from .analyze_kidney_spleen_reassessment import load_data
from .build_asgsr_paper import ROOT
from .run_kidney_spleen_seed_sensitivity import CONFIGS, OUTPUT_DIR, SEEDS


# Each entry was reviewed as a complete Reactome label in tissue context. The
# mapping is intentionally explicit so that every top-decile decision is auditable.
MANUAL_REVIEW = {
    "kidney": {
        "R-MMU-114608_PLATELET_DEGRANULATION": ("vascular_hemostatic", "plausible_context", "supporting", "Repeatable vascular or blood-cell response, but not evidence of kidney-cell platelet activation.", "Siew2024"),
        "R-MMU-9711123_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS": ("stress_response", "sensitivity_dependent", "exclude", "Lower primary score is attenuated by composition adjustment and lacks complete method support.", "Siew2024"),
        "R-MMU-1430728_METABOLISM": ("metabolic", "broad_redundant", "supporting", "A broad parent term that is directionally heterogeneous with its amino-acid, oxidation, and lipid children.", "Finch2025;Siew2024"),
        "R-MMU-381426_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS": ("matrix_growth_repair", "coherent_complementary", "main_candidate", "Higher IGF transport connects growth-factor availability to the matrix-repair response; this is a pathway-level extension rather than a demonstrated increase in IGF protein activity.", "Finch2025;Siew2024"),
        "R-MMU-6798695_NEUTROPHIL_DEGRANULATION": ("immune", "seed_inconsistent", "exclude", "The corrected seed-2022 direction reverses and preranked GSEA points the other way.", ""),
        "R-MMU-196071_METABOLISM_OF_STEROID_HORMONES": ("metabolic_endocrine", "plausible_context", "supporting", "Renal endocrine metabolism is relevant, but conventional support and project consistency are incomplete.", "Siew2024"),
        "R-MMU-9764274_REGULATION_OF_EXPRESSION_AND_FUNCTION_OF_TYPE_I_CLASSICAL_CADHERINS": ("epithelial_adhesion", "plausible_complementary", "supporting", "Higher cadherin regulation is compatible with epithelial structural adaptation, but ssGSEA does not support the direction.", "Finch2025;Siew2024"),
        "R-MMU-76005_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2": ("vascular_hemostatic", "broad_redundant", "supporting", "Redundant with the platelet-degranulation signal and strongly attenuated after composition-proxy adjustment.", "Siew2024"),
        "R-MMU-3000178_ECM_PROTEOGLYCANS": ("matrix_growth_repair", "coherent_aligned", "main_candidate", "Higher ECM proteoglycan activity agrees with renal matrix dysregulation and profibrotic structural adaptation reported across spaceflight kidney datasets.", "Finch2025;Siew2024"),
        "R-MMU-418594_G_ALPHA_I_SIGNALLING_EVENTS": ("gpcr", "broad_redundant", "exclude", "A broad signaling child without a kidney-specific ligand or downstream process and with incomplete GSEA support.", ""),
        "R-MMU-8978868_FATTY_ACID_METABOLISM": ("metabolic", "seed_inconsistent", "exclude", "The primary positive direction becomes negative in both added seeds, matching known strain-dependent lipid heterogeneity rather than one pooled direction.", "Finch2025"),
        "R-MMU-168256_IMMUNE_SYSTEM": ("immune", "composition_sensitive", "exclude", "The broad immune parent is strongly attenuated after composition-proxy adjustment and does not identify a renal mechanism.", "Finch2025"),
        "R-MMU-174824_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE": ("metabolic", "plausible_context", "supporting", "Lipoprotein handling is relevant to known renal lipid dysregulation, but seed and conventional support are incomplete.", "Finch2025"),
        "R-MMU-500792_GPCR_LIGAND_BINDING": ("gpcr", "broad_redundant", "exclude", "The parent label is too broad to assign a renal biological interpretation.", ""),
        "R-MMU-1266738_DEVELOPMENTAL_BIOLOGY": ("growth_repair", "broad_redundant", "exclude", "The label combines unrelated developmental systems and its median effect is substantially smaller than the primary estimate.", ""),
        "R-MMU-1474228_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX": ("matrix_growth_repair", "composition_sensitive", "supporting", "Matrix degradation is literature-relevant and seed-stable, but only 5% of the unadjusted magnitude remains after composition-proxy adjustment.", "Finch2025;Siew2024"),
        "R-MMU-2173782_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS": ("vascular_immune", "weak_support", "exclude", "Only two of five directional checks support the lower score.", ""),
        "R-MMU-1280218_ADAPTIVE_IMMUNE_SYSTEM": ("immune", "composition_sensitive", "exclude", "The broad adaptive parent is attenuated after composition adjustment and lacks GSEA support.", "Finch2025"),
        "R-MMU-9757110_PREDNISONE_ADME": ("drug_metabolism", "seed_inconsistent", "exclude", "The drug-specific annotation reverses across seeds and is near zero in the median run.", ""),
        "R-MMU-71291_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES": ("metabolic_capacity", "latent_gene_discordant", "exclude", "The latent score is lower across seeds and five of six projects, but raw member genes, ssGSEA, preranked GSEA, and prior kidney GSEA point higher. It is not interpreted as reduced tubular metabolic capacity.", "Hammond2018;Siew2024"),
        "R-MMU-445355_SMOOTH_MUSCLE_CONTRACTION": ("vascular_hemostatic", "composition_sensitive", "exclude", "Likely reflects renal vascular or stromal content and falls below the retained-magnitude threshold after composition adjustment.", "Siew2024"),
        "R-MMU-392499_METABOLISM_OF_PROTEINS": ("metabolic", "broad_redundant", "exclude", "The broad protein-metabolism parent is composition-sensitive and does not resolve a specific process.", "Finch2025"),
        "R-MMU-6805567_KERATINIZATION": ("epithelial_annotation", "tissue_incongruent", "exclude", "Keratinization is not a coherent whole-kidney mechanism and likely reflects shared epithelial genes or sampling composition.", ""),
        "R-MMU-2672351_STIMULI_SENSING_CHANNELS": ("ion_sensing", "seed_inconsistent", "exclude", "The added seeds reverse the primary direction and the median is near zero.", ""),
        "R-MMU-446728_CELL_JUNCTION_ORGANIZATION": ("epithelial_adhesion", "composition_sensitive", "supporting", "Cell-junction change is literature-relevant, but the expiMap effect is almost eliminated by composition-proxy adjustment.", "Finch2025"),
        "R-MMU-168249_INNATE_IMMUNE_SYSTEM": ("immune", "composition_sensitive", "exclude", "The broad innate parent retains only a small fraction of its unadjusted magnitude.", "Finch2025"),
        "R-MMU-2142753_ARACHIDONATE_METABOLISM": ("lipid_inflammatory", "weak_support", "exclude", "Only two directional checks support this isolated lipid-inflammatory child.", "Finch2025"),
        "R-MMU-8856688_GOLGI_TO_ER_RETROGRADE_TRANSPORT": ("intracellular_traffic", "seed_inconsistent", "exclude", "The seed median is small and one added seed reverses direction.", ""),
        "R-MMU-211859_BIOLOGICAL_OXIDATIONS": ("metabolic_capacity", "heterogeneous_subprogram", "supporting", "The latent score is lower across seeds and five of six projects and several high-weight CYP, UGT, and ACSM genes are lower, but the full member set, ssGSEA, and preranked GSEA point higher. This is retained only as heterogeneous enzyme regulation, not a general reduction in renal oxidation.", "Hammond2018;Siew2024"),
        "R-MMU-112315_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES": ("shared_receptor_annotation", "tissue_incongruent", "exclude", "The neuronal label is not interpreted in whole kidney and is largely removed by composition adjustment.", ""),
        "R-MMU-195721_SIGNALING_BY_WNT": ("matrix_growth_repair", "coherent_complementary", "main_candidate", "Higher aggregate WNT program scores across all six projects extend prior gene-level reports of WNT dysregulation; they do not imply that every WNT ligand, including Wnt11, is higher.", "Finch2025;Siew2024"),
    },
    "spleen": {
        "R-MMU-6798695_NEUTROPHIL_DEGRANULATION": ("innate_effector", "coherent_complementary", "main_candidate", "Lower splenic degranulation-program scores are consistent across all five unconfounded projects, seeds, ssGSEA, GSEA, and composition adjustment. This extends prior evidence of altered neutrophil regulation but is not a direct assay of neutrophil degranulation or suppressive function.", "Buchheim2026"),
        "R-MMU-202403_TCR_SIGNALING": ("adaptive_activation", "coherent_aligned", "main_candidate", "Lower T-cell receptor signaling is the strongest adaptive signal and matches prior reduced splenic T-cell abundance and activation after flight.", "Gridley2009;Martinez2015;Hwang2015"),
        "R-MMU-198933_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL": ("immune_coordination", "method_discordant", "supporting", "Lower scores across projects and seeds suggest impaired immune-cell coordination, but preranked GSEA points in the opposite direction and fewer than half of measured member genes are lower. It remains supporting context.", "Gridley2009;Hwang2015"),
        "R-MMU-9709957_SENSORY_PERCEPTION": ("shared_receptor_annotation", "tissue_incongruent", "exclude", "A sensory parent is not a spleen mechanism and lacks conventional support.", ""),
        "R-MMU-373076_CLASS_A_1_RHODOPSIN_LIKE_RECEPTORS": ("gpcr", "tissue_incongruent", "exclude", "The receptor superfamily is statistically repeatable but too broad and off-context to interpret as sensory biology in spleen.", ""),
        "R-MMU-983231_FACTORS_INVOLVED_IN_MEGAKARYOCYTE_DEVELOPMENT_AND_PLATELET_PRODUCTION": ("hematopoietic", "seed_inconsistent", "supporting", "A plausible splenic hematopoietic signal, but direction is not stable across complete training runs.", "Horie2019"),
        "R-MMU-388841_REGULATION_OF_T_CELL_ACTIVATION_BY_CD28_FAMILY": ("adaptive_activation", "method_discordant", "supporting", "Lower CD28-linked scores are stable across all projects and seeds and match prior anti-CD3/CD28 activation deficits, but preranked GSEA points in the opposite direction and fewer than half of measured member genes are lower. It supports, but does not define, the adaptive result.", "Martinez2015;Hwang2015"),
        "R-MMU-2132295_MHC_CLASS_II_ANTIGEN_PRESENTATION": ("adaptive_activation", "context_sensitive", "supporting", "Lower pooled antigen presentation supports the adaptive story, but only three of five unconfounded projects share the direction.", "Hwang2015"),
        "R-MMU-983705_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR": ("adaptive_activation", "seed_inconsistent", "exclude", "The positive primary score collapses to near zero and reverses across added seeds.", "Horie2019"),
        "R-MMU-9716542_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3": ("cytoskeletal_mitochondrial", "project_sensitive", "exclude", "The effect is seed-stable but fails held-out-project direction support and the parent label is broad.", ""),
        "R-MMU-163125_POST_TRANSLATIONAL_MODIFICATION_SYNTHESIS_OF_GPI_ANCHORED_PROTEINS": ("surface_protein_processing", "plausible_complementary", "supporting", "Lower GPI-anchor synthesis could affect immune-surface proteins, but the pathway is general processing rather than a specific immune mechanism.", ""),
        "R-MMU-2990846_SUMOYLATION": ("protein_regulation", "project_sensitive", "exclude", "Only one of five held-out projects supports the primary positive direction and seeds are inconsistent.", ""),
        "R-MMU-2980736_PEPTIDE_HORMONE_METABOLISM": ("endocrine_metabolic", "project_sensitive", "exclude", "The label is broad and does not reach held-out-project support.", ""),
        "R-MMU-5669034_TNFS_BIND_THEIR_PHYSIOLOGICAL_RECEPTORS": ("inflammatory_signaling", "seed_inconsistent", "supporting", "Lower TNF-receptor engagement is method-supported and project-consistent, but one added seed reverses the effect.", "Gridley2009;Hwang2015"),
        "R-MMU-76002_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION": ("hematopoietic", "plausible_context", "supporting", "Lower platelet activation may reflect splenic hematologic change, but conventional support is incomplete and the median effect is modest.", "Horie2019"),
        "R-MMU-422475_AXON_GUIDANCE": ("shared_cytoskeletal_annotation", "tissue_incongruent", "exclude", "Axon-guidance genes have immune and cytoskeletal roles, but this parent cannot be interpreted as neuronal remodeling in spleen.", ""),
        "R-MMU-1474290_COLLAGEN_FORMATION": ("stromal_matrix", "small_after_multiseed", "supporting", "A possible stromal response, but the three-seed median is small and no direct spleen-spaceflight validation establishes the mechanism.", ""),
        "R-MMU-9020702_INTERLEUKIN_1_SIGNALING": ("inflammatory_signaling", "context_sensitive", "supporting", "Lower IL-1 signaling is seed-stable but only three of five projects support the direction.", "Gridley2009;Hwang2015"),
        "R-MMU-556833_METABOLISM_OF_LIPIDS": ("immune_metabolism", "plausible_complementary", "supporting", "Lower lipid metabolism is internally repeatable and may indicate reduced metabolic support for immune activity, but GSEA does not support the direction.", "Wu2024"),
        "R-MMU-2172127_DAP12_INTERACTIONS": ("innate_activation", "small_after_multiseed", "exclude", "The primary lower effect collapses to approximately zero across seeds despite nominal directional checks.", ""),
        "R-MMU-5621481_C_TYPE_LECTIN_RECEPTORS_CLRS": ("innate_sensing", "coherent_complementary", "main_candidate", "Lower C-type lectin sensing is consistent across all five projects and all five checks, extending the known immune dysfunction phenotype to pathogen-recognition circuitry.", "Hwang2015;Wu2024"),
        "R-MMU-372790_SIGNALING_BY_GPCR": ("gpcr", "seed_inconsistent", "exclude", "The broad GPCR parent reverses across seeds and lacks conventional support.", ""),
        "R-MMU-9748784_DRUG_ADME": ("drug_metabolism", "weak_support", "exclude", "The drug-processing parent has only one of five supporting checks and is not a spleen-specific mechanism.", ""),
        "R-MMU-2454202_FC_EPSILON_RECEPTOR_FCERI_SIGNALING": ("innate_adaptive_receptor", "seed_inconsistent", "exclude", "The positive primary direction is not stable across full training runs.", ""),
        "R-MMU-388396_GPCR_DOWNSTREAM_SIGNALLING": ("gpcr", "broad_redundant", "exclude", "The broad negative GPCR child conflicts with the positive GPCR parent and lacks held-out support.", ""),
        "R-MMU-174824_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE": ("immune_metabolism", "plausible_context", "supporting", "Lower lipoprotein handling is seed-stable but lacks GSEA and composition-independent evidence as a distinct spleen process.", ""),
        "R-MMU-975634_RETINOID_METABOLISM_AND_TRANSPORT": ("immune_metabolism", "seed_inconsistent", "exclude", "The primary positive effect collapses to zero and changes direction across seeds.", ""),
        "R-MMU-418555_G_ALPHA_S_SIGNALLING_EVENTS": ("gpcr", "broad_redundant", "exclude", "The isolated GPCR child is not supported by either conventional method.", ""),
        "R-MMU-9958863_SLC_MEDIATED_TRANSPORT_OF_AMINO_ACIDS": ("immune_metabolism", "plausible_complementary", "supporting", "Higher amino-acid transport is fully triangulated but isolated from the dominant lower immune programs and requires targeted metabolic validation.", ""),
        "R-MMU-9013408_RHOG_GTPASE_CYCLE": ("cytoskeletal", "small_after_multiseed", "exclude", "The primary effect collapses to zero and has only one supporting check.", ""),
        "R-MMU-375276_PEPTIDE_LIGAND_BINDING_RECEPTORS": ("receptor_signaling", "weak_support", "exclude", "Only one directional check supports this broad receptor parent.", ""),
        "R-MMU-936837_ION_TRANSPORT_BY_P_TYPE_ATPASES": ("ion_transport", "small_after_multiseed", "exclude", "The primary effect collapses to zero and has only one supporting check.", ""),
        "R-MMU-6791226_MAJOR_PATHWAY_OF_RRNA_PROCESSING_IN_THE_NUCLEOLUS_AND_CYTOSOL": ("biosynthetic_capacity", "seed_inconsistent", "exclude", "Lower ribosome processing is plausible but not stable across all complete training runs.", ""),
        "R-MMU-6811442_INTRA_GOLGI_AND_RETROGRADE_GOLGI_TO_ER_TRAFFIC": ("intracellular_traffic", "plausible_context", "supporting", "Higher secretory trafficking is seed-stable but lacks held-out-project support as a central spleen response.", ""),
        "R-MMU-68886_M_PHASE": ("cell_cycle", "small_after_multiseed", "exclude", "The primary lower score shrinks substantially and lacks seed stability.", ""),
    },
}


LITERATURE_SOURCES = (
    ("Finch2025", "Finch et al. Spaceflight causes strain-dependent gene expression changes in the kidneys of mice.", "https://doi.org/10.1038/s41526-025-00465-0"),
    ("Siew2024", "Siew et al. Cosmic kidney disease: an integrated pan-omic, physiological and morphological study into spaceflight-induced renal dysfunction.", "https://doi.org/10.1038/s41467-024-49212-1"),
    ("Hammond2018", "Hammond et al. Effects of space flight on mouse liver versus kidney: gene pathway analyses.", "https://doi.org/10.3390/ijms19124106"),
    ("Horie2019", "Horie et al. Down-regulation of GATA1-dependent erythrocyte-related genes in the spleens of mice exposed to space travel.", "https://doi.org/10.1038/s41598-019-44067-9"),
    ("Gridley2009", "Gridley et al. Spaceflight effects on T lymphocyte distribution, function and gene expression.", "https://doi.org/10.1152/japplphysiol.91126.2008"),
    ("Martinez2015", "Martinez et al. Spaceflight and simulated microgravity cause a significant reduction of key gene expression in early T-cell activation.", "https://doi.org/10.1152/ajpregu.00449.2014"),
    ("Hwang2015", "Hwang et al. Post-spaceflight mouse splenocytes demonstrate altered activation properties and surface molecule expression.", "https://doi.org/10.1371/journal.pone.0124380"),
    ("Wu2024", "Wu et al. Single-cell analysis identifies conserved features of immune dysfunction in simulated microgravity and spaceflight.", "https://doi.org/10.1038/s41467-023-42013-y"),
    ("Buchheim2026", "Buchheim et al. Spaceflight alters the immune regulatory functions of neutrophil granulocytes on T lymphocytes.", "https://doi.org/10.1016/j.isci.2025.114380"),
)


def decode(values) -> list[str]:
    return [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def gene_symbol_lookup() -> dict[str, str]:
    with h5py.File(ROOT / "assets/archs4/mouse_gene_v2.5.h5", "r") as handle:
        genes = decode(handle["/meta/genes/ensembl_gene"][:])
        symbols = decode(handle["/meta/genes/symbol"][:])
    return dict(zip(genes, symbols))


def build_manual_review(matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tissue, mapping in MANUAL_REVIEW.items():
        observed = set(
            matrix.loc[
                matrix["tissue"].eq(tissue) & matrix["primary_top_decile"], "term"
            ]
        )
        expected = set(mapping)
        if observed != expected:
            raise RuntimeError(
                f"Manual {tissue} review is incomplete: missing={sorted(observed - expected)}, "
                f"extra={sorted(expected - observed)}"
            )
        for term, (family, category, decision, rationale, citations) in mapping.items():
            rows.append(
                {
                    "tissue": tissue,
                    "term": term,
                    "process_family": family,
                    "manual_category": category,
                    "decision": decision,
                    "manual_rationale": rationale,
                    "literature_keys": citations,
                }
            )
    review = matrix.merge(pd.DataFrame(rows), on=["tissue", "term"], how="inner")
    return review.sort_values(
        ["tissue", "decision", "robustness_support_count", "effect_seed2020"],
        ascending=[True, True, False, False],
    )


def gene_effects_for_data(data) -> pd.DataFrame:
    expression = pd.DataFrame(data.log2cpm, index=data.obs.index, columns=data.genes)
    condition = data.obs["condition_inferred"].astype(str)
    rows = []
    for accession, indexes in data.obs.groupby("id.accession", observed=True).indices.items():
        indexes = np.asarray(indexes)
        local = condition.iloc[indexes]
        flight = local.eq("flight").to_numpy()
        ground = local.eq("ground_control").to_numpy()
        if flight.any() and ground.any():
            effect = (
                expression.iloc[indexes].iloc[flight].mean(axis=0)
                - expression.iloc[indexes].iloc[ground].mean(axis=0)
            )
            rows.append(effect.rename(str(accession)))
    return pd.DataFrame(rows)


def decoder_weight_matrix(config, seed: int) -> pd.DataFrame:
    import scarches as sca

    query_dir = config.query_dir(seed)
    mapped = ad.read_h5ad(query_dir / "mapped_query_with_scores.h5ad")
    model = sca.models.EXPIMAP.load(query_dir / "query_model", adata=mapped)
    model.latent_directions(method="sum", adata=model.adata)
    terms = list(map(str, model.adata.uns["terms"]))
    directions = np.asarray(model.adata.uns["directions"], dtype=float)
    state = model.model.state_dict()
    weights = state["decoder.L0.expr_L.weight"].detach().cpu().numpy()
    mask = state["decoder.L0.expr_L.mask"].detach().cpu().numpy()
    oriented = weights * mask * directions[np.newaxis, :]
    return pd.DataFrame(
        oriented,
        index=mapped.var_names.astype(str),
        columns=terms,
    )


def build_gene_support(
    review: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    symbols = gene_symbol_lookup()
    detail_rows = []
    seed_summary_rows = []
    for config in CONFIGS:
        data = load_data(config)
        accession_effects = gene_effects_for_data(data)
        balanced = accession_effects.mean(axis=0)
        selected = review.loc[review["tissue"].eq(config.tissue)]
        for seed in SEEDS:
            decoder = decoder_weight_matrix(config, seed)
            effect_column = f"effect_seed{seed}"
            for item in selected.itertuples(index=False):
                pathway_effect = float(getattr(item, effect_column))
                if not np.isfinite(pathway_effect) or item.term not in decoder:
                    continue
                genes = [
                    gene
                    for gene in data.gene_sets[item.term]
                    if gene in balanced and gene in decoder.index
                ]
                effects = balanced.loc[genes]
                weights = decoder.loc[genes, item.term]
                predicted = weights * pathway_effect
                pathway_match = np.sign(effects) == np.sign(pathway_effect)
                decoder_match = np.sign(effects) == np.sign(predicted)
                absolute_weight = weights.abs()
                weight_total = float(absolute_weight.sum())
                weighted_match = (
                    float(absolute_weight.loc[decoder_match].sum() / weight_total)
                    if weight_total > 0
                    else np.nan
                )
                correlation = stats.spearmanr(predicted, effects).statistic

                rank_order = effects.abs().sort_values(ascending=False).index
                ranks = pd.Series(
                    np.arange(1, len(rank_order) + 1), index=rank_order
                )
                for gene in genes:
                    detail_rows.append(
                        {
                            "tissue": config.tissue,
                            "term": item.term,
                            "display_term": item.display_term,
                            "decision": item.decision,
                            "seed": seed,
                            "gene_rank_by_absolute_effect": int(ranks.loc[gene]),
                            "gene_id": gene,
                            "gene_symbol": symbols.get(gene, ""),
                            "project_balanced_log2cpm_effect": float(effects.loc[gene]),
                            "oriented_decoder_weight": float(weights.loc[gene]),
                            "predicted_decoder_logit_contribution": float(
                                predicted.loc[gene]
                            ),
                            "same_direction_as_pathway_score": bool(
                                pathway_match.loc[gene]
                            ),
                            "matches_decoder_predicted_direction": bool(
                                decoder_match.loc[gene]
                            ),
                        }
                    )

                concordant = pd.DataFrame(
                    {
                        "effect": effects,
                        "weight": weights,
                        "predicted": predicted,
                        "match": decoder_match,
                    }
                )
                concordant = concordant.loc[concordant["match"]].sort_values(
                    "predicted", key=lambda values: values.abs(), ascending=False
                )
                top = "; ".join(
                    f"{symbols.get(gene, gene)} (raw {row.effect:+.2f}, weight {row.weight:+.2f})"
                    for gene, row in concordant.head(8).iterrows()
                )
                seed_summary_rows.append(
                    {
                        "tissue": config.tissue,
                        "term": item.term,
                        "display_term": item.display_term,
                        "decision": item.decision,
                        "manual_category": item.manual_category,
                        "seed": seed,
                        "pathway_effect": pathway_effect,
                        "pathway_genes_in_hvg_model": int(len(genes)),
                        "member_gene_same_direction_fraction": float(
                            pathway_match.mean()
                        ),
                        "median_member_gene_log2cpm_effect": float(effects.median()),
                        "decoder_direction_match_fraction": float(
                            decoder_match.mean()
                        ),
                        "decoder_abs_weight_direction_match_fraction": weighted_match,
                        "decoder_gene_effect_spearman": float(correlation),
                        "top_decoder_concordant_member_genes": top,
                    }
                )

    detail = pd.DataFrame(detail_rows)
    by_seed = pd.DataFrame(seed_summary_rows)
    group_columns = [
        "tissue",
        "term",
        "display_term",
        "decision",
        "manual_category",
    ]
    summary = (
        by_seed.groupby(group_columns, as_index=False, observed=True)
        .agg(
            pathway_seed_median_effect=("pathway_effect", "median"),
            pathway_genes_in_hvg_model=("pathway_genes_in_hvg_model", "max"),
            member_gene_same_direction_fraction=(
                "member_gene_same_direction_fraction",
                "median",
            ),
            median_member_gene_log2cpm_effect=(
                "median_member_gene_log2cpm_effect",
                "median",
            ),
            decoder_direction_match_fraction_median=(
                "decoder_direction_match_fraction",
                "median",
            ),
            decoder_abs_weight_direction_match_fraction_median=(
                "decoder_abs_weight_direction_match_fraction",
                "median",
            ),
            decoder_abs_weight_direction_match_fraction_minimum=(
                "decoder_abs_weight_direction_match_fraction",
                "min",
            ),
            decoder_gene_effect_spearman_median=(
                "decoder_gene_effect_spearman",
                "median",
            ),
        )
    )
    primary_top = by_seed.loc[
        by_seed["seed"].eq(SEEDS[0]),
        group_columns + ["top_decoder_concordant_member_genes"],
    ]
    summary = summary.merge(primary_top, on=group_columns, how="left")
    return detail, by_seed, summary


def write_review_notes(review: pd.DataFrame, gene_summary: pd.DataFrame) -> None:
    lines = [
        "# Biological review of corrected kidney and spleen models",
        "",
        "Every pathway in the primary top decile was manually reviewed as a complete Reactome program. Magnitude, direction across projects and seeds, conventional enrichment, composition sensitivity, tissue fit, redundancy, and prior literature were considered separately. A favorable numerical score alone was not enough to promote a pathway.",
        "",
        "## Kidney",
        "",
        "The coherent kidney result is a higher structural and growth-factor-response axis: ECM proteoglycans, aggregate WNT signaling, and IGF transport are higher in flight. These programs agree with prior reports of renal ECM dysregulation, fibrosis-related signaling, nephron remodeling, and WNT involvement, while their joint expression adds a pathway-level repair or maladaptive-growth hypothesis. All three have strong member-gene directional support, but none passes conventional preranked-GSEA FDR below 0.05 and each is attenuated by broad composition-proxy adjustment, so the axis remains a corroborating or complementary hypothesis rather than a statistically confirmed pathway discovery. Lower amino-acid-metabolism and biological-oxidation latent scores are not promoted as broad decreases. Amino-acid metabolism is rejected because raw genes, ssGSEA, preranked GSEA, and prior kidney GSEA point higher. Biological oxidation is retained only as a heterogeneous enzyme subset because several influential CYP, UGT, and ACSM genes are lower while the full gene set points higher. Platelet degranulation is numerically robust but remains a vascular or blood-composition signal rather than a kidney-cell claim.",
        "",
        "## Spleen",
        "",
        "The coherent spleen result is lower adaptive activation plus lower innate pathogen-sensing and degranulation programs. T-cell receptor signaling, neutrophil degranulation, and C-type lectin receptor signaling are lower in all five unconfounded projects, all three seeds, both conventional methods, and composition-adjusted analyses; each also passes conventional preranked-GSEA FDR below 0.05 and has majority member-gene support. This recovers prior T-cell activation deficits and adds a coordinated innate-effector layer. The transcriptomic degranulation score does not directly measure neutrophil function. CD28-linked activation and lymphoid-to-nonlymphoid interaction scores are lower but disagree with preranked GSEA, while MHC class II antigen presentation is directionally mixed across projects; these remain supporting context rather than central claims.",
        "",
        "## Scope recommendation",
        "",
        "Spleen should no longer be described as a junk or null tissue. Its corrected HVG reference-query model has a coherent immune-suppression pattern across five unconfounded projects and is a stronger manuscript candidate than soleus. Kidney also has a coherent three-program structural and growth-factor axis, but it is composition-sensitive and lacks conventional pathway FDR support; it is suitable as a secondary corroborating and pathway-integration result with conservative language, not as an equal-strength discovery claim. Lung and retina remain the actual null tissues from the original screen.",
        "",
        "## Main-candidate member-gene direction",
        "",
    ]
    for tissue in ("kidney", "spleen"):
        lines.extend([f"### {tissue.title()}", ""])
        selected = gene_summary.loc[
            gene_summary["tissue"].eq(tissue)
            & gene_summary["decision"].eq("main_candidate")
        ]
        for row in selected.itertuples(index=False):
            lines.append(
                f"- **{row.display_term}:** {row.member_gene_same_direction_fraction:.0%} of measured member genes move with the pathway score. Across the three trained decoders, a median {row.decoder_abs_weight_direction_match_fraction_median:.0%} of absolute decoder weight predicts the observed member-gene direction (minimum {row.decoder_abs_weight_direction_match_fraction_minimum:.0%}). Largest concordant seed-2020 effects: {row.top_decoder_concordant_member_genes}."
            )
        lines.append("")
    lines.extend(["## Kidney metabolic-label audit", ""])
    for term in (
        "R-MMU-71291_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
        "R-MMU-211859_BIOLOGICAL_OXIDATIONS",
    ):
        row = gene_summary.loc[
            gene_summary["tissue"].eq("kidney") & gene_summary["term"].eq(term)
        ].iloc[0]
        lines.append(
            f"- **{row.display_term}:** the median latent shift is {row.pathway_seed_median_effect:+.3f}, but only {row.member_gene_same_direction_fraction:.0%} of measured genes move in that same direction. The median decoder-weighted gene-direction match is {row.decoder_abs_weight_direction_match_fraction_median:.0%} across seeds, and the median decoder-predicted versus observed gene-effect correlation is {row.decoder_gene_effect_spearman_median:.2f}."
        )
    lines.append("")
    lines.extend(["## Literature", ""])
    for key, citation, url in LITERATURE_SOURCES:
        lines.append(f"- **{key}:** [{citation}]({url})")
    lines.append("")
    (OUTPUT_DIR / "biological_review.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def plot_main_candidates(review: pd.DataFrame) -> None:
    label_overrides = {
        "R-MMU-195721_SIGNALING_BY_WNT": "WNT signaling",
        "R-MMU-3000178_ECM_PROTEOGLYCANS": "ECM proteoglycans",
        "R-MMU-381426_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS": "IGF transport and uptake",
        "R-MMU-202403_TCR_SIGNALING": "T-cell receptor signaling",
        "R-MMU-6798695_NEUTROPHIL_DEGRANULATION": "Neutrophil degranulation program",
        "R-MMU-5621481_C_TYPE_LECTIN_RECEPTORS_CLRS": "C-type lectin receptor signaling",
    }
    category_colors = {
        "coherent_aligned": "#2f855a",
        "coherent_complementary": "#2b6cb0",
    }
    accession = pd.read_csv(OUTPUT_DIR / "seed_accession_effects.tsv.gz", sep="\t")
    accession = accession.loc[
        accession["seed"].eq(SEEDS[0])
        & ~(
            accession["tissue"].eq("spleen")
            & accession["accession"].astype(str).eq("OSD-288")
        )
    ]
    project = (
        accession.groupby(["tissue", "project", "term"], as_index=False)["effect"]
        .mean()
    )
    selected = review.loc[review["decision"].eq("main_candidate")].copy()
    order = {
        "kidney": [
            "R-MMU-3000178_ECM_PROTEOGLYCANS",
            "R-MMU-195721_SIGNALING_BY_WNT",
            "R-MMU-381426_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS",
        ],
        "spleen": [
            "R-MMU-202403_TCR_SIGNALING",
            "R-MMU-6798695_NEUTROPHIL_DEGRANULATION",
            "R-MMU-5621481_C_TYPE_LECTIN_RECEPTORS_CLRS",
        ],
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.2), constrained_layout=True)
    for axis, tissue in zip(axes, ("kidney", "spleen")):
        tissue_review = selected.set_index("term").loc[order[tissue]].reset_index()
        positions = np.arange(len(tissue_review))[::-1]
        all_values = []
        for position, item in zip(positions, tissue_review.itertuples(index=False)):
            effects = project.loc[
                project["tissue"].eq(tissue) & project["term"].eq(item.term),
                "effect",
            ].to_numpy(dtype=float)
            all_values.extend(effects.tolist())
            color = category_colors[item.manual_category]
            axis.scatter(
                effects,
                np.full(len(effects), position),
                color="#737b83",
                edgecolor="white",
                linewidth=0.7,
                s=54,
                alpha=0.85,
                zorder=2,
            )
            axis.hlines(
                position,
                item.seed_effect_minimum,
                item.seed_effect_maximum,
                color=color,
                linewidth=3,
                zorder=3,
            )
            axis.scatter(
                item.seed_effect_median,
                position,
                marker="D",
                color=color,
                edgecolor="white",
                linewidth=0.9,
                s=92,
                zorder=4,
            )
            axis.text(
                0.99,
                position,
                f"GSEA q={item.gsea_fdr:.3f}",
                transform=axis.get_yaxis_transform(),
                ha="right",
                va="bottom",
                fontsize=9,
                color="#4a5157",
            )
        axis.axvline(0, color="#343a40", linewidth=1)
        axis.set_yticks(positions)
        axis.set_yticklabels(
            [label_overrides[term] for term in tissue_review["term"]], fontsize=11
        )
        for tick, category in zip(
            axis.get_yticklabels(), tissue_review["manual_category"]
        ):
            tick.set_color(category_colors[category])
            tick.set_fontweight("semibold")
        axis.set_title(tissue.title(), fontsize=16, fontweight="bold", pad=12)
        axis.set_xlabel("Flight minus ground expiMap pathway shift", fontsize=11)
        axis.grid(axis="x", color="#d9dde1", linewidth=0.8)
        axis.set_axisbelow(True)
        values = np.asarray(
            all_values
            + tissue_review["seed_effect_minimum"].tolist()
            + tissue_review["seed_effect_maximum"].tolist(),
            dtype=float,
        )
        span = max(float(np.ptp(values)), 0.1)
        axis.set_xlim(float(values.min() - 0.12 * span), float(values.max() + 0.38 * span))
        axis.spines[["top", "right"]].set_visible(False)

    legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#737b83",
            markeredgecolor="white",
            markersize=8,
            label="OSDR project, seed 2020",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="#4a5157",
            linewidth=2,
            markerfacecolor="#4a5157",
            markersize=7,
            label="Median and range across three complete trainings",
        ),
        Line2D([0], [0], color=category_colors["coherent_aligned"], linewidth=4, label="Prior-literature aligned"),
        Line2D([0], [0], color=category_colors["coherent_complementary"], linewidth=4, label="Complementary process hypothesis"),
    ]
    fig.legend(
        handles=legend,
        loc="outside lower center",
        ncol=2,
        frameon=False,
        fontsize=10,
    )
    fig.suptitle(
        "Corrected kidney and spleen HVG models: manually retained pathway programs",
        fontsize=18,
        fontweight="bold",
    )
    for suffix in ("png", "pdf"):
        fig.savefig(
            OUTPUT_DIR / f"curated_main_pathway_project_seed_evidence.{suffix}",
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    matrix = pd.read_csv(OUTPUT_DIR / "pathway_evidence_matrix.tsv", sep="\t")
    review = build_manual_review(matrix)
    detail, by_seed, summary = build_gene_support(review)
    review.to_csv(OUTPUT_DIR / "top_decile_manual_review.tsv", sep="\t", index=False)
    detail.to_csv(
        OUTPUT_DIR / "top_decile_member_gene_effects.tsv", sep="\t", index=False
    )
    by_seed.to_csv(
        OUTPUT_DIR / "top_decile_member_gene_support_by_seed.tsv",
        sep="\t",
        index=False,
    )
    summary.to_csv(
        OUTPUT_DIR / "top_decile_member_gene_support.tsv", sep="\t", index=False
    )
    main_terms = review.loc[
        review["decision"].eq("main_candidate"), ["tissue", "term"]
    ]
    main_detail = detail.merge(main_terms, on=["tissue", "term"], how="inner")
    main_summary = summary.loc[summary["decision"].eq("main_candidate")]
    main_detail.to_csv(
        OUTPUT_DIR / "main_candidate_member_gene_effects.tsv", sep="\t", index=False
    )
    main_summary.to_csv(
        OUTPUT_DIR / "main_candidate_member_gene_support.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        LITERATURE_SOURCES, columns=["key", "citation", "url"]
    ).to_csv(OUTPUT_DIR / "literature_sources.tsv", sep="\t", index=False)
    write_review_notes(review, summary)
    plot_main_candidates(review)
    print(
        review.groupby(["tissue", "decision"], observed=True)
        .size()
        .rename("pathways")
        .to_string(),
        flush=True,
    )


if __name__ == "__main__":
    main()
