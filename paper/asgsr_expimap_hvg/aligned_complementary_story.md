# Aligned and Complementary Program Audit

> **Scope update (July 14, 2026):** This document records the original thymus, skin, liver, and soleus audit. The revised manuscript promotes corrected spleen as a main positive tissue, adds kidney as secondary exploratory evidence, and moves soleus to supplementary sensitivity analysis. Use `tissue_selection_audit.md`, `source_data/table_2_retained_pathway_evidence.tsv`, and the manuscript for the current interpretation. The historical discussion below is retained as provenance rather than silently rewritten.

## Bottom line

There is a clear alternative perspective, but it is not that one new pathway explains spaceflight across all tissues.

**The defensible cross-tissue story is that spaceflight changes the multicellular coordination and structural-maintenance state of tissues, in addition to producing the established intracellular stress and metabolic phenotypes.** The direction and biological form of that change are tissue specific:

- **Thymus:** an adaptive, proliferative, and migratory thymocyte-niche state is lower, while innate and stromal-matrix programs are higher.
- **Skin:** communication, regenerative-niche, antioxidant, and barrier-support programs move together in a lower direction.
- **Liver:** established metabolic dysfunction is accompanied by lower adaptive-immune and cytoskeletal-regulatory scores; two previously blue liver terms do not survive deeper review.
- **Soleus:** matrix degradation and damage-response scores suggest an extracellular structural-injury layer beyond contractile atrophy, but this is the weakest and most mission-sensitive component.

This should be presented as a **testable tissue-state and cell-composition model**, not as proof that the named pathways causally regulate one another.

## Audit method

Every program that was labeled aligned or complementary before this review was evaluated using:

1. decoder-oriented flight-minus-ground score direction;
2. accession agreement and restricted sensitivity;
3. member-gene differential-expression support;
4. direct tissue-specific spaceflight evidence;
5. related spaceflight evidence outside the tissue;
6. non-spaceflight tissue biology needed to connect the Reactome term to the proposed phenotype;
7. gene-set overlap and within-accession pathway-score correlation.

The complete quantitative audit is in `source_data/table_s6_program_story_audit.tsv`. Pairwise score correlations and Reactome gene overlaps are in `source_data/table_s7_program_pairwise_structure.tsv`.

Evidence tiers:

- **Direct:** tissue-specific spaceflight work supports the process and compatible direction.
- **Process-supported:** spaceflight work supports the biological process, but not the exact Reactome term or direction.
- **Mechanistic only:** tissue biology makes the interpretation plausible, but direct spaceflight validation is absent.
- **Direction conflict or ambiguity:** existing evidence prevents a clean complementary claim.

## Why these should be treated as tissue-state axes

Many selected programs covary strongly even when their Reactome gene sets barely overlap:

- Thymus mitotic cell cycle and DNA repair: `r=0.996`, with gene-set Jaccard `0.123`.
- Thymus DNA repair and lymphoid-stromal interactions: `r=0.937`, with no shared genes.
- Skin chromatin modifiers and gap-junction trafficking: `r=0.957`, with no shared genes.
- Skin phase II detoxification and sphingolipid metabolism: `r=0.931`, with no shared genes.
- Liver MHC class II presentation and Rho-family cycle: `r=0.865`, with four shared genes.
- Soleus extracellular-matrix degradation and DNA repair: `r=0.900`, with no shared genes.

These are correlations of sample-level pathway scores after centering within accession. They support coordinated tissue states, but could also reflect cell composition, global condition effects, or shared model structure. They do not establish molecular coupling.

## Thymus

### Clear story

Prior work establishes thymic shrinkage and lower proliferative activity. The expiMap result adds a spatial and organizational interpretation: **flight thymus shifts away from a proliferating, T-cell-producing, migratory niche and toward an innate and stromal response**.

The lower cell-cycle, DNA-repair-associated, T-cell receptor, lymphoid-stromal, and RHOA scores form a highly correlated axis. Higher innate Toll-like receptor and extracellular-matrix organization scores oppose that axis. All complementary thymus directions survive exclusion of strain-confounded OSD-289.

This does not prove that thymic epithelial cells remodel matrix or that TLR signaling causes involution. In bulk tissue, the same pattern could arise from loss of proliferating thymocytes plus a larger relative contribution from stromal or innate cells. That composition hypothesis is itself biologically meaningful and directly testable.

| Program | Direction | Evidence | Story decision | Interpretation |
| --- | --- | --- | --- | --- |
| Mitotic cell cycle | Lower | Direct | Core aligned anchor | Reduced thymocyte proliferation; may partly reflect fewer cycling thymocytes. |
| DNA repair | Lower | Process-supported | Core | Reduced repair-associated tissue state; not proof that repair kinetics are defective. |
| T-cell receptor signaling | Lower | Direct process | Core aligned anchor | Reduced adaptive thymocyte state or abundance. |
| Innate TLR signaling | Higher | Indirect spaceflight | Core complementary | Candidate innate counter-response; direct thymus TLR validation is needed. |
| Lymphoid-stromal interactions | Lower | Mechanistic only | Core complementary | Weaker niche coordination or loss of interacting populations. |
| RHOA cytoskeletal cycle | Lower | Mechanistic only | Supporting | Supports reduced motility and adhesion; not an independent mechanism. |
| Extracellular-matrix organization | Higher | Process-supported | Core complementary | Stromal structural response; not equivalent to fibrosis or increased matrix mass. |

Key sources: [Horie et al. 2019](https://www.nature.com/articles/s41598-019-56432-9), [Gridley et al. 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075097), [Lebsack et al. 2010](https://doi.org/10.1002/jcb.22547), [Shimizu et al. 2023](https://www.nature.com/articles/s42003-023-05251-w), [Grandke et al. 2026](https://www.nature.com/articles/s41467-026-68737-1), and [thymic ECM and migration biology](https://pubmed.ncbi.nlm.nih.gov/11097218/).

## Skin

### Clear story

Prior literature emphasizes inflammation, oxidative stress, DNA damage, dermal thinning, and barrier disruption. The alternative expiMap perspective is more specific: **flight skin generally occupies a lower epidermal-differentiation and coordinated-maintenance state involving keratinization, direct cell communication, lipid barrier support, chromatin regulation, and follicular repair**.

Keratinization is lower in five of six accessions (study-balanced effect -1.22), with 13 of 70 retained member genes passing pooled gene-level FDR. Cope et al. independently reported lower filaggrin and CASP14 in most murine skin subsets, making this an aligned barrier-differentiation result rather than junk. Reviewed skin program pairs correlate at `r=0.83-0.99`, including pathway pairs with no shared genes. The finding is therefore stronger as one coordinated state than as separate pathway mechanisms.

Protocol stratification narrows the claim. Both MHU-2 microgravity contrasts are lower than Earth 1 g across all eight reviewed programs. Onboard artificial 1 g remains lower in dorsal skin but is mostly higher in femoral skin, so it is not a general rescue. RR-5 dorsal skin after 30 days of Earth recovery is the sole positive keratinization accession. RR-6 live-return and ISS-terminal directions agree, while RR-7 shows a duration and strain interaction. The general story is lower under several exposure-phase contexts, modified by site, recovery, strain, and duration.

The aligned cell-junction result provides the strongest structural anchor. Gap-junction trafficking extends the communication phenotype. Hedgehog and sphingolipid metabolism provide biologically specific regenerative and barrier hypotheses, but direct flight measurements of those pathways are absent. Lower chromatin-modifier score is compatible with altered epigenetic plasticity, but existing flight-skin studies do not establish a global loss of chromatin-modifying activity. Broad phase-II conjugation is lower while its fully nested glutathione-conjugation child is higher; their strongly anticorrelated latent scores do not identify a uniform detoxification direction.

| Program | Direction | Evidence | Story decision | Interpretation |
| --- | --- | --- | --- | --- |
| Keratinization | Lower | Direct | Core aligned anchor | Reduced epidermal differentiation and cornified-barrier support; not uniform downregulation of every keratin gene. |
| Cell-cell junction organization | Lower | Direct | Core aligned anchor | Reduced barrier and tissue coordination. |
| Chromatin-modifying enzymes | Lower | Direct process, direction unresolved | Supporting | Candidate reduction in regulatory plasticity; not established globally. |
| Gap-junction trafficking | Lower | Process-supported | Core complementary | More specific cell-communication hypothesis within known junction disruption. |
| Phase II detoxification | Lower broad node; higher nested glutathione-conjugation node | Direct process, direction unresolved | Context-sensitive; not core | Overlapping latent nodes oppose one another, so detoxification capacity cannot be assigned a direction. |
| Hedgehog signaling | Lower | Mechanistic only | Supporting | Candidate loss of follicular or regenerative-niche activity. |
| Sphingolipid metabolism | Lower | Mechanistic only | Supporting | Candidate reduction in ceramide and barrier-lipid support; metabolite abundance was not measured here. |

Key sources: [Park et al. 2024](https://www.nature.com/articles/s41467-024-48625-2), [Cope et al. 2024](https://www.nature.com/articles/s43856-024-00532-9), [Mao et al. 2014](https://pubmed.ncbi.nlm.nih.gov/24796731/), [Neutelings et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC5515501/), [flight-induced skin epigenetic plasticity](https://pmc.ncbi.nlm.nih.gov/articles/PMC11647166/), [Hedgehog-dependent follicular regeneration](https://pmc.ncbi.nlm.nih.gov/articles/PMC6249328/), and [sphingolipid barrier biology](https://www.jlr.org/article/S0022-2275%2820%2941978-9/fulltext).

## Liver

### Clear but narrower story

Prior literature already establishes lipid, xenobiotic, sulfur, mitochondrial, and insulin-related disruption. The complementary result that survives deeper review is: **flight liver adds a lower adaptive-immune and cytoskeletal-regulatory layer to the established metabolic phenotype**.

The primary liver query was remapped after removing OSD-164 and OSD-168, leaving 197 samples from 10 independent cohort sources. OSD-168 duplicates RR-1 and RR-3 cohorts represented by OSD-48 and OSD-137 and includes ERCC technical variants; OSD-164 overlaps OSD-47 animals. The de-duplicated remap preserves every reviewed pathway direction and strengthens the MHC II, T-cell receptor, and Rho-family effects. Cytochrome P450 remains split 5 positive and 5 negative, so de-duplication does not manufacture a clean metabolic direction.

MHC class II antigen presentation and T-cell receptor signaling should be collapsed into one immune axis. They likely report altered resident or infiltrating immune-cell abundance as much as regulation within individual cells. Lower Rho-family score suggests altered cytoskeletal or mechanical regulation, but does not prove reduced force sensing or identify shear as the cause.

Glutathione conjugation does not support a clean blue claim. Direct liver studies report lower glutathione pools and reducibility, which conflicts with a simple “higher detoxification” interpretation. Extracellular-matrix organization remains context-sensitive, but the de-duplicated result now supports a narrower hypothesis: coordinated matrix organization or maintenance is predominantly lower even though total matrix abundance, turnover, and fibrosis are unresolved.

| Program | Direction | Evidence | Story decision | Interpretation |
| --- | --- | --- | --- | --- |
| Regulation of insulin secretion | Lower | Direct | Core aligned anchor | Lower insulin-regulatory state; the term does not imply hepatocyte insulin secretion. |
| Glutathione conjugation | Higher | Direction conflict | Not core; context-sensitive | Compensation is possible, but direct metabolite evidence and study heterogeneity prevent a clear claim. |
| MHC class II antigen presentation | Lower | Indirect spaceflight | Core complementary | Lower hepatic immune communication or fewer antigen-presenting cells. |
| T-cell receptor signaling | Lower | Indirect spaceflight | Core complementary | Same immune-composition axis as MHC II, not a separate mechanism. |
| Rho-family GTPase cycle | Lower | Process-supported | Supporting | Candidate cytoskeletal-mechanical regulation layer. |
| Extracellular-matrix organization | Predominantly lower | Process-supported, measurement-limited | Context-sensitive hypothesis | Suggests reduced coordinated matrix organization or maintenance; does not establish lower matrix abundance or fibrosis. |

Key sources: [Mathyk et al. 2024](https://www.nature.com/articles/s42003-023-05213-2), [Kurosawa et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8575787/), [Jonscher et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4838331/), [spaceflight antigen-presentation dysfunction](https://pmc.ncbi.nlm.nih.gov/articles/PMC4430214/), and [Li et al. 2026 mechanotransduction study](https://pubmed.ncbi.nlm.nih.gov/41672847/).

## Soleus

### Follow-up hypothesis, not a settled story

Prior literature emphasizes atrophy, slow-to-fast fiber transition, mitochondrial dysfunction, and altered contraction. The complementary expiMap hypothesis is: **spaceflight soleus may also enter a matrix-disassembly and cellular-damage state**.

Higher extracellular-matrix degradation, lower glycosaminoglycan metabolism, and higher DNA repair all retain direction across three accessions and after the OSD-714 restriction. Broad immune signaling is positive in ISS-terminal OSD-104 and live-return OSD-714 but negative in live-return OSD-770, so collection endpoint alone does not explain it. Only two unconfounded accessions remain, member-gene support is weak, and broad contractile, metabolic, and immune directions are mission sensitive. The matrix interpretation is therefore useful for experiment design but is not a general conclusion.

| Program | Direction | Evidence | Story decision | Interpretation |
| --- | --- | --- | --- | --- |
| Extracellular-matrix degradation | Higher | Process-supported | Supporting | Candidate matrix turnover; compatible with reduced ECM proteins in astronaut soleus. |
| Glycosaminoglycan metabolism | Lower | Mechanistic only | Exploratory | Candidate loss of matrix hydration and organization; no direct soleus flight validation. |
| DNA repair | Higher | Indirect spaceflight | Exploratory | Candidate cellular-damage response; no soleus-specific validation. |

Key sources: [Gambara et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5226721/), [Tascher et al. 2017](https://pubs.acs.org/doi/10.1021/acs.jproteome.7b00201), [human soleus ECM remodeling](https://pmc.ncbi.nlm.nih.gov/articles/PMC9962627/), and [Murgia et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11153545/).

## Conference-ready alternative narrative

> Existing studies describe the outcomes of spaceflight in each organ: thymic atrophy, skin inflammation and barrier damage, liver metabolic dysfunction, and soleus muscle loss. Our expiMap results suggest an additional organizational layer. In thymus, the state supporting proliferating and migrating T cells is lower while innate and stromal programs are higher. In skin, keratinization, communication, regenerative, and barrier-lipid programs generally fall together, with gravity, site, recovery, strain, and duration modifying the magnitude; nested conjugation programs do not support a single detoxification direction. De-duplicated liver data add a lower immune and cytoskeletal-regulatory axis to the known metabolic response, while soleus provides preliminary evidence of matrix disassembly. Together, these results suggest that spaceflight affects not only how individual cells handle stress, but also how tissues coordinate cells with their structural environment.

## Falsifiable predictions

1. **Thymus:** spatial or single-cell data should show reduced cycling thymocytes and thymocyte-stromal contact signatures together with relatively increased innate or matrix-producing stromal states.
2. **Skin:** flight samples should show lower functional junctional coupling, altered epidermal ceramides, and reduced follicular regenerative signaling in the same studies; direct metabolite and enzyme assays should determine which conjugation branch accounts for the opposing latent scores.
3. **Liver:** lower MHC II and T-cell receptor scores should localize to altered immune-cell abundance or state; hepatocyte-only data should weaken those signals. Rho-associated changes should correlate with cytoskeletal or shear-response measurements.
4. **Soleus:** flight should increase matrix protease activity or matrix-fragment products and reduce sulfated glycosaminoglycan content; these measurements should vary by mission and strain.

Failure of these predictions would favor bulk-composition or latent-model explanations over coordinated pathway remodeling.
