# Blue / Complementary expiMap Pathway Shifts: Presentation Notes

These notes summarize the blue/complementary FLT-GC pathway shifts from the four literature-reviewed HVG expiMap models. Blue means the pathway shift is biologically plausible and complementary to prior literature, but not a direct, direction-specific replication.

Wording note for slides: "higher in flight" and "lower in flight" refer to the expiMap pathway score / latent pathway shift in FLT relative to GC. This does not mean every gene in the Reactome pathway moves in the same direction.

Citation keys refer to the per-tissue source files in this folder:

- `liver_hvg_sources.md`
- `skin_hvg_sources.md`
- `soleus_hvg_sources.md`
- `thymus_hvg_sources.md`

## Liver HVG

Prior literature frame:

- Mouse liver spaceflight studies repeatedly report lipid dysregulation, lipotoxicity, CYP/xenobiotic metabolism changes, sulfur/cofactor metabolism changes, insulin/estrogen signaling changes, mitochondrial stress, and immune dysfunction.
- Key citation families: `Jonscher2016`, `Beheshti2019`, `Moskaleva2015`, `Kurosawa2021`, `Mathyk2024`, `Vinken2022`, `daSilveira2020`, `Crucian2018`, `Gridley2009`, `Kim2024`, `Li2026`.

Representative blue pathways:

- `R-MMU-9012999_RHO_GTPASE_CYCLE`, `R-MMU-194315_SIGNALING_BY_RHO_GTPASES`, `R-MMU-8980692_RHOA_GTPASE_CYCLE`, `R-MMU-9013149_RAC1_GTPASE_CYCLE`, and `R-MMU-9716542_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3`.
- `R-MMU-2132295_MHC_CLASS_II_ANTIGEN_PRESENTATION`, `R-MMU-983169_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION`, `R-MMU-202403_TCR_SIGNALING`, and `R-MMU-913531_INTERFERON_SIGNALING`.
- `R-MMU-1474228_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX`, `R-MMU-216083_INTEGRIN_CELL_SURFACE_INTERACTIONS`, `R-MMU-418990_ADHERENS_JUNCTIONS_INTERACTIONS`, `R-MMU-421270_CELL_CELL_JUNCTION_ORGANIZATION`, `R-MMU-195721_SIGNALING_BY_WNT`, `R-MMU-9006936_SIGNALING_BY_TGFB_FAMILY_MEMBERS`, `R-MMU-186797_SIGNALING_BY_PDGF`, and `R-MMU-6806834_SIGNALING_BY_MET`.
- `R-MMU-425407_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT`, `R-MMU-9958863_SLC_MEDIATED_TRANSPORT_OF_AMINO_ACIDS`, `R-MMU-196849_METABOLISM_OF_WATER_SOLUBLE_VITAMINS_AND_COFACTORS`, `R-MMU-156590_GLUTATHIONE_CONJUGATION`, and `R-MMU-8978868_FATTY_ACID_METABOLISM`.

Direction-aware presentation wording:

- Rho GTPase cycle is higher in flight, meaning the liver model detects a cytoskeletal/mechanotransduction shift consistent with altered gravity and hepatocyte stress, even though prior liver papers do not directly establish Rho-cycle activation.
- MHC class II antigen presentation, TCR signaling, class I MHC processing, and interferon signaling are lower in flight, meaning the liver blue signal may reflect immune suppression or immune-cell redistribution rather than simple inflammatory activation.
- Extracellular-matrix degradation and SLC-mediated membrane transport are higher in flight, meaning the pathway shift points toward tissue remodeling and altered metabolite/transporter handling downstream of liver injury.
- Fatty-acid metabolism, amino-acid transport, water-soluble vitamin/cofactor metabolism, and glutathione conjugation are lower in flight, meaning the complementary signal is compatible with depleted or suppressed metabolic handling rather than only lipotoxic activation.

Complementary interpretation:

- The green liver story is mostly "spaceflight disrupts hepatic lipid, CYP/ADME, endocrine, and immune biology." The blue pathway shifts add a mechanistic layer: altered gravity may also be changing hepatic mechanotransduction, cytoskeletal state, immune-cell crosstalk, and tissue remodeling.
- The immune blue pathway shifts are especially useful because they suggest a suppression/composition-change interpretation, not simply inflammation. Several antigen-presentation and TCR-related shifts are FLT_down, consistent with immune dysfunction but not specific enough to call a direct liver replication.
- The ECM, junction, WNT, TGF-beta, PDGF, and MET pathway shifts frame spaceflight liver effects as a repair/remodeling process downstream of metabolic stress. That complements the lipotoxicity literature without claiming that these specific pathway shifts were already proven in mouse liver flight studies.

Slide phrasing:

> Liver blue pathway shifts suggest that known spaceflight lipid/ADME injury is coupled to mechanotransduction, immune suppression, and tissue-remodeling programs.

## Skin HVG

Prior literature frame:

- Spaceflight skin literature reports barrier disruption, altered keratin/cornification, inflammatory/KRAS signaling, oxidative stress, altered metabolism, dermal atrophy, cutaneous muscle changes, hair-follicle changes, and wound-healing risk.
- Key citation families: `Cope2024`, `Park2024`, `Mao2014`, `Neutelings2015`, `Radstake2022`, `Riwaldt2021`, `CuboMateo2021`, `Piipponen2020`, `Pfisterer2021`, `Afshinnekoo2020`.

Representative blue pathways:

- `R-MMU-3247509_CHROMATIN_MODIFYING_ENZYMES`, `R-MMU-4839726_CHROMATIN_ORGANIZATION`, `R-MMU-2990846_SUMOYLATION`, and `R-MMU-1632852_MACROAUTOPHAGY`.
- `R-MMU-190828_GAP_JUNCTION_TRAFFICKING`, `R-MMU-1500931_CELL_CELL_COMMUNICATION`, `R-MMU-421270_CELL_CELL_JUNCTION_ORGANIZATION`, and `R-MMU-1638074_KERATAN_SULFATE_KERATIN_METABOLISM`.
- `R-MMU-9012999_RHO_GTPASE_CYCLE`, `R-MMU-194315_SIGNALING_BY_RHO_GTPASES`, `R-MMU-195258_RHO_GTPASE_EFFECTORS`, `R-MMU-388396_GPCR_DOWNSTREAM_SIGNALLING`, `R-MMU-9006934_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES`, `R-MMU-195721_SIGNALING_BY_WNT`, and `R-MMU-5358351_SIGNALING_BY_HEDGEHOG`.
- `R-MMU-156590_GLUTATHIONE_CONJUGATION`, `R-MMU-156580_PHASE_II_CONJUGATION_OF_COMPOUNDS`, `R-MMU-428157_SPHINGOLIPID_METABOLISM`, `R-MMU-191273_CHOLESTEROL_BIOSYNTHESIS`, `R-MMU-8978868_FATTY_ACID_METABOLISM`, `R-MMU-70326_GLUCOSE_METABOLISM`, and `R-MMU-71291_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES`.
- `R-MMU-8957322_METABOLISM_OF_STEROIDS`, `R-MMU-2980736_PEPTIDE_HORMONE_METABOLISM`, and `R-MMU-381426_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS`.

Direction-aware presentation wording:

- Chromatin-modifying enzymes, chromatin organization, SUMOylation, and macroautophagy are higher in flight, meaning skin shows a stress-regulatory pathway shift consistent with DNA damage, oxidative stress, and transcriptional adaptation.
- Gap-junction trafficking is lower in flight, meaning the model points toward weaker intercellular communication during barrier stress.
- Glutathione conjugation and phase-II conjugation are lower in flight, meaning detox/antioxidant handling may be impaired or depleted rather than simply activated.
- Sphingolipid metabolism, GPCR downstream signaling, RTK signaling, WNT signaling, Hedgehog signaling, and Rho GTPase effectors are higher in flight, meaning the complementary skin signal supports repair, hair-follicle, vascular/receptor, and cell-migration remodeling.
- Cholesterol biosynthesis, fatty-acid metabolism, amino-acid metabolism/transport, and carbohydrate metabolism are lower in flight, meaning the barrier phenotype may include reduced lipid/metabolic support for repair.

Complementary interpretation:

- The direct skin literature emphasizes barrier damage, inflammation, and repair. Blue pathway shifts add a regulatory and communication layer: chromatin/SUMO/autophagy, gap junctions, GPCR/RTK/WNT/Hedgehog, and local hormone/IGF signaling.
- The metabolic blue pathway shifts complement the barrier story. Sphingolipid, cholesterol, fatty-acid, amino-acid, glucose, glutathione, and phase-II conjugation shifts suggest that skin barrier dysfunction may be coupled to altered detox capacity and repair energetics.
- The Rho/GPCR/RTK/WNT/Hedgehog pathway shifts provide a wound-healing and hair-follicle/regeneration perspective. These are plausible in flight-exposed skin but are not direct direction-specific replications.

Slide phrasing:

> Skin blue pathway shifts extend the barrier/inflammation story toward wound-healing coordination: chromatin stress memory, detox and barrier-lipid metabolism, cell-cell communication, and repair signaling.

## Soleus HVG

Prior literature frame:

- Soleus and skeletal muscle spaceflight literature strongly supports atrophy, fiber-type and contractile remodeling, mitochondrial/metabolic changes, proteome remodeling, collagen/ECM changes, immune involvement, and altered mechanotransduction.
- Key citation families: `Gambara2017`, `Sandona2012`, `Okada2021`, `Tascher2017`, `Murgia2024`, `Fitts2010`, `Chopard2009`, `Bonaldo2013`, `Miller2001`, `Tidball2002`, `Dumont2007`.

Representative blue pathways:

- `R-MMU-1430728_METABOLISM`, `R-MMU-211859_BIOLOGICAL_OXIDATIONS`, `R-MMU-382551_TRANSPORT_OF_SMALL_MOLECULES`, `R-MMU-425407_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT`, and `R-MMU-372790_SIGNALING_BY_GPCR`.
- `R-MMU-1474228_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX`, `R-MMU-1442490_COLLAGEN_DEGRADATION`, `R-MMU-1474290_COLLAGEN_FORMATION`, `R-MMU-1650814_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES`, `R-MMU-216083_INTEGRIN_CELL_SURFACE_INTERACTIONS`, and `R-MMU-421270_CELL_CELL_JUNCTION_ORGANIZATION`.
- `R-MMU-9012999_RHO_GTPASE_CYCLE`, `R-MMU-194315_SIGNALING_BY_RHO_GTPASES`, `R-MMU-9013149_RAC1_GTPASE_CYCLE`, `R-MMU-195258_RHO_GTPASE_EFFECTORS`, and `R-MMU-9716542_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3`.
- `R-MMU-1640170_CELL_CYCLE`, `R-MMU-1266738_DEVELOPMENTAL_BIOLOGY`, `R-MMU-73894_DNA_REPAIR`, `R-MMU-5684996_MAPK1_MAPK3_SIGNALING`, `R-MMU-168898_TOLL_LIKE_RECEPTOR_CASCADES`, `R-MMU-5668541_TNFR2_NON_CANONICAL_NF_KB_PATHWAY`, and `R-MMU-449147_SIGNALING_BY_INTERLEUKINS`.
- `R-MMU-5576891_CARDIAC_CONDUCTION` is blue only in an off-label sense: the Reactome label is cardiac, but overlapping ion-channel and excitation-coupling genes can be relevant to skeletal muscle.

Direction-aware presentation wording:

- Broad metabolism, biological oxidations, small-molecule transport, and SLC-mediated transport are lower in flight, meaning soleus shows reduced oxidative/metabolic and transporter pathway capacity, consistent with disuse and atrophy biology.
- Collagen degradation and extracellular-matrix degradation are lower in flight, while collagen formation and collagen-modifying enzymes are higher or mixed, meaning the pathway shift suggests matrix turnover/remodeling rather than a one-direction collagen story.
- Rho GTPase signaling, RAC1 cycle, cell-cell junction organization, MAPK1/3 signaling, and TNFR2/NF-kB are higher in flight, meaning the model detects mechanotransduction, adhesion, and inflammatory-remodeling signals on top of the canonical atrophy phenotype.
- Cell cycle and developmental biology are higher in flight, meaning the complementary signal may represent myogenic repair/reprogramming or nonmyofiber cell activity rather than mature muscle-fiber growth.
- Interleukin signaling and broad RTK signaling are lower in flight, meaning some anabolic or immune-activation pathways may be suppressed even while other inflammatory/remodeling pathways rise.

Complementary interpretation:

- The direct soleus literature supports atrophy and metabolic suppression. Blue pathway shifts add the systems context: atrophy is coupled to transport, ECM/collagen remodeling, integrin/Rho mechanotransduction, immune signaling, and possible myogenic/developmental reprogramming.
- The metabolism and biological-oxidation blue shifts are not novel by themselves, but they connect expiMap pathway shifts to a broader loss of oxidative and transport capacity.
- ECM/collagen/integrin/junction pathway shifts are presentation-relevant because they shift the framing from "muscle fibers shrink" to "the muscle niche and load-sensing matrix remodel."
- Rho/RAC/MIRO and MAPK pathway shifts provide a plausible mechanotransduction and mitochondrial-cytoskeletal bridge, but direct pathway-level direction in soleus flight studies is limited.

Slide phrasing:

> Soleus blue pathway shifts suggest that spaceflight atrophy is coupled to metabolic transport loss, ECM/integrin mechanotransduction, immune remodeling, and partial myogenic repair signals.

## Thymus HVG

Prior literature frame:

- Thymus and immune spaceflight literature supports thymic involution, reduced thymopoiesis, loss of T-cell precursors, altered thymic epithelial/stromal organization, immune dysregulation, mitochondrial stress, apoptosis/stress signaling, and broad epigenetic effects.
- Key citation families: `Gridley2013`, `Horie2019`, `Benjamin2016`, `Woods2003`, `Akiyama2020`, `Crucian2018`, `Muramatsu2025`, `Han2023`, `Novoselova2015`, `DaSilveira2020`, `Afshinnekoo2020`, `Okamura2024`, `Winer2025`.

Representative blue pathways:

- `R-MMU-8980692_RHOA_GTPASE_CYCLE`, `R-MMU-9012999_RHO_GTPASE_CYCLE`, `R-MMU-9013149_RAC1_GTPASE_CYCLE`, `R-MMU-195258_RHO_GTPASE_EFFECTORS`, and `R-MMU-194315_SIGNALING_BY_RHO_GTPASES`.
- `R-MMU-198933_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL`, `R-MMU-1500931_CELL_CELL_COMMUNICATION`, `R-MMU-421270_CELL_CELL_JUNCTION_ORGANIZATION`, `R-MMU-1474244_EXTRACELLULAR_MATRIX_ORGANIZATION`, `R-MMU-1474228_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX`, `R-MMU-1474290_COLLAGEN_FORMATION`, `R-MMU-170834_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX`, `R-MMU-195721_SIGNALING_BY_WNT`, and `R-MMU-201681_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT`.
- `R-MMU-1169410_ANTIMICROBIAL_MECHANISM_OF_IFN_STIMULATED_GENES`, `R-MMU-6798695_NEUTROPHIL_DEGRANULATION`, `R-MMU-977606_REGULATION_OF_COMPLEMENT_CASCADE`, `R-MMU-5668541_TNFR2_NON_CANONICAL_NF_KB_PATHWAY`, and `R-MMU-2142753_ARACHIDONATE_METABOLISM`.
- `R-MMU-3247509_CHROMATIN_MODIFYING_ENZYMES`, `R-MMU-4839726_CHROMATIN_ORGANIZATION`, `R-MMU-9932444_ATP_DEPENDENT_CHROMATIN_REMODELERS`, `R-MMU-1428517_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT`, `R-MMU-5684996_MAPK1_MAPK3_SIGNALING`, `R-MMU-5673001_RAF_MAP_KINASE_CASCADE`, and `R-MMU-1257604_PIP3_ACTIVATES_AKT_SIGNALING`.

Direction-aware presentation wording:

- RHOA/Rho/RAC GTPase cycle and cell-cell junction organization are lower in flight, meaning thymus shows a pathway shift consistent with impaired thymocyte/stromal motility, adhesion, and microenvironment organization.
- Lymphoid/non-lymphoid immunoregulatory interactions and broad cell-cell communication are lower in flight, meaning the complementary signal supports disrupted thymocyte-stromal crosstalk during involution.
- ECM organization, ECM degradation, collagen formation, TGF-beta receptor signaling, and TCF/WNT signaling are higher in flight, meaning the model points toward stromal niche remodeling rather than only thymocyte loss.
- Immune system, IFN-stimulated antimicrobial mechanisms, neutrophil degranulation, TNFR2/NF-kB, and arachidonate metabolism are higher in flight, meaning thymus has an innate/inflammatory stress component alongside impaired adaptive thymopoiesis.
- Aerobic respiration / respiratory electron transport, MAPK1/3, RAF-MAPK, and ATP-dependent chromatin remodelers are lower in flight, meaning the pathway shift is compatible with mitochondrial dysfunction and reduced growth/activation-state signaling.
- Chromatin organization and PIP3-AKT signaling are higher in flight, meaning some regulatory and survival/stress-response programs may be compensatory rather than simply suppressed.

Complementary interpretation:

- The direct thymus literature supports involution and impaired thymopoiesis. Blue pathway shifts suggest the mechanism is not just cell loss: the thymic niche may be remodeling through stromal ECM, junctions, thymocyte-stromal communication, and cytoskeletal/migration pathways.
- Rho/RAC/RHOA and junction/ECM pathway shifts add a physical and microenvironmental interpretation: altered gravity may affect thymocyte migration, adhesion, and epithelial/stromal architecture.
- IFN-stimulated gene, neutrophil degranulation, complement, TNFR2/NF-kB, and arachidonate pathway shifts add an innate/antiviral inflammatory layer that complements adaptive immune suppression.
- Chromatin, PTM, MAPK/RAF, PI3K/AKT, and respiratory-chain pathway shifts add regulatory and metabolic stress context, consistent with spaceflight multi-omics but not thymus-specific enough to call green.

Slide phrasing:

> Thymus blue pathway shifts reframe thymic involution as niche remodeling: impaired thymocyte-stromal communication, cytoskeletal migration changes, innate immune stress, and epigenetic/metabolic reprogramming.

## Cross-Tissue Takeaway

The blue pathway shifts are useful for presentation because they are not the headline replicated findings. They show what expiMap adds beyond standard DGE/pathway enrichment:

- Liver: metabolic injury plus mechanotransduction, immune suppression, and remodeling.
- Skin: barrier/inflammation plus chromatin stress memory, detox capacity, wound-healing signaling, and cell communication.
- Soleus: atrophy plus ECM/integrin mechanotransduction, transport loss, immune remodeling, and myogenic repair.
- Thymus: involution plus stromal niche remodeling, migration/adhesion, innate immune stress, and epigenetic/metabolic regulation.
