import ast
import hashlib
import json
from pathlib import Path
import re
import sys
import tomllib
import unittest
from unittest import mock
import zipfile

import pandas as pd
import yaml

from nasa_mouse_diffusion.paper_parity import (
    annotate_importance_literature,
    build_synthetic_guided_paper,
)
from nasa_mouse_expimap import integrate_reassessed_tissues_paper


ROOT = Path(__file__).resolve().parents[1]


class HandoffIntegrityTests(unittest.TestCase):
    def test_single_requirements_entrypoint(self):
        environment = yaml.safe_load(
            (ROOT / "environment.yml").read_text(encoding="utf-8")
        )
        pip_dependencies = next(
            item["pip"]
            for item in environment["dependencies"]
            if isinstance(item, dict) and "pip" in item
        )

        self.assertEqual(pip_dependencies, ["-r requirements.txt"])
        self.assertTrue((ROOT / "requirements.txt").is_file())
        self.assertFalse((ROOT / "requirements-nasa-mouse-glare.txt").exists())
        self.assertFalse((ROOT / "requirements-nasa-mouse-generative.txt").exists())

    def test_archs4_candidate_dumps_are_not_versioned_data(self):
        candidate_dumps = list(
            (ROOT / "data/archs4").glob("archs4_mouse_*_candidate_samples.tsv")
        )
        self.assertEqual(candidate_dumps, [])

    def test_public_reference_recovery_is_documented(self):
        artifact_text = (ROOT / "ARTIFACTS.md").read_text(encoding="utf-8")
        handoff_text = (ROOT / "MENTOR_HANDOFF.md").read_text(encoding="utf-8")
        for text in (artifact_text, handoff_text):
            self.assertIn("prepare-references", text)
            self.assertIn("public", text.lower())
        self.assertIn("mouse_gene_v2.5.h5", artifact_text)
        self.assertIn("be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad", artifact_text)

    def test_project_license_metadata_agrees(self):
        license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
        self.assertTrue(license_text.startswith("MIT License\n"))
        self.assertIn("Copyright (c) 2026 Jason Trinh", license_text)

        with (ROOT / "pyproject.toml").open("rb") as handle:
            project = tomllib.load(handle)["project"]
        citation = yaml.safe_load(
            (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        )

        self.assertEqual(project["license"], "MIT")
        self.assertIn("LICENSE", project["license-files"])
        self.assertIn("THIRD_PARTY_NOTICES.md", project["license-files"])
        self.assertEqual(citation["license"], "MIT")

    def test_internship_report_manifest_inputs_exist(self):
        manifest = pd.read_csv(
            ROOT / "paper/slstp_internship_report/source_data/source_manifest.tsv",
            sep="\t",
        )
        missing = [path for path in manifest["path"] if not (ROOT / path).exists()]
        self.assertEqual(missing, [])

    def test_frozen_figure_sources_are_complete(self):
        tables = build_synthetic_guided_paper.load_frozen_figure_tables()
        expected = {
            "arch_summary",
            "locked_repeats",
            "naive_utility",
            "thymus_core",
            "thymus_reactome",
            "muscle_summary",
            "soleus_genes",
            "muscle_reactome",
            "tissue_summary",
            "model_screen",
            "matched_utility",
            "matched_candidates",
            "matched_consensus_comparison",
            "development_highlights",
        }
        self.assertEqual(set(tables), expected)
        self.assertEqual(len(tables["development_highlights"]), 10)

        missing_root = ROOT / "does-not-exist"
        with mock.patch.object(
            integrate_reassessed_tissues_paper,
            "REASSESSMENT_DIR",
            missing_root,
        ):
            effects = integrate_reassessed_tissues_paper.project_effects(
                "spleen",
                ("R-MMU-202403_TCR_SIGNALING",),
            )
        self.assertFalse(effects.empty)
        self.assertEqual(effects["heldout_project"].nunique(), 5)

    def test_final_generative_tissue_evidence_matches_manuscript(self):
        evidence = pd.read_csv(
            ROOT
            / "paper/synthetic_guided_spaceflight/source_data/"
            "table_6_tissue_evidence.tsv",
            sep="\t",
        )
        self.assertEqual(
            evidence.columns.tolist(),
            [
                "tissue",
                "matched_all_gene_result",
                "secondary_consensus_result",
                "interpretation",
            ],
        )
        self.assertEqual(len(evidence), 8)
        self.assertEqual(evidence.iloc[0]["tissue"], "Thymus")
        self.assertIn("15 genes", evidence.iloc[0]["matched_all_gene_result"])
        self.assertIn(
            "Four shared-importance",
            evidence.iloc[1]["matched_all_gene_result"],
        )
        self.assertNotIn("guided delta", "\n".join(evidence.astype(str).stack()))

    def test_selected_feature_comparison_bundle_is_complete(self):
        bundle = ROOT / "outputs/comparison/selected_features"
        pathways = pd.read_csv(bundle / "pathway_crosswalk.tsv", sep="\t")
        genes = pd.read_csv(bundle / "gene_crosswalk.tsv", sep="\t")
        stable = pd.read_csv(
            bundle / "generative_selected_arm_stable_features.tsv",
            sep="\t",
        )
        all_arm_stable = pd.read_csv(
            bundle / "generative_all_arm_stable_features.tsv.gz",
            sep="\t",
        )
        coverage = pd.read_csv(bundle / "generative_analysis_coverage.tsv", sep="\t")
        matched = pd.read_csv(bundle / "generative_matched_genes.tsv", sep="\t")
        consensus = pd.read_csv(bundle / "generative_consensus_genes.tsv", sep="\t")

        self.assertEqual(len(pathways), 26)
        self.assertEqual((pathways["method"] == "expiMap").sum(), 16)
        self.assertEqual(
            (pathways["method"] == "conditional_DDIM_classifier").sum(),
            10,
        )
        self.assertTrue(pathways["pathway_id"].str.fullmatch(r"R-MMU-\d+").all())
        self.assertTrue(
            pathways.apply(
                lambda row: row["pathway_term"].startswith(row["pathway_id"]),
                axis=1,
            ).all()
        )
        self.assertEqual(len(stable), 679)
        self.assertEqual(len(all_arm_stable), 3262)
        self.assertEqual(
            len(
                all_arm_stable[
                    ["analysis_scope", "tissue", "gene"]
                ].drop_duplicates()
            ),
            1307,
        )
        self.assertEqual(
            len(all_arm_stable[["analysis_scope", "tissue"]].drop_duplicates()),
            27,
        )
        self.assertEqual(len(coverage), 27)
        self.assertEqual(
            coverage["selected_arm_feature_comparison_available"].sum(),
            22,
        )
        unavailable = coverage.loc[
            ~coverage["selected_arm_feature_comparison_available"],
            ["scope", "tissue"],
        ]
        self.assertEqual(
            set(map(tuple, unavailable.itertuples(index=False, name=None))),
            {
                ("tissue", "cecum"),
                ("tissue", "colon"),
                ("tissue", "liver"),
                ("muscle_group", "edl"),
                ("muscle_group", "quadriceps"),
            },
        )
        self.assertEqual(
            stable["minimum_selection_frequency"].unique().tolist(),
            [0.5],
        )
        self.assertEqual(
            stable["minimum_coefficient_sign_agreement"].unique().tolist(),
            [0.75],
        )
        self.assertEqual(len(matched), 21)
        self.assertEqual(len(consensus), 49)
        self.assertEqual(coverage["full_selected_arm_feature_count"].sum(), 4475)
        self.assertEqual(coverage["all_arm_stable_feature_row_count"].sum(), 3262)
        self.assertEqual(coverage["all_arm_stable_unique_gene_count"].sum(), 1307)
        self.assertEqual(coverage["stable_selected_feature_count"].sum(), 679)
        self.assertEqual(coverage["matched_primary_gene_count"].sum(), 21)
        self.assertEqual(coverage["consensus_secondary_gene_count"].sum(), 49)
        self.assertEqual(coverage["grouped_pathway_count"].sum(), 10)
        self.assertEqual(genes["expimap_pathway_member"].sum(), 743)
        self.assertEqual(genes["generative_any_arm_stable_feature"].sum(), 1307)
        self.assertEqual(
            genes["generative_selected_arm_stable_feature"].sum(),
            679,
        )
        self.assertEqual(genes["matched_primary"].sum(), 21)
        self.assertEqual(genes["consensus_secondary"].sum(), 49)
        expimap_symbols = genes.loc[
            genes["expimap_pathway_member"],
            "gene_symbol",
        ]
        self.assertFalse(expimap_symbols.str.startswith("ENSMUSG").any())
        expimap_pathway_ids = genes.loc[
            genes["expimap_pathway_member"],
            "retained_pathway_ids",
        ].str.split(";")
        self.assertTrue(
            expimap_pathway_ids.explode().str.fullmatch(r"R-MMU-\d+").all()
        )
        self.assertFalse(
            genes.duplicated(["analysis_scope", "tissue", "gene_id"]).any()
        )

        workbook = pd.ExcelFile(bundle / "selected_feature_comparison.xlsx")
        self.assertEqual(
            set(workbook.sheet_names),
            {
                "guide",
                "gene_crosswalk",
                "pathway_crosswalk",
                "expimap_pathways",
                "expimap_genes",
                "expimap_members",
                "gen_stable_features",
                "gen_all_arm_stable",
                "gen_matched_genes",
                "gen_consensus_genes",
                "gen_grouped_pathways",
                "gen_analysis_coverage",
            },
        )

        manifest = json.loads((bundle / "manifest.json").read_text())
        for filename, metadata in manifest["outputs"].items():
            path = bundle / filename
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                metadata["sha256"],
            )
        self.assertEqual(
            hashlib.sha256(
                (bundle / manifest["workbook"]["path"]).read_bytes()
            ).hexdigest(),
            manifest["workbook"]["sha256"],
        )

        comparison_readme = (bundle / "README.md").read_text(encoding="utf-8")
        for pathway_id in pathways["pathway_id"].unique():
            self.assertIn(pathway_id, comparison_readme)
        for symbol in pd.concat(
            [matched["symbol"], consensus["symbol"]],
            ignore_index=True,
        ).unique():
            self.assertIn(f"`{symbol}`", comparison_readme)
        for filename in (
            "gene_crosswalk.tsv",
            "generative_all_arm_stable_features.tsv.gz",
            "generative_selected_arm_stable_features.tsv",
        ):
            self.assertIn(filename, comparison_readme)

    def test_literature_annotation_sources_resolve(self):
        expimap_dir = ROOT / "paper/asgsr_expimap_hvg/source_data"
        review_dir = expimap_dir / "literature_review"
        expimap_source_keys = set()

        for tissue in ("thymus", "skin", "liver", "soleus"):
            annotations = pd.read_csv(review_dir / "final" / f"{tissue}.tsv", sep="\t")
            source_text = (review_dir / "sources" / f"{tissue}.md").read_text(
                encoding="utf-8"
            )
            tissue_source_keys = set(
                re.findall(r"^- ([^:]+):", source_text, flags=re.MULTILINE)
            )
            expimap_source_keys.update(tissue_source_keys)

            required = [
                "literature_alignment",
                "direction_assessment",
                "review_rationale",
            ]
            self.assertFalse(annotations[required].isna().any(axis=None))
            used_keys = {
                key
                for value in annotations["citations"].fillna("")
                for key in str(value).split(";")
                if key
            }
            self.assertEqual(used_keys - tissue_source_keys, set())

        reassessment_sources = pd.read_csv(
            expimap_dir / "table_s30_kidney_spleen_literature_sources.tsv",
            sep="\t",
        )
        expimap_source_keys.update(reassessment_sources["key"].astype(str))
        retained = pd.read_csv(
            expimap_dir / "table_2_retained_pathway_evidence.tsv", sep="\t"
        )
        self.assertFalse(
            retained[["manual_rationale", "literature_keys"]].isna().any(axis=None)
        )
        retained_keys = {
            key
            for value in retained["literature_keys"]
            for key in str(value).split(";")
            if key
        }
        self.assertEqual(retained_keys - expimap_source_keys, set())

        generative_dir = ROOT / "paper/synthetic_guided_spaceflight/source_data"
        annotation_sets = [
            (
                "consensus genes",
                "table_s16_promoted_gene_literature_annotations.tsv",
                "table_s17_promoted_gene_literature_sources.tsv",
                49,
            ),
            (
                "matched genes",
                "table_s22_matched_gene_literature_annotations.tsv",
                "table_s24_importance_literature_sources.tsv",
                21,
            ),
            (
                "grouped pathways",
                "table_s23_grouped_pathway_literature_annotations.tsv",
                "table_s24_importance_literature_sources.tsv",
                10,
            ),
        ]
        required = [
            "literature_classification",
            "evidence_scope",
            "evidence_relationship",
            "source_ids",
            "literature_summary",
            "interpretation",
        ]
        for name, annotation_file, source_file, expected_rows in annotation_sets:
            with self.subTest(annotation_set=name):
                annotations = pd.read_csv(generative_dir / annotation_file, sep="\t")
                sources = pd.read_csv(generative_dir / source_file, sep="\t")
                self.assertEqual(len(annotations), expected_rows)
                self.assertFalse(annotations[required].isna().any(axis=None))
                source_ids = set(sources["source_id"].astype(str))
                used_ids = {
                    source_id
                    for value in annotations["source_ids"]
                    for source_id in str(value).split(";")
                    if source_id
                }
                self.assertEqual(used_ids - source_ids, set())

    def test_annotation_prompts_and_frozen_inputs_are_preserved(self):
        prompt_path = ROOT / "docs/annotation_prompts.md"
        prompt_text = prompt_path.read_text(encoding="utf-8")
        provenance_text = (ROOT / "docs/annotation_provenance.md").read_text(
            encoding="utf-8"
        )

        for protocol_id in (
            "expimap-pathway-review-v1",
            "synthetic-feature-review-v2",
        ):
            self.assertIn(protocol_id, prompt_text)
            self.assertIn(protocol_id, provenance_text)
        for label in ("aligning", "complementary", "ambiguous", "unmatched"):
            self.assertIn(f"- {label}:", prompt_text)

        frozen = pd.read_csv(
            annotate_importance_literature.FROZEN_GROUPED_INPUT,
            sep="\t",
        )
        committed = pd.read_csv(
            ROOT
            / "paper/synthetic_guided_spaceflight/source_data/"
            "table_s23_grouped_pathway_literature_annotations.tsv",
            sep="\t",
        )
        self.assertEqual(len(frozen), 10)
        self.assertEqual(
            set(map(tuple, frozen[["tissue", "term"]].itertuples(index=False, name=None))),
            set(
                map(
                    tuple,
                    committed[["tissue", "term"]].itertuples(
                        index=False, name=None
                    ),
                )
            ),
        )

        missing_root = ROOT / "does-not-exist"
        with (
            mock.patch.object(
                annotate_importance_literature,
                "GROUPED_INPUT",
                missing_root / "eligible.tsv.gz",
            ),
            mock.patch.object(
                annotate_importance_literature,
                "NONREDUNDANT_INPUT",
                missing_root / "nonredundant.tsv",
            ),
        ):
            clone_only = annotate_importance_literature._collapse_grouped()
        self.assertEqual(len(clone_only), 10)

    def test_editorial_review_covers_final_manuscripts(self):
        record = (ROOT / "docs/editorial_review.md").read_text(encoding="utf-8")
        self.assertIn("https://github.com/blader/humanizer", record)
        self.assertIn("version 2.9.1", record)
        self.assertIn("no-fabrication", record)

        reviewed_paths = [
            ROOT / "paper/slstp_internship_report/manuscript.md",
            ROOT / "paper/asgsr_expimap_hvg/manuscript.md",
            ROOT / "paper/synthetic_guided_spaceflight/manuscript.md",
        ]
        for path in reviewed_paths:
            with self.subTest(path=path.relative_to(ROOT)):
                text = path.read_text(encoding="utf-8")
                self.assertNotRegex(text, r"[\u2013\u2014]")
                self.assertNotRegex(
                    text,
                    r"\[(?:add|confirm|insert|review|todo|tbd)[^]]*\]",
                )

    def test_presentation_deliverables_are_grouped(self):
        presentation_root = ROOT / "presentation"
        self.assertEqual(
            {path.name for path in presentation_root.iterdir()},
            {"README.md", "poster", "midpoint", "final"},
        )

        deliverables = [
            "poster/asgsr_expimap_poster.pptx",
            "poster/asgsr_expimap_poster.pdf",
            "midpoint/SLSTP_2026_Midpoint_Presentation.pptx",
            "midpoint/SLSTP_2026_Midpoint_Presentation.pdf",
            "final/SLSTP_2026_Generative_Transcriptomics.pptx",
            "final/SLSTP_2026_Generative_Transcriptomics.pdf",
            "final/speaker_notes.md",
        ]
        missing = [
            path
            for path in deliverables
            if not (presentation_root / path).is_file()
        ]
        self.assertEqual(missing, [])

        final_deck = (
            presentation_root / "final/SLSTP_2026_Generative_Transcriptomics.pptx"
        )
        with zipfile.ZipFile(final_deck) as archive:
            names = archive.namelist()
        slides = [
            name
            for name in names
            if re.fullmatch(r"ppt/slides/slide\d+\.xml", name)
        ]
        notes = [
            name
            for name in names
            if re.fullmatch(r"ppt/notesSlides/notesSlide\d+\.xml", name)
        ]
        self.assertEqual(len(slides), 29)
        self.assertEqual(len(notes), len(slides))

    def test_documentation_links_resolve(self):
        missing = []
        link_pattern = re.compile(r"\[[^]]+\]\(([^)]+)\)")

        markdown_paths = [
            ROOT / "README.md",
            ROOT / "MENTOR_HANDOFF.md",
            ROOT / "ARTIFACTS.md",
            ROOT / "REPRODUCIBILITY.md",
            ROOT / "assets/README.md",
            ROOT / "tests/README.md",
            *sorted((ROOT / "docs").rglob("*.md")),
            *sorted((ROOT / "src").rglob("README.md")),
            *sorted((ROOT / "configs").rglob("README.md")),
            *sorted((ROOT / "data").rglob("README.md")),
            *sorted((ROOT / "outputs").rglob("README.md")),
            *sorted((ROOT / "paper").rglob("*.md")),
            *sorted((ROOT / "presentation").rglob("*.md")),
        ]
        for markdown_path in markdown_paths:
            for target in link_pattern.findall(
                markdown_path.read_text(encoding="utf-8")
            ):
                target = target.split("#", 1)[0]
                if not target or "://" in target or target.startswith("mailto:"):
                    continue
                resolved = (markdown_path.parent / target).resolve()
                if not resolved.exists():
                    missing.append(f"{markdown_path.relative_to(ROOT)} -> {target}")

        self.assertEqual(missing, [])

    def test_source_packages_and_modules_are_documented(self):
        source_root = ROOT / "src"
        for package_dir in sorted(source_root.glob("nasa_mouse_*")):
            if not package_dir.is_dir() or not (package_dir / "__init__.py").is_file():
                continue
            with self.subTest(package=package_dir.name):
                self.assertTrue((package_dir / "README.md").is_file())

            for module_path in sorted(package_dir.rglob("*.py")):
                if module_path.name == "__init__.py":
                    continue
                with self.subTest(module=module_path.relative_to(ROOT)):
                    tree = ast.parse(module_path.read_text(encoding="utf-8"))
                    self.assertTrue(ast.get_docstring(tree))

                    readme_dir = module_path.parent
                    while not (readme_dir / "README.md").is_file():
                        self.assertNotEqual(readme_dir, source_root)
                        readme_dir = readme_dir.parent
                    module_name = module_path.relative_to(readme_dir).as_posix()
                    readme_text = (readme_dir / "README.md").read_text(
                        encoding="utf-8"
                    )
                    self.assertIn(f"`{module_name}`", readme_text)

            for script_path in sorted(package_dir.rglob("*.R")):
                with self.subTest(script=script_path.relative_to(ROOT)):
                    readme_text = (package_dir / "README.md").read_text(
                        encoding="utf-8"
                    )
                    self.assertIn(f"`{script_path.name}`", readme_text)

    def test_paper_paths_use_current_repository_layout(self):
        expected_grouped_dir = (
            ROOT
            / "outputs/generative/benchmark/analyses"
            / "grouped_pathway_importance_osdr_disjoint_v1"
        )
        self.assertEqual(
            annotate_importance_literature.GROUPED_DIR, expected_grouped_dir
        )

        legacy_paths = [
            "outputs/generative_benchmark",
            "paper/asgsr_expimap_hvg/poster",
            "presentation/generative_slstp_2026",
            "presentation/SLSTP_2026_Generative_Transcriptomics",
            "presentation/glare",
            "src/nasa_mouse_rna_diffusion",
            "src/expiMap_scarches",
        ]
        paper_text = "\n".join(
            path.read_text(encoding="utf-8")
            for path in sorted((ROOT / "paper").rglob("*.md"))
        )
        for legacy_path in legacy_paths:
            self.assertNotIn(legacy_path, paper_text)

    def test_generative_render_only_skips_analysis_refresh(self):
        arguments = ["build_synthetic_guided_paper", "--render-only"]
        with (
            mock.patch.object(sys, "argv", arguments),
            mock.patch.object(
                build_synthetic_guided_paper, "render_paper_documents"
            ) as render,
            mock.patch.object(
                build_synthetic_guided_paper, "build_source_tables"
            ) as refresh,
        ):
            build_synthetic_guided_paper.main()

        render.assert_called_once_with()
        refresh.assert_not_called()


if __name__ == "__main__":
    unittest.main()
