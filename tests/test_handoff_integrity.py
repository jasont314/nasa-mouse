import ast
from pathlib import Path
import re
import sys
import tomllib
import unittest
from unittest import mock

import pandas as pd
import yaml

from nasa_mouse_diffusion.paper_parity import (
    annotate_importance_literature,
    build_synthetic_guided_paper,
)


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
        ]
        missing = [
            path
            for path in deliverables
            if not (presentation_root / path).is_file()
        ]
        self.assertEqual(missing, [])

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
            if not package_dir.is_dir():
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
