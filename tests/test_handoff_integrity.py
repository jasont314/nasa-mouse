from pathlib import Path
import re
import sys
import tomllib
import unittest
from unittest import mock

import pandas as pd
import yaml

from nasa_mouse_diffusion.paper_parity import build_synthetic_guided_paper


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

        for markdown_path in sorted((ROOT / "docs").glob("*.md")):
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
