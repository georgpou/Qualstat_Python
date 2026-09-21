from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from tests.support import REPO_ROOT, normalize_input_line, run_python


class QualStatExampleTests(unittest.TestCase):
    def run_copied_example(self, name: str, expected_placeholder: str):
        source = REPO_ROOT / "examples" / name
        with tempfile.TemporaryDirectory() as directory:
            copied = Path(directory) / name
            shutil.copytree(source, copied)
            expected = (copied / "expected.log").read_text(encoding="utf-8")
            completed = run_python(
                [str(REPO_ROOT / "scripts" / "qualstat.py"), str(copied / "settings.yaml")],
                cwd=REPO_ROOT,
            )
            actual = (copied / "expected.log").read_text(encoding="utf-8")
            return completed, expected, actual

    def test_abfe_example_runs_and_matches_reference_log(self):
        completed, expected, actual = self.run_copied_example(
            "qualstat_abfe", "<ABSOLUTE_PATH_TO_ABFE_DATA>"
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Analysis: ABFE", actual)
        self.assertIn("Metric", actual)
        self.assertEqual(
            normalize_input_line(actual, "<ABSOLUTE_PATH_TO_ABFE_DATA>"),
            expected,
        )

    def test_rbfe_example_runs_and_matches_reference_log(self):
        completed, expected, actual = self.run_copied_example(
            "qualstat_rbfe", "<ABSOLUTE_PATH_TO_RBFE_DATA>"
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Analysis: RBFE", actual)
        self.assertIn("RBFE cycle analysis", actual)
        self.assertIn("Cycle #1", actual)
        self.assertEqual(
            normalize_input_line(actual, "<ABSOLUTE_PATH_TO_RBFE_DATA>"),
            expected,
        )

    def test_example_readme_documents_both_qualstat_commands(self):
        readme = (REPO_ROOT / "examples" / "README.md").read_text(encoding="utf-8")

        self.assertIn("python scripts/qualstat.py examples/qualstat_abfe/settings.yaml", readme)
        self.assertIn("python scripts/qualstat.py examples/qualstat_rbfe/settings.yaml", readme)
        self.assertIn("python scripts/rbfe_to_abfe.py", readme)
        self.assertIn("DDG(A -> B) = G(B) - G(A)", readme)
        self.assertIn("not independent absolute-binding-free-energy simulations", readme)
        self.assertIn("Input:", readme)

    def test_project_docs_and_ci_publish_the_canonical_workflow(self):
        readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
        qualstat_docs = (REPO_ROOT / "docs" / "qualstat_docs.md").read_text(encoding="utf-8")
        rbfe_docs = (REPO_ROOT / "docs" / "rbfe_to_abfe_docs.md").read_text(encoding="utf-8")
        environment = (REPO_ROOT / "environment.yml").read_text(encoding="utf-8")
        workflow = (REPO_ROOT / ".github" / "workflows" / "tests.yml").read_text(encoding="utf-8")

        self.assertIn("python -m unittest discover -s tests -v", readme)
        self.assertIn("examples/README.md", readme)
        self.assertIn("conda activate qualstat", readme)
        self.assertIn("name: qualstat", environment)
        self.assertNotIn("qualstat-python", readme + qualstat_docs + rbfe_docs)
        self.assertIn("environment.yml", workflow)
        self.assertIn("python -m unittest discover -s tests -v", workflow)
        self.assertIn("cinnabar=0.6.1", environment)

    def test_rbfe_to_abfe_example_runs_and_matches_committed_tables(self):
        source = REPO_ROOT / "examples" / "rbfe_to_abfe"
        with tempfile.TemporaryDirectory() as directory:
            copied = Path(directory) / "rbfe_to_abfe"
            shutil.copytree(source, copied)
            run_log = (copied / "run.log").read_text(encoding="utf-8")
            expected_main = pd.read_csv(copied / "calculated_abfe.csv")
            expected_cycles = pd.read_csv(copied / "calculated_abfe_cycle_closure.csv")
            expected_edges = pd.read_csv(copied / "calculated_abfe_cycle_closure_edges.csv")
            completed = run_python(
                [
                    str(REPO_ROOT / "scripts" / "rbfe_to_abfe.py"),
                    str(copied / "network.csv"),
                    str(copied / "experimental.csv"),
                    "-o",
                    str(copied / "calculated_abfe.csv"),
                    "--cycle-closure",
                    "4",
                ],
                cwd=REPO_ROOT,
            )
            actual_main = pd.read_csv(copied / "calculated_abfe.csv")
            actual_cycles = pd.read_csv(copied / "calculated_abfe_cycle_closure.csv")
            actual_edges = pd.read_csv(copied / "calculated_abfe_cycle_closure_edges.csv")

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Cycles detected: 7", completed.stdout)
        self.assertIn("Largest normalized closure error: 0.331295", run_log)
        pd.testing.assert_frame_equal(actual_main, expected_main, check_exact=False, rtol=1e-10, atol=1e-10)
        pd.testing.assert_frame_equal(actual_cycles, expected_cycles, check_exact=False, rtol=1e-10, atol=1e-10)
        pd.testing.assert_frame_equal(actual_edges, expected_edges, check_exact=False, rtol=1e-10, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
