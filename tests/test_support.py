from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from tests.support import (
    REPO_ROOT,
    load_script_module,
    normalize_input_line,
    run_python,
    write_csv,
)


class SupportTests(unittest.TestCase):
    def test_loads_both_scripts_without_installing_package(self):
        qualstat = load_script_module("qualstat.py", "qualstat_support_test")
        rbfe_to_abfe = load_script_module("rbfe_to_abfe.py", "rbfe_to_abfe_support_test")

        self.assertTrue(hasattr(qualstat, "evaluate_metric"))
        self.assertTrue(hasattr(rbfe_to_abfe, "analyze"))

    def test_write_csv_creates_readable_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nested" / "data.csv"
            write_csv(path, ["ligand", "value"], [["A", "1.25"], ["B", "2.50"]])

            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))

        self.assertEqual(rows, [["ligand", "value"], ["A", "1.25"], ["B", "2.50"]])

    def test_run_python_uses_repository_root_by_default(self):
        completed = run_python(["scripts/qualstat.py", "--help"])

        self.assertEqual(completed.returncode, 0)
        self.assertIn("ABFE:", completed.stdout)
        self.assertEqual(completed.stderr, "")
        self.assertEqual(Path.cwd(), REPO_ROOT)

    def test_normalize_input_line_changes_only_path(self):
        report = "Header\nInput: /machine/specific/path/data.csv\nMetric: 1\n"

        normalized = normalize_input_line(report, "<ABSOLUTE_PATH>")

        self.assertEqual(
            normalized,
            "Header\nInput: <ABSOLUTE_PATH>\nMetric: 1\n",
        )


if __name__ == "__main__":
    unittest.main()
