from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from tests.support import REPO_ROOT, run_python, write_csv


NETWORK_HEADER = ["ligand_A", "ligand_B", "DDG_kJ_mol", "DDG_uncertainty_kJ_mol"]
EXPERIMENTAL_HEADER = ["ligand", "DG_exp_kJ_mol"]


class RbfeToAbfeCliTests(unittest.TestCase):
    def write_inputs(self, directory: Path, rows=None, experimental=None):
        if rows is None:
            rows = [["A", "B", 2.0, 0.3], ["B", "C", -3.0, 0.4], ["C", "A", 1.5, 0.5]]
        if experimental is None:
            experimental = [["A", -10.0], ["B", -8.0], ["C", -11.0]]
        network = write_csv(directory / "network.csv", NETWORK_HEADER, rows)
        experiment = write_csv(directory / "experimental.csv", EXPERIMENTAL_HEADER, experimental)
        return network, experiment

    def run_cli(self, directory: Path, *extra):
        network, experiment = self.write_inputs(directory)
        output = directory / "calculated_abfe.csv"
        completed = run_python(
            [
                str(REPO_ROOT / "scripts" / "rbfe_to_abfe.py"),
                str(network),
                str(experiment),
                "-o",
                str(output),
                *map(str, extra),
            ],
            cwd=directory,
        )
        return completed, output

    def test_valid_cli_serializes_kj_output_without_index(self):
        with tempfile.TemporaryDirectory() as directory:
            completed, output = self.run_cli(Path(directory), "--cycle-closure", 3)
            cycle_output = Path(directory) / "calculated_abfe_cycle_closure.csv"
            edge_output = Path(directory) / "calculated_abfe_cycle_closure_edges.csv"

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("Cycles detected: 1", completed.stdout)
            self.assertIn("Largest normalized closure error", completed.stdout)
            self.assertTrue(output.exists())
            self.assertTrue(cycle_output.exists())
            self.assertTrue(edge_output.exists())

            result = pd.read_csv(output)
            self.assertEqual(
                list(result.columns),
                [
                    "ligand",
                    "DG_exp_kJ_mol",
                    "DG_calc_kJ_mol",
                    "DG_calc_uncertainty_kJ_mol",
                ],
            )
            self.assertNotIn("Unnamed: 0", result.columns)
            self.assertTrue(result.select_dtypes(include="number").notna().all().all())
            self.assertIn("-10.000000", output.read_text(encoding="utf-8"))

    def test_cycle_length_below_three_returns_exit_code_two(self):
        with tempfile.TemporaryDirectory() as directory:
            completed, output = self.run_cli(Path(directory), "--cycle-closure", 2)

        self.assertEqual(completed.returncode, 2)
        self.assertIn("--cycle-closure must be at least 3", completed.stderr)
        self.assertFalse(output.exists())

    def test_late_cycle_failure_preserves_existing_output_and_leaves_no_partial_files(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            rows = [
                ["A", "B", 1.0, 0.2],
                ["B", "A", -1.0, 0.2],
                ["B", "C", 2.0, 0.3],
            ]
            network, experiment = self.write_inputs(
                directory,
                rows=rows,
                experimental=[["A", -10.0], ["B", -9.0], ["C", -7.0]],
            )
            output = directory / "calculated_abfe.csv"
            output.write_text("SENTINEL\n", encoding="utf-8")

            completed = run_python(
                [
                    str(REPO_ROOT / "scripts" / "rbfe_to_abfe.py"),
                    str(network),
                    str(experiment),
                    "-o",
                    str(output),
                    "--cycle-closure",
                    "3",
                ],
                cwd=directory,
            )

            self.assertEqual(completed.returncode, 2)
            self.assertIn("one aggregated edge per unordered ligand pair", completed.stderr)
            self.assertEqual(output.read_text(encoding="utf-8"), "SENTINEL\n")
            self.assertFalse((directory / "calculated_abfe_cycle_closure.csv").exists())
            self.assertFalse((directory / "calculated_abfe_cycle_closure_edges.csv").exists())

    def test_tree_cli_succeeds_and_reports_no_cycles(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            network, experiment = self.write_inputs(
                directory,
                rows=[["A", "B", 2.0, 0.3], ["B", "C", -3.0, 0.4]],
            )
            output = directory / "calculated_abfe.csv"
            completed = run_python(
                [
                    str(REPO_ROOT / "scripts" / "rbfe_to_abfe.py"),
                    str(network),
                    str(experiment),
                    "-o",
                    str(output),
                    "--cycle-closure",
                    "3",
                ],
                cwd=directory,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("No cycles containing 3-3 ligands were found", completed.stdout)
            self.assertTrue((directory / "calculated_abfe_cycle_closure.csv").read_text(encoding="utf-8").startswith("source,cycle"))


if __name__ == "__main__":
    unittest.main()
