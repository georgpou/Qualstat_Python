from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import yaml

from tests.support import REPO_ROOT, run_python


class QualStatCliTests(unittest.TestCase):
    def test_all_help_aliases_show_schemas_metrics_and_descriptions(self):
        aliases = ("-h", "-help", "--help")
        metrics = (
            "MAD", "MADtr", "r2", "PI", "RMSD", "MSD", "Median", "MQ", "Q",
            "slope", "inter", "multi", "AbsMed", "tau", "regMAD", "ROCar",
            "taux", "taur", "taurx", "r22", "slope2", "max", "R", "rho", "rho2",
        )

        for alias in aliases:
            with self.subTest(alias=alias):
                completed = run_python(["scripts/qualstat.py", alias])
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertIn("ABFE: ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty", completed.stdout)
                self.assertIn("RBFE: ligand_a,ligand_b,calculated,calculated_uncertainty,experimental,experimental_uncertainty", completed.stdout)
                for metric in metrics:
                    self.assertIn(metric, completed.stdout)
                self.assertIn("Mean absolute deviation.", completed.stdout)

    def test_template_generation_has_documented_defaults_and_all_switches(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for args, expected_path in (
                ([str(REPO_ROOT / "scripts/qualstat.py"), "--write-template"], directory / "qualstat_template.yaml"),
                ([str(REPO_ROOT / "scripts/qualstat.py"), "--write-template", "custom.yaml"], directory / "custom.yaml"),
            ):
                with self.subTest(expected_path=expected_path):
                    completed = run_python(args, cwd=directory)
                    self.assertEqual(completed.returncode, 0, completed.stderr)
                    self.assertTrue(expected_path.exists())
                    document = yaml.safe_load(expected_path.read_text(encoding="utf-8"))
                    self.assertEqual(document["bootstrap_rounds"], 1000)
                    self.assertTrue(document["output_file"].endswith(".log"))
                    self.assertTrue(document["cycle_analysis"])
                    self.assertEqual(
                        set(document["statistics"]),
                        {
                            "MAD", "MADtr", "r2", "PI", "RMSD", "MSD", "Median", "MQ", "Q",
                            "slope", "inter", "multi", "AbsMed", "tau", "regMAD", "ROCar",
                            "taux", "taur", "taurx", "r22", "slope2", "max", "R", "rho", "rho2",
                        },
                    )
                    self.assertTrue(all(type(value) is bool for value in document["statistics"].values()))
                    self.assertIn("# RBFE only: enumerate unique simple cycles", expected_path.read_text(encoding="utf-8"))
                    self.assertEqual(
                        {
                            name
                            for name, enabled in document["statistics"].items()
                            if enabled
                        },
                        {"MAD", "RMSD", "taurx", "r22"},
                    )
                    self.assertIn("s / sqrt(n)", expected_path.read_text(encoding="utf-8"))

    def test_template_generation_requires_force_to_replace_an_existing_file(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            target = directory / "custom.yaml"
            target.write_text("SENTINEL\n", encoding="utf-8")
            command = [
                str(REPO_ROOT / "scripts/qualstat.py"),
                "--write-template",
                str(target),
            ]

            refused = run_python(command, cwd=directory)
            self.assertEqual(refused.returncode, 2)
            self.assertIn("already exists", refused.stderr)
            self.assertIn("--force", refused.stderr)
            self.assertEqual(target.read_text(encoding="utf-8"), "SENTINEL\n")

            replaced = run_python([*command, "--force"], cwd=directory)
            self.assertEqual(replaced.returncode, 0, replaced.stderr)
            self.assertIn("input_file: results.csv", target.read_text(encoding="utf-8"))

    def test_invalid_configuration_returns_exit_code_two(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            settings = directory / "invalid.yaml"
            settings.write_text("analysis_type: abfe\n", encoding="utf-8")

            completed = run_python(["scripts/qualstat.py", str(settings)])

        self.assertEqual(completed.returncode, 2)
        self.assertIn("Error:", completed.stderr)
        self.assertIn("missing required YAML setting", completed.stderr)


if __name__ == "__main__":
    unittest.main()
