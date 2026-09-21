from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tests.support import load_script_module, write_text


class QualStatCsvTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qs = load_script_module("qualstat.py", "qualstat_csv_test")

    def config_for(self, path: Path, analysis_type: str = "abfe"):
        return self.qs.Config(
            input_file=path,
            output_file=path.with_name("report.log"),
            analysis_type=analysis_type,
            energy_unit="kJ/mol",
            bootstrap_rounds=10,
            random_seed=1,
            significance_multiplier=1.645,
            statistics=("MAD",),
            cycle_analysis=False,
        )

    def read(self, directory: Path, text: str, analysis_type: str = "abfe"):
        path = directory / "data.csv"
        write_text(path, text)
        return self.qs.read_dataset(self.config_for(path, analysis_type))

    def test_valid_abfe_rows_and_zero_uncertainty_are_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            data = self.read(
                Path(directory),
                "ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty\n"
                "A,1.0,0.0,1.2,0.0\n",
            )

        self.assertEqual(data.ligand_a, ("A",))
        self.assertIsNone(data.ligand_b)
        self.assertEqual(data.calculated_uncertainty, (0.0,))
        self.assertEqual(data.experimental_uncertainty, (0.0,))

    def test_headers_must_match_exactly(self):
        headers = (
            "calculated,ligand,calculated_uncertainty,experimental,experimental_uncertainty",
            "ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty_typo",
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for header in headers:
                with self.subTest(header=header):
                    with self.assertRaisesRegex(ValueError, "CSV header must be exactly"):
                        self.read(directory, header + "\nA,1,0,1,0\n")

    def test_empty_and_malformed_rows_are_rejected(self):
        cases = (
            ("", "CSV header must be exactly"),
            (
                "ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty\n",
                "CSV contains no data records",
            ),
            (
                "ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty\n"
                "A,1,0,1\n",
                "expected 5 columns, got 4",
            ),
            (
                "ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty\n"
                "A,not-a-number,0,1,0\n",
                "calculated must be numeric",
            ),
            (
                "ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty\n"
                "A,nan,0,1,0\n",
                "calculated must be finite",
            ),
            (
                "ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty\n"
                "A,inf,0,1,0\n",
                "calculated must be finite",
            ),
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for text, message in cases:
                with self.subTest(message=message):
                    with self.assertRaisesRegex(ValueError, message):
                        self.read(directory, text)

    def test_blank_ligand_and_negative_uncertainties_are_rejected(self):
        header = "ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty\n"
        cases = (
            (" ,1,0,1,0\n", "ligand must not be blank"),
            ("A,1,-0.1,1,0\n", "calculated_uncertainty must be non-negative"),
            ("A,1,0,1,-0.1\n", "experimental_uncertainty must be non-negative"),
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for row, message in cases:
                with self.subTest(message=message):
                    with self.assertRaisesRegex(ValueError, message):
                        self.read(directory, header + row)

    def test_rbfe_blank_self_and_duplicate_pairs_are_rejected(self):
        header = "ligand_a,ligand_b,calculated,calculated_uncertainty,experimental,experimental_uncertainty\n"
        cases = (
            ("A, ,1,0,1,0\n", "ligand_b must not be blank"),
            ("A,A,1,0,1,0\n", "RBFE self-edge A -> A is invalid"),
            (
                "A,B,1,0,1,0\nA,B,2,0,2,0\n",
                "duplicate RBFE pair A, B",
            ),
            (
                "A,B,1,0,1,0\nB,A,-1,0,-1,0\n",
                "duplicate RBFE pair B, A",
            ),
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for rows, message in cases:
                with self.subTest(message=message):
                    with self.assertRaisesRegex(ValueError, message):
                        self.read(directory, header + rows, "rbfe")


if __name__ == "__main__":
    unittest.main()
