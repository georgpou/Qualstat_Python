from __future__ import annotations

import math
import unittest
from pathlib import Path

from tests.support import load_script_module


class QualStatCycleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qs = load_script_module("qualstat.py", "qualstat_cycles_test")

    def rbfe_data(self, ligand_a, ligand_b, values, uncertainties):
        n_records = len(values)
        return self.qs.Dataset(
            ligand_a=tuple(ligand_a),
            ligand_b=tuple(ligand_b),
            calculated=tuple(values),
            calculated_uncertainty=tuple(uncertainties),
            experimental=tuple(values),
            experimental_uncertainty=tuple(uncertainties),
        )

    def config(self, analysis_type: str, cycle_analysis: bool):
        return self.qs.Config(
            input_file=Path("data.csv"),
            output_file=Path("report.log"),
            analysis_type=analysis_type,
            energy_unit="kJ/mol",
            bootstrap_rounds=10,
            random_seed=1,
            significance_multiplier=1.645,
            statistics=("MAD",),
            cycle_analysis=cycle_analysis,
        )

    def test_reversed_stored_edge_changes_sign_and_uses_root_sum_square(self):
        data = self.rbfe_data(
            ["A", "C", "C"],
            ["B", "B", "A"],
            [2.0, 3.5, 1.0],
            [0.3, 0.4, 0.5],
        )

        cycles = self.qs.find_cycles(data)

        self.assertEqual(len(cycles), 1)
        self.assertEqual(cycles[0].ligands, ("A", "B", "C"))
        self.assertAlmostEqual(cycles[0].closure, -0.5)
        self.assertAlmostEqual(cycles[0].uncertainty, math.sqrt(0.3**2 + 0.4**2 + 0.5**2))

    def test_cycles_are_unique_and_deterministically_ordered(self):
        data = self.rbfe_data(
            ["A", "B", "C", "A", "D"],
            ["B", "C", "A", "D", "B"],
            [2.0, -3.0, 1.0, 3.0, -1.0],
            [0.1] * 5,
        )

        first = self.qs.find_cycles(data)
        second = self.qs.find_cycles(data)

        self.assertEqual(first, second)
        self.assertEqual([cycle.ligands for cycle in first], sorted({cycle.ligands for cycle in first}))
        self.assertEqual(len({cycle.ligands for cycle in first}), len(first))
        self.assertEqual(
            [cycle.ligands for cycle in first],
            [("A", "B", "C"), ("A", "B", "D"), ("A", "C", "B", "D")],
        )

    def test_no_cycle_report_is_explicit(self):
        data = self.rbfe_data(
            ["A", "B"],
            ["B", "C"],
            [1.0, 2.0],
            [0.1, 0.1],
        )
        report = self.qs.render_report(
            self.config("rbfe", True),
            data,
            (self.qs.MetricResult("MAD", 0.0, 0.0, 10),),
            self.qs.find_cycles(data),
        )

        self.assertIn("RBFE cycle analysis", report)
        self.assertIn("No cycles found.", report)

    def test_cycle_analysis_is_ignored_for_abfe_data(self):
        data = self.qs.Dataset(
            ligand_a=("A", "B"),
            ligand_b=None,
            calculated=(1.0, 2.0),
            calculated_uncertainty=(0.1, 0.1),
            experimental=(1.0, 2.0),
            experimental_uncertainty=(0.1, 0.1),
        )
        report = self.qs.render_report(
            self.config("abfe", True),
            data,
            (self.qs.MetricResult("MAD", 0.0, 0.0, 10),),
            (),
        )

        self.assertIn("Note: cycle_analysis ignored for ABFE data.", report)


if __name__ == "__main__":
    unittest.main()
