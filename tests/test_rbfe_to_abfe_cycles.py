from __future__ import annotations

import math
import unittest

import pandas as pd

from tests.support import load_script_module


class RbfeToAbfeCycleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rbfe = load_script_module("rbfe_to_abfe.py", "rbfe_to_abfe_cycles_test")

    def frame(self, rows):
        return pd.DataFrame(
            rows,
            columns=["ligand_A", "ligand_B", "DDG_kJ_mol", "DDG_uncertainty_kJ_mol"],
        )

    def test_triangle_reports_hand_calculated_closure_and_uncertainty(self):
        network = self.frame(
            [
                ["A", "B", 2.0, 0.3],
                ["B", "C", -3.0, 0.4],
                ["C", "A", 1.5, 0.5],
            ]
        )

        cycles, edges = self.rbfe.cycle_closure_diagnostics(network, max_cycle_length=3)

        self.assertEqual(len(cycles), 1)
        self.assertEqual(tuple(cycles.loc[0, "cycle"]), ("A", "B", "C"))
        self.assertAlmostEqual(cycles.loc[0, "cc_kJ_mol"], 0.5)
        self.assertAlmostEqual(cycles.loc[0, "cc_per_sqrt_edge_kJ_mol"], 0.5 / math.sqrt(3.0))
        self.assertAlmostEqual(cycles.loc[0, "cc_unc_normalized"], 0.5 / math.sqrt(0.3**2 + 0.4**2 + 0.5**2))
        self.assertTrue(edges["mean_cc_per_sqrt_edge_kJ_mol"].map(math.isfinite).all())
        self.assertTrue(edges["max_cc_per_sqrt_edge_kJ_mol"].map(math.isfinite).all())

    def test_complete_four_ligand_graph_reports_all_three_to_four_ligand_cycles(self):
        network = self.frame(
            [
                ["A", "B", 2.0, 0.2],
                ["B", "C", -3.0, 0.3],
                ["C", "A", 1.0, 0.25],
                ["A", "D", 3.0, 0.4],
                ["D", "B", -1.0, 0.35],
                ["C", "D", 4.0, 0.3],
            ]
        )

        cycles, _ = self.rbfe.cycle_closure_diagnostics(network, max_cycle_length=4)

        # K4 has four triangles and three distinct four-node simple cycles.
        self.assertEqual(len(cycles), 7)
        self.assertEqual(len({tuple(cycle) for cycle in cycles["cycle"]}), 7)
        self.assertTrue((cycles["cc_kJ_mol"].map(math.isfinite)).all())
        self.assertTrue((cycles["cc_unc_normalized"].map(math.isfinite)).all())

    def test_tree_returns_empty_reports_with_stable_headers(self):
        network = self.frame(
            [
                ["A", "B", 2.0, 0.3],
                ["B", "C", -3.0, 0.4],
            ]
        )

        cycles, edges = self.rbfe.cycle_closure_diagnostics(network, max_cycle_length=3)

        self.assertTrue(cycles.empty)
        self.assertTrue(edges.empty)
        self.assertEqual(
            list(cycles.columns),
            ["source", "cycle", "cc_unc_normalized", "cc_kJ_mol", "cc_per_sqrt_edge_kJ_mol"],
        )
        self.assertEqual(
            list(edges.columns),
            [
                "source",
                "ligandA",
                "ligandB",
                "n_cycles",
                "mean_cc_per_sqrt_edge_kJ_mol",
                "max_cc_per_sqrt_edge_kJ_mol",
            ],
        )


if __name__ == "__main__":
    unittest.main()
