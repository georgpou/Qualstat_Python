from __future__ import annotations

import math
import unittest

import pandas as pd

from tests.support import load_script_module


try:
    from cinnabar.femap import FEMap  # noqa: F401
    from openff.units import unit  # noqa: F401
except ImportError:
    HAS_CINNABAR = False
else:
    HAS_CINNABAR = True


class RbfeToAbfeReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rbfe = load_script_module("rbfe_to_abfe.py", "rbfe_to_abfe_reconstruction_test")

    def network(self):
        # Hidden ligand energies are L1=0, L2=2, L3=-1, and L4=3.
        return pd.DataFrame(
            [
                ["L1", "L2", 2.0, 0.20],
                ["L2", "L3", -3.0, 0.30],
                ["L3", "L1", 1.0, 0.25],
                ["L2", "L4", 1.0, 0.40],
                ["L4", "L1", -3.0, 0.35],
                ["L3", "L4", 4.0, 0.30],
            ],
            columns=["ligand_A", "ligand_B", "DDG_kJ_mol", "DDG_uncertainty_kJ_mol"],
        )

    def experimental(self, offset: float = -20.0):
        hidden = {"L1": 0.0, "L2": 2.0, "L3": -1.0, "L4": 3.0}
        return pd.DataFrame(
            [[ligand, value + offset] for ligand, value in hidden.items()],
            columns=["ligand", "DG_exp_kJ_mol"],
        )

    def test_consistent_network_reproduces_pairwise_differences_and_alignment(self):
        network = self.network()
        result = self.rbfe.analyze(network, self.experimental())
        values = result.set_index("ligand")["DG_calc_kJ_mol"].to_dict()

        self.assertEqual(list(result["ligand"]), ["L1", "L2", "L3", "L4"])
        self.assertEqual(
            list(result.columns),
            [
                "ligand",
                "DG_exp_kJ_mol",
                "DG_calc_kJ_mol",
                "DG_calc_uncertainty_kJ_mol",
            ],
        )
        for row in network.itertuples(index=False):
            self.assertAlmostEqual(values[row.ligand_B] - values[row.ligand_A], row.DDG_kJ_mol, places=10)

        self.assertAlmostEqual(
            result["DG_calc_kJ_mol"].mean(),
            result["DG_exp_kJ_mol"].mean(),
            places=10,
        )
        self.assertTrue(result["DG_calc_uncertainty_kJ_mol"].map(math.isfinite).all())
        self.assertTrue((result["DG_calc_uncertainty_kJ_mol"] > 0.0).all())

    def test_explicit_numpy_backend_is_reported(self):
        result = self.rbfe.analyze(
            self.network(), self.experimental(), backend="numpy"
        )

        self.assertEqual(result.attrs["estimator"], "NumPy weighted least squares")

    @unittest.skipUnless(HAS_CINNABAR, "OpenFreeEnergy Cinnabar is not installed")
    def test_explicit_cinnabar_backend_is_exercised_when_installed(self):
        result = self.rbfe.analyze(
            self.network(), self.experimental(), backend="cinnabar"
        )

        self.assertEqual(result.attrs["estimator"], "OpenFreeEnergy Cinnabar")

    def test_positive_a_to_b_edge_has_positive_reconstructed_difference(self):
        result = self.rbfe.analyze(self.network(), self.experimental())
        values = result.set_index("ligand")["DG_calc_kJ_mol"]

        # The input convention is DDG(A -> B) = G(B) - G(A).
        self.assertGreater(values["L2"] - values["L1"], 0.0)
        self.assertAlmostEqual(values["L2"] - values["L1"], 2.0, places=10)

    def test_common_experimental_offset_shifts_only_the_aligned_values(self):
        network = self.network()
        first = self.rbfe.analyze(network, self.experimental(-20.0))
        shifted = self.rbfe.analyze(network, self.experimental(-10.0))

        for column in ("DG_calc_kJ_mol",):
            for left, right in zip(first[column], shifted[column]):
                self.assertAlmostEqual(right - left, 10.0, places=10)
        for left, right in zip(
            first["DG_calc_uncertainty_kJ_mol"],
            shifted["DG_calc_uncertainty_kJ_mol"],
        ):
            self.assertAlmostEqual(left, right, places=12)

        first_values = first.set_index("ligand")["DG_calc_kJ_mol"]
        shifted_values = shifted.set_index("ligand")["DG_calc_kJ_mol"]
        for row in network.itertuples(index=False):
            first_difference = first_values[row.ligand_B] - first_values[row.ligand_A]
            shifted_difference = shifted_values[row.ligand_B] - shifted_values[row.ligand_A]
            self.assertAlmostEqual(first_difference, shifted_difference, places=12)

    def test_reversing_edge_storage_direction_does_not_change_the_fit(self):
        network = self.network()
        reversed_network = network.copy()
        rows_to_reverse = [0, 3, 5]
        for index in rows_to_reverse:
            ligand_a = reversed_network.loc[index, "ligand_A"]
            reversed_network.loc[index, "ligand_A"] = reversed_network.loc[index, "ligand_B"]
            reversed_network.loc[index, "ligand_B"] = ligand_a
            reversed_network.loc[index, "DDG_kJ_mol"] *= -1.0

        original = self.rbfe.analyze(network, self.experimental())
        reversed_result = self.rbfe.analyze(reversed_network, self.experimental())

        pd.testing.assert_frame_equal(original, reversed_result, atol=1e-12, rtol=1e-12)

    def test_scaling_all_edge_uncertainties_scales_only_output_uncertainties(self):
        network = self.network()
        scaled_network = network.copy()
        scaled_network["DDG_uncertainty_kJ_mol"] *= 2.0

        original = self.rbfe.analyze(network, self.experimental())
        scaled = self.rbfe.analyze(scaled_network, self.experimental())

        self.assertTrue(
            (original["DG_calc_kJ_mol"] - scaled["DG_calc_kJ_mol"])
            .abs()
            .lt(1e-12)
            .all()
        )
        for first, second in zip(
            original["DG_calc_uncertainty_kJ_mol"],
            scaled["DG_calc_uncertainty_kJ_mol"],
        ):
            self.assertAlmostEqual(second, 2.0 * first, places=12)


if __name__ == "__main__":
    unittest.main()
