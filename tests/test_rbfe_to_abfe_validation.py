from __future__ import annotations

import unittest

import pandas as pd

from tests.support import load_script_module


class RbfeToAbfeValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rbfe = load_script_module("rbfe_to_abfe.py", "rbfe_to_abfe_validation_test")

    def network(self, rows=None, **columns):
        if rows is not None:
            return pd.DataFrame(
                rows,
                columns=[
                    "ligand_A",
                    "ligand_B",
                    "DDG_kJ_mol",
                    "DDG_uncertainty_kJ_mol",
                ],
            )
        return pd.DataFrame(columns)

    def experimental(self, rows=None, **columns):
        if rows is not None:
            return pd.DataFrame(rows, columns=["ligand", "DG_exp_kJ_mol"])
        return pd.DataFrame(columns)

    def valid_network(self):
        return self.network(
            [
                ["A", "B", 2.0, 0.3],
                ["B", "C", -3.0, 0.4],
                ["C", "A", 1.0, 0.5],
            ]
        )

    def valid_experimental(self):
        return self.experimental([["A", -10.0], ["B", -8.0], ["C", -11.0]])

    def assert_inputs_rejected(self, network, experimental, message):
        with self.assertRaisesRegex((ValueError, RuntimeError), message):
            self.rbfe._validate_inputs(network, experimental)

    def test_missing_columns_are_rejected(self):
        self.assert_inputs_rejected(
            pd.DataFrame({"ligand_A": ["A"]}),
            self.valid_experimental(),
            "Network CSV is missing columns",
        )
        self.assert_inputs_rejected(
            self.valid_network(),
            pd.DataFrame({"ligand": ["A"]}),
            "Experimental CSV is missing columns",
        )

    def test_empty_inputs_are_rejected(self):
        self.assert_inputs_rejected(
            self.network(
                ligand_A=pd.Series(dtype=str),
                ligand_B=pd.Series(dtype=str),
                DDG_kJ_mol=pd.Series(dtype=float),
                DDG_uncertainty_kJ_mol=pd.Series(dtype=float),
            ),
            self.valid_experimental(),
            "Network CSV must not be empty",
        )
        self.assert_inputs_rejected(
            self.valid_network(),
            self.experimental(
                ligand=pd.Series(dtype=str),
                DG_exp_kJ_mol=pd.Series(dtype=float),
            ),
            "Experimental CSV must not be empty",
        )

    def test_nonfinite_and_nonnumeric_values_are_rejected(self):
        for bad_value in ("not-a-number", float("nan"), float("inf"), float("-inf")):
            with self.subTest(bad_value=bad_value):
                network = self.valid_network()
                if isinstance(bad_value, str):
                    network["DDG_kJ_mol"] = network["DDG_kJ_mol"].astype(object)
                network.loc[0, "DDG_kJ_mol"] = bad_value
                self.assert_inputs_rejected(network, self.valid_experimental(), "finite numbers")

        experimental = self.valid_experimental()
        experimental["DG_exp_kJ_mol"] = experimental["DG_exp_kJ_mol"].astype(object)
        experimental.loc[0, "DG_exp_kJ_mol"] = "not-a-number"
        self.assert_inputs_rejected(self.valid_network(), experimental, "finite numbers")

    def test_zero_and_negative_rbfe_uncertainties_are_rejected(self):
        # Unlike QualStat, which intentionally accepts zero uncertainties,
        # reconstruction needs strictly positive edge weights for Cinnabar.
        for uncertainty in (0.0, -0.1):
            with self.subTest(uncertainty=uncertainty):
                network = self.valid_network()
                network.loc[0, "DDG_uncertainty_kJ_mol"] = uncertainty
                self.assert_inputs_rejected(
                    network,
                    self.valid_experimental(),
                    "greater than zero",
                )

    def test_self_transformations_are_rejected(self):
        network = self.valid_network()
        network.loc[0, "ligand_B"] = "A"
        self.assert_inputs_rejected(network, self.valid_experimental(), "Self-transformations are not allowed")

    def test_duplicate_experimental_ligands_are_rejected(self):
        experimental = pd.concat([self.valid_experimental(), pd.DataFrame([["A", -9.0]], columns=experimental_columns())], ignore_index=True)
        self.assert_inputs_rejected(self.valid_network(), experimental, "Experimental ligand names must be unique")

    def test_missing_and_extra_experimental_ligands_are_rejected(self):
        missing = self.experimental([["A", -10.0], ["B", -8.0]])
        self.assert_inputs_rejected(self.valid_network(), missing, r"missing experimental=\['C'\]")

        extra = self.experimental([["A", -10.0], ["B", -8.0], ["C", -11.0], ["D", -7.0]])
        self.assert_inputs_rejected(self.valid_network(), extra, r"extra experimental=\['D'\]")

    def test_disconnected_network_is_rejected(self):
        network = self.network(
            [
                ["A", "B", 1.0, 0.2],
                ["C", "D", 1.0, 0.2],
            ]
        )
        experimental = self.experimental([["A", -1.0], ["B", -2.0], ["C", -3.0], ["D", -4.0]])

        with self.assertRaisesRegex(ValueError, "disconnected"):
            self.rbfe.analyze(network, experimental)

    def test_repeated_unordered_pair_is_rejected_for_cycle_diagnostics(self):
        network = self.network(
            [
                ["A", "B", 1.0, 0.2],
                ["B", "A", -1.0, 0.2],
                ["B", "C", 2.0, 0.2],
            ]
        )

        with self.assertRaisesRegex(ValueError, "one aggregated edge per unordered ligand pair"):
            self.rbfe.cycle_closure_diagnostics(network, max_cycle_length=3)


def experimental_columns():
    return ["ligand", "DG_exp_kJ_mol"]


if __name__ == "__main__":
    unittest.main()
