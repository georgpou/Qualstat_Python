from __future__ import annotations

import unittest

from tests.support import load_script_module


class QualStatBootstrapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qs = load_script_module("qualstat.py", "qualstat_bootstrap_test")

    def dataset(self, calculated_uncertainty, experimental_uncertainty):
        return self.qs.Dataset(
            ligand_a=("A", "B", "C"),
            ligand_b=None,
            calculated=(1.0, 2.0, 3.0),
            calculated_uncertainty=tuple(calculated_uncertainty),
            experimental=(1.1, 1.9, 3.2),
            experimental_uncertainty=tuple(experimental_uncertainty),
        )

    def test_zero_uncertainty_is_a_deterministic_bootstrap_oracle(self):
        data = self.dataset([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

        results = self.qs.run_statistics(
            data,
            metric_names=("MAD", "R", "rho", "tau", "slope", "max"),
            rounds=20,
            seed=123,
            multiplier=1.645,
        )

        for result in results:
            self.assertEqual(result.valid_bootstraps, 20, result.name)
            self.assertEqual(result.bootstrap_sd, 0.0, result.name)

    def test_fixed_seed_reproduces_nonzero_bootstrap_and_other_seed_changes_sd(self):
        data = self.dataset([0.2, 0.3, 0.4], [0.1, 0.2, 0.3])
        metric_names = ("MAD", "R", "rho")

        first = self.qs.run_statistics(data, metric_names, 40, 2026, 1.645)
        repeated = self.qs.run_statistics(data, metric_names, 40, 2026, 1.645)
        different = self.qs.run_statistics(data, metric_names, 40, 2027, 1.645)

        self.assertEqual(first, repeated)
        self.assertTrue(
            any(
                left.bootstrap_sd != right.bootstrap_sd
                for left, right in zip(first, different)
            )
        )


if __name__ == "__main__":
    unittest.main()
