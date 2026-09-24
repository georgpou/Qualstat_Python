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

    def test_zero_uncertainty_still_gives_sampling_variation(self):
        data = self.dataset([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

        results = self.qs.run_statistics(
            data,
            metric_names=("MAD",),
            rounds=100,
            seed=123,
            multiplier=1.645,
        )

        self.assertEqual(results[0].valid_bootstraps, 100)
        self.assertGreater(results[0].bootstrap_sd, 0.0)

    def test_each_sample_keeps_calculated_and_experimental_rows_paired(self):
        data = self.qs.Dataset(
            ligand_a=("A", "B", "C"),
            ligand_b=None,
            calculated=(1.0, 4.0, 9.0),
            calculated_uncertainty=(0.2, 0.4, 0.9),
            experimental=(2.0, 5.0, 10.0),
            experimental_uncertainty=(0.3, 0.5, 1.0),
        )

        result, = self.qs.run_statistics(data, ("MAD",), 100, 123, 1.645)

        self.assertEqual(result.estimate, 1.0)
        self.assertEqual(result.valid_bootstraps, 100)
        self.assertEqual(result.bootstrap_sd, 0.0)

    def test_duplicate_rows_are_excluded_from_tau_pairs_in_each_round(self):
        data = self.qs.Dataset(
            ligand_a=("A", "B"),
            ligand_b=None,
            calculated=(1.0, 2.0),
            calculated_uncertainty=(0.0, 0.0),
            experimental=(1.0, 2.0),
            experimental_uncertainty=(0.0, 0.0),
        )

        result, = self.qs.run_statistics(data, ("tau",), 20, 123, 1.645)

        self.assertEqual(result.estimate, 1.0)
        self.assertEqual(result.valid_bootstraps, 10)
        self.assertEqual(result.bootstrap_sd, 0.0)

    def test_resampled_rows_keep_names_values_and_uncertainties_together(self):
        data = self.qs.Dataset(
            ligand_a=("A", "B", "C"),
            ligand_b=("B", "C", "D"),
            calculated=(1.0, 2.0, 3.0),
            calculated_uncertainty=(0.1, 0.2, 0.3),
            experimental=(4.0, 5.0, 6.0),
            experimental_uncertainty=(0.4, 0.5, 0.6),
        )

        sampled = self.qs._resample_rows(data, (2, 2, 0))

        self.assertEqual(sampled.ligand_a, ("C", "C", "A"))
        self.assertEqual(sampled.ligand_b, ("D", "D", "B"))
        self.assertEqual(sampled.calculated, (3.0, 3.0, 1.0))
        self.assertEqual(sampled.calculated_uncertainty, (0.3, 0.3, 0.1))
        self.assertEqual(sampled.experimental, (6.0, 6.0, 4.0))
        self.assertEqual(sampled.experimental_uncertainty, (0.6, 0.6, 0.4))

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
