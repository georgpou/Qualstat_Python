from __future__ import annotations

import math
import unittest

from tests.support import load_script_module


class QualStatMetricTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qs = load_script_module("qualstat.py", "qualstat_metrics_test")

    def dataset(
        self,
        calculated,
        experimental,
        calculated_uncertainty=None,
        experimental_uncertainty=None,
    ):
        n_records = len(calculated)
        if calculated_uncertainty is None:
            calculated_uncertainty = [0.0] * n_records
        if experimental_uncertainty is None:
            experimental_uncertainty = [0.0] * n_records
        return self.qs.Dataset(
            ligand_a=tuple(f"L{i}" for i in range(n_records)),
            ligand_b=None,
            calculated=tuple(calculated),
            calculated_uncertainty=tuple(calculated_uncertainty),
            experimental=tuple(experimental),
            experimental_uncertainty=tuple(experimental_uncertainty),
        )

    def metric(self, name, calculated, experimental, **uncertainties):
        data = self.dataset(calculated, experimental, **uncertainties)
        context = self.qs.build_metric_context(data, 1.0)
        return self.qs.evaluate_metric(name, data.calculated, data.experimental, context)

    def test_perfect_prediction_covers_all_supported_metrics(self):
        values = [1.0, 2.0, 3.0, 4.0]
        data = self.dataset(values, values)
        context = self.qs.build_metric_context(data, 1.0)
        results = {
            name: self.qs.evaluate_metric(name, data.calculated, data.experimental, context)
            for name in self.qs.METRIC_DESCRIPTIONS
        }

        for name in (
            "MAD",
            "MADtr",
            "RMSD",
            "MSD",
            "Median",
            "AbsMed",
            "Q",
            "max",
            "regMAD",
        ):
            self.assertEqual(results[name], 0.0, name)
        for name in ("R", "r2", "rho", "tau", "PI", "taux", "taur", "taurx", "r22", "slope2"):
            self.assertAlmostEqual(results[name], 1.0, msg=name)
        self.assertAlmostEqual(results["MQ"], 1.0)
        self.assertAlmostEqual(results["slope"], 1.0)
        self.assertAlmostEqual(results["inter"], 0.0)
        self.assertAlmostEqual(results["multi"], 2.5)

        self.assertEqual(set(results), set(self.qs.METRIC_DESCRIPTIONS))

    def test_inverse_order_has_negative_signed_correlations(self):
        calculated = [4.0, 3.0, 2.0, 1.0]
        experimental = [1.0, 2.0, 3.0, 4.0]

        self.assertAlmostEqual(self.metric("R", calculated, experimental), -1.0)
        self.assertAlmostEqual(self.metric("r2", calculated, experimental), -1.0)
        self.assertAlmostEqual(self.metric("rho", calculated, experimental), -1.0)
        self.assertAlmostEqual(self.metric("tau", calculated, experimental), -1.0)

    def test_even_sample_medians_use_lower_middle_value(self):
        calculated = [2.0, 3.0, 5.0, 9.0]
        experimental = [1.0, 1.0, 1.0, 1.0]

        self.assertEqual(self.metric("Median", calculated, experimental), 2.0)
        self.assertEqual(self.metric("AbsMed", calculated, experimental), 2.0)

    def test_spearman_uses_average_ranks_for_ties(self):
        calculated = [1.0, 1.0, 3.0, 4.0]
        experimental = [1.0, 2.0, 3.0, 4.0]

        # Ranks are (1.5, 1.5, 3, 4) and (1, 2, 3, 4); r = 4.5/sqrt(4.5*5).
        expected = 4.5 / math.sqrt(22.5)
        self.assertAlmostEqual(self.metric("rho", calculated, experimental), expected)

    def test_tau_excludes_experimental_ties_but_counts_calculated_ties_as_discordant(self):
        self.assertAlmostEqual(
            self.metric("tau", [1.0, 2.0, 3.0], [1.0, 1.0, 3.0]),
            1.0,
        )
        self.assertAlmostEqual(
            self.metric("tau", [1.0, 1.0, 3.0], [1.0, 2.0, 3.0]),
            1.0 / 3.0,
        )

    def test_taux_uses_strict_uncertainty_threshold(self):
        equal_cutoff = self.dataset(
            [0.0, 1.0],
            [0.0, 2.0],
            calculated_uncertainty=[0.0, 1.0],
            experimental_uncertainty=[0.0, 0.0],
        )
        context = self.qs.build_metric_context(equal_cutoff, 1.0)
        self.assertEqual(context.taux_pairs, ())
        self.assertEqual(
            self.qs.evaluate_metric("taux", equal_cutoff.calculated, equal_cutoff.experimental, context),
            0.0,
        )

        just_above = self.dataset(
            [0.0, 1.0001],
            [0.0, 2.0],
            calculated_uncertainty=[0.0, 1.0],
            experimental_uncertainty=[0.0, 0.0],
        )
        context = self.qs.build_metric_context(just_above, 1.0)
        self.assertEqual(context.taux_pairs, ((0, 1),))
        self.assertEqual(
            self.qs.evaluate_metric("taux", just_above.calculated, just_above.experimental, context),
            1.0,
        )

    def test_taur_and_taurx_select_significant_nonzero_values(self):
        data = self.dataset(
            [1.0, 0.5, -2.0],
            [1.0, -1.0, 0.0],
            calculated_uncertainty=[0.0, 0.5, 0.0],
            experimental_uncertainty=[0.0, 0.0, 2.0],
        )
        context = self.qs.build_metric_context(data, 1.0)

        self.assertEqual(context.taur_indices, (0, 1))
        self.assertEqual(context.taurx_indices, (0,))
        self.assertEqual(
            self.qs.evaluate_metric("taur", data.calculated, data.experimental, context),
            0.0,
        )
        self.assertEqual(
            self.qs.evaluate_metric("taurx", data.calculated, data.experimental, context),
            1.0,
        )

    def test_pi_gives_no_numerator_credit_to_calculated_ties(self):
        # Ordered experimental separation sums to 8; non-tied calculated pairs
        # contribute 2+2+1+1 = 6, so PI is 6/8.
        self.assertAlmostEqual(
            self.metric("PI", [1.0, 1.0, 3.0], [1.0, 2.0, 3.0]),
            0.75,
        )

    def test_undefined_metrics_are_nan(self):
        self.assertTrue(math.isnan(self.metric("MQ", [1.0, 2.0], [0.0, 1.0])))
        self.assertTrue(math.isnan(self.metric("R", [1.0, 1.0], [0.0, 1.0])))
        self.assertTrue(math.isnan(self.metric("r2", [1.0, 1.0], [0.0, 1.0])))

    def test_origin_symmetric_metrics_have_independent_values(self):
        calculated = [1.0, 2.0, 4.0]
        experimental = [2.0, 3.0, 5.0]

        # Doubling with sign-negated values gives zero means. Covariance is 56,
        # calculated sum of squares is 42, and experimental sum is 76.
        self.assertAlmostEqual(self.metric("r22", calculated, experimental), 3136.0 / 3192.0)
        self.assertAlmostEqual(self.metric("slope2", calculated, experimental), 14.0 / 19.0)

    def test_legacy_rocar_class_selection_and_fallback(self):
        normal = self.dataset([3.0, 2.0, 1.0, 0.0], [0.0, 0.0, 1.0, 1.0])
        normal_context = self.qs.build_metric_context(normal, 1.0)
        self.assertEqual(normal_context.roc_active, (False, False, True, True))
        # Twenty-one thresholds produce the hand-checkable trapezoid area 1.0.
        self.assertAlmostEqual(
            self.qs.evaluate_metric("ROCar", normal.calculated, normal.experimental, normal_context),
            1.0,
        )

        tied = self.dataset([3.0, 2.0, 1.0, 0.0], [1.0, 2.0, 1.0, 2.0])
        tied_context = self.qs.build_metric_context(tied, 1.0)
        self.assertEqual(tied_context.roc_active, (False, True, False, True))

        unique = self.dataset([1.0, 3.0, 2.0, 0.0], [10.0, 20.0, 30.0, 40.0])
        unique_context = self.qs.build_metric_context(unique, 1.0)
        self.assertEqual(unique_context.roc_active, (False, True, True, True))
        # All labels are unique, so the first experimental label is inactive;
        # the threshold path has area 1/3 for this hand-calculated ordering.
        self.assertAlmostEqual(
            self.qs.evaluate_metric("ROCar", unique.calculated, unique.experimental, unique_context),
            1.0 / 3.0,
        )


if __name__ == "__main__":
    unittest.main()
