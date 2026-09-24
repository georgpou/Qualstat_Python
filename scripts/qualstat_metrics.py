"""Statistics used by QualStat. Each named metric has one function."""

from __future__ import annotations

from collections import Counter
import math
from dataclasses import dataclass
from typing import Callable, Sequence


METRIC_DESCRIPTIONS = {
    "MAD": "Mean absolute deviation.",
    "MADtr": "Mean absolute deviation after removing the mean signed error.",
    "r2": "Signed squared Pearson correlation (negative correlation gives negative r2).",
    "PI": "Legacy pair-ordering index weighted by experimental separation; calculated ties add no numerator.",
    "RMSD": "Root mean square deviation.",
    "MSD": "Mean signed deviation, calculated minus experimental.",
    "Median": "Median signed deviation; lower middle value for an even sample.",
    "MQ": "Mean calculated-to-experimental quotient.",
    "Q": "Squared-error sum divided by the squared experimental-value sum.",
    "slope": "Regression slope with experimental values on the x axis.",
    "inter": "Intercept of the same regression used for slope.",
    "multi": "Regression slope multiplied by the mean experimental value.",
    "AbsMed": "Median absolute prediction error; lower middle value for an even sample.",
    "tau": "Legacy Kendall tau; prediction ties count as discordant.",
    "regMAD": "MAD from a regression of experimental on calculated values.",
    "ROCar": "Legacy 21-threshold ROC area; the most frequent experimental value is the inactive class.",
    "taux": "Legacy Kendall tau after excluding insignificant pairs.",
    "taur": "Sign agreement statistic for relative free energies.",
    "taurx": "Relative-energy sign agreement after excluding insignificant values.",
    "r22": "Signed squared Pearson correlation after adding sign-negated data.",
    "slope2": "Regression slope after adding sign-negated data.",
    "max": "Maximum absolute deviation.",
    "R": "Pearson correlation coefficient.",
    "rho": "Spearman rank correlation with average ranks for ties.",
    "rho2": "Origin-symmetric Spearman correlation after adding sign-negated data (custom RBFE metric).",
}


@dataclass(frozen=True)
class Dataset:
    """CSV rows kept as aligned columns."""

    ligand_a: tuple[str, ...]
    ligand_b: tuple[str, ...] | None
    calculated: tuple[float, ...]
    calculated_uncertainty: tuple[float, ...]
    experimental: tuple[float, ...]
    experimental_uncertainty: tuple[float, ...]


@dataclass(frozen=True)
class MetricContext:
    """Pair and row selections for one dataset or bootstrap sample."""

    tau_pairs: tuple[tuple[int, int], ...]
    taux_pairs: tuple[tuple[int, int], ...]
    taur_indices: tuple[int, ...]
    taurx_indices: tuple[int, ...]
    roc_active: tuple[bool, ...]


def _safe_divide(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator != 0.0 else math.nan


def _mean(values: Sequence[float]) -> float:
    return math.fsum(values) / len(values)


def _lower_median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    return ordered[(len(ordered) - 1) // 2]


def _pearson(x: Sequence[float], y: Sequence[float]) -> float:
    mean_x = _mean(x)
    mean_y = _mean(y)
    sum_xx = math.fsum((value - mean_x) ** 2 for value in x)
    sum_yy = math.fsum((value - mean_y) ** 2 for value in y)
    numerator = math.fsum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    return _safe_divide(numerator, math.sqrt(sum_xx * sum_yy))


def _slope(calculated: Sequence[float], experimental: Sequence[float]) -> float:
    mean_calculated = _mean(calculated)
    mean_experimental = _mean(experimental)
    covariance = math.fsum(
        (a - mean_calculated) * (b - mean_experimental)
        for a, b in zip(calculated, experimental)
    )
    variance_experimental = math.fsum(
        (value - mean_experimental) ** 2 for value in experimental
    )
    return _safe_divide(covariance, variance_experimental)


def _average_ranks(values: Sequence[float]) -> tuple[float, ...]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = ((start + 1) + end) / 2.0
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return tuple(ranks)


def build_metric_context(data: Dataset, multiplier: float) -> MetricContext:
    """Select eligible comparisons using the current rows."""
    n_records = len(data.calculated)
    tau_pairs = tuple(
        (i, j)
        for i in range(n_records)
        for j in range(i + 1, n_records)
        if data.experimental[i] != data.experimental[j]
    )
    taux_pairs = []
    for i in range(n_records):
        for j in range(i + 1, n_records):
            calculated_cutoff = multiplier * math.hypot(
                data.calculated_uncertainty[i], data.calculated_uncertainty[j]
            )
            experimental_cutoff = multiplier * math.hypot(
                data.experimental_uncertainty[i], data.experimental_uncertainty[j]
            )
            if (
                abs(data.calculated[i] - data.calculated[j]) > calculated_cutoff
                and abs(data.experimental[i] - data.experimental[j]) > experimental_cutoff
            ):
                taux_pairs.append((i, j))

    taur_indices = tuple(i for i, value in enumerate(data.experimental) if value != 0.0)
    taurx_indices = tuple(
        i
        for i in range(n_records)
        if abs(data.calculated[i]) > multiplier * data.calculated_uncertainty[i]
        and abs(data.experimental[i]) > multiplier * data.experimental_uncertainty[i]
    )

    counts = Counter(data.experimental)
    inactive = max(
        counts,
        key=lambda value: (counts[value], -data.experimental.index(value)),
    )
    roc_active = tuple(value != inactive for value in data.experimental)
    return MetricContext(
        tau_pairs=tau_pairs,
        taux_pairs=tuple(taux_pairs),
        taur_indices=taur_indices,
        taurx_indices=taurx_indices,
        roc_active=roc_active,
    )


def _ordering_statistic(
    calculated: Sequence[float],
    experimental: Sequence[float],
    pairs: Sequence[tuple[int, int]],
) -> float:
    if not pairs:
        return math.nan
    score = 0.0
    for i, j in pairs:
        concordant = (
            calculated[i] < calculated[j] and experimental[i] < experimental[j]
        ) or (
            calculated[i] > calculated[j] and experimental[i] > experimental[j]
        )
        score += 1.0 if concordant else -1.0
    return score / len(pairs)


def _relative_ordering_statistic(
    calculated: Sequence[float],
    experimental: Sequence[float],
    indices: Sequence[int],
) -> float:
    if not indices:
        return math.nan
    score = 0.0
    for index in indices:
        same_sign = (
            calculated[index] < 0.0 and experimental[index] < 0.0
        ) or (
            calculated[index] > 0.0 and experimental[index] > 0.0
        )
        score += 1.0 if same_sign else -1.0
    return score / len(indices)


def _legacy_roc(
    calculated: Sequence[float],
    experimental: Sequence[float],
    active: Sequence[bool],
) -> float:
    maximum = max(calculated)
    minimum = experimental[0]
    for value in calculated:
        if minimum > value:
            minimum = value
    step = (maximum - minimum) / 20.0
    points = []
    for index in range(21):
        threshold = minimum + index * step
        true_positive = false_positive = false_negative = true_negative = 0.0
        for value, is_active in zip(calculated, active):
            if value < threshold and is_active:
                true_positive += 1.0
            elif value < threshold and not is_active:
                false_positive += 1.0
            elif value > threshold and is_active:
                false_negative += 1.0
            elif value > threshold and not is_active:
                true_negative += 1.0
        false_rate = _safe_divide(false_positive, false_positive + true_negative)
        true_rate = _safe_divide(true_positive, true_positive + false_negative)
        if not math.isfinite(false_rate) or not math.isfinite(true_rate):
            return math.nan
        points.append((false_rate, true_rate))
    area = math.fsum(
        (points[i][0] - points[i - 1][0]) * (points[i][1] + points[i - 1][1])
        for i in range(1, len(points))
    )
    return area * 0.5


def _residuals(calculated, experimental):
    """Return calculated minus experimental values."""
    return tuple(a - b for a, b in zip(calculated, experimental))


def _mirrored(calculated, experimental):
    """Add a sign-reversed copy of each pair."""
    return (
        tuple(calculated) + tuple(-value for value in calculated),
        tuple(experimental) + tuple(-value for value in experimental),
    )


def metric_mad(calculated, experimental, context):
    """Mean absolute prediction error."""
    return _mean(tuple(abs(value) for value in _residuals(calculated, experimental)))


def metric_madtr(calculated, experimental, context):
    """Mean absolute error after removing the mean signed error."""
    residuals = _residuals(calculated, experimental)
    bias = _mean(residuals)
    return _mean(tuple(abs(value - bias) for value in residuals))


def metric_r(calculated, experimental, context):
    """Pearson correlation."""
    return _pearson(calculated, experimental)


def metric_r2(calculated, experimental, context):
    """Legacy signed squared Pearson correlation."""
    correlation = _pearson(calculated, experimental)
    return correlation * abs(correlation)


def metric_pi(calculated, experimental, context):
    """Legacy ordering index weighted by experimental separation."""
    weighted_sum = 0.0
    concordance_sum = 0.0
    for i in range(len(calculated)):
        for j in range(len(calculated)):
            expected_difference = experimental[j] - experimental[i]
            calculated_difference = calculated[j] - calculated[i]
            if calculated_difference != 0.0:
                contribution = abs(expected_difference)
                concordance_sum += (
                    -contribution
                    if expected_difference / calculated_difference < 0.0
                    else contribution
                )
            weighted_sum += abs(expected_difference)
    return _safe_divide(concordance_sum, weighted_sum)


def metric_rmsd(calculated, experimental, context):
    """Root mean squared prediction error."""
    residuals = _residuals(calculated, experimental)
    return math.sqrt(_mean(tuple(value * value for value in residuals)))


def metric_msd(calculated, experimental, context):
    """Mean signed prediction error."""
    return _mean(_residuals(calculated, experimental))


def metric_median(calculated, experimental, context):
    """Lower median of the signed errors."""
    return _lower_median(_residuals(calculated, experimental))


def metric_mq(calculated, experimental, context):
    """Mean calculated-to-experimental quotient."""
    if any(value == 0.0 for value in experimental):
        return math.nan
    return _mean(tuple(a / b for a, b in zip(calculated, experimental)))


def metric_q(calculated, experimental, context):
    """Squared-error sum divided by squared experimental sum."""
    return _safe_divide(
        math.fsum(value * value for value in _residuals(calculated, experimental)),
        math.fsum(value * value for value in experimental),
    )


def metric_slope(calculated, experimental, context):
    """Slope of calculated values against experimental values."""
    return _slope(calculated, experimental)


def metric_inter(calculated, experimental, context):
    """Intercept of calculated values against experimental values."""
    slope = _slope(calculated, experimental)
    return _mean(calculated) - slope * _mean(experimental)


def metric_multi(calculated, experimental, context):
    """Regression slope times mean experimental value."""
    return _slope(calculated, experimental) * _mean(experimental)


def metric_absmed(calculated, experimental, context):
    """Lower median of absolute prediction errors."""
    residuals = _residuals(calculated, experimental)
    return _lower_median(tuple(abs(value) for value in residuals))


def metric_tau(calculated, experimental, context):
    """Legacy pair-ordering score."""
    return _ordering_statistic(calculated, experimental, context.tau_pairs)


def metric_regmad(calculated, experimental, context):
    """Mean absolute error after regressing experimental on calculated."""
    slope = _slope(experimental, calculated)
    intercept = _mean(experimental) - slope * _mean(calculated)
    return _mean(
        tuple(
            abs(expected - (predicted * slope + intercept))
            for predicted, expected in zip(calculated, experimental)
        )
    )


def metric_rocar(calculated, experimental, context):
    """Legacy 21-threshold ROC quantity."""
    return _legacy_roc(calculated, experimental, context.roc_active)


def metric_taux(calculated, experimental, context):
    """Pair-ordering score for significant pairs."""
    return _ordering_statistic(calculated, experimental, context.taux_pairs)


def metric_taur(calculated, experimental, context):
    """Sign agreement for nonzero experimental values."""
    return _relative_ordering_statistic(calculated, experimental, context.taur_indices)


def metric_taurx(calculated, experimental, context):
    """Sign agreement for significant values."""
    return _relative_ordering_statistic(calculated, experimental, context.taurx_indices)


def metric_r22(calculated, experimental, context):
    """Signed squared Pearson correlation on mirrored pairs."""
    doubled_calculated, doubled_experimental = _mirrored(calculated, experimental)
    correlation = _pearson(doubled_calculated, doubled_experimental)
    return correlation * abs(correlation)


def metric_slope2(calculated, experimental, context):
    """Regression slope on mirrored pairs."""
    doubled_calculated, doubled_experimental = _mirrored(calculated, experimental)
    return _slope(doubled_calculated, doubled_experimental)


def metric_max(calculated, experimental, context):
    """Largest absolute prediction error."""
    return max(abs(value) for value in _residuals(calculated, experimental))


def metric_rho(calculated, experimental, context):
    """Spearman correlation with average ranks for ties."""
    return _pearson(_average_ranks(calculated), _average_ranks(experimental))


def metric_rho2(calculated, experimental, context):
    """Spearman correlation on mirrored pairs."""
    doubled_calculated, doubled_experimental = _mirrored(calculated, experimental)
    return _pearson(
        _average_ranks(doubled_calculated),
        _average_ranks(doubled_experimental),
    )


# Each function accepts (calculated, experimental, context).
MetricFunction = Callable[[Sequence[float], Sequence[float], MetricContext], float]
METRICS: dict[str, MetricFunction] = {
    "MAD": metric_mad,
    "MADtr": metric_madtr,
    "r2": metric_r2,
    "PI": metric_pi,
    "RMSD": metric_rmsd,
    "MSD": metric_msd,
    "Median": metric_median,
    "MQ": metric_mq,
    "Q": metric_q,
    "slope": metric_slope,
    "inter": metric_inter,
    "multi": metric_multi,
    "AbsMed": metric_absmed,
    "tau": metric_tau,
    "regMAD": metric_regmad,
    "ROCar": metric_rocar,
    "taux": metric_taux,
    "taur": metric_taur,
    "taurx": metric_taurx,
    "r22": metric_r22,
    "slope2": metric_slope2,
    "max": metric_max,
    "R": metric_r,
    "rho": metric_rho,
    "rho2": metric_rho2,
}


def evaluate_metric(
    name: str,
    calculated: Sequence[float],
    experimental: Sequence[float],
    context: MetricContext,
) -> float:
    """Calculate a metric selected by its configuration name."""
    try:
        metric = METRICS[name]
    except KeyError as exc:
        raise ValueError(f"unknown statistic: {name}") from exc
    return metric(calculated, experimental, context)
