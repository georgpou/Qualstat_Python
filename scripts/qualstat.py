#!/usr/bin/env python3
"""Legacy-compatible QualStat analysis with YAML/CSV input."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import math
import random
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import yaml


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

_CANONICAL_METRICS = {name.lower(): name for name in METRIC_DESCRIPTIONS}

ABFE_HEADER = (
    "ligand",
    "calculated",
    "calculated_uncertainty",
    "experimental",
    "experimental_uncertainty",
)
RBFE_HEADER = (
    "ligand_a",
    "ligand_b",
    "calculated",
    "calculated_uncertainty",
    "experimental",
    "experimental_uncertainty",
)

TEMPLATE_TEXT = """# CSV file containing ABFE values or directed RBFE transformations.
input_file: results.csv

# Plain-text analysis report. A .log suffix is recommended.
output_file: qualstat_results.log

# Select abfe or rbfe; this determines the required CSV columns.
analysis_type: rbfe

# Unit label printed in the report; input values are not converted.
energy_unit: kJ/mol

# Each CSV uncertainty is the one-sigma uncertainty of the reported value.
# For a mean of n independent repeats estimated from their sample SD s, use
# the standard error s / sqrt(n); for triplicates, use s / sqrt(3).

# Number of legacy parametric uncertainty-propagation samples.
bootstrap_rounds: 1000

# Integer seed for reproducible results, or null for a random seed.
random_seed: 2026

# Multiplier used by taux and taurx (1.645 = 90%; 1.96 = 95%).
significance_multiplier: 1.645

# Every supported metric is listed. Change values to true or false.
statistics:
  MAD: true
  MADtr: false
  r2: false
  PI: false
  RMSD: true
  MSD: false
  Median: false
  MQ: false
  Q: false
  slope: false
  inter: false
  multi: false
  AbsMed: false
  tau: false
  regMAD: false
  ROCar: false
  taux: false
  taur: false
  taurx: true
  r22: true
  slope2: false
  max: false
  R: false
  rho: false
  rho2: false

# RBFE only: enumerate unique simple cycles and calculate closure errors.
cycle_analysis: true
"""


@dataclass(frozen=True)
class Config:
    input_file: Path
    output_file: Path
    analysis_type: str
    energy_unit: str
    bootstrap_rounds: int
    random_seed: int | None
    significance_multiplier: float
    statistics: tuple[str, ...]
    cycle_analysis: bool


@dataclass(frozen=True)
class Dataset:
    ligand_a: tuple[str, ...]
    ligand_b: tuple[str, ...] | None
    calculated: tuple[float, ...]
    calculated_uncertainty: tuple[float, ...]
    experimental: tuple[float, ...]
    experimental_uncertainty: tuple[float, ...]


@dataclass(frozen=True)
class MetricContext:
    calculated_uncertainty: tuple[float, ...]
    experimental_uncertainty: tuple[float, ...]
    tau_pairs: tuple[tuple[int, int], ...]
    taux_pairs: tuple[tuple[int, int], ...]
    taur_indices: tuple[int, ...]
    taurx_indices: tuple[int, ...]
    roc_active: tuple[bool, ...]


@dataclass(frozen=True)
class MetricResult:
    name: str
    estimate: float
    bootstrap_sd: float
    valid_bootstraps: int


@dataclass(frozen=True)
class CycleResult:
    ligands: tuple[str, ...]
    closure: float
    uncertainty: float


def write_template(
    path: Path = Path("qualstat_template.yaml"), *, force: bool = False
) -> Path:
    path = Path(path)
    try:
        with path.open("w" if force else "x", encoding="utf-8") as handle:
            handle.write(TEMPLATE_TEXT)
    except FileExistsError as exc:
        raise ValueError(f"{path} already exists; use --force to overwrite it") from exc
    return path


def _require_mapping(document: object) -> dict:
    if not isinstance(document, dict):
        raise ValueError("YAML root must be a mapping")
    return document


def _required(document: dict, key: str):
    if key not in document:
        raise ValueError(f"missing required YAML setting: {key}")
    return document[key]


def load_config(path: Path) -> Config:
    path = Path(path).expanduser().resolve()
    with path.open(encoding="utf-8") as handle:
        document = _require_mapping(yaml.safe_load(handle))

    analysis_type = str(_required(document, "analysis_type")).lower()
    if analysis_type not in {"abfe", "rbfe"}:
        raise ValueError("analysis_type must be 'abfe' or 'rbfe'")

    rounds = _required(document, "bootstrap_rounds")
    if type(rounds) is not int or rounds < 2:
        raise ValueError("bootstrap_rounds must be an integer of at least 2")

    raw_statistics = _required(document, "statistics")
    if not isinstance(raw_statistics, dict):
        raise ValueError("statistics must be a mapping of metric names to booleans")
    selected = []
    seen = set()
    for raw_name, enabled in raw_statistics.items():
        key = str(raw_name).lower()
        if key not in _CANONICAL_METRICS:
            raise ValueError(f"unknown statistic: {raw_name}")
        canonical = _CANONICAL_METRICS[key]
        if canonical in seen:
            raise ValueError(f"duplicate statistic: {raw_name}")
        seen.add(canonical)
        if type(enabled) is not bool:
            raise ValueError(f"statistics.{raw_name} must be true or false")
        if enabled:
            selected.append(canonical)
    if not selected:
        raise ValueError("at least one statistic must be true")

    random_seed = document.get("random_seed", 2026)
    if random_seed is not None and type(random_seed) is not int:
        raise ValueError("random_seed must be an integer or null")

    multiplier = document.get("significance_multiplier", 1.645)
    if type(multiplier) not in {int, float} or not math.isfinite(multiplier) or multiplier <= 0:
        raise ValueError("significance_multiplier must be a positive finite number")

    energy_unit = document.get("energy_unit", "kJ/mol")
    if not isinstance(energy_unit, str) or not energy_unit.strip():
        raise ValueError("energy_unit must be a non-empty string")

    cycle_analysis = document.get("cycle_analysis", False)
    if type(cycle_analysis) is not bool:
        raise ValueError("cycle_analysis must be true or false")

    base = path.parent
    input_file = (base / str(_required(document, "input_file"))).expanduser().resolve()
    output_file = (base / str(_required(document, "output_file"))).expanduser().resolve()
    if _same_file(output_file, input_file) or _same_file(output_file, path):
        raise ValueError("output_file must differ from the input CSV and YAML configuration")
    return Config(
        input_file=input_file,
        output_file=output_file,
        analysis_type=analysis_type,
        energy_unit=energy_unit.strip(),
        bootstrap_rounds=rounds,
        random_seed=random_seed,
        significance_multiplier=float(multiplier),
        statistics=tuple(selected),
        cycle_analysis=cycle_analysis,
    )


def _same_file(first: Path, second: Path) -> bool:
    if first == second:
        return True
    try:
        return first.samefile(second)
    except OSError:
        return False


def _parse_number(value: str, row_number: int, field: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"row {row_number}: {field} must be numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"row {row_number}: {field} must be finite")
    if field.endswith("uncertainty") and number < 0:
        raise ValueError(f"row {row_number}: {field} must be non-negative")
    return number


def read_dataset(config: Config) -> Dataset:
    expected_header = ABFE_HEADER if config.analysis_type == "abfe" else RBFE_HEADER
    ligand_a = []
    ligand_b = [] if config.analysis_type == "rbfe" else None
    calculated = []
    calculated_uncertainty = []
    experimental = []
    experimental_uncertainty = []
    pairs = set()
    abfe_ligands = set()

    with config.input_file.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header is None or tuple(header) != expected_header:
            raise ValueError(
                "CSV header must be exactly: " + ",".join(expected_header)
            )
        for row_number, row in enumerate(reader, start=2):
            if not row or all(not cell.strip() for cell in row):
                continue
            if len(row) != len(expected_header):
                raise ValueError(
                    f"row {row_number}: expected {len(expected_header)} columns, got {len(row)}"
                )
            first = row[0].strip()
            if not first:
                raise ValueError(f"row {row_number}: {expected_header[0]} must not be blank")
            ligand_a.append(first)
            offset = 1
            if ligand_b is not None:
                second = row[1].strip()
                if not second:
                    raise ValueError(f"row {row_number}: ligand_b must not be blank")
                if first == second:
                    raise ValueError(f"row {row_number}: RBFE self-edge {first} -> {second} is invalid")
                pair = frozenset((first, second))
                if pair in pairs:
                    raise ValueError(f"row {row_number}: duplicate RBFE pair {first}, {second}")
                pairs.add(pair)
                ligand_b.append(second)
                offset = 2
            elif first in abfe_ligands:
                raise ValueError(f"row {row_number}: duplicate ABFE ligand {first}")
            else:
                abfe_ligands.add(first)
            numeric_fields = expected_header[offset:]
            values = [
                _parse_number(value.strip(), row_number, field)
                for value, field in zip(row[offset:], numeric_fields)
            ]
            calculated.append(values[0])
            calculated_uncertainty.append(values[1])
            experimental.append(values[2])
            experimental_uncertainty.append(values[3])

    if not calculated:
        raise ValueError("CSV contains no data records")
    return Dataset(
        ligand_a=tuple(ligand_a),
        ligand_b=tuple(ligand_b) if ligand_b is not None else None,
        calculated=tuple(calculated),
        calculated_uncertainty=tuple(calculated_uncertainty),
        experimental=tuple(experimental),
        experimental_uncertainty=tuple(experimental_uncertainty),
    )


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
        calculated_uncertainty=data.calculated_uncertainty,
        experimental_uncertainty=data.experimental_uncertainty,
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


def evaluate_metric(
    name: str,
    calculated: Sequence[float],
    experimental: Sequence[float],
    context: MetricContext,
) -> float:
    residuals = tuple(a - b for a, b in zip(calculated, experimental))
    if name == "MAD":
        return _mean(tuple(abs(value) for value in residuals))
    if name == "MADtr":
        systematic_error = _mean(residuals)
        return _mean(tuple(abs(value - systematic_error) for value in residuals))
    if name == "R":
        return _pearson(calculated, experimental)
    if name == "r2":
        correlation = _pearson(calculated, experimental)
        return correlation * abs(correlation)
    if name == "PI":
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
    if name == "RMSD":
        return math.sqrt(_mean(tuple(value * value for value in residuals)))
    if name == "MSD":
        return _mean(residuals)
    if name == "Median":
        return _lower_median(residuals)
    if name == "MQ":
        if any(value == 0.0 for value in experimental):
            return math.nan
        return _mean(tuple(a / b for a, b in zip(calculated, experimental)))
    if name == "Q":
        return _safe_divide(
            math.fsum(value * value for value in residuals),
            math.fsum(value * value for value in experimental),
        )
    if name == "slope":
        return _slope(calculated, experimental)
    if name == "inter":
        slope = _slope(calculated, experimental)
        return _mean(calculated) - slope * _mean(experimental)
    if name == "multi":
        return _slope(calculated, experimental) * _mean(experimental)
    if name == "AbsMed":
        return _lower_median(tuple(abs(value) for value in residuals))
    if name == "tau":
        return _ordering_statistic(calculated, experimental, context.tau_pairs)
    if name == "regMAD":
        slope = _slope(experimental, calculated)
        intercept = _mean(experimental) - slope * _mean(calculated)
        return _mean(
            tuple(
                abs(expected - (predicted * slope + intercept))
                for predicted, expected in zip(calculated, experimental)
            )
        )
    if name == "ROCar":
        return _legacy_roc(calculated, experimental, context.roc_active)
    if name == "taux":
        return _ordering_statistic(calculated, experimental, context.taux_pairs)
    if name == "taur":
        return _relative_ordering_statistic(
            calculated, experimental, context.taur_indices
        )
    if name == "taurx":
        return _relative_ordering_statistic(
            calculated, experimental, context.taurx_indices
        )
    if name in {"r22", "slope2", "rho2"}:
        doubled_calculated = tuple(calculated) + tuple(-value for value in calculated)
        doubled_experimental = tuple(experimental) + tuple(-value for value in experimental)
        if name == "r22":
            correlation = _pearson(doubled_calculated, doubled_experimental)
            return correlation * abs(correlation)
        if name == "slope2":
            return _slope(doubled_calculated, doubled_experimental)
        return _pearson(
            _average_ranks(doubled_calculated),
            _average_ranks(doubled_experimental),
        )
    if name == "max":
        return max(abs(value) for value in residuals)
    if name == "rho":
        return _pearson(_average_ranks(calculated), _average_ranks(experimental))
    raise ValueError(f"unknown statistic: {name}")


def run_statistics(
    data: Dataset,
    metric_names: Sequence[str],
    rounds: int,
    seed: int | None,
    multiplier: float,
) -> tuple[MetricResult, ...]:
    context = build_metric_context(data, multiplier)
    estimates = {
        name: evaluate_metric(name, data.calculated, data.experimental, context)
        for name in metric_names
    }
    samples = {name: [] for name in metric_names}
    rng = random.Random(seed)
    for _ in range(rounds):
        calculated = tuple(
            rng.gauss(value, error)
            for value, error in zip(
                data.calculated, data.calculated_uncertainty
            )
        )
        experimental = tuple(
            rng.gauss(value, error)
            for value, error in zip(
                data.experimental, data.experimental_uncertainty
            )
        )
        for name in metric_names:
            value = evaluate_metric(name, calculated, experimental, context)
            if math.isfinite(value):
                samples[name].append(value)

    return tuple(
        MetricResult(
            name=name,
            estimate=estimates[name],
            bootstrap_sd=(
                statistics.stdev(samples[name])
                if len(samples[name]) >= 2
                else math.nan
            ),
            valid_bootstraps=len(samples[name]),
        )
        for name in metric_names
    )


def canonical_cycle(nodes: Sequence[str]) -> tuple[str, ...]:
    sequence = tuple(nodes)
    rotations = [sequence[index:] + sequence[:index] for index in range(len(sequence))]
    reversed_sequence = tuple(reversed(sequence))
    rotations.extend(
        reversed_sequence[index:] + reversed_sequence[:index]
        for index in range(len(reversed_sequence))
    )
    return min(rotations)


def find_cycles(data: Dataset) -> tuple[CycleResult, ...]:
    if data.ligand_b is None:
        raise ValueError("cycle analysis requires RBFE data")

    adjacency: dict[str, set[str]] = {}
    edges = {}
    for ligand_a, ligand_b, value, uncertainty in zip(
        data.ligand_a,
        data.ligand_b,
        data.calculated,
        data.calculated_uncertainty,
    ):
        adjacency.setdefault(ligand_a, set()).add(ligand_b)
        adjacency.setdefault(ligand_b, set()).add(ligand_a)
        edges[frozenset((ligand_a, ligand_b))] = (
            ligand_a,
            ligand_b,
            value,
            uncertainty,
        )

    cycles: set[tuple[str, ...]] = set()

    for start in sorted(adjacency):
        path = [start]
        visited = {start}
        stack = [(start, iter(sorted(adjacency[start])))]
        while stack:
            _, neighbors = stack[-1]
            try:
                neighbor = next(neighbors)
            except StopIteration:
                stack.pop()
                departed = path.pop()
                if departed != start:
                    visited.remove(departed)
                continue
            if neighbor == start:
                if len(path) >= 3:
                    cycles.add(canonical_cycle(path))
                continue
            if neighbor in visited:
                continue
            path.append(neighbor)
            visited.add(neighbor)
            stack.append((neighbor, iter(sorted(adjacency[neighbor]))))

    results = []
    for cycle in sorted(cycles):
        closure_terms = []
        errors = []
        pairs = zip(cycle, cycle[1:] + cycle[:1])
        for source, target in pairs:
            stored_source, stored_target, value, uncertainty = edges[
                frozenset((source, target))
            ]
            closure_terms.append(
                value if (source, target) == (stored_source, stored_target) else -value
            )
            errors.append(uncertainty)
        results.append(
            CycleResult(
                ligands=cycle,
                closure=math.fsum(closure_terms),
                uncertainty=math.sqrt(math.fsum(error * error for error in errors)),
            )
        )
    return tuple(results)


def _format_number(value: float) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.10g}"


def render_report(
    config: Config,
    data: Dataset,
    results: Sequence[MetricResult],
    cycles: Sequence[CycleResult],
) -> str:
    seed = "random" if config.random_seed is None else str(config.random_seed)
    lines = [
        "Python QualStat",
        f"Analysis: {config.analysis_type.upper()}",
        f"Input: {config.input_file}",
        f"Energy unit: {config.energy_unit}",
        f"Records: {len(data.calculated)}",
        f"Propagation rounds: {config.bootstrap_rounds}",
        f"Seed: {seed}",
        "Uncertainty method: legacy independent Gaussian uncertainty propagation "
        "(records are not resampled)",
        "",
        f"{'Metric':<8}{'Estimate':>16}{'Propagation SD':>16}{'Valid rounds':>14}",
        "-" * 54,
    ]
    for result in results:
        lines.append(
            f"{result.name:<8}"
            f"{_format_number(result.estimate):>16}"
            f"{_format_number(result.bootstrap_sd):>16}"
            f"{result.valid_bootstraps:>14}"
        )

    undefined = [
        result.name
        for result in results
        if not math.isfinite(result.estimate) or not math.isfinite(result.bootstrap_sd)
    ]
    if undefined:
        lines.extend(
            (
                "",
                "Note: undefined value(s) reported as nan for: " + ", ".join(undefined),
            )
        )

    if config.analysis_type == "abfe" and config.cycle_analysis:
        lines.extend(("", "Note: cycle_analysis ignored for ABFE data."))
    elif config.analysis_type == "rbfe" and config.cycle_analysis:
        lines.extend(("", "RBFE cycle analysis", "==================="))
        if not cycles:
            lines.append("No cycles found.")
        for number, cycle in enumerate(cycles, start=1):
            path = " -> ".join(cycle.ligands + cycle.ligands[:1])
            lines.extend(
                (
                    "",
                    f"Cycle #{number}",
                    "",
                    f"{path} = {cycle.closure:.3f} +/- "
                    f"{cycle.uncertainty:.3f} {config.energy_unit}",
                )
            )
    return "\n".join(lines) + "\n"


def run(config_path: Path) -> str:
    config = load_config(config_path)
    data = read_dataset(config)
    results = run_statistics(
        data=data,
        metric_names=config.statistics,
        rounds=config.bootstrap_rounds,
        seed=config.random_seed,
        multiplier=config.significance_multiplier,
    )
    cycles = (
        find_cycles(data)
        if config.analysis_type == "rbfe" and config.cycle_analysis
        else ()
    )
    report = render_report(config, data, results, cycles)
    config.output_file.write_text(report, encoding="utf-8")
    print(report, end="")
    return report


def _help_epilog() -> str:
    metrics = "\n".join(
        f"  {name:<7} {description}" for name, description in METRIC_DESCRIPTIONS.items()
    )
    return f"""CSV schemas:
  ABFE: ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty
  RBFE: ligand_a,ligand_b,calculated,calculated_uncertainty,experimental,experimental_uncertainty

Metrics:
{metrics}

RBFE cycle analysis enumerates all simple cycles; their number can grow rapidly
for dense graphs.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="qualstat.py",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Legacy-compatible QualStat analysis for ABFE and RBFE data.",
        epilog=_help_epilog(),
    )
    parser.add_argument("-h", "-help", "--help", action="help", help="show this help message")
    parser.add_argument(
        "--write-template",
        nargs="?",
        const="qualstat_template.yaml",
        metavar="PATH",
        help="write a commented YAML template and exit",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="allow --write-template to replace an existing file",
    )
    parser.add_argument("config", nargs="?", help="YAML configuration file")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.force and args.write_template is None:
        parser.error("--force requires --write-template")
    if args.config is None:
        if args.write_template is None:
            parser.error("a YAML configuration file is required")
    try:
        if args.write_template is not None:
            path = write_template(Path(args.write_template), force=args.force)
            print(f"Wrote {path}")
            return 0
        run(Path(args.config))
        return 0
    except (OSError, ValueError, yaml.YAMLError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
