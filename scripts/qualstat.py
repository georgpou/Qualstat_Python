#!/usr/bin/env python3
"""Legacy-compatible QualStat analysis with YAML/CSV input."""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import yaml

from qualstat_metrics import (
    Dataset,
    METRIC_DESCRIPTIONS,
    MetricContext,
    build_metric_context,
    evaluate_metric,
)


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

# Number of rounds; each round samples complete CSV records with replacement.
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


def _resample_rows(data: Dataset, indices: Sequence[int]) -> Dataset:
    """Copy whole paired records in the sampled order."""

    def pick(values: Sequence):
        return tuple(values[index] for index in indices)

    return Dataset(
        ligand_a=pick(data.ligand_a),
        ligand_b=pick(data.ligand_b) if data.ligand_b is not None else None,
        calculated=pick(data.calculated),
        calculated_uncertainty=pick(data.calculated_uncertainty),
        experimental=pick(data.experimental),
        experimental_uncertainty=pick(data.experimental_uncertainty),
    )


def run_statistics(
    data: Dataset,
    metric_names: Sequence[str],
    rounds: int,
    seed: int | None,
    multiplier: float,
) -> tuple[MetricResult, ...]:
    """Estimate metrics and their spread across replacement samples."""
    context = build_metric_context(data, multiplier)
    estimates = {
        name: evaluate_metric(name, data.calculated, data.experimental, context)
        for name in metric_names
    }
    samples = {name: [] for name in metric_names}
    rng = random.Random(seed)
    n_records = len(data.calculated)
    for _ in range(rounds):
        indices = rng.choices(range(n_records), k=n_records)
        sampled = _resample_rows(data, indices)
        sampled_context = build_metric_context(sampled, multiplier)
        for name in metric_names:
            value = evaluate_metric(
                name, sampled.calculated, sampled.experimental, sampled_context
            )
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
        f"Bootstrap rounds: {config.bootstrap_rounds}",
        f"Seed: {seed}",
        "Uncertainty method: paired records sampled with replacement",
        "",
        f"{'Metric':<8}{'Estimate':>16}{'Bootstrap SD':>16}{'Valid rounds':>14}",
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
