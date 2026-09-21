#!/usr/bin/env python3
"""Reconstruct mean-aligned ligand binding free energies from an RBFE network."""

from __future__ import annotations

import argparse
from importlib.util import find_spec
import math
import os
import tempfile
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


NETWORK_COLUMNS = {
    "ligand_A",
    "ligand_B",
    "DDG_kJ_mol",
    "DDG_uncertainty_kJ_mol",
}
EXPERIMENTAL_COLUMNS = {"ligand", "DG_exp_kJ_mol"}
KCAL_TO_KJ = 4.184
CYCLE_COLUMNS = [
    "source",
    "cycle",
    "cc_unc_normalized",
    "cc_kJ_mol",
    "cc_per_sqrt_edge_kJ_mol",
]
CYCLE_EDGE_COLUMNS = [
    "source",
    "ligandA",
    "ligandB",
    "n_cycles",
    "mean_cc_per_sqrt_edge_kJ_mol",
    "max_cc_per_sqrt_edge_kJ_mol",
]


class CinnabarUnavailableError(RuntimeError):
    """Raised only when the optional OpenFreeEnergy backend cannot be imported."""


def _require_columns(data: pd.DataFrame, required: set[str], name: str) -> None:
    """Fail early when an input CSV is missing required columns."""
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"{name} CSV is missing columns: {', '.join(sorted(missing))}")


def _validate_inputs(network: pd.DataFrame, experimental: pd.DataFrame) -> None:
    """Check reconstruction and alignment inputs."""
    _validate_network(network)
    _require_columns(experimental, EXPERIMENTAL_COLUMNS, "Experimental")

    if experimental.empty:
        raise ValueError("Experimental CSV must not be empty")
    if any(pd.isna(value) or not str(value).strip() for value in experimental["ligand"]):
        raise ValueError("Experimental ligand names must be non-empty")
    if experimental["ligand"].duplicated().any():
        raise ValueError("Experimental ligand names must be unique")
    if not np.isfinite(pd.to_numeric(experimental["DG_exp_kJ_mol"], errors="coerce")).all():
        raise ValueError("Experimental energies must contain finite numbers")

    # Mean alignment requires one experimental value for every network ligand.
    network_ligands = set(network["ligand_A"]) | set(network["ligand_B"])
    experimental_ligands = set(experimental["ligand"])
    if network_ligands != experimental_ligands:
        missing_exp = sorted(network_ligands - experimental_ligands)
        extra_exp = sorted(experimental_ligands - network_ligands)
        raise ValueError(
            "Ligand mismatch between CSV files; "
            f"missing experimental={missing_exp}, extra experimental={extra_exp}"
        )


def _validate_network(network: pd.DataFrame) -> None:
    """Check the RBFE network before fitting or cycle analysis."""
    _require_columns(network, NETWORK_COLUMNS, "Network")

    if network.empty:
        raise ValueError("Network CSV must not be empty")
    if any(
        pd.isna(value) or not str(value).strip()
        for column in (network["ligand_A"], network["ligand_B"])
        for value in column
    ):
        raise ValueError("Network ligand names must be non-empty")
    if (network["ligand_A"] == network["ligand_B"]).any():
        raise ValueError("Self-transformations are not allowed")

    # Energies must be finite; edge uncertainties must also be positive.
    numeric_columns = [
        network["DDG_kJ_mol"],
        network["DDG_uncertainty_kJ_mol"],
    ]
    if any(not np.isfinite(pd.to_numeric(column, errors="coerce")).all() for column in numeric_columns):
        raise ValueError("Energy and uncertainty columns must contain finite numbers")
    if (pd.to_numeric(network["DDG_uncertainty_kJ_mol"]) <= 0).any():
        raise ValueError("Every RBFE uncertainty must be greater than zero")

    seen_pairs: set[frozenset[object]] = set()
    repeated_pairs: set[frozenset[object]] = set()
    for row in network.itertuples(index=False):
        pair = frozenset((row.ligand_A, row.ligand_B))
        if pair in seen_pairs:
            repeated_pairs.add(pair)
        seen_pairs.add(pair)
    if repeated_pairs:
        labels = ["-".join(sorted(map(str, pair))) for pair in repeated_pairs]
        raise ValueError(
            "Use one aggregated edge per unordered ligand pair; "
            f"repeated pairs: {', '.join(sorted(labels))}"
        )


def _network_ligands(network: pd.DataFrame) -> list[object]:
    """Return stable, deterministically ordered ligand labels."""
    return sorted(set(network["ligand_A"]) | set(network["ligand_B"]), key=str)


def _check_connected(network: pd.DataFrame) -> None:
    """Require one weakly connected component without relying on Cinnabar."""
    ligands = _network_ligands(network)
    adjacency = {ligand: set() for ligand in ligands}
    for row in network.itertuples(index=False):
        adjacency[row.ligand_A].add(row.ligand_B)
        adjacency[row.ligand_B].add(row.ligand_A)

    visited = set()
    pending = [ligands[0]]
    while pending:
        ligand = pending.pop()
        if ligand in visited:
            continue
        visited.add(ligand)
        pending.extend(adjacency[ligand] - visited)
    if len(visited) != len(ligands):
        raise ValueError("The RBFE network is disconnected")


def _build_femap(network: pd.DataFrame):
    """Build and connectivity-check a Cinnabar FEMap."""
    if find_spec("cinnabar") is None or find_spec("openff") is None:
        raise CinnabarUnavailableError(
            "OpenFreeEnergy Cinnabar is required. Install it from conda-forge; "
            "do not use the unrelated PyPI package named 'cinnabar'."
        )

    # Once the optional packages are present, import/runtime errors indicate a
    # broken installation and must not be disguised as an unavailable backend.
    from cinnabar.femap import FEMap
    from openff.units import unit

    # Build the directed RBFE network in the units supplied by the CSV.
    femap = FEMap()
    for row in network.itertuples(index=False):
        femap.add_relative_calculation(
            labelA=row.ligand_A,
            labelB=row.ligand_B,
            value=float(row.DDG_kJ_mol) * unit.kilojoule_per_mole,
            uncertainty=float(row.DDG_uncertainty_kJ_mol) * unit.kilojoule_per_mole,
            source="ATM-RBFE",
        )

    # Cinnabar cannot reconstruct a common scale from disconnected components.
    if not femap.check_weakly_connected():
        raise ValueError("The RBFE network is disconnected")
    return femap


def _fit_network_numpy(network: pd.DataFrame) -> pd.DataFrame:
    """Fit zero-mean ligand energies by Gaussian weighted least squares."""
    ligands = _network_ligands(network)
    ligand_index = {ligand: index for index, ligand in enumerate(ligands)}
    design = np.zeros((len(network), len(ligands)), dtype=float)
    for row_index, row in enumerate(network.itertuples(index=False)):
        design[row_index, ligand_index[row.ligand_A]] = -1.0
        design[row_index, ligand_index[row.ligand_B]] = 1.0

    differences = network["DDG_kJ_mol"].to_numpy(dtype=float)
    uncertainties = network["DDG_uncertainty_kJ_mol"].to_numpy(dtype=float)
    weights = 1.0 / uncertainties**2
    information = design.T @ (weights[:, None] * design)
    rhs = design.T @ (weights * differences)

    # A connected difference network has one null mode: adding the same
    # constant to every ligand. The Moore-Penrose inverse fixes the zero-mean
    # gauge and gives the covariance in that same gauge.
    covariance = np.linalg.pinv(information, hermitian=True)
    energies = covariance @ rhs
    fitted_uncertainties = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    return pd.DataFrame(
        {
            "ligand": ligands,
            "DG_calc_kJ_mol": energies,
            "DG_calc_uncertainty_kJ_mol": fitted_uncertainties,
        }
    )


def _fit_network_cinnabar(network: pd.DataFrame) -> pd.DataFrame:
    """Fit ligand energies with the OpenFreeEnergy Cinnabar MLE backend."""
    femap = _build_femap(network)
    femap.generate_absolute_values()
    calculated = femap.get_absolute_dataframe()
    calculated = calculated[calculated["computational"]].copy()

    # Cinnabar's dataframe API reports kcal/mol, so convert back to kJ/mol.
    calculated = calculated.rename(columns={"label": "ligand"})
    calculated["DG_calc_kJ_mol"] = calculated["DG (kcal/mol)"] * KCAL_TO_KJ
    calculated["DG_calc_uncertainty_kJ_mol"] = (
        calculated["uncertainty (kcal/mol)"] * KCAL_TO_KJ
    )
    return calculated


def _canonical_cycle(nodes: Sequence[object]) -> tuple[object, ...]:
    """Canonicalize an undirected simple cycle up to rotation and reversal."""
    sequence = tuple(nodes)
    candidates = [
        sequence[index:] + sequence[:index] for index in range(len(sequence))
    ]
    reversed_sequence = tuple(reversed(sequence))
    candidates.extend(
        reversed_sequence[index:] + reversed_sequence[:index]
        for index in range(len(reversed_sequence))
    )
    return min(candidates, key=lambda cycle: tuple(map(str, cycle)))


def _find_simple_cycles(
    network: pd.DataFrame, max_cycle_length: int
) -> tuple[tuple[object, ...], ...]:
    """Enumerate unique undirected simple cycles within the requested limit."""
    ligands = _network_ligands(network)
    adjacency = {ligand: set() for ligand in ligands}
    for row in network.itertuples(index=False):
        adjacency[row.ligand_A].add(row.ligand_B)
        adjacency[row.ligand_B].add(row.ligand_A)

    cycles: set[tuple[object, ...]] = set()
    for start in ligands:
        path = [start]
        visited = {start}
        stack = [(start, iter(sorted(adjacency[start], key=str)))]
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
                if 3 <= len(path) <= max_cycle_length:
                    cycles.add(_canonical_cycle(path))
                continue
            if neighbor in visited or len(path) >= max_cycle_length:
                continue
            path.append(neighbor)
            visited.add(neighbor)
            stack.append((neighbor, iter(sorted(adjacency[neighbor], key=str))))
    return tuple(sorted(cycles, key=lambda cycle: (len(cycle), tuple(map(str, cycle)))))


def _cycle_closure_numpy(
    network: pd.DataFrame, max_cycle_length: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate Cinnabar-compatible cycle and edge summaries in kJ/mol."""
    edge_data = {}
    for row in network.itertuples(index=False):
        edge_data[frozenset((row.ligand_A, row.ligand_B))] = (
            row.ligand_A,
            row.ligand_B,
            float(row.DDG_kJ_mol),
            float(row.DDG_uncertainty_kJ_mol),
        )

    cycle_rows = []
    edge_values: dict[frozenset[object], list[float]] = {
        pair: [] for pair in edge_data
    }
    for cycle in _find_simple_cycles(network, max_cycle_length):
        signed_values = []
        uncertainties = []
        cycle_pairs = []
        for source, target in zip(cycle, cycle[1:] + cycle[:1]):
            pair = frozenset((source, target))
            stored_source, stored_target, value, uncertainty = edge_data[pair]
            signed_values.append(
                value if (source, target) == (stored_source, stored_target) else -value
            )
            uncertainties.append(uncertainty)
            cycle_pairs.append(pair)
        closure = abs(math.fsum(signed_values))
        propagated_uncertainty = math.sqrt(
            math.fsum(value * value for value in uncertainties)
        )
        scaled_closure = closure / math.sqrt(len(cycle))
        cycle_rows.append(
            {
                "source": "ATM-RBFE",
                "cycle": cycle,
                "cc_unc_normalized": closure / propagated_uncertainty,
                "cc_kJ_mol": closure,
                "cc_per_sqrt_edge_kJ_mol": scaled_closure,
            }
        )
        for pair in cycle_pairs:
            edge_values[pair].append(scaled_closure)

    cycles = pd.DataFrame(cycle_rows, columns=CYCLE_COLUMNS)
    edge_rows = []
    for pair, values in edge_values.items():
        if not values:
            continue
        ligand_a, ligand_b = sorted(pair, key=str)
        edge_rows.append(
            {
                "source": "ATM-RBFE",
                "ligandA": ligand_a,
                "ligandB": ligand_b,
                "n_cycles": len(values),
                "mean_cc_per_sqrt_edge_kJ_mol": math.fsum(values) / len(values),
                "max_cc_per_sqrt_edge_kJ_mol": max(values),
            }
        )
    edges = pd.DataFrame(edge_rows, columns=CYCLE_EDGE_COLUMNS)
    if not edges.empty:
        edges = edges.sort_values(["ligandA", "ligandB"], kind="stable").reset_index(drop=True)
    return cycles, edges


def cycle_closure_diagnostics(
    network: pd.DataFrame, max_cycle_length: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return cycle- and edge-level closure diagnostics in kJ/mol."""
    _validate_network(network)
    if max_cycle_length < 3:
        raise ValueError("Maximum cycle length must be at least 3")

    return _cycle_closure_numpy(network, max_cycle_length)


def analyze(
    network: pd.DataFrame, experimental: pd.DataFrame, backend: str = "auto"
) -> pd.DataFrame:
    """Fit network energies and align their mean to the experimental mean."""
    _validate_inputs(network, experimental)
    _check_connected(network)
    if backend not in {"auto", "numpy", "cinnabar"}:
        raise ValueError("backend must be 'auto', 'numpy', or 'cinnabar'")

    if backend == "numpy":
        calculated = _fit_network_numpy(network)
        estimator = "NumPy weighted least squares"
    elif backend == "cinnabar":
        calculated = _fit_network_cinnabar(network)
        estimator = "OpenFreeEnergy Cinnabar"
    else:
        try:
            calculated = _fit_network_cinnabar(network)
        except CinnabarUnavailableError:
            calculated = _fit_network_numpy(network)
            estimator = "NumPy weighted least squares"
        else:
            estimator = "OpenFreeEnergy Cinnabar"

    # Join by ligand name, then fit the single undetermined additive constant.
    result = experimental[["ligand", "DG_exp_kJ_mol"]].merge(
        calculated[["ligand", "DG_calc_kJ_mol", "DG_calc_uncertainty_kJ_mol"]],
        on="ligand",
        validate="one_to_one",
    )
    shift = (result["DG_exp_kJ_mol"] - result["DG_calc_kJ_mol"]).mean()
    result["DG_calc_kJ_mol"] += shift

    # Stable ordering makes output comparisons reproducible.
    result = result.sort_values("ligand", kind="stable").reset_index(drop=True)
    result.attrs["mean_shift_kJ_mol"] = float(shift)
    result.attrs["estimator"] = estimator
    return result


def _write_csvs_atomically(outputs: list[tuple[pd.DataFrame, Path]]) -> None:
    """Stage every CSV before replacing any requested output file."""
    staged: list[tuple[Path, Path]] = []
    try:
        for data, destination in outputs:
            destination = destination.resolve()
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="",
                prefix=f".{destination.name}.",
                suffix=".tmp",
                dir=destination.parent,
                delete=False,
            ) as handle:
                data.to_csv(handle, index=False, float_format="%.6f")
                staged.append((Path(handle.name), destination))

        # Each replacement is atomic, and calculation/serialization finishes first.
        for temporary, destination in staged:
            os.replace(temporary, destination)
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)


def _same_file(first: Path, second: Path) -> bool:
    """Compare paths safely even when a requested output does not exist yet."""
    first = first.resolve()
    second = second.resolve()
    if first == second:
        return True
    try:
        return first.samefile(second)
    except OSError:
        return False


def _validate_output_paths(inputs: Sequence[Path], outputs: Sequence[Path]) -> None:
    """Prevent any main or derived output from replacing an input CSV."""
    if any(_same_file(output, input_path) for output in outputs for input_path in inputs):
        raise ValueError("Output files must differ from both input CSV files")


def main() -> None:
    """Read CSV inputs, run the analysis, and write the requested CSV files."""
    parser = argparse.ArgumentParser(
        description="Reconstruct and mean-align ligand free energies from an RBFE network."
    )
    parser.add_argument("network_csv", type=Path, help="RBFE edge CSV in kJ/mol")
    parser.add_argument("experimental_csv", type=Path, help="Experimental DG CSV in kJ/mol")
    parser.add_argument("-o", "--output", type=Path, default=Path("calculated_abfe.csv"))
    parser.add_argument(
        "--backend",
        choices=("auto", "numpy", "cinnabar"),
        default="auto",
        help="network estimator (default: auto; prefer Cinnabar when installed)",
    )
    parser.add_argument(
        "--cycle-closure",
        type=int,
        metavar="N",
        help="Assess all cycles containing at most N ligands (N must be at least 3)",
    )
    args = parser.parse_args()

    try:
        if args.cycle_closure is not None and args.cycle_closure < 3:
            raise ValueError("--cycle-closure must be at least 3")
        cycle_path = args.output.with_name(f"{args.output.stem}_cycle_closure.csv")
        edge_path = args.output.with_name(
            f"{args.output.stem}_cycle_closure_edges.csv"
        )
        requested_paths = [args.output]
        if args.cycle_closure is not None:
            requested_paths.extend((cycle_path, edge_path))
        _validate_output_paths(
            (args.network_csv, args.experimental_csv), requested_paths
        )
        network = pd.read_csv(args.network_csv)
        experimental = pd.read_csv(args.experimental_csv)
        result = analyze(network, experimental, backend=args.backend)
        outputs = [(result, args.output)]

        if args.cycle_closure is not None:
            cycles, cycle_edges = cycle_closure_diagnostics(
                network, max_cycle_length=args.cycle_closure
            )
            outputs.extend([(cycles, cycle_path), (cycle_edges, edge_path)])

        # Do not touch output files until every requested calculation has succeeded.
        _write_csvs_atomically(outputs)
    except Exception as exc:
        parser.error(str(exc))

    print(f"Ligands: {len(result)}")
    print(f"Estimator: {result.attrs['estimator']}")
    print(f"Mean alignment shift: {result.attrs['mean_shift_kJ_mol']:.6f} kJ/mol")
    print(f"Output: {args.output}")
    if args.cycle_closure is not None:
        print(f"Cycles detected: {len(cycles)}")
        if cycles.empty:
            print(f"No cycles containing 3-{args.cycle_closure} ligands were found")
        else:
            largest = cycles["cc_unc_normalized"].max()
            print(f"Largest normalized closure error: {largest:.6f}")
        print(f"Cycle output: {cycle_path}")
        print(f"Cycle-edge output: {edge_path}")


if __name__ == "__main__":
    main()
