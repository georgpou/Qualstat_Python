#!/usr/bin/env python3
"""Reconstruct mean-aligned ligand binding free energies from an RBFE network."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

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


def _require_columns(data: pd.DataFrame, required: set[str], name: str) -> None:
    """Fail early when an input CSV is missing required columns."""
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"{name} CSV is missing columns: {', '.join(sorted(missing))}")


def _validate_inputs(network: pd.DataFrame, experimental: pd.DataFrame) -> None:
    """Check the inputs before passing them to Cinnabar."""
    _validate_network(network)
    _require_columns(experimental, EXPERIMENTAL_COLUMNS, "Experimental")

    if experimental.empty:
        raise ValueError("Experimental CSV must not be empty")
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
    """Check the RBFE network before passing it to Cinnabar."""
    _require_columns(network, NETWORK_COLUMNS, "Network")

    if network.empty:
        raise ValueError("Network CSV must not be empty")
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


def _build_femap(network: pd.DataFrame):
    """Build and connectivity-check a Cinnabar FEMap."""
    try:
        from cinnabar.femap import FEMap
        from openff.units import unit
    except ImportError as exc:
        raise RuntimeError(
            "OpenFreeEnergy Cinnabar is required. Install it from conda-forge; "
            "do not use the unrelated PyPI package named 'cinnabar'."
        ) from exc

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


def cycle_closure_diagnostics(
    network: pd.DataFrame, max_cycle_length: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return cycle- and edge-level closure diagnostics in kJ/mol."""
    _validate_network(network)
    if max_cycle_length < 3:
        raise ValueError("Maximum cycle length must be at least 3")

    # Cinnabar's cycle tables require one value for each undirected edge.
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
            "Cycle closure requires one aggregated edge per unordered ligand pair; "
            f"repeated pairs: {', '.join(sorted(labels))}"
        )

    femap = _build_femap(network)
    cycles = femap.get_cycle_closure_dataframe(max_cycle_length=max_cycle_length)
    if cycles.empty:
        # Cinnabar 0.6.1 raises a KeyError for edge statistics with no cycles.
        edges = pd.DataFrame(
            columns=[
                "source",
                "ligandA",
                "ligandB",
                "n_cycles",
                "mean_cc_per_edge (kcal/mol)",
                "max_cc_per_edge (kcal/mol)",
            ]
        )
    else:
        edges = femap.get_cycle_closure_edge_statistics_dataframe(
            max_cycle_length=max_cycle_length
        )

    # Cinnabar reports dataframes in kcal/mol; normalized closure is unitless.
    cycles["cc_kJ_mol"] = cycles.pop("cc (kcal/mol)") * KCAL_TO_KJ
    cycles["cc_per_sqrt_edge_kJ_mol"] = (
        cycles.pop("cc_per_edge (kcal/mol)") * KCAL_TO_KJ
    )
    edges["mean_cc_per_sqrt_edge_kJ_mol"] = (
        edges.pop("mean_cc_per_edge (kcal/mol)") * KCAL_TO_KJ
    )
    edges["max_cc_per_sqrt_edge_kJ_mol"] = (
        edges.pop("max_cc_per_edge (kcal/mol)") * KCAL_TO_KJ
    )
    return cycles, edges


def analyze(network: pd.DataFrame, experimental: pd.DataFrame) -> pd.DataFrame:
    """Run Cinnabar MLE and align the calculated mean to the experimental mean."""
    _validate_inputs(network, experimental)
    femap = _build_femap(network)

    # MLE combines all edges and cycles into one internally consistent solution.
    femap.generate_absolute_values()
    calculated = femap.get_absolute_dataframe()
    calculated = calculated[calculated["computational"]].copy()

    # Cinnabar's dataframe API reports kcal/mol, so convert back to kJ/mol.
    calculated = calculated.rename(columns={"label": "ligand"})
    calculated["DG_calc_kJ_mol"] = calculated["DG (kcal/mol)"] * KCAL_TO_KJ
    calculated["DG_calc_uncertainty_kJ_mol"] = (
        calculated["uncertainty (kcal/mol)"] * KCAL_TO_KJ
    )

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


def main() -> None:
    """Read CSV inputs, run the analysis, and write the requested CSV files."""
    parser = argparse.ArgumentParser(
        description="Reconstruct and mean-align ligand free energies with Cinnabar."
    )
    parser.add_argument("network_csv", type=Path, help="RBFE edge CSV in kJ/mol")
    parser.add_argument("experimental_csv", type=Path, help="Experimental DG CSV in kJ/mol")
    parser.add_argument("-o", "--output", type=Path, default=Path("calculated_abfe.csv"))
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
        network = pd.read_csv(args.network_csv)
        experimental = pd.read_csv(args.experimental_csv)
        result = analyze(network, experimental)
        outputs = [(result, args.output)]

        if args.cycle_closure is not None:
            cycles, cycle_edges = cycle_closure_diagnostics(
                network, max_cycle_length=args.cycle_closure
            )
            cycle_path = args.output.with_name(
                f"{args.output.stem}_cycle_closure.csv"
            )
            edge_path = args.output.with_name(
                f"{args.output.stem}_cycle_closure_edges.csv"
            )
            outputs.extend([(cycles, cycle_path), (cycle_edges, edge_path)])

        # Do not touch output files until every requested calculation has succeeded.
        _write_csvs_atomically(outputs)
    except Exception as exc:
        parser.error(str(exc))

    print(f"Ligands: {len(result)}")
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
