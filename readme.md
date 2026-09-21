# QualStat analysis tools

This repository contains two small command-line programs for analysing alchemical free-energy results:

- `scripts/qualstat.py` calculates statistical quality metrics for calculated versus experimental ABFE or RBFE data. It supports a parametric bootstrap, a collection of legacy QualStat metrics, Spearman correlation, and optional RBFE cycle-closure analysis.
- `scripts/rbfe_to_abfe.py` uses OpenFreeEnergy Cinnabar to reconstruct one relative-energy-derived binding free energy per ligand from a connected RBFE network, aligns the calculated energies to the experimental mean, and can generate Cinnabar cycle-closure diagnostics.

More detailed documentation is in:

- `docs/qualstat_docs.md`
- `docs/rbfe_to_abfe_docs.md`

## Repository layout

```text
.
├── .gitignore
├── environment.yml
├── readme.md
├── docs/
│   ├── qualstat_docs.md
│   └── rbfe_to_abfe_docs.md
├── examples/
│   ├── README.md
│   ├── qualstat_abfe/
│   ├── qualstat_rbfe/
│   └── rbfe_to_abfe/
├── tests/
│   └── ...
└── scripts/
    ├── qualstat.py
    └── rbfe_to_abfe.py
```

`environment.yml` is included in addition to the requested source/docs layout so that the shared Conda environment is reproducible.

## Environment setup

The shared environment is called `qualstat` and uses only the `conda-forge` channel.

Using Mamba:

```bash
mamba env create -f environment.yml
mamba activate qualstat
```

Using Conda:

```bash
conda env create -f environment.yml
conda activate qualstat
```

To update an existing environment after `environment.yml` changes:

```bash
mamba env update -n qualstat -f environment.yml --prune
```

or:

```bash
conda env update -n qualstat -f environment.yml --prune
```

### Why these versions are constrained

The two programs have very different dependency footprints.

`qualstat.py` uses the Python standard library plus PyYAML. The main Python requirement is therefore a modern Python version; Python 3.11 is used here as a conservative common version.

`rbfe_to_abfe.py` directly uses NumPy and pandas and imports `FEMap` from OpenFreeEnergy Cinnabar together with `openff.units`. Cinnabar is pinned to `0.6.1` because the script relies on the FEMap API and on the names/structure of the cycle-closure dataframes used by that release. NumPy and pandas are constrained to current major-version ranges to reduce the chance that a future major release changes behaviour underneath this fixed Cinnabar version.

Do not install `cinnabar` with `pip install cinnabar`. The package with that name on PyPI is not the OpenFreeEnergy Cinnabar used by this script. Install Cinnabar from `conda-forge`, as done by `environment.yml`.

`openff-units` is listed explicitly even though it may also be installed as a dependency of the Cinnabar stack, because `rbfe_to_abfe.py` imports it directly. This makes the repository's direct requirements clear and prevents accidental reliance on an undeclared transitive dependency.

## Quick checks after installation

From the repository root:

```bash
python scripts/qualstat.py --help
python scripts/rbfe_to_abfe.py --help
```

You can also confirm the important imports and versions:

```bash
python - <<'PY'
import sys
import yaml
import numpy
import pandas
import cinnabar
import openff.units

print("Python:", sys.version.split()[0])
print("PyYAML:", yaml.__version__)
print("NumPy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("Cinnabar:", cinnabar.__version__)
print("openff-units: import OK")
PY
```

## Quick use: QualStat

Create a commented YAML settings template:

```bash
python scripts/qualstat.py --write-template
```

Then edit `qualstat_template.yaml` and run:

```bash
python scripts/qualstat.py qualstat_template.yaml
```

The YAML selects ABFE or RBFE analysis, the input CSV, output report, bootstrap settings, requested statistics, and optional cycle analysis. Paths inside the YAML are resolved relative to the YAML file.

For the exact CSV columns, metric definitions, bootstrap procedure, and QualStat cycle-closure behaviour, see `docs/qualstat_docs.md`.

## Quick use: RBFE network reconstruction

Basic reconstruction:

```bash
python scripts/rbfe_to_abfe.py network.csv experimental.csv \
    -o calculated_abfe.csv
```

With cycle-closure diagnostics for cycles containing at most five ligands:

```bash
python scripts/rbfe_to_abfe.py network.csv experimental.csv \
    -o calculated_abfe.csv \
    --cycle-closure 5
```

The script fits one energy per ligand using Cinnabar's maximum-likelihood network estimator, then applies one common shift so that the calculated mean matches the experimental mean. The optional cycle analysis is diagnostic; it does not alter or reject fitted ligand energies.

For the required CSV schemas, sign convention, interpretation of the reconstructed energies and uncertainties, and cycle-closure outputs, see `docs/rbfe_to_abfe_docs.md`.

## Examples and tests

The `examples/README.md` file describes three small, reproducible workflows:
QualStat ABFE, QualStat RBFE with cycle analysis, and Cinnabar RBFE-to-ABFE
reconstruction. Run them from the repository root with the commands shown in
that file.

Run the complete test suite with:

```bash
python -m unittest discover -s tests -v
```

The tests check metric definitions, bootstrap reproducibility, input validation,
CLI behavior, Cinnabar reconstruction, cycle diagnostics, output safety, and
the committed examples. They do not prove that a molecular simulation is
physically converged or that a force field is scientifically correct. Those
questions still require independent replicas, sampling checks, setup review,
and scientific interpretation.

## Important: the two CSV interfaces are not identical

The programs can live in the same environment, but their file formats are currently independent.

`rbfe_to_abfe.py` expects an RBFE network with columns:

```text
ligand_A,ligand_B,DDG_kJ_mol,DDG_uncertainty_kJ_mol
```

and an experimental file with:

```text
ligand,DG_exp_kJ_mol
```

Its reconstructed output contains:

```text
ligand,DG_exp_kJ_mol,DG_calc_kJ_mol,DG_calc_uncertainty_kJ_mol
```

QualStat ABFE mode instead requires:

```text
ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty
```

QualStat RBFE mode requires:

```text
ligand_a,ligand_b,calculated,calculated_uncertainty,experimental,experimental_uncertainty
```

Therefore the Cinnabar reconstruction output cannot currently be passed directly to QualStat without a small conversion step. The key missing quantity is `experimental_uncertainty`: `rbfe_to_abfe.py` does not read or emit it, while QualStat requires it. Do not fill this column with zero unless zero genuinely represents the intended experimental standard deviation; supply the experimental uncertainty used for the statistical analysis.

## Dependency and compatibility notes

The main dependency risks are:

1. **Wrong Cinnabar package.** PyPI's `cinnabar` name refers to an unrelated package. Use the conda-forge package.
2. **Cinnabar API drift.** `rbfe_to_abfe.py` expects Cinnabar's `FEMap`, `generate_absolute_values()`, `get_absolute_dataframe()`, `get_cycle_closure_dataframe()`, and `get_cycle_closure_edge_statistics_dataframe()` interfaces and specific dataframe column names. Keeping `cinnabar=0.6.1` avoids an unnoticed API change.
3. **Future pandas/NumPy major releases.** The script itself uses ordinary NumPy/pandas operations, but Cinnabar also uses these libraries internally. The environment therefore avoids unbounded future major versions while retaining normal compatible updates within the selected ranges.
4. **Implicit OpenFF dependency.** `openff.units` is imported directly by the script, so it is declared explicitly in `environment.yml` rather than relying only on Cinnabar to pull it in.
5. **Different data schemas.** This is not a package conflict, but it is the main workflow incompatibility between the two programs. See the section above.

## Reproducibility

For a calculation used in a benchmark or publication, save the exact solved environment in addition to `environment.yml`:

```bash
conda env export -n qualstat --from-history > environment-history.yml
conda list -n qualstat --explicit > environment-explicit.txt
```

`environment.yml` records the intended direct dependencies. `environment-explicit.txt` records the exact package builds used on that machine and is useful when reproducing an analysis later.
