# QualStat Python

Two command-line tools for analysing alchemical binding free energies:

| Program | Purpose |
|---|---|
| `scripts/qualstat.py` | Compare calculated and experimental ABFE or RBFE values with legacy QualStat metrics, Spearman statistics, replacement bootstrap, and optional cycle closure. |
| `scripts/rbfe_to_abfe.py` | Reconstruct one mean-aligned binding free energy per ligand from a connected RBFE network. |

The project preserves the numerical definitions of the original
[QualStat](https://signe.teokem.lu.se/ulf/Methods/qual-stat.html) program where
compatibility matters. Legacy metrics are labelled explicitly in the detailed
documentation.

## Installation

### Standard Python environment

With Python 3.10 or newer, this route runs QualStat and the built-in NumPy
network estimator:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate with `.venv\Scripts\activate`.

### Conda/Mamba environment with Cinnabar

The supplied environment also installs OpenFreeEnergy Cinnabar 0.6.1:

```bash
mamba env create -f environment.yml
mamba activate qualstat
```

The equivalent Conda commands are:

```bash
conda env create -f environment.yml
conda activate qualstat
```

Do not run `pip install cinnabar`; that PyPI name belongs to an unrelated
package. The default `--backend auto` uses OpenFreeEnergy Cinnabar when it is
available and otherwise uses the included mathematically equivalent NumPy
weighted-least-squares implementation. Select `--backend cinnabar` or
`--backend numpy` when a workflow must require one implementation. The terminal
summary reports which estimator was used.

## Quick start

Create a fully commented QualStat configuration:

```bash
python scripts/qualstat.py --write-template
python scripts/qualstat.py qualstat_template.yaml
```

Template generation will not replace an existing file. Add `--force` only when
you intentionally want to overwrite it.

Reconstruct ligand energies from an RBFE network:

```bash
python scripts/rbfe_to_abfe.py network.csv experimental.csv \
  -o calculated_abfe.csv --cycle-closure 5
```

Paths in a QualStat YAML file are resolved relative to that YAML file.

## Input schemas

QualStat ABFE:

```text
ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty
```

QualStat RBFE:

```text
ligand_a,ligand_b,calculated,calculated_uncertainty,experimental,experimental_uncertainty
```

RBFE network reconstruction:

```text
ligand_A,ligand_B,DDG_kJ_mol,DDG_uncertainty_kJ_mol
```

with experimental values:

```text
ligand,DG_exp_kJ_mol
```

For both RBFE interfaces, the sign convention is

```text
DDG(A -> B) = G(B) - G(A)
```

## Uncertainty convention

Every uncertainty column must contain the **one-sigma uncertainty of the value
reported in the corresponding energy column**.

If the reported energy is the mean of `n` independent repeats and `s` is their
sample standard deviation, supply the standard error of the mean:

```text
SEM = s / sqrt(n)
```

For triplicates, this is `s / sqrt(3)`. If the value supplied by the simulation
or assay software is already the uncertainty of the reported mean, do not
divide it again.

QualStat samples complete CSV rows with replacement in each bootstrap round.
Calculated and experimental values from the same row stay paired. `Bootstrap SD`
is the sample standard deviation of the metric across valid rounds. The
uncertainty columns are used by metrics such as `taux` and `taurx` and by RBFE
cycle analysis; the bootstrap does not draw Gaussian noise from them.

`scripts/qualstat.py` handles configuration, input, resampling, and reports.
`scripts/qualstat_metrics.py` defines each metric in its own small function.
Run the entry point with `python scripts/qualstat.py <settings.yaml>`.

## Recommended RBFE metrics

The default RBFE template enables the four metrics commonly used by this
project:

- `MAD` — mean absolute error;
- `RMSD` — root mean square error;
- `taurx` — legacy sign agreement after excluding values that are not
  significant relative to their uncertainties;
- `r22` — legacy origin-symmetric signed squared Pearson statistic.

Ordinary Pearson `R`, Spearman `rho`, and legacy `tau` depend on the arbitrary
direction assigned to individual RBFE edges. They should not be reported for
raw RBFE edge sets.

`rho2` is an optional, orientation-independent extension. It computes standard
Spearman correlation after adding the sign-negated copy of every edge. This
makes it invariant to reversing any individual transformation while retaining
sign information. It is a custom project metric, not a standard named
statistic; report its definition whenever it is used.

## Documentation, examples, and tests

- [QualStat usage and metric definitions](docs/qualstat_docs.md)
- [RBFE-to-ABFE reconstruction and cycle closure](docs/rbfe_to_abfe_docs.md)
- [Runnable dummy datasets and expected outputs](examples/README.md)

Run the complete suite from the repository root:

```bash
python -m unittest discover -s tests -v
```

The tests verify formulas, orientation behavior, replacement sampling,
input validation, network reconstruction, cycle closure, CLI behavior, and the
committed examples. They do not establish physical convergence of a molecular
simulation.
