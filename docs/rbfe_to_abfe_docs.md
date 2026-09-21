# RBFE network reconstruction and cycle closure

`scripts/rbfe_to_abfe.py` fits one internally consistent free energy per ligand
from a connected RBFE network, then applies one common shift so that the mean
calculated energy equals the mean experimental energy.

These are reconstructed relative-network energies. They are **not** independent
ABFE simulations, and the experimental values supply the otherwise
undetermined overall offset.

## Installation and estimators

For the built-in NumPy estimator:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

For OpenFreeEnergy Cinnabar 0.6.1:

```bash
mamba env create -f environment.yml
mamba activate qualstat
```

Do not install the unrelated PyPI package named `cinnabar`.

The default `--backend auto` uses Cinnabar when its `FEMap` and `openff.units`
imports are available; otherwise it uses the included NumPy
weighted-least-squares implementation. Use `--backend cinnabar` or `--backend
numpy` to require one implementation. An explicitly requested backend never
silently falls back, and an installed-but-broken Cinnabar stack raises an
error. The terminal summary identifies the estimator used. Both implement the
Gaussian maximum-likelihood model documented by
[Cinnabar](https://cinnabar.openfree.energy/en/stable/concepts/estimators.html).

## Input files

### RBFE network

```csv
ligand_A,ligand_B,DDG_kJ_mol,DDG_uncertainty_kJ_mol
LIG1,LIG2,2.40,0.50
LIG2,LIG3,-1.80,0.60
LIG3,LIG1,-0.30,0.55
```

The direction convention is

```text
DDG(A -> B) = G(B) - G(A)
```

The graph must be connected. Self-transformations, non-finite values,
non-positive uncertainties, and repeated unordered pairs are rejected. Combine
repeats into one row before running the program.

### Experimental values

```csv
ligand,DG_exp_kJ_mol
LIG1,-25.20
LIG2,-22.90
LIG3,-27.10
```

Every network ligand must occur exactly once. Missing, extra, or duplicate
ligands are rejected.

All input and output energies are in kJ/mol.

## Uncertainty convention

`DDG_uncertainty_kJ_mol` is the **one-sigma uncertainty of the reported mean
RBFE**, because that value controls the statistical weight of its edge.

For `n` independent repeats with sample standard deviation `s`, use

```text
SEM = s / sqrt(n)
```

For triplicates, use `s / sqrt(3)`. Do not divide again if the upstream result
already reports the standard error or uncertainty of the mean. Shared
trajectories, common references, and other correlations are not represented by
this one-column independent-error model.

## Run the calculation

```bash
python scripts/rbfe_to_abfe.py network.csv experimental.csv \
  -o calculated_abfe.csv
```

To require Cinnabar, for example in a reproducibility workflow:

```bash
python scripts/rbfe_to_abfe.py network.csv experimental.csv \
  -o calculated_abfe.csv --backend cinnabar
```

Add cycle diagnostics for cycles containing at most five ligands:

```bash
python scripts/rbfe_to_abfe.py network.csv experimental.csv \
  -o calculated_abfe.csv --cycle-closure 5
```

The main output contains:

```text
ligand,DG_exp_kJ_mol,DG_calc_kJ_mol,DG_calc_uncertainty_kJ_mol
```

## Statistical model

For an edge from ligand `i` to ligand `j`, the model is

```text
DDG_ij ~ Normal(G_j - G_i, sigma_ij^2)
```

Writing all edge equations as `d = A G + error`, the fit minimizes

```text
sum_i ((d_i - (A G)_i) / sigma_i)^2
```

Thus smaller-uncertainty edges receive more weight. The information matrix is
`A.T W A`, with `W_ii = 1 / sigma_i^2`. Because adding a constant to every
ligand energy changes no relative difference, this matrix has one null mode.
The NumPy estimator uses its Moore-Penrose inverse, which selects the zero-mean
solution and supplies the fitted covariance. This reproduces Cinnabar's
zero-centred MLE values and uncertainties for the tested networks.

The script then applies

```text
shift = mean(DG_exp - DG_calc)
```

to all fitted energies. This shift changes neither pairwise differences nor
their ranking.

### What the output uncertainty includes

`DG_calc_uncertainty_kJ_mol` propagates the independent uncertainties assigned
to the RBFE edges through the network fit.

It does not include:

- experimental uncertainty or uncertainty in the mean-alignment shift;
- force-field or assay systematic error;
- correlated edge errors;
- pose, protonation, mapping, or sampling problems absent from the input
  uncertainties.

Because experiment supplies the common offset, MSD/bias of the mean-aligned
values is constrained and must not be interpreted as an independent ABFE
validation result.

## Cycle closure

For a cycle `A -> B -> C -> A`, the signed edge sum should be zero. The program
reports its absolute value:

```text
cc = abs(DDG_AB + DDG_BC + DDG_CA)
```

Assuming independent edges:

```text
cycle_uncertainty = sqrt(sigma_AB^2 + sigma_BC^2 + sigma_CA^2)
cc_unc_normalized = cc / cycle_uncertainty
cc_per_sqrt_edge = cc / sqrt(number_of_edges)
```

The cycle table contains:

| Column | Meaning |
|---|---|
| `cycle` | Canonical ligand sequence for the unique simple cycle. |
| `cc_kJ_mol` | Absolute closure error. |
| `cc_per_sqrt_edge_kJ_mol` | Closure error scaled by cycle length. |
| `cc_unc_normalized` | Closure error divided by combined edge uncertainty. |

The companion edge table reports how many assessed cycles contain each edge,
plus the mean and maximum scaled closure errors for those cycles.

A normalized value below one means only that the closure error is smaller than
the supplied combined uncertainty. It does not prove convergence; errors can
cancel, and the supplied uncertainties may themselves be inaccurate. Large
values identify cycles and shared edges to inspect, not a uniquely guilty edge.

Cycle diagnostics do not alter or reject fitted ligand energies.

## Validation and output safety

The script validates all inputs and completes every requested calculation
before replacing output files. If cycle analysis fails, an existing main output
is preserved rather than being partially overwritten.

Empty cycle results still produce CSV files with stable headers. Simple-cycle
enumeration can grow rapidly for dense networks, so a maximum length of four or
five is a practical starting point.

## Example and tests

The directory `examples/rbfe_to_abfe/` contains a connected four-ligand network,
mean-aligned output, and both cycle tables.

```bash
python -m unittest discover -s tests -v
```

The tests cover sign conventions, analytical reconstruction cases, uncertainty
scaling, arbitrary edge reversals, Cinnabar-reference output, cycle formulas,
validation, and atomic file handling. The Conda CI job preflights Cinnabar and
runs a backend-specific integration test, so it cannot pass by silently using
NumPy instead.

Further reading:

- [Cinnabar estimators](https://cinnabar.openfree.energy/en/stable/concepts/estimators.html)
- [Cinnabar cycle-closure tutorial](https://cinnabar.openfree.energy/en/latest/tutorials/cycle-closure.html)
- [Binding-affinity benchmark best practices](https://pubmed.ncbi.nlm.nih.gov/36382113/)
