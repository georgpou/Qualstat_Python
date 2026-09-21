# RBFE network to mean-aligned binding free energies

## Purpose

`rbfe_to_abfe.py` converts a connected network of relative binding free-energy
(RBFE) results into one calculated binding free energy for each ligand. It uses
the maximum-likelihood estimator in
[Cinnabar](https://cinnabar.openfree.energy/en/latest/) and then moves the
calculated values onto the experimental energy scale by matching their mean.

The script also has an optional cycle-closure analysis. This checks whether
different calculation paths through the ligand network agree with each other.

All input and output energies are in **kJ/mol**.

> Important: these are binding free energies reconstructed from an RBFE
> network. They are not independent absolute binding free-energy simulations.
> Experimental results supply the one overall offset that RBFE calculations
> cannot determine by themselves.

## Installation

Install the OpenFreeEnergy version of Cinnabar from conda-forge:

```bash
mamba create -n cinnabar -c conda-forge python=3.11 cinnabar=0.6.1 pandas
mamba activate cinnabar
```

Do not use `pip install cinnabar`. The package with that name on PyPI is an
unrelated program and does not contain the required `FEMap` class.

## Input files

### 1. RBFE network

The network CSV must contain these four columns:

```csv
ligand_A,ligand_B,DDG_kJ_mol,DDG_uncertainty_kJ_mol
LIG1,LIG2,2.40,0.50
LIG2,LIG3,-1.80,0.60
LIG3,LIG1,-0.30,0.55
```

| Column | Meaning |
|---|---|
| `ligand_A` | Ligand at the beginning of the transformation |
| `ligand_B` | Ligand at the end of the transformation |
| `DDG_kJ_mol` | Calculated RBFE difference in kJ/mol |
| `DDG_uncertainty_kJ_mol` | Uncertainty of that RBFE result in kJ/mol |

The direction convention is:

```text
DDG(A -> B) = G(B) - G(A)
```

A sign mistake will produce a valid-looking but incorrect result, so check that
the ligand order agrees with the convention used by the program that generated
the RBFE values.

When three independent repeats are used, put one combined row in this file for
each ligand pair. Use the mean RBFE result as `DDG_kJ_mol` and the standard
deviation across the repeats as `DDG_uncertainty_kJ_mol`.

For cycle closure, do not include the same ligand pair more than once, even in
the opposite direction. For example, `LIG1 -> LIG2` and `LIG2 -> LIG1` count as
the same pair. Combine repeats before running this script. This rule prevents
Cinnabar's network fit and cycle table from silently using different versions
of the same edge.

### 2. Experimental binding free energies

The experimental CSV must contain:

```csv
ligand,DG_exp_kJ_mol
LIG1,-25.20
LIG2,-22.90
LIG3,-27.10
```

Every ligand in the RBFE network must appear exactly once in this file. Extra
or missing ligand names are treated as errors.

## Basic use

```bash
python rbfe_to_abfe.py network.csv experimental.csv \
    -o calculated_abfe.csv
```

The output contains:

```csv
ligand,DG_exp_kJ_mol,DG_calc_kJ_mol,DG_calc_uncertainty_kJ_mol
```

| Column | Meaning |
|---|---|
| `ligand` | Ligand name |
| `DG_exp_kJ_mol` | Experimental binding free energy |
| `DG_calc_kJ_mol` | Mean-aligned calculated binding free energy |
| `DG_calc_uncertainty_kJ_mol` | Network-based uncertainty from Cinnabar |

## What the calculation does

### Step 1: Build the ligand network

Ligands are the nodes of the network. Each calculated transformation is an
edge connecting two ligands. The complete network must be connected: it must
be possible to reach every ligand by following one or more calculated edges.

### Step 2: Fit one energy per ligand

Each RBFE result describes a difference between two unknown ligand energies:

```text
calculated DDG(A -> B) approximately equals G(B) - G(A)
```

Cinnabar finds the set of ligand energies that best explains all edges at the
same time. Edges with smaller reported uncertainties receive more weight.
Edges with larger uncertainties receive less weight.

If the network contains cycles that do not close perfectly, no set of ligand
energies can reproduce every edge exactly. Cinnabar therefore finds the best
global compromise. The final fitted ligand energies are internally consistent,
even when the original RBFE edges were not.

### Step 3: Align the calculated mean to experiment

RBFE calculations determine differences, but not the overall zero point. For
example, changing every calculated energy by +10 kJ/mol leaves every RBFE
difference unchanged.

The script calculates one common shift:

```text
shift = mean(experimental DG - calculated DG)
```

It adds this shift to every calculated ligand energy. This makes the mean of
the calculated values equal to the mean of the experimental values.

Experimental values therefore set only the common offset. They do not change
the calculated differences or ligand ranking.

### Meaning of the calculated uncertainty

Cinnabar propagates the uncertainties supplied for the RBFE edges through the
network fit. The output uncertainty describes how well the network determines
each fitted ligand energy under Cinnabar's statistical assumptions.

It does **not** include:

- uncertainty in the experimental measurements;
- uncertainty in the mean-alignment shift;
- force-field systematic error;
- errors caused by incorrect protonation states, binding poses or atom maps;
- hidden sampling problems not represented by the supplied edge uncertainties.

## Cycle-closure analysis

Run cycle closure with:

```bash
python rbfe_to_abfe.py network.csv experimental.csv \
    -o calculated_abfe.csv \
    --cycle-closure 5
```

Here, `5` is the largest number of ligands allowed in an assessed cycle.
Cinnabar examines all detected simple cycles containing 3, 4 or 5 ligands. It
does not examine only the largest cycle. Two-edge forward/backward comparisons
are not included.

"All" means every distinct simple loop that Cinnabar finds within the selected
size limit. A simple loop does not visit a ligand twice, apart from returning
to its starting ligand. The same loop is not counted again merely because it
is read from a different starting ligand or in the opposite direction.

Increasing the maximum can find longer cycles, but the number of possible
cycles can grow rapidly in a dense network. A value of 4 or 5 is a practical
starting point for many ligand networks.

### What should happen around a cycle?

For the cycle:

```text
A -> B -> C -> A
```

the signed sum should ideally be zero:

```text
DDG(A -> B) + DDG(B -> C) + DDG(C -> A) = 0
```

For example:

```text
A -> B     +2.0 kJ/mol
B -> C     -4.0 kJ/mol
C -> A     +2.0 kJ/mol
sum         0.0 kJ/mol
```

This cycle closes exactly. A non-zero sum means that at least one edge is not
fully consistent with the others. The cause may be insufficient sampling,
underestimated uncertainties, different binding modes or another setup issue.

## Cycle-level output

With `--cycle-closure N`, the script writes:

```text
calculated_abfe_cycle_closure.csv
```

The file contains one row for every detected cycle:

| Column | Meaning |
|---|---|
| `source` | Label assigned to the computational results |
| `cycle` | Ligands forming the cycle |
| `cc_kJ_mol` | Absolute cycle-closure error in kJ/mol |
| `cc_per_sqrt_edge_kJ_mol` | Closure error divided by the square root of the number of edges |
| `cc_unc_normalized` | Closure error divided by the combined edge uncertainty |

The raw closure error is:

```text
cc = absolute value of the signed sum of DDG values around the cycle
```

The propagated cycle uncertainty is approximately:

```text
sqrt(uncertainty_1^2 + uncertainty_2^2 + ...)
```

The normalized result is:

```text
cc_unc_normalized = cc / propagated cycle uncertainty
```

### Qualitative interpretation

The normalized closure error is usually the easiest value to interpret:

| `cc_unc_normalized` | Practical interpretation |
|---:|---|
| Below 1 | Closure error is smaller than the reported combined uncertainty |
| 1 to 2 | Mild tension; inspect together with replicate and sampling results |
| 2 to 3 | Suspicious inconsistency; inspect the edges in this cycle |
| Above 3 | Strong inconsistency relative to the reported uncertainties |

These ranges are guides, not universal acceptance rules. A high value can mean
that an edge is poorly converged, but it can also mean that the edge
uncertainties are too optimistic. A small value does not prove convergence,
because errors can cancel around a loop.

The raw `cc_kJ_mol` value is still useful because it shows the energetic size
of the disagreement. The `cc_per_sqrt_edge_kJ_mol` value makes cycles of
different lengths easier to compare.

## Edge-level cycle output

The script also writes:

```text
calculated_abfe_cycle_closure_edges.csv
```

This file summarizes the cycles involving each calculated edge:

| Column | Meaning |
|---|---|
| `source` | Label assigned to the computational results |
| `ligandA`, `ligandB` | Transformation being summarized |
| `n_cycles` | Number of detected cycles containing the edge |
| `mean_cc_per_sqrt_edge_kJ_mol` | Mean scaled closure error for those cycles |
| `max_cc_per_sqrt_edge_kJ_mol` | Largest scaled closure error involving the edge |

An edge that appears repeatedly in high-error cycles deserves attention. This
does not prove that the edge is wrong: several cycles can share the same
problem elsewhere. Cycle analysis identifies where to investigate, not which
single ligand or edge is guilty.

## Does cycle closure change the fitted energies?

The `--cycle-closure` flag only creates diagnostic reports. It does not apply a
separate correction and does not reject any ligands.

The maximum-likelihood fit already uses all edges, including redundant edges
in cycles. The cycle report simply shows how consistent the original RBFE
measurements were before they were combined into one set of ligand energies.

## Is there a convergence cutoff?

Not inside this reconstruction step. Cinnabar's MLE calculation solves a
weighted system of equations using linear algebra; it is not a molecular
simulation that becomes converged after a certain number of iterations. There
is therefore no Cinnabar convergence threshold that certifies the result as
trustworthy.

A successful calculation only means that Cinnabar could solve the supplied
network. Scientific confidence still comes from several checks together:

- agreement among independent simulation repeats;
- reasonable edge uncertainties;
- acceptable cycle-closure errors;
- stable results when simulations are extended or divided into time blocks;
- careful inspection of poses, protonation states, mappings and sampling.

Do not use one universal cycle-closure cutoff as an automatic pass/fail rule.
For example, a normalized closure error below 1 is reassuring relative to the
uncertainties you supplied, but it cannot reveal a systematic error shared by
every edge.

## Log files and keeping a record

Cinnabar 0.6.1 does not automatically create a detailed log file for this
workflow. The script prints a short summary to the terminal and writes the CSV
results. To keep the terminal messages, redirect them to a text file:

```bash
python rbfe_to_abfe.py network.csv experimental.csv \
    -o calculated_abfe.csv --cycle-closure 5 \
    > rbfe_to_abfe.log 2>&1
```

This log records the ligand count, alignment shift, number of cycles, largest
normalized closure error and output paths. It is a run summary, not a detailed
iteration history. For reproducibility, keep the two input CSV files, all
output CSV files, this log and the Cinnabar version together.

## What if no cycles are found?

If no cycle containing 3 to `N` ligands exists, the script prints a message and
creates empty cycle-report CSV files with column headers. The binding-energy
calculation can still succeed for a connected network without cycles, but
there is then no independent path available for a cycle-closure check.

## Errors checked by the script

The script stops with an explanatory message when it finds:

- missing required columns;
- an empty input file;
- non-numeric, missing or infinite energy values;
- zero or negative RBFE uncertainties;
- a transformation from a ligand to itself;
- a repeated unordered ligand pair when cycle closure is requested;
- duplicate ligand names in the experimental file;
- different ligand sets in the two input files;
- a disconnected RBFE network;
- a `--cycle-closure` value below 3.

The script cannot automatically detect a reversed sign, an incorrect ligand
name that happens to be used consistently, an unsuitable protonation state or
an unconverged simulation with an unrealistically small reported uncertainty.

If the script fails, read the last error message first. Most failures can be
fixed by correcting a column name, removing a missing value, aggregating repeat
rows, giving every edge a positive uncertainty, or adding transformations that
connect separate parts of the network. An import error usually means that the
wrong `cinnabar` package was installed; create the conda environment shown
above. A completed run with a large closure error is not a software failure—it
is a warning to inspect the transformations that contribute to that loop.

## Recommended workflow

1. Combine the independent repeats for every transformation.
2. Check the ligand names, edge directions and units.
3. Run the script with cycle closure enabled.
4. Inspect cycles with the largest `cc_unc_normalized` values.
5. Use the edge-level report to find transformations shared by several poor
   cycles.
6. Compare those transformations with replicate agreement, sampling and setup
   checks.
7. Treat the reconstructed ligand energies as one part of the overall
   validation, not as proof that every calculation converged.

## Further reading

- [Cinnabar MLE explanation](https://cinnabar.openfree.energy/en/latest/concepts/estimators.html)
- [Cinnabar cycle-closure tutorial](https://cinnabar.openfree.energy/en/latest/tutorials/cycle-closure.html)
- [Cinnabar FEMap API](https://cinnabar.openfree.energy/en/latest/generated/cinnabar.femap.FEMap.html)
