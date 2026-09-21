# Python QualStat

This is a small Python replacement for the a QualStat Fortran program created by Ulf Ryde's group https://signe.teokem.lu.se/ulf/Methods/qual-stat.html. 
It preserves QualStat's legacy metric definitions and parametric bootstrap, accepts CSV/YAML input, adds Spearman's rho, and can calculate RBFE cycle-closure errors.

## Installation

Only PyYAML is required outside the Python standard library. Either install it into an existing Python 3.10+ environment or create the supplied Conda environment:

```bash
conda env create -f environment.yml
conda activate qualstat
```

## Commands

```bash
python scripts/qualstat.py --write-template
python scripts/qualstat.py --write-template my_settings.yaml
python scripts/qualstat.py my_settings.yaml
python scripts/qualstat.py -h
python scripts/qualstat.py -help
python scripts/qualstat.py --help
```

The first command creates `qualstat_template.yaml`. The template contains brief comments, every supported statistic as a `true`/`false` option, 1,000 bootstrap rounds, and a `.log` output name.

Input and output paths in YAML are resolved relative to the YAML file, not the current shell directory.

## CSV input

Headers and column order are exact. Uncertainties must be non-negative standard deviations.

ABFE:

```csv
ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty
Lig1,-25.1,0.8,-24.3,0.4
```

RBFE:

```csv
ligand_a,ligand_b,calculated,calculated_uncertainty,experimental,experimental_uncertainty
Lig1,Lig2,-3.0,0.5,-2.6,0.3
```

An RBFE value is interpreted as the directed difference `ligand_a -> ligand_b`. Self-edges and repeated ligand pairs—even if reversed—are rejected.

## Bootstrap

For every bootstrap round, the script independently draws each calculated and experimental value from a normal distribution centered on the input value with the input uncertainty as its standard deviation. It then recomputes the selected metrics. The report gives the original-data estimate, sample standard deviation of valid bootstrap values, and number of valid bootstrap rounds.

Use an integer `random_seed` for reproducible Python results or `null` for a random seed. Python and Fortran use different random-number generators, so their individual bootstrap samples are not expected to match.

## Metrics

Run `python scripts/qualstat.py --help` for a concise explanation of every metric. The following legacy details matter when comparing to other software:

- `r2` and `r22` retain the sign of Pearson correlation after squaring.
- `Median` and `AbsMed` use the lower middle value for an even sample.
- `tau` excludes tied experimental pairs, while a prediction tie counts as discordant.
- `taux` and `taurx` default to the legacy 1.645 significance multiplier.
- `PI` weights pair ordering by experimental separation; calculated ties add no numerator contribution.
- `ROCar` reproduces the legacy 21-threshold calculation rather than a modern exact AUC routine. It treats the most frequent experimental value as the inactive class; ties use first appearance, and if every label is unique, the first label is selected deterministically.
- `rho` is the added standard Spearman rank correlation with average ranks for ties.

Mathematically undefined values are printed as `nan` without stopping other metrics. For example, a correlation is undefined when either input has zero variance.

## RBFE cycle analysis

With `analysis_type: rbfe` and `cycle_analysis: true`, the script finds every unique simple cycle of length three or greater. Traversing an edge in its stored `ligand_a -> ligand_b` direction adds its calculated free energy; traversing it backward subtracts it. Independent edge uncertainties are combined by root-sum-square.

```text
Cycle #1

Lig1 -> Lig2 -> Lig3 -> Lig1 = 0.300 +/- 0.877 kJ/mol
```

Equivalent rotations and reversals are printed once in deterministic order. Enumerating all simple cycles can become expensive for a dense network because their number can grow exponentially.

## Examples and tests

```bash
python scripts/qualstat.py examples/qualstat_abfe/settings.yaml
python scripts/qualstat.py examples/qualstat_rbfe/settings.yaml
python -m unittest discover -s tests -v
```

The example reports are written into `examples/` with `.log` suffixes. The
tests validate software behavior and statistical definitions; they do not
establish physical convergence of an underlying molecular simulation.
