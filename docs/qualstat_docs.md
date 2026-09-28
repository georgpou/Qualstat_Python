# QualStat: usage and statistical definitions

`scripts/qualstat.py` compares calculated and experimental ABFE or RBFE values.
It preserves legacy QualStat definitions, adds ordinary Spearman `rho` and the
orientation-independent extension `rho2`, supports two bootstrap methods, and can
inspect RBFE cycle closure.

## Run the program

```bash
python scripts/qualstat.py --write-template
python scripts/qualstat.py qualstat_template.yaml
python scripts/qualstat.py --help
```

The template lists every metric. Input and output paths are relative to the
YAML file, not the shell's current directory.

Template generation refuses to replace an existing file. To replace one
intentionally, run `python scripts/qualstat.py --write-template --force`, or
place `--force` after an explicit template path.

## Input CSV files

Headers and their order are exact.

ABFE:

```csv
ligand,calculated,calculated_uncertainty,experimental,experimental_uncertainty
Lig1,-25.1,0.8,-24.3,0.4
```

Ligand names must be unique in ABFE mode.

RBFE:

```csv
ligand_a,ligand_b,calculated,calculated_uncertainty,experimental,experimental_uncertainty
Lig1,Lig2,-3.0,0.5,-2.6,0.3
```

Each RBFE row is directed from `ligand_a` to `ligand_b`:

```text
DDG(A -> B) = G(B) - G(A)
```

Self-edges and repeated unordered ligand pairs are rejected. Thus `A -> B` and
`B -> A` cannot both occur in one input file.

Calculated and experimental values must use the same unit. The program prints
the `energy_unit` label but does not convert values.

## What belongs in an uncertainty column?

Use the **one-sigma standard uncertainty of the reported energy**, not
automatically the raw spread across repeats.

If an energy is the mean of `n` independent repeats and `s` is their sample
standard deviation, the standard error of that mean is

```text
SEM = s / sqrt(n)
```

For three repeats, use `s / sqrt(3)`. If an upstream program already reports
the uncertainty or standard error of the mean, use it unchanged. This is the
same SD-versus-SEM distinction described by
[SciPy's `sem` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sem.html).

This conversion assumes independent repeats. Correlated repeats contain less
independent information, so `s / sqrt(n)` may then be too optimistic.

## Bootstrap sampling

Choose `bootstrap_method: replacement` (the default) or
`bootstrap_method: parametric` in the YAML file. Both use `bootstrap_rounds`
rounds and report the method in the output log.

- **Replacement:** Draw as many complete CSV rows as the input contains, with
  replacement. For example, three rows might yield `(2, 2, 0)`. Names, values,
  and uncertainties remain paired. Rebuild metric pair/significance selections
  for each sample.
- **Parametric (original Fortran method):** Keep each row exactly once. For each
  row, independently draw calculated and experimental energies from
  `Normal(reported energy, reported uncertainty)` using the CSV uncertainties
  as one-sigma standard deviations. Zero uncertainty leaves that energy fixed.
  Metric pair/significance selections are determined from the original values
  and held fixed over the draws, as in the Fortran `init()` routine. Zero
  experimental uncertainty freezes those draws, but significance selections
  still use the experimental values (including ties and zeros).

The `Estimate` column uses the original rows under either method.
`Bootstrap SD` is the sample standard deviation across finite round values;
`Valid rounds` counts rounds where that metric was defined. If fewer than two
rounds have finite values, the SD is `nan`.

With replacement, uncertainty columns do not perturb energies; they still
determine significant pairs for `taux` and `taurx` and contribute to RBFE cycle
uncertainties. With either method, cycle closure uses the original input values
and uncertainties. `random_seed` reproduces the same draws; set it to `null`
for a random seed.

## Recommended RBFE statistics

Accuracy metrics should be the main RBFE summaries. The default template uses
`MAD`, `RMSD`, `taurx`, and `r22`, matching this project's established
practice.

| Metric class | Metrics | Orientation behavior |
|---|---|---|
| Recommended accuracy | `MAD`, `RMSD` | Invariant |
| Preferred legacy RBFE | `taurx`, `r22` | Invariant |
| Optional project extension | `rho2` | Invariant |
| Other invariant summaries | `MQ`, `Q`, `AbsMed`, `taur`, `slope2`, `max` | Invariant |
| Do not use on raw RBFE edges | `MADtr`, `r2`, `PI`, `MSD`, `Median`, `slope`, `inter`, `multi`, `tau`, `taux`, `regMAD`, `ROCar`, `R`, `rho` | Depends on chosen edge directions |

Hahn et al. show that ordinary correlation values for an RBFE edge set can be
changed simply by reversing edge definitions and recommend accuracy measures
such as RMSE and MUE instead
([Living Journal of Computational Molecular Science, 2022](https://pubmed.ncbi.nlm.nih.gov/36382113/)).

### `rho2`: origin-symmetric Spearman correlation

There is no widely accepted standard orientation-independent Spearman
coefficient for RBFE edges. This repository therefore provides a clearly
labelled custom extension:

1. Start with every pair `(calculated_i, experimental_i)`.
2. Add `(-calculated_i, -experimental_i)`.
3. Calculate ordinary Spearman correlation on the doubled data.

Reversing one stored edge only swaps its two mirrored points, leaving the
doubled dataset—and therefore `rho2`—unchanged. Unlike Spearman correlation of
absolute magnitudes, `rho2` retains whether the calculated and experimental
signs agree. Because it is not a standard statistic, publications should state
this definition explicitly rather than calling it ordinary Spearman `rho`.

## Metric definitions

Let `c_i` be calculated values, `e_i` experimental values, and
`d_i = c_i - e_i`.

| Name | Definition and interpretation |
|---|---|
| `MAD` | Mean `abs(d_i)`; often called MAE or MUE elsewhere. |
| `RMSD` | Square root of mean `d_i^2`; often called RMSE. |
| `MSD` | Mean signed error, mean `d_i`. |
| `max` | Maximum `abs(d_i)`. |
| `R` | Standard Pearson correlation. |
| `rho` | Standard Spearman correlation using average ranks for ties. |
| `rho2` | Spearman correlation after adding sign-negated data; custom orientation-independent extension. |
| `r2` | Legacy signed `R * abs(R)`, not conventional non-negative `R^2`. |
| `r22` | Legacy `r2` after adding sign-negated data. |
| `taur` | Legacy sign-agreement score over non-zero experimental RBFEs. |
| `taurx` | `taur` restricted to edges significant relative to both supplied uncertainties. |
| `tau` | Legacy ordering score; experimental ties are omitted and calculated ties count as discordant. It is not Kendall tau-b. |
| `taux` | Legacy `tau` restricted to significant pairs. |
| `MADtr` | MAD after subtracting the mean signed error. |
| `Median` | Lower median of signed errors for an even-sized sample. |
| `AbsMed` | Lower median of absolute prediction errors; not median absolute deviation about a sample median. |
| `PI` | Legacy pair-ordering index weighted by experimental separation. |
| `MQ` | Mean calculated-to-experimental quotient. |
| `Q` | Sum of squared errors divided by sum of squared experimental values. |
| `slope`, `inter` | Regression of calculated values on experimental values. |
| `multi` | Legacy slope multiplied by the mean experimental value. |
| `slope2` | Regression slope after adding sign-negated data. |
| `regMAD` | MAD around the regression of experimental on calculated values. |
| `ROCar` | Compatibility-only 21-threshold legacy ROC quantity. Do not use it for continuous ABFE/RBFE energies. |

`tau`, `taux`, `taur`, and `taurx` return `nan` when no eligible comparisons
exist. Other mathematically undefined values, such as correlation with a
constant vector, also appear as `nan` without stopping the remaining metrics.

## RBFE cycle analysis

With `analysis_type: rbfe` and `cycle_analysis: true`, the program enumerates
unique simple cycles containing at least three ligands. Traversing a stored edge
forward adds its value; traversing it backward subtracts it. With independent
edge uncertainties, the cycle uncertainty is

```text
sqrt(sigma_1^2 + sigma_2^2 + ...)
```

Cycle enumeration can become expensive in dense networks because the number of
simple cycles may grow exponentially.
