# Examples

These small files show the two command-line programs with deterministic dummy
data. They are demonstrations of the file formats and report layout, not
evidence that a molecular simulation has converged physically.

Run the examples from the repository root after installing the environment in
`environment.yml`:

```bash
python scripts/qualstat.py examples/qualstat_abfe/settings.yaml
python scripts/qualstat.py examples/qualstat_rbfe/settings.yaml
python scripts/rbfe_to_abfe.py \
  examples/rbfe_to_abfe/network.csv \
  examples/rbfe_to_abfe/experimental.csv \
  -o examples/rbfe_to_abfe/calculated_abfe.csv \
  --cycle-closure 4
```

Each QualStat settings file uses paths relative to its own directory. The
committed `expected.log` files are reference reports generated with the fixed
seed in the corresponding YAML file. Only the machine-specific `Input:` path
is replaced by a placeholder in those committed logs; the numerical report is
otherwise unchanged. Running the commands in the checkout writes the report to
`expected.log`, so copy the example directory first if you want to preserve the
reference file exactly.

Every uncertainty in these examples is the one-sigma uncertainty of the value
reported on the same row. If your row is the mean of `n` independent repeats
and you have their sample standard deviation `s`, supply `s / sqrt(n)`. For
triplicates, this is `s / sqrt(3)`.

## QualStat ABFE

`qualstat_abfe/` contains six independent ligand binding free energies with
non-zero calculated and experimental uncertainties. It demonstrates the ABFE
CSV schema, replacement bootstrap, and common error/correlation metrics.

## QualStat RBFE

`qualstat_rbfe/` contains four ligands and five directed transformations. The
edges form several simple cycles, so the report also shows directed cycle
closure and root-sum-square uncertainty. Its selected statistics are the
project's preferred orientation-independent set: MAD, RMSD, taurx, and r22.

## Input conventions

QualStat expects calculated and experimental values in the same energy unit.
The RBFE convention is the directed difference from `ligand_a` to
`ligand_b`; reversed duplicate pairs are not allowed. Zero uncertainty is
valid for QualStat and is useful for deterministic tests, although the examples
use non-zero uncertainties.

## RBFE-to-ABFE reconstruction

`rbfe_to_abfe/` uses a connected four-ligand network with cycles and a small
deliberate closure inconsistency. The input sign convention is:

```text
DDG(A -> B) = G(B) - G(A)
```

The three output CSV files contain network-reconstructed values, their network
uncertainties, and cycle diagnostics from Cinnabar or the numerically
equivalent NumPy fallback. The terminal output states which
estimator was used. The calculated energies are mean-aligned to the
experimental values; they are reconstructed relative-network energies and are
not independent absolute-binding-free-energy simulations. Cycle closure is a
diagnostic and does not change or reject the fitted values.
