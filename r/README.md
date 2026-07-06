# r — the validation oracle

Validation-independent ground truth for the project's numerical estimators: the
same inputs run through external R implementations (Kisung You's packages,
base R, and Ripser via TDAstats), so the C# in `src/**` is checked against an
*external* implementation in a different language — not against itself. Sibling
of `tests/`, modelled on `lean/`: a self-contained sub-project whose **toolchain
lives outside** and whose only in-tree state is the pinned package set.

## Layout

| Path | Role |
| --- | --- |
| `oracles/*.R` | reference-value generators (read a fixture, emit JSON) |
| `oracles/_common.R` | shared fixture I/O helpers |
| `DESCRIPTION` | explicit oracle dependency manifest for `renv::snapshot()` |
| `scripts/bootstrap.R` | one-time provisioning: install packages + pin `renv.lock` |
| `scripts/r-session.ps1` | project-local R resolver + aliases (`rver`, `rnr`, `rns`, `rnsnap`, `rs`) |
| `scripts/oracle-ci.ps1` | health/repro gate — turn R on, check renv in sync, smoke an oracle |
| `.Renviron` | project R environment; disables renv's Windows sandbox for stable CLI/test activation |
| `renv.lock` | the pinned package set (reproducible via `renv::restore`) |

The R **toolchain is external** — normally the PDenv install under
`~/PDenv/rlang` (R 4.6.1), or any shell/user-level `R_HOME` / `PATH` that points
at `Rscript.exe`. R is **opt-in**: it is not assumed to be in the default env
(niche, validation-only). Only this project's *package library*
(`renv/library/`, gitignored, pinned by `renv.lock`) lives in-tree.

## One-time provisioning (network)

```powershell
Set-Location <repo>/r
. scripts/r-session.ps1     # resolves R and adds project-local aliases
Initialize-ROracleSession -SetAliases
Rscript scripts/bootstrap.R    # installs the oracle package set and writes renv.lock
```

Thereafter `rnr` (`renv::restore`) reproduces the exact library from `renv.lock`;
`rnsnap` (`renv::snapshot`) re-pins after adding a package.

## How the C# tests reach it

A parity test resolves `Rscript.exe` from the live/user environment or the PDenv
toolchain, runs it with **this directory as the working dir** (so `.Rprofile`
auto-activates the pinned library), and compares the emitted JSON within
tolerance. The fixture matrix is generated once and fed to *both* sides — R's
Mersenne-Twister and C#'s Xoshiro can't share a stream, so we compare *outputs*,
not RNG. Eigenvector / subspace comparisons are sign-agnostic + subspace-distance;
note `prcomp` eigenvalues use the (n−1) denominator vs the C# MLE (n) — rescale,
or compare the subspaces. These tests are **skipped when R is absent** (toolchain
is opt-in), so the normal `dotnet test` run is unaffected.

## Oracle map (grows with the parity work)

| C# | oracle | reference |
| --- | --- | --- |
| `Maths.LinAlg.Pca` | `pca_oracle.R` | base `prcomp` |
| `TDA.Ph.FullRips` + `PersistentHomology` | `tda_oracle.R` | Ripser via `TDAstats::calculate_homology` |
| `Maths.Estimators.MxPbf` | `mxpbf_oracle.R` *(todo)* | transcribe paper §2.2 / §3.1 |
| `Maths.Estimators.RobustDistributedPCA` | `mom_oracle.R` *(todo)* | `maotai` geometric median + the §5.3 rule |
| `TDA.DimReduction.Spred` | *(composition: component oracles + downstream ISOLET benchmark)* | PH via `tda_oracle.R`; median oracle pending |
