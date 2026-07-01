# r — the validation oracle

Validation-independent ground truth for the project's numerical estimators: the
same inputs run through Kisung You's reference R packages (Rdimtools, maotai,
T4cluster) and base R, so the C# in `src/maths/**` is checked against an
*external* implementation in a different language — not against itself. Sibling
of `tests/`, modelled on `lean/`: a self-contained sub-project whose **toolchain
lives outside** and whose only in-tree state is the pinned package set.

## Layout

| Path | Role |
| --- | --- |
| `oracles/*.R` | reference-value generators (read a fixture, emit JSON) |
| `oracles/_common.R` | shared fixture I/O helpers |
| `scripts/bootstrap.R` | one-time provisioning: install packages + pin `renv.lock` |
| `scripts/oracle-ci.ps1` | health/repro gate — turn R on, check renv in sync, smoke an oracle |
| `renv.lock` | the pinned package set (reproducible via `renv::restore`) |

The R **toolchain is external** — the portable install at `$PORTABLE_ROOT/rlang`
(R 4.6.0), turned on by sourcing `ps.core.bootstrap/helpers/env-Rlang.ps1`. R is
**opt-in**: it is not in the default env (niche, validation-only). Only this
project's *package library* (`renv/library/`, gitignored, pinned by `renv.lock`)
lives in-tree.

## One-time provisioning (network)

```powershell
. "$env:PORTABLE_ROOT/UserGithub/PowerShellCore/ps.core.bootstrap/helpers/env-Rlang.ps1"  # turn R on
Set-Location <repo>/r
Rscript scripts/bootstrap.R    # installs jsonlite + Rdimtools/maotai/T4cluster, writes renv.lock
```

Thereafter `rnr` (`renv::restore`) reproduces the exact library from `renv.lock`;
`rnsnap` (`renv::snapshot`) re-pins after adding a package.

## How the C# tests reach it

A parity test resolves `Rscript.exe` from the toolchain (`$PORTABLE_ROOT/rlang`),
runs it with **this directory as the working dir** (so `.Rprofile` auto-activates
the pinned library), and compares the emitted JSON within tolerance. The fixture
matrix is generated once and fed to *both* sides — R's Mersenne-Twister and C#'s
Xoshiro can't share a stream, so we compare *outputs*, not RNG. Eigenvector /
subspace comparisons are sign-agnostic + subspace-distance; note `prcomp`
eigenvalues use the (n−1) denominator vs the C# MLE (n) — rescale, or compare the
subspaces. These tests are **skipped when R is absent** (toolchain is opt-in), so
the normal `dotnet test` run is unaffected.

## Oracle map (grows with the parity work)

| C# | oracle | reference |
| --- | --- | --- |
| `Maths.LinAlg.Pca` | `pca_oracle.R` | base `prcomp` |
| `Maths.Estimators.MxPbf` | `mxpbf_oracle.R` *(todo)* | transcribe paper §2.2 / §3.1 |
| `Maths.Estimators.RobustDistributedPCA` | `mom_oracle.R` *(todo)* | `maotai` geometric median + the §5.3 rule |
| `Maths.LinAlg.Spred` | *(downstream: ISOLET benchmark, not R-parity)* | — |
