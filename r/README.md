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
| `GeometricMedian<GrassmannManifold>` | `mom_oracle.R` | `Riemann::riem.median` on `wrap.grassmann` |
| `Maths.Estimators.MxPbf` | `mxpbf_oracle.R` *(todo)* | transcribe paper §2.2 / §3.1 |
| `Maths.Estimators.RobustDistributedPCA` | *(todo: extend `mom_oracle.R`)* | `maotai` geometric median + the §5.3 rule |
| `TDA.DimReduction.Spred` | *(composition: component oracles + downstream ISOLET benchmark)* | PH via `tda_oracle.R`; Grassmann median via `mom_oracle.R`; full SPRED composition still indirect |
| `Maths.Samplers.Ensemble.ChainDiagnostics` (rank-normalized split R̂, bulk and tail ESS, ESS of the mean) + `Maths.Regression.PredictiveChecks.PsisLoo` (k̂, r_eff, the elpd's Monte Carlo error) | `diagnostics_oracle.R` | `posterior` and `loo` (OQ9, N13); the tail diagnostics' split-chain convention is emitted beside posterior's |
| `Maths.Regression.Spline.Bass.BassBasis` (the case-side and output-side factors, R's scaling and normalization; stage 4's S2) | `bass_basis_oracle.R` | R BASS's own `makeBasis`, `makeBasisCat` and `scale_range` at the pinned commit; run by `BassBasisParityTests` |
| BASS's collapsed g-prior marginal over τ (`BassTarget`, stage 4's S5; checked through the posterior mode, whose closed forms the oracle checks against this one) | `bass_marginal_oracle.R` | base `integrate` on log τ; the port brief's D03 propriety conditions |
| BASS's order-free admission gate, under R's intercept convention and the Liang default (centred, and the direct-functional centred Hadamard form) *(port stage 4; not built)* | `bass_admission_oracle.R` | base `chol`; QR of the explicitly centred design for the centred modes; R BASS's `getQf` rule as a defect fixture |
| BASS's truncation residuals, in-sample and leave-one-case-out *(port stage 5; not built)* | `bass_truncation_oracle.R` | base `svd`; explicit refits are the reference, the case-coordinate computation is checked against them |
| BASS's τ prior, the unit-information Gamma family: Zellner–Siow against Francom 2018's shape 1/N (D03, F-O8) *(port stage 4; not built)* | `bass_tau_prior_oracle.R` | quadrature on ln τ against ψ and ψ′; Bayes factors against the intercept-only structure by quadrature |
| BASS's intercept conventions: Liang's flat intercept against level centring, with R's zero centre beside them (D03 and its proper-defaults alternative, F-O9) *(port stage 4; not built)* | `bass_intercept_oracle.R` | quadrature of the g-form Bayes factors; on an explicit design, R's τ form against the g-form identity, at three offsets of y |
| `Maths.Distributions.Wallenius`, BASS's variable-set proposal probabilities (D01, F-O5; stage 4's S3) | `bass_wallenius_oracle.R` | brute-force sums over draw orders, run by `WalleniusParityTests`; R BASS's `dmwnchBass` beside them as a defect fixture (its inclusion–exclusion loses digits under skewed weights) |
| BASS's tempering, swap and birth-ratio count laws (D01, D05, F-O6); the count, order and λ laws and the boundary ratios under D01 and D02 as decided (F-O3, F-O4) *(port stage 4; not built)* | `bass_enumeration_oracle.R` | enumeration of a small structure space and each kernel's transition matrix; R BASS's kernels, `lp` and fixed thirds as defect fixtures |
| BASS's exact posterior over an enumerable structure space with data (stage 4's S6): D03's collapsed marginal under the proper defaults (level centring, σ² scaled to the response, Zellner–Siow τ) and its options (fixed g, Liang's flat intercept, R's zero centre, the Jeffreys σ² with the exact-fit exclusion), with categorical predictors beside the continuous ones (D01's categorical part, both subset measures); `BassTarget`'s marginals, the structure prior and the exact-fit exclusion run by `BassTargetParityTests` (S5); the chain's count law, τ and inclusions by `BassChainParityTests` (S6) | `bass_enumeration_oracle.R` (`posterior`) | enumeration of F with τ by Gauss–Legendre on log τ; each closed form against quadrature over σ² of the explicit Gaussian integral, and against `bass_marginal_oracle.R` and `bass_intercept_oracle.R` |
| BASS's tempered laws and ladder evidence on that space, the integrated and linear forms (D05; stage 4's S8): log Z_b and the slope's mean from `BassTarget` run by `BassTargetParityTests` (S5), the ladder in S8 | `bass_enumeration_oracle.R` (`posterior`, `rungs`) | the exact law and log Z_b per rung; the slope's mean against d log Z_b/db by central differences |
| `Maths.Regression.Spline.Bass.BassDesign`'s uncentred Gram and the collapsed target's statistics, scalar and direct functional (X'X, X'y, y'y, qf, the factor's diagonal and log-determinant, σ²'s rate at fixed τ; stage 4's S4 and S5) | `bass_statistics_oracle.R` | R BASS's own statements in `bass`, `birth_des`, `birth_des_func`, `addBasis*`, `getQf` and `updateMCMC`, evaluated as written; the explicit design and a direct solve; `getQf`'s refusal as a flag (D06's defect fixture); the zero centre's S against `s2_rate_direct` run by `BassTargetParityTests` (S5) |
| R BASS's Jeffreys mode (D03, F-O6) *(defect fixture)* | `bass_jeffreys_oracle.R` | the coherent flat and Jeffreys coefficient priors' factors |
| BASS's output path and sensitivity (D07, D08, O-O1, O-O2) *(port stages 5–6; not built)* | `bass_output_defects_oracle.R`, `bass_sensitivity_defects_oracle.R` | corrected expectations beside R BASS's behaviour |
| BASS's hinge-product moments (D08, O-O3) *(port stage 6; not built)* | `bass_moments_oracle.R` | Gauss–Legendre quadrature and category enumeration |
| BASS calibration's defects and covariance formulations (D09, C-O1, C-O2) *(port stage 7; not built)* | `bass_calibration_defects_oracle.R`, `bass_covariance_oracle.R` | dense eigen and Cholesky references; R BASS's Woodbury helpers for parity |

`scripts/bass_stage1_checks.R` runs the BASS oracles on the port brief's fixtures and reproduces their numbers (about two minutes).
`scripts/bass_stage4_checks.R` runs the stage-4 oracle modes (the posterior, its rungs, R's statistics) on the stage-4 oracle brief's fixtures and prints PASS or FAIL for each self-check (about five minutes).

## Current caveats

On Windows, the `TDAstats`/Ripser path and Riemann's intrinsic Grassmann median
can write complete JSON and then exit `Rscript` with `0xC0000005` during native
teardown. The C# bridge tolerates that only for the known oracle scripts and
only after the output file exists and parses.
