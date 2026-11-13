# B1 stage two: ridge, random intercepts and MoMPCA

## Completed — 2026-11-12

Scope: Decision 3.2 in [the consolidation brief](spectral-solver-consolidation-brief.md). Migrate RidgeSpectrum, BassRandomIntercepts and MoMPCA to the existing element-wise relative rule. Graph callers, the absolute default and final cleanup remain follow-up work. No subagents are dispatched.

Baseline: `cb94b17a064501317b62c8c3367eb2b05ca766e2`, archived under `artifacts/spectral-b1-stage-two/baseline-source`; unrelated guidance and notes are excluded. Builds, probes and receipts stay under `artifacts/spectral-b1-stage-two`. Preserve RidgeSpectrum's statistical column normalization and each caller's validation and cutoff policies. Keep code commits targeted on main.

The three calls are migrated. Focused regression passed 39 in each build, against 37 on the archived baseline. Geometry passed 91 in each build, against 90 on the baseline. LinAlg passed 167 in each build, including all three decided kernel witnesses and stalled convergence. The full default regression gate finished with 3,281 passes and only three old pin mismatches, with no skips. The new repeated-block MoMPCA witness also passes with the pre-stage-two call; it is analytic coverage, not evidence of a former-rule failure. Code commits: `f48613f` (MoMPCA) and `672d1c2` (ridge/intercept spectra and three qualified pin updates). The rebuilt commit gate passed 2,778 tests, with zero failures or skips.

The full default regression run finished: 3,281 passed of 3,284, exactly the three documented old held-scale pin mismatches, zero skips (`candidate-source/default/test-select/20261112_162319-slice/Slice.trx`). No numerical or distributional contract test failed. The snapshot was rebuilt with the qualified pin expectations, and its commit gate passed all 2,778 cases (`candidate-source/default/test-select/20261112_183617-commit/Commit.trx`). Production sources are unchanged; `full-run-source-hashes.json` records the pre-pin source state. The full slow tier was not repeated for changed expectations alone. Graph callers and final cleanup remain follow-up slices.

## Changes

- `RidgeSpectrum.Of` selects `ElementWiseRelative` at 1e-15. Its unit-diagonal correlation matrix, column roots, mass normalization, multiplication of eigenvalues by mass, projections, positivity/admission checks and shrinkage formulas stay. This normalization defines the statistical prior and is not removed as a solver workaround. The class docstring now names the relative test.
- `RandomInterceptBlock.Of` in `BassRandomIntercepts.cs` selects the same rule at 1e-15 for HᵀDH. Basis construction, observed-level handling, positivity checks, spectral whitening, determinant and group-conditional formulas stay.
- `MoMPCA` selects the rule at 1e-12 for the median projected covariance. The caller no longer computes an absolute tolerance from the largest matrix entry. Local PCA still reads factors through ThinSvd; aggregation, principal-axis reconstruction, orientation and nonnegative variance clamp stay.

The dense kernel is unchanged. The absolute default and graph callers remain available for stage three. `ScaledEigen` and the legacy pencil helper remain for the separately commissioned final cleanup. No superseded public surface is deleted in this slice.

## Qualification fixtures and limits

### Ridge normalization and shrinkage

`RelativeJacobiStageTwoTests.Ridge_NormalizationAndWeakDirectionSurviveColumnRescaling` reads an analytic 3×3 Gram with a correlation pair rho = 1 − 2^-30 and an independent third direction. It uses mass 16 with rows 19, so the weighted normalization is distinguished from the row count. Column scales are [2^k, 1, 2^-k] at k = -300, 0, 300; all Gram entries remain representable.

The known generalized spectrum is [16(1+rho), 16, 16(1−rho)]. The smallest eigenvalue is 1.4901161193847656e-8. Its observed relative error, like the other two eigenvalues', is zero on these dyadic fixtures in both builds. This is not a universal error bound for nearly singular Grams.

Cross products put energy into the weak direction. The analytic determinant, shrinkage and their slopes are checked over tau = {1e-10,1,1e10}, b = {0,.3,1}, with a 2e-12 scaled bound. `BassRidgeTests` separately checks dense marginals under every intercept convention, column rescaling, the orthogonal g-form identity, admission and bordered candidates against fresh decompositions.

### Random intercepts

`UnbalancedIntercepts_WhitenTheIndependentDenseCovariance` uses counts [1,2,5,20], multiplied by 1 and 3 (28 and 84 observations), and s = {1e-3,1,1e3}. It forms Pi = I − nnᵀ/(nᵀn) directly, then V = I + s²ZPiZᵀ, independently of the production contrast basis. Applying the returned whitening to coordinate vectors gives W; the test reads W V W − I and compares the determinant against DenseCholesky. It also checks that W preserves the constant.

Both builds give the same printed results. The largest entrywise whitening residual is 2.042e-10, and the largest absolute log-determinant gap is 1.170e-9, both on the 84-observation, s=1000 case. Bounds are 2e-12 at s≤1 and 1e-7 at s=1000, where the dense covariance is more ill-conditioned. Full per-case values are in `artifacts/spectral-b1-stage-two/numerical-summary.json`. Existing conditional-covariance, joint-whitening and dense-target checks also run in the focused gate.

### MoMPCA

The existing scale-equivariance test covers two nontrivial aggregation fixtures, with scale factors 1e-12, 1e-9, 1e-6, 1e-3, 1e2, 1e3 and 1e6. It checks mean and variance scaling, basis agreement and the calibrated alpha at its existing tolerances. All existing geometry tests pass in both builds.

The new repeated-block witness uses six symmetric points per block, five identical blocks, and the full three-dimensional subspace. Its covariance is a unit direction beside 1e-20·[[1,.5],[.5,1]], at data scale factors {1e-60,1,1e60}. Expected variances are scale²·[1,1.5e-20,.5e-20]; each is checked to 2e-11 relative error, with orthogonality at 2e-12. This is analytic weak-variance coverage through the public pipeline. It is not, by itself, proof that the old MoMPCA stopping rule fails: local PCA can align the weak axes before the final dense solve. The explicit old-call comparison passed in the default build (one test, `old-mom-witness/old-mom-witness.trx`). That snapshot contains the unmodified pre-stage-two production code plus the new test. The repeated-block witness therefore does not establish changed behavior or an old-rule regression. No such failure is claimed; the local PCA alignment noted in review is consistent with this result. The nontrivial existing aggregation/rescaling fixtures supply separate equivariance coverage.

## Commands and receipts

TRX paths below are relative to `artifacts/spectral-b1-stage-two/`.

| Gate | Outcome | Receipt |
|---|---|---|
| Baseline regression focus | 37 passed | `baseline-source-default-focused/focused.trx` |
| Baseline geometry | 90 passed | `baseline-source-default-geometry/geometry.trx` |
| Candidate regression focus, default / reference | 39 / 39 passed | `candidate-source-{default,reference}-focused/focused.trx` |
| Candidate full geometry, default / reference | 91 / 91 passed | `candidate-source-{default,reference}-geometry/geometry.trx` |
| Candidate LinAlg, default / reference | 167 / 167 passed | `candidate-source-{default,reference}-linalg/linalg.trx` |
| Old MoMPCA code with new witness, default | One passed; not a former-rule failure witness | `old-mom-witness/old-mom-witness.trx` |
| Baseline held-scale pins | Three passed | `pin-baseline/pin-baseline.trx` |
| Full default regression before re-pinning | 3,281 passed, three old pin mismatches, zero skips | `candidate-source/default/test-select/20261112_162319-slice/Slice.trx` |
| Rebuilt commit gate after re-pinning | 2,778 passed, zero failed or skipped | `candidate-source/default/test-select/20261112_183617-commit/Commit.trx` |
| Two-seed-per-axis comparative sweep, baseline / candidate | Passed each; all 166 compared values and outcomes unchanged | `seed-sweep-baseline/seed-sweep.trx`, `seed-sweep/seed-sweep.trx` |

All source snapshots are archived from `cb94b17a064501317b62c8c3367eb2b05ca766e2`, with only the paths in `artifacts/spectral-b1-stage-two/owned-paths.txt` overlaid for the candidate. `candidate-hashes.json` identifies those files. The old-MoMPCA comparison uses a separate baseline archive plus only the new geometry test file (`old-mom-source`), leaving the baseline and candidate unchanged.

The only source difference between the full regression run and rebuilt commit gate is `RandomEffectsTests.cs`: its three expectations and explanatory comment. `full-run-source-hashes.json` preserves the former values' source; the candidate manifest records all six final owned files. The full run's failures are retained in `full-run-pin-failures.json`, and the 11 passed exact-law/evidence cases in `full-law-checks.json`.

`artifacts/spectral-b1-stage-two/run-gate.ps1` records its exact argument arrays, settings, logs and TRX. Release builds use isolated roots, `--disable-build-servers`, MSBuild -m:1, below-normal priority, and six xUnit workers. DOTNET_PROCESSOR_COUNT is also capped at six. SourceRevisionId is explicitly set to the archive's base revision, avoiding the live-parent-HEAD metadata caveat found in stage one; candidate overlays are identified by their source hashes.

```powershell
& artifacts/spectral-b1-stage-two/run-gate.ps1 -Source baseline-source -Suite focused
& artifacts/spectral-b1-stage-two/run-gate.ps1 -Source baseline-source -Suite geometry
& artifacts/spectral-b1-stage-two/run-gate.ps1 -Suite focused
& artifacts/spectral-b1-stage-two/run-gate.ps1 -Build reference -Suite focused
& artifacts/spectral-b1-stage-two/run-gate.ps1 -Suite geometry
& artifacts/spectral-b1-stage-two/run-gate.ps1 -Build reference -Suite geometry
& artifacts/spectral-b1-stage-two/run-gate.ps1 -Suite linalg
& artifacts/spectral-b1-stage-two/run-gate.ps1 -Build reference -Suite linalg
& artifacts/spectral-b1-stage-two/run-seed-sweep.ps1 -Source baseline-source
& artifacts/spectral-b1-stage-two/run-seed-sweep.ps1 -Source candidate-source
& artifacts/spectral-b1-stage-two/pin-probe/run.ps1
$owned = Get-Content -LiteralPath artifacts/spectral-b1-stage-two/owned-paths.txt
& scripts/test-select.ps1 -SourceRoot artifacts/spectral-b1-stage-two/candidate-source -Scope Slice -Changed $owned -NoBuild -ArtifactsRoot D:/aghado01/ThermoMapper/artifacts/spectral-b1-stage-two/candidate-source/default -Threads 6
& scripts/test-select.ps1 -SourceRoot artifacts/spectral-b1-stage-two/candidate-source -Scope Commit -Changed $owned -NoBuild -ArtifactsRoot D:/aghado01/ThermoMapper/artifacts/spectral-b1-stage-two/candidate-source/default -Threads 6
```

Focused runs set BASS_QUALIFICATION=1; the full regression run leaves it unset and includes the slow and pinned tiers. Reference builds add ExtraDefineConstants=EIGEN_REFERENCE. Existing restore audit-data access warnings are retained in the logs; no numerical failure is inferred from a restore warning. No performance result is claimed on the shared machine.

The registered comparative sweep ran on baseline and candidate with SEED_SWEEP=2 (two derived seeds per axis), DOTNET_PROCESSOR_COUNT=6 and below-normal priority. `run-seed-sweep.ps1` records replayable commands. All 166 comparison values and outcomes agree numerically (`seed-sweep/comparison.json`). The same seven summary rows retain failed inequalities: five exploratory flat-histogram readings, and the two-ridges recovery claim's mean-C and mean-subspace angles on one of two chain seeds. Those are the pre-existing limits recorded in stage one. This bounded sweep is not pathwise verification of the three held-group pinned fits, nor a general stochastic-robustness claim.

## Pins and review

The unchanged baseline reproduces all three `RandomEffectsTests.AHeldScaleFit_KeepsItsDraws` expectations (`pin-baseline/pin-baseline.trx`). A reflection probe calls that test's exact private fixture and hash on each saved build and exports all 400 retained draws per fixture. It verifies these current values before the expectations are updated:

| Fixture | Baseline pin | Candidate pin | Retained structures that differ |
|---|---|---|---:|
| hinge | b9b95aee57ef459c | 36bec2bb56e6cb77 | 129 / 400 |
| legendre | f07a33dd35774c18 | dde825de6c642a96 | 0 / 400 |
| tempered hinge | 05122582e8f56f75 | 6a2e590c8c64b44e | 184 / 400 |

The hinge fits do **not** remain close draw for draw. Maximum retained-draw separations are 0.05879 and 0.10111 in tau, 0.09721 and 0.16989 in sigma², and 0.6317 and 0.8485 in group effects, for hinge and tempered hinge respectively. These are diverged sampler trajectories after the changed single-group Gram arithmetic, not numerical approximation-error bounds. The single observation-centred grouping routes through `RandomInterceptBlock` in `RandomEffectsBlock.Of`; these fixtures exercise that route. Acceptance is based on the dense covariance/conditional references, the exact-law regression tests and the paired seed sweep, not a claim of pathwise bit identity for a changed solver.

The Legendre fixture retains its structures and tau. Its maximum sigma² gap is 8.89e-16, coefficient gap 2.394e-13, and group-effect gap 2.147e-13. Raw draws, hashes and the replayable `Compare.ps1` are under `pin-probe`; `comparison.json` states explicitly that coefficient comparisons use matching term sequences only. The updated test retains its Pinned trait and records the stage-two reason. All 11 full-run cases covering the fixed-scale and drawn-scale ladder laws and SMC/evidence checks passed (`full-law-checks.json`); all other numerical/property tests passed too. The rebuilt commit gate passed all 2,778 cases after these updates. No other pin is changed.

Final self-review checked the three rule selections, preserved statistical normalization and caller cutoffs, the old-MoMPCA witness result, the three pin qualifications, matching hashes for all six owned code/test files and exclusion of unrelated changes. No outstanding code fix remains in this slice. The bounded numerical and stochastic claims above are retained; no performance gain is claimed.

## Doc changes owed

- `docs/Maths/LinearAlgebra.md`: stage-two callers now select the relative rule; absolute graph tolerance units stay until stage three.
- `docs/Maths/Regression/BASS.md`: RidgeSpectrum retains its statistical normalization but uses a relative eigen solve; the single-group random-intercept Gram does likewise.
- `docs/Maths/RiemannianGeometry.md`: the MoMPCA description should retire the locally computed largest-entry tolerance and name the direct relative solve.
- No `docs/` edits in this worker slice.
