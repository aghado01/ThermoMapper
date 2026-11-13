# B1 stage one: relative dense Jacobi and first callers

## Completed — 2026-11-12

Owner: the stage-one Sol worker, commissioned by Review WIP items. Scope is Decisions 3.1 and 4 in [the consolidation brief](spectral-solver-consolidation-brief.md): named element-wise stopping alongside the absolute rule, truthful convergence, retained reference dispatch, ScaledEigen and PriorPencil only. Later migrations and deletions remain follow-up work.

Baseline: `e3a3b7f597f98d7a71f9eaf3489fe903ad4d22b2`, archived under `artifacts/spectral-b1/baseline-source`. Existing edits to guidance and unrelated notes are excluded. Builds, qualification probes and receipts live under `artifacts/spectral-b1`, with isolated output roots and capped below-normal test runs. No timings are planned on the shared machine.

Kernel landed as `ea76fdf` (`feat(linalg): add element-wise relative Jacobi stopping`). Candidate snapshot is the baseline archive plus the paths in `artifacts/spectral-b1/owned-paths.txt`; unrelated GP-worker changes are excluded. Candidate focus passed 66 in each build, including the 200-case E8 sweep. LinAlg passed 167, focused dense/graph passed 46 and geometry passed 90 in each build. The full default regression run finished with 3,274 passes and five explained former-pin mismatches, no skips. All 166 before/after seed-sweep comparison outcomes agree. Five pins are deliberately updated, qualified below.

Concrete limit: [[1e300, 1e-300], [1e-300, 0]] requires a rotation whose tangent underflows. The pivot is preserved: one sweep, zero rotations, `Converged = false`, `Stalled = true`, on allocating and caller-owned paths. This is reported as a representability limit, not silently cleared and not treated as a reason to add a global floor.

Caller migration landed as `13ef477` (`refactor(regression): migrate E8 and prior pencils to relative Jacobi`). The rebuilt commit gate passed 2,901 tests with no skips. NewtonMetropolis passed three tests in each build. The full slow tier was not repeated for changed expected hashes alone. Later B1 migration stages and the separate final cleanup remain follow-up slices; no timings are claimed.

## Implementation and qualification

### Kernel contract

`JacobiStoppingRule` names `Absolute` and `ElementWiseRelative`. Every DenseEigen entry, EigenFast entry and the retained scalar Eigen entry accepts it. The default remains absolute for the later stages. Existing absolute arithmetic and its boundary convention stay: fast pivots rotate at equality, reference pivots do not.

The relative rule is |a_pq| ≤ tol √|a_pp a_qq|, with no matrix-norm floor and no definiteness split. Square roots are taken separately; tol multiplies the larger root before the two roots are multiplied. This avoids overflow/underflow in the diagonal product and avoids prematurely underflowing tol times the smaller root on an uneven pair. A zero pivot passes, including when a diagonal is zero. A nonzero zero-diagonal pivot fails and must rotate. Relative tol is finite and in [0, 1).

Relative rotations normalize just the pivot pair before computing a stable algebraic tangent. This avoids overflowing `aqq-app`, `2*apq` or theta² at large scales and preserves the `app-t*apq`, `aqq+t*apq` diagonal update on a graded pair. An unrepresentable tangent is not applied or counted, and its pivot is not cleared. Convergence scans the final working matrix with the same pivot test, after the last sweep, on allocating, caller-owned and eigenvalues-only paths. `Stalled` means a whole sweep rotated nothing although a final pivot failed. A cap with remaining representable work is nonconverged without that flag.

The reference build still executes the scalar reference kernel, including its mirrored matrix updates and allocating buffers. Its relative branch uses the same accuracy predicate and stable tangent as the fast branch, rather than the absolute branch's trigonometric rotation. Thus it is a build comparison with the requested relative semantics, not an independent algorithmic oracle for those shared primitives. Analytical spectra, residuals and orthogonality are the independent gates. `DecompositionLog` assertions verify the selected kernel in both build variants. DLL sizes/hashes are saved in `artifacts/spectral-b1/assembly-hashes.json`.

### Stage-one callers

All production callers of `ScaledEigen.Of` and `.Into` now call DenseEigen directly with the relative rule and 1e-14. This includes E8, whitening, matrix Poincaré and BASS sensitivity. `ScaledEigen.cs` remains unchanged for the explicit old/new qualification and final cleanup.

PriorPencil reads the whitened likelihood directly with the relative rule and 1e-15. Its prior whitening, positivity checks, statistical prior, cross products and spectral shrinkage formulas stay. Its former rescaling operation is retained as internal `DecomposeAtUnitScale` for old/new qualification; production does not call it. This helper and ScaledEigen are cleanup debt for the final B1 slice. RidgeSpectrum, BassRandomIntercepts, MoMPCA and graph tolerance units have not migrated.

### Numerical range actually checked

`RelativeJacobiTests` exercises fast, reference and DenseEigen dispatch independently, then checks allocating flat input, caller-owned columns and eigenvalues alone for identical values/reports within each dispatch:

| Fixture | Scale factors | Gate |
|---|---|---|
| Unit entry plus 1e-20·[[1, .5], [.5, 1]] | 10^k, k = -240, -120, -20, 0, 20, 120, 240 | Every eigenvalue, including .5e-20 and 1.5e-20 before scaling, within 2e-13 relative |
| [[0, 1], [1, 0]] | Same seven factors | ±scale within 2e-13 relative; nonzero pivot rotates |
| [[1, 5e-21], [5e-21, 1e-40]] | Same seven factors | Weak eigenvalue 7.5e-41·scale within 2e-13 relative |
| Equal-diagonal positive block with separate eigenvalue; 4×4 Hadamard indefinite control | 10^k, k = -280, -200, -100, 0, 100, 200, 280 | Known spectrum, residual and orthogonality |
| Dense 8×8 Hadamard spectrum [8,7,6,5,4,3,2,-1] | Unit scale | SIMD rotation paths; independently recomputed off(VᵀAV) at one sweep agrees with both reported norms within 1e-13 |

For the scale sweep, each residual component is bounded by 2e-13 times the largest absolute expected eigenvalue in the original units, computed after division by that scale to avoid overflow in the check. Orthogonality is bounded by 2e-13 per entry. Weak eigenvalue checks are separate from residuals, so a tiny independent block cannot hide behind a good global residual. The smallest decided-witness weak eigenvalue tested is 7.5e-281; this is not a claim of relative accuracy over all subnormal inputs or every matrix with that dynamic range.

Zero sweeps, one sweep, zero relative tolerance on a zero matrix, equality at the relative threshold, and a stalled unrepresentable rotation are also checked. The explicit stall case above is a counterexample to unconditional convergence, not to the rule on the three decided cases. No tested case required a definiteness split. Whole-range arbitrary SPD/indefinite accuracy is not established.

`RelativeEigenMigrationTests` compares well-conditioned scaled spectra against ScaledEigen over 10^-240 through 10^240. A nontrivial diagonal prior is whitened independently, then the old pencil solve and current pencil agree in determinant, shrinkage and conditional log evidence to 2e-13 over tau = {1e-8,.1,1,100,1e8}, b = {0,.4,1}. A separate graded pencil puts cross-product energy into the weak direction: eigenvalue .5e-20, shrinkage tau/(tau+.5e-20), and the analytic log determinant pass at tau = {1e-22,1e-20,1}. This confirms the weak direction affects the statistical terms, not just the returned spectrum. MatrixTReferenceTests also retains its independent dense/analytic reference checks.

### E8 equivariance

The existing 200-case seeded `Whitening_IsAffineEquivariant_SeededSweep` uses p = 2–5, n = 50–400, weighted and unweighted samples and affine scales from 1e-3 through 1e3. All three stated bounds are 1e-9.

| Build/state | Eigenvalue gap | QᵀQ − I | Congruence gap |
|---|---:|---:|---:|
| Archived baseline, default | 4.25e-13 | 1.76e-11 | 6.15e-14 |
| Candidate, default | 2.46e-13 | 7.75e-12 | 6.31e-14 |
| Candidate, reference | 2.42e-13 | 7.75e-12 | 6.31e-14 |

These are this fixed seeded sweep, not a universal bound. The printed summaries are in the focused TRX files and `artifacts/spectral-b1/e8-summary.json`.

## Commands and receipts

| Gate | Result | TRX receipt under `artifacts/spectral-b1/` |
|---|---|---|
| Archived baseline dense/graph selection | 236 passed, two standing Möbius skips | `baseline-default/results/baseline.trx` |
| Archived baseline regression focus | 64 passed | `baseline-source-default-focused/focused.trx` |
| Candidate LinAlg, default / reference | 167 / 167 passed | `candidate-source-{default,reference}-linalg/linalg.trx` |
| Candidate E8/pencil/readout focus, default / reference | 66 / 66 passed | `candidate-source-{default,reference}-focused/focused.trx` |
| Candidate dense/graph dispatch focus, default / reference | 46 / 46 passed | `candidate-source-{default,reference}-viz/viz.trx` |
| Candidate geometry, default / reference | 90 / 90 passed | `geometry-{default,reference}/geometry.trx` |
| NewtonMetropolis absolute caller, default / reference | 3 / 3 passed | `newton-{default,reference}/newton.trx` |
| Full default regression, including slow and pinned tiers, before re-pinning | 3,274 passed, five old pin mismatches, zero skips; 3,279 total | `candidate-source/default/test-select/20261112_151843-slice/Slice.trx` |
| Rebuilt default commit gate after deliberate re-pinning | 2,901 passed, zero failed or skipped | `candidate-source/default/test-select/20261112_154716-commit/Commit.trx` |
| Registered seed sweep, two seeds per axis, baseline / candidate | Test passed each; 166 comparisons with unchanged outcomes | `seed-sweep-baseline/seed-sweep.trx`, `seed-sweep/seed-sweep.trx` |

The full run's five failures are exactly the two modified-g and three noise pins documented below, saved separately in `full-run-pin-failures.json`. No numerical/property test failed. After changing those expectations and adding the missing pinned trait, the rebuilt commit gate passed all 2,901 selected tests, with no skips. The whole slow tier is not repeated for a change of expected hashes alone.

Source baseline is the `git archive` of `e3a3b7f597f98d7a71f9eaf3489fe903ad4d22b2`. Candidate source is that archive plus the owned files; their SHA256 hashes are in `artifacts/spectral-b1/candidate-hashes.json`. This excludes concurrent guidance changes and the GP worker's patch. .NET SDK is 10.0.201; complete runtime information is in `artifacts/spectral-b1/dotnet-info.txt`. All builds are Release, with isolated build roots, build servers disabled, MSBuild -m:1 on direct runs, below-normal priority, and xUnit capped at two threads (six for the final commit gate, the project default cap). No runtime performance claim is made under this load.

For the pre-pin full run, restore the three pin-test files from the baseline while retaining the other 14 owned paths; `full-run-source-hashes.json` identifies those exact sources. The final candidate manifest includes all 17 owned code/test files, including the updated expectations. The code/source state of the final commit gate therefore differs from the full run only in those three test files.

The replayable focused runner is `artifacts/spectral-b1/run-gate.ps1`. Its exact dotnet argument arrays and runsettings accompany every receipt. For example:

```powershell
& artifacts/spectral-b1/run-gate.ps1 -Source baseline-source -Suite focused
& artifacts/spectral-b1/run-gate.ps1 -Suite linalg
& artifacts/spectral-b1/run-gate.ps1 -Build reference -Suite linalg
& artifacts/spectral-b1/run-gate.ps1 -Suite focused
& artifacts/spectral-b1/run-gate.ps1 -Build reference -Suite focused
& artifacts/spectral-b1/run-gate.ps1 -Suite viz
& artifacts/spectral-b1/run-gate.ps1 -Build reference -Suite viz
$owned = Get-Content -LiteralPath artifacts/spectral-b1/owned-paths.txt
& scripts/test-select.ps1 -SourceRoot artifacts/spectral-b1/candidate-source -Scope Slice -Changed $owned -NoBuild -ArtifactsRoot D:/aghado01/ThermoMapper/artifacts/spectral-b1/candidate-source/default -Threads 2
& scripts/test-select.ps1 -SourceRoot artifacts/spectral-b1/candidate-source -Scope Commit -Changed $owned -NoBuild -ArtifactsRoot D:/aghado01/ThermoMapper/artifacts/spectral-b1/candidate-source/default -Threads 6
```

The reference runner adds `-p:ExtraDefineConstants=EIGEN_REFERENCE`; focused regression runs set `BASS_QUALIFICATION=1`. The whole regression assembly runs once without that flag, including the pinned and slow tiers. A direct native build precedes its `test-select -NoBuild` invocation.

Archived-source provenance caveat: .NET discovers the parent checkout's `.git` and embeds its live HEAD in `AssemblyInformationalVersion`, rather than the archive's source revision. This was verified in the generated AssemblyInfo: the focused binaries embed `e6a1e0f`, and the rebuilt default binaries embed `ea76fdf`; the actual candidate source stayed fixed apart from the three pin-test files. Thus the archive plus owned-file manifest defines the source pin, not that generated metadata. Original and post-rebuild DLL hashes are saved separately (`assembly-hashes.json`, `assembly-hashes-after-pin-rebuild.json`). A follow-up improvement to snapshot test tooling could set `SourceRevisionId` explicitly from its source manifest.

Initial `test-select` build launches exited 1 without compiler diagnostics or TRX. The same selected baseline dense/graph gate succeeded through direct dotnet with the saved runsettings and -m:1 (`artifacts/spectral-b1/baseline-direct.log`, `baseline-default/results/baseline.trx`). NoBuild test-select runs then worked. This remains an unisolated build-launch problem; it is not counted as a numerical failure or a successful baseline run. Restore emitted NU1900 access-denied audit-data warnings in both baseline and candidate; builds and tests still completed.

## Pins and review

The registered comparative sweep ran with `SEED_SWEEP=2` on the candidate and archived baseline, capped with `DOTNET_PROCESSOR_COUNT=2`, at below-normal priority. Receipts and TSVs are in `artifacts/spectral-b1/seed-sweep` and `seed-sweep-baseline`. All 166 comparison outcomes agree between builds; the largest absolute change of a compared value is 3.4313e-15 (`seed-sweep/comparison.json`). This is a two-seed-per-axis smoke qualification, not a broad stochastic robustness claim.

The same seven summary rows have failing inequalities before and after. Five concern the explicitly exploratory `WarmAgainstFlat_AlongTheLength` readings (two at the own seeds, three on the chain axis). The other two are the two-ridges E8 recovery claim's mean-C and mean-subspace angles: each holds on one of two chain seeds, while its relative-C error holds on both and all three data-axis claims hold on both. These are existing limits to the comparative claims, not newly introduced failures, and are left for the owners of those statistical tests. The sweep test reports these rates rather than asserting that every comparison holds.

`RandomEffectsTests.AHeldScaleFit_KeepsItsDraws(legendre)` deliberately moves `a1f4528886afa5f7` → `f07a33dd35774c18`, because its modified g-prior on whitened data reaches PriorPencil, whose unscaled relative solve changes rounding. The archived baseline reproduces all three former pins (`artifacts/spectral-b1/pin-baseline/pin-baseline.trx`). The hinge and tempered-hinge pins are unchanged.

The reflection probe in `artifacts/spectral-b1/pin-probe` calls the exact private pinned fixture and fit in the two saved test builds. Across all 400 retained draws, term sequences, tau and sigma² agree exactly. Maximum absolute coefficient difference is 2.8422e-14; group effects differ by at most 2.5536e-14. The old/new JSON draws and `comparison.json` retain the evidence. This is a rounding change in continuous draws; no structure or scale trajectory change was observed on this fixture. The current test keeps its pinned trait, with a dated reason next to the new expectation. The second moved pin is qualified below.

`BassBitIdentityTests.Smc_KeepsTheBaseline(khaos)` also deliberately moves `a51939ff36591784` → `3505662fb879bb6a`, for the modified-g PriorPencil arithmetic. The archived baseline reproduces both former SMC pins (`pin-baseline/smc-baseline.trx`); the hinge pin stays. The same reflection probe reconstructs the exact KHAOS SMC settings (300 particles, two replicates, three mutation steps, master seed 13). Over 600 retained draws, term sequences and tau are identical; sigma² differs by at most 4.4409e-16, coefficients by 6.5053e-14, weights by 8.1966e-17 and replicate log normalizers by 1.4211e-14. Replicate temperatures, sweep counts and genealogy errors are unchanged. JSONs and `pin-probe/smc-comparison.json` retain this evidence. Three noise pins are qualified below.

Three `RandomEffectsNoiseTests.ANoiseFitWithoutABlock_KeepsItsDraws` pins also move through the weighted PriorPencil:

| Fit | Former pin | Current pin |
|---|---|---|
| t | 5471c4bcad7d0b09 | 0d6c539efe10a3bd |
| skewed horseshoe | afd1aa3f399ea35f | f5855f7382b9826b |
| tempered normal-wald | 98b9820e14563b71 | 8410c6ea70239f29 |

The archived baseline reproduces all three former pins (`pin-baseline/noise-baseline.trx`). The reflection probe reproduces each method's noise, priors and sampling settings and exports 300 retained draws per fit. All term sequences and tau values are unchanged. Continuous latent updates amplify the initial numerical differences:

| Fit | Max absolute sigma² gap | Coefficient gap | Latent-value gap | Skew gap |
|---|---:|---:|---:|---:|
| t | 2.36e-15 | 1.61e-12 | 1.24e-12 | 0 |
| horseshoe | 2.39e-9 | 1.06e-7 | 5.91e-8 | 1.05e-8 |
| normal-wald | 4.07e-12 | 7.88e-8 | 2.38e-10 | 9.93e-11 |

These are observed differences, not a claimed rounding-level bound on entire latent trajectories. The JSONs and `pin-probe/noise-comparison.json` retain the comparison. The existing missing `Tier=Pinned` trait is added to this method, so its stored expectations enter future commit gates. Five pins in total move; none of their former values is asserted as a contract of the new solver arithmetic. The exact kernel/pencil gates, full regression contracts and before/after comparative sweep qualify the change separately from these hashes.

The new checks use analytic spectra and current-path comparisons, not stored incidental output. Self-review covered forwarding of the rule through every overload, threshold equality and zero handling, local normalization, final-state convergence and stalling, reference dispatch, eigenvalues-only parity, production removal of ScaledEigen calls, preservation of absolute caller arithmetic, and exclusion of concurrent changes. The final diff and owned-path source hashes were checked against the candidate. Code commits are `ea76fdf` and `13ef477`. No outstanding code fix remains in this slice; the representability and statistical limits above remain explicit qualifications.

## Doc changes owed

- `docs/Maths/LinearAlgebra.md`: named stopping rule, the staged absolute default, relative tolerance units, final-state convergence and `Stalled`, and the reference build's scalar relative branch and shared primitives.
- `docs/Maths/Regression/ActiveSubspaces.md`: replace the unit-scale workaround paragraph with the direct element-wise relative solve and this slice's bounded E8 evidence; final deletion of ScaledEigen remains a later slice.
- `docs/Maths/Regression/BASS.md` (the pencil) and `MultiOutput.md` (the shared generalized spectrum): name the direct relative solve; prior whitening and RidgeSpectrum's statistical column normalization remain caller policy.
- No `docs/` files were edited by this worker.
