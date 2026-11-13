# F — anisotropic tensor-spline factorization

## Decision — retain compact LDLᵀ (2026-11-12)

**Retain the current `AnisotropicTensorSpline` compact diagonal-major `BandCholesky(..., Ldlt)` path. No production code changed.** The actual `SpdFactors.ForBand` route needs dense input conversion and does not justify a blanket replacement:

- Across two warmed processes and nine layouts, complete-fit median speedups span **0.84–1.55×**, with allocated bytes **1.28–5.74× the baseline**. The 16×8 complete fit is about 15% slower in both runs; 6×16 changes direction between runs. This does not establish the proposed 2–4× end-to-end gain.
- Fixed-evidence calls span **0.76–1.75×** with allocated bytes **1.27–5.79× the baseline**. Faster factor/evidence calls do not reliably translate to a faster complete fit.
- The public fixed-smoothing API changes numerical admission at λx=λy=1e20 in three fixtures. Those points lie outside the optimizer's actual search boxes, but the API allows positive fixed smoothing parameters; the differences remain adverse qualification evidence. Neither path has an independent high-smoothing accuracy oracle here.

This retains existing behavior rather than making an unqualified public numerical-policy change for modest, layout-dependent time gains and higher allocation. It does not reject a future compact-input row-packed factor experiment; that would be a different API/implementation and needs its own measurements.

## Qualification and provenance

Source baseline: **`305a373791caf8209aa9a2617bd0593a45068442`**, archived under `artifacts/spectral-tensor/source`. The source archive excludes concurrent B1 stage-two changes. Final hash checks confirm the live tensor caller, penalty, LDLᵀ, wide/dense dispatcher and optimizer match the archived files. Unrelated guidance/drafts, issues workspace files and other workers' changes were preserved. All probe code and generated evidence remain under `artifacts/spectral-tensor/`.

- **Full baseline regression assembly:** **3,282 passed**, no failures/skips; 39m34s test time / 2,383s wrapper wall time, [Slice.trx](../../artifacts/spectral-tensor/build/test-select/20261112_163215-slice/Slice.trx). This is qualification of the unchanged baseline, not a claim that the probe-local candidate was wired into the whole assembly.
- **Focused tensor/anisotropic gates:** 7 passed, [tensor-baseline.trx](../../artifacts/spectral-tensor/gates/tensor-baseline.trx).
- **Dense/wide/strict factor gates:** 40 passed, [factor-baseline.trx](../../artifacts/spectral-tensor/gates/factor-baseline.trx). No failures/skips in either focus.
- **Candidate nominal and bounded-domain checks:** all 81 nominal fixed-smoothing cases, 225 actual-box grid cases and nine complete optima pass the stated gates. The union of both optimizer traces contains 990 distinct trial points, all cross-evaluated with zero admission mismatches. Three unidentified-nullspace optimizations and nine corresponding fixed systems refuse under both paths.

[qualification.jsonl](../../artifacts/spectral-tensor/qualification.jsonl), [qualification-summary.json](../../artifacts/spectral-tensor/qualification-summary.json) and [optimum-summary.json](../../artifacts/spectral-tensor/optimum-summary.json) preserve the detailed comparisons. Across the optima, maximum log-smoothing discrepancy is **8.0348e-7**, normalized objective discrepancy **1.2031e-14**, and coefficient discrepancy **7.9848e-8**. Both methods converge with matching bound flags. Evaluation counts can differ slightly (16×32: 112 versus 114); timings retain the actual work of each search.

The baseline timing uses the actual production class. The probe-local candidate mechanically copies that class, changes the factor return type to `ISpdFactor`, converts the accumulated compact lower band to a fresh dense lower triangle, invokes the actual `SpdFactors.ForBand`, and requires raw `MinPivot > 1e-12` to carry forward the intended floor rejection. No other evidence or optimizer formula changes. That intended threshold does not guarantee identical rounded pivots under LDLᵀ and square-root Cholesky; the observed public differences below demonstrate this.

Premeasurement numerical gates: coefficient error ≤1e-7, log determinant error ≤1e-9, REML error ≤1e-8 for fixed points; complete optima require convergence, matching bound flags, evidence error ≤1e-8, log-lambda discrepancy ≤1e-3 for interior optima and coefficient discrepancy ≤1e-3. Adversarial cases are retained separately, not silently dropped.

### Public smoothing calls versus the optimizer's box

At λx = λy = 1e20, **6×16 changes from acceptance to rejection**, and **16×16 and 32×16 change from rejection to acceptance**. These are genuine public fixed-smoothing behavior changes between LDLᵀ and the proposed strict square-root factor with the same nominal pivot threshold, and remain adverse evidence for adoption.

The actual REML search bounds come from each axis's trace-balance reference times [1e-8,1e8]. Across these nine fixtures, upper bounds are only about **1.59e7–1.73e8**, so 1e20 lies outside every actual box. The final qualification receipt records each fixture's exact lower/upper bounds (`kind=search-box`).

Inside the actual boxes, a 5×5 log-grid per fixture (225 points including every corner/edge midpoint) found **zero admission mismatches**. A replay of the exact public Nelder–Mead objective/search settings reproduced both optimizers' points, values and evaluation counts bit-for-bit; the union of their trial points contains **990 distinct points**, all cross-evaluated by both paths, also with **zero admission mismatches**. `kind=search-trace` records the reproduction checks and counts, and `kind=search-point` retains every comparison. This is bounded evidence on the grid and actual trial points, not a proof of equivalent admission throughout the continuous boxes. The outside-box public differences are not excused by the optimizer's narrower domain.

All 225 box-grid points also meet the fixed-point numerical gates: worst normalized coefficient difference **6.1484e-8**, REML difference **2.8810e-9**, log determinant difference **6.7524e-10**. At λx = λy = 1e12 (outside all boxes), both paths accept all fixtures, but worst differences grow to **4.1513e-5** (coefficients), **9.8779e-7** (REML) and **4.9855e-7** (log determinant). No independent high-smoothing oracle establishes which path is more accurate there.

All nine unidentified-nullspace fixed systems are refused by both paths. Their three optimizations all refuse with the same `ArgumentException: The start must be feasible`. The first probe stopped at this expected refusal because it did not catch it; the corrected probe records both refusals and completes. The initial receipt is retained separately. `qualification-summary.json` summarizes the complete final receipt, including the three λ = 1e20 admission mismatches.

## Fixture and measurement protocol

Every fixture uses cubic bases, second differences on each axis, open-clamped [0,1] splines, and `ν−4` interior knots at `k/(ν−3)`. The observation grid is `gx=max(20,2νx)`, `gy=max(20,2νy)`, at midpoints `x=(a+0.5)/gx`, `y=(b+0.5)/gy`. Responses are deterministic:

`sin(1.5πx) sin(5πy) + 0.1 cos(37x+53y) + 0.05 sin(89x−71y)`.

| Shape | Observations | Coefficients | Half-bandwidth | Actual proposed backend |
|---|---:|---:|---:|---|
| 4×4 | 400 | 16 | 15 | DenseCholesky |
| 6×8 | 400 | 48 | 27 | DenseCholesky |
| 6×16 | 640 | 96 | 51 | DenseCholesky |
| 8×8 | 400 | 64 | 27 | WideBandCholesky |
| 8×16 | 640 | 128 | 51 | WideBandCholesky |
| 16×8 | 640 | 128 | 27 | WideBandCholesky |
| 16×16 | 1024 | 256 | 51 | WideBandCholesky |
| 16×32 | 2048 | 512 | 99 | WideBandCholesky |
| 32×16 | 2048 | 512 | 51 | WideBandCholesky |

The baseline backend is **BandCholesky/Ldlt in every row**, reading diagonal-major compact storage directly. The proposed backend is read from the actual factor returned by `ForBand`, not inferred from the class's name. The candidate fills only the dense lower triangle because both dispatched factors read only that triangle; allocation/zeroing of the full dense array is included.

Two operations were measured independently:

- **Evidence:** `RemlLogEvidence(0.7,2.0)` on each preconstructed spline. This includes compact Gram cloning, penalty accumulation, all dense conversion in the proposed path, fresh factor allocation/setup, decomposition/log determinant, solve and evidence arithmetic.
- **Complete REML fit:** construct the spline from the shared `TensorDesign` and response, call `RemlOptimum(tolerance:1e-6,maxEvaluations:400)`, then fit coefficients at the returned optimum. Thus common Gram/RHS/spectrum setup, every trial's factor/conversion/solve/log determinant, and the final coefficient solve are included. Input tensor-design construction is common precomputation outside this operation.

Release .NET SDK 10.0.201/runtime 10.0.5; Windows 10.0.26220; Intel Core Ultra 5 125H, 18 logical processors; default tiering/PGO; below-normal probe priority. Warm each evidence operation for at least 64 alternating pairs and one second; warm each complete fit for at least two pairs and one second after the evidence loop has warmed the factor paths. Record actual warmup counts. Five trial pairs alternate execution order, with full GC before each pair. Evidence repetitions per trial: 20 through 64 coefficients, 5 through 128, 1 above that; complete fits use one. Timings use Stopwatch; allocations use `GC.GetAllocatedBytesForCurrentThread`. Two independent processes were run. Untimed `calibrationOnly` values in the qualification receipt are not benchmark evidence.

Normalized numerical errors are `|a−b|/max(1,|a|,|b|)`; coefficients take the maximum over coordinates. The whole-fit comparison also checks bound flags, convergence and the returned λ values. Both methods use the exact same input arrays, smoothing points, objective formulas, search settings and source dependencies.

## End-to-end measurements

Raw receipts: [round 1](../../artifacts/spectral-tensor/benchmark-round1.jsonl), [round 2](../../artifacts/spectral-tensor/benchmark-round2.jsonl). Each has nine fixture/backend records and 18 measurements. [summary-final.json](../../artifacts/spectral-tensor/summary-final.json) retains per-process medians, minima, warmup counts and allocation medians/minima/maxima. `probe/Summarize.ps1` computes it. The allocation-byte column below shows **run 1 medians**; every run's exact allocation ranges are in the summary/raw receipts, and small trial-to-trial variations are retained. Speedup is compact/ForBand, so values below 1 favor the retained path.

Absolute timings drift by roughly 3–5× between processes despite the same build/settings. CPU affinity, frequency and desktop background activity were not controlled; no cause is assigned to the drift. Therefore these values do not define a precise crossover or hardware-independent speedup. Within-process paired comparisons and the allocation/public-behavior findings support the retain decision; no universal time dominance is claimed.

| Shape | Operation | Compact median ms, runs 1 / 2 | ForBand median ms, runs 1 / 2 | Speedup range | Allocated bytes, run 1 median, compact → ForBand | Allocation ratio |
|---|---|---|---|---|---|---|
| 16x16 | evidence | 2.3947 / 0.5757 | 1.6512 / 0.4318 | 1.33–1.45× | 219400 → 745792 | 3.40× |
| 16x16 | reml-fit | 288.8957 / 79.4824 | 212.2622 / 55.0591 | 1.36–1.44× | 22685392 → 76405576 | 3.37× |
| 16x32 | evidence | 17.4199 / 4.6693 | 12.7065 / 3.0748 | 1.37–1.52× | 831752 → 2933056 | 3.53× |
| 16x32 | reml-fit | 2420.7127 / 642.3783 | 1738.6915 / 449.1052 | 1.39–1.43× | 95073936 → 338426336 | 3.56× |
| 16x8 | evidence | 0.3704 / 0.0876 | 0.3044 / 0.0762 | 1.15–1.22× | 60680 → 192832 | 3.18× |
| 16x8 | reml-fit | 52.6439 / 15.8891 | 60.7701 / 18.3518 | 0.87–0.87× | 7379440 → 23234832 | 3.15× |
| 32x16 | evidence | 6.2997 / 1.2069 | 3.5957 / 0.9046 | 1.33–1.75× | 438536 → 2539840 | 5.79× |
| 32x16 | reml-fit | 729.2084 / 186.9066 | 648.2957 / 140.9338 | 1.12–1.33× | 49263808 → 282613760 | 5.74× |
| 4x4 | evidence | 0.0122 / 0.0031 | 0.0091 / 0.0025 | 1.22–1.34× | 4744 → 6008 | 1.27× |
| 4x4 | reml-fit | 0.9811 / 0.3392 | 0.9525 / 0.3287 | 1.03–1.03× | 539440 → 689352 | 1.28× |
| 6x16 | evidence | 0.6464 / 0.1608 | 0.6848 / 0.2121 | 0.76–0.94× | 82440 → 154296 | 1.87× |
| 6x16 | reml-fit | 91.3309 / 21.2848 | 81.1656 / 25.2352 | 0.84–1.13× | 8529984 → 15856448 | 1.86× |
| 6x8 | evidence | 0.1040 / 0.0376 | 0.0959 / 0.0260 | 1.08–1.45× | 22920 → 40440 | 1.76× |
| 6x8 | reml-fit | 17.3922 / 7.1312 | 15.1795 / 6.5508 | 1.09–1.15× | 2773632 → 4857888 | 1.75× |
| 8x16 | evidence | 1.1469 / 0.3288 | 0.6869 / 0.1983 | 1.66–1.67× | 109832 → 241984 | 2.20× |
| 8x16 | reml-fit | 114.9565 / 31.1657 | 84.3987 / 23.7280 | 1.31–1.36× | 10806496 → 23622392 | 2.19× |
| 8x8 | evidence | 0.1769 / 0.0344 | 0.1693 / 0.0283 | 1.05–1.22× | 30472 → 63808 | 2.09× |
| 8x8 | reml-fit | 26.9329 / 5.0324 | 17.3706 / 4.1534 | 1.21–1.55× | 3498248 → 7264440 | 2.08× |

## Load and execution receipts

Heavy work used single-node builds and capped below-normal tests, with slice-specific output. B1's original full gate overlapped only untimed qualification and this slice's baseline gate. Tensor timing waited until no test hosts remained at 17:12 PST. Round 1's metadata starts **17:13:01 PST** and the receipt finished **17:14:05**. A second pre-start host command had a Nu syntax error; the corrected sample at **17:13:32** showed no tests or build drivers and only the tensor probe active, about one core, with summed process CPU about 193% (100% is one core; 18 logical processors). This is finite sampled evidence on a shared desktop, not a dedicated idle host.

The next sample was much later, at **18:35**, and showed a later B1 build. It is named `late-after-round1-first.nuon`, not presented as a during-timing sample. At **18:36** B1 started a later commit gate, so round 2 waited for it. Round 2 started about **18:43:10** and finished **18:43:39**. The pre-start check at 18:43:08 showed no tests, build drivers or active dotnet workload. The 18:43:31 during sample showed only the tensor probe and summed process CPU about 211%; the 18:43:59 post sample showed no active dotnet workload (summed CPU about 72%). No tensor timing overlaps the observed B1 heavy work. Receipts: `load-round1-start.nuon`, `late-after-round1*.nuon`, `quiet-check-round2.nuon`, `load-round2-early.nuon` and `load-round2-after.nuon`. No cross-chat messages were sent; coordination was recorded in this report and lightweight work continued.

The first qualification run stopped at the expected unidentified-nullspace optimizer refusal because the probe had not caught it. The corrected probe records both paths' refusals and completes; earlier receipts are preserved as `qualification-first.jsonl` and `qualification-before-search-box.jsonl`. No production adjustment was made to force a comparison to pass.

## Review and limits

Second-pass self-review checked the generated candidate against the archived caller, lower-band indexing and dense conversion bounds, actual size/bandwidth dispatch, threshold admission, and measurement boundaries. Reflection and backend inspection are used only outside timing. The baseline timing calls the production class directly. Candidate source differs only in its factor return type, compact-to-dense conversion, `ForBand` call/admission and an untimed inspection entry; the evidence and optimizer formulas are unchanged. No independent subagent was dispatched.

The representative family is cubic with second-difference penalties and regular observation grids. It does not establish a universal crossover for every degree, penalty order, irregular design or hardware. The public λ = 1e12/1e20 comparisons have no independent accuracy oracle; baseline acceptance at an extreme smoothing is not certification that its numerical answer is correct. These regimes expose a qualification gap, rather than permission to discard the adverse cases. The bounded-search evidence covers the listed grid and actual optimizer traces only.

## Build identity and commands

The source archive is `305a373791caf8209aa9a2617bd0593a45068442`. Because the archived projects sit below the live repository, MSBuild stamps the surrounding checkout's revision into their informational version: **1.0.0+4e9547f6a630c594a7c4c1d4713e2cd72a727f35**. This metadata is not represented as the frozen source revision. `build-versions.json`, `build-hashes.txt`, `source-hashes.txt` and `source-hashes-final.txt` retain binary and source evidence. The live caller, penalty, LDLᵀ, wide/dense dispatch and optimizer files were hash-checked against the archive and match. The probe's generated candidate and `Prepare.ps1` are retained, so its factor-only changes can be reviewed.

Commands run from `D:/aghado01/ThermoMapper` (logs under `artifacts/spectral-tensor`):

```powershell
git archive 305a373 -o artifacts/spectral-tensor/baseline.zip
Expand-Archive -LiteralPath artifacts/spectral-tensor/baseline.zip -DestinationPath artifacts/spectral-tensor/source
./artifacts/spectral-tensor/probe/Prepare.ps1
dotnet build artifacts/spectral-tensor/probe/Probe.csproj -c Release -p:ArtifactsRoot=D:/aghado01/ThermoMapper/artifacts/spectral-tensor/build/ --disable-build-servers -m:1
# Probe-only rebuilds added --no-restore -p:BuildProjectReferences=false.
dotnet artifacts/spectral-tensor/build/bin/Probe/Release/net10.0/Maths.Regression.Tests.dll
dotnet artifacts/spectral-tensor/build/bin/Probe/Release/net10.0/Maths.Regression.Tests.dll --timing
pwsh -NoProfile -File scripts/test-select.ps1 -Scope Slice -SourceRoot D:/aghado01/ThermoMapper/artifacts/spectral-tensor/source -ArtifactsRoot D:/aghado01/ThermoMapper/artifacts/spectral-tensor/build -Threads 4 -NoBuild
```

Test projects were restored and built single-node with the same artifact root; the two direct focused gates used `probe/gate.runsettings` (four xUnit threads), `--no-build`, trx logging and the class filters recorded in `baseline-tests.log` / `factor-tests.log`. All builds/runs were below normal. The initial probe build needed an explicit parent-namespace import for its renamed class; corrected probe builds succeeded with zero warnings/errors. The regression test build had 22 warnings and no errors; the LinAlg test build had one warning and no errors. Restore logs retain any NU1900 vulnerability-feed access warning. The redundant zip was removed after extraction; the source archive directory, probe sources and receipts remain.

## Pins moved

None.

## Doc changes owed

`docs/Maths/LinearAlgebra.md`, the wide-band selection paragraph: distinguish its dense-input kernel comparison against `DenseCholesky` from the measured end-to-end comparison against an already compact LDLᵀ caller. Record this caller's retain decision and conversion/allocation cost from this report. No public API or factor-dispatch change needs documenting. The coordinator owns the shared brief's completion link. `docs/` stays untouched.
