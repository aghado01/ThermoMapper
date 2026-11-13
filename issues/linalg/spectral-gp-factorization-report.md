# F — GP solve-based factorization

## Decision and verification — 2026-11-12

**Adopted.** Ordinary GP fits now keep a solve-based `DenseCholesky`; leave-out constructs selected precision blocks per unit. The historical raw-pivot admission **strictly above 1e-12** is preserved, including rejection of small positive pivots. No shared factor implementation, tensor spline, or eigensolver code changed in this slice.

**Production/test commit:** `5db2e405ac835624f5ba7f7ef79c05691e205326` — `perf(bgp): solve fits and form leave-out precision per unit`. Landed directly on shared `main`; unrelated guidance, drafts, issues workspace files and concurrent B1 changes were preserved. All three owned C# files match the qualified snapshot by SHA-256.

- **Commit gate:** 152 passed, no failures/skips; 17 selected classes plus pinned tier, [Commit.trx](../../artifacts/spectral-gp/gates/test-select/20261112_151338-commit/Commit.trx).
- **Whole regression assembly:** 3,280 passed, no failures/skips; 24m04s test time / 1,452s wrapper wall time, [Slice.trx](../../artifacts/spectral-gp/gates/test-select/20261112_151605-slice/Slice.trx).
- **Independent numerical checks:** existing wide-prior, decimal variance and subset-evidence tests passed. All six probe fixtures passed fit, prediction and singleton/grouped leave-out comparison in both final processes. Maximum relative difference: **7.8515e-14**, below the stated 1e-10 tolerance.
- **Representative n = 128/384, both final processes:** ordinary fits **2.83–3.52× faster**, **48.5–49.5% fewer allocated bytes**; public leave-out **2.17–2.88× faster**, **32.1–36.4% fewer bytes**.
- **Small n = 32 control:** fits **1.65–2.56× faster**, 43.6% fewer bytes; leave-out **1.27–2.06× faster**, about 25% fewer bytes. Every final fixture/operation median favored the candidate. Absolute small-case timings remain visibly noisier than allocation counts.

These measurements justify adoption on the measured sizes, rather than retaining the eager inverse-factor path. Timings and allocations are separated for fitting and leave-out in the table below; no claim of sampler-wide speedup is made.

## Execution and load notes

Started from `e3a3b7f597f98d7a71f9eaf3489fe903ad4d22b2` on shared main. Qualification and probes use an isolated archive plus the three GP-owned C# files, excluding concurrent B1 changes, with slice-specific output under `artifacts/spectral-gp/`.

Two normal wrapper attempts stopped during restore (exit 1, no useful diagnostic). Explicit single-node restore succeeded with NU1900 (vulnerability-feed data path denied); the test build then succeeded with 64 warnings and no errors, and the normal capped below-normal wrapper ran with `-NoBuild`. Probe builds succeeded; warnings were copied-source/pre-existing documentation warnings. One accidentally overlapping duplicate probe build was promptly cancelled. None of those build calls is timing evidence. Git's sandbox initially denied `.git/index.lock`; the production commit succeeded through the requested escalated Git action.

The first exploratory timing round started after a CPU-filtered host check showed no active dotnet work, but B1 builds began concurrently. It was stopped; partial `benchmark-round1.jsonl` is **not acceptance evidence**. Host Nushell sees workloads omitted by sandbox `Get-Process`, so subsequent checks used host `ps --long`, selecting process identity, CPU and command only. Coordination state was kept in this report; no cross-chat messages were sent and lightweight source work continued.

GP's full gate completed about 15:40 PST. B1's full gate and subsequent commit gate remained active, so GP timing waited until they exited about 15:57. Quiet preliminary rounds 2/3 then exposed insufficient small-fixture warmup; the final protocol below corrected that before the decision. **Accepted rounds 4/5 ran about 16:00–16:04 PST**, after the probe-only rebuild. Before/during/after snapshots showed no test hosts or other active .NET build/test workload; during snapshots showed only the GP probe at about one core. Normal desktop activity remained: sampled summed process CPU was roughly 140–318% (100% is one core; 18 logical processors). This is a shared desktop measurement, not an idle dedicated benchmark host. Load receipts are `quiet-before-round4.nuon`, `load-round4*.nuon`, `quiet-before-round5.nuon`, `load-round5*.nuon` and `quiet-after-round5.nuon`.

## Candidate and qualification method

`GpFactoredFit` now keeps `DenseCholesky`. The caller admits only finite positive-definite factors whose **raw** minimum pivot exceeds `1e-12`. Thus this does not silently adopt DenseCholesky's weaker positive-only admission. The two covariance solves, GLS algebra and evidence formulas are unchanged; the packed factor's back-solve changes accumulation order, so comparison uses tolerances.

`GpPrecision` is constructed only by leave-out. It copies the packed lower factor once, forward-solves the nonzero tails of the inverse-factor columns needed by one unit, and takes their dot products to form `R⁻¹_UU`. Columns are released after each unit. The partition requests each column once; off-block precision entries and a full covariance inverse are never formed. Runtime remains cubic for all-unit scoring, with reduced arithmetic and storage for small units. Ordinary fitting constructs neither inverse nor inverse-factor columns. No shared LinAlg API or factor implementation changes.

New `BgpFactorizationTests` cover:

- known and unknown scale at `BitDecrement(1e-12)`, exactly `1e-12`, and `BitIncrement(1e-12)` using a diagonal covariance;
- a nearly rank-one two-point covariance whose second positive Schur pivot lies below/above the historical boundary;
- an unordered selected precision block against solves of the covariance's coordinate vectors.

Existing `BgpConstantMeanTests`, `BgpVarianceTests` and `BgpPointwiseTests` supply independent wide-prior, decimal-elimination and full-minus-subset evidence checks. The latter covers singleton, consecutive and interleaved units under both noise models.

The standalone probe compiles the **unchanged baseline GP class** from `git show e3a3b7f:src/maths/regression/bgp/GpRegression.cs`, renamed `LegacyGpRegression`, beside the candidate in one process. Both use the same archived dependencies. It compares mean, alpha, Q, log evidence, mean/variance predictions and predictive log densities, then singleton/grouped leave-out terms at `1e-10 * max(1, |a|, |b|)` tolerance. Fixtures: n = 32, 128, 384; deterministic two-dimensional inputs `x=(u, 0.3 sin(5u))`, `u=i/(n−1)`; `y=3+sin(8u)+0.13 cos(37u)`; squared-exponential t = 0.02, known noise variance or unknown-scale nugget 0.04, no input scaling. Unit shapes: singletons, consecutive fours, interleaved fours and four consecutive quarters. GP construction and distance precomputation are outside the measured calls, identically on both sides.

Measurement protocol: Release .NET SDK 10.0.201/runtime 10.0.5, Windows 10.0.26220, Intel Core Ultra 5 125H (18 logical processors), default tiering/PGO, below-normal probe process. Each operation gets at least 64 alternating warmup pairs **and at least one second of warmup**, then nine rounds with alternating baseline/candidate order and full GC before each pair. Actual warmup pair counts are stored per operation. Repetitions per round: 60/8/1 for n = 32/128/384. Record per-call milliseconds and `GC.GetAllocatedBytesForCurrentThread` deltas. Leave-out calls include their public API's refit; the initial fit passed to that API is outside its timing. Ordinary fits are measured separately. The helper alone is not substituted for the end-to-end leave-out measurement.

Preliminary quiet rounds 2 and 3 used 64 warmup pairs only. Round 2's n = 32 known-model leave-out medians temporarily favored legacy (0.82–0.98×), while the otherwise identical covariance work later in the unknown-model rows favored the candidate (1.64–1.87×). The small fixture completes its initial warmups in milliseconds; asynchronous tiered optimization is a plausible explanation. Rather than select a favorable row, the final protocol adds elapsed-time warmup and reruns both sides. These earlier receipts are retained and labeled preliminary; adoption uses only elapsed-warmup rounds 4/5. At n = 128/384, round 2 already favored candidate fits by 2.81–3.50× and leave-out by 2.16–2.99×, with consistently lower allocations.

Probe sources and local receipts: `artifacts/spectral-gp/probe/`, `dotnet-info.txt`, `source-hashes.txt`, build logs, `gates/test-select/`, and load snapshots. Source snapshot is baseline `e3a3b7f597f98d7a71f9eaf3489fe903ad4d22b2` plus the three owned C# files. Concurrent main changes are excluded. The Slice wrapper's printed changed-path census sees the live checkout, but Slice imposes no filter and executes this already-built isolated assembly.

## Final measurements

Accepted raw receipts: [round 4](../../artifacts/spectral-gp/benchmark-round4.jsonl), [round 5](../../artifacts/spectral-gp/benchmark-round5.jsonl). Both contain all 30 operations and all six correctness passes. [summary-final.json](../../artifacts/spectral-gp/summary-final.json) retains each process's medians, minima and per-call allocation counts; `probe/Summarize.ps1` computes them. Allocation counts were constant across all nine trials and both processes for every row. The speedup range is the ratio of each process's legacy/candidate medians. Exact timings fluctuate between processes, but the direction holds throughout.

| n | Model | Call | Legacy median ms, runs 4 / 5 | Solve median ms, runs 4 / 5 | Median speedup range | Bytes/call, legacy → solve |
|---|---|---|---|---|---|---|
| 32 | known | fit | 0.0288 / 0.0300 | 0.0175 / 0.0140 | 1.65–2.15× | 26008 → 14664 |
| 32 | known | leaveout-blocks4 | 0.0580 / 0.0515 | 0.0384 / 0.0388 | 1.33–1.51× | 39640 → 29784 |
| 32 | known | leaveout-interleaved4 | 0.0783 / 0.0619 | 0.0491 / 0.0407 | 1.52–1.59× | 39640 → 29784 |
| 32 | known | leaveout-quarters | 0.0509 / 0.0519 | 0.0265 / 0.0310 | 1.68–1.92× | 39960 → 30008 |
| 32 | known | leaveout-singleton | 0.0349 / 0.0441 | 0.0276 / 0.0214 | 1.27–2.06× | 37080 → 27800 |
| 32 | unknown | fit | 0.0347 / 0.0284 | 0.0136 / 0.0144 | 1.98–2.56× | 26008 → 14664 |
| 32 | unknown | leaveout-blocks4 | 0.0469 / 0.0598 | 0.0244 / 0.0323 | 1.85–1.92× | 39640 → 29784 |
| 32 | unknown | leaveout-interleaved4 | 0.0653 / 0.0441 | 0.0390 / 0.0300 | 1.47–1.67× | 39640 → 29784 |
| 32 | unknown | leaveout-quarters | 0.0486 / 0.0436 | 0.0369 / 0.0328 | 1.32–1.33× | 39960 → 30008 |
| 32 | unknown | leaveout-singleton | 0.0476 / 0.0744 | 0.0280 / 0.0386 | 1.70–1.93× | 37080 → 27800 |
| 128 | known | fit | 1.6369 / 1.8535 | 0.5776 / 0.6235 | 2.83–2.97× | 397720 → 204744 |
| 128 | known | leaveout-blocks4 | 2.5523 / 2.5165 | 1.1193 / 1.0480 | 2.28–2.40× | 550360 → 363288 |
| 128 | known | leaveout-interleaved4 | 2.6160 / 2.2515 | 0.9099 / 1.0357 | 2.17–2.88× | 550360 → 363288 |
| 128 | known | leaveout-quarters | 2.2838 / 2.3541 | 1.0512 / 1.0347 | 2.17–2.28× | 584856 → 397112 |
| 128 | known | leaveout-singleton | 2.4809 / 2.4365 | 0.9715 / 0.9308 | 2.55–2.62× | 540120 → 355352 |
| 128 | unknown | fit | 2.0787 / 1.7280 | 0.7116 / 0.5991 | 2.88–2.92× | 397720 → 204744 |
| 128 | unknown | leaveout-blocks4 | 3.3675 / 3.7161 | 1.3046 / 1.2915 | 2.58–2.88× | 550360 → 363288 |
| 128 | unknown | leaveout-interleaved4 | 2.3850 / 2.8858 | 1.0031 / 1.1704 | 2.38–2.47× | 550360 → 363288 |
| 128 | unknown | leaveout-quarters | 2.8251 / 3.7201 | 1.2831 / 1.6064 | 2.20–2.32× | 584856 → 397112 |
| 128 | unknown | leaveout-singleton | 2.7802 / 2.9260 | 1.1014 / 1.1643 | 2.51–2.52× | 540120 → 355352 |
| 384 | known | fit | 37.4917 / 41.2408 | 11.6482 / 11.9261 | 3.22–3.46× | 3551640 → 1792968 |
| 384 | known | leaveout-blocks4 | 65.0944 / 58.7626 | 23.6136 / 23.0056 | 2.55–2.76× | 4795864 → 3054872 |
| 384 | known | leaveout-interleaved4 | 58.7246 / 60.9445 | 23.0571 / 21.8514 | 2.55–2.79× | 4795864 → 3054872 |
| 384 | known | leaveout-quarters | 70.4913 / 56.4813 | 29.0413 / 25.3656 | 2.23–2.43× | 5191832 → 3448632 |
| 384 | known | leaveout-singleton | 60.6404 / 55.3873 | 21.8687 / 19.4441 | 2.77–2.85× | 4765144 → 3031064 |
| 384 | unknown | fit | 46.2146 / 38.9611 | 13.3482 / 11.0532 | 3.46–3.52× | 3551640 → 1792968 |
| 384 | unknown | leaveout-blocks4 | 55.9740 / 54.2408 | 21.6475 / 21.0046 | 2.58–2.59× | 4795864 → 3054872 |
| 384 | unknown | leaveout-interleaved4 | 56.2222 / 61.9239 | 21.8142 / 24.4324 | 2.53–2.58× | 4795864 → 3054872 |
| 384 | unknown | leaveout-quarters | 59.4853 / 59.8496 | 26.2276 / 26.7289 | 2.24–2.27× | 5191832 → 3448632 |
| 384 | unknown | leaveout-singleton | 58.6836 / 52.5685 | 20.9872 / 19.0348 | 2.76–2.80× | 4765144 → 3031064 |

## Second-pass review

Re-read the changed GP call chain and independently checked the triangular-column algebra, tail indices, symmetry, packed-index bounds, observation-partition coverage and unknown-scale outside-unit guard. `UnitPartition` supplies unique in-range members; no column is reused across units, so keeping columns only for a current unit is valid. Both raw-pivot rejection and failed/nonfinite factor rejection remain caller policy. Prediction and joint-posterior callers continue solving through the same fit-local immutable factor, so no shared mutable cache was introduced. No code findings. This was a second-pass self-review under the no-subagent instruction, not a claim of an external independent review.

## Commands and receipts

Run from `D:/aghado01/ThermoMapper`:

```powershell
# Explicit single-node fallback after the normal wrapper's restore failures.
dotnet restore artifacts/spectral-gp/source/tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj -p:ArtifactsRoot=D:/aghado01/ThermoMapper/artifacts/spectral-gp/gates/ --disable-build-servers -m:1
dotnet build artifacts/spectral-gp/source/tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj -c Release -p:ArtifactsRoot=D:/aghado01/ThermoMapper/artifacts/spectral-gp/gates/ --disable-build-servers --no-restore -m:1
pwsh -NoProfile -File scripts/test-select.ps1 -Scope Commit -SourceRoot D:/aghado01/ThermoMapper/artifacts/spectral-gp/source -Changed src/maths/regression/bgp/GpRegression.cs,src/maths/regression/bgp/GpPrecision.cs,tests/regression/BgpFactorizationTests.cs -ArtifactsRoot D:/aghado01/ThermoMapper/artifacts/spectral-gp/gates -Threads 4 -NoBuild
pwsh -NoProfile -File scripts/test-select.ps1 -Scope Slice -SourceRoot D:/aghado01/ThermoMapper/artifacts/spectral-gp/source -ArtifactsRoot D:/aghado01/ThermoMapper/artifacts/spectral-gp/gates -Threads 4 -NoBuild
dotnet build artifacts/spectral-gp/probe/Probe.csproj -c Release -p:ArtifactsRoot=D:/aghado01/ThermoMapper/artifacts/spectral-gp/probe-build/ --disable-build-servers -m:1
# Final probe-only rebuild also used --no-restore -p:BuildProjectReferences=false.
dotnet artifacts/spectral-gp/probe-build/bin/Probe/Release/net10.0/Maths.Regression.Tests.dll
```

`commit-gate-nobuild.log`, `slice-gate.log` and `gate-summary.json` record the passing gates. `test-restore.log`, `test-build.log`, `probe-build.log` and `probe-build-final.log` and `probe-build-settled.log` record builds; `commit-gate.log` and `commit-gate-retry.log` retain the failed wrapper attempts. `environment.json` records that `SEED_SWEEP`, `BASS_QUALIFICATION`, `DOTNET_TieredCompilation` and `DOTNET_TieredPGO` were unset. The standalone single-node fallback build's inherited priority was not recorded; both passing test runs and all timing calls explicitly used below-normal priority.

## Limits

This slice measures fixed-hyperparameter fitting and the public leave-out API, not complete sampler throughput or a hyperparameter optimization. It does not establish costs beyond n = 384, a hardware-independent speedup, or subcubic all-unit leave-out. Ordinary fits remain O(n³) dense Cholesky and all-unit leave-out remains O(n³). Numerical comparison is tolerance-based, not a promise of bit-identical trajectories. Full-assembly qualification is of the isolated GP snapshot, not combined main with concurrent B1 changes. Probe code, source snapshot and generated receipts are local ignored artifacts; the report and production/test code are the reviewable committed record.

## Pins moved

None.

## Doc changes owed

[BGP.md](../../docs/Maths/Regression/BGP.md) should state that ordinary fits keep a solve-based packed dense factor, GP retains its raw-pivot `> 1e-12` admission, and leave-out forms selected precision blocks from per-unit inverse-factor columns. Its algebra and cubic all-unit complexity stay the same. Documentation is deferred for the offline batch, as commissioned; `docs/` is untouched.

[LinearAlgebra.md](../../docs/Maths/LinearAlgebra.md), the factor-admission paragraph around line 464, should distinguish BAPS's floored-factor rejection from BGP's explicit minimum-pivot admission on the strict dense factor; add BGP to that factor's consumers.
