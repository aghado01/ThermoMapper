# Performance cost of the regression engines, the samplers and the shared factorizations (chip brief)

Date: 2026-11-04. A measurement, tests only: nothing under `src/` changes and no default moves. It follows [the BASS cost brief](bass/performance-cost-brief.md), whose counted units, reference rate, quiet-machine protocol and lever table it reuses, and it covers what that brief left out. The efficiency pass for BASS is a separate chip (`bass-efficiency-pass`) that builds BASS's already priced levers; this one prices everything else. Its report goes under "Reports" below, newest first.

## Why

The idioms the repo uses elsewhere (`Vector<T>`, `Vector256`, `TensorPrimitives`, `ArrayPool`, `AggressiveInlining`, struct-constrained generics and `static abstract` as in `ISwConfig`, `IRiemannianManifold` and `IRobustLoss`) are absent from `maths/regression` and `maths/samplers`, which hold one struct-constrained generic, `Tally<TCategory>`, an enum. The inner loops of the shared factorizations are scalar too (`DenseCholesky.cs:79`, `sum -= _l[ri + k] * _l[rj + k]`). That the idioms are absent is known. Where the time goes is not, because only BASS has a cost report.

What a read of the code shows, before any timing, and what the timing must confirm or refute:

- *BARS and ABMRS.* One chain step evaluates the target once. `SplineTarget.Evaluate` builds a dense `double[m, ν]` design through `SplineBasis.Design`, which calls `BSpline.EvaluateBasis` per row, and that allocates four arrays per row. It then scans the dense design for the support gate (`HasSupport`) and again for each row's span (`new BandedDesign(design)`). `WeightedNormalModel.Resolve` clones and sorts the response to find its median on every evaluation, although the response never changes. The band Cholesky itself is O(ν·p²) at p = 3.
- *BAPS.* Every Gibbs step takes the noise block's residual sum of squares through a dense n × ν loop (`Baps.ResidualSS`), beside a band factorization of O(ν·p²).
- *BGP.* Every step builds an n × n kernel matrix with an interface call and an `exp` per pair, then `CholeskyDecomposition.Decompose`, which also forms L⁻¹ in a second O(n³/6) pass that `FitCore` never reads, and allocates three n × n arrays per fit.
- *Changepoint.* The exact DP spends three transcendentals per `LogSumExp` cell, over O(n²·k²) cells in the regression-curve loop.
- *Samplers.* `ReversibleJumpChain.Step` calls `TotalWeight` twice and `Weight` for each move again, looks the reverse move up by string key, and allocates an `ITargetEvaluation` per proposal. `SliceSampler.Step` and `RandomWalkMetropolis.Step` take a `Func<double, double>`. `SamcLearning.Update` builds a new `MulticanonicalWeights` per update. `SmcRun.Advance` forks and mutates through interfaces under `Parallel.For`; `BassSmcKernel.Move` builds a `ReversibleJumpChain` (a dictionary) per move and takes a lock (`BassSmc.cs:161`).

## The decisions

1. **Measure; build nothing.** No production file changes. A lever's price is a test-side replica compiled fully optimized, beside the production code at the same state, never a change to `src/`. A replica that changes bits says so.
2. **Counts at once, times in the quietest window available.** Allocation bytes, pause counts, analytic operation counts and realized sizes do not change under load. Wall times follow the protocol below.
3. **Parts by replica, wholes by production.** A part internal to a method is timed through the public pieces it calls, or through a verbatim test-side copy of a private loop; the whole is the production method. Parts weighted by their frequencies must sum to at most the whole, and the rest is reported.
4. **M1's timing method:** a warm-up batch, then the median over nine batches of a batch's mean per call, the parts' batches round-robin with the whole's so load falls on all alike, one thread, Release. A batch is calibrated to 5 ms (the whole's to 20 ms). A call over 100 ms is a *heavy point*: five single-call batches, round-robin all the same.
5. **Tier-0 guard:** every timed block runs once at a tenth, its rows discarded, before its record pass (the BASS cost run found its first readings partly tier-0 code).
6. **Budget:** up to 3 hours of timed machine time, waits excluded, projected from a 1/10 trial, with a declared trim order.
7. **The dynamic PGO that the runtime enables by default stays on** for every reading except the settings block, which prices it.
8. **Frame nothing around SPC.** The engine serves diverse applications; every size below is chosen from the engines' own fixtures and declared ranges.

## The measurements, fixed before any run

Declare every value in the harness classes before the first timed run. A change after a run is a deviation, reported with its reason.

### The reference rates (`machine`)

Denominators for every ns per counted unit. Each read single-threaded, warmed, at n = 1,000, 10,000, 100,000 and 1,000,000 where a length applies:
- a scalar multiply-add loop; a `Vector<double>` loop; `PairwiseSum.Dot`; `MatrixOps.Dot`; an unrolled four-accumulator `Vector256` dot;
- memory streams (a copy and a triad) at 1, 2, 4 and 8 threads;
- `Math.Exp`, `Math.Log`, `Math.Sqrt` and a division, scalar, and the in-box `Vector256.Exp` and `Log` per element (`TensorPrimitives` is built on them; the test project does not reference `System.Numerics.Tensors`, and adding it would edit its project file, which this brief does not);
- one `Func<double, double>` call, one interface call (monomorphic and bimorphic), one `Dictionary<string, T>` lookup, one `Lock` enter-exit (uncontended), one `new double[n]` for n = 8, 64, 512 and 4,096.

### The shared factorizations (`kernels:*`)

Each read as time per counted unit against the reference rate, and against a fixed-order vector reference built test-side. The counted unit of a factor is the multiply-adds its loops execute, counted by running the loop's bounds without the arithmetic (and asserted equal to the closed form where one exists). The reference is a four-lane dot with `PairwiseSum`'s fixed reduction order, so its bits do not depend on the machine's SIMD width, on a layout where the dot's operands are contiguous. Bit identity is read as: the maximum relative difference of the factor against production, and whether every entry is identical.

| Block | Points | Production | References |
| --- | --- | --- | --- |
| `kernels:band` | ν ∈ {8, 16, 32, 64, 128, 256, 512, 1,024} × p ∈ {1, 2, 3, 5, 8, 16, 32, 64} (p < ν), LDLᵀ strict and Cholesky floored | `BandCholesky.DecomposeBanded` + `SolveInPlace` | (a) scalar, contiguous row-major band; (b) compile-time p by `static abstract`, p ∈ {1, 2, 3, 5}, unrolled; (c) fixed-order 4-lane dot on the row-major band |
| `kernels:engines` | the bands the engines build: BARS ν = 8, 20, 36 at p = 3; BAPS ν = 29, 104, 504 at p = 3; BAPS tensor g × g for knots 6 and 10 at the bandwidth read from `TensorProductDesign` and `TensorPenalty` | as above | as above |
| `kernels:dense` | m ∈ {5, 10, 20, 50, 100, 200, 500, 1,000} | `DenseCholesky.Decompose`, `ForwardSolveInPlace`, `BackSolveInPlace` | (a) scalar over the same packed layout; (b) fixed-order 4-lane dot per row pair |
| `kernels:bgp` | n ∈ {50, 100, 200, 400, 800, 1,600} | `CholeskyDecomposition.Decompose` (factor and L⁻¹ passes timed apart by replica), `Solve` | `DenseCholesky` strict on the same matrix; its fixed-order vector replica |
| `kernels:wide` | (n, p) ∈ {(200, 30), (500, 60), (1,000, 100), (500, 250)} | `WideBandCholesky.Decompose`, `Solve`, with `BandCholesky` and `DenseCholesky` on the same matrix | none |
| `kernels:bordered` | M ∈ {5, 10, 20, 50, 100, 200} | `ProposeBirth`, `Accept`, `ProposeDeath`, `ProposeChange`, `Rebuild`, `Solve` | units from BASS's table (1.5M² per proposal, M³/3 per rebuild, 2M² per death); the dense row dot at the same lengths |
| `kernels:lowrank` | (m, k) ∈ {20, 50, 200} × {10, 39, 100} | `DiagonalPlusLowRank` construction (with `ThinSvd.Decompose` timed apart) and `SolveInPlace` | none (priced in BASS's calibration; here read per counted unit) |
| `kernels:khatri` | n_x ∈ {100, 500} × m ∈ {10, 50, 200, 1,000} at M = 20 factors | `KhatriRao.Gram`, `Border`, `CentredBorder`, `CrossProducts` | none |

### The engines' per-iteration split (`engine:*`)

Each engine is read at a burned state (2,000 steps from its default start on its fixture), its parts at that state, and its whole as the production method over 500 further steps. Frequencies per step come from the production counters. A part's counted unit and the reference it is read against are in its row.

**`engine:bars`** (also ABMRS). Fixture: y = sin(2π x) + 0.15·ε, x at the cell centres, Poisson(5) prior, `LocalBetaKernel(50)`, `WeightedNormalModel` defaults. ABMRS is BARS with `EbicPrior(γ, C)` on the same engine (no engine of its own exists in `src/`; `docs/Maths/BayesianRegression.md` states the mapping); the fixture reads it at γ = 1, C = 100 and C = 1,000, whole step only, plus the prior's `LogPrior`.
- Points: n ∈ {100, 500, 2,000, 10,000}, at the chain's realized k; and, at those n, forced equispaced configurations of k ∈ {3, 8, 16, 32} for the evaluation's parts.
- Parts of `SplineTarget.Evaluate`, each a unit count and a rate: the knot simplex check; the design (`SplineBasis.Design`: n·4 `EvaluateBasis` row evaluations, n·ν stores, 4 arrays allocated per row); the support gate (n·ν adds); `new BandedDesign(design)` (n·ν reads); the response's centre (`Resolve`: a clone and an n log n sort); `Accumulate` (n·(p+1)(p+2)/2 multiply-adds); the factor (ν·p²/2); admission (`DesignAdmission.Admits`: pivots and the condition estimate's solves); the centring (`Centre`: `Apply`, a row-sum scan, an array); the right-hand side and solve; the fit and its residual; `CollapsedGPrior.FromLeastSquares`, `MarginalAt` and `LogMarginal`.
- Parts of a step outside the evaluation: `TotalWeight` (×2), `SelectMove`, the move's `Propose` (birth, death, relocate separately), the reverse-key lookup, the two `Math.Log` of the selection ratio, `Temperature.LogDensity`, the uniform and the `Math.Exp`. The rest is the whole less the parts.
- Then the *run*: `Bars.Run` with one chain, burn-in 500 and 1,500 samples at n = 500, reading ns per iteration of the whole run (weights, readout and draws included), against the steps' sum.

**`engine:baps`**. Fixtures: the one-dimensional P-spline of `BapsTests` (n points on a cubic basis with equispaced knots, `DifferencePenalty(2)`), and the tensor surface of `BapsAnisotropicTests` (g × g, `TensorPenalty`).
- Points: 1-D (n, ν) ∈ {(120, 29), (1,000, 104), (5,000, 504)}; tensor (n, ν) ∈ {(900, 100), (2,500, 196)}.
- Schemes: `Collapsed` (the λ walk through `RandomWalkMetropolis` over the REML evidence), `Gibbs`, `LocallyAdaptiveGibbs` for 1-D; `AnisotropicGibbs` and `Gibbs` for the tensor.
- Whole: `Baps.Run` with one chain, burn-in 200, 600 samples; ns per iteration and bytes per iteration.
- Parts: the Gram copy; the penalty's accumulation; the band factor (units ν·p²/2); `GaussianConditional.Draw` (two solves and a draw); the noise block's `ResidualSS` (n·ν multiply-adds, and the banded `Apply` equivalent, n·(p+1)); the roughness; the conjugate draws; for the collapsed scheme the REML evidence per λ proposal (copy, accumulate, factor, solve, log-determinant) and its acceptance frequency.
- `HierarchicalPSpline` (the BAPS extension with 53 allocation sites) is read at its default hierarchical scheme on 8 replicates of the 1-D fixture at (n, ν) = (120, 29): whole ns per iteration, allocation bytes per iteration, and the shares of its blocks from `BlockStats`-style per-block timing by a test-side scan.

**`engine:bgp`**. Fixture: `BgpBenchmark`'s Swiss-roll inputs in d = 3, f = sin(2πu) + ½cos(2πv), noise 0.1, `SquaredExponentialKernel`, known σ² = 0.01, similarity scaling.
- Points: n ∈ {50, 100, 200, 400, 800, 1,600}.
- Whole: `GpRegression.Fit(t)` (one step's target) and `BgpSampler.Run` with one chain, burn-in 100 and 400 draws (ns per iteration, including its per-draw conditional variance at 100 test points, and bytes per iteration).
- Parts of `Fit`: the kernel fill (n(n+1)/2 `Kernel.Evaluate` interface calls, each an `exp`, plus 2 stores per pair), `Decompose` step 1 (n³/6 multiply-adds), step 2 (n³/6 again), the two solves (n² each), the reductions; the three n × n allocations.

**`engine:changepoint`**. Fixture: piecewise-constant means, noise sd 0.5, `GaussianSegmentModel`, with `ExactChangepoint.DefaultHyperparameters`.
- Points: (n, k) ∈ {(100, 5), (100, 20), (300, 5), (300, 10), (300, 20), (1,000, 5), (1,000, 10), (1,000, 20)}.
- Whole: `ExactChangepoint.Infer`; then `ChangepointPosterior.Sample(1,000 draws)`.
- Parts by replica: `model.Fill` (n²/2 cells), the forward DP (k·n²/2 `LogSumExp` cells), the backward DP, the break probabilities (n·k²/2 cells), the regression-curve loop (n²/2 · k²/2 cells), the range-add and readout. The unit is a `LogSumExp` cell, read against two `exp` plus one `log` from `machine`.

**`engine:gprior`**. The collapsed g-prior block both BARS and BASS use: `FromLeastSquares`, `MarginalAt` (documented to allocate nothing), `LogMarginal`, `At`, `Sample`, at M ∈ {5, 20, 50, 100, 200} columns; bytes per call.

### The samplers (`sampler:*`)

**`sampler:rj`**: `ReversibleJumpChain<TState>.Step` on two targets.
- *Cheap:* a birth-death chain on k with target ∝ λᵏ/k! (an `ITargetEvaluation` per proposal, three moves of weights 0.4, 0.4, 0.2 with the third a self-reverse relocation).
- *Expensive:* BARS's target at n = 500, k at its realized value.
- Parts, as replicas of `Step`: `TotalWeight` once and twice; `SelectMove`; the reverse-key dictionary lookup against a precomputed reverse index; the evaluation (allocation included); the acceptance arithmetic. Whole: the production `Step` with PGO on.
- A struct-constrained generic replica of the whole `Step`, on the cheap target only, with `static abstract` palette and target: what it adds over the production `Step` at the same state. The `Func<double, double>` and interface sites are read with dynamic PGO on (default) and, in the settings block, off.

**`sampler:blocks`**: `SliceSampler.Step` and `RandomWalkMetropolis.Step`, each on a cheap density (−x²/2) and two expensive ones (BAPS's REML evidence at ν = 104; `GpRegression.Fit` at n = 200), against test-side struct-generic replicas of the same algorithm. Read: ns per step, evaluations per step, the delegate's share, bytes per step.

**`sampler:pt`**: `ParallelTempering.Step` at 5, 10 and 20 rungs, on the cheap target and BARS at n = 100. Parts: the rungs' steps (the cheap/expensive chains' costs times the rungs), `Record`, `Attempt`, `TrackEnds`. The swap round's share, in ns per rung and bytes per scan.

**`sampler:smc`**: `SmcRun.Advance` on a tempered Gaussian target (prior N(0, I), likelihood N(μ, 0.5²I) in d = 2 and 10), with a cheap `Evaluate` and an expensive one (a 2,000-term `exp` sum), at N = 1,000 and 4,000 particles, `RandomWalkKernel`. The stage is repeated at one temperature, so the cost parts are separated by their options: `MutationSteps` 0, 1, 2 and `ResampleBelow` 0 (never) and 1 (every stage). Reweight, resampling (ancestors and forks) and mutation per particle; `Parallel.For`'s fan-out overhead at 1, 4 and 8 degrees of parallelism. Then `BassSmcKernel.Move` at the BASS cost harness's base point: ns per move, bytes per move, the `ReversibleJumpChain` construction timed apart, `Account.Add` timed apart, and the lock's contended cost at 4 and 8 threads (a loop of `lock (gate) Account.Add`).

**`sampler:samc`**: `SamcLearning.Update` at K ∈ {50, 200, 1,000, 5,000} nodes and 4 and 16 chains: ns per update, bytes per update (`WithLogWeights` copies), and the weights' `LogWeight` and `TryLocate` per lookup; `MulticanonicalEvaluation.Reweighed` per chain; a multicanonical chain's `Step` against the plain chain on the cheap target.

**`sampler:wham`, `sampler:ladder`, `sampler:ensemble`**:
- `Wham.Estimate` on a synthetic `TemperingPilot`, rungs ∈ {8, 16, 32} × draws per rung ∈ {2,000, 20,000, 200,000}.
- `LadderEvidence.SteppingStone` and `PathSampling` on rungs ∈ {10, 30, 60} × draws per rung ∈ {2,000, 20,000, 200,000}.
- `ChainDiagnostics.Diagnose` on 4 chains × {1,000, 10,000, 100,000} draws × 1 and 10 functionals; `ChainEnsemble.Run` of an independent-draw model (the engine's own overhead per draw) at 4 chains × 20,000 samples.

### Allocation (`alloc:*`)

Counts, run at once: bytes per iteration by `GC.GetAllocatedBytesForCurrentThread()` deltas and, in the class run alone, by `GC.GetTotalAllocatedBytes(precise: true)`; the same parts as above, so the allocations that run per iteration are told from the sites that do not (the counts of `new double[` are sites). Then the GC cost, one chain and `c` chains in parallel (c = 1, 2, 4, 8; each chain an independent instance, `Parallel.For`): `GC.GetTotalPauseDuration()` over wall time, the collection counts by generation, and the parallel speedup, for BARS (n = 500), BAPS Gibbs (1-D, ν = 104), BGP (n = 200 and 800), PT (10 rungs, cheap), SMC (4,000 particles) and the exact changepoint (n = 300, k = 10). The BASS report's bound is a pause of at most 6.5% of a sweep; the question is whether chains in parallel change it.

### The memory-bound regimes (`membound`)

A point is memory-bound when its working set leaves the caches (this machine: 2 MB L2 per P-core, 24 MB shared L3). Each is read as ns per unit from a size inside the L2 to one far past the L3, and the rate is set beside the stream rate:
- BAPS `ResidualSS` and the banded `Apply`: n·ν·8 bytes ∈ {0.1, 1, 8, 64, 512} MB.
- BARS's design build and gate at n·ν·8 from 0.1 to 256 MB (n up to 200,000 at ν = 16).
- BGP's n × n kernel fill and the Cholesky: n ∈ {200, 400, 800, 1,600, 3,200}.
- The changepoint's tables, three (n+1)² arrays: n ∈ {300, 1,000, 2,000, 3,000}.
- `KhatriRao.CrossProducts` at BASS's direct sizes.

Parallel efficiency at each size: the same instance replicated on 1, 2, 4 and 8 threads, throughput against the single-thread rate. A point whose footprint (times its threads) exceeds a third of the free physical memory is skipped and said so. This machine had 2.3 GB free of 31.5 GB when the brief was written.

### The runtime settings (`settings`)

Priced on one engine, BARS's step at n = 500 (with BGP's `Fit` at n = 400 and the BAPS Gibbs iteration at ν = 104 as two further readings of the same settings, since they allocate differently). Each setting runs in its own process, by environment variable on the `dotnet test` command; a row records the setting it ran under.

| Setting | Environment |
| --- | --- |
| default | none (workstation, concurrent GC, tiered compilation, dynamic PGO) |
| no PGO | `DOTNET_TieredPGO=0` |
| no tiering | `DOTNET_TieredCompilation=0` |
| quick jit for loops off | `DOTNET_TC_QuickJitForLoops=0` |
| non-concurrent GC | `DOTNET_gcConcurrent=0` |
| server GC | `DOTNET_gcServer=1` |
| server, non-concurrent | `DOTNET_gcServer=1`, `DOTNET_gcConcurrent=0` |
| larger gen0 budget | `DOTNET_GCgen0size=0x4000000` |

Each reads ns per iteration, bytes per iteration, pause time per iteration, and the time to the first reading at the default's rate (the warm-up's length), one chain and four chains in parallel.

## The levers' table

Each lever's row states: the site; the idiom (element-wise `Vector256`, fixed-order reduction, struct-constrained generic or `static abstract`, `ArrayPool` or reused scratch, data layout, or an algorithmic change); the share of the engine's iteration it could remove, per engine and at the axes' ends, as the measured part's share less the part's rate floor (the reference rate, where an arithmetic loop has one); and **whether it keeps bits**.

The bits rule: element-wise lanes without FMA match scalar bits; reductions and FMA change them. `PairwiseSum`'s fixed order gives the same bits on every machine; `MatrixOps.Dot`'s order depends on the machine's SIMD width. A layout change that keeps every operation and its order keeps bits. Skipping exact zeros keeps bits; a different association does not. A vectorized `exp` differs in the last place.

Candidate levers, priced from the shares and not built:

| Site | Idiom | Bits |
| --- | --- | --- |
| `SplineBasis.Design` and `BSpline.EvaluateBasis`: four arrays per row, a dense design | the band built straight into band storage, scratch reused | kept if the arithmetic is |
| `SplineTarget.HasSupport`, `BandedDesign`'s span scan: dense n·ν scans | column mass and spans from the band (n·(p+1)) | kept |
| `WeightedNormalModel.Resolve`: a median sort per evaluation | the response's level once per target | kept |
| `DesignAdmission.Admits`: the condition estimate's solves | none priced unless it is a large share | — |
| `BandCholesky` at p = 3 | compile-time bandwidth (`static abstract`), unrolled | kept |
| `BandCholesky`, `WideBandCholesky` and `DenseCholesky` at long rows | fixed-order 4-lane dot over a contiguous layout | changes (a reduction) |
| `CholeskyDecomposition`: the unused L⁻¹ pass; three n × n arrays per fit | factor only (`DenseCholesky`'s path), scratch reused | the factor kept |
| BGP's kernel fill | `TensorPrimitives.Exp` on the row | changes (last place) |
| `Baps.ResidualSS`: a dense n·ν loop | the banded `Apply` | kept |
| `ExactChangepoint`'s `LogSumExp` cells | one maximum per cell and one `log`; vectorized `exp` | changes |
| `ReversibleJumpChain.Step`: `TotalWeight` ×2, the string-keyed lookup | one pass, a precomputed reverse index | kept |
| `SliceSampler.Step`, `RandomWalkMetropolis.Step`: `Func<double, double>` | struct-constrained generic | kept |
| `SamcLearning.Update`: a new weights object per update | in-place update with a version | kept |
| `SmcRun.Advance` and `BassSmcKernel.Move`: a dictionary per move, a lock per move | a chain reused per worker; per-thread accounts merged between stages | kept |
| Allocation generally | `ArrayPool`, reused scratch | kept |
| Runtime settings | `TieredPGO`, concurrent and server GC | kept |

## The quiet-machine protocol

The BASS cost run's protocol, with the one change this machine forces.
- **Counts first.** Run every count (allocation, operation counts, realized sizes, the checks) at once.
- **Check before each timed block.** The harness records the other `testhost` processes at the block's start and end, their priority classes, the other processes' CPU in cores busy over the block (the processor-time deltas of every readable process, less its own), the free physical memory, and its own CPU over wall.
- **Two tiers, declared now.** *Q (quiet):* no other `testhost` on the machine. *B (background):* the only other test hosts run at below-normal priority (the waiting-time and basin-hold chains, which may hold 4 threads for hours) and the block's other-process load is recorded. Wait for Q with the Monitor tool's until-loop, not a sleep, for at most 30 minutes before the first block; then proceed in B. A B block is flagged, never reported as quiet, and re-timed under Q at the run's end if a Q window then exists (a sample of at least 25 points over the blocks). Single-thread readings are read under B to about ±30% (the BASS run's spread under 1.6–4.2 busy cores); **parallel-efficiency and memory-bandwidth readings under B measure a shared machine and are reported as such, or not at all**.
- **Priority.** On this machine other work holds below-normal threads on most of the 18 logical processors (the first look found the machine at about 80% CPU), so a Normal thread is timesliced against it. Each timed pass runs at `AboveNormal`, restored afterwards, and the block's rows say so. The tier and the priority are both reported; the A/B reading at `Normal` (`COST_PROFILE_PRIORITY=Normal`, its own rows) prices what the bump changes on the `machine` block and one engine.
- **If no Q window comes by the end,** the report says which blocks stand under B only, with the loads, and the Q commands are written down.
- **Record the machine:** the commit, Release, the .NET runtime, the CPU model and the logical processor count, the GC mode, tiered-compilation settings, and the active power plan (`powercfg /getactivescheme`, read only).

## The budget

- **Trial first.** A 1/10 trial, every block at a tenth, and a projection.
- **Up to 3 hours of timed machine time,** waits excluded.
- **Above that, trim in this order and record it:** (1) repeats 9 → 5; (2) the changepoint's (1,000, 20); (3) BARS n = 10,000; (4) the tensor (2,500, 196) and BGP n = 1,600; (5) the memory-bound points' largest sizes; (6) 8-thread rows; (7) the settings' `gen0size` and `quick jit` rows.
- **No need to ask:** the user is away; the brief's rules decide.

## Checks before the runs

- *The timer:* the empty part's time (the harness's overhead), reported and subtracted nowhere.
- *The allocation counter:* a `byte[1000]` reads 1,024 B; `CollapsedGPrior.MarginalAt` reads 0.
- *Counted units:* each analytic count equals the loop's executed count (the counting replica) at three sizes.
- *Replicas against production:* each replica's result equals the production result (bit-identical where the arithmetic is unchanged; the maximum relative difference where it is a different association), asserted.
- *Nothing commits:* a part's timed call leaves the production state untouched.
- *The parts against the whole:* the parts, weighted by their frequencies, sum to at most the whole within its IQR, and the rest is reported.
- *Repeatability:* the base point of `engine:bars` timed twice in one process agrees within its IQR on its parts, counts exactly.
- *Existing suites still pass* with the harness in the tree (the regression and sampler suites' fast filters).

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** No worktrees; switch no branch; push nothing.
- **Tests only.** New files under `tests/regression/` (the regression test project already references Regression, Samplers, LinAlg and Rng): `CostProfile.cs` (grid, meters, rows, the protocol), `CostProfileKernels.cs`, `CostProfileEngines.cs`, `CostProfileSamplers.cs`, `CostProfileTests.cs` (the facts), all behind `BASS_QUALIFICATION` = 1; with the switch off the facts return at once. Output goes to `COST_PROFILE_OUT` (a folder); `COST_PROFILE_SCALE` = 0.1 is the trial; `COST_PROFILE_BLOCK` names the blocks to run, separated by semicolons. Rows append and a rerun skips a block whose rows are present, so a restart loses nothing.
- **Build and test into** `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-cost-profile\`, set on each `dotnet` command itself. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-cost-profile\` when done.
- **Commit your own files only, by explicit path** (`git commit -- paths`), leaving nothing staged: the index is shared. End each message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Commit the harness before any timed run. Before each commit read `git diff --cached` (or the path-limited diff) and confirm every hunk is yours.
- **Other sessions' uncommitted work stays untouched,** including the untracked files and the deleted `issues/regression/opus-bass-completion-remainder.md` in the tree.
- **Your files:** the five test files, the per-point table `performance-cost-regression-samplers-results.csv` beside this brief, and this brief.
- **Line endings:** briefs and C# are LF.
- **No sub-agents.**

## Review before reporting

A fresh-eyes pass, fixed before the report:
- *Correctness:* each part's timing is of the work the engine does at that state (a replica matches production, asserted); allocation is read with nothing else allocating; every slope is fitted above the fixed overhead, with the points it rests on named; the analytic counts match the loops; no B block is reported as quiet.
- *Coverage:* an engine whose realized size missed its declared point; a part that is a large share somewhere but was timed at one point only; a question of the brief this run does not answer, named as missing with its reason.

## Report

Under "Reports", newest first:
1. each engine's split by part in ns per counted unit against the reference rate;
2. the factorizations against the fixed-order vector reference, with the real bandwidths and whether short bands gain from vectors at all or a compile-time bandwidth beats them;
3. the samplers' overhead on a cheap and an expensive target, with PGO on, and what a struct-constrained generic adds over it;
4. allocation (bytes per iteration, GC pause) at one chain and several in parallel, and which `new double[` sites run per iteration;
5. SMC's per-particle overhead, the account's lock included;
6. the memory-bound regimes and their parallel efficiency;
7. the runtime settings priced on one engine;
8. the lever table per engine and for the shared kernels, ranked by the share it could remove, with bits kept or not; a recommendation where the numbers make one clear; what could not be measured and why;
9. the machine, the load seen at each block, the tiers, the trims and the deviations; the checks, the review passes, and what the main session must know.

## Not in this brief

- BASS's own levers (`bass-efficiency-pass`), except where BASS exercises a shared kernel (`BorderedCholesky`, `KhatriRao`, `DiagonalPlusLowRank`, the collapsed g-prior block, `BassSmcKernel`).
- Statistical efficiency (error per unit time): costs here do not depend on mixing, and the statistical half multiplies them by an ESS.
- The rest of the repo (geometry, distributions, clustering, graphs): another day.
- Any lever's implementation, and any default.

## Reports

Newest first.

*(none yet)*
