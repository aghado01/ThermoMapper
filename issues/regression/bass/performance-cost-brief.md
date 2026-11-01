# Performance qualification, the cost half (chip brief)

Date: 2026-10-31. Status: briefed for a chip; the decisions below are the main session's recommendations, taken when the user starts the chip. This brief is the hub's performance qualification ([its section](../bass-full-port-brief.md#validation-and-performance-qualification)) for what a run costs; what it buys, statistical error per unit time and the run-length defaults, is the statistical half, briefed after Friedman's replication. Its report goes under "Reports" below, newest first. `docs/Maths` and the checklists are the main session's.

## Why

The hub asks for benchmarks by case count, predictor count, active terms, interaction order, response size and rank, and category cardinality. It asks for allocation and peak memory, factor work, basis and cross-product scans, sensitivity moment work, and statistical error per unit time.

What exists is partial:
- *S10's cost baseline* ([BayesianRegression.md](../../../docs/Maths/BayesianRegression.md), "Cost"): µs per sweep over n ∈ {200, 1,000, 5,000}, M_max ∈ {5, 20, 50}, p ∈ {2, 10}, one thread, untempered (`BassQualificationCostTests`). Its builds were not recorded, so its numbers compare only within one run.
- *M1's decomposition* (`BassMixingCostTests`, F-M2's Phase 0): where a sweep's time goes, part by part, on the same grid, with the parts' frequencies per sweep.
- *Single readings:* N15's 1.4 ms per σ² proposal at m = 20 with 39 columns (stage 7); R BASS at about 630 µs per sweep on Friedman at n = 500 against C# near 100 (S10); F-R1's per-arm seconds.

Two things split the work in two:
- *F-R1 showed that untempered BASS does not converge at the default lengths* on any of its benchmarks: R̂ medians of 2.14 and 1.74, about 5 effective draws in 8,000. An error per unit time for the default chain would measure a stalled chain. That half needs samplers that mix and a reference that stands, and Friedman's reference rests on one seed until the warm-start chip's replication lands.
- *Costs do not depend on mixing.* A sweep's time, its allocation, its factor work and its scans are properties of the state and the data, and they are what the statistical half multiplies by an ESS.

## The decisions (recommended; starting the chip takes them)

1. **Cost only.** This brief measures and builds no optimization. It changes no default and recommends none. The efficiency pass's levers are priced from the measured shares, not built.
2. **No production code.** The harness is test-side. A part internal to a move is timed as M1 timed it: through the internals the test host already sees, or by a test-side replica compiled fully optimized. If a reading truly needs a production hook, stop and say so in the report rather than adding one.
3. **Counts at once, times on a quiet machine.** Allocation bytes, retained bytes, realized M, and parts' frequencies per sweep do not change under load, so they run at once. Wall times wait for a machine with no other test host running (the protocol below).
4. **One factor at a time around a base point,** plus S10's crossed n × M, extended. A full factorial over seven axes is out of budget; interactions beyond n × M are read where a cost model predicts one.
5. **M1's timing method:** a warm-up batch, then the median over nine batches of a batch's mean per call, with the parts' batches run round-robin with the sweep's, so load falls on every part alike. One thread per chain. Release.
6. **R as context at three points,** not as a comparison. R's sweep is not the same work: its fixed ⅓ move probabilities (F-O6) and its own evaluation order. The hub's warning stands: no speedup is assumed from the language.
7. **A budget of 6 hours of timed machine time,** projected from a 1/10 trial, with a declared trim order.

## The measurements, fixed before any run

Declare every value below in the harness class before the first timed run. A change after a run is reported as a deviation, with its reason.

**The surfaces.**
- *S10's surface* (`BassQualificationCostTests.Data`): x₁ and x₂ carry the signal and the other predictors are inert. It reads the proposal's cost in p with M held by the data.
- *A rich surface* with every predictor active (for example a sum of one- and two-way sinusoids over all predictors, at noise sd 0.05), so that the cap sets M. Use it wherever the axis is M or interaction order.
- Record the realized M at every point (the mean over retained draws), and use it as the x-axis where the cap is not reached. S10 found p = 10 short of the cap at n ≤ 1,000.

**The base point:** scalar response, n = 1,000, p = 10, M_max = 20, maximum interaction order 2, D01–D06's defaults.

**The axes,** each varied with the others at the base:

| Axis | Values | Notes |
| --- | --- | --- |
| Case count n | 200, 500, 1,000, 2,000, 5,000, 10,000, 20,000, 50,000 | rich surface |
| Active terms (M_max) | 5, 10, 20, 50, 100, 200 | at n = 1,000 and n = 10,000; rich surface |
| Predictors p | 2, 5, 10, 20, 50, 100, 200 | both surfaces |
| Interaction order | 1, 2, 3, 4, 5 | rich surface |
| Category cardinality L | 2, 4, 16, 64, 256, 1,024 | two categorical predictors beside five continuous, n chosen so every level has at least 4 cases |
| Direct functional | m ∈ {10, 50, 200, 1,000} × n_x ∈ {100, 500} | p = 5, M_max = 20 |
| Basis (PCA) rank | r ∈ {1, 3, 10, 30} | n = 500, m = 200 |

Categorical subsets are level lists, not a 64-bit word, so L past 64 is in scope. Record what each subset's membership test costs per case.

**What each point reads.**
- *Time:* µs per sweep, and M1's parts: the birth's draw, ratio, column, cross products, bordered Cholesky and design; the death's and change's proposals; the target; the three commits; the strict rebuild, amortized over its cadence; τ's block with its evaluation count; the composition; the rest.
- *Frequencies per sweep:* births, deaths, changes, admitted proposals, rebuilds, τ evaluations, columns evaluated.
- *Factor work and scans:* each counted operation's analytic cost (for example O(M²) per bordered update, O(n·M) per cross-product scan, O(M³) per rebuild) beside its measured time, as time per counted unit. That shows which parts run at the machine's arithmetic rate and which are overhead-bound.
- *Memory:*
  - steady-state allocation in bytes per sweep after burn-in, by `GC.GetTotalAllocatedBytes(precise: true)` deltas, with the class run alone so nothing else allocates in the process;
  - which parts allocate, where a sweep's allocation is not zero;
  - retained bytes per fit and per retained draw, against n, M, m and the draw count (the draws keep their full state, by the retention decision);
  - peak working set at each axis's largest point, each read in its own process.

**The outputs,** per draw, at the base and along M, p and order:
- prediction per case per draw, scalar and functional;
- sensitivity: Sobol indices and effects per draw (`BassSensitivity.Of`, `BassSensitivity.Effects`), and for basis fits against r (`BassBasisSensitivity.Of`).

**The other backends,** per unit, at the base point, beside the untempered sweep:
- *tempering:* per rung-sweep and per swap round, at 5, 10 and 20 rungs;
- *SMC for fitting:* per particle mutation step, per reweighting and per resampling, at 1,000 and 4,000 particles;
- *the flat histogram:* per iteration, its weight lookup and update included.

**Calibration** (N15): `DiagonalPlusLowRank`'s construction and solve per σ² proposal at m ∈ {20, 50, 200} × k ∈ {10, 39, 100}, and the cut SMC's cost per leg at stage 7's sizes. Beside them, the operation count of N15's capacitance-Cholesky alternative, O(m·k² + k³), against the per-call Jacobi thin SVD's: priced, not built.

**Parallel fits** (context): at three points, the default fit with its chains in parallel against one chain's time × the chain count, for the machine's parallel efficiency.

**R** (context): R BASS 1.3.1 from `r/` at three points: Friedman at n = 500, p = 10; the rich surface at n = 5,000, p = 10, M_max = 20; direct functional at n_x = 100, m = 50. Read µs per sweep over R's default lengths, with `sessionInfo()` recording the BLAS in use.

## The quiet-machine protocol

- **Counts first.** Run every count (allocation, retained bytes, frequencies, realized M, the checks) at once; load does not change them.
- **Wait for quiet before each timed block.** Start one only when no `testhost` process other than your own run's is on the machine. Wait with the Monitor tool's until-loop, not a sleep inside a test.
- **Record the load around each block.** The harness records the other test hosts it can see at each block's start and end, and its own CPU over wall. A block that saw another test host is flagged and re-timed at the end of the run, once the machine is quiet.
- **What you'll be waiting on.** The warm-start chip may still be running its Friedman step when you start; its runs are the wait. The main session holds the flat histogram's loose ends and its own suite runs until your timed grid finishes.
- **If no quiet window comes,** report the counts and the trial, and leave the timed grid for a later run with its commands written down. Times taken under load are not the grid.
- **Record the machine:** the commit, Release, the .NET runtime, the CPU model, the logical processor count, the GC mode, tiered-compilation settings and the active power plan (`powercfg /getactivescheme`, read only).

## The budget

- **Trial first.** Run a 1/10 trial (a tenth of the sweeps at every point), and project the full run from it.
- **Up to 6 hours of timed machine time,** waits excluded.
- **Above that, trim in this order and record it:**
  1. the repeats from 9 to 5;
  2. n = 50,000;
  3. p = 200;
  4. L = 1,024;
  5. m = 1,000;
  6. the second M row (n = 10,000).
- **No need to ask:** the user is away overnight, and the brief's rules decide.

## Checks before the runs

- *The timer:* an empty part's time, which is the harness's own overhead, subtracted nowhere and reported.
- *The allocation counter:* a known allocation reads its size; `CollapsedGPrior.MarginalAt`, documented as allocating nothing (N12), reads 0.
- *Nothing commits:* the timed proposals leave the chain untouched (the test-side `BassMixing.Untouched`), as in M1.
- *The parts against the sweep:* the parts, weighted by their frequencies, sum to at most the sweep, and the rest is reported.
- *Repeatability:* one grid point timed twice in one process agrees within its batches' interquartile range, and its counts agree exactly.
- *The surfaces:* the rich surface reaches the cap at the M points; L's designs give every level at least 4 cases.
- *Draws unchanged:* the eight bit-identity hashes (`BassBitIdentityTests`) pass with the harness in the tree.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Tests only.** The harness is `tests/regression/BassCostQualificationTests.cs`, behind `BASS_QUALIFICATION` = 1. Test-side helpers go under `tests/regression/`. The R context script is `r/scripts/bass_cost_context.R`.
- **Other chips may work beside you,** the flat-histogram warm-start chip among them ([its brief](flat-histogram-warmstart-brief.md)), in `src/maths/samplers/flathistogram/`, `BassFlatHistogramRun.cs` and the flat-histogram tests. Build against the checkout as it is. If a commit changes an API your harness uses, adapt your harness; never edit another chip's files. If a build fails mid-change in a file you do not own, wait and retry.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-perf\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-perf\` when you finish.
- **Your files:** the harness, its helpers, the R script, the per-point table `performance-cost-results.csv` beside this brief, and this brief's Reports section.
- **Commit your own files only, by explicit path,** in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Commit the harness before any timed run. The user's untracked notes and other sessions' edits stay unstaged. Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - each part's timing is of the work the chain does at that state, not a replica's;
  - allocation is read with nothing else allocating in the process;
  - every slope is fitted above the fixed overhead, with the points it rests on named;
  - the analytic operation counts match the code's loops;
  - no flagged block is reported as quiet.
- *Coverage:*
  - an axis whose realized M missed its cap;
  - a part that is a large share somewhere but was timed at one point only;
  - a hub measure this run does not read, named for the statistical half or as missing.

## Report

Under "Reports" below, newest first:
- **the scaling first:** per axis, µs per sweep, the slope above the fixed overhead with its range, and the parts' shares at the base and at the axis's ends;
- **memory:** allocation per sweep (and by part, where nonzero), retained bytes per fit and per draw, peak working set;
- **factor work and scans:** each counted operation's frequency, analytic count and time per unit;
- **the outputs, the other backends and calibration,** N15 priced;
- **the efficiency pass's levers, ranked** by the share of a sweep each could remove, at the base point and at the axes' ends. The candidates:
  - prefix-sum cross products for a main-effect hinge (the local-move mixing brief);
  - `KhatriRao.CentredBorder`'s re-centring allocation per direct proposal (BayesianRegression's direct algebra);
  - N15;
  - the state's allocation per sweep;
  - the rebuild cadence;
  - anything the shares show that this list misses;
- **the parallel and R context;**
- **the hub's list:** each axis and measure, what covers it, and what is left for the statistical half;
- **the machine, the trims, the flagged blocks and any deviation;**
- **the checks with their tolerances, the review passes, and anything the main session must know.**

Recommend no default. Price the levers; do not build them.

## Not in this brief

- **The statistical half:** error per unit time per backend, with method-appropriate diagnostics (particle-weight ESS and autocorrelation ESS are not one currency), and the run-length defaults that D10 left as a budget. It is briefed after Friedman's replication, on targets with references that stand.
- **The efficiency pass itself,** which acts on this report's ranking.
- **The completion audit** (deferred by the user).

## Reports

Newest first.

### The cost half measured (chip, 2026-10-31)

**What landed** (on `main`, nothing pushed; tests and an R script only, no production file):
- `fa39a1d`: the harness, written before its timed runs. `BassCostQualificationTests.cs` (the facts), `BassCostQualification.cs` (`BassCost`: the surfaces, the declared grid, the burn-in rule, the meters, the rows), `BassCostSweep.cs` (M1's reader generalized), `r/scripts/bass_cost_context.R`.
- Revisions, each before the readings it affects: `c4caee4` (a restored chain's column-mean cache filled; GC per sweep), `b3ee63d` (effects on a 20-point grid), `e216f8e` (retained memory after a warm-up fit), `b484fca` (each timed block runs itself at a tenth first, discarded; the change's analytic count), `abe8a7c` and `ca729db` (the machine's reference rates).
- `50650e6` repairs my error: `fa39a1d` swept in the warm-start chip's staged `flat-histogram-warmstart-brief.md`. `50650e6` put HEAD's copy back without touching the working tree, and that chip committed its report itself (`3990c63`). I told it. Every commit since names its paths explicitly.
- `performance-cost-results.csv`, beside this brief: 83 rows (59 grid points, then the outputs', backends', calibration's, parallel and machine readings) and 288 columns. `count:` columns come from the counts pass, `time:` from the record pass, `re-run:` from the grid's second pass.

**How to run it** (Release, behind `BASS_QUALIFICATION` = 1, rows to `BASS_COST_OUT`; `BASS_COST_SCALE` = 0.1 for the trial):
```bash
BASS_QUALIFICATION=1 BASS_COST_OUT=<dir> dotnet test tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj -c Release --filter "FullyQualifiedName~BassCostQualificationTests.Counts_EveryPoint"
```
- `Timed_Blocks` reads the blocks `BASS_COST_BLOCK` names, separated by semicolons: the ten axes by name (`n`, `M_max at n = 1,000`, … `n × M_max`), then `repeat`, `outputs`, `backends`, `calibration`, `parallel` and `machine`. Wait for no other `testhost` before each.
- `Counts_RetainedMemory` runs once. `Counts_PeakWorkingSet` runs once per process, with `BASS_COST_PEAK` set to each of `n`, `M_max`, `p`, `order`, `L`, `functional` and `basis rank`.
- The R context: from `r/`, `Rscript scripts/bass_cost_context.R <BASS_COST_OUT>/r <out>`.

#### The scaling

A point is read at chain 0's burned state, started and seeded as the run does it: 2,000 sweeps, then 1,000-sweep blocks while M is below the cap, until three flat blocks, to 20,000 at most. Each part is timed by M1's method: the median of nine batch means after a warm-up batch, round-robin with the sweep's batches. The sweep is a step plus a composition, so the figures are per sweep at Thin = 1.

*µs per sweep along each axis* (the record pass; M in brackets where it fell short of the cap):

| axis | values |
| --- | --- |
| n (rich, M_max 20) | 200: 7 (M 7); 500: 15; 1,000: 22; 2,000: 23; 5,000: 57; 10,000: 217; 20,000: 284; 50,000: 1,066 |
| M_max at n = 1,000 | 5: 20; 10: 22; 20: 22; 50: 70; 100: 118; 200: 98 (M 94) |
| M_max at n = 10,000 | 5: 82; 10: 112; 20: 217; 50: 255; 100: 741; 200: 1,388 |
| p, rich | 2: 18; 5: 20; 10: 22; 20: 24; 50: 25; 100: 20 (M 11); 200: 17 (M 8) |
| p, S10 | 2: 18; 5: 20; 10: 22; 20: 23; 50: 18; 100: 17 (M 9); 200: 16 (M 4) |
| K | 1: 22; 2: 22; 3: 20; 4: 21; 5: 20 |
| L (n 4,096) | 2: 67; 4: 86; 16: 118; 64: 115; 256: 105; 1,024: 79 |
| n_x × m (direct) | 100×10: 32; 100×50: 35; 100×200: 95; 100×1,000: 122; 500×10: 47; 500×50: 69; 500×200: 164; 500×1,000: 877 |
| n × M_max (the cross's extra points) | 200×5: 7; 200×50 and ×100: 9 (M 7); 5,000×5: 48; ×50: 126; ×100: 248; 20,000×5: 177; ×50: 551; ×100: 1,381 |

*Slopes, fitted above the fixed overhead* (record pass, then the re-run):
- **n:** 21.7 and 24.7 µs per 1,000 cases at M = 20, fitted on n = 2,000–50,000. Below n ≈ 2,000 the sweep sits at its floor of 15–23 µs.
- **M:** 1.2 and 1.6 µs per term at n = 1,000 (realized M 20–100), and 6.9 and 8.6 at n = 10,000 (M 20–200).
- **Direct functional:** 1.7 and 2.6 µs per 1,000 rows of N = n_x·m, fitted on N = 20,000–500,000.
- **p and K are flat** within the spread: 17–25 µs and 20–22 µs. Past p = 50 the time is read at a lower M (below).
- **L** reads 67–118 µs at n = 4,096 and is not monotone. Its membership test rises with L, but the changes' mix moves too.
- **The fixed part at the cap** is about 5–8 µs per sweep: τ 1.5–2 µs, a composition about 2 µs, a death 1–3 µs and the target 0.3 µs.

*The parts' shares of a sweep* (record pass). "Re-run ÷ record" is the grid's second full pass over the first.

| point | µs per sweep (IQR) | re-run ÷ record | change | death | birth | target | commit | τ | composition | rest |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| base (n 1,000, p 10, M 20, K 2) | 22.1 (5.9) | 1.08 | 58% | 7% | 0% | 1% | 4% | 9% | 9% | 13% |
| n 200 (M 7) | 7.0 (1.1) | 1.72 | 18% | 5% | 16% | 3% | 3% | 22% | 13% | 19% |
| n 50,000 | 1,065.9 (254.1) | 1.12 | 84% | 0% | 0% | 0% | 2% | 0% | 0% | 14% |
| M 5, n 1,000 | 20.1 (2.3) | 1.03 | 57% | 3% | 0% | 2% | 1% | 15% | 7% | 16% |
| M 94 of 200, n 1,000 | 98.4 (13.0) | 1.85 | 38% | 11% | 25% | 0% | 4% | 2% | 14% | 5% |
| M 5, n 10,000 | 82.4 (13.1) | 0.86 | 79% | 1% | 0% | 0% | 1% | 2% | 1% | 16% |
| M 200, n 10,000 | 1,388.4 (280.8) | 1.17 | 79% | 5% | 5% | 0% | 10% | 0% | 3% | −2% |
| p 2 | 17.5 (5.8) | 1.44 | 61% | 8% | 0% | 2% | 2% | 11% | 13% | 3% |
| p 200 (M 8) | 16.8 (8.5) | 1.29 | 28% | 2% | 28% | 1% | 4% | 9% | 5% | 22% |
| K 1 | 22.2 (4.5) | 1.64 | 54% | 7% | 0% | 1% | 5% | 9% | 10% | 15% |
| K 5 | 19.7 (4.0) | 1.32 | 57% | 8% | 0% | 2% | 3% | 9% | 11% | 11% |
| L 2 | 67.4 (8.6) | 1.62 | 69% | 3% | 0% | 0% | 12% | 3% | 3% | 9% |
| L 1,024 | 78.8 (18.5) | 1.30 | 75% | 2% | 0% | 0% | 4% | 2% | 3% | 14% |
| direct 100 × 10 | 31.7 (9.6) | 1.38 | 56% | 9% | 0% | 1% | 6% | 10% | 11% | 8% |
| direct 500 × 1,000 | 877.2 (94.0) | 1.51 | 53% | 0% | 0% | 0% | 56% | 0% | 0% | −9% |

*What the shares show:*
- **At the cap, the change proposal is the sweep.** The chain proposes no births there: half its moves are deaths and half changes. The change is 54–61% of a sweep at n ≈ 1,000 and 79–96% from n = 10,000, plus its commit when accepted. A change is a column evaluation plus that column's cross products: at the base, the column costs 10.9 µs, the cross products 9.3 µs, and the whole change 25.5 µs per call.
- **The composition is no longer the sweep's largest part.** M1 read it at 34–94% at Thin = 1. That was before `9ed600a` (2026-10-29), which caches each term's column mean. A composition now costs about 2 µs a draw at M = 20 and 44 µs at M = 200. With τ, the composition is a fixed 18% of the base sweep, 36% at n = 200, and under 3% from n = 10,000.
- **In direct functional fits, the commit takes 30–56% of a sweep from N = 20,000.** Refit rebuilds the residual explicitly, at O(N·M) per accepted move: 19.7 ms per commit at 500 × 1,000.
- **Below the cap, births take 16–28% of a sweep:** at n = 200, at p ≥ 100, and at M_max ≥ 100 with n = 1,000.
- **The rest stays within the sweep's IQR at all 59 points,** −18% to +23% of the sweep. The parts, weighed by their frequencies, account for the sweep. *(The main session's recount, 2026-10-31: the rest exceeds the IQR at 8 points, all positive and at sweeps of 7–70 µs, up to 25% at n = 200, M_max = 50; the fitting checklist's report.)*

#### Memory

- **Allocation per sweep** (the thread's counter; the process's agrees with it to a median of 1.1%, at most 5.8%): 9.6 KB at the base. It grows as the change's column, 8n bytes per change: 5.4 KB at n = 200 and 201 KB at n = 50,000. Elsewhere: 43–73 KB at n = 10,000 over M = 5–200, 15–132 KB in the direct fits, 21–24 KB at the L points.
- **By part,** in bytes per call at the base:
  - a column 8,152 (its n doubles), a change about 10.4 KB, a death 1.5 KB;
  - τ's step 1.2 KB, a composition 2.1–2.2 KB;
  - a commit 350–370 B, a draw 163 B, a ratio 136 B, a target evaluation 112 B.
  - Nothing: the bordered factor, the rebuild, the scalar cross products, and the design's refusal of a birth at the cap.
  - Direct functional: `KhatriRao.CentredBorder`'s re-centred copies, 17 KB per call at 100 × 10 and 219 KB at 500 × 1,000.
- **GC:** pauses take at most 6.5% of a sweep anywhere: 0–1.4 µs at M = 20 for n ≤ 2,000, and 24 µs at most (n = 10,000, M = 100). Gen-2 collections run at most 4.4 per 1,000 sweeps.
- **Retained per draw:** about 120 B plus 13 B per term (a coefficient and a term id). That is 193 B at M = 5, 365 B at M = 20 and 740–757 B at M ≈ 48, whatever n and m. The per-fit intercept is within ±65 KB of zero everywhere except the base point, 0.1 MB of one-time state left after the warm-up. A fit keeps no O(n) state.
- **Peak working set**, each read in its own process (two chains, a 10,000-sweep burn-in, 1,000 draws, one thread), over a host baseline of 85–86 MiB:

  | largest point | peak |
  | --- | --- |
  | n = 50,000 | 171 MiB |
  | M = 200 at n = 10,000 | 171 MiB |
  | direct 500 × 1,000 | 135 MiB |
  | basis r = 30 | 123 MiB |
  | p = 200 | 110 MiB |
  | L = 1,024 | 109 MiB |
  | K = 5 | 103 MiB |

#### Factor work and scans

*The machine's reference rate* (`machine`): the production `PairwiseSum.Dot` runs at 0.28–0.46 ns per multiply-add over n = 1,000–50,000. That equals a `Vector<double>` loop (0.30–0.42) and is 2–3× faster than a plain scalar loop (0.8–1.0). Memory streams at 10.4 GB/s on one thread, 0.77 ns per double.

*µs per call (ns per counted unit)*:

| point | M | column, n(ō + 4) | cross products, n(M + 3) | bordered proposal, 1.5M² | commit, NM | rebuild, M³/3 | death, 2M² | change | composition, M²/2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| base | 20 | 10.9 (2.1) | 9.3 (0.41) | 1.04 (1.7) | 24.3 (1.2) | 6.1 (2.3) | 3.1 (3.8) | 25.5 (0.86) | 2.0 (9.8) |
| n 200 | 7 | 2.0 (1.9) | 0.53 (0.27) | 0.22 (3.0) | 2.0 (1.4) | 0.67 (5.8) | 0.94 (9.5) | 3.9 (1.2) | 0.92 (38) |
| n 50,000 | 20 | 676 (2.5) | 968 (0.84) | 1.08 (1.8) | 1,410 (1.4) | 6.0 (2.2) | 3.6 (4.5) | 1,820 (1.3) | 2.2 (11) |
| M 5, n 10,000 | 5 | 116 (2.1) | 37 (0.46) | 0.25 (6.6) | 51.5 (1.0) | 0.40 (9.7) | 0.83 (16) | 132 (0.97) | 0.86 (69) |
| M 200, n 10,000 | 200 | 115 (2.1) | 1,900 (0.94) | 74.6 (1.2) | 2,870 (1.4) | 3,510 (1.3) | 137 (1.7) | 2,350 (1.0) | 44 (2.2) |
| M 100, n 1,000 | 100 | 12.8 (2.4) | 38.8 (0.38) | 20.6 (1.4) | 145 (1.4) | 470 (1.4) | 41 (2.1) | 132 (0.89) | 17 (3.4) |
| L 1,024 | 20 | 178 (7.5) | 37.1 (0.39) | 0.95 (1.6) | 91 (1.1) | 5.6 (2.1) | 3.1 (3.9) | 116 (0.97) | 2.2 (11) |
| direct 100 × 10 | 20 | 1.8 (5.1) | 16.7 (1.9) | 1.32 (2.2) | 62 (3.1) | 9.6 (3.6) | 5.4 (6.7) | 35.6 (3.3) | 3.5 (17) |
| direct 500 × 1,000 | 20 | 12.7 (3.0) | 891 (0.82) | 1.05 (1.7) | 19,700 (2.0) | 7.6 (2.8) | 5.0 (6.2) | 944 (0.86) | 2.8 (14) |

*The counted units,* as the code's loops run:
- The direct cross products are 3M(n_x + m) + 2(N + n_x).
- A change counts its column, its cross products and 4M².
- A death's copies are left out of its 2M².

*Frequencies per sweep:* changes 0.33–0.51, accepted moves 0.013–0.12, rebuilds 0.0001–0.0012 (one per 100 accepted updates), τ evaluations 5.7–6.1 per step. The columns evaluated per sweep equal the births with a term plus the changes, which are in the CSV.

*What runs at the machine's rate, and what is overhead:*
- **Cross products** run at the machine's vector rate in cache, 0.27–0.46 ns per multiply-add. Where n·M·8 bytes leave the caches (n = 50,000, or n ≥ 10,000 with M ≥ 100), they run at the DRAM rate, 0.8–0.94 ns.
- **Columns** take 1.8–2.5 ns per element pass, and 5–7.5 ns at the categorical points: 3–7× the machine's rate. They are overhead-bound per element: each kernel makes its own pass with its own hinge arithmetic.
- **Commits** (1.0–1.4 ns per N·M, scalar; 2–3 ns direct) and **factor work** (bordered 1.0–1.8, rebuild 1.1–1.4, death 1.7–4.5 ns per unit from M = 20) run near the scalar rate.
- **The composition and the death are overhead-bound below M ≈ 50,** at a fixed 1–3 µs.
- **Categorical membership costs 9.4 ns per case at L = 2 and 65 ns at L = 1,024** (16.5, 19.8, 24.1 and 37.5 between). It is a binary search of the sorted level list, read here for a subset of L/2 levels, against about 2 ns per case for a hinge factor.
- **The timer's own overhead** (the empty part) is 4–7 ns per call, reported and subtracted nowhere.

#### The outputs, the other backends and calibration

*The outputs,* per draw, from two chains of 250 draws after a 2,000-sweep burn-in. At M_max ≥ 100 these fits held M ≈ 46–49.
- **Prediction:** 0.06–0.13 µs per case per draw at p ≤ 100, and 0.24 at p = 200, flat in M and K (latent and predictive draws together). Direct fits: 0.06 µs per case per location per draw. Basis fits: 0.02 (r = 1) to 0.12 (r = 30).
- **Sobol per draw:** 2.9 µs at M = 5, 7.7 at M = 20 and 34–39 at M ≈ 46–49. Along p: 7.8–38; along K: 7–10; direct fits: 8–54.
- **Basis Sobol (`BassBasisSensitivity.Of`) is superlinear in r:** 20 µs per draw at r = 1, 156 at r = 3, 766 at r = 10 and 1,459 at r = 30.
- **Effects:** 4.5–9.6 ns per cell per draw at the cap points, and up to 27 where cells are few. At R's 100-point grid the base's 61,000 cells would cost about 0.35 ms per draw.
- **A basis component's fit:** 10.7–14.6 µs per component sweep at n = 500, flat in r.

*The other backends,* at the base point. The untempered step (moves and τ) costs 12.6 µs there.
- **Tempering,** with rungs geometric from b = 1 to the self-sized hot end b = 0.001:
  - a rung-sweep costs 13.8–17.4 µs;
  - a swap round costs 8.7, 44 and 69 µs at 5, 10 and 20 rungs: 11%, 24% and 17% of a scan;
  - a hand-over costs 4.9–8.6 µs, two per accepted swap at most.
- **SMC for fitting,** one thread, D05's defaults:
  - 1,000 particles took 205 stages (19 resampled), 27.6 s; 4,000 took 224 stages (21 resampled), 145 s;
  - a particle's mutation step costs 24.7–26.4 µs at the run's end, twice the untempered step;
  - a reweighting costs 18–83 µs;
  - a resampling costs 17.8 and 73 ms: the ancestors' draw is 28–99 µs, and about half the particles are forked at 36–37 µs each, so a fork costs more than a mutation step;
  - the stage's β search costs 5.2 and 27 ms;
  - the mutation is 97% of a stage's advance.
- **The flat histogram,** four chains on one thread, a flat start over a window about ℓ: an iteration costs 104 µs, 26 µs per chain step. Its own work is a weight lookup (0.02 µs), the update (0.48 µs) and a reweigh (0.02 µs): under 1% of an iteration. Its chains sat near M ≈ 8 in the window, a different state from the base's.

*Calibration (N15)*, `DiagonalPlusLowRank` per σ² proposal (construction and solve), in µs. The record pass ran under the run's heaviest background load (4.2 cores busy); its first pass read 1.8–2.3× faster, given in brackets.

| m \ k | 10 | 39 | 100 |
| --- | --- | --- | --- |
| 20 | 83 (40) | 421 (187) | 445 (240) |
| 50 | 97 (43) | 2,458 (1,031) | 4,601 (2,073) |
| 200 | 214 (86) | 3,264 (1,350) | 28,130 (14,160) |

- **The construction is 97–99.8% of a proposal:** a solve costs 0.8–57 µs.
- **The alternative is priced, not built.** One-sided Jacobi costs k(k − 1)/2·(7m + 4k) multiply-adds per sweep, 2.2·10⁵ at (20, 39) and 8.9·10⁶ at (200, 100), and the timings imply 2–8 sweeps at 0.4–1 ns each. N15's capacitance Cholesky costs m·k(k + 1)/2 + k³/6 + 2mk + k² per proposal: 2.9·10⁴ at (20, 39) and 1.2·10⁶ at (200, 100), 7.7× and 7.3× below a single Jacobi sweep. Where k > m (20 × 100, 50 × 100), Σ's own m × m Cholesky is cheaper still: m³/6 + m²k ≈ 4·10⁴ at (20, 100).
- **Stage 7's single reading** was 1.4 ms per σ² proposal at m = 20 with 39 columns; here the same size costs 0.19–0.42 ms. Stage 7's figure was taken inside a calibration run.
- **The cut SMC** at stage 7's pooled-posterior fixture (default spec, 1,000 particles, three legs, 50 stages, one thread): 3.1 s per leg and 185 ms per stage (first pass: 2.0 s and 123 ms). The fixture has m = 8, so N15's share of a calibration run at larger m is not read here.

#### The efficiency pass's levers, ranked

Each lever is priced as the share of a sweep it could remove, from the record pass:

| lever | base | n 200 | n 50,000 | M 5, n 10k | M 200, n 10k | p 200 | K 5 | L 1,024 | direct 100×10 | direct 500×1,000 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1. Prefix sums for a main-effect hinge: the cross products | 16% | 3% | 26% | 10% | 40% | 10% | ≤ 11% | — | — | — |
| … and the proposal's column, which the same sums would skip | 35% | 14% | 45% | 42% | 43% | 37% | ≤ 33% | — | — | — |
| 2. The column kernels (the column's share; it runs at 3–7× the machine's rate) | 25% | 19% | 31% | 69% | 4% | 46% | 31% | ≤ 100% | 3% | 1% |
| 3. The commit's explicit residual (Refit, O(N·M)) | 4% | 3% | 2% | 1% | 10% | 4% | 3% | 4% | 6% | 56% |
| 4. The composition and τ, fixed per sweep | 18% | 36% | 0% | 3% | 3% | 14% | 20% | 5% | 21% | 1% |
| 5. The state's allocation (GC pause; bytes per sweep) | 0% (9.6 KB) | 0% (5.4 KB) | 0.5% (201 KB) | 4.7% (43 KB) | 0.5% (73 KB) | 4.2% (9.9 KB) | 5.8% (9.4 KB) | 3.6% (24 KB) | 4.9% (15 KB) | 1.5% (132 KB) |
| 6. `KhatriRao.CentredBorder`'s re-centring (and KB allocated per call) | — | — | — | — | — | — | — | — | 13% (17) | 2.7% (219) |
| 7. The rebuild cadence | 0.01% | 0.01% | 0% | 0% | 0.12% | 0% | 0.01% | 0% | 0.01% | 0% |

*How each lever is priced:*
- **1. Prefix sums.** The share is the one-way share of proposals (2 − the mean proposed order: exact at K ≤ 2, at most (K − ō)/(K − 1) past it) times the births' and changes' cross products. The draft prices the cross products alone. The same suffix sums would also skip the proposal's column, which is needed only on acceptance. The categorical and direct points hold no one-way hinge proposals to price. This lever cuts part of a change too, which M1 named.
- **2. The column kernels** are the column's whole share. The lever removes at most the part above the machine's rate, about 60–85% of it. Its share at L = 1,024 is an upper bound: the birth pool's categorical columns cost more than the changes' (115% of a sweep read that way).
- **3. The commit's explicit residual** is the one the list misses. It is 30–56% of a sweep in direct fits from N = 20,000, and 10–12% at M = 200 and at L = 2.
- **4. The composition and τ** are not on the list. They are the largest share where the sweep is small: 36% at n = 200, 18% at the base.
- **5. The state's allocation** is priced at its measured GC pause, which bounds it: at most 6.5% anywhere. The allocation's own cost is not separable from the parts that allocate.
- **6. `CentredBorder`'s re-centring** is priced as the share of the direct cross products that re-centres every current factor, 2M(n_x + m) of the units. It matters only in small direct fits.
- **7. The rebuild cadence** is negligible: at most 0.12% of a sweep.
- **N15** sits outside the sweep. It lives in calibration's σ² proposals, where the construction is the whole cost (above). The capacitance Cholesky's count is 7–8× below a single Jacobi sweep, and the timings imply 2–8 sweeps per proposal.
- **The categorical membership test** (9–65 ns per case against a hinge's 2) sits inside lever 2 at the L points.
- **SMC's per-particle cost** is twice an untempered step, and a fork costs more than a mutation step. That is the SMC kernel's own overhead (a chain built per move, its account under a lock), outside the sweep.
- **A vectorized dot product is worth nothing:** the production dot already runs at the vector rate.

The ranking by share at the base: lever 1 with its column skip (35%), lever 2 (25%), lever 4 (18%), lever 1's cross products alone (16%), lever 3 (4%), then levers 5–7. Levers 1 and 2 overlap on the column, so their shares do not add. At large n and M, lever 1 rises to 40–45%. In direct fits, lever 3 is first at 56%.

#### The parallel and R context

- **Parallel fits** (the default four chains at the default lengths; `MaxDegreeOfParallelism` 1 against unbounded): speedups of 3.4 at the base, 2.0 at n = 10,000 and 3.0 on the direct 100 × 200. Efficiency over four cores is 0.51–0.84; the first pass read 2.5–2.7. The draws were identical across the two. The CPU is a hybrid with 4 P-cores, so four chains need E-cores or SMT siblings.
- **R BASS 1.3.1,** run from the pinned source under R 4.6.1 with R's reference BLAS ("Matrix products: default", LAPACK 3.12.1), over its default lengths (10,000 iterations, 9,000 burn-in):
  - Friedman at n = 500, p = 10, R's defaults: 713 µs per sweep (S10 read about 630), mean M 11;
  - the rich surface at n = 5,000, p = 10, maxBasis 20, maxInt 2: 764 µs, M 20;
  - the direct functional at n_x = 100, m = 50: 392 µs, M 20.
- **The C# grid reads 57 and 35 µs per sweep** at the last two points, at the burned state. R's figure covers its whole run, the burn-in at growing M included, with its own per-iteration updates and fixed ⅓ move probabilities. It is context, not a comparison.
- **`truncdist` is not installed.** In the sourced environment its truncated gamma draw is replaced by the same law's inverse CDF, as the fitting checklist's R runs did.

#### The hub's list

| hub measure | covered by | left |
| --- | --- | --- |
| case count | n 200–50,000; n × M_max | — |
| predictor count | p 2–200 on both surfaces | past p = 50 at a lower M (below) |
| active terms | M_max 5–200 at n = 1,000 and 10,000; realized M recorded | M 94 of 200 at n = 1,000 |
| interaction order | K 1–5 | the chain holds mostly one-way terms (mean held order 1.05 at the base) |
| response size and rank | direct n_x × m; basis r 1–30 (the fit, prediction, Sobol) | — |
| category cardinality | L 2–1,024 with membership per case | — |
| allocation and peak memory | bytes per sweep and per part, GC, retained per fit and per draw, peak per process | — |
| factor work | bordered proposal, commit, rebuild, death, change per counted unit | — |
| basis and cross-product scans | column and cross products per unit, against the machine's rate | — |
| sensitivity moment work | inside Sobol and effects per draw | the moment work alone (`BassMoments.On`) is not isolated |
| statistical error per unit time | — | the statistical half |

#### The machine, the trims, the flagged blocks and the deviations

- **The machine:**
  - commit `b484fca`, with `abe8a7c`/`ca729db` for `machine`; Release; .NET 10.0.5; Windows 11 build 26220;
  - Intel Core Ultra 5 125H (Family 6 Model 170 Stepping 4): 14 cores (4 P, 8 E, 2 LP-E), 18 logical processors;
  - workstation GC, Interactive latency, concurrent by default; tiered compilation and dynamic PGO by default (`DOTNET_TieredCompilation` unset);
  - power plan Balanced (`381b4222-f694-41f0-9685-ff5bb260df2e`).
- **No trims.** The trial ran in 4 min 55 s and projected about 50 minutes. The timed machine time actually used was about 51 minutes, every timed pass and the trial included. Waits were short: no other test host ran during any timed block.
- **No block was flagged:** every block saw no other test host at its start and end. The machine was never idle, though:
  - other processes kept 1.6–4.2 cores busy in the record pass and up to 6.9 in the grid's re-run (a long-running `python`, about one core throughout, plus IDE and desktop apps);
  - the two full grid passes, about 25 minutes apart, differ by a median 1.29× in the sweep (quartiles 0.98–1.51, range 0.70–2.60);
  - the change's share holds to a median of 4.3 points (90th percentile 15);
  - read absolute times to about ±30% and shares to about ±5 points. Thread placement on the hybrid CPU is not controlled.
- **Deviations,** each a declaration changed after a run, with its reason:
  1. **The rich surface:** amplitudes taper as 10/(v + 1) past the tenth predictor, unchanged for p ≤ 10. With them flat, p = 100 held 5 terms of 20 and p = 200 none (first counts pass, before any timed run).
  2. **The burn-in:** it now stops after three flat blocks rather than one (second counts pass, before any timed run).
  3. **The trial's projection:** the trial ran each point's burn-in whole and a tenth of its timed sweeps, because the burn-in sets the state the times are read at.
  4. **Sweep batches:** they fill a restored chain's column-mean cache before timing, and read GC (after the trial, before the grid).
  5. **Effects:** read on a 20-point grid with a 5·10⁷ cell budget. R's 100-point default passes its budget at the full draw count; the cost per cell carries to any grid.
  6. **Retained memory:** read after a warm-up fit.
  7. **The record:** every timed block now runs itself at a tenth first, its rows discarded. The first full pass's first readings in each block were partly tier-0 code: n = 500 read slower than n = 1,000, and the outputs' first point read 3–4× high. That pass is superseded. The record is the second, and the third (the grid only) is the re-run.
  8. **Smaller ones:** the change's analytic count was corrected (above); `BASS_COST_BLOCK` splits on semicolons, because axis names hold commas; the peak fits burn 10,000 sweeps.

#### The checks, the review passes, and what the main session should know

*Checks:*
- **The timer:** the empty part costs 4–7 ns per call, reported and subtracted nowhere.
- **The allocation counter:** a `byte[1000]` reads 1,024 B. `CollapsedGPrior.MarginalAt` × 1,000 reads 0 B. Both are asserted.
- **Nothing commits:** `BassMixing.Untouched` passes at all 59 points, in both passes.
- **The parts against the sweep:** the rest is never below −IQR (−18% to +23% of a sweep). The test-side replicas match the design's own proposal: design − column − cross − Cholesky is within ±3.3 µs, and 9 of 73 µs at n = 1,000, M_max = 200.
- **Repeatability:** the base point read twice in one process agrees within its IQR on 16 of 18 parts; the rebuild and the death's commit differ by 0.5 and 2.5 µs against IQRs of 0.25 and 2.1. Its counts agree exactly (asserted).
- **The surfaces:**
  - the rich surface reaches the cap at every M point but n = 1,000, M_max = 200 (94);
  - it misses at n = 200 (7 terms for caps ≥ 20) and at p = 100 and 200 (11 and 8);
  - S10's surface sets its own M (4–9 past p = 50);
  - L's designs give every level at least 4 cases (asserted, at n = 4,096).
- **Draws unchanged:** the 8 `BassBitIdentityTests` hashes pass with the harness in the tree, and the four new facts return at once with the switch off.

*Review passes and what they fixed:* the effects budget; the retained-memory warm-up; the tier-0 first readings; the change's analytic count, against `ProposeChange`'s loops; the machine's reference rate, found missing; and that rate's own tier-0 reading of the production dot (3 ns per multiply-add before its warm-up). No flagged block is reported as quiet.

*Coverage that remains:*
- **M missed its cap** on the axes named above.
- **Parts read at one point only:** the backends at the base (as briefed); the SMC fork's cost, which grows with the design it copies; the flat histogram's chain-step cost, read at its chains' own states. N15's share inside a calibration run at m > 8 is not read.

*For the main session:*
- **The cost model at the cap:** sweep ≈ 5–8 µs + ½ × change + commits. A change costs a column, n(ō + 4) element passes at 2–2.5 ns, plus cross products, n(M + 3) multiply-adds at 0.3 ns in cache and 0.9 ns beyond it. That is what the statistical half multiplies by its ESS.
- **M1's composition-dominance reading is superseded** by `9ed600a`'s column-mean cache. Lever 1's price now rests on the change, not the composition.
- **No default is recommended,** and nothing was built.
