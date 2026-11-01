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
