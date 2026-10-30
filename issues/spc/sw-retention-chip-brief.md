# SW retention, the SPC-side items (chip brief)

Date: 2026-10-29. Builds the items of the [SW retention brief](sw-retention-brief.md) that stay inside SPC: S1 with S3, S2's Binder view, S4, S6 and S7, with a Binder-crossing check. Nothing here reaches into `src/maths` (the samplers, regression) or the BASS work.

- *The frame* (the user, 2026-10-29): accumulation is the run's state, not a resume mechanism. The per-sweep record is the state; the running sums are its views. Resume was one use of that state, and SPC does not use it: the SPCX file caches finished tasks, and the engine's `Restore` has no caller.
- *Decided* (the user, 2026-10-29, the retention brief's SQ1 and SQ2):
  - the per-sweep record is always collected, as the running sums and the histogram are, with no switch and no thinning;
  - it lives on `Accumulator`, so a task loaded from its SPCX file carries it;
  - S4 (the burn-in trace) and S7 (the exact histogram) are `AccumulationSpec` options, off by default.
- *What exists.*
  - The cycle, `SwendsenWang.Draw` (`src/clustering/graphical/spc/runtime/core/sampler/`): Pass 1 walks the edges, forms the bonds and sums `cycleEnergy` over the current spins; Pass 2 recolours each bond cluster; Pass 3a–3e accumulate Σ|c|² (in `long`), the largest cluster `maxSize`, the Blatt magnetization `m`, the log-binned histogram (64-bit since `1c42741`), the landscapes and co-membership. `BurnIn` runs `Draw` and then `ResetAccumulators`.
  - `Accumulator`, `ISwEngine`, `GetCheckpoint` and `Restore`; `AccumulatorSerializer` (`SwRunSerialization.cs`, SPCX v8): fixed fields, then each optional array behind a presence byte.
  - The Potts reductions in `src/graphs/models/potts/observables/` (`SpecificHeat.Reduce` and its siblings), assembled over temperature by `SweepProfile.From` and `SweepCurves.ByTemperature` (`src/clustering/graphical/spc/profiling/`). `ClusterSizeEntropy` (`src/graphs/observables/`) reduces a size histogram.
- *Not in scope,* because each reaches into the samplers or waits on a decision:
  - R-hat, ESS and τ_int (`ChainDiagnostics`), the convergence verdict, temperature integration (`LadderEvidence`), Ferrenberg–Swendsen reweighting, and C_v's error bar, which needs τ_int;
  - S5, the per-sweep partitions (persistence and SQ5);
  - wiring `Restore` so a cached task can be extended;
  - exporting the series (SQ3).

Read first: [AGENTS.md](../../AGENTS.md) and [CONTRIBUTING.md](../../CONTRIBUTING.md); the retention brief in full; [SpcEngine.md](../../docs/Clustering/SpcEngine.md), its SPCX layout and configuration table.

## Working rules (the shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-sw-retention\` to every `dotnet build` and `dotnet test`. Run the VizCore tests in Release. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-sw-retention\` when you finish.
- **Other sessions are running** on the BASS work in `src/maths` and `issues/regression`. Expect their commits between yours, and leave their files alone, tracked or untracked.
- **Your files.**
  - `src/clustering/graphical/spc/runtime/core/sampler/` and `src/clustering/graphical/spc/profiling/`.
  - New reductions in `src/graphs/models/potts/observables/` (the Binder cumulant) and `src/graphs/observables/` (S7's cluster functionals).
  - `docs/Clustering/SpcEngine.md`: the SPCX layout, the options, and what each recorded value describes.
  - Tests in `tests/VizCore.Tests/`.
  - This brief's Reports section.
  - Leave other production files alone. If you need a change in one, make the smallest behaviour-free change, name it in your report, and rerun the suites that cover it.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Search with `rg` through the shell;** the repo's hooks block the Grep and Glob tools.
- **Line endings:** briefs and C# are LF; `docs/` files keep theirs (`SpcEngine.md` is CRLF). Check `git diff --stat` after each doc edit.
- **Notation:** use the symbols SPC already uses (T, E, m, Σ|c|², c_max). A split you notice goes in [the notation ledger](../notation-ledger.md); rename nothing.
- **Decisions are not yours to settle.** Where this brief leaves something open, implement a reversible recommendation and put the question in your report with your reasons.

## What to build

Tasks in order, one commit each.

1. **The probe, before any change.** At the commit you start from, record the accumulator's outputs for fixed-seed runs: the running sums, the histogram's counts, the spins, the RNG words, the per-edge counts and the landscapes. Use a triangle and a small random graph, three temperatures (one cold, one near the transition, one hot), and `AccumulationSpec.None` and `Currencies` with both landscapes. Keep them as test expectations. Every later task leaves them bit-identical: retention must not change the RNG stream or any existing output.
2. **S1 with S3, the per-sweep record.**
   - One record per measured sweep: E (`cycleEnergy`), m, Σ|c|² and c_max. Σ|c|² without the largest cluster is Σ|c|² − c_max², a view, not a fifth value.
   - Plain arrays on `Accumulator` and the `ISwEngine` read surface, with no `Maths.Samplers` type. A later adapter to `ChainTrace` is a view at that seam. The storage shape and names are yours; name your choice in the report.
   - The record resets with the running sums in `ResetAccumulators`, so the measured series starts after burn-in. `GetCheckpoint` copies it. `Restore` restores it and checks its length against `DrawCount`, so a future extension appends.
   - The engine grows the series; a capacity taken from the run's budget avoids copies.
   - **What each value describes.** E, Σ|c|² and c_max belong to the spins entering the sweep (and its bonds); m belongs to the spins leaving it, since Pass 3b reads them after Pass 2's recolouring. Record this in the docstrings and in `SpcEngine.md`. A view that pairs E with m for the same spins pairs `m[t]` with `E[t+1]`; no such view is built here.
   - SPCX v9: the record follows the fixed fields, always present. Update `SpcEngine.md`'s layout.
3. **S2, the Binder cumulant.**
   - A Potts reduction beside `SpecificHeat`: U = 1 − ⟨m⁴⟩/(3⟨m²⟩²) from the series.
   - Assembled per temperature in `SweepProfile` from the replicas' pooled moments, as its magnetization variance is, not by averaging per-replica U.
   - **The check:** the q = 2 Potts model on L × L periodic lattices with unit couplings, where the crossing sits at the exact T_c = 1/ln(1 + √2) ≈ 1.1346 (H = −Σ J δ, J = 1). No lattice builder exists; write one as a test helper from `CsrGraph.FromEdges`. Use L = 8, 16 and 32 and temperatures bracketing T_c. Fix the tolerance on each pairwise crossing before you run it, and report the crossings, the tolerance and the runtime.
4. **S4, the burn-in trace.** An `AccumulationSpec` option, off by default. When on, the burn-in sweeps' records are kept in their own arrays before `ResetAccumulators` clears the measured ones, and they are serialized behind a presence byte (SPCX v10).
5. **S6, the replica spread.**
   - `SweepCurves.ByTemperature`, or a sibling beside it, also returns each temperature's standard deviation across replicas and the replica count. The deviation is NaN for a single replica (the default is 1).
   - `SweepProfile`'s curves stay numerically unchanged.
6. **S7, the exact cluster-size histogram.**
   - An `AccumulationSpec` option, off by default: 64-bit counts at each size 1..N, filled in Pass 3c beside the log-binned bins, reset with them, serialized behind a presence byte (the same or the next version bump), and validated by `Restore` as the other optional arrays are.
   - Reductions in `Graphs.Observables`, beside `ClusterSizeEntropy`:
     - S₂ = ⟨Σ_k c_k²⟩/N²;
     - the size-floored count F_n = ⟨Σ_k (c_k − n)₊⟩/n at any floor n;
     - ⟨Σ_k c_k ln c_k⟩.

## Checks

- *The probe* (task 1) after every task.
- *The series against the sums:* each series, summed in order, equals its running sum bit for bit, on every probe run.
- *The exact histogram against what is already kept:* Σ_k k·n_k = N × draws; Σ_k k²·n_k equals Σ|c|²'s running sum; rebinned to [2^k, 2^(k+1)), it equals the log-binned histogram.
- *SPCX round trips* for every new array, present and absent, and a v8 file refused with its version in the message.
- *The Binder crossing* (task 3).
- *Burn-in:* the trace's length is the burn-in count, and the measured series is unchanged by turning the trace on.

## Review before reporting

Before you write the report, take a fresh-eyes pass and fix what it finds:
- *Is it correct?* Which spins each recorded value describes; the reset between burn-in and measurement; the series and the sums agreeing; the presence bytes and versions.
- *What is missing?* An item above with no test, an option without its documentation, a doc line the format change left stale.

Name the passes and their fixes in the report.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the tests and suite counts, including a full VizCore run at the end;
- any change outside your files;
- the questions you left as reversible defaults, with your reasons;
- the review passes;
- anything the main session must know, including stale lines in the retention brief.

## Reports

### 2026-10-29: all six tasks landed

**What landed.** The record's bits were never regenerated: every later commit leaves the 12 probe cases bit-identical, with the options off and on.

| Task | Commit | What |
|---|---|---|
| 1. The probe | `6de243a` | `SwAccumulatorProbeTests`: a triangle (q = 4; T = 0.3, 1, 3) and a 16-node random graph (q = 20; T = 0.1, 0.6, 2), under `None` and `Currencies` with both landscapes, 50 + 400 draws. Recorded at 6233f93: the six sums and four RNG words as bit patterns, plus a SHA-256 of the histogram, spins, per-edge counts and landscapes. The random graph's middle temperature sits at its χ_m and C_v peak (T ≈ 0.6–0.65). |
| 2. S1 with S3 | `bc70bf2` | The per-draw record on `Accumulator` and the `ISwEngine` read surface. Also: SPCX v9, the `Restore` length check, and `InternalsVisibleTo` for the tests. |
| 3. S2 | `f542951` | `BinderCumulant` (`Moments`, `Reduce`) in the Potts observables, and the `BinderCumulant` channel in `SweepProfile`, reduced once from moments pooled across replicas. Includes the lattice crossing check. |
| 4. S4 | `1322e47` | `AccumulationSpec.BurnInTrace` and SPCX v10. |
| 5. S6 | `96de69d` | `SweepCurves.ByTemperatureWithSpread` → `SweepCurve(Temperatures, Means, StandardDeviations, Replicas)`. |
| 6. S7 | `d7269ea` | `AccumulationSpec.ExactClusterSizeHistogram`, SPCX v11, and `Graphs.Observables.ClusterSizeSums` (`Mean` for any g, `SecondMoment` = S₂, `FlooredCount` = F_n, `SizeLogSize`). |
| Review | `8c4c1ae`, `564a1ab` | The fixes below, and one notation-ledger entry. |

**The Binder crossing.** The tolerance was fixed before the first run, from a timing pass alone: |T× − T_c| ≤ 0.02 for each pair. The setup:
- T ∈ {1.09, 1.11, 1.13, 1.15, 1.17, 1.19};
- 4 replicas per point, 1,000 burn-in draws, and 25,000 / 12,500 / 6,250 measured draws for L = 8 / 16 / 32;
- seeds from `SeedTree` per L, then `SpcSeedHelper` per (T, replica);
- the crossing is the linear interpolation at the single sign change of U_larger − U_smaller.

The crossings sit 0.0022–0.0023 below the exact T_c = 1/ln(1 + √2) = 1.13459:

| Pair | Crossing |
|---|---|
| 8 × 16 | 1.1324 |
| 8 × 32 | 1.1323 |
| 16 × 32 | 1.1323 |

U near the crossing is about 0.616; the 2D Ising U* is 0.611. The check takes 8.1 s alone (`Parallel.For` over 72 runs) and about 4 s warm.

**Tests.** 93 new facts, all passing:

| Class | Facts | What they check |
|---|---|---|
| `SwAccumulatorProbeTests` | 24 | 12 cases, with the options off and on |
| `SwDrawRecordTests` | 25 | the series against the sums on every probe run; the burn-in reset; m on the leaving spins and E on the entering spins, bit for bit; v11 round trips, present and empty; v8, v9 and v10 refused with the version in the message; growth; `Restore` appends and refuses a short record |
| `BinderCumulantTests` | 7 | the reduction, pooling, NaN without a record, the crossing |
| `SwBurnInTraceTests` | 10 | the trace's length is the burn-in count; it is the head of the unburned stream; the measured record is unchanged; round trips; `Restore` |
| `SweepCurvesSpreadTests` | 3 | the means match `ByTemperature` bit for bit; the sample SD; empty input |
| `SwExactClusterSizeHistogramTests` | 24 | Σ s·n_s = N × draws, Σ s²·n_s = Σ\|c\|²'s sum, and the rebinning equals the log bins, on every probe run; the reset; round trips; `Restore`; the reductions by hand; S₂ = χ_FK/N on a run |

The full VizCore Release runs, both built into `artifacts\chip-sw-retention\`:
- after task 6: 357 passed, 2 skipped (the Möbius golden drafts), 0 failed;
- at the end, after a clean (`--no-incremental`) build: 357 passed, 2 skipped, 0 failed, in 53 s. The build raised no warnings in the changed files; three CS1574 cref warnings in the sampler folder (`IReduction`, `PottsModel{TConfig}`) predate this work.

**Changes outside my files.** All are behaviour-free, and the full VizCore run covers each.
- `projects/Clustering.Graphical.SPC/Clustering.Graphical.SPC.csproj`: `InternalsVisibleTo VizCore.Tests`, so tests can drive the internal engine's `Restore` and bare `Draw`.
- `docs/Clustering/SpcAnalysis.md`: the Binder bullet; "seven series" corrected to eight, with `BinderCumulant` added; the `ByTemperatureWithSpread` paragraph and its row in the surface table.
- `docs/Graphs/GraphSubstrate.md`: `ClusterSizeSums` in the observables row, with a short paragraph.
- `issues/notation-ledger.md`: "One SW update: draw, cycle and sweep".

**Reversible defaults, open for the user.**
- *The record's shape.* Four parallel plain arrays, named "PerDraw" after `DrawCount` and the existing `…PerNode` arrays: `EnergyPerDraw` and `MagPerDraw` (double), `SumSqClusterSizesPerDraw` (long, the cycle's own integer) and `MaxClusterSizePerDraw` (int). That is 28 bytes a draw, and every value is exact. The alternative is one record type of the four arrays, shared by the record and the trace.
- *Not `required`.* The four arrays default to empty, so the seven hand-built frames in the tests stay valid. A consumer reads "no record" as length ≠ `DrawCount`. The SPCX reader accepts any length; `Restore` enforces it.
- *Capacity.* `Run(n)` and `BurnIn(n)` reserve their own counts, so the budget arrives through the call and no constructor parameter is needed. Bare `Draw`s double the capacity.
- *The trace's scope.* It keeps every draw a `BurnIn` reset discards, including bare draws before it, and appends across repeated `BurnIn` calls. With the option on and no burn-in, it is present and empty.
- *Versions.* One bump per task (v9, v10, v11), since S4 and S7 landed in separate commits.
- *Binder's degenerate cases.* NaN when ⟨m²⟩ ≤ 0, or when a frame lacks the record, rather than the house 0: U = 0 would read as "disordered". In `SweepProfile`, ⟨m²⟩ comes from the running sums, which match the series' mean bit for bit.
- *Binder for q > 2.* m is Blatt's scalar order parameter, which is \|M\|/N for q = 2. For q > 2 this U is a scalar cumulant, not the vector-Potts one.
- *Not built.* The E³ and E⁴ views of S2, and a fourth-moment channel.
- *S6.* A sibling method with the sample SD (n − 1, so NaN for one replica). `ByTemperature` is untouched rather than delegating, so `SweepProfile` cannot move; delegating later would remove the duplicated grouping. The spreads are not `SweepProfile` channels.
- *S7.* Length N + 1, indexed by size, with entry 0 zero. F_n takes any real floor > 0. The reductions return 0 for no draws, the house style.
- *The crossing's cost.* The pre-registered budget is kept. The crossings agree to 10⁻⁴ with a tenfold margin, so a smaller budget would still pass if the suite's time matters.

**Review passes.**
1. *Correctness,* a fresh read of the engine, serializer and profiling diffs since 6233f93:
   - which spins each value describes: docstrings, `SpcEngine.md`, and the two alignment tests;
   - the reset: the record count, the exact histogram, and the trace appended before the reset;
   - series = sums: bit equality on all 12 probe runs;
   - the section order in writer and reader, and the version refusals.
   - *Found and fixed:* growing the record after burn-in copied the burn-in's stale capacity (`Array.Resize`), contradicting "a budgeted run never copies". Growth now copies only the live entries.
2. *Coverage and documentation:*
   - *Found and fixed:* `BinderCumulant` was missing from `SpcEngine.md`'s Potts surface and observables tables, and from its pooling and NaN notes.
   - *Found and fixed:* the memory costs lacked the record, trace and exact histogram, and the Accumulator concept line did not name the options.
   - *Found and fixed:* stale docstrings: `ISwEngine.GetCheckpoint`, `SwRunResult` and the serializer remark (review commit); `AccumulationSpec` and `SwRunSpec` (fixed with task 4).
   - *Found and fixed during the tasks:* `SpcAnalysis.md`'s "seven series" (task 3), and an unescaped `|` in a table cell (task 4).
   - Every check in this brief has a test; none was left without one.

**For the main session.**
- *Old caches throw.* SPCX is v11, and `SpcxDiskFrameSink.TryLoad` does not catch the version refusal, so a checkpoint directory holding v8–v10 files throws `InvalidDataException`. Delete old checkpoint directories. A cache also loads without its options' arrays (null) if it was written without them; the "Resume works per task" caveat covers this.
- *The CLI.* `SpcCommand.ParseAccumulationSpec` (user-repl) has no tokens for `BurnInTrace` or `ExactClusterSizeHistogram`. That is a small user-repl follow-up.
- *A new export column.* `SweepProfile`'s tabular export gains a `BinderCumulant` column.
- *Stale lines in the retention brief:*
  - "Where things stand": the accumulator list lacks the per-draw record and the two options, and the serialization is now SPCX v11, not v8.
  - S1's cost: the record is 28 bytes a draw (two doubles, a long, an int), not four doubles: 0.28 MB per run at 10⁴ draws, about 54 MB over 48 temperatures × 4 replicas.
  - S2: "if S1 ends up optional" is moot. Binder is built; E³ and E⁴ are views over `EnergyPerDraw`, not built.
  - S3, S4, S6 and S7 are built, so Sequencing steps 1–4 are done.
  - The Report section still says "Nothing yet".
