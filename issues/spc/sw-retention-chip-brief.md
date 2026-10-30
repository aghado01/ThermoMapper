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
