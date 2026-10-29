# BASS F-M2 Phase 0: the measurements (chip brief)

Date: 2026-10-28. Follows the F-M2 design draft, [local-move-mixing-brief.md](local-move-mixing-brief.md). The user decided on 2026-10-28 to run its Phase 0 as one chip; that is this brief. Phase 0 measures and changes nothing: its readings feed the user's design pass on the levers, and Phase 1 builds nothing until the user has read this chip's report.
- *What exists.* Stage 4 is on `main`, S10's qualification included (`tests/regression/BassQualification*Tests.cs`, behind `BASS_QUALIFICATION` = 1). Births and deaths accept about 1.5% of attempts on Friedman, and a tempered Friedman ladder holds one steep band at b ≈ 0.16–0.32. The design draft's "What S10 found" has the numbers.
- *This chip builds:* the harnesses for M1–M5, each a qualification test, and a report of their readings.
- *Not in scope.*
  - Any lever: no new kernel, schedule, option or default.
  - The open questions in the design draft, which are the user's.
  - Phase 1 or Phase 2.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the design draft in full, above all its frame (two factors) and Phase 0, which this brief makes concrete;
- the fitting brief's D01 (term representation, proposals, move schedule) and D04 (the sweep), and "Exploratory support in births";
- the fitting checklist's reports "Circle-back: argon's notes on S7–S10" and "S10: the qualification", for the runs these measurements continue;
- `docs/Maths/BayesianRegression.md`, sections "BASS: the run", "BASS: tempering and the evidence" and "BASS: qualification".

In the code:
- *The moves:* `BassTermProposal.cs` (`DrawBirth`, `DrawChange`, `LogBirthRatio`, `LogDeathRatio`) and `BassSweep.cs` (`BassStructureMoves` and its `Birth`, `Death` and `Change`; `BassChain`; `BassTauBlock`).
- *The design and target:* `BassDesign.cs` (`ProposeBirth`, `ProposeDeath`, `ProposeChange`, `Snapshot`) and `BassTarget.cs` (`Evaluate`, `LogStructurePrior`).
- *Tempering:* `BassLadder.cs`, `BassRun.cs`, `BassDiagnostics.cs` (`BassTemperingReadout`) and the shared `src/maths/samplers/rjmcmc/ParallelTempering.cs` and `TemperingLadder.cs`.
- *The tests to build on:*
  - `BassQualificationRecoveryTests.cs`: the fixtures, `Config(seed, tempering)`, and the per-pair profile;
  - `BassQualificationCostTests.cs`: the cost surface and grid;
  - `BassQualification.cs`: the switch and helpers;
  - `BassSweepTests.cs`: how a chain is built and stepped in a test.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-mixing\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-mixing\` when you finish.
- **Run the heavy tests filtered,** with `BASS_QUALIFICATION=1`, one at a time, and run timing measurements (M1) with nothing else heavy running on the machine.
- **Your files.**
  - New test files `tests/regression/BassMixing*Tests.cs`, plus a helper file of the same prefix if you need one.
  - This brief's Reports section.
  - Leave everything else alone, the design draft and the checklist included. If you need a change in one, name it in your report.
- **Production code stays as it is.** The tests can see the internals.
  - Where a measurement cannot be taken without a change, the only one allowed is an internal accessor that changes no behaviour: no new state, no change to any draw.
  - Name each such change in your report, and rerun the Regression, Samplers and Oracle (`--filter "FullyQualifiedName~Bass"`) suites to show the defaults unchanged.
- **The harness changes nothing it measures.** A proposal leaves its candidate pending in the chain's design, and the design commits it when the chain next shows that candidate's structure.
  - When you propose from a fixed state without stepping, check that the design and counts come back as they were. Compare the design's statistics before and after, bit for bit, and the next proposal's fit against a freshly loaded design's.
  - Assert it in the test.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per measurement or a coherent pair. Use the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF.
- **Read, do not judge.** Each harness prints its readings and asserts only its own soundness: the no-commit check, finite ratios, counts that add up. It asserts no pre-registered criterion, as `BassQualificationCostTests` does not.
- **Decisions are not yours.** Where a reading bears on a lever, say what it shows and how strongly. Recommending a lever's adoption is the pass's job (Phase 2), not this chip's.

## The measurements

The fixtures are S10's Friedman (seed 53), mixed categorical (57) and direct functional (61), with `BassQualificationRecoveryTests`' data and `Config`. Each fixed state is a chain's state at the end of a 5,000-sweep burn-in, four chains per fixture. Friedman's inputs 5–9 are inert, and its truth tags a term as inert when any of its continuous factors uses one.

**M1. Where a sweep's cost goes** (`BassMixingCostTests`).
- *The grid:* `BassQualificationCostTests`' surface, n ∈ {200, 1,000, 5,000} × M_max ∈ {5, 20, 50} × p ∈ {2, 10}. At each point, take a state after burn-in.
- *Time each part,* median of repeats after a warm-up, one thread, Release:
  - a birth's term draw;
  - the new column's evaluation;
  - the design's proposal, split if you can into the cross products Zᵀz and the bordered Cholesky;
  - the target's evaluation of the candidate;
  - the commit on an accept, including the explicit residual (time it on proposals forced to commit, then restore);
  - a death's and a change's proposal;
  - one τ block step, with its count of the slice's evaluations;
  - one composition draw;
  - the rest of a sweep: the block scan's bookkeeping and the state's allocation. Take this as a whole sweep's time less the parts, so the parts add up.
- *Tempered:* on Friedman tempered (S10's), one scan's swaps and hand-overs per pair, against its replicas' sweeps.
- *Report:* a table per part over the grid in µs, and each part's share of a sweep. Say what the shares imply for several structure moves per sweep (lever 1): which share a second move would add to and which it would not.

**M2. Which births fail, and why** (`BassMixingBirthTests`).
- *The harness.* From each fixed state, propose 10,000 births and score each with the target, committing none.
- *Record each proposal's outcome:*
  - an empty draw;
  - a refusal, by reason;
  - or the log Metropolis–Hastings ratio split into its parts:
    - Δℓ, the log marginal's change;
    - the count prior's ratio;
    - the term prior's log, with the within-move term −ln(M + 1) − ln q(t | S);
    - the selection ratio d_{M+1}/b_M.
  - Check that the parts sum to the ratio the chain would use.
- *Tag each proposed term* by its orders, and on Friedman whether it is inert. Record its knots' quantiles within their variables' distinct values, and its signs.
- *Read, per fixture:*
  - the share of births lost to each cause;
  - the acceptance probability min(1, e^{ratio}) by tag (order; inert or active);
  - Δℓ's quantiles for active and for inert terms;
  - the largest Δℓ among the proposals, against the threshold the other parts set.
- *Deaths:* from the same states, remove each term of the state in turn and record its Δℓ and its death's full ratio.
- *Report:* tables per fixture, then what bounds the gain of use-weighted main-effect variables (lever 3a): the active and inert shares of main-effect births and their acceptance.

**M3. Changes by knot displacement** (`BassMixingChangeTests`).
- *The harness.* From the same fixed states, propose 10,000 changes each and score them as in M2.
- *Group hinge changes:*
  - by sign flip or not;
  - by the knot's displacement |j′ − j| in the variable's distinct-value index, in buckets 0, 1, 2–3, 4–7, 8–15, 16–31, 32 and over.
- *Group categorical changes* by the subset's size change.
- *Read, per group:* its share of proposals, the refusal share, and the mean acceptance probability min(1, e^{Δℓ}).
- *Report:* the table per fixture, and what the small displacements' acceptance against the redraws' says about a local knot step (lever 2).

**M4. The chains' modes** (`BassMixingModeTests`).
- *The runs:* S10's untempered Friedman fit (seed 53), plus the mixed categorical and direct functional fits.
- *Per chain over the sampling phase:*
  - the unnormalized log posterior (the log marginal plus the log structure prior), and its quantiles;
  - M's band;
  - the structures it holds, by the variables and orders of their terms, with the share of draws in its three most frequent.
- *Across chains:*
  - the differences of their median log posterior;
  - whether any chain's states reach another chain's log-posterior range.
- *Report:* the table, and whether the chains hold modes of comparable posterior or one sits many nats below.

**M5. The steep band** (`BassMixingBandTests`).
- *Part 1, S10's tempered Friedman,* self-sized as in S10. Along each chain's frozen ladder, record each rung's mean and sd of M and of the log likelihood ℓ over the sampling phase.
  - The readout keeps neither per rung, so read them from the replicas: through `ParallelTempering`'s trace, or by stepping the ladder in the test. An internal accessor is allowed under the rule above.
  - Report the rows beside the per-pair rejections.
- *Part 2, a dense fixed ladder,* adaptation off, through the band. For example:
  - rungs from 1 down to 0.4, spaced as the self-sized ladder spaced them there;
  - 31 rungs from 0.40 to 0.10 at Δb = 0.01;
  - then the self-sized spacing down to the hot end 1/N.
  - The same fixture, seed and lengths; record the per-pair rejections and the same per-rung rows.
  - Report whether the band's pairs fall below 0.9 at Δb = 0.01.
- *The Gaussian reading's increment* h(r) = √(2/π)·Φ⁻¹((1 + r)/2), the design draft's lever 4a, for each band pair of both ladders:
  - floor 1 − r̂ at 1/(2·scans);
  - report the band's Σh against its Σr;
  - and, from Part 2, whether Σh read on Part 1's coarse pairs predicted Part 2's finer Σr.
- *Report:* whether the band is a finite barrier that the sizing under-reads, or a regime change whose rungs stay saturated at Δb = 0.01, with the per-rung M as the evidence.

## Review before reporting

Before you write the report, take a fresh-eyes pass over the harnesses and the readings, and fix what it finds:
- *Correctness:* the no-commit check holds; the ratio's parts sum to the chain's ratio; each table's counts add up; the timings sum to the sweep.
- *Coverage:* a fixture or bucket left empty; a reading that depends on one seed or one chain where the others disagree.

Name the passes and their fixes in the report.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- each measurement's tables and its reading, kept apart: the numbers, then what they show, and how strongly;
- the run times;
- any internal accessor added, and the suites rerun;
- the review passes;
- anything the main session must know.

## Reports

Newest first.

### Phase 0 measured: M1–M5 — 2026-10-28

*Landed* (tests only: no production change, no internal accessor):
- `fc88448`: the shared helper `BassMixing.cs` and M2, `BassMixingBirthTests`.
- `6b365ff`: M3, `BassMixingChangeTests`. `98381db`: M4, `BassMixingModeTests`. `11f1204`: M5, `BassMixingBandTests`. `a13eb80`: M1, `BassMixingCostTests`.
- Each runs behind `BASS_QUALIFICATION` = 1, reads, and asserts only its own soundness.
- *The helper.* S10's fixtures, a copy of `BassQualificationRecoveryTests`' generators, which are private there. Each chain or ladder is started, seeded and burned as `BassEnsemble` does it. `Score` splits a proposal's log ratio into its parts and checks them against the ratio `ReversibleJumpChain.Step` forms. `Untouched` is the no-commit check, bit for bit: the design's terms, statistics, β, residual and factor; the next birth, change and death against a design restored from the copy taken before, and against a fresh load to 1e-9; and the counts, through births drawn by the moves and by fresh counts from twin streams.
- *The fixed states are S10's.* M4 asserts that its chains' M equals the run's retained draws, draw for draw, on all three fixtures. M5's Part 1 reproduces S10's steep pairs exactly.

#### M1. Where a sweep's cost goes

*Per call, µs* (run c of two on the final code; one thread, Release; the median of nine batch means, each part's batches run round-robin with the sweep's). At M = M_max the design refuses a birth before its column and the chain proposes none, so there the birth's residuals are not read (—). τ's column gives its evaluations per step in brackets. The sweep is a step plus a composition.

| n | p | M_max | M | draw | ratio | column | Zᵀz | Cholesky | design other | candidate | death | change | target | commit b/d/c | rebuild | τ (evals) | compose | sweep |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 200 | 2 | 5 | 5 | 1.30 | 0.99 | 8.45 | 0.63 | 0.92 | — | — | 5.70 | 20.19 | 1.00 | 8.8/6.2/9.4 | 3.1 | 6.33 (6.2) | 32.9 | 62.7 |
| 200 | 2 | 20 | 20 | 1.22 | 0.81 | 8.18 | 2.41 | 4.39 | — | — | 11.15 | 32.45 | 1.50 | 21.9/20.3/22.6 | 42.3 | 10.10 (5.8) | 155.6 | 218.2 |
| 200 | 2 | 50 | 20 | 0.55 | 0.75 | 3.71 | 2.58 | 1.29 | 0.31 | 2.10 | 5.51 | 12.05 | 1.67 | 20.2/17.6/18.7 | 8.8 | 9.59 (5.7) | 67.5 | 103.7 |
| 200 | 10 | 5 | 4 | 0.40 | 0.69 | 4.00 | 0.63 | 0.22 | -0.05 | 0.45 | 1.33 | 5.71 | 1.60 | 2.4/1.5/2.0 | 0.5 | 11.18 (6.3) | 13.0 | 37.8 |
| 200 | 10 | 20 | 4 | 0.43 | 0.44 | 4.08 | 0.72 | 0.20 | 0.10 | 0.51 | 1.18 | 5.85 | 0.37 | 2.4/1.6/2.2 | 0.5 | 3.25 (6.3) | 11.9 | 23.6 |
| 200 | 10 | 50 | 4 | 0.43 | 0.44 | 5.46 | 0.61 | 0.19 | -1.85 | 0.71 | 1.22 | 5.78 | 0.40 | 2.3/1.5/2.0 | 0.5 | 2.92 (6.3) | 10.9 | 20.6 |
| 1000 | 2 | 5 | 5 | 0.26 | 0.25 | 17.04 | 4.27 | 0.28 | — | — | 1.14 | 25.88 | 0.36 | 10.1/8.0/9.9 | 0.6 | 2.53 (6.2) | 63.3 | 84.4 |
| 1000 | 2 | 20 | 20 | 0.15 | 0.14 | 10.15 | 8.60 | 0.77 | — | — | 3.16 | 23.91 | 0.21 | 18.8/19.8/20.4 | 6.2 | 1.40 (5.7) | 134.8 | 149.9 |
| 1000 | 2 | 50 | 48 | 0.19 | 0.15 | 10.98 | 18.61 | 4.03 | 3.71 | -1.25 | 12.24 | 57.89 | 0.40 | 61.9/66.7/69.0 | 64.6 | 1.78 (5.7) | 397.5 | 431.8 |
| 1000 | 10 | 5 | 5 | 0.29 | 0.31 | 11.85 | 3.17 | 0.22 | — | — | 0.67 | 15.27 | 0.21 | 5.6/3.8/5.8 | 0.3 | 1.86 (6.2) | 34.9 | 46.2 |
| 1000 | 10 | 20 | 7 | 0.26 | 0.23 | 11.68 | 3.43 | 0.21 | -1.48 | -1.07 | 0.80 | 13.20 | 0.19 | 7.1/5.5/6.1 | 0.5 | 1.15 (6.0) | 38.7 | 41.8 |
| 1000 | 10 | 50 | 7 | 0.23 | 0.21 | 8.81 | 2.48 | 0.18 | 3.15 | -2.75 | 0.81 | 13.15 | 0.19 | 7.1/4.9/5.8 | 0.5 | 1.17 (6.0) | 37.8 | 40.9 |
| 5000 | 2 | 5 | 5 | 0.13 | 0.13 | 44.57 | 13.54 | 0.17 | — | — | 0.54 | 62.29 | 0.16 | 21.7/18.0/21.4 | 0.3 | 1.21 (6.2) | 146.8 | 187.5 |
| 5000 | 2 | 20 | 20 | 0.20 | 0.19 | 59.12 | 42.07 | 0.98 | — | — | 3.60 | 115.21 | 0.24 | 133.3/96.2/107.7 | 5.9 | 1.44 (5.7) | 778.4 | 987.3 |
| 5000 | 2 | 50 | 50 | 0.20 | 0.17 | 47.03 | 93.70 | 4.41 | — | — | 9.72 | 176.57 | 0.34 | 212.5/202.1/222.4 | 51.2 | 1.00 (5.8) | 1680.1 | 1782.7 |
| 5000 | 10 | 5 | 5 | 0.24 | 0.21 | 45.78 | 11.12 | 0.15 | — | — | 0.60 | 61.54 | 0.17 | 19.7/16.4/19.8 | 0.3 | 1.21 (6.1) | 143.8 | 184.7 |
| 5000 | 10 | 20 | 14 | 0.26 | 0.24 | 52.50 | 28.65 | 0.58 | -0.21 | 0.41 | 1.54 | 80.17 | 0.23 | 70.2/52.3/63.1 | 2.3 | 0.99 (5.8) | 425.8 | 627.5 |
| 5000 | 10 | 50 | 14 | 0.25 | 0.23 | 46.06 | 25.06 | 0.43 | 0.42 | 0.34 | 1.52 | 78.08 | 0.18 | 59.2/49.8/58.2 | 2.1 | 1.25 (5.8) | 422.9 | 579.4 |

*Share of a sweep,* each part weighed by its frequency per sweep over the timed batches. "Structure move" is the births, deaths, changes, the target's evaluations and the commits. "Rest" is the sweep less the parts.

| n | p | M_max | M | sweep µs | birth (column, Zᵀz, Cholesky) | death | change | target | commit | τ | composition | rest | structure move | moves/sweep b/d/c |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 200 | 2 | 5 | 5 | 62.7 | 0.0% | 4.8% | 15.3% | 1.5% | 0.3% | 10.1% | 52.5% | 15.6% | 21.8% | 0.00/0.52/0.48 |
| 200 | 2 | 20 | 20 | 218.2 | 0.8% (0.4, 0.1, 0.2) | 2.4% | 6.5% | 0.6% | 0.7% | 4.6% | 71.3% | 13.1% | 10.9% | 0.10/0.46/0.44 |
| 200 | 2 | 50 | 20 | 103.7 | 3.6% (1.2, 0.8, 0.4) | 1.8% | 4.0% | 1.4% | 0.7% | 9.2% | 65.1% | 14.3% | 11.4% | 0.32/0.34/0.34 |
| 200 | 10 | 5 | 4 | 37.8 | 4.3% (2.7, 0.4, 0.1) | 1.4% | 5.5% | 3.7% | 0.2% | 29.6% | 34.4% | 21.0% | 15.0% | 0.25/0.38/0.36 |
| 200 | 10 | 20 | 4 | 23.6 | 8.2% (5.2, 0.9, 0.3) | 1.7% | 8.9% | 1.4% | 0.5% | 13.8% | 50.4% | 15.1% | 20.7% | 0.30/0.34/0.36 |
| 200 | 10 | 50 | 4 | 20.6 | 8.6% (7.9, 0.9, 0.3) | 2.0% | 10.1% | 1.7% | 0.6% | 14.2% | 52.7% | 10.2% | 23.0% | 0.30/0.34/0.36 |
| 1000 | 2 | 5 | 5 | 84.4 | 0.0% | 0.7% | 15.3% | 0.4% | 0.2% | 3.0% | 75.0% | 5.4% | 16.5% | 0.00/0.50/0.50 |
| 1000 | 2 | 20 | 20 | 149.9 | 0.0% | 1.1% | 7.8% | 0.1% | 0.3% | 0.9% | 89.9% | -0.3% | 9.4% | 0.00/0.50/0.49 |
| 1000 | 2 | 50 | 48 | 431.8 | 2.6% (0.8, 1.3, 0.3) | 1.1% | 4.4% | 0.1% | 2.2% | 0.4% | 92.1% | -2.8% | 10.3% | 0.30/0.37/0.33 |
| 1000 | 10 | 5 | 5 | 46.2 | 0.0% | 0.7% | 16.6% | 0.4% | 0.2% | 4.0% | 75.5% | 2.4% | 18.0% | 0.00/0.50/0.50 |
| 1000 | 10 | 20 | 7 | 41.8 | 10.5% (9.3, 2.7, 0.2) | 0.7% | 10.2% | 0.4% | 1.5% | 2.8% | 92.6% | -18.8% | 23.4% | 0.33/0.34/0.32 |
| 1000 | 10 | 50 | 7 | 40.9 | 10.1% (7.2, 2.0, 0.1) | 0.7% | 10.4% | 0.4% | 1.5% | 2.9% | 92.4% | -18.2% | 23.0% | 0.33/0.34/0.32 |
| 5000 | 2 | 5 | 5 | 187.5 | 0.0% | 0.2% | 15.9% | 0.1% | 0.2% | 0.6% | 78.3% | 4.8% | 16.3% | 0.00/0.52/0.48 |
| 5000 | 2 | 20 | 20 | 987.3 | 0.2% (0.1, 0.1, 0.0) | 0.2% | 5.8% | 0.0% | 0.7% | 0.1% | 78.8% | 14.1% | 6.9% | 0.02/0.49/0.49 |
| 5000 | 2 | 50 | 50 | 1782.7 | 0.4% (0.1, 0.3, 0.0) | 0.3% | 4.7% | 0.0% | 1.0% | 0.1% | 94.2% | -0.7% | 6.4% | 0.05/0.48/0.47 |
| 5000 | 10 | 5 | 5 | 184.7 | 0.0% | 0.2% | 16.8% | 0.1% | 0.6% | 0.7% | 77.9% | 3.9% | 17.6% | 0.00/0.50/0.50 |
| 5000 | 10 | 20 | 14 | 627.5 | 4.2% (2.7, 1.5, 0.0) | 0.1% | 4.4% | 0.0% | 1.5% | 0.2% | 67.9% | 21.7% | 10.3% | 0.32/0.33/0.34 |
| 5000 | 10 | 50 | 14 | 579.4 | 4.2% (2.7, 1.5, 0.0) | 0.1% | 4.6% | 0.0% | 1.5% | 0.2% | 73.0% | 16.4% | 10.4% | 0.34/0.33/0.34 |

*Tempered Friedman,* chain 0's self-sized ladder (20 rungs), with the two runs:

| part | run c, µs | share of a scan | run d, µs | share |
| --- | --- | --- | --- | --- |
| a scan | 251 | 100% | 477 | 100% |
| a pair's swap ratio, 19 pairs every scan | 0.17 | 1.3% | 0.61 | 2.4% |
| a hand-over, two per accepted swap at most (4.96 accepted per scan) | 6.6 | 26.1% | 12.2 | 25.3% |
| the replicas' sweeps, per replica (the scan less the above) | 9.1 | 72.6% | 17.2 | 72.3% |
| an untempered sweep at b = 1 (moves and τ) | 8.3 | — | 13.6 | — |

*What M1 shows.*
- *The machine was loaded,* at 27–44% CPU from other processes: two python processes near a core each, Defender and OneDrive. Absolute times moved 1.5–2× between the two runs on the final code. The structure-move share held within 7 points at every grid point, and the composition share within 21.
- *The decomposition is good to about ±20% of a sweep.* "Rest" runs from −19% to +22% in run c, and from −16% to +33% in run d. Read the shares to that precision.
- *The composition dominates* at Thin = 1: 34–94% of a sweep, above 65% wherever n ≥ 1,000. Each retained draw re-evaluates every term's column for its mean (`BassStructureMoves.ColumnMean`), O(N·M) per draw.
- *The structure move is 6–23% of a sweep.* τ's block is 0.1–30%: above 5% only at n = 200, and 14–30% where p = 10 too and a sweep costs 20–40 µs. Its slice takes about six evaluations per step.
- *Lever 1.* A second structure move per sweep adds the structure share, a sweep's 6–23%. It adds none of the composition, τ or the scan's rest. So k moves cost about 1 + (k − 1)·s sweeps, with s between 0.06 and 0.23. This holds at Thin = 1. Under thinning, or with the column means cached, s grows toward the whole of the moves' cost.
- *Within a birth,* the column's evaluation and Zᵀz dominate at large n: 47 and 94 µs at n = 5,000, M = 50. The bordered Cholesky is at most 4–9 µs (run c, run d). A change costs about a birth's design work plus the reduced model (176 µs there). On grid points at the cap, the chain makes half deaths and half changes.
- *The prefix-sum cross products the draft names are priced at* the birth's Zᵀz, 0–2.7% of a sweep. They would also cut part of a change's 4–17%. Cost does not bound the levers; acceptance does.
- *Tempered,* the swap ratios are 1–2.4% of a scan. The hand-overs are 25–26%: about two per accepted swap, each 0.7 of a replica's sweep, which runs near the untempered sweep's cost. Several moves per replica sweep (lever 1, tempered) would shrink the hand-overs' share in proportion.
- *How strongly:* the shares' ordering (composition ≫ structure move > τ, except at the smallest fits) holds in both runs at every grid point. The single percentages are good to the ±20% rest.

#### M2. Which births fail, and why

*Births lost by cause,* share of 40,000 draws per fixture (four chains × 10,000). Cap, variance-inflation and exact-fit refusals were 0 everywhere.

| cause | Friedman | mixed categorical | direct functional |
| --- | --- | --- | --- |
| empty draw | 0.00% | 8.62% | 11.32% |
| refused: support | 15.29% | 24.01% | 17.62% |
| refused: a repeated term | 0.03% | 1.25% | 0.43% |
| refused: angle | 0.00% | 1.11% | 0.00% |
| rejected, Σ(1 − α) | 83.08% | 64.02% | 69.32% |
| accepted, Σα | 1.60% | 1.00% | 1.30% |
| α per non-empty attempt, chains 0–3 | 1.71, 2.08, 1.04, 1.58% | 1.30, 1.12, 1.22, 0.73% | 0.50, 3.95, 0.50, 0.93% |

*By tag, Friedman,* pooled (K is the continuous order; a term is inert if any factor uses x6–x10):

| tag | births | refused | mean α | Δℓ median | 95% | max |
| --- | --- | --- | --- | --- | --- | --- |
| K1 active | 8,960 | 4.5% | 0.0163 | −3.17 | −0.72 | 2.45 |
| K1 inert | 9,001 | 4.1% | 0.0141 | −3.04 | −0.79 | 0.69 |
| K2 active | 4,648 | 16.5% | 0.0109 | −3.14 | −1.07 | 3.52 |
| K2 inert | 9,513 | 16.8% | 0.0194 | −3.15 | −0.92 | 4.80 |
| K3 active | 1,130 | 36.9% | 0.0067 | −3.20 | −1.61 | 0.69 |
| K3 inert | 6,748 | 38.2% | 0.0186 | −3.20 | −1.28 | 3.43 |

*Δℓ's quantiles over the admitted births,* with the threshold the other parts set (the Δℓ at which the ratio is 0):

| fixture, set | 5% | 25% | median | 75% | 95% | 99% | max | threshold median (5–95%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman, active | −3.55 | −3.45 | −3.16 | −2.41 | −0.81 | 0.22 | 3.52 | 2.08 (1.84–2.55) |
| Friedman, inert | −3.54 | −3.36 | −3.11 | −2.52 | −0.89 | 0.33 | 4.80 | 1.98 (0.97–2.24) |
| mixed, all | −3.26 | −2.84 | −2.61 | −2.38 | −1.90 | −1.19 | 2.61 | 1.71 (1.12–2.23) |
| functional, all | −4.40 | −4.17 | −3.79 | −3.10 | −1.31 | 1.10 | 5.02 | 1.82 (1.12–2.35) |

*Per chain, the largest Δℓ against its threshold:* Friedman 3.52 (x1+0.14 · x2−0.65; threshold 3.01, ratio 0.51), 4.80 (x6 · x8, inert; 1.30), 2.65 (x6 · x7 · x8, inert; 0.47), 3.43 (x1 · x5 · x8, inert; 1.56). Births with α > 0.05 number 476, 770, 225 and 675 of about 8,450 admitted per chain. Mixed: 1.12, 2.61, 2.19, 0.07. Functional: 0.08, 5.02 (chain 1 holds 244 births with α > 0.5), 2.33, 0.18.

*Main effects and lever 3a,* Friedman, K = 1, about 1,800 births per variable:
- Per-variable mean α: x1 0.021, x2 0.024, x3 0.015, x4 0.0043, x5 0.017; inert x6 0.0065, x7 0.0047, x8 0.028, x9 0.023, x10 0.008.
- The use weights z_u ∝ 5 + appearances move the active share from 0.50 to 0.59–0.62.
- Applied to each chain's per-variable α, they move the main effects' mean α from 0.0191, 0.0178, 0.0070, 0.0178 to 0.0236, 0.0191, 0.0070, 0.0166: −7% to +24%.
- On the mixed fixture the use weights change the continuous main effects' α by at most 3%; on the functional fixture by −18% to +23%. Their single categorical and output variables get the whole weight already.

*Active main effects by knot,* Friedman, by the knot's quantile among its variable's distinct values and its sign:
- Mean α runs from 0.008 to 0.045 across the ten cells: (−, 0.6–0.8) 0.0445; (−, 0.4–0.6) 0.019; (+, 0.6–0.8) 0.018.
- Edge cells whose hinge points outward lose about 21% to support.

*Deaths,* each term of each state removed in turn:
- *Friedman:* 38 of 41 terms have death α < 10⁻⁶, and the median Δℓ is −234. The exceptions are chain 3's inert x6 · x8 (Δℓ −2.95, α 0.20) and an x1 · x2 (0.175), and chain 1's x1 · x2 (7.6 · 10⁻⁴).
- *Mixed:* 24 of 25 terms below 10⁻⁴. Chain 1 holds an x1 · x2 · x3 term whose removal gains 2.15 (α 1).
- *Functional:* 24 of 25 below 10⁻⁴. Chain 3's t main effect has α 0.32.

*What M2 shows.*
- *The births' loss is the Occam cost, not the inputs.* An admitted birth's Δℓ concentrates near −½ ln N: −3.1 at N = 500, −2.6 at 300, −3.8 at 1,200. That is what an uninformative column costs under the unit-information g-prior. The prior and proposal parts then ask about 2 nats more (thresholds 1.7–2.1). So a random birth carries a ratio near −5, and only the 3–9% of births with α > 0.05 (Friedman) carry any weight. Strong: four chains and 26,000–34,000 scored births per fixture agree, and every ratio's parts sum to the chain's.
- *Inert inputs are not where births are lost.* Active and inert main effects accept alike (0.0163 against 0.0141). Inert interactions accept more than active ones: an inert factor on an active hinge can still fit residual structure.
- *Lever 3a's gain is small.* Use weights move main-effect acceptance by −7% to +24% per chain. Every main-effect birth on an active input would gain ×1.07 over uniform, pooled. The draft's "at most double" overstates it, because the doubling assumed inert births never accept.
- *Knot placement within active inputs varies acceptance about five-fold,* the reading 3b's condition asks for. The best cell still accepts under 5%, and these are births into states whose active inputs are already fitted.
- *Every term is load-bearing.* Deaths fail because removing a term costs tens to hundreds of nats. Only a just-born surplus term or an inert interaction comes back out. This is the local mode the draft named. Strong: 91 terms over 12 states.
- *One-chain readings:* the functional fixture's pooled birth acceptance rests on chain 1 (3.95% against 0.50–0.93%). That chain's state admits x0 · t interactions with Δℓ up to 5.

#### M3. Changes by knot displacement

*Friedman,* 40,000 changes, every one a continuous hinge. The kept-sign bucket 0 is empty by construction: a redraw of the same sign and knot is the same term, read on the row "same factor".

| sign | displacement | changes | share | refused | mean α (se) | Δℓ median |
| --- | --- | --- | --- | --- | --- | --- |
| kept | 1 | 81 | 0.20% | 0% | 0.963 (0.010) | 0.00 |
| kept | 2–3 | 147 | 0.37% | 0% | 0.928 (0.010) | 0.02 |
| kept | 4–7 | 350 | 0.88% | 0% | 0.800 (0.013) | −0.09 |
| kept | 8–15 | 651 | 1.63% | 0% | 0.656 (0.014) | −0.30 |
| kept | 16–31 | 1,170 | 2.92% | 0% | 0.418 (0.012) | −1.37 |
| kept | ≥ 32 | 17,593 | 43.98% | 6.9% | 0.041 (0.001) | −47.0 |
| flipped | 0 | 41 | 0.10% | 14.6% | 0.0001 | −175.5 |
| flipped | 1–31 | 2,349 | 5.87% | 5–18% | 0.005–0.008 | −100 to −160 |
| flipped | ≥ 32 | 17,568 | 43.92% | 5.1% | 0.021 (0.001) | −82.4 |
| same factor | — | 50 | 0.12% | 0% | 1 (a null move) | 0 |
| all | — | 40,000 | 100% | 5.9% | 0.064 (0.001) | −56.7 |

*Mixed categorical,* kept-sign hinges: 1: 0.969 (n = 84); 2–3: 0.935 (146); 4–7: 0.804 (305); 8–15: 0.591 (641); 16–31: 0.328 (1,148); ≥ 32: 0.046 (9,571). Flipped, at any displacement: 0.02–0.08. Subsets by size change: −1: 0.097; 0: 5.6 · 10⁻¹¹; +1: 0.047. Redraws equal to the current factor: 6.6% of changes, mostly subsets under R's measure.

*Direct functional,* kept-sign hinges:
- On x0 and x1 (60 distinct values): 1: 0.62 (395); 2–3: 0.23 (716); 4–7: 0.074; 8–15: 0.008; ≥ 32: 10⁻⁸.
- On t (20 values): 1: 0.37 (632); 2–3: 0.042; 4–7: 0.005.
- Flips: ≤ 0.08. Equal redraws: 1.5%.

*Per chain, kept-sign displacement 1–3:* Friedman 0.921, 0.934, 0.952, 0.953 (n ≈ 55 each); mixed 0.934–0.961; functional 0.17–0.29.

*What M3 shows.*
- A knot step of 1–3 distinct values with the sign kept accepts 0.93–0.97 on Friedman and the mixed fit. Today's redraw, at its common displacement (≥ 32), accepts 0.04.
- The redraw lands within three steps, sign kept, in 0.6% of Friedman's changes. Acceptance falls off with distance: about ½ at 16–31 steps, 3–6% of n_v.
- On the functional fit the falloff is far steeper: 0.62 at one step of x's 60 values, 0.37 at one step of t's 20. A step there is a coarse move.
- Sign flips accept ≤ 0.02 (Friedman) at any displacement, so the local step's keeping the sign costs nothing.
- Strength: consistent over four chains per fixture, with se ≤ 0.015 in the small buckets. These are acceptances, not mixing: a one-step move barely moves the knot. The rows at 4–15 (0.6–0.8 on Friedman and mixed) price wider windows.
- S10's "changes accept 5%" includes the null moves of equal redraws, 6.6% of the mixed fit's changes.

#### M4. The chains' modes

Unnormalized log posterior: the log marginal plus the log structure prior, τ's density left out.

| fixture | chain | mean M (range) | median log post (5–95%) | below the best chain's median | structures | top three's share |
| --- | --- | --- | --- | --- | --- | --- |
| Friedman (R̂ of M 2.35) | 0 | 10.28 (10–12) | −933.9 (−959.3, −930.8) | 41.5 | 17 | 83.3% |
| | 1 | 8.21 (8–10) | −892.4 (−912.7, −890.5) | 0 | 7 | 95.3% |
| | 2 | 10.77 (10–12) | −926.2 (−942.2, −913.8) | 33.8 | 6 | 97.6% |
| | 3 | 9.17 (9–11) | −916.0 (−935.4, −913.1) | 23.6 | 15 | 88.9% |
| mixed (1.54) | 0 | 5.66 (5–7) | −309.9 (−313.5, −308.3) | 0 | 5 | 99.5% |
| | 1 | 6.08 (6–7) | −312.5 (−332.0, −309.2) | 2.5 | 6 | 97.8% |
| | 2 | 6.02 (6–7) | −331.9 (−333.9, −329.9) | 22.0 | 4 | 99.9% |
| | 3 | 5.50 (5–7) | −309.9 (−324.8, −307.1) | 0.0 | 7 | 92.7% |
| functional (2.55) | 0 | 5.04 (5–6) | 181.1 (177.1, 183.5) | 0 | 5 | 99.6% |
| | 1 | 7.06 (7–8) | 152.7 (147.7, 155.6) | 28.4 | 5 | 97.8% |
| | 2 | 5.12 (5–6) | 168.6 (160.6, 176.3) | 12.5 | 5 | 98.6% |
| | 3 | 7.16 (7–9) | 151.7 (135.2, 154.3) | 29.4 | 7 | 95.6% |

*The structures held,* by variables and orders:
- Friedman chain 1's top structure (79.6%): x1 + x2 + 2 × x1·x2 + 2 × x3 + x4 + x5.
- Chain 0's (75.0%) holds five x1·x2 interactions and no x2 main effect.
- Chain 3's top (82.7%) is chain 1's second, with three x1·x2, at −915.8.
- Chain 2 holds two x1 and two x2 main effects.
- Mixed chains 0, 1 and 3 share x1 + x2·z + 3 × x3·z (+ z). Chain 2 holds two x2·x3·z terms that none of the others visits.

*Reach,* the share of chain i's states at or above chain j's median:
- Friedman: chains 0, 2 and 3 never reach chain 1's (0.0%). Chain 1's states lie inside chain 0's range 0.1% of the time.
- Functional: chains 1 and 3 never reach chain 0's or chain 2's medians.
- Mixed: chain 2 never reaches the others' medians, and they hold 0–6.8% of their states in its range.

*What M4 shows.*
- *The chains do not hold modes of comparable posterior.* On Friedman and the functional fit, three of four chains sit 12–42 nats below the best chain's median. They sit in structures that differ in which interactions carry the surface, and never reach the best chain's level over 2,000 draws. On the mixed fit, three chains share one mode and one sits 22 nats below.
- *These are local modes, not one diffuse mode.* Together with M2's deaths (every term load-bearing), leaving such a mode means replacing a load-bearing term with its substitute at once. Single births and deaths each pay tens of nats on the way.
- *Strength:* strong for the gaps, which are decades of nats on 2,000 draws per chain.
- *A caveat:* this compares the densities of the states held, not the modes' masses. A mode's mass also scales with its volume (its knots' freedom), which this does not measure. A 24-nat gap would need about e²⁴ times the volume to close.

#### M5. The steep band

*Part 1 reproduces S10.* Each chain's pairs above 0.9 are S10's exactly: chain 0 [0.222, 0.319] 1.000 and [0.169, 0.222] 0.937; chain 1 [0.195, 0.232] 0.991; chain 2 [0.156, 0.192] 0.997; chain 3 [0.205, 0.251] 0.999.

*The two ladders* (round trips per scan over the sampling phase, against Syed et al.'s ideal from the pairs' rejections):

| chain | Part 1 rungs | Λ̂ = Σr | Σh | ideal | observed | Part 2 rungs | Σr | Σh | ideal | observed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 20 | 8.64 | 11.63 | 3.3 · 10⁻⁷ | 0 | 45 | 9.87 | 10.25 | 0.030 | 0 |
| 1 | 22 | 8.93 | 10.60 | 0.0039 | 0 | 44 | 10.17 | 10.53 | 0.030 | 0 |
| 2 | 20 | 8.98 | 11.02 | 0.0013 | 0 | 43 | 9.59 | 9.98 | 0.031 | 0 |
| 3 | 23 | 8.64 | 10.75 | 0.00057 | 0 | 45 | 9.20 | 9.50 | 0.034 | 0 |

*Per-rung rows through the band, chain 0* (M and ℓ, the untempered log marginal, as mean ± sd over 2,000 scans; r is the pair to the next hotter rung):

| Part 1: b | M | ℓ | r | | Part 2: b | M | ℓ | r |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.438 | 9.20 ± 0.75 | −808.9 ± 6.4 | 0.229 | | 0.40 | 8.29 ± 0.81 | −805.2 ± 8.6 | 0.053 |
| 0.379 | 9.07 ± 0.78 | −811.8 ± 7.4 | 0.330 | | 0.32 | 8.16 ± 0.86 | −818.7 ± 17.3 | 0.107 |
| 0.319 | 8.73 ± 0.75 | −816.6 ± 9.4 | **1.000** | | 0.26 | 7.92 ± 0.87 | −868.0 ± 47.2 | 0.206 |
| 0.222 | 7.07 ± 0.35 | −968.9 ± 14.5 | **0.937** | | 0.24 | 7.59 ± 0.94 | −906.6 ± 60.9 | 0.341 |
| 0.169 | 5.94 ± 0.72 | −1045.2 ± 38.0 | 0.761 | | 0.23 | 7.27 ± 0.92 | −942.9 ± 62.7 | 0.293 |
| 0.150 | 4.50 ± 0.79 | −1130.3 ± 40.3 | 0.382 | | 0.20 | 6.56 ± 0.80 | −1011.6 ± 40.3 | 0.160 |
| 0.122 | 4.15 ± 0.59 | −1156.2 ± 32.2 | 0.603 | | 0.17 | 6.24 ± 0.56 | −1033.9 ± 22.1 | 0.089 |
| 0.099 | 3.48 ± 0.97 | −1221.7 ± 63.1 | 0.481 | | 0.13 | 5.33 ± 1.11 | −1093.4 ± 55.7 | 0.325 |

*The other chains in Part 2* show the same shape:
- M falls from 8.3–8.6 at b = 0.40 to 6.1–6.2 at b = 0.16.
- ℓ's sd peaks at 59–66 at b ≈ 0.23–0.25, against 9–11 at 0.40 and 18–34 at 0.16.
- The band's largest pair rejection is 0.34, 0.28, 0.46 and 0.36 for chains 0–3.
- In Part 1, each chain's steep pair straddles the jump: chain 2's 0.192 rung holds M 8.5 and ℓ −862, its 0.156 rung M 6.1 and ℓ −1041.

*The Gaussian reading* (h with 1 − r floored at 1/4,000, so h(1.000) = 2.92; a pair straddling a window's edge counted by its share of Δb inside):

| chain | steep span | Part 1 Σr, Σh over it | Part 2 Σr over it | Part 1 Σr, Σh over [0.10, 0.40] | Part 2 Σr over [0.10, 0.40] |
| --- | --- | --- | --- | --- | --- |
| 0 | [0.169, 0.319] | 1.94, 4.40 | 2.90 | 4.06, 6.80 | 4.98 |
| 1 | [0.195, 0.232] | 0.99, 2.09 | 0.91 | 4.34, 5.77 | 5.15 |
| 2 | [0.156, 0.192] | 1.00, 2.39 | 0.31 | 4.06, 5.66 | 4.82 |
| 3 | [0.205, 0.251] | 1.00, 2.59 | 1.19 | 3.99, 5.81 | 4.68 |

*What M5 shows.*
- *The band is a finite barrier that the sizing under-reads,* not a regime boundary whose rungs stay saturated. At Δb = 0.01 no pair in [0.10, 0.40] rejects above 0.46 on any chain.
- *It is a steep, resolvable regime change:* M falls by about two terms and ℓ by about 210 over b ≈ 0.30–0.16. ℓ's spread peaks at b ≈ 0.24 (a specific-heat-like peak).
- *The coarse ladders put one rung on each side.* Over the window, the dense ladder reads Σr 4.7–5.2 where Part 1's Λ̂ read 4.0–4.3. The whole ladder rises from 8.6–9.0 to 9.2–10.2.
- *The Gaussian increment does not predict the finer Σr better than Σr itself.* Over the window it over-reads by 0.6–1.8 where Σr under-reads by 0.7–0.9. Over the steep spans it over-reads by 1.2–2.1 on every chain, where Σr misses by −1.0 to +0.7. Its saturated pairs' value is the floor's (2.92 at 2,000 scans), not the barrier's.
- *The window is the fair comparison.* The dense runs are fresh trajectories, and the jump sits at slightly different b: chain 2's moved from [0.156, 0.192] to [0.22, 0.25], so its Part 2 Σr over Part 1's span (0.31) misses its band.
- *The draft's chain 0 band,* "a barrier of at least 4.2 where Λ̂ counts 1.94", reads 2.90 at Δb = 0.01.
- *With the band opened, the local moves bound the run.* Part 2 raised the ideal round-trip rate to 0.030–0.034 per scan, about 60 round trips in 2,000 scans. All four chains completed none. This is the draft's efficiency factor, read directly.
- *Strength:* four chains agree on every point above. The dense ladders are fresh runs, not replays of Part 1: chain 0's cold rung held M 10.65 against Part 1's 9.62.

#### Run times

Release, four chains in parallel where a measurement has them; M1 single-threaded.
- *M1:* the grid 27–49 s per run (four runs, and one with tiered compilation off, for the diagnosis below); the tempered part 1–5 s.
- *M2:* about 3 s for the three fixtures. *M3:* about 2 s. *M4:* about 2 s, with each fixture's S10 fit for the reproduction check.
- *M5:* 45 s: Part 1 7.3 s, Part 2 38.2 s.

#### Accessors and suites

- No internal accessor was needed and no production file changed, so the defaults are untouched by construction.
- The Regression suite, rerun without the switch (the new tests return at once): 828/828 (816 before, plus the 12 new cases).

#### The review passes and their fixes

1. *M2's first read.*
   - The per-attempt α counted empty draws as accepted: an empty draw's default score has no refusal, so its α read as 1. Fixed: the empty draws are left out of the per-attempt row (Friedman had none; the other two fixtures did).
   - Terms printed with 0-based indices beside the 1-based predictor names. They now print by name, with signs and knot quantiles.
   - Order tags became K/Z/T, so they do not read as variables.
2. *M1's no-commit check at the cap* found no admissible birth to probe with. The check now probes a change and a death too, and a birth only below the cap.
3. *M1's bordered proposals.* A pool term's bordered proposal can be indefinite, and the design then goes strict, so the replica takes the quadratic form only from a positive definite proposal. Repeated terms are refused before their column, as in the design.
4. *M1's birth commit at the cap* cannot be taken at the state, yet the chain leaves the cap during sweeps. It is taken from one term below, a death committed first, untimed.
5. *M1's rebuild cadence.* A forced commit can fall on the factor's rebuild. It is detected per call (a rebuild records its drift; a strict accept clears it) and replaced by the rebuild's amortized share.
6. *M1's first full tables* showed "design other" negative below the cap, and "rest" up to about 50%. There were three causes:
   - The test-side replicas were timed as tier-0 code against the engine's optimized hot paths; a run with tiering off confirmed it. Fixed with `AggressiveOptimization` on the replicas.
   - The sweeps ran on from the state while M drifted upward. Each sweep batch now restarts from the state's snapshot, with its own streams.
   - The parts were timed minutes apart from the sweep under varying load. They now run round-robin.

   After the fixes, |rest| ≤ 22% in run c and ≤ 33% in run d.
7. *Correctness held throughout.*
   - No commit: the design, its factor and the counts came back bit for bit, after every M1, M2 and M3 harness.
   - Every scored ratio's parts summed to the chain's within 1e-10 relative.
   - Births' causes added to one and matched the moves' tallies.
   - The death harness matched the death move's proposals bit for bit.
   - M2–M4 reran on the final helper with identical readings.
8. *Coverage.*
   - Empty buckets: the kept-sign displacement 0 (by construction), and t's ≥ 32 buckets (20 values).
   - One-chain readings: the functional fixture's birth acceptance rests on chain 1.
   - M1 reads one chain per grid point, and one ladder (chain 0's) for its tempered part.
   - M1 ran on a loaded machine. A quiet rerun would tighten its absolute times.

#### For the main session

- *A duplicate.* The fixtures in `BassMixing.cs` copy `BassQualificationRecoveryTests`' private generators. M4's reproduction assertion guards the copy. Making the generators internal there would remove it, a change to a file not mine.
- *Readings that bear on the design draft,* which is not mine to edit:
  - Lever 3a's measured gain is −7% to +24% per chain, and ×1.07 pooled with every main effect on an active input, not ×2 (M2).
  - The band's barrier at chain 0 is 2.9, not "at least 4.2" (M5).
  - With the band opened, the round trips are bounded by the local moves (M5's Part 2).
  - Composition, not the structure move, dominates a sweep at Thin = 1, which sets lever 1's price (M1).
- *Not a recommendation:* caching each term's column mean in the composition belongs to the efficiency pass.
