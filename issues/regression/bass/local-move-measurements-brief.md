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
