# BASS F-M2 Phase 2: the pre-registered pass (chip brief)

Date: 2026-10-29. Follows Phase 1, the local-move options, built and verified the same day ([its brief and report](local-move-options-brief.md#reports)). The user decided the pass on 2026-10-29: the design draft's section "Phase 1's decisions" in [local-move-mixing-brief.md](local-move-mixing-brief.md) fixes its arms, fixtures, seeds and rule, and this brief fixes the rest before any run. The user said "Start Phase 2" the same day.
- *What it does:* runs the eight arms on four fixtures, five master seeds each, and applies the rule exactly as written below.
- *What it does not do:*
  - change any default, or any production code; the main session verifies the verdict and then carries it out;
  - size the run lengths or the ladder, which follow the pass;
  - judge beyond the rule: other readings are context, kept apart from the verdict.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the design draft's sections "Phase 1's decisions" and "Phase 2";
- Phase 1's report, for the options' configuration, their tallies (`MoveKinds`) and its smoke harness;
- S10's qualification in the fitting checklist's reports, for the fixtures and its readings.

In the code:
- *S10's fixtures and configuration:* `tests/regression/BassQualificationRecoveryTests.cs` (`Config`, `Friedman`, `Mixed`, `Functional`).
- *Phase 1's smoke harness:* `tests/regression/BassLocalMoveSmokeTests.cs`, behind `BASS_QUALIFICATION` = 1, which times a fit on one thread.
- *Phase 0's harness:* `tests/regression/BassMixing.cs` (`LogPosterior`, `Loaded`).
- *The diagnostics:* `BassDiagnostics` (rank-normalized R̂ and bulk ESS, the monitored points, `MoveKinds`, the tempering readout), and the shared diagnostics functions it calls.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-pass\` to every `dotnet build` and `dotnet test`. Build in Release for the timed runs, and state the configuration. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-pass\` when you finish.
- **Your files.**
  - A new harness, `tests/regression/BassMovePassTests.cs`, behind `BASS_QUALIFICATION` = 1.
  - In `BassQualificationRecoveryTests.cs`, a sample-size parameter for `Friedman`. The n = 500 data must stay bit for bit; check it.
  - The per-fit results, `issues/regression/bass/local-move-pass-results.csv`.
  - This brief's Reports section.
  - Leave everything else alone, production code included. If a reading cannot be taken without a production change, stop and say so in the report.
- **Commit your own files only, by explicit path,** in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked files under `issues/regression/` (`argon-bass-qa.md`, `argon-bass-question-suggestions.md`, `bass-questions.md`) stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs, CSV and C# are LF.
- **Timing needs a quiet machine.** Run the pass in the background, and run nothing else heavy while it runs. The main session will not run builds or tests meanwhile.

## The pass, fixed before any run

**Arms.** Each is S10's configuration (`Config`: four chains, `BurnIn` 5,000, `Samples` 2,000, `Thin` 1) with these changes:

| arm | `Moves` | `LocalSteps` | `FactorMoves` | `Tempering` |
| --- | --- | --- | --- | --- |
| A0, the baseline | `EqualOverAvailable` | null | null | null |
| A1, local steps | `EqualOverAvailable` | `new()` | null | null |
| A2, the fixed-M family | `EqualOverAvailable` | null | `new()` | null |
| A3, both | `EqualOverAvailable` | `new()` | `new()` | null |
| A4, prior-cancelling | `PriorCancelling()` | null | null | null |
| A5, prior-cancelling with both | `PriorCancelling()` | `new()` | `new()` | null |
| A6, tempered | `EqualOverAvailable` | null | null | `new TemperingSpec { Split = new SplitRule() }` |
| A7, tempered with both | `EqualOverAvailable` | `new()` | `new()` | `new TemperingSpec { Split = new SplitRule() }` |

The tempered arms count `BurnIn` and `Samples` in scans, as S10's tempered run did. Equal work is the metric's: ESS per CPU-second puts a scan's cost and a sweep's on one scale.

**Fixtures.**
- S10's three: Friedman (n = 500, p = 10), the mixed categorical fit (n = 300) and the direct functional fit (n = 60, m = 20).
- Friedman at n = 2,000, from `Friedman`'s generator with its data seed 51 and n = 2,000.

**Seeds.** Five master seeds per fixture, s₀ + 1,000·k for k = 0, …, 4. s₀ is S10's for its three fixtures (53, 57, 61), and 65 for Friedman at n = 2,000. Every arm uses the same five seeds, so the arms are paired. Starts are the fit's own exact prior draws (S7).

**Running.**
- Each fit runs with `MaxDegreeOfParallelism` = 1, one fit at a time, in Release. Run the harness alone (`--filter` on its class), so xUnit runs nothing beside it.
- Within each (fixture, seed), the eight arms run in an order shuffled by a fixed seed, so slow drift in the machine falls on every arm alike.
- The time is the whole `Bass.Fit`, burn-in and sampling. Record the sampling phase's share too, if the fit exposes it without a production change.
- *A timing control:* for one seed per fixture, run A0 a second time at the end, and report the two times' ratio.

**Readings, per fit.**
- The rank-normalized bulk ESS, and R̂, of:
  - M;
  - σ²;
  - the log posterior, the collapsed log π(structure, u | y) of each retained draw at b = 1, from its structure and τ;
  - the fitted value at each monitored point (the default points).
- *The fit's minimum ESS,* over those quantities. A quantity constant over every draw has no ESS; leave it out of the minimum and report it.
- *min-ESS/s:* the minimum ESS over the fit's time in seconds.
- *Context, outside the rule:*
  - the acceptance per move kind, from `MoveKinds`;
  - each chain's median M;
  - whether σ's 95% interval covers the truth (Friedman: 1);
  - for the tempered arms, each chain's settled rung count, the steep pairs, the round trips per scan, and observed over ideal.

**The rule, verbatim from the decisions.**
- Per fixture f, an arm's ratio is r_f = exp(mean over seeds of ln(min-ESS/s of the arm) − mean over seeds of ln(min-ESS/s of A0)). The seeds are paired, so this is also the geometric mean of the per-seed ratios.
- The arm's score is the geometric mean of r_f over the four fixtures.
- An arm clears the rule when its score is at least 1.5 and no r_f lies below 0.8.
- Among arms that clear it, the highest score wins, and its changes become the defaults. If none clears, every option stays an option and the defaults stay D01's.
- No arm, fixture or seed is dropped after the runs.

**Context for the verdict,** reported beside it and not part of it:
- the spread of the per-seed log ratios: their standard deviation over seeds, and the standard error of r_f's log;
- each arm's score with any one fixture left out.

## Checks before the runs

- *The arms are what they say:* each fit's resolved configuration (`BassResolved`: the step windows, the factor shares, the split threshold) and the fingerprint differ across arms exactly as the table does.
- *The pairing holds:* the eight arms of one (fixture, seed) start from the same exact prior draws. Compare the starting structures across arms.
- *The log posterior is right:* for a few retained draws, the computed value equals the chain's own log density at that state.
- *Friedman's n = 500 data is unchanged* by the new parameter.
- Run the whole pass once at 1/10 of the lengths on one seed, and check that every reading is produced and the rule's arithmetic runs, before the timed run.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - the rule computed exactly as written: the logs, the paired seeds, the geometric means, the floor;
  - the ESS taken from the shared rank-normalized estimator;
  - the log posterior at b = 1;
  - no fit timed with another running beside it;
  - the CSV agreeing with the report's tables.
- *Coverage:*
  - a fit whose minimum ESS is near zero, or whose quantity is constant;
  - a verdict that rests on one fixture (the leave-one-out scores);
  - the timing control's spread, against the margins the rule turns on.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the verdict first: each arm's score, its r_f per fixture, and which arm clears, if any;
- then the context tables: per arm and fixture, min-ESS/s (geometric mean over seeds), the quantity that set the minimum, R̂ of each quantity, the acceptance per move kind, σ's coverage, and the tempered readings;
- the timing control, the run times, and the build configuration;
- the review passes and their fixes;
- anything the main session must know.

## Reports

Newest first.
