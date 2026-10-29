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

### Phase 2 run: no arm clears — 2026-10-29

*Landed,* on `main`, nothing pushed:
- `6f1e99a`: the harness, `tests/regression/BassMovePassTests.cs`, and `Friedman`'s sample size in `BassQualificationRecoveryTests`, committed before the timed run.
- `e45cdc0`: the per-fit results, `local-move-pass-results.csv` (160 fits and 4 timing controls, one row each), and a check of the minimum's edge cases, added after the run (it reads nothing of the pass).
- This report.

**The verdict.** Under the rule as written, over all 160 fits, **no arm clears**. No arm, fixture or seed was dropped. Every option stays an option, and the defaults stay D01's.

| arm | score | r_f Friedman | r_f mixed | r_f functional | r_f Friedman, n = 2,000 | least r_f | clears |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A1, local steps | 0.929 | 0.823 | 0.953 | 1.016 | 0.933 | 0.823 | no |
| A2, fixed-M family | 1.031 | 1.003 | 0.986 | 0.961 | 1.191 | 0.961 | no |
| A3, both | 0.998 | 0.944 | 1.180 | 0.891 | 0.998 | 0.891 | no |
| A4, prior-cancelling | 0.864 | 0.815 | 0.793 | 0.830 | 1.040 | 0.793 | no |
| A5, prior-cancelling with both | 0.984 | 0.826 | 1.023 | 0.980 | 1.131 | 0.826 | no |
| A6, tempered | 0.078 | 0.071 | 0.156 | 0.059 | 0.057 | 0.057 | no |
| A7, tempered with both | 0.067 | 0.060 | 0.128 | 0.059 | 0.043 | 0.043 | no |

- The highest score is A2's, 1.031: its mean ln r_f is 0.03, where a score of 1.5 needs 0.405.
- One r_f lies below the 0.8 floor: A4's on the mixed fit, 0.793. A4 also fails on its score, 0.864.
- An independent script re-scored the CSV (from each fit's minimum ESS and seconds, the minimum re-derived from the quantities). It reproduced every r_f, score, standard deviation, standard error and leave-one-out score below.

*Context for the verdict,* not part of it. The per-seed log ratios' standard deviation, and the standard error of ln r_f:

| arm | Friedman: sd, SE | mixed: sd, SE | functional: sd, SE | Friedman, n = 2,000: sd, SE |
| --- | --- | --- | --- | --- |
| A1, local steps | 0.228, 0.102 | 0.433, 0.194 | 0.434, 0.194 | 0.199, 0.089 |
| A2, fixed-M family | 0.159, 0.071 | 0.550, 0.246 | 0.364, 0.163 | 0.107, 0.048 |
| A3, both | 0.216, 0.097 | 0.388, 0.174 | 0.479, 0.214 | 0.138, 0.062 |
| A4, prior-cancelling | 0.141, 0.063 | 0.295, 0.132 | 0.237, 0.106 | 0.191, 0.086 |
| A5, prior-cancelling with both | 0.137, 0.061 | 0.201, 0.090 | 0.339, 0.152 | 0.110, 0.049 |
| A6, tempered | 0.248, 0.111 | 0.742, 0.332 | 0.474, 0.212 | 0.148, 0.066 |
| A7, tempered with both | 0.197, 0.088 | 0.214, 0.096 | 0.474, 0.212 | 0.425, 0.190 |

Each arm's score with one fixture left out:

| arm | without Friedman | without mixed | without functional | without Friedman, n = 2,000 |
| --- | --- | --- | --- | --- |
| A1, local steps | 0.967 | 0.921 | 0.901 | 0.927 |
| A2, fixed-M family | 1.041 | 1.047 | 1.056 | 0.983 |
| A3, both | 1.016 | 0.943 | 1.036 | 0.997 |
| A4, prior-cancelling | 0.881 | 0.890 | 0.876 | 0.813 |
| A5, prior-cancelling with both | 1.043 | 0.971 | 0.985 | 0.939 |
| A6, tempered | 0.081 | 0.062 | 0.086 | 0.087 |
| A7, tempered with both | 0.069 | 0.054 | 0.070 | 0.077 |

- No verdict rests on one fixture. With any one fixture left out, the untempered arms score 0.81–1.06 and the tempered arms 0.05–0.09.

*What the metric read.* In the untempered fits, the minimum ESS sits at the estimator's floor.
- Over the 120 untempered fits, the median minimum ESS is 4.96 (range 4.25–10.8). 119 of them read below 10, and 101 below 6. M's R̂ has a median of 1.96–2.45 per arm, and reaches 4.34.
- The chains hold their own bands of M, as S10 found, so the 8 split chains sit at different levels. Stan's bulk ESS of split chains that disagree tends to about their count over two, 4 here, which is where these minima sit.
- For the untempered arms, the rule's ratio is therefore mostly a ratio of times. The table splits each r_f into the arm's minimum-ESS ratio over A0's and its time ratio (r_f is the first over the second; geometric means over seeds, from the CSV):

| arm | Friedman | mixed | functional | Friedman, n = 2,000 |
| --- | --- | --- | --- | --- |
| A1, local steps | 1.04 / 1.26 | 1.06 / 1.11 | 1.02 / 1.00 | 1.01 / 1.08 |
| A2, fixed-M family | 1.03 / 1.03 | 1.10 / 1.12 | 0.88 / 0.92 | 1.09 / 0.92 |
| A3, both | 1.06 / 1.12 | 1.14 / 0.97 | 0.87 / 0.98 | 0.97 / 0.97 |
| A4, prior-cancelling | 0.87 / 1.06 | 0.94 / 1.18 | 0.82 / 0.99 | 0.96 / 0.92 |
| A5, prior-cancelling with both | 0.91 / 1.10 | 1.14 / 1.11 | 0.98 / 1.00 | 0.96 / 0.85 |
| A6, tempered | 2.01 / 28.42 | 3.76 / 24.16 | 1.85 / 31.23 | 2.30 / 40.28 |
| A7, tempered with both | 1.79 / 29.62 | 3.17 / 24.70 | 1.83 / 31.17 | 1.85 / 42.64 |

- Tempering lifts the minimum ESS 1.8–3.8-fold and brings M's R̂ to a median of 1.23–1.25 (at most 1.55), at 24–43 times A0's time.
- *A0 reproduces S10.* At S10's seeds (k = 0), A0 reads M's R̂ of 2.35, 1.54 and 2.55. σ's intervals are [1.013, 1.161], [0.472, 0.558] and [0.188, 0.203], S10's to the digit. Friedman's chain medians are 10, 8, 11 and 9, as in Phase 1's smoke reading.

*The context tables.* Below, lp is the log posterior, and f1–f3 are the fitted values at the default monitored points.

Min-ESS/s, the geometric mean over seeds, with the quantity that set each seed's minimum in brackets. No quantity was constant in any fit, and M never froze.

| arm | Friedman | mixed | functional | Friedman, n = 2,000 |
| --- | --- | --- | --- | --- |
| A0, baseline | 22.59 (M 4, lp 1) | 25.46 (M 3, lp 2) | 26.71 (M 3, lp 2) | 7.13 (M 2, lp 2, f1 1) |
| A1, local steps | 18.60 (M 2, f1 2, lp 1) | 24.27 (M 3, lp 2) | 27.13 (lp 3, M 2) | 6.65 (M 4, lp 1) |
| A2, fixed-M family | 22.66 (M 2, f3 1, lp 1, f2 1) | 25.10 (M 2, f1 2, lp 1) | 25.66 (M 3, lp 2) | 8.49 (M 4, lp 1) |
| A3, both | 21.33 (M 3, lp 1, f1 1) | 30.04 (M 3, f1 2) | 23.78 (M 4, lp 1) | 7.12 (M 4, f3 1) |
| A4, prior-cancelling | 18.42 (lp 4, M 1) | 20.19 (M 4, lp 1) | 22.16 (M 4, lp 1) | 7.42 (lp 4, M 1) |
| A5, prior-cancelling with both | 18.65 (lp 4, f1 1) | 26.04 (M 3, lp 2) | 26.18 (M 3, lp 2) | 8.07 (lp 4, M 1) |
| A6, tempered | 1.60 (M 4, f1 1) | 3.96 (M 5) | 1.58 (lp 5) | 0.41 (M 3, f3 1, lp 1) |
| A7, tempered with both | 1.36 (M 3, f1 1, f3 1) | 3.27 (f1 2, M 2, lp 1) | 1.57 (lp 4, M 1) | 0.31 (f1 2, lp 1, f2 1, M 1) |

The minimum ESS itself, the geometric mean over seeds [the least, the largest]:

| arm | Friedman | mixed | functional | Friedman, n = 2,000 |
| --- | --- | --- | --- | --- |
| A0, baseline | 5.4 [5.0, 6.1] | 5.1 [4.4, 6.6] | 5.6 [4.4, 8.4] | 4.9 [4.4, 5.4] |
| A1, local steps | 5.6 [4.9, 7.0] | 5.4 [4.3, 9.0] | 5.7 [4.3, 10.8] | 4.9 [4.2, 5.8] |
| A2, fixed-M family | 5.6 [5.0, 6.2] | 5.6 [4.4, 8.0] | 5.0 [4.3, 6.1] | 5.3 [5.1, 5.5] |
| A3, both | 5.7 [4.9, 7.0] | 5.8 [4.3, 8.8] | 4.9 [4.3, 6.8] | 4.7 [4.5, 5.1] |
| A4, prior-cancelling | 4.7 [4.4, 5.0] | 4.8 [4.4, 5.6] | 4.6 [4.3, 4.9] | 4.7 [4.4, 5.0] |
| A5, prior-cancelling with both | 4.9 [4.7, 5.2] | 5.8 [5.2, 6.7] | 5.5 [4.9, 6.5] | 4.7 [4.4, 4.9] |
| A6, tempered | 10.8 [7.6, 14.8] | 19.1 [8.4, 47.0] | 10.4 [7.5, 21.6] | 11.2 [9.1, 12.5] |
| A7, tempered with both | 9.6 [8.1, 12.4] | 16.1 [10.1, 22.5] | 10.3 [6.6, 20.3] | 9.0 [7.1, 14.7] |

R̂ by quantity: the diagnostics' own, the less favourable of bulk and tail. Each cell gives the median over seeds, with the largest in brackets.

| arm | fixture | M | sigma2 | logpost | f1 | f2 | f3 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A0, baseline | Friedman | 1.97 [2.35] | 1.08 [1.14] | 1.99 [2.05] | 1.73 [1.84] | 1.32 [1.50] | 1.33 [1.96] |
| A0, baseline | mixed | 2.50 [3.49] | 1.01 [1.01] | 2.03 [2.65] | 1.45 [1.71] | 1.02 [1.07] | 1.04 [1.09] |
| A0, baseline | functional | 2.15 [3.40] | 1.00 [1.02] | 2.02 [2.39] | 1.13 [1.18] | 1.48 [1.50] | 1.12 [1.14] |
| A0, baseline | Friedman, n = 2,000 | 2.20 [3.36] | 1.10 [1.22] | 2.32 [2.75] | 1.58 [1.98] | 1.49 [1.82] | 1.21 [1.31] |
| A1, local steps | Friedman | 2.02 [2.17] | 1.06 [1.16] | 1.82 [2.34] | 1.72 [1.90] | 1.28 [1.48] | 1.42 [1.75] |
| A1, local steps | mixed | 2.40 [3.78] | 1.02 [1.04] | 1.89 [2.29] | 1.27 [1.56] | 1.04 [1.16] | 1.14 [1.29] |
| A1, local steps | functional | 1.98 [4.34] | 1.01 [1.01] | 1.84 [2.77] | 1.05 [1.19] | 1.36 [1.58] | 1.05 [1.26] |
| A1, local steps | Friedman, n = 2,000 | 2.05 [4.31] | 1.07 [1.09] | 2.40 [3.09] | 1.56 [1.80] | 1.51 [1.66] | 1.46 [1.52] |
| A2, fixed-M family | Friedman | 1.74 [2.32] | 1.06 [1.23] | 1.91 [2.01] | 1.68 [1.93] | 1.30 [1.72] | 1.50 [1.72] |
| A2, fixed-M family | mixed | 1.67 [3.42] | 1.03 [1.04] | 1.41 [1.98] | 1.88 [2.15] | 1.08 [1.29] | 1.02 [1.17] |
| A2, fixed-M family | functional | 2.61 [4.15] | 1.01 [1.03] | 2.17 [2.99] | 1.15 [1.19] | 1.43 [1.48] | 1.07 [1.12] |
| A2, fixed-M family | Friedman, n = 2,000 | 2.00 [2.22] | 1.04 [1.18] | 1.62 [2.19] | 1.41 [1.65] | 1.55 [1.62] | 1.55 [1.85] |
| A3, both | Friedman | 1.73 [2.45] | 1.13 [1.33] | 1.73 [2.08] | 1.76 [1.85] | 1.32 [1.50] | 1.40 [1.88] |
| A3, both | mixed | 1.73 [3.81] | 1.01 [1.02] | 1.35 [2.12] | 1.37 [1.88] | 1.04 [1.55] | 1.11 [1.23] |
| A3, both | functional | 3.73 [4.05] | 1.00 [1.02] | 2.20 [2.86] | 1.18 [1.41] | 1.37 [1.73] | 1.11 [1.36] |
| A3, both | Friedman, n = 2,000 | 2.85 [3.12] | 1.07 [1.14] | 2.08 [2.66] | 1.36 [1.97] | 1.46 [1.60] | 1.35 [2.23] |
| A4, prior-cancelling | Friedman | 2.32 [3.23] | 1.74 [2.63] | 2.77 [3.40] | 1.91 [2.43] | 1.41 [2.00] | 1.97 [2.10] |
| A4, prior-cancelling | mixed | 2.43 [3.65] | 1.03 [1.76] | 2.17 [2.96] | 1.69 [2.48] | 1.03 [1.21] | 1.03 [1.31] |
| A4, prior-cancelling | functional | 2.89 [4.09] | 1.01 [1.02] | 2.07 [2.76] | 1.08 [1.14] | 1.18 [1.30] | 1.14 [1.16] |
| A4, prior-cancelling | Friedman, n = 2,000 | 2.35 [3.34] | 1.68 [1.76] | 2.81 [3.16] | 1.51 [1.69] | 1.77 [1.88] | 1.66 [1.79] |
| A5, prior-cancelling with both | Friedman | 1.75 [2.68] | 1.48 [1.71] | 2.29 [2.74] | 1.87 [2.66] | 1.45 [1.68] | 1.33 [1.95] |
| A5, prior-cancelling with both | mixed | 1.81 [2.10] | 1.02 [1.10] | 1.64 [2.17] | 1.30 [1.43] | 1.03 [1.11] | 1.07 [1.09] |
| A5, prior-cancelling with both | functional | 1.65 [2.48] | 1.01 [1.05] | 1.79 [2.17] | 1.08 [1.17] | 1.09 [1.50] | 1.07 [1.11] |
| A5, prior-cancelling with both | Friedman, n = 2,000 | 2.45 [2.79] | 1.49 [1.61] | 2.66 [3.39] | 1.64 [1.92] | 1.79 [1.91] | 1.74 [1.80] |
| A6, tempered | Friedman | 1.23 [1.49] | 1.04 [1.10] | 1.17 [1.26] | 1.16 [1.23] | 1.10 [1.19] | 1.10 [1.24] |
| A6, tempered | mixed | 1.16 [1.41] | 1.00 [1.01] | 1.11 [1.37] | 1.08 [1.13] | 1.00 [1.01] | 1.01 [1.02] |
| A6, tempered | functional | 1.02 [1.37] | 1.00 [1.00] | 1.36 [1.49] | 1.01 [1.02] | 1.03 [1.09] | 1.01 [1.06] |
| A6, tempered | Friedman, n = 2,000 | 1.24 [1.37] | 1.01 [1.02] | 1.26 [1.38] | 1.14 [1.27] | 1.11 [1.19] | 1.09 [1.35] |
| A7, tempered with both | Friedman | 1.34 [1.43] | 1.04 [1.10] | 1.29 [1.35] | 1.15 [1.40] | 1.09 [1.34] | 1.15 [1.31] |
| A7, tempered with both | mixed | 1.15 [1.50] | 1.00 [1.01] | 1.14 [1.29] | 1.10 [1.30] | 1.02 [1.04] | 1.01 [1.05] |
| A7, tempered with both | functional | 1.18 [1.45] | 1.00 [1.01] | 1.38 [1.61] | 1.01 [1.02] | 1.07 [1.20] | 1.02 [1.02] |
| A7, tempered with both | Friedman, n = 2,000 | 1.20 [1.55] | 1.01 [1.03] | 1.19 [1.30] | 1.28 [1.38] | 1.12 [1.45] | 1.10 [1.35] |

Acceptance by move kind, pooled over seeds: accepted over attempted, with accepted over admitted (neither null nor refused) in brackets.

| arm | fixture | Birth | Death | Redraw | LocalStep | Toggle | Prune | Grow | VariableChange |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A0, baseline | Friedman | 0.0205 (0.0243) | 0.0204 (0.0204) | 0.0602 (0.0637) |  |  |  |  |  |
| A0, baseline | mixed | 0.0106 (0.0148) | 0.0101 (0.0101) | 0.1394 (0.1517) |  |  |  |  |  |
| A0, baseline | functional | 0.0085 (0.0107) | 0.0079 (0.0079) | 0.0457 (0.0509) |  |  |  |  |  |
| A0, baseline | Friedman, n = 2,000 | 0.0158 (0.0167) | 0.0150 (0.0150) | 0.0436 (0.0442) |  |  |  |  |  |
| A1, local steps | Friedman | 0.0188 (0.0220) | 0.0190 (0.0190) | 0.0566 (0.0599) | 0.7944 (0.8178) |  |  |  |  |
| A1, local steps | mixed | 0.0122 (0.0172) | 0.0114 (0.0114) | 0.1476 (0.1610) | 0.8186 (0.8434) | 0.0066 (0.0100) |  |  |  |
| A1, local steps | functional | 0.0076 (0.0095) | 0.0069 (0.0069) | 0.0419 (0.0466) | 0.3735 (0.4157) |  |  |  |  |
| A1, local steps | Friedman, n = 2,000 | 0.0120 (0.0128) | 0.0123 (0.0123) | 0.0402 (0.0408) | 0.7167 (0.7424) |  |  |  |  |
| A2, fixed-M family | Friedman | 0.0204 (0.0240) | 0.0203 (0.0203) | 0.0646 (0.0681) |  |  | 0.0072 (0.0235) | 0.0076 (0.0085) | 0.0076 (0.0080) |
| A2, fixed-M family | mixed | 0.0121 (0.0170) | 0.0109 (0.0109) | 0.1277 (0.1388) |  |  | 0.0038 (0.0066) | 0.0041 (0.0091) | 0.0013 (0.0034) |
| A2, fixed-M family | functional | 0.0128 (0.0161) | 0.0118 (0.0118) | 0.0421 (0.0473) |  |  | 0.0019 (0.0034) | 0.0027 (0.0048) | 0.0000 (0.0000) |
| A2, fixed-M family | Friedman, n = 2,000 | 0.0174 (0.0184) | 0.0181 (0.0181) | 0.0556 (0.0564) |  |  | 0.0026 (0.0067) | 0.0059 (0.0062) | 0.0034 (0.0034) |
| A3, both | Friedman | 0.0186 (0.0219) | 0.0191 (0.0191) | 0.0554 (0.0584) | 0.8001 (0.8208) |  | 0.0095 (0.0312) | 0.0080 (0.0090) | 0.0082 (0.0087) |
| A3, both | mixed | 0.0110 (0.0156) | 0.0097 (0.0097) | 0.1585 (0.1720) | 0.8462 (0.8661) | 0.0134 (0.0206) | 0.0042 (0.0090) | 0.0058 (0.0128) | 0.0036 (0.0097) |
| A3, both | functional | 0.0062 (0.0078) | 0.0059 (0.0059) | 0.0515 (0.0576) | 0.3734 (0.4125) |  | 0.0027 (0.0051) | 0.0036 (0.0062) | 0.0004 (0.0010) |
| A3, both | Friedman, n = 2,000 | 0.0142 (0.0149) | 0.0147 (0.0147) | 0.0493 (0.0501) | 0.7079 (0.7296) |  | 0.0062 (0.0173) | 0.0041 (0.0043) | 0.0036 (0.0037) |
| A4, prior-cancelling | Friedman | 0.0882 (0.1036) | 0.0145 (0.0145) | 0.0723 (0.0766) |  |  |  |  |  |
| A4, prior-cancelling | mixed | 0.0654 (0.0946) | 0.0122 (0.0122) | 0.1453 (0.1596) |  |  |  |  |  |
| A4, prior-cancelling | functional | 0.0314 (0.0404) | 0.0063 (0.0063) | 0.0466 (0.0525) |  |  |  |  |  |
| A4, prior-cancelling | Friedman, n = 2,000 | 0.0775 (0.0823) | 0.0111 (0.0111) | 0.0558 (0.0568) |  |  |  |  |  |
| A5, prior-cancelling with both | Friedman | 0.0861 (0.1005) | 0.0140 (0.0140) | 0.0637 (0.0675) | 0.7884 (0.8124) |  | 0.0061 (0.0228) | 0.0063 (0.0071) | 0.0094 (0.0100) |
| A5, prior-cancelling with both | mixed | 0.0528 (0.0757) | 0.0102 (0.0102) | 0.1518 (0.1642) | 0.8369 (0.8585) | 0.0025 (0.0037) | 0.0075 (0.0144) | 0.0073 (0.0164) | 0.0023 (0.0062) |
| A5, prior-cancelling with both | functional | 0.0295 (0.0373) | 0.0063 (0.0063) | 0.0435 (0.0491) | 0.3412 (0.3891) |  | 0.0017 (0.0034) | 0.0021 (0.0038) | 0.0012 (0.0027) |
| A5, prior-cancelling with both | Friedman, n = 2,000 | 0.0768 (0.0816) | 0.0111 (0.0111) | 0.0558 (0.0566) | 0.7031 (0.7245) |  | 0.0047 (0.0145) | 0.0090 (0.0093) | 0.0051 (0.0052) |
| A6, tempered | Friedman | 0.0178 (0.0209) | 0.0179 (0.0179) | 0.0569 (0.0599) |  |  |  |  |  |
| A6, tempered | mixed | 0.0124 (0.0172) | 0.0124 (0.0124) | 0.1384 (0.1495) |  |  |  |  |  |
| A6, tempered | functional | 0.0065 (0.0082) | 0.0057 (0.0057) | 0.0421 (0.0469) |  |  |  |  |  |
| A6, tempered | Friedman, n = 2,000 | 0.0130 (0.0137) | 0.0123 (0.0123) | 0.0445 (0.0452) |  |  |  |  |  |
| A7, tempered with both | Friedman | 0.0186 (0.0219) | 0.0198 (0.0198) | 0.0743 (0.0787) | 0.7980 (0.8221) |  | 0.0035 (0.0103) | 0.0049 (0.0055) | 0.0084 (0.0088) |
| A7, tempered with both | mixed | 0.0114 (0.0160) | 0.0097 (0.0097) | 0.1456 (0.1573) | 0.8482 (0.8658) | 0.0029 (0.0043) | 0.0018 (0.0036) | 0.0044 (0.0098) | 0.0025 (0.0063) |
| A7, tempered with both | functional | 0.0075 (0.0093) | 0.0061 (0.0061) | 0.0512 (0.0575) | 0.3213 (0.3585) |  | 0.0022 (0.0041) | 0.0023 (0.0040) | 0.0005 (0.0010) |
| A7, tempered with both | Friedman, n = 2,000 | 0.0114 (0.0120) | 0.0117 (0.0117) | 0.0515 (0.0524) | 0.7341 (0.7539) |  | 0.0026 (0.0071) | 0.0033 (0.0034) | 0.0041 (0.0041) |

σ's 95% interval covers the truth:

| arm | Friedman | mixed | functional | Friedman, n = 2,000 |
| --- | --- | --- | --- | --- |
| A0, baseline | 2 of 5 | 5 of 5 | 5 of 5 | 4 of 5 |
| A1, local steps | 0 of 5 | 5 of 5 | 5 of 5 | 3 of 5 |
| A2, fixed-M family | 2 of 5 | 5 of 5 | 5 of 5 | 4 of 5 |
| A3, both | 1 of 5 | 5 of 5 | 5 of 5 | 5 of 5 |
| A4, prior-cancelling | 0 of 5 | 5 of 5 | 5 of 5 | 0 of 5 |
| A5, prior-cancelling with both | 0 of 5 | 5 of 5 | 5 of 5 | 0 of 5 |
| A6, tempered | 5 of 5 | 5 of 5 | 5 of 5 | 5 of 5 |
| A7, tempered with both | 5 of 5 | 5 of 5 | 5 of 5 | 5 of 5 |

The tempered arms' ladders over seeds:
- each chain's settled rungs;
- its steep pairs, those whose rejection exceeds 0.9;
- the round trips per scan, and the ideal rate;
- observed over ideal;
- the largest pair rejection.

No tempered fit completed a round trip.

| arm | fixture | rungs per chain [range] | steep pairs per chain (mean) | round trips/scan (mean) | ideal (mean) | observed/ideal (mean) | largest pair rejection |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A6, tempered | Friedman | 20.4 [18, 24] | 1.10 | 0 | 0.0067 | 0 | 1.000 |
| A6, tempered | mixed | 15.8 [14, 17] | 0.10 | 0 | 0.039 | 0 | 0.911 |
| A6, tempered | functional | 20.0 [17, 24] | 0.20 | 0 | 0.0274 | 0 | 0.974 |
| A6, tempered | Friedman, n = 2,000 | 30.3 [27, 38] | 1.15 | 0 | 0.00532 | 0 | 1.000 |
| A7, tempered with both | Friedman | 20.4 [19, 23] | 0.85 | 0 | 0.00784 | 0 | 1.000 |
| A7, tempered with both | mixed | 16.0 [14, 18] | 0.10 | 0 | 0.0386 | 0 | 0.952 |
| A7, tempered with both | functional | 21.0 [19, 26] | 0.15 | 0 | 0.0276 | 0 | 0.998 |
| A7, tempered with both | Friedman, n = 2,000 | 29.9 [26, 35] | 1.25 | 0 | 0.0047 | 0 | 1.000 |

Each chain's median M is in the CSV (`median_M`).

*Times and the build.*
- *The build.* Release (`-c Release`), net10.0, 18 logical processors, workstation GC.
- *How each fit ran.* Each fit ran alone on one thread (`MaxDegreeOfParallelism` = 1), under xUnit filtered to `BassMovePassTests.Pass`. Within each (fixture, seed), the arms ran in an order shuffled from seed 20261029.
- *The warm-up.* Before any timed fit, every arm ran untimed on every fixture at 1/20 of the lengths.
- *What was timed.* The time is the whole `Bass.Fit`.
- *The sampling share* is read from the run's progress reports, from the last burn-in sweep to the last sweep, including the ensemble's own diagnostics. It is 0.30–0.52 untempered and 0.25–0.36 tempered.
- *Duration.* The pass took 9.0 minutes in all, 8.6 of them in timed fits. The dry run took 17 seconds.

Seconds per fit, the mean over seeds, with the sampling share's and CPU over wall time's ranges over the arm's 20 fits:

| arm | Friedman | mixed | functional | Friedman, n = 2,000 | sampling share | CPU/wall |
| --- | --- | --- | --- | --- | --- | --- |
| A0, baseline | 0.24 | 0.20 | 0.21 | 0.69 | 0.30–0.45 | 0.83–1.17 |
| A1, local steps | 0.31 | 0.22 | 0.21 | 0.74 | 0.31–0.45 | 0.80–1.14 |
| A2, fixed-M family | 0.25 | 0.23 | 0.19 | 0.63 | 0.31–0.48 | 0.77–1.10 |
| A3, both | 0.27 | 0.19 | 0.21 | 0.67 | 0.31–0.47 | 0.91–1.05 |
| A4, prior-cancelling | 0.25 | 0.24 | 0.21 | 0.64 | 0.32–0.50 | 0.90–1.04 |
| A5, prior-cancelling with both | 0.27 | 0.22 | 0.21 | 0.58 | 0.33–0.52 | 0.78–1.14 |
| A6, tempered | 6.79 | 4.83 | 6.68 | 27.79 | 0.25–0.36 | 0.87–0.98 |
| A7, tempered with both | 7.10 | 4.99 | 6.61 | 29.96 | 0.29–0.36 | 0.90–0.99 |

The timing control, A0 at each fixture's first seed, run again at the end:

| fixture | first (s) | again (s) | again/first | draws identical |
| --- | --- | --- | --- | --- |
| Friedman | 0.248 | 0.225 | 0.910 | yes |
| mixed | 0.198 | 0.184 | 0.930 | yes |
| functional | 0.187 | 0.172 | 0.920 | yes |
| Friedman, n = 2,000 | 0.640 | 0.555 | 0.868 | yes |

- *Drift is by block, about ±10%.* Seed block k = 3 ran 8.6% above its cells' geometric means, and k = 4 ran 10.5% below. The arms of a block share their block's drift: they are shuffled within it, and the rule pairs the seeds.
- *Descheduling.* 11 of the 160 fits ran at CPU over wall below 0.9 (0.77 at least), spread across seven arms. A `dllhost` process was using about one core of the 18 when the machine was checked before the run.
- *Against the margins, no verdict turns on timing.*
  - The best score would need to rise 1.45-fold.
  - The floor's one miss, A4 on the mixed fit at 0.793, lies within the timing control's spread of 0.8, but A4 fails on its score regardless.

*Checks before the runs,* all passed:
- *The arms* (`Arms_AreWhatTheySay`, run always):
  - each arm's schedule, local steps, family and tempering are its row of the table;
  - its resolved echo follows: windows max(1, ⌈0.02·n_v⌉) on each variable's knots, shares ⅓, ⅙, ⅙, ⅓, threshold 0.9, each null where the option is off;
  - its lengths and seed are S10's;
  - the eight fingerprints of each fixture are distinct, and each equals S10's once its four fields are set back.
- *The pairing.* In all 20 (fixture, seed) blocks, the eight arms' starts hash alike: each chain's terms, τ and the prior draw's attempts. The pass asserted it block by block.
- *The log posterior* (`LogPosterior_IsTheChainsOwn`, gated): A0, A5 and A7 on each fixture, two chains each.
  - Each chain was replayed as the run takes it, the ladder burned and scanned for A7, and read at three draws.
  - The structure and τ match the draw's bit for bit.
  - The chain's own b = 1 density equals `Temperature.LogPosterior` at the state bit for bit.
  - The value read from the draw lies within 1.2·10⁻¹¹ of it.
- *Friedman's n = 500 data* (`Friedman_AtN500_IsS10s`, run always) hashes `6b52ce789967ee0c` before and after the parameter. The first 500 cases at n = 2,000 are S10's bit for bit. `BassBitIdentityTests`' eight pins hold.
- *The shared estimator.* In every fit, the fit's own diagnostics for M and σ² equal `ChainDiagnostics.Diagnose` on the retained draws; the pass asserted it.
- *The dry run,* at 1/10 of the lengths on one seed, produced every reading, and the rule's arithmetic ran.

*Choices the brief left open,* each settled before the timed run and stated in the code:
- *A frozen M* (each split chain constant, at different values) would read the ESS of its mean (`ChainDiagnostics.MeanEss`). It never occurred.
- *The untimed warm-up.*
- *The timing control's seed* is k = 0 for every fixture, and its draws are asserted identical.
- *R̂* is the diagnostics' `RHat`.
- *The filter* is the pass alone, a subset of its class.

*The review passes and their fixes:*
- *Correctness:*
  - the rule was re-scored independently, as above;
  - the ESS comes from the shared estimator, asserted in every fit;
  - the log posterior is at b = 1, checked against the chains' own densities;
  - no fit was timed beside another: one test, one fit at a time, on one thread;
  - the compact tables here are generated from the CSV.
- *One gap, fixed.* The dry run exercised neither a constant quantity nor a frozen M, so the minimum's two edge paths had no check. `Minimum_LeavesAConstantOut_AndReadsAFrozenMByItsMean` (synthetic chains, run always) was added in `e45cdc0`, and the class summary names it.
- *Coverage:*
  - near-zero minima: the floor reading above;
  - constant quantities: none;
  - a verdict resting on one fixture: none, by the leave-one-out scores;
  - the timing spread: above.

*For the main session:*
- *The verdict to verify:* no arm clears, so no default changes.
- *The metric's floor.* At S10's lengths, the pre-registered metric reads every untempered arm at the bulk-ESS floor, so its ratio there mostly compares times. An untempered option could clear only by lifting the chains off that floor, and none does at these lengths.
- *The run lengths.* The draft's next step sizes them "from the winner's ESS per sweep", and there is no winner. A0's ESS per sweep at these lengths is floor-bound too. Where the sizing starts from is the main session's call.
- *Tempering, as context.* The tempered arms cover σ in all 40 of their fits and hold M's R̂ to at most 1.55, against 10–16 of each untempered arm's 20 fits and R̂ up to 4.34. Yet none of their 40 fits completed a round trip, and they took 24–43 times as long.
- *Prior-cancelling (A4, A5)* misses σ on both Friedman fixtures in all 20 of those fits. Its births accept 3.5–6.2 times as often as A0's, fixture by fixture, while its deaths accept at 0.006–0.015 against A0's 0.008–0.020.
- *Times.* Here A0 takes about 7–9 µs per chain sweep at n ≤ 500, and about 25 at n = 2,000 (the whole fit over its sweeps). Phase 1's smoke reading gave 35–62, and S10's cost baseline about 100 on Friedman. Those runs' build configurations are not recorded, which matters to any sizing by time.
- `artifacts\chip-pass\` is deleted.
