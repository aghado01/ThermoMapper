# The flat-histogram warm start, robust across seeds, and Friedman replicated (chip brief)

Date: 2026-10-30. N17 Phase 2 found Friedman's posterior on one seed ([its brief and reports](flat-histogram-phase2-brief.md#reports)). On master seed 8017, from the WHAM warm start, every chain tunnelled, the chains agreed within 0.34 nats, and ln Ẑ(1) = −841.4 ± 0.2 reached the −864.6 bound, with the posterior on a1 (four x1·x2 terms).
- *A second seed (8217) failed at its start.* SAMC's fixed point refused the pilot's ĝ over [−1552.56, −773.559], K = 128, b_lo = 0. No node was empty, but a node's occupancy stood 125 from flat, against the 10⁻³ the solve must reach.
- *A flat start there, at the same budget, learned nothing:* no round trip in 100,000 iterations, b = 1 Kish ESS 8 of 500,000.
- *What this chip does:*
  - make the warm start hold across seeds, without changing a start the fixed point already accepts;
  - replicate the Friedman reading on seed 8217 as pre-registered, and on a third seed;
  - read whether a flat start can learn Friedman's window with a larger gain scale.
- *Not in scope:*
  - any default change;
  - the window's top, which a pilot that never found a1 sets below the posterior's ℓ₁ (a later option);
  - SPC (Phase 3);
  - `docs/Maths` and the checklists, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- [the design](flat-histogram-design.md), with its "Phase 1, as built" and "Phase 2, as built";
- [Phase 1's brief](flat-histogram-phase1-brief.md), its reports (the `FixedPoint` fix `9f561d0` above all);
- [Phase 2's brief](flat-histogram-phase2-brief.md), its reports and the main session's verification;
- `docs/Maths/Samplers.md`, "Flat-histogram sampling", and `docs/Maths/BayesianRegression.md`, "BASS: the flat-histogram sampler".

In the code:
- `src/maths/samplers/flathistogram/SamcLearning.cs` (`FixedPoint`: the free nodes, the empty-node rule, the start fitted per node cell, Newton, the stall rule and the 10⁻³ refusal);
- `Wham.cs` (`Estimate`, `WarmStart`, `Window`, `Nodes`) and `TemperingPilot.cs`;
- `src/maths/regression/spline/bass/BassFlatHistogramRun.cs` (the pilot, `Pooled`, the start's shape, the refusal's message);
- `tests/regression/BassFlatHistogramQualificationTests.cs`: arm (i), `Friedman_ArmOne_SecondSeed_TwiceTheSweeps` (seed 8217, `Thin` 400, its criteria written before its run) and `Friedman_ArmOne_SecondSeed_FlatStart`.

## The steps

Each step lands with its checks and ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next. Name the passes and their fixes in the report.

1. **Diagnose seed 8217's refusal.**
   - Rebuild its pilot (the pilot is deterministic from the master seed; it takes about a minute), then WHAM's ĝ and the start's shape.
   - Find the node that departs and why:
     - its hat's mass and where in its cell that mass sits;
     - its neighbours' masses;
     - the pilot's draws per segment;
     - Newton's path: the iterations, the step sizes, where it stalled or ran out.
   - Say which of these it is, with the numbers:
     - an ill-posed limit, a hat seeing mass only at an edge of its support, so θ* is extreme;
     - a well-posed limit Newton fails to reach from its start;
     - something else.
   - Do the same for seeds 7801–7803, where Phase 2's exploration also saw refusals on fixture A.
2. **Make the start robust.**
   - *Build* the remedy the diagnosis supports. The candidates:
     - widen the empty-node rule to thin nodes, a hat whose effective pilot mass falls below a stated threshold, interpolated as empty nodes are;
     - interpolate ln ĝ log-linearly across segments the pilot left empty, before the fixed point;
     - coarsen K until the fixed point flattens;
     - a better Newton start or line search, if the limit is well posed.
   - Choose by measured robustness, and state the rule and its threshold.
   - *The invariant:* a start the fixed point accepts today is unchanged, bit for bit. Check it on seed 8017's Friedman pilot, by the start's weights' hash taken before the change, and on every Phase 1 and Phase 2 test.
     - So the remedy runs where the fixed point refuses, or under a rule that no accepted start triggers.
     - Seed 8017's accepted pilot already had nodes under 10 draws (0–2 and 153–155), so a thin-node threshold must not reach them, or must apply only on refusal.
   - *Check:*
     - on a sweep of Friedman pilots (at least 16 master seeds, declared before the sweep, 8017 and 8217 among them), the refusal rate before and after;
     - on the toys and fixture A, the start's accuracy against SAMC's limit where it was accepted, unchanged, and where it was refused, now within Phase 1's tolerance (0.25 RMS);
     - the refusal kept for a ĝ with no usable mass, with its message.
3. **Replicate Friedman** (`BASS_QUALIFICATION` = 1; S10's Friedman fixture; Release; one fit at a time).
   - *Seed 8217* by `Friedman_ArmOne_SecondSeed_TwiceTheSweeps` as written: the WHAM start, the default gain start, `Thin` 400, its criteria as committed in `56c66f3`, unchanged.
   - *A third seed,* declared before any run and not used by any harness: arm (i) at the same settings, with the same criteria, in a test written and committed before it runs.
   - *What would count,* for each seed, as committed:
     - every production chain tunnels;
     - the chains' ln Ẑ_c(1) SD is at most 1 nat, and each labelled basin at 5% or more passes the χ² at p ≥ 0.01;
     - the pooled ln Ẑ(1) is at or above −864.6;
     - against seed 8017, ln Ẑ(1) lies within 1 nat plus both runs' errors of −841.36 ± 0.17, and a1's pooled share within 5 points of 65.9%.
   - *A reading, not ruled:* the flat start on seed 8217 over its pilot's window with t₀ = 100·K (Liang et al.'s range runs to 100K), the rest as `Friedman_ArmOne_SecondSeed_FlatStart`. It reads whether a flat start can learn the window when the warm start is unavailable.
   - *Read per run,* by the harness's `Read`: the window, K, t₀ and the gain's start; the pilot and the fixed point; the learning's round trips, flatness and trace; production's round trips, flatness and b = 1 Kish ESS; ln Ẑ(1) per chain and pooled; the basins with their χ²; σ's interval; D10's verdict; the time.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Other sessions may be working in this checkout,** an SPC session among them. If a build fails in a file you do not own, wait and retry; never edit it.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-fh3\` to every `dotnet build` and `dotnet test`. Build in Release for the runs. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-fh3\` when you finish.
- **Your files:**
  - `src/maths/samplers/flathistogram/` (the warm start's code) and its tests under `tests/samplers/`;
  - `BassFlatHistogramRun.cs`, only where the start's shape or message must follow the change;
  - `tests/regression/BassFlatHistogramQualificationTests.cs` and `BassFlatHistogramTests.cs`;
  - this brief's Reports section.
- **Commit your own files only, by explicit path,** one commit per step, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked notes under `issues/regression/` and other sessions' edits stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
  - Commit each run's test before running it.
- **Line endings:** briefs and C# are LF.
- **Bit-identity:** `BassBitIdentityTests`' hashes hold, and so does every start the fixed point accepts today.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - the remedy leaves an accepted start bit for bit as it was;
  - a thin or interpolated node's weight is a start only, which SAMC's learning corrects;
  - the refusal kept where nothing usable remains;
  - each run's criteria committed before it ran.
- *Coverage:*
  - a pilot with an empty segment beside a drawn node;
  - a thin node at the window's edge;
  - a pilot whose ĝ spans more than 745 nats;
  - a seed whose replication fails a criterion, reported as it falls.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the diagnosis, with its numbers;
- the remedy and its rule, with the refusal rates before and after;
- each replication: what would count as committed, then its tables, then what they show, kept apart;
- the flat-start reading;
- the suites with their counts (Samplers, Regression, the BASS oracle classes);
- the review passes and their fixes;
- anything the main session must know.

Recommend no default.

## Reports

Newest first.
