# BASS F-M2: pricing the mode exits (chip brief)

Date: 2026-10-29. Follows Phase 0 ([its brief and report](local-move-measurements-brief.md#reports)) and the design draft's "After Phase 0: proposals" section in [local-move-mixing-brief.md](local-move-mixing-brief.md). The user decided on 2026-10-29 to price these hypotheses before choosing Phase 1's levers; this brief is that work. Like Phase 0, it measures and changes nothing.
- *What Phase 0 found.*
  - The chains hold local modes 12–42 nats apart (M4). Each worse mode carries co-adapted terms that are each load-bearing: on Friedman, 38 of 41 terms have death α below 10⁻⁶ (M2).
  - A local knot step accepts about 0.95 (M3).
  - With the steep band opened, the ladder's ideal round-trip rate reached 0.03 per scan and no round trip completed (M5).
- *The hypotheses this chip prices,* P1–P4 below:
  - co-adapted pairs of terms that are jointly redundant;
  - structural moves at fixed M that re-represent a term;
  - a drift path along which sliding knots make a term redundant;
  - mid-ladder rungs frozen by b·Δℓ.
- *Not in scope.*
  - Building any move into production code.
  - Any default.
  - The design draft's questions, which are the user's.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the design draft's sections "Phase 0's readings" and "After Phase 0: proposals", which define the moves and their ratios;
- Phase 0's brief and report in full: its harness, its no-commit check, and its fixed states;
- the fitting brief's D01 (the term prior, the birth proposal, the move schedule).

In the code:
- *Phase 0's harness,* `tests/regression/BassMixing.cs`: fixtures, burned chains and ladders, `Score`, and the no-commit check. Its tests `BassMixing*Tests.cs` show how each was used.
- *The moves and the design:* `BassTermProposal.cs` (`LogBirth` through `LogBirthRatio` and `LogDeathRatio`), `BassSweep.cs` (`BassStructureMoves`), `BassDesign.cs` (`ProposeChange`, `ProposeDeath`, `Load`) and `BassTarget.cs` (`Evaluate`, `LogStructurePrior`; the tempered forms' likelihood at b).

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-exits\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-exits\` when you finish.
- **The main session works beside you,** caching the composition's column means in `src/maths/regression/spline/bass/BassSweep.cs`. That change is bit-identical and does not touch the structure moves. Leave the file alone. If a build breaks while its edit is in progress, wait and retry.
- **Your files.**
  - New test files `tests/regression/BassMixingExit*Tests.cs`.
  - Additions to `tests/regression/BassMixing.cs`, without changing what Phase 0's tests read: rerun M2–M4 and check that their readings are unchanged.
  - This brief's Reports section.
  - Leave everything else alone, the design draft included; name any change you need in your report.
- **Production code stays as it is.** The candidate moves are test-side scorers, not moves. Where a measurement cannot be taken without a change, the only one allowed is an internal accessor that changes no behaviour; name it, and rerun the Regression suite.
- **The harness changes nothing it measures.** Apply Phase 0's no-commit check after every scorer.
- **Commit your own files only, by explicit path,** one commit per measurement or a coherent pair. Use the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF.
- **Read, do not judge.** Each harness prints its readings and asserts only its own soundness. Recommend no lever; the readings feed the user's design pass.

## The measurements

- *Fixed states.* Phase 0's: each chain's state at the end of a 5,000-sweep burn-in, on S10's Friedman, mixed categorical and direct functional fits, four chains each.
- *Better and worse modes,* by M4's medians: Friedman's best is chain 1, and chains 0, 2 and 3 sit 24–42 nats below; on the mixed fit chain 2 sits 22 below the rest; on the functional fit chains 1 and 3 sit about 29 below chains 0 and 2.
- *First, check those states are in their chains' modes.* Compare each state's log posterior with its chain's M4 median. Report it, and if a state lies outside its chain's sampling range, use the chain's most frequent structure instead.
- *Report everything per chain,* marking which are the worse modes.

**P1. Pair deaths** (`BassMixingExitPairTests`).
- *For every unordered pair of terms in each state,* record:
  - Δℓ for removing both;
  - the pair death's full log ratio. That is the target's ratio (log marginal and log structure prior, through the target) plus the proposal part ln q_pair(t₁, t₂ | S′) + ln C(M, 2), where S′ is the reduced structure and q_pair(t₁, t₂ | S′) = q(t₁ | S′)·q(t₂ | S′ + t₁) + q(t₂ | S′)·q(t₁ | S′ + t₂), each q the D01 birth probability (`LogBirth`). Take the move-selection ratio as 1, and say so.
- *Against the single deaths:* for each pair, the two terms' single-death Δℓ. A pair whose joint removal costs far less than either alone is co-adapted.
- *Where M ≤ 12, also every triple's Δℓ* for the states whose best pair costs under 20 nats.
- *Report:* each state's best pairs (Δℓ, full ratio, the terms by name) and the count of pairs within 10 and 20 nats. Say whether the worse modes hold jointly redundant pairs the best mode does not.

**P2. Structural moves at fixed M** (`BassMixingExitFactorTests`).
- *Prune:* every factor of every term with K ≥ 2 in its part. Record Δℓ and the full log ratio, which the design draft derives as Δℓ alone.
  - Assert that the target's term-prior ratio and the proposal ratio cancel, to 1e-12 relative. This checks the derivation, so report it if it fails.
  - A prune onto a duplicate term, or below the support count, is a refusal: report its share.
- *Grow and variable change:* 5,000 proposals each per state, drawn as the draft defines them (grow: the variable uniform among the term's unused ones, the sign and knot from the prior; change: v → v′ uniform among the unused ones, the knot from v′'s prior). Record Δℓ, the cancellation check, and refusals.
- *Group:*
  - by part (continuous, categorical, functional);
  - on Friedman, by whether the factor's variable is active or inert;
  - by the term's order before and after.
- *Two-step profile:* after each prune or variable change costing under 20 nats, record the best single death and the best pair death in the new state. The question is whether a cheap re-representation opens a cheap removal.
- *Report:* acceptance (min(1, e^{ratio})) and Δℓ quantiles per group; each worse-mode state's ten cheapest exits by name; the two-step profile.

**P3. Profile deaths** (`BassMixingExitProfileTests`).
- *For each term of each state:* remove it, then re-optimize the remaining hinges' knots greedily.
  - In passes over the factors, try j ± δ for δ ∈ {1, …, 8}. Keep any step that raises the log marginal plus the log structure prior; the knot prior is uniform, so this is the marginal.
  - Stop after 5 passes or a pass without gain.
  - Categorical subsets stay fixed.
- *Record:* the plain death's Δℓ; the profile's Δℓ (after re-optimization, against the original state); the passes taken; and the knots' total displacement in distinct-value steps.
- *Report:* per state, plain against profile Δℓ for each term, and how many terms' profile Δℓ lies within 10 and 20 nats.
  - This is an optimization diagnostic: it bounds where a drift path could end, not that one exists or how long it takes.
  - Say so in the reading.

**P4. Structure moves along the dense ladder** (`BassMixingExitLadderTests`).
- *The ladder:* M5's Part 2 dense fixed ladder on Friedman, each chain, at the end of its burn-in.
- *At each rung's state, scored under that rung's tempered target:*
  - births (2,000 per rung);
  - every term's death;
  - change redraws (2,000);
  - local knot steps, δ ∈ {±1, ±2, ±3}, for every hinge factor;
  - and P2's prunes.
  - Use the ladder's own tempered form (the integrated form by default), through the target's evaluation at b, with the no-commit check per rung.
- *Report:* per rung (b, M), the mean α of each move kind, and the deaths' median Δℓ scaled by b. Say at which b deaths stop being negligible (α above 10⁻³), and whether the band's rungs are frozen.

## Review before reporting

Before you write the report, take a fresh-eyes pass over the harnesses and the readings, and fix what it finds:
- *Correctness:*
  - the no-commit check holds;
  - P1's q_pair sums both draw orders;
  - P2's cancellation check passes, or its failure is reported;
  - P4 scores at the rung's b and not at 1;
  - the counts add up.
- *Coverage:*
  - a state not in its chain's mode;
  - a reading that rests on one chain;
  - P3's greedy search stopping at its pass limit rather than at a local optimum.

Name the passes and their fixes in the report.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- each measurement's tables and its reading, kept apart: the numbers, then what they show, and how strongly;
- the run times;
- any accessor added, and the suites rerun;
- the review passes;
- anything the main session must know.

## Reports

Newest first.
