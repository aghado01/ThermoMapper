# BASS F-M2 Phase 1: the local moves as options (chip brief)

Date: 2026-10-29. Follows Phase 0 ([its brief and report](local-move-measurements-brief.md#reports)) and the mode-exit pricing ([its brief and report](mode-exit-pricing-brief.md#reports)). The user decided Phase 1 on 2026-10-29; the decisions are the design draft's section "Phase 1's decisions" in [local-move-mixing-brief.md](local-move-mixing-brief.md). This chip builds them as options, off by default, and proves each one exact.
- *What it builds:*
  - lever 2, local steps inside a change;
  - the fixed-M family: a factor's prune and grow, and a factor's variable change;
  - lever 4's split rule, 4b, in the shared `TemperingLadder`.
- *Not in scope:*
  - the pass (Phase 2), which follows once this report is verified, and any default;
  - levers 1 and 3, the pair death and F-X4's seam, all withdrawn or deferred by the decisions;
  - `docs/Maths`, the design draft and the checklists, which the main session updates after verifying;
  - the R oracles. Under a constant likelihood the target is the prior on F whatever the kernel, so F-O3's laws serve unchanged.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then, in the design draft:
- "Lever 2" and "Lever 4" (the kernels and 4b);
- "After Phase 0: proposals" (the fixed-M family's ratios) and "The mode-exit pricing's readings";
- "Phase 1's decisions", which govern where they differ from the earlier sections, and "Phase 1" (the checks).

Also [the fitting brief](port-fitting-brief.md)'s D01 (the term prior, the birth proposal, the move schedule) and D05 (tempering), and the pricing report's choices (its matched selection, and the "prune, part emptied" case).

In the code:
- *The moves:* `BassTermProposal.cs` (`DrawChange`, `DrawHinge`, the birth and death ratios), `BassSweep.cs` (`BassStructureMoves`: `Weight`, `Birth`, `Death`, `Change`, `Candidate`, the `Support` tally), `BassStructurePrior.cs` (`LogPrior(term)`, each part's order law, `LogFactors`, the subset measure's `LogSubset` and `DrawSubset`), `BassDesign.cs` (`ProposeChange`) and `BassTarget.cs` (the tempered forms).
- *The configuration:* `BassConfig.cs` (`BassSampling.Moves`, `TemperingSpec`, validation), `BassResolved.cs`, `BassDiagnostics.cs`.
- *The ladder:* `src/maths/samplers/rjmcmc/TemperingLadder.cs` (`Barrier`, `Size`, `Place`, `Sources`), `BassLadder.cs` (`Adapt`, `Resize`) and `ParallelTempering.cs`. BARS's `BarsEnsemble.cs` calls the same ladder.
- *The tests to extend:*
  - `BassSweepTests.cs`: `EnumerateF`, the 2,263-state fixture on `Data()` with three terms at most; `ConstantLikelihood_SamplesThePriorOnF`; `Schedules_AtTheBoundaries`; `Design_FollowsTheChain`; `BatchError`;
  - `BassProposalTests.cs`, `BassPriorSamplerTests.cs`, `BassTemperingTests.cs`;
  - `tests/samplers/TemperingLadderTests.cs` and `ParallelTemperingResizeTests.cs`;
  - Phase 0's and the pricing's harness, `tests/regression/BassMixing.cs`, whose `Replaced`, `With`, `Without` and `Stepped` scorers draw the moves as the pricing read them.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-moves\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-moves\` when you finish.
- **Your files.**
  - In `src/maths/regression/spline/bass/`: `BassSweep.cs`, `BassTermProposal.cs`, `BassConfig.cs`, `BassResolved.cs`, `BassRun.cs`, `BassDiagnostics.cs`, `BassLadder.cs`, and new files for the moves if they read better apart. `BassStructurePrior.cs` and `BassDesign.cs` for accessors only; name each in the report.
  - `src/maths/samplers/rjmcmc/TemperingLadder.cs`, generic. `ParallelTempering.cs` only if 4b cannot be placed without it; name the change.
  - The tests above, new `tests/regression/Bass*Tests.cs` files, and additions to `BassMixing.cs` that leave Phase 0's and the pricing's readings unchanged.
  - This brief's Reports section.
  - Leave everything else alone.
- **Commit your own files only, by explicit path,** one commit per step or a coherent pair, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked files under `issues/regression/` (`argon-bass-qa.md`, `argon-bass-question-suggestions.md`, `bass-questions.md`) stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.
- **Every option is off by default, and off means bit-identical.** With an option off, the sweep draws no extra random numbers and takes no new branch that changes arithmetic.

## The steps

Each kernel's exactness checks (the next section) land with its step. Each step ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next step. Name the passes and their fixes in the report.

**Step 0. The bit-identity baseline**, before any production edit.
- Hash the retained draws (M, σ², τ, λ, and the terms) of a set of fits at the current `HEAD`: S10's Friedman, mixed and functional fixtures untempered, on both move schedules; Friedman tempered with a self-sizing ladder; and one tempered BARS ensemble fit. Short runs are enough.
- Keep it as a test, or as a harness with the hashes recorded, so every later commit rechecks it with the options off.

**Step 1. The configuration.**
- Each option is a nullable record on its natural owner, `null` meaning off: the local steps (ρ = ½, ω = 0.02) and the fixed-M family (its weights) on `BassSampling`, and the split rule (threshold 0.9) on `TemperingSpec`.
  - Names follow the methods and stay true across each option's whole range (CONTRIBUTING's naming).
  - The record's defaults are the decided values.
- Validated in `BassConfig`'s checks with the existing message style; echoed by the resolved configuration; round-tripped through `ToJson`.

**Step 2. Lever 2: local steps inside a change.** As the lever's section defines them, with the decided values.
- *A hinge factor* (continuous or functional), with probability ρ: j′ = j + δ on the variable's distinct-value index, sign kept. |δ| is uniform on {1, …, w} with w = max(1, ⌈ω·n_v⌉), and its sign ±1 with ½ each. A step outside {0, …, n_v − 1} is a null move.
- *A categorical factor,* with probability ρ: toggle one level, uniform among the n_v. An empty or full subset is a null move. The log ratio is the subset measure's, ln π_v(D′) − ln π_v(D), taken from the prior's own `LogSubset`, under both subset measures.
- *With probability 1 − ρ,* today's redraw, unchanged.
- *A null move* is an attempted, rejected move that leaves the state and the design untouched, as the change's same-state path does today. Duplicate terms and the support rule are refused by the design, as for any change.

**Step 3. The fixed-M family.** As "Phase 1's decisions" states it.
- *Their slot.* With the family on, the change's slot splits in thirds: a factor's redraw or local step (lever 2 inside), a prune or grow, a variable change. The schedule's birth and death weights stay as they are on both schedules. Under the prior-cancelling schedule at M = 0, every fixed-M move is a null move.
- *Prune or grow.* A term uniformly, prune or grow with ½ each, then a part uniformly among the fit's P parts.
  - A prune takes a factor uniformly among the part's K. From an empty part it is a null move; a prune that would empty the whole term is refused.
  - A grow takes its variable uniformly among the part's unused variables, and its sign and knot, or its subset, from the prior. Into a full part (no unused variable, or the part's order at its cap) it is a null move.
- *A variable change.* A term uniformly, a part uniformly, a factor uniformly among the part's K, then v′ uniformly among the part's unused variables, with its sign and knot, or its subset, from v′'s prior. An empty part, or one with no unused variable, gives a null move.
- *The ratio* is the exact log ratio of the term prior and of the proposals, both ways, plus Δℓ.
  - Compute it from the prior (`LogPrior(term)`) and the proposal probabilities, not from the expected cancellation.
  - Where the draft derives the cancellation (a prune or grow that leaves K ≥ 1 in the part, and the variable change), assert it in a test.
  - A part's K = 1 → 0 depends on the part's order law at 0. Derive it and state it in the report.

**Step 4. 4b, the split rule,** in `TemperingLadder`, generic, and called from `BassLadder`'s sizing when the option is on.
- *The rule:* each adaptation round, a pair whose rejection exceeds the threshold gains a rung halfway along the pair's barrier, under `Place`'s own interpolation, up to `MaxLevels`.
- *How it combines with `Size` and `Place`* is yours to design. One reading: while any pair exceeds the threshold, a round splits rather than re-placing, and once none does, sizing resumes on the finer, less censored readings. Keep the ladder's contract: it settles a round before burn-in ends, so the last round equalizes it. State the design and its reason in the report.
- *Checks:*
  - on synthetic rejection profiles with a known steep band: the split resolves it, never exceeds `MaxLevels`, and leaves a ladder with no steep pair unchanged;
  - BARS's and every option-off BASS ladder stay bit-identical (Step 0).

**Step 5. Diagnostics.** Per move kind: attempted, accepted, null and refused, for the local step, the toggle, the prune, the grow and the variable change, beside today's birth, death and change. The pass reads them.

## Exactness

Ground truth comes from outside the sampler: enumeration of F, and the moves' definitions as the decisions state them.

- *A constant likelihood.* With each new option on alone and all together, on both schedules, the chain samples the prior on F: the count law and the order law against `EnumerateF`, as `ConstantLikelihood_SamplesThePriorOnF` does.
- *Detailed balance with the likelihood, enumerated,* at a fixed τ (`FixedG`), on the 2,263-state fixture.
  - π(S | τ, y) comes from enumeration: the log structure prior plus the collapsed log marginal, normalized over F.
  - The test builds each new kernel's transition matrix K(S, S′). It enumerates every proposal from its own statement of the selection rules (the term, part, factor, variable, sign, knot or δ, with their probabilities), and takes the acceptance from the production code's log ratio plus Δℓ. Internal accessors that expose a proposal's log ratio are allowed; name them.
  - Assert π(S)K(S, S′) = π(S′)K(S′, S) for every pair to rounding (state the tolerance; 10⁻⁹ relative at most), and πK = π for each kernel and for the whole sweep's structure kernel with births and deaths.
  - Repeat at one tempered b (0.3) against the tempered target, under the default path and form, so the ratio's prior and proposal part stays untempered and only the likelihood takes b.
- *The other parts.* The 2,263-state fixture is continuous only. Build small enumerable fixtures holding a categorical part and a functional part, and run the same detailed-balance check there. Run it under both subset measures. Report each fixture's |F|.
- *Visit frequencies.* The real chain at fixed τ, all options on, against the exact posterior on the 2,263-state fixture, within batch-means Monte Carlo error. This exercises the code path the enumeration does not.
- *Bit-identity* with the options off, against Step 0, after every commit.
- *Suites:* Regression, Samplers, and the BASS oracle tests (`tests/oracle`, the `Bass` classes), all passing.

## A smoke reading

This is exploration, not the pass. One master seed, four chains, S10's three fixtures untempered, at S10's run lengths. Run each option alone and all together, plus Friedman tempered with 4b.
- Read per move kind the acceptance, each chain's band of M, the time per sweep, and on the tempered run the settled count, the steep pairs and the ideal round-trip rate.
- Compare with M3's step acceptances and with Phase 0's Δb = 0.01 profile (Σr 4.7–5.2 over [0.10, 0.40]).
- Read, do not judge: recommend no default.

## Review before reporting

Take a fresh-eyes pass over the whole change and fix what it finds:
- *Correctness:*
  - every option off is bit-identical, with no extra random draw;
  - a null move leaves the state, the design and the tallies consistent;
  - the design follows the chain through the new moves' accepts and rejects (as `Design_FollowsTheChain` does);
  - the tempered ratio takes b on Δℓ only;
  - the split rule respects `MaxLevels` and the settle-before-burn-in contract;
  - the configuration validates, resolves and round-trips.
- *Coverage:*
  - each part kind, and both subset measures;
  - a prune that empties a part while the term keeps another, and the grow back into it;
  - a grow at the order cap;
  - the knot index's edges, and w = 1;
  - both schedules, and the prior-cancelling schedule at M = 0;
  - duplicate-term and support refusals.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the configuration's shape;
- each kernel's ratio as built, with the K = 1 → 0 derivation;
- 4b's design and its reason;
- the exactness checks, their tolerances and the fixtures' |F|;
- the bit-identity hashes;
- the smoke reading, kept apart from what it shows;
- the run times, and the suites with their counts;
- the review passes and their fixes;
- anything the main session must know, such as a choice this brief left open.

## Reports

Newest first.
