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

### Phase 1 options built: lever 2, the fixed-M family, 4b — 2026-10-29

*Landed,* on `main`, nothing pushed:
- `cf431c0`: Step 0, the bit-identity baseline (`BassBitIdentityTests`).
- `7a08c37`: Step 1, the configuration (`BassMoveOptionsTests`).
- `afafe03`: Step 2, lever 2 (`BassKernelEnumeration`, `BassLocalStepTests`, rows in `BassSweepTests`).
- `368d7c3`: Step 3, the fixed-M family (`BassFactorMoveTests`).
- `a99c02f`: Step 4, 4b in `TemperingLadder` (tests in `TemperingLadderTests`, `ParallelTemperingResizeTests`, `BassTemperingTests`).
- `c167af3`: Step 5, the moves by kind in the diagnostics.
- `aaa8034`: the whole-change review's fixes. `856ecfd`: the smoke harness, `BassLocalMoveSmokeTests`, behind `BASS_QUALIFICATION` = 1. `8217eba`: a pinned JSON shape the full suite caught.
- Production files: `BassConfig`, `BassResolved`, `BassTermProposal`, `BassSweep`, `BassLadder`, `BassDiagnostics`, `BassRun`, and `TemperingLadder`. Unchanged: `BassStructurePrior`, `BassDesign`, `BassModel` and `ParallelTempering`. `BassMixing.cs` is unchanged too, so Phase 0's and the pricing's readings stand.

*The configuration.* Each option is a nullable record, null meaning off; the records' defaults are the decided values.
- `BassSampling.LocalSteps`: `LocalSteps(Probability = 0.5, WindowFraction = 0.02)`. `Window(n)` = max(1, ⌈ω·n⌉), with ω·n taken in decimal so that 0.02·150 stays 3.
- `BassSampling.FactorMoves`: `FactorMoves(Change = 1, PruneOrGrow = 1, VariableChange = 1)`, relative weights within the change's slot, thirds by default. The prune and the grow take half of their weight each.
- `TemperingSpec.Split`: `SplitRule(Threshold = 0.9)`, refused unless the ladder sizes itself (`Levels` null).
- Validation names each field (`Sampling.LocalSteps.Probability`, …). `Change` may be 0, a family without the change, except with the local steps on, since they run inside the change.
- The resolved echo: `BassResolved.StepWindows`, each continuous and functional variable's w from its knots; `FactorShares`, the slot's shares; `SplitThreshold`. With the steps on, `Resolve` prepares the inputs once to count each variable's knots as the prior does, which leaves `BassModel` untouched. The JSON round-trips, and the fingerprint moves with each option.
- The palette with the family on: birth, death, change, then `prune`, `grow` and `variable`. Each new move's weight is the change slot times its share, so the births' and deaths' weights are the schedule's own. The prune and the grow reverse each other at equal weights, and the variable change reverses itself.

*Each kernel's ratio as built.* The within-move ratio is ln q(reverse) − ln q(forward), from the proposals both ways. The target adds the term prior's ratio through `LogPrior`. Every pick named below is the same both ways, and the term's 1/M is matched by the selection.
- *The local step:* 0. |δ| is 1/w and its sign ½ both ways, and the factor pick is the same. The sign is kept and the knot prior is uniform, so the acceptance is Δℓ.
- *The toggle:* 0. The level is 1/n_v both ways. The target's prior ratio carries the measure's ln π_v(D′) − ln π_v(D): ln C(n, |D|) − ln C(n, |D′|) under R's measure, 0 under the paper's.
- *The prune*, factor f of a part with K factors:
  - the forward pick: −ln P − ln K;
  - the grow back: −ln P − ln(p − K + 1) + ln q_v(f);
  - q_v(f) is 1/(2n_v) for a hinge and π_v(D) for a subset, the law the grow draws from.
- *The grow*, into a part with K factors:
  - the forward pick: −ln P − ln(p − K) + ln q_v(f);
  - the prune back: −ln P − ln(K + 1).
- *The variable change:* ln q_v(old) − ln q_{v′}(new). The part, the factor and the p − K unused variables are the same both ways.
- *The cancellation.* Over every term and every move of each fixture, the ratio plus the prior's ratio is 0 to within 1.8·10⁻¹⁵, against a tolerance of 10⁻¹²:
  - continuous: 72 prunes with K ≥ 1 left, 72 grows, 72 variable changes;
  - categorical, both measures: 192 prunes that empty a part, 192 grows into an empty part, 392 variable changes;
  - functional: 120, 120 and 336.
- *A part's K = 1 → 0, derived.* With several parts, each part's order is uniform on {0, …, cap} (D01, decided 2026-10-28), so the order term ln 1/(cap + 1) is the same at K = 1 and K = 0.
  - π(t) holds the part's set term −ln C(p, 1) = −ln p and the factor's ln π_v(f); π(t′) holds neither.
  - So ln π(t′) − ln π(t) = ln p − ln π_v(f).
  - The proposals: the prune is (1/P)(1/1). The grow back into the empty part is (1/P)(1/p)·π_v(f), all p variables being unused. Their log ratio is −ln p + ln π_v(f).
  - The sum is 0, so this case cancels as well. With one part, K = 1 → 0 empties the term, which is the intercept, outside F, and is refused.
- *Null moves and refusals.* A null move proposes the state itself with a ratio of −∞, so the chain counts an attempt, draws its uniform and rejects; the state and the design are untouched.
  - Null moves: a step off the knots; a toggle to an empty or full subset; a prune from an empty part; a grow into a part at its cap or with no unused variable; a variable change with none unused; and the whole family at M = 0, which only the prior-cancelling schedule gives weight.
  - A prune of a term's only factor is refused by the move, with no `BassRefusal` reason.
  - Duplicates and support are refused by the design, as for any change.
  - The change at M = 0 stays unavailable, as today.

*4b's design.* `TemperingLadder.Split(ladder, rejection, threshold, maximum)` is generic.
- For each pair rejecting past the threshold it adds a rung halfway along the pair's barrier, on `Place`'s monotone cubic through the cumulative reading, inverted at the pair's midpoint level.
- Every current rung is kept. The steepest pairs split first within the maximum, the colder on a tie.
- `BassLadder.Resize` resizes to the split ladder when the rule is on and the split adds a rung. Otherwise it sizes as before: `Size` then `Place`, or `AdaptLadder` at an unchanged count.
- `Sources` maps each kept rung to its own replica (at a barrier distance of 0), so `ParallelTempering.Resize` keeps every replica untouched. Each split rung starts from a neighbour of its pair. `ParallelTempering.cs` is unchanged.
- The run's last round adapts and never resizes (`BassRun.Burn`, unchanged), so a split ladder is equalized before it freezes: the contract holds.
- *The reason:* a censored pair cannot say how much barrier it hides, so a split measures the band instead of extrapolating an under-read barrier through the cubic. Keeping every rung keeps the other pairs' readings and replicas intact. Once no pair is steep, `Size` and `Place` work from finer, less censored readings.

*The exactness checks.*
- *Detailed balance by enumeration* at fixed τ = 1/N.
  - π comes from the prior's own `LogCollectionWeight` and `LogPrior` plus the collapsed log marginal.
  - K(S, S′) is built from the moves' statement, with the production palette's weights and the production ratios.
  - Tolerance: 10⁻⁹ in logs, which is relative. The worst gap was 5.7·10⁻¹⁴, for pairwise balance and for πK = π per kernel and over the whole sweep with births and deaths.
  - Both schedules, b = 1 and b = 0.3 on the default path and form, and the steps alone, the family alone and both together.

  | fixture | \|F\| | births/deaths | redraws | local steps | prunes/grows | variable changes |
  | --- | --- | --- | --- | --- | --- | --- |
  | `Data()`, continuous, M ≤ 3 | 2,263 | 12,928 | 29,600 | 9,792 | 15,248 | 7,208 |
  | categorical (x0; z1 of 2 levels, z2 of 4), caps 1, M ≤ 2, either measure | 2,918 | 11,516 | 67,156 | 18,200 | 16,624 | 18,488 |
  | functional (x0; outputs t1 and t2 over 4 locations), caps 1, M ≤ 2 | 593 | 2,300 | 5,066 | 1,576 | 3,112 | 2,616 |

  The counts are transitions S → S′ ≠ S. The continuous fixture also ran with ω = 0.7 (w = 3 over 3 knots, ±3 always off) and the functional one with ω = 0.6.
- *Negative controls* fail the check:
  - a step ratio 0.05 off;
  - a toggle missing the subset measure's prior ratio under R's measure;
  - a part-emptying prune and its grow scored with an extra ln 2.
- *The draws against the statement.* `ProposeChange` (w = 3, ρ = 0.4, 400,000 draws), toggles under both measures (300,000), and the family's prunes, grows and variable changes. Every outcome is within 5 standard errors of its stated probability.
- *Constant likelihood,* on both schedules with the steps, the family, and all options on: the count law and the order law against `EnumerateF`, π(1)/π(0), τ's prior and λ (`BassSweepTests`).
- *Visit frequencies.* The real chain at fixed τ, 4 chains × 150,000 sweeps, with the steps, the family, and all options on. The count law and the ten likeliest structures are each within 5 batch-means errors + 0.002. All options on: M = 1 0.8145 against 0.8153 (± 0.0020), M = 2 0.1696 against 0.1688 (± 0.0017), M = 3 0.0159 against 0.0160 (± 0.0005).
- *The design follows the chain* through every new move's accepts and rejects. The checks are `Design_FollowsTheChain` with all options on, scalar and functional, and every ladder replica's design (`Ladder_EveryReplicasDesignFollowsItsState` with all options on and the reference rung).
- *4b on synthetic profiles:*
  - exact midpoints on a linear barrier;
  - the steepest first, within the maximum, colder on ties;
  - no steep pair leaves the ladder bit for bit;
  - a known steep band (a step of 6 in Λ near b = 0.2, censored as r = 1 − e^(−ΔΛ)) is resolved in 4 rounds to 15 rungs, 5 of them in the band, the largest rejection 0.83. At a maximum of 12 it stops there.
- *4b on BASS,* on the N = 175 functional fixture where S8's 8 rungs rejected 0.96–1.00, every kept replica stayed the same object with its design:
  - at a maximum of 24, 5–6 split rounds and 1–2 sizing rounds, settling at 14 and 16 rungs;
  - at a maximum of 11, 2 of each, settling at 11.
- *The tallies.* The kinds of each palette move add up exactly to the chain's own attempts and acceptances, chain by chain: untempered, and at a splitting ladder's b = 1, on all three fixtures. Refusals by kind equal the reasons' tally once prunes are off. A duplicate term (variable change) and a support refusal (local step) are refused as their kind, the design untouched. A resumed run keeps the tallies.

*The bit-identity hashes* (`BassBitIdentityTests`) were recorded at `cf431c0`, before any production edit. Each covers the retained draws, the terms, each chain's move counts, refusals and empty births, and the frozen rungs. Every step's commit rechecked them with the options off, beside the older pins (`BassTemperingSizingTests.FixedLevels_RunAsS8Did`, `ParallelTemperingResizeTests.NeverResized_IsBitIdenticalToBefore`), which also hold.

| fit (2 chains, 300 burn-in, 100 draws) | hash |
| --- | --- |
| Friedman, equal / prior-cancelling | `ee7d60a2738b56d8` / `414e30ea867f9058` |
| mixed, equal / prior-cancelling | `137dad3cbeb33acb` / `72fda364ee08f1f6` |
| functional, equal / prior-cancelling | `3e5f43ed894f60bc` / `e57e8cbe0c39d339` |
| Friedman tempered, self-sizing (255 burn-in, 50 draws; settled at 14 and 16 rungs) | `0d067c924599ed93` |
| BARS, tempered, 4 levels (fit, sd, ladder readout) | `7edb59818f141da9` |

*The smoke reading.* One master seed per fixture (S10's: 53, 57, 61), four chains, 5,000 + 2,000 sweeps, each fit on one thread. The options-off rows are S10's chains. α̃ is the acceptance of admitted proposals, those neither null nor refused.

| fixture, options | µs/sweep | chains' M (median, [band]) | M's R̂ | local step α̃ (null) | toggle α̃ (null) | prune α̃ (null, refused) | grow α̃ (null, refused) | variable change α̃ (null, refused) | redraw α̃ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Friedman, off | 55 | 10 [10,12], 8 [8,10], 11 [10,12], 9 [9,11] | 2.35 | | | | | | 0.056 |
| Friedman, steps | 57 | 9 [9,12], 10 [10,13], 12 [11,13], 11 [10,14] | 2.02 | 0.830 (0.03) | | | | | 0.072 |
| Friedman, family | 62 | 10 [10,12], 10 [10,13], 9 [9,11], 9 [9,11] | 1.74 | | | 0.007 (0, 0.69) | 0.009 (0.01, 0.09) | 0.009 (0, 0.05) | 0.056 |
| Friedman, all | 51 | 8 [8,10], 10 [10,12], 9 [9,10], 11 [11,13] | 2.45 | 0.806 (0.05) | | 0.029 (0, 0.71) | 0.012 (0.01, 0.12) | 0.005 (0, 0.04) | 0.044 |
| mixed, off | 44 | 6 [5,7], 6 [6,7], 6 [6,7], 5 [5,7] | 1.54 | | | | | | 0.153 |
| mixed, steps | 39 | 6 [6,7], 7 [7,10], 5 [5,7], 5 [5,6] | 2.93 | 0.835 (0.03) | 0.008 (0.29) | | | | 0.167 |
| mixed, family | 35 | 5 [5,6], 7 [7,8], 5 [5,6], 6 [6,8] | 3.42 | | | 0.010 (0.11, 0.22) | 0.005 (0.39, 0.15) | 0 (0.51, 0.07) | 0.120 |
| mixed, all | 35 | 10 [10,11], 9 [9,11], 5 [5,6], 6 [6,7] | 3.81 | 0.834 (0.01) | 0.078 (0.36) | 0 (0.12, 0.38) | 0.011 (0.43, 0.15) | 0.003 (0.56, 0.06) | 0.228 |
| functional, off | 42 | 5 [5,6], 7 [7,8], 5 [5,6], 7 [7,9] | 2.55 | | | | | | 0.057 |
| functional, steps | 37 | 5 [5,6], 5 [5,6], 5 [5,6], 5 [5,7] | 1.02 | 0.403 (0.10) | | | | | 0.037 |
| functional, family | 37 | 5 [5,5], 5 [5,6], 6 [5,8], 5 [5,7] | 1.66 | | | 0.005 (0.21, 0.28) | 0.011 (0.30, 0.10) | 0 (0.55, 0.04) | 0.048 |
| functional, all | 46 | 8 [6,9], 7 [6,8], 5 [5,5], 6 [6,7] | 3.85 | 0.392 (0.07) | | 0.013 (0.22, 0.25) | 0.008 (0.28, 0.13) | 0.003 (0.53, 0.06) | 0.042 |

Friedman tempered, self-sizing, without and with 4b:

| | µs/scan | rungs per chain | steep pairs (> 0.9) per chain | ideal round trips/scan | observed | rungs in [0.10, 0.40]; Σr there | M's R̂ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| without | 1,072 (50 per replica sweep) | 20, 22, 20, 23 | 2, 1, 1, 1 | 0.0014 | 0 | 6, 9, 8, 9; 3.41, 4.09, 3.45, 3.64 | 1.45 |
| with 4b | 1,003 (51 per replica sweep) | 20, 21, 19, 19 | 0, 1, 2, 2 | 0.0091 | 0 | 7, 9, 6, 8; 3.20, 3.19, 3.02, 3.62 | 1.49 |

A probe of 4b's rounds on Friedman's chain 0, run once and not committed:
- the steep pair moved between rounds as the replicas' states moved, from (1.0, 0.25) at round 8 to (0.60, 0.37) at round 1,024;
- after the last sizing round split it, the equalizing round read no steep pair, and the sampling phase none (Σr 9.55 over the ladder).
- The other chains' steep pairs in sampling are therefore pairs that steepened after burn-in, not splits that the equalization undid.

*What the smoke reading shows,* read, not judged, and recommending no default:
- *The local steps accept as M3 predicted on Friedman and the mixed fit:* 0.83 admitted, pooled over |δ| ≤ 10 and ≤ 6, against M3's 0.93–0.97 for steps of 1–3 and 0.6–0.8 for 4–15.
  - On the functional fit they accept 0.39–0.40: w is 1–2 there, and one step is 2–5% of a variable's range.
  - Toggles accept 0.008–0.078, and 29–36% of them are null, on z's three levels.
- *The family accepts rarely*, at most 0.03 of admitted proposals, as P2's prices predicted.
  - On the categorical and functional fits 51–56% of its variable changes are null, parts with no unused variable.
  - On Friedman 69–71% of its prunes are refused. Every prune of a main effect is (it would leave the intercept), beside the design's refusals.
- *M's R̂ and the chains' bands move both ways* with one seed:
  - with the steps, the functional fit's four chains share M = 5 (R̂ 1.02, against 2.55 off), while the mixed fit's R̂ rises from 1.54 to 2.93;
  - with all options on, two mixed chains sit at M ≈ 10, a mode the off run did not visit.
  - These are single draws of chain behaviour; the pass reads five seeds.
- *The time per sweep* stays within −20% to +13% of the off run on one thread, where M differs too; the family's moves cost a change's design update.
- *4b settles at about the same count* (19–21 rungs against 20–23), with as many steep pairs in all (5 against 5), though chain 0's two are gone.
  - Its ideal round-trip rate is 6.5 times higher (0.0091 against 0.0014 per scan).
  - Neither ladder completes a round trip in 2,000 scans.
  - Σr over [0.10, 0.40] reads 3.0–3.6 with 4b and 3.4–4.1 without, both below Phase 0's Δb = 0.01 profile (4.7–5.2) over 6–9 rungs.

*Run times.* The step tests: `BassLocalStepTests` about 15 s, `BassFactorMoveTests` about 22 s, `BassMoveOptionsTests` about 10 s, the split tests about 5 s. The smoke reading: 1.2 minutes for the 12 untempered fits and 58 s for the two tempered ones.

*The suites,* at `8217eba`, all passing: Regression 1,030 (3 min 42 s), Samplers 194 (1 min 11 s), and the BASS oracle classes (`tests/oracle`, `Bass*`) 54 (2 min 17 s).

*The review passes and their fixes:*
- *Step 0.* The hash covered the draws only. It now covers each chain's move counts, refusals and empty births too, re-recorded before any production edit.
- *Step 1.* None found.
- *Step 2.* The enumeration had no negative control, so two were added. `Design_FollowsTheChain` did not run with the steps on; rows were added.
- *Step 3.* Two test-side bugs, fixed:
  - the enumeration built a variable change by removing the term's only factor first, which `BassTerm` refuses;
  - the ratio test's skip of a refused prune also skipped that factor's variable changes.
  - The direct duplicate and support refusal test went to Step 5's tallies.
- *Step 4.* None found.
- *Step 5.* The null-move test's palette held a grow without its prune, which the chain refuses; the class summary was stale.
- *The whole change.* The proposal's and the tempering spec's summaries now name the new options, and the ladder's design-following check runs with every new move on.
- *The full suite* caught one pinned JSON shape, `BassTemperingTests.Spec_HasD05sDefaults_RoundTrips_AndMovesTheFingerprint`: a tempering spec now serializes `"Split":null`. Updated in `8217eba`; the filtered runs of Steps 1–5 had not included it.

*For the main session:*
- *Choices the brief left open,* each stated in the code:
  - the null move as a ratio of −∞, attempted and rejected, which draws the chain's uniform, with the options on only;
  - the move keys `prune`, `grow` and `variable`;
  - `FactorMoves.Change` = 0 allowed without the local steps;
  - the fixed-M part draw takes one word even when P = 1;
  - 4b splits in place of sizing while a pair is steep.
- *Accessors,* all internal on `BassTermProposal`: `Stepped`, `Toggled`, `LogRedrawRatio`, `Pruned`, `Grown` (hinge, subset), `Exchanged` (hinge, subset), `LogPruneRatio`, `LogGrowRatio` and `LogVariableChangeRatio`. The sweep computes its ratios through the same methods. No accessor was added to `BassStructurePrior` or `BassDesign`; the enumeration reads their existing internals.
- *Shapes that changed:*
  - `BassDiagnostics` and `BassChainAccount` gained `MoveKinds` (positional, built by `BassRun` only);
  - `BassMovesSnapshot` gained `Kinds`;
  - `BassTermProposal`'s constructor gained an optional `LocalSteps`.
  - `ResetTallies()` keeps its signature, since three Phase 0 tests call it. A pending proposal's later acceptance is left out of the new phase, which keeps the kinds equal to the chain's counts.
- *D01's tilt readout* (`Support`) counts the family's candidate terms with the family on, so its refusal rates by order mix proposals then.
- *The smoke reading's R̂ and bands* differ across options in both directions on one seed. Phase 2's five seeds and its pre-registered rule decide.
