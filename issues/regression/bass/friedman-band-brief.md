# Friedman's band and modes: four readings (chip brief)

Date: 2026-10-29. The user asked for this measurement after the review notes in `issues/regression/argon-friedman-etc.md` (argon's, and the Opus chat's read of them) and the main session's response. It measures, and changes nothing in production. Its readings decide what, if anything, is built next for Friedman's worse modes:
- more particles, or an SMC schedule that keeps lineages;
- a compound replacement move;
- or nothing.

- *What is known.*
  - F-M2's Phase 0 (M4) found untempered chains in local modes 24–42 nats apart. M5 found a steep band at b ≈ 0.16–0.32 where ℓ's spread peaks and M falls by about two terms. At Δb = 0.01, no pair in [0.10, 0.40] rejects above 0.46.
  - Stage 8's C4 ran fitting SMC at 1,000 particles only, varying the data size. On Friedman at S10's size, the replicates disagree: the evidence's between-replicate error is 4.2 nats, and each replicate ends descended from one or two lineages.
  - D10's claim that the worse modes carry about e⁻¹² of the mass compared densities, not masses (corrected 2026-10-29).
- *Not in scope:* any production change, any default, and building any move. Recommend nothing; the readings go to the user's call.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- `issues/regression/argon-friedman-etc.md`, whole;
- F-M2's Phase 0 report ([local-move-measurements-brief.md](local-move-measurements-brief.md#reports)), M4 and M5;
- C4's report ([stage8-fitting-smc-brief.md](stage8-fitting-smc-brief.md#reports)).

In the code:
- *Phase 0's harness,* `tests/regression/BassMixing.cs`, with `BassMixingModeTests` (M4) and `BassMixingBandTests` (M5's dense fixed ladder).
- *C4's qualification,* `tests/regression/BassSmcQualificationTests.cs`.
- *BASS's refusal tallies by reason,* `BassDiagnostics` and `BassRefusal`.
- `BassStructureSpec.AngleFloor`, 10⁻¹⁰ by default (D06), and the support rule.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-friedman\` to every `dotnet build` and `dotnet test`. Build in Release for the runs. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-friedman\` when you finish.
- **Your files:**
  - a new harness, `tests/regression/BassFriedmanBandTests.cs`, behind `BASS_QUALIFICATION` = 1;
  - additions to `BassMixing.cs` that leave Phase 0's readings unchanged;
  - this brief's Reports section.
  - Leave production code alone. If a reading cannot be taken without a change, stop and say so.
- **Commit your own files only, by explicit path,** in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked notes under `issues/regression/` stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.
- **Read, do not judge.** Each harness prints its readings and asserts only its own soundness.

## The readings

All on S10's Friedman fixture (n = 500, p = 10, data seed 51), with S10's master seed 53 unless stated.

**A structural key,** used by R2 and R4: a structure's multiset of its terms' variable sets, ignoring signs and knots. For example, {x1·x2, x1·x2, x3, x4, x5} for a structure carrying two x1·x2 terms.
- Report each sample's basin under it. Group the rare keys (under 1% everywhere) as "other".
- Name the basins M4's best and worse chains hold.

**R1. ℓ at fixed b inside the band: coexistence or crossover.** Run M5's dense fixed ladder: Δb = 0.01 over [0.10, 0.40], its four chains, M5's lengths. At every rung, record each sample's untempered log-likelihood ℓ₁ and M.
- Report per rung:
  - the histogram of ℓ₁ (pooled over chains, and per chain);
  - its mean, sd, skewness and Sarle's bimodality coefficient;
  - a two-component Gaussian mixture against one, by BIC, with the components' weights and their separation in sd units;
  - the same for M.
- The question: at some b in the band, is ℓ₁ bimodal, two phases coexisting (a first-order transition), or unimodal and broad (a continuous crossover)?
- Also report Λ, the ladder's Σr over the band and over the whole ladder, for OQ2's record.

**R2. SMC against the particle count.** Fitting SMC at its defaults, apart from the particle count: four replicates each of 1,000, 4,000 and 16,000 particles, with S10's master seed and the replicates' seeds from the seed tree. Per run:
- each replicate's ln Ẑ with its genealogy error, and the replicates' spread (their sd and the pooled between-run error);
- each replicate's basin fractions at b = 1 under the structural key, weighted;
- each replicate's families, and the stage count;
- the mean log posterior per replicate (M4's definition, for comparison with M4);
- σ's 95% interval;
- the time.

Since ln Ẑ_basin = ln Ẑ + ln(basin fraction), also report each basin's log mass per replicate. That is the comparison of masses, not densities.

**R3. Mutation against particles.** At 1,000 particles, four replicates with `MutationSteps` 20 against C4's 5, with everything else as in R2. The same readings, plus the kernel's acceptance per move kind.

**R4. The support's share of the barrier.**
- *First, the refusals.* On M4's four untempered chains at their lengths, report each proposal kind's refusals by reason (angle, support, variance inflation, numerical, cap) as fractions of its proposals, over the sampling phase and per chain. Report the same along R1's ladder, by rung.
- *Then, only if the angle floor refuses more than 1% of any kind's proposals,* rerun M4's chains with the floor at 10⁻¹⁴, and then 10⁻¹⁸. Report each chain's log-posterior median, M's R̂, the log posterior's R̂, and the basin each chain holds.
  - A lower floor changes F, the target's support. Say so beside the readings: they measure the floor's share of the stall, not the same posterior.
- *Do the same for the support rule's minimum,* if it refuses more than 1%: half its count, the floor unchanged.

**MCMC's basins for comparison.** On M4's four chains, report the basin fractions under the same key, per chain. Where S10's tempered Friedman chains are cheap to rerun, give theirs as well.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - ℓ₁ is the untempered log-likelihood at every rung;
  - the structural key ignores knots and signs;
  - the basins' log masses add up to ln Ẑ within rounding;
  - the weights are normalized;
  - Phase 0's readings are unchanged where its helpers were touched.
- *Coverage:*
  - a reading that rests on one chain or one replicate;
  - the bimodality test's sensitivity at the band's edges;
  - any run that hits its stage budget.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- each reading's tables, and what they show, kept apart;
- the run times;
- the review passes and their fixes;
- anything the main session must know.

Recommend no route.

## Reports

Newest first.
