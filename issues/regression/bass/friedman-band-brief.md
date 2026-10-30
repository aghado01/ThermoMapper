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

### R1–R4 measured — 2026-10-29

*Landed* (tests only: no production change, no internal accessor, so the defaults are untouched by construction):
- `5673978`: the harness `tests/regression/BassFriedmanBandTests.cs`, behind `BASS_QUALIFICATION` = 1: R1 `Band_LogLikelihoodAtFixedB`, R2 and R3 `Smc_AgainstParticlesAndMutation`, R4 `Refusals_OnM4sChains`. One addition to `BassMixing.cs`: `Fixture.Basin`, the structural key.
- Each test reads and asserts only its own soundness (the review passes below).

*The basins.* The key is the brief's: the multiset of the terms' variable sets, signs and knots ignored. The harness labels the basins M4's chains hold (1% or more of a chain's states) A, B, …, the best chain's first. The lower-case labels are this report's, for basins M4's chains hold under 1%:

| label | basin | M | held by |
| --- | --- | --- | --- |
| A | x1 + 2 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 8 | M4's best chain (1), 79.6% |
| B | x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 9 | M4's chain 3, 82.7%; chain 1, 11.7% |
| G | x1 + 5 × x1·x2 + 2 × x3 + x4 + x5 | 10 | M4's chain 0, 75.0% |
| S | 2 × x1 + 3 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5 | 11 | M4's chain 2, 61.3% |
| T | x1 + 3 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5 | 10 | M4's chain 2, 30.7% |
| AC | B + x6·x8 | 10 | M4's chain 3, 1.3% |
| a1 | 4 × x1·x2 + 2 × x3 + x4 + x5 | 8 | SMC 16,000, replicate 0 |
| a2 | 4 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 9 | SMC 16,000, replicate 2 |
| a3 | 3 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 8 | SMC 16,000, replicate 2 |
| a4 | x1 + 3 × x1·x2 + x2 + 3 × x3 + x4 + x5 | 10 | SMC, several replicates |
| a5 | x1 + 2 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5 | 9 | SMC, several replicates |
| a6 | 2 × x1 + 2 × x1·x2 + 3 × x2 + 2 × x3 + x4 + x5 | 11 | SMC 4,000, replicate 2 |
| a7 | 3 × x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 11 | SMC 4,000, replicate 3 |
| a8 | x1 + 4 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 10 | SMC 4,000, replicate 1 |
| a9 | 2 × x1 + 3 × x1·x2 + 2 × x1·x2·x8 + x2 + 2 × x3 + x4 + x5 | 12 | S10's tempered chain 2, 63.4% |
| a10 | 2 × x1 + 2 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5 | 10 | S10's tempered chain 3, 27.5%; R4's rerun, chain 3, 51.0% |
| a11 | 2 × x1 + 3 × x1·x2 + 2 × x1·x2·x5 + 2 × x2 + 2 × x3 + x4 + x5 | 13 | R4's rerun, chain 2, 87.0% |
| a12 | 2 × x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5 | 10 | S10's tempered chain 0, 6.7%; R4's rerun, chain 1, 12.2% |

The key is fine-grained: M4's chains alone hold 29 basins at 1% or more, and 28 more reach 5% in some sample of the tempered chains, the rerun or SMC. The harness prints every basin at 0.1% or more of each sample; the tables here keep the largest. Beside each it prints a coarser view, the share by the count of x1·x2 terms, the count that sets M4's modes apart.

#### R1. ℓ₁ at fixed b

M5's dense ladders, rebuilt as M5 builds them: Part 1 is S10's self-sized ladder per chain; Part 2 takes its rungs above 0.40, 31 rungs from 0.40 to 0.10 at Δb = 0.01, and its rungs below 0.10. Four chains, 5,000 scans of burn-in, 2,000 sampled. ℓ₁ is each replica's untempered log marginal after the scan's swaps. The rows reproduce M5's (chain 0 at b = 0.40: M 8.29, ℓ −805.2 ± 8.6; at 0.24: ℓ −906.6 ± 60.9).

*Λ̂ = Σr* (a pair straddling an edge counted by its share of Δb inside):

| chain | Part 1 rungs | band [0.16, 0.32] | window [0.10, 0.40] | ladder | Part 2 rungs | band | window | ladder | Part 2's largest pair in the window |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 20 | 2.31 | 4.06 | 8.64 | 45 | 2.99 | 4.98 | 9.87 | [0.23, 0.24] 0.341 |
| 1 | 22 | 2.41 | 4.34 | 8.93 | 44 | 2.99 | 5.15 | 10.17 | [0.12, 0.13] 0.282 |
| 2 | 20 | 2.04 | 4.06 | 8.98 | 43 | 2.69 | 4.82 | 9.59 | [0.24, 0.25] 0.456 |
| 3 | 23 | 2.10 | 3.99 | 8.64 | 45 | 2.79 | 4.68 | 9.20 | [0.12, 0.13] 0.359 |

No ladder completed a round trip in its 2,000 scans; Part 2's ideal rate is 0.030–0.034 per scan (M5's).

*ℓ₁ at the window's rungs, the four chains pooled* (8,000 states):
- BC is Sarle's bimodality coefficient, above 5/9 for a bimodal or strongly skewed sample.
- D is Ashman's separation of the two-Gaussian fit's components, above 2 for a clean split; "means" are theirs, and "minority w" the smaller one's weight.
- ESS is the four chains' for the mean.
- Each histogram spans the rung's pooled range in 24 bins, `·` marking an empty one.

| b | ℓ₁ mean | sd | skew | BC | D | minority w | means | ESS | histogram |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.40 | −809.6 | 10.2 | −1.07 | 0.473 | 1.51 | 0.312 | −818.9, −805.3 | 42 | `▁▁▁▁▁▁▁▁▁▁▂▂▂▃▄▆▇███▇▃▁▁` |
| 0.39 | −810.7 | 10.9 | −1.02 | 0.480 | 1.55 | 0.326 | −820.5, −806.0 | 41 | `▁▁▁▁▁▁▁▁▁▁▂▂▃▄▅▆▇█▇█▇▃▁▁` |
| 0.38 | −811.9 | 11.5 | −0.96 | 0.484 | 1.58 | 0.351 | −821.9, −806.5 | 39 | `▁▁▁▁▁▁▁▁▂▂▂▂▃▄▆▇▇█▇█▆▃▁▁` |
| 0.37 | −813.3 | 12.2 | −0.90 | 0.489 | 1.61 | 0.395 | −823.0, −807.0 | 35 | `▁▁▁▁▁▁▁▁▂▂▂▃▃▅▅▆▇█▇█▇▃▁▁` |
| 0.36 | −814.9 | 12.9 | −0.88 | 0.490 | 1.65 | 0.384 | −825.5, −808.3 | 33 | `▁▁▁▁▁▁▁▁▂▂▂▃▃▅▅▆▇█▇▇▆▃▁▁` |
| 0.35 | −816.5 | 13.5 | −0.79 | 0.494 | 1.73 | 0.415 | −827.2, −809.0 | 35 | `▁▁▁▁▁▁▁▂▂▃▃▄▄▅▆▇▇██▇▆▃▁▁` |
| 0.34 | −818.7 | 14.5 | −0.89 | 0.386 | 1.64 | 0.500 | −827.9, −809.5 | 39 | `▁····▁···▁▁▁▁▁▂▂▄▄▆██▇▄▁` |
| 0.33 | −820.8 | 15.1 | −0.77 | 0.384 | 1.65 | 0.458 | −829.5, −810.4 | 38 | `▁····▁···▁▁▁▁▁▂▃▅▅▇██▇▄▁` |
| 0.32 | −822.9 | 16.0 | −0.88 | 0.351 | 1.54 | 0.418 | −831.0, −811.8 | 41 | `▁▁▁··▁▁▁▁▁▁▁▁▂▂▃▆▆▇██▆▃▁` |
| 0.31 | −825.5 | 17.1 | −0.97 | 0.349 | 1.45 | 0.389 | −833.0, −813.7 | 46 | `▁▁▁▁▁▁▁▁▁▁▁▁▁▂▃▄▇▆███▆▃▁` |
| 0.30 | −829.1 | 19.8 | −1.54 | 0.383 | 1.91 | 0.035 | −880.5, −827.2 | 50 | `▁·▁▁▁▁▁▁▁▁▁▁▁▁▁▃▅▇▇██▆▃▁` |
| 0.29 | −833.5 | 22.7 | −1.63 | 0.446 | 1.91 | 0.066 | −882.8, −830.0 | 51 | `▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▃▆▇███▆▄▁` |
| 0.28 | −840.3 | 29.8 | −1.98 | 0.570 | 2.16 | 0.106 | −901.2, −833.1 | 50 | `▁··▁▁▁▁▁▁▁▁▁▁▁▁▁▂▄▇██▇▄▁` |
| 0.27 | −851.9 | 39.8 | −1.52 | 0.652 | 2.80 | 0.188 | −919.3, −836.3 | 40 | `▁▁▁▁▁▁▁▁▁▁▂▂▁▁▂▂▃▅▇█▇▆▃▁` |
| 0.26 | −873.2 | 52.4 | −0.90 | 0.616 | 3.03 | 0.359 | −932.6, −840.0 | 20 | `▁▁▁▁▁▁▁▂▃▃▂▄▃▂▁▃▃▅▇█▇▆▃▁` |
| 0.25 | −900.3 | 61.0 | −0.41 | 0.536 | 2.90 | 0.437 | −943.8, −844.2 | 14 | `▁▁▃▂▁▂▁▃▇▇▄▇▅▃▂▃▄███▇▅▃▁` |
| 0.24 | −934.1 | 65.4 | 0.08 | 0.523 | 3.28 | 0.280 | −966.3, −851.2 | 11 | `▁▁▁▄▇▂▂▃▂▇█▅▄▃▃▂▂▅▅▄▄▂▁▁` |
| 0.23 | −963.5 | 60.1 | 0.58 | 0.549 | 2.27 | 0.252 | −1028.7, −941.6 | 11 | `▁▁▁▄█▄▂▃▂▄▇▅▂▂▂▁▂▂▂▂▂▂▁▁` |
| 0.22 | −986.2 | 51.3 | 0.95 | 0.564 | 2.03 | 0.384 | −1029.5, −959.2 | 12 | `▁··▁▁▁▄█▂▂▂▃▄▃▁▁▁▁▁▁▁▁▁▁` |
| 0.21 | −1001.9 | 43.5 | 1.30 | 0.607 | 1.91 | 0.469 | −1030.5, −969.4 | 14 | `▁··▁▁▂▆█▁▂▂▂▃▃▁▁▁▁▁▁▁▁▁▁` |
| 0.20 | −1014.6 | 34.7 | 1.63 | 0.622 | 1.61 | 0.347 | −1031.6, −982.6 | 20 | `▁·▁▁▁▂▅█▃▁▁▁▂▂▁▁▁▁▁▁▁▁▁▁` |
| 0.19 | −1023.2 | 27.3 | 1.71 | 0.545 | 1.11 | 0.269 | −1032.2, −998.8 | 40 | `▁▁▁▁▁▁▁▃▇█▁▁▁▁▁▁▁▁▁▁▁▁·▁` |
| 0.18 | −1029.4 | 23.1 | 1.15 | 0.229 | 0.54 | 0.216 | −1033.1, −1016.0 | 179 | `▁▁▁▁▁▁▁▁▂▇█▁▁▁▁▁▁▁▁▁▁▁·▁` |
| 0.17 | −1033.8 | 22.3 | −0.33 | 0.094 | 0.05 | 0.186 | −1034.1, −1032.3 | 443 | `▁▁▁▁▁▁▁▁▁▁▃▆█▂▁▁▁▁▁▁▁▁▁▁` |
| 0.16 | −1039.1 | 25.0 | −2.01 | 0.408 | 0.61 | 0.167 | −1058.4, −1035.3 | 591 | `▁··▁▁▁▁▁▁▁▁▁▁▂▆█▃▁▁▁▁▁▁▁` |
| 0.15 | −1048.2 | 34.5 | −1.94 | 0.710 | 1.42 | 0.214 | −1091.1, −1036.6 | 56 | `▁·▁▁▁▁▁▁▁▁▁▁▁▁▂▅█▅▁▁▁▁▁▁` |
| 0.14 | −1064.8 | 47.5 | −1.19 | 0.762 | 4.85 | 0.249 | −1141.2, −1039.3 | 31 | `▁▁▁▁▁▁▁▁▂▃▁▁▁▁▁▂▅█▆▁▁▁▁▁` |
| 0.13 | −1093.6 | 57.3 | −0.28 | 0.628 | 4.92 | 0.488 | −1147.9, −1041.7 | 30 | `▁▁▁▁▁▁▁▂▄▇▇▁▁▁▂▂▄██▂▁▁▁▁` |
| 0.12 | −1129.0 | 55.7 | 0.42 | 0.470 | 4.59 | 0.239 | −1155.9, −1043.3 | 53 | `▁▁▁▁▁▁▁▂▃▄▆█▁▁▁▁▁▂▂▃▂▁▁▁` |
| 0.11 | −1156.4 | 50.9 | 0.31 | 0.263 | 0.29 | 0.409 | −1162.0, −1148.4 | 113 | `▁··▁▁▁▁▁▁▂▃▂▃▅██▁▁▁▁▂▂▁▁` |
| 0.10 | −1181.6 | 51.5 | −0.52 | 0.258 | 1.18 | 0.304 | −1195.9, −1148.9 | 389 | `▁▁▁▁▁▁▁▁▁▁▂▄▄▃▄▆█▃▁▁▁▁▁▁` |

*Per chain,* both BC > 5/9 and D > 2 hold at:
- chain 0: b = 0.24–0.27, 0.21 and 0.13–0.14;
- chain 1: b = 0.25–0.28 and 0.13–0.15;
- chain 2: b = 0.23, 0.24, 0.26 and 0.13–0.14;
- chain 3: b = 0.27 and 0.13, and D > 2 with BC 0.52–0.56 at 0.26 and 0.12.

*Each chain's ESS* for ℓ₁'s mean is 2–19 at b = 0.22–0.27 on every chain, and 81–1,075 at 0.16–0.18. Above 0.30 it depends on the chain: 7–21 on chain 1, 14–45 on chain 2, 32–1,671 on chains 0 and 3.

*M at the window's rungs, pooled* (counts of 8,000):

| b | M's counts | mean | BC of M |
| --- | --- | --- | --- |
| 0.30 | 6:198 7:1259 8:3760 9:2084 10:643 11:54 12:2 | 8.24 | 0.322 |
| 0.27 | 6:292 7:1258 8:3918 9:1976 10:508 11:45 12:3 | 8.16 | 0.295 |
| 0.25 | 6:689 7:1773 8:3777 9:1426 10:294 11:39 12:2 | 7.87 | 0.309 |
| 0.24 | 6:1388 7:2486 8:2951 9:972 10:178 11:25 | 7.52 | 0.381 |
| 0.23 | 6:2253 7:2899 8:2191 9:540 10:107 11:10 | 7.17 | 0.444 |
| 0.22 | 5:2 6:3175 7:2890 8:1581 9:306 10:44 11:2 | 6.89 | 0.513 |
| 0.20 | 5:1 6:4971 7:2145 8:743 9:134 10:6 | 6.51 | 0.669 |
| 0.17 | 4:56 5:55 6:6190 7:1406 8:253 9:36 10:4 | 6.23 | 0.365 |
| 0.14 | 3:34 4:1192 5:472 6:5179 7:990 8:115 9:18 | 5.79 | 0.392 |
| 0.13 | 3:141 4:2490 5:896 6:3662 7:704 8:92 9:15 | 5.33 | 0.465 |
| 0.12 | 3:461 4:3954 5:1232 6:1883 7:385 8:76 9:9 | 4.76 | 0.533 |

*ℓ₁ given M, pooled* (histograms over the rung's pooled range: [−1053, −792] at b = 0.26, [−1054, −794] at 0.25, [−1079, −794] at 0.24):

| b | M | share | ℓ₁ mean ± sd | histogram |
| --- | --- | --- | --- | --- |
| 0.26 | 6 | 4.7% | −895.9 ± 79.2 | `▁▁▄▁··········▁▁▁▄▃█▁···` |
| | 7 | 17.8% | −898.6 ± 68.7 | `▁▁▁▁▁▂▁▃▇█▁···▁·▁▁▅▃▅█▂▁` |
| 0.25 | 6 | 8.6% | −957.5 ± 87.2 | `▁▃█▄···········▁▁▄▂▅▁···` |
| | 7 | 22.2% | −932.0 ± 61.8 | `▁▁▁▁▂▂▁▃██▂···▁▁▁▁▃▁▂▃▁▁` |
| | 8 | 47.2% | −889.1 ± 48.2 | `▁▁·▁▁▁▁▂▃▃▃█▆▃▂▃▃▇▆▆▅▃▂▁` |
| | 9 | 17.8% | −871.6 ± 41.4 | `·····▁▁▁▁▂▂▅▃▄▂▃▄▆▆█▇▄▂▁` |
| 0.24 | 6 | 17.3% | −999.8 ± 66.0 | `▁▁▁▄█▂·········▁▁▁▁▁▁···` |
| | 7 | 31.1% | −958.5 ± 45.6 | `·▁▁▁▂▁▂▂▂▇█▃···▁▁▁▁▁▁▂▁·` |
| | 8 | 36.9% | −903.1 ± 49.9 | `··▁▁▁▁▁▂▂▅▅▄█▇▄▃▃▇█▅▅▂▂▁` |

From b = 0.19 to 0.16, M = 6's ℓ₁ is one lump at −1033 to −1037 (sd 10–14). At b = 0.12–0.15 the lumps are M classes instead: M = 4 at −1148 to −1156 (sd 12–20), with M = 5 beside it, against M = 6 at −1040 to −1057.

*The tests' sensitivity*: iid mixtures of the pooled b = 0.32 phase, N(−822.9, 16.0²), and the b = 0.16 phase, N(−1039.1, 25.0²), at each minority weight; 100 samples per cell; the share flagged.

| n | minority phase | its weight | BC > 5/9 | ΔBIC > 10 and D > 2 |
| --- | --- | --- | --- | --- |
| 50 | b = 0.16 | 0.01, 0.02, 0.05, 0.10 | 38%, 67%, 88%, 98% | 35%, 66%, 86%, 98% |
| 50 | b = 0.32 | 0.01, 0.02, 0.05, 0.10 | 26%, 52%, 76%, 98% | 41%, 62%, 84%, 98% |
| 2,000 | b = 0.16 | 0.01, 0.02 | 87%, 100% | 100%, 100% |
| 2,000 | b = 0.32 | 0.01, 0.02 | 0%, 72% | 100%, 100% |
| 8,000 | b = 0.16 | 0.01, 0.02 | 98%, 100% | 100%, 100% |
| 8,000 | b = 0.32 | 0.01, 0.02 | 0%, 88% | 100%, 100% |

At weights of 0.05 and above every test flags every sample at n = 2,000 and 8,000; at 0.20 and above, at n = 50 too.

ΔBIC alone is not specific here. It exceeds 10 at every rung of the window, b = 0.40 included (ΔBIC 1,318 with D 1.51 and BC 0.47), since two Gaussians also fit one skewed lump better than one.

*What R1 shows.*
- *ℓ₁ is bimodal at fixed b inside the band.* At b = 0.23–0.28 the pooled ℓ₁ has BC 0.52–0.65 and D 2.2–3.3, and its histogram holds separated lumps (three at b = 0.24, near −1026, −954 and −865). Each chain shows the same at some rung of b = 0.23–0.28, so the phases coexist within a chain; this is not chains sitting in different phases.
- *The two phases are not M classes.* At b = 0.24–0.26, ℓ₁ given M = 6 or M = 7 has two lumps 120–190 nats apart, with empty bins between them. M itself is unimodal over the integers at every rung: its weight slides from 8 to 6 through 7, and M = 7 is never depleted.
- *The lower phase holds its place while its weight grows.* Its lump sits at ℓ₁ ≈ −1029 to −1033 (sd 7–9) from b = 0.23 to 0.18, while its weight in the two-Gaussian fit rises from 25% to 78%.
- *Outside the band ℓ₁ is one lump:* skewed, with a lower tail, at b ≥ 0.29 (skew −0.8 to −1.6); sharp at 0.16–0.18 (BC 0.09–0.41, D under 0.7).
- *A second split at b ≈ 0.12–0.15* divides M = 6 (ℓ₁ ≈ −1040) from M = 4 and 5 (ℓ₁ ≈ −1150), with D 2.4–5.9 per chain. There the phases are M classes.
- *Of the brief's two readings,* these fit the first: two phases coexisting at fixed b, not a unimodal broad crossover.
- *Strength.*
  - The lumps and their positions repeat on all four chains.
  - Their weights at a rung rest on 2–19 effective samples per chain, so the weights are not resolved.
  - At that ESS a 5% minority phase is flagged 76–88% of the time, and a 1% one 26–41%. At the band's upper edge (b = 0.29–0.30: pooled minority weight 3.5–6.6%, D 1.9) the other phase may be present at weights these tests miss.

#### R2. SMC against the particle count

Fitting SMC at its defaults (conditional ESS 0.99, 5 mutation sweeps, resampling below 0.5, the integrated form, a budget of 1,000 stages), four replicates from S10's master seed 53. The 1,000-particle run reproduces C4's exactly (ln Ẑ −875.76 ± 4.21; replicate means of the log posterior −904.0, −914.3, −917.5, −914.1).

| particles | s | stages | ln Ẑ per replicate (genealogy error 0.60–0.83) | sd (range) | pooled ± error (within, between) | log post mean per replicate (M4's) | σ's 95% interval |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1,000 | 29–31 | 179, 182, 197, 198 | −885.57, −865.22, −874.54, −877.70 | 8.42 (20.35) | −875.76 ± 4.21 (0.42, 4.21) | −904.0, −914.3, −917.5, −914.1 | [1.0033, 1.1506]: missed |
| 4,000 | 107–127 | 181, 199, 186, 198 | −868.57, −875.51, −869.47, −883.93 | 7.08 (15.36) | −874.37 ± 3.54 (0.39, 3.54) | −910.4, −915.0, −908.8, −923.7 | [0.9971, 1.1425]: covered |
| 16,000 | 561–661 | 182, 189, 189, 169 | −855.22, −868.86, −865.66, −857.50 | 6.49 (13.64) | −861.81 ± 3.25 (0.41, 3.25) | −895.9, −905.3, −905.3, −907.2 | [0.9948, 1.1431]: covered |

*Families* (the initial particles with a descendant alive) at the first stage reaching each b, then the b where they first fall to two:

| particles | replicate | 0.10 | 0.16 | 0.24 | 0.32 | 1 | ≤ 2 at b |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1,000 | 0 | 113 | 45 | 15 | 4 | 1 | 0.466 |
| | 1 | 104 | 23 | 1 | 1 | 1 | 0.215 |
| | 2 | 121 | 34 | 2 | 2 | 1 | 0.238 |
| | 3 | 71 | 32 | 3 | 3 | 2 | 0.545 |
| 4,000 | 0 | 566 | 124 | 6 | 1 | 1 | 0.298 |
| | 1 | 484 | 139 | 22 | 5 | 1 | 0.664 |
| | 2 | 528 | 80 | 13 | 3 | 2 | 0.398 |
| | 3 | 550 | 126 | 24 | 7 | 5 | — |
| 16,000 | 0 | 2,181 | 475 | 9 | 5 | 1 | 0.334 |
| | 1 | 2,175 | 304 | 54 | 10 | 4 | — |
| | 2 | 1,149 | 245 | 32 | 7 | 3 | — |
| | 3 | 1,279 | 280 | 1 | 1 | 1 | 0.209 |

*The kernel*, at every particle count, accepts 8.3–9.0% of births, 10.8–12.4% of deaths and 23.6–25.9% of redraws. The support rule refuses 16.8–17.3% of births and 4.9–5.2% of redraws.

*Each replicate's two largest basins at b = 1* (weighted fraction, then the log mass ln Ẑ + ln f), and its share by the count of x1·x2 terms:

| particles | replicate | ln Ẑ | largest | second | x1·x2 terms: 2, 3, 4 |
| --- | --- | --- | --- | --- | --- |
| 1,000 | 0 | −885.57 | A 58.6% (−886.1) | B 6.4% (−888.3) | 85%, 15%, 0% |
| | 1 | −865.22 | B 49.5% (−865.9) | a4 17.0% (−867.0) | 0%, 99%, 1% |
| | 2 | −874.54 | a4 35.5% (−875.6) | B 13.5% (−876.5) | 0%, 98%, 2% |
| | 3 | −877.70 | B 54.4% (−878.3) | T 8.9% (−880.1) | 0%, 92%, 8% |
| 4,000 | 0 | −868.57 | B 47.5% (−869.3) | T 21.2% (−870.1) | 4%, 91%, 5% |
| | 1 | −875.51 | B 30.6% (−876.7) | a8 26.1% (−876.9) | 1%, 53%, 40% (5: 6%) |
| | 2 | −869.47 | a5 37.7% (−870.5) | a6 21.9% (−871.0) | 93%, 6%, 0% |
| | 3 | −883.93 | a7 32.1% (−885.1) | B 19.7% (−885.6) | 8%, 87%, 5% |
| 16,000 | 0 | −855.22 | a1 68.8% (−855.6) | a1 + x6·x8 6.7% (−857.9) | 0%, 1%, 99% |
| | 1 | −868.86 | B 64.1% (−869.3) | a4 5.2% (−871.8) | 15%, 85%, 0% |
| | 2 | −865.66 | a2 49.5% (−866.4) | a3 23.7% (−867.1) | 0%, 34%, 66% |
| | 3 | −857.50 | A 38.3% (−858.5) | B 20.3% (−859.1) | 53%, 33%, 14% |

*What R2 shows.*
- *The replicates' evidence does not agree at any N.* Across replicates, ln Ẑ's sd is 8.4, 7.1 and 6.5 nats at 1,000, 4,000 and 16,000 particles (ranges 20.4, 15.4, 13.6). The between-run error is 4.2, 3.5 and 3.3, against a within-run (genealogy) error of 0.4. The review notes put agreement at about a nat.
- *ln Ẑ rises with N.* The pooled estimate is −875.8, −874.4, −861.8, the largest replicate −865.2, −868.6, −855.2. It has not settled by 16,000 particles.
- *The lineages collapse inside the band at every N.* The families alive at b = 0.16 grow with N (23–45, 80–139, 245–475). By b = 0.32 only 1–4, 1–7 and 1–10 remain, and 1–5 at b = 1. One 16,000-particle replicate is down to a single family at b = 0.21.
- *The replicates hold different basins.* At 16,000, replicate 0, with the largest ln Ẑ (−855.2), holds a1 at 68.8%. a1 has four x1·x2 terms and no x1 or x2 main effect, and no MCMC chain or other replicate holds it at 0.1% or more. Replicate 3 (−857.5) holds A and B, replicate 1 B, replicate 2 a2 and a3.
- *A basin's mass does not agree across replicates.* Its ln Ẑ + ln f moves with its replicate's ln Ẑ. A's reads −858.5, −872.8 and −875.2 in the three 16,000-particle replicates that hold it, B's −859.1, −869.3 and −870.0. Within a replicate the fractions are the surviving lineages' local exploration, so they inherit ln Ẑ's spread.
- *The log posterior (M4's),* against M4's best chain's median of −892.4: replicate means from −895.9 to −923.7. This compares densities, not masses.
- *σ:* the pooled 95% interval covers 1 at 4,000 and 16,000 (lower ends 0.997 and 0.995) and misses at 1,000 (1.003), as C4's did.
- *Strength.* Four replicates per N, so the spread is itself rough. a1 rests on one replicate.

#### R3. Mutation against particles

At 1,000 particles, 20 mutation sweeps per stage against R2's 5, all else as R2:

| sweeps | s | stages | ln Ẑ per replicate (genealogy error) | sd (range) | pooled ± error (within, between) | families at b = 0.16, 0.24, 0.32, 1 | log post mean | σ's 95% interval |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5 | 29–31 | 179, 182, 197, 198 | −885.57, −865.22, −874.54, −877.70 (0.83 each) | 8.42 (20.35) | −875.76 ± 4.21 (0.42, 4.21) | 23–45, 1–15, 1–4, 1–2 | −904.0 to −917.5 | [1.0033, 1.1506] |
| 20 | 64 | 185, 181, 184, 186 | −872.86 ± 0.74, −873.38 ± 0.49, −858.96 ± 0.81, −862.81 ± 0.82 | 7.24 (14.42) | −867.00 ± 3.62 (0.36, 3.62) | 46–103, 3–36, 3–10, 2–5 | −900.4 to −909.0 | [0.9951, 1.1442] |

At 20 sweeps, B is the largest basin in every replicate: 48.4% (−873.6), 27.2% (−874.7), 60.9% (−859.5) and 41.7% (−863.7). The second largest are A 16.1%, a5 17.4%, AC 7.0% and a4 25.9%.

*The kernel's acceptance by kind* is unchanged: births 8.5–8.7%, deaths 11.3–11.6%, redraws 23.9–24.5%, against 8.4–9.0%, 10.8–12.1% and 23.8–25.5% at 5 sweeps.

*What R3 shows.*
- At the same per-sweep acceptance, the extra sweeps add moves, not better moves.
- ln Ẑ's mean rises from −875.8 to −867.0, and its sd falls from 8.4 to 7.2 (range 20.4 to 14.4).
- The families at b = 0.32 rise from 1–4 to 3–10, about the 16,000-particle run's 1–10, in 64 s against 561–661 s. At b = 1 they rise from 1–2 to 2–5.
- The replicates still disagree by 14 nats.
- *Strength:* one seed, four replicates.

#### R4. The support's share of the barrier

*The refusals* on M4's four untempered chains over their sampling phase. S10's support minimum is 20 nonzero points and the angle floor 10⁻¹⁰. Each cell is the share of that kind's proposals; deaths are never refused.

| chain | births | refused: support | redraws | refused: support | deaths |
| --- | --- | --- | --- | --- | --- |
| 0 | 687 | 16.2% | 658 | 5.9% | 655 |
| 1 | 640 | 15.2% | 650 | 4.3% | 710 |
| 2 | 678 | 14.9% | 660 | 5.6% | 662 |
| 3 | 642 | 16.2% | 651 | 7.1% | 707 |
| pooled | 2,647 | 15.6% | 2,619 | 5.7% | 2,734 |

- *Every refusal is the support rule's.* The angle floor, repeated terms, failed factors, the cap, variance inflation and exact fits refuse none.
- *The refused terms' supports:* of the 563, 56 have no nonzero point and 270 have 1–9. The other 237 (42%) have 10–19, which half the minimum would admit.
- *Along R1's dense ladders* (four chains, every rung): births 15.9% support and 0.03% repeated terms; redraws 4.7% support and 0.03% repeats; the angle floor none. The shares are flat in b: births 13.6–17.1% and redraws 3.8–5.4% at each window rung, 18.7% and 5.7% below 0.10.

*The rule's outcome.*
- The angle floor refuses 0% of every kind's proposals, so the reruns at 10⁻¹⁴ and 10⁻¹⁸ were not taken. Counting the repeated terms and failed factors with it changes nothing on M4's chains.
- The support rule refuses 15.6% of births, so M4's chains were rerun with its minimum at 10, the floor unchanged.
- *This rerun changes F:* terms with 10–19 nonzero points join the support, so its chains sample a different posterior. It measures the rule's share of the stall, not the same posterior.

| minimum | chain | median log post (5–95%) | below the best | M mean (range) | basin most held | x1·x2 terms | R̂ of M | R̂ of the log post |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20 (S10's) | 0 | −933.9 (−959.3, −930.8) | 41.5 | 10.28 (10–12) | G, 75.0% | 5: 100% | 2.35 | 2.05 |
| | 1 | −892.4 (−912.7, −890.5) | 0 | 8.21 (8–10) | A, 79.6% | 2: 87%, 3: 13% | | |
| | 2 | −926.2 (−942.2, −913.8) | 33.8 | 10.77 (10–12) | S, 61.3% | 3: 95%, 4: 6% | | |
| | 3 | −916.0 (−935.4, −913.1) | 23.6 | 9.17 (9–11) | B, 82.7% | 3: 100% | | |
| 10 | 0 | −892.9 (−918.5, −890.1) | 0 | 9.22 (9–11) | B, 79.2% | 3: 100% | 2.13 | 1.99 |
| | 1 | −905.4 (−930.5, −902.2) | 12.5 | 9.38 (9–12) | B, 67.5% | 3: 99% | | |
| | 2 | −981.6 (−993.4, −978.8) | 88.6 | 13.06 (12–15) | a11, 87.0% | 3: 96%, 2: 4% | | |
| | 3 | −903.3 (−927.4, −897.6) | 10.4 | 9.65 (9–11) | a10, 51.0%; a5, 37.2% | 2: 100% | | |

The R̂s are the run's (`Bass.Fit` on the same config), whose log posterior includes τ's prior density; the medians are M4's. Under the minimum of 10, the support rule still refuses 8.4% of births and 2.5% of redraws.

*What R4 shows.*
- *The angle floor has no share in the stall on these chains.* It refuses nothing, at the cold rung or along the ladder. A replacement's intermediate state (the old term with its substitute) is never refused by the floor here.
- *The support rule refuses one birth in six, flat in b,* and 42% of the terms it refuses have 10–19 nonzero points.
- *Halving the minimum does not unstick the chains.* R̂ of M is 2.13 against 2.35, and the medians spread 88.6 nats against 41.5. Chains 0 and 1 now share B; chain 2 sits 88.6 nats below in a 13-term basin with two x1·x2·x5 terms; chain 3 holds two x1·x2 terms with doubled x1 and x2 main effects.
- *Strength:* the refusal shares rest on 2,600–2,700 proposals per kind and agree across chains within 3 points. The rerun is one run of four chains.

#### MCMC's basins for comparison

*M4's untempered chains* (S10's, 2,000 sweeps each):
- chain 0: G 75.0%, then eleven basins at 1.0–5.3%, each G with one or two more terms;
- chain 1: A 79.6%, B 11.7%, A + x1·x6 4.0%;
- chain 2: S 61.3%, T 30.7%, S with a fourth x1·x2 term 5.5%;
- chain 3: B 82.7%, then seven basins at 1.3–3.4%, each B with one more term.

Only B, above 1%, is held by two chains (1 and 3).

*S10's tempered chains* (R1's Part 1, the cold rung, 2,000 scans each), each with its share by the count of x1·x2 terms:
- chain 0: B 50.2%, T 10.4%, 2 × x1 + 4 × x1·x2 + x2 + 2 × x3 + x4 + x5 8.9%, AC 8.1%, a12 6.7%; x1·x2 terms 2: 1%, 3: 89%, 4: 10%;
- chain 1: B 55.6%, a8 11.5%, a8 + x6·x8 7.7%; 2: 1%, 3: 70%, 4: 29%;
- chain 2: a9 63.4%, T 9.8%; 2: 3%, 3: 98%;
- chain 3: T 33.2%, a10 27.5%, A 12.9%, B 8.9%; 2: 50%, 3: 50%.

*What these show.*
- Under the structural key the untempered chains hold four disjoint sets of basins.
- The tempered chains each mix several, the most held at 33–63%, and B is the most held on two of them.
- The basins split on the x1·x2 count: two in A, a5 and a10; three in B, S and T; four in a1, a2 and a8; five in G.

#### Run times

Release; SMC's replicates in parallel within the machine's bound, the chains in parallel elsewhere. The second run's 16,000-particle fit overlapped a concurrent session's work in this checkout, so the times are rough.
- *R1:* 1 min 23 s – 1 min 48 s: Part 1 6.1–8.7 s plus S10's tempered fit; Part 2 43.6–57.0 s, with chain 0's unrecorded twin beside it; the readings' fits the rest.
- *R2:* 29–31 s, 1 min 47 s – 2 min 7 s and 9 min 21 s – 11 min 1 s at 1,000, 4,000 and 16,000 particles. *R3:* 64 s. The ranges are two runs each.
- *R4:* 1–6 s, with its two `Bass.Fit` checks.

#### The review passes and their fixes

1. *Correctness, as the brief lists it.*
   - ℓ₁ is each replica's untempered log marginal at every rung. The harness asserts it equals `LogLikelihood(1)` bit for bit at each rung.
   - The structural key ignores knots and signs. The harness asserts that a retained structure keeps its basin with every hinge's sign flipped and knot moved, and leaves it when one factor moves to another variable (added in this pass).
   - The basins' log masses add up to ln Ẑ within 10⁻⁹ relative, asserted per replicate. Each replicate's weights sum to 1/R within 10⁻¹², and its fractions to one within 10⁻⁹.
   - Phase 0's readings are unchanged. `BassMixing.cs` only gained `Fixture.Basin`. `BassMixingModeTests` and `BassMixingBandTests` reran with M4's and M5's readings unchanged (M4's Friedman table, M5's Σr and Σh on every chain). R4's chains reproduce M4's medians and R̂, R1's Part 2 M5's rows, and R2's 1,000-particle run C4's.
2. *The harness's own checks.*
   - Its chains are the runs': M4's chains hold the untempered fit's M draw for draw (also under the halved minimum), and Part 1's cold rungs hold the tempered fit's.
   - The refusal recorder stands in for each chain's target and leaves the path unchanged: chain 0's dense ladder without it holds the same M and ℓ₁ at every rung and scan.
   - Each sweep's refusal is read once, with its reason matching the candidate's, and each numerical refusal is split by a whole load.
3. *Fixes after the first read.*
   - ΔBIC flagged every rung, b = 0.40 included. It is reported beside BC and D, and the readings rest on those and the histograms.
   - M was unimodal where ℓ₁ was not. ℓ₁ given M was added to reconcile the two, and shows the split within M.
   - The sensitivity check was taken at the chains' length only. It now includes n = 50, near the band's ESS.
   - The key split the samples into dozens of basins, so the x1·x2 count was added beside it as a coarser view.
4. *Coverage.*
   - *One-replicate or one-chain readings:* a1 (16,000 particles, replicate 0); R4's rerun is one run of four chains; R3 is one seed.
   - *The bimodality test at the band's edges:* above, in R1's strength.
   - *The stage budget:* no run reached it (at most 199 of 1,000 stages), and every replicate reached b = 1.
   - *Repeatability:* a second run of the final harness gave the first run's readings to the printed digit, apart from the sensitivity table, whose streams changed when n = 50 was added.

#### For the main session

- *argon's mechanism is not in these chains.* The angle floor refuses nothing on M4's chains or along the dense ladders, so a replacement's intermediate state is never refused by the floor. The support rule refuses a sixth of births, but halving its minimum leaves the chains apart (R̂ of M 2.13).
- *The Opus chat's premise, rejection near 0.95 at Δb = 0.01,* is M5's already-corrected reading: the dense pairs reject at most 0.28–0.46. Its coexistence reading holds in ℓ₁ nonetheless (R1).
- *Λ̂ for OQ2's record:* Σr is 2.0–2.4 over the band and 8.6–9.0 over the self-sized ladder; 2.7–3.0 and 9.2–10.2 on the dense one.
- *SMC's evidence has not converged in N by 16,000 particles.* One replicate finds a basin, a1, with the largest ln Ẑ, which no chain visits. This bears on C4's and D10's readings of the modes' masses.
- *A concurrent session* began editing BASS production and test files in this checkout at 21:06, after this chip's last build (20:52). Every run used that build, so no reading depends on those edits, and none of them is in this chip's commits.
- No route is recommended.
- `artifacts\chip-friedman\` is deleted.
