# Friedman: a path from a pilot reference (exploration chip brief)

Date: 2026-10-29. The user asked for this test after the Friedman measurement ([its report](friedman-band-brief.md#reports)). That measurement found the tempering band at b ≈ 0.23–0.28 to be a first-order transition, its two phases coexisting at fixed b. Fitting SMC from the prior collapses its lineages there at every particle count read, and its replicates' ln Ẑ disagrees by 6.5–8.4 nats up to 16,000 particles.

The hypothesis is that a path starting inside the fitted phase never crosses the transition. That is Surjanovic et al. 2022's variational reference for tempering, here in its simplest discrete form. This chip tests it, and changes nothing in production.
- *The question:* does SMC along a path from a pilot-built reference reach replicate agreement, ln Ẑ within about a nat, where SMC from the prior does not?
- *Not in scope:* building the reference into BASS, any production change, any default, and any recommendation of a route. The readings go to the user.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the Friedman measurement's brief and report;
- C4's report ([stage8-fitting-smc-brief.md](stage8-fitting-smc-brief.md#reports));
- D05 in [port-fitting-brief.md](port-fitting-brief.md) (the likelihood path, the normalizer as evidence).

In the code:
- *the SMC core,* `src/maths/samplers/smc/`;
- *C4's adapter,* `BassSmc.cs` (`BassParticle`, `BassParticles`, `BassParticleTarget`, `BassSmcKernel`) and `BassSmcRun.cs`;
- *the Friedman measurement's harness,* `tests/regression/BassFriedmanBandTests.cs` (its basins, `Fixture.Basin`, and its SMC runs);
- *Phase 0's harness,* `tests/regression/BassMixing.cs` (`Score`, `Loaded`, `Proposed`, `LogPosterior`), which scores BASS's moves outside the chain;
- `BassPriorSampler` (its exact draws on F and its attempts) and `BassStructurePrior` (the unconstrained term measure's log density).

## The design

**The reference** q is a distribution over (structure, τ) that can be drawn exactly and whose density is exact.
- q(S, τ) = [ε·p₀(S) + (1 − ε)·e(S)]·p(τ).
  - p₀ is D01's unconstrained structure measure, normalized on its own space, which contains F.
  - e is the empirical distribution of the pilot runs' structures (below), each distinct structure weighted by its pooled frequency.
  - p(τ) is τ's prior.
- A draw from p₀ that falls outside F gets zero weight at the first step with b > 0. That is exact, since the posterior is zero there.
- ε is 0.05 and 0.5, each a separate arm.

**The path** is geometric from q to the posterior π(S, τ) ∝ p_F(S)·p(τ)·L(S, τ), where p_F is the prior on F:
- log γ_b = (1 − b)·log q + b·log π̃, as a target on the likelihood path with base log q and "likelihood" Δ = log π̃ − log q;
- stages by conditional ESS at 0.99, multinomial resampling below half the population, C4's defaults otherwise;
- the mutation: BASS's structure moves and τ's slice step, with their acceptance evaluated against γ_b. The proposals are D01's, unchanged; only the target differs.

**The evidence.** The run's normalizer estimates ∫ π̃ / ∫ q, and q is normalized. First read how `BassTarget` states the structure prior.
- If it is the unconstrained measure restricted to F, ln p(y) under the prior on F is the normalizer less ln Z_F, where Z_F is F's share of the unconstrained measure. Estimate Z_F from `BassPriorSampler`'s acceptance over many attempts, and carry its error into the reported one.
- If it is already normalized on F, nothing is subtracted.
- Say which it is.

**The pilots,** run with seeds disjoint from the test's replicates:
- M4's four untempered chains;
- S10's tempered Friedman fit (its cold rungs);
- fitting SMC from the prior at 4,000 particles, four replicates (the measurement's R2 setting).

Their pooled draws make e. A second reference leaves out one source, the SMC pilot, to read how much the result depends on what the pilots saw.

**Where it runs.** In tests only, with no production change. The mutation kernel for γ_b is the harness's own: an `ISmcKernel` over `BassParticle` that proposes with `BassTermProposal` and the particle's design, as `BassMixing`'s scorers do, and accepts against γ_b. If the kernel cannot be written without a production seam, stop and report the seam needed, without adding it: a ridge chip is changing BASS's production files in this checkout.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **A ridge chip works beside you** on BASS's production files (`BassConfig.cs`, `BassTarget.cs`, `BassSweep.cs`, `BassDesign.cs` and others), replacing `SampledG` and `FixedG` with a two-axis coefficient prior. Build against the checkout as it is. If its commits change an API your harness uses, adapt your harness; never edit its files. If a build fails mid-change in a file you do not own, wait and retry.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-reference\` to every `dotnet build` and `dotnet test`. Build in Release for the runs. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-reference\` when you finish.
- **Your files:**
  - a new harness, `tests/regression/BassFriedmanReferenceTests.cs`, behind `BASS_QUALIFICATION` = 1;
  - additions to `BassMixing.cs` that leave Phase 0's readings unchanged;
  - this brief's Reports section.
- **Commit your own files only, by explicit path,** in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked notes and the ridge chip's uncommitted edits stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.

## Checks before the runs

- *The kernel is exact for γ_b.* On fixture A's enumerable F (`BassSweepTests.EnumerateF`), run it as a chain at b = 0.3 and 0.7 from a small reference, and compare the visit frequencies with γ_b computed by enumeration.
- *The reference is exact:* q's draws against its density, within 5 SE, over p₀'s and e's components separately and mixed.
- *The evidence:* on fixtures B and C, SMC from a reference reproduces the oracle's exact ln Z (−28.589 and −30.903) within its error, Z_F's included. The reference here is built from a short pilot on the same fixture.
- *At ε = 1* the reference is p₀, and the run is SMC from the unconstrained prior: check it against C4's from-prior SMC on fixture B.

## The runs

On S10's Friedman fixture (n = 500, seed 51; master seed 53 for the replicates), for each reference (all pilots, and all but the SMC pilot) and each ε, at 1,000 and 4,000 particles, with four replicates and C4's other defaults, read:
- each replicate's ln Ẑ (Z_F's correction included), with its genealogy error, and the replicates' sd, range and pooled between-run error;
- each replicate's basin fractions under the measurement's structural key, its two largest basins with their log masses, and its share by the count of x1·x2 terms;
- the families at b = 0.25, 0.5, 0.75 and 1, and the stage count;
- σ's 95% interval;
- the time, against the measurement's from-prior runs at the same particle count (R2).

Name whether any replicate reaches the four-x1·x2 basin (a1, a2) the measurement's 16,000-particle replicate found, and whether the references' pilots held it.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - the kernel's acceptance against γ_b, with D01's proposal ratios unchanged;
  - q's density matching its sampler;
  - Z_F's correction and its error;
  - the pilots' seeds disjoint from the replicates';
  - Phase 0's readings unchanged.
- *Coverage:*
  - a reading resting on one replicate;
  - a run whose families collapse anyway;
  - draws from p₀ outside F.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the checks with their tolerances;
- the runs' tables, then what they show, kept apart;
- the times;
- the review passes;
- anything the main session must know.

Recommend no route.

## Reports

Newest first.

### The reference path measured — 2026-10-29

*Landed* (tests only; no production change, and no seam was needed):
- `abdcd4d`: the harness `tests/regression/BassFriedmanReferenceTests.cs`, behind `BASS_QUALIFICATION` = 1.
  - Its checks: `Kernel_SamplesTheBridgeOnFixtureA`, `Reference_DrawsMatchItsDensityOnFixtureA`, `Evidence_MatchesTheExactLogZ` and `PriorReference_AgreesWithC4sSmcOnFixtureB`.
  - Its runs: `Friedman_PilotsAndTheirReference` and `Friedman_PathFromPilotReference`.
  - `BassMixing.cs` is untouched.
- *No seam.* `ReversibleJumpChain` and `BassTauBlock` take any `ISamplingTarget<BassState>`, so the kernel is `BassSmcKernel`'s sweep with γ_b in the model target's place. The harness uses existing internals only: `BassParticle`'s constructor, `BassState.TermArray` and `BassSmcAccount.Add`.

*The construction.*
- *γ_b* runs on the likelihood path as the SMC core takes it. ln q is the untempered prior and Δ = ln π̃ − ln q the likelihood, so the stage density is ln q + b·Δ.
- *How `BassTarget` states the structure prior.* It is D02's count weight on F, w = Γ(h₁ + M)/(Γ(h₁)(h₂ + 1)^M)·Π π(term). That fits neither of the brief's two cases exactly.
  - It is the unconstrained measure restricted to F up to a constant, C = Σ_{M ≤ M_max} Γ(h₁ + M)/(Γ(h₁) M! (h₂ + 1)^M), the normalizer of `BassPriorSampler`'s count law. On Friedman, ln C = 0.9531, which is h₁ ln(1 + 1/h₂) since the cap lies far out.
  - The harness's π̃ takes p₀ = w/C on F, the law of `BassPriorSampler`'s proposal over sets of distinct terms. So the run's normalizer estimates Z_F·p(y), and ln p(y) = normalizer − ln Z_F.
  - Z_F is `BassPriorSampler`'s acceptance over 40,000 accepted draws. Its error, √((1 − Z_F)/accepted), joins each replicate's error and the pooled one in quadrature.
  - Equivalently, against w itself, ln p(y) = normalizer − ln C − ln Z_F.
- *q* is ε·p₀ + (1 − ε)·e, with u drawn from τ's prior.
  - e pools the pilots' draws: an MCMC draw has mass 1, an SMC particle N·R·w, so each replicate's sum is N.
  - A p₀ draw outside F (an intercept term, a repeat, a support or angle refusal) is a particle with no design, whose γ_b is 0 at every b > 0. The kernel asserts at every stage that none holds weight, and none is ever moved.
- *The kernel* is C4's sweep. The structure move runs on D01's schedule through a reversible-jump chain over the particle's own `BassStructureMoves`, then `BassTauBlock` steps on the bridge. q's τ factor is p(u), so γ_b's τ conditional is p(u)·e^{b·ℓ₁}.
- *Seeds.* The replicates' four seeds sit under `SeedTree.Derive(53, 9)[8]`, a child of master seed 53 that no pilot derives (a fit takes the first 2·chains children, fitting SMC the first R); the harness asserts it. Z_F and q's Friedman checks take later children. The arms share the four seeds (common random numbers).

*The checks* (all pass; C4's tolerances where C4 had one):

| check | tolerance | reading |
| --- | --- | --- |
| the kernel as a chain on fixture A (2,263 structures), reference ε = 0.5 over a 12-structure pilot, 4 × 30,000 sweeps, at b = 0.3 and 0.7 | M's law 5 SE + 0.004; each term's inclusion (at least 0.5%) and the ten heaviest structures 5 SE + 0.005; E[u] 5 SE + 0.001 | largest \|z\| 2.21 (b = 0.3) and 2.46 (b = 0.7). γ_b sits at total variation 0.251 and 0.080 from the linear form's path from p₀, so the check sees q |
| q's draws on fixture A: p₀, e, and mixed at ε = 0.5, 400,000 each | binomial 5 SE | M's law in F, the share outside F (0.34575 ± 0.00075 against 1 − Z_F = 0.34505), the 15 heaviest structures, u's quartiles and mean: largest \|z\| 2.86. q sums to one over F and its complement to 10⁻⁹ |
| the exact ln Z, by this harness's enumeration with u by quadrature (1,601 nodes) | 0.0015 against the oracle | −28.5886 (B) and −30.9026 (C), against −28.589 and −30.903 |
| Z_F from the sampler's acceptance | 5 SE against enumeration | B −0.30633 ± 0.00257 against −0.30761; C −0.34798 ± 0.00271 against −0.34566 |
| SMC from q on B and C (2,000 particles, 4 replicates; pilot 2 chains × 500 draws) | each replicate and the pool within 4 errors, Z_F's included, + 0.01 | B ε 0.5: −28.5861 ± 0.0318 (z 0.09). B ε 0.05: −28.6518 ± 0.0310 (z −2.03; replicates −2.57, −0.26, −0.23, −2.18). C ε 0.05: −30.9145 ± 0.0119 (z −0.97). C ε 0.5: −30.9203 ± 0.0203 (z −0.85) |
| ε = 1 on B against C4's from-prior SMC (2,000 particles, 4 replicates, seed 2026) | 4 combined errors + 0.01; E[M] 5 SE + 0.01 | −28.5987 ± 0.0295 against −28.5974 ± 0.0266 (integrated) and −28.5926 ± 0.0316 (linear). E[M] 1.9572 against 1.9564 and 1.9572. p₀'s draws in F 0.7379 ± 0.0049 against Z_F 0.7352 |
| q on Friedman | 5 SE | this harness's p₀ in F 0.82070 ± 0.00257 against the sampler's Z_F 0.82364 (z −1.14); e's ten heaviest structures' draw frequencies, \|z\| at most 1.63 |
| the pilots against the measurement | to the printed digit | M4's medians (−892.4, −933.9, −926.2, −916.0) and basin shares. From-prior SMC at 4,000 and 1,000 particles gives R2's ln Ẑ per replicate, so R2 stands on the current code, the ridge chip's commits included |

On fixture B the pilot saw three structures, holding 72.0% of the exact posterior (ln Z + ln P(e) = −28.917). At ε = 0.05 the pool reads 0.063 below the exact ln Z (two errors), and at ε = 0.5 0.003 above. The final populations hold e's structures at 73–82% and 69–73%. On C, e held 86.9% (ln Z + ln P(e) = −31.043), and the pools read 0.012 and 0.018 below.

#### The pilots

S10's Friedman fixture, master seed 53:
- M4's four untempered chains, 2,000 draws each after burn-in;
- S10's tempered fit's cold rungs, 2,000 each;
- fitting SMC from the prior at 4,000 particles (R2's run), weighted.

The labels are the measurement's, from M4's chains. a5 is that report's "x1 + 2 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5".

| | M4 | S10 tempered | SMC 4,000 | e, all pilots | e, no SMC pilot |
| --- | --- | --- | --- | --- | --- |
| mass (draws) | 8,000 | 8,000 | 16,000 | 32,000 | 16,000 |
| distinct structures | 198 | 616 | 11,926 | 12,740 | 814 |
| A | 19.9% | 3.5% | 1.0% | 6.4% | 11.7% |
| B | 23.6% | 28.8% | 25.2% | 25.7% | 26.2% |
| G | 18.8% | 0.0% | 0.0% | 4.7% | 9.4% |
| S | 15.3% | 0.4% | 0.3% | 4.1% | 7.9% |
| T | 7.7% | 13.3% | 9.0% | 9.7% | 10.5% |
| x1·x2 terms 2, 3, 4, 5 | 22, 52, 1, 25% | 14, 77, 10, 0% | 27, 59, 12, 1% | 22, 62, 9, 7% | 18, 64, 6, 13% |
| a1, a2 | 0, 0 | 0, 0 | 0, 0.42% | 0, 0.21% | 0, 0 |

- *Z_F* = 0.82364 (ln −0.19402 ± 0.00210), from 40,000 of 48,565 attempts.
- *e's height over ε·p₀ on its own structures,* h = ln[(1 − ε)·e(S)/(ε·p₀(S))], weighted by e: its 10%, 50% and 90% points are 101.6, 121.2 and 147.8 (all pilots, ε 0.05); 98.7, 118.3 and 144.8 (all pilots, 0.5); 97.5, 123.8 and 174.7 (no SMC pilot, 0.05); and 94.5, 120.9 and 171.8 (no SMC pilot, 0.5). A move off e's structures pays (1 − b)·h at b.
- *The evidence over e's own structures,* ln Σ_{S∈e} p₀(S)·∫ p(u) L(S, u) du − ln Z_F by quadrature in u, with no move:
  - all pilots (12,740 structures): −887.15. By basin: A (359 structures) −887.20, B (3,478) −890.65, a5 (1,664) −891.11.
  - no SMC pilot (814 structures): −887.59. By basin: A (61) −887.60, B (192) −893.14.
  - The best single structure, in A: −889.55.

#### The runs

Four replicates per arm, C4's defaults. ln p(y) includes Z_F's correction, and each replicate's genealogy error is 0.25–0.76. The pooled error is the larger of the within- and between-run errors, with Z_F's (0.002) in quadrature. Families start at N at the first stage, where no arm resamples. Times are two runs each.

| reference | ε | N | stages | ln p(y) per replicate | sd (range) | pooled ± error (within, between) | families at b = 0.25, 0.5, 0.75, 1 | s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all pilots | 0.05 | 1,000 | 49–63 | −887.66, −887.94, −886.64, −886.37 | 0.76 (1.57) | −887.15 ± 0.38 (0.28, 0.38) | 219–479, 122–218, 101–137, 26–48 | 20, 29 |
| all pilots | 0.5 | 1,000 | 60–65 | −887.99, −889.15, −886.75, −885.82 | 1.45 (3.32) | −887.43 ± 0.73 (0.30, 0.73) | 228–258, 109–128, 52–121, 20–63 | 40, 38 |
| no SMC pilot | 0.05 | 1,000 | 49–56 | −888.70, −886.74, −887.84, −886.82 | 0.93 (1.96) | −887.52 ± 0.47 (0.21, 0.47) | 183–209, 125–183, 125–133, 39–98 | 41, 34 |
| no SMC pilot | 0.5 | 1,000 | 61–68 | −888.30, −887.01, −887.75, −886.64 | 0.74 (1.66) | −887.42 ± 0.37 (0.30, 0.37) | 121–257, 74–119, 59–74, 29–71 | 36, 29 |
| all pilots | 0.05 | 4,000 | 59–60 | −887.13, −887.09, −886.50, −886.60 | 0.33 (0.63) | −886.83 ± 0.19 (0.19, 0.16) | 873–920, 493–898, 449–522, 118–174 | 161, 165 |
| all pilots | 0.5 | 4,000 | 71–78 | −886.87, −886.73, −886.32, −885.57 | 0.58 (1.30) | −886.37 ± 0.29 (0.27, 0.29) | 1,011–1,021, 464–510, 257–281, 60–100 | 203, 288 |
| no SMC pilot | 0.05 | 4,000 | 56–64 | −887.27, −886.68, −887.08, −887.18 | 0.26 (0.59) | −887.05 ± 0.16 (0.16, 0.13) | 778–866, 462–525, 462–525, 92–172 | 217, 272 |
| no SMC pilot | 0.5 | 4,000 | 64–68 | −887.73, −887.35, −886.82, −887.13 | 0.38 (0.91) | −887.25 ± 0.19 (0.17, 0.19) | 412–969, 386–418, 240–265, 79–149 | 226, 296 |
| R2, from the prior | — | 1,000 | 179–198 | −885.57, −865.22, −874.54, −877.70 | 8.42 (20.35) | −875.76 ± 4.21 (0.42, 4.21) | at b = 1: 1–2 | 29–31 (R2); 38–56 here |
| R2, from the prior | — | 4,000 | 181–199 | −868.57, −875.51, −869.47, −883.93 | 7.08 (15.36) | −874.37 ± 3.54 (0.39, 3.54) | at b = 1: 1–5 | 107–127 (R2); 176–365 here |
| R2, from the prior | — | 16,000 | 169–189 | −855.22, −868.86, −865.66, −857.50 | 6.49 (13.64) | −861.81 ± 3.25 (0.41, 3.25) | at b = 1: 1–4 | 561–661 (R2) |

*The final populations* (per replicate ranges):

| reference | ε | N | A's share | A's log mass | x1·x2 terms = 2 | on e's own structures | distinct structures | mean log post (M4's) | σ's 95% interval, pooled |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all pilots | 0.05 | 1,000 | 82.5–91.2% | −888.03 to −886.51 | 92.5–99.4% | 48.0–86.8% | 153–326 | −891.7 to −890.1 | [1.0199, 1.1592] |
| all pilots | 0.5 | 1,000 | 65.2–94.0% | −889.57 to −885.88 | 87.3–99.5% | 41.2–78.1% | 219–391 | −895.1 to −889.3 | [1.0169, 1.1592] |
| no SMC pilot | 0.05 | 1,000 | 87.2–98.1% | −888.84 to −886.76 | 97.9–100% | 53.0–93.5% | 84–310 | −892.5 to −888.8 | [1.0223, 1.1595] |
| no SMC pilot | 0.5 | 1,000 | 88.3–93.4% | −888.39 to −886.74 | 98.4–99.7% | 56.1–71.3% | 258–312 | −891.2 to −889.8 | [1.0215, 1.1578] |
| all pilots | 0.05 | 4,000 | 78.1–84.4% | −887.34 to −886.75 | 96.7–98.6% | 39.8–45.8% | 1,230–1,446 | −892.8 to −891.8 | [1.0189, 1.1575] |
| all pilots | 0.5 | 4,000 | 77.0–86.7% | −887.02 to −885.71 | 94.1–99.2% | 19.3–44.7% | 1,149–1,480 | −892.3 to −891.3 | [1.0204, 1.1582] |
| no SMC pilot | 0.05 | 4,000 | 84.1–90.5% | −887.39 to −886.85 | 99.0–99.6% | 30.4–57.0% | 1,154–1,374 | −892.6 to −890.8 | [1.0221, 1.1589] |
| no SMC pilot | 0.5 | 4,000 | 88.5–92.6% | −887.85 to −886.89 | 98.8–99.3% | 49.3–67.6% | 682–1,164 | −891.2 to −890.0 | [1.0219, 1.1581] |

- *The other basins.* Most replicates' second basin is a5 (0.7–13.7%) or B (0.3–9.7%); the rest hold A with one term added or doubled (0.6–4.1%). Four x1·x2 terms hold at most 0.19% of any replicate, and a1 and a2 hold 0.00% in every replicate.
- *Along the path.* The weight on e's own structures is 100.0% at b = 0.25, 0.5 and 0.75 in all 96 such readings (32 replicates, 3 marks). At the first stage it is 96–97% at ε = 0.05 and 55–58% at ε = 0.5; the rest are p₀'s draws in F.
  - A's share starts under 14% at the first stage: e's share, diluted at ε = 0.5 by p₀'s draws. It rises to 43–59%, 59–85% and 69–97% at b = 0.25, 0.5 and 0.75 over all pilots, and to 68–85%, 88–96% and 97–100% without the SMC pilot.
  - B starts at 13–28% and falls to 2–13% (all pilots) and 0.5–3% (no SMC pilot) at 0.75.
- *Origins.* No final weight descends from p₀'s draws, in any replicate. 75–100% descends from draws in A.
- *The kernel's acceptance.*

  | | births | deaths | redraws |
  | --- | --- | --- | --- |
  | ε = 0.05 | 0.58–0.79% | 0.52–0.70% | 0.63–1.08% |
  | ε = 0.5 | 3.67–4.68% | 4.02–5.23% | 3.98–4.95% |
  | R2's from-prior runs | 8.3–9.0% | 10.8–12.4% | 23.6–25.9% |

  The support rule refuses 14.9–15.9% of births and 4.7–5.3% of redraws, as in R4.

*What the runs show.*
- *Replicate agreement within about a nat is reached at 4,000 particles in all four arms.* ln p(y)'s sd is 0.26–0.58 (ranges 0.59–1.30), and the between-run error 0.13–0.29 sits beside the within-run error 0.16–0.27. At 1,000 particles the sd is 0.74–1.45 (ranges 1.57–3.32). From the prior it was 8.4, 7.1 and 6.5 at 1,000, 4,000 and 16,000.
- *All eight arms read one value,* pooled −886.4 to −887.5. Leaving out the SMC pilot moves the pool by 0.2–0.9 nats; ε moves it by 0.1–0.5.
- *That value lies below every from-prior replicate.* Every reference replicate (−889.15 to −885.57) sits at or below every one of R2's twelve (−885.57 to −855.22).
  - Ẑ is unbiased for the evidence, up to the adaptive schedule's small bias, so by Markov's inequality P(ln Ẑ ≥ ln p(y) + t) ≤ e^{−t} for each replicate.
  - Over R2's twelve replicates, −855.22 then puts ln p(y) above −864.6 at 99.9% (t = ln 12,000 = 9.4). The arms' pools lie 21.8–22.9 nats below that bound.
  - The same holds for basin A alone, since Ẑ times a basin's fraction is unbiased for its mass. Here A's log mass reads −885.7 to −889.6; R2's 16,000-particle replicate 3 read −858.5. The shortfall lies inside A, not between basins.
- *The runs read the mass of e's own structures.*
  - The evidence over e's structures, with no move, is −887.15 (all pilots) and −887.59 (no SMC pilot). The arms' pools read from 0.3 nats below it to 0.8 above.
  - The population holds e's structures through b = 0.75 (100% of the weight). e's height over ε·p₀ there is 95–175 nats, so leaving them costs (1 − b)·h, about 30 nats at b = 0.75.
  - It leaves them only in the last stages, ending with 19–94% of its weight still on them, over 84–1,480 distinct structures. Those neighbours add at most a factor of about two to the mass.
  - On B and C, whose F holds a few thousand structures, the same last stages recover the posterior mass off e. On Friedman they do not.
- *No transition is crossed along this path.* The basins' weights move smoothly, and the families never fall to two: 20–174 at b = 1, against 1–5 from the prior. The stages number 49–78, against 179–199.
- *The basins.* Every replicate holds A most (65–98%), M4's best chain's basin, with two x1·x2 terms at 87–100%. The weighted mean log posterior, −888.8 to −895.1, sits at M4's best chain's median (−892.4). No replicate reaches a1 or a2. The pilots held a1 nowhere, and a2 only in the SMC pilot (0.42% of its weight; 0.21% of e over all pilots, none without the SMC pilot).
- *σ.* Every arm's pooled 95% interval misses 1, its lower end at 1.017–1.022, as from the prior at 1,000. R2's 4,000- and 16,000-particle intervals covered.
- *Strength.*
  - Four replicates per arm.
  - The arms share their seeds, so their agreement with one another is not independent.
  - The bound on ln p(y) rests on R2's replicates being SMC normalizers; its largest, −855.22, is one replicate.

#### Run times

Release, with the four replicates in parallel. Concurrent sessions loaded the machine throughout: the same from-prior SMC at 4,000 particles took 176, 251 and 365 s in this chip's three runs, against R2's 107–127 s. The times are rough.
- The arms at 4,000 particles took 161–296 s, against 176–365 s from the prior in the same sessions.
- The arms at 1,000 took 20–41 s, against 38–56 s.
- Stages number 49–78 against 179–199, but each costs more, since the particles carry M of about 8–9 terms from the start.
- The checks took about 2 minutes; the pilots 4–7 minutes, almost all of it the from-prior SMC.

#### The review passes and their fixes

1. *Correctness, as the brief lists it.*
   - The kernel's acceptance against γ_b, with D01's ratios unchanged: the enumeration check, which sees q (total variation 0.25 and 0.08 from the path from p₀). Every 97th sweep it also asserts that the move's evaluation equals its state's.
   - q's density matching its sampler: the draws check, with normalization to 10⁻⁹.
   - Z_F's correction and its error: exact on B and C, the estimate within 1 SE there, its error carried.
   - The pilots' seeds disjoint from the replicates': asserted.
   - Phase 0's readings unchanged: `BassMixing.cs` is untouched.
2. *Fixes in the build.*
   - `Math.Log1p` does not exist in .NET, so it became `Math.Log(1 − ε)`.
   - An arithmetic seed (4202 + 10ε) was replaced by `SeedTree` children.
   - ε = 0 was admitted, for e's own check.
   - A meaningless assertion was replaced by the evaluation identity above.
   - The families-at-b lookup was rewritten as loops.
   - z prints "—" where the SE is 0.
3. *After the first read of the runs.* Agreement on a value below the from-prior runs pointed to the population being held on the pilots' structures, but that was an inference. The harness gained:
   - the final weight on e's own structures, and the same at each mark;
   - the distinct structures;
   - the mean log posterior;
   - e's height over ε·p₀;
   - e's own evidence on Friedman;
   - on B and C, the exact posterior's mass on e's structures.

   The reading was then refined: the population does leave e's structures, but only after b = 0.75, and what it finds there adds under a nat.
4. *Coverage.*
   - *Readings on one replicate:* none here. The comparison with a1 rests on R2's one replicate, as the measurement said.
   - *Families collapsing anyway:* no, 20–174 remain at b = 1.
   - *p₀'s draws outside F:* 17–18% of p₀'s draws on Friedman. Counted, zero-weighted from the first stage (asserted at every stage), never moved.
   - *Repeatability:* the rerun with the added readings reproduced every ln Ẑ, family count and basin fraction to the printed digit, since the readings draw nothing.

#### For the main session

- *The brief's question, as read:* replicate agreement within a nat is reached at 4,000 particles in every arm. It is reached on a value 22–23 nats below the 99.9% lower bound that R2's own from-prior replicates put on ln p(y). The value is the evidence of the pilots' structures, read within 0.8 nats. So the agreement measures how well the reference covers the posterior, not ln p(y). This is a reading; no route is recommended.
- *`BassTarget`'s structure prior is the count weight w,* normalized neither by C nor on F. An evidence from any start other than `BassPriorSampler`'s draws must subtract ln C + ln Z_F; from those draws both cancel, as in C4's runs.
- *No production seam was needed;* the kernel uses existing internals only.
- *The ridge chip committed during this chip* (`e225c22`). The harness built against it unchanged, and the from-prior pilots reproduced R2 to the printed digit.
- *Another session's staged files* (`src/maths/samplers/flathistogram/Reweighting.cs` and three sampler tests) were in the index at commit time. `abdcd4d` and this report's commit were made by pathspec, so those files remain staged as that session left them.
- `artifacts\chip-reference\` is deleted.
