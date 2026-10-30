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
