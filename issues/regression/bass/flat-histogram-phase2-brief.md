# The flat-histogram sampler, Phase 2: BASS's adapter and the Friedman qualification (chip brief)

Date: 2026-10-30. Builds Phase 2 of N17 on Phase 1's engine, which was built and verified on 2026-10-30. The design's decisions are settled ([flat-histogram-design.md](flat-histogram-design.md), including its "Phase 1, as built"). This brief fixes the build.
- *What it builds:*
  - BASS's adapter, `BassSampling.FlatHistogram`, null by default: a multicanonical fit whose draws are weighted at b = 1, with its evidence and its diagnostics;
  - the enumeration checks on fixtures A, B and C;
  - the Friedman qualification, which reads and sets no default.
- *Not in scope:*
  - any default change: MCMC stays the default, and every configuration without the new spec draws bit-identically;
  - changes to `Maths.Samplers.FlatHistogram`, SMC or `ParallelTempering`, unless one is needed (then name it in the report before making it);
  - SPC's adapter (Phase 3);
  - `docs/Maths` and the checklists, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the design, whole;
- [Phase 1's brief](flat-histogram-phase1-brief.md), its reports above all, their "For Phase 2" lists included;
- `docs/Maths/Samplers.md`, "Flat-histogram sampling (N17 Phase 1)";
- the Friedman readings this qualification answers to: [the band](friedman-band-brief.md#reports) and [the reference path](friedman-reference-brief.md#reports). In `docs/Maths/BayesianRegression.md` they are under "BASS: SMC", "Friedman's band, measured" and "A path from a pilot reference";
- D05 in [port-fitting-brief.md](port-fitting-brief.md) (the linear and integrated forms, the reference rung) and D10 in [port-numerics-brief.md](port-numerics-brief.md) (adaptation in burn-in only, the verdict).

In the code:
- *The engine,* `src/maths/samplers/flathistogram/`:
  - `MulticanonicalTarget`, `MulticanonicalRun`;
  - `IMulticanonicalModel` and `IMulticanonicalChain`, with `MulticanonicalChain` as the reversible-jump example;
  - `Reweighting`, `TemperingPilot.FromLadders`, `Wham` and `FlatHistogramSpec`.
- *BASS:*
  - `Bass.cs`, whose `Fit` dispatches to SMC on `BassSampling.Smc`, the pattern for the new branch;
  - `BassSmcRun.cs` (`BassSmcFit.Run`, `Assemble`): how a weighted fit is put together;
  - `BassFit.cs` (`DrawWeights`, which every view reads) and `BassEvidence.cs`;
  - `BassSweep.cs`: `BassChain` takes any `ISamplingTarget<BassState>`, so it runs on the wrapper. `BassTauBlock` and the structure moves are the kernels;
  - `BassLadder.cs`: `TraceRungs(traceSlopes: true)` under `TemperedForm.Linear`, which the pilot reads;
  - `BassConfig.cs` (`BassSampling`, `TemperingSpec`, validation) and `BassRunLengths`.
- *The harnesses this one draws on:*
  - `tests/regression/BassMixing.cs` (`BassMixing.Get(BassMixing.Friedman)`, `Key`);
  - `BassFriedmanBandTests.cs` (the basins' `Legend`);
  - `BassFriedmanReferenceTests.cs` (the pilots, Z_F, the evidence of e's own structures);
  - `BassSweepTests.EnumerateF` and the enumeration oracle's fixtures A, B and C.

## The design, as Phase 2 builds it

These are the brief's shape calls. Revise one only with a reason, and name it in the report.

**The configuration.** `BassSampling.FlatHistogram`, a values-only record, null by default, serialized in `BassConfig`'s JSON:
- the engine's `FlatHistogramSpec`;
- for a WHAM start, the pilot: a `TemperingSpec` and its lengths (chains, burn-in, traced scans). The pilot's form must be `TemperedForm.Linear`, since WHAM needs rung laws prior·e^{b·ℓ₁}; validation refuses the integrated form there. The reference rung is on by default, so the hottest rung is b = 0 and b_lo = 0.
- *Validation refuses:*
  - a flat-histogram spec together with `BassSampling.Smc`;
  - an improper τ prior (hyper-g with a ≤ 2) wherever b_lo = 0, since the prior over (structure, τ) is then improper, as a reference rung and SMC refuse it.
- *Under Jeffreys σ² or the flat intercept* the prior over (structure, τ) stays proper (D05), so the b_lo = 0 window runs, but the fit reports no evidence.
- `BassRunLengths` reports the spec, as it reports SMC's.

**The chain.** A BASS `IMulticanonicalChain` over `BassChain` on the `MulticanonicalTarget<BassState>`:
- the structure moves and τ's slice step run on π_w unchanged, since they read the target;
- `Reweigh` re-reads ln w from the current evaluation's cached ℓ₁, with no new marginal;
- `SaveState` and restore go through `BassChain`'s snapshots.
- ℓ₁ is the collapsed log marginal at b = 1: with τ sampled, the state is (structure, u), and ℓ₁ = ln p(y | structure, τ).

**The fit.** The production chains' draws, composed as an MCMC fit's are, so each keeps its full state (D04, retain-don't-discard).
- `DrawWeights` are the b = 1 weights e^{ℓ₁}/w(ℓ₁), self-normalized, in `DrawStore`'s pooled order, so every view reads them as it reads an SMC fit's.
- `Traces()` refuses, as for SMC: the draws' order is π_w's, not the posterior's. The flat-histogram diagnostics take the traces' place.
- *The diagnostics* carry:
  - the engine's occupancy, flatness and round trips per chain;
  - the learning's weight trace;
  - the b = 1 Kish ESS, pooled and per chain;
  - D10's verdict on π_w's functionals;
  - the window, K, b_lo and t₀ as resolved.

**The evidence.** `LogEvidenceCurve` on a grid of b, with ln Z(1) as the fit's evidence and its error from the chains' spread.
- *It needs no correction for the structure prior:* the ratio estimator normalizes against the prior itself (Phase 1's correction), so it gives ln p(y) under the prior on F directly.
- *With b_lo > 0* the denominator reaches the prior only through the b_lo-tempered law below the window, a one-step importance ratio. So the evidence is reported only when b_lo = 0. With b_lo > 0 the fit reports ln Z(1) − ln Z(b_lo) (`LogRatio`), flagged as a ratio.
- `BassEvidencePriors` records the spec, as for SMC.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Other sessions may be working in this checkout,** an SPC session among them. If a build fails in a file you do not own, wait and retry; never edit it.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-fh2\` to every `dotnet build` and `dotnet test`. Build in Release for the qualification. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-fh2\` when you finish.
- **Your files:**
  - in `src/maths/regression/spline/bass/`: `BassConfig.cs`, `Bass.cs`, `BassFit.cs`, `BassEvidence.cs`, `BassDiagnostics.cs`, `BassRunLengths`'s file, the presentation files that must refuse or report the new fit, and new files for the adapter;
  - tests under `tests/regression/`, with additions to `BassMixing.cs` that leave Phase 0's readings unchanged;
  - this brief's Reports section.
- **Commit your own files only, by explicit path,** one commit per step, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked notes under `issues/regression/` and other sessions' edits stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.
- **Bit-identity:** `BassBitIdentityTests`' hashes hold, and every configuration without the new spec draws exactly as before.

## The steps

Each step lands with its checks and ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next step. Name the passes and their fixes in the report.

1. **The configuration and the chain.**
   - *Check:*
     - the JSON round trip, and every refusal (the integrated pilot, SMC together, an improper τ prior where b_lo = 0);
     - bit-identity without the spec;
     - the chain on the wrapper reproduces π_w on fixture A's enumerable F at fixed weights, against enumeration, with held τ and with sampled τ by quadrature over u;
     - `Reweigh` makes no new marginal.
2. **The fit, its evidence and its views.**
   - *Check on fixture A* (enumerable F), under τ held and sampled, and under both coefficient forms:
     - the exact density of states, from each structure's ℓ₁ and prior mass (with sampled τ, quadrature cells over u), so `FixedPoint`'s θ* is known;
     - the learned weights against θ*;
     - the reweighted structure posterior (the count law, each term's inclusion, the heaviest structures, E[τ]) against the oracle's exact posterior, within 5 SE of the chains' spread;
     - ln Z(b) on a grid against enumeration.
   - *Check on fixtures B and C:* ln Z(1) against the oracle's exact −28.589 and −30.903, within 4 errors + 0.01 (C4's tolerance).
   - *End to end:* a scalar fit and a functional fit, every view running on the weighted draws: summary, metadata, both pointwise contents, prediction, sensitivity, effects and the residual check. `Traces()` refuses.
   - *Resume,* bit for bit, across learning and production.
3. **The warm start on BASS.**
   - *Build:* the pilot as a linear-form, reference-rung, self-sized ladder, traced, read by `TemperingPilot.FromLadders`.
   - *Check on fixture A:* WHAM's ln ĝ against the exact density of states within the pilot's errors, and the warm start reaching the flat start's accuracy in fewer iterations (report both).
4. **The Friedman qualification** (`BASS_QUALIFICATION` = 1; S10's Friedman fixture, n = 500, from master seeds no pilot or earlier harness uses).
   - *Write down before running* what would count, so the reading is not fitted to its numbers:
     - *the chains tunnel:* every production chain makes round trips between the window's ends;
     - *they agree:* the chains' ln Z(1) lie within about a nat of one another, and their reweighted basin fractions agree within their spread;
     - *they reach the bound:* the pooled ln Z(1) sits at or above −864.6, the 99.9% lower bound that R2's from-prior normalizers set by Markov's inequality, against the −887 of the pilots' own structures.
   - *The arms:*
     - (i) the engine's window from b_lo = 0, the reference rung, to b = 1: the full evidence;
     - (ii) a window from b_lo just below the lower split (b ≈ 0.12–0.15), to b = 1, its pilot a ladder with its hot end at that b_lo and no reference rung: the posterior and ln Z(1) − ln Z(b_lo) only, narrower and cheaper. If fitting SMC from the prior can stop at that b_lo without a production change, add its ln Z(b_lo) to read ln Z(1) through it. Otherwise name the seam and leave it.
     - Each arm from the WHAM warm start at the default gain start. Arm (i) also from the same warm start with the gain counted from 0, since the default's trust in the pilot is regime-dependent.
   - *Size the runs* for at least 20 round trips per production chain where the chains tunnel at all, within about 30 minutes per arm. Report the budget each took.
   - *Read, per arm and per chain:*
     - the window, K, b_lo, t₀ and the gain's start;
     - whether `FixedPoint` accepted the pilot's ĝ, and where ĝ had no mass, since the pilot does not cross the band at b ≈ 0.23–0.28;
     - the learning's trace, and whether the weights settled;
     - production's occupancy and flatness, and the round trips;
     - the b = 1 Kish ESS;
     - ln Z(1) (or the ratio), pooled with its error;
     - the reweighted basin fractions under the band measurement's structural key (A, B, G, S, T, a1, a2, a5), each with its log mass, and the share by the count of x1·x2 terms;
     - σ's 95% interval;
     - D10's verdict on π_w's functionals.
   - *Against,* at matched or stated cost, from their recorded readings: R2's from-prior SMC (1,000, 4,000 and 16,000 particles), the reference arms (−886.4 to −887.5), and M4's MCMC chains.
   - Name whether any chain visits a1 or a2, the four-x1·x2 basins that only R2's 16,000-particle replicate found.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - the moves' acceptance on π_w, with D01's proposal ratios unchanged and the structure prior read from the inner target;
  - τ's slice on p(u)·w(ℓ₁(S, u));
  - `Reweigh` after every installed weight set;
  - the weights frozen in production;
  - `DrawWeights` aligned with `DrawStore`'s order;
  - the evidence only where b_lo = 0;
  - nothing existing changed in its draws.
- *Coverage:*
  - held and sampled τ;
  - both coefficient forms;
  - Jeffreys σ² and the flat intercept, sampled from b_lo = 0 with no evidence reported;
  - an improper τ prior, refused at b_lo = 0;
  - a functional response;
  - a chain that never tunnels, which the diagnostics flag;
  - resume mid-learning;
  - a pilot whose ĝ has a gap.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the API and the configuration, with any shape call revised;
- each check with its tolerance;
- the qualification: what would count as written before the run, then the tables, then what they show, kept apart;
- the times;
- the suites with their counts (Regression, Samplers, the BASS oracle classes);
- the review passes and their fixes;
- anything the main session must know.

Recommend no default.

## Reports

Newest first.
