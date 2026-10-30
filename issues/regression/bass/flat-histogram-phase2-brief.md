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

### Verified by the main session, with a second seed — 2026-10-30

- *The build:*
  - the adapter read (`BassMulticanonical.cs`, `BassFlatHistogramRun.cs`): the chain on π_w, `Reweigh`, the weights against the stores, the evidence's conditions and the pilot's linear form;
  - the suites rerun in Release: Regression 1,331, Samplers 262, the BASS oracle classes 90, all passing;
  - the nine nullable warnings cleared (`93ffbef`).
- *The below-split reading* checked. At b ≤ 0.03 the curve, SMC and exact importance sampling agree to 0.01 nats, and an error in the denominator would shift every b alike, so the 5.5-nat gap at b = 0.10 is SMC's.
- *A second seed for arm (i),* written before its run (`56c66f3`): master seed 8217, the sweeps doubled through `Thin` = 400, and the same criteria, plus agreement with seed 8017 (ln Ẑ(1) within a nat plus both errors of −841.36, a1 within 5 points of 65.9%).
  - *The run was refused at its start* (48 s). SAMC's fixed point could not flatten the pilot's ĝ over [−1552.56, −773.559], K = 128, b_lo = 0: no node was empty, but a node's occupancy stood 125 from flat.
  - *A rerun from a flat start* over that window, K and b_lo, with the gain from 0 and the same lengths, was written before it ran (`699bd8b`). It took 22 minutes, at 67,000 sweeps a second a chain.
    - Learning: no round trip on any chain in 100,000 iterations, flatness 0.99, the weights still moving 1.5 nats (RMS) between 10% and 50%.
    - Production: no round trip, flatness 0.98, no draw above the window, b = 1 Kish ESS 8 of 500,000.
    - ln Ẑ(1) −843.64 ± 0.46, chains −842.9 to −844.9; σ's interval [1.012, 1.104].
    - a1 at 74.7% pooled, the chains at 81, 79, 47 and 46%, failing the χ² (p < 0.001).
    - Its ln Z(b) curve matches seed 8017's to within a nat up to b = 0.8 (−126.1, −531.0, −687.7 at 0.1, 0.6, 0.8).
  - *Against what would count:* the first criterion fails, so the rerun reads nothing about the posterior. What it shows agrees with seed 8017 without confirming it.
- *So:* the Friedman reading stands on one seed. Replication needs a warm start that holds across seeds, whether a fixed point that does not refuse a ĝ with no empty node, or a learning budget a flat start can cross the window in.

### Phase 2 built, and the Friedman qualification read (chip, 2026-10-30)

**What landed**, on `main`, nothing pushed:
- `33e067e`, step 1: the configuration (`BassSampling.FlatHistogram`, `BassFlatHistogramSpec`, `BassFlatHistogramPilot`, validation, JSON) and BASS's chain on π_w (`BassMulticanonical.cs`: `BassMulticanonicalModel`, `BassMulticanonicalChain`). `Bass.Fit` refused the spec until step 2.
- `9d78484`, step 2: the fit (`BassFlatHistogramRun.cs`: `BassFlatHistogramFit`), its evidence (`BassFlatHistogramEvidence`; `BassEvidence.FlatHistogram`), its readout (`BassFlatHistogramReadout`, `BassDiagnostics.FlatHistogram`), `Traces()` refused, the incomplete fit refused by the posterior views, and `BassRunLengths.FlatHistogram`.
- `820982b`, step 3: the pilot and its warm start (`BassFlatHistogramFit.Pilot`, `Pooled`, `BassFlatHistogramPilotReadout`), with `BassFlatHistogramSweep.cs` registering the warm-against-flat claim.
- `68ee777`: an iteration reads `Sampling.Thin` (below, a shape addition).
- `6a656d8`, step 4: the qualification's harness, `BassFlatHistogramQualificationTests.cs`, with what would count written into it before the run.
- `acf9d2d`: the reading below the lower split, added to the harness after the run to locate the arms' difference.
- `78708b3`: `BassMovePassTests`' data hash re-pinned (the suites, below).
- No change to `Maths.Samplers.FlatHistogram`, SMC or `ParallelTempering`.

**The API and the configuration.**
- `BassSampling.FlatHistogram`, a `BassFlatHistogramSpec`, null by default and serialized in the config's JSON (`"FlatHistogram":null` beside `"Smc":null`), so every fingerprint moved and no draw did:
  - `Run`, the engine's `FlatHistogramSpec`;
  - `Pilot`, a `BassFlatHistogramPilot`: `Tempering` (by default the linear form with the reference rung, self-sizing), `Chains` (4), `BurnIn` and `TracedScans` (both required, no default chosen);
  - `MaxIterations`, a per-call budget (below).
- *Validation refuses:* the spec beside `Sampling.Smc`, beside a top-level `Sampling.Tempering`, or with an R̂ target; a pilot missing where the engine's spec needs one (`NeedsPilot`), or present where nothing reads it; the pilot's integrated form or powered path; the pilot's lengths; the engine's own fields, named under `Sampling.FlatHistogram.Run`; and an improper τ prior wherever b_lo = 0 (a declared `LowSlope` of 0, or b_lo left to a reference-rung pilot). The tempering checks moved into one helper shared by `Sampling.Tempering` and the pilot, their messages unchanged.
- `Bass.Fit` dispatches the spec to `BassFlatHistogramFit.Run`; `Bass.Resume` continues an incomplete run.
- *The chain* (`BassMulticanonicalChain`): `BassChain` on `MulticanonicalTarget<BassState>` over `model.Target`, so the structure moves and τ's slice step read π_w through the target. `Reweigh` is `SetState(Current, evaluation.Reweighed(target.Weights))`; `SaveState` is the `BassChain` snapshot plus the composition stream, term table, draws, start and the move counts when production began. A production draw is composed as an MCMC draw is, and reduced to M, τ, the log posterior (the inner prior plus the cached ℓ₁), σ² and the fitted values, in D10's three groups.
- *The fit:* one `DrawStore` per production chain, terms merged in chain order; `DrawWeights` = `Reweighting.Weigh(draws, 1).Weights`, chain-major in draw order.
- *The readout* (`BassFlatHistogramReadout`): the weights as resolved (window, K, b_lo), t₀ and the gain's start, the learning's iterations, visits and weight trace, production's visits pooled and per chain, the b = 1 Kish ESS pooled and per chain, the production draws (`MulticanonicalDraws`, for reweighting to any b), the evidence curve, the pilot's readout and `Complete`. The convergence diagnostics and `Verdict` read π_w's functionals over production (D10's rules on the chains' own target); `ModeSeparation` reads the chains' reweighted means at b = 1; `MassAtCap` is weighted; the accounts cover production.
- *The evidence* (`BassFlatHistogramEvidence`): `LogEvidenceCurve` on 101 points from b_lo to 1, `LogRatio` = ln Ẑ(1) − ln Ẑ(b_lo) with its chains'-spread error, `IsEvidence` at b_lo = 0. Computed for a complete run under a proper σ² prior and intercept when asked for; `BassFit.Evidence` holds it (with `BassEvidencePriors`) only from b_lo = 0 under a proper τ prior.
- *The pilot* (`BassFlatHistogramFit.Pilot`): `Chains` ladders on the linear form, each from an exact draw of the prior on F, burned as a tempered fit burns them (adapted after rounds of 1, 2, 4, … scans, sized after every round but the last), frozen and traced. `BassFlatHistogramPilotReadout`: the ladders' `BassTemperingReadout` over the traced scans, the pooled rungs and their inefficiencies, the draws, and the pilot's draws under each node's hat (`NodeDraws`, `EmptyNodes`).

**Shape calls revised or added**, each with its reason:
1. *`MaxIterations`, a per-call budget* (added). The engine's lengths are totals, and `Bass.Resume` needs a call that can stop short of them, across learning and production. It mirrors `SmcSpec.MaxStages`: a fit short of its lengths is flagged (`Complete` false), reports no evidence or residual check, refuses the posterior views, and its checkpoint continues it.
2. *A complete run carries no checkpoint* (revised from SMC's pattern). The engine refuses a complete checkpoint, and its production records are private to it, so a complete run's resume cannot give the fit again as SMC's does. `BassFit.Checkpoint` is null once the lengths are reached.
3. *`Sampling.Thin` is read* (added, `68ee777`): each iteration, a learning update or a production draw, follows `Thin` sweeps of every chain, MCMC's meaning. Production keeps every draw it takes; Friedman needed 25M production sweeps a chain, 100M draws at one per sweep, beyond memory.
4. *The pilot's ladders are read one at a time* (revised from "read by `TemperingPilot.FromLadders`"). `FromLadders` needs every ladder on the same rungs, and self-sizing ladders each freeze their own. Each ladder is read by `FromLadders`, and the rungs pool by b: a rung several ladders share (b = 1, b = 0) holds their draws, its inefficiency Σn / Σ(n/g), their draws over their summed effective draws.
5. *Above b_lo = 0 the ratio sits in the readout, not in `BassFit.Evidence`* (revised from "the fit reports … LogRatio, flagged as a ratio"). `BassEvidence.LogEvidence` is ln p(y) for every other sampler; a ratio there would read as one. `BassFit.Evidence` is null, and `Diagnostics.FlatHistogram.Evidence` holds the ratio with `IsEvidence` false.
6. *`BassEvidence.SteppingStone` is now nullable*: a flat-histogram fit's evidence has no stepping-stone estimate, and `LogEvidence`/`MonteCarloError` read its `FlatHistogram` estimate instead. Nine reads of `evidence.SteppingStone.…` on ladder and SMC fits, in seven harnesses this chip does not own, now compile with nullable warnings (listed for the main session below).
7. *A WHAM start the fixed point refuses throws, naming the window, K, b_lo and the empty nodes*, so a qualification reads the refusal.
8. *Production's per-chain visits* are computed by the adapter from each chain's levels, since the engine merges them; the round trips are the engine's.

**The checks, with their tolerances** (`BassFlatHistogramTests`). Fixture A is `BassSweepTests`' F (N = 18, M ≤ 3, 2,263 structures); `BassFlatHistogramSpace` enumerates it with ℓ₁ on an 801-node trapezoid grid in u, which reproduces the oracle's ln Z on B and C to 0.0004.
- *Step 1.*
  - The JSON round trip, equal and stable, `NeedsPilot` kept out; an unknown field and a missing required one refused; the fingerprint moves with the spec.
  - Every refusal above, each naming its field; improper τ accepted from b_lo > 0; Jeffreys σ² and the flat intercept accepted from b_lo = 0.
  - The chain on π_w at fixed weights (SAMC's limit on the exact density of states, K = 16, b_lo = 0.2), 4 × 30,000 sweeps, against π_w by enumeration within 5 batch-means SE + a slack (0.004 on M's law, 0.005 on inclusions and the ten heaviest structures, 0.001 on E[u], 0.01 on E[ℓ₁]):
    - τ held: 24 rows, worst |z| 4.47 on P(M = 0) = 0.0034 (read 0.0020 ± 0.0003), inside its slack. A probe at 120,000 sweeps on three seeds read 0.0029–0.0034 for every chain: the batch SE under-reads a rare state's bursts. Removed after reading.
    - τ sampled: 40 rows, worst |z| 2.44; E[u] −4.4588 ± 0.0173 against −4.4628.
    - Every 97th sweep, the chain's density equals the inner prior plus ln w(ℓ₁), and its cached ℓ₁ the state's own marginal.
  - `Reweigh` under moved weights: no inner evaluation, the same inner evaluation object, ln w re-read from the cached ℓ₁, the chain's density updated; under the same weights the evaluation is kept.
  - Bit-identity: `BassBitIdentityTests`' seven hashes hold, run at steps 1 and 2, with the `Thin` change and in the final suite.
- *Step 2*, on fixture A, 8 chains × 30,000 learning + 30,000 production iterations from flat, K = 12 over the prior's 1% to the posterior's 99% quantile of ℓ₁, t₀ = 120:

  | τ, form, b_lo | learned ln w to θ* (≤ 0.15) | production flatness | round trips | b = 1 Kish ESS of 240,000 | posterior rows, worst \|z\| | ln Z(b) − ln Z(b_lo), worst \|z\| of 5 |
  | --- | --- | --- | --- | --- | --- | --- |
  | held, g, 0 | 0.081 | 0.069 | 1,781 | 15,840 | 29, 4.49 | 1.03 |
  | sampled, g, 0 | 0.034 | 0.084 | 2,693 | 26,415 | 26, 3.33 | 0.63 |
  | held, ridge, 0 | 0.100 | 0.150 | 651 | 25,272 | 35, 4.09 | 0.20 |
  | sampled, ridge, 0 | 0.029 | 0.062 | 2,830 | 23,814 | 27, 2.87 | 1.03 |
  | sampled, g, 0.3 | 0.032 | 0.069 | 2,660 | 26,598 | 26, 3.70 | 1.07 |

  - *The posterior rows* (M's law, each term's inclusion at 0.5% or more, the ten heaviest structures, E[τ] when sampled) are the reweighted estimates against the exact posterior, within 5 SE of the chains' spread + 0.002 (0.5% of E[τ]). Every row beyond 3 SE is a probability under 3% (4.49 on an inclusion of 0.026, 4.09 on a structure of 0.007), where eight chains' spread under-reads; each sat inside its slack. With τ held every draw holds the held τ.
  - *ln Z(b)* at b = 0.1, 0.25, 0.5, 0.75, 1 of [b_lo, 1], against the enumeration within 5 errors + 0.01. At b_lo = 0 the curve's first value is 0 exactly, and `BassFit.Evidence` is ln Ẑ(1) with no stepping-stone; at b_lo = 0.3 `BassFit.Evidence` is null and the readout holds the ratio, flagged.
  - *Fixtures B and C*, 8 chains × 30,000 + 30,000 from flat, K = 12: ln Ẑ(1) −28.529 ± 0.047 against the oracle's −28.589 (z 1.26; the enumeration −28.5886) and −30.863 ± 0.054 against −30.903 (z 0.74; −30.9026), within C4's 4 errors + 0.01.
  - *End to end*, a scalar fit (N = 75) and a direct functional fit (16 cases × 5 locations), 4 chains × 4,000 + 1,000 from a rough window: the weights sum to 1 and equal `Reweighting.Weigh`'s; each sampled draw's weight, its ℓ₁ rebuilt from its structure and τ, is e^{ℓ₁}/w(ℓ₁) to one normalizer (8 decimals); summary (σ²'s and M's means against the weighted draws by hand to 1e−12), metadata (the spec, `Thin`, no burn-in or samples), both pointwise contents (carrying the weights), prediction, sensitivity (its decomposition's weights), effects, the training fit (its mean by hand to 1e−9) and, on the functional fit, the residual check, all run; `Traces()` refuses; the verdict has its three groups; a complete fit has no checkpoint.
  - *Resume:* at two sweeps an iteration, a budget of 170 stops mid-learning (no draws, every posterior view and `Traces()` refused, the metadata given), then 40 draws into production, then completes; bit for bit with the uninterrupted run (draws, accounts and move counts by `BassBitIdentityTests.Hash`, weights, the learned weights and their trace, both phases' occupancy, round trips and end states, the evidence, R̂). The accounts cover production's sweeps only. Another config's checkpoint and an MCMC one are refused.
  - *The improper options* from b_lo = 0: Jeffreys σ² and the flat intercept run to completion with no evidence, in the fit or the readout.
  - *A chain that never tunnels:* three production draws a chain; both chains flagged in `ChainsWithoutRoundTrip`.
- *Step 3*, on fixture A, τ sampled, K = 12 over the exact window from b_lo = 0:
  - *The pilot*, two ladders burned 255 scans and traced 2,000: each sized itself to 5 rungs (1, 0.325, 0.184, 0.106, 0 and 1, 0.368, 0.183, 0.104, 0), pooled to 8 rungs with b = 1 and 0 shared, 20,000 draws, inefficiencies 3.0–5.4, 0.073 round trips per scan; WHAM's residual under 1e−8. The pooling unit test reads Σn / Σ(n/g) on shared rungs to 1e−15.
  - *WHAM against the exact density of states*, eight pilots' starts against SAMC's limit θ*: per-pilot centred RMS 0.097 (tolerance 0.25), mean z² over the nodes 0.80 (tolerance 3; t₇ gives 1.4), Phase 1's measure.
  - *Warm against flat*, the mean of three learning seeds, four chains, the warm start's gain from the pilot's draws over the chains (5,000):

    | iterations | 0 | 2k | 4k | 8k | 16k | 32k |
    | --- | --- | --- | --- | --- | --- | --- |
    | flat | — | 0.205 | 0.124 | 0.086 | 0.059 | 0.029 |
    | warm, default gain start | 0.091 | 0.110 | 0.072 | 0.078 | 0.046 | 0.050 |
    | warm, the gain from 0 | 0.091 | 0.157 | 0.105 | 0.094 | 0.050 | 0.045 |

    The pilot alone is worth about 8,000 flat iterations, and the warm start reaches the flat start's accuracy in fewer iterations up to 8,000 (warm 4,000 0.072 against flat 8,000 0.086, asserted and registered). Past 16,000 flat learning overtakes it: from the default gain start the warm start's learning plateaus near 0.05, the pilot's residual error not worked off at a gain counted from 5,000, as Phase 1 found for the double well. The length of the asserted comparison was chosen after an exploratory reading of this table, said here; the 16k and 32k columns are that exploration's (three seeds, the same pilot).
  - *A pilot whose ĝ has a gap*: three fixed rungs (1, 0.9 and the reference rung's 0), τ held so ℓ₁ is discrete, K = 20. On seed 7803 the pilot drew nothing under nodes 13–15, and the start reads them off the line through nodes 12 and 16 (to 1e−9); the fit completes with its evidence. On seed 7802 SAMC's fixed point refused the ĝ (a node's occupancy 13.9 past flat, nodes 13 and 14 empty), and the fit throws naming the window, K, b_lo and those nodes.

**The qualification** (`BassFlatHistogramQualificationTests`, behind `BASS_QUALIFICATION` = 1; S10's Friedman fixture, n = 500; Release).

*What would count, as written into the harness and committed (`6a656d8`) before the run:*
1. *The chains tunnel:* every production chain completes round trips between the window's end bins, the runs sized for at least 20 per chain; a chain with none has not crossed, and its arm's agreement and evidence then read nothing about Friedman's posterior.
2. *They agree:* the chains' own ln Ẑ_c(1) (arm i; arm ii's ratios) within about a nat, read as an SD over chains of at most 1 nat; and for each labelled basin holding 5% or more of the pooled weight, the chains' reweighted fractions pass a χ² test of equality against their own batch-means errors (20 batches per chain) at p ≥ 0.01.
3. *They reach the bound:* the pooled ln Ẑ(1) at or above −864.6, the 99.9% lower bound R2's from-prior normalizers put on ln p(y), against the −887 of the pilots' own structures; arm (ii)'s ln Ẑ(1) through SMC's ln Ẑ(0.10), errors in quadrature.
4. Agreement below the bound would repeat the reference path's reading; agreement at or above it with every chain tunnelling would be the first sampler that computes Friedman's structural posterior.

*The arms*, as fixed in the harness before the run:
- (i) b_lo = 0: the pilot's reference rung, master seed 8017, from the WHAM start at the default gain start and, on the same start, with the gain from 0;
- (ii) b_lo = 0.10, just below the lower split (b ≈ 0.12–0.15 on the band measurement's ladders), its pilot's hot end there with no reference rung, master seed 8027. Fitting SMC from the prior stops at b = 0.10 with no production change: the harness drives the SMC core on the linear form (`SmcRun` with `BassParticleTarget`, `BassParticles`, `BassSmcKernel`, C4's defaults, the next b capped at 0.10), 4,000 particles × 4 replicates from master seed 8039, and adds its ln Ẑ(0.10) to arm (ii)'s ratio.
- Each pilot: four linear-form, self-sizing ladders, 5,000 scans of burn-in and 2,000 traced (R1's lengths).
- *The lengths*, the same for every arm: 4 chains × (100,000 learning + 125,000 production iterations) × 200 sweeps, 45M sweeps a chain, every 200th kept (125,000 draws a chain).
- *Their sizing*, on master seed 8101, which no arm uses: 20,000 iterations at one sweep made no round trip; 100,000 at ten sweeps (1M a chain) made one or two per chain in learning and none in 200,000 production sweeps, at about 30,000 sweeps a second a chain. Half an hour an arm was then about 45M sweeps a chain. The sizing's own readings: the pilot's ĝ accepted, ln Ẑ(1) −847.7 ± 4.0 (default gain) and −851.5 ± 3.3 (from 0), the chains 3–17 nats apart.

*The arms' tables.*

| | (i), default gain start | (i), the gain from 0 | (ii), default gain start |
| --- | --- | --- | --- |
| window (width), K, b_lo, t₀ | [−1551.0, −778.3] (773), 204, 0, 2,040 | the same | [−1223.4, −776.0] (447), 88, 0.10, 880 |
| the gain counted from | 42,000 | 0 | 32,500 |
| pilot: rungs per ladder; pooled; draws | 21, 21, 21, 21; 78; 168,000 | the same pilot | 16, 16, 17, 16; 59; 130,000 |
| pilot: round trips per traced scan; barrier; worst pair rejection; g | 0; 8.30; 0.975; 1.1–1,345 | | 0; 5.80; 0.998; 1.9–1,449 |
| SAMC's fixed point on ĝ | accepted, no node empty; under 10 draws at nodes 0–2 (ℓ₁ −1551 to −1543) and 153–155 (−969 to −961) | | accepted, no node empty; under 10 at nodes 0–4 (−1223 to −1203) |
| learning: round trips per chain; flatness; visits above the window | 5, 9, 5, 4; 2.80; 29% | 9, 11, 8, 10; 0.15; 8% | 10, 11, 8, 10; 2.41; 15% |
| the weights' centred RMS change, start–10%, 10–50%, 50–90%, 90–100% | 0.88, 3.42, 0.15, 0.093 | 4.11, 0.24, 0.21, 0.136 | 3.26, 5.26, 0.42, 0.069 |
| production: round trips per chain (all) | 17, 9, 12, 9 (47) | 7, 9, 4, 7 (27) | 16, 17, 13, 13 (59) |
| production flatness, pooled (per chain) | 0.39 (0.43, 0.42, 0.84, 0.50) | 0.41 (0.52, 0.44, 0.85, 0.52) | 0.71 (0.69–0.73) |
| production above the window, per chain of 125,000 | 8,500; 4,889; 8,970; 9,726 | 10,124; 16,377; 1,516; 7,454 | 8,431; 9,538; 7,191; 8,673 |
| b = 1 Kish ESS of 500,000 (per chain) | 34,256 (9,186; 5,242; 9,487; 10,341) | 37,563 (10,754; 17,204; 1,656; 7,949) | 39,702 (9,919; 10,938; 8,736; 10,111) |
| per chain ln Ẑ_c(1), or the ratio | −841.30, −841.90, −841.20, −841.19 | −840.84, −840.21, −842.77, −841.19 | −715.13, −715.15, −715.41, −715.33 |
| their SD (range) | 0.34 (0.71) | 1.09 (2.56) | 0.14 (0.28) |
| pooled | ln Ẑ(1) −841.36 ± 0.17 | ln Ẑ(1) −840.96 ± 0.54 | ln Ẑ(1) − ln Ẑ(0.10) −715.25 ± 0.07; through SMC's −131.35 ± 0.21, ln Ẑ(1) −846.60 ± 0.22 |
| ln Ẑ(b) at 0.1, 0.2, 0.25, 0.3, 0.4, 0.6, 0.8 | −125.9, −210.3, −251.4, −292.0, −372.3, −530.8, −686.7 | −126.0, −210.3, −251.4, −291.9, −372.2, −530.5, −686.3 | ratio only |
| σ's 95% interval | [0.9757, 1.1139] | [0.9768, 1.1129] | [0.9767, 1.1136] |
| D10's verdict on π_w's functionals | converged, every group | not converged: σ² and the log posterior at R̂ 1.01; structural converged | converged, every group |
| time: pilot, learning, production, all | 67, 337, 375, 779 s | 38, 335, 485, 858 s | 33, 463, 400, 896 s |
| sweeps a second a chain | 63,000 | 55,000 | 52,000 |

*SMC from the prior to b = 0.10* (linear form, 4,000 particles): ln Ẑ(0.10) −131.30 ± 0.31, −131.70 ± 0.31, −131.11 ± 0.56, −131.27 ± 0.45; pooled −131.35 ± 0.21 (within 0.21, between 0.13); 63–66 stages; 178–210 families at the end; 42 s.

*The basins*, reweighted at b = 1 under the band measurement's key (labels A, B, G, S, T from M4's chains, a1, a2, a5 from R2's SMC), each with its pooled fraction, its log mass ln Ẑ(1) + ln f (arm ii through SMC's ln Ẑ(0.10)), each chain's fraction, and the χ² p of their agreement where the pooled fraction is 5% or more:

| basin | (i), default | (i), from 0 | (ii) |
| --- | --- | --- | --- |
| a1, 4 × x1·x2 + 2 × x3 + x4 + x5 | 65.9% (−841.8); 61.5, 67.0, 68.2, 67.1%; p 0.917 | 67.6% (−841.4); 66.7, 66.7, 68.9, 70.5%; p 0.976 | 67.7% (−847.0); 67.6, 67.4, 67.0, 68.5%; p 0.999 |
| a1 + x6·x8 (unlabelled) | 10.6% (−843.6); 13.2, 9.1, 9.5, 10.0%; p 0.197 | 10.7% (−843.2); p 0.537 | 10.9% (−848.8); p 0.833 |
| B | 1.3% (−845.7) | 0.2% (−847.4) | 0.3% (−852.3) |
| a2 | 0.2% (−847.6) | 0.3% (−846.9) | 0.2% (−852.6) |
| A | under 0.05% (−850.9), 7,261 draws visiting | under 0.05% (−851.8), 4,604 | under 0.05% (−857.1), 10,172 |
| G, S, T, a5 | under 0.05% each (G 1 draw, S none, T 26, a5 565) | under 0.05% (S, G none) | under 0.05% |
| x1·x2 terms 3, 4, 5 | 3.8, 96.0, 0.2% (each chain 94.8–96.7% at 4) | 1.4, 98.4, 0.2% (96.9–98.6%) | 2.2, 97.6, 0.2% (96.7–98.1%) |

The rest of the weight sits on a1 with one term added or doubled (x2·x9·x10, a third x3, x8, a second x5, x9·x10, each 0.7–1.8%). Every chain of every arm visits a1 and a2.

*Against what would count:*
- (1) every production chain tunnels in all three arms, fewest 9, 4 and 13 round trips; none reached the 20 sized for.
- (2) the chains' SD of ln Ẑ_c is 0.34 and 0.14 nats in arm (i) at the default gain start and arm (ii), and 1.09 in arm (i) from 0, just past the threshold; the one labelled basin at 5% or more, a1, passes the χ² in every arm (p 0.92–1.00).
- (3) ln Ẑ(1) is −841.36, −840.96 and −846.60, all above −864.6.

*Added after the run, to locate the arms' difference* (`acf9d2d`, `Friedman_ArmOneAgainstSmcAndThePrior_BelowTheLowerSplit`; written once the difference was seen, so not pre-registered). The arms agree on ln Z(1) − ln Z(0.10) (arm i's curve −715.5, arm ii −715.25 ± 0.07) and part on ln Z(0.10). Three estimates of ln Z(b) on the linear form, below it:

| b | importance sampling, 200,000 exact prior draws | SMC from the prior, a stage forced at each b (4 × 4,000) | arm (i)'s curve, the same run | arm (i) per chain | arm (i)'s Kish ESS at b |
| --- | --- | --- | --- | --- | --- |
| 0.01 | −15.13 ± 0.00 | −15.13 ± 0.00 | −15.13 ± 0.00 | −15.13 each | 9,648 |
| 0.02 | −30.17 ± 0.00 | −30.18 ± 0.00 | −30.18 ± 0.00 | −30.17 to −30.18 | 13,500 |
| 0.03 | −45.00 ± 0.01 | −45.02 ± 0.01 | −45.01 ± 0.00 | −45.00 to −45.02 | 26,593 |
| 0.05 | −73.02 ± 0.21 | −73.10 ± 0.07 | −73.02 ± 0.01 | −72.98 to −73.04 | 165,678 |
| 0.07 | −98.84 ± 0.50 | −97.87 ± 0.16 | −97.61 ± 0.01 | −97.58 to −97.64 | 136,703 |
| 0.10 | −136.52 ± 0.85 | −131.42 ± 0.31 | −125.93 ± 0.06 | −125.78 to −126.02 | 49,007 |

- The prior's ℓ₁ (the 200,000 draws): min −1560.1, 1% −1532.1, median −1521.5, 99% −1414.7, max −1245.6; the intercept alone −1518.4, inside the window. The denominator's Kish ESS is 8,163 of 500,000; production put 1,513 draws below the window, and its lowest bin's occupancy is 0.73–0.76 of flat on every chain.
- SMC with the forced stages reads ln Ẑ(0.10) −131.42 ± 0.31 (66, 64, 66, 67 stages; 172–220 families), against the pre-registered run's −131.35 ± 0.21. The rerun of arm (i) reproduced its ln Ẑ(1), −841.36 (asserted), in 898 s against the first run's 779.

*What the tables show* (kept apart from them):
- *The flat histogram crosses Friedman's band.* Every production chain of every arm made round trips across a window of 447–773 nats that the pilots' ladders, at 16–21 rungs, never crossed (0 round trips per traced scan). The trips were fewer than sized: 4–17 a chain against the 20 planned, and the arms took 13–15 minutes rather than 30, the sweeps running about twice as fast as the sizing measured.
- *The chains agree with one another.* Arm (i) at the default gain start and arm (ii) put their chains' ln Ẑ_c within 0.34 and 0.14 nats (SD); arm (i) from 0 within 1.09, its chain 2 the one apart (4 round trips, b = 1 ESS 1,656). Every chain of every arm gives a1 61–71% of its weight and four x1·x2 terms 94.8–98.6%.
- *The posterior they agree on is a1's family*, 4 × x1·x2 + 2 × x3 + x4 + x5 at 66–68% with a1 + x6·x8 at 11%. Only R2's 16,000-particle replicate 0 found a1 (68.8% there, with a1 + x6·x8 at 6.7%); no MCMC chain, no other SMC replicate and no reference arm held it. A, B, G, S and T, the basins M4's chains and the reference path held, carry 1.3% or less together; A, M4's best chain's basin and 65–98% of every reference arm, carries under 0.05%, though 4,600–10,200 production draws visit it.
- *Their evidence sits far above every earlier reading.* ln Ẑ(1) −841.4 ± 0.2 and −841.0 ± 0.5 from b_lo = 0 lie 23 nats above the 99.9% bound and 45–46 above the reference arms. a1's log mass, −841.4 to −841.8, against R2's −855.6 for it.
- *The two ln Z(1) readings part by 5.4–5.6 nats, all of it below b = 0.10.* Below b = 0.05 the flat histogram, SMC and exact prior importance sampling agree to a tenth of a nat, so the flat histogram's prior-mass denominator reads the prior correctly. Above it they part in the order a sampler missing rare high-likelihood mass would show: importance sampling falls first, SMC next. Each estimator is unbiased in Z and reads low in ln Z where it misses mass (R2's Markov argument). So SMC from the prior, on the linear form, under-reads ln Z(0.10) by about 5.5 nats, and arm (ii)'s −846.6 through it is low by as much. Arm (ii)'s ratio agrees with arm (i)'s. This is a reading of three estimators' order, not a proof that the flat histogram's is right. The linear form's own structure below 0.10 was not mapped (the band measurement's splits are on S10's ladders, the integrated form).
- *σ* covers 1 in every arm, as R2's 4,000- and 16,000-particle runs did and the reference arms did not.
- *The gain's start.* From the default start (t = 42,000, the pilot's 168,000 draws over four chains) the learning left 29% of its visits above the window and a histogram 2.80 from flat, and its weights still moved 3.4 nats (RMS) between 10% and 50% of the learning. From 0 the learning's histogram was 0.15 from flat, with 8% above the window. Production from the default start made more round trips (47 against 27) and its chains agree more closely. One seed each, so neither start is read as the better.
- *The window's top* came from a pilot that never found a1, so ℓ_hi (−778.3, −776.0) lies below the posterior's ℓ₁. Production put 1–13% of each chain's draws above it, where the slope-1 extension holds the posterior's own law, and the b = 1 weights read them exactly. The learning spent 8–29% of its visits there, which SAMC's update leaves inert.
- *D10's verdict on π_w's functionals* reads converged at 1.01 for arm (i) at the default start and arm (ii); from 0 two functionals sit at R̂ 1.01.
- *Strength.* One master seed per arm (arm i's two runs share their pilot and seeds). Each chain's trips are 4–17, so the agreement rests on few independent crossings. The chains' spread is the evidence's error, 3 degrees of freedom. The bound rests on R2's twelve replicates.

*Against the recorded readings, at stated cost* (each with its seed's single run):

| sampler | cost | ln Ẑ(1) | replicates' or chains' SD | where its weight sits |
| --- | --- | --- | --- | --- |
| R2, SMC from the prior, 1,000 / 4,000 / 16,000 particles × 4 | 29–31 s / 1 min 47 s – 2 min 7 s / 9–11 min | −875.8 ± 4.2 / −874.4 ± 3.5 / −861.8 ± 3.3 | 8.4 / 7.1 / 6.5 | A, B, a4, a5–a8 by replicate; a1 in one 16,000 replicate (−855.2) |
| reference path, 8 arms at 1,000 and 4,000 particles | 20–41 s / 161–296 s | −886.4 to −887.5 | 0.26–1.45 | A 65–98% |
| M4, 4 untempered chains of 7,000 sweeps | seconds | none | R̂ of M 2.35 | G, A, S, B: one basin set per chain |
| flat histogram, arm (i), default start | 13 min (45M sweeps a chain) | −841.36 ± 0.17 | 0.34 | a1 66% |
| flat histogram, arm (i), the gain from 0 | 14 min | −840.96 ± 0.54 | 1.09 | a1 68% |
| flat histogram, arm (ii), b_lo = 0.10 | 15 min + SMC to 0.10 (42 s) | −846.60 ± 0.22 through SMC; ratio −715.25 ± 0.07 | 0.14 | a1 68% |

- *a1 and a2:* every chain of every arm visits both; a1 holds 66–68% and a2 0.2–0.3%.

**The times** (Release, built into `artifacts\chip-fh2\`; other sessions may have loaded the machine):
- The flat-histogram suite, `BassFlatHistogramTests`, 23 tests: about 1.5 minutes; fixture A's five fits 5–8 s each, B and C 5–6 s.
- The sizing on seed 8101: 33 s and 15 s at 20,000 iterations, 73 s and 52 s at 100,000 × 10 sweeps.
- The qualification: 43 minutes for its four tests (SMC to 0.10 42 s; the arms 13 min 3 s, 14 min 25 s, 15 min 1 s); the reading below the split 17 minutes (prior draws, SMC with marks 119 s, arm (i) again 898 s).

**The suites**, at `78708b3` with Release builds into `artifacts\chip-fh2\`:
- *Regression:* 1,331 passed, 0 failed, in 4 min 41 s. The new: `BassFlatHistogramTests` 22, `BassFlatHistogramSweep` 1 (it runs only under `SEED_SWEEP`) and `BassFlatHistogramQualificationTests` 5 (they run only under `BASS_QUALIFICATION`). The first full run failed one, `BassMovePassTests.Friedman_AtN500_IsS10s`, whose data hash reads the default config's JSON; it is re-pinned in `78708b3` with its history in the comment, the data unchanged, as N14, stage 8 and the coefficient prior's two axes re-pinned it.
- *Samplers:* 262 passed, 0 failed, in 33 s, as at Phase 1's last report.
- *The BASS oracle classes* (`Maths.Oracle.Tests.Bass*`, 13 classes): 90 passed, 0 failed, in 4 min 42 s.
- `BassBitIdentityTests`' seven hashes hold.

**The review passes and their fixes.**
- *Step 1.*
  - The held-τ chain read P(M = 0) at z −4.47, inside its slack; a probe at 120,000 sweeps on three seeds read 0.0029–0.0034 against 0.0034, so the batch SE under-reads a rare state, and the check was left as written. The probe was removed.
  - A doubled blank line in the test file, fixed in step 2's commit.
- *Step 2.*
  - The resume check compared tuples of arrays, which compare by reference: rewritten field by field.
  - `Effects` over its default grid exceeded its value budget at 4,000 draws: the check runs at 10 grid points.
  - The residual check is a functional fit's: the scalar fit's null default is asserted instead of a call that throws.
  - *Coverage added:* each draw's weight checked against e^{ℓ₁}/w(ℓ₁) of its own state, ℓ₁ rebuilt from its structure and τ, so the alignment of `DrawWeights` with the stores is shown independently of `Reweighting.Weigh`; with τ held, every draw's τ asserted the held one rather than an E[τ] row whose error is 0.
- *Step 3.*
  - An overlong guard reflowed.
  - The exploration found SAMC's fixed point refusing a sparse pilot's ĝ, some with no node empty (seeds 7801–7803, K = 20–60); the adapter's message now names the window, K, b_lo and the empty nodes, and the gap check covers both outcomes.
  - The warm-against-flat comparison's length (8,000) was chosen after its exploratory table, said in the test and here.
- *Step 4.* The harness was smoke-run at tiny lengths on seeds 8101–8105 before its constants were committed; the run found nothing to fix in it.
- *Final, against the brief's list:*
  - *Correctness:* the moves' acceptance on π_w with D01's ratios and the inner prior (enumeration, step 1); τ's slice on p(u)·w(ℓ₁) (sampled-τ enumeration, E[u]); `Reweigh` after every installed weight set (the engine reweighs every chain after each update and on restore; `Reweigh` itself checked); the weights frozen in production (the readout's weights and the draws' are one object); `DrawWeights` aligned with the stores (rebuilt per draw); the evidence only from b_lo = 0 (b_lo = 0.3 has none); nothing existing changed in its draws (`BassBitIdentityTests`' hashes).
  - *Coverage:* τ held and sampled; both coefficient forms; Jeffreys σ² and the flat intercept from b_lo = 0 with no evidence; an improper τ prior refused at b_lo = 0 and accepted above; a functional response; a chain that never tunnels, flagged; resume mid-learning and across into production; a pilot whose ĝ has a gap, interpolated and refused.

**For the main session.**
- *Three new estimates of ln p(y) on Friedman, all above the bound:* −841.4 ± 0.2 and −841.0 ± 0.5 from b_lo = 0, and −846.6 ± 0.2 through SMC's ln Ẑ(0.10), which the reading below the split puts about 5.5 nats low. Each is one seed's run. Before the posterior of a1 is taken as Friedman's, the arms want a second master seed and the 20 round trips a chain they fell short of; at the measured 52,000–63,000 sweeps a second a chain, 30 minutes buys about twice these lengths.
- *SMC from the prior on the linear form loses mass between b = 0.05 and 0.10*, where its lineages fall to 172–220 of 4,000, below the band measurement's splits (which were read on the integrated form). D05's path for SMC is the integrated form; this reading is the linear form's only.
- *The window's top from a pilot that never found the posterior's basin* lies below it; the reweighting stays exact through the slope-1 extension, but the learning's visits there are inert. A window read from the run's own production, or a margin, is a later option, not built.
- *The fixed point's refusals.* SAMC's fixed point refuses a ĝ it cannot flatten, as where the pilot left a segment empty beside a drawn node, not only where a node is empty. Friedman's pilots were accepted (no node empty, a few under 10 draws), but a sparser pilot could be refused there. This is the engine's documented rule, not changed here.
- *`BassEvidence.SteppingStone` is nullable now*, so nine reads of `evidence.SteppingStone.…` on ladder and SMC fits warn (CS8602) in files this chip does not own: `tests/oracle/BassSmcParityTests.cs:142`, `BassTemperingParityTests.cs:276`, `tests/regression/BassFriedmanBandTests.cs:515`, `BassSmcQualificationTests.cs:124`, `BassSmcTests.cs:336, 353`, `BassTemperingSizingTests.cs:75, 319`, `BassTemperingTests.cs:644`. A `!` each clears them.
- *Not logged:* `BassResolved` (not this chip's file) logs SMC's unread settings; the flat histogram's (`Chains`, `BurnIn`, `Samples`, `BatchSize`, `MinEssPerChain`, `StoppingScope`) are documented on `BassSampling.FlatHistogram` but not logged.
- *`BassFlatHistogramSweep.cs`* holds the warm-against-flat claim for `SEED_SWEEP`, beside Phase 1's `FlatHistogramSweep.cs`; both could fold into `ComparativeSweep`.
- `docs/Maths` and the checklists are untouched. `artifacts\chip-fh2\` is deleted.
- No default is recommended.
