# Stage 8, C5: cut-calibration SMC (chip brief)

Date: 2026-10-29. Builds the calibration checklist's C5, cut-calibration SMC with ordering and induced-density bridges. D09 (in [port-calibration-brief.md](port-calibration-brief.md), decided 2026-10-29) fixes its shape, following Mathews, Gopalan, Gattiker, Smith and Francom, *Sequential Monte Carlo for Cut-Bayesian Posterior Computation* (arXiv 2406.07555v2; a local copy is in `D:\aghado01\codex-scientiae\supellex-private\gauntlet\bayesian\2406.07555v2`):
- carry particles through μ_s = p(θ | y, e_s) over independent emulator draws e_s, reweighting by the ratio of conditional densities, then resampling and mutating with the inner kernel;
- order the e_s so that successive conditionals stay close;
- bridge successive conditionals through their induced densities, p(θ | y, e_{s−1})^{1−t}·p(θ | y, e_s)^t, because emulator draws are discrete. Mathews et al.'s linear tempering interpolates a continuous cut parameter and does not apply directly.

Stage 7's nested reference sampler (C3) is the reference this SMC is checked against. The choices this brief makes beyond D09 are in "Decided for stage 8" below.
- *Not in scope:*
  - fitting SMC, the parallel C4 chip's ([its brief](stage8-fitting-smc-brief.md));
  - semi-modular calibration (the extensions brief's E14), which can later reuse these bridges;
  - `docs/Maths` and the checklists, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- D09 whole, with its *(2026-10-29)* amendments and its as-built notes;
- the calibration segment's "SMC integration";
- stage 7's report ([stage7-calibration-brief.md](stage7-calibration-brief.md#reports)): its API, the conditional sampler's blocks, and the nested reference's inner-run error;
- Mathews et al., on the ordering and its finite-sample bounds.

In the code:
- *The SMC core,* `src/maths/samplers/smc/`. A cut posterior "walks through different conditional targets" (`SmcStage`'s summary). See also `SmcRun`, `SmcWeights` and `NextInverseTemperature`.
- *Calibration,* `src/maths/regression/spline/bass/calibration/`:
  - `ICalibrationEmulator`;
  - `CalibrationModel` and `CalibrationTarget` (the likelihood path, ℓ₁);
  - `CalibrationBlocks` (the discrepancy's Gaussian draw, the variance step, the categorical block);
  - `CalibrationSampler` (the θ walk and its block scan);
  - `NestedCalibration` (C3's reference and its pooling);
  - `CalibrationConfig`.

## Decided for stage 8

These are the main session's recommendations, taken by the user on starting this chip.
- **The first population** is tempered from the prior to μ_1 on D05's likelihood path, with stages by conditional ESS, as fitting SMC does. Every default prior here is proper.
- **The ordering:** a tour over the S emulator draws, nearest neighbour then 2-opt, on the L₂ distance between the draws' predictions at a fixed space-filling design of 64 points in the unit θ box: a Sobol sequence if the repo has one, else a seeded Latin hypercube, with each categorical combination where there are few.
  - Record the alternative, the distance read at the first population's particles, which adapts to the posterior but depends on the path.
  - The random order is the comparison.
- **The bridges:** between e_{s−1} and e_s, the geometric bridge (1 − t)·log p̃_{s−1} + t·log p̃_s. t is placed by conditional ESS, at the same target as the first leg.
  - The incremental weight at a step from t to t′ is (t′ − t)·(log p̃_s − log p̃_{s−1}), evaluated at the particle.
  - The unnormalized densities' constants cancel in the normalized weights.
- **The mutation** at a bridge stage is the conditional sampler's block scan, aimed at the bridge. Each block must leave the bridge invariant:
  - the θ walk and the variance steps are Metropolis on the bridged density;
  - the categorical block enumerates the bridged conditional;
  - the discrepancy's draw uses the bridged Gaussian, whose precision and mean combine the two conditionals' (it stays Gaussian, since each is quadratic in v).
  - The θ walk's scale is fixed per stage from the population (`Prepare`), which is adaptation from the past (D10).
- **The cut posterior** pools each μ_s's population, read when t reaches 1 at e_s, with weight 1/S. With weighted emulator draws, use their weights instead. The populations are recorded per e_s.
- **`CalibrationConfig.Smc`:** a `CutSmcSpec?`, null meaning C3's nested sampler. Its values are provisional, read in the checks and not tuned:
  - `Particles` 1,000;
  - `TargetConditionalEss` 0.99;
  - `MutationSteps` 5 block-scan sweeps per stage;
  - `ResampleBelow` 0.5;
  - multinomial resampling;
  - `OuterDraws` as C3's.
- **The diagnostics:** per bridge, the stage count, ESS, conditional ESS, families and acceptance; across the tour, the total stages and the families left. Keep them apart from chain R̂ (the SMC integration section).

## Working rules (a shared checkout, with a parallel chip)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **The C4 chip works beside you,** on BASS's fitting files (`Bass.cs`, `BassConfig.cs`, `BassRun.cs`, `BassDiagnostics.cs`, `BassFit.cs`, `BassPointwise.cs`, `BassEvidence.cs`, `BassResolved.cs`) and on `DrawSample` in `src/maths/samplers/ensemble/Readouts.cs`.
  - Leave those files alone. If a build fails in a file you do not own, wait and retry; never edit it.
  - Weighted quantiles change under C4. Where your checks compare quantiles, compare means and variances, or quantiles at equal weights.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-c5\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-c5\` when you finish.
- **Your files:**
  - `src/maths/regression/spline/bass/calibration/`, existing and new;
  - `src/maths/samplers/smc/`, for accessors only; name each one;
  - new tests `tests/regression/BassCalibration*Tests.cs`;
  - this brief's Reports section.
- **Commit your own files only, by explicit path,** one commit per step, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's working-copy edits (`issues/regression/argon-bass-qa.md`, `bass-questions.md`, `issues/spc/opus-spin-glass-threads.md`) stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.
- **Off means unchanged:** with `Smc` null, stage 7's nested sampler and its tests are as they were.

## The steps

Each step lands with its checks and ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next step. Name the passes and their fixes in the report.

**Step 1. The bridged target and its kernel.** A target for the bridge between two emulator draws at t, and the block scan aimed at it.
- *Checks:*
  - the bridged density against (1 − t)·log p̃_a + t·log p̃_b, to rounding;
  - at t = 0 and 1, the scan reproduces the single-draw sampler's law;
  - each block's invariance for the bridge at t = 0.3 and 0.7, by quadrature on a small fixture, the discrepancy's bridged Gaussian against the closed form included;
  - the categorical block's bridged enumeration against brute force.

**Step 2. The ordering.**
- *Checks:*
  - the tour visits every draw once;
  - 2-opt never lengthens it;
  - it is deterministic in the seed;
  - with weighted draws, it orders the resampled set.

**Step 3. The cut SMC.** The first leg from the prior, then the tour's bridges, each run on `SmcRun` over its stages. The pooled cut posterior, with the populations kept per e_s.
- *Checks:*
  - against C3's nested reference on stage 7's fixtures (a basis fit and a direct fit, a discrepancy, two variance groups, a categorical parameter): the pooled posterior means and variances of θ, σ² and v within their Monte Carlo errors, those errors read from replicate runs;
  - against quadrature of the exact equal mixture, on stage 7's pooled-posterior fixture (E[θ] 0.5340);
  - each μ_s's population against C3's inner run on the same e_s;
  - weights summing to 1, families reported, and the cut kept (no stage reweights the emulator);
  - the config's JSON round trip, and its refusals.

**Step 4. The comparison.** This is exploration: it reads and recommends nothing. On stage 7's recovery fixture:
- the tour against a random order: total bridge stages, the least ESS, families left;
- cut SMC against C3 at matched error: the time each takes to reach the same Monte Carlo error on E[θ];
- θ*'s coverage over five seeds, as stage 7 read it.

Report what the readings show, kept apart from the numbers.

## Review before reporting

Take a fresh-eyes pass over the whole change and fix what it finds:
- *Correctness:*
  - the cut kept;
  - each kernel invariant for its bridge;
  - the incremental weights;
  - the first leg starts from a proper prior;
  - the pooling's weights;
  - forks independent;
  - every change evaluated afresh (stage 7's cache rule).
- *Coverage:*
  - both formulations;
  - with and without the discrepancy;
  - a categorical parameter below the cap and one above;
  - weighted emulator draws;
  - S = 1;
  - a bridge whose conditional ESS collapses.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the API and the configuration;
- each check with its tolerance;
- the comparison's tables, then what they show;
- the suites with their counts (Regression, the BASS and calibration oracle classes);
- the review passes and their fixes;
- anything the main session must know.

## Reports

Newest first.
