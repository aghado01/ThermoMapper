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
- **The ordering:** a tour over the S emulator draws, nearest neighbour then 2-opt, on d(s, s′), the spread of the log-likelihood difference ℓ(θ; e_s) − ℓ(θ; e_s′) over the first leg's final particles, weighted.
  - This is the quantity a bridge's conditional ESS falls with, so the tour minimizes what the bridges pay.
  - It costs one S × N evaluation of those particles, kept as an S × N matrix, then S² × N arithmetic.
  - Any order gives a valid estimator, since each μ_s is targeted exactly, so letting the ordering depend on the population changes efficiency only.
  - *The fallback,* recorded: the L₂ distance between the draws' predictions at a fixed space-filling design of 64 points in the unit θ box. It is cheaper, but it measures closeness where the posterior has no mass.
  - The random order is the comparison.
- **The bridges:** between e_{s−1} and e_s, the geometric bridge (1 − t)·log p̃_{s−1} + t·log p̃_s. t is placed by conditional ESS, at the same target as the first leg.
  - The incremental weight at a step from t to t′ is (t′ − t)·(log p̃_s − log p̃_{s−1}), evaluated at the particle.
  - The unnormalized densities' constants cancel in the normalized weights.
  - *The bridge as a target on the likelihood path:* log π(θ) + ℓ_{s−1} + t·Δ, with Δ = ℓ_s − ℓ_{s−1}. Its base is the previous conditional, and Δ plays the likelihood, so `SmcSampler.TemperToPosterior` walks t to 1 on it (argon's note, 2026-10-29).
  - *Its first step* moves to a new target, and `SmcRun.ConditionalEss` re-evaluates the whole population on every call for a new target, so a bisection through it would cost about 60N evaluations. Instead:
    - evaluate e_s on the particles once;
    - bisect t₁ with `SmcWeights.ConditionalEss` on the cached t·Δ;
    - `Advance` to the bridge at t₁;
    - then `TemperToPosterior`, which reuses the evaluations.
  - *One `SmcRun` across the first leg and every bridge,* so `Families` tracks the initial lineages over the whole tour. Its `SmcOptions.MaxStages`, whose stage seeds are derived at `Start`, scales with S: a per-bridge cap times S, plus the first leg's.
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
  - `src/maths/samplers/smc/` stays as it is: its public API (`SmcRun.Start`, `Stage`, `ConditionalEss`, `Advance`, `SmcSampler.RunTempered`, `TemperToPosterior`, `NextInverseTemperature`, `SmcWeights`) is enough. Name any further need in the report rather than edit it;
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
- C3's inner runs under the default gate (R̂ 1.01, ESS 100 per chain, `MaxSamples` 4,000): how many converge. This is stage 7's open item;
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

### C5 built, Steps 1–4 — 2026-10-29

From the C5 chip. All four steps are built and committed on `main`, and nothing is pushed. `src/maths/samplers/smc/` is untouched. With `Smc` null, stage 7 is unchanged bit for bit: a fingerprint of its nested and tempered samplers' draws (three fixtures) was identical before step 1 and after steps 1 and 3.

**Commits.**
- `86c5969`: Step 1, the bridged target and its block-scan kernel.
- `2e77db3`: Step 2, the ordering.
- `eb05485`: Step 3, the cut SMC.
- `cb0ecef`: Step 4, the comparison (exploration, behind `CUT_SMC_COMPARISON=1`).
- This report, in its own commit after these.

**API** (`src/maths/regression/spline/bass/calibration/`, namespace `Maths.Regression.Spline.Bass.Calibration`).
- *The bridge* (`CalibrationBridge.cs`): `CalibrationBridge(from, to)` is an `ISamplingTarget<CalibrationState>` on D05's likelihood path.
  - Its base is log π + ℓ_a and its likelihood Δ = ℓ_b − ℓ_a, so the stage at t is log π + ℓ_a + t·Δ.
  - `BridgeEvaluation` carries both draws' `CalibrationEvaluation`s and `Difference` (Δ).
- *The kernel* (`CutSmcKernel.cs`): `CutSmcKernel(model)` is the `ISmcKernel<CalibrationState>`, one block-scan sweep per move, aimed at a stage: one draw at rung b (the first leg) or a bridge at t.
  - `StageRates` gives each stage's `CutSmcBlockRates`: the walk's, the variance steps' and the single-site level steps' attempts and acceptances.
  - `CutSmcParticles(model)` draws exactly from the prior and forks.
- *The ordering* (`CutSmcOrdering.cs`): `LogLikelihoods` (the S × N matrix, on fresh copies), `Distances`, `NearestNeighbour`, `TwoOpt`, `Tour`, `Length` and `Random`, with `CutSmcOrder { Tour, Random }`.
- *The run* (`CutSmc.cs`): `CutSmcSampler.Run(model, maxDegreeOfParallelism = −1)` returns a `CutSmcCalibration`:
  - `OuterDraws`, `Order` (the visits, as positions), `Distances` (a tour's, else null) and `TourLength`;
  - `Populations`, one `CutSmcPopulation` per visit: the draw, the visit, `CalibrationDraw`s in physical coordinates, normalized weights, `LogEvidence` ± `LogEvidenceError` (log p(y | e_s) as the run estimates it, never read by the cut), and its `CutSmcLeg`;
  - `CutSmcLeg`: from (null for the first leg), to, the tour distance, the first step t₁, the leg's stage reports and its block rates, with `StageCount`, `MinEss`, `MinConditionalEss`, `Families` and `Acceptance`;
  - `Stages`, `TotalStages`, `Families` and `Pooled()`.
  - All diagnostics are SMC's own (ESS, conditional ESS, families, acceptance); none is a chain R̂.
- *Changes inside `calibration/` to stage 7's files:*
  - `CalibrationState` caches the prediction and the covariance factor for the last two targets that read it (was one), so a bridge reads both draws from the cache. The internal `Fork()` makes a distinct state sharing the values and caches, so an offspring's later cache writes stay its own.
  - `CalibrationModel.DrawPrior(rng)` is an exact prior draw, v ~ N(0, s_v²·I). `Start` keeps v at 0 for chains, as before.
  - `DiscrepancyBlock.Conditional` and `CategoricalBlock.Enumerate` take weighted likelihood terms (`LikelihoodTerm`): precision I/s_v² + Σ_k w_k·DᵀΣ_k⁻¹D, and Σ_k w_k·ℓ_k over the level combinations. Stage 7's single-term path is arithmetic-identical.
  - `CalibrationConfig.Smc` was added. `BassCalibrationRecoveryTests`' fits and helpers went from private to internal, so the new tests reuse its BASS fits.

**Configuration.** `CalibrationConfig.Smc` is a `CutSmcSpec?`, null by default, which leaves the cut posterior to C3's `NestedCalibrationSampler`. `CutSmcSampler.Run` refuses a null spec and points there. The defaults are provisional, read in the checks and not tuned:
- `Particles` 1,000;
- `TargetConditionalEss` 0.99, for the first leg's stages and every bridge's;
- `MutationSteps` 5 sweeps;
- `ResampleBelow` 0.5;
- `Resampling` multinomial;
- `OuterDraws` null, with C3's meaning;
- `Order` Tour;
- `MaxStagesPerLeg` 1,000, so the run's `MaxStages` is that times S + 1.

The run reads `Sampling.MasterSeed`. Its outer draws come from C3's stream for the same seed, so the two samplers share them; its own streams sit beside C3's inner ones. Each field is refused by name when out of range, and the config round-trips through JSON. A config written before stage 8 reads with `Smc` null.

**Checks, Step 1** (`BassCalibrationBridgeTests`, 16 tests).
- *The bridged density* against (1 − t)·log p̃_a + t·log p̃_b at t = 0, 0.3, 0.7 and 1: worst relative gap 2.4·10⁻¹⁶ (tolerance 10⁻¹²). A step from t to t′ reweights by (t′ − t)·Δ, and the slope is Δ. Off the support the bridge is −∞.
- *The caches,* counted by emulator calls: a bridge evaluates two predictions once. A change to v keeps both draws' predictions and factors, a change to σ² keeps both predictions, and a change to θ keeps both factors. A fork carries the caches, and its later writes stay its own.
- *The prior draw:* θ, the levels, v ~ N(0, 0.09) and the IG group's E[ln σ²], each within 5 SE over 40,000 draws. An improper group is refused.
- *The scan's law,* by quadrature (θ × ln σ² × level, v integrated out), within 5 MCSE plus the grid:
  - at the bridge's ends, against the single-draw conditional: t = 0 gives E[θ] 0.5153 against 0.5164 (MCSE 0.0028), and t = 1 gives 0.5549 against 0.5506 (MCSE 0.0048);
  - on the first leg, one draw at b = 0.5: 0.5727 against 0.5720 (MCSE 0.0049);
  - inside the bridge, against its own quadrature: t = 0.3 gives 0.5297 against 0.5304 (MCSE 0.0031), and t = 0.7 gives 0.5374 against 0.5433 (MCSE 0.0046). E[ln σ²] and the level frequencies agree too;
  - the bridged quadrature reproduces the sampler's at t = 1 and at a tempered end to 10⁻⁹.
- *Each block alone* at t = 0.3 and 0.7, everything else held:
  - the walk against quadrature over θ: 0.5150 against 0.5153 (MCSE 0.0006) and 0.5272 against 0.5278;
  - the variance step against quadrature over ln σ² with its Jacobian: −5.0135 against −5.0183 (MCSE 0.0048) and −4.9955 against −5.0003;
  - the levels' bridged enumeration against brute force with the dense Σ's to 10⁻¹², and its draws within 5 binomial SE;
  - v's bridged Gaussian against its closed form: the mean to 10⁻¹⁰, the covariance to 10⁻¹¹, and the draws' means within 5 SE;
  - single-site level steps above a cap of 8 against the bridged enumeration, within 5 MCSE plus 0.002 (acceptance 0.08 at t = 0.3, 0.05 at 0.7);
  - under the multivariate formulation, and functional without a discrepancy on σ²_e·I: the bridge's density to 10⁻¹², and the walk and the variance step at t = 0.5 against quadrature (for example, multivariate: 0.5025 against 0.5020, MCSE 0.0005).
- *Refusals:* the powered path, another model's target, and a target that is neither a draw nor a bridge. A bridge's ends, and a bridge from a draw to itself, reduce to one term.

**Checks, Step 2** (`BassCalibrationOrderingTests`, 5 tests).
- Every draw is visited once from the first leg's, at 1 to 120 nodes.
- 2-opt never lengthens a path (200 instances, from nearest neighbour and from random paths). Against brute force over every open path at seven nodes, the tour was optimal in 29 of 30 instances (mean ratio 1.0001, worst 1.0019). On a line from an end it is the sorted order.
- The distances match the weighted spread of the differences computed directly, to 10⁻⁹ relative, under rows offset by 10⁶. A constant added to one draw's log-likelihood (its normalizer) moves nothing.
- The tour and the random order are deterministic in their inputs and seed, and the random order is uniform (within 0.015).
- With weighted draws, the tour orders C3's resampled set, each duplicate beside its twin, by the triangle inequality on d. The S × N matrix leaves the population's caches untouched.

**Checks, Step 3** (`BassCalibrationCutSmcTests`, 14 tests). Errors are read from paired replicates: the same seed for both samplers, so the same outer draws. For each mean and variance, the paired differences must average within 6 SE of 0, the SE read from their own spread over 6 replicates. C3 runs at a fixed inner length (an unreachable ESS gate, every run pooled).
- *Against quadrature,* stage 7's pooled-posterior fixture, N = 500, 6 replicates (6 SE plus the grid):
  - pooled E[θ] 0.5359 against the exact mixture's 0.5340 (SE 0.0013, 1.5 SE);
  - E[ln σ²] −5.1531 against −5.1447 (−1.0 SE);
  - the level probabilities within 0.9 SE;
  - each μ_s's E[θ] against its own conditional: −0.7, 1.8 and 1.1 SE;
  - each draw's log p(y | e_s), the run's normalizer after its leg, against quadrature with every prior's normalization: 2.5565 against 2.5809, 2.5957 against 2.6249 and 2.6438 against 2.6538, each within 1.6 SE. This checks the first leg's increments and every bridge's.
- *Against C3,* pooled means and variances of θ, the level indicators, σ² and v:
  - a discrepancy, two groups (IG and half-Cauchy) and a level (stage 7's conditional fixture), N = 400 against C3 with a 500-scan burn-in and 1,000 draws a chain: every moment within 2.0 SE. Each μ_s's θ against C3's inner run on the same e_s: within 1.7 SE;
  - BASS's basis fit and direct fit (stage 7's recovery fixture, six outer draws of 200), the default N against C3 with a 2,000-scan burn-in and 4,000 draws a chain: within 2.0 SE (basis) and 1.5 SE (direct);
  - the multivariate formulation, the default N against C3 with a 500-scan burn-in and 1,500 draws a chain: within 1.0 SE;
  - two categoricals above a cap of 8 (single-site steps in both samplers), N = 300 against the same C3: within 1.8 SE.
- *The account* (on every run): weights summing to 1, every visit's mass exactly 1/S, populations on their draws, the order a permutation from the first leg's draw, families in [1, N], and the stages adding up.
- *The cut kept:* a draw that misfits by ±0.15 out of the discrepancy's span has log p(y | e) −2.28 against 2.64 and 2.67. A joint posterior would weigh it at 0.36%; the cut pools it at exactly 1/3.
- *S = 1:* the first leg alone, resampled at every stage (every offspring after the first a fork), 34 families left: E[θ] within 6 SE, the SE bounded by the families' count.
- *A bridge whose conditional ESS collapses:* draw 1's conditional moved about five posterior sds.
  - Its bridge took 45 stages (the near bridge 9), first t 0.028, every stage at the 0.99·N target, resampled three times, 196 families left.
  - Its population matches its conditional's quadrature: 0.4923 against 0.4994, SE ≤ 0.024.
  - With `MaxStagesPerLeg` 8 the run stops with a message naming the leg.
- *Weighted draws* (0.2, 0.5, 0.3), eight outer draws resampled by weight, 6 replicates: the outer draws are C3's, each draw's copies are visited in one run, and each draw's pooled mass is its share of the visits. Pooled E[θ] 0.5351 against the weighted mixture's 0.5359 (SE 0.0015, the outer resampling included).
- *Determinism:* identical draws, weights, order and evidence at 1 and 6 threads; another seed differs. A random order keeps the same first leg, takes no distances, and is the same at 2 and 6 threads.
- *The config:* the JSON round trip with every field off its default; each field's refusal by name; an unknown order refused.
- *Refusals:* a null spec (pointing to C3), an improper prior, a parallelism of 0.
- *A bridge's first step:* the largest t₁ whose conditional ESS stays at the target, as `NextInverseTemperature` bisects. It is 1 between equal conditionals.

**The comparison** (Step 4, `BassCalibrationCutSmcComparisonTests`, behind `CUT_SMC_COMPARISON=1`). This is exploration: it reads and recommends nothing. It runs on stage 7's recovery fixture, field seed 1 unless stated. The machine was shared with the C4 chip, so wall and processor times moved between runs; the SMC's basis run read 26 s and then 11.5 s wall. The ratios within one run are the readings.

*The tour against a random order.* Twenty outer draws of the basis fit's 200, at the default spec, over five seeds. Each seed's two runs share the outer draws and the first leg, and both paths are measured on the tour's distances.

| Seed | Order | Path length | Bridge stages | Least ESS (bridges) | Families left | Wall s |
|---|---|---|---|---|---|---|
| 1 | Tour | 3.83 | 52 | 496 | 333 | 37.7 |
| 1 | Random | 9.12 | 113 | 453 | 214 | 76.3 |
| 2 | Tour | 4.12 | 54 | 497 | 345 | 48.4 |
| 2 | Random | 10.72 | 128 | 496 | 213 | 91.0 |
| 3 | Tour | 3.50 | 54 | 495 | 356 | 49.7 |
| 3 | Random | 8.38 | 116 | 490 | 270 | 74.4 |
| 4 | Tour | 4.66 | 53 | 494 | 339 | 49.7 |
| 4 | Random | 13.23 | 125 | 493 | 262 | 68.7 |
| 5 | Tour | 5.14 | 54 | 500 | 331 | 39.5 |
| 5 | Random | 11.46 | 108 | 487 | 244 | 56.9 |
| Mean | Tour | 4.25 | 53.4 | 497 | 341 | 45.0 |
| Mean | Random | 10.58 | 118.0 | 484 | 241 | 73.5 |

Over the bridges, the median was 12.0 stages per unit of distance.

*Cut SMC against C3 at matched error on E[θ].* Six fixed draws of the 200, each taken once, so the replicates' spread is the inner error alone. Five replicates each; the SMC at its defaults, C3 under the default gate with every inner run pooled for the error. "C3 to SMC's error" scales C3's time by (its per-run sd over the SMC's)². The table is the second run's:

| Fit | θ | SMC E[θ] ± per-run sd | C3 E[θ] ± per-run sd | SMC wall / CPU s | C3 to SMC's error, wall / CPU s |
|---|---|---|---|---|---|
| basis | θ₁ | 1.29137 ± 0.00035 | 1.29139 ± 0.00048 | 11.5 / 54 | 12.1 / 20 |
| basis | θ₂ | 0.39832 ± 0.00011 | 0.39835 ± 0.000039 | 11.5 / 54 | 0.9 / 1 |
| direct | θ₁ | 1.42238 ± 0.0018 | 1.42450 ± 0.0011 | 4.7 / 11 | 3.1 / 3 |
| direct | θ₂ | 0.40756 ± 0.00073 | 0.40711 ± 0.00069 | 4.7 / 11 | 8.0 / 7 |

The first run, with the C4 chip busier, read the SMC at 26.2 s wall and 97 s CPU on the basis fit, against C3's 23.8 s and 34 s to the SMC's θ₁ error. Five replicates estimate each per-run sd to about ±35%, so each ratio is good to a factor of about 2.

*C3's inner runs under the default gate:* R̂ 1.01, ESS 100 per chain, 4 chains, at most 4,000 draws per chain, after the default 2,000-scan burn-in. In the same runs:
- basis: 30 of 30 converged, with a median of 1,500 draws per chain (at most 2,500);
- direct: 30 of 30, median 1,500 (at most 2,750).

*θ*'s coverage over five seeds:* stage 7's four cases and seeds, six outer draws, the SMC at its defaults. The intervals come from the recovery test's own weighted quantile, not C4's `DrawSample`.

| Case | θ₁ | θ₂ | Stage 7's C3 | SMC stages |
|---|---|---|---|---|
| basis | 4/5 | 5/5 | 5/5, 5/5 | 85–94 |
| direct | 4/5 | 5/5 | 4/5, 5/5 | 62–71 |
| discrepancy | 5/5 | 5/5 | 5/5, 5/5 | 103–111 |
| two groups | 5/5 | 5/5 | 5/5, 5/5 | 91–105 |

The basis miss is seed 3: θ₁ 1.341 [1.304, 1.377], with θ* = 1.3 just below.

*What the readings show.*
- *The tour does what it is for.* Against a random order it halves the bridge stages (53 against 118 at S = 20) and keeps more lineages (341 against 241 families), and wall time falls by about 40%. The least ESS is about the same, since the adaptive schedule holds each step's conditional ESS whatever the order.
- *The distance ranks a bridge's cost but does not size it.* A bridge took about 12 stages per unit of d, where the small-step Gaussian rule, 1/√(2(1 − 0.99)), gives 7.1. d is read at μ_1's particles, and Δ spreads further at later conditionals and over larger steps.
- *At S = 6 on these fits, C3 under the default gate reaches the SMC's error on E[θ] with less processor time:* about a third on the basis fit's θ₁, far less on its θ₂, and comparable on the direct fit.
  - The SMC's cost is mostly its first leg, about 70 of the basis fit's 85–94 stages.
  - Each further draw cost about 2.8 bridge stages at S = 20, where C3 pays a whole inner run. Larger S is untested here.
- *Every C3 inner run converged under the default gate,* on both fits, within at most 2,750 draws per chain. Stage 7's open item reads converged at the default budget.
- *On the direct fit, C3 under the default gate reads E[θ₁] 0.0021 above the SMC,* about 4 SE over five replicates, as C3 at 1,000 fixed draws did in Step 3's review. C3 at 8,000 draws matched the SMC to 0.1 SE, so this is a short-run effect of C3 that its R̂ and ESS gate does not see.
- *Coverage matches stage 7's C3 reading,* except one basis seed where θ* sits just below the interval.

**Suites.**
- Regression (`Maths.Regression.Tests`, whole): 1,169 passed, 0 failed (14 min 23 s, on a machine shared with the C4 chip's runs).
- BASS and calibration oracle classes (`Maths.Oracle.Tests`, filter `Bass`): 71 passed, 0 failed (8 min 26 s).
- The comparison's three readings ran with `CUT_SMC_COMPARISON=1` (23 min); they assert only that they run. The matched-error reading was then rerun with more digits (2 min 37 s).
- No R oracle was touched, so R's stage-1 checks were not rerun.

**Review passes and their fixes.**
- *Step 1:*
  - added the bridge under the multivariate formulation and on σ²_e·I without a discrepancy, which no check covered;
  - the first leg's rung chain ran twice inside the bridge-end theory, and is now its own test;
  - corrected my cache test's count: a third draw displaces one slot, which costs two predictions on the bridge's next read, not one.
- *Step 2:* nothing to fix.
- *Step 3:*
  - *The cut check's premise was too strong.* The misfit draw's evidence fell 5 nats, not 10, as σ² widens to absorb it. The check now states what a joint posterior would do (0.36%) against the cut's 1/3.
  - *The basis fit's var θ₁ read −5.2 SE at N = 400.* A diagnostic (6 replicates; C3 with burn-in 600 and 3,000, and SMC at three settings) put it on the SMC's finite population:
    - at N = 400 and 5 sweeps, the pooled var θ₁ ran 3% under C3's (per draw 0.93–1.03 of C3's long runs);
    - at N = 1,000, or 20 sweeps, it matched within 1–3%;
    - C3's burn-in length moved it by 1.4%.
    The BASS-fit and multivariate checks now run at the default population.
  - *The direct fit's E[θ₁] read −3.9 SE against C3 at 1,000 draws a chain.* Rerun against C3 at 8,000 it read 0.1 SE, and C3 itself moved by 2.8 SE between the two lengths: C3's fixed-length bias, not the SMC's. The BASS-fit references now run 4,000 draws a chain after a 2,000-scan burn-in.
  - *Coverage gaps:* added the evidence check by quadrature, which tests every incremental weight end to end, and a random-order run.
- *Whole change,* the brief's list:
  - the cut is kept: no stage reweights an emulator draw, and each visit pools at 1/S whatever its evidence;
  - each kernel is invariant for its bridge, block by block and as a scan;
  - the incremental weights are (t′ − t)·Δ, and the evidences match quadrature;
  - the first leg starts from an exact draw of a proper prior, and an improper one is refused;
  - the pooling's weights are checked on every run;
  - forks are distinct states (S = 1 resamples at every stage);
  - every change is evaluated afresh: states are immutable, and caches are keyed by target.
  - Coverage: both formulations; with and without the discrepancy; levels below and above the cap; weighted draws; S = 1; a collapsing bridge.

**For the main session.**
- *Weighted emulator draws, a reading of the brief.* The decided text's "with weight 1/S. With weighted emulator draws, use their weights instead." sits beside Step 2's "with weighted draws, it orders the resampled set". I read them together: the outer draws are C3's, resampled by weight, and each visit pools at 1/S. A draw resampled k times thus carries k/S, its weight's estimate. Pooling the resampled set by the weights again would count them twice. If the intent was every retained draw once, pooled by its weight with no resampling, that is a second path through `OuterDraws`, and I did not build it.
- *`MaxStages`:* one field, `MaxStagesPerLeg`, serves as both the per-bridge cap and the first leg's, so `MaxStages` = cap·(S + 1). It caps the run's total, since `TemperToPosterior` runs a leg to its end; a run that exhausts it stops with a message naming the leg.
- *The kernel's scales:* the brief named the walk's per-stage scale. The variance steps get theirs the same way, 2.38 times ln σ²_g's weighted spread, as adaptation from the past.
- *R's decorrelation step is left out of the SMC sweep.* Its cadence (every 100 scans) has no counterpart in five sweeps per stage, and the population already spans θ. Every remaining block keeps its stage invariant.
- *Resampling is rare at the defaults.* At a conditional ESS of 0.99, the ESS falls by about 1% a stage. The fixtures' first legs (35–70 stages) and short bridges often end above 0.5·N without resampling, so the families stay at N. Mutation still moves every particle. The collapsing bridge, the basis fit (families 132–206 of 400) and S = 1 do resample.
- *The covariance factor dominates the sweep.* At m = 8 with seven low-rank columns, `FunctionalCovariance.Factor` (Jacobi SVD through `DiagonalPlusLowRank`) costs 164 µs against 1.6 µs for a prediction. A bridge sweep builds about four factors (two draws, a variance step per group), about 0.66 of its 0.92 ms. This is stage 7's finding again; C-O2's Cholesky route through I + UᵀD⁻¹U would cut it. I left the linalg and the covariance's factor route alone, since changing them would move stage 7's numerics. Parallel scaling: 5.8× at 6 threads, 7.2× at 18.
- *C3 at short runs is biased on the direct fit, and its gate does not see it.* At 1,000 fixed draws a chain, E[θ₁] sat 2.8 SE from its own 8,000-draw value. Under the default gate, with every inner run converged, it read 0.0021 above the SMC, about 4 SE. This bears on reading stage 7's short reference runs and on D10's gate, not on stage 7's tests.
- *The SMC core:* its public API sufficed. Two notes:
  - `SmcRun` does not expose its `MaxStages`, so the driver carries it;
  - a bisection on cached increments (`FirstStep`, the same rule as `NextInverseTemperature`'s on a weights-and-increments pair) would sit naturally in `SmcWeights`, beside `ConditionalEss`.
- *`LogEvidenceError`* is the core's genealogy estimate. It reads 0 when that estimate goes negative, as at S = 1 resampled at every stage with 34 families.
- *Test runtime:* the three new classes take 6 min 24 s wall run together on 18 cores. The summed test time is about 4 min in `BassCalibrationBridgeTests` and 6.4 min in `BassCalibrationCutSmcTests`, the latter mostly the paired replicates. The comparison runs only with `CUT_SMC_COMPARISON=1`.
- *Not done:* `docs/Maths` and the checklists, which are the main session's.
