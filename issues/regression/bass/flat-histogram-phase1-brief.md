# The flat-histogram sampler, Phase 1: the engine (chip brief)

Date: 2026-10-29. Builds Phase 1 of N17, the flat-histogram sampler, whose nine decisions the user took on 2026-10-29: [flat-histogram-design.md](flat-histogram-design.md). Read the design first, whole. Its decisions are settled; this brief fixes the build.
- *What it builds,* in `Maths.Samplers.FlatHistogram`:
  - the multicanonical weights;
  - SAMC's learning, with weights shared across chains;
  - the target wrapper;
  - the run (learning in burn-in, then frozen production);
  - the reweighting and the ln Z(b) curve;
  - the WHAM warm start.
  It is validated on toy targets in the Samplers host.
- *Not in scope:*
  - BASS's adapter and the Friedman qualification (Phase 2);
  - SPC's adapter (Phase 3);
  - any change to `Temperature`, `ParallelTempering`, SMC or any engine;
  - `docs/Maths` and the checklists, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then the sources, held locally (the design's "Sources, held locally"):
- Liang, Liu and Carroll 2007, on SAMC: its algorithm, gain sequence and convergence theorem;
- Belardinelli and Pereyra's 1/t schedule;
- Bornn et al.'s PAWL, for chains that share weights.

In the code:
- *The target contract,* `src/maths/samplers/targets/SamplingTarget.cs`: `ISamplingTarget`, `ITargetEvaluation` (`LogPrior`, `LogLikelihood(b)`), `Temperature`, and `SamplingTarget`'s function-built targets.
- *The chains and the shell:*
  - `src/maths/samplers/rjmcmc/` (`ReversibleJumpChain`, `IRjMove`, `MoveStep`);
  - `src/maths/samplers/ensemble/` (`ChainEnsemble`, `IEnsembleModel`, `ConvergenceVerdict`, `DrawSample`, `RunningMoments`);
  - `src/maths/samplers/blocks/` (`BlockScan`, `RandomWalkMetropolis`, `SliceSampler`).
- *The tempering pilot WHAM reads:* `ParallelTempering` (`TraceRungs`, `RungTraces`) and the self-sized ladder's sizing in `TemperingLadder`.
- *The evidence's shapes,* `src/maths/samplers/evidence/LadderEvidence.cs`, which the ln Z(b) curve sits beside.

## The design, as Phase 1 builds it

- **The target wrapper** (`MulticanonicalTarget`): wraps any `ISamplingTarget` with a weight function.
  - Its evaluation keeps the inner `LogPrior`, and gives ln w(ℓ₁) as `LogLikelihood(1)`, where ℓ₁ is the inner evaluation's `LogLikelihood(1)`. Other b are refused.
  - So π_w ∝ prior·w(ℓ₁) is sampled at `Temperature.Posterior` by any kernel that runs on a target.
  - The evaluation also exposes ℓ₁, for the learning and the reweighting.
- **The weights** (`MulticanonicalWeights`): immutable.
  - A window [ℓ_lo, ℓ_hi] with K nodes, and ln w piecewise linear in ℓ₁ between them.
  - Outside the window, ln w continues with slope b_lo below and slope 1 above, continuous at the edges.
- **SAMC's learning** (`SamcLearning`).
  - At each step, each chain's visit moves the nodes' log weights toward equal occupancy, with the gain γ_t = t₀/max(t₀, t), Liang et al.'s form.
  - With piecewise-linear weights, the visit's indicator is spread over its two nodes by their interpolation weights. Derive the update, state it, and check that it converges on the toys.
  - The chains share one weight function, and their visits update it in a fixed order, so the result does not depend on scheduling.
  - It runs in burn-in only; then the weights freeze (D10's policy). Its state rides the checkpoint.
- **The run** (`MulticanonicalRun`), over the ensemble shell.
  - The burn-in learns, from the warm start or from flat.
  - Production runs independent chains on the frozen weights, with D10's verdict on π_w's functionals.
  - The caller supplies a factory for its chain on a given target, so BASS's chain plugs in at Phase 2.
- **The reweighting.**
  - Each production draw carries the importance weight e^{b·ℓ₁}/w(ℓ₁), self-normalized at any b in [b_lo, 1].
  - Report the weights as a weighted sample at b = 1, its Kish ESS, and the weighted laws at chosen b.
- **The ln Z(b) curve.**
  - Z(b) = E_w[e^{b·ℓ₁}/w(ℓ₁)] / E_w[1/w(ℓ₁)], valid when the prior is normalized, on a grid of b.
  - Its Monte Carlo error is the production chains' spread.
  - It is a record beside `LadderEvidence`'s estimates.
- **The warm start** (`Wham`): multiple-histogram reweighting (Ferrenberg and Swendsen 1989) of a tempering pilot's rung samples (b_i and their ℓ₁ draws) into ln g on the nodes, iterated to its fixed point. It also sets the window from the pilot's ℓ₁ quantiles at its hot end and at b = 1, with margins.
- **The diagnostics:**
  - the histogram's flatness over the window;
  - tunnelling: the round trips between the window's lowest and highest bins, per chain;
  - the learning's trace of the weights;
  - the b = 1 weights' Kish ESS.
- **Configuration** is a values-only record (`FlatHistogramSpec`): the window or the warm start's source, K, t₀, the learning and production lengths, and b_lo, validated.

## Working rules

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Two chips work beside you.** The ridge chip is in BASS's production files and `src/maths/regression/gprior/`; the Friedman reference test is in `tests/regression/`. Your files do not overlap theirs. If a build fails in a file you do not own, wait and retry; never edit it.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-fh\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-fh\` when you finish.
- **Your files:** new files under `src/maths/samplers/flathistogram/` (namespace `Maths.Samplers.FlatHistogram`), new tests under `tests/samplers/`, and this brief's Reports section. Change no existing file. If you need a change outside these, name it in the report and do not make it.
- **Commit your own files only, by explicit path,** one commit per step, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's untracked notes and the other chips' uncommitted edits stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.

## The steps

Each step lands with its checks and ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next step. Name the passes and their fixes in the report.

1. **The weights and the wrapper.**
   - *Check:* ln w's continuity at the nodes and at the window's edges; the wrapper's density equal to the inner prior plus ln w(ℓ₁), with −∞ off the support; refusal at b ≠ 1.
2. **SAMC's learning and the run.**
   - *Toy 1,* a double well: a uniform prior on an interval and ℓ₁(x) = −c·(x² − 1)², whose g is known by quadrature. Check:
     - the learned ln g against the exact one, up to a constant, within a stated tolerance after the learning;
     - the occupancy flat within its tolerance;
     - convergence improving with the learning's length, as SAMC's theorem says.
   - The weights shared across chains give the same result under any parallelism bound.
3. **The reweighting and ln Z(b).**
   - *Check on toy 1:* the reweighted laws at three b against quadrature; ln Z(b) on a grid against quadrature, within the production chains' errors.
   - *Toy 2,* a small q-state Potts lattice: periodic, prior uniform over the qᴺ configurations, ℓ₁ = J·(the number of agreeing neighbour pairs), and single-site Metropolis moves.
     - Pick q and the size so that the exact energy distribution, by enumeration (qᴺ up to about 10⁷), is bimodal near the transition, and report the bimodality.
     - Check the learned ln g against the exact density of states, ln Z(b) against the exact one at every b on the grid, and the reweighted energy distribution at the transition's b against the exact one.
   - *Report, as exploration:* the tunnelling rate on toy 2, against a parallel tempering ladder's round-trip rate at equal work.
4. **The warm start.**
   - *Check:* on both toys, WHAM's ln g from a parallel-tempering pilot against the exact one within the pilot's errors; learning from the warm start reaching the flat-start's accuracy in fewer sweeps (report both).
   - *Then* the configuration record, its validation and JSON round trip, and resume across the learning and production phases, bit for bit.

## Review before reporting

Take a fresh-eyes pass and fix what it finds:
- *Correctness:*
  - the wrapper's density;
  - SAMC's update with interpolated nodes;
  - the gain sequence;
  - the weights frozen in production;
  - the reweighting's self-normalization;
  - Z(b)'s ratio estimator, which needs a normalized prior (assert it on the toys);
  - nothing existing changed.
- *Coverage:*
  - a window too narrow, which the chain leaves by the Boltzmann extensions;
  - a window too wide;
  - one chain;
  - a chain that never tunnels within its run, which the diagnostics flag;
  - resume mid-learning.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the API and the configuration;
- SAMC's update with interpolated nodes, as derived;
- each check with its tolerance;
- toy 2's exploration readings;
- the suites with their counts (Samplers, and Regression to show nothing moved);
- the review passes and their fixes;
- anything Phase 2 must know.

## Addition (2026-10-30): a mean-field Potts toy

Asked by the user after toy 2's tunnelling reading, and sent to this chip by message; it runs after the run weights and the normalization correction. Toy 2's 2× in tempering's favour is uninformative at nine sites. Tempering's crossing cost at a first-order transition grows exponentially with size, while the multicanonical cost should grow only polynomially (Berg and Neuhaus 1992), and at N = 9 neither shows. Torpid mixing of parallel and simulated tempering is proven for mean-field Potts with q ≥ 3 (Bhatnagar and Randall, SODA 2004; Woodard, Schmidler and Huber, *EJP* 2009).
- *The model.*
  - Curie–Weiss Potts with q = 3 (and 4 if cheap), the prior uniform over configurations, and ℓ₁ = (J/(2N))·(Σ_k n_k² − N).
  - The state is the occupation vector n. Single-site Metropolis lumps onto it exactly; derive and check that the Hastings ratio reduces to the likelihood ratio.
  - The exact ln g sums the multinomial over the n with one Σ n_k².
  - J places the transition near b_t = 0.6. b_t at each N is the exact law's equal-weight point, checked against the derived closed form βJ = 2(q − 1)·ln(q − 1)/(q − 2).
- *The checks,* in the suite: at two or three N, the learned ln g, the reweighted law at b_t and ln Z(b), each against the exact values within a stated tolerance.
- *The reading* (exploration, behind an environment variable):
  - N from about 50 to a few thousand.
  - The flat-histogram run, a self-sized tempering ladder and a fixed-b chain at b_t, at equal single-site proposals.
  - One crossing measure for all: round trips of ℓ₁ between the exact law's two peaks at b_t, per 10⁶ proposals.
  - Per N: the peak-to-dip ratio, the crossing rates with their spread, and the fitted growth (exponential against power law).
  - What would distinguish the hypotheses is written down before the run.
- *Not in scope:* a 2D lattice past enumeration, BASS, defaults and `docs/Maths`.

## Reports

Newest first.

### Phase 1 built (chip, 2026-10-29)

**What landed**, on `main`, nothing pushed, every file new under `src/maths/samplers/flathistogram/` and `tests/samplers/`; no existing file changed:
- `d7f4ab8`, step 1: `MulticanonicalWeights`, `MulticanonicalTarget`, `MulticanonicalEvaluation`.
- `3627a13`, step 2: `SamcLearning` (with `FixedPoint`), `IMulticanonicalChain`/`IMulticanonicalModel`/`MulticanonicalChain`, `MulticanonicalRun`, `MulticanonicalDraws`, `FlatHistogramDiagnostics`, `FlatHistogramSpec`.
- `54f417b`, step 3: `Reweighting` (`WeightedDraws`, `ReweightedEstimate`, `LogEvidenceCurve`).
- `2fa3bc0`, step 4: `TemperingPilot`, `Wham` (`WhamEstimate`), the run's pilot, the spec's validation, JSON and resume.
- `4c34930`, after the source check below: WHAM weighs each rung by its effective draws, as Ferrenberg and Swendsen.

**The API.**
- `MulticanonicalWeights(low, high, lowSlope, logWeights)`: ln w at K equally spaced nodes on [ℓ_lo, ℓ_hi], piecewise linear, slope b_lo below and 1 above; `LogWeight(ℓ)`, `TryLocate(ℓ, out segment, out fraction)`, `FlatOccupancy()` (π: 1/(K−1) inside, half at the edges), `Flat(...)`, `WithLogWeights(...)`. Immutable.
- `MulticanonicalTarget<TState>(inner, weights)`: `Evaluate` keeps the inner `LogPrior` and gives ln w(ℓ₁) as `LogLikelihood(1)`, refusing other b (`NotSupportedException`); `MulticanonicalEvaluation` exposes `UntemperedLogLikelihood` (ℓ₁), `LogWeight`, `Weights`, and `Reweighed(weights)`, which re-reads ln w from the cached ℓ₁. The run installs each new weight set on the target (internal), and each chain reweighs its cached evaluation.
- `SamcLearning(start, gainScale t₀, gainStart)`: `Update(levels)` (one iteration, the chains' ℓ₁ in chain order), `Gain(t₀, t)`, `NextGain`, `Save()`/`SamcLearningState`; `FixedPoint(shape, levels, logMasses, out empty)`, SAMC's limit for a known density of states (Newton on the convex G below, from the histogram start).
- The chain contract: `IMulticanonicalModel<TState>` (functional declarations as `IEnsembleModel`, `StartChain(target, chainSeed, readoutSeed)`, `Fingerprint`, `RestoreChain(target, state)`) and `IMulticanonicalChain<TState> : IChain<TState>` (`CurrentEvaluation`, `Reweigh()`, `Accumulate(draw, functionals)`, `SaveState()`). `MulticanonicalChain<TState>` adapts a `ReversibleJumpChain` on the target with a `DrawReduction<TState>`.
- `MulticanonicalRun.Run(spec, target, model, masterSeed, pilot?, resumeFrom?, iterations?, checkpoint, parallel?)` → `MulticanonicalResult<TState>`: `Weights`, `Learning`, `LearningDiagnostics`, `WeightTrace` (start, each hundredth, end), `Production` (the shell's `EnsembleRun`, D10's verdict on π_w's functionals), `Draws`, `ProductionDiagnostics`, `Chains`, `Checkpoint`, `Complete`. Learning: every chain steps at the current weights (in parallel, or serially at a bound of 1), then one update from their visits in chain order, then every chain reweighs. Production: `ChainEnsemble.Run` with burn 0, the learned chains found by their seeds (the shell's `SeedTree` split), so the shell's checkpoint carries production.
- `FlatHistogramDiagnostics`: `Occupancy`, `Desired`, `InWindow`/`BelowWindow`/`AboveWindow`, `RoundTrips` per chain (lowest bin = node 0's half-spacing cell and below, highest likewise; a trip leaves the lowest, reaches the highest, returns), `LastEnd`; `Flatness`, `MaxFlatnessDeviation`, `ChainsWithoutRoundTrip`, `TotalRoundTrips`.
- `Reweighting`: `Weigh(draws, b)` → pooled self-normalized weights e^{b·ℓ₁}/w, chain-major, with pooled and per-chain Kish ESS; `Expectation(draws, b, functional | DrawFunction)` → pooled estimate, error sd(chain estimates)/√C, chain estimates; `LogEvidence(draws, grid)` → `LogEvidenceCurve` (pooled ratio estimator, chains' spread, per-chain curves, per-b Kish ESS, the denominator's ESS, `LogRatio(from, to)`).
- `TemperingPilot(b[], levels[][], inefficiencies?)` (g_i ≥ 1 per rung, null for 1) or `FromLadders(ladders)` (traced for the same scans, likelihood path, slopes; every hotter rung checked for log ratio = (b_c − b)·slope, so a target not linear in b is refused; g_i = max(1, n/ESS) from each rung's ℓ₁ chains); `Wham.Estimate(pilot)` → `WhamEstimate` (distinct pooled ℓ₁, ln ĝ at each, f_i, iterations, residual, effective draws n_i/g_i); `Wham.WarmStart(pilot, shape)` = `FixedPoint` on ĝ; `Wham.Window(pilot, q, margin)`; `Wham.Nodes(pilot, window)`.

**The configuration**, `FlatHistogramSpec` (values only, `Validate()`, JSON round trip, `NeedsPilot` computed and not serialized): `Chains` (4), `LearningIterations` (required, ≥ 0), `ProductionDraws` (required, ≥ 1), `Start` (`Flat` | `Wham`), `Window` (null: from the pilot), `LowSlope` b_lo (null: the pilot's hottest rung), `Nodes` K (null: ⌈width/sd₁⌉ + 1, sd₁ the b = 1 rung's ℓ₁ sd), `GainScale` t₀ (null: 10·K, inside Liang et al.'s 2K–100K), `GainStart` (null: 0 from flat, the pilot's draws over the chains from WHAM), `WindowQuantile` q (0.01), `WindowMargin` (1 sd). A pilot window runs from the hottest rung's q-quantile less a margin to the b = 1 rung's (1 − q)-quantile plus a margin, held inside the pilot's pooled extreme draws (beyond them the extensions carry the hot end's and the posterior's own laws). Lengths are totals; `iterations` takes a call's share, and `resumeFrom` continues.

**SAMC's update with interpolated nodes, as derived.** With θ = −ln w, a visit at ℓ inside the window, in segment j at fraction u, has hat values φ_j = 1 − u, φ_{j+1} = u, and the iteration's update over C chains is θ_k ← θ_k + γ_t·(1/C)·Σ_c 1_in(ℓ_c)·(φ_k(ℓ_c) − π_k), γ_t = t₀/max(t₀, t), t counted from `GainStart` + 1. A visit outside the window contributes nothing. Since Σφ = Σπ = 1, each update conserves Σθ. The mean field is P_in(θ)·(E_θ[φ | in] − π) = −P_in·∇G for the convex G(θ) = ln ∫_in prior·e^{−θ·φ(ℓ)} + π·θ, so the limit is the unique (to a constant) θ* with E_θ*[φ_k | in] = π_k: a function of g over the window only, whatever the extensions, equal to ln g up to a constant when ln g is piecewise linear on the nodes. Two alternatives were rejected: counting outside visits on the edge nodes, whose limit would depend on the extensions, and leaving them in with zero indicator, which drifts every θ. Empty nodes (no mass under their hat) fall without bound while the rest take π_k + d, Liang et al.'s d. For a lattice energy on the nodes, θ* = ln g − ln(π + d) exactly; the edges' half share is the trapezoid's.

**The checks, with their tolerances.** Toy 1 is the double well (x uniform on [−2, 2], ℓ₁ = −10·(x² − 1)², window [−30, −0.02], b_lo = 0.05, K = 31, t₀ = 310, a Gaussian walk of step 0.25). Its references are a 400 000-cell midpoint quadrature of the prior and, for ln g, the closed form. Toy 2 is the 6-state Potts model on a periodic 3 × 3 lattice, J = 1.5: 6⁹ ≈ 1.0·10⁷ configurations enumerated with site 0's colour fixed, a node on each level n = 0…18, b_lo = 0, single-site Metropolis.
- *Step 1:*
  - ln w equals its node values and is continuous at every node and both edges (to 1e−6 at h = 1e−9), with slopes b_lo and 1 outside.
  - ℓ = −∞ gives −∞ when b_lo > 0, and ln w_0 at b_lo = 0.
  - The wrapper's density is the inner prior + ln w(ℓ₁); off the support it is −∞ with no likelihood call. b = 0.5, b = 0 and the slope are refused.
  - `Reweighed` makes no inner evaluation.
  - A reversible-jump chain on the wrapper over ten states reproduces prior·w(ℓ₁) (400 000 draws, within 6 sd at an inflation of 20).
- *Step 2:*
  - The gain's values.
  - The update by hand: spread, averaged, outside visits inert, Σ ln w conserved to 1e−14.
  - `FixedPoint`: recovers a piecewise-linear ln g to 3e−12 (tolerance 1e−6); on levels at the nodes it equals ln g − ln(π + d) to 1e−9, with the empty node interpolated.
  - Double well, learned against θ*: centred RMS 0.170 / 0.090 / 0.043 at 4k / 16k / 64k iterations. That is about 1/√L, and the ordering is asserted and registered for the seed sweep. Tolerance 0.15 at 64k.
  - Against the closed-form ln g at the 24 nodes more than 2.5 spacings from ℓ = −c and 1.5 from the top: θ* 0.008 (tolerance 0.1), learned 0.064 (tolerance 0.15). g is singular at the barrier's top and at the wells' floors, where θ* departs from ln g at the nodes by up to 2.4 nats.
  - Production on the frozen weights, 4 × 40 000 draws:
    - the learned weights' own occupancy is flat to 0.114 (exact by quadrature, tolerance 0.2);
    - production's visits match that occupancy to 0.088 (tolerance 0.2; about one over the root of the ~630 round trips);
    - every chain tunnels (154–164 round trips).
  - The same weights, draws and means under parallelism bounds 1, 2 and unbounded, bit for bit.
  - One chain: RMS 0.081 (tolerance 0.3).
  - A chain stepping 1e−7 is flagged in both phases.
  - A narrow window [−12, −2] is left on both sides by the extensions (production 22 648 in, 20 095 below, 37 257 above), learns its own θ* to 0.024 and is flat to 0.056.
  - A wide window [−30, 6] has its six empty nodes 1–6 rising without bound (ln w 39 against −9), the rest at θ* with d to 0.051, and the histogram unflat (> 0.99).
  - A scripted chain gives two round trips exactly.
  - The trace holds the start, every iteration and the end.
- *Step 3:*
  - The weights by hand (to 1e−14) and Kish's ESS.
  - Double well, 8 chains × 40 000 draws, every bound |z| < 4.5 in the chains' spread errors (7 dof):
    - E[ℓ₁], E[x²] and P(x > 0) at b = 0.1, 0.5 and 1, against quadrature: |z| ≤ 1.74. P(x > 0) at b = 1 is 0.5004 ± 0.0179, the two wells held equally across a 10-nat barrier.
    - ln Z(b) on b = 0, 0.02, …, 1 against quadrature: |z| ≤ 3.99. The prior's masses sum to 1, asserted.
    - That seed's 3.99 led to a 12-seed calibration: RMS z 1.19 for the denominator, 1.10 for ln Z(0.05) and 1.40 for ln Z(1), against t₇'s 1.18. There is no bias; the deviation's sign flips across seeds.
    - The narrow window's laws stay exact through its extensions.
    - One chain gives NaN errors.
    - b = 1 Kish ESS 7 347 of 320 000.
  - Potts:
    - The enumeration sums to 6⁹, with n = 18 holding 6 configurations, n = 14 holding 270, and n = 13, 15–17 empty.
    - At equal peak heights, b_c = 0.6146 (βJ = 0.9218): P(6) = 0.1117 = P(14), dip P(12) = 0.0471, so peak/dip is 2.37. The ordered mass is 0.211, with the gap 13, 15–17 between.
    - Learned against θ*: 0.040 RMS over the 15 occupied levels (tolerance 0.1). θ* equals ln g − ln(π + d) to 3e−15.
    - ln Z(b) on 13 b against the enumeration: |z| ≤ 0.94. The energy law at b_c: worst |z| 1.41 over the occupied levels, exactly 0 at the empty ones.
    - b = 1 Kish ESS 45 636 of 800 000.
- *Step 4:*
  - WHAM on exact rung samples (4 rungs, 10⁵ draws each in exact proportion) recovers ln Z(b_i) − ln Z(0) to the printed 1e−4 (tolerance 2e−3), residual < 1e−8.
  - Eight independent pilots of 8 rungs × 20 000 scans each give warm starts against θ*:
    - double well: per-pilot RMS 0.152 (tolerance 0.25), mean z² over nodes 2.37 (tolerance 3; t₇ gives 1.4);
    - Potts: 0.056 (tolerance 0.25), mean z² 1.24.
  - Warm against flat, the mean of 3 learning seeds:
    - double well: flat 0.139 / 0.120 / 0.098 at 8k / 16k / 32k; warm 0.116 / 0.093 / 0.073 at 0 / 8k / 16k;
    - Potts: flat 0.123 / 0.127 / 0.093 at 25k / 50k / 100k; warm 0.050 / 0.092 / 0.053 at 0 / 25k / 50k.
    - Asserted: warm at half the length beats flat at the full length. Registered for the seed sweep.
  - A run from the pilot alone (2 ladders):
    - window [−49.58, −0.000], K = 68, b_lo = 0.05, t₀ = 680, gain from t = 40 000;
    - ln Z at b = 0.05, 0.5, 1 within 4.5 errors (−0.4566 ± 0.0070 against −0.4600, −1.5712 ± 0.0162 against −1.5677, −1.9504 ± 0.0211 against −1.9440).
  - Refused: a pilot whose ℓ_b = ln(1 + b(e^{ℓ₁} − 1)), an untraced or powered ladder, no rung at b = 1, duplicate rungs, a NaN draw, a one-draw rung.
  - The spec: each of 15 ill-formed fields refused, by `Validate` and by `Run`; JSON round trip equal and stable; `NeedsPilot` kept out.
  - Resume: a run taken at 700 + 1 300 learning iterations, then 400 production draws, then the rest, ends bit for bit where the whole run ends. That covers weights, trace, both phases' diagnostics with `LastEnd`, draws, functionals, means, variances, acceptance and verdict. Another seed, spec or model fingerprint is refused, and so is a finished checkpoint.

**Toy 2's exploration readings.** At 800 000 single-site proposals each:
- the flat-histogram run (100 000 learning and 100 000 production iterations × 4 chains) made 760 round trips per 10⁶ production proposals between the window's end bins (n = 0 and n = 18), 304 in all;
- an 8-rung likelihood-path ladder from b = 1 to 0, adapted in its burn-in to 1.00, 0.82, 0.73, 0.65, 0.55, 0.41, 0.23, 0.00, made 1 640 per 10⁶, 656 hot-to-cold-and-back in all, at swap acceptance 0.68.

Tempering wins by about 2× here: the transition at N = 9 is weak (peak/dip 2.4), the window spans the whole energy range at one site per step, and the two counts' ends differ, so it is a reading, not a claim. The bimodality that defeats tempering needs a lattice past enumeration.

**The suites**, at `2fa3bc0` on a clean tree, built into `artifacts\chip-fh\`:
- *Samplers:* 251 passed, 0 failed, in 1.3 min. That is 206 existing tests and 45 new: weights 9, learning and run 15, reweighting 5, Potts 4, WHAM 7, spec and resume 4, and the sweep's 1, which runs only under `SEED_SWEEP`. At `4c34930`: 253 passed, 0 failed, the 2 new being WHAM's.
- *Regression:* 1 303 passed, 0 failed, in 9 m 54 s. Nothing moved.

**The review passes and their fixes.**
- *Step 1:* nothing in the code. A test reached the slope's default refusal without the interface cast. The target's internal `Install` waited for step 2's runs to be covered.
- *Step 2:*
  - `FixedPoint`'s Hessian was asymmetric: the rank-one fill of later rows overwrote the tridiagonal cross terms, so Newton returned θ = 0. The fill order is fixed.
  - Its line search spun at rounding level for 200 iterations, 112 s on the narrow window. It now starts from the histogram estimate ln(M_k/π_k) and stops when a step no longer lowers G; 0.14 s.
  - The wide window has 6 empty nodes, not 5: node 0's hat sees ℓ just below 0.
  - The flatness check was split into the learned weights' exact occupancy and production's sampling of it, after the single reading sat at 0.143 against 0.15.
  - Added: the closed-form ln g check, the scripted round trips, the trace's shape, the adapter's refusals, and `FixedPoint`'s param tags.
- *Step 3:*
  - The Potts transition read at equal phase masses put b_c (0.718) where the disordered bump is barely a peak (peak/dip 1.17). It is read at equal peak heights now (2.37).
  - The double well's ln Z at small b sat 3–4 errors high on its seed. The 12-seed calibration above found no bias and calibrated errors; the seed was kept.
  - The Potts test now asserts the prior's normalization.
- *Step 4:*
  - The window test checked the hot rung's own minimum, but the window holds inside the pilot's pooled extremes.
  - Warm against flat failed as single-seed comparisons at the default gain start. They are now seed-averaged, at half the length, and registered; see the gain-start note below.
  - `NeedsPilot` serialized into the JSON: now `[JsonIgnore]`, asserted.
  - A fan-out per learning iteration dominated the toys' time: the run steps serially at a bound of 1, and the toys run so, except where the bound is the subject.
  - WHAM collapses identical levels, which took Potts' pilots from 66 s to 2 s.
  - The enum's cref to `Wham` resolved to the enum member: now qualified.
- *Final, against the brief's list:*
  - Correctness, each covered above: the wrapper's density, the interpolated update, the gain, frozen production weights (every production chain's evaluation and the draws hold the learning's final weights, the same object), self-normalization, and the ratio estimator's normalized prior (asserted on both toys).
  - Coverage: narrow and wide windows, one chain, a chain that never tunnels, and resume mid-learning, each with its test.
  - Nothing existing changed: the four commits add 19 files and modify none.

**WHAM checked against Ferrenberg and Swendsen 1989.** Checked after step 4 against the local copies: the letter at `codex-scientiae/supellex-private/PDFs/SPC/FerrenbergSwendsen1989_physrevlett.63.1195.pdf` and the longer paper at `PDFs/SPC/FerrenbergSwendsen1.4822862.pdf`.
- *What matches:* `Wham.Estimate` solves their multiple-histogram equations (3) and (4) in the list form the letter allows for, with draws in place of histograms. It finds the f_n by the derivative-guided iteration they recommend, here Newton on the convex objective whose stationary point is (3)–(4).
- *One departure, now corrected (`4c34930`):* they weigh run n by 1/g_n, with g_n = 1 + 2τ_n its statistical inefficiency (the letter's (3), the longer paper's (10)–(15)). That combination minimizes the error in W(S). Step 4's `Wham` had set g_n = 1, as for independent draws.
  - `Wham` now counts each draw 1/g of its rung and enters each rung with n/g effective draws.
  - `FromLadders` reads g = max(1, n/ESS), ESS the mean's effective sample size over the rung's ℓ₁ chains (`ChainDiagnostics.MeanEss`). A pilot built from arrays takes g = 1 unless given its own.
  - *Checks:* exact rung histograms stay exact under unequal g (to 1e−3); a common g on every rung changes no f (to 1e−9); `FromLadders`' g matches `MeanEss` directly.
- *What it bought on the toys, a reading:* the double well's g run 1.8–2.4, and 4.3–5.4 at the hot rung. Over the same 8 pilots, the per-pilot RMS to θ* is 0.164 with FS's weights against 0.152 with g = 1 (mean z² 2.41 against 2.37). On Potts it is 0.055 against 0.056 (z² 1.16 against 1.24).
  - The warm-start comparisons still pass at the new default: double well warm 0.141 / 0.107 / 0.080 at 0 / 8k / 16k against flat 0.098 at 32k; Potts warm 0.039 / 0.083 / 0.069 at 0 / 25k / 50k against flat 0.093 at 100k.
  - The weights did not lower the error here. One g per rung, read from ℓ₁'s mean, is FS's own summary of a run, not each bin's. The gain they promise needs rungs that mix at very different rates, which Friedman's band, where tempering stalls, may supply; Phase 2 can read it there.

**For Phase 2.**
- *BASS's structure prior is unnormalized.* The Friedman reference hand-over gives D02's count weight, with ln C 0.9531 and ln Z_F −0.194 on Friedman. The ratio estimator then reads ln Z(b) off by the prior's log mass, a constant in b. Add the known normalizer, or read `LogEvidenceCurve.LogRatio`, which cancels it.
- *The denominator's tail.* E_w[1/w]'s variance is finite only if E_prior[e^{−b_lo·ℓ₁}] is finite over the low extension. With b_lo > 0 and ℓ₁ unbounded below, read `NormalizerEffectiveSampleSize`; b_lo = 0 removes the question.
- *The pilot must be on the linear form.* `TemperingPilot.FromLadders` refuses a ladder whose log ratios are not (b_c − b)·slope, which is D05's integrated form; BASS's pilot runs on ℓ_b = b·ℓ₁. The window is read at the pilot's hottest rung, which also sets b_lo.
- *BASS's chain.* BASS's chain is a `BlockScan` over a `ReversibleJumpChain`, so it needs its own `IMulticanonicalChain`:
  - `Reweigh` is `SetState(Current, evaluation.Reweighed(target.Weights))`;
  - `SaveState` carries the scan's snapshot;
  - every block must leave π_w invariant, that is, evaluate through the chain's target at b = 1; a τ block on BASS's own density would not.
- *The gain's start is regime-dependent.* t_w = pilot draws / chains is visit-count parity. On the toys it over-trusts the double well's pilot (whose error matches flat SAMC's at about t = 9k, not 40k) and under-trusts the Potts pilot (worth about t = 115k).
  - A gain too large for a good start degrades it first: Potts 0.050 → 0.092 at 25k → 0.053 at 50k.
  - A gain too small cannot fill what the pilot never sampled, Friedman's dip between phases above all; filling a δ-nat hole costs about t_w·(e^{δK/t₀} − 1) iterations.
  - The qualification should read `GainStart` with t₀ (10·K, from Liang et al.'s range, not tuned).
- *Overhead, measured.* b = 1 Kish ESS is 2.3% of production draws on the double well (31 nodes) and 5.7% on Potts (19). Production's per-node error runs at about one over the root of the round trips.
- *Memory.* Production keeps every draw's ℓ₁ and all of its functionals, for the reweighted laws at any b. BASS's per-grid functionals may want a monitored subset.
- *The posterior's weights.* `Reweighting.Weigh(draws, 1).Weights` is chain-major in draw order, the pooled order a `DrawWeights` array would need. Align it with BASS's `DrawStore`.
- *D10's verdict.* On the double well's production (4 × 40 000 draws) the verdict reads NotConverged: x mixes only by tunnelling, about 160 round trips per chain.
- *A file outside my set I would change:* `tests/samplers/ComparativeSweep.cs`, to hold this chip's three registered claims. They sit in the new `FlatHistogramSweep.cs`, rerun by the same `SeedSweep` under `SEED_SWEEP`.
