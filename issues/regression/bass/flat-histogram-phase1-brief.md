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

## Reports

Newest first.
