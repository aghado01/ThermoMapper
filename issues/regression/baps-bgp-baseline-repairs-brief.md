# BAPS and BGP baseline repairs — chip brief

Date: 2026-10-19. Status: chip brief, not started. Start once step 0.1 of the [BASS preparation brief](bass-preparation-brief.md) has landed on `main`. Append the report below when done.

## Context

ThermoMapper is preparing to add BASS to `Maths.Regression`. Before BASS is integrated, the existing engines must be repaired and qualified, so that pre-existing defects cannot be mistaken for BASS defects. The [adaptation prerequisites review](bass-adaptation-prerequisites-review.md) records the findings and groups them into batches B01–B06.

This chip owns **B04 (BAPS)** and **B05 (BGP)**. Neither has a code path into BASS, so this work runs in parallel with the main preparation track, which owns B01, B02, B03 and B06. The chip's closure report is an input to the main track's baseline qualification (its step 0.6).

Read first:

- [AGENTS.md](../../AGENTS.md) and [CONTRIBUTING.md](../../CONTRIBUTING.md).
- From the prerequisites review: R02, R03, R06 and R10 in the findings table; the sections "BAPS tail and Gibbs conditioning" and "BGP scaling qualification"; rows B04 and B05 of the remediation register, with its evidence paragraph.

## Why step 0.1 comes first

Step 0.1 changes how non-finite R̂ is treated in [ChainEnsemble](../../src/maths/samplers/ensemble/ChainEnsemble.cs) and [ChainDiagnostics](../../src/maths/samplers/ensemble/ChainDiagnostics.cs). [BgpSampler](../../src/maths/regression/bgp/BgpSampler.cs#L92) runs on `ChainEnsemble`, and [BapsSampler](../../src/maths/regression/spline/baps/BapsSampler.cs#L160) reports `ChainDiagnostics.RHat`. Starting after that change means the chip's checks run against the repaired diagnostics.

## Scope

### B05: BGP

1. **R06 protection (test only).** Add concurrent-versus-serial checks for [GpRegression.Fit](../../src/maths/regression/bgp/GpRegression.cs#L80) with immutable inputs and a stateless kernel. The local scratch fix is already in place; the check keeps it there.
2. **R10 support.**
   - Reproduce the empty support: two one-dimensional inputs at 0 and 100 with the default `gamma1 = 0.25` give a lower bandwidth bound of 2500 against the fixed upper bound of 1 in [EmpiricalBayesBandwidthPrior](../../src/maths/regression/bgp/EmpiricalBayesBandwidthPrior.cs#L67).
   - Then reject empty or non-finite support before the chain starts, with a clear error.
   - Define and document the input units, and one fitted transform shared by training and prediction.
   - Include the response in that policy: [GpRegression.Fit](../../src/maths/regression/bgp/GpRegression.cs#L98) fits a zero-mean GP to the response as given, so an uncentered response pulls predictions toward 0 away from the data and biases the bandwidth. Put response centering (a fitted offset applied in training and restored in prediction) in the transform-policy memo.
   - Cover boundary and constant-coordinate inputs.
3. **R10 wording.** In the BGP section of [BayesianRegression.md](../../docs/Maths/BayesianRegression.md) and in [bgp-twcd-lemmas.md](../../lean/Protolemmata/BARS/bgp-twcd-lemmas.md#L67), replace the claim that coordinate scaling breaks the approximation theorem:
   - A fixed nonsingular affine map preserves the rates and changes only the constants.
   - Empirical min-max transforms depend on sample size and need further control.

   Claim no new proof.

### B04: BAPS

4. **R03 hierarchical scan.**
   - Reproduce the stale statistic in [HierarchicalPSpline](../../src/maths/regression/spline/baps/HierarchicalPSpline.cs#L118): the deviation sum is formed with the old population coefficients, those coefficients are redrawn at line 133, and the deviation variance is then updated from the stale sum at line 136.
   - Repair by recomputing the deviation statistic after the population draw, or by another order verified against the full joint density.
   - Check each full conditional against its closed form on a tiny model, and the stationary behavior against an independent reference.
5. **R02 `MarginalEvidence` tail.**
   - Reproduce the improper target: on a small design, the profiled evidence of [PenalizedSpline](../../src/maths/regression/spline/baps/PenalizedSpline.cs#L81) approaches a constant as `log λ` grows, and [BapsSampler](../../src/maths/regression/spline/baps/BapsSampler.cs#L197) samples it under a flat prior on `log λ`.
   - Then write a decision memo and stop for the user's choice. Options:
     - a declared proper prior on `log λ`;
     - a declared bounded support, stated as part of the prior;
     - a separate evidence-optimization mode beside the Bayesian sampler.
   - Implement the chosen option. Because the target is one-dimensional in `log λ`, check it against quadrature.
   - `Gibbs`, `AdaptiveGibbs` and `MarginalEvidence` use different priors, so agreement between them is not a criterion.
6. **B01 policy in BAPS display.** Move [BapsSampler](../../src/maths/regression/spline/baps/BapsSampler.cs#L160) to the moments form of R̂: a `RunningMoments` per chain and `ChainDiagnostics.RHat` returning an `RHatDiagnostic`, so `RHatLogLambda` carries its state. The sums-based overload it calls today stays only until then. Step 0.1 already changed `BgpResult.BandwidthRHat` to `RHatDiagnostic`.

### Related findings outside this chip's repairs

Added 2026-10-19 from the [additional-bugs review](gemini-additional-bugs.md). None of these is a reviewed baseline defect, and none touches BASS. Record what you observe about them in the report, but don't repair them here:

- `HierarchicalPSpline` runs a single chain with no R̂ or ESS (already documented). Multi-chain support would need the ensemble shell; it is a candidate after this chip.
- `BapsSampler` and `HierarchicalPSpline` accept only `DifferencePenalty`, not the `IBandPenalty` extension interface.
- `AnisotropicTensorSpline` has no λ optimizer (already documented).
- `MarginalEvidence` mode reassembles and refactors the banded system on every recorded step. If the R02 repair touches this path, keep any caching change in a separate commit.
- BGP's acceptance rate includes burn-in. The fix belongs to the shell (main-track step 3.4), and BGP inherits it.

## Order of work

1. Write the failing checks and reproductions for items 2, 4 and 5, plus the R06 protection test, and record each failure at its parent revision.
2. Send the R02 decision memo. Include the R10 transform policy as well if the public behavior changes.
3. Make the repairs for R03 (item 4) and R10 (items 2 and 3), and apply item 6.
4. Repair R02 (item 5) once the decision is made.
5. Run the suite with `scripts/fact-harness.ps1 -Project tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj`, then write the report.

## Boundaries

**May edit:**

- `src/maths/regression/spline/baps/**` and `src/maths/regression/bgp/**`.
- Their tests under `tests/regression/`. New test files take a `Baps`, `Bgp`, `PenalizedSpline` or `HierarchicalPSpline` prefix.
- The BAPS and BGP sections of `docs/Maths/BayesianRegression.md`, including their bullets under "Status and caveats".
- `lean/Protolemmata/BARS/bgp-twcd-lemmas.md`.

**Leave untouched** (owned by the main track):

- `src/maths/samplers/**`, `linalg/**`, `distributions/**`, `rng/**`, `regression/spline/bars/**` and the shared spline files.
- `docs/Maths/Samplers.md`.

**Factor calls.** Keep the existing `BandCholesky` and `CholeskyDecomposition` calls as they are. The main track adds a status API beside them and migrates these call sites in its cleanup sweep. If a repair needs a change outside these boundaries, record it in the report as a request instead.

## Working rules

- **Check first.** Each repair begins with an independent check or counterexample, and the check and the fix land in one commit.
- **Commits.** One concern per commit, on `main`. The main track works in the same checkout, so stage files by explicit path and commit nothing outside these boundaries.
- **Semantics before tests.** Settle statistical semantics with the user before changing any test expectation or public behavior. Known-incorrect outputs become defect fixtures, not gold standards.
- **Stay in scope.** No BASS work, no unrelated subsystem fixes, and no subagent dispatch.

## Done when

- B04 and B05 meet the closure evidence in the remediation register.
- The regression suite passes at a recorded revision.
- The report below is complete.

## Decision memo: R02 sampling target and R10 transform policy

Written 2026-10-19 by the chip, at parent revision `6f481ea`. Two decisions are needed before the R02 repair and the R10 transform land.

### R02: what `MarginalEvidence` samples

**Reproduction.** On a cubic P-spline with n = 30, ν = 12 and r = 2, the profiled evidence ℓ_p levels off as λ grows:

| log₁₀ λ | 2 | 4 | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|---|---|
| ℓ_p, linear + noise | 8.3394 | 8.4282 | 8.42876 | 8.42860 | 8.4338 | 8.6385 | NaN |
| ℓ_p, sine + noise | −33.8998 | −34.6106 | −34.6175 | −34.6176 | −34.6176 | −34.6165 | −34.4152 |

Beyond about λ = 10⁸ the values drift because of rounding, and at 10¹⁴ the linear case returns NaN. On the linear data, the `MarginalEvidence` chains (4 chains, burn 500) drift upward with run length. The 97.5% quantile of log₁₀ λ is 7.55 after 1500 samples, 14.01 after 6000 and 14.01 after 24000. The chains stop near 10¹⁴ only because the evidence there is NaN, and a NaN comparison rejects the move: the effective upper bound is a floating-point failure. R̂ is 1.06 at 24000 samples, so it does not flag the drift. On the sine data no drift is visible, because the interior maximum is about 19 nats above the plateau, but the target is still improper.

**Options.**

- **A1. Uniform prior on the effective degrees of freedom (recommended).** edf(λ) = tr[(ZᵀZ + λP)⁻¹ZᵀZ] runs monotonically from ν down to r. A uniform prior on edf induces the density on u = log λ of |d edf/du| / (ν − r), with |d edf/du| = tr S − tr S² for S = A⁻¹ZᵀZ. The prior is proper, has no hyperparameter, and is unchanged by rescaling Z, y or P, which only shifts λ. Its upper tail decays like 1/λ, so the posterior normalizes; the lower tail is already controlled by ℓ_p. This is the multivariate form of a uniform shrinkage prior (Natarajan & Kass 2000), with edf as the complexity measure (Hodges & Sargent 2001). Cost: S comes from ν solves with the existing band factor, O(ν²·bw) per step. No linalg change is needed.
- **A2. Normal prior on log λ, mean and sd declared by the caller.** The simplest code. There is no honest default, because the scale of λ depends on the units of Z and P.
- **A3. PC prior** (Simpson et al. 2017). It shrinks toward the polynomial base model (λ = ∞) and needs the caller's tail statement (U, α). Principled, but needs more machinery than A1.
- **B. Declared bounded support**, uniform on log λ over [L, U]. A small code change, but for data near the null space the posterior piles up at U, so summaries depend on U directly. A scale-free U would itself be stated in edf units, which needs the A1 machinery.
- **C. Separate evidence optimization.** Remove `MarginalEvidence` from `BapsLambdaUpdate` (a clean break) and add an explicit REML optimizer that returns λ̂, β̂(λ̂), σ̂² and a flag for when the maximum is the polynomial limit. There is no prior to defend, but the λ uncertainty under the marginal model is lost; only `Gibbs` and `AdaptiveGibbs` would report λ intervals.

**Recommendation: A1.** It keeps the mode's purpose (a λ posterior with β integrated out) without inventing a scale. The check is one-dimensional: posterior mean and quantiles of log λ from the sampler against quadrature of ℓ_p(e^u) + log(tr S − tr S²), using dense S computed independently in the test, on both the linear and the sine data. The mode's prior stays different from the `Gibbs` prior, which is already the documented position. Under `Gibbs`, the IG(10⁻³, 10⁻³) priors make the λ posterior proper only through the 10⁻³ terms; that is outside this chip and is recorded as an observation only.

### R10: input units and the fitted transform

**Current state.** Inputs are used exactly as given. The fixed upper bound t ≤ 1 assumes squared distances of order 1 (TWCD2025 assumes X ⊂ [0,1]^D). Prediction uses the same raw coordinates. The response is not centered, so the zero-mean GP pulls predictions toward 0 away from the data, and the bandwidth absorbs the offset. In a probe (n = 60 on [0,1], y = c + sin 2πx + noise 0.1), adding an offset c moved t̂ from 0.074 (c = 0) to 0.132 (c = 5) and 0.282 (c = 20). In-sample MSE rose from 0.0015 to 0.0101, and the c = 20 prediction at x = 2 was 5.8 instead of about 20.

**Settled already, and implemented under item 2 regardless of this choice.** The prior's constructor rejects inputs whose support is empty or not finite in the sampled coordinate log t: γ₁T_n² ≥ 1, a non-finite T_n (NaN, infinite or overflowing inputs), T_n = 0 (every point duplicated, which leaves log t unbounded below), and fewer than two inputs.

**Options for the transform.**

- **T1. Caller-owned.** `GpRegression` keeps using coordinates and response as given. The documentation states the contract: coordinates scaled so the domain lies in [0,1]^D, and a centered response. A public fitted-transform type is offered, which callers apply to training and test inputs and use to restore predictions. Public behavior is unchanged apart from the new type.
- **T2. Instrument-owned (recommended).** `GpRegression` fits one transform from the training data at construction, applies it in `Fit`, `KernelAffinity` and `PredictMean`, and restores the response offset in predictions. The prior reads T_n from the GP's own distances, so its separate `x` argument goes away. `BgpResult` reports bandwidths in the fitted frame, with the frame's scale for conversion (t_original = s²·t). Public behavior changes: raw test inputs are mapped by the training transform, and bandwidths are in fitted units.

**Input map (either option).**

- **M1. Isotropic similarity (recommended).** x′ = (x − lo)/s, where lo is the per-coordinate training minimum and s the largest coordinate span. The data map into [0,1]^D with their shape preserved. A constant coordinate needs no special case, and the map is nonsingular unless all inputs coincide (already rejected, because then T_n = 0).
- **M2. Per-coordinate min–max.** This is what the current tests do by hand. It is anisotropic, changes the geometry the kernel sees (the BGP-2 question), and needs a rule for zero-span coordinates.

**Response.** Center on the training mean: the fitted offset is subtracted in training and added back in prediction. The response is not rescaled, so σ² stays in response units.

**Recommendation: T2 + M1 + centering,** with a `None` choice for callers whose inputs are already in the documented units. Both M1 and M2 are empirical transforms whose constants depend on n, which the corrected BGP-2 wording records. The tests that normalize by hand (`BgpSamplerTests`, `BgpBenchmark`) produce data that M1 maps by the identity: their training minimum is 0 and their largest span is 1.

**Decisions requested:** (1) R02: A1, A2, A3, B or C. (2) R10: T1 or T2, M1 or M2, and confirmation of response centering.

## Report

For each item, record:

- the reproduction or counterexample;
- the chosen specification;
- affected paths;
- commits;
- independent checks with results;
- remaining limitations;
- requests for the main track.
