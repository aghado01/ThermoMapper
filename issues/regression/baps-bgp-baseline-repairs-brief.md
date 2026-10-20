# BAPS and BGP baseline repairs — chip brief

Date: 2026-10-19. Status: complete. The report is below, and the suite passes at `ed60885`. Started after step 0.1 of the [BASS preparation brief](bass-preparation-brief.md) landed on `main`.

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

**Consequence for the reported mean.** A1 and A3 put positive prior density at the polynomial limit, so the posterior density of u = log λ decays only like e^(−u) (A3: e^(−u/2)). Then E[λ] = ∫ e^u p(u | y) du is infinite on any data, and `LambdaMean`, the mean of the λ draws, has no population value. On data with a clear interior maximum, the tail weight is tiny (about e^(−19) for the sine case), so the estimate is stable in practice; near the null space it is not. A2 and B keep E[λ] finite. With A1 I would report the posterior median `LambdaMedian` in place of `LambdaMean`, for all three modes. It sits with `LambdaLo` and `LambdaHi` as quantiles of the same pooled draws, and it reads the same on λ and log λ. That is a clean break in `BapsResult`, and the `BapsTests` assertions on `LambdaMean` would move to the median.

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

**Decisions requested:** (1) R02: A1, A2, A3, B or C, and with A1 or A3, whether `LambdaMean` becomes `LambdaMedian`. (2) R10: T1 or T2, M1 or M2, and confirmation of response centering.

**Decided by the user, 2026-10-19:** R02 = A1, with `LambdaMean` replaced by `LambdaMedian` for all modes. R10 = T2 + M1 with a `None` option, plus response centering.

## Report

For each item, record:

- the reproduction or counterexample;
- the chosen specification;
- affected paths;
- commits;
- independent checks with results;
- remaining limitations;
- requests for the main track.

### Summary

Completed 2026-10-19.

- **Qualified revision:** `ed60885`, built with .NET SDK 10.0.201 in Release.
- **Regression suite:** `scripts/fact-harness.ps1 -Project tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj` passes 117 of 117 (`artifacts/test-runs/TestHarness.Runner/20261019_190012452`).
- **Benchmark:** the env-gated `BgpBenchmark` (`BGP_BENCH=1`) also passes.
- **Parent revision for the reproductions:** `6f481ea`, whose code is identical to `f5afb7a`, unless an item names another.

| Commit | Content |
|---|---|
| `0ea0329` | R03: hierarchical scan fix and its checks |
| `4745e00` | R10: support rejection and warm start |
| `c2f0a1d` | R06: concurrency protection test |
| `9f9a0b3` | R10: theorem wording (docs and BGP-2) |
| `3d57bab` | Item 6: BAPS R̂ in moments form |
| `dc75590`, `e7acce6` | Decision memo and its addendum |
| `1d1cada` | `BapsResult.LambdaMedian` replaces `LambdaMean` |
| `5fe826e` | R02: uniform-edf prior for `MarginalEvidence` |
| `1f10ba0` | R10: fitted transform (similarity inputs, centered response) |
| `ed60885` | Enum theory data, so the fact harness can run each case |

**Closure against the register.**

- **B04 meets its closure evidence.** The sampling target is normalizable, with a documented prior (A1), and is checked against quadrature. The hierarchical update order is valid. Small reference conditional and posterior checks run with Monte Carlo tolerances. No agreement is required between intentionally different priors.
- **B05 meets its closure evidence.** Empty or non-finite support is rejected early. The input units and the transform policy are documented. Boundary and constant-coordinate cases are covered. Concurrent and serial fits are checked against each other. The theoretical wording is corrected without claiming a new proof.

### Item 1: R06 protection (B05)

- **Reproduction.** This is a historical defect, so the check passes at the parent, as it should. To show the check has teeth, the shared scratch matrix was put back temporarily (not committed): the check then failed, with a concurrent log-marginal of 67.20 against the serial −43.91.
- **Specification.** `GpRegression.Fit` is pure over immutable inputs and a stateless kernel. Concurrent results must equal serial results bit for bit.
- **Paths.** `tests/regression/BgpConcurrencyTests.cs`.
- **Commit.** `c2f0a1d`.
- **Checks.** 192 concurrent fits and predictions over 24 bandwidths equal the serial ones bit for bit, and the inputs are unchanged afterwards. Pass.
- **Limitations.** Interleaving needs at least two cores. A future mutable kernel is not covered.

### Item 2: R10 support and transform (B05)

**Reproductions.** At `6f481ea` unless stated:

- **Empty support.** Inputs {0, 100} gave a lower bound of 2500. The chains froze at t = 0.999, where the prior density is zero, with acceptance 0 and a constant R̂ that satisfies a stopping rule, and a prediction was still returned. The same happened with a lower bound of exactly 1 ({0, 2}) and with T_n = ∞ ({0, 10²⁰⁰}).
- **Narrow support.** For {0, 1.9995}, the warm start was clamped to 0.999, below the support (0.9995, 1], so the chains froze there.
- **Degenerate inputs.** Duplicated inputs gave T_n = 0, and the chains drifted to t ≈ 10⁻¹¹ with acceptance 1. A NaN input made every output NaN, silently. A single input failed with an unclear `Math.Clamp` message.
- **Uncentered response** (at `5fe826e`). Adding 20 to the response moved the log-evidence at t = 0.05 from 36.29 to −861.62, and t̂ from 0.092 to 0.542. In the memo probe, offsets 0, 5 and 20 gave t̂ = 0.074, 0.132 and 0.282.
- **Units.** Rescaling the inputs to 100x − 7 emptied the support (lower bound 5.95).

**Specification.**

- **Support rule.** The prior rejects fewer than two inputs, a non-finite γ₁T_n², γ₁T_n² ≥ 1, and T_n = 0. Non-finite inputs and responses are rejected when the transform is fitted. Chains start at √(γ₁T_n²); an explicit start outside the support is rejected.
- **Transform** (user decision: T2 + M1 + centering, with `None`). `GpRegression` fits a `GpTransform`. The default input map is `Similarity`, (x − lo)/s with s the largest span; `None` keeps the inputs as given. The response is centered on the training mean and restored in prediction, and it is not rescaled.
- **Units.** Bandwidths, the support and T_n are in the mapped frame. `BgpResult.InputScale` reports s, and t_original = s²·t. The prior's constructor is now `(gp, …)`, reading T_n from the GP's own distances.

**Paths.**

- `src/maths/regression/bgp/`: `GpTransform.cs` (new), `GpRegression.cs`, `EmpiricalBayesBandwidthPrior.cs` and `BgpSampler.cs`.
- Tests: `BgpSupportTests.cs`, `BgpTransformTests.cs` (new), and the constructor calls in `BgpTests`, `BgpSamplerTests` and `BgpBenchmark`.
- The BGP docs section.

**Commits.** `4745e00` and `1f10ba0`.

**Checks.**

- **Support** (`BgpSupportTests`):
  - empty support;
  - lower bounds of 1, 0.9025 and 0.9995 (the narrow case, where the chains now move inside);
  - three non-finite cases, duplicated inputs and a single input;
  - a constant coordinate, under both scalings;
  - `Similarity` makes {0, 100} valid, but not two opposite corners of a 4-cube.
- **Transform** (`BgpTransformTests`):
  - an offset of 20 leaves the evidence and the bandwidth posterior unchanged and shifts predictions by exactly 20, all to 10⁻⁹;
  - 100x − 7 gives the same T_n, support, evidence, predictions and posterior;
  - the s² bandwidth conversion;
  - inputs already in the unit frame map by the identity, so `None` and `Similarity` agree bit for bit;
  - a NaN response is rejected.
- **End to end.**
  - The Swiss-roll test gives t̂ = 0.184 [0.156, 0.220] and a test MSE of 0.0032 (Var f* = 0.555).
  - In the benchmark, the EB/oracle RMSE ratio is 0.997, 0.997, 0.990 and 0.991 at n = 50, 100, 200 and 400, with RMSE ∼ n^−0.628.

**Limitations.**

- The similarity map is fitted on the sample, so its constants depend on n. The theory for that case stays open (BGP-2).
- Very few points in four or more dimensions can still produce an empty support, which is rejected.
- Anisotropic data keep their anisotropy by design.
- `PredictMean` does not check that test inputs are finite.
- The tests and benchmark still apply a per-coordinate min–max (M2) by hand, so the instrument's similarity map is the identity on their data.

### Item 3: R10 wording (B05)

- **Specification.** A fixed nonsingular affine map, including a diagonal one, keeps the rate exponents. It changes only the constants, which Lemma 4.1 and Prop. 4.2 tie to curvature, second fundamental form, volume, diameter, injectivity radius and reach. Maps fitted on the sample (the similarity or a min–max) depend on n: carrying the rate over needs φₙ → φ with constants uniform near φ, and a degenerating span breaks that. No new proof is claimed.
- **Paths.** The "Coordinate scaling and the theory" paragraph in `docs/Maths/BayesianRegression.md`, and BGP-2 in `lean/Protolemmata/BARS/bgp-twcd-lemmas.md`. The numbering (Lemma 4.1 for the operator expansion, Prop. 4.2 for the RKHS approximation) was checked against arXiv 2407.09286v3 §4.1.
- **Commits.** `9f9a0b3`, updated by `1f10ba0` for the similarity map.
- **Checks.** None apply to the wording; the Lean blocks remain uncompiled sketches.
- **Limitations.** The fixed-affine half is a lemma-tier candidate and is not formalized.

### Item 4: R03 hierarchical scan (B04)

- **Reproduction** (`6f481ea`). On the τ_u² conditional, the first sweep gave mean z = −85.6 (variance z = −16.3), and sweep 5 gave mean z = 4.6 (variance z = 49.3). Against quadrature, the stationary E[σ] was off by z = 8.35 and E[τ_u] by z = −4.58, and the β_0 entries by up to |z| = 4.1.
- **Specification.** The deviation statistic is recomputed after the β_0 draw. Each sweep draws β_r, then β_0, then σ², τ_u² and τ_0², and every draw conditions on the current values of the others.
- **Paths.** `HierarchicalPSpline.cs`, `tests/regression/HierarchicalPSplineConditionalTests.cs`, and the BAPS docs (the scan is described and the caveat removed).
- **Commit.** `0ea0329`.
- **Checks.** On a tiny model (ν = 5, two replicates, eight points, IG(2, 0.2)):
  - the whitened Gaussian conditionals over 4000 draws are within |z| ≤ 3;
  - the IG conditionals at sweeps 1 and 5 are within |z| ≤ 2.1;
  - against quadrature with every coefficient integrated out (16 chains of 20000 draws), every stationary mean is within |z| ≤ 1.3;
  - the existing borrow-strength test still passes (MSE ratio 0.26).
- **Limitations.** The sampler still runs a single chain without R̂ or ESS. The tolerances are 5σ on the conditionals and 4.5 standard errors on the stationary means.

### Item 5: R02 `MarginalEvidence` target (B04)

- **Reproduction.** At `6f481ea`, see the memo: ℓ_p plateaus, and on linear data the chains drift to the NaN wall near λ = 10¹⁴, with R̂ = 1.06. At `1d1cada` (flat prior), the quadrature check fails on linear data: the sampler's ln λ median is 19.2 ± 3.2 against 1.23, and its 97.5% quantile 30.5 against 5.80.
- **Specification** (user decision: A1 with `LambdaMedian`).
  - The mode places a uniform prior on edf, which induces the density |d edf/d ln λ| on ln λ.
  - `PenalizedSpline.EffectiveDegreesOfFreedom` and `EdfSlope` read the generalized spectrum of (P, ZᵀZ + P). It is computed once, lazily, from a band LDLᵀ factor via `SampleInnovation` and `DenseEigen`, with the r null-space eigenvalues set to zero.
  - `BapsResult` reports `LambdaMedian` for all modes.
- **Paths.** `PenalizedSpline.cs`, `BapsSampler.cs`, `tests/regression/BapsMarginalEvidenceTests.cs` (new), `BapsTests.cs`, and the BAPS docs.
- **Commits.** `1d1cada`, `5fe826e` and `ed60885`.
- **Checks.**
  - edf and its slope match a dense reference from λ = 10⁻³ to 10⁹.
  - ∫ |d edf/du| du = 10.000000 = ν − r, so the prior is proper.
  - The sampler matches quadrature (10 runs × 4 chains × 4000 draws). Linear data give ln λ quantiles −1.491, 1.211 and 5.787 against −1.471, 1.225 and 5.801; sine data give −3.869, −2.423 and −1.194 against −3.861, −2.416 and −1.184.
  - On the `BapsTests` sine data, the median is 3.74 [1.25, 9.30] against the REML λ = 5.41.
- **Limitations.**
  - E[λ] is infinite under A1, which is why the median is reported.
  - ℓ_p is formed as yᵀy − β̂ᵀZᵀy, so it loses precision above about λ = 10¹⁰ and returns NaN near 10¹⁴. The prior keeps the chain far below that range.
  - The spectrum requires ZᵀZ + P to be nonsingular. Otherwise the first `MarginalEvidence` evaluation throws `InvalidOperationException`, which `Parallel.For` wraps in an `AggregateException`.
  - The sine case also passes against the old flat target, because the prior varies slowly over the bulk; the linear case is the one that discriminates.
  - Observation only: under `Gibbs`, the IG(10⁻³, 10⁻³) priors make the λ posterior proper only through the 10⁻³ terms.

### Item 6: B01 policy in BAPS display

- **Reproduction** (`6f481ea`). `RHatLogLambda` was a bare double from the sums-based overload: NaN for one chain, with the constant and frozen cases blurred by cancellation.
- **Specification.** One `RunningMoments` per chain, reduced by `ChainDiagnostics.RHat(moments)` to an `RHatDiagnostic`.
- **Paths.** `BapsSampler.cs`, `BapsTests.cs`.
- **Commit.** `3d57bab`.
- **Checks.** One chain reports `Unavailable` (value NaN) and four chains report `Estimated`, under both `Gibbs` and `AdaptiveGibbs`. The mixing assertions now require an estimated state.
- **Limitations.** `HierarchicalPSpline` still reports no R̂.

### Related findings (observed, not repaired)

- **Single-chain `HierarchicalPSpline`.** Confirmed. The new checks take their standard errors from independent seeded runs. A candidate for the ensemble shell.
- **`DifferencePenalty` only.** Confirmed: `BapsSampler` and `HierarchicalPSpline` accept only `DifferencePenalty`. `PenalizedSpline`'s edf spectrum works for any `IBandPenalty`.
- **No λ optimizer in `AnisotropicTensorSpline`.** Confirmed, unchanged.
- **`MarginalEvidence` refactorization.** Unchanged: one factorization per proposal for ℓ_p, and two more per recorded step. The new spectrum could give ℓ_p, β̂ and the edf terms in O(ν²) per step, without refactorizing and with a stable RSS. That would be a separate caching and stability change.
- **BGP acceptance rate includes burn-in.** Confirmed, unchanged; it belongs to main-track step 3.4.

### Requests for the main track

1. **Delete the sums overload (step 3.6).** `ChainDiagnostics.RHat(double[], double[], int)` has no callers left. Its doc line "Remaining caller: `BapsSampler`" is stale.
2. **Factor-status migration (step 3.6).** There is one new `BandCholesky` call site: `PenalizedSpline.PenaltySpectrum` (`DecomposeBanded`, `HitFloor`, `SampleInnovation`). The existing call sites are unchanged: `BapsSampler`, `HierarchicalPSpline`, `PenalizedSpline.Solve`, and `CholeskyDecomposition` in `GpRegression.Fit`.
3. **Baseline (step 0.6).** These public changes alter BAPS and BGP outputs by design, so any baseline recorded before them will differ:
   - `BapsResult.LambdaMedian` replaces `LambdaMean`, and `BapsResult.RHatLogLambda` is an `RHatDiagnostic`;
   - `MarginalEvidence` now uses the edf prior;
   - `GpRegression` gains a `scaling` argument (default `Similarity`) and centers the response;
   - `EmpiricalBayesBandwidthPrior(gp, …)` no longer takes `x`;
   - `BgpResult` gains `InputScale`.
4. **Fact harness.** `scripts/fact-harness.ps1` fails any theory case whose display name contains a string argument, with MSBuild error MSB4177, because the quotes break the generated `dotnet test` arguments. Escalated as a candidate harness fix, per AGENTS.md; the chip's tests now use enum data.
