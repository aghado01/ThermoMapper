# BASS adaptation prerequisites review

Date: 2026-10-18. Status: source review complete; repairs and runtime qualification pending. This document authorizes no source-code changes. R13 and R14 added 2026-10-19; the work is sequenced in the [preparation brief](bass-preparation-brief.md) and the [BAPS and BGP chip brief](baps-bgp-baseline-repairs-brief.md).

Establish a sound existing regression engine before integrating BASS. Pre-existing defects must have independent reproductions, repairs or explicit statistical resolutions, and verification before the new capability is used to assess the shared infrastructure. Keep existing-engine repairs as separately scoped work, with their own evidence, so a BASS discrepancy can be attributed to the new implementation rather than an unresolved baseline defect.

This review checks the [regression bugs digest](discussion/gemini-regression-bugs-digest.md) against the current implementation and reconciles its modernization claims with the BASS source. It supplies the prerequisite register for the [full BASS port brief](bass-full-port-brief.md). The [public inventory and backend review](bass-public-inventory-and-backend-review.md) remains the detailed account of BASS functionality and R implementation issues.

## Evidence and scope

- ThermoMapper source baseline: `c82c6368b151d28d14dfa1f2076de633efd125b1`.
- BASS reference: version 1.3.1, commit `26710146c623a258854162367a939041efedd091`, checkout `D:/aghado01/packages/BASS`.
- Assessment method: source tracing, mathematical derivation, and the cited primary literature. Source anchors refer to these baselines.
- No new executable probes completed for this review: creation of an isolated check in the designated `.codex` scratch directory was blocked. No probe files were created. The closure checks below are requirements, not reported passing results.

Eight of the digest's ten existing-engine entries describe behavior still present. Some are implementation defects; others are statistical-contract discrepancies or acknowledged approximations. The GP scratch-buffer defect has been repaired. The BGP entry combines a real bandwidth-support problem with an unsupported claim that coordinate scaling breaks a theorem.

The baseline scope is the reviewed Regression, Samplers, LinAlg, Distributions and RNG paths, together with callers affected by their repairs. This is a bounded engine qualification exercise, not a requirement to repair every unrelated ThermoMapper subsystem.

## Existing engine findings

| ID | Finding and assessment | Source and consequence |
|---|---|---|
| R01 | BARS count prior: confirmed measure discrepancy if `p(k)` is intended to be the actual knot-count prior. | [SplineTarget](../../src/maths/regression/spline/bars/SplineTarget.cs#L51) adds only the count prior; [birth](../../src/maths/regression/spline/bars/KnotMoves.cs#L50) includes `-log(k+1)`. Under a constant likelihood, the implied count marginal is proportional to `p(k)/k!`. Resolve the conditional location density before treating BARS as a qualified statistical reference. |
| R02 | BAPS `MarginalEvidence`: confirmed improper upper tail under the conditions below. | [PenalizedSpline](../../src/maths/regression/spline/baps/PenalizedSpline.cs#L81) supplies a profile objective that approaches a finite limit; [BapsSampler](../../src/maths/regression/spline/baps/BapsSampler.cs#L197) samples it against a flat prior on log lambda. Specify a proper sampling target or separate evidence optimization from Bayesian sampling. |
| R03 | Hierarchical P-spline: confirmed stale conditional statistic. | [HierarchicalPSpline](../../src/maths/regression/spline/baps/HierarchicalPSpline.cs#L118) accumulates deviations from the old population coefficients, redraws those coefficients at line 133, then updates the deviation variance from the stale sum at line 136. Repair the scan and verify the joint target. |
| R04 | Robust BARS weights: confirmed approximation, already acknowledged in the main regression documentation. | [WeightResampling](../../src/maths/regression/spline/bars/WeightResampling.cs#L48) uses conditional-mean residuals and a MAD scale. These are plug-in values, not the conditioning variables of an exact joint posterior draw. An exact robust path needs a joint-model derivation; retained approximations need explicit semantics. |
| R05 | Shared adaptive kernel: confirmed for the supplied-instance path. | [BarsEnsemble](../../src/maths/regression/spline/bars/BarsEnsemble.cs#L146) shares `_kernel` unless `adaptKernelTau` creates private instances; [AdaptiveBurn](../../src/maths/regression/spline/bars/BarsEnsemble.cs#L369) adapts any adaptive-kernel instance. Parallel chains can mutate the same object. Supplying the instance also bypasses the explicit adaptive-tempering exclusion. |
| R06 | GP shared scratch: historical defect, fixed in the inspected code. | [GpRegression.Fit](../../src/maths/regression/bgp/GpRegression.cs#L80) allocates its scratch matrix and factor locally. Protect this behavior with concurrent-versus-serial checks using immutable inputs and a stateless kernel. Do not schedule a repair for a scratch buffer that is no longer shared. |
| R07 | Weight resampling disrupts adaptation: confirmed for the default schedule. | [StepOnce](../../src/maths/regression/spline/bars/BarsEnsemble.cs#L295) replaces the chain every 25 steps; [AdaptiveBurn](../../src/maths/regression/spline/bars/BarsEnsemble.cs#L385) reads live segment counters at 200-step boundaries. Complete default windows therefore see zero counters. Other schedules can produce invalid differences across replaced segments. |
| R08 | Singular-design handling: confirmed inconsistent admission. | [WeightedNormalModel](../../src/maths/regression/spline/bars/WeightedNormalModel.cs#L61) rejects a factor that hits its floor. [PoissonModel](../../src/maths/regression/spline/bars/PoissonModel.cs#L180), [WhittleModel](../../src/maths/regression/spline/bars/WhittleModel.cs#L69), and [GibbsLossModel](../../src/maths/regression/spline/bars/GibbsLossModel.cs#L85) continue through floored factors and can return finite scores. A floor flag is a numerical warning, not a complete rank-revealing analysis. |
| R09 | R-hat stopping: confirmed false-convergence path in shared code. | [ChainDiagnostics](../../src/maths/samplers/ensemble/ChainDiagnostics.cs#L40) correctly returns infinity for frozen chains with different means. [ChainEnsemble](../../src/maths/samplers/ensemble/ChainEnsemble.cs#L131) converts infinity and NaN to zero for stopping. [BARS readout](../../src/maths/regression/spline/bars/BarsEnsemble.cs#L206) applies the same conversion to its displayed fit diagnostics. |
| R10 | BGP scaling: confirmed support/validation issue; theorem-failure claim rejected as stated. | [EmpiricalBayesBandwidthPrior](../../src/maths/regression/bgp/EmpiricalBayesBandwidthPrior.cs#L67) computes a lower bound without checking it is below the fixed upper bound 1. The [local theory note](../../lean/Protolemmata/BARS/bgp-twcd-lemmas.md#L67) poses a transport question; it does not establish failure of the approximation rate. |
| R13 | BARS Normal marginal: the `a ≤ 0` rejection detects rounding error, not a model condition. Added 2026-10-19. | [WeightedNormalModel](../../src/maths/regression/spline/bars/WeightedNormalModel.cs#L74) computes `a = yᵀWy − (m/(m+1))·Q`, which loses digits to cancellation, then rejects nonpositive values as a singular fit. In exact arithmetic `a = ‖r‖²_W + Q/(m+1) > 0`. Compute that form from an explicit residual, and state any near-interpolation limit as an admission rule under R08. See the [derivations review §5](bass-derivations-and-integration-review.md#5-d06-proposal-admission-by-angle-with-the-residual-as-shared-state). |
| R14 | BARS span coverage: π(T) is documented as a probability but counts covering spans. Added 2026-10-19 from the [additional-bugs review](gemini-additional-bugs.md). | The coverage loop in [BarsEnsemble](../../src/maths/regression/spline/bars/BarsEnsemble.cs) adds 1 per covering significant span, so the overlapping spans of nested peaks give π(T) > 1 and `PeakModes`' 0.5 cutoff stops being a probability threshold. Count each grid point once per draw, or rename the quantity. A BARS readout repair; it does not gate BASS. |

### Count priors and conditional location measures

For ordered continuous knots on the unit interval,

\[
S_k=\{0<\xi_1<\cdots<\xi_k<1\},\qquad \operatorname{vol}(S_k)=1/k!.
\]

The current density is proportional to `L(y | k, xi) p(k)` with respect to Lebesgue measure on each simplex. Consequently,

\[
p_{\mathrm{current}}(k\mid y)\propto p(k)\int_{S_k}L(y\mid k,\xi)\,d\xi.
\]

Only for a constant likelihood does this reduce to `p(k)/k!`. Drawing `k ~ p(k)` and then sorting independent uniform locations instead gives conditional location density `k!` on `S_k`. That model needs the corresponding factor in its target. Other location restrictions require their own normalizers.

This is a BARS model-resolution task. The generic RJ engine's move-selection correction does not supply a missing model density. Existing tests on a Poisson target over an integer dimension do not establish the count law for continuous sorted knots. BASS needs a separate derivation for its own term representation and measure; copying a factorial correction mechanically would be another error.

### BAPS tail and Gibbs conditioning

Let `nu` be coefficient dimension, `r` the penalty nullity, and `A = Z'Z + lambda P`. If the design identifies the null space and the limiting penalized residual sum of squares is positive, then

\[
\log|A|=(\nu-r)\log\lambda+O(1),\qquad
\ell_p(\lambda)\longrightarrow C.
\]

For `u = log(lambda)`, the sampled density therefore approaches a positive constant as `u` increases and cannot normalize. The defect is not confined to nearly polynomial data. A declared proper prior or bounded support can define a different, proper target; an undocumented numerical cap does not resolve the model specification. Gibbs and AdaptiveGibbs use different variance-component priors, so matching their lambda distributions is not a correctness criterion for the current marginal mode.

For R03, using an old value is not inherently an invalid Gibbs step. The actual problem is that the population coefficient draw uses the old deviation variance, while the subsequent deviation-variance draw uses residuals formed with the old population coefficients. These two dependent updates are not the advertised alternating full-conditionals. Recomputing the deviation statistic after the population draw is one repair direction; a different ordering must be checked against the full joint density.

For R04, drawing coefficients alone is insufficient because the scale is also plugged in. Derive the complete joint model, including any dependence of the coefficient prior on the weights, before asserting a Gamma or mixture-responsibility conditional. A correctly derived collapsed transition is also possible; collapse is not itself an approximation. Preserve the documented approximation only with an explicit name, purpose and uncertainty interpretation, and keep it out of exact-posterior reference comparisons.

### State ownership and diagnostic failure

R05 is a race across ensemble chains. Replicas within one tempering ladder advance sequentially, although they also share the supplied kernel. Per-chain or per-replica ownership must be chosen explicitly, and the adaptive-tempering rule must apply regardless of how a kernel is supplied. Store tuning parameters in declarative configuration and construct live mutable kernels at the execution boundary.

R07 needs cumulative counters or correctly reconciled segment statistics. The implementation already banks move statistics for reporting; adaptation currently ignores those banks. Test complete and partial adaptation windows and resampling periods that both divide and do not divide the window. Do not infer successful adaptation solely from a nonzero final acceptance rate.

R09 also affects a single chain or a first batch containing too few draws: the diagnostic returns NaN, which the stopping rule accepts as zero. Invalid or unavailable diagnostics must not certify convergence. A truly constant functional requires an explicit policy; it must not be conflated with chains frozen at different values. This applies to both stopping and reported diagnostics. More sophisticated split/rank-normalized diagnostics are a subsequent improvement to a correct basic failure policy.

### BGP scaling qualification

The support issue is concrete. With two one-dimensional inputs at 0 and 100 and the default `gamma1 = 0.25`, the code computes `Tn = 100`, so the lower bound is 2500 while the upper bound remains 1. No bandwidth has positive prior density. Define and validate coordinate units, the shared training/prediction transform, and a nonempty finite support before starting the chain.

The mathematical assessment of fixed affine normalization is different. A fixed nonsingular affine map sends a smooth compact embedded manifold to another smooth compact embedded manifold and preserves smoothness order under composition. Applying the approximation result to that transformed manifold preserves the powers of bandwidth and intrinsic dimension while changing geometric constants. This is a deduction from the theorem's hypotheses, not a claim that the paper explicitly analyzes the repository's preprocessing. [Tang, Wu, Cheng and Dunson, Proposition 4.2](https://arxiv.org/html/2407.09286v3#S4.SS1).

Empirical min-max transforms can depend on sample size, and degenerating coordinate spans require further control. Record that qualification; do not state that anisotropy alone breaks the theorem. The current theory sketches are not compiled formal verification.

## Additional constraints on reuse

| ID | Finding | Required treatment |
|---|---|---|
| R11 | [Dense Cholesky](../../src/maths/linalg/Cholesky.cs#L58) silently floors pivots and eagerly constructs the inverse triangular factor. | This behavior is unsuitable as the unchanged factor backend for BASS structural acceptance. Add a strict solve-oriented path with explicit failure status, then incremental maintenance. Audit existing callers' numerical contracts; preserve intentional regularization where it belongs. Neither all existing calls nor all rank-deficient Bayesian models are thereby declared invalid. |
| R12 | [IObservationModel.SampleCoefficients](../../src/maths/regression/spline/bars/IObservationModel.cs#L35) defaults to the conditional mean, inherited by Weighted Normal. [BARS readouts](../../src/maths/regression/spline/bars/BarsEnsemble.cs#L312) use it for curves and extrema. | Current default intervals reflect knot variation through conditional-mean curves, omitting conditional coefficient uncertainty. The mean remains a valid Rao-Blackwellized estimator of the linear fitted value under a correct target. Nonlinear readouts such as argmax and quantiles require conditional draws or a valid integration method. Resolve existing result semantics and require full BASS posterior uncertainty. |

R08 and R11 require a clean ownership split: LinAlg reports numerical facts; Regression specifies whether the selected prior and likelihood admit the model. A proper coefficient prior can make a rank-deficient design well-defined. Adding an unreported ridge or pivot floor to an unregularized g-prior does not do so coherently.

## Corrections to BASS modernization claims

| Claim from the digest | Verified disposition | Port requirement |
|---|---|---|
| Full factorization costs `O(M^3)` at every proposal. | [getQf](../../../packages/BASS/R/mcmc_funcs.R) in the external checkout performs full Cholesky on candidates that reach it. The R implementation already updates relevant Gram entries incrementally. | Keep a strict full-factorization reference; then add `O(M^2)` bordered insertions and Givens deletions, with rebuild checks. A border extension changes dimension and is not a same-size rank-one update. Basis evaluation and cross-products retain observation-dependent costs. |
| Schur complement at or below `1e-12` must always be rejected. | Exact dependence is a structural problem for the unregularized g-prior; a universal absolute tolerance is not scale-aware. The R code instead catches Cholesky failure and uses a diagonal-ratio heuristic. | Specify relative rank/conditioning checks, finite-value checks, and rebuild-before-rejection where numerical drift is plausible. Rejected work must leave accepted state unchanged. |
| Scalar hinge branches mispredict about half the time, requiring `Vector256.Max`. | The actual `pos` helper in `R/basis_funcs.R:13` uses `(abs(vec)+vec)/2`. The stated branch rate is unsupported. | Treat SIMD as a measured optimization with scalar fallback and hardware dispatch, not a statistical prerequisite or mandatory platform restriction. |
| Replace categorical membership with a `ulong`. | The current implementation uses `%in%`; a word mask is a useful representation for up to 64 encoded levels. | Support larger category sets, stable encoding and declared unseen-category behavior. Preserve full categorical capability. |
| Uniform proposals yield less than 1% acceptance for at least 20 predictors. | The threshold is unsupported. `R/rjmcmc_all.R:11` weights interaction orders and multi-variable selections using current counts; the one-variable branch is uniform. `R/mcmc_funcs.R:12` computes weighted unordered subset probabilities through combinatorial sums. | Benchmark actual proposal families and subset-probability cost. Compare posterior error and effective exploration per unit time, not acceptance alone. |

The external BASS source paths above are resolved through the checkout link in the references below; their line numbers are pinned by the reviewed commit.

The proposed statistical enhancements need separate treatment:

- Residual-informed moves can preserve the chosen target when their complete forward and reverse probabilities, selection probabilities and support are correct. Marginal `X' residual` screening can miss interaction-only signals, so preserve exploratory support. Locally balanced proposals require an actual normalized proposal construction; residual weighting alone is not the construction. [Zanella](https://arxiv.org/abs/1711.07424).
- Sparse Dirichlet priors change the model. Linero's construction concerns regression-tree splitting proportions; applying it to BASS subsets, especially selection without replacement, requires a new derivation. Proposal weights and prior probabilities must remain distinct. [Linero](https://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1264957).
- Intrinsic dimension does not bound interaction order in the original predictor coordinates. A hard gate on `maxInt` restricts model support. Treat it as an explicit model extension or use dimension information only to guide proposals without eliminating valid terms.

These enhancements are outside baseline repair. None is a cure for an incorrect target, stale conditional, or invalid convergence decision.

## Existing engine remediation register

The following work belongs to an existing-engine repair track, separate from BASS feature implementation. Separation preserves ownership and attribution; it does not defer the reviewed baseline defects until after BASS integration. All batches are pending. R06 is a historical repair to protect, and the theorem-failure assertion in R10 requires correction rather than a speculative kernel redesign.

| Batch | Findings and owner | Required closure evidence |
|---|---|---|
| B01 Diagnostics | R09; Samplers owns diagnostic/stopping mechanics, Regression owns its result presentation. | Independent fixtures for differing frozen chains, identical constants, insufficient draws, one chain and nonfinite observations; no invalid diagnostic can trigger convergence; reported values preserve the distinction. Run affected sampler and regression checks. |
| B02 Mutable state and adaptation | R05 and R07; BARS constructs its private kernels, Samplers provides general execution/state contracts. | Explicit ownership per chain/replica; consistent adaptive-tempering policy; deterministic seeded serial/concurrent comparisons; counters and tuning remain valid across multiple rebuild schedules, including the default 25/200 combination. |
| B03 BARS statistical contracts | R01, R04, R08, R12, R13 and R14; Regression. | Flat-likelihood count-law oracle; derived robust target and conditionals or explicitly separated approximation; model-specific rank policy; fixed-structure coefficient mean/covariance checks and qualified uncertainty for nonlinear readouts. Statistical choices and public semantics must be settled before changing their tests. |
| B04 BAPS statistical repairs | R02 and R03; Regression. | A normalizable sampling target with a documented prior, or a separately exposed optimization mode; valid hierarchical update order; small reference posterior/conditional checks and Monte Carlo tolerances. Do not require agreement between intentionally different priors. |
| B05 BGP support and regression protection | R06 and R10; Regression. | Early rejection of empty/nonfinite support; documented input units and transform policy; boundary and constant-coordinate cases; concurrent/serial fit checks protecting the local scratch fix; corrected theoretical wording without claiming a new proof. |
| B06 Existing numerical qualification | R08 and R11; LinAlg reports status, Regression owns admission. | Inventory the factor paths and callers actually reused or changed. Verify supported SPD solves, log determinants, precision-draw covariance and singular/near-singular behavior against independent references. Record intentional regularization explicitly. The new strict/incremental dense backend is completed separately in port stage 2 and must pass before BASS uses it. |

For every batch, record the reproduction or mathematical counterexample, chosen specification, affected paths, repair revision, independent checks, and remaining limitations. A passing legacy suite alone is insufficient when it does not exercise the reported defect. A statistical-contract resolution must make the code, exposed results and documentation agree; merely relabeling an incorrect exact-posterior claim does not qualify that path for exact BASS comparisons.

### Baseline completion gate

1. Establish independent failing checks or counterexamples for the confirmed current defects before introducing BASS or refactoring their shared contracts.
2. Resolve B01-B06 as separately scoped repairs/qualification work. Keep hypothetical model enhancements and unrelated subsystem work outside these batches.
3. Run the affected existing-engine suites and the new independent checks. Record a qualified source revision, toolchain and test results. Until then, the baseline gate is open.
4. Develop and verify new numerical and general-sampler facilities against that baseline. Rerun the affected baseline checks after shared-contract migration. Changes in expected statistical behavior must be explained by the separately recorded repairs.
5. Integrate and test BASS after the existing baseline and the specific shared facilities it consumes are qualified. Use the pinned R implementation only where its output is qualified; use independent mathematical references for known R defects.

The user has explicitly requested this order to prevent pre-existing bugs from confounding the BASS work. Independent specification and oracle design can proceed while repairs are being planned. This is not a claim that the whole repository is bug-free, nor permission to start implementation in this documentation-only task.

## Implications for SMC and dependency boundaries

SMC can reuse general state, weight, resampling and execution machinery after qualification. It cannot normalize an improper final target, restore a missing prior factor, or repair an invalid mutation kernel. BASS owns the target sequence, hyperparameter conditionals and calibration feedback policy. Samplers owns the population driver and method-appropriate diagnostics; LinAlg owns numerical operations. No reverse dependency on Regression or SPC is needed.

When population resampling is added, cloned particles need independent mutable workspaces and RNG streams. State exchange and target changes must refresh cached densities coherently. Particle weight ESS and ancestry diversity remain distinct from independent-chain R-hat and autocorrelation ESS. These are port-stage requirements, not evidence that the current code already supports SMC.

The existing-engine gate precedes integration. The port's D01-D10 decisions still need their own model-specific resolutions; correcting BARS or BAPS does not by itself validate BASS's prior measure, collapsed sweep, tempering path or cut-calibration target.

## References

- [Full BASS port brief](bass-full-port-brief.md).
- [BASS public inventory and backend review](bass-public-inventory-and-backend-review.md).
- [Regression bugs digest](discussion/gemini-regression-bugs-digest.md), [initial adaptation notes](discussion/gemini-bass-adaptation-notes.md), and [BASS 2.0 ideas](discussion/gemini-bass-2.0.md). Preserve these discussion sources as historical inputs.
- [Bayesian Regression documentation](../../docs/Maths/BayesianRegression.md) and [Samplers documentation](../../docs/Maths/Samplers.md).
- External checkout: [BASS source](../../../packages/BASS); specific inspected helpers: [mcmc_funcs.R](../../../packages/BASS/R/mcmc_funcs.R), [basis_funcs.R](../../../packages/BASS/R/basis_funcs.R), [rjmcmc_all.R](../../../packages/BASS/R/rjmcmc_all.R), and [bass.R](../../../packages/BASS/R/bass.R).

Append repair and qualification results to this register or link separately scoped repair briefs/reports here. Keep source-confirmed findings, implemented repairs and passing runtime evidence distinct.
