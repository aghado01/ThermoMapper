
Here is a synthesis of the scatter notes, observations, recorded bugs, and architectural edge cases concerning the **Bayesian Regression Engine** (`Maths.Regression`) and its underlying samplers (`Maths.Samplers`), gathered from `issues/`, `docs/`, `lean/`, and the test suites in [ThermoMapper](file:///d:/aghado01/ThermoMapper).

---

### 1. Theoretical & Stationary Distribution Mismatches

* **The $p(k)/k!$ Knot-Count Prior Mismatch in BARS**
  * **Source**: [`docs/Maths/BayesianRegression.md:129`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L129-L131), [`issues/grok-spc-applications-crypto.md:38`](file:///d:/aghado01/ThermoMapper/issues/grok-spc-applications-crypto.md#L38-L41).
  * **Observation**: The [SplineTarget](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/SplineTarget.cs) log-target adds only $\log p(k)$ without an explicit $\xi \mid k$ density. However, [KnotBirthMove](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/KnotMoves.cs)'s Hastings auxiliary ratio includes $-\log(k+1)$, which assumes Lebesgue measure on the ordered knot simplex. Because that simplex has volume $1/k!$, the stationary knot-count distribution is actually proportional to $p(k)/k!$ instead of $p(k)$. If an ordered i.i.d.-uniform knot prior was intended, the target is missing $+\log k!$.

* **BAPS `MarginalEvidence` Improper Prior at $\lambda \to \infty$**
  * **Source**: [`docs/Maths/BayesianRegression.md:174`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L174-L188), [`docs/Maths/BayesianRegression.md:272`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L272-L274).
  * **Observation**: In [PenalizedSpline](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/PenalizedSpline.cs), the profiled log evidence $\ell_p(\lambda)$ levels off to an asymptote as $\lambda \to \infty$. Consequently, a flat prior on $\ln \lambda$ cannot be normalized. On nearly linear or polynomial data, the random-walk Metropolis chain on $\ln \lambda$ can drift upward indefinitely.
  * **Discrepancy**: `MarginalEvidence` and `Gibbs` / `AdaptiveGibbs` modes target different posteriors and priors, so reported smoothing parameters ($\lambda$) are not directly comparable.

* **Non-Strict Gibbs Scan in `HierarchicalPSpline`**
  * **Source**: [`docs/Maths/BayesianRegression.md:189`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L189-L191), [`docs/Maths/BayesianRegression.md:274`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L274-L275).
  * **Observation**: In [HierarchicalPSpline.cs](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/HierarchicalPSpline.cs), the variance parameter $\tau_u^2$ update conditions on the population curve $\beta_0$ from *before* that sweep's redraw, resulting in a simultaneous update rather than a valid alternating Gibbs scan.

* **Approximate Gibbs Conditioning in Robust Weight Resampling**
  * **Source**: [`docs/Maths/BayesianRegression.md:147`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L147-L149).
  * **Observation**: [StudentTWeights](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WeightResampling.cs) and [ContaminatedNormalWeights](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WeightResampling.cs) resample weights from the residuals of the conditional posterior mean rather than drawing jointly with coefficients, creating an approximate rather than exact Gibbs block.

---

### 2. Concurrency, State Sharing, & Thread Safety

* **Shared Adaptive Kernel Mutation Across Threads**
  * **Source**: [`docs/Maths/BayesianRegression.md:118`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L118-L121), [`docs/Maths/BayesianRegression.md:275`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L275-L276).
  * **Bug**: If an instance of [AdaptiveLocalBetaKernel](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/AdaptiveLocalBetaKernel.cs) is passed via `BarsConfig.Kernel`, all parallel chains and tempering replicas in [ChainEnsemble](file:///d:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainEnsemble.cs) share and mutate that single kernel instance across threads. Only `BarsEnsemble`'s dedicated `adaptKernelTau` parameter properly derives independent per-chain kernel instances.

* **The `GpRegression.Fit` Shared-Scratch Incident**
  * **Source**: [`issues/spred/spred-retrospective.md:29`](file:///d:/aghado01/ThermoMapper/issues/spred/spred-retrospective.md#L29), [`issues/spred/spred-pilot-thread.md:40`](file:///d:/aghado01/ThermoMapper/issues/spred/spred-pilot-thread.md#L40).
  * **Observation**: Recorded as an infamous regression bug where static/shared mutable scratch buffers in [GpRegression](file:///d:/aghado01/ThermoMapper/src/maths/regression/bgp/GpRegression.cs) produced race conditions and corrupted fits when executed under multi-threaded parallel blocks. Serves as a standing warning against static scratch in mathematical routines.

* **Freezing Adaptation via Weight Resampling**
  * **Source**: [`docs/Maths/BayesianRegression.md:119`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L119-L121).
  * **Bug**: [AdaptiveLocalBetaKernel](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/AdaptiveLocalBetaKernel.cs) tunes $\log \tau$ based on chain relocate attempt counters. When weight resampling is active, the chain is rebuilt every `WeightEvery` steps (default 25), resetting the counters. As a result, the 200-step adaptation window sees no new attempts, and $\tau$ remains frozen at its initial value.

---

---

### 4. Model Invariants & Diagnostics Caveats

* **Singular Design Handling Discrepancy**
  * **Source**: [`docs/Maths/BayesianRegression.md:133`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L133-L146).
  * **Observation**: [WeightedNormalModel](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WeightedNormalModel.cs) detects rank deficiency and band Cholesky floor violations, returning $-\infty$. Conversely, Laplace/BIC models ([PoissonModel](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/PoissonModel.cs), [WhittleModel](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WhittleModel.cs), [GibbsLossModel](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/GibbsLossModel.cs)) only guard against $m \le \nu$; when the Gram matrix is singular, they return finite estimates instead of $-\infty$.

* **Gelman–Rubin ($\hat{R}$) Early Stopping Pass on Divergent Chains**
  * **Source**: [`docs/Maths/Samplers.md:118`](file:///d:/aghado01/ThermoMapper/docs/Maths/Samplers.md#L118-L119), [`docs/Maths/Samplers.md:168`](file:///d:/aghado01/ThermoMapper/docs/Maths/Samplers.md#L168-L169).
  * **Bug**: In [ChainEnsemble.cs](file:///d:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainEnsemble.cs), early stopping treats any non-finite $\hat{R}$ as $0$ (converged). If chains freeze at different local modes ($W \approx 0$ while $B > 0$), $\hat{R} = +\infty$, which mistakenly triggers the early stopping check.

* **BGP Coordinate Normalization & Bandwidth Boundaries**
  * **Source**: [`docs/Maths/BayesianRegression.md:206`](file:///d:/aghado01/ThermoMapper/docs/Maths/BayesianRegression.md#L206-L215), [`lean/Protolemmata/BARS/bgp-twcd-lemmas.md:67`](file:///d:/aghado01/ThermoMapper/lean/Protolemmata/BARS/bgp-twcd-lemmas.md#L67-L107).
  * **Observation**: [EmpiricalBayesBandwidthPrior](file:///d:/aghado01/ThermoMapper/src/maths/regression/bgp/EmpiricalBayesBandwidthPrior.cs) hardcodes an upper bandwidth bound of $1$, assuming inputs are min-max scaled. However, min-max scaling is anisotropic per coordinate, which warps the isotropic squared-exponential kernel into an anisotropic distance metric and breaks the theoretical RKHS approximation rate on manifolds.

---

### 5. BASS Modernization ("BASS 2.0") Notes & Pitfalls

* **Source**: [`issues/gemini-bass-adaptation-notes.md`](file:///d:/aghado01/ThermoMapper/issues/gemini-bass-adaptation-notes.md), [`issues/gemini-bass-2.0.md`](file:///d:/aghado01/ThermoMapper/issues/gemini-bass-2.0.md).
* **Deficiencies in Legacy BASS (R Package) to avoid in C# port**:
  1. **$O(M^3)$ Matrix Factorization**: Legacy BASS calls full `chol(X'X)` at each proposal. The C# engine must use $O(M^2)$ bordered Cholesky rank-1 updates on births and Givens downdates on deaths.
  2. **Collinear Basis Detection**: Proposals with $a^T a - r^T r \le 10^{-12}$ cause numerical breakdown and must be rejected speculatively.
  3. **Hinge Branch Mispredictions**: Evaluating $[s(x-t)]_+$ with scalar conditional jumps creates pipeline stalls on ~50% mispredictions; branchless SIMD (`Vector256.Max`) is required.
  4. **Categorical Evaluation**: Replacing linear search/%in% branching with 64-bit integer bitmasks (`ulong`).
  5. **Combinatorial Stalls**: Uniform knot and variable proposals lead to $<1\%$ acceptance in $p \ge 20$. Proposed modernization uses Zanella (2020) locally balanced proposals ($X^T r$ residual correlations), intrinsic dimension gating for `maxInt`, and Linero (2018) sparse Dirichlet priors $\text{Dirichlet}(\alpha/p)$.