Viewed baps-bgp-baseline-repairs-brief.md:1-100
Viewed bass-preparation-brief.md:1-256
Viewed Samplers.review.md:1-35
Viewed BayesianRegression.review.md:1-51
Viewed Samplers.md:1-183
Viewed BayesianRegression.md:1-250
Viewed BayesianRegression.md:251-292
Listed directory regression
Viewed bass-adaptation-prerequisites-review.md:1-120
Viewed bass-adaptation-prerequisites-review.md:121-156
Viewed ExactChangepoint.cs:170-185
Viewed BarsEnsemble.cs:330-345
Viewed BarsEnsemble.cs:420-460
Viewed PoissonModel.cs:100-145
Viewed ChainDiagnostics.cs:80-115
Viewed ChainDiagnostics.cs:116-172
Viewed ReversibleJumpChain.cs:45-75

A comparison of the four deliverable documents ([`Samplers.review.md`](file:///d:/aghado01/capstone/runs/run_1/deliverable-review/Samplers.review.md), [`BayesianRegression.review.md`](file:///d:/aghado01/capstone/runs/run_1/deliverable-review/BayesianRegression.review.md), [`BayesianRegression.md`](file:///d:/aghado01/capstone/runs/run_1/deliverables/docs/BayesianRegression.md), and [`Samplers.md`](file:///d:/aghado01/capstone/runs/run_1/deliverables/docs/Samplers.md)) against the two briefs ([`baps-bgp-baseline-repairs-brief.md`](file:///d:/aghado01/ThermoMapper/issues/regression/baps-bgp-baseline-repairs-brief.md) and [`bass-preparation-brief.md`](file:///d:/aghado01/ThermoMapper/issues/regression/bass-preparation-brief.md)) reveals several bugs, statistical nuances, and design considerations that are **not covered** in the briefs.

While the briefs deliberately restrict their scope to the prerequisite repair batches B01–B06 (findings R01–R13) required for BASS integration, the deliverables and their reviews surface substantive issues across other areas of `Maths.Regression` and `Maths.Samplers`.

---

### 1. Exact Changepoint Regression (Completely Omitted from Both Briefs)

Neither brief touches [`ExactChangepoint`](file:///d:/aghado01/ThermoMapper/src/maths/regression/changepoint/ExactChangepoint.cs) or [`SegmentModels`](file:///d:/aghado01/ThermoMapper/src/maths/regression/changepoint/SegmentModels.cs), which reside under `Maths.Regression`:

* **Default $\sigma^2$ Formula Discrepancy ([`BayesianRegression.review.md:17`](file:///d:/aghado01/capstone/runs/run_1/deliverable-review/BayesianRegression.review.md#L17)):** 
  In [`ExactChangepoint.DefaultHyperparameters`](file:///d:/aghado01/ThermoMapper/src/maths/regression/changepoint/ExactChangepoint.cs#L175-L183), the default noise variance is calculated as the lower median of squared first differences:
  $$\sigma^2 = \max\left(\frac{\operatorname{med}_i[(y_{i+1}-y_i)^2]}{0.9098},\, 10^{-6}\right)$$
  The documentation originally claimed it was the square of the median difference $\operatorname{med}(\Delta y)^2 / 0.9098$ (which collapses to zero whenever differences change signs).
* **Dominant Quadrature Complexity ([`BayesianRegression.review.md:23`](file:///d:/aghado01/capstone/runs/run_1/deliverable-review/BayesianRegression.review.md#L23)):** 
  [`CauchySegmentModel`](file:///d:/aghado01/ThermoMapper/src/maths/regression/changepoint/SegmentModels.cs#L97) evaluates segment evidence via $G = 512$ quadrature nodes on $\nu \pm 6\rho$. Filling the $\frac{n(n-1)}{2}$ segment tables costs $O(n^2 G)$ operations, dominating the $O(Kn^2)$ dynamic program unless $K^2 > G$ ($K \gtrsim 23$).
* **Quadrature Mass Truncation:** 
  The Cauchy quadrature abruptly truncates posterior mass outside $\nu \pm 6\rho$ without reporting boundary truncation or integrating tails.

---

### 2. BARS Extrema and Readouts ([`BarsEnsemble`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/BarsEnsemble.cs), [`SplineExtrema`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/SplineExtrema.cs))

* **Span Coverage $\pi(T) > 1.0$ from Nested Peaks ([`BayesianRegression.review.md:21`](file:///d:/aghado01/capstone/runs/run_1/deliverable-review/BayesianRegression.review.md#L21)):**
  In [`BarsEnsemble.cs:330-333`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/BarsEnsemble.cs#L330-L333), `SpanCoverageSum` is incremented by 1 for *each* covering significant peak in a draw. If a dominant peak's drop-level span contains a secondary significant peak whose span covers the same grid points, the draw contributes $+2$ (or more), causing the average $\pi(T)$ to exceed 1.0:
  * In code comments and docs, $\pi(T)$ was documented as a "posterior probability".
  * In [`BarsEnsemble.PeakModes`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/BarsEnsemble.cs#L436-L440), mode bounds are scanned outward while $\pi(T) \ge 0.5$. When nested spans artificially inflate $\pi(T)$, the 0.5 threshold no longer reflects a 50% posterior probability mass cutoff.
  * Neither brief addresses this nesting artifact or clarifies how `PeakModes` should extract spans from multi-count coverages.
* **Algorithmic Discontinuity in Root Finding:**
  [`SplineExtrema`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/SplineExtrema.cs#L57-L71) switches algorithms based on degree: degree $\le 3$ uses exact 4-point cubic recovery with closed-form root formulas, whereas degree $> 3$ switches to Vandermonde interpolation and Durand–Kerner iterative polynomial root solving with fixed stopping tolerances.

---

### 3. Non-Normal Observation Models ([`WhittleModel`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WhittleModel.cs), [`GibbsLossModel`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/GibbsLossModel.cs), [`PoissonModel`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/PoissonModel.cs))

* **Missing Posterior Coefficient Draws for Whittle & GibbsLoss:**
  In [`bass-preparation-brief.md:106`](file:///d:/aghado01/ThermoMapper/issues/regression/bass-preparation-brief.md#L106), step 0.5c (R12) repairs [`WeightedNormalModel.SampleCoefficients`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WeightedNormalModel.cs) to draw $\sigma^2$ and $\beta$. However:
  * [`WhittleModel`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WhittleModel.cs) and [`GibbsLossModel`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/GibbsLossModel.cs) inherit the default from [`IObservationModel.SampleCoefficients`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/IObservationModel.cs#L35), which returns the conditional mean $\hat\beta$.
  * Nonlinear readouts (peaks, quantiles, intervals) under Whittle and GibbsLoss will continue to represent knot uncertainty only, leaving their interval coverage inconsistent with Normal.
* **Arbitrary Ridge and Acceptance Threshold in `PoissonModel` ([`BayesianRegression.review.md:20`](file:///d:/aghado01/capstone/runs/run_1/deliverable-review/BayesianRegression.review.md#L20)):**
  [`PoissonModel.SampleCoefficients`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/PoissonModel.cs#L100-L128) accepts the first Laplace draw if the log MH ratio exceeds $-10$ (`suspectThreshold`); otherwise it runs up to 5 independence MH steps and can fall back to returning $\hat\beta$ directly. Furthermore, its log posterior target applies an undocumented $N(0, 100)$ ridge penalty (`PoissonModel.cs:141`).
* **IRLS Divergence / Failure Modes:**
  While step 0.5b (R08) addresses rank admission from Cholesky status for all four models, the briefs do not specify policies for non-convergence or step-halving divergence during IRLS fits in [`PoissonModel`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/PoissonModel.cs) and [`WhittleModel`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WhittleModel.cs).

---

### 4. BAPS Architecture and Usability Constraints

* **`IBandPenalty` Extension Seam Incompleteness ([`BayesianRegression.review.md:22`](file:///d:/aghado01/capstone/runs/run_1/deliverable-review/BayesianRegression.review.md#L22)):**
  [`IBandPenalty`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/IBandPenalty.cs) is exposed as an extension interface, but only [`PenalizedSpline`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/PenalizedSpline.cs) accepts it. [`BapsSampler`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/BapsSampler.cs) and [`HierarchicalPSpline`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/HierarchicalPSpline.cs) accept only concrete [`DifferencePenalty`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/DifferencePenalty.cs) instances because they access `.Order`, `.Roughness`, and `.SquaredDifferencesInto`. Neither brief considers unblocking this seam.
* **`HierarchicalPSpline` Lacks Multi-Chain and Convergence Diagnostics:**
  While [`baps-bgp-baseline-repairs-brief.md:38`](file:///d:/aghado01/ThermoMapper/issues/regression/baps-bgp-baseline-repairs-brief.md#L38) repairs the stale conditional update (R03), [`HierarchicalPSpline`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/HierarchicalPSpline.cs) still runs as a single chain with no $\hat R$ or ESS convergence diagnostics, unlike `BapsSampler`, `BarsEnsemble`, and `BgpSampler`.
* **Missing Optimizer in `AnisotropicTensorSpline` ([`BayesianRegression.md:277`](file:///d:/aghado01/capstone/runs/run_1/deliverables/docs/BayesianRegression.md#L277)):**
  [`AnisotropicTensorSpline`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/AnisotropicTensorSpline.cs) calculates REML log evidence over the eigen-spectrum of 2D penalties but provides no optimizer or search routine for $(\lambda_x, \lambda_y)$.
* **Precomputation and Allocation Asymmetry in `MarginalEvidence` ([`BayesianRegression.review.md:25`](file:///d:/aghado01/capstone/runs/run_1/deliverable-review/BayesianRegression.review.md#L25)):**
  While Gibbs modes in `BapsSampler` cache $Z^\top Z$ and $Z^\top y$, `MarginalEvidence` mode re-evaluates via `PenalizedSpline`, re-assembling and re-factoring the banded system 3 times per recorded step.

---

### 5. Bayesian GP (BGP) Considerations

* **Zero-Mean GP on Uncentered Response:**
  [`GpRegression`](file:///d:/aghado01/ThermoMapper/src/maths/regression/bgp/GpRegression.cs#L80-L103) implements a zero-mean GP with SE kernel $k_t(d^2) = \exp(-d^2/2t)$ directly on raw $y$ without centering. When $y$ has a non-zero mean, out-of-sample predictions revert to 0 away from the data, and bandwidth likelihood estimation is biased toward bandwidths that compensate for non-centered responses. While the brief addresses input scaling (R10), it does not address response centering.
* **Burn-in In Acceptance Rate:**
  Like `ChainEnsemble`, [`BgpResult.AcceptanceRate`](file:///d:/aghado01/ThermoMapper/src/maths/regression/bgp/BgpSampler.cs) is computed from cumulative counters that include burn-in proposals.

---

### 6. General Samplers Infrastructure ([`Maths.Samplers`](file:///d:/aghado01/ThermoMapper/src/maths/samplers/))

* **$O(Cn^2)$ Quadratic Bottleneck in `ChainDiagnostics.Ess` ([`Samplers.md:141`](file:///d:/aghado01/capstone/runs/run_1/deliverables/docs/Samplers.md#L141)):**
  In [`ChainDiagnostics.cs:150-157`](file:///d:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainDiagnostics.cs#L150-L157), the variogram loop computes $\hat\rho_t$ for *all* lags $t = 1 \dots n-1$ across all chains before checking for truncation. Even if the first pair is negative ($\hat\rho_1 + \hat\rho_2 < 0$), all $n(n-1)/2$ differences are computed. Phase 3 (step 3.4) of the preparation brief expands ESS to bulk/tail and constant handling, but does not address this quadratic cost or propose FFT-based autocovariances.
* **Acceptance Rate Burn-In Contamination ([`Samplers.md:169`](file:///d:/aghado01/capstone/runs/run_1/deliverables/docs/Samplers.md#L169)):**
  [`ChainEnsemble.Run`](file:///d:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainEnsemble.cs#L123) reports acceptance as $\sum\text{Accepted}/\sum\text{Attempts}$ from the chains' cumulative counters without subtracting burn-in steps.
* **Move Involutivity Not Validated in `ReversibleJumpChain` ([`Samplers.md:88`](file:///d:/aghado01/capstone/runs/run_1/deliverables/docs/Samplers.md#L88)):**
  [`ReversibleJumpChain`](file:///d:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/ReversibleJumpChain.cs#L62-L64) validates that each move's `ReverseKey` exists in the palette, but does not check reciprocal pairing (`moves[m.ReverseKey].ReverseKey == m.Key`).
* **`SpectralNoiseField` Boundary Periodicity ([`Samplers.md:172`](file:///d:/aghado01/capstone/runs/run_1/deliverables/docs/Samplers.md#L172)):**
  [`SpectralNoiseField.Generate`](file:///d:/aghado01/ThermoMapper/src/maths/samplers/SpectralNoiseField.cs) uses a periodic FFT grid, but [`Sample`](file:///d:/aghado01/ThermoMapper/src/maths/samplers/SpectralNoiseField.cs#L146) interpolates nodes $0 \dots N-1$ across the bounding box, breaking periodicity across the domain boundaries.

---

### Summary Checklist of Uncovered Issues

| Subsystem | File / Component | Uncovered Issue / Design Consideration | Impact / Severity |
|---|---|---|---|
| **Changepoint** | [`ExactChangepoint.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/changepoint/ExactChangepoint.cs) | Entire module omitted from both briefs; default $\sigma^2$ formula had documentation error. | Documentation / scope omission |
| **Changepoint** | [`SegmentModels.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/changepoint/SegmentModels.cs) | Cauchy table fill is $O(n^2 G)$ with $G=512$; truncates mass outside $\nu \pm 6\rho$. | Performance / unmonitored truncation |
| **BARS** | [`BarsEnsemble.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/BarsEnsemble.cs#L330-L333) | Span coverage $\pi(T) > 1.0$ when spans of nested significant peaks overlap. | Statistical semantic flaw / distorts `PeakModes` |
| **BARS** | [`WhittleModel.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WhittleModel.cs) & [`GibbsLossModel.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/GibbsLossModel.cs) | Omitted from R12 coefficient draws; return only conditional mean $\hat\beta$. | Uncertainty underestimation for nonlinear readouts |
| **BARS** | [`PoissonModel.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/PoissonModel.cs#L100-L143) | Hardcoded Laplace acceptance threshold ($-10$), 5 MH steps fallback to $\hat\beta$, $N(0, 100)$ ridge. | Heuristic inconsistency |
| **BAPS** | [`BapsSampler.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/BapsSampler.cs), [`HierarchicalPSpline.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/HierarchicalPSpline.cs) | Reject general [`IBandPenalty`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/IBandPenalty.cs); locked to [`DifferencePenalty`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/DifferencePenalty.cs). | Extension point breakage |
| **BAPS** | [`HierarchicalPSpline.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/HierarchicalPSpline.cs) | Runs single chain with zero convergence diagnostics ($\hat R$, ESS). | Reliability / diagnostic blind spot |
| **BAPS** | [`AnisotropicTensorSpline.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/AnisotropicTensorSpline.cs) | No optimizer provided for $(\lambda_x, \lambda_y)$. | Usability gap |
| **BGP** | [`GpRegression.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/bgp/GpRegression.cs) | Zero-mean GP fitted to uncentered response $y$; predictions pull to 0. | Statistical bias away from data |
| **Samplers** | [`ChainDiagnostics.cs`](file:///d:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainDiagnostics.cs#L150) | `Ess` eagerly computes all $n-1$ lags ($O(C n^2)$ cost). | Performance bottleneck on long chains |
| **Samplers** | [`ChainEnsemble.cs`](file:///d:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainEnsemble.cs#L123) | Reported acceptance rates conflate burn-in and post-burn draws. | Metric accuracy |
| **Samplers** | [`ReversibleJumpChain.cs`](file:///d:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/ReversibleJumpChain.cs#L62) | Does not validate reciprocal reverse key involution in move palettes. | Detailed balance vulnerability |