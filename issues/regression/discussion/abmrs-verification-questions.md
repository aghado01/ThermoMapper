Verify in the Biometrics text first. The paper is not yet in the compendium; the preprint is, as HYK2024.md. The external summary of the differences was wrong on one checkable point: it said the death formula changed, but the preprint already has the same formula. So check each item against the paper:

#	What to verify	Where it matters
1	Whether the prior keeps the preprint's form: π(k) ∝ C(C,k)^(1−γ) with π(ξ | k) = 1/C(C,k) (equation 3)	The γ_HYK = 1 + γ mapping, and so this step
2	γ's range (γ ≥ 0?), the three regimes and the recommended [0.5, 2]	Whether γ > 1 is their convention, to cite, or our extension, to label; constructor bounds; docs
3	Guidance on C ("typically > m"), or anything on k past C/2	Support [0, C] with guidance, or truncation
4	The b_k and d_k formulas; c = 0.3 and any stated reason or mixing claim	The mixing benchmark's arms (BARS mixing)
5	Coincident knots (multiplicity above 1)	The tensor-product extension's handling of jumps (port brief)
6	The non-normal treatment (still an EBIC plug-in?)	The docs' comparison with our Laplace marginals
7	How peak inference is defined and computed (draws or a plug-in)	The readout contract for the tensor-product extension
8	Whether the γ comparison still favors γ = 1, and whether it separates the prior from mixing	The confound noted under BARS mixing
9	Title, method name (ABMRS?), volume, issue, article ID, DOI	The citation
Answer: 

Junhui He, Ying Yang, Jian Kang, Adaptive Bayesian multivariate spline knot inference with prior specifications on model complexity, Biometrics, Volume 82, Issue 3, September 2026, ujag121, https://doi.org/10.1093/biomtc/ujag121
Method name: ABMRS (Adaptive Bayesian Multivariate Regression Splines)


--- 

# Report 

An investigation of the 9 verification questions against the official supplementary materials package ([`ABMRS_Revision_V2_Supplementary_Materials`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials)), including the Web Appendix PDF, C++ engine, R package, and demo scripts, yields the following findings:

---

### 1. Prior Form
* **Finding:** The prior in the published *Biometrics* version keeps the preprint's exact form:
  $$\pi(k) \propto \binom{C}{k}^{1-\gamma}, \quad \text{with } \pi(\xi \mid k) = \frac{1}{\binom{C}{k}}, \quad \implies \pi(k, \xi) \propto \binom{C}{k}^{-\gamma}$$
* **Evidence:**
  * In [Web Appendix A.2](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/web_appendix_extracted.txt#L51-L68) (Proof of Lemma 3):
    $$\pi(k, \xi) = \pi(\xi \mid k)\pi(k) \propto \tau(\mathcal{M}_k)^{-\gamma} = \binom{n}{k}^{-\gamma}$$
    For a birth step ($k' = k + 1$), the prior ratio is $\left(\frac{n-k}{k+1}\right)^{-\gamma}$.
  * The identity $\pi(k) b_k = \pi(k+1) d_{k+1}$ in Lemma 3 and the C++ code ([`ebars.cpp:10-18`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/src/ebars.cpp#L10-L18)) implement:
    $$\frac{b_k}{d_{k+1}} = \left(\frac{n-k}{k+1}\right)^{1-\gamma} \implies \frac{\pi(k+1)}{\pi(k)} = \left(\frac{n-k}{k+1}\right)^{1-\gamma}$$
  * *Note on earlier discussion summaries:* An external note claimed $\pi(k) \propto \binom{C}{k}^{-\gamma}$ and $\pi(\xi, k) \propto \binom{C}{k}^{-1-\gamma}$. That was algebraically mistaken; the author's text and code keep $\pi(k) \propto \binom{C}{k}^{1-\gamma}$, so $\gamma_{\text{HYK}} = 1 + \gamma_{\text{TM}}$.

---

### 2. $\gamma$'s Range and Regimes
* **Finding:** $\gamma > 1$ and the recommended range $[0.5, 2.0]$ are the authors' official convention and documentation, not an external extension.
* **Evidence:**
  * Package docstrings in [`wrapper.R:13-14, 140-141`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/R/wrapper.R#L13-L14):
    > `@param gamma double, the exponent in extended BIC, default value is 1.0. It should be constrained to the range [0.5, 2.0].`
  * The three regimes operate as:
    1. $\gamma < 1$ ($1 - \gamma > 0$): Prior $\binom{C}{k}^{1-\gamma}$ peaks at $k = C/2$, favoring higher model complexity (less regularization).
    2. $\gamma = 1$ ($1 - \gamma = 0$): Flat count prior ($\binom{C}{k}^0 = 1$); proposals become constant $b_k = c, d_k = c$.
    3. $\gamma > 1$ ($1 - \gamma < 0$): Penalizes knot count, favoring sparsity (heavy regularization).
  * Web Appendix B (p. 4, line 114) explicitly notes that $\gamma = 1.2$ over-regularizes on oscillatory curves at small sample sizes:
    > *"the proposed method with $\gamma = 1.2$ yields a high RMSE for Curve 3 when the sample size is $m = 200$, indicating that $\gamma > 1$ can lead to overly regularized models for highly oscillatory curves."*
  * Demos ([`knot_inference.Rmd:123`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/demo/knot_inference.Rmd#L123), [`triceps.Rmd:52`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/demo/triceps.Rmd#L52)) evaluate both $\gamma = 1.0$ and $\gamma = 1.5$.

---

### 3. Guidance on $C$ and Support past $C/2$
* **Finding:** Support in the engine is nominally $[0, C]$ with no truncation at $\lfloor C/2 \rfloor$. Guidance is $C \ge 2m$ or $3m$.
* **Evidence:**
  * In [`wrapper.R:17-18, 236-237`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/R/wrapper.R#L17-L18): `times = 3.0` by default ($n = m \times \text{times} = 3m$). Demos set `times = 2` ($n = 2m$).
  * In [`ebars.cpp:25-32, 117-126`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/src/ebars.cpp#L25-L32): The candidate set is an equally spaced grid of $n$ points. There is no truncation condition at $n/2$.
  * Because $C \ge 2m$ (e.g. $n = 1000$ for $m = 500$), $C/2 = 500$ is far outside the active posterior knot support ($k \le 10$). Moreover, [`utils.cpp:20-27`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/src/utils.cpp#L20-L27) uses `mpart = min(20, floor(0.1*m))` to discard basis columns with insufficient data points, naturally preventing the chain from climbing into dense knot regimes.

---

### 4. Move Schedules ($b_k, d_k$) and $c = 0.3$
* **Finding:** The formulas are:
  $$b_k = c \min\left(1, \left[\frac{n-k}{k+1}\right]^{1-\gamma}\right), \quad d_k = c \min\left(1, \left[\frac{k}{n-k+1}\right]^{1-\gamma}\right), \quad r_k = 1 - b_k - d_k$$
* **Evidence:**
  * Confirmed in [`ebars.cpp:10-23`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/src/ebars.cpp#L10-L23) and [`mebars.cpp:8-21`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/src/mebars.cpp#L8-L21).
  * The death formula is identical to the preprint.
  * Setting $c = 0.3$ guarantees $b_k + d_k \le 2c = 0.60$, enforcing a strict minimum relocation floor of $r_k \ge 0.40$ (40%) across all states.
  * No separate mixing theorem is proven for $c = 0.3$ in the Web Appendix; it is treated as an empirical tuning constant ensuring sufficient relocation mobility.

---

### 5. Coincident Knots (Knot Multiplicity)
* **Finding:** Coincident knots (multiplicity $> 1$ at the same coordinate) are **NOT supported** in the ABMRS sampler.
* **Evidence:**
  * In [`ebars.cpp:85-101, 117-148`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/src/ebars.cpp#L85-L101), all knots are sampled **without replacement** from the distinct grid points $j/(n+1)$.
  * The previous claim that ABMRS had an explicit kernel for "knot multiplicity layers" was a hallucination in the external notes.
  * In the numerical simulations with jumps (e.g., [`README.Rmd:43`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/README.Rmd#L43) and [`knot_inference.Rmd:40`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/demo/knot_inference.Rmd#L40)), true synthetic data was generated using coincident knots in `splines::bs(..., knots = c(0.2, 0.2, 0.5, 0.7))`, but ABMRS fits them by placing distinct knots arbitrarily close on the fine candidate grid.

---

### 6. Non-Normal Treatment
* **Finding:** The non-normal treatment remains an asymptotic EBIC plug-in.
* **Evidence:**
  * In [`ABMRS/README.md:24-25`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/ABMRS/README.md#L24-L25):
    > *"In the non-normal cases, we utilize the extended Bayesian information criterion (EBIC) to approximate the posterior density."*
  * In [Web Appendix A.2](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/web_appendix_extracted.txt#L78-L79):
    > *"For the EBIC-based approximation, we substitute $\hat{p}$ for the corresponding $p$."*
  * The C++ engine implements only Gaussian regression; non-normal extensions rely on the asymptotic EBIC deviance plug-in rather than Laplace marginal integration.

---

### 7. Peak Inference
* **Finding:** There is **no analytical gradient or Hessian root-finding machinery** in ABMRS.
* **Evidence:**
  * A search across the entire supplementary code and documentation turns up zero instances of gradients, Hessians, or critical-point solvers.
  * The package exposes posterior draws of knots (`knots()`), regression coefficients (`coefs()`), and curve/surface predictive samples (`predict()`).
  * In [`triceps.Rmd:81-96`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/demo/triceps.Rmd#L81-L96), knot inference is performed as a **two-stage plug-in**: the authors compute the posterior mean knot locations `node = rowMeans(nodes_mat)`, plug `node` into a truncated basis, and fit a frequentist `glm`.
  * For surfaces ([`surface_fitting.Rmd:106-114`](file:///d:/aghado01/ThermoMapper/issues/regression/ABMRS_Revision_V2_Supplementary_Materials/demo/surface_fitting.Rmd#L106-L114)), predictions and credible bands are computed pointwise across MCMC draws (`rowMeans(pred)` and `apply(pred, 1, quantile)`). Extrema are simply read off the predictive draws.

---

### 8. $\gamma$ Comparison and Confounding with Mixing
* **Finding:** The comparison does **not** favor $\gamma = 1$, and it does **not** separate the prior from mixing.
* **Evidence:**
  * In Web Table 1, for oscillatory functions at $m = 200$, $\gamma = 0.8$ outperforms $\gamma = 1.0$ (RMSE $0.680$ vs $0.751$), while $\gamma = 1.2$ severely over-regularizes (RMSE $1.533$). At larger sample sizes, results across $\gamma \in [0.8, 1.2]$ are virtually identical; the authors describe ABMRS as "highly robust to the choice of $\gamma$."
  * Because the move schedule sets $b_k = c \min(1, [(n-k)/(k+1)]^{1-\gamma})$, the prior parameter $\gamma$ is directly embedded in the proposal probabilities. Varying $\gamma$ alters the prior and the dimension-jumping proposal schedule simultaneously. The authors do not separate the two.

---

### 9. Bibliographic Details & Citation
* **Article Title:** *Adaptive Bayesian multivariate spline knot inference with prior specifications on model complexity*
* **Authors:** Junhui He, Ying Yang, and Jian Kang
* **Method Name:** ABMRS (Adaptive Bayesian Multivariate Regression Splines); R package `ABMRS` (classes `EBARS` and `MEBARS`)
* **Journal:** *Biometrics*, Volume 82, Issue 3 (September 2026), Article `ujag121`
* **DOI:** [`10.1093/biomtc/ujag121`](https://doi.org/10.1093/biomtc/ujag121)
* **Code Repositories:** [GitHub ABMRS](https://github.com/junhuihe2000/ABMRS) · [GitHub exampleABMRS](https://github.com/junhuihe2000/exampleABMRS)