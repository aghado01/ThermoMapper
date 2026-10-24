
# argon comments 


# so i've been thinking about stopping conditions for my application of BARS in SPC where ive been wanting to use a progressive bisecting strategy of measuring evenly space points along the susceptibility curve and using bayesian adaptive splines to estimate a joint posterior on an estimate of susceptibility until i have recovered the curve confidently. but its been tricky thinking about when I actually would be confident, and i just discovered some things that might be useful

The core insight in that note is strong: **BARS posterior variance alone cannot serve as a stopping condition** because nonparametric spline posteriors are blind to aliasing. If a sharp susceptibility spike sits strictly between two grid points $\beta_i$ and $\beta_{i+1}$, BARS will fit a flat spline across the gap, place zero knots there, and report high confidence. Using the internal physics of the Gibbs measure to certify that the grid has resolved the underlying smooth manifold before trusting the spline is the right architecture.

However, if you try to implement the note's prescriptions literally, two pieces will break in practice—and fixing them reveals a much cleaner, exact stopping criterion.

---

### 1. Where to Call BS on the Note (The Practical Traps)

1. **"Higher cumulants come for free" is a Monte Carlo mirage:**
In theory, $\beta$-derivatives of thermodynamic averages are cumulants, and the cubic spline error scales as $h^4 \Vert{}\chi^{(4)}\Vert{}_\infty$. In practice, $\chi(\beta)$ is already a 2nd cumulant of the order parameter (or mean squared Fortuin-Kasteleyn cluster size). Its 4th derivative $\chi^{(4)}(\beta)$ is a **6th-order cross-cumulant** between cluster sizes and bond energy (and even for heat capacity $C_V$, the 4th derivative requires the 6th energy cumulant $\kappa_6(E)$). In Swendsen-Wang runs near a pseudo-critical transition, sample estimators of 5th and 6th cumulants have catastrophic variance due to heavy tails and multimodal energy distributions. You would need $10^6\text{--}10^7$ sweeps per point to stabilize $\hat{\chi}^{(4)}$, destroying your sample economy.


2. **The Endpoint Blind Spot in Point-Wise $\Delta \beta \le c / \hat{\sigma}_E$:**
The note suggests stopping when every interval satisfies $\Delta \beta \le c / \hat{\sigma}_E$ using the sample variance $\text{Var}(E) = \sigma_E^2$ at the grid points. **This fails during early bisection stages.** Suppose a sharp cluster-splitting transition occurs at $\beta_c \in (\beta_i, \beta_{i+1})$. At $\beta_i$, the cluster is fused and $\sigma_E(\beta_i)$ is small; at $\beta_{i+1}$, the cluster is fully split and $\sigma_E(\beta_{i+1})$ is also small. Checking $\hat{\sigma}_E$ only at the endpoints $\beta_i$ and $\beta_{i+1}$ will make the interval look flat and "resolved," causing you to stop right over an undiscovered phase transition.


3. **Fisher Zeros ($\text{Im }\beta_1$) are Diagnostic, Not Operational:**
While the strip of analyticity $\vert{}\text{Im }\beta_1\vert{}$ explains *why* the spline converges, extracting $\text{Im }\beta_1$ from real-axis MCMC data requires analytic continuation (searching for complex zeros of the partition function $Z(\beta + i\tau)$). You cannot reliably estimate $\text{Im }\beta_1$ until *after* you have already resolved the peak on the real axis.



---

### 2. The Mathematical Fix: Exact Interval Certifiers from 1st & 2nd Moments

You can make the thermodynamic length idea immune to the endpoint blind spot and avoid noisy high-order cumulants by exploiting two exact identities of the Potts exponential family $P_\beta(x) = Z(\beta)^{-1} e^{-\beta E(x)}$.

#### A. The Fundamental Theorem of Calculus Trick (Zero Blind Spots)

Because the Fisher information is the exact negative derivative of the mean energy,


$$\sigma_E^2(\beta) = \text{Var}_\beta(E) = -\frac{\partial \langle E \rangle_\beta}{\partial \beta} \ge 0$$


the mean energy $\langle E \rangle_\beta$ is **strictly monotonic** in $\beta$. Integrating across any grid interval $[\beta_i, \beta_{i+1}]$ (with $\Delta \beta_i = \beta_{i+1} - \beta_i > 0$) yields:


$$\int_{\beta_i}^{\beta_{i+1}} \sigma_E^2(\beta) \, d\beta = \langle E \rangle_{\beta_i} - \langle E \rangle_{\beta_{i+1}} \equiv \vert{}\Delta \langle E \rangle_i\vert{}$$


Look at what this buys you:

* **No transition can hide between grid points:** Every SPC cluster split releases a finite latent bond energy $\vert{}\Delta \langle E \rangle_i\vert{}$. Even if a razor-thin spike in $\sigma_E^2(\beta)$ sits strictly inside $(\beta_i, \beta_{i+1})$ so that both endpoint variances $\hat{\sigma}_E^2(\beta_i)$ and $\hat{\sigma}_E^2(\beta_{i+1})$ are near zero, the **exact integral** of that spike is captured by the drop in first moments $\langle E \rangle_{\beta_i} - \langle E \rangle_{\beta_{i+1}}$.
* **Exact Jeffreys Divergence:** Multiplying by $\Delta \beta_i$ gives the exact symmetrized KL (Jeffreys) divergence between the Gibbs measures at the endpoints, which bounds the thermodynamic length $L_i = \int_{\beta_i}^{\beta_{i+1}} \sigma_E(\beta) d\beta$ via Cauchy-Schwarz:



$$L_i^2 \le \Delta \beta_i \int_{\beta_i}^{\beta_{i+1}} \sigma_E^2(\beta) \, d\beta = \Delta \beta_i \left(\langle E \rangle_{\beta_i} - \langle E \rangle_{\beta_{i+1}}\right) = D_J(P_{\beta_i} \Vert{} P_{\beta_{i+1}})$$


* **The Quadrature Discrepancy Detector:** Since your SW runs give you both $\langle E \rangle$ and $\sigma_E^2$ at the endpoints, compare the exact integral $\vert{}\Delta \langle E \rangle_i\vert{}$ against the trapezoidal rule estimate from the endpoints:



$$\mathcal{R}_i = \frac{\frac{\Delta \beta_i}{2}\left(\hat{\sigma}_E^2(\beta_i) + \hat{\sigma}_E^2(\beta_{i+1})\right)}{\langle E \rangle_{\beta_i} - \langle E \rangle_{\beta_{i+1}}}$$



If a peak is hiding inside $(\beta_i, \beta_{i+1})$, the numerator is small while the denominator is large ($\mathcal{R}_i \ll 1$). If the interval is sitting on the shoulder of a peak, $\mathcal{R}_i \neq 1$ measures the integrated curvature $\int (\sigma_E^2)'' d\beta$ without ever computing a noisy 4th sample moment. When $\mathcal{R}_i \approx 1$ and $D_J \le c^2$, the interval is provably resolved.

#### B. The Direct Bridge from $\sigma_E$ to Susceptibility $\chi(\beta)$

The note treats thermodynamic length ($\sigma_E d\beta$) and the susceptibility curve $\chi(\beta)$ as separate concepts, but they are directly coupled. Let $\hat{\chi}(x)$ be your single-sweep Fortuin-Kasteleyn susceptibility observable, so $\chi(\beta) = \langle \hat{\chi} \rangle_\beta$. Its exact $\beta$-derivative is the cross-covariance with energy:


$$\frac{\partial \chi(\beta)}{\partial \beta} = -\text{Cov}_\beta(\hat{\chi}, E)$$


By Cauchy-Schwarz, the rate of change of the true susceptibility curve is bounded by the thermodynamic length metric:


$$\left\vert{} \frac{d\chi(\beta)}{d\beta} \right\vert{} \le \sigma_{\hat{\chi}}(\beta) \, \sigma_E(\beta) \implies \frac{\vert{}d\chi(\beta)\vert{}}{\sigma_{\hat{\chi}}(\beta)} \le \sigma_E(\beta) \, d\beta = dL_{\text{thermo}}$$


Enforcing $\Delta L_i \le c$ guarantees that **the true susceptibility $\chi(\beta)$ cannot change across the interval $[\beta_i, \beta_{i+1}]$ by more than $c$ single-sweep standard deviations $\sigma_{\hat{\chi}}$**.

Even better: $-\widehat{\text{Cov}}(\hat{\chi}, E)$ is only a **2nd moment**, so it is statistically stable and truly comes for free from your SW sweeps. It gives you a direct empirical measurement of the slope $\chi'(\beta_i)$ at every grid point!

---

### 3. A Concrete 3-Gate Stopping Rule for Progressive Bisection + BARS

When running progressive bisection (or dyadic adaptive refinement, which you should strongly consider over global uniform bisection so you don't waste sweeps bisecting flat plateaus), stop bisecting an interval $[\beta_i, \beta_{i+1}]$ at stage $m$ when it passes three gates:

1. **Gate 1 — Thermodynamic Resolution (No Hidden Features):**
The exact Jeffreys divergence across the interval satisfies:

$$D_{J,i} = \Delta \beta_i \left(\hat{E}(\beta_i) - \hat{E}(\beta_{i+1})\right) \le c_J^2 \quad (\text{typically } c_J \in [1, 2])$$

AND, for any interval with non-trivial energy drop $\vert{}\Delta \hat{E}_i\vert{}$ above noise, the quadrature ratio satisfies $\vert{}\mathcal{R}_i - 1\vert{} \le \tau_{\text{quad}}$ (e.g., $\tau_{\text{quad}} \approx 0.15$), certifying that $\sigma_E^2(\beta)$ is approximately linear across the sub-interval and no narrow transition spike sits inside $(\beta_i, \beta_{i+1})$.
2. **Gate 2 — Tangent / Derivative Consistency (Physics vs. Spline):**
At the grid points, the BARS posterior spline derivative $s'(\beta_i)$ agrees within posterior + MCMC credible bounds with the empirical cross-covariance slope $-\widehat{\text{Cov}}_{\beta_i}(\hat{\chi}, E)$. (If BARS is smoothing over a sharp corner, the empirical tangent $-\widehat{\text{Cov}}(\hat{\chi}, E)$ will systematically exceed the spline derivative).
3. **Gate 3 — Cross-Stage BARS Posterior Stability:**
Between bisection stage $m-1$ and stage $m$, the functional change in the BARS posterior mean $\hat{\chi}_{\text{BARS}}^{(m)}(\beta)$ (or the Supremum / $L_2$ norm scaled by the posterior standard deviation band) is below tolerance:

$$\sup_{\beta} \frac{\left\vert{} \mathbb{E}[\chi(\beta) \mid \mathcal{D}_m] - \mathbb{E}[\chi(\beta) \mid \mathcal{D}_{m-1}] \right\vert{}}{\text{SD}[\chi(\beta) \mid \mathcal{D}_m]} \le \tau_{\text{BARS}}$$

*(Bonus lookahead)*: Once Gate 1 holds ($D_{J,i} \lesssim 2$), neighbouring energy histograms overlap. Before even running SW sweeps at the next bisection midpoint $\beta_{i+1/2}$, you can use **Ferrenberg-Swendsen multi-histogram reweighting** from $\beta_i$ and $\beta_{i+1}$ to predict $\hat{\chi}_{\text{FS}}(\beta_{i+1/2})$ for zero extra MCMC cost and check whether it already sits inside your BARS credible band.

### 4. What to Actually Target in Lean 4

Instead of formalizing complex-plane Fisher zero Cauchy bounds—where the constant $\vert{}\text{Im }\beta_1\vert{}$ is unknown on arbitrary data graphs anyway—the load-bearing theorem for your harness is the **Exponential Family Interval Certifier**:

* For any finite state space $\Omega$, energy $E: \Omega \to \mathbb{R}$, and observable $A: \Omega \to \mathbb{R}$ under the Gibbs measure $P_\beta(\omega) = Z(\beta)^{-1} e^{-\beta E(\omega)}$:
1. $\partial_\beta \langle E \rangle_\beta = -\text{Var}_\beta(E)$ and $\int_{\beta_1}^{\beta_2} \text{Var}_\beta(E) d\beta = \langle E \rangle_{\beta_1} - \langle E \rangle_{\beta_2}$.
2. $\vert{}\partial_\beta \langle A \rangle_\beta\vert{} = \vert{}\text{Cov}_\beta(A, E)\vert{} \le \sqrt{\text{Var}_\beta(A)\text{Var}_\beta(E)}$.
3. Combining (1) and (2) via Cauchy-Schwarz yields the exact excursion bound across any interval $[\beta_1, \beta_2]$:

$$\left\vert{} \langle A \rangle_{\beta_2} - \langle A \rangle_{\beta_1} \right\vert{} \le \sqrt{\int_{\beta_1}^{\beta_2} \text{Var}_\beta(A) \, d\beta} \cdot \sqrt{\langle E \rangle_{\beta_1} - \langle E \rangle_{\beta_2}}$$

That last inequality is completely free of unknown lattice scaling constants, uses only finite sums and standard calculus in Mathlib, and bounds the maximum unobserved excursion of $\chi(\beta)$ directly in terms of the measurable energy drop $\Delta \langle E \rangle$.

Are you committed to bisecting the entire temperature grid uniformly at each volley, or does your worker-volley orchestrator allow **local dyadic bisection** (only splitting intervals where $D_{J,i}$ or the BARS posterior variance exceeds tolerance)?



# i start with some number of points on the interval [0,1], always including 0 and 1 initially. then i take the temperatures in between each of the original points. then i do it again. after each iteration, i fit bars on the entire set of points up to that time, so there isn't dependency from one iteration to the next. i dont like adaptive strategy because i dont know what the system's shape looks like, i would rather fit on an even grid and see the shape reveal itself as the density increases

Refitting BARS from scratch on a globally bisected uniform grid on $T \in [0, 1]$ is actually a major statistical advantage for stopping—not a limitation.

Adaptive schemes introduce sequential experimental design bias and can easily tunnel-vision on the first broad peak while ignoring a narrow, decoupled hierarchical split in a "flat" region. By doubling a uniform grid $\mathcal{G}_{m-1} \to \mathcal{G}_m$ (with spacing $h_m = h_0 2^{-m}$) and refitting BARS independently at each stage $m$, you get **two structural invariants** that an adaptive scheme never has:

1. **Every new volley of midpoints $\mathcal{M}_m = \mathcal{G}_m \setminus \mathcal{G}_{m-1}$ is an unbiased, space-filling out-of-sample holdout set** for your Stage $m-1$ BARS posterior.
2. **Every interval $[T_j, T_{j+1}]$ from Stage $m-1$ is split into two exact halves $[T_j, T_{j+1/2}]$ and $[T_{j+1/2}, T_{j+1}]$**, allowing you to check dyadic conservation of energy and heat capacity across scales.

To turn your note's Fisher-zero and thermodynamic-length concepts into an exact stopping rule for this global bisection loop, you need to combine three checks at the end of Stage $m$: **(1) Out-of-Sample Midpoint Surprise**, **(2) $T$-Space Thermodynamic Resolution (so nothing hides at Stage $m+1$)**, and **(3) Knot Posterior Saturation**.

---

### 1. The Free Holdout Test: Midpoint Predictive Surprise (Stage $m-1 \to m$)

The note mentions using conformal Bayes to check calibration. Because your grid at Stage $m$ consists of the Stage $m-1$ grid $\mathcal{G}_{m-1}$ plus the exact midpoints $\mathcal{M}_m = \{T_{j+1/2}\}$, **you get an exact out-of-sample calibration check at every bisection step before you even refit BARS on $\mathcal{G}_m$.**

Let $\hat{\chi}_{m-1}(T)$ and $s_{m-1}^2(T)$ be the BARS posterior mean and predictive variance (spline posterior variance + MCMC observation noise variance $\hat{\sigma}_{\bar{\chi}}^2$) fitted *only* on $\mathcal{G}_{m-1}$. When Stage $m$ finishes running SW sweeps at the midpoints $T_{j+1/2} \in \mathcal{M}_m$, compute the **standardized midpoint innovations**:


$$z_{j+1/2}^{(m)} = \frac{\bar{\chi}_{\text{SW}}(T_{j+1/2}) - \hat{\chi}_{m-1}(T_{j+1/2})}{s_{m-1}(T_{j+1/2})}$$

* **While under-resolved ($h_m > \vert{}\text{Im } T_1\vert{}$):** Whenever the new midpoints hit an previously unresolved peak or steep shoulder, $\vert{}z_{j+1/2}^{(m)}\vert{} \gg 3$. The Stage $m-1$ spline bridged across the feature, and the midpoint measurement exposes the aliasing bias.


* **When resolved ($h_m \ll \vert{}\text{Im } T_1\vert{}$):** The spline error $O(h_{m-1}^4 \Vert{}\chi^{(4)}\Vert{}_\infty)$ drops below the Monte Carlo noise floor. The innovations $\{z_{j+1/2}^{(m)}\}$ across all $N_{m-1}-1$ midpoints behave like standard normal noise ($\max_j \vert{}z_{j+1/2}^{(m)}\vert{} \le \Phi^{-1}(1 - \alpha / 2\vert{}\mathcal{M}_m\vert{})$ and empirical $(1-\alpha)$ credible interval coverage matches nominal).



---

### 2. Certifying That Nothing Hides at Stage $m+1$ (In $T \in [0, 1]$ Coordinates)

Passingly passing the midpoint test between Stage $m-1$ and Stage $m$ proves that Stage $m$ didn't uncover a surprise—but how do you know a super-narrow transition isn't still hiding at the *quarter-points* (which wouldn't be sampled until Stage $m+1$)?

This is where the physics note comes in, provided we translate its $\beta$-space equations into your uniform $T \in [0, 1]$ grid to handle the Jacobian $dT = -T^2 d\beta$ and the $T=0$ boundary without singularities:

* **Exact $T$-Space Energy Identity:** Since $\beta = 1/T$, the chain rule gives:

$$\frac{\partial \langle E \rangle_T}{\partial T} = \frac{\partial \langle E \rangle}{\partial \beta} \frac{d\beta}{dT} = (-\sigma_E^2)\left(-\frac{1}{T^2}\right) = \frac{\sigma_E^2(T)}{T^2} = C_V(T) \ge 0$$



Integrating across any uniform grid interval $[T_j, T_{j+1}]$ of width $h_m = T_{j+1} - T_j$ gives:

$$\int_{T_j}^{T_{j+1}} C_V(T) \, dT = \langle E \rangle_{T_{j+1}} - \langle E \rangle_{T_j} \equiv \Delta \langle E \rangle_j$$


* **Exact $T$-Space Thermodynamic Length Bound:** The Fisher information with respect to $T$ is $\mathcal{I}_T(T) = \sigma_E^2(T) / T^4$, so the thermodynamic length element $ds = \sigma_E \vert{}d\beta\vert{} = \frac{\sigma_E(T)}{T^2} dT$ across $[T_j, T_{j+1}]$ ($T_j > 0$) satisfies Cauchy-Schwarz:



$$L_j = \int_{T_j}^{T_{j+1}} \left(\frac{\sigma_E(T)}{T}\right)\left(\frac{1}{T}\right) dT \le \sqrt{\int_{T_j}^{T_{j+1}} \frac{\sigma_E^2(T)}{T^2} dT \int_{T_j}^{T_{j+1}} \frac{dT}{T^2}} = \sqrt{\frac{h_m}{T_j T_{j+1}} \, \Delta \langle E \rangle_j}$$



Notice how $\frac{h_m}{T_j T_{j+1}} \Delta \langle E \rangle_j$ is the **exact Jeffreys divergence** $D_{J,j}$ across your uniform $T$-interval, automatically accounting for the $T^2$ compression near $T \to 0$! (For the boundary interval $[0, h_m]$, once $h_m$ is below the minimum bond-breaking gap so $\langle E \rangle_{h_m} = E_{\text{ground}}$, $\Delta \langle E \rangle_0 = 0$ and the interval is trivially resolved).



Because every Stage $m-1$ interval $[T_j, T_{j+1}]$ was just bisected at $T_{j+1/2}$, you can check two exact physical invariants at Stage $m$:

1. **Uniform Thermodynamic Overlap:** Every active sub-interval on $\mathcal{G}_m$ satisfies:

$$D_{J} = \frac{h_m}{T_k T_{k+1}} \left(\langle \hat{E} \rangle_{T_{k+1}} - \langle \hat{E} \rangle_{T_k}\right) \le c_J^2 \quad (\text{with } c_J \approx 1.5\text{--}2)$$



which guarantees overlapping energy distributions everywhere on the grid.


2. **Dyadic Quadrature Closure (No Hidden Spike):** On every bisected parent interval $[T_j, T_{j+1}]$, compare the exact energy drop $\Delta \langle \hat{E} \rangle_j = \langle \hat{E} \rangle_{T_{j+1}} - \langle \hat{E} \rangle_{T_j}$ against Simpson's / trapezoidal quadrature of $\hat{C}_V(T) = \hat{\sigma}_E^2(T)/T^2$ using the endpoints and the new midpoint $T_{j+1/2}$:

$$\left\vert{} \Delta \langle \hat{E} \rangle_j - \frac{h_m}{2}\left(\hat{C}_V(T_j) + 2\hat{C}_V(T_{j+1/2}) + \hat{C}_V(T_{j+1})\right) \right\vert{} \le \tau_{\text{quad}} \, \text{SE}(\Delta \langle \hat{E} \rangle_j)$$



If a narrow phase transition is hiding inside $[T_j, T_{j+1/2}]$ or $[T_{j+1/2}, T_{j+1}]$, the energy drop $\Delta \langle \hat{E} \rangle_j$ will be significantly larger than the heat-capacity quadrature at the grid points.

---

### 3. BARS Structural Convergence: Knot Saturation ($p(k \mid \mathcal{D}_m)$)

Because you refit BARS from scratch on the entire set of $N_m$ points at each stage, the posterior distribution over the **number of knots $k$** is a direct indicator of whether the curve's geometry has been fully discovered:

* **Discovery Regime (Early Stages):** As $h_m$ halves, new peaks and inflection points become visible. To fit them, the BARS reversible-jump MCMC must add knots, so the posterior mean knot count $\mathbb{E}[k \mid \mathcal{D}_m]$ increases with $m$.
* **Asymptotic / Resolved Regime:** Once $h_m$ is smaller than the narrowest Fisher zero strip width $\vert{}\text{Im } T_1\vert{}$, doubling the grid density adds zero new geometric features to $\chi(T)$. Because BARS places a complexity penalty on $k$ (via the marginal likelihood / BIC-like phase space volume of each knot), **doubling the number of data points $N_m = 2 N_{m-1} - 1$ stops increasing the number of knots $k$**. Instead, $\mathbb{E}[k \mid \mathcal{D}_m] \approx \mathbb{E}[k \mid \mathcal{D}_{m-1}]$, and the new midpoints only shrink the credible band width around the existing spline by $\sim 1/\sqrt{2}$.



---

### Summary of the Stopping Protocol at Bisection Stage $m$

After collecting the midpoint SW runs $\mathcal{M}_m$ and refitting BARS on the full uniform grid $\mathcal{G}_m$, stop bisecting when all four conditions hold simultaneously:

| Gate | Metric Computed at Stage $m$ | What It Certifies |
| --- | --- | --- |
| **1. Thermodynamic Density** | $\max_{k} \frac{h_m}{T_k T_{k+1}} (\langle \hat{E} \rangle_{T_{k+1}} - \langle \hat{E} \rangle_{T_k}) \le c_J^2$ | Adjacent temperatures have overlapping Boltzmann energy histograms across the worst transition.

 |
| **2. Quadrature Closure** | $\frac{h_m}{2}(\hat{C}_V(T_k) + \hat{C}_V(T_{k+1})) \approx \langle \hat{E} \rangle_{T_{k+1}} - \langle \hat{E} \rangle_{T_k}$ | No unobserved transition spike is hiding strictly between grid points $T_k$ and $T_{k+1}$. |
| **3. Midpoint Holdout** | $\max_{T \in \mathcal{M}_m} \frac{\vert{}\bar{\chi}_{\text{SW}}(T) - \hat{\chi}_{m-1}(T)\vert{}}{s_{m-1}(T)} \le z_{1-\alpha}$ | The new midpoints $\mathcal{M}_m$ produced zero out-of-sample surprises relative to Stage $m-1$'s BARS fit.

 |
| **4. Knot & Curve Stability** | $\vert{}\mathbb{E}[k \mid \mathcal{D}_m] - \mathbb{E}[k \mid \mathcal{D}_{m-1}]\vert{} < 1$ and $\sup_T \frac{\vert{}\hat{\chi}_m(T) - \hat{\chi}_{m-1}(T)\vert{}}{\text{SD}[\chi(T) \mid \mathcal{D}_m]} \le \tau$ | BARS has stopped discovering new geometric features and the posterior curve has stabilized. |

How many points are you starting with at Stage $0$ (e.g., $N_0 = 9$ or $11$), and are you fitting BARS only to susceptibility $\chi(T)$ right now, or jointly across your other lenses (like cyclomatic complexity and compressibility) on the same grid?



# initially it would only be susceptibility because that's well understood for the signal. from what i've seen the other observables add complexity and interpretation burden and the relationship to phase behavior isn't always clear 

Dropping the auxiliary observables (like compressibility and cyclomatic complexity) and focusing strictly on susceptibility $\chi(T)$ is the right call—not just for engineering simplicity, but mathematically. Those topological and information-theoretic heuristics are not canonical expectations of the Gibbs measure $P_\beta$, which means the entire analytic framework in your note (Fisher zeros of $Z(\beta)$, cumulant-derivative identities, and thermodynamic length) does not even apply to them.

For $\chi(T)$, the physics is clean and exact. Focusing solely on recovering $\chi(T)$ via uniform bisection + BARS clarifies the exact problem your stopping condition has to solve: **$\chi(T)$ is non-monotonic.**

---

### 1. The Core Vulnerability of a $\chi(T)$-Only Grid: Peak Aliasing vs. Finite-Size Scaling

Because $\chi(T)$ rises at a cluster-splitting transition $T_c$ and falls back down on both sides, a sharp transition sitting strictly inside an interval $(T_k, T_{k+1})$ leaves both endpoint measurements $\bar{\chi}(T_k)$ and $\bar{\chi}(T_{k+1})$ near baseline. If you only look at $\chi(T)$ at the grid points, an interval containing a razor-sharp spike looks indistinguishable from a flat plateau.

This is where the two ideas in your note—Fisher zeros ($\text{Im }\beta_1$) and thermodynamic length ($\sigma_E d\beta$)—lock together via a **Finite-Size Scaling Duality** specific to hierarchical SPC:

* Suppose a subcluster of $n_s$ nodes splits at some critical temperature $T_c$. By finite-size scaling, the local partition function zeros associated with that subcluster pinch toward the real axis at distance $\vert{}\text{Im }\beta_1\vert{} \sim n_s^{-1/\nu}$ (or $\sim n_s^{-1}$ for first-order Potts splits).


* **Small Splits ($n_s$ small) $\implies$ Wide Strip, Smooth $\chi(T)$:** When a small peripheral group of nodes detaches, $n_s$ is small, so $\vert{}\text{Im }\beta_1\vert{}$ is **large**. The resulting feature in $\chi(T)$ is broad and smooth, meaning a coarse bisection grid already hits it and BARS resolves it early.


* **Macro Splits ($n_s$ large) $\implies$ Narrow Strip, Spike in $\chi(T)$:** The transitions that are dangerously narrow in $T$ (tiny $\vert{}\text{Im }\beta_1\vert{}$, where a uniform grid can step completely over the peak) are the **large cluster splits** ($n_s \gg 1$).


* **The Duality:** While a large-$n_s$ split produces a *narrow, non-monotonic spike* in $\chi(T)$ that can hide between $T_k$ and $T_{k+1}$, it releases a *massive, monotonic step* $\Delta \langle E \rangle \propto n_s$ in the mean bond energy $\langle E \rangle_T$.

In short: **the transitions narrow enough to hide from $\chi(T)$ are loud in $\Delta \langle E \rangle_T$, and the transitions quiet in $\Delta \langle E \rangle_T$ are wide enough to be seen directly in $\chi(T)$.**

---

### 2. Two Stopping Architectures for Recovering $\chi(T)$

Depending on whether you want your SW workers to return *only* susceptibility statistics or include one passive scalar (mean bond count/energy) alongside $\chi(T)$, here is how the stopping criteria compare:

| Feature | Option A: Pure $\chi(T)$ Only (Zero Extra Observables) | Option B: $\chi(T)$ BARS + Passive $\Delta \bar{E}$ Guardrail |
| --- | --- | --- |
| **Worker Payload per $T_k$** | $(\bar{\chi}_k, \, s_{\hat{\chi}, k}^2)$ | $(\bar{\chi}_k, \, s_{\hat{\chi}, k}^2, \, \bar{E}_k)$ |
| **What BARS Fits** | 1D spline on $\bar{\chi}(T)$ weighted by $s_{\hat{\chi}, k}^2 / M_{\text{eff}}$ | 1D spline on $\bar{\chi}(T)$ weighted by $s_{\hat{\chi}, k}^2 / M_{\text{eff}}$ |
| **Noise / Calibration Check** | Out-of-sample midpoint $z$-scores on $\mathcal{M}_m$ + BARS knot saturation $\mathbb{E}[k \mid \mathcal{D}_m]$ | Out-of-sample midpoint $z$-scores on $\mathcal{M}_m$ + BARS knot saturation $\mathbb{E}[k \mid \mathcal{D}_m]$ |
| **Anti-Aliasing Guarantee** | **Empirical (Two-Stage Quiet Rule):** Must wait for two consecutive bisections ($m-1, m$) with zero midpoint surprises to reduce the chance of quarter-point aliasing. | **Exact (Single-Stage Certificate):** $\frac{h_m}{T_k T_{k+1}}(\bar{E}_{k+1} - \bar{E}_k) \le c_J^2$ guarantees no macro-split is hiding inside $(T_k, T_{k+1})$.

 |

#### How Option A Works (Strictly $\chi(T)$ Only)

If you want the pipeline to touch nothing but $\chi(T)$, you already collect two numbers at each temperature $T_k$: the sample mean $\bar{\chi}(T_k)$ and the across-sweep variance $s_{\hat{\chi}}^2(T_k)$ (needed for the heteroskedastic noise variance in BARS).

Even without energy, $s_{\hat{\chi}}^2(T_k) = \widehat{\text{Var}}_{T_k}(\hat{\chi})$ is a **4th-order order-parameter cumulant** that spikes even more sharply near $T_c$ than $\bar{\chi}(T_k)$ itself. In a pure-$\chi$ setup, you stop at Stage $m$ when:

1. **Midpoint Predictive Coverage:** Every new midpoint $T_{j+1/2} \in \mathcal{M}_m$ falls inside the Stage $m-1$ BARS predictive band:

$$\max_{j} \frac{\left\vert{}\bar{\chi}(T_{j+1/2}) - \hat{\chi}_{m-1}(T_{j+1/2})\right\vert{}}{\sqrt{\text{Var}_{\text{BARS}}[\chi(T_{j+1/2}) \mid \mathcal{D}_{m-1}] + s_{\hat{\chi}}^2(T_{j+1/2})/M_{\text{eff}}}} \le z_{1 - \alpha}$$


2. **Variance-Peak Resolution:** No single grid point $T_k$ carries an isolated spike in sweep variance $s_{\hat{\chi}}^2(T_k)$ without at least one neighboring grid point $T_{k\pm 1}$ also showing elevated variance (confirming the peak width is $> 2 h_m$).
3. **Knot & Functional Stabilization:** The BARS posterior knot count stops growing ($\left\vert{}\mathbb{E}[k \mid \mathcal{D}_m] - \mathbb{E}[k \mid \mathcal{D}_{m-1}]\right\vert{} < 1$) and the $L_\infty$ shift in the posterior mean curve is within a fraction of the posterior SD band.

#### Why Option B's Single Extra Scalar ($\bar{E}_k$) Is Worth Considering

You don't fit a spline to $\bar{E}_k$ or interpret it as a clustering signal. In Swendsen-Wang, when you iterate over graph edges to activate Fortuin-Kasteleyn bonds with probability $p_{ij} = 1 - e^{-J_{ij}/T}$, the sum of unsatisfied edge weights $E(s) = \sum_{\langle i,j \rangle} J_{ij}(1 - \delta_{s_i, s_j})$ is already sitting in register memory.

Accumulating its scalar mean $\bar{E}_k$ alongside $\bar{\chi}_k$ gives you the exact integral of the Fisher information across every interval, $\int_{T_k}^{T_{k+1}} \frac{\sigma_E^2(T)}{T^2} dT = \bar{E}_{k+1} - \bar{E}_k$. That turns the thermodynamic length rule from `opus-fisher-etc.md` into a single division-free inequality check on adjacent grid points before you trust Option A's BARS convergence.

In your current SW kernel, are you computing $\hat{\chi}$ from the Fortuin-Kasteleyn cluster size distribution ($\sum c_k^2$) or from the spin magnetization variance ($N_{\max}$) at each sweep? (That choice slightly changes how the within-temperature variance $s_{\hat{\chi}}^2(T_k)$ behaves on the high-$T$ side of a split.)


Historically, Blatt, Wiseman, and Domany (1996/1997) actually used **spin magnetization** $m = \frac{(N_{\max}/N)q - 1}{q - 1}$ (where $N_{\max} = \max_\sigma N_\sigma$ is the largest spin-color group) and defined $\chi(T) = \frac{N}{T}(\langle m^2 \rangle - \langle m \rangle^2)$, even though they used Swendsen-Wang (which builds Fortuin-Kasteleyn clusters) to sample the states.

Switching from spin magnetization to the **FK cluster representation** is a strict mathematical improvement—it is an exact Rao-Blackwellization of the Potts model via the Edwards-Sokal coupling. Furthermore, your memory that magnetization looks "smoother" is physically real, and understanding *why* it looks smoother explains why FK is much better suited for your BARS bisection pipeline.

---

### 1. Why Magnetization Looks Smoother (And Why It's an Artifact)

In each Swendsen-Wang sweep, you first form connected FK components of sizes $\{c_1, c_2, \dots, c_K\}$ using bond probabilities $p_{ij} = 1 - e^{-J_{ij}/T}$, and then you assign each component one of $q$ colors uniformly at random. Moving from $\{c_k\}$ to spin magnetization $N_{\max}$ introduces two effects that artificially smooth the curve:

1. **$1/q$ Accidental Color Collisions:**
Suppose a cluster splits into two disconnected FK components $c_1$ and $c_2$ just above $T_c$. In the FK representation, they are strictly separate ($0$ active bonds between them). In the spin representation, on **$1/q$ of all sweeps**, $c_1$ and $c_2$ are randomly assigned the *exact same spin color*, causing the magnetization estimator to jump back up to $c_1 + c_2$. This stochastic mixing between the fused and split states blurs the transition over a wider temperature window and injects Multinomial($1/q$) coloring variance into every sweep.
2. **Giant-Component Masking of Hierarchical Splits:**
As your note points out, hierarchical SPC transitions create multiple near-real Fisher zeros corresponding to different subcluster splits. Suppose your graph has a large dense core $A$ ($300$ nodes) that stays fused until high $T$, and a smaller peripheral cluster $B$ ($80$ nodes) that splits into $40 + 40$ at a lower temperature $T_{c,B}$.


* **$N_{\max}$ (Magnetization) is nearly blind to $B$ splitting:** Because $N_{\max}$ only tracks the *single largest* spin state ($A$'s color), the $80 \to 40+40$ split of $B$ is invisible on $(q-1)/q$ of the sweeps, only appearing as tiny $1/q$ collision noise when $B$'s pieces happen to draw $A$'s color. That makes the $\chi_{\text{mag}}(T)$ curve look artificially smooth across $T_{c,B}$ by erasing secondary hierarchical peaks.
* **FK $\sum_k c_k^2$ sees every split equally:** Whenever *any* cluster splits into $c_1 + c_2$ anywhere in the graph, $\sum_k c_k^2$ drops by exactly:

$$(c_1 + c_2)^2 - (c_1^2 + c_2^2) = 2 c_1 c_2$$



completely independent of whether a larger cluster $A$ exists elsewhere.



---

### 2. How Choice of FK Estimator Affects BARS and Your Stopping Rule

In the Edwards-Sokal coupling, the marginal partition function of the FK bond measure $P_{\text{FK}}(\omega) \propto q^{C(\omega)} \prod_{e \in \omega} (e^{\beta J_e} - 1)$ is the **exact same partition function $Z(\beta)$** as the spin Potts model. That means the Fisher zero strip width $\vert{}\text{Im }\beta_1\vert{}$ and the thermodynamic length $ds = \sigma_E d\beta$ from your note apply identically to FK observables.

However, depending on which of the three standard FK definitions of $\chi(T)$ you are computing, the curve geometry and noise model for BARS behave very differently:

| FK Estimator Definition | Curve Shape on $T \in [0, 1]$ | Captures Shadowed Hierarchical Splits? | MCMC Noise Structure for BARS |
| --- | --- | --- | --- |
| **1. Unsubtracted FK Mean**<br>

<br>$\chi_1(T) = \left\langle \frac{1}{N} \sum_k c_k^2 \right\rangle_T$ | **Monotonic step function** from $N$ (at $T=0$) down to $1$ (at $T=1$). Peaks are in $-d\chi_1 / dT$. | **Yes** (every split drops the curve by $2 c_1 c_2 / N$). | **Sample Mean** (Gaussian CLT errors, lowest variance, exact Rao-Blackwellized $\sum_{i,j} \langle \delta_{s_i, s_j} \rangle$). |
| **2. Percolation Subtracted Mean**<br>

<br>$\chi_2(T) = \left\langle \frac{1}{N} \left(\sum_k c_k^2 - c_{\max}^2\right) \right\rangle_T$ | **Peaked curve** ($0$ at $T=0$, peaks at transitions, $\to 1$ at $T=1$). | **Partially** (splits of the largest cluster peak sharply; splits of smaller clusters drop $\sum_{k \neq \max} c_k^2$). | **Sample Mean**, but has a sharp derivative cusp where the identity of $c_{\max}$ swaps. |
| **3. FK Susceptibility Variance**<br>

<br>$\chi_3(T) = \frac{1}{NT} \text{Var}_T\left(\frac{1}{N}\sum_k c_k^2\right)$<br>

<br>*(or $\frac{1}{NT}\text{Var}_T(c_{\max})$)* | **Multi-peaked curve** ($0$ on plateaus, positive spikes at every transition $T_c$). | **Yes** for $\text{Var}(\sum c_k^2)$ (spikes whenever *any* subcluster fluctuates between fused and split). | **Sample Variance** (Chi-square / Gamma error tails; requires $\log$-transform or bootstrap SEs for BARS). |

### 3. The Key Implication for Your Stopping Condition

If you use **Estimator 1** ($\langle \frac{1}{N}\sum_k c_k^2 \rangle_T$) or **Estimator 3** ($\text{Var}_T(\frac{1}{N}\sum_k c_k^2)$) together:

* Notice that $\frac{1}{N}\sum_k c_k^2$ is **monotonic in expectation** across transitions (just like mean energy $\langle E \rangle_T$), whereas its across-sweep variance $\text{Var}_T(\frac{1}{N}\sum_k c_k^2)$ forms the sharp transition spikes!
* That means **even if you don't track energy $E$ at all**, the drop in the unsubtracted FK mean $\Delta \left\langle \frac{1}{N}\sum_k c_k^2 \right\rangle$ across a grid interval $[T_j, T_{j+1}]$ acts as an **anti-aliasing certifier** for the peaks in $\text{Var}_T(\frac{1}{N}\sum_k c_k^2)$! No cluster split of size $2 c_1 c_2$ can hide its variance spike between $T_j$ and $T_{j+1}$ without producing a permanent drop of $2 c_1 c_2 / N$ between the endpoint means $\langle \frac{1}{N}\sum_k c_k^2 \rangle_{T_j}$ and $\langle \frac{1}{N}\sum_k c_k^2 \rangle_{T_{j+1}}$.
