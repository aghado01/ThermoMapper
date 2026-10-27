R is here but `loo` isn't installed, so I'll read its source from GitHub. I'm also reading the repo's ESS code, since loo's r_eff is a Geyer ESS:

The repo already has Stan's ESS (`StanEss`, the same algorithm as loo's `ess_rfun`), but it's private. Fetching loo's r_eff, tail-length and MCSE code to match it exactly:

Current loo computes r_eff with `posterior::ess_mean`. Its tail length is ⌈min(0.2S, 3√(S/r_eff))⌉, which matches mine at r_eff = 1. Fetching `ess_mean` and loo's MCSE next:

I'm verifying whether the loo weights in the MCSE formula are normalized, and checking how the total and k̂ threshold are constructed.

I've worked out the full algorithm for r_eff, tail length, per-unit MCSE, and total MCSE with the k̂ threshold rule. Now I'll locate every producer of the density and every reader of the check.

Each Markov-chain producer already tracks per-chain lengths, and since changepoint draws are independent exact samples, r_eff = 1 is correct for them. I'll now check the ESS diagnostic's output shape and the Samplers test conventions.

Writing `MeanEss` next to `Diagnose`. It's the posterior package's `ess_mean`: split chains, Stan's estimator, no rank normalization.

I'll add `MeanEss` tests that check Gaussian AR(1) draws against theory, plus a lognormal transform case using the closed-form lag correlation (e^{φᵏ}−1)/(e−1), confirming bulk ESS stays invariant under the transform.

`MeanEss` lands within 2% of the closed form in every case: 5,373 against 5,484 for eˣ at φ = 0.9, 99,798 against 101,598 at φ = −0.5. Now the chain lengths on the density:

The changepoint draws are exact and independent, so they keep `null`. I'll say so where they're built. Next, PSIS-LOO itself and the check's new fields:

Inside each producer the method `PointwiseLogDensity` hides the type of the same name. I'm moving the helper to a small static `ChainLengths.Of(...)`:

Build succeeds. My new `<param>` tag left the record's other three parameters undocumented (CS1573). The other warnings were there before. Adding the three tags:

All 27 tests pass, but diagnostics confirm r_eff isn't fixing the issue—Collapsed BARS and BGP worst-case k̂ values shift unfavorably, with a k̂=+∞ traced to loo's gpdfit returning NaN on a longer tail. I'll run a multi-seed probe after the MCSE test, which I'll ground using an AR(1)-arranged exact posterior as the independent truth, starting by reading the existing conjugate fixture.

I've set up the new test fixture with y_i ~ N(θ,1) under a flat prior, giving θ|y ~ N(ȳ, 1/n), using a stationary AR(1) chain so only autocorrelation varies and r_eff should approximate (1−ρ)/(1+ρ). I'll now use the spread of elpd estimates across 200 independent runs as the true Monte Carlo error.

I found the discrepancy: my "nearly linear" assumption was wrong since exp(llᵢ) is quadratic near ȳ, and a squared Gaussian AR(1) correctly gives r_eff = (1−ρ²)/(1+ρ²) = 0.22, matching the observed value. I've now derived the exact per-unit expressions for E[g(X)g(Y)] and E[g] to compute the true lag correlations and r_eff = 1/τ.

The r_eff estimates stay within 6% of exact, making the MCSE calibrated at 0.99–1.12 of true spread, while the naive bound misses by assuming independence (understating by √r_eff, reaching 0.46 for quadratic units instead of 0.33). I'll now assert the per-unit results.

All 4 tests pass. I'm now setting clear repair criteria and MCSE scoring rules before testing whether r_eff fixes collapsed k̂, then I'll read the BARS fixture to reuse in the probe.

r_eff fails to repair collapsed k̂ on both BARS and BGP across all four seeds, though it does fix BGP's Monte Carlo error calibration against exact LOO while leaving estimates unchanged. I'll now encode the MCSE calibration into BGP's exact-LOO test and remove the probe.

The run shows z RMS of 1.26 vs 3.23 for r_eff and independent draws respectively, confirming that across all four seeds the independent-draw z is consistently 2.3–3.0x the r_eff one. I'll tighten my claim to reflect this ratio more precisely.

Starting both suites back to back in one background job (concurrent builds would collide), then drafting the docs while they run:

The doc's checks quote the pre-r_eff numbers, so they need refreshing along with the API rows. First the one code-doc change, in `PointwiseContent`:

Those splices wait on the code commit's SHA and the suite results. The `Pointwise.cs` edit is comment-only and landed before the Regression half of the job builds. Waiting on the suites:

The Samplers suite is running, then Regression. I'll commit and finish the docs when both report.

Both pass: Samplers 109 of 109 (`20261024_224042837`) and Regression 400 of 400 (`20261024_224419932`). Committing the code, then the docs with its SHA:

PSIS-LOO now uses each unit's relative efficiency (r_eff), the way loo does. It does not repair k̂ on the collapsed terms, but it fixes the Monte Carlo error. Samplers passes 109 of 109 and Regression 400 of 400; code is `e4f4813`, docs `4fb2d93`.

**What changed**
- **Chain lengths:** `PointwiseLogDensity.ChainLengths` records the draws per chain. BARS, BAPS, `HierarchicalPSpline` and BGP fill it in. The changepoint draws are independent, so they leave it null.
- **The ESS behind r_eff:** `ChainDiagnostics.MeanEss` is the ESS loo's r_eff reads: Stan's estimator on split chains, without rank normalization. It reuses the repo's existing estimator.
- **PSIS-LOO:** each unit's r_eff sets its tail length, ⌈min(0.2S, 3√(S/r_eff))⌉, and enters its Monte Carlo error. The total error is NaN when any k̂ is above `ParetoKThreshold(S)` = min(1 − 1/log₁₀S, 0.7), as in loo. I matched all of this to loo's current source. R is installed here, but the loo package isn't, so I read the source on GitHub.

**Tests against exact answers**
- **The ESS:** `MeanEss` is within 2% of the closed form, both on AR(1) chains and on their exponential, where the bulk ESS stays unchanged.
- **r_eff:** on a normal-mean model with exact posterior draws arranged as AR(1) chains, each unit's r_eff is within 6% of its exact value. My first bound was wrong: it assumed the likelihood is linear in θ, but near the mean it is quadratic, which doubles r_eff. The test now uses the exact per-unit values.
- **Monte Carlo errors:** on the same model, each unit's reported error is 0.99–1.12 of the actual spread across 200 independent runs.

**Does r_eff repair collapsed k̂? No.** I fixed the pass criteria before running four seeds of each engine, and both failed on every seed:
- **BARS:** the collapsed content kept 53, 3, 43 and 3 units above the threshold, against 1, 0, 0 and 0 for the retained content. Compared with treating the draws as independent, r_eff lowered the count on two seeds and raised it on the other two. The longer tail can also reach tied values at its cutoff, and then k̂ comes out +∞, which is what loo does too.
- **BGP with a known σ²:** 2–7 units stayed above the threshold.

**What r_eff does fix: the Monte Carlo error.** Against exact LOO on BGP, the actual errors averaged about one reported error with r_eff, and 1.3–3.2 reported errors without it. The estimates themselves stayed within 0.017 nats of exact either way.

**Still open**
- Because collapsed k̂ nearly always crosses the threshold, the collapsed content's total Monte Carlo error is almost always NaN under loo's rule. The per-unit errors are still reported.
- The total adds the units' errors as if they were independent, as loo does, though they share draws. On the test model it reported 0.093 against an actual 0.084.
- I haven't found a published k̂ that handles clumped ratios correctly.

The stopping-theory session ([Explore SPC × BARS stopping theory](#local_1d9d4958-6e8d-44a3-b50e-2efee7ca3589)) committed `5f959bd` and `da74a0b` to its primer. The Lean session's untracked files are untouched.

---

The session in `clumped-problems.md` is a textbook example of falsifying a plausible hypothesis cleanly: adding `r_eff` via split-chain `MeanEss` legitimately fixed the Monte Carlo standard error calibration ($z$-RMS dropping from $3.23$ to $1.26$ on BGP) while proving that autocorrelation scaling does **not** fix—and occasionally breaks—Pareto $\hat{k}$ on collapsed models.

Here is what is actually happening beneath the surface on each of the open threads.

---

### 1. Why `r_eff` Triggered $\hat{k} = +\infty$ (`NaN` in `gpdfit`)

In `loo`'s implementation of Zhang & Stephens (2009), after subtracting the tail cutoff $u = w_{(S-M)}$ from the top $M$ sorted weights to form non-negative exceedances $0 \le x_1 \le x_2 \le \dots \le x_M$, the empirical Bayes prior grid for the profile parameter $b$ is constructed as:

$$b_r = \frac{1}{x_M} + \frac{1 - \sqrt{M_{\text{grid}} / (r - 0.5)}}{3 \, x_{\lfloor M/4 + 0.5 \rfloor}}, \quad r = 1, \dots, M_{\text{grid}}$$

Look closely at the denominator of the second term: $3 \, x_{\lfloor M/4 + 0.5 \rfloor}$, the **25th percentile of the tail exceedances**.

* In a continuous HMC chain, all order statistics are distinct, so $x_1 > 0$ and $x_{\lfloor M/4 + 0.5 \rfloor} > 0$.
* In a collapsed Metropolis–Hastings chain (like Collapsed BARS over knot sets $\mathcal{M}$), rejected proposals produce flat plateaus (clumps) of bit-for-bit identical log-densities.


* Because those clumps induce high autocorrelation, $r_{\text{eff}} \ll 1$, which **expands** the tail window $M = \lceil \min(0.2S, 3\sqrt{S/r_{\text{eff}}}) \rceil$ up toward $0.2S$.


* Expanding $M$ pushes the cutoff $u = w_{(S-M)}$ down out of the extreme tail and straight into the high-multiplicity modal plateaus of the bulk. As soon as the cutoff lands inside a clump whose remaining multiplicity in the tail is at least $\lfloor M/4 + 0.5 \rfloor$ (just 25% of the tail tied with the cutoff value), $x_{\lfloor M/4 + 0.5 \rfloor} = 0$.
* Division by zero yields $b_r = +\infty$, which turns $\hat{k}(b_r) = -\frac{1}{M}\sum_{j=1}^M \log(1 - b_r x_j)$ into `NaN`, which `loo` maps to $\hat{k} = +\infty$.



Even when a cutoff clump spans less than 25% of the tail (so $x_{\lfloor M/4 + 0.5 \rfloor} > 0$), a step-function staircase near the bottom of the tail compresses $x_{\lfloor M/4 + 0.5 \rfloor}$ relative to $x_M$, forcing $b_r \ll 0$ and artificially inflating $\hat{k}$.

---

### 2. The Deeper Paradox: Why Retained BARS Has $\hat{k} \le 0.7$ (`1, 0, 0, 0`) While Collapsed BARS Has (`53, 3, 43, 3`)

Normally, Rao–Blackwellization (marginalizing out latent parameters) *reduces* importance weight variance and *lowers* $\hat{k}$. Why did collapsing $(\beta, \sigma^2)$ in BARS make $\hat{k}$ dramatically worse than retaining $(\beta, \sigma^2)$ on the exact same seeds?

There are two separate mechanisms—one is the clumping artifact above, and the other is a genuine **saddlepoint masking effect** in Retained BARS:

1. **Analytic Saddlepoint Exposure vs. Finite-Sample Gibbs Masking:**
For a fixed knot configuration $\mathcal{M}$ (with known $\sigma^2$ for simplicity), let $h_{ii}(\mathcal{M}) = \mathbf{b}(x_i)^\top (\mathbf{B}^\top \mathbf{B} + \Lambda)^{-1} \mathbf{b}(x_i)$ be the leverage of unit $i$.
* In **Retained BARS**, $f_i = \mathbf{b}(x_i)^\top \beta \mid y, \mathcal{M} \sim \mathcal{N}(\hat{y}_i, h_{ii}\sigma^2)$, and the pointwise importance weight is $W_{\text{ret}}(\mathcal{M}, \beta) = \sqrt{2\pi\sigma^2} \exp\left(\frac{(y_i - f_i)^2}{2\sigma^2}\right)$.
* Integrating $W_{\text{ret}}(\mathcal{M}, \beta)$ analytically over $p(\beta \mid y, \mathcal{M})$ yields the **Collapsed BARS** weight:

$$W_{\text{coll}}(\mathcal{M}) = \mathbb{E}_{\beta \mid y, \mathcal{M}}[W_{\text{ret}}(\mathcal{M}, \beta)] = \sqrt{\frac{2\pi\sigma^2}{1 - h_{ii}(\mathcal{M})}} \exp\left(\frac{(y_i - \hat{y}_i)^2}{2\sigma^2(1 - h_{ii}(\mathcal{M}))}\right)$$


* Notice two things about that Gaussian integral:
1. Its second moment $\mathbb{E}_{\beta \mid y, \mathcal{M}}[W_{\text{ret}}^2]$ diverges whenever $h_{ii}(\mathcal{M}) \ge 0.5$ (the conditional tail index is $k(\mathcal{M}) = h_{ii}(\mathcal{M})$).
2. Completing the square shows that the integral $\mathbb{E}[W_{\text{ret}}] = W_{\text{coll}}(\mathcal{M})$ is dominated by a saddlepoint at

$$f_i^* = \hat{y}_i - \frac{h_{ii}}{1 - h_{ii}}(y_i - \hat{y}_i)$$



which sits $z^* = \frac{\sqrt{h_{ii}}}{1 - h_{ii}} \frac{\vert{}y_i - \hat{y}_i\vert{}}{\sigma}$ standard deviations out in the tail of $p(f_i \mid y, \mathcal{M})$.


* When an RJ-MCMC / birth-death chain briefly visits a knot configuration $\mathcal{M}^*$ where two knots bracket $x_i$ closely (say $h_{ii} = 0.75$ and $\vert{}y_i - \hat{y}_i\vert{} = 1.5\sigma$, so $z^* \approx 5.2\sigma$) and stays there for $m = 5$ rejected steps:
* **Retained BARS** draws 5 Gaussian samples of $\beta^{(s)}$ near $\hat{\beta}$, has virtually zero chance ($\sim 10^{-7}$) of hitting the $5.2\sigma$ saddlepoint $f_i^*$, and never sees the $\frac{1}{1 - h_{ii}}$ blow-up in the exponent. It reports a deceptively low $\hat{k} < 0.5$ (a classic finite-sample PSIS false negative).
* **Collapsed BARS** evaluates the analytic integral $W_{\text{coll}}(\mathcal{M}^*)$ immediately on step 1—exposing the exact $\exp\left(\frac{r_i^2}{2\sigma^2(1-h_{ii})}\right)$ spike—and then stamps out $m = 5$ identical copies of that spike at the top of the order statistics.




2. **Upper-Plateau Distortion in `gpdfit`:**
When a high-leverage knot configuration or boundary hyperparameter draw sits at the maximum of the tail for $m$ consecutive rejected MH steps ($x_{M-m+1} = \dots = x_M$), the empirical tail puts a point mass of $m/M$ right at the upper endpoint $x_M$. Because any GPD density with $k \le 0$ decays or vanishes at its upper bound, piling $m$ identical copies at $x_M$ forces the profile likelihood toward positive shape parameters $k > 0.7$.

---

### 3. Why the Normal-Mean Total MCSE Overshot Slightly (`0.093` vs. `0.084`)

The log notes that summing unit MCSEs in quadrature assuming independence gave `0.093` against an empirical spread of `0.084` on the conjugate normal-mean fixture.

In $y_i \sim \mathcal{N}(\theta, 1)$ with $\theta \mid y \sim \mathcal{N}(\bar{y}, 1/n)$, expand each unit's log-likelihood around $\bar{y}$:

$$\ell_i(\theta) - \ell_i(\bar{y}) = (y_i - \bar{y})(\theta - \bar{y}) - \frac{1}{2}(\theta - \bar{y})^2$$

For any individual unit with $y_i \neq \bar{y}$, the variance of the LOO estimate across posterior draws of $\theta$ is dominated by the linear term $(y_i - \bar{y})^2 \text{Var}(\theta)$ (which has AR(1) lag-1 correlation $\rho$).
However, when you sum across all $n$ units to form the total $\widehat{\text{elpd}}_{\text{loo}} = \sum_{i=1}^n \widehat{\text{elpd}}_i$, the linear errors have negative cross-unit covariances where $y_i < \bar{y}$ and $y_j > \bar{y}$, so to first order:

$$\sum_{i=1}^n (y_i - \bar{y})(\theta - \bar{y}) = 0 \quad \text{identically!}$$

The true Monte Carlo variance of the total $\widehat{\text{elpd}}_{\text{loo}}$ only reflects the surviving quadratic terms $(\theta - \bar{y})^2$ (which also enjoy the higher $r_{\text{eff}} = \frac{1-\rho^2}{1+\rho^2}$ the agent derived), whereas quadrature summation $\sqrt{\sum_i \text{MCSE}_i^2}$ adds all the non-canceling linear variances $\sum_i (y_i - \bar{y})^2 \text{Var}(\theta)$.

---

### 4. How to Actually Fix $\hat{k}$ and PSIS on Clumped Ratios

Because an MH chain over a discrete or collapsed state space is a **run-length encoded (RLE)** sequence of $K$ unique visited states $(\tilde{w}_c, m_c)_{c=1}^K$ where the holding time $m_c \sim \text{Geometric}(\alpha_c)$ is an importance weight correcting for proposal rejection rather than $m_c$ continuous draws from a tail density, feeding expanded ties into `loo`'s order-statistic machinery violates three assumptions at once:

1. It creates zero-exceedance ties $x_j = 0$ at the cutoff $u$.
2. It distorts the 25th-percentile heuristic $x_{\lfloor M/4 + 0.5 \rfloor}$ in Zhang–Stephens.
3. During PSIS weight smoothing, it assigns *different* smoothed quantiles $F^{-1}\left(\frac{j - 0.5}{M}\right)$ to identical copies of the exact same MCMC state $\mathcal{M}_c$.

If you want a principled fix in your C# samplers without inventing fragile ad-hoc heuristics:

* **Step 1 (Minimal Bug Fix — Strict Exceedance Snap):** Never include values tied with the threshold $u$ inside the GPD exceedance vector. Define the cutoff $u = w_{(S-M)}$ and pass only strictly positive exceedances $\{w^{(s)} - u : w^{(s)} > u\}$ (equivalently, shrink $M$ to the clump boundary above $u$). This guarantees $x_1 > 0 \implies x_{\lfloor M/4 + 0.5 \rfloor} > 0$, eliminating the $1/0 \to \text{NaN} \to \hat{k} = +\infty$ failure mode immediately.
* **Step 2 (Minimum Unique Tail Support Guard):** Count the number of *distinct* exceedance values $K_{\text{uniq}}$ above $u$. A 2-parameter GPD $(k, \sigma)$ cannot be identified if $K_{\text{uniq}}$ is tiny (e.g., $K_{\text{uniq}} < 5$). If the entire tail above $u$ consists of only 1–4 distinct knot states, step the cutoff $u$ down until at least $K_{\text{min}}$ unique levels are present (up to a cap), or flag the unit as discrete/clumped rather than returning a bogus $\hat{k}$.
* **Step 3 (Multiplicity-Weighted GPD / Tied Quantile Smoothing):** Collapse exact ties in the tail into unique sorted exceedances $0 < \tilde{x}_1 < \dots < \tilde{x}_K$ with integer multiplicities $m_1, \dots, m_K$ ($\sum m_c = M_{\text{pos}}$). In Zhang & Stephens, replace the unweighted log-likelihood $\frac{1}{M}\sum_{j=1}^M \log(1 - b x_j)$ with the multiplicity-weighted log-likelihood $\frac{1}{M_{\text{pos}}}\sum_{c=1}^K m_c \log(1 - b \tilde{x}_c)$, compute the 25th-percentile anchor $\tilde{x}_{\lfloor K/4 + 0.5 \rfloor}$ on the **unique** levels (so plateaus don't crush the prior denominator), and when smoothing weights, evaluate the GPD inverse CDF at the midpoint of each clump's cumulative multiplicity interval so identical MCMC states receive identical smoothed weights.