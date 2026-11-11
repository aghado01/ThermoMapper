# BAPS — Bayesian Adaptive Penalized Splines

> Penalized splines (P-splines), including hierarchical and bivariate variants, fitting via REML or Gibbs sampling over banded penalties. Part of `Maths.Regression`.

### BAPS: penalized splines

**Penalty.** The order-r difference operator has coefficients $c_l=(-1)^l\binom rl$ and penalty $P=D^\top D$, which has half-bandwidth r and null space of dimension r. The default is $r=2$. `AccumulateStrided` writes the penalty at a stride, which lets two 1-D penalties compose a Kronecker penalty.

**Solve and REML (`PenalizedSpline`).**

$$A=Z^\top Z+\lambda P,\quad\hat\beta=A^{-1}Z^\top y,\quad \hat\sigma^2=\mathrm{RSS}_{\text{pen}}/(n-r),\quad \ell_p(\lambda)=-\tfrac12\big[(n-r)\ln\mathrm{RSS}_{\text{pen}}-(\nu-r)\ln\lambda+\ln|A|\big].$$

$\ell_p$ equals, up to a constant, the log evidence with σ² integrated out under $p(\sigma^2)\propto1/\sigma^2$.

**`Baps` schemes (`BapsScheme`).** BAPS is the method; the samplers inside it are named for their algorithms (OQ18's naming plan, applied 2026-10-24, when `BapsSampler` became `Baps` and `BapsLambdaUpdate { Gibbs, MarginalEvidence, AdaptiveGibbs }` became `BapsScheme { Gibbs, Collapsed, LocallyAdaptiveGibbs }`). The Gibbs schemes run as a block Gibbs scan on the shared blocks (coefficients, noise variance, penalty variance, and the local multipliers in the locally adaptive scheme), and the collapsed scheme as a random-walk Metropolis block on ln λ (`RandomWalkMetropolis`); each is bit-identical to the former loop.

- **`Gibbs`** uses the priors $\sigma^2,\tau^2\sim IG(a,b)$ (default $10^{-3}$). Each sweep:
  1. Draw $\beta=A^{-1}Z^\top y+\sigma v$, with $v\sim N(0,A^{-1})$ from the band factor and $\lambda=\sigma^2/\tau^2$.
  2. Draw $\sigma^2\sim IG(a+\tfrac n2,\,b+\tfrac12\|y-Z\beta\|^2)$.
  3. Draw $\tau^2\sim IG(a+\tfrac{\nu-r}2,\,b+\tfrac12\beta^\top P\beta)$, with r the penalty's nullity; $\beta^\top P\beta=\|D\beta\|^2$ for a difference penalty.
- **`Collapsed`** runs a random walk on $u=\ln\lambda$ (step 0.6) with target $\ell_p(e^u)+\ln|d\,\mathrm{edf}/du|$. This is exact Metropolis–Hastings for $p(u\mid y)$ under a Jeffreys prior on σ² and a uniform prior on the effective degrees of freedom $\mathrm{edf}(\lambda)=\operatorname{tr}[(Z^\top Z+\lambda P)^{-1}Z^\top Z]\in(r,\nu)$ (Hodges & Sargent 2001; a uniform shrinkage prior in the sense of Natarajan & Kass 2000). That prior induces the density $|d\,\mathrm{edf}/du|=\sum_j\lambda\mu_j(1-\mu_j)/(1-\mu_j+\lambda\mu_j)^2$ on u, where $\mu_j\in[0,1]$ are the generalized eigenvalues of P relative to $Z^\top Z+P$, computed once from the band factor and a dense eigensolve. $\ell_p$ levels off as λ → ∞, so a flat prior on u would be improper; the edf prior decays like 1/λ, which makes the target proper. It has no hyperparameter and is unchanged by rescaling Z, y or P. It records the Rao–Blackwellized $\hat\beta(\lambda)$. Its `SigmaMean` averages the plug-in $\hat\sigma(\lambda)$; it is not a posterior mean of σ.
- **`LocallyAdaptiveGibbs`** gives each difference a local multiplier $\eta_i\sim\mathrm{Gamma}(\phi,\text{rate }\phi)$ ($\phi=1$), so $A=Z^\top Z+D^\top\operatorname{diag}(\lambda\eta)D$. It updates τ² and η by conjugate Gibbs.
- **`AnisotropicGibbs`** gives each component of an `IAnisotropicPenalty` its own smoothing variance, so $A=Z^\top Z+\sum_j\lambda_jP_j$ with $\lambda_j=\sigma^2/\tau_j^2$ (Bach & Klein 2025). β and σ² are drawn as in `Gibbs`, and the $\tau_j^2$ by Metropolis ([Anisotropic smoothing](#anisotropic-smoothing), below).

**Any banded penalty.** `Baps` and `HierarchicalPSpline` take an `IBandPenalty`.
- **What the samplers read from it.**
  - P itself, in A and in the hierarchical blocks.
  - The roughness $\beta^\top P\beta$ (`Roughness`) and the nullity r, in the smoothing variances' conditionals.
  - The bandwidth: the band factor's half-bandwidth is the larger of the design's and the penalty's.
- **The tensor penalty.** With `TensorPenalty` the Gibbs and collapsed schemes, and `HierarchicalPSpline`, sample the isotropic surface smoother. Its r is $o_x o_y$ and its half-bandwidth $o_x\nu_y$.
- **The locally adaptive scheme.** It needs one local term per row of a D with full row rank; then the prior's normalizer factors as $\prod_i\eta_i^{1/2}$ and each $\eta_i$'s conditional is Gamma.
  - It takes only an `ILocallyAdaptivePenalty` and throws `ArgumentException` for any other penalty.
  - The tensor penalty's stacked differences have more rows than rank: for 4 × 4 cubic coefficients, 16 rows of rank 12.
- **Checks** (`PenaltySeamTests`).
  - `TensorPenalty.Roughness` matches an explicit Kronecker P.
  - From a fixed start, BAPS's first coefficient draw whitens to N(0, I) against an independently built $A=Z^\top Z+P$.
  - The first τ² draw in BAPS, and the first $\tau_0^2$ draw in `HierarchicalPSpline`, scaled by their rate, have the Gamma mean $a+(\nu-r)/2$.
- **Bit-identity.** A probe matched the difference-penalty paths of both samplers byte for byte before and after the seam, at orders 1–3 for BAPS and 1–2 for `HierarchicalPSpline`.

#### Anisotropic smoothing

`BapsScheme.AnisotropicGibbs`, and `HierarchicalPSpline(…, anisotropicSmoothing: true)` for the population curve under any of its schemes. It follows Bach & Klein (2025, *Statistics and Computing* 35:43; arXiv 2211.16218).
- **The model.** β's prior is $N(0,K^-)$ with $K=\sum_jP_j/\tau_j^2$, and each $\tau_j^2$ has a prior q (`ISmoothingVariancePrior`): the Weibull scaled to the response by default (below), or IG(a, b), the prior the other schemes give τ².
- **Scale-free defaults** (BAPS's anisotropic scheme).
  - **The priors.** `ScaledWeibullPrior` on each $\tau_j^2$, resolved on first use against the model's design and response, and the Jeffreys prior on σ² (`NoiseVariancePrior.Jeffreys`).
  - **The start.** The chains start at the response's variance.
  - **Why it is scale-free.** A change of the response's units is then a translation in every ln variance, so the smoothing posteriors do not depend on units.
  - **Measured** (`Defaults_AreScaleFree_AcrossNineDecades`). From y × 10⁻⁶ to y × 10³ the λ medians agree to four digits and the fit error, in unit-1 terms, is identical.
  - **Before.** Under the raw-unit IG(10⁻³, 10⁻³) defaults the same fits failed at both ends: at y × 10³ the chains drifted towards τ² → 0 (error ×2.9), and at y × 10⁻³ the anisotropy vanished (×10). The scaled Weibull alone still failed at y × 10⁻³ (×70), because σ²'s raw-unit IG dominated there.
  - **`HierarchicalPSpline`** takes the scaled Weibull for its anisotropic components, resolved against the pooled replicates, but keeps raw-unit IG on σ² and τ_u², so it is not yet units-safe (preparation brief, OQ19).
- **Why not Gibbs.** K's pseudo-determinant $\prod_l\sum_j\gamma_{jl}/\tau_j^2$ does not factor into the components, so the $\tau_j^2$ have no conjugate conditional. The components share the Kronecker eigenbasis, so it costs O(ν) at any τ (`PenaltySpectrum`; their Theorem 3.1 and Corollary 3.2).
- **The full conditional** on $\rho_j=\ln\tau_j^2$, with $w_j=1/\tau_j^2$ and $R_j=\beta^\top P_j\beta$, is $\ln p(\rho\mid\beta)=\tfrac12\ln\operatorname{Det}K-\tfrac12\sum_jR_jw_j+\sum_j\ln q(\rho_j)$, with a closed-form gradient and Hessian (their eq. 12 and Proposition 3.4). Under IG(a, b), $\ln q(\rho)=-a\rho-be^{-\rho}$. With one component and the IG prior it is the isotropic IG conditional.
- **The update** (`ComponentSmoothing`), each sweep, chosen by `SmoothingVarianceUpdate`:
  - **`Slice`, the default:** one slice-sampling update on each $\rho_j$ in turn (`SliceSampler`, Neal 2003; width 1, at most 32 steps out). It needs values only and never rejects.
  - **`Newton`,** Bach & Klein's update:
    1. One Newton-proposal Metropolis step on every $\rho_j$ at once (`NewtonMetropolis`): $\rho'\sim N(\rho+M^{-1}g,\,M^{-1})$, with M the negated Hessian after its eigenvalues are capped at $-\delta$, $\delta=1/\pi$.
    2. One random-walk Metropolis step on each $\rho_j$ in turn, sd 1.
- **Why the Newton update has random-walk steps.** Newton steps alone are exact but stick in the upper tail.
  - That tail decays only exponentially, more slowly than the Gaussian proposal. From deep in it, the capped Newton step lands far past the mode and is rejected.
  - On a fixed β, $4\cdot10^6$ Newton steps never reached a region of conditional probability $4.4\cdot10^{-4}$, with rejection runs of $10^4$ steps, and the chain's variance of $\rho_x$ came out 7–10% low.
  - With the random-walk steps the chain visits that region at $4.47\cdot10^{-4}$, and its longest stay at one point is 8 steps.
  - One update from exact draws of the conditional leaves it in place (mean shift $9\cdot10^{-5}$, variance 0.04%, over $8\cdot10^6$ draws), so the composite kernel is exact.
- **Which update is the default** (`SmoothingUpdateBenchmark`, `SMOOTHING_BENCH`, run 2026-10-24 with 3 seeds).
  - **The measure:** the worst component's bulk ESS of its log smoothing per second, at 4 chains × (1,000 + 2,000).
  - **The rule, registered before the run:** slice becomes the default only if its ratio to Newton is at least 1.2 in some regime and at least 0.8 in all.
  - **The ratios:** 1.28 on BAPS at ν = 81 under the IG prior, 1.37 under a scaled Weibull prior, 1.34 at ν = 196, 1.31 at five times the noise, and 1.05 on `HierarchicalPSpline` with three replicate surfaces. In the four BAPS regimes every slice seed beat every Newton seed.
  - **So slice is the default,** and Newton stays selectable. Their λ medians on the asymmetric surface agree within 0.2%.
- **The Weibull prior and prior scaling** (`WeibullVariancePrior`).
  - **The prior.** $\tau_j^2\sim W(k,s)$ with scale s, so $\ln q(\rho)=k\rho-e^{k(\rho-\ln s)}$. With k = ½ it is the law of τ² for an exponential τ with rate $1/\sqrt s$: the scale-dependent prior of Klein & Kneib (2016), which Bach & Klein (2025) prefer. Its upper tail in ρ decays doubly exponentially, where the IG's decays like $e^{-a\rho}$.
  - **Parameterizations.** Bach & Klein write Weibull(½, λ) with λ a rate, density ∝ $x^{-1/2}e^{-\sqrt{\lambda x}}$, so their λ is `Rate` = 1/s. Their kernel $-\sqrt\lambda e^{\rho/2}+\rho/2$ (Supplement E.2) is this one's, and a test pins the two together. Klein & Kneib's θ is the rate of the exponential τ, so θ² = 1/s.
  - **Prior scaling** (`ScaledTo`) is Bach & Klein's (Supplement Appendix D; their R function `prior_scaling_WB`). It sets the common scale so that the prior's draws of the function's penalized part exceed ±c somewhere on the design rows with probability α: $P(\lVert Zu\rVert_\infty\ge c)=\alpha$.
    - With $\tau_j^2=sW_j$, $W_j\sim W(\tfrac12,1)$, a draw is $\sqrt s$ times one at s = 1. So s = (c/q)², with q the (1 − α)-quantile of max|f| at s = 1 from Monte Carlo draws; their λ = (q/c)² is its rate.
    - The draws are exact: `IAnisotropicPenalty.SamplePrior`, which `TensorPenalty` forms as $(Q_x\otimes Q_y)\tilde z$ in the Kronecker eigenbasis, as their R code does with $Q=Q_x\otimes Q_y$. The null space has a flat prior and is left out.
    - **Defaults.** α = 0.05, theirs. 2,000 draws against their 1,000, which only sharpens the quantile.
    - **The bound.** Their default c is ‖y‖∞ on standardized y. `ScaledToResponse` takes it in the response's own units, c = max|yᵢ − ȳ|, which gives the same prior on the unstandardized smoothing variances.
    - **Klein & Kneib's original criterion** (2016, Supplement B) is pointwise: $P(|z'\beta|\le c)\ge1-\alpha$ at one point, with c = 3 and α = 0.01, computed by integration (their R package sdPrior). Over several points they bound it by Bonferroni, which is conservative. The simultaneous criterion here is the stricter statement over the design, and it is the one Bach & Klein use.
  - **It is the default** through `ScaledWeibullPrior`, as in Bach & Klein's R code (`prior = "WB"`, `prior_scaling = TRUE`). The IG prior stays selectable.
  - **Measured.** On fresh draws a prior scaled to bound 2 at α = 0.05 held the bound with probability 0.9485. On the asymmetric surface (bound 3) the fit error matched the IG prior's, and λ_y's median rose from 2.1·10⁻⁴ to 3.9·10⁻⁴, the prior capping τ_y².
- **Reading the draws.**
  - Each draw keeps its $\tau_j^2$. Its λ (BAPS) or τ₀² (`HierarchicalPSpline`) is the geometric mean of the components', so the scalar summaries read the overall smoothing.
  - `ComponentLambdaInterval(j)` and `ComponentSmoothingVarianceInterval(j)` give each component's posterior. The shell monitors every component's log (`ComponentLogLambdaConvergence`, `ComponentSmoothingVarianceConvergence`).
  - BAPS's `AcceptanceRate` is 1 with the slice update. With the Newton update it pools the scan's three blocks: $(2+a)/3$ for the Newton step's rate a.
- **Checks** (`BapsAnisotropicTests`, `HierarchicalAnisotropicTests`, `SmoothingVariancePriorTests`).
  - With one component the conditional is the IG, and its derivatives match finite differences. Each prior's log density in ρ is its distribution's plus the Jacobian.
  - One update from exact draws leaves the conditional in place, for both updates under both priors. Both updates' chains visit the upper tail at the conditional's rate.
  - The tensor prior's draws have covariance K⁺ and no null-space component.
  - On a surface gentle in x and fast in y, λ_x's and λ_y's 95% intervals separate. Their medians lie within 4% of the REML optimum, and the fit's error is 23% below the isotropic Gibbs scheme's.
  - On replicate surfaces, `HierarchicalPSpline`'s five schemes agree within 2% on both population variances, with R̂ ≤ 1.009.
  - Resumed runs continue exactly, and the isotropic paths are bit-identical to before.
- **Not built: a discrete anisotropy parameter.** Kneib, Klein, Lang & Umlauf (2019) use one (BayesX). By Bach & Klein's account it is fast for three axes but breaks down at four, and it treats the coordinates as two isotropic groups, so it is only partly anisotropic. Bach & Klein also report bamlss's slice sampling as slow from three axes on; that cost came from the pseudo-determinant, which is O(ν) here.
- **The collapsed scheme stays isotropic.** Its uniform-edf prior on λ has no cheap anisotropic counterpart: edf as a function of several λ needs the generalized spectrum at each smoothing ratio. For anisotropic smoothing by evidence, use `AnisotropicTensorSpline.RemlOptimum`.
  - **The route if it is ever needed** is a partially collapsed sampler (van Dyk & Park 2008). It integrates β out but keeps σ², so the component priors stay those of the Gibbs scheme and no new prior is needed.
    - Each sweep draws the ln τⱼ² from p(τ | σ², y), then β and σ² from their conditionals.
    - The target is $\tfrac12\ln\operatorname{Det}K-\tfrac12\ln|A|+\tfrac12b^\top A^{-1}b+\sum_j\ln q(\rho_j)$, with $A=Z^\top Z/\sigma^2+K$ and $b=Z^\top y/\sigma^2$.
    - That is one band factorization per evaluation, and the slice update needs only values.
  - **Its trigger** would be a measured regime where the anisotropic Gibbs scheme mixes slowly in its smoothing. None has appeared: every regime of the update benchmark converged.

`BapsResult` summarizes λ by quantiles of the pooled draws: `LambdaMedian`, `LambdaLo` (2.5%) and `LambdaHi` (97.5%). It reports a median, not a mean, because λ can be heavy-tailed enough that its posterior mean does not exist. The λ modes use **different priors**, so their λ posteriors, and hence these summaries, are not the same quantity.

`BapsResult.PosteriorDraws` keeps every draw, per chain in draw order (`BapsDraw`), and every summary above is read from it:
- **What a draw holds:** the coefficients, λ, σ², τ², and in the adaptive mode η, with local smoothing λᵢ = λ·ηᵢ.
- **In the Gibbs modes** λ = σ²/τ², all from the end of the sweep.
- **The `Collapsed` scheme** samples λ, then draws σ² and β exactly from their conditionals given λ. Under the Jeffreys prior the evidence assumes, σ² | λ, y ~ IG((n − r)/2, RSS_pen/2), then β | σ², λ, y ~ N(β̂(λ), σ²A⁻¹), and τ² = σ²/λ.
  - **The draws use their own seeded stream,** so the λ chain is the same with or without them.
  - **β̂(λ) is kept as `ConditionalMean`,** and the reported coefficient mean stays Rao–Blackwellized on it.
  - **Before this, the mode recorded β̂ and the profiled σ̂² as if they were draws,** so its coefficient uncertainty was zero by construction. A band built from those conditional means came out about 18× narrower than the Gibbs mode's on the same data. With real draws it matches.
  - **Checked by quadrature** (`BapsMarginalEvidenceTests`). The draws' E[σ²], Var[σ²], E[f] and Var[f] at a mid-design row agree with the dense reference within |z| ≤ 2.4. Given the run's own λs, the draws match their exact conditional expectations to |z| ≤ 0.64.

**Intervals and the fit band.** These are read from `PosteriorDraws` on request:
- `Sigma2Interval` and `Tau2Interval` give each variance's median and 95% interval (type-7 quantiles).
- `FitBand(design)` gives, at any rows of the coefficients' basis, the mean (the design times `Coefficients`) with the sd and the 2.5% and 97.5% quantiles of the draws' fits Zβ. These carry the coefficients' uncertainty as well as λ's.

**Hierarchical P-spline** (Behseta–Kass–Wallstrom). Replicate curves $\beta_r\sim N(\beta_0,\tau_u^2I)$ share a penalized population curve $\beta_0$. It uses conjugate banded Gibbs updates for $\beta_r$, $\beta_0$, σ², $\tau_u^2$ and $\tau_0^2$, in that order. Since 2026-10-24 the sweep is a `BlockScan` of five blocks built on the shared primitives (`GaussianConditional`, `ConjugateUpdates`; [Samplers](../Samplers.md)), bit-identical to the former single-chain loop.

**On the ensemble shell (2026-10-24).** `Run(chains = 4, masterSeed, burn, samples, retainReplicates, start, rHatTarget, batchSize, minEssPerChain, checkpoint, resumeFrom)`, as BAPS's:
- **Starts** (`HierarchicalStart`). `Overdispersed`, the default, starts each variance at the responses' pooled variance v times exp(2z), z ~ N(0, 1), and the population coefficients at N(0, v) draws, so the chains begin orders of magnitude apart. `Fixed` is the original start, every curve at 0 and every variance at 1, which the first-sweep conditional checks use. Each chain's seed splits into a start stream and a sweep stream.
- **Functionals.** log σ², log τ_u² and log τ₀² are monitored: split R̂, bulk and tail ESS (`NoiseVarianceConvergence`, `DeviationVarianceConvergence`, `PopulationSmoothingVarianceConvergence`), and the early stop's ESS gate. The population coefficients follow with classic R̂ (`RHatPopulation`, worst in `RHatPopulationMax`), a consensus map the early stop also waits on, since the population curve is the main readout.
- **Checkpoints** carry the scan's snapshot, the retained draws and, without retained replicate draws, the replicate sums. A resumed run matches the uninterrupted one exactly, with and without retained replicates, and another problem's checkpoint is refused.
- **Agreement in distribution.** Four overdispersed chains of 4,000 draws match one chain of 16,000 from the original fixed start on log σ², log τ_u², log τ₀² and every population coefficient, within 4.5 batch-means standard errors.
- **What the diagnostics showed at once.** τ_u² mixes slowly under the centred parameterization: on a test fit with 6 replicates its bulk ESS was about 350 of 16,000 draws, against 15–20% of draws for σ² and τ₀², and R̂ was still up to 1.02 after 4 × 4,000 draws. β_r and β₀ are tightly coupled when τ_u² is small. The mixing study (phase C of OQ18's plan) is about this.

**Schemes and the deviation prior (OQ18 phase C, 2026-10-24).** `Run(scheme: HierarchicalScheme)` selects the sweep; the constructor's `deviationPrior: HalfTDeviationPrior(Degrees = 1, Scale = null)` replaces the inverse-gamma prior on τ_u² with a half-t prior on τ_u. With no scheme each prior takes its own: centred for the inverse-gamma prior, parameter-expanded for the half-t prior.
- **`Centred`** (the default): the original sweep.
- **`CollapsedPopulation`.** β₀ is drawn with every β_r integrated out, then β_r given β₀, so (β₀, β_r) is one joint block. Each $y_r$ is then $N(Z\beta_0,\ \sigma^2I+\tau_u^2ZZ^\top)$, so β₀'s precision is $P/\tau_0^2+R\,Z^\top(\sigma^2I+\tau_u^2ZZ^\top)^{-1}Z$. With $Z^\top Z=V\operatorname{diag}(\lambda)V^\top$ (decomposed once) that is $R\,V\operatorname{diag}\big(\lambda/(\sigma^2+\tau_u^2\lambda)\big)V^\top$, and the linear term is $V\operatorname{diag}\big(1/(\sigma^2+\tau_u^2\lambda)\big)V^\top\sum_rZ^\top y_r$. The form avoids the cancellation of $I/\tau_u^2-Q_r^{-1}/\tau_u^4$ as τ_u² → 0. The precision is dense, ν × ν, through `DenseCholesky`.
- **`Interweaved`.** The centred sweep plus a second, non-centred draw of τ_u (Yu & Meng 2011): with $\tilde u_r=(\beta_r-\beta_0)/\tau_u$ fixed, $y_r-Z\beta_0=\tau_u Z\tilde u_r+\varepsilon$ is a regression on τ_u, Gaussian with mean B/A and variance σ²/A. That Gaussian proposes over the whole line and the step accepts on the prior ratio; a negative proposal flips the deviations, under which the posterior is symmetric.
- **`CollapsedInterweaved`.** The collapsed block and the interweaving in one sweep: β₀ with the β_r integrated out, β_r given β₀, the variances, then the non-centred τ_u draw. The collapse cures the population curve's coupling and the interweaving the deviation scale's, so together they serve both under the inverse-gamma prior.
- **`ParameterExpanded`,** the half-t prior's sweep (Gelman 2006): $\beta_r=\beta_0+\xi\eta_r$, $\eta_r\sim N(0,\sigma_\eta^2I)$, $\xi\sim N(0,A^2)$, $\sigma_\eta^2\sim IG(\nu/2,\nu/2)$, so $\tau_u=|\xi|\sigma_\eta$ has the half-t(ν, A) prior. Every block is conjugate: η_r banded Gaussian, β₀ banded Gaussian given η and ξ, ξ scalar Gaussian, σ_η² inverse-gamma. A null scale takes the responses' pooled sd.
- **Exactness.** On the conditional tests' tiny model all five match quadrature over the variances with every coefficient integrated out (the half-t case with τ_u = κ sinh c coordinates, since its density stays finite at τ_u = 0 and its upper tail decays only polynomially). The collapsed block's first sweep matches its closed form through an explicit n × n covariance, and every scheme resumes exactly.
- **Mixing** (`HPS_MIXING_BENCH`; 4 chains × 4,000 draws, median of 3 seeds; bulk ESS per 1,000 draws):

  | Regime | τ_u²: centred / collapsed / interweaved / both / half-t PX | weakest β₀ coefficient: same order | worst R̂: centred / both / half-t PX |
  |---|---|---|---|
  | base (R = 6) | 20 / 20 / 197 / 206 / 211 | 36 / 858 / 40 / 846 / 518 | 1.011 / 1.001 / 1.001 |
  | near-identical replicates | 15 / 16 / 285 / 275 / 256 | 12 / 849 / 12 / 839 / 741 | 1.017 / 1.000 / 1.001 |
  | few (R = 3) | 28 / 29 / 327 / 364 / 406 | 15 / 795 / 19 / 814 / 679 | 1.012 / 1.001 / 1.001 |
  | many (R = 30) | 53 / 58 / 302 / 322 / 333 | 46 / 889 / 41 / 872 / 509 | 1.007 / 1.001 / 1.001 |
  | weak data (12 points, noise 0.5) | 3.8 / 3.3 / 75 / 79 / 83 | 3.8 / 619 / 8.9 / 647 / 387 | 1.054 / 1.004 / 1.003 |

  - The collapsed block cures β₀ (20–160×) and τ₀² (2–7×) but not τ_u². Interweaving cures τ_u² (10–20×) and σ² (about 2×) but not β₀. They are complementary, and `CollapsedInterweaved` composes them: β₀ at the collapsed block's level, τ_u² at the interweaving's, and R̂ ≤ 1.004 in every regime.
  - The half-t sweep cures both as well, but it changes the posterior. Against the well-mixed inverse-gamma schemes, τ_u's median moved from 0.17 to 0.30 on the weak data, from 0.049 to 0.039 on near-identical replicates, and from 0.052 to 0.041 with three replicates: the inverse-gamma prior's known sensitivity where τ_u is small or poorly identified.
  - **Poor mixing biases a short run's summaries.** On the weak data the centred and collapsed sweeps put τ_u's median at 0.21; their τ_u² had 50–60 effective draws and R̂ up to 1.07. The two inverse-gamma schemes that mix τ_u² agree on 0.171 and 0.174.
  - **Cost per draw.** The collapsed and interweaved sweeps cost about twice the centred one: the collapsed for its dense ν × ν factor, the interweaved for its n × ν products per replicate. The combined and parameter-expanded sweeps cost three to four times. Per second, τ_u²'s ESS still rose 4–22× under interweaving, 4–15× under the combined sweep and 3–17× under the half-t sweep (timings vary between runs; the ratios are within one run).
- **Defaults unchanged.** These are single data sets: exploration, not a study for defaults. A default changes only after one pre-registered comparison. Each draw conditions on the current values of the others; in particular, the $\tau_u^2$ statistic $\sum_r\|\beta_r-\beta_0\|^2$ uses the $\beta_0$ drawn in the same sweep.

`HierarchicalResult.PosteriorDraws` keeps every draw (`HierarchicalDraw`), and every summary is read from it:
- **What a draw holds:** the population curve, each replicate's curve, σ², τ_u² and τ₀², and the conditional means E[β₀ | ·] and E[β_r | ·] the sweep computes on the way.
- **The opt-out.** Replicate curves dominate memory, at R × ν numbers twice over per draw. `Run(retainReplicates: false)` keeps only the population curve and the variances; the replicate summaries then accumulate as each chain runs, to the same values.
- **Rao–Blackwellized means.** `RaoBlackwellizedPopulationCoefficients` and `RaoBlackwellizedReplicateCoefficients` average the conditional means, beside the current averages of the drawn curves. Both estimate the same posterior means. The Rao–Blackwellized ones drop the conditional draw's own noise; for a scan of more than two blocks that usually lowers the Monte Carlo error but is not guaranteed (Liu, Wong & Kong 1994 prove it for two). On the test data they lowered the across-seed variance of the coefficient estimates by about 16%, since most of the spread comes from how the variance components mix.
- **Intervals and bands,** read from the draws on request:
  - `NoiseVarianceInterval`, `DeviationVarianceInterval` and `PopulationSmoothingVarianceInterval` give each variance's median and 95% interval.
  - `PopulationBand(design)` and `ReplicateBand(r, design)` give pointwise bands, built the same way as BAPS's `FitBand`. A replicate band needs the replicate draws.

**Bivariate smoothing.**

- **Design.** `TensorProductDesign` builds the dense row-wise Kronecker design, with columns indexed $j\nu_y+k$.
- **Isotropic.** $P=P_x\otimes I+I\otimes P_y$ with one λ, via `PenalizedSpline` + `TensorPenalty`. `Baps` and `HierarchicalPSpline` sample it with the same penalty.
- **Anisotropic.** $A=Z^\top Z+\lambda_x(P_x\otimes I)+\lambda_y(I\otimes P_y)$. It evaluates REML from the 1-D penalty spectra: the combined penalty's eigenvalues are the pairwise sums $\lambda_xa_i+\lambda_yb_j$, so its pseudo-determinant is their product over the non-zero pairs.
- **The REML optimum.** `RemlOptimum` maximizes the anisotropic evidence over $(\ln\lambda_x,\ln\lambda_y)$.
  - **The search.** Nelder–Mead (the internal `NelderMead`, with Lagarias et al. 1998 coefficients) starts from each axis's balance point $\operatorname{tr}(Z^\top Z)/\operatorname{tr}(P_{\text{axis}})$ and searches within $10^{\pm8}$ of it.
  - **Edges.** The evidence levels off as an axis's λ grows when the data have no curvature along that axis, so an axis that ends at the box's edge is flagged (`LambdaXAtBound`, `LambdaYAtBound`). A trial point that is not numerically positive definite counts as infeasible.
  - **Checks** (`AnisotropicRemlTests`). On a surface gentle in x and fast in y it converged in 141 evaluations and beat every point of a 25 × 25 grid around it. On a surface constant in x, λ_x came out $1.6\cdot10^5$ times λ_y.
- **Cost.** The factor's half-bandwidth is $\max(\approx\deg_x\nu_y+\deg_y,\ r_x\nu_y,\ r_y)$, and the design term dominates for cubic bases. Each evaluation costs $O(\nu_x\nu_y\cdot\text{bw}^2)$.

**Factor checks (step 3.6).** BAPS's and BGP's systems are positive definite whenever their designs identify the penalty's null space: at least r distinct inputs for a difference penalty of order r. The hierarchical blocks' identity terms and BGP's nugget σ²I make theirs so unconditionally.
- **The rule.** A factor that is not finite, not positive definite, or clamped a pivot at the 1e-12 floor means degenerate input. The call throws `InvalidOperationException`, naming the matrix and its `FactorStatus`, because the floored solves would describe a different matrix.
- **Where it applies:**
  - `PenalizedSpline`'s solve and penalty spectrum;
  - both `Baps` Gibbs schemes;
  - both `HierarchicalPSpline` blocks;
  - `AnisotropicTensorSpline`'s REML factor;
  - `GpRegression.Fit`.
- **Why BARS differs.** BARS treats an inadmissible factor as outside the model's support (design admission, above), since a knot configuration is a proposal its chain can refuse.

