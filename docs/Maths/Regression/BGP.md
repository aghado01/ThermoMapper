# BGP — Bayesian Gaussian Process Bandwidth

> Gaussian process regression sampling kernel bandwidth under an empirical-Bayes prior, with known or adaptive unknown noise scale. Part of `Maths.Regression`.

### Bayesian GP bandwidth (BGP)

**Input units and the fitted transform.** `GpRegression` fits one `GpTransform` on the training inputs, and training and prediction share it:

- **Terms** (`BgpSampler.PointwiseLogDensity(result, units)`). The collapsed content only, since f is integrated and never drawn.
  - Per draw, log p(y_U | y_rest, t, η) with f, the constant mean and, for an unknown scale, s² integrated out (`GpRegression.LeaveOutLogDensities(fit, units)`).
  - The flat mean leaves the projected precision P = R⁻¹ − R⁻¹11ᵀR⁻¹/1ᵀR⁻¹1, with PY = α and YᵀPY = Q. A known σ² gives a Gaussian; an unknown scale under 1/s² gives a Student t over n − 1 effective observations.
  - A draw that repeats its predecessor's (t, η) reuses its row, so a run costs one O(n³) refit per accepted move.

- **Inputs.** The default, `GpInputScaling.Similarity`, is the isotropic map $\phi(x)=(x-\mathrm{lo})/s$. Here lo is the per-coordinate training minimum and s the largest coordinate span, so the data fill the unit cube with their shape preserved. A constant coordinate maps to 0. `GpInputScaling.None` uses the inputs as given, for callers whose coordinates already have squared distances of order 1.
- **Units.** The sampler, the prior's support and $T_n$ work in the mapped frame. A bandwidth t there is $s^2t$ in the original units. `BgpResult` reports the bandwidth posterior mean and its 2.5% and 97.5% quantiles in both frames, as `Bandwidth` and `BandwidthOriginal` (both `BandwidthSummary`), together with s as `InputScale`.
- **Retained draws.** `BgpResult.BandwidthDraws` keeps every draw, per chain in draw order (`BgpDraw`). A draw holds the bandwidth t, the GP fit the chain held there (its log marginal, α and mean), and that fit's conditional variance at the test points. Given t the GP posterior is fixed, so this is the draw's whole state. A rejected proposal shares the previous fit and its variance. The posterior mean, the posterior sd and the bandwidth summary are read from these draws.
- **Posterior sd.** `PosteriorSd` is the sd of f(x*) by the law of total variance: the spread of the draws' conditional means plus their mean conditional variance.
  - **The conditional variance** integrates the constant mean out under its flat prior (universal kriging; Rasmussen & Williams 2006, §2.7): $k(x^*,x^*)-k_*^\top K^{-1}k_*+(1-\mathbf 1^\top K^{-1}k_*)^2/\mathbf 1^\top K^{-1}\mathbf 1$. The last term is the mean's own uncertainty.
  - **It is the uncertainty of f,** not of a new observation, which adds σ².
  - **When it is computed.** The sampler computes it during the run, once per distinct fit, while that fit's factor of K exists. The factor costs n² numbers, so it is never retained. `GpRegression.PredictVariance(fit, x)` gives the variance at other points by refactoring at the fit's bandwidth.
  - **Checked** against an independent route (`BgpVarianceTests`): a proper N(0, 10⁸) prior on the mean, conditioned as one joint Gaussian and solved in `decimal`, agrees to 10⁻⁹, inside, at the edges and outside the data.
- **Checkpoints and the early stop.**
  - `Run(checkpoint: true)` returns the restart state. `Run(resumeFrom: …)` continues it exactly as the uninterrupted run would have.
  - The chain's state is its stream, its log bandwidth, the log target there, the bandwidth its current fit was made at (refitted on restore, giving the same fit), its counters and its retained draws.
  - The fingerprint covers the GP's data, noise variance, kernel and input map, the bandwidth prior, the test points, the starting bandwidth and the proposal scale.
  - `rHatTarget` stops the run once the bandwidth's split R̂ and bulk and tail ESS pass, as for the other shell clients.
- **Validation.** Non-finite inputs or responses are rejected, and so are all-identical inputs under `Similarity`.

**Noise model (2026-10-23).** `GpRegression.NoiseModel` says how the variances are treated.
- **`Known`, the default and TWCD2025's algorithm.** σ² is given and the kernel's signal variance is 1: K = Σ_t + σ²I. The response must then be on a scale where a unit signal variance is sensible.
- **`UnknownScale`** (`GpRegression.WithUnknownScale`). K = s²(Σ_t + ηI): the signal variance s² and the noise variance σ² = s²η are both unknown, and η is the noise-to-signal ratio.
  - **The marginal.** With the constant mean integrated out under its flat prior and s² under 1/s², $\log p(Y\mid t,\eta)=-\tfrac12\log|R|-\tfrac12\log\mathbf 1^\top R^{-1}\mathbf 1-\tfrac{n-1}{2}\log Q+\log\Gamma(\tfrac{n-1}{2})-\tfrac{n-1}{2}\log\pi$. Here R = Σ_t + ηI and Q = (Y − m̂1)ᵀR⁻¹(Y − m̂1).
  - **s² given (t, η):** s² | t, η, y ~ IG((n − 1)/2, Q/2). It is drawn exactly for every retained draw, on the readout stream, so the (t, η) chain does not depend on it.
  - **What s² changes.** The conditional mean m̂ + k*ᵀR⁻¹(Y − m̂1) does not depend on s². The conditional variance scales with it: each draw's `ConditionalVariance` is s² times the unit-scale variance (`PredictVariance`).
  - **The sampler.** It walks (log t, log η) jointly, one refactor per step.
    - **The proposal** is `BgpProposal.Adaptive` by default (adaptive Metropolis: Haario, Saksman & Tamminen 2001). In burn-in it learns a correlated step: its shape is (2.38²/2)·Σ̂ from the draws after the first tenth, updated every 50 steps, and its global scale takes a vanishing Robbins–Monro step toward an acceptance rate of 0.35, the two-dimensional optimum (Andrieu & Thoms 2008, algorithm 4; Roberts & Rosenthal 2001). The step is then frozen, so sampling is ordinary Metropolis–Hastings, and the frozen step travels with checkpoints.
    - **`BgpProposal.Fixed`** keeps independent steps (`proposalSd`, `noiseRatioProposalSd`, 0.4 each by default), which also start the adaptive one.
    - **The prior on η.** η has a log-uniform prior over [`noiseRatioLower`, `noiseRatioUpper`], 10⁻⁶ to 10² by default; its lower end also floors the nugget, which keeps R well conditioned. η is a second monitored functional (`NoiseRatioConvergence`), and `NoiseVarianceInterval`, `SignalVarianceInterval` and `NoiseRatioInterval` summarize the draws.
- **Alternatives considered:**
  - *σ² only, signal variance kept at 1:* a prior on σ² and a joint walk on (log t, log σ²). It is the closest match to the "+ unknown σ" row of TWCD2025's Lucky Cat table, which does not state its prior. It leaves the scale mismatch in place: a response far from unit scale still has a mis-scaled prior, and σ² absorbs the difference.
  - *Plug-in estimates* by maximizing the marginal likelihood: estimation, not inference. The paper lists empirical Bayes for σ as future work.
- **Trade-offs and implications.** Measured on single test data sets: exploration, not a study for defaults.
  - **Scale invariance.** Rescaling the response by c shifts the marginal by −(n − 1) log c, and rescales s² and σ² by c², with t and η unchanged.
  - **It fixes mis-scaling.** On a sine of amplitude 6 with noise sd 0.5, the known model, even given the true σ², shrank the bandwidth to compensate (t ≈ 0.025 against 0.13). Its fit MSE was 0.054 with 82% pointwise 95% coverage; the unknown-scale model reached 0.020 and 100%. On unit-scale data, where the known model is well specified, the two were alike (MSE 0.022 against 0.024; coverage 0.92 against 0.96).
  - **Mixing.** t and η trade off in the posterior: a longer bandwidth with more noise against a shorter one with less.
    - **With the fixed walk** this cost about tenfold. t's ESS fell from 1,567 under a known σ² to 174 (unit scale), and from 3,746 to 372 (amplitude 6), over 16,000 draws; η's ESS was 137–279.
    - **The adaptive proposal recovers most of it.** At equal work (4 chains × 4,000 draws, burn 1,000, two data sets × three seeds), t's bulk ESS was 1,695–2,098 against the fixed walk's 160–371, and η's was 1,593–2,130 against 129–269. R̂ was ≤ 1.003 against up to 1.025, with acceptance 0.32–0.37 and the same posterior.
    - **Both proposals target the right posterior.** They match a two-dimensional quadrature of the (log t, log η) posterior means: adaptive within 0.4 and 0.9 SE, fixed within 0.2 and 0.1 SE, with the fixed walk's SEs four times larger.
    - **A short burn-in leaves the adaptive step at the fixed one.** It needs some hundreds of burn-in steps to learn.
  - **s² is weakly identified:** median 87 with a 95% interval [16, 1263] on the amplitude-6 data. For squared-exponential kernels the signal variance and the bandwidth trade off, and only a combination of them is well determined (compare Zhang 2004 for Matérn variance and range). The fit is unaffected, but the bands can run conservative: the 100% coverage above is a sign of it.
  - **The η range matters at its ends.** A posterior piled against either bound means the range is too narrow.
  - **Outside the paper's theory, with an extension proved.** TWCD2025's contraction analysis assumes σ is known; the unknown-scale model is a full-Bayes extension it mentions but does not analyse. The extension is derived in [PosteriorContraction](../PosteriorContraction.md) and checked in Lean (BGP-4), conditional on the cited general contraction theorem and on TWCD2025's framework inequalities.
    - **The result:** with a proper, compactly supported prior on s (for example 1/s² restricted to a box) and σ₀ inside the induced σ range, the known-σ rate holds jointly for f and σ.
    - **The argument:** the authors' Remark A.1 rescaling carries the prior-mass bound to f*/s; their sieve and entropy bounds carry over with O(log n) extra terms; the tests for ‖f − f*‖ₙ + |σ − σ₀| are built from Gaussian tails. The paper's (60) must hold at an exponent above the prior-mass one by a σ₀-dependent margin, which its sieve can be tuned to give.
    - **Open:** a heavy-tailed amplitude prior, the posterior-mean analogue of the paper's Theorem 3.2, and a σ Bernstein–von Mises result (de Jonge & van Zanten 2013).
    - **The bridge to the improper 1/s²:** the posterior mass of s² outside such a box is measurable from a run's own draws.
- **Checked** (`BgpUnknownScaleTests`):
  - the closed-form marginal agrees with quadrature over s² against a dense reference to 10⁻⁸;
  - the rescaling shift is exact to 10⁻⁹;
  - the s² draws match their conditional, z = −0.08 against E[s² | ·] = Q/(n − 3);
  - resume is exact;
  - a known σ² reports s² = 1 and σ² as given.
- **Known σ² is unchanged:** it is bit-identical to before, down to each draw's log marginal.

The response is not rescaled, so σ² stays in response units.

**GP model.** $y=m\mathbf1+g+\varepsilon$, where g is a zero-mean GP with squared-exponential kernel $k_t(d^2)=\exp(-d^2/2t)$ on mapped distances, the noise σ² is known, and the constant mean m has a flat prior and is integrated out. With $K=K_t+\sigma^2I$, the evidence is

$$\log\int p(y\mid m,t)\,dm=-\tfrac12\Big[y^\top K^{-1}y-\frac{(\mathbf1^\top K^{-1}y)^2}{\mathbf1^\top K^{-1}\mathbf1}\Big]-\tfrac12\log|K|-\tfrac12\log\mathbf1^\top K^{-1}\mathbf1-\tfrac{n-1}2\log2\pi.$$

The term $-\tfrac12\log\mathbf1^\top K^{-1}\mathbf1$ depends on t, so it is kept. The GLS mean $\hat m=\mathbf1^\top K^{-1}y/\mathbf1^\top K^{-1}\mathbf1$ is the posterior mean of m. With $\alpha=K^{-1}(y-\hat m\mathbf1)$, predictions return the mean only, $\hat f(x_*)=\hat m+\sum_jk_t(\|\phi(x_*)-\phi(x_j)\|^2)\alpha_j$, and revert to $\hat m$ away from the data. A constant shift of y moves $\hat m$ and the predictions by the same amount and leaves the evidence unchanged. Each fit makes one extra solve, for $K^{-1}\mathbf1$, on the same factor. Internally, the training mean ȳ is subtracted before the solves for numerical conditioning only; the results do not depend on it.

**Empirical-Bayes prior** (Tang, Wu, Cheng & Dunson 2025):

$$\log p(t)=-a_0\ln t-\frac{b_0}{\hat v_n(t)},\qquad t\in(\gamma_1T_n^2,\,1],$$

- $\hat v_n(t)$ is the mean off-diagonal kernel value.
- $T_n$ is the mean distance to the K-th nearest neighbor, where the neighbor list counts the point itself first. So $K=2$ means the nearest *other* point.
- $K=2$ for $n<200$, and $\lceil\gamma_2(\ln n)^2\rceil$ otherwise.
- The defaults are $a_0=b_0=1$ and $\gamma_1=\gamma_2=0.25$.
- The upper bound 1 is fixed, which assumes squared distances of order 1 in the mapped frame; TWCD2025 take the domain to be $[0,1]^D$. The default similarity map provides that, but very few points in four or more dimensions can still give $\gamma_1T_n^2\ge1$.
- The support must be nonempty and finite in $\ln t$. Construction rejects fewer than two inputs, non-finite inputs or $T_n$, $\gamma_1T_n^2\ge1$ (empty support), and $T_n=0$ (every point's nearest others are duplicates, so $\ln t$ is unbounded below).

**Coordinate scaling and the theory.** A fixed nonsingular affine map of the inputs, including per-coordinate (anisotropic) scales, keeps the TWCD2025 rates, which depend only on the intrinsic dimension and the smoothness. This is a deduction from the theorem's hypotheses applied to the transformed manifold, not a new proof. The map changes only the constants: the curvature, reach, volume and diameter of the transformed manifold. An empirical transform fitted on the sample, such as the default similarity map or a per-coordinate min–max, depends on n, and a degenerating span needs further control; neither the paper nor a proof here covers that case. See BGP-2 in [the BGP lemma notes](../../../lean/Protolemmata/BARS/bgp-twcd-lemmas.md).

**Sampler.** A random-walk Metropolis–Hastings on $u=\ln t$ with target $\log p(y\mid e^u)+\log p(e^u)+u$. The last term is the Jacobian of $t=e^u$.

- The default proposal standard deviation is 0.4, with 2 chains, 200 burn-in steps and 600 draws.
- The chains start at the geometric midpoint $\sqrt{\gamma_1T_n^2}$ of the support; an explicit start outside the support is rejected.
- A proposal outside the prior's support is rejected without fitting.
- The posterior mean prediction averages $\hat f_t$ over all pooled draws.
- The bandwidth interval uses the type-7 2.5% and 97.5% quantiles of the pooled draws.
- Each step costs $O(n^3)$.

