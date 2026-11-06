# Probability and Information

> ThermoMapper's stochastic primitives: a checkpointable xoshiro256++ generator and a seed-expansion tree that make every Monte Carlo component reproducible from one integer, allocation-free samplers, quantiles and log-densities for the continuous distributions, truncated draws that hold far in a tail, a log-gamma function, and entropy and divergence measures. · **Projects:** `Maths.Rng`, `Maths.Distributions`, `Maths.Information` (namespaces of the same names) · **Maturity:** core; the `IDivergence<T>` abstraction is partial

## Overview

These three small projects sit at the bottom of the [architecture](./Architecture.md) (layers 0–1). They are consumed widely:

| Type | Consumers |
|---|---|
| `Xoshiro256PlusPlus` | distributions, regression, [samplers](./Samplers.md), SPC, geometry, synthetic noise fields |
| `SeedTree` | the ensemble shell, parallel tempering, BAPS, GMM K-sweeps, distributed SPRED |
| `Beta` | BARS knot-locality proposals (regression) |
| `Gamma` | regression (BAPS local-multiplier draws, BARS robust weights); the sampler behind `Beta`, `InverseGamma` and `StudentT` |
| `InverseGamma` | regression (the BARS σ² draw, BAPS and hierarchical P-spline variance components) |
| `Normal` | regression (standard draws behind coefficient and Laplace draws) |
| `Weibull` | regression (the Weibull smoothing-variance prior and its prior scaling) |
| `StudentT`, `Truncated<T>`, `Uniform`, `Mixture`, `Empirical` | BASS sensitivity's input measures, through their interval moments (regression, stage 6); calibration priors later |
| `HalfCauchy` | BASS calibration's heavy-tailed variance-group prior, C⁺(0, s) on σ_g (regression, stage 7) |
| `IIntervalMoments`, `GaussLegendre` | BASS sensitivity's moment service; `GaussLegendre` behind `StudentT`'s interval moments |
| `SpecialFunctions.LogGamma` | geometry, regression (changepoint and BARS priors) |
| `SpecialFunctions` incomplete gamma and beta | the distributions' CDFs |
| `Shannon` | SPC profiling, graph observables |

The KL, Jensen–Shannon and `IDivergence` types have no in-repository callers yet.

## Key concepts

- **The project die.** The library's standard generator is `Xoshiro256PlusPlus`. Every sampler takes a generator as an explicit argument, and none holds hidden global state. A few components elsewhere use `System.Random` instead (see [Concepts](./Concepts.md#reproducibility-contract)).
- **Seed tree.** A deterministic fan-out from one master seed to N child seeds, so a parallel ensemble of chains or workers can be reproduced exactly.
- **Checkpointable state.** The generator's four 64-bit state words can be saved and restored bit for bit. The SPC checkpoint format uses this.
- **Nats.** Every entropy and divergence here uses the natural logarithm.

## Design

- **Value types with named factories.** Each distribution is a readonly struct built by a factory that names its parameterization: `Normal.MeanSd`, `Gamma.ShapeRate` or `Gamma.ShapeScale`, `InverseGamma.ShapeScale`, `Weibull.ShapeScale`, `StudentT.LocationScale`, `Beta.Shapes`.
  - Conjugate updates produce rates, so a caller passes a rate as a rate rather than converting it at the call site. That conversion is where R's BASS τ prior slipped ([BASS derivations review §6](../../issues/regression/bass-derivations-and-integration-review.md#6-d03-finding-rs-default-τ-hyperprior)).
  - The structs allocate nothing and take the generator explicitly. `Normal.SampleStandard` and `Gamma.SampleStandard` are the unit draws.
  - `SpecialFunctions`, `KLDivergence`, `JSDistance` and `Shannon` stay static.
  - Reproducibility is purely a seeding question: one generator per chain, seeded from the seed tree.
- **Log-space tails.** CDFs, survival functions and interval probabilities come from whichever tail represents them directly, in log space, so a probability far below the double range keeps its digits as a log. `Truncated<T>` normalizes by its parent's interval probability.
- **A log-concave coordinate per family.** Each density is log-concave in one coordinate at every parameter value, and its truncated draw and its quantile both work there, with no branch on the parameter ([below](#truncated-draws-and-quantiles)).
- **Seeding convention.** `SeedTree.Derive(master, count)` has the prefix property: `Derive(m, C)` is a prefix of `Derive(m, 2C)`. The ensemble shell uses this to give chain c a *kernel* seed from children $[0, C)$ and a separate *readout* seed from children $[C, 2C)$. See [Samplers](./Samplers.md).
- **Dispatch surface for divergences (partial).** `IDivergence<T>` is a dispatch interface over a domain T. It has one adapter, `KlDivergence`, over `ReadOnlyMemory<double>`, and no consumers. The Jensen–Shannon distance mentioned in its documentation does not implement it.

## API at a glance

| Area | Entry points |
|---|---|
| Generator | `Xoshiro256PlusPlus(int? seed)`, `NextUInt64`, `NextDouble`, `NextInt(q)`, `SaveState`, `LoadState`, `Restore` |
| Seeds | `SeedTree.Derive(master, count)` |
| Distributions | `Normal.Standard`, `Normal.MeanSd(mean, sd)`, `Normal.MeanVariance`, `Normal.SampleStandard(rng)`, `Normal.StandardQuantile(p)` (Wichura's AS 241, about 1 part in 10¹⁶; checked against the log-space CDF to p = 10⁻³⁰⁰), `Normal.StandardInverseLogCdf(logP)`; `Gamma.ShapeRate(shape, rate)`, `Gamma.ShapeScale(shape, scale)`, `Gamma.SampleStandard(rng, shape)`, `Gamma.SampleLog(rng)`; `InverseGamma.ShapeScale(shape, scale)`; `Weibull.ShapeScale(shape, scale)` with `Median` and `Mode`; `StudentT.Standard(ν)`, `StudentT.LocationScale(ν, location, scale)`; `Beta.Shapes(a, b)`; `HalfCauchy.WithScale(s)`; `Uniform.On(lower, upper)`; `Mixture.Of(components, weights)`; `Truncated.Of(parent, lower, upper)`. Each is an `IContinuousDistribution`: `Sample(rng)`, `SampleOn(rng, lower, upper)`, `LogPdf`, `Cdf`, `Survival`, `LogCdf`, `LogSurvival`, `IntervalProbability`, `LogIntervalProbability`, `Quantile(p)`, `InverseLogCdf(logP)`, `InverseLogSurvival(logS)`, and its moments where it has them. `Empirical.Of(sample)`; `IIntervalMoments.IntervalMoments(lower, upper, centre, logScale, moments)`, implemented by `Normal`, `StudentT`, `Truncated<T>`, `Uniform`, `Mixture` and `Empirical`; `GaussLegendre.Of(n)` |
| Special functions | `LogGamma`, `LogBeta`; `RegularizedGammaP`, `RegularizedGammaQ` and their logs; `RegularizedBeta`, `LogRegularizedBeta`, `LogRegularizedBetaComplement` |
| Information | `KLDivergence.Forward(p, q, eps)`, `KLDivergence.Symmetric(p, q, eps)`; `IDivergence<T>`, `KlDivergence(eps)`; `JSDistance.Distance(p, q)`; `Shannon.EntropyNats` (counts or probabilities) |

## Algorithms

### xoshiro256++ and seeding

**Output and state update.** The state is four `ulong` words $s_0..s_3$. Each call returns $\text{out}=\operatorname{rotl}_{23}(s_0+s_3)+s_0$ and then advances the state:

$$t=s_1\ll17;\quad s_2\mathrel{\oplus}=s_0;\ s_3\mathrel{\oplus}=s_1;\ s_1\mathrel{\oplus}=s_2;\ s_0\mathrel{\oplus}=s_3;\ s_2\mathrel{\oplus}=t;\ s_3=\operatorname{rotl}_{45}(s_3).$$

This is the published xoshiro256++ reference algorithm.

**Conversions.**

- `NextDouble` is $\lfloor x/2^{11}\rfloor\cdot2^{-53}\in[0,1)$.
- `NextInt(q)` is a multiply-shift, $\lfloor(\lfloor x/2^{32}\rfloor\cdot q)/2^{32}\rfloor$, with no rejection step, so a relative bias of about $q/2^{32}$ remains. It throws for $q\le0$.
- Consumers convert uniforms to $(0,1]$ as `1.0 - NextDouble()` before taking a logarithm.

**Seeding.**

- An `int` seed is expanded by four successive SplitMix64 calls, one per state word. Negative seeds sign-extend.
- A `null` seed draws 32 bytes from the OS `RandomNumberGenerator`, which is **not reproducible**.
- If seeding yields the all-zero state, $s_0$ is set to the golden-ratio constant `0x9e3779b97f4a7c15`; the restore paths reject that state instead.

**Checkpointing.** `SaveState` returns the four words. `LoadState` restores them in place, and `Restore` builds a new generator from them. Both reject the all-zero state.

**`SeedTree.Derive`.** Child i is the low 32 bits of the i-th SplitMix64 output, starting from state `(ulong)master`, using the same constants as the generator. Each child gets a distinct *index*. The 32-bit child *values* are not guaranteed to be distinct, so collisions follow the usual birthday bound.

### Distributions

- **Normal.** Box–Muller, cosine branch only: $z=\sqrt{-2\ln u_1}\cos(2\pi u_2)$ with $u_i\in(0,1]$. Each draw consumes two uniforms.
- **Gamma(α, rate β)**, scale θ = 1/β, uses Marsaglia–Tsang (2000):
  1. Set $d=a-\tfrac13$ and $c=1/\sqrt{9d}$.
  2. Draw $x\sim N(0,1)$ until $1+cx>0$, then set $v=(1+cx)^3$.
  3. Accept $d\,v\,\theta$ if the squeeze $u<1-0.0331x^4$ holds, or if $\ln u<\tfrac12x^2+d(1-v+\ln v)$.

  For shape $a<1$ it uses the boost $G(a)=G(a+1)\,u^{1/a}$, which recurses at most once. The log-density is $(\alpha-1)\ln x-\beta x+\alpha\ln\beta-\ln\Gamma(\alpha)$, and $-\infty$ for $x\le0$. The CDF is $P(\alpha,\beta x)$.
- **Beta(a, b)** is sampled as $G_a/(G_a+G_b)$ with unit-scale gammas. If the sum underflows it returns 0.5. The log-density is $(a-1)\ln x+(b-1)\ln(1-x)-\ln B(a,b)$ on $(0,1)$.
- **Inverse gamma IG(α, β)** is $1/G$ with $G\sim\mathrm{Gamma}(\alpha,\text{rate }\beta)$; β is a scale of X and the rate of 1/X. The draw uses the same arithmetic as the call sites it replaced, so seeded streams are unchanged. The CDF is $Q(\alpha,\beta/x)$.
- **Student-t** draws $\mu+s\,Z\sqrt{\nu/V}$ with $V\sim\chi^2_\nu=2\,\mathrm{Gamma}(\nu/2)$. Its tail is $P(T>t)=\tfrac12 I_x(\nu/2,\tfrac12)$ with $x=\nu/(\nu+t^2)$. Both $x$ and $1-x$ are formed from whichever of $t^2/\nu$ and $\nu/t^2$ is below 1, with their logs, so neither overflows.
- **LogGamma** uses the Lanczos approximation with $g=7$ and nine coefficients for $x\ge 0.5$, and the reflection $\ln(\pi/\sin\pi x)-\ln\Gamma(1-x)$ below that. The documented domain is $x>0$.
- **Incomplete gamma.** $P(a,x)$ comes from the power series for $x<a+1$, and $Q(a,x)$ from the Lentz continued fraction otherwise (Numerical Recipes §6.2). The prefactor $x^a e^{-x}/\Gamma(a)$ is kept in log space. The other tail is the accurate complement $\ln(1-e^u)$ (Mächler 2012). Each loop runs at most 100,000 iterations and throws `ArithmeticException` if it does not converge.
- **Incomplete beta.** $I_x(a,b)$ comes from the Lentz continued fraction (Numerical Recipes §6.4), switching to $I_{1-x}(b,a)$ above $x=(a+1)/(a+b+2)$.
- **Normal tails.** $\Phi(-z)=\tfrac12 Q(\tfrac12, z^2/2)$.
- **Interval probabilities.**
  - An interval above the median is $S(a)-S(b)$ and one below is $F(b)-F(a)$, each computed in log space, so neither subtracts from 1.
  - An interval across the median is $1-F(a)-S(b)$ in general. The normal uses $\tfrac12[\operatorname{erf}(|\alpha|/\sqrt2)+\operatorname{erf}(\beta/\sqrt2)]$ instead, and the Student-t $\tfrac12[P(|T|<|\alpha|)+P(|T|<\beta)]$. Both stay accurate for an interval narrow around the centre.
- **Interval moments** (`IIntervalMoments`, 2026-10-29, for BASS's D08): e^(−logScale)·∫ (x − c)ᵏ dF over (lower, upper], with an error estimate.
  - *The interval is half-open,* as the interval probabilities are, so consecutive intervals count an atom once.
  - *The log scale* lets a truncation pass its log mass down, so moments far in a tail keep their digits.
  - *A centre inside the interval* keeps the powers from cancelling.
  - *`Uniform`* has the closed form; *`Empirical`* takes compensated sums over its atoms; *`Mixture`* sums its components' moments by weight.
  - *`Normal`* integrates the density's Taylor series on stretches with h·(1 + |z|) ≤ 1, summed to convergence: analytic, with no quadrature error. The Φ/φ recurrence was rejected because on a narrow piece it loses a factor of about 1/h² per power, which costs every digit by degree 2–3 for knots 0.0015 apart under σ = 0.2. Stretches below e⁻⁹² of the interval's peak density are dropped, and an interval beyond about 10⁸ sd is refused.
  - *`StudentT`* uses Gauss–Legendre at every ν, on a finite interval (truncated first), with panels graded geometrically from the location so that each stays three half-widths or more from the poles at location ± i·scale·√ν. It returns the 40-node values, and their difference from the 20-node ones is the error estimate.
  - *`Truncated<T>`* clips the interval and passes its log mass on as the scale. `Truncated.Of(mixture)` is R's mixture convention (weights before truncation); a mixture of truncated components is the conditional one.
- **Truncation.** `Truncated<T>` has density $f(x)/P(\text{lower}\le X\le\text{upper})$ on the interval, with the normaliser taken from the parent's log interval probability. N(0, 1) on $[10,\infty)$ therefore has density $\varphi(x)/Q(10)$, where $1-\Phi(10)$ is exactly 0 in double precision.
- **Accuracy**, checked against mpmath at 50–60 digits (`DistributionReferenceTests`):
  - **LogGamma** is within $2\times10^{-14}\max(1,|\ln\Gamma|)$ at ten points from $10^{-8}$ to $10^{5}$.
  - **Incomplete gamma and beta** are limited by their log prefactor, a sum of terms of size $S$ ($|\ln\Gamma(a)|$, $|a\ln x|$, …). Their relative error is about $\varepsilon S$: about $10^{-13}$ for moderate shapes and $10^{-12}$ at shape 1000.
  - **Far tails** keep their digits in log space, for example $\ln\Phi(-40)=-804.6084\ldots$ and $\ln P(T>10^{200})$ for ν = 1.
  - **A same-tail interval** loses relative accuracy as it narrows: about $\varepsilon\max(1,|\alpha|)/(\beta-\alpha)$ for the normal.

### Truncated draws and quantiles

Each family's density is log-concave in one coordinate at every parameter value, and both services work there:

| Family | Coordinate y | Density in y | Peak |
|---|---|---|---|
| Normal | $z=(x-\mu)/\sigma$ | $e^{-z^2/2}$ | $x=\mu$ |
| Gamma(α, rate β) | $\ln(\beta x)$ | $e^{\alpha y-e^{y}}$ | $x=\alpha/\beta$ |
| Inverse gamma(α, β) | $\ln(\beta/x)$ | $e^{\alpha y-e^{y}}$ | $x=\beta/\alpha$ |
| Beta(a, b) | $\operatorname{logit}x$ | $e^{ay}(1+e^{y})^{-a-b}$ | $x=a/(a+b)$ |
| Student-t(ν) | $\operatorname{asinh}(t/\sqrt\nu)$ | $\cosh(y)^{-\nu}$ | $t=0$ |

In x itself a gamma below shape 1 is log-convex, an inverse gamma's upper tail is log-convex and a t's tails are polynomial. Here none of them needs a branch on its parameter.

**Truncated draws.** `SampleOn(rng, lower, upper)` on every family; `Truncated<T>.Sample(rng)` is its parent's on its bounds. Each draws from the one `Xoshiro256PlusPlus` stream, deterministically.

- *The hat* (`LogConcaveEnvelope`). The anchor is the interval's point closest to the peak, and $h(\delta)=\ln f(\text{anchor}+\delta)-\ln f(\text{anchor})$ is concave with maximum 0. Each side gets two tangents: at the anchor, flat if the peak is interior, and where h first reaches −1, or at the interval's end if it never does. A proposal is an exponential cut to one tangent's span, three uniforms each, rejected against h. This is adaptive rejection's hat (Gilks and Wild 1992, from memory) with fixed design points, as in transformed density rejection (Hörmann 1995).
- *Where it holds.*
  - *The far tail:* h is a difference, $-\delta(z+\delta/2)$ for the normal, with no normalizer or CDF in it. N(0, 1) on $[10,\infty)$ accepts 99.6%, and on $[1000,\infty)$ runs the same arithmetic.
  - *A narrow interval:* a draw is the anchor plus a small offset, so it keeps the anchor's digits, and the two tangents nearly coincide: 100% on $[10,10.001]$.
  - *Far from the mode:* Gamma(1000, 1) on $[0,1]$ and IG(50, 10⁻³) on $[1,\infty)$ are nearly their power laws, $x^{\alpha-1}$ and $x^{-\alpha-1}$, which the far tangent matches.
- *The bound.* With $h(t)=-1$, a flat piece to $k$ and an exponential of rate $r$ beyond it have area $k+1/r=t$. By concavity $h\ge-x/t$ on $[0,t]$, so the density holds at least $t(1-1/e)$, and at least $1-1/e\approx63\%$ of proposals are accepted; an interval's end cuts both sides, and the tangent at the anchor only tightens the hat. Newton stops at $|h+1|\le\tfrac14$, which proves $1-e^{-3/4}\approx53\%$. Measured over a grid of intervals, the floor is 88.6% for the normal, 86.8% for the gamma, 88.6% for the t and 88.7% for the beta.
- *The guard.* Newton from a reach far past the point where h falls by 1 steps back about one unit of h's growth at a time, so it can run out of steps where h is astronomically low, and a far tangent there is so steep that rounding in its cut swamps the hat or sinks it below h. Where Newton's search ends short of its point, a guarded search runs from the reach: Newton where h ≥ −2²⁰, and below it Newton on ln(−h) against ln t, which lands at t·|h|^(−1/k) with k = t·h′/h ≥ 1 by concavity, exact for a power of t and close for an exponential, or a halving in ln t where h has overflowed. Every search that converged is unchanged bit for bit, so any reach is valid and a good one only saves steps. *Checked* (`LogConcaveEnvelopeTests`): every family's own log ratio, 30 cases each over far ranges, drawn with reaches from the root to 10¹⁰⁰ times it and +∞, against its law by quadrature. Before, acceptance fell to 0, draws threw, and Beta's and t's failed Kolmogorov–Smirnov at 10¹⁰⁰ with the hat below h; now every law holds, with acceptance 0.867 or more. A pin over 1,311 cases of the families' draws moved only 31, all the modified half-normal's whose own reach ran Newton out of steps; the GIG's truncated draw on far upper tails at small |p|, a and b, which threw in 4 of 4,000 cases, now draws its law.
- *Ratio-of-uniforms* (`RatioOfUniforms`), the untruncated draw where a family has it: (U, V) uniform on [u₋, u₊] × (0, 1], X = U/V + m, accepted when V ≤ √(f(X)/f(m)) (Kinderman and Monahan 1977; the shift by the mode m, Dagpunar 1989 and Lehner 1989). With u₋ and u₊ the extremes of (x − m)√(f(x)/f(m)) the rectangle is minimal; for a T₋½-concave f the region is convex and fills at least half of it. A family passes the sides in closed form, as the GIG does, or takes them from `RatioOfUniforms.Bounds`, golden section on each side. *Checked* (`GigRatioOfUniformsTests`): N(0, 1)'s sides found at ±√(2/e) to 7 digits and its draws, and Gamma(3.5)'s by search, on their laws; the GIG's three regions on 18 chosen cases and 300 seeded far-range cases (|p| to 10³, a and b over 16 decades), worst √n·D 1.76 against the cut 2.6.
- *Closed forms.* The Weibull's $(x/\lambda)^k$ is a unit exponential and the half-Cauchy is uniform in angle, so each inverts on the interval with one uniform and no rejection. `Mixture` picks a component by weight and mass on the interval, which reads its tails.
- *The shape.* `IContinuousDistribution` gains `Sample` and `SampleOn`, so `Truncated<T>` dispatches by a constrained call with no boxing and a family cannot lack a draw.
  - A separate sampler interface tested at run time, as `IIntervalMoments` is, would only be needed if some family could not draw. All can.
  - Static per-family classes would hide the draw from `Truncated<T>`.
  - `SampleOn` reads no log mass, so a loop whose bounds or means move calls the parent, as `Normal.MeanSd(μᵢ, 1).SampleOn(rng, 0, ∞)` for probit, and skips `Truncated.Of`, whose normalizer costs an incomplete gamma.

The alternatives, with acceptance per proposal and rough nanoseconds per draw (the 5th percentile of short rounds on one shared machine, not a benchmark):

| Sampler | Acceptance | ns per draw | Reach |
|---|---|---|---|
| **The tangent hat (chosen)** | 88.5% on $[0,\infty)$, 99.6% on $[10,\infty)$, 100% on $[10,10.1]$ | normal: 60–115 on a half-line, 150–430 finite. Gamma and inverse gamma: 300–800. Beta and t: 500–950 | Every family and interval |
| Robert (1995): a translated exponential at his optimal rate, a uniform on a short interval, normal rejection across a wide one | 76% on $[0,\infty)$, 99.5% on $[10,\infty)$, 63% on $[10,10.1]$ | 30–45 on a half-line, 20–115 finite | The normal only |
| Botev's `trandn` in one dimension: a Rayleigh tail beyond 0.66, after Marsaglia (1964), normal rejection across a wide interval, inversion across a short central one | 31% on $[0.5,\infty)$, 99% on $[10,\infty)$ | 13–30 in the tail, 300–640 central with this module's CDF | The normal only; the central inversion is a difference of CDFs |
| Minimax tilting (Botev 2017) | Not measured | Not measured | Multivariate. In one dimension there is no dependence to tilt away, and it reduces to the row above, so it belongs to correlated latents (E22), not here |
| Inversion through the log tails, `Truncated<T>.Quantile(U)` | No rejection | 2–9× the hat | Every family. A narrow interval is a small difference of log tails, which resolve only to ε·\|ln S\| |
| Dagpunar (1978); Philippe (1997), for a truncated gamma | Not measured | Not measured | Part of the parameter space: shape above 1; integer shape or a beta mixture |
| Rejection from the parent | The interval's mass | | Fails in a tail: $e^{-50}$ on $[10,\infty)$ |

*The choice.* One sampler carries every family: one bound, one proof, one set of checks, and no regime to tune. *What it costs:* about three of Robert's draws on a normal half-line, probit's case (1.4 at the mean, 2.5–9 on a finite interval), for more setup per call. Robert's acceptance is below the hat's everywhere but across a wide interval. If E20's latent block ever profiles here, a closed half-line path is an addition behind `Normal.SampleOn` that changes no caller; nothing has measured that need.

**Quantiles.** `Quantile(p)`, `InverseLogCdf(logP)` and `InverseLogSurvival(logS)` on every family. The log forms exist because a lower tail below 10⁻³⁰⁸ and an upper tail below 10⁻¹⁶ have no p.

- *The solver* (`TailInversion`). Newton on the smaller log tail, in the table's coordinate. A log-concave density has log-concave tails (Prékopa 1973, from memory), so the log tail is concave in y and a step from one side of the root stays on that side, with bisection inside a bracket when a step leaves it. A step below 10⁻⁸ ends the solve and comes back as a correction, which is second-order small, so no evaluation goes on the last digits.
- *Starts.* The tail's asymptote in y: $g^\alpha/\Gamma(\alpha+1)$ for the gamma, $x^a/(aB)$ for the beta and the t's straight line. It bounds the root from the side Newton keeps, and is the answer itself where x underflows. In the bulk, Wilson–Hilferty for the gamma, the logit's normal approximation for the beta and Cornish–Fisher for the t. A start only sets the number of evaluations.
- *The normal* keeps AS 241 where $e^{\text{logP}}$ is in range, and is Newton below that.
- *The t near its centre* solves on $P(|T|<t)$, so `Quantile(½ + 10⁻¹²)` keeps relative digits.
- *The beta* solves in the logit and forms x and 1 − x apart. An x within 2⁻⁵³ of 1 reads as 1, and `Beta.Shapes(b, a).InverseLogCdf(logS)` is the complement at full precision.
- *`Truncated<T>`* reads its parent's: $F(x)=F(\text{lower})+p\cdot\text{mass}$ or $S(x)=S(\text{upper})+(1-p)\cdot\text{mass}$, sums in log space that cancel nothing.
- *`Mixture`* solves between its components' quantiles, where the sum must cross, and bisects across a gap in its support.
- *The alternative.* Approximations tuned for one family, as in R (AS 91 for `qgamma`, AS 109 for `qbeta`, Hill 1970 for `qt`, all from memory), are faster per call and stop short in logs, as below. This solver's accuracy is the tails', and the tails are already checked. It costs 0.3–5 µs per truncated inversion, which suits everything but a hot loop, and there `SampleOn` reads no tail.
- *Checks that do not rest on this code.* Closed forms far into both tails: the exponential, IG(1, β), Beta(a, 1), Beta(1, b), the arcsine, the Cauchy, t(2), χ²₁ and Gamma(3). Round trips through this module's log tails. And base R 4.6.1 over 510 cases: 418 hold their tail through R's `p*` within $10^{-12}\max(1,|\text{tail}|)+2\times10^{-14}S$, and the rest are an end of the double's range or an x rounded onto 1. Against R's own `q*`, 361 gaps, the worst 1.4×10⁻¹⁰, at ν = 10⁵, where this side's ln B(ν/2, ½) cancels.
- *A sweep at parameters nobody chose* (`DistributionSweepTests`). Seeded, and drawn over many orders of magnitude: 600 truncated cases by Kolmogorov–Smirnov against the tails, which the draws do not read, and 3,000 cases each of quantile round trips and of the density against a thin slice of the tail. Its first pass found what the chosen cases had not, each now corrected: the half-Cauchy's draws collapsing onto a bound past $10^{16}$ scales, where $\pi/2-\arctan a$ cancels and $ab$ overflows; its log density reading $-\infty$ past $10^{154}$ scales; the Weibull's lower tail and its quantile underflowing with the root still in range; and a mixture returning a finite number for a root beyond the double range. That last one now returns the range's end, as `StudentT` does.
- *A tiny shape's mass at the doubles' end.* Under a Gamma of shape a, P(X < x) is (Rate·x)ᵃ/Γ(a + 1) near 0, so at a = 10⁻³ about half of X lies below the doubles. `SampleLog` gives ln X on `Sample`'s stream, the boost's ln Y + (ln U)/a, which stands where X would underflow. `LogCdf` and `LogSurvival` read that form in logs where Rate·x falls under 10⁻³⁰⁰, where the product had underflowed and they read −∞ and 0 with x still in range. Both checks are in `DistributionSweepTests`, `SampleLog` by Kolmogorov–Smirnov over shapes from 10⁻⁶ to 10³.
- *Where R falls short, corrected rather than ported:*
  - `qt` stops short in the far tail: 1% in t at ν = 1.5 and tail −700, where R's own tail reads −700.015, and 8×10⁻⁹ at ν = 3.
  - `qt` overflows to ∞ at ν = ½ and tail −230, where t = 1.03×10¹⁹⁹ holds it.
  - `qbeta` returns a floor of 2⁻¹⁰²³ where the root lies below the double range, as at Beta(0.01, 0.01) and tail −23, whose tail at that floor is −7.8. This side returns 0.

### Normal variance-mean mixtures: the mixing law (E1)

A noise that is Gaussian given a latent per observation, εᵢ | vᵢ ~ N(m·vᵢ, s²·vᵢ), is a normal variance-mean mixture, and its mixing law is v's (`IMixingLaw`, in `Maths.Regression`, since BASS and BARS share it; [the BASS use](./VarianceMeanMixtures.md)). It supplies v's prior draw, its conditional given a residual r at a likelihood power b, ∝ p(v)·N(r | m·v, s²·v)^b, the marginal ln ∫ N(r | m·v, s²·v)p(v)dv, and E[v] and Var[v], +∞ where infinite.
- *Student t* (`StudentTMixing`): v ~ IG(ν/2, ν/2), GIG(−ν/2, 0, ν) in `GigParameters`' form, so at m = 0 the noise is s·t_ν. Its conditional is the inverse gamma IG((ν + b)/2, (ν + b·r²/s²)/2), drawn as the reciprocal of the precision's Gamma; its marginal is the t density. ln Γ((ν + 1)/2) − ln Γ(ν/2) is taken from the two log gammas below ν/2 = 10³ and from Stirling's series term by term above, ½ln x + x(ln(1 + h) − h) + (1/(12(x + ½)) − 1/(12x)) − (1/(360(x + ½)³) − 1/(360x³)), x = ν/2, h = 1/(2x), whose dropped terms are under 10⁻¹⁸ there; and ln(1 + r²/(νs²)) by `SpecialFunctions.Log1p`, public for it. The direct forms, as BARS's `StudentTWeights` takes them, read up to 2.0·10⁻⁵ off, relative, at ν = 3·10¹⁰ on the sweep below.
- *Checked* (`BassStudentTTests`): the marginal against quadrature of the mixture over ln v, on 11 chosen cases and a seeded sweep of 400 with ν from 0.05 to 10¹¹ and s² over 24 decades, worst relative gap 1.4·10⁻¹³. The quadrature's own inverse-gamma normalizer needed Stirling's form above ν/2 = 10³, or it read 2.5·10⁻¹⁰ low at ν = 10⁶, the cancellation the law avoids. The conditional's draws against the Gamma's mean and log mean at b = 1, 0.4 and 0, within 5 standard errors.
- *Members to come:* the asymmetric Laplace's exponential mixing and the Normal-Wald's inverse Gaussian, which draw from a GIG with a > 0 (the GIG sampler, E1's slices 2 and 3), and generalized beta prime mixing, which is not conjugate (slice 4). A drift m ≠ 0 is refused until then.

### The GIG, the modified half-normal and the asymmetric Laplace (E1)

Three laws for normal variance-mean mixtures, built with E1's slice 2 ([the BASS use](./VarianceMeanMixtures.md#the-asymmetric-laplace-quantile-regression-e1s-slice-2)).
- *`GeneralizedInverseGaussian`*, GIG(p, a, b) ∝ x^(p−1)·exp(−(a·x + b/x)/2), its algebra and limits in `GigParameters`. **Draws.** In u = ln x the log density is p·u − ω·cosh(u − ln s), ω = √(ab), s = √(b/a), whose second derivative is negative: every GIG is log-concave in ln x. So the tangent hat behind every truncated draw ([above](#truncated-draws-and-quantiles)) samples it exactly for every parameter, far ranges and an interval included, about the peak u* = ln s + asinh(p/ω). The log ratio there is −A·g(δ) − B·g(−δ), g(x) = eˣ − 1 − x, A = (R + p)/2, B = (R − p)/2, R = √(p² + ω²), the smaller of A and B taken as ω²/(2(R + |p|)): two non-positive terms, so nothing cancels when ω ≪ |p|, where p·(δ − sinh δ) − 2R·sinh²(δ/2), the first form tried, lost its digits. The hat's reaches start where each term alone falls by 1: before the envelope's guard a reach far past that point put its far tangent where h is −4·10¹⁶, and the rounding in that tangent's cut, some 500 nats, left the hat with no acceptance; now they only save it a search. The limits b = 0 and a = 0 draw the gamma and the inverse gamma. An untruncated draw is Hörmann and Leydold's (2014), GIGrvg's, where `DrawsByRatio` takes it: in the standard form y ~ GIG(λ, ω, ω), λ = |p| by y ↦ 1/y and x = √(b/a)·y, ratio-of-uniforms with a shift by the mode where λ > 2 or ω > 3, the rectangle's sides at a cubic's roots by Cardano's rule; without the shift where the density is T₋½-concave, λ ≥ 1 − 9ω²/4, or ω > 0.2; and below both a hat of three pieces, constant to ω/(1 − λ), then a power, then an exponential. It needs no search for a design point and is 2.5 to 5.6 times faster a draw on the parameters fits produce; the hat keeps truncated draws and the range past |p| ≤ 10³, 10⁻⁸ ≤ ω ≤ 10⁸ ([the measurement](./VarianceMeanMixtures.md#the-gigs-draw-ratio-of-uniforms-where-it-is-faster)). Its normalizer and moments are [below](#log-bessel-k-the-gigs-normalizer-and-the-generalized-hyperbolic-e1).
- *`ModifiedHalfNormal`*, MHN(α, β, γ) ∝ x^(α−1)·exp(−βx² + γx) on x > 0 (Sun, Kong and Pal 2021), the law of σ⁻¹ for a mixture's global scale under a drift. Three exact draws: for α ≥ 1 it is log-concave in x, and the hat's ratio about the mode is (α − 1)(ln(1 + t) − t) − βδ², t = δ/x*, the mode's identity cancelling the linear terms; for α < 1 with γ ≤ 0 it is log-concave in ln x, with the ratio βz*²(2δ − (e^(2δ) − 1)) + γz*((e^δ − 1) − δ); for α < 1 with γ > 0 neither, and e^(γx) expanded gives a mixture, k ∝ (γ/√β)ᵏΓ((α + k)/2)/k!, then x² ~ Gamma((α + k)/2, β). The paper's own samplers are rejection schemes per branch.
- *`AsymmetricLaplace`*, ALD(μ, σ, q) ∝ exp(−ρ_q((x − μ)/σ)) with the check loss ρ_q, the quantile likelihood whose q-quantile is μ (Yu and Moyeed 2001). Two exponential tails meet at μ, so its density, tails, quantiles and draws, on an interval too, are closed forms in logs. It is the normal mixture μ + σ(θv + τ√v·ζ), v ~ Exp(1), θ = (1 − 2q)/(q(1 − q)), τ² = 2/(q(1 − q)).
- *Exponential mixing* (`ExponentialMixing`, the [mixing law](#normal-variance-mean-mixtures-the-mixing-law-e1) of the asymmetric Laplace): v's conditional is GIG(1 − b/2, 2 + b·m²/s², b·r²/s²), the gamma limit at r = 0, and its marginal is closed at any drift, (m² + 2s²)^(−½)·exp(m·r/s² − |r|·√(m² + 2s²)/s²).

*Checked* (`BassQuantileTests`, `GbassPrimitiveParityTests`). The GIG on 134 cases, 14 chosen and 120 swept over far ranges, 4,000 draws each, worst decile z 3.01 against its law by quadrature, and on four intervals; the modified half-normal on 134 cases across its three branches, worst decile z 3.32; against `gbass_primitives_oracle.R`'s law, 11 GIG and 9 modified half-normal cases at 50,000 draws, worst |z| 2.56 and 2.55. The asymmetric Laplace against the exponential mixture to 12 digits, its tails and quantiles round-tripped, its draws on and off intervals; exponential mixing's marginal against quadrature over 305 cases, worst relative gap 8.7·10⁻¹⁵. The sweep's first pass found the cancellation and the reach above, and the half-normal's reach landing on x = 0, each corrected.

### Log Bessel K, the GIG's normalizer and the generalized hyperbolic (E1)

Built with E1's slice 3 ([the BASS use](./VarianceMeanMixtures.md#the-normal-wald-the-skewed-t-and-the-generalized-hyperbolic-e1s-slice-3)).
- *`SpecialFunctions.LogBesselK`, `LogBesselKScaled`*: ln K_ν(x) and ln(eˣK_ν(x)) for real ν and x ≥ 0, in logs over far ranges. For the fractional order μ = ν − n, |μ| ≤ ½, Temme's series for x ≤ 2 and Steed's continued fraction CF2 beyond (Temme 1975; Thompson and Barnett 1986; Numerical Recipes §6.7); then the upward recurrence to ν, stable for K, carried in the scaled ratio s = x·K_{μ+k+1}/K_{μ+k} = 2(μ + k) + x²/s, which neither overflows nor underflows; at ν ≥ 300 Debye's uniform expansion to four terms, νη written as ν(asinh(1/z) − 1/(z + √(1 + z²))) so neither end of z = x/ν cancels. *Alternatives:* Amos's algorithm, R's, in linear space with scaling; the integral representation by quadrature, which is the check here.
- *`GeneralizedInverseGaussian.LogNormalizer`, `LogPdf`, `Mean`, `Variance`*: ln ∫ x^(p−1)e^(−(ax + b/x)/2) dx = ln 2 + (p/2)·ln(b/a) + ln K_p(√(ab)), closed at the gamma and inverse-gamma limits and +∞ where it diverges; the moments as Bessel ratios in the scaled logs, the variance's K_{p+2}K_p/K²_{p+1} − 1 by expm1, so it keeps its digits at large √(ab).
- *`GeneralizedHyperbolic`*, x = m·v + s·√v·ζ with v ~ GIG(p, a, b) (Barndorff-Nielsen 1977). Its density is closed, m·x/s² − ½ln(2πs²) + L(p − ½, a + m²/s², b + x²/s²) − L(p, a, b) with L the GIG's log normalizer. Its tails are F(x) = ∫ Φ(z(u))·g(u) du over u = ln v, z(u) = (x − m·eᵘ)/(s·e^(u/2)), taken in logs by adaptive Gauss–Kronrod on panels widening from the integrand's peak, each tail on its own side (`IntegratedLogCdf`, `IntegratedLogSurvival`). A reader takes them from a table of the normal score z(x) = Φ⁻¹(F(x)) at points walked out from the median to the 10⁻¹² quantiles, each step moving z by at most ¼ by its exact slope f/φ(z), joined by cubic Hermite and read back as ln Φ(±z); the integral beyond. The table because a tail by the integral costs 110 to 230 µs, each of its hundreds of normal log tails an incomplete gamma, and a predictive mixture's quantile reads hundreds of tails per draw.
  - *Alternatives tried:* a mixture of normals at fixed nodes of v, which held the body to 10⁻⁵ but went wrong far out under a strong drift, where Φ's step in ln v is narrower than the nodes; Hermite on ln F itself, 10⁻⁴ off near the median. *Not tried:* a faster normal log tail in the integrand, which would retire the table.

*Checked* (`BesselKTests`, `GigNormalizerTests`, `GigMixingTests`, `BassNormalWaldTests`, `GbassPrimitiveParityTests`). Log Bessel K against its integral eˣK_ν(x) = ∫ e^(−x(cosh t − 1))·cosh(νt) dt by quadrature on 2,000 cases, ν to 10⁴ and x from 10⁻³⁰⁰ to 10³⁰⁰, to 10⁻¹² relative; the half-integer orders in closed form to 10⁻¹³; base R's `besselK` at the oracle's nine points to 4·10⁻¹⁵. The GIG's normalizer, density and moments against quadrature of its kernel on 1,500 cases, both limits included. The generalized hyperbolic's density, as GIG mixing's marginal, against quadrature of the mixture itself on 600 cases to 10⁻⁸; its tails by the integral against its density's mass at its quantiles 0.001 to 0.999 to 1.3·10⁻⁹ on six laws (the Normal-Wald with and without a drift, the skewed t, the variance gamma, two general GIGs); the table to 10⁻⁵ relative in the body and 3·10⁻⁴ at a power tail's 10⁻⁹, its quantiles' levels to 4·10⁻⁵.

### Information measures

| Measure | Definition as implemented | Conventions |
|---|---|---|
| Forward KL | $\mathrm{KL}(p\Vert q)=\sum_{p_i>0}p_i\ln(p_i/\tilde q_i)$, $\tilde q_i=\max(q_i,\varepsilon)$, $\varepsilon=10^{-12}$ | Throws on a length mismatch; entries with $p_i\le0$ are skipped; nothing is normalized, so non-negativity holds only for normalized inputs |
| Jeffreys (symmetric) KL | $\mathrm{KL}(p\Vert q)+\mathrm{KL}(q\Vert p)$ | Also the default body of `IDivergence<T>.Symmetric` |
| Jensen–Shannon distance | $\sqrt{\tfrac12(S_p+S_q)}$ with $S_p=\sum_ip_i\ln\max(p_i/m_i,10^{-15})$ and $m_i=\max(\tfrac12(p_i+q_i),10^{-15})$ | Square root of the JS divergence, in nats, at most $\sqrt{\ln2}$ for normalized inputs; branch-free SIMD via `TensorPrimitives`; the length is taken from `p`: a longer `q` is silently truncated, a shorter `q` throws `ArgumentOutOfRangeException` from an internal slice, and an empty `p` returns 0 |
| Shannon entropy | $H=-\sum_{c_i>0}\hat p_i\ln\hat p_i$ with $\hat p_i=c_i/\sum_{c_j>0}c_j$ | Counts or unnormalized probabilities; non-positive entries are ignored; returns 0 for a non-positive total |

**`JSDistance` scratch memory** depends on the input length:

| Length | Scratch source |
|---|---|
| ≤ 64 | stack (`stackalloc`) |
| 65–512 | a per-thread pair of 512-element arrays, allocated once per thread |
| > 512 | rented from `ArrayPool` |

## Configuration

| Setting | Default | Meaning |
|---|---|---|
| `Xoshiro256PlusPlus(seed)` | `null` | Seed from OS entropy (not reproducible) |
| Incomplete gamma and beta iterations | 100,000 | Cap before `ArithmeticException` |
| `LogConcaveEnvelope` attempts | 100,000 | Cap before `ArithmeticException`; with over half accepted, only a NaN density reaches it |
| `RatioOfUniforms` attempts | 100,000 | Cap before `ArithmeticException`; with half the rectangle accepted, only a wrong density reaches it |
| `LogConcaveEnvelope` design point | 32 Newton steps; the guard below h = −2²⁰, 64 steps | The guard runs only where Newton's search ends short of its point; its steps count apart from Newton's |
| `TailInversion` iterations | 100 | Cap, returning the best root found |
| `KLDivergence` eps | 1e-12 | Floor applied to q |
| `KlDivergence(eps)` | 1e-12 | Floor applied to q; `default(KlDivergence)` and `new KlDivergence()` both have eps = 0 |
| `JSDistance` clamps | 1e-15 (fixed) | Hard-coded floor on $m_i$ and on the ratio; not configurable |

## Status and caveats

- **`default(KlDivergence)` and `new KlDivergence()` have eps = 0.** The parameterless form bypasses the optional-argument constructor, so a zero entry in q (with $p_i>0$) gives $+\infty$ rather than being floored. Construct with an explicit eps, e.g. `new KlDivergence(1e-12)`.
- **Unbounded rejection loop.** The Gamma rejection loop has no iteration cap. This is standard for Marsaglia–Tsang, whose acceptance rate is above 95%.
- **Integer draws are slightly biased.** `NextInt` has no rejection step. The bias is negligible for small q.
- **Quantiles are as accurate as the tails they invert.** A prefactor's ε·S carries over, about 10⁻¹⁰ relative at ν = 10⁵, and a solve in ln S carries ε·|ln S| of rounding into x, 3×10⁻¹³ at tail −700. A shape above about 10⁸ exhausts the incomplete gamma's series, so its quantile throws. Truncated draws read no tail and hold there (shape 10¹², ν = 10¹⁴).
- **A truncated draw has setup per call.** It costs about three of Robert's on a normal half-line ([the comparison](#truncated-draws-and-quantiles)).
- **`Empirical` has no draw or quantile.** It is a measure, not an `IContinuousDistribution`.
- **`Normal.StandardQuantile` returns NaN outside [0, 1]**, where `Quantile(p)` on every family refuses. Its callers read it as a plain function.
- **Large shapes.** Accuracy falls as $\varepsilon S$ with the prefactor's size (see Algorithms). For a Student-t with very large ν, use `Normal`.
- **`default` structs are invalid.** `default(Gamma)` and the other `default` values have zero parameters. Build distributions with their factories.
- **Partial divergence abstraction.** `IDivergence<T>` is a partial abstraction with a single adapter and no consumers.

## Source map

| Role | Files |
|---|---|
| Generator and seeds | `ThermoMapper/src/maths/rng/Xoshiro256PlusPlus.cs`, `SeedTree.cs` |
| Distributions | `ThermoMapper/src/maths/distributions/Normal.cs`, `Gamma.cs`, `InverseGamma.cs`, `Weibull.cs`, `HalfCauchy.cs`, `StudentT.cs`, `Beta.cs`, `Uniform.cs`, `Mixture.cs`, `Truncated.cs`, `IContinuousDistribution.cs`, `LogConcaveEnvelope.cs`, `TailInversion.cs`, `SpecialFunctions.cs`, `BesselK.cs` (log Bessel K), `GeneralizedInverseGaussian.cs`, `GeneralizedHyperbolic.cs`, `RatioOfUniforms.cs` |
| Information | `ThermoMapper/src/maths/information/Divergences.cs` (`KLDivergence`), `IDivergence.cs`, `JSDistance.cs`, `Shannon.cs` |

**See also:** [Samplers](./Samplers.md) (ensemble seeding) · [BayesianRegression](./BayesianRegression.md) · [DesignPrinciples](./DesignPrinciples.md#6-reproducible-randomness)
