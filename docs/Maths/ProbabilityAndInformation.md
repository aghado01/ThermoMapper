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
| Distributions | `Normal.Standard`, `Normal.MeanSd(mean, sd)`, `Normal.MeanVariance`, `Normal.SampleStandard(rng)`, `Normal.StandardQuantile(p)` (Wichura's AS 241, about 1 part in 10¹⁶; checked against the log-space CDF to p = 10⁻³⁰⁰), `Normal.StandardInverseLogCdf(logP)`; `Gamma.ShapeRate(shape, rate)`, `Gamma.ShapeScale(shape, scale)`, `Gamma.SampleStandard(rng, shape)`; `InverseGamma.ShapeScale(shape, scale)`; `Weibull.ShapeScale(shape, scale)` with `Median` and `Mode`; `StudentT.Standard(ν)`, `StudentT.LocationScale(ν, location, scale)`; `Beta.Shapes(a, b)`; `HalfCauchy.WithScale(s)`; `Uniform.On(lower, upper)`; `Mixture.Of(components, weights)`; `Truncated.Of(parent, lower, upper)`. Each is an `IContinuousDistribution`: `Sample(rng)`, `SampleOn(rng, lower, upper)`, `LogPdf`, `Cdf`, `Survival`, `LogCdf`, `LogSurvival`, `IntervalProbability`, `LogIntervalProbability`, `Quantile(p)`, `InverseLogCdf(logP)`, `InverseLogSurvival(logS)`, and its moments where it has them. `Empirical.Of(sample)`; `IIntervalMoments.IntervalMoments(lower, upper, centre, logScale, moments)`, implemented by `Normal`, `StudentT`, `Truncated<T>`, `Uniform`, `Mixture` and `Empirical`; `GaussLegendre.Of(n)` |
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
- *Where R falls short, corrected rather than ported:*
  - `qt` stops short in the far tail: 1% in t at ν = 1.5 and tail −700, where R's own tail reads −700.015, and 8×10⁻⁹ at ν = 3.
  - `qt` overflows to ∞ at ν = ½ and tail −230, where t = 1.03×10¹⁹⁹ holds it.
  - `qbeta` returns a floor of 2⁻¹⁰²³ where the root lies below the double range, as at Beta(0.01, 0.01) and tail −23, whose tail at that floor is −7.8. This side returns 0.

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
| Distributions | `ThermoMapper/src/maths/distributions/Normal.cs`, `Gamma.cs`, `InverseGamma.cs`, `Weibull.cs`, `HalfCauchy.cs`, `StudentT.cs`, `Beta.cs`, `Uniform.cs`, `Mixture.cs`, `Truncated.cs`, `IContinuousDistribution.cs`, `LogConcaveEnvelope.cs`, `TailInversion.cs`, `SpecialFunctions.cs` |
| Information | `ThermoMapper/src/maths/information/Divergences.cs` (`KLDivergence`), `IDivergence.cs`, `JSDistance.cs`, `Shannon.cs` |

**See also:** [Samplers](./Samplers.md) (ensemble seeding) · [BayesianRegression](./BayesianRegression.md) · [DesignPrinciples](./DesignPrinciples.md#6-reproducible-randomness)
