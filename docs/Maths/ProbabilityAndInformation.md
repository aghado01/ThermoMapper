# Probability and Information

> ThermoMapper's stochastic primitives: a checkpointable xoshiro256++ generator and a seed-expansion tree that make every Monte Carlo component reproducible from one integer, allocation-free samplers and log-densities for the Normal, Gamma and Beta distributions, a log-gamma function, and entropy and divergence measures. · **Projects:** `Maths.Rng`, `Maths.Distributions`, `Maths.Information` (namespaces of the same names) · **Maturity:** core; the `IDivergence<T>` abstraction is partial

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
| `StudentT`, `Truncated<T>` | no callers yet; built for BASS sensitivity input measures and calibration priors |
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

- **Value types with named factories.** Each distribution is a readonly struct built by a factory that names its parameterization: `Normal.MeanSd`, `Gamma.ShapeRate` or `Gamma.ShapeScale`, `InverseGamma.ShapeScale`, `StudentT.LocationScale`, `Beta.Shapes`.
  - Conjugate updates produce rates, so a caller passes a rate as a rate rather than converting it at the call site. That conversion is where R's BASS τ prior slipped ([BASS derivations review §6](../../issues/regression/bass-derivations-and-integration-review.md#6-d03-finding-rs-default-τ-hyperprior)).
  - The structs allocate nothing and take the generator explicitly. `Normal.SampleStandard` and `Gamma.SampleStandard` are the unit draws.
  - `SpecialFunctions`, `KLDivergence`, `JSDistance` and `Shannon` stay static.
  - Reproducibility is purely a seeding question: one generator per chain, seeded from the seed tree.
- **Log-space tails.** CDFs, survival functions and interval probabilities come from whichever tail represents them directly, in log space, so a probability far below the double range keeps its digits as a log. `Truncated<T>` normalizes by its parent's interval probability.
- **Seeding convention.** `SeedTree.Derive(master, count)` has the prefix property: `Derive(m, C)` is a prefix of `Derive(m, 2C)`. The ensemble shell uses this to give chain c a *kernel* seed from children $[0, C)$ and a separate *readout* seed from children $[C, 2C)$. See [Samplers](./Samplers.md).
- **Dispatch surface for divergences (partial).** `IDivergence<T>` is a dispatch interface over a domain T. It has one adapter, `KlDivergence`, over `ReadOnlyMemory<double>`, and no consumers. The Jensen–Shannon distance mentioned in its documentation does not implement it.

## API at a glance

| Area | Entry points |
|---|---|
| Generator | `Xoshiro256PlusPlus(int? seed)`, `NextUInt64`, `NextDouble`, `NextInt(q)`, `SaveState`, `LoadState`, `Restore` |
| Seeds | `SeedTree.Derive(master, count)` |
| Distributions | `Normal.Standard`, `Normal.MeanSd(mean, sd)`, `Normal.MeanVariance`, `Normal.SampleStandard(rng)`, `Normal.StandardQuantile(p)` and `Quantile(p)` (Wichura's AS 241, about 1 part in 10¹⁶; checked against the log-space CDF to p = 10⁻³⁰⁰); `Gamma.ShapeRate(shape, rate)`, `Gamma.ShapeScale(shape, scale)`, `Gamma.SampleStandard(rng, shape)`; `InverseGamma.ShapeScale(shape, scale)`; `StudentT.Standard(ν)`, `StudentT.LocationScale(ν, location, scale)`; `Beta.Shapes(a, b)`. Each has `Sample(rng)`, `LogPdf`, `Cdf`, `Survival`, `LogCdf`, `LogSurvival`, `IntervalProbability`, `LogIntervalProbability` and its moments. `Truncated.Of(parent, lower, upper)`; `IContinuousDistribution` |
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
- **Truncation.** `Truncated<T>` has density $f(x)/P(\text{lower}\le X\le\text{upper})$ on the interval, with the normaliser taken from the parent's log interval probability. N(0, 1) on $[10,\infty)$ therefore has density $\varphi(x)/Q(10)$, where $1-\Phi(10)$ is exactly 0 in double precision.
- **Accuracy**, checked against mpmath at 50–60 digits (`DistributionReferenceTests`):
  - **LogGamma** is within $2\times10^{-14}\max(1,|\ln\Gamma|)$ at ten points from $10^{-8}$ to $10^{5}$.
  - **Incomplete gamma and beta** are limited by their log prefactor, a sum of terms of size $S$ ($|\ln\Gamma(a)|$, $|a\ln x|$, …). Their relative error is about $\varepsilon S$: about $10^{-13}$ for moderate shapes and $10^{-12}$ at shape 1000.
  - **Far tails** keep their digits in log space, for example $\ln\Phi(-40)=-804.6084\ldots$ and $\ln P(T>10^{200})$ for ν = 1.
  - **A same-tail interval** loses relative accuracy as it narrows: about $\varepsilon\max(1,|\alpha|)/(\beta-\alpha)$ for the normal.

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
| `KLDivergence` eps | 1e-12 | Floor applied to q |
| `KlDivergence(eps)` | 1e-12 | Floor applied to q; `default(KlDivergence)` and `new KlDivergence()` both have eps = 0 |
| `JSDistance` clamps | 1e-15 (fixed) | Hard-coded floor on $m_i$ and on the ratio; not configurable |

## Status and caveats

- **`default(KlDivergence)` and `new KlDivergence()` have eps = 0.** The parameterless form bypasses the optional-argument constructor, so a zero entry in q (with $p_i>0$) gives $+\infty$ rather than being floored. Construct with an explicit eps, e.g. `new KlDivergence(1e-12)`.
- **Unbounded rejection loop.** The Gamma rejection loop has no iteration cap. This is standard for Marsaglia–Tsang, whose acceptance rate is above 95%.
- **Integer draws are slightly biased.** `NextInt` has no rejection step. The bias is negligible for small q.
- **No truncated sampler and no quantiles.** `Truncated<T>` cannot draw, and no distribution has a quantile function. No caller needs either yet.
- **No interval moments.** BASS sensitivity needs truncated means and higher raw moments (port brief D08). They are not implemented.
- **Large shapes.** Accuracy falls as $\varepsilon S$ with the prefactor's size (see Algorithms). For a Student-t with very large ν, use `Normal`.
- **`default` structs are invalid.** `default(Gamma)` and the other `default` values have zero parameters. Build distributions with their factories.
- **Partial divergence abstraction.** `IDivergence<T>` is a partial abstraction with a single adapter and no consumers.

## Source map

| Role | Files |
|---|---|
| Generator and seeds | `ThermoMapper/src/maths/rng/Xoshiro256PlusPlus.cs`, `SeedTree.cs` |
| Distributions | `ThermoMapper/src/maths/distributions/Normal.cs`, `Gamma.cs`, `InverseGamma.cs`, `StudentT.cs`, `Beta.cs`, `Truncated.cs`, `IContinuousDistribution.cs`, `SpecialFunctions.cs` |
| Information | `ThermoMapper/src/maths/information/Divergences.cs` (`KLDivergence`), `IDivergence.cs`, `JSDistance.cs`, `Shannon.cs` |

**See also:** [Samplers](./Samplers.md) (ensemble seeding) · [BayesianRegression](./BayesianRegression.md) · [DesignPrinciples](./DesignPrinciples.md#6-reproducible-randomness)
