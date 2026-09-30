# Probability and Information

> ThermoMapper's stochastic primitives: a checkpointable xoshiro256++ generator and a seed-expansion tree that make every Monte Carlo component reproducible from one integer, allocation-free samplers and log-densities for the Normal, Gamma and Beta distributions, a log-gamma function, and entropy and divergence measures. · **Projects:** `Maths.Rng`, `Maths.Distributions`, `Maths.Information` (namespaces of the same names) · **Maturity:** core; the `IDivergence<T>` abstraction is partial

## Overview

These three small projects sit at the bottom of the [architecture](./Architecture.md) (layers 0–1). They are consumed widely:

| Type | Consumers |
|---|---|
| `Xoshiro256PlusPlus` | distributions, regression, [samplers](./Samplers.md), SPC, geometry, synthetic noise fields |
| `SeedTree` | the ensemble shell, parallel tempering, BAPS, GMM K-sweeps, distributed SPRED |
| `Beta` | BARS knot-locality proposals (regression) |
| `Gamma` | regression (BAPS inverse-gamma and local-multiplier draws, BARS robust weights) |
| `SpecialFunctions.LogGamma` | geometry, regression (changepoint and BARS priors) |
| `Shannon` | SPC profiling, graph observables |

The KL, Jensen–Shannon and `IDivergence` types have no in-repository callers yet.

## Key concepts

- **The project die.** The library's standard generator is `Xoshiro256PlusPlus`. Every sampler takes a generator as an explicit argument, and none holds hidden global state. A few components elsewhere use `System.Random` instead (see [Concepts](./Concepts.md#reproducibility-contract)).
- **Seed tree.** A deterministic fan-out from one master seed to N child seeds, so a parallel ensemble of chains or workers can be reproduced exactly.
- **Checkpointable state.** The generator's four 64-bit state words can be saved and restored bit for bit. The SPC checkpoint format uses this.
- **Nats.** Every entropy and divergence here uses the natural logarithm.

## Design

- **Static, stateless primitives that take an explicit RNG.** `Normal`, `Gamma`, `Beta`, `SpecialFunctions`, `KLDivergence`, `JSDistance` and `Shannon` are static and allocate nothing on their hot paths. Reproducibility is purely a seeding question: one generator per chain, seeded from the seed tree.
- **Seeding convention.** `SeedTree.Derive(master, count)` has the prefix property: `Derive(m, C)` is a prefix of `Derive(m, 2C)`. The ensemble shell uses this to give chain c a *kernel* seed from children $[0, C)$ and a separate *readout* seed from children $[C, 2C)$. See [Samplers](./Samplers.md).
- **Dispatch surface for divergences (partial).** `IDivergence<T>` is a dispatch interface over a domain T. It has one adapter, `KlDivergence`, over `ReadOnlyMemory<double>`, and no consumers. The Jensen–Shannon distance mentioned in its documentation does not implement it.

## API at a glance

| Area | Entry points |
|---|---|
| Generator | `Xoshiro256PlusPlus(int? seed)`, `NextUInt64`, `NextDouble`, `NextInt(q)`, `SaveState`, `LoadState`, `Restore` |
| Seeds | `SeedTree.Derive(master, count)` |
| Distributions | `Normal.Sample(rng[, mean, sd])`, `Normal.LogPdf`; `Gamma.Sample(rng, shape, scale = 1)`, `Gamma.LogPdf`; `Beta.Sample(rng, a, b)`, `Beta.LogPdf`; `SpecialFunctions.LogGamma` |
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
- **Gamma(a, θ)** uses Marsaglia–Tsang (2000):
  1. Set $d=a-\tfrac13$ and $c=1/\sqrt{9d}$.
  2. Draw $x\sim N(0,1)$ until $1+cx>0$, then set $v=(1+cx)^3$.
  3. Accept $d\,v\,\theta$ if the squeeze $u<1-0.0331x^4$ holds, or if $\ln u<\tfrac12x^2+d(1-v+\ln v)$.

  For shape $a<1$ it uses the boost $G(a)=G(a+1)\,u^{1/a}$, which recurses at most once. The log-density is $(a-1)\ln x-x/\theta-a\ln\theta-\ln\Gamma(a)$, and $-\infty$ for $x\le0$.
- **Beta(a, b)** is sampled as $G_a/(G_a+G_b)$ with unit-scale gammas. If the sum underflows it returns 0.5. The log-density is $(a-1)\ln x+(b-1)\ln(1-x)-\ln B(a,b)$ on $(0,1)$.
- **LogGamma** uses the Lanczos approximation with $g=7$ and nine coefficients for $x\ge 0.5$, and the reflection $\ln(\pi/\sin\pi x)-\ln\Gamma(1-x)$ below that. The documented domain is $x>0$.

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
| `Gamma.Sample` / `LogPdf` scale | 1.0 | Scale θ (mean aθ) |
| `KLDivergence` eps | 1e-12 | Floor applied to q |
| `KlDivergence(eps)` | 1e-12 | Floor applied to q; `default(KlDivergence)` and `new KlDivergence()` both have eps = 0 |
| `JSDistance` clamps | 1e-15 (fixed) | Hard-coded floor on $m_i$ and on the ratio; not configurable |

## Status and caveats

- **`default(KlDivergence)` and `new KlDivergence()` have eps = 0.** The parameterless form bypasses the optional-argument constructor, so a zero entry in q (with $p_i>0$) gives $+\infty$ rather than being floored. Construct with an explicit eps, e.g. `new KlDivergence(1e-12)`.
- **Unbounded rejection loop.** The Gamma rejection loop has no iteration cap. This is standard for Marsaglia–Tsang, whose acceptance rate is above 95%.
- **Integer draws are slightly biased.** `NextInt` has no rejection step. The bias is negligible for small q.
- **Unchecked documentation claims.** The ~15-digit accuracy claimed for `LogGamma` is the author's statement and was not independently checked.
- **Partial divergence abstraction.** `IDivergence<T>` is a partial abstraction with a single adapter and no consumers.

## Source map

| Role | Files |
|---|---|
| Generator and seeds | `ThermoMapper/src/maths/rng/Xoshiro256PlusPlus.cs`, `SeedTree.cs` |
| Distributions | `ThermoMapper/src/maths/distributions/Normal.cs`, `Gamma.cs`, `Beta.cs`, `SpecialFunctions.cs` |
| Information | `ThermoMapper/src/maths/information/Divergences.cs` (`KLDivergence`), `IDivergence.cs`, `JSDistance.cs`, `Shannon.cs` |

**See also:** [Samplers](./Samplers.md) (ensemble seeding) · [BayesianRegression](./BayesianRegression.md) · [DesignPrinciples](./DesignPrinciples.md#6-reproducible-randomness)
