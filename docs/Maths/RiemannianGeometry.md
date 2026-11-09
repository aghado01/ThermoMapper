# Riemannian Geometry

> ThermoMapper's statistics on curved spaces. One manifold contract covers the sphere, two hyperbolic models, the Grassmannian and Stiefel manifolds, and product and scaled composites. On top of that contract a single IRLS engine computes Fréchet means, geometric medians and scatter matrices on any of these spaces. The project also provides heat-kernel bandwidth calibration on hyperbolic and spherical space, robust median-of-means PCA, simulated annealing over subspaces (the engine behind SPRED), and a two-sample test. · **Project:** `Maths.Geometry` · **Maturity:** core; `StiefelManifold` and `MxPbf` are partial

## Overview

`Maths.Geometry` sits in layer 3 of the [architecture](./Architecture.md). Its dependencies:

- the geodesic distance kernels in [`Maths.Distance`](./DistanceMetrics.md);
- [`Maths.LinAlg`](./LinearAlgebra.md) for eigen, PCA and Gram–Schmidt;
- `Maths.Distributions` for log-gamma;
- `Maths.Rng`;
- `System.Numerics.Tensors`.

Its consumers:

| Consumer | Uses |
|---|---|
| [GraphConstruction](./GraphConstruction.md) | manifolds and bandwidth calibrators |
| [Spred](./Spred.md) | `SubspaceAnnealer`, `GrassmannManifold`, `GeometricMedian` |
| [GaussianMixtures](./GaussianMixtures.md) | `EuclideanMedian` for robust initialization |
| [SyntheticDatasets](./SyntheticDatasets.md) | `PoincareBallManifold` |

**Namespaces.**

| Namespace | Contents |
|---|---|
| `Maths.Geometry` | manifolds and the contract |
| `Maths.Geometry.Solver` | IRLS, losses, options |
| `Maths.Geometry.Bandwidth` | heat-kernel calibration |
| `Maths.Geometry.Estimators.Intrinsic` / `.Ambient` / `.Calibration` | location and scatter estimators |
| `Maths.Geometry.DimReduction` | MoMPCA, SubspaceAnnealer |
| `Maths.Geometry.Inference` | mxPBF |

```mermaid
graph TD
  IRM[IRiemannianManifold] --- E[EuclideanVectorManifold] & S[SphericalManifold] & P[PoincareBallManifold] & L[LorentzHyperboloidManifold] & G[GrassmannManifold] & St[StiefelManifold] & Pr["ProductManifold&lt;TA,TB&gt;"] & Sc["ScaledManifold&lt;T&gt;"]
  GM[GeometricMean] -->|L2Loss| Irls["Irls.Solve&lt;TManifold,TLoss&gt;"]
  GMed[GeometricMedian] -->|L1Loss| Irls
  Irls --> IRM
  GM --> KS[KarcherScatter] --> SA[ScatterAccumulator]
  GMed --> WS[WeiszfeldScatter] --> SA
  CF[ConsistencyFactors] -.c.-> SA
  MoM[MoMPCA] --> GMed & Pr & Sc
  Ann[SubspaceAnnealer] --> G
  Hyp[Bandwidth.Hyperbolic] --> HK[HeatKernel]
  Sph[Bandwidth.Spherical] -.->|β converted by caller| HK
```

## Key concepts

- **Log and Exp maps.** $\mathrm{Log}_p(q)$ is the tangent vector at p pointing along the geodesic to q, with length $d(p,q)$. $\mathrm{Exp}_p(v)$ is the point reached by following v from p.
- **`Dimension` vs `IntrinsicDimension`.** `Dimension` is the **buffer length**, meaning the number of ambient coordinates: n for the sphere and hyperboloid, n·r for Grassmann and Stiefel. `IntrinsicDimension` is the manifold dimension: n−1, r(n−r), and nr − r(r+1)/2. For embedded manifolds, `Dimension` is larger than the tangent-space dimension. The code refers to keeping the two apart as an "off-by-one ward".
- **"Intrinsic" and "ambient" estimators.** These names mean *computed through Log/Exp* versus *computed in plain ℝᴰ coordinates*. They do not refer to intrinsic-dimension estimation; this module has no intrinsic-dimension estimators.
- **IRLS weight.** For a loss ρ with $\psi=\rho'$, the weight is $w(r)=\psi(r)/r$. That is 1 for the L2 loss and $1/r$ for the L1 loss.
- **Warm start in `destination`.** Every iterative location estimator (`GeometricMedian`, `GeometricMean` on curved manifolds, `EuclideanMedian`) uses the `destination` buffer as both its initial iterate and its output. `EuclideanMean` and the flat-space `GeometricMean` simply overwrite it.

## Design

**Compile-time dispatch through generic struct parameters.** The manifold and loss contracts use C# static abstract members: `IsFlat` on manifolds, and `IsClosedForm`, `IsSingularAtZero` and `Weight` on losses. Consumers constrain `TManifold : struct, IRiemannianManifold` and `TLoss : struct, IRobustLoss`. The JIT therefore specializes each value-type instantiation, and checks such as `if (TManifold.IsFlat)` become constants. The code calls this "three-axis dispatch": manifold × loss × options.

**Composable manifolds.** `ProductManifold<TA,TB>` and `ScaledManifold<TInner>` are themselves struct manifolds, so they nest. MoMPCA, for example, runs the whole solver on `ProductManifold<ScaledManifold<EuclideanVectorManifold>, ScaledManifold<GrassmannManifold>>`.

**Spans as currency.**

- Points and tangents are `ReadOnlySpan<double>` / `Span<double>` of length `Dimension`.
- Solver scratch comes from `ArrayPool<double>.Shared`.
- The sphere and the Poincaré ball use `stackalloc` for D ≤ 256.
- Grassmann and Stiefel operations allocate small matrices on each call.

**Declarative options with pure engines.** `IrlsOptions` is a readonly struct with a `Default`, and `Irls` is a static engine. `SubspaceAnnealerOptions` is a data record whose `Validate()` is shared by every consumer.

**One scatter primitive.** Karcher (L2) and Weiszfeld (L1) scatter differ only in the weights they pass to the internal `ScatterAccumulator`.

**Extension points.**

- New manifold: implement `IRiemannianManifold` as a struct.
- New robust loss: implement `IRobustLoss`. Only `L1Loss` and `L2Loss` exist today.
- New subspace objective: pass a `SubspaceObjectiveFunction` delegate.

## API at a glance

| Entry point | Essentials |
|---|---|
| `IRiemannianManifold` | `static abstract bool IsFlat`; `Dimension`; `Distance(p, q)`; `LogMap(p, q, dst)`; `ExpMap(p, v, dst)`; `AddScaled(dst, v, s)`; `Norm(p, v)` |
| Manifolds | `EuclideanVectorManifold(d)`, `SphericalManifold(n ≥ 2)`, `PoincareBallManifold(d ≥ 1)`, `LorentzHyperboloidManifold(n ≥ 2)`, `GrassmannManifold(n, r)`, `StiefelManifold(n, r)`, `ProductManifold<TA,TB>(a, b)`, `ScaledManifold<T>(inner, scale > 0)` |
| Solver | `Irls.Solve<TManifold,TLoss>(manifold, data, weights, destination, IrlsOptions, finalIrlsWeights = default)`; `IrlsOptions`; `HybridMode`, `SingularityPolicy`, `ConvergenceCriterion`; `L1Loss`, `L2Loss` |
| Intrinsic estimators | `GeometricMean.Compute` / `ComputeWithScatter`; `GeometricMedian.Compute` / `ComputeWithScatter` / `MedoidIndex`; `KarcherScatter.Compute`; `WeiszfeldScatter.Compute` |
| Ambient estimators | `EuclideanMean`, `EuclideanMedian` (over an index subset of `double[][]`); `ConsistencyFactors.Gaussian(D)`, `ConsistencyFactors.Laplace(D)` |
| Bandwidth | `Bandwidth.Hyperbolic.MatchBeta(m̂₂, d)` / `SecondMoment`; `Bandwidth.Spherical.MatchBeta(m̂₂, m)` / `SecondMoment`; `HeatKernel.BandwidthFromBeta(β)` |
| Subspaces | `MoMPCA.ComputeMoM(data, kBlocks, numComponents)` → `PcaResult`; `SubspaceAnnealer.Compute(data, targetDim, objective, maxIters, seed, options, ct)` → `SubspaceAnnealerResult(Projection, Objective)` |
| Inference | `MxPbf.TwoSampleMean(x, y, alpha)`, `MxPbf.TwoSampleCovariance(x, y, alpha, a0, b0)` → log mxPBF |

## Algorithms

### Manifold catalog

In the table, $\hat p=p/\|p\|$, and $\langle a,b\rangle_L$ is the Lorentz product with the time-like coordinate **last**. No small-angle series are used; every guard is a hard threshold or a clamp.

| Manifold | Points / tangents | Distance | Log | Exp | Guards |
|---|---|---|---|---|---|
| Euclidean | ℝᵈ | $\Vert q-p\Vert$ | $q-p$ | $p+v$ | none |
| Sphere $S^{n-1}$ | any nonzero vector, normalized internally / vectors ⟂ p̂ (not enforced) | $\arccos\langle\hat p,\hat q\rangle$ | $\frac{\theta}{\sin\theta}(\hat q-\langle\hat p,\hat q\rangle\hat p)$ | $\mathrm{normalize}(\cos\Vert v\Vert\,\hat p+\frac{\sin\Vert v\Vert}{\Vert v\Vert}v)$ | Log returns **0** for $\theta<10^{-12}$ or near-antipodal $\sin\theta<10^{-12}$; Log/Exp normalize a zero vector to $e_0$, but `Distance` returns π/2 for a zero vector |
| Poincaré ball $H^d$ | $\Vert p\Vert<1$ | $\operatorname{arcosh}\!\big(1+\frac{2\Vert a-b\Vert^2}{(1-\Vert a\Vert^2)(1-\Vert b\Vert^2)}\big)$ | $\frac{2}{\lambda_p}\operatorname{artanh}\Vert\delta\Vert\frac{\delta}{\Vert\delta\Vert}$, $\delta=(-p)\oplus q$ | $p\oplus\tanh(\tfrac{\lambda_p\Vert v\Vert}{2})\frac{v}{\Vert v\Vert}$, then clamped to the ball | inputs clamped to radius $1-10^{-5}$; $\lambda_p=2/\max(10^{-12},1-\Vert p\Vert^2)$ |
| Lorentz hyperboloid $H^d$ | $\langle p,p\rangle_L=-1$ | $\operatorname{arcosh}(-\langle p,q\rangle_L)$ | $\frac{\operatorname{arcosh}\alpha}{\sqrt{\alpha^2-1}}(q-\alpha p)$ | $\cosh\Vert v\Vert_L\,p+\frac{\sinh\Vert v\Vert_L}{\Vert v\Vert_L}v$, time coordinate re-projected | Log returns 0 for $\alpha\le1+10^{-12}$ |
| Grassmann $Gr(r,n)$ | n×r orthonormal (column-major) / horizontal ($Y^\top\Delta=0$, not enforced) | $\big(\sum_i\arccos^2\sigma_i(Y^\top Z)\big)^{1/2}$ | $(U_2\Theta)\,U_1^\top$, where $Y^\top Z=U_1\cos\Theta\,V_1^\top$ and column i of $U_2\Theta$ is $\frac{\theta_i}{\sin\theta_i}\big((ZV_1)_i-c_i(YU_1)_i\big)$ | $\mathrm{MGS}\big((YV\cos\Sigma+U\sin\Sigma)V^\top\big)$ | a principal angle near π/2 ($\sigma_i=\cos\theta_i<10^{-12}$) leaves a zero column in $U_1$, so that direction is **dropped from Log**, giving a shortened tangent, not a zero vector; angles with $\sin\theta_i<10^{-12}$ ($\theta_i\approx0$) contribute nothing |
| Product | [A-part \| B-part] | $\sqrt{d_A^2+d_B^2}$ | factor-wise | factor-wise | `IsFlat` = both factors flat |
| Scaled | as inner | $\sqrt s\,d$ | unchanged | unchanged | $s$ finite and $>0$ |

**Grassmann SVD.** Grassmann and Stiefel take a thin SVD from [`ThinSvd`](./LinearAlgebra.md#thin-svd-thinsvd), which never forms the r×r Gram matrix. They zero the left vectors below σ = 10⁻¹², a rule the log maps rely on at the cut locus. The Grassmann Log matches Bendokat–Zimmermann–Absil and the Exp matches Edelman–Arias–Smith.

**Cost.** Sphere, ball and hyperboloid operations are $O(D)$. Grassmann operations are $O(nr^2)$ plus Jacobi sweeps on the r×r triangular factor.

**Stiefel (partial).** The type's doc comments promise orientation-sensitive canonical (Edelman–Arias–Smith) geometry. The code does something different:

- **Distance** depends only on the column spans, so $d(Y,YQ)=0$ for any orthogonal Q.
- **Log** is an approximation.
- **Exp** is a Gram–Schmidt retraction $\mathrm{MGS}(Y+\Delta)$.

### The IRLS engine (`Irls.Solve`)

`Irls.Solve` minimizes $F(p)=\sum_iw_i\rho(d(p,x_i))$. With L2 the minimizer is the Fréchet/Karcher mean. With L1 it is the geometric median. The manifold enters only through `Distance`, `LogMap`, `ExpMap`, `AddScaled`, `Norm` and `IsFlat`. The distance-valued constants are read in a unit λ, by default the data's spread σ ([below](#the-unit-of-the-constants-irlsoptionsscale)), so the flowchart writes them as τλ, ελ and θλ.

```mermaid
flowchart TD
  A[start: p = destination] --> B{loss closed-form? L2}
  B -->|yes, flat| C[weighted average — done]
  B -->|yes, curved| D["Karcher flow: p ← Exp_p(Σw Log_p x / Σw) until shift < τλ"]
  B -->|no, L1| E["radii r_i (‖Log_p x_i‖ under Spread, d(p,x_i) under Data); λ; ω_i = 0 if r_i ≤ ελ (coincident) else w_i/r_i; minDist includes coincident points"]
  E --> F{"coincidence and policy = OptimalityCheck?"}
  F -->|check passes| G[return p — p is the optimal data point]
  F -->|fails / not run| H{"SubgradientOnly, or Hybrid and minDist < θλ?"}
  H -->|yes| I["subgradient step: p ← Exp_p(η_k Σ w_i Log_p x_i / r_i), skipping coincident points"]
  H -->|no| J["Weiszfeld step: p ← Exp_p(Σ ω_i Log_p x_i / Σ ω_i)"]
  I --> K{"shift < τλ, or the step left p unchanged?"}
  J --> K
  K -->|no| E
  K -->|yes| L[return p; final weights]
```

**How each step works.**

- **Optimality check.** This is the Weiszfeld / Vardi–Zhang condition at a coincident data point c: $\big\|\sum_{j\ne c}\omega_j\mathrm{Log}_px_j\big\|\le w_c$. If it holds, p is optimal and the solver returns.
- **Coincident points always trigger the subgradient step under `Hybrid`.** A coincident point counts toward `minDist`, which is then at most ελ. That is below θλ whenever ε < θ, as in `Default` (1e-10 < 1e-4). `WeiszfeldOnly` never takes a subgradient step; `SubgradientOnly` takes only subgradient steps and skips the optimality check. So whenever the optimality check does not return, the next step is a subgradient step.
- **The subgradient step.** Under `Spread` it is $\eta_0\sigma/\sqrt{k+1}$ times the weighted mean $\sum_iw_i\,\mathrm{Log}_px_i/r_i\,/\sum_iw_i$ of the unit directions, over the points outside the coincidence radius ελ, where k is the global iteration count. Under `Data` it keeps its old form: the skip radius is a fixed $10^{-12}$, not ε (a point with $10^{-12}\le r\le\varepsilon$ is excluded from the Weiszfeld weights but receives weight $w/r$ in the subgradient direction), and the step, $\eta_0/\sqrt{k+1}$ times the sum, is normalized neither by $\sum w$ nor by the data scale.
- **Convergence.** `Absolute` stops when shift < τλ. `RelativeToNorm` stops when shift < τ(λ + ‖p‖₂), using the ambient norm; under `Data` that is τ(1 + ‖p‖₂). Under `Spread` the shift is the norm of the tangent step, and the solver also stops when the step moves the iterate by at most $4\varepsilon_{\text{mach}}\|p\|_2$ in ambient coordinates. Under `Data` the shift is $d(p_k,p_{k+1})$.

**Final weights** (`finalIrlsWeights`) depend on how the solver exits:

| Exit | Final weights |
|---|---|
| L1, normal exit | recomputed at the returned point |
| L1, optimality-check exit | that iteration's weights, with $\omega_c=0$ |
| L2 | the outer weights |

**Cost per L1 iteration.** Under `Spread`: n Log maps and n norms (the log vectors are kept, n·D doubles, and reused by the step), one sort of n radii, and one Exp. Under `Data`: n distances, n Log maps and one Exp; subgradient steps add n norms. When the optimality check runs, up to $2n-1$ more Log maps.

#### The unit of the constants (`IrlsOptions.Scale`)

The four distance-valued constants, `Tolerance` τ, `Epsilon` ε, `SubgradientThreshold` θ and `Eta0` η₀, are read in a unit λ:

- **`Spread`** (the default, and the enum's zero value): λ = σ, the weighted lower median of the radii from the current iterate to the data, over the points of positive weight, recomputed each iteration. When half the weight or more sits on the iterate, σ is their weighted mean radius instead; when all of it does, the iterate is returned.
- **`Data`**: λ = 1, the constants in the data's units. This is the solver as it was before the spread unit, reproduced bit for bit (pinned in `IrlsScaleTests`).

**Why.** With λ = 1, the median of a cluster whose spread is below θ = 1e-4 takes a subgradient step of size up to η₀ = 1 at every iteration, runs to `MaxIterations` and lands about 5% of the spread off, on every manifold ([the Grassmann discussion's probe](../../issues/regression/discussion/grassmann-stiefel-subspaces-brief.md#report)). Under `Spread` every comparison the solver makes is between two quantities in the same unit, so the median and mean of c + s·u are c + s·(the estimate for u) to a relative accuracy that does not depend on s.

**The choice of σ.** The scale has to resist outliers, since the median exists to resist them.

| Scale | Cost | Robustness | Behaviour |
|---|---|---|---|
| Median radius from the iterate, each iteration (chosen) | one sort of n radii per iteration; the radii are computed anyway | ends as the median distance from the median (50% breakdown), whatever the warm start | the thresholds move with the iterate early on, while the shifts are still large |
| Median radius from the warm start, once | one sort | inherits the warm start: a warm start pulled toward outliers, such as an arithmetic mean under contamination, inflates σ up to its distance from the bulk, and the relative accuracy loses that factor | stationary thresholds |
| Median pairwise distance | O(n²) distances | about 29% breakdown, independent of the warm start | too costly for GMM's clusters, and the per-iteration median reaches the same robustness |
| Data units (`Data`) | none | not a scale | not equivariant; kept as the option for the old semantics |

**What else `Spread` changes, and why.**

- **Radii from the Log map.** The Weiszfeld weights must resolve τσ. On the sphere and the Grassmannian `Distance` is an arccos, with relative error about $\varepsilon_{\text{mach}}/2\theta^2$ (5e-7 at θ = 1e-5), which stalls the median below a spread of 1e-4; $\|\mathrm{Log}_px\|$ holds to about 1e-10 down to 1e-7.
- **The shift and the stall stop.** The shift is the norm of the tangent step. The sphere's and the Grassmannian's Exp maps drop tangent steps below 1e-12 (an absolute guard), and an iterate far from the origin cannot move by less than its own ulps; in both cases the step stays above τσ while the iterate stops moving, so the solver stops when the step moves it by at most $4\varepsilon_{\text{mach}}\|p\|_2$.
- **Flat steps through the tangent round trip.** $p+\sum_i\omega_i(x_i-p)/\sum_i\omega_i$ instead of $\sum_i\omega_ix_i/\sum_i\omega_i$, so a centre far from the origin costs only the rounding of p, where the direct form carries n roundings of the centre. The flat mean takes the same form from the warm start.
- **The subgradient step divided by Σw**, so it is invariant to the weights' scale as well as the data's.

**Measured** (the seeded sweeps of `IrlsScaleTests`, and the chip's probe in [the brief's report](../../issues/geometry/irls-scale-brief.md#report)):

- On ℝᵈ the median's relative error is below 1e-7 and the mean's below 1e-14, the same at every s from 1e-12 to 1e6, for centres at 0 or at 100s. A centre 10 from the origin limits the error to a few hundred ulps of the centre (8e-4 of the spread for the mean at s = 1e-12, where the data themselves are represented to 2e-3).
- On S⁹, Gr(2, 10) and Gr(3, 20) the median holds about 1e-8 down to s = 1e-4 in 7 to 12 iterations, then about 1e-12/s, the Exp guard (1e-6 at s = 1e-6). Below about 1e-6 the arccos in the Log map drops principal angles under 1e-8 and the iteration counts are no longer bounded: that is the manifolds' own resolution, a separate item. The mean holds 1e-8 down to 1e-7.
- Under `Data` the same curved cases run to 200 iterations below s = 1e-4 (1e-3 on S²) and miss by 2% to 11% of the spread.

### Location, scatter and consistency

- **`GeometricMean`** runs IRLS with `L2Loss`, and **`GeometricMedian`** runs it with `L1Loss`. `GeometricMean` is the Fréchet mean, not $\exp\circ\text{mean}\circ\log$. `MedoidIndex` returns $\arg\min_i\sum_jd(x_i,x_j)$ as a robust warm start.
- **Scatter.** $\hat\Sigma=\frac{c}{\sum_i\omega_i}\sum_i\omega_iv_iv_i^\top$ with $v_i=\mathrm{Log}_\mu x_i$, as a row-major D×D matrix. On embedded manifolds this matrix is ambient-sized and rank-deficient.
- **Consistency factors.** These make the 1/r-weighted (Weiszfeld) scatter consistent for the identity ($\hat\Sigma\to I$ as $n\to\infty$) under a reference distribution: $c_D=D\,\mathbb E[1/r]/\mathbb E[r]$. They require D ≥ 2.
  - Gaussian reference: $c_D=D/(D-1)$.
  - Spherical Laplace reference: $c_D=1/(D-1)$.
  - Both apply to **L1 (Weiszfeld)** scatter.
- **Ambient estimators.** `EuclideanMean` and `EuclideanMedian` take an index subset of a `double[][]`, so they are not drop-in replacements for the intrinsic estimators. `EuclideanMedian` also takes explicit `maxIterations`/`tolerance`/`epsilon` and a `scale` (`IrlsScale.Spread` by default), read as in the [IRLS engine](#the-unit-of-the-constants-irlsoptionsscale).
  - `EuclideanMedian` is a regularized Weiszfeld, $\omega_i=w_i/\max(r_i,\varepsilon\lambda)$. It has no optimality check. Under `Spread` it steps by the weighted mean of $x_i-p$ and stops when the step is below $\max(\tau\sigma,4\varepsilon_{\text{mach}}\|p\|_2)$; under `Data` it stops at τ(1 + ‖p‖₂) of the previous iterate, as before.
  - Unlike `IrlsOptions`, both validate finite, non-negative weights. `EuclideanMedian` also requires `maxIterations` ≥ 1 and a finite, positive tolerance and ε.

### Heat-kernel bandwidth calibration

The goal is a Gaussian heat-kernel parameter β, with bandwidth $\delta=1/\sqrt{2\beta}$ (`HeatKernel.BandwidthFromBeta`), such that the model's nearest-neighbor second moment matches an observed value $\hat m_2$.

**Radial model** (Van Vleck-corrected short-time heat kernel times the polar volume element):

$$H^d:\ w(r)=e^{-\beta r^2}(r\sinh r)^{(d-1)/2},\qquad S^m:\ w(r)=e^{-\beta r^2}(r\sin r)^{(m-1)/2},\ 0<r<\pi.$$

**Nearest-neighbor reference.** A Poisson(λ = 32) population is drawn from $w/W$. Conditioned on at least one point, the nearest distance has density

$$f_{NN}(r)=\frac{\lambda}{W}w(r)e^{-\lambda C(r)/W}\Big/(1-e^{-\lambda}),\qquad C(r)=\int_0^rw.$$

`SecondMoment` evaluates $M_2(\beta)=\int r^2f_{NN}$ by trapezoidal quadrature with 1,024 intervals on $[0,r_{\max}]$, where $r_{\max}=\max\big(8,\tfrac{(d-1)\delta^2}{2}+12\delta+d\big)$ on $H^d$ and $\pi-10^{-9}$ on $S^m$.

**Root finding (`MatchBeta`).**

1. Initial guess: $d/(2\hat m_2)$, the Euclidean identity.
2. Expand the bracket geometrically.
3. Bisect 48 times.

On $S^m$ the result saturates at `BetaFloor` = 1e-6 when $\hat m_2$ exceeds what the sphere can produce. If the bracket cannot be established, the solver returns silently.

**Dimension arguments.** Pass the **intrinsic** dimension: `SphericalManifold.IntrinsicDimension` for $S^{n-1}\subset\mathbb R^n$, and d for $H^d$. The hyperbolic parameter is named `ambientDimension` but is used as d, so a Lorentz caller must pass n−1.

### Median-of-means PCA (`MoMPCA`)

"MoM" is median-of-means: block estimates aggregated by a geometric median.

1. Split the data into K **contiguous** blocks (no shuffling) and run a local PCA on each, giving a mean $\mu_k$ and a frame $U_k$.
2. Compute preliminary medians: a Euclidean median of the $\mu_k$ and a Grassmann median of the $U_k$, the latter warm-started at the medoid.
3. Calibrate scales from the **square of the median** distances: $\hat\tau_\mu=(\mathrm{med}_k\|\mu_k-\tilde\mu\|)^2/p$ and $\hat\tau_U=(\mathrm{med}_kd_{Gr}(U_k,\tilde U))^2/[r(p-r)]$. Then set $\hat\alpha=\mathrm{clamp}(2\hat\tau_U/(\hat\tau_\mu+\hat\tau_U),10^{-6},2-10^{-6})$.
4. Compute a joint geometric median on $\mathrm{Scaled}(\mathbb R^p,\alpha)\times\mathrm{Scaled}(Gr,2-\alpha)$. α scales every distance and norm. Under the default `Spread` unit the solver's thresholds scale with them, so α moves only the joint median itself, through the balance of the two factors.

**Rescaling the data.** The calibration keeps $\hat\alpha\hat\tau_\mu=(2-\hat\alpha)\hat\tau_U$, and rescaling the data by s multiplies $\hat\tau_\mu$ by s² while $\hat\tau_U$ (in radians) stays put, so the two calibrated metrics differ by a constant factor and MoMPCA is scale-equivariant, as long as the clamp $[10^{-6},2-10^{-6}]$ does not bind. It binds once $\hat\tau_U/\hat\tau_\mu$ leaves about $[5\times10^{-7},2\times10^6]$: on a fixture with $\hat\alpha\approx0.15$, s = 1e-3 and 1e2 reproduce the scale-1 result to 1e-14, while s = 1e-6 and 1e3 move the mean by about 1% and the subspace by about 0.01 rad.

**Output (`PcaResult`).**

| Field | Contents |
|---|---|
| `Mean` | the robust mean |
| `Components` | the r median-frame columns: an arbitrary orthonormal basis of the median subspace |
| `Eigenvalues` | per-column variances over the full data (not eigenvalues, and possibly unsorted) |
| `ExplainedVarianceRatio` | relative to the sum of those variances |
| `Whitened` | false |

### Subspace annealing (`SubspaceAnnealer`)

`SubspaceAnnealer` minimizes a black-box objective $f(P)$ over k-dimensional subspaces, where $P$ is a k×d orthonormal projection with **row** vectors.

- **Start.** The PCA subspace.
- **Temperature.** Cooling is geometric: $T=T_0c^{\,\text{iter}}$.
- **Moves.** Each iteration draws one move type:
  - **Isotropic** (Grassmann Exp of a horizontal tangent with an isotropic Gaussian direction and Frobenius length exactly s);
  - **Paired** (a random in-span mixture $m=\cos\varphi\,y_a+\sin\varphi\,y_{a+1}$ of two eigen-adjacent columns is rotated toward a random direction orthogonal to the span, while the partner $-\sin\varphi\,y_a+\cos\varphi\,y_{a+1}$ stays fixed; requires $k\ge2$);
  - **Single-column Givens rotation** by angle $\theta=s\,\mathcal N(0,1)$ (otherwise). This equals the Grassmann Exp of a rank-1 horizontal tangent.
- **Acceptance.** Metropolis: $f'<f$ or $U<e^{(f-f')/T}$. The best accepted value is tracked.
- **Step control.** Each move coordinate has its own acceptance-adaptive step size (one per retained column for single-column Givens moves, one for isotropic moves and one for paired moves, $k+2$ in all, each clamped to [`StepFloor`, `StepCeiling`]): $s\leftarrow s\,e^{g(1-p^*)}$ on acceptance and $s\,e^{-gp^*}$ on rejection, with $g=0.1$ and target acceptance $p^*=0.25$.

SPRED plugs a persistent-homology objective into this engine; see [Spred](./Spred.md).

### mxPBF two-sample tests (partial)

These compute the log of the **maximum pairwise Bayes factor** (Lee, You & Lin) for mean differences, per coordinate, or for covariance differences, per ordered coordinate pair. The multiplicity is absorbed by a prior $\gamma=\max(n,p)^{-\alpha}$. They return only the statistic, with no threshold or decision. The covariance prior constants are marked provisional in the code.

## Configuration

| Setting | Default | Meaning |
|---|---|---|
| `IrlsOptions.MaxIterations` | 200 | Iteration cap. **0 is a sentinel** that replaces the whole struct with `Default` |
| `IrlsOptions.Scale` | `Spread` | The unit λ of the four constants below: the spread σ, or 1 under `Data`. Also the enum's zero value |
| `IrlsOptions.Tolerance` / `Epsilon` | 1e-8 / 1e-10 | Shift tolerance; coincidence radius (in units of λ) |
| `IrlsOptions.HybridMode` / `SubgradientThreshold` / `Eta0` | `Hybrid` / 1e-4 / 1.0 | Switch to subgradient near data (threshold and step in units of λ) |
| `IrlsOptions.SingularityPolicy` / `ConvergenceCriterion` | `OptimalityCheck` / `Absolute` | In `Default`. The enum's zero value is `Regularise` |
| `EuclideanMedian` maxIterations / tolerance / epsilon / scale | 200 / 1e-8 / 1e-10 / `Spread` | tolerance and epsilon in units of λ |
| `consistencyFactor` | 1.0 | Scatter multiplier |
| Bandwidth reference population / grid / bisections | 32 / 1024 / 48 | — |
| `SubspaceAnnealer` maxIters / seed | 1000 / null | A null seed uses OS entropy |
| `SubspaceAnnealerOptions` IsotropicFraction / PairedFraction / TargetAcceptance | 0.1 / 0.0 / 0.25 | Move mixture; step-controller target |
| `SubspaceAnnealerOptions` InitialStep / StepFloor / StepCeiling | 0.1 / 1e-3 / π/2 | Geodesic step |
| `SubspaceAnnealerOptions` InitialTemperature / CoolingRate | 1.0 / 0.99 | Geometric cooling |
| `MxPbf` alpha / a0 / b0 | 2 / 1 / 1 | a0 and b0 are provisional |

## Status and caveats

- **`SingularityPolicy.Regularise` does not do what its documentation says.** It is documented as flooring distances at ε. In practice, coincident points already get weight 0, so the floor never applies, and Regularise effectively means "drop coincident points".
- **Partially set `IrlsOptions`.** If you set only some fields and leave `MaxIterations` nonzero, the others stay at zero, and `Scale` at `Spread`. With `Tolerance = 0`, `Absolute` convergence never triggers (under `Spread` the stall stop still ends the run). Start from `IrlsOptions.Default` instead.
- **`StiefelManifold` is partial.** Its geometry depends only on the column span, not on orientation, and its Exp is a retraction.
- **`MxPbf` is provisional.** The covariance test's prior constants `a0`/`b0` are provisional. Its pooled Gram is also centered incorrectly: it subtracts $n\bar z\bar z^\top$ from already group-centered Grams instead of adding the between-group terms $n_1(\bar x-\bar z)(\bar x-\bar z)^\top+n_2(\bar y-\bar z)(\bar y-\bar z)^\top$, so `TwoSampleCovariance` is biased unless both group means are zero.
- **Precision near zero distance.** Sphere and Grassmann distances use arccos, which resolves only about 1e-8 near zero. Under `Spread` the solver reads its radii from the Log map, which in turn drops principal angles below about 1e-8, and the Exp maps drop tangent steps below 1e-12; together they bound the median's relative accuracy by about 1e-12/s and its range at about s = 1e-6 ([measured](#the-unit-of-the-constants-irlsoptionsscale)). Under `Data` the 1e-8 tolerance sits at the arccos floor.
- **`MoMPCA` is scale-equivariant only inside its α̂ clamp** ([above](#median-of-means-pca-mompca)).

## Source map

| Role | Files (`ThermoMapper/src/maths/geometry/`) |
|---|---|
| Contract and manifolds | `IRiemannianManifold.cs`, `manifolds/EuclideanVectorManifold.cs`, `SphericalManifold.cs`, `PoincareBallManifold.cs`, `LorentzHyperboloidManifold.cs`, `GrassmannManifold.cs`, `StiefelManifold.cs`, `RiemannianProductManifold.cs`, `ScaledManifold.cs` |
| Solver | `solver/irls.cs`, `irls-options.cs`, `IRobustLoss.cs` |
| Estimators | `estimators/intrinsic/ManifoldMean.cs`, `ManifoldMedian.cs`, `KarcherScatter.cs`, `WeiszfeldScatter.cs`, `ScatterAccumulator.cs`; `estimators/ambient/EuclideanMean.cs`, `EuclideanMedian.cs`; `estimators/calibration/ConsistencyFactors.cs` |
| Bandwidth | `bandwidth/HeatKernel.cs`, `Hyperbolic.cs`, `Spherical.cs` |
| Subspaces | `dim-reduction/MoMPCA.cs`, `SubspaceAnnealer.cs`, `SubspaceAnnealerOptions.cs` |
| Inference | `inference/MxPbf.cs` |

**See also:** [DistanceMetrics](./DistanceMetrics.md) (geodesic kernels) · [GraphConstruction](./GraphConstruction.md#bandwidth-estimation) · [Spred](./Spred.md) · [DesignPrinciples](./DesignPrinciples.md)
