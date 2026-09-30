# Synthetic Datasets

> Generators for labeled point clouds whose ground truth is known exactly. Each one isolates a single structural feature:
>
> - density voids
> - anisotropy
> - multi-scale hierarchy
> - non-convexity
> - non-orientable or self-intersecting topology
> - sparse supports
> - simplex-valued data
> - curved native geometry
>
> They are used to test metrics, proximity graphs and clusterers, SPC in particular, on controlled structure.
>
> **Project:** `Synthetic` (namespaces `Synthetic`, `Synthetic.Euclidean`, `Synthetic.Manifolds`)
>
> **Maturity:** core, except `GaussianManifold`, which is partial.

## Overview

`Synthetic` sits in layer 4 of the [architecture](./Architecture.md). It depends on:

- `Maths.Rng` and `Maths.Samplers`, for the eye generator's colored-noise background;
- `Maths.Geometry`, for the Poincaré-ball manifold (see [Riemannian Geometry](./RiemannianGeometry.md#manifold-catalog)).

No other in-scope library consumes it. Its consumers are the out-of-scope REPL and visualization modules, and the tests.

The source says that no generator prescribes a configuration. A generator describes a *shape*, and the caller decides which metric, graph and clusterer to test on it. Metadata tags such as `SuggestedMetric` and `FutureMetric` record the author's intent.

## Key concepts

- **Dataset.** A `SyntheticDataset` contains:
  - `Features`: a jagged `double[n][d]` array;
  - `Labels`: the finest ground-truth level;
  - `ClusterCount`;
  - `LabelsByLevel`: one label array per hierarchy level, ordered coarse → fine;
  - `Parameters`: an echo of the generator's inputs;
  - `Metadata`;
  - optionally, `ClusterGeometries`.
- **Contiguous label blocks.** Points are emitted cluster by cluster and never shuffled. Permute them yourself if you need i.i.d. order.
- **Background.** Some datasets include unstructured points:
  - The eye generators give background the label `ClusterCount`, so background is *excluded* from the count.
  - `Bwd1995Toy` labels its background 0 and *includes* it in `ClusterCount = 4`.
- **Metadata** (`SyntheticDatasetMeta`):
  - `GeneratorName`;
  - `GeometryClass`: "Euclidean", or "Manifold" for `Simplex`, `GaussianManifold` and the three hyperbolic generators;
  - `TopologyTag` and `HierarchyTag`;
  - `GTNumClusters`, which always equals `ClusterCount`;
  - `AmbientDimensionality`;
  - `LiteratureReference`, `SuggestedMetric` and `FutureMetric`, which are free-form strings rather than enums.
- **Cluster geometry.** Analytic descriptions of clusters, used for visualization and as oracles:
  - `EllipsoidGeometry(Center, Covariance)`;
  - `ArcGeometry(SpineSamples, Radius, NoiseScale)`;
  - `MobiusTubeGeometry`, with spine samples and [T, N, B] frames;
  - `ManifoldGeometry`, which is declared but not produced by any generator.

  `ClusterGeometry` is an empty base class that consumers downcast.
- **Eye skeleton.** A coordinate-free description of an "eye": strokes are geodesic arcs that carry tubes, plus an optional pupil. A skeleton can be realized in flat space or on the Poincaré ball.

## Design

- **Pure static generators returning a mutable DTO.** Each `Generate(...)`:
  1. creates one `System.Random(seed)`, with default seed 42;
  2. draws all its randomness from it in a fixed order;
  3. writes points into pre-sized arrays;
  4. fills `Parameters`, `Metadata` and, optionally, geometries.

  There is no registry or interface. A new generator is simply another static class.
- **Reproducibility contract.**
  - Identical parameters give identical output.
  - A parameter change that alters the number of random draws shifts every point generated after it.
  - Generators share no state, so concurrent calls are safe.

  This module uses `System.Random` rather than the library's [`SeedTree`](./DesignPrinciples.md#6-reproducible-randomness) discipline. The exception is the eye's noise field, which comes from a separate Xoshiro stream (seed XOR `0x9E3779B9`), so the field depends only on the seed and the noise settings. Structure points are drawn before any background draw, so background settings never move them. Background candidates still come from the main `System.Random`, so changing structure settings does move the background.
- **Shared toolkit.** An internal `SyntheticData` class provides:
  - normal, Gamma and Dirichlet samplers;
  - placement of centroids on a sphere;
  - random orthogonal matrices;
  - covariance assembly;
  - Euler rotations;
  - a regularized Cholesky factorization.
- **Skeleton → realizer (the eye family).**

  ```mermaid
  flowchart LR
    cfg[EyeTorusToyConfig] --> BS[BuildSkeleton] --> SK[EyeSkeleton]
    SK --> SL[SampleLocal: tangent vectors + fine/coarse labels]
    SL --> ER[flat realizer: c + s·v, radial cap] --> BG1[colored-noise background] --> DS1[SyntheticDataset]
    SL --> HR[Poincaré realizer: exp₀ of rescaled v] --> BG2[sinh²ρ background] --> DS2[SyntheticDataset]
  ```

  Both realizers seed identically and call `SampleLocal` first. Structure point i is therefore the *same* tangent sample in the flat and the hyperbolic dataset, and only the realization differs. This is how the family separates ambient from intrinsic clustering.
- **Composition by enums.** In the two-cluster topology generators, presets select independent aspects:
  - `EllipsoidPlacement` and `MobiusPlacement`: center and rotation;
  - `EllipsoidShellMode`: radial profile;
  - `TubeCrossSection`: tube cross-section;
  - `MobiusSpineShape`: spine.

  The eye family has its own per-stroke `CrossSectionShape` (`Shell`, `Solid`, `Ribbon`).

## API at a glance

Every entry point is `public static SyntheticDataset Generate(...)`, and all parameters have defaults (see [Catalog](#catalog)).

| Namespace | Generators and helpers |
|---|---|
| `Synthetic.Euclidean` | `SpatialBlobs`, `AnisotropicGaussian`, `BlattThreeCluster`, `BlattHierarchy`, `Bwd1995Toy`, `TwoMoons`, `SparseSupports`, `CrescentAndEllipsoid`, `MobiusAndEllipsoid` (+ `Project4DTo3D`), `EyeTorusToy` (+ `BuildSkeleton`, `SampleLocal`, `FlattenToPlane`, nested `EyeTorusToyConfig`); eye types `EyeStroke`, `EyePupil`, `EyeSkeleton`, `CrossSectionShape`; `MobiusTubeGeometry`; enums `EllipsoidPlacement`, `EllipsoidShellMode`, `MobiusPlacement`, `MobiusSpineShape`, `TubeCrossSection` |
| `Synthetic.Manifolds` | `HyperbolicBlobs`, `HyperbolicHierarchy`, `HyperbolicEyeTorus` (+ `DefaultMaxGeodesicRadius`), `Simplex`, `GaussianManifold` |
| `Synthetic` | `SyntheticDataset`, `SyntheticDatasetMeta`, `ClusterGeometry`, `EllipsoidGeometry`, `ArcGeometry`, `ManifoldGeometry` |

## Algorithms

### Catalog

Defaults are shown. n is the total number of points at those defaults.

| Generator | Native space (default dim) | Ground truth | Construction | Intended use |
|---|---|---|---|---|
| `SpatialBlobs` | $\mathbb R^d$ (2) | 4 clusters, n = 6000 | $x=\mu_c+\sigma z$, centroids on a sphere of radius `separation` = 5, σ = 0.5 | easy baseline |
| `AnisotropicGaussian` | $\mathbb R^d$, d ≥ 2 (2) | 3, n = 6000 | $x=\mu_c+R_c(s\odot z)$, $s=(a,1,\dots,1)$, a = 5, fresh rotation per cluster; means at radius 2 | covariance-aware metrics; heavy overlap |
| `BlattThreeCluster` | x–y plane (2) | 3, n = 6000 | Gaussians (σ = 1) at (−6, 0), (6, 0), (0, 6) | canonical SPC validation |
| `BlattHierarchy` | x–y plane (2) | 24 fine; levels [2, 6, 24], n = 6000 | coarse centers at radius 20; medium at 4 and fine at 0.8 around their parents, at evenly spaced angles plus jitter; leaf σ = 0.15 | temperature-dependent multi-scale resolution (Blatt–Wiseman–Domany 1997) |
| `Bwd1995Toy` | fixed 2-D | 3 stripes + background, n = 3200 | 800 uniform points on [−10, 10]×[−3, 3] plus three 14.5 × 0.8 stripes of 800 each; about 77 background points fall inside each stripe | SPC plateau test (BWD 1995 Fig. 1) |
| `TwoMoons` | $\mathbb R^d$ (2) | 2, n = 6000 | `make_moons` geometry with deterministic parameter spacing, noise 0.1 | non-convex shapes |
| `SparseSupports` | $\mathbb R^{500}$ | 4, n = 6000 | cluster c occupies features $[cs, cs+S)$ with stride $s=S-O$ (S = 50, O = 10); values in [0.5, 1.5), exactly 0 elsewhere | distance concentration; sparse supports |
| `CrescentAndEllipsoid` | $\mathbb R^3$ | 2, n = 7500 | tapered arc (R 3, width 0.4, half-angle 2.04 rad) + placed ellipsoid | false bridges in k-NN graphs |
| `MobiusAndEllipsoid` | $\mathbb R^3$ or $\mathbb R^4$ | 2, n = 7500 | twisted tube on a figure-eight or circle spine + ellipsoid | non-orientable, self-intersecting topology |
| `EyeTorusToy` | $\mathbb R^D$, D ≥ 3 (3) | 3 strokes (+ pupil) + background; levels [signal/background, fine], n = 11200 | skeleton → flat realizer, pink-noise background | ambient vs intrinsic clustering (toroidal lift of the BWD 1996 toy) |
| `HyperbolicEyeTorus` | Poincaré ball, D ≥ 3 | as the eye | skeleton → exp₀ realizer | curvature distortion |
| `HyperbolicBlobs` | Poincaré ball (2) | 3, n = 6000 | wrapped isotropic normals around exp₀ centers | Poincaré-distance clustering |
| `HyperbolicHierarchy` | Poincaré ball (2) | 27 leaves; levels [3, 9, 27], n = 6000 | recursive geodesic branching | hierarchy in hyperbolic space |
| `Simplex` | $\Delta^{K-1}$ (K = 20) | 5, n = 6000 | Dirichlet clusters with disjoint (or overlapping) supports, κ = 20 | KL/JSD saturation |
| `GaussianManifold` | (μ, log σ) plane (2) | 3, n = 6000 | isotropic blobs in (μ, log σ) coordinates | Fisher–Rao vs Euclidean (intent only) |

### Shared samplers

- **Normal.** The cosine branch of Box–Muller, $z=\sqrt{-2\ln u_1}\cos(2\pi u_2)$ with $u\in(0,1]$. Each normal consumes two uniforms.
- **Gamma(α, 1).** Marsaglia–Tsang rejection for α ≥ 1, with $d=\alpha-\tfrac13$ and $c=1/\sqrt{9d}$, a squeeze, then a log test. For α < 1 it boosts via $G_{\alpha+1}U^{1/\alpha}$.
- **Dirichlet.** Normalized Gamma draws.
- **Centroids on a sphere.** $\mu_c=r\,g/\lVert g\rVert$, independently for each cluster. r is a *radius*, not a guaranteed gap (see [Status](#status-and-caveats)).
- **Random orthogonal matrix.**
  - d = 2: a uniform rotation.
  - d ≥ 3: modified Gram–Schmidt on a Gaussian matrix. The result is Haar on O(d), so it may be a reflection. A degenerate column (norm below 1e-12) is left unnormalized; the divisor becomes 1.0.
- **Covariance and Cholesky.** $\Sigma=R\,\mathrm{diag}(s^2)R^\top$. The Cholesky factorization is regularized with floors of 1e-14.
- **Euler rotation.** $R=R_z R_y R_x$, i.e. extrinsic x → y → z.

### Topology stress tests

**Crescent.** Let A be `arcHalfAngle`.

- The arc angle is $t=\pi+2A\tilde u$, with $\tilde u=\mathrm{clamp}(0.35z,-\tfrac12,\tfrac12)$. Points concentrate at the elbow $(-R,0,0)$, and the open face points toward +x.
- The width tapers as $\cos(\pi\tilde u)$.
- About 15% of samples clamp to the tip angles, where the taper is zero. These form two needle-like concentrations.
- The default span is 2A ≈ 234°, more than a half annulus.

**Ellipsoid.** The ellipsoid is $x=c+L(r\hat g)$ with $LL^\top=Q\,\mathrm{diag}(a^2)Q^\top$, where a is the vector of semi-axes. The shell mode sets r:

| Mode | r | Effect |
|---|---|---|
| `Gaussian` | $\lVert g\rVert$ | exactly $N(c,\Sigma)$ |
| `Solid` | $U^{1/3}$ | uniform in volume |
| `Hollow` | 1 | a linearly warped sphere, *not* uniform in area |
| `Annular` | $0.75+0.25U$ | uniform in radius, not volume |

`EllipsoidGeometry.Covariance` always stores Σ, which equals the sample covariance only in `Gaussian` mode. The sample covariance is Σ/3 for `Hollow` and Σ/5 for `Solid`.

Placement presets (`NearOpenFace`, `OrthogonalElbowIntersect`, `IntersectUpperTip`, `IntersectLowerTip`, `Manual`) fix the center and the direction of the long axis. An explicit Euler input overrides the rotation.

**Möbius tube.** The twist half-angle is $h=k\theta/2$. With k = 1 the frame flips after one circuit, which gives a Möbius band; with even k the band is orientable.

- **Spines** (`MobiusSpineShape`).
  - `Circle`: $R(\cos\theta,\sin\theta,0)$, with θ density $\propto1+0.4|\cos\theta|$.
  - `FigureEight` (the default): $(R\sin\theta,\tfrac R2\sin2\theta,0)$. The spine crosses itself at the origin, at θ = 0 and θ = π, and its width narrows toward the crossing by the splay factor.
- **Cross-sections.** A ribbon, or a disk drawn from Gaussian, uniform-disk or annular radial profiles. Isotropic noise is then added.
- **4-D lift.** $x_4\propto(\text{in-section offset})\cdot\sin h$. It separates most of the crossing sheets. It does *not* lift points near the spine, nor the θ ≈ 0 branch at the figure-eight crossing, so it only partially resolves the self-intersection.
- **Projection back to 3-D.** `Project4DTo3D` rotates the (z, x₄) plane by α.

### The eye family

The skeleton has three strokes:

- a central full ring (label 0);
- upper and lower arcs (labels 1 and 2), spanning 216° and 288° by default, optionally tapered;
- an optional pupil (label 3).

Their sweep planes are:

- D = 3: all strokes share axes (0, 1, 2).
- D ≥ 4: central (0, 1, 2), upper (0, 2, 3) and lower (0, 3, 1). These planes share axis 0, so they are not independent.

**Tangent sampling.** For each stroke:

1. $\theta\sim U[\theta_0,\theta_1]$. The minor radius is scaled by $(1-\tau)+\tau\sin(\pi t)$, where $t$ is the arc fraction and τ = `HalfArcTaper` for the two arcs (default 0, i.e. untapered; always 0 for the central ring).
2. The cross-section angle φ is biased toward the inner edge by rejection with acceptance $(1-g)+g\,\tfrac{1-\cos\phi}{2}$. The mean acceptance is $1-g/2$.
3. The point is $v=(R+a)(\cos\theta\,e_A+\sin\theta\,e_B)+b\,e_{ax}+o+\sigma\xi$. Here $(a,b)=(r\cos\phi,\ z\,r\sin\phi)$, with $z$ = `ZThickness` (0.3) and $r$ the (tapered) minor radius for `Shell` (the default) or $r\sqrt U$ for `Solid`. `Ribbon` skips step 2 and draws $a\sim U[-r,r]$ and $b\sim N(0,\,\text{RibbonThickness}^2)$ (0.04).

**Flat realizer.**

- **Map:** $x=s\,v+c$, with an optional radial cap $\rho_{max}$.
- **Background count:** $\lfloor0.12\,N_{struct}\rfloor$ points.
- **Box:** centered at the origin, not at `Center`.
- **Density:** candidates are accepted with probability $e^{\kappa(f-\max f)}$ against a spectral noise field f (pink noise by default; see [Samplers](./Samplers.md)). This gives density $\propto e^{\kappa f}$ for `BackgroundContrast` κ ≥ 0. A negative κ silently yields a uniform background.
- **Guard:** after $1000N_{bg}+1000$ trials, the remaining background points are filled uniformly.
- **Flattening:** `FlattenToPlane` keeps the first two coordinates and sets `AmbientDimensionality` = 2.

### Poincaré-ball conventions

The hyperbolic generators return Euclidean coordinates inside the open unit ball, with curvature −1.

- **Exponential map.** `HyperbolicEyeTorus` and `HyperbolicHierarchy` use `PoincareBallManifold.ExpMap`, $\exp_p(v)=p\oplus\tanh(\lambda_p\lVert v\rVert/2)\,\hat v$ with $\lambda_p=2/(1-\lVert p\rVert^2)$, which clamps its result to radius $1-10^{-5}$. `HyperbolicBlobs` inlines the same origin map and its own Möbius addition, without the clamp. At the origin λ = 2, so a tangent norm η lands at geodesic radius 2η.
- **Möbius addition.**
  $$c\oplus y=\frac{(1+2\langle c,y\rangle+\lVert y\rVert^2)c+(1-\lVert c\rVert^2)y}{1+2\langle c,y\rangle+\lVert c\rVert^2\lVert y\rVert^2}.$$

**`HyperbolicEyeTorus`.** Let $s$ be the normalized radius of each tangent sample. Blend a "faithful" tangent norm $\eta_f=s\rho_{max}/2$ with a "cosmetic" one $\eta_c=\operatorname{artanh}(s\tanh\tfrac{\rho_{max}}2)$, using `WarpStrength` w.

- w = 1 puts geodesic radius exactly $s\rho_{max}$.
- w = 0 gives a linear Euclidean picture.
- $\rho_{max}$ defaults to 2.5, so every point lies within Euclidean radius 0.848.
- The background draws ρ with density $\propto\sinh^2\rho$. This is uniform in hyperbolic volume only for D = 3.

**`HyperbolicBlobs`.**

- Centers are at $\tanh(S)\hat g$, i.e. geodesic radius 2S, with S = `separation`.
- Points are $c\oplus\tanh(\lVert v\rVert)\hat v$. Since left Möbius translation is an isometry, each cluster is a wrapped isotropic normal.
- The library's Poincaré distance clamps inputs to radius $1-10^{-5}$. As a result, 2S is faithful only for S ≲ 6.1. Beyond that, a center's measured geodesic radius saturates near 12.2, and pairwise distances saturate near 24.4. They do not become infinite.

**`HyperbolicHierarchy`.**

- The root is at the origin. Children are placed by an exact geodesic step δ in a *uniformly random* direction, so there is no outward bias.
- Leaves are $\exp_c(v)$ with $v\sim N(0,(\ell^2+\nu^2)I)$ and $\ell=0.8\gamma^L$.
- To first order the Euclidean scatter is about $\lVert v\rVert$ at every depth. The geodesic scatter $(1+\cosh\rho_c)\lVert v\rVert$ grows sharply toward the boundary. At the defaults, for a center at geodesic radius 2, it is already about 1.9, comparable to the sibling separation.
- The defaults split 6000 → 3 × 2000 → 9 × ~667 → 27 × ~222 points.
- A leaf that stops early reuses its last label at the deeper levels, and that label can coincide with an unrelated node's label.

### Simplex and Gaussian-manifold data

**`Simplex`.** Concentrations are $\alpha_k=\kappa m_k+10^{-3}$.

- In disjoint mode, cluster c puts mass $1/q$ on its own q categories and $10^{-6}$ elsewhere.
- Off-support Gamma draws are tiny, and many underflow to exactly 0. This drives KL divergence toward saturation, as intended.
- Leftover categories (K mod C) belong to no support.

**`GaussianManifold`** places isotropic blobs in (μ, log σ) coordinates at alternating ±0.5 offsets. No Fisher–Rao geometry is computed: the manifold intent exists only as the `FutureMetric = "FisherRaoHalfPlane"` tag.

## Configuration

All generators take `seed = 42` and the defaults in the catalog. The exceptions are the two eye generators: `Generate(cfg = null, overrideSeed = null)` seeds from `EyeTorusToyConfig.Seed` (42) unless `overrideSeed` is given. The less obvious knobs:

| Setting | Default | Meaning |
|---|---|---|
| `CrescentAndEllipsoid` ellipsoid | 2500 points, axes (3, 0.5, 0.5), `Hollow`, `NearOpenFace` | Semi-axes, shell, preset |
| `MobiusAndEllipsoid` tube | R 2.5, half-width 1.1, half-thickness 0.12, σ 0.06, twist 1, `Annular`, `FigureEight`, splay 0.7 | Tube and spine |
| `MobiusAndEllipsoid` ellipsoid / dims | `Solid`, `CenterCrossOrtho` / 3 (3 or 4) | — |
| `gapScale` (both two-cluster generators) | 1 | Accepted and echoed but **unused** |
| `EyeTorusToyConfig` radii (major/minor) | central 2.5/0.6, upper 2.0/0.4, lower 3.8/0.9 | Stroke tori |
| `EyeTorusToyConfig` points | 4000 / 2500 / 3500; noise σ 0.08 | Per stroke |
| `EyeTorusToyConfig` background | ratio 0.12, pink exponent 1.0, contrast 1.0, grid 32 | Colored-noise blanket |
| `EyeTorusToyConfig` other | density gradient 0.4, scale 1.0, dimension 3, cap ∞, warp 1.0 | — |
| `HyperbolicEyeTorus.DefaultMaxGeodesicRadius` | 2.5 | Used when the cap is ∞ |
| `HyperbolicHierarchy` | depth 3, branches 3, step 1.0, decay 0.65, noise 0.22 | `basePointsPerLeaf` (250) acts only as a stopping threshold |
| `Simplex` | K = 20, disjoint, κ = 20 | — |

**Memory.** The eye's noise grid has `grid^D` cells, about $10^6$ in 4-D. `SparseSupports` is dense (6000 × 500).

## Status and caveats

- **Cluster separation is not guaranteed.** Centroids are placed independently on a sphere, so two clusters can land arbitrarily close. This affects `SpatialBlobs`, `AnisotropicGaussian`, the coarse level of `BlattHierarchy`, and `HyperbolicBlobs`. The `SpatialBlobs` comment promises gaps larger than the spread.
- **Background is counted inconsistently in `ClusterCount`.** `Bwd1995Toy` counts its background and the eye generators do not. Check `Metadata` before scoring.
- **`GaussianManifold` is partial.** It generates Euclidean blobs in parameter coordinates. The Fisher–Rao comparison it names is not implemented.
- **Some settings are inert, and some comments overstate the geometry.**
  - `gapScale` is unused.
  - `HyperbolicEyeTorus` ignores the `Center` coordinates, `GlobalScale` and the colored-noise settings.
  - `GaussianIsotropic` and `GaussianAnisotropic` Möbius cross-sections share one code path.
  - The comments overstate the Möbius 4-D lift and the uniformity of the `Annular` shell ("uniform in that band"), and say hierarchy branches "peel outward".
- **Placement names assume different spines.** With the default figure-eight spine, `NearSeam` sits at a lobe tip, not at a seam. With the circle spine, `CenterCrossOrtho` sits on the ring at θ = π and `CenterCrossCoPlanar` in the ring's empty center. Check the placement against the spine shape.

## Source map

| Role | Files (`ThermoMapper/src/synthetic/`) |
|---|---|
| Dataset types and shared samplers | `SyntheticData.cs` |
| Euclidean generators | `euclidean/SpatialBlobs.cs`, `AnisotropicGaussian.cs`, `BlattThreeCluster.cs`, `BlattHierarchy.cs`, `Bwd1995Toy.cs`, `TwoMoons.cs`, `SparseSupports.cs`, `CrescentEllipsoid.cs`, `MobiusEllipsoid.cs` |
| Eye family | `euclidean/EyeSkeleton.cs`, `euclidean/EyeTorusToy.cs`, `manifolds/HyperbolicEyeTorus.cs` |
| Manifold generators | `manifolds/HyperbolicBlobs.cs`, `HyperbolicHierarchy.cs`, `Simplex.cs`, `GaussianManifold.cs` |

**See also:** [RiemannianGeometry](./RiemannianGeometry.md) · [SpcEngine](./SpcEngine.md) · [DistanceMetrics](./DistanceMetrics.md) · [Archivory](./Archivory.md)
