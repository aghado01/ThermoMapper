# Active subspaces

> E8 on the BASS core: the matrix C = E_ρ[∇f ∇fᵀ] of a fitted surface under an input law ρ, per posterior draw, and its top eigenvectors, the directions along which the surface changes most on average. This page grows slice by slice with [the build brief](../../issues/regression/bass/active-subspace-brief.md). Slice 1 has landed: the moment service's derivative moments and C per draw on the exact route.

## Overview

Rumsey et al. (*Discovering Active Subspaces for High-Dimensional Computer Models*, arXiv 2307.11241) compute C for a BASS fit in closed form. Each entry reduces to one-variable integrals of the terms' factors and their derivatives, and C is computed for each posterior draw. Their code is `coactivity`'s `C_bass` (pinned at 35e4573), the oracle. Slice 0 confirmed five defects in it, corrected on the engine's side ([Corrections to the sources](./BayesianRegression.md#corrections-to-the-sources)).

C finds linear combinations of the inputs along which f moves, which Sobol indices, ranking inputs and their subsets, cannot. Its uses are reducing a calibration's parameter space, choosing directions to refit along, comparing two models' responses (co-active subspaces, slice 3), and reading which output patterns move with the inputs.

| Piece | Where | Slice |
|---|---|---|
| Reference functions with C in closed form, a quadrature for any gradient | `tests/regression/ActiveSubspaceReferences.cs` | 0 |
| The oracle and its defects | `r/oracles/active_subspace_oracle.R` | 0 |
| A factor's slope moments E[φ′], E[φ′ψ], E[φ′ψ′] | `DesignFactor.SlopeMean`, `SlopePair` | 1 |
| A term pair's gradient moments, C for a linear combination of terms | `BassMoments.GradientPair`, `GradientOuterProduct` | 1 |
| The readout: both routes, the metric, E[C] and the subspaces' mean | `BassFit.ActiveSubspace` | 2 (not built) |

## C on the exact route (slice 1)

A term is B = κ·Πᵥ φᵥ(xᵥ), and its partial derivative in a continuous predictor i differentiates the one factor on i: ∂ᵢB = κ·φ′ᵢ·Π_{v≠i} φᵥ. Under a product law, E[∂ᵢB_a·∂ⱼB_b] is therefore a product over the variables of either term:

- E[φ′_{a,i} φ′_{b,i}] at v = i = j;
- for i ≠ j, E[φ′_{a,i} φ_{b,i}] at v = i and E[φ_{a,j} φ′_{b,j}] at v = j, each the slope's mean E[φ′] where the other term lacks the variable;
- E[φ_{a,v} φ_{b,v}] on every other variable both terms hold, and the one factor's mean on a variable only one holds.

Then C = Σ_{m,n} b_m b_n·E[∇B_m ∇B_nᵀ] for a draw's coefficients b, the intercept left out. These are the paper's I₁, I₂ and I₃ (its eq. Cij). The engine takes them at every hinge degree and for Legendre factors, under every law the value moments take: uniform, truncated normal and Student mixtures, and empirical marginals.

- **The factors.** A hinge's slope q·s·(s(x − t))₊^{q−1} is a piecewise polynomial of degree q − 1, integrated as the value moments are (`PiecewisePolynomial.HingeProduct` with slope flags), piece by piece through the law's interval moments. A Legendre factor's slope differentiates its expansion about each piece's centre. That expansion's coefficients grow with the degree, so digits are lost as they are for the value moments: ≤ 4·10⁻¹⁴ of the integrand's size to degree six, ≤ 7·10⁻¹² to degree ten.
- **The kink convention under an empirical marginal.** At degree one a slope jumps at its knot, and the empirical law's atoms sit on the knots. The pieces close on the right, so they give such an atom the left limit; the engine instead takes the derivative's mean of one-sided limits, s/2, which is what the gradient readout reads there (F-X3, `Hinge.Derivative`). The correction is the atom's mass times the difference (G3's "C reads the kink convention"). Under a continuous law the knots carry no mass, so the question does not arise.
- **Categorical inputs** take no part in C (decision 3): a level has no derivative. A term's categorical factors enter as cross moments under the level probabilities, so each continuous input's gradient is taken within its level and averaged over the levels. `C_bass` reads no categorical factor, so for it every level indicator is one.
- **Units.** C is in the fit's transformed units, the box scaled to [0, 1], the default metric of decision 2. For x = Lz, C_z = LᵀC_xL (Rumsey et al.'s Theorem 1), so the raw units' C is D⁻¹CD⁻¹ with D the boxes' widths. Slice 2's readout states its metric.
- **Caching.** Slope means and slope pairs are cached by factor, and gradient moments by term pair, as the value moments are. A draw that repeats a structure costs one pass over its term pairs: 7 to 9 µs for the fixture's 8-term draws with the caches warm, 21 to 44 µs cold (indicative, thread cycles, `ActiveSubspaceCostTests`).

### What was checked (slice 1)

| Gate | Result |
|---|---|
| Each factor pair against one-variable quadrature of the factors' own values and derivatives: hinges of degree one to three, every sign pair, knots inside, at the ends of and outside the support, under a uniform law, a normal–Student mixture and an empirical marginal with atoms on the knots | ≤ 5·10⁻¹⁵ of the integrand's size |
| Legendre factors to degree ten under the same laws, and on [0, 1] against ∫ψ′_mψ_n and ∫ψ′_mψ′_n in closed form | ≤ 4·10⁻¹⁴ to degree six, ≤ 7·10⁻¹² to degree ten |
| Seeded sweep, 400 hinge pairs: knots anywhere in [−1, 2], mixtures with a component up to 40 scales outside the interval, empirical marginals with knots on atoms | ≤ 3·10⁻¹⁴ |
| Each term pair, hinge products with and without categorical factors, degrees one to three, against a tensor quadrature of the basis's own derivative columns | ≤ 5.1·10⁻¹⁴ |
| C per draw of fitted hinge models (degrees one and two, with a categorical input) and of a Legendre fit, against the gradient readout averaged over the law by a rule exact for the integrand, under the uniform law, a mixture and an empirical marginal | ≤ 8.1·10⁻¹⁴ of C's largest entry |
| SG-4(b): E[x(1 − x)ψ′_mψ′_n] = n(n + 1)·δ_mn through the slope moments under the Jacobi weight; on a Legendre fit's draws, E[x_v(1 − x_v)(∂_v f)²] assembled from them against the coefficient readout | ≤ 7·10⁻¹⁴ to degree six (3·10⁻¹¹ to ten); 5·10⁻¹⁴ on the draws |
| C per draw of an engine fit against `C_bass` on the same draws (`ActiveSubspaceParityTests`) | uniform law: as written, 9·10⁻¹⁶; a truncated normal: corrected, 1.1·10⁻¹⁵, and 6.7% from it as written |

## The decisions (the brief's, taken in commissioning)

1. **Two routes, one readout, chosen per query.** The exact route serves a product law over separable terms (built in slice 1). The sampled route averages the gradient readout over the law's points otherwise: dependent inputs, non-separable kernels, or a check. The readout says which route ran. *Alternatives:* sampled always, which is simple but never exact; exact only, which excludes dependent inputs and non-separable kernels.
2. **The input metric is declared.** By default the box scaling, as above. *Options:* the inputs' own units, or whitening by the law's covariance, which is affine-equivariant and natural under a dependent law. Under a product uniform law, box scaling and whitening agree up to a constant.
3. **Categorical inputs take no part in C,** as above.
4. **Two posterior summaries, reported apart:** E[C] with its eigendecomposition, and the subspaces' projection mean P̄ = Σ w_s W_s W_sᵀ, whose eigenvalues are mean squared cosines (G3, SG-5). Slice 2.
5. **Numerics:** each draw's subspace from the thin SVD of its gradient matrix on the sampled route; the draws share the law's points. Slice 2.
6. **Per-draw consistency checks,** reported and not enforced: Poincaré's total(v) ≤ ((b − a)²/π²)·C_vv on a uniform input, and SG-4(b)'s weighted bound on a Legendre fit. Slice 2.

## Notation

This page takes the module's symbols: B a term, κ its scale, φ a factor and φ′ its derivative, b a draw's coefficients, q a hinge's degree. The paper writes h_{im} for a factor, γ_m for a coefficient and I₁, I₂, I₃ for the factor integrals. The code calls the factor integrals of derivatives "slope" moments (`SlopeMean`, `SlopePair`) beside the pointwise `Derivative`. C is E[∇f∇fᵀ] throughout, which the code names `GradientOuterProduct`, and `BassGradient` is the pointwise readout. These crossings are in [the notation ledger](../../issues/notation-ledger.md).
