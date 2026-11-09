# Active subspaces

> E8 on the BASS core: the matrix C = E_ρ[∇f ∇fᵀ] of a fitted surface under an input law ρ, per posterior draw, and its top eigenvectors, the directions along which the surface changes most on average. This page grows slice by slice with [the build brief](../../issues/regression/bass/active-subspace-brief.md). Slices 1 to 3 have landed: the moment service's derivative moments and C per draw on the exact route; the readout with both routes, the declared metric, the posterior's two summaries and the consistency checks; then functional outputs under a declared output metric, and co-active pairs.

## Overview

Rumsey et al. (*Discovering Active Subspaces for High-Dimensional Computer Models*, arXiv 2307.11241) compute C for a BASS fit in closed form. Each entry reduces to one-variable integrals of the terms' factors and their derivatives, and C is computed for each posterior draw. Their code is `coactivity`'s `C_bass` (pinned at 35e4573), the oracle. Slice 0 confirmed five defects in it, corrected on the engine's side ([Corrections to the sources](./BayesianRegression.md#corrections-to-the-sources)).

C finds linear combinations of the inputs along which f moves, which Sobol indices, ranking inputs and their subsets, cannot. Its uses are reducing a calibration's parameter space, choosing directions to refit along, comparing two models' responses (co-active subspaces, slice 3), and reading which output patterns move with the inputs.

| Piece | Where | Slice |
|---|---|---|
| Reference functions with C in closed form, a quadrature for any gradient | `tests/regression/ActiveSubspaceReferences.cs` | 0 |
| The oracle and its defects | `r/oracles/active_subspace_oracle.R` | 0 |
| A factor's slope moments E[φ′], E[φ′ψ], E[φ′ψ′] | `DesignFactor.SlopeMean`, `SlopePair` | 1 |
| A term pair's gradient moments, C for a linear combination of terms | `BassMoments.GradientPair`, `GradientOuterProduct` | 1 |
| The readout: both routes, the metric, the checks | `BassFit.ActiveSubspace`, `ActiveSubspaceSpec`, `ActiveSubspaceAnalysis` | 2 |
| The posterior over subspaces from per-draw C, any engine's | `Maths.Regression.Subspaces.SubspacePosterior`, `DrawSubspace` | 2 |
| Functional outputs: C = E[JᵀW_yJ] under a declared W_y, for direct and basis fits | `OutputMetric`, `ActiveSubspaceSpec.OutputMetric`, `BassBasisFit.ActiveSubspace` | 3 |
| Co-active pairs C_fg = E[∇f∇gᵀ] and their posterior | `BassFit.CoActiveSubspace`, `BassBasisFit.CoActiveSubspace`, `CoActiveSpec`, `Maths.Regression.Subspaces.CoActivePosterior` | 3 |
| A pair of term lists' cross moments, Σ w_mn E[∇B_m∇B_nᵀ] | `BassMoments.GradientCross` | 3 |

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

## The readout (slice 2)

`BassFit.ActiveSubspace(data, spec)` reads a scalar fit's draws under `ActiveSubspaceSpec`: the predictors' laws (as the Sobol readout takes them), the metric, the dimension k, the route, the sampled route's points and seed, and the draws.

- **The routes.** `Auto` takes the exact route, which every term type and every product law built today allows; slice 4's joint samples will take the sampled one. The sampled route draws N points from the law once, shared by every draw, reads the gradient readout there in batches of draws, and keeps each draw's N × p gradient matrix in box units over √N as a factor F with C = FᵀF. Its subspace comes from F's thin SVD, so its gaps are not squared, and each entry of C carries its standard error.
- **The metric.** A map L from metric to box coordinates, C = LᵀC_boxL, reported as `Map`. Under a product law L is diagonal: 1 for the box, one over the box's width for raw units, and the law's standard deviation in box units, read from its interval moments, for whitening. Subspaces are read in the metric, so they move under any non-orthogonal change of units, as decision 2 says they must.
- **The summaries** (`SubspacePosterior`, decision 4), under the draw weights:
  - E[C] and its eigensystem: directions weighted by their gradient energy, a draw by its size;
  - the projection mean, the top-k eigenvectors of P̄ = Σ w_s W_s W_sᵀ, each draw counted once, with P̄'s eigenvalues as mean squared cosines and its eigenvectors in order;
  - each draw's eigenvalues, gap λ_k − λ_{k+1} and principal angles to the mean subspace;
  - the principal angles between E[C]'s top-k eigenspace and the mean subspace, which say how far apart the two summaries sit.
  Outlying draws show in their angles; nothing is down-weighted. Angles are arccos of a thin SVD's cosines, good to about 1.5·10⁻⁸ rad (decision 5).
- **The checks** (`PoincareCheck`, decision 6), per draw and continuous input:
  - *Poincaré:* total(v) against ((b − a)²/π²)·C_vv for an input under a uniform law on [a, b], with the totals from the moment route under the same law and draws;
  - *SG-4(b):* on a Legendre fit under the product uniform law on the boxes, total(v) against ½·E[x_v(1 − x_v)(∂_v f)²] from the coefficients;
  - a ratio above one beyond rounding, or beyond three Monte Carlo errors on the sampled route, counts as a violation, which is a defect.
- **Refused:** k outside 1 to the continuous predictors. Categorical predictors are named in `Excluded`. (A functional response was refused until slice 3.)

### What was checked (slice 2)

| Gate | Result |
|---|---|
| SG-5's three draws (λ, v) = (1, a), (1, a), (100, b): with a = (1, 0), b = (3/5, 4/5), 3E[C] = [[38, 48], [48, 64]] and 3P̄ = [[59/25, 12/25], [12/25, 16/25]]; with a and b at 1.2 rad, the summaries' leading directions | exact matrices; 37.99° and 54.29° apart, each the closed form's to 10⁻¹² |
| Weighted draws against replicated ones, in the summary and on a fit (`fit with { DrawWeights }` against `Draws` repeated) | 2·10⁻¹⁴, 2·10⁻¹⁵ |
| The exact and sampled routes on a hinge fit with a categorical input, 8,192 points, under the uniform law, a mixture and an empirical marginal | every entry of every draw's C within 1.3 standard errors |
| C under a declared change of metric on a fit whose box is [−3, 5]: raw units and whitening under a truncated normal, against LᵀC_boxL, with the law's sd by quadrature | 5·10⁻¹⁵ |
| Decision 6: Poincaré on hinge (degree 2) and Legendre draws, SG-4(b) on Legendre draws; on a degree-one Legendre fit, SG-4(b)'s equality | no violation (worst ratios 0.81, 0.21 and 0.38); equality to 3·10⁻¹⁶ on 88 draw-inputs |
| Recovery at N = 100, 300, 900 (G3) on functions from outside the model: a ridge sin(2.5wᵀx) on [0, 1]³ and two ridges on [0, 1]⁴ with one inert input, against their spans and C in closed form | the angles of E[C]'s subspace and of the mean subspace and E[C]'s error all fall with N (registered in `ComparativeSweep`); at 900, 0.010 rad on both functions, within the Davis–Kahan bound from E[C]'s own error at every N |

The recovery gate's first rank-two case, a polynomial of two linear forms, had a relative gap of 0.0045 between its second eigenvalue and zero. Its second direction was then identified only as well as Davis–Kahan allows, an angle near 0.1 rad with E[C] within 0.9%. The gate now uses two ridges of comparable gradient energy (relative gap 0.37); the lesson stands in the readout: the gaps say how well each subspace is identified.

**Cost** (`ActiveSubspaceCostTests.ReadoutCostRow`, thread cycles, indicative): the exact readout at 3,000 draws (60 structures repeated), 15 µs per draw with its checks and 4 µs for the summary alone; the sampled readout, 0.6 ms per draw at 4,096 points.

## Functional outputs and co-active pairs (slice 3)

**A draw as terms and H.** Every readout now reads a draw as its terms B_a and a matrix H, each term's coefficient at each output location: one row for a scalar fit; H[ℓ, a] = c_a·h_a(t_ℓ) for a direct functional fit, h_a the term's output-side factors at the training output location t_ℓ with their normalizer; and H[ℓ, (k, a)] = G_ℓk·c_ka for a basis fit, over its joint draws, with G = diag(y_s)K the reconstruction's map and c_ka component k's coefficients. The Jacobian of the outputs at x is then J = H·∇B, so every quantity below is a quadratic form in the term pairs' moments E[∇B_a∇B_bᵀ], which `BassMoments.GradientCross` sums for two term lists and a weight matrix.

**Functional C.** C = E[JᵀW_yJ] = Σ_ab M_ab·E[∇B_a∇B_bᵀ] with M = HᵀW_yH, and W_y the declared output metric (`OutputMetric`, F-X6's Q) over the training output locations, the identity by default. W_y says which output units and correlations count: e_ℓe_ℓᵀ reads location ℓ's own C, a diagonal of quadrature weights integrates C over the output grid, and C is linear in W_y. W_y must be m × m, finite, symmetric and positive semidefinite, on either route. A basis fit's C is Σ_kk′ (GᵀW_yG)_kk′·E[∇w_k∇w_k′ᵀ] in the same form, on the exact route only: its sampled route waits for slice 4's joint samples. The sampled route of a direct fit maps each point's Jacobian by W_y's root R (W_y = RᵀR, from its eigendecomposition), so each draw's factor has N·rank(W_y) rows. The readout reads at the training output locations, as the gradient readout and the conditioned Sobol form do; a basis fit exists only there.

**The checks** keep slice 2's form: the totals are Σ_ab M_ab·T_ab(v) over the term pairs' totals, which with W_y = RᵀR are the sums of the outputs Ry's totals, each meeting its own bound, so the sum meets C's. The scalar readout now takes its totals the same way, M = bbᵀ.

**Co-active pairs** (Rumsey, Hardy, Ahrens and Vander Wiel, arXiv 2311.18146): C_fg = E[∇f∇gᵀ], with C_f and C_g beside it, per draw pair and in the readout's metric (LᵀC_fgL). Each side is a scalar fit, or a functional fit read through an output functional, f = a_fᵀy and g = a_gᵀy (`CoActiveSpec.Output`, `OtherOutput`); coactivity reads a functional fit per location, a = e_ℓ.
- *Two fits* on the same inputs (the same predictors, boxes and training values) at one hinge degree pair their draws by position, each pair weighted by the product of the two fits' weights, which samples the posteriors' product. coactivity's `Cfg_bass` pairs by its `mcmc.use` rows the same way.
- *Two outputs of one fit*, direct or basis, pair each draw with itself: the joint posterior of the two outputs.
- *The exact route* reads C_fg as `GradientCross` of the two draws' terms with weights u_m v_n, u = H_fᵀa_f, v = H_gᵀa_g (a scalar side's u is its b). *The sampled route* averages both gradient readouts at shared points and gives every entry of C_f, C_g and C_fg its standard error (`CoActiveError`); a basis fit's waits for slice 4.

**The posterior** (`CoActivePosterior`, fit-agnostic): E[C_f], E[C_g] and E[C_fg] under the weights; V = (E[C_fg] + E[C_fg]ᵀ)/2, whose eigenvectors are the co-active directions, each maximizing E[(wᵀ∇f)(wᵀ∇g)] orthogonal to the earlier ones. V may be indefinite, so its eigenpairs are ordered by |λ|, as the paper advises for finding gaps. The contributions π_i = λ_i/√(E[t_f]E[t_g]), with t = tr C, sum to the concordance of the means E[t_fg]/√(E[t_f]E[t_g]); each pair's own concordance is reported too, since the concordance of the means is not the mean of the concordances. The co-activity scores are α_i(q) = Σ_{j≤q} λ_j w_ij², or with |λ_j| unsigned.

### What was checked (slice 3)

| Gate | Result |
|---|---|
| The brief's gate: a basis fit's C (2 components, 15 locations) against the direct computation on its reconstructed outputs, J = G·J_w by each component's gradient readout and Σ w·JᵀW_yJ by tensor quadrature split at the knots, under the identity and a full W_y | ≤ 9.4·10⁻¹⁶ |
| A direct functional fit's C (7 locations) against quadrature of the gradient readout's Jacobian, under the identity, trapezoid weights and a full W_y | ≤ 4.6·10⁻¹⁵ |
| W_y's readings: e_ℓe_ℓᵀ against that location's quadrature, linearity in W_y, no W_y as the identity | 5.4·10⁻¹⁶ |
| The direct fit's exact and sampled routes under a full W_y, 8,192 points, 12 draws | within 1.2 standard errors |
| 2311.18146's simple example, f₁ = x₁² + x₁x₂ and f₂ = f₁ + βx₂³ on [0, 1]²: C₁, C₂, C₁₂ by the polynomial closed form against the paper's matrices and against quadrature, at β = 0, ½, 3, −12; t, the concordance, π summing to it, ordering by \|λ\|, the scores | 10⁻¹⁴; the paper's formulas hold |
| The paper's quoted readings | 0.944 at β = ½ and 0.551 at β = 3 hold; three slips, below |
| The paper's estimate from BASS fits at β = 3, n = 200 and 1,000 Latin hypercube points: E[C₁₂]'s Frobenius distance from the closed form | 0.061 and 0.015 (0.87% and 0.22% of ‖C₁₂‖_F; the paper's single fits through `concordance`: 0.050 and 0.0086); the fall registered in `ComparativeSweep` |
| Two fits: swapping them transposes C_fg; a fit paired with itself gives its own C and concordance 1 | 3.8·10⁻¹⁶, 2.0·10⁻¹⁶ |
| Two fits against quadrature of both gradient readouts; exact against sampled | 6.1·10⁻¹⁵; within 2.1 standard errors |
| Two outputs of a direct fit: swap, a_f = a_g against C under W_y = a_fa_fᵀ, quadrature, the two-fit route with the fit paired with itself, the routes | 2.4·10⁻¹⁶, 6.0·10⁻¹⁶, 8.1·10⁻¹⁵, equal, within 1.8 standard errors |
| A functional fit through a_fᵀy paired with a scalar fit, against quadrature | 5.5·10⁻¹⁵ |
| Two outputs of a basis fit: a_f = a_g against C under a_fa_fᵀ; quadrature of the reconstructed a_fᵀJ, a_gᵀJ | 7.8·10⁻¹⁶, 1.1·10⁻¹⁵ |
| Two engine fits' pairs against `Cfg_bass` on the same draw pairs (`CoActiveParityTests`), and the mode's quadrature | uniform law: as written, 1.2·10⁻¹⁵; a truncated normal: corrected, 1.5·10⁻¹⁵, and 27% from it as written; quadrature 3.3·10⁻¹⁴ |
| A direct functional fit's C at each training location (W_y = e_ℓe_ℓᵀ) against coactivity's per-location reading, `bassfunc2bass_fixed_t` then `C_bass_univariate` with slice 0's corrections, and the mode's quadrature of the draw folded at t_ℓ: 4 draws at 4 of 9 locations | 3.5·10⁻¹⁵; 2.8·10⁻¹⁴ |

**The paper's slips** (2311.18146 v2, §3.1.1, each against its own definitions; the tests pin what each reads): the concordance −0.131 quoted for β = −12 is the value at β = −15 (at −12 it is −0.1095); the first co-active direction quoted at β = ½, (0.907, 0.422), is V's at β = 0.4 (at ½ it is (0.8971, 0.4418)); the contributions quoted at β = ½, (0.9518, −0.0077), are the unsymmetrized C₁₂'s eigenvalues over √(t₁t₂), where V's, which its eq. (contrib) defines, give (0.9684, −0.0242). Both pairs sum to the concordance, since tr V = tr C₁₂. Its script (`concordance`'s `inst/CoASM/simple_poly_concordance.R`, read and not run) symmetrizes as the text says.

**Cost** (`ActiveSubspaceCostTests.FunctionalAndCoActiveCostRow`, thread cycles, indicative), per draw at 3,000 draws: a direct functional fit's C over 7 locations, 39 µs under the identity and 44 µs under a full W_y; a basis fit's over 15 locations and 2 components, 43 µs; two outputs of the functional fit, 48 µs; per pair, two scalar fits, 15 µs exact and 0.93 ms sampled at 4,096 points. Slice 2's readout in the same run: 10 µs exact, 0.98 ms sampled. Against slice 2's own build, interleaved under one load, the sampled readout reads the same and the exact one about a quarter slower, its checks' totals now taken over every term pair through M.

**Slice 3's choices,** beyond the brief:
- *W_y is declared and checked,* the identity by default. *Alternatives:* diagonal weights only, which loses the outputs' correlations; a W_y built by the readout from the output grid, which would choose the units for the user.
- *Two fits pair draws by position* with product weights. *Alternatives:* every pair of draws, the exact product expectation at quadratic cost; random pairing, which adds its own noise.
- *A functional side is read through an output functional,* which holds coactivity's per-location reading (a = e_ℓ) and the integrated one. *Alternative:* per-location matrices only, as coactivity returns them.
- *The co-active posterior is E[C_f], E[C_g], E[C_fg] and V's eigensystem,* with each pair's concordance. *Not built:* a projection mean of the pairs' co-active subspaces, decision 4's second summary; V's eigenvalues are signed, so its top-q would need the magnitude order to be stable across pairs.
- *Not built:* a basis fit paired with another fit, and output locations beyond the training ones for a direct fit. Both fit the terms-and-H form, so each is a new H.

## The decisions (the brief's, taken in commissioning)

1. **Two routes, one readout, chosen per query.** The exact route serves a product law over separable terms (built in slice 1). The sampled route averages the gradient readout over the law's points otherwise: dependent inputs, non-separable kernels, or a check. The readout says which route ran (built in slice 2). *Alternatives:* sampled always, which is simple but never exact; exact only, which excludes dependent inputs and non-separable kernels.
2. **The input metric is declared.** By default the box scaling, as above. *Options:* the inputs' own units, or whitening by the law's covariance, which is affine-equivariant and natural under a dependent law. Under a product uniform law, box scaling and whitening agree up to a constant.
3. **Categorical inputs take no part in C,** as above.
4. **Two posterior summaries, reported apart:** E[C] with its eigendecomposition, and the subspaces' projection mean P̄ = Σ w_s W_s W_sᵀ, whose eigenvalues are mean squared cosines (G3, SG-5). Built in slice 2. *Alternatives:* the Karcher mean or the geometric median of the subspaces (Geometry's IRLS), at most labelled options; a robust median is not the default, since an outlying draw is posterior mass.
5. **Numerics:** each draw's subspace from the thin SVD of its gradient matrix on the sampled route; the draws share the law's points. Built in slice 2.
6. **Per-draw consistency checks,** reported and not enforced: Poincaré's total(v) ≤ ((b − a)²/π²)·C_vv on a uniform input, and SG-4(b)'s weighted bound on a Legendre fit. Built in slice 2.

## Notation

This page takes the module's symbols: B a term, κ its scale, φ a factor and φ′ its derivative, b a draw's coefficients, q a hinge's degree; from slice 3, H a draw's coefficients at the output locations, W_y the output metric, a_f and a_g output functionals. The paper writes h_{im} for a factor, γ_m for a coefficient and I₁, I₂, I₃ for the factor integrals. The code calls the factor integrals of derivatives "slope" moments (`SlopeMean`, `SlopePair`) beside the pointwise `Derivative`. C is E[∇f∇fᵀ] throughout, which the code names `GradientOuterProduct`, and `BassGradient` is the pointwise readout. These crossings are in [the notation ledger](../../issues/notation-ledger.md).
