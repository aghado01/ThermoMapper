# Active subspaces

> E8 on the BASS core: the matrix C = E_ρ[∇f ∇fᵀ] of a fitted surface under an input law ρ, per posterior draw, and its top eigenvectors, the directions along which the surface changes most on average. This page grows slice by slice with [the build brief](../../../issues/regression/bass/active-subspace-brief.md). Slices 1 to 5 have landed: the moment service's derivative moments and C per draw on the exact route; the readout with both routes, the declared metric, the posterior's two summaries and the consistency checks; functional outputs under a declared output metric, and co-active pairs; dependent inputs, a law given as a joint sample, with whitening by its full covariance; then the Grassmann discussion's settled items: per-draw ranks, the output side, the matrix Poincaré check, basis adequacy and stable angles.

## Overview

Rumsey et al. (*Discovering Active Subspaces for High-Dimensional Computer Models*, arXiv 2307.11241) compute C for a BASS fit in closed form. Each entry reduces to one-variable integrals of the terms' factors and their derivatives, and C is computed for each posterior draw. Their code is `coactivity`'s `C_bass` (pinned at 35e4573), the oracle. Slice 0 confirmed five defects in it, corrected on the engine's side ([Corrections to the sources](./Regression.md#corrections-to-the-sources)).

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
| Dependent inputs: a joint sample as the law, the sampled route at its points | `JointSample`, `ActiveSubspaceSpec.Sample`, `CoActiveSpec.Sample` | 4 |
| Whitening by a law's full covariance, the symmetric root; standardizing, its diagonal | `InputMetric.Whitened`, `InputMetric.Standardized`, `Maths.Regression.Subspaces.Whitening` | 4 |
| A basis fit's sampled route, product law or sample | `BassBasisFit.ActiveSubspace`, `CoActiveSubspace` | 4 |
| Per-draw k by a declared rule, mixed-rank P̄, the rank summary | `RankRule`, `ActiveSubspaceSpec.Rank`, `Maths.Regression.Subspaces.RankPosterior` | 5 |
| The output side G = E[J W_x Jᵀ] | `ActiveSubspaceSpec.OutputPatterns`, `Maths.Regression.Subspaces.OutputPosterior` | 5 |
| The matrix Poincaré check | `ActiveSubspaceAnalysis.MatrixPoincare`, `Maths.Regression.Subspaces.MatrixPoincare` | 5 |
| Basis adequacy beside a direct fit | `BassBasisFit.Adequacy`, `Maths.Regression.Subspaces.BasisAdequacy` | 5 |
| Principal angles through LinAlg's stable routine | `SubspacePosterior.PrincipalAngles` → `Maths.LinAlg.PrincipalAngles.Between` | 5 |

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
  Outlying draws show in their angles; nothing is down-weighted. Slice 2 took the angles as arccos of a thin SVD's cosines, good to about 1.5·10⁻⁸ rad (decision 5); from slice 5 they come from LinAlg's stable routine, good to about √(p·k)·ε (below).
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

## Dependent inputs (slice 4)

**A joint sample is the law.** `ActiveSubspaceSpec.Sample` (and `CoActiveSpec.Sample`) takes a `JointSample`: points in the inputs' own units, one column per predictor, categorical ones included, with optional weights, normalized. The law is the sample's own, as a calibration's posterior draws or an SMC run's weighted particles give it. The readout runs the sampled route at those points, every draw sharing them. Under a sample:
- the exact route, which needs a product law, is refused; so are `Measures` beside it; `Points` and `Seed` are not read;
- the points must lie in the boxes (the declared domains, or the training ranges), which is where a product law's readout integrates too; beyond them the fit extrapolates;
- the analysis has no product measure and no consistency check (`Measure` and `Poincare` are null): Poincaré's bound is a marginal's.

C is Σ_r w_r J_rᵀW_yJ_r over the points, and each entry's standard error comes from Σw²(v − m)²/(1 − Σw²), which with equal weights is the sample variance over N, slices 2 and 3's formula.

**Whitening.** `InputMetric.Whitened` is L = Σ_ρ^{1/2}, the symmetric root of the law's covariance in box units, so C_u = LᵀCL is the generalized problem C v = λ Σ_ρ⁻¹ v (decision 2).
- Under a product law Σ_ρ is diagonal and L = diag(σ): slice 2's readout, unchanged.
- Under a joint sample Σ_ρ is the sample's weighted covariance, the sample read as its own law (no N − 1), and L mixes the inputs.
- `InputMetric.Standardized` names the diagonal form, L = diag(σ) under either law. It keeps the input axes and is not affine-equivariant under a sample.

The fit-agnostic primitive `Whitening` (`Moments`, `SymmetricRoot`) refuses a covariance whose smallest eigenvalue is below 10⁻¹² of its largest: a degenerate law has no whitened coordinates.

**Equivariance.** Under any invertible affine map z = Ax + c of the inputs and their law, the whitened eigenvalues are invariant, and the whitened frames differ by the orthogonal Q = L_x⁻¹A⁻¹L_z. A whitened direction v is the input functional L⁻ᵀv, and those map to the same functions. The box and raw metrics' spectra move.

**Where the metric matters.** An exactly rank-k C's top-k subspace is the same set of input functionals in every invertible metric: LᵀCL has range LᵀG for G the span of the gradients, and L⁻ᵀ maps it back to G. A metric decides only among the directions a C of rank above k leaves.
- 2307.11241's triangle law (§4.1.1) has C near rank one, λ₁ carrying 99.96% of tr C. Its box and whitened leading functionals agree to 2·10⁻⁴ rad.
- On a band law with correlation 0.91, f = 16(x₁ − x₂)² + (x₁ + x₂)² changes fastest across the band, so the box metric leads with x₁ − x₂. The law barely moves that way, so whitening, which weighs each direction by how far the law moves along it, leads with x₁ + x₂. The two are perpendicular, and the readout reports both.

**A basis fit's sampled route,** product law or sample: per joint draw, each component's gradient readout at the shared points. With J = G·J_w, JᵀW_yJ = J_wᵀ(GᵀW_yG)J_w, so the output metric's root is taken in component space, at most K rows. Two outputs read a_fᵀJ = (Gᵀa_f)ᵀJ_w.

**Eigendecompositions at unit scale.** `DenseEigen`'s Jacobi sweeps stop rotating once every off-diagonal is below an absolute 10⁻¹², so for a matrix of small entries they leave off-diagonals as large as its smallest eigenvalues. E8 now decomposes every symmetric matrix divided by its largest entry, with a 10⁻¹⁴ threshold (`ScaledEigen`): each draw's C, P̄, V, W_y's root and whitening's root. The equivariance sweep found it, with QᵀQ − I at 3.8·10⁻³ at small scales, and the root primitive, with R² against Σ at 1.1·10⁻⁶. After the change they read 1.8·10⁻¹¹ and 6.8·10⁻¹⁵. LinAlg itself is unchanged; its other consumers are a follow-up.

### What was checked (slice 4)

| Gate | Result |
|---|---|
| Against a product law the exact route already gives: a joint sample drawn from it (a truncated normal, a uniform, uniform levels), on a scalar fit with a categorical input, 16,384 points | every entry of 12 draws' C within 2.45 standard errors |
| The same for a direct functional fit under a full W_y; a basis fit on its new sampled route, by the law's points and by a sample; two scalar fits' pairs; a basis fit's two outputs by the law's points and by a sample | within 1.73; 0.30 and 1.14; 1.81; 1.86 and 1.54 standard errors; no violation on the basis fit's checks |
| Whitening under a weighted triangle-law sample: the map against the covariance's root by Denman–Beavers iteration; each draw's C against L̂ᵀC_boxL̂, and standardizing against diag(σ̂); integer weights against repeated points | 9.6·10⁻¹⁵; 1.6·10⁻¹⁴; 1.1·10⁻¹⁴ |
| The whitened spectrum in the inputs' own units: on boxes [−3, 5] and [10, 12] under a correlated sample, against L_rawᵀC_rawL_raw | 4.9·10⁻¹⁵; the box and raw leading eigenvalues differ by up to 303% |
| Affine equivariance, a seeded sweep of 200: p to 5, 50 to 400 points (weighted in half), reference gradient fields, random maps with scales 10⁻³ to 10³ | whitened eigenvalues 4.3·10⁻¹³; QᵀQ − I 1.8·10⁻¹¹; C_u′ against QᵀC_uQ 6.2·10⁻¹⁴; the box spectrum moves by up to 1.7·10⁸ |
| The primitives: the root of random covariances with condition numbers to 10⁸ and scales 10⁻⁶ to 10⁶; weighted moments against repeated points; the refusals | R² against Σ 6.8·10⁻¹⁵, exactly symmetric; 6.3·10⁻¹⁸ |
| 2307.11241's triangle law: the paper's C_ρ₂ = [[1710, 741], [741, 322]]/540 against quadrature of the triangle; a BASS fit on 500 Latin hypercube points of the square, read under 20,000 points of the law | the closed form to 10⁻¹⁴; E[C] within 0.38% (box) and 0.36% (whitened, against L̂ᵀC_ρ₂L̂); leading direction 3·10⁻⁴ rad from the truth's, the two metrics' 1.7·10⁻⁴ rad apart |
| The band law, \|x₁ − x₂\| < 0.2, f = 16(x₁ − x₂)² + (x₁ + x₂)²: truths by quadrature of the band; the readouts under 20,000 points | the truths lead with x₁ − x₂ (box) and x₁ + x₂ (whitened); the readouts 0.036 and 0.051 rad from them, 1.48 rad apart |
| The refusals: a sample beside Measures, the exact route, points outside the boxes, bad weights, a missing column, whitening a degenerate sample, standardizing a constant input | each refused; the box metric reads the degenerate sample |

**Cost** (`ActiveSubspaceCostTests.JointSampleCostRow`, thread cycles, indicative): under a joint sample, 0.53 ms per draw at 4,096 points and 1.8 ms at 16,384, whitened alike; the basis fit's sampled route, 0.92 ms per joint draw at 4,096 points; whitening itself at 16,384 points, 0.23 ms for 2 inputs and 2.4 ms for 10.

**Slice 4's choices,** beyond the brief:
- *The symmetric root* for whitening. *Alternatives:* a Cholesky factor, which depends on the inputs' order; the eigenvector rotation (PCA whitening), which loses the input axes. All give the same eigenvalues.
- *The sample as its own law,* its covariance divided by the weights' sum. *Alternative:* N − 1, the unbiased estimate of a population's, which the sample is not here.
- *Points outside the boxes refused.* *Alternatives:* extrapolating with a count reported; truncating silently, which changes the law.
- *`Standardized` added;* `Whitened` keeps its name and, under a product law, its values.
- *`Measure` and `Poincare` nullable* in the analyses, null under a sample. *Alternative:* NaN-filled checks, which would read as checks that ran.
- *Unit-scale eigendecompositions in E8's code,* not a change to `DenseEigen`, whose other consumers would move.

## The Grassmann discussion's settled items (slice 5)

The user settled the [Grassmann discussion](../../../issues/regression/discussion/grassmann-stiefel-subspaces-brief.md)'s open items on 2026-11-08, and they form this slice.

**Per-draw k, opt-in.** `ActiveSubspaceSpec.Rank` takes a declared `RankRule`, with no default. Without one, the fixed-k readout is what it was, bit for bit, apart from the angles' new routine.
- *The rules.* A trace share: the smallest k whose leading eigenvalues carry the share of tr C_s. A gap: the k up to k_max with the largest (λ_k − λ_{k+1})/λ₁, with λ_{p+1} = 0 and k = p allowed. A gap against λ_k itself, or on a log scale, would read any trailing zero eigenvalue (an input no term holds) as the largest gap. Against λ₁ it favours the first gap when the spectrum decays geometrically: on two ridges with λ₂ = 0.37λ₁ it reads rank 1, where the trace share at 0.95 reads 2. A draw whose C is zero has k_s = 0.
- *The posterior over ranks* (`RankPosterior`): each k_s and their weighted distribution; P̄ = Σ w_s P_s over mixed ranks, whose trace is E[k_s] and whose ordered eigenvectors are the flag mean for mixed dimensions (Draper et al. 2014, cited from memory).
- *The rank summary.* Over projections of every rank, Σ_s w_s‖P − P_s‖²_F = tr P + E[k_s] − 2 tr(P P̄), which for P = UUᵀ is E[k_s] + Σ_i (1 − 2u_iᵀP̄u_i). It is least at P̄'s eigenvectors whose eigenvalue exceeds ½ (Ky Fan); on coordinate subspaces that is the median probability model (Barbieri and Berger 2004, cited from memory). It is reported beside the mean rank, not used as a default k.
- *Each draw's distance to the summary,* d² = ½‖P_s − P*‖²_F = |k_s − k*|/2 + Σ sin²θ_i over the min(k_s, k*) principal angles (Ye and Lim 2016, cited from memory), with the rank term and the containment defect reported apart.

**The output side.** `ActiveSubspaceSpec.OutputPatterns` (0 by default) reads G = E[J W_x Jᵀ] for a direct or basis fit: the m × m matrix of the output patterns that move with the inputs, which this page's Why names.
- W_x = LLᵀ is the input metric, and the output metric enters by its symmetric root S, so the patterns stay over the locations: G_W = S·G·S, and tr G_W = tr C under the same metrics.
- *Exact route:* G = S·H·K·Hᵀ·S with K_ab = Σ_ij (W_x)_ij E[∂ᵢB_a ∂ⱼB_b] (`BassMoments.GradientGram`). *Sampled route:* each point's S·J·L·√w side by side, m × Np. A basis fit's J is diag(y_s)K·J_w.
- G is never formed. Each draw's spectrum and leading patterns come from its factor's thin SVD, and E[G]'s from the thin SVD of the draws' compressed factors side by side (`OutputPosterior`). A pattern's sign is set so that its largest entry is positive.

**The matrix Poincaré check** (`ActiveSubspaceAnalysis.MatrixPoincare`). Under a product law, Efron–Stein and each input's Poincaré inequality give Cov_ρ(f) ≼ E[J D Jᵀ], with D the inputs' constants, (b − a)²/π² for a uniform input.
- Per draw it reads the largest generalized eigenvalue of the covariance against the bound, both under the output metric's root, from the moment route: Cov = H·Σ_B·Hᵀ with Σ_B the terms' covariances (`CovarianceGram`), and the bound H·K_D·Hᵀ.
- At most 1 holds; above 1 + 10⁻⁹ is a violation. It is decision 6's check in matrix form, and it makes "low gradient energy, low variance" exact per draw.
- *Opt-in* (`ActiveSubspaceSpec.MatrixPoincare`, off by default): it multiplies the exact readout's cost per draw by 1.3 to 3.1 on scalar, functional and basis fits alike. `MatrixPoincareStatus` says whether it was read, not asked, or not applicable: asked under a categorical input, a law that is not uniform, or a joint sample, which has no such bound. Slice 2's per-input check stays on.

**Basis adequacy** (`BassBasisFit.Adequacy`). Beside a direct-mode fit of the same cases, it reads the share of the direct fit's G outside the basis's span, tr((I − P)G)/tr G with the span diag(y_s)K under S, and the containment defect of G's leading patterns in that span, per draw and for E[G]. A basis fit's own G lies in its span by construction, so the direct fit answers D07's question: did the truncation drop structure that moves with the inputs? E37 asks it from the truncation residuals instead.

**The angles.** `SubspacePosterior.PrincipalAngles` delegates to LinAlg's `PrincipalAngles.Between`, which pairs the cosines of aᵀb with the sines of the residual through atan2. Small angles keep their accuracy, about √(p·k)·ε, where slice 2's arccos resolved about 1.5·10⁻⁸. It reads bases of different ranks, as the mixed-rank distances need. Regression keeps no copy of it.

### What was checked (slice 5)

| Gate | Result |
|---|---|
| Without a rule: the fixed-k readout with and without ranks and output patterns asked, both routes | bit for bit, every draw's C, eigenvalues, basis and angles, E[C], P̄ and their eigensystems |
| The rules on chosen spectra, trailing zeros and k = p included; refusals | as stated |
| Sums of r ridges recover r: BASS fits of one ridge on [0, 1]³ and two on [0, 1]⁴ (an inert input), N = 900, the trace share at 0.95 | summary rank 1 and 2, every draw at r; at 0.99 one more (each fit puts a few percent of its gradient energy off the span); the gap rule reads the two ridges as 1 |
| tr P̄ against E[k_s]; the rank summary's loss against Σ w‖P* − P_s‖² computed directly; rank term plus containment against ½‖P_s − P*‖²; the summary against random projections of every rank and P̄'s leading k, a seeded sweep of 150 posteriors | 1.8·10⁻¹⁵; 1.2·10⁻¹⁴; 7.6·10⁻¹⁵; never beaten |
| tr G = tr C under the same metrics: direct and basis fits, box, raw and whitened, W_y = I and full, both routes, 96 draws | 5.5·10⁻¹⁴ |
| G against quadrature of the Jacobian: a basis fit's reconstructed J by each component's gradient readout, a direct fit's by its own; E[G]'s spectrum | 7.0·10⁻¹⁵ |
| The matrix check on (cos πx₁, cos πx₂), the extremal case; D halved; a mixed field; covariance outside the bound's range | 1 to 10⁻¹⁵; 2; 0.94; ∞ |
| The matrix check on fitted draws, asked: a scalar fit, a direct functional fit under a full W_y, a basis fit, 280 draws; by default; asked where it does not apply | no violation, worst 0.867; not asked; not applicable under a categorical input, a normal law and a sample |
| Adequacy on chosen factors: inside the span; F = [2e₁, 0.7e₄] against span(e₁, e₂); a non-orthonormal basis of the span | 0; b²/(a² + b²) to 10⁻¹⁵, the leading pattern contained, the second's defect 1; the same |
| Adequacy on fits: curves a(x)φ₁ + b(x)φ₂ + noise, a direct fit beside PCA bases of rank 1 and 2, against the truth's share for each fitted basis's own span | 0.2251 against 0.2251 (span(φ₁) itself: 0.2296); 0.0011 against 0.0001 |
| The angles against slice 2's arccos on 300 random pairs of mixed ranks, above 10⁻⁶ rad; rotations by 10⁻¹² to 10⁻⁶ rad | within 0.11 of arccos's own error; 2·10⁻¹⁶ relative, where arccos returns 0 below 10⁻⁸ |

**Cost** (`ActiveSubspaceCostTests.SettledItemsCostRow`, thread cycles, indicative), per draw at 3,000 draws:
- the exact readout with the matrix check and without it: scalar 20.4 and 8.6 µs, functional 35.9 and 22.9 µs, basis 59.2 and 19.1 µs (a third reading: 13.4 and 7.1, 28.9 and 15.6, 58.2 and 31.6; a fourth, on the wall clock with the final code: 16.6 and 12.9, 28.4 and 14.1, 65.2 and 30.0); per-draw ranks add less than the readings resolve (0.3 µs in the first row, −2.8 µs in the fourth);
- a functional fit's exact readout, 25 µs, and 43 µs with two output patterns;
- the sampled output side's own cost, 2.9 ms at 4,096 points, the thin SVD of the m × Np stacked Jacobians;
- basis adequacy, 0.16 ms per direct draw at 12 locations.

The matrix check is therefore off by default: it multiplies the exact readout's cost by 1.3 to 3.1 over four readings.

**Slice 5's choices,** beyond the brief:
- *The gap rule against λ₁.* *Alternatives:* against λ_k, or on a log scale, which read a trailing zero eigenvalue as the largest gap.
- *G's output metric by its symmetric root,* so patterns stay over the locations. *Alternative:* slice 3's eigen-root, whose rows are W_y's eigenvectors.
- *Output patterns opt-in,* the count given. *Alternative:* always computed, which costs every functional readout a thin SVD per draw.
- *The matrix check opt-in,* off by default, since it about doubles the exact readout on every fit kind: the user's standing rule that a payload which makes a path expensive becomes opt-in for that path and keeps its output. *Alternatives:* on by default, as slice 2's per-input check stays; caching its term-pair Grams per structure, which would cut its cost and is a follow-up.
- *The checks read from the moment route even on the sampled route,* as slice 2's do, so they carry no Monte Carlo error.

## The decisions (the brief's, taken in commissioning)

1. **Two routes, one readout, chosen per query.** The exact route serves a product law over separable terms (built in slice 1). The sampled route averages the gradient readout over the law's points otherwise: dependent inputs, non-separable kernels, or a check. The readout says which route ran (built in slice 2); a joint sample takes the sampled route (built in slice 4). *Alternatives:* sampled always, which is simple but never exact; exact only, which excludes dependent inputs and non-separable kernels.
2. **The input metric is declared.** By default the box scaling, as above. *Options:* the inputs' own units, or whitening by the law's covariance, which is affine-equivariant and natural under a dependent law. Under a product uniform law, box scaling and whitening agree up to a constant. Built in slice 2 for product laws; whitening by a joint sample's full covariance, and `Standardized` for the diagonal form, in slice 4.
3. **Categorical inputs take no part in C,** as above.
4. **Two posterior summaries, reported apart:** E[C] with its eigendecomposition, and the subspaces' projection mean P̄ = Σ w_s W_s W_sᵀ, whose eigenvalues are mean squared cosines (G3, SG-5). Built in slice 2. *Alternatives:* the Karcher mean or the geometric median of the subspaces (Geometry's IRLS), at most labelled options; a robust median is not the default, since an outlying draw is posterior mass.
5. **Numerics:** each draw's subspace from the thin SVD of its gradient matrix on the sampled route; the draws share the law's points. Built in slice 2.
6. **Per-draw consistency checks,** reported and not enforced: Poincaré's total(v) ≤ ((b − a)²/π²)·C_vv on a uniform input, and SG-4(b)'s weighted bound on a Legendre fit. Built in slice 2.

## Notation

This page takes the module's symbols: B a term, κ its scale, φ a factor and φ′ its derivative, b a draw's coefficients, q a hinge's degree; from slice 3, H a draw's coefficients at the output locations, W_y the output metric, a_f and a_g output functionals; from slice 4, Σ_ρ the law's covariance in box units and L = Σ_ρ^{1/2} whitening's map, x = μ + Lu; from slice 5, k_s a draw's own rank, k* the rank summary's, G = E[J W_x Jᵀ] the output side and S the output metric's symmetric root. Slice 3 wrote G for the reconstruction's map, which slice 5 writes out as diag(y_s)K. The paper writes h_{im} for a factor, γ_m for a coefficient and I₁, I₂, I₃ for the factor integrals. The code calls the factor integrals of derivatives "slope" moments (`SlopeMean`, `SlopePair`) beside the pointwise `Derivative`. C is E[∇f∇fᵀ] throughout, which the code names `GradientOuterProduct`, and `BassGradient` is the pointwise readout. These crossings are in [the notation ledger](../../../issues/notation-ledger.md).
