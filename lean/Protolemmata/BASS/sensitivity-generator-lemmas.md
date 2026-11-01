# BASS — sensitivity as a Markov generator

Grounds the variance decomposition ([BassMoments](../../../docs/Maths/BayesianRegression.md)), active-subspace analysis (E8) and KHAOS's coefficient prior (E5) in one frame: Hoeffding components are the eigenfunctions of resampling dynamics. These are the statements behind gaps G3, G4 and G5 of the [lineage-gaps register](../../../issues/regression/bass/lineage-gaps-register.md).

**Status (2026-11-01).** Posed. [`lineage-gap-checks.py`](lineage-gap-checks.py) checks each item numerically, on the labels below. The Lean is uncompiled: Mathlib's build output is absent from `lean/.lake`, so names and modules were checked against the Mathlib source only.

**Prior art.** The parts are known separately: the Efron–Stein and Hoeffding decomposition, DGSM's Poincaré link to total indices (Lamboni et al., arXiv 1202.0943; Roustant et al., arXiv 1612.03689), eigenfunction lower bounds (arXiv 1906.09883), weighted Poincaré inequalities (arXiv 2412.04918, 2605.30979). The Boolean version of SG-1 is in O'Donnell's *Analysis of Boolean Functions* (from memory, to confirm). The searches found no sensitivity paper stating the whole frame. It may be folklore, so this is synthesis, not a claim of discovery.

---

## SG-1 — the Sobol decomposition is a spectrum

**Setting.** Finitely many variables v, each on a finite set with a probability vector μ_v; μ = ⊗μ_v; ⟨f, g⟩ = E_μ[fg]. E_v averages out v. Write total(v) = E[Var(f | x_{−v})] for an input's total, and V_u for a subset's main variance, as in the docs.

**Claims.**
- (a) Each E_v is idempotent and self-adjoint, and they commute. Commutativity is where independence is used.
- (b) P_u = Π_{v∈u}(I − E_v)·Π_{v∉u}E_v is a complete orthogonal family of projections. f_u = P_u f is the Hoeffding component and V_u = ‖P_u f‖².
- (c) With L = Σ_v(E_v − I), L·P_u = −|u|·P_u.
- (d) Σ_v total(v) = −⟨f, Lf⟩ = Σ_u |u|·V_u, a Dirichlet energy. Efron–Stein is the spectral gap of 1: Var f ≤ Σ_v total(v), with equality exactly when f is additive.
- (e) U_ρ = Π_v(ρI + (1−ρ)E_v) = Σ_u ρ^{|u|}·P_u. For product-form a = Π a_v and b = Π b_v,
  ⟨a, U_ρ b⟩ = Π_v(E a_v·E b_v + ρ·Cov(a_v, b_v)).

**Why (e) matters.** The variance at each interaction order is a coefficient of ρ in Σ_{j,l} a_j a_l ⟨B_j, U_ρ B_l⟩. Its factors are `BassMoments`' cached Eφ_v and Cov(φ_v, φ′_v), so no new integrals are needed. A pair's polynomial costs O(s²) in its s shared variables, against 2^s subsets. At BASS's usual orders that saving is small. It matters at KHAOS's higher orders, and the readout (variance by order) is new either way.

**Transport.**
- *To continuous inputs.* (a)–(d) use only that the E_v are commuting self-adjoint idempotents, so an abstract statement over such a family covers finite and continuous spaces at once; Fubini supplies the family.
- *To other input kinds.* The frame takes a generator per variable. Resampling gives the total index, and Langevin gives DGSM, with C = E[∇f∇fᵀ] as its matrix. A path graph on levels gives an ordinal input's discrete gradient. Their Dirichlet forms are comparable through each generator's own Poincaré constant.
- *To dependent inputs.* Replace E_v by E[· | x_{−v}]. These are still self-adjoint idempotents, but they stop commuting (0.92 in the check), so (b), (c) and the spectrum are lost. The Dirichlet form survives as the generator of random-scan Gibbs dynamics: Var f ≤ (1/γ_μ)·Σ_v total(v), where γ_μ is that chain's gap. So γ_μ measures how far totals undercount under dependence: 1 under independence, and 0 when two inputs are identical, where both totals vanish while the variance stays. The check reads γ_μ = 0.063, Var 0.446 ≤ 4.65.

**Lean.** (c) as the abstract step; (a) and Efron–Stein on the finite product.

```lean
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Real.Basic

open Finset

theorem generator_eigen {E ι : Type*} [AddCommGroup E] [Module ℝ E] [Fintype ι] [DecidableEq ι]
    (P : ι → E →ₗ[ℝ] E) (u : Finset ι) (f : E)
    (h₀ : ∀ i ∈ u, P i f = 0) (h₁ : ∀ i ∉ u, P i f = f) :
    (∑ i, (P i f - f)) = -(u.card : ℝ) • f := by
  sorry

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {X : ι → Type*} [∀ i, Fintype (X i)]
  [∀ i, DecidableEq (X i)]

noncomputable def coordMean (μ : ∀ i, X i → ℝ) (i : ι) (f : (∀ j, X j) → ℝ) :
    (∀ j, X j) → ℝ :=
  fun x => ∑ y : X i, μ i y * f (Function.update x i y)

noncomputable def meanOf (μ : ∀ i, X i → ℝ) (f : (∀ j, X j) → ℝ) : ℝ :=
  ∑ x : ∀ j, X j, (∏ i, μ i (x i)) * f x

theorem coordMean_comm (μ : ∀ i, X i → ℝ) (i j : ι) (f : (∀ j, X j) → ℝ) :
    coordMean μ i (coordMean μ j f) = coordMean μ j (coordMean μ i f) := by
  sorry

theorem efron_stein (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y) (hμ : ∀ i, ∑ y, μ i y = 1)
    (f : (∀ j, X j) → ℝ) :
    meanOf μ (fun x => f x ^ 2) - meanOf μ f ^ 2 ≤
      ∑ i, meanOf μ (fun x => f x * (f x - coordMean μ i f x)) := by
  sorry
```

**Fishable:** no Efron–Stein in Mathlib (`rg -i efron` finds nothing), so this is self-contained. `coordMean_comm` is `Finset.sum_comm` with `Function.update_comm`. **Guards:** i ≠ j in the update swap; μ a probability vector for idempotence.

---

## SG-2 — Sobol and active subspaces have different standing under contraction

**Claim (a).** |√V_u(f) − √V_u(g)| ≤ ‖P_u(f − g)‖ ≤ ‖f − g‖, since P_u is linear and an orthogonal projection. So a contraction rate for f in L²(μ) passes to √V_u at the same rate. A ratio index needs Var f away from zero.

**Claim (b).** C = E[∇f∇fᵀ] has no such bound. A triangle wave of fixed amplitude and n teeth keeps its L² size while its gradient grows with n (the check: L² 0.00577 throughout, gradient 0.02, 0.2, 2.0 at n = 1, 10, 100). It takes about 2n hinges, so it lies in the model class.

**Limits and open pieces.**
- *Which norm.* Contraction theorems here report ‖·‖ₙ at the design points ([PosteriorContraction](../../../docs/Maths/PosteriorContraction.md)), while Sobol reads L²(μ). Passing between them needs a random-design argument and μ comparable to the design. So (a) is conditional on that step.
- *What would give C a rate.* An inverse inequality on a sieve with a lower bound on knot spacing, or interpolation against a smoother truth. BASS's knots sit at data values and can be adjacent, so the spacing bound is not free. Yoo & Ghosal (arXiv 1411.6716) get derivative rates for tensor-product B-splines, and Odin, Bachoc & Lagnoux (arXiv 2403.03540) subspace recovery for GP priors. The searches found nothing for gradient outer products under an adaptive hinge posterior. Both papers known from abstracts only.
- *Kinks under atoms.* Knots sit at data values, and `EmpiricalMeasure` puts atoms there. So under an empirical measure C reads the kink convention (F-X3: the mean of the one-sided limits); under a continuous measure the kinks have mass zero and it does not. E8 should report which it computed.

**Triage.** (a) formalize. (b) as a rate: validate empirically first, by C recovery on functions with known C, at growing N.

---

## SG-3 — a joint mixture of products keeps the pair integrals closed

**Claim.** For μ = Σ_c π_c·⊗_v μ_{c,v} and product-form terms,

E_μ[B_j B_l] = Σ_c π_c·Π_v ∫φ_{jv}φ_{lv} dμ_{c,v},

and the same for the derivative pairs that give C. One component is `ProductMeasure`. N components with delta marginals are a joint empirical sample, whose direct form (evaluate f at the points) is the cheaper path to the same number. So closed-form moments and Monte Carlo over inputs are one formula, with latent-class inputs, product-kernel smoothing and (from memory: Sancetta & Satchell 2004, to confirm) Bernstein copulas, mixtures of Beta products, between them.

**What it keeps and loses.**
- *Kept, exact per draw:* means, pair moments, Var f, C. The law of total variance, Var f = Σ_c π_c Var_c f + Var_c(E_c f), keeps the centred pair shares per component, so `BassMoments`' digits survive.
- *Lost:* Hoeffding's meaning, since inputs are dependent (SG-1). Each component has its own decomposition, and the between-component term reads as a latent label c, which is an extra input of the measure but not of f.
- *Reduced, not closed:* Shapley effects. E_μ[f | x_w] = Σ_c ω_c(x_w)·E_c[f | x_w] with ω_c ∝ π_c·Π_{v∈w} p_{c,v}(x_v) is a ratio, so one outer loop over x_w remains, where the general case needs a double loop.
- *Worse:* the Poincaré constant, as components separate, which loosens active-subspace error bounds.

**In the code.** `MixtureMeasure` mixes one variable's marginal inside `ProductMeasure`, so the joint form is new. `IIntervalMoments` already supplies each component's one-variable integrals. **Literature:** the searches found nothing on this for sensitivity analysis; Iooss & Prieur (arXiv 1707.01334) gives Shapley's analytic forms for Gaussian inputs.

```lean
theorem mixture_product_integral {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    {X : ι → Type*} [∀ i, Fintype (X i)] [∀ i, DecidableEq (X i)]
    (w : κ → ℝ) (μ : κ → ∀ i, X i → ℝ) (a : ∀ i, X i → ℝ) :
    ∑ x : ∀ j, X j, (∑ c, w c * ∏ i, μ c i (x i)) * ∏ i, a i (x i) =
      ∑ c, w c * ∏ i, ∑ y, μ c i y * a i y := by
  sorry
```

**Fishable:** `Finset.prod_univ_sum` for the finite form; `MeasureTheory.integral_fintype_prod_eq_prod` (`Mathlib/MeasureTheory/Integral/Pi.lean`) for the continuous one.

---

## SG-4 — KHAOS's basis is a generator's eigenbasis, which suggests a prior

**Claim (a).** The shifted Legendre polynomials are eigenfunctions of the Jacobi operator: [x(1−x)P′_n]′ = −n(n+1)·P_n.

**Proof route.** With P_n = Σ_k c_k x^k, the operator sends x^k to k²x^{k−1} − k(k+1)x^k. So the claim is (k+1)²c_{k+1} = −(n−k)(n+k+1)c_k, which `coeff_shiftedLegendre` gives directly.

**Claim (b).** Integrating by parts (the boundary term vanishes, since x(1−x) is 0 at both ends) gives orthogonality with no separate proof, and for orthonormal ψ_n, E[x(1−x)ψ′_mψ′_n] = n(n+1)·δ_mn under the uniform measure on [0, 1]. So for f = Σ_α a_α Π_v ψ_{α_v},

Σ_v E[x_v(1−x_v)(∂_v f)²] = Σ_α a_α²·Σ_v α_v(α_v+1),

and for each variable total(v) ≤ ½·E[x_v(1−x_v)(∂_v f)²], with ½ sharp at degree 1. Both come straight from the coefficients.

**The prior this suggests.** KHAOS's g_m = (1 + q(d+q−2))^(−ζ/2) sees only the order q and total degree d. α = (5,1) and (3,3) both give 13; the generator's eigenvalues give 32 and 24. An option beside it:

g_m = (1 + Σ_v α_v(α_v+1))^(−ζ/2).

**Limits.**
- For a main effect of degree d, KHAOS's factor is d^(−ζ/2) and this one about d^(−ζ), so the same ζ means twice the decay. Compare tuned, not at equal ζ.
- The Sobolev reading (Σ_α(1 + λ_α)^s a_α² finite) is plain for a dense expansion. Under sparse reversible-jump selection every draw has finitely many terms, so the prior's smoothness acts through term selection and coefficient size. Any contraction argument has to handle that, for which random-series priors are the nearest theory (Shen & Ghosal, arXiv 1403.0625; abstract only).
- The weight x(1−x) vanishes at the edges, so this functional reads sensitivity less there. The diagonal identity does not make a weighted C diagonal.
- KHAOS reports its ridge best on noise-free simulators and its g-prior on most noisy ones. A new prior has to earn its place predictively, beside those.

R2D2 on polynomial chaos already exists (Bürkner et al., arXiv 2204.06043): that is a prior on variance shares. The smoothness-matched prior is what the searches did not find.

```lean
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.RingTheory.Polynomial.ShiftedLegendre

open Polynomial

theorem shiftedLegendre_jacobi (n : ℕ) :
    derivative ((X * (1 - X)) * derivative (shiftedLegendre n)) =
      -C ((n : ℤ) * (n + 1)) * shiftedLegendre n := by
  sorry
```

**Fishable:** `shiftedLegendre` lives in ℤ[X], with `coeff_shiftedLegendre`, `degree_shiftedLegendre` and a Rodrigues form (`factorial_mul_shiftedLegendre_eq`). No orthogonality or differential equation was seen in that module. **Guards:** Mathlib's sign is P_n(1 − 2x), which the eigenrelation does not see; normalization to ψ_n needs ∫P_n² = 1/(2n+1), for which the Rodrigues form is the route.

---

## SG-5 — E[C]'s top eigenspace is not the mean of the draws' subspaces

**Claim.** Over rank-r projectors P, the top eigenspace of E[C] maximizes tr(P·E[C]) = E_j Σ_k λ_{jk}‖P v_{jk}‖², which weights each draw's directions by its eigenvalues. The chordal mean of the draws' leading subspaces P_j maximizes tr(P·E[P_j]), weighting each draw equally. The two differ: with draws (λ, v) = (1, a), (1, a), (100, b) and a, b at 1.2 radians, they sit 54° apart.

**Reading.** Both are legitimate and mean different things, so offer both, with the Grassmann Fréchet mean a third. E[C] answers "which directions carry the average squared gradient"; the subspace mean answers "which subspace do draws tend to agree on". A draw with a large gradient dominates the first and counts once in the second.

**Fishable:** r = 1 is the Rayleigh quotient, which Mathlib has. Ky Fan for general r was not found (`rg -i 'ky.?fan'` finds nothing). The counterexample needs only explicit 2 × 2 matrices.
