# BASS — walking the contraction chain

A counterfactual proof-walk of [PosteriorContraction](../../../docs/Maths/PosteriorContraction.md) (the Bayesian GP's unknown-scale result, UV-1 to UV-6) against BASS's prior. Gap G2 of the [lineage-gaps register](../../../issues/regression/bass/lineage-gaps-register.md). The searches found no posterior contraction result for Bayesian MARS; BART has them.

**Status (2026-11-01).** Walked against Ghosal & van der Vaart (2017), Chapter 8, by reading, and *MARS via LASSO* read. The local cover, an M-term bound and the normalizer's bound are proved in Lean (`BassContractionCover`). No contraction theorem for BASS is proved. The prior-mass step is no longer wholly conditional: for a truth of bounded variation it has a slow rate from the M-term bound; the fast rate stays open. Book references below are to Chapter 8 (`9781139029834.009.pdf`).

**The model.** Fixed design x₁…x_N, y = f*(x) + w with w ~ N(0, σ₀²). A structure S in F, chosen by the count weight and term measure (D01, D02); f = Σ_j a_j B_j on S's columns; a | S, σ², τ from the g-form; σ² ~ IG(c, c·s²) scaled to the response (D03); the level-centred intercept.

---

## The walk

| Step | The GP result's form | BASS | Verdict, by reading |
|---|---|---|---|
| UV-1, Gaussian KL and log-ratio variance | About the likelihood only | Same likelihood | **Transports, with one small lemma** (CW-1). The one-observation identities carry over, and the book's Lemma 2.7 gives the known-σ form. But `slab_subset_klBall` packages them through a sup-norm slab, and BASS's prior-mass step works in ‖·‖ₙ |
| UV-5, tests for ‖f − f′‖ₙ + \|σ − σ′\| | Centred χ² statistics, on the whole parameter set, with σ in a compact range | Same, given compact σ. f unbounded is fine | **Transports**: `BgpGaussianTails`. `exists_tests` quantifies over all of E × [σ_min, σ_max]. The book's own Lemma 8.27(ii) needs less: its exponent is K[‖θ₀ − θ₁‖² + n\|σ₀ − σ₁\|²]/(σ₀²σ₁²), and the remark before it asks only bounded support above, "[0, σ̄]" |
| UV-2, prior mass of the KL ball | TWCD2025's small-ball (59) for a GP | Below: the structure's prior times the coefficient mass | **New.** It needs an M-term rate, which *MARS via LASSO* does not state. `exists_mTerm_approx` gives a slow one for its class (CW-3). The g-form's step holds as written only for g near N (CW-2) |
| UV-3, remaining mass | TWCD2025's (60) | Θₙ = {M ≤ M̄ₙ}: only the count prior's tail | **New, and simple.** The count law is negative binomial, ∝ Γ(h₁ + M)/(M!·(h₂ + 1)^M), tail ratio 1/11 at the defaults. Theorem 8.20(iii′) asks Π(Θₙᶜ)/Π(B) = o(e^{−2nε²}), so M̄ₙ ≈ C·Nε²/ln(h₂ + 1): the log in the rate, as expected. The default cap min(1000, N − 1) has to grow for the asymptotics, so a theorem reads the cap as N − 1 |
| UV-4, entropy | A global sup-norm cover, TWCD2025's (61), with a polynomial radius | Below: a union of finite-dimensional subspaces in ‖·‖ₙ, locally | **New, and this is where it differs. Proved:** `exists_dictionary_local_cover` |

---

## Where it differs: local entropy needs no norm bound

The GP result needs compact amplitude ranges because its entropy is counted in sup norm, where the amplitude multiplies, and a heavy tail forces an exponentially large cutoff. That is the failure PosteriorContraction records for heavy-tailed amplitudes.

**That proof step need not transport here, and Lean now says so.** Under a structure S, f lies in the column space V_S of dimension at most M(S). Theorem 1 of Ghosal & van der Vaart (2007), condition (2.4), asks for a *local* cover: {θ ∈ Θₙ : d_n(θ, θ₀) < r} by balls of radius ξr/4.

`exists_dictionary_local_cover` covers it. For structures of at most M terms from a dictionary D, the part of the union times ℝ within r of (f₀, σ₀), in ‖f − f′‖ + |σ − σ′|, is covered at radius ξr/4 by at most

(|D| + 1)^M · (1 + 32/ξ)^M · (⌈8/ξ⌉ + 1)

balls. The count reads neither r, nor any bound on the coefficients, nor a range for σ. `card_le_of_pairwise_dist` is the volume argument in any finite-dimensional normed space, so it does not use that ‖·‖ₙ is Euclidean. `log_dictionary_count_le` is the walk's form:

log cover ≤ M̄ₙ·log(|D| + 1) + M̄ₙ·log(1 + 32/ξ) + log(⌈8/ξ⌉ + 1) ≲ M̄ₙ·(K·log(pN)),

with |D| at most K_max·C(p, K)·(2N)^K for hinge terms of order up to K. Categorical subsets multiply |D| by their count, which adds their levels to the log.

So Zellner–Siow's heavy tail in g = 1/τ costs nothing in entropy.

**Which route can use this, confirmed against the book.**
- *Theorem 8.19(b) as printed asks only for local entropy.* Its (8.21) is sup_{ε>ε_n} log N(ξε, {θ ∈ Θ_{n,1} : d_n(θ, θ_{n,0}) < 2ε}, e_n) ≤ nε_n², and (b) concludes O_P(e^{−nε_n²}) for sufficiently large M, given (8.20) on B_{n,k} with k > 1 and nε_n² → ∞.
- *The global cover in the fixed-constant axiom is a strengthening,* and the book names it as one: after Lemma 8.21, "weaker, but simpler theorems" replace the local entropy by the global, (8.22) and (8.23) implying (8.20) and (8.21).
- *So a local restatement stays a consequence of the printed theorem.* The hypothesis to put in place of `hentropy` in `fixed_constant_contraction`:

  ```lean
  (hentropy_local : ∀ᶠ n in atTop, ∀ r, εbar n < r → ∃ N : Finset (E × ℝ),
    (N.card : ℝ) ≤ Real.exp (n * εbar n ^ 2) ∧
      {θ ∈ Θ n | d n θ θ₀ < 2 * r} ⊆ ⋃ θ' ∈ N, {θ | d n θ θ' ≤ ξ * r / 2})
  ```

  It sits on the strong side of (8.21) at ε = r, by the same half-radius device the file already uses: a cover by closed balls of radius ξr/2 with arbitrary centres converts to one at radius ξr with centres in the set, where the tests apply. The current global `hentropy` implies it, since ξε̄ₙ/2 ≤ ξr/2 for r > ε̄ₙ. So the restated axiom gives the existing theorem, and widens it. **Proposed only;** `BgpUnknownVarianceTheorem` is unchanged.
- *Theorem 1 of the 2007 paper* takes the local cover as its axiom stands, so the mean-M_n route needs no restatement.

## The prior-mass step

Given an M-term structure S_M with ‖f_M − f*‖ₙ ≤ ε/2, the g-form is what makes the coefficient mass easy. With a ~ N(0, gσ²(BᵀB)⁻¹), Ba is Gaussian on V_S, isotropic in the norm the metric uses, with no design constants. Writing z for its coordinates in an orthonormal basis of V_S, z ~ N(0, gσ²I_M) and ‖f − f_M‖ₙ = ‖z − z_M‖/√N, so

−log Π(‖f − f_M‖ₙ ≤ ε/2 | S_M, σ, g) ≲ (M/2)·log(M g σ²/(Nε²)) + N‖f_M‖ₙ²/(2gσ²).

**CW-2, the role of g.** The second term is linear in N unless g grows with it. At g ≈ N the bound is the first draft's (M/2)·log(Mσ²/ε²) + ‖f_M‖ₙ²/(2σ²). At a constant g it is of order N, which exceeds Nε_n² and breaks the condition. The requirement is g_n ≳ ‖f_M‖ₙ²/ε_n², so **"fixed g" below means a deterministic sequence g_n of order N, not a constant.** SC-0's fixed τ for mixing is a different use of the phrase. Zellner–Siow puts constant mass on g ∈ [N, 2N], which is what the sampled τ needs.

Adding the structure's prior, −log Π(S_M) ≤ M·(K·log(2pN) + log K_max + log(h₂ + 1)) + log W_F, with W_F below, gives the condition M_n·log(pN) ≲ Nε_n².

**So the theorem's shape, given an approximation rate ε(M):** ε_n² ≈ M_n·log(pN)/N where ε(M_n) ≈ ε_n, then a further log from the sieve, so ε̄_n² ≈ M_n·log²(pN)/N. The two rates compose as they do in form 2 of UV-6: Lemma 1 removes the sieve's complement at ε_n, and Theorem 1 runs on the sieve at ε̄_n. p enters only through the structure's prior and the count, and a truth depending on few variables keeps M_n small.

### What *MARS via LASSO* gives (arXiv 2111.11694v5, read)

**The class.** On [0, 1]^d, F_∞^{d,s} holds a₀ + Σ_{|α|≤s} ∫Π_{j∈S(α)}(x_j − t_j)₊ dν_α(t), with each ν_α a finite signed measure, so products of hinges with each variable at most once, as BASS's terms are. Its variation is V_mars = Σ_α |ν_α|([0, 1)^{|α|} ∖ {0}). The atom at t = 0 is excluded, so products of linear functions are unpenalized.

**The results.**
- *Fixed lattice design*, with spacings at least ρ/n_k and sub-Gaussian errors: risk ≤ C_{ρ,d}(σ²V^{1/2}/n)^{4/5}[log(2 + Vn^{1/2}/σ)]^{3(2s−1)/5} + C_{ρ,d}(σ²/n)log²n, so n^{−4/5}(log n)^{3(2s−1)/5}.
- *Random design*, density bounded by B and finite L^{5,1} errors: O_p(n^{−4/5}(log n)^{8(s−1)/5}).
- *Entropy:* log N(ε, D_m, ‖·‖₂) ≤ C_m ε^{−1/2}[log(1/ε)]^{3(2m−1)/4}, through integrated Brownian sheet small-ball bounds.
- *Reduction (Lemma `reduction-to-discrete-measures`):* for any design, f agrees at the design points with a discrete measure on the data-value lattice, with no larger V. That is BASS's knot set, plus 0.
- *Lattice approximation (Lemma `discrete-measures-approximation`):* ‖f_μ − f_ν‖_∞ ≤ (2/N)V and ‖f_μ − f_ν‖_{p₀,2} ≤ C_d(B/N³)^{1/2}V, on a pre-set lattice with up to N^{|α|} atoms per α.

**Two corrections to the first summary.**
- **It states no M-term approximation rate.** The rate is the least squares estimator's, by entropy, not by counting terms. The lattice lemma is a non-adaptive grid, so it gives ε(M) ≈ V·M^{−1/s} in sup norm at best, polynomial in the number of variables, and its C_d depends on d.
- **"Depends only logarithmically on dimension" means less than it reads.** The paper's sense is that n's exponent is free of d and the log's power depends on s alone ("to some extent", in its words). Its constants C_{ρ,d} and C_m depend on d, with no claim that the dependence is logarithmic. So it is a fixed-d statement, not a high-dimensional one. In this walk, p's logarithmic entry comes from the structure's prior and the count, not from this paper.

**CW-3, the M-term step, as far as it goes.**
- *Slow rate, proved:* `exists_mTerm_approx` is Maurey's empirical method in ‖·‖ₙ: a combination with coefficient norm V and atoms of size at most R has an M-term approximation with ‖f − f_M‖ₙ² ≤ V²R²/M, reading neither the dictionary's size nor the design. On [0, 1]^p the hinge products have R ≤ 1, and the reduction lemma puts f* on BASS's knots with the same V. So ε(M)² = V²/M, which gives M_n ≈ V·(N/log(pN))^{1/2}, ε_n² ≈ V·(log(pN)/N)^{1/2}, and ε̄_n² ≈ ε_n²·log(pN), for every K and every p, with constants free of p. This is the theorem's conditional shape made concrete, at a rate slower than least squares'.
- *Fast rate at K = 1, by hand and not in Lean:* a piecewise-linear interpolant has error at most |I|·|ν|(I)/4 on an interval I, and Cauchy–Schwarz on a balanced partition gives ε(M) ≤ V/M² in sup norm with M free knots, per variable. So M_n ≈ (N/log)^{1/5} and ε_n² ≈ (log(pN)/N)^{4/5}, polynomial in the number q of active variables. That recovers the paper's exponent, which is why the additive case looks settled.
- *Fast rate at K ≥ 2: open.* The entropy bound ε^{−1/2} is consistent with ε(M) ≈ M^{−2} up to logs, but entropy does not imply an M-term rate. The paper points to the mixed-smoothness literature (Temlyakov 2018; Dũng, Temlyakov & Ullrich 2018) as where such a result would sit. Not read.
- *Limits on the slow rate.*
  - V excludes the atoms at 0, so products of linear functions need counting separately: include their coefficients in V, or add up to Σ_{k≤K} C(p, k) terms.
  - The class is on [0, 1]^p, and V is not scale-invariant, so inputs are taken as scaled.
  - The approximant has to be in F. The support count costs little: atoms below the threshold each sit on at most min(20, 0.1N) cases, so dropping them costs at most V²·min(20, 0.1N)/N in ‖·‖ₙ², below either rate. D06's angle floor is the part still to say, since the draws can land on near-collinear columns.

---

## Where it will take work

Each of these is a place the walk stops, in the order I would take them.

1. **The response-scaled σ² prior.** IG(c, c·s²) reads s² from y, so the prior is data-dependent, and it is not compact. Route: on the event that s² lies in a fixed interval, whose probability tends to 1, argue uniformly over the hyperparameter. Two refinements from the reading:
   - *The upper tail is the hard one, and one test removes it.* At c = 10⁻³ the tail is t^(−c), so a sieve cut at σ̄ₙ would need σ̄ₙ² exponentially large, and UV-5's K, which divides by σ_max⁴, would collapse. No sieve is needed there. The region {σ ≥ σ̄} for a fixed σ̄ above σ₀ is testable at once: S_{f*}(Y)/n separates it, whatever the alternative's f, at e^{−cN}, as UV-5's scale-up case does. That goes through Theorem 8.22's testing form, or Theorem 8.12 in the i.i.d. section, not through entropy.
   - *The lower tail is exponentially light,* Π(σ² < t) ≈ exp(−c·s²/t), so an ordinary sieve handles it.
   - Yoo & Ghosal treat both a hierarchical IG and an empirical-Bayes σ for tensor B-splines, and Rousseau & Szabó (2017) is the empirical-Bayes reference PosteriorContraction names.
2. **The level-centred intercept** also reads y. PosteriorContraction records the flat constant mean as a separate delta, a candidate BGP-5.
3. **The structure prior's normalization: a constant, and no real obstacle.** `BassTarget` uses D02's weight unnormalized over F, with W_F the sum over F of Γ(h₁ + M)/(Γ(h₁)(h₂ + 1)^M)·Π_t c(t). The intercept-only structure is in F with weight 1, so W_F ≥ 1. The c(t) sum to at most 1, so the unordered M-term collections sum to at most 1/M! (`sum_powersetCard_prod_le`), and W_F ≤ ((h₂ + 1)/h₂)^{h₁}, which is 1.1¹⁰ ≈ 2.59 at the defaults. So the normalizer costs under one nat, free of N and of the design. F-O3's enumerated count law agrees: its 0.6037 at M = 0 is 1/W_F, inside [0.386, 1].
4. **The admission gate (D06) and the support count.** Support refuses near-singular structures, which depends on the design only. The approximating S_M has to be admissible, a condition on the design. The first draft missed the support threshold: a term needs min(20, 0.1·n) nonzero values, so a hinge near an edge is outside F. At 20 of N cases that is a vanishing part of the range, so it should cost only at the boundary, and the approximation step has to say so.
5. **From ‖·‖ₙ to L²(μ).** Sobol's Lipschitz bound ([SG-2](sensitivity-generator-lemmas.md)) reads L²(μ). Random design, the analogue of TWCD2025's Theorem 3.3, is the step PosteriorContraction leaves unformalized too. Knots at data values make the dictionary depend on the design, which fixed design hides and random design does not.

**CW-1, the small lemma.** With ‖f − f*‖ₙ ≤ σ₀ε/6 and |σ − σ₀| ≤ σ₀ε/6, (f, σ) lies in the KL ball of radius ε. UV-1's per-observation bounds K_i ≤ 10D + 2Δ_i and V_i ≤ 18D + 16Δ_i are linear in Δ_i = δ_i²/σ₀², and they need only r ∈ [½, 2], so their sum reads n‖δ‖ₙ²/σ₀². The Lean's `klDiv_le_and_evariance_le_of_slab` bounds each coordinate by the slab, so it wants unpacking one level. Nothing new is needed, and this is posed rather than done because it sits against `BgpGaussianKL`'s internals.

**Start at:** g_n of order N, compact σ, the Theorem-1 route, with the local cover and W_F's bound in hand; then widen to the sampled τ, the point where the walk's claim about local entropy gets used, and to σ's upper tail by the test above.

**Lean.** `BgpGaussianKL` and `BgpGaussianTails` carry over, with CW-1 between. `BassContractionCover` holds the cover, the M-term bound and the normalizer's bound, beside `BgpCovering`, whose interval cover it uses. Theorem 1 stays cited as an axiom there.
