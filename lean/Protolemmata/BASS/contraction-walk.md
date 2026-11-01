# BASS — walking the contraction chain

A counterfactual proof-walk of [PosteriorContraction](../../../docs/Maths/PosteriorContraction.md) (the Bayesian GP's unknown-scale result, UV-1 to UV-6) against BASS's prior. Gap G2 of the [lineage-gaps register](../../../issues/regression/bass/lineage-gaps-register.md). The searches found no posterior contraction result for Bayesian MARS; BART has them.

**Status (2026-11-01).** A walk on paper, with every verdict an expectation. Nothing here is proved, and the approximation paper below is known only from its abstract and is not yet on disk.

**The model.** Fixed design x₁…x_N, y = f*(x) + w with w ~ N(0, σ₀²). A structure S in F, chosen by the count weight and term measure (D01, D02); f = Σ_j a_j B_j on S's columns; a | S, σ², τ from the g-form; σ² ~ IG(c, c·s²) scaled to the response (D03); the level-centred intercept.

---

## The walk

| Step | The GP result's form | BASS | Expected verdict |
|---|---|---|---|
| UV-1, Gaussian KL and log-ratio variance | About the likelihood only | Same likelihood | **Transports** as it stands: `BgpGaussianKL` |
| UV-5, tests for ‖f − f′‖ₙ + \|σ − σ′\| | Centred χ² statistics, on the whole parameter set, with σ in a compact range | Same, given compact σ. f unbounded is fine | **Transports**: `BgpGaussianTails` |
| UV-2, prior mass of the KL ball | TWCD2025's small-ball (59) for a GP | Below: the structure's prior times the coefficient mass | **New.** Needs an approximation rate |
| UV-3, remaining mass | TWCD2025's (60) | Θₙ = {M ≤ M̄ₙ}: only the count prior's tail | **New, and simple**: negative binomial tails are geometric, which costs a log in the rate |
| UV-4, entropy | A global sup-norm cover, TWCD2025's (61), with a polynomial radius | Below: a union of finite-dimensional subspaces in ‖·‖ₙ, locally | **New, and this is where it differs** |

---

## Where it differs: local entropy needs no norm bound

The GP result needs compact amplitude ranges because its entropy is counted in sup norm, where the amplitude multiplies, and a heavy tail forces an exponentially large cutoff. That is the failure PosteriorContraction records for heavy-tailed amplitudes.

**That proof step need not transport here.** Under a structure S, f lies in the column space V_S of dimension M(S), and ‖·‖ₙ is the Euclidean norm on it, scaled. So in the route through Ghosal & van der Vaart (2007), Theorem 1, condition (2.4) asks for a *local* cover: {θ ∈ Θₙ : d_n(θ, θ₀) < r} by balls of radius ξr/4. For one structure, that set lies in a ball of radius r in V_S (times an interval for σ), and is covered by a constant-to-the-M number of balls, with no dependence on r or on any norm bound on the coefficients. Over the union,

log cover ≲ log|{S : M(S) ≤ M̄ₙ}| + M̄ₙ·log(c/ξ) ≲ M̄ₙ·(K·log(pN)).

So Zellner–Siow's heavy tail in g = 1/τ may cost nothing in entropy, and it puts constant mass on g near N for the prior-mass step. **Which route can use this:** Theorem 1 as its axiom stands. The fixed-constant axiom in `BgpUnknownVarianceTheorem` asks for a global cover, a strengthening made for the GP case, which needs a norm bound. The book's Theorem 8.19(b) itself asks for local entropy, so restating that axiom on the printed side would serve this route too.

## The prior-mass step

Given an M-term structure S_M with ‖f_M − f*‖ₙ ≤ ε/2, the g-form is what makes the coefficient mass easy. With a ~ N(0, gσ²(BᵀB)⁻¹), Ba is Gaussian on V_S in the same norm the metric uses, so with g near N,

−log Π(‖f − f_M‖ₙ ≤ ε/2 | S_M, σ, g) ≲ (M/2)·log(Mσ²/ε²) + ‖f_M‖ₙ²/(2σ²),

with no design constants. Adding the structure's prior, −log Π(S_M) ≲ M·K·log(pN), gives the condition M_n·log(pN) ≲ Nε_n².

**So the theorem's shape, conditional on an approximation rate ε(M):** ε_n² ≈ M_n·log(pN)/N where ε(M_n) ≈ ε_n, then a further log from the sieve, so ε̄_n² ≈ M_n·log²(pN)/N. The two rates compose as they do in form 2 of UV-6: Lemma 1 removes the sieve's complement at ε_n, and Theorem 1 runs on the sieve at ε̄_n. p enters only logarithmically, and a truth depending on few variables keeps M_n small.

**The approximation class.** Ki, Fang & Guntuboyina, *MARS via LASSO* (*Annals of Statistics* 2024, arXiv 2111.11694) study least squares over infinite combinations of MARS basis functions under a variation constraint, with a rate depending on dimension only logarithmically. That is the natural class for ε(M), and an M-term approximation result is what to look for in it. Abstract only so far.

---

## Where it will take work

Each of these is a place the walk stops, in the order I would take them.

1. **The response-scaled σ² prior.** IG(c, c·s²) reads s² from y, so the prior is data-dependent, and it is not compact. Route: on the event that s² lies in a fixed interval, whose probability tends to 1, argue uniformly over the hyperparameter. PosteriorContraction notes that inverse-gamma tails on the noise alone can be handled on a log scale; its tests assume a compact σ, so start there. Rousseau & Szabó (2017) is the empirical-Bayes reference it names.
2. **The level-centred intercept** also reads y. PosteriorContraction records the flat constant mean as a separate delta, a candidate BGP-5.
3. **The structure prior's normalization.** `BassTarget` uses D02's count weight unnormalized over F. The theorem needs a proper Π, so the lower bound has to carry the normalizer's dependence on N.
4. **The admission gate (D06).** Support refuses near-singular structures, which depends on the design only. The approximating S_M has to be admissible, a condition on the design.
5. **From ‖·‖ₙ to L²(μ).** Sobol's Lipschitz bound ([SG-2](sensitivity-generator-lemmas.md)) reads L²(μ). Random design, the analogue of TWCD2025's Theorem 3.3, is the step PosteriorContraction leaves unformalized too. Knots at data values make the dictionary depend on the design, which fixed design hides and random design does not.

**Start at:** fixed g, compact σ, a proper normalized structure prior, the Theorem-1 route; then widen to the sampled τ, the point where this walk's claim about local entropy gets tested.

**Lean.** `BgpGaussianKL` and `BgpGaussianTails` carry over. The new piece is a cover of a union of Euclidean balls in subspaces, beside `BgpCovering`'s lemmas, with Theorem 1 cited as an axiom as there.
