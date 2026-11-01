# BASS — the structure chain

Grounds what the mixing results rest on, and the checks a new reversible-jump kernel owes before it lands (F-X4). Gap G1 of the [lineage-gaps register](../../../issues/regression/bass/lineage-gaps-register.md), and the KHAOS pieces of E4 and E5.

**Status (2026-11-01).** Posed. [`lineage-gap-checks.py`](lineage-gap-checks.py) checks SC-1, SC-2 and SC-3 numerically. The Lean is uncompiled: Mathlib's build output is absent from `lean/.lake`, so names were checked against the Mathlib source only.

---

## SC-0 — at fixed τ, the structure chain is a finite reversible chain

**Claim.** With knots at distinct data values and a capped count and order, F is finite for given data. With the coefficients and σ² collapsed and τ held, each structure move is a Metropolis–Hastings kernel on F, reversible with respect to π(S | τ, y). So Aldous & Fill's machinery applies directly, which it does not to BARS's continuous knots. A sampled τ adds one continuous coordinate; a sweep's block scan composes kernels, and a composition need not be reversible even when its parts are.

**Already in the code.** `BassKernelEnumeration` (tests) builds K(S, S′) and π(S | τ, y) at a fixed τ, from the moves' statements, on three fixtures, the continuous one with 2,263 states. It checks detailed balance and πK = π. Nothing computes a gap.

**Triage:** validate empirically. Compute 1 − λ₂ of the symmetrized kernel per move mix, and SC-2's bound over the basins. The fixtures are small and F grows fast with N, so these readings show scaling, not Friedman's barrier itself.

---

## SC-1 — reciprocity, and what KHAOS's printed ratios do

**Claim.** On any set with positive π and positive q both ways, r(s, t) = π(t)q(t, s) / (π(s)q(s, t)) satisfies r(s, t)·r(t, s) = 1, and the kernel q(s, t)·min(1, r(s, t)) satisfies detailed balance. So a printed pair whose product is not 1 is not a Metropolis–Hastings pair for any π and q. Reciprocity is cheap to assert, and belongs in F-X4's contract as a property test: for every enumerated move and its reverse, the product of the two log ratios is 0.

**Reading of arXiv 2510.25036v1 (manuscript only), in the paper's symbols, since this checks its displays.** The acceptance is log α = log-likelihood ratio + log prior ratio + log A_X, so A_X must be reverse over forward. KHAOS's birth draws q₀ with weight q₀^(−s_q)/c_q, coins with product Π_j η_j(q₀)^{χ_j}(1−η_j(q₀))^{1−χ_j} =: H(q₀), a degree with weight d^(−s_d)/c_d, and one of binom(d−1, q−1) partitions. With H̄ = (1/c_q)Σ_{q₀} q₀^(−s_q) H(q₀), the Metropolis–Hastings forms are:

| | Metropolis–Hastings | Printed |
|---|---|---|
| A_Birth(M) | P_D·binom·d^{s_d}c_d / (P_B·(M+1)·H̄) | P_D·binom·H(q₀) / (P_B·(M+1)·q₀^{s_q}c_q·d^{s_d}c_d) |
| A_Death(M) | P_B·M·H̄ / (P_D·binom·d^{s_d}c_d) | P_B·d^{s_d}c_d·H̄ / (P_D·M·binom) |
| A_Mutate1 | d_cand^{s_d}·binom_cand / (d_curr^{s_d}·binom_curr) | d_curr^{s_d}·binom_curr / (d_cand^{s_d}·binom_cand) |

So, as printed:
- the birth's coin product and q₀ weight are inverted, and q₀ is conditioned on where the death marginalizes it;
- the degree weight is inverted in birth, death and Mutate1, and the partition count in Mutate1;
- the death's M is inverted.

The printed pair's product is H(q₀)·H̄ / ((M+1)²·q₀^{s_q}c_q), not 1 (the check matches this form to 16 digits). The (M+1) double count with the prior ratio is already recorded in the port hub's register.

**Also in the manuscript.** P₁(x) = 1 in the recurrence (x is meant). q = Σ_{j=1}^M χ_j (j ≤ p is meant). The partition probability binom(d−1, d−1)⁻¹ (d−1, q−1 is meant). Degrees sampled up to p, where the prior's set stops at d_max. The g₀² conditional has no data term (the extensions brief notes it).

**What this is not.** Not a claim about the `khaos` R package, whose source has not been read, nor about the supplement's Appendix C (the η construction and delayed rejection), which was not read either. η reads the current model's activity, so the reverse birth must read it from the state after the death, which the printed form leaves implicit. Read the code before calling any of this more than a manuscript reading. By the no-parity rule, a KHAOS port derives its ratios and tests reciprocity, whatever the print says.

```lean
import Mathlib.Data.Real.Basic

noncomputable def mhRatio {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S) : ℝ :=
  π t * q t s / (π s * q s t)

theorem mhRatio_mul_swap {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S)
    (hs : π s ≠ 0) (ht : π t ≠ 0) (hst : q s t ≠ 0) (hts : q t s ≠ 0) :
    mhRatio π q s t * mhRatio π q t s = 1 := by
  sorry

theorem mh_detailed_balance {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S)
    (hs : 0 < π s) (ht : 0 < π t) (hst : 0 < q s t) (hts : 0 < q t s) :
    π s * (q s t * min 1 (mhRatio π q s t)) = π t * (q t s * min 1 (mhRatio π q t s)) := by
  sorry
```

---

## SC-2 — a gap is bounded by any set's flow

**Claim.** For a kernel P reversible with respect to π, and any γ with γ·Var f ≤ E(f, f) for all f, every set A has

γ·π(A)(1 − π(A)) ≤ flow(A) = Σ_{s∈A, t∉A} π(s)P(s, t).

**Proof.** Take f = 1_A. Var f = π(A)(1 − π(A)), and reversibility gives E(1_A, 1_A) = ½Σ π(s)P(s, t)(1_A(s) − 1_A(t))² = flow(A).

**Why it matters.** It is the easy half of Cheeger, in the same Dirichlet-form language as [SG-1](sensitivity-generator-lemmas.md). So relaxation time is at least π(A)(1 − π(A)) / flow(A), and any basin gives a lower bound computable on an enumerated kernel. The check has a true gap of 0.295, with flow/variance above it on 20,000 random sets.

```lean
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Real.Basic

open Finset

variable {S : Type*} [Fintype S] [DecidableEq S]

noncomputable def dirichlet (π : S → ℝ) (P : S → S → ℝ) (f : S → ℝ) : ℝ :=
  (1 / 2) * ∑ s, ∑ t, π s * P s t * (f s - f t) ^ 2

noncomputable def varianceOf (π : S → ℝ) (f : S → ℝ) : ℝ :=
  ∑ s, π s * f s ^ 2 - (∑ s, π s * f s) ^ 2

noncomputable def cutFlow (π : S → ℝ) (P : S → S → ℝ) (A : Finset S) : ℝ :=
  ∑ s ∈ A, ∑ t ∈ Aᶜ, π s * P s t

theorem gap_le_cut (π : S → ℝ) (P : S → S → ℝ) (hπ : ∑ s, π s = 1)
    (hrev : ∀ s t, π s * P s t = π t * P t s) (γ : ℝ)
    (hgap : ∀ f : S → ℝ, γ * varianceOf π f ≤ dirichlet π P f) (A : Finset S) :
    γ * ((∑ s ∈ A, π s) * (1 - ∑ s ∈ A, π s)) ≤ cutFlow π P A := by
  sorry
```

**Fishable:** no Cheeger, spectral-gap or conductance result turned up in Mathlib (`rg -il` on those words finds nothing), so this is self-contained. The proof is two finite-sum computations.

---

## SC-3 — KHAOS's posterior precision

**Claim, in the paper's symbols.** With S₀ = g₀²D(g)(ΨᵀΨ)⁻¹D(g) and G_ml = (g₀²g_m g_l + 1)/(g₀²g_m g_l), G∘ΨᵀΨ = ΨᵀΨ + S₀⁻¹ (checked to 1e-13). So the paper's Σ_n = (G∘ΨᵀΨ)⁻¹ is the ordinary conjugate form, and positive definite as a sum of two positive definite matrices.

**Why it matters.** Unless g ∝ 1, this precision is neither `GForm`'s scalar times the Gram nor `RidgeForm`'s Gram plus a diagonal on the Gram's spectrum. It is a third `CoefficientForm`, with its own fast path to find.

```lean
import Mathlib.Data.Real.Basic
import Mathlib.LinearAlgebra.Matrix.Defs

theorem khaos_precision_entry {m : Type*} (A : Matrix m m ℝ) (g : m → ℝ) (g₀ : ℝ) (i j : m)
    (hg₀ : g₀ ≠ 0) (hi : g i ≠ 0) (hj : g j ≠ 0) :
    (g₀ * g i * g j + 1) / (g₀ * g i * g j) * A i j = A i j + A i j / (g₀ * g i * g j) := by
  sorry
```

---

## SC-4 — proof-walks to run

Neither paper has been read past its abstract; both are on disk (`gauntlet/bayesian/2210.09352v1`, `PDFs/Bayesian/Woodard2009.pdf`). The [chip brief](lean-leads-chip-brief.md) takes these walks. Each walk's first step that does not transport is the lemma to pose.

- **Ronen, Saarinen, Tan, Duncan & Yu (arXiv 2210.09352).** A mixing lower bound for a simplified BART (one tree, a subset of moves), exponential in the data size. *Question for BASS.* Leaving a basin here means losing a load-bearing term, which costs about N times the variance per case it explains, over 2σ², and F-M2 measured every term load-bearing and births losing to an Occam cost near ½ ln N. So the barrier should grow linearly in N, which suggests the same exponential form. SC-2 plus SC-0's N sweep tests that before any proof.
- **Woodard, Schmidler & Huber (EJP 2009).** Persistence as a sufficient condition for tempering's torpid mixing, explaining Bhatnagar & Randall's mean-field Potts result (q ≥ 3). *Questions.* Is the steep band a persistence phenomenon, and can a flat-histogram run's ln g measure it? Separately, the [mean-field Potts toy](../../../docs/Maths/Samplers.md) fits tempering replicas at N^3 for q = 3, polynomial where that theorem says torpid, so either the sizes run are too short or the toy's route differs from the theorem's setting. Read that before relying on the match at q = 4.
