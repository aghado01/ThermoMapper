# BGP — transport deltas from TWCD2025

The `maths/regression/bgp/` instrument implements Tang–Wu–Cheng–Dunson 2025 (TWCD2025),
*Adaptive Bayesian Regression on Data with Low Intrinsic Dimensionality*. The paper's own
results — the contraction rate (Thm. in §3.3), the manifold RKHS approximation (Lemma 4.1, Prop. 4.2),
the empirical-prior validity (Prop. 4.4) — are **given**. We do not re-prove them; a step that
rests on one of them is a *citable* `sorry`.

What needs rigorous treatment is only the **delta**: the handful of places where the
implementation extrapolates past, or bends, the paper's stated setup. The counterfactual
proof-walk over the build surfaced three. Posed below for Lean, each with its citation boundary
marked and a note on whether I expect it to land as a finished **lemma** (complete up to citable
sorries) or to sit as an **enthymeme** (a forward boundary I'm not yet satisfied with).

The Lean blocks are sketches — statements with proof-strategy comments, not compiled.

---

### BGP-1: The sampler targets the true bandwidth posterior

`BgpSampler` does not sample `t` directly. It runs random-walk Metropolis on `u = log t` against
the target `g(u) = logMarginal(eᵘ) + logPrior(eᵘ) + u`, and the `+ u` is the change-of-variables
Jacobian (`dt = eᵘ du`). Two implementation facts have to hold for the draws to be the *true*
`t`-posterior `π(t) ∝ L(t)·p(t)`:

1. **(reparametrization)** the pushforward through `t = eᵘ` of the `u`-chain's stationary law is
   exactly `π(t)` — i.e. the `+u` term is right;
2. **(unnormalized invariance)** dropping the prior normalizer `Ẑₙ` (we only ever form ratios)
   does not change the chain.

Both are standard, but they are *our* correctness obligations, not the paper's. (1) is a clean
change-of-variables fact; (2) is detailed-balance invariance under positive rescaling of the
target. **Lemma-tier candidate** — complete up to citing mathlib's MH/COV machinery.

```lean
import Mathlib.MeasureTheory.Measure.Map
import Mathlib.Analysis.SpecialFunctions.Log.Basic

open MeasureTheory

variable (L p : ℝ → ℝ)                      -- marginal likelihood and (unnormalized) EB prior, on t > 0

/-- The t-posterior density (up to the global normalizer), supported on t > 0. -/
noncomputable def piT (t : ℝ) : ℝ := L t * p t

/-- The log-space target the sampler actually uses: π evaluated at eᵘ, times the Jacobian eᵘ. -/
noncomputable def gU (u : ℝ) : ℝ := piT L p (Real.exp u) * Real.exp u

/-- (1) Pushing `gU` forward through `t = eᵘ` recovers `piT`.  Citation boundary: the mathlib
    change-of-variables for `Measure.map` under `exp` (a diffeomorphism (−∞,∞) → (0,∞)). -/
theorem pushforward_logspace_target_eq_posterior :
    -- (Measure.map exp (volume.withDensity gU)) = volume.withDensity piT  (on (0,∞))
    sorry

/-- (2) The Metropolis chain is invariant under `p ↦ c • p` for any `c > 0`; hence dropping `Ẑₙ`
    is harmless.  Citation boundary: Metropolis–Hastings detailed balance. -/
theorem mh_invariant_under_target_rescale (c : ℝ) (hc : 0 < c) :
    sorry
```

*Formalization question:* is it cleaner to carry the acceptance ratio symbolically and prove
ratio-invariance directly, or to lean on a mathlib `Kernel`/detailed-balance API if one exists at
the version we pin? The first keeps it self-contained; the second risks API drift.

---

### BGP-2: Coordinate transforms and the adaptive rate

TWCD2025 assume the data domain `X ⊂ [0,1]^D`. The instrument uses its inputs as given; the BGP
tests and benchmark rescale each coordinate by its training min–max, `φ(x)_k = (x_k − loₖ)/spanₖ`,
before fitting. Corrected 2026-10-19: an earlier version of this note treated the anisotropy of `φ`
as a possible failure of the approximation rate. Anisotropy alone does not break it; the open
question is the data dependence of a fitted transform.

**Fixed affine maps (a deduction, citable).** Let `φ(x) = A x + c` with `A` nonsingular and fixed,
not estimated from the sample. If `M` is a `d`-dimensional smooth closed manifold embedded in
`ℝ^D`, then `φ(M)` is another one. With its induced metric it satisfies Assumption 3, and
`f* ∘ φ⁻¹` is Hölder of the same order `s` on `φ(M)`. Apply the manifold RKHS approximation
(Lemma 4.1 and Prop. 4.2) and the contraction theorem to `φ(M)`. The squared-exponential kernel then
sees the ambient Euclidean distance of the new embedding, which is exactly the theorem's setting.
The powers of the bandwidth and of `n` depend only on `d` and `s`, so they are unchanged. The
constants change: Lemma 4.1's constants depend on the curvature, second fundamental form, volume,
diameter, injectivity radius and reach of the embedded manifold, and `φ` moves all of them.
Per-coordinate scales are a diagonal `A`, so they are covered. This follows from the theorem's
hypotheses; the paper does not analyze preprocessing, and no new proof is claimed here.
**Lemma-tier candidate**, complete up to citable sorries.

**Fitted transforms (the open part).** The min–max map `φₙ` is estimated from the sample, so its
constants `(loₖ, spanₖ)` depend on `n` and on the data, and the transformed design is not a fixed
manifold sampled i.i.d. Carrying the rate over needs at least `φₙ → φ` for a fixed nonsingular `φ`,
and Lemma 4.1's constants uniform over maps near `φ`. A coordinate whose span shrinks with `n` makes
`φₙ` degenerate, and then the constants are not uniform. The empirical-Bayes prior (`T_n`, `v̂_n`)
is computed after the transform, so Prop. 4.4 needs the same conditional argument. **Enthymeme**:
this is a forward boundary, not a citable step.

```lean
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse

variable {D : ℕ}

/-- A fixed affine map of the ambient coordinates. -/
def affineMap (A : Matrix (Fin D) (Fin D) ℝ) (c : Fin D → ℝ) (x : Fin D → ℝ) : Fin D → ℝ :=
  A.mulVec x + c

/-- Citable half.  A fixed nonsingular affine map sends a smooth closed embedded d-manifold to
    another, and intrinsic Hölder-s functions to intrinsic Hölder-s functions, so the contraction
    rate keeps its exponent and only the geometric constants move.  Citation boundary: TWCD2025
    Lemma 4.1, Prop. 4.2 and the manifold contraction theorem, applied to φ(M). -/
theorem fixed_affine_preserves_rate (A : Matrix (Fin D) (Fin D) ℝ) (hA : IsUnit A.det) :
    -- rate(φ(M), f* ∘ φ⁻¹) = rate(M, f*) up to constants depending on A
    sorry := by
  sorry

/-- Open half, NOT citable.  For the sample-fitted min–max map φₙ: if φₙ → φ with φ nonsingular
    and Lemma 4.1's constants are uniform near φ, the rate carries over.  Degenerating spans
    violate the uniformity. -/
theorem fitted_minmax_preserves_rate :
    sorry := by
  sorry
```

*Resolved for fixed maps:* a kernel with fixed per-coordinate scales is the isotropic kernel after a
fixed diagonal affine map, so it falls under the citable half and needs no kernel redesign. What
remains open is the data dependence of `φₙ`.

---

### BGP-3: The natural-log neighbor count still satisfies (A3)

The empirical-Bayes prior uses `k = ⌈γ₂ (ln n)²⌉`. The paper writes `k = ⌈γ₂ log²(n)⌉` with the
log base unpinned, so the natural-log choice is an extrapolation of an ambiguous spec. The claim
to discharge: Prop. 4.4 (the prior satisfies condition (A3) with probability `≥ 1 − 2n^{−10}`) is
unaffected.

I expect this one to **dissolve on a careful read of Prop. 4.4's quantifier** rather than need a
real proof: a base change is a constant factor on `k`, i.e. a rescaling of `γ₂`, and Prop. 4.4 is
stated for *fixed positive* `γ₂` with existentially-quantified constants `c₁, c₂, K₁, …`. If the
quantifier is "for any γ₂ > 0 there exist constants," the base is absorbed and there is no delta —
the proof-walk terminates immediately at the citation. Recording it anyway so the check is on the
ledger, not silently assumed.

```lean
/-- (A3) robustness to the log base.  Citation boundary: Prop. 4.4 with γ₂' = γ₂ · (ln b)²,
    PROVIDED the paper quantifies "for any fixed γ₂ > 0."  If so this is a one-line cite. -/
theorem a3_holds_under_natural_log_k :
    -- prior with k = ⌈γ₂ (ln n)²⌉ satisfies (A3)  ⇐  Prop. 4.4 at γ₂' = γ₂ (ln b)²
    sorry
```

*To confirm against the paper:* read Prop. 4.4's statement of the `γ₂` quantifier (§4.3 / Appendix
B.2). If `γ₂` is fixed-but-arbitrary, mark BGP-3 *not a delta* and drop it.
