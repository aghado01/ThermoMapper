# BGP — transport deltas from TWCD2025

The `maths/regression/bgp/` instrument implements Tang–Wu–Cheng–Dunson 2025 (TWCD2025),
*Adaptive Bayesian Regression on Data with Low Intrinsic Dimensionality*. The paper's own
results — the contraction rate (Thm. in §3.3), the manifold RKHS approximation (Lemma 4.1),
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

### BGP-2: Normalization preserves the adaptive rate (the real interface)

TWCD2025 *assumes* the data domain `X ⊂ [0,1]^D`. The implementation *imposes* it: each ambient
coordinate is min–max rescaled, `φ(x)_k = (x_k − loₖ)/spanₖ`. This is a **per-coordinate** (hence
**anisotropic**) affine map, and that anisotropy is the thing I'm not yet satisfied transports.

The benign part is citable: `φ` with `0 < spanₖ < ∞` is bi-Lipschitz, and the covering-number /
Minkowski-dimension condition (A1) is a bi-Lipschitz invariant — so the intrinsic dimension `ϱ`
survives, constants aside. The part that opens a forward boundary is the **RKHS approximation**
(A2 / Lemma 4.1): the squared-exponential kernel is isotropic in `‖x−x'‖²`, but after `φ` it sees
the *weighted* distance `Σₖ (xₖ−x'ₖ)²/spanₖ²`. Does `F^ε`'s approximation of `f*` hold at the same
`s`-rate under that warped metric, or does the anisotropy interact with the manifold's second
fundamental form in a way the isotropic analysis didn't cover? That is the unsatisfied step.
**Enthymeme** until the proof-walk through Lemma 4.1's approximation closes — the `sorry` there is
*not* citable; it's the gap I opened.

```lean
import Mathlib.Topology.MetricSpace.Lipschitz
import Mathlib.Topology.MetricSpace.HausdorffDimension

variable {D : ℕ}
variable (lo span : Fin D → ℝ) (hspan : ∀ k, 0 < span k)

/-- Per-coordinate min–max normalization: anisotropic affine. -/
def normalize (x : Fin D → ℝ) : Fin D → ℝ := fun k => (x k - lo k) / span k

/-- Benign half — citable.  `normalize` is bi-Lipschitz, so it preserves the covering-number
    dimension of the data domain.  Citation boundary: bi-Lipschitz invariance of Minkowski
    dimension (standard; Falconer / mathlib `dimH`-style results). -/
theorem normalize_preserves_intrinsic_dimension :
    AntilipschitzWith sorry (normalize lo span) ∧ LipschitzWith sorry (normalize lo span) := by
  sorry

/-- Open half — the forward boundary, NOT citable.  Under the anisotropic warp the SE kernel's
    RKHS approximation (A2 / Lemma 4.1) of an intrinsic Hölder `f*` must still achieve the
    `ε^{s/2}` sup-error.  This is what I'm exploring: does the isotropic approximation argument
    transport, or does coordinate-anisotropy × curvature break the rate? -/
theorem normalized_rkhs_approximation_holds :
    -- ∀ ε < ε₀, ∃ F ∈ Hε(φ X), ‖F − f*‖∞ ≤ ν₁ ε^{s/2} ∧ ‖F‖²_{Hε} ≤ ν₂ ε^{−ℓ/2}
    sorry
```

*Open question:* the cleanest escape may be to make the kernel anisotropic-by-construction (absorb
`spanₖ` into per-coordinate bandwidths) and show that is exactly an isotropic fit in the original
coordinates — collapsing the delta. Worth checking before grinding the warped-metric approximation.

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
