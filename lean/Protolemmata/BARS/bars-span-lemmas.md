# BARS — span lemmas

The continuation of [`bars-multipeak-lemmas.md`](bars-multipeak-lemmas.md) into the FWHM-span readout (the
"MP-5/6" forward-reference there). The span trio is the **level-set analogue** of the peak trio: where
MP-1/MP-2 are about the *critical* set (`f' = 0`) and its argmax, SP-1..3 are about the *level* set
(`f = ℓ`) and its width. They ground the landed `SplineExtrema.SignificantPeakSpans` (exact closed-form
crossings) + `BarsResult.SpanCoverage` π(T).

**Two widths (load-bearing).** STRUCTURAL = the FWHM of the curve itself, a per-draw level-set width;
EPISTEMIC = the credible interval, the spread of the *location* posterior. SP-1..3 are about the
**structural** span — computed *within each draw's curve* — which is exactly why they are level-set lemmas,
not posterior-spread lemmas. Same `DOMAIN-PREMISE` line as the parent note: the crossing arithmetic is
owned here; the drop-fraction/level and the meaning of the span are the consumer's.

---

## SP-1 — span endpoints are the closed-form level-set roots (level-set analogue of MP-1)

The endpoints of a peak's span are the nearest roots of `(f − ℓ)` on each side, and that level set is finite
and exactly enumerable per span (a cubic root per piece) — no scan, the level-set twin of MP-1's critical-set
enumeration. Basin-clamped: the relevant crossing lies between the peak and its bounding valley, and the level
sits strictly above the valley floor (`ℓ = h − dropFraction·prominence`, `dropFraction < 1`), so the crossing
is guaranteed to exist there unless the domain boundary truncates it (the `clip` flag).

```lean
import Mathlib.Algebra.Polynomial.Roots

open Set Polynomial

/-- The ℓ-level set of a non-constant polynomial span is finite (it injects into the roots of `p − C ℓ`).
    The exact, finite candidate set the span endpoints are drawn from — no grid scan. -/
theorem level_set_finite (p : Polynomial ℝ) (ℓ : ℝ) (hp : p ≠ C ℓ) (a b : ℝ) :
    {x ∈ Icc a b | p.eval x = ℓ}.Finite := by
  -- {x | p.eval x = ℓ} = (p - C ℓ).roots (a finite multiset, since p - C ℓ ≠ 0); intersect with Icc.
  sorry
```

**Fishable:** yes — `(p - C ℓ).roots` finiteness is first-class in Mathlib, the same finiteness machinery MP-1
uses for `p.derivative`. **Guards (where they bite):** `p ≠ C ℓ` — a span pinned flat *at* the level has an
infinite crossing set (the genuinely degenerate/clip case the type system should force you to name); and the
cubic-span (`Degree == 3`) hypothesis the implementation's 4-eval recovery already needs.

---

## SP-2 — exact unique crossing per draw ⇒ zero-slop pushforward (argmax-exactness analogue)

Inside a monotone bracket straddling `ℓ` there is a **unique** crossing (IVT existence + strict-monotone
uniqueness), and the implementation returns *that* root in closed form. So each draw's span endpoint carries
zero scan/optimizer slop, and π(T) is an **exact** pushforward of the curve posterior — its width reflects the
genuine structural-plus-epistemic spread and nothing the solver injected. This is the span counterpart of
`argmax_in_closed_form_set`'s zero-slop certificate, and it is the lemma the deterministic `SpanCrossingTests`
already witnesses numerically (edge value = level to 6 digits, via an independent `SplineBasis.Evaluate`).

```lean
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Order.Monotone.Basic

open Set

/-- A monotone bracket straddling ℓ has exactly one crossing — existence by IVT, uniqueness by strict
    monotonicity. The span endpoint IS this root, so the per-draw reduction injects no slop. -/
theorem unique_crossing_in_monotone_bracket
    (f : ℝ → ℝ) (a b ℓ : ℝ) (hab : a ≤ b)
    (hcont : ContinuousOn f (Icc a b)) (hmono : StrictMonoOn f (Icc a b))
    (hlo : f a ≤ ℓ) (hhi : ℓ ≤ f b) :
    ∃! x, x ∈ Icc a b ∧ f x = ℓ := by
  -- existence: intermediate_value_Icc; uniqueness: hmono.injOn on the Icc.
  sorry
```

**Fishable:** yes — `intermediate_value_Icc` + `StrictMonoOn.injOn` are both in Mathlib; no theory gap. This is
the cleanest fish of the trio (and it is already empirically true in the test suite).

---

## SP-3 — FWHM-of-mean ≠ mean-of-FWHM (MP-2 analogue)

The structural width does not survive averaging: the FWHM of the pooled mean curve is **not** the mean of the
per-draw FWHMs — averaging shifts and reshapes the level set. So the span, like the count, **must** be computed
within each draw (which the code does) and pooled, never read off the mean fit `r.Fit`. One-witness existence,
the structural-width counterpart of MP-2's `count_of_mean_undercounts`.

```lean
open Set

/-- Width of the super-level set of `f` at level `ℓ` on `[a,b]` (the FWHM when `ℓ` is half the peak). -/
noncomputable def levelWidth (f : ℝ → ℝ) (a b ℓ : ℝ) : ℝ :=
  (MeasureTheory.volume {x ∈ Icc a b | ℓ ≤ f x}).toReal

/-- FWHM does not commute with averaging: two curves whose mean's width differs from the mean of their widths.
    Certifies span-must-be-per-draw. Witness (deferred): two equal-height humps offset so the mean broadens. -/
theorem fwhm_of_mean_ne_mean_of_fwhm :
    ∃ (f g : ℝ → ℝ) (a b ℓ : ℝ), a < b ∧
      levelWidth (fun x => (f x + g x) / 2) a b ℓ ≠ (levelWidth f a b ℓ + levelWidth g a b ℓ) / 2 := by
  sorry
```

**Fishable:** the witness is explicit (offset humps), but `levelWidth` via `MeasureTheory.volume` makes the
discharge heavier than MP-2's `ncard` witness. Alternative: state the width as `right − left` of the two
crossings (reusing SP-2's unique crossing) to keep it in pure real analysis. Either way one-witness — no theory
gap, just a definitional choice for `levelWidth`.

---

## Fishability ladder & sequencing

| Lemma | Certifies | Fishable now? |
|---|---|---|
| **SP-2** | exact unique per-draw crossing → zero-slop π(T) (already test-witnessed) | ✅ IVT + strict-mono |
| **SP-1** | span endpoints = closed-form finite level set | ✅ polynomial root finiteness |
| **SP-3** | span must be per-draw (width doesn't commute with averaging) | ✅ one witness (def-choice for width) |

**Recommendation:** lead with **SP-2** — it is the smallest true statement, already corroborated numerically by
`SpanCrossingTests`, and a clean IVT + `StrictMonoOn.injOn` close; graduating it to `Enthymemata/BARS.lean`
alongside MP-2 would pair the count and span non-commutation certificates. **SP-1** is the enumeration behind it;
**SP-3** is the per-draw mandate (define `levelWidth` as `right − left` of the SP-2 crossings to keep it pure
real analysis). All three are the structural-width analogues of the peak trio, not posterior-spread claims.
