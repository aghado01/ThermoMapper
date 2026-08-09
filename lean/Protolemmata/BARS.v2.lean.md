import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.Analysis.Calculus.LocalExtr.Basic
import Mathlib.Data.Set.Card
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import Mathlib.Tactic.NormNum
import Mathlib.Topology.Order.IntermediateValue

/-
  PROTOLEMMA STATUS (2026-08-08)

  Preserved here rather than under Enthymemata because a standalone Lean 4.32.2
  check does not yet elaborate. Two calls use the obsolete
  `Polynomial.setOf_isRoot_finite` name, and the measure rewrites in
  `peak_width_eq_bracket_roots` and `levelWidth_not_affine` leave goals open.
  The statements and remaining obligations also await the Fable 5 semantic
  review. This is therefore revision source, not an active formal module.
-/

/-!
# BARS — argmax, multi-peak & span lemmas (revised)

Formalization targets for the BARS posterior-curve readout. Statements faithfully translate the
protolemmata (work from those, not memory): `fable-BARS-lemma.md` (single peak),
`bars-multipeak-lemmas.md` (MP-1..4), `bars-span-lemmas.md` (SP-1..3).

## Revisions in this pass

* **MP-1 is interval-free.** `IsLocalMax` is the two-sided filter statement
  `∀ᶠ x in 𝓝 t, f x ≤ f t`; Fermat needs no bracket. The old `Ioo` hypothesis was dead weight
  and blocked applying MP-1 at the `Icc` endpoints inside `finite_local_max`.
* **SP-2 gains its antitone dual.** A peak has two shoulders; FWHM needs the up-crossing and
  the down-crossing.
* **Per-peak width is the bracket-root difference.** `superlevel_of_peak_bracket` shows the
  super-level set of a single-peak bracket is exactly `Icc xL xR`, so `levelWidth` restricted
  to the bracket computes `xR − xL`. Globally, `levelWidth` is *occupancy* — total measure
  above the level, summed across every hump — and is a peak width only on a bracket.
* **SP-3 splits.** `levelWidth_not_affine` keeps the free-level statement under an honest name
  (constants witness it; it certifies only non-affinity). The FWHM claim proper is
  `fwhm_of_mean_ne_mean_of_fwhm`, with the level *draw-relative* (each curve's own half-max) —
  that is the statement that indicts pooled-mean reduction.
* **MP-2′ graduates in witness form.** Two single-peak draws whose pointwise mean is bimodal —
  the overcount direction. Jointly, MP-2 (undercount instance) and MP-2′ (overcount instance)
  show the peak-count functional admits **no signed bound under mixing**: it is not
  conservatively biased in a correctable direction, it is directionless. Reduce-per-draw is
  justified by the *absence* of any one-sided theorem, not by a chosen direction. (The
  measure-carrier version of MP-2′ — posterior multimodality proper — stays pending.)

## Contracts

* `peakCount` uses full-topology `IsLocalMax` (weak `≤`), so: (i) boundary maxima of the
  *restricted* curve are not counted — interior features only; (ii) any plateau makes the set
  infinite and `ncard` junks to `0`. Safe on non-constant polynomial spans; document at the
  call sites feeding `SignificantPeaks` / `PeakCountMean`.
* Everything here is per-span. Globalization keystone (pending spline infra): a spline draw's
  candidate set is the finite union of per-span root sets, exhaustive at knots *because* `C¹`
  gluing forces knot extrema to be critical points of both adjacent pieces. Below `C¹`, knots
  join the candidate set explicitly, like range endpoints.

Proof status: mechanical proofs are attempted against Mathlib from memory (drafted, not
compiled — expect name-level fixes, flagged inline); structural `sorry`s carry their intended
route in comments.
-/

open Set

namespace BARS

/-! ### Peaks (MP-1, MP-2, MP-2′) -/

/-- **MP-1** (interval-free). Fermat: a local maximiser of a polynomial is a root of its
derivative — the closed-form candidate set behind `SignificantPeaks`. `IsLocalMax` is
two-sided, so no bracket hypothesis belongs here; this is what lets `finite_local_max` use it
at `Icc` endpoints. -/
theorem local_max_is_critical (p : Polynomial ℝ) (t : ℝ)
    (hmax : IsLocalMax (fun x => p.eval x) t) :
    (Polynomial.derivative p).eval t = 0 := by
  -- `IsLocalMax.deriv_eq_zero` is unconditional (the junk value `deriv = 0` at
  -- non-differentiable points absorbs the case split); `Polynomial.deriv` rewrites
  -- `deriv (p.eval ·)` to the formal derivative.
  simpa [Polynomial.deriv] using hmax.deriv_eq_zero

/-- The local-max set on a non-constant span is finite: by MP-1 it injects into the root set of
`derivative p ≠ 0` — the finiteness `SignificantPeaks` / `PeakCountMean` rely on. -/
theorem finite_local_max (p : Polynomial ℝ) (hp : Polynomial.derivative p ≠ 0) (a b : ℝ) :
    {t ∈ Icc a b | IsLocalMax (fun x => p.eval x) t}.Finite := by
  have hroots : {x : ℝ | (Polynomial.derivative p).IsRoot x}.Finite :=
    Polynomial.setOf_isRoot_finite hp  -- NOTE: check current name (`finite_setOf_root`?)
  refine hroots.subset ?_
  rintro t ⟨-, hmax⟩
  exact local_max_is_critical p t hmax

/-- Number of local maxima of `f` on `[a, b]`. Contract: meaningful only when the peak set is
finite (e.g. non-constant polynomial spans); plateaus make the set infinite and `ncard` junks
to `0`; boundary maxima of the restricted curve are excluded by the full-topology predicate. -/
noncomputable def peakCount (f : ℝ → ℝ) (a b : ℝ) : ℕ :=
  {x ∈ Icc a b | IsLocalMax f x}.ncard

/-- Witness span for MP-2: the antiderivative of `−(X−1)(X−2)(X−3)`. The derivative's sign
pattern `+,−,+,−` across roots `1, 2, 3` gives local maxima exactly `{1, 3}` (and a minimum
at `2`). -/
noncomputable def mp2F : Polynomial ℝ :=
  Polynomial.C (-(1 : ℝ)/4) * Polynomial.X ^ 4 + Polynomial.C 2 * Polynomial.X ^ 3
    + Polynomial.C (-(11 : ℝ)/2) * Polynomial.X ^ 2 + Polynomial.C 6 * Polynomial.X

/-- Reflection witness: `mp2G.eval x = mp2F.eval (−x)` (via `Polynomial.eval_comp`);
maxima exactly `{−3, −1}`. -/
noncomputable def mp2G : Polynomial ℝ := mp2F.comp (-Polynomial.X)

/-- **MP-2.** Peak *count* of the mean can undercount: two curves each with two maxima whose
pointwise mean has one. Count analogue of `argmax_expectation_noncommute`; teeth for
reduce-per-draw. NB the *direction* is not general — see MP-2′.

Witness: `f := mp2F.eval`, `g := mp2G.eval` on `[−4, 4]`. The mean's derivative is the odd part
of `f′`, which collapses to `−x(x² + 11)` — one sign change, hence exactly one maximum, at
`0`. -/
theorem count_of_mean_undercounts :
    ∃ (f g : ℝ → ℝ) (a b : ℝ), a < b ∧
      peakCount f a b = 2 ∧ peakCount g a b = 2 ∧
      peakCount (fun x => (f x + g x) / 2) a b = 1 := by
  -- Route: `use ⟨mp2F.eval, mp2G.eval, −4, 4⟩`. Factor the derivatives; `nlinarith` the factor
  -- signs on each interval; `strictMonoOn_of_deriv_pos` / `strictAntiOn_of_deriv_neg` on the
  -- pieces; assemble `IsLocalMax` at each junction from mono-left / anti-right on an `Ioo`
  -- neighbourhood (`Filter.eventually_of_mem (Ioo_mem_nhds ..)` with a case split).
  -- Exactness (no other maxima): interval-free MP-1 excludes non-critical points; `2` (resp.
  -- `−2`) is a minimum of `f` (resp. `g`); the mean's only critical point is `0`. Mechanical,
  -- ~150 lines.
  sorry

/-- **MP-2′** (witness form). The overcount direction: two curves each with a *single* maximum
whose pointwise mean has two — the mean manufactures structure no draw has. With MP-2 this
shows peak count admits no signed bound under mixing. The measure-carrier version (pointwise
mean of a genuinely multimodal posterior) stays pending.

Witness: `f x = −|x + 1| + x²/4`, `g x = −|x − 1| + x²/4` on `[−3, 3]`. Piecewise-quadratic
with no plateaus (the `x²/4` tilt is what dodges the `ncard` junk that kills bare tents):
`f` peaks only at the corner `−1`, `g` at `1`; the mean is `−1 + x²/4` on `[−1, 1]`, rising
into both corners and falling beyond them, so its maxima are exactly `{−1, 1}`. -/
theorem count_of_mean_overcounts :
    ∃ (f g : ℝ → ℝ) (a b : ℝ), a < b ∧
      peakCount f a b = 1 ∧ peakCount g a b = 1 ∧
      peakCount (fun x => (f x + g x) / 2) a b = 2 := by
  -- Route: corner `IsLocalMax` from strict mono/anti on the adjacent quadratic pieces;
  -- endpoint exclusion is automatic (full topology — the curve keeps moving outside the
  -- window); interior exactness from explicit derivative signs (`±1 + x/2`, `−1 + x/2`,
  -- and `x/2` on the middle piece of the mean).
  sorry

/-! ### Spans (SP-1, SP-2 + dual, widths, SP-3) -/

/-- **SP-1.** The ℓ-level set of a non-constant polynomial span is finite — the closed-form
candidate set the span endpoints are drawn from (the level-set analogue of MP-1). -/
theorem level_set_finite (p : Polynomial ℝ) (ℓ : ℝ) (hp : p ≠ Polynomial.C ℓ) (a b : ℝ) :
    {x ∈ Icc a b | p.eval x = ℓ}.Finite := by
  have h0 : p - Polynomial.C ℓ ≠ 0 := sub_ne_zero.mpr hp
  have hroots : {x : ℝ | (p - Polynomial.C ℓ).IsRoot x}.Finite :=
    Polynomial.setOf_isRoot_finite h0  -- NOTE: same name check as above
  refine hroots.subset ?_
  rintro x ⟨-, hx⟩
  simp [Polynomial.IsRoot, Polynomial.eval_sub, Polynomial.eval_C, hx]

/-- **SP-2** (rising shoulder). A monotone bracket straddling `ℓ` has exactly one crossing —
IVT for existence, strict monotonicity for uniqueness. The span endpoint is *this* root, so
the per-draw reduction injects zero scan slop and `π(T)` is an exact pushforward; witnessed by
`SpanCrossingTests`. -/
theorem unique_crossing_in_monotone_bracket
    (f : ℝ → ℝ) (a b ℓ : ℝ) (hab : a ≤ b)
    (hcont : ContinuousOn f (Icc a b)) (hmono : StrictMonoOn f (Icc a b))
    (hlo : f a ≤ ℓ) (hhi : ℓ ≤ f b) :
    ∃! x, x ∈ Icc a b ∧ f x = ℓ := by
  obtain ⟨x, hx, hfx⟩ := intermediate_value_Icc hab hcont ⟨hlo, hhi⟩
  refine ⟨x, ⟨hx, hfx⟩, ?_⟩
  rintro y ⟨hy, hfy⟩
  exact hmono.injOn hy hx (hfy.trans hfx.symm)

/-- **SP-2, antitone dual** (falling shoulder). Same content for the down-crossing; a peak's
width needs both shoulders. -/
theorem unique_crossing_in_antitone_bracket
    (f : ℝ → ℝ) (a b ℓ : ℝ) (hab : a ≤ b)
    (hcont : ContinuousOn f (Icc a b)) (hanti : StrictAntiOn f (Icc a b))
    (hhi : ℓ ≤ f a) (hlo : f b ≤ ℓ) :
    ∃! x, x ∈ Icc a b ∧ f x = ℓ := by
  obtain ⟨x, hx, hfx⟩ := intermediate_value_Icc' hab hcont ⟨hlo, hhi⟩
  refine ⟨x, ⟨hx, hfx⟩, ?_⟩
  rintro y ⟨hy, hfy⟩
  exact hanti.injOn hy hx (hfy.trans hfx.symm)

/-- Measure of the super-level set `{x ∈ [a,b] | ℓ ≤ f x}` — *occupancy*: total length above
the level, summed across every hump. It equals a peak's width only on a single-peak bracket
(`peak_width_eq_bracket_roots`); never read it as an FWHM on a multi-peak window. -/
noncomputable def levelWidth (f : ℝ → ℝ) (a b ℓ : ℝ) : ℝ :=
  (MeasureTheory.volume {x ∈ Icc a b | ℓ ≤ f x}).toReal

/-- On a single-peak bracket — rising on `[a, t]`, falling on `[t, b]` — the super-level set is
exactly the interval between the shoulder crossings. Order-theoretic; no continuity needed. -/
theorem superlevel_of_peak_bracket
    (f : ℝ → ℝ) (a t b ℓ xL xR : ℝ)
    (hmono : StrictMonoOn f (Icc a t)) (hanti : StrictAntiOn f (Icc t b))
    (hxL : xL ∈ Icc a t) (hfxL : f xL = ℓ)
    (hxR : xR ∈ Icc t b) (hfxR : f xR = ℓ) :
    {x ∈ Icc a b | ℓ ≤ f x} = Icc xL xR := by
  -- ⊆: split on `x ≤ t`. A point of `[a, xL)` sits strictly below `ℓ` by `hmono` against
  --    `xL`; `(xR, b]` symmetrically by `hanti` against `xR`.
  -- ⊇: `[xL, t]` is `≥ ℓ` by `hmono` from `xL` (equality at `xL`); `[t, xR]` by `hanti`
  --    into `xR`. Membership in `Icc a b` from `hxL.1` and `hxR.2`.
  sorry

/-- Per-peak width: on a single-peak bracket, `levelWidth` computes the difference of the two
shoulder roots delivered by the bracket lemmas — the FWHM when `ℓ` is half the peak. This is
the reconciliation of the occupancy functional with the per-draw span readout. -/
theorem peak_width_eq_bracket_roots
    (f : ℝ → ℝ) (a t b ℓ xL xR : ℝ)
    (hmono : StrictMonoOn f (Icc a t)) (hanti : StrictAntiOn f (Icc t b))
    (hxL : xL ∈ Icc a t) (hfxL : f xL = ℓ)
    (hxR : xR ∈ Icc t b) (hfxR : f xR = ℓ) :
    levelWidth f a b ℓ = xR - xL := by
  have hset := superlevel_of_peak_bracket f a t b ℓ xL xR hmono hanti hxL hfxL hxR hfxR
  have hle : xL ≤ xR := le_trans hxL.2 hxR.1
  simp [levelWidth, hset, Real.volume_Icc, ENNReal.toReal_ofReal (sub_nonneg.mpr hle)]

/-- **SP-3a** (renamed from the old SP-3). `levelWidth` is not affine in the curve: the mean's
occupancy differs from the mean of occupancies. Constants witness it — this certifies
non-affinity of the functional and *nothing about FWHM*. Kept as the zero-cost sanity lemma;
the pooled-mean indictment is SP-3 below. -/
theorem levelWidth_not_affine :
    ∃ (f g : ℝ → ℝ) (a b ℓ : ℝ), a < b ∧
      levelWidth (fun x => (f x + g x) / 2) a b ℓ ≠
        (levelWidth f a b ℓ + levelWidth g a b ℓ) / 2 := by
  refine ⟨fun _ => 1, fun _ => -1, 0, 1, 0, one_pos, ?_⟩
  have h1 : {x ∈ Icc (0 : ℝ) 1 | (0 : ℝ) ≤ 1} = Icc (0 : ℝ) 1 := by ext x; norm_num
  have h2 : {x ∈ Icc (0 : ℝ) 1 | (0 : ℝ) ≤ -1} = (∅ : Set ℝ) := by ext x; norm_num
  have h3 : {x ∈ Icc (0 : ℝ) 1 | (0 : ℝ) ≤ (1 + -1) / 2} = Icc (0 : ℝ) 1 := by
    ext x; norm_num
  -- widths: mean ↦ 1, f ↦ 1, g ↦ 0; and 1 ≠ (1 + 0) / 2.
  simp [levelWidth, h1, h2, h3, Real.volume_Icc]
  norm_num

/-- Draw-relative level: half of the curve's supremum on the window. FWHM's defining feature —
the level depends on the draw, which is what makes the functional *doubly* nonlinear. -/
noncomputable def halfMax (f : ℝ → ℝ) (a b : ℝ) : ℝ := sSup (f '' Icc a b) / 2

/-- FWHM as occupancy at the draw's own half-max. On a single-peak window this equals the
shoulder-root difference via `peak_width_eq_bracket_roots`. -/
noncomputable def fwhm (f : ℝ → ℝ) (a b : ℝ) : ℝ := levelWidth f a b (halfMax f a b)

/-- **SP-3.** FWHM does not commute with averaging *even for genuine single-peak draws with the
level tied to each curve*: the mean's FWHM (at the mean's own half-max) differs from the mean
of the draws' FWHMs. Certifies the span must be computed per draw, not on the pooled mean fit.

Witness on `[−2, 2]`: tent `f x = max 0 (1 − |x|)` (half-max `1/2`, width `1`; its plateau is
fine here — `fwhm` is measure-based, not `peakCount`-based) and parabola `g x = 2 − x²/2`
(half-max `1`, width `2√2`). The mean has supremum `3/2` at `0`, half-max `3/4`, and
super-level set exactly `Icc (−1) 1` (inner piece: `x² + 2|x| − 3 = (|x|+3)(|x|−1) ≤ 0`),
so width `2`; and `2 ≠ (1 + 2√2)/2` since `2√2 < 3`. -/
theorem fwhm_of_mean_ne_mean_of_fwhm :
    ∃ (f g : ℝ → ℝ) (a b : ℝ), a < b ∧
      ContinuousOn f (Icc a b) ∧ ContinuousOn g (Icc a b) ∧
      fwhm (fun x => (f x + g x) / 2) a b ≠ (fwhm f a b + fwhm g a b) / 2 := by
  -- Route: the three suprema are `1`, `2`, `3/2`, attained at `0` (`IsGreatest.csSup_eq`);
  -- the three super-level sets are `Icc (−1/2) (1/2)`, `Icc (−√2) √2`, `Icc (−1) 1`;
  -- widths via `Real.volume_Icc`; conclude with `Real.sq_sqrt` / `nlinarith` from `8 < 9`.
  sorry

/-! ### Pending (spline infra / measure carrier)

* Globalization keystone: a spline draw's candidate set is the finite union of per-span root
  sets; exhaustive at knots iff the gluing is `C¹` (knot extrema are then critical points of
  both adjacent pieces). This is the bridge from the per-span lemmas above to what
  `SignificantPeaks` actually scans.
* MP-2′, measure form: pointwise-mean curve of a genuinely multimodal posterior vs draws.
* MP-3: prominence-count pushforward.  MP-4: point-process intensity of peak locations.

Notes:

Drop-in revision, everything promised plus MP-2′ graduated early. Status split: five proofs are attempted in full (MP-1 as a two-liner, `finite_local_max`, SP-1, both SP-2 shoulders, SP-3a, and the width corollary), four sorries remain with their routes in comments (MP-2 with the `mp2F`/reflection witnesses now defined in the file, MP-2′, the bracket set-equality, and draw-relative SP-3 with the tent-vs-parabola witness worked out in the docstring).

One design note worth knowing: the natural MP-2′ witness — two tents — dies on the `ncard` junk value, because a tent's flat tails are all weak local maxima and the peak set goes infinite. The witness in the file adds an `x²/4` tilt so every piece is strictly monotone: `f x = −|x+1| + x²/4`, `g` its mirror, mean `−1 + x²/4` on the middle piece with corners at ±1. Counts 1/1/2, no plateaus, and endpoint exclusion comes free from the full-topology predicate.

The Lean 4.32.2 check confirms the two flagged name/proof gaps: the root-finiteness lemma
(`Polynomial.setOf_isRoot_finite` vs. its current replacement) and the simp choreography in the
two measure proofs (`Real.volume_Icc` / `toReal_ofReal`). The logic of every attempted proof is
sound even though those declarations do not yet elaborate.

-/


end BARS
