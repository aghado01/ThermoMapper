import Mathlib

/-!
# BARS — argmax, multi-peak & span lemmas

Formalization targets for the BARS posterior-curve readout. Statements faithfully translate the
protolemmata (work from those, not memory): `fable-BARS-lemma.md` (single peak),
`bars-multipeak-lemmas.md` (MP-1..4), `bars-span-lemmas.md` (SP-1..3).

Graduated here (statements compile, proofs apologize):
* **MP-1** `local_max_is_critical`, `finite_local_max` — the local-max set is the closed-form set.
* **MP-2** `count_of_mean_undercounts` — peak *count* of the mean undercounts.
* **SP-1** `level_set_finite` — the level set is finite (level-set analogue of MP-1).
* **SP-2** `unique_crossing_in_monotone_bracket` — unique exact crossing; the zero-slop endpoint.
* **SP-3** `fwhm_of_mean_ne_mean_of_fwhm` — FWHM doesn't commute with averaging.

Pending (need spline infra / a measure carrier — a placeholder for now): MP-2′ (curve vs posterior
multimodality), MP-3 (prominence-count pushforward), MP-4 (point-process intensity).
-/

open Set

namespace BARS

/-! ### Peaks (MP-1, MP-2) -/

/-- **MP-1.** Fermat: an interior local maximiser of a polynomial span is a root of its
derivative — the closed-form candidate set behind `SignificantPeaks`. -/
theorem local_max_is_critical (p : Polynomial ℝ) (a b t : ℝ) (ht : t ∈ Ioo a b)
    (hmax : IsLocalMax (fun x => p.eval x) t) :
    (Polynomial.derivative p).eval t = 0 := by
  sorry

/-- The local-max set on a non-constant span is finite (it injects into the roots of
`derivative p ≠ 0`) — the finiteness `SignificantPeaks` / `PeakCountMean` rely on. -/
theorem finite_local_max (p : Polynomial ℝ) (hp : Polynomial.derivative p ≠ 0) (a b : ℝ) :
    {t ∈ Icc a b | IsLocalMax (fun x => p.eval x) t}.Finite := by
  sorry

/-- Number of local maxima of `f` on `[a, b]` — finite `ncard` for piecewise-polynomial draws. -/
noncomputable def peakCount (f : ℝ → ℝ) (a b : ℝ) : ℕ :=
  {x ∈ Icc a b | IsLocalMax f x}.ncard

/-- **MP-2.** Peak *count* of the mean undercounts: two curves each with two significant maxima
can have a pointwise mean with only one. The count analogue of `argmax_expectation_noncommute`,
the teeth for reduce-per-draw (peak-detecting the pooled mean curve undercounts the posterior). -/
theorem count_of_mean_undercounts :
    ∃ (f g : ℝ → ℝ) (a b : ℝ), a < b ∧
      peakCount f a b = 2 ∧ peakCount g a b = 2 ∧
      peakCount (fun x => (f x + g x) / 2) a b = 1 := by
  sorry

/-! ### Spans (SP-1, SP-2, SP-3) -/

/-- **SP-1.** The ℓ-level set of a non-constant polynomial span is finite — the closed-form
candidate set the span endpoints are drawn from (the level-set analogue of MP-1). -/
theorem level_set_finite (p : Polynomial ℝ) (ℓ : ℝ) (hp : p ≠ Polynomial.C ℓ) (a b : ℝ) :
    {x ∈ Icc a b | p.eval x = ℓ}.Finite := by
  sorry

/-- **SP-2.** A monotone bracket straddling `ℓ` has exactly one crossing — IVT for existence,
strict monotonicity for uniqueness. The span endpoint is *this* root, so the per-draw reduction
injects zero scan slop and π(T) is an exact pushforward (the span counterpart of
`argmax_in_closed_form_set`); witnessed by `SpanCrossingTests`. -/
theorem unique_crossing_in_monotone_bracket
    (f : ℝ → ℝ) (a b ℓ : ℝ) (hab : a ≤ b)
    (hcont : ContinuousOn f (Icc a b)) (hmono : StrictMonoOn f (Icc a b))
    (hlo : f a ≤ ℓ) (hhi : ℓ ≤ f b) :
    ∃! x, x ∈ Icc a b ∧ f x = ℓ := by
  sorry

/-- Width of the super-level set `{x ∈ [a,b] | ℓ ≤ f x}` (the FWHM when ℓ is half the peak). -/
noncomputable def levelWidth (f : ℝ → ℝ) (a b ℓ : ℝ) : ℝ :=
  (MeasureTheory.volume {x ∈ Icc a b | ℓ ≤ f x}).toReal

/-- **SP-3.** FWHM does not commute with averaging: two curves whose mean's width differs from the
mean of their widths. Certifies the span must be computed per draw, not on the pooled mean fit. -/
theorem fwhm_of_mean_ne_mean_of_fwhm :
    ∃ (f g : ℝ → ℝ) (a b ℓ : ℝ), a < b ∧
      levelWidth (fun x => (f x + g x) / 2) a b ℓ ≠
        (levelWidth f a b ℓ + levelWidth g a b ℓ) / 2 := by
  sorry

end BARS
