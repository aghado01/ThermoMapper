import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic.NormNum

/-!
# Cumulative fields — exponential cut coordinate

The reusable analytic kernel for cumulative-field clustering. A nonnegative
edge score is mapped through the exponential survival link, then thresholded.
The two apparent controls, temperature and affinity threshold, factor through
one effective cut coordinate.
-/

open Real

namespace CumulativeField

/-- Exponential affinity associated with score `H` at positive scale `T`. -/
noncomputable def affinity (H T : ℝ) : ℝ := 1 - exp (-(H / T))

/-- The score multiplier induced by an affinity threshold `theta`. -/
noncomputable def cutScale (theta : ℝ) : ℝ := -log (1 - theta)

/-- Thresholding exponential affinity is exactly thresholding the underlying
score at `T * cutScale theta`. The algebraic equivalence needs only `theta < 1`;
the usual probability-threshold interpretation additionally assumes
`0 ≤ theta`. -/
theorem affinity_gt_iff (H T theta : ℝ) (hT : 0 < T) (htheta : theta < 1) :
    affinity H T > theta ↔ H > T * cutScale theta := by
  have hOneSub : 0 < 1 - theta := sub_pos.mpr htheta
  unfold affinity cutScale
  rw [gt_iff_lt, gt_iff_lt, lt_sub_comm,
    ← lt_log_iff_exp_lt hOneSub]
  constructor
  · intro h
    have hdiv : -log (1 - theta) < H / T := by
      have hneg := neg_lt_neg h
      simpa only [neg_neg] using hneg
    have hmul : -log (1 - theta) * T < H :=
      (lt_div_iff₀ hT).mp hdiv
    simpa only [mul_comm] using hmul
  · intro h
    have hmul : -log (1 - theta) * T < H := by
      simpa only [mul_comm] using h
    have hdiv : -log (1 - theta) < H / T :=
      (lt_div_iff₀ hT).mpr hmul
    have hneg := neg_lt_neg hdiv
    simpa only [neg_neg] using hneg

/-- The conventional half-affinity threshold has cut scale `log 2`. -/
theorem cutScale_half : cutScale (1 / 2 : ℝ) = log 2 := by
  rw [cutScale, show (1 : ℝ) - 1 / 2 = 1 / 2 by norm_num,
    one_div, Real.log_inv, neg_neg]

/-- Half-affinity is the familiar `H > T * log 2` slice. -/
theorem affinity_gt_half_iff (H T : ℝ) (hT : 0 < T) :
    affinity H T > 1 / 2 ↔ H > T * log 2 := by
  simpa only [cutScale_half] using
    affinity_gt_iff H T (1 / 2) hT (by norm_num)

end CumulativeField
