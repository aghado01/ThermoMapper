import Enthymemata.ActiveSubspace.ContinuousPoincare
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.MeasureTheory.Integral.DivergenceTheorem
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Probability.ConditionalProbability
import Mathlib.Probability.Distributions.Gaussian.Real

/-!
# Active subspaces — weighted Poincaré on an interval (AS-2(e), AS-2(g))

The truncated normal's constant in E8's matrix check, `min(σ², ((b − a)/π)²)`, from two bounds for
a law `∝ e^{−V}` on `[a, b]`:
* Bakry–Émery: with `V″ ≥ κ > 0`, `Var h ≤ κ⁻¹ E[h′²]`;
* Payne–Weinberger: with `V″ ≥ 0`, a log-concave law, `Var h ≤ ((b − a)/π)² E[h′²]`.

Both proofs are one-dimensional and use no ODE theory. They pass through the dual variable: with
`g = h − m` (`m` the weighted mean), `W(x) = ∫ₐˣ g e^{−V}` vanishes at both ends, and
`W′² e^{V} = g² e^{−V}`.
* `exists_dual`: `W`.
* `poincare_of_dual`, the duality: an inequality with zero ends for `W`,
  `∫ W² e^{V} ≤ C ∫ g² e^{−V}`, gives `∫ g² e^{−V} ≤ C ∫ h′² e^{−V}`. By parts, kinks allowed,
  `∫ g² e^{−V} = −∫ W h′`; and pointwise `0 ≤ (W e^{V} + C h′)² e^{−V}`.
* `integral_bochner`: the identity through `Ψ = V′ W² e^{V}`, which vanishes at both ends,
  `∫ V″ W² e^{V} + ∫ (2 V′ W g + V′² W² e^{V}) = 0`.
* `poincare_weighted`, Bakry–Émery: `κ ∫ W² e^{V} ≤ ∫ V″ W² e^{V} ≤ ∫ g² e^{−V}`, the second
  since the difference is `∫ (g + V′ W e^{V})² e^{−V}`.
* `poincare_logConcave`, Payne–Weinberger: `z = e^{V/2} W` is zero at both ends, and
  `∫ z′² = ∫ g² e^{−V} − ½ ∫ V″ W² e^{V} − ¼ ∫ V′² W² e^{V} ≤ ∫ g² e^{−V}`. Wirtinger with
  zero ends (`wirtinger_dirichlet`) gives `∫ W² e^{V} = ∫ z² ≤ ((b − a)/π)² ∫ z′²`.
* `poincare_truncNormal`: `V = (x − μ)²/(2σ²)`, so the constant is `min(σ², ((b − a)/π)²)`.
* `hasPoincareOn_of_density`: a weighted inequality as a law's, `HasPoincareOn ν a b c` when
  `ν ∝ e^{−V}` on `[a, b]`.
* `truncNormal μ v a b`, `gaussianReal μ v` conditioned on `[a, b]`, with
  `hasPoincareOn_truncNormal`: its constant is `min(v, ((b − a)/π)²)`, so its coordinates join
  uniform ones in `matrix_poincare_prod_box` (`ae_mem_truncNormal`).
* `matrix_poincare_box_mixed`: the matrix check on E8's laws, each input uniform or a truncated
  normal, with `Dᵢ` its constant.
-/

open MeasureTheory Set Real

namespace ActiveSubspace

section Weighted

variable {a b : ℝ}

/-- The dual variable `W(x) = ∫ₐˣ (h − m) e^{−V}`: continuous on `[a, b]`, with
`W′ = (h − m) e^{−V}` inside, and zero at both ends when `m` is the weighted mean. -/
theorem exists_dual (hab : a < b) {V : ℝ → ℝ} (hVc : Continuous V) {h : ℝ → ℝ}
    (hc : ContinuousOn h (Icc a b)) (m : ℝ) (hmean : ∫ x in a..b, (h x - m) * exp (-V x) = 0) :
    ∃ W : ℝ → ℝ, ContinuousOn W (Icc a b) ∧
      (∀ x ∈ Ioo a b, HasDerivAt W ((h x - m) * exp (-V x)) x) ∧ W a = 0 ∧ W b = 0 := by
  have hgωc : ContinuousOn (fun x => (h x - m) * exp (-V x)) (Icc a b) :=
    (hc.sub continuousOn_const).mul (Real.continuous_exp.comp hVc.neg).continuousOn
  refine ⟨fun x => ∫ t in a..x, (h t - m) * exp (-V t), ?_, fun x hx => ?_,
    intervalIntegral.integral_same, hmean⟩
  · have hp := intervalIntegral.continuousOn_primitive_interval (μ := volume) (a := a) (b := b)
      (f := fun t => (h t - m) * exp (-V t))
      (by rw [uIcc_of_le hab.le]; exact hgωc.integrableOn_Icc)
    rwa [uIcc_of_le hab.le] at hp
  · exact intervalIntegral.integral_hasDerivAt_right
      ((hgωc.mono (Icc_subset_Icc le_rfl hx.2.le)).intervalIntegrable_of_Icc hx.1.le)
      ((hgωc.mono Ioo_subset_Icc_self).stronglyMeasurableAtFilter isOpen_Ioo x hx)
      (hgωc.continuousAt (Icc_mem_nhds hx.1 hx.2))

/-- **The one-dimensional duality.** For `W` with `W′ = (h − m) e^{−V}` inside `[a, b]` and zero
at both ends, so that `W′² e^{V} = (h − m)² e^{−V}`: an inequality with zero ends for `W`,
`∫ W² e^{V} ≤ C ∫ (h − m)² e^{−V}`, gives the one with no end condition for `h`, with the same
`C`: `∫ (h − m)² e^{−V} ≤ C ∫ h′² e^{−V}`. -/
theorem poincare_of_dual (hab : a < b) {V : ℝ → ℝ} (hVc : Continuous V)
    {h h' : ℝ → ℝ} {s : Set ℝ} (hs : s.Countable) (hc : ContinuousOn h (Icc a b))
    (hd : ∀ x ∈ Ioo a b \ s, HasDerivAt h (h' x) x) (hm : Measurable h') {K : ℝ}
    (hK : ∀ x ∈ Icc a b, |h' x| ≤ K) (m : ℝ) {W : ℝ → ℝ} (hWc : ContinuousOn W (Icc a b))
    (hWd : ∀ x ∈ Ioo a b, HasDerivAt W ((h x - m) * exp (-V x)) x) (hWa : W a = 0)
    (hWb : W b = 0) {C : ℝ} (hC : 0 < C)
    (hW : ∫ x in a..b, W x ^ 2 * exp (V x) ≤
      C * ∫ x in a..b, (h x - m) ^ 2 * exp (-V x)) :
    ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) ≤ C * ∫ x in a..b, h' x ^ 2 * exp (-V x) := by
  have hωc : Continuous fun x => exp (-V x) := Real.continuous_exp.comp hVc.neg
  have heVc : Continuous fun x => exp (V x) := Real.continuous_exp.comp hVc
  have hωe : ∀ x, exp (-V x) * exp (V x) = 1 := fun x => by
    rw [← Real.exp_add, neg_add_cancel, Real.exp_zero]
  have hgc : ContinuousOn (fun x => h x - m) (Icc a b) := hc.sub continuousOn_const
  have hIoc : ∀ x ∈ uIoc a b, x ∈ Icc a b := fun x hx => by
    rw [uIoc_of_le hab.le] at hx
    exact Ioc_subset_Icc_self hx
  have hh'i : IntervalIntegrable h' volume a b :=
    (intervalIntegrable_const (c := K)).mono_fun' hm.aestronglyMeasurable
      (ae_restrict_of_forall_mem measurableSet_uIoc fun x hx => by
        simpa [Real.norm_eq_abs] using hK x (hIoc x hx))
  have hh'2i : IntervalIntegrable (fun x => h' x ^ 2) volume a b :=
    (intervalIntegrable_const (c := K ^ 2)).mono_fun' (hm.pow_const 2).aestronglyMeasurable
      (ae_restrict_of_forall_mem measurableSet_uIoc fun x hx => by
        simp only [Real.norm_eq_abs, abs_pow]
        exact pow_le_pow_left₀ (abs_nonneg _) (hK x (hIoc x hx)) 2)
  have hiWh : IntervalIntegrable (fun x => W x * h' x) volume a b :=
    hh'i.continuousOn_mul (by rw [uIcc_of_le hab.le]; exact hWc)
  have hih2 : IntervalIntegrable (fun x => h' x ^ 2 * exp (-V x)) volume a b :=
    hh'2i.mul_continuousOn hωc.continuousOn
  have hig2 : IntervalIntegrable (fun x => (h x - m) ^ 2 * exp (-V x)) volume a b :=
    ((hgc.pow 2).mul hωc.continuousOn).intervalIntegrable_of_Icc hab.le
  have hiW2 : IntervalIntegrable (fun x => W x ^ 2 * exp (V x)) volume a b :=
    ((hWc.pow 2).mul heVc.continuousOn).intervalIntegrable_of_Icc hab.le
  -- by parts, kinks allowed: `∫ (h − m)² e^{−V} = −∫ W h′`
  have hparts : ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) = -∫ x in a..b, W x * h' x := by
    have hftc := integral_eq_of_hasDerivAt_off_countable_of_le (fun x => W x * (h x - m))
      (fun x => (h x - m) ^ 2 * exp (-V x) + W x * h' x) hab.le hs (hWc.mul hgc)
      (fun x hx => ((hWd x hx.1).mul ((hd x hx).sub_const m)).congr_deriv (by ring))
      (hig2.add hiWh)
    simp only [hWa, hWb, zero_mul, sub_self] at hftc
    rw [intervalIntegral.integral_add hig2 hiWh] at hftc
    linarith
  -- `0 ≤ ∫ (W e^{V} + C h′)² e^{−V} = ∫ W² e^{V} + 2C ∫ W h′ + C² ∫ h′² e^{−V}`
  have hsq : ∀ x, (W x * exp (V x) + C * h' x) ^ 2 * exp (-V x) =
      W x ^ 2 * exp (V x) + 2 * C * (W x * h' x) + C ^ 2 * (h' x ^ 2 * exp (-V x)) := by
    intro x
    have := hωe x
    linear_combination (W x ^ 2 * exp (V x) + 2 * C * (W x * h' x)) * this
  have hnn : 0 ≤ ∫ x in a..b, (W x * exp (V x) + C * h' x) ^ 2 * exp (-V x) :=
    intervalIntegral.integral_nonneg hab.le fun x _ =>
      mul_nonneg (sq_nonneg _) (Real.exp_pos _).le
  rw [intervalIntegral.integral_congr (fun x _ => hsq x),
    intervalIntegral.integral_add (hiW2.add (hiWh.const_mul (2 * C))) (hih2.const_mul (C ^ 2)),
    intervalIntegral.integral_add hiW2 (hiWh.const_mul (2 * C)),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul] at hnn
  -- `0 ≤ C G − 2 C G + C² H`
  have hCI : C * ∫ x in a..b, W x * h' x =
      -(C * ∫ x in a..b, (h x - m) ^ 2 * exp (-V x)) := by
    rw [hparts]
    ring
  have hCG : C * ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) ≤
      C * (C * ∫ x in a..b, h' x ^ 2 * exp (-V x)) := by
    nlinarith [hnn, hW, hCI]
  exact le_of_mul_le_mul_left hCG hC

/-- **The Bochner step,** as one identity through `Ψ = V′ W² e^{V}`, which vanishes at both
ends: `∫ V″ W² e^{V} + ∫ (2 V′ W (h − m) + V′² W² e^{V}) = 0` for `W` as in `poincare_of_dual`. -/
theorem integral_bochner (hab : a < b) {V V' V'' : ℝ → ℝ} (hV : ∀ x, HasDerivAt V (V' x) x)
    (hV' : ∀ x, HasDerivAt V' (V'' x) x) (hV''c : Continuous V'') {h : ℝ → ℝ}
    (hc : ContinuousOn h (Icc a b)) (m : ℝ) {W : ℝ → ℝ} (hWc : ContinuousOn W (Icc a b))
    (hWd : ∀ x ∈ Ioo a b, HasDerivAt W ((h x - m) * exp (-V x)) x) (hWa : W a = 0)
    (hWb : W b = 0) :
    (∫ x in a..b, V'' x * W x ^ 2 * exp (V x)) +
      ∫ x in a..b, (2 * V' x * W x * (h x - m) + V' x ^ 2 * W x ^ 2 * exp (V x)) = 0 := by
  have hVc : Continuous V := continuous_iff_continuousAt.2 fun x => (hV x).continuousAt
  have hV'c : Continuous V' := continuous_iff_continuousAt.2 fun x => (hV' x).continuousAt
  have heVc : Continuous fun x => exp (V x) := Real.continuous_exp.comp hVc
  have hωe : ∀ x, exp (-V x) * exp (V x) = 1 := fun x => by
    rw [← Real.exp_add, neg_add_cancel, Real.exp_zero]
  have hgc : ContinuousOn (fun x => h x - m) (Icc a b) := hc.sub continuousOn_const
  have hiQ : IntervalIntegrable (fun x => V'' x * W x ^ 2 * exp (V x)) volume a b :=
    ((hV''c.continuousOn.mul (hWc.pow 2)).mul heVc.continuousOn).intervalIntegrable_of_Icc hab.le
  have hiR : IntervalIntegrable
      (fun x => 2 * V' x * W x * (h x - m) + V' x ^ 2 * W x ^ 2 * exp (V x)) volume a b :=
    ((((continuousOn_const.mul hV'c.continuousOn).mul hWc).mul hgc).add
      (((hV'c.continuousOn.pow 2).mul (hWc.pow 2)).mul heVc.continuousOn))
      |>.intervalIntegrable_of_Icc hab.le
  have hftc := integral_eq_of_hasDerivAt_off_countable_of_le
    (fun x => V' x * W x ^ 2 * exp (V x))
    (fun x => V'' x * W x ^ 2 * exp (V x) +
      (2 * V' x * W x * (h x - m) + V' x ^ 2 * W x ^ 2 * exp (V x)))
    hab.le (s := ∅) countable_empty
    ((hV'c.continuousOn.mul (hWc.pow 2)).mul heVc.continuousOn)
    (fun x hx => (((hV' x).mul ((hWd x hx.1).fun_pow 2)).mul (hV x).exp).congr_deriv (by
      have := hωe x
      simp only [Pi.mul_apply]
      linear_combination (2 * V' x * W x * (h x - m)) * this))
    (hiQ.add hiR)
  simp only [hWa, hWb, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero,
    zero_mul, sub_self] at hftc
  rw [intervalIntegral.integral_add hiQ hiR] at hftc
  exact hftc

/-- **AS-2(e), Bakry–Émery on an interval.** For `V` twice differentiable with `V″` continuous and
at least `κ > 0` on `[a, b]`, and `h` continuous on `[a, b]` with a measurable derivative off a
countable set, bounded on `[a, b]`: if `m` is the `e^{−V}`-weighted mean, `∫ (h − m) e^{−V} = 0`,
then `∫ (h − m)² e^{−V} ≤ κ⁻¹ ∫ h′² e^{−V}`. -/
theorem poincare_weighted (hab : a < b) {V V' V'' : ℝ → ℝ} {κ : ℝ} (hκ : 0 < κ)
    (hV : ∀ x, HasDerivAt V (V' x) x) (hV' : ∀ x, HasDerivAt V' (V'' x) x)
    (hV''c : Continuous V'') (hV''κ : ∀ x ∈ Icc a b, κ ≤ V'' x)
    {h h' : ℝ → ℝ} {s : Set ℝ} (hs : s.Countable) (hc : ContinuousOn h (Icc a b))
    (hd : ∀ x ∈ Ioo a b \ s, HasDerivAt h (h' x) x) (hm : Measurable h') {K : ℝ}
    (hK : ∀ x ∈ Icc a b, |h' x| ≤ K) (m : ℝ) (hmean : ∫ x in a..b, (h x - m) * exp (-V x) = 0) :
    ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) ≤ κ⁻¹ * ∫ x in a..b, h' x ^ 2 * exp (-V x) := by
  have hVc : Continuous V := continuous_iff_continuousAt.2 fun x => (hV x).continuousAt
  have hV'c : Continuous V' := continuous_iff_continuousAt.2 fun x => (hV' x).continuousAt
  have hωc : Continuous fun x => exp (-V x) := Real.continuous_exp.comp hVc.neg
  have heVc : Continuous fun x => exp (V x) := Real.continuous_exp.comp hVc
  have hωe : ∀ x, exp (-V x) * exp (V x) = 1 := fun x => by
    rw [← Real.exp_add, neg_add_cancel, Real.exp_zero]
  obtain ⟨W, hWc, hWd, hWa, hWb⟩ := exists_dual hab hVc hc m hmean
  have hΨ := integral_bochner hab hV hV' hV''c hc m hWc hWd hWa hWb
  have hgc : ContinuousOn (fun x => h x - m) (Icc a b) := hc.sub continuousOn_const
  have hig2 : IntervalIntegrable (fun x => (h x - m) ^ 2 * exp (-V x)) volume a b :=
    ((hgc.pow 2).mul hωc.continuousOn).intervalIntegrable_of_Icc hab.le
  have hiQ : IntervalIntegrable (fun x => V'' x * W x ^ 2 * exp (V x)) volume a b :=
    ((hV''c.continuousOn.mul (hWc.pow 2)).mul heVc.continuousOn).intervalIntegrable_of_Icc hab.le
  have hiR : IntervalIntegrable
      (fun x => 2 * V' x * W x * (h x - m) + V' x ^ 2 * W x ^ 2 * exp (V x)) volume a b :=
    ((((continuousOn_const.mul hV'c.continuousOn).mul hWc).mul hgc).add
      (((hV'c.continuousOn.pow 2).mul (hWc.pow 2)).mul heVc.continuousOn))
      |>.intervalIntegrable_of_Icc hab.le
  have hiW2 : IntervalIntegrable (fun x => W x ^ 2 * exp (V x)) volume a b :=
    ((hWc.pow 2).mul heVc.continuousOn).intervalIntegrable_of_Icc hab.le
  -- `∫ (h − m)² e^{−V} − ∫ V″ W² e^{V} = ∫ (h − m + V′ W e^{V})² e^{−V} ≥ 0`
  have hsq : ∀ x, (h x - m + V' x * W x * exp (V x)) ^ 2 * exp (-V x) =
      (h x - m) ^ 2 * exp (-V x) +
        (2 * V' x * W x * (h x - m) + V' x ^ 2 * W x ^ 2 * exp (V x)) := by
    intro x
    have := hωe x
    linear_combination (2 * V' x * W x * (h x - m) + V' x ^ 2 * W x ^ 2 * exp (V x)) * this
  have hnn : 0 ≤ ∫ x in a..b, (h x - m + V' x * W x * exp (V x)) ^ 2 * exp (-V x) :=
    intervalIntegral.integral_nonneg hab.le fun x _ =>
      mul_nonneg (sq_nonneg _) (Real.exp_pos _).le
  rw [intervalIntegral.integral_congr (fun x _ => hsq x),
    intervalIntegral.integral_add hig2 hiR] at hnn
  -- `κ ∫ W² e^{V} ≤ ∫ V″ W² e^{V}`
  have hκW : κ * ∫ x in a..b, W x ^ 2 * exp (V x) ≤
      ∫ x in a..b, V'' x * W x ^ 2 * exp (V x) := by
    rw [← intervalIntegral.integral_const_mul]
    refine intervalIntegral.integral_mono_on hab.le (hiW2.const_mul κ) hiQ fun x hx => ?_
    calc κ * (W x ^ 2 * exp (V x)) ≤ V'' x * (W x ^ 2 * exp (V x)) :=
          mul_le_mul_of_nonneg_right (hV''κ x hx) (mul_nonneg (sq_nonneg _) (Real.exp_pos _).le)
      _ = V'' x * W x ^ 2 * exp (V x) := by ring
  refine poincare_of_dual hab hVc hs hc hd hm hK m hWc hWd hWa hWb (inv_pos.2 hκ) ?_
  calc ∫ x in a..b, W x ^ 2 * exp (V x)
      = κ⁻¹ * (κ * ∫ x in a..b, W x ^ 2 * exp (V x)) := by
        rw [← mul_assoc, inv_mul_cancel₀ hκ.ne', one_mul]
    _ ≤ κ⁻¹ * ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) :=
        mul_le_mul_of_nonneg_left (by linarith) (inv_pos.2 hκ).le

/-- **AS-2(g), Payne–Weinberger on an interval.** For `V` convex on `[a, b]` (`V″ ≥ 0` there),
with `V` twice differentiable and `V″` continuous, and `h` as in `poincare_weighted`: if `m` is the
`e^{−V}`-weighted mean, then `∫ (h − m)² e^{−V} ≤ ((b − a)/π)² ∫ h′² e^{−V}`. A log-concave law
on `[a, b]` has the uniform law's constant, whatever `V`. -/
theorem poincare_logConcave (hab : a < b) {V V' V'' : ℝ → ℝ}
    (hV : ∀ x, HasDerivAt V (V' x) x) (hV' : ∀ x, HasDerivAt V' (V'' x) x)
    (hV''c : Continuous V'') (hV''0 : ∀ x ∈ Icc a b, 0 ≤ V'' x)
    {h h' : ℝ → ℝ} {s : Set ℝ} (hs : s.Countable) (hc : ContinuousOn h (Icc a b))
    (hd : ∀ x ∈ Ioo a b \ s, HasDerivAt h (h' x) x) (hm : Measurable h') {K : ℝ}
    (hK : ∀ x ∈ Icc a b, |h' x| ≤ K) (m : ℝ) (hmean : ∫ x in a..b, (h x - m) * exp (-V x) = 0) :
    ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) ≤
      ((b - a) / π) ^ 2 * ∫ x in a..b, h' x ^ 2 * exp (-V x) := by
  have hVc : Continuous V := continuous_iff_continuousAt.2 fun x => (hV x).continuousAt
  have hV'c : Continuous V' := continuous_iff_continuousAt.2 fun x => (hV' x).continuousAt
  have hωc : Continuous fun x => exp (-V x) := Real.continuous_exp.comp hVc.neg
  have heVc : Continuous fun x => exp (V x) := Real.continuous_exp.comp hVc
  have hEc : Continuous fun x => exp (V x / 2) := Real.continuous_exp.comp (hVc.div_const 2)
  have hωe : ∀ x, exp (-V x) * exp (V x) = 1 := fun x => by
    rw [← Real.exp_add, neg_add_cancel, Real.exp_zero]
  have hE2 : ∀ x, exp (V x / 2) ^ 2 = exp (V x) := fun x => by
    rw [sq, ← Real.exp_add, add_halves]
  obtain ⟨W, hWc, hWd, hWa, hWb⟩ := exists_dual hab hVc hc m hmean
  have hΨ := integral_bochner hab hV hV' hV''c hc m hWc hWd hWa hWb
  have hgc : ContinuousOn (fun x => h x - m) (Icc a b) := hc.sub continuousOn_const
  -- `z = e^{V/2} W`, zero at both ends
  set z : ℝ → ℝ := fun x => exp (V x / 2) * W x with hzdef
  set z' : ℝ → ℝ := fun x => exp (V x / 2) * (V' x / 2 * W x + (h x - m) * exp (-V x))
    with hz'def
  have hzd : ∀ x ∈ Ioo a b \ (∅ : Set ℝ), HasDerivAt z (z' x) x := fun x hx =>
    (((hV x).div_const 2).exp.mul (hWd x hx.1)).congr_deriv (by simp only [hz'def]; ring)
  have hzc : ContinuousOn z (Icc a b) := hEc.continuousOn.mul hWc
  have hz'c : ContinuousOn z' (Icc a b) :=
    hEc.continuousOn.mul (((hV'c.div_const 2).continuousOn.mul hWc).add
      (hgc.mul hωc.continuousOn))
  have hD := wirtinger_dirichlet hab countable_empty hzc hzd hz'c
    (by simp only [hzdef, hWa, mul_zero]) (by simp only [hzdef, hWb, mul_zero])
  -- `∫ z² = ∫ W² e^{V}`
  have hz2 : ∀ x, z x ^ 2 = W x ^ 2 * exp (V x) := fun x => by
    simp only [hzdef]
    rw [mul_pow, hE2]
    ring
  -- `∫ z′² = ∫ (h − m)² e^{−V} + ½ ∫ (2 V′ W (h − m) + V′² W² e^{V}) − ¼ ∫ V′² W² e^{V}`
  have hz'2 : ∀ x, z' x ^ 2 = (h x - m) ^ 2 * exp (-V x) +
      (1 / 2) * (2 * V' x * W x * (h x - m) + V' x ^ 2 * W x ^ 2 * exp (V x)) -
      (1 / 4) * (V' x ^ 2 * W x ^ 2 * exp (V x)) := fun x => by
    simp only [hz'def]
    have h1 := hE2 x
    have h2 := hωe x
    linear_combination (V' x / 2 * W x + (h x - m) * exp (-V x)) ^ 2 * h1 +
      (V' x * W x * (h x - m) + (h x - m) ^ 2 * exp (-V x)) * h2
  have hig2 : IntervalIntegrable (fun x => (h x - m) ^ 2 * exp (-V x)) volume a b :=
    ((hgc.pow 2).mul hωc.continuousOn).intervalIntegrable_of_Icc hab.le
  have hiR : IntervalIntegrable
      (fun x => 2 * V' x * W x * (h x - m) + V' x ^ 2 * W x ^ 2 * exp (V x)) volume a b :=
    ((((continuousOn_const.mul hV'c.continuousOn).mul hWc).mul hgc).add
      (((hV'c.continuousOn.pow 2).mul (hWc.pow 2)).mul heVc.continuousOn))
      |>.intervalIntegrable_of_Icc hab.le
  have hiP : IntervalIntegrable (fun x => V' x ^ 2 * W x ^ 2 * exp (V x)) volume a b :=
    (((hV'c.continuousOn.pow 2).mul (hWc.pow 2)).mul heVc.continuousOn).intervalIntegrable_of_Icc
      hab.le
  have hint : ∫ x in a..b, z' x ^ 2 = (∫ x in a..b, (h x - m) ^ 2 * exp (-V x)) +
      (1 / 2) * (∫ x in a..b, (2 * V' x * W x * (h x - m) + V' x ^ 2 * W x ^ 2 * exp (V x))) -
      (1 / 4) * ∫ x in a..b, V' x ^ 2 * W x ^ 2 * exp (V x) := by
    rw [intervalIntegral.integral_congr (fun x _ => hz'2 x),
      intervalIntegral.integral_sub (hig2.add (hiR.const_mul (1 / 2))) (hiP.const_mul (1 / 4)),
      intervalIntegral.integral_add hig2 (hiR.const_mul (1 / 2)),
      intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul]
  have hQnn : 0 ≤ ∫ x in a..b, V'' x * W x ^ 2 * exp (V x) :=
    intervalIntegral.integral_nonneg hab.le fun x hx =>
      mul_nonneg (mul_nonneg (hV''0 x hx) (sq_nonneg _)) (Real.exp_pos _).le
  have hPnn : 0 ≤ ∫ x in a..b, V' x ^ 2 * W x ^ 2 * exp (V x) :=
    intervalIntegral.integral_nonneg hab.le fun x _ =>
      mul_nonneg (mul_nonneg (sq_nonneg _) (sq_nonneg _)) (Real.exp_pos _).le
  have hz'le : ∫ x in a..b, z' x ^ 2 ≤ ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) := by
    rw [hint]
    linarith
  rw [intervalIntegral.integral_congr (fun x _ => hz2 x)] at hD
  have hC : 0 < ((b - a) / π) ^ 2 := pow_pos (div_pos (sub_pos.2 hab) Real.pi_pos) 2
  refine poincare_of_dual hab hVc hs hc hd hm hK m hWc hWd hWa hWb hC ?_
  exact hD.trans (mul_le_mul_of_nonneg_left hz'le hC.le)

/-- **AS-2(e) and (g), the truncated normal.** For the normal `N(μ, σ²)` restricted to `[a, b]`
(`V = (x − μ)²/(2σ²)`, so `V″ = σ⁻²`), the weighted variance is at most
`min(σ², ((b − a)/π)²)` times the weighted mean squared derivative: Bakry–Émery's `σ²`, or
Payne–Weinberger's `((b − a)/π)²`. -/
theorem poincare_truncNormal (hab : a < b) {μ σ : ℝ} (hσ : 0 < σ) {h h' : ℝ → ℝ} {s : Set ℝ}
    (hs : s.Countable) (hc : ContinuousOn h (Icc a b))
    (hd : ∀ x ∈ Ioo a b \ s, HasDerivAt h (h' x) x) (hm : Measurable h') {K : ℝ}
    (hK : ∀ x ∈ Icc a b, |h' x| ≤ K) (m : ℝ)
    (hmean : ∫ x in a..b, (h x - m) * exp (-((x - μ) ^ 2 / (2 * σ ^ 2))) = 0) :
    ∫ x in a..b, (h x - m) ^ 2 * exp (-((x - μ) ^ 2 / (2 * σ ^ 2))) ≤
      min (σ ^ 2) (((b - a) / π) ^ 2) *
        ∫ x in a..b, h' x ^ 2 * exp (-((x - μ) ^ 2 / (2 * σ ^ 2))) := by
  have hσ2 : 0 < σ ^ 2 := pow_pos hσ 2
  have hV : ∀ x, HasDerivAt (fun x => (x - μ) ^ 2 / (2 * σ ^ 2)) ((x - μ) / σ ^ 2) x := by
    intro x
    have h1 : HasDerivAt (fun x => (x - μ) ^ 2) (2 * (x - μ)) x := by
      simpa using ((hasDerivAt_id x).sub_const μ).fun_pow 2
    refine (h1.div_const (2 * σ ^ 2)).congr_deriv ?_
    field_simp
    try ring
  have hV' : ∀ x, HasDerivAt (fun x => (x - μ) / σ ^ 2) (1 / σ ^ 2) x := by
    intro x
    simpa using ((hasDerivAt_id x).sub_const μ).div_const (σ ^ 2)
  rcases min_choice (σ ^ 2) (((b - a) / π) ^ 2) with hmin | hmin <;> rw [hmin]
  · have h := poincare_weighted hab (κ := 1 / σ ^ 2) (by positivity) hV hV' continuous_const
      (fun _ _ => le_rfl) hs hc hd hm hK m hmean
    rwa [one_div, inv_inv] at h
  · exact poincare_logConcave hab hV hV' continuous_const (fun _ _ => (one_div_pos.2 hσ2).le) hs
      hc hd hm hK m hmean

end Weighted

section Laws

open ProbabilityTheory
open scoped NNReal

variable {a b : ℝ}

/-- **A weighted inequality as a law's.** A probability law `ν` with `∫ φ dν = C ∫ₐᵇ φ e^{−V}` for
every `φ` (`ν ∝ e^{−V}` on `[a, b]`) has the one-variable constant `c` of any weighted inequality
`∫ (h − m)² e^{−V} ≤ c ∫ h′² e^{−V}` that holds at the weighted mean `m`: Bakry–Émery's `κ⁻¹`
(`poincare_weighted`), or Payne–Weinberger's `((b − a)/π)²` (`poincare_logConcave`). -/
theorem hasPoincareOn_of_density (hab : a < b) {V : ℝ → ℝ} (hVc : Continuous V) {c : ℝ}
    (hw : ∀ ⦃h h' : ℝ → ℝ⦄ ⦃s : Set ℝ⦄ ⦃K : ℝ⦄, s.Countable → ContinuousOn h (Icc a b) →
      (∀ x ∈ Ioo a b \ s, HasDerivAt h (h' x) x) → Measurable h' →
      (∀ x ∈ Icc a b, |h' x| ≤ K) → ∀ m : ℝ, ∫ x in a..b, (h x - m) * exp (-V x) = 0 →
      ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) ≤ c * ∫ x in a..b, h' x ^ 2 * exp (-V x))
    (ν : Measure ℝ) [IsProbabilityMeasure ν] {C : ℝ}
    (hν : ∀ φ : ℝ → ℝ, ∫ x, φ x ∂ν = C * ∫ x in a..b, φ x * exp (-V x)) :
    HasPoincareOn ν a b c := by
  intro h h' s K hs hc hd hm hK
  have hωc : Continuous fun x => exp (-V x) := Real.continuous_exp.comp hVc.neg
  have hZ : 1 = C * ∫ x in a..b, exp (-V x) := by
    have h1 := hν fun _ => 1
    simp only [integral_const, probReal_univ, smul_eq_mul, mul_one, one_mul] at h1
    exact h1
  have hZnn : 0 ≤ ∫ x in a..b, exp (-V x) :=
    intervalIntegral.integral_nonneg hab.le fun x _ => (Real.exp_pos _).le
  have hC0 : 0 ≤ C := by nlinarith
  set m := ∫ x, h x ∂ν with hmdef
  have hhi : IntervalIntegrable (fun x => h x * exp (-V x)) volume a b :=
    (hc.mul hωc.continuousOn).intervalIntegrable_of_Icc hab.le
  have hmean : ∫ x in a..b, (h x - m) * exp (-V x) = 0 := by
    have e1 : ∀ x, (h x - m) * exp (-V x) = h x * exp (-V x) - m * exp (-V x) := fun x => by
      ring
    simp only [e1]
    rw [intervalIntegral.integral_sub hhi ((hωc.intervalIntegrable a b).const_mul m),
      intervalIntegral.integral_const_mul]
    have hm' : m = C * ∫ x in a..b, h x * exp (-V x) := hν h
    rw [hm']
    linear_combination (∫ x in a..b, h x * exp (-V x)) * hZ
  have hW := hw hs hc hd hm (fun x _ => hK x) m hmean
  have hvar : varM ν h = C * ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) :=
    hν fun x => (h x - m) ^ 2
  have hen : ∫ x, h' x ^ 2 ∂ν = C * ∫ x in a..b, h' x ^ 2 * exp (-V x) :=
    hν fun x => h' x ^ 2
  rw [hvar, hen]
  calc C * ∫ x in a..b, (h x - m) ^ 2 * exp (-V x)
      ≤ C * (c * ∫ x in a..b, h' x ^ 2 * exp (-V x)) := mul_le_mul_of_nonneg_left hW hC0
    _ = c * (C * ∫ x in a..b, h' x ^ 2 * exp (-V x)) := by ring

/-- The normal law `N(μ, v)` truncated to `[a, b]`: `gaussianReal μ v` conditioned on `[a, b]`. -/
noncomputable def truncNormal (μ : ℝ) (v : ℝ≥0) (a b : ℝ) : Measure ℝ :=
  ProbabilityTheory.cond (gaussianReal μ v) (Icc a b)

theorem truncNormal_isProbabilityMeasure (hab : a < b) (μ : ℝ) {v : ℝ≥0} (hv : v ≠ 0) :
    IsProbabilityMeasure (truncNormal μ v a b) := by
  unfold truncNormal
  refine cond_isProbabilityMeasure fun h0 => ?_
  have h1 := gaussianReal_absolutelyContinuous' μ hv h0
  rw [Real.volume_Icc, ENNReal.ofReal_eq_zero] at h1
  linarith

/-- The truncated normal lives on `[a, b]`. -/
theorem ae_mem_truncNormal (μ : ℝ) (v : ℝ≥0) (a b : ℝ) :
    ∀ᵐ t ∂(truncNormal μ v a b), t ∈ Icc a b := by
  unfold truncNormal ProbabilityTheory.cond
  exact Measure.ae_smul_measure (ae_restrict_mem measurableSet_Icc) _

/-- Integrals under the truncated normal, as `e^{−(x−μ)²/(2v)}`-weighted integrals over `[a, b]`. -/
theorem integral_truncNormal (hab : a < b) (μ : ℝ) {v : ℝ≥0} (hv : v ≠ 0) (φ : ℝ → ℝ) :
    ∫ x, φ x ∂(truncNormal μ v a b) =
      (gaussianReal μ v (Icc a b)).toReal⁻¹ * (Real.sqrt (2 * π * v))⁻¹ *
        ∫ x in a..b, φ x * exp (-((x - μ) ^ 2 / (2 * (v : ℝ)))) := by
  have hpdf : ∀ x, (gaussianPDF μ v x).toReal * φ x =
      (Real.sqrt (2 * π * v))⁻¹ * (φ x * exp (-((x - μ) ^ 2 / (2 * (v : ℝ))))) := fun x => by
    rw [toReal_gaussianPDF, gaussianPDFReal, neg_div]
    ring
  rw [truncNormal, ProbabilityTheory.cond, integral_smul_measure, ENNReal.toReal_inv,
    gaussianReal_of_var_ne_zero μ hv,
    setIntegral_withDensity_eq_setIntegral_toReal_smul (measurable_gaussianPDF μ v)
      (ae_of_all _ fun _ => ENNReal.ofReal_lt_top) φ measurableSet_Icc,
    integral_Icc_eq_integral_Ioc, ← intervalIntegral.integral_of_le hab.le]
  simp only [smul_eq_mul, hpdf, intervalIntegral.integral_const_mul]
  ring

/-- **AS-2(e) and (g), the truncated normal as a law.** `N(μ, v)` truncated to `[a, b]` has the
one-variable Poincaré constant `min(v, ((b − a)/π)²)` on `[a, b]`: Bakry–Émery's `v`
(`V″ = 1/v`), or Payne–Weinberger's `((b − a)/π)²` (`V` convex). -/
theorem hasPoincareOn_truncNormal (hab : a < b) (μ : ℝ) {v : ℝ≥0} (hv : v ≠ 0) :
    HasPoincareOn (truncNormal μ v a b) a b (min (v : ℝ) (((b - a) / π) ^ 2)) := by
  haveI := truncNormal_isProbabilityMeasure hab μ hv
  have hv0 : (0 : ℝ) < v := NNReal.coe_pos.2 (pos_iff_ne_zero.2 hv)
  have hV : ∀ x, HasDerivAt (fun x => (x - μ) ^ 2 / (2 * (v : ℝ))) ((x - μ) / v) x := by
    intro x
    have h1 : HasDerivAt (fun x => (x - μ) ^ 2) (2 * (x - μ)) x := by
      simpa using ((hasDerivAt_id x).sub_const μ).fun_pow 2
    refine (h1.div_const (2 * (v : ℝ))).congr_deriv ?_
    field_simp
    try ring
  have hV' : ∀ x, HasDerivAt (fun x => (x - μ) / (v : ℝ)) (1 / v) x := by
    intro x
    simpa using ((hasDerivAt_id x).sub_const μ).div_const (v : ℝ)
  have hVc : Continuous fun x => (x - μ) ^ 2 / (2 * (v : ℝ)) :=
    continuous_iff_continuousAt.2 fun x => (hV x).continuousAt
  refine hasPoincareOn_of_density hab hVc ?_ (truncNormal μ v a b) (integral_truncNormal hab μ hv)
  intro h h' s K hs hc hd hm hK m hmean
  rcases min_choice (v : ℝ) (((b - a) / π) ^ 2) with hmin | hmin <;> rw [hmin]
  · have h := poincare_weighted hab (κ := 1 / v) (one_div_pos.2 hv0) hV hV' continuous_const
      (fun _ _ => le_rfl) hs hc hd hm hK m hmean
    rwa [one_div, inv_inv] at h
  · exact poincare_logConcave hab hV hV' continuous_const (fun _ _ => (one_div_pos.2 hv0).le) hs
      hc hd hm hK m hmean

/-- **AS-2(d), (e) and (g), on E8's laws.** Each input is uniform on `[aᵢ, bᵢ]`, or, for
`i ∈ N`, the normal `N(μᵢ, vᵢ)` truncated to it. Then `Cov(f) ≼ E[J D Jᵀ]` with
`Dᵢ = ((bᵢ − aᵢ)/π)²` for a uniform input and `Dᵢ = min(vᵢ, ((bᵢ − aᵢ)/π)²)` for a truncated
normal, every hypothesis asked only on the box. -/
theorem matrix_poincare_box_mixed {o : Type*} [Finite o] {n : ℕ} (a b : Fin n → ℝ)
    (hab : ∀ i, a i < b i) (N : Finset (Fin n)) (μ : Fin n → ℝ) (v : Fin n → ℝ≥0)
    (hv : ∀ i ∈ N, v i ≠ 0)
    {f : (Fin n → ℝ) → o → ℝ} (hf : ∀ p, Measurable fun x => f x p) {C : ℝ}
    (hC : ∀ x, (∀ i, x i ∈ Icc (a i) (b i)) → ∀ p, |f x p| ≤ C)
    (J : (Fin n → ℝ) → Matrix o (Fin n) ℝ) (hJ : ∀ p i, Measurable fun x => J x p i) {K : ℝ}
    (hK : ∀ x, (∀ i, x i ∈ Icc (a i) (b i)) → ∀ p i, |J x p i| ≤ K)
    (s : Fin n → (Fin n → ℝ) → Set ℝ) (hs : ∀ i x, (s i x).Countable)
    (hcont : ∀ p i x, (∀ j, x j ∈ Icc (a j) (b j)) →
      ContinuousOn (fun t => f (Function.update x i t) p) (Icc (a i) (b i)))
    (hder : ∀ p i x, (∀ j, x j ∈ Icc (a j) (b j)) → ∀ t ∈ Ioo (a i) (b i) \ s i x,
      HasDerivAt (fun t => f (Function.update x i t) p) (J (Function.update x i t) p i) t) :
    (energyM (Measure.pi fun i =>
        if i ∈ N then truncNormal (μ i) (v i) (a i) (b i) else uniformIcc (a i) (b i)) J
        (fun i => if i ∈ N then min (v i : ℝ) (((b i - a i) / π) ^ 2)
          else ((b i - a i) / π) ^ 2) -
      covM (Measure.pi fun i =>
        if i ∈ N then truncNormal (μ i) (v i) (a i) (b i) else uniformIcc (a i) (b i)) f
      ).PosSemidef := by
  haveI : ∀ i, IsProbabilityMeasure
      (if i ∈ N then truncNormal (μ i) (v i) (a i) (b i) else uniformIcc (a i) (b i)) := by
    intro i
    by_cases hi : i ∈ N
    · simp only [hi, ↓reduceIte]
      exact truncNormal_isProbabilityMeasure (hab i) (μ i) (hv i hi)
    · simp only [hi, ↓reduceIte]
      exact uniformIcc_isProbabilityMeasure (hab i)
  refine matrix_poincare_prod_box _ a b _ hab (fun i => ?_) (fun i => ?_) hf hC J hJ hK s hs
    hcont hder
  · by_cases hi : i ∈ N
    · simp only [hi, ↓reduceIte]
      exact hasPoincareOn_truncNormal (hab i) (μ i) (hv i hi)
    · simp only [hi, ↓reduceIte]
      exact hasPoincareOn_uniformIcc (hab i)
  · by_cases hi : i ∈ N
    · simp only [hi, ↓reduceIte]
      exact ae_mem_truncNormal (μ i) (v i) (a i) (b i)
    · simp only [hi, ↓reduceIte]
      exact ae_mem_uniformIcc (a i) (b i)

end Laws

end ActiveSubspace
