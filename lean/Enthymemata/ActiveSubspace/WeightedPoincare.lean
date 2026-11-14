import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.MeasureTheory.Integral.DivergenceTheorem
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# Active subspaces — Bakry–Émery on an interval (AS-2(e))

The truncated normal's constant in E8's matrix check, `σ²`, as Bakry–Émery's bound for a law
`∝ e^{−V}` on `[a, b]` with `V″ ≥ κ > 0`: `Var h ≤ κ⁻¹ E[h′²]`.

The proof is one-dimensional and uses no ODE theory. With `g = h − m` (`m` the weighted mean) and
`W(x) = ∫ₐˣ g e^{−V}`, which vanishes at both ends:
* integration by parts, kinks allowed: `∫ g² e^{−V} = −∫ W h′`;
* the Bochner step, as one identity for `Ψ = V′ W² e^{V}`:
  `∫ g² e^{−V} − ∫ V″ W² e^{V} = ∫ (g + V′ W e^{V})² e^{−V} ≥ 0`;
* pointwise, `−W h′ ≤ ½ V″ W² e^{V} + ½ h′² e^{−V} / V″`.
Together, `∫ g² e^{−V} ≤ κ⁻¹ ∫ h′² e^{−V}`.

* `poincare_weighted`: the inequality, with the weighted mean given by `∫ (h − m) e^{−V} = 0`.
* `poincare_truncNormal`: `V = (x − μ)²/(2σ²)`, so the constant is `σ²`.
-/

open MeasureTheory Set Real

namespace ActiveSubspace

section Weighted

variable {a b : ℝ}

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
  set g : ℝ → ℝ := fun x => h x - m
  have hgc : ContinuousOn g (Icc a b) := hc.sub continuousOn_const
  have hgωc : ContinuousOn (fun x => g x * exp (-V x)) (Icc a b) := hgc.mul hωc.continuousOn
  -- the primitive `W`
  set W : ℝ → ℝ := fun x => ∫ t in a..x, g t * exp (-V t)
  have hWc : ContinuousOn W (Icc a b) := by
    have h := intervalIntegral.continuousOn_primitive_interval (μ := volume) (a := a) (b := b)
      (f := fun t => g t * exp (-V t)) (by rw [uIcc_of_le hab.le]; exact hgωc.integrableOn_Icc)
    rwa [uIcc_of_le hab.le] at h
  have hWd : ∀ x ∈ Ioo a b, HasDerivAt W (g x * exp (-V x)) x := fun x hx =>
    intervalIntegral.integral_hasDerivAt_right
      ((hgωc.mono (Icc_subset_Icc le_rfl hx.2.le)).intervalIntegrable_of_Icc hx.1.le)
      ((hgωc.mono Ioo_subset_Icc_self).stronglyMeasurableAtFilter isOpen_Ioo x hx)
      (hgωc.continuousAt (Icc_mem_nhds hx.1 hx.2))
  have hWa : W a = 0 := intervalIntegral.integral_same
  have hWb : W b = 0 := hmean
  -- integrability
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
  have hig2 : IntervalIntegrable (fun x => g x ^ 2 * exp (-V x)) volume a b :=
    ((hgc.pow 2).mul hωc.continuousOn).intervalIntegrable_of_Icc hab.le
  have hiQ : IntervalIntegrable (fun x => V'' x * W x ^ 2 * exp (V x)) volume a b :=
    ((hV''c.continuousOn.mul (hWc.pow 2)).mul heVc.continuousOn).intervalIntegrable_of_Icc hab.le
  -- Step 1: `∫ g² e^{−V} = −∫ W h′`, by parts with the kinks allowed
  have hstep1 : ∫ x in a..b, g x ^ 2 * exp (-V x) = -∫ x in a..b, W x * h' x := by
    have hftc := integral_eq_of_hasDerivAt_off_countable_of_le (fun x => W x * g x)
      (fun x => g x ^ 2 * exp (-V x) + W x * h' x) hab.le hs (hWc.mul hgc)
      (fun x hx => ((hWd x hx.1).mul ((hd x hx).sub_const m)).congr_deriv (by ring))
      (hig2.add hiWh)
    simp only [hWa, hWb, zero_mul, sub_self] at hftc
    rw [intervalIntegral.integral_add hig2 hiWh] at hftc
    linarith
  -- Step 2: the Bochner identity, through `Ψ = V′ W² e^{V}`
  have hstep2 : ∫ x in a..b, V'' x * W x ^ 2 * exp (V x) ≤ ∫ x in a..b, g x ^ 2 * exp (-V x) := by
    have hiR : IntervalIntegrable
        (fun x => 2 * V' x * W x * g x + V' x ^ 2 * W x ^ 2 * exp (V x)) volume a b :=
      ((((continuousOn_const.mul hV'c.continuousOn).mul hWc).mul hgc).add
        (((hV'c.continuousOn.pow 2).mul (hWc.pow 2)).mul heVc.continuousOn))
        |>.intervalIntegrable_of_Icc hab.le
    have hftc := integral_eq_of_hasDerivAt_off_countable_of_le
      (fun x => V' x * W x ^ 2 * exp (V x))
      (fun x => V'' x * W x ^ 2 * exp (V x) +
        (2 * V' x * W x * g x + V' x ^ 2 * W x ^ 2 * exp (V x)))
      hab.le (s := ∅) countable_empty
      ((hV'c.continuousOn.mul (hWc.pow 2)).mul heVc.continuousOn)
      (fun x hx => (((hV' x).mul ((hWd x hx.1).fun_pow 2)).mul (hV x).exp).congr_deriv (by
        have := hωe x
        simp only [Pi.mul_apply]
        linear_combination (2 * V' x * W x * g x) * this))
      (hiQ.add hiR)
    simp only [hWa, hWb, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero,
      zero_mul, sub_self] at hftc
    rw [intervalIntegral.integral_add hiQ hiR] at hftc
    -- `∫ g² e^{−V} − ∫ V″ W² e^V = ∫ (g + V′ W e^V)² e^{−V} ≥ 0`
    have hsq : ∀ x, (g x + V' x * W x * exp (V x)) ^ 2 * exp (-V x) =
        g x ^ 2 * exp (-V x) + (2 * V' x * W x * g x + V' x ^ 2 * W x ^ 2 * exp (V x)) := by
      intro x
      have := hωe x
      linear_combination (2 * V' x * W x * g x + V' x ^ 2 * W x ^ 2 * exp (V x)) * this
    have hnn : 0 ≤ ∫ x in a..b, (g x + V' x * W x * exp (V x)) ^ 2 * exp (-V x) :=
      intervalIntegral.integral_nonneg hab.le fun x _ =>
        mul_nonneg (sq_nonneg _) (Real.exp_pos _).le
    rw [intervalIntegral.integral_congr (fun x _ => hsq x),
      intervalIntegral.integral_add hig2 hiR] at hnn
    linarith
  -- Step 3: pointwise, `−W h′ ≤ ½ V″ W² e^V + ½ κ⁻¹ h′² e^{−V}` on the interval
  have hpt : ∀ x ∈ Icc a b, -(W x * h' x) ≤
      (1 / 2) * (V'' x * W x ^ 2 * exp (V x)) + (1 / 2) * κ⁻¹ * (h' x ^ 2 * exp (-V x)) := by
    intro x hx
    have hV'' : κ ≤ V'' x := hV''κ x hx
    have hA : 0 < V'' x * exp (V x) := mul_pos (hκ.trans_le hV'') (Real.exp_pos _)
    have hω := hωe x
    have key : 0 ≤ (V'' x * exp (V x) * W x + h' x) ^ 2 / (V'' x * exp (V x)) :=
      div_nonneg (sq_nonneg _) hA.le
    have hexp : (V'' x * exp (V x) * W x + h' x) ^ 2 / (V'' x * exp (V x)) =
        V'' x * W x ^ 2 * exp (V x) + 2 * (W x * h' x) + h' x ^ 2 * exp (-V x) / V'' x := by
      have hE : exp (V x) ≠ 0 := (Real.exp_pos _).ne'
      have hV0 : V'' x ≠ 0 := (hκ.trans_le hV'').ne'
      rw [Real.exp_neg]
      field_simp
      ring
    have hκ' : h' x ^ 2 * exp (-V x) / V'' x ≤ κ⁻¹ * (h' x ^ 2 * exp (-V x)) := by
      rw [div_eq_mul_inv, mul_comm]
      exact mul_le_mul_of_nonneg_right (inv_anti₀ hκ hV'')
        (mul_nonneg (sq_nonneg _) (Real.exp_pos _).le)
    nlinarith [key, hexp, hκ']
  have hstep3 : -∫ x in a..b, W x * h' x ≤
      (1 / 2) * (∫ x in a..b, V'' x * W x ^ 2 * exp (V x)) +
        (1 / 2) * κ⁻¹ * ∫ x in a..b, h' x ^ 2 * exp (-V x) := by
    have hmono : ∫ x in a..b, -(W x * h' x) ≤ ∫ x in a..b,
        ((1 / 2) * (V'' x * W x ^ 2 * exp (V x)) + (1 / 2) * κ⁻¹ * (h' x ^ 2 * exp (-V x))) :=
      intervalIntegral.integral_mono_on hab.le hiWh.neg
        ((hiQ.const_mul _).add (hih2.const_mul _)) hpt
    have hA : ∫ x in a..b, (1 / 2) * (V'' x * W x ^ 2 * exp (V x)) =
        (1 / 2) * ∫ x in a..b, V'' x * W x ^ 2 * exp (V x) :=
      intervalIntegral.integral_const_mul _ _
    have hB : ∫ x in a..b, (1 / 2) * κ⁻¹ * (h' x ^ 2 * exp (-V x)) =
        (1 / 2) * κ⁻¹ * ∫ x in a..b, h' x ^ 2 * exp (-V x) :=
      intervalIntegral.integral_const_mul _ _
    have hsum := intervalIntegral.integral_add (hiQ.const_mul (1 / 2))
      (hih2.const_mul ((1 / 2) * κ⁻¹))
    rw [hsum, hA, hB, intervalIntegral.integral_neg] at hmono
    exact hmono
  have hg2 : ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) = ∫ x in a..b, g x ^ 2 * exp (-V x) := rfl
  rw [hg2]
  set G := ∫ x in a..b, g x ^ 2 * exp (-V x)
  set Q := ∫ x in a..b, V'' x * W x ^ 2 * exp (V x)
  set H := ∫ x in a..b, h' x ^ 2 * exp (-V x)
  set I := ∫ x in a..b, W x * h' x
  set T := κ⁻¹ * H with hT
  have e : (1 / 2) * κ⁻¹ * H = (1 / 2) * T := by rw [hT]; ring
  rw [e] at hstep3
  linarith [hstep1, hstep2, hstep3]

/-- **AS-2(e), the truncated normal.** For the normal `N(μ, σ²)` restricted to `[a, b]`
(`V = (x − μ)²/(2σ²)`, so `V″ = σ⁻²`), the weighted variance is at most `σ²` times the weighted
mean squared derivative. -/
theorem poincare_truncNormal (hab : a < b) {μ σ : ℝ} (hσ : 0 < σ) {h h' : ℝ → ℝ} {s : Set ℝ}
    (hs : s.Countable) (hc : ContinuousOn h (Icc a b))
    (hd : ∀ x ∈ Ioo a b \ s, HasDerivAt h (h' x) x) (hm : Measurable h') {K : ℝ}
    (hK : ∀ x ∈ Icc a b, |h' x| ≤ K) (m : ℝ)
    (hmean : ∫ x in a..b, (h x - m) * exp (-((x - μ) ^ 2 / (2 * σ ^ 2))) = 0) :
    ∫ x in a..b, (h x - m) ^ 2 * exp (-((x - μ) ^ 2 / (2 * σ ^ 2))) ≤
      σ ^ 2 * ∫ x in a..b, h' x ^ 2 * exp (-((x - μ) ^ 2 / (2 * σ ^ 2))) := by
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
  have h := poincare_weighted hab (κ := 1 / σ ^ 2) (by positivity) hV hV' continuous_const
    (fun _ _ => le_rfl) hs hc hd hm hK m hmean
  rwa [one_div, inv_inv] at h

end Weighted

end ActiveSubspace
