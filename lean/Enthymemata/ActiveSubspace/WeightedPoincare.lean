import Enthymemata.ActiveSubspace.ContinuousPoincare
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.MeasureTheory.Integral.DivergenceTheorem
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Probability.ConditionalProbability
import Mathlib.Probability.Distributions.Gaussian.Real

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
* `hasPoincareOn_of_density`: the same for a law, `HasPoincareOn ν a b κ⁻¹` when `ν ∝ e^{−V}` on
  `[a, b]`.
* `truncNormal μ v a b`, `gaussianReal μ v` conditioned on `[a, b]`, with
  `hasPoincareOn_truncNormal`: its constant is `v`, so its coordinates join uniform ones in
  `matrix_poincare_prod_box` (`ae_mem_truncNormal`).
* `matrix_poincare_box_mixed`: the matrix check on E8's laws, each input uniform or a truncated
  normal, with `Dᵢ` its constant.
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

section Laws

open ProbabilityTheory
open scoped NNReal

variable {a b : ℝ}

/-- **AS-2(e) as a law.** A probability law `ν` with `∫ φ dν = C ∫ₐᵇ φ e^{−V}` for every `φ`
(`ν ∝ e^{−V}` on `[a, b]`), with `V″ ≥ κ > 0` on `[a, b]`, has the one-variable constant `κ⁻¹`. -/
theorem hasPoincareOn_of_density (hab : a < b) {V V' V'' : ℝ → ℝ} {κ : ℝ} (hκ : 0 < κ)
    (hV : ∀ x, HasDerivAt V (V' x) x) (hV' : ∀ x, HasDerivAt V' (V'' x) x)
    (hV''c : Continuous V'') (hV''κ : ∀ x ∈ Icc a b, κ ≤ V'' x)
    (ν : Measure ℝ) [IsProbabilityMeasure ν] {C : ℝ}
    (hν : ∀ φ : ℝ → ℝ, ∫ x, φ x ∂ν = C * ∫ x in a..b, φ x * exp (-V x)) :
    HasPoincareOn ν a b κ⁻¹ := by
  intro h h' s K hs hc hd hm hK
  have hVc : Continuous V := continuous_iff_continuousAt.2 fun x => (hV x).continuousAt
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
  have hW := poincare_weighted hab hκ hV hV' hV''c hV''κ hs hc hd hm (fun x _ => hK x) m hmean
  have hvar : varM ν h = C * ∫ x in a..b, (h x - m) ^ 2 * exp (-V x) :=
    hν fun x => (h x - m) ^ 2
  have hen : ∫ x, h' x ^ 2 ∂ν = C * ∫ x in a..b, h' x ^ 2 * exp (-V x) :=
    hν fun x => h' x ^ 2
  rw [hvar, hen]
  calc C * ∫ x in a..b, (h x - m) ^ 2 * exp (-V x)
      ≤ C * (κ⁻¹ * ∫ x in a..b, h' x ^ 2 * exp (-V x)) := mul_le_mul_of_nonneg_left hW hC0
    _ = κ⁻¹ * (C * ∫ x in a..b, h' x ^ 2 * exp (-V x)) := by ring

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

/-- **AS-2(e), the truncated normal as a law.** `N(μ, v)` truncated to `[a, b]` has the
one-variable Poincaré constant `v` on `[a, b]`: Bakry–Émery with `V″ = 1/v`. -/
theorem hasPoincareOn_truncNormal (hab : a < b) (μ : ℝ) {v : ℝ≥0} (hv : v ≠ 0) :
    HasPoincareOn (truncNormal μ v a b) a b v := by
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
  have h := hasPoincareOn_of_density hab (κ := 1 / v) (one_div_pos.2 hv0) hV hV'
    continuous_const (fun _ _ => le_rfl) (truncNormal μ v a b) (integral_truncNormal hab μ hv)
  rwa [one_div, inv_inv] at h

/-- **AS-2(d) and (e), on E8's laws.** Each input is uniform on `[aᵢ, bᵢ]`, or, for `i ∈ N`, the
normal `N(μᵢ, vᵢ)` truncated to it. Then `Cov(f) ≼ E[J D Jᵀ]` with `Dᵢ = ((bᵢ − aᵢ)/π)²` for a
uniform input and `Dᵢ = vᵢ` for a truncated normal, every hypothesis asked only on the box. -/
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
        (fun i => if i ∈ N then (v i : ℝ) else ((b i - a i) / π) ^ 2) -
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
