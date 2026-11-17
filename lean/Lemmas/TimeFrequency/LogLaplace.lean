import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.MeasureTheory.Integral.Prod

/-! # Logarithms from Laplace transforms on the positive half-line -/

open MeasureTheory Real Set Filter

namespace TimeFrequency

noncomputable def logKernel (x t : ℝ) := (exp (-t) - exp (-x * t)) / t

/-- Frullani's formula with absolute integrability, first for ordered endpoints. -/
theorem frullani_of_le {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) :
    IntegrableOn (fun t ↦ (exp (-a * t) - exp (-b * t)) / t) (Ioi 0) ∧
      (∫ t in Ioi 0, (exp (-a * t) - exp (-b * t)) / t) = log (b / a) := by
  let μ := volume.restrict (Ioc a b)
  let ν := volume.restrict (Ioi (0 : ℝ))
  let f := fun p : ℝ × ℝ ↦ exp (-p.1 * p.2)
  have hinner {u : ℝ} (hu : 0 < u) : ∫ t, f (u,t) ∂ν = u⁻¹ := by
    simpa [f, ν] using integral_exp_mul_Ioi (neg_lt_zero.mpr hu) 0
  have hf : Integrable f (μ.prod ν) := by
    apply (integrable_prod_iff (by dsimp [f]; fun_prop)).mpr
    constructor
    · filter_upwards [ae_restrict_mem measurableSet_Ioc] with u hu
      exact integrableOn_exp_mul_Ioi (neg_lt_zero.mpr (ha.trans hu.1)) 0
    · have huinv : IntegrableOn (fun u : ℝ ↦ u⁻¹) (Ioc a b) :=
        (continuousOn_inv₀.mono fun u hu ↦ ne_of_gt (ha.trans_le hu.1)).integrableOn_Icc.mono_set
          Ioc_subset_Icc_self
      apply huinv.congr
      filter_upwards [ae_restrict_mem measurableSet_Ioc] with u hu
      simpa [f, Real.norm_eq_abs, abs_of_pos (exp_pos _)] using (hinner (ha.trans hu.1)).symm
  have ht : ∀ᵐ t ∂ν, (∫ u, f (u,t) ∂μ) = (exp (-a*t) - exp (-b*t)) / t := by
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    change (∫ u in Ioc a b, exp (-u * t)) = _
    rw [← intervalIntegral.integral_of_le hab]
    simp_rw [neg_mul, ← mul_neg]
    rw [intervalIntegral.integral_comp_mul_right exp (neg_ne_zero.mpr (ne_of_gt ht)), integral_exp]
    simp only [smul_eq_mul]
    ring
  refine ⟨hf.integral_prod_right.congr ht, ?_⟩
  rw [← integral_congr_ae ht, ← integral_prod_symm _ hf, integral_prod _ hf]
  calc
    _ = ∫ u in Ioc a b, u⁻¹ := by
      apply integral_congr_ae
      filter_upwards [ae_restrict_mem measurableSet_Ioc] with u hu
      exact hinner (ha.trans hu.1)
    _ = _ := by rw [← intervalIntegral.integral_of_le hab, integral_inv_of_pos ha (ha.trans_le hab)]

theorem logKernel_spec {x : ℝ} (hx : 0 < x) :
    IntegrableOn (logKernel x) (Ioi 0) ∧ (∫ t in Ioi 0, logKernel x t) = log x := by
  rcases le_total 1 x with h | h
  · unfold logKernel
    simpa using frullani_of_le zero_lt_one h
  · obtain ⟨hi, he⟩ := frullani_of_le hx h
    have hf : logKernel x = fun t ↦ -((exp (-x*t) - exp (-(1:ℝ)*t)) / t) := by
      funext t
      simp [logKernel]
      ring
    rw [hf, integral_neg, he]
    exact ⟨hi.neg, by simp⟩

theorem integral_norm_logKernel {x : ℝ} (hx : 0 < x) :
    (∫ t in Ioi 0, ‖logKernel x t‖) = |log x| := by
  obtain ⟨_, he⟩ := logKernel_spec hx
  rcases le_total 1 x with h | h
  · rw [abs_of_nonneg (log_nonneg h), ← he]
    apply integral_congr_ae
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    apply Real.norm_of_nonneg
    exact div_nonneg (sub_nonneg.mpr (exp_le_exp.mpr (by nlinarith [show 0 < t from ht])))
      (le_of_lt ht)
  · rw [abs_of_nonpos (log_nonpos hx.le h), ← he, ← integral_neg]
    apply integral_congr_ae
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    apply Real.norm_of_nonpos
    exact div_nonpos_of_nonpos_of_nonneg (sub_nonpos.mpr
      (exp_le_exp.mpr (by nlinarith [show 0 < t from ht]))) (le_of_lt ht)

#print axioms frullani_of_le
#print axioms logKernel_spec
#print axioms integral_norm_logKernel

noncomputable abbrev positiveLine : Measure ℝ := volume.restrict (Ioi 0)

/-- Fubini for the log-power product. Absolute integrability follows from the actual
log-product moment because each one-dimensional kernel has constant sign in its scale. -/
theorem log_product_laplace {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω} [IsFiniteMeasure μ]
    {X Y : Ω → ℝ} (hX : Measurable X) (hY : Measurable Y) (hpX : ∀ᵐ ω ∂μ, 0 < X ω)
    (hpY : ∀ᵐ ω ∂μ, 0 < Y ω) (hi : Integrable (fun ω ↦ log (X ω) * log (Y ω)) μ) :
    Integrable (fun p : ℝ × ℝ ↦ ∫ ω, logKernel (X ω) p.1 * logKernel (Y ω) p.2 ∂μ)
        (positiveLine.prod positiveLine) ∧
      (∫ p, (∫ ω, logKernel (X ω) p.1 * logKernel (Y ω) p.2 ∂μ)
        ∂positiveLine.prod positiveLine) = ∫ ω, log (X ω) * log (Y ω) ∂μ := by
  let F := fun z : Ω × (ℝ × ℝ) ↦ logKernel (X z.1) z.2.1 * logKernel (Y z.1) z.2.2
  have hF : Integrable F (μ.prod (positiveLine.prod positiveLine)) := by
    apply (integrable_prod_iff (by dsimp [F, logKernel]; fun_prop)).mpr
    constructor
    · filter_upwards [hpX, hpY] with ω hx hy
      exact (logKernel_spec hx).1.mul_prod (logKernel_spec hy).1
    · apply hi.norm.congr
      filter_upwards [hpX, hpY] with ω hx hy
      simp only [F, norm_mul]
      rw [integral_prod_mul (fun t ↦ ‖logKernel (X ω) t‖) (fun t ↦ ‖logKernel (Y ω) t‖),
        integral_norm_logKernel hx, integral_norm_logKernel hy]
      rfl
  refine ⟨hF.integral_prod_right, (integral_prod_symm F hF).symm.trans ?_⟩
  rw [integral_prod F hF]
  apply integral_congr_ae
  filter_upwards [hpX, hpY] with ω hx hy
  simp only [F]
  rw [integral_prod_mul (logKernel (X ω)) (logKernel (Y ω)), (logKernel_spec hx).2,
    (logKernel_spec hy).2]

#print axioms log_product_laplace

end TimeFrequency
