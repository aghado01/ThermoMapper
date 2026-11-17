import Lemmas.TimeFrequency.PeakOfExponentials
import Mathlib.MeasureTheory.Integral.Layercake
import Mathlib.Probability.Distributions.Exponential
import Mathlib.Probability.Independence.Basic

/-! # The maximum of a sample of independent unit exponentials -/

open MeasureTheory ProbabilityTheory Real Set

namespace TimeFrequency

variable {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω} [IsProbabilityMeasure μ]

noncomputable def sampleMaximum {n : ℕ} (X : Fin (n + 1) → Ω → ℝ) (ω : Ω) : ℝ :=
  Finset.univ.sup' Finset.univ_nonempty (fun i ↦ X i ω)

theorem sampleMaximum_measurable {n : ℕ} {X : Fin (n + 1) → Ω → ℝ}
    (hm : ∀ i, Measurable (X i)) : Measurable (sampleMaximum X) := by
  have he : sampleMaximum X = Finset.univ.sup' Finset.univ_nonempty X := by
    funext ω
    simp [sampleMaximum, Finset.sup'_apply]
  rw [he]
  exact Finset.measurable_sup' Finset.univ_nonempty fun i _ ↦ hm i

omit [IsProbabilityMeasure μ] in
theorem exponentialMaximum_cdf {n : ℕ} {X : Fin (n + 1) → Ω → ℝ}
    (hm : ∀ i, Measurable (X i)) (law : ∀ i, μ.map (X i) = expMeasure 1)
    (ind : iIndepFun X μ) {t : ℝ} (ht : 0 ≤ t) :
    μ.real {ω | sampleMaximum X ω ≤ t} = (1-exp (-t))^(n + 1) := by
  letI := isProbabilityMeasure_expMeasure (r := 1) zero_lt_one
  have he : {ω | sampleMaximum X ω ≤ t} = ⋂ i, X i ⁻¹' Iic t := by
    ext ω
    simp [sampleMaximum, Finset.sup'_le_iff]
  have hi (i) : μ.real (X i ⁻¹' Iic t) = 1-exp (-t) := by
    rw [Measure.real, ← Measure.map_apply (hm i) measurableSet_Iic, law,
      ← Measure.real, ← cdf_eq_real, cdf_expMeasure_eq zero_lt_one]
    simp [ht]
  rw [he, Measure.real, ind.meas_iInter (fun i ↦ ⟨Iic t, measurableSet_Iic, rfl⟩),
    ENNReal.toReal_prod]
  change (∏ i, μ.real (X i ⁻¹' Iic t)) = _
  simp only [hi, Finset.prod_const, Finset.card_univ, Fintype.card_fin]

/-- S12 for every nonempty sample, including the single-bin case. -/
theorem expected_exponentialMaximum {n : ℕ} {X : Fin (n + 1) → Ω → ℝ}
    (hm : ∀ i, Measurable (X i)) (law : ∀ i, μ.map (X i) = expMeasure 1)
    (ind : iIndepFun X μ) : ∫ ω, sampleMaximum X ω ∂μ = (harmonic (n + 1) : ℝ) := by
  have hmM := sampleMaximum_measurable hm
  have hc {t : ℝ} (ht : 0 ≤ t) := exponentialMaximum_cdf hm law ind ht
  have hn : ∀ᵐ ω ∂μ, 0 ≤ sampleMaximum X ω := by
    have hz : μ {ω | sampleMaximum X ω ≤ 0} = 0 := by
      apply (measureReal_eq_zero_iff (measure_ne_top _ _)).mp
      simpa using hc (t := 0) le_rfl
    have hpos : ∀ᵐ ω ∂μ, 0 < sampleMaximum X ω := by simpa [ae_iff, not_lt] using hz
    exact hpos.mono fun _ h ↦ h.le
  have ht (t : ℝ) (h : 0 ≤ t) :
      μ.real {ω | t < sampleMaximum X ω} = 1-(1-exp (-t))^(n + 1) := by
    have he : {ω | t < sampleMaximum X ω} = {ω | sampleMaximum X ω ≤ t}ᶜ := by
      ext; simp
    rw [he, measureReal_compl (measurableSet_le hmM measurable_const), probReal_univ, hc h]
  rw [integral_eq_lintegral_of_nonneg_ae hn hmM.aestronglyMeasurable,
    lintegral_eq_lintegral_meas_lt μ hn hmM.aemeasurable,
    ← integral_one_sub_pow_one_sub_exp (n + 1)]
  rw [integral_eq_lintegral_of_nonneg_ae (by
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht'
    rw [← ht t (le_of_lt ht')]
    positivity) (by apply Measurable.aestronglyMeasurable; fun_prop)]
  congr 1
  apply lintegral_congr_ae
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht'
  rw [← ht t (le_of_lt ht'), Measure.real, ENNReal.ofReal_toReal (measure_ne_top _ _)]

theorem exponentialMaximum_integrable {n : ℕ} {X : Fin (n + 1) → Ω → ℝ}
    (hm : ∀ i, Measurable (X i)) (law : ∀ i, μ.map (X i) = expMeasure 1)
    (ind : iIndepFun X μ) : Integrable (sampleMaximum X) μ := by
  by_contra hn
  have h := expected_exponentialMaximum hm law ind
  rw [integral_undef hn] at h
  have hpos : 0 < harmonic (n + 1) := by
    rw [harmonic_succ]
    have : 0 ≤ harmonic n := Finset.sum_nonneg fun i _ ↦ by positivity
    positivity
  exact (ne_of_gt (show 0 < (harmonic (n + 1) : ℝ) by exact_mod_cast hpos)) h.symm

/-- The two independent coordinates of the product exponential law have expected maximum 3/2. -/
theorem exponentialMaximum_two :
    ∫ z : Fin 2 → ℝ, sampleMaximum (fun i z ↦ z i) z ∂Measure.pi (fun _ ↦ expMeasure 1) =
      3/2 := by
  letI := isProbabilityMeasure_expMeasure (r := 1) zero_lt_one
  rw [expected_exponentialMaximum (n := 1) (fun i ↦ measurable_pi_apply i)
    (fun i ↦ (measurePreserving_eval (fun _ : Fin 2 ↦ expMeasure 1) i).map_eq)
    (iIndepFun_pi (X := fun _ : Fin 2 ↦ id) (fun _ ↦ aemeasurable_id))]
  norm_num [harmonic]

#print axioms sampleMaximum_measurable
#print axioms exponentialMaximum_cdf
#print axioms expected_exponentialMaximum
#print axioms exponentialMaximum_integrable
#print axioms exponentialMaximum_two

end TimeFrequency
