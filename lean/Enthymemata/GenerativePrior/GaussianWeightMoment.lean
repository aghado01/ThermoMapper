import Enthymemata.GenerativePrior.WeightCondition
import Mathlib.MeasureTheory.Integral.Gamma
import Mathlib.Probability.Distributions.Gaussian.Real

/-! # The standard planar Gaussian realizes the weight-condition matrix (S18) -/

open Real Matrix Set MeasureTheory MeasureTheory.Measure ProbabilityTheory

namespace GenerativePrior

/-- A row's contribution to the second moment with both input gates open. -/
noncomputable def gatedEntry (θ : ℝ) (i j : Fin 2) (p : ℝ × ℝ) : ℝ :=
  if 0 < p.1 ∧ 0 < p.1 * cos θ + p.2 * sin θ then
    ![p.1, p.2] i * ![p.1, p.2] j else 0

theorem integrable_gatedEntry (θ : ℝ) (i j : Fin 2) :
    Integrable (gatedEntry θ i j) ((gaussianReal 0 1).prod (gaussianReal 0 1)) := by
  have hcoord (k : Fin 2) : MemLp (fun p : ℝ × ℝ ↦ ![p.1, p.2] k) 2
      ((gaussianReal 0 1).prod (gaussianReal 0 1)) := by
    have h := memLp_id_gaussianReal' (μ := 0) (v := 1) 2 (by norm_num)
    fin_cases k
    · exact h.comp_measurePreserving measurePreserving_fst
    · exact h.comp_measurePreserving measurePreserving_snd
  have hm : MeasurableSet {p : ℝ × ℝ | 0 < p.1 ∧ 0 < p.1 * cos θ + p.2 * sin θ} :=
    (measurableSet_lt measurable_const measurable_fst).inter
      (measurableSet_lt measurable_const
        (show Measurable (fun p : ℝ × ℝ ↦ p.1 * cos θ + p.2 * sin θ) by fun_prop))
  apply (((hcoord i).integrable_mul (hcoord j)).indicator hm).congr
  apply ae_of_all
  intro p
  by_cases h : 0 < p.1 ∧ 0 < p.1 * cos θ + p.2 * sin θ <;> simp [gatedEntry, indicator, h]

theorem standardPlane_integral (f : ℝ × ℝ → ℝ) :
    ∫ p, f p ∂(gaussianReal 0 1).prod (gaussianReal 0 1) =
      ∫ p, (gaussianPDFReal 0 1 p.1 * gaussianPDFReal 0 1 p.2) * f p := by
  have hm := measurable_gaussianPDF 0 1
  have hm₂ : Measurable (fun p : ℝ × ℝ ↦ gaussianPDF 0 1 p.1 * gaussianPDF 0 1 p.2) :=
    (hm.comp measurable_fst).mul (hm.comp measurable_snd)
  rw [gaussianReal_of_var_ne_zero _ one_ne_zero, prod_withDensity hm hm,
    integral_withDensity_eq_integral_toReal_smul hm₂ (ae_of_all _ fun _ ↦
      ENNReal.mul_lt_top gaussianPDF_lt_top gaussianPDF_lt_top)]
  simp only [ENNReal.toReal_mul, toReal_gaussianPDF, smul_eq_mul, volume_eq_prod]

theorem standardPlane_polar_density (r φ : ℝ) :
    gaussianPDFReal 0 1 (r * cos φ) * gaussianPDFReal 0 1 (r * sin φ) =
      exp (-(1 / 2) * r ^ 2) / (2 * π) := by
  simp only [gaussianPDFReal, NNReal.coe_one, mul_one, sub_zero]
  rw [mul_mul_mul_comm, ← exp_add, ← mul_inv, mul_self_sqrt (by positivity)]
  have h : -(r * cos φ) ^ 2 / 2 + -(r * sin φ) ^ 2 / 2 = -(1 / 2) * r ^ 2 := by
    nlinarith only [congrArg (fun t : ℝ ↦ r ^ 2 * t) (cos_sq_add_sin_sq φ)]
  rw [h]
  ring

theorem gaussian_radial_moment :
    (∫ r in Ioi (0 : ℝ), r ^ 3 * exp (-(1 / 2) * r ^ 2) / (2 * π)) = 1 / π := by
  rw [integral_div]
  have h := integral_rpow_mul_exp_neg_mul_rpow (p := 2) (q := 3) (b := 1 / 2)
    (by norm_num) (by norm_num) (by norm_num)
  norm_num [Gamma_add_one, Gamma_one] at h
  simp_rw [neg_mul, h]
  field_simp

theorem polar_gatedEntry {r φ θ : ℝ} (hr : 0 < r) (hφ : φ ∈ Ioo (-π) π)
    (hθ : θ ∈ Icc 0 π) (i j : Fin 2) :
    gatedEntry θ i j (r * cos φ, r * sin φ) =
      r ^ 2 * (Ioo (θ - π / 2) (π / 2)).indicator
        (fun t ↦ planeDirection t i * planeDirection t j) φ := by
  have hgate : 0 < r * cos φ ∧ 0 < r * cos φ * cos θ + r * sin φ * sin θ ↔
      φ ∈ Ioo (θ - π / 2) (π / 2) := by
    rw [mul_assoc, mul_assoc r (sin φ), ← mul_add]
    simpa only [mul_pos_iff_of_pos_left hr] using both_gates_iff hφ hθ
  have hcoord (k : Fin 2) : ![r * cos φ, r * sin φ] k = r * planeDirection φ k := by
    fin_cases k <;> rfl
  by_cases ha : φ ∈ Ioo (θ - π / 2) (π / 2)
  · simp [gatedEntry, hgate.mpr ha, indicator_of_mem ha, hcoord]
    ring
  · simp [gatedEntry, mt hgate.mp ha, indicator_of_notMem ha]

/-- The expectation under independent standard Gaussian coordinates, including both endpoints. -/
theorem gaussian_gatedMoment_eq_weightQ {θ : ℝ} (hθ : θ ∈ Icc 0 π) :
    (fun i j ↦ ∫ p, gatedEntry θ i j p ∂(gaussianReal 0 1).prod (gaussianReal 0 1)) =
      weightQ θ := by
  rw [← angularMoment_eq_weightQ]
  ext i j
  let arc := Ioo (θ - π / 2) (π / 2)
  let kernel := fun t ↦ planeDirection t i * planeDirection t j
  have hsub : arc ⊆ Ioo (-π) π := fun t ht ↦
    ⟨by linarith [ht.1, hθ.1, pi_pos], by linarith [ht.2, pi_pos]⟩
  rw [standardPlane_integral, ← integral_comp_polarCoord_symm]
  calc
    _ = ∫ p in Ioi (0 : ℝ) ×ˢ Ioo (-π) π,
        (p.1 ^ 3 * exp (-(1 / 2) * p.1 ^ 2) / (2 * π)) * arc.indicator kernel p.2 := by
      refine setIntegral_congr_fun (measurableSet_Ioi.prod measurableSet_Ioo) ?_
      rintro ⟨r, φ⟩ ⟨hr, hφ⟩
      change r * (gaussianPDFReal 0 1 (r * cos φ) * gaussianPDFReal 0 1 (r * sin φ) *
        gatedEntry θ i j (r * cos φ, r * sin φ)) = _
      rw [standardPlane_polar_density, polar_gatedEntry hr hφ hθ]
      dsimp only [Prod.fst, Prod.snd, arc, kernel]
      ring
    _ = angularMoment θ i j := by
      rw [volume_eq_prod, setIntegral_prod_mul
        (fun r : ℝ ↦ r ^ 3 * exp (-(1 / 2) * r ^ 2) / (2 * π)) (arc.indicator kernel),
        gaussian_radial_moment,
        setIntegral_indicator measurableSet_Ioo, inter_eq_right.mpr hsub]
      unfold angularMoment
      rw [intervalIntegral.integral_of_le (by linarith [hθ.2])]
      exact congrArg ((1 / π) * ·) (setIntegral_congr_set Ioo_ae_eq_Ioc)

#print axioms integrable_gatedEntry
#print axioms standardPlane_integral
#print axioms standardPlane_polar_density
#print axioms gaussian_radial_moment
#print axioms polar_gatedEntry
#print axioms gaussian_gatedMoment_eq_weightQ

end GenerativePrior
