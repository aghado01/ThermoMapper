import Lemmas.TimeFrequency.GaussianFourthMoment
import Lemmas.TimeFrequency.QuadraticFormVariance

/-! # Isserlis for the coordinates of a Gaussian vector, and the variance of its quadratic form -/

open MeasureTheory ProbabilityTheory Matrix
open scoped ProbabilityTheory

namespace TimeFrequency

variable {Ω n : Type*} [MeasurableSpace Ω] [Fintype n] {μ : Measure Ω} {x : Ω → n → ℝ}

noncomputable def secondMoments (μ : Measure Ω) (x : Ω → n → ℝ) : Matrix n n ℝ :=
  fun i j ↦ ∫ ω, x ω i*x ω j ∂μ

theorem gaussian_coordinate_pairing (hx : HasGaussianLaw x μ) (h0 : ∫ ω, x ω ∂μ = 0) :
    Pairing (secondMoments μ x) (fun p q r s ↦ ∫ ω, x ω p*x ω q*(x ω r*x ω s) ∂μ) := by
  let ν := μ.map x
  letI : IsGaussian ν := hx.isGaussian_map
  let A (i : n) : StrongDual ℝ (n → ℝ) := ContinuousLinearMap.proj i
  have hm : ∫ v, v ∂ν = 0 := (integral_map hx.aemeasurable (by fun_prop)).trans h0
  have hK (i j) : covarianceBilinDual ν (A i) (A j) = secondMoments μ x i j := by
    rw [covarianceBilinDual_apply' (IsGaussian.hasGaussianLaw_id.memLp_two)]
    change ∫ v, A i (v - ∫ z, z ∂ν) * A j (v - ∫ z, z ∂ν) ∂ν = _
    simp only [hm, sub_zero]
    exact integral_map hx.aemeasurable (by fun_prop : AEStronglyMeasurable
      (fun v ↦ A i v * A j v) ν)
  intro p q r s
  have h := gaussian_four_projections ν hm (A p) (A q) (A r) (A s)
  simp only [hK] at h
  refine Eq.trans ?_ h
  symm
  exact integral_map hx.aemeasurable (by fun_prop : AEStronglyMeasurable
    (fun v ↦ A p v * A q v * (A r v * A s v)) ν)

/-- S15, allowing a singular Gaussian covariance and any symmetric real matrix. -/
theorem gaussian_quadratic_variance (hx : HasGaussianLaw x μ) (h0 : ∫ ω, x ω ∂μ = 0)
    {H : Matrix n n ℝ} (hH : H.IsSymm) :
    (∫ ω, (x ω ⬝ᵥ (H *ᵥ x ω))^2 ∂μ) - (∫ ω, x ω ⬝ᵥ (H *ᵥ x ω) ∂μ)^2 =
      2*((H * secondMoments μ x)*(H * secondMoments μ x)).trace := by
  have h2 (i j : n) : MemLp (fun ω ↦ x ω i*x ω j) 2 μ := by
    exact gaussian_product_memLp_two (hx.eval i) (hx.eval j)
  apply variance_quadForm μ x hH (secondMoments μ x)
    (fun i j ↦ ⟨(hx.eval i).memLp_two.integrable_mul (hx.eval j).memLp_two, rfl⟩)
    (fun p q r s ↦ ∫ ω, x ω p*x ω q*(x ω r*x ω s) ∂μ)
    (fun p q r s ↦ ⟨(h2 p q).integrable_mul (h2 r s), rfl⟩)
  exact gaussian_coordinate_pairing hx h0

theorem variance_gaussian_quadForm (hx : HasGaussianLaw x μ) (h0 : ∫ ω, x ω ∂μ = 0)
    {H : Matrix n n ℝ} (hH : H.IsSymm) :
    Var[fun ω ↦ x ω ⬝ᵥ (H *ᵥ x ω); μ] =
      2*((H * secondMoments μ x)*(H * secondMoments μ x)).trace := by
  letI := hx.isProbabilityMeasure
  have hM : MemLp (fun ω ↦ x ω ⬝ᵥ (H *ᵥ x ω)) 2 μ := by
    simp_rw [dotProduct_mulVec_eq_sum]
    exact memLp_finsetSum _ fun i _ ↦ memLp_finsetSum _ fun j _ ↦
      (gaussian_product_memLp_two (hx.eval i) (hx.eval j)).const_mul (H i j)
  rw [variance_eq_sub hM]
  exact gaussian_quadratic_variance hx h0 hH

/-- The standard-normal one-coordinate quadratic form has variance two. -/
theorem standardNormal_square_variance :
    (∫ z, (z^2)^2 ∂gaussianReal 0 1) - (∫ z, z^2 ∂gaussianReal 0 1)^2 = 2 := by
  have h2 := variance_eq_integral (μ := gaussianReal 0 1) measurable_id.aemeasurable
  change Var[id; gaussianReal 0 1] = ∫ z, (z - ∫ u, u ∂gaussianReal 0 1)^2 ∂gaussianReal 0 1 at h2
  simp only [variance_id_gaussianReal, integral_id_gaussianReal, sub_zero, NNReal.coe_one] at h2
  rw [← h2]
  simp only [← pow_mul]
  rw [gaussianReal_fourth]
  norm_num

#print axioms gaussian_coordinate_pairing
#print axioms gaussian_quadratic_variance
#print axioms variance_gaussian_quadForm
#print axioms standardNormal_square_variance

end TimeFrequency
