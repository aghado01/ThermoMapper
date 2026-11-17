import Mathlib.Probability.Distributions.Gaussian.HasGaussianLaw.Basic
import Mathlib.Probability.Moments.CovarianceBilinDual

/-! # The Gaussian fourth moment and its polarization -/

open MeasureTheory ProbabilityTheory Real
open scoped NNReal ENNReal ProbabilityTheory

namespace TimeFrequency

/-- The fourth moment also covers the point mass at variance zero. -/
theorem gaussianReal_fourth (v : ℝ≥0) : ∫ x, x^4 ∂gaussianReal 0 v = 3*(v:ℝ)^2 := by
  let a := fun t : ℝ ↦ (v:ℝ)*t
  let f := fun t : ℝ ↦ exp ((v:ℝ)*t^2/2)
  have da (t : ℝ) : HasDerivAt a (v:ℝ) t := by
    simpa [a] using (hasDerivAt_id t).const_mul (v:ℝ)
  have df (t : ℝ) : HasDerivAt f (a t * f t) t := by
    convert ((((hasDerivAt_id t).pow 2).const_mul (v:ℝ)).div_const 2).exp using 1
    all_goals simp [f,a]
    all_goals ring
  have d2 (t : ℝ) : deriv (fun t ↦ a t*f t) t = ((v:ℝ)+(a t)^2)*f t := by
    refine ((da t).mul (df t)).deriv.trans ?_
    ring
  have d3 (t : ℝ) : deriv (fun t ↦ ((v:ℝ)+(a t)^2)*f t) t =
      (3*(v:ℝ)*a t+(a t)^3)*f t := by
    refine ((((da t).pow 2).const_add (v:ℝ)).mul (df t)).deriv.trans ?_
    simp only [Pi.pow_apply]
    ring
  have d4 := (((((da 0).const_mul (3*(v:ℝ))).add ((da 0).pow 3))).mul (df 0)).deriv
  change (∫ x, ((fun t : ℝ ↦ t)^4) x ∂gaussianReal 0 v) = _
  rw [← iteratedDeriv_mgf_zero (by simp), mgf_fun_id_gaussianReal]
  simp only [zero_mul, zero_add]
  change iteratedDeriv 4 f 0 = _
  rw [iteratedDeriv_succ, iteratedDeriv_succ, iteratedDeriv_succ, iteratedDeriv_one,
    show deriv f = _ from funext fun t ↦ (df t).deriv, funext d2, funext d3]
  refine d4.trans ?_
  simp [a, f]
  ring

theorem gaussian_fourth_eq_variance {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω}
    {Y : Ω → ℝ} (hY : HasGaussianLaw Y μ) (h0 : ∫ ω, Y ω ∂μ = 0) :
    ∫ ω, (Y ω)^4 ∂μ = 3*Var[Y; μ]^2 := by
  have hmap := hY.map_eq_gaussianReal
  rw [h0] at hmap
  calc ∫ ω, (Y ω)^4 ∂μ = ∫ y, y^4 ∂μ.map Y :=
        (integral_map hY.aemeasurable (by fun_prop : AEStronglyMeasurable (fun y : ℝ ↦ y^4)
          (μ.map Y))).symm
    _ = _ := by rw [hmap, gaussianReal_fourth, Real.coe_toNNReal _ (variance_nonneg _ _)]

section Polarization

theorem gaussian_product_memLp_two {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    {Y Z : Ω → ℝ} (hY : HasGaussianLaw Y P) (hZ : HasGaussianLaw Z P) :
    MemLp (fun ω ↦ Y ω * Z ω) 2 P := by
  letI : ENNReal.HolderTriple 4 4 2 := ⟨by
    apply (ENNReal.toReal_eq_toReal_iff' (by finiteness) (by finiteness)).mp
    rw [ENNReal.toReal_add (by finiteness) (by finiteness)]
    norm_num [ENNReal.toReal_inv]⟩
  exact (hZ.memLp (p := 4) (by norm_num)).mul' (r := 2) (hY.memLp (p := 4) (by norm_num))

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
  [SecondCountableTopology E] [MeasurableSpace E] [BorelSpace E]
  (μ : Measure E) [IsGaussian μ]

omit [CompleteSpace E] [SecondCountableTopology E] in
private theorem projection_law (A : StrongDual ℝ E) : HasGaussianLaw A μ :=
  (IsGaussian.hasGaussianLaw_id (μ := μ)).map_fun A

private theorem projection_mean (h0 : (∫ x, x ∂μ) = 0) (A : StrongDual ℝ E) :
    (∫ x, A x ∂μ) = 0 := by
  refine (A.integral_comp_comm (IsGaussian.hasGaussianLaw_id (μ := μ)).integrable).trans ?_
  simp [h0]

private theorem projection_fourth (h0 : (∫ x, x ∂μ) = 0) (A : StrongDual ℝ E) :
    ∫ x, (A x)^4 ∂μ = 3*(covarianceBilinDual μ A A)^2 := by
  rw [covarianceBilinDual_self_eq_variance (IsGaussian.hasGaussianLaw_id.memLp_two),
    gaussian_fourth_eq_variance (projection_law μ A) (projection_mean μ h0 A)]

omit [CompleteSpace E] [SecondCountableTopology E] in
private theorem integrable_projection_fourth (A : StrongDual ℝ E) :
    Integrable (fun x ↦ (A x)^4) μ := by
  have h := (projection_law μ A).memLp (p := 4) (by norm_num)
  simpa [Real.norm_eq_abs, Even.pow_abs (by decide : Even 4)] using
    (h.integrable_norm_pow (p := 4) (by decide))

omit [CompleteSpace E] [SecondCountableTopology E] in
private theorem integrable_square_product (A B : StrongDual ℝ E) :
    Integrable (fun x ↦ (A x)^2*(B x)^2) μ := by
  have h (L : StrongDual ℝ E) : MemLp (fun x ↦ (L x)^2) 2 μ := by
    apply (memLp_two_iff_integrable_sq (by fun_prop)).mpr
    simpa only [← pow_mul] using integrable_projection_fourth μ L
  exact (h A).integrable_mul (h B)

/-- Polarizing the one-dimensional fourth moment gives the mixed square moment. -/
theorem gaussian_square_product (h0 : (∫ x, x ∂μ) = 0) (A B : StrongDual ℝ E) :
    ∫ x, (A x)^2*(B x)^2 ∂μ =
      covarianceBilinDual μ A A * covarianceBilinDual μ B B + 2*(covarianceBilinDual μ A B)^2 := by
  have he (x : E) : 12*((A x)^2*(B x)^2) =
      ((A+B) x)^4 + ((A-B) x)^4 - 2*(A x)^4 - 2*(B x)^4 := by
    simp only [add_apply, sub_apply]
    ring
  have h := integral_congr_ae (μ := μ) (ae_of_all _ he)
  have hi := integrable_projection_fourth μ
  erw [integral_const_mul, integral_sub (((hi (A+B)).add (hi (A-B))).sub ((hi A).const_mul 2))
    ((hi B).const_mul 2), integral_sub ((hi (A+B)).add (hi (A-B))) ((hi A).const_mul 2),
    integral_add (hi (A+B)) (hi (A-B)), integral_const_mul, integral_const_mul] at h
  simp only [projection_fourth μ h0] at h
  simp only [map_add, map_sub, add_apply, sub_apply, covarianceBilinDual_comm B A] at h
  linarith only [h]

/-- Isserlis for any four continuous real linear observations of a centered Gaussian law. -/
theorem gaussian_four_projections (h0 : (∫ x, x ∂μ) = 0) (A B C D : StrongDual ℝ E) :
    ∫ x, A x * B x * (C x * D x) ∂μ =
      covarianceBilinDual μ A B * covarianceBilinDual μ C D +
      covarianceBilinDual μ A C * covarianceBilinDual μ B D +
      covarianceBilinDual μ A D * covarianceBilinDual μ B C := by
  have he (x : E) : 16*(A x*B x*(C x*D x)) =
      ((A+B) x)^2*((C+D) x)^2 - ((A-B) x)^2*((C+D) x)^2 -
      ((A+B) x)^2*((C-D) x)^2 + ((A-B) x)^2*((C-D) x)^2 := by
    simp only [add_apply, sub_apply]
    ring
  have h := integral_congr_ae (μ := μ) (ae_of_all _ he)
  have hi := integrable_square_product μ
  erw [integral_const_mul, integral_add
    (((hi (A+B) (C+D)).sub (hi (A-B) (C+D))).sub (hi (A+B) (C-D))) (hi (A-B) (C-D)),
    integral_sub ((hi (A+B) (C+D)).sub (hi (A-B) (C+D))) (hi (A+B) (C-D)),
    integral_sub (hi (A+B) (C+D)) (hi (A-B) (C+D))] at h
  simp only [gaussian_square_product μ h0] at h
  simp only [map_add, map_sub, add_apply, sub_apply, covarianceBilinDual_comm B A,
    covarianceBilinDual_comm D C] at h
  linarith only [h]

end Polarization

#print axioms gaussianReal_fourth
#print axioms gaussian_fourth_eq_variance
#print axioms gaussian_product_memLp_two
#print axioms gaussian_square_product
#print axioms gaussian_four_projections

end TimeFrequency
