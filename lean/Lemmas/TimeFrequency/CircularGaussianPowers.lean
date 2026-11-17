import Lemmas.TimeFrequency.LogPowerLaplaceLaw
import Mathlib.Probability.Distributions.Gaussian.Real
import Mathlib.Probability.HasLaw

/-! # The powers of a correlated unit circular complex Gaussian pair -/

open Real Set MeasureTheory ProbabilityTheory
open scoped NNReal ProbabilityTheory

namespace TimeFrequency

noncomputable abbrev halfGaussian := gaussianReal 0 (1/2)
noncomputable abbrev gaussianPair := halfGaussian.prod halfGaussian
noncomputable abbrev circularPairLaw := gaussianPair.prod gaussianPair

/-- The Cartesian coordinates of two independent circular complex innovations are interleaved
so their real and imaginary blocks can be integrated independently. -/
noncomputable def complexSample (p : (ℝ × ℝ) × (ℝ × ℝ)) : ℂ := ⟨p.1.1, p.2.1⟩

noncomputable def correlatedSample (c d : ℝ) (p : (ℝ × ℝ) × (ℝ × ℝ)) : ℂ :=
  ⟨c*p.1.1 + d*p.1.2, c*p.2.1 + d*p.2.2⟩

theorem halfGaussian_density (f : ℝ → ℝ) :
    ∫ x, f x ∂halfGaussian = ∫ x, (√π)⁻¹ * exp (-x^2) * f x := by
  rw [integral_gaussianReal_eq_integral_smul (by norm_num : (1/2 : ℝ≥0) ≠ 0)]
  norm_num [gaussianPDFReal, smul_eq_mul]
  congr 1
  funext x
  field_simp

/-- Integrating a shifted Gaussian with quadratic tilt, including zero loading. -/
theorem halfGaussian_tilt (a d t : ℝ) (ht : 0 ≤ t) :
    (∫ v, exp (-t*(a+d*v)^2) ∂halfGaussian) =
      (√(1+t*d^2))⁻¹ * exp (-(t*a^2/(1+t*d^2))) := by
  let B := 1+t*d^2
  have hB : 0 < B := by dsimp [B]; positivity
  rw [halfGaussian_density]
  have hcomp (v : ℝ) : (√π)⁻¹ * exp (-v^2) * exp (-t*(a+d*v)^2) =
      ((√π)⁻¹ * exp (-(t*a^2/B))) * exp (-B*(v + t*a*d/B)^2) := by
    rw [mul_assoc, ← exp_add, mul_assoc, ← exp_add]
    congr 2
    dsimp [B]
    field_simp
    ring
  simp_rw [hcomp]
  rw [integral_const_mul, show (∫ v : ℝ, exp (-B*(v + t*a*d/B)^2)) = √(π/B) by
    rw [← integral_gaussian]
    exact integral_add_right_eq_self (fun v : ℝ ↦ exp (-B*v^2)) _]
  rw [Real.sqrt_div pi_nonneg]
  change _ = (√B)⁻¹ * exp (-(t*a^2/B))
  field_simp

theorem realBlock_laplace {c d : ℝ} (hcd : c ^ 2 + d ^ 2 = 1) {s t : ℝ} (hs : 0 ≤ s) (ht : 0 ≤ t) :
    (∫ p : ℝ × ℝ, exp (-s*p.1^2) * exp (-t*(c*p.1+d*p.2)^2) ∂gaussianPair) =
      (√((1+s)*(1+t)-c^2*s*t))⁻¹ := by
  let B := 1+t*d^2
  let A := s+t*c^2/B
  have hB : 0 < B := by dsimp [B]; positivity
  have hA : 0 ≤ A := by dsimp [A]; positivity
  have hdet : B*(1+A) = (1+s)*(1+t)-c^2*s*t := by
    dsimp [B,A]
    field_simp
    nlinarith only [hcd, congrArg (fun z : ℝ ↦ s*t*z) hcd, congrArg (fun z : ℝ ↦ t*z) hcd]
  have hi : Integrable (fun p : ℝ × ℝ ↦ exp (-s*p.1^2) * exp (-t*(c*p.1+d*p.2)^2))
      gaussianPair := by
    apply Integrable.of_bound (by fun_prop) 1
    apply ae_of_all
    intro p
    rw [← exp_add, Real.norm_of_nonneg (exp_pos _).le, exp_le_one_iff]
    nlinarith [mul_nonneg hs (sq_nonneg p.1), mul_nonneg ht (sq_nonneg (c*p.1+d*p.2))]
  rw [integral_prod _ hi]
  simp_rw [integral_const_mul, halfGaussian_tilt _ _ _ ht]
  have he (u : ℝ) : exp (-s*u^2) * ((√B)⁻¹ * exp (-(t*(c*u)^2/B))) =
      (√B)⁻¹ * exp (-A*u^2) := by
    rw [mul_left_comm, ← exp_add]
    congr 2
    dsimp [A]
    ring
  change ∫ u, exp (-s*u^2) * ((√B)⁻¹ * exp (-(t*(c*u)^2/B))) ∂halfGaussian = _
  simp_rw [he]
  rw [integral_const_mul]
  have hv := halfGaussian_tilt 0 1 A hA
  simp only [one_mul, zero_add, one_pow, mul_one, zero_pow (by decide : 2 ≠ 0), mul_zero,
    zero_div, neg_zero, exp_zero] at hv
  rw [hv, ← mul_inv, ← Real.sqrt_mul hB.le, hdet]

theorem circularPower_laplace {c d : ℝ} (hcd : c ^ 2 + d ^ 2 = 1) {s t : ℝ}
    (hs : 0 ≤ s) (ht : 0 ≤ t) :
    (∫ p, exp (-Complex.normSq (complexSample p)*s) *
      exp (-Complex.normSq (correlatedSample c d p)*t) ∂circularPairLaw) =
      1 / ((1+s)*(1+t)-c^2*s*t) := by
  let f := fun p : ℝ × ℝ ↦ exp (-s*p.1^2) * exp (-t*(c*p.1+d*p.2)^2)
  have he (p) : exp (-Complex.normSq (complexSample p)*s) *
      exp (-Complex.normSq (correlatedSample c d p)*t) = f p.1 * f p.2 := by
    simp only [complexSample, correlatedSample, Complex.normSq_apply, f, ← exp_add]
    congr 1
    ring
  simp_rw [he]
  rw [integral_prod_mul f f, show (∫ p, f p ∂gaussianPair) = _ from realBlock_laplace hcd hs ht,
    ← mul_inv, Real.mul_self_sqrt]
  · simp
  · nlinarith only [hs, ht, mul_nonneg hs ht, mul_nonneg (mul_nonneg hs ht) (sq_nonneg d),
      congrArg (fun z : ℝ ↦ s*t*z) hcd]

/-- The pair has an arbitrary correlation phase. Its squared coherence is `c²`. -/
noncomputable def gaussianSamples (u : ℂ) (c d : ℝ) (p : (ℝ × ℝ) × (ℝ × ℝ)) : ℂ × ℂ :=
  (complexSample p, u * correlatedSample c d p)

theorem gaussianSamples_measurable (u : ℂ) (c d : ℝ) : Measurable (gaussianSamples u c d) := by
  unfold gaussianSamples complexSample correlatedSample
  simp only [Complex.mk_eq_add_mul_I]
  apply Continuous.measurable
  fun_prop

/-- S10c on any probability space carrying the jointly circular Gaussian sampling law.
The joint law is a linear image of four independent real Gaussians of variance one half. -/
theorem circularGaussian_logPower_covariance {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω}
    [IsProbabilityMeasure μ] {a b : Ω → ℂ} (ha : Measurable a) (hb : Measurable b)
    {c d : ℝ} (hcd : c ^ 2 + d ^ 2 = 1) {u : ℂ} (hu : Complex.normSq u = 1)
    (h : HasLaw (fun ω ↦ (a ω, b ω)) (circularPairLaw.map (gaussianSamples u c d)) μ) :
    HasSum (fun n : ℕ ↦ (c^2)^(n+1)/((n:ℝ)+1)^2)
      (cov[fun ω ↦ log (Complex.normSq (a ω)), fun ω ↦ log (Complex.normSq (b ω)); μ]) := by
  apply covariance_log_of_laplace ⟨sq_nonneg _, by nlinarith [sq_nonneg d]⟩
    (Complex.continuous_normSq.measurable.comp ha) (Complex.continuous_normSq.measurable.comp hb)
    (ae_of_all _ fun _ ↦ Complex.normSq_nonneg _) (ae_of_all _ fun _ ↦ Complex.normSq_nonneg _)
  intro s t hs ht
  let f := fun z : ℂ × ℂ ↦ exp (-Complex.normSq z.1*s) * exp (-Complex.normSq z.2*t)
  have hf : Measurable f := by dsimp [f]; fun_prop
  have hm := gaussianSamples_measurable u c d
  refine (h.integral_comp hf.aestronglyMeasurable).trans ?_
  rw [integral_map hm.aemeasurable hf.aestronglyMeasurable]
  simpa only [f, gaussianSamples, Complex.normSq_mul, hu, one_mul] using
    circularPower_laplace hcd hs ht

theorem canonical_logPower_covariance {c d : ℝ} (hcd : c ^ 2 + d ^ 2 = 1) :
    HasSum (fun n : ℕ ↦ (c^2)^(n+1)/((n:ℝ)+1)^2)
      (cov[fun p ↦ log (Complex.normSq (complexSample p)),
        fun p ↦ log (Complex.normSq (correlatedSample c d p)); circularPairLaw]) := by
  apply covariance_log_of_laplace ⟨sq_nonneg _, by nlinarith [sq_nonneg d]⟩
    (by simp only [complexSample, Complex.normSq_apply]; fun_prop)
    (by simp only [correlatedSample, Complex.normSq_apply]; fun_prop)
    (ae_of_all _ fun _ ↦ Complex.normSq_nonneg _) (ae_of_all _ fun _ ↦ Complex.normSq_nonneg _)
  exact fun _ _ ↦ circularPower_laplace hcd

/-- At perfect coherence the two powers agree pointwise. Their log-power variance is ζ(2). -/
theorem circularGaussian_logPower_variance :
    Var[fun p ↦ log (Complex.normSq (complexSample p)); circularPairLaw] = π^2/6 := by
  have he : correlatedSample 1 0 = complexSample := by
    funext p
    simp [correlatedSample, complexSample]
  have h := canonical_logPower_covariance (c := 1) (d := 0) (by norm_num)
  rw [he, covariance_self (by
    apply Measurable.aemeasurable
    simp only [complexSample, Complex.normSq_apply]
    fun_prop), one_pow] at h
  exact h.unique dilog_one

example := canonical_logPower_covariance (c := 3/5) (d := 4/5) (by norm_num)

#print axioms gaussianSamples_measurable
#print axioms circularGaussian_logPower_covariance
#print axioms canonical_logPower_covariance
#print axioms circularGaussian_logPower_variance
#print axioms halfGaussian_density
#print axioms halfGaussian_tilt
#print axioms realBlock_laplace
#print axioms circularPower_laplace

end TimeFrequency
