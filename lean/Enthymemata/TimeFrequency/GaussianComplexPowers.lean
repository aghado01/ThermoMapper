import Enthymemata.TimeFrequency.GaussianFourthMoment
import Enthymemata.TimeFrequency.LogPowerCovariance

/-! # Power covariance of a centered, potentially improper complex Gaussian pair -/

open MeasureTheory ProbabilityTheory Complex
open scoped ComplexConjugate ProbabilityTheory

namespace TimeFrequency

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
  [SecondCountableTopology E] [MeasurableSpace E] [BorelSpace E]

private noncomputable def part (L : E →L[ℝ] ℂ) : Fin 2 → StrongDual ℝ E :=
  ![Complex.reCLM.comp L, Complex.imCLM.comp L]

private noncomputable def phaseFactor : Fin 2 → ℂ := ![1, I]

omit [CompleteSpace E] [SecondCountableTopology E] [MeasurableSpace E] [BorelSpace E] in
private theorem parts_sum (L : E →L[ℝ] ℂ) (x : E) :
    L x = ∑ i, phaseFactor i * (part L i x : ℂ) := by
  simp [Fin.sum_univ_two, part, phaseFactor, mul_comm I, Complex.re_add_im]

variable (μ : Measure E) [IsGaussian μ]

private theorem real_second (h0 : (∫ x, x ∂μ) = 0) (A B : StrongDual ℝ E) :
    ∫ x, A x*B x ∂μ = covarianceBilinDual μ A B := by
  rw [covarianceBilinDual_apply' IsGaussian.hasGaussianLaw_id.memLp_two]
  change _ = ∫ x, A (x - ∫ z, z ∂μ) * B (x - ∫ z, z ∂μ) ∂μ
  simp [h0]

private noncomputable def contraction (A B : E →L[ℝ] ℂ) (i j : Fin 2) : ℂ :=
  phaseFactor i * phaseFactor j * (covarianceBilinDual μ (part A i) (part B j) : ℂ)

private theorem complex_second (h0 : (∫ x, x ∂μ) = 0) (A B : E →L[ℝ] ℂ) :
    (∫ x, A x * B x ∂μ) = ∑ i, ∑ j, contraction μ A B i j := by
  let f (i j : Fin 2) (x : E) := (phaseFactor i * phaseFactor j) *
    ((part A i x * part B j x : ℝ) : ℂ)
  have he (x : E) : A x*B x = ∑ i, ∑ j, f i j x := by
    rw [parts_sum A x, parts_sum B x, Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro i _
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro j _
    simp only [f, Complex.ofReal_mul]
    ring
  have hi (i j) : Integrable (f i j) μ := by
    have h (L : StrongDual ℝ E) : MemLp (fun x ↦ L x) 2 μ :=
      (IsGaussian.hasGaussianLaw_id.map_fun L).memLp_two
    exact (Complex.ofRealCLM.integrable_comp ((h _).integrable_mul (h _))).const_mul _
  simp_rw [he]
  rw [integral_finsetSum _ fun i _ ↦ integrable_finsetSum _ fun j _ ↦ hi i j]
  simp_rw [integral_finsetSum _ fun j _ ↦ hi _ j, f, integral_const_mul, integral_complex_ofReal,
    real_second μ h0, contraction]

/-- Complexification of the real Gaussian fourth tensor. All four maps are real-linear,
so conjugate-linear observations and singular realifications are included. -/
theorem gaussian_complex_fourth (h0 : (∫ x, x ∂μ) = 0) (A B C D : E →L[ℝ] ℂ) :
    (∫ x, A x*B x*(C x*D x) ∂μ) =
      (∫ x, A x*B x ∂μ)*(∫ x, C x*D x ∂μ) +
      (∫ x, A x*C x ∂μ)*(∫ x, B x*D x ∂μ) +
      (∫ x, A x*D x ∂μ)*(∫ x, B x*C x ∂μ) := by
  let f (i j k l : Fin 2) (x : E) := (phaseFactor i*phaseFactor j*phaseFactor k*phaseFactor l) *
    ((part A i x * part B j x * (part C k x * part D l x) : ℝ) : ℂ)
  have he (x : E) : A x*B x*(C x*D x) = ∑ i, ∑ j, ∑ k, ∑ l, f i j k l x := by
    rw [parts_sum A x, parts_sum B x, parts_sum C x, parts_sum D x]
    simp only [Fin.sum_univ_two, f, Complex.ofReal_mul]
    ring
  have hi (i j k l) : Integrable (f i j k l) μ := by
    have h (U V : StrongDual ℝ E) : MemLp (fun x ↦ U x * V x) 2 μ :=
      gaussian_product_memLp_two (IsGaussian.hasGaussianLaw_id.map_fun U)
        (IsGaussian.hasGaussianLaw_id.map_fun V)
    exact (Complex.ofRealCLM.integrable_comp ((h _ _).integrable_mul (h _ _))).const_mul _
  have hi3 (i j k) := integrable_finsetSum Finset.univ fun l _ ↦ hi i j k l
  have hi2 (i j) := integrable_finsetSum Finset.univ fun k _ ↦ hi3 i j k
  have hi1 (i) := integrable_finsetSum Finset.univ fun j _ ↦ hi2 i j
  simp_rw [he, integral_finsetSum _ fun i _ ↦ hi1 i, integral_finsetSum _ fun j _ ↦ hi2 _ j,
    integral_finsetSum _ fun k _ ↦ hi3 _ _ k, integral_finsetSum _ fun l _ ↦ hi _ _ _ l,
    f, integral_const_mul, integral_complex_ofReal, gaussian_four_projections μ h0,
    complex_second μ h0]
  simp only [Fin.sum_univ_two, contraction, Complex.ofReal_add, Complex.ofReal_mul]
  ring

/-- S11 for any jointly Gaussian realification, with no circularity assumption. -/
theorem gaussian_power_covariance {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    {a b : Ω → ℂ} (h : HasGaussianLaw (fun ω ↦ (a ω, b ω)) P)
    (h0 : (∫ ω, (a ω, b ω) ∂P) = 0) :
    (∫ ω, ‖a ω‖^2*‖b ω‖^2 ∂P) - (∫ ω, ‖a ω‖^2 ∂P)*(∫ ω, ‖b ω‖^2 ∂P) =
      ‖∫ ω, a ω*conj (b ω) ∂P‖^2 + ‖∫ ω, a ω*b ω ∂P‖^2 := by
  let x := fun ω ↦ (a ω, b ω)
  let ν := P.map x
  letI : IsGaussian ν := h.isGaussian_map
  have hm : (∫ v, v ∂ν) = 0 := (integral_map h.aemeasurable (by fun_prop)).trans h0
  let A := ContinuousLinearMap.fst ℝ ℂ ℂ
  let B := ContinuousLinearMap.snd ℝ ℂ ℂ
  let J := Complex.conjCLE.toContinuousLinearMap
  have hm2 (U V : (ℂ × ℂ) →L[ℝ] ℂ) : (∫ v, U v*V v ∂ν) = ∫ ω, U (x ω)*V (x ω) ∂P :=
    integral_map h.aemeasurable (by fun_prop)
  have hp := gaussian_complex_fourth ν hm A (J.comp A) B (J.comp B)
  simp only [hm2] at hp
  rw [add_right_comm] at hp
  apply cov_normSq_pairing
  refine Eq.trans ?_ hp
  exact (integral_map h.aemeasurable (by fun_prop : AEStronglyMeasurable
    (fun v ↦ A v * (J.comp A) v * (B v*(J.comp B) v)) ν)).symm

theorem covariance_gaussian_powers {Ω : Type*} [MeasurableSpace Ω] {P : Measure Ω}
    {a b : Ω → ℂ} (h : HasGaussianLaw (fun ω ↦ (a ω, b ω)) P)
    (h0 : (∫ ω, (a ω, b ω) ∂P) = 0) :
    cov[fun ω ↦ ‖a ω‖^2, fun ω ↦ ‖b ω‖^2; P] =
      ‖∫ ω, a ω*conj (b ω) ∂P‖^2 + ‖∫ ω, a ω*b ω ∂P‖^2 := by
  letI := h.isProbabilityMeasure
  have hp {Y : Ω → ℂ} (hg : HasGaussianLaw Y P) : MemLp (fun ω ↦ ‖Y ω‖^2) 2 P := by
    apply (memLp_two_iff_integrable_sq (hg.memLp_two.aestronglyMeasurable.norm.pow 2)).mpr
    simpa only [Pi.pow_apply, ← pow_mul] using
      ((hg.memLp (p := 4) (by norm_num)).integrable_norm_pow (p := 4) (by decide))
  rw [covariance_eq_sub (hp h.fst) (hp h.snd)]
  exact gaussian_power_covariance h h0

/-- A standard real Gaussian embedded twice in ℂ has pseudo-covariance one and
power covariance two. This witnesses the contribution of the improper term. -/
theorem realGaussian_complex_witness :
    (∫ z, (z : ℂ)*(z : ℂ) ∂gaussianReal 0 1) = 1 ∧
    (∫ z, ‖(z : ℂ)‖^2*‖(z : ℂ)‖^2 ∂gaussianReal 0 1) -
      (∫ z, ‖(z : ℂ)‖^2 ∂gaussianReal 0 1)^2 = 2 := by
  let L := Complex.ofRealCLM.prod Complex.ofRealCLM
  have hg := (IsGaussian.hasGaussianLaw_id (μ := gaussianReal 0 1)).map_fun L
  have hm : (∫ z, L z ∂gaussianReal 0 1) = 0 := by
    refine (L.integral_comp_comm IsGaussian.hasGaussianLaw_id.integrable).trans ?_
    simp
  have hp := gaussian_power_covariance (a := fun z : ℝ ↦ (z : ℂ)) (b := fun z ↦ (z : ℂ)) hg hm
  have hs := variance_eq_integral (μ := gaussianReal 0 1) measurable_id.aemeasurable
  change Var[id; gaussianReal 0 1] = ∫ z, (z - ∫ u, u ∂gaussianReal 0 1)^2 ∂gaussianReal 0 1 at hs
  simp only [variance_id_gaussianReal, integral_id_gaussianReal, sub_zero, NNReal.coe_one] at hs
  have h2 : (∫ z, (z : ℂ)*(z : ℂ) ∂gaussianReal 0 1) = 1 := by
    simp_rw [← Complex.ofReal_mul, ← sq, integral_complex_ofReal, ← hs, Complex.ofReal_one]
  refine ⟨h2, ?_⟩
  simpa [sq, h2, one_add_one_eq_two] using hp

#print axioms gaussian_complex_fourth
#print axioms gaussian_power_covariance
#print axioms covariance_gaussian_powers
#print axioms realGaussian_complex_witness

end TimeFrequency
