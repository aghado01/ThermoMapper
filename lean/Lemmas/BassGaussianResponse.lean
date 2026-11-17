import Lemmas.BassWhitenedMarginal
import Mathlib.Probability.Distributions.Gaussian.HasGaussianLaw.Independence
import Mathlib.Probability.Independence.Integration

/-! # The Gaussian law of the BASS linear response -/

open Matrix MeasureTheory ProbabilityTheory
open scoped ProbabilityTheory

namespace BASS

variable {ι κ Ω : Type*} [Fintype ι] [Fintype κ] [MeasurableSpace Ω]
  {μ : Measure Ω} {β : Ω → κ → ℝ} {ε : Ω → ι → ℝ}

omit [Fintype ι] in
/-- A linear image of the joint Gaussian remains Gaussian, including singular laws. -/
theorem gaussian_linear_response [Finite ι] (Z : Matrix ι κ ℝ)
    (h : HasGaussianLaw (fun ω ↦ (β ω, ε ω)) μ) :
    HasGaussianLaw (fun ω ↦ Z *ᵥ β ω + ε ω) μ := by
  letI := Fintype.ofFinite ι
  let L := Z.mulVecLin.toContinuousLinearMap
  exact h.map_fun (L.comp (.fst ℝ (κ → ℝ) (ι → ℝ)) + .snd ℝ (κ → ℝ) (ι → ℝ))

private theorem centered_coordinate {j : Type*} [Fintype j] {X : Ω → j → ℝ}
    (hX : HasGaussianLaw X μ) (h0 : ∫ ω, X ω ∂μ = 0) (i : j) :
    (∫ ω, X ω i ∂μ) = 0 := by
  let L : (j → ℝ) →L[ℝ] ℝ := ContinuousLinearMap.proj i
  exact (L.integral_comp_comm hX.integrable).trans (by rw [h0, map_zero])

/-- S14/S4 with independent centered Gaussian coefficients and noise and the specified
marginal covariance matrices. The response is centered and its covariance is S14's matrix. -/
theorem gaussian_marginal_response [DecidableEq ι] [DecidableEq κ]
    (Z : Matrix ι κ ℝ) (Om : Matrix ι ι ℝ) (σ2 τ : ℝ)
    (hβ : HasGaussianLaw β μ) (hε : HasGaussianLaw ε μ) (hind : IndepFun β ε μ)
    (mβ : ∫ ω, β ω ∂μ = 0) (mε : ∫ ω, ε ω ∂μ = 0)
    (cβ : ∀ a b, cov[fun ω ↦ β ω a, fun ω ↦ β ω b; μ] =
      ((σ2 * τ⁻¹) • (Zᵀ * Om⁻¹ * Z)⁻¹) a b)
    (cε : ∀ i j, cov[fun ω ↦ ε ω i, fun ω ↦ ε ω j; μ] = (σ2 • Om) i j) :
    HasGaussianLaw (fun ω ↦ Z *ᵥ β ω + ε ω) μ ∧ (∫ ω, Z *ᵥ β ω + ε ω ∂μ) = 0 ∧
      ∀ i j, cov[fun ω ↦ (Z *ᵥ β ω + ε ω) i, fun ω ↦ (Z *ᵥ β ω + ε ω) j; μ] =
        (σ2 • (Om + τ⁻¹ • (Z * (Zᵀ * Om⁻¹ * Z)⁻¹ * Zᵀ))) i j := by
  letI := hβ.isProbabilityMeasure
  let L := Z.mulVecLin.toContinuousLinearMap
  have hR := gaussian_linear_response Z (hind.hasGaussianLaw hβ hε)
  have hM : (∫ ω, Z *ᵥ β ω + ε ω ∂μ) = 0 := by
    change (∫ ω, L (β ω) + ε ω ∂μ) = 0
    rw [integral_add (hβ.map_fun L).integrable hε.integrable, L.integral_comp_comm hβ.integrable,
      mβ, mε, map_zero, add_zero]
  refine ⟨hR,hM,fun i j ↦ ?_⟩
  rw [covariance_eq_sub (hR.eval i).memLp_two (hR.eval j).memLp_two,
    centered_coordinate hR hM i, zero_mul, sub_zero]
  apply marginal_covariance μ Z Om β ε
  · intro a b
    have ha := (hβ.eval a).memLp_two
    have hb := (hβ.eval b).memLp_two
    refine ⟨ha.integrable_mul hb, ?_⟩
    simpa only [covariance_eq_sub ha hb, Pi.mul_apply, centered_coordinate hβ mβ a,
      zero_mul, sub_zero]
      using cβ a b
  · intro i j
    have hi := (hε.eval i).memLp_two
    have hj := (hε.eval j).memLp_two
    refine ⟨hi.integrable_mul hj, ?_⟩
    simpa only [covariance_eq_sub hi hj, Pi.mul_apply, centered_coordinate hε mε i,
      zero_mul, sub_zero]
      using cε i j
  · intro a i
    have ha := (hβ.eval a).memLp_two
    have hi := (hε.eval i).memLp_two
    refine ⟨ha.integrable_mul hi, ?_⟩
    have h := (hind.comp (measurable_pi_apply a) (measurable_pi_apply i)
      ).integral_mul_eq_mul_integral ha.aestronglyMeasurable hi.aestronglyMeasurable
    exact h.trans (by simp [centered_coordinate hβ mβ a])

/-- One coefficient, identity design and independent unit Gaussian noise give N(0,2). -/
theorem unit_gaussian_response :
    ((gaussianReal 0 1).prod (gaussianReal 0 1)).map (fun p : ℝ × ℝ ↦
      (((1 : Matrix Unit Unit ℝ) *ᵥ (fun _ ↦ p.1) + (fun _ : Unit ↦ p.2)) : Unit → ℝ) ()) =
      gaussianReal 0 2 := by
  simp only [Matrix.one_mulVec, Pi.add_apply]
  have hi : IndepFun (Prod.fst : ℝ × ℝ → ℝ) Prod.snd
      ((gaussianReal 0 1).prod (gaussianReal 0 1)) := indepFun_prod measurable_id measurable_id
  have h := gaussianReal_add_gaussianReal_of_indepFun hi
    measurePreserving_fst.map_eq measurePreserving_snd.map_eq
  change ((gaussianReal 0 1).prod (gaussianReal 0 1)).map
    ((Prod.fst : ℝ × ℝ → ℝ) + Prod.snd) = gaussianReal 0 2
  simpa only [zero_add, one_add_one_eq_two] using h

#print axioms gaussian_linear_response
#print axioms gaussian_marginal_response
#print axioms unit_gaussian_response

end BASS
