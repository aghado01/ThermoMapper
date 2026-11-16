import Enthymemata.ActiveSubspace.ProjectionMean
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.Normed.Lp.Matrix
import Mathlib.MeasureTheory.Integral.Bochner.SumMeasure
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# Isotropic noise energy under a fixed orthogonal projection (S21)

The raw second moment matrix is the expectation of the outer product. The coordinate products
are integrable, so expectation commutes with the finite sums in a quadratic form.
-/

open Matrix MeasureTheory
open scoped ENNReal

namespace GenerativePrior

variable {Ω ι : Type*} [MeasurableSpace Ω] [Fintype ι] [DecidableEq ι]
  {μ : Measure Ω} [IsProbabilityMeasure μ] {η : Ω → EuclideanSpace ℝ ι} {v : ℝ}

/-- Each entry of the raw second moment matrix is integrable and its expectation is
the corresponding entry of `v I`, where `v = σ²` is the second moment per coordinate. -/
structure IsotropicSecondMoments (μ : Measure Ω) [IsProbabilityMeasure μ]
    (η : Ω → EuclideanSpace ℝ ι) (v : ℝ) : Prop where
  integrable (i j : ι) : Integrable (fun ω ↦ η ω i * η ω j) μ
  moment (i j : ι) : ∫ ω, η ω i * η ω j ∂μ = if i = j then v else 0

/-- Integrating the quadratic form against the second moments contracts it to a trace. -/
theorem IsotropicSecondMoments.integral_quadratic (h : IsotropicSecondMoments μ η v)
    (M : Matrix ι ι ℝ) :
    ∫ ω, (WithLp.ofLp (η ω)) ⬝ᵥ (M *ᵥ WithLp.ofLp (η ω)) ∂μ = v * M.trace := by
  have term (i j : ι) : Integrable (fun ω ↦ M i j * (η ω i * η ω j)) μ :=
    (h.integrable i j).const_mul _
  calc
    _ = ∫ ω, ∑ i, ∑ j, M i j * (η ω i * η ω j) ∂μ := by
      apply integral_congr_ae
      filter_upwards with ω
      simp only [dotProduct, mulVec, Finset.mul_sum]
      exact Finset.sum_congr rfl fun i _ ↦ Finset.sum_congr rfl fun j _ ↦ by
        change η ω i * (M i j * η ω j) = _
        ring
    _ = ∑ i, ∑ j, M i j * (if i = j then v else 0) := by
      rw [integral_finsetSum _ fun i _ ↦ integrable_finsetSum _ fun j _ ↦ term i j]
      apply Finset.sum_congr rfl
      intro i _
      rw [integral_finsetSum _ fun j _ ↦ term i j]
      simp_rw [integral_const_mul, h.moment]
    _ = _ := by simp [mul_ite, Matrix.trace, mul_comm v, Finset.sum_mul]

/-- The mean squared Euclidean norm of the noise in the ambient space. -/
theorem IsotropicSecondMoments.ambient_energy (h : IsotropicSecondMoments μ η v) :
    ∫ ω, ‖η ω‖ ^ 2 ∂μ = v * Fintype.card ι := by
  simpa only [one_mulVec, trace_one, dotProduct, ← sq, EuclideanSpace.real_norm_sq_eq]
    using h.integral_quadratic (1 : Matrix ι ι ℝ)

/-- The mean squared Euclidean norm after a fixed orthogonal projection of rank `P.rank`. -/
theorem IsotropicSecondMoments.projected_energy (h : IsotropicSecondMoments μ η v)
    {P : Matrix ι ι ℝ} (hP : ActiveSubspace.IsOrthProj P) :
    ∫ ω, ‖Matrix.toLpLin 2 2 P (η ω)‖ ^ 2 ∂μ = v * P.rank := by
  rw [hP.rank_eq_trace, ← h.integral_quadratic P]
  apply integral_congr_ae
  filter_upwards with ω
  rw [EuclideanSpace.real_norm_sq_eq]
  change ∑ i, ((P *ᵥ WithLp.ofLp (η ω)) i) ^ 2 = _
  simp only [sq]
  change (P *ᵥ _) ⬝ᵥ (P *ᵥ _) = _
  rw [← dotProduct_transpose_mulVec, mulVec_mulVec, hP.1, hP.2]

/-- In positive ambient dimension, a rank-`k` projection retains `k/n` of the mean energy. -/
theorem IsotropicSecondMoments.projected_fraction [Nonempty ι] (h : IsotropicSecondMoments μ η v)
    {P : Matrix ι ι ℝ} (hP : ActiveSubspace.IsOrthProj P) :
    (∫ ω, ‖Matrix.toLpLin 2 2 P (η ω)‖ ^ 2 ∂μ) =
      ((P.rank : ℝ) / Fintype.card ι) * ∫ ω, ‖η ω‖ ^ 2 ∂μ := by
  rw [h.projected_energy hP, h.ambient_energy]
  have hc : (Fintype.card ι : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  field_simp [hc]

/-! ### A rank-one projection of four equally likely signed coordinate vectors -/

namespace SignedAxes

/-- Each of the four outcomes has probability one quarter. -/
noncomputable def law : Measure (Fin 4) := (4 : ℝ≥0∞)⁻¹ • Measure.count

instance law_isProbabilityMeasure : IsProbabilityMeasure law := by
  constructor
  norm_num [law, Measure.smul_apply, Measure.count_univ, ENat.card_eq_coe_fintype_card]
  exact ENNReal.inv_mul_cancel (by norm_num) (by norm_num)

/-- The outcomes are the two signed coordinate axes, each with length two. -/
def noise (ω : Fin 4) : EuclideanSpace ℝ (Fin 2) :=
  WithLp.toLp 2 (![![2, 0], ![-2, 0], ![0, 2], ![0, -2]] ω)

theorem isotropic : IsotropicSecondMoments law noise 2 := by
  constructor
  · exact fun _ _ ↦ .of_finite
  · intro i j
    rw [law, integral_smul_measure, integral_count]
    fin_cases i <;> fin_cases j <;>
      norm_num [noise, PiLp.toLp_apply, Fin.sum_univ_succ, smul_eq_mul]

/-- Projection onto the first coordinate. -/
def proj : Matrix (Fin 2) (Fin 2) ℝ := diagonal ![1, 0]

theorem orthogonal : ActiveSubspace.IsOrthProj proj := by
  constructor
  · exact diagonal_transpose _
  · rw [proj, diagonal_mul_diagonal]
    congr 1
    funext i
    fin_cases i <;> norm_num

theorem rank_one : proj.rank = 1 := by
  have hr := orthogonal.rank_eq_trace
  norm_num [proj, Matrix.trace, Fin.sum_univ_succ] at hr
  exact_mod_cast hr

/-- The rank-one projection retains mean energy two of the ambient four. -/
theorem energies :
    (∫ ω, ‖Matrix.toLpLin 2 2 proj (noise ω)‖ ^ 2 ∂law) = 2 ∧
      (∫ ω, ‖noise ω‖ ^ 2 ∂law) = 4 := by
  constructor
  · simpa [rank_one] using isotropic.projected_energy orthogonal
  · rw [isotropic.ambient_energy]
    norm_num

end SignedAxes

#print axioms IsotropicSecondMoments.integral_quadratic
#print axioms IsotropicSecondMoments.ambient_energy
#print axioms IsotropicSecondMoments.projected_energy
#print axioms IsotropicSecondMoments.projected_fraction
#print axioms SignedAxes.law_isProbabilityMeasure
#print axioms SignedAxes.isotropic
#print axioms SignedAxes.orthogonal
#print axioms SignedAxes.rank_one
#print axioms SignedAxes.energies

end GenerativePrior
