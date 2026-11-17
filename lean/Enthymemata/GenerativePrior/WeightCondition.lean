import Enthymemata.GenerativePrior.AngleMap
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Geometry.Euclidean.Angle.Unoriented.Basic
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-! # The angular second moment of two active ReLU gates (S18) -/

open Real Matrix Set MeasureTheory InnerProductGeometry
open scoped RealInnerProductSpace

namespace GenerativePrior

noncomputable def planeDirection (φ : ℝ) : Fin 2 → ℝ := ![cos φ, sin φ]

/-- The planar formula for the weight-condition matrix in the basis of the two directions. -/
noncomputable def weightQ (θ : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  ((π - θ) / (2 * π)) • 1 + (sin θ / (2 * π)) • !![cos θ, sin θ; sin θ, -cos θ]

/-- The integral of the outer product, taken entry by entry. -/
noncomputable def angularMoment (θ : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := fun i j ↦
  (1 / π) * ∫ φ in θ - π / 2..π / 2, planeDirection φ i * planeDirection φ j

/-- Within the polar chart, both dot products are positive exactly on the open arc. -/
theorem both_gates_iff {φ θ : ℝ} (hφ : φ ∈ Ioo (-π) π) (hθ : θ ∈ Icc 0 π) :
    0 < cos φ ∧ 0 < cos φ * cos θ + sin φ * sin θ ↔ φ ∈ Ioo (θ - π / 2) (π / 2) := by
  rw [← cos_sub]
  obtain ⟨hl, hr⟩ := hφ
  obtain ⟨ht, hT⟩ := hθ
  constructor
  · rintro ⟨hc, hd⟩
    have hupper : φ < π / 2 := by
      by_contra! hn
      exact hc.not_ge (cos_nonpos_of_pi_div_two_le_of_le hn (by linarith [pi_pos]))
    have hlower : -(π / 2) < φ := by
      by_contra! hn
      have h := cos_nonpos_of_pi_div_two_le_of_le (x := -φ) (by linarith) (by linarith [pi_pos])
      exact hc.not_ge (by simpa using h)
    refine ⟨?_, hupper⟩
    by_contra! hn
    have h := cos_nonpos_of_pi_div_two_le_of_le (x := -(φ - θ)) (by linarith) (by linarith)
    exact hd.not_ge (by rwa [cos_neg] at h)
  · rintro ⟨ha, hb⟩
    exact ⟨cos_pos_of_mem_Ioo ⟨by linarith, hb⟩, cos_pos_of_mem_Ioo ⟨by linarith, by linarith⟩⟩

theorem angularMoment_eq_weightQ (θ : ℝ) : angularMoment θ = weightQ θ := by
  have hcc := integral_cos_sq (a := θ - π / 2) (b := π / 2)
  have hss := integral_sin_sq (a := θ - π / 2) (b := π / 2)
  have hsc := integral_sin_mul_cos₂ (a := θ - π / 2) (b := π / 2)
  simp only [cos_sub_pi_div_two, sin_sub_pi_div_two, cos_pi_div_two, sin_pi_div_two,
    zero_mul, mul_zero, zero_pow (by decide : 2 ≠ 0)] at hcc hss hsc
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [angularMoment, weightQ, planeDirection, ← sq, hcc, hss, hsc,
      mul_comm (cos _) (sin _)] <;> ring

theorem weightQ_zero : weightQ 0 = (1 / 2 : ℝ) • 1 := by
  simp [weightQ]
  field_simp

theorem weightQ_pi : weightQ π = 0 := by simp [weightQ]

theorem weightQ_cross (θ : ℝ) :
    planeDirection 0 ⬝ᵥ (weightQ θ *ᵥ planeDirection θ) = reluAngleCos θ / 2 := by
  simp [weightQ, planeDirection, Matrix.mulVec, dotProduct, Fin.sum_univ_two, Matrix.one_apply,
    reluAngleCos]
  field_simp
  nlinarith only [congrArg (fun t : ℝ ↦ sin θ * t) (sin_sq_add_cos_sq θ)]

theorem weightQ_self : planeDirection 0 ⬝ᵥ (weightQ 0 *ᵥ planeDirection 0) = 1 / 2 := by
  rw [weightQ_cross, reluAngleCos_zero]

/-- Exact unit-input second moments give the S17 output angle. -/
theorem exact_layer_angle {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    {θ : ℝ} {u v : E} (hu : ‖u‖ ^ 2 = 1 / 2) (hv : ‖v‖ ^ 2 = 1 / 2)
    (huv : ⟪u, v⟫ = planeDirection 0 ⬝ᵥ (weightQ θ *ᵥ planeDirection θ)) :
    angle u v = reluAngle θ := by
  have he : ‖u‖ = ‖v‖ := (sq_eq_sq₀ (norm_nonneg _) (norm_nonneg _)).mp (hu.trans hv.symm)
  have hp : ‖u‖ * ‖v‖ = 1 / 2 := by rwa [← he, ← sq]
  rw [angle, huv, weightQ_cross, hp]
  unfold reluAngle
  congr 1
  ring

section FiniteLayer

variable {ι : Type*} [Fintype ι] (w : ι → Fin 2 → ℝ)

noncomputable def planarRelu (x : Fin 2 → ℝ) : EuclideanSpace ℝ ι :=
  WithLp.toLp 2 (fun i ↦ max 0 (w i ⬝ᵥ x))

/-- Sum of the row outer products for which both input gates are open. -/
noncomputable def gatedGram (x y : Fin 2 → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  ∑ i, if 0 < w i ⬝ᵥ x ∧ 0 < w i ⬝ᵥ y then vecMulVec (w i) (w i) else 0

theorem planarRelu_inner (x y : Fin 2 → ℝ) :
    ⟪planarRelu w x, planarRelu w y⟫ = x ⬝ᵥ (gatedGram w x y *ᵥ y) := by
  rw [PiLp.inner_apply, gatedGram, sum_mulVec, dotProduct_sum]
  apply Finset.sum_congr rfl
  intro i _
  change max 0 (w i ⬝ᵥ y) * max 0 (w i ⬝ᵥ x) = _
  by_cases hx : 0 < w i ⬝ᵥ x <;> by_cases hy : 0 < w i ⬝ᵥ y
  · rw [if_pos ⟨hx, hy⟩, vecMulVec_mulVec, op_smul_eq_smul, dotProduct_smul, smul_eq_mul,
      max_eq_right hx.le, max_eq_right hy.le, dotProduct_comm x (w i)]
  · rw [if_neg (fun h ↦ hy h.2), max_eq_left (le_of_not_gt hy), zero_mul,
      zero_mulVec, dotProduct_zero]
  all_goals
    rw [if_neg (fun h ↦ hx h.1), max_eq_left (le_of_not_gt hx), mul_zero,
      zero_mulVec, dotProduct_zero]

theorem planeDirection_self (θ : ℝ) :
    planeDirection θ ⬝ᵥ planeDirection θ = 1 := by
  simp [planeDirection, dotProduct, Fin.sum_univ_two, ← sq]

/-- The three exact WDC matrices for this input pair imply the S17 output angle. -/
theorem planarRelu_angle_of_exact_gram (θ : ℝ)
    (hxx : gatedGram w (planeDirection 0) (planeDirection 0) = weightQ 0)
    (hyy : gatedGram w (planeDirection θ) (planeDirection θ) = weightQ 0)
    (hxy : gatedGram w (planeDirection 0) (planeDirection θ) = weightQ θ) :
    angle (planarRelu w (planeDirection 0)) (planarRelu w (planeDirection θ)) = reluAngle θ := by
  apply exact_layer_angle
  · rw [← real_inner_self_eq_norm_sq, planarRelu_inner, hxx, weightQ_self]
  · rw [← real_inner_self_eq_norm_sq, planarRelu_inner, hyy, weightQ_zero, smul_mulVec,
      one_mulVec, dotProduct_smul, planeDirection_self, smul_eq_mul, mul_one]
  · rw [planarRelu_inner, hxy]

end FiniteLayer

/-- At right angles, gating selects one quadrant and the off-diagonal entries survive. -/
theorem angularMoment_right_angle :
    angularMoment (π / 2) = !![(1 / 4 : ℝ), 1 / (2 * π); 1 / (2 * π), 1 / 4] := by
  rw [angularMoment_eq_weightQ]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [weightQ] <;> field_simp <;> ring

#print axioms both_gates_iff
#print axioms angularMoment_eq_weightQ
#print axioms weightQ_zero
#print axioms weightQ_pi
#print axioms weightQ_cross
#print axioms weightQ_self
#print axioms exact_layer_angle
#print axioms planarRelu_inner
#print axioms planeDirection_self
#print axioms planarRelu_angle_of_exact_gram
#print axioms angularMoment_right_angle

end GenerativePrior
