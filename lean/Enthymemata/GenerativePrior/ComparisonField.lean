import Enthymemata.GenerativePrior.AngleMap
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Algebra.Order.BigOperators.GroupWithZero.Finset
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Geometry.Euclidean.Angle.Unoriented.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# Zeros of the deterministic comparison field (S20)

Hand and Voroninski, proof sketch and Lemma `lemma:Seps`. The angle at layer `i` is
`reluAngle^[i] θ`. The reference code is the last argument of `comparisonField`.
-/

open Real Set Function Filter InnerProductGeometry NormedSpace
open scoped Topology

namespace GenerativePrior

noncomputable def angleWeight (θ : ℝ) : ℝ := (π - θ) / π

/-- The coefficient ξ multiplying the reference direction. -/
noncomputable def angleProduct (θ : ℝ) (d : ℕ) : ℝ :=
  ∏ i ∈ Finset.range d, angleWeight (reluAngle^[i] θ)

/-- The coefficient ζ. The tail product starts strictly after layer i. -/
noncomputable def angleSum (θ : ℝ) (d : ℕ) : ℝ :=
  ∑ i ∈ Finset.range d, sin (reluAngle^[i] θ) / π *
    ∏ j ∈ Finset.Ico (i + 1) d, angleWeight (reluAngle^[j] θ)

noncomputable def antipodalRadius (d : ℕ) : ℝ := angleSum π d

theorem angleSum_succ (θ : ℝ) (d : ℕ) :
    angleSum θ (d + 1) = angleSum θ d * angleWeight (reluAngle^[d] θ) +
      sin (reluAngle^[d] θ) / π := by
  simp only [angleSum, Finset.sum_range_succ, Finset.Ico_self, Finset.prod_empty, mul_one]
  rw [Finset.sum_mul]
  congr 1
  apply Finset.sum_congr rfl
  intro i hi
  rw [Finset.prod_Ico_succ_top (Finset.mem_range.mp hi), mul_assoc]

/-- Expansion of the affine cosine recurrence, including its initial-value term. -/
theorem cos_iterate_eq_angleProduct_add_angleSum {θ : ℝ} (hθ : θ ∈ Icc 0 π) (d : ℕ) :
    cos (reluAngle^[d] θ) = angleProduct θ d * cos θ + angleSum θ d := by
  induction d with
  | zero => simp [angleProduct, angleSum]
  | succ d ih =>
    rw [iterate_succ_apply', cos_reluAngle (reluAngle_mapsTo.iterate d hθ), angleSum_succ]
    rw [show angleProduct θ (d + 1) = angleProduct θ d * angleWeight (reluAngle^[d] θ)
      from Finset.prod_range_succ ..]
    unfold reluAngleCos angleWeight
    rw [ih]
    ring

theorem reluAngle_iterate_zero (n : ℕ) : reluAngle^[n] 0 = 0 :=
  (show IsFixedPt reluAngle 0 from reluAngle_zero).iterate n

theorem angleProduct_at_zero (d : ℕ) : angleProduct 0 d = 1 := by
  simp [angleProduct, reluAngle_iterate_zero, angleWeight, pi_ne_zero]

theorem angleSum_at_zero (d : ℕ) : angleSum 0 d = 0 := by
  simp [angleSum, reluAngle_iterate_zero]

theorem angleProduct_pos {θ : ℝ} (hθ : θ ∈ Ico 0 π) (d : ℕ) : 0 < angleProduct θ d := by
  apply Finset.prod_pos
  intro i _
  have hi : reluAngle^[i] θ ≤ θ := reluAngle_iterate_antitone ⟨hθ.1, hθ.2.le⟩ (Nat.zero_le i)
  exact div_pos (sub_pos.mpr (hi.trans_lt hθ.2)) pi_pos

theorem angleProduct_at_pi {d : ℕ} (hd : 1 ≤ d) : angleProduct π d = 0 :=
  Finset.prod_eq_zero (Finset.mem_range.mpr (show 0 < d from hd)) (by simp [angleWeight])

theorem antipodalRadius_eq_cos {d : ℕ} (hd : 1 ≤ d) :
    antipodalRadius d = cos (reluAngle^[d] π) := by
  rw [cos_iterate_eq_angleProduct_add_angleSum (right_mem_Icc.mpr pi_pos.le), angleProduct_at_pi hd]
  simp [antipodalRadius]

theorem antipodalRadius_two : antipodalRadius 2 = 1 / π := by
  rw [antipodalRadius_eq_cos (by decide), show reluAngle^[2] π = reluAngle (reluAngle π) by rfl,
    reluAngle_pi, cos_reluAngle ⟨(half_pos pi_pos).le, half_le_self pi_pos.le⟩]
  simp [reluAngleCos]

theorem antipodalRadius_pos {d : ℕ} (hd : 2 ≤ d) : 0 < antipodalRadius d := by
  have hπ := right_mem_Icc.mpr pi_pos.le
  have h := cos_le_cos_of_nonneg_of_le_pi (reluAngle_mapsTo.iterate d hπ).1
    (reluAngle_mapsTo.iterate 2 hπ).2 (reluAngle_iterate_antitone hπ hd)
  rw [← antipodalRadius_eq_cos (by decide : 1 ≤ 2), antipodalRadius_two,
    ← antipodalRadius_eq_cos (by omega : 1 ≤ d)] at h
  exact (one_div_pos.mpr pi_pos).trans_le h

theorem antipodalRadius_tendsto_one : Tendsto antipodalRadius atTop (𝓝 1) := by
  have h := continuous_cos.continuousAt.tendsto.comp
    (reluAngle_iterate_tendsto_zero (right_mem_Icc.mpr pi_pos.le))
  rw [cos_zero] at h
  apply h.congr'
  filter_upwards [eventually_ge_atTop 1] with d hd
  exact (antipodalRadius_eq_cos hd).symm

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The paper's unnormalized vector field h, with the two terms parallel to x collected. -/
noncomputable def comparisonField (d : ℕ) (x y : E) : E :=
  let θ := angle x y
  ((2 : ℝ) ^ d)⁻¹ • ((1 - angleSum θ d * (‖y‖ / ‖x‖)) • x - angleProduct θ d • y)

/-- The unit-direction formula in Lemma Seps with its original relative norm and signs. -/
theorem comparisonField_normalized (d : ℕ) {x y : E} (hx : x ≠ 0) (hy : y ≠ 0) :
    ‖y‖⁻¹ • comparisonField d x y = ((2 : ℝ) ^ d)⁻¹ •
      ((‖x‖ / ‖y‖ - angleSum (angle x y) d) • normalize x -
        angleProduct (angle x y) d • normalize y) := by
  have hx' := norm_ne_zero_iff.mpr hx
  have hy' := norm_ne_zero_iff.mpr hy
  rw [comparisonField, smul_comm ‖y‖⁻¹, smul_sub]
  simp only [NormedSpace.normalize, smul_smul]
  congr 2 <;> congr 1 <;> field_simp

theorem comparisonField_zero_iff_balance (d : ℕ) (x y : E) :
    comparisonField d x y = 0 ↔
      (1 - angleSum (angle x y) d * (‖y‖ / ‖x‖)) • x = angleProduct (angle x y) d • y := by
  simp [comparisonField, sub_eq_zero]

theorem comparisonField_zero_iff_of_angle_zero {d : ℕ} {x y : E} (hθ : angle x y = 0) :
    comparisonField d x y = 0 ↔ x = y := by
  simp [comparisonField_zero_iff_balance, hθ, angleProduct_at_zero, angleSum_at_zero]

theorem comparisonField_zero_iff_of_neg_smul {d : ℕ} (hd : 1 ≤ d) {y : E} (hy : y ≠ 0)
    {r : ℝ} (hr : r < 0) : comparisonField d (r • y) y = 0 ↔ r = -antipodalRadius d := by
  have hθ : angle (r • y) y = π := by simp [angle_smul_left_of_neg _ _ hr, hy]
  have hratio : ‖y‖ / ‖r • y‖ = (-r)⁻¹ := by
    rw [norm_smul, Real.norm_eq_abs, abs_of_neg hr]
    field_simp [norm_ne_zero_iff.mpr hy, hr.ne]
  rw [comparisonField_zero_iff_balance, hθ, angleProduct_at_pi hd, zero_smul, hratio,
    smul_eq_zero_iff_left (smul_ne_zero hr.ne hy), ← div_eq_mul_inv, sub_eq_zero,
    eq_div_iff (neg_ne_zero.mpr hr.ne), one_mul]
  change -r = antipodalRadius d ↔ _
  exact neg_eq_iff_eq_neg

theorem comparisonField_ne_zero_of_interior {d : ℕ} {x y : E} (hx : x ≠ 0) (hy : y ≠ 0)
    (hθ : angle x y ∈ Ioo 0 π) : comparisonField d x y ≠ 0 := by
  intro hzero
  let c := 1 - angleSum (angle x y) d * (‖y‖ / ‖x‖)
  have heq : c • x = angleProduct (angle x y) d • y :=
    (comparisonField_zero_iff_balance d x y).mp hzero
  have hp := angleProduct_pos ⟨hθ.1.le, hθ.2⟩ d
  have hne : c ≠ 0 := by
    intro h
    rw [h, zero_smul] at heq
    exact smul_ne_zero hp.ne' hy heq.symm
  have ha : angle x (c • x) = angle x y := by rw [heq, angle_smul_right_of_pos _ _ hp]
  rcases hne.lt_or_gt with hc | hc
  · rw [angle_smul_right_of_neg _ _ hc, angle_self_neg_of_nonzero hx] at ha
    exact hθ.2.ne ha.symm
  · rw [angle_smul_right_of_pos _ _ hc, angle_self hx] at ha
    exact hθ.1.ne ha

/-- For nonzero codes and depth at least two the comparison field has exactly two zeros. -/
theorem comparisonField_zero_iff {d : ℕ} (hd : 2 ≤ d) {x y : E} (hx : x ≠ 0) (hy : y ≠ 0) :
    comparisonField d x y = 0 ↔ x = y ∨ x = (-antipodalRadius d) • y := by
  have hd' : 1 ≤ d := by omega
  constructor
  · intro h
    rcases (angle_nonneg x y).eq_or_lt with hz | hp
    · exact Or.inl ((comparisonField_zero_iff_of_angle_zero hz.symm).mp h)
    rcases (angle_le_pi x y).lt_or_eq with hi | hπ
    · exact (comparisonField_ne_zero_of_interior hx hy ⟨hp, hi⟩ h).elim
    · obtain ⟨_, r, hr, rfl⟩ := angle_eq_pi_iff.mp ((angle_comm y x).trans hπ)
      have hradius := (comparisonField_zero_iff_of_neg_smul hd' hy hr).mp h
      exact Or.inr (congrArg (fun t : ℝ ↦ t • y) hradius)
  · rintro (rfl | rfl)
    · exact (comparisonField_zero_iff_of_angle_zero (angle_self hx)).mpr rfl
    · have hneg := neg_lt_zero.mpr (antipodalRadius_pos hd)
      exact (comparisonField_zero_iff_of_neg_smul hd' hy hneg).mpr rfl

/-- The two distinct nonzero roots for a unit reference on the real line at depth two. -/
theorem comparisonField_two_roots_witness :
    -(1 / π) < (0 : ℝ) ∧ comparisonField 2 (1 : ℝ) 1 = 0 ∧
      comparisonField 2 (-(1 / π) : ℝ) 1 = 0 := by
  have hz := neg_lt_zero.mpr (one_div_pos.mpr pi_pos)
  refine ⟨hz, (comparisonField_zero_iff (by decide) one_ne_zero one_ne_zero).mpr (Or.inl rfl), ?_⟩
  apply (comparisonField_zero_iff (by decide) hz.ne one_ne_zero).mpr
  exact Or.inr (by simp [antipodalRadius_two])

/-- Orthogonal coordinate vectors give a nonzero field in the plane. -/
theorem comparisonField_orthogonal_witness :
    comparisonField 2 (EuclideanSpace.single (0 : Fin 2) (1 : ℝ))
      (EuclideanSpace.single (1 : Fin 2) (1 : ℝ)) ≠ 0 := by
  apply comparisonField_ne_zero_of_interior (by simp [EuclideanSpace.single])
    (by simp [EuclideanSpace.single])
  have hθ := (inner_eq_zero_iff_angle_eq_pi_div_two
    (EuclideanSpace.single (0 : Fin 2) (1 : ℝ)) (EuclideanSpace.single (1 : Fin 2) (1 : ℝ))).mp
    (by norm_num [EuclideanSpace.inner_single_left, EuclideanSpace.single])
  rw [hθ]
  exact ⟨half_pos pi_pos, half_lt_self pi_pos⟩

#print axioms angleSum_succ
#print axioms cos_iterate_eq_angleProduct_add_angleSum
#print axioms reluAngle_iterate_zero
#print axioms angleProduct_at_zero
#print axioms angleSum_at_zero
#print axioms angleProduct_pos
#print axioms angleProduct_at_pi
#print axioms antipodalRadius_eq_cos
#print axioms antipodalRadius_two
#print axioms antipodalRadius_pos
#print axioms antipodalRadius_tendsto_one
#print axioms comparisonField_normalized
#print axioms comparisonField_zero_iff_balance
#print axioms comparisonField_zero_iff_of_angle_zero
#print axioms comparisonField_zero_iff_of_neg_smul
#print axioms comparisonField_ne_zero_of_interior
#print axioms comparisonField_zero_iff
#print axioms comparisonField_two_roots_witness
#print axioms comparisonField_orthogonal_witness

end GenerativePrior
