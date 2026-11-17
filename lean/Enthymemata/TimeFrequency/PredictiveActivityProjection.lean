import Enthymemata.TimeFrequency.PredictiveActivityCore

/-! # Conditional projection on a frozen finite encoder -/

namespace TimeFrequency.PredictiveActivity

open Finset hiding expect
open ContextTree.DepthEvidence

variable {X Z : Type*} [Fintype X] [DecidableEq Z]

noncomputable def fiberSum (P : Law X) (f : X → Z) (a : X → ℝ) (z : Z) : ℝ :=
  ∑ x, if f x = z then P.mass x * a x else 0

noncomputable def codeMass (P : Law X) (f : X → Z) (z : Z) : ℝ :=
  fiberSum P f (fun _ ↦ 1) z

noncomputable def conditionalMean (P : Law X) (f : X → Z) (a : X → ℝ) (z : Z) : ℝ :=
  fiberSum P f a z / codeMass P f z

theorem code_mass_pos (P : Law X) (f : X → Z) (hf : Function.Surjective f) (z : Z) :
    0 < codeMass P f z := by
  obtain ⟨x, hx⟩ := hf z
  unfold codeMass fiberSum
  have hn (y : X) : 0 ≤ if f y = z then P.mass y * 1 else 0 := by
    split_ifs
    · exact mul_nonneg (P.positive y).le (by norm_num)
    · exact le_rfl
  apply sum_pos' (fun y _ ↦ hn y) ⟨x, mem_univ _, ?_⟩
  simpa only [hx, if_true, mul_one] using P.positive x

theorem fiber_encoded (P : Law X) (f : X → Z) (a : Z → ℝ) (z : Z) :
    fiberSum P f (fun x ↦ a (f x)) z = codeMass P f z * a z := by
  unfold codeMass fiberSum
  rw [sum_mul]
  apply sum_congr rfl
  intro x _
  by_cases h : f x = z <;> simp [h]

theorem fiber_sub (P : Law X) (f : X → Z) (a b : X → ℝ) (z : Z) :
    fiberSum P f (fun x ↦ a x - b x) z =
      fiberSum P f a z - fiberSum P f b z := by
  unfold fiberSum
  rw [← sum_sub_distrib]
  apply sum_congr rfl
  intro x _
  split_ifs <;> ring

theorem regroup [Fintype Z] (P : Law X) (f : X → Z) (a : X → ℝ) (b : Z → ℝ) :
    expect P (fun x ↦ a x * b (f x)) = ∑ z, fiberSum P f a z * b z := by
  classical
  simp only [fiberSum, sum_mul, expect]
  rw [sum_comm]
  apply sum_congr rfl
  intro x _
  simp only [ite_mul, zero_mul]
  simp [mul_assoc]

theorem conditional_orthogonal [Finite Z] (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (a : X → ℝ) (b : Z → ℝ) :
    expect P (fun x ↦ (a x - conditionalMean P f a (f x)) * b (f x)) = 0 := by
  letI := Fintype.ofFinite Z
  rw [regroup]
  apply sum_eq_zero
  intro z _
  rw [fiber_sub, fiber_encoded]
  unfold conditionalMean
  field_simp [ne_of_gt (code_mass_pos P f hf z)]
  ring

theorem projection_identity [Finite Z] (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (a : X → ℝ) (head : Z → ℝ) :
    expect P (fun x ↦ (a x - head (f x)) ^ 2) =
      expect P (fun x ↦ (a x - conditionalMean P f a (f x)) ^ 2) +
      expect P (fun x ↦ (conditionalMean P f a (f x) - head (f x)) ^ 2) := by
  have ho := conditional_orthogonal P f hf a
    (fun z ↦ conditionalMean P f a z - head z)
  have he :
      expect P (fun x ↦ (a x - head (f x)) ^ 2) =
        expect P (fun x ↦ (a x - conditionalMean P f a (f x)) ^ 2) +
        expect P (fun x ↦ (conditionalMean P f a (f x) - head (f x)) ^ 2) +
        2 * expect P (fun x ↦ (a x - conditionalMean P f a (f x)) *
          (conditionalMean P f a (f x) - head (f x))) := by
    simp only [expect, mul_sum, ← sum_add_distrib]
    apply sum_congr rfl
    intro x _
    ring
  rwa [ho, mul_zero, add_zero] at he

theorem projection_error_le [Finite Z] (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (a : X → ℝ) (head : Z → ℝ) :
    expect P (fun x ↦ (a x - conditionalMean P f a (f x)) ^ 2) ≤
      expect P (fun x ↦ (a x - head (f x)) ^ 2) := by
  rw [projection_identity P f hf a head]
  apply le_add_of_nonneg_right
  exact sum_nonneg (fun x _ ↦ mul_nonneg (P.positive x).le (sq_nonneg _))

theorem encoded_expect [Finite Z] (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (a : X → ℝ) (b : Z → ℝ) :
    expect P (fun x ↦ a x * b (f x)) =
      expect P (fun x ↦ conditionalMean P f a (f x) * b (f x)) := by
  have ho := conditional_orthogonal P f hf a b
  have he : expect P (fun x ↦ (a x - conditionalMean P f a (f x)) * b (f x)) =
      expect P (fun x ↦ a x * b (f x)) -
      expect P (fun x ↦ conditionalMean P f a (f x) * b (f x)) := by
    simp only [expect, ← sum_sub_distrib]
    apply sum_congr rfl
    intro x _
    ring
  exact sub_eq_zero.mp (he.symm.trans ho)

#print axioms encoded_expect
#print axioms code_mass_pos
#print axioms fiber_encoded
#print axioms fiber_sub
#print axioms regroup
#print axioms conditional_orthogonal
#print axioms projection_identity
#print axioms projection_error_le

end TimeFrequency.PredictiveActivity
