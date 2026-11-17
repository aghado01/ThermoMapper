import Enthymemata.TimeFrequency.PredictiveActivityProjection

/-! # Three-head population excess controls the encoded joint activity law -/

namespace TimeFrequency.PredictiveActivity

open Finset hiding expect
open ContextTree.DepthEvidence

variable {X Z : Type*} [Fintype X] [Finite Z] [DecidableEq Z]

omit [Finite Z] in
theorem conditional_add (P : Law X) (f : X → Z) (a b : X → ℝ) (z : Z) :
    conditionalMean P f (fun x ↦ a x + b x) z =
      conditionalMean P f a z + conditionalMean P f b z := by
  unfold conditionalMean fiberSum
  rw [← add_div, ← sum_add_distrib]
  congr 1
  apply sum_congr rfl
  intro x _
  split_ifs <;> ring

noncomputable def encodedLaw (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (p : X → ActivityLaw) (z : Z) : ActivityLaw where
  prob y := conditionalMean P f (fun x ↦ (p x).prob y) z
  nonneg y := div_nonneg
    (sum_nonneg (fun x _ ↦ by
      split_ifs
      · exact mul_nonneg (P.positive x).le ((p x).nonneg y)
      · exact le_rfl)) (code_mass_pos P f hf z).le
  total := by
    simp only [conditionalMean, ← sum_div, fiberSum]
    rw [sum_comm]
    have h (x : X) :
        (∑ y, if f x = z then P.mass x * (p x).prob y else 0) =
          if f x = z then P.mass x else 0 := by
      split_ifs with hx
      · simp only [← mul_sum, (p x).total, mul_one]
      · simp
    simp only [h]
    have hm : (∑ x, if f x = z then P.mass x else 0) = codeMass P f z := by
      simp [codeMass, fiberSum]
    rw [hm, div_self (ne_of_gt (code_mass_pos P f hf z))]

omit [Finite Z] in
theorem encoded_target_means (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (p : X → ActivityLaw) (z : Z) (j : Fin 3) :
    targetMeans (encodedLaw P f hf p z) j =
      conditionalMean P f (fun x ↦ targetMeans (p x) j) z := by
  fin_cases j <;> simp [targetMeans, encodedLaw, conditional_add]

noncomputable def headExcess (P : Law X) (f : X → Z) (p : X → ActivityLaw)
    (head : Z → Fin 3 → ℝ) : ℝ :=
  ∑ j, expect P (fun x ↦ binaryRisk (targetMeans (p x) j) (head (f x) j) -
    binaryRisk (targetMeans (p x) j) (targetMeans (p x) j))

omit [Finite Z] [DecidableEq Z] in
theorem head_excess_squares (P : Law X) (f : X → Z) (p : X → ActivityLaw)
    (head : Z → Fin 3 → ℝ) :
    headExcess P f p head =
      ∑ j, expect P (fun x ↦ (targetMeans (p x) j - head (f x) j) ^ 2) := by
  simp only [headExcess, binary_excess, sub_sq_comm]

theorem joint_retention_bound (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (p : X → ActivityLaw) (head : Z → Fin 3 → ℝ) :
    expect P (fun x ↦ ∑ y, ((p x).prob y - (encodedLaw P f hf p (f x)).prob y) ^ 2)
      ≤ 8 * headExcess P f p head := by
  calc
    _ ≤ expect P (fun x ↦ 8 * ∑ j,
        (targetMeans (p x) j - targetMeans (encodedLaw P f hf p (f x)) j) ^ 2) := by
      exact sum_le_sum (fun x _ ↦ mul_le_mul_of_nonneg_left
        (activity_squared_bound _ _) (P.positive x).le)
    _ = 8 * ∑ j, expect P (fun x ↦
        (targetMeans (p x) j - conditionalMean P f (fun x ↦ targetMeans (p x) j)
          (f x)) ^ 2) := by
      simp only [encoded_target_means, expect, mul_sum]
      rw [sum_comm]
      apply sum_congr rfl
      intro j _
      apply sum_congr rfl
      intro x _
      ring
    _ ≤ 8 * headExcess P f p head := by
      rw [head_excess_squares]
      exact mul_le_mul_of_nonneg_left (sum_le_sum (fun j _ ↦
        projection_error_le P f hf (fun x ↦ targetMeans (p x) j)
          (fun z ↦ head z j))) (by norm_num)

theorem joint_retention_exact (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (p : X → ActivityLaw) (head : Z → Fin 3 → ℝ)
    (he : headExcess P f p head = 0) :
    ∀ x, (p x).prob = (encodedLaw P f hf p (f x)).prob := by
  have hb := joint_retention_bound P f hf p head
  rw [he, mul_zero] at hb
  intro x
  have hn (x : X) : 0 ≤ ∑ y,
      ((p x).prob y - (encodedLaw P f hf p (f x)).prob y) ^ 2 :=
    sum_nonneg (fun _ _ ↦ sq_nonneg _)
  have hx := (single_le_sum (fun a _ ↦ mul_nonneg (P.positive a).le (hn a))
    (mem_univ x)).trans hb
  have hz : (∑ y, ((p x).prob y - (encodedLaw P f hf p (f x)).prob y) ^ 2) = 0 :=
    le_antisymm (le_of_not_gt (fun hs ↦ (not_lt_of_ge hx) (mul_pos (P.positive x) hs))) (hn x)
  funext y
  have hy := (single_le_sum (fun _ _ ↦ sq_nonneg _) (mem_univ y)).trans_eq hz
  exact sub_eq_zero.mp (sq_eq_zero_iff.mp (le_antisymm hy (sq_nonneg _)))

#print axioms conditional_add
#print axioms encoded_target_means
#print axioms head_excess_squares
#print axioms joint_retention_bound
#print axioms joint_retention_exact

end TimeFrequency.PredictiveActivity
