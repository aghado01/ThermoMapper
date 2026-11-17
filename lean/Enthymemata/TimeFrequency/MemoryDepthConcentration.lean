import Enthymemata.TimeFrequency.MemoryDepthUpper
import Mathlib.Analysis.Complex.Exponential

/-! # Finite concentration tools for history-conditioned transition estimates -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset MeasureTheory ProbabilityTheory

variable {X H : Type*} [Fintype X] [Fintype H]

open scoped Classical in
noncomputable def probability (P : Law X) (E : X → Prop) : ℝ :=
  expect P (fun x ↦ if E x then 1 else 0)

theorem expect_mono (P : Law X) {f g : X → ℝ} (h : ∀ x, f x ≤ g x) :
    expect P f ≤ expect P g :=
  sum_le_sum (fun x _ ↦ mul_le_mul_of_nonneg_left (h x) (P.positive x).le)

theorem probability_mono (P : Law X) {E F : X → Prop} (h : ∀ x, E x → F x) :
    probability P E ≤ probability P F := by
  classical
  apply expect_mono
  intro x
  by_cases he : E x <;> by_cases hf : F x <;> simp_all

theorem probability_union (P : Law X) (E F : X → Prop) :
    probability P (fun x ↦ E x ∨ F x) ≤ probability P E + probability P F := by
  classical
  have h := expect_mono P (f := fun x ↦ if E x ∨ F x then 1 else 0)
    (g := fun x ↦ (if E x then 1 else 0) + (if F x then 1 else 0)) (by
      intro x; by_cases he : E x <;> by_cases hf : F x <;> simp [he, hf])
  simp only [probability, expect, mul_add, sum_add_distrib] at *
  convert h using 1
  apply sum_congr rfl
  intro x _
  split_ifs <;> rfl

theorem probability_exists (P : Law X) (E : H → X → Prop) :
    probability P (fun x ↦ ∃ h, E h x) ≤ ∑ h, probability P (E h) := by
  classical
  have hb (x : X) : (if ∃ h, E h x then (1 : ℝ) else 0) ≤
      ∑ h, if E h x then (1 : ℝ) else 0 := by
    split_ifs with hx
    · obtain ⟨h, hh⟩ := hx
      exact (by simp [hh] : (1 : ℝ) = if E h x then 1 else 0).le.trans
        (single_le_sum (f := fun h ↦ if E h x then (1 : ℝ) else 0)
          (fun _ _ ↦ by split_ifs <;> norm_num) (mem_univ h))
    · exact sum_nonneg (fun _ _ ↦ by split_ifs <;> norm_num)
  have h := expect_mono P hb
  simp only [probability, expect, mul_sum] at *
  rw [sum_comm]
  convert h using 1
  apply sum_congr rfl
  intro x _
  split_ifs <;> rfl

theorem iid_exponential (P : Law X) (f : X → ℝ) (n : ℕ) :
    expect (iid n P) (fun x ↦ Real.exp (∑ i, f (x i))) =
      (expect P (fun x ↦ Real.exp (f x))) ^ n := by
  simp only [expect, iid, Real.exp_sum, ← prod_mul_distrib]
  rw [← Fintype.prod_sum (f := fun (_ : Fin n) x ↦ P.mass x * Real.exp (f x))]
  simp

theorem finite_chernoff (P : Law X) (f : X → ℝ) (n : ℕ) {t b : ℝ} (ht : 0 < t) :
    probability (iid n P) (fun x ↦ b ≤ ∑ i, f (x i)) ≤
      Real.exp (-t * b) * (expect P (fun x ↦ Real.exp (t * f x))) ^ n := by
  classical
  have hp (x : Fin n → X) :
      (if b ≤ ∑ i, f (x i) then (1 : ℝ) else 0) ≤
        Real.exp (-t * b) * Real.exp (∑ i, t * f (x i)) := by
    rw [← Real.exp_add, ← mul_sum]
    split_ifs with hx
    · apply Real.one_le_exp_iff.mpr
      nlinarith
    · exact (Real.exp_pos _).le
  have h := expect_mono (iid n P) hp
  have he : expect (iid n P)
      (fun x ↦ Real.exp (-t * b) * Real.exp (∑ i, t * f (x i))) =
      Real.exp (-t * b) * expect (iid n P) (fun x ↦ Real.exp (∑ i, t * f (x i))) := by
    simp only [expect, mul_sum]; apply sum_congr rfl; intros; ring
  rw [he, iid_exponential P (fun x ↦ t * f x) n] at h
  exact h

theorem bounded_mgf [MeasurableSpace X] [MeasurableSingletonClass X]
    (P : Law X) (f : X → ℝ) (hf : ∀ x, f x ∈ Set.Icc (0 : ℝ) 1) (t : ℝ) :
    expect P (fun x ↦ Real.exp (t * (f x - expect P f))) ≤ Real.exp (t ^ 2 / 8) := by
  have h := hasSubgaussianMGF_of_mem_Icc (μ := P.measure)
    (measurable_of_finite f).aemeasurable (Filter.Eventually.of_forall hf)
  have hm := h.mgf_le t
  simp only [mgf, integral_expect] at hm
  norm_num at hm
  convert hm using 1
  congr 1
  ring

theorem exp_decrement {t : ℝ} (ht : 0 ≤ t) (hu : t ≤ 1 / 2) :
    Real.exp (-2 * t ^ 2) ≤ 1 - t ^ 2 := by
  have h := Real.add_one_le_exp (2 * t ^ 2)
  have hpos := Real.exp_pos (2 * t ^ 2)
  have he : Real.exp (-2 * t ^ 2) = 1 / Real.exp (2 * t ^ 2) := by
    rw [show -2 * t ^ 2 = -(2 * t ^ 2) by ring, Real.exp_neg, one_div]
  rw [he]
  apply (div_le_iff₀ hpos).mpr
  have hsq : t ^ 2 ≤ 1 / 4 := by nlinarith
  have hm := mul_le_mul_of_nonneg_left h (show 0 ≤ 1 - t ^ 2 by linarith)
  nlinarith [sq_nonneg t]

theorem exp_neg_one_bound : Real.exp (-1 : ℝ) ≤ 3 / 8 := by
  have h := Real.sum_le_exp_of_nonneg (show (0 : ℝ) ≤ 1 by norm_num) 4
  norm_num [sum_range_succ] at h
  rw [Real.exp_neg, ← one_div]
  apply (div_le_iff₀ (Real.exp_pos 1)).mpr
  linarith

#print axioms expect_mono
#print axioms probability_mono
#print axioms probability_union
#print axioms probability_exists
#print axioms iid_exponential
#print axioms finite_chernoff
#print axioms bounded_mgf
#print axioms exp_decrement
#print axioms exp_neg_one_bound

end TimeFrequency.ContextTree.DepthEvidence
