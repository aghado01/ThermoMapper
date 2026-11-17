import Enthymemata.TimeFrequency.PredictiveActivityProjection

/-! # Strictly positive predictive information from a lost conditional distinction -/

namespace TimeFrequency.PredictiveActivity

open Finset hiding expect
open ContextTree.DepthEvidence

theorem scalar_kl_lower {p q : ℝ} (hp : 0 < p) (hq : 0 < q) :
    p - q ≤ p * Real.log (p / q) := by
  have h := mul_le_mul_of_nonneg_left
    (Real.log_le_sub_one_of_pos (div_pos hq hp)) hp.le
  have he : Real.log (q / p) = -Real.log (p / q) := by
    rw [Real.log_div hq.ne' hp.ne', Real.log_div hp.ne' hq.ne']
    ring
  have hc : p * (q / p - 1) = q - p := by field_simp [hp.ne']
  rw [he, hc] at h
  linarith

theorem scalar_kl_strict {p q : ℝ} (hp : 0 < p) (hq : 0 < q) (hne : p ≠ q) :
    p - q < p * Real.log (p / q) := by
  have hn : q / p ≠ 1 := by
    intro he
    exact hne ((div_eq_one_iff_eq hp.ne').mp he).symm
  have h := mul_lt_mul_of_pos_left
    (Real.log_lt_sub_one_of_pos (div_pos hq hp) hn) hp
  have he : Real.log (q / p) = -Real.log (p / q) := by
    rw [Real.log_div hq.ne' hp.ne', Real.log_div hp.ne' hq.ne']
    ring
  have hc : p * (q / p - 1) = q - p := by field_simp [hp.ne']
  rw [he, hc] at h
  linarith

theorem divergence_nonnegative {Y : Type*} [Fintype Y] (p q : Law Y) :
    0 ≤ divergence p q := by
  have h := sum_le_sum (fun y (_ : y ∈ univ) ↦ scalar_kl_lower (p.positive y) (q.positive y))
  simpa only [sum_sub_distrib, p.total, q.total, sub_self, divergence, expect,
    likelihood] using h

theorem divergence_positive {Y : Type*} [Fintype Y] (p q : Law Y)
    (hne : ∃ y, p.mass y ≠ q.mass y) : 0 < divergence p q := by
  obtain ⟨y, hy⟩ := hne
  have h := sum_lt_sum
    (fun a (_ : a ∈ univ) ↦ scalar_kl_lower (p.positive a) (q.positive a))
    ⟨y, mem_univ _, scalar_kl_strict (p.positive y) (q.positive y) hy⟩
  simpa only [sum_sub_distrib, p.total, q.total, sub_self, divergence, expect,
    likelihood] using h

variable {X Z : Type*} [Fintype X] [DecidableEq Z]

theorem conditional_range (P : Law X) (f : X → Z) (hf : Function.Surjective f)
    (a : X → ℝ) (ha : ∀ x, a x ∈ Set.Ioo 0 1) (z : Z) :
    conditionalMean P f a z ∈ Set.Ioo 0 1 := by
  obtain ⟨x, hx⟩ := hf z
  have hl : 0 < fiberSum P f a z := by
    have hn (y : X) : 0 ≤ if f y = z then P.mass y * a y else 0 := by
      split_ifs
      · exact (mul_pos (P.positive y) (ha y).1).le
      · exact le_rfl
    apply sum_pos' (fun y _ ↦ hn y) ⟨x, mem_univ _, ?_⟩
    simp only [hx, if_true]
    exact mul_pos (P.positive x) (ha x).1
  have hu : fiberSum P f a z < codeMass P f z := by
    apply sum_lt_sum
    · intro y _
      by_cases hy : f y = z
      · simp only [hy, if_true, mul_one]
        exact mul_le_of_le_one_right (P.positive y).le (ha y).2.le
      · simp [hy]
    · refine ⟨x, mem_univ _, ?_⟩
      simp only [hx, if_true, mul_one]
      exact mul_lt_of_lt_one_right (P.positive x) (ha x).2
  exact ⟨div_pos hl (code_mass_pos P f hf z),
    (div_lt_one (code_mass_pos P f hf z)).mpr hu⟩

noncomputable def binaryInformation (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (a : X → ℝ) (ha : ∀ x, a x ∈ Set.Ioo 0 1) : ℝ :=
  expect P (fun x ↦ divergence (bernoulli (a x) (ha x))
    (bernoulli (conditionalMean P f a (f x)) (conditional_range P f hf a ha (f x)))) /
    Real.log 2

theorem information_positive (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (a : X → ℝ) (ha : ∀ x, a x ∈ Set.Ioo 0 1)
    (hd : ∃ x x', f x = f x' ∧ a x ≠ a x') :
    0 < binaryInformation P f hf a ha := by
  have hw : ∃ x, a x ≠ conditionalMean P f a (f x) := by
    by_contra h
    push Not at h
    obtain ⟨x, x', he, hn⟩ := hd
    exact hn ((h x).trans (by rw [he, ← h x']))
  obtain ⟨x, hx⟩ := hw
  apply div_pos ?_ (Real.log_pos (by norm_num))
  apply sum_pos' (fun y _ ↦ mul_nonneg (P.positive y).le (divergence_nonnegative _ _))
  refine ⟨x, mem_univ _, mul_pos (P.positive x) (divergence_positive _ _ ?_)⟩
  exact ⟨true, by simpa [bernoulli, bitMass] using hx⟩

theorem conditional_constant (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (c : ℝ) (z : Z) :
    conditionalMean P f (fun _ ↦ c) z = c := by
  unfold conditionalMean
  have he := fiber_encoded P f (fun _ : Z ↦ c) z
  change fiberSum P f (fun _ ↦ c) z = codeMass P f z * c at he
  rw [he]
  exact mul_div_cancel_left₀ c (code_mass_pos P f hf z).ne'

theorem information_constant (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) {c : ℝ} (hc : c ∈ Set.Ioo 0 1) :
    binaryInformation P f hf (fun _ ↦ c) (fun _ ↦ hc) = 0 := by
  simp only [binaryInformation, conditional_constant P f hf]
  simp [divergence, expect, likelihood, bernoulli, bitMass]

#print axioms scalar_kl_lower
#print axioms scalar_kl_strict
#print axioms divergence_nonnegative
#print axioms divergence_positive
#print axioms conditional_range
#print axioms information_positive
#print axioms conditional_constant
#print axioms information_constant

end TimeFrequency.PredictiveActivity
