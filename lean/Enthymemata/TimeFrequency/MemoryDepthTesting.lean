import Enthymemata.TimeFrequency.MemoryDepthLikelihood

/-! # Necessary taper-bit budgets, including randomized decisions -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset

variable {X Y : Type*} [Fintype X] [Fintype Y]

theorem log_tangent {p q z : ℝ} (hp : 0 < p) (hq : 0 < q) (hz : 0 < z) :
    p * Real.log z + p - z * q ≤ p * Real.log (p / q) := by
  have h := Real.log_le_sub_one_of_pos (div_pos (mul_pos hz hq) hp)
  have he : Real.log (z * q / p) = Real.log z - Real.log (p / q) := by
    rw [Real.log_div (mul_ne_zero hz.ne' hq.ne') hp.ne',
      Real.log_mul hz.ne' hq.ne', Real.log_div hp.ne' hq.ne']
    ring
  rw [he] at h
  have hm := mul_le_mul_of_nonneg_left h hp.le
  have hf : p * (z * q / p - 1) = z * q - p := by field_simp
  rw [hf] at hm
  linarith

theorem decision_tangent (P Q : Law X) (t : X → ℝ) (ht : ∀ x, t x ∈ Set.Icc 0 1)
    {z w : ℝ} (hz : 0 < z) (hw : 0 < w) :
    expect P t * Real.log z + (1 - expect P t) * Real.log w + 1 -
      z * expect Q t - w * (1 - expect Q t) ≤ divergence P Q := by
  have hc (R : Law X) : (∑ x, R.mass x * (1 - t x)) = 1 - expect R t := by
    simp [expect, mul_sub, sum_sub_distrib, R.total]
  have hp (x : X) :
      (P.mass x * t x) * Real.log z + (P.mass x * (1 - t x)) * Real.log w +
        P.mass x - z * (Q.mass x * t x) - w * (Q.mass x * (1 - t x)) ≤
        P.mass x * Real.log (P.mass x / Q.mass x) := by
    have h1 := mul_le_mul_of_nonneg_left
      (log_tangent (P.positive x) (Q.positive x) hz) (ht x).1
    have h0 := mul_le_mul_of_nonneg_left
      (log_tangent (P.positive x) (Q.positive x) hw) (sub_nonneg.mpr (ht x).2)
    nlinarith only [h1, h0]
  have hs := sum_le_sum (fun x (_ : x ∈ univ) ↦ hp x)
  simpa only [sum_sub_distrib, sum_add_distrib, ← sum_mul, ← mul_sum,
    P.total, hc, divergence, likelihood, expect] using hs

theorem binary_data_processing (P Q : Law X) (t : X → ℝ)
    (ht : ∀ x, t x ∈ Set.Icc 0 1)
    (hp : expect P t ∈ Set.Ioo 0 1) (hq : expect Q t ∈ Set.Ioo 0 1) :
    kl (expect P t) (expect Q t) ≤ divergence P Q := by
  have h := decision_tangent P Q t ht (div_pos hp.1 hq.1)
    (div_pos (sub_pos.mpr hp.2) (sub_pos.mpr hq.2))
  have h1 : (expect P t / expect Q t) * expect Q t = expect P t :=
    div_mul_cancel₀ _ hq.1.ne'
  have h0 : ((1 - expect P t) / (1 - expect Q t)) * (1 - expect Q t) =
      1 - expect P t := div_mul_cancel₀ _ (sub_pos.mpr hq.2).ne'
  rw [h1, h0] at h
  unfold kl
  linarith

theorem testing_divergence (P Q : Law X) (t : X → ℝ)
    (ht : ∀ x, t x ∈ Set.Icc 0 1) {δ : ℝ} (hδ : 0 < δ) (hu : δ < 1 / 2)
    (hP : 1 - δ ≤ expect P t) (hQ : expect Q t ≤ δ) :
    kl (1 - δ) δ ≤ divergence P Q := by
  have hδ1 : δ < 1 := by linarith
  let z := (1 - δ) / δ
  let w := δ / (1 - δ)
  have hz : 0 < z := div_pos (sub_pos.mpr hδ1) hδ
  have hw : 0 < w := div_pos hδ (sub_pos.mpr hδ1)
  have hz1 : 1 ≤ z := (le_div_iff₀ hδ).mpr (by linarith)
  have hw1 : w ≤ 1 := (div_le_iff₀ (sub_pos.mpr hδ1)).mpr (by linarith)
  have hlz : 0 ≤ Real.log z := Real.log_nonneg hz1
  have hlw : Real.log w = -Real.log z := by
    dsimp [z, w]
    rw [Real.log_div hδ.ne' (sub_pos.mpr hδ1).ne',
      Real.log_div (sub_pos.mpr hδ1).ne' hδ.ne']
    ring
  have hez : z * δ = 1 - δ := div_mul_cancel₀ _ hδ.ne'
  have hew : w * (1 - δ) = δ := div_mul_cancel₀ _ (sub_pos.mpr hδ1).ne'
  have h := decision_tangent P Q t ht hz hw
  have h1 := mul_nonneg (sub_nonneg.mpr hP)
    (show 0 ≤ Real.log z - Real.log w by rw [hlw]; linarith)
  have h0 := mul_nonneg (sub_nonneg.mpr hQ) (show 0 ≤ z - w by linarith)
  have he : kl (1 - δ) δ = (1 - δ) * Real.log z + δ * Real.log w := by
    simp [kl, z, w]
  rw [he]
  nlinarith only [h, h1, h0, hez, hew]

theorem bernoulli_quadratic_upper {p q : ℝ}
    (hp : p ∈ Set.Ioo 0 1) (hq : q ∈ Set.Ioo 0 1) :
    kl p q ≤ (p - q) ^ 2 / (q * (1 - q)) := by
  have h1 := mul_le_mul_of_nonneg_left
    (Real.log_le_sub_one_of_pos (div_pos hp.1 hq.1)) hp.1.le
  have h0 := mul_le_mul_of_nonneg_left
    (Real.log_le_sub_one_of_pos (div_pos (sub_pos.mpr hp.2) (sub_pos.mpr hq.2)))
    (sub_nonneg.mpr hp.2.le)
  have he : p * (p / q - 1) + (1 - p) * ((1 - p) / (1 - q) - 1) =
      (p - q) ^ 2 / (q * (1 - q)) := by
    field_simp [hq.1.ne', (sub_pos.mpr hq.2).ne']
    ring
  unfold kl
  linarith

theorem taper_divergence_upper {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2) :
    kl ((1 + a) / 8) (1 / 8) ≤ a ^ 2 / 7 := by
  have h := bernoulli_quadratic_upper (taper_parameters ha hu).2 (taper_parameters ha hu).1
  have he : (((1 + a) / 8 - 1 / 8) ^ 2) / ((1 / 8 : ℝ) * (1 - 1 / 8)) =
      a ^ 2 / 7 := by ring
  rwa [he] at h

theorem necessary_budget {a δ : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (hδ : 0 < δ) (hδu : δ < 1 / 2) (n : ℕ) (R : Law Y)
    (t : Y × Word n → ℝ) (ht : ∀ x, t x ∈ Set.Icc 0 1)
    (hP : 1 - δ ≤ expect
      (product R (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2))) t)
    (hQ : expect
      (product R (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1))) t ≤ δ) :
    7 * kl (1 - δ) δ / a ^ 2 ≤ n := by
  have h := testing_divergence _ _ t ht hδ hδu hP hQ
  rw [expected_log_evidence ha hu] at h
  have hb := mul_le_mul_of_nonneg_left (taper_divergence_upper ha hu) (Nat.cast_nonneg n)
  apply (div_le_iff₀ (sq_pos_of_pos ha)).mpr
  nlinarith

theorem necessary_budget_rounded {a δ : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (hδ : 0 < δ) (hδu : δ < 1 / 2) (n : ℕ) (R : Law Y)
    (t : Y × Word n → ℝ) (ht : ∀ x, t x ∈ Set.Icc 0 1)
    (hP : 1 - δ ≤ expect
      (product R (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2))) t)
    (hQ : expect
      (product R (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1))) t ≤ δ) :
    ⌈7 * kl (1 - δ) δ / a ^ 2⌉₊ ≤ n := by
  exact Nat.ceil_le.mpr (necessary_budget ha hu hδ hδu n R t ht hP hQ)

theorem decision_divergence_formula {δ : ℝ} (hδ : 0 < δ) (hu : δ < 1 / 2) :
    kl (1 - δ) δ = (1 - 2 * δ) * Real.log ((1 - δ) / δ) := by
  have hd : 0 < 1 - δ := by linarith
  have he : Real.log (δ / (1 - δ)) = -Real.log ((1 - δ) / δ) := by
    rw [Real.log_div hδ.ne' hd.ne', Real.log_div hd.ne' hδ.ne']
    ring
  simp only [kl, sub_sub_cancel, he]
  ring

theorem decision_divergence_pos {δ : ℝ} (hδ : 0 < δ) (hu : δ < 1 / 2) :
    0 < kl (1 - δ) δ := by
  rw [decision_divergence_formula hδ hu]
  apply mul_pos (by linarith)
  apply Real.log_pos
  exact (lt_div_iff₀ hδ).mpr (by linarith)

#print axioms log_tangent
#print axioms decision_tangent
#print axioms binary_data_processing
#print axioms testing_divergence
#print axioms bernoulli_quadratic_upper
#print axioms taper_divergence_upper
#print axioms necessary_budget
#print axioms necessary_budget_rounded
#print axioms decision_divergence_formula
#print axioms decision_divergence_pos

end TimeFrequency.ContextTree.DepthEvidence
