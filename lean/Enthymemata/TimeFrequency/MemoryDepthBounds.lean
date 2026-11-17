import Enthymemata.TimeFrequency.MemoryDepthCounts

/-! # Uniform finite-sample transition bounds with rare-history dependence -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset

variable {H : Type*} [Fintype H] [DecidableEq H]

theorem ratio_tail (P : Law (H × Bool)) (h : H) (n : ℕ) {t s : ℝ}
    (ht : 0 < t) (hu : t ≤ 1 / 2) (hs : s ^ 2 = 1) :
    probability (iid n P)
        (fun x ↦ t * visits x h ≤ s * (positives x h - conditional P h * visits x h)) ≤
      Real.exp (-(n : ℝ) * historyMass P h * t ^ 2) := by
  let f := fun z : H × Bool ↦
    if z.1 = h then s * (bitValue z.2 - conditional P h) - t else 0
  have hc := finite_chernoff P f n (b := 0) (show 0 < 4 * t by positivity)
  have he (x : Fin n → H × Bool) : (0 ≤ ∑ i, f (x i)) ↔
      t * visits x h ≤ s * (positives x h - conditional P h * visits x h) := by
    dsimp [f]
    rw [masked_sum]
    constructor <;> intro hx <;> linarith
  simp_rw [he] at hc
  simp only [mul_zero, Real.exp_zero, one_mul] at hc
  have hm := ratio_mgf P h ht.le hu hs
  have hn : 0 ≤ expect P (fun z ↦ Real.exp (4 * t * f z)) :=
    sum_nonneg (fun _ _ ↦ mul_nonneg (P.positive _).le (Real.exp_pos _).le)
  have hp := pow_le_pow_left₀ hn hm n
  calc _ ≤ (expect P (fun z ↦ Real.exp (4 * t * f z))) ^ n := hc
    _ ≤ (Real.exp (-historyMass P h * t ^ 2)) ^ n := hp
    _ = _ := by rw [← Real.exp_nat_mul]; congr 1; ring

theorem transition_error (P : Law (H × Bool)) (h : H) (n : ℕ) {t : ℝ}
    (ht : 0 < t) (hu : t ≤ 1 / 2) :
    probability (iid n P) (fun x ↦ t ≤ |estimate x h - conditional P h|) ≤
      2 * Real.exp (-(n : ℝ) * historyMass P h * t ^ 2) := by
  let E := fun x : Fin n → H × Bool ↦
    t * visits x h ≤ positives x h - conditional P h * visits x h
  let F := fun x : Fin n → H × Bool ↦
    t * visits x h ≤ -(positives x h - conditional P h * visits x h)
  have hsub : ∀ x, t ≤ |estimate x h - conditional P h| → E x ∨ F x := by
    intro x hx
    by_cases hn : visits x h = 0
    · left
      simp only [E, hn, mul_zero, sub_zero]
      exact positives_nonneg x h
    · have hnpos := lt_of_le_of_ne (visits_nonneg x h) (Ne.symm hn)
      rw [estimate, if_neg hn] at hx
      rcases (le_abs).mp hx with hp | hp
      · left
        dsimp [E]
        have hm := mul_le_mul_of_nonneg_right hp hnpos.le
        have he : (positives x h / visits x h - conditional P h) * visits x h =
            positives x h - conditional P h * visits x h := by field_simp [hn]
        rwa [he] at hm
      · right
        dsimp [F]
        have hm := mul_le_mul_of_nonneg_right hp hnpos.le
        have he : (-(positives x h / visits x h - conditional P h)) * visits x h =
            -(positives x h - conditional P h * visits x h) := by field_simp [hn]
        rwa [he] at hm
  have hU := ratio_tail P h n ht hu (show (1 : ℝ) ^ 2 = 1 by norm_num)
  have hL := ratio_tail P h n ht hu (show (-1 : ℝ) ^ 2 = 1 by norm_num)
  simp only [one_mul] at hU
  simp only [neg_one_mul] at hL
  exact (probability_mono _ hsub).trans ((probability_union _ E F).trans (by
    dsimp [E, F]
    linarith))

theorem simultaneous_transition_error (P : Law (H × Bool)) (n : ℕ) {t π : ℝ}
    (ht : 0 < t) (hu : t ≤ 1 / 2) (hπ : ∀ h, π ≤ historyMass P h) :
    probability (iid n P) (fun x ↦ ∃ h, t ≤ |estimate x h - conditional P h|) ≤
      2 * Fintype.card H * Real.exp (-(n : ℝ) * π * t ^ 2) := by
  have hbound (h : H) := transition_error P h n ht hu
  have hb (h : H) : Real.exp (-(n : ℝ) * historyMass P h * t ^ 2) ≤
      Real.exp (-(n : ℝ) * π * t ^ 2) := by
    apply Real.exp_le_exp.mpr
    have hm := mul_le_mul_of_nonneg_left (hπ h)
      (show 0 ≤ (n : ℝ) * t ^ 2 by positivity)
    nlinarith
  have hs := sum_le_sum (fun h (_ : h ∈ univ) ↦
    (hbound h).trans (mul_le_mul_of_nonneg_left (hb h) (by norm_num : (0 : ℝ) ≤ 2)))
  have hp := probability_exists (iid n P)
    (fun h x ↦ t ≤ |estimate x h - conditional P h|)
  simp only [sum_const, card_univ, nsmul_eq_mul] at hs
  exact hp.trans (by nlinarith only [hs])

theorem transition_budget (P : Law (H × Bool)) (n : ℕ) {t π δ : ℝ}
    (ht : 0 < t) (hu : t ≤ 1 / 2) (hπ : 0 < π) (hm : ∀ h, π ≤ historyMass P h)
    (hδ : 0 < δ) (hH : 0 < Fintype.card H)
    (hn : ⌈Real.log (2 * Fintype.card H / δ) / (π * t ^ 2)⌉₊ ≤ n) :
    probability (iid n P) (fun x ↦ ∃ h, t ≤ |estimate x h - conditional P h|) ≤ δ := by
  have hr := (Nat.le_ceil (Real.log (2 * Fintype.card H / δ) / (π * t ^ 2))).trans
    (Nat.cast_le.mpr hn)
  have hmul := (div_le_iff₀ (mul_pos hπ (sq_pos_of_pos ht))).mp hr
  have hc : (0 : ℝ) < 2 * Fintype.card H := by exact_mod_cast (by omega : 0 < 2 * Fintype.card H)
  have he : Real.exp (-(n : ℝ) * π * t ^ 2) ≤ δ / (2 * Fintype.card H) := by
    apply (Real.le_log_iff_exp_le (div_pos hδ hc)).mp
    have hlog : Real.log (δ / (2 * Fintype.card H)) =
        -Real.log (2 * Fintype.card H / δ) := by
      rw [Real.log_div hδ.ne' hc.ne', Real.log_div hc.ne' hδ.ne']
      ring
    rw [hlog]
    nlinarith only [hmul]
  have hb := simultaneous_transition_error P n ht hu hm
  have hfinal := mul_le_mul_of_nonneg_left he hc.le
  have heq : (2 * Fintype.card H : ℝ) * (δ / (2 * Fintype.card H)) = δ := by
    field_simp [hc.ne']
  rw [heq] at hfinal
  exact hb.trans hfinal

theorem too_few_visits (P : Law (H × Bool)) (h : H) (n : ℕ) :
    probability (iid n P) (fun x ↦ visits x h ≤ n * historyMass P h / 2) ≤
      Real.exp (-(n : ℝ) * historyMass P h / 8) := by
  let f := fun z : H × Bool ↦ if z.1 = h then (-1 : ℝ) else 0
  have hc := finite_chernoff P f n (t := 1)
    (b := -(n : ℝ) * historyMass P h / 2) (by norm_num)
  have he (x : Fin n → H × Bool) : (∑ i, f (x i)) = -visits x h := by
    simp only [visits, ← sum_neg_distrib]
    apply sum_congr rfl
    intro i _
    dsimp [f]
    split_ifs <;> ring
  simp_rw [he] at hc
  have hevent (x : Fin n → H × Bool) :
      (-(n : ℝ) * historyMass P h / 2 ≤ -visits x h) ↔
        visits x h ≤ (n : ℝ) * historyMass P h / 2 := by
      constructor <;> intro hx <;> linarith
  simp_rw [hevent] at hc
  have hm : expect P (fun z ↦ Real.exp (1 * f z)) ≤
      Real.exp (-5 * historyMass P h / 8) := by
    simp only [one_mul, f]
    rw [masked_exponential P h (fun _ ↦ -1)]
    have hconst : expect (bernoulli (conditional P h) (conditional_parameters P h))
        (fun _ ↦ Real.exp (-1)) = Real.exp (-1) := by
      simp only [expect, ← sum_mul, Law.total, one_mul]
    rw [hconst]
    have h1 := mul_le_mul_of_nonneg_left exp_neg_one_bound (historyMass_pos P h).le
    have h2 := Real.add_one_le_exp (-5 * historyMass P h / 8)
    linarith
  have hn0 : 0 ≤ expect P (fun z ↦ Real.exp (1 * f z)) :=
    sum_nonneg (fun _ _ ↦ mul_nonneg (P.positive _).le (Real.exp_pos _).le)
  have hp := pow_le_pow_left₀ hn0 hm n
  calc _ ≤ Real.exp (-(1 : ℝ) * (-(n : ℝ) * historyMass P h / 2)) *
        (expect P (fun z ↦ Real.exp (1 * f z))) ^ n := hc
    _ ≤ Real.exp (-(1 : ℝ) * (-(n : ℝ) * historyMass P h / 2)) *
        (Real.exp (-5 * historyMass P h / 8)) ^ n :=
      mul_le_mul_of_nonneg_left hp (Real.exp_pos _).le
    _ = _ := by rw [← Real.exp_nat_mul, ← Real.exp_add]; congr 1; ring

theorem any_too_few_visits (P : Law (H × Bool)) (n : ℕ) {π : ℝ}
    (hπ : ∀ h, π ≤ historyMass P h) :
    probability (iid n P)
        (fun x ↦ ∃ h, visits x h ≤ n * historyMass P h / 2) ≤
      Fintype.card H * Real.exp (-(n : ℝ) * π / 8) := by
  have hb (h : H) : Real.exp (-(n : ℝ) * historyMass P h / 8) ≤
      Real.exp (-(n : ℝ) * π / 8) := by
    apply Real.exp_le_exp.mpr
    have hm := mul_le_mul_of_nonneg_left (hπ h) (Nat.cast_nonneg n)
    linarith
  have hs := sum_le_sum (fun h (_ : h ∈ univ) ↦ (too_few_visits P h n).trans (hb h))
  have hp := probability_exists (iid n P)
    (fun h x ↦ visits x h ≤ n * historyMass P h / 2)
  simp only [sum_const, card_univ, nsmul_eq_mul] at hs
  exact hp.trans hs

#print axioms too_few_visits
#print axioms any_too_few_visits

#print axioms ratio_tail
#print axioms transition_error
#print axioms simultaneous_transition_error
#print axioms transition_budget

end TimeFrequency.ContextTree.DepthEvidence
