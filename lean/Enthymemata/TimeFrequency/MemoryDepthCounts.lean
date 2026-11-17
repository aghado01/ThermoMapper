import Enthymemata.TimeFrequency.MemoryDepthConcentration

/-! # Counts, visits, and transition errors in independent history-response records -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset

variable {H : Type*} [Fintype H] [DecidableEq H]

noncomputable def historyMass (P : Law (H × Bool)) (h : H) :=
  P.mass (h, false) + P.mass (h, true)

noncomputable def conditional (P : Law (H × Bool)) (h : H) :=
  P.mass (h, true) / historyMass P h

omit [DecidableEq H] in
theorem historyMass_pos (P : Law (H × Bool)) (h : H) : 0 < historyMass P h :=
  add_pos (P.positive _) (P.positive _)

omit [DecidableEq H] in
theorem conditional_parameters (P : Law (H × Bool)) (h : H) :
    conditional P h ∈ Set.Ioo 0 1 := by
  constructor
  · exact div_pos (P.positive _) (historyMass_pos P h)
  · apply (div_lt_one (historyMass_pos P h)).mpr
    dsimp [historyMass]
    linarith [P.positive (h, false)]

omit [DecidableEq H] in
theorem conditional_masses (P : Law (H × Bool)) (h : H) :
    historyMass P h * conditional P h = P.mass (h, true) ∧
    historyMass P h * (1 - conditional P h) = P.mass (h, false) := by
  have he : historyMass P h * conditional P h = P.mass (h, true) := by
    unfold conditional
    field_simp [(historyMass_pos P h).ne']
  exact ⟨he, by dsimp [historyMass] at *; nlinarith only [he]⟩

theorem history_expect (P : Law (H × Bool)) (h : H) (f : Bool → ℝ) (c : ℝ) :
    expect P (fun z ↦ if z.1 = h then f z.2 else c) =
      c + ∑ b, P.mass (h, b) * (f b - c) := by
  have he (z : H × Bool) : (if z.1 = h then f z.2 else c) =
      c + if z.1 = h then f z.2 - c else 0 := by
    split_ifs <;> ring
  simp only [expect, he, mul_add, sum_add_distrib]
  rw [← sum_mul, P.total, one_mul]
  congr 1
  simp [Fintype.sum_prod_type, mul_ite]

theorem masked_exponential (P : Law (H × Bool)) (h : H) (f : Bool → ℝ) :
    expect P (fun z ↦ Real.exp (if z.1 = h then f z.2 else 0)) =
      1 - historyMass P h + historyMass P h *
        expect (bernoulli (conditional P h) (conditional_parameters P h))
          (fun b ↦ Real.exp (f b)) := by
  have he (z : H × Bool) : Real.exp (if z.1 = h then f z.2 else 0) =
      if z.1 = h then Real.exp (f z.2) else 1 := by split_ifs <;> simp
  simp_rw [he]
  rw [history_expect P h (fun b ↦ Real.exp (f b)) 1]
  simp only [Fintype.sum_bool, expect, bernoulli, bitMass, Bool.false_eq_true,
    if_false, if_true]
  have hm := conditional_masses P h
  have h1 := congrArg (fun z : ℝ ↦ z * Real.exp (f true)) hm.1
  have h0 := congrArg (fun z : ℝ ↦ z * Real.exp (f false)) hm.2
  dsimp [historyMass] at *
  nlinarith only [h1, h0]

noncomputable def visits {n : ℕ} (x : Fin n → H × Bool) (h : H) : ℝ :=
  ∑ i, if (x i).1 = h then 1 else 0

noncomputable def positives {n : ℕ} (x : Fin n → H × Bool) (h : H) : ℝ :=
  ∑ i, if (x i).1 = h then bitValue (x i).2 else 0

noncomputable def estimate {n : ℕ} (x : Fin n → H × Bool) (h : H) :=
  if visits x h = 0 then 1 / 2 else positives x h / visits x h

omit [Fintype H] in
theorem visits_nonneg {n : ℕ} (x : Fin n → H × Bool) (h : H) : 0 ≤ visits x h :=
  sum_nonneg (fun _ _ ↦ by split_ifs <;> norm_num)

omit [Fintype H] in
theorem positives_nonneg {n : ℕ} (x : Fin n → H × Bool) (h : H) : 0 ≤ positives x h :=
  sum_nonneg (fun i _ ↦ by
    by_cases hi : (x i).1 = h
    · simp only [if_pos hi, bitValue]; split_ifs <;> norm_num
    · simp [hi])

omit [Fintype H] in
theorem masked_sum {n : ℕ} (x : Fin n → H × Bool) (h : H) (q t s : ℝ) :
    (∑ i, if (x i).1 = h then s * (bitValue (x i).2 - q) - t else 0) =
      s * (positives x h - q * visits x h) - t * visits x h := by
  simp only [positives, visits, mul_sum, ← sum_sub_distrib]
  apply sum_congr rfl
  intro i _
  split_ifs <;> ring

theorem ratio_mgf (P : Law (H × Bool)) (h : H) {t s : ℝ}
    (ht : 0 ≤ t) (hu : t ≤ 1 / 2) (hs : s ^ 2 = 1) :
    expect P (fun z ↦ Real.exp (4 * t *
        (if z.1 = h then s * (bitValue z.2 - conditional P h) - t else 0))) ≤
      Real.exp (-historyMass P h * t ^ 2) := by
  let q := conditional P h
  let B := bernoulli q (conditional_parameters P h)
  have hb : expect B bitValue = q := by simp [B, expect, bernoulli, bitMass, bitValue]
  have hmgf := bounded_mgf B bitValue (by intro b; cases b <;> norm_num [bitValue]) (4 * t * s)
  rw [hb] at hmgf
  have hexp : (4 * t * s) ^ 2 / 8 = 2 * t ^ 2 := by nlinarith [hs]
  rw [hexp] at hmgf
  have he (b : Bool) : 4 * t * (s * (bitValue b - q) - t) =
      -4 * t ^ 2 + (4 * t * s) * (bitValue b - q) := by ring
  have hm : expect B (fun b ↦ Real.exp (4 * t * (s * (bitValue b - q) - t))) ≤
      Real.exp (-2 * t ^ 2) := by
    simp_rw [he, Real.exp_add]
    have hf : expect B (fun b ↦ Real.exp (-4 * t ^ 2) *
        Real.exp (4 * t * s * (bitValue b - q))) =
        Real.exp (-4 * t ^ 2) * expect B
          (fun b ↦ Real.exp (4 * t * s * (bitValue b - q))) := by
      simp only [expect, mul_sum]; apply sum_congr rfl; intros; ring
    rw [hf]
    calc _ ≤ Real.exp (-4 * t ^ 2) * Real.exp (2 * t ^ 2) :=
          mul_le_mul_of_nonneg_left hmgf (Real.exp_pos _).le
      _ = _ := by rw [← Real.exp_add]; congr 1; ring
  have hez (z : H × Bool) : 4 * t *
      (if z.1 = h then s * (bitValue z.2 - q) - t else 0) =
      if z.1 = h then 4 * t * (s * (bitValue z.2 - q) - t) else 0 := by
    split_ifs <;> ring
  change expect P (fun z ↦ Real.exp (4 * t *
    (if z.1 = h then s * (bitValue z.2 - q) - t else 0))) ≤ _
  simp_rw [hez]
  rw [masked_exponential P h (fun b ↦ 4 * t * (s * (bitValue b - q) - t))]
  have hπ := (historyMass_pos P h).le
  have h1 := mul_le_mul_of_nonneg_left hm hπ
  have h2 := mul_le_mul_of_nonneg_left (exp_decrement ht hu) hπ
  have h3 := Real.add_one_le_exp (-historyMass P h * t ^ 2)
  dsimp [B] at h1
  nlinarith only [h1, h2, h3]

#print axioms historyMass_pos
#print axioms conditional_parameters
#print axioms conditional_masses
#print axioms history_expect
#print axioms masked_exponential
#print axioms visits_nonneg
#print axioms positives_nonneg
#print axioms masked_sum
#print axioms ratio_mgf

end TimeFrequency.ContextTree.DepthEvidence
