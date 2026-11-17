import Enthymemata.TimeFrequency.MemoryDepthEmpirical
import Enthymemata.TimeFrequency.MemoryDepthReconstruction

/-! # The six-view depth-four sampling experiment -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset DepthFour

noncomputable def sourceBlockLaw (P : Source) (hp : ∀ x, 0 < P.law 5 x) : Law (Word 5) where
  mass := P.law 5
  positive := hp
  total := P.total 5

def powerStatistic (j : Fin 8) (x : Word 5) : ℚ :=
  ![(dc (slice x 0 2 (by decide)) : ℚ) / 4,
    (dc (slice x 0 3 (by decide)) : ℚ) / 9,
    (dc (slice x 0 4 (by decide)) : ℚ) / 16,
    (dc x : ℚ) / 25,
    if dc (slice x 0 4 (by decide)) = 16 then 1 else 0,
    if dc x = 25 then 1 else 0,
    if taper 3 x = 16 then 1 else 0,
    if taper 2 x = 16 then 1 else 0] j

def realPowerStatistic (j : Fin 8) (x : Word 5) : ℝ := (powerStatistic j x : ℝ)

private theorem integer_power_ranges (x : Word 5) :
    dc (slice x 0 2 (by decide)) ∈ Set.Icc (0 : ℤ) 4 ∧
    dc (slice x 0 3 (by decide)) ∈ Set.Icc (0 : ℤ) 9 ∧
    dc (slice x 0 4 (by decide)) ∈ Set.Icc (0 : ℤ) 16 ∧
    dc x ∈ Set.Icc (0 : ℤ) 25 := by decide +revert

private theorem normalized_power_range {p c : ℤ} (hc : 0 < c)
    (hp : p ∈ Set.Icc 0 c) : (p : ℚ) / c ∈ Set.Icc (0 : ℚ) 1 := by
  have hcQ : (0 : ℚ) < c := by exact_mod_cast hc
  constructor
  · exact div_nonneg (by exact_mod_cast hp.1) hcQ.le
  · exact (div_le_one hcQ).mpr (by exact_mod_cast hp.2)

theorem power_statistic_range (j : Fin 8) (x : Word 5) :
    powerStatistic j x ∈ Set.Icc (0 : ℚ) 1 := by
  have h := integer_power_ranges x
  fin_cases j
  · exact normalized_power_range (by norm_num : (0 : ℤ) < 4) h.1
  · exact normalized_power_range (by norm_num : (0 : ℤ) < 9) h.2.1
  · exact normalized_power_range (by norm_num : (0 : ℤ) < 16) h.2.2.1
  · exact normalized_power_range (by norm_num : (0 : ℤ) < 25) h.2.2.2
  · change (if dc (slice x 0 4 (by decide)) = 16 then (1 : ℚ) else 0) ∈ Set.Icc 0 1
    split_ifs <;> norm_num
  · change (if dc x = 25 then (1 : ℚ) else 0) ∈ Set.Icc 0 1
    split_ifs <;> norm_num
  · change (if taper 3 x = 16 then (1 : ℚ) else 0) ∈ Set.Icc 0 1
    split_ifs <;> norm_num
  · change (if taper 2 x = 16 then (1 : ℚ) else 0) ∈ Set.Icc 0 1
    split_ifs <;> norm_num

theorem real_power_statistic_range (j : Fin 8) (x : Word 5) :
    realPowerStatistic j x ∈ Set.Icc (0 : ℝ) 1 := by
  have h := power_statistic_range j x
  change 0 ≤ powerStatistic j x ∧ powerStatistic j x ≤ 1 at h
  change 0 ≤ (powerStatistic j x : ℝ) ∧ (powerStatistic j x : ℝ) ≤ 1
  exact_mod_cast h

private theorem prefix_power_expect (P : Source) (hp : ∀ x, 0 < P.law 5 x)
    (w : ℕ) (hw : 0 + w ≤ 5) (c : ℝ) :
    expect (sourceBlockLaw P hp) (fun x ↦ (dc (slice x 0 w hw) : ℝ) / c) =
      powerMean P w / c := by
  simp only [expect, sourceBlockLaw, ← mul_div_assoc, ← sum_div]
  change mean (P.law 5) (fun x ↦ dc (slice x 0 w hw)) / c = _
  rw [P.translation]
  rfl

private theorem prefix_event_expect (P : Source) (hp : ∀ x, 0 < P.law 5 x)
    (w : ℕ) (hw : 0 + w ≤ 5) (v : ℤ) :
    expect (sourceBlockLaw P hp)
      (fun x ↦ if dc (slice x 0 w hw) = v then 1 else 0) = bin (P.law w) dc v := by
  have h := P.translation 5 0 w hw (fun x ↦ hit (dc x) v)
  simpa [expect, sourceBlockLaw, mean, hit, bin, mul_ite] using h

private theorem event_expect (P : Source) (hp : ∀ x, 0 < P.law 5 x)
    (f : Word 5 → ℤ) (v : ℤ) :
    expect (sourceBlockLaw P hp) (fun x ↦ if f x = v then 1 else 0) =
      bin (P.law 5) f v := by
  simp [expect, sourceBlockLaw, bin, mul_ite]

private theorem cast_indicator (A : Prop) [Decidable A] :
    (((if A then (1 : ℚ) else 0) : ℚ) : ℝ) = if A then 1 else 0 := by
  split_ifs <;> norm_num

theorem power_statistic_expect (P : Source) (hp : ∀ x, 0 < P.law 5 x) (j : Fin 8) :
    expect (sourceBlockLaw P hp) (realPowerStatistic j) = powerSummaries P j := by
  fin_cases j
  · change expect (sourceBlockLaw P hp)
      (fun x ↦ (((dc (slice x 0 2 (by decide)) : ℚ) / 4 : ℚ) : ℝ)) = powerMean P 2 / 4
    push_cast
    exact prefix_power_expect P hp 2 (by decide) 4
  · change expect (sourceBlockLaw P hp)
      (fun x ↦ (((dc (slice x 0 3 (by decide)) : ℚ) / 9 : ℚ) : ℝ)) = powerMean P 3 / 9
    push_cast
    exact prefix_power_expect P hp 3 (by decide) 9
  · change expect (sourceBlockLaw P hp)
      (fun x ↦ (((dc (slice x 0 4 (by decide)) : ℚ) / 16 : ℚ) : ℝ)) = powerMean P 4 / 16
    push_cast
    exact prefix_power_expect P hp 4 (by decide) 16
  · change expect (sourceBlockLaw P hp)
      (fun x ↦ (((dc x : ℚ) / 25 : ℚ) : ℝ)) = powerMean P 5 / 25
    push_cast
    simp only [expect, sourceBlockLaw, ← mul_div_assoc, ← sum_div, powerMean, mean]
  · change expect (sourceBlockLaw P hp)
      (fun x ↦ (((if dc (slice x 0 4 (by decide)) = 16 then (1 : ℚ) else 0) : ℚ) : ℝ)) =
        bin (P.law 4) dc 16
    simp_rw [cast_indicator]
    exact prefix_event_expect P hp 4 (by decide) 16
  · change expect (sourceBlockLaw P hp)
      (fun x ↦ (((if dc x = 25 then (1 : ℚ) else 0) : ℚ) : ℝ)) = bin (P.law 5) dc 25
    simp_rw [cast_indicator]
    exact event_expect P hp dc 25
  · change expect (sourceBlockLaw P hp)
      (fun x ↦ (((if taper 3 x = 16 then (1 : ℚ) else 0) : ℚ) : ℝ)) = bin (P.law 5) (taper 3) 16
    simp_rw [cast_indicator]
    exact event_expect P hp (taper 3) 16
  · change expect (sourceBlockLaw P hp)
      (fun x ↦ (((if taper 2 x = 16 then (1 : ℚ) else 0) : ℚ) : ℝ)) = bin (P.law 5) (taper 2) 16
    simp_rw [cast_indicator]
    exact event_expect P hp (taper 2) 16

/-- Means and top-bin frequencies share the width-four and width-five views. -/
def summaryView : Fin 8 → Fin 6 := ![0, 1, 2, 3, 2, 3, 4, 5]

noncomputable def recordedSummaries {n : ℕ} (x : Fin 6 → Fin n → Word 5) (j : Fin 8) :=
  sampleMean (realPowerStatistic j) (x (summaryView j))

theorem unpaired_summary_error (P : Source) (hp : ∀ x, 0 < P.law 5 x)
    (n : ℕ) (hn : 0 < n) {ε : ℝ} (he : 0 < ε) :
    probability (iid 6 (iid n (sourceBlockLaw P hp))) (fun x ↦ ∃ j,
      ε ≤ |recordedSummaries x j - powerSummaries P j|) ≤
      16 * Real.exp (-2 * n * ε ^ 2) := by
  simpa only [recordedSummaries, power_statistic_expect, Fintype.card_fin, Nat.cast_ofNat,
    show (2 : ℝ) * 8 = 16 from by norm_num] using
    separate_summary_error (sourceBlockLaw P hp) realPowerStatistic
      real_power_statistic_range 6 n summaryView hn he

theorem unpaired_summary_budget (P : Source) (hp : ∀ x, 0 < P.law 5 x)
    (n : ℕ) {ε δ : ℝ} (he : 0 < ε) (hδ : 0 < δ) (hδu : δ < 1)
    (hn : ⌈Real.log (16 / δ) / (2 * ε ^ 2)⌉₊ ≤ n) :
    probability (iid 6 (iid n (sourceBlockLaw P hp))) (fun x ↦ ∃ j,
      ε ≤ |recordedSummaries x j - powerSummaries P j|) ≤ δ := by
  have hcard : (2 * Fintype.card (Fin 8) : ℝ) = 16 := by norm_num
  have h := separate_summary_budget (sourceBlockLaw P hp) realPowerStatistic
    real_power_statistic_range 6 n summaryView he hδ hδu (by decide) (by rwa [hcard])
  simpa only [recordedSummaries, power_statistic_expect] using h

#print axioms power_statistic_range
#print axioms real_power_statistic_range
#print axioms power_statistic_expect
#print axioms unpaired_summary_error
#print axioms unpaired_summary_budget

end TimeFrequency.ContextTree.DepthEvidence
