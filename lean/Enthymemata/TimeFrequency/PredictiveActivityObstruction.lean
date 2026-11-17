import Enthymemata.TimeFrequency.PredictiveActivityEntropy

/-! # A one-bit obstruction to marginal future-activity objectives -/

namespace TimeFrequency.PredictiveActivity

open Finset hiding expect
open ContextTree.DepthEvidence

noncomputable def fairHistory : Law Bool := bernoulli (1 / 2) (by constructor <;> norm_num)

noncomputable def parityLaw (older : Bool) : ActivityLaw where
  prob := if older then ![0, 1 / 2, 1 / 2, 0] else ![1 / 2, 0, 0, 1 / 2]
  nonneg y := by cases older <;> fin_cases y <;> norm_num
  total := by cases older <;> norm_num [Fin.sum_univ_succ]

noncomputable def uniformActivity : ActivityLaw where
  prob _ := 1 / 4
  nonneg _ := by norm_num
  total := by norm_num

def constantCode (_ : Bool) : Unit := ()

theorem constant_code_surjective : Function.Surjective constantCode := by
  intro z
  cases z
  exact ⟨false, rfl⟩

theorem constant_encoded_prob (z : Unit) (y : Fin 4) :
    (encodedLaw fairHistory constantCode constant_code_surjective parityLaw z).prob y = 1 / 4 := by
  cases z
  fin_cases y <;> norm_num [encodedLaw, conditionalMean, fiberSum, codeMass, fairHistory,
    bernoulli, bitMass, constantCode, parityLaw, Fintype.sum_bool]

theorem marginal_means (older : Bool) :
    targetMeans (parityLaw older) 0 = 1 / 2 ∧
    targetMeans (parityLaw older) 1 = 1 / 2 := by
  cases older <;> norm_num [targetMeans, parityLaw, Matrix.cons_val_two, Matrix.cons_val_three]

theorem constant_marginal_optimum (older : Bool) (j : Fin 2) :
    binaryRisk (targetMeans (parityLaw older) (Fin.castLE (by decide) j)) (1 / 2) -
      binaryRisk (targetMeans (parityLaw older) (Fin.castLE (by decide) j))
        (targetMeans (parityLaw older) (Fin.castLE (by decide) j)) = 0 := by
  have hm := marginal_means older
  fin_cases j
  · change binaryRisk (targetMeans (parityLaw older) 0) (1 / 2) -
      binaryRisk (targetMeans (parityLaw older) 0) (targetMeans (parityLaw older) 0) = 0
    rw [hm.1, sub_self]
  · change binaryRisk (targetMeans (parityLaw older) 1) (1 / 2) -
      binaryRisk (targetMeans (parityLaw older) 1) (targetMeans (parityLaw older) 1) = 0
    rw [hm.2, sub_self]

theorem parity_discards_one_bit (older : Bool) :
    klNats (parityLaw older) uniformActivity / Real.log 2 = 1 := by
  have hl : Real.log 2 ≠ 0 := ne_of_gt (Real.log_pos (by norm_num))
  cases older <;> norm_num [klNats, parityLaw, uniformActivity, Fin.sum_univ_succ]
  all_goals field_simp; ring

theorem constant_encoder_discards_one_bit :
    expect fairHistory (fun older ↦ klNats (parityLaw older)
      (encodedLaw fairHistory constantCode constant_code_surjective parityLaw
        (constantCode older))) / Real.log 2 = 1 := by
  have he (older : Bool) : klNats (parityLaw older)
      (encodedLaw fairHistory constantCode constant_code_surjective parityLaw
        (constantCode older)) = klNats (parityLaw older) uniformActivity := by
    simp only [klNats, constant_encoded_prob, uniformActivity]
  simp only [expect, he, sum_div, mul_div_assoc, parity_discards_one_bit, mul_one,
    fairHistory.total]

theorem one_bit_identity_retains (older : Bool) (y : Fin 4) :
    (encodedLaw fairHistory id Function.surjective_id parityLaw older).prob y =
      (parityLaw older).prob y := by
  simp [encodedLaw, conditionalMean, codeMass, fiberSum, fairHistory.positive older |>.ne']

theorem eligible_code_sizes :
    Fintype.card Unit ≤ 2 ∧ Fintype.card Bool ≤ 2 := by decide

def pulseWave {n : ℕ} (pulse : Fin n → ℝ) (active : Bool) : Fin n → ℝ :=
  fun i ↦ if active then pulse i else 0

noncomputable def energy {n : ℕ} (wave : Fin n → ℝ) : ℝ := ∑ i, wave i ^ 2

noncomputable def activity {n : ℕ} (threshold : ℝ) (wave : Fin n → ℝ) : Bool :=
  decide (threshold < energy wave)

theorem pulse_activity {n : ℕ} (pulse : Fin n → ℝ) {threshold : ℝ}
    (ht : 0 < threshold) (hp : threshold < energy pulse) (b : Bool) :
    activity threshold (pulseWave pulse b) = b := by
  cases b
  · simp [activity, pulseWave, energy, not_lt_of_ge ht.le]
  · change decide (threshold < energy pulse) = true
    exact decide_eq_true hp

/-- The older segment and the two future channels use separate, ordered slots.
The recent slot is common silence. Future labels share the same future slot. -/
def waveformRecord {n : ℕ} (pulse : Fin n → ℝ) (older a b : Bool) :
    (Fin n → ℝ) × (Fin n → ℝ) × ((Fin n → ℝ) × (Fin n → ℝ)) :=
  (pulseWave pulse older, fun _ ↦ 0, pulseWave pulse a, pulseWave pulse b)

theorem waveform_labels {n : ℕ} (pulse : Fin n → ℝ) {threshold : ℝ}
    (ht : 0 < threshold) (hp : threshold < energy pulse) (older a b : Bool) :
    activity threshold (waveformRecord pulse older a b).1 = older ∧
    activity threshold (waveformRecord pulse older a b).2.2.1 = a ∧
    activity threshold (waveformRecord pulse older a b).2.2.2 = b := by
  exact ⟨pulse_activity pulse ht hp older, pulse_activity pulse ht hp a,
    pulse_activity pulse ht hp b⟩

theorem constant_encoder_entropy_gap :
    conditionalEntropyBits fairHistory (fun older ↦
      encodedLaw fairHistory constantCode constant_code_surjective parityLaw
        (constantCode older)) - conditionalEntropyBits fairHistory parityLaw = 1 := by
  rw [← information_discarded_identity fairHistory constantCode constant_code_surjective
    parityLaw (fun z y ↦ by rw [constant_encoded_prob]; norm_num)]
  exact constant_encoder_discards_one_bit

def slotTime : Fin 3 → ℤ := ![-2, -1, 1]

theorem causal_slots : slotTime 0 < 0 ∧ slotTime 1 < 0 ∧ 0 < slotTime 2 := by
  change (-2 : ℤ) < 0 ∧ (-1 : ℤ) < 0 ∧ 0 < (1 : ℤ)
  norm_num

/-- A concrete nonzero pulse and threshold realize every activity label. -/
theorem concrete_waveform_realization (older a b : Bool) :
    activity (1 / 2) (waveformRecord (fun _ : Fin 1 ↦ (1 : ℝ)) older a b).1 = older ∧
    activity (1 / 2) (waveformRecord (fun _ : Fin 1 ↦ (1 : ℝ)) older a b).2.2.1 = a ∧
    activity (1 / 2) (waveformRecord (fun _ : Fin 1 ↦ (1 : ℝ)) older a b).2.2.2 = b := by
  apply waveform_labels
  · norm_num
  · norm_num [energy]

#print axioms causal_slots
#print axioms concrete_waveform_realization
#print axioms constant_encoder_entropy_gap
#print axioms constant_code_surjective
#print axioms constant_encoded_prob
#print axioms marginal_means
#print axioms constant_marginal_optimum
#print axioms parity_discards_one_bit
#print axioms constant_encoder_discards_one_bit
#print axioms one_bit_identity_retains
#print axioms eligible_code_sizes
#print axioms pulse_activity
#print axioms waveform_labels

end TimeFrequency.PredictiveActivity
