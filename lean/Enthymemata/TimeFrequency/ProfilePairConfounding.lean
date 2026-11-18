import Enthymemata.TimeFrequency.ProfileAcousticPairs

/-! # Holding the speaker fixed while crossing a binary recording condition -/

namespace TimeFrequency.LinkedProfiles
open MeasureTheory Set
open scoped ENNReal
noncomputable section

/-- A binary factor selects the high-fraction emission probability. -/
def emissionChance (label : Bool) : ℝ := if label then 3/4 else 1/4

private theorem chance_mem (b : Bool) : emissionChance b ∈ Icc 0 1 := by
  cases b <;> norm_num [Set.mem_Icc, emissionChance]

def conditionalPair (b c : Bool) :=
  (bitLaw (emissionChance b)).prod (bitLaw (emissionChance c))

/-- Speaker and first recording condition are independent fair bits. The speaker returns.
The protocol chooses the second condition, and the mechanism chooses which factor drives
each emission. -/
def recordingExperiment (driver : Bool → Bool → Bool) (secondCondition : Bool → Bool) :=
  let record := fun s c ↦ conditionalPair (driver s c) (driver s (secondCondition c))
  balanced (balanced (record true true) (record true false))
    (balanced (record false true) (record false false))

private theorem balanced_self {X : Type*} [MeasurableSpace X] (p : Measure X) :
    balanced p p = p := by
  ext A
  simp only [balanced, Measure.smul_apply, Measure.add_apply, smul_eq_mul, ← two_mul, ← mul_assoc]
  rw [ENNReal.inv_mul_cancel (by norm_num) (by norm_num), one_mul]

private theorem balanced_real {X : Type*} [MeasurableSpace X] (p q : Measure X)
    [IsFiniteMeasure p] [IsFiniteMeasure q] (A : Set X) :
    (balanced p q).real A = (p.real A + q.real A) / 2 := by
  rw [balanced, measureReal_ennreal_smul_apply, measureReal_add_apply]
  norm_num
  ring

theorem shared_recording_law :
    recordingExperiment (fun s _ ↦ s) id = recordingExperiment (fun _ c ↦ c) id := by
  simp [recordingExperiment, balanced_self]

private theorem conditional_high (b c : Bool) :
    (conditionalPair b c).real {(true,true)} = emissionChance b * emissionChance c := by
  rw [conditionalPair, ← singleton_prod_singleton, measureReal_prod_prod,
    bitLaw_mass (chance_mem b), bitLaw_mass (chance_mem c)]
  rfl

private instance (b c : Bool) : IsFiniteMeasure (conditionalPair b c) := by
  unfold conditionalPair
  infer_instance

theorem crossed_recording_prediction :
    (readPairs waveletReadout (recordingExperiment (fun s _ ↦ s) Bool.not)).real
      {(waveletReadout true,waveletReadout true)} = 5/16 ∧
    (readPairs waveletReadout (recordingExperiment (fun _ c ↦ c) Bool.not)).real
      {(waveletReadout true,waveletReadout true)} = 3/16 := by
  simp only [readPairs_cell _ waveletReadout_injective, recordingExperiment, Bool.not_true,
    Bool.not_false, balanced_self]
  simp only [balanced_real, conditional_high]
  norm_num [emissionChance]

/-- Independently redrawing the recording condition gives independent fair marginal draws
under the condition-only mechanism. -/
theorem independent_recording_prediction :
    (balanced (recordingExperiment (fun _ c ↦ c) id)
      (recordingExperiment (fun _ c ↦ c) Bool.not)).real {(true,true)} = 1/4 := by
  simp only [recordingExperiment, id_eq, Bool.not_true, Bool.not_false, balanced_self]
  simp only [balanced_real, conditional_high]
  norm_num [emissionChance]

#print axioms shared_recording_law
#print axioms crossed_recording_prediction
#print axioms independent_recording_prediction

end
end TimeFrequency.LinkedProfiles
