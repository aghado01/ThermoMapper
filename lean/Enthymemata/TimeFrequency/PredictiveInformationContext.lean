import Enthymemata.TimeFrequency.PredictiveInformationPositive
import Enthymemata.TimeFrequency.ContextTreeFourier

/-! # The all-depth context-tree collision hides positive incremental information -/

namespace TimeFrequency.PredictiveActivity

open Finset hiding expect
open ContextTree ContextTree.Collision ContextTree.DepthEvidence

noncomputable def historyLaw (m : ℕ) : Law (Word (D m)) where
  mass _ := uniform (D m)
  positive _ := by unfold uniform; positivity
  total := uniform_total (D m)

def recentHistory (m : ℕ) (h : Word (D m)) : Word (m + 3) := Fin.tail h

theorem recent_surjective (m : ℕ) : Function.Surjective (recentHistory m) := by
  intro h
  exact ⟨Fin.cons false h, rfl⟩

theorem deep_transition_range (m : ℕ) {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (h : Word (D m)) : transition m a h true ∈ Set.Ioo 0 1 := by
  have hs := transition_stochastic m ha hu h
  have ht := hs.2
  simp only [Fintype.sum_bool] at ht
  exact ⟨hs.1 true, by linarith [hs.1 false]⟩

theorem deep_incremental_information (m : ℕ) {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2) :
    0 < binaryInformation (historyLaw m) (recentHistory m) (recent_surjective m)
      (fun h ↦ transition m a h true) (deep_transition_range m ha hu) := by
  apply information_positive
  exact exact_memory m ha

theorem fair_incremental_information (m : ℕ) :
    binaryInformation (historyLaw m) (recentHistory m) (recent_surjective m)
      (fun _ ↦ 1 / 2) (fun _ ↦ by constructor <;> norm_num) = 0 :=
  information_constant _ _ _ (by constructor <;> norm_num)

theorem rectangular_information_obstruction (m : ℕ) {a : ℝ}
    (ha : 0 < a) (hu : a < 1 / 2) :
    (∀ w o (ho : o + w ≤ L m) (v : Fin w → ℝ),
      bin (weight m a) (fun x ↦ powerVector (slice x o w ho)) v =
        bin (fun _ ↦ uniform w) powerVector v) ∧
    binaryInformation (historyLaw m) (recentHistory m) (recent_surjective m)
      (fun _ ↦ 1 / 2) (fun _ ↦ by constructor <;> norm_num) = 0 ∧
    0 < binaryInformation (historyLaw m) (recentHistory m) (recent_surjective m)
      (fun h ↦ transition m a h true) (deep_transition_range m ha hu) := by
  exact ⟨fun _ _ h _ ↦ all_width_power_laws h a _, fair_incremental_information m,
    deep_incremental_information m ha hu⟩

#print axioms rectangular_information_obstruction
#print axioms recent_surjective
#print axioms deep_transition_range
#print axioms deep_incremental_information
#print axioms fair_incremental_information

end TimeFrequency.PredictiveActivity
