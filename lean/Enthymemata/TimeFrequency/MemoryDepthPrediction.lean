import Enthymemata.TimeFrequency.MemoryDepthConfidence
import Enthymemata.TimeFrequency.DepthFourPrediction

/-! # Prediction tolerance despite uncertainty about deeper memory -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset

noncomputable def sourceExcess (P : Source) (qhat : Word 4 → ℝ) : ℝ :=
  ∑ h, P.law 4 h * (qhat h - nextKernel 4 (P.law 5) h true) ^ 2

theorem uniform_prediction_tolerance (P : Source) (qhat : Word 4 → ℝ) {t : ℝ}
    (herr : ∀ h, |qhat h - nextKernel 4 (P.law 5) h true| ≤ t) :
    sourceExcess P qhat ≤ t ^ 2 := by
  have hb (h : Word 4) :
      (qhat h - nextKernel 4 (P.law 5) h true) ^ 2 ≤ t ^ 2 := by
    have hp := pow_le_pow_left₀ (abs_nonneg _) (herr h) 2
    rwa [sq_abs] at hp
  have hs := sum_le_sum (fun h (_ : h ∈ univ) ↦
    mul_le_mul_of_nonneg_left (hb h) (P.nonneg 4 h))
  simpa only [sourceExcess, ← sum_mul, P.total, one_mul] using hs

theorem fair_prediction_loss (a : ℝ) :
    DepthFour.excess (fun _ ↦ 1 / 2) (DepthFour.predict a 0) = a ^ 2 / 2 := by
  have hz : DepthFour.predict 0 0 = fun _ ↦ (1 / 2 : ℝ) := by
    funext h
    simp [DepthFour.predict, DepthFour.prediction_formula]
  have he := DepthFour.prediction_distance 0 0 a 0
  rw [hz] at he
  convert he using 1
  ring

theorem arbitrarily_weak_deep_prediction {τ : ℝ} (hτ : 0 < τ) :
    ∃ a : ℝ, 0 < a ∧ a < 1 / 2 ∧ DepthFour.Admissible a 0 ∧
      (∃ h h' : Word 4, Fin.tail h = Fin.tail h' ∧
        DepthFour.kernel a 0 h true ≠ DepthFour.kernel a 0 h' true) ∧
      (∀ h : Word 3, (DepthFour.predict a 0 (Fin.cons false h) +
        DepthFour.predict a 0 (Fin.cons true h)) / 2 = 1 / 2) ∧
      (∀ w o (ho : o + w ≤ 5) (v : Fin w → ℝ),
        bin (DepthFour.weight a 0) (fun x ↦ powerVector (slice x o w ho)) v =
        bin (DepthFour.weight 0 0) (fun x ↦ powerVector (slice x o w ho)) v) ∧
      DepthFour.excess (fun _ ↦ 1 / 2) (DepthFour.predict a 0) ≤ τ := by
  let a := min (1 / 4 : ℝ) τ
  have ha : 0 < a := lt_min (by norm_num) hτ
  have hau : a ≤ 1 / 4 := min_le_left _ _
  have hat : a ≤ τ := min_le_right _ _
  refine ⟨a, ha, by linarith, ?_, DepthFour.exact_four (Or.inl ha.ne'),
    DepthFour.suffix_prediction a 0, ?_, ?_⟩
  · simp only [DepthFour.Admissible, abs_of_pos ha, abs_zero, add_zero]
    linarith
  · intro w o ho v
    rw [DepthFour.all_fourier_laws _ _ ho, DepthFour.all_fourier_laws _ _ ho]
  · rw [fair_prediction_loss]
    have hs := mul_le_mul_of_nonneg_left hau ha.le
    nlinarith

#print axioms uniform_prediction_tolerance
#print axioms fair_prediction_loss
#print axioms arbitrarily_weak_deep_prediction

end TimeFrequency.ContextTree.DepthEvidence
