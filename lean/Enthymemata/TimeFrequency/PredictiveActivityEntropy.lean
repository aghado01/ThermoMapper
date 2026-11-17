import Enthymemata.TimeFrequency.PredictiveActivityBits

/-! # Conditional entropy bookkeeping for a frozen representation -/

namespace TimeFrequency.PredictiveActivity

open Finset hiding expect
open ContextTree.DepthEvidence

variable {X Z : Type*} [Fintype X] [Finite Z] [DecidableEq Z]

noncomputable def conditionalEntropyBits (P : Law X) (p : X → ActivityLaw) : ℝ :=
  -(expect P (fun x ↦ ∑ y, (p x).prob y * Real.log ((p x).prob y))) / Real.log 2

theorem weighted_log_div {p q : ℝ} (hq : 0 < q) :
    p * Real.log (p / q) = p * Real.log p - p * Real.log q := by
  by_cases hz : p = 0
  · simp [hz]
  · rw [Real.log_div hz hq.ne']
    ring

theorem encoded_log_expect (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (p : X → ActivityLaw) :
    expect P (fun x ↦ ∑ y, (p x).prob y *
      Real.log ((encodedLaw P f hf p (f x)).prob y)) =
    expect P (fun x ↦ ∑ y, (encodedLaw P f hf p (f x)).prob y *
      Real.log ((encodedLaw P f hf p (f x)).prob y)) := by
  simp only [expect, mul_sum]
  conv_lhs => rw [sum_comm]
  conv_rhs => rw [sum_comm]
  apply sum_congr rfl
  intro y _
  exact encoded_expect P f hf (fun x ↦ (p x).prob y)
    (fun z ↦ Real.log ((encodedLaw P f hf p z).prob y))

theorem information_discarded_identity (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (p : X → ActivityLaw)
    (hq : ∀ z y, 0 < (encodedLaw P f hf p z).prob y) :
    expect P (fun x ↦ klNats (p x) (encodedLaw P f hf p (f x))) / Real.log 2 =
      conditionalEntropyBits P (fun x ↦ encodedLaw P f hf p (f x)) -
      conditionalEntropyBits P p := by
  have he : expect P (fun x ↦ klNats (p x) (encodedLaw P f hf p (f x))) =
      expect P (fun x ↦ ∑ y, (p x).prob y * Real.log ((p x).prob y)) -
      expect P (fun x ↦ ∑ y, (encodedLaw P f hf p (f x)).prob y *
        Real.log ((encodedLaw P f hf p (f x)).prob y)) := by
    rw [← encoded_log_expect P f hf p]
    simp only [klNats, weighted_log_div (hq _ _),
      sum_sub_distrib, expect, mul_sub]
  rw [he]
  unfold conditionalEntropyBits
  ring

theorem information_loss_from_three_heads (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (p : X → ActivityLaw) (head : Z → Fin 3 → ℝ)
    {β : ℝ} (hb : 0 < β) (hq : ∀ z y, β ≤ (encodedLaw P f hf p z).prob y) :
    conditionalEntropyBits P (fun x ↦ encodedLaw P f hf p (f x)) -
      conditionalEntropyBits P p ≤ 8 * headExcess P f p head / (β * Real.log 2) := by
  rw [← information_discarded_identity P f hf p (fun z y ↦ hb.trans_le (hq z y))]
  exact discarded_bits_bound P f hf p head hb hq

#print axioms weighted_log_div
#print axioms encoded_log_expect
#print axioms information_discarded_identity
#print axioms information_loss_from_three_heads

end TimeFrequency.PredictiveActivity
