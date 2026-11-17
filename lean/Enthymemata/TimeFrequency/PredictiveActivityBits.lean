import Enthymemata.TimeFrequency.PredictiveActivityRetention

/-! # Joint activity loss in bits, with an encoded-law probability floor -/

namespace TimeFrequency.PredictiveActivity

open Finset hiding expect
open ContextTree.DepthEvidence

noncomputable def klNats (p q : ActivityLaw) : ℝ :=
  ∑ y, p.prob y * Real.log (p.prob y / q.prob y)

noncomputable def entropyBits (p : ActivityLaw) : ℝ :=
  -(∑ y, p.prob y * Real.log (p.prob y)) / Real.log 2

theorem scalar_kl_upper {p q : ℝ} (hp : 0 ≤ p) (hq : 0 < q) :
    p * Real.log (p / q) ≤ p - q + (p - q) ^ 2 / q := by
  by_cases hz : p = 0
  · subst p
    simp only [zero_mul]
    have he : (0 - q + (0 - q) ^ 2 / q) = 0 := by
      field_simp [hq.ne']
      ring
    rw [he]
  · have hp' : 0 < p := lt_of_le_of_ne hp (Ne.symm hz)
    have hl := mul_le_mul_of_nonneg_left
      (Real.log_le_sub_one_of_pos (div_pos hp' hq)) hp
    have he : p * (p / q - 1) = p - q + (p - q) ^ 2 / q := by
      field_simp [hq.ne']
      ring
    rwa [he] at hl

theorem kl_floor_upper (p q : ActivityLaw) {β : ℝ} (hb : 0 < β)
    (hq : ∀ y, β ≤ q.prob y) :
    klNats p q ≤ (∑ y, (p.prob y - q.prob y) ^ 2) / β := by
  have hs : (∑ y, (p.prob y - q.prob y)) = 0 := by
    rw [sum_sub_distrib, p.total, q.total, sub_self]
  calc
    _ ≤ ∑ y, (p.prob y - q.prob y + (p.prob y - q.prob y) ^ 2 / β) := by
      apply sum_le_sum
      intro y _
      exact (scalar_kl_upper (p.nonneg y) (hb.trans_le (hq y))).trans
        (add_le_add le_rfl (div_le_div_of_nonneg_left (sq_nonneg _) hb (hq y)))
    _ = _ := by rw [sum_add_distrib, hs, zero_add, sum_div]

variable {X Z : Type*} [Fintype X] [Finite Z] [DecidableEq Z]

theorem discarded_bits_bound (P : Law X) (f : X → Z)
    (hf : Function.Surjective f) (p : X → ActivityLaw) (head : Z → Fin 3 → ℝ)
    {β : ℝ} (hb : 0 < β) (hq : ∀ z y, β ≤ (encodedLaw P f hf p z).prob y) :
    expect P (fun x ↦ klNats (p x) (encodedLaw P f hf p (f x))) / Real.log 2 ≤
      8 * headExcess P f p head / (β * Real.log 2) := by
  have he : expect P (fun x ↦ klNats (p x) (encodedLaw P f hf p (f x))) ≤
      expect P (fun x ↦ ∑ y,
        ((p x).prob y - (encodedLaw P f hf p (f x)).prob y) ^ 2) / β := by
    unfold expect
    rw [sum_div]
    simp only [mul_div_assoc]
    exact sum_le_sum (fun x _ ↦ mul_le_mul_of_nonneg_left
      (kl_floor_upper _ _ hb (hq (f x))) (P.positive x).le)
  have ht := he.trans (div_le_div_of_nonneg_right (joint_retention_bound P f hf p head) hb.le)
  have hl : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hr := div_le_div_of_nonneg_right ht hl.le
  simpa only [div_div] using hr

#print axioms scalar_kl_upper
#print axioms kl_floor_upper
#print axioms discarded_bits_bound

end TimeFrequency.PredictiveActivity
