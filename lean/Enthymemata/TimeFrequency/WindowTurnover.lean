import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Complex.ExponentialBounds
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity

/-!
# The rectangular window's turnover (S13)

A rectangular window of length `w` on a resonance decaying as `e^{−j/τ}` reads power
`τ·φ(w/τ)` at zero offset, with `φ(u) = (1 − e^{−u})²/u`. The turnover in window length is
`φ`'s maximizer:

* `turnover_spec`: `φ` increases on `(0, u⋆]` and decreases on `[u⋆, ∞)`. Here `u⋆` is the unique
  positive root of `eᵘ = 1 + 2u`, and `1.25 < u⋆ < 1.26`. The numerical chip read `1.255·τ`.

`φ′(u) = (1 − e^{−u})·(e^{−u}(1 + 2u) − 1)/u²`. The middle factor has the sign of
`ψ(u) = 1 + 2u − eᵘ`, which is zero at `0`, increases to `ln 2` and decreases after, so it has one
positive root. The bracket comes from `e^{5/4} < 3.5` and `e^{1.26} > 3.52`, by five Taylor terms
and Mathlib's nine-digit bounds on `e`.
-/

open Real Set

namespace TimeFrequency

/-- `φ(u) = (1 − e^{−u})²/u`. -/
noncomputable def turnover (u : ℝ) : ℝ :=
  (1 - exp (-u)) ^ 2 / u

/-- `ψ(u) = 1 + 2u − eᵘ`, whose sign is `φ′`'s. -/
noncomputable def turnoverGap (u : ℝ) : ℝ :=
  1 + 2 * u - exp u

theorem hasDerivAt_turnoverGap (u : ℝ) : HasDerivAt turnoverGap (2 - exp u) u := by
  have h := (((hasDerivAt_id u).const_mul 2).const_add 1).sub (hasDerivAt_exp u)
  exact h.congr_deriv (by simp)

theorem hasDerivAt_turnover {u : ℝ} (hu : u ≠ 0) :
    HasDerivAt turnover ((1 - exp (-u)) * (exp (-u) * (1 + 2 * u) - 1) / u ^ 2) u := by
  have hq : HasDerivAt (fun x => 1 - exp (-x)) (exp (-u)) u := by
    simpa using ((hasDerivAt_neg u).exp).const_sub 1
  have h := (hq.pow 2).div (hasDerivAt_id u) hu
  refine h.congr_deriv ?_
  simp only [id]
  norm_num
  ring

/-- `e^{−u}(1 + 2u) − 1 = e^{−u}·ψ(u)`. -/
theorem exp_neg_mul_sub_one (u : ℝ) :
    exp (-u) * (1 + 2 * u) - 1 = exp (-u) * turnoverGap u := by
  have h : exp (-u) * exp u = 1 := by rw [← exp_add, neg_add_cancel, exp_zero]
  unfold turnoverGap
  linear_combination h

theorem exp_quarter_lt : exp (1 / 4) < 1.2841 := by
  have h := Real.exp_bound (x := 1 / 4) (by rw [abs_of_pos (by norm_num)]; norm_num)
    (n := 5) (by norm_num)
  rw [abs_le] at h
  norm_num [Finset.sum_range_succ, Nat.factorial] at h
  linarith [h.2]

theorem exp_thirteen_fiftieths_gt : 1.2969 < exp (13 / 50) := by
  have h := Real.exp_bound (x := 13 / 50) (by rw [abs_of_pos (by norm_num)]; norm_num)
    (n := 5) (by norm_num)
  rw [abs_le] at h
  norm_num [Finset.sum_range_succ, Nat.factorial] at h
  linarith [h.1]

theorem turnoverGap_five_quarters_pos : 0 < turnoverGap (5 / 4) := by
  have he : exp (5 / 4) = exp 1 * exp (1 / 4) := by rw [← exp_add]; norm_num
  have h1 := exp_one_lt_d9
  have h2 := exp_quarter_lt
  have h3 : exp 1 * exp (1 / 4) < 2.7182818286 * 1.2841 :=
    mul_lt_mul'' h1 h2 (exp_pos _).le (exp_pos _).le
  unfold turnoverGap
  rw [he]
  norm_num at h3 ⊢
  linarith

theorem turnoverGap_sixty_three_fiftieths_neg : turnoverGap (63 / 50) < 0 := by
  have he : exp (63 / 50) = exp 1 * exp (13 / 50) := by rw [← exp_add]; norm_num
  have h1 := exp_one_gt_d9
  have h2 := exp_thirteen_fiftieths_gt
  have h3 : 2.7182818283 * 1.2969 < exp 1 * exp (13 / 50) :=
    mul_lt_mul'' h1 h2 (by norm_num) (by norm_num)
  unfold turnoverGap
  rw [he]
  norm_num at h3 ⊢
  linarith

theorem turnoverGap_strictMonoOn : StrictMonoOn turnoverGap (Icc 0 (log 2)) := by
  refine strictMonoOn_of_deriv_pos (convex_Icc 0 (log 2))
    (fun x _ => (hasDerivAt_turnoverGap x).continuousAt.continuousWithinAt) fun x hx => ?_
  rw [interior_Icc] at hx
  rw [(hasDerivAt_turnoverGap x).deriv, sub_pos]
  calc exp x < exp (log 2) := exp_lt_exp.2 hx.2
    _ = 2 := exp_log (by norm_num)

theorem turnoverGap_strictAntiOn : StrictAntiOn turnoverGap (Ici (log 2)) := by
  refine strictAntiOn_of_deriv_neg (convex_Ici (log 2))
    (fun x _ => (hasDerivAt_turnoverGap x).continuousAt.continuousWithinAt) fun x hx => ?_
  rw [interior_Ici] at hx
  rw [(hasDerivAt_turnoverGap x).deriv, sub_neg]
  calc (2 : ℝ) = exp (log 2) := (exp_log (by norm_num)).symm
    _ < exp x := exp_lt_exp.2 hx

theorem turnoverGap_zero : turnoverGap 0 = 0 := by
  simp [turnoverGap]

/-- **S13.** `φ(u) = (1 − e^{−u})²/u` increases on `(0, u⋆]` and decreases on `[u⋆, ∞)`, where
`u⋆` is the unique positive root of `eᵘ = 1 + 2u`, and `1.25 < u⋆ < 1.26`. -/
theorem turnover_spec : ∃ u₀ : ℝ, 0 < u₀ ∧ exp u₀ = 1 + 2 * u₀ ∧
    (∀ u, 0 < u → exp u = 1 + 2 * u → u = u₀) ∧ 5 / 4 < u₀ ∧ u₀ < 63 / 50 ∧
    StrictMonoOn turnover (Ioc 0 u₀) ∧ StrictAntiOn turnover (Ici u₀) := by
  have hcont : Continuous turnoverGap :=
    continuous_iff_continuousAt.2 fun x => (hasDerivAt_turnoverGap x).continuousAt
  -- a root between 1.25 and 1.26
  obtain ⟨u₀, ⟨hlo, hhi⟩, hroot⟩ : ∃ u₀ ∈ Ioo (5 / 4 : ℝ) (63 / 50), turnoverGap u₀ = 0 :=
    intermediate_value_Ioo' (by norm_num) hcont.continuousOn
      ⟨turnoverGap_sixty_three_fiftieths_neg, turnoverGap_five_quarters_pos⟩
  have hlog : log 2 < 5 / 4 := by linarith [log_two_lt_d9]
  have hu₀ : log 2 < u₀ := hlog.trans hlo
  have hlog0 : 0 < log 2 := log_pos (by norm_num)
  -- the sign of `ψ` on either side of `u₀`
  have hpos : ∀ u ∈ Ioo 0 u₀, 0 < turnoverGap u := by
    intro u hu
    by_cases h : u ≤ log 2
    · have := turnoverGap_strictMonoOn ⟨le_rfl, hlog0.le⟩ ⟨hu.1.le, h⟩ hu.1
      rwa [turnoverGap_zero] at this
    · have := turnoverGap_strictAntiOn (not_le.1 h).le hu₀.le hu.2
      rwa [hroot] at this
  have hneg : ∀ u, u₀ < u → turnoverGap u < 0 := fun u hu => by
    have := turnoverGap_strictAntiOn hu₀.le (hu₀.trans hu).le hu
    rwa [hroot] at this
  have hu₀pos : 0 < u₀ := hlog0.trans hu₀
  refine ⟨u₀, hu₀pos, ?_, fun u hu he => ?_, hlo, hhi, ?_, ?_⟩
  · unfold turnoverGap at hroot
    linarith
  · have hgap : turnoverGap u = 0 := by unfold turnoverGap; linarith
    by_contra hne
    rcases lt_or_gt_of_ne hne with h | h
    · exact (hpos u ⟨hu, h⟩).ne' hgap
    · exact (hneg u h).ne hgap
  · refine strictMonoOn_of_deriv_pos (convex_Ioc 0 u₀)
      (fun x hx => (hasDerivAt_turnover hx.1.ne').continuousAt.continuousWithinAt)
      fun x hx => ?_
    rw [interior_Ioc] at hx
    rw [(hasDerivAt_turnover hx.1.ne').deriv, exp_neg_mul_sub_one]
    have h1 : 0 < 1 - exp (-x) := sub_pos.2 (exp_lt_one_iff.2 (neg_lt_zero.2 hx.1))
    have h2 : 0 < exp (-x) * turnoverGap x := mul_pos (exp_pos _) (hpos x hx)
    exact div_pos (mul_pos h1 h2) (pow_pos hx.1 2)
  · refine strictAntiOn_of_deriv_neg (convex_Ici u₀)
      (fun x hx => (hasDerivAt_turnover (hu₀pos.trans_le hx).ne').continuousAt.continuousWithinAt)
      fun x hx => ?_
    rw [interior_Ici] at hx
    have hx0 : 0 < x := hu₀pos.trans hx
    rw [(hasDerivAt_turnover hx0.ne').deriv, exp_neg_mul_sub_one]
    have h1 : 0 < 1 - exp (-x) := sub_pos.2 (exp_lt_one_iff.2 (neg_lt_zero.2 hx0))
    have h2 : exp (-x) * turnoverGap x < 0 := mul_neg_of_pos_of_neg (exp_pos _) (hneg x hx)
    exact div_neg_of_neg_of_pos (mul_neg_of_pos_of_neg h1 h2) (by positivity)

/-- The maximizer: `φ(u) ≤ φ(u⋆)` for every `u > 0`. -/
theorem turnover_le : ∃ u₀ : ℝ, 5 / 4 < u₀ ∧ u₀ < 63 / 50 ∧
    ∀ u, 0 < u → turnover u ≤ turnover u₀ := by
  obtain ⟨u₀, hu₀, -, -, hlo, hhi, hmono, hanti⟩ := turnover_spec
  refine ⟨u₀, hlo, hhi, fun u hu => ?_⟩
  rcases le_total u u₀ with h | h
  · exact hmono.monotoneOn ⟨hu, h⟩ ⟨hu₀, le_rfl⟩ h
  · exact hanti.antitoneOn self_mem_Ici h h

end TimeFrequency
