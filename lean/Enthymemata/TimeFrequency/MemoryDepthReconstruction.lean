import Enthymemata.TimeFrequency.MemoryDepthConfidence
import Enthymemata.TimeFrequency.DepthFourCompletion

/-! # Explicit depth-four inversion from eight normalized power summaries -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset DepthFour

private theorem sum_five (f : Fin 5 → ℝ) :
    ∑ i, f i = f 0 + f 1 + f 2 + f 3 + f 4 := by
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]
  change f 0 + (f 1 + (f 2 + (f 3 + f 4))) = _
  ring

private theorem moment_add {n : ℕ} (p : Word n → ℝ) (f g : Word n → ℤ) :
    mean p (fun x ↦ f x + g x) = mean p f + mean p g := by
  simp [mean, Int.cast_add, mul_add, sum_add_distrib]

private theorem moment_scale {n : ℕ} (p : Word n → ℝ) (c : ℤ) (f : Word n → ℤ) :
    mean p (fun x ↦ c * f x) = c * mean p f := by
  simp only [mean, Int.cast_mul, mul_sum]
  apply sum_congr rfl
  intros
  ring

private theorem moment_sum {n : ℕ} {I : Type*} [Fintype I]
    (p : Word n → ℝ) (f : I → Word n → ℤ) :
    mean p (fun x ↦ ∑ i, f i x) = ∑ i, mean p (f i) := by
  simp only [mean, Int.cast_sum, mul_sum]
  rw [sum_comm]

private theorem moment_one (P : Source) (n : ℕ) : mean (P.law n) (fun _ ↦ 1) = 1 := by
  simp [mean, P.total]

private theorem moment_event {n : ℕ} (p : Word n → ℝ) (f : Word n → ℤ) (v : ℤ) :
    mean p (fun x ↦ hit (f x) v) = bin p f v := by
  simp [mean, hit, bin, mul_ite]

private theorem moment_atom {n : ℕ} (p : Word n → ℝ) (t : Word n) :
    mean p (fun x ↦ hit x t) = p t := by
  classical
  simp [mean, hit, mul_ite]

def lagCoefficients (t : Word 5) : Fin 4 → ℤ :=
  ![sign (t 0) * sign (t 1) + sign (t 1) * sign (t 2) +
      sign (t 2) * sign (t 3) + sign (t 3) * sign (t 4),
    sign (t 0) * sign (t 2) + sign (t 1) * sign (t 3) + sign (t 2) * sign (t 4),
    sign (t 0) * sign (t 3) + sign (t 1) * sign (t 4),
    sign (t 0) * sign (t 4)]

theorem pair_moments (P : Source) (t : Word 5) :
    mean (P.law 5) (pairSum t) =
      (lagCoefficients t 0 : ℝ) * lag P 1 + (lagCoefficients t 1 : ℝ) * lag P 2 +
      (lagCoefficients t 2 : ℝ) * lag P 3 + (lagCoefficients t 3 : ℝ) * lag P 4 := by
  have hf : pairSum t = fun x ↦ (sign (t 0) * sign (t 1)) * (sign (x 0) * sign (x 1)) +
      (sign (t 0) * sign (t 2)) * (sign (x 0) * sign (x 2)) +
      (sign (t 0) * sign (t 3)) * (sign (x 0) * sign (x 3)) +
      (sign (t 0) * sign (t 4)) * (sign (x 0) * sign (x 4)) +
      (sign (t 1) * sign (t 2)) * (sign (x 1) * sign (x 2)) +
      (sign (t 1) * sign (t 3)) * (sign (x 1) * sign (x 3)) +
      (sign (t 1) * sign (t 4)) * (sign (x 1) * sign (x 4)) +
      (sign (t 2) * sign (t 3)) * (sign (x 2) * sign (x 3)) +
      (sign (t 2) * sign (t 4)) * (sign (x 2) * sign (x 4)) +
      (sign (t 3) * sign (t 4)) * (sign (x 3) * sign (x 4)) := by
    funext x
    decide +revert
  rw [hf]
  simp only [moment_add, moment_scale, mean_pair P (0 : Fin 5) 1 (by decide),
    mean_pair P (0 : Fin 5) 2 (by decide),
    mean_pair P (0 : Fin 5) 3 (by decide),
    mean_pair P (0 : Fin 5) 4 (by decide),
    mean_pair P (1 : Fin 5) 2 (by decide),
    mean_pair P (1 : Fin 5) 3 (by decide),
    mean_pair P (1 : Fin 5) 4 (by decide),
    mean_pair P (2 : Fin 5) 3 (by decide),
    mean_pair P (2 : Fin 5) 4 (by decide),
    mean_pair P (3 : Fin 5) 4 (by decide)]
  dsimp [lagCoefficients]
  push_cast
  ring

theorem selected_pair_moments (P : Source) (r : Fin 5) :
    mean (P.law 5) (selectedPairs r) =
      ∑ i : Fin 5, ∑ j : Fin 5,
        if i < j ∧ i ≠ r ∧ j ≠ r then lag P (j.val - i.val) else 0 := by
  have he (i j : Fin 5) : mean (P.law 5) (fun x ↦
      if i < j ∧ i ≠ r ∧ j ≠ r then sign (x i) * sign (x j) else 0) =
      if i < j ∧ i ≠ r ∧ j ≠ r then lag P (j.val - i.val) else 0 := by
    by_cases hij : i < j ∧ i ≠ r ∧ j ≠ r
    · simp only [if_pos hij, mean_pair P i j hij.1]
    · simp [hij, mean]
  change mean (P.law 5) (fun x ↦ ∑ i, ∑ j,
    if i < j ∧ i ≠ r ∧ j ≠ r then sign (x i) * sign (x j) else 0) = _
  rw [moment_sum]
  simp_rw [moment_sum, he]

theorem taper_event_moment (P : Source) (r : Fin 5) :
    8 * bin (P.law 5) (taper r) 16 =
      1 + mean (P.law 5) (selectedPairs r) + mean (P.law 5) (quartic r) := by
  have he := congrArg (mean (P.law 5)) (funext (tapered_fourth r))
  simpa only [moment_scale, moment_event, moment_add, moment_one, Int.cast_ofNat] using he

theorem top_event_moment (P : Source) :
    16 * bin (P.law 5) dc 25 =
      1 + mean (P.law 5) (pairSum (fun _ ↦ true)) + ∑ r, mean (P.law 5) (quartic r) := by
  have he := congrArg (mean (P.law 5)) (funext five_top)
  simpa only [moment_scale, moment_event, moment_add, moment_one, moment_sum,
    Int.cast_ofNat] using he

theorem endpoint_moments (P : Source) :
    mean (P.law 5) (quartic 0) = mean (P.law 4) character ∧
    mean (P.law 5) (quartic 4) = mean (P.law 4) character := by
  have h0 : quartic 0 = fun x ↦ character (slice x 1 4 (by decide)) :=
    funext (fun x ↦ (endpoint_quartics x).2)
  have h4 : quartic 4 = fun x ↦ character (slice x 0 4 (by decide)) :=
    funext (fun x ↦ (endpoint_quartics x).1)
  rw [h0, h4, P.translation, P.translation]
  exact ⟨rfl, rfl⟩

noncomputable def powerSummaries (P : Source) : Fin 8 → ℝ :=
  ![powerMean P 2 / 4, powerMean P 3 / 9, powerMean P 4 / 16, powerMean P 5 / 25,
    bin (P.law 4) dc 16, bin (P.law 5) dc 25,
    bin (P.law 5) (taper 3) 16, bin (P.law 5) (taper 2) 16]

def decodePower {R : Type*} [Field R] (m : Fin 8 → R) : Fin 8 → R :=
  let r1 := (4 * m 0 - 2) / 2
  let r2 := (9 * m 1 - 8 * m 0 + 1) / 2
  let r3 := (16 * m 2 - 18 * m 1 + 4 * m 0) / 2
  let r4 := (25 * m 3 - 32 * m 2 + 9 * m 1) / 2
  let q := 8 * m 4 - 1 - 3 * r1 - 2 * r2 - r3
  let a := 8 * m 6 - 1 - 2 * r1 - 2 * r2 - r3 - r4
  let b := 8 * m 7 - 1 - 2 * r1 - r2 - 2 * r3 - r4
  let c := 16 * m 5 - 1 - 4 * r1 - 3 * r2 - 2 * r3 - r4 - 2 * q - a - b
  ![r1, r2, r3, r4, q, a, b, c]

def reconstructPower {R : Type*} [Field R] (m : Fin 8 → R) (t : Word 5) : R :=
  let r := decodePower m
  (1 + (lagCoefficients t 0 : R) * r 0 + (lagCoefficients t 1 : R) * r 1 +
    (lagCoefficients t 2 : R) * r 2 + (lagCoefficients t 3 : R) * r 3 +
    ((quartic 0 t : R) + (quartic 4 t : R)) * r 4 + (quartic 3 t : R) * r 5 +
    (quartic 2 t : R) * r 6 + (quartic 1 t : R) * r 7) / 32

theorem decoded_moments (P : Source) : decodePower (powerSummaries P) =
    ![lag P 1, lag P 2, lag P 3, lag P 4, mean (P.law 4) character,
      mean (P.law 5) (quartic 3), mean (P.law 5) (quartic 2),
      mean (P.law 5) (quartic 1)] := by
  have h2 := power_mean_formula P 2
  have h3 := power_mean_formula P 3
  have h4 := power_mean_formula P 4
  have h5 := power_mean_formula P 5
  norm_num [sum_range_succ] at h2 h3 h4 h5
  have hq := fourth_from_dc P
  have ha := taper_event_moment P 3
  have hb := taper_event_moment P 2
  have hs3 : mean (P.law 5) (selectedPairs 3) =
      2 * lag P 1 + 2 * lag P 2 + lag P 3 + lag P 4 := by
    have hf : selectedPairs 3 = fun x ↦ (sign (x 0) * sign (x 1)) +
        (sign (x 0) * sign (x 2)) +
        (sign (x 0) * sign (x 4)) +
        (sign (x 1) * sign (x 2)) +
        (sign (x 1) * sign (x 4)) +
        (sign (x 2) * sign (x 4)) := by
      funext x
      decide +revert
    rw [hf]
    simp only [moment_add, mean_pair P (0 : Fin 5) 1 (by decide),
      mean_pair P (0 : Fin 5) 2 (by decide),
      mean_pair P (0 : Fin 5) 4 (by decide),
      mean_pair P (1 : Fin 5) 2 (by decide),
      mean_pair P (1 : Fin 5) 4 (by decide),
      mean_pair P (2 : Fin 5) 4 (by decide)]
    change lag P 1 + lag P 2 + lag P 4 + lag P 1 + lag P 3 + lag P 2 = _
    ring
  have hs2 : mean (P.law 5) (selectedPairs 2) =
      2 * lag P 1 + lag P 2 + 2 * lag P 3 + lag P 4 := by
    have hf : selectedPairs 2 = fun x ↦ (sign (x 0) * sign (x 1)) +
        (sign (x 0) * sign (x 3)) +
        (sign (x 0) * sign (x 4)) +
        (sign (x 1) * sign (x 3)) +
        (sign (x 1) * sign (x 4)) +
        (sign (x 3) * sign (x 4)) := by
      funext x
      decide +revert
    rw [hf]
    simp only [moment_add, mean_pair P (0 : Fin 5) 1 (by decide),
      mean_pair P (0 : Fin 5) 3 (by decide),
      mean_pair P (0 : Fin 5) 4 (by decide),
      mean_pair P (1 : Fin 5) 3 (by decide),
      mean_pair P (1 : Fin 5) 4 (by decide),
      mean_pair P (3 : Fin 5) 4 (by decide)]
    change lag P 1 + lag P 3 + lag P 4 + lag P 2 + lag P 3 + lag P 1 = _
    ring
  rw [hs3] at ha
  rw [hs2] at hb
  have hc := top_event_moment P
  rw [pair_moments] at hc
  rw [sum_five] at hc
  dsimp [lagCoefficients, sign] at hc
  norm_num only [Int.cast_ofNat, Int.cast_one] at hc
  rw [(endpoint_moments P).1, (endpoint_moments P).2] at hc
  funext i
  fin_cases i <;> simp [decodePower, powerSummaries] <;> linarith

theorem five_atom_moments (P : Source) (t : Word 5) :
    32 * P.law 5 t =
      1 + mean (P.law 5) (pairSum t) +
        ∑ r, (quartic r t : ℝ) * mean (P.law 5) (quartic r) := by
  have he := congrArg (mean (P.law 5)) (funext (five_inversion t))
  simp only [moment_scale, moment_add, moment_atom, P.symmetry, moment_one, moment_sum,
    Int.cast_ofNat] at he
  linarith

theorem reconstruct_power_law (P : Source) (t : Word 5) :
    reconstructPower (powerSummaries P) t = P.law 5 t := by
  have he := five_atom_moments P t
  rw [pair_moments] at he
  rw [sum_five] at he
  rw [(endpoint_moments P).1, (endpoint_moments P).2] at he
  simp only [reconstructPower, decoded_moments]
  change (1 + (lagCoefficients t 0 : ℝ) * lag P 1 +
    (lagCoefficients t 1 : ℝ) * lag P 2 + (lagCoefficients t 2 : ℝ) * lag P 3 +
    (lagCoefficients t 3 : ℝ) * lag P 4 +
    ((quartic 0 t : ℝ) + (quartic 4 t : ℝ)) * mean (P.law 4) character +
    (quartic 3 t : ℝ) * mean (P.law 5) (quartic 3) +
    (quartic 2 t : ℝ) * mean (P.law 5) (quartic 2) +
    (quartic 1 t : ℝ) * mean (P.law 5) (quartic 1)) / 32 = _
  linarith

#print axioms pair_moments
#print axioms selected_pair_moments
#print axioms taper_event_moment
#print axioms top_event_moment
#print axioms endpoint_moments
#print axioms decoded_moments
#print axioms five_atom_moments
#print axioms reconstruct_power_law

end TimeFrequency.ContextTree.DepthEvidence
