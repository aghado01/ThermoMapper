import Enthymemata.TimeFrequency.MemoryDepthConcentration

/-! # Sampling errors for bounded power summaries on separate independent views -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset MeasureTheory

variable {X J : Type*} [Fintype X] [Fintype J]

noncomputable def sampleMean {n : ℕ} (f : X → ℝ) (x : Fin n → X) :=
  (∑ i, f (x i)) / n

theorem probability_coordinate (P : Law X) (n : ℕ) (i : Fin n) (E : X → Prop) :
    probability (iid n P) (fun x ↦ E (x i)) = probability P E := by
  classical
  unfold probability
  exact expect_coordinate n P i (fun x ↦ if E x then 1 else 0)

theorem empirical_direction [MeasurableSpace X] [MeasurableSingletonClass X]
    (P : Law X) (f : X → ℝ) (hf : ∀ x, f x ∈ Set.Icc (0 : ℝ) 1)
    (n : ℕ) (hn : 0 < n) {ε s : ℝ} (he : 0 < ε) (hs : s ^ 2 = 1) :
    probability (iid n P) (fun x ↦ ε ≤ s * (sampleMean f x - expect P f)) ≤
      Real.exp (-2 * n * ε ^ 2) := by
  let g := fun x ↦ s * (f x - expect P f)
  have hc := finite_chernoff P g n (b := n * ε) (show 0 < 4 * ε by positivity)
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have hsum (x : Fin n → X) :
      (∑ i, g (x i)) = n * (s * (sampleMean f x - expect P f)) := by
    simp only [g, sampleMean, mul_sub, sum_sub_distrib, ← mul_sum, sum_const,
      card_univ, Fintype.card_fin, nsmul_eq_mul]
    field_simp [hnR.ne']
  have hev (x : Fin n → X) :
      ((n : ℝ) * ε ≤ n * (s * (sampleMean f x - expect P f))) ↔
        ε ≤ s * (sampleMean f x - expect P f) := by
    constructor <;> intro hx <;> nlinarith
  simp_rw [hsum, hev] at hc
  have hm := bounded_mgf P f hf (4 * ε * s)
  have hexp : (4 * ε * s) ^ 2 / 8 = 2 * ε ^ 2 := by nlinarith [hs]
  rw [hexp] at hm
  have hefun : (fun x ↦ Real.exp (4 * ε * g x)) =
      (fun x ↦ Real.exp (4 * ε * s * (f x - expect P f))) := by
    funext x
    congr 1
    dsimp [g]
    ring
  rw [← hefun] at hm
  have hn0 : 0 ≤ expect P (fun x ↦ Real.exp (4 * ε * g x)) :=
    sum_nonneg (fun _ _ ↦ mul_nonneg (P.positive _).le (Real.exp_pos _).le)
  have hp := pow_le_pow_left₀ hn0 hm n
  calc _ ≤ Real.exp (-(4 * ε) * (n * ε)) *
        (expect P (fun x ↦ Real.exp (4 * ε * g x))) ^ n := hc
    _ ≤ Real.exp (-(4 * ε) * (n * ε)) * (Real.exp (2 * ε ^ 2)) ^ n :=
      mul_le_mul_of_nonneg_left hp (Real.exp_pos _).le
    _ = _ := by rw [← Real.exp_nat_mul, ← Real.exp_add]; congr 1; ring

theorem empirical_error [MeasurableSpace X] [MeasurableSingletonClass X]
    (P : Law X) (f : X → ℝ) (hf : ∀ x, f x ∈ Set.Icc (0 : ℝ) 1)
    (n : ℕ) (hn : 0 < n) {ε : ℝ} (he : 0 < ε) :
    probability (iid n P) (fun x ↦ ε ≤ |sampleMean f x - expect P f|) ≤
      2 * Real.exp (-2 * n * ε ^ 2) := by
  have hU := empirical_direction P f hf n hn he (show (1 : ℝ) ^ 2 = 1 by norm_num)
  have hL := empirical_direction P f hf n hn he (show (-1 : ℝ) ^ 2 = 1 by norm_num)
  simp only [one_mul] at hU
  simp only [neg_one_mul] at hL
  have hp := probability_union (iid n P)
    (fun x ↦ ε ≤ sampleMean f x - expect P f)
    (fun x ↦ ε ≤ -(sampleMean f x - expect P f))
  have hd : (fun x : Fin n → X ↦ ε ≤ |sampleMean f x - expect P f|) =
      (fun x ↦ ε ≤ sampleMean f x - expect P f ∨
        ε ≤ -(sampleMean f x - expect P f)) := by
    funext x
    exact propext le_abs
  rw [hd]
  linarith

theorem separate_summary_error [MeasurableSpace X] [MeasurableSingletonClass X]
    (P : Law X) (f : J → X → ℝ) (hf : ∀ j x, f j x ∈ Set.Icc (0 : ℝ) 1)
    (r n : ℕ) (view : J → Fin r) (hn : 0 < n) {ε : ℝ} (he : 0 < ε) :
    probability (iid r (iid n P)) (fun x ↦ ∃ j,
      ε ≤ |sampleMean (f j) (x (view j)) - expect P (f j)|) ≤
      2 * Fintype.card J * Real.exp (-2 * n * ε ^ 2) := by
  have hb (j : J) : probability (iid r (iid n P))
      (fun x ↦ ε ≤ |sampleMean (f j) (x (view j)) - expect P (f j)|) ≤
      2 * Real.exp (-2 * n * ε ^ 2) := by
    rw [probability_coordinate (iid n P) r (view j)
      (fun y ↦ ε ≤ |sampleMean (f j) y - expect P (f j)|)]
    exact empirical_error P (f j) (hf j) n hn he
  have hs := sum_le_sum (fun j (_ : j ∈ univ) ↦ hb j)
  have hu := probability_exists (iid r (iid n P))
    (fun j x ↦ ε ≤ |sampleMean (f j) (x (view j)) - expect P (f j)|)
  simp only [sum_const, card_univ, nsmul_eq_mul] at hs
  exact hu.trans (by nlinarith only [hs])

theorem exponential_budget {c δ rate : ℝ} (hc : 0 < c) (hδ : 0 < δ)
    (hr : 0 < rate) (n : ℕ) (hn : ⌈Real.log (c / δ) / rate⌉₊ ≤ n) :
    c * Real.exp (-(n : ℝ) * rate) ≤ δ := by
  have hnR := (Nat.le_ceil (Real.log (c / δ) / rate)).trans (Nat.cast_le.mpr hn)
  have hm := (div_le_iff₀ hr).mp hnR
  have he : Real.exp (-(n : ℝ) * rate) ≤ δ / c := by
    apply (Real.le_log_iff_exp_le (div_pos hδ hc)).mp
    rw [Real.log_div hc.ne' hδ.ne'] at hm
    rw [Real.log_div hδ.ne' hc.ne']
    linarith
  have h := mul_le_mul_of_nonneg_left he hc.le
  have heq : c * (δ / c) = δ := by field_simp
  rwa [heq] at h

theorem separate_summary_budget [MeasurableSpace X] [MeasurableSingletonClass X]
    (P : Law X) (f : J → X → ℝ) (hf : ∀ j x, f j x ∈ Set.Icc (0 : ℝ) 1)
    (r n : ℕ) (view : J → Fin r) {ε δ : ℝ} (he : 0 < ε)
    (hδ : 0 < δ) (hδu : δ < 1) (hJ : 0 < Fintype.card J)
    (hn : ⌈Real.log (2 * Fintype.card J / δ) / (2 * ε ^ 2)⌉₊ ≤ n) :
    probability (iid r (iid n P)) (fun x ↦ ∃ j,
      ε ≤ |sampleMean (f j) (x (view j)) - expect P (f j)|) ≤ δ := by
  have hc : (1 : ℝ) ≤ 2 * Fintype.card J := by
    exact_mod_cast (by omega : 1 ≤ 2 * Fintype.card J)
  have hl : 0 < Real.log (2 * Fintype.card J / δ) :=
    Real.log_pos ((lt_div_iff₀ hδ).mpr (by linarith))
  have hnR := (Nat.le_ceil (Real.log (2 * Fintype.card J / δ) / (2 * ε ^ 2))).trans
    (Nat.cast_le.mpr hn)
  have hn0 : 0 < n := by
    have hpos : (0 : ℝ) < n := (div_pos hl (by positivity)).trans_le hnR
    exact_mod_cast hpos
  have hb := separate_summary_error P f hf r n view hn0 he
  have heq : -2 * (n : ℝ) * ε ^ 2 = -(n : ℝ) * (2 * ε ^ 2) := by ring
  rw [heq] at hb
  exact hb.trans (exponential_budget (by linarith) hδ (by positivity) n hn)

#print axioms separate_summary_error
#print axioms exponential_budget
#print axioms separate_summary_budget

#print axioms probability_coordinate
#print axioms empirical_direction
#print axioms empirical_error

end TimeFrequency.ContextTree.DepthEvidence
