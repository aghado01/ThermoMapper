import Enthymemata.TimeFrequency.MemoryDepthReconstruction

/-! # Quantified error propagation through the depth-four inversion -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset DepthFour

def reconstructionBase (t : Word 5) : ℚ := reconstructPower (fun _ ↦ 0) t
def reconstructionRow (t : Word 5) (j : Fin 8) : ℚ :=
  reconstructPower (fun k ↦ if k = j then 1 else 0) t - reconstructionBase t

def reconstructionIntegerRow (t : Word 5) (j : Fin 8) : ℤ :=
  let l := lagCoefficients t
  let q := quartic 0 t + quartic 4 t
  let a := quartic 3 t
  let b := quartic 2 t
  let c := quartic 1 t
  ![4 * l 0 - 8 * l 1 + 4 * l 2 + 4 * a - 8 * b + 4 * c,
    9 * l 1 - 18 * l 2 + 9 * l 3 - 9 * a + 18 * b - 9 * c,
    16 * l 2 - 32 * l 3 - 16 * q + 16 * a + 16 * c,
    25 * l 3 - 25 * a - 25 * b + 25 * c,
    16 * q - 32 * c, 32 * c, 16 * a - 16 * c, 16 * b - 16 * c] j

theorem reconstruction_integer_bound (t : Word 5) :
    ∑ j, (reconstructionIntegerRow t j).natAbs ≤ 1540 := by
  decide +revert

theorem reconstruction_row_integer (t : Word 5) (j : Fin 8) :
    reconstructionRow t j = (reconstructionIntegerRow t j : ℚ) / 64 := by
  fin_cases j <;>
    simp [reconstructionRow, reconstructionBase, reconstructPower, decodePower,
      reconstructionIntegerRow] <;> ring

theorem reconstruction_row_bound (t : Word 5) :
    ∑ j, |reconstructionRow t j| ≤ (385 / 16 : ℚ) := by
  have hI : (∑ j, |reconstructionIntegerRow t j| : ℤ) ≤ 1540 := by
    have hc : ((∑ j, (reconstructionIntegerRow t j).natAbs : ℕ) : ℤ) ≤ 1540 := by
      exact_mod_cast reconstruction_integer_bound t
    simpa only [Nat.cast_sum, Int.natCast_natAbs] using hc
  have hQ : (∑ j, |(reconstructionIntegerRow t j : ℚ)|) ≤ 1540 := by
    exact_mod_cast hI
  simp_rw [reconstruction_row_integer, abs_div]
  simp only [abs_of_pos (by norm_num : (0 : ℚ) < 64)]
  rw [← sum_div]
  apply (div_le_iff₀ (by norm_num : (0 : ℚ) < 64)).mpr
  norm_num
  exact hQ

theorem reconstruction_affine (m : Fin 8 → ℝ) (t : Word 5) :
    reconstructPower m t = (reconstructionBase t : ℝ) +
      ∑ j, (reconstructionRow t j : ℝ) * m j := by
  simp [reconstructionBase, reconstructionRow, reconstructPower, decodePower,
    Fin.sum_univ_succ]
  ring

theorem reconstruction_stability {u v : Fin 8 → ℝ} {ε : ℝ} (he : 0 ≤ ε)
    (hm : ∀ j, |u j - v j| ≤ ε) (t : Word 5) :
    |reconstructPower u t - reconstructPower v t| ≤ (385 / 16 : ℝ) * ε := by
  have hr : (∑ j, |(reconstructionRow t j : ℝ)|) ≤ 385 / 16 := by
    have hq : (16 : ℚ) * (∑ j, |reconstructionRow t j|) ≤ 385 := by
      linarith [reconstruction_row_bound t]
    have hR : (16 : ℝ) * (∑ j, |(reconstructionRow t j : ℝ)|) ≤ 385 := by
      exact_mod_cast hq
    linarith
  rw [reconstruction_affine, reconstruction_affine]
  have hd : (reconstructionBase t : ℝ) + ∑ j, (reconstructionRow t j : ℝ) * u j -
      ((reconstructionBase t : ℝ) + ∑ j, (reconstructionRow t j : ℝ) * v j) =
      ∑ j, (reconstructionRow t j : ℝ) * (u j - v j) := by
    simp only [mul_sub, sum_sub_distrib]
    ring
  rw [hd]
  calc _ ≤ ∑ j, |(reconstructionRow t j : ℝ) * (u j - v j)| := abs_sum_le_sum_abs _ _
    _ ≤ ∑ j, |(reconstructionRow t j : ℝ)| * ε := by
      apply sum_le_sum
      intro j _
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_left (hm j) (abs_nonneg _)
    _ = (∑ j, |(reconstructionRow t j : ℝ)|) * ε := by rw [sum_mul]
    _ ≤ _ := mul_le_mul_of_nonneg_right hr he

theorem source_block_stability (P Q : Source) {ε : ℝ} (he : 0 ≤ ε)
    (hm : ∀ j, |powerSummaries P j - powerSummaries Q j| ≤ ε) (t : Word 5) :
    |P.law 5 t - Q.law 5 t| ≤ (385 / 16 : ℝ) * ε := by
  simpa only [reconstruct_power_law] using reconstruction_stability he hm t

theorem binary_ratio_stability {p0 p1 q0 q1 η ε : ℝ}
    (hp0 : 0 ≤ p0) (hp1 : 0 ≤ p1)
    (hη : 0 < η) (hp : η ≤ p0 + p1) (hq : η ≤ q0 + q1)
    (he0 : |p0 - q0| ≤ ε) (he1 : |p1 - q1| ≤ ε) :
    |p1 / (p0 + p1) - q1 / (q0 + q1)| ≤ ε / η := by
  have hpp : 0 < p0 + p1 := hη.trans_le hp
  have hqq : 0 < q0 + q1 := hη.trans_le hq
  have he : 0 ≤ ε := (abs_nonneg _).trans he0
  have h0 := abs_le.mp he0
  have h1 := abs_le.mp he1
  have hu : p1 * q0 - q1 * p0 ≤ ε * (p0 + p1) := by
    nlinarith [mul_le_mul_of_nonneg_left h0.1 hp1, mul_le_mul_of_nonneg_left h1.2 hp0]
  have hl : -ε * (p0 + p1) ≤ p1 * q0 - q1 * p0 := by
    nlinarith [mul_le_mul_of_nonneg_left h0.2 hp1, mul_le_mul_of_nonneg_left h1.1 hp0]
  have hd : p1 / (p0 + p1) - q1 / (q0 + q1) =
      (p1 * q0 - q1 * p0) / ((p0 + p1) * (q0 + q1)) := by field_simp [hpp.ne', hqq.ne']; ring
  rw [hd, abs_div, abs_of_pos (mul_pos hpp hqq)]
  apply (div_le_iff₀ (mul_pos hpp hqq)).mpr
  have hn : |p1 * q0 - q1 * p0| ≤ ε * (p0 + p1) := abs_le.mpr ⟨by linarith, hu⟩
  have heq : (ε / η) * ((p0 + p1) * (q0 + q1)) =
      (p0 + p1) * ((ε / η) * (q0 + q1)) := by ring
  rw [heq]
  have hs : ε ≤ (ε / η) * (q0 + q1) := by
    have hm := mul_le_mul_of_nonneg_left hq (div_nonneg he hη.le)
    have hc : ε / η * η = ε := div_mul_cancel₀ _ hη.ne'
    rwa [hc] at hm
  exact hn.trans (by nlinarith [mul_le_mul_of_nonneg_left hs hpp.le])

theorem source_kernel_stability (P Q : Source) {η ε : ℝ}
    (hη : 0 < η) (he : 0 ≤ ε)
    (hp : ∀ h : Word 4, η ≤ P.law 4 h) (hq : ∀ h : Word 4, η ≤ Q.law 4 h)
    (hm : ∀ j, |powerSummaries P j - powerSummaries Q j| ≤ ε) (h : Word 4) :
    |nextKernel 4 (P.law 5) h true - nextKernel 4 (Q.law 5) h true| ≤
      (385 / 16 : ℝ) * ε / η := by
  apply binary_ratio_stability (P.nonneg _ _) (P.nonneg _ _)
    hη (by rw [prefix_mass]; exact hp h) (by rw [prefix_mass]; exact hq h)
  · exact source_block_stability P Q he hm _
  · exact source_block_stability P Q he hm _

#print axioms reconstruction_integer_bound
#print axioms reconstruction_row_integer
#print axioms reconstruction_row_bound
#print axioms reconstruction_affine
#print axioms reconstruction_stability
#print axioms source_block_stability
#print axioms binary_ratio_stability
#print axioms source_kernel_stability

end TimeFrequency.ContextTree.DepthEvidence
