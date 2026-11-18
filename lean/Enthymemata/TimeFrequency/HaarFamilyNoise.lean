import Enthymemata.TimeFrequency.HaarFamilies
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Analysis.InnerProductSpace.PiL2

/-! # Perturbing a raw block before computing normalized Haar energies -/

namespace TimeFrequency.Haar
noncomputable section

open Finset

/-- Root mean square over all 24 circular samples. -/
def blockRms (x : Block) := ‖(WithLp.toLp 2 x : EuclideanSpace ℝ Site)‖ / Real.sqrt 24

theorem blockRms_nonneg (x : Block) : 0 ≤ blockRms x := by unfold blockRms; positivity

theorem blockRms_sq (x : Block) : blockRms x ^ 2 = meanSquare x := by
  rw [blockRms, div_pow, EuclideanSpace.real_norm_sq_eq, Real.sq_sqrt (by norm_num)]
  rfl

theorem blockRms_change (x y : Block) : |blockRms x - blockRms y| ≤ blockRms (x-y) := by
  unfold blockRms
  rw [← sub_div, abs_div, abs_of_nonneg (Real.sqrt_nonneg _), WithLp.toLp_sub]
  exact div_le_div_of_nonneg_right (abs_norm_sub_norm_le _ _) (Real.sqrt_nonneg _)

private theorem filtered_sub (b : ℕ) (x y : Block) :
    filtered b (x-y) = filtered b x - filtered b y := by
  funext i
  simp only [filtered, Pi.sub_apply, sub_sub_sub_comm (x _) (y _), sum_sub_distrib, sub_div]

/-- Cauchy–Schwarz on the taps, followed by circular shift invariance of block energy. -/
theorem filtered_energy_le (b : ℕ) (hb : 0 < b) (x : Block) :
    meanSquare (filtered b x) ≤ (2*b) * meanSquare x := by
  have hpos : (0:ℝ) < 2*b := by positivity
  have tap (i : Site) : (filtered b x i)^2 ≤
      ∑ d ∈ range b, ((x (i-d))^2 + (x (i-(d+b : ℕ)))^2) := by
    rw [filtered, div_pow, Real.sq_sqrt hpos.le, div_le_iff₀ hpos]
    have hc := sum_mul_sq_le_sq_mul_sq (range b) (fun _ ↦ (1:ℝ))
      (fun d ↦ x (i-d) - x (i-(d+b : ℕ)))
    simp only [one_mul,one_pow,sum_const,card_range,nsmul_eq_mul,mul_one] at hc
    have ht : (∑ d ∈ range b, (x (i-d) - x (i-(d+b : ℕ)))^2) ≤
        2 * ∑ d ∈ range b, ((x (i-d))^2 + (x (i-(d+b : ℕ)))^2) := by
      rw [mul_sum]
      exact sum_le_sum fun d _ ↦ by nlinarith [sq_nonneg (x (i-d) + x (i-(d+b : ℕ)))]
    exact hc.trans ((mul_le_mul_of_nonneg_left ht (Nat.cast_nonneg b)).trans_eq (by ring))
  have shift (k : Site) : ∑ i, (x (i-k))^2 = ∑ i, (x i)^2 :=
    (Equiv.subRight k).sum_comp (fun i ↦ (x i)^2)
  have hs := sum_le_sum (s := univ) (fun i _ ↦ tap i)
  rw [sum_comm] at hs
  simp only [sum_add_distrib,shift,sum_const,card_range,nsmul_eq_mul] at hs
  unfold meanSquare
  linarith only [hs]

private theorem filtered_rms_le (b : ℕ) (hb : 0 < b) (x : Block) :
    blockRms (filtered b x) ≤ Real.sqrt (2*b) * blockRms x := by
  have he := filtered_energy_le b hb x
  rw [← blockRms_sq, ← blockRms_sq] at he
  have hs : (blockRms (filtered b x))^2 ≤ (Real.sqrt (2*b)*blockRms x)^2 := by
    rwa [mul_pow, Real.sq_sqrt (by positivity)]
  exact le_of_abs_le ((sq_le_sq.mp hs).trans_eq (abs_of_nonneg
    (mul_nonneg (Real.sqrt_nonneg _) (blockRms_nonneg _))))

/-- A Lipschitz bound for one energy divided by raw block power. -/
theorem normalized_band_change (b : ℕ) (hb : 0 < b) (x y : Block)
    (hx : 0 < blockRms x) (hy : 0 < blockRms y) :
    |meanSquare (filtered b x)/meanSquare x - meanSquare (filtered b y)/meanSquare y| ≤
      8*b*(blockRms (x-y)/blockRms x) := by
  let k : ℝ := Real.sqrt (2*b)
  let e := blockRms (x-y)
  let A := blockRms (filtered b x)/blockRms x
  let B := blockRms (filtered b y)/blockRms y
  have hk : 0 ≤ k := Real.sqrt_nonneg _
  have he : 0 ≤ e := blockRms_nonneg _
  have hA : 0 ≤ A := div_nonneg (blockRms_nonneg _) hx.le
  have hB : 0 ≤ B := div_nonneg (blockRms_nonneg _) hy.le
  have hAk : A ≤ k := (div_le_iff₀ hx).mpr (filtered_rms_le b hb x)
  have hBk : B ≤ k := (div_le_iff₀ hy).mpr (filtered_rms_le b hb y)
  have hband : |blockRms (filtered b x)-blockRms (filtered b y)| ≤ k*e := by
    exact (blockRms_change _ _).trans ((filtered_sub b x y) ▸ filtered_rms_le b hb (x-y))
  have hraw : |blockRms y-blockRms x| ≤ e := by rw [abs_sub_comm]; exact blockRms_change x y
  have hdiff : |A-B| ≤ 2*k*(e/blockRms x) := by
    have hidentity : (A-B)*blockRms x = (blockRms (filtered b x)-blockRms (filtered b y)) +
        B*(blockRms y-blockRms x) := by
      dsimp [A,B]
      field_simp
      ring
    have h := abs_add_le (blockRms (filtered b x) - blockRms (filtered b y))
      (B * (blockRms y - blockRms x))
    rw [← hidentity, abs_mul, abs_of_pos hx, abs_mul, abs_of_nonneg hB] at h
    have hh := mul_le_mul hBk hraw (abs_nonneg _) hk
    rw [← mul_div_assoc]
    apply (le_div_iff₀ hx).mpr
    linarith
  have hsq : |A^2-B^2| ≤ (2*k*(e/blockRms x))*(2*k) := by
    rw [sq_sub_sq,abs_mul,abs_of_nonneg (add_nonneg hA hB)]
    rw [mul_comm (2*k*(e/blockRms x))]
    exact mul_le_mul (by linarith : A+B ≤ 2*k) hdiff (abs_nonneg _) (by positivity)
  have hks : k^2 = 2*b := Real.sq_sqrt (by positivity)
  calc
    _ = |A^2-B^2| := by simp only [A,B,div_pow,blockRms_sq]
    _ ≤ _ := hsq
    _ = 8*b*(blockRms (x-y)/blockRms x) := by
      calc
        _ = 4*k^2*(e/blockRms x) := by ring
        _ = _ := by rw [hks]; ring

/-- Euclidean error in the three normalized observations, in units of relative raw RMS error. -/
theorem normalized_waveform_error (x y : Block) (hx : 0 < blockRms x)
    (hsmall : blockRms (x - y) < blockRms x) :
    distance3 (normalized x) (normalized y) ≤ 8*Real.sqrt 21*(blockRms (x-y)/blockRms x) := by
  have hy : 0 < blockRms y := by
    have h := (abs_sub_le_iff.mp (blockRms_change x y)).1
    linarith
  let ρ := blockRms (x-y)/blockRms x
  have hρ : 0 ≤ ρ := div_nonneg (blockRms_nonneg _) hx.le
  have band (b : ℕ) (hb : 0 < b) :
      (meanSquare (filtered b x)/meanSquare x - meanSquare (filtered b y)/meanSquare y)^2 ≤
        (8*b*ρ)^2 := by
    apply sq_le_sq.mpr
    rw [abs_of_nonneg (show (0:ℝ) ≤ 8*b*ρ by positivity)]
    exact normalized_band_change b hb x y hx hy
  have h1 := band 1 (by decide)
  have h2 := band 2 (by decide)
  have h4 := band 4 (by decide)
  norm_num only [Nat.cast_ofNat] at h1 h2 h4
  change Real.sqrt _ ≤ 8*Real.sqrt 21*ρ
  apply Real.sqrt_le_iff.mpr
  refine ⟨by positivity, ?_⟩
  let δ (b : ℕ) := meanSquare (filtered b x) / meanSquare x -
    meanSquare (filtered b y) / meanSquare y
  change δ 1 * δ 1 + δ 2 * δ 2 + δ 4 * δ 4 ≤ _
  change (δ 1)^2 ≤ _ at h1
  change (δ 2)^2 ≤ _ at h2
  change (δ 4)^2 ≤ _ at h4
  have hs := Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 21)
  nlinarith only [h1,h2,h4,show (8*Real.sqrt 21*ρ)^2 = 64*21*ρ^2 by rw [mul_pow,mul_pow,hs]; ring]

/-- Raw waveform tolerance for the explicit collision intervals with margin 1/16. -/
theorem collision_family_waveform_tolerance (x y : Block) (hx : 0 < blockRms x)
    (h : 256 * Real.sqrt 21 * blockRms (x - y) < blockRms x) :
    distance3 (normalized y) (normalized x) < 1/32 := by
  have hk : 1 ≤ Real.sqrt 21 := Real.le_sqrt_of_sq_le (by norm_num)
  have he := blockRms_nonneg (x-y)
  have hsmall : blockRms (x-y) < blockRms x := by
    nlinarith [mul_nonneg (sub_nonneg.mpr hk) he]
  rw [distance3_comm]
  apply (normalized_waveform_error x y hx hsmall).trans_lt
  rw [← mul_div_assoc, div_lt_iff₀ hx]
  linarith

#print axioms collision_family_waveform_tolerance
#print axioms blockRms_nonneg
#print axioms blockRms_sq
#print axioms blockRms_change
#print axioms filtered_energy_le
#print axioms normalized_band_change
#print axioms normalized_waveform_error

end
end TimeFrequency.Haar
