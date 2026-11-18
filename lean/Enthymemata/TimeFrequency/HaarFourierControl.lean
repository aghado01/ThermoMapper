import Enthymemata.TimeFrequency.HaarGeometry

/-! # Full-block rectangular Fourier control for the two-scale collision -/

namespace TimeFrequency.Haar

open Finset Complex
open scoped ComplexConjugate
noncomputable section

def quadrature (x : Site → ℝ) (w : ℂ) := (∑ i, x i * tone w 3 i)/24
def binThreePower (x : Site → ℝ) := (quadrature x 1)^2 + (quadrature x I)^2

/-- Native unnormalized DFT, unit rectangular window, only bin three retained. -/
theorem fourier_control_eq (x : Site → ℝ) :
    binThreePower x = normSq (ZMod.dft (fun i ↦ (x i : ℂ)) 3)/24^2 := by
  have hr (i : Site) : (ZMod.stdAddChar (-(i*3)) * (x i : ℂ)).re = x i * tone 1 3 i := by
    rw [AddChar.map_neg_eq_conj]
    rw [mul_comm i (3:Site)]
    simp only [tone, carrier, one_mul, mul_re, conj_re, ofReal_re, ofReal_im, mul_zero, sub_zero]
    ring
  have hi (i : Site) : (ZMod.stdAddChar (-(i*3)) * (x i : ℂ)).im = x i * tone I 3 i := by
    rw [AddChar.map_neg_eq_conj]
    rw [mul_comm i (3:Site)]
    simp only [tone, carrier, mul_re, mul_im, conj_im, ofReal_re, ofReal_im, I_re, I_im,
      zero_mul, mul_zero, one_mul, zero_sub, zero_add]
    ring
  simp only [ZMod.dft_apply, smul_eq_mul, normSq_apply, re_sum, im_sum, hr, hi,
    binThreePower, quadrature]
  ring

private theorem control_projection (p : Pitch) (a b w : ℂ) :
    quadrature (fun i ↦ tone a (bin p) i + tone b (2*bin p) i) w =
      (if p = .three then (a*conj w).re else 0)/2 := by
  have h1 := tone_inner a w (bin p) 3
  have h2 := tone_inner b w (2*bin p) 3
  rw [if_neg (by cases p <;> decide)] at h1
  rw [if_neg (by cases p <;> decide), if_neg (by cases p <;> decide)] at h2
  simp only [quadrature, add_mul, sum_add_distrib]
  cases p
  · rw [if_pos (by decide)] at h1
    simp only [if_true]
    linarith
  · rw [if_neg (by decide)] at h1
    simp
    linarith

theorem bin_three_signal (p : Pitch) {t : ℝ} (ht : t ∈ Set.Icc 0 1) (φ ψ : ℝ) :
    binThreePower (signal p t φ ψ) = if p = .three then (1-t)/2 else 0 := by
  let a := phasor (Real.sqrt (2*(1-t))) φ
  have ha : normSq a = 2*(1-t) := by
    rw [show a = phasor _ _ from rfl, phasor_power, Real.sq_sqrt (by linarith [ht.2])]
  have hq (w : ℂ) : quadrature (signal p t φ ψ) w =
      (if p = .three then (a*conj w).re else 0)/2 :=
    control_projection p a (phasor (Real.sqrt (2*t)) ψ) w
  rw [binThreePower, hq, hq]
  cases p
  · simp only [if_true, map_one, mul_one, Complex.conj_I, mul_neg, Complex.neg_re,
      Complex.mul_I_re, neg_neg]
    change (a.re/2)^2 + (a.im/2)^2 = _
    simp only [normSq_apply] at ha
    nlinarith
  · norm_num [show Pitch.four ≠ Pitch.three by decide]

theorem fourier_distinguishes_collision (φ ψ α β : ℝ) :
    binThreePower (signal .four collisionU α β) < binThreePower (signal .three collisionT φ ψ) := by
  obtain ⟨h1,h2,h3,h4⟩ := collision_ranges
  rw [bin_three_signal _ ⟨by linarith,by linarith⟩, bin_three_signal _ ⟨by linarith,by linarith⟩]
  simp
  linarith

#print axioms fourier_control_eq
#print axioms bin_three_signal
#print axioms fourier_distinguishes_collision

end
end TimeFrequency.Haar
