import Mathlib.Analysis.Fourier.ZMod
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-! # Two harmonics in a circular block of 24 samples -/

namespace TimeFrequency.Haar

open Finset Complex
open scoped ComplexConjugate
noncomputable section

abbrev Site := ZMod 24

inductive Pitch | three | four deriving Fintype, DecidableEq

def bin : Pitch → Site | .three => 3 | .four => 4

def carrier (j i : Site) : ℂ := ZMod.stdAddChar (j*i)
def tone (a : ℂ) (j i : Site) : ℝ := (a * carrier j i).re

def phasor (a φ : ℝ) : ℂ := a * exp (φ * I)

theorem phasor_power (a φ : ℝ) : normSq (phasor a φ) = a^2 := by
  simp [phasor, normSq_eq_norm_sq]

theorem tone_cos (a φ : ℝ) (j i : ℕ) :
    tone (phasor a φ) j i = a * Real.cos (φ + 2*Real.pi*j*i/24) := by
  have hc : carrier j i = exp (((2*Real.pi*j*i/24 : ℝ) : ℂ)*I) := by
    unfold carrier
    rw [← Nat.cast_mul]
    have h := ZMod.stdAddChar_coe (N := 24) (j*i : ℤ)
    push_cast at h ⊢
    convert h using 2
    ring
  rw [tone, phasor, hc, mul_assoc, ← Complex.exp_add, ← add_mul]
  rw [← Complex.ofReal_add]
  simp only [Complex.mul_re, ofReal_re, ofReal_im, zero_mul, sub_zero,
    Complex.exp_re, Complex.mul_im, Complex.I_im, Complex.I_re, mul_zero, add_zero,
    mul_one, Real.exp_zero, one_mul]

theorem carrier_sum (j : Site) :
    ∑ i, carrier j i = if j = 0 then 24 else 0 := by
  split_ifs with h
  · simp [carrier, h, ZMod.card]
  · exact AddChar.sum_eq_zero_of_ne_one (ZMod.isPrimitive_stdAddChar 24 h)

private theorem sum_tone (a : ℂ) (j : Site) :
    ∑ i, tone a j i = if j = 0 then 24*a.re else 0 := by
  simp only [tone, ← re_sum, ← mul_sum, carrier_sum]
  split_ifs <;> simp
  ring

private theorem product_tone (a b : ℂ) (j k i : Site) :
    2*tone a j i*tone b k i =
      tone (a*b) (j+k) i + tone (a*conj b) (j-k) i := by
  have hc : carrier (j-k) i = carrier j i * conj (carrier k i) := by
    simp only [carrier, sub_eq_add_neg, add_mul, neg_mul, AddChar.map_add_eq_mul,
      AddChar.map_neg_eq_conj]
  have hd : carrier (j+k) i = carrier j i * carrier k i := by
    simp [carrier, add_mul, AddChar.map_add_eq_mul]
  simp only [tone, hc, hd, mul_re, mul_im, conj_re, conj_im]
  ring

theorem tone_inner (a b : ℂ) (j k : Site) :
    2 * (∑ i, tone a j i * tone b k i) =
      (if j+k = 0 then 24*(a*b).re else 0) +
      (if j=k then 24*(a*conj b).re else 0) := by
  simp only [mul_sum, ← mul_assoc, product_tone, sum_add_distrib, sum_tone, sub_eq_zero]

def meanSquare (x : Site → ℝ) := (∑ i, x i ^ 2) / 24

theorem two_tone_energy (p : Pitch) (a b : ℂ) :
    meanSquare (fun i ↦ tone a (bin p) i + tone b (2*bin p) i) =
      (normSq a + normSq b)/2 := by
  have h (j : Site) (z : ℂ) (hj : j+j ≠ 0) :
      2*(∑ i, tone z j i * tone z j i) = 24 * normSq z := by
    simpa [hj, mul_conj] using tone_inner z z j j
  have h1 := h (bin p) a (by cases p <;> decide)
  have h2 := h (2*bin p) b (by cases p <;> decide)
  have h3 := tone_inner a b (bin p) (2*bin p)
  rw [if_neg (by cases p <;> decide), if_neg (by cases p <;> decide)] at h3
  simp only [meanSquare, add_sq, sum_add_distrib, mul_assoc, ← mul_sum]
  simp only [sq]
  linarith

def signal (p : Pitch) (t φ ψ : ℝ) (i : Site) :=
  tone (phasor (Real.sqrt (2*(1-t))) φ) (bin p) i +
  tone (phasor (Real.sqrt (2*t)) ψ) (2*bin p) i

theorem signal_energy (p : Pitch) {t : ℝ} (ht : t ∈ Set.Icc 0 1) (φ ψ : ℝ) :
    meanSquare (signal p t φ ψ) = 1 := by
  unfold signal
  rw [two_tone_energy, phasor_power, phasor_power, Real.sq_sqrt (by linarith [ht.2]),
    Real.sq_sqrt (by linarith [ht.1])]
  ring

#print axioms phasor_power
#print axioms carrier_sum
#print axioms tone_cos
#print axioms tone_inner
#print axioms two_tone_energy
#print axioms signal_energy

end
end TimeFrequency.Haar
