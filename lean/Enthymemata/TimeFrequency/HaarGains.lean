import Enthymemata.TimeFrequency.HaarMeasurements

/-! # Exact responses at bins 3, 6, 4, and 8 -/

namespace TimeFrequency.Haar

open Complex
open scoped ComplexConjugate
noncomputable section

/-- Coordinates in order of increasing half-width, 1, 2, 4. -/
def gains (j : Site) : Fin 3 → ℝ := ![gain 1 j, gain 2 j, gain 4 j]

private theorem double_real (z : ℂ) (h : normSq z = 1) :
    (z^2).re = 2*z.re^2-1 := by
  simp only [pow_two, mul_re, normSq_apply] at *
  nlinarith

private theorem real_part (j : ℕ) :
    (carrier j (-1)).re = Real.cos (2*Real.pi*j/24) := by
  have hh := tone_cos 1 0 j 1
  simp only [tone, phasor, ofReal_one, Complex.ofReal_zero, zero_mul, Complex.exp_zero,
    one_mul, Nat.cast_one, mul_one, zero_add] at hh
  simpa [carrier, AddChar.map_neg_eq_conj] using hh

private theorem polynomial_gains (j : Site) (x : ℝ) (hx : (carrier j (-1)).re = x) :
    gains j = ![1-x, (1+x)*(1-(2*x^2-1)),
      (1+x)*(1+(2*x^2-1))*(1-(2*(2*x^2-1)^2-1))] := by
  let z := carrier j (-1)
  have hu : normSq z = 1 := Circle.normSq_coe _
  have h2 := double_real z hu
  have h4 := double_real (z^2) (by rw [map_pow, hu, one_pow])
  rw [← pow_mul] at h4
  obtain ⟨g1,g2,g4⟩ := gain_formulas j
  change z.re = x at hx
  change gain 1 j = 1 - z.re at g1
  change gain 2 j = (1+z.re)*(1-(z^2).re) at g2
  change gain 4 j = (1+z.re)*(1+(z^2).re)*(1-(z^4).re) at g4
  unfold gains
  rw [g1,g2,g4,h4,h2,hx]

theorem bin_gains :
    gains 3 = ![1-Real.sqrt 2/2,1+Real.sqrt 2/2,2+Real.sqrt 2] ∧
    gains 6 = ![1,2,0] ∧ gains 4 = ![1/2,9/4,9/8] ∧ gains 8 = ![3/2,3/4,3/8] := by
  have r := Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 2)
  have angle (j : ℕ) := real_part j
  have h3 := angle 3
  have h4 := angle 4
  have h6 := angle 6
  have h8 := angle 8
  norm_num only [Nat.cast_ofNat] at h3 h4 h6 h8
  rw [show 2*Real.pi*3/24 = Real.pi/4 by ring, Real.cos_pi_div_four] at h3
  rw [show 2*Real.pi*4/24 = Real.pi/3 by ring, Real.cos_pi_div_three] at h4
  rw [show 2*Real.pi*6/24 = Real.pi/2 by ring, Real.cos_pi_div_two] at h6
  rw [show 2*Real.pi*8/24 = Real.pi-(Real.pi/3) by ring, Real.cos_pi_sub,
    Real.cos_pi_div_three] at h8
  rw [polynomial_gains _ _ h3, polynomial_gains _ _ h6, polynomial_gains _ _ h4,
    polynomial_gains _ _ h8]
  norm_num [div_pow, r]
  all_goals ring

def features (p : Pitch) (t : ℝ) (j : Fin 3) :=
  (1-t)*gains (bin p) j + t*gains (2*bin p) j

def measured (x : Site → ℝ) : Fin 3 → ℝ :=
  ![meanSquare (filtered 1 x), meanSquare (filtered 2 x), meanSquare (filtered 4 x)]

theorem measured_signal (p : Pitch) {t : ℝ} (ht : t ∈ Set.Icc 0 1) (φ ψ : ℝ) :
    measured (signal p t φ ψ) = features p t := by
  funext j
  fin_cases j
  · exact filtered_energy p ht φ ψ (b := 1) (by decide)
  · exact filtered_energy p ht φ ψ (b := 2) (by decide)
  · exact filtered_energy p ht φ ψ (b := 4) (by decide)

#print axioms bin_gains
#print axioms measured_signal

end
end TimeFrequency.Haar
