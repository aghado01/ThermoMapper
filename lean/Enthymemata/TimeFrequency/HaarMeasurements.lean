import Enthymemata.TimeFrequency.HaarHarmonics

/-! # Translation-averaged Haar measurements and their frequency gains -/

namespace TimeFrequency.Haar

open Finset Complex
open scoped ComplexConjugate
noncomputable section

/-- The first half of the taps has positive sign, the second half negative sign. -/
def filtered (b : ℕ) (x : Site → ℝ) (i : Site) : ℝ :=
  (∑ d ∈ range b, (x (i-d) - x (i-(d+b : ℕ)))) / Real.sqrt (2*b)

def response (b : ℕ) (j : Site) : ℂ :=
  ∑ d ∈ range b, (carrier j (-d) - carrier j (-(d+b : ℕ)))

def gain (b : ℕ) (j : Site) := normSq (response b j) / (2*b)

private theorem delay (a : ℂ) (j i d : Site) :
    tone a j (i-d) = (a * carrier j (-d) * carrier j i).re := by
  rw [tone, carrier, sub_eq_add_neg, mul_add, AddChar.map_add_eq_mul]
  unfold carrier
  congr 1
  ring

private theorem filtered_tone (a : ℂ) (j : Site) (b : ℕ) (i : Site) :
    filtered b (tone a j) i = tone (a*response b j) j i / Real.sqrt (2*b) := by
  simp only [filtered, delay]
  simp only [tone, response, mul_sum, sum_mul, re_sum]
  congr 1
  apply sum_congr rfl
  intro d _
  simp [mul_sub, sub_mul]

theorem filtered_energy (p : Pitch) {t : ℝ} (ht : t ∈ Set.Icc 0 1)
    (φ ψ : ℝ) {b : ℕ} (hb : 0 < b) :
    meanSquare (filtered b (signal p t φ ψ)) =
      (1-t)*gain b (bin p) + t*gain b (2*bin p) := by
  let a := phasor (Real.sqrt (2*(1-t))) φ
  let c := phasor (Real.sqrt (2*t)) ψ
  let f := fun i ↦ tone (a*response b (bin p)) (bin p) i +
    tone (c*response b (2*bin p)) (2*bin p) i
  have hf (i : Site) : filtered b (signal p t φ ψ) i = f i / Real.sqrt (2*b) := by
    have hadd : filtered b (signal p t φ ψ) i =
        filtered b (tone a (bin p)) i + filtered b (tone c (2*bin p)) i := by
      simp only [filtered, signal, add_sub_add_comm, sum_add_distrib, add_div, a, c]
    rw [hadd, filtered_tone, filtered_tone, ← add_div]
  have he : meanSquare (filtered b (signal p t φ ψ)) =
      meanSquare f / (Real.sqrt (2*b))^2 := by
    simp only [meanSquare, hf, div_pow, ← sum_div]
    ring
  rw [he, Real.sq_sqrt (by positivity)]
  change meanSquare (fun i ↦ _ + _) / _ = _
  rw [two_tone_energy, normSq_mul, normSq_mul]
  simp only [a, c, phasor_power, Real.sq_sqrt (by linarith [ht.2] : 0 ≤ 2*(1-t)),
    Real.sq_sqrt (by linarith [ht.1] : 0 ≤ 2*t), gain]
  ring

private theorem carrier_pow (j : Site) (n : ℕ) :
    carrier j (-(n : Site)) = carrier j (-1) ^ n := by
  induction n with
  | zero => simp [carrier]
  | succ n ih =>
    rw [Nat.cast_succ, neg_add, carrier, mul_add, AddChar.map_add_eq_mul, pow_succ]
    exact congrArg (fun z ↦ z * carrier j (-1)) ih

private theorem one_plus_power (z : ℂ) (hz : normSq z = 1) :
    normSq (1+z) = 2*(1+z.re) := by
  simp only [normSq_apply, add_re, add_im, one_re, one_im] at *
  nlinarith

private theorem one_minus_power (z : ℂ) (hz : normSq z = 1) :
    normSq (1-z) = 2*(1-z.re) := by
  simp only [normSq_apply, sub_re, sub_im, one_re, one_im] at *
  nlinarith

theorem gain_formulas (j : Site) :
    let z := carrier j (-1)
    gain 1 j = 1-z.re ∧
    gain 2 j = (1+z.re)*(1-(z^2).re) ∧
    gain 4 j = (1+z.re)*(1+(z^2).re)*(1-(z^4).re) := by
  let z := carrier j (-1)
  have hu : normSq z = 1 := Circle.normSq_coe _
  have hu2 : normSq (z^2) = 1 := by rw [map_pow, hu, one_pow]
  have hu4 : normSq (z^4) = 1 := by rw [map_pow, hu, one_pow]
  have response_poly (b : ℕ) : response b j = ∑ d ∈ range b, (z^d-z^(d+b)) := by
    unfold response
    apply sum_congr rfl
    intro d _
    rw [carrier_pow j d, carrier_pow j (d+b)]
  have h1 : response 1 j = 1-z := by
    rw [response_poly]
    norm_num [sum_range_succ]
  have h2 : response 2 j = (1+z)*(1-z^2) := by
    rw [response_poly]
    norm_num [sum_range_succ]
    ring
  have h4 : response 4 j = (1+z)*(1+z^2)*(1-z^4) := by
    rw [response_poly]
    norm_num [sum_range_succ]
    ring
  dsimp only
  change gain 1 j = 1-z.re ∧ gain 2 j = _ ∧ gain 4 j = _
  simp only [gain, h1, h2, h4, normSq_mul, one_plus_power z hu, one_plus_power _ hu2,
    one_minus_power z hu, one_minus_power _ hu2, one_minus_power _ hu4]
  norm_num
  constructor <;> ring

#print axioms filtered_energy
#print axioms gain_formulas

end
end TimeFrequency.Haar
