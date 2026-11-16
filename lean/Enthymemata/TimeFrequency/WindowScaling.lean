import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Analysis.Complex.Trigonometric
import Mathlib.Analysis.Real.Sqrt
import Mathlib.MeasureTheory.Integral.Bochner.ContinuousLinearMap
import Mathlib.MeasureTheory.Integral.Bochner.SumMeasure
import Mathlib.MeasureTheory.Measure.Count
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# Window scaling (S2)

The power a window of `w` samples reads from a tone, an impulse and white noise. Sample `j` sits at
time `tⱼ`, and the windowed transform at frequency `ω` is `X(ω) = Σⱼ gⱼ·xⱼ·e^{−iωtⱼ}` (`windowed`).

* `norm_windowed_tone_sq`: a tone `xⱼ = A·e^{iωtⱼ}` reads `|A|²·|Σⱼ gⱼ|²` at its own frequency.
* `norm_windowed_impulse_sq`: an impulse of amplitude `a` at sample `j₀` reads `|a|²·|g_{j₀}|²`
  at every frequency.
* `noisePower_white`: noise with second moments `E[εⱼ·conj εₖ] = σ²·δⱼₖ` has expected power
  `σ²·Σⱼ|gⱼ|²` at every frequency. `noisePower` is the expected power as a function of the
  second-moment matrix alone; `integral_norm_windowed_sq` says it is the expectation for any
  random samples with those second moments.
* Under the energy normalization `Σⱼ|gⱼ|² = 1`, the rectangular window `gⱼ = 1/√w` (`rect`) reads
  `|A|²·w`, `|a|²/w` and `σ²`: slopes `+1`, `−1` and `0` in log power against log window length.
* `norm_windowed_tone_sq_le`: no window of `w` samples with unit energy reads more than `|A|²·w`
  from a tone, the rectangular window's power.
-/

open Complex MeasureTheory

namespace TimeFrequency

/-- The phase factor `e^{−iωτ}`. -/
noncomputable def phase (ω τ : ℝ) : ℂ :=
  exp (↑(-(ω * τ)) * I)

theorem norm_phase (ω τ : ℝ) : ‖phase ω τ‖ = 1 :=
  norm_exp_ofReal_mul_I _

variable {ι : Type*} [Fintype ι]

/-- The windowed transform at frequency `ω`, `Σⱼ gⱼ·xⱼ·e^{−iωtⱼ}`, sample `j` at time `tⱼ`. -/
noncomputable def windowed (g x : ι → ℂ) (t : ι → ℝ) (ω : ℝ) : ℂ :=
  ∑ j, g j * x j * phase ω (t j)

/-- A tone of amplitude `A` and frequency `ω`, sampled at the times `t`. -/
noncomputable def tone (A : ℂ) (ω : ℝ) (t : ι → ℝ) : ι → ℂ :=
  fun j => A * exp (↑(ω * t j) * I)

theorem windowed_tone (g : ι → ℂ) (t : ι → ℝ) (A : ℂ) (ω : ℝ) :
    windowed g (tone A ω t) t ω = A * ∑ j, g j := by
  simp only [windowed, tone, Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  have h : exp (↑(ω * t j) * I) * phase ω (t j) = 1 := by
    rw [phase, ← exp_add,
      show (↑(ω * t j) * I + ↑(-(ω * t j)) * I : ℂ) = 0 by push_cast; ring, exp_zero]
  calc g j * (A * exp (↑(ω * t j) * I)) * phase ω (t j)
      = A * g j * (exp (↑(ω * t j) * I) * phase ω (t j)) := by ring
    _ = A * g j := by rw [h, mul_one]

/-- **S2, the tone.** At its own frequency a tone reads `|A|²·|Σⱼ gⱼ|²`. -/
theorem norm_windowed_tone_sq (g : ι → ℂ) (t : ι → ℝ) (A : ℂ) (ω : ℝ) :
    ‖windowed g (tone A ω t) t ω‖ ^ 2 = ‖A‖ ^ 2 * ‖∑ j, g j‖ ^ 2 := by
  rw [windowed_tone, norm_mul, mul_pow]

theorem windowed_impulse [DecidableEq ι] (g : ι → ℂ) (t : ι → ℝ) (a : ℂ) (j₀ : ι) (ω : ℝ) :
    windowed g (Pi.single j₀ a) t ω = g j₀ * a * phase ω (t j₀) := by
  unfold windowed
  rw [Finset.sum_eq_single j₀ (fun j _ hj => by simp [Pi.single_eq_of_ne hj])
    (fun h => absurd (Finset.mem_univ _) h), Pi.single_eq_same]

/-- **S2, the impulse.** An impulse of amplitude `a` at sample `j₀` reads `|a|²·|g_{j₀}|²` at
every frequency. -/
theorem norm_windowed_impulse_sq [DecidableEq ι] (g : ι → ℂ) (t : ι → ℝ) (a : ℂ) (j₀ : ι)
    (ω : ℝ) : ‖windowed g (Pi.single j₀ a) t ω‖ ^ 2 = ‖a‖ ^ 2 * ‖g j₀‖ ^ 2 := by
  rw [windowed_impulse, norm_mul, norm_mul, norm_phase]
  ring

/-- The expected power at frequency `ω` of samples whose second moments are
`C j k = E[εⱼ·conj εₖ]`: `Σⱼₖ uⱼ·C j k·conj uₖ` with `uⱼ = gⱼ·e^{−iωtⱼ}`. -/
noncomputable def noisePower (g : ι → ℂ) (t : ι → ℝ) (C : ι → ι → ℂ) (ω : ℝ) : ℂ :=
  ∑ j, ∑ k, g j * phase ω (t j) * C j k * (starRingEnd ℂ) (g k * phase ω (t k))

/-- **S2, the noise, on its second moments.** White second moments `σ²·δⱼₖ` give expected power
`σ²·Σⱼ|gⱼ|²` at every frequency. -/
theorem noisePower_white [DecidableEq ι] (g : ι → ℂ) (t : ι → ℝ) (s : ℝ) (ω : ℝ) :
    noisePower g t (fun j k => if j = k then (s : ℂ) else 0) ω =
      ((s * ∑ j, ‖g j‖ ^ 2 : ℝ) : ℂ) := by
  unfold noisePower
  simp only [mul_ite, ite_mul, mul_zero, zero_mul, Finset.sum_ite_eq, Finset.mem_univ, if_true]
  push_cast
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  calc g j * phase ω (t j) * (s : ℂ) * (starRingEnd ℂ) (g j * phase ω (t j))
      = (s : ℂ) * (g j * phase ω (t j) * (starRingEnd ℂ) (g j * phase ω (t j))) := by ring
    _ = (s : ℂ) * (‖g j‖ : ℂ) ^ 2 := by
      rw [mul_conj', norm_mul, norm_phase, mul_one]

/-- **The expectation.** For random samples `ε` whose products `εⱼ·conj εₖ` are integrable with
means `C j k`, the expected power of the windowed transform is `noisePower g t C ω`. -/
theorem integral_norm_windowed_sq {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) (g : ι → ℂ)
    (t : ι → ℝ) (ε : ι → Ω → ℂ) (C : ι → ι → ℂ)
    (hint : ∀ j k, Integrable (fun w => ε j w * (starRingEnd ℂ) (ε k w)) μ)
    (hmom : ∀ j k, ∫ w, ε j w * (starRingEnd ℂ) (ε k w) ∂μ = C j k) (ω : ℝ) :
    ((∫ w, ‖windowed g (fun j => ε j w) t ω‖ ^ 2 ∂μ : ℝ) : ℂ) = noisePower g t C ω := by
  rw [← integral_complex_ofReal]
  have h : ∀ w, ((‖windowed g (fun j => ε j w) t ω‖ ^ 2 : ℝ) : ℂ) =
      ∑ j, ∑ k, (g j * phase ω (t j) * (starRingEnd ℂ) (g k * phase ω (t k))) *
        (ε j w * (starRingEnd ℂ) (ε k w)) := by
    intro w
    push_cast
    rw [← mul_conj', windowed, map_sum, Finset.sum_mul_sum]
    refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
    simp only [map_mul]
    ring
  simp_rw [h]
  rw [integral_finsetSum _ fun j _ => integrable_finsetSum _ fun k _ => (hint j k).const_mul _]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [integral_finsetSum _ fun k _ => (hint j k).const_mul _]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [integral_const_mul, hmom]
  ring

/-- **S2, the noise.** Random samples with white second moments `σ²·δⱼₖ`, through a window of
unit energy, have expected power `σ²` at every frequency. -/
theorem integral_norm_windowed_sq_white [DecidableEq ι] {Ω : Type*} [MeasurableSpace Ω]
    (μ : Measure Ω) (g : ι → ℂ) (hg : ∑ j, ‖g j‖ ^ 2 = 1) (t : ι → ℝ) (ε : ι → Ω → ℂ) (s : ℝ)
    (hint : ∀ j k, Integrable (fun w => ε j w * (starRingEnd ℂ) (ε k w)) μ)
    (hmom : ∀ j k, ∫ w, ε j w * (starRingEnd ℂ) (ε k w) ∂μ = if j = k then (s : ℂ) else 0)
    (ω : ℝ) : ∫ w, ‖windowed g (fun j => ε j w) t ω‖ ^ 2 ∂μ = s := by
  have h := integral_norm_windowed_sq μ g t ε _ hint hmom ω
  rw [noisePower_white, hg, mul_one] at h
  exact_mod_cast h

/-- Cauchy–Schwarz: `|Σⱼ gⱼ|² ≤ w·Σⱼ|gⱼ|²` over `w` samples. -/
theorem norm_sum_sq_le (g : ι → ℂ) : ‖∑ j, g j‖ ^ 2 ≤ Fintype.card ι * ∑ j, ‖g j‖ ^ 2 := by
  have h := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (fun j => ‖g j‖) (fun _ => (1 : ℝ))
  simp only [mul_one, one_pow, Finset.sum_const, Finset.card_univ, nsmul_eq_mul] at h
  calc ‖∑ j, g j‖ ^ 2 ≤ (∑ j, ‖g j‖) ^ 2 := by
        gcongr
        exact norm_sum_le _ _
    _ ≤ (∑ j, ‖g j‖ ^ 2) * Fintype.card ι := h
    _ = _ := mul_comm _ _

/-- No window of `w` samples with unit energy reads more than `|A|²·w` from a tone. -/
theorem norm_windowed_tone_sq_le (g : ι → ℂ) (hg : ∑ j, ‖g j‖ ^ 2 = 1) (t : ι → ℝ) (A : ℂ)
    (ω : ℝ) : ‖windowed g (tone A ω t) t ω‖ ^ 2 ≤ ‖A‖ ^ 2 * Fintype.card ι := by
  rw [norm_windowed_tone_sq]
  have := norm_sum_sq_le g
  rw [hg, mul_one] at this
  exact mul_le_mul_of_nonneg_left this (sq_nonneg _)

/-! ### The rectangular window -/

/-- The rectangular window of `w` samples with unit energy, `gⱼ = 1/√w`. -/
noncomputable def rect (w : ℕ) : Fin w → ℂ :=
  fun _ => (((Real.sqrt w)⁻¹ : ℝ) : ℂ)

theorem norm_rect_sq {w : ℕ} (j : Fin w) : ‖rect w j‖ ^ 2 = (w : ℝ)⁻¹ := by
  rw [rect, norm_real, Real.norm_eq_abs, sq_abs, inv_pow, Real.sq_sqrt (Nat.cast_nonneg _)]

theorem rect_energy {w : ℕ} (hw : 0 < w) : ∑ j, ‖rect w j‖ ^ 2 = 1 := by
  simp only [norm_rect_sq, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  exact mul_inv_cancel₀ (Nat.cast_ne_zero.2 hw.ne')

theorem norm_sum_rect_sq {w : ℕ} (hw : 0 < w) : ‖∑ j, rect w j‖ ^ 2 = w := by
  have hw0 : (w : ℝ) ≠ 0 := Nat.cast_ne_zero.2 hw.ne'
  simp only [rect, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  rw [norm_mul, mul_pow, norm_natCast, norm_real, Real.norm_eq_abs, sq_abs, inv_pow,
    Real.sq_sqrt (Nat.cast_nonneg _), sq, mul_assoc, mul_inv_cancel₀ hw0, mul_one]

/-- **S2, rectangular corollary: the tone** reads `|A|²·w`. -/
theorem rect_tone {w : ℕ} (hw : 0 < w) (t : Fin w → ℝ) (A : ℂ) (ω : ℝ) :
    ‖windowed (rect w) (tone A ω t) t ω‖ ^ 2 = ‖A‖ ^ 2 * w := by
  rw [norm_windowed_tone_sq, norm_sum_rect_sq hw]

/-- **S2, rectangular corollary: the impulse** reads `|a|²/w` at every frequency. -/
theorem rect_impulse {w : ℕ} (t : Fin w → ℝ) (a : ℂ) (j₀ : Fin w) (ω : ℝ) :
    ‖windowed (rect w) (Pi.single j₀ a) t ω‖ ^ 2 = ‖a‖ ^ 2 / w := by
  rw [norm_windowed_impulse_sq, norm_rect_sq, div_eq_mul_inv]

/-- **S2, rectangular corollary: the noise** reads `σ²` at every frequency. -/
theorem rect_noise {w : ℕ} (hw : 0 < w) (t : Fin w → ℝ) (s : ℝ) (ω : ℝ) :
    noisePower (rect w) t (fun j k => if j = k then (s : ℂ) else 0) ω = s := by
  rw [noisePower_white, rect_energy hw, mul_one]

/-! ### Witnesses -/

/-- A window of four samples at times `0, 1, 2, 3` meets every hypothesis: the tone reads
`|A|²·4`. -/
example : ‖windowed (rect 4) (tone 2 1 fun j => (j : ℝ)) (fun j => (j : ℝ)) 1‖ ^ 2 = 16 := by
  rw [rect_tone (by norm_num)]
  norm_num

/-- The impulse at sample `1` of a four-sample rectangular window reads `|a|²/4`. -/
example : ‖windowed (rect 4) (Pi.single 1 2) (fun j => (j : ℝ)) 3‖ ^ 2 = 1 := by
  rw [rect_impulse]
  norm_num

/-- The white-noise expectation's hypotheses are met on a concrete space: two outcomes under
the counting measure, sample `j` equal to `1` on outcome `j` and `0` elsewhere, so
`E[εⱼ·conj εₖ] = δⱼₖ`. -/
example : ∫ w, ‖windowed (rect 2) (fun j => if j = w then (1 : ℂ) else 0) (fun j => (j : ℝ)) 5‖ ^ 2
    ∂(Measure.count : Measure (Fin 2)) = 1 :=
  integral_norm_windowed_sq_white _ _ (rect_energy (by norm_num)) _
    (fun j w => if j = w then (1 : ℂ) else 0) 1 (fun _ _ => Integrable.of_finite)
    (fun j k => by
      rw [integral_fintype Integrable.of_finite]
      fin_cases j <;> fin_cases k <;> simp) 5

end TimeFrequency
