import Mathlib.Algebra.Ring.GeomSum
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.MeasureTheory.Integral.IntegralEqImproper
import Mathlib.NumberTheory.Harmonic.Defs
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.NormNum

/-!
# The peak of `n` independent unit exponentials (S12)

A single frame's largest bin is not flat for noise: the largest of `n` independent unit
exponentials has survival function `1 − (1 − e^{−x})ⁿ`, and its mean is the harmonic number:

* `integral_one_sub_pow_one_sub_exp`: `∫₀^∞ (1 − (1 − e^{−x})ⁿ) dx = H_n`.

By the geometric sum, `1 − qⁿ = (1 − q)·Σ_{k<n} qᵏ` with `q = 1 − e^{−x}`, so the integrand is the
derivative of `Σ_{k<n} q^{k+1}/(k + 1)`, which runs from `0` at `x = 0` to `H_n` at infinity.
That the maximum's mean is this integral (the layer-cake formula with independence) is not
formalized here.
-/

open Real Set Filter Topology MeasureTheory

namespace TimeFrequency

/-- The antiderivative `Σ_{k<n} (1 − e^{−x})^{k+1}/(k + 1)`. -/
noncomputable def peakPrimitive (n : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ Finset.range n, (1 - exp (-x)) ^ (k + 1) / (k + 1)

theorem hasDerivAt_peakPrimitive (n : ℕ) (x : ℝ) :
    HasDerivAt (peakPrimitive n) (1 - (1 - exp (-x)) ^ n) x := by
  have hq : HasDerivAt (fun x => 1 - exp (-x)) (exp (-x)) x := by
    have := ((hasDerivAt_neg x).exp).const_sub 1
    simpa using this
  have hterm : ∀ k ∈ Finset.range n, HasDerivAt (fun x => (1 - exp (-x)) ^ (k + 1) / (k + 1))
      ((1 - exp (-x)) ^ k * exp (-x)) x := by
    intro k _
    have := (hq.pow (k + 1)).div_const ((k : ℝ) + 1)
    refine this.congr_deriv ?_
    have hk : ((k : ℝ) + 1) ≠ 0 := by positivity
    push_cast
    field_simp
  have hsum := HasDerivAt.fun_sum hterm
  refine hsum.congr_deriv ?_
  have h := geom_sum_mul_neg (1 - exp (-x)) n
  rw [show (1 : ℝ) - (1 - exp (-x)) = exp (-x) by ring] at h
  rw [← Finset.sum_mul, h]

/-- **S12.** `∫₀^∞ (1 − (1 − e^{−x})ⁿ) dx = H_n`. -/
theorem integral_one_sub_pow_one_sub_exp (n : ℕ) :
    ∫ x in Ioi (0 : ℝ), (1 - (1 - exp (-x)) ^ n) = (harmonic n : ℝ) := by
  have hlim : Tendsto (peakPrimitive n) atTop
      (𝓝 (∑ k ∈ Finset.range n, (1 - 0) ^ (k + 1) / ((k : ℝ) + 1))) := by
    unfold peakPrimitive
    refine tendsto_finsetSum _ fun k _ => ?_
    exact ((tendsto_exp_neg_atTop_nhds_zero.const_sub 1).pow (k + 1)).div_const _
  rw [integral_Ioi_of_hasDerivAt_of_nonneg
    ((continuous_iff_continuousAt.1
      (continuous_iff_continuousAt.2 fun x => (hasDerivAt_peakPrimitive n x).continuousAt))
        0).continuousWithinAt
    (fun x _ => hasDerivAt_peakPrimitive n x) (fun x hx => ?_) hlim]
  · simp [peakPrimitive, harmonic]
  · have h0 : 0 ≤ 1 - exp (-x) := by
      have : exp (-x) ≤ 1 := exp_le_one_iff.2 (by linarith [hx.out])
      linarith
    have h1 : 1 - exp (-x) ≤ 1 := by linarith [exp_pos (-x)]
    linarith [pow_le_one₀ h0 h1 (n := n)]

/-- The witness at `n = 2`: the larger of two unit exponentials has mean `3/2`. -/
example : ∫ x in Ioi (0 : ℝ), (1 - (1 - exp (-x)) ^ 2) = 3 / 2 := by
  rw [integral_one_sub_pow_one_sub_exp]
  norm_num [harmonic]

end TimeFrequency
