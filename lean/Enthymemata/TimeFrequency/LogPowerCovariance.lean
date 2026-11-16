import Mathlib.Algebra.BigOperators.Field
import Mathlib.Analysis.MellinTransform
import Mathlib.Analysis.SpecialFunctions.Gamma.Deriv
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Data.Rat.BigOperators
import Mathlib.Analysis.Real.Sqrt
import Mathlib.Data.Real.Basic
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.MeasureTheory.Integral.Bochner.ContinuousLinearMap
import Mathlib.MeasureTheory.Integral.Bochner.SumMeasure
import Mathlib.MeasureTheory.Measure.Count
import Mathlib.NumberTheory.Harmonic.Defs
import Mathlib.NumberTheory.Harmonic.GammaDeriv
import Mathlib.NumberTheory.ZetaValues
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# The covariance of log powers (S10) and powers of a non-circular pair (S11)

The covariance of `ln|a|²` and `ln|b|²` for jointly circular complex Gaussians with squared
coherence `κ` is the dilogarithm `Li₂(κ) = Σ_{n≥1} κⁿ/n²`.

* `alternating_choose_harmonic` (S10a): `Σ_{k=0}^{n} (−1)ᵏ·C(n, k)·H_k = −1/n` for `n ≥ 1`.
* `integral_pow_mul_log_mul_exp`: `∫₀^∞ xᵏ·ln x·e^{−x} dx = k!·(H_k − γ)`. This is `Γ′(k + 1)`,
  read from Mathlib's `Complex.hasDerivAt_Gamma_nat` and the derivative of the Gamma integral.
  `integrableOn_pow_mul_log_mul_exp` supplies integrability from the Mellin derivative.
* `integral_log_mul_laguerre` (S10b): `∫₀^∞ ln x·L_n(x)·e^{−x} dx = −1/n` for `n ≥ 1`. The `γ`
  terms cancel because `Σ_k (−1)ᵏ·C(n, k) = 0`.
* `logPower_covariance` (S10c). The premise is `LaguerreProductPremise`: the bivariate expansion
  `e^{−x−y}·Σ κⁿ·L_n(x)·L_n(y)` of the pair's joint density, integrated term by term against
  `ln x·ln y`. Then `Cov(ln|a|², ln|b|²) = Σ_{n≥1} κⁿ/n²`. `dilog_one`: at `κ = 1` the series is
  `π²/6`.
* `cov_normSq_pairing` (S11): with the complex pairing rule as premise,
  `Cov(|a|², |b|²) = |E[a·conj b]|² + |E[a·b]|²`. The second term is the pseudo-Gram entry.
-/

open Finset Real Set Filter Topology MeasureTheory Asymptotics ComplexConjugate

namespace TimeFrequency

/-- `Σ_{k=0}^{m} (−1)ᵏ·C(m, k)/(k + 1) = 1/(m + 1)`, the integral of `(1 − x)ᵐ` over `[0, 1]`. -/
theorem sum_alternating_choose_div_succ (m : ℕ) :
    ∑ k ∈ range (m + 1), (-1 : ℚ) ^ k * m.choose k / (k + 1) = 1 / (m + 1) := by
  have hterm : ∀ k, (-1 : ℚ) ^ k * m.choose k / (k + 1) =
      (-1) ^ k * (m + 1).choose (k + 1) / (m + 1) := by
    intro k
    have h' : ((m : ℚ) + 1) * m.choose k = (m + 1).choose (k + 1) * (k + 1) := by
      exact_mod_cast Nat.add_one_mul_choose_eq m k
    have hk1 : (k : ℚ) + 1 ≠ 0 := by positivity
    have hm1 : (m : ℚ) + 1 ≠ 0 := by positivity
    rw [div_eq_div_iff hk1 hm1]
    linear_combination (-1 : ℚ) ^ k * h'
  simp_rw [hterm]
  rw [← Finset.sum_div]
  have halt : (∑ i ∈ range (m + 1 + 1), (-1 : ℚ) ^ i * ((m + 1).choose i : ℚ)) = 0 := by
    exact_mod_cast Int.alternating_sum_range_choose_of_ne (n := m + 1) (Nat.succ_ne_zero m)
  rw [Finset.sum_range_succ'] at halt
  have hneg : ∑ k ∈ range (m + 1), (-1 : ℚ) ^ (k + 1) * ((m + 1).choose (k + 1) : ℚ) =
      -∑ k ∈ range (m + 1), (-1 : ℚ) ^ k * ((m + 1).choose (k + 1) : ℚ) := by
    rw [← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    ring
  rw [hneg] at halt
  simp only [pow_zero, Nat.choose_zero_right, Nat.cast_one, mul_one] at halt
  have hs : ∑ k ∈ range (m + 1), (-1 : ℚ) ^ k * ((m + 1).choose (k + 1) : ℚ) = 1 := by linarith
  rw [hs]

/-- **S10a.** For `n ≥ 1`, `Σ_{k=0}^{n} (−1)ᵏ·C(n, k)·H_k = −1/n`, written at `n + 1`. -/
theorem alternating_choose_harmonic (n : ℕ) :
    ∑ k ∈ range (n + 1 + 1), (-1 : ℚ) ^ k * (n + 1).choose k * harmonic k = -1 / (n + 1) := by
  -- the same sum one row down, read with its index shifted
  have h3 : ∑ j ∈ range (n + 1), (-1 : ℚ) ^ (j + 1) * n.choose (j + 1) * harmonic (j + 1) =
      ∑ k ∈ range (n + 1), (-1 : ℚ) ^ k * n.choose k * harmonic k := by
    have e1 := Finset.sum_range_succ' (fun k => (-1 : ℚ) ^ k * n.choose k * harmonic k) (n + 1)
    have e2 := Finset.sum_range_succ (fun k => (-1 : ℚ) ^ k * n.choose k * harmonic k) (n + 1)
    simp only [harmonic_zero, mul_zero, add_zero, Nat.choose_succ_self, Nat.cast_zero,
      zero_mul] at e1 e2
    rw [← e1, e2]
  -- Pascal's rule and `H_{j+1} = H_j + 1/(j + 1)`, term by term
  have key : ∀ j, (-1 : ℚ) ^ (j + 1) * (n + 1).choose (j + 1) * harmonic (j + 1) =
      -((-1 : ℚ) ^ j * n.choose j * harmonic j) - (-1 : ℚ) ^ j * n.choose j / (j + 1) +
        (-1 : ℚ) ^ (j + 1) * n.choose (j + 1) * harmonic (j + 1) := by
    intro j
    have hH : harmonic (j + 1) = harmonic j + 1 / ((j : ℚ) + 1) := by
      rw [harmonic_succ]
      push_cast
      ring
    rw [Nat.choose_succ_succ', Nat.cast_add]
    linear_combination (-((-1 : ℚ) ^ j * n.choose j)) * hH
  rw [Finset.sum_range_succ', harmonic_zero, mul_zero, add_zero,
    Finset.sum_congr rfl fun j _ => key j, Finset.sum_add_distrib, Finset.sum_sub_distrib,
    Finset.sum_neg_distrib, h3, sum_alternating_choose_div_succ]
  ring

/-- S10a over `ℝ`. -/
theorem alternating_choose_harmonic_real (n : ℕ) :
    ∑ k ∈ range (n + 1 + 1), (-1 : ℝ) ^ k * (n + 1).choose k * (harmonic k : ℝ) =
      -1 / (n + 1) := by
  have h := congrArg (fun q : ℚ => (q : ℝ)) (alternating_choose_harmonic n)
  simp only [Rat.cast_sum, Rat.cast_mul, Rat.cast_pow, Rat.cast_neg, Rat.cast_one,
    Rat.cast_natCast, Rat.cast_div, Rat.cast_add] at h
  exact h

/-! ### S10b, the Laguerre coefficient of the logarithm -/

/-- `∫₀^∞ xᵏ·ln x·e^{−x} dx = k!·(H_k − γ)`, which is `Γ′(k + 1)`. -/
theorem integral_pow_mul_log_mul_exp (k : ℕ) :
    ∫ x in Ioi (0 : ℝ), x ^ k * (log x * exp (-x)) =
      k.factorial * ((harmonic k : ℝ) - eulerMascheroniConstant) := by
  have hs : 0 < ((k : ℂ) + 1).re := by
    simp only [Complex.add_re, Complex.natCast_re, Complex.one_re]
    positivity
  have h1 := Complex.hasDerivAt_GammaIntegral hs
  have hev : Complex.Gamma =ᶠ[𝓝 ((k : ℂ) + 1)] Complex.GammaIntegral := by
    filter_upwards [(Complex.continuous_re.isOpen_preimage _ isOpen_Ioi).mem_nhds hs] with s hs'
    exact Complex.Gamma_eq_integral hs'
  have heq := (h1.congr_of_eventuallyEq hev).unique (Complex.hasDerivAt_Gamma_nat k)
  have hint : (∫ t in Ioi (0 : ℝ), (t : ℂ) ^ ((k : ℂ) + 1 - 1) *
      ((Real.log t : ℂ) * (Real.exp (-t) : ℂ))) =
      ((∫ x in Ioi (0 : ℝ), x ^ k * (log x * exp (-x)) : ℝ) : ℂ) := by
    rw [← integral_complex_ofReal]
    refine setIntegral_congr_fun measurableSet_Ioi fun t _ => ?_
    simp only [add_sub_cancel_right, Complex.cpow_natCast]
    push_cast
    ring
  rw [hint] at heq
  have h2 : ((∫ x in Ioi (0 : ℝ), x ^ k * (log x * exp (-x)) : ℝ) : ℂ) =
      ((k.factorial * ((harmonic k : ℝ) - eulerMascheroniConstant) : ℝ) : ℂ) := by
    rw [heq]
    push_cast
    ring
  exact_mod_cast h2

/-- `xᵏ·ln x·e^{−x}` is integrable on `(0, ∞)`, from the Mellin derivative of `e^{−x}`. -/
theorem integrableOn_pow_mul_log_mul_exp (k : ℕ) :
    IntegrableOn (fun x : ℝ => x ^ k * (log x * exp (-x))) (Ioi 0) := by
  have hs : 0 < ((k : ℂ) + 1).re := by
    simp only [Complex.add_re, Complex.natCast_re, Complex.one_re]
    positivity
  have hmc := (mellin_hasDerivAt_of_isBigO_rpow (E := ℂ)
    (f := fun t : ℝ => ((Real.exp (-t) : ℝ) : ℂ)) (a := ((k : ℂ) + 1).re + 1) (b := 0)
    (s := (k : ℂ) + 1) ?_ ?_ (lt_add_one _) ?_ hs).1
  · have hc : IntegrableOn (fun t : ℝ => (t : ℂ) ^ ((k : ℂ) + 1 - 1) •
        (Real.log t • ((Real.exp (-t) : ℝ) : ℂ))) (Ioi 0) := hmc
    refine IntegrableOn.congr_fun (Integrable.re hc) (fun t _ => ?_) measurableSet_Ioi
    have e : (t : ℂ) ^ ((k : ℂ) + 1 - 1) • (Real.log t • ((Real.exp (-t) : ℝ) : ℂ)) =
        ((t ^ k * (log t * exp (-t)) : ℝ) : ℂ) := by
      simp only [add_sub_cancel_right, Complex.cpow_natCast, smul_eq_mul, Complex.real_smul]
      push_cast
      ring
    simp only [e]
    exact Complex.ofReal_re _
  · exact (Continuous.continuousOn (Complex.continuous_ofReal.comp
      (Real.continuous_exp.comp continuous_neg))).locallyIntegrableOn measurableSet_Ioi
  · rw [← isBigO_norm_left]
    simp_rw [Complex.norm_real, isBigO_norm_left]
    simpa only [neg_one_mul] using (isLittleO_exp_neg_mul_rpow_atTop zero_lt_one _).isBigO
  · simp_rw [neg_zero, Real.rpow_zero]
    refine isBigO_const_of_tendsto (?_ : Tendsto _ _ (𝓝 (1 : ℂ))) one_ne_zero
    rw [(by simp : (1 : ℂ) = Real.exp (-0))]
    exact (Complex.continuous_ofReal.comp
      (Real.continuous_exp.comp continuous_neg)).continuousWithinAt

/-- The Laguerre polynomial `L_n(x) = Σ_k C(n, k)·(−x)ᵏ/k!`. -/
noncomputable def laguerre (n : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ range (n + 1), (n.choose k : ℝ) * (-x) ^ k / k.factorial

/-- **S10b.** For `n ≥ 1`, `∫₀^∞ ln x·L_n(x)·e^{−x} dx = −1/n`, written at `n + 1`. -/
theorem integral_log_mul_laguerre (n : ℕ) :
    ∫ x in Ioi (0 : ℝ), log x * laguerre (n + 1) x * exp (-x) = -1 / (n + 1) := by
  have hpt : ∀ x, log x * laguerre (n + 1) x * exp (-x) =
      ∑ k ∈ range (n + 1 + 1), ((n + 1).choose k * (-1) ^ k / k.factorial : ℝ) *
        (x ^ k * (log x * exp (-x))) := by
    intro x
    unfold laguerre
    rw [Finset.mul_sum, Finset.sum_mul]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [neg_pow]
    ring
  simp_rw [hpt]
  rw [integral_finsetSum _ fun k _ => Integrable.const_mul (integrableOn_pow_mul_log_mul_exp k) _]
  simp_rw [integral_const_mul, integral_pow_mul_log_mul_exp]
  have h1 := alternating_choose_harmonic_real n
  have h0 : ∑ k ∈ range (n + 1 + 1), (-1 : ℝ) ^ k * ((n + 1).choose k : ℝ) = 0 := by
    exact_mod_cast Int.alternating_sum_range_choose_of_ne (n := n + 1) (Nat.succ_ne_zero n)
  have hsplit : ∑ k ∈ range (n + 1 + 1), ((n + 1).choose k * (-1) ^ k / k.factorial : ℝ) *
      (k.factorial * ((harmonic k : ℝ) - eulerMascheroniConstant)) =
      ∑ k ∈ range (n + 1 + 1), (-1 : ℝ) ^ k * (n + 1).choose k * (harmonic k : ℝ) -
        eulerMascheroniConstant * ∑ k ∈ range (n + 1 + 1), (-1 : ℝ) ^ k * (n + 1).choose k := by
    rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    have : (k.factorial : ℝ) ≠ 0 := by positivity
    field_simp
  rw [hsplit, h1, h0]
  ring

theorem laguerre_zero (x : ℝ) : laguerre 0 x = 1 := by
  simp [laguerre]

/-! ### S10c, the dilogarithm -/

/-- **S10c's premise, from memory.** The pair `(|a|², |b|²)` of unit-variance jointly circular
complex Gaussians with squared coherence `κ < 1` has joint density
`e^{−x−y}·Σ_{n≥0} κⁿ·L_n(x)·L_n(y)` (Hille–Hardy). Integrated term by term against `ln x·ln y`,
that gives `E[ln X·ln Y] = Σ_n κⁿ·(∫₀^∞ ln x·L_n(x)·e^{−x} dx)²`. That consequence is the
premise. -/
def LaguerreProductPremise (κ : ℝ) {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω)
    (X Y : Ω → ℝ) : Prop :=
  HasSum (fun n : ℕ => κ ^ n * (∫ x in Ioi (0 : ℝ), log x * laguerre n x * exp (-x)) ^ 2)
    (∫ w, log (X w) * log (Y w) ∂μ)

/-- **S10c.** Under the premise, and with each marginal unit exponential
(`E[ln X] = ∫₀^∞ ln x·e^{−x} dx`), `Cov(ln X, ln Y) = Σ_{n≥1} κⁿ/n²`. -/
theorem logPower_covariance {κ : ℝ} {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω}
    {X Y : Ω → ℝ} (hprem : LaguerreProductPremise κ μ X Y)
    (hX : ∫ w, log (X w) ∂μ = ∫ x in Ioi (0 : ℝ), log x * exp (-x))
    (hY : ∫ w, log (Y w) ∂μ = ∫ x in Ioi (0 : ℝ), log x * exp (-x)) :
    HasSum (fun n : ℕ => κ ^ (n + 1) / ((n : ℝ) + 1) ^ 2)
      (∫ w, log (X w) * log (Y w) ∂μ - (∫ w, log (X w) ∂μ) * ∫ w, log (Y w) ∂μ) := by
  have h := (hasSum_nat_add_iff' 1).2 hprem
  have hc0 : (∫ x in Ioi (0 : ℝ), log x * laguerre 0 x * exp (-x)) =
      ∫ x in Ioi (0 : ℝ), log x * exp (-x) := by
    simp_rw [laguerre_zero, mul_one]
  simp only [Finset.sum_range_one, pow_zero, one_mul, hc0] at h
  rw [hX, hY, ← sq]
  have hf : (fun n : ℕ => κ ^ (n + 1) / ((n : ℝ) + 1) ^ 2) = fun n : ℕ =>
      κ ^ (n + 1) * (∫ x in Ioi (0 : ℝ), log x * laguerre (n + 1) x * exp (-x)) ^ 2 := by
    funext n
    rw [integral_log_mul_laguerre n, div_pow, neg_one_sq]
    ring
  rw [hf]
  exact h

/-- **S10c's endpoint.** `Σ_{n≥1} 1/n² = π²/6`. -/
theorem dilog_one : HasSum (fun n : ℕ => (1 : ℝ) ^ (n + 1) / ((n : ℝ) + 1) ^ 2) (π ^ 2 / 6) := by
  have h := (hasSum_nat_add_iff' 1).2 hasSum_zeta_two
  simp only [Finset.sum_range_one, Nat.cast_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true,
    zero_pow, div_zero, sub_zero, Nat.cast_add, Nat.cast_one] at h
  simpa using h

/-! ### S11, powers of a non-circular pair -/

/-- **S11.** If the fourth moment of `a` and `b` obeys the complex pairing rule,
`E[aā·bb̄] = E[aā]·E[bb̄] + E[ab̄]·E[āb] + E[ab]·E[āb̄]`, then
`Cov(|a|², |b|²) = |E[a·conj b]|² + |E[a·b]|²`. -/
theorem cov_normSq_pairing {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) (a b : Ω → ℂ)
    (hpair : ∫ w, a w * conj (a w) * (b w * conj (b w)) ∂μ =
      (∫ w, a w * conj (a w) ∂μ) * (∫ w, b w * conj (b w) ∂μ) +
      (∫ w, a w * conj (b w) ∂μ) * (∫ w, conj (a w) * b w ∂μ) +
      (∫ w, a w * b w ∂μ) * (∫ w, conj (a w) * conj (b w) ∂μ)) :
    ∫ w, ‖a w‖ ^ 2 * ‖b w‖ ^ 2 ∂μ - (∫ w, ‖a w‖ ^ 2 ∂μ) * (∫ w, ‖b w‖ ^ 2 ∂μ) =
      ‖∫ w, a w * conj (b w) ∂μ‖ ^ 2 + ‖∫ w, a w * b w ∂μ‖ ^ 2 := by
  have hA : ((∫ w, ‖a w‖ ^ 2 ∂μ : ℝ) : ℂ) = ∫ w, a w * conj (a w) ∂μ := by
    rw [← integral_complex_ofReal]
    congr 1
    funext w
    rw [Complex.mul_conj']
    push_cast
    ring
  have hB : ((∫ w, ‖b w‖ ^ 2 ∂μ : ℝ) : ℂ) = ∫ w, b w * conj (b w) ∂μ := by
    rw [← integral_complex_ofReal]
    congr 1
    funext w
    rw [Complex.mul_conj']
    push_cast
    ring
  have hAB : ((∫ w, ‖a w‖ ^ 2 * ‖b w‖ ^ 2 ∂μ : ℝ) : ℂ) =
      ∫ w, a w * conj (a w) * (b w * conj (b w)) ∂μ := by
    rw [← integral_complex_ofReal]
    congr 1
    funext w
    rw [Complex.mul_conj', Complex.mul_conj']
    push_cast
    ring
  have hc1 : ∫ w, conj (a w) * b w ∂μ = conj (∫ w, a w * conj (b w) ∂μ) := by
    rw [← integral_conj]
    congr 1
    funext w
    simp [mul_comm]
  have hc2 : ∫ w, conj (a w) * conj (b w) ∂μ = conj (∫ w, a w * b w ∂μ) := by
    rw [← integral_conj]
    congr 1
    funext w
    simp
  apply Complex.ofReal_injective
  rw [Complex.ofReal_add, Complex.ofReal_sub, Complex.ofReal_mul, hAB, hA, hB, hpair, hc1, hc2,
    Complex.mul_conj', Complex.mul_conj']
  push_cast
  ring

/-- The witness at `n = 2`: `H₀ − 2·H₁ + H₂ = 0 − 2 + 3/2 = −1/2`. -/
example : ∑ k ∈ range (1 + 1 + 1), (-1 : ℚ) ^ k * (1 + 1).choose k * harmonic k = -1 / 2 := by
  rw [alternating_choose_harmonic 1]
  norm_num

/-- S10b at `n = 1`: `∫₀^∞ ln x·(1 − x)·e^{−x} dx = −1`. -/
example : ∫ x in Ioi (0 : ℝ), log x * laguerre (0 + 1) x * exp (-x) = -1 := by
  rw [integral_log_mul_laguerre 0]
  norm_num

/-- `∫₀^∞ ln x·e^{−x} dx = −γ`. -/
theorem integral_log_mul_exp : ∫ x in Ioi (0 : ℝ), log x * exp (-x) = -eulerMascheroniConstant := by
  have h := integral_pow_mul_log_mul_exp 0
  simp only [pow_zero, one_mul, Nat.factorial_zero, Nat.cast_one, harmonic_zero, Rat.cast_zero,
    zero_sub] at h
  exact h

/-- S10c's premise is satisfiable: at `κ = 0` a law with `ln X = ln Y = −γ` meets it. It is
degenerate, not the Gaussian pair, but it shows the premise is no contradiction. -/
example := logPower_covariance (κ := 0) (μ := Measure.dirac ())
  (X := fun _ : Unit => exp (-eulerMascheroniConstant))
  (Y := fun _ : Unit => exp (-eulerMascheroniConstant))
  (by
    unfold LaguerreProductPremise
    rw [integral_dirac]
    simp only [log_exp]
    convert hasSum_single (f := fun n : ℕ => (0 : ℝ) ^ n *
      (∫ x in Ioi (0 : ℝ), log x * laguerre n x * exp (-x)) ^ 2) 0
      (fun n hn => by simp [zero_pow hn]) using 1
    simp_rw [pow_zero, one_mul, laguerre_zero, mul_one, integral_log_mul_exp]
    ring)
  (by rw [integral_dirac, log_exp, integral_log_mul_exp])
  (by rw [integral_dirac, log_exp, integral_log_mul_exp])

/-- The six-point law: `−√3` at `0`, `√3` at `5`, `0` elsewhere; under `(1/6)·count` it has
`E x² = 1` and `E x⁴ = 3`. -/
noncomputable def sixPointR (w : Fin 6) : ℝ :=
  if w = 0 then -√3 else if w = 5 then √3 else 0

theorem sixPointR_moments :
    ∫ w, sixPointR w * sixPointR w ∂((1 / 6 : NNReal) • Measure.count) = 1 ∧
      ∫ w, sixPointR w * sixPointR w * (sixPointR w * sixPointR w)
        ∂((1 / 6 : NNReal) • Measure.count) = 3 := by
  constructor
  · rw [integral_smul_nnreal_measure, integral_fintype Integrable.of_finite]
    simp [sixPointR, Fin.sum_univ_succ, Measure.real, NNReal.smul_def]
    norm_num
  · rw [integral_smul_nnreal_measure, integral_fintype Integrable.of_finite]
    simp [sixPointR, Fin.sum_univ_succ, Measure.real, NNReal.smul_def]
    norm_num

/-- S11's premise is met by the six-point law taken as a real complex pair `a = b`: the fourth
moment `3` equals `1·1 + 1·1 + 1·1`, and `Var(x²) = 2 = 1² + 1²`. -/
example := cov_normSq_pairing ((1 / 6 : NNReal) • Measure.count)
  (fun w => (sixPointR w : ℂ)) (fun w => (sixPointR w : ℂ)) (by
    simp only [Complex.conj_ofReal, ← Complex.ofReal_mul, integral_complex_ofReal,
      sixPointR_moments.1, sixPointR_moments.2]
    norm_num)

end TimeFrequency
