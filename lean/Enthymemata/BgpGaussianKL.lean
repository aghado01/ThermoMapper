import Enthymemata.BgpUnknownVariance
import Mathlib.InformationTheory.KullbackLeibler.Basic
import Mathlib.Probability.Distributions.Gaussian.Real
import Mathlib.Probability.Moments.MGFAnalytic
import Mathlib.Probability.Moments.Variance
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
# BGP-4, UV-1c — the Gaussian KL divergence and log-ratio variance in closed form

For one observation with law `N(μ₀, σ₀²)` under the truth and `N(μ, σ²)` under the model, with
`δ = μ − μ₀` and `r = σ₀/σ` (`Protolemmata/BARS/bgp-unknown-variance-lemmas.md`, UV-1):

* `klDiv_gaussianReal`: `KL(N(μ₀,σ₀²) ‖ N(μ,σ²)) = −log r + (r² − 1)/2 + δ² r²/(2σ₀²)`, against
  Mathlib's `InformationTheory.klDiv`.
* `variance_llr_gaussianReal`: `Var_{N(μ₀,σ₀²)}[log p₀/p] = (r² − 1)²/2 + δ² r⁴/σ₀²`, against
  Mathlib's `variance` of the log-likelihood ratio `llr`; `evariance_llr_gaussianReal` is the same
  identity in `ℝ≥0∞`.

The route is the protolemma's: the log ratio is a quadratic in `y = x − μ₀`, so its mean needs
`E y` and `E y²`, and its variance also `E y³ = 0` and `E y⁴ = 3σ₀⁴`. Mathlib has no Gaussian fourth
moment at the pinned version, so `integral_pow_gaussianReal_zero` reads the moments of `N(0, v)` off
its moment-generating function `exp(v t²/2)` (`iteratedDeriv_mgf_zero`), differentiated four times.

Variances enter as `(σ ^ 2).toNNReal`, so σ₀ and σ stay real and the identities compose with UV-1a,
UV-1b and UV-1d (`BgpUnknownVariance.lean`) in the same `(μ₀, σ₀, μ, σ)` parameterization.
`klDiv_le_and_evariance_le_of_slab` is that composition: the protolemma's inclusion of the
`|μ − μ₀|`, `|σ − σ₀|` slab of radius `σ₀ε/6` in the KL neighbourhood `{K ≤ ε², V ≤ ε²}`.
-/

open MeasureTheory ProbabilityTheory InformationTheory Real
open scoped ENNReal NNReal

namespace BGP

/-! ### Moments of `N(0, v)` from its moment-generating function -/

/-- The derivative of `t ↦ exp(c t²/2)`. -/
theorem hasDerivAt_exp_mul_sq_div_two (c t : ℝ) :
    HasDerivAt (fun t => rexp (c * t ^ 2 / 2)) (c * t * rexp (c * t ^ 2 / 2)) t := by
  exact ((((hasDerivAt_pow 2 t).const_mul c).div_const 2).exp).congr_deriv (by norm_num; ring)

/-- Differentiating `p(t)·exp(c t²/2)` gives `(p′(t) + c t p(t))·exp(c t²/2)`. -/
theorem deriv_mul_exp_mul_sq_div_two {p p' : ℝ → ℝ} (hp : ∀ t, HasDerivAt p (p' t) t) (c : ℝ) :
    deriv (fun t => p t * rexp (c * t ^ 2 / 2))
      = fun t => (p' t + c * t * p t) * rexp (c * t ^ 2 / 2) := by
  funext t
  refine (((hp t).mul (hasDerivAt_exp_mul_sq_div_two c t)).deriv).trans ?_
  ring

/-- The first four derivatives of `g(t) = exp(c t²/2)`: `g⁽⁴⁾(t) = (3c² + 6c³t² + c⁴t⁴) g(t)`, and
the lower ones in the same form. -/
theorem iteratedDeriv_exp_mul_sq_div_two (c : ℝ) :
    iteratedDeriv 1 (fun t => rexp (c * t ^ 2 / 2)) = (fun t => c * t * rexp (c * t ^ 2 / 2)) ∧
    iteratedDeriv 2 (fun t => rexp (c * t ^ 2 / 2))
      = (fun t => (c + c ^ 2 * t ^ 2) * rexp (c * t ^ 2 / 2)) ∧
    iteratedDeriv 3 (fun t => rexp (c * t ^ 2 / 2))
      = (fun t => (3 * c ^ 2 * t + c ^ 3 * t ^ 3) * rexp (c * t ^ 2 / 2)) ∧
    iteratedDeriv 4 (fun t => rexp (c * t ^ 2 / 2))
      = (fun t => (3 * c ^ 2 + 6 * c ^ 3 * t ^ 2 + c ^ 4 * t ^ 4) * rexp (c * t ^ 2 / 2)) := by
  have d1 : deriv (fun t => rexp (c * t ^ 2 / 2)) = fun t => c * t * rexp (c * t ^ 2 / 2) := by
    funext t
    exact (hasDerivAt_exp_mul_sq_div_two c t).deriv
  have d2 : deriv (fun t => c * t * rexp (c * t ^ 2 / 2))
      = fun t => (c + c ^ 2 * t ^ 2) * rexp (c * t ^ 2 / 2) := by
    rw [deriv_mul_exp_mul_sq_div_two (p := fun t => c * t) (p' := fun _ => c)
      (fun t => by simpa using (hasDerivAt_id t).const_mul c) c]
    funext t
    ring
  have d3 : deriv (fun t => (c + c ^ 2 * t ^ 2) * rexp (c * t ^ 2 / 2))
      = fun t => (3 * c ^ 2 * t + c ^ 3 * t ^ 3) * rexp (c * t ^ 2 / 2) := by
    rw [deriv_mul_exp_mul_sq_div_two (p := fun t => c + c ^ 2 * t ^ 2)
      (p' := fun t => 2 * c ^ 2 * t)
      (fun t => (((hasDerivAt_pow 2 t).const_mul (c ^ 2)).const_add c).congr_deriv
        (by norm_num; ring)) c]
    funext t
    ring
  have d4 : deriv (fun t => (3 * c ^ 2 * t + c ^ 3 * t ^ 3) * rexp (c * t ^ 2 / 2))
      = fun t => (3 * c ^ 2 + 6 * c ^ 3 * t ^ 2 + c ^ 4 * t ^ 4) * rexp (c * t ^ 2 / 2) := by
    rw [deriv_mul_exp_mul_sq_div_two (p := fun t => 3 * c ^ 2 * t + c ^ 3 * t ^ 3)
      (p' := fun t => 3 * c ^ 2 + 3 * c ^ 3 * t ^ 2)
      (fun t => (((hasDerivAt_id t).const_mul (3 * c ^ 2)).add
        ((hasDerivAt_pow 3 t).const_mul (c ^ 3))).congr_deriv (by norm_num; ring)) c]
    funext t
    ring
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [iteratedDeriv_succ', iteratedDeriv_zero, d1]
  · rw [iteratedDeriv_succ', iteratedDeriv_succ', iteratedDeriv_zero, d1, d2]
  · rw [iteratedDeriv_succ', iteratedDeriv_succ', iteratedDeriv_succ', iteratedDeriv_zero, d1, d2,
      d3]
  · rw [iteratedDeriv_succ', iteratedDeriv_succ', iteratedDeriv_succ', iteratedDeriv_succ',
      iteratedDeriv_zero, d1, d2, d3, d4]

/-- Every power of `x` is integrable under `N(μ, v)`: its moment-generating function is finite
everywhere. -/
theorem integrable_pow_gaussianReal (μ : ℝ) (v : ℝ≥0) (n : ℕ) :
    Integrable (fun x : ℝ => x ^ n) (gaussianReal μ v) := by
  have h : (0 : ℝ) ∈ interior (integrableExpSet id (gaussianReal μ v)) := by simp
  simpa using integrable_pow_of_mem_interior_integrableExpSet h n

/-- The moments of `N(0, v)` are the derivatives at `0` of `t ↦ exp(v t²/2)`. -/
theorem integral_pow_gaussianReal_zero (v : ℝ≥0) (n : ℕ) :
    ∫ x, x ^ n ∂gaussianReal 0 v = iteratedDeriv n (fun t => rexp (v * t ^ 2 / 2)) 0 := by
  have h : (0 : ℝ) ∈ interior (integrableExpSet id (gaussianReal 0 v)) := by simp
  have hm := iteratedDeriv_mgf_zero h n
  rw [mgf_id_gaussianReal] at hm
  simp only [zero_mul, zero_add] at hm
  rw [hm]
  rfl

/-- The first four moments of `N(0, v)`: `0, v, 0, 3v²`. -/
theorem moments_gaussianReal_zero (v : ℝ≥0) :
    ∫ x, x ∂gaussianReal 0 v = 0 ∧ ∫ x, x ^ 2 ∂gaussianReal 0 v = v ∧
    ∫ x, x ^ 3 ∂gaussianReal 0 v = 0 ∧ ∫ x, x ^ 4 ∂gaussianReal 0 v = 3 * (v : ℝ) ^ 2 := by
  obtain ⟨d1, d2, d3, d4⟩ := iteratedDeriv_exp_mul_sq_div_two (v : ℝ)
  refine ⟨?_, ?_, ?_, ?_⟩
  · have := integral_pow_gaussianReal_zero v 1
    simp only [pow_one] at this
    rw [this, d1]
    simp
  · rw [integral_pow_gaussianReal_zero, d2]
    simp
  · rw [integral_pow_gaussianReal_zero, d3]
    simp
  · rw [integral_pow_gaussianReal_zero, d4]
    simp

/-- A quartic in `y` integrated against `N(0, v)`: `a₀ + a₂ v + 3 a₄ v²`. -/
theorem integral_quartic_gaussianReal_zero (v : ℝ≥0) (a₀ a₁ a₂ a₃ a₄ : ℝ) :
    ∫ y, (a₀ + a₁ * y + a₂ * y ^ 2 + a₃ * y ^ 3 + a₄ * y ^ 4) ∂gaussianReal 0 v
      = a₀ + a₂ * v + 3 * a₄ * (v : ℝ) ^ 2 := by
  obtain ⟨m1, m2, m3, m4⟩ := moments_gaussianReal_zero v
  have i1 : Integrable (fun y : ℝ => y) (gaussianReal 0 v) := by
    simpa using integrable_pow_gaussianReal 0 v 1
  have i2 := integrable_pow_gaussianReal 0 v 2
  have i3 := integrable_pow_gaussianReal 0 v 3
  have i4 := integrable_pow_gaussianReal 0 v 4
  have j1 : Integrable (fun y : ℝ => a₀ + a₁ * y) (gaussianReal 0 v) :=
    (integrable_const a₀).add (i1.const_mul a₁)
  have j2 : Integrable (fun y : ℝ => a₀ + a₁ * y + a₂ * y ^ 2) (gaussianReal 0 v) :=
    j1.add (i2.const_mul a₂)
  have j3 : Integrable (fun y : ℝ => a₀ + a₁ * y + a₂ * y ^ 2 + a₃ * y ^ 3) (gaussianReal 0 v) :=
    j2.add (i3.const_mul a₃)
  rw [integral_add j3 (i4.const_mul a₄), integral_add j2 (i3.const_mul a₃),
    integral_add j1 (i2.const_mul a₂), integral_add (integrable_const a₀) (i1.const_mul a₁),
    integral_const, integral_const_mul, integral_const_mul, integral_const_mul,
    integral_const_mul, m1, m2, m3, m4]
  simp
  ring

/-- Recentring: `∫ F(x − μ₀) dN(μ₀, v) = ∫ F dN(0, v)`. -/
theorem integral_comp_sub_gaussianReal (μ₀ : ℝ) (v : ℝ≥0) {F : ℝ → ℝ} (hF : Continuous F) :
    ∫ x, F (x - μ₀) ∂gaussianReal μ₀ v = ∫ y, F y ∂gaussianReal 0 v := by
  have hmap : (gaussianReal μ₀ v).map (fun x => x - μ₀) = gaussianReal 0 v := by
    rw [gaussianReal_map_sub_const, sub_self]
  rw [← hmap, integral_map (measurable_sub_const μ₀).aemeasurable hF.aestronglyMeasurable]

/-- Recentring preserves integrability. -/
theorem integrable_comp_sub_gaussianReal (μ₀ : ℝ) (v : ℝ≥0) {F : ℝ → ℝ}
    (hF : Integrable F (gaussianReal 0 v)) :
    Integrable (fun x => F (x - μ₀)) (gaussianReal μ₀ v) := by
  have hmap : (gaussianReal μ₀ v).map (fun x => x - μ₀) = gaussianReal 0 v := by
    rw [gaussianReal_map_sub_const, sub_self]
  rw [← hmap] at hF
  exact (integrable_map_measure hF.aestronglyMeasurable
    (measurable_sub_const μ₀).aemeasurable).1 hF

/-- A quartic in `y` is integrable under `N(0, v)`. -/
theorem integrable_quartic_gaussianReal_zero (v : ℝ≥0) (a₀ a₁ a₂ a₃ a₄ : ℝ) :
    Integrable (fun y : ℝ => a₀ + a₁ * y + a₂ * y ^ 2 + a₃ * y ^ 3 + a₄ * y ^ 4)
      (gaussianReal 0 v) := by
  have i1 : Integrable (fun y : ℝ => y) (gaussianReal 0 v) := by
    simpa using integrable_pow_gaussianReal 0 v 1
  exact ((((integrable_const a₀).add (i1.const_mul a₁)).add
    ((integrable_pow_gaussianReal 0 v 2).const_mul a₂)).add
    ((integrable_pow_gaussianReal 0 v 3).const_mul a₃)).add
    ((integrable_pow_gaussianReal 0 v 4).const_mul a₄)

/-! ### The log-likelihood ratio of two Gaussians -/

/-- `log p₀(x)/p(x)` for `p₀ = N(μ₀, σ₀²)` and `p = N(μ, σ²)`, written out. -/
noncomputable def gaussianLogRatio (μ₀ σ₀ μ σ x : ℝ) : ℝ :=
  -Real.log (σ₀ / σ) - (x - μ₀) ^ 2 / (2 * σ₀ ^ 2) + (x - μ) ^ 2 / (2 * σ ^ 2)

theorem continuous_gaussianLogRatio (μ₀ σ₀ μ σ : ℝ) : Continuous (gaussianLogRatio μ₀ σ₀ μ σ) := by
  unfold gaussianLogRatio
  fun_prop

/-- The log of the Gaussian density, with the variance written as `σ²`. -/
theorem log_gaussianPDFReal (μ : ℝ) {σ : ℝ} (hσ : 0 < σ) (x : ℝ) :
    Real.log (gaussianPDFReal μ (σ ^ 2).toNNReal x)
      = -Real.log (√(2 * π)) - Real.log σ - (x - μ) ^ 2 / (2 * σ ^ 2) := by
  have hpi : (0 : ℝ) ≤ 2 * π := by positivity
  have hsq : √(2 * π * σ ^ 2) = √(2 * π) * σ := by
    rw [Real.sqrt_mul hpi, Real.sqrt_sq hσ.le]
  have hs0 : √(2 * π) ≠ 0 := (Real.sqrt_pos.2 (by positivity)).ne'
  rw [gaussianPDFReal, Real.coe_toNNReal _ (sq_nonneg σ), hsq,
    Real.log_mul (inv_ne_zero (mul_ne_zero hs0 hσ.ne')) (Real.exp_pos _).ne', Real.log_inv,
    Real.log_mul hs0 hσ.ne', Real.log_exp]
  ring

/-- The log ratio of the two densities is `gaussianLogRatio`. -/
theorem log_gaussianPDFReal_div (μ₀ μ : ℝ) {σ₀ σ : ℝ} (h0 : 0 < σ₀) (h : 0 < σ) (x : ℝ) :
    Real.log (gaussianPDFReal μ₀ (σ₀ ^ 2).toNNReal x / gaussianPDFReal μ (σ ^ 2).toNNReal x)
      = gaussianLogRatio μ₀ σ₀ μ σ x := by
  have hv₀ : (σ₀ ^ 2).toNNReal ≠ 0 := by simp [h0.ne']
  have hv : (σ ^ 2).toNNReal ≠ 0 := by simp [h.ne']
  rw [Real.log_div (gaussianPDFReal_pos _ _ _ hv₀).ne' (gaussianPDFReal_pos _ _ _ hv).ne',
    log_gaussianPDFReal μ₀ h0, log_gaussianPDFReal μ h, gaussianLogRatio,
    Real.log_div h0.ne' h.ne']
  ring

/-- `llr N(μ₀,σ₀²) N(μ,σ²)` is `gaussianLogRatio` almost everywhere under the truth. -/
theorem llr_gaussianReal_ae (μ₀ μ : ℝ) {σ₀ σ : ℝ} (h0 : 0 < σ₀) (h : 0 < σ) :
    llr (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) (gaussianReal μ (σ ^ 2).toNNReal)
      =ᵐ[gaussianReal μ₀ (σ₀ ^ 2).toNNReal] gaussianLogRatio μ₀ σ₀ μ σ := by
  have hv₀ : (σ₀ ^ 2).toNNReal ≠ 0 := by simp [h0.ne']
  have hv : (σ ^ 2).toNNReal ≠ 0 := by simp [h.ne']
  have hP₀vol := gaussianReal_absolutelyContinuous μ₀ hv₀
  have hPvol := gaussianReal_absolutelyContinuous μ hv
  have hP₀P := hP₀vol.trans (gaussianReal_absolutelyContinuous' μ hv)
  filter_upwards [hP₀P.ae_le (Measure.rnDeriv_eq_div hP₀vol hPvol),
    hP₀vol.ae_le (rnDeriv_gaussianReal μ₀ (σ₀ ^ 2).toNNReal),
    hP₀vol.ae_le (rnDeriv_gaussianReal μ (σ ^ 2).toNNReal)] with x h1 h2 h3
  rw [llr_def]
  simp only
  rw [h1, h2, h3, ENNReal.toReal_div, toReal_gaussianPDF, toReal_gaussianPDF,
    log_gaussianPDFReal_div μ₀ μ h0 h]

/-- The log ratio as a quadratic in `y = x − μ₀`, with `δ = μ − μ₀`. -/
theorem gaussianLogRatio_eq (μ₀ σ₀ μ σ x : ℝ) :
    gaussianLogRatio μ₀ σ₀ μ σ x
      = (-Real.log (σ₀ / σ) + (μ - μ₀) ^ 2 / (2 * σ ^ 2)) + (-(μ - μ₀) / σ ^ 2) * (x - μ₀)
        + (1 / (2 * σ ^ 2) - 1 / (2 * σ₀ ^ 2)) * (x - μ₀) ^ 2 + 0 * (x - μ₀) ^ 3
        + 0 * (x - μ₀) ^ 4 := by
  unfold gaussianLogRatio
  ring

/-- The mean of the log ratio under the truth, before the algebra into `r = σ₀/σ`. -/
theorem integral_gaussianLogRatio (μ₀ μ : ℝ) {σ₀ σ : ℝ} :
    ∫ x, gaussianLogRatio μ₀ σ₀ μ σ x ∂gaussianReal μ₀ (σ₀ ^ 2).toNNReal
      = (-Real.log (σ₀ / σ) + (μ - μ₀) ^ 2 / (2 * σ ^ 2))
        + (1 / (2 * σ ^ 2) - 1 / (2 * σ₀ ^ 2)) * σ₀ ^ 2 := by
  have hF : Continuous fun y : ℝ =>
      (-Real.log (σ₀ / σ) + (μ - μ₀) ^ 2 / (2 * σ ^ 2)) + (-(μ - μ₀) / σ ^ 2) * y
        + (1 / (2 * σ ^ 2) - 1 / (2 * σ₀ ^ 2)) * y ^ 2 + 0 * y ^ 3 + 0 * y ^ 4 := by
    fun_prop
  simp_rw [gaussianLogRatio_eq]
  rw [integral_comp_sub_gaussianReal μ₀ _ hF, integral_quartic_gaussianReal_zero,
    Real.coe_toNNReal _ (sq_nonneg σ₀)]
  ring

theorem integrable_gaussianLogRatio (μ₀ σ₀ μ σ : ℝ) :
    Integrable (gaussianLogRatio μ₀ σ₀ μ σ) (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) := by
  have h := integrable_comp_sub_gaussianReal μ₀ (σ₀ ^ 2).toNNReal
    (integrable_quartic_gaussianReal_zero (σ₀ ^ 2).toNNReal
      (-Real.log (σ₀ / σ) + (μ - μ₀) ^ 2 / (2 * σ ^ 2)) (-(μ - μ₀) / σ ^ 2)
      (1 / (2 * σ ^ 2) - 1 / (2 * σ₀ ^ 2)) 0 0)
  refine h.congr (ae_of_all _ fun x => ?_)
  exact (gaussianLogRatio_eq μ₀ σ₀ μ σ x).symm

/-! ### UV-1c: the closed forms -/

/-- **UV-1c (KL).** `KL(N(μ₀,σ₀²) ‖ N(μ,σ²)) = −log r + (r² − 1)/2 + δ² r²/(2σ₀²)` with
`r = σ₀/σ` and `δ = μ − μ₀`. The log ratio is integrable and, by `gaussianLogRatio_eq`, a quadratic
in `x − μ₀`; its mean needs only the first two moments of `N(0, σ₀²)`. -/
theorem klDiv_gaussianReal (μ₀ μ : ℝ) {σ₀ σ : ℝ} (h0 : 0 < σ₀) (h : 0 < σ) :
    klDiv (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) (gaussianReal μ (σ ^ 2).toNNReal)
      = ENNReal.ofReal (-Real.log (σ₀ / σ) + ((σ₀ / σ) ^ 2 - 1) / 2
          + (μ - μ₀) ^ 2 * (σ₀ / σ) ^ 2 / (2 * σ₀ ^ 2)) := by
  have hv : (σ ^ 2).toNNReal ≠ 0 := by simp [h.ne']
  have hv₀ : (σ₀ ^ 2).toNNReal ≠ 0 := by simp [h0.ne']
  have hP₀P := (gaussianReal_absolutelyContinuous μ₀ hv₀).trans
    (gaussianReal_absolutelyContinuous' μ hv)
  have hae := llr_gaussianReal_ae μ₀ μ h0 h
  have hint : Integrable (llr (gaussianReal μ₀ (σ₀ ^ 2).toNNReal)
      (gaussianReal μ (σ ^ 2).toNNReal)) (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) :=
    (integrable_gaussianLogRatio μ₀ σ₀ μ σ).congr hae.symm
  rw [klDiv_of_ac_of_integrable hP₀P hint, integral_congr_ae hae, integral_gaussianLogRatio]
  simp only [probReal_univ, add_sub_cancel_right]
  congr 1
  field_simp
  ring

/-- The variance of the log ratio under the truth, before the algebra into `r = σ₀/σ`. -/
theorem variance_gaussianLogRatio (μ₀ μ : ℝ) {σ₀ σ : ℝ} :
    Var[gaussianLogRatio μ₀ σ₀ μ σ; gaussianReal μ₀ (σ₀ ^ 2).toNNReal]
      = ((μ - μ₀) / σ ^ 2) ^ 2 * σ₀ ^ 2
        + 2 * (1 / (2 * σ ^ 2) - 1 / (2 * σ₀ ^ 2)) ^ 2 * (σ₀ ^ 2) ^ 2 := by
  set A := -Real.log (σ₀ / σ) + (μ - μ₀) ^ 2 / (2 * σ ^ 2)
  set B := -(μ - μ₀) / σ ^ 2
  set C := 1 / (2 * σ ^ 2) - 1 / (2 * σ₀ ^ 2)
  set v : ℝ := σ₀ ^ 2 with hv_def
  have hcoe : ((σ₀ ^ 2).toNNReal : ℝ) = v := Real.coe_toNNReal _ (sq_nonneg σ₀)
  have hsq : ∀ x, (gaussianLogRatio μ₀ σ₀ μ σ x - (A + C * v)) ^ 2
      = (C ^ 2 * v ^ 2) + (-2 * B * C * v) * (x - μ₀) + (B ^ 2 - 2 * C ^ 2 * v) * (x - μ₀) ^ 2
        + (2 * B * C) * (x - μ₀) ^ 3 + C ^ 2 * (x - μ₀) ^ 4 := by
    intro x
    rw [gaussianLogRatio_eq]
    ring
  have hF : Continuous fun y : ℝ =>
      (C ^ 2 * v ^ 2) + (-2 * B * C * v) * y + (B ^ 2 - 2 * C ^ 2 * v) * y ^ 2
        + (2 * B * C) * y ^ 3 + C ^ 2 * y ^ 4 := by
    fun_prop
  have hmean : ∫ x, gaussianLogRatio μ₀ σ₀ μ σ x ∂gaussianReal μ₀ (σ₀ ^ 2).toNNReal
      = A + C * v := by
    rw [integral_gaussianLogRatio]
  rw [variance_eq_integral (continuous_gaussianLogRatio μ₀ σ₀ μ σ).measurable.aemeasurable, hmean]
  simp_rw [hsq]
  rw [integral_comp_sub_gaussianReal μ₀ _ hF, integral_quartic_gaussianReal_zero, hcoe]
  have hB : ((μ - μ₀) / σ ^ 2) ^ 2 = B ^ 2 := by ring
  rw [hB]
  ring

/-- **UV-1c (variance).** `Var_{N(μ₀,σ₀²)}[log p₀/p] = (r² − 1)²/2 + δ² r⁴/σ₀²` with `r = σ₀/σ` and
`δ = μ − μ₀`, for Mathlib's log-likelihood ratio `llr`. The variance of a quadratic in `x − μ₀`
needs the third and fourth moments of `N(0, σ₀²)`, `0` and `3σ₀⁴`. -/
theorem variance_llr_gaussianReal (μ₀ μ : ℝ) {σ₀ σ : ℝ} (h0 : 0 < σ₀) (h : 0 < σ) :
    Var[llr (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) (gaussianReal μ (σ ^ 2).toNNReal);
        gaussianReal μ₀ (σ₀ ^ 2).toNNReal]
      = ((σ₀ / σ) ^ 2 - 1) ^ 2 / 2 + (μ - μ₀) ^ 2 * (σ₀ / σ) ^ 4 / σ₀ ^ 2 := by
  rw [variance_congr (llr_gaussianReal_ae μ₀ μ h0 h), variance_gaussianLogRatio]
  field_simp
  ring

/-- The log ratio has a finite second moment under the truth. -/
theorem memLp_two_gaussianLogRatio (μ₀ σ₀ μ σ : ℝ) :
    MemLp (gaussianLogRatio μ₀ σ₀ μ σ) 2 (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) := by
  rw [memLp_two_iff_integrable_sq
    (continuous_gaussianLogRatio μ₀ σ₀ μ σ).aestronglyMeasurable]
  have h := integrable_comp_sub_gaussianReal μ₀ (σ₀ ^ 2).toNNReal
    (integrable_quartic_gaussianReal_zero (σ₀ ^ 2).toNNReal
      ((-Real.log (σ₀ / σ) + (μ - μ₀) ^ 2 / (2 * σ ^ 2)) ^ 2)
      (2 * (-Real.log (σ₀ / σ) + (μ - μ₀) ^ 2 / (2 * σ ^ 2)) * (-(μ - μ₀) / σ ^ 2))
      ((-(μ - μ₀) / σ ^ 2) ^ 2 + 2 * (-Real.log (σ₀ / σ) + (μ - μ₀) ^ 2 / (2 * σ ^ 2))
        * (1 / (2 * σ ^ 2) - 1 / (2 * σ₀ ^ 2)))
      (2 * (-(μ - μ₀) / σ ^ 2) * (1 / (2 * σ ^ 2) - 1 / (2 * σ₀ ^ 2)))
      ((1 / (2 * σ ^ 2) - 1 / (2 * σ₀ ^ 2)) ^ 2))
  refine h.congr (ae_of_all _ fun x => ?_)
  dsimp only
  rw [gaussianLogRatio_eq]
  ring

/-- **UV-1c (variance, in `ℝ≥0∞`).** The same identity for `evariance`, the form the
KL neighbourhood `{K ≤ ε², V ≤ ε²}` takes without junk values. -/
theorem evariance_llr_gaussianReal (μ₀ μ : ℝ) {σ₀ σ : ℝ} (h0 : 0 < σ₀) (h : 0 < σ) :
    evariance (llr (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) (gaussianReal μ (σ ^ 2).toNNReal))
        (gaussianReal μ₀ (σ₀ ^ 2).toNNReal)
      = ENNReal.ofReal (((σ₀ / σ) ^ 2 - 1) ^ 2 / 2 + (μ - μ₀) ^ 2 * (σ₀ / σ) ^ 4 / σ₀ ^ 2) := by
  rw [← variance_llr_gaussianReal μ₀ μ h0 h, evariance_congr (llr_gaussianReal_ae μ₀ μ h0 h),
    variance_congr (llr_gaussianReal_ae μ₀ μ h0 h),
    ofReal_variance (memLp_two_gaussianLogRatio μ₀ σ₀ μ σ)]

/-! ### UV-1 assembled: the slab lies in the KL neighbourhood -/

/-- **UV-1 (the inclusion).** A slab `|μ − μ₀| ≤ σ₀ε/6`, `|σ − σ₀| ≤ σ₀ε/6` with `0 ≤ ε ≤ 3` lies
in the KL neighbourhood: `K ≤ ε²` and `V ≤ ε²`. UV-1d puts `r = σ₀/σ` in `[1/2, 2]` with
`(r − 1)² ≤ 4(σ − σ₀)²/σ₀²`; then UV-1a and UV-1b bound the closed forms of UV-1c by
`K ≤ 10(σ − σ₀)²/σ₀² + 2δ²/σ₀² ≤ ε²/3` and `V ≤ 18(σ − σ₀)²/σ₀² + 16δ²/σ₀² ≤ 17ε²/18`. -/
theorem klDiv_le_and_evariance_le_of_slab {μ₀ μ σ₀ σ ε : ℝ} (h0 : 0 < σ₀) (hε0 : 0 ≤ ε)
    (hε : ε ≤ 3) (hδ : |μ - μ₀| ≤ σ₀ * ε / 6) (hσ : |σ - σ₀| ≤ σ₀ * ε / 6) :
    klDiv (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) (gaussianReal μ (σ ^ 2).toNNReal)
        ≤ ENNReal.ofReal (ε ^ 2) ∧
      evariance (llr (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) (gaussianReal μ (σ ^ 2).toNNReal))
        (gaussianReal μ₀ (σ₀ ^ 2).toNNReal) ≤ ENNReal.ofReal (ε ^ 2) := by
  have hσ2 : |σ - σ₀| ≤ σ₀ / 2 := hσ.trans (by nlinarith)
  obtain ⟨hr, hr1⟩ := scale_ratio_bridge h0 hσ2
  have hσpos : 0 < σ := by
    have := (abs_le.1 hσ2).1
    linarith
  have hs0 : 0 < σ₀ ^ 2 := by positivity
  have hD : (σ - σ₀) ^ 2 / σ₀ ^ 2 ≤ ε ^ 2 / 36 := by
    rw [div_le_iff₀ hs0]
    have h := pow_le_pow_left₀ (abs_nonneg _) hσ 2
    rw [sq_abs] at h
    nlinarith
  have hΔ : (μ - μ₀) ^ 2 / σ₀ ^ 2 ≤ ε ^ 2 / 36 := by
    rw [div_le_iff₀ hs0]
    have h := pow_le_pow_left₀ (abs_nonneg _) hδ 2
    rw [sq_abs] at h
    nlinarith
  have hΔ0 : 0 ≤ (μ - μ₀) ^ 2 / σ₀ ^ 2 := by positivity
  have hr1' : (σ₀ / σ - 1) ^ 2 ≤ 4 * ((σ - σ₀) ^ 2 / σ₀ ^ 2) := by rwa [mul_div_assoc] at hr1
  obtain ⟨hrlo, hrhi⟩ := hr
  have hr2 : (σ₀ / σ) ^ 2 ≤ 4 := by nlinarith
  have hr4 : (σ₀ / σ) ^ 4 ≤ 16 := by nlinarith
  have hK1 := gaussian_kl_scale_part_le (σ₀ / σ) hrlo
  have hV1 := gaussian_v_scale_part_le (σ₀ / σ) ⟨hrlo, hrhi⟩
  refine ⟨?_, ?_⟩
  · rw [klDiv_gaussianReal μ₀ μ h0 hσpos]
    apply ENNReal.ofReal_le_ofReal
    have e : (μ - μ₀) ^ 2 * (σ₀ / σ) ^ 2 / (2 * σ₀ ^ 2)
        = (σ₀ / σ) ^ 2 / 2 * ((μ - μ₀) ^ 2 / σ₀ ^ 2) := by
      field_simp
    have hδK : (σ₀ / σ) ^ 2 / 2 * ((μ - μ₀) ^ 2 / σ₀ ^ 2) ≤ 2 * ((μ - μ₀) ^ 2 / σ₀ ^ 2) :=
      mul_le_mul_of_nonneg_right (by linarith) hΔ0
    rw [e]
    nlinarith
  · rw [evariance_llr_gaussianReal μ₀ μ h0 hσpos]
    apply ENNReal.ofReal_le_ofReal
    have e : (μ - μ₀) ^ 2 * (σ₀ / σ) ^ 4 / σ₀ ^ 2 = (σ₀ / σ) ^ 4 * ((μ - μ₀) ^ 2 / σ₀ ^ 2) := by
      field_simp
    have hδV : (σ₀ / σ) ^ 4 * ((μ - μ₀) ^ 2 / σ₀ ^ 2) ≤ 16 * ((μ - μ₀) ^ 2 / σ₀ ^ 2) :=
      mul_le_mul_of_nonneg_right hr4 hΔ0
    rw [e]
    nlinarith

end BGP
