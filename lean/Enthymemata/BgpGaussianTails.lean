import Mathlib.MeasureTheory.Integral.Lebesgue.Markov
import Mathlib.Probability.Distributions.Gaussian.Real
import Mathlib.Probability.Independence.Integration
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# BGP-4, UV-5 — Gaussian tail bounds for the tests

The analytic core of the tests in `Protolemmata/BARS/bgp-unknown-variance-lemmas.md` (UV-5), for
`n` independent observations `Yᵢ ~ N(mᵢ, v)`, as `Measure.pi`:

* `lintegral_exp_quadratic_gaussianReal`: completing the square,
  `E exp(a(Y − c)² + b(Y − c)) = (1 − 2av)^(−1/2) exp((a(m − c)² + b(m − c) + b²v/2)/(1 − 2av))`
  for `2av < 1`, from the density identity and `lintegral_gaussianPDF_eq_one`;
* `lintegral_pi_exp_sum`: under `Measure.pi`, the moment-generating function of a coordinate-wise
  sum factors (`iIndepFun_pi`);
* `measure_ge_le_exp_mul_lintegral`: Chernoff's bound, from Markov's inequality in `ℝ≥0∞`;
* `pi_gaussian_linear_tail`: `P(Σ cᵢ(Yᵢ − mᵢ) ≥ t) ≤ exp(−t²/(2v Σ cᵢ²))`, for the mean tests;
* `pi_gaussian_sq_upper_tail`: with `A = nv + Σ(mᵢ − cᵢ)²` and `nτ > A`,
  `P(Σ(Yᵢ − cᵢ)² ≥ nτ) ≤ exp(−(nτ − A)²/(16 v nτ))`;
* `pi_gaussian_sq_lower_tail`: for `τ < v`, `P(Σ(Yᵢ − cᵢ)² ≤ nτ) ≤ exp(−n(v − τ)²/(8v²))`.

The last two are the χ² tests of the residual scale. Lebesgue integrals keep every step free of
integrability side conditions.
-/

open MeasureTheory ProbabilityTheory Real
open scoped ENNReal NNReal

namespace BGP

/-! ### One Gaussian coordinate: completing the square -/

/-- The density identity behind `lintegral_exp_quadratic_gaussianReal`: tilting `N(m, v)` by
`exp(a(y − c)² + b(y − c))` gives a constant times another Gaussian density. -/
theorem gaussianPDFReal_mul_exp_quadratic (m c a b : ℝ) {v : ℝ} (hv : 0 < v)
    (h : 2 * a * v < 1) (y : ℝ) :
    gaussianPDFReal m v.toNNReal y * rexp (a * (y - c) ^ 2 + b * (y - c))
      = (1 / √(1 - 2 * a * v)
          * rexp ((a * (m - c) ^ 2 + b * (m - c) + b ^ 2 * v / 2) / (1 - 2 * a * v)))
        * gaussianPDFReal (c + ((m - c) + b * v) / (1 - 2 * a * v))
            (v / (1 - 2 * a * v)).toNNReal y := by
  have hw : 0 < 1 - 2 * a * v := by linarith
  have hv' : 0 < v / (1 - 2 * a * v) := div_pos hv hw
  simp only [gaussianPDFReal, Real.coe_toNNReal _ hv.le, Real.coe_toNNReal _ hv'.le]
  have hexp : rexp (-(y - m) ^ 2 / (2 * v)) * rexp (a * (y - c) ^ 2 + b * (y - c))
      = rexp ((a * (m - c) ^ 2 + b * (m - c) + b ^ 2 * v / 2) / (1 - 2 * a * v))
        * rexp (-(y - (c + ((m - c) + b * v) / (1 - 2 * a * v))) ^ 2
            / (2 * (v / (1 - 2 * a * v)))) := by
    rw [← Real.exp_add, ← Real.exp_add]
    congr 1
    obtain ⟨w, hW⟩ : ∃ w, 1 - 2 * a * v = w := ⟨_, rfl⟩
    have hw0 : w ≠ 0 := hW ▸ hw.ne'
    have hv0 : v ≠ 0 := hv.ne'
    rw [hW]
    field_simp
    rw [← hW]
    ring
  have hsq : √(2 * π * (v / (1 - 2 * a * v))) = √(2 * π * v) / √(1 - 2 * a * v) := by
    rw [show 2 * π * (v / (1 - 2 * a * v)) = (2 * π * v) / (1 - 2 * a * v) by ring,
      Real.sqrt_div' _ hw.le]
  rw [hsq, mul_assoc, hexp]
  obtain ⟨s, hs⟩ : ∃ s, √(2 * π * v) = s := ⟨_, rfl⟩
  obtain ⟨t, ht⟩ : ∃ t, √(1 - 2 * a * v) = t := ⟨_, rfl⟩
  have hs0 : s ≠ 0 := hs ▸ (Real.sqrt_pos.2 (by positivity)).ne'
  have ht0 : t ≠ 0 := ht ▸ (Real.sqrt_pos.2 hw).ne'
  rw [hs, ht]
  set E₁ := rexp ((a * (m - c) ^ 2 + b * (m - c) + b ^ 2 * v / 2) / (1 - 2 * a * v))
  set E₂ := rexp (-(y - (c + ((m - c) + b * v) / (1 - 2 * a * v))) ^ 2
    / (2 * (v / (1 - 2 * a * v))))
  rw [inv_div]
  field_simp

/-- **Completing the square.** For `2av < 1`,
`E exp(a(Y − c)² + b(Y − c)) = (1 − 2av)^(−1/2) exp((a(m − c)² + b(m − c) + b²v/2)/(1 − 2av))`
under `Y ~ N(m, v)`. -/
theorem lintegral_exp_quadratic_gaussianReal (m c a b : ℝ) {v : ℝ} (hv : 0 < v)
    (h : 2 * a * v < 1) :
    ∫⁻ y, ENNReal.ofReal (rexp (a * (y - c) ^ 2 + b * (y - c))) ∂gaussianReal m v.toNNReal
      = ENNReal.ofReal (1 / √(1 - 2 * a * v)
          * rexp ((a * (m - c) ^ 2 + b * (m - c) + b ^ 2 * v / 2) / (1 - 2 * a * v))) := by
  have hw : 0 < 1 - 2 * a * v := by linarith
  have hv0 : v.toNNReal ≠ 0 := by simp [hv]
  have hv0' : (v / (1 - 2 * a * v)).toNNReal ≠ 0 := by simp [div_pos hv hw]
  have hC : 0 ≤ 1 / √(1 - 2 * a * v)
      * rexp ((a * (m - c) ^ 2 + b * (m - c) + b ^ 2 * v / 2) / (1 - 2 * a * v)) := by
    positivity
  have hmeas : Measurable fun y => ENNReal.ofReal (rexp (a * (y - c) ^ 2 + b * (y - c))) :=
    ENNReal.measurable_ofReal.comp (Real.measurable_exp.comp (by fun_prop))
  rw [gaussianReal_of_var_ne_zero _ hv0,
    lintegral_withDensity_eq_lintegral_mul _ (measurable_gaussianPDF _ _) hmeas]
  have hpt : ∀ y, (gaussianPDF m v.toNNReal * fun y =>
        ENNReal.ofReal (rexp (a * (y - c) ^ 2 + b * (y - c)))) y
      = ENNReal.ofReal (1 / √(1 - 2 * a * v)
          * rexp ((a * (m - c) ^ 2 + b * (m - c) + b ^ 2 * v / 2) / (1 - 2 * a * v)))
        * gaussianPDF (c + ((m - c) + b * v) / (1 - 2 * a * v))
            (v / (1 - 2 * a * v)).toNNReal y := by
    intro y
    simp only [Pi.mul_apply, gaussianPDF]
    rw [← ENNReal.ofReal_mul (gaussianPDFReal_nonneg _ _ _), ← ENNReal.ofReal_mul hC,
      gaussianPDFReal_mul_exp_quadratic m c a b hv h y]
  rw [lintegral_congr hpt, lintegral_const_mul _ (measurable_gaussianPDF _ _),
    lintegral_gaussianPDF_eq_one _ hv0', mul_one]

/-! ### Independent coordinates -/

/-- Under `Measure.pi` of Gaussians, `E exp(Σᵢ gᵢ(Yᵢ))` factors into the coordinate integrals. -/
theorem lintegral_pi_exp_sum {n : ℕ} (m : Fin n → ℝ) (v : ℝ≥0) (g : Fin n → ℝ → ℝ)
    (hg : ∀ i, Measurable (g i)) :
    ∫⁻ y, ENNReal.ofReal (rexp (∑ i, g i (y i))) ∂Measure.pi (fun i => gaussianReal (m i) v)
      = ∏ i, ∫⁻ x, ENNReal.ofReal (rexp (g i x)) ∂gaussianReal (m i) v := by
  have hpt : ∀ y : Fin n → ℝ, ENNReal.ofReal (rexp (∑ i, g i (y i)))
      = ∏ i, ENNReal.ofReal (rexp (g i (y i))) := by
    intro y
    rw [Real.exp_sum, ENNReal.ofReal_prod_of_nonneg (fun i _ => (Real.exp_pos _).le)]
  have hmeas : ∀ i, Measurable fun x => ENNReal.ofReal (rexp (g i x)) :=
    fun i => ENNReal.measurable_ofReal.comp (Real.measurable_exp.comp (hg i))
  have hind := lintegral_prod_eq_prod_lintegral_of_indepFun Finset.univ
    (fun i (ω : Fin n → ℝ) => ENNReal.ofReal (rexp (g i (ω i))))
    (iIndepFun_pi (μ := fun i => gaussianReal (m i) v)
      (X := fun i x => ENNReal.ofReal (rexp (g i x))) (fun i => (hmeas i).aemeasurable))
    (fun i => (hmeas i).comp (measurable_pi_apply i))
  rw [lintegral_congr hpt, hind]
  refine Finset.prod_congr rfl fun i _ => ?_
  exact (measurePreserving_eval (fun i => gaussianReal (m i) v) i).lintegral_comp (hmeas i)

/-! ### Chernoff -/

/-- **Chernoff's bound.** `P(X ≥ t) ≤ e^{−κt} E e^{κX}` for `κ ≥ 0`, by Markov's inequality applied
to `e^{κX}`. -/
theorem measure_ge_le_exp_mul_lintegral {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    {X : Ω → ℝ} (hX : Measurable X) {κ : ℝ} (hκ : 0 ≤ κ) (t : ℝ) :
    P {ω | t ≤ X ω} ≤ ENNReal.ofReal (rexp (-(κ * t)))
      * ∫⁻ ω, ENNReal.ofReal (rexp (κ * X ω)) ∂P := by
  have hm : Measurable fun ω => ENNReal.ofReal (rexp (κ * X ω)) :=
    ENNReal.measurable_ofReal.comp (Real.measurable_exp.comp (hX.const_mul κ))
  have hmk := mul_meas_ge_le_lintegral₀ (hm.aemeasurable (μ := P)) (ENNReal.ofReal (rexp (κ * t)))
  have hsub : {ω | t ≤ X ω}
      ⊆ {ω | ENNReal.ofReal (rexp (κ * t)) ≤ ENNReal.ofReal (rexp (κ * X ω))} := by
    intro ω hω
    exact ENNReal.ofReal_le_ofReal (Real.exp_le_exp.2 (mul_le_mul_of_nonneg_left hω hκ))
  have hone : ENNReal.ofReal (rexp (-(κ * t))) * ENNReal.ofReal (rexp (κ * t)) = 1 := by
    rw [← ENNReal.ofReal_mul (Real.exp_pos _).le, ← Real.exp_add, neg_add_cancel, Real.exp_zero,
      ENNReal.ofReal_one]
  calc P {ω | t ≤ X ω}
      = ENNReal.ofReal (rexp (-(κ * t))) * (ENNReal.ofReal (rexp (κ * t)) * P {ω | t ≤ X ω}) := by
        rw [← mul_assoc, hone, one_mul]
    _ ≤ ENNReal.ofReal (rexp (-(κ * t)))
        * (ENNReal.ofReal (rexp (κ * t))
          * P {ω | ENNReal.ofReal (rexp (κ * t)) ≤ ENNReal.ofReal (rexp (κ * X ω))}) :=
        mul_le_mul_right (mul_le_mul_right (measure_mono hsub) _) _
    _ ≤ _ := mul_le_mul_right hmk _

/-! ### Elementary bounds on the square-root factor -/

/-- `(1 − x)^(−1/2) ≤ exp((x/2)/(1 − x))` for `x < 1`, from `log z ≤ z − 1` at
`z = 1/(1 − x)`. -/
theorem one_div_sqrt_one_sub_le_exp {x : ℝ} (hx1 : x < 1) :
    1 / √(1 - x) ≤ rexp (x / 2 / (1 - x)) := by
  have hw : 0 < 1 - x := by linarith
  have hz : 0 < 1 / (1 - x) := by positivity
  have h1 : 1 / √(1 - x) = √(1 / (1 - x)) := by
    rw [Real.sqrt_div' _ hw.le, Real.sqrt_one]
  rw [h1, ← Real.exp_log (Real.sqrt_pos.2 hz), Real.log_sqrt hz.le]
  apply Real.exp_le_exp.2
  have h2 := Real.log_le_sub_one_of_pos hz
  have h3 : 1 / (1 - x) - 1 = x / (1 - x) := by
    rw [div_sub_one hw.ne', sub_sub_cancel]
  have h4 : x / 2 / (1 - x) = (x / (1 - x)) / 2 := by ring
  rw [h4]
  linarith

/-- `(1 + x)^(−1/2) ≤ exp(−(x/2 − x²/2))` for `x ≥ 0`, from
`log(1 + x) ≥ 1 − 1/(1 + x) ≥ x − x²`. -/
theorem one_div_sqrt_one_add_le_exp {x : ℝ} (hx0 : 0 ≤ x) :
    1 / √(1 + x) ≤ rexp (-(x / 2 - x ^ 2 / 2)) := by
  have hw : 0 < 1 + x := by linarith
  have h1 : 1 / √(1 + x) = rexp (-(Real.log (1 + x) / 2)) := by
    rw [Real.exp_neg, ← Real.log_sqrt hw.le, Real.exp_log (Real.sqrt_pos.2 hw), one_div]
  rw [h1]
  apply Real.exp_le_exp.2
  have h2 := Real.one_sub_inv_le_log_of_pos hw
  obtain ⟨u, hu⟩ : ∃ u, (1 + x)⁻¹ = u := ⟨_, rfl⟩
  have hu1 : u * (1 + x) = 1 := by rw [← hu]; exact inv_mul_cancel₀ hw.ne'
  have hu0 : 0 ≤ u := by rw [← hu]; positivity
  rw [hu] at h2
  -- x − x² ≤ 1 − u, since (1 − u − (x − x²))(1 + x) = x³ ≥ 0
  have h3 : x - x ^ 2 ≤ 1 - u := by
    have e : (1 - u - (x - x ^ 2)) * (1 + x) = x ^ 3 := by
      have : u * (1 + x) = 1 := hu1
      nlinarith
    have hx3 : 0 ≤ x ^ 3 := by positivity
    have : 0 ≤ 1 - u - (x - x ^ 2) := by
      by_contra hneg
      have := mul_neg_of_neg_of_pos (not_le.mp hneg) hw
      linarith
    linarith
  linarith

/-! ### The three tails -/

variable {n : ℕ}

/-- **The mean test's tail.** `P(Σ cᵢ(Yᵢ − mᵢ) ≥ t) ≤ exp(−t²/(2v Σ cᵢ²))` for independent
`Yᵢ ~ N(mᵢ, v)` and `t ≥ 0`: Chernoff at `κ = t/(v Σ cᵢ²)`. -/
theorem pi_gaussian_linear_tail (m c : Fin n → ℝ) {v t : ℝ} (hv : 0 < v)
    (hQ : 0 < ∑ i, c i ^ 2) (ht : 0 ≤ t) :
    Measure.pi (fun i => gaussianReal (m i) v.toNNReal) {y | t ≤ ∑ i, c i * (y i - m i)}
      ≤ ENNReal.ofReal (rexp (-(t ^ 2 / (2 * v * ∑ i, c i ^ 2)))) := by
  obtain ⟨Q, hQdef⟩ : ∃ Q, ∑ i, c i ^ 2 = Q := ⟨_, rfl⟩
  rw [hQdef] at hQ ⊢
  obtain ⟨κ, hκ⟩ : ∃ κ, t / (v * Q) = κ := ⟨_, rfl⟩
  have hκ0 : 0 ≤ κ := hκ ▸ div_nonneg ht (mul_pos hv hQ).le
  have hX : Measurable fun y : Fin n → ℝ => ∑ i, c i * (y i - m i) :=
    Finset.measurable_sum _ fun i _ => (measurable_const.mul
      ((measurable_pi_apply i).sub measurable_const))
  refine (measure_ge_le_exp_mul_lintegral _ hX hκ0 t).trans (le_of_eq ?_)
  have hsum : ∀ y : Fin n → ℝ, κ * ∑ i, c i * (y i - m i)
      = ∑ i, (0 * (y i - m i) ^ 2 + κ * c i * (y i - m i)) := by
    intro y
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  simp_rw [hsum]
  rw [lintegral_pi_exp_sum m v.toNNReal (fun i x => 0 * (x - m i) ^ 2 + κ * c i * (x - m i))
    (fun i => by fun_prop)]
  have hfac : ∀ i, ∫⁻ x, ENNReal.ofReal (rexp (0 * (x - m i) ^ 2 + κ * c i * (x - m i)))
      ∂gaussianReal (m i) v.toNNReal = ENNReal.ofReal (rexp (κ ^ 2 * v / 2 * c i ^ 2)) := by
    intro i
    rw [lintegral_exp_quadratic_gaussianReal (m i) (m i) 0 (κ * c i) hv (by linarith)]
    congr 1
    simp only [mul_zero, zero_mul, sub_zero, sub_self, Real.sqrt_one, div_one, one_mul,
      zero_add]
    ring_nf
  rw [Finset.prod_congr rfl fun i _ => hfac i,
    ← ENNReal.ofReal_prod_of_nonneg (fun i _ => (Real.exp_pos _).le), ← Real.exp_sum,
    ← ENNReal.ofReal_mul (Real.exp_pos _).le, ← Real.exp_add, ← Finset.mul_sum, hQdef]
  congr 2
  have hvQ : v * Q ≠ 0 := (mul_pos hv hQ).ne'
  rw [← hκ]
  field_simp
  ring

/-- **The residual scale's upper tail.** For independent `Yᵢ ~ N(mᵢ, v)`, with
`A = nv + Σ(mᵢ − cᵢ)²` and `nτ > A`, `P(Σ(Yᵢ − cᵢ)² ≥ nτ) ≤ exp(−(nτ − A)²/(16 v nτ))`: Chernoff at
`κ = G/(8v(A + G))` with `G = nτ − A`, so that `2κv ≤ 1/4`. -/
theorem pi_gaussian_sq_upper_tail (m c : Fin n → ℝ) {v τ : ℝ} (hv : 0 < v)
    (hG : n * v + ∑ i, (m i - c i) ^ 2 < n * τ) :
    Measure.pi (fun i => gaussianReal (m i) v.toNNReal) {y | n * τ ≤ ∑ i, (y i - c i) ^ 2}
      ≤ ENNReal.ofReal (rexp (-((n * τ - (n * v + ∑ i, (m i - c i) ^ 2)) ^ 2
          / (16 * v * (n * τ))))) := by
  obtain ⟨M, hM⟩ : ∃ M, ∑ i, (m i - c i) ^ 2 = M := ⟨_, rfl⟩
  have hM0 : 0 ≤ M := hM ▸ Finset.sum_nonneg fun i _ => sq_nonneg _
  rw [hM] at hG ⊢
  obtain ⟨A, hA⟩ : ∃ A, n * v + M = A := ⟨_, rfl⟩
  obtain ⟨D, hD⟩ : ∃ D, (n : ℝ) * τ = D := ⟨_, rfl⟩
  rw [hA, hD] at hG ⊢
  have hA0 : 0 ≤ A := hA ▸ by positivity
  have hGpos : 0 < D - A := by linarith
  have hD0 : 0 < D := by linarith
  obtain ⟨κ, hκ⟩ : ∃ κ, (D - A) / (8 * v * D) = κ := ⟨_, rfl⟩
  have hκ0 : 0 < κ := hκ ▸ by positivity
  obtain ⟨x, hx⟩ : ∃ x, 2 * κ * v = x := ⟨_, rfl⟩
  have hx_eq : x = (D - A) / (4 * D) := by
    rw [← hx, ← hκ]
    field_simp
    ring
  have hx0 : 0 ≤ x := by rw [hx_eq]; positivity
  have hx4 : x ≤ 1 / 4 := by
    rw [hx_eq, div_le_iff₀ (by positivity)]
    linarith
  have hX : Measurable fun y : Fin n → ℝ => ∑ i, (y i - c i) ^ 2 :=
    Finset.measurable_sum _ fun i _ => ((measurable_pi_apply i).sub measurable_const).pow_const 2
  rw [← hD]
  refine (measure_ge_le_exp_mul_lintegral _ hX hκ0.le (n * τ)).trans ?_
  rw [hD]
  have hsum : ∀ y : Fin n → ℝ, κ * ∑ i, (y i - c i) ^ 2
      = ∑ i, (κ * (y i - c i) ^ 2 + 0 * (y i - c i)) := by
    intro y
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  simp_rw [hsum]
  rw [lintegral_pi_exp_sum m v.toNNReal (fun i x => κ * (x - c i) ^ 2 + 0 * (x - c i))
    (fun i => by fun_prop)]
  have hw : 0 < 1 - x := by linarith
  -- each factor is at most exp(κ(v + (mᵢ − cᵢ)²)/(1 − x))
  have hfac : ∀ i, ∫⁻ y, ENNReal.ofReal (rexp (κ * (y - c i) ^ 2 + 0 * (y - c i)))
      ∂gaussianReal (m i) v.toNNReal
      ≤ ENNReal.ofReal (rexp (κ * (v + (m i - c i) ^ 2) / (1 - x))) := by
    intro i
    rw [lintegral_exp_quadratic_gaussianReal (m i) (c i) κ 0 hv (by linarith)]
    apply ENNReal.ofReal_le_ofReal
    have h1 := one_div_sqrt_one_sub_le_exp (by linarith : x < 1)
    have e1 : 1 - 2 * κ * v = 1 - x := by rw [hx]
    have e2 : κ * (m i - c i) ^ 2 + 0 * (m i - c i) + 0 ^ 2 * v / 2 = κ * (m i - c i) ^ 2 := by
      ring
    rw [e1, e2]
    calc 1 / √(1 - x) * rexp (κ * (m i - c i) ^ 2 / (1 - x))
        ≤ rexp (x / 2 / (1 - x)) * rexp (κ * (m i - c i) ^ 2 / (1 - x)) :=
          mul_le_mul_of_nonneg_right h1 (Real.exp_pos _).le
      _ = rexp (κ * (v + (m i - c i) ^ 2) / (1 - x)) := by
          rw [← Real.exp_add, ← hx]
          congr 1
          ring
  have hsumeq : ∑ i, κ * (v + (m i - c i) ^ 2) / (1 - x) = κ * A / (1 - x) := by
    rw [← Finset.sum_div, ← Finset.mul_sum, Finset.sum_add_distrib, Finset.sum_const,
      Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, hM, hA]
  calc ENNReal.ofReal (rexp (-(κ * D)))
        * ∏ i, ∫⁻ y, ENNReal.ofReal (rexp (κ * (y - c i) ^ 2 + 0 * (y - c i)))
          ∂gaussianReal (m i) v.toNNReal
      ≤ ENNReal.ofReal (rexp (-(κ * D)))
        * ∏ i, ENNReal.ofReal (rexp (κ * (v + (m i - c i) ^ 2) / (1 - x))) :=
        mul_le_mul_right (Finset.prod_le_prod' fun i _ => hfac i) _
    _ = ENNReal.ofReal (rexp (-(κ * D) + κ * A / (1 - x))) := by
        rw [← ENNReal.ofReal_prod_of_nonneg (fun i _ => (Real.exp_pos _).le), ← Real.exp_sum,
          ← ENNReal.ofReal_mul (Real.exp_pos _).le, ← Real.exp_add, hsumeq]
    _ ≤ ENNReal.ofReal (rexp (-((D - A) ^ 2 / (16 * v * D)))) := by
        apply ENNReal.ofReal_le_ofReal
        apply Real.exp_le_exp.2
        -- 1/(1 − x) ≤ 1 + 2x on [0, 1/2]
        have hκA : 0 ≤ κ * A := by positivity
        have hinv : κ * A / (1 - x) ≤ κ * A * (1 + 2 * x) := by
          rw [div_le_iff₀ hw]
          have h12 : 0 ≤ x * (1 - 2 * x) := mul_nonneg hx0 (by linarith)
          nlinarith [mul_nonneg hκA h12]
        -- −κD + κA(1 + 2x) = −κ(D − A) + 2xκA, and that is −(D − A)²/(16vD) − (D − A)³/(16vD²)
        have hkey : -(κ * D) + κ * A * (1 + 2 * x)
            = -((D - A) ^ 2 / (16 * v * D)) - (D - A) ^ 3 / (16 * v * D ^ 2) := by
          rw [hx_eq, ← hκ]
          field_simp
          ring
        have hcube : 0 ≤ (D - A) ^ 3 / (16 * v * D ^ 2) := by positivity
        linarith

/-- **The residual scale's lower tail.** For independent `Yᵢ ~ N(mᵢ, v)` and `τ < v`,
`P(Σ(Yᵢ − cᵢ)² ≤ nτ) ≤ exp(−n(v − τ)²/(8v²))`: Chernoff for `−Σ(Yᵢ − cᵢ)²` at
`κ = (v − τ)/(4v²)`, dropping the means, which only shift the sum up. -/
theorem pi_gaussian_sq_lower_tail (m c : Fin n → ℝ) {v τ : ℝ} (hv : 0 < v) (hτ : τ < v) :
    Measure.pi (fun i => gaussianReal (m i) v.toNNReal) {y | ∑ i, (y i - c i) ^ 2 ≤ n * τ}
      ≤ ENNReal.ofReal (rexp (-(n * (v - τ) ^ 2 / (8 * v ^ 2)))) := by
  obtain ⟨κ, hκ⟩ : ∃ κ, (v - τ) / (4 * v ^ 2) = κ := ⟨_, rfl⟩
  have hκ0 : 0 < κ := hκ ▸ div_pos (by linarith) (by positivity)
  have hX : Measurable fun y : Fin n → ℝ => -∑ i, (y i - c i) ^ 2 :=
    (Finset.measurable_sum _ fun i _ =>
      ((measurable_pi_apply i).sub measurable_const).pow_const 2).neg
  have hset : {y : Fin n → ℝ | ∑ i, (y i - c i) ^ 2 ≤ n * τ}
      = {y | -(n * τ) ≤ -∑ i, (y i - c i) ^ 2} := by
    ext y
    simp only [Set.mem_setOf_eq, neg_le_neg_iff]
  rw [hset]
  refine (measure_ge_le_exp_mul_lintegral _ hX hκ0.le (-(n * τ))).trans ?_
  have hsum : ∀ y : Fin n → ℝ, κ * -∑ i, (y i - c i) ^ 2
      = ∑ i, (-κ * (y i - c i) ^ 2 + 0 * (y i - c i)) := by
    intro y
    rw [mul_neg, ← neg_mul, Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  simp_rw [hsum]
  rw [lintegral_pi_exp_sum m v.toNNReal (fun i x => -κ * (x - c i) ^ 2 + 0 * (x - c i))
    (fun i => by fun_prop)]
  obtain ⟨x, hx⟩ : ∃ x, 2 * κ * v = x := ⟨_, rfl⟩
  have hx0 : 0 ≤ x := hx ▸ by positivity
  have hfac : ∀ i, ∫⁻ y, ENNReal.ofReal (rexp (-κ * (y - c i) ^ 2 + 0 * (y - c i)))
      ∂gaussianReal (m i) v.toNNReal
      ≤ ENNReal.ofReal (rexp (-(x / 2 - x ^ 2 / 2))) := by
    intro i
    rw [lintegral_exp_quadratic_gaussianReal (m i) (c i) (-κ) 0 hv (by nlinarith)]
    apply ENNReal.ofReal_le_ofReal
    have e1 : 1 - 2 * -κ * v = 1 + x := by rw [← hx]; ring
    rw [e1]
    have h1 := one_div_sqrt_one_add_le_exp hx0
    have h2 : rexp ((-κ * (m i - c i) ^ 2 + 0 * (m i - c i) + 0 ^ 2 * v / 2) / (1 + x)) ≤ 1 := by
      rw [Real.exp_le_one_iff]
      apply div_nonpos_of_nonpos_of_nonneg _ (by linarith)
      nlinarith [sq_nonneg (m i - c i)]
    calc 1 / √(1 + x) * rexp ((-κ * (m i - c i) ^ 2 + 0 * (m i - c i) + 0 ^ 2 * v / 2) / (1 + x))
        ≤ 1 / √(1 + x) * 1 := mul_le_mul_of_nonneg_left h2 (by positivity)
      _ ≤ rexp (-(x / 2 - x ^ 2 / 2)) := by rw [mul_one]; exact h1
  calc ENNReal.ofReal (rexp (-(κ * -(n * τ))))
        * ∏ i, ∫⁻ y, ENNReal.ofReal (rexp (-κ * (y - c i) ^ 2 + 0 * (y - c i)))
          ∂gaussianReal (m i) v.toNNReal
      ≤ ENNReal.ofReal (rexp (-(κ * -(n * τ))))
        * ∏ _i : Fin n, ENNReal.ofReal (rexp (-(x / 2 - x ^ 2 / 2))) :=
        mul_le_mul_right (Finset.prod_le_prod' fun i _ => hfac i) _
    _ = ENNReal.ofReal (rexp (κ * (n * τ) - n * (x / 2 - x ^ 2 / 2))) := by
        rw [← ENNReal.ofReal_prod_of_nonneg (fun i _ => (Real.exp_pos _).le), ← Real.exp_sum,
          ← ENNReal.ofReal_mul (Real.exp_pos _).le, ← Real.exp_add, Finset.sum_const,
          Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
        congr 2
        ring
    _ = ENNReal.ofReal (rexp (-(n * (v - τ) ^ 2 / (8 * v ^ 2)))) := by
        congr 2
        rw [← hx, ← hκ]
        have hv0 : v ≠ 0 := hv.ne'
        field_simp
        ring

end BGP
