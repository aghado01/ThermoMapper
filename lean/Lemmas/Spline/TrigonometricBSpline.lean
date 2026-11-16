import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity

/-!
# The order-2 trigonometric B-spline and the damped oscillator (S8)

After Schumaker § 10.3 and § 10.8. With carrier `ω > 0`, the pieces live in
`span{cos ωx, sin ωx}`, the null space of `D² + ω²`; the posed case is `ω = ½`.

* `trigBSpline ω y₀ y₁ y₂`: (10.103) applied to (10.99), the piece
  `sin(ω(x − y₀))/(sin(ω(y₁ − y₀))·sin(ω(y₂ − y₀)))` on `[y₀, y₁)` and
  `sin(ω(y₂ − x))/(sin(ω(y₂ − y₁))·sin(ω(y₂ − y₀)))` on `[y₁, y₂)`, zero elsewhere
  (`trigBSpline_eq_left`, `trigBSpline_eq_right`, `trigBSpline_eq_zero`).
* `sin_shift_mem`, `hasDerivAt_trig`, `trig_ode`: each piece is in the span, and the span solves
  `f″ + ω²f = 0`.
* `continuous_trigBSpline` and `trigBSpline_pos`: when `ω(y₂ − y₀) < π`, that is
  `y₂ − y₀ < π/ω`, half the carrier's period, the spline is continuous on `ℝ` and positive on
  `(y₀, y₂)`. For `ω = ½` the condition is `y₂ − y₀ < 2π`.
* Sharpness. `collocation_det` says the Tchebycheff determinant of the span at `s < t` is
  `sin(ω(t − s))`, zero at `t − s = π/ω`. There, `twoKnotBump` (`sin(ω(x − y₀))` on `[y₀, y₂]`,
  zero outside) is continuous, positive inside, zero outside, with no interior knot. That is the
  local-support element the corrected Lemma 10.16 excludes when its hypothesis holds.
* The damped oscillator.
  - `nullDamped γ ω A B = e^{−γx}(A cos ωx + B sin ωx)` spans the null space of
    `D² + 2γD + (γ² + ω²)` (`nullDamped_ode`).
  - The Green's function `green γ ω` is `e^{−γx} sin(ωx)/ω` for `x ≥ 0` and `0` for `x < 0`. It is
    continuous, `green 0 = 0`, it has right derivative `1` at `0`, and it solves the equation
    for `x > 0`.
  - Theorem 10.15 at order 2 with simple knots, `τ 0 < τ 1 < … < τ (k + 1)`: the null space with
    `green(· − τ i)` spans the splines whose pieces lie in the null space (`spline_span`), and
    they are linearly independent (`spline_indep`).
-/

open Set Real

namespace Spline

/-! ### The span of `cos ωx`, `sin ωx` -/

/-- `sin(ω(x − c))` is in `span{cos ωx, sin ωx}`. -/
theorem sin_shift_mem (ω c : ℝ) :
    ∃ A B : ℝ, ∀ x, sin (ω * (x - c)) = A * cos (ω * x) + B * sin (ω * x) :=
  ⟨-sin (ω * c), cos (ω * c), fun x => by rw [mul_sub, sin_sub]; ring⟩

/-- `sin(ω(c − x))` is in `span{cos ωx, sin ωx}`. -/
theorem sin_shift_mem' (ω c : ℝ) :
    ∃ A B : ℝ, ∀ x, sin (ω * (c - x)) = A * cos (ω * x) + B * sin (ω * x) :=
  ⟨sin (ω * c), -cos (ω * c), fun x => by rw [mul_sub, sin_sub]; ring⟩

theorem hasDerivAt_linear (ω x : ℝ) : HasDerivAt (fun x => ω * x) ω x :=
  ((hasDerivAt_id' x).const_mul ω).congr_deriv (mul_one ω)

theorem hasDerivAt_trig (ω A B x : ℝ) :
    HasDerivAt (fun x => A * cos (ω * x) + B * sin (ω * x))
      ((ω * B) * cos (ω * x) + (-(ω * A)) * sin (ω * x)) x := by
  have hc : HasDerivAt (fun x => cos (ω * x)) (-sin (ω * x) * ω) x :=
    (hasDerivAt_cos (ω * x)).comp x (hasDerivAt_linear ω x)
  have hs : HasDerivAt (fun x => sin (ω * x)) (cos (ω * x) * ω) x :=
    (hasDerivAt_sin (ω * x)).comp x (hasDerivAt_linear ω x)
  exact ((hc.const_mul A).add (hs.const_mul B)).congr_deriv (by ring)

/-- **The null space of `D² + ω²`.** Every `A cos ωx + B sin ωx` solves `f″ + ω²f = 0`. -/
theorem trig_ode (ω A B x : ℝ) :
    deriv (deriv fun x => A * cos (ω * x) + B * sin (ω * x)) x +
      ω ^ 2 * (A * cos (ω * x) + B * sin (ω * x)) = 0 := by
  have h1 : deriv (fun x => A * cos (ω * x) + B * sin (ω * x)) =
      fun x => (ω * B) * cos (ω * x) + (-(ω * A)) * sin (ω * x) :=
    deriv_eq fun x => hasDerivAt_trig ω A B x
  rw [h1, (hasDerivAt_trig ω (ω * B) (-(ω * A)) x).deriv]
  ring

/-! ### The order-2 trigonometric B-spline -/

/-- The order-2 trigonometric B-spline with carrier `ω` on knots `y₀, y₁, y₂`, written so each
branch is continuous; `trigBSpline_eq_left`, `_right` and `_zero` give the posed half-open
pieces. -/
noncomputable def trigBSpline (ω y₀ y₁ y₂ : ℝ) (x : ℝ) : ℝ :=
  if x ≤ y₁ then
    (if x ≤ y₀ then 0 else sin (ω * (x - y₀)) / (sin (ω * (y₁ - y₀)) * sin (ω * (y₂ - y₀))))
  else
    (if x ≤ y₂ then sin (ω * (y₂ - x)) / (sin (ω * (y₂ - y₁)) * sin (ω * (y₂ - y₀))) else 0)

section Knots

variable {ω y₀ y₁ y₂ : ℝ}

theorem trigBSpline_eq_left {x : ℝ} (hx : x ∈ Ico y₀ y₁) :
    trigBSpline ω y₀ y₁ y₂ x =
      sin (ω * (x - y₀)) / (sin (ω * (y₁ - y₀)) * sin (ω * (y₂ - y₀))) := by
  unfold trigBSpline
  rw [if_pos hx.2.le]
  rcases hx.1.eq_or_lt with h | h
  · subst h
    simp
  · rw [if_neg (not_le.2 h)]

theorem trigBSpline_eq_right (h01 : y₀ < y₁)
    (hs1 : sin (ω * (y₁ - y₀)) ≠ 0) (hs2 : sin (ω * (y₂ - y₁)) ≠ 0) {x : ℝ} (hx : x ∈ Ico y₁ y₂) :
    trigBSpline ω y₀ y₁ y₂ x =
      sin (ω * (y₂ - x)) / (sin (ω * (y₂ - y₁)) * sin (ω * (y₂ - y₀))) := by
  unfold trigBSpline
  rcases hx.1.eq_or_lt with h | h
  · subst h
    rw [if_pos le_rfl, if_neg (not_le.2 h01)]
    rw [div_mul_cancel_left₀ hs1, div_mul_cancel_left₀ hs2]
  · rw [if_neg (not_le.2 h), if_pos hx.2.le]

theorem trigBSpline_eq_zero (h01 : y₀ < y₁) (h12 : y₁ < y₂) {x : ℝ} (hx : x ∉ Ico y₀ y₂) :
    trigBSpline ω y₀ y₁ y₂ x = 0 := by
  unfold trigBSpline
  rcases not_and_or.1 hx with h | h
  · rw [not_le] at h
    rw [if_pos (h.le.trans h01.le), if_pos h.le]
  · rw [not_lt] at h
    rcases h.eq_or_lt with h' | h'
    · subst h'
      rw [if_neg (not_le.2 h12), if_pos le_rfl, sub_self, mul_zero, sin_zero, zero_div]
    · rw [if_neg (not_le.2 (h12.trans h')), if_neg (not_le.2 h')]

/-- **S8, continuity.** With the sines at the knot gaps nonzero, which `ω(y₂ − y₀) < π` gives,
the spline is continuous on `ℝ`. -/
theorem continuous_trigBSpline (h01 : y₀ < y₁) (h12 : y₁ < y₂)
    (hs1 : sin (ω * (y₁ - y₀)) ≠ 0) (hs2 : sin (ω * (y₂ - y₁)) ≠ 0) :
    Continuous (trigBSpline ω y₀ y₁ y₂) := by
  unfold trigBSpline
  refine Continuous.if_le ?_ ?_ continuous_id continuous_const ?_
  · refine Continuous.if_le continuous_const (by fun_prop) continuous_id continuous_const ?_
    intro x hx
    rw [hx, sub_self, mul_zero, sin_zero, zero_div]
  · refine Continuous.if_le (by fun_prop) continuous_const continuous_id continuous_const ?_
    intro x hx
    rw [hx, sub_self, mul_zero, sin_zero, zero_div]
  · intro x hx
    subst hx
    rw [if_neg (not_le.2 h01), if_pos h12.le, div_mul_cancel_left₀ hs1,
      div_mul_cancel_left₀ hs2]

theorem sin_pos_of_lt (hω : 0 < ω) {z : ℝ} (hz : 0 < z) (hzπ : ω * z < π) : 0 < sin (ω * z) :=
  sin_pos_of_pos_of_lt_pi (mul_pos hω hz) hzπ

/-- **S8, positivity.** When `ω(y₂ − y₀) < π`, the spline is positive on `(y₀, y₂)`. -/
theorem trigBSpline_pos (hω : 0 < ω) (h01 : y₀ < y₁) (h12 : y₁ < y₂) (hπ : ω * (y₂ - y₀) < π)
    {x : ℝ} (hx : x ∈ Ioo y₀ y₂) : 0 < trigBSpline ω y₀ y₁ y₂ x := by
  have hlt : ∀ z, z ≤ y₂ - y₀ → ω * z < π := fun z hz =>
    lt_of_le_of_lt (mul_le_mul_of_nonneg_left hz hω.le) hπ
  have hD₀ : 0 < sin (ω * (y₂ - y₀)) := sin_pos_of_lt hω (by linarith) hπ
  unfold trigBSpline
  by_cases hx1 : x ≤ y₁
  · rw [if_pos hx1, if_neg (not_le.2 hx.1)]
    exact div_pos (sin_pos_of_lt hω (by linarith [hx.1]) (hlt _ (by linarith [hx.2])))
      (mul_pos (sin_pos_of_lt hω (by linarith) (hlt _ (by linarith))) hD₀)
  · rw [if_neg hx1, if_pos hx.2.le]
    exact div_pos (sin_pos_of_lt hω (by linarith [hx.2]) (hlt _ (by linarith [hx.1])))
      (mul_pos (sin_pos_of_lt hω (by linarith) (hlt _ (by linarith))) hD₀)

/-- **S8 as posed** (`ω = ½`): continuous on `ℝ` and positive on `(y₀, y₂)` when
`y₂ − y₀ < 2π`. -/
theorem trigBSpline_half (h01 : y₀ < y₁) (h12 : y₁ < y₂) (h2π : y₂ - y₀ < 2 * π) :
    Continuous (trigBSpline (1 / 2) y₀ y₁ y₂) ∧
      ∀ x ∈ Ioo y₀ y₂, 0 < trigBSpline (1 / 2) y₀ y₁ y₂ x := by
  have hω : (0 : ℝ) < 1 / 2 := by norm_num
  have hπ : 1 / 2 * (y₂ - y₀) < π := by linarith
  have hgap : ∀ z, 0 < z → z < y₂ - y₀ → sin (1 / 2 * z) ≠ 0 := fun z hz hz' =>
    (sin_pos_of_lt hω hz (by linarith)).ne'
  exact ⟨continuous_trigBSpline h01 h12 (hgap _ (by linarith) (by linarith))
      (hgap _ (by linarith) (by linarith)),
    fun x hx => trigBSpline_pos hω h01 h12 hπ hx⟩

end Knots

/-! ### Sharpness at half the period -/

/-- The Tchebycheff determinant of `span{cos ωx, sin ωx}` at `s, t` is `sin(ω(t − s))`: nonzero
for `0 < ω(t − s) < π`, zero at `ω(t − s) = π`. -/
theorem collocation_det (ω s t : ℝ) :
    cos (ω * s) * sin (ω * t) - sin (ω * s) * cos (ω * t) = sin (ω * (t - s)) := by
  rw [mul_sub, sin_sub]
  ring

/-- `sin(ω(x − y₀))` on `[y₀, y₂]`, zero outside. -/
noncomputable def twoKnotBump (ω y₀ y₂ : ℝ) (x : ℝ) : ℝ :=
  if x ≤ y₀ then 0 else (if x ≤ y₂ then sin (ω * (x - y₀)) else 0)

/-- **S8, sharpness.** At `ω(y₂ − y₀) = π`, the bump is continuous on `ℝ`, positive on
`(y₀, y₂)`, zero off `(y₀, y₂)`, and on `[y₀, y₂]` it is one function of the span: a
local-support element on two knots. -/
theorem twoKnotBump_spec {ω y₀ y₂ : ℝ} (hω : 0 < ω) (hπ : ω * (y₂ - y₀) = π) :
    Continuous (twoKnotBump ω y₀ y₂) ∧ (∀ x ∈ Ioo y₀ y₂, 0 < twoKnotBump ω y₀ y₂ x) ∧
      (∀ x ∉ Ioo y₀ y₂, twoKnotBump ω y₀ y₂ x = 0) ∧
      ∃ A B : ℝ, ∀ x ∈ Icc y₀ y₂, twoKnotBump ω y₀ y₂ x = A * cos (ω * x) + B * sin (ω * x) := by
  have h02 : y₀ < y₂ := by
    by_contra h
    have : ω * (y₂ - y₀) ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hω.le (by linarith)
    linarith [pi_pos]
  refine ⟨?_, fun x hx => ?_, fun x hx => ?_, ?_⟩
  · unfold twoKnotBump
    refine Continuous.if_le continuous_const ?_ continuous_id continuous_const ?_
    · refine Continuous.if_le (by fun_prop) continuous_const continuous_id continuous_const ?_
      intro x hx
      rw [hx, hπ, sin_pi]
    · intro x hx
      rw [hx, if_pos h02.le, sub_self, mul_zero, sin_zero]
  · unfold twoKnotBump
    rw [if_neg (not_le.2 hx.1), if_pos hx.2.le]
    refine sin_pos_of_lt hω (by linarith [hx.1]) ?_
    rw [← hπ]
    exact mul_lt_mul_of_pos_left (by linarith [hx.2]) hω
  · unfold twoKnotBump
    rcases not_and_or.1 hx with h | h
    · rw [if_pos (not_lt.1 h)]
    · rw [not_lt] at h
      rcases h.eq_or_lt with h' | h'
      · subst h'
        rw [if_neg (not_le.2 h02), if_pos le_rfl, hπ, sin_pi]
      · rw [if_neg (not_le.2 (h02.trans h')), if_neg (not_le.2 h')]
  · obtain ⟨A, B, hAB⟩ := sin_shift_mem ω y₀
    refine ⟨A, B, fun x hx => ?_⟩
    unfold twoKnotBump
    rcases hx.1.eq_or_lt with h | h
    · subst h
      rw [if_pos le_rfl, ← hAB, sub_self, mul_zero, sin_zero]
    · rw [if_neg (not_le.2 h), if_pos hx.2, hAB]

/-! ### The damped oscillator -/

/-- The null space of `D² + 2γD + (γ² + ω²)`: `e^{−γx}(A cos ωx + B sin ωx)`. -/
noncomputable def nullDamped (γ ω A B : ℝ) (x : ℝ) : ℝ :=
  exp (-γ * x) * (A * cos (ω * x) + B * sin (ω * x))

/-- Differentiation maps the null space to itself:
`(A, B) ↦ (−γA + ωB, −ωA − γB)`. -/
theorem hasDerivAt_nullDamped (γ ω A B x : ℝ) :
    HasDerivAt (nullDamped γ ω A B) (nullDamped γ ω (-γ * A + ω * B) (-ω * A - γ * B) x) x := by
  have he : HasDerivAt (fun x => exp (-γ * x)) (exp (-γ * x) * (-γ)) x :=
    (hasDerivAt_exp (-γ * x)).comp x (hasDerivAt_linear (-γ) x)
  have hc : HasDerivAt (fun x => cos (ω * x)) (-sin (ω * x) * ω) x :=
    (hasDerivAt_cos (ω * x)).comp x (hasDerivAt_linear ω x)
  have hs : HasDerivAt (fun x => sin (ω * x)) (cos (ω * x) * ω) x :=
    (hasDerivAt_sin (ω * x)).comp x (hasDerivAt_linear ω x)
  exact (he.mul ((hc.const_mul A).add (hs.const_mul B))).congr_deriv
    (by unfold nullDamped; simp only [Pi.add_apply]; ring)

theorem deriv_nullDamped (γ ω A B : ℝ) :
    deriv (nullDamped γ ω A B) = nullDamped γ ω (-γ * A + ω * B) (-ω * A - γ * B) :=
  deriv_eq fun x => hasDerivAt_nullDamped γ ω A B x

/-- **The null space.** Every `e^{−γx}(A cos ωx + B sin ωx)` solves
`f″ + 2γf′ + (γ² + ω²)f = 0`. -/
theorem nullDamped_ode (γ ω A B x : ℝ) :
    deriv (deriv (nullDamped γ ω A B)) x + 2 * γ * deriv (nullDamped γ ω A B) x +
      (γ ^ 2 + ω ^ 2) * nullDamped γ ω A B x = 0 := by
  rw [deriv_nullDamped, deriv_nullDamped]
  unfold nullDamped
  ring

/-- The damped oscillator's Green's function: `e^{−γx} sin(ωx)/ω` for `x ≥ 0`, `0` for `x < 0`. -/
noncomputable def green (γ ω : ℝ) (x : ℝ) : ℝ :=
  if 0 ≤ x then nullDamped γ ω 0 ω⁻¹ x else 0

theorem green_of_nonneg {γ ω x : ℝ} (hx : 0 ≤ x) :
    green γ ω x = exp (-γ * x) * sin (ω * x) / ω := by
  rw [green, if_pos hx, nullDamped]
  ring

theorem green_zero (γ ω : ℝ) : green γ ω 0 = 0 := by
  rw [green_of_nonneg le_rfl]
  simp

theorem green_eq_zero_of_nonpos {γ ω x : ℝ} (hx : x ≤ 0) : green γ ω x = 0 := by
  rcases hx.eq_or_lt with h | h
  · rw [h, green_zero]
  · rw [green, if_neg (not_le.2 h)]

/-- **S8, the Green's function.** It is continuous, vanishes at `0`, and has right derivative
`1` there. -/
theorem green_spec (γ : ℝ) {ω : ℝ} (hω : ω ≠ 0) :
    Continuous (green γ ω) ∧ green γ ω 0 = 0 ∧ HasDerivWithinAt (green γ ω) 1 (Ici 0) 0 := by
  refine ⟨?_, green_zero γ ω, ?_⟩
  · unfold green
    refine Continuous.if_le ?_ continuous_const continuous_const continuous_id ?_
    · unfold nullDamped
      fun_prop
    · intro x hx
      rw [← hx]
      unfold nullDamped
      simp
  · have h := (hasDerivAt_nullDamped γ ω 0 ω⁻¹ 0).hasDerivWithinAt (s := Ici 0)
    have h1 : nullDamped γ ω (-γ * 0 + ω * ω⁻¹) (-ω * 0 - γ * ω⁻¹) 0 = 1 := by
      unfold nullDamped
      simp [hω]
    rw [h1] at h
    exact h.congr (fun x hx => if_pos hx) (if_pos le_rfl)

/-- **S8, the equation.** For `x > 0`, `g″ + 2γg′ + (γ² + ω²)g = 0`. -/
theorem green_ode (γ ω : ℝ) {x : ℝ} (hx : 0 < x) :
    deriv (deriv (green γ ω)) x + 2 * γ * deriv (green γ ω) x +
      (γ ^ 2 + ω ^ 2) * green γ ω x = 0 := by
  have hev : green γ ω =ᶠ[nhds x] nullDamped γ ω 0 ω⁻¹ := by
    filter_upwards [Ioi_mem_nhds hx] with y hy using if_pos (le_of_lt hy)
  rw [hev.deriv.deriv_eq, hev.deriv_eq, hev.eq_of_nhds]
  exact nullDamped_ode γ ω 0 ω⁻¹ x

/-! ### Theorem 10.15 at order 2 with simple knots -/

/-- On `[τ, ∞)` the Green's function shifted to `τ` is in the null space. -/
theorem green_shift_mem (γ ω τ : ℝ) :
    ∃ P Q : ℝ, ∀ t, τ ≤ t → green γ ω (t - τ) = nullDamped γ ω P Q t := by
  refine ⟨-(exp (γ * τ) * sin (ω * τ) / ω), exp (γ * τ) * cos (ω * τ) / ω, fun t ht => ?_⟩
  rw [green_of_nonneg (sub_nonneg.2 ht), nullDamped,
    show -γ * (t - τ) = -γ * t + γ * τ by ring, exp_add, mul_sub, sin_sub]
  ring

/-- A null-space function vanishing at `τ` is a multiple of `green(· − τ)` on `[τ, ∞)`. -/
theorem nullDamped_eq_green {γ ω A B τ : ℝ} (hω : ω ≠ 0) (h0 : nullDamped γ ω A B τ = 0) :
    ∃ c : ℝ, ∀ t, τ ≤ t → nullDamped γ ω A B t = c * green γ ω (t - τ) := by
  have hz : A * cos (ω * τ) + B * sin (ω * τ) = 0 := by
    unfold nullDamped at h0
    rcases mul_eq_zero.1 h0 with h | h
    · exact absurd h (exp_pos _).ne'
    · exact h
  refine ⟨ω * (exp (-γ * τ) * (B * cos (ω * τ) - A * sin (ω * τ))), fun t ht => ?_⟩
  have key : ∀ W X : ℝ, ω * W * (X / ω) = W * X := fun W X => by field_simp
  rw [green_of_nonneg (sub_nonneg.2 ht), key, nullDamped,
    show ω * t = ω * (t - τ) + ω * τ by ring, cos_add, sin_add,
    show -γ * t = -γ * (t - τ) + -γ * τ by ring, exp_add]
  linear_combination (exp (-γ * (t - τ)) * exp (-γ * τ) * cos (ω * (t - τ))) * hz

theorem nullDamped_add (γ ω A B C D t : ℝ) :
    nullDamped γ ω A B t + nullDamped γ ω C D t = nullDamped γ ω (A + C) (B + D) t := by
  unfold nullDamped
  ring

theorem nullDamped_zero (γ ω t : ℝ) : nullDamped γ ω 0 0 t = 0 := by
  unfold nullDamped
  ring

/-- **S8, Theorem 10.15's span at order 2.** On knots `τ 0 < … < τ (k + 1)`, a function whose
pieces on `[τ i, τ (i + 1)]` lie in the null space is a null-space function plus
`Σᵢ dᵢ·green(· − τ i)` over the interior knots. -/
theorem spline_span {γ ω : ℝ} (hω : ω ≠ 0) (τ : ℕ → ℝ) (hτ : StrictMono τ) (s : ℝ → ℝ) :
    ∀ k : ℕ, (∀ i ≤ k, ∃ A B : ℝ, ∀ t ∈ Icc (τ i) (τ (i + 1)), s t = nullDamped γ ω A B t) →
      ∃ (A B : ℝ) (d : ℕ → ℝ), ∀ t ∈ Icc (τ 0) (τ (k + 1)),
        s t = nullDamped γ ω A B t + ∑ i ∈ Finset.Icc 1 k, d i * green γ ω (t - τ i) := by
  intro k
  induction k with
  | zero =>
    intro hs
    obtain ⟨A, B, hAB⟩ := hs 0 le_rfl
    exact ⟨A, B, 0, fun t ht => by simpa using hAB t ht⟩
  | succ k ih =>
    intro hs
    obtain ⟨A, B, d, hrep⟩ := ih fun i hi => hs i (by omega)
    obtain ⟨A', B', hpiece⟩ := hs (k + 1) le_rfl
    -- on `[τ (k + 1), ∞)` the current representation is one null-space function
    have hshift : ∀ i, ∃ P Q : ℝ, ∀ t, τ i ≤ t → green γ ω (t - τ i) = nullDamped γ ω P Q t :=
      fun i => green_shift_mem γ ω (τ i)
    choose P Q hPQ using hshift
    set A'' := A + ∑ i ∈ Finset.Icc 1 k, d i * P i
    set B'' := B + ∑ i ∈ Finset.Icc 1 k, d i * Q i
    have hcur : ∀ t, τ (k + 1) ≤ t →
        nullDamped γ ω A B t + ∑ i ∈ Finset.Icc 1 k, d i * green γ ω (t - τ i) =
          nullDamped γ ω A'' B'' t := by
      intro t ht
      have hsum : ∑ i ∈ Finset.Icc 1 k, d i * green γ ω (t - τ i) =
          nullDamped γ ω (∑ i ∈ Finset.Icc 1 k, d i * P i) (∑ i ∈ Finset.Icc 1 k, d i * Q i) t := by
        rw [Finset.sum_congr rfl fun i hi => by
          rw [hPQ i t ((hτ.monotone (by simp at hi; omega)).trans ht)]]
        unfold nullDamped
        rw [Finset.sum_mul, Finset.sum_mul, ← Finset.sum_add_distrib, Finset.mul_sum]
        refine Finset.sum_congr rfl fun i _ => ?_
        ring
      rw [hsum, nullDamped_add]
    -- the jump at `τ (k + 1)` is a multiple of the shifted Green's function
    have hk1 : τ (k + 1) ∈ Icc (τ 0) (τ (k + 1)) := ⟨hτ.monotone (Nat.zero_le _), le_rfl⟩
    have hk1' : τ (k + 1) ∈ Icc (τ (k + 1)) (τ (k + 1 + 1)) :=
      ⟨le_rfl, (hτ (Nat.lt_succ_self _)).le⟩
    have hvanish : nullDamped γ ω (A' + -A'') (B' + -B'') (τ (k + 1)) = 0 := by
      have e1 := hrep _ hk1
      have e2 := hpiece _ hk1'
      rw [hcur _ le_rfl] at e1
      rw [← nullDamped_add, ← e2, e1]
      unfold nullDamped
      ring
    obtain ⟨c, hc⟩ := nullDamped_eq_green hω hvanish
    refine ⟨A, B, Function.update d (k + 1) c, fun t ht => ?_⟩
    rw [Finset.sum_Icc_succ_top (by omega), Function.update_self]
    have hupd : ∑ i ∈ Finset.Icc 1 k, Function.update d (k + 1) c i * green γ ω (t - τ i) =
        ∑ i ∈ Finset.Icc 1 k, d i * green γ ω (t - τ i) :=
      Finset.sum_congr rfl fun i hi => by
        rw [Function.update_of_ne (by simp at hi; omega)]
    rw [hupd]
    by_cases htk : t ≤ τ (k + 1)
    · rw [green_eq_zero_of_nonpos (by linarith), mul_zero, add_zero]
      exact hrep t ⟨ht.1, htk⟩
    · have htk' : τ (k + 1) ≤ t := (not_le.1 htk).le
      rw [← add_assoc, hcur t htk', ← hc t htk', nullDamped_add, hpiece t ⟨htk', ht.2⟩]
      congr 1 <;> ring

/-- **S8, Theorem 10.15's independence at order 2.** On knots `τ 0 < … < τ (k + 1)`, a null-space
function plus `Σᵢ dᵢ·green(· − τ i)` that vanishes on `[τ 0, τ (k + 1)]` has every coefficient
zero. -/
theorem spline_indep {γ ω : ℝ} (hω : 0 < ω) (τ : ℕ → ℝ) (hτ : StrictMono τ) (k : ℕ)
    (A B : ℝ) (d : ℕ → ℝ)
    (h : ∀ t ∈ Icc (τ 0) (τ (k + 1)),
      nullDamped γ ω A B t + ∑ i ∈ Finset.Icc 1 k, d i * green γ ω (t - τ i) = 0) :
    A = 0 ∧ B = 0 ∧ ∀ i ∈ Finset.Icc 1 k, d i = 0 := by
  -- a step inside each piece on which `sin` stays positive
  have hstep : ∀ i, ∃ δ, 0 < δ ∧ τ i + δ ≤ τ (i + 1) ∧ ω * δ < π := fun i => by
    refine ⟨min (τ (i + 1) - τ i) (π / (2 * ω)), lt_min (by linarith [hτ (Nat.lt_succ_self i)])
      (by positivity), by linarith [min_le_left (τ (i + 1) - τ i) (π / (2 * ω))], ?_⟩
    calc ω * min (τ (i + 1) - τ i) (π / (2 * ω)) ≤ ω * (π / (2 * ω)) :=
          mul_le_mul_of_nonneg_left (min_le_right _ _) hω.le
      _ = π / 2 := by field_simp
      _ < π := by linarith [pi_pos]
  have hτmono := hτ.monotone
  -- off its knot's right, every Green's term vanishes
  have hzero_after : ∀ t i, t ≤ τ i → green γ ω (t - τ i) = 0 := fun t i ht =>
    green_eq_zero_of_nonpos (by linarith)
  -- on the first piece only the null-space part remains
  have hfirst : ∀ t ∈ Icc (τ 0) (τ 1), nullDamped γ ω A B t = 0 := by
    intro t ht
    have := h t ⟨ht.1, ht.2.trans (hτmono (by omega))⟩
    rwa [Finset.sum_eq_zero fun i hi => by
      rw [hzero_after t i (ht.2.trans (hτmono (by simp at hi; omega))), mul_zero], add_zero] at this
  obtain ⟨δ, hδ, hδle, hδπ⟩ := hstep 0
  have e1 := hfirst (τ 0) ⟨le_rfl, (hτ (by omega)).le⟩
  have e2 := hfirst (τ 0 + δ) ⟨by linarith, hδle⟩
  unfold nullDamped at e1 e2
  have e1' := (mul_eq_zero.1 e1).resolve_left (exp_pos _).ne'
  have e2' := (mul_eq_zero.1 e2).resolve_left (exp_pos _).ne'
  have hsinδ : 0 < sin (ω * δ) := sin_pos_of_pos_of_lt_pi (mul_pos hω hδ) hδπ
  have hdet := collocation_det ω (τ 0) (τ 0 + δ)
  rw [show τ 0 + δ - τ 0 = δ by ring] at hdet
  have hA : A * sin (ω * δ) = 0 := by
    rw [← hdet]
    linear_combination sin (ω * (τ 0 + δ)) * e1' - sin (ω * τ 0) * e2'
  have hB : B * sin (ω * δ) = 0 := by
    rw [← hdet]
    linear_combination -(cos (ω * (τ 0 + δ)) * e1') + cos (ω * τ 0) * e2'
  have hA0 : A = 0 := (mul_eq_zero.1 hA).resolve_right hsinδ.ne'
  have hB0 : B = 0 := (mul_eq_zero.1 hB).resolve_right hsinδ.ne'
  refine ⟨hA0, hB0, ?_⟩
  -- the knots' coefficients, left to right
  have hall : ∀ m, ∀ i ∈ Finset.Icc 1 k, i ≤ m → d i = 0 := by
    intro m
    induction m with
    | zero => intro i hi him; simp at hi; omega
    | succ m ihm =>
      intro i hi him
      rcases Nat.lt_or_ge i (m + 1) with hlt | hge
      · exact ihm i hi (by omega)
      have hieq : i = m + 1 := by omega
      subst hieq
      obtain ⟨δ', hδ', hδ'le, hδ'π⟩ := hstep (m + 1)
      have hmk : m + 1 ≤ k := (Finset.mem_Icc.1 hi).2
      set t := τ (m + 1) + δ'
      have ht : t ∈ Icc (τ 0) (τ (k + 1)) :=
        ⟨by linarith [hτmono (Nat.zero_le (m + 1))], hδ'le.trans (hτmono (by omega))⟩
      have hval := h t ht
      rw [hA0, hB0, nullDamped_zero, zero_add,
        Finset.sum_eq_single (m + 1) (fun j hj hjne => ?_) (fun hno => absurd hi hno)] at hval
      · have hg : 0 < green γ ω (t - τ (m + 1)) := by
          rw [show t - τ (m + 1) = δ' by simp [t], green_of_nonneg hδ'.le]
          exact div_pos (mul_pos (exp_pos _) (sin_pos_of_pos_of_lt_pi (mul_pos hω hδ') hδ'π)) hω
        exact (mul_eq_zero.1 hval).resolve_right hg.ne'
      · rcases Nat.lt_or_ge j (m + 1) with hjlt | hjge
        · rw [ihm j hj (by omega), zero_mul]
        · rw [hzero_after t j (hδ'le.trans (hτmono (by omega))), mul_zero]
  exact fun i hi => hall i i hi le_rfl

/-! ### Witnesses -/

/-- The posed spline on knots `0, 1, 2` (`2 < 2π`) is continuous and positive inside. -/
example : Continuous (trigBSpline (1 / 2) 0 1 2) ∧
    ∀ x ∈ Ioo (0 : ℝ) 2, 0 < trigBSpline (1 / 2) 0 1 2 x :=
  trigBSpline_half (by norm_num) (by norm_num) (by linarith [two_le_pi])

/-- The two-knot bump at `ω = ½` on `[0, 2π]`. -/
example : Continuous (twoKnotBump (1 / 2) 0 (2 * π)) :=
  (twoKnotBump_spec (by norm_num) (by ring)).1

/-- The Green's function at `γ = 1`, `ω = 2`. -/
example : Continuous (green 1 2) ∧ green 1 2 0 = 0 ∧ HasDerivWithinAt (green 1 2) 1 (Ici 0) 0 :=
  green_spec 1 (by norm_num)

/-- Independence on knots `0, 1, 2` at `γ = 1`, `ω = 2`. -/
example (A B : ℝ) (d : ℕ → ℝ)
    (h : ∀ t ∈ Icc ((fun n : ℕ => (n : ℝ)) 0) ((fun n : ℕ => (n : ℝ)) 2),
      nullDamped 1 2 A B t +
        ∑ i ∈ Finset.Icc 1 1, d i * green 1 2 (t - (fun n : ℕ => (n : ℝ)) i) = 0) :
    A = 0 ∧ B = 0 ∧ ∀ i ∈ Finset.Icc 1 1, d i = 0 :=
  spline_indep (by norm_num) (fun n : ℕ => (n : ℝ)) Nat.strictMono_cast 1 A B d h

end Spline
