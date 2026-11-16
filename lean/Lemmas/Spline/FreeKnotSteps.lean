import Mathlib.Algebra.Order.Floor.Semiring
import Mathlib.Analysis.Real.Sqrt
import Mathlib.NumberTheory.Real.Irrational
import Mathlib.Order.Interval.Set.OrdConnected
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity

/-!
# Free knots against fixed knots, for piecewise constants (S9)

The order-1 case of Schumaker's Theorem 7.12: with free knots the sup error of piecewise
constants is set by how much a monotone function rises, not by where it rises. Fixed uniform knots
cannot follow a jump.

* `PiecewiseConst s a b n`: `s` is constant on each of `n` intervals that cover `[a, b]`. That is
  `n − 1` free knots, the pieces open or closed at either end.
* `exists_piecewiseConst_monotone` (free): a monotone `f` on `[a, b]` and `k ≥ 0` admit `s` with
  `k + 1` pieces and `sup |f − s| ≤ (f b − f a) / (2(k + 1))`. Put `s` at the midpoint of the band
  of height `(f b − f a)/(k + 1)` that holds `f x`; the bands' preimages are intervals because
  `f` is monotone.
* `exists_monotone_far_from_uniform` (fixed): one monotone `f` with `f b − f a = 1` keeps
  `sup |f − s| ≥ ½` against every `s` constant on the open pieces of the `k` uniform breakpoints,
  for every `k`. Its jump sits at `a + (b − a)(√2 − 1)`, which is never a uniform breakpoint
  because `(k + 1)(√2 − 1)` is irrational.
-/

open Set

namespace Spline

/-- `s` is piecewise constant on `[a, b]` with at most `n` pieces: `n` intervals inside `[a, b]`
cover it, and `s` is constant on each. -/
def PiecewiseConst (s : ℝ → ℝ) (a b : ℝ) (n : ℕ) : Prop :=
  ∃ I : Fin n → Set ℝ, (∀ j, (I j).OrdConnected) ∧ (∀ j, I j ⊆ Icc a b) ∧
    Icc a b ⊆ ⋃ j, I j ∧ ∀ j, ∀ x ∈ I j, ∀ y ∈ I j, s x = s y

/-- **S9, free knots.** A function monotone on `[a, b]` is within `(f b − f a)/(2(k + 1))` of a
piecewise constant with `k + 1` pieces. -/
theorem exists_piecewiseConst_monotone {f : ℝ → ℝ} {a b : ℝ} (hab : a ≤ b)
    (hf : MonotoneOn f (Icc a b)) (k : ℕ) :
    ∃ s : ℝ → ℝ, PiecewiseConst s a b (k + 1) ∧
      ∀ x ∈ Icc a b, |f x - s x| ≤ (f b - f a) / (2 * (k + 1)) := by
  have ha : a ∈ Icc a b := left_mem_Icc.2 hab
  have hb : b ∈ Icc a b := right_mem_Icc.2 hab
  have hk : (0 : ℝ) < k + 1 := by positivity
  set h := (f b - f a) / (k + 1) with hh
  have hV : 0 ≤ f b - f a := sub_nonneg.2 (hf ha hb hab)
  have h0 : 0 ≤ h := div_nonneg hV hk.le
  -- the band of `f x`, capped at the top band `k`
  set q : ℝ → ℕ := fun x => min ⌊(f x - f a) / h⌋₊ k with hq
  have hqmono : ∀ x ∈ Icc a b, ∀ y ∈ Icc a b, x ≤ y → q x ≤ q y := fun x hx y hy hxy =>
    min_le_min_right _ (Nat.floor_le_floor (div_le_div_of_nonneg_right
      (sub_le_sub_right (hf hx hy hxy) _) h0))
  have hqk : ∀ x, q x ≤ k := fun x => min_le_right _ _
  refine ⟨fun x => f a + (q x + 1 / 2) * h, ?_, ?_⟩
  · refine ⟨fun j => {x | x ∈ Icc a b ∧ q x = j}, fun j => ⟨?_⟩, fun j x hx => hx.1, ?_, ?_⟩
    · rintro x ⟨hx, hxj⟩ z ⟨hz, hzj⟩ y ⟨hxy, hyz⟩
      have hy : y ∈ Icc a b := ⟨hx.1.trans hxy, hyz.trans hz.2⟩
      exact ⟨hy, le_antisymm (hzj ▸ hqmono y hy z hz hyz) (hxj ▸ hqmono x hx y hy hxy)⟩
    · intro x hx
      exact mem_iUnion.2 ⟨⟨q x, Nat.lt_succ_of_le (hqk x)⟩, hx, rfl⟩
    · rintro j x ⟨-, hxj⟩ y ⟨-, hyj⟩
      simp only [hxj, hyj]
  · intro x hx
    have hfx : f a ≤ f x := hf ha hx hx.1
    have hfb : f x ≤ f b := hf hx hb hx.2
    rcases h0.eq_or_lt with hzero | hpos
    · -- no rise: `f` is constant on `[a, b]`
      have hc : f b - f a = 0 := by
        rw [hh] at hzero
        rcases div_eq_zero_iff.1 hzero.symm with h' | h'
        · exact h'
        · exact absurd h' hk.ne'
      rw [← hzero, hc]
      simp only [mul_zero, add_zero, zero_div]
      rw [abs_nonpos_iff, sub_eq_zero]
      linarith
    · set u := (f x - f a) / h with hu
      have hu0 : 0 ≤ u := div_nonneg (sub_nonneg.2 hfx) h0
      have hlo : (q x : ℝ) ≤ u :=
        (Nat.cast_le.2 (min_le_left _ _)).trans (Nat.floor_le hu0)
      have hhi : u ≤ q x + 1 := by
        rcases le_total ⌊u⌋₊ k with hle | hle
        · have : q x = ⌊u⌋₊ := min_eq_left hle
          rw [this]
          exact (Nat.lt_floor_add_one u).le
        · have : q x = k := min_eq_right hle
          have e : ((k : ℝ) + 1) * ((f b - f a) / (k + 1)) = f b - f a := by field_simp
          rw [this, hu, div_le_iff₀ hpos, hh, e]
          linarith
      have hfx' : f x - f a = u * h := by
        rw [hu, div_mul_cancel₀ _ hpos.ne']
      have hgoal : (f b - f a) / (2 * (k + 1)) = h / 2 := by
        rw [hh]
        field_simp
      rw [hgoal, abs_le]
      constructor <;> nlinarith [mul_le_mul_of_nonneg_right hlo hpos.le,
        mul_le_mul_of_nonneg_right hhi hpos.le]

/-- `s` is constant on each open piece `(t_j, t_{j+1})` of the `k` uniform breakpoints
`t_j = a + j·(b − a)/(k + 1)`; its values at the breakpoints are free. -/
def StepOnUniform (s : ℝ → ℝ) (a b : ℝ) (k : ℕ) : Prop :=
  ∀ j : ℕ, j ≤ k → ∀ x y,
    x ∈ Ioo (a + j * ((b - a) / (k + 1))) (a + (j + 1) * ((b - a) / (k + 1))) →
    y ∈ Ioo (a + j * ((b - a) / (k + 1))) (a + (j + 1) * ((b - a) / (k + 1))) → s x = s y

/-- **S9, fixed knots.** One monotone `f` with `f b − f a = 1` stays at sup distance at least `½`
from every piecewise constant on the `k` uniform breakpoints, for every `k`. -/
theorem exists_monotone_far_from_uniform {a b : ℝ} (hab : a < b) :
    ∃ f : ℝ → ℝ, Monotone f ∧ f b - f a = 1 ∧
      ∀ k : ℕ, ∀ s : ℝ → ℝ, StepOnUniform s a b k → ∃ x ∈ Icc a b, 1 / 2 ≤ |f x - s x| := by
  have hs1 : 1 < √2 := by
    rw [show (1 : ℝ) = √1 by simp]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  have hs2 : √2 < 2 := by
    rw [show (2 : ℝ) = √4 by rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  set r := √2 - 1 with hr
  have hr0 : 0 < r := by linarith
  have hr1 : r < 1 := by linarith
  set c := a + (b - a) * r with hc
  have hca : a < c := by nlinarith
  have hcb : c < b := by nlinarith
  refine ⟨fun x => if c ≤ x then 1 else 0, ?_, ?_, ?_⟩
  · intro x y hxy
    by_cases hx : c ≤ x
    · simp [hx, hx.trans hxy]
    · simp only [hx, if_false]
      split_ifs <;> norm_num
  · simp [hcb.le, not_le.2 hca]
  · intro k s hs
    have hk : (0 : ℝ) < k + 1 := by positivity
    set δ := (b - a) / (k + 1) with hδ
    have hδ0 : 0 < δ := div_pos (by linarith) hk
    set ρ := ((k : ℝ) + 1) * r with hρ
    have hρirr : Irrational ρ := by
      have := (irrational_sqrt_two.sub_natCast 1).natCast_mul (Nat.succ_ne_zero k)
      simpa [hρ, hr] using this
    set j := ⌊ρ⌋₊ with hj
    have hjρ : (j : ℝ) < ρ :=
      lt_of_le_of_ne (Nat.floor_le (by positivity)) (fun h => hρirr.ne_nat j h.symm)
    have hρj : ρ < j + 1 := Nat.lt_floor_add_one ρ
    have hjk : j ≤ k := by
      have : ρ < k + 1 := by rw [hρ]; nlinarith
      exact Nat.le_of_lt_succ (by exact_mod_cast hjρ.trans this)
    have hcρ : c = a + ρ * δ := by
      rw [hc, hρ, hδ]
      field_simp
    have hlo : a + j * δ < c := by rw [hcρ]; nlinarith
    have hhi : c < a + (j + 1) * δ := by rw [hcρ]; nlinarith
    set x₁ := (a + j * δ + c) / 2 with hx₁
    have hx₁c : x₁ < c := by rw [hx₁]; linarith
    have hx₁lo : a + j * δ < x₁ := by rw [hx₁]; linarith
    have hseq : s x₁ = s c := hs j hjk x₁ c ⟨hx₁lo, hx₁c.trans hhi⟩ ⟨hlo, hhi⟩
    have hjδ : 0 ≤ (j : ℝ) * δ := by positivity
    have hx₁a : a ≤ x₁ := by linarith
    have hf₁ : (if c ≤ x₁ then (1 : ℝ) else 0) = 0 := if_neg (not_le.2 hx₁c)
    have hfc : (if c ≤ c then (1 : ℝ) else 0) = 1 := if_pos le_rfl
    by_cases hv : 1 / 2 ≤ |(1 : ℝ) - s c|
    · exact ⟨c, ⟨hca.le, hcb.le⟩, by beta_reduce; rw [hfc]; exact hv⟩
    · refine ⟨x₁, ⟨hx₁a, (hx₁c.trans hcb).le⟩, ?_⟩
      beta_reduce
      rw [hf₁, hseq]
      rw [not_le, abs_lt] at hv
      rw [abs_sub_comm, sub_zero]
      rw [le_abs]
      left
      linarith [hv.1]

/-! ### Witnesses -/

/-- The free-knot statement's hypotheses hold for `f = id` on `[0, 1]`: two pieces reach `¼`. -/
example : ∃ s : ℝ → ℝ, PiecewiseConst s 0 1 (1 + 1) ∧
    ∀ x ∈ Icc (0 : ℝ) 1, |id x - s x| ≤ (id 1 - id 0) / (2 * ((1 : ℕ) + 1)) :=
  exists_piecewiseConst_monotone (by norm_num) (monotone_id.monotoneOn _) 1

/-- The fixed-knot statement on `[0, 1]`. -/
example : ∃ f : ℝ → ℝ, Monotone f ∧ f 1 - f 0 = 1 ∧
    ∀ k : ℕ, ∀ s : ℝ → ℝ, StepOnUniform s 0 1 k → ∃ x ∈ Icc (0 : ℝ) 1, 1 / 2 ≤ |f x - s x| :=
  exists_monotone_far_from_uniform (by norm_num)

end Spline
