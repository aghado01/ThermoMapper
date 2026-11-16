import Mathlib.Analysis.Normed.Module.Basic
import Mathlib.Tactic.Abel
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# Contraction from a regularity inequality (S6)

The deterministic half of the dense phase-retrieval argument (arXiv 2011.00288, Lemma 1 and
Theorem 2) in a single basin. In a real normed space, a step direction `v` that is within half
the error plus `2η` of the error `x − x⋆`, on the ball of radius `r` about `x⋆`, drives the
iteration `x_{t+1} = x_t − α·v(x_t)` to within `4η` at the rate `1 − α/2`:

* `norm_descent_step_le`: one step contracts, `‖x' − x⋆‖ ≤ (1 − α/2)·‖x − x⋆‖ + 2αη`;
* `norm_descent_le`: from `‖x₀ − x⋆‖ ≤ r`, every iterate stays in the ball and
  `‖x_t − x⋆‖ ≤ (1 − α/2)ᵗ·‖x₀ − x⋆‖ + 4η·(1 − (1 − α/2)ᵗ)`;
* `norm_descent_le_posed`: S6 as posed, `‖x_t − x⋆‖ ≤ (1 − α/2)ᵗ·‖x₀ − x⋆‖ + 4η`.

The step needs `0 ≤ α ≤ 1`; staying in the ball needs `4η ≤ r`.
-/

namespace TimeFrequency

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- The iteration `x_{t+1} = x_t − α·v(x_t)` from `x₀`. -/
def descent (v : E → E) (α : ℝ) (x₀ : E) : ℕ → E
  | 0 => x₀
  | t + 1 => descent v α x₀ t - α • v (descent v α x₀ t)

/-- **One step.** On the ball, `‖x − α·v(x) − x⋆‖ ≤ (1 − α/2)·‖x − x⋆‖ + 2αη`. -/
theorem norm_descent_step_le {v : E → E} {xs : E} {η r α : ℝ} (hα0 : 0 ≤ α) (hα1 : α ≤ 1)
    (hreg : ∀ x, ‖x - xs‖ ≤ r → ‖v x - (x - xs)‖ ≤ ‖x - xs‖ / 2 + 2 * η) {x : E}
    (hx : ‖x - xs‖ ≤ r) :
    ‖x - α • v x - xs‖ ≤ (1 - α / 2) * ‖x - xs‖ + 2 * α * η := by
  have e : x - α • v x - xs = (1 - α) • (x - xs) - α • (v x - (x - xs)) := by
    simp only [sub_smul, one_smul, smul_sub]
    abel
  have hreg' := mul_le_mul_of_nonneg_left (hreg x hx) hα0
  rw [e]
  calc ‖(1 - α) • (x - xs) - α • (v x - (x - xs))‖
      ≤ ‖(1 - α) • (x - xs)‖ + ‖α • (v x - (x - xs))‖ := norm_sub_le _ _
    _ = (1 - α) * ‖x - xs‖ + α * ‖v x - (x - xs)‖ := by
      rw [norm_smul, norm_smul, Real.norm_of_nonneg (by linarith), Real.norm_of_nonneg hα0]
    _ ≤ (1 - α / 2) * ‖x - xs‖ + 2 * α * η := by linarith

/-- **S6, sharpened.** From `‖x₀ − x⋆‖ ≤ r`, every iterate stays in the ball, and
`‖x_t − x⋆‖ ≤ (1 − α/2)ᵗ·‖x₀ − x⋆‖ + 4η·(1 − (1 − α/2)ᵗ)`. -/
theorem norm_descent_le {v : E → E} {xs x₀ : E} {η r α : ℝ} (hηr : 4 * η ≤ r) (hα0 : 0 ≤ α)
    (hα1 : α ≤ 1) (hreg : ∀ x, ‖x - xs‖ ≤ r → ‖v x - (x - xs)‖ ≤ ‖x - xs‖ / 2 + 2 * η)
    (h₀ : ‖x₀ - xs‖ ≤ r) (t : ℕ) :
    ‖descent v α x₀ t - xs‖ ≤
        (1 - α / 2) ^ t * ‖x₀ - xs‖ + 4 * η * (1 - (1 - α / 2) ^ t) ∧
      ‖descent v α x₀ t - xs‖ ≤ r := by
  set ρ := 1 - α / 2 with hρ
  have hρ0 : 0 ≤ ρ := by linarith
  have hρ1 : ρ ≤ 1 := by linarith
  induction t with
  | zero => exact ⟨by simp [descent], by simpa [descent] using h₀⟩
  | succ t ih =>
    obtain ⟨hb, hr⟩ := ih
    have hstep := norm_descent_step_le hα0 hα1 hreg hr
    have hpow0 : 0 ≤ ρ ^ (t + 1) := pow_nonneg hρ0 _
    have hpow1 : ρ ^ (t + 1) ≤ 1 := pow_le_one₀ hρ0 hρ1
    have hbound : ‖descent v α x₀ (t + 1) - xs‖ ≤
        ρ ^ (t + 1) * ‖x₀ - xs‖ + 4 * η * (1 - ρ ^ (t + 1)) := by
      have hmul := mul_le_mul_of_nonneg_left hb hρ0
      simp only [descent]
      calc ‖descent v α x₀ t - α • v (descent v α x₀ t) - xs‖
          ≤ ρ * ‖descent v α x₀ t - xs‖ + 2 * α * η := hstep
        _ ≤ ρ * (ρ ^ t * ‖x₀ - xs‖ + 4 * η * (1 - ρ ^ t)) + 2 * α * η := by linarith
        _ = ρ ^ (t + 1) * ‖x₀ - xs‖ + 4 * η * (1 - ρ ^ (t + 1)) := by
          rw [hρ]
          ring
    refine ⟨hbound, hbound.trans ?_⟩
    have h1 : ρ ^ (t + 1) * ‖x₀ - xs‖ ≤ ρ ^ (t + 1) * r := mul_le_mul_of_nonneg_left h₀ hpow0
    have h2 : 4 * η * (1 - ρ ^ (t + 1)) ≤ r * (1 - ρ ^ (t + 1)) :=
      mul_le_mul_of_nonneg_right hηr (by linarith)
    linarith

/-- **S6 as posed.** With `η ≥ 0`, `4η ≤ r` and `0 ≤ α ≤ 1`, from `‖x₀ − x⋆‖ ≤ r` every iterate
stays in the ball and `‖x_t − x⋆‖ ≤ (1 − α/2)ᵗ·‖x₀ − x⋆‖ + 4η`. -/
theorem norm_descent_le_posed {v : E → E} {xs x₀ : E} {η r α : ℝ} (hη : 0 ≤ η)
    (hηr : 4 * η ≤ r) (hα0 : 0 ≤ α) (hα1 : α ≤ 1)
    (hreg : ∀ x, ‖x - xs‖ ≤ r → ‖v x - (x - xs)‖ ≤ ‖x - xs‖ / 2 + 2 * η)
    (h₀ : ‖x₀ - xs‖ ≤ r) (t : ℕ) :
    ‖descent v α x₀ t - xs‖ ≤ (1 - α / 2) ^ t * ‖x₀ - xs‖ + 4 * η ∧
      ‖descent v α x₀ t - xs‖ ≤ r := by
  obtain ⟨hb, hr⟩ := norm_descent_le hηr hα0 hα1 hreg h₀ t
  refine ⟨hb.trans ?_, hr⟩
  have : 0 ≤ (1 - α / 2) ^ t := pow_nonneg (by linarith) _
  nlinarith

/-! ### Witness -/

/-- On `ℝ` with `x⋆ = 0`, the direction `v(x) = 3x/2 + 1/10` meets the regularity inequality with
`η = 1/20` on the ball of radius `1`, and with `α = 1/2` the iterates contract at rate `3/4` to
within `1/5`. -/
example (x₀ : ℝ) (h₀ : ‖x₀‖ ≤ 1) (t : ℕ) :
    ‖descent (fun x : ℝ => 3 / 2 * x + 1 / 10) (1 / 2) x₀ t‖ ≤ (3 / 4) ^ t * ‖x₀‖ + 1 / 5 := by
  have hreg : ∀ x : ℝ, ‖x - 0‖ ≤ 1 →
      ‖(3 / 2 * x + 1 / 10) - (x - 0)‖ ≤ ‖x - 0‖ / 2 + 2 * (1 / 20) := fun x _ => by
    simp only [sub_zero, Real.norm_eq_abs]
    have e : 3 / 2 * x + 1 / 10 - x = x / 2 + 1 / 10 := by ring
    rw [e]
    calc |x / 2 + 1 / 10| ≤ |x / 2| + |1 / 10| := abs_add_le _ _
      _ = |x| / 2 + 2 * (1 / 20) := by
        rw [abs_div, abs_of_pos (by norm_num : (0 : ℝ) < 1 / 10), abs_two]
        ring
  have h := (norm_descent_le_posed (v := fun x : ℝ => 3 / 2 * x + 1 / 10) (xs := 0) (x₀ := x₀)
    (η := 1 / 20) (r := 1) (α := 1 / 2) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    hreg (by simpa using h₀) t).1
  have e1 : (1 - (1 : ℝ) / 2 / 2) = 3 / 4 := by norm_num
  rw [e1] at h
  simp only [sub_zero] at h
  linarith

end TimeFrequency
