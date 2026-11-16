import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.Real.Sqrt
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# What the range restricted isometry condition buys (S19)

Hand and Voroninski (arXiv 1705.07576), "Main Results". `A` satisfies the range restricted
isometry condition (RRIC) with constant `ε` on the range of `G` when, for all codes,
`|⟨A(Gx₁ − Gx₂), A(Gx₃ − Gx₄)⟩ − ⟨Gx₁ − Gx₂, Gx₃ − Gx₄⟩| ≤ ε·‖Gx₁ − Gx₂‖·‖Gx₃ − Gx₄‖`.

* `rric_sandwich`: on the diagonal `x₃ = x₁`, `x₄ = x₂`,
  `(1 − ε)·‖Gx₁ − Gx₂‖² ≤ ‖A(Gx₁ − Gx₂)‖² ≤ (1 + ε)·‖Gx₁ − Gx₂‖²`.
* `rric_injective`: for `ε < 1`, `A∘G` separates every pair that `G` separates.
* `rric_objective_sandwich`: `(1 − ε)·f_I ≤ f_A ≤ (1 + ε)·f_I` for the compressed and uncompressed
  objectives `f_A(x) = ½‖A(Gx − Gx₀)‖²` and `f_I(x) = ½‖Gx − Gx₀‖²`.
-/

open RealInnerProductSpace

namespace GenerativePrior

variable {X E F : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [NormedAddCommGroup F]
  [InnerProductSpace ℝ F]

/-- The range restricted isometry condition with constant `ε`. -/
def RangeRestrictedIsometry (A : E →ₗ[ℝ] F) (G : X → E) (ε : ℝ) : Prop :=
  ∀ x₁ x₂ x₃ x₄ : X, |⟪A (G x₁ - G x₂), A (G x₃ - G x₄)⟫ - ⟪G x₁ - G x₂, G x₃ - G x₄⟫| ≤
    ε * ‖G x₁ - G x₂‖ * ‖G x₃ - G x₄‖

/-- **S19, the sandwich.** -/
theorem rric_sandwich {A : E →ₗ[ℝ] F} {G : X → E} {ε : ℝ} (h : RangeRestrictedIsometry A G ε)
    (x₁ x₂ : X) :
    (1 - ε) * ‖G x₁ - G x₂‖ ^ 2 ≤ ‖A (G x₁ - G x₂)‖ ^ 2 ∧
      ‖A (G x₁ - G x₂)‖ ^ 2 ≤ (1 + ε) * ‖G x₁ - G x₂‖ ^ 2 := by
  have h' := h x₁ x₂ x₁ x₂
  rw [real_inner_self_eq_norm_sq, real_inner_self_eq_norm_sq, abs_le] at h'
  constructor <;> nlinarith [h'.1, h'.2]

/-- **S19, injectivity.** For `ε < 1`, `A∘G` separates every pair that `G` separates. -/
theorem rric_injective {A : E →ₗ[ℝ] F} {G : X → E} {ε : ℝ} (h : RangeRestrictedIsometry A G ε)
    (hε : ε < 1) {x₁ x₂ : X} (hG : G x₁ ≠ G x₂) : A (G x₁) ≠ A (G x₂) := by
  intro hA
  have hs := (rric_sandwich h x₁ x₂).1
  rw [map_sub, hA, sub_self, norm_zero] at hs
  have hpos : 0 < ‖G x₁ - G x₂‖ := norm_pos_iff.2 (sub_ne_zero.2 hG)
  nlinarith [mul_pos (sub_pos.2 hε) (pow_pos hpos 2)]

/-- Compression preserves injectivity on any set of codes on which the generator is injective. -/
theorem rric_injOn {A : E →ₗ[ℝ] F} {G : X → E} {ε : ℝ} (h : RangeRestrictedIsometry A G ε)
    (hε : ε < 1) {s : Set X} (hG : Set.InjOn G s) : Set.InjOn (A ∘ G) s := by
  intro x hx y hy hA
  exact hG hx hy (by_contra fun hne ↦ rric_injective h hε hne hA)

/-- The compressed objective `f_A(x) = ½‖A(Gx − Gx₀)‖²`. -/
noncomputable def compressedObjective (A : E →ₗ[ℝ] F) (G : X → E) (x₀ x : X) : ℝ :=
  ‖A (G x - G x₀)‖ ^ 2 / 2

/-- The uncompressed objective `f_I(x) = ½‖Gx − Gx₀‖²`. -/
noncomputable def uncompressedObjective (G : X → E) (x₀ x : X) : ℝ :=
  ‖G x - G x₀‖ ^ 2 / 2

/-- **S19, the objectives.** `(1 − ε)·f_I ≤ f_A ≤ (1 + ε)·f_I`. -/
theorem rric_objective_sandwich {A : E →ₗ[ℝ] F} {G : X → E} {ε : ℝ}
    (h : RangeRestrictedIsometry A G ε) (x₀ x : X) :
    (1 - ε) * uncompressedObjective G x₀ x ≤ compressedObjective A G x₀ x ∧
      compressedObjective A G x₀ x ≤ (1 + ε) * uncompressedObjective G x₀ x := by
  obtain ⟨h1, h2⟩ := rric_sandwich h x x₀
  unfold compressedObjective uncompressedObjective
  constructor <;> linarith

/-! ### Witness -/

/-- A scaled identity `√(3/2)·I` on `ℝ` meets the condition with `ε = ½` for any `G`, so the
injectivity statement's hypotheses are satisfiable with `0 < ε < 1`. -/
theorem rric_scaled_id (G : X → ℝ) :
    RangeRestrictedIsometry (√(3 / 2) • LinearMap.id : ℝ →ₗ[ℝ] ℝ) G (1 / 2) := by
  intro x₁ x₂ x₃ x₄
  have hc : √(3 / 2) * √(3 / 2) = (3 / 2 : ℝ) := Real.mul_self_sqrt (by norm_num)
  simp only [LinearMap.smul_apply, LinearMap.id_coe, id, real_inner_smul_left,
    real_inner_smul_right, ← mul_assoc, hc]
  have hcs := abs_real_inner_le_norm (G x₁ - G x₂) (G x₃ - G x₄)
  rw [show 3 / 2 * ⟪G x₁ - G x₂, G x₃ - G x₄⟫ - ⟪G x₁ - G x₂, G x₃ - G x₄⟫ =
    1 / 2 * ⟪G x₁ - G x₂, G x₃ - G x₄⟫ by ring, abs_mul, abs_of_pos (by norm_num : (0 : ℝ) < 1 / 2)]
  nlinarith [hcs]

theorem rric_scaled_id_separates :
    (√(3 / 2) • LinearMap.id : ℝ →ₗ[ℝ] ℝ) 1 ≠ (√(3 / 2) • LinearMap.id : ℝ →ₗ[ℝ] ℝ) 2 :=
  rric_injective (rric_scaled_id id) (by norm_num) (show (1 : ℝ) ≠ 2 by norm_num)

#print axioms rric_sandwich
#print axioms rric_injective
#print axioms rric_injOn
#print axioms rric_objective_sandwich
#print axioms rric_scaled_id
#print axioms rric_scaled_id_separates

end GenerativePrior
