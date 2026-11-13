import Enthymemata.BassSensitivity
import Mathlib.Algebra.Order.Star.Real
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Active subspaces — the matrix Poincaré check (AS-2)

The statements of `Protolemmata/BASS/active-subspace-lemmas.md`, AS-2, on a finite product law.

* **(a) The reduction.** `loewner_iff`: `A ≼ B` exactly when `vᵀAv ≤ vᵀBv` for every `v`.
  `covForm_eq_var`: `vᵀ Cov(f) v = Var(vᵀf)`.
  `energyForm_eq`: `vᵀ E[J D Jᵀ] v = Σᵢ cᵢ E[(Σₐ vₐJₐᵢ)²]`.
* **(b) The product step.** `meanOf_mul_sub_coordMean`: Efron–Stein's summand is `E[Varᵢ g]`.
  `var_le_of_coord_poincare`: Efron–Stein (`BASS.efron_stein`) and each coordinate's own inequality
  `Varᵢ g ≤ cᵢ eᵢ` give `Var g ≤ Σᵢ cᵢ E[eᵢ]`, each coordinate keeping its own constant;
  `coordwise_le_scalar`: never looser than the largest constant taken for all.
  `matrix_poincare`: so `Cov(f) ≼ E[J D Jᵀ]` when the one-variable inequality holds for every linear
  combination `vᵀf`, whose partials are `vᵀJ`.
* **(c)** `energy_box`: `J diag(wᵢ² r) Jᵀ = r · (J diag w)(J diag w)ᵀ`, so at `r = 1/π²` the
  uniform's
  bound is `G/π²` in the box metric.

The one-variable inequality is a hypothesis (`hcoord`): Mathlib carries no Poincaré or Wirtinger
inequality, and a finite product has no derivative, so the partials enter as data `J`. The uniform's
constant `(b − a)²/π²` and the truncated normal's `min(σ², (b − a)²/π²)` rest on the numeric checks.
-/

open Matrix Finset BASS

namespace ActiveSubspace

/-! ### (a) The reduction -/

section Reduction

variable {o : Type*} [Fintype o]

/-- **AS-2(a).** For symmetric `A` and `B`, `A ≼ B` exactly when every quadratic form is ordered. -/
theorem loewner_iff (A B : Matrix o o ℝ) (hA : A.IsHermitian) (hB : B.IsHermitian) :
    (B - A).PosSemidef ↔ ∀ v : o → ℝ, v ⬝ᵥ (A *ᵥ v) ≤ v ⬝ᵥ (B *ᵥ v) := by
  rw [posSemidef_iff_dotProduct_mulVec]
  simp only [star_trivial, sub_mulVec, dotProduct_sub, sub_nonneg]
  exact ⟨fun h => h.2, fun h => ⟨hB.sub hA, h⟩⟩

end Reduction

section Finite

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {X : ι → Type*} [∀ i, Fintype (X i)]
variable {o : Type*} [Fintype o]

theorem meanOf_add (μ : ∀ i, X i → ℝ) (f g : (∀ j, X j) → ℝ) :
    meanOf μ (fun x => f x + g x) = meanOf μ f + meanOf μ g := by
  simp only [meanOf, mul_add, Finset.sum_add_distrib]

theorem meanOf_sub (μ : ∀ i, X i → ℝ) (f g : (∀ j, X j) → ℝ) :
    meanOf μ (fun x => f x - g x) = meanOf μ f - meanOf μ g := by
  simp only [meanOf, mul_sub, Finset.sum_sub_distrib]

theorem meanOf_const_mul (μ : ∀ i, X i → ℝ) (c : ℝ) (f : (∀ j, X j) → ℝ) :
    meanOf μ (fun x => c * f x) = c * meanOf μ f := by
  simp only [meanOf, Finset.mul_sum]
  exact Finset.sum_congr rfl fun _ _ => by ring

theorem meanOf_sum {κ : Type*} (s : Finset κ) (μ : ∀ i, X i → ℝ) (f : κ → (∀ j, X j) → ℝ) :
    meanOf μ (fun x => ∑ a ∈ s, f a x) = ∑ a ∈ s, meanOf μ (f a) := by
  simp only [meanOf, Finset.mul_sum]
  exact Finset.sum_comm

theorem meanOf_mono (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y) {f g : (∀ j, X j) → ℝ}
    (h : ∀ x, f x ≤ g x) : meanOf μ f ≤ meanOf μ g :=
  Finset.sum_le_sum fun x _ =>
    mul_le_mul_of_nonneg_left (h x) (Finset.prod_nonneg fun i _ => h₀ i _)

/-- The variance under the product law. -/
noncomputable def varOf (μ : ∀ i, X i → ℝ) (g : (∀ j, X j) → ℝ) : ℝ :=
  meanOf μ (fun x => g x ^ 2) - meanOf μ g ^ 2

/-- The covariance matrix of `m` outputs under the product law. -/
noncomputable def covOf (μ : ∀ i, X i → ℝ) (f : (∀ j, X j) → o → ℝ) : Matrix o o ℝ :=
  fun a b => meanOf μ (fun x => f x a * f x b) -
    meanOf μ (fun x => f x a) * meanOf μ (fun x => f x b)

/-- `E[J D Jᵀ]`, with `D = diag c` and the partials `J x a i` given as data. -/
noncomputable def energyOf (μ : ∀ i, X i → ℝ) (J : (∀ j, X j) → Matrix o ι ℝ) (c : ι → ℝ) :
    Matrix o o ℝ :=
  fun a b => ∑ i, c i * meanOf μ (fun x => J x a i * J x b i)

/-- `Varᵢ g`: the variance over coordinate `i` with the others held. -/
noncomputable def condVar (μ : ∀ i, X i → ℝ) (i : ι) (g : (∀ j, X j) → ℝ) (x : ∀ j, X j) : ℝ :=
  coordMean μ i (fun z => g z ^ 2) x - coordMean μ i g x ^ 2

omit [Fintype o] in
theorem covOf_isHermitian (μ : ∀ i, X i → ℝ) (f : (∀ j, X j) → o → ℝ) :
    (covOf μ f).IsHermitian := by
  ext a b
  simp only [conjTranspose_apply, star_trivial, covOf, mul_comm]

omit [Fintype o] in
theorem energyOf_isHermitian (μ : ∀ i, X i → ℝ) (J : (∀ j, X j) → Matrix o ι ℝ) (c : ι → ℝ) :
    (energyOf μ J c).IsHermitian := by
  ext a b
  simp only [conjTranspose_apply, star_trivial, energyOf, mul_comm (J _ a _)]

omit [DecidableEq ι] [∀ i, Fintype (X i)] in
/-- A quadratic form as a double sum. -/
theorem quadForm_eq (M : Matrix o o ℝ) (v : o → ℝ) :
    v ⬝ᵥ (M *ᵥ v) = ∑ a, ∑ b, v a * v b * M a b := by
  simp only [dotProduct, mulVec, Finset.mul_sum]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

/-- `E[(Σₐ vₐ hₐ)²] = Σₐ Σ_b vₐ v_b E[hₐ h_b]`. -/
theorem meanOf_sq_lin (μ : ∀ i, X i → ℝ) (h : o → (∀ j, X j) → ℝ) (v : o → ℝ) :
    meanOf μ (fun x => (∑ a, v a * h a x) ^ 2) =
      ∑ a, ∑ b, v a * v b * meanOf μ (fun x => h a x * h b x) := by
  have hsq : (fun x => (∑ a, v a * h a x) ^ 2) =
      fun x => ∑ a, ∑ b, v a * v b * (h a x * h b x) := by
    funext x
    rw [sq, Finset.sum_mul_sum]
    exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  rw [hsq, meanOf_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [meanOf_sum]
  exact Finset.sum_congr rfl fun b _ => meanOf_const_mul μ _ _

/-- **AS-2(a).** `vᵀ Cov(f) v = Var(vᵀf)`. -/
theorem covForm_eq_var (μ : ∀ i, X i → ℝ) (f : (∀ j, X j) → o → ℝ) (v : o → ℝ) :
    v ⬝ᵥ (covOf μ f *ᵥ v) = varOf μ (fun x => v ⬝ᵥ f x) := by
  have hm : meanOf μ (fun x => v ⬝ᵥ f x) = ∑ a, v a * meanOf μ (fun x => f x a) := by
    simp only [dotProduct]
    rw [meanOf_sum]
    exact Finset.sum_congr rfl fun a _ => meanOf_const_mul μ (v a) (fun x => f x a)
  have h1 : meanOf μ (fun x => (v ⬝ᵥ f x) ^ 2) =
      ∑ a, ∑ b, v a * v b * meanOf μ (fun x => f x a * f x b) :=
    meanOf_sq_lin μ (fun a x => f x a) v
  have h2 : (∑ a, v a * meanOf μ (fun x => f x a)) ^ 2 =
      ∑ a, ∑ b, v a * v b * (meanOf μ (fun x => f x a) * meanOf μ (fun x => f x b)) := by
    rw [sq, Finset.sum_mul_sum]
    exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  rw [quadForm_eq, varOf, hm, h2, h1, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun b _ => ?_
  simp only [covOf]
  ring

/-- **AS-2(a).** `vᵀ E[J D Jᵀ] v = Σᵢ cᵢ E[(Σₐ vₐ Jₐᵢ)²]`, the partials of `vᵀf`. -/
theorem energyForm_eq (μ : ∀ i, X i → ℝ) (J : (∀ j, X j) → Matrix o ι ℝ) (c : ι → ℝ)
    (v : o → ℝ) :
    v ⬝ᵥ (energyOf μ J c *ᵥ v) = ∑ i, c i * meanOf μ (fun x => (∑ a, v a * J x a i) ^ 2) := by
  have h : ∀ i, meanOf μ (fun x => (∑ a, v a * J x a i) ^ 2) =
      ∑ a, ∑ b, v a * v b * meanOf μ (fun x => J x a i * J x b i) :=
    fun i => meanOf_sq_lin μ (fun a x => J x a i) v
  have hL : v ⬝ᵥ (energyOf μ J c *ᵥ v) =
      ∑ a, ∑ b, ∑ i, v a * v b * (c i * meanOf μ (fun x => J x a i * J x b i)) := by
    rw [quadForm_eq]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    simp only [energyOf, Finset.mul_sum]
  have hR : ∑ i, c i * meanOf μ (fun x => (∑ a, v a * J x a i) ^ 2) =
      ∑ i, ∑ a, ∑ b, v a * v b * (c i * meanOf μ (fun x => J x a i * J x b i)) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [h i, Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun b _ => by ring
  have hswap : ∀ a, ∑ b, ∑ i, v a * v b * (c i * meanOf μ (fun x => J x a i * J x b i)) =
      ∑ i, ∑ b, v a * v b * (c i * meanOf μ (fun x => J x a i * J x b i)) :=
    fun a => Finset.sum_comm
  rw [hL, hR]
  simp only [hswap]
  exact Finset.sum_comm

/-! ### (b) The product step -/

/-- Efron–Stein's summand `E[g (g − Eᵢ g)]` is `E[Varᵢ g]`. -/
theorem meanOf_mul_sub_coordMean (μ : ∀ i, X i → ℝ) (hμ : ∀ i, ∑ y, μ i y = 1) (i : ι)
    (g : (∀ j, X j) → ℝ) :
    meanOf μ (fun x => g x * (g x - coordMean μ i g x)) = meanOf μ (condVar μ i g) := by
  have hl : meanOf μ (fun x => g x * (g x - coordMean μ i g x)) =
      meanOf μ (fun x => g x ^ 2) - meanOf μ (fun x => g x * coordMean μ i g x) := by
    rw [← meanOf_sub]
    congr 1
    funext x
    ring
  have hr : meanOf μ (condVar μ i g) =
      meanOf μ (fun x => g x ^ 2) - meanOf μ (fun x => coordMean μ i g x * coordMean μ i g x) := by
    have hc : condVar μ i g =
        fun x => coordMean μ i (fun z => g z ^ 2) x - coordMean μ i g x ^ 2 := rfl
    rw [hc, meanOf_sub, meanOf_coordMean μ hμ i (fun z => g z ^ 2)]
    simp only [sq]
  have hsa := coordMean_selfAdjoint μ hμ i g (coordMean μ i g)
  rw [coordMean_idem μ hμ i g] at hsa
  rw [hl, hr, hsa]

/-- **AS-2(b), the product step.** If each coordinate satisfies its own one-variable inequality
`Varᵢ g ≤ cᵢ · eᵢ` (with the other coordinates held), then `Var g ≤ Σᵢ cᵢ E[eᵢ]`: Efron–Stein puts
each coordinate's constant on its own term. -/
theorem var_le_of_coord_poincare (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y)
    (hμ : ∀ i, ∑ y, μ i y = 1) (g : (∀ j, X j) → ℝ) (c : ι → ℝ) (e : ι → (∀ j, X j) → ℝ)
    (hcoord : ∀ i x, condVar μ i g x ≤ c i * e i x) :
    varOf μ g ≤ ∑ i, c i * meanOf μ (e i) := by
  have hes := efron_stein μ h₀ hμ g
  simp only [meanOf_mul_sub_coordMean μ hμ] at hes
  refine hes.trans (Finset.sum_le_sum fun i _ => ?_)
  rw [← meanOf_const_mul]
  exact meanOf_mono μ h₀ (hcoord i)

/-- **AS-2(b).** The coordinate-wise bound is never looser than one constant `C ≥ cᵢ` for all. -/
theorem coordwise_le_scalar (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y) (c : ι → ℝ) (C : ℝ)
    (hc : ∀ i, c i ≤ C) (e : ι → (∀ j, X j) → ℝ) (he : ∀ i x, 0 ≤ e i x) :
    ∑ i, c i * meanOf μ (e i) ≤ C * ∑ i, meanOf μ (e i) := by
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_right (hc i) ?_
  have h := meanOf_mono μ h₀ (f := fun _ => (0 : ℝ)) (g := e i) (he i)
  simpa [meanOf] using h

/-- **AS-2, the matrix check.** If the one-variable inequality holds coordinate by coordinate for
every linear combination `vᵀf`, whose partial along `i` is `Σₐ vₐ Jₐᵢ`, then
`Cov(f) ≼ E[J D Jᵀ]`. -/
theorem matrix_poincare (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y) (hμ : ∀ i, ∑ y, μ i y = 1)
    (f : (∀ j, X j) → o → ℝ) (J : (∀ j, X j) → Matrix o ι ℝ) (c : ι → ℝ)
    (hcoord : ∀ (v : o → ℝ) i x,
      condVar μ i (fun z => v ⬝ᵥ f z) x ≤ c i * (∑ a, v a * J x a i) ^ 2) :
    (energyOf μ J c - covOf μ f).PosSemidef := by
  rw [loewner_iff _ _ (covOf_isHermitian μ f) (energyOf_isHermitian μ J c)]
  intro v
  rw [covForm_eq_var, energyForm_eq]
  exact var_le_of_coord_poincare μ h₀ hμ _ c (fun i x => (∑ a, v a * J x a i) ^ 2) (hcoord v)

end Finite

/-! ### (c) The box metric -/

/-- **AS-2(c).** `J diag(wᵢ² r) Jᵀ = r · (J diag w)(J diag w)ᵀ`: with `r = 1/π²` and `w` the box
widths, the uniform's bound `E[J D Jᵀ]` is `G/π²` with `G` read in the box metric. -/
theorem energy_box {o ι : Type*} [Fintype ι] [DecidableEq ι] (J : Matrix o ι ℝ) (w : ι → ℝ)
    (r : ℝ) :
    J * diagonal (fun i => w i ^ 2 * r) * Jᵀ = r • ((J * diagonal w) * (J * diagonal w)ᵀ) := by
  ext a b
  rw [Matrix.mul_apply, Matrix.smul_apply, Matrix.mul_apply, smul_eq_mul, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp only [mul_diagonal, transpose_apply]
  ring

end ActiveSubspace
