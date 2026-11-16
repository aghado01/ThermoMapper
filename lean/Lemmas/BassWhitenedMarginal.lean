import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Data.Matrix.Block
import Mathlib.Data.Matrix.ColumnRowPartitioned
import Mathlib.LinearAlgebra.Matrix.Hadamard
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.SchurComplement
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.MeasureTheory.Integral.Bochner.SumMeasure
import Mathlib.MeasureTheory.Measure.Count
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Module
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# The whitened marginal (S4)

The specification for letting BASS's Gaussian core serve an output inner product `Q ≠ I`
(`WorkingResponse`'s contract: noise covariance `σ²(I ⊗ R)`, `Q = R⁻¹`). With `Ω = L·Lᵀ`, write
`ỹ = L⁻¹y` and `Z̃ = L⁻¹Z`.

* `residual_whiten`: `yᵀΩ⁻¹y − (ZᵀΩ⁻¹y)ᵀ(ZᵀΩ⁻¹Z)⁻¹(ZᵀΩ⁻¹y) = ỹᵀỹ − (Z̃ᵀỹ)ᵀ(Z̃ᵀZ̃)⁻¹(Z̃ᵀỹ)`,
  for any `L`.
* `gaussLogDensity_whiten`: the Gaussian log density of `y` under mean `Zβ` and covariance
  `σ²Ω` is that of `ỹ` under mean `Z̃β` and covariance `σ²I`, less `½·log det Ω`.
* `log_det_one_kronecker`, `inv_one_kronecker`: for `Ω = I_n ⊗ R`, `log det Ω = n·log det R` and
  `Ω⁻¹ = I_n ⊗ R⁻¹`. By `quadForm_one_kronecker`, `ỹᵀỹ = Σᵢ yᵢᵀQyᵢ`.
* `gPrior_logMarginal_whiten`: the engine's terms. Under noise `σ²Ω` and the coefficient prior
  `β ~ N(0, σ²(τ·ZᵀΩ⁻¹Z)⁻¹)` on the whitened Gram, the marginal covariance is
  `σ²(Ω + τ⁻¹·Z(ZᵀΩ⁻¹Z)⁻¹Zᵀ)`. Its log density is the unit-noise closed form on `(ỹ, Z̃)`
  (`gaussLogDensity_gPrior`, `CollapsedGPrior`'s at `b = 1`), less `½·log det Ω`:
  `−(N/2)·log(2πσ²) + (k/2)·log(τ/(τ+1)) − S/(2σ²) − ½·log det Ω`, with
  `S = ỹᵀỹ − ỹᵀP̃ỹ/(τ+1)` and `P̃` the projection onto `Z̃`'s columns.
* `level_whiten`, `level_one_kronecker`, `level_eq_mean`: the level-centred intercept under `Q`.
  The constant column `1` is whitened like any other, so its level is the generalized
  least-squares level `1ᵀΩ⁻¹y / 1ᵀΩ⁻¹1`. Under `Ω = I_n ⊗ R` that is `Σᵢ 1ᵀQyᵢ / (n·1ᵀQ1)`. It
  is the plain mean when `1ᵀQ = λ·1ᵀ` with `λ ≠ 0`.

That `y`'s marginal under the linear-Gaussian model is `N(0, σ²(Ω + τ⁻¹·Z(ZᵀΩ⁻¹Z)⁻¹Zᵀ))` is the
Gaussian convolution, stated here as the covariance and not derived from an integral.
`gaussLogDensity` is the density's formula, with `log` of the determinant. S14
(`marginal_covariance`) proves that covariance on second moments: `β` with second moments
`σ²τ⁻¹(ZᵀΩ⁻¹Z)⁻¹`, `ε` with `σ²Ω`, uncorrelated, give `y = Zβ + ε` exactly that matrix. That the sum
is Gaussian stays a premise.
-/

open Matrix MeasureTheory
open scoped Kronecker

namespace BASS

variable {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]

/-! ### Whitening -/

theorem dotProduct_transpose_mulVec {a b : Type*} [Fintype a] [Fintype b] (A : Matrix a b ℝ)
    (v : b → ℝ) (w : a → ℝ) : v ⬝ᵥ (Aᵀ *ᵥ w) = (A *ᵥ v) ⬝ᵥ w := by
  rw [dotProduct_mulVec, vecMul_transpose]

/-- `(L·Lᵀ)⁻¹ = (L⁻¹)ᵀ·L⁻¹`, for any `L`. -/
theorem inv_mul_transpose (L : Matrix ι ι ℝ) : (L * Lᵀ)⁻¹ = (L⁻¹)ᵀ * L⁻¹ := by
  rw [Matrix.mul_inv_rev, transpose_nonsing_inv]

/-- `vᵀΩ⁻¹w = (L⁻¹v)ᵀ(L⁻¹w)` for `Ω = L·Lᵀ`. -/
theorem quadForm_whiten (L : Matrix ι ι ℝ) (v w : ι → ℝ) :
    v ⬝ᵥ ((L * Lᵀ)⁻¹ *ᵥ w) = (L⁻¹ *ᵥ v) ⬝ᵥ (L⁻¹ *ᵥ w) := by
  rw [inv_mul_transpose, ← mulVec_mulVec, dotProduct_transpose_mulVec]

omit [Fintype κ] [DecidableEq κ] in
/-- The whitened Gram: `ZᵀΩ⁻¹Z = Z̃ᵀZ̃`. -/
theorem gram_whiten (L : Matrix ι ι ℝ) (Z : Matrix ι κ ℝ) :
    Zᵀ * (L * Lᵀ)⁻¹ * Z = (L⁻¹ * Z)ᵀ * (L⁻¹ * Z) := by
  simp only [inv_mul_transpose, transpose_mul, Matrix.mul_assoc]

omit [Fintype κ] [DecidableEq κ] in
/-- The whitened cross products: `ZᵀΩ⁻¹y = Z̃ᵀỹ`. -/
theorem cross_whiten (L : Matrix ι ι ℝ) (Z : Matrix ι κ ℝ) (y : ι → ℝ) :
    Zᵀ *ᵥ ((L * Lᵀ)⁻¹ *ᵥ y) = (L⁻¹ * Z)ᵀ *ᵥ (L⁻¹ *ᵥ y) := by
  simp only [inv_mul_transpose, transpose_mul, ← mulVec_mulVec]

/-- **S4, the residual forms.** For `Ω = L·Lᵀ`,
`yᵀΩ⁻¹y − (ZᵀΩ⁻¹y)ᵀ(ZᵀΩ⁻¹Z)⁻¹(ZᵀΩ⁻¹y) = ỹᵀỹ − (Z̃ᵀỹ)ᵀ(Z̃ᵀZ̃)⁻¹(Z̃ᵀỹ)`, with Mathlib's inverse
(so for any `L` and `Z`). -/
theorem residual_whiten (L Ω : Matrix ι ι ℝ) (hΩ : Ω = L * Lᵀ) (Z : Matrix ι κ ℝ) (y : ι → ℝ) :
    y ⬝ᵥ (Ω⁻¹ *ᵥ y) - (Zᵀ *ᵥ (Ω⁻¹ *ᵥ y)) ⬝ᵥ ((Zᵀ * Ω⁻¹ * Z)⁻¹ *ᵥ (Zᵀ *ᵥ (Ω⁻¹ *ᵥ y))) =
      (L⁻¹ *ᵥ y) ⬝ᵥ (L⁻¹ *ᵥ y) - ((L⁻¹ * Z)ᵀ *ᵥ (L⁻¹ *ᵥ y)) ⬝ᵥ
        (((L⁻¹ * Z)ᵀ * (L⁻¹ * Z))⁻¹ *ᵥ ((L⁻¹ * Z)ᵀ *ᵥ (L⁻¹ *ᵥ y))) := by
  subst hΩ
  rw [quadForm_whiten, gram_whiten, cross_whiten]

/-! ### The Gaussian log density -/

/-- The Gaussian log density of `y` under mean `μ` and covariance `C`,
`−(N/2)·log 2π − ½·log det C − ½·(y − μ)ᵀC⁻¹(y − μ)`. -/
noncomputable def gaussLogDensity (y μ : ι → ℝ) (C : Matrix ι ι ℝ) : ℝ :=
  -(Fintype.card ι : ℝ) / 2 * Real.log (2 * Real.pi) - Real.log C.det / 2 -
    (y - μ) ⬝ᵥ (C⁻¹ *ᵥ (y - μ)) / 2

/-- Congruence: the density under `L·M·Lᵀ` is the density of `L⁻¹y` under `M`, less
`½·log det(L·Lᵀ)`. -/
theorem gaussLogDensity_congr (L M : Matrix ι ι ℝ) (hL : L.det ≠ 0) (hM : M.det ≠ 0)
    (y μ : ι → ℝ) :
    gaussLogDensity y μ (L * M * Lᵀ) =
      gaussLogDensity (L⁻¹ *ᵥ y) (L⁻¹ *ᵥ μ) M - Real.log (L * Lᵀ).det / 2 := by
  have hdet : (L * M * Lᵀ).det = M.det * (L * Lᵀ).det := by
    rw [det_mul, det_mul, det_mul, det_transpose]
    ring
  have hLL : (L * Lᵀ).det ≠ 0 := by
    rw [det_mul, det_transpose]
    exact mul_ne_zero hL hL
  have hinv : (L * M * Lᵀ)⁻¹ = (L⁻¹)ᵀ * (M⁻¹ * L⁻¹) := by
    rw [Matrix.mul_inv_rev, Matrix.mul_inv_rev, transpose_nonsing_inv]
  unfold gaussLogDensity
  rw [hdet, Real.log_mul hM hLL, hinv, ← mulVec_mulVec, ← mulVec_mulVec,
    dotProduct_transpose_mulVec, mulVec_sub]
  ring

omit [DecidableEq κ] in
/-- **S4, the density.** The Gaussian log density of `y` under mean `Zβ` and covariance `σ²Ω`,
`Ω = L·Lᵀ`, is that of `ỹ = L⁻¹y` under mean `Z̃β` and covariance `σ²I`, less `½·log det Ω`. -/
theorem gaussLogDensity_whiten (L : Matrix ι ι ℝ) (hL : L.det ≠ 0) (Z : Matrix ι κ ℝ)
    (β : κ → ℝ) {σ2 : ℝ} (hσ : σ2 ≠ 0) (y : ι → ℝ) :
    gaussLogDensity y (Z *ᵥ β) (σ2 • (L * Lᵀ)) =
      gaussLogDensity (L⁻¹ *ᵥ y) ((L⁻¹ * Z) *ᵥ β) (σ2 • 1) - Real.log (L * Lᵀ).det / 2 := by
  have e : σ2 • (L * Lᵀ) = L * (σ2 • 1) * Lᵀ := by
    rw [Matrix.mul_smul, Matrix.mul_one, Matrix.smul_mul]
  have hM : (σ2 • (1 : Matrix ι ι ℝ)).det ≠ 0 := by
    rw [det_smul, det_one, mul_one]
    exact pow_ne_zero _ hσ
  rw [e, gaussLogDensity_congr L _ hL hM, mulVec_mulVec]

/-! ### The unit-noise g-prior marginal in closed form -/

/-- The projection onto the columns of `X`, `X(XᵀX)⁻¹Xᵀ`, with Mathlib's inverse. -/
noncomputable def proj (X : Matrix ι κ ℝ) : Matrix ι ι ℝ :=
  X * (Xᵀ * X)⁻¹ * Xᵀ

omit [DecidableEq ι] in
theorem proj_mul_self (X : Matrix ι κ ℝ) (h : IsUnit (Xᵀ * X).det) :
    proj X * proj X = proj X := by
  simp only [proj, Matrix.mul_assoc]
  rw [← Matrix.mul_assoc Xᵀ X, ← Matrix.mul_assoc (Xᵀ * X)⁻¹ (Xᵀ * X), nonsing_inv_mul _ h,
    Matrix.one_mul]

/-- `det(I + a·P) = (1 + a)ᵏ` for the projection onto `k` independent columns, by
Weinstein–Aronszajn. -/
theorem det_one_add_smul_proj (X : Matrix ι κ ℝ) (h : IsUnit (Xᵀ * X).det) (a : ℝ) :
    (1 + a • proj X).det = (1 + a) ^ Fintype.card κ := by
  have e : a • proj X = (a • X) * ((Xᵀ * X)⁻¹ * Xᵀ) := by
    rw [proj, Matrix.smul_mul, Matrix.mul_assoc]
  rw [e, det_one_add_mul_comm, Matrix.mul_smul, Matrix.mul_assoc, nonsing_inv_mul _ h,
    show (1 : Matrix κ κ ℝ) + a • 1 = (1 + a) • 1 by rw [add_smul, one_smul], det_smul, det_one,
    mul_one]

/-- `(s·(I + a·P))⁻¹ = s⁻¹·(I − a/(1 + a)·P)`. -/
theorem inv_smul_one_add_smul_proj (X : Matrix ι κ ℝ) (h : IsUnit (Xᵀ * X).det) {s a : ℝ}
    (hs : s ≠ 0) (ha : 1 + a ≠ 0) :
    (s • (1 + a • proj X))⁻¹ = s⁻¹ • (1 - (a / (1 + a)) • proj X) := by
  apply inv_eq_right_inv
  rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul, mul_inv_cancel₀ hs, one_smul]
  have hP := proj_mul_self X h
  set P := proj X
  set d := a / (1 + a) with hd_def
  have hd : a - d - a * d = 0 := by
    rw [hd_def]
    field_simp
    ring
  calc (1 + a • P) * (1 - d • P) = 1 + (a - d - a * d) • P := by
        rw [Matrix.add_mul, Matrix.one_mul, Matrix.mul_sub, Matrix.mul_one, Matrix.smul_mul,
          Matrix.mul_smul, smul_smul, hP, sub_smul, sub_smul]
        abel
    _ = 1 := by rw [hd, zero_smul, add_zero]

/-- **The unit-noise marginal.** With `σ² > 0`, `τ > 0` and `XᵀX` invertible, the log density
of `y` under covariance `σ²(I + τ⁻¹·P)` is
`−(N/2)·log(2πσ²) + (k/2)·log(τ/(τ+1)) − (yᵀy − yᵀPy/(τ+1))/(2σ²)`: `CollapsedGPrior`'s
determinant term and scale at `b = 1`. -/
theorem gaussLogDensity_gPrior (X : Matrix ι κ ℝ) (h : IsUnit (Xᵀ * X).det) {σ2 τ : ℝ}
    (hσ : 0 < σ2) (hτ : 0 < τ) (y : ι → ℝ) :
    gaussLogDensity y 0 (σ2 • (1 + τ⁻¹ • proj X)) =
      -(Fintype.card ι : ℝ) / 2 * Real.log (2 * Real.pi * σ2) +
        (Fintype.card κ : ℝ) / 2 * Real.log (τ / (τ + 1)) -
          (y ⬝ᵥ y - y ⬝ᵥ (proj X *ᵥ y) / (τ + 1)) / (2 * σ2) := by
  have h1a : (1 : ℝ) + τ⁻¹ ≠ 0 := by positivity
  have hlog : Real.log (σ2 • (1 + τ⁻¹ • proj X)).det =
      Fintype.card ι * Real.log σ2 - Fintype.card κ * Real.log (τ / (τ + 1)) := by
    rw [det_smul, det_one_add_smul_proj X h,
      Real.log_mul (pow_ne_zero _ hσ.ne') (pow_ne_zero _ h1a), Real.log_pow, Real.log_pow,
      show (1 + τ⁻¹ : ℝ) = (τ / (τ + 1))⁻¹ by field_simp, Real.log_inv]
    ring
  have hd : τ⁻¹ / (1 + τ⁻¹) = 1 / (τ + 1) := by
    field_simp
  have hquad : (y - 0) ⬝ᵥ ((σ2 • (1 + τ⁻¹ • proj X))⁻¹ *ᵥ (y - 0)) =
      (y ⬝ᵥ y - y ⬝ᵥ (proj X *ᵥ y) / (τ + 1)) / σ2 := by
    rw [sub_zero, inv_smul_one_add_smul_proj X h hσ.ne' h1a, hd, smul_mulVec, sub_mulVec,
      one_mulVec, smul_mulVec, dotProduct_smul, dotProduct_sub, dotProduct_smul]
    simp only [smul_eq_mul]
    ring
  have h2π : Real.log (2 * Real.pi * σ2) = Real.log (2 * Real.pi) + Real.log σ2 :=
    Real.log_mul (by positivity) hσ.ne'
  unfold gaussLogDensity
  rw [hlog, hquad, h2π]
  ring

/-- **S4 in the engine's terms.** Under noise `σ²Ω`, `Ω = L·Lᵀ`, and the coefficient prior on
the whitened Gram `ZᵀΩ⁻¹Z`, the collapsed marginal's log density is the unit-noise closed form
on `(ỹ, Z̃)` less `½·log det Ω`. -/
theorem gPrior_logMarginal_whiten (L : Matrix ι ι ℝ) (hL : L.det ≠ 0) (Z : Matrix ι κ ℝ)
    (h : IsUnit ((L⁻¹ * Z)ᵀ * (L⁻¹ * Z)).det) {σ2 τ : ℝ} (hσ : 0 < σ2) (hτ : 0 < τ)
    (y : ι → ℝ) :
    gaussLogDensity y 0 (σ2 • (L * Lᵀ + τ⁻¹ • (Z * (Zᵀ * (L * Lᵀ)⁻¹ * Z)⁻¹ * Zᵀ))) =
      -(Fintype.card ι : ℝ) / 2 * Real.log (2 * Real.pi * σ2) +
        (Fintype.card κ : ℝ) / 2 * Real.log (τ / (τ + 1)) -
          ((L⁻¹ *ᵥ y) ⬝ᵥ (L⁻¹ *ᵥ y) -
            (L⁻¹ *ᵥ y) ⬝ᵥ (proj (L⁻¹ * Z) *ᵥ (L⁻¹ *ᵥ y)) / (τ + 1)) / (2 * σ2) -
        Real.log (L * Lᵀ).det / 2 := by
  have hLu : IsUnit L.det := isUnit_iff_ne_zero.2 hL
  have hZ : Z = L * (L⁻¹ * Z) := by
    rw [← Matrix.mul_assoc, mul_nonsing_inv _ hLu, Matrix.one_mul]
  have hcov : σ2 • (L * Lᵀ + τ⁻¹ • (Z * (Zᵀ * (L * Lᵀ)⁻¹ * Z)⁻¹ * Zᵀ)) =
      L * (σ2 • (1 + τ⁻¹ • proj (L⁻¹ * Z))) * Lᵀ := by
    rw [gram_whiten]
    set Zt := L⁻¹ * Z
    rw [hZ]
    simp only [proj, transpose_mul, Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul,
      Matrix.smul_mul, Matrix.mul_one, Matrix.mul_assoc, smul_add]
  have hM : (σ2 • (1 + τ⁻¹ • proj (L⁻¹ * Z))).det ≠ 0 := by
    rw [det_smul, det_one_add_smul_proj _ h]
    exact mul_ne_zero (pow_ne_zero _ hσ.ne') (pow_ne_zero _ (by positivity))
  rw [hcov, gaussLogDensity_congr L _ hL hM, mulVec_zero, gaussLogDensity_gPrior _ h hσ hτ]

/-! ### The marginal covariance (S14) -/

/-- A linear image's second moments: `E[X_a·X_b] = M a b` gives
`E[(AX)_i·(AX)_j] = (A·M·Aᵀ)_{ij}`. -/
theorem integral_mulVec_mul_mulVec {m p Ω : Type*} [Fintype p] [MeasurableSpace Ω]
    (μ : Measure Ω) (A : Matrix m p ℝ) (X : Ω → p → ℝ) (M : Matrix p p ℝ)
    (h : ∀ a b, Integrable (fun w => X w a * X w b) μ ∧ ∫ w, X w a * X w b ∂μ = M a b)
    (i j : m) : ∫ w, (A *ᵥ X w) i * (A *ᵥ X w) j ∂μ = (A * M * Aᵀ) i j := by
  have hpt : ∀ w, (A *ᵥ X w) i * (A *ᵥ X w) j =
      ∑ a, ∑ b, A i a * A j b * (X w a * X w b) := by
    intro w
    simp only [mulVec, dotProduct, Finset.sum_mul, Finset.mul_sum]
    conv_lhs => rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    ring
  have hL : ∫ w, (A *ᵥ X w) i * (A *ᵥ X w) j ∂μ = ∑ a, ∑ b, A i a * A j b * M a b := by
    simp_rw [hpt]
    rw [integral_finsetSum _ fun a _ => integrable_finsetSum _ fun b _ => (h a b).1.const_mul _]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [integral_finsetSum _ fun b _ => (h a b).1.const_mul _]
    refine Finset.sum_congr rfl fun b _ => ?_
    rw [integral_const_mul, (h a b).2]
  rw [hL]
  simp only [Matrix.mul_apply, transpose_apply, Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun b _ => Finset.sum_congr rfl fun a _ => ?_
  ring

omit [DecidableEq ι] [DecidableEq κ] in
-- The statement needs no finiteness of `ι`; the proof stacks `(β, ε)` over `κ ⊕ ι`.
set_option linter.unusedFintypeInType false in
/-- **S14, on second moments.** With `E[ββᵀ] = Σ_β`, `E[εεᵀ] = Σ_ε` and `E[βεᵀ] = 0`, the response
`y = Zβ + ε` has `E[yyᵀ] = Z·Σ_β·Zᵀ + Σ_ε`. With zero means these are covariances. -/
theorem integral_secondMoment_response {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω)
    (Z : Matrix ι κ ℝ) (β : Ω → κ → ℝ) (ε : Ω → ι → ℝ) (Sβ : Matrix κ κ ℝ) (Sε : Matrix ι ι ℝ)
    (hββ : ∀ a b, Integrable (fun w => β w a * β w b) μ ∧ ∫ w, β w a * β w b ∂μ = Sβ a b)
    (hεε : ∀ i j, Integrable (fun w => ε w i * ε w j) μ ∧ ∫ w, ε w i * ε w j ∂μ = Sε i j)
    (hβε : ∀ a i, Integrable (fun w => β w a * ε w i) μ ∧ ∫ w, β w a * ε w i ∂μ = 0)
    (i j : ι) :
    ∫ w, (Z *ᵥ β w + ε w) i * (Z *ᵥ β w + ε w) j ∂μ = (Z * Sβ * Zᵀ + Sε) i j := by
  classical
  set X : Ω → κ ⊕ ι → ℝ := fun w => Sum.elim (β w) (ε w) with hX_def
  have hX : ∀ w, Z *ᵥ β w + ε w = fromCols Z 1 *ᵥ X w := fun w => by
    rw [fromCols_mulVec, one_mulVec]
    rfl
  have hM : ∀ a b, Integrable (fun w => X w a * X w b) μ ∧
      ∫ w, X w a * X w b ∂μ = fromBlocks Sβ 0 0 Sε a b := by
    rintro (a | a) (b | b)
    · simpa [hX_def] using hββ a b
    · simpa [hX_def] using hβε a b
    · have hc : (fun w => ε w a * β w b) = fun w => β w b * ε w a :=
        funext fun w => mul_comm _ _
      simp only [hX_def, Sum.elim_inl, Sum.elim_inr, fromBlocks_apply₂₁, Matrix.zero_apply, hc]
      exact hβε b a
    · simpa [hX_def] using hεε a b
  simp_rw [hX]
  rw [integral_mulVec_mul_mulVec μ _ X _ hM, transpose_fromCols, fromCols_mul_fromBlocks,
    fromCols_mul_fromRows]
  simp

/-- **S14.** With `Σ_β = σ²τ⁻¹(ZᵀΩ⁻¹Z)⁻¹` and `Σ_ε = σ²Ω`, uncorrelated, the response's second
moments are `σ²(Ω + τ⁻¹Z(ZᵀΩ⁻¹Z)⁻¹Zᵀ)`, the covariance `gPrior_logMarginal_whiten` reads. -/
theorem marginal_covariance {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) (Z : Matrix ι κ ℝ)
    (Om : Matrix ι ι ℝ) {σ2 τ : ℝ} (β : Ω → κ → ℝ) (ε : Ω → ι → ℝ)
    (hββ : ∀ a b, Integrable (fun w => β w a * β w b) μ ∧
      ∫ w, β w a * β w b ∂μ = ((σ2 * τ⁻¹) • (Zᵀ * Om⁻¹ * Z)⁻¹) a b)
    (hεε : ∀ i j, Integrable (fun w => ε w i * ε w j) μ ∧ ∫ w, ε w i * ε w j ∂μ = (σ2 • Om) i j)
    (hβε : ∀ a i, Integrable (fun w => β w a * ε w i) μ ∧ ∫ w, β w a * ε w i ∂μ = 0)
    (i j : ι) :
    ∫ w, (Z *ᵥ β w + ε w) i * (Z *ᵥ β w + ε w) j ∂μ =
      (σ2 • (Om + τ⁻¹ • (Z * (Zᵀ * Om⁻¹ * Z)⁻¹ * Zᵀ))) i j := by
  rw [integral_secondMoment_response μ Z β ε _ _ hββ hεε hβε]
  suffices hmat : Z * ((σ2 * τ⁻¹) • (Zᵀ * Om⁻¹ * Z)⁻¹) * Zᵀ + σ2 • Om =
      σ2 • (Om + τ⁻¹ • (Z * (Zᵀ * Om⁻¹ * Z)⁻¹ * Zᵀ)) by rw [hmat]
  rw [Matrix.mul_smul, Matrix.smul_mul]
  module

/-! ### The output inner product `Ω = I_n ⊗ R` -/

variable {c o : Type*} [Fintype c] [DecidableEq c] [Fintype o] [DecidableEq o]

theorem det_one_kronecker (R : Matrix o o ℝ) :
    ((1 : Matrix c c ℝ) ⊗ₖ R).det = R.det ^ Fintype.card c := by
  rw [det_kronecker, det_one, one_pow, one_mul]

/-- **S4, the determinant.** `log det(I_n ⊗ R) = n·log det R`. -/
theorem log_det_one_kronecker (R : Matrix o o ℝ) :
    Real.log ((1 : Matrix c c ℝ) ⊗ₖ R).det = Fintype.card c * Real.log R.det := by
  rw [det_one_kronecker, Real.log_pow]

/-- **S4, the inverse.** `(I_n ⊗ R)⁻¹ = I_n ⊗ R⁻¹`. -/
theorem inv_one_kronecker (R : Matrix o o ℝ) :
    ((1 : Matrix c c ℝ) ⊗ₖ R)⁻¹ = 1 ⊗ₖ R⁻¹ := by
  rw [inv_kronecker, inv_one]

omit [DecidableEq o] in
/-- A factor of `R` gives one of `I_n ⊗ R`: `(I ⊗ L_R)(I ⊗ L_R)ᵀ = I ⊗ L_R·L_Rᵀ`. -/
theorem one_kronecker_mul_transpose (LR : Matrix o o ℝ) :
    ((1 : Matrix c c ℝ) ⊗ₖ LR) * ((1 : Matrix c c ℝ) ⊗ₖ LR)ᵀ =
      (1 : Matrix c c ℝ) ⊗ₖ (LR * LRᵀ) := by
  rw [← kroneckerMap_transpose, transpose_one, ← mul_kronecker_mul, Matrix.mul_one]

omit [DecidableEq o] in
theorem one_kronecker_mulVec_apply (Q : Matrix o o ℝ) (w : c × o → ℝ) (i : c) (l : o) :
    (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ w) (i, l) = (Q *ᵥ fun l' => w (i, l')) l := by
  simp only [mulVec, dotProduct, kroneckerMap_apply, Fintype.sum_prod_type]
  rw [Finset.sum_eq_single i (fun b _ hb => by simp [one_apply, Ne.symm hb]) (by simp)]
  simp

omit [DecidableEq o] in
/-- `vᵀ(I_n ⊗ Q)w = Σᵢ vᵢᵀQwᵢ`, with `vᵢ = v (i, ·)`, the case index outer. -/
theorem dotProduct_one_kronecker_mulVec (Q : Matrix o o ℝ) (v w : c × o → ℝ) :
    v ⬝ᵥ (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ w) =
      ∑ i, (fun l => v (i, l)) ⬝ᵥ (Q *ᵥ fun l => w (i, l)) := by
  simp only [dotProduct, Fintype.sum_prod_type, one_kronecker_mulVec_apply]

/-- **S4, the quadratic term.** Whitening by `I ⊗ L_R` gives `ỹᵀỹ = Σᵢ yᵢᵀQyᵢ`,
`Q = (L_R·L_Rᵀ)⁻¹ = R⁻¹`: the contract's `Σᵢ yᵢᵀQyᵢ`. -/
theorem quadForm_one_kronecker (LR : Matrix o o ℝ) (y : c × o → ℝ) :
    (((1 : Matrix c c ℝ) ⊗ₖ LR)⁻¹ *ᵥ y) ⬝ᵥ (((1 : Matrix c c ℝ) ⊗ₖ LR)⁻¹ *ᵥ y) =
      ∑ i, (fun l => y (i, l)) ⬝ᵥ ((LR * LRᵀ)⁻¹ *ᵥ fun l => y (i, l)) := by
  rw [← quadForm_whiten, one_kronecker_mul_transpose, inv_one_kronecker,
    dotProduct_one_kronecker_mulVec]

/-- BASS's direct functional design, the column-wise Khatri–Rao product `A ⊙ H`: column `j` is
`aⱼ ⊗ hⱼ`, rows indexed `(case, location)` with the case outer, as `KhatriRao` lays them out. -/
def khatriRao (A : Matrix c κ ℝ) (H : Matrix o κ ℝ) : Matrix (c × o) κ ℝ :=
  of fun x j => A x.1 j * H x.2 j

omit [DecidableEq ι] [Fintype κ] [DecidableEq κ] [DecidableEq o] in
/-- **The Khatri–Rao Gram.** Under `W = I_n ⊗ Q`, `(A ⊙ H)ᵀW(A ⊙ H) = (AᵀA)∘(HᵀQH)`, the
Hadamard product: `KhatriRao.Gram`'s form, and the whitened Gram of the functional design for
`Q = R⁻¹`. -/
theorem khatriRao_gram (A : Matrix c κ ℝ) (H : Matrix o κ ℝ) (Q : Matrix o o ℝ) :
    (khatriRao A H)ᵀ * ((1 : Matrix c c ℝ) ⊗ₖ Q) * khatriRao A H =
      (Aᵀ * A) ⊙ (Hᵀ * Q * H) := by
  ext j j'
  have hentry : ((khatriRao A H)ᵀ * ((1 : Matrix c c ℝ) ⊗ₖ Q) * khatriRao A H) j j' =
      (fun x => khatriRao A H x j) ⬝ᵥ
        (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ fun x => khatriRao A H x j') := by
    rw [Matrix.mul_assoc]
    rfl
  have hHQH : (Hᵀ * Q * H) j j' = (fun l => H l j) ⬝ᵥ (Q *ᵥ fun l => H l j') := by
    rw [Matrix.mul_assoc]
    rfl
  have hAA : (Aᵀ * A) j j' = ∑ i, A i j * A i j' := rfl
  rw [hentry, dotProduct_one_kronecker_mulVec, hadamard_apply, hHQH, hAA, Finset.sum_mul]
  refine Finset.sum_congr rfl fun i _ => ?_
  have e1 : (fun l => khatriRao A H (i, l) j) = A i j • fun l => H l j := by
    ext l
    simp [khatriRao]
  have e2 : (fun l => khatriRao A H (i, l) j') = A i j' • fun l => H l j' := by
    ext l
    simp [khatriRao]
  rw [e1, e2, mulVec_smul, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul]
  ring

omit [DecidableEq ι] [Fintype κ] [DecidableEq κ] [DecidableEq c] [DecidableEq o] in
/-- Unit output weights: `(A ⊙ H)ᵀ(A ⊙ H) = (AᵀA)∘(HᵀH)`, the Cartesian shortcut. -/
theorem khatriRao_gram_one (A : Matrix c κ ℝ) (H : Matrix o κ ℝ) :
    (khatriRao A H)ᵀ * khatriRao A H = (Aᵀ * A) ⊙ (Hᵀ * H) := by
  classical
  have := khatriRao_gram A H 1
  rwa [one_kronecker_one, Matrix.mul_one, Matrix.mul_one] at this

/-! ### The level-centred intercept under `Q` -/

/-- The whitened constant column's least-squares level is the generalized least-squares level:
`(1̃ᵀỹ)/(1̃ᵀ1̃) = 1ᵀΩ⁻¹y / 1ᵀΩ⁻¹1`. -/
theorem level_whiten (L : Matrix ι ι ℝ) (y : ι → ℝ) :
    ((L⁻¹ *ᵥ fun _ => 1) ⬝ᵥ (L⁻¹ *ᵥ y)) / ((L⁻¹ *ᵥ fun _ => 1) ⬝ᵥ (L⁻¹ *ᵥ fun _ => 1)) =
      ((fun _ => 1) ⬝ᵥ ((L * Lᵀ)⁻¹ *ᵥ y)) / ((fun _ => 1) ⬝ᵥ ((L * Lᵀ)⁻¹ *ᵥ fun _ => 1)) := by
  rw [quadForm_whiten, quadForm_whiten]

omit [DecidableEq o] in
/-- **S4, the level under `Q`.** With `W = I_n ⊗ Q`, the level is `Σᵢ 1ᵀQyᵢ / (n·1ᵀQ1)`: the
`Q`-weighted level of the mean curve, with mass `n·1ᵀQ1`. -/
theorem level_one_kronecker (Q : Matrix o o ℝ) (y : c × o → ℝ) :
    ((fun _ => 1) ⬝ᵥ (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ y)) /
        ((fun _ => 1) ⬝ᵥ (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ fun _ => 1)) =
      (∑ i, (fun _ => 1) ⬝ᵥ (Q *ᵥ fun l => y (i, l))) /
        (Fintype.card c * ((fun _ => 1) ⬝ᵥ (Q *ᵥ fun _ => 1))) := by
  rw [dotProduct_one_kronecker_mulVec, dotProduct_one_kronecker_mulVec, Finset.sum_const,
    Finset.card_univ, nsmul_eq_mul]

omit [DecidableEq o] in
/-- **S4, when the level is the mean.** If `1ᵀQ = λ·1ᵀ` with `λ ≠ 0` (constant column sums), the
level under `W = I_n ⊗ Q` is the plain mean of all `n·m` values. -/
theorem level_eq_mean (Q : Matrix o o ℝ) {lam : ℝ} (hQ : (fun _ => (1 : ℝ)) ᵥ* Q = lam • fun _ => 1)
    (hlam : lam ≠ 0) [Nonempty c] [Nonempty o] (y : c × o → ℝ) :
    ((fun _ => 1) ⬝ᵥ (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ y)) /
        ((fun _ => 1) ⬝ᵥ (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ fun _ => 1)) =
      (∑ x, y x) / (Fintype.card c * Fintype.card o) := by
  have key : ∀ v : o → ℝ, (fun _ => (1 : ℝ)) ⬝ᵥ (Q *ᵥ v) = lam * ∑ l, v l := fun v => by
    rw [dotProduct_mulVec, hQ, smul_dotProduct, smul_eq_mul]
    simp [dotProduct]
  rw [level_one_kronecker]
  simp only [key, ← Finset.mul_sum, Fintype.sum_prod_type, Finset.sum_const, Finset.card_univ,
    nsmul_eq_mul, mul_one]
  have hc : (Fintype.card c : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Fintype.card_ne_zero
  have ho : (Fintype.card o : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Fintype.card_ne_zero
  field_simp

/-! ### Witnesses -/

/-- S14's hypotheses on a concrete space: four outcomes under the counting measure, one
coefficient and one observation, `Z = Ω = σ² = τ = 1`, `β = ½(1, 1, −1, −1)` and
`ε = ½(1, −1, 1, −1)`, uncorrelated with unit second moments. -/
example := marginal_covariance (Measure.count : Measure (Fin 4)) (1 : Matrix (Fin 1) (Fin 1) ℝ) 1
  (σ2 := 1) (τ := 1)
  (fun w _ => if w.val < 2 then (1 / 2 : ℝ) else -1 / 2)
  (fun w _ => if w.val % 2 = 0 then (1 / 2 : ℝ) else -1 / 2)
  (fun a b => ⟨Integrable.of_finite, by
    rw [integral_fintype Integrable.of_finite]
    fin_cases a; fin_cases b
    simp [Fin.sum_univ_succ, Measure.real]
    norm_num⟩)
  (fun a b => ⟨Integrable.of_finite, by
    rw [integral_fintype Integrable.of_finite]
    fin_cases a; fin_cases b
    simp [Fin.sum_univ_succ, Measure.real]
    norm_num⟩)
  (fun a b => ⟨Integrable.of_finite, by
    rw [integral_fintype Integrable.of_finite]
    simp [Fin.sum_univ_succ, Measure.real]
    norm_num⟩)

/-- The engine-terms statement's hypotheses hold on a correlated pair: `L = [[1, 0], [1, 1]]`,
so `Ω = [[1, 1], [1, 2]]`, the design `Z = (1, 1)ᵀ`, `σ² = 1` and `τ = 1`. -/
example (y : Fin 2 → ℝ) :=
  have hinv : (!![(1 : ℝ), 0; 1, 1])⁻¹ = !![1, 0; -1, 1] := by
    apply inv_eq_left_inv
    ext i j
    fin_cases i <;> fin_cases j <;> simp [Matrix.mul_apply, Fin.sum_univ_two]
  have hG : IsUnit ((!![(1 : ℝ), 0; 1, 1]⁻¹ * !![(1 : ℝ); 1])ᵀ *
      (!![(1 : ℝ), 0; 1, 1]⁻¹ * !![(1 : ℝ); 1])).det := by
    rw [hinv]
    exact isUnit_iff_ne_zero.2 (by simp [Matrix.mul_apply, Fin.sum_univ_two])
  gPrior_logMarginal_whiten !![(1 : ℝ), 0; 1, 1] (by simp [det_fin_two_of]) !![(1 : ℝ); 1] hG
    one_pos one_pos y

/-- The Khatri–Rao Gram on two cases, two locations and one column. -/
example : (khatriRao !![(1 : ℝ); 2] !![(3 : ℝ); 4])ᵀ * khatriRao !![(1 : ℝ); 2] !![(3 : ℝ); 4] =
    (!![(1 : ℝ); 2]ᵀ * !![(1 : ℝ); 2]) ⊙ (!![(3 : ℝ); 4]ᵀ * !![(3 : ℝ); 4]) :=
  khatriRao_gram_one _ _

/-- The level equals the mean under `Q = I` on two cases and two locations. -/
example (y : Fin 2 × Fin 2 → ℝ) :
    ((fun _ => 1) ⬝ᵥ (((1 : Matrix (Fin 2) (Fin 2) ℝ) ⊗ₖ (1 : Matrix (Fin 2) (Fin 2) ℝ)) *ᵥ y)) /
        ((fun _ => 1) ⬝ᵥ (((1 : Matrix (Fin 2) (Fin 2) ℝ) ⊗ₖ (1 : Matrix (Fin 2) (Fin 2) ℝ)) *ᵥ
          fun _ => 1)) =
      (∑ x, y x) / (Fintype.card (Fin 2) * Fintype.card (Fin 2)) :=
  level_eq_mean (1 : Matrix (Fin 2) (Fin 2) ℝ) (lam := 1) (by simp) one_ne_zero y

end BASS
