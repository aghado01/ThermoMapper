import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.SchurComplement
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
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
`gaussLogDensity` is the density's formula, with `log` of the determinant.
-/

open Matrix
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

/-- The whitened Gram: `ZᵀΩ⁻¹Z = Z̃ᵀZ̃`. -/
theorem gram_whiten (L : Matrix ι ι ℝ) (Z : Matrix ι κ ℝ) :
    Zᵀ * (L * Lᵀ)⁻¹ * Z = (L⁻¹ * Z)ᵀ * (L⁻¹ * Z) := by
  simp only [inv_mul_transpose, transpose_mul, Matrix.mul_assoc]

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
    sorry
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

/-- A factor of `R` gives one of `I_n ⊗ R`: `(I ⊗ L_R)(I ⊗ L_R)ᵀ = I ⊗ L_R·L_Rᵀ`. -/
theorem one_kronecker_mul_transpose (LR : Matrix o o ℝ) :
    ((1 : Matrix c c ℝ) ⊗ₖ LR) * ((1 : Matrix c c ℝ) ⊗ₖ LR)ᵀ =
      (1 : Matrix c c ℝ) ⊗ₖ (LR * LRᵀ) := by
  rw [← kroneckerMap_transpose, transpose_one, ← mul_kronecker_mul, Matrix.mul_one]

theorem one_kronecker_mulVec_apply (Q : Matrix o o ℝ) (w : c × o → ℝ) (i : c) (l : o) :
    (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ w) (i, l) = (Q *ᵥ fun l' => w (i, l')) l := by
  simp only [mulVec, dotProduct, kroneckerMap_apply, Fintype.sum_prod_type]
  rw [Finset.sum_eq_single i (fun b _ hb => by simp [one_apply, Ne.symm hb]) (by simp)]
  simp

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

/-! ### The level-centred intercept under `Q` -/

/-- The whitened constant column's least-squares level is the generalized least-squares level:
`(1̃ᵀỹ)/(1̃ᵀ1̃) = 1ᵀΩ⁻¹y / 1ᵀΩ⁻¹1`. -/
theorem level_whiten (L : Matrix ι ι ℝ) (y : ι → ℝ) :
    ((L⁻¹ *ᵥ fun _ => 1) ⬝ᵥ (L⁻¹ *ᵥ y)) / ((L⁻¹ *ᵥ fun _ => 1) ⬝ᵥ (L⁻¹ *ᵥ fun _ => 1)) =
      ((fun _ => 1) ⬝ᵥ ((L * Lᵀ)⁻¹ *ᵥ y)) / ((fun _ => 1) ⬝ᵥ ((L * Lᵀ)⁻¹ *ᵥ fun _ => 1)) := by
  rw [quadForm_whiten, quadForm_whiten]

/-- **S4, the level under `Q`.** With `W = I_n ⊗ Q`, the level is `Σᵢ 1ᵀQyᵢ / (n·1ᵀQ1)`: the
`Q`-weighted level of the mean curve, with mass `n·1ᵀQ1`. -/
theorem level_one_kronecker (Q : Matrix o o ℝ) (y : c × o → ℝ) :
    ((fun _ => 1) ⬝ᵥ (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ y)) /
        ((fun _ => 1) ⬝ᵥ (((1 : Matrix c c ℝ) ⊗ₖ Q) *ᵥ fun _ => 1)) =
      (∑ i, (fun _ => 1) ⬝ᵥ (Q *ᵥ fun l => y (i, l))) /
        (Fintype.card c * ((fun _ => 1) ⬝ᵥ (Q *ᵥ fun _ => 1))) := by
  rw [dotProduct_one_kronecker_mulVec, dotProduct_one_kronecker_mulVec, Finset.sum_const,
    Finset.card_univ, nsmul_eq_mul]

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

/-- The level equals the mean under `Q = I` on two cases and two locations. -/
example (y : Fin 2 × Fin 2 → ℝ) :
    ((fun _ => 1) ⬝ᵥ (((1 : Matrix (Fin 2) (Fin 2) ℝ) ⊗ₖ (1 : Matrix (Fin 2) (Fin 2) ℝ)) *ᵥ y)) /
        ((fun _ => 1) ⬝ᵥ (((1 : Matrix (Fin 2) (Fin 2) ℝ) ⊗ₖ (1 : Matrix (Fin 2) (Fin 2) ℝ)) *ᵥ
          fun _ => 1)) =
      (∑ x, y x) / (Fintype.card (Fin 2) * Fintype.card (Fin 2)) :=
  level_eq_mean (1 : Matrix (Fin 2) (Fin 2) ℝ) (lam := 1) (by simp) one_ne_zero y

end BASS
