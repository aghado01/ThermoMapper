import Enthymemata.ActiveSubspace.ProjectionMean
import Mathlib.Algebra.Order.Field.Basic
import Mathlib.Data.Real.Basic
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.Rank
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# Active subspaces — rotate-and-refit's non-invariance (AS-4)

The statements of `Protolemmata/BASS/active-subspace-lemmas.md`, AS-4, both readings.

**The frame.** BASS's class in the active variables `u` is a sum of products of one-variable hinge
factors (`IsHingeSum`).
* `ridgeHinge_not_separable`: for `a, b ≠ 0`, `(a u₁ + b u₂ − t)₊` is no finite sum
  `Σₘ φₘ(u₁)ψₘ(u₂)`: on a designed grid it is triangular with diagonal `½` (`ridgeGrid_det`), so of
  full rank, where `M` separable terms have rank at most `M`.
* `hinge_isHingeSum`, `hinge_comp_not_hingeSum`, `hingeSum_not_closed`: a linear map `R` of `u` with
  two nonzero entries in a row carries a hinge out of the class. For orthogonal `R`, only a signed
  permutation keeps it. `hingeFactor_neg`, `hingeFactor_scale` and `IsHingeFactor.monomial` give the
  factor level, and `IsHingeSum.comp_monomial` the class: every monomial map keeps it.

**The scale.**
* `outer_pullback`, `activeMatrix_pullback`: `C_z = LᵀC_xL` for `x = Lz`.
* `pullback_eigen_of_orthogonal`, `conformal_pullback`: an orthogonal `L` carries eigenpairs; `cQ`
  scales the eigenvalues by `c²`.
* `pullback_ker`, `pullback_range`: under every invertible `L`, `ker C_z = L⁻¹ ker C_x` and
  `range C_z = Lᵀ range C_x`, so a rank-`k` `C`'s leading `k` (its range) is carried whatever `L`
  is.
* `leading_moves_of_not_conformal`: every `L` with `LLᵀ` not a multiple of `I` moves the leading
  eigenvector of some positive semidefinite `C`, built on two eigenvectors of `LLᵀ`.
* `leading_moves`, `leading_C_coord`, `leading_pullback_coord`: `C = diag(1, ½)`, `L = diag(1, 2)`:
  `C`'s leading eigenvectors are the multiples of `e₁`, `LᵀCL`'s those of `e₂`, and `Lᵀe₁` is not
  leading for `LᵀCL`. The leading direction moves from `x₁` to `x₂`.
* `hat_invariant`: the g-form's hat `B(BᵀB)⁻¹Bᵀ` is unchanged by `B ↦ BT`, `T` invertible, so its
  marginal reads no column's scale. `rawRidge_moves`: the raw ridge's covariance does.
* `trainingRange_factor`, `trainingRange_hinge`: under the training-range rule the evaluated hinge
  product does not read the declared box; `kappa_moves`, `kappa_ratio`: `κ` does.
-/

open Matrix Finset

namespace ActiveSubspace

/-! ### The frame -/

/-- The designed grid's matrix `(j − i + ½)₊` is upper triangular with diagonal `½`. -/
theorem ridgeGrid_det (N : ℕ) :
    (Matrix.of fun i j : Fin N => max 0 (((j : ℕ) : ℝ) - ((i : ℕ) : ℝ) + 1 / 2)).det =
      (1 / 2) ^ N := by
  rw [det_of_upperTriangular]
  · simp only [of_apply, sub_self, zero_add]
    rw [max_eq_right (by norm_num), Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  · intro i j hij
    have h : ((j : ℕ) : ℝ) + 1 ≤ ((i : ℕ) : ℝ) := by exact_mod_cast hij
    simp only [of_apply]
    exact max_eq_left (by linarith)

/-- **AS-4(i).** For `a, b ≠ 0`, the ridge hinge `(a u₁ + b u₂ − t)₊` is not a finite sum of
products `φₘ(u₁) ψₘ(u₂)`. -/
theorem ridgeHinge_not_separable {a b : ℝ} (ha : a ≠ 0) (hb : b ≠ 0) (t : ℝ) (M : ℕ)
    (φ ψ : Fin M → ℝ → ℝ) : ¬ ∀ x y, max 0 (a * x + b * y - t) = ∑ m, φ m x * ψ m y := by
  intro h
  let G : Matrix (Fin (M + 1)) (Fin (M + 1)) ℝ :=
    Matrix.of fun i j => max 0 (((j : ℕ) : ℝ) - ((i : ℕ) : ℝ) + 1 / 2)
  have hfac : G = (Matrix.of fun (i : Fin (M + 1)) (m : Fin M) => φ m (-((i : ℕ) : ℝ) / a)) *
      (Matrix.of fun (m : Fin M) (j : Fin (M + 1)) => ψ m ((t + ((j : ℕ) : ℝ) + 1 / 2) / b)) := by
    ext i j
    rw [mul_apply]
    simp only [G, of_apply]
    rw [← h]
    congr 1
    field_simp
    ring
  have hrank : G.rank ≤ M := by
    rw [hfac]
    exact (rank_mul_le_left _ _).trans ((rank_le_card_width _).trans (by simp))
  have hdet : G.det ≠ 0 := by
    simp only [G, ridgeGrid_det]
    positivity
  have hfull := rank_of_isUnit G ((isUnit_iff_isUnit_det G).2 (isUnit_iff_ne_zero.2 hdet))
  rw [Fintype.card_fin] at hfull
  omega

/-- A hinge factor: `1` (its variable absent from the term) or `(s (x − t))₊^q` with `s = ±1`. -/
def IsHingeFactor (g : ℝ → ℝ) : Prop :=
  g = (fun _ => 1) ∨ ∃ (s t : ℝ) (q : ℕ), (s = 1 ∨ s = -1) ∧ g = fun x => max 0 (s * (x - t)) ^ q

/-- BASS's class in `d` variables: a finite sum of coefficients times products of hinge factors, one
per variable. The intercept is the term whose factors are all `1`. -/
def IsHingeSum {d : ℕ} (f : (Fin d → ℝ) → ℝ) : Prop :=
  ∃ (M : ℕ) (β : Fin M → ℝ) (g : Fin M → Fin d → ℝ → ℝ), (∀ m l, IsHingeFactor (g m l)) ∧
    ∀ u, f u = ∑ m, β m * ∏ l, g m l (u l)

/-- A single hinge in one variable is in the class. -/
theorem hinge_isHingeSum {d : ℕ} (i : Fin d) (t : ℝ) :
    IsHingeSum (fun u : Fin d → ℝ => max 0 (u i - t)) := by
  refine ⟨1, fun _ => 1, fun _ l => if l = i then (fun x => max 0 (1 * (x - t)) ^ 1)
    else (fun _ => 1), ?_, ?_⟩
  · intro m l
    by_cases hl : l = i
    · exact Or.inr ⟨1, t, 1, Or.inl rfl, by simp [hl]⟩
    · exact Or.inl (by simp [hl])
  · intro u
    simp only [Finset.univ_unique, Finset.sum_singleton, one_mul]
    rw [Finset.prod_eq_single i]
    · simp
    · intro l _ hl
      simp [hl]
    · simp

/-- **AS-4(iii).** If a row of `R` has two nonzero entries, the hinge `(u_i − t)₊` composed with
`R` is not in the class: on the plane of those two variables it is AS-4(i)'s ridge. -/
theorem hinge_comp_not_hingeSum {d : ℕ} (R : Matrix (Fin d) (Fin d) ℝ) {i j k : Fin d}
    (hjk : j ≠ k) (hj : R i j ≠ 0) (hk : R i k ≠ 0) (t : ℝ) :
    ¬ IsHingeSum (fun u => max 0 ((R *ᵥ u) i - t)) := by
  rintro ⟨M, β, g, -, hf⟩
  apply ridgeHinge_not_separable hj hk t M (fun m x => β m * g m j x)
    (fun m y => ∏ l ∈ univ.erase j, g m l (if l = k then y else 0))
  intro x y
  let u : Fin d → ℝ := fun l => if l = j then x else if l = k then y else 0
  have hRu : (R *ᵥ u) i = R i j * x + R i k * y := by
    have hpt : ∀ l, R i l * u l =
        (if l = j then R i l * x else 0) + (if l = k then R i l * y else 0) := by
      intro l
      by_cases hl : l = j
      · subst hl
        simp [u, hjk]
      · by_cases hl' : l = k
        · subst hl'
          simp [u, hl]
        · simp [u, hl, hl']
    simp only [mulVec, dotProduct, hpt, Finset.sum_add_distrib, Finset.sum_ite_eq',
      Finset.mem_univ, if_true]
  have hfu := hf u
  simp only [hRu] at hfu
  rw [hfu]
  refine Finset.sum_congr rfl fun m _ => ?_
  rw [← Finset.mul_prod_erase univ (fun l => g m l (u l)) (Finset.mem_univ j)]
  have huj : u j = x := by simp [u]
  rw [huj, ← mul_assoc]
  congr 1
  refine Finset.prod_congr rfl fun l hl => ?_
  have hlj : l ≠ j := Finset.ne_of_mem_erase hl
  simp [u, hlj]

/-- **AS-4(iii), the class is not closed.** For `R` with two nonzero entries in a row, some member
of the class leaves it under `u ↦ Ru`. -/
theorem hingeSum_not_closed {d : ℕ} (R : Matrix (Fin d) (Fin d) ℝ) {i j k : Fin d} (hjk : j ≠ k)
    (hj : R i j ≠ 0) (hk : R i k ≠ 0) :
    ∃ f : (Fin d → ℝ) → ℝ, IsHingeSum f ∧ ¬ IsHingeSum (fun u => f (R *ᵥ u)) :=
  ⟨fun u => max 0 (u i - 0), hinge_isHingeSum i 0, hinge_comp_not_hingeSum R hjk hj hk 0⟩

/-- A sign flip maps a hinge factor to the factor of the opposite sign at the reflected knot. -/
theorem hingeFactor_neg (s t x : ℝ) (q : ℕ) :
    max 0 (s * (-x - t)) ^ q = max 0 (-s * (x - -t)) ^ q := by
  rw [show s * (-x - t) = -s * (x - -t) by ring]

/-- A positive scale maps a hinge factor to a multiple of the factor at the scaled knot. -/
theorem hingeFactor_scale (s t x c : ℝ) (q : ℕ) (hc : 0 < c) :
    max 0 (s * (c * x - t)) ^ q = c ^ q * max 0 (s * (x - t / c)) ^ q := by
  rw [show s * (c * x - t) = c * (s * (x - t / c)) by field_simp, ← mul_pow,
    mul_max_of_nonneg _ _ hc.le, mul_zero]

/-- A hinge factor read through `x ↦ ε c x`, with `ε = ±1` and `c > 0`, is a positive multiple of a
hinge factor: the sign moves into the factor's sign and the knot, the scale into the multiple. -/
theorem IsHingeFactor.monomial {g : ℝ → ℝ} (hg : IsHingeFactor g) {ε c : ℝ}
    (hε : ε = 1 ∨ ε = -1) (hc : 0 < c) :
    ∃ κ : ℝ, 0 < κ ∧ ∃ g' : ℝ → ℝ, IsHingeFactor g' ∧ ∀ x, g (ε * c * x) = κ * g' x := by
  rcases hg with rfl | ⟨s, t, q, hs, rfl⟩
  · exact ⟨1, one_pos, fun _ => 1, Or.inl rfl, fun x => by simp⟩
  · rcases hε with rfl | rfl
    · refine ⟨c ^ q, pow_pos hc q, fun x => max 0 (s * (x - t / c)) ^ q,
        Or.inr ⟨s, t / c, q, hs, rfl⟩, fun x => ?_⟩
      simp only [one_mul]
      exact hingeFactor_scale s t x c q hc
    · refine ⟨c ^ q, pow_pos hc q, fun x => max 0 (-s * (x - -t / c)) ^ q,
        Or.inr ⟨-s, -t / c, q, by rcases hs with rfl | rfl <;> norm_num, rfl⟩, fun x => ?_⟩
      have h := hingeFactor_scale (-s) (-t) x c q hc
      rw [show -s * (c * x - -t) = s * (-1 * c * x - t) by ring] at h
      exact h

/-- **AS-4(iii), the converse.** The class is kept by every monomial map of the variables — a
permutation with signs and positive scales: the factors' multiples move into the coefficients. So
BASS's class, and with it a projected refit's model, does not see a frame's signs, order or
scale. -/
theorem IsHingeSum.comp_monomial {d : ℕ} {f : (Fin d → ℝ) → ℝ} (hf : IsHingeSum f)
    (σ : Equiv.Perm (Fin d)) (ε c : Fin d → ℝ) (hε : ∀ l, ε l = 1 ∨ ε l = -1)
    (hc : ∀ l, 0 < c l) : IsHingeSum (fun u => f (fun l => ε l * c l * u (σ l))) := by
  obtain ⟨M, β, g, hg, hfg⟩ := hf
  choose κ _hκ g' hg' hgg' using fun m l => (hg m l).monomial (hε l) (hc l)
  refine ⟨M, fun m => β m * ∏ l, κ m l, fun m k => g' m (σ.symm k),
    fun m k => hg' m (σ.symm k), fun u => ?_⟩
  dsimp only
  rw [hfg]
  refine Finset.sum_congr rfl fun m _ => ?_
  simp only [hgg', Finset.prod_mul_distrib]
  have hre : ∏ l, g' m l (u (σ l)) = ∏ k, g' m (σ.symm k) (u k) := by
    rw [← Equiv.prod_comp σ (fun k => g' m (σ.symm k) (u k))]
    simp only [Equiv.symm_apply_apply]
  rw [hre]
  ring

/-! ### The scale -/

section Scale

variable {n : Type*} [Fintype n]

/-- **AS-4(v).** For `x = Lz`, the gradient `∇_z = Lᵀ∇_x` has outer product `Lᵀ(∇_x∇_xᵀ)L`. -/
theorem outer_pullback (L : Matrix n n ℝ) (g : n → ℝ) :
    vecMulVec (Lᵀ *ᵥ g) (Lᵀ *ᵥ g) = Lᵀ * vecMulVec g g * L := by
  ext a b
  simp only [vecMulVec_apply, mulVec, dotProduct, mul_apply, transpose_apply, Finset.sum_mul,
    Finset.mul_sum]
  exact Finset.sum_congr rfl fun _ _ => Finset.sum_congr rfl fun _ _ => by ring

/-- **AS-4(v).** So `C_z = LᵀC_xL` for any weighted average of gradients. -/
theorem activeMatrix_pullback {Ω : Type*} [Fintype Ω] (p : Ω → ℝ) (L : Matrix n n ℝ)
    (g : Ω → n → ℝ) :
    ∑ ω, p ω • vecMulVec (Lᵀ *ᵥ g ω) (Lᵀ *ᵥ g ω) = Lᵀ * (∑ ω, p ω • vecMulVec (g ω) (g ω)) * L := by
  simp only [outer_pullback, Matrix.mul_sum, Matrix.sum_mul, Matrix.mul_smul, Matrix.smul_mul]

variable [DecidableEq n]

/-- **AS-4(vi).** An orthogonal `Q` carries eigenpairs: `Cv = μv` exactly when
`(QᵀCQ)(Qᵀv) = μ(Qᵀv)`. -/
theorem pullback_eigen_of_orthogonal (Q C : Matrix n n ℝ) (hQ : Q * Qᵀ = 1) (μ : ℝ) (v : n → ℝ) :
    C *ᵥ v = μ • v ↔ (Qᵀ * C * Q) *ᵥ (Qᵀ *ᵥ v) = μ • (Qᵀ *ᵥ v) := by
  have key : (Qᵀ * C * Q) *ᵥ (Qᵀ *ᵥ v) = Qᵀ *ᵥ (C *ᵥ v) := by
    rw [mulVec_mulVec, mulVec_mulVec, Matrix.mul_assoc, hQ, Matrix.mul_one]
  rw [key]
  constructor
  · intro h
    rw [h, mulVec_smul]
  · intro h
    have h2 := congrArg (Q *ᵥ ·) h
    simp only [mulVec_smul, mulVec_mulVec, ← Matrix.mul_assoc, hQ, Matrix.one_mul,
      one_mulVec] at h2
    exact h2

omit [DecidableEq n] in
/-- **AS-4(vi).** A conformal `cQ` scales the pulled-back matrix by `c²`. -/
theorem conformal_pullback (Q C : Matrix n n ℝ) (c : ℝ) :
    (c • Q)ᵀ * C * (c • Q) = c ^ 2 • (Qᵀ * C * Q) := by
  rw [transpose_smul, Matrix.smul_mul, Matrix.smul_mul, Matrix.mul_smul, smul_smul, sq]

/-- **AS-4(vii).** Under every invertible `L`, `ker(LᵀCL) = L⁻¹ ker C`. -/
theorem pullback_ker (L C : Matrix n n ℝ) (hL : IsUnit L.det) (v : n → ℝ) :
    (Lᵀ * C * L) *ᵥ v = 0 ↔ C *ᵥ (L *ᵥ v) = 0 := by
  have key : (Lᵀ * C * L) *ᵥ v = Lᵀ *ᵥ (C *ᵥ (L *ᵥ v)) := by
    rw [mulVec_mulVec, mulVec_mulVec]
  rw [key]
  constructor
  · intro h
    have hLt : IsUnit Lᵀ.det := by rwa [det_transpose]
    have h2 := congrArg ((Lᵀ)⁻¹ *ᵥ ·) h
    rw [mulVec_mulVec (M := (Lᵀ)⁻¹) (N := Lᵀ), nonsing_inv_mul _ hLt, one_mulVec,
      mulVec_zero] at h2
    exact h2
  · intro h
    rw [h, mulVec_zero]

/-- **AS-4(vii).** Under every invertible `L`, `range(LᵀCL) = Lᵀ range C`: the exact active
subspace of a rank-`k` `C`, its range, is carried whatever the metric. -/
theorem pullback_range (L C : Matrix n n ℝ) (hL : IsUnit L.det) (w : n → ℝ) :
    (∃ v, (Lᵀ * C * L) *ᵥ v = w) ↔ ∃ u, Lᵀ *ᵥ (C *ᵥ u) = w := by
  constructor
  · rintro ⟨v, rfl⟩
    exact ⟨L *ᵥ v, by rw [mulVec_mulVec, mulVec_mulVec]⟩
  · rintro ⟨u, rfl⟩
    refine ⟨L⁻¹ *ᵥ u, ?_⟩
    rw [mulVec_mulVec, Matrix.mul_assoc (Lᵀ * C), mul_nonsing_inv _ hL, Matrix.mul_one,
      ← mulVec_mulVec]

/-- `v` is a leading eigenvector of `C`: an eigenvector whose eigenvalue bounds every Rayleigh
quotient of `C`. -/
def IsLeading (C : Matrix n n ℝ) (v : n → ℝ) : Prop :=
  v ≠ 0 ∧ ∃ μ : ℝ, C *ᵥ v = μ • v ∧ ∀ w, w ⬝ᵥ (C *ᵥ w) ≤ μ * (w ⬝ᵥ w)

end Scale

section General

variable {n : Type*} [Fintype n] [DecidableEq n]

omit [DecidableEq n] in
/-- A rank-one matrix applied to a vector. -/
theorem outer_mulVec (a b w : n → ℝ) : vecMulVec a b *ᵥ w = (b ⬝ᵥ w) • a := by
  ext k
  simp only [mulVec, dotProduct, vecMulVec_apply, Pi.smul_apply, smul_eq_mul, Finset.sum_mul]
  exact Finset.sum_congr rfl fun _ _ => by ring

omit [DecidableEq n] in
theorem dot_self_nonneg (w : n → ℝ) : 0 ≤ w ⬝ᵥ w :=
  Finset.sum_nonneg fun i _ => mul_self_nonneg (w i)

omit [DecidableEq n] in
/-- Bessel's inequality for an orthonormal pair. -/
theorem bessel_pair (u₁ u₂ w : n → ℝ) (h11 : u₁ ⬝ᵥ u₁ = 1) (h22 : u₂ ⬝ᵥ u₂ = 1)
    (h12 : u₁ ⬝ᵥ u₂ = 0) : (u₁ ⬝ᵥ w) ^ 2 + (u₂ ⬝ᵥ w) ^ 2 ≤ w ⬝ᵥ w := by
  have h := dot_self_nonneg (w - (u₁ ⬝ᵥ w) • u₁ - (u₂ ⬝ᵥ w) • u₂)
  have e1 : w ⬝ᵥ u₁ = u₁ ⬝ᵥ w := dotProduct_comm _ _
  have e2 : w ⬝ᵥ u₂ = u₂ ⬝ᵥ w := dotProduct_comm _ _
  have e3 : u₂ ⬝ᵥ u₁ = 0 := by rw [dotProduct_comm, h12]
  simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, smul_eq_mul, e1,
    e2, e3, h11, h22, h12] at h
  nlinarith [h]

/-- The eigenvectors of a symmetric matrix are orthonormal as plain vectors. -/
theorem eigvec_dot {A : Matrix n n ℝ} (hA : A.IsHermitian) (i j : n) :
    (⇑(hA.eigenvectorBasis i) : n → ℝ) ⬝ᵥ ⇑(hA.eigenvectorBasis j) = if i = j then 1 else 0 := by
  have h := orthonormal_iff_ite.1 hA.eigenvectorBasis.orthonormal i j
  rw [EuclideanSpace.inner_eq_star_dotProduct] at h
  simpa [dotProduct_comm] using h

/-- **AS-4(viii), in general.** If `LLᵀ` is not a multiple of `I` (`L` is not conformal), some
positive semidefinite `C` has a leading eigenvector `v` whose image `Lᵀv` does not lead `LᵀCL`:
the leading direction moves under every non-conformal change of units. -/
theorem leading_moves_of_not_conformal (L : Matrix n n ℝ)
    (hL : ∀ c : ℝ, L * Lᵀ ≠ c • (1 : Matrix n n ℝ)) :
    ∃ (C : Matrix n n ℝ) (v : n → ℝ), C.PosSemidef ∧ IsLeading C v ∧
      ¬ IsLeading (Lᵀ * C * L) (Lᵀ *ᵥ v) := by
  have hpsd : (L * Lᵀ).PosSemidef := by
    have := posSemidef_self_mul_conjTranspose L
    rwa [conjTranspose_eq_transpose_of_trivial] at this
  have hS : (L * Lᵀ).IsHermitian := hpsd.1
  set μ := hS.eigenvalues with hμdef
  set u : n → n → ℝ := fun k => ⇑(hS.eigenvectorBasis k) with hudef
  have hu : ∀ k, (L * Lᵀ) *ᵥ u k = μ k • u k := fun k => hS.mulVec_eigenvectorBasis k
  have hdot : ∀ k l, u k ⬝ᵥ u l = if k = l then 1 else 0 := fun k l => eigvec_dot hS k l
  -- two distinct eigenvalues
  obtain ⟨i, j, hij⟩ : ∃ i j, μ j < μ i := by
    by_contra hne
    push Not at hne
    rcases isEmpty_or_nonempty n with hn | ⟨⟨k⟩⟩
    · exact hL 0 (by ext a; exact isEmptyElim a)
    · apply hL (μ k)
      have hd : diagonal μ = μ k • (1 : Matrix n n ℝ) := by
        ext a b
        by_cases hab : a = b
        · subst hab
          simp [le_antisymm (hne a k) (hne k a)]
        · simp [hab]
      conv_lhs => rw [eq_eigU hS]
      rw [← hμdef, hd, Matrix.mul_smul, Matrix.mul_one, Matrix.smul_mul, eigU_mul_transpose]
  have hji : j ≠ i := fun h => (lt_irrefl _ (h ▸ hij))
  have hμj : 0 ≤ μ j := hpsd.eigenvalues_nonneg j
  have hμi : 0 < μ i := lt_of_le_of_lt hμj hij
  have hii : u i ⬝ᵥ u i = 1 := by simp [hdot]
  have hjj : u j ⬝ᵥ u j = 1 := by simp [hdot]
  have hij0 : u i ⬝ᵥ u j = 0 := by simp [hdot, hji.symm]
  have hji0 : u j ⬝ᵥ u i = 0 := by rw [dotProduct_comm, hij0]
  set ci : ℝ := (μ i + μ j) / (2 * μ i) with hci
  have hci0 : 0 < ci := by positivity
  have hci1 : ci < 1 := by rw [hci, div_lt_one (by positivity)]; linarith
  have hcimu : μ j < μ i * ci := by
    have h : μ i * ci = (μ i + μ j) / 2 := by
      rw [hci, mul_div_assoc', div_eq_iff (by positivity)]
      ring
    rw [h]
    linarith
  set C : Matrix n n ℝ := ci • vecMulVec (u i) (u i) + vecMulVec (u j) (u j) with hC
  have hCv : ∀ w, C *ᵥ w = (ci * (u i ⬝ᵥ w)) • u i + (u j ⬝ᵥ w) • u j := by
    intro w
    rw [hC, add_mulVec, smul_mulVec, outer_mulVec, outer_mulVec, smul_smul]
  have hray : ∀ w, w ⬝ᵥ (C *ᵥ w) = ci * (u i ⬝ᵥ w) ^ 2 + (u j ⬝ᵥ w) ^ 2 := by
    intro w
    rw [hCv, dotProduct_add, dotProduct_smul, dotProduct_smul, smul_eq_mul, smul_eq_mul,
      dotProduct_comm w (u i), dotProduct_comm w (u j)]
    ring
  refine ⟨C, u j, ?_, ⟨?_, 1, ?_, fun w => ?_⟩, ?_⟩
  · have h1 := (posSemidef_vecMulVec_self_star (u i)).smul hci0.le
    have h2 := posSemidef_vecMulVec_self_star (u j)
    simp only [star_trivial] at h1 h2
    exact h1.add h2
  · intro h0
    rw [h0, dotProduct_zero] at hjj
    exact zero_ne_one hjj
  · rw [hCv, hij0, hjj, mul_zero, zero_smul, zero_add, one_smul]
  · rw [hray, one_mul]
    have hb := bessel_pair (u i) (u j) w hii hjj hij0
    nlinarith [sq_nonneg (u i ⬝ᵥ w)]
  · rintro ⟨hne, m, heig, hrayL⟩
    -- the eigenvalue at Lᵀuⱼ is μⱼ
    have hLv : (Lᵀ * C * L) *ᵥ (Lᵀ *ᵥ u j) = μ j • (Lᵀ *ᵥ u j) := by
      have hre : (Lᵀ * C * L) *ᵥ (Lᵀ *ᵥ u j) = Lᵀ *ᵥ (C *ᵥ ((L * Lᵀ) *ᵥ u j)) := by
        simp only [mulVec_mulVec, Matrix.mul_assoc]
      rw [hre, hu j, mulVec_smul, hCv, hij0, hjj, mul_zero, zero_smul, zero_add, one_smul,
        mulVec_smul]
    rw [hLv] at heig
    obtain ⟨k, hk⟩ := Function.ne_iff.1 hne
    have hm : m = μ j := by
      have := congrFun heig k
      simp only [Pi.smul_apply, smul_eq_mul, Pi.zero_apply] at this hk
      exact (mul_right_cancel₀ hk this).symm
    -- the Rayleigh quotient at Lᵀuᵢ is μᵢ cᵢ, above μⱼ
    have hq := hrayL (Lᵀ *ᵥ u i)
    have hLw : L *ᵥ (Lᵀ *ᵥ u i) = μ i • u i := by rw [mulVec_mulVec, hu i]
    have hlhs : (Lᵀ *ᵥ u i) ⬝ᵥ ((Lᵀ * C * L) *ᵥ (Lᵀ *ᵥ u i)) = μ i ^ 2 * ci := by
      rw [← mulVec_mulVec, ← mulVec_mulVec, hLw, dotProduct_mulVec, vecMul_transpose, hLw,
        mulVec_smul, dotProduct_smul, smul_dotProduct, hray, hii, hji0]
      simp only [smul_eq_mul]
      ring
    have hrhs : (Lᵀ *ᵥ u i) ⬝ᵥ (Lᵀ *ᵥ u i) = μ i := by
      rw [dotProduct_mulVec, vecMul_transpose, hLw, smul_dotProduct, hii, smul_eq_mul, mul_one]
    rw [hlhs, hrhs, hm] at hq
    nlinarith [hcimu, hμi]

end General

section Example

/-- `C = diag(1, ½)`. -/
noncomputable def Cex : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, 1 / 2]

/-- `L = diag(1, 2)`, a box twice as wide on `x₂`. -/
def Lex : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, 2]

theorem pullback_ex : Lexᵀ * Cex * Lex = !![1, 0; 0, 2] := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [Lex, Cex, mul_apply, Fin.sum_univ_two]

/-- `C`'s leading eigenvectors are the multiples of `e₁`. -/
theorem leading_C_coord (v : Fin 2 → ℝ) (hv : IsLeading Cex v) : v 1 = 0 := by
  obtain ⟨-, μ, heig, hray⟩ := hv
  have h1 : 1 ≤ μ := by simpa [Cex, mulVec, dotProduct, Fin.sum_univ_two] using hray ![1, 0]
  have h2 : 1 / 2 * v 1 = μ * v 1 := by
    simpa [Cex, mulVec, dotProduct, Fin.sum_univ_two] using congrFun heig 1
  have h3 : (μ - 1 / 2) * v 1 = 0 := by linarith
  rcases mul_eq_zero.1 h3 with h | h
  · exfalso
    linarith
  · exact h

/-- `LᵀCL`'s leading eigenvectors are the multiples of `e₂`. -/
theorem leading_pullback_coord (v : Fin 2 → ℝ) (hv : IsLeading (Lexᵀ * Cex * Lex) v) : v 0 = 0 := by
  rw [pullback_ex] at hv
  obtain ⟨-, μ, heig, hray⟩ := hv
  have h1 := hray ![0, 1]
  have h2 := congrFun heig 0
  simp [mulVec, dotProduct, Fin.sum_univ_two] at h1 h2
  nlinarith

/-- **AS-4(viii).** The leading direction moves. `e₁` is a leading eigenvector of `C = diag(1, ½)`;
carried by `Lᵀ` for `L = diag(1, 2)` it is not leading for `LᵀCL`, whose leading eigenvector is
`e₂`. So the leading direction is `x₁` in one metric and `x₂` in the other. -/
theorem leading_moves :
    IsLeading Cex ![1, 0] ∧ ¬ IsLeading (Lexᵀ * Cex * Lex) (Lexᵀ *ᵥ ![1, 0]) ∧
      IsLeading (Lexᵀ * Cex * Lex) ![0, 1] := by
  refine ⟨⟨fun h => by simpa using congrFun h 0, 1, ?_, fun w => ?_⟩, fun h => ?_,
    ⟨fun h => by simpa using congrFun h 1, 2, ?_, fun w => ?_⟩⟩
  · ext i
    fin_cases i <;> simp [Cex, mulVec, dotProduct, Fin.sum_univ_two]
  · simp [Cex, mulVec, dotProduct, Fin.sum_univ_two]
    nlinarith [sq_nonneg (w 1)]
  · have h0 := leading_pullback_coord _ h
    simp [Lex, mulVec, dotProduct, Fin.sum_univ_two] at h0
  · rw [pullback_ex]
    ext i
    fin_cases i <;> simp [mulVec, dotProduct, Fin.sum_univ_two]
  · rw [pullback_ex]
    simp [mulVec, dotProduct, Fin.sum_univ_two]
    nlinarith [sq_nonneg (w 0)]

end Example

/-! ### The projected fit's rules -/

/-- **AS-4(ix).** The g-form's hat is unchanged by `B ↦ BT` for invertible `T` (a column scale is a
diagonal `T`), so its `R²` and collapsed marginal read no column's scale. -/
theorem hat_invariant {N k : Type*} [Fintype N] [Fintype k] [DecidableEq k] (B : Matrix N k ℝ)
    (T : Matrix k k ℝ) (hT : IsUnit T.det) :
    (B * T) * ((B * T)ᵀ * (B * T))⁻¹ * (B * T)ᵀ = B * (Bᵀ * B)⁻¹ * Bᵀ := by
  have hTt : IsUnit Tᵀ.det := by rwa [det_transpose]
  rw [transpose_mul, show Tᵀ * Bᵀ * (B * T) = Tᵀ * (Bᵀ * B) * T by
      simp only [Matrix.mul_assoc], Matrix.mul_inv_rev, Matrix.mul_inv_rev]
  simp only [Matrix.mul_assoc]
  rw [mul_nonsing_inv_cancel_left _ _ hT, nonsing_inv_mul_cancel_left _ _ hTt]

/-- **AS-4(ix).** The raw ridge's marginal covariance `σ²I + τ bbᵀ` moves when a column is scaled
by `c` with `c² ≠ 1`. -/
theorem rawRidge_moves {N : Type*} [DecidableEq N] (b : N → ℝ) (i : N) (hb : b i ≠ 0)
    (σ2 τ c : ℝ) (hτ : τ ≠ 0) (hc : c ^ 2 ≠ 1) :
    σ2 • (1 : Matrix N N ℝ) + τ • vecMulVec (c • b) (c • b) ≠ σ2 • 1 + τ • vecMulVec b b := by
  intro h
  have hii := congrFun (congrFun h i) i
  simp only [Matrix.add_apply, Matrix.smul_apply, Matrix.one_apply_eq, vecMulVec_apply,
    Pi.smul_apply, smul_eq_mul] at hii
  have h2 : τ * b i ^ 2 * (c ^ 2 - 1) = 0 := by linear_combination hii
  rcases mul_eq_zero.1 h2 with h3 | h3
  · exact mul_ne_zero hτ (pow_ne_zero 2 hb) h3
  · exact hc (by linarith)

/-- **AS-4(x).** A factor read in a box of width `W`, `((z/W)₊)^q`, normalized by its peak in the
same units, `(P/W)^q`, does not read `W`. -/
theorem trainingRange_factor (z P W : ℝ) (q : ℕ) (hW : 0 < W) :
    max 0 (z / W) ^ q / (P / W) ^ q = max 0 z ^ q / P ^ q := by
  rw [show max 0 (z / W) = max 0 z / W by rw [← max_div_div_right hW.le, zero_div], div_pow,
    div_pow, div_div_div_cancel_right₀ (pow_ne_zero q hW.ne')]

/-- **AS-4(x).** So a term's evaluated hinge product, under the training-range rule (`P` the
training-range peak `U − t` or `t − L`, `z = ±(x − t)`), is the same under any declared box. -/
theorem trainingRange_hinge {ι : Type*} (s : Finset ι) (z P W : ι → ℝ) (q : ℕ)
    (hW : ∀ l, 0 < W l) :
    ∏ l ∈ s, max 0 (z l / W l) ^ q / (P l / W l) ^ q = ∏ l ∈ s, max 0 (z l) ^ q / P l ^ q :=
  Finset.prod_congr rfl fun l _ => trainingRange_factor _ _ _ q (hW l)

/-- **AS-4(x).** `κ`, the peak in the box's units, does move with the box. -/
theorem kappa_moves (P W W' : ℝ) (q : ℕ) (hq : q ≠ 0) (hP : 0 < P) (hW : 0 < W) (hW' : 0 < W')
    (h : W ≠ W') : (P / W) ^ q ≠ (P / W') ^ q := by
  intro h'
  rw [pow_left_inj₀ (by positivity) (by positivity) hq, div_eq_div_iff hW.ne' hW'.ne'] at h'
  exact h (mul_left_cancel₀ hP.ne' h').symm

/-- **AS-4(x).** By the factor `(W′/W)^q` per factor. -/
theorem kappa_ratio (P W W' : ℝ) (q : ℕ) (hP : P ≠ 0) (hW : W ≠ 0) (hW' : W' ≠ 0) :
    (P / W) ^ q / (P / W') ^ q = (W' / W) ^ q := by
  rw [← div_pow]
  congr 1
  field_simp

end ActiveSubspace
