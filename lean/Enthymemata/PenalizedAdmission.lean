import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Algebra.Order.Star.Real
import Mathlib.Algebra.QuadraticDiscriminant
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.Analysis.Real.Sqrt
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv

/-!
# A penalized system's admission (PA-5 to PA-7)

BAPS factors `A(c) = Σₖ cₖ Sₖ` at every sampler step, every `Sₖ` positive semidefinite and fixed
by the design. A gate admits a factor that is positive definite with `κ₁` of its unit-diagonal
scaling `E(c)` at most `10¹⁴`. The statements, checks and prior art are in
`Protolemmata/BAPS/penalized-admission-lemmas.md`.

* **PA-5, identification.** For `c > 0`, `ker A(c) = ⋂ₖ ker Sₖ` (`comb_mulVec_eq_zero_iff`), so
  `A(c)` is positive definite at one `c > 0` exactly when at every one (`comb_posDef_iff_of_pos`).
  BAPS's instance: `penalized_det_eq_zero_iff`.
* **PA-6, the bound.** With `m·c′ₖ ≤ cₖ ≤ M·c′ₖ`: the form (`form_comb_bound`), the extreme
  Rayleigh quotients (`lamMin_comb_ge`, `lamMax_comb_le`), and `κ₂(A(c)) ≤ (M/m)·κ₂(A(c′))`
  (`cond2_comb_le`). The smallest instance attains `M/m` (`cond2_diagPair`).
* **PA-7, the chain to the gate.** The diagonal (`diag_comb_bound`), the scaled matrix's least
  Rayleigh quotient (`lamMin_unitScale_ge`), the norm lines (`norm1_le_card`,
  `norm1_le_sqrt_card`, `mulVec_sq_le_norm1`), and
  `κ₁(E(c)) ≤ q·((1 + √n)/2)·(M/m)·‖E(c′)⁻¹‖₁` (`cond1_unitScale_le_sharp`), whose column step
  `‖P‖₁ ≤ λ_max(P)·(1 + √n)/2` (`norm1_le_half_one_add_sqrt`) is attained (`norm1_columnWitness`);
  the first chain's `q·√n` (`cond1_unitScale_le`) follows.

`λ_min` and `λ_max` are the infimum and supremum of the Rayleigh quotient (`lamMin`, `lamMax`), and
the results above use no spectral theorem. For a symmetric matrix they are its least and greatest
eigenvalues (`lamMin_eq_iInf_eigenvalues`, `lamMax_eq_iSup_eigenvalues`), through Mathlib's
`Matrix.IsHermitian.spectral_theorem`. `‖·‖₂` enters only through bounds `‖Ax‖ ≤ β‖x‖`, and no
operator norm is defined. A coefficient vector's ratio to the reference enters as a bracket
`m ≤ ρₖ ≤ M`; the least and greatest `ρₖ` are its tightest choice.
-/

open Matrix Finset

namespace PenalizedAdmission

variable {n ι : Type*} [Fintype n]

/-! ### Forms -/

section Forms

variable {A : Matrix n n ℝ}

theorem dot_self_nonneg (x : n → ℝ) : 0 ≤ x ⬝ᵥ x :=
  Finset.sum_nonneg fun i _ => mul_self_nonneg (x i)

theorem dot_self_pos {x : n → ℝ} (hx : x ≠ 0) : 0 < x ⬝ᵥ x :=
  lt_of_le_of_ne (dot_self_nonneg x) fun h => hx (dotProduct_self_eq_zero.1 h.symm)

theorem form_symm (hs : Aᵀ = A) (x y : n → ℝ) : y ⬝ᵥ A *ᵥ x = x ⬝ᵥ A *ᵥ y := by
  rw [dotProduct_mulVec, dotProduct_comm, ← mulVec_transpose, hs]

/-- **Cauchy–Schwarz for a positive semidefinite form.** -/
theorem form_sq_le (hs : Aᵀ = A) (hp : ∀ x, 0 ≤ x ⬝ᵥ A *ᵥ x) (x y : n → ℝ) :
    (y ⬝ᵥ A *ᵥ x) ^ 2 ≤ (y ⬝ᵥ A *ᵥ y) * (x ⬝ᵥ A *ᵥ x) := by
  have h : ∀ t : ℝ, 0 ≤ (y ⬝ᵥ A *ᵥ y) * (t * t) + (2 * (y ⬝ᵥ A *ᵥ x)) * t + (x ⬝ᵥ A *ᵥ x) := by
    intro t
    have := hp (x + t • y)
    simp only [mulVec_add, mulVec_smul, add_dotProduct, dotProduct_add, smul_dotProduct,
      dotProduct_smul, smul_eq_mul] at this
    rw [← form_symm hs x y] at this
    nlinarith [this]
  have := discrim_le_zero h
  unfold discrim at this
  nlinarith [this]

/-- A positive semidefinite form vanishes at `x` only where `Ax = 0`. -/
theorem mulVec_eq_zero_of_form (hs : Aᵀ = A) (hp : ∀ x, 0 ≤ x ⬝ᵥ A *ᵥ x) {x : n → ℝ}
    (hx : x ⬝ᵥ A *ᵥ x = 0) : A *ᵥ x = 0 := by
  have h := form_sq_le hs hp x (A *ᵥ x)
  rw [hx, mul_zero] at h
  exact dotProduct_self_eq_zero.1 (pow_eq_zero_iff two_ne_zero |>.1 (le_antisymm h (sq_nonneg _)))

theorem form_nonneg {M : Matrix n n ℝ} (h : M.PosSemidef) (x : n → ℝ) : 0 ≤ x ⬝ᵥ M *ᵥ x := by
  simpa using h.dotProduct_mulVec_nonneg x

omit [Fintype n] in
theorem transpose_eq {M : Matrix n n ℝ} (h : M.IsHermitian) : Mᵀ = M := by
  rw [← conjTranspose_eq_transpose_of_trivial]; exact h.eq

theorem abs_form_le (A : Matrix n n ℝ) (x : n → ℝ) :
    |x ⬝ᵥ A *ᵥ x| ≤ (∑ i, ∑ j, |A i j|) * (x ⬝ᵥ x) := by
  have hxx : ∀ i, x i ^ 2 ≤ x ⬝ᵥ x := fun i => by
    have := Finset.single_le_sum (f := fun j => x j * x j) (fun j _ => mul_self_nonneg (x j))
      (Finset.mem_univ i)
    simpa [sq, dotProduct] using this
  calc |x ⬝ᵥ A *ᵥ x| = |∑ i, ∑ j, x i * A i j * x j| := by
        simp [dotProduct, mulVec, Finset.mul_sum, mul_assoc]
    _ ≤ ∑ i, ∑ j, |x i * A i j * x j| :=
        (abs_sum_le_sum_abs _ _).trans (sum_le_sum fun i _ => abs_sum_le_sum_abs _ _)
    _ ≤ ∑ i, ∑ j, |A i j| * (x ⬝ᵥ x) := by
        gcongr with i _ j _
        rw [abs_mul, abs_mul]
        have : |x i| * |x j| ≤ x ⬝ᵥ x := by
          nlinarith [hxx i, hxx j, sq_nonneg (|x i| - |x j|), sq_abs (x i), sq_abs (x j)]
        nlinarith [abs_nonneg (A i j), abs_nonneg (x i), abs_nonneg (x j)]
    _ = _ := by simp [Finset.sum_mul]

end Forms

/-! ### Rayleigh bounds -/

section Rayleigh

/-- The least Rayleigh quotient, `⨅ x ≠ 0, xᵀAx / xᵀx`: `λ_min` for a symmetric `A`. -/
noncomputable def lamMin (A : Matrix n n ℝ) : ℝ :=
  ⨅ x : {x : n → ℝ // x ≠ 0}, x.1 ⬝ᵥ A *ᵥ x.1 / (x.1 ⬝ᵥ x.1)

/-- The greatest Rayleigh quotient, `⨆ x ≠ 0, xᵀAx / xᵀx`: `λ_max` for a symmetric `A`. -/
noncomputable def lamMax (A : Matrix n n ℝ) : ℝ :=
  ⨆ x : {x : n → ℝ // x ≠ 0}, x.1 ⬝ᵥ A *ᵥ x.1 / (x.1 ⬝ᵥ x.1)

/-- `κ₂ = λ_max / λ_min`. -/
noncomputable def cond2 (A : Matrix n n ℝ) : ℝ := lamMax A / lamMin A

variable {A : Matrix n n ℝ}

theorem bddBelow_rayleigh (A : Matrix n n ℝ) :
    BddBelow (Set.range fun x : {x : n → ℝ // x ≠ 0} => x.1 ⬝ᵥ A *ᵥ x.1 / (x.1 ⬝ᵥ x.1)) := by
  refine ⟨-(∑ i, ∑ j, |A i j|), ?_⟩
  rintro _ ⟨x, rfl⟩
  rw [le_div_iff₀ (dot_self_pos x.2)]
  have := abs_form_le A x.1
  linarith [neg_abs_le (x.1 ⬝ᵥ A *ᵥ x.1)]

theorem bddAbove_rayleigh (A : Matrix n n ℝ) :
    BddAbove (Set.range fun x : {x : n → ℝ // x ≠ 0} => x.1 ⬝ᵥ A *ᵥ x.1 / (x.1 ⬝ᵥ x.1)) := by
  refine ⟨∑ i, ∑ j, |A i j|, ?_⟩
  rintro _ ⟨x, rfl⟩
  rw [div_le_iff₀ (dot_self_pos x.2)]
  have := abs_form_le A x.1
  linarith [le_abs_self (x.1 ⬝ᵥ A *ᵥ x.1)]

theorem lamMin_mul_le (A : Matrix n n ℝ) (x : n → ℝ) : lamMin A * (x ⬝ᵥ x) ≤ x ⬝ᵥ A *ᵥ x := by
  by_cases hx : x = 0
  · subst hx; simp
  · have : lamMin A ≤ x ⬝ᵥ A *ᵥ x / (x ⬝ᵥ x) := ciInf_le (bddBelow_rayleigh A) ⟨x, hx⟩
    exact (le_div_iff₀ (dot_self_pos hx)).1 this

theorem form_le_lamMax_mul (A : Matrix n n ℝ) (x : n → ℝ) : x ⬝ᵥ A *ᵥ x ≤ lamMax A * (x ⬝ᵥ x) := by
  by_cases hx : x = 0
  · subst hx; simp
  · have : x ⬝ᵥ A *ᵥ x / (x ⬝ᵥ x) ≤ lamMax A := le_ciSup (bddAbove_rayleigh A) ⟨x, hx⟩
    exact (div_le_iff₀ (dot_self_pos hx)).1 this

omit [Fintype n] in
theorem nonempty_ne_zero [Nonempty n] : Nonempty {x : n → ℝ // x ≠ 0} :=
  ⟨⟨fun _ => 1, fun h => by simpa using congrFun h (Classical.arbitrary n)⟩⟩

theorem le_lamMin [Nonempty n] {μ : ℝ} (h : ∀ x, μ * (x ⬝ᵥ x) ≤ x ⬝ᵥ A *ᵥ x) : μ ≤ lamMin A :=
  haveI := nonempty_ne_zero (n := n)
  le_ciInf fun x => (le_div_iff₀ (dot_self_pos x.2)).2 (h x.1)

theorem lamMax_le [Nonempty n] {μ : ℝ} (h : ∀ x, x ⬝ᵥ A *ᵥ x ≤ μ * (x ⬝ᵥ x)) : lamMax A ≤ μ :=
  haveI := nonempty_ne_zero (n := n)
  ciSup_le fun x => (div_le_iff₀ (dot_self_pos x.2)).2 (h x.1)

theorem lamMin_nonneg [Nonempty n] (hp : ∀ x, 0 ≤ x ⬝ᵥ A *ᵥ x) : 0 ≤ lamMin A :=
  le_lamMin fun x => by simpa using hp x

theorem lamMax_nonneg [Nonempty n] (hp : ∀ x, 0 ≤ x ⬝ᵥ A *ᵥ x) : 0 ≤ lamMax A := by
  have h1 := form_le_lamMax_mul A (fun _ => 1)
  have h2 := hp (fun _ => 1)
  have h3 : (0 : ℝ) < (fun _ : n => (1 : ℝ)) ⬝ᵥ (fun _ => 1) :=
    dot_self_pos fun h => by simpa using congrFun h (Classical.arbitrary n)
  nlinarith

/-- A symmetric matrix with a positive least Rayleigh quotient is positive definite. -/
theorem posDef_of_lamMin_pos (hs : Aᵀ = A) (h : 0 < lamMin A) : A.PosDef := by
  refine posDef_iff_dotProduct_mulVec.2 ⟨?_, fun x hx => ?_⟩
  · rw [IsHermitian, conjTranspose_eq_transpose_of_trivial, hs]
  · simpa using lt_of_lt_of_le (mul_pos h (dot_self_pos hx)) (lamMin_mul_le A x)

end Rayleigh

/-! ### Inverses -/

section Inverse

variable [DecidableEq n] {A : Matrix n n ℝ}

theorem inv_transpose_eq (hs : Aᵀ = A) : (A⁻¹)ᵀ = A⁻¹ := by rw [transpose_nonsing_inv, hs]

theorem form_inv_nonneg (hu : IsUnit A.det) (hp : ∀ x, 0 ≤ x ⬝ᵥ A *ᵥ x) (x : n → ℝ) :
    0 ≤ x ⬝ᵥ A⁻¹ *ᵥ x := by
  have := hp (A⁻¹ *ᵥ x)
  rwa [mulVec_mulVec, mul_nonsing_inv _ hu, one_mulVec, dotProduct_comm] at this

omit [DecidableEq n] in
/-- `‖Py‖² ≤ β · yᵀPy` for a positive semidefinite `P` whose form is at most `β · ‖·‖²`. -/
theorem sq_norm_le_form {P : Matrix n n ℝ} (hs : Pᵀ = P) (hp : ∀ x, 0 ≤ x ⬝ᵥ P *ᵥ x) {β : ℝ}
    (hβ0 : 0 ≤ β) (hβ : ∀ z, z ⬝ᵥ P *ᵥ z ≤ β * (z ⬝ᵥ z)) (y : n → ℝ) :
    (P *ᵥ y) ⬝ᵥ (P *ᵥ y) ≤ β * (y ⬝ᵥ P *ᵥ y) := by
  have hcs := form_sq_le hs hp y (P *ᵥ y)
  have h2 := hβ (P *ᵥ y)
  have hy := hp y
  rcases (dot_self_nonneg (P *ᵥ y)).eq_or_lt with h0 | hpos
  · rw [← h0]; exact mul_nonneg hβ0 hy
  · have : (P *ᵥ y) ⬝ᵥ (P *ᵥ y) * ((P *ᵥ y) ⬝ᵥ (P *ᵥ y)) ≤
        (P *ᵥ y) ⬝ᵥ (P *ᵥ y) * (β * (y ⬝ᵥ P *ᵥ y)) := by
      nlinarith [mul_le_mul_of_nonneg_right h2 hy]
    exact le_of_mul_le_mul_left this hpos

/-- If `A⁻¹`'s form is at most `β · ‖·‖²`, then `‖x‖² ≤ β · xᵀAx`: `λ_min(A) ≥ 1/β`. -/
theorem dot_le_form_of_inv (hs : Aᵀ = A) (hu : IsUnit A.det) (hp : ∀ x, 0 ≤ x ⬝ᵥ A *ᵥ x)
    {β : ℝ} (hβ0 : 0 ≤ β) (hβ : ∀ z, z ⬝ᵥ A⁻¹ *ᵥ z ≤ β * (z ⬝ᵥ z)) (x : n → ℝ) :
    x ⬝ᵥ x ≤ β * (x ⬝ᵥ A *ᵥ x) := by
  have := sq_norm_le_form (inv_transpose_eq hs) (form_inv_nonneg hu hp) hβ0 hβ (A *ᵥ x)
  rwa [mulVec_mulVec, nonsing_inv_mul _ hu, one_mulVec, dotProduct_comm (A *ᵥ x)] at this

omit [DecidableEq n] in
/-- A positive definite matrix has a positive least Rayleigh quotient. -/
theorem lamMin_pos [Nonempty n] (hA : A.PosDef) : 0 < lamMin A := by
  classical
  have hs := transpose_eq hA.isHermitian
  have hu : IsUnit A.det := (isUnit_iff_isUnit_det A).1 hA.isUnit
  have hp : ∀ x, 0 ≤ x ⬝ᵥ A *ᵥ x := form_nonneg hA.posSemidef
  set K := ∑ i, ∑ j, |A⁻¹ i j|
  have hK0 : 0 ≤ K := sum_nonneg fun i _ => sum_nonneg fun j _ => abs_nonneg _
  have hle := dot_le_form_of_inv hs hu hp hK0
    (fun z => (le_abs_self _).trans (abs_form_le A⁻¹ z))
  have hKpos : 0 < K := by
    rcases hK0.eq_or_lt with h | h
    · have h1 := hle (fun _ => 1)
      rw [← h, zero_mul] at h1
      have : (0 : ℝ) < (fun _ : n => (1 : ℝ)) ⬝ᵥ (fun _ => 1) :=
        dot_self_pos fun h => by simpa using congrFun h (Classical.arbitrary n)
      linarith
    · exact h
  refine lt_of_lt_of_le (one_div_pos.2 hKpos) (le_lamMin fun x => ?_)
  rw [one_div, inv_mul_le_iff₀ hKpos]
  exact hle x

end Inverse

/-! ### PA-5 — identification does not read the coefficients -/

section Identification

variable [Fintype ι]

/-- The penalized system `A(c) = Σₖ cₖ Sₖ`. -/
def comb (S : ι → Matrix n n ℝ) (c : ι → ℝ) : Matrix n n ℝ := ∑ k, c k • S k

variable {S : ι → Matrix n n ℝ} {c c' : ι → ℝ}

theorem form_comb (S : ι → Matrix n n ℝ) (c : ι → ℝ) (x : n → ℝ) :
    x ⬝ᵥ comb S c *ᵥ x = ∑ k, c k * (x ⬝ᵥ S k *ᵥ x) := by
  simp only [comb, sum_mulVec, smul_mulVec, dotProduct_sum, dotProduct_smul, smul_eq_mul]

theorem comb_mulVec (S : ι → Matrix n n ℝ) (c : ι → ℝ) (x : n → ℝ) :
    comb S c *ᵥ x = ∑ k, c k • (S k *ᵥ x) := by
  simp only [comb, sum_mulVec, smul_mulVec]

omit [Fintype n] in
theorem comb_transpose (hS : ∀ k, (S k).PosSemidef) (c : ι → ℝ) : (comb S c)ᵀ = comb S c := by
  simp only [comb, transpose_sum, transpose_smul, fun k => transpose_eq (hS k).isHermitian]

theorem comb_form_nonneg (hS : ∀ k, (S k).PosSemidef) (hc : ∀ k, 0 ≤ c k) (x : n → ℝ) :
    0 ≤ x ⬝ᵥ comb S c *ᵥ x := by
  rw [form_comb]
  exact sum_nonneg fun k _ => mul_nonneg (hc k) (form_nonneg (hS k) x)

/-- **PA-5(a), the form.** For `c > 0`, `xᵀA(c)x = 0` exactly when `Sₖx = 0` for every `k`. -/
theorem form_comb_eq_zero_iff (hS : ∀ k, (S k).PosSemidef) (hc : ∀ k, 0 < c k) (x : n → ℝ) :
    x ⬝ᵥ comb S c *ᵥ x = 0 ↔ ∀ k, S k *ᵥ x = 0 := by
  constructor
  · intro h k
    rw [form_comb] at h
    have h0 := (sum_eq_zero_iff_of_nonneg fun k _ =>
      mul_nonneg (hc k).le (form_nonneg (hS k) x)).1 h k (mem_univ k)
    have hq : x ⬝ᵥ S k *ᵥ x = 0 := (mul_eq_zero.1 h0).resolve_left (hc k).ne'
    exact mulVec_eq_zero_of_form (transpose_eq (hS k).isHermitian) (form_nonneg (hS k)) hq
  · intro h
    rw [form_comb]
    exact sum_eq_zero fun k _ => by rw [h k, dotProduct_zero, mul_zero]

/-- **PA-5(a), the kernel.** For `c > 0`, `ker A(c) = ⋂ₖ ker Sₖ`. -/
theorem comb_mulVec_eq_zero_iff (hS : ∀ k, (S k).PosSemidef) (hc : ∀ k, 0 < c k) (x : n → ℝ) :
    comb S c *ᵥ x = 0 ↔ ∀ k, S k *ᵥ x = 0 := by
  constructor
  · intro h
    exact (form_comb_eq_zero_iff hS hc x).1 (by rw [h, dotProduct_zero])
  · intro h
    rw [comb_mulVec]
    exact sum_eq_zero fun k _ => by rw [h k, smul_zero]

/-- **PA-5(b).** For `c > 0`, `A(c)` is positive definite exactly when the `Sₖ` share no null
vector but `0`. -/
theorem comb_posDef_iff (hS : ∀ k, (S k).PosSemidef) (hc : ∀ k, 0 < c k) :
    (comb S c).PosDef ↔ ∀ x, (∀ k, S k *ᵥ x = 0) → x = 0 := by
  rw [posDef_iff_dotProduct_mulVec]
  have hherm : (comb S c).IsHermitian := by
    rw [IsHermitian, conjTranspose_eq_transpose_of_trivial, comb_transpose hS]
  constructor
  · rintro ⟨-, h⟩ x hx
    by_contra hne
    have := h hne
    rw [star_trivial, (form_comb_eq_zero_iff hS hc x).2 hx] at this
    exact lt_irrefl _ this
  · intro h
    refine ⟨hherm, fun x hx => ?_⟩
    rw [star_trivial]
    refine lt_of_le_of_ne (comb_form_nonneg hS (fun k => (hc k).le) x) fun h0 => hx ?_
    exact h x ((form_comb_eq_zero_iff hS hc x).1 h0.symm)

omit [Fintype n] in
/-- **PA-5(b).** Positive definiteness at one `c > 0` is positive definiteness at every `c > 0`. -/
theorem comb_posDef_iff_of_pos [Finite n] (hS : ∀ k, (S k).PosSemidef) (hc : ∀ k, 0 < c k)
    (hc' : ∀ k, 0 < c' k) : (comb S c).PosDef ↔ (comb S c').PosDef := by
  have := Fintype.ofFinite n
  rw [comb_posDef_iff hS hc, comb_posDef_iff hS hc']

/-- **PA-5(c), BAPS's instance.** With every `λⱼ > 0` and every `Pⱼ` positive semidefinite,
`ZᵀZ + Σⱼ λⱼPⱼ` is singular exactly when some `x ≠ 0` has `Zx = 0` and `Pⱼx = 0` for every `j`. -/
theorem penalized_det_eq_zero_iff {m J : Type*} [Fintype m] [Fintype J] [DecidableEq n]
    (Z : Matrix m n ℝ) (P : J → Matrix n n ℝ) (hP : ∀ j, (P j).PosSemidef) (lam : J → ℝ)
    (hlam : ∀ j, 0 < lam j) :
    (Zᵀ * Z + ∑ j, lam j • P j).det = 0 ↔ ∃ x ≠ 0, Z *ᵥ x = 0 ∧ ∀ j, P j *ᵥ x = 0 := by
  let S : Option J → Matrix n n ℝ := fun k => k.elim (Zᵀ * Z) P
  let c : Option J → ℝ := fun k => k.elim 1 lam
  have hS : ∀ k, (S k).PosSemidef := by
    rintro (_ | j)
    · simpa [S, conjTranspose_eq_transpose_of_trivial] using posSemidef_conjTranspose_mul_self Z
    · exact hP j
  have hc : ∀ k, 0 < c k := by
    rintro (_ | j)
    · exact one_pos
    · exact hlam j
  have hA : Zᵀ * Z + ∑ j, lam j • P j = comb S c := by
    simp [comb, Fintype.sum_option, S, c]
  have hZ : ∀ x, (Zᵀ * Z) *ᵥ x = 0 ↔ Z *ᵥ x = 0 := fun x => by
    constructor
    · intro h
      have : x ⬝ᵥ (Zᵀ * Z) *ᵥ x = (Z *ᵥ x) ⬝ᵥ (Z *ᵥ x) := by
        rw [← mulVec_mulVec, dotProduct_mulVec, ← mulVec_transpose, transpose_transpose,
          dotProduct_comm]
      rw [h, dotProduct_zero] at this
      exact dotProduct_self_eq_zero.1 this.symm
    · intro h; rw [← mulVec_mulVec, h, mulVec_zero]
  rw [hA, ← exists_mulVec_eq_zero_iff]
  refine exists_congr fun x => and_congr_right fun _ => ?_
  rw [comb_mulVec_eq_zero_iff hS hc]
  constructor
  · intro h; exact ⟨(hZ x).1 (h none), fun j => h (some j)⟩
  · rintro ⟨h0, h1⟩ (_ | j)
    · exact (hZ x).2 h0
    · exact h1 j

end Identification

/-! ### PA-6 — the bound -/

section Bound

variable [Fintype ι] {S : ι → Matrix n n ℝ} {c c' : ι → ℝ} {m M : ℝ}

/-- **PA-6(a).** With `m·c′ₖ ≤ cₖ ≤ M·c′ₖ` for every `k` (so `m ≤ min ρ` and `max ρ ≤ M`, where
`ρₖ = cₖ / c′ₖ`): `m·xᵀA(c′)x ≤ xᵀA(c)x ≤ M·xᵀA(c′)x`. -/
theorem form_comb_bound (hS : ∀ k, (S k).PosSemidef) (hm : ∀ k, m * c' k ≤ c k)
    (hM : ∀ k, c k ≤ M * c' k) (x : n → ℝ) :
    m * (x ⬝ᵥ comb S c' *ᵥ x) ≤ x ⬝ᵥ comb S c *ᵥ x ∧
      x ⬝ᵥ comb S c *ᵥ x ≤ M * (x ⬝ᵥ comb S c' *ᵥ x) := by
  rw [form_comb, form_comb, mul_sum, mul_sum]
  refine ⟨sum_le_sum fun k _ => ?_, sum_le_sum fun k _ => ?_⟩
  · rw [← mul_assoc]; exact mul_le_mul_of_nonneg_right (hm k) (form_nonneg (hS k) x)
  · rw [← mul_assoc]; exact mul_le_mul_of_nonneg_right (hM k) (form_nonneg (hS k) x)

/-- **PA-6(b), the smallest eigenvalue.** `λ_min(A(c)) ≥ m·λ_min(A(c′))`. -/
theorem lamMin_comb_ge [Nonempty n] (hS : ∀ k, (S k).PosSemidef) (hm0 : 0 ≤ m)
    (hm : ∀ k, m * c' k ≤ c k) (hM : ∀ k, c k ≤ M * c' k) :
    m * lamMin (comb S c') ≤ lamMin (comb S c) :=
  le_lamMin fun x => by
    have h1 := lamMin_mul_le (comb S c') x
    have h2 := (form_comb_bound hS hm hM x).1
    calc m * lamMin (comb S c') * (x ⬝ᵥ x) = m * (lamMin (comb S c') * (x ⬝ᵥ x)) := by ring
      _ ≤ m * (x ⬝ᵥ comb S c' *ᵥ x) := mul_le_mul_of_nonneg_left h1 hm0
      _ ≤ _ := h2

/-- **PA-6(b), the largest eigenvalue.** `λ_max(A(c)) ≤ M·λ_max(A(c′))`. -/
theorem lamMax_comb_le [Nonempty n] (hS : ∀ k, (S k).PosSemidef) (hM0 : 0 ≤ M)
    (hm : ∀ k, m * c' k ≤ c k) (hM : ∀ k, c k ≤ M * c' k) :
    lamMax (comb S c) ≤ M * lamMax (comb S c') :=
  lamMax_le fun x => by
    have h1 := form_le_lamMax_mul (comb S c') x
    have h2 := (form_comb_bound hS hm hM x).2
    calc _ ≤ M * (x ⬝ᵥ comb S c' *ᵥ x) := h2
      _ ≤ M * (lamMax (comb S c') * (x ⬝ᵥ x)) := mul_le_mul_of_nonneg_left h1 hM0
      _ = _ := by ring

/-- **PA-6(c).** Where `A(c′)` is positive definite, `κ₂(A(c)) ≤ (M/m)·κ₂(A(c′))`. -/
theorem cond2_comb_le [Nonempty n] (hS : ∀ k, (S k).PosSemidef) (hc : ∀ k, 0 ≤ c k)
    (hm0 : 0 < m) (hM0 : 0 ≤ M) (hm : ∀ k, m * c' k ≤ c k) (hM : ∀ k, c k ≤ M * c' k)
    (hpd : (comb S c').PosDef) :
    cond2 (comb S c) ≤ M / m * cond2 (comb S c') := by
  classical
  have hl' := lamMin_pos hpd
  have hl := lamMin_comb_ge hS hm0.le hm hM
  have hu := lamMax_comb_le hS hM0 hm hM
  have hx0 : 0 ≤ lamMax (comb S c) := lamMax_nonneg (comb_form_nonneg hS hc)
  unfold cond2
  calc lamMax (comb S c) / lamMin (comb S c)
      ≤ M * lamMax (comb S c') / (m * lamMin (comb S c')) :=
        div_le_div₀ (hx0.trans hu) hu (mul_pos hm0 hl') hl
    _ = M / m * (lamMax (comb S c') / lamMin (comb S c')) := by
        field_simp

omit [Fintype ι] in
/-- **PA-6(e).** A coefficient held at `1` in both enters with ratio `1`, so `m ≤ 1 ≤ M`. -/
theorem bracket_one (hm : ∀ k, m * c' k ≤ c k) (hM : ∀ k, c k ≤ M * c' k) {k₀ : ι}
    (h : c k₀ = 1) (h' : c' k₀ = 1) : m ≤ 1 ∧ 1 ≤ M := by
  have := hm k₀; have := hM k₀; rw [h, h'] at *; constructor <;> linarith

end Bound

/-! ### PA-6(d) — the constant is attained -/

section Sharp

/-- The smallest instance, `S₁ = diag(1, 0)` and `S₂ = diag(0, 1)`. -/
def diagPair : Fin 2 → Matrix (Fin 2) (Fin 2) ℝ :=
  ![diagonal ![1, 0], diagonal ![0, 1]]

theorem comb_diagPair (c : Fin 2 → ℝ) : comb diagPair c = diagonal ![c 0, c 1] := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [comb, diagPair, Fin.sum_univ_two]

theorem form_diagPair (c : Fin 2 → ℝ) (x : Fin 2 → ℝ) :
    x ⬝ᵥ comb diagPair c *ᵥ x = c 0 * x 0 ^ 2 + c 1 * x 1 ^ 2 := by
  rw [comb_diagPair]
  simp [dotProduct, mulVec_diagonal, Fin.sum_univ_two]
  ring

theorem diagPair_posSemidef : ∀ k, (diagPair k).PosSemidef := by
  intro k
  fin_cases k <;> simp [diagPair, posSemidef_diagonal_iff, Fin.forall_fin_two]

/-- The smallest instance meets the hypotheses: each `Sₖ` is singular, yet `A(c)` is positive
definite at every `c > 0`. -/
theorem comb_diagPair_posDef {c : Fin 2 → ℝ} (hc : ∀ k, 0 < c k) : (comb diagPair c).PosDef := by
  rw [comb_posDef_iff diagPair_posSemidef hc]
  intro x hx
  have h0 := congrFun (hx 0) 0
  have h1 := congrFun (hx 1) 1
  simp [diagPair, mulVec_diagonal] at h0 h1
  funext i
  fin_cases i <;> simp [h0, h1]

theorem lamMin_diagPair (c : Fin 2 → ℝ) : lamMin (comb diagPair c) = min (c 0) (c 1) := by
  apply le_antisymm
  · refine le_min ?_ ?_
    · have := lamMin_mul_le (comb diagPair c) ![1, 0]
      rw [form_diagPair] at this
      simpa [dotProduct, Fin.sum_univ_two] using this
    · have := lamMin_mul_le (comb diagPair c) ![0, 1]
      rw [form_diagPair] at this
      simpa [dotProduct, Fin.sum_univ_two] using this
  · refine le_lamMin fun x => ?_
    rw [form_diagPair]
    simp only [dotProduct, Fin.sum_univ_two]
    nlinarith [min_le_left (c 0) (c 1), min_le_right (c 0) (c 1), sq_nonneg (x 0),
      sq_nonneg (x 1)]

theorem lamMax_diagPair (c : Fin 2 → ℝ) : lamMax (comb diagPair c) = max (c 0) (c 1) := by
  apply le_antisymm
  · refine lamMax_le fun x => ?_
    rw [form_diagPair]
    simp only [dotProduct, Fin.sum_univ_two]
    nlinarith [le_max_left (c 0) (c 1), le_max_right (c 0) (c 1), sq_nonneg (x 0),
      sq_nonneg (x 1)]
  · refine max_le ?_ ?_
    · have := form_le_lamMax_mul (comb diagPair c) ![1, 0]
      rw [form_diagPair] at this
      simpa [dotProduct, Fin.sum_univ_two] using this
    · have := form_le_lamMax_mul (comb diagPair c) ![0, 1]
      rw [form_diagPair] at this
      simpa [dotProduct, Fin.sum_univ_two] using this

/-- **PA-6(d).** On the smallest instance with `c′ = (1, 1)`, `κ₂(A(c′)) = 1` and
`κ₂(A(c)) = max ρ / min ρ`, where `ρ = c`: the constant `M/m` is attained. -/
theorem cond2_diagPair (c : Fin 2 → ℝ) :
    cond2 (comb diagPair 1) = 1 ∧
      cond2 (comb diagPair c) = max (c 0) (c 1) / min (c 0) (c 1) * cond2 (comb diagPair 1) := by
  have h1 : cond2 (comb diagPair 1) = 1 := by
    simp [cond2, lamMin_diagPair, lamMax_diagPair]
  refine ⟨h1, ?_⟩
  rw [h1, mul_one, cond2, lamMin_diagPair, lamMax_diagPair]

end Sharp

/-! ### PA-7 — from the bound to what the gate reads -/

section Norms

/-- `‖A‖₁`, the largest column absolute sum. -/
noncomputable def norm1 (A : Matrix n n ℝ) : ℝ := ⨆ j, ∑ i, |A i j|

/-- `κ₁(A) = ‖A‖₁ ‖A⁻¹‖₁`. -/
noncomputable def cond1 [DecidableEq n] (A : Matrix n n ℝ) : ℝ := norm1 A * norm1 A⁻¹

theorem colSum_le_norm1 (A : Matrix n n ℝ) (j : n) : ∑ i, |A i j| ≤ norm1 A :=
  le_ciSup (f := fun j => ∑ i, |A i j|) (Set.finite_range _).bddAbove j

theorem norm1_le [Nonempty n] {A : Matrix n n ℝ} {β : ℝ} (h : ∀ j, ∑ i, |A i j| ≤ β) :
    norm1 A ≤ β :=
  ciSup_le h

theorem norm1_nonneg (A : Matrix n n ℝ) : 0 ≤ norm1 A :=
  Real.iSup_nonneg fun _ => sum_nonneg fun _ _ => abs_nonneg _

omit [Fintype n] in
theorem entry_symm {A : Matrix n n ℝ} (hs : Aᵀ = A) (i j : n) : A i j = A j i := by
  rw [← transpose_apply A j i, hs]

theorem single_form [DecidableEq n] (A : Matrix n n ℝ) (i j : n) :
    Pi.single i 1 ⬝ᵥ A *ᵥ Pi.single j 1 = A i j := by
  rw [single_dotProduct, one_mul]
  simp [mulVec, dotProduct_single]

theorem dot_sq_le (u v : n → ℝ) : (u ⬝ᵥ v) ^ 2 ≤ (u ⬝ᵥ u) * (v ⬝ᵥ v) := by
  classical
  have := form_sq_le (A := 1) transpose_one (fun x => by simpa using dot_self_nonneg x) v u
  simpa using this

/-- **PA-7(c3), Schur's test.** For a symmetric `A`, `‖Ax‖ ≤ ‖A‖₁ ‖x‖` (in squares), so
`‖A‖₂ ≤ ‖A‖₁`. -/
theorem mulVec_sq_le_norm1 {A : Matrix n n ℝ} (hs : Aᵀ = A) (x : n → ℝ) :
    (A *ᵥ x) ⬝ᵥ (A *ᵥ x) ≤ norm1 A ^ 2 * (x ⬝ᵥ x) := by
  set N := norm1 A
  have hN : 0 ≤ N := norm1_nonneg A
  have hrow : ∀ i, ∑ j, |A i j| ≤ N := fun i => by
    rw [show ∑ j, |A i j| = ∑ j, |A j i| from sum_congr rfl fun j _ => by rw [entry_symm hs]]
    exact colSum_le_norm1 A i
  have hi : ∀ i, (A *ᵥ x) i ^ 2 ≤ N * ∑ j, |A i j| * x j ^ 2 := fun i => by
    have hcs := sum_sq_le_sum_mul_sum_of_sq_le_mul univ (r := fun j => A i j * x j)
      (f := fun j => |A i j|) (g := fun j => |A i j| * x j ^ 2)
      (fun j _ => abs_nonneg _) (fun j _ => mul_nonneg (abs_nonneg _) (sq_nonneg _))
      (fun j _ => le_of_eq (by rw [mul_pow, ← sq_abs (A i j)]; ring))
    calc (A *ᵥ x) i ^ 2 = (∑ j, A i j * x j) ^ 2 := rfl
      _ ≤ (∑ j, |A i j|) * ∑ j, |A i j| * x j ^ 2 := hcs
      _ ≤ N * ∑ j, |A i j| * x j ^ 2 :=
          mul_le_mul_of_nonneg_right (hrow i)
            (sum_nonneg fun j _ => mul_nonneg (abs_nonneg _) (sq_nonneg _))
  calc (A *ᵥ x) ⬝ᵥ (A *ᵥ x) = ∑ i, (A *ᵥ x) i ^ 2 := sum_congr rfl fun i _ => (sq _).symm
    _ ≤ ∑ i, N * ∑ j, |A i j| * x j ^ 2 := sum_le_sum fun i _ => hi i
    _ = N * ∑ j, x j ^ 2 * ∑ i, |A i j| := by
        rw [← mul_sum, sum_comm]
        congr 1
        exact sum_congr rfl fun j _ => by rw [mul_sum]; exact sum_congr rfl fun i _ => by ring
    _ ≤ N * ∑ j, x j ^ 2 * N :=
        mul_le_mul_of_nonneg_left
          (sum_le_sum fun j _ => mul_le_mul_of_nonneg_left (colSum_le_norm1 A j) (sq_nonneg _)) hN
    _ = N ^ 2 * (x ⬝ᵥ x) := by
        rw [← sum_mul, show ∑ j, x j ^ 2 = x ⬝ᵥ x from sum_congr rfl fun j _ => sq _]
        ring

/-- For a symmetric `P`, `zᵀPz ≤ ‖P‖₁ ‖z‖²`. -/
theorem form_le_norm1 {P : Matrix n n ℝ} (hs : Pᵀ = P) (z : n → ℝ) :
    z ⬝ᵥ P *ᵥ z ≤ norm1 P * (z ⬝ᵥ z) := by
  have h1 := dot_sq_le z (P *ᵥ z)
  have h2 := mulVec_sq_le_norm1 hs z
  have hz := dot_self_nonneg z
  have hN := norm1_nonneg P
  have h3 : (z ⬝ᵥ P *ᵥ z) ^ 2 ≤ (norm1 P * (z ⬝ᵥ z)) ^ 2 := by
    calc _ ≤ (z ⬝ᵥ z) * ((P *ᵥ z) ⬝ᵥ (P *ᵥ z)) := h1
      _ ≤ (z ⬝ᵥ z) * (norm1 P ^ 2 * (z ⬝ᵥ z)) := mul_le_mul_of_nonneg_left h2 hz
      _ = _ := by ring
  exact abs_le_of_sq_le_sq' h3 (mul_nonneg hN hz) |>.2

/-- **PA-7(c2).** If `‖Ax‖ ≤ β‖x‖` for every `x` (in squares), then `‖A‖₁ ≤ √n · β`; so
`‖A‖₁ ≤ √n ‖A‖₂`. -/
theorem norm1_le_sqrt_card [Nonempty n] {A : Matrix n n ℝ} {β : ℝ} (hβ : 0 ≤ β)
    (h : ∀ x, (A *ᵥ x) ⬝ᵥ (A *ᵥ x) ≤ β ^ 2 * (x ⬝ᵥ x)) :
    norm1 A ≤ √(Fintype.card n) * β := by
  classical
  refine norm1_le fun j => ?_
  have he : ∀ i, (A *ᵥ Pi.single j 1) i = A i j := fun i => by simp [mulVec, dotProduct_single]
  have h1 := h (Pi.single j 1)
  rw [show (A *ᵥ Pi.single j 1) ⬝ᵥ (A *ᵥ Pi.single j 1) = ∑ i, A i j ^ 2 from
      sum_congr rfl fun i _ => by rw [he i, sq],
    show (Pi.single j (1 : ℝ)) ⬝ᵥ Pi.single j 1 = 1 by simp, mul_one] at h1
  have hcs : (∑ i, |A i j|) ^ 2 ≤ (Fintype.card n : ℝ) * ∑ i, A i j ^ 2 := by
    have := sum_mul_sq_le_sq_mul_sq univ (fun _ => (1 : ℝ)) (fun i => |A i j|)
    simpa [sq_abs] using this
  have hsq : (∑ i, |A i j|) ^ 2 ≤ (Fintype.card n : ℝ) * β ^ 2 :=
    hcs.trans (mul_le_mul_of_nonneg_left h1 (Nat.cast_nonneg _))
  calc ∑ i, |A i j| = √((∑ i, |A i j|) ^ 2) :=
        (Real.sqrt_sq (sum_nonneg fun i _ => abs_nonneg _)).symm
    _ ≤ √((Fintype.card n : ℝ) * β ^ 2) := Real.sqrt_le_sqrt hsq
    _ = √(Fintype.card n) * β := by rw [Real.sqrt_mul (Nat.cast_nonneg _), Real.sqrt_sq hβ]

/-- A positive semidefinite `E` with a unit diagonal has every entry in `[-1, 1]`. -/
theorem abs_entry_le_one {E : Matrix n n ℝ} (hs : Eᵀ = E)
    (hp : ∀ x, 0 ≤ x ⬝ᵥ E *ᵥ x) (hd : ∀ i, E i i = 1) (i j : n) : |E i j| ≤ 1 := by
  classical
  have := form_sq_le hs hp (Pi.single j 1) (Pi.single i 1)
  rw [single_form, single_form, single_form, hd, hd, mul_one] at this
  exact (sq_le_one_iff_abs_le_one _).1 this

/-- **PA-7(c1).** A positive semidefinite `E` with a unit diagonal and at most `q` nonzeros in any
column has `‖E‖₁ ≤ q`. -/
theorem norm1_le_card [Nonempty n] {E : Matrix n n ℝ} (hs : Eᵀ = E)
    (hp : ∀ x, 0 ≤ x ⬝ᵥ E *ᵥ x) (hd : ∀ i, E i i = 1) {q : ℕ}
    (hq : ∀ j, (univ.filter fun i => E i j ≠ 0).card ≤ q) : norm1 E ≤ q := by
  refine norm1_le fun j => ?_
  calc ∑ i, |E i j| = ∑ i ∈ univ.filter (fun i => E i j ≠ 0), |E i j| :=
        (sum_filter_of_ne fun i _ h => abs_ne_zero.1 h).symm
    _ ≤ ∑ i ∈ univ.filter (fun i => E i j ≠ 0), (1 : ℝ) :=
        sum_le_sum fun i _ => abs_entry_le_one hs hp hd i j
    _ = (univ.filter fun i => E i j ≠ 0).card := by simp
    _ ≤ q := by exact_mod_cast hq j

/-- `‖A⁻¹y‖ ≤ ‖y‖ / λ_min(A)`, in squares, for a symmetric invertible `A` with `λ_min(A) > 0`. -/
theorem inv_mulVec_sq_le [DecidableEq n] {A : Matrix n n ℝ} (hu : IsUnit A.det)
    (hl : 0 < lamMin A) (y : n → ℝ) :
    (A⁻¹ *ᵥ y) ⬝ᵥ (A⁻¹ *ᵥ y) ≤ (1 / lamMin A) ^ 2 * (y ⬝ᵥ y) := by
  set u := A⁻¹ *ᵥ y
  have hAu : A *ᵥ u = y := by rw [mulVec_mulVec, mul_nonsing_inv _ hu, one_mulVec]
  have h1 : lamMin A * (u ⬝ᵥ u) ≤ u ⬝ᵥ y := by
    have := lamMin_mul_le A u; rwa [hAu] at this
  have h0 : 0 ≤ lamMin A * (u ⬝ᵥ u) := mul_nonneg hl.le (dot_self_nonneg u)
  have h3 : (lamMin A * (u ⬝ᵥ u)) ^ 2 ≤ (u ⬝ᵥ u) * (y ⬝ᵥ y) :=
    (pow_le_pow_left₀ h0 h1 2).trans (dot_sq_le u y)
  rw [div_pow, one_pow, one_div_mul_eq_div, le_div_iff₀ (pow_pos hl 2)]
  rcases (dot_self_nonneg u).eq_or_lt with h | hpos
  · rw [← h, zero_mul]; exact dot_self_nonneg y
  · refine le_of_mul_le_mul_left ?_ hpos
    calc u ⬝ᵥ u * (u ⬝ᵥ u * lamMin A ^ 2) = (lamMin A * (u ⬝ᵥ u)) ^ 2 := by ring
      _ ≤ _ := h3

/-- `zᵀA⁻¹z ≤ ‖z‖² / λ_min(A)` for an invertible `A` with `λ_min(A) > 0`. -/
theorem form_inv_le [DecidableEq n] {A : Matrix n n ℝ} (hu : IsUnit A.det) (hl : 0 < lamMin A)
    (z : n → ℝ) : z ⬝ᵥ A⁻¹ *ᵥ z ≤ 1 / lamMin A * (z ⬝ᵥ z) := by
  have h1 := dot_sq_le z (A⁻¹ *ᵥ z)
  have h2 := inv_mulVec_sq_le hu hl z
  have hz := dot_self_nonneg z
  have h3 : (z ⬝ᵥ A⁻¹ *ᵥ z) ^ 2 ≤ (1 / lamMin A * (z ⬝ᵥ z)) ^ 2 := by
    calc _ ≤ (z ⬝ᵥ z) * ((A⁻¹ *ᵥ z) ⬝ᵥ (A⁻¹ *ᵥ z)) := h1
      _ ≤ (z ⬝ᵥ z) * ((1 / lamMin A) ^ 2 * (z ⬝ᵥ z)) := mul_le_mul_of_nonneg_left h2 hz
      _ = _ := by ring
  exact (abs_le_of_sq_le_sq' h3 (mul_nonneg (one_div_pos.2 hl).le hz)).2

/-- **PA-7(c2′), the column step, sharp.** A symmetric positive semidefinite `P` whose form is at
most `μ‖·‖²` has `‖P‖₁ ≤ μ (1 + √n) / 2`. A column `u = Peⱼ` has `‖u‖² ≤ μ uⱼ`, and the rest of
its absolute sum is at most `√(n−1) · √(μuⱼ − uⱼ²)`. A rank-one `μvvᵀ` attains the constant. -/
theorem norm1_le_half_one_add_sqrt [Nonempty n] {P : Matrix n n ℝ} (hs : Pᵀ = P)
    (hp : ∀ x, 0 ≤ x ⬝ᵥ P *ᵥ x) {μ : ℝ} (hμ0 : 0 ≤ μ) (hμ : ∀ z, z ⬝ᵥ P *ᵥ z ≤ μ * (z ⬝ᵥ z)) :
    norm1 P ≤ μ * (1 + √(Fintype.card n)) / 2 := by
  classical
  refine norm1_le fun j => ?_
  set u : n → ℝ := P *ᵥ Pi.single j 1 with hudef
  have hu : ∀ i, u i = P i j := fun i => by simp [hudef]
  have ha : u j = Pi.single j 1 ⬝ᵥ P *ᵥ Pi.single j 1 := by rw [single_form, hu]
  have hsq : u ⬝ᵥ u ≤ μ * u j := by
    have := sq_norm_le_form hs hp hμ0 hμ (Pi.single j 1)
    rwa [← ha] at this
  have ha0 : 0 ≤ u j := by rw [ha]; exact hp _
  have hN1 : (1 : ℝ) ≤ Fintype.card n := by exact_mod_cast Fintype.card_pos
  have haa : u j ^ 2 ≤ u ⬝ᵥ u := by
    have := Finset.single_le_sum (f := fun i => u i * u i) (fun i _ => mul_self_nonneg (u i))
      (Finset.mem_univ j)
    simpa [sq, dotProduct] using this
  have hsplit : ∑ i, |P i j| = u j + ∑ i ∈ univ.erase j, |u i| := by
    rw [← add_sum_erase _ _ (mem_univ j)]
    simp only [← hu, abs_of_nonneg ha0]
  have hrest : ∑ i ∈ univ.erase j, u i ^ 2 = u ⬝ᵥ u - u j ^ 2 := by
    rw [eq_sub_iff_add_eq, show u ⬝ᵥ u = ∑ i, u i * u i from rfl, ← add_sum_erase _ _ (mem_univ j)]
    simp only [sq]
    ring
  have hcard : ((univ.erase j).card : ℝ) = Fintype.card n - 1 := by
    rw [card_erase_of_mem (mem_univ j), card_univ, Nat.cast_sub Fintype.card_pos, Nat.cast_one]
  set T := ∑ i ∈ univ.erase j, |u i|
  have hT0 : 0 ≤ T := sum_nonneg fun i _ => abs_nonneg _
  have hT2 : T ^ 2 ≤ (Fintype.card n - 1) * (μ * u j - u j ^ 2) := by
    have := sum_mul_sq_le_sq_mul_sq (univ.erase j) (fun _ => (1 : ℝ)) (fun i => |u i|)
    simp only [one_mul, one_pow, sum_const, nsmul_eq_mul, mul_one, sq_abs] at this
    rw [hcard, hrest] at this
    exact this.trans (mul_le_mul_of_nonneg_left (by linarith) (by linarith))
  have hpos : 0 ≤ μ * u j - u j ^ 2 := by nlinarith
  set s := √(Fintype.card n - 1 : ℝ)
  set r := √(μ * u j - u j ^ 2)
  have hs2 : s ^ 2 = Fintype.card n - 1 := Real.sq_sqrt (by linarith)
  have hr2 : r ^ 2 = μ * u j - u j ^ 2 := Real.sq_sqrt hpos
  have hs0 : 0 ≤ s := Real.sqrt_nonneg _
  have hr0 : 0 ≤ r := Real.sqrt_nonneg _
  have hTsr : T ≤ s * r := by
    have h : T ^ 2 ≤ (s * r) ^ 2 := by rw [mul_pow, hs2, hr2]; exact hT2
    exact (abs_le_of_sq_le_sq' h (mul_nonneg hs0 hr0)).2
  set w := √(Fintype.card n : ℝ)
  have hw2 : w ^ 2 = Fintype.card n := Real.sq_sqrt (by linarith)
  have hw0 : 0 ≤ w := Real.sqrt_nonneg _
  have he : (w * μ / 2) ^ 2 - (u j - μ / 2 + s * r) ^ 2 = (s * (u j - μ / 2) - r) ^ 2 := by
    linear_combination (μ ^ 2 / 4) * hw2 - (μ ^ 2 / 4) * hs2 - (1 + s ^ 2) * hr2
  have hkey : (u j - μ / 2 + s * r) ^ 2 ≤ (w * μ / 2) ^ 2 := by
    nlinarith [sq_nonneg (s * (u j - μ / 2) - r)]
  have hkey' := (abs_le_of_sq_le_sq' hkey (by positivity)).2
  rw [hsplit]
  linarith

end Norms

section Scaling

/-- `E = D^{-1/2} A D^{-1/2}` with `D = diag A`: the unit-diagonal scaling. -/
noncomputable def unitScale (A : Matrix n n ℝ) : Matrix n n ℝ :=
  of fun i j => (√(A i i))⁻¹ * A i j * (√(A j j))⁻¹

theorem form_unitScale (A : Matrix n n ℝ) (x : n → ℝ) :
    x ⬝ᵥ unitScale A *ᵥ x =
      (fun i => (√(A i i))⁻¹ * x i) ⬝ᵥ A *ᵥ (fun i => (√(A i i))⁻¹ * x i) := by
  simp only [unitScale, dotProduct, mulVec, of_apply, mul_sum]
  exact sum_congr rfl fun i _ => sum_congr rfl fun j _ => by ring

omit [Fintype n] in
theorem unitScale_transpose {A : Matrix n n ℝ} (hs : Aᵀ = A) : (unitScale A)ᵀ = unitScale A := by
  ext i j
  simp only [transpose_apply, unitScale, of_apply, entry_symm hs j i]
  ring

omit [Fintype n] in
theorem unitScale_diag {A : Matrix n n ℝ} {i : n} (h : 0 < A i i) : unitScale A i i = 1 := by
  have h1 := Real.sqrt_pos.2 h
  have h2 := Real.mul_self_sqrt h.le
  simp only [unitScale, of_apply]
  field_simp
  linarith

omit [Fintype n] in
theorem unitScale_ne_zero {A : Matrix n n ℝ} {i j : n} (h : unitScale A i j ≠ 0) : A i j ≠ 0 :=
  fun h0 => h (by simp [unitScale, h0])

/-- A positive definite `A` has a positive definite unit-diagonal scaling. -/
theorem unitScale_form_pos {A : Matrix n n ℝ} (hA : A.PosDef) {x : n → ℝ} (hx : x ≠ 0) :
    0 < x ⬝ᵥ unitScale A *ᵥ x := by
  rw [form_unitScale]
  have hw : (fun i => (√(A i i))⁻¹ * x i) ≠ 0 := by
    intro h
    apply hx
    funext i
    have := congrFun h i
    simp only [Pi.zero_apply, mul_eq_zero, inv_eq_zero] at this
    exact this.resolve_left (Real.sqrt_pos.2 hA.diag_pos).ne'
  simpa using hA.dotProduct_mulVec_pos hw

variable [Fintype ι] {S : ι → Matrix n n ℝ} {c c' : ι → ℝ} {m M : ℝ}

omit [Fintype n] in
/-- **PA-7(a), the diagonal.** `m·A(c′)ᵢᵢ ≤ A(c)ᵢᵢ ≤ M·A(c′)ᵢᵢ`. -/
theorem diag_comb_bound [Finite n] (hS : ∀ k, (S k).PosSemidef) (hm : ∀ k, m * c' k ≤ c k)
    (hM : ∀ k, c k ≤ M * c' k) (i : n) :
    m * comb S c' i i ≤ comb S c i i ∧ comb S c i i ≤ M * comb S c' i i := by
  classical
  have := Fintype.ofFinite n
  have := form_comb_bound hS hm hM (Pi.single i 1)
  rwa [single_form, single_form] at this

/-- **PA-7(b).** With `A(c′)` positive definite, `λ_min(E(c)) ≥ (m/M)·λ_min(E(c′))`. -/
theorem lamMin_unitScale_ge [Nonempty n] (hS : ∀ k, (S k).PosSemidef) (hm0 : 0 < m)
    (hM0 : 0 < M) (hm : ∀ k, m * c' k ≤ c k) (hM : ∀ k, c k ≤ M * c' k)
    (hpd : (comb S c').PosDef) :
    m / M * lamMin (unitScale (comb S c')) ≤ lamMin (unitScale (comb S c)) := by
  set A := comb S c
  set A' := comb S c'
  have hd' : ∀ i, 0 < A' i i := fun i => hpd.diag_pos
  have hdiag := diag_comb_bound hS hm hM
  have hd : ∀ i, 0 < A i i := fun i => lt_of_lt_of_le (mul_pos hm0 (hd' i)) (hdiag i).1
  have hl0 : 0 ≤ lamMin (unitScale A') :=
    lamMin_nonneg fun x => by rw [form_unitScale]; exact form_nonneg hpd.posSemidef _
  refine le_lamMin fun x => ?_
  set w : n → ℝ := fun i => (√(A i i))⁻¹ * x i
  set v : n → ℝ := fun i => √(A' i i) * w i
  have hvw : (fun i => (√(A' i i))⁻¹ * v i) = w := funext fun i => by
    have := (Real.sqrt_pos.2 (hd' i)).ne'
    simp only [v]
    field_simp
  have e1 : x ⬝ᵥ unitScale A *ᵥ x = w ⬝ᵥ A *ᵥ w := form_unitScale A x
  have e2 : v ⬝ᵥ unitScale A' *ᵥ v = w ⬝ᵥ A' *ᵥ w := by rw [form_unitScale, hvw]
  have h3 : m * (w ⬝ᵥ A' *ᵥ w) ≤ w ⬝ᵥ A *ᵥ w := (form_comb_bound hS hm hM w).1
  have h4 := lamMin_mul_le (unitScale A') v
  have h5 : x ⬝ᵥ x ≤ M * (v ⬝ᵥ v) := by
    simp only [dotProduct, mul_sum]
    refine sum_le_sum fun i _ => ?_
    have hv : v i * v i =
        (√(A' i i) * √(A' i i)) * ((√(A i i))⁻¹ * (√(A i i))⁻¹) * (x i * x i) := by
      simp only [v, w]; ring
    rw [Real.mul_self_sqrt (hd' i).le, ← mul_inv, Real.mul_self_sqrt (hd i).le] at hv
    have hr : 1 ≤ M * A' i i * (A i i)⁻¹ := by
      rw [le_mul_inv_iff₀ (hd i), one_mul]; exact (hdiag i).2
    rw [hv]
    nlinarith [mul_self_nonneg (x i)]
  calc m / M * lamMin (unitScale A') * (x ⬝ᵥ x)
      = m * lamMin (unitScale A') * ((x ⬝ᵥ x) / M) := by ring
    _ ≤ m * lamMin (unitScale A') * (v ⬝ᵥ v) :=
        mul_le_mul_of_nonneg_left ((div_le_iff₀' hM0).2 h5) (mul_nonneg hm0.le hl0)
    _ = m * (lamMin (unitScale A') * (v ⬝ᵥ v)) := by ring
    _ ≤ m * (v ⬝ᵥ unitScale A' *ᵥ v) := mul_le_mul_of_nonneg_left h4 hm0.le
    _ = m * (w ⬝ᵥ A' *ᵥ w) := by rw [e2]
    _ ≤ w ⬝ᵥ A *ᵥ w := h3
    _ = _ := e1.symm

/-- `λ_min(E) ≥ 1/‖E⁻¹‖₁ > 0` for the unit-diagonal scaling of a positive definite symmetric `A`. -/
theorem inv_norm1_le_lamMin [Nonempty n] [DecidableEq n] {A : Matrix n n ℝ} (hs : Aᵀ = A)
    (hA : A.PosDef) :
    0 < norm1 (unitScale A)⁻¹ ∧ 1 / norm1 (unitScale A)⁻¹ ≤ lamMin (unitScale A) := by
  set E := unitScale A
  have hsE : Eᵀ = E := unitScale_transpose hs
  have hpdE : E.PosDef := posDef_iff_dotProduct_mulVec.2
    ⟨by rw [IsHermitian, conjTranspose_eq_transpose_of_trivial, hsE],
     fun x hx => by simpa using unitScale_form_pos hA hx⟩
  have huE : IsUnit E.det := (isUnit_iff_isUnit_det E).1 hpdE.isUnit
  have hpE : ∀ x, 0 ≤ x ⬝ᵥ E *ᵥ x := form_nonneg hpdE.posSemidef
  set N := norm1 E⁻¹
  have hle : ∀ x, x ⬝ᵥ x ≤ N * (x ⬝ᵥ E *ᵥ x) := dot_le_form_of_inv hsE huE hpE (norm1_nonneg _)
    (form_le_norm1 (inv_transpose_eq hsE))
  have hNpos : 0 < N := by
    rcases (show (0 : ℝ) ≤ N from norm1_nonneg _).eq_or_lt with h | h
    · have h1 := hle (fun _ => 1)
      rw [← h, zero_mul] at h1
      have : (0 : ℝ) < (fun _ : n => (1 : ℝ)) ⬝ᵥ (fun _ => 1) :=
        dot_self_pos fun h => by simpa using congrFun h (Classical.arbitrary n)
      linarith
    · exact h
  refine ⟨hNpos, le_lamMin fun x => ?_⟩
  rw [one_div, inv_mul_le_iff₀ hNpos]
  exact hle x

/-- **PA-7(d′), composed, with the sharp column step.** With `A(c′)` positive definite and at most
`q` nonzeros in any column of `A(c)`: `κ₁(E(c)) ≤ q · ((1 + √n)/2) · (M/m) · ‖E(c′)⁻¹‖₁`. -/
theorem cond1_unitScale_le_sharp [Nonempty n] [DecidableEq n] (hS : ∀ k, (S k).PosSemidef)
    (hm0 : 0 < m) (hM0 : 0 < M) (hm : ∀ k, m * c' k ≤ c k) (hM : ∀ k, c k ≤ M * c' k)
    (hpd : (comb S c').PosDef) {q : ℕ}
    (hq : ∀ j, (univ.filter fun i => comb S c i j ≠ 0).card ≤ q) :
    cond1 (unitScale (comb S c)) ≤
      q * ((1 + √(Fintype.card n)) / 2) * (M / m) * norm1 (unitScale (comb S c'))⁻¹ := by
  set A := comb S c
  set A' := comb S c'
  set E := unitScale A
  set N' := norm1 (unitScale A')⁻¹
  have hd' : ∀ i, 0 < A' i i := fun i => hpd.diag_pos
  have hdiag := diag_comb_bound hS hm hM
  have hd : ∀ i, 0 < A i i := fun i => lt_of_lt_of_le (mul_pos hm0 (hd' i)) (hdiag i).1
  have hsA : Aᵀ = A := comb_transpose hS c
  have hsE : Eᵀ = E := unitScale_transpose hsA
  have hApsd : ∀ x, 0 ≤ x ⬝ᵥ A *ᵥ x := fun x =>
    le_trans (mul_nonneg hm0.le (form_nonneg hpd.posSemidef x)) (form_comb_bound hS hm hM x).1
  have hpE : ∀ x, 0 ≤ x ⬝ᵥ E *ᵥ x := fun x => by rw [form_unitScale]; exact hApsd _
  -- `λ_min(E(c′)) ≥ 1/‖E(c′)⁻¹‖₁` and `λ_min(E(c)) ≥ (m/M)·λ_min(E(c′))`
  obtain ⟨hN'pos, hl'⟩ := inv_norm1_le_lamMin (comb_transpose hS c') hpd
  have hl := lamMin_unitScale_ge hS hm0 hM0 hm hM hpd
  have hlow : 0 < m / M * (1 / N') := mul_pos (div_pos hm0 hM0) (one_div_pos.2 hN'pos)
  have hlE : m / M * (1 / N') ≤ lamMin E :=
    (mul_le_mul_of_nonneg_left hl' (div_pos hm0 hM0).le).trans hl
  have hlEpos : 0 < lamMin E := hlow.trans_le hlE
  have huE : IsUnit E.det :=
    (isUnit_iff_isUnit_det E).1 (posDef_of_lamMin_pos hsE hlEpos).isUnit
  -- the two norms
  have hnE : norm1 E ≤ q := norm1_le_card hsE hpE (fun i => unitScale_diag (hd i))
    (fun j => (card_le_card fun i hi => by
      simp only [mem_filter] at hi ⊢; exact ⟨hi.1, unitScale_ne_zero hi.2⟩).trans (hq j))
  have hninv : norm1 E⁻¹ ≤ 1 / lamMin E * (1 + √(Fintype.card n)) / 2 :=
    norm1_le_half_one_add_sqrt (inv_transpose_eq hsE) (form_inv_nonneg huE hpE)
      (one_div_pos.2 hlEpos).le (form_inv_le huE hlEpos)
  have hrecip : 1 / lamMin E ≤ M / m * N' := by
    calc 1 / lamMin E ≤ 1 / (m / M * (1 / N')) := one_div_le_one_div_of_le hlow hlE
      _ = M / m * N' := by field_simp
  calc cond1 E = norm1 E * norm1 E⁻¹ := rfl
    _ ≤ q * (1 / lamMin E * (1 + √(Fintype.card n)) / 2) :=
        mul_le_mul hnE hninv (norm1_nonneg _) (Nat.cast_nonneg _)
    _ ≤ q * (M / m * N' * (1 + √(Fintype.card n)) / 2) := by
        gcongr
    _ = _ := by ring

/-- **PA-7(d), composed.** With `A(c′)` positive definite and at most `q` nonzeros in any column of
`A(c)`: `κ₁(E(c)) ≤ q · √n · (M/m) · ‖E(c′)⁻¹‖₁`. It follows from `cond1_unitScale_le_sharp`,
since `(1 + √n)/2 ≤ √n`. -/
theorem cond1_unitScale_le [Nonempty n] [DecidableEq n] (hS : ∀ k, (S k).PosSemidef)
    (hm0 : 0 < m) (hM0 : 0 < M) (hm : ∀ k, m * c' k ≤ c k) (hM : ∀ k, c k ≤ M * c' k)
    (hpd : (comb S c').PosDef) {q : ℕ}
    (hq : ∀ j, (univ.filter fun i => comb S c i j ≠ 0).card ≤ q) :
    cond1 (unitScale (comb S c)) ≤
      q * √(Fintype.card n) * (M / m) * norm1 (unitScale (comb S c'))⁻¹ := by
  have h := cond1_unitScale_le_sharp hS hm0 hM0 hm hM hpd hq
  have h1 : (1 : ℝ) ≤ √(Fintype.card n) := by
    rw [show (1 : ℝ) = √1 from Real.sqrt_one.symm]
    exact Real.sqrt_le_sqrt (by exact_mod_cast Fintype.card_pos)
  have h0 : 0 ≤ M / m * norm1 (unitScale (comb S c'))⁻¹ :=
    mul_nonneg (div_pos hM0 hm0).le (norm1_nonneg _)
  refine h.trans ?_
  calc (q : ℝ) * ((1 + √(Fintype.card n)) / 2) * (M / m) * norm1 (unitScale (comb S c'))⁻¹
      = q * ((1 + √(Fintype.card n)) / 2) * (M / m * norm1 (unitScale (comb S c'))⁻¹) := by ring
    _ ≤ q * √(Fintype.card n) * (M / m * norm1 (unitScale (comb S c'))⁻¹) :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left (by linarith) (Nat.cast_nonneg _)) h0
    _ = _ := by ring

end Scaling

/-! ### The Rayleigh extremes are the extreme eigenvalues -/

section Spectrum

variable [DecidableEq n] {A : Matrix n n ℝ}

/-- The spectral theorem in coordinates: with `U` the eigenvector matrix and `y = Uᵀx`,
`xᵀAx = Σᵢ λᵢ yᵢ²` and `yᵀy = xᵀx`. -/
theorem form_eq_sum_eigenvalues (hA : A.IsHermitian) (x : n → ℝ) :
    x ⬝ᵥ A *ᵥ x =
        ∑ i, hA.eigenvalues i * ((hA.eigenvectorUnitary : Matrix n n ℝ)ᵀ *ᵥ x) i ^ 2 ∧
      ((hA.eigenvectorUnitary : Matrix n n ℝ)ᵀ *ᵥ x) ⬝ᵥ
          ((hA.eigenvectorUnitary : Matrix n n ℝ)ᵀ *ᵥ x) = x ⬝ᵥ x := by
  set U : Matrix n n ℝ := ↑hA.eigenvectorUnitary
  have hstar : star U = Uᵀ := by rw [star_eq_conjTranspose, conjTranspose_eq_transpose_of_trivial]
  have hUU : U * Uᵀ = 1 := by rw [← hstar]; exact Unitary.coe_mul_star_self _
  have hspec : A = U * diagonal hA.eigenvalues * Uᵀ := by
    have h := hA.spectral_theorem
    rw [Unitary.conjStarAlgAut_apply] at h
    rw [← hstar]
    simpa [Unitary.coe_star, RCLike.ofReal_real_eq_id] using h
  set y := Uᵀ *ᵥ x
  have hxU : x ᵥ* U = y := by rw [← mulVec_transpose]
  constructor
  · conv_lhs => rw [hspec, ← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec, hxU]
    simp only [dotProduct, mulVec_diagonal]
    exact sum_congr rfl fun i _ => by ring
  · calc y ⬝ᵥ y = x ᵥ* U ⬝ᵥ y := by rw [hxU]
      _ = x ⬝ᵥ U *ᵥ (Uᵀ *ᵥ x) := (dotProduct_mulVec x U y).symm
      _ = x ⬝ᵥ x := by rw [mulVec_mulVec, hUU, one_mulVec]

theorem eigenvector_facts (hA : A.IsHermitian) (j : n) :
    ⇑(hA.eigenvectorBasis j) ⬝ᵥ A *ᵥ ⇑(hA.eigenvectorBasis j) = hA.eigenvalues j ∧
      ⇑(hA.eigenvectorBasis j) ⬝ᵥ ⇑(hA.eigenvectorBasis j) = 1 := by
  have hstar : star (hA.eigenvectorUnitary : Matrix n n ℝ) =
      (hA.eigenvectorUnitary : Matrix n n ℝ)ᵀ := by
    rw [star_eq_conjTranspose, conjTranspose_eq_transpose_of_trivial]
  have hy := (form_eq_sum_eigenvalues hA ⇑(hA.eigenvectorBasis j)).2
  rw [← hstar, hA.star_eigenvectorUnitary_mulVec j] at hy
  have h1 : ⇑(hA.eigenvectorBasis j) ⬝ᵥ ⇑(hA.eigenvectorBasis j) = 1 := by
    rw [← hy]; simp
  refine ⟨?_, h1⟩
  rw [hA.mulVec_eigenvectorBasis j, dotProduct_smul, h1, smul_eq_mul, mul_one]

/-- **The bridge.** For a symmetric `A`, the least Rayleigh quotient is the least eigenvalue. -/
theorem lamMin_eq_iInf_eigenvalues [Nonempty n] (hA : A.IsHermitian) :
    lamMin A = ⨅ i, hA.eigenvalues i := by
  apply le_antisymm
  · refine le_ciInf fun j => ?_
    have h := lamMin_mul_le A ⇑(hA.eigenvectorBasis j)
    rwa [(eigenvector_facts hA j).1, (eigenvector_facts hA j).2, mul_one] at h
  · refine le_lamMin fun x => ?_
    obtain ⟨h1, h2⟩ := form_eq_sum_eigenvalues hA x
    rw [h1, ← h2, dotProduct, mul_sum]
    refine sum_le_sum fun i _ => ?_
    rw [← sq]
    exact mul_le_mul_of_nonneg_right (ciInf_le (Set.finite_range _).bddBelow i) (sq_nonneg _)

/-- **The bridge.** For a symmetric `A`, the greatest Rayleigh quotient is the greatest
eigenvalue. -/
theorem lamMax_eq_iSup_eigenvalues [Nonempty n] (hA : A.IsHermitian) :
    lamMax A = ⨆ i, hA.eigenvalues i := by
  apply le_antisymm
  · refine lamMax_le fun x => ?_
    obtain ⟨h1, h2⟩ := form_eq_sum_eigenvalues hA x
    rw [h1, ← h2, dotProduct, mul_sum]
    refine sum_le_sum fun i _ => ?_
    rw [← sq]
    exact mul_le_mul_of_nonneg_right (le_ciSup (Set.finite_range _).bddAbove i) (sq_nonneg _)
  · refine ciSup_le fun j => ?_
    have h := form_le_lamMax_mul A ⇑(hA.eigenvectorBasis j)
    rwa [(eigenvector_facts hA j).1, (eigenvector_facts hA j).2, mul_one] at h

end Spectrum

/-- The column step's witness at `n = 4`: `P = vvᵀ/12` with `v = (3, 1, 1, 1)`. -/
noncomputable def columnWitness : Matrix (Fin 4) (Fin 4) ℝ :=
  (1 / 12 : ℝ) • vecMulVec ![3, 1, 1, 1] ![3, 1, 1, 1]

theorem form_columnWitness (z : Fin 4 → ℝ) :
    z ⬝ᵥ columnWitness *ᵥ z = (3 * z 0 + z 1 + z 2 + z 3) ^ 2 / 12 := by
  simp [columnWitness, vecMulVec, dotProduct, mulVec, Fin.sum_univ_four]
  ring

/-- **PA-7(c2′) is sharp.** `vvᵀ/12` with `v = (3, 1, 1, 1)` is symmetric positive semidefinite with
form at most `‖·‖²`, and its first column's absolute sum is `3/2 = (1 + √4)/2`. -/
theorem norm1_columnWitness :
    columnWitnessᵀ = columnWitness ∧ (∀ z, 0 ≤ z ⬝ᵥ columnWitness *ᵥ z) ∧
      (∀ z, z ⬝ᵥ columnWitness *ᵥ z ≤ 1 * (z ⬝ᵥ z)) ∧
      norm1 columnWitness = 1 * (1 + √(Fintype.card (Fin 4))) / 2 := by
  have hs : columnWitnessᵀ = columnWitness := by
    ext i j; simp [columnWitness, vecMulVec, mul_comm]
  have hp : ∀ z, 0 ≤ z ⬝ᵥ columnWitness *ᵥ z := fun z => by
    rw [form_columnWitness]; positivity
  have hb : ∀ z, z ⬝ᵥ columnWitness *ᵥ z ≤ 1 * (z ⬝ᵥ z) := fun z => by
    rw [form_columnWitness]
    simp only [dotProduct, Fin.sum_univ_four, one_mul]
    nlinarith [sq_nonneg (z 0 - 3 * z 1), sq_nonneg (z 0 - 3 * z 2), sq_nonneg (z 0 - 3 * z 3),
      sq_nonneg (z 1 - z 2), sq_nonneg (z 1 - z 3), sq_nonneg (z 2 - z 3)]
  have h4 : √(Fintype.card (Fin 4) : ℝ) = 2 := by
    rw [Fintype.card_fin, show ((4 : ℕ) : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]
  refine ⟨hs, hp, hb, le_antisymm (norm1_le_half_one_add_sqrt hs hp zero_le_one hb) ?_⟩
  rw [h4]
  have := colSum_le_norm1 columnWitness 0
  have hc : ∑ i, |columnWitness i 0| = 3 / 2 := by
    simp [columnWitness, vecMulVec, Fin.sum_univ_four, abs_of_nonneg]
    norm_num
  linarith

/-! ### The inverse part is unbounded over every family -/

section Dual

/-- **A dual lower bound on `‖E⁻¹‖₁`.** For a symmetric invertible `E` and any `y` with every
`|(Ey)ₖ| ≤ β`: `|yⱼ| ≤ ‖E⁻¹‖₁ · β`. -/
theorem abs_le_norm1_inv_mul [DecidableEq n] {E : Matrix n n ℝ} (hs : Eᵀ = E)
    (hu : IsUnit E.det) {y : n → ℝ} {β : ℝ} (hβ : ∀ k, |(E *ᵥ y) k| ≤ β) (j : n) :
    |y j| ≤ norm1 E⁻¹ * β := by
  have hy : y = E⁻¹ *ᵥ (E *ᵥ y) := by rw [mulVec_mulVec, nonsing_inv_mul _ hu, one_mulVec]
  have hsi : (E⁻¹)ᵀ = E⁻¹ := inv_transpose_eq hs
  have hβ0 : 0 ≤ β := (abs_nonneg _).trans (hβ j)
  calc |y j| = |(E⁻¹ *ᵥ (E *ᵥ y)) j| := by rw [← hy]
    _ = |∑ k, E⁻¹ j k * (E *ᵥ y) k| := rfl
    _ ≤ ∑ k, |E⁻¹ j k| * β := (abs_sum_le_sum_abs _ _).trans
        (sum_le_sum fun k _ => by
          rw [abs_mul]; exact mul_le_mul_of_nonneg_left (hβ k) (abs_nonneg _))
    _ = (∑ k, |E⁻¹ k j|) * β := by
        rw [sum_mul]; exact sum_congr rfl fun k _ => by rw [entry_symm hsi j k]
    _ ≤ norm1 E⁻¹ * β := mul_le_mul_of_nonneg_right (colSum_le_norm1 _ j) hβ0

/-- `E (D^{1/2} w) = D^{-1/2} (A w)` for the unit-diagonal scaling of an `A` with positive
diagonal. -/
theorem unitScale_mulVec {A : Matrix n n ℝ} (hd : ∀ i, 0 < A i i) (w : n → ℝ) :
    unitScale A *ᵥ (fun i => √(A i i) * w i) = fun i => (√(A i i))⁻¹ * (A *ᵥ w) i := by
  funext i
  simp only [unitScale, mulVec, dotProduct, of_apply, mul_sum]
  refine sum_congr rfl fun j _ => ?_
  have := (Real.sqrt_pos.2 (hd j)).ne'
  field_simp

/-! ### The inverse part is unbounded

On `r² + 1` coordinates, `v = (3/5, 4/(5r), …, 4/(5r))` is a unit vector. The family
`S₁ = I − vvᵀ`, `S₂ = vvᵀ` has `A(c′) = I` at `c′ = (1, 1)` and `A(c) = 2I − vvᵀ` at `c = (2, 1)`,
so `M/m = 2` and `‖E(c′)⁻¹‖₁ = 1`, while `‖E(c)⁻¹‖₁ ≥ 6r/25`. -/

section Unbounded

variable (r : ℕ)

/-- The unit vector `(3/5, 4/(5r), …, 4/(5r))` on `r² + 1` coordinates. -/
noncomputable def spikeVec : Fin (r ^ 2 + 1) → ℝ := fun i => if i = 0 then 3 / 5 else 4 / (5 * r)

/-- The two-term family `S₁ = I − vvᵀ`, `S₂ = vvᵀ`. -/
noncomputable def spikeFamily : Fin 2 → Matrix (Fin (r ^ 2 + 1)) (Fin (r ^ 2 + 1)) ℝ :=
  ![1 - vecMulVec (spikeVec r) (spikeVec r), vecMulVec (spikeVec r) (spikeVec r)]

variable {r}

theorem spikeVec_dot (hr : 1 ≤ r) : spikeVec r ⬝ᵥ spikeVec r = 1 := by
  have hr0 : (r : ℝ) ≠ 0 := by positivity
  simp only [dotProduct, Fin.sum_univ_succ, spikeVec, Fin.succ_ne_zero, if_true, if_false,
    sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
  push_cast
  field_simp
  ring

theorem vecMulVec_mulVec' (v w : Fin (r ^ 2 + 1) → ℝ) :
    vecMulVec v v *ᵥ w = (v ⬝ᵥ w) • v := by
  funext i
  simp [vecMulVec, mulVec, dotProduct, mul_sum, mul_comm, mul_left_comm]

theorem spikeFamily_posSemidef (hr : 1 ≤ r) : ∀ k, (spikeFamily r k).PosSemidef := by
  have hv := spikeVec_dot hr
  rw [Fin.forall_fin_two]
  constructor
  · change (1 - vecMulVec (spikeVec r) (spikeVec r)).PosSemidef
    refine posSemidef_iff_dotProduct_mulVec.2 ⟨?_, fun x => ?_⟩
    · rw [IsHermitian, conjTranspose_eq_transpose_of_trivial, transpose_sub, transpose_one,
        transpose_vecMulVec]
    · rw [star_trivial, sub_mulVec, one_mulVec, vecMulVec_mulVec', dotProduct_sub,
        dotProduct_smul, smul_eq_mul]
      have := dot_sq_le (spikeVec r) x
      rw [hv, one_mul] at this
      rw [dotProduct_comm x (spikeVec r)]
      nlinarith [this]
  · change (vecMulVec (spikeVec r) (spikeVec r)).PosSemidef
    simpa using posSemidef_vecMulVec_self_star (spikeVec r)

theorem comb_spike_one : comb (spikeFamily r) ![1, 1] = 1 := by
  simp only [comb, Fin.sum_univ_two, spikeFamily, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_fin_one, one_smul]
  abel

theorem comb_spike_two :
    comb (spikeFamily r) ![2, 1] =
      (2 : ℝ) • (1 : Matrix _ _ ℝ) - vecMulVec (spikeVec r) (spikeVec r) := by
  simp only [comb, Fin.sum_univ_two, spikeFamily, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_fin_one, one_smul, smul_sub]
  rw [two_smul, two_smul]
  abel

theorem unitScale_one : unitScale (1 : Matrix (Fin (r ^ 2 + 1)) (Fin (r ^ 2 + 1)) ℝ) = 1 := by
  ext i j
  by_cases h : i = j
  · subst h; simp [unitScale]
  · simp [unitScale, one_apply_ne h]

theorem norm1_one : norm1 (1 : Matrix (Fin (r ^ 2 + 1)) (Fin (r ^ 2 + 1)) ℝ) = 1 := by
  have h : ∀ j : Fin (r ^ 2 + 1),
      ∑ i, |(1 : Matrix (Fin (r ^ 2 + 1)) (Fin (r ^ 2 + 1)) ℝ) i j| = 1 :=
    fun j => by
      rw [Finset.sum_eq_single j (fun b _ hb => by simp [one_apply_ne hb]) (by simp)]
      simp
  simp only [norm1, h, ciSup_const]

/-- **The spike family's inverse part.** `‖E(c)⁻¹‖₁ ≥ 6r/25` at `c = (2, 1)`. The certificate is
`y = D^{1/2} w` with `Aw = (0, 5θ/(4r), …, 5θ/(4r))`, so `|(Ey)ₖ| ≤ 1` and `y₀ ≥ 6r/25`. -/
theorem norm1_inv_spike (hr : 1 ≤ r) :
    6 * (r : ℝ) / 25 ≤ norm1 (unitScale (comb (spikeFamily r) ![2, 1]))⁻¹ := by
  classical
  have hv := spikeVec_dot hr
  have hr0 : (0 : ℝ) < r := by exact_mod_cast hr
  have hr1 : (1 : ℝ) ≤ r := by exact_mod_cast hr
  set v := spikeVec r with hvdef
  set A := comb (spikeFamily r) ![2, 1] with hAdef
  have hA : A = (2 : ℝ) • (1 : Matrix _ _ ℝ) - vecMulVec v v := comb_spike_two
  have hAii : ∀ i, A i i = 2 - v i ^ 2 := fun i => by rw [hA]; simp [vecMulVec, sq]
  have hv2 : ∀ i, v i ^ 2 ≤ 1 := fun i => by
    by_cases h : i = 0
    · simp only [hvdef, spikeVec, h, if_true]; norm_num
    · simp only [hvdef, spikeVec, h, if_false]
      rw [div_pow, div_le_one (by positivity)]
      nlinarith
  have hd : ∀ i, 1 ≤ A i i := fun i => by rw [hAii]; linarith [hv2 i]
  have hd0 : ∀ i, 0 < A i i := fun i => by linarith [hd i]
  have hsA : Aᵀ = A := by
    rw [hA, transpose_sub, transpose_smul, transpose_one, transpose_vecMulVec]
  have hAmul : ∀ w, A *ᵥ w = (2 : ℝ) • w - (v ⬝ᵥ w) • v := fun w => by
    rw [hA, sub_mulVec, smul_mulVec, one_mulVec, vecMulVec_mulVec']
  have hApd : A.PosDef := by
    refine posDef_iff_dotProduct_mulVec.2 ⟨?_, fun x hx => ?_⟩
    · rw [IsHermitian, conjTranspose_eq_transpose_of_trivial, hsA]
    · rw [star_trivial, hAmul, dotProduct_sub, dotProduct_smul, dotProduct_smul, smul_eq_mul,
        smul_eq_mul]
      have h1 := dot_sq_le v x
      rw [hv, one_mul] at h1
      rw [dotProduct_comm x v]
      have h2 := dot_self_pos hx
      nlinarith
  set E := unitScale A
  have hsE : Eᵀ = E := unitScale_transpose hsA
  have hpdE : E.PosDef := posDef_iff_dotProduct_mulVec.2
    ⟨by rw [IsHermitian, conjTranspose_eq_transpose_of_trivial, hsE],
     fun x hx => by simpa using unitScale_form_pos hApd hx⟩
  have huE : IsUnit E.det := (isUnit_iff_isUnit_det E).1 hpdE.isUnit
  -- the certificate
  set Δ : ℝ := 2 - (4 / (5 * r)) ^ 2 with hΔdef
  have hΔ : 1 ≤ Δ := by
    rw [hΔdef, div_pow]
    have : (4 : ℝ) ^ 2 / (5 * r) ^ 2 ≤ 1 := by rw [div_le_one (by positivity)]; nlinarith
    linarith
  have hsΔ : 1 ≤ √Δ := by
    rw [show (1 : ℝ) = √1 from Real.sqrt_one.symm]; exact Real.sqrt_le_sqrt hΔ
  set θ : ℝ := 4 * r * √Δ / 5 with hθ
  set w : Fin (r ^ 2 + 1) → ℝ :=
    fun i => if i = 0 then 3 * θ / 10 else 41 * θ / (40 * r) with hw
  have hvw : v ⬝ᵥ w = θ := by
    simp only [dotProduct, Fin.sum_univ_succ, hvdef, spikeVec, hw, Fin.succ_ne_zero, if_true,
      if_false, sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
    push_cast
    field_simp
    ring
  have hAw : ∀ i, (A *ᵥ w) i = if i = 0 then 0 else 5 * θ / (4 * r) := fun i => by
    rw [hAmul, hvw]
    by_cases h : i = 0
    · simp only [h, if_true, Pi.sub_apply, Pi.smul_apply, smul_eq_mul, hw, hvdef, spikeVec]; ring
    · simp only [h, if_false, Pi.sub_apply, Pi.smul_apply, smul_eq_mul, hw, hvdef, spikeVec]
      field_simp
      ring
  have hEy : ∀ k, |(E *ᵥ fun i => √(A i i) * w i) k| ≤ 1 := fun k => by
    rw [unitScale_mulVec hd0]
    change |(√(A k k))⁻¹ * (A *ᵥ w) k| ≤ 1
    rw [hAw k]
    by_cases h : k = 0
    · simp [h]
    · rw [if_neg h, hAii k]
      have hvk : v k = 4 / (5 * r) := by simp [hvdef, spikeVec, h]
      rw [hvk, ← hΔdef]
      have hpos : 0 < √Δ := by linarith
      have : 5 * θ / (4 * r) = √Δ := by rw [hθ]; field_simp
      rw [this, inv_mul_cancel₀ hpos.ne', abs_one]
  have hy0 := abs_le_norm1_inv_mul hsE huE hEy 0
  have hy0v : √(A 0 0) * w 0 = √(A 0 0) * (3 * θ / 10) := by simp [hw]
  have hs0 : 1 ≤ √(A 0 0) := by
    rw [show (1 : ℝ) = √1 from Real.sqrt_one.symm]; exact Real.sqrt_le_sqrt (hd 0)
  have hθ0 : 0 ≤ 3 * θ / 10 := by rw [hθ]; positivity
  have hθlow : 6 * (r : ℝ) / 25 ≤ 3 * θ / 10 := by
    rw [hθ]; nlinarith
  calc 6 * (r : ℝ) / 25 ≤ 3 * θ / 10 := hθlow
    _ ≤ √(A 0 0) * (3 * θ / 10) := le_mul_of_one_le_left hθ0 hs0
    _ ≤ |√(A 0 0) * w 0| := by rw [hy0v]; exact le_abs_self _
    _ ≤ norm1 E⁻¹ * 1 := hy0
    _ = norm1 E⁻¹ := mul_one _

/-- **The conjecture is false.** No constant bounds `‖E(c)⁻¹‖₁ / ((M/m)·‖E(c′)⁻¹‖₁)` over every
family: for every `C`, the spike family on `r² + 1` coordinates, `r > 25C/3`, has `m = 1`, `M = 2`,
`A(c′) = I` positive definite and `‖E(c)⁻¹‖₁ > C·(M/m)·‖E(c′)⁻¹‖₁`. -/
theorem inverse_part_unbounded (C : ℝ) :
    ∃ (k : ℕ) (S : Fin 2 → Matrix (Fin (k + 1)) (Fin (k + 1)) ℝ) (c c' : Fin 2 → ℝ),
      (∀ i, (S i).PosSemidef) ∧ (∀ i, 1 * c' i ≤ c i) ∧ (∀ i, c i ≤ 2 * c' i) ∧
      (comb S c').PosDef ∧
      C * (2 / 1) * norm1 (unitScale (comb S c'))⁻¹ < norm1 (unitScale (comb S c))⁻¹ := by
  set r := ⌈25 * C / 3⌉₊ + 1 with hrdef
  have hr : 1 ≤ r := Nat.le_add_left 1 _
  refine ⟨r ^ 2, spikeFamily r, ![2, 1], ![1, 1], spikeFamily_posSemidef hr, ?_, ?_, ?_, ?_⟩
  · rw [Fin.forall_fin_two]; norm_num
  · rw [Fin.forall_fin_two]; norm_num
  · rw [comb_spike_one]; exact PosDef.one
  · rw [comb_spike_one, unitScale_one, inv_one, norm1_one, mul_one]
    have h1 := norm1_inv_spike hr
    have h2 : 25 * C / 3 < r := by
      rw [hrdef]; push_cast
      exact lt_of_le_of_lt (Nat.le_ceil _) (lt_add_one _)
    linarith

end Unbounded

end Dual

/-! ### The smallest instance

`S₁ = diag(1, 0)` and `S₂ = diag(0, 1)`, each singular, with `c′ = (1, 1)`, `c = (2, 3)`, `m = 2`,
`M = 3` and `q = 1`, meet every main theorem's hypotheses. BAPS's form is met by `Z = (1 0)` and
`P = diag(0, 1)`. -/

section Instance

theorem card_filter_diagonal_le_one (d : n → ℝ) [DecidableEq n] (j : n) :
    (univ.filter fun i => diagonal d i j ≠ 0).card ≤ 1 :=
  Finset.card_le_one.2 fun a ha b hb => by
    simp only [mem_filter, diagonal_apply] at ha hb
    split_ifs at ha hb with h1 h2 <;> simp_all

theorem bracket_low : ∀ k, (2 : ℝ) * (1 : Fin 2 → ℝ) k ≤ ![2, 3] k := by
  intro k; fin_cases k <;> norm_num

theorem bracket_high : ∀ k, (![2, 3] : Fin 2 → ℝ) k ≤ 3 * (1 : Fin 2 → ℝ) k := by
  intro k; fin_cases k <;> norm_num

example : (comb diagPair ![2, 3]).PosDef ↔ (comb diagPair 1).PosDef :=
  comb_posDef_iff_of_pos diagPair_posSemidef (fun k => by fin_cases k <;> norm_num)
    fun _ => one_pos

example : (comb diagPair 1).PosDef := comb_diagPair_posDef fun _ => one_pos

example (x : Fin 2 → ℝ) :
    2 * (x ⬝ᵥ comb diagPair 1 *ᵥ x) ≤ x ⬝ᵥ comb diagPair ![2, 3] *ᵥ x ∧
      x ⬝ᵥ comb diagPair ![2, 3] *ᵥ x ≤ 3 * (x ⬝ᵥ comb diagPair 1 *ᵥ x) :=
  form_comb_bound diagPair_posSemidef bracket_low bracket_high x

example : cond2 (comb diagPair ![2, 3]) ≤ 3 / 2 * cond2 (comb diagPair 1) :=
  cond2_comb_le diagPair_posSemidef (fun k => by fin_cases k <;> norm_num) (by norm_num)
    (by norm_num) bracket_low bracket_high (comb_diagPair_posDef fun _ => one_pos)

example :
    2 / 3 * lamMin (unitScale (comb diagPair 1)) ≤ lamMin (unitScale (comb diagPair ![2, 3])) :=
  lamMin_unitScale_ge diagPair_posSemidef (by norm_num) (by norm_num) bracket_low bracket_high
    (comb_diagPair_posDef fun _ => one_pos)

example : cond1 (unitScale (comb diagPair ![2, 3])) ≤
    (1 : ℕ) * √(Fintype.card (Fin 2)) * (3 / 2) * norm1 (unitScale (comb diagPair 1))⁻¹ :=
  cond1_unitScale_le diagPair_posSemidef (by norm_num) (by norm_num) bracket_low bracket_high
    (comb_diagPair_posDef fun _ => one_pos) fun j => by
      rw [comb_diagPair]; exact card_filter_diagonal_le_one _ j

/-- BAPS's form on the smallest instance: `Z = (1 0)`. -/
def rowZ : Matrix (Fin 1) (Fin 2) ℝ := !![1, 0]

/-- BAPS's form on the smallest instance: one penalty, `P = diag(0, 1)`. -/
noncomputable def penaltyP : Fin 1 → Matrix (Fin 2) (Fin 2) ℝ := ![diagonal ![0, 1]]

theorem penaltyP_posSemidef : ∀ j, (penaltyP j).PosSemidef := by
  intro j
  fin_cases j
  simp [penaltyP, posSemidef_diagonal_iff, Fin.forall_fin_two]

example := penalized_det_eq_zero_iff rowZ penaltyP penaltyP_posSemidef (1 : Fin 1 → ℝ)
  fun _ => one_pos

end Instance

end PenalizedAdmission
