import Mathlib.Data.Real.Basic
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.LinearAlgebra.Matrix.Vec

/-!
# The Kronecker Gram (S5)

On a complete grid a product design is the Kronecker product `Z = A ⊗ B` of its axes' designs,
and its Gram and projection factor by axis. That is the BAMRS brief's full-grid shortcut, which
reads the fitted energy as `‖P_A·Y·P_B‖²_F`. Scattered rows give only the row-wise
(face-splitting) product, whose Gram `Σᵣ (aᵣaᵣᵀ) ⊗ (bᵣbᵣᵀ)` does not factor. BASS's functional
design pairs columns instead (column-wise Khatri–Rao), and its Gram is a Hadamard product,
`BASS.khatriRao_gram`.

* `kronecker_gram`: `(A ⊗ B)ᵀ(A ⊗ B) = (AᵀA) ⊗ (BᵀB)`.
* `proj_kronecker`: `P_{A ⊗ B} = P_A ⊗ P_B` for `P_A = A(AᵀA)⁻¹Aᵀ` (`proj`). The identity needs no
  invertibility; `P_A` is the orthogonal projection onto `A`'s columns when `AᵀA` is invertible
  (`proj_transpose`, `proj_mul_self`).
* The vec convention. `y (i, j) = Y i j`, the `A`-index outer: Mathlib's `vec Yᵀ`.
  `rowMajor_energy`: `yᵀ(P_A ⊗ P_B)y = ‖P_A·Y·P_B‖²_F`. Under Mathlib's column-major `vec Y` the
  factors swap: `colMajor_energy`, `vec Yᵀ(P_B ⊗ P_A) vec Y = ‖P_A·Y·P_B‖²_F`.
-/

open Matrix
open scoped Kronecker

namespace BAMRS

variable {m n p q : Type*} [Fintype m] [Fintype n] [Fintype p] [Fintype q]

/-- The squared Frobenius norm as a trace: `Σᵢⱼ Xᵢⱼ² = tr(XᵀX)`. -/
theorem sum_sq_eq_trace (X : Matrix m n ℝ) : ∑ i, ∑ j, X i j ^ 2 = (Xᵀ * X).trace := by
  simp only [trace, diag, mul_apply, transpose_apply, sq]
  exact Finset.sum_comm

omit [Fintype p] [Fintype q] in
/-- **S5, the Gram.** `(A ⊗ B)ᵀ(A ⊗ B) = (AᵀA) ⊗ (BᵀB)`. -/
theorem kronecker_gram (A : Matrix m p ℝ) (B : Matrix n q ℝ) :
    (A ⊗ₖ B)ᵀ * (A ⊗ₖ B) = (Aᵀ * A) ⊗ₖ (Bᵀ * B) := by
  rw [← kroneckerMap_transpose, ← mul_kronecker_mul]

variable [DecidableEq p] [DecidableEq q]

/-- The projection onto the columns of `A`, `A(AᵀA)⁻¹Aᵀ`, with Mathlib's inverse (zero when `AᵀA`
is singular). -/
noncomputable def proj (A : Matrix m p ℝ) : Matrix m m ℝ :=
  A * (Aᵀ * A)⁻¹ * Aᵀ

/-- **S5, the projection.** `P_{A ⊗ B} = P_A ⊗ P_B`. -/
theorem proj_kronecker (A : Matrix m p ℝ) (B : Matrix n q ℝ) :
    proj (A ⊗ₖ B) = proj A ⊗ₖ proj B := by
  rw [proj, kronecker_gram, inv_kronecker, ← kroneckerMap_transpose, ← mul_kronecker_mul,
    ← mul_kronecker_mul]
  rfl

theorem proj_transpose (A : Matrix m p ℝ) : (proj A)ᵀ = proj A := by
  simp only [proj, transpose_mul, transpose_transpose, transpose_nonsing_inv, Matrix.mul_assoc]

theorem proj_mul_self (A : Matrix m p ℝ) (h : IsUnit (Aᵀ * A).det) :
    proj A * proj A = proj A := by
  simp only [proj, Matrix.mul_assoc]
  rw [← Matrix.mul_assoc Aᵀ A, ← Matrix.mul_assoc (Aᵀ * A)⁻¹ (Aᵀ * A), nonsing_inv_mul _ h,
    Matrix.one_mul]

/-- **S5, the energy, column-major.** With Mathlib's `vec Y` (columns stacked, `B`-index outer),
`vec Yᵀ(P_B ⊗ P_A) vec Y = ‖P_A·Y·P_B‖²_F`. -/
theorem colMajor_energy (A : Matrix m p ℝ) (B : Matrix n q ℝ) (hA : IsUnit (Aᵀ * A).det)
    (hB : IsUnit (Bᵀ * B).det) (Y : Matrix m n ℝ) :
    vec Y ⬝ᵥ ((proj B ⊗ₖ proj A) *ᵥ vec Y) = ∑ i, ∑ j, (proj A * Y * proj B) i j ^ 2 := by
  rw [kronecker_mulVec_vec, proj_transpose, vec_dotProduct_vec, sum_sq_eq_trace]
  have hX : (proj A * Y * proj B)ᵀ * (proj A * Y * proj B) =
      proj B * (Yᵀ * proj A * Y * proj B) := by
    simp only [transpose_mul, proj_transpose, Matrix.mul_assoc]
    rw [← Matrix.mul_assoc (proj A) (proj A), proj_mul_self A hA]
  rw [hX, trace_mul_comm (proj B), Matrix.mul_assoc (Yᵀ * proj A * Y) (proj B) (proj B),
    proj_mul_self B hB]
  simp only [Matrix.mul_assoc]

/-- **S5, the energy, as posed.** For `y (i, j) = Y i j` (the `A`-index outer, Mathlib's
`vec Yᵀ`), `yᵀ(P_A ⊗ P_B)y = ‖P_A·Y·P_B‖²_F`. -/
theorem rowMajor_energy (A : Matrix m p ℝ) (B : Matrix n q ℝ) (hA : IsUnit (Aᵀ * A).det)
    (hB : IsUnit (Bᵀ * B).det) (Y : Matrix m n ℝ) :
    vec Yᵀ ⬝ᵥ ((proj A ⊗ₖ proj B) *ᵥ vec Yᵀ) = ∑ i, ∑ j, (proj A * Y * proj B) i j ^ 2 := by
  rw [colMajor_energy B A hB hA Yᵀ]
  have : proj B * Yᵀ * proj A = (proj A * Y * proj B)ᵀ := by
    simp only [transpose_mul, proj_transpose, Matrix.mul_assoc]
  rw [this]
  simp only [transpose_apply]
  exact Finset.sum_comm

/-! ### Witness -/

/-- The energy's hypotheses hold for the one-column design `(1, 1)ᵀ` on both axes. -/
example (Y : Matrix (Fin 2) (Fin 2) ℝ) :
    vec Yᵀ ⬝ᵥ ((proj !![(1 : ℝ); 1] ⊗ₖ proj !![(1 : ℝ); 1]) *ᵥ vec Yᵀ) =
      ∑ i, ∑ j, (proj !![(1 : ℝ); 1] * Y * proj !![(1 : ℝ); 1]) i j ^ 2 :=
  have h : IsUnit ((!![(1 : ℝ); 1])ᵀ * !![(1 : ℝ); 1]).det :=
    isUnit_iff_ne_zero.2 (by simp [Matrix.mul_apply, Fin.sum_univ_two])
  rowMajor_energy _ _ h h Y

end BAMRS
