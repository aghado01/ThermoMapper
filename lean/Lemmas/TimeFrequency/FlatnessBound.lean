import Mathlib.Analysis.Real.Sqrt
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.LinearAlgebra.Matrix.Hermitian
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.LinearAlgebra.Matrix.Symmetric
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# The flatness bound (S3)

How flat a quadratic `xᵀHx` can be made over a time–frequency cell of fixed area. For
`x ~ N(0, Σ)`, `Var(xᵀHx) = 2·tr((HΣ)²)` (S15, `TimeFrequency.variance_quadForm`, under the
pairing rule), and for any factor `Σ = S·Sᵀ` that trace is `‖SᵀHS‖²_F` (`trace_sq_eq_frobSq`). In
two dimensions:

* `two_mul_abs_det_le_frobSq`: a real symmetric 2×2 `B` has `‖B‖²_F ≥ 2·|det B|`, with equality
  exactly when `B` is a multiple of the identity or has trace zero
  (`frobSq_eq_two_mul_abs_det_iff`).
* `frobSq_congr_ge`: `‖SᵀHS‖²_F ≥ 2·(det S)²·|det H|` for any real 2×2 `S` and symmetric `H`;
  `trace_sq_ge`: `tr((HΣ)²) ≥ 2·det Σ·|det H|` for `Σ = S·Sᵀ`.
* `exists_posDef_frobSq_eq`: for `det H ≠ 0` and any `d > 0`, a positive definite `S` with
  `det S = d` attains the bound. So over cells of fixed area the least variance is fixed by
  `|det H|` alone, whatever the cell's shape.

`frobSq` is the squared Frobenius norm written as a sum, so the statements need no normed
instance on matrices.
-/

open Matrix

namespace TimeFrequency

/-- The squared Frobenius norm, `Σᵢⱼ Bᵢⱼ²`. -/
def frobSq {m n : Type*} [Fintype m] [Fintype n] (B : Matrix m n ℝ) : ℝ :=
  ∑ i, ∑ j, B i j ^ 2

theorem frobSq_fin_two (B : Matrix (Fin 2) (Fin 2) ℝ) :
    frobSq B = B 0 0 ^ 2 + B 0 1 ^ 2 + B 1 0 ^ 2 + B 1 1 ^ 2 := by
  simp [frobSq, Fin.sum_univ_two]
  ring

theorem frobSq_smul {m n : Type*} [Fintype m] [Fintype n] (c : ℝ) (B : Matrix m n ℝ) :
    frobSq (c • B) = c ^ 2 * frobSq B := by
  simp only [frobSq, Matrix.smul_apply, smul_eq_mul, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  ring

/-- A congruence `SᵀHS` of a symmetric `H` is symmetric. -/
theorem isSymm_transpose_mul_mul {n : Type*} [Fintype n] {H : Matrix n n ℝ} (hH : H.IsSymm)
    (S : Matrix n n ℝ) : (Sᵀ * H * S).IsSymm := by
  rw [Matrix.IsSymm, transpose_mul, transpose_mul, transpose_transpose, hH.eq, Matrix.mul_assoc]

/-- For symmetric `B`, `tr(B·B) = ‖B‖²_F`. -/
theorem trace_mul_self_eq_frobSq {n : Type*} [Fintype n] {B : Matrix n n ℝ} (hB : B.IsSymm) :
    (B * B).trace = frobSq B := by
  simp only [trace, diag, mul_apply, frobSq]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  rw [hB.apply j i, sq]

/-- **The variance link.** For symmetric `H` and any factor `Σ = S·Sᵀ`,
`tr((HΣ)²) = ‖SᵀHS‖²_F`; so `Var(xᵀHx) = 2‖SᵀHS‖²_F` for `x ~ N(0, Σ)`, whatever the factor. -/
theorem trace_sq_eq_frobSq {n : Type*} [Fintype n] {H : Matrix n n ℝ}
    (hH : H.IsSymm) (S : Matrix n n ℝ) :
    ((H * (S * Sᵀ)) * (H * (S * Sᵀ))).trace = frobSq (Sᵀ * H * S) := by
  rw [← trace_mul_self_eq_frobSq (isSymm_transpose_mul_mul hH S)]
  have : H * (S * Sᵀ) * (H * (S * Sᵀ)) = (H * S * Sᵀ * H * S) * Sᵀ := by
    simp only [Matrix.mul_assoc]
  rw [this, trace_mul_comm]
  simp only [Matrix.mul_assoc]

/-- **S3.** A real symmetric 2×2 matrix has `‖B‖²_F ≥ 2·|det B|`. -/
theorem two_mul_abs_det_le_frobSq {B : Matrix (Fin 2) (Fin 2) ℝ} (hB : B.IsSymm) :
    2 * |B.det| ≤ frobSq B := by
  rw [frobSq_fin_two, det_fin_two, hB.apply 0 1]
  rcases abs_cases (B 0 0 * B 1 1 - B 0 1 * B 0 1) with ⟨h, _⟩ | ⟨h, _⟩ <;> rw [h] <;>
    nlinarith [sq_nonneg (B 0 0 - B 1 1), sq_nonneg (B 0 0 + B 1 1), sq_nonneg (B 0 1)]

/-- **S3, the equality case.** `‖B‖²_F = 2·|det B|` exactly when `B` is a multiple of the
identity (then `det B ≥ 0`) or has trace zero (then `det B ≤ 0`). -/
theorem frobSq_eq_two_mul_abs_det_iff {B : Matrix (Fin 2) (Fin 2) ℝ} (hB : B.IsSymm) :
    frobSq B = 2 * |B.det| ↔ (∃ a : ℝ, B = a • (1 : Matrix (Fin 2) (Fin 2) ℝ)) ∨ B.trace = 0 := by
  have h10 := hB.apply 0 1
  rw [frobSq_fin_two, det_fin_two, trace_fin_two, h10]
  constructor
  · intro h
    rcases abs_cases (B 0 0 * B 1 1 - B 0 1 * B 0 1) with ⟨ha, _⟩ | ⟨ha, _⟩ <;> rw [ha] at h
    · left
      have h1 : B 0 1 = 0 := by nlinarith [sq_nonneg (B 0 0 - B 1 1), sq_nonneg (B 0 1)]
      have h2 : B 0 0 = B 1 1 := by nlinarith [sq_nonneg (B 0 0 - B 1 1), sq_nonneg (B 0 1)]
      refine ⟨B 0 0, ?_⟩
      ext i j
      fin_cases i <;> fin_cases j <;> simp [h1, h2, h10]
    · right
      nlinarith [sq_nonneg (B 0 0 + B 1 1)]
  · rintro (⟨a, rfl⟩ | h)
    · simp
      ring
    · have hc : B 1 1 = -B 0 0 := by linarith
      rw [hc, abs_of_nonpos (by nlinarith [sq_nonneg (B 0 0), sq_nonneg (B 0 1)])]
      ring

/-- **S3, the congruence form.** For any real 2×2 `S` and symmetric `H`,
`‖SᵀHS‖²_F ≥ 2·(det S)²·|det H|`. -/
theorem frobSq_congr_ge (S : Matrix (Fin 2) (Fin 2) ℝ) {H : Matrix (Fin 2) (Fin 2) ℝ}
    (hH : H.IsSymm) : 2 * S.det ^ 2 * |H.det| ≤ frobSq (Sᵀ * H * S) := by
  have h := two_mul_abs_det_le_frobSq (isSymm_transpose_mul_mul hH S)
  have hd : (Sᵀ * H * S).det = S.det ^ 2 * H.det := by
    rw [det_mul, det_mul, det_transpose]
    ring
  rw [hd, abs_mul, abs_of_nonneg (sq_nonneg _)] at h
  linarith

/-- **S3 in the covariance's terms.** For symmetric `H` and `Σ = S·Sᵀ`,
`tr((HΣ)²) ≥ 2·det Σ·|det H|`, so `Var(xᵀHx) ≥ 4·det Σ·|det H|` for `x ~ N(0, Σ)`. -/
theorem trace_sq_ge (S : Matrix (Fin 2) (Fin 2) ℝ) {H : Matrix (Fin 2) (Fin 2) ℝ}
    (hH : H.IsSymm) :
    2 * (S * Sᵀ).det * |H.det| ≤ ((H * (S * Sᵀ)) * (H * (S * Sᵀ))).trace := by
  rw [trace_sq_eq_frobSq hH, det_mul, det_transpose, ← sq]
  exact frobSq_congr_ge S hH

/-- Sylvester's criterion in two dimensions: a symmetric `S` with `S₀₀ > 0` and `det S > 0` is
positive definite. -/
theorem posDef_fin_two {S : Matrix (Fin 2) (Fin 2) ℝ} (hS : S.IsSymm) (h0 : 0 < S 0 0)
    (hd : 0 < S.det) : S.PosDef := by
  refine PosDef.of_dotProduct_mulVec_pos (isHermitian_iff_isSymm.2 hS) fun x hx => ?_
  rw [det_fin_two, hS.apply 0 1] at hd
  have hQ : star x ⬝ᵥ (S *ᵥ x) =
      x 0 * (S 0 0 * x 0 + S 0 1 * x 1) + x 1 * (S 0 1 * x 0 + S 1 1 * x 1) := by
    simp only [star_trivial, dotProduct, mulVec, Fin.sum_univ_two, hS.apply 0 1]
  rw [hQ]
  have key : S 0 0 * (x 0 * (S 0 0 * x 0 + S 0 1 * x 1) + x 1 * (S 0 1 * x 0 + S 1 1 * x 1)) =
      (S 0 0 * x 0 + S 0 1 * x 1) ^ 2 + (S 0 0 * S 1 1 - S 0 1 * S 0 1) * x 1 ^ 2 := by ring
  by_cases h1 : x 1 = 0
  · have h0' : x 0 ≠ 0 := fun h => hx (funext fun i => by fin_cases i <;> simp [h, h1])
    rw [h1]
    simp only [mul_zero, add_zero, zero_mul]
    nlinarith [sq_pos_of_ne_zero h0', mul_pos h0 (sq_pos_of_ne_zero h0')]
  · have : 0 < S 0 0 * (x 0 * (S 0 0 * x 0 + S 0 1 * x 1) + x 1 * (S 0 1 * x 0 + S 1 1 * x 1)) := by
      rw [key]
      nlinarith [sq_nonneg (S 0 0 * x 0 + S 0 1 * x 1), mul_pos hd (sq_pos_of_ne_zero h1)]
    exact pos_of_mul_pos_right this h0.le

/-- A witness of the equality case scales to any determinant `d > 0`. -/
theorem exists_posDef_of_witness {H S₀ : Matrix (Fin 2) (Fin 2) ℝ} (hH : H.IsSymm)
    (hS : S₀.IsSymm) (h0 : 0 < S₀ 0 0) (hdet : 0 < S₀.det)
    (heq : (∃ κ : ℝ, S₀ᵀ * H * S₀ = κ • (1 : Matrix (Fin 2) (Fin 2) ℝ)) ∨
      (S₀ᵀ * H * S₀).trace = 0)
    {d : ℝ} (hd : 0 < d) :
    ∃ S : Matrix (Fin 2) (Fin 2) ℝ,
      S.PosDef ∧ S.det = d ∧ frobSq (Sᵀ * H * S) = 2 * d ^ 2 * |H.det| := by
  have hq : 0 < d / S₀.det := div_pos hd hdet
  set t := Real.sqrt (d / S₀.det) with ht_def
  have ht : 0 < t := Real.sqrt_pos.2 hq
  have ht2 : t ^ 2 = d / S₀.det := Real.sq_sqrt hq.le
  have hdetS : (t • S₀).det = d := by
    rw [det_smul, Fintype.card_fin, ht2]
    field_simp
  have h00 : 0 < (t • S₀) 0 0 := by
    rw [Matrix.smul_apply, smul_eq_mul]
    exact mul_pos ht h0
  refine ⟨t • S₀, posDef_fin_two (hS.smul t) h00 (hdetS ▸ hd), hdetS, ?_⟩
  have hcong : (t • S₀)ᵀ * H * (t • S₀) = (t * t) • (S₀ᵀ * H * S₀) := by
    simp only [transpose_smul, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  have hB := isSymm_transpose_mul_mul hH (t • S₀)
  rw [(frobSq_eq_two_mul_abs_det_iff hB).2 ?_]
  · rw [det_mul, det_mul, det_transpose, hdetS, abs_mul, abs_mul, abs_of_pos hd]
    ring
  · rw [hcong]
    rcases heq with ⟨κ, hκ⟩ | htr
    · exact Or.inl ⟨t * t * κ, by rw [hκ, smul_smul]⟩
    · exact Or.inr (by rw [trace_smul, htr, smul_zero])

theorem isSymm_fin_two (x y z : ℝ) : (!![x, y; y, z]).IsSymm := by
  ext i j
  fin_cases i <;> fin_cases j <;> rfl

/-- **S3, attainment.** For symmetric `H` with `det H ≠ 0` and any `d > 0`, some positive
definite `S` with `det S = d` attains `‖SᵀHS‖²_F = 2·d²·|det H|`. A definite `H` takes `S`
with `S² ∝ H⁻¹`, making `SᵀHS` scalar; an indefinite one takes `S² ∝ ‖H‖²_F·I − tr H·H`,
making it traceless. -/
theorem exists_posDef_frobSq_eq {H : Matrix (Fin 2) (Fin 2) ℝ} (hH : H.IsSymm) (hdet : H.det ≠ 0)
    {d : ℝ} (hd : 0 < d) :
    ∃ S : Matrix (Fin 2) (Fin 2) ℝ,
      S.PosDef ∧ S.det = d ∧ frobSq (Sᵀ * H * S) = 2 * d ^ 2 * |H.det| := by
  obtain ⟨a, b, c, rfl⟩ : ∃ a b c : ℝ, H = !![a, b; b, c] :=
    ⟨H 0 0, H 0 1, H 1 1, by ext i j; fin_cases i <;> fin_cases j <;> simp [hH.apply 0 1]⟩
  rw [det_fin_two_of] at hdet
  rcases lt_or_gt_of_ne hdet with hneg | hpos
  · -- indefinite: S² ∝ M = ‖H‖²_F·I − tr H·H, with entries p, q, r
    have hp : 0 < 2 * b ^ 2 + c ^ 2 - a * c := by nlinarith [sq_nonneg c, sq_nonneg b]
    have hr : 0 < a ^ 2 + 2 * b ^ 2 - a * c := by nlinarith [sq_nonneg a, sq_nonneg b]
    have hgap : 0 < (a - c) ^ 2 + 4 * b ^ 2 := by
      by_contra h
      nlinarith [sq_nonneg (a + c), sq_nonneg (a - c), sq_nonneg b]
    have hD : 0 < (2 * b ^ 2 + c ^ 2 - a * c) * (a ^ 2 + 2 * b ^ 2 - a * c) -
        (-(a + c) * b) ^ 2 := by
      have e : (2 * b ^ 2 + c ^ 2 - a * c) * (a ^ 2 + 2 * b ^ 2 - a * c) - (-(a + c) * b) ^ 2 =
          (b * b - a * c) * ((a - c) ^ 2 + 4 * b ^ 2) := by ring
      rw [e]
      exact mul_pos (by linarith) hgap
    set δ := Real.sqrt ((2 * b ^ 2 + c ^ 2 - a * c) * (a ^ 2 + 2 * b ^ 2 - a * c) -
      (-(a + c) * b) ^ 2) with hδ_def
    have hδ0 : 0 < δ := Real.sqrt_pos.2 hD
    have hδ : δ ^ 2 = (2 * b ^ 2 + c ^ 2 - a * c) * (a ^ 2 + 2 * b ^ 2 - a * c) -
        (-(a + c) * b) ^ 2 := Real.sq_sqrt hD.le
    refine exists_posDef_of_witness (S₀ := !![2 * b ^ 2 + c ^ 2 - a * c + δ, -(a + c) * b;
      -(a + c) * b, a ^ 2 + 2 * b ^ 2 - a * c + δ]) hH (isSymm_fin_two _ _ _) ?_ ?_ (Or.inr ?_) hd
    · simp only [of_apply, cons_val', cons_val_zero, empty_val', cons_val_fin_one]
      linarith
    · rw [det_fin_two_of]
      nlinarith [mul_pos hδ0 (add_pos hp hr)]
    · simp [Matrix.mul_apply, Fin.sum_univ_two, trace_fin_two]
      linear_combination (a + c) * hδ
  · have ha0 : a ≠ 0 := by
      rintro rfl
      nlinarith [sq_nonneg b]
    rcases lt_or_gt_of_ne ha0 with ha | ha
    · -- negative definite: S² ∝ −adj H, and SᵀHS = Δ(a + c − 2δ)·I
      have hc : c < 0 := by nlinarith [sq_nonneg b]
      set δ := Real.sqrt (a * c - b * b) with hδ_def
      have hδ0 : 0 < δ := Real.sqrt_pos.2 hpos
      have hδ : δ ^ 2 = a * c - b * b := Real.sq_sqrt hpos.le
      refine exists_posDef_of_witness (S₀ := !![-c + δ, b; b, -a + δ]) hH (isSymm_fin_two _ _ _)
        ?_ ?_ (Or.inl ⟨(a * c - b * b) * (a + c - 2 * δ), ?_⟩) hd
      · simp only [of_apply, cons_val', cons_val_zero, empty_val', cons_val_fin_one]
        linarith
      · rw [det_fin_two_of]
        nlinarith [mul_pos hδ0 (show 0 < -(a + c) by linarith)]
      · ext i j
        fin_cases i <;> fin_cases j
        · simp [Matrix.mul_apply, Fin.sum_univ_two]
          linear_combination a * hδ
        · simp [Matrix.mul_apply, Fin.sum_univ_two]
          linear_combination b * hδ
        · simp [Matrix.mul_apply, Fin.sum_univ_two]
          linear_combination b * hδ
        · simp [Matrix.mul_apply, Fin.sum_univ_two]
          linear_combination c * hδ
    · -- positive definite: S² ∝ adj H, and SᵀHS = Δ(a + c + 2δ)·I
      have hc : 0 < c := by nlinarith [sq_nonneg b]
      set δ := Real.sqrt (a * c - b * b) with hδ_def
      have hδ0 : 0 < δ := Real.sqrt_pos.2 hpos
      have hδ : δ ^ 2 = a * c - b * b := Real.sq_sqrt hpos.le
      refine exists_posDef_of_witness (S₀ := !![c + δ, -b; -b, a + δ]) hH (isSymm_fin_two _ _ _)
        ?_ ?_ (Or.inl ⟨(a * c - b * b) * (a + c + 2 * δ), ?_⟩) hd
      · simp only [of_apply, cons_val', cons_val_zero, empty_val', cons_val_fin_one]
        linarith
      · rw [det_fin_two_of]
        nlinarith [mul_pos hδ0 (show 0 < a + c by linarith)]
      · ext i j
        fin_cases i <;> fin_cases j
        · simp [Matrix.mul_apply, Fin.sum_univ_two]
          linear_combination a * hδ
        · simp [Matrix.mul_apply, Fin.sum_univ_two]
          linear_combination b * hδ
        · simp [Matrix.mul_apply, Fin.sum_univ_two]
          linear_combination b * hδ
        · simp [Matrix.mul_apply, Fin.sum_univ_two]
          linear_combination c * hδ

/-! ### Witnesses -/

/-- The bound's hypotheses hold on a concrete pair, and its equality case is reached: `S = I`
and the traceless `H = diag(1, −1)`. -/
example : frobSq ((1 : Matrix (Fin 2) (Fin 2) ℝ)ᵀ * !![1, 0; 0, -1] * 1) =
    2 * (1 : Matrix (Fin 2) (Fin 2) ℝ).det ^ 2 * |(!![(1 : ℝ), 0; 0, -1]).det| := by
  simp [frobSq_fin_two, det_fin_two_of]
  norm_num

/-- The 2×2 bound applied to `B = diag(1, 2)`, where it is strict: `‖B‖²_F = 5 > 4 = 2·|det B|`. -/
example : 2 * |(!![(1 : ℝ), 0; 0, 2]).det| ≤ frobSq !![(1 : ℝ), 0; 0, 2] :=
  two_mul_abs_det_le_frobSq (isSymm_fin_two _ _ _)

example : 2 * |(!![(1 : ℝ), 0; 0, 2]).det| < frobSq !![(1 : ℝ), 0; 0, 2] := by
  simp [frobSq_fin_two, det_fin_two_of]
  norm_num

/-- The attainment's hypotheses are met: `H = diag(1, −1)` and `d = 2`. -/
example : ∃ S : Matrix (Fin 2) (Fin 2) ℝ, S.PosDef ∧ S.det = 2 ∧
    frobSq (Sᵀ * !![1, 0; 0, -1] * S) = 2 * 2 ^ 2 * |(!![(1 : ℝ), 0; 0, -1]).det| :=
  exists_posDef_frobSq_eq (isSymm_fin_two _ _ _) (by simp [det_fin_two_of]) (by norm_num)

end TimeFrequency
