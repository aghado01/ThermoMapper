import Mathlib.Analysis.Real.Sqrt
import Mathlib.LinearAlgebra.Matrix.Symmetric
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.MeasureTheory.Integral.Bochner.SumMeasure
import Mathlib.MeasureTheory.Measure.Count
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# The variance of a Gaussian quadratic form (S15)

This closes S3's reading. For a real zero-mean vector `x` with second moments `S` whose fourth
moments obey the Gaussian pairing rule (Isserlis), `Var(xᵀHx) = 2·tr((HS)²)` for symmetric `H`.

* `Pairing S m4`: `E[xₚx_qx_rx_s] = S_pq·S_rs + S_pr·S_qs + S_ps·S_qr`, the premise. Mathlib
  does not supply Isserlis's theorem for Gaussian vectors, so the rule enters by name.
* `pairing_variance`: the algebra. The three pairings give `(tr HS)²` and twice `tr((HS)²)`.
* `variance_quadForm`: the expectation form. Any random vector with integrable moments to order
  four, second moments `S` and paired fourth moments has
  `E[(xᵀHx)²] − (E[xᵀHx])² = 2·tr((HS)²)`.
* The witness meets the premise with a six-point law in one dimension: `±√3` with weight `1/6`
  each, `0` with weight `2/3`. So `E x² = 1` and `E x⁴ = 3`.
-/

open Matrix Finset MeasureTheory

namespace TimeFrequency

variable {n : Type*} [Fintype n]

/-- The Gaussian pairing rule for fourth moments `m4` against second moments `S`. -/
def Pairing (S : Matrix n n ℝ) (m4 : n → n → n → n → ℝ) : Prop :=
  ∀ p q r s, m4 p q r s = S p q * S r s + S p r * S q s + S p s * S q r

theorem trace_mul_mul_eq_sum (H S : Matrix n n ℝ) :
    ((H * S) * (H * S)).trace =
      ∑ p, ∑ a, ∑ c, ∑ b, H p b * S b a * (H a c * S c p) := by
  rw [Matrix.trace]
  simp only [Matrix.diag, mul_apply, Finset.sum_mul, Finset.mul_sum]

/-- One cross pairing, `Σ H_pq·H_rs·S_pr·S_qs = tr((HS)²)`. -/
theorem sum_pairing_cross₁ {H S : Matrix n n ℝ} (hH : H.IsSymm) (hS : S.IsSymm) :
    ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * (S p r * S q s) = ((H * S) * (H * S)).trace := by
  rw [trace_mul_mul_eq_sum]
  refine Finset.sum_congr rfl fun p _ => ?_
  symm
  conv_lhs => arg 2; ext a; rw [Finset.sum_comm]
  conv_lhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun q _ => ?_
  conv_lhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun r _ => Finset.sum_congr rfl fun s _ => ?_
  rw [hH.apply r s, hS.apply p r]
  ring

/-- The other cross pairing, `Σ H_pq·H_rs·S_ps·S_qr = tr((HS)²)`. -/
theorem sum_pairing_cross₂ {H S : Matrix n n ℝ} (hS : S.IsSymm) :
    ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * (S p s * S q r) = ((H * S) * (H * S)).trace := by
  rw [trace_mul_mul_eq_sum]
  refine Finset.sum_congr rfl fun p _ => ?_
  symm
  conv_lhs => arg 2; ext a; rw [Finset.sum_comm]
  conv_lhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun q _ => Finset.sum_congr rfl fun r _ =>
    Finset.sum_congr rfl fun s _ => ?_
  rw [hS.apply p s]
  ring

/-- The direct pairing, `Σ H_pq·H_rs·S_pq·S_rs = (Σ H_pq·S_pq)²`. -/
theorem sum_pairing_direct (H S : Matrix n n ℝ) :
    ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * (S p q * S r s) = (∑ p, ∑ q, H p q * S p q) ^ 2 := by
  rw [sq, Finset.sum_mul]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun q _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun r _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun s _ => ?_
  ring

/-- **S15, the algebra.** Under the pairing rule and symmetric `H`, `S`,
`Σ H_pq·H_rs·m4_pqrs − (Σ H_pq·S_pq)² = 2·tr((HS)²)`. -/
theorem pairing_variance {H S : Matrix n n ℝ} (hH : H.IsSymm) (hS : S.IsSymm)
    {m4 : n → n → n → n → ℝ} (hp : Pairing S m4) :
    ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * m4 p q r s - (∑ p, ∑ q, H p q * S p q) ^ 2 =
      2 * ((H * S) * (H * S)).trace := by
  have hsplit : ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * m4 p q r s =
      ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * (S p q * S r s) +
        ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * (S p r * S q s) +
          ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * (S p s * S q r) := by
    simp only [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun q _ =>
      Finset.sum_congr rfl fun r _ => Finset.sum_congr rfl fun s _ => ?_
    rw [hp]
    ring
  rw [hsplit, sum_pairing_direct, sum_pairing_cross₁ hH hS, sum_pairing_cross₂ hS]
  ring

/-- `xᵀHx` as a double sum. -/
theorem dotProduct_mulVec_eq_sum (H : Matrix n n ℝ) (x : n → ℝ) :
    x ⬝ᵥ (H *ᵥ x) = ∑ p, ∑ q, H p q * (x p * x q) := by
  simp only [dotProduct, mulVec, Finset.mul_sum]
  refine Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun q _ => ?_
  ring

/-- **S15.** A random vector with integrable moments to order four, second moments `S` and paired
fourth moments has `E[(xᵀHx)²] − (E[xᵀHx])² = 2·tr((HS)²)` for symmetric `H`. -/
theorem variance_quadForm {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) (x : Ω → n → ℝ)
    {H : Matrix n n ℝ} (hH : H.IsSymm) (S : Matrix n n ℝ)
    (h2 : ∀ p q, Integrable (fun w => x w p * x w q) μ ∧ ∫ w, x w p * x w q ∂μ = S p q)
    (m4 : n → n → n → n → ℝ)
    (h4 : ∀ p q r s, Integrable (fun w => x w p * x w q * (x w r * x w s)) μ ∧
      ∫ w, x w p * x w q * (x w r * x w s) ∂μ = m4 p q r s)
    (hp : Pairing S m4) :
    ∫ w, (x w ⬝ᵥ (H *ᵥ x w)) ^ 2 ∂μ - (∫ w, x w ⬝ᵥ (H *ᵥ x w) ∂μ) ^ 2 =
      2 * ((H * S) * (H * S)).trace := by
  have hS : S.IsSymm := by
    ext p q
    rw [transpose_apply, ← (h2 q p).2, ← (h2 p q).2]
    congr 1
    funext w
    ring
  have hmean : ∫ w, x w ⬝ᵥ (H *ᵥ x w) ∂μ = ∑ p, ∑ q, H p q * S p q := by
    simp_rw [dotProduct_mulVec_eq_sum]
    rw [integral_finsetSum _ fun p _ => integrable_finsetSum _ fun q _ => (h2 p q).1.const_mul _]
    refine Finset.sum_congr rfl fun p _ => ?_
    rw [integral_finsetSum _ fun q _ => (h2 p q).1.const_mul _]
    refine Finset.sum_congr rfl fun q _ => ?_
    rw [integral_const_mul, (h2 p q).2]
  have hsq : ∀ w, (x w ⬝ᵥ (H *ᵥ x w)) ^ 2 =
      ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * (x w p * x w q * (x w r * x w s)) := by
    intro w
    rw [dotProduct_mulVec_eq_sum, sq, Finset.sum_mul]
    refine Finset.sum_congr rfl fun p _ => ?_
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun q _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun r _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun s _ => ?_
    ring
  have hsecond : ∫ w, (x w ⬝ᵥ (H *ᵥ x w)) ^ 2 ∂μ =
      ∑ p, ∑ q, ∑ r, ∑ s, H p q * H r s * m4 p q r s := by
    simp_rw [hsq]
    rw [integral_finsetSum _ fun p _ => integrable_finsetSum _ fun q _ =>
      integrable_finsetSum _ fun r _ => integrable_finsetSum _ fun s _ =>
        (h4 p q r s).1.const_mul _]
    refine Finset.sum_congr rfl fun p _ => ?_
    rw [integral_finsetSum _ fun q _ => integrable_finsetSum _ fun r _ =>
      integrable_finsetSum _ fun s _ => (h4 p q r s).1.const_mul _]
    refine Finset.sum_congr rfl fun q _ => ?_
    rw [integral_finsetSum _ fun r _ => integrable_finsetSum _ fun s _ =>
      (h4 p q r s).1.const_mul _]
    refine Finset.sum_congr rfl fun r _ => ?_
    rw [integral_finsetSum _ fun s _ => (h4 p q r s).1.const_mul _]
    refine Finset.sum_congr rfl fun s _ => ?_
    rw [integral_const_mul, (h4 p q r s).2]
  rw [hsecond, hmean]
  exact pairing_variance hH hS hp

/-! ### Witness -/

/-- The six-point law: `−√3` at `0`, `√3` at `5`, `0` elsewhere. -/
noncomputable def sixPoint (w : Fin 6) : Fin 1 → ℝ :=
  fun _ => if w = 0 then -√3 else if w = 5 then √3 else 0

/-- The pairing premise is met: under `(1/6)·count`, `E x² = 1` and `E x⁴ = 3 = 3·1²`, so the
quadratic form `x²` has variance `2·tr((1·1)²) = 2`. -/
example : ∫ w, (sixPoint w ⬝ᵥ ((1 : Matrix (Fin 1) (Fin 1) ℝ) *ᵥ sixPoint w)) ^ 2
      ∂((1 / 6 : NNReal) • Measure.count) -
    (∫ w, sixPoint w ⬝ᵥ ((1 : Matrix (Fin 1) (Fin 1) ℝ) *ᵥ sixPoint w)
      ∂((1 / 6 : NNReal) • Measure.count)) ^ 2 =
      2 * (((1 : Matrix (Fin 1) (Fin 1) ℝ) * 1) * (1 * 1)).trace := by
  have h3 : √3 ^ 2 = (3 : ℝ) := Real.sq_sqrt (by norm_num)
  refine variance_quadForm _ sixPoint (Matrix.isSymm_one) 1 (fun p q => ⟨Integrable.of_finite, ?_⟩)
    (fun _ _ _ _ => 3) (fun p q r s => ⟨Integrable.of_finite, ?_⟩) ?_
  · rw [integral_smul_nnreal_measure, integral_fintype Integrable.of_finite]
    fin_cases p
    fin_cases q
    simp [sixPoint, Fin.sum_univ_succ, Measure.real, NNReal.smul_def]
    norm_num
  · rw [integral_smul_nnreal_measure, integral_fintype Integrable.of_finite]
    simp [sixPoint, Fin.sum_univ_succ, Measure.real, NNReal.smul_def]
    norm_num
  · intro p q r s
    fin_cases p; fin_cases q; fin_cases r; fin_cases s
    norm_num

end TimeFrequency
