import Mathlib.Algebra.Order.Star.Real
import Mathlib.Analysis.Matrix.PosDef
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# Active subspaces — the projection mean over mixed ranks (AS-1)

The statements of `Protolemmata/BASS/active-subspace-lemmas.md`, AS-1. Projections are orthogonal
(symmetric and idempotent), of any rank; the weights are nonnegative and sum to 1.

* **(a)** `IsOrthProj.posSemidef`, `projMean_posSemidef`, `one_sub_projMean_posSemidef`,
  `trace_projMean`: `0 ≼ P̄ ≼ I` and `tr P̄ = Σ w_s k_s`.
* **(b)** `projLoss_eq`: `Σ w_s ‖Q − P_s‖²_F = tr Q + Σ w_s tr P_s − 2 tr(Q P̄)`.
  `trace_sub_two_trace_mul_ge` and `summaryProj_trace_sub`: over projections of every rank,
  `tr Q − 2 tr(Q P̄)` is at least `Σᵢ min 0 (1 − 2μᵢ)`, attained by the projection onto any span of
  eigenvectors holding every eigenvalue above ½ and none below. `summaryLoss_le` puts them together:
  the rank summary is least, and an eigenvalue of exactly ½ is indifferent.
* **(c)** `half_frobSq_cross_rank`: `½‖P_s − UUᵀ‖²_F = (tr P_s − k)/2 + Σᵢ (1 − cos²θᵢ)`, with the
  `cos²θᵢ` the eigenvalues of `UᵀP_sU`, which lie in `[0, 1]` (`compress_eigenvalues_mem`).
-/

open Matrix Finset

namespace ActiveSubspace

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- An orthogonal projection: symmetric and idempotent. -/
def IsOrthProj (P : Matrix n n ℝ) : Prop := Pᵀ = P ∧ P * P = P

/-- The squared Frobenius norm, entry by entry. -/
noncomputable def frobSq {m : Type*} [Fintype m] (A : Matrix m n ℝ) : ℝ := ∑ i, ∑ j, A i j ^ 2

omit [DecidableEq n] in
theorem frobSq_eq_trace {m : Type*} [Fintype m] (A : Matrix m n ℝ) : frobSq A = (Aᵀ * A).trace := by
  simp only [frobSq, Matrix.trace, Matrix.diag_apply, mul_apply, transpose_apply, sq]
  exact Finset.sum_comm

omit [DecidableEq n] in
theorem IsOrthProj.isHermitian {P : Matrix n n ℝ} (hP : IsOrthProj P) : P.IsHermitian := by
  rw [IsHermitian, conjTranspose_eq_transpose_of_trivial, hP.1]

theorem IsOrthProj.one_sub {P : Matrix n n ℝ} (hP : IsOrthProj P) : IsOrthProj (1 - P) := by
  refine ⟨by rw [transpose_sub, transpose_one, hP.1], ?_⟩
  rw [sub_mul, mul_sub, mul_sub, one_mul, mul_one, one_mul, hP.2]
  abel

omit [DecidableEq n] in
/-- A projection is positive semidefinite: `P = PᵀP`. -/
theorem IsOrthProj.posSemidef {P : Matrix n n ℝ} (hP : IsOrthProj P) : P.PosSemidef := by
  have h : P = Pᴴ * P := by rw [conjTranspose_eq_transpose_of_trivial, hP.1, hP.2]
  rw [h]
  exact posSemidef_conjTranspose_mul_self P

omit [DecidableEq n] in
/-- A projection's diagonal lies in `[0, 1]`: `Pᵢᵢ = Σⱼ Pᵢⱼ² ≥ Pᵢᵢ²`. -/
theorem IsOrthProj.diag_mem {P : Matrix n n ℝ} (hP : IsOrthProj P) (i : n) :
    0 ≤ P i i ∧ P i i ≤ 1 := by
  have hsq : P i i = ∑ j, P i j ^ 2 := by
    conv_lhs => rw [← hP.2]
    rw [mul_apply]
    refine Finset.sum_congr rfl fun j _ => ?_
    have : P j i = P i j := by
      have := congrFun (congrFun hP.1 i) j
      rwa [transpose_apply] at this
    rw [this, sq]
  have h0 : 0 ≤ P i i := hsq ▸ Finset.sum_nonneg fun j _ => sq_nonneg _
  have hle : P i i ^ 2 ≤ ∑ j, P i j ^ 2 :=
    Finset.single_le_sum (f := fun j => P i j ^ 2) (fun j _ => sq_nonneg _) (Finset.mem_univ i)
  rw [← hsq] at hle
  exact ⟨h0, by nlinarith⟩

omit [DecidableEq n] in
/-- `‖P − Q‖²_F = tr P + tr Q − 2 tr(PQ)` for orthogonal projections. -/
theorem IsOrthProj.frobSq_sub {P Q : Matrix n n ℝ} (hP : IsOrthProj P) (hQ : IsOrthProj Q) :
    frobSq (P - Q) = P.trace + Q.trace - 2 * (P * Q).trace := by
  rw [frobSq_eq_trace, transpose_sub, hP.1, hQ.1, sub_mul, mul_sub, mul_sub, hP.2, hQ.2,
    trace_sub, trace_sub, trace_sub, trace_mul_comm Q P]
  ring

section Mean

variable {ι : Type*} [Fintype ι]

omit [DecidableEq n] in
/-- **AS-1(a).** `P̄ = Σ w_s P_s` is positive semidefinite. -/
theorem projMean_posSemidef (P : ι → Matrix n n ℝ) (hP : ∀ s, IsOrthProj (P s)) (w : ι → ℝ)
    (hw : ∀ s, 0 ≤ w s) : (∑ s, w s • P s).PosSemidef :=
  posSemidef_sum _ fun s _ => (hP s).posSemidef.smul (hw s)

/-- **AS-1(a).** `I − P̄ = Σ w_s (I − P_s)` is positive semidefinite. -/
theorem one_sub_projMean_posSemidef (P : ι → Matrix n n ℝ) (hP : ∀ s, IsOrthProj (P s))
    (w : ι → ℝ) (hw : ∀ s, 0 ≤ w s) (hsum : ∑ s, w s = 1) : (1 - ∑ s, w s • P s).PosSemidef := by
  have h : 1 - ∑ s, w s • P s = ∑ s, w s • (1 - P s) := by
    simp only [smul_sub, Finset.sum_sub_distrib, ← Finset.sum_smul, hsum, one_smul]
  rw [h]
  exact projMean_posSemidef (fun s => 1 - P s) (fun s => (hP s).one_sub) w hw

omit [DecidableEq n] in
/-- **AS-1(a).** `tr P̄ = Σ w_s tr P_s`, the weighted mean rank. -/
theorem trace_projMean (P : ι → Matrix n n ℝ) (w : ι → ℝ) :
    (∑ s, w s • P s).trace = ∑ s, w s * (P s).trace := by
  rw [trace_sum]
  simp only [trace_smul, smul_eq_mul]

omit [DecidableEq n] in
/-- **AS-1(b), the loss.** For a projection `Q` of any rank,
`Σ w_s ‖Q − P_s‖²_F = tr Q + Σ w_s tr P_s − 2 tr(Q P̄)`. -/
theorem projLoss_eq (P : ι → Matrix n n ℝ) (hP : ∀ s, IsOrthProj (P s)) (w : ι → ℝ)
    (hsum : ∑ s, w s = 1) (Q : Matrix n n ℝ) (hQ : IsOrthProj Q) :
    ∑ s, w s * frobSq (Q - P s) =
      Q.trace + ∑ s, w s * (P s).trace - 2 * (Q * ∑ s, w s • P s).trace := by
  simp only [hQ.frobSq_sub (hP _), Matrix.mul_sum, trace_sum, Matrix.mul_smul, trace_smul,
    smul_eq_mul]
  have e : ∀ s, w s * (Q.trace + (P s).trace - 2 * (Q * P s).trace) =
      w s * Q.trace + w s * (P s).trace - 2 * (w s * (Q * P s).trace) := fun s => by ring
  rw [Finset.sum_congr rfl fun s _ => e s, Finset.sum_sub_distrib, Finset.sum_add_distrib,
    ← Finset.sum_mul, hsum, one_mul, ← Finset.mul_sum]

end Mean

section Spectral

variable {A : Matrix n n ℝ} (hA : A.IsHermitian)

/-- The eigenvector matrix of a real symmetric matrix, as a plain matrix. -/
noncomputable abbrev eigU : Matrix n n ℝ := (hA.eigenvectorUnitary : Matrix n n ℝ)

theorem eigU_transpose_mul : (eigU hA)ᵀ * eigU hA = 1 := by
  rw [← conjTranspose_eq_transpose_of_trivial, ← star_eq_conjTranspose]
  exact Unitary.coe_star_mul_self _

theorem eigU_mul_transpose : eigU hA * (eigU hA)ᵀ = 1 := by
  rw [← conjTranspose_eq_transpose_of_trivial, ← star_eq_conjTranspose]
  exact Unitary.coe_mul_star_self _

/-- The real spectral theorem: `A = U diag(μ) Uᵀ`. -/
theorem eq_eigU : A = eigU hA * diagonal hA.eigenvalues * (eigU hA)ᵀ := by
  conv_lhs => rw [hA.spectral_theorem]
  rw [Unitary.conjStarAlgAut_apply, star_eq_conjTranspose, conjTranspose_eq_transpose_of_trivial,
    RCLike.ofReal_real_eq_id, Function.id_comp]

/-- `tr Q − 2 tr(QA) = Σᵢ (1 − 2μᵢ)(UᵀQU)ᵢᵢ` for any `Q`. -/
theorem trace_sub_two_trace_mul_eq (Q : Matrix n n ℝ) :
    Q.trace - 2 * (Q * A).trace =
      ∑ i, (1 - 2 * hA.eigenvalues i) * ((eigU hA)ᵀ * Q * eigU hA) i i := by
  have h1 : Q.trace = ((eigU hA)ᵀ * Q * eigU hA).trace := by
    rw [Matrix.mul_assoc, trace_mul_comm, Matrix.mul_assoc, eigU_mul_transpose, Matrix.mul_one]
  have h2 : (Q * A).trace = ((eigU hA)ᵀ * Q * eigU hA * diagonal hA.eigenvalues).trace := by
    conv_lhs => rw [eq_eigU hA]
    rw [← Matrix.mul_assoc, ← Matrix.mul_assoc, trace_mul_comm, ← Matrix.mul_assoc,
      ← Matrix.mul_assoc]
  rw [h1, h2]
  simp only [Matrix.trace, Matrix.diag_apply, mul_diagonal]
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- **AS-1(b), the bound.** For a projection `Q` of any rank,
`tr Q − 2 tr(QA) ≥ Σᵢ min 0 (1 − 2μᵢ)`. -/
theorem trace_sub_two_trace_mul_ge (Q : Matrix n n ℝ) (hQ : IsOrthProj Q) :
    ∑ i, min 0 (1 - 2 * hA.eigenvalues i) ≤ Q.trace - 2 * (Q * A).trace := by
  have hR : IsOrthProj ((eigU hA)ᵀ * Q * eigU hA) := by
    refine ⟨by rw [transpose_mul, transpose_mul, transpose_transpose, hQ.1, Matrix.mul_assoc], ?_⟩
    calc (eigU hA)ᵀ * Q * eigU hA * ((eigU hA)ᵀ * Q * eigU hA)
        = (eigU hA)ᵀ * Q * (eigU hA * (eigU hA)ᵀ) * Q * eigU hA := by simp only [Matrix.mul_assoc]
      _ = (eigU hA)ᵀ * Q * eigU hA := by
          rw [eigU_mul_transpose, Matrix.mul_one, Matrix.mul_assoc ((eigU hA)ᵀ), hQ.2]
  rw [trace_sub_two_trace_mul_eq hA]
  refine Finset.sum_le_sum fun i _ => ?_
  obtain ⟨h0, h1⟩ := hR.diag_mem i
  rcases le_total 0 (1 - 2 * hA.eigenvalues i) with hc | hc
  · rw [min_eq_left hc]; exact mul_nonneg hc h0
  · rw [min_eq_right hc]; nlinarith

/-- The projection onto the span of the eigenvectors indexed by `S`. -/
noncomputable def summaryProj (S : Finset n) : Matrix n n ℝ :=
  eigU hA * diagonal (fun i => if i ∈ S then (1 : ℝ) else 0) * (eigU hA)ᵀ

theorem summaryProj_isOrthProj (S : Finset n) : IsOrthProj (summaryProj hA S) := by
  refine ⟨?_, ?_⟩
  · rw [summaryProj, transpose_mul, transpose_mul, transpose_transpose, diagonal_transpose,
      Matrix.mul_assoc]
  · simp only [summaryProj]
    calc eigU hA * diagonal (fun i => if i ∈ S then (1 : ℝ) else 0) * (eigU hA)ᵀ *
          (eigU hA * diagonal (fun i => if i ∈ S then (1 : ℝ) else 0) * (eigU hA)ᵀ)
        = eigU hA * diagonal (fun i => if i ∈ S then (1 : ℝ) else 0) * ((eigU hA)ᵀ * eigU hA) *
            diagonal (fun i => if i ∈ S then (1 : ℝ) else 0) * (eigU hA)ᵀ := by
          simp only [Matrix.mul_assoc]
      _ = _ := by
          rw [eigU_transpose_mul, Matrix.mul_one, Matrix.mul_assoc (eigU hA), diagonal_mul_diagonal]
          congr 3
          funext i
          split_ifs <;> norm_num

/-- **AS-1(b), attained.** For any `S` holding every eigenvalue above ½ and none below, the
projection onto those eigenvectors meets the bound: an eigenvalue of exactly ½ is indifferent. -/
theorem summaryProj_trace_sub (S : Finset n) (hS₁ : ∀ i, 1 / 2 < hA.eigenvalues i → i ∈ S)
    (hS₂ : ∀ i ∈ S, 1 / 2 ≤ hA.eigenvalues i) :
    (summaryProj hA S).trace - 2 * (summaryProj hA S * A).trace =
      ∑ i, min 0 (1 - 2 * hA.eigenvalues i) := by
  have hD : (eigU hA)ᵀ * summaryProj hA S * eigU hA =
      diagonal (fun i => if i ∈ S then (1 : ℝ) else 0) := by
    simp only [summaryProj, ← Matrix.mul_assoc, eigU_transpose_mul, Matrix.one_mul]
    rw [Matrix.mul_assoc, eigU_transpose_mul, Matrix.mul_one]
  rw [trace_sub_two_trace_mul_eq hA, hD]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp only [diagonal_apply_eq]
  by_cases hi : i ∈ S
  · rw [if_pos hi, mul_one, min_eq_right (by linarith [hS₂ i hi])]
  · rw [if_neg hi, mul_zero, min_eq_left (by linarith [not_lt.1 (mt (hS₁ i) hi)])]

/-- **AS-1(b), the rank summary is least.** Over projections `Q` of every rank, the expected
Frobenius loss to the draws is least at the projection onto `P̄`'s eigenvectors with eigenvalue
above ½ (with any of those at exactly ½ added or not). -/
theorem summaryLoss_le {ι : Type*} [Fintype ι] (P : ι → Matrix n n ℝ) (hP : ∀ s, IsOrthProj (P s))
    (w : ι → ℝ) (hsum : ∑ s, w s = 1) (hAdef : A = ∑ s, w s • P s) (S : Finset n)
    (hS₁ : ∀ i, 1 / 2 < hA.eigenvalues i → i ∈ S) (hS₂ : ∀ i ∈ S, 1 / 2 ≤ hA.eigenvalues i)
    (Q : Matrix n n ℝ) (hQ : IsOrthProj Q) :
    ∑ s, w s * frobSq (summaryProj hA S - P s) ≤ ∑ s, w s * frobSq (Q - P s) := by
  rw [projLoss_eq P hP w hsum _ (summaryProj_isOrthProj hA S), projLoss_eq P hP w hsum Q hQ,
    ← hAdef]
  have h1 := summaryProj_trace_sub hA S hS₁ hS₂
  have h2 := trace_sub_two_trace_mul_ge hA Q hQ
  linarith

end Spectral

section CrossRank

variable {k : ℕ}

omit [DecidableEq n] in
theorem IsOrthProj.isHermitian_compress {P : Matrix n n ℝ} (hP : IsOrthProj P)
    (U : Matrix n (Fin k) ℝ) : (Uᵀ * P * U).IsHermitian := by
  rw [IsHermitian, conjTranspose_eq_transpose_of_trivial, transpose_mul, transpose_mul,
    transpose_transpose, hP.1, Matrix.mul_assoc]

omit [DecidableEq n] in
/-- The squared cosines of the principal angles lie in `[0, 1]`. -/
theorem compress_eigenvalues_mem {P : Matrix n n ℝ} (hP : IsOrthProj P) (U : Matrix n (Fin k) ℝ)
    (hU : Uᵀ * U = 1) (i : Fin k) :
    0 ≤ (hP.isHermitian_compress U).eigenvalues i ∧
      (hP.isHermitian_compress U).eigenvalues i ≤ 1 := by
  classical
  have hpsd : (Uᵀ * P * U).PosSemidef := by
    have := hP.posSemidef.conjTranspose_mul_mul_same U
    rwa [conjTranspose_eq_transpose_of_trivial] at this
  have hpsd' : (1 - Uᵀ * P * U).PosSemidef := by
    have h : 1 - Uᵀ * P * U = Uᵀ * (1 - P) * U := by
      rw [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_one, hU]
    rw [h]
    have := hP.one_sub.posSemidef.conjTranspose_mul_mul_same U
    rwa [conjTranspose_eq_transpose_of_trivial] at this
  refine ⟨hpsd.eigenvalues_nonneg i, ?_⟩
  set hM := hP.isHermitian_compress U
  set v : Fin k → ℝ := ⇑(hM.eigenvectorBasis i)
  have hev : hM.eigenvalues i = v ⬝ᵥ ((Uᵀ * P * U) *ᵥ v) := by
    rw [hM.eigenvalues_eq i]
    simp [v]
  have hvv : v ⬝ᵥ v = 1 := by
    have h := hM.eigenvectorBasis.orthonormal.1 i
    have h2 : (inner ℝ (hM.eigenvectorBasis i) (hM.eigenvectorBasis i) : ℝ) = 1 := by
      rw [real_inner_self_eq_norm_sq, h, one_pow]
    rw [EuclideanSpace.inner_eq_star_dotProduct] at h2
    simpa [v, dotProduct_comm] using h2
  have hnn := hpsd'.dotProduct_mulVec_nonneg v
  simp only [star_trivial, sub_mulVec, one_mulVec, dotProduct_sub, hvv] at hnn
  linarith

omit [DecidableEq n] in
/-- **AS-1(c), the cross-rank distance.** For a projection `P_s` and an orthonormal `U` of `k`
columns, `½‖P_s − UUᵀ‖²_F = (tr P_s − k)/2 + Σᵢ (1 − cos²θᵢ)`, with `cos²θᵢ` the eigenvalues of
`UᵀP_sU`. At `k ≤ rank P_s` these are the `min(k_s, k)` principal angles' squared cosines, and the
right side is the rank term plus the containment defect. -/
theorem half_frobSq_cross_rank {P : Matrix n n ℝ} (hP : IsOrthProj P) (U : Matrix n (Fin k) ℝ)
    (hU : Uᵀ * U = 1) :
    (1 / 2) * frobSq (P - U * Uᵀ) =
      (P.trace - k) / 2 + ∑ i, (1 - (hP.isHermitian_compress U).eigenvalues i) := by
  have hQ : IsOrthProj (U * Uᵀ) := by
    refine ⟨by rw [transpose_mul, transpose_transpose], ?_⟩
    rw [Matrix.mul_assoc, ← Matrix.mul_assoc Uᵀ, hU, Matrix.one_mul]
  have htrQ : (U * Uᵀ).trace = k := by
    rw [trace_mul_comm, hU, trace_one, Fintype.card_fin]
  have htrPQ : (P * (U * Uᵀ)).trace = ∑ i, (hP.isHermitian_compress U).eigenvalues i := by
    have h := (hP.isHermitian_compress U).trace_eq_sum_eigenvalues
    simp only [RCLike.ofReal_real_eq_id, id] at h
    rw [← h, ← Matrix.mul_assoc, trace_mul_comm, ← Matrix.mul_assoc]
  rw [hP.frobSq_sub hQ, htrQ, htrPQ, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ,
    Fintype.card_fin, nsmul_eq_mul, mul_one]
  ring

end CrossRank

end ActiveSubspace
