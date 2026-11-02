import Mathlib.Algebra.Field.ZMod
import Mathlib.LinearAlgebra.Matrix.Block

/-!
# Persistence relative to a subcomplex that comes first

Take a filtered boundary matrix `D` whose first cells form a subcomplex `K₀`, so `K₀`'s cells
precede every other cell and their faces stay in `K₀`. The persistence of the pair `H(K_δ, K₀)`
reduces the matrix of the remaining cells, with `K₀`'s rows and columns deleted. That matrix needs
no second reduction: any reduction of `D` restricts to one of it, so one absolute run determines
the relative barcode.

| Absolute bar, `K₀` first | Relative bar |
|---|---|
| Born outside `K₀` | The same bar |
| Born in `K₀`, killed outside it at `j` | Essential, born at `j`, one degree up |
| Born in `K₀`, killed within it, or essential | None |

`relative_barcode` is the statement, for any reduction of the relative matrix and any of `D`, and
`low_degree` gives the degree. The last row holds because its cells are not in the relative index.

This is the analogue, for a fixed `K₀`, of Prop 2.4 of de Silva, Morozov and Vejdemo-Johansson's
*Dualities in persistent (co)homology*, whose relative module `H(Xₙ, Xᵢ)` moves with the
filtration. Stated for any field; the engine's case is `ZMod 2`. A bar is a pair `(i, j)` with
`IsLow R i j`, and `isLow_iff_of_reductions` makes that independent of the reduction, so the
statements hold for whichever one an algorithm produces. The code this covers is
`TDA.Ph.RelativeFiltration`, checked empirically by
`RelativeBarcode_IsARereadingOfTheAbsoluteRunWithK0First`.
-/

namespace RelativePersistence

open Matrix

variable {ι F : Type*} [LinearOrder ι] [Field F]

/-- Row `i` holds the last nonzero entry of column `j`: cell `i` is born and cell `j` kills it. -/
def IsLow (M : Matrix ι ι F) (i j : ι) : Prop :=
  M i j ≠ 0 ∧ ∀ r, i < r → M r j = 0

/-- No two columns share a low. -/
def Reduced (M : Matrix ι ι F) : Prop :=
  ∀ i j₁ j₂, IsLow M i j₁ → IsLow M i j₂ → j₁ = j₂

/-- Cell `i` starts a class that nothing kills: its column is zero and it is no column's low. -/
def Essential (M : Matrix ι ι F) (i : ι) : Prop :=
  (∀ r, M r i = 0) ∧ ∀ j, ¬ IsLow M i j

/-- `K` is a subcomplex of `D` that comes strictly first: it is an initial segment of the order, and
the faces of its cells stay in it. -/
structure IsFirstSubcomplex (D : Matrix ι ι F) (K : ι → Prop) : Prop where
  initial : ∀ a b, a ≤ b → K b → K a
  closed : ∀ i l, K l → ¬ K i → D i l = 0

/-- The matrix of the cells outside `K`: `K`'s rows and columns deleted. -/
def rel (K : ι → Prop) (M : Matrix ι ι F) : Matrix {x // ¬ K x} {x // ¬ K x} F :=
  M.submatrix Subtype.val Subtype.val

theorem IsLow.unique {M : Matrix ι ι F} {i i' j : ι} (h : IsLow M i j) (h' : IsLow M i' j) :
    i = i' := by
  rcases lt_trichotomy i i' with hlt | heq | hgt
  · exact absurd (h.2 i' hlt) h'.1
  · exact heq
  · exact absurd (h'.2 i hgt) h.1

/-- A nonzero column has a low. -/
theorem exists_isLow [Finite ι] {M : Matrix ι ι F} {j : ι} (h : ∃ r, M r j ≠ 0) :
    ∃ i, IsLow M i j := by
  classical
  cases nonempty_fintype ι
  obtain ⟨r₀, hr₀⟩ := h
  let s : Finset ι := Finset.univ.filter fun r => M r j ≠ 0
  have hs : s.Nonempty := ⟨r₀, by simp [s, hr₀]⟩
  refine ⟨s.max' hs, ?_, fun r hr => ?_⟩
  · simpa [s] using s.max'_mem hs
  · by_contra hne
    exact absurd (s.le_max' r (by simp [s, hne])) (not_le.mpr hr)

variable [Fintype ι]

/-- A reduction of `D`: `D = R * U` with `R` reduced and `U` upper triangular and invertible. -/
structure IsReduction (D R U : Matrix ι ι F) : Prop where
  factor : D = R * U
  reduced : Reduced R
  upper : ∀ i j, j < i → U i j = 0
  diag : ∀ i, U i i ≠ 0

/-- A combination of columns of a reduced matrix, using at least one nonzero column, is nonzero,
and its last nonzero entry sits at the low of a column it uses. Lows being distinct, the largest
one among the columns used cannot cancel. -/
theorem combination_top {R : Matrix ι ι F} (hR : Reduced R) (y : ι → F)
    (hy : ∃ l, y l ≠ 0 ∧ ∃ r, R r l ≠ 0) :
    ∃ i l, y l ≠ 0 ∧ IsLow R i l ∧ (∑ s, R i s * y s) ≠ 0 ∧
      ∀ r, i < r → ∑ s, R r s * y s = 0 := by
  classical
  obtain ⟨l₀, hl₀, hc₀⟩ := hy
  obtain ⟨i₀, hi₀⟩ := exists_isLow hc₀
  let T : Finset ι := Finset.univ.filter fun i => ∃ l, y l ≠ 0 ∧ IsLow R i l
  have hT : T.Nonempty := ⟨i₀, by simpa [T] using ⟨l₀, hl₀, hi₀⟩⟩
  obtain ⟨l, hl, hlow⟩ : ∃ l, y l ≠ 0 ∧ IsLow R (T.max' hT) l := by
    simpa [T] using T.max'_mem hT
  have hother : ∀ s r, s ≠ l → T.max' hT ≤ r → R r s * y s = 0 := by
    intro s r hs hr
    by_cases hys : y s = 0
    · simp [hys]
    by_cases hcs : ∃ q, R q s ≠ 0
    · obtain ⟨q, hq⟩ := exists_isLow hcs
      have hqT : q ≤ T.max' hT := T.le_max' q (by simpa [T] using ⟨s, hys, hq⟩)
      have hne : q ≠ T.max' hT := fun h => hs (hR _ s l (h ▸ hq) hlow)
      simp [hq.2 r (lt_of_lt_of_le (lt_of_le_of_ne hqT hne) hr)]
    · push Not at hcs
      simp [hcs r]
  refine ⟨T.max' hT, l, hl, hlow, ?_, fun r hr => ?_⟩
  · rw [Finset.sum_eq_single l (fun s _ hs => hother s _ hs le_rfl) (by simp)]
    exact mul_ne_zero hlow.1 hl
  · refine Finset.sum_eq_zero fun s _ => ?_
    by_cases hs : s = l
    · subst hs
      simp [hlow.2 r hr]
    · exact hother s r hs hr.le

/-- An upper triangular matrix with a nonzero diagonal is invertible. -/
theorem isUnit_det {U : Matrix ι ι F} (hu : ∀ i j, j < i → U i j = 0) (hd : ∀ i, U i i ≠ 0) :
    IsUnit U.det := by
  have hb : U.BlockTriangular id := fun i j hij => hu i j hij
  rw [Matrix.det_of_upperTriangular hb]
  exact isUnit_iff_ne_zero.mpr (Finset.prod_ne_zero_iff.mpr fun i _ => hd i)

/-- A cell that is some column's low has a zero column: no cell both kills a class and is born.
This is where `D * D = 0` enters. -/
theorem column_zero_of_isLow {D R U : Matrix ι ι F} (hD : D * D = 0) (h : IsReduction D R U)
    {i j : ι} (hlow : IsLow R i j) : ∀ r, R r i = 0 := by
  classical
  have hdet := isUnit_det h.upper h.diag
  have hR : R = D * U⁻¹ := by
    rw [h.factor, Matrix.mul_assoc, Matrix.mul_nonsing_inv U hdet, Matrix.mul_one]
  have h1 : D * R = 0 := by
    calc D * R = D * (D * U⁻¹) := by rw [← hR]
      _ = 0 := by rw [← Matrix.mul_assoc, hD, Matrix.zero_mul]
  have hzero : R * (U * R) = 0 := by
    rw [← Matrix.mul_assoc, ← h.factor, h1]
  have hyi : (U * R) i j ≠ 0 := by
    rw [Matrix.mul_apply, Finset.sum_eq_single i]
    · exact mul_ne_zero (h.diag i) hlow.1
    · intro s _ hs
      rcases lt_or_gt_of_ne hs with hlt | hgt
      · simp [h.upper i s hlt]
      · simp [hlow.2 s hgt]
    · simp
  by_contra hne
  push Not at hne
  obtain ⟨r, hr⟩ := hne
  obtain ⟨q, _, _, _, hq, _⟩ :=
    combination_top h.reduced (fun l => (U * R) l j) ⟨i, hyi, r, hr⟩
  apply hq
  have := congrFun (congrFun hzero q) j
  simpa [Matrix.mul_apply] using this

/-- `R`'s columns inherit `D`'s closure: a column in `K` has no entry outside `K`. -/
theorem reduction_closed {D R U : Matrix ι ι F} {K : ι → Prop} (hK : IsFirstSubcomplex D K)
    (h : IsReduction D R U) (i : ι) (hi : ¬ K i) (l : ι) (hl : K l) : R i l = 0 := by
  induction l using WellFoundedLT.induction with
  | ind l ih =>
    have hsum : D i l = R i l * U l l := by
      rw [h.factor, Matrix.mul_apply, Finset.sum_eq_single l]
      · intro s _ hs
        rcases lt_or_gt_of_ne hs with hlt | hgt
        · simp [ih s hlt (hK.initial s l hlt.le hl)]
        · simp [h.upper s l hgt]
      · simp
    have := hK.closed i l hl hi
    rw [hsum] at this
    exact (mul_eq_zero.mp this).resolve_right (h.diag l)

/-- A bar born outside `K` is killed outside `K`. -/
theorem not_mem_of_isLow {D R U : Matrix ι ι F} {K : ι → Prop} (hK : IsFirstSubcomplex D K)
    (h : IsReduction D R U) {i j : ι} (hi : ¬ K i) (hlow : IsLow R i j) : ¬ K j :=
  fun hj => hlow.1 (reduction_closed hK h i hi j hj)

omit [Fintype ι] in
/-- The relative matrix's lows are exactly the absolute lows outside `K`: finite bars born outside
`K` carry over unchanged, and nothing else appears. -/
theorem relative_isLow_iff {D : Matrix ι ι F} {R : Matrix ι ι F} {K : ι → Prop}
    (hK : IsFirstSubcomplex D K) (i j : {x // ¬ K x}) :
    IsLow (rel K R) i j ↔ IsLow R i.1 j.1 := by
  constructor
  · refine fun h => ⟨h.1, fun r hr => ?_⟩
    have hrK : ¬ K r := fun hk => i.2 (hK.initial _ _ hr.le hk)
    exact h.2 ⟨r, hrK⟩ hr
  · exact fun h => ⟨h.1, fun r hr => h.2 r.1 hr⟩

variable {K : ι → Prop}

/-- Any reduction of `D` restricts to a reduction of the relative matrix, so the relative barcode
needs no second reduction. -/
theorem relative_isReduction [DecidablePred K] {D R U : Matrix ι ι F}
    (hK : IsFirstSubcomplex D K) (h : IsReduction D R U) :
    IsReduction (rel K D) (rel K R) (rel K U) where
  factor := by
    ext i j
    have hsplit := Fintype.sum_subtype_add_sum_subtype K (fun s => R i.1 s * U s j.1)
    have hin : ∑ s : {x // K x}, R i.1 s.1 * U s.1 j.1 = 0 :=
      Finset.sum_eq_zero fun s _ => by simp [reduction_closed hK h i.1 i.2 s.1 s.2]
    rw [hin, zero_add] at hsplit
    simp only [rel, Matrix.submatrix_apply, Matrix.mul_apply]
    rw [h.factor, Matrix.mul_apply, hsplit]
  reduced := fun i j₁ j₂ h₁ h₂ =>
    Subtype.ext (h.reduced i.1 j₁.1 j₂.1 ((relative_isLow_iff hK i j₁).1 h₁)
      ((relative_isLow_iff hK i j₂).1 h₂))
  upper := fun i j hij => h.upper i.1 j.1 hij
  diag := fun i => h.diag i.1

omit [Fintype ι] in
/-- Every cell in `K` comes before every cell outside it. -/
theorem lt_of_mem_of_not_mem {D : Matrix ι ι F} (hK : IsFirstSubcomplex D K) {l r : ι} (hl : K l)
    (hr : ¬ K r) : l < r := by
  by_contra hle
  exact hr (hK.initial _ _ (not_lt.mp hle) hl)

/-- The relative essential classes are the absolute ones outside `K`, together with one for each
cell outside `K` that kills a class born in `K`. -/
theorem relative_essential_iff {D R U : Matrix ι ι F} (hD : D * D = 0)
    (hK : IsFirstSubcomplex D K) (h : IsReduction D R U) (i : {x // ¬ K x}) :
    Essential (rel K R) i ↔ Essential R i.1 ∨ ∃ l, K l ∧ IsLow R l i.1 := by
  constructor
  · rintro ⟨hz, hn⟩
    by_cases hc : ∃ r, R r i.1 ≠ 0
    · obtain ⟨l, hl⟩ := exists_isLow hc
      by_cases hlK : K l
      · exact Or.inr ⟨l, hlK, hl⟩
      · exact absurd (hz ⟨l, hlK⟩) hl.1
    · push Not at hc
      refine Or.inl ⟨hc, fun j hj => ?_⟩
      exact hn ⟨j, not_mem_of_isLow hK h i.2 hj⟩ ((relative_isLow_iff hK i _).2 hj)
  · rintro (⟨hz, hn⟩ | ⟨l, hlK, hl⟩)
    · exact ⟨fun r => hz r.1, fun j hj => hn j.1 ((relative_isLow_iff hK i j).1 hj)⟩
    · refine ⟨fun r => hl.2 r.1 (lt_of_mem_of_not_mem hK hlK r.2), fun j hj => ?_⟩
      · exact hl.1 (column_zero_of_isLow hD h ((relative_isLow_iff hK i j).1 hj) l)

/-- Reductions of graded matrices stay graded, so every bar `(i, j)` has `j` one degree above `i`.
A bar born in `K` and killed outside it therefore moves up one degree in the relative barcode. -/
theorem low_degree {D R U : Matrix ι ι F} (h : IsReduction D R U) (dim : ι → ℕ)
    (hDg : ∀ i j, D i j ≠ 0 → dim j = dim i + 1) (hUg : ∀ i j, U i j ≠ 0 → dim i = dim j)
    (i j : ι) (hij : R i j ≠ 0) : dim j = dim i + 1 := by
  classical
  induction j using WellFoundedLT.induction with
  | ind j ih =>
    by_cases hd : D i j = 0
    · have hsum : ∑ s, R i s * U s j = 0 := by
        rw [← Matrix.mul_apply, ← h.factor, hd]
      rw [← Finset.add_sum_erase _ _ (Finset.mem_univ j)] at hsum
      have hrest : ∑ s ∈ Finset.univ.erase j, R i s * U s j ≠ 0 := by
        intro hzero
        rw [hzero, add_zero] at hsum
        exact mul_ne_zero hij (h.diag j) hsum
      obtain ⟨s, hs, hne⟩ := Finset.exists_ne_zero_of_sum_ne_zero hrest
      have hsj : s ≠ j := (Finset.mem_erase.mp hs).1
      have hRs : R i s ≠ 0 := left_ne_zero_of_mul hne
      have hUs : U s j ≠ 0 := right_ne_zero_of_mul hne
      have hlt : s < j := by
        rcases lt_or_gt_of_ne hsj with hlt | hgt
        · exact hlt
        · exact absurd (h.upper s j hgt) hUs
      rw [← hUg s j hUs]
      exact ih s hlt hRs
    · exact hDg i j hd

/-- If `R₁ = R₂ * T` with `T` upper triangular and `R₂` reduced, each low of `R₁`'s column `j` is
the low of some column of `R₂` at or before `j`. -/
theorem exists_isLow_of_factor {R₁ R₂ T : Matrix ι ι F} (hR : Reduced R₂) (hf : R₁ = R₂ * T)
    (hT : ∀ i j, j < i → T i j = 0) {i j : ι} (hlow : IsLow R₁ i j) :
    ∃ l, l ≤ j ∧ IsLow R₂ i l := by
  classical
  have hcol : ∀ r, R₁ r j = ∑ s, R₂ r s * T s j := fun r => by rw [hf, Matrix.mul_apply]
  have hy : ∃ l, T l j ≠ 0 ∧ ∃ r, R₂ r l ≠ 0 := by
    by_contra hnone
    push Not at hnone
    apply hlow.1
    rw [hcol]
    refine Finset.sum_eq_zero fun s _ => ?_
    by_cases hs : T s j = 0
    · simp [hs]
    · simp [hnone s hs i]
  obtain ⟨q, l, hl, hql, hne, htop⟩ := combination_top hR (fun s => T s j) hy
  have hq : IsLow R₁ q j :=
    ⟨by rw [hcol]; exact hne, fun r hr => by rw [hcol]; exact htop r hr⟩
  have : i = q := hlow.unique hq
  subst this
  refine ⟨l, ?_, hql⟩
  by_contra hgt
  exact hl (hT l j (not_le.mp hgt))

/-- Two reductions of one matrix have the same lows, so a barcode belongs to the matrix and not to
the algorithm that reduced it. -/
theorem isLow_iff_of_reductions {D R₁ U₁ R₂ U₂ : Matrix ι ι F} (h₁ : IsReduction D R₁ U₁)
    (h₂ : IsReduction D R₂ U₂) (i j : ι) : IsLow R₁ i j ↔ IsLow R₂ i j := by
  classical
  have hdet₁ := isUnit_det h₁.upper h₁.diag
  have hdet₂ := isUnit_det h₂.upper h₂.diag
  have hb₁ : U₁.BlockTriangular id := fun a b hab => h₁.upper a b hab
  have hb₂ : U₂.BlockTriangular id := fun a b hab => h₂.upper a b hab
  letI := Matrix.invertibleOfIsUnitDet U₁ hdet₁
  letI := Matrix.invertibleOfIsUnitDet U₂ hdet₂
  have hT : (U₂ * U₁⁻¹).BlockTriangular id :=
    hb₂.mul (Matrix.blockTriangular_inv_of_blockTriangular hb₁)
  have hS : (U₁ * U₂⁻¹).BlockTriangular id :=
    hb₁.mul (Matrix.blockTriangular_inv_of_blockTriangular hb₂)
  have hf₁ : R₁ = R₂ * (U₂ * U₁⁻¹) := by
    rw [← Matrix.mul_assoc, ← h₂.factor, h₁.factor, Matrix.mul_assoc,
      Matrix.mul_nonsing_inv U₁ hdet₁, Matrix.mul_one]
  have hf₂ : R₂ = R₁ * (U₁ * U₂⁻¹) := by
    rw [← Matrix.mul_assoc, ← h₁.factor, h₂.factor, Matrix.mul_assoc,
      Matrix.mul_nonsing_inv U₂ hdet₂, Matrix.mul_one]
  induction j using WellFoundedLT.induction generalizing i with
  | ind j ih =>
    constructor
    · intro hlow
      obtain ⟨l, hle, hl⟩ := exists_isLow_of_factor h₂.reduced hf₁ (fun a b hab => hT hab) hlow
      rcases lt_or_eq_of_le hle with hlt | rfl
      · exact absurd (h₁.reduced i l j ((ih l hlt i).2 hl) hlow) hlt.ne
      · exact hl
    · intro hlow
      obtain ⟨l, hle, hl⟩ := exists_isLow_of_factor h₁.reduced hf₂ (fun a b hab => hS hab) hlow
      rcases lt_or_eq_of_le hle with hlt | rfl
      · exact absurd (h₂.reduced i l j ((ih l hlt i).1 hl) hlow) hlt.ne
      · exact hl

omit [Fintype ι] in
/-- A column is zero exactly when it has no low, so the lows alone decide which classes are
essential. -/
theorem essential_iff_of_isLow_iff [Finite ι] {M N : Matrix ι ι F}
    (hMN : ∀ i j, IsLow M i j ↔ IsLow N i j) (i : ι) : Essential M i ↔ Essential N i := by
  have hzero : ∀ {A : Matrix ι ι F}, (∀ r, A r i = 0) ↔ ∀ r, ¬ IsLow A r i := by
    intro A
    constructor
    · exact fun h r hl => hl.1 (h r)
    · intro h r
      by_contra hne
      obtain ⟨q, hq⟩ := exists_isLow ⟨r, hne⟩
      exact h q hq
  unfold Essential
  rw [hzero, hzero]
  simp only [hMN]

/-- One absolute run with `K` first determines the relative barcode. For any reduction `R` of `D`
and any reduction `R'` of the relative matrix, the relative bars are the absolute bars born outside
`K`, and the relative essential classes are the absolute ones outside `K` together with one for
each cell that kills a class born in `K`. -/
theorem relative_barcode [DecidablePred K] {D R U : Matrix ι ι F}
    {R' U' : Matrix {x // ¬ K x} {x // ¬ K x} F} (hD : D * D = 0) (hK : IsFirstSubcomplex D K)
    (h : IsReduction D R U) (h' : IsReduction (rel K D) R' U') (i j : {x // ¬ K x}) :
    (IsLow R' i j ↔ IsLow R i.1 j.1) ∧
      (Essential R' i ↔ Essential R i.1 ∨ ∃ l, K l ∧ IsLow R l i.1) := by
  have hlow := isLow_iff_of_reductions h' (relative_isReduction hK h)
  exact ⟨(hlow i j).trans (relative_isLow_iff hK i j),
    (essential_iff_of_isLow_iff hlow i).trans (relative_essential_iff hD hK h i)⟩

section Example

/-! The hypotheses can be met, and the case the file is about fires. Two vertices and the edge
joining them over `ZMod 2`, with `K₀` the two vertices: absolutely the edge kills vertex `1`'s
class, and relatively it is an essential class, one degree up. -/

/-- Two vertices, cells `0` and `1`, and the edge joining them, cell `2`. -/
private def edge : Matrix (Fin 3) (Fin 3) (ZMod 2) := !![0, 0, 1; 0, 0, 1; 0, 0, 0]

private theorem edge_isReduction : IsReduction edge edge 1 where
  factor := by simp
  reduced := by unfold Reduced IsLow; decide
  upper := fun _ _ hij => Matrix.one_apply_ne hij.ne'
  diag := fun i => by simp

private theorem edge_first : IsFirstSubcomplex edge (· ≤ (1 : Fin 3)) where
  initial := fun _ _ hab hb => hab.trans hb
  closed := by decide

private theorem edge_boundary : edge * edge = 0 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [edge, Matrix.mul_apply, Fin.sum_univ_succ]

example : Essential (rel (· ≤ (1 : Fin 3)) edge) ⟨2, by decide⟩ :=
  (relative_essential_iff edge_boundary edge_first edge_isReduction _).2
    (Or.inr ⟨1, by decide, by unfold IsLow; decide⟩)

end Example

end RelativePersistence
