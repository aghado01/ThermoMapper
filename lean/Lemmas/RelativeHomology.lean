import Lemmas.RelativePersistence
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.ToLin
import Mathlib.LinearAlgebra.Quotient.Basic

/-!
# The relative barcode is the persistent homology of the pair

`RelativePersistence` proves a table about matrices, taking a bar to be whatever reduction reads
off. This file proves what that reading means. It defines the relative homology of a filtered pair
directly, `H_p(S, K) = Z_p(S, K) / B_p(S, K)`, with no reduction in sight, and shows that the rank
of `H_p(S, K) → H_p(T, K)` counts exactly the bars of the table alive from `S` through `T`
(`finrank_persistence`, read by `alive_iff_table`). With `K` empty this is the classical statement
that reduction computes persistent homology (`finrank_persistence_absolute`).

Two more results complete the picture. `exists_gradedReduction` shows that every filtered complex
has a reduction, so the statements are never vacuous. `finrank_quotient_persistence` shows that the
complex with `K`'s cells deleted, which is what `TDA.Ph.RelativeFiltration` reduces, has the ranks
of the pair, with no reduction in the statement. `SimplicialFiltration` instantiates all of this
for ℤ/2 simplicial chains, the engine's case.

The result here is stated in ranks. `IntervalDecomposition` proves the module is the sum of those
intervals, and that the ranks determine them. The C# is not covered: `TDA.Ph.RelativeFiltration` is
checked by tests. The statement is folklore in substance, not a new result: see
`issues/ph/relative-barcode-provenance.md`.
-/

namespace RelativePersistence

open Module

section Counts

variable {α β : Type*}

/-- The number of cells satisfying `A`. -/
noncomputable abbrev cellCount (A : α → Prop) : ℕ := Nat.card {i // A i}

theorem cellCount_congr {A B : α → Prop} (h : ∀ i, A i ↔ B i) : cellCount A = cellCount B :=
  Nat.card_congr (Equiv.subtypeEquivRight h)

theorem cellCount_bij {A : α → Prop} {B : β → Prop} (f : α → β) (hmaps : ∀ i, A i → B (f i))
    (hinj : ∀ i j, A i → A j → f i = f j → i = j) (hsurj : ∀ b, B b → ∃ a, A a ∧ f a = b) :
    cellCount A = cellCount B := by
  refine Nat.card_congr (Equiv.ofBijective (fun x => ⟨f x.1, hmaps _ x.2⟩) ⟨?_, ?_⟩)
  · intro x y hxy
    exact Subtype.ext (hinj _ _ x.2 y.2 (congrArg Subtype.val hxy))
  · intro y
    obtain ⟨a, ha, hfa⟩ := hsurj y.1 y.2
    exact ⟨⟨a, ha⟩, Subtype.ext hfa⟩

theorem cellCount_eq_zero {A : α → Prop} (h : ∀ i, ¬ A i) : cellCount A = 0 := by
  haveI : IsEmpty {i // A i} := ⟨fun i => h i.1 i.2⟩
  exact Nat.card_of_isEmpty

/-- Cells split by any further condition. -/
theorem cellCount_split [Finite α] (A Q : α → Prop) :
    cellCount (fun i => A i ∧ Q i) + cellCount (fun i => A i ∧ ¬ Q i) = cellCount A := by
  classical
  cases nonempty_fintype α
  simp only [cellCount, Nat.card_eq_fintype_card, Fintype.card_subtype]
  rw [← Finset.card_filter_add_card_filter_not (s := Finset.univ.filter A) Q]
  congr 2 <;> (ext i; simp)

end Counts

section Chains

variable {ι : Type*} (F : Type*) [Field F]

/-- Chains supported on the cells satisfying `A`. -/
def chains (A : ι → Prop) : Submodule F (ι → F) where
  carrier := {x | ∀ i, ¬ A i → x i = 0}
  add_mem' := fun hx hy i hi => by simp [hx i hi, hy i hi]
  zero_mem' := fun _ _ => rfl
  smul_mem' := fun c _ hx i hi => by simp [hx i hi]

variable {F}

theorem chains_mono {A B : ι → Prop} (h : ∀ i, A i → B i) : chains F A ≤ chains F B :=
  fun _ hx i hi => hx i (mt (h i) hi)

/-- The chains on a set of cells have one dimension per cell. -/
theorem finrank_chains [Finite ι] (A : ι → Prop) : finrank F (chains F A) = cellCount A := by
  classical
  cases nonempty_fintype ι
  let e : chains F A ≃ₗ[F] ({i // A i} → F) :=
    { toFun := fun x i => x.1 i.1
      map_add' := fun _ _ => rfl
      map_smul' := fun _ _ => rfl
      invFun := fun y => ⟨fun i => if h : A i then y ⟨i, h⟩ else 0, fun i hi => by simp [hi]⟩
      left_inv := fun x => Subtype.ext (funext fun i => by
        by_cases h : A i
        · simp [h]
        · simpa [h] using (x.2 i h).symm)
      right_inv := fun y => funext fun i => by simp [i.2] }
  rw [e.finrank_eq, finrank_fintype_fun_eq_card, cellCount, Nat.card_eq_fintype_card]

end Chains

section Low

variable {ι F : Type*} [LinearOrder ι] [Field F]

open scoped Classical in
/-- The low of a column, where it has one. -/
noncomputable def lowOf (R : Matrix ι ι F) (j : ι) : ι :=
  if h : ∃ i, IsLow R i j then h.choose else j

theorem lowOf_eq {R : Matrix ι ι F} {i j : ι} (h : IsLow R i j) : lowOf R j = i := by
  have hex : ∃ i, IsLow R i j := ⟨i, h⟩
  rw [lowOf, dif_pos hex]
  exact hex.choose_spec.unique h

end Low

variable {ι F : Type*} [LinearOrder ι] [Fintype ι] [Field F]

/-- A set of cells is an initial segment of the filtration order. -/
def IsInitial (S : ι → Prop) : Prop := ∀ a b, a ≤ b → S b → S a

/-- A based chain complex with its cells in filtration order. Every boundary is made of earlier
cells, so every initial segment is a subcomplex. -/
structure FilteredComplex (ι F : Type*) [LinearOrder ι] [Fintype ι] [Field F] where
  /-- The boundary map in the basis of cells. -/
  boundary : Matrix ι ι F
  /-- The degree of a cell. -/
  deg : ι → ℕ
  boundary_sq : boundary * boundary = 0
  boundary_lt : ∀ i j, boundary i j ≠ 0 → i < j
  boundary_deg : ∀ i j, boundary i j ≠ 0 → deg j = deg i + 1

/-- Column `j` has its low outside `X`. -/
def LowOutside (R : Matrix ι ι F) (X : ι → Prop) (j : ι) : Prop := ∃ i, IsLow R i j ∧ ¬ X i

section Reduced

variable {R : Matrix ι ι F} {X : ι → Prop}

omit [Fintype ι] in
theorem not_lowOutside_iff {j : ι} : ¬ LowOutside R X j ↔ ∀ l, IsLow R l j → X l := by
  simp [LowOutside]

omit [Fintype ι] in
/-- A column whose low is in an initial segment lies in that segment. -/
theorem column_eq_zero [Finite ι] (hX : IsInitial X) {j : ι} (hj : ¬ LowOutside R X j) (r : ι)
    (hr : ¬ X r) : R r j = 0 := by
  by_contra hne
  obtain ⟨q, hq⟩ := exists_isLow ⟨r, hne⟩
  have hrq : r ≤ q := by
    by_contra hlt
    exact hne (hq.2 r (not_le.mp hlt))
  exact hj ⟨q, hq, fun hqX => hr (hX r q hrq hqX)⟩

/-- A combination of columns reaches outside an initial segment exactly when it uses a column
whose low is outside it: lows being distinct, the largest cannot cancel. -/
theorem mulVec_mem_chains_iff (hR : Reduced R) (hX : IsInitial X) {c : ι → F} :
    R.mulVec c ∈ chains F X ↔ ∀ j, c j ≠ 0 → ¬ LowOutside R X j := by
  constructor
  · intro hm
    by_contra hbad
    push Not at hbad
    classical
    let c₁ : ι → F := fun j => if LowOutside R X j then c j else 0
    have hc₂ : R.mulVec (c - c₁) ∈ chains F X := by
      intro r hr
      refine Finset.sum_eq_zero fun s _ => ?_
      by_cases hl : LowOutside R X s
      · simp [c₁, hl]
      · simp [column_eq_zero hX hl r hr]
    have hc₁ : R.mulVec c₁ ∈ chains F X := by
      have hsplit : R.mulVec c₁ = R.mulVec c - R.mulVec (c - c₁) := by
        rw [Matrix.mulVec_sub, sub_sub_cancel]
      rw [hsplit]
      exact Submodule.sub_mem _ hm hc₂
    obtain ⟨l₀, hl₀, i₀, hi₀, hiX⟩ := hbad
    have hl₀O : LowOutside R X l₀ := ⟨i₀, hi₀, hiX⟩
    obtain ⟨i, l, hl, hlow, htop, _⟩ :=
      combination_top hR c₁ ⟨l₀, by simpa [c₁, hl₀O], i₀, hi₀.1⟩
    have hlO : LowOutside R X l := by
      by_contra hn
      exact hl (by simp [c₁, hn])
    obtain ⟨q, hq, hqX⟩ := hlO
    have hiq : i = q := hlow.unique hq
    subst hiq
    exact htop (hc₁ i hqX)
  · intro h r hr
    refine Finset.sum_eq_zero fun s _ => ?_
    by_cases hs : c s = 0
    · simp [hs]
    · simp [column_eq_zero hX (h s hs) r hr]

/-- A nonzero combination of columns whose lows are all outside `X` leaves `X`. -/
theorem eq_zero_of_mulVec_mem (hR : Reduced R) {c : ι → F} (hc : ∀ j, c j ≠ 0 → LowOutside R X j)
    (hm : R.mulVec c ∈ chains F X) : c = 0 := by
  by_contra hne
  obtain ⟨l₀, hl₀⟩ := Function.ne_iff.mp hne
  obtain ⟨i₀, hi₀, _⟩ := hc l₀ hl₀
  obtain ⟨i, l, hl, hlow, htop, _⟩ := combination_top hR c ⟨l₀, hl₀, i₀, hi₀.1⟩
  obtain ⟨q, hq, hqX⟩ := hc l hl
  have hiq : i = q := hlow.unique hq
  subst hiq
  exact htop (hm i hqX)

end Reduced

/-- An invertible matrix that keeps a set of cells closed maps its chains onto themselves. -/
theorem chains_map_eq {V : Matrix ι ι F} (hV : IsUnit V.det) {B : ι → Prop}
    (hB : ∀ i j, B j → V i j ≠ 0 → B i) : (chains F B).map V.mulVecLin = chains F B := by
  have hinj : Function.Injective V.mulVecLin :=
    Matrix.mulVec_injective_iff_isUnit.2 ((Matrix.isUnit_iff_isUnit_det V).2 hV)
  refine Submodule.eq_of_le_of_finrank_eq ?_ ?_
  · rintro _ ⟨c, hc, rfl⟩ i hi
    change ∑ s, V i s * c s = 0
    refine Finset.sum_eq_zero fun s _ => ?_
    by_cases hs : B s
    · have hzero : V i s = 0 := by
        by_contra hne
        exact hi (hB i s hs hne)
      simp [hzero]
    · simp [hc s hs]
  · exact (Submodule.equivMapOfInjective _ hinj _).finrank_eq.symm

omit [Fintype ι] in
/-- Subtracting columns `J` of `A` from column `a` can clear every collision with them: take the
combination whose last nonzero entry is earliest. -/
theorem exists_column_coeffs [Finite ι] (A : Matrix ι ι F) (a : ι) (J : Finset ι) :
    ∃ c : ι → F, ∀ i l, l ∈ J → IsLow A i l →
      ¬ ((A i a - ∑ s ∈ J, c s * A i s ≠ 0) ∧ ∀ r, i < r → A r a - ∑ s ∈ J, c s * A r s = 0) := by
  classical
  cases nonempty_fintype ι
  let y : (ι → F) → ι → F := fun c r => A r a - ∑ s ∈ J, c s * A r s
  let μ : (ι → F) → ℕ := fun c => (Finset.univ.filter fun q => ∃ r, q ≤ r ∧ y c r ≠ 0).card
  refine ⟨Function.argmin μ, fun i l hl hlow ⟨hne, htop⟩ => ?_⟩
  set c := Function.argmin μ with hc
  let c' : ι → F := fun s => c s + if s = l then y c i / A i l else 0
  have hy' : ∀ r, y c' r = y c r - y c i / A i l * A r l := by
    intro r
    simp only [y, c', add_mul, Finset.sum_add_distrib, ite_mul, zero_mul, Finset.sum_ite_eq', hl,
      if_true]
    ring
  have hzero : ∀ r, i ≤ r → y c' r = 0 := by
    intro r hr
    rcases eq_or_lt_of_le hr with rfl | hlt
    · rw [hy', div_mul_cancel₀ _ hlow.1, sub_self]
    · rw [hy', show y c r = 0 from htop r hlt, hlow.2 r hlt, mul_zero, sub_zero]
  refine Function.not_lt_argmin μ c' (Finset.card_lt_card ?_)
  refine (Finset.ssubset_iff_of_subset ?_).2 ⟨i, ?_, ?_⟩
  · intro q hq
    obtain ⟨r, hqr, hr⟩ := (Finset.mem_filter.mp hq).2
    have hri : r < i := by
      by_contra hle
      exact hr (hzero r (not_lt.mp hle))
    exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, i, (hqr.trans_lt hri).le, hne⟩
  · exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, i, le_rfl, hne⟩
  · intro hq
    obtain ⟨r, hqr, hr⟩ := (Finset.mem_filter.mp hq).2
    exact hr (hzero r hqr)

namespace FilteredComplex

variable (C : FilteredComplex ι F)

/-- The `p`-cells of `S`. -/
def cells (S : ι → Prop) (p : ℕ) : ι → Prop := fun i => S i ∧ C.deg i = p

/-- `Z_p(S, K)`: the `p`-chains on `S` whose boundary lies in `K`. -/
def relCycles (K S : ι → Prop) (p : ℕ) : Submodule F (ι → F) :=
  chains F (C.cells S p) ⊓ (chains F K).comap C.boundary.mulVecLin

/-- `B_p(S, K) = ∂C_{p+1}(S) + C_p(K)`. -/
def relBoundaries (K S : ι → Prop) (p : ℕ) : Submodule F (ι → F) :=
  (chains F (C.cells S (p + 1))).map C.boundary.mulVecLin ⊔ chains F (C.cells K p)

/-- `H_p(S, K) = Z_p(S, K) / B_p(S, K)`. -/
abbrev RelHomology (K S : ι → Prop) (p : ℕ) : Type _ :=
  C.relCycles K S p ⧸ (C.relBoundaries K S p).comap (C.relCycles K S p).subtype

variable {C} {K S T : ι → Prop}

theorem relCycles_mono (hST : ∀ i, S i → T i) (p : ℕ) :
    C.relCycles K S p ≤ C.relCycles K T p :=
  inf_le_inf_right _ (chains_mono fun i h => ⟨hST i h.1, h.2⟩)

theorem relBoundaries_mono (hST : ∀ i, S i → T i) (p : ℕ) :
    C.relBoundaries K S p ≤ C.relBoundaries K T p :=
  sup_le_sup_right (Submodule.map_mono (chains_mono fun i h => ⟨hST i h.1, h.2⟩)) _

variable (C) in
/-- The map `H_p(S, K) → H_p(T, K)` that the inclusion induces. -/
def persistence (hST : ∀ i, S i → T i) (p : ℕ) : C.RelHomology K S p →ₗ[F] C.RelHomology K T p :=
  Submodule.mapQ _ _ (Submodule.inclusion (relCycles_mono hST p))
    fun _ hx => relBoundaries_mono hST p hx

/-- The boundary of a chain on an initial segment stays on it. -/
theorem boundary_mem_chains (hK : IsInitial K) {w : ι → F} (hw : w ∈ chains F K) :
    C.boundary.mulVec w ∈ chains F K := by
  intro r hr
  refine Finset.sum_eq_zero fun s _ => ?_
  by_cases hs : K s
  · have hzero : C.boundary r s = 0 := by
      by_contra hne
      exact hr (hK r s (C.boundary_lt r s hne).le hs)
    simp [hzero]
  · simp [hw s hs]

/-- The quotient is honest: on initial segments, every relative boundary is a relative cycle. -/
theorem relBoundaries_le_relCycles (hK : IsInitial K) (hS : IsInitial S) (hKS : ∀ i, K i → S i)
    (p : ℕ) : C.relBoundaries K S p ≤ C.relCycles K S p := by
  refine sup_le ?_ fun w hw => ⟨chains_mono (fun i h => ⟨hKS i h.1, h.2⟩) hw,
    boundary_mem_chains hK (chains_mono (fun i h => h.1) hw)⟩
  rintro _ ⟨z, hz, rfl⟩
  refine ⟨fun r hr => ?_, ?_⟩
  · change ∑ s, C.boundary r s * z s = 0
    refine Finset.sum_eq_zero fun s _ => ?_
    by_cases hs : C.cells S (p + 1) s
    · have hzero : C.boundary r s = 0 := by
        by_contra hne
        refine hr ⟨hS r s (C.boundary_lt r s hne).le hs.1, ?_⟩
        have := C.boundary_deg r s hne
        have := hs.2
        omega
      simp [hzero]
    · simp [hz s hs]
  · change C.boundary.mulVec (C.boundary.mulVec z) ∈ chains F K
    rw [Matrix.mulVec_mulVec, C.boundary_sq, Matrix.zero_mulVec]
    exact Submodule.zero_mem _

/-- The rank of `H_p(S, K) → H_p(T, K)` is `dim Z_p(S, K) - dim (Z_p(S, K) ∩ B_p(T, K))`. -/
theorem finrank_range_persistence (hST : ∀ i, S i → T i) (p : ℕ) :
    finrank F (LinearMap.range (C.persistence (K := K) hST p)) +
      finrank F ↥(C.relCycles K S p ⊓ C.relBoundaries K T p) = finrank F (C.relCycles K S p) := by
  set Z := C.relCycles K S p with hZ
  set B := C.relBoundaries K T p with hB
  let f : Z →ₗ[F] C.RelHomology K T p :=
    Submodule.mkQ _ ∘ₗ Submodule.inclusion (relCycles_mono hST p)
  have hrange : LinearMap.range f = LinearMap.range (C.persistence (K := K) hST p) := by
    have hcomp : f = C.persistence (K := K) hST p ∘ₗ Submodule.mkQ _ := rfl
    rw [hcomp, LinearMap.range_comp, Submodule.range_mkQ, Submodule.map_top]
  have hker : LinearMap.ker f = B.comap Z.subtype := by
    ext x
    simp [f, hB, Submodule.Quotient.mk_eq_zero]
  have hinf : (Z ⊓ B).comap Z.subtype = B.comap Z.subtype := by
    rw [Submodule.comap_inf, Submodule.comap_subtype_self, top_inf_eq]
  have hfin : finrank F (B.comap Z.subtype) = finrank F ↥(Z ⊓ B) := by
    rw [← hinf]
    exact (Submodule.comapSubtypeEquivOfLe (inf_le_left : Z ⊓ B ≤ Z)).finrank_eq
  rw [← hfin, ← hker, ← hrange]
  exact LinearMap.finrank_range_add_finrank_ker f

/-- A reduction `D * V = R` whose column operations stay within a degree. -/
structure GradedReduction (C : FilteredComplex ι F) (R V : Matrix ι ι F) : Prop where
  mul : C.boundary * V = R
  reduced : Reduced R
  upper : ∀ i j, j < i → V i j = 0
  diag : ∀ i, V i i ≠ 0
  homog : ∀ i j, V i j ≠ 0 → C.deg i = C.deg j

theorem deg_of_mul_ne_zero {V : Matrix ι ι F} (hh : ∀ i j, V i j ≠ 0 → C.deg i = C.deg j) {i j : ι}
    (h : (C.boundary * V) i j ≠ 0) : C.deg j = C.deg i + 1 := by
  rw [Matrix.mul_apply] at h
  obtain ⟨s, _, hs⟩ := Finset.exists_ne_zero_of_sum_ne_zero h
  rw [← hh s j (right_ne_zero_of_mul hs)]
  exact C.boundary_deg i s (left_ne_zero_of_mul hs)

section Graded

variable {R V : Matrix ι ι F}

theorem GradedReduction.isUnit (g : C.GradedReduction R V) : IsUnit V.det :=
  isUnit_det g.upper g.diag

theorem GradedReduction.deg (g : C.GradedReduction R V) {i j : ι} (h : R i j ≠ 0) :
    C.deg j = C.deg i + 1 :=
  deg_of_mul_ne_zero g.homog (by rwa [g.mul])

theorem GradedReduction.boundary_mul (g : C.GradedReduction R V) : C.boundary * R = 0 := by
  rw [← g.mul, ← Matrix.mul_assoc, C.boundary_sq, Matrix.zero_mul]

theorem GradedReduction.inv_upper (g : C.GradedReduction R V) (i j : ι) (hij : j < i) :
    V⁻¹ i j = 0 := by
  have hb : V.BlockTriangular id := fun a b hab => g.upper a b hab
  letI := Matrix.invertibleOfIsUnitDet V g.isUnit
  exact Matrix.blockTriangular_inv_of_blockTriangular hb hij

/-- A graded reduction is a reduction in the sense of `RelativePersistence`. -/
theorem GradedReduction.toIsReduction (g : C.GradedReduction R V) :
    IsReduction C.boundary R V⁻¹ where
  factor := by rw [← g.mul, Matrix.mul_assoc, Matrix.mul_nonsing_inv V g.isUnit, Matrix.mul_one]
  reduced := g.reduced
  upper := g.inv_upper
  diag := by
    intro i hi
    have h1 : (V⁻¹ * V) i i = 1 := by rw [Matrix.nonsing_inv_mul V g.isUnit, Matrix.one_apply_eq]
    rw [Matrix.mul_apply, Finset.sum_eq_single i] at h1
    · simp [hi] at h1
    · intro s _ hs
      rcases lt_or_gt_of_ne hs with hlt | hgt
      · simp [g.inv_upper i s hlt]
      · simp [g.upper s i hgt]
    · simp

theorem cells_closed (g : C.GradedReduction R V) (hS : IsInitial S) (p : ℕ) :
    ∀ i j, C.cells S p j → V i j ≠ 0 → C.cells S p i := by
  intro i j hj hne
  have hij : i ≤ j := by
    by_contra hlt
    exact hne (g.upper i j (not_le.mp hlt))
  exact ⟨hS i j hij hj.1, (g.homog i j hne).trans hj.2⟩

theorem V_injective (g : C.GradedReduction R V) : Function.Injective V.mulVecLin :=
  Matrix.mulVec_injective_iff_isUnit.2 ((Matrix.isUnit_iff_isUnit_det V).2 g.isUnit)

/-- The relative cycles are spanned by the columns of `V` born on `S` whose boundary lands in `K`:
those with an empty column, or with their low in `K`. -/
theorem relCycles_eq (g : C.GradedReduction R V) (hK : IsInitial K) (hS : IsInitial S) (p : ℕ) :
    C.relCycles K S p =
      (chains F fun j => C.cells S p j ∧ ¬ LowOutside R K j).map V.mulVecLin := by
  have hA := chains_map_eq (F := F) g.isUnit (cells_closed g hS p)
  ext x
  constructor
  · rintro ⟨hx, hb⟩
    rw [SetLike.mem_coe, ← hA] at hx
    obtain ⟨c, hc, rfl⟩ := hx
    have hRc : R.mulVec c ∈ chains F K := by
      have hbd : C.boundary.mulVec (V.mulVec c) ∈ chains F K := hb
      rwa [Matrix.mulVec_mulVec, g.mul] at hbd
    have hlow := (mulVec_mem_chains_iff g.reduced hK).1 hRc
    refine ⟨c, fun i hi => ?_, rfl⟩
    by_contra hne
    exact hi ⟨by_contra fun h => hne (hc i h), hlow i hne⟩
  · rintro ⟨c, hc, rfl⟩
    have hcA : c ∈ chains F (C.cells S p) := chains_mono (fun i h => h.1) hc
    refine ⟨?_, ?_⟩
    · change V.mulVecLin c ∈ chains F (C.cells S p)
      rw [← hA]
      exact ⟨c, hcA, rfl⟩
    · change C.boundary.mulVec (V.mulVec c) ∈ chains F K
      rw [Matrix.mulVec_mulVec, g.mul]
      refine (mulVec_mem_chains_iff g.reduced hK).2 fun j hj => ?_
      intro hl
      exact hj (hc j fun h => h.2 hl)

theorem finrank_relCycles (g : C.GradedReduction R V) (hK : IsInitial K) (hS : IsInitial S)
    (p : ℕ) : finrank F (C.relCycles K S p) =
      cellCount fun j => C.cells S p j ∧ ¬ LowOutside R K j := by
  rw [relCycles_eq g hK hS p, ← finrank_chains (F := F)]
  exact (Submodule.equivMapOfInjective _ (V_injective g) _).finrank_eq.symm

/-- The columns that kill, within `T`, a class of `H_p(S, K)`: of degree `p + 1`, with their low
outside `K` and on `S`. -/
def Killer (C : FilteredComplex ι F) (R : Matrix ι ι F) (K S T : ι → Prop) (p : ℕ) (j : ι) : Prop :=
  C.cells T (p + 1) j ∧ LowOutside R K j ∧ ¬ LowOutside R S j

/-- The classes of `H_p(S, K)` that are dead in `T` are spanned by the killers' columns, and `K`'s
own chains are the zero class. -/
theorem relCycles_inf_relBoundaries (g : C.GradedReduction R V) (hK : IsInitial K)
    (hS : IsInitial S) (hT : IsInitial T) (hKS : ∀ i, K i → S i) (p : ℕ) :
    C.relCycles K S p ⊓ C.relBoundaries K T p =
      (chains F (Killer C R K S T p)).map R.mulVecLin ⊔ chains F (C.cells K p) := by
  classical
  have hA' := chains_map_eq (F := F) g.isUnit (cells_closed g hT (p + 1))
  have hKA : chains F (C.cells K p) ≤ chains F (C.cells S p) :=
    chains_mono fun i h => ⟨hKS i h.1, h.2⟩
  apply le_antisymm
  · rintro x ⟨⟨hxA, -⟩, hxB⟩
    obtain ⟨y, hy, w, hw, rfl⟩ := Submodule.mem_sup.1 hxB
    obtain ⟨z, hz, rfl⟩ := hy
    rw [← hA'] at hz
    obtain ⟨c, hc, rfl⟩ := hz
    have hbd : C.boundary.mulVecLin (V.mulVecLin c) = R.mulVec c := by
      change C.boundary.mulVec (V.mulVec c) = R.mulVec c
      rw [Matrix.mulVec_mulVec, g.mul]
    rw [hbd] at hxA ⊢
    have hRc : R.mulVec c ∈ chains F (C.cells S p) := by
      have hsub := Submodule.sub_mem _ hxA (hKA hw)
      rwa [add_sub_cancel_right] at hsub
    have hcS := (mulVec_mem_chains_iff g.reduced hS).1 (chains_mono (fun i h => h.1) hRc)
    let cG : ι → F := fun j => if LowOutside R K j then c j else 0
    let cK : ι → F := c - cG
    have hsplit : R.mulVec c = R.mulVec cG + R.mulVec cK := by
      rw [← Matrix.mulVec_add]
      simp [cK]
    have hcG : cG ∈ chains F (Killer C R K S T p) := by
      intro i hi
      by_cases hl : LowOutside R K i
      · have hci : c i = 0 := by
          by_contra hne
          exact hi ⟨by_contra fun h => hne (hc i h), hl, hcS i hne⟩
        simp [cG, hci]
      · simp [cG, hl]
    have hRcK : R.mulVec cK ∈ chains F (C.cells K p) := by
      intro r hr
      refine Finset.sum_eq_zero fun s _ => ?_
      by_cases hl : LowOutside R K s
      · simp [cK, cG, hl]
      · by_cases hcs : c s = 0
        · simp [cK, cG, hl, hcs]
        · have hRrs : R r s = 0 := by
            by_contra hne
            refine hr ⟨by_contra fun hKr => hne (column_eq_zero hK hl r hKr), ?_⟩
            have hd := g.deg hne
            have hs : C.cells T (p + 1) s := by_contra fun h => hcs (hc s h)
            have := hs.2
            omega
          simp [hRrs]
    rw [hsplit, add_assoc]
    exact Submodule.mem_sup.2 ⟨R.mulVec cG, ⟨cG, hcG, rfl⟩, R.mulVec cK + w,
      Submodule.add_mem _ hRcK hw, rfl⟩
  · refine sup_le ?_ fun w hw =>
      ⟨⟨hKA hw, boundary_mem_chains hK (chains_mono (fun i h => h.1) hw)⟩,
        Submodule.mem_sup_right hw⟩
    rintro _ ⟨c, hc, rfl⟩
    refine ⟨⟨fun r hr => ?_, ?_⟩, ?_⟩
    · change ∑ s, R r s * c s = 0
      refine Finset.sum_eq_zero fun s _ => ?_
      by_cases hs : Killer C R K S T p s
      · have hzero : R r s = 0 := by
          by_contra hne
          refine hr ⟨by_contra fun hSr => hne (column_eq_zero hS hs.2.2 r hSr), ?_⟩
          have hd := g.deg hne
          have := hs.1.2
          omega
        simp [hzero]
      · simp [hc s hs]
    · change C.boundary.mulVec (R.mulVec c) ∈ chains F K
      rw [Matrix.mulVec_mulVec, g.boundary_mul, Matrix.zero_mulVec]
      exact Submodule.zero_mem _
    · refine Submodule.mem_sup_left ⟨V.mulVec c, ?_, ?_⟩
      · change V.mulVecLin c ∈ chains F (C.cells T (p + 1))
        rw [← hA']
        exact ⟨c, chains_mono (fun i h => h.1) hc, rfl⟩
      · change C.boundary.mulVec (V.mulVec c) = R.mulVec c
        rw [Matrix.mulVec_mulVec, g.mul]

theorem finrank_relCycles_inf (g : C.GradedReduction R V) (hK : IsInitial K) (hS : IsInitial S)
    (hT : IsInitial T) (hKS : ∀ i, K i → S i) (p : ℕ) :
    finrank F ↥(C.relCycles K S p ⊓ C.relBoundaries K T p) =
      cellCount (Killer C R K S T p) + cellCount (C.cells K p) := by
  rw [relCycles_inf_relBoundaries g hK hS hT hKS p]
  set M := (chains F (Killer C R K S T p)).map R.mulVecLin with hMdef
  have hKil : ∀ c ∈ chains F (Killer C R K S T p), ∀ j, c j ≠ 0 → LowOutside R K j := by
    intro c hc j hj
    exact (by_contra fun h => hj (hc j h) : Killer C R K S T p j).2.1
  have hinj : Function.Injective (R.mulVecLin.domRestrict (chains F (Killer C R K S T p))) := by
    refine (injective_iff_map_eq_zero _).2 fun x hx => Subtype.ext ?_
    refine eq_zero_of_mulVec_mem (X := K) g.reduced (hKil x.1 x.2) ?_
    have hzero : R.mulVec x.1 = 0 := hx
    rw [hzero]
    exact Submodule.zero_mem _
  have hM : finrank F M = cellCount (Killer C R K S T p) := by
    rw [hMdef, ← LinearMap.range_domRestrict, LinearMap.finrank_range_of_inj hinj,
      finrank_chains]
  have hdisj : M ⊓ chains F (C.cells K p) = ⊥ := by
    refine (Submodule.eq_bot_iff _).2 ?_
    rintro _ ⟨⟨c, hc, rfl⟩, hy⟩
    have hmem : R.mulVec c ∈ chains F K := chains_mono (fun i h => h.1) hy
    have hzero := eq_zero_of_mulVec_mem g.reduced (hKil c hc) hmem
    simp [hzero]
  have hsum := Submodule.finrank_sup_add_finrank_inf_eq M (chains F (C.cells K p))
  rw [hdisj, finrank_bot, add_zero, hM, finrank_chains] at hsum
  exact hsum

/-- Cell `i` carries a class of `H_p(S, K)` still alive in `T`. It is on `S` and outside `K`, its
column is empty or its low is in `K`, and nothing in `T` kills it. -/
def Alive (C : FilteredComplex ι F) (R : Matrix ι ι F) (K S T : ι → Prop) (p : ℕ) (i : ι) : Prop :=
  C.deg i = p ∧ S i ∧ ¬ K i ∧ (∀ l, IsLow R l i → K l) ∧ ∀ j, T j → ¬ IsLow R i j

/-- An initial segment of a filtered complex is a first subcomplex. -/
theorem isFirstSubcomplex (hK : IsInitial K) : IsFirstSubcomplex C.boundary K :=
  ⟨hK, fun i l hl hi => by_contra fun hne => hi (hK i l (C.boundary_lt i l hne).le hl)⟩

/-- The rank of `H_p(S, K) → H_p(T, K)`, read from a graded reduction. -/
theorem finrank_persistence_of_graded (g : C.GradedReduction R V) (hK : IsInitial K)
    (hS : IsInitial S) (hT : IsInitial T) (hKS : ∀ i, K i → S i) (hST : ∀ i, S i → T i) (p : ℕ) :
    finrank F (LinearMap.range (C.persistence (K := K) hST p)) =
      cellCount (Alive C R K S T p) := by
  have hrank := finrank_range_persistence (C := C) (K := K) hST p
  rw [finrank_relCycles g hK hS p, finrank_relCycles_inf g hK hS hT hKS p] at hrank
  have hred := g.toIsReduction
  have hfirst := isFirstSubcomplex (C := C) hK
  have hsplit₁ := cellCount_split (fun j => C.cells S p j ∧ ¬ LowOutside R K j) K
  have hsplit₂ := cellCount_split (fun j => (C.cells S p j ∧ ¬ LowOutside R K j) ∧ ¬ K j)
    (fun i => ∃ j, T j ∧ IsLow R i j)
  have hin : cellCount (fun j => (C.cells S p j ∧ ¬ LowOutside R K j) ∧ K j) =
      cellCount (C.cells K p) := by
    refine cellCount_congr fun i =>
      ⟨fun h => ⟨h.2, h.1.1.2⟩, fun h => ⟨⟨⟨hKS i h.1, h.2⟩, ?_⟩, h.1⟩⟩
    rintro ⟨q, hq, hqK⟩
    exact hq.1 (reduction_closed hfirst hred q hqK i h.1)
  have halive : cellCount (fun i => ((C.cells S p i ∧ ¬ LowOutside R K i) ∧ ¬ K i) ∧
      ¬ ∃ j, T j ∧ IsLow R i j) = cellCount (Alive C R K S T p) := by
    refine cellCount_congr fun i => ?_
    rw [not_lowOutside_iff]
    constructor
    · rintro ⟨⟨⟨⟨hSi, hdi⟩, hl⟩, hKi⟩, hk⟩
      exact ⟨hdi, hSi, hKi, hl, fun j hj hlow => hk ⟨j, hj, hlow⟩⟩
    · rintro ⟨hdi, hSi, hKi, hl, hk⟩
      exact ⟨⟨⟨⟨hSi, hdi⟩, hl⟩, hKi⟩, fun ⟨j, hj, hlow⟩ => hk j hj hlow⟩
  have hkill : cellCount (Killer C R K S T p) = cellCount (fun i =>
      ((C.cells S p i ∧ ¬ LowOutside R K i) ∧ ¬ K i) ∧ ∃ j, T j ∧ IsLow R i j) := by
    refine cellCount_bij (lowOf R) ?_ ?_ ?_
    · rintro j ⟨⟨hTj, hdj⟩, ⟨i, hlow, hiK⟩, hnS⟩
      rw [lowOf_eq hlow]
      have hd := g.deg hlow.1
      refine ⟨⟨⟨⟨by_contra fun hSi => hnS ⟨i, hlow, hSi⟩, by omega⟩, ?_⟩, hiK⟩, j, hTj, hlow⟩
      rintro ⟨q, hq, _⟩
      exact hq.1 (column_zero_of_isLow C.boundary_sq hred hlow q)
    · rintro j₁ j₂ ⟨_, ⟨i₁, hl₁, _⟩, _⟩ ⟨_, ⟨i₂, hl₂, _⟩, _⟩ hlow
      rw [lowOf_eq hl₁, lowOf_eq hl₂] at hlow
      subst hlow
      exact g.reduced i₁ j₁ j₂ hl₁ hl₂
    · rintro i ⟨⟨⟨⟨hSi, hdi⟩, _⟩, hiK⟩, j, hTj, hlow⟩
      have hd := g.deg hlow.1
      refine ⟨j, ⟨⟨hTj, by omega⟩, ⟨i, hlow, hiK⟩, ?_⟩, lowOf_eq hlow⟩
      rintro ⟨q, hq, hqS⟩
      exact hqS ((hq.unique hlow) ▸ hSi)
  omega

end Graded

/-- Every filtered complex has a reduction, by the usual column algorithm: each column in turn
subtracts earlier columns of its own degree until its low is unused. -/
theorem exists_gradedReduction (C : FilteredComplex ι F) : ∃ R V, C.GradedReduction R V := by
  classical
  suffices h : ∀ s : Finset ι, ∃ V : Matrix ι ι F, (∀ i j, j < i → V i j = 0) ∧ (∀ i, V i i ≠ 0) ∧
      (∀ i j, V i j ≠ 0 → C.deg i = C.deg j) ∧ ∀ i j₁ j₂, j₁ ∈ s → j₂ ∈ s →
        IsLow (C.boundary * V) i j₁ → IsLow (C.boundary * V) i j₂ → j₁ = j₂ by
    obtain ⟨V, hu, hd, hh, hr⟩ := h Finset.univ
    exact ⟨_, V, rfl, fun i j₁ j₂ => hr i j₁ j₂ (Finset.mem_univ _) (Finset.mem_univ _), hu, hd, hh⟩
  intro s
  induction s using Finset.induction_on_max with
  | empty =>
    refine ⟨1, fun i j hij => Matrix.one_apply_ne hij.ne', fun i => by simp, fun i j h => ?_,
      by simp⟩
    have hij : i = j := by
      by_contra hne
      exact h (Matrix.one_apply_ne hne)
    rw [hij]
  | insert a s hlt ih =>
    obtain ⟨V, hu, hd, hh, hr⟩ := ih
    obtain ⟨c, hc⟩ := exists_column_coeffs (C.boundary * V) a (s.filter fun l => C.deg l = C.deg a)
    set J := s.filter fun l => C.deg l = C.deg a with hJ
    have hJlt : ∀ l ∈ J, l ∈ s ∧ C.deg l = C.deg a := fun l hl => Finset.mem_filter.1 hl
    let V' : Matrix ι ι F := fun i j => if j = a then V i a - ∑ l ∈ J, c l * V i l else V i j
    have hcol : ∀ i, V' i a = V i a - ∑ l ∈ J, c l * V i l := fun i => by simp [V']
    have hne : ∀ i j, j ≠ a → V' i j = V i j := fun i j hj => by simp [V', hj]
    have hu' : ∀ i j, j < i → V' i j = 0 := by
      intro i j hij
      by_cases hj : j = a
      · rw [hj] at hij ⊢
        rw [hcol, hu i a hij, zero_sub, neg_eq_zero]
        refine Finset.sum_eq_zero fun l hl => ?_
        rw [hu i l ((hlt l (hJlt l hl).1).trans hij), mul_zero]
      · rw [hne i j hj]
        exact hu i j hij
    have hd' : ∀ i, V' i i ≠ 0 := by
      intro i
      by_cases hi : i = a
      · rw [hi, hcol, Finset.sum_eq_zero fun l hl => by rw [hu a l (hlt l (hJlt l hl).1), mul_zero],
          sub_zero]
        exact hd a
      · rw [hne i i hi]
        exact hd i
    have hh' : ∀ i j, V' i j ≠ 0 → C.deg i = C.deg j := by
      intro i j hij
      by_cases hj : j = a
      · rw [hj] at hij ⊢
        by_contra hdeg
        apply hij
        rw [hcol, show V i a = 0 from by_contra fun h => hdeg (hh i a h), zero_sub, neg_eq_zero]
        refine Finset.sum_eq_zero fun l hl => ?_
        have hzero : V i l = 0 := by
          by_contra h
          exact hdeg ((hh i l h).trans (hJlt l hl).2)
        rw [hzero, mul_zero]
      · rw [hne i j hj] at hij
        exact hh i j hij
    have hmul : ∀ r, (C.boundary * V') r a =
        (C.boundary * V) r a - ∑ s ∈ J, c s * (C.boundary * V) r s := by
      intro r
      simp only [Matrix.mul_apply, hcol, mul_sub, Finset.sum_sub_distrib, Finset.mul_sum]
      congr 1
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun l _ => Finset.sum_congr rfl fun x _ => ?_
      ring
    have hmul' : ∀ j, j ≠ a → ∀ r, (C.boundary * V') r j = (C.boundary * V) r j := by
      intro j hj r
      simp only [Matrix.mul_apply, hne _ j hj]
    have hlow' : ∀ i j, j ≠ a → (IsLow (C.boundary * V') i j ↔ IsLow (C.boundary * V) i j) := by
      intro i j hj
      simp only [IsLow, hmul' j hj]
    have hcoll : ∀ i l, l ∈ s → IsLow (C.boundary * V') i a → IsLow (C.boundary * V') i l →
        False := by
      intro i l hl h₁ h₂
      have hl' := (hlow' i l (hlt l hl).ne).1 h₂
      have hdl := deg_of_mul_ne_zero hh hl'.1
      have hda := deg_of_mul_ne_zero hh' h₁.1
      refine hc i l (Finset.mem_filter.2 ⟨hl, by omega⟩) hl' ⟨?_, fun r hr => ?_⟩
      · rw [← hmul]
        exact h₁.1
      · rw [← hmul]
        exact h₁.2 r hr
    refine ⟨V', hu', hd', hh', fun i j₁ j₂ hj₁ hj₂ h₁ h₂ => ?_⟩
    rcases Finset.mem_insert.1 hj₁ with e₁ | hj₁ <;> rcases Finset.mem_insert.1 hj₂ with e₂ | hj₂
    · exact e₁.trans e₂.symm
    · exact (hcoll i j₂ hj₂ (e₁ ▸ h₁) h₂).elim
    · exact (hcoll i j₁ hj₁ (e₂ ▸ h₂) h₁).elim
    · exact hr i j₁ j₂ hj₁ hj₂ ((hlow' i j₁ (hlt _ hj₁).ne).1 h₁) ((hlow' i j₂ (hlt _ hj₂).ne).1 h₂)

variable {R U : Matrix ι ι F}

/-- **The relative barcode is the persistent homology of the pair.** For any reduction of the
boundary, the rank of `H_p(S, K) → H_p(T, K)` is the number of cells alive from `S` through `T`,
which `alive_iff_table` reads as the table's bars. -/
theorem finrank_persistence (hK : IsInitial K) (hS : IsInitial S) (hT : IsInitial T)
    (hKS : ∀ i, K i → S i) (hST : ∀ i, S i → T i) (p : ℕ) (h : IsReduction C.boundary R U) :
    finrank F (LinearMap.range (C.persistence (K := K) hST p)) =
      cellCount (Alive C R K S T p) := by
  obtain ⟨Rg, Vg, g⟩ := C.exists_gradedReduction
  rw [finrank_persistence_of_graded g hK hS hT hKS hST p]
  refine cellCount_congr fun i => ?_
  simp only [Alive, isLow_iff_of_reductions h g.toIsReduction]

/-- The table: a class alive from `S` through `T` is an absolute bar born at `i` outside `K` and not
yet killed in `T`, or a class born in `K` and killed at `i`, which starts a relative class at `i`
one degree up and is never killed. Bars born in `K` and not killed outside it leave no trace. -/
theorem alive_iff_table (hK : IsInitial K) (h : IsReduction C.boundary R U) (p : ℕ) (i : ι) :
    Alive C R K S T p i ↔ C.deg i = p ∧ S i ∧ ¬ K i ∧
      (((∀ r, R r i = 0) ∧ ∀ j, T j → ¬ IsLow R i j) ∨ ∃ l, K l ∧ IsLow R l i) := by
  have hfirst := isFirstSubcomplex (C := C) hK
  constructor
  · rintro ⟨hd, hSi, hKi, hl, hk⟩
    refine ⟨hd, hSi, hKi, ?_⟩
    by_cases hc : ∃ r, R r i ≠ 0
    · obtain ⟨l, hlow⟩ := exists_isLow hc
      exact Or.inr ⟨l, hl l hlow, hlow⟩
    · push Not at hc
      exact Or.inl ⟨hc, hk⟩
  · rintro ⟨hd, hSi, hKi, (⟨hz, hk⟩ | ⟨l, hlK, hlow⟩)⟩
    · exact ⟨hd, hSi, hKi, fun l hl => absurd (hz l) hl.1, hk⟩
    · refine ⟨hd, hSi, hKi, fun q hq => (hq.unique hlow) ▸ hlK, fun j _ hj => ?_⟩
      exact hlow.1 (column_zero_of_isLow C.boundary_sq h hj l)

/-- With `K` empty the cycles are the ordinary ones. -/
theorem relCycles_empty (S : ι → Prop) (p : ℕ) :
    C.relCycles (fun _ => False) S p =
      chains F (C.cells S p) ⊓ LinearMap.ker C.boundary.mulVecLin := by
  have hbot : chains F (fun _ : ι => False) = ⊥ := by
    refine (Submodule.eq_bot_iff _).2 fun x hx => funext fun i => hx i not_false
  rw [relCycles, hbot, LinearMap.ker]

/-- With `K` empty the boundaries are the ordinary ones. -/
theorem relBoundaries_empty (S : ι → Prop) (p : ℕ) :
    C.relBoundaries (fun _ => False) S p =
      (chains F (C.cells S (p + 1))).map C.boundary.mulVecLin := by
  have hbot : chains F (C.cells (fun _ : ι => False) p) = ⊥ := by
    refine (Submodule.eq_bot_iff _).2 fun x hx => funext fun i => hx i fun h => h.1
  rw [relBoundaries, hbot, sup_bot_eq]

/-- **Reduction computes persistent homology.** With `K` empty, the rank of `H_p(S) → H_p(T)` counts
the cells of degree `p` on `S` with an empty column that nothing in `T` kills. -/
theorem finrank_persistence_absolute (hS : IsInitial S) (hT : IsInitial T) (hST : ∀ i, S i → T i)
    (p : ℕ) (h : IsReduction C.boundary R U) :
    finrank F (LinearMap.range (C.persistence (K := fun _ => False) hST p)) =
      cellCount fun i => C.deg i = p ∧ S i ∧ (∀ l, ¬ IsLow R l i) ∧ ∀ j, T j → ¬ IsLow R i j := by
  rw [finrank_persistence (fun _ _ _ hb => hb) hS hT (fun _ hf => hf.elim) hST p h]
  exact cellCount_congr fun i =>
    ⟨fun ⟨hd, hSi, _, hl, hk⟩ => ⟨hd, hSi, hl, hk⟩,
      fun ⟨hd, hSi, hl, hk⟩ => ⟨hd, hSi, not_false, fun l hlow => hl l hlow, hk⟩⟩

section Quotient

variable [DecidablePred K]

variable (C K) in
/-- The complex with `K`'s cells deleted: the chain complex `C / C(K)` in the basis of the remaining
cells, which is what `TDA.Ph.RelativeFiltration` reduces. -/
def quotient (hK : IsInitial K) : FilteredComplex {x // ¬ K x} F where
  boundary := rel K C.boundary
  deg := fun i => C.deg i.1
  boundary_sq := by
    have hfirst := isFirstSubcomplex (C := C) hK
    ext i j
    have hsplit :=
      Fintype.sum_subtype_add_sum_subtype K (fun s => C.boundary i.1 s * C.boundary s j.1)
    have hin : ∑ s : {x // K x}, C.boundary i.1 s.1 * C.boundary s.1 j.1 = 0 :=
      Finset.sum_eq_zero fun s _ => by simp [hfirst.closed i.1 s.1 s.2 i.2]
    have hall : ∑ s, C.boundary i.1 s * C.boundary s j.1 = 0 := by
      rw [← Matrix.mul_apply, C.boundary_sq, Matrix.zero_apply]
    rw [hin, zero_add, hall] at hsplit
    simpa [rel, Matrix.mul_apply] using hsplit
  boundary_lt := fun i j h => C.boundary_lt i.1 j.1 h
  boundary_deg := fun i j h => C.boundary_deg i.1 j.1 h

/-- **The relative matrix is the right matrix.** Reducing the complex with `K`'s cells deleted gives
the ranks of the pair: its absolute persistence is the persistence of `H_p(S, K)`. No reduction
appears in the statement. -/
theorem finrank_quotient_persistence (hK : IsInitial K) (hS : IsInitial S) (hT : IsInitial T)
    (hKS : ∀ i, K i → S i) (hST : ∀ i, S i → T i) (p : ℕ) :
    finrank F (LinearMap.range ((C.quotient K hK).persistence (K := fun _ => False)
        (S := fun i => S i.1) (T := fun i => T i.1) (fun i => hST i.1) p)) =
      finrank F (LinearMap.range (C.persistence (K := K) hST p)) := by
  obtain ⟨Rg, Vg, g⟩ := C.exists_gradedReduction
  have hred := g.toIsReduction
  have hfirst := isFirstSubcomplex (C := C) hK
  have hred' : IsReduction (C.quotient K hK).boundary (rel K Rg) (rel K Vg⁻¹) :=
    relative_isReduction hfirst hred
  rw [finrank_persistence_absolute (C := C.quotient K hK) (fun a b hab hb => hS a.1 b.1 hab hb)
      (fun a b hab hb => hT a.1 b.1 hab hb) (fun i => hST i.1) p hred',
    finrank_persistence hK hS hT hKS hST p hred]
  refine cellCount_bij Subtype.val ?_ (fun i j _ _ hij => Subtype.ext hij) ?_
  · rintro i ⟨hd, hSi, hl, hk⟩
    refine ⟨hd, hSi, i.2, fun l hlow => by_contra fun hlK => ?_, fun j hj hlow => ?_⟩
    · exact hl ⟨l, hlK⟩ ((relative_isLow_iff hfirst ⟨l, hlK⟩ i).2 hlow)
    · have hjK := not_mem_of_isLow hfirst hred i.2 hlow
      exact hk ⟨j, hjK⟩ hj ((relative_isLow_iff hfirst i ⟨j, hjK⟩).2 hlow)
  · rintro b ⟨hd, hSb, hKb, hl, hk⟩
    refine ⟨⟨b, hKb⟩, ⟨hd, hSb, fun l hlow => ?_, fun j hj hlow => ?_⟩, rfl⟩
    · exact l.2 (hl l.1 ((relative_isLow_iff hfirst l ⟨b, hKb⟩).1 hlow))
    · exact hk j.1 hj ((relative_isLow_iff hfirst ⟨b, hKb⟩ j).1 hlow)

end Quotient

end FilteredComplex

/-- A filtered finite simplicial complex: its cells are distinct nonempty vertex sets, every proper
nonempty subset of a cell is a cell, and it comes first. -/
structure SimplicialFiltration (ι V : Type*) [LinearOrder ι] [DecidableEq V] where
  /-- The vertices of a cell. -/
  cell : ι → Finset V
  injective : Function.Injective cell
  nonempty : ∀ i, (cell i).Nonempty
  face : ∀ j (s : Finset V), s.Nonempty → s ⊂ cell j → ∃ i, cell i = s ∧ i < j

namespace SimplicialFiltration

section Boundary

variable {ι V : Type*} [LinearOrder ι] [DecidableEq V] (X : SimplicialFiltration ι V)

/-- The ℤ/2 boundary: each cell's codimension-one faces. -/
def boundary : Matrix ι ι (ZMod 2) := fun i j =>
  if X.cell i ⊂ X.cell j ∧ (X.cell i).card + 1 = (X.cell j).card then 1 else 0

variable {X} in
theorem boundary_ne_zero {i j : ι} : X.boundary i j ≠ 0 ↔
    X.cell i ⊂ X.cell j ∧ (X.cell i).card + 1 = (X.cell j).card := by
  unfold boundary
  split_ifs with h
  · simp [h]
  · simp [h]

/-- A codimension-two face has exactly two cells between it and the cell, so the boundary of a
boundary vanishes mod 2. -/
theorem boundary_sq [Fintype ι] : X.boundary * X.boundary = 0 := by
  classical
  ext i j
  rw [Matrix.mul_apply, Matrix.zero_apply]
  by_cases hex : ∃ l, X.boundary i l * X.boundary l j ≠ 0
  · obtain ⟨l₀, h₀⟩ := hex
    obtain ⟨hA₀, hc₁⟩ := boundary_ne_zero.1 (left_ne_zero_of_mul h₀)
    obtain ⟨hB₀, hc₂⟩ := boundary_ne_zero.1 (right_ne_zero_of_mul h₀)
    have hAB : X.cell i ⊂ X.cell j := hA₀.trans hB₀
    have hcard : (X.cell j \ X.cell i).card = 2 := by
      rw [Finset.card_sdiff_of_subset hAB.subset]
      omega
    obtain ⟨a, b, hab, hBA⟩ := Finset.card_eq_two.1 hcard
    have ha := Finset.mem_sdiff.1 (show a ∈ X.cell j \ X.cell i by simp [hBA])
    have hb := Finset.mem_sdiff.1 (show b ∈ X.cell j \ X.cell i by simp [hBA])
    have hins : ∀ x, x ∈ X.cell j → x ∉ X.cell i → insert x (X.cell i) ⊂ X.cell j := by
      intro x hxB hxA
      refine Finset.ssubset_iff_subset_ne.2 ⟨Finset.insert_subset hxB hAB.subset, fun h => ?_⟩
      have hc := Finset.card_insert_of_notMem hxA
      rw [h] at hc
      omega
    obtain ⟨la, hla, -⟩ := X.face j _ (Finset.insert_nonempty a _) (hins a ha.1 ha.2)
    obtain ⟨lb, hlb, -⟩ := X.face j _ (Finset.insert_nonempty b _) (hins b hb.1 hb.2)
    have hne : la ≠ lb := by
      intro h
      rw [h, hlb] at hla
      have hmem : a ∈ insert b (X.cell i) := by
        rw [hla]
        exact Finset.mem_insert_self _ _
      rcases Finset.mem_insert.1 hmem with h' | h'
      · exact hab h'
      · exact ha.2 h'
    have hone : ∀ x l, x ∈ X.cell j → x ∉ X.cell i → X.cell l = insert x (X.cell i) →
        X.boundary i l * X.boundary l j = 1 := by
      intro x l hxB hxA hl
      have hc := Finset.card_insert_of_notMem hxA
      have hfirst : X.cell i ⊂ X.cell l ∧ (X.cell i).card + 1 = (X.cell l).card := by
        rw [hl]
        exact ⟨Finset.ssubset_insert hxA, hc.symm⟩
      have hsecond : X.cell l ⊂ X.cell j ∧ (X.cell l).card + 1 = (X.cell j).card := by
        rw [hl]
        exact ⟨hins x hxB hxA, by omega⟩
      simp [boundary, hfirst, hsecond]
    have hterm : ∀ l, X.boundary i l * X.boundary l j = if l = la ∨ l = lb then 1 else 0 := by
      intro l
      by_cases hl : l = la ∨ l = lb
      · rw [if_pos hl]
        rcases hl with rfl | rfl
        · exact hone a _ ha.1 ha.2 hla
        · exact hone b _ hb.1 hb.2 hlb
      · rw [if_neg hl]
        by_contra hne0
        obtain ⟨hA, hc₃⟩ := boundary_ne_zero.1 (left_ne_zero_of_mul hne0)
        obtain ⟨hB, hc₄⟩ := boundary_ne_zero.1 (right_ne_zero_of_mul hne0)
        have hone' : (X.cell l \ X.cell i).card = 1 := by
          rw [Finset.card_sdiff_of_subset hA.subset]
          omega
        obtain ⟨x, hx⟩ := Finset.card_eq_one.1 hone'
        have hxmem := Finset.mem_sdiff.1 (show x ∈ X.cell l \ X.cell i by simp [hx])
        have hcell : X.cell l = insert x (X.cell i) := by
          ext y
          simp only [Finset.mem_insert]
          constructor
          · intro hy
            by_cases hyA : y ∈ X.cell i
            · exact Or.inr hyA
            · have hy' : y ∈ X.cell l \ X.cell i := Finset.mem_sdiff.2 ⟨hy, hyA⟩
              rw [hx] at hy'
              exact Or.inl (Finset.mem_singleton.1 hy')
          · rintro (rfl | hyA)
            · exact hxmem.1
            · exact hA.subset hyA
        have hxBA : x ∈ X.cell j \ X.cell i := Finset.mem_sdiff.2 ⟨hB.subset hxmem.1, hxmem.2⟩
        rw [hBA, Finset.mem_insert, Finset.mem_singleton] at hxBA
        rcases hxBA with rfl | rfl
        · exact hl (Or.inl (X.injective (hcell.trans hla.symm)))
        · exact hl (Or.inr (X.injective (hcell.trans hlb.symm)))
    rw [Finset.sum_congr rfl fun l _ => hterm l, Finset.sum_boole]
    have hf : (Finset.univ.filter fun l => l = la ∨ l = lb) = {la, lb} := by
      ext l
      simp
    rw [hf, Finset.card_pair hne]
    decide
  · push Not at hex
    exact Finset.sum_eq_zero fun l _ => hex l

/-- The filtered ℤ/2 chain complex of a simplicial filtration, a cell's degree its dimension. -/
def toComplex [Fintype ι] : FilteredComplex ι (ZMod 2) where
  boundary := X.boundary
  deg := fun i => (X.cell i).card - 1
  boundary_sq := X.boundary_sq
  boundary_lt := by
    intro i j h
    obtain ⟨hs, _⟩ := boundary_ne_zero.1 h
    obtain ⟨i', hi', hlt⟩ := X.face j _ (X.nonempty i) hs
    rwa [X.injective hi'] at hlt
  boundary_deg := by
    intro i j h
    obtain ⟨_, hc⟩ := boundary_ne_zero.1 h
    have hpos := (X.nonempty i).card_pos
    omega

end Boundary

section Example

/-! The definitions give the textbook answers on an edge `e` with endpoints `v₀ v₁`.

Against `K₀` the two endpoints, `H₁(e, ∂e)` has rank one and `H₀(e, ∂e)` is zero, where absolutely
`H₀(e)` has rank one and `H₁(e)` is zero. That relative class is the table's second row: `v₁` is
born in `K₀`, the edge kills it, and the relative class starts at the edge, one degree up.

Against `K₀ = {v₀}`, `v₁` carries a class of the pair that the edge kills, the first row, a finite
bar: the map from before the edge has rank zero into the edge and rank one into itself. The same
maps, absolutely, have ranks one and two. -/

private def edge : SimplicialFiltration (Fin 3) (Fin 2) where
  cell := ![{0}, {1}, {0, 1}]
  injective := by decide
  nonempty := by decide
  face := by decide

private theorem edge_reduction : IsReduction edge.toComplex.boundary edge.boundary 1 where
  factor := (Matrix.mul_one _).symm
  reduced := by
    simp only [Reduced, IsLow]
    decide
  upper := fun _ _ hij => Matrix.one_apply_ne hij.ne'
  diag := fun i => by simp

private theorem edge_relative_one :
    finrank (ZMod 2) (LinearMap.range (edge.toComplex.persistence (K := (· ≤ 1))
      (S := fun _ => True) (T := fun _ => True) (fun _ _ => trivial) 1)) = 1 := by
  rw [FilteredComplex.finrank_persistence (fun _ _ hab hb => hab.trans hb) (fun _ _ _ _ => trivial)
    (fun _ _ _ _ => trivial) (fun _ _ => trivial) _ 1 edge_reduction]
  have hiff : ∀ i, FilteredComplex.Alive edge.toComplex edge.boundary (· ≤ 1) (fun _ => True)
      (fun _ => True) 1 i ↔ i = 2 := by
    simp only [FilteredComplex.Alive, IsLow, true_and, forall_true_left]
    decide
  rw [cellCount_congr hiff]
  exact Nat.card_unique

private theorem edge_relative_zero :
    finrank (ZMod 2) (LinearMap.range (edge.toComplex.persistence (K := (· ≤ 1))
      (S := fun _ => True) (T := fun _ => True) (fun _ _ => trivial) 0)) = 0 := by
  rw [FilteredComplex.finrank_persistence (fun _ _ hab hb => hab.trans hb) (fun _ _ _ _ => trivial)
    (fun _ _ _ _ => trivial) (fun _ _ => trivial) _ 0 edge_reduction]
  have hiff : ∀ i, ¬ FilteredComplex.Alive edge.toComplex edge.boundary (· ≤ 1) (fun _ => True)
      (fun _ => True) 0 i := by
    simp only [FilteredComplex.Alive, IsLow, true_and, forall_true_left]
    decide
  exact cellCount_eq_zero hiff

private theorem edge_absolute_zero :
    finrank (ZMod 2) (LinearMap.range (edge.toComplex.persistence (K := fun _ => False)
      (S := fun _ => True) (T := fun _ => True) (fun _ _ => trivial) 0)) = 1 := by
  rw [FilteredComplex.finrank_persistence_absolute (fun _ _ _ _ => trivial)
    (fun _ _ _ _ => trivial) (fun _ _ => trivial) 0 edge_reduction]
  have hiff : ∀ i, (edge.toComplex.deg i = 0 ∧ True ∧ (∀ l, ¬ IsLow edge.boundary l i) ∧
      ∀ j, True → ¬ IsLow edge.boundary i j) ↔ i = 0 := by
    simp only [IsLow, true_and, forall_true_left]
    decide
  rw [cellCount_congr hiff]
  exact Nat.card_unique

private theorem edge_absolute_one :
    finrank (ZMod 2) (LinearMap.range (edge.toComplex.persistence (K := fun _ => False)
      (S := fun _ => True) (T := fun _ => True) (fun _ _ => trivial) 1)) = 0 := by
  rw [FilteredComplex.finrank_persistence_absolute (fun _ _ _ _ => trivial)
    (fun _ _ _ _ => trivial) (fun _ _ => trivial) 1 edge_reduction]
  have hiff : ∀ i, ¬ (edge.toComplex.deg i = 1 ∧ True ∧ (∀ l, ¬ IsLow edge.boundary l i) ∧
      ∀ j, True → ¬ IsLow edge.boundary i j) := by
    simp only [IsLow, true_and, forall_true_left]
    decide
  exact cellCount_eq_zero hiff

private theorem edge_absolute_before :
    finrank (ZMod 2) (LinearMap.range (edge.toComplex.persistence (K := fun _ => False)
      (S := (· ≤ 1)) (T := (· ≤ 1)) (fun _ h => h) 0)) = 2 := by
  rw [FilteredComplex.finrank_persistence_absolute (fun _ _ hab hb => hab.trans hb)
    (fun _ _ hab hb => hab.trans hb) (fun _ h => h) 0 edge_reduction]
  have hiff : ∀ i, (edge.toComplex.deg i = 0 ∧ i ≤ 1 ∧ (∀ l, ¬ IsLow edge.boundary l i) ∧
      ∀ j, j ≤ 1 → ¬ IsLow edge.boundary i j) ↔ i ≤ 1 := by
    simp only [IsLow]
    decide
  rw [cellCount_congr hiff, cellCount, Nat.card_eq_fintype_card]
  decide

private theorem edge_absolute_across :
    finrank (ZMod 2) (LinearMap.range (edge.toComplex.persistence (K := fun _ => False)
      (S := (· ≤ 1)) (T := fun _ => True) (fun _ _ => trivial) 0)) = 1 := by
  rw [FilteredComplex.finrank_persistence_absolute (fun _ _ hab hb => hab.trans hb)
    (fun _ _ _ _ => trivial) (fun _ _ => trivial) 0 edge_reduction]
  have hiff : ∀ i, (edge.toComplex.deg i = 0 ∧ i ≤ 1 ∧ (∀ l, ¬ IsLow edge.boundary l i) ∧
      ∀ j, True → ¬ IsLow edge.boundary i j) ↔ i = 0 := by
    simp only [IsLow, forall_true_left]
    decide
  rw [cellCount_congr hiff]
  exact Nat.card_unique

private theorem edge_finite_before :
    finrank (ZMod 2) (LinearMap.range (edge.toComplex.persistence (K := (· ≤ 0))
      (S := (· ≤ 1)) (T := (· ≤ 1)) (fun _ h => h) 0)) = 1 := by
  rw [FilteredComplex.finrank_persistence (fun _ _ hab hb => hab.trans hb)
    (fun _ _ hab hb => hab.trans hb) (fun _ _ hab hb => hab.trans hb)
    (fun _ hb => hb.trans (by decide)) (fun _ h => h) 0 edge_reduction]
  have hiff : ∀ i, FilteredComplex.Alive edge.toComplex edge.boundary (· ≤ 0) (· ≤ 1) (· ≤ 1) 0 i ↔
      i = 1 := by
    simp only [FilteredComplex.Alive, IsLow]
    decide
  rw [cellCount_congr hiff]
  exact Nat.card_unique

private theorem edge_finite_across :
    finrank (ZMod 2) (LinearMap.range (edge.toComplex.persistence (K := (· ≤ 0))
      (S := (· ≤ 1)) (T := fun _ => True) (fun _ _ => trivial) 0)) = 0 := by
  rw [FilteredComplex.finrank_persistence (fun _ _ hab hb => hab.trans hb)
    (fun _ _ hab hb => hab.trans hb) (fun _ _ _ _ => trivial)
    (fun _ hb => hb.trans (by decide)) (fun _ _ => trivial) 0 edge_reduction]
  have hiff : ∀ i, ¬ FilteredComplex.Alive edge.toComplex edge.boundary (· ≤ 0) (· ≤ 1)
      (fun _ => True) 0 i := by
    simp only [FilteredComplex.Alive, IsLow, forall_true_left]
    decide
  exact cellCount_eq_zero hiff

end Example

end SimplicialFiltration

end RelativePersistence
