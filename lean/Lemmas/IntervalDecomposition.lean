import Lemmas.RelativeHomology

/-!
# The persistent homology of a pair is a sum of the table's intervals

`RelativeHomology` shows that the ranks of `H_p(S, K) → H_p(T, K)` count the table's bars. This file
proves the structure behind those counts, with no appeal to the general structure theorem.

`exists_intervalDecomposition` builds it: each `H_p(S, K)` has a basis indexed by the bars alive on
`S`, and the map to `H_p(T, K)` sends a basis element to the same bar's basis element, or to zero if
the bar has died in `T`. That is an isomorphism with a direct sum of interval modules.

`barcode_unique` shows nothing else fits: any family of intervals with the module's ranks has the
table's intervals, with the table's multiplicities.

The representatives (`rep`) are the classical ones. A bar that dies is represented by the boundary
that kills it, not by the column of `V` at its birth: that column's class can outlive the bar, since
the boundary only says it is homologous to older classes, the elder rule. A bar that never dies is
represented by its column of `V`, which for a class moved up from `K` is a relative cycle and not a
cycle. This covers the mathematics, not `TDA.Ph.RelativeFiltration`, which is checked by tests.
-/

namespace RelativePersistence

open Module

section Representative

variable {ι F : Type*} [LinearOrder ι] [Field F]

open scoped Classical in
/-- The chain representing the bar born at `i`: the boundary that kills it if it dies, and its own
column of `V` if not. -/
noncomputable def rep (R V : Matrix ι ι F) (i : ι) : ι → F :=
  if h : ∃ j, IsLow R i j then fun r => R r h.choose else fun r => V r i

theorem rep_of_isLow {R V : Matrix ι ι F} (hR : Reduced R) {i j : ι} (h : IsLow R i j) :
    rep R V i = fun r => R r j := by
  have hex : ∃ j, IsLow R i j := ⟨j, h⟩
  rw [rep, dif_pos hex, hR i _ _ hex.choose_spec h]

theorem rep_of_forall_not {R V : Matrix ι ι F} {i : ι} (h : ∀ j, ¬ IsLow R i j) :
    rep R V i = fun r => V r i := by
  rw [rep, dif_neg (not_exists.2 h)]

open scoped Classical in
/-- The cell that kills the class born at `i`, if any. -/
noncomputable def death (R : Matrix ι ι F) (i : ι) : Option ι :=
  if h : ∃ j, IsLow R i j then some h.choose else none

theorem death_eq_some_iff {R : Matrix ι ι F} (hR : Reduced R) {i j : ι} :
    death R i = some j ↔ IsLow R i j := by
  constructor
  · intro h
    unfold death at h
    split_ifs at h with hex
    rw [← Option.some.inj h]
    exact hex.choose_spec
  · intro h
    have hex : ∃ j, IsLow R i j := ⟨j, h⟩
    rw [death, dif_pos hex, hR i _ _ hex.choose_spec h]

theorem death_eq_none_iff {R : Matrix ι ι F} {i : ι} : death R i = none ↔ ∀ j, ¬ IsLow R i j := by
  unfold death
  split_ifs with hex
  · simpa using hex
  · simpa using hex

end Representative

variable {ι F : Type*} [LinearOrder ι] [Fintype ι] [Field F]

namespace FilteredComplex

variable {C : FilteredComplex ι F} {K S T : ι → Prop} {p : ℕ} {R V U : Matrix ι ι F}

theorem rep_self (g : C.GradedReduction R V) (i : ι) : rep R V i i ≠ 0 := by
  by_cases hex : ∃ j, IsLow R i j
  · obtain ⟨j, hj⟩ := hex
    rw [rep_of_isLow g.reduced hj]
    exact hj.1
  · push Not at hex
    rw [rep_of_forall_not hex]
    exact g.diag i

theorem rep_above (g : C.GradedReduction R V) {i r : ι} (h : i < r) : rep R V i r = 0 := by
  by_cases hex : ∃ j, IsLow R i j
  · obtain ⟨j, hj⟩ := hex
    rw [rep_of_isLow g.reduced hj]
    exact hj.2 r h
  · push Not at hex
    rw [rep_of_forall_not hex]
    exact g.upper r i h

/-- A bar's representative is a relative cycle on any `S` it is born on. -/
theorem rep_mem_relCycles (g : C.GradedReduction R V) (hK : IsInitial K) (hS : IsInitial S)
    {i : ι} (hi : Alive C R K S T p i) : rep R V i ∈ C.relCycles K S p := by
  obtain ⟨hd, hSi, -, hl, -⟩ := hi
  refine ⟨fun r hr => ?_, ?_⟩
  · by_contra hne
    have hri : r ≤ i := by_contra fun h => hne (rep_above g (not_le.mp h))
    refine hr ⟨hS r i hri hSi, ?_⟩
    by_cases hex : ∃ j, IsLow R i j
    · obtain ⟨j, hj⟩ := hex
      rw [rep_of_isLow g.reduced hj] at hne
      have h₁ := g.deg hne
      have h₂ := g.deg hj.1
      omega
    · push Not at hex
      rw [rep_of_forall_not hex] at hne
      exact (g.homog r i hne).trans hd
  · change C.boundary.mulVec (rep R V i) ∈ chains F K
    by_cases hex : ∃ j, IsLow R i j
    · obtain ⟨j, hj⟩ := hex
      rw [rep_of_isLow g.reduced hj]
      have hzero : C.boundary.mulVec (fun r => R r j) = 0 :=
        funext fun q => congrFun (congrFun g.boundary_mul q) j
      rw [hzero]
      exact Submodule.zero_mem _
    · push Not at hex
      rw [rep_of_forall_not hex]
      have hcol : C.boundary.mulVec (fun r => V r i) = fun q => R q i :=
        funext fun q => congrFun (congrFun g.mul q) i
      rw [hcol]
      exact fun q hq => column_eq_zero hK (not_lowOutside_iff.2 hl) q hq

/-- Once its killer is on `T`, a bar's representative is a boundary there. -/
theorem rep_mem_relBoundaries (g : C.GradedReduction R V) (hT : IsInitial T) {i j : ι} (hj : T j)
    (hlow : IsLow R i j) (hd : C.deg i = p) : rep R V i ∈ C.relBoundaries K T p := by
  rw [rep_of_isLow g.reduced hlow]
  refine Submodule.mem_sup_left ⟨fun r => V r j, fun r hr => ?_, ?_⟩
  · by_contra hne
    have hdj := g.deg hlow.1
    exact hr (cells_closed g hT (p + 1) r j ⟨hj, by omega⟩ hne)
  · exact funext fun q => congrFun (congrFun g.mul q) j

/-- The class in `H_p(S, K)` of a bar alive on `S`. -/
noncomputable def barClass (g : C.GradedReduction R V) (hK : IsInitial K) (hS : IsInitial S)
    (p : ℕ) (i : {i // Alive C R K S S p i}) : C.RelHomology K S p :=
  Submodule.mkQ _ ⟨rep R V i.1, rep_mem_relCycles g hK hS i.2⟩

/-- The classes of the bars alive on `S` are independent. Each representative's last nonzero entry
is at its own birth, so a combination ends at the latest bar it uses. If the combination were a
relative boundary, that entry would be the low of a column on `S`, and the bar would be dead. -/
theorem barClass_linearIndependent (g : C.GradedReduction R V) (hK : IsInitial K)
    (hS : IsInitial S) (p : ℕ) : LinearIndependent F (barClass g hK hS p) := by
  classical
  haveI : Fintype {i // Alive C R K S S p i} := Fintype.ofFinite _
  rw [Fintype.linearIndependent_iff]
  intro c hc
  by_contra hne
  push Not at hne
  obtain ⟨i₀, hi₀⟩ := hne
  obtain ⟨m, hm, hmax⟩ := (Finset.univ.filter fun i => c i ≠ 0).exists_max_image
    (fun i : {i // Alive C R K S S p i} => i.1) ⟨i₀, by simp [hi₀]⟩
  have hcm : c m ≠ 0 := (Finset.mem_filter.1 hm).2
  have hle : ∀ i, c i ≠ 0 → i.1 ≤ m.1 := fun i hi => hmax i (by simp [hi])
  obtain ⟨-, -, hmK, -, hmk⟩ := m.2
  set z : C.relCycles K S p := ∑ i, c i • ⟨rep R V i.1, rep_mem_relCycles g hK hS i.2⟩ with hz
  have hmkQ : ((C.relBoundaries K S p).comap (C.relCycles K S p).subtype).mkQ z = 0 := by
    rw [hz, map_sum, ← hc]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [map_smul]
    rfl
  have hzB : (z : ι → F) ∈ C.relBoundaries K S p :=
    (Submodule.Quotient.mk_eq_zero
      ((C.relBoundaries K S p).comap (C.relCycles K S p).subtype)).1 hmkQ
  have hzval : ∀ r, (z : ι → F) r = ∑ i, c i * rep R V i.1 r := by
    intro r
    simp [hz, Finset.sum_apply]
  have hzm : (z : ι → F) m.1 ≠ 0 := by
    rw [hzval, Finset.sum_eq_single m]
    · exact mul_ne_zero hcm (rep_self g m.1)
    · intro i _ hi
      by_cases hci : c i = 0
      · simp [hci]
      · have hlt : i.1 < m.1 := lt_of_le_of_ne (hle i hci) fun h => hi (Subtype.ext h)
        rw [rep_above g hlt, mul_zero]
    · simp
  have hzabove : ∀ r, m.1 < r → (z : ι → F) r = 0 := by
    intro r hr
    rw [hzval]
    refine Finset.sum_eq_zero fun i _ => ?_
    by_cases hci : c i = 0
    · simp [hci]
    · rw [rep_above g ((hle i hci).trans_lt hr), mul_zero]
  obtain ⟨y, ⟨w', hw', rfl⟩, w, hw, hsum⟩ := Submodule.mem_sup.1 hzB
  rw [← chains_map_eq (F := F) g.isUnit (cells_closed g hS (p + 1))] at hw'
  obtain ⟨c', hc', rfl⟩ := hw'
  have hbd : C.boundary.mulVecLin (V.mulVecLin c') = R.mulVec c' := by
    change C.boundary.mulVec (V.mulVec c') = R.mulVec c'
    rw [Matrix.mulVec_mulVec, g.mul]
  rw [hbd] at hsum
  have hR : ∀ r, m.1 ≤ r → R.mulVec c' r = (z : ι → F) r := by
    intro r hr
    have hwr : w r = 0 := hw r fun h => hmK (hK _ _ hr h.1)
    have hval := congrFun hsum r
    rw [Pi.add_apply, hwr, add_zero] at hval
    exact hval
  have hRm : ∑ s, R m.1 s * c' s ≠ 0 := by
    have hval := hR m.1 le_rfl
    rw [← hval] at hzm
    exact hzm
  have hRabove : ∀ r, m.1 < r → ∑ s, R r s * c' s = 0 := by
    intro r hr
    have hval := hR r hr.le
    rw [hzabove r hr] at hval
    exact hval
  obtain ⟨s, _, hs⟩ := Finset.exists_ne_zero_of_sum_ne_zero hRm
  obtain ⟨q, l, hl, hql, hq, htop⟩ :=
    combination_top g.reduced c' ⟨s, right_ne_zero_of_mul hs, m.1, left_ne_zero_of_mul hs⟩
  have hqm : q = m.1 := by
    rcases lt_trichotomy q m.1 with h | h | h
    · exact absurd (htop m.1 h) hRm
    · exact h
    · exact absurd (hRabove q h) hq
  rw [hqm] at hql
  have hcl : C.cells S (p + 1) l := by_contra fun h => hl (hc' l h)
  exact hmk l hcl.1 hql

/-- The map from a space to itself is the identity. -/
theorem persistence_self (p : ℕ) :
    C.persistence (K := K) (S := S) (T := S) (fun _ h => h) p = LinearMap.id :=
  Submodule.linearMap_qext _ (LinearMap.ext fun _ => rfl)

/-- `H_p(S, K)` has one dimension per bar alive on `S`. -/
theorem finrank_relHomology (g : C.GradedReduction R V) (hK : IsInitial K) (hS : IsInitial S)
    (hKS : ∀ i, K i → S i) (p : ℕ) :
    finrank F (C.RelHomology K S p) = cellCount (Alive C R K S S p) := by
  rw [← finrank_persistence_of_graded g hK hS hS hKS (fun _ h => h) p, persistence_self,
    LinearMap.range_id, finrank_top]

/-- The basis of `H_p(S, K)` whose elements are the classes of the bars alive on `S`. -/
noncomputable def barBasis (g : C.GradedReduction R V) (hK : IsInitial K) (hS : IsInitial S)
    (hKS : ∀ i, K i → S i) (p : ℕ) :
    Basis {i // Alive C R K S S p i} F (C.RelHomology K S p) :=
  letI : Fintype {i // Alive C R K S S p i} := Fintype.ofFinite _
  Basis.mk (barClass_linearIndependent g hK hS p)
    ((barClass_linearIndependent g hK hS p).span_eq_top_of_card_eq_finrank' (by
      rw [finrank_relHomology g hK hS hKS p, cellCount, Nat.card_eq_fintype_card])).ge

theorem barBasis_apply (g : C.GradedReduction R V) (hK : IsInitial K) (hS : IsInitial S)
    (hKS : ∀ i, K i → S i) (p : ℕ) (i : {i // Alive C R K S S p i}) :
    barBasis g hK hS hKS p i = barClass g hK hS p i :=
  Basis.mk_apply _ _ _

open scoped Classical in
/-- The maps act diagonally: a bar's basis element goes to the same bar's basis element while it
lives, and to zero once it has died. -/
theorem persistence_barBasis (g : C.GradedReduction R V) (hK : IsInitial K) (hS : IsInitial S)
    (hT : IsInitial T) (hKS : ∀ i, K i → S i) (hKT : ∀ i, K i → T i) (hST : ∀ i, S i → T i)
    (p : ℕ) (i : {i // Alive C R K S S p i}) :
    C.persistence hST p (barBasis g hK hS hKS p i) =
      if h : Alive C R K T T p i.1 then barBasis g hK hT hKT p ⟨i.1, h⟩ else 0 := by
  rw [barBasis_apply]
  split_ifs with h
  · rw [barBasis_apply]
    rfl
  · obtain ⟨hd, hSi, hKi, hl, _⟩ := i.2
    have hkilled : ∃ j, T j ∧ IsLow R i.1 j := by
      by_contra hne
      push Not at hne
      exact h ⟨hd, hST _ hSi, hKi, hl, hne⟩
    obtain ⟨j, hTj, hlow⟩ := hkilled
    exact (Submodule.Quotient.mk_eq_zero _).2 (rep_mem_relBoundaries g hT hTj hlow hd)

open scoped Classical in
/-- **The persistent homology of the pair is a sum of the table's intervals.** For any reduction of
the boundary, each `H_p(S, K)` has a basis indexed by the bars alive on `S`, and the map to
`H_p(T, K)` sends a bar's basis element to the same bar's, or to zero if it has died in `T`. -/
theorem exists_intervalDecomposition (hK : IsInitial K) (p : ℕ) (h : IsReduction C.boundary R U) :
    ∃ b : ∀ S : ι → Prop, IsInitial S → (∀ i, K i → S i) →
        Basis {i // Alive C R K S S p i} F (C.RelHomology K S p),
      ∀ (S T : ι → Prop) (hS : IsInitial S) (hKS : ∀ i, K i → S i) (hT : IsInitial T)
        (hKT : ∀ i, K i → T i) (hST : ∀ i, S i → T i) (i : {i // Alive C R K S S p i}),
        C.persistence hST p (b S hS hKS i) =
          if hi : Alive C R K T T p i.1 then b T hT hKT ⟨i.1, hi⟩ else 0 := by
  obtain ⟨Rg, Vg, g⟩ := C.exists_gradedReduction
  have halive : ∀ S T i, Alive C Rg K S T p i ↔ Alive C R K S T p i := by
    intro S T i
    simp only [Alive, isLow_iff_of_reductions h g.toIsReduction]
  refine ⟨fun S hS hKS => (barBasis g hK hS hKS p).reindex
    (Equiv.subtypeEquivRight fun i => halive S S i), ?_⟩
  intro S T hS hKS hT hKT hST i
  rw [Basis.reindex_apply]
  by_cases hi : Alive C R K T T p i.1
  · rw [dif_pos hi, Basis.reindex_apply, persistence_barBasis g hK hS hT hKS hKT hST, dif_pos]
    rfl
  · rw [dif_neg hi, persistence_barBasis g hK hS hT hKS hKT hST, dif_neg]
    exact fun hg => hi ((halive T T i.1).1 hg)

section Unique

variable {B B' : Type*} [Finite B] [Finite B']

/-- How many of a family of intervals are alive from `S` through `T`: born on `S`, with no death on
`T`. -/
noncomputable def intervalRank (b : B → ι) (d : B → Option ι) (S T : ι → Prop) : ℕ :=
  cellCount fun x => S (b x) ∧ ∀ j, d x = some j → ¬ T j

omit [Fintype ι] in
/-- Moving the left end past `i` picks out the intervals born at `i`. -/
theorem intervalRank_born (b : B → ι) (d : B → Option ι) (i : ι) (T : ι → Prop) :
    intervalRank b d (· ≤ i) T = intervalRank b d (· < i) T +
      cellCount fun x => b x = i ∧ ∀ j, d x = some j → ¬ T j := by
  have hsplit := cellCount_split (fun x => b x ≤ i ∧ ∀ j, d x = some j → ¬ T j) (fun x => b x < i)
  have h₁ : cellCount (fun x => (b x ≤ i ∧ ∀ j, d x = some j → ¬ T j) ∧ b x < i) =
      cellCount fun x => b x < i ∧ ∀ j, d x = some j → ¬ T j :=
    cellCount_congr fun x => ⟨fun h => ⟨h.2, h.1.2⟩, fun h => ⟨⟨h.1.le, h.2⟩, h.1⟩⟩
  have h₂ : cellCount (fun x => (b x ≤ i ∧ ∀ j, d x = some j → ¬ T j) ∧ ¬ b x < i) =
      cellCount fun x => b x = i ∧ ∀ j, d x = some j → ¬ T j :=
    cellCount_congr fun x => ⟨fun h => ⟨le_antisymm h.1.1 (not_lt.mp h.2), h.1.2⟩,
      fun h => ⟨⟨h.1.le, h.2⟩, h.1.not_lt⟩⟩
  rw [h₁, h₂] at hsplit
  exact hsplit.symm

omit [Fintype ι] in
/-- Moving the right end past `j` picks out the intervals that die at `j`. -/
theorem intervalCount_died (b : B → ι) (d : B → Option ι) (i j : ι) :
    (cellCount fun x => b x = i ∧ ∀ q, d x = some q → ¬ q < j) =
      (cellCount fun x => b x = i ∧ ∀ q, d x = some q → ¬ q ≤ j) +
        cellCount fun x => b x = i ∧ d x = some j := by
  have hsplit := cellCount_split (fun x => b x = i ∧ ∀ q, d x = some q → ¬ q < j)
    (fun x => d x = some j)
  have h₁ : cellCount (fun x => (b x = i ∧ ∀ q, d x = some q → ¬ q < j) ∧ d x = some j) =
      cellCount fun x => b x = i ∧ d x = some j := by
    refine cellCount_congr fun x => ⟨fun h => ⟨h.1.1, h.2⟩, fun h => ⟨⟨h.1, fun q hq => ?_⟩, h.2⟩⟩
    have hqj : q = j := Option.some.inj (hq.symm.trans h.2)
    rw [hqj]
    exact lt_irrefl j
  have h₂ : cellCount (fun x => (b x = i ∧ ∀ q, d x = some q → ¬ q < j) ∧ ¬ d x = some j) =
      cellCount fun x => b x = i ∧ ∀ q, d x = some q → ¬ q ≤ j := by
    refine cellCount_congr fun x => ⟨fun h => ⟨h.1.1, fun q hq hle => ?_⟩, fun h => ⟨⟨h.1,
      fun q hq hlt => h.2 q hq hlt.le⟩, fun hd => h.2 j hd le_rfl⟩⟩
    have hge : j ≤ q := not_lt.mp (h.1.2 q hq)
    have hqj : q = j := le_antisymm hle hge
    exact h.2 (hqj ▸ hq)
  rw [h₁, h₂] at hsplit
  omega

omit [LinearOrder ι] [Fintype ι] [Finite B] in
/-- An interval with no death on the whole complex is one that never dies. -/
theorem intervalCount_essential (b : B → ι) (d : B → Option ι) (i : ι) :
    (cellCount fun x => b x = i ∧ ∀ j, d x = some j → ¬ True) =
      cellCount fun x => b x = i ∧ d x = none := by
  refine cellCount_congr fun x => and_congr_right fun _ => ?_
  rw [Option.eq_none_iff_forall_ne_some]
  exact ⟨fun h j hj => h j hj trivial, fun h j hj _ => h j hj⟩

omit [Fintype ι] in
/-- Intervals are determined by their ranks. Two families of intervals with the same ranks agree on
every interval born outside `K`, with multiplicity. An empty interval has no ranks to show it, which
is why the finite case asks that the interval run forward. -/
theorem interval_count_eq (hK : IsInitial K) {b : B → ι} {d : B → Option ι} {b' : B' → ι}
    {d' : B' → Option ι} (hr : ∀ S T : ι → Prop, IsInitial S → IsInitial T → (∀ i, K i → S i) →
      (∀ i, S i → T i) → intervalRank b d S T = intervalRank b' d' S T) {i : ι} (hi : ¬ K i) :
    (∀ j, i < j → (cellCount fun x => b x = i ∧ d x = some j) =
      cellCount fun x => b' x = i ∧ d' x = some j) ∧
    (cellCount fun x => b x = i ∧ d x = none) = cellCount fun x => b' x = i ∧ d' x = none := by
  have hlt : ∀ x, K x → x < i := fun x hx => by_contra fun h => hi (hK i x (not_lt.mp h) hx)
  have hle : ∀ a : ι, IsInitial (· ≤ a) := fun a _ _ hab hb => hab.trans hb
  have hlt' : ∀ a : ι, IsInitial (· < a) := fun a _ _ hab hb => hab.trans_lt hb
  constructor
  · intro j hij
    have e₁ := hr (· ≤ i) (· < j) (hle i) (hlt' j) (fun x hx => (hlt x hx).le)
      (fun x hx => lt_of_le_of_lt hx hij)
    have e₂ := hr (· < i) (· < j) (hlt' i) (hlt' j) hlt (fun x hx => hx.trans hij)
    have e₃ := hr (· ≤ i) (· ≤ j) (hle i) (hle j) (fun x hx => (hlt x hx).le)
      (fun x hx => le_trans hx hij.le)
    have e₄ := hr (· < i) (· ≤ j) (hlt' i) (hle j) hlt (fun x hx => (hx.trans hij).le)
    have a₁ := intervalRank_born b d i (· < j)
    have a₂ := intervalRank_born b' d' i (· < j)
    have a₃ := intervalRank_born b d i (· ≤ j)
    have a₄ := intervalRank_born b' d' i (· ≤ j)
    have d₁ := intervalCount_died b d i j
    have d₂ := intervalCount_died b' d' i j
    omega
  · have e₁ := hr (· ≤ i) (fun _ => True) (hle i) (fun _ _ _ _ => trivial)
      (fun x hx => (hlt x hx).le) (fun _ _ => trivial)
    have e₂ := hr (· < i) (fun _ => True) (hlt' i) (fun _ _ _ _ => trivial) hlt (fun _ _ => trivial)
    have a₁ := intervalRank_born b d i (fun _ => True)
    have a₂ := intervalRank_born b' d' i (fun _ => True)
    rw [intervalCount_essential b d i] at a₁
    rw [intervalCount_essential b' d' i] at a₂
    omega

/-- Cell `i` starts a bar of `H_p(·, K)`: of degree `p`, outside `K`, with its column empty or its
low in `K`. Its interval runs from `i` to the cell that kills it, if any. -/
def IsBar (C : FilteredComplex ι F) (R : Matrix ι ι F) (K : ι → Prop) (p : ℕ) (i : ι) : Prop :=
  C.deg i = p ∧ ¬ K i ∧ ∀ l, IsLow R l i → K l

/-- A bar's interval: it is alive from `S` through `T` exactly when it is born on `S` and its
killer, if it has one, is not on `T`. -/
theorem alive_iff_isBar (hR : Reduced R) {S T : ι → Prop} {i : ι} :
    Alive C R K S T p i ↔ IsBar C R K p i ∧ S i ∧ ∀ j, death R i = some j → ¬ T j := by
  constructor
  · rintro ⟨hd, hSi, hKi, hl, hk⟩
    exact ⟨⟨hd, hKi, hl⟩, hSi, fun j hj hT => hk j hT ((death_eq_some_iff hR).1 hj)⟩
  · rintro ⟨⟨hd, hKi, hl⟩, hSi, hk⟩
    exact ⟨hd, hSi, hKi, hl, fun j hj hlow => hk j ((death_eq_some_iff hR).2 hlow) hj⟩

/-- The table as a family of intervals has the module's ranks. -/
theorem intervalRank_table (hR : Reduced R) (S T : ι → Prop) :
    intervalRank (Subtype.val : {i // IsBar C R K p i} → ι) (fun x => death R x.1) S T =
      cellCount (Alive C R K S T p) := by
  refine cellCount_bij (fun x => x.1) ?_ (fun x y _ _ h => Subtype.ext h) ?_
  · rintro ⟨i, hi⟩ ⟨hSi, hk⟩
    exact (alive_iff_isBar hR).2 ⟨hi, hSi, hk⟩
  · intro i hi
    obtain ⟨hbar, hSi, hk⟩ := (alive_iff_isBar hR).1 hi
    exact ⟨⟨i, hbar⟩, ⟨hSi, hk⟩, rfl⟩

/-- **The barcode is the table, and nothing else.** Any family of intervals with the ranks of the
persistent homology of the pair has exactly the table's intervals born outside `K`. An interval
`[i, j)` occurs once if `i` starts a bar that `j` kills, and `[i, ∞)` occurs once if nothing kills
it. Neither occurs otherwise. -/
theorem barcode_unique (hK : IsInitial K) (p : ℕ) (h : IsReduction C.boundary R U) {b : B → ι}
    {d : B → Option ι} (hr : ∀ S T : ι → Prop, IsInitial S → IsInitial T → (∀ i, K i → S i) →
      ∀ hST : (∀ i, S i → T i), finrank F (LinearMap.range (C.persistence (K := K) hST p)) =
        intervalRank b d S T) {i : ι} (hi : ¬ K i) :
    (∀ j, i < j → (cellCount fun x => b x = i ∧ d x = some j) =
      cellCount fun x : ι => x = i ∧ IsBar C R K p x ∧ IsLow R x j) ∧
    (cellCount fun x => b x = i ∧ d x = none) =
      cellCount fun x : ι => x = i ∧ IsBar C R K p x ∧ ∀ j, ¬ IsLow R x j := by
  obtain ⟨hfin, hess⟩ := interval_count_eq hK (b' := (Subtype.val : {i // IsBar C R K p i} → ι))
    (d' := fun x => death R x.1)
    (fun S T hS hT hKS hST => by
      rw [← hr S T hS hT hKS hST, finrank_persistence hK hS hT hKS hST p h,
        intervalRank_table h.reduced]) hi
  constructor
  · intro j hij
    rw [hfin j hij]
    refine cellCount_bij (fun x => x.1) ?_ (fun x y _ _ e => Subtype.ext e) ?_
    · rintro ⟨x, hx⟩ ⟨hxi, hdx⟩
      exact ⟨hxi, hx, (death_eq_some_iff h.reduced).1 hdx⟩
    · rintro x ⟨hxi, hx, hlow⟩
      exact ⟨⟨x, hx⟩, ⟨hxi, (death_eq_some_iff h.reduced).2 hlow⟩, rfl⟩
  · rw [hess]
    refine cellCount_bij (fun x => x.1) ?_ (fun x y _ _ e => Subtype.ext e) ?_
    · rintro ⟨x, hx⟩ ⟨hxi, hdx⟩
      exact ⟨hxi, hx, death_eq_none_iff.1 hdx⟩
    · rintro x ⟨hxi, hx, hn⟩
      exact ⟨⟨x, hx⟩, ⟨hxi, death_eq_none_iff.2 hn⟩, rfl⟩

end Unique

end FilteredComplex

end RelativePersistence
