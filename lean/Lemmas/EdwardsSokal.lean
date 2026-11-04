import Lemmas.RandomCluster
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.SetTheory.Cardinal.Finite
import Mathlib.Tactic.FieldSimp

/-!
# The Potts model and the Edwards–Sokal two-point identity

The Potts model on a finite graph with a coupling per edge, and `edwards_sokal_twopoint`, Grimmett
2006, Thm 1.16 and (11.25): spin agreement is `1 / q` plus `1 - 1 / q` times the connection
probability under the random-cluster measure.

The proof is the coupling, done as a double sum over spins and configurations. Each pair's factor
splits as `1 - p_e + p_e δ_e`, and expanding the product gives a sum over configurations in which a
spin assignment counts only if it is constant on clusters (`pottsWeight_eq_sum`). Swapping the sums,
`n ^ k(ω)` assignments are constant on the clusters of `ω`, which is the random-cluster weight
(`card_compat`). All of them agree across a connected pair, and a fraction `1 / n` of them across a
disconnected one (`card_compat_agree_mul`).

`pottsWeight_pe` ties the weight to the physics: it is `exp (-β H)` for `H = ∑_e J_e (1 - δ_e)`.

All of this is about the exact expectation. A finite-sample estimate need not be monotone.
-/

universe u

namespace Bifiltration

variable {V : Type u}

/-! ### The Potts model -/

/-- The graph of the pairs on which `σ` agrees. -/
def agreeGraph {n : ℕ} (σ : V → Fin n) : SimpleGraph V where
  Adj x y := x ≠ y ∧ σ x = σ y
  symm := ⟨fun _ _ h => ⟨h.1.symm, h.2.symm⟩⟩
  loopless := ⟨fun _ h => h.1 rfl⟩

open Classical in
/-- The unordered pairs of distinct vertices. -/
noncomputable abbrev edgePairs (V : Type u) [Fintype V] : Finset (Sym2 V) :=
  Finset.univ.filter (fun e : Sym2 V => ¬ e.IsDiag)

open Classical in
/-- The Boltzmann weight of `σ`, with `exp (-β J_e) = 1 - p_e`: a pair contributes 1 if `σ` agrees
across it and `1 - p_e` if not. -/
noncomputable def pottsWeight [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) {n : ℕ}
    (σ : V → Fin n) : ℝ :=
  ∏ e ∈ edgePairs V, if e ∈ (agreeGraph σ).edgeSet then 1 else 1 - edgeParam hp e

open Classical in
/-- The Potts partition function. -/
noncomputable def pottsPartition [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) (n : ℕ) : ℝ :=
  ∑ σ : V → Fin n, pottsWeight hp σ

open Classical in
/-- The probability that `i` and `j` agree. -/
noncomputable def pottsAgree [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) (n : ℕ)
    (i j : V) : ℝ :=
  ∑ σ : V → Fin n, if σ i = σ j then pottsWeight hp σ / pottsPartition hp n else 0

open Classical in
/-- The Potts Hamiltonian `H(σ) = ∑_e J_e (1 - δ_e(σ))`. -/
noncomputable def pottsHamiltonian [Fintype V] (J : V → V → ℝ) (hJs : ∀ i j, J i j = J j i) {n : ℕ}
    (σ : V → Fin n) : ℝ :=
  ∑ e ∈ edgePairs V, if e ∈ (agreeGraph σ).edgeSet then 0 else Sym2.lift ⟨J, hJs⟩ e

/-- At `p_e = 1 - exp (-β J_e)` the weight is the Boltzmann factor of the Hamiltonian. -/
theorem pottsWeight_pe [Fintype V] (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j)
    (hJs : ∀ i j, J i j = J j i) {β : ℝ} (hβ : 0 ≤ β) {n : ℕ} (σ : V → Fin n) :
    pottsWeight (pe_isEdgeParameter J hJ hJs hβ) σ =
      Real.exp (-(β * pottsHamiltonian J hJs σ)) := by
  classical
  rw [pottsHamiltonian, Finset.mul_sum, ← Finset.sum_neg_distrib, Real.exp_sum, pottsWeight]
  refine Finset.prod_congr rfl fun e _ => ?_
  by_cases he : e ∈ (agreeGraph σ).edgeSet
  · simp [he]
  · induction e using Sym2.ind with
    | _ x y => simp [he, edgeParam_mk, pe]

/-! ### Expanding the product into configurations -/

open Classical in
/-- A product over the pairs of a sum of two terms is a sum over configurations: each configuration
takes the first term on its open edges and the second on the rest. -/
theorem prod_add_eq_sum_config [Fintype V] (g h : Sym2 V → ℝ) :
    ∏ e ∈ edgePairs V, (g e + h e) =
      ∑ ω : Config V, ∏ e ∈ edgePairs V, if e ∈ ω.edgeSet then g e else h e := by
  rw [Finset.prod_add]
  have hset : ∀ t ∈ (edgePairs V).powerset, ∀ e,
      e ∈ (SimpleGraph.fromEdgeSet (t : Set (Sym2 V))).edgeSet ↔ e ∈ t := by
    intro t ht e
    have hsub := Finset.mem_powerset.mp ht
    rw [SimpleGraph.edgeSet_fromEdgeSet]
    constructor
    · exact fun he => by simpa using he.1
    · intro he
      exact ⟨by simpa using he, by simpa using (Finset.mem_filter.mp (hsub he)).2⟩
  have hfilter : ∀ (ω : Config V) (e : Sym2 V),
      e ∈ (edgePairs V).filter (· ∈ ω.edgeSet) ↔ e ∈ edgePairs V ∧ e ∈ ω.edgeSet :=
    fun _ _ => Finset.mem_filter
  refine Finset.sum_nbij' (fun t => SimpleGraph.fromEdgeSet (t : Set (Sym2 V)))
    (fun ω => (edgePairs V).filter (· ∈ ω.edgeSet))
    (fun _ _ => Finset.mem_univ _)
    (fun ω _ => Finset.mem_powerset.mpr (Finset.filter_subset _ _)) ?_ ?_ ?_
  · intro t ht
    ext e
    rw [hfilter, hset t ht e]
    exact ⟨fun h => h.2, fun he => ⟨Finset.mem_powerset.mp ht he, he⟩⟩
  · intro ω _
    ext x y
    rw [SimpleGraph.fromEdgeSet_adj]
    simp only [Finset.coe_filter, Set.mem_setOf_eq, SimpleGraph.mem_edgeSet]
    constructor
    · exact fun h => h.1.2
    · intro h
      exact ⟨⟨Finset.mem_filter.mpr ⟨Finset.mem_univ _, by simpa using h.ne⟩, h⟩, h.ne⟩
  · intro t ht
    have hsub := Finset.mem_powerset.mp ht
    have hprod : (∏ i ∈ t, g i) * ∏ i ∈ edgePairs V \ t, h i =
        ∏ e ∈ edgePairs V, if e ∈ t then g e else h e := by
      rw [Finset.prod_ite, Finset.filter_mem_eq_inter, Finset.inter_eq_right.mpr hsub,
        Finset.filter_notMem_eq_sdiff]
    rw [hprod]
    refine Finset.prod_congr rfl fun e _ => ?_
    by_cases he : e ∈ t
    · rw [if_pos he, if_pos ((hset t ht e).mpr he)]
    · rw [if_neg he, if_neg fun h' => he ((hset t ht e).mp h')]

open Classical in
/-- The coupling: a spin assignment's weight is the sum, over the configurations it is constant
along, of their edge products. -/
theorem pottsWeight_eq_sum [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) {n : ℕ}
    (σ : V → Fin n) :
    pottsWeight hp σ = ∑ ω : Config V,
      (∏ e ∈ edgePairs V, edgeFactor hp ω e) * if ω ≤ agreeGraph σ then 1 else 0 := by
  have hsplit : ∀ e, (if e ∈ (agreeGraph σ).edgeSet then (1 : ℝ) else 1 - edgeParam hp e) =
      (if e ∈ (agreeGraph σ).edgeSet then edgeParam hp e else 0) + (1 - edgeParam hp e) := by
    intro e
    split_ifs <;> ring
  simp only [pottsWeight, hsplit]
  rw [prod_add_eq_sum_config]
  refine Finset.sum_congr rfl fun ω _ => ?_
  by_cases hω : ω ≤ agreeGraph σ
  · rw [if_pos hω, mul_one]
    refine Finset.prod_congr rfl fun e _ => ?_
    by_cases he : e ∈ ω.edgeSet
    · have hag : e ∈ (agreeGraph σ).edgeSet := SimpleGraph.edgeSet_mono hω he
      simp [edgeFactor, he, hag]
    · simp [edgeFactor, he]
  · rw [if_neg hω, mul_zero]
    obtain ⟨e, he, hne⟩ := Set.not_subset.mp
      fun hsub => hω (SimpleGraph.edgeSet_subset_edgeSet.mp hsub)
    refine Finset.prod_eq_zero (i := e) ?_ ?_
    · simpa using ω.not_isDiag_of_mem_edgeSet he
    · simp [he, hne]

open Classical in
/-- Swapping the sums: any function of the spins, averaged against the weight, becomes a sum over
configurations of a count over the spins constant along them. -/
theorem sum_spin_mul_weight [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) (n : ℕ)
    (A : (V → Fin n) → ℝ) :
    ∑ σ : V → Fin n, A σ * pottsWeight hp σ = ∑ ω : Config V,
      (∏ e ∈ edgePairs V, edgeFactor hp ω e) *
        ∑ σ : V → Fin n, A σ * if ω ≤ agreeGraph σ then 1 else 0 := by
  simp only [pottsWeight_eq_sum, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun ω _ => Finset.sum_congr rfl fun σ _ => ?_
  ring

/-! ### Counting the spins constant on clusters -/

theorem le_agreeGraph_iff {n : ℕ} (ω : SimpleGraph V) (σ : V → Fin n) :
    ω ≤ agreeGraph σ ↔ ∀ x y, ω.Adj x y → σ x = σ y :=
  ⟨fun h _ _ hxy => (h hxy).2, fun h x y hxy => ⟨hxy.ne, h x y hxy⟩⟩

theorem eq_of_reachable {n : ℕ} {ω : SimpleGraph V} {σ : V → Fin n} (h : ω ≤ agreeGraph σ) {x y : V}
    (hr : ω.Reachable x y) : σ x = σ y := by
  obtain ⟨w⟩ := hr
  induction w with
  | nil => rfl
  | cons hadj _ ih => exact ((le_agreeGraph_iff ω σ).mp h _ _ hadj).trans ih

/-- Spin assignments constant along `ω` are functions of its clusters. -/
def spinsEquiv {n : ℕ} (ω : SimpleGraph V) :
    (ω.ConnectedComponent → Fin n) ≃ {σ : V → Fin n // ω ≤ agreeGraph σ} where
  toFun f := ⟨f ∘ ω.connectedComponentMk, (le_agreeGraph_iff ω _).mpr fun _ _ h =>
    congrArg f (SimpleGraph.ConnectedComponent.connectedComponentMk_eq_of_adj h)⟩
  invFun σ := SimpleGraph.ConnectedComponent.lift σ.1 fun _ _ w _ => eq_of_reachable σ.2 ⟨w⟩
  left_inv f := by
    funext c
    induction c using SimpleGraph.ConnectedComponent.ind with
    | h v => rfl
  right_inv _ := rfl

/-- There are `n ^ k(ω)` spin assignments constant along `ω`. -/
theorem card_compat [Finite V] {n : ℕ} (ω : SimpleGraph V) :
    Nat.card {σ : V → Fin n // ω ≤ agreeGraph σ} = n ^ componentCount ω := by
  rw [← Nat.card_congr (spinsEquiv ω), Nat.card_fun, componentCount]
  simp

/-- Functions agreeing at two distinct points are a fraction `1 / |β|` of all functions: a value at
`c'` and a function agreeing there determine any function, and conversely. -/
theorem card_agree_mul {C β : Type*} [Finite C] [Finite β] {c c' : C} (h : c ≠ c') :
    Nat.card {f : C → β // f c = f c'} * Nat.card β = Nat.card β ^ Nat.card C := by
  classical
  rw [← Nat.card_prod, ← Nat.card_fun]
  refine Nat.card_congr
    { toFun := fun x => Function.update x.1.1 c' x.2
      invFun := fun g =>
        (⟨Function.update g c' (g c), by simp [Function.update_of_ne h]⟩, g c')
      left_inv := ?_
      right_inv := ?_ }
  · rintro ⟨⟨f, hf⟩, b⟩
    simp only [Function.update_of_ne h, Function.update_self, Function.update_idem, hf,
      Function.update_eq_self]
  · intro g
    simp [Function.update_idem]

/-- Across a disconnected pair, a fraction `1 / n` of the assignments constant along `ω` agree. -/
theorem card_compat_agree_mul [Finite V] {n : ℕ} {ω : SimpleGraph V} {i j : V}
    (hconn : ¬ ω.Reachable i j) :
    Nat.card {σ : V → Fin n // ω ≤ agreeGraph σ ∧ σ i = σ j} * n = n ^ componentCount ω := by
  have hne : ω.connectedComponentMk i ≠ ω.connectedComponentMk j := fun h =>
    hconn (SimpleGraph.ConnectedComponent.exact h)
  have hcard := card_agree_mul (β := Fin n) hne
  have hequiv : {f : ω.ConnectedComponent → Fin n //
      f (ω.connectedComponentMk i) = f (ω.connectedComponentMk j)} ≃
      {σ : V → Fin n // ω ≤ agreeGraph σ ∧ σ i = σ j} :=
    ((spinsEquiv ω).subtypeEquiv
      (q := fun x : {σ : V → Fin n // ω ≤ agreeGraph σ} => x.1 i = x.1 j) fun _ => Iff.rfl).trans
      (Equiv.subtypeSubtypeEquivSubtypeInter (fun σ : V → Fin n => ω ≤ agreeGraph σ)
        fun σ => σ i = σ j)
  rw [Nat.card_congr hequiv] at hcard
  simpa [componentCount] using hcard

theorem sum_boole_natCard {α : Type*} [Fintype α] (P : α → Prop) [DecidablePred P] :
    ∑ a, (if P a then (1 : ℝ) else 0) = Nat.card {a // P a} := by
  rw [Finset.sum_boole, Nat.card_eq_fintype_card, Fintype.card_subtype]

open Classical in
theorem sum_compat [Fintype V] (n : ℕ) (ω : Config V) :
    ∑ σ : V → Fin n, (1 : ℝ) * (if ω ≤ agreeGraph σ then 1 else 0) =
      (n : ℝ) ^ componentCount ω := by
  simp only [one_mul]
  rw [sum_boole_natCard, card_compat]
  push_cast
  rfl

open Classical in
/-- The count behind the identity: over the assignments constant along `ω`, those agreeing across
`i, j` are all of them if `i ↔ j`, and a fraction `1 / n` if not. -/
theorem sum_agree_compat [Fintype V] {n : ℕ} (hn : 1 ≤ n) (ω : Config V) (i j : V) :
    ∑ σ : V → Fin n, (if σ i = σ j then (1 : ℝ) else 0) * (if ω ≤ agreeGraph σ then 1 else 0) =
      (n : ℝ) ^ componentCount ω *
        (1 / (n : ℝ) + (1 - 1 / (n : ℝ)) * if ω ∈ Conn i j then 1 else 0) := by
  have hn₀ : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hprod : ∀ σ : V → Fin n,
      (if σ i = σ j then (1 : ℝ) else 0) * (if ω ≤ agreeGraph σ then 1 else 0) =
        if ω ≤ agreeGraph σ ∧ σ i = σ j then 1 else 0 := by
    intro σ
    by_cases h₁ : σ i = σ j <;> by_cases h₂ : ω ≤ agreeGraph σ <;> simp [h₁, h₂]
  simp only [hprod]
  rw [sum_boole_natCard]
  by_cases hconn : ω.Reachable i j
  · have hiff : ∀ σ : V → Fin n, (ω ≤ agreeGraph σ ∧ σ i = σ j) ↔ ω ≤ agreeGraph σ :=
      fun σ => ⟨fun h => h.1, fun h => ⟨h, eq_of_reachable h hconn⟩⟩
    rw [Nat.card_congr (Equiv.subtypeEquivRight hiff), card_compat,
      if_pos (show ω ∈ Conn i j from hconn)]
    push_cast
    ring
  · have hc : (Nat.card {σ : V → Fin n // ω ≤ agreeGraph σ ∧ σ i = σ j} : ℝ) * n =
        (n : ℝ) ^ componentCount ω := by
      exact_mod_cast card_compat_agree_mul hconn
    rw [if_neg (show ω ∉ Conn i j from hconn), ← hc]
    field_simp
    ring

/-! ### The identity -/

/-- The two partition functions are equal: the Potts model is the random-cluster model summed over
spins. -/
theorem pottsPartition_eq [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) (n : ℕ) :
    pottsPartition hp n = rcPartition hp n := by
  classical
  have h := sum_spin_mul_weight hp n fun _ => 1
  simp only [sum_compat] at h
  simpa only [one_mul, pottsPartition, rcPartition, rcWeight] using h

/-- **Grimmett Thm 1.16, with parameters per edge.** Spin agreement is the baseline `1 / n` plus the
rescaled connection probability. -/
theorem pottsAgree_eq [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) {n : ℕ} (hn : 1 ≤ n)
    (i j : V) :
    pottsAgree hp n i j = 1 / (n : ℝ) + (1 - 1 / (n : ℝ)) * rcProb hp n (Conn i j) := by
  classical
  have hq₀ : (0 : ℝ) < n := by exact_mod_cast hn
  have hZ := rcPartition_pos hp hq₀
  have hnum : ∑ σ : V → Fin n, (if σ i = σ j then (1 : ℝ) else 0) * pottsWeight hp σ =
      (1 / (n : ℝ)) * rcPartition hp n + (1 - 1 / (n : ℝ)) *
        ∑ ω : Config V, (if ω ∈ Conn i j then (1 : ℝ) else 0) * rcWeight hp n ω := by
    rw [sum_spin_mul_weight]
    simp only [sum_agree_compat hn]
    rw [rcPartition, Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun ω _ => ?_
    simp only [rcWeight]
    ring
  have hagree : pottsAgree hp n i j =
      (∑ σ : V → Fin n, (if σ i = σ j then (1 : ℝ) else 0) * pottsWeight hp σ) /
        rcPartition hp n := by
    rw [pottsAgree, pottsPartition_eq, Finset.sum_div]
    refine Finset.sum_congr rfl fun σ _ => ?_
    split_ifs <;> simp
  have hprob : rcProb hp n (Conn i j) =
      (∑ ω : Config V, (if ω ∈ Conn i j then (1 : ℝ) else 0) * rcWeight hp n ω) /
        rcPartition hp n := by
    rw [rcProb, Finset.sum_div]
    refine Finset.sum_congr rfl fun ω _ => ?_
    by_cases h : ω ∈ Conn i j <;> simp [h, rcMass]
  rw [hagree, hnum, hprob, add_div, mul_div_assoc, mul_div_assoc, div_self hZ.ne', mul_one]

/-- On a single edge, agreement has the value the Potts model gives by hand: `n` agreeing
assignments of weight 1 against `n (n - 1)` disagreeing ones of weight `1 - p`. This combines the
identity with `rcProb_one_edge`, so a wrong factor on either side fails here. -/
theorem pottsAgree_one_edge {p : Fin 2 → Fin 2 → ℝ} (hp : IsEdgeParameter p) {n : ℕ} (hn : 1 ≤ n) :
    pottsAgree hp n 0 1 = 1 / (p 0 1 + n * (1 - p 0 1)) := by
  have hq₀ : (0 : ℝ) < n := by exact_mod_cast hn
  have hlt := hp.lt_one 0 1
  have hge := hp.nonneg 0 1
  have hd : 0 < p 0 1 + n * (1 - p 0 1) := by nlinarith [mul_pos hq₀ (sub_pos.mpr hlt)]
  have hevent : Conn (0 : Fin 2) 1 = {ω : Config (Fin 2) | ω.Adj 0 1} := by
    ext ω
    rcases config_fin_two ω with rfl | rfl
    · simp [Conn, SimpleGraph.reachable_bot]
    · have hadj : (⊤ : SimpleGraph (Fin 2)).Adj 0 1 := by simp
      exact ⟨fun _ => hadj, fun _ => hadj.reachable⟩
  rw [pottsAgree_eq hp hn, hevent, rcProb_one_edge hp hq₀]
  field_simp
  ring

/-- `spinAgree` is the Potts two-point function at an integer `q`, and `φ` is the random-cluster
measure it is coupled to. The Potts model exists only at integer `q`. -/
def IsPottsTwoPoint [Fintype V] (spinAgree : (V → V → ℝ) → V → V → ℝ)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (q : ℝ) : Prop :=
  IsRandomCluster φ q ∧ ∃ n : ℕ, (n : ℝ) = q ∧
    ∀ (p : V → V → ℝ) (hp : IsEdgeParameter p) (i j : V), spinAgree p i j = pottsAgree hp n i j

open Classical in
/-- The Potts two-point function as a total function, zero off the parameters it covers. -/
noncomputable def pottsTwoPoint [Fintype V] (n : ℕ) (p : V → V → ℝ) (i j : V) : ℝ :=
  if hp : IsEdgeParameter p then pottsAgree hp n i j else 0

/-- So the hypothesis `IsPottsTwoPoint` can be met, at every integer `q`. -/
theorem isPottsTwoPoint_pottsTwoPoint [Fintype V] (n : ℕ) :
    IsPottsTwoPoint (pottsTwoPoint (V := V) n) (randomCluster (n : ℝ)) n :=
  ⟨isRandomCluster_randomCluster _, n, rfl, fun _ hp _ _ => dif_pos hp⟩

/-- So the symmetry hypothesis of the spin-agreement cut graph can be met too. -/
theorem pottsTwoPoint_symm [Fintype V] (n : ℕ) (p : V → V → ℝ) (i j : V) :
    pottsTwoPoint n p i j = pottsTwoPoint n p j i := by
  classical
  unfold pottsTwoPoint
  split_ifs with hp
  · exact Finset.sum_congr rfl fun σ _ => if_congr eq_comm rfl rfl
  · rfl

/-- **A2 (Grimmett Thm 1.16; with parameters per edge, (11.25)).** Spin agreement is the baseline
`1 / q` plus the rescaled connection probability. -/
theorem edwards_sokal_twopoint [Fintype V]
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ)
    (spinAgree : (V → V → ℝ) → V → V → ℝ)
    (hES : IsPottsTwoPoint spinAgree φ q)
    (p : V → V → ℝ) (hp : IsEdgeParameter p) (i j : V) :
    spinAgree p i j = 1 / q + (1 - 1 / q) * φ p (Conn i j) := by
  obtain ⟨hφ, n, rfl, hspin⟩ := hES
  rw [hspin p hp i j, hφ p hp _]
  exact pottsAgree_eq hp (by exact_mod_cast hq) i j

end Bifiltration
