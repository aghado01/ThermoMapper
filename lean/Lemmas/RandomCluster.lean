import Mathlib.Algebra.BigOperators.Field
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Combinatorics.SetFamily.FourFunctions
import Mathlib.Combinatorics.SimpleGraph.Connectivity.Connected
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.Order.Monotone.Basic
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Linarith

/-!
# The random-cluster measure and its comparison inequality

The random-cluster measure on a finite graph with a parameter per edge (Grimmett 2006, §1.4 and
§11.4), and `comparison_inequality`, his Thm 3.21: for `q ≥ 1`, raising edge parameters cannot lower
the probability of an increasing event. He states it for one `p`, and §11.4 says it stays valid per
edge. This proof covers any finite graph with a parameter per edge.

It runs through two facts. The per-edge factors satisfy the lattice inequality one edge at a time,
which is where `p₁ ≤ p₂` enters. The cluster count is supermodular, Grimmett's (3.11), which is
where `q ≥ 1` enters. `componentCount_supermodular` proves that by dimension: functions constant
along a graph's edges form a subspace whose dimension is its cluster count, and subspaces obey the
modular law. Mathlib's `holley` then needs no positivity, so `p_e = 0` needs no special case.

So connection probability is monotone in temperature outright (`connProb_mono`), and so are the
clusters of its cut graphs (`connCutGraph_reachable_mono`). Spin agreement takes one more cited
result, in `Enthymemata.Bifiltration`.

All of this is about the exact expectation. A finite-sample estimate need not be monotone.
-/

universe u

namespace Bifiltration

variable {V : Type u}

/-! ### Configurations and events -/

/-- An edge configuration is its open subgraph, ordered by inclusion. -/
abbrev Config (V : Type u) := SimpleGraph V

/-- The event that `i` and `j` share an open cluster. -/
def Conn (i j : V) : Set (Config V) := {ω | ω.Reachable i j}

theorem Conn_comm (i j : V) : Conn i j = Conn j i := by
  ext ω
  exact ⟨SimpleGraph.Reachable.symm, SimpleGraph.Reachable.symm⟩

/-- An event is increasing if adding open edges keeps it. -/
def Increasing (A : Set (Config V)) : Prop :=
  ∀ ⦃ω ω' : Config V⦄, ω ≤ ω' → ω ∈ A → ω' ∈ A

/-- **P3.** Adding open edges cannot break a path. -/
theorem connection_event_increasing (i j : V) : Increasing (Conn i j) :=
  fun _ _ hle hmem => SimpleGraph.Reachable.mono hle hmem

/-- The parameters the model covers: symmetric, and in `[0, 1)` on every pair. -/
structure IsEdgeParameter (p : V → V → ℝ) : Prop where
  nonneg : ∀ i j, 0 ≤ p i j
  lt_one : ∀ i j, p i j < 1
  symm : ∀ i j, p i j = p j i

/-! ### The cluster count is supermodular -/

/-- The number of open clusters `k(ω)`, with isolated vertices counted. -/
noncomputable def componentCount (G : SimpleGraph V) : ℕ := Nat.card G.ConnectedComponent

/-- Real functions constant along the edges of `G`. -/
def edgeConstant (G : SimpleGraph V) : Submodule ℝ (V → ℝ) where
  carrier := {f | ∀ x y, G.Adj x y → f x = f y}
  add_mem' := fun ha hb x y h => by simp [ha x y h, hb x y h]
  zero_mem' := fun _ _ _ => rfl
  smul_mem' := fun _ _ hf x y h => by simp [hf x y h]

/-- The functions constant along edges are the functions of the cluster, so there are as many
dimensions as clusters. -/
theorem finrank_edgeConstant [Finite V] (G : SimpleGraph V) :
    Module.finrank ℝ (edgeConstant G) = componentCount G := by
  letI : Fintype G.ConnectedComponent := Fintype.ofFinite _
  let L : (G.ConnectedComponent → ℝ) →ₗ[ℝ] (V → ℝ) :=
    LinearMap.funLeft ℝ ℝ G.connectedComponentMk
  have hL : Function.Injective L :=
    LinearMap.funLeft_injective_of_surjective ℝ ℝ _ Quot.mk_surjective
  have hrange : LinearMap.range L = edgeConstant G := by
    ext f
    constructor
    · rintro ⟨g, rfl⟩ x y h
      exact congrArg g (SimpleGraph.ConnectedComponent.connectedComponentMk_eq_of_adj h)
    · intro hf
      have hw : ∀ {x y : V} (w : G.Walk x y), f x = f y := by
        intro x y w
        induction w with
        | nil => rfl
        | cons h _ ih => exact (hf _ _ h).trans ih
      exact ⟨SimpleGraph.ConnectedComponent.lift f (fun _ _ w _ => hw w), rfl⟩
  rw [← hrange, LinearMap.finrank_range_of_inj hL, Module.finrank_fintype_fun_eq_card,
    componentCount, Nat.card_eq_fintype_card]

/-- **Grimmett (3.11).** `k(a) + k(b) ≤ k(a ⊓ b) + k(a ⊔ b)`. A function constant along both graphs
is constant along their join, and one constant along either is constant along their meet, so the
modular law for subspaces gives the inequality. -/
theorem componentCount_supermodular [Finite V] (a b : SimpleGraph V) :
    componentCount a + componentCount b ≤ componentCount (a ⊓ b) + componentCount (a ⊔ b) := by
  have hmod := Submodule.finrank_sup_add_finrank_inf_eq (edgeConstant a) (edgeConstant b)
  have hsup : edgeConstant a ⊔ edgeConstant b ≤ edgeConstant (a ⊓ b) :=
    sup_le (fun _ hf x y h => hf x y h.1) (fun _ hf x y h => hf x y h.2)
  have hinf : edgeConstant a ⊓ edgeConstant b ≤ edgeConstant (a ⊔ b) := fun _ hf x y h =>
    h.elim (hf.1 x y) (hf.2 x y)
  have h₁ := Submodule.finrank_mono hsup
  have h₂ := Submodule.finrank_mono hinf
  rw [finrank_edgeConstant] at h₁ h₂
  rw [finrank_edgeConstant, finrank_edgeConstant] at hmod
  omega

/-- With no edges open, every vertex is its own cluster. -/
theorem componentCount_bot : componentCount (⊥ : SimpleGraph V) = Nat.card V := by
  refine Nat.card_congr (Equiv.ofBijective (⊥ : SimpleGraph V).connectedComponentMk
    ⟨fun _ _ h => ?_, Quot.mk_surjective⟩).symm
  exact SimpleGraph.reachable_bot.mp (SimpleGraph.ConnectedComponent.exact h)

/-- With every edge open, there is one cluster. -/
theorem componentCount_top [Nonempty V] : componentCount (⊤ : SimpleGraph V) = 1 := by
  have hreach : ∀ x y : V, (⊤ : SimpleGraph V).Reachable x y := fun x y => by
    by_cases h : x = y
    · exact h ▸ SimpleGraph.Reachable.refl x
    · exact SimpleGraph.Adj.reachable h
  have hsub : Subsingleton (⊤ : SimpleGraph V).ConnectedComponent := by
    constructor
    intro c d
    induction c using SimpleGraph.ConnectedComponent.ind with
    | h x =>
      induction d using SimpleGraph.ConnectedComponent.ind with
      | h y => exact SimpleGraph.ConnectedComponent.sound (hreach x y)
  exact Nat.card_of_subsingleton (Classical.arbitrary _)

/-! ### The measure -/

/-- `p` on unordered pairs. -/
def edgeParam {p : V → V → ℝ} (hp : IsEdgeParameter p) : Sym2 V → ℝ :=
  Sym2.lift ⟨p, hp.symm⟩

theorem edgeParam_mk {p : V → V → ℝ} (hp : IsEdgeParameter p) (i j : V) :
    edgeParam hp s(i, j) = p i j :=
  Sym2.lift_mk _ i j

theorem edgeParam_nonneg {p : V → V → ℝ} (hp : IsEdgeParameter p) (e : Sym2 V) :
    0 ≤ edgeParam hp e := by
  induction e using Sym2.ind with
  | _ i j => rw [edgeParam_mk]; exact hp.nonneg i j

theorem edgeParam_lt_one {p : V → V → ℝ} (hp : IsEdgeParameter p) (e : Sym2 V) :
    edgeParam hp e < 1 := by
  induction e using Sym2.ind with
  | _ i j => rw [edgeParam_mk]; exact hp.lt_one i j

theorem edgeParam_le {p₁ p₂ : V → V → ℝ} (hp₁ : IsEdgeParameter p₁) (hp₂ : IsEdgeParameter p₂)
    (hp : ∀ i j, p₁ i j ≤ p₂ i j) (e : Sym2 V) : edgeParam hp₁ e ≤ edgeParam hp₂ e := by
  induction e using Sym2.ind with
  | _ i j => rw [edgeParam_mk, edgeParam_mk]; exact hp i j

open Classical in
/-- An edge's factor in `ω`: `p_e` if open, `1 - p_e` if closed. -/
noncomputable def edgeFactor {p : V → V → ℝ} (hp : IsEdgeParameter p) (ω : Config V) (e : Sym2 V) :
    ℝ :=
  if e ∈ ω.edgeSet then edgeParam hp e else 1 - edgeParam hp e

open Classical in
/-- The unnormalized weight `∏_e p_e^ω(e) (1 - p_e)^(1 - ω(e)) · q^k(ω)`, over the unordered pairs
of distinct vertices. An edge with `p_e = 0` is never open, so a sparser graph is the case where `p`
vanishes off it. -/
noncomputable def rcWeight [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) (q : ℝ)
    (ω : Config V) : ℝ :=
  (∏ e ∈ Finset.univ.filter (fun e : Sym2 V => ¬ e.IsDiag), edgeFactor hp ω e) *
    q ^ componentCount ω

open Classical in
/-- The partition function `Z`. -/
noncomputable def rcPartition [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) (q : ℝ) : ℝ :=
  ∑ ω : Config V, rcWeight hp q ω

/-- The probability of the configuration `ω`. -/
noncomputable def rcMass [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) (q : ℝ)
    (ω : Config V) : ℝ :=
  rcWeight hp q ω / rcPartition hp q

open Classical in
/-- The probability of the event `A`. -/
noncomputable def rcProb [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) (q : ℝ)
    (A : Set (Config V)) : ℝ :=
  ∑ ω : Config V, A.indicator (rcMass hp q) ω

/-- `φ` is the random-cluster measure with cluster weight `q`, on every parameter it covers. -/
def IsRandomCluster [Fintype V] (φ : (V → V → ℝ) → Set (Config V) → ℝ) (q : ℝ) : Prop :=
  ∀ (p : V → V → ℝ) (hp : IsEdgeParameter p) (A : Set (Config V)), φ p A = rcProb hp q A

open Classical in
/-- The random-cluster measure as a total function, zero off the parameters it covers. -/
noncomputable def randomCluster [Fintype V] (q : ℝ) (p : V → V → ℝ) (A : Set (Config V)) : ℝ :=
  if hp : IsEdgeParameter p then rcProb hp q A else 0

/-- So the hypothesis `IsRandomCluster` can be met. -/
theorem isRandomCluster_randomCluster [Fintype V] (q : ℝ) :
    IsRandomCluster (randomCluster (V := V) q) q :=
  fun _ hp _ => dif_pos hp

/-! ### It is a probability measure -/

theorem edgeFactor_nonneg {p : V → V → ℝ} (hp : IsEdgeParameter p) (ω : Config V) (e : Sym2 V) :
    0 ≤ edgeFactor hp ω e := by
  have h₀ := edgeParam_nonneg hp e
  have h₁ := edgeParam_lt_one hp e
  unfold edgeFactor
  split_ifs <;> linarith

theorem edgeFactor_prod_nonneg {p : V → V → ℝ} (hp : IsEdgeParameter p) (ω : Config V)
    (s : Finset (Sym2 V)) : 0 ≤ ∏ e ∈ s, edgeFactor hp ω e :=
  Finset.prod_nonneg fun e _ => edgeFactor_nonneg hp ω e

theorem rcWeight_nonneg [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) {q : ℝ} (hq : 0 ≤ q)
    (ω : Config V) : 0 ≤ rcWeight hp q ω :=
  mul_nonneg (edgeFactor_prod_nonneg hp ω _) (pow_nonneg hq _)

/-- The empty configuration always has positive weight, which is what `p_e < 1` is for. -/
theorem rcWeight_bot_pos [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) {q : ℝ} (hq : 0 < q) :
    0 < rcWeight hp q ⊥ := by
  refine mul_pos (Finset.prod_pos fun e _ => ?_) (pow_pos hq _)
  have h₁ := edgeParam_lt_one hp e
  have hbot : e ∉ (⊥ : SimpleGraph V).edgeSet := by simp
  simp only [edgeFactor, hbot, if_false]
  linarith

theorem rcPartition_pos [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) {q : ℝ} (hq : 0 < q) :
    0 < rcPartition hp q := by
  classical
  exact Finset.sum_pos' (fun ω _ => rcWeight_nonneg hp hq.le ω)
    ⟨⊥, Finset.mem_univ _, rcWeight_bot_pos hp hq⟩

theorem rcMass_nonneg [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) {q : ℝ} (hq : 0 < q)
    (ω : Config V) : 0 ≤ rcMass hp q ω :=
  div_nonneg (rcWeight_nonneg hp hq.le ω) (rcPartition_pos hp hq).le

theorem sum_rcMass [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) {q : ℝ} (hq : 0 < q) :
    open Classical in ∑ ω : Config V, rcMass hp q ω = 1 := by
  classical
  simp only [rcMass]
  rw [← Finset.sum_div]
  exact div_self (rcPartition_pos hp hq).ne'

theorem rcProb_univ [Fintype V] {p : V → V → ℝ} (hp : IsEdgeParameter p) {q : ℝ} (hq : 0 < q) :
    rcProb hp q Set.univ = 1 := by
  simp only [rcProb, Set.indicator_univ]
  exact sum_rcMass hp hq

/-! ### An outside check on the definition: one edge -/

theorem config_fin_two (G : SimpleGraph (Fin 2)) : G = ⊥ ∨ G = ⊤ := by
  by_cases h : G.Adj 0 1
  · right
    ext i j
    fin_cases i <;> fin_cases j <;> simp [h, h.symm]
  · left
    ext i j
    have h' : ¬ G.Adj 1 0 := fun h₁ => h h₁.symm
    fin_cases i <;> fin_cases j <;> simp [h, h']

theorem sum_config_fin_two (_ : Fintype (SimpleGraph (Fin 2))) (f : SimpleGraph (Fin 2) → ℝ) :
    ∑ G, f G = f ⊥ + f ⊤ := by
  classical
  have hne : (⊥ : SimpleGraph (Fin 2)) ≠ ⊤ := fun h => by
    have hadj : (⊤ : SimpleGraph (Fin 2)).Adj 0 1 := by simp
    rw [← h] at hadj
    exact hadj
  have huniv : (Finset.univ : Finset (SimpleGraph (Fin 2))) = {⊥, ⊤} := by
    ext G
    simpa using config_fin_two G
  rw [huniv, Finset.sum_pair hne]

theorem prod_pairs_fin_two (_ : Fintype (Sym2 (Fin 2)))
    (_ : DecidablePred fun e : Sym2 (Fin 2) => ¬ e.IsDiag) (f : Sym2 (Fin 2) → ℝ) :
    ∏ e ∈ Finset.univ.filter (fun e : Sym2 (Fin 2) => ¬ e.IsDiag), f e = f s(0, 1) := by
  classical
  have hfilter : Finset.univ.filter (fun e : Sym2 (Fin 2) => ¬ e.IsDiag) = {s(0, 1)} := by
    ext e
    induction e using Sym2.ind with
    | _ i j => fin_cases i <;> fin_cases j <;> simp [Sym2.eq_swap (a := (0 : Fin 2)) (b := 1)]
  rw [hfilter, Finset.prod_singleton]

/-- On a single edge the measure gives Grimmett's (3.3) for an edge whose endpoints are joined no
other way: open with probability `p / (p + q (1 - p))`. The formula comes from the book rather than
from the definition, so a wrong factor or cluster weight fails here. -/
theorem rcProb_one_edge {p : Fin 2 → Fin 2 → ℝ} (hp : IsEdgeParameter p) {q : ℝ} (hq : 0 < q) :
    rcProb hp q {ω | ω.Adj 0 1} = p 0 1 / (p 0 1 + q * (1 - p 0 1)) := by
  have hlt := hp.lt_one 0 1
  have hge := hp.nonneg 0 1
  have hbot : rcWeight hp q ⊥ = (1 - p 0 1) * q ^ 2 := by
    have hnot : s((0 : Fin 2), 1) ∉ (⊥ : SimpleGraph (Fin 2)).edgeSet := by simp
    simp only [rcWeight]
    rw [prod_pairs_fin_two, componentCount_bot, Nat.card_eq_fintype_card, Fintype.card_fin]
    simp only [edgeFactor, hnot, if_false, edgeParam_mk]
  have htop : rcWeight hp q ⊤ = p 0 1 * q ^ 1 := by
    have hmem : s((0 : Fin 2), 1) ∈ (⊤ : SimpleGraph (Fin 2)).edgeSet := by simp
    simp only [rcWeight]
    rw [prod_pairs_fin_two, componentCount_top]
    simp only [edgeFactor, hmem, if_true, edgeParam_mk]
  have hZ : rcPartition hp q = q * (p 0 1 + q * (1 - p 0 1)) := by
    simp only [rcPartition]
    rw [sum_config_fin_two, hbot, htop]
    ring
  have hbot_not : (⊥ : SimpleGraph (Fin 2)) ∉ {ω : Config (Fin 2) | ω.Adj 0 1} := by simp
  have htop_mem : (⊤ : SimpleGraph (Fin 2)) ∈ {ω : Config (Fin 2) | ω.Adj 0 1} := by simp
  simp only [rcProb]
  rw [sum_config_fin_two, Set.indicator_of_notMem hbot_not, Set.indicator_of_mem htop_mem, zero_add]
  simp only [rcMass]
  rw [htop, hZ, pow_one, mul_comm (p 0 1) q, mul_div_mul_left _ _ hq.ne']

/-! ### The comparison inequality -/

/-- The lattice inequality one edge at a time. Only an edge open in `a` and closed in `b` is not an
equality, and there it says `p₁ ≤ p₂`. -/
theorem edgeFactor_lattice {p₁ p₂ : V → V → ℝ} (hp₁ : IsEdgeParameter p₁)
    (hp₂ : IsEdgeParameter p₂) (hp : ∀ i j, p₁ i j ≤ p₂ i j) (a b : Config V) (e : Sym2 V) :
    edgeFactor hp₁ a e * edgeFactor hp₂ b e ≤
      edgeFactor hp₁ (a ⊓ b) e * edgeFactor hp₂ (a ⊔ b) e := by
  have hle := edgeParam_le hp₁ hp₂ hp e
  by_cases ha : e ∈ a.edgeSet <;> by_cases hb : e ∈ b.edgeSet <;>
    simp only [edgeFactor, SimpleGraph.edgeSet_inf, SimpleGraph.edgeSet_sup, Set.mem_inter_iff,
      Set.mem_union, ha, hb, and_self, and_true, and_false, or_self, or_true, or_false, if_true,
      if_false, le_refl]
  nlinarith

/-- The lattice inequality for the weights: the edges by `edgeFactor_lattice`, and the clusters by
`componentCount_supermodular`, with `1 ≤ q`. -/
theorem rcWeight_lattice [Fintype V] {p₁ p₂ : V → V → ℝ} (hp₁ : IsEdgeParameter p₁)
    (hp₂ : IsEdgeParameter p₂) (hp : ∀ i j, p₁ i j ≤ p₂ i j) {q : ℝ} (hq : 1 ≤ q)
    (a b : Config V) :
    rcWeight hp₁ q a * rcWeight hp₂ q b ≤ rcWeight hp₁ q (a ⊓ b) * rcWeight hp₂ q (a ⊔ b) := by
  have hq₀ : 0 ≤ q := zero_le_one.trans hq
  simp only [rcWeight]
  rw [mul_mul_mul_comm, mul_mul_mul_comm (∏ e ∈ _, edgeFactor hp₁ (a ⊓ b) e),
    ← Finset.prod_mul_distrib, ← Finset.prod_mul_distrib, ← pow_add, ← pow_add]
  refine mul_le_mul ?_ (pow_le_pow_right₀ hq (componentCount_supermodular a b))
    (pow_nonneg hq₀ _) (Finset.prod_nonneg fun e _ => ?_)
  · exact Finset.prod_le_prod
      (fun e _ => mul_nonneg (edgeFactor_nonneg hp₁ a e) (edgeFactor_nonneg hp₂ b e))
      (fun e _ => edgeFactor_lattice hp₁ hp₂ hp a b e)
  · exact mul_nonneg (edgeFactor_nonneg hp₁ _ e) (edgeFactor_nonneg hp₂ _ e)

/-- The comparison inequality for the measure itself: Holley's inequality on the configurations,
with the indicator of `A` as the monotone function. -/
theorem rcProb_mono [Fintype V] {p₁ p₂ : V → V → ℝ} (hp₁ : IsEdgeParameter p₁)
    (hp₂ : IsEdgeParameter p₂) (hp : ∀ i j, p₁ i j ≤ p₂ i j) {q : ℝ} (hq : 1 ≤ q)
    {A : Set (Config V)} (hA : Increasing A) : rcProb hp₁ q A ≤ rcProb hp₂ q A := by
  classical
  have hq₀ : 0 < q := zero_lt_one.trans_le hq
  have hZ₁ := rcPartition_pos hp₁ hq₀
  have hZ₂ := rcPartition_pos hp₂ hq₀
  have key : ∑ ω : Config V, A.indicator (1 : Config V → ℝ) ω * rcMass hp₁ q ω ≤
      ∑ ω : Config V, A.indicator (1 : Config V → ℝ) ω * rcMass hp₂ q ω := by
    refine holley (rcMass hp₁ q) (rcMass hp₂ q) (A.indicator 1)
      (fun ω => Set.indicator_nonneg (fun _ _ => zero_le_one) ω)
      (rcMass_nonneg hp₁ hq₀) (rcMass_nonneg hp₂ hq₀) (fun a b hab => ?_)
      ((sum_rcMass hp₁ hq₀).trans (sum_rcMass hp₂ hq₀).symm) (fun a b => ?_)
    · by_cases ha : a ∈ A
      · simp [ha, hA hab ha]
      · by_cases hb : b ∈ A <;> simp [ha, hb]
    · simp only [rcMass]
      rw [div_mul_div_comm, div_mul_div_comm]
      exact div_le_div_of_nonneg_right (rcWeight_lattice hp₁ hp₂ hp hq a b)
        (mul_pos hZ₁ hZ₂).le
  have hind : ∀ {r : Config V → ℝ} (ω : Config V),
      A.indicator r ω = A.indicator (1 : Config V → ℝ) ω * r ω := by
    intro r ω
    by_cases h : ω ∈ A <;> simp [h]
  simpa only [rcProb, hind (r := rcMass _ q)] using key

/-- **A1 (Grimmett Thm 3.21, with parameters per edge).** Raising edge parameters cannot lower the
probability of an increasing event, for `q ≥ 1` on a finite graph. -/
theorem comparison_inequality [Fintype V]
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    {p₁ p₂ : V → V → ℝ} (hp₁ : IsEdgeParameter p₁) (hp₂ : IsEdgeParameter p₂)
    (hp : ∀ i j, p₁ i j ≤ p₂ i j)
    {A : Set (Config V)} (hA : Increasing A) :
    φ p₁ A ≤ φ p₂ A := by
  rw [hφ p₁ hp₁ A, hφ p₂ hp₂ A]
  exact rcProb_mono hp₁ hp₂ hp hq hA

/-! ### Temperature -/

/-- FK edge-open probability `p_e(β) = 1 - exp (-β J_e)`. -/
noncomputable def pe (J : V → V → ℝ) (β : ℝ) (i j : V) : ℝ := 1 - Real.exp (-(β * J i j))

/-- **P1.** With non-negative couplings, `p_e` is non-decreasing in β. -/
theorem pe_mono_in_beta (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (i j : V) :
    Monotone (fun β => pe J β i j) := by
  intro β₁ β₂ hβ
  have hineq : -(β₂ * J i j) ≤ -(β₁ * J i j) := by nlinarith [hJ i j]
  have := Real.exp_le_exp.mpr hineq
  simp only [pe]; linarith

/-- **P2.** β-monotonicity in the coordinatewise order on the parameters. -/
theorem p_vector_mono (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j)
    {β₁ β₂ : ℝ} (hβ : β₁ ≤ β₂) : ∀ i j, pe J β₁ i j ≤ pe J β₂ i j :=
  fun i j => pe_mono_in_beta J hJ i j hβ

/-- Symmetric non-negative couplings at `β ≥ 0` give parameters the model covers. -/
theorem pe_isEdgeParameter (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    {β : ℝ} (hβ : 0 ≤ β) : IsEdgeParameter (pe J β) where
  nonneg i j := by
    have hle : Real.exp (-(β * J i j)) ≤ 1 :=
      Real.exp_le_one_iff.mpr (neg_nonpos.mpr (mul_nonneg hβ (hJ i j)))
    simp only [pe]; linarith
  lt_one i j := by
    have := Real.exp_pos (-(β * J i j))
    simp only [pe]; linarith
  symm i j := by
    simp only [pe]; rw [hJs i j]

/-- **P4.** Connection probability is non-decreasing in β: A1 on the event of P3. -/
theorem connProb_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) (i j : V) :
    φ (pe J β₁) (Conn i j) ≤ φ (pe J β₂) (Conn i j) :=
  comparison_inequality q hq φ hφ (pe_isEdgeParameter J hJ hJs hβ₁)
    (pe_isEdgeParameter J hJ hJs (hβ₁.trans hβ)) (p_vector_mono J hJ hβ)
    (connection_event_increasing i j)

/-! ### The cut graph of the connection probability -/

/-- The cut graph at `(β, θ)`: edges whose connection probability clears the threshold. This is what
SPC's co-membership columns estimate. -/
def connCutGraph (J : V → V → ℝ) (φ : (V → V → ℝ) → Set (Config V) → ℝ) (β θ : ℝ) :
    SimpleGraph V where
  Adj i j := i ≠ j ∧ θ ≤ φ (pe J β) (Conn i j)
  symm := ⟨fun i j h => ⟨h.1.symm, by rw [Conn_comm]; exact h.2⟩⟩
  loopless := ⟨fun _ h => h.1 rfl⟩

/-- Cooling and lowering the threshold grow the cut graph. -/
theorem connCutGraph_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) {θ₁ θ₂ : ℝ} (hθ : θ₂ ≤ θ₁) :
    connCutGraph J φ β₁ θ₁ ≤ connCutGraph J φ β₂ θ₂ :=
  fun i j hadj =>
    ⟨hadj.1, hθ.trans (hadj.2.trans (connProb_mono J hJ hJs q hq φ hφ hβ₁ hβ i j))⟩

/-- The clusters nest: two points together at `(β₁, θ₁)` are together at every colder or
lower-threshold corner, so each partition refines the next. -/
theorem connCutGraph_reachable_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) {θ₁ θ₂ : ℝ} (hθ : θ₂ ≤ θ₁) {i j : V}
    (h : (connCutGraph J φ β₁ θ₁).Reachable i j) :
    (connCutGraph J φ β₂ θ₂).Reachable i j :=
  h.mono (connCutGraph_mono J hJ hJs q hq φ hφ hβ₁ hβ hθ)

end Bifiltration
