import Mathlib.Algebra.BigOperators.Field
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Data.Finset.Card
import Mathlib.Data.Finset.Max
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Fintype.Perm
import Mathlib.Data.Prod.Lex
import Mathlib.Data.Real.Basic
import Mathlib.Order.Interval.Set.OrdConnected
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# Conformal validity (CB-1 to CB-4)

The finite statements of `Protolemmata/Conformal/conformal-lemmas.md`: prediction sets valid
whatever the model does, as `GuardedEstimator`'s means are valid whatever the emulator does.
Exchangeability is written out as the average over the orderings of a fixed data set, and an
exchangeable law on a finite space as a mixture of those averages.

Conventions, fixed in CB-1. A score is a conformity score: higher conforms better. A unit's rank
`rankLe` counts the units of the ranked set `C` scoring at most its score, itself included, and
the accept rule at level `k` accepts it when that rank exceeds `k`. `C` is every unit, or the test
unit's stratum.

* **CB-1.** `card_rankLe_le` is the counting core; `coverage_orderings` the rank guarantee over
  orderings, `coverage_orderings_eq` its exact form with no ties, `coverage_strict_le` and
  `coverage_orderings_le` the band with ties, `card_rank_eq` and `rank_uniform` the uniform rank,
  `coverage_tieBreak` random tie-breaking, `coverage_given` and `coverage_law` the guarantee
  given `w` and under any exchangeable law, and `coverage_stratum` the per-stratum form.
* **CB-2.** `split_coverage`: a score fitted on a training split covers whatever the fit is;
  `exchangeable_given_training` is the hypothesis it needs, from exchangeability of every unit;
  `cqr_ordConnected`, CQR's accept set is an interval.
* **CB-4.** `sum_jointU` and `post_addOneIn`, the add-one-in posterior; the integrated score and
  the latent-posterior score both equivariant (`integratedScore_comp`, `latentScore_comp`), the
  latter's closed form and its computation from the `n`-unit draws (`latentScore_eq`,
  `sum_jointU_addOneIn`); the mixed score covering no better than the latent one
  (`accepts_latent_of_mixed`) and covering nothing in `mixed_counterexample`.
* **CB-3**, with an error bound per unit. `accepts_or_fragile` and `acceptsInflated_of_accepts`,
  the deterministic part; `coverage_close` (bounds carried by the units) and `coverage_inflated`
  (bounds by ordering and position) over orderings; `card_fragile_le_estimated` and
  `fragileEstimated_mono`, the margin read from the estimated scores and any dominating bounds.
-/

open Finset

namespace Conformal

/-! ### CB-1: ranks and the accept rule -/

section Rank

variable {ι α : Type*} [LinearOrder α]

/-- The units of `C` scoring at most `s u`: `u`'s rank, `u` itself counted when it is in `C`. -/
def rankLe (C : Finset ι) (s : ι → α) (u : ι) : ℕ := (C.filter fun j => s j ≤ s u).card

/-- The units of `C` scoring strictly below `s u`. -/
def rankLt (C : Finset ι) (s : ι → α) (u : ι) : ℕ := (C.filter fun j => s j < s u).card

/-- **The accept rule at level `k`.** `u` is accepted when more than `k` units of `C` score at most
its score, itself included: Fong and Holmes's `π(y) > k/(n+1)`. Ties count for `u`. -/
def Accepts (k : ℕ) (C : Finset ι) (s : ι → α) (u : ι) : Prop := k < rankLe C s u

/-- The strict rule: at least `k` units of `C` score strictly below `u`. Hoff's
"`c_{n+1}` above the `k`-th order statistic", which he states without ties. -/
def AcceptsStrict (k : ℕ) (C : Finset ι) (s : ι → α) (u : ι) : Prop := k ≤ rankLt C s u

instance (k : ℕ) (C : Finset ι) (s : ι → α) (u : ι) : Decidable (Accepts k C s u) :=
  inferInstanceAs (Decidable (k < _))

instance (k : ℕ) (C : Finset ι) (s : ι → α) (u : ι) : Decidable (AcceptsStrict k C s u) :=
  inferInstanceAs (Decidable (k ≤ _))

/-- **The counting core.** Whatever the scores and their ties, at most `k` units of `C` have rank
at most `k`: take the highest-scoring of them, and every one of them counts toward its rank. -/
theorem card_rankLe_le (C : Finset ι) (s : ι → α) (k : ℕ) :
    (C.filter fun u => rankLe C s u ≤ k).card ≤ k := by
  rcases (C.filter fun u => rankLe C s u ≤ k).eq_empty_or_nonempty with h | h
  · rw [h, card_empty]
    exact Nat.zero_le k
  obtain ⟨m, hm, hmax⟩ := (C.filter fun u => rankLe C s u ≤ k).exists_max_image s h
  calc (C.filter fun u => rankLe C s u ≤ k).card ≤ rankLe C s m := by
        apply card_le_card
        intro a ha
        exact mem_filter.2 ⟨(mem_filter.1 ha).1, hmax a ha⟩
    _ ≤ k := (mem_filter.1 hm).2

/-- The units of `C` scoring strictly below `u` and those scoring at least as high partition `C`. -/
theorem rankLt_add_card_ge (C : Finset ι) (s : ι → α) (u : ι) :
    rankLt C s u + (C.filter fun j => s u ≤ s j).card = C.card := by
  have h := C.card_filter_add_card_filter_not (fun j => s j < s u)
  simpa only [rankLt, not_lt] using h

/-- The dual count: at least `k` units of `C` have fewer than `k` units strictly below them. -/
theorem le_card_rankLt_lt (C : Finset ι) (s : ι → α) {k : ℕ} (hk : k ≤ C.card) :
    k ≤ (C.filter fun u => rankLt C s u < k).card := by
  have hdual := card_rankLe_le (α := αᵒᵈ) C (fun j => OrderDual.toDual (s j)) (C.card - k)
  have hrank : ∀ u, rankLe (α := αᵒᵈ) C (fun j => OrderDual.toDual (s j)) u =
      (C.filter fun j => s u ≤ s j).card := fun u => by
    simp only [rankLe, OrderDual.toDual_le_toDual]
  simp only [hrank] at hdual
  have hfilt : (C.filter fun u => rankLt C s u < k) =
      C.filter fun u => ¬ (C.filter fun j => s u ≤ s j).card ≤ C.card - k := by
    refine filter_congr fun u _ => ?_
    have := rankLt_add_card_ge C s u
    omega
  have hsplit := C.card_filter_add_card_filter_not
    (fun u => (C.filter fun j => s u ≤ s j).card ≤ C.card - k)
  rw [hfilt]
  omega

/-- A unit of `C` counts itself, so its rank exceeds the count strictly below it. -/
theorem rankLt_lt_rankLe (C : Finset ι) (s : ι → α) {u : ι} (hu : u ∈ C) :
    rankLt C s u < rankLe C s u := by
  apply card_lt_card
  refine (ssubset_iff_of_subset fun j hj => ?_).2 ⟨u, mem_filter.2 ⟨hu, le_rfl⟩, ?_⟩
  · exact mem_filter.2 ⟨(mem_filter.1 hj).1, (mem_filter.1 hj).2.le⟩
  · exact fun h => lt_irrefl _ (mem_filter.1 h).2

/-- With no ties on `C`, the rank is one more than the count strictly below. -/
theorem rankLe_eq_rankLt_add_one {C : Finset ι} {s : ι → α} (hs : Set.InjOn s C) {u : ι}
    (hu : u ∈ C) : rankLe C s u = rankLt C s u + 1 := by
  classical
  have h : (C.filter fun j => s j ≤ s u) = insert u (C.filter fun j => s j < s u) := by
    ext j
    simp only [mem_filter, mem_insert]
    constructor
    · rintro ⟨hj, hle⟩
      rcases hle.lt_or_eq with hlt | heq
      · exact Or.inr ⟨hj, hlt⟩
      · exact Or.inl (hs hj hu heq)
    · rintro (rfl | ⟨hj, hlt⟩)
      · exact ⟨hu, le_rfl⟩
      · exact ⟨hj, hlt.le⟩
  have hnot : u ∉ C.filter fun j => s j < s u := fun h => lt_irrefl _ (mem_filter.1 h).2
  rw [rankLe, h, card_insert_of_notMem hnot, rankLt]

/-- The strict rule accepts no unit the accept rule refuses. -/
theorem accepts_of_acceptsStrict {k : ℕ} {C : Finset ι} {s : ι → α} {u : ι} (hu : u ∈ C)
    (h : AcceptsStrict k C s u) : Accepts k C s u :=
  lt_of_le_of_lt h (rankLt_lt_rankLe C s hu)

/-- With no ties the two rules agree. -/
theorem accepts_iff_acceptsStrict {k : ℕ} {C : Finset ι} {s : ι → α} (hs : Set.InjOn s C)
    {u : ι} (hu : u ∈ C) : Accepts k C s u ↔ AcceptsStrict k C s u := by
  unfold Accepts AcceptsStrict
  rw [rankLe_eq_rankLt_add_one hs hu]
  omega

/-- With no ties, exactly `k` units of `C` are refused at level `k`. -/
theorem card_rankLe_le_eq {C : Finset ι} {s : ι → α} (hs : Set.InjOn s C) {k : ℕ}
    (hk : k ≤ C.card) : (C.filter fun u => rankLe C s u ≤ k).card = k := by
  have heq : (C.filter fun u => rankLe C s u ≤ k) = C.filter fun u => rankLt C s u < k := by
    refine filter_congr fun u hu => ?_
    rw [rankLe_eq_rankLt_add_one hs hu]
    omega
  have h₁ := card_rankLe_le C s k
  have h₂ := le_card_rankLt_lt C s hk
  rw [← heq] at h₂
  omega

/-- **The uniform rank.** With no ties, each rank `1, …, |C|` is held by exactly one unit. -/
theorem card_rank_eq {C : Finset ι} {s : ι → α} (hs : Set.InjOn s C) {r : ℕ} (hr : 1 ≤ r)
    (hrC : r ≤ C.card) : (C.filter fun u => rankLe C s u = r).card = 1 := by
  classical
  have hsplit : (C.filter fun u => rankLe C s u ≤ r) =
      (C.filter fun u => rankLe C s u ≤ r - 1) ∪ C.filter fun u => rankLe C s u = r := by
    ext u
    by_cases hu : u ∈ C
    · simp only [mem_filter, mem_union, hu, true_and]
      omega
    · simp only [mem_filter, mem_union, hu, false_and, or_self]
  have hdisj : Disjoint (C.filter fun u => rankLe C s u ≤ r - 1)
      (C.filter fun u => rankLe C s u = r) := by
    rw [disjoint_filter]
    intro u _ h₁ h₂
    omega
  have h₁ := card_rankLe_le_eq hs hrC
  have h₂ := card_rankLe_le_eq hs (k := r - 1) (by omega)
  rw [hsplit, card_union_of_disjoint hdisj] at h₁
  omega

/-- **The band with ties.** The accept rule accepts at most `|C| − k` units, plus those whose ties
straddle the threshold. -/
theorem card_accepts_add_le (C : Finset ι) (s : ι → α) {k : ℕ} (hk : k ≤ C.card) :
    (C.filter fun u => Accepts k C s u).card + k ≤
      C.card + (C.filter fun u => rankLt C s u < k ∧ k < rankLe C s u).card := by
  classical
  have hsub : (C.filter fun u => Accepts k C s u) ⊆
      (C.filter fun u => AcceptsStrict k C s u) ∪
        C.filter fun u => rankLt C s u < k ∧ k < rankLe C s u := by
    intro u hu
    simp only [mem_filter, mem_union, Accepts, AcceptsStrict] at hu ⊢
    by_cases h : k ≤ rankLt C s u
    · exact Or.inl ⟨hu.1, h⟩
    · exact Or.inr ⟨hu.1, by omega, hu.2⟩
  have hstrict : (C.filter fun u => AcceptsStrict k C s u).card + k ≤ C.card := by
    have h₁ := le_card_rankLt_lt C s hk
    have h₂ := C.card_filter_add_card_filter_not (fun u => AcceptsStrict k C s u)
    have h₃ : (C.filter fun u => ¬ AcceptsStrict k C s u) = C.filter fun u => rankLt C s u < k := by
      refine filter_congr fun u _ => ?_
      simp only [AcceptsStrict, not_le]
    rw [h₃] at h₂
    omega
  have := (card_le_card hsub).trans (card_union_le _ _)
  omega

end Rank

/-! ### CB-1: orderings -/

section Orderings

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The permutations that move only units of `C`. -/
def permsOn (C : Finset ι) : Finset (Equiv.Perm ι) := univ.filter fun σ => ∀ i, i ∉ C → σ i = i

theorem one_mem_permsOn (C : Finset ι) : (1 : Equiv.Perm ι) ∈ permsOn C :=
  mem_filter.2 ⟨mem_univ _, fun _ _ => rfl⟩

theorem mem_iff_of_mem_permsOn {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (i : ι) :
    σ i ∈ C ↔ i ∈ C := by
  have hfix := (mem_filter.1 hσ).2
  constructor
  · intro h
    by_contra hi
    rw [hfix i hi] at h
    exact hi h
  · intro hi
    by_contra h
    have h' : σ (σ i) = σ i := hfix (σ i) h
    rw [σ.injective h'] at h
    exact h hi

theorem swap_mem_permsOn {C : Finset ι} {t t' : ι} (ht : t ∈ C) (ht' : t' ∈ C) :
    Equiv.swap t t' ∈ permsOn C := by
  refine mem_filter.2 ⟨mem_univ _, fun i hi => Equiv.swap_apply_of_ne_of_ne ?_ ?_⟩
  · rintro rfl
    exact hi ht
  · rintro rfl
    exact hi ht'

theorem mul_mem_permsOn {C : Finset ι} {σ τ : Equiv.Perm ι} (hσ : σ ∈ permsOn C)
    (hτ : τ ∈ permsOn C) : σ * τ ∈ permsOn C := by
  refine mem_filter.2 ⟨mem_univ _, fun i hi => ?_⟩
  rw [Equiv.Perm.mul_apply, (mem_filter.1 hτ).2 i hi, (mem_filter.1 hσ).2 i hi]

theorem card_permsOn_pos (C : Finset ι) : 0 < (permsOn C).card :=
  card_pos.2 ⟨1, one_mem_permsOn C⟩

/-- **Averaging over orderings.** Over the permutations of `C`, the test unit `t` lands on each
unit of `C` equally often: counting the orderings that put it on a unit with property `p` is
counting those units. -/
theorem card_mul_card_filter_permsOn {C : Finset ι} {t : ι} (ht : t ∈ C) (p : ι → Prop)
    [DecidablePred p] :
    C.card * ((permsOn C).filter fun σ : Equiv.Perm ι => p (σ t)).card =
      (permsOn C).card * (C.filter p).card := by
  have hshift : ∀ t' ∈ C, ((permsOn C).filter fun σ : Equiv.Perm ι => p (σ t')).card =
      ((permsOn C).filter fun σ : Equiv.Perm ι => p (σ t)).card := by
    intro t' ht'
    have hsw := swap_mem_permsOn ht ht'
    apply card_equiv (Equiv.mulRight (Equiv.swap t t'))
    intro σ
    simp only [mem_filter, Equiv.coe_mulRight, Equiv.Perm.mul_apply, Equiv.swap_apply_left]
    constructor
    · rintro ⟨hσ, hp⟩
      exact ⟨mul_mem_permsOn hσ hsw, hp⟩
    · rintro ⟨hσ, hp⟩
      refine ⟨?_, hp⟩
      have := mul_mem_permsOn hσ hsw
      rwa [mul_assoc, Equiv.swap_mul_self, mul_one] at this
  have hdouble : ∑ t' ∈ C, ((permsOn C).filter fun σ : Equiv.Perm ι => p (σ t')).card =
      ∑ _σ ∈ permsOn C, (C.filter p).card := by
    simp only [card_filter]
    rw [sum_comm]
    refine sum_congr rfl fun σ hσ => ?_
    exact sum_equiv σ (fun i => (mem_iff_of_mem_permsOn hσ i).symm) fun _ _ => rfl
  rw [sum_congr rfl hshift, sum_const, sum_const, smul_eq_mul, smul_eq_mul] at hdouble
  exact hdouble

theorem rankLe_comp {α : Type*} [LinearOrder α] {C : Finset ι} {σ : Equiv.Perm ι}
    (hσ : σ ∈ permsOn C) (s : ι → α) (u : ι) : rankLe C (s ∘ σ) u = rankLe C s (σ u) := by
  apply card_equiv σ
  intro j
  simp only [mem_filter, Function.comp_apply, mem_iff_of_mem_permsOn hσ]

theorem rankLt_comp {α : Type*} [LinearOrder α] {C : Finset ι} {σ : Equiv.Perm ι}
    (hσ : σ ∈ permsOn C) (s : ι → α) (u : ι) : rankLt C (s ∘ σ) u = rankLt C s (σ u) := by
  apply card_equiv σ
  intro j
  simp only [mem_filter, Function.comp_apply, mem_iff_of_mem_permsOn hσ]

/-- A score on data sets `z : ι → Z` is **equivariant on `C`** when reordering the units of `C`
reorders their scores with them. -/
def Equivariant {Z α : Type*} (C : Finset ι) (S : (ι → Z) → ι → α) : Prop :=
  ∀ σ ∈ permsOn C, ∀ z i, S (z ∘ σ) i = S z (σ i)

/-- The share of orderings putting the test unit on a unit with property `q` is the share of units
of `C` with it. -/
theorem ratio_permsOn {C : Finset ι} {t : ι} (ht : t ∈ C) (q : ι → Prop) [DecidablePred q] :
    (((permsOn C).filter fun σ : Equiv.Perm ι => q (σ t)).card : ℝ) / (permsOn C).card =
      ((C.filter q).card : ℝ) / C.card := by
  have hcount := card_mul_card_filter_permsOn ht q
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  rw [div_eq_div_iff hg.ne' hc.ne']
  have : ((C.card * ((permsOn C).filter fun σ : Equiv.Perm ι => q (σ t)).card : ℕ) : ℝ) =
      (((permsOn C).card * (C.filter q).card : ℕ) : ℝ) := by rw [hcount]
  push_cast at this
  linarith

/-- **Coverage over orderings is a count of units.** For an equivariant score and any property of
a unit that reordering carries along, the share of orderings of a fixed data set giving the test
unit the property is the share of units of `C` that have it. -/
theorem coverage_eq {Z α : Type*} [LinearOrder α] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → α} (hS : Equivariant C S) (z : ι → Z) (p : (ι → α) → ι → Prop)
    [∀ s, DecidablePred (p s)] (hp : ∀ σ ∈ permsOn C, ∀ s u, p (s ∘ σ) u ↔ p s (σ u)) :
    (((permsOn C).filter fun σ : Equiv.Perm ι => p (S (z ∘ σ)) t).card : ℝ) / (permsOn C).card =
      ((C.filter (p (S z))).card : ℝ) / C.card := by
  have hfilt : ((permsOn C).filter fun σ : Equiv.Perm ι => p (S (z ∘ σ)) t) =
      (permsOn C).filter fun σ : Equiv.Perm ι => p (S z) (σ t) := by
    refine filter_congr fun σ hσ => ?_
    have : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
    rw [this, hp σ hσ]
  rw [hfilt]
  exact ratio_permsOn ht (p (S z))

/-- **CB-1, the rank guarantee.** Over every ordering of a fixed data set, an equivariant score's
test unit is accepted at level `k` at least `1 − k/|C|` of the time, whatever the ties. -/
theorem coverage_orderings {Z α : Type*} [LinearOrder α] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → α} (hS : Equivariant C S) (z : ι → Z) (k : ℕ) :
    1 - (k : ℝ) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (S (z ∘ σ)) t).card : ℝ) /
        (permsOn C).card := by
  rw [coverage_eq ht hS z (fun s u => Accepts k C s u)
    (fun σ hσ s u => by simp only [Accepts, rankLe_comp hσ])]
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have hrej := card_rankLe_le C (S z) k
  have hsplit := C.card_filter_add_card_filter_not (fun u => Accepts k C (S z) u)
  have hnot : (C.filter fun u => ¬ Accepts k C (S z) u) = C.filter fun u => rankLe C (S z) u ≤ k :=
    filter_congr fun u _ => by simp only [Accepts, not_lt]
  rw [hnot] at hsplit
  have hcount : (C.card : ℝ) - k ≤ (C.filter fun u => Accepts k C (S z) u).card := by
    have : C.card ≤ (C.filter fun u => Accepts k C (S z) u).card + k := by omega
    have : (C.card : ℝ) ≤ (C.filter fun u => Accepts k C (S z) u).card + k := by exact_mod_cast this
    linarith
  rw [one_sub_div hc.ne']
  exact div_le_div_of_nonneg_right hcount hc.le

/-- The strict rule covers at most `1 − k/|C|`. -/
theorem coverage_strict_le {Z α : Type*} [LinearOrder α] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → α} (hS : Equivariant C S) (z : ι → Z) {k : ℕ} (hk : k ≤ C.card) :
    (((permsOn C).filter fun σ : Equiv.Perm ι => AcceptsStrict k C (S (z ∘ σ)) t).card : ℝ) /
      (permsOn C).card ≤
      1 - (k : ℝ) / C.card := by
  rw [coverage_eq ht hS z (fun s u => AcceptsStrict k C s u)
    (fun σ hσ s u => by simp only [AcceptsStrict, rankLt_comp hσ])]
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have h₁ := le_card_rankLt_lt C (S z) hk
  have h₂ := C.card_filter_add_card_filter_not (fun u => AcceptsStrict k C (S z) u)
  have h₃ : (C.filter fun u => ¬ AcceptsStrict k C (S z) u) =
      C.filter fun u => rankLt C (S z) u < k :=
    filter_congr fun u _ => by simp only [AcceptsStrict, not_le]
  rw [h₃] at h₂
  have hcount : ((C.filter fun u => AcceptsStrict k C (S z) u).card : ℝ) ≤ C.card - k := by
    have : (C.filter fun u => AcceptsStrict k C (S z) u).card + k ≤ C.card := by omega
    have : ((C.filter fun u => AcceptsStrict k C (S z) u).card : ℝ) + k ≤ C.card := by
      exact_mod_cast this
    linarith
  rw [one_sub_div hc.ne']
  exact div_le_div_of_nonneg_right hcount hc.le

/-- **The band with ties.** The accept rule covers at most `1 − k/|C|` plus the share of units
whose ties straddle the threshold. -/
theorem coverage_orderings_le {Z α : Type*} [LinearOrder α] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → α} (hS : Equivariant C S) (z : ι → Z) {k : ℕ} (hk : k ≤ C.card) :
    (((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (S (z ∘ σ)) t).card : ℝ) /
      (permsOn C).card ≤
      1 - (k : ℝ) / C.card +
        ((C.filter fun u => rankLt C (S z) u < k ∧ k < rankLe C (S z) u).card : ℝ) / C.card := by
  rw [coverage_eq ht hS z (fun s u => Accepts k C s u)
    (fun σ hσ s u => by simp only [Accepts, rankLe_comp hσ])]
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have h := card_accepts_add_le C (S z) hk
  have h' : ((C.filter fun u => Accepts k C (S z) u).card : ℝ) + k ≤
      C.card + (C.filter fun u => rankLt C (S z) u < k ∧ k < rankLe C (S z) u).card := by
    exact_mod_cast h
  rw [one_sub_div hc.ne', ← add_div]
  exact div_le_div_of_nonneg_right (by linarith) hc.le

/-- **CB-1 with no ties: exact.** When the scores of the data set are distinct on `C`, the test
unit is accepted on exactly `1 − k/|C|` of the orderings. -/
theorem coverage_orderings_eq {Z α : Type*} [LinearOrder α] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → α} (hS : Equivariant C S) (z : ι → Z) (hinj : Set.InjOn (S z) C) {k : ℕ}
    (hk : k ≤ C.card) :
    (((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (S (z ∘ σ)) t).card : ℝ) /
      (permsOn C).card =
      1 - (k : ℝ) / C.card := by
  rw [coverage_eq ht hS z (fun s u => Accepts k C s u)
    (fun σ hσ s u => by simp only [Accepts, rankLe_comp hσ])]
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  have h₁ := card_rankLe_le_eq hinj hk
  have h₂ := C.card_filter_add_card_filter_not (fun u => Accepts k C (S z) u)
  have h₃ : (C.filter fun u => ¬ Accepts k C (S z) u) = C.filter fun u => rankLe C (S z) u ≤ k :=
    filter_congr fun u _ => by simp only [Accepts, not_lt]
  rw [h₃, h₁] at h₂
  have : ((C.filter fun u => Accepts k C (S z) u).card : ℝ) = C.card - k := by
    have : (C.filter fun u => Accepts k C (S z) u).card = C.card - k := by omega
    rw [this, Nat.cast_sub hk]
  rw [this, one_sub_div hc.ne']

/-- **The uniform rank over orderings.** With no ties, the test unit holds each rank
`1, …, |C|` on `1/|C|` of the orderings. -/
theorem rank_uniform {Z α : Type*} [LinearOrder α] {C : Finset ι} {t : ι} (ht : t ∈ C)
    {S : (ι → Z) → ι → α} (hS : Equivariant C S) (z : ι → Z) (hinj : Set.InjOn (S z) C) {r : ℕ}
    (hr : 1 ≤ r) (hrC : r ≤ C.card) :
    (((permsOn C).filter fun σ : Equiv.Perm ι => rankLe C (S (z ∘ σ)) t = r).card : ℝ) /
      (permsOn C).card =
      1 / C.card := by
  rw [coverage_eq ht hS z (fun s u => rankLe C s u = r)
    (fun σ hσ s u => by simp only [rankLe_comp hσ]), card_rank_eq hinj hr hrC, Nat.cast_one]

/-- **Random tie-breaking.** Carry a label with each unit, distinct within `C`, and break the
score's ties by it: the score stays equivariant, has no ties, and covers exactly `1 − k/|C|`. -/
theorem coverage_tieBreak {Z L α : Type*} [LinearOrder α] [LinearOrder L] {C : Finset ι} {t : ι}
    (ht : t ∈ C) {S : (ι → Z) → ι → α} (hS : Equivariant C S) (z : ι → Z × L)
    (hlab : Set.InjOn (fun i => (z i).2) C) {k : ℕ} (hk : k ≤ C.card) :
    (((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C
        (fun i => toLex (S (fun j => ((z ∘ σ) j).1) i, ((z ∘ σ) i).2)) t).card : ℝ) /
        (permsOn C).card = 1 - (k : ℝ) / C.card := by
  let S' : (ι → Z × L) → ι → α ×ₗ L := fun w i => toLex (S (fun j => (w j).1) i, (w i).2)
  have hS' : Equivariant C S' := by
    intro σ hσ w i
    simp only [S', Function.comp_apply]
    rw [show (fun j => (w (σ j)).1) = (fun j => (w j).1) ∘ σ from rfl, hS σ hσ]
  have hinj : Set.InjOn (S' z) C := by
    intro i hi j hj h
    exact hlab hi hj (congrArg (fun x => (ofLex x).2) h)
  exact coverage_orderings_eq ht hS' z hinj hk

/-! ### CB-1 under an exchangeable law -/

/-- Reordering the units, `z ↦ z ∘ σ`, as an equivalence of data sets. -/
def reorder {Z : Type*} (σ : Equiv.Perm ι) : (ι → Z) ≃ (ι → Z) where
  toFun z := z ∘ σ
  invFun z := z ∘ σ.symm
  left_inv z := by
    funext i
    simp
  right_inv z := by
    funext i
    simp

/-- An exchangeable law is the average, over orderings, of its own reorderings: the expectation of
any `A` is the expectation of `A`'s average over the orderings of each data set. -/
theorem sum_law_average {Z : Type*} [Fintype Z] {C : Finset ι} (P : (ι → Z) → ℝ)
    (hP : ∀ σ ∈ permsOn C, ∀ z, P (z ∘ σ) = P z) (A : (ι → Z) → ℝ) :
    ∑ z, P z * A z = ∑ z, P z * ((∑ σ ∈ permsOn C, A (z ∘ σ)) / (permsOn C).card) := by
  have hone : ∀ σ ∈ permsOn C, ∑ z, P z * A (z ∘ σ) = ∑ z, P z * A z := by
    intro σ hσ
    calc ∑ z, P z * A (z ∘ σ) = ∑ z, P (z ∘ σ) * A (z ∘ σ) :=
          sum_congr rfl fun z _ => by rw [hP σ hσ]
      _ = ∑ z, P z * A z :=
          Fintype.sum_equiv (reorder σ) (fun z => P (z ∘ σ) * A (z ∘ σ)) (fun z => P z * A z)
            fun _ => rfl
  have hg : ((permsOn C).card : ℝ) ≠ 0 := by exact_mod_cast (card_permsOn_pos C).ne'
  have hsum : ∑ z, P z * ((∑ σ ∈ permsOn C, A (z ∘ σ)) / (permsOn C).card) =
      (∑ σ ∈ permsOn C, ∑ z, P z * A (z ∘ σ)) / (permsOn C).card := by
    rw [sum_comm, sum_div]
    refine sum_congr rfl fun z _ => ?_
    rw [mul_div_assoc', mul_sum, sum_div]
  rw [hsum, sum_congr rfl hone, sum_const, nsmul_eq_mul]
  field_simp

/-- **CB-1 given `w`.** For a sub-law `P` on data sets (the law given one value of whatever else
the score reads, unnormalized), with the units of `C` exchangeable under it and the score
equivariant on `C`, the test unit is accepted on at least `1 − k/|C|` of its mass. -/
theorem coverage_given {Z α : Type*} [Fintype Z] [LinearOrder α] {C : Finset ι} {t : ι}
    (ht : t ∈ C) (P : (ι → Z) → ℝ) (hP0 : ∀ z, 0 ≤ P z) (hP : ∀ σ ∈ permsOn C, ∀ z, P (z ∘ σ) = P z)
    {S : (ι → Z) → ι → α} (hS : Equivariant C S) (k : ℕ) :
    (1 - (k : ℝ) / C.card) * ∑ z, P z ≤ ∑ z, P z * if Accepts k C (S z) t then 1 else 0 := by
  rw [sum_law_average P hP, mul_sum]
  refine sum_le_sum fun z _ => ?_
  rw [mul_comm]
  refine mul_le_mul_of_nonneg_left ?_ (hP0 z)
  have h := coverage_orderings ht hS z k
  rwa [← sum_boole] at h

/-- **CB-1 under an exchangeable law.** `W` carries anything the score may also read (a training
split, a sampler's seeds), with any law; given each `w` the units of `C` are exchangeable, and the
score equivariant on `C`. Then the test unit is accepted with probability at least `1 − k/|C|`. -/
theorem coverage_law {W Z α : Type*} [Fintype W] [Fintype Z] [LinearOrder α] {C : Finset ι}
    {t : ι} (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z)
    (hP1 : ∑ w, ∑ z, P w z = 1) (hP : ∀ σ ∈ permsOn C, ∀ w z, P w (z ∘ σ) = P w z)
    (S : W → (ι → Z) → ι → α) (hS : ∀ w, Equivariant C (S w)) (k : ℕ) :
    1 - (k : ℝ) / C.card ≤ ∑ w, ∑ z, P w z * if Accepts k C (S w z) t then 1 else 0 := by
  calc 1 - (k : ℝ) / C.card = ∑ w, (1 - (k : ℝ) / C.card) * ∑ z, P w z := by
        rw [← mul_sum, hP1, mul_one]
    _ ≤ _ := sum_le_sum fun w _ =>
        coverage_given ht (P w) (hP0 w) (fun σ hσ z => hP σ hσ w z) (hS w) k

/-- **CB-1 per stratum.** Units exchangeable within strata, under every permutation that keeps each
unit in its stratum, and a score equivariant under those permutations (it may read every stratum):
ranked within the test unit's stratum, the test unit is accepted with probability at least
`1 − k/m`, `m` the stratum's size. -/
theorem coverage_stratum {W Z α κ : Type*} [Fintype W] [Fintype Z] [LinearOrder α]
    (str : ι → κ) {C : Finset ι} {t : ι} (hC : ∀ j, j ∈ C ↔ str j = str t)
    (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (hP : ∀ σ : Equiv.Perm ι, (∀ i, str (σ i) = str i) → ∀ w z, P w (z ∘ σ) = P w z)
    (S : W → (ι → Z) → ι → α)
    (hS : ∀ σ : Equiv.Perm ι, (∀ i, str (σ i) = str i) → ∀ w z i, S w (z ∘ σ) i = S w z (σ i))
    (k : ℕ) :
    1 - (k : ℝ) / C.card ≤ ∑ w, ∑ z, P w z * if Accepts k C (S w z) t then 1 else 0 := by
  have hstr : ∀ σ ∈ permsOn C, ∀ i, str (σ i) = str i := by
    intro σ hσ i
    by_cases hi : i ∈ C
    · rw [(hC _).1 ((mem_iff_of_mem_permsOn hσ i).2 hi), (hC i).1 hi]
    · rw [(mem_filter.1 hσ).2 i hi]
  exact coverage_law ((hC t).2 rfl) P hP0 hP1 (fun σ hσ => hP σ (hstr σ hσ)) S
    (fun w σ hσ => hS σ (hstr σ hσ) w) k

end Orderings

/-! ### CB-2: the split form -/

section Split

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The split form's score: a function fitted on the training split, read at each unit. -/
def splitScore {W Z α : Type*} (fit : W → Z → α) (w : W) (z : ι → Z) (i : ι) : α := fit w (z i)

/-- Any fit, read unit by unit, is equivariant on every `C`. -/
theorem splitScore_equivariant {W Z α : Type*} (fit : W → Z → α) (w : W) (C : Finset ι) :
    Equivariant C (splitScore fit w) := fun _ _ _ _ => rfl

/-- **CB-2, the split form holds for any fit.** `W` is the training split with anything else the
fit reads (a sampler's seeds, its whole unmixed chain), under any law. Given `w`, the calibration
and test units are exchangeable. Then the test unit is accepted with probability at least
`1 − k/|C|`, whatever map `fit` turns `w` into a score. -/
theorem split_coverage {W Z α : Type*} [Fintype W] [Fintype Z] [LinearOrder α] {C : Finset ι}
    {t : ι} (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z)
    (hP1 : ∑ w, ∑ z, P w z = 1) (hP : ∀ σ ∈ permsOn C, ∀ w z, P w (z ∘ σ) = P w z)
    (fit : W → Z → α) (k : ℕ) :
    1 - (k : ℝ) / C.card ≤ ∑ w, ∑ z, P w z * if Accepts k C (splitScore fit w z) t then 1 else 0 :=
  coverage_law ht P hP0 hP1 hP (splitScore fit) (fun w => splitScore_equivariant fit w C) k

omit [Fintype ι] [DecidableEq ι] in
/-- The split form's hypothesis, from exchangeability of every unit: a law on training units `T`
and the rest, invariant under every permutation of all of them, is invariant under permutations
of the rest with the training split held. -/
theorem exchangeable_given_training {T Z : Type*} (P : (T ⊕ ι → Z) → ℝ)
    (hP : ∀ τ : Equiv.Perm (T ⊕ ι), ∀ z, P (z ∘ τ) = P z) (σ : Equiv.Perm ι) (w : T → Z)
    (z : ι → Z) : P (Sum.elim w (z ∘ σ)) = P (Sum.elim w z) := by
  have h : Sum.elim w (z ∘ σ) = Sum.elim w z ∘ Equiv.sumCongr (Equiv.refl T) σ := by
    funext x
    cases x <;> rfl
  rw [h, hP]

/-- CQR's conformity score, the negative of Romano, Patterson and Candès's
`max(q_lo(x) − y, y − q_hi(x))`. -/
def cqrScore {X : Type*} (qlo qhi : X → ℝ) (p : X × ℝ) : ℝ := min (p.2 - qlo p.1) (qhi p.1 - p.2)

omit [Fintype ι] [DecidableEq ι] in
/-- Raising one unit's own score, the others held, keeps it accepted. -/
theorem accepts_mono {α : Type*} [LinearOrder α] {k : ℕ} {C : Finset ι} {s s' : ι → α} {u : ι}
    (hu : s u ≤ s' u) (hs : ∀ j, j ≠ u → s' j = s j) (h : Accepts k C s u) : Accepts k C s' u := by
  classical
  refine lt_of_lt_of_le h (card_le_card fun j hj => ?_)
  rw [mem_filter] at hj ⊢
  refine ⟨hj.1, ?_⟩
  by_cases hju : j = u
  · rw [hju]
  · rw [hs j hju]
    exact hj.2.trans hu

/-- **CQR's accept set is an interval.** With the calibration units held and the test unit at `x`,
the candidates `y` accepted at level `k` form an order-connected set. -/
theorem cqr_ordConnected {ι₀ X : Type*} [Fintype ι₀] (qlo qhi : X → ℝ)
    (cal : ι₀ → X × ℝ) (x : X) (k : ℕ) :
    Set.OrdConnected {y : ℝ |
      Accepts k univ (fun i : Option ι₀ => cqrScore qlo qhi (i.elim (x, y) cal)) none} := by
  refine ⟨fun y₁ h₁ y₃ h₃ y₂ hy => ?_⟩
  have hother : ∀ y y' : ℝ, ∀ j : Option ι₀, j ≠ none →
      cqrScore qlo qhi (j.elim (x, y') cal) = cqrScore qlo qhi (j.elim (x, y) cal) := by
    intro y y' j hj
    cases j with
    | none => exact absurd rfl hj
    | some j => rfl
  have hmid : min (cqrScore qlo qhi (x, y₁)) (cqrScore qlo qhi (x, y₃)) ≤
      cqrScore qlo qhi (x, y₂) := by
    obtain ⟨h12, h23⟩ := hy
    simp only [cqrScore]
    apply le_min
    · exact (min_le_left _ _).trans ((min_le_left _ _).trans (by linarith))
    · exact (min_le_right _ _).trans ((min_le_right _ _).trans (by linarith))
  rcases le_total (cqrScore qlo qhi (x, y₁)) (cqrScore qlo qhi (x, y₃)) with h | h
  · rw [min_eq_left h] at hmid
    exact accepts_mono hmid (hother y₁ y₂) h₁
  · rw [min_eq_right h] at hmid
    exact accepts_mono hmid (hother y₃ y₂) h₃

end Split

/-! ### CB-4: latents and the score's symmetry -/

section Latent

variable {Θ V Z : Type*} [Fintype Θ] [Fintype V]
variable (prior : Θ → ℝ) (q : Θ → V → ℝ) (f : Θ → V → Z → ℝ)

/-- A unit's density with its latent integrated, `g(z | θ) = Σ_v q(v | θ) f(z | θ, v)`. -/
noncomputable def integrated (θ : Θ) (z : Z) : ℝ := ∑ v, q θ v * f θ v z

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The unnormalized posterior of `θ` given the units `z`, each unit's latent integrated. -/
noncomputable def postU (z : ι → Z) (θ : Θ) : ℝ := prior θ * ∏ i, integrated q f θ (z i)

/-- The posterior of `θ` given the units `z`. -/
noncomputable def post (z : ι → Z) (θ : Θ) : ℝ :=
  postU prior q f z θ / ∑ θ', postU prior q f z θ'

/-- The unnormalized joint posterior of `θ` and every unit's latent. -/
noncomputable def jointU (z : ι → Z) (θ : Θ) (v : ι → V) : ℝ :=
  prior θ * ∏ i, (q θ (v i) * f θ (v i) (z i))

omit [Fintype Θ] in
/-- The latents integrate out unit by unit: the joint posterior's `θ`-marginal is the posterior
with each latent integrated, so draws of `θ` with every latent are draws of `θ`. -/
theorem sum_jointU (z : ι → Z) (θ : Θ) :
    ∑ v, jointU prior q f z θ v = postU prior q f z θ := by
  have h := Finset.prod_univ_sum (t := fun _ : ι => (univ : Finset V))
    (f := fun i v => q θ v * f θ v (z i))
  rw [Fintype.piFinset_univ] at h
  unfold jointU postU integrated
  rw [h, mul_sum]

omit [Fintype Θ] [DecidableEq ι] in
/-- **CB-4(a), add-one-in.** The candidate's unit enters the posterior as one more factor, its
integrated density. -/
theorem postU_addOneIn {ι₀ : Type*} [Fintype ι₀] (z : Option ι₀ → Z) (θ : Θ) :
    postU prior q f z θ = postU prior q f (fun j => z (some j)) θ * integrated q f θ (z none) := by
  simp only [postU, Fintype.prod_option]
  ring

omit [DecidableEq ι] in
/-- **CB-4(a), normalized.** The posterior given the `n` units and the candidate is the `n`-unit
posterior reweighted by `g(z_new | θ)`: Fong and Holmes's self-normalized add-one-in weights. -/
theorem post_addOneIn {ι₀ : Type*} [Fintype ι₀] (z : Option ι₀ → Z)
    (hN : ∑ θ, postU prior q f (fun j => z (some j)) θ ≠ 0) (θ : Θ) :
    post prior q f z θ = post prior q f (fun j => z (some j)) θ * integrated q f θ (z none) /
      ∑ θ', post prior q f (fun j => z (some j)) θ' * integrated q f θ' (z none) := by
  simp only [post, postU_addOneIn prior q f z, div_mul_eq_mul_div, ← sum_div]
  rw [div_div_div_cancel_right₀ hN]

omit [Fintype Θ] [DecidableEq ι] in
theorem postU_comp (z : ι → Z) (σ : Equiv.Perm ι) (θ : Θ) :
    postU prior q f (z ∘ σ) θ = postU prior q f z θ := by
  simp only [postU, Function.comp_apply]
  rw [Equiv.prod_comp σ (fun i => integrated q f θ (z i))]

omit [DecidableEq ι] in
theorem post_comp (z : ι → Z) (σ : Equiv.Perm ι) (θ : Θ) :
    post prior q f (z ∘ σ) θ = post prior q f z θ := by
  simp only [post, postU_comp]

/-- **CB-4(b), the integrated score**: each unit's posterior predictive density given all the
units, every latent integrated. -/
noncomputable def integratedScore (z : ι → Z) (i : ι) : ℝ :=
  ∑ θ, post prior q f z θ * integrated q f θ (z i)

omit [DecidableEq ι] in
theorem integratedScore_comp (z : ι → Z) (σ : Equiv.Perm ι) (i : ι) :
    integratedScore prior q f (z ∘ σ) i = integratedScore prior q f z (σ i) := by
  simp only [integratedScore, post_comp, Function.comp_apply]

/-- **CB-4(b).** The integrated score is equivariant on every `C`. -/
theorem integratedScore_equivariant (C : Finset ι) : Equivariant C (integratedScore prior q f) :=
  fun σ _ z i => integratedScore_comp prior q f z σ i

/-- **CB-4(c), the latent-posterior score**: the expectation of `f(z_i | θ, v_i)` under the joint
posterior of `θ` and every unit's latent, given all the units. -/
noncomputable def latentScore (z : ι → Z) (i : ι) : ℝ :=
  (∑ θ, ∑ v, jointU prior q f z θ v * f θ (v i) (z i)) / ∑ θ, postU prior q f z θ

omit [Fintype Θ] [Fintype V] [DecidableEq ι] in
theorem jointU_comp (z : ι → Z) (σ : Equiv.Perm ι) (θ : Θ) (v : ι → V) :
    jointU prior q f (z ∘ σ) θ (v ∘ σ) = jointU prior q f z θ v := by
  simp only [jointU, Function.comp_apply]
  rw [Equiv.prod_comp σ (fun i => q θ (v i) * f θ (v i) (z i))]

theorem latentScore_comp (z : ι → Z) (σ : Equiv.Perm ι) (i : ι) :
    latentScore prior q f (z ∘ σ) i = latentScore prior q f z (σ i) := by
  unfold latentScore
  have hden : ∑ θ, postU prior q f (z ∘ σ) θ = ∑ θ, postU prior q f z θ :=
    sum_congr rfl fun θ _ => postU_comp prior q f z σ θ
  have hnum : ∀ θ, ∑ v, jointU prior q f (z ∘ σ) θ v * f θ (v i) ((z ∘ σ) i) =
      ∑ v, jointU prior q f z θ v * f θ (v (σ i)) (z (σ i)) := fun θ =>
    (Fintype.sum_equiv (reorder σ) (fun v => jointU prior q f z θ v * f θ (v (σ i)) (z (σ i)))
      (fun v => jointU prior q f (z ∘ σ) θ v * f θ (v i) ((z ∘ σ) i))
      (fun v => by rw [show reorder σ v = v ∘ σ from rfl, jointU_comp]; rfl)).symm
  rw [hden, sum_congr rfl fun θ _ => hnum θ]

/-- **CB-4(c).** The latent-posterior score is equivariant on every `C`. -/
theorem latentScore_equivariant (C : Finset ι) : Equivariant C (latentScore prior q f) :=
  fun σ _ z i => latentScore_comp prior q f z σ i

/-- Every latent but unit `i`'s integrates out. -/
theorem sum_prod_mul_coord (a : ι → V → ℝ) (H : V → ℝ) (i : ι) :
    ∑ v : ι → V, (∏ j, a j (v j)) * H (v i) =
      (∑ x, a i x * H x) * ∏ j ∈ univ.erase i, ∑ x, a j x := by
  have h := Finset.prod_univ_sum (t := fun _ : ι => (univ : Finset V))
    (f := fun j x => a j x * if j = i then H x else 1)
  rw [Fintype.piFinset_univ] at h
  have hl : ∏ j, ∑ x, a j x * (if j = i then H x else 1) =
      (∑ x, a i x * H x) * ∏ j ∈ univ.erase i, ∑ x, a j x := by
    rw [← mul_prod_erase univ _ (mem_univ i)]
    simp only [if_true]
    congr 1
    refine prod_congr rfl fun j hj => ?_
    simp only [if_neg (ne_of_mem_erase hj), mul_one]
  rw [← hl, h]
  refine sum_congr rfl fun v _ => ?_
  rw [prod_mul_distrib]
  simp

/-- **CB-4(c), the closed form.** Given all the units, unit `i`'s latent-posterior score is
`Σ_θ p(θ | z) · Σ_v q f² / g`: the latent's posterior given `θ` and `z_i`, `q f / g`, averages
`f`. -/
theorem latentScore_eq (z : ι → Z) (i : ι) (hg : ∀ θ, integrated q f θ (z i) ≠ 0) :
    latentScore prior q f z i = ∑ θ, post prior q f z θ *
      ((∑ v, q θ v * f θ v (z i) ^ 2) / integrated q f θ (z i)) := by
  have hθ : ∀ θ, ∑ v, jointU prior q f z θ v * f θ (v i) (z i) = postU prior q f z θ *
      ((∑ v, q θ v * f θ v (z i) ^ 2) / integrated q f θ (z i)) := by
    intro θ
    have h := sum_prod_mul_coord (fun j x => q θ x * f θ x (z j)) (fun x => f θ x (z i)) i
    have hp : postU prior q f z θ = prior θ *
        (integrated q f θ (z i) * ∏ j ∈ univ.erase i, integrated q f θ (z j)) := by
      rw [postU, mul_prod_erase univ (fun j => integrated q f θ (z j)) (mem_univ i)]
    have hl : ∑ v, jointU prior q f z θ v * f θ (v i) (z i) =
        prior θ * ∑ v : ι → V, (∏ j, q θ (v j) * f θ (v j) (z j)) * f θ (v i) (z i) := by
      rw [mul_sum]
      exact sum_congr rfl fun v _ => by rw [jointU, mul_assoc]
    rw [hl, h, hp]
    have hsq : ∑ x, q θ x * f θ x (z i) * f θ x (z i) = ∑ v, q θ v * f θ v (z i) ^ 2 :=
      sum_congr rfl fun x _ => by ring
    rw [hsq]
    have hgθ := hg θ
    simp only [integrated] at hgθ ⊢
    field_simp
  rw [latentScore, sum_congr rfl fun θ _ => hθ θ, sum_div]
  refine sum_congr rfl fun θ _ => ?_
  rw [post]
  ring

omit [Fintype Θ] in
/-- **CB-4(c) from the `n`-unit draws.** The joint posterior given the candidate too is the
`n`-unit joint posterior with one more factor, `q(v' | θ) f(z_new | θ, v')`, which is the
add-one-in weight `g(z_new | θ)` times the candidate latent's conditional `p(v' | θ, z_new)`.
So the score is computed from the retained draws, reweighted, with the candidate's latent drawn
from its conditional. -/
theorem sum_jointU_addOneIn {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z)
    (θ : Θ) (F : (Option ι₀ → V) → ℝ) :
    ∑ v, jointU prior q f z θ v * F v =
      ∑ v₀ : ι₀ → V, jointU prior q f (fun j => z (some j)) θ v₀ *
        ∑ v', q θ v' * f θ v' (z none) * F (fun o => o.elim v' v₀) := by
  rw [← Fintype.sum_equiv Equiv.piOptionEquivProd.symm
    (fun p => jointU prior q f z θ (Equiv.piOptionEquivProd.symm p) *
      F (Equiv.piOptionEquivProd.symm p)) _ fun _ => rfl]
  rw [Fintype.sum_prod_type, sum_comm]
  refine sum_congr rfl fun v₀ _ => ?_
  rw [mul_sum]
  refine sum_congr rfl fun v' _ => ?_
  have he : (Equiv.piOptionEquivProd.symm (v', v₀) : Option ι₀ → V) = fun o => o.elim v' v₀ := by
    funext o
    cases o <;> rfl
  rw [he]
  simp only [jointU, Fintype.prod_option, Option.elim_none, Option.elim_some]
  ring

/-- **CB-4(d), the mixed score**: data units by their latents, the candidate by the integrated
density. -/
noncomputable def mixedScore {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z)
    (o : Option ι₀) : ℝ :=
  o.elim (integratedScore prior q f z none) fun j => latentScore prior q f z (some j)

omit [Fintype Θ] [DecidableEq ι] in
theorem postU_nonneg (hprior : ∀ θ, 0 ≤ prior θ) (hq : ∀ θ v, 0 ≤ q θ v)
    (hf : ∀ θ v y, 0 ≤ f θ v y) (z : ι → Z) (θ : Θ) : 0 ≤ postU prior q f z θ :=
  mul_nonneg (hprior θ) (prod_nonneg fun i _ =>
    sum_nonneg fun v _ => mul_nonneg (hq θ v) (hf θ v (z i)))

/-- **The latent score sits above the integrated one**, unit by unit: `Σ_v q f² / g − g` is the
latent's conditional variance of `f` over `g`, which vanishes only where `f` does not vary with
the latent. -/
theorem integratedScore_le_latentScore (hprior : ∀ θ, 0 ≤ prior θ) (hq : ∀ θ v, 0 ≤ q θ v)
    (hf : ∀ θ v y, 0 ≤ f θ v y) (hq1 : ∀ θ, ∑ v, q θ v = 1) (z : ι → Z) (i : ι)
    (hg : ∀ θ, 0 < integrated q f θ (z i)) :
    integratedScore prior q f z i ≤ latentScore prior q f z i := by
  rw [latentScore_eq prior q f z i fun θ => (hg θ).ne', integratedScore]
  refine sum_le_sum fun θ _ => mul_le_mul_of_nonneg_left ?_
    (div_nonneg (postU_nonneg prior q f hprior hq hf z θ)
      (sum_nonneg fun θ' _ => postU_nonneg prior q f hprior hq hf z θ'))
  rw [le_div_iff₀ (hg θ)]
  have hvar : 0 ≤ ∑ v, q θ v * (f θ v (z i) - integrated q f θ (z i)) ^ 2 :=
    sum_nonneg fun v _ => mul_nonneg (hq θ v) (sq_nonneg _)
  have hterm : ∀ v, q θ v * (f θ v (z i) - integrated q f θ (z i)) ^ 2 =
      q θ v * f θ v (z i) ^ 2 - 2 * integrated q f θ (z i) * (q θ v * f θ v (z i)) +
        integrated q f θ (z i) ^ 2 * q θ v := fun v => by ring
  have hsum : ∑ v, q θ v * (f θ v (z i) - integrated q f θ (z i)) ^ 2 =
      ∑ v, q θ v * f θ v (z i) ^ 2 - 2 * integrated q f θ (z i) * integrated q f θ (z i) +
        integrated q f θ (z i) ^ 2 * 1 := by
    rw [sum_congr rfl fun v _ => hterm v, sum_add_distrib, sum_sub_distrib, ← mul_sum, ← mul_sum,
      hq1]
    rfl
  rw [hsum] at hvar
  nlinarith [hvar]

/-- **CB-4(d): the mixed score accepts no candidate the latent score refuses.** The data units
score alike under both, and the candidate scores lower under the mixed one. -/
theorem accepts_latent_of_mixed {ι₀ : Type*} [Fintype ι₀] [DecidableEq ι₀] (z : Option ι₀ → Z)
    (hle : integratedScore prior q f z none ≤ latentScore prior q f z none) {k : ℕ}
    (h : Accepts k univ (mixedScore prior q f z) none) :
    Accepts k univ (latentScore prior q f z) none := by
  refine lt_of_lt_of_le h (card_le_card fun j hj => ?_)
  rw [mem_filter] at hj ⊢
  refine ⟨hj.1, ?_⟩
  cases j with
  | none => exact le_rfl
  | some j => exact le_trans hj.2 hle

/-- The counterexample's latent law: a fair bit. -/
noncomputable def exQ : Unit → Bool → ℝ := fun _ _ => 1 / 2

/-- The counterexample's density: `y` is the latent bit with probability `3/4`. -/
noncomputable def exF : Unit → Bool → Bool → ℝ := fun _ v y => if v = y then 3 / 4 else 1 / 4

/-- **CB-4(d), the counterexample.** One `θ`, a fair latent bit, `y` the bit with probability
`3/4`, two units. The data unit scores `5/8` by its latent and the candidate `1/2` integrated,
whatever the data, so the candidate is never accepted at `k = 1`: coverage `0`, against CB-1's
`1/2` for either symmetric score. -/
theorem mixed_counterexample (z : Option Unit → Bool) :
    ¬ Accepts 1 univ (mixedScore (fun _ => (1 : ℝ)) exQ exF z) none := by
  have hg : ∀ y, integrated exQ exF () y = 1 / 2 := by
    intro y
    cases y <;> simp [integrated, exQ, exF] <;> norm_num
  have hpost : ∀ θ, post (fun _ => (1 : ℝ)) exQ exF z θ = 1 := by
    intro θ
    have hU : postU (fun _ => (1 : ℝ)) exQ exF z θ = 1 / 4 := by
      simp [postU, hg]
      norm_num
    simp only [post, Fintype.sum_unique, hU]
    norm_num
  have hint : mixedScore (fun _ => (1 : ℝ)) exQ exF z none = 1 / 2 := by
    simp [mixedScore, integratedScore, hpost, hg]
  have hlat : ∀ u, mixedScore (fun _ => (1 : ℝ)) exQ exF z (some u) = 5 / 8 := by
    intro u
    change latentScore (fun _ => (1 : ℝ)) exQ exF z (some u) = 5 / 8
    rw [latentScore_eq _ _ _ z _ fun θ => by rw [hg]; norm_num]
    simp only [Fintype.sum_unique, hpost, hg, exQ, exF, Fintype.sum_bool]
    cases z (some u) <;> norm_num
  have hfilt : (univ.filter fun j => mixedScore (fun _ => (1 : ℝ)) exQ exF z j ≤
      mixedScore (fun _ => (1 : ℝ)) exQ exF z none) = {none} := by
    ext j
    cases j with
    | none => simp
    | some u =>
      simp only [mem_filter, mem_univ, true_and, mem_singleton, reduceCtorEq, iff_false, not_le,
        hint, hlat]
      norm_num
  simp only [Accepts, rankLe, hfilt, card_singleton, lt_irrefl, not_false_eq_true]

end Latent

/-! ### CB-3: what an imperfect score costs -/

section Margin

variable {ι : Type*} [DecidableEq ι]

/-- The units of `C` at least `δ j` below `s u`, `u` itself always counted. -/
noncomputable def rankShift (C : Finset ι) (s : ι → ℝ) (δ : ι → ℝ) (u : ι) : ℕ :=
  (C.filter fun j => j = u ∨ s j + δ j ≤ s u).card

/-- **The fragile units** under per-unit error bounds `ε`: accepted by the exact score, but with at
most `k` units of `C` at least `ε j + ε u` below it, itself counted, so estimates within their
bounds can refuse it. With every `ε j = ε` the shift is `2ε`. -/
def Fragile (k : ℕ) (C : Finset ι) (s : ι → ℝ) (ε : ι → ℝ) (u : ι) : Prop :=
  k < rankLe C s u ∧ rankShift C s (fun j => ε j + ε u) u ≤ k

/-- **The inflated rule** (Ndiaye's stable set): accept `u` when more than `k` units of `C` have an
estimate, lowered by its bound, at most `u`'s, raised by its own. With every `ε j = ε`: at most
`2ε` above `u`'s. -/
def AcceptsInflated (k : ℕ) (C : Finset ι) (s : ι → ℝ) (ε : ι → ℝ) (u : ι) : Prop :=
  k < (C.filter fun j => s j - ε j ≤ s u + ε u).card

/-- The fragile band read from estimated scores and their bounds alone: accepted by the inflated
rule, with at most `k` units at least `2ε j + 2ε u` below. With every `ε j = ε`: `4ε`. -/
def FragileEstimated (k : ℕ) (C : Finset ι) (s : ι → ℝ) (ε : ι → ℝ) (u : ι) : Prop :=
  AcceptsInflated k C s ε u ∧ rankShift C s (fun j => 2 * ε j + 2 * ε u) u ≤ k

noncomputable instance (k : ℕ) (C : Finset ι) (s : ι → ℝ) (ε : ι → ℝ) (u : ι) :
    Decidable (Fragile k C s ε u) :=
  inferInstanceAs (Decidable (_ ∧ _))

noncomputable instance (k : ℕ) (C : Finset ι) (s : ι → ℝ) (ε : ι → ℝ) (u : ι) :
    Decidable (AcceptsInflated k C s ε u) :=
  inferInstanceAs (Decidable (_ < _))

noncomputable instance (k : ℕ) (C : Finset ι) (s : ι → ℝ) (ε : ι → ℝ) (u : ι) :
    Decidable (FragileEstimated k C s ε u) :=
  inferInstanceAs (Decidable (_ ∧ _))

/-- Estimates within their bounds on `C` give each unit of `C` at least the exact score's shifted
rank. -/
theorem rankShift_le_rankLe {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) :
    rankShift C s (fun j => ε j + ε u) u ≤ rankLe C ŝ u := by
  refine card_le_card fun j hj => ?_
  rw [mem_filter] at hj ⊢
  refine ⟨hj.1, ?_⟩
  rcases hj.2 with rfl | h
  · exact le_rfl
  · have hj' := abs_le.1 (hclose j hj.1)
    have hu' := abs_le.1 (hclose u hu)
    linarith

/-- **CB-3, the deterministic part.** Estimates within their bounds refuse a unit the exact score
accepts only when that unit is fragile. -/
theorem accepts_or_fragile {k : ℕ} {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) (h : Accepts k C s u) :
    Accepts k C ŝ u ∨ Fragile k C s ε u := by
  by_cases hf : rankShift C s (fun j => ε j + ε u) u ≤ k
  · exact Or.inr ⟨h, hf⟩
  · exact Or.inl (lt_of_lt_of_le (not_le.1 hf) (rankShift_le_rankLe hclose hu))

omit [DecidableEq ι] in
/-- The inflated rule keeps every unit the exact score accepts. -/
theorem acceptsInflated_of_accepts {k : ℕ} {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) (h : Accepts k C s u) :
    AcceptsInflated k C ŝ ε u := by
  refine lt_of_lt_of_le h (card_le_card fun j hj => ?_)
  rw [mem_filter] at hj ⊢
  have hj' := abs_le.1 (hclose j hj.1)
  have hu' := abs_le.1 (hclose u hu)
  exact ⟨hj.1, by linarith [hj.2]⟩

omit [DecidableEq ι] in
/-- The other side: estimates accept a unit the exact score refuses only when that unit sits within
`ε j + ε u` of acceptance (Ndiaye's deflated set). With `accepts_or_fragile`, a decision changes
only inside the band. -/
theorem accepts_shift_of_accepts_close {k : ℕ} {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) (h : Accepts k C ŝ u) :
    k < (C.filter fun j => s j ≤ s u + (ε j + ε u)).card := by
  refine lt_of_lt_of_le h (card_le_card fun j hj => ?_)
  rw [mem_filter] at hj ⊢
  have hj' := abs_le.1 (hclose j hj.1)
  have hu' := abs_le.1 (hclose u hu)
  exact ⟨hj.1, by linarith [hj.2]⟩

/-- **The engine's reading.** A fragile unit of the exact score sits in the band read from the
estimated scores and their bounds alone. -/
theorem fragileEstimated_of_fragile {k : ℕ} {C : Finset ι} {s ŝ ε : ι → ℝ}
    (hclose : ∀ j ∈ C, |ŝ j - s j| ≤ ε j) {u : ι} (hu : u ∈ C) (h : Fragile k C s ε u) :
    FragileEstimated k C ŝ ε u := by
  refine ⟨acceptsInflated_of_accepts hclose hu h.1, le_trans (card_le_card fun j hj => ?_) h.2⟩
  rw [mem_filter] at hj ⊢
  refine ⟨hj.1, ?_⟩
  rcases hj.2 with rfl | h'
  · exact Or.inl rfl
  · have hj' := abs_le.1 (hclose j hj.1)
    have hu' := abs_le.1 (hclose u hu)
    exact Or.inr (by linarith)

/-- The estimated band only widens with the bounds, so the engine may read it with any bounds that
dominate the true ones. -/
theorem fragileEstimated_mono {k : ℕ} {C : Finset ι} {s ε ε' : ι → ℝ} (hε : ∀ j, ε j ≤ ε' j)
    {u : ι} (h : FragileEstimated k C s ε u) : FragileEstimated k C s ε' u := by
  refine ⟨lt_of_lt_of_le h.1 (card_le_card fun j hj => ?_),
    le_trans (card_le_card fun j hj => ?_) h.2⟩
  · rw [mem_filter] at hj ⊢
    exact ⟨hj.1, by linarith [hj.2, hε j, hε u]⟩
  · rw [mem_filter] at hj ⊢
    refine ⟨hj.1, ?_⟩
    rcases hj.2 with h' | h'
    · exact Or.inl h'
    · exact Or.inr (by linarith [hε j, hε u])

variable [Fintype ι]

theorem fragile_comp {k : ℕ} {C : Finset ι} {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C)
    (s ε : ι → ℝ) (u : ι) : Fragile k C (s ∘ σ) (ε ∘ σ) u ↔ Fragile k C s ε (σ u) := by
  have h : rankShift C (s ∘ σ) (fun j => (ε ∘ σ) j + (ε ∘ σ) u) u =
      rankShift C s (fun j => ε j + ε (σ u)) (σ u) := by
    apply card_equiv σ
    intro j
    simp only [mem_filter, Function.comp_apply, mem_iff_of_mem_permsOn hσ, σ.injective.eq_iff]
  simp only [Fragile, rankLe_comp hσ, h]

/-- **CB-3 over orderings.** The exact score `S` is equivariant; each unit's error bound travels
with it (`ε` indexed by the units of the fixed data set); the estimate `ŝ σ` for the ordering `σ`
is anything within those bounds on `C` (Monte Carlo error, an asymmetric computation). The estimate
covers at least `1 − (k + |M|)/|C|`, `M` the fragile units of the exact score. -/
theorem coverage_close {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C) {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant C S) (z : ι → Z) (ŝ : Equiv.Perm ι → ι → ℝ) (ε : ι → ℝ)
    (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j - S (z ∘ σ) j| ≤ ε (σ j)) (k : ℕ) :
    1 - ((k : ℝ) + (C.filter fun u => Fragile k C (S z) ε u).card) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (ŝ σ) t).card : ℝ) /
        (permsOn C).card := by
  have hfun : ∀ σ ∈ permsOn C, S (z ∘ σ) = S z ∘ σ := fun σ hσ => funext fun i => hS σ hσ z i
  have hsub : ((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (S (z ∘ σ)) t) ⊆
      ((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (ŝ σ) t) ∪
        (permsOn C).filter fun σ : Equiv.Perm ι => Fragile k C (S z) ε (σ t) := by
    intro σ hσ
    rw [mem_filter] at hσ
    rw [mem_union, mem_filter, mem_filter]
    rcases accepts_or_fragile (ε := ε ∘ σ) (hclose σ hσ.1) ht hσ.2 with h | h
    · exact Or.inl ⟨hσ.1, h⟩
    · rw [hfun σ hσ.1, fragile_comp hσ.1] at h
      exact Or.inr ⟨hσ.1, h⟩
  have hcard := (card_le_card hsub).trans (card_union_le _ _)
  have hacc := coverage_orderings ht hS z k
  have hfrag := ratio_permsOn ht (fun u => Fragile k C (S z) ε u)
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  have hcard' : (((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (S (z ∘ σ)) t).card : ℝ)
      / (permsOn C).card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (ŝ σ) t).card : ℝ) /
        (permsOn C).card +
      (((permsOn C).filter fun σ : Equiv.Perm ι => Fragile k C (S z) ε (σ t)).card : ℝ) /
        (permsOn C).card := by
    rw [← add_div]
    exact div_le_div_of_nonneg_right (by exact_mod_cast hcard) hg.le
  rw [hfrag] at hcard'
  rw [add_div, ← sub_sub]
  linarith

/-- **The inflated rule over orderings** keeps CB-1's guarantee, `1 − k/|C|`, at the price of a
wider set. The bounds may differ by ordering and by position, as Monte Carlo errors do (the
candidate's against a data unit's). -/
theorem coverage_inflated {Z : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C) {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant C S) (z : ι → Z) (ŝ ε : Equiv.Perm ι → ι → ℝ)
    (hclose : ∀ σ ∈ permsOn C, ∀ j ∈ C, |ŝ σ j - S (z ∘ σ) j| ≤ ε σ j) (k : ℕ) :
    1 - (k : ℝ) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι => AcceptsInflated k C (ŝ σ) (ε σ) t).card : ℝ) /
        (permsOn C).card := by
  have hsub : ((permsOn C).filter fun σ : Equiv.Perm ι => Accepts k C (S (z ∘ σ)) t) ⊆
      (permsOn C).filter fun σ : Equiv.Perm ι => AcceptsInflated k C (ŝ σ) (ε σ) t := by
    intro σ hσ
    rw [mem_filter] at hσ ⊢
    exact ⟨hσ.1, acceptsInflated_of_accepts (hclose σ hσ.1) ht hσ.2⟩
  have hg : (0 : ℝ) < (permsOn C).card := by exact_mod_cast card_permsOn_pos C
  exact (coverage_orderings ht hS z k).trans
    (div_le_div_of_nonneg_right (by exact_mod_cast card_le_card hsub) hg.le)

/-- **The margin read from the estimate.** For every ordering, the exact score's fragile count is
at most the count of the band read from that ordering's estimates and the units' bounds, so the
engine can report `1 − (k + m̂)/|C|` with `m̂` read from what it holds (with any dominating bounds,
by `fragileEstimated_mono`). -/
theorem card_fragile_le_estimated {Z : Type*} {C : Finset ι} {S : (ι → Z) → ι → ℝ}
    (hS : Equivariant C S) (z : ι → Z) {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (ŝ ε : ι → ℝ)
    (hclose : ∀ j ∈ C, |ŝ j - S (z ∘ σ) j| ≤ ε (σ j)) (k : ℕ) :
    (C.filter fun u => Fragile k C (S z) ε u).card ≤
      (C.filter fun u => FragileEstimated k C ŝ (ε ∘ σ) u).card := by
  have heq : (C.filter fun u => Fragile k C (S z) ε u).card =
      (C.filter fun u => Fragile k C (S (z ∘ σ)) (ε ∘ σ) u).card := by
    symm
    apply card_equiv σ
    intro u
    have : S (z ∘ σ) = S z ∘ σ := funext fun i => hS σ hσ z i
    simp only [mem_filter, this, fragile_comp hσ, mem_iff_of_mem_permsOn hσ]
  rw [heq]
  refine card_le_card fun u hu => ?_
  rw [mem_filter] at hu ⊢
  exact ⟨hu.1, fragileEstimated_of_fragile hclose hu.1 hu.2⟩

end Margin

end Conformal
