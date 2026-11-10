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
import Mathlib.GroupTheory.Perm.Fin
import Mathlib.MeasureTheory.Constructions.Pi
import Mathlib.MeasureTheory.Integral.Lebesgue.Markov
import Mathlib.Order.Interval.Set.OrdConnected
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# Conformal validity (CB-1 to CB-4, E35's M-1 and M-2, and slice 9's cross-fitted forms)

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
* **M-1** (E35, the smoothed p-values). `seqRank_uniform`, the finite core: the randomized ranks of
  the smoothed p-values are independent and uniform, `sum_tieRankProb` the tie class covering a
  rank; `seqRank_law` under any exchangeable law.
* **M-2** (E35, the betting product). `lintegral_capital_snoc` and `lintegral_capital`, the link
  (a martingale with expectation 1); `lintegral_capital_le` its supermartingale form;
  `lintegral_stoppedCapital`, the e-value at a stopping time; `ville_horizon` and `ville`, Ville's
  inequality.
* **Slice 9, the jackknife+** (Barber, Candès, Ramdas and Tibshirani). `card_refused_pairwise_le`,
  the counting core: at most `2k − 1` units are each stranger than `|C| − k` others or more. The
  cross-fitted rule covers at least `1 − (2k − 1)/|C|` over orderings
  (`coverage_pairwise_orderings`) and under a law exchangeable up to a relabeling of the fits'
  randomness (`coverage_pairwise`, `crossFit_coverage`). The jackknife+ interval refuses only a
  strange test unit (`not_acceptsPairwise_of_jackRefuses`), so `jackknifePlus_orderings`,
  `jackknifePlus_coverage` and `jackknifePlus_two_alpha`, its `1 − 2α`. The same bound holds for a
  fit that reads its units' order (`OrderReading`): `coverage_orderReading` and
  `jackknifePlus_orderReading`, by moving any unit to the front with the others kept in order
  (`acceptsPairwise_cycleRange`).
* **Slice 9, CV+ with `K` folds** (Barber et al., Theorem 4). By merging the folds' ranks (Vovk):
  `merged_le`, ranks super-uniform on `1, …, M` summing to at most `L` with probability at most
  `2L/(KM)`, through `layer_cake`, `sum_antitone_le` and `sum_max_sub_le`; `crossFold_refused_le`,
  `sum_rankLe_lt_of_cvRefuses`, `cvPlus_merged` and `cvPlus_merged_alpha`, with no symmetric fit;
  `cvPlus_oneFold`, one fold alone, at most `ℓ/(m + 1)`, above `1 − 2α` with two folds.
  By a tournament across folds: `card_refused_across_le`, at most `2ℓ + m − 2` strange units;
  `ratio_prodPerm`, averaging over folds and slots; `coverage_across`,
  `not_acceptsAcross_of_cvRefuses` and `cvPlus_tournament`.
* **Slice 9, asymmetric CQR** (Romano, Patterson and Candès, Theorem 2). `cqrAsym_coverage`, each
  tail missing at most `k/|C|` and the set covering at least `1 − (k_lo + k_hi)/|C|`;
  `cqrAsym_coverage_half` at `α/2` per tail; `cqrLo_isUpperSet`, `cqrHi_isLowerSet` and
  `cqrAsym_ordConnected`, each tail refusing on its own side of an interval.
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

/-! ### E35, M-1: the smoothed p-values' ranks are i.i.d. uniform -/

section Sequential

variable {Z α : Type*} [LinearOrder α]

/-- The probability that a unit's randomized rank in `C` is `r`: the units strictly below it, plus
a uniform draw among the units tied with it, itself included. This is the finite part of the
smoothed p-value `(L + θ T)/n`: its rank `L + ⌈θ T⌉`. -/
noncomputable def tieRankProb {ι : Type*} (C : Finset ι) (s : ι → α) (u : ι) (r : ℕ) : ℝ :=
  if rankLt C s u < r ∧ r ≤ rankLe C s u then 1 / ((rankLe C s u - rankLt C s u : ℕ) : ℝ) else 0

/-- A unit's rank splits into the units strictly below it and the units tied with it. -/
theorem rankLe_eq_rankLt_add {ι : Type*} (C : Finset ι) (s : ι → α) (u : ι) :
    rankLe C s u = rankLt C s u + (C.filter fun j => s j = s u).card := by
  classical
  unfold rankLe rankLt
  rw [← card_union_of_disjoint]
  · congr 1
    ext j
    simp only [mem_filter, mem_union]
    constructor
    · rintro ⟨hj, hle⟩
      rcases hle.lt_or_eq with h | h
      · exact Or.inl ⟨hj, h⟩
      · exact Or.inr ⟨hj, h⟩
    · rintro (⟨hj, h⟩ | ⟨hj, h⟩)
      · exact ⟨hj, h.le⟩
      · exact ⟨hj, h.le⟩
  · rw [disjoint_filter]
    intro j _ h₁ h₂
    exact (ne_of_lt h₁) h₂

/-- **The tie class covering a rank.** For each rank `1, …, |C|`, the units' probabilities of
holding it sum to 1: one tie class straddles it, and its `T` members each hold it with
probability `1/T`. -/
theorem sum_tieRankProb {ι : Type*} (C : Finset ι) (s : ι → α) {r : ℕ} (hr : 1 ≤ r)
    (hrC : r ≤ C.card) : ∑ u ∈ C, tieRankProb C s u r = 1 := by
  classical
  set A := C.filter fun u => rankLt C s u < r ∧ r ≤ rankLe C s u with hA
  obtain ⟨u₀, hu₀⟩ : A.Nonempty := by
    have h₁ := card_rankLe_le C s (r - 1)
    have h₂ := le_card_rankLt_lt C s hrC
    obtain ⟨u, hu, hnot⟩ := exists_mem_notMem_of_card_lt_card
      (s := C.filter fun u => rankLe C s u ≤ r - 1) (t := C.filter fun u => rankLt C s u < r)
      (by omega)
    refine ⟨u, mem_filter.2 ⟨(mem_filter.1 hu).1, (mem_filter.1 hu).2, ?_⟩⟩
    by_contra h
    exact hnot (mem_filter.2 ⟨(mem_filter.1 hu).1, by omega⟩)
  have hu₀' := mem_filter.1 hu₀
  have hclass : A = C.filter fun v => s v = s u₀ := by
    ext v
    simp only [hA, mem_filter]
    constructor
    · rintro ⟨hv, hlt, hle⟩
      refine ⟨hv, ?_⟩
      rcases lt_trichotomy (s v) (s u₀) with h | h | h
      · exfalso
        have : rankLe C s v ≤ rankLt C s u₀ := card_le_card fun j hj => by
          rw [mem_filter] at hj ⊢
          exact ⟨hj.1, lt_of_le_of_lt hj.2 h⟩
        omega
      · exact h
      · exfalso
        have : rankLe C s u₀ ≤ rankLt C s v := card_le_card fun j hj => by
          rw [mem_filter] at hj ⊢
          exact ⟨hj.1, lt_of_le_of_lt hj.2 h⟩
        omega
    · rintro ⟨hv, heq⟩
      have h1 : rankLt C s v = rankLt C s u₀ := by unfold rankLt; rw [heq]
      have h2 : rankLe C s v = rankLe C s u₀ := by unfold rankLe; rw [heq]
      rw [h1, h2]
      exact ⟨hv, hu₀'.2⟩
  have hT : ∀ v ∈ A, rankLe C s v - rankLt C s v = A.card := by
    intro v hv
    have hv' := mem_filter.1 (hclass ▸ hv)
    rw [rankLe_eq_rankLt_add C s v, hclass]
    have : (C.filter fun j => s j = s v) = C.filter fun j => s j = s u₀ := by rw [hv'.2]
    rw [this]
    omega
  have hpos : (0 : ℝ) < A.card := by exact_mod_cast card_pos.2 ⟨u₀, hu₀⟩
  unfold tieRankProb
  rw [← sum_filter, ← hA, sum_congr rfl fun v hv => by rw [hT v hv], sum_const, nsmul_eq_mul]
  field_simp

omit [LinearOrder α] in
theorem tieRankProb_comp {ι : Type*} [Fintype ι] [DecidableEq ι] [LinearOrder α] {C : Finset ι}
    {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (s : ι → α) (u : ι) (r : ℕ) :
    tieRankProb C (s ∘ σ) u r = tieRankProb C s (σ u) r := by
  simp only [tieRankProb, rankLe_comp hσ, rankLt_comp hσ]

/-- The probability that the randomized ranks, newest first, are `r`, for the units `v` in the
order given (`v 0` the newest). At each time the newest unit is ranked among every unit so far,
its ties broken by a fresh uniform draw; `S n` scores the first `n` units. -/
noncomputable def seqRankProb (S : (n : ℕ) → (Fin n → Z) → Fin n → α) :
    (n : ℕ) → (Fin n → Z) → List ℕ → ℝ
  | 0, _, r => if r = [] then 1 else 0
  | _ + 1, _, [] => 0
  | n + 1, v, r₀ :: r' =>
    tieRankProb univ (S (n + 1) v) 0 r₀ * seqRankProb S n (fun i => v i.succ) r'

/-- The rank sequences, newest first, that `n` units can produce: the newest of `m` units has a
rank in `1, …, m`. -/
def ValidRanks : ℕ → List ℕ → Prop
  | 0, r => r = []
  | _ + 1, [] => False
  | n + 1, r₀ :: r' => 1 ≤ r₀ ∧ r₀ ≤ n + 1 ∧ ValidRanks n r'

theorem mem_permsOn_univ {ι : Type*} [Fintype ι] [DecidableEq ι] (σ : Equiv.Perm ι) :
    σ ∈ permsOn (univ : Finset ι) :=
  mem_filter.2 ⟨mem_univ _, fun i hi => absurd (mem_univ i) hi⟩

/-- **M-1, the finite core.** For scores equivariant at every time, and over every ordering of a
fixed data set, each rank sequence that `n` units can produce has total probability 1, so
probability `1/n!` under a uniform ordering: the randomized ranks are independent, the newest of
`m` units uniform on `1, …, m`. The premises are the definition's: every earlier unit ranked,
calibration and test alike; ties broken by a fresh uniform draw; the scores at each time a
symmetric function of the units so far. -/
theorem seqRank_uniform (S : (n : ℕ) → (Fin n → Z) → Fin n → α)
    (hS : ∀ n, Equivariant (univ : Finset (Fin n)) (S n)) :
    ∀ (n : ℕ) (v : Fin n → Z) (r : List ℕ), ValidRanks n r →
      ∑ σ : Equiv.Perm (Fin n), seqRankProb S n (v ∘ σ) r = 1 := by
  intro n
  induction n with
  | zero =>
    intro v r hr
    simp only [ValidRanks] at hr
    subst hr
    simp [seqRankProb]
  | succ n ih =>
    intro v r hr
    cases r with
    | nil => exact absurd hr (by simp [ValidRanks])
    | cons r₀ r' =>
      obtain ⟨h₁, h₂, hr'⟩ := hr
      rw [← Equiv.sum_comp Equiv.Perm.decomposeFin.symm, Fintype.sum_prod_type]
      have hstep : ∀ (p : Fin (n + 1)) (e : Equiv.Perm (Fin n)),
          seqRankProb S (n + 1) (v ∘ Equiv.Perm.decomposeFin.symm (p, e)) (r₀ :: r') =
            tieRankProb univ (S (n + 1) v) p r₀ *
              seqRankProb S n ((fun i => v (Equiv.swap 0 p i.succ)) ∘ e) r' := by
        intro p e
        have hσ := mem_permsOn_univ (Equiv.Perm.decomposeFin.symm (p, e))
        have hfun : S (n + 1) (v ∘ Equiv.Perm.decomposeFin.symm (p, e)) =
            S (n + 1) v ∘ Equiv.Perm.decomposeFin.symm (p, e) :=
          funext fun i => hS (n + 1) _ hσ v i
        simp only [seqRankProb]
        rw [hfun, tieRankProb_comp hσ, Equiv.Perm.decomposeFin_symm_apply_zero]
        congr 1
        congr 1
        funext i
        simp [Equiv.Perm.decomposeFin_symm_apply_succ]
      simp only [hstep, ← mul_sum]
      rw [sum_congr rfl fun p _ => by rw [ih _ r' hr', mul_one]]
      exact sum_tieRankProb univ (S (n + 1) v) h₁ (by simpa using h₂)

/-- **M-1 under an exchangeable law.** For a law on data sets invariant under every reordering,
each rank sequence `n` units can produce has probability `1/n!`. -/
theorem seqRank_law [Fintype Z] (S : (n : ℕ) → (Fin n → Z) → Fin n → α)
    (hS : ∀ n, Equivariant (univ : Finset (Fin n)) (S n)) (n : ℕ) (P : (Fin n → Z) → ℝ)
    (hP1 : ∑ z, P z = 1) (hP : ∀ σ : Equiv.Perm (Fin n), ∀ z, P (z ∘ σ) = P z) (r : List ℕ)
    (hr : ValidRanks n r) :
    ∑ z, P z * seqRankProb S n z r = 1 / (n.factorial : ℝ) := by
  rw [sum_law_average (C := univ) P (fun σ _ z => hP σ z)]
  have hcard : ((permsOn (univ : Finset (Fin n))).card : ℝ) = n.factorial := by
    have : permsOn (univ : Finset (Fin n)) = univ :=
      eq_univ_of_forall fun σ => mem_permsOn_univ σ
    rw [this, card_univ, Fintype.card_perm, Fintype.card_fin]
  have hsum : ∀ z : Fin n → Z, ∑ σ ∈ permsOn (univ : Finset (Fin n)), seqRankProb S n (z ∘ σ) r =
      1 := by
    intro z
    have : permsOn (univ : Finset (Fin n)) = univ :=
      eq_univ_of_forall fun σ => mem_permsOn_univ σ
    rw [this]
    exact seqRank_uniform S hS n z r hr
  simp only [hsum, hcard]
  rw [← sum_mul, hP1, one_mul]

end Sequential

/-! ### E35, M-2: the betting product and Ville's inequality -/

section Betting

open MeasureTheory
open scoped ENNReal

variable (μ : Measure ℝ) (f : (n : ℕ) → (Fin n → ℝ) → ℝ → ℝ≥0∞)

/-- The capital after `n` bets on the p-values `h`: the product of the betting functions, the
`n`th chosen from the earlier p-values `Fin.init h` and read at the newest. -/
noncomputable def capital : (n : ℕ) → (Fin n → ℝ) → ℝ≥0∞
  | 0, _ => 1
  | n + 1, h => capital n (Fin.init h) * f n (Fin.init h) (h (Fin.last n))

theorem measurable_init (n : ℕ) : Measurable fun h : Fin (n + 1) → ℝ => Fin.init h :=
  measurable_pi_lambda _ fun i => measurable_pi_apply (Fin.castSucc i)

variable {f}

theorem measurable_capital (hf : ∀ n, Measurable (Function.uncurry (f n))) :
    ∀ n, Measurable (capital f n)
  | 0 => measurable_const
  | n + 1 => by
    have h₁ := (measurable_capital hf n).comp (measurable_init n)
    have h₂ := (hf n).comp ((measurable_init n).prodMk (measurable_pi_apply (Fin.last n)))
    exact h₁.mul h₂

/-- **The link, one step.** Given the history `h`, the next capital integrates over the newest
p-value to the current one: the martingale property, with the betting function chosen from `h`
and integrating to 1. -/
theorem lintegral_capital_snoc (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ = 1) (n : ℕ) (h : Fin n → ℝ) :
    ∫⁻ x, capital f (n + 1) (Fin.snoc h x) ∂μ = capital f n h := by
  simp only [capital, Fin.init_snoc, Fin.snoc_last]
  have hm : Measurable fun x => f n h x := (hf n).comp (measurable_const.prodMk measurable_id)
  rw [lintegral_const_mul _ hm, h1, mul_one]

/-- Tonelli at the newest p-value: the capital's integral at `n + 1` is the history's integral of
the next step's integral over the newest p-value. -/
theorem lintegral_capital_succ [IsProbabilityMeasure μ]
    (hf : ∀ n, Measurable (Function.uncurry (f n))) (n : ℕ) :
    ∫⁻ h, capital f (n + 1) h ∂(Measure.pi fun _ : Fin (n + 1) => μ) =
      ∫⁻ g, ∫⁻ x, capital f (n + 1) (Fin.snoc g x) ∂μ ∂(Measure.pi fun _ : Fin n => μ) := by
  set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) (Fin.last n)
  have hmp : MeasurePreserving e (Measure.pi fun _ : Fin (n + 1) => μ)
      (μ.prod (Measure.pi fun _ : Fin n => μ)) :=
    measurePreserving_piFinSuccAbove (fun _ : Fin (n + 1) => μ) (Fin.last n)
  have hG : Measurable fun z : ℝ × (Fin n → ℝ) => capital f (n + 1) (e.symm z) :=
    (measurable_capital hf (n + 1)).comp e.symm.measurable
  have hsnoc : ∀ (g : Fin n → ℝ) (x : ℝ), e.symm (x, g) = Fin.snoc g x := fun g x => by
    ext i
    simp [e, MeasurableEquiv.piFinSuccAbove_symm_apply]
  calc ∫⁻ h, capital f (n + 1) h ∂(Measure.pi fun _ : Fin (n + 1) => μ)
      = ∫⁻ h, capital f (n + 1) (e.symm (e h)) ∂(Measure.pi fun _ : Fin (n + 1) => μ) := by
        simp only [MeasurableEquiv.symm_apply_apply]
    _ = ∫⁻ z, capital f (n + 1) (e.symm z) ∂(μ.prod (Measure.pi fun _ : Fin n => μ)) :=
        hmp.lintegral_comp_emb e.measurableEmbedding fun z => capital f (n + 1) (e.symm z)
    _ = ∫⁻ g, ∫⁻ x, capital f (n + 1) (e.symm (x, g)) ∂μ ∂(Measure.pi fun _ : Fin n => μ) :=
        lintegral_prod_symm _ hG.aemeasurable
    _ = ∫⁻ g, ∫⁻ x, capital f (n + 1) (Fin.snoc g x) ∂μ ∂(Measure.pi fun _ : Fin n => μ) := by
        simp only [hsnoc]

/-- **The link.** With the p-values independent with law `μ`, the capital at every horizon has
expectation 1. -/
theorem lintegral_capital [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ = 1) :
    ∀ n, ∫⁻ h, capital f n h ∂(Measure.pi fun _ : Fin n => μ) = 1
  | 0 => by simp [capital]
  | n + 1 => by
    rw [lintegral_capital_succ μ hf n]
    simp only [lintegral_capital_snoc μ hf h1 n]
    exact lintegral_capital hf h1 n

/-- Under a bet integrating to at most 1, the next capital integrates to at most the current one:
the capital is a supermartingale. -/
theorem lintegral_capital_snoc_le (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) (n : ℕ) (h : Fin n → ℝ) :
    ∫⁻ x, capital f (n + 1) (Fin.snoc h x) ∂μ ≤ capital f n h := by
  simp only [capital, Fin.init_snoc, Fin.snoc_last]
  have hm : Measurable fun x => f n h x := (hf n).comp (measurable_const.prodMk measurable_id)
  rw [lintegral_const_mul _ hm]
  calc capital f n h * ∫⁻ x, f n h x ∂μ ≤ capital f n h * 1 := by gcongr; exact h1 n h
    _ = capital f n h := mul_one _

/-- Under bets integrating to at most 1, the capital has expectation at most 1. -/
theorem lintegral_capital_le [IsProbabilityMeasure μ]
    (hf : ∀ n, Measurable (Function.uncurry (f n))) (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) :
    ∀ n, ∫⁻ h, capital f n h ∂(Measure.pi fun _ : Fin n => μ) ≤ 1
  | 0 => by simp [capital]
  | n + 1 => by
    rw [lintegral_capital_succ μ hf n]
    calc ∫⁻ g, ∫⁻ x, capital f (n + 1) (Fin.snoc g x) ∂μ ∂(Measure.pi fun _ : Fin n => μ)
        ≤ ∫⁻ g, capital f n g ∂(Measure.pi fun _ : Fin n => μ) :=
          lintegral_mono fun g => lintegral_capital_snoc_le μ hf h1 n g
      _ ≤ 1 := lintegral_capital_le hf h1 n

/-- The times stopped at, by `n`: the rule `stop k` reads the first `k` p-values at time `k`. -/
def stoppedSet (stop : (n : ℕ) → Set (Fin n → ℝ)) : (n : ℕ) → Set (Fin n → ℝ)
  | 0 => stop 0
  | n + 1 => (fun h => Fin.init h) ⁻¹' stoppedSet stop n ∪ stop (n + 1)

open scoped Classical in
variable (f) in
/-- The bets switched off once stopped: still betting functions, still reading only the past. -/
noncomputable def stoppedBet (stop : (n : ℕ) → Set (Fin n → ℝ)) (n : ℕ) (h : Fin n → ℝ) (x : ℝ) :
    ℝ≥0∞ :=
  if h ∈ stoppedSet stop n then 1 else f n h x

theorem measurableSet_stoppedSet {stop : (n : ℕ) → Set (Fin n → ℝ)}
    (hstop : ∀ n, MeasurableSet (stop n)) : ∀ n, MeasurableSet (stoppedSet stop n)
  | 0 => hstop 0
  | n + 1 => ((measurableSet_stoppedSet hstop n).preimage (measurable_init n)).union (hstop (n + 1))

theorem measurable_stoppedBet (hf : ∀ n, Measurable (Function.uncurry (f n)))
    {stop : (n : ℕ) → Set (Fin n → ℝ)} (hstop : ∀ n, MeasurableSet (stop n)) (n : ℕ) :
    Measurable (Function.uncurry (stoppedBet f stop n)) := by
  classical
  have hS : MeasurableSet {z : (Fin n → ℝ) × ℝ | z.1 ∈ stoppedSet stop n} :=
    (measurableSet_stoppedSet hstop n).preimage measurable_fst
  exact Measurable.ite hS measurable_const (hf n)

theorem lintegral_stoppedBet [IsProbabilityMeasure μ] (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ = 1)
    (stop : (n : ℕ) → Set (Fin n → ℝ)) (n : ℕ) (h : Fin n → ℝ) :
    ∫⁻ x, stoppedBet f stop n h x ∂μ = 1 := by
  classical
  by_cases hs : h ∈ stoppedSet stop n
  · simp [stoppedBet, hs]
  · simp only [stoppedBet, hs, if_false]
    exact h1 n h

theorem lintegral_stoppedBet_le [IsProbabilityMeasure μ] (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1)
    (stop : (n : ℕ) → Set (Fin n → ℝ)) (n : ℕ) (h : Fin n → ℝ) :
    ∫⁻ x, stoppedBet f stop n h x ∂μ ≤ 1 := by
  classical
  by_cases hs : h ∈ stoppedSet stop n
  · simp [stoppedBet, hs]
  · simp only [stoppedBet, hs, if_false]
    exact h1 n h

/-- **The capital at a stopping time is an e-value.** Stopped by any measurable rule, by the
horizon `n` at the latest, the capital has expectation 1. -/
theorem lintegral_stoppedCapital [IsProbabilityMeasure μ]
    (hf : ∀ n, Measurable (Function.uncurry (f n))) (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ = 1)
    {stop : (n : ℕ) → Set (Fin n → ℝ)} (hstop : ∀ n, MeasurableSet (stop n)) (n : ℕ) :
    ∫⁻ h, capital (stoppedBet f stop) n h ∂(Measure.pi fun _ : Fin n => μ) = 1 :=
  lintegral_capital μ (measurable_stoppedBet hf hstop) (lintegral_stoppedBet μ h1 stop) n

/-- Until it stops, the stopped capital is the capital. -/
theorem capital_stoppedBet_of_notMem {stop : (n : ℕ) → Set (Fin n → ℝ)} :
    ∀ (n : ℕ) (h : Fin n → ℝ), h ∉ stoppedSet stop n →
      capital (stoppedBet f stop) n h = capital f n h
  | 0, _, _ => rfl
  | n + 1, h, hs => by
    classical
    have hinit : Fin.init h ∉ stoppedSet stop n := fun hh => hs (Or.inl hh)
    simp only [capital, stoppedBet, hinit, if_false]
    rw [capital_stoppedBet_of_notMem n (Fin.init h) hinit]

/-- Once stopped, the stopped capital stays put. -/
theorem capital_stoppedBet_of_mem {stop : (n : ℕ) → Set (Fin n → ℝ)} (n : ℕ)
    (h : Fin (n + 1) → ℝ) (hs : Fin.init h ∈ stoppedSet stop n) :
    capital (stoppedBet f stop) (n + 1) h = capital (stoppedBet f stop) n (Fin.init h) := by
  classical
  simp only [capital, stoppedBet, hs, if_true, mul_one]

/-- Stopped by the alarm `c ≤ capital`, the stopped capital is at least `c` from then on. -/
theorem le_capital_stoppedBet (c : ℝ≥0∞) :
    ∀ (n : ℕ) (h : Fin n → ℝ), h ∈ stoppedSet (fun k => {h | c ≤ capital f k h}) n →
      c ≤ capital (stoppedBet f fun k => {h | c ≤ capital f k h}) n h
  | 0, _, hs => hs
  | n + 1, h, hs => by
    classical
    by_cases hi : Fin.init h ∈ stoppedSet (fun k => {h | c ≤ capital f k h}) n
    · rw [capital_stoppedBet_of_mem n h hi]
      exact le_capital_stoppedBet c n (Fin.init h) hi
    · have hs' : c ≤ capital f (n + 1) h := by
        rcases hs with hs | hs
        · exact absurd hs hi
        · exact hs
      have heq : capital (stoppedBet f fun k => {h | c ≤ capital f k h}) (n + 1) h =
          capital f (n + 1) h := by
        simp only [capital, stoppedBet, hi, if_false]
        rw [capital_stoppedBet_of_notMem n (Fin.init h) hi]
      rw [heq]
      exact hs'

/-- **M-2 at a horizon, Ville's inequality.** With the p-values independent with law `μ` and each
bet chosen from the past and integrating to at most 1, the probability that the capital reaches
`c` by the horizon `n` is at most `1/c`. The bet is switched off at the alarm, the stopped capital
is again a betting product with expectation at most 1, and Markov's inequality finishes. -/
theorem ville_horizon [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) (c : ℝ≥0∞) (n : ℕ) :
    (Measure.pi fun _ : Fin n => μ) (stoppedSet (fun k => {h | c ≤ capital f k h}) n) ≤ c⁻¹ := by
  have hstop : ∀ k, MeasurableSet {h : Fin k → ℝ | c ≤ capital f k h} := fun k =>
    measurableSet_le measurable_const (measurable_capital hf k)
  have hg := measurable_capital (measurable_stoppedBet hf hstop) n
  have hint := lintegral_capital_le μ (measurable_stoppedBet hf hstop)
    (lintegral_stoppedBet_le μ h1 _) n
  have hmarkov := mul_meas_ge_le_lintegral₀ (μ := Measure.pi fun _ : Fin n => μ) hg.aemeasurable c
  have hsub : (Measure.pi fun _ : Fin n => μ) (stoppedSet (fun k => {h | c ≤ capital f k h}) n) ≤
      (Measure.pi fun _ : Fin n => μ)
        {h | c ≤ capital (stoppedBet f fun k => {h | c ≤ capital f k h}) n h} :=
    measure_mono fun h hh => le_capital_stoppedBet c n h hh
  rw [ENNReal.le_inv_iff_mul_le, mul_comm]
  calc c * (Measure.pi fun _ : Fin n => μ) (stoppedSet (fun k => {h | c ≤ capital f k h}) n)
      ≤ c * (Measure.pi fun _ : Fin n => μ)
          {h | c ≤ capital (stoppedBet f fun k => {h | c ≤ capital f k h}) n h} := by gcongr
    _ ≤ 1 := hmarkov.trans hint

/-- The stopped set read along a path: stopped by `N` exactly when the rule fired at some
`k ≤ N`. -/
theorem mem_stoppedSet_iff (stop : (n : ℕ) → Set (Fin n → ℝ)) (x : ℕ → ℝ) :
    ∀ N, (fun i : Fin N => x i) ∈ stoppedSet stop N ↔ ∃ k ≤ N, (fun i : Fin k => x i) ∈ stop k
  | 0 => by simp [stoppedSet]
  | N + 1 => by
    have hinit : Fin.init (fun i : Fin (N + 1) => x i) = fun i : Fin N => x i := rfl
    simp only [stoppedSet, Set.mem_union, Set.mem_preimage, hinit, mem_stoppedSet_iff stop x N]
    constructor
    · rintro (⟨k, hk, h⟩ | h)
      · exact ⟨k, by omega, h⟩
      · exact ⟨N + 1, le_rfl, h⟩
    · rintro ⟨k, hk, h⟩
      rcases Nat.lt_or_ge k (N + 1) with hlt | hge
      · exact Or.inl ⟨k, by omega, h⟩
      · obtain rfl : k = N + 1 := by omega
        exact Or.inr h

open Filter Topology in
/-- **M-2, Ville's inequality for the betting product.** Let the p-values `p i` on a space `Ω` be
independent with law `μ` at every horizon (M-1's conclusion, `μ` uniform on `[0, 1]`), and let each
bet `f n` be measurable, integrate to at most 1 under `μ`, and read only the earlier p-values.
Then the capital ever reaches `c` with probability at most `1/c`. -/
theorem ville [IsProbabilityMeasure μ] (hf : ∀ n, Measurable (Function.uncurry (f n)))
    (h1 : ∀ n h, ∫⁻ x, f n h x ∂μ ≤ 1) {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω)
    (p : ℕ → Ω → ℝ) (hp : ∀ i, Measurable (p i))
    (hlaw : ∀ n, P.map (fun ω (i : Fin n) => p i ω) = Measure.pi fun _ => μ) (c : ℝ≥0∞) :
    P {ω | ∃ n, c ≤ capital f n (fun i : Fin n => p i ω)} ≤ c⁻¹ := by
  set stop : (n : ℕ) → Set (Fin n → ℝ) := fun k => {h | c ≤ capital f k h} with hstopdef
  have hstop : ∀ k, MeasurableSet (stop k) := fun k =>
    measurableSet_le measurable_const (measurable_capital hf k)
  set E : ℕ → Set Ω := fun N => {ω | (fun i : Fin N => p i ω) ∈ stoppedSet stop N} with hEdef
  have hE : {ω | ∃ n, c ≤ capital f n (fun i : Fin n => p i ω)} = ⋃ N, E N := by
    ext ω
    simp only [Set.mem_setOf_eq, Set.mem_iUnion, hEdef, mem_stoppedSet_iff stop (fun i => p i ω)]
    constructor
    · rintro ⟨n, hn⟩
      exact ⟨n, n, le_rfl, hn⟩
    · rintro ⟨N, k, -, hk⟩
      exact ⟨k, hk⟩
  have hmono : Monotone E := by
    intro N M hNM ω hω
    simp only [hEdef, Set.mem_setOf_eq, mem_stoppedSet_iff stop (fun i => p i ω)] at hω ⊢
    obtain ⟨k, hk, h⟩ := hω
    exact ⟨k, hk.trans hNM, h⟩
  have hbound : ∀ N, P (E N) ≤ c⁻¹ := by
    intro N
    have hφ : Measurable fun ω (i : Fin N) => p i ω := measurable_pi_lambda _ fun i => hp i
    have : P (E N) = (P.map fun ω (i : Fin N) => p i ω) (stoppedSet stop N) := by
      rw [Measure.map_apply hφ (measurableSet_stoppedSet hstop N)]
      rfl
    rw [this, hlaw N]
    exact ville_horizon μ hf h1 c N
  rw [hE]
  exact le_of_tendsto' (tendsto_measure_iUnion_atTop hmono) hbound

end Betting

/-! ### Slice 9: the jackknife+ and the cross-fitted rule -/

section CrossFitted

variable {ι : Type*}

/-- The units of `C` that `u` is stranger than: those `j` with `s u j < s j u`, where `s i j` is
unit `i`'s conformity under the fit that leaves out `i` and `j`. A tie makes neither stranger. -/
def strangerThan {α : Type*} [LinearOrder α] (C : Finset ι) (s : ι → ι → α) (u : ι) : Finset ι :=
  C.filter fun j => s u j < s j u

/-- **The cross-fitted accept rule at level `k`.** `u` is refused when it is stranger than at least
`|C| − k` units of `C`. On leave-two-out scores this is the `n`-fold cross-conformal set, and the
jackknife+ interval contains it (`not_acceptsPairwise_of_jackRefuses`). -/
def AcceptsPairwise {α : Type*} [LinearOrder α] (k : ℕ) (C : Finset ι) (s : ι → ι → α) (u : ι) :
    Prop :=
  (strangerThan C s u).card + k < C.card

instance {α : Type*} [LinearOrder α] (k : ℕ) (C : Finset ι) (s : ι → ι → α) (u : ι) :
    Decidable (AcceptsPairwise k C s u) :=
  inferInstanceAs (Decidable (_ < _))

theorem notMem_strangerThan_self {α : Type*} [LinearOrder α] (C : Finset ι) (s : ι → ι → α)
    (u : ι) : u ∉ strangerThan C s u :=
  fun h => lt_irrefl _ (mem_filter.1 h).2

/-- At level `0` the cross-fitted rule accepts every unit of `C`: the set is the whole space. -/
theorem acceptsPairwise_zero {α : Type*} [LinearOrder α] {C : Finset ι} (s : ι → ι → α) {u : ι}
    (hu : u ∈ C) : AcceptsPairwise 0 C s u := by
  have h : strangerThan C s u ⊂ C :=
    (ssubset_iff_of_subset (filter_subset _ _)).2 ⟨u, hu, notMem_strangerThan_self C s u⟩
  have := card_lt_card h
  unfold AcceptsPairwise
  omega

/-- **The counting core**, Landau's argument as Barber, Candès, Ramdas and Tibshirani use it.
Whatever the comparisons and their ties, at most `2k − 1` units of `C` are refused at level `k`:
each refused unit is stranger than at least `|C| − k` units, and among the refused units each
comparison goes at most one way. -/
theorem card_refused_pairwise_le {α : Type*} [LinearOrder α] (C : Finset ι) (s : ι → ι → α)
    (k : ℕ) : (C.filter fun u => ¬ AcceptsPairwise k C s u).card ≤ 2 * k - 1 := by
  classical
  set F := C.filter fun u => ¬ AcceptsPairwise k C s u
  have hFC : F ⊆ C := filter_subset _ _
  set a : ι → ℕ := fun u => (F.filter fun j => s u j < s j u).card with ha
  have hunit : ∀ u ∈ F, C.card ≤ a u + (C \ F).card + k := by
    intro u hu
    have href : C.card ≤ (strangerThan C s u).card + k := by
      simpa only [AcceptsPairwise, not_lt] using (mem_filter.1 hu).2
    have hsplit : (strangerThan C s u).card ≤ a u + (C \ F).card := by
      calc (strangerThan C s u).card ≤ ((F.filter fun j => s u j < s j u) ∪ (C \ F)).card := by
            refine card_le_card fun j hj => ?_
            rw [strangerThan, mem_filter] at hj
            by_cases hjF : j ∈ F
            · exact mem_union_left _ (mem_filter.2 ⟨hjF, hj.2⟩)
            · exact mem_union_right _ (mem_sdiff.2 ⟨hj.1, hjF⟩)
        _ ≤ _ := card_union_le _ _
    omega
  have hcard : C.card = F.card + (C \ F).card := by
    rw [card_sdiff_of_subset hFC]
    have := card_le_card hFC
    omega
  have hsum : F.card * C.card ≤ (∑ u ∈ F, a u) + F.card * (C \ F).card + F.card * k := by
    have h := sum_le_sum hunit
    simp only [sum_const, smul_eq_mul, sum_add_distrib] at h
    linarith
  have hle : ∀ u j, ((if s u j < s j u then 1 else 0) + (if s j u < s u j then 1 else 0) +
      (if u = j then 1 else 0) : ℕ) ≤ 1 := by
    intro u j
    by_cases huj : u = j
    · subst huj
      simp
    · by_cases h1 : s u j < s j u
      · have h2 : ¬ s j u < s u j := fun h2 => lt_asymm h1 h2
        simp [h1, h2, huj]
      · simp only [h1, huj, if_false]
        split_ifs <;> omega
  have hA : ∑ u ∈ F, a u = ∑ u ∈ F, ∑ j ∈ F, if s u j < s j u then 1 else 0 := by
    simp only [ha, card_filter]
  have hB : ∑ u ∈ F, a u = ∑ u ∈ F, ∑ j ∈ F, if s j u < s u j then 1 else 0 := by
    rw [hA, sum_comm]
  have hD : F.card = ∑ u ∈ F, ∑ j ∈ F, if u = j then 1 else 0 := by
    rw [card_eq_sum_ones]
    refine sum_congr rfl fun u hu => ?_
    rw [sum_ite_eq]
    simp [hu]
  have hpairs : 2 * (∑ u ∈ F, a u) + F.card ≤ F.card * F.card := by
    calc 2 * (∑ u ∈ F, a u) + F.card
        = ∑ u ∈ F, ∑ j ∈ F, ((if s u j < s j u then 1 else 0) + (if s j u < s u j then 1 else 0) +
            (if u = j then 1 else 0)) := by
          rw [two_mul]
          nth_rw 1 [hA]
          rw [hB, hD, ← sum_add_distrib, ← sum_add_distrib]
          refine sum_congr rfl fun u _ => ?_
          rw [← sum_add_distrib, ← sum_add_distrib]
      _ ≤ ∑ u ∈ F, ∑ j ∈ F, 1 := sum_le_sum fun u _ => sum_le_sum fun j _ => hle u j
      _ = F.card * F.card := by simp
  have h1 : F.card * F.card ≤ (∑ u ∈ F, a u) + F.card * k := by
    rw [hcard, mul_add] at hsum
    linarith
  have h2 : F.card * (F.card + 1) ≤ F.card * (2 * k) := by
    have : F.card * F.card + F.card ≤ 2 * (F.card * k) := by linarith
    calc F.card * (F.card + 1) = F.card * F.card + F.card := by ring
      _ ≤ 2 * (F.card * k) := this
      _ = F.card * (2 * k) := by ring
  rcases Nat.eq_zero_or_pos F.card with h0 | hpos
  · rw [h0]
    exact Nat.zero_le _
  · have := Nat.le_of_mul_le_mul_left h2 hpos
    omega

/-- The share of units the cross-fitted rule accepts: at least `1 − (2k − 1)/|C|`. -/
theorem ratio_acceptsPairwise_ge {α : Type*} [LinearOrder α] {C : Finset ι} (hC : C.Nonempty)
    (s : ι → ι → α) (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / C.card ≤
      ((C.filter fun u => AcceptsPairwise k C s u).card : ℝ) / C.card := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast hC.card_pos
  have hrej := card_refused_pairwise_le C s k
  have hsplit := C.card_filter_add_card_filter_not (fun u => AcceptsPairwise k C s u)
  have hcount : (C.card : ℝ) - ((2 * k - 1 : ℕ) : ℝ) ≤
      (C.filter fun u => AcceptsPairwise k C s u).card := by
    have : C.card ≤ (C.filter fun u => AcceptsPairwise k C s u).card + (2 * k - 1) := by omega
    have : (C.card : ℝ) ≤ ((C.filter fun u => AcceptsPairwise k C s u).card : ℝ) +
        ((2 * k - 1 : ℕ) : ℝ) := by exact_mod_cast this
    linarith
  rw [one_sub_div hc.ne']
  exact div_le_div_of_nonneg_right hcount hc.le

/-- The level in Barber et al.'s terms: at `k ≤ α m`, the bound `1 − (2k − 1)/m` is at least
`1 − 2α`, and above it by `1/m` at `k = α m ≥ 1`. -/
theorem one_sub_two_mul_le {k m : ℕ} {α : ℝ} (hm : 0 < m) (hk : (k : ℝ) ≤ α * m) :
    1 - 2 * α ≤ 1 - ((2 * k - 1 : ℕ) : ℝ) / m := by
  have hm' : (0 : ℝ) < m := by exact_mod_cast hm
  have h1 : ((2 * k - 1 : ℕ) : ℝ) ≤ 2 * k := by
    have : 2 * k - 1 ≤ 2 * k := Nat.sub_le _ _
    exact_mod_cast this
  have : ((2 * k - 1 : ℕ) : ℝ) / m ≤ 2 * α := by
    rw [div_le_iff₀ hm']
    linarith
  linarith

/-- The jackknife+'s conformity: minus a regression fit's absolute residual. `cov` and `resp` read a
unit's covariate and response. -/
def jackFit {Z X : Type*} (cov : Z → X) (resp : Z → ℝ) (alg : Finset ι → (ι → Z) → X → ℝ)
    (T : Finset ι) (z : ι → Z) (p : Z) : ℝ :=
  -|resp p - alg T z (cov p)|

variable [Fintype ι] [DecidableEq ι]

/-- A matrix score on data sets is **pairwise-equivariant on `C`** when reordering the units of `C`
reorders the matrix's rows and columns with them. -/
def PairEquivariant {Z α : Type*} (C : Finset ι) (S : (ι → Z) → ι → ι → α) : Prop :=
  ∀ σ ∈ permsOn C, ∀ z i j, S (z ∘ σ) i j = S z (σ i) (σ j)

theorem acceptsPairwise_comp {α : Type*} [LinearOrder α] {k : ℕ} {C : Finset ι}
    {σ : Equiv.Perm ι} (hσ : σ ∈ permsOn C) (s : ι → ι → α) (u : ι) :
    AcceptsPairwise k C (fun i j => s (σ i) (σ j)) u ↔ AcceptsPairwise k C s (σ u) := by
  have h : (strangerThan C (fun i j => s (σ i) (σ j)) u).card = (strangerThan C s (σ u)).card := by
    apply card_equiv σ
    intro j
    simp only [strangerThan, mem_filter, mem_iff_of_mem_permsOn hσ]
  simp only [AcceptsPairwise, h]

/-- **The cross-fitted rule over orderings.** For a pairwise-equivariant score, over every ordering
of a fixed data set, the test unit is accepted at level `k` at least `1 − (2k − 1)/|C|` of the time,
whatever the ties. -/
theorem coverage_pairwise_orderings {Z α : Type*} [LinearOrder α] {C : Finset ι} {t : ι}
    (ht : t ∈ C) {S : (ι → Z) → ι → ι → α} (hS : PairEquivariant C S) (z : ι → Z) (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι => AcceptsPairwise k C (S (z ∘ σ)) t).card : ℝ) /
        (permsOn C).card := by
  have hfilt : ((permsOn C).filter fun σ : Equiv.Perm ι => AcceptsPairwise k C (S (z ∘ σ)) t) =
      (permsOn C).filter fun σ : Equiv.Perm ι => AcceptsPairwise k C (S z) (σ t) := by
    refine filter_congr fun σ hσ => ?_
    have : S (z ∘ σ) = fun i j => S z (σ i) (σ j) :=
      funext fun i => funext fun j => hS σ hσ z i j
    rw [this, acceptsPairwise_comp hσ]
  rw [hfilt, ratio_permsOn ht (fun u => AcceptsPairwise k C (S z) u)]
  exact ratio_acceptsPairwise_ge ⟨t, ht⟩ (S z) k

/-- **The cross-fitted rule under an exchangeable law.** `W` carries the fits' randomness, and
`act σ` relabels it when the units of `C` are reordered by `σ`. The law of the randomness and the
data together is invariant under each reordering with its relabeling, and the matrix equivariant
under the pair. With `act` the identity, the law is exchangeable given `w` and the matrix
equivariant for each `w`; a relabeling admits fits symmetric only in law, such as a sampler whose
seeds travel with the units. Then the test unit is accepted with probability at least
`1 − (2k − 1)/|C|`. -/
theorem coverage_pairwise {W Z α : Type*} [Fintype W] [Fintype Z] [LinearOrder α] {C : Finset ι}
    {t : ι} (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z)
    (hP1 : ∑ w, ∑ z, P w z = 1) (act : Equiv.Perm ι → W ≃ W)
    (hP : ∀ σ ∈ permsOn C, ∀ w z, P (act σ w) (z ∘ σ) = P w z) (S : W → (ι → Z) → ι → ι → α)
    (hS : ∀ σ ∈ permsOn C, ∀ w z i j, S (act σ w) (z ∘ σ) i j = S w z (σ i) (σ j)) (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / C.card ≤
      ∑ w, ∑ z, P w z * if AcceptsPairwise k C (S w z) t then 1 else 0 := by
  rw [← Fintype.sum_prod_type'] at hP1 ⊢
  have hG : ((permsOn C).card : ℝ) ≠ 0 := by exact_mod_cast (card_permsOn_pos C).ne'
  have hswap : ∀ σ ∈ permsOn C,
      ∑ x : W × (ι → Z), P x.1 x.2 * (if AcceptsPairwise k C (S x.1 x.2) t then (1 : ℝ) else 0) =
        ∑ x : W × (ι → Z),
          P x.1 x.2 * (if AcceptsPairwise k C (S x.1 x.2) (σ t) then 1 else 0) := by
    intro σ hσ
    rw [← Equiv.sum_comp ((act σ).prodCongr (reorder σ))]
    refine Fintype.sum_congr _ _ fun x => ?_
    simp only [Equiv.prodCongr_apply, Prod.map_fst, Prod.map_snd]
    rw [show reorder σ x.2 = x.2 ∘ σ from rfl, hP σ hσ]
    have : S (act σ x.1) (x.2 ∘ σ) = fun i j => S x.1 x.2 (σ i) (σ j) :=
      funext fun i => funext fun j => hS σ hσ x.1 x.2 i j
    rw [this]
    simp only [acceptsPairwise_comp hσ]
  have key : ∀ x : W × (ι → Z), (1 - ((2 * k - 1 : ℕ) : ℝ) / C.card) * P x.1 x.2 ≤
      P x.1 x.2 * ((∑ σ ∈ permsOn C,
        if AcceptsPairwise k C (S x.1 x.2) (σ t) then (1 : ℝ) else 0) / (permsOn C).card) := by
    intro x
    rw [mul_comm]
    refine mul_le_mul_of_nonneg_left ?_ (hP0 x.1 x.2)
    rw [sum_boole, ratio_permsOn ht (fun u => AcceptsPairwise k C (S x.1 x.2) u)]
    exact ratio_acceptsPairwise_ge ⟨t, ht⟩ _ k
  calc 1 - ((2 * k - 1 : ℕ) : ℝ) / C.card
      = ∑ x : W × (ι → Z), (1 - ((2 * k - 1 : ℕ) : ℝ) / C.card) * P x.1 x.2 := by
        rw [← mul_sum, hP1, mul_one]
    _ ≤ ∑ x : W × (ι → Z), P x.1 x.2 * ((∑ σ ∈ permsOn C,
          if AcceptsPairwise k C (S x.1 x.2) (σ t) then (1 : ℝ) else 0) / (permsOn C).card) :=
        sum_le_sum fun x _ => key x
    _ = (∑ σ ∈ permsOn C, ∑ x : W × (ι → Z),
          P x.1 x.2 * if AcceptsPairwise k C (S x.1 x.2) (σ t) then 1 else 0) /
          (permsOn C).card := by
        simp only [mul_div_assoc', mul_sum, sum_div]
        exact sum_comm
    _ = (∑ σ ∈ permsOn C, ∑ x : W × (ι → Z),
          P x.1 x.2 * if AcceptsPairwise k C (S x.1 x.2) t then 1 else 0) /
          (permsOn C).card := by
        rw [sum_congr rfl fun σ hσ => (hswap σ hσ).symm]
    _ = _ := by
        rw [sum_const, nsmul_eq_mul]
        field_simp

/-- Leaving out `i` and `j`, renamed by `σ`, leaves out `σ i` and `σ j`. -/
theorem image_sdiff_pair (σ : Equiv.Perm ι) (i j : ι) :
    (univ \ {i, j}).image σ = univ \ {σ i, σ j} := by
  ext x
  simp only [mem_image, mem_sdiff, mem_univ, true_and, mem_insert, mem_singleton]
  constructor
  · rintro ⟨a, ha, rfl⟩ h
    exact ha (h.imp (fun h => σ.injective h) (fun h => σ.injective h))
  · intro hx
    refine ⟨σ.symm x, fun h => hx ?_, σ.apply_symm_apply x⟩
    rcases h with h | h
    · left
      rw [← h, σ.apply_symm_apply]
    · right
      rw [← h, σ.apply_symm_apply]

/-- The leave-two-out score: unit `i` read by the fit to every unit but `i` and `j`. -/
def pairScore {Z β : Type*} (fit : Finset ι → (ι → Z) → Z → β) (z : ι → Z) (i j : ι) : β :=
  fit (univ \ {i, j}) z (z i)

/-- A **symmetric fit** reads its training units as a set: on the data reordered by `σ`, the fit
to the units `T` is the fit to the units `T.image σ` of the original, which hold the same data. -/
def FitSymmetric {Z β : Type*} (C : Finset ι) (fit : Finset ι → (ι → Z) → β) : Prop :=
  ∀ σ ∈ permsOn C, ∀ T z, fit T (z ∘ σ) = fit (T.image σ) z

theorem pairScore_comp {Z β : Type*} {fit fit' : Finset ι → (ι → Z) → Z → β} {σ : Equiv.Perm ι}
    {z : ι → Z} (h : ∀ T, fit' T (z ∘ σ) = fit (T.image σ) z) (i j : ι) :
    pairScore fit' (z ∘ σ) i j = pairScore fit z (σ i) (σ j) := by
  simp only [pairScore, h, image_sdiff_pair, Function.comp_apply]

/-- A symmetric fit's leave-two-out score is pairwise-equivariant. -/
theorem pairScore_equivariant {Z β : Type*} {C : Finset ι} {fit : Finset ι → (ι → Z) → Z → β}
    (h : FitSymmetric C fit) : PairEquivariant C (pairScore fit) :=
  fun σ hσ z i j => pairScore_comp (fun T => h σ hσ T z) i j

/-- **The cross-fitted rule with any score** (Vovk's `n`-fold cross-conformal set). Each unit is
scored by a fit to every unit but itself and the other unit of the pair, under any conformity score
the fit defines; the fit is symmetric up to the relabeling `act` of its randomness. Then the test
unit is accepted with probability at least `1 − (2k − 1)/|C|`. -/
theorem crossFit_coverage {W Z α : Type*} [Fintype W] [Fintype Z] [LinearOrder α] {C : Finset ι}
    {t : ι} (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z)
    (hP1 : ∑ w, ∑ z, P w z = 1) (act : Equiv.Perm ι → W ≃ W)
    (hP : ∀ σ ∈ permsOn C, ∀ w z, P (act σ w) (z ∘ σ) = P w z)
    (fit : W → Finset ι → (ι → Z) → Z → α)
    (hfit : ∀ σ ∈ permsOn C, ∀ w T z, fit (act σ w) T (z ∘ σ) = fit w (T.image σ) z) (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / C.card ≤
      ∑ w, ∑ z, P w z * if AcceptsPairwise k C (pairScore (fit w) z) t then 1 else 0 :=
  coverage_pairwise ht P hP0 hP1 act hP (fun w => pairScore (fit w))
    (fun σ hσ w z i j => pairScore_comp (fun T => hfit σ hσ w T z) i j) k

/-- **The jackknife+ refuses** the test unit `t` at level `k` when at least `|C| − k` of the upper
ends `μ̂₋ᵢ(Xₜ) + Rᵢ` lie strictly below `Yₜ`, or at least `|C| − k` of the lower ends
`μ̂₋ᵢ(Xₜ) − Rᵢ` strictly above it. Here `i` runs over the other units of `C`, `μ̂₋ᵢ` is the fit
leaving out `i` and `t`, and `Rᵢ` its absolute residual at `i`. That is `Yₜ` outside the closed
interval from the `k`-th smallest lower end to the `k`-th largest upper end, with
`k = ⌊α(n + 1)⌋` at Barber et al.'s level `α`. With a fit that reads only its training units, the
interval is computed from the training data and `Xₜ`. -/
def JackRefuses {Z X : Type*} (cov : Z → X) (resp : Z → ℝ) (k : ℕ) (C : Finset ι)
    (alg : Finset ι → (ι → Z) → X → ℝ) (z : ι → Z) (t : ι) : Prop :=
  C.card ≤ ((C.erase t).filter fun i => alg (univ \ {i, t}) z (cov (z t)) +
      |resp (z i) - alg (univ \ {i, t}) z (cov (z i))| < resp (z t)).card + k ∨
    C.card ≤ ((C.erase t).filter fun i => resp (z t) < alg (univ \ {i, t}) z (cov (z t)) -
      |resp (z i) - alg (univ \ {i, t}) z (cov (z i))|).card + k

noncomputable instance {Z X : Type*} (cov : Z → X) (resp : Z → ℝ) (k : ℕ) (C : Finset ι)
    (alg : Finset ι → (ι → Z) → X → ℝ) (z : ι → Z) (t : ι) :
    Decidable (JackRefuses cov resp k C alg z t) :=
  inferInstanceAs (Decidable (_ ∨ _))

/-- **The jackknife+ refuses only a strange test unit.** Each upper end below `Yₜ`, or lower end
above it, is a unit `i` whose residual under the fit leaving out `i` and `t` is below the test
unit's: a refusal by the interval is a refusal by the cross-fitted rule, so the jackknife+ interval
contains the cross-fitted set. -/
theorem not_acceptsPairwise_of_jackRefuses {Z X : Type*} {cov : Z → X} {resp : Z → ℝ} {k : ℕ}
    {C : Finset ι} {alg : Finset ι → (ι → Z) → X → ℝ} {z : ι → Z} {t : ι}
    (h : JackRefuses cov resp k C alg z t) :
    ¬ AcceptsPairwise k C (pairScore (jackFit cov resp alg) z) t := by
  unfold AcceptsPairwise
  rw [not_lt]
  have hpair : ∀ i, (univ \ {t, i} : Finset ι) = univ \ {i, t} := fun i => by rw [pair_comm]
  rcases h with h | h
  · refine h.trans (Nat.add_le_add_right (card_le_card fun i hi => ?_) k)
    rw [mem_filter, mem_erase] at hi
    rw [strangerThan, mem_filter]
    refine ⟨hi.1.2, ?_⟩
    simp only [pairScore, jackFit, hpair]
    rw [neg_lt_neg_iff]
    have h3 := le_abs_self (resp (z t) - alg (univ \ {i, t}) z (cov (z t)))
    linarith [hi.2]
  · refine h.trans (Nat.add_le_add_right (card_le_card fun i hi => ?_) k)
    rw [mem_filter, mem_erase] at hi
    rw [strangerThan, mem_filter]
    refine ⟨hi.1.2, ?_⟩
    simp only [pairScore, jackFit, hpair]
    rw [neg_lt_neg_iff]
    have h3 := le_abs_self (alg (univ \ {i, t}) z (cov (z t)) - resp (z t))
    rw [abs_sub_comm] at h3
    linarith [hi.2]

/-- **The jackknife+ over orderings.** With a symmetric fit, over every ordering of a fixed data
set, the jackknife+ covers the test unit at least `1 − (2k − 1)/|C|` of the time. -/
theorem jackknifePlus_orderings {Z X : Type*} {C : Finset ι} {t : ι} (ht : t ∈ C) (cov : Z → X)
    (resp : Z → ℝ) {alg : Finset ι → (ι → Z) → X → ℝ} (halg : FitSymmetric C alg) (z : ι → Z)
    (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / C.card ≤
      (((permsOn C).filter fun σ : Equiv.Perm ι =>
        ¬ JackRefuses cov resp k C alg (z ∘ σ) t).card : ℝ) / (permsOn C).card := by
  have hsym : FitSymmetric C (jackFit cov resp alg) := fun σ hσ T z => by
    funext p
    simp only [jackFit, halg σ hσ T z]
  refine (coverage_pairwise_orderings ht (pairScore_equivariant hsym) z k).trans ?_
  have hG : (0 : ℝ) ≤ (permsOn C).card := Nat.cast_nonneg _
  refine div_le_div_of_nonneg_right ?_ hG
  have hsub : ((permsOn C).filter fun σ : Equiv.Perm ι =>
      AcceptsPairwise k C (pairScore (jackFit cov resp alg) (z ∘ σ)) t) ⊆
      (permsOn C).filter fun σ : Equiv.Perm ι => ¬ JackRefuses cov resp k C alg (z ∘ σ) t := by
    intro σ hσ
    rw [mem_filter] at hσ ⊢
    exact ⟨hσ.1, fun href => not_acceptsPairwise_of_jackRefuses href hσ.2⟩
  exact_mod_cast card_le_card hsub

/-- **The jackknife+ under an exchangeable law** (Barber, Candès, Ramdas and Tibshirani, Theorem 1,
with its finite bound). Given the fits' randomness `w`, the units of `C` are exchangeable up to the
relabeling `act`, and the regression fit is symmetric up to it. Then the jackknife+ covers the test
unit with probability at least `1 − (2k − 1)/|C|`, whatever the fit. -/
theorem jackknifePlus_coverage {W Z X : Type*} [Fintype W] [Fintype Z] {C : Finset ι} {t : ι}
    (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (act : Equiv.Perm ι → W ≃ W) (hP : ∀ σ ∈ permsOn C, ∀ w z, P (act σ w) (z ∘ σ) = P w z)
    (cov : Z → X) (resp : Z → ℝ) (alg : W → Finset ι → (ι → Z) → X → ℝ)
    (halg : ∀ σ ∈ permsOn C, ∀ w T z, alg (act σ w) T (z ∘ σ) = alg w (T.image σ) z) (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / C.card ≤
      ∑ w, ∑ z, P w z * if ¬ JackRefuses cov resp k C (alg w) z t then 1 else 0 := by
  refine (crossFit_coverage ht P hP0 hP1 act hP (fun w => jackFit cov resp (alg w))
    (fun σ hσ w T z => by funext p; simp only [jackFit, halg σ hσ w T z]) k).trans ?_
  refine sum_le_sum fun w _ => sum_le_sum fun z _ => mul_le_mul_of_nonneg_left ?_ (hP0 w z)
  by_cases hacc : AcceptsPairwise k C (pairScore (jackFit cov resp (alg w)) z) t
  · have : ¬ JackRefuses cov resp k C (alg w) z t := fun h =>
      not_acceptsPairwise_of_jackRefuses h hacc
    simp [hacc, this]
  · simp only [hacc, if_false]
    split_ifs <;> norm_num

/-- **The jackknife+'s `1 − 2α`.** At `k ≤ α|C|`, `k = ⌊α(n + 1)⌋` in particular, the jackknife+
covers with probability at least `1 − 2α`. -/
theorem jackknifePlus_two_alpha {W Z X : Type*} [Fintype W] [Fintype Z] {C : Finset ι} {t : ι}
    (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (act : Equiv.Perm ι → W ≃ W) (hP : ∀ σ ∈ permsOn C, ∀ w z, P (act σ w) (z ∘ σ) = P w z)
    (cov : Z → X) (resp : Z → ℝ) (alg : W → Finset ι → (ι → Z) → X → ℝ)
    (halg : ∀ σ ∈ permsOn C, ∀ w T z, alg (act σ w) T (z ∘ σ) = alg w (T.image σ) z) {k : ℕ}
    {α : ℝ} (hk : (k : ℝ) ≤ α * C.card) :
    1 - 2 * α ≤ ∑ w, ∑ z, P w z * if ¬ JackRefuses cov resp k C (alg w) z t then 1 else 0 :=
  (one_sub_two_mul_le (card_pos.2 ⟨t, ht⟩) hk).trans
    (jackknifePlus_coverage ht P hP0 hP1 act hP cov resp alg halg k)

end CrossFitted

/-! ### Slice 9: the jackknife+ with a fit that reads its units' order -/

section OrderReading

variable {n : ℕ}

/-- A fit **reads its training units' order**: on the data reordered by `σ`, the fit to the units
`T` is the fit to the units `T.image σ` of the original whenever `σ` keeps the units of `T` in
their order. It may depend on that order, as a sampler sweeping its units in index order with a
fixed seed does, but not on where in the index its units sit. Every symmetric fit qualifies. -/
def OrderReading {Z β : Type*} (fit : Finset (Fin (n + 1)) → (Fin (n + 1) → Z) → β) : Prop :=
  ∀ σ : Equiv.Perm (Fin (n + 1)), ∀ (T : Finset (Fin (n + 1))) (z : Fin (n + 1) → Z),
    StrictMonoOn σ (↑T : Set (Fin (n + 1))) → fit T (z ∘ σ) = fit (T.image σ) z

/-- Moving unit `p` to the front, the others keeping their order: `(Fin.cycleRange p).symm` sends
`0` to `p` and is increasing on every other unit. -/
theorem strictMonoOn_cycleRange_symm (p : Fin (n + 1)) :
    StrictMonoOn (Fin.cycleRange p).symm ({0}ᶜ : Set (Fin (n + 1))) := by
  intro a ha b hb hab
  obtain ⟨a', rfl⟩ := Fin.exists_succ_eq.2 ha
  obtain ⟨b', rfl⟩ := Fin.exists_succ_eq.2 hb
  rw [Fin.cycleRange_symm_succ, Fin.cycleRange_symm_succ]
  exact Fin.strictMono_succAbove p (Fin.succ_lt_succ_iff.1 hab)

/-- With the test unit moved to the front, its comparisons are unit `p`'s. -/
theorem pairScore_cycleRange {Z β : Type*} {fit : Finset (Fin (n + 1)) → (Fin (n + 1) → Z) → Z → β}
    (h : OrderReading fit) (p : Fin (n + 1)) (z : Fin (n + 1) → Z) (j : Fin (n + 1)) :
    pairScore fit (z ∘ (Fin.cycleRange p).symm) 0 j =
        pairScore fit z p ((Fin.cycleRange p).symm j) ∧
      pairScore fit (z ∘ (Fin.cycleRange p).symm) j 0 =
        pairScore fit z ((Fin.cycleRange p).symm j) p := by
  set σ := (Fin.cycleRange p).symm with hσ
  have h0 : σ 0 = p := Fin.cycleRange_symm_zero p
  have hmono : ∀ T : Finset (Fin (n + 1)), (0 : Fin (n + 1)) ∉ T →
      StrictMonoOn σ (T : Set (Fin (n + 1))) := fun T hT =>
    (strictMonoOn_cycleRange_symm p).mono fun x hx => by
      simp only [Set.mem_compl_iff, Set.mem_singleton_iff]
      rintro rfl
      exact hT hx
  have h1 : (0 : Fin (n + 1)) ∉ (univ \ {0, j} : Finset (Fin (n + 1))) := by simp
  have h2 : (0 : Fin (n + 1)) ∉ (univ \ {j, 0} : Finset (Fin (n + 1))) := by simp
  constructor
  · simp only [pairScore, h σ _ z (hmono _ h1), image_sdiff_pair, h0, Function.comp_apply]
  · simp only [pairScore, h σ _ z (hmono _ h2), image_sdiff_pair, h0, Function.comp_apply]

theorem acceptsPairwise_cycleRange {Z α : Type*} [LinearOrder α]
    {fit : Finset (Fin (n + 1)) → (Fin (n + 1) → Z) → Z → α} (h : OrderReading fit) (k : ℕ)
    (p : Fin (n + 1)) (z : Fin (n + 1) → Z) :
    AcceptsPairwise k univ (pairScore fit (z ∘ (Fin.cycleRange p).symm)) 0 ↔
      AcceptsPairwise k univ (pairScore fit z) p := by
  have hc : (strangerThan univ (pairScore fit (z ∘ (Fin.cycleRange p).symm)) 0).card =
      (strangerThan univ (pairScore fit z) p).card := by
    apply card_equiv (Fin.cycleRange p).symm
    intro j
    simp only [strangerThan, mem_filter, mem_univ, true_and, (pairScore_cycleRange h p z j).1,
      (pairScore_cycleRange h p z j).2]
  simp only [AcceptsPairwise, hc]

/-- **The cross-fitted rule with a fit that reads its units' order.** Given the fits' randomness
`w`, the units are exchangeable, and each fit reads only its training units' order. Then the test
unit, unit `0`, is accepted with probability at least `1 − (2k − 1)/(n + 1)`: moving any unit to
the front keeps the others' order, so the test unit is strange exactly as often as any other. -/
theorem coverage_orderReading {W Z α : Type*} [Fintype W] [Fintype Z] [LinearOrder α]
    (P : W → (Fin (n + 1) → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (hP : ∀ σ : Equiv.Perm (Fin (n + 1)), ∀ w z, P w (z ∘ σ) = P w z)
    (fit : W → Finset (Fin (n + 1)) → (Fin (n + 1) → Z) → Z → α)
    (hfit : ∀ w, OrderReading (fit w)) (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / (n + 1) ≤
      ∑ w, ∑ z, P w z * if AcceptsPairwise k univ (pairScore (fit w) z) 0 then 1 else 0 := by
  have hN : (0 : ℝ) < n + 1 := by positivity
  have hmove : ∀ p : Fin (n + 1),
      ∑ w, ∑ z, P w z * (if AcceptsPairwise k univ (pairScore (fit w) z) 0 then (1 : ℝ) else 0) =
        ∑ w, ∑ z, P w z * (if AcceptsPairwise k univ (pairScore (fit w) z) p then 1 else 0) := by
    intro p
    refine sum_congr rfl fun w _ => ?_
    rw [← Equiv.sum_comp (reorder (Fin.cycleRange p).symm)]
    refine Fintype.sum_congr _ _ fun z => ?_
    rw [show reorder (Fin.cycleRange p).symm z = z ∘ (Fin.cycleRange p).symm from rfl, hP]
    simp only [acceptsPairwise_cycleRange (hfit w) k p z]
  have hratio : ∀ w z, 1 - ((2 * k - 1 : ℕ) : ℝ) / (n + 1) ≤
      (∑ p : Fin (n + 1),
        if AcceptsPairwise k univ (pairScore (fit w) z) p then (1 : ℝ) else 0) / (n + 1) := by
    intro w z
    have h := ratio_acceptsPairwise_ge (C := (univ : Finset (Fin (n + 1)))) ⟨0, mem_univ _⟩
      (pairScore (fit w) z) k
    rw [card_univ, Fintype.card_fin] at h
    push_cast at h
    rwa [sum_boole]
  calc 1 - ((2 * k - 1 : ℕ) : ℝ) / (n + 1)
      = ∑ w, ∑ z, (1 - ((2 * k - 1 : ℕ) : ℝ) / (n + 1)) * P w z := by
        simp only [← mul_sum, hP1, mul_one]
    _ ≤ ∑ w, ∑ z, P w z * ((∑ p : Fin (n + 1),
          if AcceptsPairwise k univ (pairScore (fit w) z) p then (1 : ℝ) else 0) / (n + 1)) :=
        sum_le_sum fun w _ => sum_le_sum fun z _ => by
          rw [mul_comm]
          exact mul_le_mul_of_nonneg_left (hratio w z) (hP0 w z)
    _ = (∑ p : Fin (n + 1), ∑ w, ∑ z,
          P w z * if AcceptsPairwise k univ (pairScore (fit w) z) p then 1 else 0) / (n + 1) := by
        simp only [mul_div_assoc', mul_sum, sum_div]
        rw [sum_congr rfl fun w _ => sum_comm, sum_comm]
    _ = (∑ _p : Fin (n + 1), ∑ w, ∑ z,
          P w z * if AcceptsPairwise k univ (pairScore (fit w) z) 0 then 1 else 0) / (n + 1) := by
        rw [sum_congr rfl fun p _ => (hmove p).symm]
    _ = _ := by
        rw [sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
        push_cast
        field_simp

/-- **The jackknife+ with a fit that reads its units' order.** The fit may read its training units
in index order, as a sampler sweeping them with a fixed seed does. Given its randomness, the units
are exchangeable. Then the jackknife+ covers the test unit with probability at least
`1 − (2k − 1)/(n + 1)`, the symmetric fit's bound. -/
theorem jackknifePlus_orderReading {W Z X : Type*} [Fintype W] [Fintype Z]
    (P : W → (Fin (n + 1) → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (hP : ∀ σ : Equiv.Perm (Fin (n + 1)), ∀ w z, P w (z ∘ σ) = P w z) (cov : Z → X)
    (resp : Z → ℝ) (alg : W → Finset (Fin (n + 1)) → (Fin (n + 1) → Z) → X → ℝ)
    (halg : ∀ w, OrderReading (alg w)) (k : ℕ) :
    1 - ((2 * k - 1 : ℕ) : ℝ) / (n + 1) ≤
      ∑ w, ∑ z, P w z * if ¬ JackRefuses cov resp k univ (alg w) z 0 then 1 else 0 := by
  have hjack : ∀ w, OrderReading (jackFit cov resp (alg w)) := fun w σ T z hmono => by
    funext p
    simp only [jackFit, halg w σ T z hmono]
  refine (coverage_orderReading P hP0 hP1 hP (fun w => jackFit cov resp (alg w)) hjack k).trans ?_
  refine sum_le_sum fun w _ => sum_le_sum fun z _ => mul_le_mul_of_nonneg_left ?_ (hP0 w z)
  by_cases hacc : AcceptsPairwise k univ (pairScore (jackFit cov resp (alg w)) z) 0
  · have : ¬ JackRefuses cov resp k univ (alg w) z 0 := fun h =>
      not_acceptsPairwise_of_jackRefuses h hacc
    simp [hacc, this]
  · simp only [hacc, if_false]
    split_ifs <;> norm_num

end OrderReading

/-! ### Slice 9: CV+ by merging the folds' ranks (Vovk; Barber et al., Theorem 4(a)) -/

section MergedFolds

/-- **The merging inequality.** If `K` ranks sum to at most `L`, then `L ≤ Σₖ (2L/K − rₖ)₊`. -/
theorem le_sum_max_of_sum_le {κ : Type*} [Fintype κ] (r : κ → ℝ) {L : ℝ}
    (hK : 0 < Fintype.card κ) (h : ∑ k, r k ≤ L) :
    L ≤ ∑ k, max (2 * L / Fintype.card κ - r k) 0 := by
  have hK' : (0 : ℝ) < Fintype.card κ := by exact_mod_cast hK
  have hc : ∑ _k : κ, 2 * L / Fintype.card κ = 2 * L := by
    rw [sum_const, card_univ, nsmul_eq_mul]
    field_simp
  calc L ≤ ∑ k, (2 * L / Fintype.card κ - r k) := by
        rw [sum_sub_distrib, hc]
        linarith
    _ ≤ _ := sum_le_sum fun k _ => le_max_left _ _

/-- Telescoping on `1, …, M`: for `1 ≤ r ≤ M + 1`, `h r` is `h (M + 1)` plus the drops
`h j − h (j + 1)` at every `j ≥ r`. -/
theorem layer_cake (h : ℕ → ℝ) {M r : ℕ} (hr1 : 1 ≤ r) (hrM : r ≤ M + 1) :
    h r = h (M + 1) + ∑ j ∈ Icc 1 M, (h j - h (j + 1)) * (if r ≤ j then 1 else 0) := by
  have hfilt : ∑ j ∈ Icc 1 M, (h j - h (j + 1)) * (if r ≤ j then 1 else 0) =
      ∑ j ∈ Ico r (M + 1), (h j - h (j + 1)) := by
    simp only [mul_ite, mul_one, mul_zero]
    rw [← sum_filter]
    congr 1
    ext j
    simp only [mem_filter, mem_Icc, mem_Ico]
    omega
  rw [hfilt, sum_Ico_eq_sum_range]
  have := sum_range_sub' (fun i => h (r + i)) (M + 1 - r)
  simp only [show ∀ i, r + (i + 1) = r + i + 1 from fun i => by omega] at this
  rw [this, show r + (M + 1 - r) = M + 1 by omega, add_zero]
  ring

/-- Abel's summation on `1, …, M`. -/
theorem sum_Icc_sub_mul (h : ℕ → ℝ) (M : ℕ) :
    ∑ j ∈ Icc 1 M, (h j - h (j + 1)) * j = ∑ j ∈ Icc 1 M, h j - M * h (M + 1) := by
  induction M with
  | zero => simp
  | succ M ih =>
    rw [sum_Icc_succ_top (by omega), sum_Icc_succ_top (by omega), ih]
    push_cast
    ring

/-- **The layer cake.** A rank super-uniform on `1, …, M`, `P(r ≤ j) ≤ j/M` for every `j`, keeps the
expectation of any antitone `h` at most its average over `1, …, M`. -/
theorem sum_antitone_le {W V : Type*} [Fintype W] [Fintype V] (P : W → V → ℝ)
    (hP1 : ∑ w, ∑ v, P w v = 1) {M : ℕ} (hM : 0 < M) (r : W → V → ℕ)
    (hr : ∀ w v, 1 ≤ r w v ∧ r w v ≤ M)
    (hsu : ∀ j, ∑ w, ∑ v, P w v * (if r w v ≤ j then 1 else 0) ≤ (j : ℝ) / M)
    {h : ℕ → ℝ} (hanti : Antitone h) :
    ∑ w, ∑ v, P w v * h (r w v) ≤ (∑ j ∈ Icc 1 M, h j) / M := by
  have hM' : (0 : ℝ) < M := by exact_mod_cast hM
  have hexp : ∀ w v, P w v * h (r w v) = P w v * h (M + 1) +
      ∑ j ∈ Icc 1 M, (h j - h (j + 1)) * (P w v * (if r w v ≤ j then 1 else 0)) := by
    intro w v
    rw [layer_cake h (M := M) (hr w v).1 (by have := (hr w v).2; omega), mul_add, mul_sum]
    congr 1
    exact sum_congr rfl fun j _ => by ring
  have hswap : ∑ w, ∑ v, ∑ j ∈ Icc 1 M, (h j - h (j + 1)) * (P w v * (if r w v ≤ j then 1 else 0)) =
      ∑ j ∈ Icc 1 M, (h j - h (j + 1)) * ∑ w, ∑ v, P w v * (if r w v ≤ j then 1 else 0) := by
    rw [sum_congr rfl fun w _ => sum_comm, sum_comm]
    exact sum_congr rfl fun j _ => by rw [mul_sum]; exact sum_congr rfl fun w _ => by rw [mul_sum]
  have hconst : ∑ w, ∑ v, P w v * h (M + 1) = h (M + 1) := by
    simp only [← sum_mul, hP1, one_mul]
  simp only [hexp, sum_add_distrib]
  rw [hconst, hswap]
  calc h (M + 1) + ∑ j ∈ Icc 1 M, (h j - h (j + 1)) *
        ∑ w, ∑ v, P w v * (if r w v ≤ j then 1 else 0)
      ≤ h (M + 1) + ∑ j ∈ Icc 1 M, (h j - h (j + 1)) * ((j : ℝ) / M) := by
        gcongr with j hj
        · exact sub_nonneg.2 (hanti (Nat.le_succ j))
        · exact hsu j
    _ = (∑ j ∈ Icc 1 M, h j) / M := by
        have : ∑ j ∈ Icc 1 M, (h j - h (j + 1)) * ((j : ℝ) / M) =
            (∑ j ∈ Icc 1 M, (h j - h (j + 1)) * j) / M := by
          rw [sum_div]
          exact sum_congr rfl fun j _ => by ring
        rw [this, sum_Icc_sub_mul]
        field_simp
        ring

/-- The drops of `(c − j)₊` sum to at most `c²/2`: a right Riemann sum of a decreasing function. -/
theorem sum_max_sub_le {c : ℝ} (hc : 0 ≤ c) (M : ℕ) :
    ∑ j ∈ Icc 1 M, max (c - j) 0 ≤ c ^ 2 / 2 := by
  have key : ∀ a : ℝ, max a 0 ≤ max (a + 1) 0 ^ 2 / 2 - max a 0 ^ 2 / 2 := by
    intro a
    rcases le_or_gt 0 a with ha | ha
    · rw [max_eq_left ha, max_eq_left (by linarith)]
      nlinarith
    · rw [max_eq_right ha.le]
      have : 0 ≤ max (a + 1) 0 ^ 2 / 2 := by positivity
      linarith
  have hind : ∀ M : ℕ, ∑ j ∈ Icc 1 M, max (c - j) 0 ≤
      max c 0 ^ 2 / 2 - max (c - M) 0 ^ 2 / 2 := by
    intro M
    induction M with
    | zero => simp
    | succ M ih =>
      rw [sum_Icc_succ_top (by omega)]
      have hk := key (c - ((M + 1 : ℕ) : ℝ))
      have he : c - ((M + 1 : ℕ) : ℝ) + 1 = c - M := by push_cast; ring
      rw [he] at hk
      linarith
  have h0 : 0 ≤ max (c - M) 0 ^ 2 / 2 := by positivity
  have := hind M
  rw [max_eq_left hc] at this
  linarith

/-- **The merged ranks** (Vovk and Wang's average of p-values, as a count). `K` ranks, each
super-uniform on `1, …, M` and dependent in any way, sum to at most `L` with probability at most
`2L/(KM)`. -/
theorem merged_le {W V κ : Type*} [Fintype W] [Fintype V] [Fintype κ] (P : W → V → ℝ)
    (hP0 : ∀ w v, 0 ≤ P w v) (hP1 : ∑ w, ∑ v, P w v = 1) {M : ℕ} (hM : 0 < M)
    (r : κ → W → V → ℕ) (hr : ∀ k w v, 1 ≤ r k w v ∧ r k w v ≤ M)
    (hsu : ∀ k j, ∑ w, ∑ v, P w v * (if r k w v ≤ j then 1 else 0) ≤ (j : ℝ) / M)
    (hK : 0 < Fintype.card κ) {L : ℕ} (hL : 0 < L) :
    ∑ w, ∑ v, P w v * (if ∑ k, r k w v ≤ L then 1 else 0) ≤
      2 * L / (Fintype.card κ * M) := by
  have hK' : (0 : ℝ) < Fintype.card κ := by exact_mod_cast hK
  have hL' : (0 : ℝ) < L := by exact_mod_cast hL
  have hM' : (0 : ℝ) < M := by exact_mod_cast hM
  set c : ℝ := 2 * L / Fintype.card κ with hcdef
  have hc : 0 ≤ c := by positivity
  have hpt : ∀ w v, (if ∑ k, r k w v ≤ L then (1 : ℝ) else 0) ≤
      (∑ k, max (c - r k w v) 0) / L := by
    intro w v
    split_ifs with h
    · rw [le_div_iff₀ hL', one_mul]
      have h' : ∑ k, (r k w v : ℝ) ≤ L := by exact_mod_cast h
      exact le_sum_max_of_sum_le (fun k => (r k w v : ℝ)) hK h'
    · exact div_nonneg (sum_nonneg fun k _ => le_max_right _ _) hL'.le
  have hfold : ∀ k, ∑ w, ∑ v, P w v * max (c - r k w v) 0 ≤ c ^ 2 / 2 / M := by
    intro k
    have hanti : Antitone fun j : ℕ => max (c - j) 0 := by
      intro i j hij
      have : (i : ℝ) ≤ j := by exact_mod_cast hij
      exact max_le_max (by linarith) le_rfl
    calc ∑ w, ∑ v, P w v * max (c - r k w v) 0 ≤ (∑ j ∈ Icc 1 M, max (c - j) 0) / M :=
          sum_antitone_le P hP1 hM (r k) (hr k) (hsu k) hanti
      _ ≤ c ^ 2 / 2 / M := div_le_div_of_nonneg_right (sum_max_sub_le hc M) hM'.le
  calc ∑ w, ∑ v, P w v * (if ∑ k, r k w v ≤ L then 1 else 0)
      ≤ ∑ w, ∑ v, P w v * ((∑ k, max (c - r k w v) 0) / L) :=
        sum_le_sum fun w _ => sum_le_sum fun v _ => mul_le_mul_of_nonneg_left (hpt w v) (hP0 w v)
    _ = (∑ k, ∑ w, ∑ v, P w v * max (c - r k w v) 0) / L := by
        simp only [mul_div_assoc', mul_sum, sum_div]
        rw [sum_congr rfl fun w _ => sum_comm, sum_comm]
    _ ≤ (∑ _k : κ, c ^ 2 / 2 / M) / L :=
        div_le_div_of_nonneg_right (sum_le_sum fun k _ => hfold k) hL'.le
    _ = 2 * L / (Fintype.card κ * M) := by
        rw [sum_const, card_univ, nsmul_eq_mul, hcdef]
        field_simp

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype ι] [DecidableEq ι] in
/-- A unit of `C` holds rank `1` at least and `|C|` at most. -/
theorem one_le_rankLe {α : Type*} [LinearOrder α] {C : Finset ι} (s : ι → α) {u : ι} (hu : u ∈ C) :
    1 ≤ rankLe C s u :=
  card_pos.2 ⟨u, mem_filter.2 ⟨hu, le_rfl⟩⟩

omit [Fintype ι] [DecidableEq ι] in
theorem rankLe_le_card {α : Type*} [LinearOrder α] (C : Finset ι) (s : ι → α) (u : ι) :
    rankLe C s u ≤ C.card :=
  card_le_card (filter_subset _ _)

/-- **Each fold's rank is super-uniform** (CB-1 in the fold). -/
theorem rank_superUniform {W Z α : Type*} [Fintype W] [Fintype Z] [LinearOrder α] {C : Finset ι}
    {t : ι} (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z)
    (hP1 : ∑ w, ∑ z, P w z = 1) (hP : ∀ σ ∈ permsOn C, ∀ w z, P w (z ∘ σ) = P w z)
    (S : W → (ι → Z) → ι → α) (hS : ∀ w, Equivariant C (S w)) (j : ℕ) :
    ∑ w, ∑ z, P w z * (if rankLe C (S w z) t ≤ j then 1 else 0) ≤ (j : ℝ) / C.card := by
  have h := coverage_law ht P hP0 hP1 hP S hS j
  have hsplit : ∀ w z, P w z * (if rankLe C (S w z) t ≤ j then (1 : ℝ) else 0) =
      P w z - P w z * (if Accepts j C (S w z) t then 1 else 0) := by
    intro w z
    by_cases ha : Accepts j C (S w z) t
    · have : ¬ rankLe C (S w z) t ≤ j := not_le.2 ha
      simp [ha, this]
    · have : rankLe C (S w z) t ≤ j := not_lt.1 ha
      simp [ha, this]
  simp only [hsplit, sum_sub_distrib, hP1]
  linarith

/-- **Cross-conformal and CV+ with `K` folds, by merging the folds' ranks** (Vovk; Barber, Candès,
Ramdas and Tibshirani, Theorem 4(a)). `C k` is fold `k` with the test unit `t`, `m + 1` units each.
Given `w`, the units of each `C k` are exchangeable, and fold `k`'s score is equivariant on `C k`:
a fit that does not read fold `k` or the test unit qualifies, symmetric or not. The rule refuses
when `Σₖ rₖ < ℓ + K`, `rₖ` the test unit's rank in `C k`; that has probability at most
`2(ℓ + K − 1)/(K(m + 1))`. -/
theorem crossFold_refused_le {W Z α κ : Type*} [Fintype W] [Fintype Z] [LinearOrder α] [Fintype κ]
    {C : κ → Finset ι} {t : ι} (ht : ∀ k, t ∈ C k) {m : ℕ} (hC : ∀ k, (C k).card = m + 1)
    (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (hP : ∀ k, ∀ σ ∈ permsOn (C k), ∀ w z, P w (z ∘ σ) = P w z)
    (S : κ → W → (ι → Z) → ι → α) (hS : ∀ k w, Equivariant (C k) (S k w)) (ℓ : ℕ)
    (hK : 0 < Fintype.card κ) :
    ∑ w, ∑ z, P w z * (if ∑ k, rankLe (C k) (S k w z) t < ℓ + Fintype.card κ then 1 else 0) ≤
      2 * ((ℓ : ℝ) + Fintype.card κ - 1) / (Fintype.card κ * (m + 1)) := by
  have hK' : (0 : ℝ) < Fintype.card κ := by exact_mod_cast hK
  have hr : ∀ k w z, 1 ≤ rankLe (C k) (S k w z) t ∧ rankLe (C k) (S k w z) t ≤ m + 1 :=
    fun k w z => ⟨one_le_rankLe _ (ht k), (rankLe_le_card _ _ _).trans_eq (hC k)⟩
  have hsu : ∀ k j, ∑ w, ∑ z, P w z * (if rankLe (C k) (S k w z) t ≤ j then 1 else 0) ≤
      (j : ℝ) / (m + 1 : ℕ) := fun k j => by
    have := rank_superUniform (ht k) P hP0 hP1 (hP k) (S k) (hS k) j
    rwa [hC k] at this
  rcases Nat.eq_zero_or_pos (ℓ + Fintype.card κ - 1) with hL | hL
  · -- the event is empty: every rank is at least 1
    have hnone : ∀ w z, ¬ ∑ k, rankLe (C k) (S k w z) t < ℓ + Fintype.card κ := by
      intro w z h
      have : Fintype.card κ ≤ ∑ k, rankLe (C k) (S k w z) t := by
        calc Fintype.card κ = ∑ _k : κ, 1 := by simp
          _ ≤ _ := sum_le_sum fun k _ => (hr k w z).1
      omega
    simp only [hnone, if_false, mul_zero, sum_const_zero]
    have h1 : (1 : ℝ) ≤ Fintype.card κ := by exact_mod_cast hK
    have h2 : (0 : ℝ) ≤ ℓ := Nat.cast_nonneg ℓ
    exact div_nonneg (by linarith) (by positivity)
  · have hev : ∀ w z, (∑ k, rankLe (C k) (S k w z) t < ℓ + Fintype.card κ) ↔
        ∑ k, rankLe (C k) (S k w z) t ≤ ℓ + Fintype.card κ - 1 := fun w z => by omega
    simp only [hev]
    have h := merged_le P hP0 hP1 (Nat.succ_pos m) (fun k w z => rankLe (C k) (S k w z) t)
      hr hsu hK hL
    have hcast : ((ℓ + Fintype.card κ - 1 : ℕ) : ℝ) = (ℓ : ℝ) + Fintype.card κ - 1 := by
      rw [Nat.cast_sub (by omega)]
      push_cast
      ring
    rw [hcast] at h
    push_cast at h
    exact h

/-- **CV+ refuses** the test unit `t` at level `ℓ` when at least `n + 1 − ℓ` of the upper ends
`μ̂₋ₖ(Xₜ) + Rᵢ` lie strictly below `Yₜ`, or at least `n + 1 − ℓ` of the lower ends `μ̂₋ₖ(Xₜ) − Rᵢ`
strictly above it. Here `i` runs over the training units of each fold `k`, the units of `C k` other
than `t`; `n` counts them; `μ̂₋ₖ = alg k` is fold `k`'s fit and `Rᵢ` its absolute residual at `i`.
With `ℓ = ⌊α(n + 1)⌋` this is Barber et al.'s CV+ interval, closed. -/
def CVRefuses {Z X κ : Type*} [Fintype κ] (cov : Z → X) (resp : Z → ℝ) (ℓ : ℕ) (C : κ → Finset ι)
    (alg : κ → (ι → Z) → X → ℝ) (z : ι → Z) (t : ι) : Prop :=
  (∑ k, ((C k).erase t).card) + 1 ≤ (∑ k, (((C k).erase t).filter fun i =>
      alg k z (cov (z t)) + |resp (z i) - alg k z (cov (z i))| < resp (z t)).card) + ℓ ∨
    (∑ k, ((C k).erase t).card) + 1 ≤ (∑ k, (((C k).erase t).filter fun i =>
      resp (z t) < alg k z (cov (z t)) - |resp (z i) - alg k z (cov (z i))|).card) + ℓ

noncomputable instance {Z X κ : Type*} [Fintype κ] (cov : Z → X) (resp : Z → ℝ) (ℓ : ℕ)
    (C : κ → Finset ι) (alg : κ → (ι → Z) → X → ℝ) (z : ι → Z) (t : ι) :
    Decidable (CVRefuses cov resp ℓ C alg z t) :=
  inferInstanceAs (Decidable (_ ∨ _))

omit [Fintype ι] in
/-- The rank in `C` of a unit of `C` is one more than the other units scoring at most its score. -/
theorem rankLe_eq_card_erase_add_one {α : Type*} [LinearOrder α] {C : Finset ι} (s : ι → α) {t : ι}
    (ht : t ∈ C) : rankLe C s t = ((C.erase t).filter fun i => s i ≤ s t).card + 1 := by
  have hmem : t ∈ C.filter fun i => s i ≤ s t := mem_filter.2 ⟨ht, le_rfl⟩
  rw [rankLe, filter_erase, card_erase_of_mem hmem]
  have : 0 < (C.filter fun i => s i ≤ s t).card := card_pos.2 ⟨t, hmem⟩
  omega

omit [Fintype ι] in
/-- **CV+ refuses only by the merged rule.** A unit `i` whose upper end lies below `Yₜ` (or lower
end above it) is one the test unit is stranger than under fold `k`'s fit, so it does not count
toward the test unit's rank: a refusal leaves the ranks summing below `ℓ + K`. -/
theorem sum_rankLe_lt_of_cvRefuses {Z X κ : Type*} [Fintype κ] {cov : Z → X} {resp : Z → ℝ}
    {ℓ : ℕ} {C : κ → Finset ι} {alg : κ → (ι → Z) → X → ℝ} {z : ι → Z} {t : ι}
    (ht : ∀ k, t ∈ C k) (h : CVRefuses cov resp ℓ C alg z t) :
    ∑ k, rankLe (C k) (fun i => -|resp (z i) - alg k z (cov (z i))|) t < ℓ + Fintype.card κ := by
  set s : κ → ι → ℝ := fun k i => -|resp (z i) - alg k z (cov (z i))| with hs
  have hrank : ∀ k,
      rankLe (C k) (s k) t = (((C k).erase t).filter fun i => s k i ≤ s k t).card + 1 :=
    fun k => rankLe_eq_card_erase_add_one (s k) (ht k)
  have hsum : ∑ k, rankLe (C k) (s k) t =
      ∑ k, (((C k).erase t).filter fun i => s k i ≤ s k t).card + Fintype.card κ := by
    rw [sum_congr rfl fun k _ => hrank k, sum_add_distrib]
    simp
  -- a unit beaten on either side does not score at most the test unit
  have hdisj : ∀ (k : κ) (p : ι → Prop) [DecidablePred p], (∀ i, p i → s k t < s k i) →
      (((C k).erase t).filter fun i => s k i ≤ s k t).card + (((C k).erase t).filter p).card ≤
        ((C k).erase t).card := by
    intro k p _ hp
    rw [← card_union_of_disjoint]
    · exact card_le_card (union_subset (filter_subset _ _) (filter_subset _ _))
    · rw [disjoint_filter]
      intro i _ h1 h2
      exact absurd h1 (not_le.2 (hp i h2))
  rcases h with h | h
  · have hb := sum_le_sum fun k (_ : k ∈ univ) => hdisj k
      (fun i => alg k z (cov (z t)) + |resp (z i) - alg k z (cov (z i))| < resp (z t)) (by
        intro i hi
        simp only [hs]
        rw [neg_lt_neg_iff]
        have := le_abs_self (resp (z t) - alg k z (cov (z t)))
        linarith)
    rw [sum_add_distrib] at hb
    rw [hsum]
    omega
  · have hb := sum_le_sum fun k (_ : k ∈ univ) => hdisj k
      (fun i => resp (z t) < alg k z (cov (z t)) - |resp (z i) - alg k z (cov (z i))|) (by
        intro i hi
        simp only [hs]
        rw [neg_lt_neg_iff]
        have := le_abs_self (alg k z (cov (z t)) - resp (z t))
        rw [abs_sub_comm] at this
        linarith)
    rw [sum_add_distrib] at hb
    rw [hsum]
    omega

/-- **CV+ by merging the folds** (Barber et al., Theorem 4(a), in counts). Given `w`, the units of
each fold with the test unit are exchangeable, and fold `k`'s fit `alg k w` does not change when
they are reordered: it reads neither, in any order. Then CV+ covers with probability at least
`1 − 2(ℓ + K − 1)/(K(m + 1))`. -/
theorem cvPlus_merged {W Z X κ : Type*} [Fintype W] [Fintype Z] [Fintype κ] {C : κ → Finset ι}
    {t : ι} (ht : ∀ k, t ∈ C k) {m : ℕ} (hC : ∀ k, (C k).card = m + 1) (P : W → (ι → Z) → ℝ)
    (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (hP : ∀ k, ∀ σ ∈ permsOn (C k), ∀ w z, P w (z ∘ σ) = P w z) (cov : Z → X) (resp : Z → ℝ)
    (alg : κ → W → (ι → Z) → X → ℝ)
    (halg : ∀ k, ∀ σ ∈ permsOn (C k), ∀ w z, alg k w (z ∘ σ) = alg k w z) (ℓ : ℕ)
    (hK : 0 < Fintype.card κ) :
    1 - 2 * ((ℓ : ℝ) + Fintype.card κ - 1) / (Fintype.card κ * (m + 1)) ≤
      ∑ w, ∑ z, P w z * if ¬ CVRefuses cov resp ℓ C (fun k => alg k w) z t then 1 else 0 := by
  set S : κ → W → (ι → Z) → ι → ℝ := fun k w z i => -|resp (z i) - alg k w z (cov (z i))|
  have hS : ∀ k w, Equivariant (C k) (S k w) := by
    intro k w σ hσ z i
    simp only [S, Function.comp_apply, halg k σ hσ w z]
  have href := crossFold_refused_le ht hC P hP0 hP1 hP S hS ℓ hK
  have hpt : ∀ w z,
      P w z * (if ¬ CVRefuses cov resp ℓ C (fun k => alg k w) z t then (1 : ℝ) else 0) ≥
      P w z - P w z * (if ∑ k, rankLe (C k) (S k w z) t < ℓ + Fintype.card κ then 1 else 0) := by
    intro w z
    by_cases hc : CVRefuses cov resp ℓ C (fun k => alg k w) z t
    · have := sum_rankLe_lt_of_cvRefuses (alg := fun k => alg k w) ht hc
      simp [hc, S, this]
    · have h0 := hP0 w z
      by_cases hr : ∑ k, rankLe (C k) (S k w z) t < ℓ + Fintype.card κ
      · simp [hc, hr, h0]
      · simp [hc, hr]
  calc 1 - 2 * ((ℓ : ℝ) + Fintype.card κ - 1) / (Fintype.card κ * (m + 1))
      ≤ ∑ w, ∑ z, (P w z - P w z *
          (if ∑ k, rankLe (C k) (S k w z) t < ℓ + Fintype.card κ then 1 else 0)) := by
        simp only [sum_sub_distrib, hP1]
        linarith
    _ ≤ _ := sum_le_sum fun w _ => sum_le_sum fun z _ => hpt w z

/-- **CV+ by merging, at level `α`.** At `ℓ ≤ α(n + 1)`, `n = Km`, CV+ covers with probability at
least `1 − 2α − 2(1 − α)(K − 1)/(n + K)`: Vovk's term, as Barber et al. state it. -/
theorem cvPlus_merged_alpha {K m ℓ : ℕ} {α : ℝ} (hK : 0 < K)
    (hℓ : (ℓ : ℝ) ≤ α * (K * m + 1)) :
    1 - 2 * α - 2 * (1 - α) * (K - 1) / (K * m + K) ≤
      1 - 2 * ((ℓ : ℝ) + K - 1) / (K * (m + 1)) := by
  have hK' : (0 : ℝ) < K := by exact_mod_cast hK
  have hpos : (0 : ℝ) < K * (m + 1) := by positivity
  have he : (K : ℝ) * m + K = K * (m + 1) := by ring
  rw [he]
  rw [show 1 - 2 * α - 2 * (1 - α) * (K - 1) / (K * (m + 1)) =
      1 - 2 * (α * (K * m + 1) + K - 1) / (K * (m + 1)) by field_simp; ring]
  have : 2 * ((ℓ : ℝ) + K - 1) / (K * (m + 1)) ≤ 2 * (α * (K * m + 1) + K - 1) / (K * (m + 1)) :=
    div_le_div_of_nonneg_right (by linarith) hpos.le
  linarith

/-- **CV+ by one fold.** A refusal leaves the test unit's rank at most `ℓ` in every fold, so in fold
`k₀` alone, and that fold with the test unit is a split conformal calibration set: CV+ covers with
probability at least `1 − ℓ/(m + 1)`. Only fold `k₀`'s exchangeability with the test unit and its
fit are needed. With two folds, `ℓ = ⌊α(2m + 1)⌋`, this is above `1 − 2α`. -/
theorem cvPlus_oneFold {W Z X κ : Type*} [Fintype W] [Fintype Z] [Fintype κ] {C : κ → Finset ι}
    {t : ι} (ht : ∀ k, t ∈ C k) {m : ℕ} (hC : ∀ k, (C k).card = m + 1) (P : W → (ι → Z) → ℝ)
    (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1) (k₀ : κ)
    (hP : ∀ σ ∈ permsOn (C k₀), ∀ w z, P w (z ∘ σ) = P w z) (cov : Z → X) (resp : Z → ℝ)
    (alg : κ → W → (ι → Z) → X → ℝ)
    (halg : ∀ σ ∈ permsOn (C k₀), ∀ w z, alg k₀ w (z ∘ σ) = alg k₀ w z) (ℓ : ℕ) :
    1 - (ℓ : ℝ) / (m + 1) ≤
      ∑ w, ∑ z, P w z * if ¬ CVRefuses cov resp ℓ C (fun k => alg k w) z t then 1 else 0 := by
  classical
  set S : κ → W → (ι → Z) → ι → ℝ := fun k w z i => -|resp (z i) - alg k w z (cov (z i))|
    with hSdef
  have hS : ∀ w, Equivariant (C k₀) (S k₀ w) := by
    intro w σ hσ z i
    simp only [hSdef, Function.comp_apply, halg σ hσ w z]
  have hcov := coverage_law (ht k₀) P hP0 hP1 (fun σ hσ w z => hP σ hσ w z) (S k₀) hS ℓ
  rw [hC k₀] at hcov
  push_cast at hcov
  refine hcov.trans (sum_le_sum fun w _ => sum_le_sum fun z _ =>
    mul_le_mul_of_nonneg_left ?_ (hP0 w z))
  by_cases hc : CVRefuses cov resp ℓ C (fun k => alg k w) z t
  · have hsum : ∑ k, rankLe (C k) (S k w z) t < ℓ + Fintype.card κ :=
      sum_rankLe_lt_of_cvRefuses (alg := fun k => alg k w) ht hc
    have hK : 0 < Fintype.card κ := Fintype.card_pos_iff.2 ⟨k₀⟩
    have h1 := add_sum_erase univ (fun k => rankLe (C k) (S k w z) t) (mem_univ k₀)
    have h2 : (univ.erase k₀).card ≤ ∑ k ∈ univ.erase k₀, rankLe (C k) (S k w z) t := by
      calc (univ.erase k₀).card = ∑ _k ∈ univ.erase k₀, 1 := by simp
        _ ≤ _ := sum_le_sum fun k _ => one_le_rankLe _ (ht k)
    rw [card_erase_of_mem (mem_univ _), card_univ] at h2
    have hle : rankLe (C k₀) (S k₀ w z) t ≤ ℓ := by omega
    have hna : ¬ Accepts ℓ (C k₀) (S k₀ w z) t := not_lt.2 hle
    simp [hc, hna]
  · simp only [hc, not_false_eq_true, if_true]
    split_ifs <;> norm_num

end MergedFolds

/-! ### Slice 9: CV+ by a tournament across folds (Barber et al., Theorem 4(b)) -/

section FoldTournament

variable {ι β : Type*} [DecidableEq β]

/-- The units of `C` in other folds than `u`'s that `u` is stranger than. -/
def strangerAcross {α : Type*} [LinearOrder α] (C : Finset ι) (f : ι → β) (s : ι → ι → α)
    (u : ι) : Finset ι :=
  C.filter fun j => f j ≠ f u ∧ s u j < s j u

/-- **The cross-fold rule at level `ℓ`**, with folds of `m` units: `u` is refused when it is
stranger than at least `n + 1 − ℓ` units of other folds, `n = |C| − m`. A tie makes neither
stranger. -/
def AcceptsAcross {α : Type*} [LinearOrder α] (m ℓ : ℕ) (C : Finset ι) (f : ι → β)
    (s : ι → ι → α) (u : ι) : Prop :=
  (strangerAcross C f s u).card + ℓ + m < C.card + 1

instance {α : Type*} [LinearOrder α] (m ℓ : ℕ) (C : Finset ι) (f : ι → β) (s : ι → ι → α)
    (u : ι) : Decidable (AcceptsAcross m ℓ C f s u) :=
  inferInstanceAs (Decidable (_ < _))

/-- **The counting core across folds** (Barber et al.'s Step 1 for CV+). Every fold of `C` holds
`m` units, and units play only across folds, ties allowed. Then at most `2ℓ + m − 2` units are
refused at level `ℓ`: each refused unit is stranger than at least `n + 1 − ℓ` others, and among the
refused units each comparison across folds goes at most one way. -/
theorem card_refused_across_le {α : Type*} [LinearOrder α] (C : Finset ι) (f : ι → β) {m : ℕ}
    (hm : ∀ u ∈ C, (C.filter fun j => f j = f u).card = m) (s : ι → ι → α) (ℓ : ℕ) :
    (C.filter fun u => ¬ AcceptsAcross m ℓ C f s u).card ≤ 2 * ℓ + m - 2 := by
  classical
  set F := C.filter fun u => ¬ AcceptsAcross m ℓ C f s u
  have hFC : F ⊆ C := filter_subset _ _
  set a : ι → ℕ := fun u => (F.filter fun j => f j ≠ f u ∧ s u j < s j u).card with ha
  set q : ι → ℕ := fun u => (F.filter fun j => f j = f u).card with hq
  set g : ι → ℕ := fun u => ((C \ F).filter fun j => f j ≠ f u).card with hg
  have hunit : ∀ u ∈ F, C.card + 1 ≤ a u + g u + ℓ + m := by
    intro u hu
    have href : C.card + 1 ≤ (strangerAcross C f s u).card + ℓ + m := by
      simpa only [AcceptsAcross, not_lt] using (mem_filter.1 hu).2
    have hsplit : (strangerAcross C f s u).card ≤ a u + g u := by
      calc (strangerAcross C f s u).card ≤ ((F.filter fun j => f j ≠ f u ∧ s u j < s j u) ∪
            ((C \ F).filter fun j => f j ≠ f u)).card := by
            refine card_le_card fun j hj => ?_
            rw [strangerAcross, mem_filter] at hj
            by_cases hjF : j ∈ F
            · exact mem_union_left _ (mem_filter.2 ⟨hjF, hj.2⟩)
            · exact mem_union_right _ (mem_filter.2 ⟨mem_sdiff.2 ⟨hj.1, hjF⟩, hj.2.1⟩)
        _ ≤ _ := card_union_le _ _
    omega
  have hfold : ∀ u ∈ F, g u + m + F.card = C.card + q u := by
    intro u hu
    have h1 := (C \ F).card_filter_add_card_filter_not (fun j => f j ≠ f u)
    have h2 : ((C \ F).filter fun j => ¬ f j ≠ f u) =
        (C.filter fun j => f j = f u) \ (F.filter fun j => f j = f u) := by
      ext j
      simp only [mem_filter, mem_sdiff, not_not]
      tauto
    have h3 : (F.filter fun j => f j = f u) ⊆ (C.filter fun j => f j = f u) :=
      filter_subset_filter _ hFC
    have h4 := card_le_card h3
    have h5 := card_le_card hFC
    rw [h2, card_sdiff_of_subset h3, card_sdiff_of_subset hFC] at h1
    rw [hm u (hFC hu)] at h1 h4
    simp only [hg, hq]
    omega
  have hqm : ∀ u ∈ F, q u ≤ m := fun u hu => by
    rw [← hm u (hFC hu)]
    exact card_le_card (filter_subset_filter _ hFC)
  have hle : ∀ u j, ((if f j ≠ f u ∧ s u j < s j u then 1 else 0) +
      (if f u ≠ f j ∧ s j u < s u j then 1 else 0) + (if f j = f u then 1 else 0) : ℕ) ≤ 1 := by
    intro u j
    by_cases hf : f j = f u
    · simp [hf]
    · have hf' : f u ≠ f j := fun h => hf h.symm
      by_cases h1 : s u j < s j u
      · have h2 : ¬ s j u < s u j := fun h2 => lt_asymm h1 h2
        simp [hf, hf', h1, h2]
      · simp only [hf, hf', h1, and_false, if_false, ne_eq, not_false_eq_true, true_and]
        split_ifs <;> omega
  have hA : ∑ u ∈ F, a u = ∑ u ∈ F, ∑ j ∈ F, if f j ≠ f u ∧ s u j < s j u then 1 else 0 := by
    simp only [ha, card_filter]
  have hB : ∑ u ∈ F, a u = ∑ u ∈ F, ∑ j ∈ F, if f u ≠ f j ∧ s j u < s u j then 1 else 0 := by
    rw [hA, sum_comm]
  have hQ : ∑ u ∈ F, q u = ∑ u ∈ F, ∑ j ∈ F, if f j = f u then 1 else 0 := by
    simp only [hq, card_filter]
  have hpairs : 2 * (∑ u ∈ F, a u) + ∑ u ∈ F, q u ≤ F.card * F.card := by
    calc 2 * (∑ u ∈ F, a u) + ∑ u ∈ F, q u
        = ∑ u ∈ F, ∑ j ∈ F, ((if f j ≠ f u ∧ s u j < s j u then 1 else 0) +
            (if f u ≠ f j ∧ s j u < s u j then 1 else 0) + (if f j = f u then 1 else 0)) := by
          rw [two_mul]
          nth_rw 1 [hA]
          rw [hB, hQ, ← sum_add_distrib, ← sum_add_distrib]
          refine sum_congr rfl fun u _ => ?_
          rw [← sum_add_distrib, ← sum_add_distrib]
      _ ≤ ∑ u ∈ F, ∑ j ∈ F, 1 := sum_le_sum fun u _ => sum_le_sum fun j _ => hle u j
      _ = F.card * F.card := by simp
  have hs1 : F.card * (C.card + 1) ≤
      (∑ u ∈ F, a u) + (∑ u ∈ F, g u) + F.card * ℓ + F.card * m := by
    have h := sum_le_sum hunit
    simp only [sum_const, smul_eq_mul, sum_add_distrib] at h
    linarith
  have hs3 : (∑ u ∈ F, g u) + F.card * m + F.card * F.card =
      F.card * C.card + ∑ u ∈ F, q u := by
    have h := sum_congr rfl hfold
    simp only [sum_add_distrib, sum_const, smul_eq_mul] at h
    linarith
  have hs5 : ∑ u ∈ F, q u ≤ F.card * m := by
    have h := sum_le_sum hqm
    simp only [sum_const, smul_eq_mul] at h
    linarith
  have h2 : F.card * (F.card + 2) ≤ F.card * (m + 2 * ℓ) := by nlinarith
  rcases Nat.eq_zero_or_pos F.card with h0 | hpos
  · rw [h0]
    exact Nat.zero_le _
  · have := Nat.le_of_mul_le_mul_left h2 hpos
    omega

omit [DecidableEq β] in
/-- Averaging over the permutations of `β`: each point is `π a` equally often. -/
theorem sum_perm_apply [Fintype β] [DecidableEq β] (a : β) (g : β → ℝ) :
    (∑ π : Equiv.Perm β, g (π a)) / Fintype.card (Equiv.Perm β) =
      (∑ b, g b) / Fintype.card β := by
  have hG : (0 : ℝ) < Fintype.card (Equiv.Perm β) := by exact_mod_cast Fintype.card_pos
  have hB : (0 : ℝ) < Fintype.card β := by exact_mod_cast Fintype.card_pos_iff.2 ⟨a⟩
  have huniv : permsOn (univ : Finset β) = univ := eq_univ_of_forall mem_permsOn_univ
  have hfib : ∀ b, (((univ : Finset (Equiv.Perm β)).filter fun π => π a = b).card : ℝ) =
      Fintype.card (Equiv.Perm β) / Fintype.card β := by
    intro b
    have h := ratio_permsOn (C := (univ : Finset β)) (mem_univ a) (fun x => x = b)
    rw [huniv, card_univ, card_univ] at h
    have h1 : ((univ : Finset β).filter fun x => x = b).card = 1 := by
      rw [filter_eq' univ b]
      simp
    rw [h1, Nat.cast_one, div_eq_div_iff hG.ne' hB.ne'] at h
    field_simp
    linarith
  have hrw : ∀ π : Equiv.Perm β, g (π a) = ∑ b, g b * if π a = b then 1 else 0 := by
    intro π
    simp [mul_ite]
  calc (∑ π : Equiv.Perm β, g (π a)) / Fintype.card (Equiv.Perm β)
      = (∑ b, g b * (((univ : Finset (Equiv.Perm β)).filter fun π => π a = b).card : ℝ)) /
          Fintype.card (Equiv.Perm β) := by
        simp only [hrw]
        rw [sum_comm]
        congr 1
        refine sum_congr rfl fun b _ => ?_
        rw [← mul_sum, sum_boole]
    _ = (∑ b, g b) / Fintype.card β := by
        simp only [hfib]
        rw [← sum_mul]
        field_simp

section Product

variable {γ : Type*} [Fintype β] [Fintype γ] [DecidableEq γ]

/-- Averaging over the permutations of the folds and of the slots within them,
`(b, c) ↦ (π b, ρ c)`: each unit is the image of `t` equally often. -/
theorem ratio_prodPerm (t : β × γ) (q : β × γ → Prop) [DecidablePred q] :
    (∑ p : Equiv.Perm β × Equiv.Perm γ,
        if q (Equiv.prodCongr p.1 p.2 t) then (1 : ℝ) else 0) /
        Fintype.card (Equiv.Perm β × Equiv.Perm γ) =
      ((univ.filter q).card : ℝ) / (Fintype.card β * Fintype.card γ) := by
  set Y : β → γ → ℝ := fun b c => if q (b, c) then 1 else 0 with hY
  have step1 : ∀ π : Equiv.Perm β, (∑ ρ : Equiv.Perm γ, Y (π t.1) (ρ t.2)) /
      Fintype.card (Equiv.Perm γ) = (∑ c, Y (π t.1) c) / Fintype.card γ :=
    fun π => sum_perm_apply t.2 (fun c => Y (π t.1) c)
  have step2 := sum_perm_apply t.1 (fun b => (∑ c, Y b c) / Fintype.card γ)
  have hcount : ∑ b, ∑ c, Y b c = ((univ.filter q).card : ℝ) := by
    rw [← sum_boole, ← Fintype.sum_prod_type']
  rw [Fintype.card_prod, Fintype.sum_prod_type]
  push_cast
  calc (∑ π : Equiv.Perm β, ∑ ρ : Equiv.Perm γ,
          if q (Equiv.prodCongr (π, ρ).1 (π, ρ).2 t) then (1 : ℝ) else 0) /
          (Fintype.card (Equiv.Perm β) * Fintype.card (Equiv.Perm γ))
      = (∑ π : Equiv.Perm β, (∑ ρ : Equiv.Perm γ, Y (π t.1) (ρ t.2)) /
          Fintype.card (Equiv.Perm γ)) / Fintype.card (Equiv.Perm β) := by
        rw [← sum_div, div_div, mul_comm]
        rfl
    _ = (∑ π : Equiv.Perm β, (∑ c, Y (π t.1) c) / Fintype.card γ) /
          Fintype.card (Equiv.Perm β) := by simp only [step1]
    _ = (∑ b, (∑ c, Y b c) / Fintype.card γ) / Fintype.card β := step2
    _ = ((univ.filter q).card : ℝ) / (Fintype.card β * Fintype.card γ) := by
        rw [← sum_div, div_div, mul_comm, hcount]

omit [DecidableEq γ] in
/-- Each fold of `β × γ` holds `|γ|` units. -/
theorem card_fold (u : β × γ) :
    ((univ : Finset (β × γ)).filter fun j => j.1 = u.1).card = Fintype.card γ := by
  rw [show ((univ : Finset (β × γ)).filter fun j => j.1 = u.1) = {u.1} ×ˢ univ by
    ext j
    simp only [mem_filter, mem_univ, true_and, mem_product, mem_singleton, and_true]]
  simp

omit [DecidableEq γ] in
theorem acceptsAcross_prodCongr {α : Type*} [LinearOrder α] {m ℓ : ℕ} (π : Equiv.Perm β)
    (ρ : Equiv.Perm γ) (s : β × γ → β × γ → α) (u : β × γ) :
    AcceptsAcross m ℓ univ Prod.fst
        (fun i j => s (Equiv.prodCongr π ρ i) (Equiv.prodCongr π ρ j)) u ↔
      AcceptsAcross m ℓ univ Prod.fst s (Equiv.prodCongr π ρ u) := by
  have h : (strangerAcross univ Prod.fst
        (fun i j => s (Equiv.prodCongr π ρ i) (Equiv.prodCongr π ρ j)) u).card =
      (strangerAcross univ Prod.fst s (Equiv.prodCongr π ρ u)).card := by
    apply card_equiv (Equiv.prodCongr π ρ)
    intro j
    simp only [strangerAcross, mem_filter, mem_univ, true_and, Equiv.prodCongr_apply,
      Prod.map_fst, ne_eq, π.injective.eq_iff]
  simp only [AcceptsAcross, h]

omit [DecidableEq γ] in
/-- The share of units the cross-fold rule accepts: at least `1 − (2ℓ + m − 2)/N`. -/
theorem ratio_acceptsAcross_ge [Nonempty β] [Nonempty γ] {α : Type*} [LinearOrder α]
    (s : β × γ → β × γ → α) (ℓ : ℕ) :
    1 - ((2 * ℓ + Fintype.card γ - 2 : ℕ) : ℝ) / (Fintype.card β * Fintype.card γ) ≤
      ((univ.filter fun u => AcceptsAcross (Fintype.card γ) ℓ univ Prod.fst s u).card : ℝ) /
        (Fintype.card β * Fintype.card γ) := by
  have hN : (0 : ℝ) < Fintype.card β * Fintype.card γ := by positivity
  have hrej := card_refused_across_le (univ : Finset (β × γ)) Prod.fst
    (fun u _ => card_fold u) s ℓ
  have hsplit := (univ : Finset (β × γ)).card_filter_add_card_filter_not
    (fun u => AcceptsAcross (Fintype.card γ) ℓ univ Prod.fst s u)
  rw [card_univ, Fintype.card_prod] at hsplit
  have hnat : Fintype.card β * Fintype.card γ ≤
      (univ.filter fun u => AcceptsAcross (Fintype.card γ) ℓ univ Prod.fst s u).card +
        (2 * ℓ + Fintype.card γ - 2) := by omega
  have hcount : ((Fintype.card β * Fintype.card γ : ℕ) : ℝ) ≤
      ((univ.filter fun u => AcceptsAcross (Fintype.card γ) ℓ univ Prod.fst s u).card : ℝ) +
        ((2 * ℓ + Fintype.card γ - 2 : ℕ) : ℝ) := by exact_mod_cast hnat
  push_cast at hcount
  rw [one_sub_div hN.ne']
  exact div_le_div_of_nonneg_right (by linarith) hN.le

/-- **The cross-fold rule under an exchangeable law.** Units are `(fold, slot)` in `β × γ`; the test
unit's fold is made up by further units from the same law, which the interval never reads. Given the
fits' randomness up to its relabeling `act`, the law is invariant under permuting the folds and the
slots, and the matrix equivariant. Then the test unit is accepted with probability at least
`1 − (2ℓ + m − 2)/(n + m)`, `m = |γ|` units a fold and `n + m = |β||γ|`. -/
theorem coverage_across [Nonempty β] [Nonempty γ] {W Z α : Type*} [Fintype W] [Fintype Z]
    [LinearOrder α] (t : β × γ) (P : W → (β × γ → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z)
    (hP1 : ∑ w, ∑ z, P w z = 1) (act : Equiv.Perm β × Equiv.Perm γ → W ≃ W)
    (hP : ∀ p w z, P (act p w) (z ∘ Equiv.prodCongr p.1 p.2) = P w z)
    (S : W → (β × γ → Z) → β × γ → β × γ → α)
    (hS : ∀ p w z i j, S (act p w) (z ∘ Equiv.prodCongr p.1 p.2) i j =
      S w z (Equiv.prodCongr p.1 p.2 i) (Equiv.prodCongr p.1 p.2 j)) (ℓ : ℕ) :
    1 - ((2 * ℓ + Fintype.card γ - 2 : ℕ) : ℝ) / (Fintype.card β * Fintype.card γ) ≤
      ∑ w, ∑ z, P w z *
        if AcceptsAcross (Fintype.card γ) ℓ univ Prod.fst (S w z) t then 1 else 0 := by
  set Acc : W × (β × γ → Z) → β × γ → Prop := fun x u =>
    AcceptsAcross (Fintype.card γ) ℓ univ Prod.fst (S x.1 x.2) u with hAcc
  rw [← Fintype.sum_prod_type'] at hP1 ⊢
  have hG : ((Fintype.card (Equiv.Perm β × Equiv.Perm γ) : ℕ) : ℝ) ≠ 0 := by
    exact_mod_cast Fintype.card_pos.ne'
  have hswap : ∀ p : Equiv.Perm β × Equiv.Perm γ,
      ∑ x : W × (β × γ → Z), P x.1 x.2 * (if Acc x t then (1 : ℝ) else 0) =
        ∑ x : W × (β × γ → Z),
          P x.1 x.2 * (if Acc x (Equiv.prodCongr p.1 p.2 t) then 1 else 0) := by
    intro p
    rw [← Equiv.sum_comp ((act p).prodCongr (reorder (Equiv.prodCongr p.1 p.2)))]
    refine Fintype.sum_congr _ _ fun x => ?_
    simp only [Equiv.prodCongr_apply, Prod.map_fst, Prod.map_snd, hAcc]
    rw [show reorder (Equiv.prodCongr p.1 p.2) x.2 = x.2 ∘ Equiv.prodCongr p.1 p.2 from rfl,
      hP p]
    have : S (act p x.1) (x.2 ∘ Equiv.prodCongr p.1 p.2) =
        fun i j => S x.1 x.2 (Equiv.prodCongr p.1 p.2 i) (Equiv.prodCongr p.1 p.2 j) :=
      funext fun i => funext fun j => hS p x.1 x.2 i j
    rw [this]
    simp only [acceptsAcross_prodCongr]
    rfl
  have key : ∀ x : W × (β × γ → Z),
      (1 - ((2 * ℓ + Fintype.card γ - 2 : ℕ) : ℝ) / (Fintype.card β * Fintype.card γ)) *
        P x.1 x.2 ≤ P x.1 x.2 * ((∑ p : Equiv.Perm β × Equiv.Perm γ,
          if Acc x (Equiv.prodCongr p.1 p.2 t) then (1 : ℝ) else 0) /
            Fintype.card (Equiv.Perm β × Equiv.Perm γ)) := by
    intro x
    rw [mul_comm]
    refine mul_le_mul_of_nonneg_left ?_ (hP0 x.1 x.2)
    rw [ratio_prodPerm t (Acc x)]
    exact ratio_acceptsAcross_ge (S x.1 x.2) ℓ
  calc 1 - ((2 * ℓ + Fintype.card γ - 2 : ℕ) : ℝ) / (Fintype.card β * Fintype.card γ)
      = ∑ x : W × (β × γ → Z),
          (1 - ((2 * ℓ + Fintype.card γ - 2 : ℕ) : ℝ) / (Fintype.card β * Fintype.card γ)) *
            P x.1 x.2 := by
        rw [← mul_sum, hP1, mul_one]
    _ ≤ ∑ x : W × (β × γ → Z), P x.1 x.2 * ((∑ p : Equiv.Perm β × Equiv.Perm γ,
          if Acc x (Equiv.prodCongr p.1 p.2 t) then (1 : ℝ) else 0) /
            Fintype.card (Equiv.Perm β × Equiv.Perm γ)) :=
        sum_le_sum fun x _ => key x
    _ = (∑ p : Equiv.Perm β × Equiv.Perm γ, ∑ x : W × (β × γ → Z),
          P x.1 x.2 * if Acc x (Equiv.prodCongr p.1 p.2 t) then 1 else 0) /
            Fintype.card (Equiv.Perm β × Equiv.Perm γ) := by
        simp only [mul_div_assoc', mul_sum, sum_div]
        exact sum_comm
    _ = (∑ _p : Equiv.Perm β × Equiv.Perm γ, ∑ x : W × (β × γ → Z),
          P x.1 x.2 * if Acc x t then 1 else 0) /
            Fintype.card (Equiv.Perm β × Equiv.Perm γ) := by
        rw [sum_congr rfl fun p _ => (hswap p).symm]
    _ = _ := by
        rw [sum_const, card_univ, nsmul_eq_mul]
        field_simp
        rfl

/-- The units outside the folds of `i` and `j`: the training set of the fit that leaves both
folds out. -/
def foldsOut (i j : β × γ) : Finset (β × γ) := univ.filter fun u => u.1 ≠ i.1 ∧ u.1 ≠ j.1

omit [DecidableEq γ] in
theorem foldsOut_comm (i j : β × γ) : foldsOut i j = foldsOut j i := by
  unfold foldsOut
  exact filter_congr fun u _ => and_comm

theorem image_foldsOut (π : Equiv.Perm β) (ρ : Equiv.Perm γ) (i j : β × γ) :
    (foldsOut i j).image (Equiv.prodCongr π ρ) =
      foldsOut (Equiv.prodCongr π ρ i) (Equiv.prodCongr π ρ j) := by
  ext x
  simp only [foldsOut, mem_image, mem_filter, mem_univ, true_and, Equiv.prodCongr_apply,
    Prod.map_fst]
  constructor
  · rintro ⟨a, ⟨h1, h2⟩, rfl⟩
    exact ⟨fun h => h1 (π.injective h), fun h => h2 (π.injective h)⟩
  · rintro ⟨h1, h2⟩
    refine ⟨(Equiv.prodCongr π ρ).symm x, ⟨fun h => h1 ?_, fun h => h2 ?_⟩,
      Prod.ext (π.apply_symm_apply x.1) (ρ.apply_symm_apply x.2)⟩
    · rw [← h]
      simp
    · rw [← h]
      simp

/-- **CV+ refuses**, folds as `β`, the test unit `t`'s fold made up by further units: at least
`n + 1 − ℓ` upper ends `μ̂₋ₖ(Xₜ) + Rᵢ` strictly below `Yₜ`, or as many lower ends strictly above,
`i` over the units of other folds than `t`'s, `μ̂₋ₖ` the fit to every unit outside the folds of `i`
and `t`. That fit reads only the `K` training folds, so the interval is CV+'s. -/
def CVRefusesAcross {Z X : Type*} (cov : Z → X) (resp : Z → ℝ) (ℓ : ℕ)
    (alg : Finset (β × γ) → (β × γ → Z) → X → ℝ) (z : β × γ → Z) (t : β × γ) : Prop :=
  Fintype.card β * Fintype.card γ + 1 ≤ ((univ.filter fun i : β × γ => i.1 ≠ t.1 ∧
      alg (foldsOut i t) z (cov (z t)) + |resp (z i) - alg (foldsOut i t) z (cov (z i))| <
        resp (z t)).card) + ℓ + Fintype.card γ ∨
    Fintype.card β * Fintype.card γ + 1 ≤ ((univ.filter fun i : β × γ => i.1 ≠ t.1 ∧
      resp (z t) < alg (foldsOut i t) z (cov (z t)) -
        |resp (z i) - alg (foldsOut i t) z (cov (z i))|).card) + ℓ + Fintype.card γ

noncomputable instance {Z X : Type*} (cov : Z → X) (resp : Z → ℝ) (ℓ : ℕ)
    (alg : Finset (β × γ) → (β × γ → Z) → X → ℝ) (z : β × γ → Z) (t : β × γ) :
    Decidable (CVRefusesAcross cov resp ℓ alg z t) :=
  inferInstanceAs (Decidable (_ ∨ _))

/-- The cross-fold score: unit `i` read by the fit that leaves out the folds of `i` and `j`. -/
def foldScore {Z X : Type*} (cov : Z → X) (resp : Z → ℝ)
    (alg : Finset (β × γ) → (β × γ → Z) → X → ℝ) (z : β × γ → Z) (i j : β × γ) : ℝ :=
  -|resp (z i) - alg (foldsOut i j) z (cov (z i))|

omit [DecidableEq γ] in
/-- **CV+ refuses only a strange test unit** across folds. -/
theorem not_acceptsAcross_of_cvRefuses {Z X : Type*} {cov : Z → X} {resp : Z → ℝ} {ℓ : ℕ}
    {alg : Finset (β × γ) → (β × γ → Z) → X → ℝ} {z : β × γ → Z} {t : β × γ}
    (h : CVRefusesAcross cov resp ℓ alg z t) :
    ¬ AcceptsAcross (Fintype.card γ) ℓ univ Prod.fst (foldScore cov resp alg z) t := by
  unfold AcceptsAcross
  rw [not_lt, card_univ, Fintype.card_prod]
  rcases h with h | h
  · refine h.trans (Nat.add_le_add_right (Nat.add_le_add_right
      (card_le_card fun i hi => ?_) ℓ) _)
    rw [mem_filter] at hi
    rw [strangerAcross, mem_filter]
    refine ⟨mem_univ _, hi.2.1, ?_⟩
    simp only [foldScore, foldsOut_comm t i]
    rw [neg_lt_neg_iff]
    have h3 := le_abs_self (resp (z t) - alg (foldsOut i t) z (cov (z t)))
    linarith [hi.2.2]
  · refine h.trans (Nat.add_le_add_right (Nat.add_le_add_right
      (card_le_card fun i hi => ?_) ℓ) _)
    rw [mem_filter] at hi
    rw [strangerAcross, mem_filter]
    refine ⟨mem_univ _, hi.2.1, ?_⟩
    simp only [foldScore, foldsOut_comm t i]
    rw [neg_lt_neg_iff]
    have h3 := le_abs_self (alg (foldsOut i t) z (cov (z t)) - resp (z t))
    rw [abs_sub_comm] at h3
    linarith [hi.2.2]

/-- **CV+ by the tournament across folds** (Barber, Candès, Ramdas and Tibshirani, Theorem 4(b),
in counts). The `n + 1` units come with `m − 1` further units from the same law, filling the test
unit's fold; the law of all `n + m` is invariant under permuting folds and slots, given the fit's
randomness up to its relabeling; the fit reads its training units as a set, up to that relabeling.
Then CV+ covers with probability at least `1 − (2ℓ + m − 2)/(n + m)`. -/
theorem cvPlus_tournament [Nonempty β] [Nonempty γ] {W Z X : Type*} [Fintype W] [Fintype Z]
    (t : β × γ) (P : W → (β × γ → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z)
    (hP1 : ∑ w, ∑ z, P w z = 1) (act : Equiv.Perm β × Equiv.Perm γ → W ≃ W)
    (hP : ∀ p w z, P (act p w) (z ∘ Equiv.prodCongr p.1 p.2) = P w z) (cov : Z → X)
    (resp : Z → ℝ) (alg : W → Finset (β × γ) → (β × γ → Z) → X → ℝ)
    (halg : ∀ p w T z, alg (act p w) T (z ∘ Equiv.prodCongr p.1 p.2) =
      alg w (T.image (Equiv.prodCongr p.1 p.2)) z) (ℓ : ℕ) :
    1 - ((2 * ℓ + Fintype.card γ - 2 : ℕ) : ℝ) / (Fintype.card β * Fintype.card γ) ≤
      ∑ w, ∑ z, P w z * if ¬ CVRefusesAcross cov resp ℓ (alg w) z t then 1 else 0 := by
  have hS : ∀ p w z i j,
      foldScore cov resp (alg (act p w)) (z ∘ Equiv.prodCongr p.1 p.2) i j =
        foldScore cov resp (alg w) z (Equiv.prodCongr p.1 p.2 i) (Equiv.prodCongr p.1 p.2 j) := by
    intro p w z i j
    simp only [foldScore, halg, image_foldsOut, Function.comp_apply]
  refine (coverage_across t P hP0 hP1 act hP (fun w => foldScore cov resp (alg w)) hS ℓ).trans ?_
  refine sum_le_sum fun w _ => sum_le_sum fun z _ => mul_le_mul_of_nonneg_left ?_ (hP0 w z)
  by_cases hacc : AcceptsAcross (Fintype.card γ) ℓ univ Prod.fst (foldScore cov resp (alg w) z) t
  · have : ¬ CVRefusesAcross cov resp ℓ (alg w) z t := fun h =>
      not_acceptsAcross_of_cvRefuses h hacc
    simp [hacc, this]
  · simp only [hacc, if_false]
    split_ifs <;> norm_num

end Product

end FoldTournament

/-! ### Slice 9: asymmetric CQR, a threshold per tail -/

section AsymmetricCQR

/-- CQR's lower-tail conformity `y − q_lo(x)`, the negative of Romano, Patterson and Candès's
`E_lo`. -/
def cqrLoScore {X : Type*} (qlo : X → ℝ) (p : X × ℝ) : ℝ := p.2 - qlo p.1

/-- CQR's upper-tail conformity `q_hi(x) − y`, the negative of their `E_hi`. -/
def cqrHiScore {X : Type*} (qhi : X → ℝ) (p : X × ℝ) : ℝ := qhi p.1 - p.2

/-- At level `0` the accept rule accepts every unit of `C`: a tail at level `0` has an infinite
end. -/
theorem accepts_zero {ι α : Type*} [LinearOrder α] {C : Finset ι} (s : ι → α) {u : ι}
    (hu : u ∈ C) : Accepts 0 C s u :=
  card_pos.2 ⟨u, mem_filter.2 ⟨hu, le_rfl⟩⟩

/-- **Two events at once.** Their joint probability is at least the sum of theirs, less one. -/
theorem sum_and_ge {W V : Type*} [Fintype W] [Fintype V] (P : W → V → ℝ) (hP0 : ∀ w v, 0 ≤ P w v)
    (hP1 : ∑ w, ∑ v, P w v = 1) (A B : W → V → Prop) [∀ w v, Decidable (A w v)]
    [∀ w v, Decidable (B w v)] :
    (∑ w, ∑ v, P w v * if A w v then 1 else 0) + (∑ w, ∑ v, P w v * if B w v then 1 else 0) - 1 ≤
      ∑ w, ∑ v, P w v * if A w v ∧ B w v then 1 else 0 := by
  have key : ∀ w v, P w v * (if A w v then 1 else 0) + P w v * (if B w v then 1 else 0) - P w v ≤
      P w v * (if A w v ∧ B w v then 1 else 0) := by
    intro w v
    have h := hP0 w v
    by_cases ha : A w v <;> by_cases hb : B w v
    · rw [if_pos ha, if_pos hb, if_pos ⟨ha, hb⟩]
      linarith
    · rw [if_pos ha, if_neg hb, if_neg fun h => hb h.2]
      linarith
    · rw [if_neg ha, if_pos hb, if_neg fun h => ha h.1]
      linarith
    · rw [if_neg ha, if_neg hb, if_neg fun h => ha h.1]
      linarith
  have hsum := sum_le_sum fun w (_ : w ∈ univ) => sum_le_sum fun v (_ : v ∈ univ) => key w v
  simp only [sum_sub_distrib, sum_add_distrib] at hsum
  linarith

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- **Asymmetric CQR, each tail and both** (Romano, Patterson and Candès, Theorem 2). `W` is the
training split with anything else the quantile fits read, under any law; given `w`, the calibration
and test units are exchangeable, and `pt` reads a unit as `(x, y)`. Each tail is the accept rule on
its own conformity at its own level. The test unit falls below the set with probability at most
`k_lo/|C|`, above it with probability at most `k_hi/|C|`, and in it with probability at least
`1 − (k_lo + k_hi)/|C|`, whatever the fits. -/
theorem cqrAsym_coverage {W Z X : Type*} [Fintype W] [Fintype Z] {C : Finset ι} {t : ι}
    (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (hP : ∀ σ ∈ permsOn C, ∀ w z, P w (z ∘ σ) = P w z) (pt : Z → X × ℝ) (qlo qhi : W → X → ℝ)
    (klo khi : ℕ) :
    1 - (klo : ℝ) / C.card ≤ ∑ w, ∑ z, P w z *
        (if Accepts klo C (fun i => cqrLoScore (qlo w) (pt (z i))) t then 1 else 0) ∧
      1 - (khi : ℝ) / C.card ≤ ∑ w, ∑ z, P w z *
        (if Accepts khi C (fun i => cqrHiScore (qhi w) (pt (z i))) t then 1 else 0) ∧
      1 - ((klo : ℝ) + khi) / C.card ≤ ∑ w, ∑ z, P w z *
        (if Accepts klo C (fun i => cqrLoScore (qlo w) (pt (z i))) t ∧
          Accepts khi C (fun i => cqrHiScore (qhi w) (pt (z i))) t then 1 else 0) := by
  have hlo : 1 - (klo : ℝ) / C.card ≤ ∑ w, ∑ z, P w z *
      (if Accepts klo C (fun i => cqrLoScore (qlo w) (pt (z i))) t then 1 else 0) :=
    split_coverage ht P hP0 hP1 hP (fun w p => cqrLoScore (qlo w) (pt p)) klo
  have hhi : 1 - (khi : ℝ) / C.card ≤ ∑ w, ∑ z, P w z *
      (if Accepts khi C (fun i => cqrHiScore (qhi w) (pt (z i))) t then 1 else 0) :=
    split_coverage ht P hP0 hP1 hP (fun w p => cqrHiScore (qhi w) (pt p)) khi
  have hand := sum_and_ge P hP0 hP1
    (fun w z => Accepts klo C (fun i => cqrLoScore (qlo w) (pt (z i))) t)
    (fun w z => Accepts khi C (fun i => cqrHiScore (qhi w) (pt (z i))) t)
  refine ⟨hlo, hhi, ?_⟩
  have : 1 - ((klo : ℝ) + khi) / C.card =
      (1 - (klo : ℝ) / C.card) + (1 - (khi : ℝ) / C.card) - 1 := by
    ring
  rw [this]
  exact le_trans (by linarith [hlo, hhi]) hand

/-- **Asymmetric CQR at `α/2` per tail.** At `k_lo, k_hi ≤ (α/2)|C|`, `⌊(α/2)(n₂ + 1)⌋` each, the
set covers with probability at least `1 − α`, each tail missing at most `α/2`. -/
theorem cqrAsym_coverage_half {W Z X : Type*} [Fintype W] [Fintype Z] {C : Finset ι} {t : ι}
    (ht : t ∈ C) (P : W → (ι → Z) → ℝ) (hP0 : ∀ w z, 0 ≤ P w z) (hP1 : ∑ w, ∑ z, P w z = 1)
    (hP : ∀ σ ∈ permsOn C, ∀ w z, P w (z ∘ σ) = P w z) (pt : Z → X × ℝ) (qlo qhi : W → X → ℝ)
    {klo khi : ℕ} {α : ℝ} (hlo : (klo : ℝ) ≤ α / 2 * C.card) (hhi : (khi : ℝ) ≤ α / 2 * C.card) :
    1 - α ≤ ∑ w, ∑ z, P w z *
      (if Accepts klo C (fun i => cqrLoScore (qlo w) (pt (z i))) t ∧
        Accepts khi C (fun i => cqrHiScore (qhi w) (pt (z i))) t then 1 else 0) := by
  have hc : (0 : ℝ) < C.card := by exact_mod_cast card_pos.2 ⟨t, ht⟩
  refine le_trans ?_ (cqrAsym_coverage ht P hP0 hP1 hP pt qlo qhi klo khi).2.2
  have : ((klo : ℝ) + khi) / C.card ≤ α := by
    rw [div_le_iff₀ hc]
    linarith
  linarith

/-- **The lower tail refuses below, the upper tail above.** With the calibration units held and the
test unit at `x`, the lower tail accepts an upper set of `y`. -/
theorem cqrLo_isUpperSet {ι₀ X : Type*} [Fintype ι₀] (qlo : X → ℝ) (cal : ι₀ → X × ℝ) (x : X)
    (k : ℕ) :
    IsUpperSet {y : ℝ |
      Accepts k univ (fun i : Option ι₀ => cqrLoScore qlo (i.elim (x, y) cal)) none} := by
  intro y y' hyy' hy
  refine accepts_mono ?_ (fun j hj => ?_) hy
  · simp only [Option.elim, cqrLoScore]
    linarith
  · cases j with
    | none => exact absurd rfl hj
    | some j => rfl

/-- The upper tail accepts a lower set of `y`. -/
theorem cqrHi_isLowerSet {ι₀ X : Type*} [Fintype ι₀] (qhi : X → ℝ) (cal : ι₀ → X × ℝ) (x : X)
    (k : ℕ) :
    IsLowerSet {y : ℝ |
      Accepts k univ (fun i : Option ι₀ => cqrHiScore qhi (i.elim (x, y) cal)) none} := by
  intro y y' hyy' hy
  refine accepts_mono ?_ (fun j hj => ?_) hy
  · simp only [Option.elim, cqrHiScore]
    linarith
  · cases j with
    | none => exact absurd rfl hj
    | some j => rfl

/-- **Asymmetric CQR's set is an interval**: an upper set from the lower tail meeting a lower set
from the upper tail. -/
theorem cqrAsym_ordConnected {ι₀ X : Type*} [Fintype ι₀] (qlo qhi : X → ℝ) (cal : ι₀ → X × ℝ)
    (x : X) (klo khi : ℕ) :
    Set.OrdConnected {y : ℝ |
      Accepts klo univ (fun i : Option ι₀ => cqrLoScore qlo (i.elim (x, y) cal)) none ∧
        Accepts khi univ (fun i : Option ι₀ => cqrHiScore qhi (i.elim (x, y) cal)) none} := by
  rw [Set.setOf_and]
  exact (cqrLo_isUpperSet qlo cal x klo).ordConnected.inter
    (cqrHi_isLowerSet qhi cal x khi).ordConnected

end AsymmetricCQR

end Conformal
